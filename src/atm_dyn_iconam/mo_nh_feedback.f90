!>
!!  This module contains the routines needed for nesting in the nonhydrostatic.
!!  version.
!!
!! @par Revision History
!!  Developed and tested by Guenther Zaengl, DWD (2010-02-10)
!!
!! @par Copyright
!! 2002-2006 by DWD and MPI-M
!! This software is provided for non-commercial use only.
!! See the LICENSE and the WARRANTY conditions.
!!
!! @par License
!! The use of ICON is hereby granted free of charge for an unlimited time,
!! provided the following rules are accepted and applied:
!! <ol>
!! <li> You may use or modify this code for your own non commercial and non
!!    violent purposes.
!! <li> The code may not be re-distributed without the consent of the authors.
!! <li> The copyright notice and statement of authorship must appear in all
!!    copies.
!! <li> You accept the warranty conditions (see WARRANTY).
!! <li> In case you intend to use the code commercially, we oblige you to sign
!!    an according license agreement with DWD and MPI-M.
!! </ol>
!!
!! @par Warranty
!! This code has been tested up to a certain level. Defects and weaknesses,
!! which may be included in the code, do not establish any warranties by the
!! authors.
!! The authors do not make any warranty, express or implied, or assume any
!! liability or responsibility for the use, acquisition or application of this
!! software.
!!
!!
MODULE mo_nh_feedback
!
!
USE mo_kind,                ONLY: wp
USE mo_exception,           ONLY: message_text, message
USE mo_model_domain,        ONLY: t_patch, t_grid_cells, t_grid_edges, p_patch_local_parent
USE mo_model_domain_import, ONLY: n_dom, n_dom_start
USE mo_interpolation,       ONLY: t_int_state, rbf_vec_interpol_vertex
USE mo_grf_interpolation,   ONLY: t_gridref_state, grf_velfbk, p_grf_state_local_parent
USE mo_nonhydrostatic_config, ONLY: l_masscorr_nest
USE mo_dynamics_config,     ONLY: nnow, nnew, nnow_rcf, nnew_rcf, nsav1, nsav2 
USE mo_parallel_config,     ONLY: nproma, p_test_run
USE mo_run_config,          ONLY: ltransport, iforcing, msg_level, ntracer
USE mo_nonhydro_types,      ONLY: t_nh_state, t_nh_prog, t_nh_diag
USE mo_impl_constants,      ONLY: min_rlcell, min_rledge, min_rlcell_int, min_rledge_int, &
            &                     min_rlvert_int, MAX_CHAR_LENGTH
USE mo_loopindices,         ONLY: get_indices_c, get_indices_e, get_indices_v
USE mo_impl_constants_grf,  ONLY: grf_fbk_start_c, grf_fbk_start_e,          &
                                  grf_bdywidth_c
USE mo_mpi,                 ONLY:  my_process_is_mpi_seq
USE mo_communication,       ONLY: exchange_data_mult
USE mo_sync,                ONLY: SYNC_C, SYNC_C1, SYNC_E, sync_patch_array, &
                                  global_sum_array3, sync_patch_array_mult
USE mo_physical_constants,  ONLY: rd, cvd_o_rd, p0ref

IMPLICIT NONE

PRIVATE

CHARACTER(len=*), PARAMETER :: version = '$Id$'

PUBLIC :: feedback, relax_feedback

CONTAINS


!>
!! This routine computes the feedback of the prognostic variables from the fine mesh.
!!
!! This routine computes the feedback of the prognostic variables from the fine mesh
!! to the corresponding grid point on the coarse mesh
!! jg in this case denotes the fine mesh level; output goes to parent_id(jg)
!!
!! @par Revision History
!! Developed  by Guenther Zaengl, DWD, 2008-04-15
!! @par
!! Modification by Guenther Zaengl, DWD, 2008-09-12:
!! Change feedback for cell-based variables from area-weighted averaging
!! to using fbk_wgt (see above routine)
!!
SUBROUTINE feedback(p_patch, p_nh_state, p_int_state, p_grf_state, jg, jgp, l_trac_fbk)

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_nh_feedback:feedback'


TYPE(t_patch),       TARGET, INTENT(IN)    ::  p_patch(n_dom_start:n_dom)
TYPE(t_nh_state), TARGET, INTENT(INOUT)    ::  p_nh_state(n_dom)
TYPE(t_int_state),   TARGET, INTENT(IN)    ::  p_int_state(n_dom_start:n_dom)
TYPE(t_gridref_state), TARGET, INTENT(IN)  ::  p_grf_state(n_dom_start:n_dom)

INTEGER, INTENT(IN) :: jg   ! child grid level
INTEGER, INTENT(IN) :: jgp  ! parent grid level

! Switch if feedback is done for tracers
! (when calling transport and microphysics not every dynamics time step, tracer feedback
!  should probably be restricted to transport time steps)
LOGICAL, INTENT(IN) :: l_trac_fbk

! local variables

TYPE(t_nh_prog),    POINTER     :: p_parent_prog => NULL()
TYPE(t_nh_prog),    POINTER     :: p_parent_save => NULL()
TYPE(t_nh_diag),    POINTER     :: p_parent_tend => NULL()
TYPE(t_nh_prog),    POINTER     :: p_child_prog => NULL()
TYPE(t_nh_prog),    POINTER     :: p_child_save => NULL()
TYPE(t_nh_diag),    POINTER     :: p_child_tend => NULL()
TYPE(t_nh_prog),    POINTER     :: p_parent_prog_rcf => NULL()
TYPE(t_nh_prog),    POINTER     :: p_child_prog_rcf  => NULL()
TYPE(t_grid_cells), POINTER     :: p_gcp => NULL()
TYPE(t_grid_cells), POINTER     :: p_gcc => NULL()
TYPE(t_grid_edges), POINTER     :: p_gep => NULL()
TYPE(t_grid_edges), POINTER     :: p_gec => NULL()
TYPE(t_gridref_state), POINTER  :: p_grf => NULL()
TYPE(t_int_state), POINTER      :: p_intc => NULL()
TYPE(t_patch),      POINTER     :: p_pp => NULL()
TYPE(t_patch),      POINTER     :: p_pc => NULL()

! Indices
INTEGER :: jb, jc, jk, jt, je, js, jgc, i_nchdom, i_chidx, &
           i_startblk, i_endblk, i_startidx, i_endidx, ic, i_ncd

INTEGER :: nlev_c, nlevp1_c  ! number of full and half levels (child dom)
INTEGER :: nlev_p, nlevp1_p  ! number of full and half levels (parent dom)
INTEGER :: nshift, nshift_c

REAL(wp), DIMENSION(nproma,p_patch(jg)%nlev,p_patch(jg)%nblks_v) :: z_u, z_v
REAL(wp) ::   &  ! RBF-reconstructed velocity
  &  vn_aux(nproma,p_patch(jg)%nlev,p_patch(jg)%nblks_e,2)
REAL(wp), ALLOCATABLE :: feedback_thv_tend(:,:,:)
REAL(wp), ALLOCATABLE :: feedback_rho_tend(:,:,:)
REAL(wp), ALLOCATABLE :: feedback_vn(:,:,:)
REAL(wp), ALLOCATABLE :: feedback_w_tend(:,:,:)
REAL(wp), ALLOCATABLE :: feedback_tracer_mass(:,:,:,:)
REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)   :: fbk_tend, parent_tend
REAL(wp), ALLOCATABLE, DIMENSION(:,:,:,:) :: fbk_tr_mass, parent_tr_mass
REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)   :: fbk_tr_totmass, parent_tr_totmass
REAL(wp) :: tendency_corr(p_patch(jgp)%nlev),        &
  &         aux_diff((ntracer+1)*p_patch(jgp)%nlev), &
  &         tracer_corr(ntracer*p_patch(jgp)%nlev)

REAL(wp) :: rd_o_cvd, rd_o_p0ref

INTEGER, DIMENSION(:,:,:), POINTER :: iidx, iblk, iidxv, iblkv
LOGICAL :: l_parallel
REAL(wp), DIMENSION(:,:,:), POINTER :: p_fbkwgt, p_fbkwgt_tr, p_fb_layer_thickness
REAL(wp), DIMENSION(:,:),   POINTER :: p_fbarea

!-----------------------------------------------------------------------

IF (msg_level >= 10) THEN
  WRITE(message_text,'(a,i2,a,i2)') '========= Feedback:',jg,' =>',jgp
  CALL message(TRIM(routine),message_text)
ENDIF

IF ( my_process_is_mpi_seq() ) THEN
  l_parallel = .FALSE.
ELSE
  l_parallel = .TRUE.
ENDIF


p_parent_prog    => p_nh_state(jgp)%prog(nnew(jgp))
p_parent_prog_rcf=> p_nh_state(jgp)%prog(nnew_rcf(jgp))
p_parent_save    => p_nh_state(jgp)%prog(nsav1(jgp))
p_parent_tend    => p_nh_state(jgp)%diag
p_child_prog     => p_nh_state(jg)%prog(nnow(jg))
p_child_prog_rcf => p_nh_state(jg)%prog(nnow_rcf(jg))
p_child_save     => p_nh_state(jg)%prog(nsav2(jg))
p_child_tend     => p_nh_state(jg)%diag
p_intc           => p_int_state(jg)
p_gcc            => p_patch(jg)%cells
p_gec            => p_patch(jg)%edges
p_pc             => p_patch(jg)

IF(l_parallel) THEN
  p_grf => p_grf_state_local_parent(jg)
  p_gcp => p_patch_local_parent(jg)%cells
  p_gep => p_patch_local_parent(jg)%edges
  p_pp  => p_patch_local_parent(jg)
ELSE
  p_grf => p_grf_state(jgp)
  p_gcp => p_patch(jgp)%cells
  p_gep => p_patch(jgp)%edges
  p_pp  => p_patch(jgp)
ENDIF

nlev_c   = p_pc%nlev
nlevp1_c = p_pc%nlevp1

nlev_p   = p_pp%nlev
nlevp1_p = p_pp%nlevp1
 
nshift = p_pc%nshift
js     = nshift

i_nchdom = MAX(1,p_pc%n_childdom)
i_chidx  = p_pc%parent_child_index

! R/c_v (not present in physical constants)
rd_o_cvd = 1._wp / cvd_o_rd

! R / p0ref
rd_o_p0ref = rd / p0ref

! parent_tend, parent_tr_mass, parent_tr_totmass are always calculated on the global parent
! and thus have to be allocated within global parent limits

i_startblk = p_patch(jgp)%cells%start_blk(grf_fbk_start_c,i_chidx)
i_endblk   = p_patch(jgp)%cells%end_blk(min_rlcell_int,i_chidx)

ALLOCATE(parent_tend(nproma, nlev_p, i_startblk:i_endblk))
IF (ltransport .AND. l_trac_fbk) &
   ALLOCATE(parent_tr_mass(nproma, nlev_p, i_startblk:i_endblk, ntracer), &
            parent_tr_totmass(nproma, nlev_p, i_startblk:i_endblk) )

! Allocation of storage fields
! In parallel runs the lower bound of the feedback_* arrays must be 1 for use in exchange data,
! the lower bound of the fbk_* arrays must be i_startblk for the use in global sum.

i_startblk = p_gcp%start_blk(grf_fbk_start_c,i_chidx)
i_endblk   = p_gcp%end_blk(min_rlcell_int,i_chidx)

ALLOCATE(fbk_tend(nproma, nlev_p,  i_startblk:i_endblk))
IF(ltransport .AND. l_trac_fbk) &
  ALLOCATE(fbk_tr_mass(nproma, nlev_p, i_startblk:i_endblk, ntracer), &
           fbk_tr_totmass(nproma, nlev_p, i_startblk:i_endblk)  )

IF(l_parallel) i_startblk = 1

ALLOCATE(feedback_thv_tend  (nproma, nlev_p, i_startblk:i_endblk),   &
         feedback_rho_tend  (nproma, nlev_p, i_startblk:i_endblk),   &
         feedback_w_tend    (nproma, nlevp1_p, i_startblk:i_endblk))

IF(ltransport .AND. l_trac_fbk) &
  ALLOCATE(feedback_tracer_mass(nproma, nlev_p, i_startblk:i_endblk, ntracer))

i_startblk = p_gep%start_blk(grf_fbk_start_e,i_chidx)
i_endblk   = p_gep%end_blk(min_rledge,i_chidx)

IF(l_parallel) i_startblk = 1

ALLOCATE(feedback_vn(nproma, nlev_p, i_startblk:i_endblk))


! Set pointers to index and coefficient fields for cell-based variables
iidx => p_gcp%child_idx
iblk => p_gcp%child_blk

p_fbkwgt    => p_grf%fbk_wgt_c
p_fbkwgt_tr => p_grf%fbk_wgt_ct
p_fbarea    => p_gcp%area

IF(l_parallel) THEN
  p_fb_layer_thickness => p_gcp%ddqz_z_full
ELSE
  p_fb_layer_thickness => p_nh_state(jgp)%metrics%ddqz_z_full
ENDIF

! Preparation of feedback: compute child tendencies


!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)
!
! Part 1: cell-based variables
i_startblk = p_gcc%start_blk(3,1)
i_endblk   = p_gcc%end_blk(min_rlcell_int,i_nchdom)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk,jt)
DO jb = i_startblk, i_endblk

  CALL get_indices_c(p_pc, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, 3, min_rlcell_int)

  DO jk = 1, nlev_c
    DO jc = i_startidx, i_endidx

      p_child_tend%grf_tend_rho(jc,jk,jb) = &
        p_child_prog%rho(jc,jk,jb) - p_child_save%rho(jc,jk,jb)

      p_child_tend%grf_tend_thv(jc,jk,jb) = &
        p_child_prog%theta_v(jc,jk,jb) - p_child_save%theta_v(jc,jk,jb)

      p_child_tend%grf_tend_w(jc,jk,jb) = &
        p_child_prog%w(jc,jk,jb) - p_child_save%w(jc,jk,jb)
    ENDDO
  ENDDO

  DO jc = i_startidx, i_endidx
    p_child_tend%grf_tend_w(jc,nlevp1_c,jb) = &
      p_child_prog%w(jc,nlevp1_c,jb) - p_child_save%w(jc,nlevp1_c,jb)
  ENDDO

ENDDO
!$OMP END DO

! Compute feedback tendency for density and tracers, including a layer-wise conservation correction

! Start/End block in the parent domain
i_startblk = p_gcp%start_blk(grf_fbk_start_c,i_chidx)
i_endblk   = p_gcp%end_blk(min_rlcell_int,i_chidx)


!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk,jt)
DO jb = i_startblk, i_endblk

  CALL get_indices_c(p_pp, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, grf_fbk_start_c, min_rlcell_int, i_chidx)

  fbk_tend(:,:,jb) = 0._wp

#ifdef __LOOP_EXCHANGE
  DO jc = i_startidx, i_endidx
    DO jk = 1, nlev_c
#else
!CDIR UNROLL=8
  DO jk = 1, nlev_c
    DO jc = i_startidx, i_endidx
#endif

      feedback_rho_tend(jc,jk+js,jb) =                                                &
        p_child_tend%grf_tend_rho(iidx(jc,jb,1),jk,iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
        p_child_tend%grf_tend_rho(iidx(jc,jb,2),jk,iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
        p_child_tend%grf_tend_rho(iidx(jc,jb,3),jk,iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
        p_child_tend%grf_tend_rho(iidx(jc,jb,4),jk,iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

      fbk_tend(jc,jk+js,jb) = feedback_rho_tend(jc,jk+js,jb)*p_fbarea(jc,jb)* &
                           p_fb_layer_thickness(jc,jk+js,jb)

    ENDDO
  ENDDO

  ! Tracers
  IF (ltransport .AND. l_trac_fbk) THEN

    fbk_tr_totmass(:,:,jb) = 0._wp

    DO jt = 1, ntracer

      fbk_tr_mass(:,:,jb,jt) = 0._wp

      feedback_tracer_mass(:,1:nshift,jb,jt) = 0._wp

#ifdef __LOOP_EXCHANGE
    DO jc = i_startidx, i_endidx
      DO jk = 1, nlev_c
#else
!CDIR UNROLL=8
      DO jk = 1, nlev_c
        DO jc = i_startidx, i_endidx
#endif

          feedback_tracer_mass(jc,jk+js,jb,jt) =                                        &
          p_fbkwgt_tr(jc,jb,1)*p_child_prog%rho(iidx(jc,jb,1),jk,iblk(jc,jb,1))*        &
          p_child_prog_rcf%tracer(iidx(jc,jb,1),jk,iblk(jc,jb,1),jt) +                  &
          p_fbkwgt_tr(jc,jb,2)*p_child_prog%rho(iidx(jc,jb,2),jk,iblk(jc,jb,2))*        &
          p_child_prog_rcf%tracer(iidx(jc,jb,2),jk,iblk(jc,jb,2),jt) +                  &
          p_fbkwgt_tr(jc,jb,3)*p_child_prog%rho(iidx(jc,jb,3),jk,iblk(jc,jb,3))*        &
          p_child_prog_rcf%tracer(iidx(jc,jb,3),jk,iblk(jc,jb,3),jt) +                  &
          p_fbkwgt_tr(jc,jb,4)*p_child_prog%rho(iidx(jc,jb,4),jk,iblk(jc,jb,4))*        &
          p_child_prog_rcf%tracer(iidx(jc,jb,4),jk,iblk(jc,jb,4),jt)

          fbk_tr_mass(jc,jk+js,jb,jt) = feedback_tracer_mass(jc,jk+js,jb,jt)        &
            &                       * p_fbarea(jc,jb)*p_fb_layer_thickness(jc,jk+js,jb)

          fbk_tr_totmass(jc,jk+js,jb) = fbk_tr_totmass(jc,jk+js,jb)                 &
            &                       + fbk_tr_mass(jc,jk+js,jb,jt)

        ENDDO
      ENDDO
    ENDDO

  ENDIF

ENDDO
!$OMP END DO

! Calculate the area-weighted sum of (p_parent_prog%rho-p_parent_save%rho)
! as well as the feedback area - this has always to be done in the parent patch

! Start/End block in the parent domain
i_startblk = p_patch(jgp)%cells%start_blk(grf_fbk_start_c,i_chidx)
i_endblk   = p_patch(jgp)%cells%end_blk(min_rlcell_int,i_chidx)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc)
DO jb = i_startblk, i_endblk

  parent_tend(:,:,jb) = 0._wp

  CALL get_indices_c(p_patch(jgp), jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, grf_fbk_start_c, min_rlcell_int, i_chidx)

  DO jk = 1, nlev_p
    DO jc = i_startidx, i_endidx
      parent_tend(jc,jk,jb) =  &
        (p_parent_prog%rho(jc,jk,jb) - p_parent_save%rho(jc,jk,jb))*   &
         p_patch(jgp)%cells%area(jc,jb)*p_nh_state(jgp)%metrics%ddqz_z_full(jc,jk,jb)
    ENDDO
! obsolete because halo points are now all at the end
!   WHERE(.NOT.p_patch(jgp)%cells%owner_mask(:,jb)) parent_tend(:,jk,jb) = 0._wp
  ENDDO

  IF (ltransport .AND. l_trac_fbk) THEN

    parent_tr_totmass(:,:,jb) = 0._wp

    DO jt = 1, ntracer

      parent_tr_mass(:,:,jb,jt) = 0._wp

      DO jk = 1, nlev_p
        DO jc = i_startidx, i_endidx
          parent_tr_mass(jc,jk,jb,jt) =  &
            p_parent_prog%rho(jc,jk,jb)*p_parent_prog_rcf%tracer(jc,jk,jb,jt)* &
            p_patch(jgp)%cells%area(jc,jb)*p_nh_state(jgp)%metrics%ddqz_z_full(jc,jk,jb)

          parent_tr_totmass(jc,jk,jb) = parent_tr_totmass(jc,jk,jb) + &
                                        parent_tr_mass(jc,jk,jb,jt)
        ENDDO
! obsolete because halo points are now all at the end
!       WHERE(.NOT.p_patch(jgp)%cells%owner_mask(:,jb))
!         parent_tr_mass(:,jk,jb,jt) = 0._wp
!         parent_tr_totmass(:,jk,jb) = 0._wp
!       END WHERE
      ENDDO
    ENDDO
  ENDIF

ENDDO
!$OMP END DO
!$OMP END PARALLEL

! fbk_dom_volume is now set in p_nh_state(jg)%metrics

IF ( .NOT. (ltransport .AND. l_trac_fbk)) THEN
  ! compute conservation correction for global mass only
  aux_diff(1:nlev_p) = global_sum_array3(1,.TRUE.,parent_tend,fbk_tend,diffmask=(/1/))
  DO jk = 1, nlev_p
    tendency_corr(jk) = aux_diff(jk) / p_nh_state(jg)%metrics%fbk_dom_volume(jk)
  ENDDO
ELSE IF (iforcing <= 1) THEN
  ! compute conservation correction for global mass and each tracer separately
  ! the correction is additive for density and multiplicative for tracer masses
  aux_diff = global_sum_array3(ntracer+1,.TRUE.,parent_tend,fbk_tend,f4din=parent_tr_mass,&
                               f4dd=fbk_tr_mass,diffmask=(/1,(2,jt=1,ntracer)/))
  DO jk = 1, nlev_p
    tendency_corr(jk) = aux_diff(jk) / p_nh_state(jg)%metrics%fbk_dom_volume(jk)
  ENDDO
  DO jk = 1, nlev_p*ntracer
    tracer_corr(jk) = aux_diff(nlev_p+jk)
   ENDDO
ELSE ! iforcing >= 2; tracers represent moisture variables
  ! compute conservation correction for global mass and total tracer mass
  ! the correction is additive for density and multiplicative for tracer mass
  aux_diff(1:2*nlev_p) = global_sum_array3(2,.TRUE.,parent_tend,fbk_tend,              &
                                         f3din2=parent_tr_totmass,f3dd2=fbk_tr_totmass,&
                                         diffmask=(/1,2/))
  DO jk = 1, nlev_p
    tendency_corr(jk) = aux_diff(jk) / p_nh_state(jg)%metrics%fbk_dom_volume(jk)
    tracer_corr(jk)   = aux_diff(nlev_p+jk)
  ENDDO
ENDIF

!$OMP PARALLEL PRIVATE(i_startblk,i_endblk,nshift_c,i_ncd,ic,jgc)

IF (l_masscorr_nest) THEN
  ! Add mass conservation correction to child domain in order to prevent
  ! possible inconsistencies in the mass fields; without l_nest_rcf, the
  ! correction increments are smaller, so that this correction does not appear
  ! to be needed.
  i_startblk = p_gcc%start_blk(grf_bdywidth_c+1,1)
  i_endblk   = p_gcc%end_blk(min_rlcell,i_nchdom)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk)
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(p_pc, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, grf_bdywidth_c+1, min_rlcell)

    DO jk = 1, nlev_c
      DO jc = i_startidx, i_endidx

        p_child_prog%rho(jc,jk,jb) = p_child_prog%rho(jc,jk,jb) + tendency_corr(jk)

        p_child_prog%rhotheta_v(jc,jk,jb) = p_child_prog%rho(jc,jk,jb) * &
                                            p_child_prog%theta_v(jc,jk,jb)

        p_child_prog%exner(jc,jk,jb) = &
          EXP(rd_o_cvd*LOG(rd_o_p0ref*p_child_prog%rhotheta_v(jc,jk,jb)))

      ENDDO
    ENDDO

  ENDDO
!$OMP END DO

  ! The conservation correction also needs to be applied to all nested domains
  IF (p_pc%n_childdom > 0) THEN

    DO ic = 1, p_pc%n_chd_total
      jgc = p_pc%child_id_list(ic)

      i_ncd      = MAX(1,p_patch(jgc)%n_childdom)
      i_startblk = p_patch(jgc)%cells%start_blk(grf_bdywidth_c+1,1)
      i_endblk   = p_patch(jgc)%cells%end_blk(min_rlcell,i_ncd)
      nshift_c   = p_patch(jgc)%nshift_total - p_pc%nshift_total

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk)
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(p_patch(jgc), jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, grf_bdywidth_c+1, min_rlcell)

        DO jk = 1, p_patch(jgc)%nlev
          DO jc = i_startidx, i_endidx

            p_nh_state(jgc)%prog(nnow(jgc))%rho(jc,jk,jb) = &
              p_nh_state(jgc)%prog(nnow(jgc))%rho(jc,jk,jb) + tendency_corr(jk+nshift_c)

            p_nh_state(jgc)%prog(nnow(jgc))%rhotheta_v(jc,jk,jb) = &
              p_nh_state(jgc)%prog(nnow(jgc))%rho(jc,jk,jb) * &
              p_nh_state(jgc)%prog(nnow(jgc))%theta_v(jc,jk,jb)

            p_nh_state(jgc)%prog(nnow(jgc))%exner(jc,jk,jb) = &
              EXP(rd_o_cvd*LOG(rd_o_p0ref*p_nh_state(jgc)%prog(nnow(jgc))%rhotheta_v(jc,jk,jb)))

          ENDDO
        ENDDO

      ENDDO
!$OMP END DO

    ENDDO
  ENDIF
ENDIF

i_startblk = p_gcp%start_blk(grf_fbk_start_c,i_chidx)
i_endblk   = p_gcp%end_blk(min_rlcell_int,i_chidx)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk)
DO jb = i_startblk, i_endblk

  CALL get_indices_c(p_pp, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, grf_fbk_start_c, min_rlcell_int, i_chidx)

#ifdef __LOOP_EXCHANGE
  DO jc = i_startidx, i_endidx
    DO jk = 1, nlev_c
#else
!CDIR UNROLL=4
  DO jk = 1, nlev_c
    DO jc = i_startidx, i_endidx
#endif

      feedback_rho_tend(jc,jk+js,jb) = feedback_rho_tend(jc,jk+js,jb) + tendency_corr(jk+js)

      feedback_thv_tend(jc,jk+js,jb) =                                             &
        p_child_tend%grf_tend_thv(iidx(jc,jb,1),jk,iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
        p_child_tend%grf_tend_thv(iidx(jc,jb,2),jk,iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
        p_child_tend%grf_tend_thv(iidx(jc,jb,3),jk,iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
        p_child_tend%grf_tend_thv(iidx(jc,jb,4),jk,iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

      feedback_w_tend(jc,jk+js,jb) =                                             &
        p_child_tend%grf_tend_w(iidx(jc,jb,1),jk,iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
        p_child_tend%grf_tend_w(iidx(jc,jb,2),jk,iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
        p_child_tend%grf_tend_w(iidx(jc,jb,3),jk,iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
        p_child_tend%grf_tend_w(iidx(jc,jb,4),jk,iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

    ENDDO
  ENDDO

  DO jc = i_startidx, i_endidx
    feedback_w_tend(jc,nlevp1_p,jb) =                                             &
      p_child_tend%grf_tend_w(iidx(jc,jb,1),nlevp1_c,iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
      p_child_tend%grf_tend_w(iidx(jc,jb,2),nlevp1_c,iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
      p_child_tend%grf_tend_w(iidx(jc,jb,3),nlevp1_c,iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
      p_child_tend%grf_tend_w(iidx(jc,jb,4),nlevp1_c,iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)
  ENDDO

ENDDO
!$OMP END DO

IF (.NOT. l_parallel) THEN

  ! Add feedback tendencies to prognostic variables on parent grid

!$OMP DO PRIVATE(jb,i_startidx,i_endidx)
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(p_pp, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, grf_fbk_start_c, min_rlcell_int, i_chidx)

    p_parent_prog%rho(i_startidx:i_endidx,nshift+1:nlev_p,jb) =   &
      p_parent_save%rho(i_startidx:i_endidx,nshift+1:nlev_p,jb) + &
      feedback_rho_tend(i_startidx:i_endidx,nshift+1:nlev_p,jb)

    p_parent_prog%theta_v(i_startidx:i_endidx,nshift+1:nlev_p,jb) =   &
      p_parent_save%theta_v(i_startidx:i_endidx,nshift+1:nlev_p,jb) + &
      feedback_thv_tend(i_startidx:i_endidx,nshift+1:nlev_p,jb)

    ! rhotheta_v is diagnosed from rho and theta_v
    p_parent_prog%rhotheta_v(i_startidx:i_endidx,nshift+1:nlev_p,jb) =   &
      p_parent_prog%theta_v(i_startidx:i_endidx,nshift+1:nlev_p,jb) *   &
      p_parent_prog%rho(i_startidx:i_endidx,nshift+1:nlev_p,jb)

    ! exner is diagnosed from rhotheta_v
    p_parent_prog%exner(i_startidx:i_endidx,nshift+1:nlev_p,jb) = EXP(rd_o_cvd*LOG(  &
      rd_o_p0ref*p_parent_prog%rhotheta_v(i_startidx:i_endidx,nshift+1:nlev_p,jb)))

    p_parent_prog%w(i_startidx:i_endidx,nshift+1:nlevp1_p,jb) =   &
      p_parent_save%w(i_startidx:i_endidx,nshift+1:nlevp1_p,jb) + &
      feedback_w_tend(i_startidx:i_endidx,nshift+1:nlevp1_p,jb)

    IF (ltransport .AND. l_trac_fbk) THEN ! perform tracer feedback
      IF (iforcing <= 1) THEN
        DO jt = 1, ntracer
          DO jk = nshift+1, nlev_p
            p_parent_prog_rcf%tracer(i_startidx:i_endidx,jk,jb,jt) =                 &
              feedback_tracer_mass(i_startidx:i_endidx,jk,jb,jt) *                   &
              tracer_corr(jk+(jt-1)*nlev_p)/p_parent_prog%rho(i_startidx:i_endidx,jk,jb)
          ENDDO
        ENDDO
      ELSE
        DO jt = 1, ntracer
          DO jk = nshift+1, nlev_p
            p_parent_prog_rcf%tracer(i_startidx:i_endidx,jk,jb,jt) =         &
              feedback_tracer_mass(i_startidx:i_endidx,jk,jb,jt) *           &
              tracer_corr(jk)/p_parent_prog%rho(i_startidx:i_endidx,jk,jb)
          ENDDO
        ENDDO
      ENDIF
    ENDIF

  ENDDO
!$OMP END DO

ENDIF
!$OMP END PARALLEL


IF (grf_velfbk == 2) THEN ! Interpolate velocity tendencies in child domain to vertices
  CALL rbf_vec_interpol_vertex( p_child_prog%vn, p_pc, p_intc, z_u, z_v)
ENDIF

! Set pointers to index and coefficient fields
! (this needs to be done outside a parallel section!)
iidx => p_gep%child_idx
iblk => p_gep%child_blk

p_fbkwgt => p_grf%fbk_wgt_e

iidxv => p_gec%vertex_idx
iblkv => p_gec%vertex_blk

!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)

! Velocity feedback
IF (grf_velfbk == 1) THEN ! Averaging weighted with child edge lenghts

  i_startblk = p_gep%start_blk(grf_fbk_start_e,i_chidx)
  i_endblk   = p_gep%end_blk(min_rledge_int,i_chidx)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,je,jk)
  DO jb = i_startblk, i_endblk

    CALL get_indices_e(p_pp, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, grf_fbk_start_e, min_rledge_int, i_chidx)

    IF (l_parallel) THEN
      feedback_vn(:,1:nshift,jb) = 0._wp
    ENDIF

#ifdef __LOOP_EXCHANGE
    DO je = i_startidx, i_endidx
      DO jk = 1, nlev_c
#else
!CDIR UNROLL=8
    DO jk = 1, nlev_c
      DO je = i_startidx, i_endidx
#endif

        feedback_vn(je,jk+js,jb) =                                            &
          p_child_prog%vn(iidx(je,jb,1),jk,iblk(je,jb,1))*p_fbkwgt(je,jb,1) + &
          p_child_prog%vn(iidx(je,jb,2),jk,iblk(je,jb,2))*p_fbkwgt(je,jb,2)
      ENDDO
    ENDDO

  ENDDO
!$OMP END DO

ELSE IF (grf_velfbk == 2) THEN ! Second-order interpolation of normal velocities
                               ! using RBF reconstruction to child vertices

  ! Projection of reconstructed velocities to the edge normals
  i_startblk = p_gec%start_blk(4,1)
  i_endblk   = p_gec%end_blk(min_rledge,i_nchdom)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,je,jk)
  DO jb = i_startblk, i_endblk

    CALL get_indices_e(p_pc, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, 4, min_rledge)

#ifdef __LOOP_EXCHANGE
    DO je = i_startidx, i_endidx
      DO jk = 1, nlev_c
#else
!CDIR UNROLL=6
    DO jk = 1, nlev_c
      DO je = i_startidx, i_endidx
#endif

        vn_aux(je,jk,jb,1) = z_u(iidxv(je,jb,1),jk,iblkv(je,jb,1)) &
             &             * p_gec%primal_normal_vert(je,jb,1)%v1  &
             &             + z_v(iidxv(je,jb,1),jk,iblkv(je,jb,1)) &
             &             * p_gec%primal_normal_vert(je,jb,1)%v2

        vn_aux(je,jk,jb,2) = z_u(iidxv(je,jb,2),jk,iblkv(je,jb,2)) &
             &             * p_gec%primal_normal_vert(je,jb,2)%v1  &
             &             + z_v(iidxv(je,jb,2),jk,iblkv(je,jb,2)) &
             &             * p_gec%primal_normal_vert(je,jb,2)%v2
      ENDDO
    ENDDO
  ENDDO
!$OMP END DO

  i_startblk = p_gep%start_blk(grf_fbk_start_e,i_chidx)
  i_endblk   = p_gep%end_blk(min_rledge_int,i_chidx)


!$OMP DO PRIVATE(jb,i_startidx,i_endidx,je,jk)
  DO jb = i_startblk, i_endblk

    CALL get_indices_e(p_pp, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, grf_fbk_start_e, min_rledge_int, i_chidx)

    IF (l_parallel) THEN
      feedback_vn(:,1:nshift,jb) = 0._wp
    ENDIF

#ifdef __LOOP_EXCHANGE
    DO je = i_startidx, i_endidx
      DO jk = 1, nlev_c
#else
!CDIR UNROLL=6
    DO jk = 1, nlev_c
      DO je = i_startidx, i_endidx
#endif

        feedback_vn(je,jk+js,jb) =                                             &
        ( p_fbkwgt(je,jb,1)*(p_child_prog%vn(iidx(je,jb,1),jk,iblk(je,jb,1)) + &
                             p_child_prog%vn(iidx(je,jb,2),jk,iblk(je,jb,2)))+ &
          p_fbkwgt(je,jb,3)*vn_aux(iidx(je,jb,1),jk,iblk(je,jb,1),1) +         &
          p_fbkwgt(je,jb,4)*vn_aux(iidx(je,jb,1),jk,iblk(je,jb,1),2) +         &
          p_fbkwgt(je,jb,5)*vn_aux(iidx(je,jb,2),jk,iblk(je,jb,2),1) +         &
          p_fbkwgt(je,jb,6)*vn_aux(iidx(je,jb,2),jk,iblk(je,jb,2),2) )
      ENDDO
    ENDDO
  ENDDO
!$OMP END DO

ENDIF

IF (.NOT. l_parallel) THEN

  ! Add feedback tendencies to prognostic variables on parent grid

!$OMP DO PRIVATE(jb,i_startidx,i_endidx)
  DO jb = i_startblk, i_endblk

    CALL get_indices_e(p_pp, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, grf_fbk_start_e, min_rledge_int, i_chidx)

    p_parent_prog%vn(i_startidx:i_endidx,nshift+1:nlev_p,jb) =   &
      feedback_vn(i_startidx:i_endidx,nshift+1:nlev_p,jb)

  ENDDO
!$OMP END DO

ENDIF

!$OMP END PARALLEL

IF (l_parallel) THEN

  CALL exchange_data_mult(p_pp%comm_pat_loc_to_glb_c_fbk, 3, 3*nlev_c+1,        &
                          RECV1=p_parent_prog%rho, SEND1=feedback_rho_tend,     &
                          ADD1= p_parent_save%rho, RECV2=p_parent_prog%theta_v, &
                          SEND2=feedback_thv_tend, ADD2= p_parent_save%theta_v, &
                          RECV3=p_parent_prog%w,   SEND3=feedback_w_tend,       &
                          ADD3=p_parent_save%w,    nshift=nshift )


  CALL exchange_data_mult(p_pp%comm_pat_loc_to_glb_e_fbk, 1, nlev_c, &
                          RECV1=p_parent_prog%vn, SEND1=feedback_vn, &
                          nshift=nshift)

IF (ltransport .AND. l_trac_fbk) THEN

  CALL exchange_data_mult(p_pp%comm_pat_loc_to_glb_c_fbk, ntracer, ntracer*nlev_c, &
    &                     RECV4D=p_parent_prog_rcf%tracer,                         &
    &                     SEND4D=feedback_tracer_mass, nshift=nshift)

ENDIF

  ! Compute rhotheta and exner in feedback area

  i_startblk = p_patch(jgp)%cells%start_blk(grf_fbk_start_c,i_chidx)
  i_endblk   = p_patch(jgp)%cells%end_blk(min_rlcell_int,i_chidx)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx)
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(p_patch(jgp), jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, grf_fbk_start_c, min_rlcell_int, i_chidx)

    ! rhotheta_v is diagnosed from rho and theta_v
    p_parent_prog%rhotheta_v(i_startidx:i_endidx,nshift+1:nlev_p,jb) =   &
      p_parent_prog%theta_v(i_startidx:i_endidx,nshift+1:nlev_p,jb) *   &
      p_parent_prog%rho(i_startidx:i_endidx,nshift+1:nlev_p,jb)

    ! exner is diagnosed from rhotheta_v
    p_parent_prog%exner(i_startidx:i_endidx,nshift+1:nlev_p,jb) = EXP(rd_o_cvd*LOG(  &
      rd_o_p0ref*p_parent_prog%rhotheta_v(i_startidx:i_endidx,nshift+1:nlev_p,jb)))

    ! divide tracer density (which is the feedback quantity) by air density,
    ! and apply multiplicative mass conservation correction
    IF (ltransport .AND. l_trac_fbk) THEN
      IF (iforcing <= 1) THEN
        DO jt = 1, ntracer
          DO jk = nshift+1, nlev_p
            p_parent_prog_rcf%tracer(i_startidx:i_endidx,jk,jb,jt) =                 &
              p_parent_prog_rcf%tracer(i_startidx:i_endidx,jk,jb,jt) *               &
              tracer_corr(jk+(jt-1)*nlev_p)/p_parent_prog%rho(i_startidx:i_endidx,jk,jb)
          ENDDO
        ENDDO
      ELSE
        DO jt = 1, ntracer
          DO jk = nshift+1, nlev_p
            p_parent_prog_rcf%tracer(i_startidx:i_endidx,jk,jb,jt) =         &
              p_parent_prog_rcf%tracer(i_startidx:i_endidx,jk,jb,jt) *       &
              tracer_corr(jk)/p_parent_prog%rho(i_startidx:i_endidx,jk,jb)
          ENDDO
        ENDDO
      ENDIF
    ENDIF

  ENDDO
!$OMP END DO

! Recompute rhotheta and exner also on the halo points
!$OMP DO PRIVATE(jk,jc,jb,ic)
#ifdef __LOOP_EXCHANGE
  DO ic = 1, p_nh_state(jgp)%metrics%ovlp_halo_c_dim(i_chidx)
    jc = p_nh_state(jgp)%metrics%ovlp_halo_c_idx(ic,i_chidx)
    jb = p_nh_state(jgp)%metrics%ovlp_halo_c_blk(ic,i_chidx)
    DO jk = nshift+1, nlev_p
#else
  DO jk = nshift+1, nlev_p
!CDIR NODEP,VOVERTAKE,VOB
    DO ic = 1, p_nh_state(jgp)%metrics%ovlp_halo_c_dim(i_chidx)
      jc = p_nh_state(jgp)%metrics%ovlp_halo_c_idx(ic,i_chidx)
      jb = p_nh_state(jgp)%metrics%ovlp_halo_c_blk(ic,i_chidx)
#endif

      p_parent_prog%rhotheta_v(jc,jk,jb) =   &
        p_parent_prog%theta_v(jc,jk,jb) * p_parent_prog%rho(jc,jk,jb)

      ! exner is diagnosed from rhotheta_v
      p_parent_prog%exner(jc,jk,jb) =   &
        EXP(rd_o_cvd*LOG(rd_o_p0ref*p_parent_prog%rhotheta_v(jc,jk,jb)))

    ENDDO
  ENDDO
!$OMP END DO

! Recompute tracer also on the halo points
  IF (ltransport .AND. l_trac_fbk) THEN
    IF (iforcing <= 1) THEN
!$OMP DO PRIVATE(jk,jt,jc,jb,ic)
#ifdef __LOOP_EXCHANGE
      DO ic = 1, p_nh_state(jgp)%metrics%ovlp_halo_c_dim(i_chidx)
        jc = p_nh_state(jgp)%metrics%ovlp_halo_c_idx(ic,i_chidx)
        jb = p_nh_state(jgp)%metrics%ovlp_halo_c_blk(ic,i_chidx)
        DO jt = 1, ntracer
          DO jk = nshift+1, nlev_p
#else
      DO jt = 1, ntracer
        DO jk = nshift+1, nlev_p
!CDIR NODEP,VOVERTAKE,VOB
          DO ic = 1, p_nh_state(jgp)%metrics%ovlp_halo_c_dim(i_chidx)
            jc = p_nh_state(jgp)%metrics%ovlp_halo_c_idx(ic,i_chidx)
            jb = p_nh_state(jgp)%metrics%ovlp_halo_c_blk(ic,i_chidx)
#endif

            p_parent_prog_rcf%tracer(jc,jk,jb,jt) =                 &
              p_parent_prog_rcf%tracer(jc,jk,jb,jt) *               &
              tracer_corr(jk+(jt-1)*nlev_p)/p_parent_prog%rho(jc,jk,jb)

          ENDDO
        ENDDO
      ENDDO
!$OMP END DO
    ELSE ! iforcing > 1
!$OMP DO PRIVATE(jk,jt,jc,jb,ic)
#ifdef __LOOP_EXCHANGE
      DO ic = 1, p_nh_state(jgp)%metrics%ovlp_halo_c_dim(i_chidx)
        jc = p_nh_state(jgp)%metrics%ovlp_halo_c_idx(ic,i_chidx)
        jb = p_nh_state(jgp)%metrics%ovlp_halo_c_blk(ic,i_chidx)
        DO jt = 1, ntracer
          DO jk = nshift+1, nlev_p
#else
      DO jt = 1, ntracer
        DO jk = nshift+1, nlev_p
!CDIR NODEP,VOVERTAKE,VOB
          DO ic = 1, p_nh_state(jgp)%metrics%ovlp_halo_c_dim(i_chidx)
            jc = p_nh_state(jgp)%metrics%ovlp_halo_c_idx(ic,i_chidx)
            jb = p_nh_state(jgp)%metrics%ovlp_halo_c_blk(ic,i_chidx)
#endif

            p_parent_prog_rcf%tracer(jc,jk,jb,jt) =         &
              p_parent_prog_rcf%tracer(jc,jk,jb,jt) *       &
              tracer_corr(jk)/p_parent_prog%rho(jc,jk,jb)

          ENDDO
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
  ENDIF

!$OMP END PARALLEL

ENDIF ! l_parallel


DEALLOCATE(parent_tend, fbk_tend, feedback_thv_tend, feedback_rho_tend, &
           feedback_w_tend, feedback_vn)
IF (ltransport .AND. l_trac_fbk) &
  DEALLOCATE(feedback_tracer_mass,parent_tr_mass,parent_tr_totmass,fbk_tr_mass,fbk_tr_totmass)

END SUBROUTINE feedback



!>
!! This routine computes the feedback of the prognostic variables from the fine mesh.
!!
!! This routine computes the feedback of the prognostic variables from the fine mesh
!! to the corresponding grid point on the coarse mesh
!! jg in this case denotes the fine mesh level; output goes to parent_id(jg)
!!
!! @par Revision History
!! Developed  by Guenther Zaengl, DWD, 2008-04-15
!! @par
!! Modification by Guenther Zaengl, DWD, 2008-09-12:
!! Change feedback for cell-based variables from area-weighted averaging
!! to using fbk_wgt (see above routine)
!!
SUBROUTINE relax_feedback(p_patch, p_nh_state, p_int_state, p_grf_state, jg, jgp, l_trac_fbk)

  CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_nh_feedback:relax_feedback'


  TYPE(t_patch),       TARGET, INTENT(IN)    ::  p_patch(n_dom_start:n_dom)
  TYPE(t_nh_state), TARGET, INTENT(INOUT)    ::  p_nh_state(n_dom)
  TYPE(t_int_state),   TARGET, INTENT(IN)    ::  p_int_state(n_dom_start:n_dom)
  TYPE(t_gridref_state), TARGET, INTENT(IN)  ::  p_grf_state(n_dom_start:n_dom)

  INTEGER, INTENT(IN) :: jg   ! child grid level
  INTEGER, INTENT(IN) :: jgp  ! parent grid level

  ! Switch if feedback is done for tracers
  ! (when calling transport and microphysics not every dynamics time step, tracer feedback
  !  should probably be restricted to transport time steps)
  LOGICAL, INTENT(IN) :: l_trac_fbk

  ! local variables

  TYPE(t_nh_prog),    POINTER     :: p_parent_prog => NULL()
  TYPE(t_nh_prog),    POINTER     :: p_child_prog => NULL()
  TYPE(t_nh_prog),    POINTER     :: p_parent_prog_rcf => NULL()
  TYPE(t_nh_prog),    POINTER     :: p_child_prog_rcf  => NULL()
  TYPE(t_grid_cells), POINTER     :: p_gcp => NULL()
  TYPE(t_grid_cells), POINTER     :: p_gcc => NULL()
  TYPE(t_grid_edges), POINTER     :: p_gep => NULL()
  TYPE(t_int_state),  POINTER     :: p_int => NULL()
  TYPE(t_gridref_state), POINTER  :: p_grf => NULL()
  TYPE(t_patch),      POINTER     :: p_pp => NULL()
  TYPE(t_patch),      POINTER     :: p_pc => NULL()

  ! Indices
  INTEGER :: jb, jc, jk, js, jt, je, jv, jgc, i_nchdom, i_chidx,    &
             i_startblk, i_endblk, i_startidx, i_endidx, ic, i_ncd, &
             i_rlstart_c, i_rlstart_e, i_rlend_c, i_rlend_e

  INTEGER :: nlev_c, nlevp1_c  ! number of full and half levels (child dom)
  INTEGER :: nlev_p, nlevp1_p  ! number of full and half levels (parent dom)
  INTEGER :: nshift, nshift_c
  INTEGER :: nblks_cp, nblks_cc, nblks_ep, nblks_ec

  REAL(wp), ALLOCATABLE, DIMENSION(:,:,:), TARGET :: feedback_rho, feedback_thv, &
                                                     feedback_vn, feedback_w
  REAL(wp), ALLOCATABLE, TARGET :: feedback_rhoqx(:,:,:,:)

  ! Note: as w(nlevp1) is diagnostic, it is excluded from feedback
  REAL(wp), DIMENSION(nproma,p_patch(jgp)%nlev,p_patch(jgp)%nblks_c), TARGET :: &
    parent_rho, diff_rho, parent_thv, diff_thv, parent_w, diff_w
  REAL(wp),DIMENSION(nproma,p_patch(jgp)%nlev,p_patch(jgp)%nblks_c,ntracer), TARGET :: &
    parent_rhoqx
  REAL(wp), DIMENSION(nproma,p_patch(jgp)%nlev,p_patch(jgp)%nblks_e), TARGET :: &
    parent_vn, diff_vn

  REAL(wp) :: rot_diff_vn(nproma,p_patch(jgp)%nlev,p_patch(jgp)%verts%start_blk(-1,1):&
    p_patch(jgp)%nblks_v)

  REAL(wp), DIMENSION(nproma,p_patch(jgp)%nlev,p_patch(jgp)%cells%start_blk(-1,1):&
    p_patch(jgp)%nblks_c) ::                                                      &
    div_diff_vn, diff_mass, parent_trmass, fbk_trmass

  REAL(wp) :: theta_v_pr(nproma,p_patch(jg)%nlev,p_patch(jg)%nblks_c)

  REAL(wp) :: tendency_corr(p_patch(jgp)%nlev), tracer_corr(p_patch(jgp)%nlev), &
              aux_diff(2*p_patch(jgp)%nlev), z_fbk_rho(nproma,4,p_patch(jg)%nlev)

  REAL(wp) :: rd_o_cvd, rd_o_p0ref, relfac, dcoef_vec


  INTEGER,  DIMENSION(:,:,:), POINTER :: iccidx, iccblk, iceidx, iceblk, iveidx, iveblk, &
                                         iecidx, iecblk, ievidx, ievblk, icnidx, icnblk
  REAL(wp), DIMENSION(:,:,:), POINTER :: p_fbk_rho, p_fbk_thv, p_fbk_w, p_fbk_vn
  REAL(wp), DIMENSION(:,:,:,:), POINTER :: p_fbk_rhoqx
  REAL(wp), DIMENSION(:,:,:), POINTER :: p_fbkwgt, p_fbkwgt_e
  REAL(wp), DIMENSION(:,:),   POINTER :: p_fbarea

  LOGICAL :: l_parallel
!-----------------------------------------------------------------------

IF (msg_level >= 10) THEN
  WRITE(message_text,'(a,i2,a,i2)') '========= Feedback:',jg,' =>',jgp
  CALL message(TRIM(routine),message_text)
ENDIF

IF ( my_process_is_mpi_seq() ) THEN
  l_parallel = .FALSE.
ELSE
  l_parallel = .TRUE.
ENDIF


p_parent_prog    => p_nh_state(jgp)%prog(nnew(jgp))
p_parent_prog_rcf=> p_nh_state(jgp)%prog(nnew_rcf(jgp))
p_child_prog     => p_nh_state(jg)%prog(nnow(jg))
p_child_prog_rcf => p_nh_state(jg)%prog(nnow_rcf(jg))
p_gcc            => p_patch(jg)%cells
p_pc             => p_patch(jg)
p_int            => p_int_state(jgp)

IF(l_parallel) THEN
  p_grf => p_grf_state_local_parent(jg)
  p_gcp => p_patch_local_parent(jg)%cells
  p_gep => p_patch_local_parent(jg)%edges
  p_pp  => p_patch_local_parent(jg)
ELSE
  p_grf => p_grf_state(jgp)
  p_gcp => p_patch(jgp)%cells
  p_gep => p_patch(jgp)%edges
  p_pp  => p_patch(jgp)
ENDIF

nlev_c   = p_pc%nlev
nlevp1_c = p_pc%nlevp1
nblks_cc = p_pc%nblks_c
nblks_ec = p_pc%nblks_e

nlev_p   = p_pp%nlev
nlevp1_p = p_pp%nlevp1
nblks_cp = p_pp%nblks_c
nblks_ep = p_pp%nblks_e

nshift = p_pc%nshift
js     = nshift

i_nchdom = MAX(1,p_pc%n_childdom)
i_chidx  = p_pc%parent_child_index

! start and end index levels for application of feedback relaxation
i_rlstart_c = grf_fbk_start_c   ! grf_fbk_start_c or min_rlcell_int
i_rlstart_e = grf_fbk_start_e   ! grf_fbk_start_e or min_rledge_int
i_rlend_c   = min_rlcell_int
i_rlend_e   = min_rledge_int

! R/c_v (not present in physical constants)
rd_o_cvd = 1._wp / cvd_o_rd

! R / p0ref
rd_o_p0ref = rd / p0ref

relfac = 0.075_wp
dcoef_vec  = 1._wp/12._wp


! Allocation of storage fields
! In parallel runs the lower bound of the feedback_* arrays must be 1 for use in exchange data,
! the lower bound of the fbk_* arrays must be i_startblk for the use in global sum.

i_startblk = p_gcp%start_blk(grf_fbk_start_c,i_chidx)
i_endblk   = p_gcp%end_blk(min_rlcell_int,i_chidx)

IF(l_parallel) i_startblk = 1

ALLOCATE(feedback_thv       (nproma, nlev_p, i_startblk:i_endblk),   &
         feedback_rho       (nproma, nlev_p, i_startblk:i_endblk),   &
         feedback_w         (nproma, nlev_p, i_startblk:i_endblk))

IF(ltransport .AND. l_trac_fbk) &
  ALLOCATE(feedback_rhoqx(nproma, nlev_p, i_startblk:i_endblk, ntracer))

i_startblk = p_gep%start_blk(grf_fbk_start_e,i_chidx)
i_endblk   = p_gep%end_blk(min_rledge,i_chidx)

IF(l_parallel) i_startblk = 1

ALLOCATE(feedback_vn(nproma, nlev_p, i_startblk:i_endblk))


! Set pointers to index and coefficient fields for cell-based variables
iccidx => p_gcp%child_idx
iccblk => p_gcp%child_blk

iceidx => p_gep%child_idx
iceblk => p_gep%child_blk



p_fbkwgt    => p_grf%fbk_wgt_c
p_fbkwgt_e  => p_grf%fbk_wgt_e
p_fbarea    => p_gcp%area

IF (p_test_run) THEN
  diff_rho = 0._wp
  diff_thv = 0._wp
  diff_w   = 0._wp
  diff_vn  = 0._wp
ENDIF

! 1. Feedback of child-domain variables to the parent grid

!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)

! Compute perturbation theta_v

i_startblk = p_gcc%start_blk(grf_bdywidth_c+1,1)
i_endblk   = p_gcc%end_blk(min_rlcell,i_nchdom)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk)
DO jb = i_startblk, i_endblk

  CALL get_indices_c(p_pc, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, grf_bdywidth_c+1, min_rlcell)

  DO jk = 1, nlev_c
    DO jc = i_startidx, i_endidx

      theta_v_pr(jc,jk,jb) = p_child_prog%theta_v(jc,jk,jb) - &
        p_nh_state(jg)%metrics%theta_ref_mc(jc,jk,jb) 

      ENDDO
    ENDDO

  ENDDO
!$OMP END DO

i_startblk = p_gcp%start_blk(grf_fbk_start_c,i_chidx)
i_endblk   = p_gcp%end_blk(min_rlcell_int,i_chidx)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk,jt,z_fbk_rho)
DO jb = i_startblk, i_endblk

  CALL get_indices_c(p_pp, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, grf_fbk_start_c, min_rlcell_int, i_chidx)

#ifdef __LOOP_EXCHANGE
  DO jc = i_startidx, i_endidx
    DO jk = 1, nlev_c
#else
!CDIR UNROLL=4
  DO jk = 1, nlev_c
    DO jc = i_startidx, i_endidx
#endif

      z_fbk_rho(jc,1,jk) =                                                   & 
        p_child_prog%rho(iccidx(jc,jb,1),jk,iccblk(jc,jb,1))*p_fbkwgt(jc,jb,1)
      z_fbk_rho(jc,2,jk) =                                                   &
        p_child_prog%rho(iccidx(jc,jb,2),jk,iccblk(jc,jb,2))*p_fbkwgt(jc,jb,2)
      z_fbk_rho(jc,3,jk) =                                                   &
        p_child_prog%rho(iccidx(jc,jb,3),jk,iccblk(jc,jb,3))*p_fbkwgt(jc,jb,3)
      z_fbk_rho(jc,4,jk) =                                                   &
        p_child_prog%rho(iccidx(jc,jb,4),jk,iccblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

      feedback_rho(jc,jk+js,jb) = SUM(z_fbk_rho(jc,1:4,jk)) - &
        p_nh_state(jg)%metrics%rho_ref_corr(jc,jk,jb) 

      feedback_thv(jc,jk+js,jb) =                                          &
        theta_v_pr(iccidx(jc,jb,1),jk,iccblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
        theta_v_pr(iccidx(jc,jb,2),jk,iccblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
        theta_v_pr(iccidx(jc,jb,3),jk,iccblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
        theta_v_pr(iccidx(jc,jb,4),jk,iccblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

      feedback_w(jc,jk+js,jb) =                                                &
        p_child_prog%w(iccidx(jc,jb,1),jk,iccblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
        p_child_prog%w(iccidx(jc,jb,2),jk,iccblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
        p_child_prog%w(iccidx(jc,jb,3),jk,iccblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
        p_child_prog%w(iccidx(jc,jb,4),jk,iccblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

    ENDDO
  ENDDO


  IF (ltransport .AND. l_trac_fbk) THEN ! tracer mass feedback
#ifdef __LOOP_EXCHANGE
    dO jc = i_startidx, i_endidx
      DO jt = 1, ntracer
        DO jk = 1, nlev_c
#else
    DO jt = 1, ntracer
!CDIR UNROLL=8
      DO jk = 1, nlev_c
        DO jc = i_startidx, i_endidx
#endif
          feedback_rhoqx(jc,jk+js,jb,jt) =                                    &
            z_fbk_rho(jc,1,jk) *                                              &
            p_child_prog_rcf%tracer(iccidx(jc,jb,1),jk,iccblk(jc,jb,1),jt) +  &
            z_fbk_rho(jc,2,jk) *                                              &
            p_child_prog_rcf%tracer(iccidx(jc,jb,2),jk,iccblk(jc,jb,2),jt) +  &
            z_fbk_rho(jc,3,jk) *                                              &
            p_child_prog_rcf%tracer(iccidx(jc,jb,3),jk,iccblk(jc,jb,3),jt) +  &
            z_fbk_rho(jc,4,jk) *                                              &
            p_child_prog_rcf%tracer(iccidx(jc,jb,4),jk,iccblk(jc,jb,4),jt)
        ENDDO
      ENDDO
    ENDDO
  ENDIF

ENDDO
!$OMP END DO

! Velocity feedback

  i_startblk = p_gep%start_blk(grf_fbk_start_e,i_chidx)
  i_endblk   = p_gep%end_blk(min_rledge_int,i_chidx)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,je,jk)
  DO jb = i_startblk, i_endblk

    CALL get_indices_e(p_pp, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, grf_fbk_start_e, min_rledge_int, i_chidx)

 !   IF (l_parallel) THEN
 !     feedback_vn(:,1:nshift,jb) = 0._wp
 !   ENDIF

#ifdef __LOOP_EXCHANGE
    DO je = i_startidx, i_endidx
      DO jk = 1, nlev_c
#else
!CDIR UNROLL=8
    DO jk = 1, nlev_c
      DO je = i_startidx, i_endidx
#endif

        feedback_vn(je,jk+js,jb) =                                                  &
          p_child_prog%vn(iceidx(je,jb,1),jk,iceblk(je,jb,1))*p_fbkwgt_e(je,jb,1) + &
          p_child_prog%vn(iceidx(je,jb,2),jk,iceblk(je,jb,2))*p_fbkwgt_e(je,jb,2)
      ENDDO
    ENDDO

  ENDDO
!$OMP END DO

!$OMP END PARALLEL

IF (l_parallel) THEN ! communication from local parent to parent

  CALL exchange_data_mult(p_pp%comm_pat_loc_to_glb_c_fbk, 3, 3*nlev_c, &
                          RECV1=parent_rho,     SEND1=feedback_rho,    &
                          RECV2=parent_thv,     SEND2=feedback_thv,    &
                          RECV3=parent_w,       SEND3=feedback_w,      &
                          nshift=nshift )


  CALL exchange_data_mult(p_pp%comm_pat_loc_to_glb_e_fbk, 1, nlev_c, &
                          RECV1=parent_vn, SEND1=feedback_vn,        &
                          nshift=nshift)

  IF (ltransport .AND. l_trac_fbk) &
  CALL exchange_data_mult(p_pp%comm_pat_loc_to_glb_c_fbk, ntracer, ntracer*nlev_c, &
                          RECV4D=parent_rhoqx, SEND4D=feedback_rhoqx, nshift=nshift)

  p_fbk_rho => parent_rho
  p_fbk_thv => parent_thv
  p_fbk_w   => parent_w
  p_fbk_vn  => parent_vn

  IF (ltransport) p_fbk_rhoqx => parent_rhoqx
ELSE

  p_fbk_rho => feedback_rho
  p_fbk_thv => feedback_thv
  p_fbk_w   => feedback_w
  p_fbk_vn  => feedback_vn

  IF (ltransport) p_fbk_rhoqx => feedback_rhoqx
ENDIF


! 2. Compute differences between feedback fields and corresponding parent fields
!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)

  i_startblk = p_patch(jgp)%cells%start_blk(-1,i_chidx)
  i_endblk   = p_patch(jgp)%cells%start_blk(i_rlstart_c,i_chidx)

!$OMP WORKSHARE
  diff_rho     (:,:,i_startblk:i_endblk) = 0._wp
  diff_thv     (:,:,i_startblk:i_endblk) = 0._wp
  diff_w       (:,:,i_startblk:i_endblk) = 0._wp
  diff_vn      (:,:,:)                   = 0._wp
  diff_mass    (:,:,:)                   = 0._wp 
  parent_trmass(:,:,:)                   = 0._wp
  fbk_trmass   (:,:,:)                   = 0._wp
!$OMP END WORKSHARE

  i_startblk = p_patch(jgp)%cells%start_blk(i_rlstart_c,i_chidx)
  i_endblk   = p_patch(jgp)%cells%end_blk(i_rlend_c,i_chidx)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk,jt)
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(p_patch(jgp), jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, i_rlstart_c, i_rlend_c, i_chidx)

    DO jk = nshift+1,nlev_p
      DO jc = i_startidx,i_endidx

        ! density
        diff_rho(jc,jk,jb) = p_fbk_rho(jc,jk,jb) - p_parent_prog%rho(jc,jk,jb)

        ! theta_v
        diff_thv(jc,jk,jb) = p_fbk_thv(jc,jk,jb) - p_parent_prog%theta_v(jc,jk,jb) + &
          p_nh_state(jgp)%metrics%theta_ref_mc(jc,jk,jb)

        ! w
        diff_w(jc,jk,jb) = p_fbk_w(jc,jk,jb) - p_parent_prog%w(jc,jk,jb)

        ! mass
        diff_mass(jc,jk,jb) = diff_rho(jc,jk,jb)*p_patch(jgp)%cells%area(jc,jb)*   &
         p_nh_state(jgp)%metrics%ddqz_z_full(jc,jk,jb)

      ENDDO
    ENDDO

    IF (ltransport .AND. l_trac_fbk) THEN

      DO jt = 1, ntracer
        DO jk = nshift+1,nlev_p
          DO jc = i_startidx, i_endidx

            parent_trmass(jc,jk,jb) = parent_trmass(jc,jk,jb) +                          &
              p_parent_prog%rho(jc,jk,jb)*p_parent_prog_rcf%tracer(jc,jk,jb,jt)*         &
              p_patch(jgp)%cells%area(jc,jb)*p_nh_state(jgp)%metrics%ddqz_z_full(jc,jk,jb)

            fbk_trmass(jc,jk,jb) = fbk_trmass(jc,jk,jb) +              &
              p_fbk_rhoqx(jc,jk,jb,jt)*p_patch(jgp)%cells%area(jc,jb)* &
              p_nh_state(jgp)%metrics%ddqz_z_full(jc,jk,jb)

          ENDDO
        ENDDO
      ENDDO
    ENDIF

  ENDDO
!$OMP END DO

  i_startblk = p_patch(jgp)%edges%start_blk(i_rlstart_e,i_chidx)
  i_endblk   = p_patch(jgp)%edges%end_blk(i_rlend_e,i_chidx)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,je,jk)
  DO jb = i_startblk, i_endblk

    CALL get_indices_e(p_patch(jgp), jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, i_rlstart_e, i_rlend_e, i_chidx)

    DO jk = nshift+1,nlev_p
      DO je = i_startidx,i_endidx
        diff_vn(je,jk,jb) = p_fbk_vn(je,jk,jb) - p_parent_prog%vn(je,jk,jb)
      ENDDO
    ENDDO

  ENDDO
!$OMP END DO
!$OMP END PARALLEL

 ! CALL sync_patch_array_mult(SYNC_C1,p_patch(jgp),3,diff_rho,diff_thv,diff_w)
  CALL sync_patch_array(SYNC_E,p_patch(jgp),diff_vn)

! 3. Smoothing of feedback-parent differences 

  iceidx => p_patch(jgp)%cells%edge_idx
  iceblk => p_patch(jgp)%cells%edge_blk

  iveidx => p_patch(jgp)%verts%edge_idx
  iveblk => p_patch(jgp)%verts%edge_blk

  iecidx => p_patch(jgp)%edges%cell_idx
  iecblk => p_patch(jgp)%edges%cell_blk

  ievidx => p_patch(jgp)%edges%vertex_idx
  ievblk => p_patch(jgp)%edges%vertex_blk


!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)

  i_startblk = p_patch(jgp)%verts%start_blk(i_rlstart_c,i_chidx)
  i_endblk   = p_patch(jgp)%verts%end_blk(min_rlvert_int,i_chidx)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jv,jk)
  DO jb = i_startblk, i_endblk

    CALL get_indices_v(p_patch(jgp), jb, i_startblk, i_endblk,         &
      i_startidx, i_endidx, i_rlstart_c, min_rlvert_int, i_chidx)

#ifdef __LOOP_EXCHANGE
    DO jv = i_startidx, i_endidx
      DO jk = nshift+1,nlev_p
#else
!CDIR UNROLL=6
    DO jk = nshift+1,nlev_p
      DO jv = i_startidx, i_endidx
#endif

        rot_diff_vn(jv,jk,jb) =   &
          diff_vn(iveidx(jv,jb,1),jk,iveblk(jv,jb,1))*p_int%geofac_rot(jv,1,jb) + &
          diff_vn(iveidx(jv,jb,2),jk,iveblk(jv,jb,2))*p_int%geofac_rot(jv,2,jb) + &
          diff_vn(iveidx(jv,jb,3),jk,iveblk(jv,jb,3))*p_int%geofac_rot(jv,3,jb) + &
          diff_vn(iveidx(jv,jb,4),jk,iveblk(jv,jb,4))*p_int%geofac_rot(jv,4,jb) + &
          diff_vn(iveidx(jv,jb,5),jk,iveblk(jv,jb,5))*p_int%geofac_rot(jv,5,jb) + &
          diff_vn(iveidx(jv,jb,6),jk,iveblk(jv,jb,6))*p_int%geofac_rot(jv,6,jb)

      ENDDO
    ENDDO
  ENDDO
!$OMP END DO

  i_startblk = p_patch(jgp)%cells%start_blk(i_rlstart_c+1,i_chidx)
  i_endblk   = p_patch(jgp)%cells%end_blk(min_rlcell_int-1,i_chidx)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk)
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(p_patch(jgp), jb, i_startblk, i_endblk,  &
      i_startidx, i_endidx, i_rlstart_c+1, min_rlcell_int-1, i_chidx)

#ifdef __LOOP_EXCHANGE
    DO jc = i_startidx, i_endidx
      DO jk = nshift+1,nlev_p
#else
!CDIR UNROLL=6
    DO jk = nshift+1,nlev_p
      DO jc = i_startidx, i_endidx
#endif

        div_diff_vn(jc,jk,jb) =   &
          diff_vn(iceidx(jc,jb,1),jk,iceblk(jc,jb,1))*p_int%geofac_div(jc,1,jb) + &
          diff_vn(iceidx(jc,jb,2),jk,iceblk(jc,jb,2))*p_int%geofac_div(jc,2,jb) + &
          diff_vn(iceidx(jc,jb,3),jk,iceblk(jc,jb,3))*p_int%geofac_div(jc,3,jb)

      ENDDO
    ENDDO
  ENDDO
!$OMP END DO

  i_startblk = p_patch(jgp)%edges%start_blk(i_rlstart_e,i_chidx)
  i_endblk   = p_patch(jgp)%edges%end_blk(i_rlend_e,i_chidx)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,je,jk)
  DO jb = i_startblk, i_endblk

    CALL get_indices_e(p_patch(jgp), jb, i_startblk, i_endblk,         &
      i_startidx, i_endidx, i_rlstart_e, i_rlend_e, i_chidx)

#ifdef __LOOP_EXCHANGE
    DO je = i_startidx, i_endidx
      DO jk = nshift+1,nlev_p
#else
!CDIR UNROLL=3
    DO jk = nshift+1,nlev_p
      DO je = i_startidx, i_endidx
#endif

        diff_vn(je,jk,jb) = diff_vn(je,jk,jb) +                &
          dcoef_vec * p_patch(jgp)%edges%area_edge(je,jb) *    &
          ( p_patch(jgp)%edges%system_orientation(je,jb) *     &
          ( rot_diff_vn(ievidx(je,jb,2),jk,ievblk(je,jb,2))    &
          - rot_diff_vn(ievidx(je,jb,1),jk,ievblk(je,jb,1)) )  &
          * p_patch(jgp)%edges%inv_primal_edge_length(je,jb) + &
          ( div_diff_vn(iecidx(je,jb,2),jk,iecblk(je,jb,2))    &
          - div_diff_vn(iecidx(je,jb,1),jk,iecblk(je,jb,1)) )  &
          * p_patch(jgp)%edges%inv_dual_edge_length(je,jb)     )

      END DO
    END DO

  END DO
!$OMP END DO

!$OMP WORKSHARE
div_diff_vn(:,:,:) = 0._wp
!$OMP END WORKSHARE

!$OMP END PARALLEL



IF ( .NOT. (ltransport .AND. l_trac_fbk)) THEN
  ! compute conservation correction for global mass only
  aux_diff(1:nlev_p) = global_sum_array3(1,.FALSE.,diff_mass)
  DO jk = 1, nlev_p
    tendency_corr(jk) = aux_diff(jk) / p_nh_state(jg)%metrics%fbk_dom_volume(jk)
  ENDDO
ELSE
  ! compute conservation correction for global mass and total tracer mass
  ! the correction is additive for density and multiplicative for tracer mass
  aux_diff(1:2*nlev_p) = global_sum_array3(2,.TRUE.,diff_mass,div_diff_vn,     &
                                         f3din2=parent_trmass,f3dd2=fbk_trmass,&
                                         diffmask=(/1,2/))
  DO jk = 1, nlev_p
    tendency_corr(jk) = aux_diff(jk) / p_nh_state(jg)%metrics%fbk_dom_volume(jk)
    tracer_corr(jk)   = aux_diff(nlev_p+jk)
  ENDDO
ENDIF

!$OMP PARALLEL PRIVATE(i_startblk,i_endblk,nshift_c,i_ncd,ic,jgc)

IF (l_masscorr_nest) THEN
  ! Add mass conservation correction to child domain in order to prevent
  ! possible inconsistencies in the mass fields; without l_nest_rcf, the
  ! correction increments are smaller, so that this correction does not appear
  ! to be needed.
  i_startblk = p_gcc%start_blk(grf_bdywidth_c+1,1)
  i_endblk   = p_gcc%end_blk(min_rlcell,i_nchdom)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk)
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(p_pc, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, grf_bdywidth_c+1, min_rlcell)

    DO jk = 1, nlev_c
      DO jc = i_startidx, i_endidx

        p_child_prog%rho(jc,jk,jb) = p_child_prog%rho(jc,jk,jb) - tendency_corr(jk)

        p_child_prog%rhotheta_v(jc,jk,jb) = p_child_prog%rho(jc,jk,jb) * &
                                            p_child_prog%theta_v(jc,jk,jb)

        p_child_prog%exner(jc,jk,jb) = &
          EXP(rd_o_cvd*LOG(rd_o_p0ref*p_child_prog%rhotheta_v(jc,jk,jb)))

      ENDDO
    ENDDO

  ENDDO
!$OMP END DO

  ! The conservation correction also needs to be applied to all nested domains
  IF (p_pc%n_childdom > 0) THEN

    DO ic = 1, p_pc%n_chd_total
      jgc = p_pc%child_id_list(ic)

      i_ncd      = MAX(1,p_patch(jgc)%n_childdom)
      i_startblk = p_patch(jgc)%cells%start_blk(grf_bdywidth_c+1,1)
      i_endblk   = p_patch(jgc)%cells%end_blk(min_rlcell,i_ncd)
      nshift_c   = p_patch(jgc)%nshift_total - p_pc%nshift_total

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk)
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(p_patch(jgc), jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, grf_bdywidth_c+1, min_rlcell)

        DO jk = 1, p_patch(jgc)%nlev
          DO jc = i_startidx, i_endidx

            p_nh_state(jgc)%prog(nnow(jgc))%rho(jc,jk,jb) = &
              p_nh_state(jgc)%prog(nnow(jgc))%rho(jc,jk,jb) + tendency_corr(jk+nshift_c)

            p_nh_state(jgc)%prog(nnow(jgc))%rhotheta_v(jc,jk,jb) = &
              p_nh_state(jgc)%prog(nnow(jgc))%rho(jc,jk,jb) * &
              p_nh_state(jgc)%prog(nnow(jgc))%theta_v(jc,jk,jb)

            p_nh_state(jgc)%prog(nnow(jgc))%exner(jc,jk,jb) = &
              EXP(rd_o_cvd*LOG(rd_o_p0ref*p_nh_state(jgc)%prog(nnow(jgc))%rhotheta_v(jc,jk,jb)))

          ENDDO
        ENDDO

      ENDDO
!$OMP END DO

    ENDDO
  ENDIF
ENDIF

  ! Execute relaxation

  i_startblk = p_patch(jgp)%cells%start_blk(i_rlstart_c,i_chidx)
  i_endblk   = p_patch(jgp)%cells%end_blk(i_rlend_c,i_chidx)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk,jt)
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(p_patch(jgp), jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, i_rlstart_c, i_rlend_c, i_chidx)

    DO jk = nshift+1,nlev_p
      DO jc = i_startidx,i_endidx

        ! density
        p_parent_prog%rho(jc,jk,jb) = p_parent_prog%rho(jc,jk,jb) + &
          relfac*diff_rho(jc,jk,jb)

        ! theta_v
        p_parent_prog%theta_v(jc,jk,jb) = p_parent_prog%theta_v(jc,jk,jb) + &
          relfac*diff_thv(jc,jk,jb)

        ! w
        p_parent_prog%w(jc,jk,jb) = p_parent_prog%w(jc,jk,jb) + &
          relfac*diff_w(jc,jk,jb)

        ! rhotheta_v is diagnosed from rho and theta_v
        p_parent_prog%rhotheta_v(jc,jk,jb) =   &
          p_parent_prog%theta_v(jc,jk,jb) * p_parent_prog%rho(jc,jk,jb)

        ! exner is diagnosed from rhotheta_v
        p_parent_prog%exner(jc,jk,jb) = EXP(rd_o_cvd*LOG(rd_o_p0ref* &
          p_parent_prog%rhotheta_v(jc,jk,jb)))

      ENDDO
    ENDDO

    IF (ltransport .AND. l_trac_fbk) THEN
      DO jt = 1, ntracer
        DO jk = nshift+1, nlev_p
          DO jc = i_startidx,i_endidx
            p_parent_prog_rcf%tracer(jc,jk,jb,jt) = p_fbk_rhoqx(jc,jk,jb,jt)* &
              tracer_corr(jk)/p_parent_prog%rho(jc,jk,jb)
          ENDDO
        ENDDO
      ENDDO
    ENDIF

  ENDDO
!$OMP END DO

  i_startblk = p_patch(jgp)%edges%start_blk(i_rlstart_e,i_chidx)
  i_endblk   = p_patch(jgp)%edges%end_blk(i_rlend_e,i_chidx)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,je)
  DO jb = i_startblk, i_endblk

    CALL get_indices_e(p_patch(jgp), jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, i_rlstart_e, i_rlend_e, i_chidx)

    DO jk = nshift+1,nlev_p
      DO je = i_startidx,i_endidx
        p_parent_prog%vn(je,jk,jb) = p_parent_prog%vn(je,jk,jb) + &
          relfac*diff_vn(je,jk,jb)
      ENDDO
    ENDDO

  ENDDO
!$OMP END DO 
!$OMP END PARALLEL

  CALL sync_patch_array_mult(SYNC_C,p_patch(jgp),3,p_parent_prog%rho,p_parent_prog%theta_v,   &
                             p_parent_prog%w)
  CALL sync_patch_array(SYNC_E,p_patch(jgp),p_parent_prog%vn)

  IF (ltransport .AND. l_trac_fbk) THEN
     CALL sync_patch_array_mult(SYNC_C, p_patch(jgp), ntracer, F4DIN=p_parent_prog_rcf%tracer,&
                                lpart4d=.TRUE.)
  ENDIF

IF (l_parallel) THEN ! Recompute rhotheta and exner on the halo points after sync of theta

#ifdef __LOOP_EXCHANGE
  DO ic = 1, p_nh_state(jgp)%metrics%ovlp_halo_c_dim(i_chidx)
    jc = p_nh_state(jgp)%metrics%ovlp_halo_c_idx(ic,i_chidx)
    jb = p_nh_state(jgp)%metrics%ovlp_halo_c_blk(ic,i_chidx)
    DO jk = nshift+1, nlev_p
#else
  DO jk = nshift+1, nlev_p
!CDIR NODEP,VOVERTAKE,VOB
    DO ic = 1, p_nh_state(jgp)%metrics%ovlp_halo_c_dim(i_chidx)
      jc = p_nh_state(jgp)%metrics%ovlp_halo_c_idx(ic,i_chidx)
      jb = p_nh_state(jgp)%metrics%ovlp_halo_c_blk(ic,i_chidx)
#endif

      p_parent_prog%rhotheta_v(jc,jk,jb) =   &
        p_parent_prog%theta_v(jc,jk,jb) * p_parent_prog%rho(jc,jk,jb)

      ! exner is diagnosed from rhotheta_v
      p_parent_prog%exner(jc,jk,jb) =   &
        EXP(rd_o_cvd*LOG(rd_o_p0ref*p_parent_prog%rhotheta_v(jc,jk,jb)))

    ENDDO
  ENDDO

ENDIF


DEALLOCATE(feedback_thv,feedback_rho,feedback_w,feedback_vn)
IF (ltransport .AND. l_trac_fbk) DEALLOCATE(feedback_rhoqx)

END SUBROUTINE relax_feedback

END MODULE mo_nh_feedback
