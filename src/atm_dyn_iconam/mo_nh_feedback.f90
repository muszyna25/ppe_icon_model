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
USE mo_model_domain,        ONLY: t_patch, t_grid_cells, t_grid_edges
USE mo_model_domain_import, ONLY: n_dom, n_dom_start
USE mo_interpolation,       ONLY: t_int_state, rbf_vec_interpol_vertex
USE mo_grf_interpolation,   ONLY: t_gridref_state, grf_velfbk
USE mo_nonhydrostatic_nml,  ONLY: l_nest_rcf, l_masscorr_nest
USE mo_dynamics_config,     ONLY: dynamics_config
USE mo_parallel_configuration,  ONLY: nproma
USE mo_run_config,          ONLY: ltransport, iforcing, msg_level, ntracer
USE mo_nonhydro_state,      ONLY: t_nh_state, t_nh_prog, t_nh_diag
USE mo_impl_constants,      ONLY: min_rlcell, min_rledge, min_rlcell_int, min_rledge_int, &
            &                     MAX_CHAR_LENGTH
USE mo_loopindices,         ONLY: get_indices_c, get_indices_e
USE mo_impl_constants_grf,  ONLY: grf_fbk_start_c, grf_fbk_start_e,          &
                                  grf_bdywidth_c
USE mo_mpi,                 ONLY:  my_process_is_mpi_seq
USE mo_communication,       ONLY: exchange_data, exchange_data_mult
USE mo_sync,                ONLY: SYNC_C, SYNC_E, sync_patch_array, &
                                  global_sum_array3, sync_patch_array_mult
USE mo_physical_constants,  ONLY: rd, cvd_o_rd, p0ref

USE mo_subdivision,         ONLY: p_patch_local_parent, p_int_state_local_parent, &
                                  p_grf_state_local_parent

IMPLICIT NONE

PRIVATE

CHARACTER(len=*), PARAMETER :: version = '$Id$'

PUBLIC :: feedback

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
INTEGER :: jb, jc, jk, jks, jt, je, jgc, i_nchdom, i_chidx, &
           i_startblk, i_endblk, i_startidx, i_endidx, ic

INTEGER :: nlev_c, nlevp1_c  ! number of full and half levels (child dom)
INTEGER :: nlev_p, nlevp1_p  ! number of full and half levels (parent dom)
INTEGER :: nshift, nshift_c

REAL(wp), DIMENSION(nproma,p_patch(jg)%nlev,p_patch(jg)%nblks_v) :: z_u, z_v
REAL(wp) ::   &  ! RBF-reconstructed velocity
  &  vn_aux(nproma,p_patch(jg)%nlev,p_patch(jg)%nblks_e,2)
REAL(wp), ALLOCATABLE :: feedback_thv_tend(:,:,:)
REAL(wp), ALLOCATABLE :: feedback_rho_tend(:,:,:)
REAL(wp), ALLOCATABLE :: feedback_vn_tend(:,:,:)
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

INTEGER :: nnow(n_dom), nnow_rcf(n_dom)
INTEGER :: nnew(n_dom), nnew_rcf(n_dom)
INTEGER :: nsav1(n_dom), nsav2(n_dom)

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

nnow    (:) = dynamics_config(1:n_dom)%nnow
nnow_rcf(:) = dynamics_config(1:n_dom)%nnow_rcf
nnew    (:) = dynamics_config(1:n_dom)%nnew
nnew_rcf(:) = dynamics_config(1:n_dom)%nnew_rcf
nsav1   (:) = dynamics_config(1:n_dom)%nsav1
nsav2   (:) = dynamics_config(1:n_dom)%nsav2

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

ALLOCATE(feedback_vn_tend(nproma, nlev_p, i_startblk:i_endblk))


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

! Part 2: Velocity components
i_startblk = p_gec%start_blk(4,1)
i_endblk   = p_gec%end_blk(min_rledge_int-1,i_nchdom)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,je,jk)
DO jb = i_startblk, i_endblk

  CALL get_indices_e(p_pc, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, 4, min_rledge_int-1)

  DO jk = 1, nlev_c
    DO je = i_startidx, i_endidx
      p_child_tend%grf_tend_vn(je,jk,jb) = &
        p_child_prog%vn(je,jk,jb) - p_child_save%vn(je,jk,jb)
    ENDDO
  ENDDO

ENDDO
!$OMP END DO

! Compute feedback tendency for density and tracers, including a layer-wise conservation correction

! Start/End block in the parent domain
i_startblk = p_gcp%start_blk(grf_fbk_start_c,i_chidx)
i_endblk   = p_gcp%end_blk(min_rlcell_int,i_chidx)


!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk,jks,jt)
DO jb = i_startblk, i_endblk

  CALL get_indices_c(p_pp, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, grf_fbk_start_c, min_rlcell_int, i_chidx)

  fbk_tend(:,:,jb) = 0._wp

#ifdef __LOOP_EXCHANGE
  DO jc = i_startidx, i_endidx
    DO jk = 1, nlev_c
      jks = jk + nshift
#else
!CDIR UNROLL=8
  DO jk = 1, nlev_c
    jks = jk + nshift
    DO jc = i_startidx, i_endidx
#endif

      feedback_rho_tend(jc,jks,jb) =                                             &
        p_child_tend%grf_tend_rho(iidx(jc,jb,1),jk,iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
        p_child_tend%grf_tend_rho(iidx(jc,jb,2),jk,iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
        p_child_tend%grf_tend_rho(iidx(jc,jb,3),jk,iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
        p_child_tend%grf_tend_rho(iidx(jc,jb,4),jk,iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

      fbk_tend(jc,jks,jb) = feedback_rho_tend(jc,jks,jb)*p_fbarea(jc,jb)* &
                           p_fb_layer_thickness(jc,jks,jb)

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
        jks = jk + nshift
#else
!CDIR UNROLL=8
      DO jk = 1, nlev_c
        jks = jk + nshift
        DO jc = i_startidx, i_endidx
#endif

          feedback_tracer_mass(jc,jks,jb,jt) =                                          &
          p_fbkwgt_tr(jc,jb,1)*p_child_prog%rho(iidx(jc,jb,1),jk,iblk(jc,jb,1))*        &
          p_child_prog_rcf%tracer(iidx(jc,jb,1),jk,iblk(jc,jb,1),jt) +                  &
          p_fbkwgt_tr(jc,jb,2)*p_child_prog%rho(iidx(jc,jb,2),jk,iblk(jc,jb,2))*        &
          p_child_prog_rcf%tracer(iidx(jc,jb,2),jk,iblk(jc,jb,2),jt) +                  &
          p_fbkwgt_tr(jc,jb,3)*p_child_prog%rho(iidx(jc,jb,3),jk,iblk(jc,jb,3))*        &
          p_child_prog_rcf%tracer(iidx(jc,jb,3),jk,iblk(jc,jb,3),jt) +                  &
          p_fbkwgt_tr(jc,jb,4)*p_child_prog%rho(iidx(jc,jb,4),jk,iblk(jc,jb,4))*        &
          p_child_prog_rcf%tracer(iidx(jc,jb,4),jk,iblk(jc,jb,4),jt)

          fbk_tr_mass(jc,jks,jb,jt) = feedback_tracer_mass(jc,jks,jb,jt)        &
            &                       * p_fbarea(jc,jb)*p_fb_layer_thickness(jc,jks,jb)

          fbk_tr_totmass(jc,jks,jb) = fbk_tr_totmass(jc,jks,jb)                 &
            &                       + fbk_tr_mass(jc,jks,jb,jt)

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
    WHERE(.NOT.p_patch(jgp)%cells%owner_mask(:,jb)) parent_tend(:,jk,jb) = 0._wp
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
        WHERE(.NOT.p_patch(jgp)%cells%owner_mask(:,jb))
          parent_tr_mass(:,jk,jb,jt) = 0._wp
          parent_tr_totmass(:,jk,jb) = 0._wp
        END WHERE
      ENDDO
    ENDDO
  ENDIF

ENDDO
!$OMP END DO
!$OMP END PARALLEL

! fbk_dom_volume not set in p_grf

IF ( .NOT. (ltransport .AND. l_trac_fbk)) THEN
  ! compute conservation correction for global mass only
  aux_diff(1:nlev_p) = global_sum_array3(1,.TRUE.,parent_tend,fbk_tend,diffmask=(/1/))
  DO jk = 1, nlev_p
    tendency_corr(jk) = aux_diff(jk) / p_grf_state(jgp)%fbk_dom_volume(jk,i_chidx)
  ENDDO
ELSE IF (iforcing <= 1) THEN
  ! compute conservation correction for global mass and each tracer separately
  ! the correction is additive for density and multiplicative for tracer masses
  aux_diff = global_sum_array3(ntracer+1,.TRUE.,parent_tend,fbk_tend,f4din=parent_tr_mass,&
                               f4dd=fbk_tr_mass,diffmask=(/1,(2,jt=1,ntracer)/))
  DO jk = 1, nlev_p
    tendency_corr(jk) = aux_diff(jk) / p_grf_state(jgp)%fbk_dom_volume(jk,i_chidx)
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
    tendency_corr(jk) = aux_diff(jk) / p_grf_state(jgp)%fbk_dom_volume(jk,i_chidx)
    tracer_corr(jk)   = aux_diff(nlev_p+jk)
  ENDDO
ENDIF

!$OMP PARALLEL PRIVATE(i_startblk,i_endblk,nshift_c)

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

    DO ic = 1, p_pc%n_childdom
      jgc = p_pc%child_id(ic)

      i_startblk = p_patch(jgc)%cells%start_blk(grf_bdywidth_c+1,1)
      i_endblk   = p_patch(jgc)%cells%end_blk(min_rlcell,i_nchdom)
      nshift_c   = p_patch(jgc)%nshift

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

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk,jks)
DO jb = i_startblk, i_endblk

  CALL get_indices_c(p_pp, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, grf_fbk_start_c, min_rlcell_int, i_chidx)

#ifdef __LOOP_EXCHANGE
  DO jc = i_startidx, i_endidx
    DO jk = 1, nlev_c
      jks = jk + nshift
#else
!CDIR UNROLL=4
  DO jk = 1, nlev_c
    jks = jk + nshift
    DO jc = i_startidx, i_endidx
#endif

      feedback_rho_tend(jc,jks,jb) = feedback_rho_tend(jc,jks,jb) + tendency_corr(jks)

      feedback_thv_tend(jc,jks,jb) =                                             &
        p_child_tend%grf_tend_thv(iidx(jc,jb,1),jk,iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
        p_child_tend%grf_tend_thv(iidx(jc,jb,2),jk,iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
        p_child_tend%grf_tend_thv(iidx(jc,jb,3),jk,iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
        p_child_tend%grf_tend_thv(iidx(jc,jb,4),jk,iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

      feedback_w_tend(jc,jks,jb) =                                             &
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
  CALL rbf_vec_interpol_vertex( p_child_tend%grf_tend_vn, p_pc,  &
                                p_intc, z_u, z_v)
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

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,je,jk,jks)
  DO jb = i_startblk, i_endblk

    CALL get_indices_e(p_pp, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, grf_fbk_start_e, min_rledge_int, i_chidx)

    IF (l_parallel) THEN
      feedback_vn_tend(:,1:nshift,jb) = 0._wp
    ENDIF

#ifdef __LOOP_EXCHANGE
    DO je = i_startidx, i_endidx
      DO jk = 1, nlev_c
        jks = jk + nshift
#else
!CDIR UNROLL=8
    DO jk = 1, nlev_c
      jks = jk + nshift
      DO je = i_startidx, i_endidx
#endif

        feedback_vn_tend(je,jks,jb) =                                                 &
          p_child_tend%grf_tend_vn(iidx(je,jb,1),jk,iblk(je,jb,1))*p_fbkwgt(je,jb,1) + &
          p_child_tend%grf_tend_vn(iidx(je,jb,2),jk,iblk(je,jb,2))*p_fbkwgt(je,jb,2)
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
#ifdef _URD
!CDIR UNROLL=_URD
#endif
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


!$OMP DO PRIVATE(jb,i_startidx,i_endidx,je,jk,jks)
  DO jb = i_startblk, i_endblk

    CALL get_indices_e(p_pp, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, grf_fbk_start_e, min_rledge_int, i_chidx)

    IF (l_parallel) THEN
      feedback_vn_tend(:,1:nshift,jb) = 0._wp
    ENDIF

#ifdef __LOOP_EXCHANGE
    DO je = i_startidx, i_endidx
      DO jk = 1, nlev_c
        jks = jk + nshift
#else
#ifdef _URD
!CDIR UNROLL=_URD
#endif
    DO jk = 1, nlev_c
      jks = jk + nshift
      DO je = i_startidx, i_endidx
#endif

        feedback_vn_tend(je,jks,jb) =                                             &
        ( p_fbkwgt(je,jb,1)*(p_child_tend%grf_tend_vn(iidx(je,jb,1),jk,iblk(je,jb,1)) + &
                             p_child_tend%grf_tend_vn(iidx(je,jb,2),jk,iblk(je,jb,2)))+ &
          p_fbkwgt(je,jb,3)*vn_aux(iidx(je,jb,1),jk,iblk(je,jb,1),1) +                  &
          p_fbkwgt(je,jb,4)*vn_aux(iidx(je,jb,1),jk,iblk(je,jb,1),2) +                  &
          p_fbkwgt(je,jb,5)*vn_aux(iidx(je,jb,2),jk,iblk(je,jb,2),1) +                  &
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
      p_parent_save%vn(i_startidx:i_endidx,nshift+1:nlev_p,jb) + &
      feedback_vn_tend(i_startidx:i_endidx,nshift+1:nlev_p,jb)

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


  CALL exchange_data_mult(p_pp%comm_pat_loc_to_glb_e_fbk, 1, nlev_c,      &
                          RECV1=p_parent_prog%vn, SEND1=feedback_vn_tend, &
                          ADD1=p_parent_save%vn,  nshift=nshift)

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
!$OMP END PARALLEL

ENDIF ! l_parallel

CALL sync_patch_array_mult(SYNC_C,p_patch(jgp),3,p_parent_prog%rho,p_parent_prog%theta_v,   &
                           p_parent_prog%w)
CALL sync_patch_array(SYNC_E,p_patch(jgp),p_parent_prog%vn)

IF (ltransport .AND. l_trac_fbk) THEN
     CALL sync_patch_array_mult(SYNC_C, p_patch(jgp), ntracer, F4DIN=p_parent_prog_rcf%tracer,&
                                lpart4d=.TRUE.)
ENDIF

IF (l_parallel) THEN ! Recompute rhotheta and exner on the halo points after sync of theta

!$OMP PARALLEL
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
!$OMP END PARALLEL

ENDIF

DEALLOCATE(parent_tend, fbk_tend, feedback_thv_tend, feedback_rho_tend, &
           feedback_w_tend, feedback_vn_tend)
IF (ltransport .AND. l_trac_fbk) &
  DEALLOCATE(feedback_tracer_mass,parent_tr_mass,parent_tr_totmass,fbk_tr_mass,fbk_tr_totmass)

END SUBROUTINE feedback


END MODULE mo_nh_feedback
