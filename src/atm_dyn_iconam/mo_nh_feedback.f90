!>
!!  This module contains the routines needed for nesting in the nonhydrostatic.
!!  version.
!!
!! @par Revision History
!!  Developed and tested by Guenther Zaengl, DWD (2010-02-10)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!!

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_nh_feedback
  !
  !
  USE mo_kind,                ONLY: wp, vp
  USE mo_exception,           ONLY: message_text, message
  USE mo_model_domain,        ONLY: t_patch, t_grid_cells, t_grid_edges, p_patch_local_parent
  USE mo_grid_config,         ONLY: n_dom, n_dom_start
  USE mo_intp_data_strc,      ONLY: t_int_state
  USE mo_intp_rbf,            ONLY: rbf_vec_interpol_vertex
  USE mo_grf_intp_data_strc,  ONLY: t_gridref_state, p_grf_state_local_parent
  USE mo_gridref_config,      ONLY: grf_velfbk, l_mass_consvcorr, fbk_relax_timescale, grf_scalfbk, grf_tracfbk
  USE mo_nonhydrostatic_config, ONLY: l_masscorr_nest
  USE mo_dynamics_config,     ONLY: nnow, nnew, nnow_rcf, nnew_rcf, nsav1, nsav2 
  USE mo_parallel_config,     ONLY: nproma, p_test_run
  USE mo_run_config,          ONLY: ltransport, iforcing, msg_level, ntracer
  USE mo_nonhydro_types,      ONLY: t_nh_state, t_nh_prog, t_nh_diag
  USE mo_impl_constants,      ONLY: min_rlcell, min_rledge, min_rlcell_int, min_rledge_int, &
    &                     min_rlvert_int, MAX_CHAR_LENGTH, nclass_aero
  USE mo_loopindices,         ONLY: get_indices_c, get_indices_e, get_indices_v
  USE mo_impl_constants_grf,  ONLY: grf_fbk_start_c, grf_fbk_start_e,          &
    grf_bdywidth_c
  USE mo_communication,       ONLY: exchange_data_mult, exchange_data_mult_mixprec
  USE mo_sync,                ONLY: SYNC_C, SYNC_E, sync_patch_array, &
    global_sum_array3, sync_patch_array_mult
  USE mo_physical_constants,  ONLY: rd, cvd_o_rd, p0ref
  USE mo_nwp_lnd_types,       ONLY: t_lnd_state, t_lnd_prog, t_wtr_prog
  USE mo_nwp_phy_types,       ONLY: t_nwp_phy_diag
  USE mo_lnd_nwp_config,      ONLY: ntiles_total, ntiles_water, lseaice
  USE mo_atm_phy_nwp_config,  ONLY: atm_phy_nwp_config, iprog_aero

  IMPLICIT NONE

  PRIVATE


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
  SUBROUTINE feedback(p_patch, p_nh_state, p_int_state, p_grf_state, p_lnd_state, &
    jg, jgp)

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_nh_feedback:feedback'


    TYPE(t_patch),       TARGET, INTENT(IN)    ::  p_patch(n_dom_start:n_dom)
    TYPE(t_nh_state), TARGET, INTENT(INOUT)    ::  p_nh_state(n_dom)
    TYPE(t_int_state),   TARGET, INTENT(IN)    ::  p_int_state(n_dom_start:n_dom)
    TYPE(t_gridref_state), TARGET, INTENT(IN)  ::  p_grf_state(n_dom_start:n_dom)
    TYPE(t_lnd_state), TARGET, INTENT(IN)      ::  p_lnd_state(n_dom)

    INTEGER, INTENT(IN) :: jg   ! child grid level
    INTEGER, INTENT(IN) :: jgp  ! parent grid level

    ! local variables

    TYPE(t_nh_prog),    POINTER     :: p_parent_prog => NULL()
    TYPE(t_nh_prog),    POINTER     :: p_parent_save => NULL()
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
    TYPE(t_gridref_state), POINTER  :: p_grfp => NULL()
    TYPE(t_int_state), POINTER      :: p_intc => NULL()
    TYPE(t_patch),      POINTER     :: p_pp => NULL()
    TYPE(t_patch),      POINTER     :: p_pc => NULL()
    TYPE(t_lnd_prog),   POINTER     :: p_lndp => NULL()
    TYPE(t_lnd_prog),   POINTER     :: p_lndc => NULL()
    TYPE(t_wtr_prog),   POINTER     :: p_wtrp => NULL()

    ! Indices
    INTEGER :: jb, jc, jk, jt, je, js, jgc, i_nchdom, i_chidx, &
      i_startblk, i_endblk, i_startidx, i_endidx, ic, i_ncd

    INTEGER :: nlev_c, nlevp1_c  ! number of full and half levels (child dom)
    INTEGER :: nlev_p, nlevp1_p  ! number of full and half levels (parent dom)
    INTEGER :: nshift, nshift_c
    INTEGER :: ntiles, ntiles_h2o

    REAL(wp), DIMENSION(nproma,p_patch(jg)%nlev,p_patch(jg)%nblks_v) :: z_u, z_v
    REAL(wp) ::   &  ! RBF-reconstructed velocity
      &  vn_aux(nproma,p_patch(jg)%nlev,p_patch(jg)%nblks_e,2)
    REAL(wp), ALLOCATABLE :: feedback_thv_tend(:,:,:)
    REAL(wp), ALLOCATABLE :: feedback_rho_tend(:,:,:)
    REAL(wp), ALLOCATABLE :: feedback_vn(:,:,:)
    REAL(wp), ALLOCATABLE :: feedback_tg(:,:,:)
    REAL(wp), ALLOCATABLE :: feedback_w_tend(:,:,:)
    REAL(wp), ALLOCATABLE :: feedback_tracer_mass(:,:,:,:)
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)   :: fbk_tend, parent_tend
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:,:) :: fbk_tr_mass, parent_tr_mass
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)   :: fbk_tr_totmass, parent_tr_totmass
    REAL(wp) :: tendency_corr(p_patch(jgp)%nlev),        &
      aux_diff((ntracer+1)*p_patch(jgp)%nlev), &
      tracer_corr(ntracer*p_patch(jgp)%nlev),  &
      diff_tg(nproma,p_patch(jgp)%nblks_c),    &
      tg_pr(nproma,p_patch(jg)%nblks_c),       &
      parent_tg(nproma,1,p_patch(jgp)%nblks_c)

    REAL(wp) :: rd_o_cvd, rd_o_p0ref, relfac

    INTEGER, DIMENSION(:,:,:), POINTER :: iidx, iblk, iidxv, iblkv
    REAL(wp), DIMENSION(:,:,:), POINTER :: p_fbkwgt, p_fbkwgt_tr, p_fb_layer_thickness
    REAL(wp), DIMENSION(:,:),   POINTER :: p_fbarea

    !-----------------------------------------------------------------------

    IF (msg_level >= 10) THEN
      WRITE(message_text,'(a,i2,a,i2)') '========= Feedback:',jg,' =>',jgp
      CALL message(TRIM(routine),message_text)
    ENDIF


    p_parent_prog    => p_nh_state(jgp)%prog(nnew(jgp))
    p_parent_prog_rcf=> p_nh_state(jgp)%prog(nnew_rcf(jgp))
    p_parent_save    => p_nh_state(jgp)%prog(nsav1(jgp))
    p_child_prog     => p_nh_state(jg)%prog(nnow(jg))
    p_child_prog_rcf => p_nh_state(jg)%prog(nnow_rcf(jg))
    p_child_save     => p_nh_state(jg)%prog(nsav2(jg))
    p_child_tend     => p_nh_state(jg)%diag
    p_intc           => p_int_state(jg)
    p_gcc            => p_patch(jg)%cells
    p_gec            => p_patch(jg)%edges
    p_pc             => p_patch(jg)
    p_lndp           => p_lnd_state(jgp)%prog_lnd(nnew_rcf(jgp))
    p_lndc           => p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))
    p_wtrp           => p_lnd_state(jgp)%prog_wtr(nnew_rcf(jgp))

    p_grf => p_grf_state_local_parent(jg)
    p_grfp => p_grf_state(jgp)
    p_gcp => p_patch_local_parent(jg)%cells
    p_gep => p_patch_local_parent(jg)%edges
    p_pp  => p_patch_local_parent(jg)

    nlev_c   = p_pc%nlev
    nlevp1_c = p_pc%nlevp1

    nlev_p   = p_pp%nlev
    nlevp1_p = p_pp%nlevp1

    nshift = p_pc%nshift
    js     = nshift

    IF (atm_phy_nwp_config(jgp)%inwp_surface > 0 ) THEN
      ntiles     = ntiles_total
      ntiles_h2o = ntiles_water
    ELSE
      ntiles     = 0
      ntiles_h2o = 0
    ENDIF

    i_nchdom = MAX(1,p_pc%n_childdom)
    i_chidx  = p_pc%parent_child_index

    ! R/c_v (not present in physical constants)
    rd_o_cvd = 1._wp / cvd_o_rd

    ! R / p0ref
    rd_o_p0ref = rd / p0ref

    ! Relaxation factor used for ground temperature relaxation
    relfac = 0.075_wp

    ! parent_tend, parent_tr_mass, parent_tr_totmass are always calculated on the global parent
    ! and thus have to be allocated within global parent limits

    i_startblk = p_patch(jgp)%cells%start_blk(grf_fbk_start_c,i_chidx)
    i_endblk   = p_patch(jgp)%cells%end_blk(min_rlcell_int,i_chidx)

    ALLOCATE(parent_tend(nproma, nlev_p, i_startblk:i_endblk))
    IF (ltransport) &
      ALLOCATE(parent_tr_mass(nproma, nlev_p, i_startblk:i_endblk, ntracer), &
      parent_tr_totmass(nproma, nlev_p, i_startblk:i_endblk) )

    ! Allocation of storage fields
    ! In parallel runs the lower bound of the feedback_* arrays must be 1 for use in exchange data,
    ! the lower bound of the fbk_* arrays must be i_startblk for the use in global sum.

    i_startblk = p_gcp%start_blk(grf_fbk_start_c,i_chidx)
    i_endblk   = p_gcp%end_blk(min_rlcell_int,i_chidx)

    ALLOCATE(fbk_tend(nproma, nlev_p,  i_startblk:i_endblk))
    IF(ltransport) &
      ALLOCATE(fbk_tr_mass(nproma, nlev_p, i_startblk:i_endblk, ntracer), &
      fbk_tr_totmass(nproma, nlev_p, i_startblk:i_endblk)  )

    i_startblk = 1

    ALLOCATE(feedback_thv_tend  (nproma, nlev_p, i_startblk:i_endblk),   &
      feedback_rho_tend  (nproma, nlev_p, i_startblk:i_endblk),   &
      feedback_w_tend    (nproma, nlevp1_p, i_startblk:i_endblk),  &
      feedback_tg        (nproma, 1, i_startblk:i_endblk)  )

    IF(ltransport) &
      ALLOCATE(feedback_tracer_mass(nproma, nlev_p, i_startblk:i_endblk, ntracer))

    i_startblk = p_gep%start_blk(grf_fbk_start_e,i_chidx)
    i_endblk   = p_gep%end_blk(min_rledge,i_chidx)

    i_startblk = 1

    ALLOCATE(feedback_vn(nproma, nlev_p, i_startblk:i_endblk))


    ! Set pointers to index and coefficient fields for cell-based variables
    iidx => p_gcp%child_idx
    iblk => p_gcp%child_blk

    IF (grf_scalfbk == 1) THEN
      p_fbkwgt    => p_grf%fbk_wgt_aw
    ELSE
      p_fbkwgt    => p_grf%fbk_wgt_bln
    ENDIF
    IF (grf_tracfbk == 1) THEN
      p_fbkwgt_tr => p_grf%fbk_wgt_aw
    ELSE
      p_fbkwgt_tr => p_grf%fbk_wgt_bln
    ENDIF
    p_fbarea    => p_gcp%area

    p_fb_layer_thickness => p_gcp%ddqz_z_full

    ! Preparation of feedback: compute child tendencies


!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)
    !
    ! Part 1: cell-based variables
    i_startblk = p_gcc%start_blk(3,1)
    i_endblk   = p_gcc%end_blk(min_rlcell_int,i_nchdom)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk,jt) ICON_OMP_DEFAULT_SCHEDULE
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
        tg_pr(jc,jb) = p_lndc%t_g(jc,jb) - p_nh_state(jg)%metrics%tsfc_ref(jc,jb)
      ENDDO

    ENDDO
!$OMP END DO

    ! Compute feedback tendency for density and tracers, including a layer-wise conservation correction

    ! Start/End block in the parent domain
    i_startblk = p_gcp%start_blk(grf_fbk_start_c,i_chidx)
    i_endblk   = p_gcp%end_blk(min_rlcell_int,i_chidx)


!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk,jt) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_pp, jb, i_startblk, i_endblk, &
        i_startidx, i_endidx, grf_fbk_start_c, min_rlcell_int)

      fbk_tend(:,:,jb) = 0._wp

#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jk = 1, nlev_c
#else
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
      IF (ltransport) THEN

        fbk_tr_totmass(:,:,jb) = 0._wp

        DO jt = 1, ntracer

          fbk_tr_mass(:,:,jb,jt) = 0._wp

          feedback_tracer_mass(:,1:nshift,jb,jt) = 0._wp

#ifdef __LOOP_EXCHANGE
          DO jc = i_startidx, i_endidx
            DO jk = 1, nlev_c
#else
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
    i_startblk = p_patch(jgp)%cells%start_blk(1,1)
    i_endblk   = p_patch(jgp)%cells%end_blk(min_rlcell_int,i_nchdom)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      parent_tend(:,:,jb) = 0._wp

      CALL get_indices_c(p_patch(jgp), jb, i_startblk, i_endblk, &
        i_startidx, i_endidx, 1, min_rlcell_int)

      DO jk = 1, nlev_p
        DO jc = i_startidx, i_endidx
          IF (p_grfp%mask_ovlp_c(jc,jb,i_chidx)) THEN
            parent_tend(jc,jk,jb) =  &
              (p_parent_prog%rho(jc,jk,jb) - p_parent_save%rho(jc,jk,jb))*   &
              p_patch(jgp)%cells%area(jc,jb)*p_nh_state(jgp)%metrics%ddqz_z_full(jc,jk,jb)
          ENDIF
        ENDDO
      ENDDO

      IF (ltransport) THEN

        parent_tr_totmass(:,:,jb) = 0._wp

        DO jt = 1, ntracer

          parent_tr_mass(:,:,jb,jt) = 0._wp

          DO jk = 1, nlev_p
            DO jc = i_startidx, i_endidx
              IF (p_grfp%mask_ovlp_c(jc,jb,i_chidx)) THEN
                parent_tr_mass(jc,jk,jb,jt) =  &
                  p_parent_prog%rho(jc,jk,jb)*p_parent_prog_rcf%tracer(jc,jk,jb,jt)* &
                  p_patch(jgp)%cells%area(jc,jb)*p_nh_state(jgp)%metrics%ddqz_z_full(jc,jk,jb)

                parent_tr_totmass(jc,jk,jb) = parent_tr_totmass(jc,jk,jb) + &
                  parent_tr_mass(jc,jk,jb,jt)
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDIF

    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    ! fbk_dom_volume is now set in p_nh_state(jg)%metrics

    IF (l_mass_consvcorr) THEN
      IF ( .NOT. (ltransport)) THEN
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
    ELSE ! conservation correction turned off
      tendency_corr(:) = 0.0_wp
      tracer_corr(:)   = 1.0_wp
    ENDIF

!$OMP PARALLEL PRIVATE(i_startblk,i_endblk,nshift_c,i_ncd,ic,jgc)

    IF (l_masscorr_nest) THEN
      ! Add mass conservation correction to child domain in order to prevent
      ! possible inconsistencies in the mass fields
      i_startblk = p_gcc%start_blk(grf_bdywidth_c+1,1)
      i_endblk   = p_gcc%end_blk(min_rlcell,i_nchdom)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(p_pc, jb, i_startblk, i_endblk, &
          i_startidx, i_endidx, grf_bdywidth_c+1, min_rlcell)

        DO jk = 1, nlev_c
          DO jc = i_startidx, i_endidx

            p_child_prog%rho(jc,jk,jb) = p_child_prog%rho(jc,jk,jb) + tendency_corr(jk+js)

            p_child_prog%exner(jc,jk,jb) = &
              EXP(rd_o_cvd*LOG(rd_o_p0ref*p_child_prog%rho(jc,jk,jb)*p_child_prog%theta_v(jc,jk,jb)))

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
          nshift_c   = p_patch(jgc)%nshift_total - p_pp%nshift_total

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk) ICON_OMP_DEFAULT_SCHEDULE
          DO jb = i_startblk, i_endblk

            CALL get_indices_c(p_patch(jgc), jb, i_startblk, i_endblk, &
              i_startidx, i_endidx, grf_bdywidth_c+1, min_rlcell)

            DO jk = 1, p_patch(jgc)%nlev
              DO jc = i_startidx, i_endidx

                p_nh_state(jgc)%prog(nnow(jgc))%rho(jc,jk,jb) = &
                  p_nh_state(jgc)%prog(nnow(jgc))%rho(jc,jk,jb) + tendency_corr(jk+nshift_c)

                p_nh_state(jgc)%prog(nnow(jgc))%exner(jc,jk,jb) = &
                  EXP(rd_o_cvd*LOG(rd_o_p0ref*p_nh_state(jgc)%prog(nnow(jgc))%rho(jc,jk,jb) * &
                  p_nh_state(jgc)%prog(nnow(jgc))%theta_v(jc,jk,jb)))

              ENDDO
            ENDDO

          ENDDO
!$OMP END DO

        ENDDO
      ENDIF
    ENDIF

    i_startblk = p_gcp%start_blk(grf_fbk_start_c,i_chidx)
    i_endblk   = p_gcp%end_blk(min_rlcell_int,i_chidx)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_pp, jb, i_startblk, i_endblk, &
        i_startidx, i_endidx, grf_fbk_start_c, min_rlcell_int)

#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jk = 1, nlev_c
#else
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

        feedback_tg(jc,1,jb) =                                   &
          tg_pr(iidx(jc,jb,1),iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
          tg_pr(iidx(jc,jb,2),iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
          tg_pr(iidx(jc,jb,3),iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
          tg_pr(iidx(jc,jb,4),iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)
      ENDDO

    ENDDO
!$OMP END DO

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

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,je,jk) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_e(p_pp, jb, i_startblk, i_endblk, &
          i_startidx, i_endidx, grf_fbk_start_e, min_rledge_int)

        feedback_vn(:,1:nshift,jb) = 0._wp

#ifdef __LOOP_EXCHANGE
        DO je = i_startidx, i_endidx
          DO jk = 1, nlev_c
#else
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

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,je,jk) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_e(p_pc, jb, i_startblk, i_endblk, &
          i_startidx, i_endidx, 4, min_rledge)

#ifdef __LOOP_EXCHANGE
        DO je = i_startidx, i_endidx
          DO jk = 1, nlev_c
#else
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


!$OMP DO PRIVATE(jb,i_startidx,i_endidx,je,jk) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_e(p_pp, jb, i_startblk, i_endblk, &
          i_startidx, i_endidx, grf_fbk_start_e, min_rledge_int)

        feedback_vn(:,1:nshift,jb) = 0._wp

#ifdef __LOOP_EXCHANGE
        DO je = i_startidx, i_endidx
          DO jk = 1, nlev_c
#else
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
!$OMP END PARALLEL

    CALL exchange_data_mult(p_pp%comm_pat_loc_to_glb_c_fbk, 4, 3*nlev_c+2,        &
      RECV1=p_parent_prog%rho,     SEND1=feedback_rho_tend, &
      RECV2=p_parent_prog%theta_v, SEND2=feedback_thv_tend, &
      RECV3=p_parent_prog%w,       SEND3=feedback_w_tend,   &
      RECV4=parent_tg,             SEND4=feedback_tg,       &
      nshift=nshift )


    CALL exchange_data_mult(p_pp%comm_pat_loc_to_glb_e_fbk, 1, nlev_c, &
      RECV1=p_parent_prog%vn, SEND1=feedback_vn, &
      nshift=nshift)

    IF (ltransport) THEN

      CALL exchange_data_mult(p_pp%comm_pat_loc_to_glb_c_fbk, ntracer, ntracer*nlev_c, &
        &                     RECV4D=p_parent_prog_rcf%tracer,                         &
        &                     SEND4D=feedback_tracer_mass, nshift=nshift)

    ENDIF


    i_startblk = p_patch(jgp)%cells%start_blk(1,1)
    i_endblk   = p_patch(jgp)%cells%end_blk(min_rlcell_int,i_nchdom)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk,jt) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch(jgp), jb, i_startblk, i_endblk, &
        i_startidx, i_endidx, 1, min_rlcell_int)

      ! So far, the prognostic state variables contain only the feedback increment;
      ! here, the save value from the previous time step is added

      DO jk = nshift+1, nlev_p
        DO jc = i_startidx, i_endidx
          IF (p_grfp%mask_ovlp_c(jc,jb,i_chidx)) THEN
            p_parent_prog%w(jc,jk,jb) = p_parent_prog%w(jc,jk,jb) + p_parent_save%w(jc,jk,jb)

            p_parent_prog%rho(jc,jk,jb) = p_parent_prog%rho(jc,jk,jb) + p_parent_save%rho(jc,jk,jb)

            p_parent_prog%theta_v(jc,jk,jb) = p_parent_prog%theta_v(jc,jk,jb) + p_parent_save%theta_v(jc,jk,jb)

            ! exner in feedback area is diagnosed from rho*theta_v
            p_parent_prog%exner(jc,jk,jb) = EXP(rd_o_cvd*LOG(rd_o_p0ref*p_parent_prog%theta_v(jc,jk,jb)* &
              p_parent_prog%rho(jc,jk,jb)))
          ENDIF
        ENDDO
      ENDDO

      jk = nlevp1_p
      DO jc = i_startidx, i_endidx
        IF (p_grfp%mask_ovlp_c(jc,jb,i_chidx)) THEN
          ! surface level of w
          p_parent_prog%w(jc,jk,jb) = p_parent_prog%w(jc,jk,jb) + p_parent_save%w(jc,jk,jb)

          ! ground temperature - relaxation is used here because stability problems appeared otherwise
          diff_tg(jc,jb) = parent_tg(jc,1,jb) - p_lndp%t_g(jc,jb) + p_nh_state(jgp)%metrics%tsfc_ref(jc,jb)
          p_lndp%t_g(jc,jb) = p_lndp%t_g(jc,jb) + relfac*diff_tg(jc,jb)
        ENDIF
      ENDDO

      DO jt = 1, ntiles
        DO jc = i_startidx, i_endidx
          IF (p_grfp%mask_ovlp_c(jc,jb,i_chidx)) THEN
            p_lndp%t_g_t(jc,jb,jt)    = p_lndp%t_g_t(jc,jb,jt)    +  relfac*diff_tg(jc,jb)
            p_lndp%t_s_t(jc,jb,jt)    = p_lndp%t_s_t(jc,jb,jt)    +  relfac*diff_tg(jc,jb)
            p_lndp%t_so_t(jc,1,jb,jt) = p_lndp%t_so_t(jc,1,jb,jt) +  relfac*diff_tg(jc,jb)
            p_lndp%t_so_t(jc,2,jb,jt) = p_lndp%t_so_t(jc,2,jb,jt) +  relfac*diff_tg(jc,jb)
          ENDIF
        ENDDO
      ENDDO

      DO jt = ntiles+1, ntiles+ntiles_h2o
        DO jc = i_startidx, i_endidx
          IF (p_grfp%mask_ovlp_c(jc,jb,i_chidx)) THEN
            p_lndp%t_g_t(jc,jb,jt) = p_lndp%t_g_t(jc,jb,jt) + relfac*diff_tg(jc,jb)
            p_lndp%t_s_t(jc,jb,jt) = p_lndp%t_s_t(jc,jb,jt) + relfac*diff_tg(jc,jb)
          ENDIF
        ENDDO
      ENDDO
      IF (lseaice) THEN
        DO jc = i_startidx, i_endidx
          IF (p_grfp%mask_ovlp_c(jc,jb,i_chidx)) THEN
            p_wtrp%t_ice(jc,jb)     = p_wtrp%t_ice(jc,jb)     + relfac*diff_tg(jc,jb)
            p_wtrp%t_snow_si(jc,jb) = p_wtrp%t_snow_si(jc,jb) + relfac*diff_tg(jc,jb)
          ENDIF
        ENDDO
      ENDIF

      ! divide tracer density (which is the feedback quantity) by air density,
      ! and apply multiplicative mass conservation correction
      IF (ltransport) THEN
        IF (iforcing <= 1) THEN
          DO jt = 1, ntracer
            DO jk = nshift+1, nlev_p
              DO jc = i_startidx, i_endidx
                IF (p_grfp%mask_ovlp_c(jc,jb,i_chidx)) THEN
                  p_parent_prog_rcf%tracer(jc,jk,jb,jt) =                 &
                    p_parent_prog_rcf%tracer(jc,jk,jb,jt) *               &
                    tracer_corr(jk+(jt-1)*nlev_p)/p_parent_prog%rho(jc,jk,jb)
                ENDIF
              ENDDO
            ENDDO
          ENDDO
        ELSE
          DO jt = 1, ntracer
            DO jk = nshift+1, nlev_p
              DO jc = i_startidx, i_endidx
                IF (p_grfp%mask_ovlp_c(jc,jb,i_chidx)) THEN
                  p_parent_prog_rcf%tracer(jc,jk,jb,jt) =         &
                    p_parent_prog_rcf%tracer(jc,jk,jb,jt) *       &
                    tracer_corr(jk)/p_parent_prog%rho(jc,jk,jb)
                ENDIF
              ENDDO
            ENDDO
          ENDDO
        ENDIF
      ENDIF

    ENDDO
!$OMP END DO

    ! Do the same computations also on the halo points
    ! NOTE: the feedback increments have already been set on the halo points because the
    !       communication pattern is defined accordingly
!$OMP DO PRIVATE(jk,jc,jb,ic) ICON_OMP_DEFAULT_SCHEDULE
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

        p_parent_prog%w(jc,jk,jb) = p_parent_prog%w(jc,jk,jb) + p_parent_save%w(jc,jk,jb)

        p_parent_prog%rho(jc,jk,jb) = p_parent_prog%rho(jc,jk,jb) + p_parent_save%rho(jc,jk,jb)

        p_parent_prog%theta_v(jc,jk,jb) = p_parent_prog%theta_v(jc,jk,jb) + &
          p_parent_save%theta_v(jc,jk,jb)

        ! exner is diagnosed from rho*theta_v
        p_parent_prog%exner(jc,jk,jb) =   &
          EXP(rd_o_cvd*LOG(rd_o_p0ref*p_parent_prog%theta_v(jc,jk,jb)*p_parent_prog%rho(jc,jk,jb)))

      ENDDO
    ENDDO
!$OMP END DO

!$OMP DO PRIVATE(jc,jb,ic) ICON_OMP_DEFAULT_SCHEDULE
!CDIR NODEP,VOVERTAKE,VOB
    DO ic = 1, p_nh_state(jgp)%metrics%ovlp_halo_c_dim(i_chidx)
      jc = p_nh_state(jgp)%metrics%ovlp_halo_c_idx(ic,i_chidx)
      jb = p_nh_state(jgp)%metrics%ovlp_halo_c_blk(ic,i_chidx)
      p_parent_prog%w(jc,nlevp1_p,jb) = p_parent_prog%w(jc,nlevp1_p,jb) + &
        p_parent_save%w(jc,nlevp1_p,jb)
    ENDDO
!$OMP END DO

    ! Recompute tracer also on the halo points
    IF (ltransport) THEN
      IF (iforcing <= 1) THEN
!$OMP DO PRIVATE(jk,jt,jc,jb,ic) ICON_OMP_DEFAULT_SCHEDULE
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
!$OMP END DO NOWAIT
      ELSE ! iforcing > 1
!$OMP DO PRIVATE(jk,jt,jc,jb,ic) ICON_OMP_DEFAULT_SCHEDULE
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
!$OMP END DO NOWAIT
      ENDIF
    ENDIF

!$OMP END PARALLEL

    DEALLOCATE(parent_tend, fbk_tend, feedback_thv_tend, feedback_rho_tend, &
      feedback_w_tend, feedback_vn, feedback_tg)
    IF (ltransport) &
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
  SUBROUTINE relax_feedback(p_patch, p_nh_state, p_int_state, p_grf_state, prm_diag, jg, jgp, dt_fbk)

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_nh_feedback:relax_feedback'

    TYPE(t_patch),       TARGET, INTENT(IN)    ::  p_patch(n_dom_start:n_dom)
    TYPE(t_nh_state), TARGET, INTENT(INOUT)    ::  p_nh_state(n_dom)
    TYPE(t_int_state),   TARGET, INTENT(IN)    ::  p_int_state(n_dom_start:n_dom)
    TYPE(t_gridref_state), TARGET, INTENT(IN)  ::  p_grf_state(n_dom_start:n_dom)
    TYPE(t_nwp_phy_diag), TARGET, INTENT(INOUT)::  prm_diag(n_dom)

    INTEGER, INTENT(IN) :: jg   ! child grid level
    INTEGER, INTENT(IN) :: jgp  ! parent grid level

    REAL(wp), INTENT(IN) :: dt_fbk ! time step at which feedback is called

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
    TYPE(t_gridref_state), POINTER  :: p_grfp => NULL()
    TYPE(t_patch),      POINTER     :: p_pp => NULL()
    TYPE(t_patch),      POINTER     :: p_pc => NULL()
    TYPE(t_nwp_phy_diag), POINTER   :: prm_diagp => NULL()
    TYPE(t_nwp_phy_diag), POINTER   :: prm_diagc => NULL()

    ! Indices
    INTEGER :: jb, jc, jk, js, jt, je, jv, i_nchdom, i_chidx,  &
      i_startblk, i_endblk, i_startidx, i_endidx, ic, &
      i_rlend_c, i_rlend_e, i_nchdom_p

    INTEGER :: nlev_c            ! number of full levels (child dom)
    INTEGER :: nlev_p            ! number of full levels (parent dom)
    INTEGER :: nshift, nst_fbk
    INTEGER :: ntracer_fbk

    REAL(vp), ALLOCATABLE, DIMENSION(:,:,:), TARGET :: feedback_rho, feedback_thv,        &
      feedback_vn, feedback_w
    REAL(vp), ALLOCATABLE, TARGET :: feedback_rhoqx(:,:,:,:), feedback_aero(:,:,:)

    ! Note: as w(nlevp1) is diagnostic, it is excluded from feedback
    REAL(vp), DIMENSION(nproma,p_patch(jg)%nlev,p_patch(jgp)%nblks_c), TARGET :: &
      parent_rho, parent_thv, parent_w

    REAL(vp), DIMENSION(nproma,p_patch(jg)%nlev) :: diff_rho, diff_thv, diff_w

    REAL(vp),DIMENSION(nproma,p_patch(jg)%nlev,p_patch(jgp)%nblks_c,ntracer), TARGET :: &
      parent_rhoqx
    REAL(vp),DIMENSION(nproma,nclass_aero,p_patch(jgp)%nblks_c) :: parent_aero

    REAL(vp), DIMENSION(nproma,p_patch(jg)%nlev,p_patch(jgp)%nblks_e), TARGET :: &
      parent_vn, diff_vn

#ifdef __LOOP_EXCHANGE
    REAL(vp) :: rot_diff_vn(p_patch(jg)%nlev,nproma,p_patch(jgp)%nblks_v)
    REAL(vp) :: div_diff_vn(p_patch(jg)%nlev,nproma,p_patch(jgp)%nblks_c)
#else
    REAL(vp) :: rot_diff_vn(nproma,p_patch(jg)%nlev,p_patch(jgp)%nblks_v)
    REAL(vp) :: div_diff_vn(nproma,p_patch(jg)%nlev,p_patch(jgp)%nblks_c)
#endif
    REAL(vp) :: theta_v_pr(nproma,p_patch(jg)%nlev,p_patch(jg)%nblks_c)

    REAL(wp) :: z_fbk_rho(nproma,4,p_patch(jg)%nlev)

    REAL(wp) :: rd_o_cvd, rd_o_p0ref, relfac, dcoef_vec


    INTEGER,  DIMENSION(:,:,:), POINTER :: iccidx, iccblk, iceidx, iceblk, iveidx, iveblk, &
      iecidx, iecblk, ievidx, ievblk
    REAL(vp), DIMENSION(:,:,:), POINTER :: p_fbk_rho, p_fbk_thv, p_fbk_w, p_fbk_vn
    REAL(vp), DIMENSION(:,:,:,:), POINTER :: p_fbk_rhoqx
    REAL(wp), DIMENSION(:,:,:), POINTER :: p_fbkwgt, p_fbkwgt_e

    !-----------------------------------------------------------------------

    ! write(0,*) "n_dom_start,n_dom, jg, jgp=", n_dom_start, n_dom, jg, jgp
    IF (msg_level >= 10) THEN
      WRITE(message_text,'(a,i2,a,i2)') '========= Feedback:',jg,' =>',jgp
      CALL message(TRIM(routine),message_text)
    ENDIF

    p_parent_prog    => p_nh_state(jgp)%prog(nnew(jgp))
    p_parent_prog_rcf=> p_nh_state(jgp)%prog(nnew_rcf(jgp))
    p_child_prog     => p_nh_state(jg)%prog(nnow(jg))
    p_child_prog_rcf => p_nh_state(jg)%prog(nnow_rcf(jg))
    p_gcc            => p_patch(jg)%cells
    p_pc             => p_patch(jg)
    p_int            => p_int_state(jgp)
    prm_diagp        => prm_diag(jgp)
    prm_diagc        => prm_diag(jg)

    p_grf  => p_grf_state_local_parent(jg)
    p_grfp => p_grf_state(jgp)
    p_gcp  => p_patch_local_parent(jg)%cells
    p_gep  => p_patch_local_parent(jg)%edges
    p_pp   => p_patch_local_parent(jg)

    nlev_c   = p_pc%nlev

    nlev_p   = p_pp%nlev

    nshift = p_pc%nshift
    js     = nshift

    ! Exclude scalar variables from feedback in the upper three layers if the current domain is vertically nested
    ! This is needed in order to avoid temperature biases in the interface layer region
    IF (nshift == 0) THEN
      nst_fbk = 1
    ELSE
      nst_fbk = 4
    ENDIF

    IF (iforcing > 1) THEN  ! tracers represent moisture variables
      ntracer_fbk = 3       ! take only QV, QC and QI
    ELSE
      ntracer_fbk = ntracer
    ENDIF

    i_nchdom = MAX(1,p_pc%n_childdom)
    i_chidx  = p_pc%parent_child_index
    i_nchdom_p = MAX(1,p_patch(jgp)%n_childdom)

    ! end index levels for application of feedback relaxation
    i_rlend_c   = min_rlcell_int
    i_rlend_e   = min_rledge_int

    ! R/c_v (not present in physical constants)
    rd_o_cvd = 1._wp / cvd_o_rd

    ! R / p0ref
    rd_o_p0ref = rd / p0ref

    relfac = MIN(0.075_wp, dt_fbk/fbk_relax_timescale) ! default for relaxation time scale is 3 hours
    dcoef_vec  = 1._wp/12._wp

    ! Allocation of storage fields
    ! In parallel runs the lower bound of the feedback_* arrays must be 1 for use in exchange data,
    ! the lower bound of the fbk_* arrays must be i_startblk for the use in global sum.

    i_startblk = p_gcp%start_blk(grf_fbk_start_c,i_chidx)
    i_endblk   = p_gcp%end_blk(min_rlcell_int,i_chidx)

    i_startblk = 1

    ALLOCATE(feedback_thv(nproma, nlev_c, i_startblk:i_endblk),   &
      feedback_rho       (nproma, nlev_c, i_startblk:i_endblk),   &
      feedback_w         (nproma, nlev_c, i_startblk:i_endblk))

    IF(ltransport) &
      ALLOCATE(feedback_rhoqx(nproma, nlev_c, i_startblk:i_endblk, ntracer_fbk))

    IF(ltransport .AND. iprog_aero == 1) &
      ALLOCATE(feedback_aero(nproma, nclass_aero, i_startblk:i_endblk))

    i_startblk = p_gep%start_blk(grf_fbk_start_e,i_chidx)
    i_endblk   = p_gep%end_blk(min_rledge,i_chidx)

    i_startblk = 1

    ALLOCATE(feedback_vn(nproma, nlev_c, i_startblk:i_endblk))


    ! Set pointers to index and coefficient fields for cell-based variables
    iccidx => p_gcp%child_idx
    iccblk => p_gcp%child_blk

    iceidx => p_gep%child_idx
    iceblk => p_gep%child_blk

    ! Note: for consistency, there is no distiction between scalar and tracer feedback weights in this routine
    ! Temperature and vertical wind feedback is always bilinear
    IF (grf_scalfbk == 1) THEN
      p_fbkwgt    => p_grf%fbk_wgt_aw
    ELSE
      p_fbkwgt    => p_grf%fbk_wgt_bln
    ENDIF
    p_fbkwgt_e  => p_grf%fbk_wgt_e

    IF (p_test_run) THEN
      diff_vn  = 0._wp
    ENDIF

    ! 1. Feedback of child-domain variables to the parent grid

!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)

    ! Compute perturbation theta_v

    i_startblk = p_gcc%start_blk(grf_bdywidth_c+1,1)
    i_endblk   = p_gcc%end_blk(min_rlcell,i_nchdom)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk) ICON_OMP_DEFAULT_SCHEDULE
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

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk,jt,z_fbk_rho) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_pp, jb, i_startblk, i_endblk, &
        i_startidx, i_endidx, grf_fbk_start_c, min_rlcell_int)

#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
!DIR$ IVDEP
        DO jk = 1, nlev_c
#else
      DO jk = 1, nlev_c
        DO jc = i_startidx, i_endidx
#endif

          feedback_thv(jc,jk,jb) =                                                      &
            theta_v_pr(iccidx(jc,jb,1),jk,iccblk(jc,jb,1))*p_grf%fbk_wgt_bln(jc,jb,1) + &
            theta_v_pr(iccidx(jc,jb,2),jk,iccblk(jc,jb,2))*p_grf%fbk_wgt_bln(jc,jb,2) + &
            theta_v_pr(iccidx(jc,jb,3),jk,iccblk(jc,jb,3))*p_grf%fbk_wgt_bln(jc,jb,3) + &
            theta_v_pr(iccidx(jc,jb,4),jk,iccblk(jc,jb,4))*p_grf%fbk_wgt_bln(jc,jb,4)

          z_fbk_rho(jc,1,jk) =                                                   & 
            p_child_prog%rho(iccidx(jc,jb,1),jk,iccblk(jc,jb,1))*p_fbkwgt(jc,jb,1)
          z_fbk_rho(jc,2,jk) =                                                   &
            p_child_prog%rho(iccidx(jc,jb,2),jk,iccblk(jc,jb,2))*p_fbkwgt(jc,jb,2)
          z_fbk_rho(jc,3,jk) =                                                   &
            p_child_prog%rho(iccidx(jc,jb,3),jk,iccblk(jc,jb,3))*p_fbkwgt(jc,jb,3)
          z_fbk_rho(jc,4,jk) =                                                   &
            p_child_prog%rho(iccidx(jc,jb,4),jk,iccblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

          ! The semi-empirical correction for density feedback is necessitated by the fact that
          ! increased orography resolution goes along with increased mass. This needs to be
          ! corrected for in order to avoid a systematic mass drift in two-way nesting
          feedback_rho(jc,jk,jb) =                                                               &
           (z_fbk_rho(jc,1,jk)+z_fbk_rho(jc,2,jk)+z_fbk_rho(jc,3,jk)+z_fbk_rho(jc,4,jk)) -       &
           (1.05_wp-0.005_wp*feedback_thv(jc,jk,jb))*p_nh_state(jg)%metrics%rho_ref_corr(jc,jk,jb) 

          feedback_w(jc,jk,jb) =                                                            &
            p_child_prog%w(iccidx(jc,jb,1),jk,iccblk(jc,jb,1))*p_grf%fbk_wgt_bln(jc,jb,1) + &
            p_child_prog%w(iccidx(jc,jb,2),jk,iccblk(jc,jb,2))*p_grf%fbk_wgt_bln(jc,jb,2) + &
            p_child_prog%w(iccidx(jc,jb,3),jk,iccblk(jc,jb,3))*p_grf%fbk_wgt_bln(jc,jb,3) + &
            p_child_prog%w(iccidx(jc,jb,4),jk,iccblk(jc,jb,4))*p_grf%fbk_wgt_bln(jc,jb,4)

        ENDDO
      ENDDO

      IF (ltransport) THEN ! tracer mass feedback
#ifdef __LOOP_EXCHANGE
        DO jc = i_startidx, i_endidx
          DO jt = 1, ntracer_fbk
!DIR$ IVDEP
            DO jk = 1, nlev_c
#else
        DO jt = 1, ntracer_fbk
          DO jk = 1, nlev_c
            DO jc = i_startidx, i_endidx
#endif
              feedback_rhoqx(jc,jk,jb,jt) =                                       &
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

      IF (ltransport .AND. iprog_aero == 1) THEN
#ifdef __LOOP_EXCHANGE
        DO jc = i_startidx, i_endidx
!DIR$ IVDEP
          DO jt = 1, nclass_aero
#else
        DO jt = 1, nclass_aero
          DO jc = i_startidx, i_endidx
#endif
            feedback_aero(jc,jt,jb) =                                                            &
              prm_diagc%aerosol(iccidx(jc,jb,1),jt,iccblk(jc,jb,1))*p_grf%fbk_wgt_bln(jc,jb,1) + &
              prm_diagc%aerosol(iccidx(jc,jb,2),jt,iccblk(jc,jb,2))*p_grf%fbk_wgt_bln(jc,jb,2) + &
              prm_diagc%aerosol(iccidx(jc,jb,3),jt,iccblk(jc,jb,3))*p_grf%fbk_wgt_bln(jc,jb,3) + &
              prm_diagc%aerosol(iccidx(jc,jb,4),jt,iccblk(jc,jb,4))*p_grf%fbk_wgt_bln(jc,jb,4)
          ENDDO
        ENDDO
      ENDIF

    ENDDO
!$OMP END DO

    ! Velocity feedback

    i_startblk = p_gep%start_blk(grf_fbk_start_e,i_chidx)
    i_endblk   = p_gep%end_blk(min_rledge_int,i_chidx)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,je,jk) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(p_pp, jb, i_startblk, i_endblk, &
        i_startidx, i_endidx, grf_fbk_start_e, min_rledge_int)


#ifdef __LOOP_EXCHANGE
      DO je = i_startidx, i_endidx
!DIR$ IVDEP
        DO jk = 1, nlev_c
#else
      DO jk = 1, nlev_c
        DO je = i_startidx, i_endidx
#endif

          feedback_vn(je,jk,jb) =                                                     &
            p_child_prog%vn(iceidx(je,jb,1),jk,iceblk(je,jb,1))*p_fbkwgt_e(je,jb,1) + &
            p_child_prog%vn(iceidx(je,jb,2),jk,iceblk(je,jb,2))*p_fbkwgt_e(je,jb,2)
        ENDDO
      ENDDO

    ENDDO
!$OMP END DO NOWAIT

!$OMP END PARALLEL

#ifdef __MIXED_PRECISION
    CALL exchange_data_mult_mixprec(p_pp%comm_pat_loc_to_glb_c_fbk, 0, 0, 3, 3*nlev_c, &
      RECV1_SP=parent_rho,     SEND1_SP=feedback_rho,      &
      RECV2_SP=parent_thv,     SEND2_SP=feedback_thv,      &
      RECV3_SP=parent_w,       SEND3_SP=feedback_w         )


    CALL exchange_data_mult_mixprec(p_pp%comm_pat_loc_to_glb_e_fbk, 0, 0, 1, nlev_c, &
      RECV1_SP=parent_vn, SEND1_SP=feedback_vn )

    IF (ltransport .AND. iprog_aero == 1) THEN
      CALL exchange_data_mult_mixprec(p_pp%comm_pat_loc_to_glb_c_fbk, 0, 0, ntracer_fbk+1, ntracer_fbk*nlev_c+nclass_aero, &
        RECV1_SP=parent_aero,     SEND1_SP=feedback_aero,                    &
        RECV4D_SP=parent_rhoqx(:,:,:,1:ntracer_fbk), SEND4D_SP=feedback_rhoqx)
    ELSE IF (ltransport) THEN
      CALL exchange_data_mult_mixprec(p_pp%comm_pat_loc_to_glb_c_fbk, 0, 0, ntracer_fbk, ntracer_fbk*nlev_c, &
        RECV4D_SP=parent_rhoqx(:,:,:,1:ntracer_fbk), SEND4D_SP=feedback_rhoqx)
    ENDIF
#else
    CALL exchange_data_mult(p_pp%comm_pat_loc_to_glb_c_fbk, 3, 3*nlev_c, &
      RECV1=parent_rho,     SEND1=feedback_rho,      &
      RECV2=parent_thv,     SEND2=feedback_thv,      &
      RECV3=parent_w,       SEND3=feedback_w         )


    CALL exchange_data_mult(p_pp%comm_pat_loc_to_glb_e_fbk, 1, nlev_c, &
      RECV1=parent_vn, SEND1=feedback_vn )

    IF (ltransport .AND. iprog_aero == 1) THEN
      CALL exchange_data_mult(p_pp%comm_pat_loc_to_glb_c_fbk, ntracer_fbk+1, ntracer_fbk*nlev_c+nclass_aero, &
        RECV1=parent_aero,     SEND1=feedback_aero,                    &
        RECV4D=parent_rhoqx(:,:,:,1:ntracer_fbk), SEND4D=feedback_rhoqx)
    ELSE IF (ltransport) THEN
      CALL exchange_data_mult(p_pp%comm_pat_loc_to_glb_c_fbk, ntracer_fbk, ntracer_fbk*nlev_c, &
        RECV4D=parent_rhoqx(:,:,:,1:ntracer_fbk), SEND4D=feedback_rhoqx)
    ENDIF
#endif

    p_fbk_rho => parent_rho
    p_fbk_thv => parent_thv
    p_fbk_w   => parent_w
    p_fbk_vn  => parent_vn

    IF (ltransport) p_fbk_rhoqx => parent_rhoqx


    ! 2. Compute differences between feedback velocity and corresponding parent field,
    !    smooth velocity increment and execute velocity relaxation
!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)

    i_startblk = p_patch(jgp)%edges%start_blk(1,1)
    i_endblk   = p_patch(jgp)%edges%end_blk(i_rlend_e,i_nchdom_p)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,je,jk) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      diff_vn (:,:,jb) = 0._wp

      CALL get_indices_e(p_patch(jgp), jb, i_startblk, i_endblk, i_startidx, i_endidx, 1, i_rlend_e)

#ifdef __LOOP_EXCHANGE
      DO je = i_startidx,i_endidx
        IF (p_grfp%mask_ovlp_e(je,jb,i_chidx)) THEN
!DIR$ IVDEP
          DO jk = 1, nlev_c
#else
      DO jk = 1, nlev_c
        DO je = i_startidx,i_endidx
          IF (p_grfp%mask_ovlp_e(je,jb,i_chidx)) THEN
#endif
            diff_vn(je,jk,jb) = p_fbk_vn(je,jk,jb) - p_parent_prog%vn(je,jk+js,jb)

#ifdef __LOOP_EXCHANGE
          ENDDO
        ENDIF
#else
          ENDIF
        ENDDO
#endif
      ENDDO

    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    CALL sync_patch_array(SYNC_E,p_patch(jgp),diff_vn)

    ! 2a. Smoothing of velocity feedback-parent differences 

    iceidx => p_patch(jgp)%cells%edge_idx
    iceblk => p_patch(jgp)%cells%edge_blk

    iveidx => p_patch(jgp)%verts%edge_idx
    iveblk => p_patch(jgp)%verts%edge_blk

    iecidx => p_patch(jgp)%edges%cell_idx
    iecblk => p_patch(jgp)%edges%cell_blk

    ievidx => p_patch(jgp)%edges%vertex_idx
    ievblk => p_patch(jgp)%edges%vertex_blk


!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)

    i_startblk = p_patch(jgp)%verts%start_blk(1,1)
    i_endblk   = p_patch(jgp)%verts%end_blk(min_rlvert_int-1,i_nchdom_p)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jv,jk) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_v(p_patch(jgp), jb, i_startblk, i_endblk, i_startidx, i_endidx, 1, min_rlvert_int-1)

#ifdef __LOOP_EXCHANGE
      DO jv = i_startidx, i_endidx
        IF (p_grfp%mask_ovlp_v(jv,jb,i_chidx)) THEN
          DO jk = 1, nlev_c
            rot_diff_vn(jk,jv,jb) =   &
#else
      DO jk = 1, nlev_c
        DO jv = i_startidx, i_endidx
          IF (p_grfp%mask_ovlp_v(jv,jb,i_chidx)) THEN
            rot_diff_vn(jv,jk,jb) =   &
#endif
              diff_vn(iveidx(jv,jb,1),jk,iveblk(jv,jb,1))*p_int%geofac_rot(jv,1,jb) + &
              diff_vn(iveidx(jv,jb,2),jk,iveblk(jv,jb,2))*p_int%geofac_rot(jv,2,jb) + &
              diff_vn(iveidx(jv,jb,3),jk,iveblk(jv,jb,3))*p_int%geofac_rot(jv,3,jb) + &
              diff_vn(iveidx(jv,jb,4),jk,iveblk(jv,jb,4))*p_int%geofac_rot(jv,4,jb) + &
              diff_vn(iveidx(jv,jb,5),jk,iveblk(jv,jb,5))*p_int%geofac_rot(jv,5,jb) + &
              diff_vn(iveidx(jv,jb,6),jk,iveblk(jv,jb,6))*p_int%geofac_rot(jv,6,jb)

#ifdef __LOOP_EXCHANGE
          ENDDO
        ENDIF
#else
          ENDIF
        ENDDO
#endif
      ENDDO
    ENDDO
!$OMP END DO

    i_startblk = p_patch(jgp)%cells%start_blk(1,1)
    i_endblk   = p_patch(jgp)%cells%end_blk(min_rlcell_int-1,i_nchdom_p)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch(jgp), jb, i_startblk, i_endblk, i_startidx, i_endidx, 1, min_rlcell_int-1)

#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        IF (p_grfp%mask_ovlp_ch(jc,jb,i_chidx)) THEN
          DO jk = 1, nlev_c
            div_diff_vn(jk,jc,jb) =   &
#else
      DO jk = 1, nlev_c
        DO jc = i_startidx, i_endidx
          IF (p_grfp%mask_ovlp_ch(jc,jb,i_chidx)) THEN
            div_diff_vn(jc,jk,jb) =   &
#endif
              diff_vn(iceidx(jc,jb,1),jk,iceblk(jc,jb,1))*p_int%geofac_div(jc,1,jb) + &
              diff_vn(iceidx(jc,jb,2),jk,iceblk(jc,jb,2))*p_int%geofac_div(jc,2,jb) + &
              diff_vn(iceidx(jc,jb,3),jk,iceblk(jc,jb,3))*p_int%geofac_div(jc,3,jb)

#ifdef __LOOP_EXCHANGE
          ENDDO
        ENDIF
#else
          ENDIF
        ENDDO
#endif
      ENDDO
    ENDDO
!$OMP END DO

    i_startblk = p_patch(jgp)%edges%start_blk(1,1)
    i_endblk   = p_patch(jgp)%edges%end_blk(i_rlend_e,i_nchdom_p)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,je,jk) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(p_patch(jgp), jb, i_startblk, i_endblk, i_startidx, i_endidx, 1, i_rlend_e)

#ifdef __LOOP_EXCHANGE
      DO je = i_startidx, i_endidx
        IF (p_grfp%mask_ovlp_e(je,jb,i_chidx)) THEN
          DO jk = 1, nlev_c
            diff_vn(je,jk,jb) = diff_vn(je,jk,jb) +                &
              dcoef_vec * p_patch(jgp)%edges%area_edge(je,jb) *    &
              ( p_patch(jgp)%edges%tangent_orientation(je,jb) *     &
              ( rot_diff_vn(jk,ievidx(je,jb,2),ievblk(je,jb,2))    &
              - rot_diff_vn(jk,ievidx(je,jb,1),ievblk(je,jb,1)) )  &
              * p_patch(jgp)%edges%inv_primal_edge_length(je,jb) + &
              ( div_diff_vn(jk,iecidx(je,jb,2),iecblk(je,jb,2))    &
              - div_diff_vn(jk,iecidx(je,jb,1),iecblk(je,jb,1)) )  &
              * p_patch(jgp)%edges%inv_dual_edge_length(je,jb)     )
          ENDDO
        ENDIF
#else
      DO jk = 1, nlev_c
        DO je = i_startidx, i_endidx
          IF (p_grfp%mask_ovlp_e(je,jb,i_chidx)) THEN
            diff_vn(je,jk,jb) = diff_vn(je,jk,jb) +                &
              dcoef_vec * p_patch(jgp)%edges%area_edge(je,jb) *    &
              ( p_patch(jgp)%edges%tangent_orientation(je,jb) *     &
              ( rot_diff_vn(ievidx(je,jb,2),jk,ievblk(je,jb,2))    &
              - rot_diff_vn(ievidx(je,jb,1),jk,ievblk(je,jb,1)) )  &
              * p_patch(jgp)%edges%inv_primal_edge_length(je,jb) + &
              ( div_diff_vn(iecidx(je,jb,2),jk,iecblk(je,jb,2))    &
              - div_diff_vn(iecidx(je,jb,1),jk,iecblk(je,jb,1)) )  &
              * p_patch(jgp)%edges%inv_dual_edge_length(je,jb)     )
          ENDIF
        ENDDO
#endif
      ENDDO

      ! 2b. Execute relaxation
#ifdef __LOOP_EXCHANGE
      DO je = i_startidx,i_endidx
        IF (p_grfp%mask_ovlp_e(je,jb,i_chidx)) THEN
!DIR$ IVDEP
          DO jk = nshift+1,nlev_p
#else
      DO jk = nshift+1,nlev_p
        DO je = i_startidx,i_endidx
          IF (p_grfp%mask_ovlp_e(je,jb,i_chidx)) THEN
#endif
            p_parent_prog%vn(je,jk,jb) = p_parent_prog%vn(je,jk,jb) + &
              relfac*diff_vn(je,jk-js,jb)

#ifdef __LOOP_EXCHANGE
          ENDDO
        ENDIF
#else
          ENDIF
        ENDDO
#endif
      ENDDO
    ENDDO
!$OMP END DO


    ! 3. The same for mass-point variables

    i_startblk = p_patch(jgp)%cells%start_blk(1,1)
    i_endblk   = p_patch(jgp)%cells%end_blk(i_rlend_c,i_nchdom_p)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk,jt,diff_rho,diff_thv,diff_w) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch(jgp), jb, i_startblk, i_endblk, i_startidx, i_endidx, 1, i_rlend_c)

      ! Compute differences between feedback fields and corresponding parent fields
#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx,i_endidx
        IF (p_grfp%mask_ovlp_c(jc,jb,i_chidx)) THEN
!DIR$ IVDEP
          DO jk = nst_fbk, nlev_c
#else
      DO jk = nst_fbk, nlev_c
        DO jc = i_startidx,i_endidx
          IF (p_grfp%mask_ovlp_c(jc,jb,i_chidx)) THEN
#endif

            ! density
            diff_rho(jc,jk) = p_fbk_rho(jc,jk,jb) - p_parent_prog%rho(jc,jk+js,jb)

            ! theta_v
            diff_thv(jc,jk) = p_fbk_thv(jc,jk,jb) - p_parent_prog%theta_v(jc,jk+js,jb) + &
              p_nh_state(jgp)%metrics%theta_ref_mc(jc,jk+js,jb)

            ! w
            diff_w(jc,jk) = p_fbk_w(jc,jk,jb) - p_parent_prog%w(jc,jk+js,jb)

#ifdef __LOOP_EXCHANGE
          ENDDO
        ENDIF
#else
          ENDIF
        ENDDO
#endif
      ENDDO

      ! Relaxation of tracer variables
      IF (ltransport) THEN
#ifdef __LOOP_EXCHANGE
        DO jc = i_startidx,i_endidx
          IF (p_grfp%mask_ovlp_c(jc,jb,i_chidx)) THEN
            DO jt = 1, ntracer_fbk
!DIR$ IVDEP
              DO jk = nshift+nst_fbk, nlev_p
#else
        DO jt = 1, ntracer_fbk
          DO jk = nshift+nst_fbk, nlev_p
            DO jc = i_startidx,i_endidx
              IF (p_grfp%mask_ovlp_c(jc,jb,i_chidx)) THEN
#endif
                p_parent_prog_rcf%tracer(jc,jk,jb,jt) = p_parent_prog_rcf%tracer(jc,jk,jb,jt) + &
                  relfac * ( p_fbk_rhoqx(jc,jk-js,jb,jt) / ( p_parent_prog%rho(jc,jk,jb) +      &
                  diff_rho(jc,jk-js) ) -  p_parent_prog_rcf%tracer(jc,jk,jb,jt) )

#ifdef __LOOP_EXCHANGE
              ENDDO
            ENDDO
          ENDIF
#else
              ENDIF
            ENDDO
          ENDDO
#endif
        ENDDO
      ENDIF

      IF (ltransport .AND. iprog_aero == 1) THEN

        DO jt = 1, nclass_aero
          DO jc = i_startidx,i_endidx
            IF (p_grfp%mask_ovlp_c(jc,jb,i_chidx)) THEN
              prm_diagp%aerosol(jc,jt,jb) = prm_diagp%aerosol(jc,jt,jb) + relfac *  &
                (parent_aero(jc,jt,jb) - prm_diagp%aerosol(jc,jt,jb))
            ENDIF
          ENDDO
        ENDDO

      ENDIF

      ! Relaxation of dynamical variables
#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx,i_endidx
        IF (p_grfp%mask_ovlp_c(jc,jb,i_chidx)) THEN
!DIR$ IVDEP
          DO jk = nshift+nst_fbk,nlev_p
#else
      DO jk = nshift+nst_fbk,nlev_p
        DO jc = i_startidx,i_endidx
          IF (p_grfp%mask_ovlp_c(jc,jb,i_chidx)) THEN
#endif

            ! density
            p_parent_prog%rho(jc,jk,jb) = p_parent_prog%rho(jc,jk,jb) + relfac*diff_rho(jc,jk-js)

            ! theta_v
            p_parent_prog%theta_v(jc,jk,jb) = p_parent_prog%theta_v(jc,jk,jb) + relfac*diff_thv(jc,jk-js)

            ! w
            p_parent_prog%w(jc,jk,jb) = p_parent_prog%w(jc,jk,jb) + relfac*diff_w(jc,jk-js)

            ! exner is diagnosed from rho*theta_v
            p_parent_prog%exner(jc,jk,jb) = EXP(rd_o_cvd*LOG(rd_o_p0ref* &
              p_parent_prog%theta_v(jc,jk,jb)*p_parent_prog%rho(jc,jk,jb)))

#ifdef __LOOP_EXCHANGE
          ENDDO
        ENDIF
#else
          ENDIF
        ENDDO
#endif
      ENDDO

    ENDDO
!$OMP END DO
!$OMP END PARALLEL

    CALL sync_patch_array(SYNC_E,p_patch(jgp),p_parent_prog%vn)

    IF (ltransport .AND. iprog_aero == 1) THEN
      CALL sync_patch_array_mult(SYNC_C, p_patch(jgp), ntracer_fbk+4, p_parent_prog%rho, p_parent_prog%theta_v, &
        p_parent_prog%w, prm_diagp%aerosol, f4din=p_parent_prog_rcf%tracer(:,:,:,1:ntracer_fbk))
    ELSE IF (ltransport) THEN
      CALL sync_patch_array_mult(SYNC_C, p_patch(jgp), ntracer_fbk+3, p_parent_prog%rho, p_parent_prog%theta_v, &
        p_parent_prog%w, f4din=p_parent_prog_rcf%tracer(:,:,:,1:ntracer_fbk))
    ELSE
      CALL sync_patch_array_mult(SYNC_C,p_patch(jgp),3,p_parent_prog%rho,p_parent_prog%theta_v,   &
        p_parent_prog%w)
    ENDIF

    ! Recomputation of exner on halo points (saves communication of one field)
#ifdef __LOOP_EXCHANGE
    DO ic = 1, p_nh_state(jgp)%metrics%ovlp_halo_c_dim(i_chidx)
      jc = p_nh_state(jgp)%metrics%ovlp_halo_c_idx(ic,i_chidx)
      jb = p_nh_state(jgp)%metrics%ovlp_halo_c_blk(ic,i_chidx)
!DIR$ IVDEP
      DO jk = nshift+nst_fbk, nlev_p
#else
    DO jk = nshift+nst_fbk, nlev_p
!CDIR NODEP,VOVERTAKE,VOB
      DO ic = 1, p_nh_state(jgp)%metrics%ovlp_halo_c_dim(i_chidx)
        jc = p_nh_state(jgp)%metrics%ovlp_halo_c_idx(ic,i_chidx)
        jb = p_nh_state(jgp)%metrics%ovlp_halo_c_blk(ic,i_chidx)
#endif

        p_parent_prog%exner(jc,jk,jb) =   &
          EXP(rd_o_cvd*LOG(rd_o_p0ref*p_parent_prog%theta_v(jc,jk,jb)*p_parent_prog%rho(jc,jk,jb)))

      ENDDO
    ENDDO


    DEALLOCATE(feedback_thv,feedback_rho,feedback_w,feedback_vn)
    IF (ltransport) DEALLOCATE(feedback_rhoqx)
    IF (ltransport .AND. iprog_aero == 1) DEALLOCATE(feedback_aero)

  END SUBROUTINE relax_feedback

END MODULE mo_nh_feedback
