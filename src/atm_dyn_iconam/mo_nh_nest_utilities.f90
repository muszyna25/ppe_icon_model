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

MODULE mo_nh_nest_utilities
  !
  !
  USE mo_kind,                ONLY: wp, vp
  USE mo_exception,           ONLY: message_text, message, finish
  USE mo_model_domain,        ONLY: t_patch, t_grid_cells, t_grid_edges, p_patch_local_parent, &
    p_patch
  USE mo_grid_config,         ONLY: n_dom
  USE mo_intp_data_strc,      ONLY: t_int_state, p_int_state, p_int_state_local_parent
  USE mo_grf_intp_data_strc,  ONLY: t_gridref_state, p_grf_state, p_grf_state_local_parent
  USE mo_gridref_config,      ONLY: grf_intmethod_c, grf_intmethod_e, grf_intmethod_ct, grf_scalfbk, grf_tracfbk
  USE mo_grf_bdyintp,         ONLY: interpol_scal_grf, interpol_vec_grf, interpol2_vec_grf
  USE mo_grf_nudgintp,        ONLY: interpol_scal_nudging, interpol_vec_nudging
  USE mo_grf_ubcintp,         ONLY: interpol_scal_ubc,interpol_vec_ubc
  USE mo_dynamics_config,     ONLY: nnow, nsav1, nnow_rcf
  USE mo_parallel_config,     ONLY: nproma, p_test_run
  USE mo_run_config,          ONLY: ltransport, msg_level, ntracer, lvert_nest, iqv, iqc, iforcing
  USE mo_nonhydro_types,      ONLY: t_nh_state, t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_nonhydro_state,      ONLY: p_nh_state
  USE mo_nwp_phy_state,       ONLY: prm_diag
  USE mo_nonhydrostatic_config,ONLY: ndyn_substeps_var
  USE mo_atm_phy_nwp_config,  ONLY: iprog_aero
  USE mo_impl_constants,      ONLY: min_rlcell_int, min_rledge_int, MAX_CHAR_LENGTH, min_rlcell, min_rledge
  USE mo_loopindices,         ONLY: get_indices_c, get_indices_e
  USE mo_impl_constants_grf,  ONLY: grf_bdyintp_start_c,                       &
    grf_bdyintp_end_c,                         &
    grf_fbk_start_c,                           &
    grf_bdywidth_c, grf_bdywidth_e,            &
    grf_nudgintp_start_c, grf_nudgintp_start_e,&
    grf_nudge_start_c
  USE mo_mpi,                 ONLY: my_process_is_mpi_parallel
  USE mo_communication,       ONLY: exchange_data, exchange_data_mult
  USE mo_sync,                ONLY: SYNC_C, SYNC_E, sync_patch_array, &
    global_sum_array3, sync_patch_array_mult
  USE mo_physical_constants,  ONLY: rd, cvd_o_rd, p0ref, vtmpc1, cpd
  USE mo_limarea_config,      ONLY: latbc_config
  USE mo_initicon_types,      ONLY: t_pi_atm
  USE mo_advection_config,    ONLY: advection_config
  USE mo_nudging_config,      ONLY: nudging_config, indg_type, ithermdyn_type

  IMPLICIT NONE

  PRIVATE


  PUBLIC :: compute_tendencies, boundary_interpolation, complete_nesting_setup, &
    prep_bdy_nudging, nest_boundary_nudging, save_progvars,                     &
    prep_rho_bdy_nudging, density_boundary_nudging, limarea_bdy_nudging

CONTAINS

  !>
  !! Computes geometric information needed for the conservation correction
  !! in the feedback routine
  !!
  !! @par Revision History
  !! Developed by Guenther Zaengl, DWD, 2010-05-05
  !! (Note: since we have not modified 'src/atm_dyn_iconam/mo_nh_feedback: feedback' 
  !! for the deep atmosphere, there is no reason to modify 'fbk_dom_volume'.)
  !!
  SUBROUTINE complete_nesting_setup


    TYPE(t_patch),      POINTER     :: p_pp => NULL()
    TYPE(t_patch),      POINTER     :: p_pc => NULL()
    TYPE(t_patch),      POINTER     :: p_lp => NULL()
    TYPE(t_grid_cells), POINTER     :: p_gcp => NULL()

    REAL(wp), ALLOCATABLE :: cell_volume(:,:,:)
    REAL(vp), ALLOCATABLE :: z_rho_ref(:,:,:)

    INTEGER :: jg, ji, jgc, jb, jc, jk, jks
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom
    INTEGER :: nlev, nlev_c, nshift

    INTEGER, DIMENSION(:,:,:), POINTER :: iidx, iblk
    REAL(wp), DIMENSION(:,:,:), POINTER :: p_fbkwgt

    ALLOCATE(p_nh_state(1)%metrics%fbk_dom_volume(p_patch(1)%nlev))
    p_nh_state(1)%metrics%fbk_dom_volume(:) = 0._wp

    DO jg = 1, n_dom-1

      p_pp   => p_patch(jg)
      p_gcp  => p_pp%cells

      nlev = p_pp%nlev

      i_nchdom = p_pp%n_childdom
      IF (i_nchdom == 0) CYCLE

      DO ji = 1, i_nchdom

        jgc    =  p_pp%child_id(ji)
        p_pc   => p_patch(jgc)

        ! Note: the number of levels of fbk_dom_volume is that of the parent grid
        ALLOCATE(p_nh_state(jgc)%metrics%fbk_dom_volume(nlev))

        ! Copy layer thicknesses to local parent of current child

        p_lp => p_patch_local_parent(jgc)
        ALLOCATE(p_lp%cells%ddqz_z_full(nproma, p_lp%nlev, p_lp%n_patch_cells))
        p_lp%cells%ddqz_z_full(:,:,:) = 0._wp ! Safety only
        CALL exchange_data(p_lp%comm_pat_glb_to_loc_c, p_lp%cells%ddqz_z_full, &
          &                p_nh_state(jg)%metrics%ddqz_z_full)

        i_startblk = p_gcp%start_blk(1,1)
        i_endblk   = p_gcp%end_blk(min_rlcell_int,i_nchdom)

        ALLOCATE(cell_volume(nproma, nlev, i_endblk))

        cell_volume = 0._wp

        DO jb = i_startblk, i_endblk

          CALL get_indices_c(p_pp, jb, i_startblk, i_endblk, i_startidx, i_endidx, &
            1, min_rlcell_int)

          DO jk = 1, nlev
            DO jc = i_startidx, i_endidx
              ! Sum must be taken over owned points of nest overlap region only
              IF (p_gcp%refin_ctrl(jc,jb) <= grf_fbk_start_c .AND. p_gcp%child_id(jc,jb) == jgc .AND. &
                p_gcp%decomp_info%owner_mask(jc,jb)) &
                cell_volume(jc,jk,jb) = p_gcp%area(jc,jb)*p_nh_state(jg)%metrics%ddqz_z_full(jc,jk,jb)
            ENDDO
          ENDDO

        ENDDO

        p_nh_state(jgc)%metrics%fbk_dom_volume(:) = global_sum_array3(1,.FALSE.,cell_volume)

        DEALLOCATE(cell_volume)

      ENDDO

      ! Second part: compute correction term needed to use perturbation density for boundary nudging

      DO ji = 1, i_nchdom

        jgc    =  p_pp%child_id(ji)
        p_pc   => p_patch(jgc)

        nlev_c = p_pc%nlev
        nshift = p_pc%nshift

        p_fbkwgt => p_grf_state_local_parent(jgc)%fbk_wgt_bln
        p_gcp => p_patch_local_parent(jgc)%cells
        p_pp  => p_patch_local_parent(jgc)

        iidx  => p_gcp%child_idx
        iblk  => p_gcp%child_blk

        i_startblk = p_gcp%start_blk(grf_nudgintp_start_c+1,ji)
        i_endblk   = p_gcp%end_blk(min_rlcell_int,ji)

        ALLOCATE(p_nh_state(jgc)%metrics%rho_ref_corr(nproma, nlev_c, p_pp%nblks_c), &
          &      z_rho_ref(nproma, nlev, p_pp%nblks_c))
        z_rho_ref(:,:,:) = 0._vp
        p_nh_state(jgc)%metrics%rho_ref_corr(:,:,:) = 0._vp

        CALL exchange_data(p_pp%comm_pat_glb_to_loc_c, RECV=z_rho_ref, &
          &                SEND=p_nh_state(jg)%metrics%rho_ref_mc)

        DO jb = i_startblk, i_endblk

          CALL get_indices_c(p_pp, jb, i_startblk, i_endblk, i_startidx, i_endidx, &
            grf_nudgintp_start_c+1, min_rlcell_int)

          DO jk = 1, nlev_c
            jks = jk + nshift
            DO jc = i_startidx, i_endidx

              p_nh_state(jgc)%metrics%rho_ref_corr(jc,jk,jb) = - z_rho_ref(jc,jks,jb) + (           &
                p_nh_state(jgc)%metrics%rho_ref_mc(iidx(jc,jb,1),jk,iblk(jc,jb,1))*p_fbkwgt(jc,jb,1)+ &
                p_nh_state(jgc)%metrics%rho_ref_mc(iidx(jc,jb,2),jk,iblk(jc,jb,2))*p_fbkwgt(jc,jb,2)+ &
                p_nh_state(jgc)%metrics%rho_ref_mc(iidx(jc,jb,3),jk,iblk(jc,jb,3))*p_fbkwgt(jc,jb,3)+ &
                p_nh_state(jgc)%metrics%rho_ref_mc(iidx(jc,jb,4),jk,iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)  )

            ENDDO
          ENDDO
        ENDDO

        DEALLOCATE(z_rho_ref)

      ENDDO
    ENDDO

  END SUBROUTINE complete_nesting_setup


  !-------------------------------------------------------------------------
  !
  !>
  !! Saves the dynamic prognostic variables needed afterwards for computing
  !! the lateral boundary tendencies for nested domains
  !!
  !! @par Revision History
  !! Developed by Guenther Zaengl, DWD, 2015-10-08
  !!
  SUBROUTINE save_progvars (jg,p_nh_prog,p_nh_save)

    INTEGER,         INTENT(IN)    :: jg
    TYPE(t_nh_prog), INTENT(IN)    :: p_nh_prog
    TYPE(t_nh_prog), INTENT(INOUT) :: p_nh_save

    ! local variables
    !
    INTEGER :: ib, jb, ic, jc, ie, je, jk, jshift
    INTEGER :: nlev, nlevp1
    INTEGER :: nproma_bdyintp, nblks_bdyintp, npromz_bdyintp, nlen

    TYPE(t_gridref_state), POINTER :: p_grf

    !-----------------------------------------------------------------------


    p_grf  => p_grf_state(jg)

    ! number of vertical levels
    nlev   = p_patch(jg)%nlev
    nlevp1 = p_patch(jg)%nlevp1

    ! for dynamic nproma blocking
    nproma_bdyintp = MIN(nproma,256)

!$OMP PARALLEL PRIVATE(nblks_bdyintp,npromz_bdyintp)

    ! cell-based variables

    ! parameters for dynamic nproma blocking
    nblks_bdyintp  = INT(p_grf%npoints_bdyintp_src_c/nproma_bdyintp)
    npromz_bdyintp = MOD(p_grf%npoints_bdyintp_src_c,nproma_bdyintp)
    IF (npromz_bdyintp > 0) THEN
      nblks_bdyintp = nblks_bdyintp + 1
    ELSE
      npromz_bdyintp = nproma_bdyintp
    ENDIF

!$OMP DO PRIVATE(ib,jb,nlen,ic,jc,jk,jshift) ICON_OMP_DEFAULT_SCHEDULE
    DO ib = 1, nblks_bdyintp
      IF (ib == nblks_bdyintp) THEN
        nlen = npromz_bdyintp
      ELSE
        nlen = nproma_bdyintp
      ENDIF
      jshift = (ib-1)*nproma_bdyintp

#ifdef __LOOP_EXCHANGE
      DO ic = jshift+1, jshift+nlen
        jc = p_grf%idxlist_bdyintp_src_c(ic)
        jb = p_grf%blklist_bdyintp_src_c(ic)
!DIR$ IVDEP
        DO jk = 1, nlev
#else
      DO jk = 1, nlev
!CDIR NODEP
        DO ic = jshift+1, jshift+nlen
          jc = p_grf%idxlist_bdyintp_src_c(ic)
          jb = p_grf%blklist_bdyintp_src_c(ic)
#endif

          p_nh_save%w(jc,jk,jb)       = p_nh_prog%w(jc,jk,jb)
          p_nh_save%rho(jc,jk,jb)     = p_nh_prog%rho(jc,jk,jb)
          p_nh_save%theta_v(jc,jk,jb) = p_nh_prog%theta_v(jc,jk,jb)

        ENDDO
      ENDDO

      DO ic = jshift+1, jshift+nlen
        jc = p_grf%idxlist_bdyintp_src_c(ic)
        jb = p_grf%blklist_bdyintp_src_c(ic)
        p_nh_save%w(jc,nlevp1,jb)  = p_nh_prog%w(jc,nlevp1,jb)
      ENDDO

    ENDDO
!$OMP END DO

    ! edge-based variables

    ! parameters for dynamic nproma blocking
    nblks_bdyintp  = INT(p_grf%npoints_bdyintp_src_e/nproma_bdyintp)
    npromz_bdyintp = MOD(p_grf%npoints_bdyintp_src_e,nproma_bdyintp)
    IF (npromz_bdyintp > 0) THEN
      nblks_bdyintp = nblks_bdyintp + 1
    ELSE
      npromz_bdyintp = nproma_bdyintp
    ENDIF

!$OMP DO PRIVATE(ib,jb,nlen,ie,je,jk,jshift) ICON_OMP_DEFAULT_SCHEDULE
    DO ib = 1, nblks_bdyintp
      IF (ib == nblks_bdyintp) THEN
        nlen = npromz_bdyintp
      ELSE
        nlen = nproma_bdyintp
      ENDIF
      jshift = (ib-1)*nproma_bdyintp

#ifdef __LOOP_EXCHANGE
      DO ie = jshift+1, jshift+nlen
        je = p_grf%idxlist_bdyintp_src_e(ie)
        jb = p_grf%blklist_bdyintp_src_e(ie)
!DIR$ IVDEP
        DO jk = 1, nlev
#else
      DO jk = 1, nlev
!CDIR NODEP
        DO ie = jshift+1, jshift+nlen
          je = p_grf%idxlist_bdyintp_src_e(ie)
          jb = p_grf%blklist_bdyintp_src_e(ie)
#endif

          p_nh_save%vn(je,jk,jb) = p_nh_prog%vn(je,jk,jb)

        ENDDO
      ENDDO

    ENDDO
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE save_progvars

  !-------------------------------------------------------------------------
  !
  !
  !
  !>
  !! Computes the time tendencies of the prognostic variables needed for
  !! interpolation to the lateral boundaries of the nested domains
  !!
  !! @par Revision History
  !! Developed by Guenther Zaengl, DWD, 2010-02-10
  !!
  SUBROUTINE compute_tendencies (jg,n_new,n_now,n_new_rcf,n_now_rcf,&
    &                            rdt,rdt_mflx)


    INTEGER,  INTENT(IN) :: jg  ! domain ID
    ! Time levels from which tendencies are computed
    INTEGER,  INTENT(IN) ::  n_new,n_now
    ! Time levels from which tracer-tendencies are computed
    INTEGER,  INTENT(IN) ::  n_new_rcf,n_now_rcf

    ! Inverse value of time step needed for computing the tendencies
    REAL(wp), INTENT(IN) ::  rdt
    ! Inverse value of time step needed for computing the mass flux tendencies
    REAL(wp), INTENT(IN) ::  rdt_mflx


    ! local variables
    !
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx,       &
      ib, jb, ic, jc, ie, je, jk, jt, js, i_nchdom, nshift, jk_start, jshift
    INTEGER :: nlev, nlevp1           !< number of full and half levels
    INTEGER :: nproma_bdyintp, nblks_bdyintp, npromz_bdyintp, nlen, ntracer_bdyintp, nsubs

    REAL(wp) :: rdt_ubc, dthalf
    ! Switch to control if the child domain is vertically nested and therefore
    ! needs interpolation of upper boundary conditions
    LOGICAL :: l_child_vertnest

    TYPE(t_gridref_state), POINTER :: p_grf

    TYPE(t_nh_state), POINTER :: p_nh
    TYPE(t_nh_prog),  POINTER :: p_prog_now
    TYPE(t_nh_prog),  POINTER :: p_prog_new
    TYPE(t_nh_prog),  POINTER :: p_prog_now_rcf
    TYPE(t_nh_prog),  POINTER :: p_prog_new_rcf

    !-----------------------------------------------------------------------

    i_nchdom = MAX(1,p_patch(jg)%n_childdom)

    nsubs    = ndyn_substeps_var(jg)

    p_grf          => p_grf_state(jg)
    p_nh           => p_nh_state(jg)
    p_prog_now     => p_nh%prog(n_now)
    p_prog_new     => p_nh%prog(n_new)
    p_prog_now_rcf => p_nh%prog(n_now_rcf)
    p_prog_new_rcf => p_nh%prog(n_new_rcf)

    ! number of vertical levels
    nlev   = p_patch(jg)%nlev
    nlevp1 = p_patch(jg)%nlevp1

    IF (advection_config(jg)%iadv_tke == 1) THEN
      ntracer_bdyintp = ntracer-1
    ELSE
      ntracer_bdyintp = ntracer
    ENDIF

    ! determine if upper boundary interpolation is needed
    IF (lvert_nest .AND. (p_patch(jg)%nshift_child > 0)) THEN
      l_child_vertnest = .TRUE.
      nshift = p_patch(jg)%nshift_child + 1
    ELSE
      l_child_vertnest = .FALSE.
      nshift = 1
    ENDIF

    jk_start = nshift ! start index for tendency computation

    ! for dynamic nproma blocking
    nproma_bdyintp = MIN(nproma,256)

    rdt_ubc = rdt*REAL(nsubs,wp)/REAL(2*(MAX(1,nsubs-2)),wp)
    dthalf = 1._wp/(2._wp*rdt)

!$OMP PARALLEL PRIVATE(i_startblk,i_endblk,nblks_bdyintp,npromz_bdyintp)

    IF (l_child_vertnest) THEN ! Compute upper boundary condition and its time derivative for nested domain

      ! cell-based variables
      i_startblk = p_patch(jg)%cells%start_blk(grf_bdywidth_c+1,1)
      i_endblk   = p_patch(jg)%cells%end_blk(min_rlcell_int,i_nchdom)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk,js) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(p_patch(jg), jb, i_startblk, i_endblk, &
          i_startidx, i_endidx, grf_bdywidth_c+1, min_rlcell_int)

        p_nh%diag%dw_int         (:,jb,nsubs+1)         = 0._wp
        p_nh%diag%mflx_ic_int    (:,jb,nsubs+1:nsubs+2) = 0._wp
        p_nh%diag%dtheta_v_ic_int(:,jb,nsubs+1)         = 0._wp
        DO js = 1, nsubs
!DIR$ IVDEP
          DO jc = i_startidx, i_endidx
            p_nh%diag%dw_int(jc,jb,nsubs+1)          = p_nh%diag%dw_int(jc,jb,nsubs+1) +          &
              p_nh%diag%dw_int(jc,jb,js)
            p_nh%diag%mflx_ic_int(jc,jb,nsubs+1)     = p_nh%diag%mflx_ic_int(jc,jb,nsubs+1) +     &
              p_nh%diag%mflx_ic_int(jc,jb,js)
            p_nh%diag%dtheta_v_ic_int(jc,jb,nsubs+1) = p_nh%diag%dtheta_v_ic_int(jc,jb,nsubs+1) + &
              p_nh%diag%dtheta_v_ic_int(jc,jb,js)
          ENDDO
        ENDDO
        DO jc = i_startidx, i_endidx
          p_nh%diag%dw_int(jc,jb,nsubs+1)          = p_nh%diag%dw_int(jc,jb,nsubs+1)          / REAL(nsubs,wp)
          p_nh%diag%mflx_ic_int(jc,jb,nsubs+1)     = p_nh%diag%mflx_ic_int(jc,jb,nsubs+1)     / REAL(nsubs,wp)
          p_nh%diag%dtheta_v_ic_int(jc,jb,nsubs+1) = p_nh%diag%dtheta_v_ic_int(jc,jb,nsubs+1) / REAL(nsubs,wp)
        ENDDO
        IF (nsubs >= 3) THEN
!DIR$ IVDEP
          DO jc = i_startidx, i_endidx
            ! Compute time tendency of mass flux upper boundary condition to obtain second-order accuracy in time
            p_nh%diag%mflx_ic_int(jc,jb,nsubs+2) = (SUM(p_nh%diag%mflx_ic_int(jc,jb,nsubs-1:nsubs)) - &
              SUM(p_nh%diag%mflx_ic_int(jc,jb,1:2)))*rdt_ubc
            ! Shift time level of averaged field back to the beginning of the first dynamic substep
            p_nh%diag%mflx_ic_int(jc,jb,nsubs+1) = p_nh%diag%mflx_ic_int(jc,jb,nsubs+1) - &
              dthalf*p_nh%diag%mflx_ic_int(jc,jb,nsubs+2)
          ENDDO
        ENDIF

      ENDDO
!$OMP END DO

      ! edge-based variables
      i_startblk = p_patch(jg)%edges%start_blk(grf_bdywidth_e+1,1)
      i_endblk   = p_patch(jg)%edges%end_blk(min_rledge_int-2,i_nchdom)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,je,jk) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_e(p_patch(jg), jb, i_startblk, i_endblk, &
          i_startidx, i_endidx, grf_bdywidth_e+1, min_rledge_int-2)

!DIR$ IVDEP
        DO je = i_startidx, i_endidx
          p_nh%diag%dvn_ie_int(je,jb) = 0.5_wp*(p_nh%diag%dvn_ie_int(je,jb) + &
            p_nh%diag%vn_ie(je,nshift,jb) - p_nh%diag%vn_ie(je,nshift+1,jb))
        ENDDO

      ENDDO
!$OMP END DO

    ENDIF ! l_child_vertnest

    ! Computation of tendencies for lateral boundary interpolation

    ! cell-based variables

    ! parameters for dynamic nproma blocking
    nblks_bdyintp  = INT(p_grf%npoints_bdyintp_src_c/nproma_bdyintp)
    npromz_bdyintp = MOD(p_grf%npoints_bdyintp_src_c,nproma_bdyintp)
    IF (npromz_bdyintp > 0) THEN
      nblks_bdyintp = nblks_bdyintp + 1
    ELSE
      npromz_bdyintp = nproma_bdyintp
    ENDIF

!$OMP DO PRIVATE(ib,jb,nlen,ic,jc,jk,jt,jshift) ICON_OMP_DEFAULT_SCHEDULE
    DO ib = 1, nblks_bdyintp
      IF (ib == nblks_bdyintp) THEN
        nlen = npromz_bdyintp
      ELSE
        nlen = nproma_bdyintp
      ENDIF
      jshift = (ib-1)*nproma_bdyintp

#ifdef __LOOP_EXCHANGE
      DO ic = jshift+1, jshift+nlen
        jc = p_grf%idxlist_bdyintp_src_c(ic)
        jb = p_grf%blklist_bdyintp_src_c(ic)
!DIR$ IVDEP
        DO jk = jk_start, nlev
#else
      DO jk = jk_start, nlev
!CDIR NODEP
        DO ic = jshift+1, jshift+nlen
          jc = p_grf%idxlist_bdyintp_src_c(ic)
          jb = p_grf%blklist_bdyintp_src_c(ic)
#endif

          p_nh%diag%grf_tend_rho(jc,jk,jb) = &
            ( p_prog_new%rho(jc,jk,jb) - p_prog_now%rho(jc,jk,jb) )*rdt

          p_nh%diag%grf_tend_thv(jc,jk,jb) = &
            ( p_prog_new%theta_v(jc,jk,jb) - p_prog_now%theta_v(jc,jk,jb) )*rdt

          ! the div field carries perturbation density for use in SR boundary_interpolation
          p_nh%diag%div(jc,jk,jb) = &
            p_prog_now%rho(jc,jk,jb) - p_nh%metrics%rho_ref_mc(jc,jk,jb)

          ! the dpres_mc field carries perturbation potential temperature for use in SR boundary_interpolation
          p_nh%diag%dpres_mc(jc,jk,jb) = &
            p_prog_now%theta_v(jc,jk,jb) - p_nh%metrics%theta_ref_mc(jc,jk,jb)

          p_nh%diag%grf_tend_w(jc,jk,jb) = &
            ( p_prog_new%w(jc,jk,jb) - p_prog_now%w(jc,jk,jb) )*rdt
        ENDDO
      ENDDO

      DO ic = jshift+1, jshift+nlen
        jc = p_grf%idxlist_bdyintp_src_c(ic)
        jb = p_grf%blklist_bdyintp_src_c(ic)
        p_nh%diag%grf_tend_w(jc,nlevp1,jb) = &
          ( p_prog_new%w(jc,nlevp1,jb) - p_prog_now%w(jc,nlevp1,jb) )*rdt
      ENDDO

      IF (ltransport) THEN

#ifdef __LOOP_EXCHANGE
      DO ic = jshift+1, jshift+nlen
        jc = p_grf%idxlist_bdyintp_src_c(ic)
        jb = p_grf%blklist_bdyintp_src_c(ic)
        DO jt = 1,ntracer_bdyintp
!DIR$ IVDEP
          DO jk = jk_start, nlev
#else
      DO jt = 1,ntracer_bdyintp
        DO jk = jk_start, nlev
!CDIR NODEP
          DO ic = jshift+1, jshift+nlen
            jc = p_grf%idxlist_bdyintp_src_c(ic)
            jb = p_grf%blklist_bdyintp_src_c(ic)
#endif

              p_nh%diag%grf_tend_tracer(jc,jk,jb,jt) =                 &
                &            ( p_prog_new_rcf%tracer(jc,jk,jb,jt)      &
                &            -  p_prog_now_rcf%tracer(jc,jk,jb,jt) )*rdt
            ENDDO
          ENDDO
        ENDDO

      ENDIF

    ENDDO
!$OMP END DO

    ! edge-based variables

    ! parameters for dynamic nproma blocking
    nblks_bdyintp  = INT(p_grf%npoints_bdyintp_src_e/nproma_bdyintp)
    npromz_bdyintp = MOD(p_grf%npoints_bdyintp_src_e,nproma_bdyintp)
    IF (npromz_bdyintp > 0) THEN
      nblks_bdyintp = nblks_bdyintp + 1
    ELSE
      npromz_bdyintp = nproma_bdyintp
    ENDIF

!$OMP DO PRIVATE(ib,jb,nlen,ie,je,jk,jshift) ICON_OMP_DEFAULT_SCHEDULE
    DO ib = 1, nblks_bdyintp
      IF (ib == nblks_bdyintp) THEN
        nlen = npromz_bdyintp
      ELSE
        nlen = nproma_bdyintp
      ENDIF
      jshift = (ib-1)*nproma_bdyintp

#ifdef __LOOP_EXCHANGE
      DO ie = jshift+1, jshift+nlen
        je = p_grf%idxlist_bdyintp_src_e(ie)
        jb = p_grf%blklist_bdyintp_src_e(ie)
!DIR$ IVDEP
        DO jk = jk_start, nlev
#else
      DO jk = jk_start, nlev
!CDIR NODEP
        DO ie = jshift+1, jshift+nlen
          je = p_grf%idxlist_bdyintp_src_e(ie)
          jb = p_grf%blklist_bdyintp_src_e(ie)
#endif
          p_nh%diag%grf_tend_vn(je,jk,jb) = &
            ( p_prog_new%vn(je,jk,jb) - p_prog_now%vn(je,jk,jb) )*rdt
          p_nh%diag%grf_tend_mflx(je,jk,jb) = &
            ( p_nh%diag%mass_fl_e(je,jk,jb) - p_nh%diag%mass_fl_e_sv(je,jk,jb) )*rdt_mflx
        ENDDO
      ENDDO

    ENDDO
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE compute_tendencies

  !-------------------------------------------------------------------------
  !
  !>
  !! Interpolates time tendencies of prognostic variables to the lateral boundary
  !! of a refined mesh
  !!
  !! @par Revision History
  !! Developed  by Guenther Zaengl, DWD, 2008-07-10
  !!
  SUBROUTINE boundary_interpolation (jg,jgc,ntp_dyn,ntc_dyn,ntp_tr,ntc_tr, &
    mass_flx_p,mass_flx_c)

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_nh_nest_utilities:boundary_interpolation'


    INTEGER, INTENT(IN)     :: jg, jgc      ! domain ID of parent and child grid

    ! Parent and child time levels for dynamical variables and tracers
    INTEGER, INTENT(IN)     :: ntp_dyn, ntc_dyn, ntp_tr, ntc_tr

    ! Mass fluxes at parent and child level
    REAL(wp), DIMENSION(:,:,:), INTENT(INOUT) :: mass_flx_p, mass_flx_c

    ! local variables

    TYPE(t_nh_diag), POINTER           :: p_diagp => NULL()
    TYPE(t_nh_diag), POINTER           :: p_diagc  => NULL()
    TYPE(t_nh_prog), POINTER           :: p_nhp_dyn => NULL()
    TYPE(t_nh_prog), POINTER           :: p_nhc_dyn => NULL()
    TYPE(t_nh_prog), POINTER           :: p_nhp_tr  => NULL()
    TYPE(t_nh_prog), POINTER           :: p_nhc_tr  => NULL()
    TYPE(t_patch), POINTER             :: p_pp => NULL()
    TYPE(t_patch), POINTER             :: p_pc => NULL()
    TYPE(t_gridref_state), POINTER     :: p_grf => NULL()
    TYPE(t_grid_cells), POINTER        :: p_gcp => NULL()


    INTEGER :: i_startblk              ! start block
    INTEGER :: i_endblk                ! end index
    INTEGER :: i_startidx              ! start index
    INTEGER :: i_endidx                ! end index

    INTEGER :: jb, jc, jk, jt, ic      ! loop indices

    INTEGER :: nlev_c                  ! number of full levels (child domain)

    INTEGER :: i_chidx, i_nchdom, i_sbc, i_ebc
    INTEGER :: ntracer_bdyintp, nsubs

    REAL(wp) :: aux3dp(nproma,ntracer+4,p_patch(jg)%nblks_c), &
      aux3dc(nproma,ntracer+4,p_patch(jgc)%nblks_c), &
      theta_prc(nproma,p_patch(jgc)%nlev,p_patch(jgc)%nblks_c), &
      rho_prc(nproma,p_patch(jgc)%nlev,p_patch(jgc)%nblks_c)

    ! Switch to control if the child domain is vertically nested and therefore
    ! needs interpolation of upper boundary conditions
    LOGICAL :: l_child_vertnest

    LOGICAL :: l_limit(2*ntracer)

    !$ INTEGER :: num_threads_omp, omp_get_max_threads
    !-----------------------------------------------------------------------

    IF (msg_level >= 10) THEN
      WRITE(message_text,'(a,i2,a,i2)') '========= Interpolate:',jg,' =>',jgc
      CALL message(TRIM(routine),message_text)
    ENDIF

    !$  num_threads_omp = omp_get_max_threads()

    p_diagp       => p_nh_state(jg)%diag
    p_diagc       => p_nh_state(jgc)%diag
    p_nhp_dyn     => p_nh_state(jg)%prog(ntp_dyn)
    p_nhc_dyn     => p_nh_state(jgc)%prog(ntc_dyn)
    p_nhp_tr      => p_nh_state(jg)%prog(ntp_tr)
    p_nhc_tr      => p_nh_state(jgc)%prog(ntc_tr)
    p_grf         => p_grf_state(jg)
    p_pp          => p_patch(jg)
    p_pc          => p_patch(jgc)
    p_gcp         => p_patch(jg)%cells


    i_chidx = p_patch(jgc)%parent_child_index
    i_nchdom   = MAX(1,p_pc%n_childdom)

    nsubs    = ndyn_substeps_var(jg)

    ! number of vertical levels (child domain)
    nlev_c   = p_pc%nlev

    ! determine if upper boundary interpolation is needed
    IF (lvert_nest .AND. (p_pp%nshift_child > 0)) THEN
      l_child_vertnest = .TRUE.
    ELSE
      l_child_vertnest = .FALSE.
    ENDIF

    ! exclude TKE from boundary interpolation if it is only vertically advected
    IF (advection_config(jg)%iadv_tke == 1) THEN
      ntracer_bdyintp = ntracer-1
    ELSE
      ntracer_bdyintp = ntracer
    ENDIF

    ! Perform interpolations to upper nest boundary needed for vertical nesting
    IF (l_child_vertnest) THEN

      IF (p_test_run) THEN
        aux3dp = 0._wp
        aux3dc = 0._wp
      ENDIF

      CALL sync_patch_array(SYNC_E,p_pp,p_diagp%dvn_ie_int)

      CALL interpol_vec_ubc (p_pp, p_pc, p_grf%p_dom(i_chidx), &
        p_diagp%dvn_ie_int, p_diagc%dvn_ie_ubc)

      ! Start and end blocks for which interpolation is needed
      i_startblk = p_pp%cells%start_block(0)
      i_endblk   = p_pp%cells%end_block(min_rlcell_int)

      ! For back-copying at child level
      i_nchdom   = MAX(1,p_pc%n_childdom)
      i_sbc      = p_pc%cells%start_blk(grf_nudge_start_c,1)
      i_ebc      = p_pc%cells%end_blk(min_rlcell_int,i_nchdom)

      IF (ltransport) THEN

!$OMP PARALLEL DO PRIVATE(jb,i_startidx,i_endidx,jc,jt) ICON_OMP_DEFAULT_SCHEDULE
        DO jb =  i_startblk, i_endblk

          CALL get_indices_c(p_pp, jb, i_startblk, i_endblk, i_startidx, i_endidx, &
            0, min_rlcell_int)

          DO jc = i_startidx, i_endidx
            aux3dp(jc,1,jb)   = p_diagp%dw_int(jc,jb,nsubs+1)
            aux3dp(jc,2:3,jb) = p_diagp%mflx_ic_int(jc,jb,nsubs+1:nsubs+2)
            aux3dp(jc,4,jb)   = p_diagp%dtheta_v_ic_int(jc,jb,nsubs+1)
          ENDDO
          DO jt = 1, ntracer
            DO jc = i_startidx, i_endidx
              aux3dp(jc,4+jt,jb) = p_diagp%q_int(jc,jb,jt)
            ENDDO
          ENDDO
        ENDDO
!$OMP END PARALLEL DO

        CALL sync_patch_array(SYNC_C,p_pp,aux3dp)

        CALL interpol_scal_ubc (p_pp, p_pc, p_grf%p_dom(i_chidx),  &
          ntracer+4, aux3dp, aux3dc)

!$OMP PARALLEL DO PRIVATE(jb,i_startidx,i_endidx,jc,jt) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = i_sbc, i_ebc

          CALL get_indices_c(p_pc, jb, i_sbc, i_ebc, i_startidx, i_endidx, &
            grf_nudge_start_c, min_rlcell_int)

          DO jc = i_startidx, i_endidx
            p_diagc%dw_ubc(jc,jb)          = aux3dc(jc,1,jb)
            p_diagc%mflx_ic_ubc(jc,jb,1:2) = aux3dc(jc,2:3,jb)
            p_diagc%dtheta_v_ic_ubc(jc,jb) = aux3dc(jc,4,jb)
          ENDDO
          DO jt = 1, ntracer
            DO jc = i_startidx, i_endidx
              p_diagc%q_ubc(jc,jb,jt) = aux3dc(jc,jt+4,jb)
            ENDDO
          ENDDO
        ENDDO
!$OMP END PARALLEL DO

      ELSE

!$OMP PARALLEL DO PRIVATE(jb,i_startidx,i_endidx,jc) ICON_OMP_DEFAULT_SCHEDULE
        DO jb =  i_startblk, i_endblk

          CALL get_indices_c(p_pp, jb, i_startblk, i_endblk, i_startidx, i_endidx, &
            0, min_rlcell_int)

          DO jc = i_startidx, i_endidx
            aux3dp(jc,1,jb)   = p_diagp%dw_int(jc,jb,nsubs+1)
            aux3dp(jc,2:3,jb) = p_diagp%mflx_ic_int(jc,jb,nsubs+1:nsubs+2)
            aux3dp(jc,4,jb)   = p_diagp%dtheta_v_ic_int(jc,jb,nsubs+1)
          ENDDO
        ENDDO
!$OMP END PARALLEL DO

        CALL sync_patch_array(SYNC_C,p_pp,aux3dp)

        CALL interpol_scal_ubc(p_pp, p_pc, p_grf%p_dom(i_chidx), 4, aux3dp, aux3dc)

!$OMP PARALLEL DO PRIVATE(jb,i_startidx,i_endidx,jc) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = i_sbc, i_ebc

          CALL get_indices_c(p_pc, jb, i_sbc, i_ebc, i_startidx, i_endidx, &
            grf_nudge_start_c, min_rlcell_int)

          DO jc = i_startidx, i_endidx
            p_diagc%dw_ubc(jc,jb)          = aux3dc(jc,1,jb)
            p_diagc%mflx_ic_ubc(jc,jb,1:2) = aux3dc(jc,2:3,jb)
            p_diagc%dtheta_v_ic_ubc(jc,jb) = aux3dc(jc,4,jb)
          ENDDO
        ENDDO
!$OMP END PARALLEL DO

      ENDIF

    ENDIF

    ! Lateral boundary interpolation of cell-based dynamical variables

    IF (grf_intmethod_c == 1) THEN ! tendency copying for all cell-based variables

      CALL exchange_data(p_pc%comm_pat_interpolation_c, &
        RECV=p_diagc%grf_tend_rho,     &
        SEND=p_diagp%grf_tend_rho)

      CALL exchange_data(p_pc%comm_pat_interpolation_c, &
        RECV=p_diagc%grf_tend_thv,     &
        SEND=p_diagp%grf_tend_thv)

      ! exchange_data should also work for w because it determines the
      ! vertical dimension with UBOUND
      CALL exchange_data(p_pc%comm_pat_interpolation_c, &
        RECV=p_diagc%grf_tend_w,       &
        SEND=p_diagp%grf_tend_w)

      ! grf_intmethod_c = 2, use gradient at cell center for interpolation
    ELSE IF (grf_intmethod_c == 2) THEN

      ! Interpolation of temporal tendencies, full w, perturbation density (stored in div)
      !  and perturbationvirtual potential temperature (stored in dpres_mc)
      CALL interpol_scal_grf (p_pp, p_pc, p_grf%p_dom(i_chidx), 6,         &
        p_diagp%grf_tend_rho, p_diagc%grf_tend_rho,  &
        p_diagp%grf_tend_thv, p_diagc%grf_tend_thv,  &
        p_diagp%grf_tend_w,   p_diagc%grf_tend_w,    &
        p_nhp_dyn%w,          p_nhc_dyn%w,           &
        p_nh_state(jg)%diag%div, rho_prc,            &
        p_nh_state(jg)%diag%dpres_mc, theta_prc      )

      ! Start and end blocks for which interpolation is needed
      i_startblk = p_pc%cells%start_blk(1,1)
      i_endblk   = p_pc%cells%end_blk(grf_bdywidth_c,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk) ICON_OMP_DEFAULT_SCHEDULE
      DO jb =  i_startblk, i_endblk

        CALL get_indices_c(p_pc, jb, i_startblk, i_endblk, i_startidx, i_endidx, &
          1, grf_bdywidth_c)

        DO jk = 1, nlev_c
!DIR$ IVDEP
          DO jc = i_startidx, i_endidx
            p_nhc_dyn%rho(jc,jk,jb) = rho_prc(jc,jk,jb) + &
              p_nh_state(jgc)%metrics%rho_ref_mc(jc,jk,jb)
            p_nhc_dyn%theta_v(jc,jk,jb) = theta_prc(jc,jk,jb) + &
              p_nh_state(jgc)%metrics%theta_ref_mc(jc,jk,jb)
          ENDDO
        ENDDO
      ENDDO
!$OMP END DO

      ! The following index list contains the halo points of the lateral boundary
      ! cells. These have to be copied as well in order for rho and theta to be
      ! synchronized.
!$OMP DO PRIVATE(ic,jb,jc,jk)
      DO ic = 1, p_nh_state(jgc)%metrics%bdy_halo_c_dim

        jb = p_nh_state(jgc)%metrics%bdy_halo_c_blk(ic)
        jc = p_nh_state(jgc)%metrics%bdy_halo_c_idx(ic)
!DIR$ IVDEP
        DO jk = 1, nlev_c
          p_nhc_dyn%rho(jc,jk,jb) = rho_prc(jc,jk,jb) + &
            p_nh_state(jgc)%metrics%rho_ref_mc(jc,jk,jb)
          p_nhc_dyn%theta_v(jc,jk,jb) = theta_prc(jc,jk,jb) + &
            p_nh_state(jgc)%metrics%theta_ref_mc(jc,jk,jb)

        ENDDO
      ENDDO
!$OMP END DO
!$OMP END PARALLEL
    ENDIF

    ! Lateral boundary interpolation of cell based tracer variables
    IF (ltransport .AND. grf_intmethod_ct == 1) THEN

      ! Start and end blocks for which interpolation is needed
      i_startblk = p_gcp%start_blk(grf_bdyintp_start_c,i_chidx)
      i_endblk   = p_gcp%end_blk(grf_bdyintp_end_c,i_chidx)

      CALL exchange_data_mult(p_pc%comm_pat_interpolation_c, ntracer_bdyintp, ntracer_bdyintp*nlev_c, &
        RECV4D=p_diagc%grf_tend_tracer(:,:,:,1:ntracer_bdyintp),                                      &
        SEND4D=p_diagp%grf_tend_tracer(:,:,:,1:ntracer_bdyintp))


    ELSE IF (ltransport .AND. grf_intmethod_ct == 2) THEN

      ! Apply positive definite limiter on full tracer fields but not on tendencies
      l_limit(1:ntracer_bdyintp) = .FALSE.
      l_limit(ntracer_bdyintp+1:2*ntracer_bdyintp) = .TRUE.

      CALL interpol_scal_grf ( p_pp, p_pc, p_grf%p_dom(i_chidx), 2*ntracer_bdyintp, &
        f4din1 =  p_diagp%grf_tend_tracer(:,:,:,1:ntracer_bdyintp),                 &
        f4dout1 = p_diagc%grf_tend_tracer(:,:,:,1:ntracer_bdyintp),                 &
        f4din2  = p_nhp_tr%tracer(:,:,:,1:ntracer_bdyintp),                         &
        f4dout2 = p_nhc_tr%tracer(:,:,:,1:ntracer_bdyintp), llimit_nneg=l_limit     )

    ENDIF

    IF (ltransport .AND. iprog_aero >= 1) THEN
     CALL interpol_scal_grf (p_pp, p_pc, p_grf%p_dom(i_chidx), 1, prm_diag(jg)%aerosol,    &
                             prm_diag(jgc)%aerosol, llimit_nneg=(/.TRUE./), lnoshift=.TRUE.)
    ENDIF

    ! Lateral boundary interpolation of edge-based variables  (velocity components)
    IF (grf_intmethod_e == 1 .OR. grf_intmethod_e == 2) THEN

      CALL interpol_vec_grf (p_pp, p_pc, p_grf%p_dom(i_chidx), p_diagp%grf_tend_vn, p_diagc%grf_tend_vn)

    ELSE IF (grf_intmethod_e == 3 .OR. grf_intmethod_e == 4) THEN

      CALL interpol2_vec_grf (p_pp, p_pc, p_grf%p_dom(i_chidx), 1, p_diagp%grf_tend_vn, p_diagc%grf_tend_vn)

    ELSE IF (grf_intmethod_e == 5 .OR. grf_intmethod_e == 6) THEN

      CALL interpol2_vec_grf (p_pp, p_pc, p_grf%p_dom(i_chidx), 3,        &
        p_diagp%grf_tend_vn, p_diagc%grf_tend_vn, mass_flx_p, mass_flx_c, &
        p_diagp%grf_tend_mflx, p_diagc%grf_tend_mflx )

    ENDIF


  END SUBROUTINE boundary_interpolation


  !>
  !! This routine prepares boundary nudging for use with 1-way nesting.
  !!
  !! The following steps are executed:
  !! 1. Mapping of coarse-grid prognostic variables to intermediate grid sharing
  !!    the domain decomposition and vertical dimension with the child grid
  !! 2. Computation of differences between parent-grid values and averaged child-grid
  !!    variables
  !! 3. Interpolation of difference fields to the child grid
  !!
  !! @par Revision History
  !! Developed  by Guenther Zaengl, DWD, 2010-06-18
  SUBROUTINE prep_bdy_nudging(jgp, jg)

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_nh_nest_utilities:prep_bdy_nudging'


    INTEGER, INTENT(IN) :: jg   ! child grid level
    INTEGER, INTENT(IN) :: jgp  ! parent grid level


    ! local variables
    !
    TYPE(t_nh_prog),    POINTER     :: p_parent_prog => NULL()
    TYPE(t_nh_prog),    POINTER     :: p_child_prog  => NULL()
    TYPE(t_nh_prog),    POINTER     :: p_parent_prog_rcf => NULL()
    TYPE(t_nh_prog),    POINTER     :: p_child_prog_rcf  => NULL()
    TYPE(t_nh_diag),    POINTER     :: p_diag        => NULL()
    TYPE(t_grid_cells), POINTER     :: p_gcp => NULL()
    TYPE(t_grid_edges), POINTER     :: p_gep => NULL()
    TYPE(t_gridref_state), POINTER  :: p_grf => NULL()
    TYPE(t_gridref_state), POINTER  :: p_grfc => NULL()
    TYPE(t_int_state), POINTER      :: p_int => NULL()
    TYPE(t_patch),      POINTER     :: p_pp => NULL()
    TYPE(t_patch),      POINTER     :: p_pc => NULL()

    ! Indices
    INTEGER :: jb, jc, jk, jt, je, js, i_chidx, i_startblk, i_endblk, &
      i_startidx, i_endidx, istartblk_c, istartblk_e
    INTEGER :: nlev_c, nlev_p
    INTEGER :: nshift      !< difference between upper boundary of parent or feedback-parent
                           !< domain and upper boundary of child domain (in terms
                           !< of vertical levels)
    INTEGER :: ntracer_nudge !< number of tracers to be nudged

    ! Local arrays for interpolated parent-grid values, and difference fields. These have
    ! to be allocatable because their dimensions differ between MPI and non-MPI runs
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)   :: parent_thv, diff_thv
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)   :: parent_rho, diff_rho
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)   :: parent_vn , diff_vn
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)   :: parent_w  , diff_w
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:,:) :: parent_tr , diff_tr

    INTEGER, DIMENSION(:,:,:), POINTER :: iidx, iblk, ieidx, ieblk
    LOGICAL :: l_parallel
    REAL(wp), DIMENSION(:,:,:), POINTER :: p_fbkwgt, p_fbkwgt_tr, p_fbkwgt_v
    !-----------------------------------------------------------------------

    IF (msg_level >= 10) THEN
      WRITE(message_text,'(a,i2,a,i2)') '1-way nesting: == Boundary nudging:',jg
      CALL message(TRIM(routine),message_text)
    ENDIF

    l_parallel = my_process_is_mpi_parallel()

    IF (iforcing > 1) THEN  ! tracers represent moisture variables
      ntracer_nudge = 3     ! take only QV, QC and QI - nudging precip variables does not make much sense
    ELSE
      ntracer_nudge = ntracer
    ENDIF

    p_parent_prog     => p_nh_state(jgp)%prog(nsav1(jgp))
    p_child_prog      => p_nh_state(jg)%prog(nnow(jg))
    p_parent_prog_rcf => p_nh_state(jgp)%prog(nnow_rcf(jgp))
    p_child_prog_rcf  => p_nh_state(jg)%prog(nnow_rcf(jg))
    p_diag            => p_nh_state(jg)%diag
    p_grfc            => p_grf_state(jg)
    p_pc              => p_patch(jg)

    p_grf => p_grf_state_local_parent(jg)
    p_int => p_int_state_local_parent(jg)
    p_gcp => p_patch_local_parent(jg)%cells
    p_gep => p_patch_local_parent(jg)%edges
    p_pp  => p_patch_local_parent(jg)

    i_chidx  = p_pc%parent_child_index

    ! number of full levels of child domain
    nlev_c   = p_pc%nlev

    ! number of full/half levels of parent domain
    nlev_p   = p_pp%nlev

    ! shift between upper model boundaries
    nshift = p_pc%nshift
    js     = nshift

    ! Please note: In the parallel case
    ! - lower bound must be 1 due to synchronization calls
    ! - upper bound must be nblks_c/e to include halo cells/edges
    ! - this doesn't cost extra memory since p_patch_local_parent
    !   only includes the cells/edges really needed

    ! Value of i_startblk needed for subroutine call for parent-to-child interpolation
    istartblk_c = 1

    ALLOCATE(parent_thv  (nproma, nlev_p, p_patch_local_parent(jg)%nblks_c),  &
      diff_thv    (nproma, nlev_c, p_patch_local_parent(jg)%nblks_c),  &
      parent_rho  (nproma, nlev_p, p_patch_local_parent(jg)%nblks_c),  &
      diff_rho    (nproma, nlev_c, p_patch_local_parent(jg)%nblks_c),  &
      parent_w    (nproma, nlev_p, p_patch_local_parent(jg)%nblks_c),  &
      diff_w      (nproma, nlev_c+1, p_patch_local_parent(jg)%nblks_c) )

    IF(ltransport) &
      ALLOCATE(parent_tr(nproma, nlev_p, p_patch_local_parent(jg)%nblks_c, ntracer_nudge),&
      &        diff_tr  (nproma, nlev_c, p_patch_local_parent(jg)%nblks_c, ntracer_nudge) )

    ! Value of i_startblk needed for subroutine call for parent-to-child interpolation
    istartblk_e = 1

    ALLOCATE(parent_vn  (nproma, nlev_p, p_patch_local_parent(jg)%nblks_e), &
      &      diff_vn    (nproma, nlev_c, p_patch_local_parent(jg)%nblks_e)  )

    ! Set pointers to index and coefficient fields for cell-based variables
    iidx  => p_gcp%child_idx
    iblk  => p_gcp%child_blk
    ieidx => p_gep%child_idx
    ieblk => p_gep%child_blk

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
    p_fbkwgt_v  => p_grf%fbk_wgt_e

    ! 1st step: Copy prognostic variables from parent grid to fields on feedback-parent grid
    ! (trivial without MPI parallelization, but communication call needed for MPI)

    CALL exchange_data_mult(p_pp%comm_pat_glb_to_loc_c, 3, 3*nlev_p, &
      RECV1=parent_rho, SEND1=p_parent_prog%rho,         &
      RECV2=parent_thv, SEND2=p_parent_prog%theta_v,     &
      RECV3=parent_w,   SEND3=p_parent_prog%w            )

    CALL exchange_data(p_pp%comm_pat_glb_to_loc_e,     &
      RECV=parent_vn, SEND=p_parent_prog%vn)

    IF (ltransport) &
      CALL exchange_data_mult(p_pp%comm_pat_glb_to_loc_c, ntracer_nudge, ntracer_nudge*nlev_p, &
      RECV4D=parent_tr, SEND4D=p_parent_prog_rcf%tracer(:,:,:,1:ntracer_nudge))

    ! 2nd step: perform feedback from refined grid to intermediate grid and compute differences

!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)

    ! a) cell-based variables
    ! Start/End block in the parent domain
    i_startblk = p_gcp%start_blk(grf_nudgintp_start_c+1,i_chidx)
    i_endblk   = p_gcp%end_blk(min_rlcell_int,i_chidx)


!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk,jt) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_pp, jb, i_startblk, i_endblk, i_startidx, i_endidx, &
        grf_nudgintp_start_c+1, min_rlcell_int)

      ! initialize diff_w at surface with zero
      diff_w(:,nlev_c+1,jb) = 0._wp

#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jk = 1, nlev_c
#else
!CDIR UNROLL=5
      DO jk = 1, nlev_c
        DO jc = i_startidx, i_endidx
#endif

          diff_thv(jc,jk,jb) = parent_thv(jc,jk+js,jb) - (                           &
            p_child_prog%theta_v(iidx(jc,jb,1),jk,iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
            p_child_prog%theta_v(iidx(jc,jb,2),jk,iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
            p_child_prog%theta_v(iidx(jc,jb,3),jk,iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
            p_child_prog%theta_v(iidx(jc,jb,4),jk,iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)   )

          diff_rho(jc,jk,jb) = parent_rho(jc,jk+js,jb) - (                       &
            p_child_prog%rho(iidx(jc,jb,1),jk,iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
            p_child_prog%rho(iidx(jc,jb,2),jk,iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
            p_child_prog%rho(iidx(jc,jb,3),jk,iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
            p_child_prog%rho(iidx(jc,jb,4),jk,iblk(jc,jb,4))*p_fbkwgt(jc,jb,4))+ &
            p_nh_state(jg)%metrics%rho_ref_corr(jc,jk,jb)

          diff_w(jc,jk,jb) = parent_w(jc,jk+js,jb) - (                         &
            p_child_prog%w(iidx(jc,jb,1),jk,iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
            p_child_prog%w(iidx(jc,jb,2),jk,iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
            p_child_prog%w(iidx(jc,jb,3),jk,iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
            p_child_prog%w(iidx(jc,jb,4),jk,iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)   )

        ENDDO
      ENDDO

      ! Tracers
      IF (ltransport) THEN

        DO jt = 1, ntracer_nudge

#ifdef __LOOP_EXCHANGE
          DO jc = i_startidx, i_endidx
            DO jk = 1, nlev_c
#else
!CDIR UNROLL=8
          DO jk = 1, nlev_c
            DO jc = i_startidx, i_endidx
#endif

              diff_tr(jc,jk,jb,jt) = parent_tr(jc,jk+js,jb,jt) - (                                &
                p_child_prog_rcf%tracer(iidx(jc,jb,1),jk,iblk(jc,jb,1),jt)*p_fbkwgt_tr(jc,jb,1) + &
                p_child_prog_rcf%tracer(iidx(jc,jb,2),jk,iblk(jc,jb,2),jt)*p_fbkwgt_tr(jc,jb,2) + &
                p_child_prog_rcf%tracer(iidx(jc,jb,3),jk,iblk(jc,jb,3),jt)*p_fbkwgt_tr(jc,jb,3) + &
                p_child_prog_rcf%tracer(iidx(jc,jb,4),jk,iblk(jc,jb,4),jt)*p_fbkwgt_tr(jc,jb,4)   )
            ENDDO
          ENDDO

        ENDDO

      ENDIF

    ENDDO
!$OMP END DO

    ! b) velocity
    ! Start/End block in the parent domain
    i_startblk = p_gep%start_blk(grf_nudgintp_start_e+2,i_chidx)
    i_endblk   = p_gep%end_blk(min_rledge_int,i_chidx)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,je,jk) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(p_pp, jb, i_startblk, i_endblk, i_startidx, i_endidx, &
        grf_nudgintp_start_e+2, min_rledge_int)


#ifdef __LOOP_EXCHANGE
      DO je = i_startidx, i_endidx
        DO jk = 1, nlev_c
#else
!CDIR UNROLL=5
      DO jk = 1, nlev_c
        DO je = i_startidx, i_endidx
#endif

          diff_vn(je,jk,jb) = parent_vn(je,jk+js,jb) - (                            &
            p_child_prog%vn(ieidx(je,jb,1),jk,ieblk(je,jb,1))*p_fbkwgt_v(je,jb,1) + &
            p_child_prog%vn(ieidx(je,jb,2),jk,ieblk(je,jb,2))*p_fbkwgt_v(je,jb,2)   )

        ENDDO
      ENDDO

    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL


    ! Interpolate differences to child grid; the differences are stored in the grf_tend fields

    ! interpol_vec/scal_nudging needs all fields with full boundaries, so the arrays
    ! have to be sync'd before calling these routines.
    ! Please note that we cannot use sync_patch_array here (comparing parallel/non parallel results)
    ! since the arrays don't start with lower bound 1 in the non parallel case!

    ! Synchronization is needed after the interpolation step because the nudging tendencies are applied outside 
    ! the dynamical core. This is needed for the scalars for reasons of mass consistency, but is also done for the
    ! wind tendencies because this turns out to improve noise filtering

    IF(l_parallel) CALL exchange_data(p_pp%comm_pat_e, diff_vn)
    CALL interpol_vec_nudging (p_pp, p_pc, p_int, p_grf%p_dom(i_chidx), p_grfc,    &
      &                        i_chidx, 0, istartblk_e, diff_vn,p_diag%grf_tend_vn )
    CALL sync_patch_array(SYNC_E,p_pc,p_diag%grf_tend_vn)

    IF(l_parallel) CALL exchange_data_mult(p_pp%comm_pat_c, 3, 3*nlev_c+1, &
      recv1=diff_thv, recv2=diff_rho, recv3=diff_w)
    CALL interpol_scal_nudging (p_pp, p_int, p_grf%p_dom(i_chidx), i_chidx, 0, 3, istartblk_c, &
      &                         f3din1=diff_thv, f3dout1=p_diag%grf_tend_thv,                  &
      &                         f3din2=diff_rho, f3dout2=p_diag%grf_tend_rho,                  &
      &                         f3din3=diff_w,   f3dout3=p_diag%grf_tend_w                     )
    CALL sync_patch_array_mult(SYNC_C,p_pc,3,p_diag%grf_tend_thv,p_diag%grf_tend_rho,  &
      p_diag%grf_tend_w)

    IF (ltransport) THEN
      IF(l_parallel) CALL exchange_data_mult(p_pp%comm_pat_c, ntracer_nudge, ntracer_nudge*nlev_c, recv4d=diff_tr)

      CALL interpol_scal_nudging (p_pp, p_int, p_grf%p_dom(i_chidx), i_chidx,          &
        &                         0, ntracer_nudge, istartblk_c, f4din=diff_tr,        &
        &                         f4dout=p_diag%grf_tend_tracer(:,:,:,1:ntracer_nudge) )
      CALL sync_patch_array_mult(SYNC_C,p_pc,ntracer_nudge,f4din=p_diag%grf_tend_tracer(:,:,:,1:ntracer_nudge))
    ENDIF

    DEALLOCATE(parent_thv, diff_thv, parent_rho, diff_rho, parent_w, diff_w, parent_vn, diff_vn)

    IF(ltransport) DEALLOCATE(parent_tr, diff_tr)


  END SUBROUTINE prep_bdy_nudging

  !>
  !! This routine prepares boundary nudging for density only (for use with 2-way nesting)
  !!
  !! The following steps are executed:
  !! 1. Mapping of coarse-grid density to intermediate grid sharing
  !!    the domain decomposition and vertical dimension with the child grid
  !! 2. Computation of differences between parent-grid values and averaged child-grid
  !!    values
  !! 3. Interpolation of difference fields to the child grid
  !!
  !! @par Revision History
  !! Developed  by Guenther Zaengl, DWD, 2011-12-08
  SUBROUTINE prep_rho_bdy_nudging(jgp, jg)


    INTEGER, INTENT(IN) :: jg   ! child grid level
    INTEGER, INTENT(IN) :: jgp  ! parent grid level


    ! local variables

    TYPE(t_nh_prog),    POINTER     :: p_parent_prog => NULL()
    TYPE(t_nh_prog),    POINTER     :: p_child_prog  => NULL()
    TYPE(t_nh_diag),    POINTER     :: p_diag        => NULL()
    TYPE(t_grid_cells), POINTER     :: p_gcp => NULL()
    TYPE(t_gridref_state), POINTER  :: p_grf => NULL()
    TYPE(t_int_state), POINTER      :: p_int => NULL()
    TYPE(t_patch),      POINTER     :: p_pp => NULL()
    TYPE(t_patch),      POINTER     :: p_pc => NULL()

    ! Indices
    INTEGER :: jb, jc, jk, js, i_chidx, i_startblk, i_endblk, &
      i_startidx, i_endidx, istartblk_c
    INTEGER :: nlev_c, nlev_p
    INTEGER :: nshift      !< difference between upper boundary of parent or feedback-parent
    !< domain and upper boundary of child domain (in terms
    !< of vertical levels)

    ! Local arrays for interpolated parent-grid values, and difference fields. These have
    ! to be allocatable because their dimensions differ between MPI and non-MPI runs
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)   :: parent_rho, diff_rho

    INTEGER, DIMENSION(:,:,:), POINTER :: iidx, iblk
    LOGICAL :: l_parallel
    REAL(wp), DIMENSION(:,:,:), POINTER :: p_fbkwgt
    !-----------------------------------------------------------------------


    l_parallel = my_process_is_mpi_parallel()

    p_parent_prog     => p_nh_state(jgp)%prog(nsav1(jgp))
    p_child_prog      => p_nh_state(jg)%prog(nnow(jg))
    p_diag            => p_nh_state(jg)%diag
    p_pc              => p_patch(jg)

    p_grf => p_grf_state_local_parent(jg)
    p_int => p_int_state_local_parent(jg)
    p_gcp => p_patch_local_parent(jg)%cells
    p_pp  => p_patch_local_parent(jg)

    i_chidx  = p_pc%parent_child_index

    ! number of full levels of child domain
    nlev_c   = p_pc%nlev

    ! number of full/half levels of parent domain
    nlev_p   = p_pp%nlev

    ! shift between upper model boundaries
    nshift = p_pc%nshift
    js     = nshift

    ! Please note: In the parallel case
    ! - lower bound must be 1 due to synchronization calls
    ! - upper bound must be nblks_c/e to include halo cells/edges
    ! - this doesn't cost extra memory since p_patch_local_parent
    !   only includes the cells/edges really needed

    ! Value of i_startblk needed for subroutine call for parent-to-child interpolation
    istartblk_c = 1

    ALLOCATE(parent_rho  (nproma, nlev_p, p_patch_local_parent(jg)%nblks_c), &
      &      diff_rho    (nproma, nlev_c, p_patch_local_parent(jg)%nblks_c)  )


    ! Set pointers to index and coefficient fields for cell-based variables
    iidx  => p_gcp%child_idx
    iblk  => p_gcp%child_blk

    IF (grf_scalfbk == 1) THEN
      p_fbkwgt    => p_grf%fbk_wgt_aw
    ELSE
      p_fbkwgt    => p_grf%fbk_wgt_bln
    ENDIF

    ! 1st step: Copy prognostic variables from parent grid to fields on feedback-parent grid
    ! (trivial without MPI parallelization, but communication call needed for MPI)

    CALL exchange_data(p_pp%comm_pat_glb_to_loc_c,        &
      &                RECV=parent_rho, SEND=p_parent_prog%rho )

    ! 2nd step: perform feedback from refined grid to intermediate grid and compute differences

!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)

    ! Start/End block in the parent domain
    i_startblk = p_gcp%start_blk(grf_nudgintp_start_c+1,i_chidx)
    i_endblk   = p_gcp%end_blk(min_rlcell_int,i_chidx)


!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_pp, jb, i_startblk, i_endblk, i_startidx, i_endidx, &
        grf_nudgintp_start_c+1, min_rlcell_int)

#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jk = 1, nlev_c
#else
!CDIR UNROLL=8
      DO jk = 1, nlev_c
        DO jc = i_startidx, i_endidx
#endif

          diff_rho(jc,jk,jb) = parent_rho(jc,jk+js,jb) - (                       &
            p_child_prog%rho(iidx(jc,jb,1),jk,iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
            p_child_prog%rho(iidx(jc,jb,2),jk,iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
            p_child_prog%rho(iidx(jc,jb,3),jk,iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
            p_child_prog%rho(iidx(jc,jb,4),jk,iblk(jc,jb,4))*p_fbkwgt(jc,jb,4))+ &
            p_nh_state(jg)%metrics%rho_ref_corr(jc,jk,jb)

        ENDDO
      ENDDO

    ENDDO
!$OMP END DO NOWAIT

!$OMP END PARALLEL


    ! Interpolate differences to child grid; the differences are stored in the grf_tend fields

    ! interpol_scal_nudging needs all fields with full boundaries, so the arrays
    ! have to be sync'd before calling these routines.
    ! Please note that we cannot use sync_patch_array here (comparing parallel/non parallel results)
    ! since the arrays don't start with lower bound 1 in the non parallel case!

    ! Synchronization is needed after the interpolation step for cell-based variables because for
    ! those, the nudging tendencies are applied outside the dynamical core for reasons of mass consistency

    IF(l_parallel) CALL exchange_data(p_pp%comm_pat_c, diff_rho)
    CALL interpol_scal_nudging (p_pp, p_int, p_grf%p_dom(i_chidx), i_chidx, 0, 1, istartblk_c, &
      &                         f3din1=diff_rho, f3dout1=p_diag%grf_tend_rho                   )
    CALL sync_patch_array(SYNC_C,p_pc,p_diag%grf_tend_rho)


    DEALLOCATE(parent_rho, diff_rho)


  END SUBROUTINE prep_rho_bdy_nudging


  !>
  !! This routine executes boundary nudging for the limited-area mode.
  !!
  !! @par Revision History
  !! Developed  by Guenther Zaengl, DWD, 2013-21-10
  SUBROUTINE limarea_bdy_nudging (p_patch, p_prog, p_prog_rcf, p_metrics, p_diag, &
                                  p_int, tsrat, p_latbc_const, p_latbc_old, p_latbc_new)

    TYPE(t_patch),   INTENT(IN)    :: p_patch
    TYPE(t_nh_prog), INTENT(IN)    :: p_prog, p_prog_rcf
    TYPE(t_nh_metrics), INTENT(IN) :: p_metrics
    TYPE(t_nh_diag), INTENT(INOUT) :: p_diag
    TYPE(t_int_state), INTENT(IN)  :: p_int

    REAL(wp), INTENT(IN) :: tsrat ! Ratio between advective and dynamical time step

    ! alternative input data, either for constant or time-dependent lateral boundary conditions
    TYPE(t_nh_prog), INTENT(IN), OPTIONAL :: p_latbc_const
    TYPE(t_pi_atm),  INTENT(IN), OPTIONAL :: p_latbc_old, p_latbc_new

    ! local variables
    INTEGER :: jb, jc, jk, je, ic, nlev
    INTEGER :: nshift, i_startblk, i_endblk, i_startidx, i_endidx
    REAL(wp) :: wfac_old, wfac_new, pres, temp, qv, tempv_inc, pres_inc
    REAL(wp) :: rho_tend, thv_tend, vn_tend, qv_tend
    REAL(wp) :: rd_o_cvd, rd_o_cpd, rd_o_p0ref, nudgecoeff
    REAL(wp) :: max_nudge_coeff_vn, max_nudge_coeff_thermdyn
    LOGICAL  :: bdymask(nproma)
    LOGICAL  :: lupper_bdy_nudging, lnudge_hydro_pres_ubn


    ! number of full levels of child domain
    nlev   = p_patch%nlev

    ! R/c_v (not present in physical constants)
    rd_o_cvd = 1._wp / cvd_o_rd

    ! R/c_p
    rd_o_cpd = rd / cpd

    ! R / p0ref
    rd_o_p0ref = rd / p0ref

    IF (nudging_config%ltype(indg_type%ubn)) THEN
      ! Upper boundary nudging is switched on
      lupper_bdy_nudging = .TRUE.
      nshift = nudging_config%ilev_end
      !
      ! Check if hydrostatic or nonhydrostatic thermodynamic variables shall be used for computing nudging increments 
      lnudge_hydro_pres_ubn = nudging_config%thermdyn_type == ithermdyn_type%hydrostatic .AND. ltransport
      !
      ! Max. nudging coefficients (qv is not nudged in upper boundary zone)
      max_nudge_coeff_vn       = nudging_config%max_nudge_coeff_vn
      max_nudge_coeff_thermdyn = nudging_config%max_nudge_coeff_thermdyn
    ELSE
      nshift                   = 0
      lupper_bdy_nudging       = .FALSE.
      lnudge_hydro_pres_ubn    = .FALSE.
      max_nudge_coeff_vn       = 0._wp
      max_nudge_coeff_thermdyn = 0._wp
    ENDIF

    IF (PRESENT(p_latbc_const) .AND. (PRESENT(p_latbc_old) .OR. PRESENT(p_latbc_new))) THEN

      CALL finish('prep_limarea_bdy_nudging','conflicting arguments')

    ELSE IF (PRESENT(p_latbc_const)) THEN ! Mode for constant lateral boundary data

      ! compute differences between lateral boundary data and prognostic variables

!$OMP PARALLEL
!$OMP DO PRIVATE(jk,jc,jb,ic,rho_tend,thv_tend) ICON_OMP_DEFAULT_SCHEDULE
#ifdef __LOOP_EXCHANGE
      DO ic = 1, p_metrics%nudge_c_dim
        jc = p_metrics%nudge_c_idx(ic)
        jb = p_metrics%nudge_c_blk(ic)
!DIR$ IVDEP
        DO jk = 1, nlev
#else
      DO jk = 1, nlev
!CDIR NODEP,VOVERTAKE,VOB
        DO ic = 1, p_metrics%nudge_c_dim
          jc = p_metrics%nudge_c_idx(ic)
          jb = p_metrics%nudge_c_blk(ic)
#endif
          thv_tend = p_latbc_const%theta_v(jc,jk,jb) - p_prog%theta_v(jc,jk,jb)
          rho_tend = p_latbc_const%rho    (jc,jk,jb) - p_prog%rho    (jc,jk,jb)

          p_prog%rho(jc,jk,jb)     = p_prog%rho(jc,jk,jb)     + tsrat*p_int%nudgecoeff_c(jc,jb)*rho_tend
          p_prog%theta_v(jc,jk,jb) = p_prog%theta_v(jc,jk,jb) + tsrat*p_int%nudgecoeff_c(jc,jb)*thv_tend
          p_prog%exner(jc,jk,jb)   = EXP(rd_o_cvd*LOG(rd_o_p0ref*p_prog%rho(jc,jk,jb)*p_prog%theta_v(jc,jk,jb)))

        ENDDO
      ENDDO
!$OMP END DO

!$OMP DO PRIVATE(jk,je,jb,ic,vn_tend) ICON_OMP_DEFAULT_SCHEDULE
#ifdef __LOOP_EXCHANGE
        DO ic = 1, p_metrics%nudge_e_dim
          je = p_metrics%nudge_e_idx(ic)
          jb = p_metrics%nudge_e_blk(ic)
!DIR$ IVDEP
          DO jk = 1, nlev
#else
        DO jk = 1, nlev
!CDIR NODEP,VOVERTAKE,VOB
          DO ic = 1, p_metrics%nudge_e_dim
            je = p_metrics%nudge_e_idx(ic)
            jb = p_metrics%nudge_e_blk(ic)
#endif
            vn_tend = p_latbc_const%vn(je,jk,jb) - p_prog%vn(je,jk,jb)

            p_prog%vn(je,jk,jb) = p_prog%vn(je,jk,jb) + tsrat*p_int%nudgecoeff_e(je,jb)*vn_tend

        ENDDO
      ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    ELSE IF (PRESENT(p_latbc_old) .AND. PRESENT(p_latbc_new)) THEN ! Mode for time-dependent lateral boundary data

      ! compute differences between lateral boundary data and prognostic variables

      wfac_old = latbc_config%lc1
      wfac_new = latbc_config%lc2

!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)

      IF (latbc_config%nudge_hydro_pres .AND. ltransport) THEN
!$OMP DO PRIVATE(jk,jc,jb,ic,pres,temp,qv,tempv_inc,pres_inc,rho_tend,thv_tend) ICON_OMP_DEFAULT_SCHEDULE
#ifdef __LOOP_EXCHANGE
        DO ic = 1, p_metrics%nudge_c_dim
          jc = p_metrics%nudge_c_idx(ic)
          jb = p_metrics%nudge_c_blk(ic)
!DIR$ IVDEP
          DO jk = nshift+1, nlev
#else
        DO jk = nshift+1, nlev
!CDIR NODEP,VOVERTAKE,VOB
          DO ic = 1, p_metrics%nudge_c_dim
            jc = p_metrics%nudge_c_idx(ic)
            jb = p_metrics%nudge_c_blk(ic)
#endif
            pres = wfac_old*p_latbc_old%pres(jc,jk,jb) + wfac_new*p_latbc_new%pres(jc,jk,jb)
            temp = wfac_old*p_latbc_old%temp(jc,jk,jb) + wfac_new*p_latbc_new%temp(jc,jk,jb)
            qv   = wfac_old*p_latbc_old%qv(jc,jk,jb)   + wfac_new*p_latbc_new%qv(jc,jk,jb)

            tempv_inc = (temp-p_diag%temp(jc,jk,jb))*(1._wp+vtmpc1*qv) + &
               (qv-p_prog_rcf%tracer(jc,jk,jb,iqv))*vtmpc1*temp
            pres_inc  = pres-p_diag%pres(jc,jk,jb)

            thv_tend = tempv_inc/p_prog%exner(jc,jk,jb) - rd_o_cpd*p_prog%theta_v(jc,jk,jb)/pres*pres_inc
            rho_tend = ( pres_inc/p_diag%tempv(jc,jk,jb) - &
              tempv_inc*p_diag%pres(jc,jk,jb)/p_diag%tempv(jc,jk,jb)**2 )/rd

            p_prog%rho(jc,jk,jb)     = p_prog%rho(jc,jk,jb)     + tsrat*p_int%nudgecoeff_c(jc,jb)*rho_tend
            p_prog%theta_v(jc,jk,jb) = p_prog%theta_v(jc,jk,jb) + tsrat*p_int%nudgecoeff_c(jc,jb)*thv_tend
            p_prog%exner(jc,jk,jb)   = EXP(rd_o_cvd*LOG(rd_o_p0ref*p_prog%rho(jc,jk,jb)*p_prog%theta_v(jc,jk,jb)))

          ENDDO
        ENDDO
!$OMP END DO

      ELSE

!$OMP DO PRIVATE(jk,jc,jb,ic,rho_tend,thv_tend) ICON_OMP_DEFAULT_SCHEDULE
#ifdef __LOOP_EXCHANGE
        DO ic = 1, p_metrics%nudge_c_dim
          jc = p_metrics%nudge_c_idx(ic)
          jb = p_metrics%nudge_c_blk(ic)
!DIR$ IVDEP
          DO jk = nshift+1, nlev
#else
        DO jk = nshift+1, nlev
!CDIR NODEP,VOVERTAKE,VOB
          DO ic = 1, p_metrics%nudge_c_dim
            jc = p_metrics%nudge_c_idx(ic)
            jb = p_metrics%nudge_c_blk(ic)
#endif
            thv_tend = wfac_old*p_latbc_old%theta_v(jc,jk,jb) + wfac_new*p_latbc_new%theta_v(jc,jk,jb) - p_prog%theta_v(jc,jk,jb)
            rho_tend = wfac_old*p_latbc_old%rho(jc,jk,jb) + wfac_new*p_latbc_new%rho(jc,jk,jb)- p_prog%rho(jc,jk,jb)

            p_prog%rho(jc,jk,jb)     = p_prog%rho(jc,jk,jb)     + tsrat*p_int%nudgecoeff_c(jc,jb)*rho_tend
            p_prog%theta_v(jc,jk,jb) = p_prog%theta_v(jc,jk,jb) + tsrat*p_int%nudgecoeff_c(jc,jb)*thv_tend
            p_prog%exner(jc,jk,jb)   = EXP(rd_o_cvd*LOG(rd_o_p0ref*p_prog%rho(jc,jk,jb)*p_prog%theta_v(jc,jk,jb)))

          ENDDO
        ENDDO
!$OMP END DO

      ENDIF

      IF (ltransport) THEN ! apply QV nudging in subsaturated regions
!$OMP DO PRIVATE(jk,jc,jb,ic,qv_tend) ICON_OMP_DEFAULT_SCHEDULE
#ifdef __LOOP_EXCHANGE
        DO ic = 1, p_metrics%nudge_c_dim
          jc = p_metrics%nudge_c_idx(ic)
          jb = p_metrics%nudge_c_blk(ic)
!DIR$ IVDEP
          DO jk = 1, nlev
#else
        DO jk = 1, nlev
!CDIR NODEP,VOVERTAKE,VOB
          DO ic = 1, p_metrics%nudge_c_dim
            jc = p_metrics%nudge_c_idx(ic)
            jb = p_metrics%nudge_c_blk(ic)
#endif
            qv_tend = wfac_old*p_latbc_old%qv(jc,jk,jb) + wfac_new*p_latbc_new%qv(jc,jk,jb) - p_prog_rcf%tracer(jc,jk,jb,iqv)

            ! Suppress positive nudging tendencies in saturated (=cloudy) regions in order to avoid runaway effects
            qv_tend = MERGE(MIN(0._wp,qv_tend), qv_tend, p_prog_rcf%tracer(jc,jk,jb,iqc) > 1.e-10_wp)

            ! using a weaker nudging coefficient for QV than for thermodynamic variables turned out to have a slightly
            ! beneficial impact on forecast quality
            p_prog_rcf%tracer(jc,jk,jb,iqv) = p_prog_rcf%tracer(jc,jk,jb,iqv) + 0.5_wp*tsrat*p_int%nudgecoeff_c(jc,jb)*qv_tend

          ENDDO
        ENDDO
!$OMP END DO
      ENDIF


!$OMP DO PRIVATE(jk,je,jb,ic,vn_tend) ICON_OMP_DEFAULT_SCHEDULE
#ifdef __LOOP_EXCHANGE
        DO ic = 1, p_metrics%nudge_e_dim
          je = p_metrics%nudge_e_idx(ic)
          jb = p_metrics%nudge_e_blk(ic)
!DIR$ IVDEP
          DO jk = nshift+1, nlev
#else
        DO jk = nshift+1, nlev
!CDIR NODEP,VOVERTAKE,VOB
          DO ic = 1, p_metrics%nudge_e_dim
            je = p_metrics%nudge_e_idx(ic)
            jb = p_metrics%nudge_e_blk(ic)
#endif
            vn_tend = wfac_old*p_latbc_old%vn(je,jk,jb) + wfac_new*p_latbc_new%vn(je,jk,jb) - p_prog%vn(je,jk,jb)

            ! using a weaker nudging coefficient for vn than for thermodynamic variables turned out to have a
            ! beneficial impact on forecast quality
            p_prog%vn(je,jk,jb) = p_prog%vn(je,jk,jb) + 0.5_wp*tsrat*p_int%nudgecoeff_e(je,jb)*vn_tend
        ENDDO
      ENDDO
!$OMP END DO


    IF (lupper_bdy_nudging) THEN

    i_startblk = p_patch%cells%start_block(grf_bdywidth_c+1)
    i_endblk   = p_patch%cells%end_block(min_rlcell)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk,pres,temp,qv,tempv_inc,pres_inc,rho_tend,thv_tend,&
!$OMP            nudgecoeff,bdymask) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
        i_startidx, i_endidx, grf_bdywidth_c+1, min_rlcell)

      ! Exclude halo points of boundary interpolation zone (causes sync error otherwise)
      DO jc = i_startidx, i_endidx
        bdymask(jc) = p_patch%cells%refin_ctrl(jc,jb)>=1 .AND. p_patch%cells%refin_ctrl(jc,jb)<=grf_bdywidth_c
      ENDDO

      IF (lnudge_hydro_pres_ubn) THEN
        DO jk = 1, nshift
!DIR$ IVDEP
          DO jc = i_startidx, i_endidx

            nudgecoeff = MERGE(0._wp, tsrat*MAX(p_int%nudgecoeff_c(jc,jb), &
              &          max_nudge_coeff_thermdyn*p_metrics%nudgecoeff_vert(jk)), bdymask(jc))

            pres = wfac_old*p_latbc_old%pres(jc,jk,jb) + wfac_new*p_latbc_new%pres(jc,jk,jb)
            temp = wfac_old*p_latbc_old%temp(jc,jk,jb) + wfac_new*p_latbc_new%temp(jc,jk,jb)
            qv   = wfac_old*p_latbc_old%qv(jc,jk,jb)   + wfac_new*p_latbc_new%qv(jc,jk,jb)

            tempv_inc = (temp-p_diag%temp(jc,jk,jb))*(1._wp+vtmpc1*qv) + &
               (qv-p_prog_rcf%tracer(jc,jk,jb,iqv))*vtmpc1*temp
            pres_inc  = pres-p_diag%pres(jc,jk,jb)

            thv_tend = tempv_inc/p_prog%exner(jc,jk,jb) - rd_o_cpd*p_prog%theta_v(jc,jk,jb)/pres*pres_inc
            rho_tend = ( pres_inc/p_diag%tempv(jc,jk,jb) -                 &
              tempv_inc*p_diag%pres(jc,jk,jb)/p_diag%tempv(jc,jk,jb)**2 )/rd

            p_prog%rho(jc,jk,jb)     = p_prog%rho(jc,jk,jb)     + nudgecoeff*rho_tend
            p_prog%theta_v(jc,jk,jb) = p_prog%theta_v(jc,jk,jb) + nudgecoeff*thv_tend
            p_prog%exner(jc,jk,jb)   = MERGE(p_prog%exner(jc,jk,jb), &
              EXP(rd_o_cvd*LOG(rd_o_p0ref*p_prog%rho(jc,jk,jb)*p_prog%theta_v(jc,jk,jb))),bdymask(jc))

        ENDDO
      ENDDO
      ELSE
        DO jk = 1, nshift
!DIR$ IVDEP
          DO jc = i_startidx, i_endidx

            nudgecoeff = MERGE(0._wp, tsrat*MAX(p_int%nudgecoeff_c(jc,jb), &
              &          max_nudge_coeff_thermdyn*p_metrics%nudgecoeff_vert(jk)), bdymask(jc))

            thv_tend = wfac_old*p_latbc_old%theta_v(jc,jk,jb) + wfac_new*p_latbc_new%theta_v(jc,jk,jb) - p_prog%theta_v(jc,jk,jb)
            rho_tend = wfac_old*p_latbc_old%rho(jc,jk,jb) + wfac_new*p_latbc_new%rho(jc,jk,jb)- p_prog%rho(jc,jk,jb)

            p_prog%rho(jc,jk,jb)     = p_prog%rho(jc,jk,jb)     + nudgecoeff*rho_tend
            p_prog%theta_v(jc,jk,jb) = p_prog%theta_v(jc,jk,jb) + nudgecoeff*thv_tend
            p_prog%exner(jc,jk,jb)   = MERGE(p_prog%exner(jc,jk,jb), &
              EXP(rd_o_cvd*LOG(rd_o_p0ref*p_prog%rho(jc,jk,jb)*p_prog%theta_v(jc,jk,jb))),bdymask(jc))

        ENDDO
      ENDDO
      ENDIF
    ENDDO
!$OMP END DO

    i_startblk = p_patch%edges%start_block(grf_bdywidth_e+1)
    i_endblk   = p_patch%edges%end_block(min_rledge)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,je,jk,vn_tend,nudgecoeff,bdymask) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
        i_startidx, i_endidx, grf_bdywidth_e+1, min_rledge)

      ! Exclude halo points of boundary interpolation zone (causes sync error otherwise)
      DO je = i_startidx, i_endidx
        bdymask(je) = p_patch%edges%refin_ctrl(je,jb)>=1 .AND. p_patch%edges%refin_ctrl(je,jb)<=grf_bdywidth_e
      ENDDO

      DO jk = 1, nshift
!DIR$ IVDEP
        DO je = i_startidx, i_endidx

            nudgecoeff = MERGE(0._wp, tsrat*MAX(0.5_wp*p_int%nudgecoeff_e(je,jb), &
              &          max_nudge_coeff_vn*p_metrics%nudgecoeff_vert(jk)), bdymask(je))

            vn_tend = wfac_old*p_latbc_old%vn(je,jk,jb) + wfac_new*p_latbc_new%vn(je,jk,jb) - p_prog%vn(je,jk,jb)

            p_prog%vn(je,jk,jb) = p_prog%vn(je,jk,jb) + nudgecoeff*vn_tend
        ENDDO
      ENDDO
    ENDDO
!$OMP END DO NOWAIT

    ENDIF

!$OMP END PARALLEL

    ELSE
      CALL finish('limarea_bdy_nudging','missing arguments')
    ENDIF

  END SUBROUTINE limarea_bdy_nudging



  !>
  !! This routine executes boundary nudging for one-way nested domains
  !!
  !!
  !! @par Revision History
  !! Developed  by Guenther Zaengl, DWD, 2010-06-18
  SUBROUTINE nest_boundary_nudging(jg, nnew, nnew_rcf, rcffac)


    INTEGER, INTENT(IN)  :: jg, nnew, nnew_rcf

    REAL(wp), INTENT(IN) :: rcffac ! Ratio between advective and dynamical time step

    ! Pointers
    TYPE(t_nh_state),  POINTER ::  p_nh
    TYPE(t_int_state), POINTER ::  p_int

    ! Indices
    INTEGER :: jb, jc, je, jk, jt, ic

    INTEGER :: nlev          ! number of vertical full levels
    INTEGER :: ntracer_nudge !< number of tracers to be nudged

    REAL(wp) :: rd_o_cvd, rd_o_p0ref, upper_lim, lower_lim

    ! Set pointers
    p_nh  => p_nh_state(jg)
    p_int => p_int_state(jg)

    IF (iforcing > 1) THEN  ! tracers represent moisture variables
      ntracer_nudge = 3     ! take only QV, QC and QI - nudging precip variables does not make much sense
    ELSE
      ntracer_nudge = ntracer
    ENDIF

    ! R/c_v (not present in physical constants)
    rd_o_cvd = 1._wp / cvd_o_rd

    ! R / p0ref
    rd_o_p0ref = rd / p0ref

    ! number of vertical levels
    nlev = p_patch(jg)%nlev

!$OMP PARALLEL

!$OMP DO PRIVATE(jk,jc,jb,ic,upper_lim,lower_lim) ICON_OMP_DEFAULT_SCHEDULE
#ifdef __LOOP_EXCHANGE
    DO ic = 1, p_nh%metrics%nudge_c_dim
      jc = p_nh%metrics%nudge_c_idx(ic)
      jb = p_nh%metrics%nudge_c_blk(ic)
!DIR$ IVDEP
      DO jk = 1, nlev
#else
    DO jk = 1, nlev
!CDIR NODEP,VOVERTAKE,VOB
      DO ic = 1, p_nh%metrics%nudge_c_dim
        jc = p_nh%metrics%nudge_c_idx(ic)
        jb = p_nh%metrics%nudge_c_blk(ic)
#endif
        upper_lim = 1.0025_wp*p_nh%prog(nnew)%rho(jc,jk,jb)
        lower_lim = 0.9975_wp*p_nh%prog(nnew)%rho(jc,jk,jb)
        p_nh%prog(nnew)%rho(jc,jk,jb) =                                      &
          p_nh%prog(nnew)%rho(jc,jk,jb) + rcffac*p_int%nudgecoeff_c(jc,jb)*  &
          p_nh%diag%grf_tend_rho(jc,jk,jb)
        p_nh%prog(nnew)%rho(jc,jk,jb) = MAX(lower_lim,MIN(upper_lim,p_nh%prog(nnew)%rho(jc,jk,jb)))

        upper_lim = 1.0025_wp*p_nh%prog(nnew)%theta_v(jc,jk,jb)
        lower_lim = 0.9975_wp*p_nh%prog(nnew)%theta_v(jc,jk,jb)
        p_nh%prog(nnew)%theta_v(jc,jk,jb) =                                     &
          p_nh%prog(nnew)%theta_v(jc,jk,jb) + rcffac*p_int%nudgecoeff_c(jc,jb)* &
          p_nh%diag%grf_tend_thv(jc,jk,jb)
        p_nh%prog(nnew)%theta_v(jc,jk,jb) = MAX(lower_lim,MIN(upper_lim,p_nh%prog(nnew)%theta_v(jc,jk,jb)))

        p_nh%prog(nnew)%exner(jc,jk,jb) =                                  &
          EXP(rd_o_cvd*LOG(rd_o_p0ref*p_nh%prog(nnew)%rho(jc,jk,jb)*p_nh%prog(nnew)%theta_v(jc,jk,jb)))

        p_nh%prog(nnew)%w(jc,jk,jb) =                                      &
          p_nh%prog(nnew)%w(jc,jk,jb) + rcffac*p_int%nudgecoeff_c(jc,jb)*  &
          p_nh%diag%grf_tend_w(jc,jk,jb)

      ENDDO
    ENDDO
!$OMP END DO

!$OMP DO PRIVATE(jb,jk,je,ic) ICON_OMP_DEFAULT_SCHEDULE
#ifdef __LOOP_EXCHANGE
    DO ic = 1, p_nh%metrics%nudge_e_dim
      je = p_nh%metrics%nudge_e_idx(ic)
      jb = p_nh%metrics%nudge_e_blk(ic)
!DIR$ IVDEP
      DO jk = 1, nlev
#else
    DO jk = 1, nlev
!CDIR NODEP,VOVERTAKE,VOB
      DO ic = 1, p_nh%metrics%nudge_e_dim
        je = p_nh%metrics%nudge_e_idx(ic)
        jb = p_nh%metrics%nudge_e_blk(ic)
#endif
        p_nh%prog(nnew)%vn(je,jk,jb) = p_nh%prog(nnew)%vn(je,jk,jb)        &
          + rcffac*p_int%nudgecoeff_e(je,jb)*p_nh%diag%grf_tend_vn(je,jk,jb)
      ENDDO
    ENDDO
!$OMP END DO

    IF (ltransport) THEN
!$OMP DO PRIVATE(jk,jc,jb,jt,ic,upper_lim,lower_lim) ICON_OMP_DEFAULT_SCHEDULE
#ifdef __LOOP_EXCHANGE
      DO ic = 1, p_nh%metrics%nudge_c_dim
        jc = p_nh%metrics%nudge_c_idx(ic)
        jb = p_nh%metrics%nudge_c_blk(ic)
        DO jt = 1, ntracer_nudge
!DIR$ IVDEP
          DO jk = 1, nlev
#else
      DO jt = 1, ntracer_nudge
        DO jk = 1, nlev
!CDIR NODEP,VOVERTAKE,VOB
          DO ic = 1, p_nh%metrics%nudge_c_dim
            jc = p_nh%metrics%nudge_c_idx(ic)
            jb = p_nh%metrics%nudge_c_blk(ic)
#endif

            upper_lim = 1.01_wp*p_nh%prog(nnew_rcf)%tracer(jc,jk,jb,jt)
            lower_lim = 0.99_wp*p_nh%prog(nnew_rcf)%tracer(jc,jk,jb,jt)
            p_nh%prog(nnew_rcf)%tracer(jc,jk,jb,jt) = p_nh%prog(nnew_rcf)%tracer(jc,jk,jb,jt) + &
              rcffac*p_int%nudgecoeff_c(jc,jb)*p_nh%diag%grf_tend_tracer(jc,jk,jb,jt)
            p_nh%prog(nnew_rcf)%tracer(jc,jk,jb,jt) =  MAX(lower_lim,MIN(upper_lim,&
              p_nh%prog(nnew_rcf)%tracer(jc,jk,jb,jt)))

          ENDDO
        ENDDO
      ENDDO
!$OMP END DO NOWAIT
    ENDIF

!$OMP END PARALLEL


  END SUBROUTINE nest_boundary_nudging

  !>
  !! This routine executes boundary nudging for density (for use with 2-way nesting)
  !!
  !!
  !! @par Revision History
  !! Developed  by Guenther Zaengl, DWD, 2011-12-08
  SUBROUTINE density_boundary_nudging(jg, nnew, rcffac)


    INTEGER, INTENT(IN)  :: jg, nnew

    REAL(wp), INTENT(IN) :: rcffac ! Ratio between advective and dynamical time step

    ! Pointers
    TYPE(t_nh_state),  POINTER ::  p_nh
    TYPE(t_int_state), POINTER ::  p_int

    ! Indices
    INTEGER :: jb, jc, jk, ic

    INTEGER :: nlev  ! number of vertical full levels

    REAL(wp) :: rd_o_cvd, rd_o_p0ref

    ! Set pointers
    p_nh  => p_nh_state(jg)
    p_int => p_int_state(jg)

    ! R/c_v (not present in physical constants)
    rd_o_cvd = 1._wp / cvd_o_rd

    ! R / p0ref
    rd_o_p0ref = rd / p0ref

    ! number of vertical levels
    nlev = p_patch(jg)%nlev

!$OMP PARALLEL

!$OMP DO PRIVATE(jk,jc,jb,ic) ICON_OMP_DEFAULT_SCHEDULE
#ifdef __LOOP_EXCHANGE
    DO ic = 1, p_nh%metrics%nudge_c_dim
      jc = p_nh%metrics%nudge_c_idx(ic)
      jb = p_nh%metrics%nudge_c_blk(ic)
!DIR$ IVDEP
      DO jk = 1, nlev
#else
    DO jk = 1, nlev
!CDIR NODEP,VOVERTAKE,VOB
      DO ic = 1, p_nh%metrics%nudge_c_dim
        jc = p_nh%metrics%nudge_c_idx(ic)
        jb = p_nh%metrics%nudge_c_blk(ic)
#endif
        p_nh%prog(nnew)%rho(jc,jk,jb) = p_nh%prog(nnew)%rho(jc,jk,jb) +  &
          MIN(0.333_wp,3._wp*rcffac*p_int%nudgecoeff_c(jc,jb))*          &
          p_nh%diag%grf_tend_rho(jc,jk,jb)

        p_nh%prog(nnew)%exner(jc,jk,jb) =                                  &
          EXP(rd_o_cvd*LOG(rd_o_p0ref*p_nh%prog(nnew)%rho(jc,jk,jb)*p_nh%prog(nnew)%theta_v(jc,jk,jb)))

      ENDDO
    ENDDO
!$OMP END DO NOWAIT

!$OMP END PARALLEL

  END SUBROUTINE density_boundary_nudging


END MODULE mo_nh_nest_utilities

