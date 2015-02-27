!>
!! This module contains the driver routines needed for initializing nested
!! domains that are started sometime during the integration
!!
!! @author Guenther Zaengl, DWD
!!
!!
!! @par Revision History
!! First version by Guenther Zaengl, DWD (2012-06-06)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_nh_init_nest_utils

  USE mo_kind,                  ONLY: wp
  USE mo_model_domain,          ONLY: t_patch, p_patch, p_patch_local_parent
  USE mo_nonhydro_types,        ONLY: t_nh_metrics, t_nh_prog, t_nh_diag
  USE mo_nonhydro_state,        ONLY: p_nh_state
  USE mo_initicon_types,        ONLY: t_initicon_state
  USE mo_nwp_phy_state,         ONLY: prm_diag
  USE mo_parallel_config,       ONLY: nproma, p_test_run
  USE mo_run_config,            ONLY: ltransport, msg_level, ntracer, iforcing
  USE mo_dynamics_config,       ONLY: nnow, nnow_rcf, nnew_rcf
  USE mo_physical_constants,    ONLY: rd, cvd_o_rd, p0ref, rhoh2o
  USE mo_phyparam_soil,         ONLY: crhosminf
  USE mo_impl_constants,        ONLY: min_rlcell, min_rlcell_int, min_rledge_int, &
    &                                 min_rlvert, min_rlvert_int, MAX_CHAR_LENGTH,&
    &                                 dzsoil, inwp
  USE mo_grf_nudgintp,          ONLY: interpol_scal_nudging, interpol_vec_nudging
  USE mo_grf_bdyintp,           ONLY: interpol_scal_grf, interpol2_vec_grf
  USE mo_grid_config,           ONLY: lfeedback, ifeedback_type
  USE mo_exception,             ONLY: message, message_text, finish
  USE mo_mpi,                   ONLY: my_process_is_mpi_parallel
  USE mo_communication,         ONLY: exchange_data, exchange_data_mult
  USE mo_sync,                  ONLY: sync_patch_array, sync_patch_array_mult, &
                                      SYNC_C, SYNC_E, SYNC_V
  USE mo_intp_data_strc,        ONLY: t_int_state, p_int_state, p_int_state_local_parent
  USE mo_grf_intp_data_strc,    ONLY: t_gridref_single_state, t_gridref_state, &
                                      p_grf_state, p_grf_state_local_parent
  USE mo_loopindices,           ONLY: get_indices_c, get_indices_e, get_indices_v
  USE mo_impl_constants_grf,    ONLY: grf_bdywidth_c, grf_bdywidth_e, grf_fbk_start_c, &
                                      grf_nudgintp_start_c, grf_nudgintp_start_e
  USE mo_nwp_lnd_types,         ONLY: t_lnd_prog, t_lnd_diag, t_wtr_prog
  USE mo_lnd_nwp_config,        ONLY: ntiles_total, ntiles_water, nlev_soil, lseaice,  &
    &                                 lmulti_snow, nlev_snow, llake, isub_lake
  USE mo_nwp_lnd_state,         ONLY: p_lnd_state
  USE mo_atm_phy_nwp_config,    ONLY: atm_phy_nwp_config
  USE mo_interpol_config,       ONLY: nudge_zone_width
  USE mo_intp,                  ONLY: cells2verts_scalar
  USE mo_ext_data_state,        ONLY: ext_data
  USE mo_nh_diagnose_pres_temp, ONLY: diagnose_pres_temp
  USE mo_intp_rbf,              ONLY: rbf_vec_interpol_cell
  USE mo_seaice_nwp,            ONLY: frsi_min
  USE mo_nwp_sfc_interp,        ONLY: smi_to_wsoil, wsoil_to_smi
  USE mo_flake,                 ONLY: flake_coldinit

  IMPLICIT NONE

  PRIVATE


  REAL(wp), PARAMETER :: rd_o_cvd = 1._wp / cvd_o_rd

  PUBLIC :: initialize_nest, topo_blending_and_fbk, interpolate_increments, &
            interpolate_sfcana

  CONTAINS
  !-------------
  !>
  !! SUBROUTINE initialize_nest
  !!
  !! Driver routine for initializing a nested domain during runtime from the parent domain
  !!
  !!
  !! @par Revision History
  !! Initial version by Guenther Zaengl, DWD(2012-06-06)
  !!
  !! Note: the interpolation of the land variables does not yet include the
  !! land-water mask and adaptations to make optimal use of tiles
  !! The aggregate_landvars routine must be called before this routine
  !!
  !!
  SUBROUTINE initialize_nest(jg, jgc)

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'initialize_nest'


    INTEGER, INTENT(IN) :: jg   ! parent (source) domain ID
    INTEGER, INTENT(IN) :: jgc  ! child  (target) domain ID


    ! local pointers
    TYPE(t_nh_prog),    POINTER     :: p_parent_prog
    TYPE(t_nh_prog),    POINTER     :: p_child_prog
    TYPE(t_nh_diag),    POINTER     :: p_child_diag
    TYPE(t_nh_prog),    POINTER     :: p_parent_prog_rcf
    TYPE(t_nh_prog),    POINTER     :: p_child_prog_rcf
    TYPE(t_nh_metrics), POINTER     :: p_parent_metrics
    TYPE(t_nh_metrics), POINTER     :: p_child_metrics
    TYPE(t_lnd_prog),   POINTER     :: p_parent_lprog
    TYPE(t_lnd_prog),   POINTER     :: p_child_lprog
    TYPE(t_lnd_prog),   POINTER     :: p_child_lprog2
    TYPE(t_wtr_prog),   POINTER     :: p_parent_wprog
    TYPE(t_wtr_prog),   POINTER     :: p_child_wprog
    TYPE(t_lnd_diag),   POINTER     :: p_parent_ldiag
    TYPE(t_lnd_diag),   POINTER     :: p_child_ldiag
    TYPE(t_gridref_state), POINTER  :: p_grf => NULL()
    TYPE(t_gridref_state), POINTER  :: p_grfc => NULL()
    TYPE(t_int_state), POINTER      :: p_int => NULL()
    TYPE(t_patch),      POINTER     :: p_pp => NULL()
    TYPE(t_patch),      POINTER     :: p_pc => NULL()

    TYPE(t_lnd_diag),   POINTER     :: lnd_diag
    TYPE(t_lnd_prog),   POINTER     :: lnd_prog

    ! local variables

    ! Indices
    INTEGER :: jb, jc, jk, jk1, jt, i_nchdom, i_chidx, i_startblk, i_endblk, &
               i_startidx, i_endidx
    INTEGER :: rl_start, rl_end
    INTEGER :: nlev_c, nlev_p
    INTEGER :: nshift      ! difference between upper boundary of parent or feedback-parent 
                           ! domain and upper boundary of child domain (in terms 
                           ! of vertical levels) 
    INTEGER :: num_lndvars, num_wtrvars, num_phdiagvars

    ! Local arrays for variables living on the local parent grid in the MPI case. These have
    ! to be allocatable because their dimensions differ between MPI and non-MPI runs
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)   :: vn_lp, w_lp, thv_pr_lp, rho_pr_lp, phdiag_lp, &
                                                 lndvars_lp, wtrvars_lp
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:,:) :: tracer_lp

    ! Local arrays on the parent or child grid. These would not have to be allocatable,
    ! but the computational overhead does not matter for an initialization routine
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)   :: thv_pr_par, rho_pr_par, lndvars_par, lndvars_chi, &
                                                 wtrvars_par, wtrvars_chi, phdiag_par, phdiag_chi
    REAL(wp), ALLOCATABLE :: tsfc_ref_p(:,:), tsfc_ref_c(:,:) ! Reference temperature at lowest level

    LOGICAL :: l_parallel, l_limit(ntracer)

    INTEGER :: i_count, ic

    !-----------------------------------------------------------------------

    IF (msg_level >= 10) THEN
      WRITE(message_text,'(a,i2,a,i2)') 'Nest initialization, domain ',jg,' =>',jgc
      CALL message(TRIM(routine),message_text)
    ENDIF

    IF (.NOT. my_process_is_mpi_parallel()) THEN
      l_parallel = .FALSE.
    ELSE
      l_parallel = .TRUE.
    ENDIF

    p_parent_prog     => p_nh_state(jg)%prog(nnow(jg))
    p_child_prog      => p_nh_state(jgc)%prog(nnow(jgc))
    p_child_diag      => p_nh_state(jgc)%diag
    p_parent_prog_rcf => p_nh_state(jg)%prog(nnow_rcf(jg))
    p_child_prog_rcf  => p_nh_state(jgc)%prog(nnow_rcf(jgc))
    p_parent_metrics  => p_nh_state(jg)%metrics
    p_child_metrics   => p_nh_state(jgc)%metrics
    p_parent_lprog    => p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))
    p_parent_wprog    => p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))
    p_child_lprog     => p_lnd_state(jgc)%prog_lnd(nnow_rcf(jgc))
    p_child_lprog2    => p_lnd_state(jgc)%prog_lnd(nnew_rcf(jgc))
    p_child_wprog     => p_lnd_state(jgc)%prog_wtr(nnow_rcf(jgc))
    p_parent_ldiag    => p_lnd_state(jg)%diag_lnd
    p_child_ldiag     => p_lnd_state(jgc)%diag_lnd
    p_grfc            => p_grf_state(jgc)
    p_pc              => p_patch(jgc)

    lnd_diag          => p_lnd_state(jg)%diag_lnd
    lnd_prog          => p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))

    p_grf => p_grf_state_local_parent(jgc)
    p_int => p_int_state_local_parent(jgc)
    p_pp  => p_patch_local_parent(jgc)

    i_nchdom = MAX(1,p_patch(jg)%n_childdom)
    i_chidx  = p_pc%parent_child_index

    ! number of full levels of child domain
    nlev_c   = p_pc%nlev

    ! number of full levels of parent domain
    nlev_p   = p_patch(jg)%nlev

    ! shift between upper model boundaries
    nshift = p_pc%nshift

    ! number of land and water variables to be interpolated
    ! Remark (GZ): the multi-layer snow variables are initialized afterwards in terra_multlay_init. This 
    ! turned out to cause occasional conflicts with directly interpolating those variables here; thus
    ! the interpolation of the multi-layer snow fields has been completely removed from this routine
    num_lndvars = 2*nlev_soil+1+ &     ! multi-layer soil variables t_so and w_so (w_so_ice is initialized in terra_multlay_init)
                  5+4+1                ! single-layer prognostic variables + t_g, freshsnow, t_seasfc and qv_s + aux variable for lake temp
    num_wtrvars  = 5                   ! water state fields + fr_seaice
    num_phdiagvars = 21                ! number of physics diagnostic variables (copied from interpol_phys_grf)
    
    ALLOCATE(thv_pr_par  (nproma, nlev_p,      p_patch(jg)%nblks_c), &
             rho_pr_par  (nproma, nlev_p,      p_patch(jg)%nblks_c), &
             lndvars_par (nproma, num_lndvars, p_patch(jg)%nblks_c), &
             wtrvars_par (nproma, num_wtrvars, p_patch(jg)%nblks_c), &
             phdiag_par  (nproma, num_phdiagvars, p_patch(jg)%nblks_c), &
             lndvars_chi (nproma, num_lndvars, p_patch(jgc)%nblks_c),&
             wtrvars_chi (nproma, num_wtrvars, p_patch(jgc)%nblks_c),&
             phdiag_chi  (nproma, num_phdiagvars, p_patch(jgc)%nblks_c),&
             tsfc_ref_p  (nproma,              p_patch(jg)%nblks_c), &
             tsfc_ref_c  (nproma,              p_patch(jgc)%nblks_c) )


    ALLOCATE(vn_lp      (nproma, nlev_p,      p_pp%nblks_e),          &
             w_lp       (nproma, nlev_p+1,    p_pp%nblks_c),          &
             thv_pr_lp  (nproma, nlev_p,      p_pp%nblks_c),          &
             rho_pr_lp  (nproma, nlev_p,      p_pp%nblks_c),          &
             phdiag_lp  (nproma, num_phdiagvars,p_pp%nblks_c),        &
             tracer_lp  (nproma, nlev_p,      p_pp%nblks_c, ntracer), &
             lndvars_lp (nproma, num_lndvars, p_pp%nblks_c),          &
             wtrvars_lp (nproma, num_wtrvars ,p_pp%nblks_c)           )

    IF (p_test_run) THEN
      lndvars_par  = 0._wp
      wtrvars_par  = 0._wp
      phdiag_par   = 0._wp
      lndvars_chi  = 0._wp
      wtrvars_chi  = 0._wp
      phdiag_chi   = 0._wp
    ENDIF



    IF (atm_phy_nwp_config(jg)%inwp_surface == 1) THEN
      ! Step 0: get lake surface temperature
      !       
      ! exclude nest boundary and halo points
      rl_start = grf_bdywidth_c+1
      rl_end   = min_rlcell_int

      i_startblk = p_patch(jg)%cells%start_blk(rl_start,1)
      i_endblk   = p_patch(jg)%cells%end_blk(rl_end,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,i_count,ic,jk) ICON_OMP_GUIDED_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(p_patch(jg), jb, i_startblk, i_endblk, &
          & i_startidx, i_endidx, rl_start, rl_end)

        ! the last element of lndvars contains the source data for lake sfc temperature;
        ! a) initialize with t_g
        lndvars_par(:,num_lndvars,jb) = lnd_prog%t_g(:,jb)

        i_count = ext_data(jg)%atm%fp_count(jb) 
!CDIR NODEP,VOVERTAKE,VOB
        DO ic = 1, i_count
          jc = ext_data(jg)%atm%idx_lst_fp(ic,jb)

          ! b) take lake sfc temperature where available
          lndvars_par(jc,num_lndvars,jb) = lnd_prog%t_g_t(jc,jb,isub_lake)

        ENDDO

      ENDDO
!$OMP END DO
!$OMP END PARALLEL

    END IF



    ! Step 1: boundary interpolation
    !
    ! Note: the sepration between boundary interpolation and the interpolation
    ! of the fields in the prognostic region of the nest is needed because the
    ! existing interpolation routines are runtime-optimized for operations that
    ! are needed regularly. Boundary interpolation is executed from the parent grid
    ! to the child grid, whereas the so-called nudging interpolation works from
    ! the local parent grid to the child grid.


    ! Step 1a: prepare boundary interpolation

!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)

    ! cell-based variables
    i_startblk = p_patch(jg)%cells%start_blk(grf_bdywidth_c+1,1)
    i_endblk   = p_patch(jg)%cells%end_blk(min_rlcell,i_nchdom)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk,jk1) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch(jg), jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, grf_bdywidth_c+1, min_rlcell)

      ! Compute perturbation quantities for rho and theta at parent level
      DO jk = 1, nlev_p
        DO jc = i_startidx, i_endidx
          rho_pr_par(jc,jk,jb) = &
            p_parent_prog%rho(jc,jk,jb) - p_parent_metrics%rho_ref_mc(jc,jk,jb)

          thv_pr_par(jc,jk,jb) = &
            p_parent_prog%theta_v(jc,jk,jb) - p_parent_metrics%theta_ref_mc(jc,jk,jb)
        ENDDO
      ENDDO

      ! Reference temperature at lowest level (taken as proxy for the surface)
      DO jc = i_startidx, i_endidx
        tsfc_ref_p(jc,jb) = p_parent_metrics%theta_ref_mc(jc,nlev_p,jb) * &
                            p_parent_metrics%exner_ref_mc(jc,nlev_p,jb)
      ENDDO

      IF (iforcing == inwp) THEN
        ! Collect diagnostic physics fields (the only really important ones are the precip fields)
        DO jc = i_startidx, i_endidx
          phdiag_par(jc,1,jb) = prm_diag(jg)%tot_prec(jc,jb)
          phdiag_par(jc,2,jb) = prm_diag(jg)%rain_gsp(jc,jb)
          phdiag_par(jc,3,jb) = prm_diag(jg)%snow_gsp(jc,jb)
          phdiag_par(jc,4,jb) = prm_diag(jg)%rain_con(jc,jb)
          phdiag_par(jc,5,jb) = prm_diag(jg)%snow_con(jc,jb)
          phdiag_par(jc,6,jb) = prm_diag(jg)%rain_gsp_rate(jc,jb)
          phdiag_par(jc,7,jb) = prm_diag(jg)%snow_gsp_rate(jc,jb)
          phdiag_par(jc,8,jb) = prm_diag(jg)%rain_con_rate(jc,jb)
          phdiag_par(jc,9,jb) = prm_diag(jg)%snow_con_rate(jc,jb)
          phdiag_par(jc,10,jb) = prm_diag(jg)%gz0(jc,jb)
          phdiag_par(jc,11,jb) = prm_diag(jg)%tcm(jc,jb)
          phdiag_par(jc,12,jb) = prm_diag(jg)%tch(jc,jb)
          phdiag_par(jc,13,jb) = prm_diag(jg)%tfm(jc,jb)
          phdiag_par(jc,14,jb) = prm_diag(jg)%tfh(jc,jb)
          phdiag_par(jc,15,jb) = prm_diag(jg)%tfv(jc,jb)
          phdiag_par(jc,16,jb) = prm_diag(jg)%t_2m(jc,jb)
          phdiag_par(jc,17,jb) = prm_diag(jg)%qv_2m(jc,jb)
          phdiag_par(jc,18,jb) = prm_diag(jg)%td_2m(jc,jb)
          phdiag_par(jc,19,jb) = prm_diag(jg)%rh_2m(jc,jb)
          phdiag_par(jc,20,jb) = prm_diag(jg)%u_10m(jc,jb)
          phdiag_par(jc,21,jb) = prm_diag(jg)%v_10m(jc,jb)
        ENDDO
      ENDIF

      IF (atm_phy_nwp_config(jg)%inwp_surface == 1) THEN
        ! Collect soil variables
        DO jk = 1, nlev_soil
          jk1 = jk + nlev_soil
          DO jc = i_startidx, i_endidx
            lndvars_par(jc,jk,jb)  = p_parent_ldiag%w_so(jc,jk,jb)
            lndvars_par(jc,jk1,jb) = p_parent_ldiag%t_so(jc,jk,jb) - tsfc_ref_p(jc,jb)
          ENDDO
        ENDDO

        jk1 = 2*nlev_soil + 1
        DO jc = i_startidx, i_endidx
          lndvars_par(jc,jk1,jb)   = p_parent_ldiag%t_so(jc,nlev_soil+1,jb)-tsfc_ref_p(jc,jb)
          lndvars_par(jc,jk1+1,jb) = p_parent_lprog%t_g(jc,jb)    - tsfc_ref_p(jc,jb)
          IF (p_parent_ldiag%t_s(jc,jb) > 10._wp) THEN
            lndvars_par(jc,jk1+2,jb) = p_parent_ldiag%t_s(jc,jb)  - tsfc_ref_p(jc,jb)
          ELSE ! t_s has missing values over water - use t_g instead
            lndvars_par(jc,jk1+2,jb) = lndvars_par(jc,jk1+1,jb)
          ENDIF
          lndvars_par(jc,jk1+3,jb) = p_parent_ldiag%t_snow(jc,jb) - tsfc_ref_p(jc,jb)
          lndvars_par(jc,jk1+4,jb) = p_parent_ldiag%w_snow(jc,jb)
          lndvars_par(jc,jk1+5,jb) = p_parent_ldiag%rho_snow(jc,jb)
          lndvars_par(jc,jk1+6,jb) = p_parent_ldiag%w_i(jc,jb)
          lndvars_par(jc,jk1+7,jb) = p_parent_ldiag%freshsnow(jc,jb)
          lndvars_par(jc,jk1+8,jb) = MAX(271._wp,p_parent_ldiag%t_seasfc(jc,jb)) ! to avoid trouble with missing values
          lndvars_par(jc,jk1+9,jb) = p_parent_ldiag%qv_s(jc,jb)
        ENDDO
      ENDIF

      IF (atm_phy_nwp_config(jg)%inwp_surface == 1 .AND. lseaice) THEN
        DO jc = i_startidx, i_endidx
          IF (p_parent_wprog%t_ice(jc,jb) > 10._wp) THEN
            wtrvars_par(jc,1,jb) = p_parent_wprog%t_ice(jc,jb)
          ELSE
            wtrvars_par(jc,1,jb) = p_parent_lprog%t_g(jc,jb)
          ENDIF
          wtrvars_par(jc,2,jb) = p_parent_wprog%h_ice(jc,jb) 
          wtrvars_par(jc,3,jb) = p_parent_wprog%t_snow_si(jc,jb) 
          wtrvars_par(jc,4,jb) = p_parent_wprog%h_snow_si(jc,jb)
          wtrvars_par(jc,5,jb) = p_parent_ldiag%fr_seaice(jc,jb)
        ENDDO
      ENDIF

    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    ! Convert wsoil into SMI for interpolation
    CALL wsoil_to_smi(p_patch(jg), lndvars_par(:,1:nlev_soil,:))

    ! Step 1b: execute boundary interpolation

    CALL interpol2_vec_grf (p_patch(jg), p_pc, p_grf_state(jg)%p_dom(i_chidx), 1, &
      p_parent_prog%vn, p_child_prog%vn)

    CALL interpol_scal_grf (p_patch(jg), p_pc, p_grf_state(jg)%p_dom(i_chidx), 3, &
      rho_pr_par,      p_child_prog%rho,                                          &
      thv_pr_par,      p_child_prog%theta_v,                                      &
      p_parent_prog%w, p_child_prog%w                                             )

    IF (ltransport) THEN
      l_limit(:) = .TRUE. ! apply positive definite limiter on tracers

      CALL interpol_scal_grf ( p_patch(jg), p_pc, p_grf_state(jg)%p_dom(i_chidx), ntracer,   &
        f4din1=p_parent_prog_rcf%tracer, f4dout1=p_child_prog_rcf%tracer, llimit_nneg=l_limit)
    ENDIF

    IF (iforcing == inwp) THEN
      CALL sync_patch_array(SYNC_C,p_patch(jg),phdiag_par)
      CALL interpol_scal_grf (p_patch(jg), p_pc, p_grf_state(jg)%p_dom(i_chidx), 1, &
        phdiag_par, phdiag_chi, lnoshift=.TRUE.                 )
    ENDIF

    IF (atm_phy_nwp_config(jg)%inwp_surface == 1) THEN
      CALL sync_patch_array(SYNC_C,p_patch(jg),lndvars_par)
      CALL interpol_scal_grf (p_patch(jg), p_pc, p_grf_state(jg)%p_dom(i_chidx), 1,  &
        lndvars_par, lndvars_chi, lnoshift=.TRUE.                 )
    ENDIF

    IF (atm_phy_nwp_config(jg)%inwp_surface == 1 .AND. lseaice) THEN
      CALL sync_patch_array(SYNC_C,p_patch(jg),wtrvars_par)
      CALL interpol_scal_grf (p_patch(jg), p_pc, p_grf_state(jg)%p_dom(i_chidx), 1, &
        wtrvars_par, wtrvars_chi, lnoshift=.TRUE.                 )
    ENDIF

    ! Step 2: Interpolation of fields in the model interior

    ! Step 2a: Copy prognostic variables from parent grid to fields on feedback-parent grid
    ! (trivial without MPI parallelization, but communication call needed for MPI)

    CALL exchange_data_mult(p_pp%comm_pat_glb_to_loc_c, 3, 3*nlev_p+1, &
      &                     RECV1=rho_pr_lp, SEND1=rho_pr_par,         &
      &                     RECV2=thv_pr_lp, SEND2=thv_pr_par,         &
      &                     RECV3=w_lp,      SEND3=p_parent_prog%w     )

    CALL exchange_data(p_pp%comm_pat_glb_to_loc_e, RECV=vn_lp, SEND=p_parent_prog%vn)

    CALL exchange_data_mult(p_pp%comm_pat_glb_to_loc_c, ntracer, ntracer*nlev_p, &
      &                     RECV4D=tracer_lp, SEND4D=p_parent_prog_rcf%tracer    )

    IF (iforcing == inwp) &
      &      CALL exchange_data(p_pp%comm_pat_glb_to_loc_c, RECV=phdiag_lp, SEND=phdiag_par)

    IF (atm_phy_nwp_config(jg)%inwp_surface == 1) &
      &      CALL exchange_data(p_pp%comm_pat_glb_to_loc_c, RECV=lndvars_lp, SEND=lndvars_par)

    IF (atm_phy_nwp_config(jg)%inwp_surface == 1 .AND. lseaice) &
      &      CALL exchange_data(p_pp%comm_pat_glb_to_loc_c, RECV=wtrvars_lp, SEND=wtrvars_par)

    ! Step 2b: Perform interpolation from local parent to child grid

    ! Note: the sync routines cannot be used for the synchronization on the local
    ! parent grid

    IF(l_parallel) CALL exchange_data(p_pp%comm_pat_e, vn_lp)
    CALL interpol_vec_nudging (p_pp, p_pc, p_int, p_grf%p_dom(i_chidx), p_grfc, &
                               i_chidx, nshift, 1, vn_lp, p_child_prog%vn       )
    CALL sync_patch_array(SYNC_E,p_pc,p_child_prog%vn)

    IF(l_parallel) CALL exchange_data_mult(p_pp%comm_pat_c, 3, 3*nlev_p+1, &
                               recv1=thv_pr_lp, recv2=rho_pr_lp, recv3=w_lp)
    CALL interpol_scal_nudging (p_pp, p_int, p_grf%p_dom(i_chidx), i_chidx, nshift, 3, 1, &
                                f3din1=thv_pr_lp, f3dout1=p_child_prog%theta_v,           &
                                f3din2=rho_pr_lp, f3dout2=p_child_prog%rho,               &
                                f3din3=w_lp,      f3dout3=p_child_prog%w                  )
    CALL sync_patch_array_mult(SYNC_C,p_pc,3,p_child_prog%theta_v,p_child_prog%rho,  &
                               p_child_prog%w)

    IF (ltransport) THEN
      IF(l_parallel) CALL exchange_data_mult(p_pp%comm_pat_c, ntracer, ntracer*nlev_p, &
                                             recv4d=tracer_lp)
      l_limit(:) = .TRUE. ! apply positive definite limiter
      CALL interpol_scal_nudging (p_pp, p_int, p_grf%p_dom(i_chidx), i_chidx, &
                                  nshift, ntracer, 1, f4din=tracer_lp,        &
                                  f4dout=p_child_prog_rcf%tracer, llimit_nneg=l_limit)
      CALL sync_patch_array_mult(SYNC_C,p_pc,ntracer,f4din=p_child_prog_rcf%tracer)
    ENDIF

    IF (iforcing == inwp) THEN
      IF(l_parallel) CALL exchange_data(p_pp%comm_pat_c, phdiag_lp)
      CALL interpol_scal_nudging (p_pp, p_int, p_grf%p_dom(i_chidx), i_chidx, 0, &
                                  1, 1, f3din1=phdiag_lp, f3dout1=phdiag_chi, overshoot_fac=1.005_wp )
      CALL sync_patch_array(SYNC_C,p_pc,phdiag_chi)
    ENDIF

    IF (atm_phy_nwp_config(jg)%inwp_surface == 1) THEN
      IF(l_parallel) CALL exchange_data(p_pp%comm_pat_c, lndvars_lp)
      CALL interpol_scal_nudging (p_pp, p_int, p_grf%p_dom(i_chidx), i_chidx, 0, &
                                  1, 1, f3din1=lndvars_lp, f3dout1=lndvars_chi, overshoot_fac=1.005_wp )
      CALL sync_patch_array(SYNC_C,p_pc,lndvars_chi)
    ENDIF

    IF (atm_phy_nwp_config(jg)%inwp_surface == 1 .AND. lseaice) THEN
      IF(l_parallel) CALL exchange_data(p_pp%comm_pat_c, wtrvars_lp)
      CALL interpol_scal_nudging (p_pp, p_int, p_grf%p_dom(i_chidx), i_chidx, 0, &
                                  1, 1, f3din1=wtrvars_lp, f3dout1=wtrvars_chi, overshoot_fac=1.005_wp )
      CALL sync_patch_array(SYNC_C,p_pc,wtrvars_chi)
    ENDIF

    ! Convert SMI back to wsoil
    CALL smi_to_wsoil(p_patch(jgc), lndvars_chi(:,1:nlev_soil,:))

    ! Step 3: Add reference state to thermodynamic variables and copy land fields
    ! from the container arrays to the prognostic variables (for the time being,
    ! all tiles are initialized with the interpolated aggregated variables)
    ! Note that nest boundary points have to be included here to obtain a proper
    ! initialization of all grid points

    i_nchdom = MAX(1,p_pc%n_childdom)

!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)

    ! cell-based variables
    i_startblk = p_pc%cells%start_blk(1,1)
    i_endblk   = p_pc%cells%end_blk(min_rlcell,i_nchdom)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk,jt,jk1) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_pc, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, 1, min_rlcell)

      ! rho and theta currently contain perturbation fields; thus, the reference
      ! state still needs to be added. Then, rhotheta and exner are diagnosed
      DO jk = 1, nlev_c
        DO jc = i_startidx, i_endidx
          p_child_prog%rho(jc,jk,jb) = &
            p_child_prog%rho(jc,jk,jb) + p_child_metrics%rho_ref_mc(jc,jk,jb)

          p_child_prog%theta_v(jc,jk,jb) = &
            p_child_prog%theta_v(jc,jk,jb) + p_child_metrics%theta_ref_mc(jc,jk,jb)

          p_child_prog%exner(jc,jk,jb) = EXP(rd_o_cvd*LOG(rd/p0ref* &
            p_child_prog%rho(jc,jk,jb) * p_child_prog%theta_v(jc,jk,jb)))

          ! exner_old also needs to be initialized here
          p_child_diag%exner_old(jc,jk,jb) = p_child_prog%exner(jc,jk,jb)
        ENDDO
      ENDDO

      ! Reference temperature at lowest level (taken as proxy for the surface)
      DO jc = i_startidx, i_endidx
        tsfc_ref_c(jc,jb) = p_child_metrics%theta_ref_mc(jc,nlev_c,jb) * &
                            p_child_metrics%exner_ref_mc(jc,nlev_c,jb)
      ENDDO

      IF (iforcing == inwp) THEN
        DO jc = i_startidx, i_endidx
          prm_diag(jgc)%tot_prec(jc,jb)       = MAX(0._wp,phdiag_chi(jc,1,jb))
          prm_diag(jgc)%rain_gsp(jc,jb)       = MAX(0._wp,phdiag_chi(jc,2,jb))
          prm_diag(jgc)%snow_gsp(jc,jb)       = MAX(0._wp,phdiag_chi(jc,3,jb))
          prm_diag(jgc)%rain_con(jc,jb)       = MAX(0._wp,phdiag_chi(jc,4,jb))
          prm_diag(jgc)%snow_con(jc,jb)       = MAX(0._wp,phdiag_chi(jc,5,jb))
          prm_diag(jgc)%rain_gsp_rate(jc,jb)  = phdiag_chi(jc,6,jb)
          prm_diag(jgc)%snow_gsp_rate(jc,jb)  = phdiag_chi(jc,7,jb)
          prm_diag(jgc)%rain_con_rate(jc,jb)  = phdiag_chi(jc,8,jb)
          prm_diag(jgc)%snow_con_rate(jc,jb)  = phdiag_chi(jc,9,jb)
          prm_diag(jgc)%gz0(jc,jb)            = phdiag_chi(jc,10,jb)
          prm_diag(jgc)%tcm(jc,jb)            = phdiag_chi(jc,11,jb)
          prm_diag(jgc)%tch(jc,jb)            = phdiag_chi(jc,12,jb)
          prm_diag(jgc)%tfm(jc,jb)            = phdiag_chi(jc,13,jb)
          prm_diag(jgc)%tfh(jc,jb)            = phdiag_chi(jc,14,jb)
          prm_diag(jgc)%tfv(jc,jb)            = phdiag_chi(jc,15,jb)
          prm_diag(jgc)%t_2m(jc,jb)           = phdiag_chi(jc,16,jb)
          prm_diag(jgc)%qv_2m(jc,jb)          = phdiag_chi(jc,17,jb)
          prm_diag(jgc)%td_2m(jc,jb)          = phdiag_chi(jc,18,jb)
          prm_diag(jgc)%rh_2m(jc,jb)          = phdiag_chi(jc,19,jb)
          prm_diag(jgc)%u_10m(jc,jb)          = phdiag_chi(jc,20,jb)
          prm_diag(jgc)%v_10m(jc,jb)          = phdiag_chi(jc,21,jb)
        ENDDO
      ENDIF

      IF (atm_phy_nwp_config(jgc)%inwp_surface == 1) THEN
        ! Distribute soil variables
        DO jt = 1, ntiles_total
          DO jk = 1, nlev_soil
            jk1 = jk + nlev_soil
            DO jc = i_startidx, i_endidx
              p_child_lprog%w_so_t(jc,jk,jb,jt)  = MAX(0._wp,lndvars_chi(jc,jk,jb))
              p_child_lprog%t_so_t(jc,jk,jb,jt)  = lndvars_chi(jc,jk1,jb) + tsfc_ref_c(jc,jb)
              ! limit w_so_t to pore volume and dryness point - TERRA crashes otherwise
              IF (ext_data(jgc)%atm%soiltyp(jc,jb) == 3) THEN
                p_child_lprog%w_so_t(jc,jk,jb,jt) = MAX(dzsoil(jk)*0.012_wp, &
                  MIN(dzsoil(jk)*0.364_wp,p_child_lprog%w_so_t(jc,jk,jb,jt)))
              ELSE IF (ext_data(jgc)%atm%soiltyp(jc,jb) == 4) THEN
                p_child_lprog%w_so_t(jc,jk,jb,jt) = MAX(dzsoil(jk)*0.030_wp, &
                  MIN(dzsoil(jk)*0.445_wp,p_child_lprog%w_so_t(jc,jk,jb,jt)))
              ELSE IF (ext_data(jgc)%atm%soiltyp(jc,jb) == 5) THEN
                p_child_lprog%w_so_t(jc,jk,jb,jt) = MAX(dzsoil(jk)*0.035_wp, &
                  MIN(dzsoil(jk)*0.455_wp,p_child_lprog%w_so_t(jc,jk,jb,jt)))
              ELSE IF (ext_data(jgc)%atm%soiltyp(jc,jb) == 6) THEN
                p_child_lprog%w_so_t(jc,jk,jb,jt) = MAX(dzsoil(jk)*0.060_wp, &
                  MIN(dzsoil(jk)*0.475_wp,p_child_lprog%w_so_t(jc,jk,jb,jt)))
              ELSE IF (ext_data(jgc)%atm%soiltyp(jc,jb) == 7) THEN
                p_child_lprog%w_so_t(jc,jk,jb,jt) = MAX(dzsoil(jk)*0.065_wp, &
                  MIN(dzsoil(jk)*0.507_wp,p_child_lprog%w_so_t(jc,jk,jb,jt)))
              ELSE IF (ext_data(jgc)%atm%soiltyp(jc,jb) == 8) THEN
                p_child_lprog%w_so_t(jc,jk,jb,jt) = MAX(dzsoil(jk)*0.098_wp, &
                  MIN(dzsoil(jk)*0.863_wp,p_child_lprog%w_so_t(jc,jk,jb,jt)))
              ENDIF
              p_child_lprog2%t_so_t(jc,jk,jb,jt) = p_child_lprog%t_so_t(jc,jk,jb,jt)
              p_child_lprog2%w_so_t(jc,jk,jb,jt) = p_child_lprog%w_so_t(jc,jk,jb,jt)
            ENDDO
          ENDDO
        ENDDO

        jk1 = 2*nlev_soil + 1
        DO jt = 1, ntiles_total
          DO jc = i_startidx, i_endidx
            p_child_lprog%t_so_t(jc,nlev_soil+1,jb,jt) = lndvars_chi(jc,jk1,jb) + &
              tsfc_ref_c(jc,jb)
            ! here we need to initialize t_g and t_g_t because t_g is copied to t_g_t
            ! in nwp_surface_init for land points
            p_child_lprog%t_g(jc,jb) = lndvars_chi(jc,jk1+1,jb) + tsfc_ref_c(jc,jb)
            p_child_lprog2%t_g(jc,jb) = p_child_lprog%t_g(jc,jb)
            p_child_lprog%t_g_t(jc,jb,jt) = p_child_lprog%t_g(jc,jb)
            p_child_lprog2%t_g_t(jc,jb,jt) = p_child_lprog%t_g(jc,jb)
            p_child_lprog%t_s_t(jc,jb,jt) = lndvars_chi(jc,jk1+2,jb) + tsfc_ref_c(jc,jb)
            p_child_lprog2%t_s_t(jc,jb,jt) = p_child_lprog%t_s_t(jc,jb,jt)
            p_child_lprog%t_snow_t(jc,jb,jt) = lndvars_chi(jc,jk1+3,jb) + tsfc_ref_c(jc,jb)
            p_child_lprog%w_snow_t(jc,jb,jt) = MAX(0._wp,lndvars_chi(jc,jk1+4,jb))
            p_child_lprog%rho_snow_t(jc,jb,jt) = lndvars_chi(jc,jk1+5,jb)
            IF (p_child_lprog%rho_snow_t(jc,jb,jt) < 0.75_wp*crhosminf .OR. &
                p_child_lprog%w_snow_t(jc,jb,jt) < 1.e-6_wp) p_child_lprog%rho_snow_t(jc,jb,jt) = 250._wp
            p_child_lprog%w_i_t(jc,jb,jt) = MAX(0._wp,lndvars_chi(jc,jk1+6,jb))
            p_child_ldiag%freshsnow_t(jc,jb,jt) = MAX(0._wp,MIN(1._wp,lndvars_chi(jc,jk1+7,jb)))
            p_child_ldiag%t_seasfc(jc,jb) = lndvars_chi(jc,jk1+8,jb)
            p_child_ldiag%qv_s(jc,jb)     = lndvars_chi(jc,jk1+9,jb)
          ENDDO
        ENDDO
        DO jt = ntiles_total+1, ntiles_total+ntiles_water
          DO jc = i_startidx, i_endidx
            p_child_lprog%t_g_t(jc,jb,jt)  = p_child_lprog%t_g(jc,jb)
            p_child_lprog2%t_g_t(jc,jb,jt) = p_child_lprog%t_g(jc,jb)
            p_child_lprog%t_s_t(jc,jb,jt)  = p_child_ldiag%t_s(jc,jb)
            p_child_lprog2%t_s_t(jc,jb,jt) = p_child_ldiag%t_s(jc,jb)
            p_child_ldiag%qv_s_t(jc,jb,jt) = p_child_ldiag%qv_s(jc,jb)
          ENDDO
        ENDDO

      ENDIF


      IF (atm_phy_nwp_config(jgc)%inwp_surface == 1 .AND. lseaice) THEN
        DO jc = i_startidx, i_endidx
          p_child_wprog%t_ice(jc,jb)     = wtrvars_chi(jc,1,jb)
          p_child_wprog%h_ice(jc,jb)     = wtrvars_chi(jc,2,jb)
          p_child_wprog%t_snow_si(jc,jb) = wtrvars_chi(jc,3,jb)
          p_child_wprog%h_snow_si(jc,jb) = wtrvars_chi(jc,4,jb)
          p_child_ldiag%fr_seaice(jc,jb) = MAX(0._wp,MIN(1._wp,wtrvars_chi(jc,5,jb)))
          IF (p_child_ldiag%fr_seaice(jc,jb) < frsi_min )         p_child_ldiag%fr_seaice(jc,jb) = 0._wp
          IF (p_child_ldiag%fr_seaice(jc,jb) > (1._wp-frsi_min) ) p_child_ldiag%fr_seaice(jc,jb) = 1._wp
        ENDDO
      ENDIF

      IF (atm_phy_nwp_config(jgc)%inwp_surface == 1 .AND. llake) THEN

        CALL flake_coldinit(                                        &
          &   nflkgb      = ext_data(jgc)%atm%fp_count    (jb),     &  ! in
          &   idx_lst_fp  = ext_data(jgc)%atm%idx_lst_fp(:,jb),     &  ! in
          &   depth_lk    = ext_data(jgc)%atm%depth_lk  (:,jb),     &  ! in
          ! here, a proper estimate of the lake surface temperature is required
          &   tskin       = lndvars_chi(:,num_lndvars,jb),          &  ! in
          &   t_snow_lk_p = p_child_wprog%t_snow_lk(:,jb),          &
          &   h_snow_lk_p = p_child_wprog%h_snow_lk(:,jb),          &
          &   t_ice_p     = p_child_wprog%t_ice    (:,jb),          &
          &   h_ice_p     = p_child_wprog%h_ice    (:,jb),          &
          &   t_mnw_lk_p  = p_child_wprog%t_mnw_lk (:,jb),          &
          &   t_wml_lk_p  = p_child_wprog%t_wml_lk (:,jb),          & 
          &   t_bot_lk_p  = p_child_wprog%t_bot_lk (:,jb),          &
          &   c_t_lk_p    = p_child_wprog%c_t_lk   (:,jb),          &
          &   h_ml_lk_p   = p_child_wprog%h_ml_lk  (:,jb),          &
          &   t_b1_lk_p   = p_child_wprog%t_b1_lk  (:,jb),          &
          &   h_b1_lk_p   = p_child_wprog%h_b1_lk  (:,jb),          &
          &   t_g_lk_p    = p_child_lprog%t_g_t    (:,jb,isub_lake) )

      ENDIF

    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    CALL rbf_vec_interpol_cell(p_child_prog%vn, p_patch(jgc), p_int_state(jgc),&
                               p_child_diag%u, p_child_diag%v)

    CALL diagnose_pres_temp(p_child_metrics, p_child_prog, p_child_prog_rcf, p_child_diag, &
                            p_patch(jgc), opt_calc_temp=.TRUE., opt_calc_pres=.TRUE.,      &
                            lnd_prog=p_child_lprog)


    DEALLOCATE(thv_pr_par, rho_pr_par, lndvars_par, wtrvars_par, phdiag_par, lndvars_chi, &
      &       wtrvars_chi, phdiag_chi, tsfc_ref_p, tsfc_ref_c, vn_lp, w_lp, thv_pr_lp   , &
      &       rho_pr_lp, phdiag_lp, tracer_lp, lndvars_lp, wtrvars_lp)

  END SUBROUTINE initialize_nest


  !-------------
  !>
  !! SUBROUTINE interpolate_increments
  !!
  !! Driver routine for interpolating data assimilation increments from a parent domain to a child domain
  !!
  !!
  !! @par Revision History
  !! Initial version by Guenther Zaengl, DWD(2014-07-16)
  !!
  !!
  SUBROUTINE interpolate_increments(initicon, jg, jgc)

    TYPE(t_initicon_state), TARGET, INTENT(INOUT) :: initicon(:)

    INTEGER, INTENT(IN) :: jg   ! parent (source) domain ID
    INTEGER, INTENT(IN) :: jgc  ! child  (target) domain ID

    ! local pointers

    TYPE(t_patch),         POINTER  :: p_pp, p_pc
    TYPE(t_gridref_state), POINTER  :: p_grf
    TYPE(t_int_state),     POINTER  :: p_int

    ! Local arrays for variables living on the local parent grid in the MPI case.
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)   :: vn_lp, temp_lp, pres_lp, qv_lp

    INTEGER :: i_chidx, rl_start, rl_end, i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER :: jb, jk, je
    INTEGER :: nlev_p, nshift

    INTEGER, POINTER :: iidx(:,:,:), iblk(:,:,:)

    LOGICAL :: l_parallel


    IF (.NOT. my_process_is_mpi_parallel()) THEN
      l_parallel = .FALSE.
    ELSE
      l_parallel = .TRUE.
    ENDIF


    p_pp           => p_patch_local_parent(jgc)
    p_pc           => p_patch(jgc)
    p_grf          => p_grf_state_local_parent(jgc)
    p_int          => p_int_state_local_parent(jgc)

    iidx           => p_patch(jg)%edges%cell_idx
    iblk           => p_patch(jg)%edges%cell_blk

    i_chidx  = p_pc%parent_child_index

    ! number of full levels of parent domain
    nlev_p   = p_patch(jg)%nlev

    ! shift between upper model boundaries
    nshift = p_pc%nshift

    IF (msg_level >= 7) THEN
      WRITE(message_text,'(a,i2,a,i2)') 'Interpolation of DA increments, domain ',jg,' =>',jgc
      CALL message('interpolate_increments',message_text)
    ENDIF

    ALLOCATE(vn_lp   (nproma, nlev_p, p_pp%nblks_e),  &
             temp_lp (nproma, nlev_p, p_pp%nblks_c),  &
             pres_lp (nproma, nlev_p, p_pp%nblks_c),  &
             qv_lp   (nproma, nlev_p, p_pp%nblks_c)   )


    ! Preparation: interpolate u and v increments to vn
    !
    rl_start = 2
    rl_end   = min_rledge_int - 2
    i_startblk = p_patch(jg)%edges%start_block(rl_start)
    i_endblk   = p_patch(jg)%edges%end_block(rl_end)

!$OMP PARALLEL DO PRIVATE(jb,jk,je,i_startidx,i_endidx)
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(p_patch(jg), jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)

      DO jk = 1, nlev_p
        DO je = i_startidx, i_endidx

          initicon(jg)%atm_inc%vn(je,jk,jb) = p_int_state(jg)%c_lin_e(je,1,jb)       &
            &               *(initicon(jg)%atm_inc%u(iidx(je,jb,1),jk,iblk(je,jb,1)) &
            &               * p_patch(jg)%edges%primal_normal_cell(je,jb,1)%v1       &
            &               + initicon(jg)%atm_inc%v(iidx(je,jb,1),jk,iblk(je,jb,1)) &
            &               * p_patch(jg)%edges%primal_normal_cell(je,jb,1)%v2)      &
            &               + p_int_state(jg)%c_lin_e(je,2,jb)                       &
            &               *(initicon(jg)%atm_inc%u(iidx(je,jb,2),jk,iblk(je,jb,2)) & 
            &               * p_patch(jg)%edges%primal_normal_cell(je,jb,2)%v1       &
            &               + initicon(jg)%atm_inc%v(iidx(je,jb,2),jk,iblk(je,jb,2)) &
            &               * p_patch(jg)%edges%primal_normal_cell(je,jb,2)%v2  )

        ENDDO  ! je
      ENDDO  ! jk

    ENDDO  ! jb
!$OMP END PARALLEL DO

    CALL sync_patch_array(SYNC_E,p_patch(jg),initicon(jg)%atm_inc%vn)

    ! Parent-to-child interpolation of atmospheric increments
    !
    ! Step 1: boundary interpolation
    !
    CALL interpol2_vec_grf (p_patch(jg), p_pc, p_grf_state(jg)%p_dom(i_chidx), 1, &
      initicon(jg)%atm_inc%vn,        initicon(jgc)%atm_inc%vn                    )

    CALL interpol_scal_grf (p_patch(jg), p_pc, p_grf_state(jg)%p_dom(i_chidx), 3, &
      initicon(jg)%atm_inc%temp,      initicon(jgc)%atm_inc%temp,                 &
      initicon(jg)%atm_inc%pres,      initicon(jgc)%atm_inc%pres,                 &
      initicon(jg)%atm_inc%qv,        initicon(jgc)%atm_inc%qv                    )

    ! Step 2: Interpolation of fields in the model interior

    ! Step 2a: Copy prognostic variables from parent grid to fields on local parent grid
    ! (trivial without MPI parallelization, but communication call needed for MPI)

    CALL exchange_data(p_pp%comm_pat_glb_to_loc_e, RECV=vn_lp, SEND=initicon(jg)%atm_inc%vn)

    CALL exchange_data_mult(p_pp%comm_pat_glb_to_loc_c, 3, 3*nlev_p,         &
      &                     RECV1=temp_lp, SEND1=initicon(jg)%atm_inc%temp,  &
      &                     RECV2=pres_lp, SEND2=initicon(jg)%atm_inc%pres,  &
      &                     RECV3=qv_lp,   SEND3=initicon(jg)%atm_inc%qv     )

    ! Step 2b: Synchronize variables on local parent grids. Note that the sync routines cannot be used in this case
    IF (l_parallel) THEN
      CALL exchange_data(p_pp%comm_pat_e, vn_lp)
      CALL exchange_data_mult(p_pp%comm_pat_c, 3, 3*nlev_p,            &
                              recv1=temp_lp, recv2=pres_lp, recv3=qv_lp)
    ENDIF

    ! Step 2c: Perform interpolation from local parent to child grid
    CALL interpol_vec_nudging (p_pp, p_pc, p_int, p_grf%p_dom(i_chidx), p_grf_state(jgc), &
                               i_chidx, nshift, 1, vn_lp, initicon(jgc)%atm_inc%vn       )
    CALL sync_patch_array(SYNC_E,p_pc,initicon(jgc)%atm_inc%vn)

    CALL interpol_scal_nudging (p_pp, p_int, p_grf%p_dom(i_chidx), i_chidx, nshift, 3, 1, &
                                f3din1=temp_lp, f3dout1=initicon(jgc)%atm_inc%temp,       &
                                f3din2=pres_lp, f3dout2=initicon(jgc)%atm_inc%pres,       &
                                f3din3=qv_lp,   f3dout3=initicon(jgc)%atm_inc%qv          )
    CALL sync_patch_array_mult(SYNC_C,p_pc,3,initicon(jgc)%atm_inc%temp,initicon(jgc)%atm_inc%pres,  &
                               initicon(jgc)%atm_inc%qv)

    DEALLOCATE(vn_lp, temp_lp, pres_lp, qv_lp)

  END SUBROUTINE interpolate_increments

  !-------------
  !>
  !! SUBROUTINE interpolate_sfcana
  !!
  !! Driver routine for interpolating surface analysis data from a parent domain to a child domain
  !! The routine is supposed to work incombination with incremental analysis update only;
  !! it processes sst (i.e. t_so(0) over sea points only), w_so increments, fr_ice, w_snow, rho_snow,
  !! h_snow and freshsnow
  !!
  !!
  !! @par Revision History
  !! Initial version by Guenther Zaengl, DWD(2014-12-12)
  !!
  !!
  SUBROUTINE interpolate_sfcana(initicon, jg, jgc)

    TYPE(t_initicon_state), TARGET, INTENT(INOUT) :: initicon(:)

    INTEGER, INTENT(IN) :: jg   ! parent (source) domain ID
    INTEGER, INTENT(IN) :: jgc  ! child  (target) domain ID

    ! local pointers
    TYPE(t_patch),         POINTER  :: p_pp, p_pc
    TYPE(t_gridref_state), POINTER  :: p_grf
    TYPE(t_int_state),     POINTER  :: p_int

    TYPE(t_lnd_prog),   POINTER     :: p_parent_lprog
    TYPE(t_lnd_prog),   POINTER     :: p_child_lprog
    TYPE(t_lnd_diag),   POINTER     :: p_parent_ldiag
    TYPE(t_lnd_diag),   POINTER     :: p_child_ldiag

    ! local variables

    ! Indices
    INTEGER :: jb, jc, jk, jk1, jt, i_chidx, i_startblk, i_endblk, &
               i_startidx, i_endidx
    INTEGER :: rl_start, rl_end

    INTEGER :: num_lndvars

    ! Local arrays for variables living on the local parent grid in the MPI case. These have
    ! to be allocatable because their dimensions differ between MPI and non-MPI runs
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)   :: lndvars_lp

    ! Local arrays on the parent or child grid. These would not have to be allocatable,
    ! but the computational overhead does not matter for an initialization routine
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)   :: lndvars_par, lndvars_chi

    INTEGER :: i_count, ic

    LOGICAL :: l_parallel

    !-----------------------------------------------------------------------


    IF (msg_level >= 7) THEN
      WRITE(message_text,'(a,i2,a,i2)') 'Interpolation of surface analysis data, domain ',jg,' =>',jgc
      CALL message('interpolate_sfcana',message_text)
    ENDIF

    IF (.NOT. my_process_is_mpi_parallel()) THEN
      l_parallel = .FALSE.
    ELSE
      l_parallel = .TRUE.
    ENDIF

    p_parent_lprog    => p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))
    p_child_lprog     => p_lnd_state(jgc)%prog_lnd(nnow_rcf(jgc))
    p_parent_ldiag    => p_lnd_state(jg)%diag_lnd
    p_child_ldiag     => p_lnd_state(jgc)%diag_lnd

    p_pc              => p_patch(jgc)

    p_grf => p_grf_state_local_parent(jgc)
    p_int => p_int_state_local_parent(jgc)
    p_pp  => p_patch_local_parent(jgc)

    i_chidx  = p_pc%parent_child_index

    ! Number of variables to be interpolated
    num_lndvars = nlev_soil + 5
    
    ALLOCATE(lndvars_par (nproma, num_lndvars, p_patch(jg)%nblks_c), &
             lndvars_chi (nproma, num_lndvars, p_patch(jgc)%nblks_c) )

    ALLOCATE(lndvars_lp (nproma, num_lndvars, p_pp%nblks_c) )

    IF (p_test_run) THEN
      lndvars_par  = 0._wp
      lndvars_chi  = 0._wp
    ENDIF


    ! Step 1: boundary interpolation
    !
    ! Note: the sepration between boundary interpolation and the interpolation
    ! of the fields in the prognostic region of the nest is needed because the
    ! existing interpolation routines are runtime-optimized for operations that
    ! are needed regularly. Boundary interpolation is executed from the parent grid
    ! to the child grid, whereas the so-called nudging interpolation works from
    ! the local parent grid to the child grid.


    ! Step 1a: prepare boundary interpolation

!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)

    ! cell-based variables
    i_startblk = p_patch(jg)%cells%start_block(grf_bdywidth_c+1)
    i_endblk   = p_patch(jg)%cells%end_block(min_rlcell)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk,jk1) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch(jg), jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, grf_bdywidth_c+1, min_rlcell)

      DO jk = 1, nlev_soil
        DO jc = i_startidx, i_endidx
          lndvars_par(jc,jk,jb)  = initicon(jg)%sfc_inc%w_so(jc,jk,jb)
        ENDDO
      ENDDO

      jk1 = nlev_soil
      DO jc = i_startidx, i_endidx
        lndvars_par(jc,jk1+1,jb) = initicon(jg)%sfc%sst(jc,jb)
        lndvars_par(jc,jk1+2,jb) = p_parent_ldiag%fr_seaice(jc,jb)
        lndvars_par(jc,jk1+3,jb) = p_parent_lprog%w_snow_t(jc,jb,1)
        lndvars_par(jc,jk1+4,jb) = MAX(crhosminf,p_parent_lprog%rho_snow_t(jc,jb,1))
        lndvars_par(jc,jk1+5,jb) = p_parent_ldiag%freshsnow_t(jc,jb,1)
      ENDDO

    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL


    ! Step 1b: execute boundary interpolation
    CALL sync_patch_array(SYNC_C,p_patch(jg),lndvars_par)
    CALL interpol_scal_grf (p_patch(jg), p_pc, p_grf_state(jg)%p_dom(i_chidx), 1,  &
      lndvars_par, lndvars_chi, lnoshift=.TRUE.                 )


    ! Step 2: Interpolation of fields in the model interior

    ! Step 2a: Copy variables from parent grid to fields on feedback-parent grid
    ! (trivial without MPI parallelization, but communication call needed for MPI)

    CALL exchange_data(p_pp%comm_pat_glb_to_loc_c, RECV=lndvars_lp, SEND=lndvars_par)

    ! Step 2b: Perform interpolation from local parent to child grid

    ! Note: the sync routines cannot be used for the synchronization on the local
    ! parent grid

    IF(l_parallel) CALL exchange_data(p_pp%comm_pat_c, lndvars_lp)
    CALL interpol_scal_nudging (p_pp, p_int, p_grf%p_dom(i_chidx), i_chidx, 0, &
                                1, 1, f3din1=lndvars_lp, f3dout1=lndvars_chi, overshoot_fac=1.005_wp )
    CALL sync_patch_array(SYNC_C,p_pc,lndvars_chi)



    ! Step 3: Add reference state to thermodynamic variables and copy land fields
    ! from the container arrays to the prognostic variables (for the time being,
    ! all tiles are initialized with the interpolated aggregated variables)
    ! Note that nest boundary points have to be included here to obtain a proper
    ! initialization of all grid points


!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)

    ! cell-based variables
    i_startblk = p_pc%cells%start_block(1)
    i_endblk   = p_pc%cells%end_block(min_rlcell)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk,jt,jk1) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_pc, jb, i_startblk, i_endblk, i_startidx, i_endidx, 1, min_rlcell)

      DO jk = 1, nlev_soil
        DO jc = i_startidx, i_endidx
          initicon(jgc)%sfc_inc%w_so(jc,jk,jb) = lndvars_chi(jc,jk,jb)
        ENDDO
      ENDDO

      jk1 = nlev_soil

      DO jc = i_startidx, i_endidx
        initicon(jgc)%sfc%sst(jc,jb)       = lndvars_chi(jc,jk1+1,jb)
        p_child_ldiag%fr_seaice(jc,jb)     = lndvars_chi(jc,jk1+2,jb)
        p_child_lprog%w_snow_t(jc,jb,1)    = lndvars_chi(jc,jk1+3,jb)
        p_child_lprog%rho_snow_t(jc,jb,1)  = lndvars_chi(jc,jk1+4,jb)
        p_child_ldiag%freshsnow_t(jc,jb,1) = lndvars_chi(jc,jk1+5,jb)

        ! diagnose h_snow after interpolation
        p_child_ldiag%h_snow_t(jc,jb,1) = p_child_lprog%w_snow_t(jc,jb,1)/p_child_lprog%rho_snow_t(jc,jb,1)*rhoh2o

        ! set limits and remove small fractions of fr_seaice
        p_child_ldiag%fr_seaice(jc,jb) = MAX(0._wp,MIN(1._wp,p_child_ldiag%fr_seaice(jc,jb)))
        IF (p_child_ldiag%fr_seaice(jc,jb) < frsi_min )         p_child_ldiag%fr_seaice(jc,jb) = 0._wp
        IF (p_child_ldiag%fr_seaice(jc,jb) > (1._wp-frsi_min) ) p_child_ldiag%fr_seaice(jc,jb) = 1._wp
        p_child_ldiag%freshsnow_t(jc,jb,1) = MIN(1._wp,p_child_ldiag%freshsnow_t(jc,jb,1))
      ENDDO

    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    DEALLOCATE(lndvars_par, lndvars_chi, lndvars_lp)

  END SUBROUTINE interpolate_sfcana


  RECURSIVE SUBROUTINE topo_blending_and_fbk(jg)

    INTEGER, INTENT(IN) :: jg

    INTEGER :: jgc, jn


    ! Loop over nested domains
    DO jn = 1, p_patch(jg)%n_childdom

      jgc = p_patch(jg)%child_id(jn)

      CALL topography_blending(p_patch(jg), p_patch(jgc), p_grf_state(jg)%p_dom(jn), jn, &
               ext_data(jg)%atm%topography_c, ext_data(jgc)%atm%topography_c  )

      IF (p_patch(jgc)%n_childdom > 0) &
        CALL topo_blending_and_fbk(jgc)

    ENDDO

    DO jn = 1, p_patch(jg)%n_childdom

      jgc = p_patch(jg)%child_id(jn)

      IF (lfeedback(jgc) .AND. ifeedback_type == 1) THEN
        CALL topography_feedback(p_patch(jg), jn, ext_data(jg)%atm%topography_c, ext_data(jgc)%atm%topography_c )
      ENDIF

    ENDDO

  END SUBROUTINE topo_blending_and_fbk

  !---------------------------------------------------------------------------
  !>
  !! Computes topography blending for nested domains initialized with real 
  !! topography data
  !!
  !! @par Revision History
  !! Initial release by Guenther Zaengl, DWD, (2011-07-05)
  !!
  SUBROUTINE topography_blending(p_pp, p_pc, p_grf, i_chidx, &
               topo_cp, topo_cc)

    ! patch at parent and child level
    TYPE(t_patch),                TARGET, INTENT(IN) :: p_pp, p_pc
    ! grf state is needed at parent level only
    TYPE(t_gridref_single_state), TARGET, INTENT(IN) :: p_grf

    ! child domain index
    INTEGER, INTENT(IN) :: i_chidx

    ! topography on cells: p = parent level, c = child level
    REAL(wp), INTENT(IN),     DIMENSION(:,:) :: topo_cp
    REAL(wp), INTENT(INOUT),  DIMENSION(:,:) :: topo_cc

    ! auxiliary fields needed for SR calls (SR's expect 3D fields)
    REAL(wp), TARGET :: z_topo_cp(nproma,1,p_pp%nblks_c)
    REAL(wp)         :: z_topo_cc(nproma,1,p_pc%nblks_c)

    ! another auxiliary needed to handle MPI parallelization, and a related pointer
    REAL(wp), ALLOCATABLE, TARGET :: z_topo_clp (:,:,:)
    REAL(wp), POINTER             :: ptr_topo_cp(:,:,:)

    TYPE(t_gridref_single_state), POINTER :: ptr_grf
    TYPE(t_patch),                POINTER :: ptr_pp
    TYPE(t_int_state),            POINTER :: ptr_int

    INTEGER  :: jgc, jb, jc, nblks_c
    INTEGER  :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom, i_rlstart, i_rlend
    REAL(wp) :: wfac

    !-------------------------------------------------------------------------

    jgc = p_pc%id

    ptr_grf    => p_grf_state_local_parent(jgc)%p_dom(i_chidx)
    ptr_pp     => p_patch_local_parent(jgc)
    ptr_int    => p_int_state_local_parent(jgc)

    ALLOCATE(z_topo_clp(nproma,1,ptr_pp%nblks_c))
    
    IF (msg_level >= 10) THEN
      WRITE(message_text,'(2(a,i2))') 'topography blending, domain ',&
        p_pp%id,'  => domain ',p_pc%id
      CALL message('', TRIM(message_text))
    ENDIF

    z_topo_cc(:,:,:) = 0._wp

    nblks_c = p_pp%nblks_c

    ! 1.(a) Copy input topography to auxiliary 3D field
    DO jb = 1, nblks_c
      CALL get_indices_c(p_pp, jb, 1, nblks_c, i_startidx, i_endidx, 1)

      DO jc = i_startidx, i_endidx
        z_topo_cp(jc,1,jb) = topo_cp(jc,jb)
      ENDDO
    ENDDO

    ! 1.(b) Copy this auxiliary field to the local parent in case of MPI parallelization
    
    CALL exchange_data(ptr_pp%comm_pat_glb_to_loc_c, RECV=z_topo_clp, SEND=z_topo_cp)
    ptr_topo_cp => z_topo_clp

    IF (my_process_is_mpi_parallel()) THEN
      ! synchronization (CALL sync does not work on local parent)
      CALL exchange_data(ptr_pp%comm_pat_c, ptr_topo_cp)
    END IF

    ! 1.(c) Interpolate coarse topography on fine mesh

    ! Comment: one could avoid the awkward procedure of interpolating boundary zone and
    ! prognostic zone separately by introducing appropriate interpolation routines and communication
    ! patterns (both of which are optimized for runtime tasks). However, those routines
    ! would not be needed anywhere else, and terrain blending is not runtime-critical

    ! Lateral boundary zone
    CALL interpol_scal_grf (p_pp, p_pc, p_grf, 1, f3din1=z_topo_cp, f3dout1=z_topo_cc, lnoshift=.TRUE.)

    ! Prognostic part of the model domain
    ! Note: in contrast to boundary interpolation, nudging expects the input on the local parent grid
    CALL interpol_scal_nudging (ptr_pp, ptr_int, ptr_grf, i_chidx, 0, 1, 1, &
                                f3din1=ptr_topo_cp, f3dout1=z_topo_cc       )

    ! 2. Apply terrain blending to cell points
    ! In the boundary interpolation zone (cell rows 1-4), the interpolated values
    ! are taken; in the subsequent 8 rows, linear blending is applied


    i_nchdom = MAX(1,p_pc%n_childdom)

    ! 2.(a) Pure interpolation in boundary interpolation zone
    i_rlstart  = 1
    i_rlend    = grf_bdywidth_c
    i_startblk = p_pc%cells%start_blk(i_rlstart,1)
    i_endblk   = p_pc%cells%end_blk(i_rlend,i_nchdom)

    DO jb = i_startblk, i_endblk
      CALL get_indices_c(p_pc, jb, i_startblk, i_endblk,         &
                         i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jc = i_startidx, i_endidx
        topo_cc(jc,jb) = z_topo_cc(jc,1,jb)
      ENDDO
    ENDDO

    ! 2.(b) Blending in prognostic part of the nested domain
    i_rlstart  = grf_bdywidth_c + 1
    i_rlend    = min_rlcell
    i_startblk = p_pc%cells%start_blk(i_rlstart,1)
    i_endblk   = p_pc%cells%end_blk(i_rlend,i_nchdom)

    DO jb = i_startblk, i_endblk
      CALL get_indices_c(p_pc, jb, i_startblk, i_endblk,         &
                         i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jc = i_startidx, i_endidx

        IF (p_pc%cells%refin_ctrl(jc,jb) > 0 .AND. &
            p_pc%cells%refin_ctrl(jc,jb) <= grf_bdywidth_c + nudge_zone_width) THEN

          wfac = REAL(p_pc%cells%refin_ctrl(jc,jb) - grf_bdywidth_c,wp)/ &
            REAL(nudge_zone_width+1,wp)
          topo_cc(jc,jb) = wfac*topo_cc(jc,jb) + (1._wp-wfac)*z_topo_cc(jc,1,jb)

        ENDIF ! on the interior grid points, topo_cc remains unchanged

      ENDDO
    ENDDO

    CALL sync_patch_array(SYNC_C,p_pc,topo_cc)

    DEALLOCATE(z_topo_clp)

  END SUBROUTINE topography_blending


  !---------------------------------------------------------------------------
  !>
  !! Computes topography feedback for nested domains initialized with real 
  !! topography data
  !!
  !! @par Revision History
  !! Initial release by Guenther Zaengl, DWD, (2011-07-05)
  !!
  SUBROUTINE topography_feedback(p_pp, i_chidx, topo_cp, topo_cc)

    ! patch at parent level
    TYPE(t_patch),                TARGET, INTENT(IN) :: p_pp

    ! child domain index
    INTEGER, INTENT(IN) :: i_chidx

    ! topography on cells: p = parent level, c = child level
    REAL(wp), INTENT(IN),    DIMENSION(:,:)         :: topo_cc
    REAL(wp), INTENT(INOUT), DIMENSION(:,:), TARGET :: topo_cp

    ! auxiliary field needed for SR calls (SR's expect 3D fields)
    REAL(wp) :: z_topo_cp(nproma,1,p_pp%nblks_c)

    ! another auxiliary to handle MPI parallelization, and related pointer
    REAL(wp), ALLOCATABLE, TARGET :: z_topo_clp(:,:)
    REAL(wp), POINTER             :: ptr_topo_cp(:,:)

    TYPE(t_gridref_state), POINTER :: ptr_grf
    TYPE(t_patch),         POINTER :: ptr_pp

    INTEGER  :: jgc, jb, jc
    INTEGER  :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom, i_rlstart, i_rlend

    INTEGER, POINTER, DIMENSION(:,:,:) :: iidx, iblk

    !-------------------------------------------------------------------------

    jgc = p_pp%child_id(i_chidx) ! child domain ID

    ALLOCATE(z_topo_clp(nproma,p_patch_local_parent(jgc)%nblks_c))
    ptr_pp      => p_patch_local_parent(jgc)
    ptr_grf     => p_grf_state_local_parent(jgc)
    ptr_topo_cp => z_topo_clp

    IF (msg_level >= 10) THEN
      WRITE(message_text,'(2(a,i2))') 'topography feedback, domain ',&
        p_pp%child_id(i_chidx),'  => domain ',p_pp%id
        CALL message('', TRIM(message_text))
    ENDIF

    ! Start/End block in the parent domain
    i_rlstart  = grf_fbk_start_c
    i_rlend    = min_rlcell_int

    i_startblk = ptr_pp%cells%start_blk(i_rlstart,i_chidx)
    i_endblk   = ptr_pp%cells%end_blk(i_rlend,i_chidx)

    iidx => ptr_pp%cells%child_idx
    iblk => ptr_pp%cells%child_blk

    ! Perform feedback for cell points in the prognostic part of the nest overlap area
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(ptr_pp, jb, i_startblk, i_endblk,       &
                         i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jc = i_startidx, i_endidx
        ptr_topo_cp(jc,jb) =                                                  &
          topo_cc(iidx(jc,jb,1),iblk(jc,jb,1))*ptr_grf%fbk_wgt_bln(jc,jb,1) + &
          topo_cc(iidx(jc,jb,2),iblk(jc,jb,2))*ptr_grf%fbk_wgt_bln(jc,jb,2) + &
          topo_cc(iidx(jc,jb,3),iblk(jc,jb,3))*ptr_grf%fbk_wgt_bln(jc,jb,3) + &
          topo_cc(iidx(jc,jb,4),iblk(jc,jb,4))*ptr_grf%fbk_wgt_bln(jc,jb,4)
      ENDDO
    ENDDO

    ! Attention: the feedback communication pattern is defined on the local parent (ptr_pp), ...
    CALL exchange_data(ptr_pp%comm_pat_loc_to_glb_c_fbk, RECV=topo_cp, SEND=ptr_topo_cp, &
      &                l_recv_exists=.TRUE.)

    ! ... whereas the subsequent halo communication needs to be done on the parent grid (p_pp)
    CALL sync_patch_array(SYNC_C,p_pp,topo_cp)


    ! Copy fed-back topography to z_topo_cp to prepare interpolation to vertices
    i_nchdom = MAX(1,p_pp%n_childdom)
    i_rlstart  = 1
    i_rlend    = min_rlcell
    i_startblk = p_pp%cells%start_blk(i_rlstart,1)
    i_endblk   = p_pp%cells%end_blk(i_rlend,i_nchdom)

    DO jb = i_startblk, i_endblk
      CALL get_indices_c(p_pp, jb, i_startblk, i_endblk,         &
                         i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jc = i_startidx, i_endidx

        z_topo_cp(jc,1,jb) = topo_cp(jc,jb) 

      ENDDO
    ENDDO

    DEALLOCATE(z_topo_clp)

  END SUBROUTINE topography_feedback


END MODULE mo_nh_init_nest_utils

