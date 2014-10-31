!>
!! mo_solve_nonhydro
!!
!! This module contains the nonhydrostatic dynamical core for the triangular version
!! Its routines were previously contained in mo_divergent_modes and mo_vector_operations
!! but have been extracted for better memory efficiency
!!
!! @author Guenther Zaengl, DWD
!!
!! @par Revision History
!! Initial release by Guenther Zaengl (2010-10-13) based on earlier work
!! by Almut Gassmann, MPI-M
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

MODULE mo_solve_nonhydro

  USE mo_kind,                 ONLY: wp, vp
  USE mo_nonhydrostatic_config,ONLY: itime_scheme,iadv_rhotheta, igradp_method, l_open_ubc, &
                                     kstart_moist, lhdiff_rcf, divdamp_fac, divdamp_order,  &
                                     divdamp_type, rayleigh_type, iadv_rcf, rhotheta_offctr,&
                                     veladv_offctr, divdamp_fac_o2
  USE mo_dynamics_config,   ONLY: idiv_method
  USE mo_parallel_config,   ONLY: nproma, p_test_run, itype_comm, use_dycore_barrier, &
    & use_icon_comm
  USE mo_run_config,        ONLY: ltimer, timers_level, lvert_nest
  USE mo_model_domain,      ONLY: t_patch
  USE mo_grid_config,       ONLY: l_limited_area
  USE mo_gridref_config,    ONLY: grf_intmethod_e
  USE mo_interpol_config,   ONLY: nudge_max_coeff
  USE mo_intp_data_strc,    ONLY: t_int_state
  USE mo_intp,              ONLY: cells2verts_scalar
  USE mo_nonhydro_types,    ONLY: t_nh_state
  USE mo_physical_constants,ONLY: cpd, rd, cvd, cvd_o_rd, grav, rd_o_cpd, p0ref
  USE mo_math_gradients,    ONLY: grad_green_gauss_cell
  USE mo_velocity_advection,ONLY: velocity_tendencies
  USE mo_math_constants,    ONLY: pi, dbl_eps
  USE mo_math_divrot,       ONLY: div_avg
  USE mo_vertical_grid,     ONLY: nrdmax, nflat_gradp
  USE mo_nh_init_utils,     ONLY: nflatlev
  USE mo_loopindices,       ONLY: get_indices_c, get_indices_e
  USE mo_impl_constants,    ONLY: min_rlcell_int, min_rledge_int, min_rlvert_int, &
    &                             min_rlcell, RAYLEIGH_CLASSIC, RAYLEIGH_KLEMP
  USE mo_impl_constants_grf,ONLY: grf_bdywidth_c, grf_bdywidth_e
  USE mo_advection_hflux,   ONLY: upwind_hflux_miura3
  USE mo_advection_traj,    ONLY: btraj
  USE mo_sync,              ONLY: SYNC_E, SYNC_C, sync_patch_array, sync_patch_array_mult
  USE mo_mpi,               ONLY: my_process_is_mpi_all_seq, work_mpi_barrier
  USE mo_timer,             ONLY: timer_solve_nh, timer_barrier, timer_start, timer_stop,       &
                                  timer_solve_nh_cellcomp, timer_solve_nh_edgecomp,             &
                                  timer_solve_nh_vnupd, timer_solve_nh_vimpl, timer_solve_nh_exch
  USE mo_icon_comm_lib,     ONLY: icon_comm_sync
  USE mo_vertical_coord_table,ONLY: vct_a
  USE mo_nh_prepadv_types,  ONLY: t_prepare_adv
  USE mo_initicon_config,   ONLY: is_iau_active, iau_wgt_dyn
  IMPLICIT NONE

  PRIVATE


  REAL(wp), PARAMETER :: rd_o_cvd = 1._wp / cvd_o_rd
  REAL(wp), PARAMETER :: cpd_o_rd = 1._wp / rd_o_cpd
  REAL(wp), PARAMETER :: rd_o_p0ref = rd / p0ref
  REAL(wp), PARAMETER :: grav_o_cpd = grav / cpd

  PUBLIC :: solve_nh

  CONTAINS


  !>
  !! solve_nh
  !!
  !! Main solver routine for nonhydrostatic dynamical core
  !!
  !! @par Revision History
  !! Development started by Guenther Zaengl on 2010-02-03
  !!
  SUBROUTINE solve_nh (p_nh, p_patch, p_int, prep_adv, nnow, nnew, l_init, l_recompute, lsave_mflx, &
                       lprep_adv, lclean_mflx, idyn_timestep, jstep, l_bdy_nudge, dtime)

    TYPE(t_nh_state),  TARGET, INTENT(INOUT) :: p_nh
    TYPE(t_int_state), TARGET, INTENT(IN)    :: p_int
    TYPE(t_patch),     TARGET, INTENT(IN)    :: p_patch
    TYPE(t_prepare_adv),       INTENT(INOUT) :: prep_adv

    ! Initialization switch that has to be .TRUE. at the initial time step only (not for restart)
    LOGICAL,                   INTENT(INOUT) :: l_init
    ! Switch to recompute velocity tendencies after a physics call irrespective of the time scheme option
    LOGICAL,                   INTENT(IN)    :: l_recompute
    ! Switch if mass flux needs to be saved for nest boundary interpolation tendency computation
    LOGICAL,                   INTENT(IN)    :: lsave_mflx
    ! Switch if preparations for tracer advection shall be computed
    LOGICAL,                   INTENT(IN)    :: lprep_adv
    ! Switch if mass fluxes computed for tracer advection need to be reinitialized
    LOGICAL,                   INTENT(IN)    :: lclean_mflx
    ! Counter of dynamics time step within a large time step (ranges from 1 to iadv_rcf)
    INTEGER,                   INTENT(IN)    :: idyn_timestep
    ! Time step count since last boundary interpolation (ranges from 0 to 2*iadv_rcf-1)
    INTEGER,                   INTENT(IN)    :: jstep
    ! Switch to determine if boundary nudging is executed
    LOGICAL,                   INTENT(IN)    :: l_bdy_nudge
    ! Time levels
    INTEGER,                   INTENT(IN)    :: nnow, nnew
    ! Time step
    REAL(wp),                  INTENT(IN)    :: dtime

    ! Local variables
    INTEGER  :: jb, jk, jc, je, jks
    INTEGER  :: nlev, nlevp1              !< number of full levels
    INTEGER  :: i_startblk, i_endblk, i_startidx, i_endidx, ishift
    INTEGER  :: rl_start, rl_end, istep, ntl1, ntl2, nvar, nshift, nshift_total
    INTEGER  :: ic, ie, ilc0, ibc0, ikp1, ikp2

    REAL(wp) :: z_theta_v_fl_e  (nproma,p_patch%nlev  ,p_patch%nblks_e), &
                z_theta_v_e     (nproma,p_patch%nlev  ,p_patch%nblks_e), &
                z_rho_e         (nproma,p_patch%nlev  ,p_patch%nblks_e), &
                z_mass_fl_div   (nproma,p_patch%nlev  ,p_patch%nblks_c), & ! used for idiv_method=2 only
                z_theta_v_fl_div(nproma,p_patch%nlev  ,p_patch%nblks_c), & ! used for idiv_method=2 only
                z_dwdz_dd       (nproma,p_patch%nlev,  p_patch%nblks_c), &
                z_exner_pr      (nproma,p_patch%nlev  ,p_patch%nblks_c), &
                z_theta_v_v     (nproma,p_patch%nlev  ,p_patch%nblks_v), & ! used for iadv_rhotheta=1 only
                z_rho_v         (nproma,p_patch%nlev  ,p_patch%nblks_v)    ! used for iadv_rhotheta=1 only

#ifndef __LOOP_EXCHANGE
    REAL(vp) :: z_distv_bary  (nproma,p_patch%nlev  ,p_patch%nblks_e,2)
    INTEGER ::  z_cell_idx    (nproma,p_patch%nlev  ,p_patch%nblks_e)
    INTEGER ::  z_cell_blk    (nproma,p_patch%nlev  ,p_patch%nblks_e)
#endif

    ! The data type vp (variable precision) is by default the same as wp but reduces
    ! to single precision when the __MIXED_PRECISION cpp flag is set at compile time
    REAL(vp) :: z_th_ddz_exner_c(nproma,p_patch%nlev,p_patch%nblks_c), &
                z_dexner_dz_c (2,nproma,p_patch%nlev,p_patch%nblks_c), &
                z_vt_ie         (nproma,p_patch%nlev,p_patch%nblks_e), &
                z_kin_hor_e     (nproma,p_patch%nlev,p_patch%nblks_e), &
                z_exner_ex_pr (nproma,p_patch%nlevp1,p_patch%nblks_c), & ! nlevp1 is intended here
                z_gradh_exner   (nproma,p_patch%nlev,p_patch%nblks_e), &
                z_rth_pr      (2,nproma,p_patch%nlev,p_patch%nblks_c), &
                z_grad_rth    (4,nproma,p_patch%nlev,p_patch%nblks_c), &
                z_w_concorr_me  (nproma,p_patch%nlev,p_patch%nblks_e)

    ! This field in addition has reversed index order (vertical first) for optimization
#ifdef __LOOP_EXCHANGE
    REAL(vp) :: z_graddiv_vn    (p_patch%nlev,nproma,p_patch%nblks_e)
#else
    REAL(vp) :: z_graddiv_vn    (nproma,p_patch%nlev,p_patch%nblks_e)
#endif

    REAL(wp) :: z_w_expl        (nproma,p_patch%nlevp1),          &
                z_thermal_exp   (nproma,p_patch%nblks_c),         &
                z_vn_avg        (nproma,p_patch%nlev  ),          &
                z_mflx_top      (nproma,p_patch%nblks_c),         &
                z_contr_w_fl_l  (nproma,p_patch%nlevp1),          &
                z_rho_expl      (nproma,p_patch%nlev  ),          &
                z_exner_expl    (nproma,p_patch%nlev  ),          &
                z_theta_tavg    (nproma,p_patch%nlev  ),          &
                z_rho_tavg      (nproma,p_patch%nlev  )

    ! The data type vp (variable precision) is by default the same as wp but reduces
    ! to single precision when the __MIXED_PRECISION cpp flag is set at compile time
    REAL(vp) :: z_alpha         (nproma,p_patch%nlevp1),          &
                z_beta          (nproma,p_patch%nlev  ),          &
                z_gamma         (nproma,p_patch%nlev  ),          &
                z_a             (nproma,p_patch%nlevp1),          &
                z_b             (nproma,p_patch%nlevp1),          &
                z_c             (nproma,p_patch%nlev  ),          &
                z_g             (nproma,p_patch%nlev  ),          &
                z_q             (nproma,p_patch%nlev  ),          &
                z_w_backtraj    (nproma,p_patch%nlevp1),          &
                z_graddiv2_vn   (nproma,p_patch%nlev  ),          &
                z_theta_v_pr_mc (nproma,p_patch%nlev  ),          &
                z_theta_v_pr_ic (nproma,p_patch%nlevp1),          &
                z_exner_ic      (nproma,p_patch%nlevp1),          &
                z_w_concorr_mc  (nproma,p_patch%nlev  ),          &
                z_flxdiv_mass   (nproma,p_patch%nlev  ),          &
                z_flxdiv_theta  (nproma,p_patch%nlev  ),          &
                z_hydro_corr    (nproma,p_patch%nblks_e)


    REAL(wp):: z_theta1, z_theta2, z_raylfac, wgt_nnow_vel, wgt_nnew_vel,     &
               dt_shift, wgt_nnow_rth, wgt_nnew_rth, dthalf,                  &
               z_ntdistv_bary(2), distv_bary(2), r_iadv_rcf, scal_divdamp_o2, &
               z1_div3d, z2_div3d, zf

    REAL(wp), DIMENSION(p_patch%nlev) :: scal_divdamp, bdy_divdamp, enh_divdamp_fac, scal_div3d

    INTEGER :: nproma_gradp, nblks_gradp, npromz_gradp, nlen_gradp, jk_start
    LOGICAL :: lcompute, lcleanup, lvn_only, lvn_pos

    ! Local variables to control vertical nesting
    LOGICAL :: l_vert_nested, l_child_vertnest

    ! Pointers
    INTEGER, POINTER   &
#ifdef _CRAYFTN
      , CONTIGUOUS     &
#endif
      ::               &
      ! to cell indices
      icidx(:,:,:), icblk(:,:,:), &
      ! to edge indices
      ieidx(:,:,:), ieblk(:,:,:), &
      ! to vertex indices
      ividx(:,:,:), ivblk(:,:,:), &
      ! to vertical neighbor indices for pressure gradient computation
      ikidx(:,:,:,:),             &
      ! to quad edge indices
      iqidx(:,:,:), iqblk(:,:,:), &
      ! for igradp_method = 3
      iplev(:), ipeidx(:), ipeblk(:)


     !-------------------------------------------------------------------
    IF (use_dycore_barrier) THEN
      CALL timer_start(timer_barrier)
      CALL work_mpi_barrier()
      CALL timer_stop(timer_barrier)
    ENDIF
    !-------------------------------------------------------------------
    
    IF (ltimer) CALL timer_start(timer_solve_nh)

    IF (lvert_nest .AND. (p_patch%nshift_total > 0)) THEN  
      l_vert_nested = .TRUE.
      nshift_total  = p_patch%nshift_total
    ELSE
      l_vert_nested = .FALSE.
      nshift_total  = 0
    ENDIF
    IF (lvert_nest .AND. p_patch%n_childdom > 0 .AND.              &
      (p_patch%nshift_child > 0 .OR. p_patch%nshift_total > 0)) THEN
      l_child_vertnest = .TRUE.
      nshift = p_patch%nshift_child + 1
    ELSE
      l_child_vertnest = .FALSE.
      nshift = 0
    ENDIF
    dthalf  = 0.5_wp*dtime

    ! Inverse value of iadv_rcf for tracer advection precomputations
    r_iadv_rcf = 1._wp/REAL(iadv_rcf,wp)

    ! number of vertical levels
    nlev   = p_patch%nlev
    nlevp1 = p_patch%nlevp1

    ! Set pointers to neighbor cells
    icidx => p_patch%edges%cell_idx
    icblk => p_patch%edges%cell_blk

    ! Set pointers to neighbor edges
    ieidx => p_patch%cells%edge_idx
    ieblk => p_patch%cells%edge_blk

    ! Set pointers to vertices of an edge
    ividx => p_patch%edges%vertex_idx
    ivblk => p_patch%edges%vertex_blk

    ! Set pointer to vertical neighbor indices for pressure gradient
    ikidx => p_nh%metrics%vertidx_gradp

    ! Set pointers to quad edges
    iqidx => p_patch%edges%quad_idx
    iqblk => p_patch%edges%quad_blk

    ! scaling factor for second-order divergence damping: divdamp_fac_o2*delta_x**2
    ! delta_x**2 is approximated by the mean cell area
    scal_divdamp_o2 = divdamp_fac_o2 * p_patch%geometry_info%mean_cell_area 

    ! Fourth-order divergence damping
    !
    ! Impose a minimum value to divergence damping factor that, starting at 20 km, increases linearly 
    ! with height to a value of 0.004 (= the namelist default) at 40 km
    DO jk = 1, nlev
      jks = jk + nshift_total
      zf = 0.5_wp*(vct_a(jks)+vct_a(jks+1))
      IF (divdamp_order == 24) THEN
        enh_divdamp_fac(jk) = MAX( 0._wp, -0.25_wp*divdamp_fac_o2 + MAX(divdamp_fac, &
        MIN(0.004_wp,0.004_wp*(zf-20000._wp)/20000._wp)) )
      ELSE
        enh_divdamp_fac(jk) = MAX(divdamp_fac,MIN(0.004_wp,0.004_wp*(zf-20000._wp)/20000._wp))
      ENDIF
    ENDDO

    scal_divdamp(:) = - enh_divdamp_fac(:) * p_patch%geometry_info%mean_cell_area**2

    ! Settings for using vertical component of divergence for divdamp_type=3/32
    IF (divdamp_type == 3) THEN ! Use 3D divergence everywhere
      scal_div3d(:) = 1._wp
    ELSE IF (divdamp_type == 32) THEN ! Apply blending between 3D divergence damping in the troposphere
      z1_div3d = 12500._wp            ! and 2D divergence damping in the stratosphere
      z2_div3d = 17500._wp
      DO jk = 1, nlev
        jks = jk + nshift_total
        zf = 0.5_wp*(vct_a(jks)+vct_a(jks+1))
        IF (zf >= z2_div3d) THEN
          scal_div3d(jk) = 0._wp
        ELSE IF (zf >= z1_div3d) THEN
          scal_div3d(jk) = (z2_div3d-zf)/(z2_div3d-z1_div3d)
        ELSE
          scal_div3d(jk) = 1._wp
        ENDIF
      ENDDO
    ENDIF

    ! Time increment for backward-shifting of lateral boundary mass flux 
    dt_shift = dtime*REAL(2*iadv_rcf-1,wp)/2._wp

    ! Coefficient for reduced fourth-order divergence damping along nest boundaries
    bdy_divdamp(:) = 0.75_wp/(nudge_max_coeff + dbl_eps)*ABS(scal_divdamp(:))


    IF (p_test_run) THEN
      z_rho_e     = 0._wp
      z_theta_v_e = 0._wp
      z_dwdz_dd   = 0._wp
    ENDIF

    ! Set time levels of ddt_adv fields for call to velocity_tendencies
    IF (itime_scheme >= 4) THEN ! Velocity advection averaging nnow and nnew tendencies
      ntl1 = nnow
      ntl2 = nnew
    ELSE                        ! Velocity advection is taken at nnew only
      ntl1 = 1
      ntl2 = 1
    ENDIF

    ! Weighting coefficients for velocity advection if tendency averaging is used
    ! The off-centering specified here turned out to be beneficial to numerical
    ! stability in extreme situations
    wgt_nnow_vel = 0.5_wp - veladv_offctr ! default value for veladv_offctr is 0.25
    wgt_nnew_vel = 0.5_wp + veladv_offctr

    ! Weighting coefficients for rho and theta at interface levels in the corrector step
    ! This empirically determined weighting minimizes the vertical wind off-centering
    ! needed for numerical stability of vertical sound wave propagation
    wgt_nnew_rth = 0.5_wp + rhotheta_offctr ! default value for rhotheta_offctr is -0.1
    wgt_nnow_rth = 1._wp - wgt_nnew_rth

    DO istep = 1, 2

      IF (istep == 1) THEN ! predictor step
        IF (itime_scheme >= 6 .OR. l_init .OR. l_recompute) THEN
          IF (itime_scheme < 6 .AND. .NOT. l_init) THEN
            lvn_only = .TRUE. ! Recompute only vn tendency
          ELSE
            lvn_only = .FALSE.
          ENDIF
          CALL velocity_tendencies(p_nh%prog(nnow),p_patch,p_int,p_nh%metrics,p_nh%diag,z_w_concorr_me, &
                                   z_kin_hor_e,z_vt_ie,ntl1,istep,lvn_only,dtime)
        ENDIF
        nvar = nnow
      ELSE                 ! corrector step
        lvn_only = .FALSE.
        CALL velocity_tendencies(p_nh%prog(nnew),p_patch,p_int,p_nh%metrics,p_nh%diag,z_w_concorr_me, &
                                 z_kin_hor_e,z_vt_ie,ntl2,istep,lvn_only,dtime)
        nvar = nnew
      ENDIF

      l_init = .FALSE. ! should be .TRUE. only at initial predictor step

    ! Preparations for igradp_method = 3/5 (reformulated extrapolation below the ground)
    IF (istep == 1 .AND. (igradp_method == 3 .OR. igradp_method == 5)) THEN
      
      iplev  => p_nh%metrics%pg_vertidx
      ipeidx => p_nh%metrics%pg_edgeidx
      ipeblk => p_nh%metrics%pg_edgeblk

      nproma_gradp = MIN(nproma,256)
      nblks_gradp  = INT(p_nh%metrics%pg_listdim/nproma_gradp)
      npromz_gradp = MOD(p_nh%metrics%pg_listdim,nproma_gradp)
      IF (npromz_gradp > 0) THEN
        nblks_gradp = nblks_gradp + 1
      ELSE
        npromz_gradp = nproma_gradp
      ENDIF

    ENDIF

    IF (timers_level > 5) CALL timer_start(timer_solve_nh_cellcomp)

    ! Computations on mass points
!$OMP PARALLEL PRIVATE (rl_start,rl_end,i_startblk,i_endblk)

    rl_start = 3
    IF (istep == 1) THEN
      rl_end = min_rlcell_int - 1
    ELSE ! halo points are not needed in step 2
      rl_end = min_rlcell_int
    ENDIF

    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

    ! initialize nest boundary points of z_rth_pr with zero
    IF (istep == 1 .AND. (p_patch%id > 1 .OR. l_limited_area)) THEN
!$OMP WORKSHARE
      z_rth_pr(:,:,:,1:i_startblk) = 0._wp
!$OMP END WORKSHARE
    ENDIF

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,jc,z_exner_ic,z_theta_v_pr_mc,z_theta_v_pr_ic,z_w_backtraj,&
!$OMP            z_rho_tavg,z_theta_tavg) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)

      IF (istep == 1) THEN ! to be executed in predictor step only

        DO jk = 1, nlev
!DIR$ IVDEP
          DO jc = i_startidx, i_endidx
            ! temporally extrapolated perturbation Exner pressure (used for horizontal gradients only)
            z_exner_ex_pr(jc,jk,jb) = - p_nh%metrics%exner_ref_mc(jc,jk,jb) +              &
              (1._wp + p_nh%metrics%exner_exfac(jc,jk,jb))*p_nh%prog(nnow)%exner(jc,jk,jb) &
                     - p_nh%metrics%exner_exfac(jc,jk,jb) *p_nh%diag%exner_old(jc,jk,jb)

            ! non-extrapolated perturbation Exner pressure
            z_exner_pr(jc,jk,jb) = p_nh%prog(nnow)%exner(jc,jk,jb) - &
              p_nh%metrics%exner_ref_mc(jc,jk,jb)

            ! Now save current time level in exner_old
            p_nh%diag%exner_old(jc,jk,jb) = p_nh%prog(nnow)%exner(jc,jk,jb)

          ENDDO
        ENDDO

        ! The purpose of the extra level of exner_pr is to simplify coding for
        ! igradp_method=4/5. It is multiplied with zero and thus actually not used
        z_exner_ex_pr(:,nlevp1,jb) = 0._wp

        IF (l_open_ubc .AND. .NOT. l_vert_nested) THEN
          ! Compute contribution of thermal expansion to vertical wind at model top
          ! Isothermal expansion is assumed
          z_thermal_exp(:,jb) = 0._wp
          DO jk = 1, nlev
!DIR$ IVDEP
            DO jc = i_startidx, i_endidx
              z_thermal_exp(jc,jb) = z_thermal_exp(jc,jb) + cvd_o_rd                      &
                * p_nh%diag%ddt_exner_phy(jc,jk,jb)                                       &
                /  (p_nh%prog(nnow)%exner(jc,jk,jb)*p_nh%metrics%inv_ddqz_z_full(jc,jk,jb))
            ENDDO
          ENDDO
        ENDIF

        IF (igradp_method <= 3) THEN
          ! Perturbation Exner pressure on bottom half level
!DIR$ IVDEP
          DO jc = i_startidx, i_endidx
            z_exner_ic(jc,nlevp1) =                                         &
              p_nh%metrics%wgtfacq_c(jc,1,jb)*z_exner_ex_pr(jc,nlev  ,jb) + &
              p_nh%metrics%wgtfacq_c(jc,2,jb)*z_exner_ex_pr(jc,nlev-1,jb) + &
              p_nh%metrics%wgtfacq_c(jc,3,jb)*z_exner_ex_pr(jc,nlev-2,jb)
          ENDDO

          DO jk = nlev, MAX(2,nflatlev(p_patch%id)), -1
!DIR$ IVDEP
            DO jc = i_startidx, i_endidx
              ! Exner pressure on remaining half levels for metric correction term
              z_exner_ic(jc,jk) =                                              &
                p_nh%metrics%wgtfac_c(jc,jk,jb)*z_exner_ex_pr(jc,jk,jb) +      &
                (1._vp-p_nh%metrics%wgtfac_c(jc,jk,jb))*z_exner_ex_pr(jc,jk-1,jb)

              ! First vertical derivative of perturbation Exner pressure
              z_dexner_dz_c(1,jc,jk,jb) =                    &
               (z_exner_ic(jc,jk) - z_exner_ic(jc,jk+1)) *   &
                p_nh%metrics%inv_ddqz_z_full(jc,jk,jb)
            ENDDO
          ENDDO

          IF (nflatlev(p_patch%id) == 1) THEN
            ! Perturbation Exner pressure on top half level
!DIR$ IVDEP
            DO jc = i_startidx, i_endidx
            z_exner_ic(jc,1) =                                          &
              p_nh%metrics%wgtfacq1_c(jc,1,jb)*z_exner_ex_pr(jc,1,jb) + &
              p_nh%metrics%wgtfacq1_c(jc,2,jb)*z_exner_ex_pr(jc,2,jb) + &
              p_nh%metrics%wgtfacq1_c(jc,3,jb)*z_exner_ex_pr(jc,3,jb)

              ! First vertical derivative of perturbation Exner pressure
              z_dexner_dz_c(1,jc,1,jb) =                    &
               (z_exner_ic(jc,1) - z_exner_ic(jc,2)) *   &
                p_nh%metrics%inv_ddqz_z_full(jc,1,jb)
            ENDDO
          ENDIF

        ENDIF

        z_rth_pr(1,i_startidx:i_endidx,1,jb) =  p_nh%prog(nnow)%rho(i_startidx:i_endidx,1,jb) - &
          p_nh%metrics%rho_ref_mc(i_startidx:i_endidx,1,jb)

        z_rth_pr(2,i_startidx:i_endidx,1,jb) =  p_nh%prog(nnow)%theta_v(i_startidx:i_endidx,1,jb) - &
          p_nh%metrics%theta_ref_mc(i_startidx:i_endidx,1,jb)

        z_theta_v_pr_mc(i_startidx:i_endidx,1) = z_rth_pr(2,i_startidx:i_endidx,1,jb)

        DO jk = 2, nlev
!DIR$ IVDEP
          DO jc = i_startidx, i_endidx
            ! density at interface levels for vertical flux divergence computation
            p_nh%diag%rho_ic(jc,jk,jb) = p_nh%metrics%wgtfac_c(jc,jk,jb)*p_nh%prog(nnow)%rho(jc,jk,jb) + &
              (1._wp-p_nh%metrics%wgtfac_c(jc,jk,jb))*p_nh%prog(nnow)%rho(jc,jk-1,jb)

            ! perturbation density and virtual potential temperature at main levels for horizontal flux divergence term
            ! (needed in the predictor step only)
            z_rth_pr(1,jc,jk,jb) =  p_nh%prog(nnow)%rho(jc,jk,jb)     - p_nh%metrics%rho_ref_mc(jc,jk,jb)
            z_rth_pr(2,jc,jk,jb) =  p_nh%prog(nnow)%theta_v(jc,jk,jb) - p_nh%metrics%theta_ref_mc(jc,jk,jb)

            ! perturbation virtual potential temperature at main levels (needed in both substeps - see branch below)
            z_theta_v_pr_mc(jc,jk) = z_rth_pr(2,jc,jk,jb)

            ! perturbation virtual potential temperature at interface levels
            z_theta_v_pr_ic(jc,jk) = &
              p_nh%metrics%wgtfac_c(jc,jk,jb)*z_theta_v_pr_mc(jc,jk) +       &
              (1._vp-p_nh%metrics%wgtfac_c(jc,jk,jb))*z_theta_v_pr_mc(jc,jk-1)

            ! virtual potential temperature at interface levels
            p_nh%diag%theta_v_ic(jc,jk,jb) = &
              p_nh%metrics%wgtfac_c(jc,jk,jb)*p_nh%prog(nnow)%theta_v(jc,jk,jb) +       &
              (1._wp-p_nh%metrics%wgtfac_c(jc,jk,jb))*p_nh%prog(nnow)%theta_v(jc,jk-1,jb)

            ! vertical pressure gradient * theta_v
            z_th_ddz_exner_c(jc,jk,jb) = p_nh%metrics%vwind_expl_wgt(jc,jb)* &
              p_nh%diag%theta_v_ic(jc,jk,jb) * (z_exner_pr(jc,jk-1,jb)-      &
              z_exner_pr(jc,jk,jb)) / p_nh%metrics%ddqz_z_half(jc,jk,jb) +   &
              z_theta_v_pr_ic(jc,jk)*p_nh%metrics%d_exner_dz_ref_ic(jc,jk,jb)
          ENDDO
        ENDDO

      ELSE  ! istep = 2 - in this step, an upwind-biased discretization is used for rho_ic and theta_v_ic
            ! in order to reduce the numerical dispersion errors

        z_theta_tavg(i_startidx:i_endidx,1) = wgt_nnow_rth*p_nh%prog(nnow)%theta_v(i_startidx:i_endidx,1,jb) +  &
          wgt_nnew_rth*p_nh%prog(nvar)%theta_v(i_startidx:i_endidx,1,jb)

        z_theta_v_pr_mc(i_startidx:i_endidx,1) = z_theta_tavg(i_startidx:i_endidx,1) - &
          p_nh%metrics%theta_ref_mc(i_startidx:i_endidx,1,jb)

        z_rho_tavg(i_startidx:i_endidx,1) = wgt_nnow_rth*p_nh%prog(nnow)%rho(i_startidx:i_endidx,1,jb) +  &
          wgt_nnew_rth*p_nh%prog(nvar)%rho(i_startidx:i_endidx,1,jb)

        DO jk = 2, nlev
!DIR$ IVDEP
          DO jc = i_startidx, i_endidx
            ! backward trajectory - use w(nnew) in order to be at the same time level as w_concorr
            z_w_backtraj(jc,jk) = - (p_nh%prog(nnew)%w(jc,jk,jb) - p_nh%diag%w_concorr_c(jc,jk,jb)) * &
              dtime*0.5_wp/p_nh%metrics%ddqz_z_half(jc,jk,jb)

            ! temporally averaged density and virtual potential temperature depending on rhotheta_offctr
            ! (see pre-computation above)
            z_rho_tavg(jc,jk) = wgt_nnow_rth*p_nh%prog(nnow)%rho(jc,jk,jb) + &
              wgt_nnew_rth*p_nh%prog(nvar)%rho(jc,jk,jb)
            z_theta_tavg(jc,jk) = wgt_nnow_rth*p_nh%prog(nnow)%theta_v(jc,jk,jb) + &
              wgt_nnew_rth*p_nh%prog(nvar)%theta_v(jc,jk,jb)

            ! density at interface levels for vertical flux divergence computation
            p_nh%diag%rho_ic(jc,jk,jb) = p_nh%metrics%wgtfac_c(jc,jk,jb)*z_rho_tavg(jc,jk) + &
              (1._wp-p_nh%metrics%wgtfac_c(jc,jk,jb))*z_rho_tavg(jc,jk-1) + &
              z_w_backtraj(jc,jk)*(z_rho_tavg(jc,jk-1)-z_rho_tavg(jc,jk))

            ! perturbation virtual potential temperature at main levels
            z_theta_v_pr_mc(jc,jk) = z_theta_tavg(jc,jk) - p_nh%metrics%theta_ref_mc(jc,jk,jb)

            ! perturbation virtual potential temperature at interface levels
            z_theta_v_pr_ic(jc,jk) = &
              p_nh%metrics%wgtfac_c(jc,jk,jb)*z_theta_v_pr_mc(jc,jk) +       &
              (1._vp-p_nh%metrics%wgtfac_c(jc,jk,jb))*z_theta_v_pr_mc(jc,jk-1)

            ! virtual potential temperature at interface levels
            p_nh%diag%theta_v_ic(jc,jk,jb) = p_nh%metrics%wgtfac_c(jc,jk,jb)*z_theta_tavg(jc,jk) + &
              (1._wp-p_nh%metrics%wgtfac_c(jc,jk,jb))*z_theta_tavg(jc,jk-1) +                      &
              z_w_backtraj(jc,jk)*(z_theta_tavg(jc,jk-1)-z_theta_tavg(jc,jk))

            ! vertical pressure gradient * theta_v
            z_th_ddz_exner_c(jc,jk,jb) = p_nh%metrics%vwind_expl_wgt(jc,jb)* &
              p_nh%diag%theta_v_ic(jc,jk,jb) * (z_exner_pr(jc,jk-1,jb)-      &
              z_exner_pr(jc,jk,jb)) / p_nh%metrics%ddqz_z_half(jc,jk,jb) +   &
              z_theta_v_pr_ic(jc,jk)*p_nh%metrics%d_exner_dz_ref_ic(jc,jk,jb)
          ENDDO
        ENDDO

      ENDIF ! istep = 1/2

      ! rho and theta at top level (in case of vertical nesting, upper boundary conditions 
      !                             are set in the vertical solver loop)
      IF (l_open_ubc .AND. .NOT. l_vert_nested) THEN
!DIR$ IVDEP
        DO jc = i_startidx, i_endidx
          p_nh%diag%theta_v_ic(jc,1,jb) = p_nh%metrics%theta_ref_ic(jc,1,jb) + &
            p_nh%metrics%wgtfacq1_c(jc,1,jb)*z_theta_v_pr_mc(jc,1) +           &
            p_nh%metrics%wgtfacq1_c(jc,2,jb)*z_theta_v_pr_mc(jc,2) +           &
            p_nh%metrics%wgtfacq1_c(jc,3,jb)*z_theta_v_pr_mc(jc,3)
          p_nh%diag%rho_ic(jc,1,jb) =  wgt_nnow_rth*(                        &
            p_nh%metrics%wgtfacq1_c(jc,1,jb)*p_nh%prog(nnow)%rho(jc,1,jb) +  &
            p_nh%metrics%wgtfacq1_c(jc,2,jb)*p_nh%prog(nnow)%rho(jc,2,jb) +  &
            p_nh%metrics%wgtfacq1_c(jc,3,jb)*p_nh%prog(nnow)%rho(jc,3,jb))+  &
            wgt_nnew_rth * (                                                 &
            p_nh%metrics%wgtfacq1_c(jc,1,jb)*p_nh%prog(nvar)%rho(jc,1,jb) +  &
            p_nh%metrics%wgtfacq1_c(jc,2,jb)*p_nh%prog(nvar)%rho(jc,2,jb) +  &
            p_nh%metrics%wgtfacq1_c(jc,3,jb)*p_nh%prog(nvar)%rho(jc,3,jb) )
        ENDDO
      ENDIF

      IF (istep == 1) THEN

        ! Perturbation theta at top and surface levels
!DIR$ IVDEP
        DO jc = i_startidx, i_endidx
          z_theta_v_pr_ic(jc,1)      = 0._wp
          z_theta_v_pr_ic(jc,nlevp1) =                                   &
            p_nh%metrics%wgtfacq_c(jc,1,jb)*z_theta_v_pr_mc(jc,nlev) +   &
            p_nh%metrics%wgtfacq_c(jc,2,jb)*z_theta_v_pr_mc(jc,nlev-1) + &
            p_nh%metrics%wgtfacq_c(jc,3,jb)*z_theta_v_pr_mc(jc,nlev-2)
            p_nh%diag%theta_v_ic(jc,nlevp1,jb) =                                  &
              p_nh%metrics%theta_ref_ic(jc,nlevp1,jb) + z_theta_v_pr_ic(jc,nlevp1)
        ENDDO

        IF (igradp_method <= 3) THEN
          DO jk = nflat_gradp(p_patch%id), nlev
!DIR$ IVDEP
            DO jc = i_startidx, i_endidx
              ! Second vertical derivative of perturbation Exner pressure (hydrostatic approximation)
              z_dexner_dz_c(2,jc,jk,jb) = -0.5_vp *                                &
               ((z_theta_v_pr_ic(jc,jk) - z_theta_v_pr_ic(jc,jk+1)) *              &
                p_nh%metrics%d2dexdz2_fac1_mc(jc,jk,jb) + z_theta_v_pr_mc(jc,jk)*  &
                p_nh%metrics%d2dexdz2_fac2_mc(jc,jk,jb))
            ENDDO
          ENDDO
        ENDIF

      ENDIF ! istep == 1

    ENDDO
!$OMP END DO NOWAIT

    IF (istep == 1) THEN
      ! Add computation of z_grad_rth (perturbation density and virtual potential temperature at main levels)
      ! at outer halo points: needed for correct calculation of the upwind gradients for Miura scheme
      rl_start = min_rlcell_int - 2
      rl_end   = min_rlcell_int - 2

      i_startblk = p_patch%cells%start_block(rl_start)
      i_endblk   = p_patch%cells%end_block(rl_end)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,jc) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)

        DO jk = 1, nlev
!DIR$ IVDEP
          DO jc = i_startidx, i_endidx

            z_rth_pr(1,jc,jk,jb) =  p_nh%prog(nnow)%rho(jc,jk,jb)     - p_nh%metrics%rho_ref_mc(jc,jk,jb)
            z_rth_pr(2,jc,jk,jb) =  p_nh%prog(nnow)%theta_v(jc,jk,jb) - p_nh%metrics%theta_ref_mc(jc,jk,jb)

          ENDDO
        ENDDO
      ENDDO
!$OMP END DO NOWAIT
    ENDIF
!$OMP END PARALLEL

    IF (timers_level > 5) THEN
      CALL timer_stop(timer_solve_nh_cellcomp)
      CALL timer_start(timer_solve_nh_vnupd)
    ENDIF

    ! Compute rho and theta at edges for horizontal flux divergence term
    IF (istep == 1) THEN
      IF (iadv_rhotheta == 1) THEN ! Simplified Miura scheme

        ! Compute density and potential temperature at vertices
        CALL cells2verts_scalar(p_nh%prog(nnow)%rho,p_patch, p_int%cells_aw_verts, &
          z_rho_v, opt_rlend=min_rlvert_int-1)
        CALL cells2verts_scalar(p_nh%prog(nnow)%theta_v,p_patch, p_int%cells_aw_verts, &
          z_theta_v_v, opt_rlend=min_rlvert_int-1)

      ELSE IF (iadv_rhotheta == 2) THEN ! Miura second-order upwind scheme

#ifndef __LOOP_EXCHANGE
        ! Compute backward trajectory - code is inlined for cache-based machines (see below)
        CALL btraj(p_patch, p_int, p_nh%prog(nnow)%vn, p_nh%diag%vt, &
                   0.5_wp*dtime, z_cell_idx, z_cell_blk, z_distv_bary, &
                   opt_rlstart=7, opt_rlend=min_rledge_int-1 )
#endif

        ! Compute Green-Gauss gradients for rho and theta
        CALL grad_green_gauss_cell(z_rth_pr, p_patch, p_int, z_grad_rth,    &
                                   opt_rlstart=3, opt_rlend=min_rlcell_int-1)

      ELSE IF (iadv_rhotheta == 3) THEN ! Third-order Miura scheme (does not perform well yet)

        lcompute =.TRUE.
        lcleanup =.FALSE.
        ! First call: compute backward trajectory with wind at time level nnow
        CALL upwind_hflux_miura3(p_patch, p_nh%prog(nnow)%rho, p_nh%prog(nnow)%vn, &
                                p_nh%prog(nnow)%vn, p_nh%diag%vt, dtime, p_int,    &
                                lcompute, lcleanup, 0, z_rho_e,                    &
                                opt_rlstart=7, opt_lout_edge=.TRUE. )

        ! Second call: compute only reconstructed value for flux divergence
        lcompute =.FALSE.
        lcleanup =.TRUE.
        CALL upwind_hflux_miura3(p_patch, p_nh%prog(nnow)%theta_v, p_nh%prog(nnow)%vn, &
                                p_nh%prog(nnow)%vn, p_nh%diag%vt, dtime, p_int,        &
                                lcompute, lcleanup, 0, z_theta_v_e,                    &
                                opt_rlstart=7, opt_lout_edge=.TRUE. )

      ENDIF
    ENDIF ! istep = 1

!$OMP PARALLEL PRIVATE (rl_start,rl_end,i_startblk,i_endblk)

    IF (istep == 1) THEN 
      ! Compute 'edge values' of density and virtual potential temperature for horizontal
      ! flux divergence term; this is included in upwind_hflux_miura3 for option 3
      IF (iadv_rhotheta <= 2) THEN

        ! Initialize halo edges with zero in order to avoid access of uninitialized array elements
        i_startblk = p_patch%edges%start_block(min_rledge_int-2)
        i_endblk   = p_patch%edges%end_block  (min_rledge_int-3)

!$OMP WORKSHARE
        z_rho_e    (:,:,i_startblk:i_endblk) = 0._wp
        z_theta_v_e(:,:,i_startblk:i_endblk) = 0._wp
!$OMP END WORKSHARE

        rl_start = 7
        rl_end   = min_rledge_int-1

        i_startblk = p_patch%edges%start_block(rl_start)
        i_endblk   = p_patch%edges%end_block  (rl_end)

        ! initialize also nest boundary points with zero
        IF (p_patch%id > 1 .OR. l_limited_area) THEN
!$OMP WORKSHARE
          z_rho_e    (:,:,1:i_startblk) = 0._wp
          z_theta_v_e(:,:,1:i_startblk) = 0._wp
!$OMP END WORKSHARE
        ENDIF

!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx,ilc0,ibc0,lvn_pos,z_ntdistv_bary,distv_bary) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = i_startblk, i_endblk

          CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx, rl_start, rl_end)

          IF (iadv_rhotheta == 2) THEN
            ! Operations from upwind_hflux_miura are inlined in order to process both
            ! fields in one step

#ifdef __LOOP_EXCHANGE
            ! For cache-based machines, also the back-trajectory computation is inlined to improve efficiency
            DO je = i_startidx, i_endidx
!DIR$ IVDEP, PREFERVECTOR
              DO jk = 1, nlev

                lvn_pos = p_nh%prog(nnow)%vn(je,jk,jb) >= 0._wp

                ! line and block indices of upwind neighbor cell
                ilc0 = MERGE(p_patch%edges%cell_idx(je,jb,1),p_patch%edges%cell_idx(je,jb,2),lvn_pos)
                ibc0 = MERGE(p_patch%edges%cell_blk(je,jb,1),p_patch%edges%cell_blk(je,jb,2),lvn_pos)

                ! distances from upwind mass point to the end point of the backward trajectory
                ! in edge-normal and tangential directions
                z_ntdistv_bary(1) =  - ( p_nh%prog(nnow)%vn(je,jk,jb) * dthalf +    &
                  MERGE(p_int%pos_on_tplane_e(je,jb,1,1), p_int%pos_on_tplane_e(je,jb,2,1),lvn_pos))

                z_ntdistv_bary(2) =  - ( p_nh%diag%vt(je,jk,jb) * dthalf +    &
                  MERGE(p_int%pos_on_tplane_e(je,jb,1,2), p_int%pos_on_tplane_e(je,jb,2,2),lvn_pos))

                ! rotate distance vectors into local lat-lon coordinates:
                !
                ! component in longitudinal direction
                distv_bary(1) =                                                                  &
                   z_ntdistv_bary(1)*MERGE(p_patch%edges%primal_normal_cell(je,jb,1)%v1,         &
                                           p_patch%edges%primal_normal_cell(je,jb,2)%v1,lvn_pos) &
                 + z_ntdistv_bary(2)*MERGE(p_patch%edges%dual_normal_cell(je,jb,1)%v1,           &
                                           p_patch%edges%dual_normal_cell(je,jb,2)%v1,lvn_pos)

                ! component in latitudinal direction
                distv_bary(2) =                                                                  &
                   z_ntdistv_bary(1)*MERGE(p_patch%edges%primal_normal_cell(je,jb,1)%v2,         &
                                           p_patch%edges%primal_normal_cell(je,jb,2)%v2,lvn_pos) &
                 + z_ntdistv_bary(2)*MERGE(p_patch%edges%dual_normal_cell(je,jb,1)%v2,           &
                                           p_patch%edges%dual_normal_cell(je,jb,2)%v2,lvn_pos)


                ! Calculate "edge values" of rho and theta_v
                ! Note: z_rth_pr contains the perturbation values of rho and theta_v,
                ! and the corresponding gradients are stored in z_grad_rth.
                z_rho_e(je,jk,jb) = p_nh%metrics%rho_ref_me(je,jk,jb)     &
                  +                            z_rth_pr(1,ilc0,jk,ibc0)   &
                  + distv_bary(1) * z_grad_rth(1,ilc0,jk,ibc0) &
                  + distv_bary(2) * z_grad_rth(2,ilc0,jk,ibc0)

                z_theta_v_e(je,jk,jb) = p_nh%metrics%theta_ref_me(je,jk,jb) &
                  +                            z_rth_pr(2,ilc0,jk,ibc0)     &
                  + distv_bary(1) * z_grad_rth(3,ilc0,jk,ibc0)   &
                  + distv_bary(2) * z_grad_rth(4,ilc0,jk,ibc0)

#else
            DO jk = 1, nlev
              DO je = i_startidx, i_endidx

                ilc0 = z_cell_idx(je,jk,jb)
                ibc0 = z_cell_blk(je,jk,jb)  

                ! Calculate "edge values" of rho and theta_v
                ! Note: z_rth_pr contains the perturbation values of rho and theta_v,
                ! and the corresponding gradients are stored in z_grad_rth.
                z_rho_e(je,jk,jb) = p_nh%metrics%rho_ref_me(je,jk,jb)     &
                  +                            z_rth_pr(1,ilc0,jk,ibc0)   &
                  + z_distv_bary(je,jk,jb,1) * z_grad_rth(1,ilc0,jk,ibc0) &
                  + z_distv_bary(je,jk,jb,2) * z_grad_rth(2,ilc0,jk,ibc0)

                z_theta_v_e(je,jk,jb) = p_nh%metrics%theta_ref_me(je,jk,jb) &
                  +                            z_rth_pr(2,ilc0,jk,ibc0)     &
                  + z_distv_bary(je,jk,jb,1) * z_grad_rth(3,ilc0,jk,ibc0)   &
                  + z_distv_bary(je,jk,jb,2) * z_grad_rth(4,ilc0,jk,ibc0)
#endif

              ENDDO ! loop over edges
            ENDDO   ! loop over vertical levels

          ELSE ! iadv_rhotheta = 1

#ifdef __LOOP_EXCHANGE
            DO je = i_startidx, i_endidx
!DIR$ IVDEP
              DO jk = 1, nlev
#else
            DO jk = 1, nlev
              DO je = i_startidx, i_endidx
#endif

                ! Compute upwind-biased values for rho and theta starting from centered differences
                ! Note: the length of the backward trajectory should be 0.5*dtime*(vn,vt) in order to arrive
                ! at a second-order accurate FV discretization, but twice the length is needed for numerical
                ! stability
                z_rho_e(je,jk,jb) =                                                                          &
                  p_int%c_lin_e(je,1,jb)*p_nh%prog(nnow)%rho(icidx(je,jb,1),jk,icblk(je,jb,1)) +             &
                  p_int%c_lin_e(je,2,jb)*p_nh%prog(nnow)%rho(icidx(je,jb,2),jk,icblk(je,jb,2)) -             &
                  dtime * (p_nh%prog(nnow)%vn(je,jk,jb)*p_patch%edges%inv_dual_edge_length(je,jb)*           &
                 (p_nh%prog(nnow)%rho(icidx(je,jb,2),jk,icblk(je,jb,2)) -                                    &
                  p_nh%prog(nnow)%rho(icidx(je,jb,1),jk,icblk(je,jb,1)) ) + p_nh%diag%vt(je,jk,jb) *         &
                  p_patch%edges%inv_primal_edge_length(je,jb) * p_patch%edges%system_orientation(je,jb) *    &
                 (z_rho_v(ividx(je,jb,2),jk,ivblk(je,jb,2)) - z_rho_v(ividx(je,jb,1),jk,ivblk(je,jb,1)) ) )

                z_theta_v_e(je,jk,jb) =                                                                          &
                  p_int%c_lin_e(je,1,jb)*p_nh%prog(nnow)%theta_v(icidx(je,jb,1),jk,icblk(je,jb,1)) +             &
                  p_int%c_lin_e(je,2,jb)*p_nh%prog(nnow)%theta_v(icidx(je,jb,2),jk,icblk(je,jb,2)) -             &
                  dtime * (p_nh%prog(nnow)%vn(je,jk,jb)*p_patch%edges%inv_dual_edge_length(je,jb)*               &
                 (p_nh%prog(nnow)%theta_v(icidx(je,jb,2),jk,icblk(je,jb,2)) -                                    &
                  p_nh%prog(nnow)%theta_v(icidx(je,jb,1),jk,icblk(je,jb,1)) ) + p_nh%diag%vt(je,jk,jb) *         &
                  p_patch%edges%inv_primal_edge_length(je,jb) * p_patch%edges%system_orientation(je,jb) *        &
                 (z_theta_v_v(ividx(je,jb,2),jk,ivblk(je,jb,2)) - z_theta_v_v(ividx(je,jb,1),jk,ivblk(je,jb,1)) ))

              ENDDO ! loop over edges
            ENDDO   ! loop over vertical levels
          ENDIF

        ENDDO
!$OMP END DO

      ENDIF

    ELSE IF (istep == 2 .AND. lhdiff_rcf .AND. divdamp_type >= 3) THEN ! apply div damping on 3D divergence

      ! add dw/dz contribution to divergence damping term

      rl_start = 7
      rl_end   = min_rledge_int-2

      i_startblk = p_patch%edges%start_block(rl_start)
      i_endblk   = p_patch%edges%end_block  (rl_end)

!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, rl_start, rl_end)

#ifdef __LOOP_EXCHANGE
        DO je = i_startidx, i_endidx
!DIR$ IVDEP, PREFERVECTOR
          DO jk = 1, nlev
            z_graddiv_vn(jk,je,jb) = z_graddiv_vn(jk,je,jb) + scal_div3d(jk) * p_patch%edges%inv_dual_edge_length(je,jb)* &
             ( z_dwdz_dd(icidx(je,jb,2),jk,icblk(je,jb,2)) - z_dwdz_dd(icidx(je,jb,1),jk,icblk(je,jb,1)) )
#else
        DO jk = 1, nlev
          DO je = i_startidx, i_endidx
            z_graddiv_vn(je,jk,jb) = z_graddiv_vn(je,jk,jb) + scal_div3d(jk) * p_patch%edges%inv_dual_edge_length(je,jb)* &
             ( z_dwdz_dd(icidx(je,jb,2),jk,icblk(je,jb,2)) - z_dwdz_dd(icidx(je,jb,1),jk,icblk(je,jb,1)) )
#endif
          ENDDO
        ENDDO
      ENDDO
!$OMP END DO

    ENDIF ! istep = 1/2

    ! Remaining computations at edge points

    rl_start = grf_bdywidth_e + 1   ! boundary update follows below
    rl_end   = min_rledge_int

    i_startblk = p_patch%edges%start_block(rl_start)
    i_endblk   = p_patch%edges%end_block(rl_end)

    IF (istep == 1) THEN
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,je,z_theta1,z_theta2,ikp1,ikp2) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, rl_start, rl_end)

        ! Store values at nest interface levels
        IF (idyn_timestep == 1 .AND. l_child_vertnest) THEN
!DIR$ IVDEP
          DO je = i_startidx, i_endidx
            p_nh%diag%dvn_ie_int(je,jb) = p_nh%diag%vn_ie(je,nshift,jb) - &
                                          p_nh%diag%vn_ie(je,nshift+1,jb)
          ENDDO
        ENDIF

#ifdef __LOOP_EXCHANGE
        DO je = i_startidx, i_endidx
          DO jk = 1, nflatlev(p_patch%id)-1
#else
        DO jk = 1, nflatlev(p_patch%id)-1
          DO je = i_startidx, i_endidx
#endif
            ! horizontal gradient of Exner pressure where coordinate surfaces are flat
            z_gradh_exner(je,jk,jb) = p_patch%edges%inv_dual_edge_length(je,jb)* &
             (z_exner_ex_pr(icidx(je,jb,2),jk,icblk(je,jb,2)) -                  &
              z_exner_ex_pr(icidx(je,jb,1),jk,icblk(je,jb,1)) )
          ENDDO
        ENDDO

        IF (igradp_method <= 3) THEN
#ifdef __LOOP_EXCHANGE
          DO je = i_startidx, i_endidx
!DIR$ IVDEP
            DO jk = nflatlev(p_patch%id), nflat_gradp(p_patch%id)
#else
          DO jk = nflatlev(p_patch%id), nflat_gradp(p_patch%id)
            DO je = i_startidx, i_endidx
#endif
              ! horizontal gradient of Exner pressure, including metric correction
              z_gradh_exner(je,jk,jb) = p_patch%edges%inv_dual_edge_length(je,jb)*         &
               (z_exner_ex_pr(icidx(je,jb,2),jk,icblk(je,jb,2)) -                          &
                z_exner_ex_pr(icidx(je,jb,1),jk,icblk(je,jb,1)) ) -                        &
                p_nh%metrics%ddxn_z_full(je,jk,jb) *                                       &
               (p_int%c_lin_e(je,1,jb)*z_dexner_dz_c(1,icidx(je,jb,1),jk,icblk(je,jb,1)) + &
                p_int%c_lin_e(je,2,jb)*z_dexner_dz_c(1,icidx(je,jb,2),jk,icblk(je,jb,2)))

            ENDDO
          ENDDO

#ifdef __LOOP_EXCHANGE
          DO je = i_startidx, i_endidx
!DIR$ IVDEP, PREFERVECTOR
            DO jk = nflat_gradp(p_patch%id)+1, nlev
#else
          DO jk = nflat_gradp(p_patch%id)+1, nlev
            DO je = i_startidx, i_endidx
#endif
              ! horizontal gradient of Exner pressure, Taylor-expansion-based reconstruction
              z_gradh_exner(je,jk,jb) = p_patch%edges%inv_dual_edge_length(je,jb)*         &
               (z_exner_ex_pr(icidx(je,jb,2),ikidx(2,je,jk,jb),icblk(je,jb,2)) +           &
                p_nh%metrics%zdiff_gradp(2,je,jk,jb)*                                      &
               (z_dexner_dz_c(1,icidx(je,jb,2),ikidx(2,je,jk,jb),icblk(je,jb,2)) +         &
                p_nh%metrics%zdiff_gradp(2,je,jk,jb)*                                      &
                z_dexner_dz_c(2,icidx(je,jb,2),ikidx(2,je,jk,jb),icblk(je,jb,2))) -        &
               (z_exner_ex_pr(icidx(je,jb,1),ikidx(1,je,jk,jb),icblk(je,jb,1)) +           &
                p_nh%metrics%zdiff_gradp(1,je,jk,jb)*                                      &
               (z_dexner_dz_c(1,icidx(je,jb,1),ikidx(1,je,jk,jb),icblk(je,jb,1)) +         &
                p_nh%metrics%zdiff_gradp(1,je,jk,jb)*                                      &
                z_dexner_dz_c(2,icidx(je,jb,1),ikidx(1,je,jk,jb),icblk(je,jb,1)))))

            ENDDO
          ENDDO
        ELSE IF (igradp_method == 4 .OR. igradp_method == 5) THEN
#ifdef __LOOP_EXCHANGE
          DO je = i_startidx, i_endidx
            DO jk = nflatlev(p_patch%id), nlev
#else
          DO jk = nflatlev(p_patch%id), nlev
            DO je = i_startidx, i_endidx
#endif
              ! horizontal gradient of Exner pressure, cubic/quadratic interpolation
              z_gradh_exner(je,jk,jb) = p_patch%edges%inv_dual_edge_length(je,jb)*   &
               (z_exner_ex_pr(icidx(je,jb,2),ikidx(2,je,jk,jb)-1,icblk(je,jb,2)) *   &
                p_nh%metrics%coeff_gradp(5,je,jk,jb) +                               &
                z_exner_ex_pr(icidx(je,jb,2),ikidx(2,je,jk,jb)  ,icblk(je,jb,2)) *   &
                p_nh%metrics%coeff_gradp(6,je,jk,jb) +                               &
                z_exner_ex_pr(icidx(je,jb,2),ikidx(2,je,jk,jb)+1,icblk(je,jb,2)) *   &
                p_nh%metrics%coeff_gradp(7,je,jk,jb) +                               &
                z_exner_ex_pr(icidx(je,jb,2),ikidx(2,je,jk,jb)+2,icblk(je,jb,2)) *   &
                p_nh%metrics%coeff_gradp(8,je,jk,jb) -                               &
               (z_exner_ex_pr(icidx(je,jb,1),ikidx(1,je,jk,jb)-1,icblk(je,jb,1)) *   &
                p_nh%metrics%coeff_gradp(1,je,jk,jb) +                               &
                z_exner_ex_pr(icidx(je,jb,1),ikidx(1,je,jk,jb)  ,icblk(je,jb,1)) *   &
                p_nh%metrics%coeff_gradp(2,je,jk,jb) +                               &
                z_exner_ex_pr(icidx(je,jb,1),ikidx(1,je,jk,jb)+1,icblk(je,jb,1)) *   &
                p_nh%metrics%coeff_gradp(3,je,jk,jb) +                               &
                z_exner_ex_pr(icidx(je,jb,1),ikidx(1,je,jk,jb)+2,icblk(je,jb,1)) *   &
                p_nh%metrics%coeff_gradp(4,je,jk,jb)) )

            ENDDO
          ENDDO
        ENDIF

        ! compute hydrostatically approximated correction term that replaces downward extrapolation
        IF (igradp_method == 3) THEN

          DO je = i_startidx, i_endidx

            z_theta1 = &
              p_nh%prog(nnow)%theta_v(icidx(je,jb,1),ikidx(1,je,nlev,jb),icblk(je,jb,1)) +  &
              p_nh%metrics%zdiff_gradp(1,je,nlev,jb)*                                       &
             (p_nh%diag%theta_v_ic(icidx(je,jb,1),ikidx(1,je,nlev,jb),  icblk(je,jb,1)) -   &
              p_nh%diag%theta_v_ic(icidx(je,jb,1),ikidx(1,je,nlev,jb)+1,icblk(je,jb,1))) *  &
              p_nh%metrics%inv_ddqz_z_full(icidx(je,jb,1),ikidx(1,je,nlev,jb),icblk(je,jb,1))

            z_theta2 = &
              p_nh%prog(nnow)%theta_v(icidx(je,jb,2),ikidx(2,je,nlev,jb),icblk(je,jb,2)) +  &
              p_nh%metrics%zdiff_gradp(2,je,nlev,jb)*                                       &
             (p_nh%diag%theta_v_ic(icidx(je,jb,2),ikidx(2,je,nlev,jb),  icblk(je,jb,2)) -   &
              p_nh%diag%theta_v_ic(icidx(je,jb,2),ikidx(2,je,nlev,jb)+1,icblk(je,jb,2))) *  &
              p_nh%metrics%inv_ddqz_z_full(icidx(je,jb,2),ikidx(2,je,nlev,jb),icblk(je,jb,2))

            z_hydro_corr(je,jb) = grav_o_cpd*p_patch%edges%inv_dual_edge_length(je,jb)*    &
              (z_theta2-z_theta1)*4._wp/(z_theta1+z_theta2)**2

          ENDDO
        ELSE IF (igradp_method == 5) THEN

          DO je = i_startidx, i_endidx

            ikp1 = MIN(nlev,ikidx(1,je,nlev,jb)+2)
            ikp2 = MIN(nlev,ikidx(2,je,nlev,jb)+2)

            z_theta1 =                                                                       &
              p_nh%prog(nnow)%theta_v(icidx(je,jb,1),ikidx(1,je,nlev,jb)-1,icblk(je,jb,1)) * &
              p_nh%metrics%coeff_gradp(1,je,nlev,jb) +                                         &
              p_nh%prog(nnow)%theta_v(icidx(je,jb,1),ikidx(1,je,nlev,jb)  ,icblk(je,jb,1)) * &
              p_nh%metrics%coeff_gradp(2,je,nlev,jb) +                                         &
              p_nh%prog(nnow)%theta_v(icidx(je,jb,1),ikidx(1,je,nlev,jb)+1,icblk(je,jb,1)) * &
              p_nh%metrics%coeff_gradp(3,je,nlev,jb) +                                         &
              p_nh%prog(nnow)%theta_v(icidx(je,jb,1),ikp1                 ,icblk(je,jb,1)) * &
              p_nh%metrics%coeff_gradp(4,je,nlev,jb) 

            z_theta2 =                                                                       &
              p_nh%prog(nnow)%theta_v(icidx(je,jb,2),ikidx(2,je,nlev,jb)-1,icblk(je,jb,2)) * &
              p_nh%metrics%coeff_gradp(5,je,nlev,jb) +                                         &
              p_nh%prog(nnow)%theta_v(icidx(je,jb,2),ikidx(2,je,nlev,jb)  ,icblk(je,jb,2)) * &
              p_nh%metrics%coeff_gradp(6,je,nlev,jb) +                                         &
              p_nh%prog(nnow)%theta_v(icidx(je,jb,2),ikidx(2,je,nlev,jb)+1,icblk(je,jb,2)) * &
              p_nh%metrics%coeff_gradp(7,je,nlev,jb) +                                         &
              p_nh%prog(nnow)%theta_v(icidx(je,jb,2),ikp2                 ,icblk(je,jb,2)) * &
              p_nh%metrics%coeff_gradp(8,je,nlev,jb) 

            z_hydro_corr(je,jb) = grav_o_cpd*p_patch%edges%inv_dual_edge_length(je,jb)*    &
              (z_theta2-z_theta1)*4._wp/(z_theta1+z_theta2)**2

          ENDDO
        ENDIF

      ENDDO
!$OMP END DO
    ENDIF ! istep = 1

    IF (istep == 1 .AND. (igradp_method == 3 .OR. igradp_method == 5)) THEN

!$OMP DO PRIVATE(jb,je,ie,nlen_gradp,ishift) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = 1, nblks_gradp
        IF (jb == nblks_gradp) THEN
          nlen_gradp = npromz_gradp
        ELSE
          nlen_gradp = nproma_gradp
        ENDIF
        ishift = (jb-1)*nproma_gradp
!CDIR NODEP,VOVERTAKE,VOB
        DO je = 1, nlen_gradp
          ie = ishift+je

          z_gradh_exner(ipeidx(ie),iplev(ie),ipeblk(ie))  =              &
            z_gradh_exner(ipeidx(ie),iplev(ie),ipeblk(ie)) +             &
            p_nh%metrics%pg_exdist(ie)*z_hydro_corr(ipeidx(ie),ipeblk(ie))

        ENDDO
      ENDDO
!$OMP END DO
    ENDIF

    ! Update horizontal velocity field: advection (including Coriolis force) and pressure-gradient term
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,je,z_graddiv2_vn) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)

      IF ((itime_scheme >= 4) .AND. istep == 2) THEN ! use temporally averaged velocity advection terms
        DO jk = 1, nlev
!DIR$ IVDEP
          DO je = i_startidx, i_endidx
            p_nh%prog(nnew)%vn(je,jk,jb) = p_nh%prog(nnow)%vn(je,jk,jb)+ dtime                  &
            & *(wgt_nnow_vel*p_nh%diag%ddt_vn_adv(je,jk,jb,ntl1)                                &
            & + wgt_nnew_vel*p_nh%diag%ddt_vn_adv(je,jk,jb,ntl2)+p_nh%diag%ddt_vn_phy(je,jk,jb) &
            & -cpd*z_theta_v_e(je,jk,jb)*z_gradh_exner(je,jk,jb))
          ENDDO
        ENDDO
      ELSE
        DO jk = 1, nlev
!DIR$ IVDEP
          DO je = i_startidx, i_endidx
            p_nh%prog(nnew)%vn(je,jk,jb) = p_nh%prog(nnow)%vn(je,jk,jb)+ dtime     &
            & *(p_nh%diag%ddt_vn_adv(je,jk,jb,ntl1)+p_nh%diag%ddt_vn_phy(je,jk,jb) &
            & -cpd*z_theta_v_e(je,jk,jb)*z_gradh_exner(je,jk,jb))
          ENDDO
        ENDDO
      ENDIF

      IF (lhdiff_rcf .AND. istep == 2 .AND. ANY( (/24,4/) == divdamp_order)) THEN ! fourth-order divergence damping

        ! Compute gradient of divergence of gradient of divergence for fourth-order divergence damping
#ifdef __LOOP_EXCHANGE
        DO je = i_startidx, i_endidx
!DIR$ IVDEP
          DO jk = 1, nlev
            z_graddiv2_vn(je,jk) = p_int%geofac_grdiv(je,1,jb)*z_graddiv_vn(jk,je,jb)      &
              + p_int%geofac_grdiv(je,2,jb)*z_graddiv_vn(jk,iqidx(je,jb,1),iqblk(je,jb,1)) &
              + p_int%geofac_grdiv(je,3,jb)*z_graddiv_vn(jk,iqidx(je,jb,2),iqblk(je,jb,2)) &
              + p_int%geofac_grdiv(je,4,jb)*z_graddiv_vn(jk,iqidx(je,jb,3),iqblk(je,jb,3)) &
              + p_int%geofac_grdiv(je,5,jb)*z_graddiv_vn(jk,iqidx(je,jb,4),iqblk(je,jb,4))
#else
        DO jk = 1, nlev
          DO je = i_startidx, i_endidx
            z_graddiv2_vn(je,jk) = p_int%geofac_grdiv(je,1,jb)*z_graddiv_vn(je,jk,jb)      &
              + p_int%geofac_grdiv(je,2,jb)*z_graddiv_vn(iqidx(je,jb,1),jk,iqblk(je,jb,1)) &
              + p_int%geofac_grdiv(je,3,jb)*z_graddiv_vn(iqidx(je,jb,2),jk,iqblk(je,jb,2)) &
              + p_int%geofac_grdiv(je,4,jb)*z_graddiv_vn(iqidx(je,jb,3),jk,iqblk(je,jb,3)) &
              + p_int%geofac_grdiv(je,5,jb)*z_graddiv_vn(iqidx(je,jb,4),jk,iqblk(je,jb,4))
#endif

           ENDDO
        ENDDO
      ENDIF

      IF (lhdiff_rcf .AND. istep == 2) THEN
        ! apply divergence damping if diffusion is not called every sound-wave time step
        IF (divdamp_order == 2 .OR. (divdamp_order == 24 .AND. scal_divdamp_o2 > 1.e-6_wp) ) THEN ! second-order divergence damping
          DO jk = 1, nlev
!DIR$ IVDEP
            DO je = i_startidx, i_endidx
              p_nh%prog(nnew)%vn(je,jk,jb) = p_nh%prog(nnew)%vn(je,jk,jb) + scal_divdamp_o2* &
#ifdef __LOOP_EXCHANGE
                z_graddiv_vn(jk,je,jb)
#else
                z_graddiv_vn(je,jk,jb)
#endif
            ENDDO
          ENDDO
        ENDIF
        IF (divdamp_order == 4 .OR. (divdamp_order == 24 .AND. divdamp_fac_o2 <= 4._wp*divdamp_fac) ) THEN
          IF (l_limited_area .OR. p_patch%id > 1) THEN 
            ! fourth-order divergence damping with reduced damping coefficient along nest boundary
            ! (scal_divdamp is negative whereas bdy_divdamp is positive; decreasing the divergence
            ! damping along nest boundaries is beneficial because this reduces the interference
            ! with the increased diffusion applied in nh_diffusion)
            DO jk = 1, nlev
!DIR$ IVDEP
              DO je = i_startidx, i_endidx
                p_nh%prog(nnew)%vn(je,jk,jb) = p_nh%prog(nnew)%vn(je,jk,jb)                         &
                  + (scal_divdamp(jk)+bdy_divdamp(jk)*p_int%nudgecoeff_e(je,jb))*z_graddiv2_vn(je,jk)
              ENDDO
            ENDDO
          ELSE ! fourth-order divergence damping
            DO jk = 1, nlev
!DIR$ IVDEP
              DO je = i_startidx, i_endidx
                p_nh%prog(nnew)%vn(je,jk,jb) = p_nh%prog(nnew)%vn(je,jk,jb)  &
                  + scal_divdamp(jk)*z_graddiv2_vn(je,jk)
              ENDDO
            ENDDO
          ENDIF
        ENDIF
      ENDIF

      IF (is_iau_active) THEN ! add analysis increment from data assimilation
        DO jk = 1, nlev
!DIR$ IVDEP
          DO je = i_startidx, i_endidx
            p_nh%prog(nnew)%vn(je,jk,jb) = p_nh%prog(nnew)%vn(je,jk,jb) +  &
              iau_wgt_dyn*p_nh%diag%vn_incr(je,jk,jb)
          ENDDO
        ENDDO
      ENDIF

      ! Classic Rayleigh damping mechanism for vn (requires reference state !!)
      !
      IF ( rayleigh_type == RAYLEIGH_CLASSIC ) THEN
        DO jk = 1, nrdmax(p_patch%id)
!DIR$ IVDEP
          DO je = i_startidx, i_endidx
            p_nh%prog(nnew)%vn(je,jk,jb) = p_nh%prog(nnew)%vn(je,jk,jb)       &
              &                          - dtime*p_nh%metrics%rayleigh_vn(jk) &
              &                          * (p_nh%prog(nnew)%vn(je,jk,jb)      &
              &                          - p_nh%ref%vn_ref(je,jk,jb))
          ENDDO
        ENDDO
      ENDIF
    ENDDO
!$OMP END DO

    IF (istep == 2 .AND. l_bdy_nudge) THEN ! apply boundary nudging if requested
!$OMP DO PRIVATE(jb,jk,je,ic) ICON_OMP_DEFAULT_SCHEDULE
      DO ic = 1, p_nh%metrics%nudge_e_dim
        je = p_nh%metrics%nudge_e_idx(ic)
        jb = p_nh%metrics%nudge_e_blk(ic)
!DIR$ IVDEP
        DO jk = 1, nlev
          p_nh%prog(nnew)%vn(je,jk,jb) = p_nh%prog(nnew)%vn(je,jk,jb)  &
            + p_int%nudgecoeff_e(je,jb)*p_nh%diag%grf_tend_vn(je,jk,jb)
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF


    ! Boundary update of horizontal velocity
    IF (istep == 1 .AND. (l_limited_area .OR. p_patch%id > 1)) THEN
      rl_start = 1
      rl_end   = grf_bdywidth_e

      i_startblk = p_patch%edges%start_block(rl_start)
      i_endblk   = p_patch%edges%end_block(rl_end)


!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,je) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, rl_start, rl_end)

        DO jk = 1, nlev
!DIR$ IVDEP
          DO je = i_startidx, i_endidx
            p_nh%prog(nnew)%vn(je,jk,jb) = p_nh%prog(nnow)%vn(je,jk,jb) + &
              dtime*p_nh%diag%grf_tend_vn(je,jk,jb)
          ENDDO
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF

    ! Preparations for nest boundary interpolation of mass fluxes from parent domain
    IF (p_patch%id > 1 .AND. grf_intmethod_e >= 5 .AND. idiv_method == 1 .AND. jstep == 0 .AND. istep == 1) THEN
!$OMP DO PRIVATE(ic,je,jb,jk) ICON_OMP_DEFAULT_SCHEDULE
      DO ic = 1, p_nh%metrics%bdy_mflx_e_dim
        je = p_nh%metrics%bdy_mflx_e_idx(ic)
        jb = p_nh%metrics%bdy_mflx_e_blk(ic)
!DIR$ IVDEP
        DO jk = 1, nlev
          p_nh%diag%grf_bdy_mflx(jk,ic,2) = p_nh%diag%grf_tend_mflx(je,jk,jb)
          p_nh%diag%grf_bdy_mflx(jk,ic,1) = prep_adv%mass_flx_me(je,jk,jb) - dt_shift*p_nh%diag%grf_bdy_mflx(jk,ic,2)
        ENDDO

      ENDDO
!$OMP END DO
    ENDIF

!$OMP END PARALLEL

    !-------------------------
    ! communication phase
    IF (timers_level > 5) THEN
      CALL timer_stop(timer_solve_nh_vnupd)
      CALL timer_start(timer_solve_nh_exch)
    ENDIF

    IF (use_icon_comm) THEN 
      IF (istep == 1) THEN
        CALL icon_comm_sync(p_nh%prog(nnew)%vn, z_rho_e, p_patch%sync_edges_not_owned, &
          & name="solve_step1_vn")
      ELSE
        CALL icon_comm_sync(p_nh%prog(nnew)%vn, p_patch%sync_edges_not_owned, &
            & name="solve_step2_vn")
      ENDIF
    ELSE IF (itype_comm == 1) THEN
      IF (istep == 1) THEN
        CALL sync_patch_array_mult(SYNC_E,p_patch,2,p_nh%prog(nnew)%vn,z_rho_e)
      ELSE
        CALL sync_patch_array(SYNC_E,p_patch,p_nh%prog(nnew)%vn)
      ENDIF
    ENDIF

    IF (idiv_method == 2 .AND. istep == 1) CALL sync_patch_array(SYNC_E,p_patch,z_theta_v_e)

    IF (timers_level > 5) THEN
      CALL timer_stop(timer_solve_nh_exch)
      CALL timer_start(timer_solve_nh_edgecomp)
    ENDIF
    ! end communication phase
    !-------------------------

!$OMP PARALLEL PRIVATE (rl_start,rl_end,i_startblk,i_endblk)

    rl_start = 5
    rl_end   = min_rledge_int - 2

    i_startblk = p_patch%edges%start_block(rl_start)
    i_endblk   = p_patch%edges%end_block(rl_end)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,je,z_vn_avg) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)

      IF (istep == 1) THEN

#ifdef __LOOP_EXCHANGE
        DO je = i_startidx, i_endidx
!DIR$ IVDEP
          DO jk = 1, nlev
#else
        DO jk = 1, nlev
          DO je = i_startidx, i_endidx
#endif
            ! Average normal wind components in order to get nearly second-order accurate divergence
            z_vn_avg(je,jk) = p_int%e_flx_avg(je,1,jb)*p_nh%prog(nnew)%vn(je,jk,jb)           &
              + p_int%e_flx_avg(je,2,jb)*p_nh%prog(nnew)%vn(iqidx(je,jb,1),jk,iqblk(je,jb,1)) &
              + p_int%e_flx_avg(je,3,jb)*p_nh%prog(nnew)%vn(iqidx(je,jb,2),jk,iqblk(je,jb,2)) &
              + p_int%e_flx_avg(je,4,jb)*p_nh%prog(nnew)%vn(iqidx(je,jb,3),jk,iqblk(je,jb,3)) &
              + p_int%e_flx_avg(je,5,jb)*p_nh%prog(nnew)%vn(iqidx(je,jb,4),jk,iqblk(je,jb,4))

            ! Compute gradient of divergence of vn for divergence damping
#ifdef __LOOP_EXCHANGE
            z_graddiv_vn(jk,je,jb) = p_int%geofac_grdiv(je,1,jb)*p_nh%prog(nnew)%vn(je,jk,jb)    &
#else
            z_graddiv_vn(je,jk,jb) = p_int%geofac_grdiv(je,1,jb)*p_nh%prog(nnew)%vn(je,jk,jb)    &
#endif
              + p_int%geofac_grdiv(je,2,jb)*p_nh%prog(nnew)%vn(iqidx(je,jb,1),jk,iqblk(je,jb,1)) &
              + p_int%geofac_grdiv(je,3,jb)*p_nh%prog(nnew)%vn(iqidx(je,jb,2),jk,iqblk(je,jb,2)) &
              + p_int%geofac_grdiv(je,4,jb)*p_nh%prog(nnew)%vn(iqidx(je,jb,3),jk,iqblk(je,jb,3)) &
              + p_int%geofac_grdiv(je,5,jb)*p_nh%prog(nnew)%vn(iqidx(je,jb,4),jk,iqblk(je,jb,4))

            ! RBF reconstruction of tangential wind component
            p_nh%diag%vt(je,jk,jb) = p_int%rbf_vec_coeff_e(1,je,jb)  &
              * p_nh%prog(nnew)%vn(iqidx(je,jb,1),jk,iqblk(je,jb,1)) &
              + p_int%rbf_vec_coeff_e(2,je,jb)                       &
              * p_nh%prog(nnew)%vn(iqidx(je,jb,2),jk,iqblk(je,jb,2)) &
              + p_int%rbf_vec_coeff_e(3,je,jb)                       &
              * p_nh%prog(nnew)%vn(iqidx(je,jb,3),jk,iqblk(je,jb,3)) &
              + p_int%rbf_vec_coeff_e(4,je,jb)                       &
              * p_nh%prog(nnew)%vn(iqidx(je,jb,4),jk,iqblk(je,jb,4))

           ENDDO
        ENDDO

      ELSE IF (itime_scheme >= 5) THEN

#ifdef __LOOP_EXCHANGE
        DO je = i_startidx, i_endidx
!DIR$ IVDEP
          DO jk = 1, nlev
#else
        DO jk = 1, nlev
          DO je = i_startidx, i_endidx
#endif
            ! Average normal wind components in order to get nearly second-order accurate divergence
            z_vn_avg(je,jk) = p_int%e_flx_avg(je,1,jb)*p_nh%prog(nnew)%vn(je,jk,jb)           &
              + p_int%e_flx_avg(je,2,jb)*p_nh%prog(nnew)%vn(iqidx(je,jb,1),jk,iqblk(je,jb,1)) &
              + p_int%e_flx_avg(je,3,jb)*p_nh%prog(nnew)%vn(iqidx(je,jb,2),jk,iqblk(je,jb,2)) &
              + p_int%e_flx_avg(je,4,jb)*p_nh%prog(nnew)%vn(iqidx(je,jb,3),jk,iqblk(je,jb,3)) &
              + p_int%e_flx_avg(je,5,jb)*p_nh%prog(nnew)%vn(iqidx(je,jb,4),jk,iqblk(je,jb,4))

            ! RBF reconstruction of tangential wind component
            p_nh%diag%vt(je,jk,jb) = p_int%rbf_vec_coeff_e(1,je,jb)  &
              * p_nh%prog(nnew)%vn(iqidx(je,jb,1),jk,iqblk(je,jb,1)) &
              + p_int%rbf_vec_coeff_e(2,je,jb)                       &
              * p_nh%prog(nnew)%vn(iqidx(je,jb,2),jk,iqblk(je,jb,2)) &
              + p_int%rbf_vec_coeff_e(3,je,jb)                       &
              * p_nh%prog(nnew)%vn(iqidx(je,jb,3),jk,iqblk(je,jb,3)) &
              + p_int%rbf_vec_coeff_e(4,je,jb)                       &
              * p_nh%prog(nnew)%vn(iqidx(je,jb,4),jk,iqblk(je,jb,4))

           ENDDO
        ENDDO

      ELSE

#ifdef __LOOP_EXCHANGE
        DO je = i_startidx, i_endidx
!DIR$ IVDEP
          DO jk = 1, nlev
#else
        DO jk = 1, nlev
          DO je = i_startidx, i_endidx
#endif
            ! Average normal wind components in order to get nearly second-order accurate divergence
            z_vn_avg(je,jk) = p_int%e_flx_avg(je,1,jb)*p_nh%prog(nnew)%vn(je,jk,jb)           &
              + p_int%e_flx_avg(je,2,jb)*p_nh%prog(nnew)%vn(iqidx(je,jb,1),jk,iqblk(je,jb,1)) &
              + p_int%e_flx_avg(je,3,jb)*p_nh%prog(nnew)%vn(iqidx(je,jb,2),jk,iqblk(je,jb,2)) &
              + p_int%e_flx_avg(je,4,jb)*p_nh%prog(nnew)%vn(iqidx(je,jb,3),jk,iqblk(je,jb,3)) &
              + p_int%e_flx_avg(je,5,jb)*p_nh%prog(nnew)%vn(iqidx(je,jb,4),jk,iqblk(je,jb,4))
           
           ENDDO
        ENDDO

      ENDIF

      IF (idiv_method == 1) THEN  ! Compute fluxes at edges using averaged velocities
                                  ! corresponding computation for idiv_method=2 follows later
        DO jk = 1,nlev
!DIR$ IVDEP
          DO je = i_startidx, i_endidx

            p_nh%diag%mass_fl_e(je,jk,jb) = z_rho_e(je,jk,jb) *        &
              z_vn_avg(je,jk) * p_nh%metrics%ddqz_z_full_e(je,jk,jb)
            z_theta_v_fl_e(je,jk,jb) = p_nh%diag%mass_fl_e(je,jk,jb) * &
              z_theta_v_e(je,jk,jb)

          ENDDO
        ENDDO

        IF (lsave_mflx .AND. istep == 2) THEN ! store mass flux for nest boundary interpolation
          DO je = i_startidx, i_endidx
            IF (p_patch%edges%refin_ctrl(je,jb) <= -4 .AND. p_patch%edges%refin_ctrl(je,jb) >= -6) THEN
!DIR$ IVDEP
              p_nh%diag%mass_fl_e_sv(je,:,jb) = p_nh%diag%mass_fl_e(je,:,jb)
            ENDIF
          ENDDO
        ENDIF

        IF (lprep_adv .AND. istep == 2) THEN ! Preprations for tracer advection
          IF (lclean_mflx) THEN
            prep_adv%mass_flx_me(:,:,jb) = 0._wp
            prep_adv%vn_traj    (:,:,jb) = 0._wp
          ENDIF
          DO jk = 1, nlev
            DO je = i_startidx, i_endidx
              prep_adv%vn_traj(je,jk,jb)     = prep_adv%vn_traj(je,jk,jb)     + r_iadv_rcf*z_vn_avg(je,jk)
              prep_adv%mass_flx_me(je,jk,jb) = prep_adv%mass_flx_me(je,jk,jb) + r_iadv_rcf*p_nh%diag%mass_fl_e(je,jk,jb)
            ENDDO
          ENDDO
        ENDIF

      ENDIF

      IF (istep == 1 .OR. itime_scheme >= 5) THEN
        ! Compute contravariant correction for vertical velocity at full levels
        DO jk = nflatlev(p_patch%id), nlev
!DIR$ IVDEP
          DO je = i_startidx, i_endidx
            z_w_concorr_me(je,jk,jb) =                                          &
              p_nh%prog(nnew)%vn(je,jk,jb)*p_nh%metrics%ddxn_z_full(je,jk,jb) + &
              p_nh%diag%vt(je,jk,jb)      *p_nh%metrics%ddxt_z_full(je,jk,jb)
          ENDDO
        ENDDO
      ENDIF

      IF (istep == 1) THEN
        ! Interpolate vn to interface levels and compute horizontal part of kinetic energy on edges
        ! (needed in velocity tendencies called at istep=2)
        DO jk = 2, nlev
!DIR$ IVDEP
          DO je = i_startidx, i_endidx
            p_nh%diag%vn_ie(je,jk,jb) =                                             &
              p_nh%metrics%wgtfac_e(je,jk,jb)*p_nh%prog(nnew)%vn(je,jk,jb) +        &
             (1._wp - p_nh%metrics%wgtfac_e(je,jk,jb))*p_nh%prog(nnew)%vn(je,jk-1,jb)
            z_vt_ie(je,jk,jb) =                                               &
              p_nh%metrics%wgtfac_e(je,jk,jb)*p_nh%diag%vt(je,jk,jb) +        &
             (1._wp - p_nh%metrics%wgtfac_e(je,jk,jb))*p_nh%diag%vt(je,jk-1,jb)
            z_kin_hor_e(je,jk,jb) = 0.5_wp*(p_nh%prog(nnew)%vn(je,jk,jb)**2 + p_nh%diag%vt(je,jk,jb)**2)
          ENDDO
        ENDDO

        IF (.NOT. l_vert_nested) THEN
          ! Top and bottom levels
!DIR$ IVDEP
          DO je = i_startidx, i_endidx
            ! Quadratic extrapolation at the top turned out to cause numerical instability in pathological cases,
            ! thus we use a no-gradient condition in the upper half layer
            p_nh%diag%vn_ie(je,1,jb) = p_nh%prog(nnew)%vn(je,1,jb)
            ! vt_ie(jk=1) is actually unused, but we need it for convenience of implementation
            z_vt_ie(je,1,jb) = p_nh%diag%vt(je,1,jb)
            !
            z_kin_hor_e(je,1,jb) = 0.5_wp*(p_nh%prog(nnew)%vn(je,1,jb)**2 + p_nh%diag%vt(je,1,jb)**2)
            p_nh%diag%vn_ie(je,nlevp1,jb) =                           &
              p_nh%metrics%wgtfacq_e(je,1,jb)*p_nh%prog(nnew)%vn(je,nlev,jb) +   &
              p_nh%metrics%wgtfacq_e(je,2,jb)*p_nh%prog(nnew)%vn(je,nlev-1,jb) + &
              p_nh%metrics%wgtfacq_e(je,3,jb)*p_nh%prog(nnew)%vn(je,nlev-2,jb)
          ENDDO
        ELSE
          ! vn_ie(jk=1) is extrapolated using parent domain information in this case
!DIR$ IVDEP
          DO je = i_startidx, i_endidx
            p_nh%diag%vn_ie(je,1,jb) = p_nh%diag%vn_ie(je,2,jb) + p_nh%diag%dvn_ie_ubc(je,jb)
            ! vt_ie(jk=1) is actually unused, but we need it for convenience of implementation
            z_vt_ie(je,1,jb) = p_nh%diag%vt(je,1,jb)
            !
            z_kin_hor_e(je,1,jb) = 0.5_wp*(p_nh%prog(nnew)%vn(je,1,jb)**2 + p_nh%diag%vt(je,1,jb)**2)
            p_nh%diag%vn_ie(je,nlevp1,jb) =                           &
              p_nh%metrics%wgtfacq_e(je,1,jb)*p_nh%prog(nnew)%vn(je,nlev,jb) +   &
              p_nh%metrics%wgtfacq_e(je,2,jb)*p_nh%prog(nnew)%vn(je,nlev-1,jb) + &
              p_nh%metrics%wgtfacq_e(je,3,jb)*p_nh%prog(nnew)%vn(je,nlev-2,jb)
          ENDDO
        ENDIF
      ENDIF

    ENDDO
!$OMP END DO

    ! Apply mass fluxes across lateral nest boundary interpolated from parent domain
    IF (p_patch%id > 1 .AND. grf_intmethod_e >= 5 .AND. idiv_method == 1) THEN
!$OMP DO PRIVATE(ic,je,jb,jk) ICON_OMP_DEFAULT_SCHEDULE
      DO ic = 1, p_nh%metrics%bdy_mflx_e_dim
        je = p_nh%metrics%bdy_mflx_e_idx(ic)
        jb = p_nh%metrics%bdy_mflx_e_blk(ic)

        ! This is needed for tracer mass consistency along the lateral boundaries
        IF (lprep_adv .AND. istep == 2) THEN ! subtract mass flux added previously...
          DO jk = 1, nlev
            prep_adv%mass_flx_me(je,jk,jb) = prep_adv%mass_flx_me(je,jk,jb) - r_iadv_rcf*p_nh%diag%mass_fl_e(je,jk,jb)
            prep_adv%vn_traj(je,jk,jb)     = prep_adv%vn_traj(je,jk,jb) - r_iadv_rcf*p_nh%diag%mass_fl_e(je,jk,jb) / &
                                             (z_rho_e(je,jk,jb) * p_nh%metrics%ddqz_z_full_e(je,jk,jb))
          ENDDO
        ENDIF

!DIR$ IVDEP
        DO jk = 1, nlev
          p_nh%diag%mass_fl_e(je,jk,jb) = p_nh%diag%grf_bdy_mflx(jk,ic,1) + &
            REAL(jstep,wp)*dtime*p_nh%diag%grf_bdy_mflx(jk,ic,2)
          z_theta_v_fl_e(je,jk,jb) = p_nh%diag%mass_fl_e(je,jk,jb) * z_theta_v_e(je,jk,jb)
        ENDDO

        IF (lprep_adv .AND. istep == 2) THEN ! ... and add the corrected one again
          DO jk = 1, nlev
            prep_adv%mass_flx_me(je,jk,jb) = prep_adv%mass_flx_me(je,jk,jb) + r_iadv_rcf*p_nh%diag%mass_fl_e(je,jk,jb)
            prep_adv%vn_traj(je,jk,jb)     = prep_adv%vn_traj(je,jk,jb) + r_iadv_rcf*p_nh%diag%mass_fl_e(je,jk,jb) / &
                                             (z_rho_e(je,jk,jb) * p_nh%metrics%ddqz_z_full_e(je,jk,jb))
          ENDDO
        ENDIF

      ENDDO
!$OMP END DO
    ENDIF


    ! It turned out that it is sufficient to compute the contravariant correction in the
    ! predictor step at time level n+1; repeating the calculation in the corrector step
    ! has negligible impact on the results except in very-high resolution runs with extremely steep mountains
    IF (istep == 1 .OR. itime_scheme >= 5) THEN

      rl_start = 3
      rl_end = min_rlcell_int - 1

      i_startblk = p_patch%cells%start_block(rl_start)
      i_endblk   = p_patch%cells%end_block(rl_end)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,jc,z_w_concorr_mc) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, rl_start, rl_end)

        ! Interpolate contravariant correction to cell centers...
#ifdef __LOOP_EXCHANGE
        DO jc = i_startidx, i_endidx
!DIR$ IVDEP
          DO jk = nflatlev(p_patch%id), nlev
#else
        DO jk = nflatlev(p_patch%id), nlev
          DO jc = i_startidx, i_endidx
#endif

            z_w_concorr_mc(jc,jk) =  &
              p_int%e_bln_c_s(jc,1,jb)*z_w_concorr_me(ieidx(jc,jb,1),jk,ieblk(jc,jb,1)) + &
              p_int%e_bln_c_s(jc,2,jb)*z_w_concorr_me(ieidx(jc,jb,2),jk,ieblk(jc,jb,2)) + &
              p_int%e_bln_c_s(jc,3,jb)*z_w_concorr_me(ieidx(jc,jb,3),jk,ieblk(jc,jb,3))

          ENDDO
        ENDDO

        ! ... and to interface levels
        DO jk = nflatlev(p_patch%id)+1, nlev
!DIR$ IVDEP
          DO jc = i_startidx, i_endidx
            p_nh%diag%w_concorr_c(jc,jk,jb) =                                &
              p_nh%metrics%wgtfac_c(jc,jk,jb)*z_w_concorr_mc(jc,jk) +        &
             (1._vp - p_nh%metrics%wgtfac_c(jc,jk,jb))*z_w_concorr_mc(jc,jk-1) 
          ENDDO
        ENDDO
!DIR$ IVDEP
        DO jc = i_startidx, i_endidx
          p_nh%diag%w_concorr_c(jc,nlevp1,jb) =                         &
            p_nh%metrics%wgtfacq_c(jc,1,jb)*z_w_concorr_mc(jc,nlev) +   &
            p_nh%metrics%wgtfacq_c(jc,2,jb)*z_w_concorr_mc(jc,nlev-1) + &
            p_nh%metrics%wgtfacq_c(jc,3,jb)*z_w_concorr_mc(jc,nlev-2)
        ENDDO

      ENDDO
!$OMP END DO
    ENDIF
 
    IF (idiv_method == 2) THEN ! Compute fluxes at edges from original velocities
      rl_start = 7
      rl_end = min_rledge_int - 3

      i_startblk = p_patch%edges%start_block(rl_start)
      i_endblk   = p_patch%edges%end_block(rl_end)

      IF (p_patch%id > 1 .OR. l_limited_area) THEN
!$OMP WORKSHARE
        z_theta_v_fl_e(:,:,p_patch%edges%start_block(5):i_startblk) = 0._wp
!$OMP END WORKSHARE
      ENDIF

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,je) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, rl_start, rl_end)

        DO jk = 1,nlev
!DIR$ IVDEP
          DO je = i_startidx, i_endidx

            p_nh%diag%mass_fl_e(je,jk,jb) = z_rho_e(je,jk,jb)         &
              * p_nh%prog(nnew)%vn(je,jk,jb) * p_nh%metrics%ddqz_z_full_e(je,jk,jb)
            z_theta_v_fl_e(je,jk,jb)= p_nh%diag%mass_fl_e(je,jk,jb)   &
              * z_theta_v_e(je,jk,jb)

          ENDDO
        ENDDO

      ENDDO
!$OMP END DO
    ENDIF  ! idiv_method = 2

!$OMP END PARALLEL

    IF (timers_level > 5) THEN
      CALL timer_stop(timer_solve_nh_edgecomp)
      CALL timer_start(timer_solve_nh_vimpl)
    ENDIF

    IF (idiv_method == 2) THEN ! use averaged divergence - idiv_method=1 is inlined for better cache efficiency

      ! horizontal divergences of rho and rhotheta are processed in one step for efficiency
      CALL div_avg(p_nh%diag%mass_fl_e, p_patch, p_int, p_int%c_bln_avg, z_mass_fl_div, &
                   opt_in2=z_theta_v_fl_e, opt_out2=z_theta_v_fl_div, opt_rlstart=4,    &
                   opt_rlend=min_rlcell_int)
    ENDIF

!$OMP PARALLEL PRIVATE (rl_start,rl_end,i_startblk,i_endblk,jk_start)

    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

    IF (l_vert_nested) THEN
      jk_start = 2
    ELSE
      jk_start = 1
    ENDIF

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,jc,z_w_expl,z_contr_w_fl_l,z_rho_expl,z_exner_expl, &
!$OMP   z_a,z_b,z_c,z_g,z_q,z_alpha,z_beta,z_gamma,ic,z_raylfac,z_flxdiv_mass,z_flxdiv_theta  ) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)

      IF (idiv_method == 1) THEN
        ! horizontal divergences of rho and rhotheta are inlined and processed in one step for efficiency
#ifdef __LOOP_EXCHANGE
        DO jc = i_startidx, i_endidx
!DIR$ IVDEP, PREFERVECTOR
          DO jk = 1, nlev
#else
        DO jk = 1, nlev
          DO jc = i_startidx, i_endidx
#endif

            z_flxdiv_mass(jc,jk) =  &
              p_nh%diag%mass_fl_e(ieidx(jc,jb,1),jk,ieblk(jc,jb,1)) * p_int%geofac_div(jc,1,jb) + &
              p_nh%diag%mass_fl_e(ieidx(jc,jb,2),jk,ieblk(jc,jb,2)) * p_int%geofac_div(jc,2,jb) + &
              p_nh%diag%mass_fl_e(ieidx(jc,jb,3),jk,ieblk(jc,jb,3)) * p_int%geofac_div(jc,3,jb)

            z_flxdiv_theta(jc,jk) =  &
              z_theta_v_fl_e(ieidx(jc,jb,1),jk,ieblk(jc,jb,1)) * p_int%geofac_div(jc,1,jb) + &
              z_theta_v_fl_e(ieidx(jc,jb,2),jk,ieblk(jc,jb,2)) * p_int%geofac_div(jc,2,jb) + &
              z_theta_v_fl_e(ieidx(jc,jb,3),jk,ieblk(jc,jb,3)) * p_int%geofac_div(jc,3,jb)

          END DO
        END DO

      ELSE ! idiv_method = 2 - just copy values to local 2D array

        DO jk = 1, nlev
          DO jc = i_startidx, i_endidx
            z_flxdiv_mass(jc,jk)  = z_mass_fl_div(jc,jk,jb)
            z_flxdiv_theta(jc,jk) = z_theta_v_fl_div(jc,jk,jb)
          END DO
        END DO

      ENDIF

      ! upper boundary conditions for rho_ic and theta_v_ic in the case of vertical nesting
      IF (l_vert_nested .AND. istep == 1) THEN
!DIR$ IVDEP
        DO jc = i_startidx, i_endidx
          p_nh%diag%theta_v_ic(jc,1,jb) = p_nh%diag%theta_v_ic(jc,2,jb) + &
            p_nh%diag%dtheta_v_ic_ubc(jc,jb)
          z_mflx_top(jc,jb)             = p_nh%diag%mflx_ic_ubc(jc,jb,1) + &
            REAL(jstep,wp)*dtime*p_nh%diag%mflx_ic_ubc(jc,jb,2)
          p_nh%diag%rho_ic(jc,1,jb) =  0._wp ! not used in dynamical core in this case, will be set for tracer interface later
        ENDDO
      ELSE IF (l_vert_nested .AND. istep == 2) THEN
!DIR$ IVDEP
        DO jc = i_startidx, i_endidx
          p_nh%diag%theta_v_ic(jc,1,jb) = p_nh%diag%theta_v_ic(jc,2,jb) + &
            p_nh%diag%dtheta_v_ic_ubc(jc,jb)
          z_mflx_top(jc,jb)             = p_nh%diag%mflx_ic_ubc(jc,jb,1) + &
            (REAL(jstep,wp)+0.5_wp)*dtime*p_nh%diag%mflx_ic_ubc(jc,jb,2)
        ENDDO
      ENDIF

      ! Start of vertically implicit solver part for sound-wave terms;
      ! advective terms and gravity-wave terms are treated explicitly
      !
      IF (istep == 2 .AND. (itime_scheme >= 4)) THEN
        DO jk = 2, nlev
!DIR$ IVDEP
          DO jc = i_startidx, i_endidx

            ! explicit part for w - use temporally averaged advection terms for better numerical stability
            ! the explicit weight for the pressure-gradient term is already included in z_th_ddz_exner_c
            z_w_expl(jc,jk) = p_nh%prog(nnow)%w(jc,jk,jb) + dtime *   &
              (wgt_nnow_vel*p_nh%diag%ddt_w_adv(jc,jk,jb,ntl1) +      &
               wgt_nnew_vel*p_nh%diag%ddt_w_adv(jc,jk,jb,ntl2)        &
              -cpd*z_th_ddz_exner_c(jc,jk,jb) )

            ! contravariant vertical velocity times density for explicit part
            z_contr_w_fl_l(jc,jk) = p_nh%diag%rho_ic(jc,jk,jb)*(-p_nh%diag%w_concorr_c(jc,jk,jb) &
              + p_nh%metrics%vwind_expl_wgt(jc,jb)*p_nh%prog(nnow)%w(jc,jk,jb) )

          ENDDO
        ENDDO
      ELSE
        DO jk = 2, nlev
!DIR$ IVDEP
          DO jc = i_startidx, i_endidx

            ! explicit part for w
            z_w_expl(jc,jk) = p_nh%prog(nnow)%w(jc,jk,jb) + dtime *             &
              (p_nh%diag%ddt_w_adv(jc,jk,jb,ntl1)-cpd*z_th_ddz_exner_c(jc,jk,jb))

            ! contravariant vertical velocity times density for explicit part
            z_contr_w_fl_l(jc,jk) = p_nh%diag%rho_ic(jc,jk,jb)*(-p_nh%diag%w_concorr_c(jc,jk,jb) &
              + p_nh%metrics%vwind_expl_wgt(jc,jb)*p_nh%prog(nnow)%w(jc,jk,jb) )

          ENDDO
        ENDDO
      ENDIF

      ! Solver coefficients
      DO jk = 1, nlev
!DIR$ IVDEP
        DO jc = i_startidx, i_endidx
          z_beta(jc,jk)=dtime*rd*p_nh%prog(nnow)%exner(jc,jk,jb) /                 &
           (cvd*p_nh%prog(nnow)%rho(jc,jk,jb)*p_nh%prog(nnow)%theta_v(jc,jk,jb)) * &
            p_nh%metrics%inv_ddqz_z_full(jc,jk,jb)

          z_alpha(jc,jk)= p_nh%metrics%vwind_impl_wgt(jc,jb)*         &
            &  p_nh%diag%theta_v_ic(jc,jk,jb)*p_nh%diag%rho_ic(jc,jk,jb)
        ENDDO
      ENDDO

      z_alpha(:,nlevp1) = 0.0_wp
      DO jk = 2, nlev
!DIR$ IVDEP
        DO jc = i_startidx, i_endidx
          z_gamma(jc,jk) = dtime*cpd*p_nh%metrics%vwind_impl_wgt(jc,jb)*    &
            p_nh%diag%theta_v_ic(jc,jk,jb)/p_nh%metrics%ddqz_z_half(jc,jk,jb)

          z_a(jc,jk) = -z_gamma(jc,jk)*z_beta(jc,jk-1)*z_alpha(jc,jk-1)
          z_c(jc,jk) = -z_gamma(jc,jk)*z_beta(jc,jk  )*z_alpha(jc,jk+1)
          z_b(jc,jk) = 1.0_vp+z_gamma(jc,jk)*z_alpha(jc,jk) &
            *(z_beta(jc,jk-1)+z_beta(jc,jk))
        ENDDO
      ENDDO
!DIR$ IVDEP
      DO jc = i_startidx, i_endidx
        z_q(jc,2) = -z_c(jc,2)/z_b(jc,2)
      ENDDO

      DO jk = 3, nlev
!DIR$ IVDEP
        DO jc = i_startidx, i_endidx
          z_g(jc,jk) = 1.0_vp/(z_b(jc,jk)+z_a(jc,jk)*z_q(jc,jk-1))
          z_q(jc,jk) = - z_c(jc,jk)*z_g(jc,jk)
        ENDDO
      ENDDO

      ! upper boundary condition for w (interpolated from parent domain in case of vertical nesting)
      ! Note: the upper b.c. reduces to w(1) = 0 in the absence of diabatic heating
      IF (l_open_ubc .AND. .NOT. l_vert_nested) THEN
!DIR$ IVDEP
        DO jc = i_startidx, i_endidx
          p_nh%prog(nnew)%w(jc,1,jb) = z_thermal_exp(jc,jb)
          z_contr_w_fl_l(jc,1) = p_nh%diag%rho_ic(jc,1,jb)*p_nh%prog(nnow)%w(jc,1,jb)   &
            * p_nh%metrics%vwind_expl_wgt(jc,jb)
        ENDDO
      ELSE IF (.NOT. l_open_ubc .AND. .NOT. l_vert_nested) THEN
        p_nh%prog(nnew)%w(:,1,jb) = 0._wp
        z_contr_w_fl_l(:,1)       = 0._wp
      ELSE  ! l_vert_nested
!DIR$ IVDEP
        DO jc = i_startidx, i_endidx
          z_contr_w_fl_l(jc,1) = z_mflx_top(jc,jb) * p_nh%metrics%vwind_expl_wgt(jc,jb)
        ENDDO
      ENDIF

      ! lower boundary condition for w, consistent with contravariant correction
!DIR$ IVDEP
      DO jc = i_startidx, i_endidx
        p_nh%prog(nnew)%w(jc,nlevp1,jb) = p_nh%diag%w_concorr_c(jc,nlevp1,jb)
        z_contr_w_fl_l(jc,nlevp1)       = 0.0_wp
      ENDDO


      ! Explicit parts of density and Exner pressure
      !
      ! Top level first
!DIR$ IVDEP
      DO jc = i_startidx, i_endidx
        z_rho_expl(jc,1)=        p_nh%prog(nnow)%rho(jc,1,jb) &
        &        -dtime*p_nh%metrics%inv_ddqz_z_full(jc,1,jb) &
        &                            *(z_flxdiv_mass(jc,1)    &
        &                            +z_contr_w_fl_l(jc,1   ) &
        &                            -z_contr_w_fl_l(jc,2   ))

        z_exner_expl(jc,1)=             z_exner_pr(jc,1,jb)    &
        &      -z_beta (jc,1)*(z_flxdiv_theta(jc,1)            &
        & +p_nh%diag%theta_v_ic(jc,1,jb)*z_contr_w_fl_l(jc,1)  &
        & -p_nh%diag%theta_v_ic(jc,2,jb)*z_contr_w_fl_l(jc,2)) &
        & +dtime*p_nh%diag%ddt_exner_phy(jc,1,jb)
      ENDDO

      ! Other levels
      DO jk = 2, nlev
!DIR$ IVDEP
        DO jc = i_startidx, i_endidx
          z_rho_expl(jc,jk)=       p_nh%prog(nnow)%rho(jc,jk  ,jb) &
          &        -dtime*p_nh%metrics%inv_ddqz_z_full(jc,jk  ,jb) &
          &                            *(z_flxdiv_mass(jc,jk     ) &
          &                            +z_contr_w_fl_l(jc,jk     ) &
          &                             -z_contr_w_fl_l(jc,jk+1   ))

          z_exner_expl(jc,jk)=          z_exner_pr(jc,jk,jb) - z_beta(jc,jk) &
          &                             *(z_flxdiv_theta(jc,jk)              &
          &   +p_nh%diag%theta_v_ic(jc,jk  ,jb)*z_contr_w_fl_l(jc,jk  )      &
          &   -p_nh%diag%theta_v_ic(jc,jk+1,jb)*z_contr_w_fl_l(jc,jk+1))     &
          &   +dtime*p_nh%diag%ddt_exner_phy(jc,jk,jb)

          p_nh%prog(nnew)%w(jc,jk,jb) = z_w_expl(jc,jk) - z_gamma(jc,jk)  &
          &      *(z_exner_expl(jc,jk-1)-z_exner_expl(jc,jk))
        ENDDO
      ENDDO

      ! Solve tridiagonal matrix for w
!DIR$ IVDEP
      DO jc = i_startidx, i_endidx
        p_nh%prog(nnew)%w(jc,2,jb)= p_nh%prog(nnew)%w(jc,2,jb)/z_b(jc,2)
      ENDDO

      DO jk = 3, nlev
!DIR$ IVDEP
        DO jc = i_startidx, i_endidx
          p_nh%prog(nnew)%w(jc,jk,jb) = (p_nh%prog(nnew)%w(jc,jk,jb)  &
            -z_a(jc,jk)*p_nh%prog(nnew)%w(jc,jk-1,jb))*z_g(jc,jk)
        ENDDO
      ENDDO

      DO jk = nlev-1, 2, -1
!DIR$ IVDEP
        DO jc = i_startidx, i_endidx
          p_nh%prog(nnew)%w(jc,jk,jb) = p_nh%prog(nnew)%w(jc,jk,jb)&
          &             +p_nh%prog(nnew)%w(jc,jk+1,jb)*z_q(jc,jk)
        ENDDO
      ENDDO

      IF (l_vert_nested) THEN
!DIR$ IVDEP
        DO jc = i_startidx, i_endidx
          p_nh%prog(nnew)%w(jc,1,jb) = p_nh%prog(nnew)%w(jc,2,jb) + p_nh%diag%dw_ubc(jc,jb)
        ENDDO
      ENDIF

      ! Rayleigh damping mechanism (Klemp,Dudhia,Hassiotis: MWR136,pp.3987-4004)
      !
      IF ( rayleigh_type == RAYLEIGH_KLEMP ) THEN
        DO jk = 2, nrdmax(p_patch%id)
          z_raylfac = 1.0_wp/(1.0_wp+dtime*p_nh%metrics%rayleigh_w(jk))
!DIR$ IVDEP
          DO jc = i_startidx, i_endidx
            p_nh%prog(nnew)%w(jc,jk,jb) = z_raylfac*p_nh%prog(nnew)%w(jc,jk,jb) +    &
                                          (1._wp-z_raylfac)*p_nh%prog(nnew)%w(jc,1,jb)
          ENDDO
        ENDDO
      ! Classic Rayleigh damping mechanism for w (requires reference state !!)
      !
      ELSE IF ( rayleigh_type == RAYLEIGH_CLASSIC ) THEN
        DO jk = 2, nrdmax(p_patch%id)
!DIR$ IVDEP
          DO jc = i_startidx, i_endidx
            p_nh%prog(nnew)%w(jc,jk,jb) = p_nh%prog(nnew)%w(jc,jk,jb)       &
              &                         - dtime*p_nh%metrics%rayleigh_w(jk) &
              &                         * ( p_nh%prog(nnew)%w(jc,jk,jb)     &
              &                         - p_nh%ref%w_ref(jc,jk,jb) )
          ENDDO
        ENDDO
      ENDIF

      IF (is_iau_active) THEN ! add analysis increments from data assimilation to density and exner pressure
        DO jk = 1, nlev
!DIR$ IVDEP
          DO jc = i_startidx, i_endidx
            z_rho_expl(jc,jk)   = z_rho_expl(jc,jk)   + iau_wgt_dyn*p_nh%diag%rho_incr(jc,jk,jb)
            z_exner_expl(jc,jk) = z_exner_expl(jc,jk) + iau_wgt_dyn*p_nh%diag%exner_incr(jc,jk,jb)
          ENDDO
        ENDDO
      ENDIF

      ! Results for thermodynamic variables
      DO jk = jk_start, nlev
!DIR$ IVDEP
        DO jc = i_startidx, i_endidx

          ! density
          p_nh%prog(nnew)%rho(jc,jk,jb) = z_rho_expl(jc,jk)              &
            - p_nh%metrics%vwind_impl_wgt(jc,jb)*dtime                   &
            * p_nh%metrics%inv_ddqz_z_full(jc,jk,jb)                     &
            *(p_nh%diag%rho_ic(jc,jk  ,jb)*p_nh%prog(nnew)%w(jc,jk  ,jb) &
            - p_nh%diag%rho_ic(jc,jk+1,jb)*p_nh%prog(nnew)%w(jc,jk+1,jb))

          ! exner
          p_nh%prog(nnew)%exner(jc,jk,jb) = z_exner_expl(jc,jk) &
            + p_nh%metrics%exner_ref_mc(jc,jk,jb)-z_beta(jc,jk) &
            *(z_alpha(jc,jk  )*p_nh%prog(nnew)%w(jc,jk  ,jb)    &
            - z_alpha(jc,jk+1)*p_nh%prog(nnew)%w(jc,jk+1,jb))

          ! theta
          p_nh%prog(nnew)%theta_v(jc,jk,jb) = p_nh%prog(nnow)%rho(jc,jk,jb)*p_nh%prog(nnow)%theta_v(jc,jk,jb) &
            *( (p_nh%prog(nnew)%exner(jc,jk,jb)/p_nh%prog(nnow)%exner(jc,jk,jb)-1.0_wp) * cvd_o_rd+1.0_wp   ) &
            / p_nh%prog(nnew)%rho(jc,jk,jb)

        ENDDO
      ENDDO

      ! Special treatment of uppermost layer in the case of vertical nesting
      IF (l_vert_nested) THEN
!DIR$ IVDEP
        DO jc = i_startidx, i_endidx

          ! density
          p_nh%prog(nnew)%rho(jc,1,jb) = z_rho_expl(jc,1)                             &
            - p_nh%metrics%vwind_impl_wgt(jc,jb)*dtime                                &
            * p_nh%metrics%inv_ddqz_z_full(jc,1,jb)                                   &
            *(z_mflx_top(jc,jb) - p_nh%diag%rho_ic(jc,2,jb)*p_nh%prog(nnew)%w(jc,2,jb))

          ! exner
          p_nh%prog(nnew)%exner(jc,1,jb) = z_exner_expl(jc,1)                  &
            + p_nh%metrics%exner_ref_mc(jc,1,jb)-z_beta(jc,1)                  &
            *(p_nh%metrics%vwind_impl_wgt(jc,jb)*p_nh%diag%theta_v_ic(jc,1,jb) &
            * z_mflx_top(jc,jb) - z_alpha(jc,2)*p_nh%prog(nnew)%w(jc,2,jb))

          ! theta
          p_nh%prog(nnew)%theta_v(jc,1,jb) = p_nh%prog(nnow)%rho(jc,1,jb)*p_nh%prog(nnow)%theta_v(jc,1,jb) &
            *( (p_nh%prog(nnew)%exner(jc,1,jb)/p_nh%prog(nnow)%exner(jc,1,jb)-1.0_wp) * cvd_o_rd+1.0_wp  ) &
            /p_nh%prog(nnew)%rho(jc,1,jb)

        ENDDO
      ENDIF

      IF (istep == 2 .AND. l_vert_nested) THEN
        ! Diagnose rho_ic(jk=1) for tracer transport, and rediagnose appropriate w(jk=1)
!DIR$ IVDEP
        DO jc = i_startidx, i_endidx
          p_nh%diag%rho_ic(jc,1,jb) =  wgt_nnow_rth*(                        &
            p_nh%metrics%wgtfacq1_c(jc,1,jb)*p_nh%prog(nnow)%rho(jc,1,jb) +  &
            p_nh%metrics%wgtfacq1_c(jc,2,jb)*p_nh%prog(nnow)%rho(jc,2,jb) +  &
            p_nh%metrics%wgtfacq1_c(jc,3,jb)*p_nh%prog(nnow)%rho(jc,3,jb))+  &
            wgt_nnew_rth * (                                                 &
            p_nh%metrics%wgtfacq1_c(jc,1,jb)*p_nh%prog(nvar)%rho(jc,1,jb) +  &
            p_nh%metrics%wgtfacq1_c(jc,2,jb)*p_nh%prog(nvar)%rho(jc,2,jb) +  &
            p_nh%metrics%wgtfacq1_c(jc,3,jb)*p_nh%prog(nvar)%rho(jc,3,jb) )
          p_nh%prog(nnew)%w(jc,1,jb) = z_mflx_top(jc,jb)/p_nh%diag%rho_ic(jc,1,jb)
        ENDDO
      ENDIF

      ! compute dw/dz for divergence damping term
      IF (lhdiff_rcf .AND. istep == 1 .AND. divdamp_type >= 3) THEN
        DO jk = 1, nlev
!DIR$ IVDEP
          DO jc = i_startidx, i_endidx
            z_dwdz_dd(jc,jk,jb) = p_nh%metrics%inv_ddqz_z_full(jc,jk,jb) *          &
              ( (p_nh%prog(nnew)%w(jc,jk,jb)-p_nh%prog(nnew)%w(jc,jk+1,jb)) -       &
                (p_nh%diag%w_concorr_c(jc,jk,jb)-p_nh%diag%w_concorr_c(jc,jk+1,jb)) )
          ENDDO
        ENDDO
      ENDIF

      ! Preparations for tracer advection
      IF (lprep_adv .AND. istep == 2) THEN
        IF (lclean_mflx) prep_adv%mass_flx_ic(:,:,jb) = 0._wp
        DO jk = 1, nlev
          DO jc = i_startidx, i_endidx
            prep_adv%mass_flx_ic(jc,jk,jb) = prep_adv%mass_flx_ic(jc,jk,jb) + r_iadv_rcf * ( z_contr_w_fl_l(jc,jk) + &
              p_nh%diag%rho_ic(jc,jk,jb) * p_nh%metrics%vwind_impl_wgt(jc,jb) * p_nh%prog(nnew)%w(jc,jk,jb) )
          ENDDO
        ENDDO
      ENDIF

      IF (istep == 2) THEN
        ! store dynamical part of exner time increment in exner_dyn_incr
        ! the conversion into a temperature tendency is done in the NWP interface
        DO jk = kstart_moist(p_patch%id), nlev
!DIR$ IVDEP
          DO jc = i_startidx, i_endidx
            p_nh%diag%exner_dyn_incr(jc,jk,jb) = p_nh%diag%exner_dyn_incr(jc,jk,jb) + &
              p_nh%prog(nnew)%exner(jc,jk,jb) - p_nh%prog(nnow)%exner(jc,jk,jb) -     &
              dtime*p_nh%diag%ddt_exner_phy(jc,jk,jb)
          ENDDO
        ENDDO
      ENDIF

      IF (istep == 2 .AND. l_child_vertnest) THEN
        ! Store values at nest interface levels
!DIR$ IVDEP
        DO jc = i_startidx, i_endidx

          p_nh%diag%dw_int(jc,jb,idyn_timestep) =                                         &
            0.5_wp*(p_nh%prog(nnow)%w(jc,nshift,jb)   + p_nh%prog(nnew)%w(jc,nshift,jb) - &
                   (p_nh%prog(nnow)%w(jc,nshift+1,jb) + p_nh%prog(nnew)%w(jc,nshift+1,jb)))

          p_nh%diag%mflx_ic_int(jc,jb,idyn_timestep) = p_nh%diag%rho_ic(jc,nshift,jb) * &
            (p_nh%metrics%vwind_expl_wgt(jc,jb)*p_nh%prog(nnow)%w(jc,nshift,jb) + &
             p_nh%metrics%vwind_impl_wgt(jc,jb)*p_nh%prog(nnew)%w(jc,nshift,jb))

          p_nh%diag%dtheta_v_ic_int(jc,jb,idyn_timestep) = p_nh%diag%theta_v_ic(jc,nshift,jb) - &
            p_nh%diag%theta_v_ic(jc,nshift+1,jb)
        ENDDO
      ENDIF

    ENDDO
!$OMP END DO

    ! Boundary update in case of nesting
    IF (l_limited_area .OR. p_patch%id > 1) THEN

      rl_start = 1
      rl_end   = grf_bdywidth_c

      i_startblk = p_patch%cells%start_block(rl_start)
      i_endblk   = p_patch%cells%end_block(rl_end)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,jc) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, rl_start, rl_end)

        ! non-MPI-parallelized (serial) case
        IF (istep == 1 .AND. my_process_is_mpi_all_seq() ) THEN

          DO jk = 1, nlev
!DIR$ IVDEP
            DO jc = i_startidx, i_endidx

              p_nh%prog(nnew)%rho(jc,jk,jb) = p_nh%prog(nnow)%rho(jc,jk,jb) + &
                dtime*p_nh%diag%grf_tend_rho(jc,jk,jb)

              p_nh%prog(nnew)%theta_v(jc,jk,jb) = p_nh%prog(nnow)%theta_v(jc,jk,jb) + &
                dtime*p_nh%diag%grf_tend_thv(jc,jk,jb)

              ! Diagnose exner from rho*theta
              p_nh%prog(nnew)%exner(jc,jk,jb) = EXP(rd_o_cvd*LOG(rd_o_p0ref* &
                p_nh%prog(nnew)%rho(jc,jk,jb)*p_nh%prog(nnew)%theta_v(jc,jk,jb)))

              p_nh%prog(nnew)%w(jc,jk,jb) = p_nh%prog(nnow)%w(jc,jk,jb) + &
                dtime*p_nh%diag%grf_tend_w(jc,jk,jb)

            ENDDO
          ENDDO

          DO jc = i_startidx, i_endidx
            p_nh%prog(nnew)%w(jc,nlevp1,jb) = p_nh%prog(nnow)%w(jc,nlevp1,jb) + &
              dtime*p_nh%diag%grf_tend_w(jc,nlevp1,jb)
          ENDDO

        ELSE IF (istep == 1 ) THEN

          ! In the MPI-parallelized case, only rho and w are updated here,
          ! and theta_v is preliminarily stored on exner in order to save
          ! halo communications

          DO jk = 1, nlev
!DIR$ IVDEP
            DO jc = i_startidx, i_endidx

              p_nh%prog(nnew)%rho(jc,jk,jb) = p_nh%prog(nnow)%rho(jc,jk,jb) + &
                dtime*p_nh%diag%grf_tend_rho(jc,jk,jb)

              ! *** Storing theta_v on exner is done to save MPI communications ***
              ! DO NOT TOUCH THIS!
              p_nh%prog(nnew)%exner(jc,jk,jb) = p_nh%prog(nnow)%theta_v(jc,jk,jb) + &
                dtime*p_nh%diag%grf_tend_thv(jc,jk,jb)

              p_nh%prog(nnew)%w(jc,jk,jb) = p_nh%prog(nnow)%w(jc,jk,jb) + &
                dtime*p_nh%diag%grf_tend_w(jc,jk,jb)

            ENDDO
          ENDDO

          DO jc = i_startidx, i_endidx
            p_nh%prog(nnew)%w(jc,nlevp1,jb) = p_nh%prog(nnow)%w(jc,nlevp1,jb) + &
              dtime*p_nh%diag%grf_tend_w(jc,nlevp1,jb)
          ENDDO

        ENDIF

        ! compute dw/dz for divergence damping term
        IF (lhdiff_rcf .AND. istep == 1 .AND. divdamp_type >= 3) THEN
          DO jk = 1, nlev
!DIR$ IVDEP
            DO jc = i_startidx, i_endidx
              z_dwdz_dd(jc,jk,jb) = p_nh%metrics%inv_ddqz_z_full(jc,jk,jb) *          &
                ( (p_nh%prog(nnew)%w(jc,jk,jb)-p_nh%prog(nnew)%w(jc,jk+1,jb)) -       &
                  (p_nh%diag%w_concorr_c(jc,jk,jb)-p_nh%diag%w_concorr_c(jc,jk+1,jb)) )
            ENDDO
          ENDDO
        ENDIF

        ! Preparations for tracer advection
        IF (lprep_adv .AND. istep == 2) THEN
        IF (lclean_mflx) prep_adv%mass_flx_ic(i_startidx:i_endidx,:,jb) = 0._wp
          DO jk = 1, nlev
!DIR$ IVDEP
            DO jc = i_startidx, i_endidx
              prep_adv%mass_flx_ic(jc,jk,jb) = prep_adv%mass_flx_ic(jc,jk,jb) + r_iadv_rcf*p_nh%diag%rho_ic(jc,jk,jb) * &
                (p_nh%metrics%vwind_expl_wgt(jc,jb)*p_nh%prog(nnow)%w(jc,jk,jb) +                                       &
                 p_nh%metrics%vwind_impl_wgt(jc,jb)*p_nh%prog(nnew)%w(jc,jk,jb) - p_nh%diag%w_concorr_c(jc,jk,jb) )
            ENDDO
          ENDDO
        ENDIF

      ENDDO
!OMP END DO

    ENDIF

!$OMP END PARALLEL

    !-------------------------
    ! communication phase

    IF (timers_level > 5) THEN
      CALL timer_stop(timer_solve_nh_vimpl)
      CALL timer_start(timer_solve_nh_exch)
    ENDIF

    IF (use_icon_comm) THEN 
      IF (istep == 1 .AND. lhdiff_rcf .AND. divdamp_type >= 3) THEN
        CALL icon_comm_sync(p_nh%prog(nnew)%w, z_dwdz_dd, p_patch%sync_cells_not_owned, &
            & name="solve_step1_w")
      ELSE IF (istep == 1) THEN ! Only w is updated in the predictor step
        CALL icon_comm_sync(p_nh%prog(nnew)%w, p_patch%sync_cells_not_owned, &
            & name="solve_step1_w")
      ELSE IF (istep == 2) THEN
        ! Synchronize all prognostic variables
        CALL icon_comm_sync(p_nh%prog(nnew)%rho, p_nh%prog(nnew)%exner, p_nh%prog(nnew)%w, &
          & p_patch%sync_cells_not_owned, name="solve_step2_w")
      ENDIF
    ELSE IF (itype_comm == 1) THEN
      IF (istep == 1) THEN
        IF (lhdiff_rcf .AND. divdamp_type >= 3) THEN
          ! Synchronize w and vertical contribution to divergence damping
          CALL sync_patch_array_mult(SYNC_C,p_patch,2,p_nh%prog(nnew)%w,z_dwdz_dd)
        ELSE 
          ! Only w needs to be synchronized
          CALL sync_patch_array(SYNC_C,p_patch,p_nh%prog(nnew)%w)
        ENDIF
      ELSE ! istep = 2: synchronize all prognostic variables
        CALL sync_patch_array_mult(SYNC_C,p_patch,3,p_nh%prog(nnew)%rho, &
                                   p_nh%prog(nnew)%exner,p_nh%prog(nnew)%w)
      ENDIF
    ENDIF

    IF (timers_level > 5) CALL timer_stop(timer_solve_nh_exch)

    ! end communication phase
    !-------------------------

    ENDDO ! istep-loop

    ! The remaining computations are needed for MPI-parallelized applications only
    IF (my_process_is_mpi_all_seq() ) THEN
      IF (ltimer) CALL timer_stop(timer_solve_nh)
      RETURN
    ENDIF 

! OpenMP directives are commented for the NEC because the overhead is too large
#ifndef __SX__
!$OMP PARALLEL PRIVATE(rl_start,rl_end,i_startblk,i_endblk)
#endif
    IF (l_limited_area .OR. p_patch%id > 1) THEN

      ! Index list over halo points lying in the boundary interpolation zone
      ! Note: this list typically contains at most 10 grid points 
#ifndef __SX__
!$OMP DO PRIVATE(jb,ic,jk,jc) ICON_OMP_DEFAULT_SCHEDULE
#endif
      DO ic = 1, p_nh%metrics%bdy_halo_c_dim

        jb = p_nh%metrics%bdy_halo_c_blk(ic)
        jc = p_nh%metrics%bdy_halo_c_idx(ic)
!DIR$ IVDEP
        DO jk = 1, nlev
          p_nh%prog(nnew)%theta_v(jc,jk,jb) = p_nh%prog(nnew)%exner(jc,jk,jb)

          ! Diagnose exner from rho*theta
          p_nh%prog(nnew)%exner(jc,jk,jb) = EXP(rd_o_cvd*LOG(rd_o_p0ref* &
            p_nh%prog(nnew)%rho(jc,jk,jb)*p_nh%prog(nnew)%theta_v(jc,jk,jb)))

        ENDDO
      ENDDO
#ifndef __SX__
!$OMP END DO
#endif

      rl_start = 1
      rl_end   = grf_bdywidth_c

      i_startblk = p_patch%cells%start_block(rl_start)
      i_endblk   = p_patch%cells%end_block(rl_end)

#ifndef __SX__
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,jc) ICON_OMP_DEFAULT_SCHEDULE
#endif
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, rl_start, rl_end)

        DO jk = 1, nlev
!DIR$ IVDEP
          DO jc = i_startidx, i_endidx

            p_nh%prog(nnew)%theta_v(jc,jk,jb) = p_nh%prog(nnew)%exner(jc,jk,jb)

            ! Diagnose exner from rhotheta
            p_nh%prog(nnew)%exner(jc,jk,jb) = EXP(rd_o_cvd*LOG(rd_o_p0ref* &
              p_nh%prog(nnew)%rho(jc,jk,jb)*p_nh%prog(nnew)%theta_v(jc,jk,jb)))

          ENDDO
        ENDDO
      ENDDO
#ifndef __SX__
!$OMP END DO
#endif
    ENDIF

    rl_start = min_rlcell_int - 1
    rl_end   = min_rlcell

    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

#ifndef __SX__
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,jc) ICON_OMP_DEFAULT_SCHEDULE
#endif
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)

#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        IF (p_nh%metrics%mask_prog_halo_c(jc,jb)) THEN
!DIR$ IVDEP
          DO jk = 1, nlev
#else
      DO jk = 1, nlev
        DO jc = i_startidx, i_endidx
          IF (p_nh%metrics%mask_prog_halo_c(jc,jb)) THEN
#endif
            p_nh%prog(nnew)%theta_v(jc,jk,jb) = p_nh%prog(nnow)%rho(jc,jk,jb)*p_nh%prog(nnow)%theta_v(jc,jk,jb) &
              *( (p_nh%prog(nnew)%exner(jc,jk,jb)/p_nh%prog(nnow)%exner(jc,jk,jb)-1.0_wp) * cvd_o_rd+1.0_wp   ) &
              / p_nh%prog(nnew)%rho(jc,jk,jb)

#ifdef __LOOP_EXCHANGE
          ENDDO
        ENDIF
#else
          ENDIF
        ENDDO
#endif
      ENDDO

    ENDDO
#ifndef __SX__
!$OMP END DO NOWAIT
!$OMP END PARALLEL
#endif

   IF (ltimer) CALL timer_stop(timer_solve_nh)

  END SUBROUTINE solve_nh

END MODULE mo_solve_nonhydro
