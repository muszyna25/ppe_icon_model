!>
!!
!! @brief prepares and postprocesses the fields for and from nwp physics
!!
!! Depending on the action item different sets of physics will be called:
!! this is
!! 1. condensation only so to apply saturation adjustment where needed
!! 2. 'slow physics' means up to now the whole physical package beside
!!     microphysics.
!! 3. turbulence, microphysics and condensation are considered fast physical package
!! 4. Updating the moist tracers in synchrone time intervalls to
!!    advection and saturation adjustment
!!
!!
!!
!! @par Revision History
!!  first implementation by Kristina Froehlich, DWD (2009-06-12)
!!  Call nwp_diagnosis with ih_clch, ih_clcm by Helmut Frank, DWD (2013-01-18)
!!  Calculate gusts in 6 hours               by Helmut Frank, DWD (2013-03-13)
!!
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

! Workaround note: With Cray Fortran 8.5.5, a segmentation fault occurred in
! the SUBROUTINE radheat (source line 1892) when accessing the dummy array
! "pqv".  the workaround here is to copy the "pqv" dummy array to a temporary
! pqv=prm_diag%tot_cld(:,:,jb,iqv)
#if _CRAYFTN == 1 && ( _RELEASE == 8 && _RELEASE_MINOR == 5)
#define __CRAY8_5_5_WORKAROUND
#endif

MODULE mo_nh_interface_nwp

  USE mtime,                      ONLY: datetime
  USE mo_util_mtime,              ONLY: getElapsedSimTimeInSeconds
  USE mo_kind,                    ONLY: wp

  USE mo_timer
  USE mo_exception,               ONLY: message, message_text, finish
  USE mo_impl_constants,          ONLY: itconv, itccov, itrad, itgscp,                        &
    &                                   itsatad, itturb, itsfc, itradheat,                    &
    &                                   itsso, itgwd, itfastphy, icosmo, igme, iedmf,         &
    &                                   min_rlcell_int, min_rledge_int, min_rlcell
  USE mo_impl_constants_grf,      ONLY: grf_bdywidth_c, grf_bdywidth_e
  USE mo_loopindices,             ONLY: get_indices_c, get_indices_e
  USE mo_intp_rbf,                ONLY: rbf_vec_interpol_cell
  USE mo_model_domain,            ONLY: t_patch
  USE mo_intp_data_strc,          ONLY: t_int_state
  USE mo_nonhydro_types,          ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_nonhydrostatic_config,   ONLY: kstart_moist, lhdiff_rcf, ih_clch, ih_clcm, &
    &                                   lcalc_dpsdt
  USE mo_nwp_lnd_types,           ONLY: t_lnd_prog, t_wtr_prog, t_lnd_diag
  USE mo_ext_data_types,          ONLY: t_external_data
  USE mo_nwp_phy_types,           ONLY: t_nwp_phy_diag, t_nwp_phy_tend
  USE mo_parallel_config,         ONLY: nproma, p_test_run, use_icon_comm, use_physics_barrier
  USE mo_diffusion_config,        ONLY: diffusion_config
  USE mo_initicon_config,         ONLY: is_iau_active
  USE mo_run_config,              ONLY: ntracer, iqv, iqc, iqi, iqs, iqtvar, iqtke,  &
    &                                   msg_level, ltimer, timers_level, lart, ldass_lhn
  USE mo_grid_config,             ONLY: l_limited_area
  USE mo_physical_constants,      ONLY: rd, rd_o_cpd, vtmpc1, p0ref, rcvd, cvd, cvv, tmelt, grav

  USE mo_nh_diagnose_pres_temp,   ONLY: diagnose_pres_temp, diag_pres, diag_temp

  USE mo_atm_phy_nwp_config,      ONLY: atm_phy_nwp_config, iprog_aero
  USE mo_util_phys,               ONLY: tracer_add_phytend, iau_update_tracer
  USE mo_lnd_nwp_config,          ONLY: ntiles_total, ntiles_water
  USE mo_cover_koe,               ONLY: cover_koe
  USE mo_satad,                   ONLY: satad_v_3D
  USE mo_aerosol_util,            ONLY: prog_aerosol_2D
  USE mo_radiation,               ONLY: radheat, pre_radiation_nwp
  USE mo_radiation_config,        ONLY: irad_aero
  USE mo_nwp_gw_interface,        ONLY: nwp_gwdrag
  USE mo_nwp_gscp_interface,      ONLY: nwp_microphysics
  USE mo_nwp_turbtrans_interface, ONLY: nwp_turbtrans
  USE mo_nwp_turbdiff_interface,  ONLY: nwp_turbdiff
  USE mo_nwp_sfc_interface,       ONLY: nwp_surface
  USE mo_nwp_conv_interface,      ONLY: nwp_convection
  USE mo_nwp_rad_interface,       ONLY: nwp_radiation
  USE mo_sync,                    ONLY: sync_patch_array, sync_patch_array_mult, SYNC_E,      &
                                        SYNC_C, SYNC_C1
  USE mo_mpi,                     ONLY: my_process_is_mpi_all_parallel, work_mpi_barrier
  USE mo_nwp_diagnosis,           ONLY: nwp_statistics, nwp_diag_output_1, nwp_diag_output_2
  USE mo_icon_comm_lib,           ONLY: new_icon_comm_variable,                               &
    &                                   icon_comm_sync_all, is_ready, until_sync
  USE mo_art_diagnostics_interface,ONLY: art_diagnostics_interface

  USE mo_art_washout_interface,   ONLY: art_washout_interface
  USE mo_art_reaction_interface,  ONLY: art_reaction_interface
  USE mo_linked_list,             ONLY: t_var_list
  USE mo_ls_forcing_nml,          ONLY: is_ls_forcing
  USE mo_ls_forcing,              ONLY: apply_ls_forcing
  USE mo_advection_config,        ONLY: advection_config
  USE mo_o3_util,                 ONLY: calc_o3_gems
  USE mo_edmf_param,              ONLY: edmf_conf
  USE mo_nh_supervise,            ONLY: compute_dpsdt

  USE mo_radar_data_state,        ONLY: radar_data, lhn_fields
  USE mo_latent_heat_nudging,     ONLY: organize_lhn
  USE mo_assimilation_config,     ONLY: assimilation_config

  IMPLICIT NONE

  PRIVATE


  REAL(wp), PARAMETER :: rd_o_p0ref = rd / p0ref
  REAL(wp), PARAMETER :: cpd_o_rd = 1._wp / rd_o_cpd

  PUBLIC :: nwp_nh_interface

CONTAINS
  !
  !-----------------------------------------------------------------------
  !
  SUBROUTINE nwp_nh_interface(lcall_phy_jg, linit, lredgrid,       & !input
                            & dt_loc, dt_phy_jg,                   & !input
                            & mtime_datetime,                      & !input
                            & pt_patch, pt_int_state, p_metrics,   & !input
                            & pt_par_patch,                        & !input
                            & ext_data,                            & !input
                            & pt_prog,                             & !inout
                            & pt_prog_now_rcf, pt_prog_rcf,        & !in/inout
                            & pt_diag ,                            & !inout
                            & prm_diag, prm_nwp_tend, lnd_diag,    & !inout
                            & lnd_prog_now, lnd_prog_new,          & !inout
                            & wtr_prog_now, wtr_prog_new,          & !inout
                            & p_prog_list                          ) !in

    !>
    ! !INPUT PARAMETERS:

    LOGICAL, INTENT(IN)          ::   &             !< physics package time control (switches)
         &                          lcall_phy_jg(:) !< for domain jg
    LOGICAL, INTENT(IN)          :: linit           !< .TRUE. if initialization call (this switch is currently used
                                                    !  to call turbtran in addition to the slow-physics routines
    LOGICAL, INTENT(IN)          :: lredgrid        !< use reduced grid for radiation
    REAL(wp),INTENT(in)          :: dt_loc          !< (advective) time step applicable to local grid level
    REAL(wp),INTENT(in)          :: dt_phy_jg(:)    !< time interval for all physics
                                                    !< packages on domain jg
    TYPE(datetime), POINTER                :: mtime_datetime !< date/time information (in)
    TYPE(t_patch),     TARGET,INTENT(inout):: pt_patch     !<grid/patch info.
    TYPE(t_patch),        TARGET,INTENT(in):: pt_par_patch !<grid/patch info (parent grid)

    TYPE(t_int_state),    TARGET,INTENT(in):: pt_int_state      !< interpolation state
    TYPE(t_nh_metrics)   ,       INTENT(in):: p_metrics
    TYPE(t_external_data),       INTENT(inout):: ext_data

    TYPE(t_nh_diag), TARGET, INTENT(inout) :: pt_diag     !<the diagnostic variables
    TYPE(t_nh_prog), TARGET, INTENT(inout) :: pt_prog     !<the prognostic variables
    TYPE(t_nh_prog), TARGET, INTENT(inout) :: pt_prog_now_rcf !<old state for tke
    TYPE(t_nh_prog), TARGET, INTENT(inout) :: pt_prog_rcf !<the prognostic variables (with
                                                          !< red. calling frequency for tracers!
    TYPE(t_nwp_phy_diag),       INTENT(inout) :: prm_diag
    TYPE(t_nwp_phy_tend),TARGET,INTENT(inout) :: prm_nwp_tend
    TYPE(t_lnd_prog),           INTENT(inout) :: lnd_prog_now, lnd_prog_new
    TYPE(t_wtr_prog),           INTENT(inout) :: wtr_prog_now, wtr_prog_new
    TYPE(t_lnd_diag),           INTENT(inout) :: lnd_diag

    TYPE(t_var_list), INTENT(in) :: p_prog_list !current prognostic state list


    ! !OUTPUT PARAMETERS:            !<variables induced by the whole physics
    ! Local array bounds:

    INTEGER :: nlev, nlevp1            !< number of full levels
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk    !> blocks
    INTEGER :: i_startidx, i_endidx    !< slices

    ! Local scalars:

    INTEGER :: jc,jk,jb,jce      !loop indices
    INTEGER :: jg,jgc            !domain id

    LOGICAL :: ltemp, lpres, ltemp_ifc, l_any_fastphys, l_any_slowphys
    LOGICAL :: lcall_lhn, lcall_lhn_v, lapply_lhn, lcall_lhn_c  !< switches for latent heat nudging
    LOGICAL :: lcompute_tt_lheat                                !< TRUE: store temperature tendency
                                                                ! due to grid scale microphysics 
                                                                ! and satad for latent heat nudging

    INTEGER,  POINTER ::  iidx(:,:,:), iblk(:,:,:)

    REAL(wp), TARGET :: &                                              !> temporal arrays for
      & z_ddt_u_tot (nproma,pt_patch%nlev,pt_patch%nblks_c),&
      & z_ddt_v_tot (nproma,pt_patch%nlev,pt_patch%nblks_c),& !< hor. wind tendencies
      & z_ddt_temp  (nproma,pt_patch%nlev)   !< Temperature tendency

    REAL(wp) :: z_exner_sv(nproma,pt_patch%nlev,pt_patch%nblks_c), z_tempv, sqrt_ri(nproma), n2, dvdz2, &
      zddt_u_raylfric(nproma,pt_patch%nlev), zddt_v_raylfric(nproma,pt_patch%nlev), convfac, wfac

    !< vertical interfaces

    REAL(wp) :: zsct ! solar constant (at time of year)
    REAL(wp) :: zcosmu0 (nproma,pt_patch%nblks_c), cosmu0_slope(nproma,pt_patch%nblks_c)
#ifdef __CRAY8_5_5_WORKAROUND
    REAL(wp) :: pqv(nproma,pt_patch%nlev)
#endif

    REAL(wp) :: z_qsum(nproma,pt_patch%nlev)  !< summand of virtual increment
    REAL(wp) :: z_ddt_alpha(nproma,pt_patch%nlev)  !< tendency of virtual increment

    ! auxiliaries for Rayleigh friction computation
    REAL(wp) :: vabs, rfric_fac, ustart, uoffset_q, ustart_q, max_relax

    ! Variables for EDMF DUALM
    REAL(wp) :: qtvar(nproma,pt_patch%nlev)

    ! communication ids, these do not need to be different variables,
    ! since they are not treated individualy
    INTEGER :: ddt_u_tot_comm, ddt_v_tot_comm, z_ddt_u_tot_comm, z_ddt_v_tot_comm, &
      & tracers_comm, tempv_comm, exner_pr_comm, w_comm

    INTEGER :: ntracer_sync

    ! Pointer to IDs of tracers which contain prognostic condensate.
    ! Required for computing the water loading term 
    INTEGER, POINTER :: condensate_list(:)

    REAL(wp) :: p_sim_time      !< elapsed simulation time on this grid level


    IF (ltimer) CALL timer_start(timer_physics)

    ! calculate elapsed simulation time in seconds (local time for
    ! this domain!)
    p_sim_time = getElapsedSimTimeInSeconds(mtime_datetime) 

    ! local variables related to the blocking

    jg        = pt_patch%id

    IF (pt_patch%n_childdom > 0) THEN
      jgc = pt_patch%child_id(jg)
    ELSE
      jgc = jg
    ENDIF

    ! number of vertical levels
    nlev   = pt_patch%nlev
    nlevp1 = pt_patch%nlevp1

    !define pointers
    iidx  => pt_patch%edges%cell_idx
    iblk  => pt_patch%edges%cell_blk


    IF (lcall_phy_jg(itsatad) .OR. lcall_phy_jg(itgscp) .OR. &
        lcall_phy_jg(itturb)  .OR. lcall_phy_jg(itsfc)) THEN
      l_any_fastphys = .TRUE.
    ELSE
      l_any_fastphys = .FALSE.
    ENDIF


    IF (lcall_phy_jg(itrad) .OR.  lcall_phy_jg(itconv) .OR. lcall_phy_jg(itccov)  &
       .OR. lcall_phy_jg(itsso) .OR. lcall_phy_jg(itgwd)) THEN
      l_any_slowphys = .TRUE.
    ELSE
      l_any_slowphys = .FALSE.
    ENDIF

    ! condensate tracer IDs
    condensate_list => advection_config(jg)%trHydroMass%list


    ! Check whether latent heat nudging is active
    !
    IF (ldass_lhn .AND. .NOT. linit) THEN
      !
      IF ( jg == jgc )  CALL assimilation_config(jg)%dass_g%reinitEvents()
      lcall_lhn   = assimilation_config(jg)%dass_lhn%isActive(mtime_datetime)
      lcall_lhn_v = assimilation_config(jg)%dass_lhn_verif%isActive(mtime_datetime)
      IF (msg_level >= 15) CALL assimilation_config(jg)%dass_g%printStatus(mtime_datetime)
      !
      lcompute_tt_lheat = lcall_lhn .OR. lcall_lhn_v
    ELSE
      lcompute_tt_lheat = .FALSE.
    ENDIF



    IF ( lcall_phy_jg(itturb) .OR. lcall_phy_jg(itconv) .OR.           &
         lcall_phy_jg(itsso)  .OR. lcall_phy_jg(itgwd) .OR. linit ) THEN

      !-------------------------------------------------------------------------
      !>
      !!   Interpolation from v_n onto u,v =>  Reconstruct u and v
      !!   This is needed for turbulence, convection and SSO/GWdrag
      !!
      !-------------------------------------------------------------------------

      IF (msg_level >= 15) &
           & CALL message('mo_nh_interface_nwp:', 'reconstruct u/v')

      IF (timers_level > 3) CALL timer_start(timer_phys_u_v)

      CALL rbf_vec_interpol_cell(pt_prog%vn,            & !< normal wind comp.
        &                        pt_patch,              & !< patch
        &                        pt_int_state,          & !< interpolation state
        &                        pt_diag%u, pt_diag%v,  & !<  reconstr. u,v wind
        &                        opt_rlend=min_rlcell_int )

      IF (timers_level > 3) CALL timer_stop(timer_phys_u_v)

    ENDIF ! diagnose u/v

    IF (msg_level >= 18) THEN

      ! Diagnose temperature needed for debugging output
      CALL diagnose_pres_temp (p_metrics, pt_prog, pt_prog_rcf,    &
           &                              pt_diag, pt_patch,       &
           &                              opt_calc_temp=.TRUE.,    &
           &                              opt_calc_pres=.FALSE.,   &
           &                              opt_rlend=min_rlcell_int )

      ! Write extensive debugging output
      CALL nwp_diag_output_1(pt_patch, pt_diag, pt_prog_rcf)

    ENDIF


    !-------------------------------------------------------------------------
    !>  Update the slow-physics tendencies on the tracer fields,
    !!  afterwards perform saturation adjustment
    !-------------------------------------------------------------------------

    IF (msg_level >= 15) CALL message('mo_nh_interface_nwp:', 'satad')
    IF (timers_level > 2) CALL timer_start(timer_satad_v_3D)

!$OMP PARALLEL PRIVATE (rl_start,rl_end,i_startblk,i_endblk)

    ! Diagnose pressure and temperature in nest boundary zone
    ! (needed for the reduced radiation grid)
    rl_start = 1
    rl_end   = grf_bdywidth_c

    i_startblk = pt_patch%cells%start_block(rl_start)
    i_endblk   = pt_patch%cells%end_block(rl_end)


    IF (jg > 1 .OR. l_limited_area) THEN
!$OMP DO PRIVATE(jb,i_startidx,i_endidx)
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk,  &
                             i_startidx, i_endidx, rl_start, rl_end)

        CALL diag_temp (pt_prog, pt_prog_rcf, condensate_list, pt_diag,    &
                        jb, i_startidx, i_endidx, 1, kstart_moist(jg), nlev)

        CALL diag_pres (pt_prog, pt_diag, p_metrics,     &
                        jb, i_startidx, i_endidx, 1, nlev)

      ENDDO
!$OMP END DO NOWAIT
    ENDIF

    ! Save Exner pressure field on halo points (prognostic points are treated below)
    rl_start = min_rlcell_int-1
    rl_end   = min_rlcell

    i_startblk = pt_patch%cells%start_block(rl_start)
    i_endblk   = pt_patch%cells%end_block(rl_end)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk,  &
                         i_startidx, i_endidx, rl_start, rl_end)

      z_exner_sv(i_startidx:i_endidx,:,jb) = pt_prog%exner(i_startidx:i_endidx,:,jb)
    ENDDO
!$OMP END DO NOWAIT


    ! computations on prognostic points
    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = pt_patch%cells%start_block(rl_start)
    i_endblk   = pt_patch%cells%end_block(rl_end)

!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,z_qsum,z_tempv) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
        & i_startidx, i_endidx, rl_start, rl_end)


      IF (.NOT. linit) THEN

        IF (is_iau_active) THEN
          ! add analysis increments from data assimilation to qv (during IAU phase)
          CALL iau_update_tracer( pt_prog     = pt_prog,     & !in
           &                      p_metrics   = p_metrics,   & !in
           &                      pt_diag     = pt_diag,     & !inout
           &                      pt_prog_rcf = pt_prog_rcf, & !inout tracer
           &                      jg          = jg,          & !in
           &                      jb          = jb,          & !in
           &                      i_startidx  = i_startidx,  & !in
           &                      i_endidx    = i_endidx,    & !in
           &                      kend        = nlev         ) !in
        ENDIF


        ! The provisional "new" tracer state, resulting from the advection 
        ! step, still needs to be updated with the SLOW-physics tracer tendencies 
        ! computed at the end of the last physics call for the then final 
        ! "new" state. The corresponding update for the dynamics variables has 
        ! already happened in the dynamical core.
        !
        CALL tracer_add_phytend( pt_prog      = pt_prog,              & !in density
          &                      prm_nwp_tend = prm_nwp_tend,         & !in
          &                      pdtime       = dt_phy_jg(itfastphy), & !in
          &                      prm_diag     = prm_diag,             & !inout phyfields
          &                      pt_prog_rcf  = pt_prog_rcf,          & !inout tracer
          &                      pt_diag      = pt_diag,              & !inout
          &                      p_metrics    = p_metrics,            & !in
          &                      jg           = jg,                   & !in
          &                      jb           = jb,                   & !in
          &                      i_startidx   = i_startidx,           & !in
          &                      i_endidx     = i_endidx,             & !in
          &                      kend         = nlev                  ) !in

      ENDIF  ! linit


      IF (l_any_fastphys .OR. linit) THEN  ! diagnose temperature
        CALL diag_temp (pt_prog, pt_prog_rcf, condensate_list, pt_diag,    &
                        jb, i_startidx, i_endidx, 1, kstart_moist(jg), nlev)
      ENDIF

      ! Save Exner pressure field (this is needed for a correction to reduce sound-wave generation by latent heating)
      z_exner_sv(i_startidx:i_endidx,:,jb) = pt_prog%exner(i_startidx:i_endidx,:,jb)

      !!-------------------------------------------------------------------------
      !> Initial saturation adjustment (a second one follows at the end of the microphysics)
      !!-------------------------------------------------------------------------

      IF (lcall_phy_jg(itsatad)) THEN

        ! initialize tt_lheat to be in used LHN
        IF (lcompute_tt_lheat) THEN
          prm_diag%tt_lheat (:,:,jb) = - pt_diag%temp   (:,:,jb)
        ENDIF

        IF ( atm_phy_nwp_config(jg)%inwp_turb == iedmf ) THEN   ! EDMF DUALM: no satad in PBL

          CALL satad_v_3D( &
               & maxiter  = 10                             ,& !> IN
               & tol      = 1.e-3_wp                       ,& !> IN
               & te       = pt_diag%temp       (:,:,jb)    ,& !> INOUT
               & qve      = pt_prog_rcf%tracer (:,:,jb,iqv),& !> INOUT
               & qce      = pt_prog_rcf%tracer (:,:,jb,iqc),& !> INOUT
               & rhotot   = pt_prog%rho        (:,:,jb)    ,& !> IN
               & qtvar    = pt_prog_rcf%tracer (:,:,jb,iqtvar) ,& !> IN
               & idim     = nproma                         ,& !> IN
               & kdim     = nlev                           ,& !> IN
               & ilo      = i_startidx                     ,& !> IN
               & iup      = i_endidx                       ,& !> IN
               & klo      = kstart_moist(jg)               ,& !> IN
               & kup      = nlev                            & !> IN
               )
        ELSE

          CALL satad_v_3D( &
               & maxiter  = 10                             ,& !> IN
               & tol      = 1.e-3_wp                       ,& !> IN
               & te       = pt_diag%temp       (:,:,jb)    ,& !> INOUT
               & qve      = pt_prog_rcf%tracer (:,:,jb,iqv),& !> INOUT
               & qce      = pt_prog_rcf%tracer (:,:,jb,iqc),& !> INOUT
               & rhotot   = pt_prog%rho        (:,:,jb)    ,& !> IN
               & idim     = nproma                         ,& !> IN
               & kdim     = nlev                           ,& !> IN
               & ilo      = i_startidx                     ,& !> IN
               & iup      = i_endidx                       ,& !> IN
               & klo      = kstart_moist(jg)               ,& !> IN
               & kup      = nlev                            & !> IN
               )

        ENDIF

        DO jk = kstart_moist(jg), nlev
          DO jc = i_startidx, i_endidx

            ! calculate virtual temperature from condens' output temperature
            ! taken from SUBROUTINE update_tempv_geopot in hydro_atmos/mo_ha_update_diag.f90
            z_qsum(jc,jk) = SUM(pt_prog_rcf%tracer (jc,jk,jb,condensate_list))
          ENDDO
        ENDDO


        DO jk = kstart_moist(jg), nlev
!DIR$ IVDEP
          DO jc = i_startidx, i_endidx

            z_tempv                 = pt_diag%tempv(jc,jk,jb)

            pt_diag%tempv(jc,jk,jb) = pt_diag%temp(jc,jk,jb)                         &
              &                   * ( 1._wp +  vtmpc1                                &
              &                   * pt_prog_rcf%tracer(jc,jk,jb,iqv) - z_qsum(jc,jk) )

            pt_prog%exner(jc,jk,jb) = pt_prog%exner(jc,jk,jb) *        &
              & (1._wp+rd_o_cpd*(pt_diag%tempv(jc,jk,jb)/z_tempv-1._wp))

          ENDDO
        ENDDO

        ! initialize tt_lheat to be in used LHN
        IF (lcompute_tt_lheat) THEN
          prm_diag%tt_lheat (:,:,jb) = prm_diag%tt_lheat (:,:,jb) + pt_diag%temp   (:,:,jb)
        ENDIF
      ENDIF

      IF (lcall_phy_jg(itgscp) .OR. lcall_phy_jg(itturb) .OR. lcall_phy_jg(itsfc)) THEN
        ! diagnose pressure for subsequent fast-physics parameterizations
        CALL diag_pres (pt_prog, pt_diag, p_metrics, jb, i_startidx, i_endidx, 1, nlev)
      ENDIF

      IF (iprog_aero == 1 .AND. .NOT. linit) THEN
        CALL prog_aerosol_2D (nproma,i_startidx,i_endidx,dt_loc,                                         &
                              prm_diag%aerosol(:,:,jb),prm_diag%aercl_ss(:,jb),prm_diag%aercl_or(:,jb),  &
                              prm_diag%aercl_bc(:,jb),prm_diag%aercl_su(:,jb),prm_diag%aercl_du(:,jb),   &
                              prm_diag%dyn_gust(:,jb),prm_diag%con_gust(:,jb),ext_data%atm%soiltyp(:,jb),&
                              ext_data%atm%plcov_t(:,jb,:),ext_data%atm%frac_t(:,jb,:),                  &
                              lnd_prog_now%w_so_t(:,1,jb,:),lnd_prog_now%t_so_t(:,1,jb,:),               &
                              lnd_diag%h_snow_t(:,jb,:)                                                  )
      ENDIF

    ENDDO ! nblks

!$OMP END DO NOWAIT
!$OMP END PARALLEL

    IF (timers_level > 2) CALL timer_stop(timer_satad_v_3D)


    !!-------------------------------------------------------------------------
    !>  turbulent transfer and diffusion and microphysics
    !!
    !!  Because we consider the following physical processes as fast ones
    !!  we allow here the update of prognostic variables inside the subroutines
    !!  This means that the conversion back to the ICON-prognostic variables
    !!  has to be done afterwards
    !!-------------------------------------------------------------------------


    IF ( (lcall_phy_jg(itturb) .OR. linit) .AND. atm_phy_nwp_config(jg)%inwp_turb==iedmf ) THEN

      IF (timers_level > 1) CALL timer_start(timer_nwp_turbulence)

      ! compute turbulent transfer coefficients (atmosphere-surface interface)
      CALL nwp_turbtrans  ( dt_phy_jg(itfastphy),             & !>in
                          & pt_patch, p_metrics,              & !>in
                          & ext_data,                         & !>in
                          & pt_prog_rcf,                      & !>inout
                          & pt_diag,                          & !>inout
                          & prm_diag,                         & !>inout
                          & wtr_prog_now,                     & !>in
                          & lnd_prog_now,                     & !>inout
                          & lnd_diag                          ) !>inout

      IF (timers_level > 1) CALL timer_stop(timer_nwp_turbulence)
    ENDIF !lcall(itturb)


    !For turbulence schemes NOT including the call to the surface scheme.
    !nwp_surface must even be called in inwp_surface = 0 because the
    !the lower boundary conditions for the turbulence scheme
    !are not set otherwise

    IF ( l_any_fastphys .AND. ( ANY( (/icosmo,igme/)==atm_phy_nwp_config(jg)%inwp_turb ) &
                  & .OR. ( edmf_conf==2  .AND. iedmf==atm_phy_nwp_config(jg)%inwp_turb ) ) ) THEN
      IF (timers_level > 2) CALL timer_start(timer_nwp_surface)

       !> as pressure is needed only for an approximate adiabatic extrapolation
       !! of the temperature at the lowest model level towards ground level,
       !! a recalculation is not required
       CALL nwp_surface    (  dt_phy_jg(itfastphy),              & !>input
                             & pt_patch,                         & !>input
                             & ext_data,                         & !>input
                             & pt_prog_rcf,                      & !>in/inout rcf=reduced calling freq.
                             & pt_diag ,                         & !>inout
                             & prm_diag,                         & !>inout
                             & lnd_prog_now, lnd_prog_new,       & !>inout
                             & wtr_prog_now, wtr_prog_new,       & !>inout
                             & lnd_diag                          ) !>input

      IF (timers_level > 2) CALL timer_stop(timer_nwp_surface)
    END IF


    !Call to turbulent parameterization schemes
    IF (  lcall_phy_jg(itturb) ) THEN

      IF (timers_level > 1) CALL timer_start(timer_nwp_turbulence)

      SELECT CASE (atm_phy_nwp_config(jg)%inwp_turb)

      !Turbulence schemes NOT including the call to the surface scheme
      CASE(icosmo,igme,iedmf)

        ! compute turbulent diffusion (atmospheric column)
        CALL nwp_turbdiff   (  dt_phy_jg(itfastphy),              & !>in
                              & pt_patch, p_metrics,              & !>in
                              & ext_data,                         & !>in
                              & pt_prog,                          & !>in
                              & pt_prog_now_rcf, pt_prog_rcf,     & !>in/inout
                              & pt_diag,                          & !>inout
                              & prm_diag, prm_nwp_tend,           & !>inout
                              & wtr_prog_now,                     & !>in
                              & lnd_prog_now,                     & !>in
                              & lnd_diag                          ) !>in


      CASE DEFAULT

        CALL finish('mo_nh_interface_nwp:','unknown choice of turbulence scheme')

      END SELECT

      IF (timers_level > 1) CALL timer_stop(timer_nwp_turbulence)

    END IF


    !-------------------------------------------------------------------------
    !  prognostic microphysic and precipitation scheme
    !-------------------------------------------------------------------------
    IF ( lcall_phy_jg(itgscp)) THEN

      IF (msg_level >= 15) &
        & CALL message('mo_nh_interface_nwp:', 'microphysics')

      !> temperature and tracers have been updated by turbulence;
      !! an update of the pressure field is not needed because pressure
      !! is not needed at high accuracy in the microphysics scheme
      !! note: after the microphysics the second call to SATAD is within 
      !!       the nwp_microphysics routine (first one is above)

      IF (timers_level > 1) CALL timer_start(timer_nwp_microphysics)

      CALL nwp_microphysics ( dt_phy_jg(itfastphy),             & !>input
                            & lcall_phy_jg(itsatad),            & !>input
                            & pt_patch, p_metrics,              & !>input
                            & pt_prog,                          & !>inout
                            & pt_prog_rcf,                      & !>inout
                            & pt_diag ,                         & !>inout
                            & prm_diag, prm_nwp_tend,           & !>inout
                            & lcompute_tt_lheat                 ) !>in

      IF (timers_level > 1) CALL timer_stop(timer_nwp_microphysics)

    ENDIF

    IF (lart) THEN
      CALL calc_o3_gems(pt_patch,mtime_datetime,pt_diag,prm_diag,ext_data)

      IF (.NOT. linit) THEN
        CALL art_reaction_interface(ext_data,              & !> in
                &                   pt_patch,              & !> in
                &                   mtime_datetime,        & !> in
                &                   dt_phy_jg(itfastphy),  & !> in
                &                   p_prog_list,           & !> in
                &                   pt_prog,               & !> in
                &                   p_metrics,             & !> in
                &                   pt_diag,               & !> inout
                &                   pt_prog_rcf%tracer,    & !>
                &                   prm_diag = prm_diag)     !> optional
                
      END IF

      CALL art_washout_interface(pt_prog,pt_diag,              & !>in
                &          dt_phy_jg(itfastphy),               & !>in
                &          pt_patch,                           & !>in
                &          prm_diag,                           & !>in
                &          p_metrics,                          & !>in
                &          pt_prog_rcf%tracer)                   !>inout
    ENDIF !lart


    !!------------------------------------------------------------------
    !> Latent heat nudging (optional)
    !!------------------------------------------------------------------
    IF (ldass_lhn .AND. .NOT. linit) THEN

      IF (msg_level >= 15) CALL message('mo_nh_interface_nwp:', 'applying LHN')
      IF (timers_level > 1) CALL timer_start(timer_datass)

      IF (lcall_lhn .OR. lcall_lhn_v) THEN
         CALL organize_lhn (   &
                               & dt_loc,                           & !>input
                               & p_sim_time,                       & ! in
                               & pt_patch, p_metrics,              & !>input
                               & pt_int_state,                     & !>input
                               & pt_prog_rcf,                      & !>inout
                               & pt_diag ,                         & !>inout
                               & prm_diag,                         & !>inout
                               & lhn_fields(jg),                   & !>inout
                               & radar_data(jg),                   & 
                               & prm_nwp_tend,                     &
                               & mtime_datetime,                   &
                               & lcall_lhn, lcall_lhn_v            )
      ELSE IF (msg_level >= 15) THEN
        WRITE (message_text,'(a,f10.2,i5)') 'LHN not running because of specified times!',p_sim_time,jg
        CALL message('mo_nh_interface_nwp:', message_text)
      ENDIF


      lcall_lhn_c = assimilation_config(jgc)%dass_lhn%isActive(mtime_datetime)
      lapply_lhn  = lcall_lhn .OR. lcall_lhn_c

      IF (lapply_lhn) THEN

        rl_start = grf_bdywidth_c+1
        rl_end   = min_rlcell_int
  
        i_startblk = pt_patch%cells%start_block(rl_start)
        i_endblk   = pt_patch%cells%end_block(rl_end)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = i_startblk, i_endblk
  
          CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
            & i_startidx, i_endidx, rl_start, rl_end )
  
          ! update prognostic variables
          DO jk = 1, nlev
            DO jc = i_startidx, i_endidx
              pt_diag%temp(jc,jk,jb) = pt_diag%temp(jc,jk,jb) + lhn_fields(jg)%ttend_lhn(jc,jk,jb) * dt_loc
              pt_prog_rcf%tracer(jc,jk,jb,iqv) = pt_prog_rcf%tracer(jc,jk,jb,iqv) + lhn_fields(jg)%qvtend_lhn(jc,jk,jb) * dt_loc
            ENDDO
          ENDDO
          !-------------------------------------------------------------------------
          !   call the saturation adjustment
          !-------------------------------------------------------------------------
  
          CALL satad_v_3D( &
               & maxiter  = 10                             ,& !> IN
               & tol      = 1.e-3_wp                       ,& !> IN
               & te       = pt_diag%temp       (:,:,jb)    ,& !> INOUT
               & qve      = pt_prog_rcf%tracer (:,:,jb,iqv),& !> INOUT
               & qce      = pt_prog_rcf%tracer (:,:,jb,iqc),& !> INOUT
               & rhotot   = pt_prog%rho        (:,:,jb)    ,& !> IN
               & idim     = nproma                         ,& !> IN
               & kdim     = nlev                           ,& !> IN
               & ilo      = i_startidx                     ,& !> IN
               & iup      = i_endidx                       ,& !> IN
               & klo      = kstart_moist(jg)               ,& !> IN
               & kup      = nlev                            ) !> IN
  
        ENDDO ! nblks
  
!$OMP END DO NOWAIT
!$OMP END PARALLEL

      ENDIF

      IF (timers_level > 1) CALL timer_stop(timer_datass)

    ENDIF ! ldass_lhn


    IF (timers_level > 1) CALL timer_start(timer_fast_phys)

    ! Remark: in the (unusual) case that satad is used without any other physics,
    ! recalculation of the thermodynamic variables is duplicated here. However,
    ! this is the easiest way to combine minimization of halo communications
    ! with a failsafe flow control

    IF (msg_level >= 15) &
      & CALL message('mo_nh_interface_nwp:', 'recalculate thermodynamic variables')


    ! exclude boundary interpolation zone of nested domains
    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = pt_patch%cells%start_block(rl_start)
    i_endblk   = pt_patch%cells%end_block(rl_end)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx, i_endidx, z_qsum) ICON_OMP_DEFAULT_SCHEDULE

    DO jb = i_startblk, i_endblk
      CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
        & i_startidx, i_endidx, rl_start, rl_end )

      IF (lcall_phy_jg(itsatad) .OR. lcall_phy_jg(itgscp) .OR. lcall_phy_jg(itturb)) THEN
        !-------------------------------------------------------------------------
        !>
        !! re-calculate scalar prognostic variables out of physics variables!
        !!
        !-------------------------------------------------------------------------

        IF (kstart_moist(jg) > 1) z_qsum(:,1:kstart_moist(jg)-1) = 0._wp

        DO jk = kstart_moist(jg), nlev
          DO jc = i_startidx, i_endidx

            z_qsum(jc,jk) = SUM(pt_prog_rcf%tracer (jc,jk,jb,condensate_list))

          ENDDO
        ENDDO


        DO jk = 1, nlev
!DIR$ IVDEP
          DO jc =  i_startidx, i_endidx

            pt_diag%tempv(jc,jk,jb) =  pt_diag%temp(jc,jk,jb)          &
&                                  * ( 1._wp +  vtmpc1                 &
&                                  *  pt_prog_rcf%tracer(jc,jk,jb,iqv) &
&                                   - z_qsum(jc,jk) )

            pt_prog%exner(jc,jk,jb) = EXP(rd_o_cpd*LOG(rd_o_p0ref                   &
              &                     * pt_prog%rho(jc,jk,jb)*pt_diag%tempv(jc,jk,jb)))

            pt_diag%exner_pr(jc,jk,jb) = pt_diag%exner_pr(jc,jk,jb) + &
              pt_prog%exner(jc,jk,jb) - z_exner_sv(jc,jk,jb)

            pt_prog%theta_v(jc,jk,  jb) = pt_diag%tempv(jc,jk,jb) &
&                                       / pt_prog%exner(jc,jk,jb)

          ENDDO
        ENDDO

        ! compute dynamical temperature tendency from increments of Exner function and density
        ! the virtual increment is neglected here because this tendency is used only as
        ! input for the convection scheme, which is rather insensitive against this quantity
        IF ( lcall_phy_jg(itconv) ) THEN
          DO jk = kstart_moist(jg), nlev
!DIR$ IVDEP
            DO jc =  i_startidx, i_endidx

              pt_diag%ddt_temp_dyn(jc,jk,jb) = pt_diag%temp(jc,jk,jb)/dt_phy_jg(itfastphy) * &
                ( cpd_o_rd/pt_prog%exner(jc,jk,jb)*pt_diag%exner_dyn_incr(jc,jk,jb) -        &
                ( pt_diag%airmass_new(jc,jk,jb)-pt_diag%airmass_now(jc,jk,jb) ) /            &
                pt_diag%airmass_new(jc,jk,jb) )

            ENDDO
          ENDDO
        ENDIF

      ENDIF ! recalculation

      IF (lcall_phy_jg(itturb) .OR. linit .OR. l_any_slowphys) THEN
        ! rediagnose pressure
        CALL diag_pres (pt_prog, pt_diag, p_metrics,     &
                        jb, i_startidx, i_endidx, 1, nlev)
      ENDIF

    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    IF (timers_level > 1) CALL timer_stop(timer_fast_phys)

    IF ( (lcall_phy_jg(itturb) .OR. linit) .AND. ANY( (/icosmo,igme/)==atm_phy_nwp_config(jg)%inwp_turb ) ) THEN

      IF (timers_level > 1) CALL timer_start(timer_nwp_turbulence)

      ! compute turbulent transfer coefficients (atmosphere-surface interface)
      CALL nwp_turbtrans  ( dt_phy_jg(itfastphy),             & !>in
                          & pt_patch, p_metrics,              & !>in
                          & ext_data,                         & !>in
                          & pt_prog_rcf,                      & !>inout
                          & pt_diag,                          & !>inout
                          & prm_diag,                         & !>inout
                          & wtr_prog_new,                     & !>in
                          & lnd_prog_new,                     & !>inout
                          & lnd_diag                          ) !>inout

      IF (timers_level > 1) CALL timer_stop(timer_nwp_turbulence)
    ENDIF !lcall(itturb)



    !!-------------------------------------------------------------------------
    !!  slow physics part
    !!-------------------------------------------------------------------------

    IF (l_any_slowphys) THEN

      IF (msg_level >= 15) &
         CALL message('mo_nh_interface', 'diagnose pres/temp for slow physics')

      ! If slow physics is called without fast physics (which should happen
      ! at the initial time step only), temperature needs to be calculated
      ! Otherwise, temperature is up to date
      IF ( .NOT. (lcall_phy_jg(itgscp) .OR. lcall_phy_jg(itturb))) THEN
        ltemp = .TRUE.
      ELSE
        ltemp = .FALSE.
      ENDIF


      ! Pressure has already been updated at the end of the fast physics part
      lpres = .FALSE.

      ! Temperature at interface levels is needed if irad_aero = 5 or 6
      ! or if Ritter-Geleyn radiation is called
      IF ( lcall_phy_jg(itrad) .AND. ( irad_aero == 5 .OR. irad_aero == 6 &
           .OR. irad_aero == 9 .OR. atm_phy_nwp_config(jg)%inwp_radiation == 2 ) ) THEN
        ltemp_ifc = .TRUE.
      ELSE
        ltemp_ifc = .FALSE.
      ENDIF

      CALL diagnose_pres_temp (p_metrics, pt_prog, pt_prog_rcf,   &
        &                      pt_diag, pt_patch,                 &
        &                      opt_calc_temp     = ltemp,         &
        &                      opt_calc_pres     = lpres,         &
        &                      lnd_prog          = lnd_prog_new,  &
        &                      opt_calc_temp_ifc = ltemp_ifc,     &
        &                      opt_rlend         = min_rlcell_int )

    ENDIF


    !-------------------------------------------------------------------------
    !> Convection
    !-------------------------------------------------------------------------

    IF ( lcall_phy_jg(itconv)  ) THEN


      IF (msg_level >= 15) &
&           CALL message('mo_nh_interface', 'convection')

      IF (timers_level > 2) CALL timer_start(timer_nwp_convection)
      CALL nwp_convection (  dt_phy_jg(itconv),                 & !>input
                            & pt_patch, p_metrics,              & !>input
                            & ext_data,                         & !>input
                            & pt_prog,                          & !>input
                            & pt_prog_rcf,                      & !>input
                            & pt_diag,                          & !>inout
                            & prm_diag, prm_nwp_tend            ) !>inout
      IF (timers_level > 2) CALL timer_stop(timer_nwp_convection)

    ENDIF! convection


    !-------------------------------------------------------------------------
    !> Cloud cover
    !-------------------------------------------------------------------------

    IF ( lcall_phy_jg(itccov) ) THEN

      ! When using a reduced grid for radiation, part of the boundary points need
      ! to be included in order to compute spatial gradients for back-interpolation
      IF (lredgrid) THEN
        rl_start = grf_bdywidth_c-1
      ELSE
        rl_start = grf_bdywidth_c+1
      ENDIF
      rl_end   = min_rlcell_int

      i_startblk = pt_patch%cells%start_block(rl_start)
      i_endblk   = pt_patch%cells%end_block(rl_end)

      IF (msg_level >= 15) &
        &  CALL message('mo_nh_interface', 'cloud cover')

      IF (timers_level > 2) CALL timer_start(timer_cover_koe)


      !-------------------------------------------------------------------------
      !> Cloud water distribution: cloud cover, cloud water, cloud ice
      !  inwp_cldcover =
      !  (0) no clouds
      !  (1) diagnostic cloud cover
      !  (2) prognostic total water variance (not yet started)
      !  (3) clouds as in COSMO
      !  (4) clouds as in turbulence
      !  (5) grid-scale cloud cover [1 or 0]
      !-------------------------------------------------------------------------

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,qtvar) ICON_OMP_GUIDED_SCHEDULE
      DO jb = i_startblk, i_endblk
        !
        CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
&                       i_startidx, i_endidx, rl_start, rl_end)

        IF ( atm_phy_nwp_config(jg)%inwp_turb == iedmf ) THEN
          qtvar(:,:) = pt_prog_rcf%tracer(:,:,jb,iqtvar)        ! EDMF DUALM
        ELSE
          qtvar(:,:) = 0.0_wp                                   ! other turb schemes
        ENDIF

        CALL cover_koe &
&             (kidia  = i_startidx ,   kfdia  = i_endidx  ,       & !! in:  horizonal begin, end indices
&              klon = nproma,  kstart = kstart_moist(jg)  ,       & !! in:  horiz. and vert. vector length
&              klev   = nlev,                                     &
&              icldscheme = atm_phy_nwp_config(jg)%inwp_cldcover ,& !! in:  cloud cover option
&              inwp_turb  = atm_phy_nwp_config(jg)%inwp_turb,     & !! in:  turbulence scheme number
&              tt     = pt_diag%temp         (:,:,jb)     ,       & !! in:  temperature at full levels
&              pp     = pt_diag%pres         (:,:,jb)     ,       & !! in:  pressure at full levels
&              ps     = pt_diag%pres_sfc     (:,jb)       ,       & !! in:  surface pressure at full levels
&              t_g    = lnd_prog_new%t_g     (:,jb)       ,       & !! in:  surface temperature
&              pgeo   = p_metrics%geopot_agl (:,:,jb)     ,       & !! in:  geopotential height
&              rho    = pt_prog%rho          (:,:,jb  )   ,       & !! in:  density
&              rcld   = prm_diag%rcld        (:,:,jb)     ,       & !! in:  standard deviation of saturation deficit
&              ldland = ext_data%atm%llsm_atm_c (:,jb)    ,       & !! in:  land/sea mask
&              ldcum  = prm_diag%locum       (:,jb)       ,       & !! in:  convection on/off
&              kcbot  = prm_diag%mbas_con    (:,jb)       ,       & !! in:  convective cloud base
&              kctop  = prm_diag%mtop_con    (:,jb)       ,       & !! in:  convective cloud top
&              ktype  = prm_diag%ktype       (:,jb)       ,       & !! in:  convection type
&              pmfude_rate = prm_diag%con_udd(:,:,jb,3)   ,       & !! in:  convective updraft detrainment rate
&              plu         = prm_diag%con_udd(:,:,jb,7)   ,       & !! in:  updraft condensate
&              qc_tend= prm_nwp_tend%ddt_tracer_pconv(:,:,jb,iqc),& !! in:  convective qc tendency
&              qv     = pt_prog_rcf%tracer   (:,:,jb,iqv) ,       & !! in:  spec. humidity
&              qc     = pt_prog_rcf%tracer   (:,:,jb,iqc) ,       & !! in:  cloud water
&              qi     = pt_prog_rcf%tracer   (:,:,jb,iqi) ,       & !! in:  cloud ice
&              qs     = pt_prog_rcf%tracer   (:,:,jb,iqs) ,       & !! in:  snow
&              qtvar  = qtvar                             ,       & !! in:  qtvar
&              cc_tot = prm_diag%clc         (:,:,jb)     ,       & !! out: cloud cover
&              qv_tot = prm_diag%tot_cld     (:,:,jb,iqv) ,       & !! out: qv       -"-
&              qc_tot = prm_diag%tot_cld     (:,:,jb,iqc) ,       & !! out: clw      -"-
&              qi_tot = prm_diag%tot_cld     (:,:,jb,iqi) )         !! out: ci       -"-


      ENDDO

!$OMP END DO NOWAIT
!$OMP END PARALLEL

      IF (timers_level > 2) CALL timer_stop(timer_cover_koe)

    ENDIF! cloud cover

    !-------------------------------------------------------------------------
    !> Radiation
    !-------------------------------------------------------------------------

    IF ( lcall_phy_jg(itrad) ) THEN

      IF (ltimer) CALL timer_start(timer_nwp_radiation)
      CALL nwp_radiation (lredgrid,              & ! in
           &              p_sim_time,            & ! in
           &              mtime_datetime,        & ! in
           &              pt_patch,              & ! in
           &              pt_par_patch,          & ! in
           &              ext_data,              & ! in
           &              lnd_diag,              & ! in
           &              pt_prog,               & ! inout
           &              pt_diag,               & ! inout
           &              prm_diag,              & ! inout
           &              lnd_prog_new,          & ! in
           &              wtr_prog_new           ) ! in
      IF (ltimer) CALL timer_stop(timer_nwp_radiation)

    ENDIF


    IF ( lcall_phy_jg(itradheat) ) THEN

      IF (msg_level >= 15) &
&           CALL message('mo_nh_interface', 'radiative heating')


      IF (timers_level > 10) CALL timer_start(timer_pre_radiation_nwp)

      CALL pre_radiation_nwp (                      &
        & kbdim      = nproma,                      &
        & p_inc_rad  = dt_phy_jg(itfastphy),        &
        & p_sim_time = p_sim_time,                  &
        & pt_patch   = pt_patch,                    &
        & zsmu0      = zcosmu0,                     &
        & zsct       = zsct,                        &
        & slope_ang  = p_metrics%slope_angle,       &
        & slope_azi  = p_metrics%slope_azimuth,     &
        & cosmu0_slp = cosmu0_slope                 )

      IF (timers_level > 10) CALL timer_stop(timer_pre_radiation_nwp)

      ! exclude boundary interpolation zone of nested domains
      rl_start = grf_bdywidth_c+1
      rl_end   = min_rlcell_int

      i_startblk = pt_patch%cells%start_block(rl_start)
      i_endblk   = pt_patch%cells%end_block(rl_end)

      IF (timers_level > 2) CALL timer_start(timer_radheat)
#ifdef __CRAY8_5_5_WORKAROUND
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,pqv) ICON_OMP_DEFAULT_SCHEDULE
#else
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
#endif
!
      DO jb = i_startblk, i_endblk
        !
        CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
&                       i_startidx, i_endidx, rl_start, rl_end)

        zcosmu0 (i_startidx:i_endidx,jb) &
          = 0.5_wp * (ABS(zcosmu0(i_startidx:i_endidx,jb)) &
          &           + zcosmu0(i_startidx:i_endidx,jb))

        !calculate solar incoming flux at TOA
        prm_diag%flxdwswtoa(i_startidx:i_endidx,jb) = zcosmu0(i_startidx:i_endidx,jb) &
          &                                         * zsct                 !zsct by pre_radiation

        prm_diag%swflxsfc (:,jb)=0._wp
        prm_diag%lwflxsfc (:,jb)=0._wp
        prm_diag%swflxtoa (:,jb)=0._wp

#ifdef __CRAY8_5_5_WORKAROUND
        ! workaround for Cray Fortran 8.5.5
        pqv=prm_diag%tot_cld(:,:,jb,iqv)
#endif

        IF (atm_phy_nwp_config(jg)%inwp_surface >= 1) THEN

          prm_diag%swflxsfc_t (:,jb,:)=0._wp
          prm_diag%lwflxsfc_t (:,jb,:)=0._wp

          CALL radheat (                   &
          !
          ! input
          ! -----
          !
          & jcs=i_startidx                         ,&! in     start index of inner do loop
          & jce=i_endidx                           ,&! in     end index of inner do loop
          & kbdim=nproma                           ,&! in     loop length and dimension size
          & klev=nlev                              ,&! in     vertical dimension size
          & klevp1=nlevp1                          ,&! in     vertical dimension size
          & ntiles=ntiles_total                    ,&! in     number of tiles of sfc flux fields
          & ntiles_wtr=ntiles_water                ,&! in     number of extra tiles for ocean and lakes
          & pmair=pt_diag%airmass_new(:,:,jb)      ,&! in     layer air mass             [kg/m2]
#ifdef __CRAY8_5_5_WORKAROUND
          & pqv=pqv                                ,&! in     specific moisture           [kg/kg]
#else
          & pqv=prm_diag%tot_cld(:,:,jb,iqv)       ,&! in     specific moisture           [kg/kg]
#endif

          & pcd=cvd                                ,&! in     specific heat of dry air  [J/kg/K]
          & pcv=cvv                                ,&! in     specific heat of vapor    [J/kg/K]
          & pi0=prm_diag%flxdwswtoa(:,jb)          ,&! in     solar incoming flux at TOA  [W/m2]
          & pemiss=ext_data%atm%emis_rad(:,jb)     ,&! in     lw sfc emissivity
          & pqc=prm_diag%tot_cld    (:,:,jb,iqc)   ,&! in     specific cloud water        [kg/kg]
          & pqi=prm_diag%tot_cld    (:,:,jb,iqi)   ,&! in     specific cloud ice          [kg/kg]
          & ppres_ifc=pt_diag%pres_ifc(:,:,jb)     ,&! in     pressure at layer boundaries [Pa]
          & albedo=prm_diag%albdif(:,jb)           ,&! in     grid-box average shortwave albedo
          & albedo_t=prm_diag%albdif_t(:,jb,:)     ,&! in     tile-specific shortwave albedo
          & lp_count=ext_data%atm%lp_count(jb)     ,&! in     number of land points
          & gp_count_t=ext_data%atm%gp_count_t(jb,:),&! in    number of land points per tile
          & spi_count =ext_data%atm%spi_count(jb)  ,&! in     number of seaice points
          & fp_count  =ext_data%atm%fp_count(jb)   ,&! in     number of (f)lake points
          & idx_lst_lp=ext_data%atm%idx_lst_lp(:,jb), &! in   index list of land points
          & idx_lst_t=ext_data%atm%idx_lst_t(:,jb,:), &! in   index list of land points per tile
          & idx_lst_spi=ext_data%atm%idx_lst_spi(:,jb),&! in  index list of seaice points
          & idx_lst_fp=ext_data%atm%idx_lst_fp(:,jb),&! in    index list of (f)lake points
          & cosmu0=zcosmu0(:,jb)                   ,&! in     cosine of solar zenith angle (w.r.t. plain surface)
          & cosmu0_slp=cosmu0_slope(:,jb)          ,&! in     slope-dependent cosine of solar zenith angle
          & opt_nh_corr=.TRUE.                     ,&! in     switch for NH mode
          & ptsfc=lnd_prog_new%t_g(:,jb)           ,&! in     surface temperature         [K]
          & ptsfc_t=lnd_prog_new%t_g_t(:,jb,:)     ,&! in     tile-specific surface temperature         [K]
          & ptsfctrad=prm_diag%tsfctrad(:,jb)      ,&! in     sfc temp. used for pflxlw   [K]
          & ptrmsw=prm_diag%trsolall (:,:,jb)      ,&! in     shortwave net tranmissivity []
          & pflxlw=prm_diag%lwflxall (:,:,jb)      ,&! in     longwave net flux           [W/m2]
          & lwflx_up_sfc_rs=prm_diag%lwflx_up_sfc_rs(:,jb), &! in longwave upward flux at surface [W/m2]
          & trsol_up_toa=prm_diag%trsol_up_toa(:,jb),   & ! in shortwave upward transm. at the top of the atmosphere
          & trsol_up_sfc=prm_diag%trsol_up_sfc(:,jb),   & ! in shortwave upward transm. at the surface
          & trsol_par_sfc=prm_diag%trsol_par_sfc(:,jb), & ! in photosynthetically active downward transm. at the surface
          & trsol_dn_sfc_diff=prm_diag%trsol_dn_sfc_diff(:,jb),&! in shortwave diffuse downward transm. at the surface
          & trsol_clr_sfc=prm_diag%trsolclr_sfc(:,jb),  & ! in clear-sky net transmissivity at surface
          & use_trsolclr_sfc=atm_phy_nwp_config(jg)%inwp_radiation/=2, &
          !
          ! output
          ! ------
          !
          & pdtdtradsw=prm_nwp_tend%ddt_temp_radsw(:,:,jb),&! out    rad. heating by SW        [K/s]
          & pdtdtradlw=prm_nwp_tend%ddt_temp_radlw(:,:,jb),&! out    rad. heating by lw        [K/s]
          & pflxsfcsw =prm_diag%swflxsfc (:,jb)   ,&        ! out shortwave surface net flux [W/m2]
          & pflxsfclw =prm_diag%lwflxsfc (:,jb)   ,&        ! out longwave surface net flux  [W/m2]
          & pflxsfcsw_t=prm_diag%swflxsfc_t (:,jb,:)   ,&   ! out tile-specific shortwave surface net flux [W/m2]
          & pflxsfclw_t=prm_diag%lwflxsfc_t (:,jb,:)   ,&   ! out tile-specific longwave surface net flux  [W/m2]
          & pflxtoasw =prm_diag%swflxtoa (:,jb)        ,&   ! out shortwave toa net flux     [W/m2]
          & lwflx_up_sfc=prm_diag%lwflx_up_sfc(:,jb)   ,&   ! out longwave upward flux at surface [W/m2]
          & swflx_up_toa=prm_diag%swflx_up_toa(:,jb)   ,&   ! out shortwave upward flux at the TOA [W/m2]
          & swflx_up_sfc=prm_diag%swflx_up_sfc(:,jb)   ,&   ! out shortwave upward flux at the surface [W/m2]
          & swflx_par_sfc=prm_diag%swflx_par_sfc(:,jb) ,&   ! out shortwave upward flux at the surface [W/m2]
          & swflx_dn_sfc_diff=prm_diag%swflx_dn_sfc_diff(:,jb) ) ! out shortwave diffuse downward flux at the surface [W/m2]

        ELSE
          CALL radheat (                   &
          !
          ! input
          ! -----
          !
          & jcs=i_startidx                         ,&! in     start index of inner do loop
          & jce=i_endidx                           ,&! in     end index of inner do loop
          & kbdim=nproma                           ,&! in     loop length and dimension size
          & klev=nlev                              ,&! in     vertical dimension size
          & klevp1=nlevp1                          ,&! in     vertical dimension size
          & ntiles=1                               ,&! in     number of tiles of sfc flux fields
          & ntiles_wtr=0                           ,&! in     number of extra tiles for ocean and lakes
          & pmair=pt_diag%airmass_new(:,:,jb)      ,&! in     layer air mass             [kg/m2]
#ifdef __CRAY8_5_5_WORKAROUND
          & pqv=pqv                                ,&! in     specific moisture           [kg/kg]
#else
          & pqv=prm_diag%tot_cld(:,:,jb,iqv)       ,&! in     specific moisture           [kg/kg]
#endif
          & pcd=cvd                                ,&! in     specific heat of dry air  [J/kg/K]
          & pcv=cvv                                ,&! in     specific heat of vapor    [J/kg/K]
          & pi0=prm_diag%flxdwswtoa(:,jb)          ,&! in     solar incoming flux at TOA  [W/m2]
          & pemiss=ext_data%atm%emis_rad(:,jb)     ,&! in     lw sfc emissivity
          & pqc=prm_diag%tot_cld    (:,:,jb,iqc)   ,&! in     specific cloud water        [kg/kg]
          & pqi=prm_diag%tot_cld    (:,:,jb,iqi)   ,&! in     specific cloud ice          [kg/kg]
          & ppres_ifc=pt_diag%pres_ifc(:,:,jb)     ,&! in     pressure at layer boundaries [Pa]
          & cosmu0=zcosmu0(:,jb)                   ,&! in     cosine of solar zenith angle
          & cosmu0_slp=cosmu0_slope(:,jb)          ,&! in     slope-dependent cosine of solar zenith angle
          & opt_nh_corr=.TRUE.                     ,&! in     switch for NH mode
          & ptsfc=lnd_prog_new%t_g(:,jb)           ,&! in     surface temperature         [K]
          & ptsfctrad=prm_diag%tsfctrad(:,jb)      ,&! in     sfc temp. used for pflxlw   [K]
          & ptrmsw=prm_diag%trsolall (:,:,jb)      ,&! in     shortwave net tranmissivity []
          & pflxlw=prm_diag%lwflxall (:,:,jb)      ,&! in     longwave net flux           [W/m2]
          & lwflx_up_sfc_rs=prm_diag%lwflx_up_sfc_rs(:,jb), &! in longwave upward flux at surface [W/m2]
          & trsol_up_toa=prm_diag%trsol_up_toa(:,jb),   & ! in shortwave upward transm. at the top of the atmosphere
          & trsol_up_sfc=prm_diag%trsol_up_sfc(:,jb),   & ! in shortwave upward transm. at the surface
          & trsol_par_sfc=prm_diag%trsol_par_sfc(:,jb), & ! in photosynthetically active downward transm. at the surface
          & trsol_dn_sfc_diff=prm_diag%trsol_dn_sfc_diff(:,jb),&! in shortwave diffuse downward transm. at the surface
          !
          ! output
          ! ------
          !
          & pdtdtradsw=prm_nwp_tend%ddt_temp_radsw(:,:,jb),&! out    rad. heating by SW        [K/s]
          & pdtdtradlw=prm_nwp_tend%ddt_temp_radlw(:,:,jb),&! out    rad. heating by lw        [K/s]
          & pflxsfcsw =prm_diag%swflxsfc (:,jb)   ,&        ! out shortwave surface net flux [W/m2]
          & pflxsfclw =prm_diag%lwflxsfc (:,jb)   ,&        ! out longwave surface net flux  [W/m2]
          & pflxtoasw =prm_diag%swflxtoa (:,jb)   ,&        ! out shortwave toa net flux     [W/m2]
          & lwflx_up_sfc=prm_diag%lwflx_up_sfc(:,jb)   ,&   ! out longwave upward flux at surface [W/m2]
          & swflx_up_toa=prm_diag%swflx_up_toa(:,jb)   ,&   ! out shortwave upward flux at the TOA [W/m2]
          & swflx_up_sfc=prm_diag%swflx_up_sfc(:,jb)   ,&   ! out shortwave upward flux at the surface [W/m2]
          & swflx_par_sfc=prm_diag%swflx_par_sfc(:,jb) ,&   ! out shortwave upward flux at the surface [W/m2]
          & swflx_dn_sfc_diff=prm_diag%swflx_dn_sfc_diff(:,jb) ) ! out shortwave diffuse downward flux at the surface [W/m2]
        ENDIF

      ENDDO ! blocks

!$OMP END DO NOWAIT
!$OMP END PARALLEL

      IF (timers_level > 2) CALL timer_stop(timer_radheat)

    ENDIF  ! inwp_radiation


    !-------------------------------------------------------------------------
    !  Gravity waves drag: orographic and non-orographic
    !-------------------------------------------------------------------------

    IF (lcall_phy_jg(itsso) .OR. lcall_phy_jg(itgwd)) THEN

      IF (msg_level >= 15) &
        &  CALL message('mo_nh_interface', 'gravity waves')

      IF (timers_level > 3) CALL timer_start(timer_sso)

      ! GZ: use fast-physics time step instead of dt_phy_jg(itsso) in order to avoid calling-frequency dependence of low-level blocking
      CALL nwp_gwdrag ( dt_loc,                    & !>input
        &               lcall_phy_jg(itsso),       & !>input
        &               dt_phy_jg(itgwd),          & !>input
        &               lcall_phy_jg(itgwd),       & !>input
        &               pt_patch, p_metrics,       & !>input
        &               ext_data,                  & !>input
        &               pt_diag,                   & !>inout
        &               prm_diag, prm_nwp_tend     ) !>inout

      IF (timers_level > 3) CALL timer_stop(timer_sso)
    ENDIF ! inwp_sso
    !-------------------------------------------------------------------------


    !-------------------------------------------------------------------------
    ! Anurag Dipankar MPIM (2013-May-29)
    ! Large-scale forcing is to be applied at the end of all physics so that
    ! the most updated variable is used. Ideally it should be "next" timestep
    ! variable. Also note that its not actually a part of physics (sub-grid
    ! activity). It is called here to take advantage of u,v.
    !
    ! These LS forcing act as slow process so the tendencies from them are
    ! accumulated with the slow physics tendencies next
    !
    !(2013-25-June) LS forcing is called every physics step
    !-------------------------------------------------------------------------
    IF(is_ls_forcing)THEN

      IF (msg_level >= 15) &
        &  CALL message('mo_nh_interface:', 'LS forcing')

      IF (timers_level > 3) CALL timer_start(timer_ls_forcing)

      ! exclude boundary interpolation zone of nested domains
      rl_start = grf_bdywidth_c+1
      rl_end   = min_rlcell_int

      ! Modified by Christopher Moseley:
      ! Call to apply_ls_forcing
      CALL apply_ls_forcing ( pt_patch,          &  !>in
        &                     p_metrics,         &  !>in
        &                     p_sim_time,        &  !>in
        &                     pt_prog,           &  !>in
        &                     pt_diag,           &  !>in
        &                     pt_prog_rcf%tracer(:,:,:,iqv),  & !>in
        &                     rl_start,                       & !>in
        &                     rl_end,                         & !>in
        &                     prm_nwp_tend%ddt_u_ls,          & !>out
        &                     prm_nwp_tend%ddt_v_ls,          & !>out
        &                     prm_nwp_tend%ddt_temp_ls,       & !>out
        &                     prm_nwp_tend%ddt_tracer_ls(:,iqv),& !>out
        &                     prm_nwp_tend%ddt_temp_subs_ls,    & !>output
        &                     prm_nwp_tend%ddt_qv_subs_ls,      & !>output
        &                     prm_nwp_tend%ddt_temp_adv_ls,     & !>output
        &                     prm_nwp_tend%ddt_qv_adv_ls,       & !>output
        &                     prm_nwp_tend%ddt_temp_nud_ls,     & !>output
        &                     prm_nwp_tend%ddt_qv_nud_ls,       & !>output
        &                     prm_nwp_tend%wsub)                  !>output

      IF (timers_level > 3) CALL timer_stop(timer_ls_forcing)

    ENDIF


    IF (timers_level > 2) CALL timer_start(timer_phys_acc)
    !-------------------------------------------------------------------------
    !>  accumulate tendencies of slow_physics: Not called when LS focing is ON
    !-------------------------------------------------------------------------
    IF( (l_any_slowphys .OR. lcall_phy_jg(itradheat)) .OR. is_ls_forcing) THEN

      IF (p_test_run) THEN
        z_ddt_u_tot = 0._wp
        z_ddt_v_tot = 0._wp
      ENDIF

      IF (timers_level > 10) CALL timer_start(timer_phys_acc_1)

      ! Coefficients for extra Rayleigh friction
      ustart    = atm_phy_nwp_config(jg)%ustart_raylfric
      uoffset_q = ustart + 40._wp
      ustart_q  = ustart + 50._wp
      max_relax = -1._wp/atm_phy_nwp_config(jg)%efdt_min_raylfric

      ! exclude boundary interpolation zone of nested domains
      rl_start = grf_bdywidth_c+1
      rl_end   = min_rlcell_int

      i_startblk = pt_patch%cells%start_block(rl_start)
      i_endblk   = pt_patch%cells%end_block(rl_end)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,z_qsum,z_ddt_temp,z_ddt_alpha,vabs, &
!$OMP  rfric_fac,zddt_u_raylfric,zddt_v_raylfric,convfac,sqrt_ri,n2,dvdz2,wfac) ICON_OMP_DEFAULT_SCHEDULE
!
      DO jb = i_startblk, i_endblk
!
        CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
&                       i_startidx, i_endidx, rl_start, rl_end)


        ! artificial Rayleigh friction: active if GWD or SSO scheme is active
        IF (atm_phy_nwp_config(jg)%inwp_sso > 0 .OR. atm_phy_nwp_config(jg)%inwp_gwd > 0) THEN
          DO jk = 1, nlev
!DIR$ IVDEP
            DO jc = i_startidx, i_endidx
              vabs = SQRT(pt_diag%u(jc,jk,jb)**2 + pt_diag%v(jc,jk,jb)**2)
              rfric_fac = MAX(0._wp, 8.e-4_wp*(vabs-ustart))
              IF (vabs > ustart_q) THEN
                rfric_fac = MIN(1._wp,4.e-4_wp*(vabs-uoffset_q)**2)
              ENDIF
              zddt_u_raylfric(jc,jk) = max_relax*rfric_fac*pt_diag%u(jc,jk,jb)
              zddt_v_raylfric(jc,jk) = max_relax*rfric_fac*pt_diag%v(jc,jk,jb)
            ENDDO
          ENDDO
        ELSE
          zddt_u_raylfric(:,:) = 0._wp
          zddt_v_raylfric(:,:) = 0._wp
        ENDIF

        ! SQRT of Richardson number between the two lowest model levels
        ! This is used below to reduce frictional heating near the surface under very stable conditions
        DO jc = i_startidx, i_endidx
          n2 = 2._wp*grav/(pt_prog%theta_v(jc,nlev,jb)+pt_prog%theta_v(jc,nlev-1,jb)) * MAX(1.e-4_wp,        &
               (pt_prog%theta_v(jc,nlev-1,jb)-pt_prog%theta_v(jc,nlev,jb))/p_metrics%ddqz_z_half(jc,nlev,jb) )
          dvdz2 = MAX(1.e-6_wp, ( (pt_diag%u(jc,nlev-1,jb)-pt_diag%u(jc,nlev,jb))**2 +                      &
                  (pt_diag%v(jc,nlev-1,jb)-pt_diag%v(jc,nlev,jb))**2 )/p_metrics%ddqz_z_half(jc,nlev,jb)**2 )
          sqrt_ri(jc) = MAX(1._wp, SQRT(n2/dvdz2))
        ENDDO

        ! heating related to momentum deposition by SSO, GWD and Rayleigh friction
        DO jk = 1, nlev
!DIR$ IVDEP
          DO jc = i_startidx, i_endidx
            wfac = MIN(1._wp, 0.004_wp*p_metrics%geopot_agl(jc,jk,jb)/grav)
            prm_nwp_tend%ddt_temp_drag(jc,jk,jb) = -rcvd*(pt_diag%u(jc,jk,jb)*             &
                                                   (prm_nwp_tend%ddt_u_sso(jc,jk,jb)+      &
                                                    prm_nwp_tend%ddt_u_gwd(jc,jk,jb)+      &
                                                    zddt_u_raylfric(jc,jk))                &
                                                   +      pt_diag%v(jc,jk,jb)*             &
                                                   (prm_nwp_tend%ddt_v_sso(jc,jk,jb)+      &
                                                    prm_nwp_tend%ddt_v_gwd(jc,jk,jb)+      &
                                                    zddt_v_raylfric(jc,jk)) ) /            &
                                                    ((1._wp-wfac)*sqrt_ri(jc) + wfac)
          ENDDO
        ENDDO
#ifdef __INTEL_COMPILER
        DO jk = 1, nlev
!DIR$ IVDEP
           DO jc = i_startidx, i_endidx
        z_ddt_temp(jc,jk) =                                                      &
   &                                       prm_nwp_tend%ddt_temp_radsw(jc,jk,jb) &
   &                                    +  prm_nwp_tend%ddt_temp_radlw(jc,jk,jb) &
   &                                    +  prm_nwp_tend%ddt_temp_drag (jc,jk,jb) &
   &                                    +  prm_nwp_tend%ddt_temp_pconv(jc,jk,jb)
          ENDDO
        ENDDO
#else
!DIR$ IVDEP
        z_ddt_temp(i_startidx:i_endidx,:) =                                                      &
   &                                       prm_nwp_tend%ddt_temp_radsw(i_startidx:i_endidx,:,jb) &
   &                                    +  prm_nwp_tend%ddt_temp_radlw(i_startidx:i_endidx,:,jb) &
   &                                    +  prm_nwp_tend%ddt_temp_drag (i_startidx:i_endidx,:,jb) &
   &                                    +  prm_nwp_tend%ddt_temp_pconv(i_startidx:i_endidx,:,jb)
#endif

        IF (kstart_moist(jg) > 1) THEN
          z_qsum(:,1:kstart_moist(jg)-1)      = 0._wp
          z_ddt_alpha(:,1:kstart_moist(jg)-1) = 0._wp
        ENDIF

        DO jk = kstart_moist(jg), nlev
          DO jc = i_startidx, i_endidx

            ! summand of virtual increment
            z_qsum(jc,jk) = SUM(pt_prog_rcf%tracer (jc,jk,jb,condensate_list))

            ! tendency of virtual increment
            z_ddt_alpha(jc,jk) = vtmpc1 * prm_nwp_tend%ddt_tracer_pconv(jc,jk,jb,iqv) &
             &                 - prm_nwp_tend%ddt_tracer_pconv(jc,jk,jb,iqc)          &
             &                 - prm_nwp_tend%ddt_tracer_pconv(jc,jk,jb,iqi)
          ENDDO
        ENDDO


        ! Convert temperature tendency into Exner function tendency
        DO jk = 1, nlev
!DIR$ IVDEP
          DO jc = i_startidx, i_endidx

            pt_diag%ddt_exner_phy(jc,jk,jb) = rd_o_cpd / pt_prog%theta_v(jc,jk,jb)           &
              &                             * (z_ddt_temp(jc,jk)                             &
              &                             *(1._wp + vtmpc1*pt_prog_rcf%tracer(jc,jk,jb,iqv)&
              &                             - z_qsum(jc,jk))                                 &
              &                             + pt_diag%temp(jc,jk,jb) * z_ddt_alpha(jc,jk))
          ENDDO
        ENDDO



        ! Accumulate wind tendencies of slow physics
        ! Strictly spoken, this would not be necessary if only radiation was called
        ! in the current time step, but the radiation time step should be a multiple
        ! of the convection time step anyway in order to obtain up-to-date cloud cover fields
        IF (l_any_slowphys) THEN
#ifdef __INTEL_COMPILER
          DO jk = 1, nlev
!DIR$ IVDEP
            DO jc = i_startidx, i_endidx
              z_ddt_u_tot(jc,jk,jb) =                   &
   &              prm_nwp_tend%ddt_u_gwd     (jc,jk,jb) &
   &            + zddt_u_raylfric            (jc,jk)    &
   &            + prm_nwp_tend%ddt_u_sso     (jc,jk,jb) &
   &            + prm_nwp_tend%ddt_u_pconv   (jc,jk,jb)
              z_ddt_v_tot(jc,jk,jb) =                   &
   &              prm_nwp_tend%ddt_v_gwd     (jc,jk,jb) &
   &            + zddt_v_raylfric            (jc,jk)    &
   &            + prm_nwp_tend%ddt_v_sso     (jc,jk,jb) &
   &            + prm_nwp_tend%ddt_v_pconv   (jc,jk,jb)
          ENDDO
        ENDDO
#else
!DIR$ IVDEP
          z_ddt_u_tot(i_startidx:i_endidx,:,jb) =                   &
   &          prm_nwp_tend%ddt_u_gwd     (i_startidx:i_endidx,:,jb) &
   &        + zddt_u_raylfric            (i_startidx:i_endidx,:)    &
   &        + prm_nwp_tend%ddt_u_sso     (i_startidx:i_endidx,:,jb) &
   &        + prm_nwp_tend%ddt_u_pconv  ( i_startidx:i_endidx,:,jb)
!DIR$ IVDEP
          z_ddt_v_tot(i_startidx:i_endidx,:,jb) =                   &
   &          prm_nwp_tend%ddt_v_gwd     (i_startidx:i_endidx,:,jb) &
   &        + zddt_v_raylfric            (i_startidx:i_endidx,:)    &
   &        + prm_nwp_tend%ddt_v_sso     (i_startidx:i_endidx,:,jb) &
   &        + prm_nwp_tend%ddt_v_pconv  ( i_startidx:i_endidx,:,jb)
#endif
        ELSE IF (is_ls_forcing) THEN
          z_ddt_u_tot(i_startidx:i_endidx,:,jb) = 0._wp
          z_ddt_v_tot(i_startidx:i_endidx,:,jb) = 0._wp
        ENDIF



        !-------------------------------------------------------------------------
        !>  accumulate tendencies of slow_physics when LS forcing is ON
        !-------------------------------------------------------------------------
        IF (is_ls_forcing) THEN

          DO jk = 1, nlev
            DO jc = i_startidx, i_endidx
              z_ddt_temp(jc,jk) = z_ddt_temp(jc,jk)               &
                                +  prm_nwp_tend%ddt_temp_ls(jk)


              ! Convert temperature tendency into Exner function tendency
              z_ddt_alpha(jc,jk) = z_ddt_alpha(jc,jk)                          &
                &                + vtmpc1 * prm_nwp_tend%ddt_tracer_ls(jk,iqv) &
                &                - prm_nwp_tend%ddt_tracer_ls(jk,iqc)          &
                &                - prm_nwp_tend%ddt_tracer_ls(jk,iqi)

              pt_diag%ddt_exner_phy(jc,jk,jb) = rd_o_cpd / pt_prog%theta_v(jc,jk,jb)           &
                &                             * (z_ddt_temp(jc,jk)                             &
                &                             *(1._wp + vtmpc1*pt_prog_rcf%tracer(jc,jk,jb,iqv)&
                &                             - z_qsum(jc,jk))                                 &
                &                             + pt_diag%temp(jc,jk,jb) * z_ddt_alpha(jc,jk) )


              ! add u/v forcing tendency here
              z_ddt_u_tot(jc,jk,jb) = z_ddt_u_tot(jc,jk,jb) &
                &                   + prm_nwp_tend%ddt_u_ls(jk)

              z_ddt_v_tot(jc,jk,jb) = z_ddt_v_tot(jc,jk,jb) &
                &                   + prm_nwp_tend%ddt_v_ls(jk)

            END DO  ! jc
          END DO  ! jk

        ENDIF ! END of LS forcing tendency accumulation

        ! combine convective and EDMF rain and snow
        IF ( atm_phy_nwp_config(jg)%inwp_turb == iedmf ) THEN
#ifdef __INTEL_COMPILER
!DIR$ IVDEP
          DO jc = i_startidx, i_endidx
            prm_diag%rain_con_rate          (jc,       jb) = &
              &   prm_diag%rain_con_rate_3d (jc,nlevp1,jb)   &
              & + prm_diag%rain_edmf_rate_3d(jc,nlevp1,jb)
            prm_diag%snow_con_rate          (jc,       jb) = &
              &   prm_diag%snow_con_rate_3d (jc,nlevp1,jb)   &
              & + prm_diag%snow_edmf_rate_3d(jc,nlevp1,jb)
          ENDDO
#else
!DIR$ IVDEP
          prm_diag%rain_con_rate          (i_startidx:i_endidx,       jb) = &
            &   prm_diag%rain_con_rate_3d (i_startidx:i_endidx,nlevp1,jb)   &
            & + prm_diag%rain_edmf_rate_3d(i_startidx:i_endidx,nlevp1,jb)
!DIR$ IVDEP
          prm_diag%snow_con_rate          (i_startidx:i_endidx,       jb) = &
            &   prm_diag%snow_con_rate_3d (i_startidx:i_endidx,nlevp1,jb)   &
            & + prm_diag%snow_edmf_rate_3d(i_startidx:i_endidx,nlevp1,jb)
#endif
        ELSE IF (lcall_phy_jg(itconv)) THEN
!DIR$ IVDEP
          DO jc = i_startidx, i_endidx
            ! rain-snow conversion factor to avoid 'snow showers' at temperatures when they don't occur in practice
            convfac = MIN(1._wp,MAX(0._wp,pt_diag%temp(jc,prm_diag%k950(jc,jb),jb)-tmelt)* &
              MAX(0._wp,prm_diag%t_2m(jc,jb)-(tmelt+1.5_wp)) )
            prm_diag%rain_con_rate(jc,jb) = prm_diag%rain_con_rate_3d(jc,nlevp1,jb) + &
              convfac*prm_diag%snow_con_rate_3d(jc,nlevp1,jb)
            prm_diag%snow_con_rate(jc,jb) = (1._wp-convfac)*prm_diag%snow_con_rate_3d(jc,nlevp1,jb)
          ENDDO
        ENDIF

      ENDDO  ! jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL
      IF (timers_level > 10) CALL timer_stop(timer_phys_acc_1)

    END IF  !END OF slow physics tendency accumulation


    !--------------------------------------------------------
    ! Final section: Synchronization of updated prognostic variables,
    !                interpolation of u/v tendendies to edge points,
    !                and diagnostic computations
    !--------------------------------------------------------

    ! Synchronize tracers if any of the updating (fast-physics) processes was active.
    ! In addition, tempv needs to be synchronized, and in case of lhdiff_rcf, also exner_pr
    IF (advection_config(jg)%iadv_tke == 1) THEN
      ! TKE does not need to be synchronized if it is advected only vertically
      ntracer_sync = ntracer-1
    ELSE
      ntracer_sync = ntracer
    ENDIF

    IF (l_any_fastphys) THEN

      IF (timers_level > 10) CALL timer_start(timer_phys_sync_tracers)

      IF (use_icon_comm) THEN ! use communication library

        tracers_comm = new_icon_comm_variable(pt_prog_rcf%tracer, pt_patch%sync_cells_not_in_domain,  &
          & status=is_ready, scope=until_sync, name="pt_prog_rcf%tracer")
        tempv_comm = new_icon_comm_variable(pt_diag%tempv, pt_patch%sync_cells_not_in_domain, &
          & status=is_ready, scope=until_sync, name="pt_diag%tempv")

        IF (lhdiff_rcf) THEN
          exner_pr_comm = new_icon_comm_variable(pt_diag%exner_pr, &
            & pt_patch%sync_cells_not_in_domain, &
            & status=is_ready, scope=until_sync, name="pt_diag%exner_pr")
          IF (diffusion_config(jg)%lhdiff_w) &
            w_comm = new_icon_comm_variable(pt_prog%w, &
              & pt_patch%sync_cells_not_in_domain, &
              & status=is_ready, scope=until_sync, name="pt_prog%w")
        ENDIF

      ELSE
        IF (lhdiff_rcf .AND. diffusion_config(jg)%lhdiff_w .AND. iprog_aero == 1) THEN
          CALL sync_patch_array_mult(SYNC_C, pt_patch, ntracer_sync+4, pt_diag%tempv, pt_prog%w, &
                                     pt_diag%exner_pr, prm_diag%aerosol,                         &
                                     f4din=pt_prog_rcf%tracer(:,:,:,1:ntracer_sync))
        ELSE IF (lhdiff_rcf .AND. diffusion_config(jg)%lhdiff_w) THEN
          CALL sync_patch_array_mult(SYNC_C, pt_patch, ntracer_sync+3, pt_diag%tempv, pt_prog%w, &
                                     pt_diag%exner_pr, f4din=pt_prog_rcf%tracer(:,:,:,1:ntracer_sync))
        ELSE IF (lhdiff_rcf) THEN
          CALL sync_patch_array_mult(SYNC_C, pt_patch, ntracer_sync+2, pt_diag%tempv, &
                                     pt_diag%exner_pr, f4din=pt_prog_rcf%tracer(:,:,:,1:ntracer_sync))
        ELSE
          CALL sync_patch_array_mult(SYNC_C, pt_patch, ntracer_sync+1, pt_diag%tempv, &
                                     f4din=pt_prog_rcf%tracer(:,:,:,1:ntracer_sync))
        ENDIF

      ENDIF

      IF (timers_level > 10) THEN
        CALL timer_stop(timer_phys_sync_tracers)
      ENDIF
    ELSE IF (linit .AND. advection_config(jg)%iadv_tke == 2) THEN
      CALL sync_patch_array(SYNC_C, pt_patch, pt_prog_rcf%tracer(:,:,:,iqtke))
    ENDIF

    !------------------------------------------------------------
    ! sync here the slowphys for aggregation
    !-------------------------------------------------------------------
    IF (use_physics_barrier) THEN
      CALL timer_start(timer_barrier)
      CALL work_mpi_barrier()
      CALL timer_stop(timer_barrier)
    ENDIF
    !-------------------------------------------------------------------
    IF (timers_level > 10) CALL timer_start(timer_phys_sync_ddt_u)
    IF (use_icon_comm) THEN

      IF (lcall_phy_jg(itturb) ) THEN
        ddt_u_tot_comm = new_icon_comm_variable(prm_nwp_tend%ddt_u_turb, &
          & pt_patch%sync_cells_one_edge_in_domain, status=is_ready, scope=until_sync, &
          & name="prm_nwp_tend%ddt_u_turb")
        ddt_v_tot_comm = new_icon_comm_variable(prm_nwp_tend%ddt_v_turb, &
          & pt_patch%sync_cells_one_edge_in_domain, status=is_ready, scope=until_sync, &
          & name="prm_nwp_tend%ddt_v_turb")

        IF ( l_any_slowphys .OR. is_ls_forcing ) THEN
          z_ddt_u_tot_comm = new_icon_comm_variable(z_ddt_u_tot, &
            & pt_patch%sync_cells_one_edge_in_domain, &
            & status=is_ready, scope=until_sync, name="z_ddt_u_tot")
          z_ddt_v_tot_comm = new_icon_comm_variable(z_ddt_v_tot, &
            & pt_patch%sync_cells_one_edge_in_domain, &
            & status=is_ready, scope=until_sync, name="z_ddt_v_tot")
        ENDIF
      ENDIF

       ! sync everything here
      CALL icon_comm_sync_all()

    ELSE

      IF ( (is_ls_forcing .OR. l_any_slowphys) .AND. lcall_phy_jg(itturb) ) THEN

        CALL sync_patch_array_mult(SYNC_C1, pt_patch, 4, z_ddt_u_tot, z_ddt_v_tot, &
                                 prm_nwp_tend%ddt_u_turb, prm_nwp_tend%ddt_v_turb)

      ELSE IF (lcall_phy_jg(itturb) ) THEN

        CALL sync_patch_array_mult(SYNC_C1, pt_patch, 2, prm_nwp_tend%ddt_u_turb, &
                                 prm_nwp_tend%ddt_v_turb)
      ENDIF
    ENDIF

    IF (timers_level > 10) CALL timer_stop(timer_phys_sync_ddt_u)
    !------------------------------------------------------------


    !------------------------------------------------------------
    ! compute on the halos
    IF (timers_level > 10) CALL timer_start(timer_phys_acc_par)
    IF (l_any_fastphys) THEN
      IF (my_process_is_mpi_all_parallel() ) THEN

        rl_start = min_rlcell_int-1
        rl_end   = min_rlcell

        i_startblk = pt_patch%cells%start_block(rl_start)
        i_endblk   = pt_patch%cells%end_block(rl_end)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx, i_endidx) ICON_OMP_DEFAULT_SCHEDULE

        DO jb = i_startblk, i_endblk
          CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
            & i_startidx, i_endidx, rl_start, rl_end )

          IF (lhdiff_rcf) THEN
            DO jk = 1, nlev
!DIR$ IVDEP
              DO jc =  i_startidx, i_endidx

                IF (p_metrics%mask_prog_halo_c(jc,jb)) THEN
                  pt_prog%exner(jc,jk,jb) = EXP(rd_o_cpd*LOG(rd_o_p0ref                   &
                    &                     * pt_prog%rho(jc,jk,jb)*pt_diag%tempv(jc,jk,jb)))

                  pt_prog%theta_v(jc,jk,  jb) = pt_diag%tempv(jc,jk,jb) &
&                                             / pt_prog%exner(jc,jk,jb)

                ENDIF

              ENDDO
            ENDDO
          ELSE
            DO jk = 1, nlev
!DIR$ IVDEP
              DO jc =  i_startidx, i_endidx

                IF (p_metrics%mask_prog_halo_c(jc,jb)) THEN
                  pt_prog%exner(jc,jk,jb) = EXP(rd_o_cpd*LOG(rd_o_p0ref                   &
                    &                     * pt_prog%rho(jc,jk,jb)*pt_diag%tempv(jc,jk,jb)))

                  pt_diag%exner_pr(jc,jk,jb) = pt_diag%exner_pr(jc,jk,jb) + &
                    pt_prog%exner(jc,jk,jb) - z_exner_sv(jc,jk,jb)

                  pt_prog%theta_v(jc,jk,  jb) = pt_diag%tempv(jc,jk,jb) &
&                                             / pt_prog%exner(jc,jk,jb)

                ENDIF

              ENDDO
            ENDDO
          ENDIF
        ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

      ENDIF ! my_process_is_mpi_all_parallel
    ENDIF ! fast-physics synchronization
    IF (timers_level > 10) CALL timer_stop(timer_phys_acc_par)


    !-------------------------------------------------------------------------
    !>
    !!    @par Interpolation from  u,v onto v_n
    !!      ddt_vn_phy  =interpol(ddt_u_tot)+interpol(ddt_v_tot)
    !!      Calculate normal velocity at edge midpoints
    !-------------------------------------------------------------------------

    IF (timers_level > 10)  CALL timer_start(timer_phys_acc_2)
!$OMP PARALLEL PRIVATE(rl_start,rl_end,i_startblk,i_endblk)

    ! exclude boundary interpolation zone of nested domains
    rl_start = grf_bdywidth_e+1
    rl_end   = min_rledge_int

    i_startblk = pt_patch%edges%start_block(rl_start)
    i_endblk   = pt_patch%edges%end_block(rl_end)


!$OMP DO PRIVATE(jb,jk,jce,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE

    DO jb = i_startblk, i_endblk

      CALL get_indices_e(pt_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)

      IF ( (is_ls_forcing .OR. l_any_slowphys) .AND. lcall_phy_jg(itturb) ) THEN

#ifdef __LOOP_EXCHANGE
        DO jce = i_startidx, i_endidx
!DIR$ IVDEP
          DO jk = 1, nlev
#else
!CDIR UNROLL=5
        DO jk = 1, nlev
          DO jce = i_startidx, i_endidx
#endif

            pt_diag%ddt_vn_phy(jce,jk,jb) =   pt_int_state%c_lin_e(jce,1,jb)           &
&                                 * ( z_ddt_u_tot(iidx(jce,jb,1),jk,iblk(jce,jb,1))    &
&                                   *  pt_patch%edges%primal_normal_cell(jce,jb,1)%v1  &
&                                   + z_ddt_v_tot(iidx(jce,jb,1),jk,iblk(jce,jb,1))    &
&                                   *  pt_patch%edges%primal_normal_cell(jce,jb,1)%v2 )&
&                                                 + pt_int_state%c_lin_e(jce,2,jb)     &
&                                 * ( z_ddt_u_tot(iidx(jce,jb,2),jk,iblk(jce,jb,2))    &
&                                    * pt_patch%edges%primal_normal_cell(jce,jb,2)%v1  &
&                                  +  z_ddt_v_tot(iidx(jce,jb,2),jk,iblk(jce,jb,2))    &
&                                   *  pt_patch%edges%primal_normal_cell(jce,jb,2)%v2 )

            pt_prog%vn(jce,jk,jb) = pt_prog%vn(jce,jk,jb) + dt_loc * (                 &
                                              pt_int_state%c_lin_e(jce,1,jb)           &
&                     * ( prm_nwp_tend%ddt_u_turb(iidx(jce,jb,1),jk,iblk(jce,jb,1))    &
&                                   *  pt_patch%edges%primal_normal_cell(jce,jb,1)%v1  &
&                       + prm_nwp_tend%ddt_v_turb(iidx(jce,jb,1),jk,iblk(jce,jb,1))    &
&                                   *  pt_patch%edges%primal_normal_cell(jce,jb,1)%v2 )&
&                                                 + pt_int_state%c_lin_e(jce,2,jb)     &
&                     * ( prm_nwp_tend%ddt_u_turb(iidx(jce,jb,2),jk,iblk(jce,jb,2))    &
&                                    * pt_patch%edges%primal_normal_cell(jce,jb,2)%v1  &
&                      +  prm_nwp_tend%ddt_v_turb(iidx(jce,jb,2),jk,iblk(jce,jb,2))    &
&                                  *  pt_patch%edges%primal_normal_cell(jce,jb,2)%v2 ) )

          ENDDO
        ENDDO

      ELSE IF (lcall_phy_jg(itturb) ) THEN
#ifdef __LOOP_EXCHANGE
        DO jce = i_startidx, i_endidx
!DIR$ IVDEP
          DO jk = 1, nlev
#else
!CDIR UNROLL=8
        DO jk = 1, nlev
          DO jce = i_startidx, i_endidx
#endif

            pt_prog%vn(jce,jk,jb) = pt_prog%vn(jce,jk,jb) + dt_loc * (                 &
                                              pt_int_state%c_lin_e(jce,1,jb)           &
&                     * ( prm_nwp_tend%ddt_u_turb(iidx(jce,jb,1),jk,iblk(jce,jb,1))    &
&                                   *  pt_patch%edges%primal_normal_cell(jce,jb,1)%v1  &
&                       + prm_nwp_tend%ddt_v_turb(iidx(jce,jb,1),jk,iblk(jce,jb,1))    &
&                                   *  pt_patch%edges%primal_normal_cell(jce,jb,1)%v2 )&
&                                                 + pt_int_state%c_lin_e(jce,2,jb)     &
&                     * ( prm_nwp_tend%ddt_u_turb(iidx(jce,jb,2),jk,iblk(jce,jb,2))    &
&                                    * pt_patch%edges%primal_normal_cell(jce,jb,2)%v1  &
&                      +  prm_nwp_tend%ddt_v_turb(iidx(jce,jb,2),jk,iblk(jce,jb,2))    &
&                                  *  pt_patch%edges%primal_normal_cell(jce,jb,2)%v2 ) )

          ENDDO
        ENDDO

      ENDIF

    ENDDO
!$OMP END DO
!$OMP END PARALLEL

    IF (timers_level > 10) CALL timer_stop(timer_phys_acc_2)


    IF (timers_level > 10) CALL timer_start(timer_phys_dpsdt)
    !
    ! dpsdt diagnostic
    IF (lcalc_dpsdt) THEN
      CALL compute_dpsdt (pt_patch = pt_patch, &
        &                 dt       = dt_loc,   &
        &                 pt_diag  = pt_diag   )
    ENDIF
    IF (timers_level > 10) CALL timer_stop(timer_phys_dpsdt)


    IF (timers_level > 10) CALL timer_start(timer_phys_sync_vn)
    IF (lcall_phy_jg(itturb)) CALL sync_patch_array(SYNC_E, pt_patch, pt_prog%vn)
    IF (timers_level > 10) CALL timer_stop(timer_phys_sync_vn)
    IF (timers_level > 2) CALL timer_stop(timer_phys_acc)



    IF (lcall_phy_jg(itturb) .AND. msg_level >= 18) THEN ! extended diagnostic for turbulence quantities
      CALL nwp_diag_output_2(pt_patch, pt_prog_rcf, prm_nwp_tend)
    ENDIF


    ! time averages, accumulations and vertical integrals
    CALL nwp_statistics(lcall_phy_jg,                    & !in
                        & dt_phy_jg,p_sim_time,          & !in
                        & kstart_moist(jg),              & !in
                        & ih_clch(jg), ih_clcm(jg),      & !in
                        & pt_patch, p_metrics,           & !in
                        & pt_prog, pt_prog_rcf,          & !in
                        & pt_diag,                       & !inout
                        & prm_diag                       ) !inout

    ! Call the ART diagnostics
    CALL art_diagnostics_interface(pt_patch,               &
      &                            pt_prog%rho,            &
      &                            pt_diag%pres,           &
      &                            pt_prog_now_rcf%tracer, &
      &                            p_metrics%ddqz_z_full,  &
      &                            p_metrics%z_mc, jg,     &
      &                            dt_phy_jg, p_sim_time)

    IF (ltimer) CALL timer_stop(timer_physics)


  END SUBROUTINE nwp_nh_interface

  !-------------------------------------------------------------------------
  !-------------------------------------------------------------------------


END MODULE mo_nh_interface_nwp

