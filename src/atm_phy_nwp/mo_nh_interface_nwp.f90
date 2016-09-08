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

MODULE mo_nh_interface_nwp

  USE mo_datetime,                ONLY: t_datetime
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
  USE mo_nonhydrostatic_config,   ONLY: kstart_moist, lhdiff_rcf, ih_clch, ih_clcm
  USE mo_nwp_lnd_types,           ONLY: t_lnd_prog, t_wtr_prog, t_lnd_diag
  USE mo_ext_data_types,          ONLY: t_external_data
  USE mo_nwp_phy_types,           ONLY: t_nwp_phy_diag, t_nwp_phy_tend
  USE mo_parallel_config,         ONLY: nproma, p_test_run, use_icon_comm, use_physics_barrier
  USE mo_diffusion_config,        ONLY: diffusion_config
  USE mo_run_config,              ONLY: ntracer, iqv, iqc, iqi, iqs, iqtvar, iqtke,  &
    &                                   msg_level, ltimer, timers_level, lart
  USE mo_grid_config,             ONLY: l_limited_area
  USE mo_physical_constants,      ONLY: rd, rd_o_cpd, vtmpc1, p0ref, rcvd, cvd, cvv

  USE mo_nh_diagnose_pres_temp,   ONLY: diagnose_pres_temp, diag_pres, diag_temp

  USE mo_atm_phy_nwp_config,      ONLY: atm_phy_nwp_config, iprog_aero
  USE mo_util_phys,               ONLY: nh_update_tracer_phy
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
  USE mo_nwp_turb_sfc_interface,  ONLY: nwp_turbulence_sfc
  USE mo_nwp_sfc_interface,       ONLY: nwp_surface
  USE mo_nwp_conv_interface,      ONLY: nwp_convection
  USE mo_nwp_rad_interface,       ONLY: nwp_radiation
  USE mo_sync,                    ONLY: sync_patch_array, sync_patch_array_mult, SYNC_E,      &
                                        SYNC_C, SYNC_C1, global_sum_array
  USE mo_mpi,                     ONLY: my_process_is_mpi_all_parallel, work_mpi_barrier,     &
    &                                   process_mpi_stdio_id, my_process_is_stdio
  USE mo_nwp_diagnosis,           ONLY: nwp_statistics, nwp_diag_output_1, nwp_diag_output_2
  USE mo_icon_comm_lib,           ONLY: new_icon_comm_variable,                               &
    &                                   icon_comm_sync_all, is_ready, until_sync
  USE mo_art_washout_interface,   ONLY: art_washout_interface
  USE mo_art_reaction_interface,  ONLY: art_reaction_interface
  USE mo_linked_list,             ONLY: t_var_list
  USE mo_ls_forcing_nml,          ONLY: is_ls_forcing
  USE mo_ls_forcing,              ONLY: apply_ls_forcing
  USE mo_advection_config,        ONLY: advection_config
  USE mo_o3_util,                 ONLY: calc_o3_gems

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
                            & p_sim_time, datetime,                & !input
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
    REAL(wp),INTENT(in)          :: p_sim_time

    TYPE(t_datetime),            INTENT(in):: datetime
    TYPE(t_patch),        TARGET,INTENT(in):: pt_patch     !<grid/patch info.
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
    INTEGER :: i_nchdom                !< domain index

    ! Local scalars:

    INTEGER :: jc,jk,jb,jce      !loop indices
    INTEGER :: jg                !domain id

    LOGICAL :: ltemp, lpres, ltemp_ifc, l_any_fastphys, l_any_slowphys

    INTEGER,  POINTER ::  iidx(:,:,:), iblk(:,:,:)

    REAL(wp), TARGET :: &                                              !> temporal arrays for
      & z_ddt_u_tot (nproma,pt_patch%nlev,pt_patch%nblks_c),&
      & z_ddt_v_tot (nproma,pt_patch%nlev,pt_patch%nblks_c),& !< hor. wind tendencies
      & z_ddt_temp  (nproma,pt_patch%nlev)   !< Temperature tendency

    REAL(wp) :: z_exner_sv(nproma,pt_patch%nlev,pt_patch%nblks_c), z_tempv, &
      zddt_u_raylfric(nproma,pt_patch%nlev), zddt_v_raylfric(nproma,pt_patch%nlev)

    !< vertical interfaces

    REAL(wp) :: zsct ! solar constant (at time of year)
    REAL(wp) :: zcosmu0 (nproma,pt_patch%nblks_c)

    REAL(wp) :: z_qsum(nproma,pt_patch%nlev)  !< summand of virtual increment
    REAL(wp) :: z_ddt_qsum                    !< summand of tendency of virtual increment

    ! auxiliaries for Rayleigh friction computation
    REAL(wp) :: vabs, rfric_fac, ustart, uoffset_q, ustart_q, max_relax

    ! Variables for dpsdt diagnostic
    REAL(wp) :: dps_blk(pt_patch%nblks_c), dpsdt_avg
    INTEGER  :: npoints_blk(pt_patch%nblks_c), npoints

    ! Variables for EDMF DUALM
    REAL(wp) :: qtvar(nproma,pt_patch%nlev)

    ! communication ids, these do not need to be different variables,
    ! since they are not treated individualy
    INTEGER :: ddt_u_tot_comm, ddt_v_tot_comm, z_ddt_u_tot_comm, z_ddt_v_tot_comm, &
      & tracers_comm, tempv_comm, exner_old_comm, w_comm

    INTEGER :: ntracer_sync

    ! Pointer to IDs of tracers which contain prognostic condensate.
    ! Required for computing the water loading term 
    INTEGER, POINTER :: condensate_list(:)



    IF (ltimer) CALL timer_start(timer_physics)

    ! local variables related to the blocking

    i_nchdom  = MAX(1,pt_patch%n_childdom)
    jg        = pt_patch%id

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
    condensate_list => advection_config(jg)%ilist_hydroMass



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

        ! The provisional "new" tracer state, resulting from the advection 
        ! step, still needs to be updated with the SLOW-physics tracer tendencies 
        ! computed at the end of the last physics call for the then final 
        ! "new" state. The corresponding update for the dynamics variables has 
        ! already happened in the dynamical core.
        !
        ! In addition: Update qv with DA increment (during IAU phase)
        CALL nh_update_tracer_phy(pt_patch          ,& !in
           &                  dt_phy_jg(itfastphy)  ,& !in
           &                  pt_diag               ,& !in
           &                  p_metrics             ,& !in
           &                  prm_nwp_tend          ,& !in
           &                  prm_diag              ,& !inout phyfields
           &                  pt_prog_rcf           ,& !inout tracer
           &                  pt_prog               ,& !in density
           &                  jb, i_startidx, i_endidx ) !in
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



    !For turbulence schemes NOT including the call to the surface scheme.
    !nwp_surface must even be called in inwp_surface = 0 because the
    !the lower boundary conditions for the turbulence scheme
    !are not set otherwise

    IF ( l_any_fastphys .AND. ANY( (/icosmo,igme/)==atm_phy_nwp_config(jg)%inwp_turb ) ) THEN
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
      CASE(icosmo,igme)

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

      !Turbulence schemes including the call to the surface scheme
      CASE(iedmf)

        CALL nwp_turbulence_sfc (  dt_phy_jg(itfastphy),              & !>input
                                  & pt_patch, p_metrics,              & !>input
                                  & ext_data,                         & !>input
                                  & pt_prog,                          & !>inout
                                  & pt_prog_now_rcf, pt_prog_rcf,     & !>in/inout
                                  & pt_diag ,                         & !>inout
                                  & prm_diag, prm_nwp_tend,           & !>inout
                                  & lnd_prog_now, lnd_prog_new,       & !>inout
                                  & wtr_prog_now, wtr_prog_new,       & !>inout
                                  & lnd_diag                          ) !>inout

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

      IF (timers_level > 1) CALL timer_start(timer_nwp_microphysics)

      CALL nwp_microphysics ( dt_phy_jg(itfastphy),             & !>input
                            & pt_patch, p_metrics,              & !>input
                            & pt_prog,                          & !>inout
                            & pt_prog_rcf,                      & !>inout
                            & pt_diag ,                         & !>inout
                            & prm_diag                          ) !>inout

      IF (timers_level > 1) CALL timer_stop(timer_nwp_microphysics)

    ENDIF

    IF (lart) THEN
      CALL calc_o3_gems(pt_patch,datetime,pt_diag,prm_diag,ext_data)

      IF (.NOT. linit) THEN
        CALL art_reaction_interface(ext_data,              & !> in
                &                   pt_patch,              & !> in
                &                   datetime,              & !> in
                &                   dt_phy_jg(itfastphy),  & !> in
                &                   p_prog_list,           & !> in
                &                   pt_prog,               & !> in
                &                   p_metrics,             & !> in
                &                   prm_diag,              & !> in
                &                   pt_diag,               & !> inout
                &                   pt_prog_rcf%tracer)
      END IF

      CALL art_washout_interface(pt_prog,pt_diag,              & !>in
                &          dt_phy_jg(itfastphy),               & !>in
                &          pt_patch,                           & !>in
                &          prm_diag,                           & !>in
                &          pt_prog_rcf%tracer)                   !>inout
    ENDIF !lart


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

    i_startblk = pt_patch%cells%start_blk(rl_start,1)
    i_endblk   = pt_patch%cells%end_blk(rl_end,i_nchdom)

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

            pt_diag%exner_old(jc,jk,jb) = pt_diag%exner_old(jc,jk,jb) + &
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

        ! reset dynamical exner increment to zero
        ! (it is accumulated over one advective time step in solve_nh)
        pt_diag%exner_dyn_incr(:,kstart_moist(jg):nlev,jb) = 0._wp

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

      i_startblk = pt_patch%cells%start_blk(rl_start,1)
      i_endblk   = pt_patch%cells%end_blk(rl_end,i_nchdom)

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
&              pmfude_rate = prm_diag%con_udd(:,:,jb,3)   ,       & !! in:  convective updraft detrainment rate
&              plu         = prm_diag%con_udd(:,:,jb,7)   ,       & !! in:  updraft condensate
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
           &              datetime,              & ! in
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
        & zsct       = zsct )
      IF (timers_level > 10) CALL timer_stop(timer_pre_radiation_nwp)

      ! exclude boundary interpolation zone of nested domains
      rl_start = grf_bdywidth_c+1
      rl_end   = min_rlcell_int

      i_startblk = pt_patch%cells%start_blk(rl_start,1)
      i_endblk   = pt_patch%cells%end_blk(rl_end,i_nchdom)

      IF (timers_level > 2) CALL timer_start(timer_radheat)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
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
          & pqv=prm_diag%tot_cld(:,:,jb,iqv)       ,&! in     specific moisture           [kg/kg]
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
          & cosmu0=zcosmu0(:,jb)                   ,&! in     cosine of solar zenith angle
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
          & pqv=prm_diag%tot_cld(:,:,jb,iqv)       ,&! in     specific moisture           [kg/kg]
          & pcd=cvd                                ,&! in     specific heat of dry air  [J/kg/K]
          & pcv=cvv                                ,&! in     specific heat of vapor    [J/kg/K]
          & pi0=prm_diag%flxdwswtoa(:,jb)          ,&! in     solar incoming flux at TOA  [W/m2]
          & pemiss=ext_data%atm%emis_rad(:,jb)     ,&! in     lw sfc emissivity
          & pqc=prm_diag%tot_cld    (:,:,jb,iqc)   ,&! in     specific cloud water        [kg/kg]
          & pqi=prm_diag%tot_cld    (:,:,jb,iqi)   ,&! in     specific cloud ice          [kg/kg]
          & ppres_ifc=pt_diag%pres_ifc(:,:,jb)     ,&! in     pressure at layer boundaries [Pa]
          & cosmu0=zcosmu0(:,jb)                   ,&! in     cosine of solar zenith angle
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

      CALL nwp_gwdrag ( dt_phy_jg(itsso),          & !>input
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

      CALL apply_ls_forcing ( pt_patch,          &  !>in
        &                     p_metrics,         &  !>in
        &                     pt_prog,           &  !>in
        &                     pt_diag,           &  !>in
        &                     pt_prog_rcf%tracer(:,:,:,iqv),  & !>in
        &                     rl_start,                       & !>in
        &                     rl_end,                         & !>in
        &                     prm_nwp_tend%ddt_u_ls,          & !>out
        &                     prm_nwp_tend%ddt_v_ls,          & !>out
        &                     prm_nwp_tend%ddt_temp_ls,       & !>out
        &                     prm_nwp_tend%ddt_tracer_ls(:,iqv) ) !>out

      IF (timers_level > 3) CALL timer_stop(timer_ls_forcing)

    ENDIF


    IF (timers_level > 2) CALL timer_start(timer_phys_acc)
    !-------------------------------------------------------------------------
    !>  accumulate tendencies of slow_physics
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

      i_startblk = pt_patch%cells%start_blk(rl_start,1)
      i_endblk   = pt_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,z_qsum,z_ddt_temp,z_ddt_qsum,vabs, &
!$OMP  rfric_fac,zddt_u_raylfric,zddt_v_raylfric) ICON_OMP_DEFAULT_SCHEDULE
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

        ! heating related to momentum deposition by SSO, GWD and Rayleigh friction
        DO jk = 1, nlev
!DIR$ IVDEP
          DO jc = i_startidx, i_endidx
            prm_nwp_tend%ddt_temp_drag(jc,jk,jb) = -rcvd*(pt_diag%u(jc,jk,jb)*             &
                                                   (prm_nwp_tend%ddt_u_sso(jc,jk,jb)+      &
                                                    prm_nwp_tend%ddt_u_gwd(jc,jk,jb)+      &
                                                    zddt_u_raylfric(jc,jk))                &
                                                   +      pt_diag%v(jc,jk,jb)*             &
                                                   (prm_nwp_tend%ddt_v_sso(jc,jk,jb)+      &
                                                    prm_nwp_tend%ddt_v_gwd(jc,jk,jb)+      &
                                                    zddt_v_raylfric(jc,jk))                )
          ENDDO
        ENDDO
!DIR$ IVDEP
        z_ddt_temp(i_startidx:i_endidx,:) =                                                      &
   &                                       prm_nwp_tend%ddt_temp_radsw(i_startidx:i_endidx,:,jb) &
   &                                    +  prm_nwp_tend%ddt_temp_radlw(i_startidx:i_endidx,:,jb) &
   &                                    +  prm_nwp_tend%ddt_temp_drag (i_startidx:i_endidx,:,jb) &
   &                                    +  prm_nwp_tend%ddt_temp_pconv(i_startidx:i_endidx,:,jb)



        IF (kstart_moist(jg) > 1) z_qsum(:,1:kstart_moist(jg)-1) = 0._wp

        DO jk = kstart_moist(jg), nlev
          DO jc = i_startidx, i_endidx

            z_qsum(jc,jk) = SUM(pt_prog_rcf%tracer (jc,jk,jb,condensate_list))

          ENDDO
        ENDDO


        ! Convert temperature tendency into Exner function tendency
        DO jk = 1, nlev
!DIR$ IVDEP
          DO jc = i_startidx, i_endidx

            z_ddt_qsum =   prm_nwp_tend%ddt_tracer_pconv(jc,jk,jb,iqc) &
              &          + prm_nwp_tend%ddt_tracer_pconv(jc,jk,jb,iqi)

            pt_diag%ddt_exner_phy(jc,jk,jb) = rd_o_cpd / pt_prog%theta_v(jc,jk,jb)           &
              &                             * (z_ddt_temp(jc,jk)                             &
              &                             *(1._wp + vtmpc1*pt_prog_rcf%tracer(jc,jk,jb,iqv)&
              &                             - z_qsum(jc,jk)) + pt_diag%temp(jc,jk,jb)        &
              &           * (vtmpc1 * prm_nwp_tend%ddt_tracer_pconv(jc,jk,jb,iqv)-z_ddt_qsum ))

          ENDDO
        ENDDO

        ! Accumulate wind tendencies of slow physics
        ! Strictly spoken, this would not be necessary if only radiation was called
        ! in the current time step, but the radiation time step should be a multiple
        ! of the convection time step anyway in order to obtain up-to-date cloud cover fields
        IF (l_any_slowphys) THEN
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
              z_ddt_qsum =   prm_nwp_tend%ddt_tracer_pconv(jc,jk,jb,iqc) &
                &          + prm_nwp_tend%ddt_tracer_pconv(jc,jk,jb,iqi) &
                &          + prm_nwp_tend%ddt_tracer_ls(jk,iqc)          &
                &          + prm_nwp_tend%ddt_tracer_ls(jk,iqi)

              pt_diag%ddt_exner_phy(jc,jk,jb) = rd_o_cpd / pt_prog%theta_v(jc,jk,jb)           &
                &                             * (z_ddt_temp(jc,jk)                             &
                &                             *(1._wp + vtmpc1*pt_prog_rcf%tracer(jc,jk,jb,iqv)&
                &                             - z_qsum(jc,jk)) + pt_diag%temp(jc,jk,jb)        &
                &           * (vtmpc1 * (prm_nwp_tend%ddt_tracer_pconv(jc,jk,jb,iqv) +         &
                &              prm_nwp_tend%ddt_tracer_ls(jk,iqv) ) - z_ddt_qsum ) )


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
!DIR$ IVDEP
          prm_diag%rain_con_rate          (i_startidx:i_endidx,       jb) = &
            &   prm_diag%rain_con_rate_3d (i_startidx:i_endidx,nlevp1,jb)   &
            & + prm_diag%rain_edmf_rate_3d(i_startidx:i_endidx,nlevp1,jb)
!DIR$ IVDEP
          prm_diag%snow_con_rate          (i_startidx:i_endidx,       jb) = &
            &   prm_diag%snow_con_rate_3d (i_startidx:i_endidx,nlevp1,jb)   &
            & + prm_diag%snow_edmf_rate_3d(i_startidx:i_endidx,nlevp1,jb)
        ELSE IF (lcall_phy_jg(itconv)) THEN
!DIR$ IVDEP
          prm_diag%rain_con_rate          (i_startidx:i_endidx,       jb) = &
            &   prm_diag%rain_con_rate_3d (i_startidx:i_endidx,nlevp1,jb)
!DIR$ IVDEP
          prm_diag%snow_con_rate          (i_startidx:i_endidx,       jb) = &
            &   prm_diag%snow_con_rate_3d (i_startidx:i_endidx,nlevp1,jb)
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
    ! In addition, tempv needs to be synchronized, and in case of lhdiff_rcf, also exner_old
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
          exner_old_comm = new_icon_comm_variable(pt_diag%exner_old, &
            & pt_patch%sync_cells_not_in_domain, &
            & status=is_ready, scope=until_sync, name="pt_diag%exner_old")
          IF (diffusion_config(jg)%lhdiff_w) &
            w_comm = new_icon_comm_variable(pt_prog%w, &
              & pt_patch%sync_cells_not_in_domain, &
              & status=is_ready, scope=until_sync, name="pt_prog%w")
        ENDIF

      ELSE
        IF (lhdiff_rcf .AND. diffusion_config(jg)%lhdiff_w .AND. iprog_aero == 1) THEN
          CALL sync_patch_array_mult(SYNC_C, pt_patch, ntracer_sync+4, pt_diag%tempv, pt_prog%w, &
                                     pt_diag%exner_old, prm_diag%aerosol,                        &
                                     f4din=pt_prog_rcf%tracer(:,:,:,1:ntracer_sync))
        ELSE IF (lhdiff_rcf .AND. diffusion_config(jg)%lhdiff_w) THEN
          CALL sync_patch_array_mult(SYNC_C, pt_patch, ntracer_sync+3, pt_diag%tempv, pt_prog%w, &
                                     pt_diag%exner_old, f4din=pt_prog_rcf%tracer(:,:,:,1:ntracer_sync))
        ELSE IF (lhdiff_rcf) THEN
          CALL sync_patch_array_mult(SYNC_C, pt_patch, ntracer_sync+2, pt_diag%tempv, &
                                     pt_diag%exner_old, f4din=pt_prog_rcf%tracer(:,:,:,1:ntracer_sync))
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

        i_startblk = pt_patch%cells%start_blk(rl_start,1)
        i_endblk   = pt_patch%cells%end_blk(rl_end,i_nchdom)

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

                  pt_diag%exner_old(jc,jk,jb) = pt_diag%exner_old(jc,jk,jb) + &
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


    ! Initialize fields for runtime diagnostics
    ! In case that average ABS(dpsdt) is diagnosed
    IF (msg_level >= 11) THEN
      dps_blk(:)   = 0._wp
      npoints_blk(:) = 0
    ENDIF

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

    i_startblk = pt_patch%edges%start_blk(rl_start,1)
    i_endblk   = pt_patch%edges%end_blk(rl_end,i_nchdom)


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


    ! Diagnosis of ABS(dpsdt) if msg_level >= 11
    IF (msg_level >= 11) THEN

      rl_start = grf_bdywidth_c+1
      rl_end   = min_rlcell_int

      i_startblk = pt_patch%cells%start_blk(rl_start,1)
      i_endblk   = pt_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, rl_start, rl_end)

        DO jc = i_startidx, i_endidx
          ! Note: division by time step follows below
          dps_blk(jb) = dps_blk(jb) + &
            ABS(pt_diag%pres_sfc(jc,jb)-pt_diag%pres_sfc_old(jc,jb))
          npoints_blk(jb) = npoints_blk(jb) + 1
          pt_diag%pres_sfc_old(jc,jb) = pt_diag%pres_sfc(jc,jb)
        ENDDO
      ENDDO
!$OMP END DO NOWAIT
    ENDIF

!$OMP END PARALLEL

    IF (timers_level > 10) CALL timer_stop(timer_phys_acc_2)
    IF (timers_level > 10) CALL timer_start(timer_phys_sync_vn)
    IF (lcall_phy_jg(itturb)) CALL sync_patch_array(SYNC_E, pt_patch, pt_prog%vn)
    IF (timers_level > 10) CALL timer_stop(timer_phys_sync_vn)
    IF (timers_level > 2) CALL timer_stop(timer_phys_acc)


    ! dpsdt diagnostic - omitted in the case of a parallization test (p_test_run) because this
    ! is a purely diagnostic quantitiy, for which it does not make sense to implement an order-invariant
    ! summation
    IF (.NOT. p_test_run .AND. msg_level >= 11) THEN
      dpsdt_avg = SUM(dps_blk)
      npoints   = SUM(npoints_blk)
      dpsdt_avg = global_sum_array(dpsdt_avg, opt_iroot=process_mpi_stdio_id)
      npoints   = global_sum_array(npoints  , opt_iroot=process_mpi_stdio_id)
      IF (my_process_is_stdio()) THEN
        dpsdt_avg = dpsdt_avg/(REAL(npoints,wp)*dt_loc)
        ! Exclude initial time step where pres_sfc_old is zero
        IF (dpsdt_avg < 10000._wp/dt_loc) THEN
          WRITE(message_text,'(a,f12.6,a,i3)') 'average |dPS/dt| =',dpsdt_avg,' Pa/s in domain',jg
          CALL message('nwp_nh_interface: ', TRIM(message_text))
        ENDIF
      ENDIF
    ENDIF

    IF (msg_level >= 18) THEN ! extended diagnostic
      CALL nwp_diag_output_2(pt_patch, pt_prog_rcf, prm_nwp_tend, lcall_phy_jg(itturb))
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

    IF (ltimer) CALL timer_stop(timer_physics)


  END SUBROUTINE nwp_nh_interface

  !-------------------------------------------------------------------------
  !-------------------------------------------------------------------------


END MODULE mo_nh_interface_nwp

