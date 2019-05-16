!>
!! @brief branch for the non-hydrostatic ICON workflow
!!
!! @author Kristina Froehlich, MPI-M (2011-07-19)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_atmo_nonhydrostatic

USE mo_kind,                 ONLY: wp
USE mo_exception,            ONLY: message, finish, print_value
USE mtime,                   ONLY: datetimeToString, OPERATOR(>)
USE mo_fortran_tools,        ONLY: copy, init
USE mo_impl_constants,       ONLY: SUCCESS, max_dom, inwp, iecham
USE mo_timer,                ONLY: timers_level, timer_start, timer_stop, timer_init_latbc, &
  &                                timer_model_init, timer_init_icon, timer_read_restart
USE mo_master_config,        ONLY: isRestart
USE mo_time_config,          ONLY: time_config
USE mo_load_restart,         ONLY: read_restart_files
USE mo_restart_attributes,   ONLY: t_RestartAttributeList, getAttributesForRestarting
USE mo_io_config,            ONLY: configure_io
USE mo_parallel_config,      ONLY: nproma, num_prefetch_proc
USE mo_nh_pzlev_config,      ONLY: configure_nh_pzlev
USE mo_advection_config,     ONLY: configure_advection
USE mo_art_config,           ONLY: configure_art
USE mo_assimilation_config,  ONLY: configure_lhn
USE mo_lnd_nwp_config,       ONLY: groups_smi
USE mo_run_config,           ONLY: dtime,                & !    namelist parameter
  &                                ltestcase,            &
  &                                ldynamics,            &
  &                                ltransport,           &
  &                                iforcing,             & !    namelist parameter
  &                                output_mode,          &
  &                                lvert_nest, ntracer,  &
  &                                ldass_lhn, msg_level, &
  &                                iqc, iqt,             &
  &                                ico2, io3,            &
  &                                number_of_grid_used
USE mo_nh_testcases,         ONLY: init_nh_testcase
USE mo_ls_forcing_nml,       ONLY: is_ls_forcing, is_nudging
USE mo_ls_forcing,           ONLY: init_ls_forcing
USE mo_dynamics_config,      ONLY: iequations, nnow, nnow_rcf, nnew, nnew_rcf, idiv_method, &
  &                                ldeepatmo
! Horizontal grid
USE mo_model_domain,         ONLY: p_patch
USE mo_grid_config,          ONLY: n_dom, start_time, end_time, &
     &                             is_plane_torus, l_limited_area, n_dom_start
USE mo_intp_data_strc,       ONLY: p_int_state
USE mo_intp_lonlat_types,    ONLY: lonlat_grids
USE mo_grf_intp_data_strc,   ONLY: p_grf_state
! Vertical grid
USE mo_vertical_grid,        ONLY: set_nh_metrics
! Grid nesting
USE mo_nh_nest_utilities,    ONLY: complete_nesting_setup
! NH-namelist state
USE mo_nonhydrostatic_config,ONLY: kstart_moist, kend_qvsubstep, l_open_ubc, &
  &                                itime_scheme, ndyn_substeps

USE mo_atm_phy_nwp_config,   ONLY: configure_atm_phy_nwp, atm_phy_nwp_config
USE mo_ensemble_pert_config, ONLY: configure_ensemble_pert, compute_ensemble_pert
USE mo_synsat_config,        ONLY: configure_synsat
! NH-Model states
USE mo_nonhydro_state,       ONLY: p_nh_state, p_nh_state_lists,               &
  &                                construct_nh_state, destruct_nh_state
USE mo_opt_diagnostics,      ONLY: construct_opt_diag, destruct_opt_diag,      &
  &                                compute_lonlat_area_weights
USE mo_nwp_phy_state,        ONLY: prm_diag, prm_nwp_tend,                     &
  &                                construct_nwp_phy_state, destruct_nwp_phy_state
USE mo_nwp_lnd_state,        ONLY: p_lnd_state, construct_nwp_lnd_state,       &
  &                                destruct_nwp_lnd_state
USE mo_interface_les,        ONLY: init_les_phy_interface
! Time integration
USE mo_nh_stepping,          ONLY: perform_nh_stepping
! Initialization with real data
USE mo_initicon,            ONLY: init_icon
USE mo_ext_data_state,      ONLY: ext_data
USE mo_ext_data_init,       ONLY: init_index_lists
! meteogram output
USE mo_meteogram_output,    ONLY: meteogram_init, meteogram_finalize
USE mo_meteogram_config,    ONLY: meteogram_output_config
USE mo_name_list_output_config,   ONLY: first_output_name_list, &
  &                               is_variable_in_output
USE mo_name_list_output_init, ONLY:  init_name_list_output,        &
  &                                  parse_variable_groups,        &
  &                                  collect_requested_ipz_levels, &
  &                                  output_file, create_vertical_axes
USE mo_level_selection,     ONLY: create_mipz_level_selections
USE mo_name_list_output,    ONLY: close_name_list_output
USE mo_pp_scheduler,        ONLY: pp_scheduler_init, pp_scheduler_finalize

! ECHAM physics
USE mo_echam_phy_memory,    ONLY: construct_echam_phy_state
USE mo_psrad_forcing_memory, ONLY: construct_psrad_forcing_list
USE mo_physical_constants,  ONLY: amd, amco2
USE mo_echam_phy_config,    ONLY: echam_phy_tc, dt_zero, echam_phy_config
USE mo_echam_rad_config,    ONLY: echam_rad_config
USE mo_echam_phy_init,      ONLY: init_echam_phy_params, init_echam_phy_external, &
   &                              init_echam_phy_field, init_o3_lcariolle
USE mo_echam_phy_cleanup,   ONLY: cleanup_echam_phy
#ifndef __NO_JSBACH__
  USE mo_jsb_model_init,    ONLY: jsbach_init_after_restart
#endif

USE mo_util_mtime,          ONLY: getElapsedSimTimeInSeconds
USE mo_output_event_types,  ONLY: t_sim_step_info
USE mo_action,              ONLY: ACTION_RESET, reset_act
USE mo_turbulent_diagnostic,ONLY: init_les_turbulent_output, close_les_turbulent_output
USE mo_limarea_config,      ONLY: latbc_config
USE mo_async_latbc_types,   ONLY: t_latbc_data
USE mo_async_latbc,         ONLY: init_prefetch, close_prefetch
USE mo_sync_latbc,          ONLY: deallocate_latbc_data
USE mo_radar_data_state,    ONLY: radar_data, init_radar_data, construct_lhn, lhn_fields, destruct_lhn
USE mo_rttov_interface,     ONLY: rttov_finalize, rttov_initialize
USE mo_synsat_config,       ONLY: lsynsat
USE mo_derived_variable_handling, ONLY: init_statistics_streams, finish_statistics_streams
USE mo_mpi,                 ONLY: my_process_is_stdio
USE mo_var_list,            ONLY: print_group_details
USE mo_sync,                ONLY: sync_patch_array, sync_c
USE mo_initicon_config,     ONLY: init_mode
USE mo_sleve_config,        ONLY: flat_height
USE mo_vertical_coord_table, ONLY: vct_a
USE mo_upatmo_config,       ONLY: configure_upatmo, destruct_upatmo

!-------------------------------------------------------------------------
#ifdef HAVE_CDI_PIO
  USE mo_impl_constants,      ONLY: pio_type_cdipio
  USE mo_parallel_config,     ONLY: pio_type
  USE mo_cdi,                 ONLY: namespaceGetActive, namespaceSetActive
  USE mo_cdi_pio_interface,         ONLY: nml_io_cdi_pio_namespace
#endif

IMPLICIT NONE
PRIVATE

PUBLIC :: atmo_nonhydrostatic
PUBLIC :: construct_atmo_nonhydrostatic, destruct_atmo_nonhydrostatic


CONTAINS

  !---------------------------------------------------------------------
  SUBROUTINE atmo_nonhydrostatic
    TYPE(t_latbc_data) :: latbc !< data structure for async latbc prefetching

!!$    CHARACTER(*), PARAMETER :: routine = "atmo_nonhydrostatic"

    CALL construct_atmo_nonhydrostatic(latbc)

    !------------------------------------------------------------------
    ! Now start the time stepping:
    ! The special initial time step for the three time level schemes
    ! is executed within process_grid_level
    !------------------------------------------------------------------

    CALL perform_nh_stepping( time_config%tc_current_date, latbc )

    !---------------------------------------------------------------------
    ! 6. Integration finished. Clean up.
    !---------------------------------------------------------------------
    CALL destruct_atmo_nonhydrostatic(latbc)

  END SUBROUTINE atmo_nonhydrostatic
  !---------------------------------------------------------------------

  !---------------------------------------------------------------------
  SUBROUTINE construct_atmo_nonhydrostatic(latbc)
    TYPE(t_latbc_data), INTENT(INOUT) :: latbc !< data structure for async latbc prefetching

    CHARACTER(*), PARAMETER :: routine = "construct_atmo_nonhydrostatic"

    INTEGER :: jg, jt, ist, jgroup
    LOGICAL :: l_pres_msl(n_dom) !< Flag. TRUE if computation of mean sea level pressure desired
    LOGICAL :: l_omega(n_dom)    !< Flag. TRUE if computation of vertical velocity desired
    LOGICAL :: l_rh(n_dom)       !< Flag. TRUE if computation of relative humidity desired
    LOGICAL :: l_pv(n_dom)       !< Flag. TRUE if computation of potential vorticity desired
    LOGICAL :: l_smi(n_dom)      !< Flag. TRUE if computation of soil moisture index desired
    TYPE(t_sim_step_info) :: sim_step_info  
    INTEGER :: jstep0
    INTEGER :: n_now, n_new, n_now_rcf, n_new_rcf
    REAL(wp) :: sim_time
    TYPE(t_RestartAttributeList), POINTER :: restartAttributes

    IF (timers_level > 3) CALL timer_start(timer_model_init)

    IF(iforcing == inwp) THEN

      CALL configure_ensemble_pert(ext_data, time_config%tc_exp_startdate)

      ! - generate index lists for tiles (land, ocean, lake)
      ! index lists for ice-covered and non-ice covered ocean points
      ! are initialized in init_nwp_phy
      CALL init_index_lists (p_patch(1:), ext_data)

      CALL configure_atm_phy_nwp(n_dom, p_patch(1:), dtime)

      CALL configure_synsat()

     ! initialize number of chemical tracers for convection
     DO jg = 1, n_dom
       CALL configure_art(jg)
     ENDDO

    ENDIF

    ! upper atmosphere
    CALL configure_upatmo( n_dom_start=n_dom_start, n_dom=n_dom,                 & 
      & p_patch=p_patch(n_dom_start:), ldeepatmo=ldeepatmo, init_mode=init_mode, &
      & iforcing=iforcing, dtime=dtime, ndyn_substeps=ndyn_substeps,             &
      & flat_height=flat_height, msg_level=msg_level, vct_a=vct_a                ) 

    ! initialize ldom_active flag if this is not a restart run

    ! calculate elapsed simulation time in seconds
    sim_time = getElapsedSimTimeInSeconds(time_config%tc_current_date) 

    DO jg=1, n_dom
      IF (jg > 1 .AND. start_time(jg) > sim_time .OR. end_time(jg) <= sim_time) THEN
        p_patch(jg)%ldom_active = .FALSE. ! domain not active
      ELSE
        p_patch(jg)%ldom_active = .TRUE.
      ENDIF
    ENDDO
    
    !---------------------------------------------------------------------
    ! 4.c Non-Hydrostatic / NWP
    !---------------------------------------------------------------------

    ALLOCATE (p_nh_state(n_dom), stat=ist)
    IF (ist /= success) THEN
      CALL finish(TRIM(routine),'allocation for p_nh_state failed')
    ENDIF

    ALLOCATE (p_nh_state_lists(n_dom), stat=ist)
    IF (ist /= success) THEN
      CALL finish(TRIM(routine),'allocation for p_nh_state_lists failed')
    ENDIF

    ! Note(GZ): Land state now needs to be allocated even if physics is turned
    ! off because ground temperature is included in feedback since r8133
    ! However, setting inwp_surface = 0 effects that only a few 2D fields are allocated
    ALLOCATE (p_lnd_state(n_dom), stat=ist)
    IF (ist /= success) THEN
      CALL finish(TRIM(routine),'allocation for p_lnd_state failed')
    ENDIF

    IF(iforcing /= inwp) atm_phy_nwp_config(:)%inwp_surface = 0



    ! Now allocate memory for the states
    DO jg=1,n_dom
      l_pres_msl(jg) = is_variable_in_output(first_output_name_list, var_name="pres_msl") .OR. &
        &              is_variable_in_output(first_output_name_list, var_name="psl_m")
      l_omega(jg)    = is_variable_in_output(first_output_name_list, var_name="omega")    .OR. &
        &              is_variable_in_output(first_output_name_list, var_name="wap_m")
    END DO
    CALL construct_nh_state(p_patch(1:), p_nh_state, p_nh_state_lists, n_timelevels=2, &
      &                     l_pres_msl=l_pres_msl, l_omega=l_omega)

    ! Add optional diagnostic variable lists (might remain empty)
    CALL construct_opt_diag(p_patch(1:), .TRUE.)

    IF(iforcing == inwp) THEN
      DO jg=1,n_dom
        l_rh(jg)  = is_variable_in_output(first_output_name_list, var_name="rh")
        l_pv(jg)  = is_variable_in_output(first_output_name_list, var_name="pv")
        l_smi(jg) = is_variable_in_output(first_output_name_list, var_name="smi")
        ! Check for special case: SMI is not in one of the output lists but it is part of a output group.
        ! In this case, the group can not be checked, as the connection between SMI and the group will be
        ! established during the add_var call. However, add_var for SMI will only be called if l_smi =.true.
        ! As a crutch, a character array containing the output groups of SMI from mo_lnd_nwp_config is used
        ! here and also at the add_var call.
        ! The check loops through the output groups. It has to be checked if l_smi is already .true., to not
        ! overwrite an existing .true. with a false. 
        IF(.not.l_smi(jg) ) THEN 
          ! Check for output groups containing SMI
          DO jgroup = 1,SIZE(groups_smi)
            IF(.not.l_smi(jg) ) THEN
              l_smi(jg) = is_variable_in_output(first_output_name_list, var_name='group:'//TRIM(groups_smi(jgroup)))
            END IF
          END DO
        END IF
      END DO
    END IF

    IF (iforcing == inwp) THEN
      CALL construct_nwp_phy_state( p_patch(1:), l_rh, l_pv )
      CALL construct_nwp_lnd_state( p_patch(1:), p_lnd_state, l_smi, n_timelevels=2 )
      CALL compute_ensemble_pert  ( p_patch(1:), ext_data, prm_diag, time_config%tc_current_date)
    END IF

#ifdef MESSY
    CALL messy_init_memory(n_dom)
#endif

    ! Due to the required ability to overwrite advection-Namelist settings
    ! via add_ref/add_tracer_ref for ICON-ART, configure_advection is called
    ! AFTER the nh_state is created. Otherwise, potential modifications of the
    ! advection-Namelist can not be taken into account properly.
    ! Unfortunately this conflicts with our trying to call the config-routines
    ! as early as possible.
    DO jg =1,n_dom
     CALL configure_advection( jg, p_patch(jg)%nlev, p_patch(1)%nlev,  &
       &                      iequations, iforcing, iqc, iqt,          &
       &                      kstart_moist(jg), kend_qvsubstep(jg),    &
       &                      lvert_nest, l_open_ubc, ntracer,         &
       &                      idiv_method, itime_scheme,               &
       &                      p_nh_state_lists(jg)%tracer_list(:)  )
    ENDDO

   IF (ldass_lhn) THEN 
     ALLOCATE (radar_data(n_dom), STAT=ist)
     IF (ist /= SUCCESS) THEN
          CALL finish(TRIM(routine),'allocation for radar_data failed')
     ENDIF
     ALLOCATE (lhn_fields(n_dom), STAT=ist)
     IF (ist /= SUCCESS) THEN
          CALL finish(TRIM(routine),'allocation for lhn_fields failed')
     ENDIF
     CALL message(TRIM(routine),'configure_lhn')
     DO jg =1,n_dom
       CALL configure_lhn(jg)
     ENDDO 

     CALL init_radar_data(p_patch(1:), radar_data)
     CALL construct_lhn(lhn_fields,p_patch(1:))
   ENDIF

    IF (iforcing == iecham) THEN
      CALL init_echam_phy_params       ( p_patch(1:) )
      CALL construct_echam_phy_state   ( p_patch(1:), ntracer )
      CALL construct_psrad_forcing_list( p_patch(1:) )
    END IF

    !------------------------------------------------------------------
    ! Prepare for time integration
    !------------------------------------------------------------------

    CALL set_nh_metrics(p_patch(1:)     ,&
         &              p_nh_state      ,&
         &              p_int_state(1:) ,&
         &              ext_data        )

    IF (n_dom > 1) THEN
      CALL complete_nesting_setup()
    END IF

    ! init LES
    DO jg = 1 , n_dom
      IF(atm_phy_nwp_config(jg)%is_les_phy) THEN
        CALL init_les_phy_interface(jg, p_patch(jg)       ,&
           &                        p_int_state(jg)       ,&
           &                        p_nh_state(jg)%metrics)
      END IF
    END DO

    !------------------------------------------------------------------
    ! Prepare initial conditions for time integration.
    !------------------------------------------------------------------
    !
    IF (isRestart()) THEN
      !
      ! This is a resumed integration. Read model state from restart file(s).
      !
      IF (timers_level > 5) CALL timer_start(timer_read_restart)
      !
      DO jg = 1,n_dom
        IF (p_patch(jg)%ldom_active) THEN
          CALL read_restart_files( p_patch(jg), n_dom)
        END IF
      END DO
      !
      CALL message(TRIM(routine),'normal exit from read_restart_files')
      !
      IF (timers_level > 5) CALL timer_stop(timer_read_restart)
      !
#ifndef __NO_JSBACH__
      DO jg = 1,n_dom
        IF (.NOT. p_patch(jg)%ldom_active) CYCLE
        IF (echam_phy_config(jg)%ljsb) THEN
          CALL jsbach_init_after_restart(jg)
        END IF
      END DO
#endif
      !
    ELSE
      !
      ! This is a new integration.
      !
      IF (ltestcase) THEN
        !
        ! Initialize testcase analytically
        !
        CALL init_nh_testcase(p_patch(1:)     ,&
          &                   p_nh_state      ,&
          &                   p_int_state(1:) ,&
          &                   p_lnd_state(1:) ,&
          &                   ext_data        ,&
          &                   ntl=2           )
        !
        IF(is_ls_forcing .OR. is_nudging) &
          CALL init_ls_forcing(p_nh_state(1)%metrics)
        !
      ELSE
        !
        ! Initialize with real atmospheric data
        !
        IF (timers_level > 5) CALL timer_start(timer_init_icon)
        !
        IF (iforcing == inwp) THEN
          !
          ! Initialize atmosphere, surface and land
          !
          CALL init_icon (p_patch(1:)     ,&
            &             p_int_state(1:) ,&
            &             p_grf_state(1:) ,&
            &             p_nh_state(1:)  ,&
            &             ext_data(1:)    ,&
            &             prm_diag(1:)    ,&
            &             p_lnd_state(1:) )
          !
        ELSE ! iforcing == iecham, inoforcing, ...
          !
          ! Initialize the atmosphere only
          !
          CALL init_icon (p_patch(1:)     ,&
            &             p_int_state(1:) ,&
            &             p_grf_state(1:) ,&
            &             p_nh_state(1:)  ,&
            &             ext_data(1:)    )
          !
        END IF ! iforcing
        !
        IF (timers_level > 5) CALL timer_stop(timer_init_icon)
        !
      END IF ! ltestcase
      !
      ! Initialize tracers fields jt=iqt to jt=ntracer, which are not available in the analysis file,
      ! but may be used with ECHAM physics, for real cases or test cases.
      !
      IF (iforcing == iecham ) THEN
        DO jg = 1,n_dom
          IF (.NOT. p_patch(jg)%ldom_active) CYCLE
          !
          ! CO2 tracer
          IF ( iqt <= ico2 .AND. ico2 <= ntracer) THEN
!$OMP PARALLEL
            CALL init(p_nh_state(jg)%prog(nnow_rcf(jg))%tracer(:,:,:,ico2),echam_rad_config(jg)% vmr_co2*amco2/amd)
!$OMP END PARALLEL
            CALL print_value(TRIM(routine)//': CO2 tracer initialized with constant vmr', &
              &              echam_rad_config(jg)% vmr_co2*amco2/amd)
          END IF
          !
          ! O3 tracer
          IF ( iqt <= io3 .AND. io3 <= ntracer) THEN
            IF (echam_phy_tc(jg)%dt_car > dt_zero) THEN
              CALL init_o3_lcariolle( time_config%tc_current_date                          ,&
                &                     p_patch(jg)                                          ,&
                &                     p_nh_state(jg)%diag% pres               (:,:,:)      ,&
                &                     p_nh_state(jg)%prog(nnow_rcf(jg))%tracer(:,:,:,io3)  )
              CALL sync_patch_array ( sync_c,p_patch(jg)                                   ,&
                &                     p_nh_state(jg)%prog(nnow_rcf(jg))%tracer(:,:,:,io3)  )
              CALL message(TRIM(routine),'o3 tracer is initialized by the Cariolle lin. o3 scheme')
            ELSE
!$OMP PARALLEL
              CALL init(p_nh_state(jg)%prog(nnow_rcf(jg))%tracer(:,:,:,io3),0.0_wp)
!$OMP END PARALLEL
              CALL message(TRIM(routine),'o3 tracer is initialized to zero, check setup')
            END IF
          END IF
          !
        END DO
        !
      END IF
      !
    END IF ! isRestart()


    ! Now set up ECHAM physics fields
    !
    IF ( iforcing == iecham ) THEN
      !
      ! read external data for real case
      IF (.NOT. ltestcase) THEN 
        CALL init_echam_phy_external( p_patch(1:)                 ,&
           &                          time_config%tc_current_date )
      END IF
      !
      ! prepare fields of the physics state, real and test case
      DO jg = 1,n_dom
        CALL init_echam_phy_field( p_patch(jg)                                        ,&
          &                        ext_data  (jg)% atm%topography_c       (:,  :)     ,&
          &                        p_nh_state(jg)% metrics% z_ifc         (:,:,:)     ,&
          &                        p_nh_state(jg)% metrics% z_mc          (:,:,:)     ,&
          &                        p_nh_state(jg)% metrics% ddqz_z_full   (:,:,:)     ,&
          &                        p_nh_state(jg)% metrics% geopot_agl_ifc(:,:,:)     ,&
          &                        p_nh_state(jg)% metrics% geopot_agl    (:,:,:)     ,&
          &                        p_nh_state(jg)% diag% temp             (:,:,:)     )
      END DO
      !
    END IF


    ! Copy prognostic variables, for which no tendencies are computed,
    ! from time level "now" to time level "new", so that they are correctly set
    ! at odd and even time steps.
    !
    IF (.NOT.ldynamics)  THEN ! copy prognostic variables of the dynamics
      DO jg = 1,n_dom
         IF (.NOT. p_patch(jg)%ldom_active) CYCLE
         n_now = nnow(jg)
         n_new = nnew(jg)
!$OMP PARALLEL
         CALL copy(p_nh_state(jg)%prog(n_now)%vn     , p_nh_state(jg)%prog(n_new)%vn     )
         CALL copy(p_nh_state(jg)%prog(n_now)%w      , p_nh_state(jg)%prog(n_new)%w      )
         CALL copy(p_nh_state(jg)%prog(n_now)%theta_v, p_nh_state(jg)%prog(n_new)%theta_v)
         CALL copy(p_nh_state(jg)%prog(n_now)%exner  , p_nh_state(jg)%prog(n_new)%exner  )
         CALL copy(p_nh_state(jg)%prog(n_now)%rho    , p_nh_state(jg)%prog(n_new)%rho    )
!$OMP END PARALLEL
      END DO
    END IF

    IF (.NOT.ltransport) THEN ! copy prognostic variables of the transport
      DO jg = 1,n_dom
         IF (.NOT. p_patch(jg)%ldom_active) CYCLE
         n_now_rcf = nnow_rcf(jg)
         n_new_rcf = nnew_rcf(jg)
         DO jt = 1,ntracer
!$OMP PARALLEL
            CALL copy(p_nh_state(jg)%prog(n_now_rcf)%tracer(:,:,:,jt), p_nh_state(jg)%prog(n_new_rcf)%tracer(:,:,:,jt) )
!$OMP END PARALLEL
         END DO
      END DO
    END IF


    !------------------------------------------------------------------
    ! Asynchronous pre-fetching
    !------------------------------------------------------------------

    ! If async prefetching is in effect, init_prefetch is a collective call
    ! with the prefetching processor and effectively starts async prefetching
    IF ((num_prefetch_proc == 1) .AND. (latbc_config%itype_latbc > 0)) THEN
      IF (timers_level > 5) CALL timer_start(timer_init_latbc)
      CALL init_prefetch(latbc)
      IF (timers_level > 5) CALL timer_stop(timer_init_latbc)
    ENDIF

    !------------------------------------------------------------------
    ! Prepare output file
    !------------------------------------------------------------------

    CALL configure_io()   ! set n_chkpt and n_diag, which control
                          ! writing of restart files and tot_int diagnostics.

    ! Add a special metrics variable containing the area weights of
    ! the regular lon-lat grid.
    CALL compute_lonlat_area_weights(lonlat_grids, p_nh_state_lists)

    ! Map the variable groups given in the output namelist onto the
    ! corresponding variable subsets:
    ! ATTENTION: all add_vars must be finished before calling this routine.
    IF (output_mode%l_nml) THEN
      CALL parse_variable_groups()
    END IF

    ! if output on z and/or p-levels is required do some config
    !
    ! Note that on the compute PEs we must call this *before* the
    ! initialization of the vertical interpolation,
    ! "pp_scheduler_init"
    IF (output_mode%l_nml) THEN
      CALL collect_requested_ipz_levels()
      DO jg = 1, n_dom
        CALL configure_nh_pzlev(jg, nproma, p_patch(jg)%npromz_c,  &
          &                     p_patch(jg)%nblks_c)
      ENDDO
    END IF

    ! setup of post-processing job queue, e.g. setup of optional
    ! diagnostic quantities like pz-level interpolation
    CALL pp_scheduler_init( (iforcing == inwp) )

    ! setup of RTTOV interface (assumes expanded variable groups)
    IF (ANY(lsynsat(:)))  CALL rttov_initialize()

    ! If async IO is in effect, init_name_list_output is a collective call
    ! with the IO procs and effectively starts async IO
    IF (output_mode%l_nml) THEN
      ! compute sim_start, sim_end
      CALL datetimeToString(time_config%tc_exp_startdate, sim_step_info%sim_start)
      CALL datetimeToString(time_config%tc_exp_stopdate, sim_step_info%sim_end)
      CALL datetimeToString(time_config%tc_startdate, sim_step_info%run_start)
      CALL datetimeToString(time_config%tc_stopdate, sim_step_info%restart_time)

      sim_step_info%dtime      = dtime
      jstep0 = 0

      restartAttributes => getAttributesForRestarting()
      IF (ASSOCIATED(restartAttributes)) THEN
        ! get start counter for time loop from restart file:
        jstep0 = restartAttributes%getInteger("jstep")
      END IF
      sim_step_info%jstep0    = jstep0
      CALL init_statistics_streams
      CALL init_name_list_output(sim_step_info)

      !---------------------------------------------------------------------
      !     Setup of meteogram output
      !---------------------------------------------------------------------
      DO jg =1,n_dom
        IF (meteogram_output_config(jg)%lenabled) THEN
          ! For dry test cases: do not sample moist variables
          ! (but allow for TORUS moist runs; see also mo_mtgrm_output.F90)
          IF (ltestcase .and. .not. is_plane_torus) THEN
            CALL meteogram_init(meteogram_output_config(jg), jg, p_patch(jg), &
              &                ext_data(jg), p_nh_state(jg),                  &
              &                p_lnd_state=p_lnd_state(jg), iforcing=iforcing,&
              &                grid_uuid=p_patch(jg)%grid_uuid,               &
              &                number_of_grid_used=number_of_grid_used(jg) )
          ELSE
            CALL meteogram_init(meteogram_output_config(jg), jg, p_patch(jg), &
              &                ext_data(jg), p_nh_state(jg), prm_diag(jg),    &
              &                p_lnd_state(jg), prm_nwp_tend(jg), iforcing,                     &
              &                grid_uuid=p_patch(jg)%grid_uuid,               &
              &                number_of_grid_used=number_of_grid_used(jg) )
          END IF
        END IF
      END DO

      CALL create_mipz_level_selections(output_file)
      CALL create_vertical_axes(output_file)
    END IF

#ifdef MESSY
    CALL messy_init_coupling
    CALL messy_init_tracer
#endif

    ! Determine if temporally averaged vertically integrated moisture quantities need to be computed

    IF (iforcing == inwp) THEN
        atm_phy_nwp_config(1:n_dom)%lcalc_moist_integral_avg = &
        is_variable_in_output(first_output_name_list, var_name="clct_avg")        .OR. &
        is_variable_in_output(first_output_name_list, var_name="tracer_vi_avg01") .OR. &
        is_variable_in_output(first_output_name_list, var_name="tracer_vi_avg02") .OR. &
        is_variable_in_output(first_output_name_list, var_name="tracer_vi_avg03") .OR. &
        is_variable_in_output(first_output_name_list, var_name="avg_qv")          .OR. &
        is_variable_in_output(first_output_name_list, var_name="avg_qc")          .OR. &
        is_variable_in_output(first_output_name_list, var_name="avg_qi")

        atm_phy_nwp_config(1:n_dom)%lcalc_extra_avg = &
        is_variable_in_output(first_output_name_list, var_name="astr_u_sso")      .OR. &
        is_variable_in_output(first_output_name_list, var_name="accstr_u_sso")    .OR. &
        is_variable_in_output(first_output_name_list, var_name="astr_v_sso")      .OR. &
        is_variable_in_output(first_output_name_list, var_name="accstr_v_sso")    .OR. &
        is_variable_in_output(first_output_name_list, var_name="adrag_u_grid")    .OR. &
        is_variable_in_output(first_output_name_list, var_name="adrag_v_grid")
     ENDIF

    !Anurag Dipankar, MPIM (2015-08-01): always call this routine
    !for LES simulation
    DO jg = 1 , n_dom
      atm_phy_nwp_config(jg)%lcalc_moist_integral_avg &
           = atm_phy_nwp_config(jg)%lcalc_moist_integral_avg &
           .OR. atm_phy_nwp_config(jg)%is_les_phy
    END DO

    !----------------------!
    !  Initialize actions  !
    !----------------------!
    !

    ! Initialize reset-Action, i.e. assign variables to action object
    CALL reset_act%initialize(ACTION_RESET)


    !Anurag Dipankar, MPIM (2014-01-14)
    !Special 1D and 0D output for LES runs till we get add_var/nml_out working
    !Only for Torus runs with single domain
    DO jg = 1 , n_dom

      IF(atm_phy_nwp_config(jg)%is_les_phy .AND. is_plane_torus) &
           CALL init_les_turbulent_output(p_patch(jg), p_nh_state(jg)%metrics, &
           &                              time_config%tc_startdate, l_rh(jg), ldelete=(.NOT. isRestart()))

    END DO


    !-------------------------------------------------------!
    !  (Optional) detailed print-out of some variable info  !
    !-------------------------------------------------------!

    ! variable group information
    IF (my_process_is_stdio() .AND. (msg_level >= 15)) THEN
      CALL print_group_details(idom=1,                            &
        &                      opt_latex_fmt           = .TRUE., &
        &                      opt_reduce_trailing_num = .TRUE.,  &
        &                      opt_skip_trivial        = .TRUE.)
    END IF

    IF (timers_level > 3) CALL timer_stop(timer_model_init)

  END SUBROUTINE construct_atmo_nonhydrostatic

  !---------------------------------------------------------------------
  SUBROUTINE destruct_atmo_nonhydrostatic(latbc)
    TYPE(t_latbc_data), INTENT(INOUT) :: latbc !< data structure for async latbc prefetching

    CHARACTER(*), PARAMETER :: routine = "destruct_atmo_nonhydrostatic"


    INTEGER :: jg, ist
    
#ifdef HAVE_CDI_PIO
    INTEGER :: prev_cdi_namespace
#endif

    !---------------------------------------------------------------------
    ! 6. Integration finished. Clean up.
    !---------------------------------------------------------------------

    CALL message(TRIM(routine),'start to clean up')

#ifdef MESSY
    CALL messy_free_memory
#endif

    ! Destruction of post-processing job queue
    CALL pp_scheduler_finalize()

    ! Destruction of some RTTOV data structures  (if enabled)
    IF (ANY(lsynsat(:)))  CALL rttov_finalize()

    ! Delete optional diagnostics
    CALL destruct_opt_diag()

    ! Delete state variables

    CALL destruct_nh_state( p_nh_state, p_nh_state_lists )
    DEALLOCATE (p_nh_state, STAT=ist)
    IF (ist /= SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for p_nh_state failed')
    ENDIF
    DEALLOCATE (p_nh_state_lists, STAT=ist)
    IF (ist /= SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for p_nh_state_lists failed')
    ENDIF

    IF (iforcing == inwp) THEN
      CALL destruct_nwp_phy_state
      CALL destruct_nwp_lnd_state( p_lnd_state )
      DO jg = 1, n_dom
        CALL atm_phy_nwp_config(jg)%finalize()
      ENDDO
    ENDIF

    IF (iforcing == iecham) THEN
      CALL cleanup_echam_phy
    ENDIF

    CALL destruct_upatmo()

    ! call close name list prefetch
    IF (l_limited_area .AND. latbc_config%itype_latbc > 0) THEN
      IF (num_prefetch_proc >= 1) THEN
        CALL close_prefetch()
        CALL latbc%finalize()
      ELSE
        CALL deallocate_latbc_data()
      END IF
    END IF

    ! Delete output variable lists
    IF (output_mode%l_nml) THEN
      CALL close_name_list_output
      CALL finish_statistics_streams
    END IF
#ifdef HAVE_CDI_PIO
    IF (pio_type == pio_type_cdipio) THEN
      prev_cdi_namespace = namespaceGetActive()
      CALL namespaceSetActive(nml_io_cdi_pio_namespace)
      CALL pioFinalize
      CALL namespaceSetActive(prev_cdi_namespace)
    END IF
#endif
    ! finalize meteogram output
    IF (output_mode%l_nml) THEN
      DO jg = 1, n_dom
        IF (meteogram_output_config(jg)%lenabled) THEN
          CALL meteogram_finalize(jg)
        END IF
      END DO
      DO jg = 1, max_dom
        DEALLOCATE(meteogram_output_config(jg)%station_list)
      END DO
    END IF

    !Close LES diag files
    DO jg = 1 , n_dom
     IF(atm_phy_nwp_config(jg)%is_les_phy .AND. is_plane_torus) &
       CALL close_les_turbulent_output(jg)
    END DO

    IF (ldass_lhn) THEN 
      ! deallocate ext_data array
      DEALLOCATE(radar_data, STAT=ist)
      IF (ist /= SUCCESS) THEN
        CALL finish(TRIM(routine), 'deallocation of radar_data')
      ENDIF
      CALL destruct_lhn (lhn_fields)
    ENDIF
 
    CALL message(TRIM(routine),'clean-up finished')

  END SUBROUTINE destruct_atmo_nonhydrostatic
  !---------------------------------------------------------------------

END MODULE mo_atmo_nonhydrostatic

