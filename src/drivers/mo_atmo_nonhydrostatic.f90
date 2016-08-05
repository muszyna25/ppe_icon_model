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
USE mo_exception,            ONLY: message, finish
USE mo_fortran_tools,        ONLY: copy
USE mo_impl_constants,       ONLY: SUCCESS, max_dom, inwp, iecham
USE mo_timer,                ONLY: timers_level, timer_start, timer_stop, &
  &                                timer_model_init, timer_init_icon, timer_read_restart
USE mo_master_config,        ONLY: isRestart
USE mo_time_config,          ONLY: time_config      ! variable
USE mo_restart,              ONLY: read_restart_files
USE mo_restart_attributes,   ONLY: t_RestartAttributeList, getAttributesForRestarting
USE mo_io_config,            ONLY: configure_io
USE mo_parallel_config,      ONLY: nproma, num_prefetch_proc
USE mo_nh_pzlev_config,      ONLY: configure_nh_pzlev
USE mo_advection_config,     ONLY: configure_advection
USE mo_art_config,           ONLY: configure_art
USE mo_run_config,           ONLY: dtime,                & !    namelist parameter
  &                                ltestcase,            &
  &                                ldynamics,            &
  &                                ltransport,           &
  &                                iforcing,             & !    namelist parameter
  &                                output_mode,          &
  &                                lvert_nest, ntracer,  &
  &                                nlev,                 &
  &                                iqv, iqc, iqt,        &
  &                                number_of_grid_used
USE mo_dynamics_config,      ONLY: iequations, nnow, nnow_rcf, nnew, nnew_rcf, idiv_method
! Horizontal grid
USE mo_model_domain,         ONLY: p_patch
USE mo_grid_config,          ONLY: n_dom, start_time, end_time, is_plane_torus
USE mo_intp_data_strc,       ONLY: p_int_state
USE mo_grf_intp_data_strc,   ONLY: p_grf_state
! NH-namelist state
USE mo_nonhydrostatic_config,ONLY: kstart_moist, kend_qvsubstep, l_open_ubc, &
  &                                itime_scheme

USE mo_atm_phy_nwp_config,   ONLY: configure_atm_phy_nwp, atm_phy_nwp_config
USE mo_ensemble_pert_config, ONLY: configure_ensemble_pert, compute_ensemble_pert
USE mo_synsat_config,        ONLY: configure_synsat
! NH-Model states
USE mo_nonhydro_state,       ONLY: p_nh_state, p_nh_state_lists,               &
  &                                construct_nh_state, destruct_nh_state
USE mo_opt_diagnostics,      ONLY: construct_opt_diag, destruct_opt_diag
USE mo_nwp_phy_state,        ONLY: prm_diag, construct_nwp_phy_state,          &
  &                                destruct_nwp_phy_state
USE mo_nwp_lnd_state,        ONLY: p_lnd_state, construct_nwp_lnd_state,       &
  &                                destruct_nwp_lnd_state
! Time integration
USE mo_nh_stepping,          ONLY: prepare_nh_integration, perform_nh_stepping
! Initialization with real data
USE mo_initicon,            ONLY: init_icon
USE mo_initicon_config,     ONLY: timeshift
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
  &                                  output_file
USE mo_name_list_output_zaxes, ONLY: create_mipz_level_selections
USE mo_name_list_output,    ONLY: close_name_list_output
USE mo_pp_scheduler,        ONLY: pp_scheduler_init, pp_scheduler_finalize
USE mo_intp_lonlat,         ONLY: compute_lonlat_area_weights

! LGS - for the implementation of ECHAM physics in iconam
USE mo_echam_phy_init,      ONLY: init_echam_phy, initcond_echam_phy
USE mo_echam_phy_cleanup,   ONLY: cleanup_echam_phy
USE mo_vertical_coord_table,ONLY: vct_a, vct_b
USE mo_nh_testcases_nml,    ONLY: nh_test_name

USE mo_master_config,       ONLY: tc_exp_startdate, tc_exp_stopdate, tc_startdate, tc_stopdate
USE mtime,                  ONLY: datetimeToString
USE mo_mtime_extensions,    ONLY: get_datetime_string
USE mo_output_event_types,  ONLY: t_sim_step_info
USE mo_action,              ONLY: ACTION_RESET, reset_act
USE mo_turbulent_diagnostic,ONLY: init_les_turbulent_output, close_les_turbulent_output
USE mo_limarea_config,      ONLY: latbc_config
USE mo_async_latbc,         ONLY: init_prefetch, close_prefetch

USE mo_rttov_interface,     ONLY: rttov_finalize, rttov_initialize
USE mo_synsat_config,       ONLY: lsynsat
USE mo_derived_variable_handling, ONLY: init_mean_stream, finish_mean_stream
!-------------------------------------------------------------------------

IMPLICIT NONE
PRIVATE

PUBLIC :: atmo_nonhydrostatic
PUBLIC :: construct_atmo_nonhydrostatic, destruct_atmo_nonhydrostatic


CONTAINS

  !---------------------------------------------------------------------
  SUBROUTINE atmo_nonhydrostatic

!!$    CHARACTER(*), PARAMETER :: routine = "mo_atmo_nonhydrostatic"

    CALL construct_atmo_nonhydrostatic()

    !------------------------------------------------------------------
    ! Now start the time stepping:
    ! The special initial time step for the three time level schemes
    ! is executed within process_grid_level
    !------------------------------------------------------------------

    CALL perform_nh_stepping( time_config%cur_datetime )

    !---------------------------------------------------------------------
    ! 6. Integration finished. Clean up.
    !---------------------------------------------------------------------
    CALL destruct_atmo_nonhydrostatic()

  END SUBROUTINE atmo_nonhydrostatic
  !---------------------------------------------------------------------

  !---------------------------------------------------------------------
  SUBROUTINE construct_atmo_nonhydrostatic

    CHARACTER(*), PARAMETER :: routine = "construct_atmo_nonhydrostatic"

    INTEGER :: jg, jt, ist
    LOGICAL :: l_pres_msl(n_dom) !< Flag. TRUE if computation of mean sea level pressure desired
    LOGICAL :: l_omega(n_dom)    !< Flag. TRUE if computation of vertical velocity desired
    LOGICAL :: l_rh(n_dom)       !< Flag. TRUE if computation of relative humidity desired
    LOGICAL :: l_pv(n_dom)       !< Flag. TRUE if computation of potential vorticity desired
    TYPE(t_sim_step_info) :: sim_step_info  
    INTEGER :: jstep0
    INTEGER :: n_now, n_new, n_now_rcf, n_new_rcf
    REAL(wp) :: sim_time
    TYPE(t_RestartAttributeList), POINTER :: restartAttributes

    IF (timers_level > 3) CALL timer_start(timer_model_init)

    IF(iforcing == inwp) THEN

      CALL configure_ensemble_pert(ext_data)

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

    ! initialize ldom_active flag if this is not a restart run
    restartAttributes => getAttributesForRestarting()
    IF (.NOT. ASSOCIATED(restartAttributes)) THEN
      DO jg=1, n_dom
        IF (jg > 1 .AND. start_time(jg) - timeshift%dt_shift > 0._wp) THEN
          p_patch(jg)%ldom_active = .FALSE. ! domain not active from the beginning
        ELSE
          p_patch(jg)%ldom_active = .TRUE.
        ENDIF
      ENDDO
    ELSE
      sim_time = restartAttributes%getReal("sim_time_DOM01")
      DO jg=1, n_dom
        IF (jg > 1 .AND. start_time(jg) > sim_time .OR. end_time(jg) <= sim_time) THEN
          p_patch(jg)%ldom_active = .FALSE. ! domain not active at restart time
        ELSE
          p_patch(jg)%ldom_active = .TRUE.
        ENDIF
      ENDDO
    ENDIF

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
        l_rh(jg) = is_variable_in_output(first_output_name_list, var_name="rh")
        l_pv(jg) = is_variable_in_output(first_output_name_list, var_name="pv")
      END DO
    END IF

    IF (iforcing == inwp) THEN
      CALL construct_nwp_phy_state( p_patch(1:), l_rh, l_pv )
      CALL construct_nwp_lnd_state( p_patch(1:),p_lnd_state,n_timelevels=2 )
      CALL compute_ensemble_pert  ( p_patch(1:), ext_data, prm_diag)
    END IF

#ifdef MESSY
    CALL messy_init_memory(n_dom)
#endif

    ! Due to the required ability to overwrite advection-Namelist settings
    ! via add_ref/add_tracer_ref for ICON-ART, configure_advection is called
    ! AFTER the nh_state is created. Otherwise, potential modifications of the
    ! advection-Namelist can not be taken into account properly.
    ! Unfortunatley this conflicts with our trying to call the config-routines
    ! as early as possible.
    DO jg =1,n_dom
     CALL configure_advection( jg, p_patch(jg)%nlev, p_patch(1)%nlev,  &
       &                      iequations, iforcing, iqc, iqt,          &
       &                      kstart_moist(jg), kend_qvsubstep(jg),    &
       &                      lvert_nest, l_open_ubc, ntracer,         &
       &                      idiv_method, itime_scheme,               &
       &                      p_nh_state_lists(jg)%tracer_list(:)  )
    ENDDO


    !---------------------------------------------------------------------
    ! 5. Perform time stepping
    !---------------------------------------------------------------------
      !------------------------------------------------------------------
      ! Prepare for time integration
      !------------------------------------------------------------------


    ! CALL prepare_nh_integration(p_patch(1:), p_nh_state, p_int_state(1:), p_grf_state(1:))

    CALL prepare_nh_integration( )

! LGS
!
    IF ( iforcing == iecham ) THEN
      CALL init_echam_phy( p_patch(1:), nh_test_name, &
        & nlev, vct_a, vct_b, time_config%cur_datetime )
      !! many of the initial conditions for the echam 'field' are set here
      !! Note: it is not certain that p_nh_state(jg)%diag%temp has been initialized at this point in time.
      !!       initcond_echam_phy should therefore not rely on the fact that this has been properly set.
      !!       It is the case for some testcases, but not e.g. for coupled or AMIP runs if the atmosphere
      !!       is initialized with IFS analyses.
      DO jg = 1,n_dom
        CALL initcond_echam_phy( jg                                               ,&
          &                      p_patch(jg)                                      ,&
          &                      p_nh_state(jg)%diag          % temp  (:,:,:)     ,&
          &                      p_nh_state(jg)%prog(nnow(jg))% tracer(:,:,:,iqv) ,&
          &                      nh_test_name                                     )
      END DO
    END IF
!---

    !
    ! Read restart files (if necessary)
    !
    IF (isRestart()) THEN
      ! This is a resumed integration. Read model state from restart file(s).

      IF (timers_level > 5) CALL timer_start(timer_read_restart)
#ifdef NOMPI
      ! TODO : Non-MPI mode does not work for multiple domains
      DO jg = 1,n_dom
        IF (p_patch(jg)%ldom_active) THEN
          CALL read_restart_files( p_patch(jg), n_dom)
        ENDIF
      END DO
#else
      DO jg = 1,n_dom
        IF (p_patch(jg)%ldom_active) THEN
          CALL read_restart_files( p_patch(jg), n_dom )
        ENDIF
      END DO
#endif
      CALL message(TRIM(routine),'normal exit from read_restart_files')
      IF (timers_level > 5) CALL timer_stop(timer_read_restart)

    ENDIF



    !------------------------------------------------------------------
    ! Prepare initial conditions for time integration.
    !------------------------------------------------------------------
    !
    ! Initialize model with real atmospheric data if appropriate switches are set
    !
    IF (.NOT. ltestcase .AND. .NOT. isRestart() ) THEN

      IF (iforcing == inwp) THEN

        ! Initialize atmosphere, surface and land

        IF (timers_level > 5) CALL timer_start(timer_init_icon)
        CALL init_icon (p_patch(1:)     ,&
          &             p_int_state(1:) ,&
          &             p_grf_state(1:) ,&
          &             p_nh_state(1:)  ,&
          &             ext_data(1:)    ,&
          &             prm_diag(1:)    ,&
          &             p_lnd_state(1:) )
        IF (timers_level > 5) CALL timer_stop(timer_init_icon)

      ELSE

        ! Initialize the atmosphere only

        IF (timers_level > 5) CALL timer_start(timer_init_icon)
        CALL init_icon (p_patch(1:)     ,&
          &             p_int_state(1:) ,&
          &             p_grf_state(1:) ,&
          &             p_nh_state(1:)  ,&
          &             ext_data(1:)    )
        IF (timers_level > 5) CALL timer_stop(timer_init_icon)

      END IF

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

    ! If async prefetching is in effect, init_prefetch is a collective call
    ! with the prefetching processor and effectively starts async prefetching
    IF ((num_prefetch_proc == 1) .AND. (latbc_config%itype_latbc > 0)) &
       CALL init_prefetch

    !------------------------------------------------------------------
    ! Prepare output file
    !------------------------------------------------------------------

    CALL configure_io()   ! set n_chkpt and n_diag, which control
                          ! writing of restart files and tot_int diagnostics.

    ! Add a special metrics variable containing the area weights of
    ! the regular lon-lat grid.
    CALL compute_lonlat_area_weights()

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
#ifdef USE_MTIME_LOOP
      CALL datetimeToString(tc_exp_startdate, sim_step_info%sim_start)
      CALL datetimeToString(tc_exp_stopdate, sim_step_info%sim_end)
      CALL datetimeToString(tc_startdate, sim_step_info%run_start)
      CALL datetimeToString(tc_stopdate, sim_step_info%restart_time)
#else
      CALL get_datetime_string(sim_step_info%sim_start, time_config%ini_datetime)
      CALL get_datetime_string(sim_step_info%sim_end,   time_config%end_datetime)
      CALL get_datetime_string(sim_step_info%restart_time,  time_config%cur_datetime, &
        &                      INT(time_config%dt_restart))
      CALL get_datetime_string(sim_step_info%run_start, time_config%cur_datetime)
#endif
      sim_step_info%dtime      = dtime
      jstep0 = 0
      IF (ASSOCIATED(restartAttributes) .AND. .NOT. time_config%is_relative_time) THEN
        ! get start counter for time loop from restart file:
        jstep0 = restartAttributes%getInteger("jstep")
      END IF
      sim_step_info%jstep0    = jstep0
      CALL init_mean_stream(p_patch(1))
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
              &                p_lnd_state(jg), iforcing,                     &
              &                grid_uuid=p_patch(jg)%grid_uuid,               &
              &                number_of_grid_used=number_of_grid_used(jg) )
          END IF
        END IF
      END DO

      CALL create_mipz_level_selections(output_file)
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
           time_config%sim_time(jg), l_rh(jg), ldelete=(.NOT. isRestart()))

    END DO

    IF (timers_level > 3) CALL timer_stop(timer_model_init)

  END SUBROUTINE construct_atmo_nonhydrostatic

  !---------------------------------------------------------------------
  SUBROUTINE destruct_atmo_nonhydrostatic

    CHARACTER(*), PARAMETER :: routine = "destruct_atmo_nonhydrostatic"


    INTEGER :: jg, ist

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
    ENDIF

    IF (iforcing == iecham) THEN
      CALL cleanup_echam_phy
    ENDIF

    ! call close name list prefetch
    IF((num_prefetch_proc == 1) .AND. (latbc_config%itype_latbc > 0)) &
       CALL close_prefetch

    ! Delete output variable lists
    IF (output_mode%l_nml) THEN
      CALL close_name_list_output
      CALL finish_mean_stream()
    END IF

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

    CALL message(TRIM(routine),'clean-up finished')

  END SUBROUTINE destruct_atmo_nonhydrostatic
  !---------------------------------------------------------------------

END MODULE mo_atmo_nonhydrostatic

