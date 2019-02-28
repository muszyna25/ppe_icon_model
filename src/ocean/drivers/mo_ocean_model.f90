!>
!! @page pagecontrolmodelf901 Main program for the ICON ocean model
!!
!! @par Revision History
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_ocean_model

  USE mo_exception,           ONLY: message, finish
  USE mo_master_config,       ONLY: isRestart
  USE mo_parallel_config,     ONLY: p_test_run, l_test_openmp, num_io_procs, &
       &                            pio_type, num_test_pe
  USE mo_mpi,                 ONLY: set_mpi_work_communicators, process_mpi_io_size, &
       &                            stop_mpi, my_process_is_io, my_process_is_mpi_test,   &
       &                            set_mpi_work_communicators, process_mpi_io_size
  USE mo_timer,               ONLY: init_timer, timer_start, timer_stop, print_timer, &
       &                            timer_model_init
  USE mo_memory_log,              ONLY: memory_log_terminate
  USE mtime,                  ONLY: MAX_DATETIME_STR_LEN, datetimeToString
  USE mo_name_list_output_init, ONLY: init_name_list_output, parse_variable_groups, &
    &                                 create_vertical_axes, output_file
  USE mo_derived_variable_handling, ONLY: init_mean_stream, finish_mean_stream
  USE mo_name_list_output,    ONLY: close_name_list_output, name_list_io_main_proc
  USE mo_name_list_output_config,  ONLY: use_async_name_list_io
  USE mo_level_selection, ONLY: create_mipz_level_selections
  USE mo_dynamics_config,     ONLY: configure_dynamics
  USE mo_zaxis_type,          ONLY: zaxisTypeList, t_zaxisTypeList

  !  USE mo_advection_config,    ONLY: configure_advection
  USE mo_run_config,          ONLY: configure_run, output_mode
  USE mo_gribout_config,      ONLY: configure_gribout

  ! Control parameters: run control, dynamics, i/o
  !
  USE mo_run_config,          ONLY: &
    & test_mode,              &
    & dtime,                  & !    :
    & ltimer,                 & !    :
    & num_lev,                &
    & nshift,                 &
    & grid_generatingcenter,  & ! grid generating center
    & grid_generatingsubcenter  ! grid generating subcenter

  USE mo_ocean_nml_crosscheck,   ONLY: ocean_crosscheck
  USE mo_ocean_nml,              ONLY: i_sea_ice, no_tracer, use_omip_forcing, lhamocc, &
    & initialize_fromRestart, ncheckpoints

  USE mo_model_domain,        ONLY: t_patch_3d, p_patch_local_parent

  ! Horizontal grid
  !
  USE mo_grid_config,         ONLY: n_dom, use_dummy_cell_closure

  USE mo_build_decomposition, ONLY: build_decomposition
  USE mo_complete_subdivision,ONLY: setup_phys_patches
  USE mtime,                  ONLY: datetimeToString

  USE mo_ocean_ext_data,      ONLY: ext_data, construct_ocean_ext_data, destruct_ocean_ext_data
  USE mo_bgc_bcond,           ONLY: construct_bgc_ext_data, destruct_bgc_ext_data, ext_data_bgc
  USE mo_ocean_types,           ONLY: t_hydro_ocean_state, &
    & t_operator_coeff, t_solverCoeff_singlePrecision
  USE mo_ocean_state,           ONLY:  v_base, &
    & construct_hydro_ocean_base, &! destruct_hydro_ocean_base, &
    & construct_hydro_ocean_state, destruct_hydro_ocean_state, &
    & construct_patch_3d, destruct_patch_3d, ocean_default_list, ocean_restart_list, construct_ocean_nudge, &
    & construct_ocean_var_lists
  USE mo_ocean_initialization, ONLY: init_ho_base, &
    & init_ho_basins, init_coriolis_oce, init_patch_3d,   &
    & init_patch_3d
  USE mo_hamocc_output,        ONLY: construct_hamocc_var_lists, construct_hamocc_state, &
    &                                destruct_hamocc_state         
  USE mo_ocean_initial_conditions,  ONLY:  apply_initial_conditions, init_ocean_bathymetry
  USE mo_ocean_check_tools,     ONLY: init_oce_index
  USE mo_util_dbg_prnt,       ONLY: init_dbg_index
  USE mo_ext_data_types,      ONLY: t_external_data
  USE mo_ocean_physics_types,  ONLY: t_ho_params, construct_ho_params, v_params, &
                                   & destruct_ho_params
  USE mo_ocean_physics,          ONLY: init_ho_params                                 
  USE mo_operator_ocean_coeff_3d,ONLY: construct_operators_coefficients, &
    & destruct_operators_coefficients

  USE mo_hydro_ocean_run,     ONLY: perform_ho_stepping, &
    & prepare_ho_stepping, write_initial_ocean_timestep, &
    & end_ho_stepping
  USE mo_sea_ice_types,       ONLY: t_atmos_fluxes, t_sea_ice, v_sea_ice
  USE mo_ice_init_thermo,     ONLY: ice_init, construct_sea_ice, &
                                    &  construct_atmos_fluxes, destruct_sea_ice
  USE mo_ocean_surface_types, ONLY: t_ocean_surface, v_oce_sfc, t_atmos_for_ocean

  USE mo_ocean_forcing,       ONLY: construct_ocean_surface, destruct_ocean_forcing, &
                                    & construct_atmos_for_ocean, destruct_atmos_for_ocean
  USE mo_ocean_forcing,       ONLY: init_ocean_forcing
  USE mo_impl_constants,      ONLY: success

  USE mo_ocean_nudging,       ONLY: ocean_nudge

  USE mo_alloc_patches,        ONLY: destruct_patches, destruct_comm_patterns
  USE mo_ocean_read_namelists, ONLY: read_ocean_namelists
  USE mo_load_restart,         ONLY: read_restart_header, read_restart_files
  USE mo_restart_attributes,   ONLY: t_RestartAttributeList, getAttributesForRestarting
  USE mo_ocean_patch_setup,    ONLY: complete_ocean_patch
  USE mo_icon_comm_interface,  ONLY: construct_icon_communication, destruct_icon_communication
  USE mo_output_event_types,   ONLY: t_sim_step_info
  USE mo_grid_tools,           ONLY: create_dummy_cell_closure
  USE mo_ocean_diagnostics,    ONLY: construct_oce_diagnostics, destruct_oce_diagnostics
  USE mo_ocean_testbed,        ONLY: ocean_testbed
  USE mo_ocean_postprocessing, ONLY: ocean_postprocess
  USE mo_io_config,            ONLY: restartWritingParameters
  USE mo_bgc_icon_comm,        ONLY: hamocc_state
  USE mo_ocean_time_events,    ONLY: init_ocean_time_events, getCurrentDate_to_String
  !-------------------------------------------------------------
  ! For the coupling
  USE mo_ocean_coupling,      ONLY: construct_ocean_coupling, destruct_ocean_coupling
  USE mo_coupling_config,     ONLY: is_coupled_run

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: ocean_model
    PUBLIC :: construct_ocean_model, destruct_ocean_model
    PUBLIC :: ocean_patch_3d, ocean_state, operators_coefficients

    TYPE(t_patch_3d), POINTER                       :: ocean_patch_3d
    TYPE(t_atmos_for_ocean)                         :: p_as
    TYPE(t_atmos_fluxes)                            :: atmos_fluxes
    TYPE(t_operator_coeff)                          :: operators_coefficients
    TYPE(t_solverCoeff_singlePrecision)             :: solverCoefficients_sp
    TYPE(t_hydro_ocean_state), ALLOCATABLE, TARGET  :: ocean_state(:)
    
  !  TYPE(t_oce_timeseries), POINTER :: oce_ts

  CONTAINS


  !--------------------------------------------------------------------------
  !>
  SUBROUTINE ocean_model(oce_namelist_filename,shr_namelist_filename)

    CHARACTER(LEN=*), INTENT(in) :: oce_namelist_filename,shr_namelist_filename

    CHARACTER(*), PARAMETER :: method_name = "mo_ocean_model:ocean_model"

    !-------------------------------------------------------------------
    IF (isRestart()) THEN
      CALL read_restart_header("oce")
    END IF

    !-------------------------------------------------------------------
    ! initialize dynamic list of vertical axes
    !-------------------------------------------------------------------

    zaxisTypeList = t_zaxisTypeList()

    !-------------------------------------------------------------------
    CALL construct_ocean_model(oce_namelist_filename,shr_namelist_filename)

    !-------------------------------------------------------------------
    IF (isRestart() .OR. initialize_fromRestart) THEN
      ! This is an resumed integration. Read model state from restart file(s).
      CALL read_restart_files( ocean_patch_3d%p_patch_2d(1) )
      CALL message(TRIM(method_name),'normal exit from read_restart_files')
      !ELSE
      !  Prepare the initial conditions:
      !  forcing is part of the restart file
    END IF ! isRestart()
    !------------------------------------------------------------------
    ! Now start the time stepping:
    ! The special initial time step for the three time level schemes
    ! is executed within process_grid_level
    !------------------------------------------------------------------

    !------------------------------------------------------------------
    ! Initialize output file if necessary;
    ! Write out initial conditions.
    !------------------------------------------------------------------

    CALL prepare_output()

    CALL prepare_ho_stepping(ocean_patch_3d, operators_coefficients, &
      & ocean_state(1), v_oce_sfc, p_as, v_sea_ice, ext_data(1), isRestart(), solverCoefficients_sp)

    !------------------------------------------------------------------
    ! write initial state
    !------------------------------------------------------------------
    IF (output_mode%l_nml .and. .true.) THEN
      CALL write_initial_ocean_timestep(ocean_patch_3d,ocean_state(1),v_oce_sfc,v_sea_ice,hamocc_state, operators_coefficients)
    ENDIF
    !------------------------------------------------------------------
    SELECT CASE (test_mode)
      CASE (0)  !  ocean model
        CALL perform_ho_stepping( ocean_patch_3d, ocean_state, &
          & ext_data,                                          &
          & v_oce_sfc, v_params, p_as, atmos_fluxes,v_sea_ice, &
          & hamocc_state,                                      &
          & operators_coefficients,                            &
          & solverCoefficients_sp)

      CASE (1 : 1999) !
        CALL ocean_testbed( oce_namelist_filename,shr_namelist_filename,   &
          & ocean_patch_3d, ocean_state,                                   &
          & ext_data,                                                      &
          & v_oce_sfc, v_params, p_as, atmos_fluxes, v_sea_ice,            &
          & operators_coefficients,                                        &
          & solverCoefficients_sp)

      CASE (2000 : 3999) !
        CALL ocean_postprocess( oce_namelist_filename,shr_namelist_filename, &
          & ocean_patch_3d, ocean_state,                                     &
          & ext_data,                                                        &
          & operators_coefficients,                                          &
          & solverCoefficients_sp)

      CASE DEFAULT
        CALL finish(method_name, "Unknown test_mode")

    END SELECT

    !------------------------------------------------------------------

    CALL print_timer()

    !------------------------------------------------------------------
    !  cleaning up process
    CALL end_ho_stepping()
    CALL destruct_ocean_model()

  END SUBROUTINE ocean_model
  !--------------------------------------------------------------------------

  !--------------------------------------------------------------------------
  !>
  !!
!<Optimize:inUse>
  SUBROUTINE destruct_ocean_model()

    CHARACTER(*), PARAMETER :: method_name = "mo_ocean_model:destruct_ocean_model"

    INTEGER :: error_status

    !------------------------------------------------------------------
    !  cleaning up process
    !------------------------------------------------------------------
    CALL message(TRIM(method_name),'start to clean up')

    CALL destruct_oce_diagnostics()
    !------------------------------------------------------------------
    ! destruct ocean physics and forcing
    ! destruct ocean state is in control_model
    !------------------------------------------------------------------
!    CALL finalise_ho_integration(ocean_state, v_params, &
!      & p_as, atmos_fluxes, v_sea_ice, v_oce_sfc)
    CALL destruct_hydro_ocean_state(ocean_state)
    !CALL destruct_hydro_ocean_base(v_base)
    CALL destruct_ho_params(v_params)

    IF(no_tracer>0) CALL destruct_ocean_forcing(v_oce_sfc)
    CALL destruct_sea_ice(v_sea_ice)

    IF(lhamocc) CALL destruct_hamocc_state(hamocc_state)
    CALL destruct_atmos_for_ocean(p_as)
    !CALL destruct_atmos_fluxes(atmos_fluxes)

    !---------------------------------------------------------------------
    ! 13. Integration finished. Carry out the shared clean-up processes
    !---------------------------------------------------------------------
    ! Destruct external data state
    CALL destruct_ocean_ext_data
    IF(lhamocc)CALL destruct_bgc_ext_data

    ! deallocate ext_data array
    DEALLOCATE(ext_data, stat=error_status)
    IF (error_status/=success) THEN
      CALL finish(TRIM(method_name), 'deallocation of ext_data')
    ENDIF


    ! Destruct communication patterns
    CALL destruct_comm_patterns( ocean_patch_3d%p_patch_2d, p_patch_local_parent )

    !The 3D-ocean version of previous calls
    CALL destruct_patches( ocean_patch_3d%p_patch_2d )
    CALL destruct_patches( p_patch_local_parent )
    NULLIFY( ocean_patch_3d%p_patch_2d )
    CALL destruct_patch_3d( ocean_patch_3d )


    ! Delete variable lists

    IF (output_mode%l_nml) THEN
      CALL close_name_list_output
      CALL finish_mean_stream()
    ENDIF

    CALL destruct_icon_communication()
    CALL destruct_ocean_coupling ()

    CALL destruct_operators_coefficients(operators_coefficients, solverCoefficients_sp)

    ! close memory logging files
    CALL memory_log_terminate

    CALL message(TRIM(method_name),'clean-up finished')

  END SUBROUTINE destruct_ocean_model
  !--------------------------------------------------------------------------

  !--------------------------------------------------------------------------
  !>
  !!
  !! It does not include the restart processes, these are called from the calling method_name ocean_model
  !!
!<Optimize:inUse>
  SUBROUTINE construct_ocean_model(oce_namelist_filename,shr_namelist_filename)

    CHARACTER(LEN=*), INTENT(in) :: oce_namelist_filename,shr_namelist_filename

    CHARACTER(*), PARAMETER :: method_name = "mo_ocean_model:construct_ocean_model"
    INTEGER :: ist, error_status, dedicatedRestartProcs
    !-------------------------------------------------------------------

    !---------------------------------------------------------------------
    ! 1.1 Read namelists (newly) specified by the user; fill the
    !     corresponding sections of the configuration states.
    !---------------------------------------------------------------------

    CALL read_ocean_namelists(oce_namelist_filename,shr_namelist_filename)

    !---------------------------------------------------------------------
    ! 1.2 Cross-check namelist setups
    !---------------------------------------------------------------------
    CALL ocean_crosscheck()

    !---------------------------------------------------------------------
    ! 2. Call configure_run to finish filling the run_config state.
    !    This needs to be done very early (but anyway after atm_crosscheck)
    !    because some component of the state, e.g., num_lev, may be
    !    modified in this subroutine which affect the following CALLs.
    !---------------------------------------------------------------------
    CALL configure_run

    !-------------------------------------------------------------------
    CALL init_ocean_time_events()

    !-------------------------------------------------------------------
    ! 3.1 Initialize the mpi work groups
    !-------------------------------------------------------------------
    CALL restartWritingParameters(opt_dedicatedProcCount = dedicatedRestartProcs)
    CALL set_mpi_work_communicators(p_test_run, l_test_openmp, num_io_procs, &
      &                             dedicatedRestartProcs, num_test_pe, pio_type)

    !-------------------------------------------------------------------
    ! 3.2 Initialize various timers
    !-------------------------------------------------------------------
    CALL init_timer

    IF (ltimer) CALL timer_start(timer_model_init)

    !-------------------------------------------------------------------
    ! 4. Setup IO procs
    !-------------------------------------------------------------------
    ! If we belong to the I/O PEs just call xxx_io_main_proc before
    ! reading patches.  This routine will never return
    CALL init_io_processes()

    ! 4. Import patches
    !-------------------------------------------------------------------
    CALL build_decomposition(num_lev,nshift, is_ocean_decomposition =.TRUE., &
      & patch_3d=ocean_patch_3d)
    CALL construct_icon_communication(ocean_patch_3d%p_patch_2d(:), n_dom=1)
    CALL complete_ocean_patch(ocean_patch_3d%p_patch_2d(1))
    ! Setup the information for the physical patches
    CALL setup_phys_patches

    ! we need the nnow info
    CALL configure_dynamics ( n_dom )

    CALL construct_ocean_var_lists(ocean_patch_3d%p_patch_2d(1))
    
    IF(lhamocc) CALL construct_hamocc_var_lists(ocean_patch_3d%p_patch_2d(1))
    !------------------------------------------------------------------
    ! step 5b: allocate state variables
    !------------------------------------------------------------------
    ALLOCATE (ocean_state(n_dom), stat=ist)
    IF (ist /= success) THEN
      CALL finish(TRIM(method_name),'allocation for ocean_state failed')
    ENDIF
    !if(lhamocc)then
    !ALLOCATE (hamocc_state(n_dom), stat=ist)
    !IF (ist /= success) THEN
    !  CALL finish(TRIM(method_name),'allocation for hamocc_state failed')
    !ENDIF
    !ENDIF
    !---------------------------------------------------------------------
    ! 9. Horizontal and vertical grid(s) are now defined.
    !    Assign values to derived variables in the configuration states
    !---------------------------------------------------------------------

    CALL configure_gribout(grid_generatingcenter, grid_generatingsubcenter, n_dom)

    !    DO jg =1,n_dom
    !      !The 3D-ocean version of previous calls
    !      CALL configure_advection( jg, ocean_patch_3d%p_patch_2D(jg)%nlev, ocean_patch_3d%p_patch_2D(1)%nlev, &
    !        &                      iequations, iforcing, iqc, iqi, iqr, iqs, iqni, iqni_nuc, iqg, &
    !        &                      0, 1, .false., .true., ntracer )
    !    ENDDO

    !------------------------------------------------------------------
    ! 10. Create and optionally read external data fields
    !------------------------------------------------------------------
    ALLOCATE (ext_data(n_dom), stat=error_status)
    IF (error_status /= success) THEN
      CALL finish(TRIM(method_name),'allocation for ext_data failed')
    ENDIF

    ! allocate memory for oceanic external data and
    ! optionally read those data from netCDF file.
    CALL construct_ocean_ext_data(ocean_patch_3d%p_patch_2d(1:), ext_data)
    IF(lhamocc)CALL construct_bgc_ext_data(ocean_patch_3d%p_patch_2d(1:), ext_data,ext_data_bgc)
    ! initial analytic bathymetry via namelist
    CALL init_ocean_bathymetry(patch_3d=ocean_patch_3d,  &
      & cells_bathymetry=ext_data(1)%oce%bathymetry_c(:,:))

    !---------------------------------------------------------------------
    ! Prepare time integration
    CALL construct_ocean_states(ocean_patch_3d, ocean_state, ext_data, &
      & v_params, p_as, atmos_fluxes, v_sea_ice, v_oce_sfc, operators_coefficients, solverCoefficients_sp)!,p_int_state(1:))

    !---------------------------------------------------------------------
    IF (use_dummy_cell_closure) CALL create_dummy_cell_closure(ocean_patch_3D)
    

    ! initialize ocean indices for debug output (including 3-dim lsm)
    CALL init_oce_index(ocean_patch_3d%p_patch_2d,ocean_patch_3d, ocean_state, ext_data )

    CALL init_ho_params(ocean_patch_3d, v_params, p_as%fu10)

!    IF (.not. isRestart()) &
    CALL apply_initial_conditions(ocean_patch_3d, ocean_state(1), ext_data(1), operators_coefficients)
      
    ! initialize forcing after the initial conditions, since it may require knowledge
    ! of the initial conditions
    CALL init_ocean_forcing(ocean_patch_3d%p_patch_2d(1),  &
      &                     ocean_patch_3d,                &
      &                     ocean_state(1),         &
      &                     v_oce_sfc,            &
      &                     p_as%fu10)
      
    IF (i_sea_ice >= 1) &
      &   CALL ice_init(ocean_patch_3D, ocean_state(1), v_sea_ice, v_oce_sfc%cellThicknessUnderIce)

    IF (ltimer) CALL timer_stop(timer_model_init)

  END SUBROUTINE construct_ocean_model
  !--------------------------------------------------------------------------

  !--------------------------------------------------------------------------
  !>
  !! Simple method_name for preparing hydrostatic ocean model.
  !!
  !! @par Revision History
  !! Initial release by Stephan Lorenz, MPI-M (2010-07)
!<Optimize:inUse>
  SUBROUTINE construct_ocean_states(patch_3d, ocean_state, external_data, &
    & p_phys_param, p_as,&
    & atmos_fluxes, p_ice, p_oce_sfc, operators_coefficients, solverCoeff_sp)

    TYPE(t_patch_3d ),TARGET,   INTENT(inout)  :: patch_3d
    TYPE(t_hydro_ocean_state),  INTENT(inout)  :: ocean_state(n_dom)
    TYPE(t_external_data),      INTENT(inout)  :: external_data(n_dom)
    TYPE(t_ho_params),          INTENT(inout)  :: p_phys_param
    TYPE(t_atmos_for_ocean ),   INTENT(inout)  :: p_as
    TYPE(t_atmos_fluxes ),      INTENT(inout)  :: atmos_fluxes
    TYPE(t_sea_ice),            INTENT(inout)  :: p_ice
    TYPE(t_ocean_surface),      INTENT(inout)  :: p_oce_sfc
    TYPE(t_operator_coeff),     INTENT(inout), TARGET  :: operators_coefficients
    TYPE(t_solverCoeff_singlePrecision), INTENT(inout) :: solverCoeff_sp

    ! local variables
    CHARACTER(LEN=MAX_DATETIME_STR_LEN)         :: datestring
    INTEGER, PARAMETER :: kice = 1
    CHARACTER(LEN=*), PARAMETER :: &
      & method_name = 'mo_ocean_model:construct_ocean_states'

    CALL message (TRIM(method_name),'start')
    !------------------------------------------------------------------
    ! no grid refinement allowed here so far
    !------------------------------------------------------------------

    IF (n_dom > 1 ) THEN
      CALL finish(TRIM(method_name), ' N_DOM > 1 is not allowed')
    END IF

    !------------------------------------------------------------------
    ! construct ocean state and physics
    !------------------------------------------------------------------
    ! initialize ocean indices for debug output (before ocean state, no 3-dim)
    CALL init_dbg_index(patch_3d%p_patch_2d(1))!(patch_2D(1))

    ! hydro_ocean_base contains the 3-dimensional structures for the ocean state

    CALL construct_patch_3d(patch_3d)

    CALL construct_hydro_ocean_base(patch_3d%p_patch_2d(1), v_base)
    CALL init_ho_base (patch_3d%p_patch_2d(1), external_data(1), v_base)
    IF (use_omip_forcing .or. is_coupled_run()) CALL init_ho_basins(patch_3d%p_patch_2d(1), v_base)
    CALL init_coriolis_oce(patch_3d%p_patch_2d(1) )
    CALL init_patch_3d    (patch_3d,                external_data(1), v_base)
    !CALL init_patch_3D(patch_3D, v_base)

    CALL construct_operators_coefficients     ( patch_3d, operators_coefficients, solverCoeff_sp, ocean_default_list)
    !------------------------------------------------------------------
    ! construct ocean state and physics
    !------------------------------------------------------------------

    ! patch_2D and ocean_state have dimension n_dom
    CALL construct_hydro_ocean_state(patch_3d, ocean_state, ncheckpoints)
    ocean_state(1)%operator_coeff => operators_coefficients

    CALL construct_ocean_nudge(patch_3d%p_patch_2d(1),  ocean_nudge)

    if(lhamocc) CALL construct_hamocc_state(patch_3d%p_patch_2d, hamocc_state)


    CALL construct_ho_params(patch_3d%p_patch_2d(1), p_phys_param, ocean_restart_list)

    !------------------------------------------------------------------
    ! construct ocean initial conditions and forcing
    !------------------------------------------------------------------

    CALL construct_sea_ice(patch_3d, p_ice, kice)
    CALL construct_atmos_for_ocean(patch_3d%p_patch_2d(1), p_as)
    CALL construct_atmos_fluxes(patch_3d%p_patch_2d(1), atmos_fluxes, kice)

    CALL construct_ocean_surface(patch_3d, p_oce_sfc)
    CALL construct_ocean_coupling(ocean_patch_3d)

    !------------------------------------------------------------------
    datestring = getCurrentDate_to_String()
    CALL construct_oce_diagnostics( ocean_patch_3d, ocean_state(1), datestring)

    !------------------------------------------------------------------
    CALL message (TRIM(method_name),'end')

  END SUBROUTINE construct_ocean_states
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE init_io_processes()
    USE mo_time_config,         ONLY: time_config

    TYPE(t_sim_step_info)   :: sim_step_info
    INTEGER                 :: jstep0
    TYPE(t_RestartAttributeList), POINTER :: restartAttributes
    CHARACTER(LEN=*), PARAMETER :: &
      & method_name = 'mo_ocean_model:init_io_processes'
    
    IF (process_mpi_io_size < 1) THEN
      IF (output_mode%l_nml) THEN
        ! -----------------------------------------
        ! non-asynchronous I/O (performed by PE #0)
        ! -----------------------------------------
        CALL message(method_name,'synchronous namelist I/O scheme is enabled.')
      ENDIF
      ! nothing to do
      RETURN
    ENDIF

    ! Decide whether async vlist or name_list IO is to be used,
    ! only one of both may be enabled!

    IF (output_mode%l_nml) THEN
      ! -----------------------------------------
      ! asynchronous I/O
      ! -----------------------------------------
      !
      use_async_name_list_io = .TRUE.
      CALL message(method_name,'asynchronous namelist I/O scheme is enabled.')
      ! consistency check
      IF (my_process_is_io() .AND. (.NOT. my_process_is_mpi_test())) THEN

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
!         CALL name_list_io_main_proc(sim_step_info, isample=1)
        CALL name_list_io_main_proc(sim_step_info)
      END IF
    ELSE IF (my_process_is_io() .AND. (.NOT. my_process_is_mpi_test())) THEN
      ! Shut down MPI
      CALL stop_mpi
      STOP
    ENDIF
    
  END SUBROUTINE init_io_processes
  !-------------------------------------------------------------------

  !--------------------------------------------------------------------------
  SUBROUTINE prepare_output()
    USE mo_time_config,         ONLY: time_config

    CHARACTER(*), PARAMETER :: method_name = "mo_ocean_model:prepare_output"

    TYPE(t_sim_step_info)               :: sim_step_info
    INTEGER                             :: jstep0
    TYPE(t_RestartAttributeList), POINTER :: restartAttributes

    !------------------------------------------------------------------
    ! Initialize output file if necessary;
    ! Write out initial conditions.
    !------------------------------------------------------------------

    IF (output_mode%l_nml) THEN
!       WRITE(0,*)'process_mpi_io_size:',process_mpi_io_size
!       IF (process_mpi_io_size > 0) use_async_name_list_io = .TRUE.
      CALL parse_variable_groups()
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
      CALL init_mean_stream(ocean_patch_3d%p_patch_2d(1))
      CALL init_name_list_output(sim_step_info, opt_lprintlist=.TRUE.,opt_l_is_ocean=.TRUE.)
      CALL create_mipz_level_selections(output_file)
      CALL create_vertical_axes(output_file)
    ENDIF
  
  END SUBROUTINE prepare_output
  !--------------------------------------------------------------------------

END MODULE mo_ocean_model

