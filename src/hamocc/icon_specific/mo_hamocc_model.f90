!>
!! @page pagecontrolmodelf901 Main program for the ICON hamocc model
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
MODULE mo_hamocc_model

  USE mo_exception,           ONLY: message, message_text, finish
  USE mo_master_config,       ONLY: isRestart
  USE mo_parallel_config,     ONLY: p_test_run, l_test_openmp, num_io_procs, &
       &                            pio_type, num_test_pe
  USE mo_mpi,                 ONLY: set_mpi_work_communicators, process_mpi_io_size, &
       &                            stop_mpi, my_process_is_io, my_process_is_mpi_test,   &
       &                            process_mpi_io_size
  USE mo_timer,               ONLY: init_timer, timer_start, timer_stop, print_timer, &
       &  timer_model_init, timer_total
  USE mo_memory_log,              ONLY: memory_log_terminate
  USE mtime,                  ONLY: MAX_DATETIME_STR_LEN, datetimeToString
  USE mo_name_list_output_init, ONLY: init_name_list_output, parse_variable_groups, &
    &                                 create_vertical_axes, output_file
  USE mo_derived_variable_handling, ONLY: init_statistics_streams, finish_statistics_streams
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

  USE mo_ocean_ext_data,      ONLY: ext_data, construct_ocean_ext_data, destruct_ocean_ext_data
  USE mo_bgc_bcond,           ONLY: construct_bgc_ext_data, destruct_bgc_ext_data, ext_data_bgc
  USE mo_ocean_types,           ONLY: t_hydro_ocean_state, &
    & t_operator_coeff, t_solverCoeff_singlePrecision
  USE mo_ocean_state,           ONLY:  v_base, &
    & construct_hydro_ocean_base, &! destruct_hydro_ocean_base, &
    & construct_hydro_ocean_state, destruct_hydro_ocean_state, &
    & construct_patch_3d, destruct_patch_3d, ocean_default_list, ocean_restart_list, &
    & construct_ocean_var_lists, construct_ocean_nudge
    
  USE mo_ocean_initialization, ONLY: init_ho_base, &
    & init_ho_basins, init_coriolis_oce, init_patch_3d,   &
    & init_patch_3d
  USE mo_hamocc_output,        ONLY: construct_hamocc_var_lists, construct_hamocc_state, &
    &                                destruct_hamocc_state         
  USE mo_ocean_initial_conditions,  ONLY:  init_ocean_bathymetry
  USE mo_ocean_check_tools,     ONLY: init_oce_index
  USE mo_util_dbg_prnt,       ONLY: init_dbg_index
  USE mo_ext_data_types,      ONLY: t_external_data
  USE mo_operator_ocean_coeff_3d,ONLY: construct_operators_coefficients, &
    & destruct_operators_coefficients

   USE mo_impl_constants,      ONLY: success
  
  USE mo_alloc_patches,        ONLY: destruct_patches, destruct_comm_patterns
  USE mo_ocean_read_namelists, ONLY: read_ocean_namelists
  USE mo_load_restart,         ONLY: read_restart_header, read_restart_files
  USE mo_restart_attributes,   ONLY: t_RestartAttributeList, getAttributesForRestarting, ocean_initFromRestart_OVERRIDE
  USE mo_ocean_patch_setup,    ONLY: complete_ocean_patch
  USE mo_icon_comm_interface,  ONLY: construct_icon_communication, destruct_icon_communication
  USE mo_output_event_types,   ONLY: t_sim_step_info
  USE mo_grid_tools,           ONLY: create_dummy_cell_closure
  USE mo_io_config,            ONLY: restartWritingParameters
  USE mo_bgc_icon_comm,        ONLY: hamocc_state
  USE mo_ocean_time_events,    ONLY: init_ocean_time_events, getCurrentDate_to_String, &
    & ocean_time_nextStep, isCheckpoint, isEndOfThisRun, newNullDatetime
  USE mtime,                     ONLY: datetime, datetimeToString, deallocateDatetime

  !-------------------------------------------------------------
  USE mo_construct_icon_hamocc, ONLY: construct_icon_hamocc, destruct_icon_hamocc, init_icon_hamocc
  USE mo_ocean_hamocc_couple_state, ONLY:   t_ocean_to_hamocc_state, t_hamocc_to_ocean_state, &
    & t_ocean_transport_state, t_hamocc_ocean_state, &
    & hamocc_ocean_state, hamocc_ocean_state_list, &
    & construct_hamocc_ocean_state, destruct_hamocc_ocean_state
  USE mo_hamocc_ocean_physics,   ONLY: tracer_biochemistry_transport

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: hamocc_model
  PUBLIC :: patch_3d, ocean_operators_coefficients

  TYPE(t_patch_3d), POINTER                       :: patch_3d
  TYPE(t_operator_coeff), TArGET                  :: ocean_operators_coefficients ! needed for running the transport
  TYPE(t_solverCoeff_singlePrecision), TARGET     :: solverCoefficients_sp

  CONTAINS


  !--------------------------------------------------------------------------
  !>
  SUBROUTINE hamocc_model(hamocc_namelist_filename,shr_namelist_filename)

    CHARACTER(LEN=*), INTENT(in) :: hamocc_namelist_filename,shr_namelist_filename

    CHARACTER(*), PARAMETER :: method_name = "mo_hamocc_model:hamocc_model"

    !-------------------------------------------------------------------
    IF (isRestart()) THEN
      CALL read_restart_header("oce")
    END IF

    !-------------------------------------------------------------------
    ! initialize dynamic list of vertical axes
    !-------------------------------------------------------------------

    zaxisTypeList = t_zaxisTypeList()

    !-------------------------------------------------------------------
    CALL construct_hamocc_model(hamocc_namelist_filename,shr_namelist_filename)

    !-------------------------------------------------------------------
    IF (isRestart() .OR. initialize_fromRestart) THEN
      ocean_initFromRestart_OVERRIDE = initialize_fromRestart
      ! This is an resumed integration. Read model state from restart file(s).
      CALL read_restart_files( patch_3d%p_patch_2d(1) )
      CALL message(TRIM(method_name),'normal exit from read_restart_files')
      !ELSE
      !  Prepare the initial conditions:
      !  forcing is part of the restart file
    END IF ! isRestart()
    

    CALL init_icon_hamocc(hamocc_ocean_state)
    !------------------------------------------------------------------
    ! Initialize output file if necessary;
    ! Write out initial conditions.
    !------------------------------------------------------------------
    CALL prepare_output()

    !------------------------------------------------------------------
    ! Now start the time stepping:
    ! The special initial time step for the three time level schemes
    ! is executed within process_grid_level
    !------------------------------------------------------------------
    CALL hamocc_timestep()
    
    
    !------------------------------------------------------------------
    ! This is the end...
    !------------------------------------------------------------------   
    CALL print_timer()
    CALL destruct_hamocc_model()
     
  END SUBROUTINE hamocc_model
  !--------------------------------------------------------------------------
  
  !--------------------------------------------------------------------------
  SUBROUTINE hamocc_timestep()
 
    TYPE(datetime), POINTER        :: current_time     => NULL()
    CHARACTER(LEN=32)              :: datestring
    INTEGER :: jstep,  jstep0
    CHARACTER(*), PARAMETER :: method_name = "mo_hamocc_model:hamocc_timestep"
    
    current_time => newNullDatetime()
    jstep0 = 0
    jstep = jstep0
    
    CALL timer_start(timer_total)

    ! timestep ...
    !-------------------------------------------------------------------------
    TIME_LOOP: DO
      jstep = jstep + 1    
      current_time = ocean_time_nextStep()

      CALL datetimeToString(current_time, datestring)
      WRITE(message_text,'(a,i10,2a)') '  Begin of timestep =',jstep,'  datetime:  ', datestring
      CALL message (TRIM(method_name), message_text)

      CALL tracer_biochemistry_transport(hamocc_ocean_state, ocean_operators_coefficients, current_time)
    !-------------------------------------------------------------------------
      IF (isEndOfThisRun()) THEN
        ! leave time loop
        EXIT TIME_LOOP
      END IF

    ENDDO TIME_LOOP
    !-------------------------------------------------------------------------
     
    CALL timer_stop(timer_total)
    
    CALL deallocateDatetime(current_time)   


  END SUBROUTINE hamocc_timestep
  !--------------------------------------------------------------------------

  !--------------------------------------------------------------------------
  !>
  !!
  SUBROUTINE destruct_hamocc_model()

    CHARACTER(*), PARAMETER :: method_name = "destruct_hamocc_model"

    INTEGER :: error_status

    !------------------------------------------------------------------
    !  cleaning up process
    !------------------------------------------------------------------
    CALL message(TRIM(method_name),'start to clean up')
    
    CALL destruct_hamocc_ocean_state()
    CALL destruct_icon_hamocc()

    ! Destruct communication patterns
    CALL destruct_comm_patterns( patch_3d%p_patch_2d, p_patch_local_parent )

    !The 3D-ocean version of previous calls
    CALL destruct_patches( patch_3d%p_patch_2d )
    CALL destruct_patches( p_patch_local_parent )
    NULLIFY( patch_3d%p_patch_2d )
    CALL destruct_patch_3d( patch_3d )


    ! Delete variable lists

    IF (output_mode%l_nml) THEN
      CALL close_name_list_output
      CALL finish_statistics_streams
    ENDIF

    CALL destruct_icon_communication()

    CALL destruct_operators_coefficients(ocean_operators_coefficients, solverCoefficients_sp)

    ! close memory logging files
    CALL memory_log_terminate

    CALL message(TRIM(method_name),'clean-up finished')

  END SUBROUTINE destruct_hamocc_model
  !--------------------------------------------------------------------------
  
  !--------------------------------------------------------------------------
  !>
  !!
  !! It does not include the restart processes, these are called from the calling method_name ocean_model
  !!
!<Optimize:inUse>
  SUBROUTINE construct_hamocc_model(hamocc_namelist_filename,shr_namelist_filename)

    CHARACTER(LEN=*), INTENT(in) :: hamocc_namelist_filename,shr_namelist_filename

    CHARACTER(*), PARAMETER :: method_name = "mo_ocean_model:construct_ocean_model"
    INTEGER :: ist, error_status, dedicatedRestartProcs
    !-------------------------------------------------------------------

    !---------------------------------------------------------------------
    ! 1.1 Read namelists (newly) specified by the user; fill the
    !     corresponding sections of the configuration states.
    !---------------------------------------------------------------------
    ! for the moment we just read the hamocc nemalis
    CALL read_ocean_namelists(hamocc_namelist_filename,shr_namelist_filename)

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
      & patch_3d=patch_3d)
    CALL construct_icon_communication(patch_3d%p_patch_2d(:), n_dom=1)
    CALL complete_ocean_patch(patch_3d%p_patch_2d(1))
    ! Setup the information for the physical patches
    CALL setup_phys_patches

    ! we need the nnow info
    CALL configure_dynamics ( n_dom )

    !------------------------------------------------------------------
    CALL construct_ocean_var_lists(patch_3d%p_patch_2d(1))
    !------------------------------------------------------------------
    ! step 5b: allocate state variables
    !------------------------------------------------------------------
!     ALLOCATE (ocean_state(n_dom), stat=ist)
!     IF (ist /= success) THEN
!       CALL finish(TRIM(method_name),'allocation for ocean_state failed')
!     ENDIF

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

    !------------------------------------------------------------------
    ! 10. Create and optionally read external data fields
    !------------------------------------------------------------------
     ALLOCATE (ext_data(1), stat=error_status)
     IF (error_status /= success) THEN
       CALL finish(TRIM(method_name),'allocation for ext_data failed')
    ENDIF
    
! 
!     ! allocate memory for oceanic external data and
!     ! optionally read those data from netCDF file.
    CALL construct_ocean_ext_data(patch_3d%p_patch_2d(1:), ext_data)
    ! initial analytic bathymetry via namelist
    CALL init_ocean_bathymetry(patch_3d=patch_3d,  &
      & cells_bathymetry=ext_data(1)%oce%bathymetry_c(:,:))    


    CALL construct_patch_3d(patch_3d)

    CALL construct_hydro_ocean_base(patch_3d%p_patch_2d(1), v_base)
    CALL init_ho_base (patch_3d%p_patch_2d(1), ext_data(1), v_base)
    CALL init_ho_basins(patch_3d%p_patch_2d(1), v_base) ! This initializes the wet_c,..., for all cells ! unbelievable !
    CALL init_coriolis_oce(patch_3d%p_patch_2d(1) )
    CALL init_patch_3d    (patch_3d, ext_data(1), v_base)
    CALL construct_operators_coefficients     ( patch_3d, ocean_operators_coefficients, solverCoefficients_sp)

    IF (use_dummy_cell_closure) CALL create_dummy_cell_closure(patch_3d)
    !------------------------------------------------------------------
    
    CALL construct_icon_hamocc(patch_3d, ext_data(1))
    CALL construct_hamocc_ocean_state(patch_3d)   
    
    !---------------------------------------------------------------------

    IF (ltimer) CALL timer_stop(timer_model_init)

  END SUBROUTINE construct_hamocc_model
  !--------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE init_io_processes()
    USE mo_time_config,         ONLY: time_config

    TYPE(t_sim_step_info)   :: sim_step_info
    INTEGER                 :: jstep0
    TYPE(t_RestartAttributeList), POINTER :: restartAttributes
    CHARACTER(LEN=*), PARAMETER :: &
      & method_name = 'mo_hamocc_model:init_io_processes'
    
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

    CHARACTER(*), PARAMETER :: method_name = "mo_hamocc_model:prepare_output"

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
      CALL init_statistics_streams
      CALL init_name_list_output(sim_step_info, opt_lprintlist=.TRUE.,opt_l_is_ocean=.TRUE.)
      CALL create_mipz_level_selections(output_file)
      CALL create_vertical_axes(output_file)
    ENDIF
  
  END SUBROUTINE prepare_output
  !--------------------------------------------------------------------------

END MODULE mo_hamocc_model

