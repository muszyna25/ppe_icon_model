!>
!! @page pagecontrolmodelf901 Main program for the ICON ocean model
!!
!! @par Revision History
!!
!! @par Copyright
!! 2002-2011 by DWD and MPI-M
!! This software is provided for non-commercial use only.
!! See the LICENSE and the WARRANTY conditions.
!!
!! @par License
!! The use of ICON is hereby granted free of charge for an unlimited time,
!! provided the following rules are accepted and applied:
!! <ol>
!! <li> You may use or modify this code for your own non commercial and non
!!    violent purposes.
!! <li> The code may not be re-distributed without the consent of the authors.
!! <li> The copyright notice and statement of authorship must appear in all
!!    copies.
!! <li> You accept the warranty conditions (see WARRANTY).
!! <li> In case you intend to use the code commercially, we oblige you to sign
!!    an according license agreement with DWD and MPI-M.
!! </ol>
!!
!! @par Warranty
!! This code has been tested up to a certain level. Defects and weaknesses,
!! which may be included in the code, do not establish any warranties by the
!! authors.
!! The authors do not make any warranty, express or implied, or assume any
!! liability or responsibility for the use, acquisition or application of this
!! software.
!!
MODULE mo_ocean_model

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message, finish
  USE mo_master_control,      ONLY: is_restart_run, get_my_process_name, get_my_model_no
  USE mo_parallel_config,     ONLY: p_test_run, l_test_openmp, num_io_procs , num_restart_procs
  USE mo_mpi,                 ONLY: my_process_is_io,set_mpi_work_communicators,p_pe_work, process_mpi_io_size
  USE mo_timer,               ONLY: init_timer, timer_start, timer_stop, print_timer, timer_model_init
  USE mo_datetime,            ONLY: t_datetime
  USE mo_name_list_output_init, ONLY: init_name_list_output, parse_variable_groups
  USE mo_name_list_output,    ONLY: close_name_list_output
  USE mo_name_list_output_config,  ONLY: use_async_name_list_io
  USE mo_dynamics_config,     ONLY: iequations, configure_dynamics

  !  USE mo_advection_config,    ONLY: configure_advection
  USE mo_run_config,          ONLY: configure_run, output_mode
  USE mo_gribout_config,      ONLY: configure_gribout

  ! Control parameters: run control, dynamics, i/o
  !
  USE mo_run_config,          ONLY: &
    & dtime,                  & !    :
    & nsteps,                 & !    :
    & ltimer,                 & !    :
    & iforcing,               & !
    & num_lev, num_levp1,     &
    & nlev, nlevp1,           &
    & iqc, iqi, iqr, iqs,     &
    & iqni, iqni_nuc, iqg,    &
    & nshift, ntracer,        &
    & grid_generatingcenter,  & ! grid generating center
    & grid_generatingsubcenter  ! grid generating subcenter

  USE mo_ocean_nml_crosscheck,   ONLY: oce_crosscheck
  USE mo_ocean_nml,              ONLY: i_sea_ice, use_new_forcing

  USE mo_model_domain,        ONLY: t_patch, t_patch_3d, p_patch_local_parent

  ! Horizontal grid
  !
  USE mo_grid_config,         ONLY: n_dom, use_dummy_cell_closure

  USE mo_build_decomposition, ONLY: build_decomposition
  USE mo_complete_subdivision,ONLY: setup_phys_patches

  USE mo_ocean_ext_data,      ONLY: ext_data, construct_ocean_ext_data, destruct_ocean_ext_data
  USE mo_oce_types,           ONLY: t_hydro_ocean_state, &
    & t_hydro_ocean_acc, t_hydro_ocean_diag, &
    & t_hydro_ocean_prog
  USE mo_oce_state,           ONLY:  v_base, &
    & construct_hydro_ocean_base, &! destruct_hydro_ocean_base, &
    & construct_hydro_ocean_state, destruct_hydro_ocean_state, &
    & construct_patch_3d, destruct_patch_3d, ocean_default_list
  USE mo_ocean_initialization,    ONLY:    setup_ocean_namelists,  init_ho_base, &
    & init_ho_basins, init_coriolis_oce, init_oce_config,  init_patch_3d,   &
    & init_patch_3d, setup_ocean_namelists
  USE mo_ocean_initial_conditions,  ONLY:  apply_initial_conditions
  USE mo_oce_check_tools,     ONLY: init_oce_index
  USE mo_util_dbg_prnt,       ONLY: init_dbg_index
  USE mo_ext_data_types,      ONLY: t_external_data
  USE mo_oce_physics,         ONLY: t_ho_params, construct_ho_params, init_ho_params, v_params
  USE mo_operator_ocean_coeff_3d,ONLY: t_operator_coeff, construct_operators_coefficients, &
    & destruct_operators_coefficients

  USE mo_hydro_ocean_run,     ONLY: perform_ho_stepping,&
    & prepare_ho_stepping,  finalise_ho_integration
  USE mo_sea_ice_types,       ONLY: t_atmos_fluxes, t_atmos_for_ocean, &
    & v_sfc_flx, v_sea_ice, t_sfc_flx, t_sea_ice
  USE mo_sea_ice,             ONLY: ice_init, &
    & construct_atmos_for_ocean, construct_atmos_fluxes, construct_sea_ice
  USE mo_oce_forcing,         ONLY: construct_ocean_forcing, init_ocean_forcing
  USE mo_impl_constants,      ONLY: max_char_length, success

  USE mo_alloc_patches,       ONLY: destruct_patches
  USE mo_ocean_read_namelists, ONLY: read_ocean_namelists
  USE mo_io_restart,          ONLY: read_restart_info_file, read_restart_files
  USE mo_io_restart_namelist, ONLY: read_restart_namelists
  USE mo_io_restart_attributes,ONLY: read_restart_attributes, get_restart_attribute
  USE mo_oce_patch_setup,     ONLY: complete_ocean_patch
  USE mo_time_config,         ONLY: time_config
  USE mo_icon_comm_interface, ONLY: construct_icon_communication, destruct_icon_communication
  USE mo_mtime_extensions,    ONLY: get_datetime_string
  USE mo_output_event_types,  ONLY: t_sim_step_info
  USE mtime,                  ONLY: setcalendar, proleptic_gregorian
  USE mo_grid_tools,          ONLY: create_dummy_cell_closure

  !-------------------------------------------------------------
  ! For the coupling
  USE mo_ocean_coupling,      ONLY: construct_ocean_coupling, destruct_ocean_coupling

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: ocean_model
  PUBLIC :: construct_ocean_model, destruct_ocean_model
  PUBLIC :: ocean_patch_3d, ocean_state, operators_coefficients

  TYPE(t_patch_3d), POINTER :: ocean_patch_3d
  TYPE(t_atmos_for_ocean)                         :: p_as
  TYPE(t_atmos_fluxes)                            :: p_atm_f
  TYPE(t_operator_coeff)                          :: operators_coefficients
  TYPE(t_hydro_ocean_state), ALLOCATABLE, TARGET :: ocean_state(:)
  TYPE(t_datetime)                                :: start_datetime


CONTAINS


  !--------------------------------------------------------------------------
  !>
  !!
  SUBROUTINE ocean_model(oce_namelist_filename,shr_namelist_filename)

    CHARACTER(LEN=*), INTENT(in) :: oce_namelist_filename,shr_namelist_filename

    CHARACTER(*), PARAMETER :: routine = "mo_ocean_model:ocean_model"

    INTEGER :: jg, ist
    INTEGER :: error_status
    TYPE(t_sim_step_info) :: sim_step_info
    INTEGER :: jstep0

    !-------------------------------------------------------------------
    IF (is_restart_run()) THEN
      CALL read_restart_ocean_namelists('restart_oce.nc')
    END IF

    !-------------------------------------------------------------------
    CALL construct_ocean_model(oce_namelist_filename,shr_namelist_filename)

    !-------------------------------------------------------------------
    IF (is_restart_run()) THEN
      ! This is an resumed integration. Read model state from restart file(s).
#ifdef NOMPI
      CALL read_restart_files
#else
      jg = 1 !no nesting
      !DO jg = ,n_dom
      CALL read_restart_files( ocean_patch_3d%p_patch_2d(jg) )
      !END DO
#endif
      CALL message(TRIM(routine),'normal exit from read_restart_files')
      !ELSE
      !  Prepare the initial conditions:
      !  forcing is part of the restart file
    END IF ! is_restart_run()
    !------------------------------------------------------------------
    ! Now start the time stepping:
    ! The special initial time step for the three time level schemes
    ! is executed within process_grid_level
    !------------------------------------------------------------------

    !------------------------------------------------------------------
    ! Initialize output file if necessary;
    ! Write out initial conditions.
    !------------------------------------------------------------------

    IF (output_mode%l_nml) THEN
      WRITE(0,*)'process_mpi_io_size:',process_mpi_io_size
      IF (process_mpi_io_size > 0) use_async_name_list_io = .TRUE.
      CALL parse_variable_groups()
      CALL setcalendar(proleptic_gregorian)
      ! compute sim_start, sim_end
      CALL get_datetime_string(sim_step_info%sim_start, time_config%ini_datetime)
      CALL get_datetime_string(sim_step_info%sim_end,   time_config%end_datetime)
      CALL get_datetime_string(sim_step_info%restart_time,  time_config%cur_datetime, &
        & INT(time_config%dt_restart))
      CALL get_datetime_string(sim_step_info%run_start, time_config%cur_datetime)
      sim_step_info%dtime      = dtime
      sim_step_info%iadv_rcf   = 1
      jstep0 = 0
      IF (is_restart_run() .AND. .NOT. time_config%is_relative_time) THEN
        ! get start counter for time loop from restart file:
        CALL get_restart_attribute("jstep", jstep0)
      END IF
      sim_step_info%jstep0    = jstep0
      CALL init_name_list_output(sim_step_info, opt_lprintlist=.TRUE.,opt_l_is_ocean=.TRUE.)
    ENDIF

    CALL prepare_ho_stepping(ocean_patch_3d,operators_coefficients,ocean_state(1),v_params, is_restart_run())

    !------------------------------------------------------------------
    CALL perform_ho_stepping( ocean_patch_3d, ocean_state,                    &
      & ext_data, start_datetime,                     &
      & (nsteps == INT(time_config%dt_restart/dtime)),&
      & v_sfc_flx,                                    &
      & v_params, p_as, p_atm_f,v_sea_ice,operators_coefficients)
    !------------------------------------------------------------------

    CALL print_timer()

    !------------------------------------------------------------------
    !  cleaning up process
    CALL destruct_ocean_model()

  END SUBROUTINE ocean_model
  !--------------------------------------------------------------------------

  !--------------------------------------------------------------------------
  !>
  !!
  SUBROUTINE destruct_ocean_model()

    CHARACTER(*), PARAMETER :: routine = "mo_ocean_model:destruct_ocean_model"

    INTEGER :: error_status

    !------------------------------------------------------------------
    !  cleaning up process
    !------------------------------------------------------------------
    CALL message(TRIM(routine),'start to clean up')

    CALL finalise_ho_integration(ocean_state, v_params, &
      & p_as, p_atm_f, v_sea_ice, v_sfc_flx)

    !---------------------------------------------------------------------
    ! 13. Integration finished. Carry out the shared clean-up processes
    !---------------------------------------------------------------------
    ! Destruct external data state
    CALL destruct_ocean_ext_data

    ! deallocate ext_data array
    DEALLOCATE(ext_data, stat=error_status)
    IF (error_status/=success) THEN
      CALL finish(TRIM(routine), 'deallocation of ext_data')
    ENDIF

    !The 3D-ocean version of previous calls
    CALL destruct_patches( ocean_patch_3d%p_patch_2d )
    CALL destruct_patches( p_patch_local_parent )
    NULLIFY( ocean_patch_3d%p_patch_2d )
    CALL destruct_patch_3d( ocean_patch_3d )


    ! Delete variable lists

    IF (output_mode%l_nml) THEN
      CALL close_name_list_output
    ENDIF

    CALL destruct_icon_communication()

    CALL destruct_ocean_coupling ()

    CALL destruct_operators_coefficients(operators_coefficients)

    CALL message(TRIM(routine),'clean-up finished')

  END SUBROUTINE destruct_ocean_model
  !--------------------------------------------------------------------------

  !--------------------------------------------------------------------------
  !>
  !!
  !! It does not include the restart processes, these are called from the calling routine ocean_model
  !!
  SUBROUTINE construct_ocean_model(oce_namelist_filename,shr_namelist_filename)

    CHARACTER(LEN=*), INTENT(in) :: oce_namelist_filename,shr_namelist_filename

    CHARACTER(*), PARAMETER :: routine = "mo_ocean_model:construct_ocean_model"

    INTEGER :: ist

    INTEGER :: error_status
    !-------------------------------------------------------------------

    !---------------------------------------------------------------------
    ! 1.1 Read namelists (newly) specified by the user; fill the
    !     corresponding sections of the configuration states.
    !---------------------------------------------------------------------

    CALL read_ocean_namelists(oce_namelist_filename,shr_namelist_filename)

    !---------------------------------------------------------------------
    ! 1.2 Cross-check namelist setups
    !---------------------------------------------------------------------
    CALL oce_crosscheck()

    !---------------------------------------------------------------------
    ! 2. Call configure_run to finish filling the run_config state.
    !    This needs to be done very early (but anyway after atm_crosscheck)
    !    because some component of the state, e.g., num_lev, may be
    !    modified in this subroutine which affect the following CALLs.
    !---------------------------------------------------------------------
    CALL configure_run

    !-------------------------------------------------------------------
    ! 3.1 Initialize the mpi work groups
    !-------------------------------------------------------------------
    CALL set_mpi_work_communicators(p_test_run, l_test_openmp, num_io_procs, num_restart_procs)

    !-------------------------------------------------------------------
    ! 3.2 Initialize various timers
    !-------------------------------------------------------------------
    CALL init_timer

    IF (ltimer) CALL timer_start(timer_model_init)

    !-------------------------------------------------------------------
    ! Initialize date and time
    !-------------------------------------------------------------------
    start_datetime = time_config%cur_datetime

    !-------------------------------------------------------------------
    ! 4. Import patches
    !-------------------------------------------------------------------
    CALL build_decomposition(num_lev,num_levp1,nshift, is_ocean_decomposition =.TRUE., &
      & patch_3d=ocean_patch_3d)
    CALL construct_icon_communication(ocean_patch_3d%p_patch_2d(:), n_dom=1)
    CALL complete_ocean_patch(ocean_patch_3d%p_patch_2d(1))

    !--------------------------------------------
    ! Setup the information for the physical patches
    CALL setup_phys_patches

    CALL setup_ocean_namelists(ocean_patch_3d%p_patch_2d(1))
    !------------------------------------------------------------------
    ! step 5b: allocate state variables
    !------------------------------------------------------------------
    ALLOCATE (ocean_state(n_dom), stat=ist)
    IF (ist /= success) THEN
      CALL finish(TRIM(routine),'allocation for ocean_state failed')
    ENDIF

    !---------------------------------------------------------------------
    ! 9. Horizontal and vertical grid(s) are now defined.
    !    Assign values to derived variables in the configuration states
    !---------------------------------------------------------------------

    CALL configure_dynamics ( n_dom )

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
      CALL finish(TRIM(routine),'allocation for ext_data failed')
    ENDIF

    ! allocate memory for oceanic external data and
    ! optionally read those data from netCDF file.
    CALL construct_ocean_ext_data(ocean_patch_3d%p_patch_2d(1:), ext_data)

    ! Prepare time integration
    CALL construct_ocean_states(ocean_patch_3d, ocean_state, ext_data, v_sfc_flx, &
      & v_params, p_as, p_atm_f, v_sea_ice,operators_coefficients)!,p_int_state(1:))


  END SUBROUTINE construct_ocean_model
  !--------------------------------------------------------------------------

  !--------------------------------------------------------------------------
  !>
  !! Simple routine for preparing hydrostatic ocean model.
  !!
  !! @par Revision History
  !! Initial release by Stephan Lorenz, MPI-M (2010-07)
  SUBROUTINE construct_ocean_states(patch_3d, p_os, p_ext_data, p_sfc_flx, &
    & p_phys_param, p_as,&
    & p_atm_f, p_ice, p_op_coeff)

    TYPE(t_patch_3d ),TARGET,   INTENT(inout)  :: patch_3d
    TYPE(t_hydro_ocean_state),  INTENT(inout)  :: p_os(n_dom)
    TYPE(t_external_data),      INTENT(inout)  :: p_ext_data(n_dom)
    TYPE(t_sfc_flx),            INTENT(inout)  :: p_sfc_flx
    TYPE(t_ho_params),          INTENT(inout)  :: p_phys_param
    TYPE(t_atmos_for_ocean ),   INTENT(inout)  :: p_as
    TYPE(t_atmos_fluxes ),      INTENT(inout)  :: p_atm_f
    TYPE(t_sea_ice),            INTENT(inout)  :: p_ice
    TYPE(t_operator_coeff),     INTENT(inout)  :: p_op_coeff

    ! local variables
    INTEGER, PARAMETER :: kice = 1
    INTEGER :: jg
    CHARACTER(LEN=max_char_length), PARAMETER :: &
      & routine = 'mo_test_hydro_ocean:construct_ocean_states'

    CALL message (TRIM(routine),'start')
    !------------------------------------------------------------------
    ! no grid refinement allowed here so far
    !------------------------------------------------------------------

    IF (n_dom > 1 ) THEN
      CALL finish(TRIM(routine), ' N_DOM > 1 is not allowed')
    END IF
    jg = n_dom

    !------------------------------------------------------------------
    ! construct ocean state and physics
    !------------------------------------------------------------------
    CALL init_oce_config

    ! initialize ocean indices for debug output (before ocean state, no 3-dim)
    CALL init_dbg_index(patch_3d%p_patch_2d(jg))!(patch_2D(jg))

    ! hydro_ocean_base contains the 3-dimensional structures for the ocean state

    CALL construct_patch_3d(patch_3d)

    CALL construct_hydro_ocean_base(patch_3d%p_patch_2d(jg), v_base)
    CALL init_ho_base     (patch_3d%p_patch_2d(jg), p_ext_data(jg), v_base)
    CALL init_ho_basins   (patch_3d%p_patch_2d(jg),                 v_base)
    CALL init_coriolis_oce(patch_3d%p_patch_2d(jg) )
    CALL init_patch_3d    (patch_3d,                p_ext_data(jg), v_base)
    !CALL init_patch_3D(patch_3D, v_base)

    CALL construct_operators_coefficients     ( patch_3d, p_op_coeff, ocean_default_list)
    !------------------------------------------------------------------
    ! construct ocean state and physics
    !------------------------------------------------------------------

    ! patch_2D and p_os have dimension n_dom
    CALL construct_hydro_ocean_state(patch_3d%p_patch_2d, p_os)

    CALL construct_ho_params(patch_3d%p_patch_2d(jg), p_phys_param)

    !------------------------------------------------------------------
    ! construct ocean initial conditions and forcing
    !------------------------------------------------------------------

    CALL construct_sea_ice(patch_3d, p_ice, kice)
    CALL construct_atmos_for_ocean(patch_3d%p_patch_2d(jg), p_as)
    CALL construct_atmos_fluxes(patch_3d%p_patch_2d(jg), p_atm_f, kice)

    CALL construct_ocean_forcing(patch_3d%p_patch_2d(jg),p_sfc_flx, ocean_default_list)
    CALL construct_ocean_coupling(ocean_patch_3d)

    !------------------------------------------------------------------
    CALL init_ho_params(patch_3d, p_phys_param)

    CALL apply_initial_conditions(patch_3d%p_patch_2d(jg),patch_3d, p_os(jg), p_ext_data(jg), p_op_coeff)
    ! initialize forcing after the initial conditions, since it may require knowledge
    ! of the initial conditions
    CALL init_ocean_forcing(patch_3d%p_patch_2d(jg)%cells%All, &
      & patch_3d%lsm_c(:,1,:), &
      & p_sfc_flx)

    IF (i_sea_ice >= 1) &
      &   CALL ice_init(patch_3D, p_os(jg), p_ice)

    ! initialize ocean indices for debug output (including 3-dim lsm)
    CALL init_oce_index( patch_3d%p_patch_2d,patch_3d, p_os, p_ext_data )

    IF (use_dummy_cell_closure) CALL create_dummy_cell_closure(patch_3D)

    CALL message (TRIM(routine),'end')

  END SUBROUTINE construct_ocean_states
  !-------------------------------------------------------------------------


  !--------------------------------------------------------------------------
  !>
  !!
  SUBROUTINE read_restart_ocean_namelists(restart_filename)

    CHARACTER(LEN=*), INTENT(in) :: restart_filename

    CHARACTER(*), PARAMETER :: routine = "mo_ocean_model:read_restart_namelists"

    CHARACTER(LEN=max_char_length) :: grid_file_name
    LOGICAL :: lsuccess
    !-------------------------------------------------------------------

    ! First read the restart master file (ASCII format) to find out
    ! which NetCDF files the model should read.
    ! Comment by Hui Wan:
    ! The namelist variable atmo_restart_info_filename should be
    ! an input argument of the subroutine read_restart_info_file.
    CALL read_restart_info_file(grid_file_name, lsuccess) ! out, out

    IF (lsuccess) THEN
      CALL message( TRIM(routine),                          &
        & 'Running model in restart mode. '       &
        & //'Horizontal grid should be read from '&
        & //TRIM(grid_file_name) )
    ELSE
      CALL finish(TRIM(routine),'Failed to read restart info file')
    END IF

    ! Read all namelists used in the previous run
    ! and store them in a buffer. These values will overwrite the
    ! model default, and will later be overwritten if the user has
    ! specified something different for this integraion.
    CALL read_restart_namelists(restart_filename)
    CALL message(TRIM(routine), 'read namelists from ocean restart file')

    ! Read all global attributs in the restart file and store them in a buffer.

    CALL read_restart_attributes(restart_filename)
    CALL message(TRIM(routine), 'read global attributes from ocean restart file')

  END SUBROUTINE read_restart_ocean_namelists
  !--------------------------------------------------------------------------


END MODULE mo_ocean_model

