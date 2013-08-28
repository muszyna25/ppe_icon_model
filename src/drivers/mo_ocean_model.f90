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
  USE mo_output,              ONLY: init_output_files, write_output_oce, close_output_files
  USE mo_name_list_output_init, ONLY: init_name_list_output, parse_variable_groups
  USE mo_name_list_output,    ONLY: close_name_list_output
  USE mo_name_list_output_config,  ONLY: use_async_name_list_io
  USE mo_grid_config,         ONLY: n_dom
  USE mo_dynamics_config,     ONLY: iequations
  
  USE mo_io_async,            ONLY: vlist_io_main_proc
  
  !  USE mo_advection_config,    ONLY: configure_advection
  USE mo_dynamics_config,     ONLY: configure_dynamics  ! subroutine
  USE mo_run_config,          ONLY: configure_run, output_mode
  USE mo_gribout_config,      ONLY: configure_gribout
  
  ! Control parameters: run control, dynamics, i/o
  !
  USE mo_io_config,           ONLY:  lwrite_initial,n_ios
  USE mo_run_config,          ONLY: &
    & dtime,                  & !    :
    & nsteps,                 & !    :
    & ltimer,                 & !    :
  !& ldump_states,           & ! flag if states should be dumped
    & lrestore_states,        & ! flag if states should be restored
    & iforcing,               & !
    & num_lev, num_levp1,     &
    & nlev, nlevp1,           &
    & iqc, iqi, iqr, iqs,     &
    & iqni, iqni_nuc, iqg,    &
    & nshift, ntracer,        &
    & grid_generatingcenter,  & ! grid generating center
    & grid_generatingsubcenter  ! grid generating subcenter
  
  USE mo_ocean_nml_crosscheck,   ONLY: oce_crosscheck
  
  USE mo_model_domain,        ONLY: t_patch, t_patch_3d, p_patch_local_parent
  
  ! Horizontal grid
  !
  USE mo_grid_config,         ONLY: n_dom
  
  USE mo_oce_state,           ONLY: t_hydro_ocean_state, setup_ocean_namelists
  USE mo_build_decomposition, ONLY: build_decomposition
  USE mo_complete_subdivision,ONLY: setup_phys_patches
  
  USE mo_impl_constants,      ONLY: success !, ihs_ocean
  
  ! External data
  USE mo_ocean_ext_data,      ONLY: ext_data, construct_ocean_ext_data, destruct_ocean_ext_data
  
  USE mo_hydro_ocean_run,     ONLY: perform_ho_stepping,&
    & prepare_ho_integration,&
    & finalise_ho_integration
  USE mo_operator_ocean_coeff_3d, ONLY: t_operator_coeff
  USE mo_oce_physics,         ONLY: v_params!, t_ho_params, t_ho_physics
  USE mo_sea_ice_types,       ONLY: t_atmos_fluxes, t_atmos_for_ocean, &
    & v_sfc_flx, v_sea_ice!, t_sfc_flx
  USE mo_impl_constants,      ONLY: max_char_length

  USE mo_alloc_patches,       ONLY: destruct_patches
  USE mo_ocean_read_namelists, ONLY: read_ocean_namelists
  USE mo_io_restart,          ONLY: read_restart_info_file, read_restart_files
  USE mo_io_restart_namelist, ONLY: read_restart_namelists
  USE mo_io_restart_attributes,ONLY: read_restart_attributes, get_restart_attribute

  USE mo_time_config,         ONLY: time_config
  USE mo_icon_comm_interface, ONLY: construct_icon_communication, destruct_icon_communication
  
  !-------------------------------------------------------------
  ! For the coupling
#ifndef __ICON_OCEAN_ONLY__
  USE mo_icon_cpl_init,       ONLY: icon_cpl_init
  USE mo_icon_cpl_init_comp,  ONLY: icon_cpl_init_comp
  USE mo_coupling_config,     ONLY: is_coupled_run, config_debug_coupler_level
  USE mo_icon_cpl_def_grid,   ONLY: icon_cpl_def_grid, icon_cpl_def_location
  USE mo_icon_cpl_def_field,  ONLY: icon_cpl_def_field
  USE mo_icon_cpl_search,     ONLY: icon_cpl_search
  USE mo_icon_cpl_finalize,   ONLY: icon_cpl_finalize
#endif
  !-------------------------------------------------------------
  
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC :: ocean_model
  PUBLIC :: construct_ocean_model, destruct_ocean_model
  PUBLIC :: ocean_patch_3d, ocean_state, operators_coefficients

  LOGICAL :: l_have_output
  TYPE(t_patch_3d), POINTER :: ocean_patch_3d
  TYPE(t_atmos_for_ocean)                         :: p_as
  TYPE(t_atmos_fluxes)                            :: p_atm_f
  TYPE(t_operator_coeff)                          :: operators_coefficients
  TYPE(t_hydro_ocean_state), ALLOCATABLE, TARGET  :: ocean_state(:)
  TYPE(t_datetime)                                :: start_datetime

  
CONTAINS


  !--------------------------------------------------------------------------
  !>
  !!
  SUBROUTINE ocean_model(oce_namelist_filename,shr_namelist_filename)
    
    CHARACTER(LEN=*), INTENT(in) :: oce_namelist_filename,shr_namelist_filename
    
    CHARACTER(*), PARAMETER :: routine = "mo_ocean_model:ocean_model"
    
    INTEGER :: n_io, jg, jfile, ist
    INTEGER :: error_status

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
    n_io    = n_ios()   ! number of: write output

    !------------------------------------------------------------------
    ! Initialize output file if necessary;
    ! Write out initial conditions.
    !------------------------------------------------------------------

    IF (output_mode%l_nml) THEN
      WRITE(0,*)'process_mpi_io_size:',process_mpi_io_size
      IF (process_mpi_io_size > 0) use_async_name_list_io = .TRUE.
      CALL parse_variable_groups()
      CALL init_name_list_output(lprintlist=.TRUE.,l_is_ocean=.TRUE.)
    ENDIF

    IF (.NOT.is_restart_run()) THEN
      ! Initialize the first output file which will contain also the
      ! initial conditions.
      jfile = 1
      CALL init_output_files(jfile, lclose=.FALSE.,p_patch_2d=ocean_patch_3d%p_patch_2d)
      IF (lwrite_initial)       CALL init_output_files(jfile, lclose=.FALSE.,p_patch_2d=ocean_patch_3d%p_patch_2d)
      IF (lwrite_initial) THEN
        !IF (output_mode%l_nml) THEN
        !  CALL write_name_list_output( time_config%cur_datetime, 0._wp, .FALSE. )
        !ENDIF
        IF (output_mode%l_vlist) THEN
          CALL write_output_oce( time_config%cur_datetime,z_sim_time=(/1.0_wp/),&
            & p_patch_3d=ocean_patch_3d, p_os=ocean_state)
        ENDIF
      ENDIF
      l_have_output = .TRUE.
    ELSE
      ! No need to write out the initial condition, thus no output
      ! during the first integration step. This run will produce
      ! output if n_io <= integration_length.

      CALL get_restart_attribute('next_output_file',jfile)
      !      IF (n_io.le.nsteps) THEN
      CALL init_output_files(jfile, lclose=.FALSE.,p_patch_2d=ocean_patch_3d%p_patch_2d)
      l_have_output = .TRUE.
      !      ELSE
      !         l_have_output = .FALSE.
      !      END IF

    END IF ! (not) is_restart_run()

    !------------------------------------------------------------------
    CALL perform_ho_stepping( ocean_patch_3d, ocean_state,                    &
      & ext_data, start_datetime, n_io,                     &
      & jfile,                                        &
      & (nsteps == INT(time_config%dt_restart/dtime)),&
      & v_sfc_flx,                                    &
      & v_params, p_as, p_atm_f,v_sea_ice,operators_coefficients,&
      & l_have_output)
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

    ! Delete variable lists

    IF (output_mode%l_nml) THEN
      CALL close_name_list_output
    ENDIF

    IF (output_mode%l_vlist) THEN
      IF (l_have_output) THEN
        CALL close_output_files
      ENDIF
    ENDIF

    CALL destruct_icon_communication()

#ifndef __ICON_OCEAN_ONLY__
    IF ( is_coupled_run() ) CALL icon_cpl_finalize
#endif

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
    ! If we belong to the I/O PEs just call vlist_io_main_proc before reading patches.
    ! This routine will never return

    IF (my_process_is_io()) CALL vlist_io_main_proc

    CALL build_decomposition(num_lev,num_levp1,nshift, is_ocean_decomposition =.TRUE., &
      & l_restore_states=lrestore_states, patch_3d=ocean_patch_3d)
    CALL construct_icon_communication(ocean_patch_3d%p_patch_2d(:), n_dom=1)

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

#ifndef __ICON_OCEAN_ONLY__
    IF ( is_coupled_run() ) THEN
      CALL construct_ocean_coupling()
    ENDIF
#endif

    ! Prepare time integration
    CALL prepare_ho_integration(ocean_patch_3d, ocean_state, ext_data, v_sfc_flx, &
      & v_params, p_as, p_atm_f, v_sea_ice,operators_coefficients)!,p_int_state(1:))


  END SUBROUTINE construct_ocean_model
  !--------------------------------------------------------------------------

  !--------------------------------------------------------------------------
#ifndef __ICON_OCEAN_ONLY__
  !------------------------------------------------------------------
  ! Prepare the coupling
  !
  ! For the time being this could all go into a subroutine which is
  ! common to atmo and ocean. Does this make sense if the setup deviates
  ! too much in future.
  !------------------------------------------------------------------
  SUBROUTINE construct_ocean_coupling()

    INTEGER, PARAMETER :: no_of_fields = 10
    CHARACTER(LEN=max_char_length) ::  field_name(no_of_fields)
    INTEGER :: field_id(no_of_fields)
    INTEGER :: grid_id
    INTEGER :: grid_shape(2)
    INTEGER :: field_shape(3)
    INTEGER :: i, error_status

    INTEGER :: patch_no

    !------------------------------------------------------------
    CALL icon_cpl_init(debug_level=config_debug_coupler_level)
    ! Inform the coupler about what we are
    CALL icon_cpl_init_comp ( get_my_process_name(), get_my_model_no(), error_status )
    ! split the global_mpi_communicator into the components
    !------------------------------------------------------------
    patch_no      = 1

    grid_shape(1) = 1
    grid_shape(2) = ocean_patch_3d%p_patch_2d(patch_no)%n_patch_cells

    CALL icon_cpl_def_grid ( &
      & grid_shape, ocean_patch_3d%p_patch_2d(patch_no)%cells%glb_index, & ! input
      & grid_id, error_status )                          ! output

    ! Marker for internal and halo points, a list which contains the
    ! rank where the native cells are located.
    CALL icon_cpl_def_location ( &
      & grid_id, grid_shape, ocean_patch_3d%p_patch_2d(patch_no)%cells%owner_local, & ! input
      & p_pe_work,  & ! this owner id
      & error_status )                                            ! output

    field_name(1) = "TAUX"
    field_name(2) = "TAUY"
    field_name(3) = "SFWFLX" ! bundled field containing two components
    field_name(4) = "SFTEMP"
    field_name(5) = "THFLX"  ! bundled field containing two components
    field_name(6) = "ICEATM" ! bundled field containing four components
    field_name(7) = "SST"
    field_name(8) = "OCEANU"
    field_name(9) = "OCEANV"
    field_name(10) = "ICEOCE" ! bundled field containing four components

    field_shape(1:2) = grid_shape(1:2)

    DO i = 1, no_of_fields
      IF ( i == 3 .OR. i == 5 ) THEN
        field_shape(3) = 2
      ELSE IF ( i == 6 .OR. i == 10 ) THEN
        field_shape(3) = 4
      ELSE
        field_shape(3) = 1
      ENDIF
      CALL icon_cpl_def_field ( field_name(i), grid_id, field_id(i), &
        & field_shape, error_status )
    ENDDO

    CALL icon_cpl_search

  END SUBROUTINE construct_ocean_coupling
#endif
  !--------------------------------------------------------------------------

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

