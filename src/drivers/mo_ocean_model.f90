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

  USE mo_exception,           ONLY: message, finish  ! use always
  USE mo_master_control,      ONLY: is_restart_run, get_my_process_name, get_my_model_no
  USE mo_parallel_config,     ONLY: p_test_run, l_test_openmp, num_io_procs
  USE mo_mpi,                 ONLY: p_stop, &
    & my_process_is_io,  my_process_is_mpi_seq, my_process_is_mpi_test, &
    & my_process_is_mpi_parallel,                                       &
    & set_mpi_work_communicators, set_comm_input_bcast, null_comm_type, &
    & p_pe_work
  USE mo_sync,                ONLY: enable_sync_checks, disable_sync_checks
  USE mo_timer,               ONLY: init_timer, timer_start, timer_stop, print_timer, &
    &                               timer_model_init
  USE mo_datetime,            ONLY: t_datetime
  USE mo_output,              ONLY: init_output_files, write_output, close_output_files
  USE mo_grid_config,         ONLY: n_dom, n_dom_start, global_cell_type !, &
!                                  dynamics_parent_grid_id
  USE mo_dynamics_config,     ONLY: iequations

  USE mo_io_async,            ONLY: vlist_io_main_proc            ! main procedure for I/O PEs

  USE mo_interpol_config,     ONLY: configure_interpolation
  USE mo_advection_config,    ONLY: configure_advection
  USE mo_dynamics_config,     ONLY: configure_dynamics  ! subroutine
  USE mo_run_config,          ONLY: configure_run

  ! Control parameters: run control, dynamics, i/o
  !
  USE mo_io_config,           ONLY:  dt_data,dt_file,dt_diag,lwrite_initial
  USE mo_io_config,           ONLY:  n_diags, n_checkpoints,n_files,n_ios
  USE mo_run_config,          ONLY: &
    & dtime,                & !    :
    & nsteps,                & !    :
!    & ltransport,           & !    :
    & ltimer,               & !    :
    & ldump_states,         & ! flag if states should be dumped
    & lrestore_states,      &  ! flag if states should be restored
    & iforcing,             & !  
    & nlev, nlevp1,         & !
    & num_lev, num_levp1,   &
    & iqv, nshift,          &
    & ntracer
  USE mo_nml_crosscheck,    ONLY: oce_crosscheck

!  USE mo_advection_nml,       ONLY: transport_nml_setup,  & ! process transport
!    & setup_transport         ! control parameters

  USE mo_subdivision,         ONLY: decompose_domain,         &
    & complete_parallel_setup,       &
    & finalize_decomposition,        &
    & copy_processor_splitting,      &
    & set_patch_communicators
  USE mo_dump_restore,        ONLY: dump_patch_state_netcdf,       &
    & restore_patches_netcdf,        &
    & restore_interpol_state_netcdf

  USE mo_icoham_dyn_memory,   ONLY: p_hydro_state
  USE mo_model_domain,        ONLY: t_patch, p_patch, p_patch_local_parent
  USE mo_intp_data_strc,      ONLY: p_int_state, p_int_state_local_parent
  USE mo_grf_intp_data_strc,  ONLY: p_grf_state, p_grf_state_local_parent

  ! Horizontal grid
  !
  USE mo_grid_config,         ONLY: &!  grid_nml_setup,          & ! process grid control parameters
    & n_dom,                & !    :
    & n_dom_start
  USE mo_model_domimp_patches,ONLY: &
    & import_basic_patches, & !
    & complete_patches

  ! Horizontal interpolation
  !
!  USE mo_interpol_nml,        ONLY: interpol_nml_setup   ! process interpol. ctl. params.
  USE mo_intp_state,          ONLY: construct_2d_interpol_state, &
    & destruct_2d_interpol_state, transfer_interpol_state

  USE mo_oce_state,           ONLY: t_hydro_ocean_state, v_ocean_state

!0  USE mo_mpiom_phy_state,     ONLY: construct_mpiom_phy_state, &
!0    &                               destruct_mpiom_phy_state

  USE mo_impl_constants,      ONLY: success !, ihs_ocean

  ! External data
  USE mo_ext_data_state,       ONLY: ext_data, init_ext_data, destruct_ext_data

  USE mo_hydro_ocean_run,      ONLY: perform_ho_stepping,&
    & prepare_ho_integration,&
    & finalise_ho_integration

!   USE mo_oce_forcing,         ONLY: t_sfc_flx, t_atmos_fluxes, t_atmos_for_ocean, &
!     &                               v_sfc_flx
  USE mo_oce_physics,         ONLY: t_ho_params, v_params!, t_ho_physics
! #
  USE mo_sea_ice_types,       ONLY: t_sfc_flx, t_atmos_fluxes, t_atmos_for_ocean, &
    &                               v_sfc_flx, v_sea_ice
  ! For the coupling
  USE mo_icon_cpl_init,      ONLY: icon_cpl_init
  USE mo_icon_cpl_init_comp, ONLY: icon_cpl_init_comp
  USE mo_impl_constants,      ONLY: CELLS, MAX_CHAR_LENGTH
  USE mo_coupling_config,     ONLY : is_coupled_run, config_debug_coupler_level
  USE mo_icon_cpl_def_grid,   ONLY : ICON_cpl_def_grid, ICON_cpl_def_location
  USE mo_icon_cpl_def_field,  ONLY : ICON_cpl_def_field
  USE mo_icon_cpl_search,     ONLY : ICON_cpl_search
  USE mo_icon_cpl_finalize,   ONLY: icon_cpl_finalize
  USE mo_alloc_patches,       ONLY : destruct_patches

  !-------------------------------------------------------------
  USE mo_read_namelists,       ONLY: read_ocean_namelists
  USE mo_io_restart,           ONLY: read_restart_info_file, read_restart_files
  USE mo_io_restart_namelist,  ONLY: read_restart_namelists
  USE mo_io_restart_attributes,ONLY: read_restart_attributes, get_restart_attribute

  USE mo_time_config,         ONLY: time_config

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: ocean_model

CONTAINS
  !>
  !! 
  SUBROUTINE ocean_model(oce_namelist_filename,shr_namelist_filename)

    CHARACTER(LEN=*), INTENT(in) :: oce_namelist_filename,shr_namelist_filename

    CHARACTER(*), PARAMETER :: routine = "mo_ocean_model:ocean_model"

    TYPE(t_atmos_for_ocean)                         :: p_as
    TYPE(t_atmos_fluxes)                            :: p_atm_f

    TYPE(t_datetime)                                :: datetime

    INTEGER :: n_io, jg, jgp, jfile, ist

    TYPE(t_patch), ALLOCATABLE :: p_patch_global(:)

    ! For the coupling

    INTEGER, PARAMETER :: no_of_fields = 8

    CHARACTER(LEN=MAX_CHAR_LENGTH) ::  field_name(no_of_fields)
    CHARACTER(LEN=MAX_CHAR_LENGTH) :: grid_file_name
    INTEGER :: field_id(no_of_fields)
    INTEGER :: grid_id
    INTEGER :: grid_shape(2) 
    INTEGER :: field_shape(3) 
    INTEGER :: i, error_status
!    INTEGER :: no_of_entities
    INTEGER :: patch_no
!    INTEGER, POINTER :: grid_glob_index(:)
    LOGICAL :: lsuccess, l_have_output

    !-------------------------------------------------------------------

    IF (is_restart_run()) THEN

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

      CALL read_restart_namelists('restart_oce.nc')
      CALL message(TRIM(routine), 'read namelists from ocean restart file')

      ! Read all global attributs in the restart file and store them in a buffer.

      CALL read_restart_attributes('restart_oce.nc')
      CALL message(TRIM(routine), 'read global attributes from ocean restart file')

    END IF ! is_restart_run()

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
    CALL set_mpi_work_communicators(p_test_run, l_test_openmp, num_io_procs)

    !-------------------------------------------------------------------
    ! 3.2 Initialize various timers
    !-------------------------------------------------------------------
    IF (ltimer) CALL init_timer

    IF (ltimer) CALL timer_start(timer_model_init)

    !-------------------------------------------------------------------
    ! Initialize date and time
    !-------------------------------------------------------------------
    datetime = time_config%cur_datetime

    !-------------------------------------------------------------------
    ! 4. Import patches
    !-------------------------------------------------------------------
    ! If we belong to the I/O PEs just call vlist_io_main_proc before reading patches.
    ! This routine will never return

    IF (my_process_is_io()) CALL vlist_io_main_proc

    ! Check patch allocation status
    IF ( ALLOCATED(p_patch)) THEN
      CALL finish(TRIM(routine), 'patch already allocated')
    END IF

    ! Allocate patch array to start patch construction
    ALLOCATE(p_patch(n_dom_start:n_dom), stat=error_status)
    IF (error_status/=success) THEN
      CALL finish(TRIM(routine), 'allocation of patch failed')
    ENDIF

    IF(lrestore_states) THEN
      ! Before the restore set p_comm_input_bcast to null
      CALL set_comm_input_bcast(null_comm_type)
      IF( .NOT. my_process_is_mpi_test()) THEN
        CALL restore_patches_netcdf( p_patch, .TRUE. )
        CALL set_patch_communicators(p_patch)
      ELSE
        CALL import_basic_patches(p_patch,nlev,nlevp1,num_lev,num_levp1,nshift)
        CALL disable_sync_checks
        CALL complete_patches( p_patch )
        CALL enable_sync_checks
      ENDIF
      ! After the restore is done set p_comm_input_bcast in the
      ! same way as it would be set in parallel_nml_setup when
      ! no restore is wanted:
      CALL set_comm_input_bcast()

    ELSE

      ! Please note: ldump_dd/lread_dd not (yet?) implemented
      IF(my_process_is_mpi_parallel()) THEN
        ALLOCATE(p_patch_global(n_dom_start:n_dom))
        CALL import_basic_patches(p_patch_global,nlev,nlevp1,num_lev,num_levp1,nshift)      
        CALL decompose_domain(p_patch_global)
        DEALLOCATE(p_patch_global)
        CALL complete_parallel_setup()
      ELSE
        CALL import_basic_patches(p_patch,nlev,nlevp1,num_lev,num_levp1,nshift)
      ENDIF

      ! Complete information which is not yet read or calculated
      CALL complete_patches( p_patch )

    ENDIF

    ! In case of a test run: Copy processor splitting to test PE
    IF(p_test_run) CALL copy_processor_splitting(p_patch)

    !--------------------------------------------------------------------------------
    ! 5. Construct interpolation state, compute interpolation coefficients.
    !--------------------------------------------------------------------------------

    CALL configure_interpolation( global_cell_type, n_dom, p_patch(1:)%level )

    ! Allocate array for interpolation state

    ALLOCATE( p_int_state(n_dom_start:n_dom), &
            & p_grf_state(n_dom_start:n_dom),STAT=error_status)
    IF (error_status /= SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for ptr_int_state failed')
    ENDIF

    IF(my_process_is_mpi_parallel()) THEN
      ALLOCATE( p_int_state_local_parent(n_dom_start+1:n_dom), &
              & p_grf_state_local_parent(n_dom_start+1:n_dom), &
              & STAT=error_status)
      IF (error_status /= SUCCESS) &
        CALL finish(TRIM(routine),'allocation for local parents failed')
    ENDIF

    IF(lrestore_states) THEN !TODOram
      ! Interpolation state is read from NetCDF
      ! On the test PE it is constructed at the same time to be able to check
      ! the state read in with the constructed state
      IF( .NOT. my_process_is_mpi_test()) THEN
        CALL restore_interpol_state_netcdf(p_patch, p_int_state)
      ELSE
        ! construct_2d_interpol_state makes sync calls for checking the
        ! results on the parallel PEs to the results on the test PE.
        ! These checks must be disabled here!
        CALL disable_sync_checks
        CALL construct_2d_interpol_state(p_patch, p_int_state)
        CALL enable_sync_checks
      ENDIF
    ELSE
      ! Construct interpolation state
      ! Please note that for pararllel runs the divided state is constructed here
      CALL construct_2d_interpol_state(p_patch, p_int_state)
      IF(my_process_is_mpi_parallel()) THEN
        ! Transfer interpolation state to local parent
        DO jg = n_dom_start+1, n_dom
          jgp = p_patch(jg)%parent_id
          CALL transfer_interpol_state(p_patch(jgp),p_patch_local_parent(jg), &
                                    &  p_int_state(jgp), p_int_state_local_parent(jg))
        ENDDO
      ENDIF
    ENDIF


    !------------------------------------------------------------------
    ! step 5b: allocate state variables
    !------------------------------------------------------------------
    ALLOCATE (v_ocean_state(n_dom), stat=ist)
    IF (ist /= success) THEN
      CALL finish(TRIM(routine),'allocation for v_ocean_state failed')
    ENDIF

    !-------------------------------------------------------------------
    ! 7. Finalize domain decomposition
    !-------------------------------------------------------------------   
    IF (my_process_is_mpi_parallel() .AND. .NOT.lrestore_states) THEN

      CALL finalize_decomposition()

    ENDIF

    IF(ldump_states)THEN

      ! Dump divided patches with interpolation and grf state to NetCDF file and exit

      CALL message(TRIM(routine),'ldump_states is set: '//&
                  'dumping patches+states and finishing')

      IF(.NOT. my_process_is_mpi_test()) THEN
        DO jg = n_dom_start, n_dom
          CALL dump_patch_state_netcdf(p_patch(jg),p_int_state(jg),p_grf_state(jg))
        ENDDO
      ENDIF

      CALL p_stop
      STOP

    ENDIF

    !---------------------------------------------------------------------
    ! 9. Horizontal and vertical grid(s) are now defined.
    !    Assign values to derived variables in the configuration states
    !---------------------------------------------------------------------

    CALL configure_dynamics ( n_dom )
!     CALL configure_diffusion( n_dom, dynamics_parent_grid_id, &
!                             & nlev, vct_a, vct_b, apzero      )

    DO jg =1,n_dom
      CALL configure_advection( jg, p_patch(jg)%nlev, p_patch(1)%nlev, &
        &                      iequations, iforcing, iqv, 0, &
        &                      0, 1, .false., .true., ntracer )
    ENDDO


    !------------------------------------------------------------------
    ! 10. Create and optionally read external data fields
    !------------------------------------------------------------------
    ALLOCATE (ext_data(n_dom), STAT=error_status)
    IF (error_status /= SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for ext_data failed')
    ENDIF

    ! allocate memory for oceanic external data and
    ! optionally read those data from netCDF file.
    CALL init_ext_data (p_patch(1:), p_int_state(1:), ext_data)

    !------------------------------------------------------------------
    ! Prepare the coupling
    !
    ! For the time being this could all go into a subroutine which is
    ! common to atmo and ocean. Does this make sense if the setup deviates
    ! too much in future.
    !------------------------------------------------------------------

    IF ( is_coupled_run() ) THEN

     !------------------------------------------------------------
      CALL icon_cpl_init(debug_level=config_debug_coupler_level)
      ! Inform the coupler about what we are
      CALL icon_cpl_init_comp ( get_my_process_name(), get_my_model_no(), error_status )
      ! split the global_mpi_communicator into the components
     !------------------------------------------------------------
      patch_no      = 1

      grid_shape(1) = 1
      grid_shape(2) = p_patch(patch_no)%n_patch_cells

      ! CALL get_patch_global_indexes ( patch_no, CELLS, no_of_entities, grid_glob_index )
      ! should grid_glob_index become a pointer in ICON_cpl_def_grid as well?
      CALL ICON_cpl_def_grid ( &
        & grid_shape, p_patch(patch_no)%cells%glb_index, & ! input
        & grid_id, error_status )                          ! output


      ! Marker for internal and halo points, a list which contains the
      ! rank where the native cells are located.
      ! p_patch(patch_no)%cells%owner_local(:) = 0 ! the ocean run sequentially
      CALL ICON_cpl_def_location ( &
        & grid_id, grid_shape, p_patch(patch_no)%cells%owner_local, & ! input
        & p_pe_work,  & ! this owner id
        & error_status )                                            ! output

      field_name(1) = "TAUX"
      field_name(2) = "TAUY"
      field_name(3) = "SFWFLX" ! bundled field containing two components
      field_name(4) = "SFTEMP"
      field_name(5) = "THFLX"  ! bundled field containing four components
      field_name(6) = "SST"
      field_name(7) = "OCEANU"
      field_name(8) = "OCEANV"

      field_shape(1:2) = grid_shape(1:2)

      DO i = 1, no_of_fields
      IF ( i == 3 ) THEN
         field_shape(3) = 2
      ELSE IF ( i == 5 ) THEN
         field_shape(3) = 4
      ELSE
         field_shape(3) = 1
      ENDIF
         CALL ICON_cpl_def_field ( field_name(i), grid_id, field_id(i), &
    &                               field_shape, error_status )
      ENDDO

      CALL ICON_cpl_search

    ENDIF

    ! Prepare time integration
    CALL prepare_ho_integration(p_patch(1:), v_ocean_state, ext_data, v_sfc_flx, &
      &                         v_params, p_as, p_atm_f, v_sea_ice)

    !------------------------------------------------------------------
    ! Daniel: Suggestion for point 5 of Feature #333
    ! (5. Subroutine to setup model components depending on this
    ! namelist, to be called after all namelists have been read, and
    ! a synoptic check has been done)
    !------------------------------------------------------------------
    ! set dependent variables/model components, depending on this (transport)
    ! namelist and potentially others
!    IF (ltransport) THEN
!      CALL setup_transport( ihs_ocean )
!    ENDIF
!
    !------------------------------------------------------------------
    ! Set initial conditions for time integration.
    !------------------------------------------------------------------

    IF (is_restart_run()) THEN
    ! This is an resumed integration. Read model state from restart file(s).
#ifdef NOMPI
      CALL read_restart_files
#else
      jg = 1 !no nesting
     !DO jg = n_dom_start,n_dom
     CALL read_restart_files( p_patch(jg) )
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

    IF (.NOT.is_restart_run()) THEN
      ! Initialize the first output file which will contain also the
      ! initial conditions.
      jfile = 1
      CALL init_output_files(jfile, lclose=.FALSE.)
      IF (lwrite_initial) CALL write_output( time_config%cur_datetime )
      l_have_output = .TRUE.
    ELSE
    ! No need to write out the initial condition, thus no output
    ! during the first integration step. This run will produce
    ! output if n_io <= integration_length.

      CALL get_restart_attribute('next_output_file',jfile)
      IF (n_io.le.nsteps) THEN
         CALL init_output_files(jfile, lclose=.FALSE.)
         l_have_output = .TRUE.
      ELSE
         l_have_output = .FALSE.
      END IF

    END IF ! (not) is_restart_run()


    IF (ltimer) CALL timer_stop(timer_model_init)

    CALL perform_ho_stepping( p_patch(1:), v_ocean_state,                   &
      &                       ext_data, datetime, n_io,                     &
      &                       jfile,                                        &
      &                       (nsteps == INT(time_config%dt_restart/dtime)),&
      &                       p_int_state(1:),                              &
      &                       v_sfc_flx,                                    &
      &                       v_params, p_as, p_atm_f,v_sea_ice,  &
      &                       l_have_output)

    IF (ltimer) CALL print_timer

    !------------------------------------------------------------------
    !  cleaning up process
    !------------------------------------------------------------------
    CALL message(TRIM(routine),'start to clean up')

    CALL finalise_ho_integration(v_ocean_state, v_params, &
      &                          p_as, p_atm_f, v_sea_ice, v_sfc_flx)



    !---------------------------------------------------------------------
    ! 13. Integration finished. Carry out the shared clean-up processes
    !---------------------------------------------------------------------
    ! Destruct external data state
    CALL destruct_ext_data

    ! deallocate ext_data array
    DEALLOCATE(ext_data, stat=error_status)
    IF (error_status/=success) THEN
      CALL finish(TRIM(routine), 'deallocation of ext_data')
    ENDIF

    ! Deallocate interpolation fields
    ! interpolation state not used for ocean model
    ! #slo# - temporarily switched on for comparison with rbf-reconstruction
    CALL destruct_2d_interpol_state( p_int_state )
    DEALLOCATE (p_int_state, stat=ist)
    IF (ist /= success) THEN
      CALL finish(TRIM(routine),'deallocation for ptr_int_state failed')
    ENDIF

    ! Deallocate grid patches

    CALL destruct_patches( p_patch )

    DEALLOCATE( p_patch, stat=ist )
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocate for patch array failed')
    ENDIF

    ! Delete variable lists
    CALL close_output_files

    IF ( is_coupled_run() ) CALL ICON_cpl_finalize
    
    CALL message(TRIM(routine),'clean-up finished')

  END SUBROUTINE ocean_model

END MODULE mo_ocean_model

