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
  USE mo_mpi,                 ONLY: p_stop, p_pe, p_io, p_nprocs
  USE mo_timer,               ONLY: init_timer, print_timer
  USE mo_namelist,            ONLY: open_nml,  close_nml, open_nml_output, close_nml_output
  USE mo_datetime,            ONLY: t_datetime
  USE mo_output,              ONLY: init_output_files
  USE mo_io_vlist,            ONLY: write_vlist_oce, destruct_vlist_oce
  
  USE mo_parallel_configuration,        ONLY:   & ! process parallel run ctl. params.
    & p_comm_work_test, p_comm_input_bcast, & ! communicators
    & p_test_pe,            & !    internal parameter
    & p_comm_work,          &
    & p_test_run,           &
    & p_io_pe0                ! Number of first I/O PE
  
  USE mo_io_async,            ONLY: io_main_proc            ! main procedure for I/O PEs
  
  
  ! Control parameters: run control, dynamics, i/o
  !
  USE mo_global_variables,    ONLY: setup_physics        ! process forcing control parameters
!0  &                               impiom,            & !    :
  USE mo_ocean_nml,           ONLY: setup_ocean_nml
  USE mo_dynamics_nml,        ONLY: dynamics_nml_setup,   &
    & cleanup_dyn_params 
  USE mo_diffusion_nml,       ONLY: diffusion_nml_setup
  USE mo_io_nml,              ONLY: io_nml_setup,         & ! process I/O
    & dt_data,              & !    :
    & dt_file,              & !    :
    & dt_diag,              & !    :
    & lprepare_output         ! internal parameter
  USE mo_run_nml,             ONLY: run_nml_setup,            & ! process run control parameters
    & dtime,                & !    :
    & ltransport,           & !    :
    & lforcing,             & !    :
    & ltimer,               & !    :
!0  & locean,               & !    :
!0  & iequations,           & !    internal parameters
!0  & ihs_ocean,            & !
!0  & impiom,               & !    :
    & ldump_states,         & ! flag if states should be dumped
    & lrestore_states         ! flag if states should be restored
  
  USE mo_advection_nml,       ONLY: transport_nml_setup,  & ! process transport
    & setup_transport         ! control parameters
  
  USE mo_subdivision,         ONLY: decompose_atmo_domain,         &
    & copy_processor_splitting,      &
    & set_patch_communicators
  USE mo_dump_restore,        ONLY: dump_patch_state_netcdf,       &
    & restore_patches_netcdf,        &
    & restore_interpol_state_netcdf
  
  USE mo_icoham_dyn_memory,   ONLY: p_hydro_state
  USE mo_atmo_control,        ONLY: p_patch_global, p_patch_subdiv, p_patch,             &
    & p_int_state_global, p_int_state_subdiv, p_int_state, &
    & p_grf_state_global, p_grf_state
  
  ! Horizontal grid
  !
  USE mo_model_domain_import, ONLY: &!  grid_nml_setup,          & ! process grid control parameters
    & n_dom,                & !    :
    & n_dom_start,          & !    :
    & parent_id,            & !    :
    & import_patches,       & !
    & destruct_patches        !
  
  ! Horizontal interpolation
  !
  USE mo_interpol_nml,        ONLY: interpol_nml_setup   ! process interpol. ctl. params.
  USE mo_intp_state,          ONLY: construct_2d_interpol_state, &
    & destruct_2d_interpol_state
  
  USE mo_oce_state,           ONLY: t_hydro_ocean_state, destruct_hydro_ocean_state  
  USE mo_gmt_output,          ONLY: setup_gmt_output
  
!0  USE mo_mpiom_phy_state,     ONLY: construct_mpiom_phy_state, &
!0    &                               destruct_mpiom_phy_state
  
  USE mo_impl_constants,      ONLY: success
  
  ! External data
  USE mo_ext_data,            ONLY: ext_data
  
!   USE mo_test_hydro_ocean,    ONLY: prepare_ho_integration, finalise_ho_integration!, &
!     &                               test_hydro_ocean
  USE mo_hydro_ocean_run,      ONLY: perform_ho_stepping,&
    & prepare_ho_integration,&
    & finalise_ho_integration
  
  USE mo_oce_forcing,         ONLY: t_ho_sfc_flx
  USE mo_oce_physics,         ONLY: t_ho_params, t_ho_physics

  !-------------------------------------------------------------

  USE mo_time_config, ONLY: time_config

  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC :: ocean_model
  
CONTAINS
  !>
  !! 
  SUBROUTINE ocean_model(namelist_filename)
    
    CHARACTER(LEN=*), INTENT(in) :: namelist_filename

    CHARACTER(*), PARAMETER :: routine = "mo_ocean_model:ocean_model"

    ! Ad-hoc declaration for the hydrostatic ocean state variables.
    ! This declaration and the declarations for the atmospheric model states, see
    ! "USE mo_atmo_control" above, should eventually be organized in the same way.
    TYPE(t_hydro_ocean_state), TARGET, ALLOCATABLE, SAVE :: p_hyoce_state(:)
    TYPE(t_ho_sfc_flx)                              :: p_sfc_flx
    TYPE (t_ho_params)                              :: p_phys_param
    TYPE(t_ho_physics)                              :: p_physics_oce
    TYPE(t_datetime)                                :: datetime
    
    INTEGER :: n_io, jg, jfile, n_file, ist, n_diag
    !-------------------------------------------------------------------
    
    CALL message(TRIM(routine),'start model initialization.')
    
    !-------------------------------------------------------------------
    ! step 1: Open the namelist file and read most basic namelists.
    !-------------------------------------------------------------------
    
    CALL open_nml(TRIM(namelist_filename))
    
    ! Create a new file in which all the namelist variables and their
    ! actual values used in the model run will be stored.
    
    IF(p_pe == p_io) CALL open_nml_output('NAMELIST_ICON_output')
    
    ! Read namelists 'run_ctl' and 'ocean_ctl' 
    
    CALL run_nml_setup
    CALL setup_ocean_nml
    
    !-------------------------------------------------------------------
    ! Initialize timers
    !-------------------------------------------------------------------
    IF (ltimer) CALL init_timer
    
    ! parallel_nml_setup must be called after setup_run since it needs
    ! some variables read in setup_run
    
!     CALL parallel_nml_setup
    
    !-------------------------------------------------------------------
    ! Initialize date and time
    !-------------------------------------------------------------------
    datetime = time_config%ini_datetime
    
    !-------------------------------------------------------------------
    ! step 2: import computational domain of the model
    !-------------------------------------------------------------------
    ! 2a) read namelist 'grid_ctl' from the namelist file (already
    !     opened above) and set up the grid/patch configuration.
    !-------------------------------------------------------------------
  ! CALL grid_nml_setup
    
    !-------------------------------------------------------------------
    ! 2b) patch import
    !-------------------------------------------------------------------
    ! If we belong to the I/O PEs just call io_main_proc before reading patches.
    ! This routine will never return
    
    IF(p_pe >= p_io_pe0) CALL io_main_proc
    !-------------------------------------------------------------------
    
    !check patch allocation status
    IF ( ALLOCATED(p_patch_global)) THEN
      CALL finish(TRIM(routine), 'patch already allocated')
    END IF
    !
    ! allocate patch array to start patch construction
    ALLOCATE(p_patch_global(n_dom_start:n_dom), stat=ist)
    IF (ist/=success) THEN
      CALL finish(TRIM(routine), 'allocation of patch failed')
    ENDIF
    
    IF(lrestore_states .AND. p_pe /= p_test_pe) THEN
      CALL restore_patches_netcdf( p_patch_global )
    ELSE
      CALL import_patches( p_patch_global )
    ENDIF
    
    IF(lrestore_states) THEN
      ! After the restore is done set p_comm_input_bcast in the
      ! same way as it would be set in parallel_nml_setup when
      ! no restore is wanted:
      IF(p_test_run) THEN
        ! Test PE reads and broadcasts to workers
        p_comm_input_bcast = p_comm_work_test
      ELSE
        ! PE 0 reads and broadcasts
        p_comm_input_bcast = p_comm_work
      ENDIF
    ENDIF
    
    !------------------------------------------------------------------
    ! step 3a: Read namelist 'dynamics_ctl'. This has to be done
    !          before the interpolation state is computed.
    !------------------------------------------------------------------
    CALL dynamics_nml_setup(n_dom)
    
    !------------------------------------------------------------------
    ! step 3b: Read namelist 'transport_ctl',
    !------------------------------------------------------------------
    IF (ltransport) THEN
      CALL transport_nml_setup
    END IF
    
    !------------------------------------------------------------------
    ! step 4: construct interpolation and the grid refinement state
    !------------------------------------------------------------------
    ! - Allocate memory for the interpolation state;
    ! - Calculate interpolation coefficients.
    !------------------------------------------------------------------
    ! interpolation state not used for ocean model
    ! #slo# - temporarily switched on for comparison with rbf-reconstruction
    !IF (iequations /= ihs_ocean) THEN
    
    CALL interpol_nml_setup(p_patch_global)
    
    ! allocate type for interpolation state
    
    ALLOCATE (p_int_state_global(n_dom_start:n_dom), &
      & p_grf_state_global(n_dom_start:n_dom),stat=ist)
    IF (ist /= success) THEN
      CALL finish(TRIM(routine),'allocation for ptr_int_state failed')
    ENDIF
    
    IF(lrestore_states .AND. p_pe /= p_test_pe) THEN
      ! Read interpolation state from NetCDF
      CALL restore_interpol_state_netcdf(p_patch_global, p_int_state_global)
    ELSE
      ! Interpolation state is constructed for
      ! the full domain on every PE and divided later
      CALL construct_2d_interpol_state(p_patch_global, p_int_state_global)
    ENDIF
    
    !END IF
    
    !-------------------------------------------------------------------------
    ! set up horizontal diffusion after the vertical coordinate is configured
    !-------------------------------------------------------------------------
    CALL diffusion_nml_setup(n_dom,parent_id)
    
    !------------------------------------------------------------------
    ! step 5b: allocate state variables
    !------------------------------------------------------------------
    
    ALLOCATE (p_hyoce_state(n_dom), stat=ist)
    IF (ist /= success) THEN
      CALL finish(TRIM(routine),'allocation for p_hyoce_state failed')
    ENDIF
    ! #slo# - temporarily switched on for comparison with rbf-reconstruction
    ALLOCATE (p_hydro_state(n_dom), stat=ist)
    IF (ist /= success) THEN
      CALL finish(TRIM(routine),'allocation for p_hydro_state failed')
    ENDIF

    !------------------------------------------------------------------
    !  Divide patches and interpolation states for parallel runs.
    !  This is only done if the model runs really in parallel.
    !------------------------------------------------------------------
    IF(p_nprocs == 1 .OR. p_pe == p_test_pe .OR. lrestore_states) THEN
      
      ! This is a verification run or a run on a single processor
      ! or the divided states have been read, just set pointers
      
      p_patch => p_patch_global
      ! #slo# - temporarily switched on for comparison with rbf-reconstruction
      !IF (iequations /= ihs_ocean) THEN
      p_int_state => p_int_state_global
      p_grf_state => p_grf_state_global
      !ENDIF
      
      IF(p_nprocs == 1 .OR. p_pe == p_test_pe) THEN
        p_patch(:)%comm = p_comm_work
      ELSE
        CALL set_patch_communicators(p_patch)
      ENDIF
      
    ELSE
      
      CALL decompose_atmo_domain()
      
    ENDIF
    
    ! In case of a test run: Copy processor splitting to test PE
    IF(p_test_run) CALL copy_processor_splitting(p_patch)
    
    IF(ldump_states)THEN
      
      ! Dump divided patches with interpolation and grf state to NetCDF file and exit
      
      CALL message(TRIM(routine),'ldump_states is set: dumping patches+states and finishing')
      
      IF(p_pe /= p_test_pe) THEN
        DO jg = n_dom_start, n_dom
          CALL dump_patch_state_netcdf(p_patch(jg),p_int_state(jg),p_grf_state(jg))
        ENDDO
      ENDIF
      
      CALL p_stop
      STOP
      
    ENDIF
    
    !------------------------------------------------------------------
    ! step i: Setup physics
    !------------------------------------------------------------------
    IF ( lforcing ) THEN
      CALL setup_physics
    END IF
    
    !------------------------------------------------------------------
    ! step 6: initialize output
    !------------------------------------------------------------------
    
    CALL io_nml_setup
    CALL setup_gmt_output(p_patch(n_dom)%nlev)
    
    ! The model produces output files for all grid levels
    
    ALLOCATE(lprepare_output(n_dom),stat=ist)
    IF (ist /= success) THEN
      CALL finish(TRIM(routine),'allocation for lprepare_output failed')
    ENDIF
    
    !------------------------------------------------------------------
    ! Create and optionally read external data fields
    !------------------------------------------------------------------
    ALLOCATE (ext_data(n_dom), stat=ist)
    IF (ist /= success) THEN
      CALL finish(TRIM(routine),'allocation for ext_data failed')
    ENDIF
    
    !------------------------------------------------------------------
    ! Prepare raw data output file
    !------------------------------------------------------------------
    ! The index of the output file starts from 1.
    jfile   = 1
    
    ! Set up global attributes and variable lists for NetCDF output file

    CALL init_output_files(jfile, lclose=.FALSE.)
    
    ! Prepare time integration

    CALL prepare_ho_integration( p_patch(1:), p_hyoce_state,           &
                               & p_sfc_flx, p_phys_param, p_physics_oce)
    
    !------------------------------------------------------------------
    ! Daniel: Suggestion for point 5 of Feature #333
    ! (5. Subroutine to setup model components depending on this
    ! namelist, to be called after all namelists have been read, and
    ! a synoptic check has been done)
    !------------------------------------------------------------------
    ! set dependent variables/model components, depending on this (transport)
    ! namelist and potentially others
    IF (ltransport) THEN
      CALL setup_transport
    ENDIF
    
    !------------------------------------------------------------------
    ! Close the namelist files.
    !------------------------------------------------------------------
    CALL close_nml
    IF (p_pe == p_io) THEN
      CALL close_nml_output
    END IF
    
    !------------------------------------------------------------------
    ! Write out the inital values
    !------------------------------------------------------------------
    CALL write_vlist_oce( p_patch(1:), p_hyoce_state, p_sfc_flx, datetime )
    
    !------------------------------------------------------------------
    ! Now start the time stepping:
    ! The special initial time step for the three time level schemes
    ! is executed within process_grid_level
    !------------------------------------------------------------------
                                         ! Number of time steps to: 
    n_io   = NINT(dt_data/dtime)         ! - write output
    n_file = NINT(dt_file/dtime)         ! - trigger output file
    n_diag = MAX(1,NINT(dt_diag/dtime))  ! - compute global intervals
    
      CALL perform_ho_stepping( p_patch(1:), p_hyoce_state,&
        & datetime, n_io, n_file,    &
        & p_int_state(1:),           &
        & p_sfc_flx,                 &
        & p_phys_param)
    
    IF (ltimer) CALL print_timer
    
    !------------------------------------------------------------------
    !  cleaning up process
    !------------------------------------------------------------------
    CALL message(TRIM(routine),'start to clean up')
    
    CALL finalise_ho_integration(p_hyoce_state, p_phys_param, p_physics_oce)
    
    
    ! Delete variable lists

    DO jg = 1, n_dom
      CALL destruct_vlist_oce( jg )
    END DO
    
    ! Deallocate interpolation fields
    ! interpolation state not used for ocean model
    ! #slo# - temporarily switched on for comparison with rbf-reconstruction
    !IF (iequations /= ihs_ocean) THEN

    CALL destruct_2d_interpol_state( p_int_state )
    IF(p_nprocs == 1 .OR. p_pe == p_test_pe .OR. lrestore_states) THEN
      DEALLOCATE (p_int_state_global, stat=ist)
    ELSE
      DEALLOCATE (p_int_state_subdiv, stat=ist)
    ENDIF
    IF (ist /= success) THEN
      CALL finish(TRIM(routine),'deallocation for ptr_int_state failed')
    ENDIF
    
    !ENDIF
    
    ! Deallocate grid patches

    CALL destruct_patches( p_patch )

    IF(p_nprocs == 1 .OR. p_pe == p_test_pe .OR. lrestore_states) THEN
      DEALLOCATE( p_patch_global, stat=ist )
    ELSE
      DEALLOCATE( p_patch_subdiv, stat=ist )
    ENDIF
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocate for patch array failed')
    ENDIF
    
    ! deallocate output switches

    DEALLOCATE(lprepare_output,stat=ist)
    IF (ist /= success) THEN
      CALL finish(TRIM(routine),'deallocation of lprepare_output failed')
    ENDIF
    
    ! deallocate time level variables

    CALL cleanup_dyn_params  ! Time level variables 
    
    CALL message(TRIM(routine),'clean-up finished')
    
  END SUBROUTINE ocean_model
  
END MODULE mo_ocean_model

