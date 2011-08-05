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
  USE mo_master_control,      ONLY: is_restart_run, get_my_couple_id
  USE mo_parallel_config,     ONLY: p_test_run, l_test_openmp, num_io_procs
  USE mo_mpi,                 ONLY: p_stop, &
    & my_process_is_io,  my_process_is_mpi_seq, my_process_is_mpi_test, &
    & set_mpi_work_communicators, set_comm_input_bcast, null_comm_type
  USE mo_timer,               ONLY: init_timer, print_timer
  USE mo_datetime,            ONLY: t_datetime
  USE mo_output,              ONLY: init_output_files
  USE mo_io_vlist,            ONLY: write_vlist_oce, destruct_vlist_oce
  USE mo_grid_config,         ONLY: n_dom, n_dom_start, global_cell_type !, &
!                                  dynamics_parent_grid_id
  USE mo_dynamics_config,     ONLY: iequations
  
  USE mo_io_async,            ONLY: io_main_proc            ! main procedure for I/O PEs
  
  USE mo_interpol_config,     ONLY: configure_interpolation
  USE mo_advection_config,    ONLY: configure_advection
  USE mo_dynamics_config,     ONLY: configure_dynamics  ! subroutine
  
  ! Control parameters: run control, dynamics, i/o
  !
  USE mo_io_config,           ONLY:  dt_data,dt_file,dt_diag!,dt_checkpoint
  USE mo_run_config,          ONLY: &
    & dtime,                & !    :
!    & ltransport,           & !    :
    & ltimer,               & !    :
    & ldump_states,         & ! flag if states should be dumped
    & lrestore_states,      &  ! flag if states should be restored
    & iforcing,             & !  
    & nlev, nlevp1,         & !
    & num_lev, num_levp1,   &
    & iqv, nshift,          &
    & ntracer
  
!  USE mo_advection_nml,       ONLY: transport_nml_setup,  & ! process transport
!    & setup_transport         ! control parameters
  
  USE mo_subdivision,         ONLY: decompose_atmo_domain,         &
    & copy_processor_splitting,      &
    & set_patch_communicators
  USE mo_dump_restore,        ONLY: dump_patch_state_netcdf,       &
    & restore_patches_netcdf,        &
    & restore_interpol_state_netcdf
  
  USE mo_icoham_dyn_memory,   ONLY: p_hydro_state
  USE mo_atmo_control,        ONLY: p_patch_global, p_patch_subdiv, p_patch
  USE mo_intp_data_strc,      ONLY: p_int_state_global, p_int_state_subdiv, p_int_state
  USE mo_grf_intp_data_strc,  ONLY: p_grf_state_global, p_grf_state
  
  ! Horizontal grid
  !
  USE mo_model_domain_import, ONLY: &!  grid_nml_setup,          & ! process grid control parameters
    & n_dom,                & !    :
    & n_dom_start,          & !    :
    & import_patches,       & !
    & destruct_patches        !
  
  ! Horizontal interpolation
  !
!  USE mo_interpol_nml,        ONLY: interpol_nml_setup   ! process interpol. ctl. params.
  USE mo_intp_state,          ONLY: construct_2d_interpol_state, &
    & destruct_2d_interpol_state
  
  USE mo_oce_state,           ONLY: t_hydro_ocean_state 
  
!0  USE mo_mpiom_phy_state,     ONLY: construct_mpiom_phy_state, &
!0    &                               destruct_mpiom_phy_state
  
  USE mo_impl_constants,      ONLY: success !, ihs_ocean
  
  ! External data
  USE mo_ext_data,            ONLY: ext_data, init_ext_data, destruct_ext_data
  
!  USE mo_test_hydro_ocean,    ONLY: prepare_ho_integration, finalise_ho_integration!, &
!    &                               test_hydro_ocean
  USE mo_hydro_ocean_run,      ONLY: perform_ho_stepping,&
    & prepare_ho_integration,&
    & finalise_ho_integration
  
  USE mo_oce_forcing,         ONLY: t_ho_sfc_flx
  USE mo_oce_physics,         ONLY: t_ho_params, t_ho_physics

  ! For the coupling
  USE mo_impl_constants,      ONLY: CELLS, MAX_CHAR_LENGTH
  USE mo_master_control,      ONLY : ocean_process, is_coupled_run
  USE mo_icon_cpl_init_comp,  ONLY : get_my_local_comp_id
  USE mo_icon_cpl_def_grid,   ONLY : ICON_cpl_def_grid
  USE mo_icon_cpl_def_field,  ONLY : ICON_cpl_def_field
  USE mo_icon_cpl_search,     ONLY : ICON_cpl_search
  USE mo_model_domain_import, ONLY : get_patch_global_indexes

  !-------------------------------------------------------------
  USE mo_read_namelists,      ONLY: read_ocean_namelists

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

    ! Ad-hoc declaration for the hydrostatic ocean state variables.
    ! This declaration and the declarations for the atmospheric model states, see
    ! "USE mo_atmo_control" above, should eventually be organized in the same way.
    TYPE(t_hydro_ocean_state), TARGET, ALLOCATABLE, SAVE :: p_hyoce_state(:)
    TYPE(t_ho_sfc_flx)                              :: p_sfc_flx
    TYPE (t_ho_params)                              :: p_phys_param
    TYPE(t_ho_physics)                              :: p_physics_oce
    TYPE(t_datetime)                                :: datetime
    
    INTEGER :: n_io, jg, jfile, n_file, ist, n_diag

    ! For the coupling
    
    INTEGER, PARAMETER :: no_of_fields = 12

    CHARACTER(LEN=MAX_CHAR_LENGTH) ::  field_name(no_of_fields)
    INTEGER :: field_id(no_of_fields)
    INTEGER :: comp_id
    INTEGER :: grid_id
    INTEGER :: grid_shape(2) 
    INTEGER :: field_shape(3) 
    INTEGER :: i, error_status
    INTEGER :: no_of_entities
    INTEGER :: patch_no
    INTEGER, POINTER :: grid_glob_index(:)

    !-------------------------------------------------------------------
    
    !---------------------------------------------------------------------
    ! 1.1 Read namelists (newly) specified by the user; fill the 
    !     corresponding sections of the configuration states.
    !---------------------------------------------------------------------

    CALL read_ocean_namelists(oce_namelist_filename,shr_namelist_filename)

    !---------------------------------------------------------------------
    ! 1.2 Cross-check namelist setups
    !---------------------------------------------------------------------

!     CALL atm_crosscheck

    !---------------------------------------------------------------------
    ! 2. Call configure_run to finish filling the run_config state.
    !    This needs to be done very early (but anyway after atm_crosscheck)
    !    because some component of the state, e.g., num_lev, may be 
    !    modified in this subroutine which affect the following CALLs.
    !---------------------------------------------------------------------
!     CALL configure_run

    !-------------------------------------------------------------------
    ! 3.1 Initialize the mpi work groups
    !-------------------------------------------------------------------
    CALL set_mpi_work_communicators(p_test_run, l_test_openmp, num_io_procs)

    !-------------------------------------------------------------------
    ! 3.2 Initialize various timers
    !-------------------------------------------------------------------
    IF (ltimer) CALL init_timer

    
    !-------------------------------------------------------------------
    ! Initialize date and time
    !-------------------------------------------------------------------
    datetime = time_config%ini_datetime
    
    !-------------------------------------------------------------------
    ! 4. Import patches
    !-------------------------------------------------------------------
    ! If we belong to the I/O PEs just call io_main_proc before reading patches.
    ! This routine will never return
    
    IF (my_process_is_io()) CALL io_main_proc
    
    ! Check patch allocation status

    IF ( ALLOCATED(p_patch_global)) THEN
      CALL finish(TRIM(routine), 'patch already allocated')
    END IF
     
    ! Allocate patch array to start patch construction

    ALLOCATE(p_patch_global(n_dom_start:n_dom), stat=error_status)
    IF (error_status/=success) THEN
      CALL finish(TRIM(routine), 'allocation of patch failed')
    ENDIF
    
    IF(lrestore_states) THEN
      ! Before the restore set p_comm_input_bcast to null
      CALL set_comm_input_bcast(null_comm_type)
      IF( .NOT. my_process_is_mpi_test()) THEN
        CALL restore_patches_netcdf( p_patch_global )
      ELSE
        CALL import_patches( p_patch_global,                       &
                            nlev,nlevp1,num_lev,num_levp1,nshift )
      ENDIF
      ! After the restore is done set p_comm_input_bcast in the
      ! same way as it would be set in parallel_nml_setup when
      ! no restore is wanted:
      CALL set_comm_input_bcast()
    ELSE
        CALL import_patches( p_patch_global,                       &
                            nlev,nlevp1,num_lev,num_levp1,nshift )      
    ENDIF

    !--------------------------------------------------------------------------------
    ! 5. Construct interpolation state, compute interpolation coefficients.
    !--------------------------------------------------------------------------------
    
    CALL configure_interpolation( global_cell_type, n_dom, p_patch_global(1:)%level )

    ! Allocate array for interpolation state
    
    ALLOCATE( p_int_state_global(n_dom_start:n_dom), &
            & p_grf_state_global(n_dom_start:n_dom),STAT=error_status)
    IF (error_status /= SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for ptr_int_state failed')
    ENDIF
    
    IF(lrestore_states .AND. .NOT. my_process_is_mpi_test()) THEN
      ! Read interpolation state from NetCDF
      CALL restore_interpol_state_netcdf(p_patch_global, p_int_state_global)
    ELSE
      ! Interpolation state is constructed for
      ! the full domain on every PE and divided later
      CALL construct_2d_interpol_state(p_patch_global, p_int_state_global)
    ENDIF
    
    
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

    
    !-------------------------------------------------------------------
    ! 7. Domain decomposition: 
    !    Divide patches and interpolation states for parallel runs.
    !    This is only done if the model runs really in parallel.
    !-------------------------------------------------------------------   
    IF (my_process_is_mpi_seq()  &
      &  .OR. lrestore_states) THEN
      
      ! This is a verification run or a run on a single processor
      ! or the divided states have been read, just set pointers
      
      p_patch => p_patch_global
      p_int_state => p_int_state_global
      p_grf_state => p_grf_state_global
      
!       IF (my_process_is_mpi_seq()) THEN
!         p_patch(:)%comm = p_comm_work
!       ELSE
        CALL set_patch_communicators(p_patch)
!       ENDIF
      
    ELSE
      
      CALL decompose_atmo_domain( )
      
    ENDIF

    ! Note: fro this point the p_patch is used
    ! In case of a test run: Copy processor splitting to test PE
    IF(p_test_run) CALL copy_processor_splitting(p_patch)
    
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
      CALL configure_advection( jg, p_patch(jg)%nlev, p_patch(1)%nlev,      &
        &                      iequations, iforcing, iqv, 0, &
        &                      0, .false., .true., ntracer )
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
      CALL init_ext_data (p_patch(1:), p_int_state, ext_data)
    

    !------------------------------------------------------------------
    ! Prepare raw data output file
    !------------------------------------------------------------------
    ! The index of the output file starts from 1.
    jfile   = 1
    
    ! Set up global attributes and variable lists for NetCDF output file
     n_io    = NINT(dt_data/dtime)        ! write output
     n_file  = NINT(dt_file/dtime)        ! trigger new output file
!      n_chkpt = NINT(dt_checkpoint/dtime)  ! write restart files
     n_diag  = MAX(1,NINT(dt_diag/dtime)) ! diagnose of total integrals

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
!    IF (ltransport) THEN
!      CALL setup_transport( ihs_ocean )
!    ENDIF
    
    !------------------------------------------------------------------
    ! Write out the inital values
    !------------------------------------------------------------------
    CALL write_vlist_oce( p_patch(1:), p_hyoce_state, p_sfc_flx, datetime )


    !------------------------------------------------------------------
    ! Prepare the coupling
    !
    ! For the time being this could all go into a subroutine which is
    ! common to atmo and ocean. Does this make sense if the setup deviates
    ! too much in future.
    !------------------------------------------------------------------

    IF ( is_coupled_run() ) THEN
 
 
      comp_id = get_my_couple_id ()
      patch_no = 1

      grid_shape(1)=1
      grid_shape(2)=p_patch(patch_no)%n_patch_cells

      ! CALL get_patch_global_indexes ( patch_no, CELLS, no_of_entities, grid_glob_index )
      ! should grid_glob_index become a pointer in ICON_cpl_def_grid as well?
      CALL ICON_cpl_def_grid ( &
        & comp_id, grid_shape, p_patch(patch_no)%cells%glb_index, & ! input
        & grid_id, error_status )                                   ! output
  
      field_name(1) = "SST"
      field_name(2) = "TAUX"
      field_name(3) = "TAUY"
      field_name(4) = ""
      field_name(5) = ""
      field_name(6) = ""
      field_name(7) = ""
      field_name(8) = ""
 
!       DO i = 1, no_of_fields
!          CALL ICON_cpl_def_field ( field_name(i), comp_id, grid_id, field_id(i), &
!    &                               field_shape, error_status )
!       ENDDO

      CALL ICON_cpl_search

    ENDIF

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
    
    ! Delete variable lists

    DO jg = 1, n_dom
      CALL destruct_vlist_oce( jg )
    END DO
    
    ! Deallocate interpolation fields
    ! interpolation state not used for ocean model
    ! #slo# - temporarily switched on for comparison with rbf-reconstruction
    !IF (iequations /= ihs_ocean) THEN

    CALL destruct_2d_interpol_state( p_int_state )
    IF(my_process_is_mpi_seq() .OR. lrestore_states) THEN
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

    IF(my_process_is_mpi_seq() .OR. lrestore_states) THEN
      DEALLOCATE( p_patch_global, stat=ist )
    ELSE
      DEALLOCATE( p_patch_subdiv, stat=ist )
    ENDIF
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocate for patch array failed')
    ENDIF
    
    CALL message(TRIM(routine),'clean-up finished')

  END SUBROUTINE ocean_model
  
END MODULE mo_ocean_model

