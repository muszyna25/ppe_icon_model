!>
!! @brief Main program for the ICON atmospheric model
!!
!! @author
!!  Leonidas Linardakis (MPI-M)
!!  Hui Wan             (MPI-M)
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
MODULE mo_atmo_model

USE mo_exception,           ONLY: message, finish
USE mo_mpi,                 ONLY: p_stop, p_pe, p_io, p_nprocs, &
& p_comm_work_test, p_comm_input_bcast, p_comm_work, &
& my_process_is_io,  my_process_is_mpi_seq, my_process_is_mpi_test, &
& my_process_is_stdio
USE mo_timer,               ONLY: init_timer, print_timer
USE mo_master_nml,          ONLY: lrestart
USE mo_namelist,            ONLY: open_nml,  close_nml, open_nml_output, close_nml_output
USE mo_output,              ONLY: init_output_files, close_output_files, write_output

USE mo_parallel_configuration, ONLY: p_test_run

USE mo_io_async,            ONLY: io_main_proc            ! main procedure for I/O PEs

! Control parameters: run control, dynamics, i/o
!
USE mo_global_variables,    ONLY: setup_physics           ! process forcing control parameters
USE mo_nonhydrostatic_nml,  ONLY: ivctype,              & ! type of vertical coordinate
& nonhydrostatic_nml_setup
USE mo_dynamics_nml,        ONLY: dynamics_nml_setup
USE mo_diffusion_nml,       ONLY: diffusion_nml_setup
USE mo_io_nml,              ONLY: io_nml_setup !,         & ! process I/O
!& dt_data,              & !    :
!& dt_file,              & !    :
!& dt_diag,              & !    :
!& dt_checkpoint 
USE mo_io_config,         ONLY:  dt_data,dt_file,dt_diag,dt_checkpoint
USE mo_dynamics_config,   ONLY: iequations
USE mo_run_nml,           ONLY: run_nml_setup
USE mo_run_config,        ONLY: &
& dtime,                & !    namelist parameter
& nsteps,               & !    :
& ltransport,           & !    :
& lforcing,             & !    :
& ltestcase,            & !    :
& ltimer,               & !    :
& iforcing,             & !    namelist parameter
& ldump_states,         & ! flag if states should be dumped
& lrestore_states,      & ! flag if states should be restored
& nlev,nlevp1,          &
& num_lev,num_levp1,nshift

USE mo_impl_constants, ONLY:&
    & ihs_atm_temp,         & !    :
    & ihs_atm_theta,        & !    :
    & inh_atmosphere,       & !    :
    & ishallow_water,       & !    :
    & ildf_dry,             & !    :
    & ildf_echam,           & !    :
    & inoforcing,           & !    :
    & iheldsuarez,          & !    :
    & iecham,               & !    :
    & inwp

USE mo_advection_nml,       ONLY: transport_nml_setup,  & ! process transport
& setup_transport         ! control parameters

! For the coupling
USE mo_impl_constants, ONLY: CELLS
USE mo_master_control, ONLY : atmo_process, is_coupled_run
USE mo_icon_cpl_init_comp, ONLY : get_my_local_comp_id
USE mo_icon_cpl_def_grid, ONLY : ICON_cpl_def_grid
USE mo_icon_cpl_def_field, ONLY : ICON_cpl_def_field
USE mo_icon_cpl_search, ONLY : ICON_cpl_search
USE mo_model_domain_import, ONLY : get_patch_global_indexes

! Test cases
!
USE mo_hydro_testcases,     ONLY: setup_testcase          ! process hyd. atm. tests ctl. params.
USE mo_nh_testcases,        ONLY: read_nh_testcase_namelist ! process non-hyd. atm. test ctl. par.

! Memory
!
USE mo_subdivision,         ONLY: decompose_atmo_domain,         &
& copy_processor_splitting,      &
& set_patch_communicators
USE mo_dump_restore,        ONLY: dump_patch_state_netcdf,       &
& restore_patches_netcdf,        &
& restore_interpol_state_netcdf, &
& restore_gridref_state_netcdf

USE mo_icoham_dyn_memory,   ONLY: p_hydro_state
USE mo_nonhydro_state,      ONLY: p_nh_state
USE mo_atmo_control,        ONLY: p_patch_global, p_patch_subdiv, p_patch

! Horizontal grid
!
USE mo_grid_nml,            ONLY: read_grid_namelist
USE mo_model_domain_import, ONLY: &  !grid_nml_setup,          & ! process grid control parameters
& n_dom,                & !    :
& n_dom_start,          & !    :
!    & parent_id,            & !    :
& import_patches,       & !
& destruct_patches        !

USE mo_grid_configuration,   ONLY: parent_id

! Horizontal interpolation
!
USE mo_interpol_nml,        ONLY: interpol_nml_setup   ! process interpol. ctl. params.
USE mo_intp_state,          ONLY: construct_2d_interpol_state, &
& destruct_2d_interpol_state
USE mo_interpolation,       ONLY: rbf_vec_interpol_cell,       &
& edges2cells_scalar
USE mo_gridref_nml,         ONLY: gridref_nml_setup
USE mo_grf_interpolation,   ONLY: construct_2d_gridref_state,  &
& destruct_2d_gridref_state

! Vertical grid
!
USE mo_vertical_coord_table,ONLY: init_vertical_coord_table
USE mo_vertical_grid,       ONLY: init_hybrid_coord, init_sleve_coord

! State variables
!
USE mo_icoham_dyn_memory,   ONLY: destruct_icoham_dyn_state
USE mo_nonhydro_state,      ONLY: destruct_nh_state


! Parameterized forcing
!
USE mo_echam_phy_memory,    ONLY: destruct_echam_phy_state
USE mo_echam_phy_setup,     ONLY: setup_echam_phy
USE mo_echam_phy_init,      ONLY: prepare_echam_phy, initcond_echam_phy, &
			  & additional_restart_init
USE mo_echam_phy_cleanup,   ONLY: cleanup_echam_phy
USE mo_gmt_output,          ONLY: setup_gmt_output
USE mo_nwp_phy_state,       ONLY: construct_nwp_phy_state,   &
& destruct_nwp_phy_state
USE mo_atm_phy_nwp_config, ONLY: atm_phy_nwp_config,setup_atm_nwp_phy
USE mo_lnd_nwp_nml,         ONLY: setup_nwp_lnd
USE mo_nwp_lnd_state,       ONLY: construct_nwp_lnd_state,   &
& destruct_nwp_lnd_state, p_lnd_state

USE mo_impl_constants,      ONLY: SUCCESS, MAX_CHAR_LENGTH

! Time integration
!
USE mo_ha_stepping,         ONLY: prepare_ha_dyn, initcond_ha_dyn, perform_ha_stepping
USE mo_nh_stepping,         ONLY: prepare_nh_integration, perform_nh_stepping
! External data
USE mo_ext_data,            ONLY: ext_data, init_ext_data, destruct_ext_data

!  USE mo_nwp_phy_init,          ONLY: init_nwp_phy
!!$  USE mo_gscp_cosmo,          ONLY: hydci_pp_init


!-------------------------------------------------------------------------
! to break circular dependency

USE mo_intp_data_strc,      ONLY: p_int_state_global, p_int_state_subdiv, p_int_state
USE mo_grf_intp_data_strc,  ONLY: p_grf_state_global, p_grf_state_subdiv, p_grf_state

!-------------------------------------------------------------------------
USE mo_io_restart,           ONLY: read_restart_info_file, read_restart_files
USE mo_io_restart_namelist,  ONLY: read_restart_namelists
USE mo_io_restart_attributes,ONLY: read_restart_attributes, get_restart_attribute

USE mo_atmo_setup_configuration, ONLY: read_atmo_namelists
USE mo_atm_nml_crosscheck,       ONLY: atm_crosscheck

USE mo_time_config,     ONLY: time_config      ! variable
USE mo_dynamics_config, ONLY: config_dynamics  ! subroutine


!-------------------------------------------------------------------------
IMPLICIT NONE
PRIVATE

PUBLIC :: atmo_model

CONTAINS
!>
!!
SUBROUTINE atmo_model(namelist_filename)

CHARACTER(LEN=*), INTENT(in) :: namelist_filename

CHARACTER(LEN=MAX_CHAR_LENGTH) :: grid_file_name 
CHARACTER(*), PARAMETER :: routine = "mo_atmo_model:atmo_model"
LOGICAL :: lsuccess

! For the coupling

INTEGER, PARAMETER :: no_of_fields = 12

CHARACTER(LEN=MAX_CHAR_LENGTH) ::  field_name(no_of_fields)
INTEGER :: field_id(no_of_fields)
INTEGER :: comp_id
INTEGER :: grid_id
INTEGER :: grid_shape(2) 
INTEGER :: field_shape(3) 
INTEGER :: i, ierror
INTEGER :: no_of_entities
INTEGER :: patch_no
INTEGER, POINTER :: grid_glob_index(:)

    !---------------------------------------------------------------------
    ! 0. If this is a resumed or warm-start run...
    !---------------------------------------------------------------------
    IF (lrestart) THEN

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

      CALL read_restart_namelists('restart_atm.nc')
      CALL message(TRIM(routine), 'read namelists from atm restart file')

      ! Read all global attributs in the restart file and store them in a buffer.

      CALL read_restart_attributes('restart_atm.nc')
      CALL message(TRIM(routine), 'read global attributes from atm restart file')

    END IF ! lrestart

    !---------------------------------------------------------------------
    ! 1. Read namelists (newly) specified by the user; fill the 
    !    corresponding sections of the configuration states.
    !---------------------------------------------------------------------
    CALL read_atmo_namelists(namelist_filename)

    !---------------------------------------------------------------------
    ! 2. Cross-check namelists
    !---------------------------------------------------------------------

    CALL atm_crosscheck

    !---------------------------------------------------------------------
    ! 3. Assign values to derived variables in the configuration states
    !---------------------------------------------------------------------
    CALL config_dynamics( lrestart, n_dom )

    !---------------------------------------------------------------------
    ! 4. Construct model states (variable lists); set initial conditions.
    !---------------------------------------------------------------------

    !---------------------------------------------------------------------
    ! 4.a Do the setup for the coupled run
    !
    ! For the time being this could all go into a subroutine which is
    ! common to atmo and ocean. Does this make sense if the setup deviates
    ! too much in future.
    !---------------------------------------------------------------------

    IF ( is_coupled_run() ) THEN
 
      comp_id = get_my_local_comp_id (atmo_process)
      patch_no = 1

      CALL get_patch_global_indexes ( patch_no, CELLS, no_of_entities, grid_glob_index )
      ! should grid_glob_index become a pointer in ICON_cpl_def_grid as well?
      CALL ICON_cpl_def_grid ( comp_id, grid_shape, grid_glob_index, grid_id, ierror )
  
      field_name(1) = "SST"
      field_name(2) = "TAUX"
      field_name(3) = "TAUY"
      field_name(4) = ""
      field_name(5) = ""
      field_name(6) = ""
      field_name(7) = ""
      field_name(8) = ""

      DO i = 1, no_of_fields
         CALL ICON_cpl_def_field ( field_name(i), comp_id, grid_id, field_id(i), &
   &                               field_shape, ierror )
      ENDDO

      CALL ICON_cpl_search

    ENDIF
    !
    !---------------------------------------------------------------------
    ! 5. Perform time stepping
    !---------------------------------------------------------------------

    !---------------------------------------------------------------------
    ! 6. Integration finished. Clean up.
    !---------------------------------------------------------------------

  END SUBROUTINE atmo_model




  !---------------------------------------------------------------------
  !---------------------------------------------------------------------
  !---------------------------------------------------------------------
  !---------------------------------------------------------------------
  !---------------------------------------------------------------------
  SUBROUTINE atmo_model_old(namelist_filename)
    
    CHARACTER(LEN=*), INTENT(in) :: namelist_filename

    CHARACTER(*), PARAMETER :: routine = "mo_atmo_model:atmo_model"

    CHARACTER(LEN=MAX_CHAR_LENGTH) :: grid_file_name 
    INTEGER :: n_io, jg, jfile, n_file, ist, n_diag, n_chkpt
    LOGICAL :: l_have_output
   
    !-------------------------------------------------------------------
    ! 1. Open the atmosphere-specific namelist file and create a 
    !    new file in which all the namelist variables and their
    !    actual values used in the model run will be stored.
    !-------------------------------------------------------------------
    CALL open_nml(TRIM(namelist_filename))
    IF(my_process_is_stdio()) CALL open_nml_output('NAMELIST_ICON_output_atm')

    ! The namelists ('run_nml' and 'testcase_ctl') are read in seperate
    ! subroutines. The validity of the user-specified configurations is
    ! checked right after each namelist is read.
    ! The two 'setup' subroutines above must be called before grid and patch
    ! import because during the import procedure the topography and land-sea
    ! mask will be initialized. They are related to the atmos/ocean switch
    ! as well as the selected test case.
    
    CALL run_nml_setup
    
    !-------------------------------------------------------------------
    ! Initialize various timers; initialize data and time 
    !-------------------------------------------------------------------
    IF (ltimer) CALL init_timer
    
    ! parallel_nml_setup must be called after setup_run since it needs
    ! some variables read in setup_run
    
!     CALL parallel_nml_setup
    
    !-------------------------------------------------------------------
    ! Initialize test case setup 
    !-------------------------------------------------------------------
    IF (ltestcase) THEN
      
      SELECT CASE (iequations)
      CASE (ishallow_water, ihs_atm_temp, ihs_atm_theta)
        CALL setup_testcase
        
      CASE (inh_atmosphere)
        CALL read_nh_testcase_namelist
        
      CASE DEFAULT
        CALL finish( TRIM(routine),' invalid value for iequations!' )
      END SELECT
      
    ENDIF
    
    !-------------------------------------------------------------------
    ! step 2: import computational domain of the model
    !-------------------------------------------------------------------
    ! 2a) read namelist 'grid_ctl' from the namelist file (already
    !     opened above) and set up the grid/patch configuration.
    !-------------------------------------------------------------------
    
!    CALL grid_nml_setup
!   CALL read_grid_namelist( )
    
    !-------------------------------------------------------------------
    ! 2b) patch import
    !-------------------------------------------------------------------
    
    !-------------------------------------------------------------------
    ! If we belong to the I/O PEs just call io_main_proc before reading patches.
    ! This routine will never return
    
    IF (my_process_is_io()) CALL io_main_proc
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
    
    IF(lrestore_states .AND. .NOT. my_process_is_mpi_test()) THEN
      CALL restore_patches_netcdf( p_patch_global )
    ELSE
      CALL import_patches( p_patch_global,nlev,nlevp1,num_lev,num_levp1,nshift )
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
    ! Read specific dynamics namelists for the nonhydrost. dynamical core
    !------------------------------------------------------------------
    SELECT CASE (iequations)
    
    CASE (inh_atmosphere)
      CALL nonhydrostatic_nml_setup
      
    CASE DEFAULT
    END SELECT
    
    !------------------------------------------------------------------
    ! step 3a-a: ! read nwp physics namelist, ...
    !------------------------------------------------------------------
    IF ( iforcing == inwp) THEN
      CALL setup_atm_nwp_phy !( p_patch_global(1:) )  ! read Namelist, ...
!      IF (inwp_surface > 0)
      CALL setup_nwp_lnd
    ENDIF
    
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
    ! KF empty SR, call needs to be shifted
    CALL gridref_nml_setup
    
    ! interpolation state not used for ocean model
    ! #slo# - temporarily switched on for comparison with rbf-reconstruction
    
    CALL interpol_nml_setup(p_patch_global)
    
    !
    ! allocate type for interpolation state
    !
    
    ALLOCATE( p_int_state_global(n_dom_start:n_dom), &
            & p_grf_state_global(n_dom_start:n_dom),STAT=ist)
    IF (ist /= SUCCESS) THEN
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
    ! step 5a: init the structure of the model equations
    !------------------------------------------------------------------
    SELECT CASE (iequations)
    
    CASE (ishallow_water)
      CALL init_vertical_coord_table(iequations, p_patch_global(1)%nlev)
      
    CASE (ihs_atm_temp, ihs_atm_theta)
      CALL init_vertical_coord_table(iequations, p_patch_global(1)%nlev)
      
    CASE (inh_atmosphere)
      IF (ivctype == 1) THEN
        CALL init_hybrid_coord(iequations, p_patch_global(1)%nlev)
      ELSE IF (ivctype == 2) THEN
        CALL init_sleve_coord(p_patch_global(1)%nlev)
      ENDIF
    CASE DEFAULT
    END SELECT
 
    ! For the NH model, the initialization routines called from
    ! construct_2d_gridref_state require the metric terms to be present
    IF (n_dom_start==0 .OR. n_dom > 1) THEN
      IF(lrestore_states .AND. .NOT. my_process_is_mpi_test()) THEN
        ! Read gridref state from NetCDF
        CALL restore_gridref_state_netcdf(p_patch_global, p_grf_state_global)
      ELSE
        CALL construct_2d_gridref_state (p_patch_global, p_grf_state_global)
      ENDIF
    ENDIF
    
    !-------------------------------------------------------------------------
    ! set up horizontal diffusion after the vertical coordinate is configured
    !-------------------------------------------------------------------------
    CALL diffusion_nml_setup(n_dom,parent_id,nlev)
    
    !------------------------------------------------------------------
    ! step 5b: allocate state variables
    !------------------------------------------------------------------
    
    SELECT CASE (iequations)
    !
    CASE (ishallow_water, ihs_atm_temp, ihs_atm_theta)
      ALLOCATE (p_hydro_state(n_dom), stat=ist)
      IF (ist /= success) THEN
        CALL finish(TRIM(routine),'allocation for p_hydro_state failed')
      ENDIF
      !
    CASE (inh_atmosphere)
!     ALLOCATE (p_nh_state(n_dom), stat=ist)
!     IF (ist /= success) THEN
!       CALL finish(TRIM(routine),'allocation for p_nh_state failed')
!     ENDIF
!     ALLOCATE (p_lnd_state(n_dom), stat=ist)
!     IF (ist /= success) THEN
!       CALL finish(TRIM(routine),'allocation for p_lnd_state failed')
!     ENDIF
      !
    END SELECT
    
    !------------------------------------------------------------------
    !  Divide patches and interpolation states for parallel runs.
    !  This is only done if the model runs really in parallel.
    !------------------------------------------------------------------
    
    IF (my_process_is_mpi_seq()  &
      &  .OR. lrestore_states) THEN
      
      ! This is a verification run or a run on a single processor
      ! or the divided states have been read, just set pointers
      
      p_patch => p_patch_global
      p_int_state => p_int_state_global
      p_grf_state => p_grf_state_global
      
      IF (my_process_is_mpi_seq()) THEN
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
      
      IF(.NOT. my_process_is_mpi_test()) THEN
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
    
    CALL io_nml_setup ! is empty and now only a place holder
    CALL setup_gmt_output(p_patch(n_dom)%nlev)
    
    ! The model produces output files for all grid levels
 
    
    !------------------------------------------------------------------
    ! Create and optionally read external data fields
    !------------------------------------------------------------------
    ALLOCATE (ext_data(n_dom), STAT=ist)
    IF (ist /= SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for ext_data failed')
    ENDIF
    
    ! allocate memory for atmospheric/oceanic external data and
    ! optionally read those data from netCDF file.
    CALL init_ext_data (p_patch(1:), p_int_state, ext_data)
    
    ! Parameterized forcing:
    ! 1. Create forcing state variables
    ! 2. Set up parameterizations
    !
    SELECT CASE (iforcing)
    CASE (inoforcing,iheldsuarez,ildf_dry)
      ! nothing to be done
      !
    CASE (iecham,ildf_echam)
      ! ECHAM forcing
      CALL setup_echam_phy
    CASE (inwp)
      ! NWP forcing
      !    CALL setup_nwp_phy( p_patch(1:) )  ! read Namelist, ... moved before setup transport
      CALL construct_nwp_phy_state( p_patch(1:) )
      CALL construct_nwp_lnd_state( p_patch(1:),p_lnd_state,n_timelevels=2 )
      !    CALL init_nwp_phy(p_patch(1:))
    CASE default
      CALL finish(TRIM(routine),'iforcing has value that is not allowed')
    END SELECT

    !------------------------------------------------------------------
    ! Prepare for time integration
    !------------------------------------------------------------------
    SELECT CASE (iequations)
    CASE (ishallow_water, ihs_atm_temp, ihs_atm_theta)

    !------------------------------------------------------------------
    ! Initialize parameters and solvers;
    ! Allocate memory for model state vectors.
    !------------------------------------------------------------------
    CALL prepare_ha_dyn( p_patch(1:) )
    IF (iforcing==IECHAM.OR.iforcing==ILDF_ECHAM) THEN
      CALL prepare_echam_phy( p_patch(1:) )
    END IF

    !------------------------------------------------------------------
    ! Set initial conditions for time integration.
    !------------------------------------------------------------------
    IF (lrestart) THEN
    ! This is an resumed integration. Read model state from restart file(s).

      CALL read_restart_files
      CALL message(TRIM(routine),'normal exit from read_restart_files')

      ! Initialize logical variables in echam physics state.
      ! This is necessary for now because logical arrays can not yet be
      ! written into restart files.

      IF (iforcing==IECHAM.OR.iforcing==ILDF_ECHAM) THEN                                       
        CALL additional_restart_init( p_patch(1:) )                                            
      END IF                                                                                   

    ELSE
    ! This is an initial run (cold start). Compute initial condition for 
    ! test cases, or read externally given initial conditions.

      CALL initcond_ha_dyn( p_patch(1:), p_int_state(1:),  &
                          & p_grf_state(1:), p_hydro_state )

      IF (iforcing==IECHAM.OR.iforcing==ILDF_ECHAM)      &
      CALL initcond_echam_phy( p_patch(1:),p_hydro_state )

    END IF ! lrestart
                                                                                               
    !--------------------
    CASE (inh_atmosphere)
      CALL prepare_nh_integration(p_patch(1:), p_nh_state, p_int_state(1:), p_grf_state(1:))
 
    CASE DEFAULT
    END SELECT
 
    !------------------------------------------------------------------
    ! Daniel: Suggestion for point 5 of Feature #333
    ! (5. Subroutine to setup model components depending on this
    ! namelist, to be called after all namelists have been read, and
    ! a synoptic check has been done)
    !------------------------------------------------------------------
    ! set dependent variables/model components, depending on this (transport)
    ! namelist and potentially others
    IF (ltransport) THEN
      CALL setup_transport( iequations )
    ENDIF
    
    !------------------------------------------------------------------
    ! Close the namelist files.
    !------------------------------------------------------------------
    
    CALL close_nml
    IF (my_process_is_stdio()) THEN
      CALL close_nml_output
    END IF
    
    !---------------------------------------------------------
    ! The most primitive event handling algorithm: 
    ! compute time step interval for taking a certain action
    !--------------------------------------------------------- 
 
    n_io    = NINT(dt_data/dtime)        ! write output
    n_file  = NINT(dt_file/dtime)        ! trigger new output file
    n_chkpt = NINT(dt_checkpoint/dtime)  ! write restart files
    n_diag  = MAX(1,NINT(dt_diag/dtime)) ! diagnose of total integrals

    !------------------------------------------------------------------
    ! Prepare output file
    !------------------------------------------------------------------
    IF (.NOT.lrestart) THEN
    ! Initialize the first output file which will contain also the 
    ! initial conditions.

      jfile = 1
      CALL init_output_files(jfile, lclose=.FALSE.)

    ELSE
    ! No need to write out the initial condition, thus no output
    ! during the first integration step. This run will produce
    ! output if n_io <= integration_length. 

      CALL get_restart_attribute('next_output_file',jfile)

      IF (n_io.le.(nsteps-1)) THEN
         CALL init_output_files(jfile, lclose=.FALSE.)
         l_have_output = .TRUE.  
      ELSE
         l_have_output = .FALSE.
      END IF

    END IF
 
    !------------------------------------------------------------------
    !  get and write out some of the inital values
    !------------------------------------------------------------------
    IF (.NOT.lrestart) THEN

    ! diagnose u and v to have meaningful initial output
    
    DO jg = 1, n_dom

      SELECT CASE (iequations)

      CASE (ishallow_water, ihs_atm_temp, ihs_atm_theta)
        SELECT CASE (p_patch(jg)%cell_type)
        CASE (3)
          CALL rbf_vec_interpol_cell(p_hydro_state(jg)%prog(1)%vn,p_patch(jg), &
            & p_int_state(jg),p_hydro_state(jg)%diag%u,p_hydro_state(jg)%diag%v)
        CASE (6)
          CALL edges2cells_scalar(p_hydro_state(jg)%prog(1)%vn,p_patch(jg), &
            & p_int_state(jg)%hex_east,p_hydro_state(jg)%diag%u)
          CALL edges2cells_scalar(p_hydro_state(jg)%prog(1)%vn,p_patch(jg), &
            & p_int_state(jg)%hex_north,p_hydro_state(jg)%diag%v)
        END SELECT

      CASE (inh_atmosphere)
        SELECT CASE (p_patch(jg)%cell_type)
        CASE (3)
          CALL rbf_vec_interpol_cell(p_nh_state(jg)%prog(1)%vn,p_patch(jg),&
            & p_int_state(jg),p_nh_state(jg)%diag%u,p_nh_state(jg)%diag%v)
        CASE (6)
          CALL edges2cells_scalar(p_nh_state(jg)%prog(1)%vn,p_patch(jg), &
            & p_int_state(jg)%hex_east,p_nh_state(jg)%diag%u)
          CALL edges2cells_scalar(p_nh_state(jg)%prog(1)%vn,p_patch(jg), &
            & p_int_state(jg)%hex_north,p_nh_state(jg)%diag%v)
        END SELECT

      CASE DEFAULT
      END SELECT
    ENDDO
    
    ! Note: here the derived output variables are not yet available
    ! (omega, divergence, vorticity)
    CALL write_output( time_config%current_datetime )
    l_have_output = .TRUE.

    END IF ! not lrestart

    !------------------------------------------------------------------
    ! Now start the time stepping:
    ! The special initial time step for the three time level schemes
    ! is executed within process_grid_level
    !------------------------------------------------------------------
    SELECT CASE (iequations)

    CASE (ishallow_water, ihs_atm_temp, ihs_atm_theta)
      CALL perform_ha_stepping( p_patch(1:), p_int_state(1:), p_grf_state(1:), &
                              & p_hydro_state, time_config%current_datetime,   &
                              & n_io, n_file, n_chkpt, n_diag, jfile,          &
                              & l_have_output                                  )

    CASE (inh_atmosphere)
      CALL perform_nh_stepping( p_patch, p_int_state, p_grf_state, p_nh_state,   &
                              & time_config%current_datetime,                    &
                              & n_io, n_file, n_chkpt, n_diag, l_have_output     )
    CASE DEFAULT
    END SELECT
 
    IF (ltimer) CALL print_timer
    
    !------------------------------------------------------------------
    !  cleaning up process
    !------------------------------------------------------------------
    
    CALL message(TRIM(routine),'start to clean up')
    
    ! Delete state variables

    SELECT CASE (iequations)

    CASE (ishallow_water, ihs_atm_temp, ihs_atm_theta)
      CALL destruct_icoham_dyn_state
      DEALLOCATE (p_hydro_state, STAT=ist)
      IF (ist /= SUCCESS) THEN
        CALL finish(TRIM(routine),'deallocation for p_hydro_state failed')
      ENDIF

    CASE (inh_atmosphere)

      CALL destruct_nh_state( p_nh_state )
      DEALLOCATE (p_nh_state, STAT=ist)
      IF (ist /= SUCCESS) THEN
        CALL finish(TRIM(routine),'deallocation for p_nh_state failed')
      ENDIF

    CASE DEFAULT
    END SELECT

    ! Delete output variable lists
    IF (l_have_output) CALL close_output_files
    
    ! Delete grids
    IF (n_dom > 1) THEN
      CALL destruct_2d_gridref_state( p_patch, p_grf_state )
    ENDIF

    IF (my_process_is_mpi_seq()  &
      & .OR. lrestore_states) THEN
      DEALLOCATE (p_grf_state_global, STAT=ist)
    ELSE
      DEALLOCATE (p_grf_state_subdiv, STAT=ist)
    ENDIF
    IF (ist /= SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for ptr_grf_state failed')
    ENDIF
      
    ! Deallocate memory for the parameterized forcing
    SELECT CASE (iforcing)

    CASE (inoforcing,iheldsuarez,ildf_dry)
      ! nothing to be done

    CASE (iecham,ildf_echam)
      CALL destruct_echam_phy_state  ! deallocate state vector
      CALL cleanup_echam_phy         ! deallocate parameter arrays

    CASE (inwp)
      CALL destruct_nwp_phy_state
      CALL destruct_nwp_lnd_state(p_lnd_state)

    CASE DEFAULT
      CALL finish(TRIM(routine),'iforcing has value that is not allowed')

    END SELECT
    
    ! Deallocate interpolation fields

    CALL destruct_2d_interpol_state( p_int_state )
    IF  (my_process_is_mpi_seq()  &
      & .OR. lrestore_states) THEN
      DEALLOCATE (p_int_state_global, STAT=ist)
    ELSE
      DEALLOCATE (p_int_state_subdiv, STAT=ist)
    ENDIF
    IF (ist /= SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for ptr_int_state failed')
    ENDIF
    
    ! Deallocate external data
    CALL destruct_ext_data

    ! deallocate ext_data array
    DEALLOCATE(ext_data, stat=ist)
    IF (ist/=success) THEN
      CALL finish(TRIM(routine), 'deallocation of ext_data')
    ENDIF
    
    ! Deallocate grid patches
    !
    CALL destruct_patches( p_patch )

    IF (my_process_is_mpi_seq()  &
      & .OR. lrestore_states) THEN
      DEALLOCATE( p_patch_global, STAT=ist )
    ELSE
      DEALLOCATE( p_patch_subdiv, STAT=ist )
    ENDIF
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocate for patch array failed')
    ENDIF
    
    ! deallocate output switches
    !
    
    CALL message(TRIM(routine),'clean-up finished')
    
  END SUBROUTINE atmo_model_old
  
END MODULE mo_atmo_model

