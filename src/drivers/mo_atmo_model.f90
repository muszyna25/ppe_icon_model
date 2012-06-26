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

USE mo_kind,                ONLY: wp, dp
USE mo_exception,           ONLY: message, finish, message_text
USE mo_mpi,                 ONLY: p_stop, &
  & my_process_is_io,  my_process_is_mpi_seq, my_process_is_mpi_test, &
  & my_process_is_mpi_parallel,                                       &
  & set_mpi_work_communicators, set_comm_input_bcast, null_comm_type, &
  & p_pe_work, get_my_mpi_all_id, p_min, p_max, p_comm_work,          &
  & process_mpi_io_size
USE mo_sync,                ONLY: enable_sync_checks, disable_sync_checks, &
                                  decomposition_statistics
USE mo_timer,               ONLY: init_timer, timer_start, timer_stop, &
  &                               timers_level, timer_model_init
USE mo_parallel_config,     ONLY: p_test_run, l_test_openmp, &
  &                               num_io_procs, nproma, use_icon_comm
USE mo_intp_lonlat,         ONLY: init_lonlat_grid_list,      &
  &                               compute_lonlat_intp_coeffs, &
  &                               destroy_lonlat_grid_list
USE mo_master_control,      ONLY: is_restart_run, get_my_process_name, get_my_model_no


USE mo_io_async,            ONLY: vlist_io_main_proc, &            ! main procedure for I/O PEs
                                  use_async_vlist_io
USE mo_name_list_output,    ONLY: name_list_io_main_proc
USE mo_name_list_output_config, ONLY: use_async_name_list_io, name_list_output_active


! Control parameters: run control, dynamics, i/o
!
USE mo_nonhydrostatic_config,ONLY: ivctype, kstart_moist, kstart_qv,     &
  &                                iadv_rcf, kend_qvsubstep, l_open_ubc, &
  &                                configure_nonhydrostatic
USE mo_dynamics_config,   ONLY: configure_dynamics, iequations
USE mo_run_config,        ONLY: configure_run, &
  & ltimer,               & !    :
  & iforcing,             & !    namelist parameter
  & ldump_states,         & ! flag if states should be dumped
  & lrestore_states,      & ! flag if states should be restored
  & ldump_dd, lread_dd,   &
  & nproc_dd,             &
  & nlev,nlevp1,          &
  & num_lev,num_levp1,    &
  & iqv, nshift,          &
  & lvert_nest, ntracer,  &
  & msg_level, dtime, output_mode

USE mo_impl_constants, ONLY:&
  & ihs_atm_temp,         & !    :
  & ihs_atm_theta,        & !    :
  & inh_atmosphere,       & !    :
  & ishallow_water,       &
  & max_dom, inwp


! For the coupling
USE mo_icon_cpl_init,      ONLY: icon_cpl_init
USE mo_icon_cpl_init_comp, ONLY: icon_cpl_init_comp
USE mo_impl_constants, ONLY: CELLS
USE mo_coupling_config,   ONLY : is_coupled_run, config_debug_coupler_level
USE mo_icon_cpl_def_grid, ONLY : ICON_cpl_def_grid, ICON_cpl_def_location
USE mo_icon_cpl_def_field, ONLY : ICON_cpl_def_field
USE mo_icon_cpl_search, ONLY : ICON_cpl_search
USE mo_alloc_patches,       ONLY : destruct_patches
USE mo_icon_cpl_finalize,   ONLY: icon_cpl_finalize

! Memory
!
USE mo_setup_subdivision,    ONLY: decompose_domain
USE mo_complete_subdivision, ONLY: &
  & complete_parallel_setup,       &
  & finalize_decomposition,        &
  & copy_processor_splitting,      &
  & set_patch_communicators,       &
  & setup_phys_patches

#ifndef NOMPI
USE mo_setup_subdivision,   ONLY:  npts_local
#endif

USE mo_dump_restore,        ONLY: dump_patch_state_netcdf,       &
& dump_domain_decomposition,     &
& restore_patches_netcdf,        &
& restore_interpol_state_netcdf, &
& restore_gridref_state_netcdf,  &
& dump_lonlat_data_netcdf,       &
& restore_lonlat_data_netcdf

USE mo_model_domain,        ONLY: t_patch, p_patch, p_patch_local_parent
USE mo_util_sysinfo,        ONLY: util_get_maxrss

! Horizontal grid
USE mo_grid_config,         ONLY: n_dom, n_dom_start, global_cell_type, &
                                  dynamics_parent_grid_id
USE mo_model_domimp_patches,ONLY: import_basic_patches, complete_patches
! Horizontal interpolation
!
USE mo_intp_state,          ONLY: construct_2d_interpol_state,  &
  &                               destruct_2d_interpol_state,   &
  &                               transfer_interpol_state
USE mo_grf_intp_state,      ONLY: construct_2d_gridref_state,   &
  &                               destruct_2d_gridref_state, transfer_grf_state

! Vertical grid
!
USE mo_vertical_coord_table,ONLY: init_vertical_coord_table, &
  &                               vct_a, vct_b, apzero
USE mo_nh_init_utils,       ONLY: init_hybrid_coord, init_sleve_coord

! Parameterized forcing
!
USE mo_impl_constants,      ONLY: SUCCESS, MAX_CHAR_LENGTH

! Time integration
!
! External data
USE mo_ext_data_state,      ONLY: ext_data, init_ext_data, init_index_lists, &
  &                               destruct_ext_data


!-------------------------------------------------------------------------
! to break circular dependency

USE mo_intp_data_strc,      ONLY: p_int_state, p_int_state_local_parent
USE mo_grf_intp_data_strc,  ONLY: p_grf_state, p_grf_state_local_parent

!-------------------------------------------------------------------------
USE mo_io_restart,           ONLY: read_restart_info_file
USE mo_io_restart_namelist,  ONLY: read_restart_namelists, &
  &                                delete_restart_namelists
USE mo_io_restart_attributes,ONLY: read_restart_attributes

USE mo_read_namelists,     ONLY: read_atmo_namelists
USE mo_nml_crosscheck,       ONLY: atm_crosscheck

!USE mo_interpol_config,    ONLY: configure_interpolation 
USE mo_interpol_config
USE mo_advection_config,   ONLY: configure_advection
USE mo_diffusion_config,   ONLY: configure_diffusion

USE mo_atmo_hydrostatic,    ONLY: atmo_hydrostatic 
USE mo_atmo_nonhydrostatic, ONLY: atmo_nonhydrostatic 
  
USE mo_icon_comm_interface,  ONLY: construct_icon_communication, &
    & destruct_icon_communication

!-------------------------------------------------------------------------
IMPLICIT NONE
PRIVATE

PUBLIC :: atmo_model, construct_atmo_model, destruct_atmo_model

CONTAINS
  
  !-------------------------------------------------------------------
  !>
  SUBROUTINE atmo_model(atm_namelist_filename,shr_namelist_filename)

    CHARACTER(LEN=*), INTENT(in) :: atm_namelist_filename
    CHARACTER(LEN=*), INTENT(in) :: shr_namelist_filename

    CHARACTER(*), PARAMETER :: routine = "mo_atmo_model:atmo_model"

#ifndef NOMPI    
#if defined(__SX__)
    INTEGER  :: maxrss
    REAL(dp) :: maxrss_gridpt,maxrss_gridpt_min, maxrss_gridpt_max
#endif
#endif


    !---------------------------------------------------------------------
    ! construct the atmo model
    CALL construct_atmo_model(atm_namelist_filename,shr_namelist_filename)
    
    !---------------------------------------------------------------------
    ! construct the coupler
    IF ( is_coupled_run() ) THEN
      CALL construct_atmo_coupler()
    ENDIF 
    
    !---------------------------------------------------------------------
    ! 12. The hydrostatic and nonhydrostatic models branch from this point
    !---------------------------------------------------------------------
    SELECT CASE(iequations)
    CASE(ishallow_water,ihs_atm_temp,ihs_atm_theta)
      CALL atmo_hydrostatic

    CASE(inh_atmosphere)
      CALL atmo_nonhydrostatic

    CASE DEFAULT
      CALL finish( TRIM(routine),'unknown choice for iequaions.')
    END SELECT


    !---------------------------------------------------------------------
    ! 13. Integration finished. Carry out the shared clean-up processes
    !---------------------------------------------------------------------
    CALL destruct_atmo_model()
  
    !---------------------------------------------------------------------
    ! destruct the coupler
    IF ( is_coupled_run() ) CALL ICON_cpl_finalize
    
    !---------------------------------------------------------------------
    ! (optional:) write resident set size from OS
#ifndef NOMPI    
#if defined(__SX__)
    IF (msg_level >= 16) THEN
      CALL util_get_maxrss(maxrss)
      PRINT  *, "PE #", get_my_mpi_all_id(), &
        &    ": MAXRSS (MiB) = ", maxrss
      ! compute memory consumption per grid point (w/out halo points)
      maxrss_gridpt = REAL(maxrss)/SUM(npts_local(n_dom_start:n_dom, 1:3))
      maxrss_gridpt_min = p_min(zfield=maxrss_gridpt, comm=p_comm_work)
      maxrss_gridpt_max = p_max(zfield=maxrss_gridpt, comm=p_comm_work)
      WRITE (message_text,'(a,a,f5.3,a,a,f5.3)')   &
        & "memory consumption (MiB/grid point): ", &
        & "min = ", maxrss_gridpt_min, " / ", &
        & "max = ", maxrss_gridpt_max
      CALL message(TRIM(routine),message_text)
    END IF
#endif
#endif

  END SUBROUTINE atmo_model
  !-------------------------------------------------------------------

  !-------------------------------------------------------------------
  !>
  SUBROUTINE construct_atmo_model(atm_namelist_filename,shr_namelist_filename)

    CHARACTER(LEN=*), INTENT(in) :: atm_namelist_filename
    CHARACTER(LEN=*), INTENT(in) :: shr_namelist_filename

    CHARACTER(*), PARAMETER :: routine = "mo_atmo_model:init_atmo_model"
    CHARACTER(LEN=MAX_CHAR_LENGTH) :: grid_file_name 
    LOGICAL :: lsuccess
    INTEGER :: jg, jgp

    INTEGER :: error_status

    TYPE(t_patch), ALLOCATABLE :: p_patch_global(:)


    ! initialize global registry of lon-lat grids
    CALL init_lonlat_grid_list()

    !---------------------------------------------------------------------
    ! 0. If this is a resumed or warm-start run...
    !---------------------------------------------------------------------
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

      CALL read_restart_namelists('restart_atm.nc')
      CALL message(TRIM(routine), 'read namelists from atm restart file')

      ! Read all global attributs in the restart file and store them in a buffer.

      CALL read_restart_attributes('restart_atm.nc')
      CALL message(TRIM(routine), 'read global attributes from atm restart file')

    END IF ! is_restart_run()

    !---------------------------------------------------------------------
    ! 1.1 Read namelists (newly) specified by the user; fill the 
    !     corresponding sections of the configuration states.
    !---------------------------------------------------------------------

    CALL read_atmo_namelists(atm_namelist_filename,shr_namelist_filename)

    !---------------------------------------------------------------------
    ! 1.2 Cross-check namelist setups
    !---------------------------------------------------------------------

    CALL atm_crosscheck

    !---------------------------------------------------------------------
    ! 2. Call configure_run to finish filling the run_config state.
    !    This needs to be done very early (but anyway after atm_crosscheck)
    !    because some component of the state, e.g., num_lev, may be 
    !    modified in this subroutine which affect the following CALLs.
    !---------------------------------------------------------------------
    SELECT CASE(iequations)
    CASE (inh_atmosphere)
      CALL configure_run( iadv_rcf )
    CASE DEFAULT
      CALL configure_run    
    END SELECT

    !-------------------------------------------------------------------
    ! 3.1 Initialize the mpi work groups
    !-------------------------------------------------------------------
    CALL set_mpi_work_communicators(p_test_run, l_test_openmp, num_io_procs)

    !-------------------------------------------------------------------
    ! 3.2 Initialize various timers
    !-------------------------------------------------------------------
    IF (ltimer) CALL init_timer
    IF (timers_level > 3) CALL timer_start(timer_model_init)

    !------------------
    ! Next, define the horizontal and vertical grids since they are aready
    ! needed for some derived control parameters. This includes
    ! - patch import
    ! - domain decompistion
    ! - vertical coordinates
    !-------------------------------------------------------------------
    ! 4. Import patches
    !-------------------------------------------------------------------
    ! If we belong to the I/O PEs just call xxx_io_main_proc before reading patches.
    ! This routine will never return

    IF (process_mpi_io_size > 0) THEN
      ! Decide whether async vlist or name_list IO is to be used,
      ! only one of both may be enabled!

      ! Currently, async name_list IO is used if at least one name_list was read
      IF(name_list_output_active) THEN
        use_async_name_list_io = .TRUE.
        CALL message('','asynchronous namelist I/O scheme is enabled.')
        IF (my_process_is_io() .AND. (.NOT. my_process_is_mpi_test())) THEN
          CALL name_list_io_main_proc
        END IF
        IF (output_mode%l_vlist) THEN
          output_mode%l_vlist = .FALSE.
          CALL message('','vlist I/O scheme has been disabled because combining it&
            & with namelist I/O is not possible for asynchronous output.')
        ENDIF
      ELSE
        use_async_vlist_io = .TRUE.
        CALL message('','asynchronous vlist I/O scheme is enabled.')
        IF (my_process_is_io()) CALL vlist_io_main_proc
      ENDIF
    ELSE
      IF(name_list_output_active) THEN
        CALL message('','synchronous namelist I/O scheme is enabled.')
      ENDIF
      IF (output_mode%l_vlist) THEN
        CALL message('','synchronous vlist I/O scheme is enabled.')
      ENDIF
    ENDIF

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

      IF(my_process_is_mpi_parallel()) THEN

        IF(lread_dd) THEN
          CALL restore_patches_netcdf( p_patch, .FALSE. )
        ELSE
          ALLOCATE(p_patch_global(n_dom_start:n_dom))
          CALL import_basic_patches(p_patch_global,nlev,nlevp1,num_lev,num_levp1,nshift)
          CALL decompose_domain(p_patch_global)
          DEALLOCATE(p_patch_global)
        ENDIF
        IF(ldump_dd) THEN
          DO jg = n_dom_start, n_dom
            IF(jg==n_dom_start) THEN
              CALL dump_domain_decomposition(p_patch(jg))
            ELSE
              CALL dump_domain_decomposition(p_patch(jg),p_patch_local_parent(jg))
            ENDIF
          ENDDO
          CALL p_stop
          STOP
        ENDIF

        CALL complete_parallel_setup()

      ELSE IF(my_process_is_mpi_test()) THEN

        ! Test run:
        ! If lread_dd is not set, import_basic_patches() goes parallel with the call above,
        ! actually the Test PE reads and broadcasts to the others.
        ! IF lread_dd is set, the read is in standalone mode.

        IF(lread_dd) CALL set_comm_input_bcast(null_comm_type)
        CALL import_basic_patches(p_patch,nlev,nlevp1,num_lev,num_levp1,nshift)
        IF(lread_dd) CALL set_comm_input_bcast()

        IF(ldump_dd) THEN
          ! ldump_dd makes not much sense in a test run, we simply stop
          CALL p_stop
          STOP
        ENDIF

      ELSE

        ! Single processor run

        IF(ldump_dd) THEN
          ! If ldump_dd is set in a single processor run, a domain decomposition for
          ! nproc_dd processors is done
          ALLOCATE(p_patch_global(n_dom_start:n_dom))
          CALL import_basic_patches(p_patch_global,nlev,nlevp1,num_lev,num_levp1,nshift)
          CALL decompose_domain(p_patch_global, nproc_dd)
          DEALLOCATE(p_patch_global)
          ! We are done, the dump is done within decompose_domain
          CALL p_stop
          STOP
        ELSE
          CALL import_basic_patches(p_patch,nlev,nlevp1,num_lev,num_levp1,nshift)
        ENDIF

      ENDIF

      ! Complete information which is not yet read or calculated
      CALL complete_patches( p_patch )

    ENDIF

    ! In case of a test run: Copy processor splitting to test PE
    IF(p_test_run) CALL copy_processor_splitting(p_patch)
    !--------------------------------------------------------------------------------
    
    !--------------------------------------------------------------------------------
    ! 5. Construct interpolation state, compute interpolation coefficients.
    !--------------------------------------------------------------------------------

    CALL configure_interpolation( global_cell_type, n_dom, p_patch(1:)%level )

    ! Allocate array for interpolation state

    ALLOCATE( p_int_state(n_dom_start:n_dom), &
            & p_grf_state(n_dom_start:n_dom), &
            & STAT=error_status)
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

    IF(lrestore_states) THEN
      ! Interpolation state is read from NetCDF
      ! On the test PE it is constructed at the same time to be able to check
      ! the state read in with the constructed state
      IF( .NOT. my_process_is_mpi_test()) THEN
        CALL restore_interpol_state_netcdf(p_patch, p_int_state )
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

    !-----------------------------------------------------------------------------
    ! 6. Construct grid refinment state, compute coefficients
    !-----------------------------------------------------------------------------
    ! For the NH model, the initialization routines called from
    ! construct_2d_gridref_state require the metric terms to be present

    IF (n_dom_start==0 .OR. n_dom > 1) THEN
      IF(lrestore_states) THEN
        ! Gridref state is read from NetCDF
        ! On the test PE it is constructed at the same time to be able to check
        ! the state read in with the constructed state
        IF( .NOT. my_process_is_mpi_test()) THEN
          CALL restore_gridref_state_netcdf(p_patch, p_grf_state)
        ELSE
          ! construct_2d_gridref_state makes sync calls for checking the
          ! results on the parallel PEs to the results on the test PE.
          ! These checks must be disabled here!
          CALL disable_sync_checks
          CALL construct_2d_gridref_state (p_patch, p_grf_state)
          CALL enable_sync_checks
        ENDIF
      ELSE
        ! Construct gridref state
        ! For non-parallel runs (either 1 PE only or on the test PE) this is done as usual,
        ! for parallel runs, the main part of the gridref state is constructed on the
        ! local parent with the following call
        CALL construct_2d_gridref_state (p_patch, p_grf_state)
        IF(my_process_is_mpi_parallel()) THEN
          ! Transfer gridref state from local parent to p_grf_state
          DO jg = n_dom_start+1, n_dom
            jgp = p_patch(jg)%parent_id
            CALL transfer_grf_state(p_patch(jgp), p_patch_local_parent(jg), &
                                 &  p_grf_state(jgp), p_grf_state_local_parent(jg), &
                                 &  p_patch(jg)%parent_child_index)
          ENDDO
        ENDIF
      ENDIF
    ENDIF

    !-------------------------------------------------------------------
    ! 7. Finalize domain decomposition
    !-------------------------------------------------------------------   

    IF(my_process_is_mpi_parallel() .AND. .NOT.lrestore_states) THEN

      CALL finalize_decomposition()

      IF(.NOT.p_test_run) THEN ! the call below hangs in test mode
        ! Print diagnostic information about domain decomposition
        DO jg = 1, n_dom
          CALL decomposition_statistics(p_patch(jg))
        ENDDO
      ENDIF

    ENDIF
    !--------------------------------------------------------------------------------
    
    !--------------------------------------------        
    ! Setup the information for the physical patches
    CALL setup_phys_patches

    !-------------------------------------------------------------------
    ! 7b. Constructing data for lon-lat interpolation
    !-------------------------------------------------------------------

    IF(lrestore_states) THEN
      ! restore data from file(s):
      CALL restore_lonlat_data_netcdf(p_patch)
    ELSE
      CALL compute_lonlat_intp_coeffs(p_patch(1:), p_int_state(1:))
    END IF

    ! Dump divided patches with interpolation and grf state to NetCDF
    ! file and exit
    IF(ldump_states)THEN
      
      CALL message(TRIM(routine),'ldump_states is set: '//&
                  'dumping patches+states and finishing')

      IF(.NOT. my_process_is_mpi_test()) THEN
        DO jg = n_dom_start, n_dom
          CALL dump_patch_state_netcdf(p_patch(jg),p_int_state(jg),p_grf_state(jg))
        ENDDO
      ENDIF

      ! dump the setup of all lon-lat interpolation processes.
      CALL dump_lonlat_data_netcdf(p_patch)

      CALL p_stop
      STOP

    ENDIF

   !---------------------------------------------------------------------
   ! 8. Import vertical grid/ define vertical coordinate
   !---------------------------------------------------------------------
    SELECT CASE (iequations)

    CASE (ishallow_water)
      CALL init_vertical_coord_table(iequations, p_patch(1)%nlev)

    CASE (ihs_atm_temp, ihs_atm_theta)
      CALL init_vertical_coord_table(iequations, p_patch(1)%nlev)

    CASE (inh_atmosphere)
      IF (ivctype == 1) THEN
        CALL init_hybrid_coord(p_patch(1)%nlev)
      ELSE IF (ivctype == 2) THEN
        CALL init_sleve_coord(p_patch(1)%nlev)
      ENDIF

    CASE DEFAULT
    END SELECT

    !---------------------------------------------------------------------
    ! 9. Horizontal and vertical grid(s) are now defined.
    !    Assign values to derived variables in the configuration states
    !---------------------------------------------------------------------

    CALL configure_dynamics ( n_dom )
    CALL configure_diffusion( n_dom, dynamics_parent_grid_id, &
      &                       nlev, vct_a, vct_b, apzero      )

    IF (iequations == inh_atmosphere) THEN
      DO jg =1,n_dom
        CALL configure_nonhydrostatic( jg, p_patch(jg)%nlev,     &
          &                            p_patch(jg)%nshift_total  )
      ENDDO
    ENDIF

    DO jg =1,n_dom
     CALL configure_advection( jg, p_patch(jg)%nlev, p_patch(1)%nlev,      &
       &                      iequations, iforcing, iqv, kstart_moist(jg), &
       &                      kstart_qv(jg), kend_qvsubstep(jg),           &
       &                      lvert_nest, l_open_ubc, ntracer ) 
    ENDDO

    !------------------------------------------------------------------
    ! 10. Create and optionally read external data fields
    !------------------------------------------------------------------
    ALLOCATE (ext_data(n_dom), STAT=error_status)
    IF (error_status /= SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for ext_data failed')
    ENDIF

    ! allocate memory for atmospheric/oceanic external data and
    ! optionally read those data from netCDF file.
    CALL init_ext_data (p_patch(1:), p_int_state(1:), ext_data)
    !
    ! generation of tiles index lists
    IF ( iequations == inh_atmosphere .AND. iforcing == inwp ) THEN
      CALL init_index_lists (p_patch(1:), ext_data)
    ENDIF

    !-------------------------------------------------------------------
    ! Initialize icon_comm_lib
    !-------------------------------------------------------------------
    IF (use_icon_comm) THEN
      CALL construct_icon_communication()
    ENDIF
    
    IF (timers_level > 3) CALL timer_stop(timer_model_init)

  END SUBROUTINE construct_atmo_model

  !-------------------------------------------------------------------
  !>
  SUBROUTINE destruct_atmo_model()
    
    CHARACTER(*), PARAMETER :: routine = "mo_atmo_model:destruct_atmo_model"

    INTEGER :: error_status
    ! Destruct external data state

    CALL destruct_ext_data
    IF (msg_level > 5) CALL message(TRIM(routine),'destruct_ext_data is done')

    ! deallocate ext_data array
    DEALLOCATE(ext_data, stat=error_status)
    IF (error_status/=success) THEN
      CALL finish(TRIM(routine), 'deallocation of ext_data')
    ENDIF

    ! Deconstruct grid refinement state

    IF (n_dom > 1) THEN
      CALL destruct_2d_gridref_state( p_patch, p_grf_state )
    ENDIF
    IF (msg_level > 5) CALL message(TRIM(routine),'destruct_2d_gridref_state is done')

    DEALLOCATE (p_grf_state, STAT=error_status)
    IF (error_status /= SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for ptr_grf_state failed')
    ENDIF

    ! Deallocate interpolation fields

    CALL destruct_2d_interpol_state( p_int_state )
    IF (msg_level > 5) CALL message(TRIM(routine),'destruct_2d_interpol_state is done')

    DEALLOCATE (p_int_state, STAT=error_status)
    IF (error_status /= SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for ptr_int_state failed')
    ENDIF

    ! Deallocate global registry for lon-lat grids
    CALL destroy_lonlat_grid_list()

    ! Deallocate grid patches
    CALL destruct_patches( p_patch )
    IF (msg_level > 5) CALL message(TRIM(routine),'destruct_patches is done')

    DEALLOCATE( p_patch, STAT=error_status )
    IF (error_status/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocate for patch array failed')
    ENDIF

    ! clear restart namelist buffer
    CALL delete_restart_namelists()
    IF (msg_level > 5) CALL message(TRIM(routine),'delete_restart_namelists is done')
    
    IF (use_icon_comm) THEN
      CALL destruct_icon_communication()
    ENDIF
    
    CALL message(TRIM(routine),'clean-up finished')
    
  END SUBROUTINE destruct_atmo_model
  !-------------------------------------------------------------------
 
  !-------------------------------------------------------------------
  !>
  SUBROUTINE construct_atmo_coupler()
    ! For the coupling

    INTEGER, PARAMETER :: no_of_fields = 8

    CHARACTER(LEN=MAX_CHAR_LENGTH) ::  field_name(no_of_fields)
    INTEGER :: field_id(no_of_fields)
    INTEGER :: grid_id
    INTEGER :: grid_shape(2) 
    INTEGER :: field_shape(3) 
    INTEGER :: i, patch_no, error_status
    
    !---------------------------------------------------------------------
    ! 11. Do the setup for the coupled run
    !
    ! For the time being this could all go into a subroutine which is
    ! common to atmo and ocean. Does this make sense if the setup deviates
    ! too much in future.
    !---------------------------------------------------------------------
    IF ( is_coupled_run() ) THEN

      !------------------------------------------------------------
      CALL icon_cpl_init(debug_level=config_debug_coupler_level)
      ! Inform the coupler about what we are
      CALL icon_cpl_init_comp ( get_my_process_name(), get_my_model_no(), error_status )
      ! split the global_mpi_communicator into the components
      !------------------------------------------------------------
      
      patch_no = 1

      grid_shape(1)=1
      grid_shape(2)=p_patch(patch_no)%n_patch_cells

      ! CALL get_patch_global_indexes ( patch_no, CELLS, no_of_entities, grid_glob_index )
      ! should grid_glob_index become a pointer in ICON_cpl_def_grid as well?
      CALL ICON_cpl_def_grid ( &
        & grid_shape, p_patch(patch_no)%cells%glb_index, & ! input
        & grid_id, error_status )                          ! output
      
      ! Marker for internal and halo points, a list which contains the
      ! rank where the native cells are located.
      CALL ICON_cpl_def_location ( &
        & grid_id, grid_shape, p_patch(patch_no)%cells%owner_local, & ! input
        & p_pe_work,  & ! this owner id
        & error_status )                                            ! output

      field_name(1) = "TAUX"
      field_name(2) = "TAUY"
      field_name(3) = "SFWFLX"   ! bundled field containing two flux components
      field_name(4) = "SFTEMP"
      field_name(5) = "THFLX"    ! bundled field containing four flux components
      field_name(6) = "SST"
      field_name(7) = "OCEANU"
      field_name(8) = "OCEANV"

      field_shape(1:2) = grid_shape(1:2)

      DO i = 1, no_of_fields
         IF ( i == 3 ) THEN
           field_shape(3) = 2
         ELSE  IF ( i == 5 ) THEN
           field_shape(3) = 4
         ELSE
           field_shape(3) = 1
         ENDIF
         CALL ICON_cpl_def_field ( field_name(i), grid_id, field_id(i), &
                                 & field_shape, error_status )
      ENDDO

      CALL ICON_cpl_search

    ENDIF
    
  END SUBROUTINE construct_atmo_coupler

END MODULE mo_atmo_model
