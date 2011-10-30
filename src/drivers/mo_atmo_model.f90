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

USE mo_kind,                ONLY: wp
USE mo_exception,           ONLY: message, finish, message_text
USE mo_mpi,                 ONLY: p_stop, &
  & my_process_is_io,  my_process_is_mpi_seq, my_process_is_mpi_test, &
  & my_process_is_mpi_parallel,                                       &
  & set_mpi_work_communicators, set_comm_input_bcast, null_comm_type, &
  & p_pe_work
USE mo_sync,                ONLY: enable_sync_checks, disable_sync_checks
USE mo_timer,               ONLY: init_timer, timer_start, timer_stop, &
  &                               timer_lonlat_setup
USE mo_parallel_config,     ONLY: p_test_run, l_test_openmp, &
  &                               num_io_procs, nproma
USE mo_lonlat_intp_config,  ONLY: configure_lonlat_intp
USE mo_master_control,      ONLY: is_restart_run, get_my_couple_id


USE mo_io_async,            ONLY: io_main_proc            ! main procedure for I/O PEs

! Control parameters: run control, dynamics, i/o
!
USE mo_nonhydrostatic_config,ONLY: ivctype, kstart_moist, kstart_qv,    &
  &                                kend_qvsubstep, l_open_ubc,          &
  &                                configure_nonhydrostatic
USE mo_dynamics_config,   ONLY: iequations
USE mo_run_config,        ONLY: configure_run, &
  & ltimer,               & !    :
  & iforcing,             & !    namelist parameter
  & ldump_states,         & ! flag if states should be dumped
  & lrestore_states,      & ! flag if states should be restored
  & nlev,nlevp1,          &
  & num_lev,num_levp1,    &
  & iqv, nshift,          &
  & lvert_nest, ntracer

USE mo_impl_constants, ONLY:&
  & ihs_atm_temp,         & !    :
  & ihs_atm_theta,        & !    :
  & inh_atmosphere,       & !    :
  & ishallow_water,       &
  & max_dom


! For the coupling
USE mo_impl_constants, ONLY: CELLS
USE mo_master_control, ONLY : is_coupled_run
USE mo_icon_cpl_init_comp, ONLY : get_my_local_comp_id
USE mo_icon_cpl_def_grid, ONLY : ICON_cpl_def_grid, ICON_cpl_def_location
USE mo_icon_cpl_def_field, ONLY : ICON_cpl_def_field
USE mo_icon_cpl_search, ONLY : ICON_cpl_search
USE mo_model_domain_import, ONLY : get_patch_global_indexes

! Memory
!
USE mo_subdivision,         ONLY: decompose_domain,         &
& finalize_decomposition, &
& copy_processor_splitting,      &
& set_patch_communicators
USE mo_dump_restore,        ONLY: dump_patch_state_netcdf,       &
& restore_patches_netcdf,        &
& restore_interpol_state_netcdf, &
& restore_gridref_state_netcdf

USE mo_model_domain,        ONLY: p_patch_global, p_patch_subdiv, p_patch, &
  &                               p_patch_local_parent

! Horizontal grid
USE mo_grid_config,         ONLY: n_dom, n_dom_start, global_cell_type, &
                                  dynamics_parent_grid_id
USE mo_model_domain_import, ONLY: import_patches, destruct_patches

! Horizontal interpolation
!
USE mo_intp_state,          ONLY: construct_2d_interpol_state,  &
  &                               destruct_2d_interpol_state,   &
  &                               allocate_int_state_lonlat,    &
  &                               rbf_setup_interpol_lonlat,    &
  &                               deallocate_int_state_lonlat,  &
  &                               transfer_interpol_state
USE mo_grf_interpolation,   ONLY: construct_2d_gridref_state,   &
  &                               destruct_2d_gridref_state
USE mo_grf_intp_state,      ONLY: transfer_grf_state
USE mo_lonlat_intp_config,  ONLY: lonlat_intp_config

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
USE mo_ext_data,            ONLY: ext_data, init_ext_data, destruct_ext_data

!  USE mo_nwp_phy_init,          ONLY: init_nwp_phy
!!$  USE mo_gscp_cosmo,          ONLY: hydci_pp_init

! meteogram output
USE mo_mtgrm_output,        ONLY: mtgrm_init, mtgrm_finalize
USE mo_mtgrm_config,        ONLY: mtgrm_output_config


!-------------------------------------------------------------------------
! to break circular dependency

USE mo_intp_data_strc,      ONLY: p_int_state, p_int_state_local_parent,     &
  &                               p_int_state_lonlat_global,                 &
  &                               p_int_state_lonlat
USE mo_grf_intp_data_strc,  ONLY: p_grf_state, p_grf_state_local_parent

!-------------------------------------------------------------------------
USE mo_io_restart,           ONLY: read_restart_info_file
USE mo_io_restart_namelist,  ONLY: read_restart_namelists, &
  &                                delete_restart_namelists
USE mo_io_restart_attributes,ONLY: read_restart_attributes

USE mo_read_namelists,     ONLY: read_atmo_namelists
USE mo_nml_crosscheck,       ONLY: atm_crosscheck

USE mo_dynamics_config,    ONLY: configure_dynamics  ! subroutine
!USE mo_interpol_config,    ONLY: configure_interpolation 
USE mo_interpol_config
USE mo_advection_config,   ONLY: configure_advection
USE mo_diffusion_config,   ONLY: configure_diffusion

USE mo_atmo_hydrostatic,    ONLY: atmo_hydrostatic 
USE mo_atmo_nonhydrostatic, ONLY: atmo_nonhydrostatic 

!-------------------------------------------------------------------------
IMPLICIT NONE
PRIVATE

PUBLIC :: atmo_model

CONTAINS
!>
!!
  SUBROUTINE atmo_model(atm_namelist_filename,shr_namelist_filename)

    CHARACTER(LEN=*), INTENT(in) :: atm_namelist_filename
    CHARACTER(LEN=*), INTENT(in) :: shr_namelist_filename

    CHARACTER(LEN=MAX_CHAR_LENGTH) :: grid_file_name 
    CHARACTER(*), PARAMETER :: routine = "mo_atmo_model:atmo_model"
    LOGICAL :: lsuccess, l_construct_lonlat_coeffs
    INTEGER :: jg, jgp

    ! For the coupling

    INTEGER, PARAMETER :: no_of_fields = 8

    CHARACTER(LEN=MAX_CHAR_LENGTH) ::  field_name(no_of_fields)
    INTEGER :: field_id(no_of_fields)
    INTEGER :: comp_id
    INTEGER :: grid_id
    INTEGER :: grid_shape(2) 
    INTEGER :: field_shape(3) 
    INTEGER :: i, error_status
    INTEGER :: patch_no

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
    CALL configure_run

    ! configure (optional) lon-lat interpolation of output variables
    CALL configure_lonlat_intp(max_dom, nproma, num_io_procs)

    !-------------------------------------------------------------------
    ! 3.1 Initialize the mpi work groups
    !-------------------------------------------------------------------
    CALL set_mpi_work_communicators(p_test_run, l_test_openmp, num_io_procs)

    !-------------------------------------------------------------------
    ! 3.2 Initialize various timers
    !-------------------------------------------------------------------
    IF (ltimer) CALL init_timer

    !------------------
    ! Next, define the horizontal and vertical grids since they are aready
    ! needed for some derived control parameters. This includes
    ! - patch import
    ! - domain decompistion
    ! - vertical coordinates
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
                            nlev,nlevp1,num_lev,num_levp1,nshift)
      ENDIF
      ! After the restore is done set p_comm_input_bcast in the
      ! same way as it would be set in parallel_nml_setup when
      ! no restore is wanted:
      CALL set_comm_input_bcast()

      p_patch => p_patch_global
    ELSE
      CALL import_patches( p_patch_global,                       &
                           nlev,nlevp1,num_lev,num_levp1,nshift)
      IF(my_process_is_mpi_parallel()) then
        CALL decompose_domain()
        p_patch => p_patch_subdiv
      ELSE
        p_patch => p_patch_global
      ENDIF
    ENDIF
    ! Note: from this point the p_patch is used

    !--------------------------------------------------------------------------------
    ! 5. Construct interpolation state, compute interpolation coefficients.
    !--------------------------------------------------------------------------------

    CALL configure_interpolation( global_cell_type, n_dom, p_patch_global(1:)%level )

    ! Allocate array for interpolation state

    ALLOCATE( p_int_state(n_dom_start:n_dom), &
            & p_grf_state(n_dom_start:n_dom), &
            & STAT=error_status)
    IF (error_status /= SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for ptr_int_state failed')
    ENDIF

    ALLOCATE(p_int_state_lonlat_global(n_dom_start:n_dom), &
      & STAT=error_status)
    IF (error_status /= SUCCESS) &
      CALL finish(TRIM(routine),'allocation for lon-lat int_state failed')

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
        CALL restore_interpol_state_netcdf(p_patch, p_int_state, &
          &                                p_int_state_lonlat_global )
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
          CALL restore_gridref_state_netcdf(p_patch_global, p_grf_state)
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
    IF (my_process_is_mpi_seq()  &
      &  .OR. lrestore_states) THEN

      CALL set_patch_communicators(p_patch)

    ELSE

      CALL finalize_decomposition()

    ENDIF

    ! In case of a test run: Copy processor splitting to test PE
    IF(p_test_run) CALL copy_processor_splitting(p_patch)

    !-------------------------------------------------------------------
    ! 7b. Constructing data for lon-lat interpolation
    !-------------------------------------------------------------------

    DO jg = 1, n_dom

      l_construct_lonlat_coeffs = &
        &  lonlat_intp_config(jg)%l_enabled  .AND. &
        &  (.NOT. lrestore_states .OR. my_process_is_mpi_test())

      LONLAT : IF (l_construct_lonlat_coeffs) THEN
        
        CALL allocate_int_state_lonlat(jg, p_int_state_lonlat_global(jg) )
        
        ! compute coefficients needed for lon-lat interpolation
        ! note: coefficients are only required on IO PEs
        IF (ltimer) CALL timer_start(timer_lonlat_setup)
        ! Note: Up to now, only the nearest neighbor search operates
        ! on the distributed patches. Therefore, we have to provide
        ! the global patch as well:
        CALL rbf_setup_interpol_lonlat( jg, p_patch(jg),               &
          &                             p_int_state_lonlat_global(jg), &
          &                             p_int_state(jg) )
        IF (ltimer) CALL timer_stop(timer_lonlat_setup)
      END IF LONLAT

    END DO

    p_int_state_lonlat => p_int_state_lonlat_global

    ! Dump divided patches with interpolation and grf state to NetCDF
    ! file and exit
    ! Note: Dump of lonlat interpolation coefficients only makes sense
    ! with asynchronous IO *disabled* because otherwise the necessary
    ! coefficients were assembled on the pure IO PEs.
    IF(ldump_states)THEN
      
      CALL message(TRIM(routine),'ldump_states is set: '//&
                  'dumping patches+states and finishing')

      IF(.NOT. my_process_is_mpi_test()) THEN
        DO jg = n_dom_start, n_dom
          CALL dump_patch_state_netcdf(jg, p_patch(jg),p_int_state(jg),p_grf_state(jg), &
            &                          p_int_state_lonlat_global(jg))
        ENDDO
      ENDIF

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
                            & nlev, vct_a, vct_b, apzero      )

    IF (iequations == inh_atmosphere) THEN
      DO jg =1,n_dom
        CALL configure_nonhydrostatic(jg, p_patch(jg)%nlev, p_patch(jg)%nshift_total)
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

    !---------------------------------------------------------------------
    ! 11. Do the setup for the coupled run
    !
    ! For the time being this could all go into a subroutine which is
    ! common to atmo and ocean. Does this make sense if the setup deviates
    ! too much in future.
    !---------------------------------------------------------------------
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
      
      ! Marker for internal and halo points, a list which contains the
      ! rank where the native cells are located.
      CALL ICON_cpl_def_location ( &
        & grid_id, grid_shape, p_patch(patch_no)%cells%owner_local, & ! input
        & p_pe_work,  & ! this owner id
        & error_status )                                            ! output

      field_name(1) = "TAUX"
      field_name(2) = "TAUY"
      field_name(3) = "SFWFLX"
      field_name(4) = "SFTEMP"
      field_name(5) = "THFLX"
      field_name(6) = "SST"
      field_name(7) = "OCEANU"
      field_name(8) = "OCEANV"

      field_shape(1:2) = grid_shape(1:2)
      field_shape(3)   = 1

      DO i = 1, no_of_fields
         CALL ICON_cpl_def_field ( field_name(i), comp_id, grid_id, field_id(i), &
                                 & field_shape, error_status )
      ENDDO

      CALL ICON_cpl_search

    ENDIF

    !---------------------------------------------------------------------
    !     Setup of meteogram output
    !---------------------------------------------------------------------

    DO jg =1,n_dom
      IF (mtgrm_output_config(jg)%lenabled) THEN
        CALL mtgrm_init(mtgrm_output_config(jg), p_patch(jg), &
          &             ext_data(jg), iforcing, jg)
      END IF
    END DO

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
    ! Destruct external data state

    CALL destruct_ext_data

    ! deallocate ext_data array
    DEALLOCATE(ext_data, stat=error_status)
    IF (error_status/=success) THEN
      CALL finish(TRIM(routine), 'deallocation of ext_data')
    ENDIF

    ! Deconstruct grid refinement state

    IF (n_dom > 1) THEN
      CALL destruct_2d_gridref_state( p_patch, p_grf_state )
    ENDIF

    DEALLOCATE (p_grf_state, STAT=error_status)
    IF (error_status /= SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for ptr_grf_state failed')
    ENDIF

    ! Deallocate interpolation fields

    CALL destruct_2d_interpol_state( p_int_state )
    IF (l_construct_lonlat_coeffs) THEN
      DO jg = 1, n_dom
        IF (lonlat_intp_config(jg)%l_enabled) &
          CALL deallocate_int_state_lonlat(p_int_state_lonlat_global(jg))
      ENDDO
    END IF

    IF  (my_process_is_mpi_seq()  &
      & .OR. lrestore_states) THEN
      DEALLOCATE (p_int_state_lonlat_global, STAT=error_status)
    ENDIF
    DEALLOCATE (p_int_state, STAT=error_status)
    IF (error_status /= SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for ptr_int_state failed')
    ENDIF

    ! Deallocate grid patches

    CALL destruct_patches( p_patch )

    IF (my_process_is_mpi_seq()  &
      & .OR. lrestore_states) THEN
      DEALLOCATE( p_patch_global, STAT=error_status )
    ELSE
      DEALLOCATE( p_patch_subdiv, STAT=error_status )
    ENDIF
    IF (error_status/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocate for patch array failed')
    ENDIF

    ! finalize meteogram output
    DO jg = 1, n_dom
      IF (mtgrm_output_config(jg)%lenabled) THEN
        CALL mtgrm_finalize(mtgrm_output_config(jg), jg)
      END IF
    END DO
    DO jg = 1, max_dom
      DEALLOCATE(mtgrm_output_config(jg)%station_list)
    END DO

    ! clear restart namelist buffer
    CALL delete_restart_namelists()

    CALL message(TRIM(routine),'clean-up finished')

  END SUBROUTINE atmo_model

END MODULE mo_atmo_model

