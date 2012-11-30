! ----------------------------------------------------------------------
! Main subroutine for the "old" prepicon tool
! Handles operation modes for
! - generating coordinate fields
! - converting IFS2ICON output
! - doing vertical interpolation
!
! Moved from the former drivers/prep_icon: F. Prill, DWD (2012-11-13)
! ----------------------------------------------------------------------
!
MODULE mo_prepicon_old

  USE mo_exception,             ONLY: message, finish

  USE mo_parallel_config,       ONLY: p_test_run, nproma
  USE mo_mpi,                   ONLY: my_process_is_mpi_seq,     &
    &                                 my_process_is_mpi_parallel
  USE mo_timer,                 ONLY: init_timer


  ! Control parameters: run control, dynamics, i/o
  !
  USE mo_nonhydrostatic_config, ONLY: ivctype
  USE mo_dynamics_config,       ONLY: iequations
  USE mo_run_config,            ONLY: ltimer,               &
    &                                 nlev,nlevp1,          &
    &                                 num_lev,num_levp1, nshift

  USE mo_impl_constants, ONLY:  inh_atmosphere

  ! Memory
  !
  USE mo_setup_subdivision,     ONLY: decompose_domain
  USE mo_complete_subdivision,  ONLY:                            &
    &                                 complete_parallel_setup,   &
    &                                 copy_processor_splitting,  &
    &                                 set_patch_communicators,   &
    &                                 finalize_decomposition


  USE mo_model_domain,        ONLY: t_patch, p_patch, p_patch_local_parent

  ! Horizontal grid
  !
  USE mo_grid_config,           ONLY: n_dom, n_dom_start, global_cell_type

  USE mo_model_domimp_patches,  ONLY: import_basic_patches, complete_patches

  USE mo_alloc_patches,  ONLY: destruct_patches
  ! Horizontal interpolation
  !
  USE mo_intp_state,            ONLY: construct_2d_interpol_state, &
    & destruct_2d_interpol_state, transfer_interpol_state

  USE mo_grf_intp_state,     ONLY: construct_2d_gridref_state,  &
    & destruct_2d_gridref_state, transfer_grf_state


  ! Vertical grid
  !
  USE mo_nh_init_utils,         ONLY: init_hybrid_coord, init_sleve_coord


  USE mo_impl_constants,        ONLY: SUCCESS


  USE mo_intp_data_strc,        ONLY: p_int_state, p_int_state_local_parent
  USE mo_grf_intp_data_strc,    ONLY: p_grf_state, p_grf_state_local_parent

  USE mo_time_config,           ONLY: time_config         ! variable
  USE mo_dynamics_config,       ONLY: configure_dynamics  ! subroutine
  USE mo_interpol_config
  USE mo_lnd_nwp_config,        ONLY: configure_lnd_nwp
  USE mo_ext_data_state,        ONLY: ext_data, init_ext_data, destruct_ext_data

  ! USE statements referring directly to prep_icon
  !
  USE mo_time_config,           ONLY: time_config 
  USE mo_prepicon_utils,        ONLY: init_prepicon, prepicon, write_prepicon_output, &
    & compute_coord_fields, init_topo_output_files, close_prepicon_output_files,      &
    & convert_variables, init_atmo_output_files, deallocate_prepicon

  USE mo_prepicon_config,       ONLY: i_oper_mode, l_zp_out
  USE mo_nh_vert_interp,        ONLY: vertical_interpolation,                         &
    &                                 intp_to_p_and_z_levels_prepicon

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: prepicon_main

CONTAINS

  ! Main subroutine for the "old" prepicon tool
  ! Handles operation modes for
  ! - generating coordinate fields
  ! - converting IFS2ICON output
  ! - doing vertical interpolation
  !
  SUBROUTINE prepicon_main()

    ! local variables
    CHARACTER(*), PARAMETER :: routine = "prep_icon"
    
    TYPE(t_patch), ALLOCATABLE :: p_patch_global(:)
    
    INTEGER :: ist, jg, jgp
    INTEGER :: error_status

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

    ! Check patch allocation status

    IF ( ALLOCATED(p_patch)) THEN
      CALL finish(TRIM(routine), 'patch already allocated')
    END IF

    ! Allocate patch array to start patch construction

    ALLOCATE(p_patch(n_dom_start:n_dom), stat=error_status)
    IF (error_status/=success) THEN
      CALL finish(TRIM(routine), 'allocation of patch failed')
    ENDIF

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


    !-----------------------------------------------------------------------------
    ! 6. Construct grid refinment state, compute coefficients
    !-----------------------------------------------------------------------------
    ! For the NH model, the initialization routines called from
    ! construct_2d_gridref_state require the metric terms to be present

    IF (n_dom_start==0 .OR. n_dom > 1) THEN
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

    !-------------------------------------------------------------------
    ! 7. Finalize domain decomposition
    !-------------------------------------------------------------------   
    IF (my_process_is_mpi_parallel()) THEN

      CALL finalize_decomposition()

    ENDIF



    !---------------------------------------------------------------------
    ! 8. Import vertical grid/ define vertical coordinate
    !---------------------------------------------------------------------
    SELECT CASE (iequations)

    CASE (inh_atmosphere)
      IF (ivctype == 1) THEN
        CALL init_hybrid_coord(p_patch(1)%nlev)
      ELSE IF (ivctype == 2) THEN
        CALL init_sleve_coord(p_patch(1)%nlev)
      ENDIF
    CASE DEFAULT
      CALL finish(TRIM(routine),'prep_icon only works for NH model')
    END SELECT


    !---------------------------------------------------------------------
    ! 9. Horizontal and vertical grid(s) are now defined.
    !    Assign values to derived variables in the configuration states
    !---------------------------------------------------------------------

    CALL configure_dynamics ( n_dom )

    CALL configure_lnd_nwp( p_patch(1:), n_dom, nproma )

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


    ! Allocate prepicon data type
    ALLOCATE (prepicon(n_dom), stat=ist)
    IF (ist /= success) THEN
      CALL finish(TRIM(routine),'allocation for prepicon failed')
    ENDIF

    ! Set all model domains to active at initial time
    DO jg = 1, n_dom
      p_patch(jg)%ldom_active = .TRUE.
    ENDDO

    ! allocate memory for topography and coordinate fields,
    ! read topo data from netCDF file, 
    ! optionally smooth topography data,
    ! and, in case of nesting, perform topography blending and feedback
    CALL init_prepicon (p_int_state(1:), p_grf_state(1:), prepicon, ext_data)

    ! Compute the 3D coordinate fields
    CALL compute_coord_fields(p_int_state(1:), prepicon)

    IF (i_oper_mode == 2) THEN

      ! Convert (hydrostatic) IFS2ICON output variables to the NH set of prognostic variables
      CALL convert_variables(p_int_state(1:), prepicon)

    ELSE IF (i_oper_mode == 3) THEN
      ! Perform vertical interpolation from intermediate IFS2ICON grid to ICON grid
      ! and convert variables to the NH set of prognostic variables
      CALL vertical_interpolation(p_patch(1:), p_int_state(1:), p_grf_state(1:), prepicon)
    ENDIF

    IF (i_oper_mode >= 2 .AND. l_zp_out) THEN
      ! Interpolate prognostic variables to pressure and height levels for diagnostic output
      CALL intp_to_p_and_z_levels_prepicon(p_patch(1:), prepicon)
    ENDIF

    !------------------------------------------------------------------
    ! Write output
    !------------------------------------------------------------------

    IF (i_oper_mode == 1) THEN

      CALL init_topo_output_files
      CALL write_prepicon_output( time_config%ini_datetime )
      CALL close_prepicon_output_files

    ELSE

      CALL init_atmo_output_files
      CALL write_prepicon_output( time_config%ini_datetime )
      CALL close_prepicon_output_files

    ENDIF

    ! Deallocate prepicon data type
    CALL deallocate_prepicon(prepicon)
    DEALLOCATE (prepicon, stat=ist)
    IF (ist /= success) THEN
      CALL finish(TRIM(routine),'deallocation for prepicon failed')
    ENDIF


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
    DEALLOCATE (p_int_state, STAT=error_status)
    IF (error_status /= SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for ptr_int_state failed')
    ENDIF

    ! Deallocate grid patches

    CALL destruct_patches( p_patch )

    DEALLOCATE( p_patch, STAT=error_status )
    IF (error_status/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocate for patch array failed')
    ENDIF

    CALL message(TRIM(routine),'clean-up finished')

  END SUBROUTINE prepicon_main
  
END MODULE mo_prepicon_old
