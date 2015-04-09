!>
!! @brief Main program for the ICON atmospheric model
!!
!! @author
!!  Leonidas Linardakis (MPI-M)
!!  Hui Wan             (MPI-M)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_atmo_model

  ! basic modules
  USE mo_exception,               ONLY: message, finish !, message_text
  USE mo_mpi,                     ONLY: stop_mpi, my_process_is_io, my_process_is_mpi_test,   &
    &                                   set_mpi_work_communicators,                           &
    &                                   p_pe_work, process_mpi_io_size,                       &
    &                                   my_process_is_restart, process_mpi_restart_size,      &
    &                                   my_process_is_pref, process_mpi_pref_size  
  USE mo_timer,                   ONLY: init_timer, timer_start, timer_stop,                  &
    &                                   timers_level, timer_model_init,                       &
    &                                   timer_domain_decomp, timer_compute_coeffs,            &
    &                                   timer_ext_data, print_timer
  USE mo_parallel_config,         ONLY: p_test_run, l_test_openmp, num_io_procs,               &
    &                                   num_restart_procs, use_async_restart_output,          &
    &                                   num_prefetch_proc
  USE mo_master_control,          ONLY: is_restart_run, get_my_process_name, get_my_model_no
#ifndef NOMPI
#if defined(__GET_MAXRSS__)
  USE mo_mpi,                     ONLY: get_my_mpi_all_id
  USE mo_util_sysinfo,            ONLY: util_get_maxrss
#endif
#endif
  USE mo_impl_constants,          ONLY: SUCCESS, MAX_CHAR_LENGTH,                             &
    &                                   ihs_atm_temp, ihs_atm_theta, inh_atmosphere,          &
    &                                   ishallow_water, inwp
  USE mo_io_restart,              ONLY: read_restart_header
  USE mo_io_restart_attributes,   ONLY: get_restart_attribute

  ! namelist handling; control parameters: run control, dynamics
  USE mo_read_namelists,          ONLY: read_atmo_namelists
  USE mo_nml_crosscheck,          ONLY: atm_crosscheck
  USE mo_nonhydrostatic_config,   ONLY: configure_nonhydrostatic
  USE mo_initicon_config,         ONLY: configure_initicon
  USE mo_lnd_nwp_config,          ONLY: configure_lnd_nwp
  USE mo_dynamics_config,         ONLY: configure_dynamics, iequations
  USE mo_run_config,              ONLY: configure_run,                                        &
    &                                   ltimer, ltestcase,                                    &
    &                                   nshift,                                               &
    &                                   num_lev,                                              &
    &                                   msg_level,                                            &
    &                                   dtime, output_mode,                                   &
    &                                   grid_generatingCenter,                                & ! grid generating center
    &                                   grid_generatingSubcenter,                             & ! grid generating subcenter
    &                                   iforcing
  USE mo_gribout_config,          ONLY: configure_gribout

  ! time stepping
  USE mo_atmo_hydrostatic,        ONLY: atmo_hydrostatic
  USE mo_atmo_nonhydrostatic,     ONLY: atmo_nonhydrostatic

  USE mo_nh_testcases,            ONLY: init_nh_testtopo

  ! coupling
# ifdef YAC_coupling
  USE mo_kind,                    ONLY: wp
  USE mo_mpi,                     ONLY: p_n_work
  USE mo_math_constants,          ONLY: pi
  USE mo_loopindices,             ONLY: get_indices_c
  USE mo_parallel_config,         ONLY: nproma
  USE mo_grid_subset,             ONLY: t_subset_range, get_index_range
  USE mo_yac_finterface,          ONLY: yac_finit, yac_fdef_comp,                    &
    &                                   yac_fdef_datetime,                           &
    &                                   yac_fdef_subdomain, yac_fconnect_subdomains, &
    &                                   yac_fdef_elements, yac_fdef_points,          &
    &                                   yac_fdef_mask, yac_fdef_field, yac_fsearch,  &
    &                                   yac_ffinalize, yac_redirstdout

  USE mo_coupling_config,         ONLY: is_coupled_run
  USE mo_model_domain,            ONLY: t_patch
# else
  USE mo_icon_cpl_init,           ONLY: icon_cpl_init
  USE mo_icon_cpl_init_comp,      ONLY: icon_cpl_init_comp
  USE mo_coupling_config,         ONLY: is_coupled_run, config_debug_coupler_level
  USE mo_icon_cpl_def_grid,       ONLY: icon_cpl_def_grid, icon_cpl_def_location
  USE mo_icon_cpl_def_field,      ONLY: icon_cpl_def_field
  USE mo_icon_cpl_search,         ONLY: icon_cpl_search
  USE mo_icon_cpl_finalize,       ONLY: icon_cpl_finalize
#endif

  USE mo_alloc_patches,           ONLY: destruct_patches

  ! horizontal grid, domain decomposition, memory
  USE mo_grid_config,             ONLY: n_dom, n_dom_start,                 &
    &                                   dynamics_parent_grid_id, n_phys_dom
  USE mo_model_domain,            ONLY: p_patch, p_patch_local_parent
  USE mo_build_decomposition,     ONLY: build_decomposition
  USE mo_complete_subdivision,    ONLY: setup_phys_patches
  USE mo_icon_comm_interface,     ONLY: construct_icon_communication,                         &
    &                                   destruct_icon_communication
  ! Vertical grid
  USE mo_vertical_coord_table,    ONLY: apzero, vct_a, vct_b, vct, allocate_vct_atmo
  USE mo_nh_init_utils,           ONLY: nflatlev
  USE mo_util_vgrid,              ONLY: construct_vertical_grid

  ! external data, physics
  USE mo_ext_data_state,          ONLY: ext_data, init_ext_data, destruct_ext_data
  USE mo_rrtm_data_interface,     ONLY: init_rrtm_model_repart, destruct_rrtm_model_repart
  USE mo_nwp_ww,                  ONLY: configure_ww

  USE mo_diffusion_config,        ONLY: configure_diffusion

  ! horizontal interpolation
  USE mo_interpol_config,         ONLY: configure_interpolation
  USE mo_intp_state,              ONLY: construct_2d_interpol_state,                          &
    &                                   destruct_2d_interpol_state, transfer_interpol_state
  USE mo_grf_intp_state,          ONLY: construct_2d_gridref_state,                           &
    &                                   destruct_2d_gridref_state, transfer_grf_state,        &
    &                                   create_grf_index_lists
  USE mo_intp_data_strc,          ONLY: p_int_state, p_int_state_local_parent
  USE mo_grf_intp_data_strc,      ONLY: p_grf_state, p_grf_state_local_parent
  USE mo_intp_lonlat,             ONLY: init_lonlat_grid_list, compute_lonlat_intp_coeffs,    &
    &                                   destroy_lonlat_grid_list

  ! I/O
  USE mo_io_restart_async,        ONLY: restart_main_proc                                       ! main procedure for Restart PEs
  USE mo_name_list_output,        ONLY: name_list_io_main_proc
  USE mo_name_list_output_config, ONLY: use_async_name_list_io
  USE mo_io_restart_namelist,     ONLY: delete_restart_namelists
  USE mo_time_config,             ONLY: time_config      ! variable
  USE mo_mtime_extensions,        ONLY: get_datetime_string
  USE mo_output_event_types,      ONLY: t_sim_step_info
  USE mtime,                      ONLY: setCalendar,          &
# ifdef YAC_coupling
       &                                MAX_DATETIME_STR_LEN, &
#endif
       &                                PROLEPTIC_GREGORIAN

  ! Prefetching  
  USE mo_async_latbc,             ONLY: prefetch_main_proc
  ! ART
  USE mo_art_init_interface,      ONLY: art_init_interface

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
      CALL finish( TRIM(routine),'unknown choice for iequations.')
    END SELECT

    ! print performance timers:
    IF (ltimer) CALL print_timer


    !---------------------------------------------------------------------
    ! 13. Integration finished. Carry out the shared clean-up processes
    !---------------------------------------------------------------------
    CALL destruct_atmo_model()

    !---------------------------------------------------------------------
    ! destruct the coupler
# ifdef YAC_coupling
    IF ( is_coupled_run() ) CALL yac_ffinalize
# else
    IF ( is_coupled_run() ) CALL icon_cpl_finalize ()
# endif

    !---------------------------------------------------------------------
    ! (optional:) write resident set size from OS
#ifndef NOMPI
#if defined(__GET_MAXRSS__)
    IF (msg_level >= 16) THEN
      CALL util_get_maxrss(maxrss)
      PRINT  *, "PE #", get_my_mpi_all_id(), &
        &    ": MAXRSS (MiB) = ", maxrss
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
    ! local variables
    CHARACTER(*), PARAMETER :: routine = "mo_atmo_model:init_atmo_model"
    INTEGER                 :: jg, jgp, jstep0, error_status
    TYPE(t_sim_step_info)   :: sim_step_info  

    ! set mtime-Calendar
    CALL setCalendar(PROLEPTIC_GREGORIAN)

    ! initialize global registry of lon-lat grids
    CALL init_lonlat_grid_list()

    !---------------------------------------------------------------------
    ! 0. If this is a resumed or warm-start run...
    !---------------------------------------------------------------------

    IF (is_restart_run()) THEN
      CALL read_restart_header("atm")
    END IF ! is_restart_run()

    !---------------------------------------------------------------------
    ! 1.1 Read namelists (newly) specified by the user; fill the
    !     corresponding sections of the configuration states.
    !---------------------------------------------------------------------

    CALL read_atmo_namelists(atm_namelist_filename,shr_namelist_filename)

    !---------------------------------------------------------------------
    ! 1.2 Cross-check namelist setups
    !---------------------------------------------------------------------

    !! DR end temporary hack !!
    CALL atm_crosscheck

#ifdef MESSY
    CALL messy_setup
#endif

    !---------------------------------------------------------------------
    ! 2. Call configure_run to finish filling the run_config state.
    !    This needs to be done very early (but anyway after atm_crosscheck)
    !    because some components of the state, e.g., num_lev, may be
    !    modified in this subroutine which affects the following CALLs.
    !---------------------------------------------------------------------
    CALL configure_run( )


    ! complete initicon config-state
    CALL configure_initicon


    !-------------------------------------------------------------------
    ! 3.1 Initialize the mpi work groups
    !-------------------------------------------------------------------
    CALL set_mpi_work_communicators(p_test_run, l_test_openmp, num_io_procs, num_restart_procs, &
               &                    num_prefetch_proc)

    !-------------------------------------------------------------------
    ! 3.2 Initialize various timers
    !-------------------------------------------------------------------
    IF (ltimer) CALL init_timer
    IF (timers_level > 3) CALL timer_start(timer_model_init)

    !-------------------------------------------------------------------
    ! 3.3 I/O initialization
    !-------------------------------------------------------------------

    ! If we belong to the Restart PEs just call restart_main_proc before reading patches.
    ! This routine will never return
    IF (process_mpi_restart_size > 0) THEN
      use_async_restart_output = .TRUE.
      CALL message('','asynchronous restart output is enabled.')
      IF (my_process_is_restart()) THEN
        CALL restart_main_proc
      ENDIF
    ENDIF

    ! If we belong to the prefetching PEs just call prefetch_main_proc before reading patches.
    ! This routine will never return
    IF (process_mpi_pref_size > 0) THEN
      num_prefetch_proc = 1
      CALL message(routine,'asynchronous input prefetching is enabled.')
      IF (my_process_is_pref() .AND. (.NOT. my_process_is_mpi_test())) THEN
        CALL prefetch_main_proc  
      ENDIF
    ENDIF
 
    ! If we belong to the I/O PEs just call xxx_io_main_proc before
    ! reading patches.  This routine will never return
    IF (process_mpi_io_size > 0) THEN
      ! Decide whether async vlist or name_list IO is to be used,
      ! only one of both may be enabled!

      IF (output_mode%l_nml) THEN
        ! -----------------------------------------
        ! asynchronous I/O
        ! -----------------------------------------
        !
        use_async_name_list_io = .TRUE.
        CALL message(routine,'asynchronous namelist I/O scheme is enabled.')
        ! consistency check
        IF (my_process_is_io() .AND. (.NOT. my_process_is_mpi_test())) THEN
          ! Stop timer which is already started but would not be stopped
          ! since xxx_io_main_proc never returns
          IF (timers_level > 3) CALL timer_stop(timer_model_init)

          ! compute sim_start, sim_end
          CALL get_datetime_string(sim_step_info%sim_start, time_config%ini_datetime)
          CALL get_datetime_string(sim_step_info%sim_end,   time_config%end_datetime)
          CALL get_datetime_string(sim_step_info%restart_time,  time_config%cur_datetime, &
            &                      INT(time_config%dt_restart))
          CALL get_datetime_string(sim_step_info%run_start, time_config%cur_datetime)
          sim_step_info%dtime      = dtime
          jstep0 = 0
          IF (is_restart_run() .AND. .NOT. time_config%is_relative_time) THEN
            ! get start counter for time loop from restart file:
            CALL get_restart_attribute("jstep", jstep0)
          END IF
          sim_step_info%jstep0    = jstep0
          CALL name_list_io_main_proc(sim_step_info)
        END IF
      ELSE IF (my_process_is_io() .AND. (.NOT. my_process_is_mpi_test())) THEN
        ! Shut down MPI
        CALL stop_mpi
        STOP
      ENDIF
    ELSE
      ! -----------------------------------------
      ! non-asynchronous I/O (performed by PE #0)
      ! -----------------------------------------
      !
      IF (output_mode%l_nml) THEN
        CALL message(routine,'synchronous namelist I/O scheme is enabled.')
      ENDIF
    ENDIF

    !------------------
    ! Next, define the horizontal and vertical grids since they are aready
    ! needed for some derived control parameters. This includes
    ! - patch import
    ! - domain decompistion
    ! - vertical coordinates
    !-------------------------------------------------------------------
    ! 4. Import patches, perform domain decomposition
    !-------------------------------------------------------------------

    IF (timers_level > 5) CALL timer_start(timer_domain_decomp)
    CALL build_decomposition(num_lev, nshift, is_ocean_decomposition = .false.)
    IF (timers_level > 5) CALL timer_stop(timer_domain_decomp)

    !--------------------------------------------------------------------------------
    ! 5. Construct interpolation state, compute interpolation coefficients.
    !--------------------------------------------------------------------------------

    IF (timers_level > 5) CALL timer_start(timer_compute_coeffs)
    CALL configure_interpolation( n_dom, p_patch(1:)%level, &
                                  p_patch(1)%geometry_info )

    ! Allocate array for interpolation state

    ALLOCATE( p_int_state(n_dom_start:n_dom), &
            & p_grf_state(n_dom_start:n_dom), &
            & STAT=error_status)
    IF (error_status /= SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for ptr_int_state failed')
    ENDIF

    ALLOCATE( p_int_state_local_parent(n_dom_start+1:n_dom), &
      &       p_grf_state_local_parent(n_dom_start+1:n_dom), &
      &       STAT=error_status)
    IF (error_status /= SUCCESS) &
      CALL finish(TRIM(routine),'allocation for local parents failed')

    ! Construct interpolation state
    ! Please note that for parallel runs the divided state is constructed here
    CALL construct_2d_interpol_state(p_patch, p_int_state)

    ! Transfer interpolation state to local parent
    DO jg = n_dom_start+1, n_dom
      jgp = p_patch(jg)%parent_id
      CALL transfer_interpol_state(p_patch(jgp),p_patch_local_parent(jg), &
        &  p_int_state(jgp), p_int_state_local_parent(jg))
    ENDDO

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

      ! Transfer gridref state from local parent to p_grf_state
      DO jg = n_dom_start+1, n_dom
        jgp = p_patch(jg)%parent_id
        CALL transfer_grf_state(p_patch(jgp), p_patch_local_parent(jg),         &
          &                     p_grf_state(jgp), p_grf_state_local_parent(jg), &
          &                     p_patch(jg)%parent_child_index)
      ENDDO
    ENDIF


    !-------------------------------------------------------------------
    ! Initialize icon_comm_lib
    !-------------------------------------------------------------------
!    IF (use_icon_comm) THEN
      CALL construct_icon_communication(p_patch, n_dom)
!    ENDIF


    !--------------------------------------------
    ! Setup the information for the physical patches
    CALL setup_phys_patches

    !-------------------------------------------------------------------
    ! 7. Constructing data for lon-lat interpolation
    !-------------------------------------------------------------------
    
    CALL compute_lonlat_intp_coeffs(p_patch(1:), p_int_state(1:))

    IF (n_dom_start==0 .OR. n_dom > 1) THEN
      CALL create_grf_index_lists(p_patch, p_grf_state, p_int_state)
    ENDIF
    IF (timers_level > 5) CALL timer_stop(timer_compute_coeffs)

   !---------------------------------------------------------------------
   ! Prepare dynamics and land 
   !---------------------------------------------------------------------

    CALL configure_dynamics ( n_dom )

    IF (iforcing == inwp) THEN ! set dimensions of tile-based variables
      CALL configure_lnd_nwp()
    ENDIF

    !------------------------------------------------------------------
    ! Create and optionally read external data fields
    !------------------------------------------------------------------
    ALLOCATE (ext_data(n_dom), STAT=error_status)
    IF (error_status /= SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for ext_data failed')
    ENDIF

    ! allocate memory for atmospheric/oceanic external data and
    ! optionally read those data from netCDF file.
    IF (timers_level > 5) CALL timer_start(timer_ext_data)
    CALL init_ext_data (p_patch(1:), p_int_state(1:), ext_data)
    IF (timers_level > 5) CALL timer_stop(timer_ext_data)

   !---------------------------------------------------------------------
   ! Import vertical grid/ define vertical coordinate
   !---------------------------------------------------------------------

    CALL allocate_vct_atmo(p_patch(1)%nlevp1)
    IF (iequations == inh_atmosphere .AND. ltestcase) THEN
      CALL init_nh_testtopo(p_patch(1:), ext_data)   ! set analytic topography
    ENDIF
    CALL construct_vertical_grid(p_patch(1:), p_int_state(1:), ext_data, &
      &                          vct_a, vct_b, vct, nflatlev)


    !---------------------------------------------------------------------
    ! Horizontal and vertical grid(s) are now defined.
    ! Assign values to derived variables in the configuration states
    !---------------------------------------------------------------------


    CALL configure_diffusion( n_dom, dynamics_parent_grid_id,       &
      &                       p_patch(1)%nlev, vct_a, vct_b, apzero )

    CALL configure_gribout(grid_generatingCenter, grid_generatingSubcenter, n_phys_dom)

    IF (iequations == inh_atmosphere) THEN
      DO jg =1,n_dom
        CALL configure_nonhydrostatic( jg, p_patch(jg)%nlev,     &
          &                            p_patch(jg)%nshift_total  )
        IF ( iforcing == inwp) THEN
          CALL configure_ww( jg, p_patch(jg)%nlev, p_patch(jg)%nshift_total)
        END IF
      ENDDO
    ENDIF



    !------------------------------------------------------------------
    ! 11. Repartitioning of radiation grid (Karteileiche?!)
    !------------------------------------------------------------------
    CALL init_rrtm_model_repart()

#ifdef MESSY
    CALL messy_initialize(n_dom)
    CALL messy_new_tracer
#endif

    !------------------------------------------------------------------
    ! 11. Create ART data fields
    !------------------------------------------------------------------
    
    CALL art_init_interface(n_dom,'construct')

    !------------------------------------------------------------------

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
    CALL destruct_patches( p_patch_local_parent )
    IF (msg_level > 5) CALL message(TRIM(routine),'destruct_patches is done')

    DEALLOCATE( p_patch, STAT=error_status )
    IF (error_status/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocate for patch array failed')
    ENDIF

    ! clear restart namelist buffer
    CALL delete_restart_namelists()
    IF (msg_level > 5) CALL message(TRIM(routine),'delete_restart_namelists is done')

    CALL destruct_rrtm_model_repart()
!    IF (use_icon_comm) THEN
      CALL destruct_icon_communication()
!    ENDIF

    ! Destruct ART data fields
    CALL art_init_interface(n_dom,'destruct')

    CALL message(TRIM(routine),'clean-up finished')

  END SUBROUTINE destruct_atmo_model
  !-------------------------------------------------------------------

  !-------------------------------------------------------------------
  !>
  SUBROUTINE construct_atmo_coupler()
    ! For the coupling

    INTEGER, PARAMETER :: no_of_fields = 10

    CHARACTER(LEN=MAX_CHAR_LENGTH) ::  field_name(no_of_fields)
    INTEGER :: field_id(no_of_fields)
    INTEGER :: i, patch_no
    INTEGER :: error_status

    !---------------------------------------------------------------------
    ! 11. Do the setup for the coupled run
    !
    ! For the time being this could all go into a subroutine which is
    ! common to atmo and ocean. Does this make sense if the setup deviates
    ! too much in future.
    !---------------------------------------------------------------------

# ifdef YAC_coupling

    INTEGER, PARAMETER :: nbr_subdomain_ids = 1
    INTEGER, PARAMETER :: CELL = 0 ! one point per cell
    ! (see definition of enum location in points.h)

    REAL(wp), PARAMETER :: deg = 180.0_wp / pi

    CHARACTER(LEN=max_char_length) :: xml_filename
    CHARACTER(LEN=max_char_length) :: xsd_filename
    CHARACTER(LEN=max_char_length) :: grid_name
    CHARACTER(LEN=max_char_length) :: comp_name

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: iso8601_ref_datetime ! ISO_8601

    INTEGER :: comp_id
    INTEGER :: comp_ids(1)
    INTEGER :: cell_point_ids(1)
    INTEGER :: cell_mask_ids(1)
    INTEGER :: domain_id
    INTEGER :: subdomain_id
    INTEGER :: subdomain_ids(nbr_subdomain_ids)

    INTEGER :: mask_checksum
    INTEGER :: i_startidx, i_endidx
    INTEGER :: nblks
    INTEGER :: BLOCK, idx, INDEX
    INTEGER :: nbr_vertices_per_cell

    REAL(wp), ALLOCATABLE :: buffer_lon(:)
    REAL(wp), ALLOCATABLE :: buffer_lat(:)
    INTEGER,  ALLOCATABLE :: buffer_c(:,:)
    INTEGER,  ALLOCATABLE :: ibuffer(:)

    TYPE(t_sim_step_info)   :: sim_step_info  
    TYPE(t_subset_range), POINTER :: cells_in_domain
    TYPE(t_patch), POINTER :: patch_horz

    IF ( .NOT. is_coupled_run() ) RETURN

    comp_name = TRIM(get_my_process_name())

    patch_no = 1

    i = LEN_TRIM(comp_name)
    CALL yac_redirstdout ( TRIM(comp_name), i, 1, p_pe_work, p_n_work, error_status )

    ! Initialise the coupler
    xml_filename = "coupling.xml"
    xsd_filename = "coupling.xsd"
    CALL yac_finit ( xml_filename, xsd_filename )

    ! Inform the coupler about what we are
    CALL yac_fdef_comp ( comp_name, comp_id )
    comp_ids(1) = comp_id

    ! Overwrite job start and end date with component data
    CALL get_datetime_string(sim_step_info%sim_start, time_config%ini_datetime)
    CALL get_datetime_string(sim_step_info%sim_end,   time_config%end_datetime)

    CALL yac_fdef_datetime ( start_datetime = TRIM(sim_step_info%sim_start),  &
      &                      end_datetime   = TRIM(sim_step_info%sim_end)   )

    ! Announce one subdomain (patch) to the coupler
    grid_name = "grid1"
    CALL yac_fdef_subdomain ( comp_id, grid_name, subdomain_id )

    subdomain_ids(1) = subdomain_id

    patch_horz => p_patch(patch_no)
    cells_in_domain  => patch_horz%cells%in_domain

    ! Extract cell information
    !
    ! cartesian coordinates of cell vertices are stored in
    ! patch_horz%verts%cartesian(:,:)%x(1:3)
    ! Here we use the longitudes and latitudes.

    nblks = max(patch_horz%nblks_c,patch_horz%nblks_v)

    ALLOCATE(buffer_lon(nproma*nblks))
    ALLOCATE(buffer_lat(nproma*nblks))
    ALLOCATE(buffer_c(3,nproma*nblks))

    nbr_vertices_per_cell = 3

    DO BLOCK = 1, patch_horz%nblks_v
      DO idx = 1, nproma
        INDEX = (BLOCK-1)*nproma+idx
        buffer_lon(INDEX) = patch_horz%verts%vertex(idx,BLOCK)%lon * deg
        buffer_lat(INDEX) = patch_horz%verts%vertex(idx,BLOCK)%lat * deg
      ENDDO
    ENDDO

    DO BLOCK = 1, patch_horz%nblks_c
      DO idx = 1, nproma
        INDEX = (BLOCK-1)*nproma+idx
        buffer_c(1,INDEX) = (patch_horz%cells%vertex_blk(idx,BLOCK,1)-1)*nproma + &
          &                  patch_horz%cells%vertex_idx(idx,BLOCK,1)
        buffer_c(2,INDEX) = (patch_horz%cells%vertex_blk(idx,BLOCK,2)-1)*nproma + &
          &                  patch_horz%cells%vertex_idx(idx,BLOCK,2)
        buffer_c(3,INDEX) = (patch_horz%cells%vertex_blk(idx,BLOCK,3)-1)*nproma + &
                             patch_horz%cells%vertex_idx(idx,BLOCK,3)
      ENDDO
    ENDDO

    ! Description of elements, here as unstructured grid
    CALL yac_fdef_elements (      &
      & subdomain_id,             &
      & patch_horz%n_patch_verts, &
      & patch_horz%n_patch_cells, &
      & nbr_vertices_per_cell,    &
      & buffer_lon,               &
      & buffer_lat,               &
      & buffer_c )

    ! Can we have two fdef_point calls for the same subdomain, i.e.
    ! one single set of cells?
    !
    ! Define cell center points (location = 0)
    !
    ! cartesian coordinates of cell centers are stored in
    ! patch_horz%cells%cartesian_center(:,:)%x(1:3)
    ! Here we use the longitudes and latitudes.

    DO BLOCK = 1, patch_horz%nblks_c
      DO idx = 1, nproma
        INDEX = (BLOCK-1)*nproma+idx
        buffer_lon(INDEX) = patch_horz%cells%center(idx,BLOCK)%lon * deg
        buffer_lat(INDEX) = patch_horz%cells%center(idx,BLOCK)%lat * deg
      ENDDO
    ENDDO

    ! center points in cells (needed e.g. for patch recovery and nearest neighbour)
    CALL yac_fdef_points (        &
      & subdomain_id,             &
      & patch_horz%n_patch_cells, &
      & CELL,                     &
      & buffer_lon,               &
      & buffer_lat,               &
      & cell_point_ids(1) )

    DEALLOCATE (buffer_lon, buffer_lat, buffer_c)

    ALLOCATE(ibuffer(nproma*patch_horz%nblks_c))

    DO idx = 1, patch_horz%n_patch_cells
       IF ( p_pe_work == patch_horz%cells%decomp_info%owner_local(idx) ) THEN
         ibuffer(idx) = -1
       ELSE
         ibuffer(idx) = patch_horz%cells%decomp_info%owner_local(idx)
       ENDIF
    ENDDO

    ! decomposition information
    CALL yac_fdef_index_location (              &
      & subdomain_id,                           &
      & patch_horz%n_patch_cells,               &
      & CELL,                                   &
      & patch_horz%cells%decomp_info%glb_index, &
      & ibuffer )

    ! Connect subdomains
    CALL yac_fconnect_subdomains ( &
      & comp_id,                   &
      & nbr_subdomain_ids,         &
      & subdomain_ids,             &
      & domain_id )

    !
    ! The land-sea mask for the ocean is available in patch_3D%surface_cell_sea_land_mask(:,:)
    !
    !          -2: inner ocean
    !          -1: boundary ocean
    !           1: boundary land
    !           2: inner land
    !
    ! The mask which is used in the ECHAM physics is prm_field(1)%lsmask(:,:).
    ! This locial mask is set to .TRUE. for land points.
    ! We can get access to via "USE mo_echam_phy_memory,ONLY: prm_field"
    !
    ! Here we use a mask which is hopefully identical to the one used by the
    ! ocean, and which works independent of the physics chosen. 
    !
    mask_checksum = 0

!rr    DO BLOCK = 1, patch_horz%nblks_c
!rr       CALL get_indices_c ( patch_horz, BLOCK, 1, patch_horz%nblks_c,  &
!rr                               i_startidx, i_endidx, 2 )
!rr       DO idx = i_startidx, i_endidx
!rr          mask_checksum = mask_checksum + ABS(ext_data(1)%atm%lsm_ctr_c(idx,BLOCK))
!rr       ENDDO
!rr    ENDDO

    DO i = 1, patch_horz%n_patch_cells
       ibuffer(i) = 0
    ENDDO

    IF ( mask_checksum > 0 ) THEN
       DO BLOCK = 1, patch_horz%nblks_c
          CALL get_indices_c ( patch_horz, BLOCK, 1, patch_horz%nblks_c,  &
                               i_startidx, i_endidx, 2 )
          DO idx = i_startidx, i_endidx
             IF ( ext_data(1)%atm%lsm_ctr_c(idx,BLOCK) < 0 ) THEN
                ibuffer((BLOCK-1)*nproma+idx) = 0
             ELSE
                ibuffer((BLOCK-1)*nproma+idx) = 1
             ENDIF
          ENDDO
       ENDDO
    ELSE
       DO i = 1, patch_horz%n_patch_cells
          ibuffer(i) = 0
       ENDDO
    ENDIF

    CALL yac_fdef_mask (           &
      & patch_horz%n_patch_cells,  &
      & ibuffer,                   &
      & cell_point_ids(1),         &
      & cell_mask_ids(1) )
    
    DEALLOCATE (ibuffer)

    field_name(1) = "surface_downward_eastward_stress"   ! bundled field containing two components
    field_name(2) = "surface_downward_northward_stress"  ! bundled field containing two components
    field_name(3) = "surface_fresh_water_flux"           ! bundled field containing three components
    field_name(4) = "surface_temperature"
    field_name(5) = "total_heat_flux"                    ! bundled field containing four components
    field_name(6) = "atmosphere_sea_ice_bundle"          ! bundled field containing four components
    field_name(7) = "sea_surface_temperature"
    field_name(8) = "eastward_sea_water_velocity"
    field_name(9) = "northward_sea_water_velocity"
    field_name(10) = "ocean_sea_ice_bundle"              ! bundled field containing five components

    DO i = 1, no_of_fields
      CALL yac_fdef_field ( &
        & field_name(i),    &
        & comp_id,          &
        & domain_id,        &
        & cell_point_ids,   &
        & cell_mask_ids,    &
        & 1,                &
        & field_id(i) )
    ENDDO

    CALL yac_fsearch ( 1, comp_ids, no_of_fields, field_id, error_status )

# else

    INTEGER :: grid_id
    INTEGER :: grid_shape(2)
    INTEGER :: field_shape(3)

    IF ( .NOT. is_coupled_run() ) RETURN

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
    ! should grid_glob_index become a pointer in icon_cpl_def_grid as well?
    CALL icon_cpl_def_grid ( &
      & grid_shape, p_patch(patch_no)%cells%decomp_info%glb_index, & ! input
      & grid_id, error_status )                                      ! output

    ! Marker for internal and halo points, a list which contains the
    ! rank where the native cells are located.
    CALL icon_cpl_def_location ( &
      & grid_id, grid_shape, p_patch(patch_no)%cells%decomp_info%owner_local, & ! input
      & p_pe_work,  &                                                           ! this owner id
      & error_status )                                                          ! output

    field_name(1) = "TAUX"   ! bundled field containing two components
    field_name(2) = "TAUY"   ! bundled field containing two components
    field_name(3) = "SFWFLX" ! bundled field containing three components
    field_name(4) = "SFTEMP"
    field_name(5) = "THFLX"  ! bundled field containing four components
    field_name(6) = "ICEATM" ! bundled field containing four components
    field_name(7) = "SST"
    field_name(8) = "OCEANU"
    field_name(9) = "OCEANV"
    field_name(10) = "ICEOCE" ! bundled field containing five components

    field_shape(1:2) = grid_shape(1:2)

    ! see equivalent atmosphere counterpart in ocean/boundary/mo_ocean_coupling.f90
    ! routine construct_ocean_coupling

    DO i = 1, no_of_fields

       IF ( i == 1 .OR. i == 2 ) THEN
         field_shape(3) = 2
       ELSE IF ( i == 3 ) THEN
         field_shape(3) = 3
       ELSE IF ( i == 5 .OR. i == 6 ) THEN
         field_shape(3) = 4
       ELSE IF ( i == 10 ) THEN
         field_shape(3) = 5
       ELSE
         field_shape(3) = 1
       ENDIF

       CALL icon_cpl_def_field ( &
         & field_name(i), grid_id, field_id(i), &
         & field_shape, error_status )

    ENDDO

    CALL icon_cpl_search

#endif

  END SUBROUTINE construct_atmo_coupler

END MODULE mo_atmo_model
