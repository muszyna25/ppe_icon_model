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
  USE mo_exception,               ONLY: message, finish, message_text
  USE mo_mpi,                     ONLY: stop_mpi, my_process_is_io, my_process_is_mpi_test,   &
    &                                   set_mpi_work_communicators,                           &
    &                                   p_pe_work, process_mpi_io_size,                       &
    &                                   my_process_is_restart, process_mpi_restart_size,      &
    &                                   my_process_is_pref, process_mpi_pref_size  
  USE mo_timer,                   ONLY: init_timer, timer_start, timer_stop,                  &
    &                                   timers_level, timer_model_init,                       &
    &                                   timer_domain_decomp, timer_compute_coeffs,            &
    &                                   timer_ext_data, print_timer
  USE mo_parallel_config,         ONLY: p_test_run, l_test_openmp, num_io_procs,              &
    &                                   num_restart_procs, use_async_restart_output,          &
    &                                   num_prefetch_proc
  USE mo_master_config,           ONLY: isRestart, setStartdate, setStopdate, setExpStopdate, &
       &                                tc_startdate, tc_stopdate, tc_dt_restart,             &
       &                                tc_exp_stopdate, tc_exp_startdate
  USE mo_master_control,          ONLY: get_my_process_name, get_my_model_no
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

  ! coupling
  USE mo_coupling_config,         ONLY: is_coupled_run
  USE mo_interface_echam_ocean,   ONLY: construct_atmo_coupler, destruct_atmo_coupler

  ! I/O
  USE mo_io_restart_async,        ONLY: restart_main_proc                                       ! main procedure for Restart PEs
  USE mo_name_list_output,        ONLY: name_list_io_main_proc
  USE mo_name_list_output_config, ONLY: use_async_name_list_io
  USE mo_io_restart_namelist,     ONLY: delete_restart_namelists
  USE mo_time_config,             ONLY: time_config      ! variable
  USE mo_output_event_types,      ONLY: t_sim_step_info
  USE mtime,                      ONLY: setCalendar, MAX_DATETIME_STR_LEN,         &
       &                                datetime, newDatetime, deallocateDatetime, &
       &                                datetimeToString, OPERATOR(<), OPERATOR(+)
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
    !
    IF ( is_coupled_run() ) THEN
      CALL construct_atmo_coupler(p_patch)
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
    CALL destruct_atmo_model ()

    !---------------------------------------------------------------------
    ! destruct the coupler
    !
    IF ( is_coupled_run() ) THEN
      CALL destruct_atmo_coupler ()
    ENDIF

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
    CHARACTER(*), PARAMETER :: routine = "mo_atmo_model:construct_atmo_model"
    INTEGER                 :: jg, jgp, jstep0, error_status
    TYPE(t_sim_step_info)   :: sim_step_info  

    TYPE(datetime), POINTER :: calculatedStopDate => NULL()
    CHARACTER(len=max_datetime_str_len) :: startDate = '', dstring

    ! initialize global registry of lon-lat grids
    CALL init_lonlat_grid_list()

    !---------------------------------------------------------------------
    ! 0. If this is a resumed or warm-start run...
    !---------------------------------------------------------------------

    IF (isRestart()) THEN
      CALL message('','Read restart file meta data ...')
      CALL read_restart_header("atm")
      CALL get_restart_attribute('tc_startdate', startDate)
    ELSE
      CALL datetimeToString(tc_exp_startdate, startDate)
    ENDIF

    !---------------------------------------------------------------------
    ! 0.1 check for start and stop dates
    !---------------------------------------------------------------------

    CALL setStartdate(startDate)
    
    IF (ASSOCIATED(tc_startdate) .AND. ASSOCIATED(tc_dt_restart)) THEN
      calculatedStopDate => newDatetime('0001-01-01T00:00:00')
      calculatedStopDate = tc_startdate + tc_dt_restart 
      IF (ASSOCIATED(tc_exp_stopdate)) THEN
        IF (tc_exp_stopdate < calculatedStopDate) THEN
          calculatedStopDate = tc_exp_stopdate
          CALL message('','Experiment stop date earlier than run stop date. '// &
           &         'Reset end to experiment stop date!')   
        ENDIF
      ENDIF
      CALL datetimeToString(calculatedStopDate, dstring)
      CALL setStopdate(dstring)
      IF (.NOT. ASSOCIATED(tc_exp_stopdate)) THEN
        CALL setExpStopdate(dstring)
      ENDIF
      CALL deallocateDatetime(calculatedStopDate)
    ELSE
#ifdef USE_MTIME_LOOP
      CALL finish('','Cannot calculate this runs stop date.')
#else
      CALL message('use_mtime_loop','Cannot calculate this runs stop date.')
#endif
    ENDIF

    IF (ASSOCIATED(tc_startdate)) THEN
      CALL datetimeToString(tc_startdate, dstring)
      WRITE(message_text,'(a,a)') 'Start date of this run   : ', dstring
      CALL message('',message_text)
    ENDIF

    IF (ASSOCIATED(tc_stopdate)) THEN
      CALL datetimeToString(tc_stopdate, dstring)
      WRITE(message_text,'(a,a)') 'Stop date of this run    : ', dstring
      CALL message('',message_text)
    ENDIF

    CALL message('','')

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
        CALL prefetch_main_proc(tc_startdate)
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
          CALL datetimeToString(tc_exp_startdate, sim_step_info%sim_start)
          CALL datetimeToString(tc_exp_stopdate, sim_step_info%sim_end)
          CALL datetimeToString(tc_startdate, sim_step_info%run_start)
          CALL datetimeToString(tc_stopdate, sim_step_info%restart_time)
          sim_step_info%dtime      = dtime
          jstep0 = 0
          IF (isRestart() .AND. .NOT. time_config%is_relative_time) THEN
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

END MODULE mo_atmo_model
