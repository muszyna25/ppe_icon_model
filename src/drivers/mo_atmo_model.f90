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
  USE mo_exception,               ONLY: message, finish
  USE mo_mpi,                     ONLY: stop_mpi, my_process_is_io,   &
    &                                   set_mpi_work_communicators, process_mpi_io_size,      &
    &                                   my_process_is_pref, process_mpi_pref_size
#ifdef HAVE_CDI_PIO
  USE mo_mpi,                     ONLY: mpi_comm_null, p_comm_work_io
#endif
  USE mo_timer,                   ONLY: init_timer, timer_start, timer_stop,                  &
    &                                   timers_level, timer_model_init,                       &
    &                                   timer_domain_decomp, timer_compute_coeffs,            &
    &                                   timer_ext_data, print_timer
  USE mo_parallel_config,         ONLY: p_test_run, num_test_pe, l_test_openmp,               &
    &                                   num_io_procs,                                         &
    &                                   num_prefetch_proc, pio_type
  USE mo_master_config,           ONLY: isRestart
  USE mo_memory_log,              ONLY: memory_log_terminate
  USE mo_impl_constants,          ONLY: pio_type_async, pio_type_cdipio
#ifdef HAVE_CDI_PIO
  USE yaxt,                       ONLY: xt_initialize, xt_initialized
  USE mo_cdi,                     ONLY: namespacegetactive
  USE mo_cdi_pio_interface,       ONLY: nml_io_cdi_pio_namespace, &
    &                                   cdi_base_namespace, &
    &                                   nml_io_cdi_pio_client_comm, &
    &                                   nml_io_cdi_pio_conf_handle
#endif
#ifndef NOMPI
#if defined(__GET_MAXRSS__)
  USE mo_mpi,                     ONLY: get_my_mpi_all_id
  USE mo_util_sysinfo,            ONLY: util_get_maxrss
#endif
#endif
  USE mo_impl_constants,          ONLY: SUCCESS,                                              &
    &                                   ihs_atm_temp, ihs_atm_theta, inh_atmosphere,          &
    &                                   ishallow_water, inwp
  USE mo_zaxis_type,              ONLY: zaxisTypeList, t_zaxisTypeList
  USE mo_load_restart,            ONLY: read_restart_header
  USE mo_restart_attributes,      ONLY: t_RestartAttributeList, getAttributesForRestarting

  ! namelist handling; control parameters: run control, dynamics
  USE mo_read_namelists,          ONLY: read_atmo_namelists
  USE mo_nml_crosscheck,          ONLY: atm_crosscheck
  USE mo_nonhydrostatic_config,   ONLY: configure_nonhydrostatic
  USE mo_initicon_config,         ONLY: configure_initicon
  USE mo_io_config,               ONLY: restartWritingParameters
  USE mo_lnd_nwp_config,          ONLY: configure_lnd_nwp, tile_list
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
#ifndef __NO_JSBACH__
  USE mo_echam_phy_config,        ONLY: echam_phy_config
  USE mo_master_control,          ONLY: master_namelist_filename
  USE mo_jsb_base,                ONLY: jsbach_setup => jsbach_setup_models, jsbach_setup_tiles
  USE mo_jsb_model_init,          ONLY: jsbach_setup_grid
  USE mo_jsb_model_final,         ONLY: jsbach_finalize
#endif
  USE mo_master_control,          ONLY: atmo_process

  ! time stepping
  USE mo_atmo_hydrostatic,        ONLY: atmo_hydrostatic
  USE mo_atmo_nonhydrostatic,     ONLY: atmo_nonhydrostatic

  USE mo_nh_testcases,            ONLY: init_nh_testtopo

  USE mo_alloc_patches,           ONLY: destruct_patches, destruct_comm_patterns

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
  USE mo_init_vgrid,              ONLY: nflatlev
  USE mo_util_vgrid,              ONLY: construct_vertical_grid

  ! external data, physics
  USE mo_ext_data_state,          ONLY: ext_data, destruct_ext_data
  USE mo_ext_data_init,           ONLY: init_ext_data
  USE mo_nwp_ww,                  ONLY: configure_ww

  USE mo_diffusion_config,        ONLY: configure_diffusion

  ! horizontal interpolation
  USE mo_interpol_config,         ONLY: configure_interpolation
  USE mo_intp_state,              ONLY: construct_2d_interpol_state,                          &
    &                                   destruct_2d_interpol_state, transfer_interpol_state
  USE mo_grf_intp_state,          ONLY: construct_2d_gridref_state,                           &
    &                                   destruct_2d_gridref_state, transfer_grf_state,        &
    &                                   create_grf_index_lists, destruct_interpol_patterns
  USE mo_intp_data_strc,          ONLY: p_int_state, p_int_state_local_parent
  USE mo_intp_lonlat_types,       ONLY: lonlat_grids
  USE mo_grf_intp_data_strc,      ONLY: p_grf_state, p_grf_state_local_parent
  USE mo_intp_lonlat,             ONLY: compute_lonlat_intp_coeffs

  ! coupling
  USE mo_coupling_config,         ONLY: is_coupled_run
  USE mo_interface_echam_ocean,   ONLY: construct_atmo_coupler, destruct_atmo_coupler

  ! I/O
  USE mo_restart,                 ONLY: detachRestartProcs
  USE mo_name_list_output,        ONLY: name_list_io_main_proc
#ifdef HAVE_CDI_PIO
  USE mo_name_list_output_init,   ONLY: init_cdipio_cb
  USE mo_name_list_output,        ONLY: write_ready_files_cdipio
#endif
  USE mo_name_list_output_config, ONLY: use_async_name_list_io
  USE mo_time_config,             ONLY: time_config      ! variable
  USE mo_output_event_types,      ONLY: t_sim_step_info
  USE mtime,                      ONLY: datetimeToString, OPERATOR(<), OPERATOR(+)
  ! Prefetching  
  USE mo_async_latbc,             ONLY: prefetch_main_proc
  ! ART
  USE mo_art_init_interface,      ONLY: art_init_interface

  !-------------------------------------------------------------------------

  IMPLICIT NONE
  PRIVATE
#ifdef HAVE_CDI_PIO
  INCLUDE 'cdipio.inc'
#endif

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
    INTEGER                 :: jg, jgp, jstep0, error_status, dedicatedRestartProcs
    TYPE(t_sim_step_info)   :: sim_step_info  
    TYPE(t_RestartAttributeList), POINTER :: restartAttributes

    ! initialize global registry of lon-lat grids
    CALL lonlat_grids%init()

    !---------------------------------------------------------------------
    ! 0. If this is a resumed or warm-start run...
    !---------------------------------------------------------------------

    restartAttributes => NULL()
    IF (isRestart()) THEN
      CALL message('','Read restart file meta data ...')
      CALL read_restart_header("atm")
    ENDIF

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
    CALL configure_initicon(dtime)


    !-------------------------------------------------------------------
    ! 3.1 Initialize the mpi work groups
    !-------------------------------------------------------------------
    CALL restartWritingParameters(opt_dedicatedProcCount = dedicatedRestartProcs)
    CALL set_mpi_work_communicators(p_test_run, l_test_openmp, &
         &                          num_io_procs, dedicatedRestartProcs, &
         &                          num_prefetch_proc, num_test_pe,      &
         &                          pio_type, opt_comp_id=atmo_process)

    !-------------------------------------------------------------------
    ! 3.2 Initialize various timers
    !-------------------------------------------------------------------
    IF (ltimer) CALL init_timer
    IF (timers_level > 3) CALL timer_start(timer_model_init)

    !-------------------------------------------------------------------
    ! initialize dynamic list of vertical axes
    !-------------------------------------------------------------------

    zaxisTypeList = t_zaxisTypeList()

#ifndef __NO_JSBACH__
    ! Setup JSBACH: read namelists, configure models for each domain
    ! This has to be after (!) the ICON zaxes have been created in the above line but
    ! before (!) the restart PEs are detached a few lines below since JSBACH
    ! adds its zaxes to zaxisTypeList
    IF (ANY(echam_phy_config(:)%ljsb)) THEN
      ! Do basic initialization of JSBACH
      CALL jsbach_setup(master_namelist_filename)
    END IF
#endif


    !-------------------------------------------------------------------
    ! 3.3 I/O initialization
    !-------------------------------------------------------------------

    ! This won't RETURN on dedicated restart PEs, starting their main loop instead.
    CALL detachRestartProcs()

    ! If we belong to the prefetching PEs just call prefetch_main_proc before reading patches.
    ! This routine will never return
    IF (process_mpi_pref_size > 0) THEN
      num_prefetch_proc = 1
      CALL message(routine,'asynchronous input prefetching is enabled.')
      IF (my_process_is_pref()) CALL prefetch_main_proc
    ENDIF

    ! If we belong to the I/O PEs just call xxx_io_main_proc before
    ! reading patches.  This routine will never return
    IF (process_mpi_io_size > 0 .AND. pio_type == pio_type_async) THEN
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
        IF (my_process_is_io()) THEN
          ! Stop timer which is already started but would not be stopped
          ! since xxx_io_main_proc never returns
          IF (timers_level > 3) CALL timer_stop(timer_model_init)

          ! compute sim_start, sim_end
          CALL datetimeToString(time_config%tc_exp_startdate, sim_step_info%sim_start)
          CALL datetimeToString(time_config%tc_exp_stopdate, sim_step_info%sim_end)
          CALL datetimeToString(time_config%tc_startdate, sim_step_info%run_start)
          CALL datetimeToString(time_config%tc_stopdate, sim_step_info%restart_time)
          sim_step_info%dtime      = dtime
          jstep0 = 0

          restartAttributes => getAttributesForRestarting()
          IF (ASSOCIATED(restartAttributes)) THEN
            ! get start counter for time loop from restart file:
            jstep0 = restartAttributes%getInteger("jstep")
          END IF
          sim_step_info%jstep0    = jstep0
          CALL name_list_io_main_proc(sim_step_info)
        END IF
      ELSE IF (my_process_is_io()) THEN
        ! Shut down MPI
        CALL stop_mpi
        STOP
      ENDIF
    ELSE IF (process_mpi_io_size > 0 .AND. pio_type == pio_type_cdipio) THEN
      ! initialize parallel output via CDI-PIO
#ifdef HAVE_CDI_PIO
      IF (.NOT. xt_initialized()) CALL xt_initialize(p_comm_work_io)
      cdi_base_namespace = namespaceGetActive()
      CALL cdiPioConfSetCallBackActions(nml_io_cdi_pio_conf_handle, &
        cdipio_callback_postcommsetup, init_cdipio_cb)
      CALL cdiPioConfSetCallBackActions(nml_io_cdi_pio_conf_handle, &
        cdipio_callback_postwritebatch, write_ready_files_cdipio)
      nml_io_cdi_pio_client_comm = &
        &   cdiPioInit(p_comm_work_io, nml_io_cdi_pio_conf_handle, &
        &              nml_io_cdi_pio_namespace)
      IF (nml_io_cdi_pio_client_comm == mpi_comm_null) THEN
        ! todo: terminate program cleanly here
        CALL stop_mpi
        STOP
      END IF
#else
      CALL finish(routine, 'CDI-PIO requested but unavailable')
#endif
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
          CALL configure_ww( time_config%tc_startdate, jg, p_patch(jg)%nlev, p_patch(jg)%nshift_total, 'ICON')
        END IF
      ENDDO
    ENDIF

#ifndef __NO_JSBACH__
    ! Setup horizontal grids and tiles for JSBACH
    DO jg=1,n_dom
      IF (echam_phy_config(jg)%ljsb) THEN 
        CALL jsbach_setup_grid( jg, p_patch(jg)) !< in
        CALL jsbach_setup_tiles(jg)
      END IF
    END DO
#endif

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

    ! destruct surface tile list
    IF (iforcing == inwp) THEN
      CALL tile_list%destruct()
    ENDIF

    ! destruct interpolation patterns generate in create_grf_index_lists
    IF (n_dom_start==0 .OR. n_dom > 1) THEN
      CALL destruct_interpol_patterns(p_patch)
    END IF

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
    CALL lonlat_grids%finalize()

    ! Destruct communication patterns
    CALL destruct_comm_patterns( p_patch, p_patch_local_parent )

    ! Deallocate grid patches
    CALL destruct_patches( p_patch )
    CALL destruct_patches( p_patch_local_parent )
    IF (msg_level > 5) CALL message(TRIM(routine),'destruct_patches is done')

    DEALLOCATE( p_patch, STAT=error_status )
    IF (error_status/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocate for patch array failed')
    ENDIF

    ! close memory logging files
    CALL memory_log_terminate

!    IF (use_icon_comm) THEN
      CALL destruct_icon_communication()
!    ENDIF

    ! Destruct ART data fields
    CALL art_init_interface(n_dom,'destruct')

#ifndef __NO_JSBACH__
    CALL jsbach_finalize()
#endif
    CALL message(TRIM(routine),'clean-up finished')

  END SUBROUTINE destruct_atmo_model
  !-------------------------------------------------------------------

END MODULE mo_atmo_model
