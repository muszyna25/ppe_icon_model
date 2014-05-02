!>
!! @brief branch for the non-hydrostatic ICON workflow
!!
!! @author Kristina Froehlich, MPI-M (2011-07-19)
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
MODULE mo_atmo_nonhydrostatic

USE mo_kind,                 ONLY: wp
USE mo_exception,            ONLY: message, finish
USE mo_impl_constants,       ONLY: SUCCESS, max_dom, inwp, iecham
USE mo_mpi,                  ONLY: my_process_is_stdio
USE mo_timer,                ONLY: timers_level, timer_start, timer_stop, &
  &                                timer_model_init, timer_init_icon, timer_read_restart
USE mo_master_control,       ONLY: is_restart_run
USE mo_var_list,             ONLY: print_var_list
USE mo_time_config,          ONLY: time_config      ! variable
USE mo_io_restart,           ONLY: read_restart_files
USE mo_io_restart_attributes,ONLY: get_restart_attribute
USE mo_io_config,            ONLY: dt_diag,dt_checkpoint
USE mo_parallel_config,      ONLY: nproma
USE mo_nh_pzlev_config,      ONLY: configure_nh_pzlev
USE mo_advection_config,     ONLY: configure_advection
USE mo_art_config,           ONLY: configure_art
USE mo_run_config,           ONLY: dtime, dtime_adv,     & !    namelist parameter
  &                                ltestcase,            &
  &                                iforcing,             & !    namelist parameter
  &                                output_mode,          &
  &                                msg_level,            & !    namelist parameter
  &                                lvert_nest, ntracer,  &
  &                                iqc, iqt
USE mo_dynamics_config,      ONLY: iequations
! Horizontal grid
USE mo_model_domain,         ONLY: p_patch
USE mo_grid_config,          ONLY: n_dom, start_time, is_plane_torus
! to break circular dependency KF???
USE mo_intp_data_strc,       ONLY: p_int_state
USE mo_grf_intp_data_strc,   ONLY: p_grf_state
! NH-namelist state
USE mo_nonhydrostatic_config,ONLY: iadv_rcf, kstart_moist, kend_qvsubstep, l_open_ubc

USE mo_atm_phy_nwp_config,   ONLY: configure_atm_phy_nwp, atm_phy_nwp_config
! NH-Model states
USE mo_nonhydro_state,       ONLY: p_nh_state, construct_nh_state, destruct_nh_state
USE mo_opt_diagnostics,      ONLY: construct_opt_diag, destruct_opt_diag
USE mo_nwp_phy_state,        ONLY: prm_diag, prm_nwp_diag_list, prm_nwp_tend_list, &
  &                                construct_nwp_phy_state,                        &
  &                                destruct_nwp_phy_state
USE mo_nwp_lnd_state,        ONLY: p_lnd_state, construct_nwp_lnd_state,       &
  &                                destruct_nwp_lnd_state
! Time integration
USE mo_nh_stepping,          ONLY: prepare_nh_integration, perform_nh_stepping
! Initialization with real data
USE mo_nh_initicon,         ONLY: init_icon
USE mo_ext_data_state,      ONLY: ext_data, init_index_lists
! meteogram output
USE mo_meteogram_output,    ONLY: meteogram_init, meteogram_finalize
USE mo_meteogram_config,    ONLY: meteogram_output_config
USE mo_name_list_output_config,   ONLY: first_output_name_list, &
  &                               is_variable_in_output
USE mo_name_list_output_init, ONLY: init_name_list_output, &
  &                               parse_variable_groups
USE mo_name_list_output,    ONLY: close_name_list_output
USE mo_pp_scheduler,        ONLY: pp_scheduler_init, pp_scheduler_finalize
USE mo_intp_lonlat,         ONLY: compute_lonlat_area_weights
USE mo_mtime_extensions,    ONLY: get_datetime_string
USE mo_output_event_types,  ONLY: t_sim_step_info
USE mo_action,              ONLY: action_init
USE mo_turbulent_diagnostic,ONLY: init_les_turbulent_output, close_les_turbulent_output

!-------------------------------------------------------------------------

IMPLICIT NONE
PRIVATE

PUBLIC :: atmo_nonhydrostatic
PUBLIC :: construct_atmo_nonhydrostatic, destruct_atmo_nonhydrostatic

! module data
INTEGER :: n_diag, n_chkpt

CONTAINS

  !---------------------------------------------------------------------
  SUBROUTINE atmo_nonhydrostatic
    
!!$    CHARACTER(*), PARAMETER :: routine = "mo_atmo_nonhydrostatic"

    CALL construct_atmo_nonhydrostatic()
    
    !------------------------------------------------------------------
    ! Now start the time stepping:
    ! The special initial time step for the three time level schemes
    ! is executed within process_grid_level
    !------------------------------------------------------------------

    CALL perform_nh_stepping( time_config%cur_datetime,       &
      &                       n_chkpt, n_diag  )
 
    !---------------------------------------------------------------------
    ! 6. Integration finished. Clean up.
    !---------------------------------------------------------------------
    CALL destruct_atmo_nonhydrostatic()
   
  END SUBROUTINE atmo_nonhydrostatic
  !---------------------------------------------------------------------
  
  !---------------------------------------------------------------------
  SUBROUTINE construct_atmo_nonhydrostatic
    
    CHARACTER(*), PARAMETER :: routine = "construct_atmo_nonhydrostatic"

    INTEGER :: jg, ist
    LOGICAL :: l_realcase
    LOGICAL :: l_pres_msl(n_dom) !< Flag. TRUE if computation of mean sea level pressure desired
    LOGICAL :: l_omega(n_dom)    !< Flag. TRUE if computation of vertical velocity desired
    LOGICAL :: l_rh(n_dom)       !< Flag. TRUE if computation of relative humidity desired
    TYPE(t_sim_step_info) :: sim_step_info  
    INTEGER :: jstep0


    IF (timers_level > 3) CALL timer_start(timer_model_init)

    IF(iforcing == inwp) THEN
      !
      ! - generate index lists for tiles (land, ocean, lake)
      ! index lists for ice-covered and non-ice covered ocean points 
      ! are initialized in init_nwp_phy
      CALL init_index_lists (p_patch(1:), ext_data)

      CALL configure_atm_phy_nwp(n_dom, p_patch(1:), dtime_adv)

     ! initialize number of chemical tracers for convection 
     DO jg = 1, n_dom
       CALL configure_art(jg)
     ENDDO

    ENDIF

    IF (.NOT. ltestcase .AND. (iforcing == inwp .OR. iforcing == iecham)) THEN
      l_realcase = .TRUE.
    ELSE
      l_realcase = .FALSE.
    ENDIF
 
    ! initialize ldom_active flag
    DO jg=1, n_dom
      IF (jg > 1 .AND. start_time(jg) > 0._wp) THEN
        p_patch(jg)%ldom_active = .FALSE. ! domain not active from the beginning
      ELSE
        p_patch(jg)%ldom_active = .TRUE.
      ENDIF
    ENDDO

    !---------------------------------------------------------------------
    ! 4.c Non-Hydrostatic / NWP
    !---------------------------------------------------------------------

    ALLOCATE (p_nh_state(n_dom), stat=ist)
    IF (ist /= success) THEN
      CALL finish(TRIM(routine),'allocation for p_nh_state failed')
    ENDIF

    ! Note(GZ): Land state now needs to be allocated even if physics is turned
    ! off because ground temperature is included in feedback since r8133
    ! However, setting inwp_surface = 0 effects that only a few 2D fields are allocated
    ALLOCATE (p_lnd_state(n_dom), stat=ist)
    IF (ist /= success) THEN
      CALL finish(TRIM(routine),'allocation for p_lnd_state failed')
    ENDIF

    IF(iforcing /= inwp) atm_phy_nwp_config(:)%inwp_surface = 0



    ! Now allocate memory for the states
    DO jg=1,n_dom
      l_pres_msl(jg) = is_variable_in_output(first_output_name_list, var_name="pres_msl")
      l_omega(jg)    = is_variable_in_output(first_output_name_list, var_name="omega")
    END DO
    CALL construct_nh_state(p_patch(1:), p_nh_state, n_timelevels=2, &
      &                     l_pres_msl=l_pres_msl, l_omega=l_omega)

    ! Add optional diagnostic variable lists (might remain empty)
    CALL construct_opt_diag(p_patch(1:), .TRUE.)

    IF(iforcing == inwp) THEN
      DO jg=1,n_dom
        l_rh(jg) = is_variable_in_output(first_output_name_list, &
          &                              var_name="rh")
      END DO
      CALL construct_nwp_phy_state( p_patch(1:), l_rh )
    ENDIF

    CALL construct_nwp_lnd_state( p_patch(1:),p_lnd_state,n_timelevels=2 )

#ifdef MESSY
    DO jg=1,n_dom
       CALL messy_init_memory(jg)
    END DO
#endif

    ! Due to the required ability to overwrite advection-Namelist settings 
    ! via add_ref/add_tracer_ref for ICON-ART, configure_advection is called 
    ! AFTER the nh_state is created. Otherwise, potential modifications of the 
    ! advection-Namelist can not be taken into account properly.
    ! Unfortunatley this conflicts with our trying to call the config-routines 
    ! as early as possible. 
    DO jg =1,n_dom
     CALL configure_advection( jg, p_patch(jg)%nlev, p_patch(1)%nlev,  &
       &                      iequations, iforcing, iqc, iqt,          &
       &                      kstart_moist(jg), kend_qvsubstep(jg),    &
       &                      lvert_nest, l_open_ubc, ntracer ) 
    ENDDO

#ifdef MESSY
    CALL messy_init_coupling
    CALL messy_init_tracer
#endif

    !---------------------------------------------------------------------
    ! 5. Perform time stepping
    !---------------------------------------------------------------------
      !------------------------------------------------------------------
      ! Prepare for time integration
      !------------------------------------------------------------------


    ! CALL prepare_nh_integration(p_patch(1:), p_nh_state, p_int_state(1:), p_grf_state(1:))

    CALL prepare_nh_integration( )

    !
    ! Read restart files (if necessary)
    !
    IF (is_restart_run()) THEN
      ! This is a resumed integration. Read model state from restart file(s).

      IF (timers_level > 5) CALL timer_start(timer_read_restart)
#ifdef NOMPI
      ! TODO : Non-MPI mode does not work for multiple domains
      DO jg = 1,n_dom
        CALL read_restart_files( p_patch(jg), n_dom)
      END DO
#else
      DO jg = 1,n_dom
        CALL read_restart_files( p_patch(jg), n_dom )
      END DO
#endif      
      CALL message(TRIM(routine),'normal exit from read_restart_files')
      IF (timers_level > 5) CALL timer_stop(timer_read_restart)

    ENDIF


    !---------------------------------------------------------------------
    !     Setup of meteogram output
    !---------------------------------------------------------------------

    IF (output_mode%l_nml) THEN
      DO jg =1,n_dom
        IF (meteogram_output_config(jg)%lenabled) THEN
          ! For dry test cases: do not sample moist variables
          ! (but allow for TORUS moist runs; see also mo_mtgrm_output.F90)
          IF (ltestcase .and. .not. is_plane_torus) THEN
             CALL meteogram_init(meteogram_output_config(jg), jg, p_patch(jg), &
              &                ext_data(jg), p_nh_state(jg),                  &
              &                p_lnd_state=p_lnd_state(jg), iforcing=iforcing)
          ELSE
            CALL meteogram_init(meteogram_output_config(jg), jg, p_patch(jg), &
              &                ext_data(jg), p_nh_state(jg), prm_diag(jg),    &
              &                p_lnd_state(jg), iforcing)
          END IF
        END IF
      END DO
    END IF



    !------------------------------------------------------------------
    ! Prepare initial conditions for time integration.
    !------------------------------------------------------------------
    !
    ! Initialize model with real atmospheric data if appropriate switches are set
    !
    IF (l_realcase .AND. .NOT. is_restart_run()) THEN

      IF (timers_level > 5) CALL timer_start(timer_init_icon)
      CALL init_icon (p_patch(1:), p_nh_state(1:), prm_diag(1:), p_lnd_state(1:), &
        &             p_int_state(1:), p_grf_state(1:), ext_data(1:))
      IF (timers_level > 5) CALL timer_stop(timer_init_icon)

    ENDIF


    !---------------------------------------------------------
    ! The most primitive event handling algorithm: 
    ! compute time step interval for taking a certain action
    !--------------------------------------------------------- 
 
    ! !!OUTDATED!! writing output is now controlled via 'istime4output' !!OUTDATED!!
    n_chkpt = NINT(dt_checkpoint/dtime)  ! write restart files
    n_diag  = MAX(1,NINT(dt_diag/dtime)) ! diagnose of total integrals


    !------------------------------------------------------------------
    ! Prepare output file
    !------------------------------------------------------------------
    !
    ! if output on z and/or p-levels is required do some config
    DO jg = 1, n_dom
      CALL configure_nh_pzlev(jg, nproma, p_patch(jg)%npromz_c,  &
        &                     p_patch(jg)%nblks_c)
    ENDDO


    ! Add a special metrics variable containing the area weights of
    ! the regular lon-lat grid.
    CALL compute_lonlat_area_weights()

    ! Map the variable groups given in the output namelist onto the
    ! corresponding variable subsets:
    ! ATTENTION: all add_vars must be finished before calling this routine.
    IF (output_mode%l_nml) THEN
      CALL parse_variable_groups()
    END IF
   

    ! setup of post-processing job queue, e.g. setup of optional
    ! diagnostic quantities like pz-level interpolation
    CALL pp_scheduler_init( (.NOT. ltestcase .OR. iforcing == inwp) )

    ! If async IO is in effect, init_name_list_output is a collective call
    ! with the IO procs and effectively starts async IO
    IF (output_mode%l_nml) THEN
      ! compute sim_start, sim_end
      CALL get_datetime_string(sim_step_info%sim_start, time_config%ini_datetime)
      CALL get_datetime_string(sim_step_info%sim_end,   time_config%end_datetime)
      CALL get_datetime_string(sim_step_info%restart_time,  time_config%cur_datetime, &
        &                      INT(time_config%dt_restart))
      CALL get_datetime_string(sim_step_info%run_start, time_config%cur_datetime)
      sim_step_info%dtime      = dtime
      sim_step_info%iadv_rcf   = iadv_rcf
      jstep0 = 0
      IF (is_restart_run() .AND. .NOT. time_config%is_relative_time) THEN
        ! get start counter for time loop from restart file:
        CALL get_restart_attribute("jstep", jstep0)
      END IF
      sim_step_info%jstep0    = jstep0
      CALL init_name_list_output(sim_step_info)
    END IF


    !----------------------!
    !  Initialize actions  !
    !----------------------!
    !
    ! assign variables to existing actions
    CALL action_init()


    ! for debug purpose: print var lists
    IF ( msg_level >=20 .AND. my_process_is_stdio() .AND. .NOT. ltestcase) THEN
      CALL print_var_list (p_nh_state(1)%prog_list(1))
      CALL print_var_list (p_nh_state(1)%diag_list)
      CALL print_var_list (p_nh_state(1)%metrics_list)
      CALL print_var_list (prm_nwp_diag_list(1))
      CALL print_var_list (prm_nwp_tend_list(1))
      CALL print_var_list (p_lnd_state(1)%lnd_prog_nwp_list(1))
      CALL print_var_list (p_lnd_state(1)%lnd_diag_nwp_list)
      CALL print_var_list (ext_data(1)%atm_list)
      CALL print_var_list (ext_data(1)%atm_td_list)
    ENDIF

    !Anurag Dipankar, MPIM (2014-01-14)
    !Special 1D and 0D output for LES runs till we get add_var/nml_out working
    IF(atm_phy_nwp_config(1)%is_les_phy .AND. is_restart_run()) &
      CALL init_les_turbulent_output(p_patch(1), p_nh_state(1)%metrics, &
                             time_config%sim_time(1), ldelete=.FALSE.)

    IF(atm_phy_nwp_config(1)%is_les_phy .AND. .NOT.is_restart_run()) &
      CALL init_les_turbulent_output(p_patch(1), p_nh_state(1)%metrics, &
                             time_config%sim_time(1), ldelete=.TRUE.)


    IF (timers_level > 3) CALL timer_stop(timer_model_init)

  END SUBROUTINE construct_atmo_nonhydrostatic
  
  !---------------------------------------------------------------------
  SUBROUTINE destruct_atmo_nonhydrostatic
    
    CHARACTER(*), PARAMETER :: routine = "destruct_atmo_nonhydrostatic"


    INTEGER :: jg, ist

    !---------------------------------------------------------------------
    ! 6. Integration finished. Clean up.
    !---------------------------------------------------------------------

    CALL message(TRIM(routine),'start to clean up')

#ifdef MESSY
    CALL messy_free_memory
#endif

    ! Destruction of post-processing job queue
    CALL pp_scheduler_finalize()

    ! Delete optional diagnostics
    CALL destruct_opt_diag()
   
    ! Delete state variables

    CALL destruct_nh_state( p_nh_state )
    DEALLOCATE (p_nh_state, STAT=ist)
    IF (ist /= SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for p_nh_state failed')
    ENDIF

    IF (iforcing == inwp) THEN
      CALL destruct_nwp_phy_state
    ENDIF
    CALL destruct_nwp_lnd_state( p_lnd_state )

    ! Delete output variable lists
    IF (output_mode%l_nml) THEN
      CALL close_name_list_output
    END IF

    ! finalize meteogram output
    IF (output_mode%l_nml) THEN
      DO jg = 1, n_dom
        IF (meteogram_output_config(jg)%lenabled) THEN
          CALL meteogram_finalize(jg)
        END IF
      END DO
      DO jg = 1, max_dom
        DEALLOCATE(meteogram_output_config(jg)%station_list)
      END DO
    END IF

    !Close LES diag files
    DO jg = 1 , n_dom
     IF(atm_phy_nwp_config(jg)%is_les_phy) &
       CALL close_les_turbulent_output(jg)
    END DO

    CALL message(TRIM(routine),'clean-up finished')
    
  END SUBROUTINE destruct_atmo_nonhydrostatic
  !---------------------------------------------------------------------

END MODULE mo_atmo_nonhydrostatic

