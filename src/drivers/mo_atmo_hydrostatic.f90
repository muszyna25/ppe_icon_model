!>
!! @brief workflow for the ICON atmospheric hydrostatic model
!!
!! @author Hui Wan (MPI-M)
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
MODULE mo_atmo_hydrostatic

  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message
  USE mo_impl_constants,    ONLY: iecham, ildf_echam
  USE mo_timer,             ONLY: print_timer

  USE mo_master_control,    ONLY: is_restart_run
  USE mo_time_config,       ONLY: time_config
  USE mo_run_config,        ONLY: dtime, nsteps, ltestcase, ltimer, iforcing, nlev, &
    &                             msg_level, output_mode, ntracer, iqc, iqt
  USE mo_dynamics_config,   ONLY: iequations
  USE mo_advection_config,  ONLY: configure_advection
  USE mo_ha_testcases,      ONLY: ctest_name
  USE mo_io_config,         ONLY: n_diags, n_checkpoints
  USE mo_grid_config,       ONLY: n_dom

  USE mo_model_domain,        ONLY: p_patch
  USE mo_intp_data_strc,      ONLY: p_int_state
  USE mo_grf_intp_data_strc,  ONLY: p_grf_state

  USE mo_vertical_coord_table,ONLY: vct_a, vct_b, ceta
  USE mo_icoham_dyn_memory,   ONLY: p_hydro_state, destruct_icoham_dyn_state
  USE mo_ha_stepping,         ONLY: prepare_ha_dyn, initcond_ha_dyn, &
                                    perform_ha_stepping

  USE mo_echam_phy_init,      ONLY: init_echam_phy, initcond_echam_phy, &
                                    additional_restart_init
  USE mo_echam_phy_cleanup,   ONLY: cleanup_echam_phy

  USE mo_io_restart,           ONLY: read_restart_files
  USE mo_io_restart_attributes,ONLY: get_restart_attribute
  USE mo_name_list_output_init, ONLY: init_name_list_output
  USE mo_name_list_output,     ONLY:  write_name_list_output, &
       &                              close_name_list_output
  USE mo_parallel_config,      ONLY: use_icon_comm
  USE mo_mtime_extensions,     ONLY: get_datetime_string
  USE mo_output_event_types,   ONLY: t_sim_step_info
  USE mtime,                   ONLY: setCalendar, PROLEPTIC_GREGORIAN


  IMPLICIT NONE
  PRIVATE

  PUBLIC :: atmo_hydrostatic
  PUBLIC :: construct_atmo_hydrostatic, destruct_atmo_hydrostatic


  INTEGER :: n_diag, n_chkpt


CONTAINS
  !>
  !!
  !!
  SUBROUTINE atmo_hydrostatic

!!$    CHARACTER(*), PARAMETER :: method_name = "atmo_hydrostatic"

    !------------------------------------------------------------------
    ! Initialize parameters and solvers;
    ! Allocate memory for model state vectors.
    !------------------------------------------------------------------

    CALL construct_atmo_hydrostatic()

    !------------------------------------------------------------------
    ! Time integraion
    !------------------------------------------------------------------

    CALL perform_ha_stepping( p_patch(1:), p_int_state(1:),                  &
                            & p_grf_state(1:),                               &
                            & p_hydro_state, time_config%cur_datetime,       &
                            & n_chkpt, n_diag )

    IF (ltimer) CALL print_timer

    !---------------------------------------------------------------------
    ! Integration finished. Start to clean up.
    !---------------------------------------------------------------------
    CALL destruct_atmo_hydrostatic()

  END SUBROUTINE atmo_hydrostatic
  !-------------------------------------------------------------------------------


  !-------------------------------------------------------------------------------
  !>
  SUBROUTINE construct_atmo_hydrostatic

    INTEGER :: jg

    CHARACTER(*), PARAMETER :: method_name = "construct_atmo_hydrostatic"
    TYPE(t_sim_step_info) :: sim_step_info  
    INTEGER :: jstep0

    !------------------------------------------------------------------
    ! Initialize parameters and solvers;
    ! Allocate memory for model state vectors.
    !------------------------------------------------------------------

    CALL prepare_ha_dyn( p_patch(1:) )


    ! Due to the required ability to overwrite advection-Namelist settings 
    ! via add_ref/add_tracer_ref for ICON-ART, configure_advection is called 
    ! AFTER the ha_state is craeted. Otherwise, potential modifications of the 
    ! advection-Namelist can not be taken into account properly.
    ! Unfortunatley this conflicts with our trying to call the config-routines 
    ! as early as possible. Input variables which are usually provided by 
    ! nonhydrostatic_nml are set to default values. 
    DO jg =1,n_dom
     CALL configure_advection( jg, p_patch(jg)%nlev, p_patch(1)%nlev,  &
       &                      iequations, iforcing, iqc, iqt,          &
       &                      1, 0, .FALSE., .FALSE., ntracer          ) 
    ENDDO


    IF (iforcing==IECHAM.OR.iforcing==ILDF_ECHAM) THEN
      CALL init_echam_phy( p_patch(1:), ltestcase, ctest_name, &
                            & nlev, vct_a, vct_b, ceta, time_config%cur_datetime )
    END IF

    !------------------------------------------------------------------
    ! Set initial conditions for time integration.
    !------------------------------------------------------------------
    ! Here constant values are also initialized,
    ! such as analytical topography.
    ! It should be caleed even in reastart mode
    CALL initcond_ha_dyn( p_patch(1:), p_int_state(1:),  &
                        & p_grf_state(1:), p_hydro_state )

    IF (iforcing==IECHAM.OR.iforcing==ILDF_ECHAM) &
      CALL initcond_echam_phy( p_patch(1:),p_hydro_state, ltestcase, ctest_name )


    !------------------------------------------------------------------
    ! The most primitive event handling algorithm:
    ! compute time step interval for taking a certain action
    !------------------------------------------------------------------

    n_chkpt = n_checkpoints()       ! number of: write restart files
    n_diag  = n_diags() ! number of: diagnose of total integrals

    !------------------------------------------------------------------
    ! Initialize output file if necessary;
    ! Write out initial conditions.
    !------------------------------------------------------------------

    IF (output_mode%l_nml) THEN
      CALL setCalendar(PROLEPTIC_GREGORIAN)
      ! compute sim_start, sim_end
      CALL get_datetime_string(sim_step_info%sim_start, time_config%ini_datetime)
      CALL get_datetime_string(sim_step_info%sim_end,   time_config%end_datetime)
      CALL get_datetime_string(sim_step_info%restart_time,  time_config%cur_datetime, &
        &                      INT(time_config%dt_restart))
      CALL get_datetime_string(sim_step_info%run_start, time_config%cur_datetime)
      sim_step_info%dtime      = dtime
      sim_step_info%iadv_rcf   = 1
      jstep0 = 0
      IF (is_restart_run() .AND. .NOT. time_config%is_relative_time) THEN
        ! get start counter for time loop from restart file:
        CALL get_restart_attribute("jstep", jstep0)
      END IF
      sim_step_info%jstep0    = jstep0
      CALL init_name_list_output(sim_step_info)
    ENDIF

    IF (.NOT.is_restart_run()) THEN
    ! Initialize the first output file which will contain also the
    ! initial conditions.

      IF (output_mode%l_nml) THEN
          ! Mis-use optional parameter last_step=.TRUE to force output of initial state.
          ! If one wants e.g. one day of 6-hrly data per output file (i.e. 4 steps) ending at 00:00,
          ! one has to set the first value of output_bounds namelist variable to 06:00, but then the
          ! initial state is now written because sim_time=0._wp is lower than start time of output_bounds.
          ! This should be replaced by an explicit handling of "first_step" similar to "last_step" in
          ! mo_name_list_output, similar to lwrite_inital for vlist output.
        CALL write_name_list_output(jstep=0)
      ENDIF

    END IF ! (not) is_restart_run()

    !------------------------------------------------------------------
    IF (is_restart_run()) THEN
    ! This is an resumed integration. Read model state from restart file(s).

#ifdef NOMPI
      CALL read_restart_files( p_patch(jg) )
#else
      jg = 1
     !DO jg = n_dom_start,n_dom
        CALL read_restart_files( p_patch(jg) )
     !END DO
#endif
      CALL message(TRIM(method_name),'normal exit from read_restart_files')

      ! Re-initialize SST, sea ice and glacier for certain experiments; 
      ! Initialize logical variables in echam physics state.
      ! The latter is necessary for now because logical arrays can not yet be
      ! written into restart files.

      ! LL: initialization has already been done
      ! IF (iforcing==IECHAM.OR.iforcing==ILDF_ECHAM) THEN
      !  CALL additional_restart_init( p_patch(1:), ltestcase, ctest_name )
      ! END IF

    END IF ! is_restart_run()
  END SUBROUTINE construct_atmo_hydrostatic
  !-------------------------------------------------------------------------------


  !-------------------------------------------------------------------------------
  !>
  SUBROUTINE destruct_atmo_hydrostatic()

    CHARACTER(*), PARAMETER :: method_name = "destruct_atmo_hydrostatic"
    !---------------------------------------------------------------------
    ! Integration finished. Start to clean up.
    !---------------------------------------------------------------------
    CALL message(TRIM(method_name),'start to clean up')


    CALL destruct_icoham_dyn_state
    
    IF (iforcing==IECHAM .OR. iforcing==ILDF_ECHAM) THEN
      CALL cleanup_echam_phy
    ENDIF

    IF (msg_level > 5) CALL message(TRIM(method_name),'echam_phy clean up is done')
    
    IF (output_mode%l_nml) THEN    
      CALL close_name_list_output
    ENDIF

  END SUBROUTINE destruct_atmo_hydrostatic


END MODULE mo_atmo_hydrostatic

