!>
!! @brief workflow for the ICON atmospheric hydrostatic model
!!
!! @author Hui Wan (MPI-M)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_atmo_hydrostatic

  USE mo_exception,         ONLY: message
  USE mo_impl_constants,    ONLY: iecham, ildf_echam

  USE mo_master_config,     ONLY: isRestart
  USE mo_time_config,       ONLY: time_config
  USE mo_run_config,        ONLY: dtime, iforcing, nlev, &
    &                             msg_level, output_mode, ntracer, iqv, iqc, iqt
  USE mo_dynamics_config,   ONLY: iequations, nnow, idiv_method
  USE mo_advection_config,  ONLY: configure_advection
  USE mo_ha_testcases,      ONLY: ctest_name
  USE mo_io_config,         ONLY: configure_io
  USE mo_grid_config,       ONLY: n_dom

  USE mo_model_domain,        ONLY: p_patch
  USE mo_intp_data_strc,      ONLY: p_int_state
  USE mo_grf_intp_data_strc,  ONLY: p_grf_state

  USE mo_vertical_coord_table,ONLY: vct_a, vct_b
  USE mo_icoham_dyn_memory,   ONLY: p_hydro_state, destruct_icoham_dyn_state
  USE mo_ha_stepping,         ONLY: prepare_ha_dyn, initcond_ha_dyn, &
                                    perform_ha_stepping
  USE mo_ha_dyn_config,       ONLY: ha_dyn_config
  USE mo_echam_phy_init,      ONLY: init_echam_phy, initcond_echam_phy !, additional_restart_init
  USE mo_echam_phy_cleanup,   ONLY: cleanup_echam_phy

  USE mo_io_restart,           ONLY: read_restart_files
  USE mo_io_restart_attributes,ONLY: get_restart_attribute
  USE mo_name_list_output_init, ONLY: init_name_list_output
  USE mo_name_list_output,     ONLY:  write_name_list_output, &
       &                              close_name_list_output
  USE mo_mtime_extensions,     ONLY: get_datetime_string
  USE mo_output_event_types,   ONLY: t_sim_step_info
  USE mtime,                   ONLY: setCalendar, PROLEPTIC_GREGORIAN


  IMPLICIT NONE
  PRIVATE

  PUBLIC :: atmo_hydrostatic
  PUBLIC :: construct_atmo_hydrostatic, destruct_atmo_hydrostatic



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
                            & p_hydro_state, time_config%cur_datetime )

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
        &                      1, 0, .FALSE., .FALSE., ntracer,         &
        &                      idiv_method, ha_dyn_config%itime_scheme ) 
    ENDDO

    IF (iforcing==IECHAM.OR.iforcing==ILDF_ECHAM) THEN
      CALL init_echam_phy( p_patch(1:), ctest_name, &
                            & nlev, vct_a, vct_b, time_config%cur_datetime )
    END IF

    !------------------------------------------------------------------
    ! Set initial conditions for time integration.
    !------------------------------------------------------------------
    ! Here constant values are also initialized,
    ! such as analytical topography.
    ! It should be called even in restart mode
    CALL initcond_ha_dyn( p_patch(1:), p_int_state(1:),  &
                        & p_grf_state(1:), p_hydro_state )

    IF (iforcing==IECHAM.OR.iforcing==ILDF_ECHAM) THEN
      DO jg = 1,n_dom
        CALL initcond_echam_phy( jg                                                  ,&
          &                      p_patch(jg)                                         ,&
          &                      p_hydro_state(jg)%prog(nnow(jg))% temp  (:,:,:)     ,&
          &                      p_hydro_state(jg)%prog(nnow(jg))% tracer(:,:,:,iqv) ,&
          &                      ctest_name                                          )
      END DO
    END IF

    !------------------------------------------------------------------
    ! The most primitive event handling algorithm:
    ! compute time step interval for taking a certain action
    !------------------------------------------------------------------

    CALL configure_io()   ! set n_chkpt and n_diag, which control 
                          ! writing of restart files and tot_int diagnostics.


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
      jstep0 = 0
      IF (isRestart() .AND. .NOT. time_config%is_relative_time) THEN
        ! get start counter for time loop from restart file:
        CALL get_restart_attribute("jstep", jstep0)
      END IF
      sim_step_info%jstep0    = jstep0
      CALL init_name_list_output(sim_step_info)
    ENDIF

    IF (.NOT.isRestart()) THEN
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

    END IF ! (not) isRestart()

    !------------------------------------------------------------------
    IF (isRestart()) THEN
    ! This is an resumed integration. Read model state from restart file(s).

      jg = 1
#ifdef NOMPI
      CALL read_restart_files( p_patch(jg) )
#else
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
      !  CALL additional_restart_init( p_patch(1:), ctest_name )
      ! END IF

    END IF ! isRestart()
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

