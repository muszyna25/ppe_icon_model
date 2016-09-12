!> This module contains routines to set start date, end date etc. 
!!
!! These routines are necessary, since (for reasons of backward
!! compatibility) different naemlist settings may be used to specify
!! this date information.
!!
!! @author F. Prill, DWD
!!
!! @par Revision History
!! Initial implementation  by  F. Prill, DWD (2015-10-06)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!! -----------------------------------------------------------------------------------
MODULE mo_time_management

  USE, INTRINSIC :: iso_c_binding, ONLY: c_int32_t
  USE mo_kind,                     ONLY: wp, i8
  USE mo_parallel_config,          ONLY: num_restart_procs
  USE mo_util_string,              ONLY: tolower, int2string
  USE mo_nonhydrostatic_config,    ONLY: divdamp_order
  USE mtime,                       ONLY: MAX_DATETIME_STR_LEN, datetime,                   &
    &                                    MAX_CALENDAR_STR_LEN,                             &
    &                                    MAX_TIMEDELTA_STR_LEN,                            &
    &                                    OPERATOR(>),OPERATOR(/=), OPERATOR(-),            &
    &                                    newDatetime, deallocateDatetime, timedelta,       &
    &                                    getPTStringFromMS, newTimedelta, min,             &
    &                                    deallocateTimedelta, OPERATOR(+), OPERATOR(==),   &
    &                                    timedeltatostring, datetimetostring,              &
    &                                    OPERATOR(*), OPERATOR(<), OPERATOR(<=),           &
    &                                    calendarToString, divideTimeDeltaInSeconds,       &
    &                                    divisionquotienttimespan,                         &
    &                                    mtime_proleptic_gregorian => proleptic_gregorian, &
    &                                    mtime_year_of_365_days => year_of_365_days,       &
    &                                    mtime_year_of_360_days => year_of_360_days,       &
    &                                    getTotalMilliSecondsTimeDelta
  USE mo_time_config,              ONLY: dt_restart,                                       &
    &                                    ini_datetime_string, is_relative_time,            &
    &                                    time_nml_icalendar => icalendar,                  &
    &                                    restart_calendar, restart_ini_datetime_string,    &
    &                                    set_calendar, set_is_relative_time,               &
    &                                    set_tc_dt_model, calendar_index2string
  USE mo_run_config,               ONLY: dtime, mtime_modelTimeStep => modelTimeStep
  USE mo_master_control,           ONLY: atmo_process, get_my_process_type
  USE mo_impl_constants,           ONLY: max_dom, IHS_ATM_TEMP, IHS_ATM_THETA,             &
    &                                    DEFAULT_RESTART_INTVL, DEFAULT_CHECKPT_INTVL,     &
    &                                    inh_atmosphere,                                   &
    &                                    dtime_proleptic_gregorian => proleptic_gregorian, &
    &                                    dtime_cly360              => cly360,              &
    &                                    dtime_julian_gregorian    => julian_gregorian
  USE mo_dynamics_config,          ONLY: iequations
  USE mo_exception,                ONLY: message, message_text, finish
  USE mo_grid_config,              ONLY: patch_weight, grid_rescale_factor, n_dom,         &
    &                                    start_time
  USE mo_io_config,                ONLY: dt_checkpoint
  USE mo_echam_phy_config,         ONLY: echam_phy_config
  USE mo_atm_phy_nwp_config,       ONLY: atm_phy_nwp_config
  USE mo_master_config,            ONLY: experimentReferenceDate,                          &
    &                                    experimentStartDate,                              &
    &                                    checkpointTimeIntval, restartTimeIntval,          &
    &                                    experimentStopDate, isRestart,                    &
    &                                    master_nml_calendar => calendar_str
  USE mo_time_config,              ONLY: set_tc_exp_refdate, set_tc_exp_startdate,         &
    &                                    set_tc_exp_stopdate, set_tc_startdate,            &
    &                                    set_tc_stopdate,                                  &
    &                                    set_tc_dt_restart, set_tc_dt_checkpoint,          &
    &                                    set_tc_current_date, set_tc_write_restart,        &
    &                                    end_datetime_string
  USE mo_io_restart_attributes,    ONLY: get_restart_attribute
  USE mo_io_restart,               ONLY: read_restart_header


  IMPLICIT NONE

  PUBLIC :: compute_timestep_settings
  PUBLIC :: compute_restart_settings
  PUBLIC :: compute_date_settings  

  PRIVATE

  !> module name
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_time_management'

CONTAINS 

  !---------------------------------------------------------------------------------------
  !> Set time step for this run.
  !
  !  This routine is necessary, since (for reasons of backward
  !  compatibility) different naemlist settings may be used to specify
  !  the date information.
  !
  !  Initial revision:  10/2015 : F. Prill, DWD
  !
  SUBROUTINE compute_timestep_settings()
    ! local variables
    CHARACTER(len=*), PARAMETER ::  routine = modname//'::compute_timestep_settings'

    TYPE(timedelta), POINTER             ::  dtime1, dtime2, zero_dt
    CHARACTER(LEN=MAX_TIMEDELTA_STR_LEN) ::  dtime_str2, dtime_string
    INTEGER(i8)                          ::  dtime_ms
    REAL(wp)                             ::  dtime_real
    INTEGER                              ::  jg
    TYPE(datetime), POINTER              ::  reference_dt

    ! --------------------------------------------------------------
    ! PART I: Collect time step as ISO8601 string
    ! --------------------------------------------------------------

    ! The time step is determined either by the "dtime" parameter or
    ! by the "modelTimeStep" parameter from the namelist "run_nml".
    !
    dtime_string = ""
    IF (TRIM(mtime_modelTimeStep) /= "")   dtime_string = TRIM(mtime_modelTimeStep)
    
    dtime_real = dtime
    IF (dtime_real > 0._wp) THEN
      dtime_ms   = NINT(dtime_real*1000, i8)
      CALL getPTStringFromMS(dtime_ms, dtime_str2)
      
      IF (dtime_string == "") THEN
        dtime_string = TRIM(dtime_str2)
      ELSE
        ! Obviously, both namelist parameters for the time step have
        ! been used. We need to test for equality.
        dtime1 => newTimedelta(dtime_string)
        dtime2 => newTimedelta(dtime_str2)
        IF (dtime1 /= dtime2) THEN
          CALL finish(routine, "Two inconsistent time steps have been specified!")
        END IF
        CALL deallocateTimedelta(dtime1)
        CALL deallocateTimedelta(dtime2)
      END IF
    END IF
    dtime1 => newTimedelta(dtime_string)

    ! Furthermore, the time step may be rescaled by the
    ! "grid_rescale_factor".
    !
    IF (grid_rescale_factor /= 1.0_wp) THEN
      dtime1 = dtime1 * grid_rescale_factor
      CALL timedeltaToString(dtime1, dtime_string)
      IF (dtime_real > 0._wp)  dtime_real = dtime_real * grid_rescale_factor

      IF (get_my_process_type() == atmo_process) THEN
        echam_phy_config%dt_rad = &
          & echam_phy_config%dt_rad * grid_rescale_factor
        
        DO jg=1,max_dom
          atm_phy_nwp_config(jg)%dt_conv = &
            atm_phy_nwp_config(jg)%dt_conv * grid_rescale_factor
          atm_phy_nwp_config(jg)%dt_rad  = &
            atm_phy_nwp_config(jg)%dt_rad  * grid_rescale_factor
          atm_phy_nwp_config(jg)%dt_sso  = &
            atm_phy_nwp_config(jg)%dt_sso  * grid_rescale_factor
          atm_phy_nwp_config(jg)%dt_gwd  = &
            atm_phy_nwp_config(jg)%dt_gwd  * grid_rescale_factor
        END DO
      END IF
    ELSE
      CALL timedeltaToString(dtime1, dtime_string)
    END IF

    ! consistency check
    zero_dt => newTimedelta("PT0S") ! mtime object for zero
    IF (dtime1 <= zero_dt) CALL finish(TRIM(routine),'"dtime" must be positive')

    CALL deallocateTimedelta( zero_dt )

    ! --------------------------------------------------------------
    ! PART II: Convert ISO8601 string into "mtime" and old REAL
    ! --------------------------------------------------------------

    CALL set_tc_dt_model(dtime_string)
    IF (dtime_real > 0._wp) THEN
      ! In case that we came from the REAL-valued namelist setting of
      ! the time step we try to avoid rounding errors in floating
      ! point arithmetic
      dtime = dtime_real
    ELSE
      ! For conversion to milliseconds, we need an anchor date,
      ! although we know that "dtime" is too small for this to be
      ! relevant.
      reference_dt => newDatetime("1980-06-01T00:00:00.000")
      dtime = REAL(getTotalMilliSecondsTimeDelta(dtime1, reference_dt),wp)/1000._wp
      CALL deallocateDatetime(reference_dt)
    END IF

    CALL deallocateTimedelta( dtime1  )
    
    ! --------------------------------------------------------------
    ! PART III: Print time step
    ! --------------------------------------------------------------

    CALL message('','')
    WRITE(message_text,'(a,a)') 'Model time step          : ', TRIM(dtime_string)
    CALL message('',message_text)

  END SUBROUTINE compute_timestep_settings


  !---------------------------------------------------------------------------------------
  !> Set restart and checkpoint interval for this run.
  !
  !  This routine is necessary, since (for reasons of backward
  !  compatibility) different naemlist settings may be used to specify
  !  the date information.
  !
  !  Initial revision:  10/2015 : F. Prill, DWD
  !
  SUBROUTINE compute_restart_settings()
    ! local variables:
    CHARACTER(len=*), PARAMETER ::  routine = modname//'::compute_restart_settings'
    CHARACTER(LEN=MAX_TIMEDELTA_STR_LEN) :: restart_intvl_string, checkpt_intvl_string, &
      &                                     checkpt_intvl2, restart_intvl2, dtime_string
    TYPE(timedelta), POINTER             :: mtime_2_5h, mtime_dt_checkpoint,            &
      &                                     mtime_dt_restart, mtime_dom_start,          &
      &                                     tmp_td1, tmp_td2
    TYPE(datetime), POINTER              :: reference_dt
    INTEGER                              :: jg


    ! --------------------------------------------------------------
    ! PART I: Collect the restart and checkpoint intervals 
    !         as ISO8601 strings
    ! --------------------------------------------------------------

    ! --- --- RESTART INTERVAL:
    !
    !         This time interval specifies when the run is supposed to
    !         save its state to a file and stop, later to be resumed.
    !
    !         The restart interval is set by the namelist parameter
    !         "restartTimeIntval" (in the namelist
    !         "master_time_control_nml") - the new way -, and
    !         "dt_restart" (in "time_nml") - the old way.
    !
    !         TODO: The restart interval needs to be multiple of the model
    !               time steps.
    !
    restart_intvl_string = ""
    IF (TRIM(restartTimeIntval) /= "")  restart_intvl_string = TRIM(restartTimeIntval)
    IF (dt_restart > 0._wp) THEN
      restart_intvl2 = "PT"//TRIM(int2string(INT(dt_restart), '(i0)'))//"S"
      IF (TRIM(restart_intvl_string) == "") THEN
        restart_intvl_string = TRIM(restart_intvl2)
      ELSE
        tmp_td1 => newTimedelta(restart_intvl_string)
        tmp_td2 => newTimedelta(restart_intvl2)        
        IF (.NOT. (tmp_td1 < tmp_td2) .AND. .NOT. (tmp_td2 < tmp_td1)) THEN
          restart_intvl_string = TRIM(restart_intvl2)
        ELSE
          CALL finish(routine, "Inconsistent setting of restart interval: " &
               &               //TRIM(restart_intvl_string)//"/"//TRIM(restart_intvl2))
        END IF
        CALL deallocateTimedelta(tmp_td1)
        CALL deallocateTimedelta(tmp_td2)
      END IF
    ELSE IF (dt_restart == 0._wp) THEN
      restart_intvl_string = 'PT0.000S'
      CALL set_tc_write_restart(.FALSE.)
    END IF
    ! if "restart_intvl_string" still unspecified: set default
    IF (TRIM(restart_intvl_string) == "") THEN
      restart_intvl_string = DEFAULT_RESTART_INTVL
    END IF



    ! --- --- CHECKPOINT INTERVAL:
    !
    !         This time interval specifies when the run is supposed to
    !         save its state to a file (but not to stop afterwards).
    !
    !         The checkpoint interval is set by the namelist parameter
    !         "checkpointTimeIntval" (in the namelist
    !         "master_time_control_nml") - the new way -, and
    !         "dt_checkpoint" (in "io_nml") - the old way.
    !
    !         TODO: The checkpoint interval needs to be multiple of the
    !               model time steps.
    !
    checkpt_intvl_string = ""
    IF (TRIM(checkpointTimeIntval) /= "")  checkpt_intvl_string = TRIM(checkpointTimeIntval)
    IF (dt_checkpoint > 0._wp) THEN
      checkpt_intvl2 = "PT"//TRIM(int2string(INT(dt_checkpoint), '(i0)'))//"S"
      IF (TRIM(checkpt_intvl_string) == "") THEN
        checkpt_intvl_string = TRIM(checkpt_intvl2)
      ELSE
        tmp_td1 => newTimedelta(checkpt_intvl_string)
        tmp_td2 => newTimedelta(checkpt_intvl2)        
        IF (.NOT. (tmp_td1 < tmp_td2) .AND. .NOT. (tmp_td2 < tmp_td1)) THEN
          checkpt_intvl_string = TRIM(checkpt_intvl2)
        ELSE
          CALL finish(routine, "Inconsistent setting of checkpoint interval: " &
               &               //TRIM(checkpt_intvl_string)//"/"//TRIM(checkpt_intvl2))
        END IF
        CALL deallocateTimedelta(tmp_td1)
        CALL deallocateTimedelta(tmp_td2)
      END IF
    END IF
    ! if "checkpt_intvl_string" still unspecified: set default
    IF (TRIM(checkpt_intvl_string) == "") THEN
      checkpt_intvl_string = DEFAULT_CHECKPT_INTVL
    END IF

    mtime_dt_restart    => newTimedelta(restart_intvl_string)
    mtime_dt_checkpoint => newTimedelta(checkpt_intvl_string)

    ! consistency checks:
    !
    ! When increased sound-wave and gravity-wave damping is chosen
    ! during the spinup phase (i.e. divdamp_order = 24),
    ! checkpointing/restarting is not allowed earlier than three hours
    ! into the integration because the results would not be
    ! bit-identical in this case
    !
    mtime_2_5h          => newTimedelta("PT02H30M")
    IF ((iequations    == inh_atmosphere) .AND. &
      & (divdamp_order == 24)             .AND. &
      & (mtime_dt_checkpoint < mtime_2_5h)) THEN
        WRITE(message_text,'(a)') &
          &  'dt_checkpoint < 2.5 hours not allowed in combination with divdamp_order = 24'
        CALL finish(routine, message_text)
    ENDIF
    CALL deallocateTimedelta(mtime_2_5h)

    ! Writing a checkpoint file exactly at the start time of a nest is
    ! not allowed:
    !
    DO jg =1,n_dom
      dtime_string = "PT"//TRIM(int2string(INT(start_time(jg)), '(i0)'))//"S"
      mtime_dom_start => newTimedelta(dtime_string)
      IF (mtime_dom_start == mtime_dt_checkpoint) THEN
        WRITE(message_text,'(a)') &
          &  'writing a checkpoint file exactly at the start time of a nest is not allowed'
        CALL finish(routine, message_text)
      END IF
      CALL deallocateTimedelta(mtime_dom_start)
    END DO


    ! If undefined, the restart interval is set to the checkpoint
    ! interval. On the other hand, the checkpoint interval is set to
    ! the restart interval, if the first is not specified:
    !
    IF (TRIM(restart_intvl_string) == "")  restart_intvl_string = checkpt_intvl_string
    IF (TRIM(checkpt_intvl_string) == "")  checkpt_intvl_string = restart_intvl_string


    ! --------------------------------------------------------------
    ! PART II: Convert ISO8601 string into "mtime" and old REAL
    ! --------------------------------------------------------------

    CALL set_tc_dt_restart     ( restart_intvl_string )
    CALL set_tc_dt_checkpoint  ( checkpt_intvl_string )

    ! For conversion to (milli-)seconds, we need an anchor date:
    !
    reference_dt => newDatetime("1980-06-01T00:00:00.000")
    dt_restart    = REAL(getTotalMilliSecondsTimeDelta(mtime_dt_restart,    reference_dt),wp)/1000._wp
    dt_checkpoint = REAL(getTotalMilliSecondsTimeDelta(mtime_dt_checkpoint, reference_dt),wp)/1000._wp
    CALL deallocateTimedelta(mtime_dt_checkpoint)
    CALL deallocateTimedelta(mtime_dt_restart)
    CALL deallocateDatetime(reference_dt)

    ! --------------------------------------------------------------
    ! PART III: Print restart and checkpoint intervals
    ! --------------------------------------------------------------

    WRITE(message_text,'(a,a)') 'Checkpoint interval      : ', TRIM(checkpt_intvl_string)
    CALL message('',message_text)
    WRITE(message_text,'(a,a)') 'Restart interval         : ', TRIM(restart_intvl_string)
    CALL message('',message_text)
    CALL message('','')

  END SUBROUTINE compute_restart_settings


  !---------------------------------------------------------------------------------------
  !> Set start date, end date etc. for this run.
  !
  !  This routine is necessary, since (for reasons of backward
  !  compatibility) different naemlist settings may be used to specify
  !  this date information.
  !
  !  This subroutine checks for contradictory namelist settings and
  !  finally creates 
  !   - data objects of type "t_datetime" 
  !   - data objects from  the mtime library.
  !   - the "nsteps" time loop count
  !
  !  As input this subroutine requires final settings dtime and
  !  dt_restart (these may have been modified in the crosscheck
  !  routine in "mo_nml_crosscheck").
  !
  !  Initial revision:  10/2015 : F. Prill, DWD
  !
  SUBROUTINE compute_date_settings(model_string, dt_restart, nsteps)
    CHARACTER(LEN=*),    INTENT(IN)    :: model_string       !< e.g.  "atm"
    REAL(wp),            INTENT(IN)    :: dt_restart         !< Length of restart cycle in seconds
    INTEGER,             INTENT(INOUT) :: nsteps
    ! local variables
    CHARACTER(len=*), PARAMETER ::  routine = modname//'::compute_date_settings'

    CHARACTER(LEN=MAX_DATETIME_STR_LEN)   ::  ini_datetime1, end_datetime1, dstring
    CHARACTER(LEN=32)                     ::  ini_datetime2, end_datetime2
    CHARACTER(LEN=MAX_DATETIME_STR_LEN)   ::  start_datetime_string,         & !< run start date
      &                                       stop_datetime_string,          & !< run stop date
      &                                       exp_start_datetime_string,     & !< experiment start date
      &                                       exp_stop_datetime_string,      & !< experiment stop date
      &                                       exp_ref_datetime_string,       & !< experiment reference date
      &                                       cur_datetime_string
    CHARACTER(LEN=MAX_TIMEDELTA_STR_LEN)  ::  td_string
    TYPE(datetime),  POINTER              ::  mtime_start, mtime_stop,       &
      &                                       mtime_exp_stop,                &
      &                                       mtime_restart_stop,            &
      &                                       mtime_nsteps_stop,             &
      &                                       tmp_dt1, tmp_dt2 
    TYPE(timedelta), POINTER              ::  mtime_dt_restart, mtime_dtime, &
      &                                       mtime_td
    TYPE(divisionquotienttimespan)        ::  mtime_quotient
    INTEGER                               ::  mtime_calendar, dtime_calendar,&
      &                                       errno
    CHARACTER(len=MAX_CALENDAR_STR_LEN)   ::  calendar1, calendar2, calendar


    ! --------------------------------------------------------------
    ! PART I: Collect all the dates as ISO8601 strings
    ! --------------------------------------------------------------

    ! --- Prelimary: set calendars in old and new date module
    !
    !     The calendar may be set by the namelist parameter "calendar"
    !     in "time_nml" or the namelist parameter "calendar" in
    !     "master_time_control_nml".

    ! Convert the calendar setting (which is an integer value for this
    ! namelist) into a string. The naming scheme is then compatible
    ! with concurrent namelist settings of the calendar (mtime):
    calendar1 = calendar_index2string(time_nml_icalendar)
    calendar2 = TRIM(master_nml_calendar)
    IF (TRIM(calendar1) /= "")  calendar = calendar1
    IF (TRIM(calendar2) /= "")  calendar = calendar2
    IF ((TRIM(calendar1) /= "") .AND. (TRIM(calendar2) /= "")) THEN
      ! both settings were used; we need to test for equality
      IF (TRIM(tolower(calendar1)) /= TRIM(tolower(calendar2)))  &
        &  CALL finish(routine, "Inconsistent setting of calendar")
    END IF
    SELECT CASE (toLower(calendar))
    CASE ('julian gregorian')
      dtime_calendar  = dtime_julian_gregorian
      mtime_calendar  = -1
      CALL finish(routine, "Julian-Gregorian calendar unsupported by mtime library!")
    CASE ('proleptic gregorian')
      dtime_calendar  = dtime_proleptic_gregorian
      mtime_calendar  = mtime_proleptic_gregorian
    CASE ('365 day year')  
      dtime_calendar  = -1
      mtime_calendar  = mtime_year_of_365_days
      !      CALL finish(routine, "Year-of-365-days calendar unsupported by datetime module!")
    CASE ('360 day year')  
      dtime_calendar  = dtime_cly360
      mtime_calendar  = mtime_year_of_360_days
    CASE default
      dtime_calendar       = dtime_proleptic_gregorian
      mtime_calendar       = mtime_proleptic_gregorian
      time_nml_icalendar   = dtime_proleptic_gregorian
      CALL message('','No calendar selected! Use default proleptic Gregorian.')
    END SELECT

    ! --- Now, we merge ISO time stamp strings from the namelists
    !     "time_nml" and "master_time_control_nml". If both namelist
    !     specify the same setting, then we compare the dates and
    !     throw an error in case of inconsistency.  We take also into
    !     account the namelist settings for "nsteps" and "dt_restart".
    !
    !    ***************************  EXPERIMENT RUN   *************
    !    *                                                         *
    !    *                    ***  JOB RUN   ***                   *
    !    *                    *                *                   *
    !  --|--------------------[----------------]-------------------|---------------> (time axis)
    !   exp_start           start            stop                exp_stop
    !                         |----------------|
    !                          restart interval
    !
    ! --- --- EXPERIMENT START DATE:
    !
    !         This is the start date for the whole experiment, which
    !         may be much earlier than the start of the current run.
    !
    !         The experiment start date is set by the namelist
    !         parameter "experimentStartDate" (in namelist
    !         "master_time_control_nml") - the new way -, and
    !         "ini_datetime_string" (in "time_nml") - the old way.
    !
    !         Note that if the namelist parameter
    !         "experimentStartDate" (in "master_time_control_nml") is
    !         unspecified, but the namelist parameter
    !         "experimentReferenceDate" (see also below) is set, then
    !         set the experiment start date to the latter.
    !
    exp_start_datetime_string = ""
    ini_datetime1 = TRIM(experimentStartDate)
    ini_datetime2 = TRIM(ini_datetime_string)
    IF (TRIM(ini_datetime1) /= "")  exp_start_datetime_string = ini_datetime1
    IF (TRIM(ini_datetime2) /= "")  exp_start_datetime_string = ini_datetime2
    IF ((TRIM(ini_datetime1) /= "") .AND. (LEN_TRIM(ini_datetime2) > 0)) THEN
      ! both settings were used; we need to test for equality
      tmp_dt1 => newDatetime(ini_datetime1)
      tmp_dt2 => newDatetime(ini_datetime2)      
      IF (.NOT. (tmp_dt1 == tmp_dt2)) THEN
        CALL finish(routine, "Inconsistent setting of experiment start date: " &
             &               //TRIM(ini_datetime1)//"/"//TRIM(ini_datetime2))
      END IF
      CALL deallocateDatetime(tmp_dt1)
      CALL deallocateDatetime(tmp_dt2)      
    END IF
    IF ( (TRIM(ini_datetime1) == "")        .AND.  &
      &  (TRIM(ini_datetime2) == "")        .AND.  &
      &  (TRIM(experimentReferenceDate) /= "")) THEN
      ! set the experiment start date to the reference date
      exp_start_datetime_string = experimentReferenceDate
    END IF
    ! throw an error, if no start date has been specified at all
    IF (TRIM(exp_start_datetime_string) == "") THEN
      CALL finish(routine, "No experiment start date has been set!")
    END IF


    ! --- --- EXPERIMENT STOP DATE:
    !
    !         This is when the whole experiment ends; note that the
    !         current run may stop before this end date has been
    !         reached (and be restarted later).
    !
    !         The experiment stop date is set by the namelist
    !         parameter "experimentStopDate" (in namelist
    !         "master_time_control_nml") - the new way -, and
    !         "end_datetime_string" (in "time_nml") - the old way.
    !
    exp_stop_datetime_string = ""
    end_datetime1 = TRIM(experimentStopDate)
    end_datetime2 = TRIM(end_datetime_string)
    IF (TRIM(end_datetime1) /= "")  exp_stop_datetime_string = end_datetime1
    IF (TRIM(end_datetime2) /= "")  exp_stop_datetime_string = end_datetime2
    IF ((TRIM(end_datetime1) /= "") .AND. (TRIM(end_datetime2) /= "")) THEN
      ! both settings were used; we need to test for equality
      tmp_dt1 => newDatetime(end_datetime1)
      tmp_dt2 => newDatetime(end_datetime2)      
      IF (.NOT. (tmp_dt1 == tmp_dt2)) THEN
        CALL finish(routine, "Inconsistent setting of experiment stop date: " &
             &               //TRIM(end_datetime1)//"/"//TRIM(end_datetime2))
      END IF
      CALL deallocateDatetime(tmp_dt1)
      CALL deallocateDatetime(tmp_dt2)      
    END IF
    ! throw an error, if no start date has been specified at all
    !
    ! Note: this error check is disabled since some scripts do not
    !       specify the experiment stop date.
    !
    !  IF (TRIM(exp_stop_datetime_string) == "") THEN
    !    CALL finish(routine, "No experiment stop date has been set!")
    !  END IF

    ! --- --- REFERENCE DATE:
    !
    !         This specifies the reference date for the calendar in
    !         use (some kind of "anchor date" on the time line).
    !
    !         If the namelist parameter "experimentReferenceDate (in
    !         "master_time_control_nml") is unspecified, then the
    !         reference date is set to the experiment start date.
    !
    IF (TRIM(experimentReferenceDate) == '') THEN
      exp_ref_datetime_string = exp_start_datetime_string
    ELSE
      exp_ref_datetime_string = experimentReferenceDate
    END IF

    ! --- --- START DATE:
    !
    !         This is the start date of the current run.
    !
    !         a) in case of restart:  start date := restart attribute "tc_startdate"
    !         b) no-restart case:     start date := experiment start date
    !
    !         Special case: In a restart situation, if the calendar or
    !         initial date/time is different from those in the restart
    !         file, we regard this integration as a new one with its
    !         own calendar.  Model time at which the previous run
    !         stopped is thus not relevant.  Simulation will start
    !         from the user-specified initial date/time, which is also
    !         the current model date/time.
    !
    IF (isRestart()) THEN

      IF ((restart_calendar /= dtime_calendar) .AND. &
        & (restart_calendar /= -1)) THEN
        
        CALL message('','Restart calendar is not matching, fallback to experiment start date: '//&
          &int2string(restart_calendar)//' /= '//int2string(dtime_calendar))
        start_datetime_string = exp_start_datetime_string

      ELSE IF (TRIM(restart_ini_datetime_string) /= TRIM(ini_datetime_string)) THEN

        CALL message('','Restart ini date is not matching, fallback to experiment start date: '//&
        &TRIM(restart_ini_datetime_string)//' /= '//TRIM(ini_datetime_string))

      ELSE
        CALL message('','Read restart file meta data ...')
        CALL read_restart_header(TRIM(model_string))
        CALL get_restart_attribute('tc_startdate', start_datetime_string)
      END IF

    ELSE
      start_datetime_string = exp_start_datetime_string
    ENDIF

    ! --- --- CURRENT DATE:
    !
    !         This is the current model date, that is stepped forward
    !         during the integration loop. We begin by setting it to
    !         the start date computed above.
    !
    cur_datetime_string = start_datetime_string

    ! --- --- STOP DATE:
    !
    !         This is the date when the current run stops time stepping.
    !
    !         The simulation stops at the next (i.e. earliest) of the
    !         following dates:
    !         a) experiment stop date
    !         b) start date + dt_restart
    !         c) start date + nsteps*dtime
    !
    NULLIFY(mtime_exp_stop)
    NULLIFY(mtime_restart_stop)
    NULLIFY(mtime_nsteps_stop)

    CALL getPTStringFromMS(INT(dtime,i8)*1000, td_string)
    mtime_dtime => newTimeDelta(td_string)
    IF (.NOT. ASSOCIATED(mtime_dtime))  CALL finish(routine, "Error in conversion of dtime to mtime!")
    mtime_start        => newDatetime(start_datetime_string, errno)
    IF (errno /= 0)  CALL finish(routine, "Error in conversion of start date: "//start_datetime_string)
    IF (TRIM(exp_stop_datetime_string) /= "") THEN
      mtime_exp_stop     => newDatetime(exp_stop_datetime_string, errno)
      IF (errno /= 0)  CALL finish(routine, "Error in conversion of exp stop date: "//exp_stop_datetime_string)
    END IF

    IF (INT(dt_restart) > 0) THEN
      mtime_dt_restart   => newTimedelta("PT"//TRIM(int2string(INT(dt_restart),'(i0)'))//"S", errno)
      IF (errno /= 0)  CALL finish(routine, "Error in conversion of dt_restart!")
      mtime_restart_stop => newDatetime(mtime_start, errno)
      IF (errno /= 0)  CALL finish(routine, "Error in initialization of restart date")
      mtime_restart_stop =  mtime_restart_stop + mtime_dt_restart
    ELSE
      IF (TRIM(exp_stop_datetime_string) /= "") THEN
        mtime_restart_stop => newDatetime(exp_stop_datetime_string, errno)
        IF (errno /= 0)  CALL finish(routine, "Error in initialization of restart date")
      END IF
    END IF
    IF (nsteps >= 0) THEN   

      ! Special treatment for the hydro atm model
      !
      ! TODO: Is this weird workaround really needed?
      IF ( (iequations == IHS_ATM_TEMP) .OR. &
        &  (iequations == IHS_ATM_THETA)     ) THEN
        
        ! If running the HYDROSTATIC version, let the model integrate
        ! one more step after the desired end of simulation in order
        ! to get the proper output. This additional step is necessary
        ! because the HYDROSTATIC model writes out values of step N
        ! after the integration from N to N+1 is finished. Also note
        ! that this additional step is done only for the regular
        ! output, and is ignored for restart.
        nsteps = nsteps + 1

        ! The additional step is not needed in the NON-hydrostatic
        ! version because in this case the model writes out values of
        ! step N after the integration from N-1 to N is finished.
      END IF
      mtime_nsteps_stop  => newDatetime(mtime_start, errno)
      IF (errno /= 0)  CALL finish(routine, "Error in initialization of nsteps stop date")
      mtime_nsteps_stop = mtime_nsteps_stop + mtime_dtime * INT(nsteps,c_int32_t)
    ELSE
      IF (TRIM(exp_stop_datetime_string) /= "") THEN
        mtime_nsteps_stop  => newDatetime(mtime_exp_stop, errno)
        IF (errno /= 0)  CALL finish(routine, "Error in initialization of nsteps  stop date")
      END IF
    END IF

    ! now, we could simply set
    !
    !   mtime_stop => newDatetime(MIN(MIN(mtime_exp_stop, mtime_restart_stop), mtime_nsteps_stop))
    !
    ! but we need to check cases where one or two of these dates have
    ! not been specified by the user...
    IF (.NOT. ASSOCIATED(mtime_exp_stop)     .AND. &
      & .NOT. ASSOCIATED(mtime_restart_stop) .AND. &
      & .NOT. ASSOCIATED(mtime_nsteps_stop)) THEN
      CALL finish(routine, "Error in initialization of stop date")
    ELSE IF (ASSOCIATED(mtime_exp_stop)      .AND. &
      & ASSOCIATED(mtime_restart_stop)       .AND. &
      & ASSOCIATED(mtime_nsteps_stop)) THEN
      mtime_stop => newDatetime(MIN(MIN(mtime_exp_stop, mtime_restart_stop), mtime_nsteps_stop))
    ELSE IF (ASSOCIATED(mtime_exp_stop)      .AND. &
      & ASSOCIATED(mtime_restart_stop)) THEN
      mtime_stop => newDatetime(MIN(mtime_exp_stop, mtime_restart_stop))
    ELSE IF (ASSOCIATED(mtime_exp_stop)      .AND. &
      & ASSOCIATED(mtime_nsteps_stop)) THEN
      mtime_stop => newDatetime(MIN(mtime_exp_stop, mtime_nsteps_stop))
    ELSE IF (ASSOCIATED(mtime_restart_stop)  .AND. &
      & ASSOCIATED(mtime_nsteps_stop)) THEN
      mtime_stop => newDatetime(MIN(mtime_restart_stop, mtime_nsteps_stop))
    END IF

    CALL datetimeToString(mtime_stop, stop_datetime_string)

    ! if it has not been specified by the user, set experiment stop
    ! date to stop date:
    IF (.NOT. ASSOCIATED(mtime_exp_stop)) THEN
      mtime_exp_stop => newDatetime(stop_datetime_string)
    END IF


    ! consistency checks:
    !
    IF (mtime_stop < mtime_start) THEN
      CALL finish(routine, 'The end date and time must not be '// &
        &                  'before the current date and time')
    END IF
    ! If a restart event occurs, check for unsupported combinations of
    ! namelist settings:
    IF (.NOT. (mtime_nsteps_stop < mtime_restart_stop)) THEN
      ! processor splitting cannot be combined with synchronous restart:
      IF ((num_restart_procs == 0) .AND. ANY(patch_weight(1:) > 0._wp)) THEN
        CALL finish(routine, "Processor splitting cannot be combined with synchronous restart!")
      END IF
    END IF
    CALL deallocateDatetime(mtime_start)
    IF (ASSOCIATED(mtime_restart_stop))  CALL deallocateDatetime(mtime_exp_stop)
    IF (ASSOCIATED(mtime_restart_stop))  CALL deallocateDatetime(mtime_restart_stop)
    IF (ASSOCIATED(mtime_restart_stop))  CALL deallocateDatetime(mtime_nsteps_stop)
    CALL deallocateDatetime(mtime_stop)
    IF (INT(dt_restart) > 0)  CALL deallocateTimedelta(mtime_dt_restart)
    CALL deallocateTimedelta(mtime_dtime)

    ! --- --- NSTEPS
    !
    !         This specifies the number of steps in the time loop
    !         (INTEGER).
    !
    !         If the namelist parameter "nsteps" has not been
    !         specified in the namelist "run_nml", then this value is
    !         computed from the stop date and the start date. On the
    !         other hand, if the namelist parameter "nsteps" is
    !         specified in the namelist "run_nml", then this defines
    !         the end date (see above).
    !
    IF (nsteps < 0) THEN   
      ! User did not specified a value, we need to compute "nsteps" as
      ! (stop date - start date)/dtime:
      mtime_start => newDatetime(start_datetime_string)
      mtime_stop  => newDatetime(stop_datetime_string)
      mtime_td    => newTimedelta("PT0S")
      mtime_td    =  mtime_stop - mtime_start
      CALL getPTStringFromMS(INT(dtime*1000._wp,i8), td_string)
      mtime_dtime => newTimedelta(td_string)
      CALL divideTimeDeltaInSeconds(mtime_td, mtime_dtime, mtime_quotient)
      nsteps = INT(mtime_quotient%quotient)
      CALL deallocateTimedelta(mtime_td)
      CALL deallocateTimedelta(mtime_dtime)
      CALL deallocateDatetime(mtime_start)
      CALL deallocateDatetime(mtime_stop)
    END IF

    ! --------------------------------------------------------------
    ! PART II: Convert ISO8601 strings into "mtime" and "t_datetime" 
    ! --------------------------------------------------------------

    ! --- Second, create date-time objects from the mtime library

    CALL set_tc_startdate    ( start_datetime_string     )
    CALL set_tc_stopdate     ( stop_datetime_string      )
    CALL set_tc_exp_startdate( exp_start_datetime_string )
    CALL set_tc_exp_stopdate ( exp_stop_datetime_string  )
    CALL set_tc_exp_refdate  ( exp_ref_datetime_string   )
    CALL set_tc_current_date ( cur_datetime_string       )

    ! --- Finally, store the same information in a "t_datetime" data
    !     structure
    CALL set_calendar        (dtime_calendar)
    CALL set_is_relative_time(is_relative_time)

    ! --------------------------------------------------------------
    ! PART III: Print all date and time components
    ! --------------------------------------------------------------

    CALL calendarToString(dstring)
    CALL message('','Calendar: '//TRIM(dstring))
    call message('','')

    WRITE(message_text,'(a,a)') 'Experiment reference date: ', TRIM(exp_ref_datetime_string)
    CALL message('',message_text)
    WRITE(message_text,'(a,a)') 'Experiment start date    : ', TRIM(exp_start_datetime_string)
    CALL message('',message_text)
    WRITE(message_text,'(a,a)') 'Experiment stop date     : ', TRIM(exp_stop_datetime_string)
    CALL message('',message_text)
    CALL message(' ',' ')

    CALL message(routine,'Initial date and time')
    CALL message('','---------------------')

    CALL message('',message_text)
    IF (isRestart()) THEN
      WRITE(message_text,'(a,a,a)') 'Start date    : ', TRIM(start_datetime_string), ' (restart run)'
    ELSE
      WRITE(message_text,'(a,a)') 'Start date    : ', TRIM(start_datetime_string)
    END IF
    CALL message('',message_text)

    CALL message(' ',' ')
    CALL message(routine,'End date and time')
    CALL message('','-----------------')
    WRITE(message_text,'(a,a)') 'Stop date     : ', TRIM(stop_datetime_string)
    CALL message('',message_text)
    CALL message(' ',' ')
    
  END SUBROUTINE compute_date_settings
 
END MODULE mo_time_management
