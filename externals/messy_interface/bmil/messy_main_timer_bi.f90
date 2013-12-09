MODULE messy_main_timer_bi

  !  Author: Astrid Kerkweg,  IPA, Uni-Mainz, 2008
  !  
  ! MAIN-TIMER is based on the date and time management of ECHAM5
  ! For the use within MESSy the basic subroutines have been heavily
  ! simplified.
  ! MAIN_TIMER comprises 3 core modules:
  ! - messy_main_timer: 
  !                       This subroutine defines the basic date and
  !                       variables needed within a model run to coordinate
  !                       the start, stop, restarts and defined EVENTS
  ! - messy_main_timer_manager: 
  !                       It determines the "heart beat" of the  model.
  !                       Based on a given start date and a time step
  !                       messy_main_timer_manager will always provide
  !                       the current date and time step number
  ! - messy_main_timer_event:
  !                       This module provides the interface routines to
  !                       determine and manage certain events (e.g. restarts,
  !                       data output or emissions)


  ! BM/MESSy
  USE messy_main_mpi_bi,       ONLY: p_parallel_io
  USE messy_main_blather_bi,   ONLY: error_bi, info_bi, warning_bi &
       , start_message_bi, end_message_bi

  ! MESSy
  USE messy_main_constants_mem, ONLY: dp, i4, OneDay              &
       , STRLEN_ULONG, STRLEN_SHORT

  ! CORE
  USE messy_main_timer_manager
  USE messy_main_timer_event
  USE messy_main_timer

  IMPLICIT NONE
  PRIVATE
  SAVE

  TYPE(time_manager), PUBLIC :: BM_time  !  base model time axis

  ! CALLED FROM CONTROL
  PUBLIC :: main_timer_setup
  PUBLIC :: main_timer_initialize
  PUBLIC :: main_timer_global_start
  !
  PUBLIC :: messy_timer_init_manager
  PUBLIC :: messy_timer_reset_time
#ifdef COSMO
  PUBLIC :: messy_timer_COSMO_reinit_time
#endif
  ! PRIVATE :: timer_read_namelist
  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------
  ! event_state
  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------
  PUBLIC :: p_bcast_event
  PUBLIC :: event_state                     ! get the state of event
  INTERFACE event_state
     !  get interval
     MODULE PROCEDURE I_event_steps          ! int<-(I:event,I:delta) in steps
     MODULE PROCEDURE I_event_seconds        ! int<-(I:event)  in seconds
     MODULE PROCEDURE I_event_seconds_next   ! int<-(I:event,I:true) future (s)
     !  evaluate next trigger, check date
     MODULE PROCEDURE L_event_trigger        ! log<-(I:event,I:time_days)
  END INTERFACE

  PUBLIC  :: event_eval        ! int<-(I:event,I:delta) fit delta into interval
  PUBLIC  :: event_print       ! (I:event[,I:format]) print event
  PUBLIC  :: timer_event_init

  ! PUBLIC HELPER ROUTINE
  PUBLIC :: timer_set_rerun_event
  PUBLIC :: timer_get_rerun_event
  PUBLIC :: timer_message
  PUBLIC :: timer_get_time_step

  ! PRIVATE :: write_date
  ! PRIVATE :: init_BM_manager
  ! PRIVATE :: reset_BM_manager


  !TYPE(time_event)     :: TEST_EVENT 
  !TYPE(io_time_event)  :: IO_TEST_EV = &
  !     io_time_event(1, TIME_INC_MONTHS,TRIG_EXACT,0) ! DEFAULT

  TYPE(time_event)     :: RERUN_EVENT 
  TYPE(io_time_event)  :: IO_RERUN_EV = &
       io_time_event(1, TIME_INC_MONTHS,TRIG_EXACT,0) ! DEFAULT
  INTEGER :: MODEL_START(6) = 0 
  INTEGER :: MODEL_STOP(6)  = 0
  TYPE(time_days) :: dummy_date    ! um_ak_20120307

  ! DUMMY INTEGER (see no_cycles)
  INTEGER         :: ccount = 1

CONTAINS

  ! -------------------------------------------------------------------------
  SUBROUTINE main_timer_setup(flag)
    
    ! MESSy/BMIL
    USE messy_main_mpi_bi,    ONLY: p_parallel_io, p_bcast, p_io
    USE messy_main_blather_bi, ONLY: error_bi
    ! MESSy/SMCL
    USE messy_main_tools,      ONLY: find_next_free_unit
#ifdef COSMO
    USE data_modelconfig,      ONLY: dt
#endif


    IMPLICIT NONE

    INTEGER, INTENT(IN) :: flag

    INTEGER :: iou
    INTEGER :: status
    LOGICAL :: lsmaller
    LOGICAL :: linit

    CHARACTER(LEN=*), PARAMETER :: substr='messy_timer_setup'

    SELECT CASE(flag)
    CASE(1)
       ! INITIALIZE CPL
       IF (p_parallel_io) THEN
          iou = find_next_free_unit(100,200)
          CALL main_timer_read_nml_cpl(status, iou)
          IF (status /= 0) CALL error_bi(' ',substr)
       ! BROADCAST RESULTS
#ifndef ECHAM5
          IF (delta_time < 0) THEN
             CALL info_bi('delta_time < 0 ',substr)
             CALL info_bi('delta_time must be given in namelist',substr)
             CALL error_bi('delta_time < 0 ',substr)
          ENDIF
#endif
       END IF

       CALL p_bcast(CAL_TYPE   , p_io)
       CALL p_bcast(delta_time , p_io)
       CALL p_bcast(lresume    , p_io)
       CALL p_bcast(NO_CYCLES  , p_io)
       CALL p_bcast(NO_DAYS    , p_io) ! um_ak_20120301
       CALL p_bcast(NO_STEPS   , p_io) ! um_ak_20120301
       CALL p_bcast(LABORT     , p_io)
#ifdef COSMO
       dt = delta_time
#endif

       ! set logicals
       lfirst_cycle = .TRUE.
       IF (.NOT. lresume) THEN
          lstart = .TRUE.
       ELSE
          lstart = .FALSE.
       ENDIF

    CASE(2)

       CALL is_init(start_date, linit)
       
       IF (.NOT. linit) THEN

          CALL p_bcast(MODEL_START, p_io)
          
          IF (SUM(MODEL_START) == 0 ) THEN
             CALL error_bi(substr, 'MODEL START DATE NOT GIVEN IN NAMELIST')
          ELSE
             YEAR_START   = MODEL_START(1)
             MONTH_START  = MODEL_START(2)
             DAY_START    = MODEL_START(3)
             HOUR_START   = MODEL_START(4)
             MINUTE_START = MODEL_START(5)
             SECOND_START = MODEL_START(6)
          
             CALL timer_set_date(status,'start',MODEL_START(1), MODEL_START(2)&
              , MODEL_START(3), MODEL_START(4), MODEL_START(5), MODEL_START(6))
             ! set current time first to start_date, will be overwritten later
             CALL timer_get_date(status,'start' &
                  ,YEAR,MONTH,DAY,HOUR,MINUTE,SECOND)
             CALL timer_set_date(status,'resume' &
                  ,YEAR,MONTH,DAY,HOUR,MINUTE,SECOND)
          ENDIF
          
       ENDIF
       
       CALL is_init(stop_date, linit)
       
       IF (.NOT. linit) THEN
          CALL p_bcast(MODEL_STOP , p_io)

          IF (SUM(MODEL_STOP) == 0 ) THEN
! um_ak_20120301+
!!$             CALL error_bi(substr, 'MODEL STOP DATE NOT GIVEN IN NAMELIST')
!!$          ELSE
!!$             CALL timer_set_date(status, 'stop', MODEL_STOP(1), MODEL_STOP(2) &
!!$                  , MODEL_STOP(3), MODEL_STOP(4), MODEL_STOP(5), MODEL_STOP(6))
             CALL warning_bi(substr, 'MODEL STOP DATE NOT GIVEN IN NAMELIST')
             IF (NO_DAYS > 0) THEN
                CALL copy_date (resume_date, dummy_date)
                IF (lstart) THEN
                   CALL add_date  (no_days, 0, dummy_date)
                ELSE
                   CALL add_date  (no_days, INT(delta_time), dummy_date)
                ENDIF
                CALL warning_bi(substr, 'Using NO_DAYS for model stop.')
             ELSE IF (NO_STEPS > 0) THEN
                CALL copy_date (resume_date, dummy_date)
                CALL add_date  ( 0, no_steps*INT(delta_time), dummy_date)
                CALL warning_bi(substr, 'Using NO_STEPS for model stop.')
             ELSE
                CALL error_bi(substr, 'MODEL STOP DATE NOT GIVEN IN NAMELIST')
             ENDIF
             CALL timer_get_date (status, dummy_date&
                  , MODEL_STOP(1), MODEL_STOP(2) &
                  , MODEL_STOP(3), MODEL_STOP(4) &
                  , MODEL_STOP(5), MODEL_STOP(6))
          ENDIF
       ENDIF
! um_ak_20120307+
       CALL timer_set_date(status, 'stop', MODEL_STOP(1), MODEL_STOP(2) &
            , MODEL_STOP(3), MODEL_STOP(4), MODEL_STOP(5), MODEL_STOP(6))
! um_ak_20120307-

       CALL if_less(stop_date,start_date,lsmaller,status)
       CALL timer_message(status,substr)
       IF (lsmaller) CALL error_bi( &
            'stop_date smaller than start_date', substr)
          
       IF (lstart) &
            CALL messy_timer_init_manager(INT(delta_time), INIT_STEP)

       ! BROADCAST RERUN EVENT
       CALL p_bcast_event(IO_RERUN_EV, p_io)

       ! BROADCAST TEST EVENT
       !CALL p_bcast_event(IO_TEST_EV, p_io)

    END SELECT

  END SUBROUTINE main_timer_setup
  ! -------------------------------------------------------------------------

  ! -------------------------------------------------------------------------
  SUBROUTINE messy_timer_init_manager (timestep, step)      

    ! time manager initialization
    INTEGER, INTENT(IN) :: timestep   ! delta time in seconds
    INTEGER, INTENT(IN) :: step       ! model time step

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'messy_timer_init_manager'
    INTEGER :: status
    CHARACTER(LEN=STRLEN_ULONG) :: m_text

    IF (p_parallel_io) &
         CALL write_date(start_date, 'start_date'//substr)

    CALL manager_init(BM_time,'base model time', start_date, timestep, m_text)
    CALL timer_message (0, substr, m_text)

    IF (lresume) THEN
       !
       CALL write_date(resume_date,'Experiment resumed at: ')
       !
       IF (lfirst_cycle) THEN
          CALL manager_init(BM_time, step+1,status)
          CALL timer_message (status, substr)
          CALL info_bi(' Set manager to resumed time step + 1.', 'TIMER')
       END IF
    ELSE
       CALL manager_state(BM_time,resume_date,step,status)
       CALL timer_message (status, substr)
    END IF

  END SUBROUTINE messy_timer_init_manager
  ! -------------------------------------------------------------------------

  ! -------------------------------------------------------------------------
  SUBROUTINE main_timer_initialize

    ! BM/MESSy
    ! op_bk_20130820+
    USE messy_main_mpi_bi,    ONLY: p_parallel_io, p_bcast, p_io

    IMPLICIT NONE
    INTRINSIC :: TRIM

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='main_timer_initialize'
    INTEGER                     :: status
    INTEGER                     :: iou    ! I/O unit
    CHARACTER(len=STRLEN_SHORT) :: ev_adjust

    CALL start_message_bi(modstr, 'INITIALIZE TIMER',' ') 

    CALL init_BM_time

    ! SET START TIME
    CALL timer_get_date(status, 'start', YEAR_START,MONTH_START,DAY_START, &
         HOUR_START,MINUTE_START,SECOND_START)

    ! um_ak_20081104+
    ! set current time first to start_date, will be overwritten later
    ! qqq needed here for NEW TIMER??
    IF (lstart) &
    CALL timer_get_date(status, 'start', YEAR,MONTH,DAY,HOUR,MINUTE,SECOND)
    ! SET CURRENT TIME STEP
    current_time_step = timer_get_time_step()
    DAYOFYEAR=NINT(YearDay(current_date))
    ! um_ak_20081104-
    CALL timer_message(status, substr)

    CALL end_message_bi(modstr, 'INITIALIZE TIMER',' ') 

    CALL start_message_bi(modstr, 'INITIALIZE RERUN',' ') 

    CALL timer_event_init(  RERUN_EVENT,  IO_RERUN_EV &
         , 'rerun_event', EV_TLEV_NEXT)

    !CALL timer_event_init(  TEST_EVENT,  IO_TEST_EV &
    !     , 'test_event', EV_TLEV_PREV)

    CALL end_message_bi(modstr, 'INITIALIZE RERUN',' ') 

  END SUBROUTINE main_timer_initialize
  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------

  SUBROUTINE main_timer_global_start

    ! MAKE TIME STEPPING FOR CLOCK AND UPDATE EVENTS

    ! evaluate events and date/time control elements at the
    ! beginning of the time loop 

    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER :: substr='main_timer_global_start'
    INTEGER         :: istep
    INTEGER         :: status
    LOGICAL         :: lsmaller, lequal !, l_test

    ! evaluate events
    IF (.NOT. lforcedtime) THEN
       l_rerun  = event_state(RERUN_EVENT,  next_date) .OR. L_TRIGGER_RESTART 
    ELSE
       l_rerun = .FALSE.
    ENDIF
    !l_test   = event_state(TEST_EVENT,  previous_date)
    !IF (p_parallel_io) write (*,*) ' TEST EVENT previous', l_test 

    CALL timer_get_date(status, 'next'        &
         ,YEAR_NEXT,MONTH_NEXT,DAY_NEXT,HOUR_NEXT,MINUTE_NEXT,SECOND_NEXT)
    CALL timer_message(status, substr)
    
    !IF (p_parallel_io) write (*,*) ' RERUN EVENT NEXT ', l_rerun &
    !     , YEAR_NEXT,MONTH_NEXT,DAY_NEXT,HOUR_NEXT,MINUTE_NEXT,SECOND_NEXT
    !CALL event_print(RERUN_EVENT, status, .FALSE., p_parallel_io)


    ! evaluate the stop of the model
    CALL if_less(stop_date,next_date,lsmaller,status)
    CALL timer_message(status,substr)
    CALL if_equal(stop_date,next_date,lequal,status)
    CALL timer_message(status,substr)
    lstop = (lsmaller .OR. lequal)

    !mz_ak20090528+
    IF (lstop) THEN
       l_rerun = .TRUE.
       L_TRIGGER_RESTART  = .TRUE.
    ENDIF
    !mz_ak20090528-

    IF (l_rerun) THEN
       IF (no_cycles > 0) THEN
          lcycbreak  = (ccount >= no_cycles)   ! break rerun cycles
          ccount = ccount + 1                  ! count rerun intervals
       ELSE
          lcycbreak  = .FALSE.                 ! break rerun cycles
       END IF
       lbreak = lcycbreak .OR. L_TRIGGER_RESTART
    ELSE
       lbreak = .FALSE.
    END IF

    ! print settings
    IF (lstop) THEN
       CALL write_date(next_date,'Stop model, last prognostic date/time is: ')
    ELSE IF (lbreak) THEN
       CALL write_date(next_date &
            ,'Interrupt model, last prognostic date/time is: ')
    END IF

    ! calculate integration interval 
    CALL manager_state(BM_time,delta_time,status)
    CALL timer_message(status,substr)

    ! SET CURRENT DATE / TIME
    CALL timer_get_date(status, 'current'     &
         ,YEAR,MONTH,DAY,HOUR,MINUTE,SECOND)
    CALL timer_message(status, substr)
    ! SET NEXT DATE / TIME
    CALL timer_get_date(status, 'next'        &
         ,YEAR_NEXT,MONTH_NEXT,DAY_NEXT,HOUR_NEXT,MINUTE_NEXT,SECOND_NEXT)
    CALL timer_message(status, substr)

    ! op_bk_20131119+
    WRITE(*,*) "MESSY: main_timer_global_start"
    WRITE(*  &
         , '(i4,a1,i2.2,a1,i2.2,a1,i2.2,a1,i2.2,a1,i2.2)')&
         YEAR_NEXT, '-', MONTH_NEXT, '-',  DAY_NEXT, ' ' &
         , HOUR_NEXT, ':', MINUTE_NEXT, ':', SECOND_NEXT
    ! op_bk_20131119-

    ! SET CURRENT TIME STEP
    current_time_step = timer_get_time_step()

    ! GET DAY OF YEAR (e.g. 1 Feb = 32)
    DAYOFYEAR=NINT(YearDay(current_date))
 
  END SUBROUTINE main_timer_global_start

  !----------------------------------------------------------------------------

  SUBROUTINE messy_timer_reset_time

    ! evaluate actions at the end of the time step loop 

    INTEGER :: istep, status
    CHARACTER(LEN=*), PARAMETER :: substr='messy_timer_reset_time'

    ! reset time manager and dates
    CALL manager_state(BM_time,previous_date,status)
    CALL timer_message(status, substr)

    IF (p_parallel_io .AND. LDEBUG) CALL write_date(previous_date &
         ,'reset_time previous date: ')   

    ! increment time step
    CALL manager_init (BM_time,1, status)       
    CALL timer_message(status, substr)

    ! redefine time window
    CALL manager_state(BM_time, istep,status)
    CALL timer_message(status, substr)
    CALL manager_state(BM_time, current_date,status)
    CALL timer_message(status, substr)
    IF (p_parallel_io .AND. LDEBUG) CALL write_date(current_date &
         ,'reset_time current date: ')   
    CALL manager_state(BM_time, next_date,istep+1,status)
    CALL timer_message(status, substr)
    IF (p_parallel_io .AND. LDEBUG) CALL write_date(next_date &
         ,'reset_time next date: ')   

    IF (p_parallel_io .AND. LDEBUG) THEN
       CALL manager_print(BM_time, status)
       CALL timer_message(status,substr)
    ENDIF

    ! reset switches 
    lstart       = .FALSE.
    lresume      = .FALSE.
    lfirst_cycle = .FALSE.

  END SUBROUTINE messy_timer_reset_time
  ! --------------------------------------------------------------------------

  FUNCTION timer_get_time_step() RESULT (istep)  !****************************

    ! provide time step from echam time manager
    !-
    IMPLICIT NONE

    INTEGER :: istep
    INTEGER :: status
    CHARACTER(LEN=*), PARAMETER :: substr='main_timer: get_time_step'

    CALL manager_state(BM_time,istep,status)
    CALL timer_message(status,substr)

  END FUNCTION timer_get_time_step

  ! ---------------------------------------------------------------------------

#ifdef COSMO
  ! --------------------------------------------------------------------------
  SUBROUTINE messy_timer_COSMO_reinit_time

    ! COSMO
    USE data_runcontrol,    ONLY: hstart,   nstart  &
         , hlastmxu,   hnextmxu,   hincmxu   & 
         , hlastmxt,   hnextmxt,   hincmxt   & 
         , nlastmxu,   nnextmxu,   nlastmxt  &
         , nnextmxt,   hnextrad,   nextrad   & 
         , hlastbound, nlastbound, nincbound &
         , hnextbound, nnextbound, hincbound &
         , yakdat1,    yakdat2, itype_calendar

    USE data_modelconfig,   ONLY: dt  
    USE data_io,            ONLY: ydate_ini
    USE utilities,          ONLY: get_utc_date
    ! MESSy
    USE messy_main_blather_bi,   ONLY: info_bi

    IMPLICIT NONE
    INTRINSIC :: NINT, REAL

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'reinit_COSMO_time'
    REAL(dp) :: timediff ! time difference in second (resume-start date) 
    INTEGER  :: ryr, rmo, rdy, rhr, rmn, rse  ! resume date components
    INTEGER  :: dummy1
    REAL(dp) :: dummy2
    INTEGER  :: status
    CHARACTER(LEN=60) :: text

    CALL timer_get_date(status, 'resume', ryr, rmo, rdy, rhr, rmn, rse)
    CALL timer_message(status, substr)
    CALL time_span_d( timediff , YEAR_START, MONTH_START,  DAY_START    &
         , HOUR_START, MINUTE_START, SECOND_START &
         , ryr, rmo,   rdy, rhr,     rmn, rse     )
    hstart = timediff * 24.0_dp ! resume hour
    nstart = NINT(timediff*86400.0_dp /dt)

    ! RE-DO some CHECKS WHICH ARE PERFORMED ORIGIALLY IN ORGANIZE_SETUP
    ! endless loop for finding the last hour (for restart runs)

    hlastmxu = 0.0_dp
    endless_u: DO
!um_ak_20090427 IF ((hlastmxu <= hstart).AND.(hstart < hlastmxu + hincmxu)) THEN
       IF ( (hlastmxu <= hstart) .AND. (hstart <= hlastmxu + hincmxu) ) THEN
          EXIT endless_u
       ENDIF
       hlastmxu = hlastmxu + hincmxu
    ENDDO endless_u
    hnextmxu = hlastmxu + hincmxu
    nlastmxu = NINT (hlastmxu * 3600.0_dp / dt)
    nnextmxu = NINT (hnextmxu * 3600.0_dp / dt)

    hlastmxt = 0.0_dp
    endless_t: DO
!um_ak_20090427 IF ((hlastmxt <= hstart).AND.(hstart < hlastmxt + hincmxt)) THEN
       IF ( (hlastmxt <= hstart) .AND. (hstart <= hlastmxt + hincmxt) ) THEN
          EXIT endless_t
       ENDIF
       hlastmxt = hlastmxt + hincmxt
    ENDDO endless_t
    hnextmxt = hlastmxt + hincmxt
    nlastmxt = NINT (hlastmxt * 3600.0_dp / dt)
    nnextmxt = NINT (hnextmxt * 3600.0_dp / dt)

    ! compute the actual date
    CALL get_utc_date(nstart, ydate_ini, dt, itype_calendar, yakdat1 &
         , yakdat2, dummy1, dummy2)

  END SUBROUTINE messy_timer_COSMO_reinit_time
  ! --------------------------------------------------------------------------
#endif
  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------
  ! EVENT_MANAGER
  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------
  SUBROUTINE timer_event_init (event, io_ev, ev_name, eval_date)

    ! the initial state of an event is calculated
    ! it is dependent on the initial date of a run
    ! at rerun the event triggers will be recalculated from beginning

    IMPLICIT NONE
    INTRINSIC :: NINT, TRIM

    TYPE(time_event),    INTENT(inout) :: event     ! will be initialized
    TYPE(io_time_event), INTENT(inout) :: io_ev     ! pass external settings
    CHARACTER(len=*)                   :: ev_name   ! submit a name
    CHARACTER(len=*)                   :: eval_date ! the evaluation date
    ! (next or present)

    TYPE(time_days)             :: my_date
    REAL(dp)                    :: zdtime, zztime
    CHARACTER(len=STRLEN_SHORT) :: newunit
    INTEGER                     :: newcount, idtime
    CHARACTER(LEN=*),PARAMETER  :: substr='timer_event_init'
    CHARACTER(len=STRLEN_ULONG) :: message
    INTEGER                     :: status
    LOGICAL                     :: lsm

    IF (p_parallel_io) write (*,*) substr, ' EVENT NAME: ', ev_name

    CALL manager_state(BM_time,zztime,status)
    CALL timer_message(status,substr)
    zdtime = zztime - 1.0_dp
    idtime = NINT(zztime)

    ! convert from steps into other units
    IF (TRIM(io_ev%unit) == TIME_INC_STEPS &
         .OR. TRIM(io_ev%unit) == TIME_INC_ALWAYS) THEN

       IF (TRIM(io_ev%unit) == TIME_INC_ALWAYS) THEN 
          ! for debugging and testing only
             io_ev%counter = 1
       END IF

       CALL convert_steps2unit(io_ev%counter,zztime,newunit,newcount &
            , status, message)
       CALL timer_message(status,substr, message)

       io_ev%unit    = TRIM(newunit)
       io_ev%counter = newcount

       SELECT CASE(io_ev%unit)
       CASE(TIME_INC_SECONDS)
          CALL timer_message(0, ' ' &
               ,'Convert event interval from STEPS into SECONDS')
       CASE(TIME_INC_MINUTES)
          CALL timer_message(0,' '&
               ,'Convert event interval from STEPS into MINUTES')
       END SELECT

    END IF

    ! initialize event with external settings
    CALL event_init(event, ev_name &
         , io_ev%counter, io_ev%unit, io_ev%adjustment, zdtime &
         , io_ev%offset, status, message)
    CALL timer_message(status,substr,message)

    IF (LDEBUG) &
    CALL event_print(event,status, .TRUE., p_parallel_io)

    ! check time stepping against event intervals
    SELECT CASE(io_ev%adjustment)
    CASE(TRIG_FIRST,TRIG_LAST)

       SELECT CASE(io_ev%unit)
       CASE(TIME_INC_SECONDS, TIME_INC_MINUTES, TIME_INC_HOURS, TIME_INC_DAYS)
          IF (event_eval(event,zztime) < 0) THEN
             IF (p_parallel_io .AND. LDEBUG) THEN
                write (*,*) substr, ' io_ev%unit', io_ev%unit
                write (*,*) substr, ' ev%unit', event_unit(event)
             ENDIF
             IF (LDEBUG) &
                  CALL event_print(event, status, .FALSE., p_parallel_io )
             CALL timer_message(3400, substr &
                  , 'Event counter mismatch with time stepping.')
          END IF

       CASE(TIME_INC_MONTHS, TIME_INC_YEARS)
          CALL timer_message(0,' ',&
               'No time stepping mismatch check defined for MONTHS and YEARS.')

       CASE default
          CALL timer_message(3400,substr,'Event unit not defined.')

       END SELECT
    END SELECT

    ! define initial date and trigger dates
    CALL manager_state(BM_time,my_date,INIT_STEP,status)
    CALL timer_message(status,substr)
    IF (p_parallel_io .AND. LDEBUG) &
         CALL write_date(my_date,'my_date: timer_event_init 1')     

!!$    IF (TRIM(eval_date) == EV_TLEV_NEXT) THEN
!!$       ! um_ak_20100525+
!!$       ! Do I really want this shift ???? (eg. 30 min => 0036)
!!$       !CALL add_date(0,idtime,my_date,status)
!!$       !CALL timer_message(status, substr)
!!$       ! um_ak_20100525-
!!$
!!$       CALL event_init(event,my_date,.FALSE.,status,message)
!!$       CALL timer_message(status, substr,message)
!!$       IF (LDEBUG) &
!!$            CALL event_print(event, status, .TRUE., p_parallel_io ) !qqq
!!$       CALL manager_state(BM_time,my_date,INIT_STEP,status)
!!$       CALL timer_message(status,substr)
!!$       CALL event_reinit(event,my_date,status, message)
!!$       CALL timer_message(status,substr, message)
!!$       IF (LDEBUG) &
!!$            CALL event_print(event, status, .TRUE., p_parallel_io ) !qqq
!!$    ELSE IF (TRIM(eval_date) == EV_TLEV_PREV) THEN
!!$       ! um_ak_20100525+
!!$       ! Do I really want this shift ???? (eg. 30 min => 0036)
!!$       !CALL add_date(0,-idtime,my_date,status)
!!$       !CALL timer_message(status, substr)
!!$       ! um_ak_20100525-
!!$       CALL event_init(event,my_date,.FALSE.,status,message)
!!$       CALL timer_message(status, substr,message)
!!$       IF (LDEBUG) &
!!$            CALL event_print(event, status, .TRUE., p_parallel_io ) !qqq
!!$       CALL manager_state(BM_time,my_date,INIT_STEP,status)
!!$       CALL timer_message(status,substr)
!!$       CALL event_reinit(event,my_date,status, message)
!!$       CALL timer_message(status,substr, message)
!!$       IF (LDEBUG) &
!!$            CALL event_print(event, status, .TRUE., p_parallel_io ) !qqq
!!$    ELSE         
       CALL event_init(event,my_date,.FALSE.,status,message)
       CALL timer_message(status,substr, message)
!!$    END IF

    ! find next trigger date starting at initial date 
    IF (lresume) THEN

       !================= preliminar
       ! the finding of the next possible trigger can take a lot of
       ! time if the present date is very fare from the start date
       !
       ! revision needed with a new I/O concept:
       !     event elements must be available in a rerun file

       DO
          CALL event_next_date(event,my_date)
          IF (p_parallel_io .AND. LDEBUG) &
               CALL write_date(my_date, 'next_date 1') 

          SELECT CASE(eval_date)
          CASE(EV_TLEV_PRES)          ! check with current date
             CALL if_less(my_date,current_date,lsm,status)
             CALL timer_message(status,substr)
             IF (.NOT.lsm) EXIT

          CASE(EV_TLEV_NEXT)          ! check with next date
             CALL if_less(my_date,next_date,lsm,status)
             CALL timer_message(status,substr)
             IF (.NOT. lsm) EXIT

          CASE(EV_TLEV_PREV)          ! check with previous date
             CALL if_less(my_date,previous_date,lsm,status)
             CALL timer_message(status,substr)
             IF (.NOT. lsm) EXIT

          END SELECT

          ! next date smaller as current date, rotate all dates
          CALL event_init(event,my_date,.TRUE.,status, message)
          CALL timer_message(status,substr, message)

       END DO
       IF (p_parallel_io .AND. LDEBUG) &
            CALL write_date(my_date, 'next_date 2') 
    END IF

    !print out event settings 
    CALL event_print(event, status, .TRUE., p_parallel_io )

  END SUBROUTINE timer_event_init

  !-----------------------------------------------------------------------------

  SUBROUTINE p_bcast_event (event, p_source, comm) 

    ! distribute events to all nodes
    !
    USE messy_main_mpi_bi,  ONLY: p_bcast

    IMPLICIT NONE
    INTRINSIC :: INT, PRESENT

    TYPE(io_time_event), INTENT(inout) :: event
    INTEGER,             INTENT(in)    :: p_source
    INTEGER, OPTIONAL,   INTENT(in)    :: comm
    ! LOCAL
    INTEGER :: isender, icomm

    isender = INT(p_source)
#ifdef COSMO
    IF (PRESENT(comm)) THEN
       icomm   = INT(comm,i4)
       CALL p_bcast(event%counter,    isender, icomm=icomm)
       CALL p_bcast(event%unit,       isender, icomm=icomm)
       CALL p_bcast(event%adjustment, isender, icomm=icomm)
       CALL p_bcast(event%offset,     isender, icomm=icomm)
    ELSE
#endif
       CALL p_bcast(event%counter,    isender)
       CALL p_bcast(event%unit,       isender)
       CALL p_bcast(event%adjustment, isender)
       CALL p_bcast(event%offset,     isender)
#ifdef COSMO
    ENDIF
#endif
  END SUBROUTINE p_bcast_event

  ! ---------------------------------------------------------------------------

  FUNCTION I_event_steps (event, delta) RESULT (ix)  !************************

    ! get out interval between two event trigger points in steps
    !
    IMPLICIT NONE
    INTRINSIC :: INT, REAL

    TYPE (time_event), INTENT(in) :: event
    REAL(dp),          INTENT(in) :: delta   ! seconds of special unit
    INTEGER                       :: ix
    INTEGER                       :: isecs

    isecs = I_event_seconds(event)
    ix = INT(REAL(isecs,dp)+0.0001_dp/delta) ! convert into steps

  END FUNCTION I_event_steps

  ! ----------------------------------------------------------------------------

  FUNCTION I_event_seconds (event) RESULT (ix) !********************************

    ! get out interval between two event trigges in seconds

    IMPLICIT NONE
    INTRINSIC :: INT

    TYPE (time_event), INTENT(in) :: event
    INTEGER                       :: ix

    !LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='main_timer event_seconds'
    INTEGER           :: iday, isec
    TYPE (time_days)  :: my_date, previous_date
    INTEGER           :: status

    ! get difference between last and present trigger date
    CALL event_previous_date(event, previous_date)
    CALL date_get(previous_date, iday, isec, status)
    CALL timer_message(status,substr)

    iday = - iday
    isec = - isec

    CALL event_current_date(event, my_date)
    CALL add_date (iday, isec, my_date, status)
    CALL timer_message(status,substr)

    ! check result
    CALL date_get(my_date,iday,isec,status)
    CALL timer_message(status,substr)

    IF (iday < 0 .OR. isec < 0) &
         CALL timer_message(3400,substr,'Event interval < 0')

    ! convert into seconds
    ix = iday*INT(OneDay) + isec !+ 0.0001_dp

    ! at initial time the interval is zero
    ! this (may be) is important for accumulation

  END FUNCTION I_event_seconds

  ! ----------------------------------------------------------------------------

  FUNCTION I_event_seconds_next (event, lnext) RESULT (ix)  !*******************

    ! return distance between present and next trigger
    !
    IMPLICIT NONE
    INTRINSIC :: INT

    TYPE (time_event), INTENT(in) :: event
    LOGICAL,           INTENT(in) :: lnext
    INTEGER                       :: ix

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='main_timer event_seconds_next'
    INTEGER           :: iday, isec
    TYPE (time_days)  :: my_date, current_date
    INTEGER           :: status

    ix = 0
    IF (lnext) THEN
       CALL event_current_date(event, current_date)
       CALL date_get(current_date,iday,isec,status)
       CALL timer_message(status,substr)

       CALL event_next_date(event, my_date)
       CALL add_date (-iday, -isec, my_date,status)
       CALL timer_message(status,substr)

       CALL date_get(my_date,iday,isec,status)
       CALL timer_message(status,substr)

       IF (iday < 0 .OR. isec < 0) &
            CALL timer_message(3400, substr, 'Event interval < 0')

       ! convert into seconds
       ix = iday*INT(OneDay) + isec !+ 0.0001_dp)

    END IF

  END FUNCTION I_event_seconds_next

  ! ----------------------------------------------------------------------------
  !
  LOGICAL FUNCTION L_event_trigger (event, date) 
    !
    IMPLICIT NONE
    INTRINSIC :: INT, REAL, TRIM

    TYPE (time_event), INTENT(inout) :: event
    TYPE (time_days),  INTENT(in)    :: date

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='main_timer L_event_trigger'
    REAL(dp)         :: zsec  ! length of event increment in seconds
    INTEGER          :: iday1, isec1, iday2, isec2
    INTEGER          :: status
    TYPE (time_days) :: adj_start, adj_stop
    LOGICAL          :: lequal1, lequal2, learly1, learly2
    CHARACTER(LEN=STRLEN_ULONG) :: m_text

    !IF (p_parallel_io) CALL event_print(event, status, .TRUE., p_parallel_io )
    L_event_trigger = .FALSE.
    m_text = ' '

    IF (.NOT. event_is_active(event))  RETURN
    ! check if the date will fit into the trigger interval
    zsec = event_delta(event)
    IF (TRIM(event_adjust(event)) == TRIM(TRIG_EXACT)) &
         zsec = event_half_delta(event)
    iday1 = INT((zsec+0.0001_dp)/OneDay)
    isec1 = INT(zsec - REAL(iday1,dp)*OneDay + 0.0001_dp)

    IF (TRIM(event_adjust(event)) == TRIM(TRIG_EXACT)) THEN
       zsec  = event_half_delta(event)-1.0_dp
       iday2 = INT((zsec+0.0001_dp)/OneDay)
       isec2 = INT(zsec - REAL(iday2,dp)*OneDay + 0.0001_dp)
    END IF

    ! detect first or last second in given unit 
    CALL event_next_date(event, adj_start)
    CALL event_next_date(event, adj_stop)

    SELECT CASE(TRIM(event_adjust(event)))
    CASE(TRIM(TRIG_FIRST))
       CALL add_date( iday1, isec1, adj_stop,status)
       CALL timer_message(status,substr)
    CASE(TRIM(TRIG_LAST))
       CALL add_date(-iday1,-isec1, adj_start,status)
       CALL timer_message(status,substr)
    CASE(TRIM(TRIG_EXACT))
       CALL add_date( iday1, isec1, adj_stop,status)
       CALL timer_message(status,substr)
       CALL add_date(-iday2,-isec2, adj_start,status)
       CALL timer_message(status,substr)
    END SELECT

    IF (p_parallel_io .AND. LDEBUG) THEN
       CALL write_date(date,      ' date1')
       CALL write_date(adj_start, ' date2')
       CALL write_date(adj_stop,  ' date3')
    ENDIF

    ! check the position of the present date 
    CALL if_less(adj_stop, date, learly1,status)
    CALL timer_message(status,substr)

    IF (learly1) THEN
       m_text = 'Warning: event trigger not longer in future <'&
            // TRIM(event_label(event)) // '>'
       CALL timer_message(0,substr,m_text)

    ELSE    ! adj_start <= my_date <= adj_stop
       CALL if_equal(adj_start,date,lequal1,status)
       CALL timer_message(status,substr)
       CALL if_equal(adj_stop,date,lequal2,status)
       CALL timer_message(status,substr)
       CALL if_less(adj_start,date,learly1,status)
       CALL timer_message(status,substr)
       CALL if_less(date,adj_stop,learly2,status)
       CALL timer_message(status,substr)

       L_event_trigger = lequal1 .OR. lequal2 .OR. (learly1 .AND. learly2) 

       IF (L_event_trigger)  THEN
          CALL event_init (event, date, .TRUE. , status, m_text)
          CALL timer_message(status,substr,m_text)
       END IF

    ENDIF

  END FUNCTION L_event_trigger

  ! ----------------------------------------------------------------------------
  !
  FUNCTION event_eval (event, delta) RESULT (slen) 

    ! check the multiple of delta in event interval
    !
    IMPLICIT NONE
    INTRINSIC :: INT, MOD, PRESENT, TRIM

    TYPE(time_event),   INTENT(in) :: event
    REAL(dp), OPTIONAL, INTENT(in) :: delta

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='main_timer event_eval'
    INTEGER                     :: slen
    INTEGER                     :: fit
    CHARACTER(LEN=STRLEN_ULONG) :: m_text

    slen = -1

    IF (.NOT. event_is_init(event)) THEN
       WRITE (*,*) ('event not initialised, no unit given.')

    ELSE
       SELECT CASE(event_unit(event))
       CASE(TIME_INC_SECONDS);   slen = event_count(event)
       CASE(TIME_INC_MINUTES);   slen = event_count(event) * 60
       CASE(TIME_INC_HOURS);     slen = event_count(event) * 60 * 60
       CASE(TIME_INC_DAYS);      slen = event_count(event) * 60 * 60 * 24
       CASE(TIME_INC_MONTHS, TIME_INC_YEARS)
          WRITE (*,*) substr,'Exact value undefined.'
          WRITE (*,*) substr &
               ,'Event trigger interval for months or years may varied.'
       CASE default
          m_text = 'Counter unit unknown ::'//TRIM(event_unit(event))
          CALL timer_message(3400, TRIM(substr), TRIM(m_text))
       END SELECT

       IF (slen > 0 .AND. PRESENT(delta)) THEN
          fit = MOD(slen,INT(delta))
          IF (fit /= 0) THEN
             WRITE(m_text,*) 'Event <',TRIM(event_label(event))&
                  ,   '> interval not a multiple of ',INT(delta)
             CALL timer_message(1,substr,m_text)
             slen = -fit
          END IF
       END IF

    END IF

  END FUNCTION event_eval

  ! ----------------------------------------------------------------------------

  SUBROUTINE event_print (event, status, short, lwrite)

    ! print event contents
    ! CALL for IO PE only

    IMPLICIT NONE
    INTRINSIC :: INT, NINT, PRESENT, TRIM

    TYPE (time_event), INTENT(in) :: event
    INTEGER,           INTENT(OUT):: status
    LOGICAL, OPTIONAL, INTENT(in) :: short
    LOGICAL, OPTIONAL, INTENT(IN) :: lwrite

    !LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='main_timer: event_print'
    TYPE (time_days)            :: previous_date, current_date, next_date
    LOGICAL                     :: lequal, llshort
    CHARACTER(LEN=STRLEN_ULONG) :: date_text

    status = 0    
    IF (PRESENT(lwrite)) THEN
       IF (.NOT. lwrite) RETURN
    ENDIF

    IF (event_is_init(event)) THEN

       IF (PRESENT(SHORT)) THEN
          llshort = .FALSE.
       ELSE
          llshort = short
       ENDIF
       IF (llshort) THEN
          WRITE(*,*) 'Event <',TRIM(event_label(event)),&
               '> : interval ',event_count(event),' '&
               ,TRIM(event_unit(event)),&
               ' : adjustment ',TRIM(event_adjust(event)) &
               ,' : offset[sec] ',event_offset(event)
       ELSE
          IF (.NOT. event_is_active(event)) THEN
             WRITE(*,*) 'Event <',TRIM(event_label(event)),'> ... not active'
          ELSE
             WRITE(*,*) ' '
             WRITE(*,*) 'State of event <',TRIM(event_label(event)) &
                  ,'> ... initialized'

             SELECT CASE(event_unit(event))
             CASE(TIME_INC_SECONDS)
                WRITE(*,*) ' trigger each ',event_count(event),' seconds'

             CASE(TIME_INC_MINUTES)
                WRITE(*,*) ' trigger each ',event_count(event),' minutes'

             CASE(TIME_INC_HOURS)
                WRITE(*,*) ' trigger each ',event_count(event),' hours'

             CASE(TIME_INC_DAYS)
                WRITE(*,*) ' trigger each ',event_count(event),' days'

             CASE(TIME_INC_MONTHS)
                WRITE(*,*) ' trigger each ',event_count(event),' months'

             CASE(TIME_INC_YEARS)
                WRITE(*,*) ' trigger each ',event_count(event),' years'

             CASE default
                WRITE(*,*) 'Counter unit unknown ::',event_unit(event)
                STATUS = 3400
                RETURN
             END SELECT

             CALL event_previous_date(event, previous_date)
             CALL print_date_components(previous_date, status &
                  , mess=date_text)
             IF (status /= 0) RETURN
             WRITE (*,*) '  last event trigger date is ...    ' // TRIM(date_text)

             CALL event_current_date(event, current_date)
             CALL print_date_components(current_date, status, mess=date_text)
             IF (status /= 0) RETURN
             CALL if_equal(current_date,previous_date,lequal, status)
             IF (status /= 0) RETURN
             IF (lequal) THEN
                WRITE(*,*) '  initial event date is ...    ' // TRIM(date_text)
             ELSE
                WRITE(*,*) '  time between two triggers: ',&
                     I_event_seconds(event),' seconds'
                WRITE(*,*) '  present trigger date is ...  ' // TRIM(date_text)
             END IF

             CALL event_next_date(event, next_date)
             CALL print_date_components(next_date,status, mess=date_text)
             IF (status /=0 ) RETURN
             WRITE(*,*)    '  next trigger date is ...         ',&
                  TRIM(date_text),' (offset of ',event_offset(event) &
                  ,' seconds included)'

             SELECT CASE(TRIM(event_adjust(event)))
             CASE(TRIM(TRIG_FIRST))
                WRITE(*,*) '  adjustment at beginning of unit,',&
                     ' interval: ( trigger : trigger + '&
                     , NINT(event_delta(event)),'s )'

             CASE(TRIM(TRIG_LAST))
                WRITE(*,*) '  adjustment at end of unit,',&
                     ' interval: ( trigger - ', NINT(event_delta(event))&
                     ,'s : trigger )'

             CASE(TRIM(TRIG_EXACT))
                WRITE(*,*) '  no adjustment,',&
                     ' interval: ( trigger - ',INT(event_half_delta(event)),&
                     's : trigger + ',INT(event_half_delta(event))-1,'s )'

             END SELECT

          END IF

       END IF

    ELSE 
       STATUS = 1
       WRITE (*,*)  substr,'no printout of uninitialized event'

    END IF


  END SUBROUTINE Event_print


  ! ------------------------------------------------------------------------
  ! ------------------------------------------------------------------------
  ! PRIVATE ROUTINES
  ! ------------------------------------------------------------------------
  ! ------------------------------------------------------------------------

  SUBROUTINE init_BM_time

    ! initialize BASE model dates
    ! CALLED FROM messy_timer_initialize
    IMPLICIT NONE
    INTRINSIC :: ABS

    INTEGER           :: istep, status
    INTEGER           :: yr, mo, dy, hr, mi, se
    REAL(dp)          :: zdtold
    LOGICAL           :: lsmaller, lequal
    TYPE(time_days)   :: date
    CHARACTER(LEN=STRLEN_ULONG) :: m_text
    CHARACTER(LEN=*), PARAMETER :: substr='main_timer: init_BM_time'

    CALL manager_state(BM_time,zdtold,status)
    CALL timer_message(status,substr)

    IF (ABS(zdtold-delta_time) > 0.0_dp) THEN

       IF (lresume) THEN
          CALL info_bi('Attention: length of timestep has been changed!', ' ')
          WRITE (m_text,*) ' Old timestep was ', zdtold, ' seconds'
          CALL info_bi(m_text, ' ')
          WRITE (m_text,*) ' New timestep is  ', delta_time, ' seconds'
          CALL info_bi(m_text,' ')
       END IF
       CALL manager_init(BM_time,delta_time,status)
       CALL timer_message(status,substr)
       ! correct to old start_date (was changed by change of delta_time)
       CALL manager_init(BM_time,start_date,status)
       CALL timer_message(status,substr)
       ! Correct the start_date components
       CALL timer_get_date(status,start_date, yr, mo, dy, hr, mi, se)
       CALL timer_message(status,substr)
    END IF

    IF (.NOT. lresume) THEN

! um_ak_20090925+
! This should not be needed anymore
!!$#if defined(ECHAM5) || defined(MPIOM)
!!$       ! in ECHAM5 the manager is initialize by an artificial date in the
!!$       ! first place, which needs to be replaced by the namelist date here
!!$       ! get start_date from BASE model
!!$       !
!!$       CALL manager_init(BM_time,start_date,status, m_text)
!!$       CALL timer_message(status,substr, m_text)
!!$       !
!!$       IF (p_parallel_io) CALL write_date(start_date &
!!$            ,'Start date replaced by namelist start date: ')
!!$       !
!!$#endif
! um_ak_20090925-
       ! set the time manager to the init_step now
       ! reset time step at beginning of experiment
       ! count backward
       ! qqq needed for COSMO ?
       istep = INIT_STEP - timer_get_time_step()
       CALL manager_init(BM_time,istep,status)
       CALL timer_message(status,substr)
    ELSE

       CALL manager_state(BM_time, date,status)
       CALL timer_message(status,substr)

       CALL if_less(date,start_date,lsmaller,status)
       CALL timer_message(status,substr)
       CALL if_equal(date,start_date,lequal,status)
       CALL timer_message(status,substr)
       IF (lsmaller) THEN
          CALL write_date(date,'Current date ...')
          CALL write_date(start_date,'Start date ...')
          CALL timer_message (3400, substr,'Start date in future')

       ELSE IF (lequal) THEN
          lresume = .FALSE.
          CALL timer_message(0,substr &
               ,'Set start date to resume date, force initial run')
       END IF

    END IF

    ! preliminary initializion of  previous date and next date
    istep = timer_get_time_step()
    CALL manager_state(BM_time,previous_date,istep-1, status)
    CALL timer_message(status, substr)
    CALL manager_state(BM_time,next_date,    istep+1, status)
    CALL timer_message(status, substr)

    CALL timer_message(0, ' ','Time step and start date evaluation done.') 

    IF (lfirst_cycle)  CALL manager_init(BM_time,.TRUE.,status,m_text) 
    CALL timer_message(status, substr,m_text)

    IF (p_parallel_io) THEN
       CALL manager_print(BM_time, status)
       CALL timer_message(status,substr)
    ENDIF

    ! define time stepping 
    istep  = timer_get_time_step()
    lstart = (istep == INIT_STEP) 
    !IF (lstart) THEN
    time_step_len = delta_time
    !ELSE
    ! will be called from data_global_start
    !   CALL timer_set_time_step_len(l2tls)
    !END IF

    ! start date of manager can be changed final definition of start date here
    CALL manager_state(BM_time,current_date, status)
    CALL timer_message(status,substr)

    CALL manager_state(BM_time,start_date,INIT_STEP,status)
    CALL timer_message(status,substr)

    ! stop of experiment evaluated only during first rerun cycle
    ! evaluation in the following order (highest priority left)
    ! NO_STEPS or NO_DAYS or DT_STOP or default

    IF (.NOT. lfirst_cycle) THEN
       CALL timer_message(0,' ','No MESSy evaluation of model stop date.') 
    END IF

    ! check date order
    CALL if_less(start_date,stop_date, lsmaller,status)
    CALL timer_message(status,substr)
    IF (.NOT. (lsmaller)) THEN
       CALL timer_message(3441, substr)
    ELSE
       CALL if_less(stop_date,current_date, lsmaller,status)
       CALL timer_message(status,substr)
       IF (lsmaller) CALL timer_message(3442, substr)
    END IF

    CALL write_date(stop_date,'Stop experiment at: ')

    ! initialize previous date and next date
    CALL manager_state(BM_time,previous_date,istep-1,status)
    CALL timer_message(status,substr)
    CALL write_date   (previous_date,'Previous date: ')

    CALL manager_state(BM_time,next_date,    istep+1,status)
    CALL timer_message(status,substr)
    CALL write_date   (next_date,    'Next date    : ')

    CALL timer_message (0,' ',' ')

  END SUBROUTINE init_BM_time
  ! ---------------------------------------------------------------------------

  ! ---------------------------------------------------------------------------
#ifdef COSMO
  SUBROUTINE check_rerun_event(event)

    IMPLICIT NONE

    TYPE(time_event), INTENT(INOUT) :: event
    !LOCAL
    LOGICAL :: le
    CHARACTER(LEN=*), PARAMETER :: substr = 'check_rerun_event'
    INTEGER                     :: status
    CHARACTER(LEN=STRLEN_ULONG) :: message

    status = 0
    IF (event_is_active(event)) THEN

       ! CHECK if restart stop date is reached 
       CALL if_less(rerun_stop_date, current_date, le, status)
       CALL timer_message(status, substr)

       ! DEACTIVE EVENT
       IF (le) THEN
          CALL event_init(event, TRIG_NONE, status, message)
          CALL timer_message(status, substr, message)
       ENDIF
    ENDIF

  END SUBROUTINE check_rerun_event
#endif
  ! ---------------------------------------------------------------------------

  ! -------------------------------------------------------------------------
  SUBROUTINE write_date(date, text) 

    ! write date in constant format to standard output
    ! input date can be different declared
    !-
    IMPLICIT NONE
    INTRINSIC :: PRESENT, TRIM

    TYPE (time_days), INTENT(in)           :: date
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: text

    ! LOCAL
    CHARACTER(len=STRLEN_ULONG) :: date_mess1, date_mess2
    INTEGER                     :: status

    IF (p_parallel_io) THEN
       CALL print_date_components(date,status,date_mess1)
       CALL timer_message(status, ' ')

       IF (PRESENT(text)) THEN
          date_mess2 = TRIM(text)//' '//TRIM(date_mess1)
       ELSE
          date_mess2 = TRIM(date_mess1)
       END IF
       CALL timer_message(0,' ',date_mess2)
    ENDIF

  END SUBROUTINE write_date
  ! --------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------
  ! PUBLIC HELPER ROUTINES
  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------

  ! ---------------------------------------------------------------------------
  SUBROUTINE timer_set_rerun_event(status, counter, unit, adjustment, offset &
       , rerun_stop)

    IMPLICIT NONE

    ! I/O
    INTEGER,  INTENT(OUT)        :: status
    INTEGER,  INTENT(IN)         :: counter    ! No. of steps in given unit
    CHARACTER(LEN=*), INTENT(IN) :: unit       ! counter unit type
    CHARACTER(LEN=*), INTENT(IN) :: adjustment ! adjustment in side the unit
    ! offset to initial date in seconds
    INTEGER, INTENT(IN)                     :: offset
    ! 
    ! number of hours after which rerun output is switched off
    INTEGER,  INTENT(IN), OPTIONAL          :: rerun_stop

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'timer_set_rerun_event'

    IO_RERUN_EV%counter    = counter 
    IO_RERUN_EV%unit       = unit
    IO_RERUN_EV%adjustment = adjustment
    IO_RERUN_EV%offset     = offset

    status = 0
    RETURN

  END SUBROUTINE timer_set_rerun_event
  ! ---------------------------------------------------------------------------

  ! ---------------------------------------------------------------------------
  SUBROUTINE timer_get_rerun_event(status, counter, unit, adjustment, offset &
       , ev_adjust)

    IMPLICIT NONE

    ! I/O
    INTEGER,  INTENT(OUT)         :: status
    INTEGER,  INTENT(OUT)         :: counter    ! No. of steps in given unit
    CHARACTER(LEN=*), INTENT(OUT) :: unit       ! counter unit type
    CHARACTER(LEN=*), INTENT(OUT) :: adjustment ! adjustment in side the unit
    ! offset to initial date in seconds
    INTEGER, INTENT(OUT)          :: offset
    ! 
    CHARACTER(LEN=*), INTENT(OUT) :: ev_adjust

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'timer_get_rerun_event'

    status = 1 

    counter    = event_count(RERUN_EVENT)
    unit       = event_unit(RERUN_EVENT)
    adjustment = event_adjust(RERUN_EVENT)
    offset     = event_offset(RERUN_EVENT)

    ev_adjust = 'next'

    status = 0

  END SUBROUTINE timer_get_rerun_event
  ! ---------------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE main_timer_read_nml_cpl(status, iou)

    ! MODULE ROUTINE (SMIL)
    !
    ! READ NAMELIST, CHECK IT, AND INITIALIZE GLOBAL VARIABLES
    !
    ! Author: Astrid Kerkweg, UNI-MZ, Jun 2009

    USE messy_main_tools,   ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE
    INTRINSIC :: TRIM

    ! I/O
    INTEGER, INTENT(OUT) :: status ! error status
    INTEGER, INTENT(IN)  :: iou    ! logical I/O unit

    NAMELIST /CPL/ CAL_TYPE, MODEL_START, MODEL_STOP, IO_RERUN_EV &
                 !, IO_TEST_EV &
                 , delta_time, lresume, NO_CYCLES, LABORT &
                 , NO_DAYS, NO_STEPS  ! um_ak_20120301

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_timer_read_nml_cpl'
    LOGICAL                     :: lex          ! file exists ?
    INTEGER                     :: fstat        ! file status

    ! INITIALIZE
    status = 1 ! ERROR

    MODEL_START(:)= 0
    MODEL_STOP(:) = 0

    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist
    
    READ(iou, NML=CPL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    ! CHECK NAMELIST
    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! no ERROR

  END SUBROUTINE main_timer_read_nml_cpl
  ! -------------------------------------------------------------------

  ! ---------------------------------------------------------------------------
  SUBROUTINE timer_message(status,substr, message)

    IMPLICIT NONE
    INTRINSIC :: PRESENT, TRIM

    INTEGER, INTENT(IN)          :: status
    CHARACTER(LEN=*), INTENT(IN) :: substr
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: message

    ! LOCAL
    CHARACTER(LEN=STRLEN_ULONG)  :: errorstr
    CHARACTER(LEN=*), PARAMETER  :: m_text= ' '

    SELECT CASE (status)
    CASE (0)
       ! no error
       IF (PRESENT(message)) THEN
          IF (.NOT. TRIM(message) == TRIM(m_text)) &
               CALL info_bi(TRIM(message), substr)
       ENDIF
       RETURN
    CASE (1)
       ! WARNING
       IF (PRESENT(message)) CALL warning_bi(TRIM(message),substr)
       RETURN
    CASE(3400)
       IF (PRESENT(message)) THEN
          CALL error_bi(TRIM(message), substr)
       ELSE
          CALL error_bi('ERROR 3400' , substr)
       ENDIF
       RETURN
    CASE(3434)
       errorstr = 'manager not initialized'
    CASE(3435)
       errorstr = 'date not initialized'
    CASE(3436)
       errorstr = 'messy_main_timer manager_reinit_days//&
            &// Time manager was not initialized, no reinit possible.'
    CASE(3437)
       errorstr = 'messy_main_timer:manager_step negative steps invalid' 
    CASE(3438)
       errorstr = 'new date before start date'
    CASE(3439)
       errorstr = 'new date do not fit to a time step'  
    CASE(3440)
       errorstr = 'combination of offset/unit invalid' 
    CASE(3441)
       errorstr = 'Start date larger/equal than stop date ....'
    CASE(3442)
       errorstr = 'Current date larger than stop date ....'
    CASE(3443)
       errorstr = 'Unknown date ....'
    CASE(3444)
       errorstr = 'Unknown calendar type ....'
    CASE(3445)
       errorstr = 'Time difference between start and resume date negative...'
    CASE(3446)
       errorstr = 'Start date in future'  
    CASE(3447)
       errorstr = 'Set start date to resume date, force initial run'
    CASE(3450)
       errorstr = 'Unit months not unambiguously convertable to seconds'
    CASE(3451)
       errorstr = 'Unit years not unambiguously convertable to seconds'
    CASE DEFAULT
       errorstr = ' '
    END SELECT

    CALL error_bi(errorstr, substr)

  END SUBROUTINE timer_message
  ! ----------------------------------------------------------------------------
END MODULE messy_main_timer_bi
