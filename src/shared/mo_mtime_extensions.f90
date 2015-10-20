!> Additional C-bindings for the mtime library.
!!
!! @author F. Prill, DWD
!!
!! @par Revision History
!! Initial implementation  by  F. Prill, DWD (2013-09-17)
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
MODULE mo_mtime_extensions

  USE, INTRINSIC :: iso_c_binding, ONLY: c_char, c_ptr, c_f_pointer, c_associated, c_bool, &
    &               c_loc, c_int64_t, c_null_char, c_int, c_double
  USE mo_kind,     ONLY: wp, i8
  USE mtime,       ONLY: datetime, event, MAX_TIMEDELTA_STR_LEN, MAX_DATETIME_STR_LEN, &
    &                    newDatetime, deallocateDatetime, datetimeToString, timedelta, &
    &                    newTimedelta, deallocateTimedelta, OPERATOR(+), OPERATOR(-),  &
    &                    getPTStringFromMS,                                            &
    &                    divisionquotienttimespan, OPERATOR(==),                       &
    &                    divideTimeDeltaInSeconds, OPERATOR(<)
  USE mo_datetime, ONLY: t_datetime
  USE mo_exception,ONLY: message,finish
  USE mo_util_string, ONLY: int2string, real2string
  IMPLICIT NONE
  
  PRIVATE
!  PUBLIC :: getPTStringFromMS
!  PUBLIC :: getTimeDeltaFromDateTime
!  PUBLIC :: get_duration_string
!  PUBLIC :: get_duration_string_real
   PUBLIC :: get_datetime_string
   PUBLIC :: convert_datetime_string_old2new
   PUBLIC :: datetime_str_equal
   PUBLIC :: timedelta_str_equal
   PUBLIC :: get_timedelta_divide_by_seconds
   PUBLIC :: getTimedeltaFromMS
!  PUBLIC :: getTriggeredPreviousEventAtDateTime

  !> module name
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_mtime_extensions'


  INTERFACE

!    SUBROUTINE my_getptstringfromms(ms, ptstr) BIND(c, name='getPTStringFromMS')
!#ifdef __SX__
!      USE, INTRINSIC :: iso_c_binding, ONLY: c_int64_t, c_ptr, c_char
!#else
!      import :: c_int64_t, c_ptr, c_char
!#endif
!      INTEGER(c_int64_t), VALUE :: ms
!      CHARACTER(c_char), DIMENSION(*) :: ptstr
!    END SUBROUTINE my_getptstringfromms

!    SUBROUTINE my_gettimedeltafromdatetime(dt1, dt2, td_return) BIND(c, name='getTimeDeltaFromDateTime')
!#ifdef __SX__
!      USE, INTRINSIC :: iso_c_binding, ONLY: c_ptr
!#else
!      IMPORT :: c_ptr
!#endif
!      TYPE(c_ptr), value :: dt1, dt2
!      TYPE(c_ptr), value :: td_return
!    END SUBROUTINE my_gettimedeltafromdatetime

!    SUBROUTINE my_gettriggeredpreviouseventatdatetime(my_event, my_datetime) BIND(c, name='getTriggeredPreviousEventAtDateTime')
!#ifdef __SX__
!      USE, INTRINSIC :: iso_c_binding, ONLY: c_ptr
!#else
!      IMPORT :: c_ptr
!#endif
!      TYPE(c_ptr), value :: my_event
!      TYPE(c_ptr), value :: my_datetime
!    END SUBROUTINE my_gettriggeredpreviouseventatdatetime
  END INTERFACE


  INTERFACE get_datetime_string
    MODULE PROCEDURE get_datetime_string
    MODULE PROCEDURE get_datetime_string_str
  END INTERFACE

!  INTERFACE getPTStringFromMS
!    MODULE PROCEDURE getPTStringFromMS_int
!    MODULE PROCEDURE getPTStringFromMS_real
!  END INTERFACE

CONTAINS 
  
  ! SUBROUTINE getPTStringFromMS_int(ms, string)
  !   INTEGER, INTENT(in) :: ms
  !   CHARACTER(len=max_timedelta_str_len), INTENT(out) :: string
  !   INTEGER :: i
  !   INTEGER(c_int64_t) :: ms_int64

  !   ms_int64 = INT(ms,c_int64_t)
  !   CALL my_getptstringfromms(ms_int64, string)
  !   char_loop: DO i = 1 , LEN(string)
  !     IF (string(i:i) == c_null_char) EXIT char_loop
  !   END DO char_loop
  !   string(i:LEN(string)) = ' '
  ! END SUBROUTINE getPTStringFromMS_int

  ! SUBROUTINE getPTStringFromMS_real(ms, string)
  !   REAL, INTENT(in) :: ms
  !   CHARACTER(len=max_timedelta_str_len), INTENT(out) :: string
  !   CALL getPTStringFromMS_int(NINT(ms), string)
  ! END SUBROUTINE getPTStringFromMS_real

  ! SUBROUTINE getTimeDeltaFromDateTime(dt1, dt2, td_return)
  !   TYPE(datetime),  INTENT(IN),    TARGET   :: dt1,dt2 
  !   TYPE(timedelta), INTENT(INOUT), TARGET   :: td_return     !< OUT

  !   CALL my_gettimedeltafromdatetime(C_LOC(dt1), C_LOC(dt2), C_LOC(td_return))
  ! END SUBROUTINE getTimeDeltaFromDateTime

!  SUBROUTINE getTriggeredPreviousEventAtDateTime(my_event, my_datetime)
!    TYPE(event)   , TARGET                ::  my_event
!    TYPE(datetime), TARGET, INTENT(INOUT) ::  my_datetime    !< OUT

!   CALL my_gettriggeredpreviouseventatdatetime(C_LOC(my_event), C_LOC(my_datetime))
! END SUBROUTINE getTriggeredPreviousEventAtDateTime


  !> Conversion of a time stamp string which is suitable for a
  !> t_datetime object to an ISO8601 conforming time stamp string.
  !
  !  This function is required for backward-compatibility of
  !  t_datetime-based namelist parameters! It has been derived from
  !  the subroutine mo_datetime::string_to_datetime.
  !
  !  Note: The "old" t_datetime time stamp strings were in fact
  !        ISO8601 conforming. This subroutine exists only to ake sure
  !        that the proper delimiters etc. are used.
  !
  FUNCTION convert_datetime_string_old2new(str_oldfmt) RESULT(ret)
    CHARACTER(LEN=MAX_DATETIME_STR_LEN)  :: ret          !  Result: ISO8601 conforming time stamp string.
    CHARACTER(LEN=*), INTENT(IN)         :: str_oldfmt   !  Time stamp string suitable for t_datetime
    ! local variables
    CHARACTER(len=32) :: datetime_format
    CHARACTER(len=16) :: year_fmt_str
    INTEGER           :: position_of_first_dash, position_of_second_dash, &
      &                  length, year_format
    INTEGER           :: second = 0
    !
    INTEGER     :: dt_year      ! [yr]
    INTEGER     :: dt_month     ! [mo]  1:12
    INTEGER     :: dt_day       ! [dy]  1:{28,29,30,31}
    INTEGER     :: dt_hour      ! [hr]  0:23
    INTEGER     :: dt_minute    ! [mn]  0:59
    REAL(wp)    :: dt_second    ! [s]   [0.,60.[

    ! --- build format string
    length = LEN_TRIM(str_oldfmt)
    position_of_first_dash = INDEX ( str_oldfmt, '-', .FALSE. )
    IF ( position_of_first_dash == 1 ) THEN
       ! we have negative years 
       position_of_second_dash = INDEX ( str_oldfmt(2:length),  '-', .FALSE. )
       year_format = position_of_second_dash - 1
       WRITE (year_fmt_str, '(A1,I1)' )  'I', year_format + 1
    ELSE
       year_format = position_of_first_dash - 1
       WRITE (year_fmt_str, '(A1,I1)' )  'I', year_format
    ENDIF
    WRITE (datetime_format,'(A1,A2,A13)') '(',year_fmt_str,',5(1X,I2),I2)'

    ! --- read from string
    READ(str_oldfmt, TRIM(datetime_format)) dt_year,   &
      &         dt_month, dt_day, dt_hour, dt_minute, second
    dt_second = REAL(second,wp)

    ! --- build an ISO8601 string
    ret = "P"//TRIM(int2string(dt_year,'(i0)'))//"-"//TRIM(int2string(dt_month,'(i0)'))//"-"//&
      &TRIM(int2string(dt_day,'(i0)'))//"T"//TRIM(int2string(dt_hour))//":"//&
      &TRIM(int2string(dt_minute))//":"//TRIM(real2string(dt_second,'(5.2f)'))//"Z"
  END FUNCTION convert_datetime_string_old2new


  !> Test, if two given ISO8601 date-time stamp strings specify the
  !  same date.
  FUNCTION datetime_str_equal(datetime_str1, datetime_str2)
    LOGICAL :: datetime_str_equal
    CHARACTER(LEN=*), INTENT(IN) :: datetime_str1, datetime_str2
    ! local variables
    TYPE(datetime), POINTER :: mtime_datetime1, mtime_datetime2

    mtime_datetime1 => newDatetime(datetime_str1)
    mtime_datetime2 => newDatetime(datetime_str2)
    datetime_str_equal = (mtime_datetime1 == mtime_datetime2)
    CALL deallocateDatetime(mtime_datetime1)
    CALL deallocateDatetime(mtime_datetime2)
  END FUNCTION datetime_str_equal


  !> Test, if two given ISO8601 date-time stamp strings specify the
  !  same date.
  FUNCTION timedelta_str_equal(timedelta_str1, timedelta_str2)
    LOGICAL :: timedelta_str_equal
    CHARACTER(LEN=*), INTENT(IN) :: timedelta_str1, timedelta_str2
    ! local variables
    TYPE(timedelta), POINTER :: mtime_timedelta1, mtime_timedelta2

    mtime_timedelta1 => newTimedelta(timedelta_str1)
    mtime_timedelta2 => newTimedelta(timedelta_str2)

    !    timedelta_str_equal = (mtime_timedelta1 == mtime_timedelta2)
    !
    !    workaround:
    timedelta_str_equal = (.NOT. (mtime_timedelta1 < mtime_timedelta2)) .AND. &
      &                   (.NOT. (mtime_timedelta2 < mtime_timedelta1))
    CALL deallocateTimedelta(mtime_timedelta1)
    CALL deallocateTimedelta(mtime_timedelta2)
  END FUNCTION timedelta_str_equal


  !> Divide a time span, given by two ISO8601 date-time stamp strings, by 
  FUNCTION get_timedelta_divide_by_seconds(start_datetime_string, &
    &                                      stop_datetime_string,  &
    &                                      dtime)  RESULT(ret)
    CHARACTER(LEN=*), INTENT(IN) :: start_datetime_string, stop_datetime_string
    REAL(wp),         INTENT(IN) :: dtime
    INTEGER :: ret
    ! local variables
    TYPE(datetime),  POINTER             :: mtime_start, mtime_stop
    TYPE(timedelta), POINTER             :: mtime_dtime, mtime_td
    CHARACTER(LEN=MAX_TIMEDELTA_STR_LEN) :: td_string
    TYPE(divisionquotienttimespan)       :: mtime_quotient

    ! User did not specified a value, we need to compute "nsteps" as
    ! (stop date - start date)/dtime:
    mtime_start => newDatetime(start_datetime_string)
    mtime_stop  => newDatetime(stop_datetime_string)
    mtime_td    => newTimedelta("PT0S")
    mtime_td    =  mtime_stop - mtime_start
    CALL getPTStringFromMS(INT(dtime*1000._wp,i8), td_string)
    mtime_dtime => newTimedelta(td_string)
    CALL divideTimeDeltaInSeconds(mtime_td, mtime_dtime, mtime_quotient)
    ret = INT(mtime_quotient%quotient)
    CALL deallocateTimedelta(mtime_td)
    CALL deallocateTimedelta(mtime_dtime)
    CALL deallocateDatetime(mtime_start)
    CALL deallocateDatetime(mtime_stop)
  END FUNCTION get_timedelta_divide_by_seconds


  !> Utility function.
  !
  !  @return Timedelta object from mtime libray representing a given
  !  span of milliseconds
  !
  FUNCTION getTimedeltaFromMS(timeMS) RESULT(ret)
    TYPE(timedelta), POINTER :: ret
    INTEGER(i8) :: timeMS !< time span given in milliseconds
    ! local variables
    CHARACTER(LEN=MAX_TIMEDELTA_STR_LEN)  ::  td_string
    CALL getPTStringFromMS(timeMS, td_string)
    ret => newTimeDelta(td_string)
  END FUNCTION getTimedeltaFromMS


  !> compute an ISO 8601 datetime string from a "t_datetime" object
  !
  SUBROUTINE get_datetime_string(datetime_string, timestamp, opt_add_seconds, opt_td_string)
    CHARACTER(LEN=MAX_DATETIME_STR_LEN), INTENT(INOUT) :: datetime_string
    TYPE(t_datetime),                    INTENT(INOUT) :: timestamp
    INTEGER, OPTIONAL,                   INTENT(IN)    :: opt_add_seconds !< additional offset
    CHARACTER(LEN=MAX_TIMEDELTA_STR_LEN), OPTIONAL     :: opt_td_string
    ! local variables
    INTEGER                  :: additional_days, iadd_days
    TYPE(datetime),  POINTER :: mtime_datetime
    TYPE(timedelta), POINTER :: mtime_td, delta_1day

    CHARACTER(LEN=MAX_DATETIME_STR_LEN)  :: result_string
    CHARACTER(LEN=MAX_TIMEDELTA_STR_LEN) :: td_string

    IF (NINT(timestamp%second) >= 60) THEN
      timestamp%second = 0
      timestamp%minute = timestamp%minute + 1
    ENDIF

    mtime_datetime => newDatetime(timestamp%year, timestamp%month, timestamp%day, &
      &                           timestamp%hour, timestamp%minute, NINT(timestamp%second), 0)

    IF (PRESENT(opt_add_seconds)) THEN
      IF (opt_add_seconds>0) THEN

        ! create a "timedelta" object
        delta_1day => newTimedelta("P01D")          ! create a time delta for 1 day
        CALL get_duration_string(opt_add_seconds, td_string, additional_days)
        DO iadd_days=1,additional_days
          mtime_datetime = mtime_datetime + delta_1day
        END DO
        mtime_td => newTimedelta(TRIM(td_string))
        mtime_datetime =   mtime_datetime + mtime_td
        CALL deallocateTimedelta(mtime_td)
        CALL deallocateTimedelta(delta_1day)
      ENDIF
    END IF
    IF (PRESENT(opt_td_string)) THEN
      mtime_td => newTimedelta(TRIM(opt_td_string))
      mtime_datetime =   mtime_datetime + mtime_td
      CALL deallocateTimedelta(mtime_td)
    ENDIF

    IF (PRESENT(opt_add_seconds) .AND. PRESENT(opt_td_string)) THEN
      CALL finish('get_datetime_string','to many optional arguments')
    ENDIF

    CALL datetimeToString(mtime_datetime, result_string)
    CALL deallocateDatetime(mtime_datetime)
    datetime_string = result_string
  END SUBROUTINE get_datetime_string


  SUBROUTINE get_duration_string(iseconds, td_string, additional_days) 
    INTEGER,                              INTENT(IN)    :: iseconds
    CHARACTER(LEN=MAX_TIMEDELTA_STR_LEN), INTENT(INOUT) :: td_string
    INTEGER,                              INTENT(OUT)   :: additional_days
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::get_duration_string"
    INTEGER :: seconds, hours, minutes

    ! additional_days: Additional days to be added to the output
    ! interval. This namelist parameter is required when the output
    ! interval has been provided as a number of seconds, e.g. which is
    ! so large that it cannot be converted into a valid ISO duration
    ! string.
    additional_days = 0

    seconds = iseconds
    ! create a "timedelta" object
    IF (seconds <= 8640) THEN
      ! for small durations: use mtime's conversion routine
      seconds = seconds*1000
      CALL getPTStringFromMS(INT(seconds,i8), td_string)
    ELSE
      ! for larger durations: convert seconds to duration string
       additional_days = seconds/86400
       seconds = seconds - 86400*additional_days
       hours   = seconds/3600
       seconds = seconds - 3600*hours
       minutes = seconds/60
       seconds = seconds - 60*minutes

       WRITE (td_string ,'(a,3(I0.2,a))') &
            & 'PT',                       &
            & hours  ,'H', minutes,'M',seconds,'S'
    END IF
  END SUBROUTINE get_duration_string


  !> compute an ISO 8601 datetime string from another date-time string
  !
  SUBROUTINE get_datetime_string_str(datetime_string, timestamp, opt_add_seconds)
    CHARACTER(LEN=*),                    INTENT(INOUT) :: datetime_string
    CHARACTER(LEN=*),                    INTENT(IN)    :: timestamp
    INTEGER, OPTIONAL,                   INTENT(IN)    :: opt_add_seconds !< additional offset
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::get_datetime_string_str"
    INTEGER                  :: iadd_days, additional_days
    TYPE(datetime),  POINTER :: mtime_datetime
    TYPE(timedelta), POINTER :: mtime_td, delta_1day
    CHARACTER(LEN=MAX_DATETIME_STR_LEN)  :: result_string
    CHARACTER(LEN=MAX_TIMEDELTA_STR_LEN) :: td_string

    mtime_datetime => newDatetime(timestamp)
    IF (PRESENT(opt_add_seconds)) THEN
      IF (opt_add_seconds>0) THEN

        ! create a "timedelta" object
        delta_1day => newTimedelta("P01D")          ! create a time delta for 1 day
        CALL get_duration_string(opt_add_seconds, td_string, additional_days)
        DO iadd_days=1,additional_days
          mtime_datetime = mtime_datetime + delta_1day
        END DO
        mtime_td => newTimedelta(TRIM(td_string))
        mtime_datetime =   mtime_datetime + mtime_td
        CALL deallocateTimedelta(mtime_td)
        CALL deallocateTimedelta(delta_1day)
      ENDIF
    END IF

    CALL datetimeToString(mtime_datetime, result_string)
    CALL deallocateDatetime(mtime_datetime)
    datetime_string = result_string
  END SUBROUTINE get_datetime_string_str

! SUBROUTINE get_duration_string_real(iseconds, td_string) 
!    REAL(wp), INTENT(IN)      :: iseconds
!    CHARACTER(LEN=MAX_TIMEDELTA_STR_LEN), INTENT(INOUT) :: td_string
!    ! local variables
!    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::get_duration_string_real"
!    REAL(wp) :: seconds
!
!    seconds = iseconds
!    ! create a "timedelta" object
!    IF (seconds <= 8640) THEN
!      ! for small durations: use mtime's conversion routine
!      seconds = seconds*1000
!      CALL getPTStringFromMS(NINT(seconds,i8), td_string)
!    ENDIF
!  END SUBROUTINE get_duration_string_real

END MODULE mo_mtime_extensions
