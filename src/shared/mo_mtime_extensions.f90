!> Additional C-bindings for the mtime library.
!!
!! @author F. Prill, DWD
!!
!! @par Revision History
!! Initial implementation  by  F. Prill, DWD (2013-09-17)
!!
!! @par Copyright
!! 2002-2010 by DWD and MPI-M
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
!! -----------------------------------------------------------------------------------
MODULE mo_mtime_extensions

  USE, INTRINSIC :: iso_c_binding, ONLY: c_char, c_ptr, c_f_pointer, c_associated, c_bool, &
    &               c_loc, c_int64_t, c_null_char, c_int
  USE mo_kind,     ONLY: wp
  USE mtime,       ONLY: datetime, event, MAX_TIMEDELTA_STR_LEN, MAX_DATETIME_STR_LEN, &
    &                    newDatetime, deallocateDatetime, datetimeToString, timedelta, &
    &                    newTimedelta, deallocateTimedelta, OPERATOR(+),               &
    &                    setCalendar, PROLEPTIC_GREGORIAN
  USE mo_datetime, ONLY: t_datetime
  USE mo_exception,ONLY: message
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: getPTStringFromMS
  PUBLIC :: getTimeDeltaFromDateTime
  PUBLIC :: get_duration_string
  PUBLIC :: get_datetime_string

  !> module name
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_mtime_extensions'


  INTERFACE

    SUBROUTINE my_getptstringfromms(ms, ptstr) BIND(c, name='getPTStringFromMS')
#ifdef __SX__
      USE, INTRINSIC :: iso_c_binding, ONLY: c_int64_t, c_ptr, c_char
#else
      import :: c_int64_t, c_ptr, c_char
#endif
      INTEGER(c_int64_t), VALUE :: ms
      CHARACTER(c_char), DIMENSION(*) :: ptstr
    END SUBROUTINE my_getptstringfromms

    SUBROUTINE my_gettimedeltafromdatetime(dt1, dt2, td_return) BIND(c, name='getTimeDeltaFromDateTime')
#ifdef __SX__
      USE, INTRINSIC :: iso_c_binding, ONLY: c_ptr
#else
      IMPORT :: c_ptr
#endif
      TYPE(c_ptr), value :: dt1, dt2
      TYPE(c_ptr), value :: td_return
    END SUBROUTINE my_gettimedeltafromdatetime
  END INTERFACE

  INTERFACE get_datetime_string
    MODULE PROCEDURE get_datetime_string
    MODULE PROCEDURE get_datetime_string_str
  END INTERFACE

CONTAINS 
  
  SUBROUTINE getPTStringFromMS(ms, string)
    INTEGER, INTENT(in) :: ms
    CHARACTER(len=max_timedelta_str_len), INTENT(out) :: string
    INTEGER :: i
    INTEGER(c_int64_t) :: ms_int64

    ms_int64 = INT(ms,c_int64_t)
    CALL my_getptstringfromms(ms_int64, string)
    char_loop: DO i = 1 , LEN(string)
      IF (string(i:i) == c_null_char) EXIT char_loop
    END DO char_loop
    string(i:LEN(string)) = ' '
  END SUBROUTINE getPTStringFromMS

  SUBROUTINE getTimeDeltaFromDateTime(dt1, dt2, td_return)
    TYPE(datetime),  INTENT(IN),    TARGET   :: dt1,dt2 
    TYPE(timedelta), INTENT(INOUT), TARGET   :: td_return     !< OUT

    CALL my_gettimedeltafromdatetime(C_LOC(dt1), C_LOC(dt2), C_LOC(td_return))
  END SUBROUTINE getTimeDeltaFromDateTime


  !> compute an ISO 8601 datetime string from a "t_datetime" object
  !
  SUBROUTINE get_datetime_string(datetime_string, timestamp, opt_add_seconds)
    CHARACTER(LEN=MAX_DATETIME_STR_LEN), INTENT(INOUT) :: datetime_string
    TYPE(t_datetime),                    INTENT(INOUT) :: timestamp
    INTEGER, OPTIONAL,                   INTENT(IN)    :: opt_add_seconds !< additional offset
    ! local variables
    INTEGER                  :: add_seconds, additional_days, iadd_days
    TYPE(datetime),  POINTER :: mtime_datetime
    TYPE(timedelta), POINTER :: mtime_td, delta_1day

    CHARACTER(LEN=MAX_DATETIME_STR_LEN)  :: result_string
    CHARACTER(LEN=MAX_TIMEDELTA_STR_LEN) :: td_string

    CALL setCalendar(PROLEPTIC_GREGORIAN)
    mtime_datetime => newDatetime(timestamp%year, timestamp%month, timestamp%day, &
      &                           timestamp%hour, timestamp%minute, NINT(timestamp%second), ims=0)

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
      CALL getPTStringFromMS(seconds, td_string)
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

    CALL setCalendar(PROLEPTIC_GREGORIAN)
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

END MODULE mo_mtime_extensions
