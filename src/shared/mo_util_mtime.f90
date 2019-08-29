!> This module contains auxiliary routines for mtime date/time objects. 
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
MODULE mo_util_mtime

  USE, INTRINSIC :: iso_c_binding, ONLY: c_int64_t
  USE mo_kind,                     ONLY: wp, i8
  USE mo_exception,                ONLY: message, message_text, finish
  USE mo_impl_constants,           ONLY: MAX_CHAR_LENGTH
  USE mtime,                       ONLY: datetime, newDatetime, timedelta, newTimeDelta,   &
    &                                    OPERATOR(-), OPERATOR(+),                         &
    &                                    juliandelta, newJulianDelta,                      &
    &                                    timeDeltaToJulianDelta, deallocateJulianDelta,    &
    &                                    deallocateTimeDelta, deallocateDatetime,          &
    &                                    getTimeDeltaFromDateTime,                         &
    &                                    getTotalMillisecondsTimedelta
  USE mo_time_config,              ONLY: time_config
  USE mo_util_string,              ONLY: t_keyword_list,                   &
                                         associate_keyword, with_keywords, &
                                         int2string

!DR Test
  USE mo_exception,                ONLY: message, message_text, finish
!DR End Test

  IMPLICIT NONE

  PUBLIC :: t_mtime_utils, mtime_utils
  PUBLIC ::  FMT_DDHHMMSS_ANNOTATED, FMT_DDDHHMMSS_ANNOTATED, &
    &        FMT_DDHHMMSS, FMT_DDDHHMMSS, FMT_DDDHH, FMT_HHH, FMT_DDHHMMSS_DAYSEP
  PUBLIC :: assumePrevMidnight
  PUBLIC :: assumeNextMidnight
  PUBLIC :: getElapsedSimTimeInSeconds
  PUBLIC :: dummyDateTime
  PUBLIC :: t_datetime_ptr
  PUBLIC :: mtime_convert_netcdf_units
  PUBLIC :: mtime_divide_timedelta

  TYPE t_datetime_ptr
    TYPE(datetime), POINTER :: ptr => NULL()
  END TYPE t_datetime_ptr

  PRIVATE

  !> module name
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_util_mtime'

  !> a collection of relative time formats
  CHARACTER(LEN=*), PARAMETER :: FMT_DDHHMMSS_ANNOTATED  = "<dd>D<hh>H<mm>M<ss><ms>S"
  CHARACTER(LEN=*), PARAMETER :: FMT_DDHHMMSS_DAYSEP     = "<dd>D <hh>H<mm>M<ss><ms>S"
  CHARACTER(LEN=*), PARAMETER :: FMT_DDDHHMMSS_ANNOTATED = "<ddd>D<hh>H<mm>M<ss><ms>S"
  CHARACTER(LEN=*), PARAMETER :: FMT_DDHHMMSS            = "<dd><hh><mm><ss><ms>"
  CHARACTER(LEN=*), PARAMETER :: FMT_DDDHHMMSS           = "<ddd><hh><mm><ss><ms>"
  CHARACTER(LEN=*), PARAMETER :: FMT_DDDHH               = "<ddd><hh>"
  CHARACTER(LEN=*), PARAMETER :: FMT_HHH                 = "<hhh>"

  TYPE t_mtime_utils
  CONTAINS
    PROCEDURE, NOPASS :: ddhhmmss => t_mtime_utils_ddhhmmss
  END TYPE t_mtime_utils

  ! global object serving as a collection of static utility routines
  TYPE(t_mtime_utils) :: mtime_utils
  
CONTAINS 

  !> @return DDHHMMSS string of a given time span.
  !
  !  Note: We return the forecast time delta as an ISO 8601-like
  !        string, where the exact format is given by the "format" argument.
  !
  FUNCTION t_mtime_utils_ddhhmmss(start_date, end_date, fmt)  RESULT(result_str)
    CHARACTER(len=MAX_CHAR_LENGTH)  :: result_str
    TYPE(datetime),    INTENT(IN)   :: start_date, end_date
    CHARACTER(LEN=*),  INTENT(IN)   :: fmt
    
    TYPE(timedelta),       POINTER  :: time_delta
    TYPE(juliandelta),     POINTER  :: julian_delta
    TYPE (t_keyword_list), POINTER  :: keywords => NULL()

    ! compute time delta:
    time_delta => newTimeDelta("P01D")
    time_delta = end_date - start_date
    ! compute the time delta in a Julian calendar, which allows more
    ! than 30 days:
    julian_delta => newJulianDelta('+', 0_c_int64_t, 0_c_int64_t)
    CALL timeDeltaToJulianDelta(time_delta,start_date,julian_delta)

    ! now collect the feasible "building-blocks"
    IF (INT(julian_delta%day) < 100) THEN
      CALL associate_keyword("<dd>",  TRIM(int2string(INT(julian_delta%day), '(i2.2)')), keywords)
    ELSE
      CALL associate_keyword("<dd>",  TRIM(int2string(INT(julian_delta%day), '(i0)')),   keywords)
    END IF
    CALL associate_keyword("<ddd>", TRIM(int2string(INT(julian_delta%day), '(i3.3)')),   keywords)
    CALL associate_keyword("<hhh>", TRIM(int2string(INT(julian_delta%day)*24+time_delta%hour, '(i3.3)')),   keywords)
    CALL associate_keyword("<hh>",  TRIM(int2string(time_delta%hour,       '(i2.2)')),   keywords)
    CALL associate_keyword("<mm>",  TRIM(int2string(time_delta%minute,     '(i2.2)')),   keywords)
    CALL associate_keyword("<ss>",  TRIM(int2string(time_delta%second,     '(i2.2)')),   keywords)
    IF (time_delta%ms /= 0) THEN
      CALL associate_keyword("<ms>",  '.'//TRIM(int2string(time_delta%second,'(i3.3)')),   keywords)
    ELSE
      CALL associate_keyword("<ms>",  "", keywords)
    END IF
    ! replace keywords in "ddhhmmss"
    result_str = TRIM(with_keywords(keywords, fmt))

    CALL deallocateTimeDelta(time_delta)
    CALL deallocateJulianDelta(julian_delta)
  END FUNCTION t_mtime_utils_ddhhmmss


  !>
  !! Returns previous 'midnight' datetime for given datetime
  !!
  !! Returns previous 'midnight' datetime for given datetime in datetime-format
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2017-12-01)
  !!
  TYPE(datetime) FUNCTION assumePrevMidnight(current_date)
    TYPE(datetime), INTENT(IN) :: current_date

    assumePrevMidnight = current_date
    !
    assumePrevMidnight%time%hour   = 0 
    assumePrevMidnight%time%minute = 0
    assumePrevMidnight%time%second = 0
    assumePrevMidnight%time%ms     = 0

  END FUNCTION assumePrevMidnight

  !>
  !! Returns nest 'midnight' datetime for given datetime
  !!
  !! Returns next 'midnight' datetime for given datetime in datetime-format
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2017-12-01)
  !!
  TYPE(datetime) FUNCTION assumeNextMidnight(current_date)
    TYPE(datetime), INTENT(IN) :: current_date

    ! local
    TYPE(timedelta), POINTER :: td_1day  => NULL()
    !---------------------------------------------------------

    td_1day =>newTimedelta("P01D")
    assumeNextMidnight = current_date + td_1day
    !
    assumeNextMidnight%time%hour   = 0 
    assumeNextMidnight%time%minute = 0
    assumeNextMidnight%time%second = 0
    assumeNextMidnight%time%ms     = 0

    CALL deallocateTimedelta(td_1day)

  END FUNCTION assumeNextMidnight


  !>
  !! Returns elapsed simulation time in seconds
  !!
  !! Elapsed simulation time is computed as the timedelta 
  !! between the current datetime and the anchor  
  !! datetime. Unless specified otherwise, the anchor date  
  !! is set to the experiment start date time_config%tc_exp_startdate.
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2017-12-01)
  !!
  REAL(wp) FUNCTION getElapsedSimTimeInSeconds(datetime_current, anchor_datetime) RESULT(sim_time)
    !
    TYPE(datetime)          , INTENT(IN) :: datetime_current
    TYPE(datetime), OPTIONAL, INTENT(IN) :: anchor_datetime   ! anchor date
    !
    ! local
    TYPE(timedelta), POINTER :: time_diff => NULL()
    TYPE(datetime)           :: anchor                        ! anchor date for 
                                                              ! timeDelta computation
    !---------------------------------------------------------

    If (PRESENT(anchor_datetime)) THEN
      anchor = anchor_datetime
    ELSE
      ! set anchor date to experiment start date
      anchor = time_config%tc_exp_startdate
    ENDIF

    time_diff  => newTimedelta("PT0S")
    time_diff  =  getTimeDeltaFromDateTime(datetime_current, anchor)
    sim_time   =  getTotalMillisecondsTimedelta(time_diff, datetime_current)*1.e-3_wp

    CALL deallocateTimedelta(time_diff)

  END FUNCTION getElapsedSimTimeInSeconds



  !>
  !! Creates a datetime object with a defined dummy date
  !!
  !! Creates a datetime object with a defined dummy date
  !! Can be used for initializing a newly defined datetime 
  !! object and for checking, whether the datetime object 
  !! at hand has already been touched and filed with a meaningful value.
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2018-12-01)
  !!
  TYPE(datetime) FUNCTION dummyDateTime()
    !
    ! local
    TYPE(datetime), POINTER :: dummy_date
    !---------------------------------------------------------

    dummy_date    => newDatetime("0-01-01T00:00:00.000")
    dummyDateTime =  dummy_date

    CALL deallocateDatetime(dummy_date)

  END FUNCTION dummyDateTime

  ! via seconds... should be solved by mtime, could be done without
  ! conversion to a prespecified time unit.
  SUBROUTINE mtime_divide_timedelta (dividend, divisor, quotient)
    TYPE(timedelta), POINTER, INTENT(in   ) :: dividend
    TYPE(timedelta), POINTER, INTENT(in   ) :: divisor
    REAL(wp),                 INTENT(  out) :: quotient
    INTEGER(i8)                             :: dd, ds

    CALL mtime_timedelta_to_seconds(dividend, dd)
    CALL mtime_timedelta_to_seconds(divisor,  ds)

    quotient = REAL(dd,wp) / REAL(ds,wp)

  END SUBROUTINE mtime_divide_timedelta

  SUBROUTINE mtime_timedelta_to_seconds (mtime_dt, seconds_dt)
    TYPE(timedelta), POINTER, INTENT(in   ) :: mtime_dt
    INTEGER(i8),              INTENT(  out) :: seconds_dt
    CHARACTER(*),                 PARAMETER :: routine = &
      & modname//"::mtime_timedelta_to_seconds"

    IF (mtime_dt%year > 0 .OR. mtime_dt%month > 0) THEN
      CALL finish(routine, "year or month nonzero, result is undefined")
    ENDIF

    seconds_dt =   86400 * INT(mtime_dt%day)    &
      &          + 3600  * INT(mtime_dt%hour)   &
      &          + 60    * INT(mtime_dt%minute) &
      &          +         INT(mtime_dt%second)

  END SUBROUTINE mtime_timedelta_to_seconds

  ! the expected input looks similar to this "seconds since 2013-04-24T00:00:00"

  SUBROUTINE mtime_convert_netcdf_units (unit_string, starttime, timeincrement)
    CHARACTER(len=*),         INTENT(in   ) :: unit_string
    TYPE(datetime),  POINTER, INTENT(  out) :: starttime
    TYPE(timeDelta), POINTER, INTENT(  out) :: timeincrement

    INTEGER           :: errno
    INTEGER           :: idx0, idx1, idx2, idx3, idx4, idx5, idx6
    INTEGER           :: month, day
    INTEGER(i8)       :: year
    CHARACTER(len=32) :: iso8601_string
    CHARACTER(len=32) :: timedelta_string

    idx1 = SCAN(unit_string,          " ")
    idx2 = SCAN(unit_string(idx1+1:), " ")
    idx2 = idx1+idx2
    idx3 = SCAN(unit_string(idx2+1:), " ")
    idx3 = idx2+idx3

    idx0 = 1
    idx4 = LEN_TRIM(unit_string)
    idx5 = SCAN(unit_string(idx2+1:idx3-1), "-")
    idx6 = SCAN(unit_string(idx2+idx5+1:idx3-1), "-")
    idx6 = idx5+idx6

    READ(unit_string(idx2+1:idx2+idx5-1),*) year
    READ(unit_string(idx2+idx5+1:idx2+idx6-1),*) month
    READ(unit_string(idx2+idx6+1:idx3-1),*) day

    WRITE(iso8601_string,'(i0,a,i2.2,a,i2.2,a,a)') &
         year, '-', month, '-', day, 'T', unit_string(idx3+1:idx4)

    IF (unit_string(idx1+1:idx2) /= "since") THEN
      CALL finish(modname, "Unit string is in unknown format"//unit_string(idx1+1:idx2))
    ENDIF

    SELECT CASE ( unit_string(:idx1-1) )
    CASE ("seconds")
      timedelta_string = "PT01S"
    CASE ("hours")
      timedelta_string = "PT01H"
    CASE ("days")
      timedelta_string = "PT01D"
    CASE default
      CALL finish(modname, "Unknown time increment: "//unit_string(:idx1+1))
    END SELECT

    starttime     => newdatetime(TRIM(iso8601_string), errno)
    IF (errno /= 0)  CALL finish(modname, "Error in conversion of date: "//TRIM(iso8601_string))

    timeincrement => newTimeDelta(timedelta_string, errno)
    IF (errno /= 0)  CALL finish(modname, "Error in conversion of unit: "//timedelta_string)

  END SUBROUTINE mtime_convert_netcdf_units

END MODULE mo_util_mtime
