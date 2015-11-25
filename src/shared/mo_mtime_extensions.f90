!> Additional functions for the mtime library.
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

  USE mo_kind,     ONLY: wp, i8
  USE mtime,       ONLY: datetime, MAX_TIMEDELTA_STR_LEN, MAX_DATETIME_STR_LEN,        &
    &                    newDatetime, deallocateDatetime, timedelta,                   &
    &                    newTimedelta, deallocateTimedelta, OPERATOR(+), OPERATOR(-),  &
    &                    getPTStringFromMS,                                            &
    &                    divisionquotienttimespan, OPERATOR(==),                       &
    &                    divideTimeDeltaInSeconds, OPERATOR(<)
  USE mo_util_string, ONLY: int2string, real2string
  IMPLICIT NONE
  
  PRIVATE
   PUBLIC :: convert_datetime_string_old2new
   PUBLIC :: datetime_str_equal
   PUBLIC :: timedelta_str_equal
   PUBLIC :: get_timedelta_divide_by_seconds
   PUBLIC :: getTimedeltaFromMS

  !> module name
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_mtime_extensions'

CONTAINS 

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

END MODULE mo_mtime_extensions
