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
   PUBLIC :: datetime_str_equal
   PUBLIC :: timedelta_str_equal

  !> module name
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_mtime_extensions'

CONTAINS 

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

END MODULE mo_mtime_extensions
