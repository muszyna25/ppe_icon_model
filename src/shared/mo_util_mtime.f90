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
  USE mo_impl_constants,           ONLY: MAX_CHAR_LENGTH
  USE mtime,                       ONLY: datetime, timedelta, newTimeDelta, OPERATOR(-),   &
    &                                    OPERATOR(+), juliandelta, newJulianDelta,         &
    &                                    timeDeltaToJulianDelta, deallocateJulianDelta,    &
    &                                    deallocateTimeDelta
  USE mo_util_string,              ONLY: t_keyword_list,                   &
                                         associate_keyword, with_keywords, &
                                         int2string


  IMPLICIT NONE

  PUBLIC :: t_mtime_utils, mtime_utils
  PUBLIC ::  FMT_DDHHMMSS_ANNOTATED, FMT_DDDHHMMSS_ANNOTATED, &
    &        FMT_DDHHMMSS, FMT_DDDHHMMSS, FMT_DDDHH, FMT_DDHHMMSS_DAYSEP
  PUBLIC :: assumePrevMidnight
  PUBLIC :: assumeNextMidnight

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


END MODULE mo_util_mtime
