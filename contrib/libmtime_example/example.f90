! libmtime example
! 2013-09-17, F. Prill (DWD)
!
MODULE libmtime_extensions

  USE, INTRINSIC :: iso_c_binding, ONLY: c_char, c_ptr, c_f_pointer, c_associated, c_bool, c_loc, c_int64_t, c_null_char
  USE mtime, ONLY: datetime, event, max_timedelta_str_len
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: getTriggerNextEventAtDatetime
  PUBLIC :: isCurrentEventActive
  PUBLIC :: getPTStringFromMS

  INTERFACE
    FUNCTION my_isCurrentEventActive(my_event, current_dt) RESULT(ret) BIND(c,name='isCurrentEventActive')
      IMPORT :: c_bool, c_ptr
      LOGICAL(c_bool) :: ret
      TYPE(c_ptr), VALUE :: my_event
      TYPE(c_ptr), VALUE :: current_dt
    END FUNCTION my_isCurrentEventActive

    FUNCTION my_getptstringfromms(ms, ptstr) RESULT(string_ptr) BIND(c, name='getPTStringFromMS')
      IMPORT :: c_int64_t, c_ptr, c_char
      TYPE(c_ptr) :: string_ptr
      INTEGER(c_int64_t), VALUE :: ms
      CHARACTER(c_char), DIMENSION(*) :: ptstr
    END FUNCTION my_getptstringfromms
    
  END INTERFACE

CONTAINS 
  
  FUNCTION getTriggerNextEventAtDatetime(my_event) RESULT(ret_datetime)
    TYPE(datetime), POINTER :: ret_datetime
    TYPE(event), POINTER :: my_event
    TYPE(c_ptr) :: c_pointer
    c_pointer = my_event%triggerNextEventDateTime
    IF (C_ASSOCIATED(c_pointer)) THEN
      CALL C_F_POINTER(c_pointer, ret_datetime)
    ELSE
      ret_datetime => NULL()
    ENDIF
  END FUNCTION getTriggerNextEventAtDatetime

  FUNCTION isCurrentEventActive(my_event, current_dt) 
    LOGICAL :: isCurrentEventActive
    TYPE(datetime),  POINTER :: current_dt
    TYPE(event),     POINTER :: my_event
    isCurrentEventActive = my_isCurrentEventActive(C_LOC(my_event), C_LOC(current_dt))
  END FUNCTION isCurrentEventActive

  SUBROUTINE getPTStringFromMS(ms, string)
    INTEGER, INTENT(in) :: ms
    CHARACTER(len=max_timedelta_str_len), INTENT(out) :: string
    TYPE(c_ptr) :: dummy_ptr
    INTEGER :: i
    INTEGER(c_int64_t) :: ms_int64

    ms_int64 = ms
    dummy_ptr = my_getptstringfromms(ms_int64, string)
    char_loop: DO i = 1 , LEN(string)
      IF (string(i:i) == c_null_char) EXIT char_loop
    END DO char_loop
    string(i:LEN(string)) = ' '
  END SUBROUTINE getPTStringFromMS

END MODULE libmtime_extensions


PROGRAM main

  USE mtime, ONLY: setCalendar, resetCalendar, PROLEPTIC_GREGORIAN,                   &
    &              max_calendar_str_len, max_datetime_str_len,                        &
    &              max_timedelta_str_len, max_eventname_str_len,                      &
    &              calendarToString, datetime, timedelta, event,                      &
    &              newDatetime, deallocateDatetime, datetimeToString,                 &
    &              newTimedelta, deallocateTimedelta, timedeltaToString,              &
    &              OPERATOR(+),                                                       &
    &              newEvent, deallocateEvent, eventToString
  USE libmtime_extensions, ONLY: getTriggerNextEventAtDatetime,                       &
    &              isCurrentEventActive, getPTStringFromMS

  IMPLICIT NONE

  TYPE(datetime),  POINTER :: date1, date2
  TYPE(timedelta), POINTER :: delta1
  TYPE(event),     POINTER :: event1
  CHARACTER(len=max_calendar_str_len)  :: calendar
  CHARACTER(len=max_datetime_str_len)  :: date1_string, date2_string
  CHARACTER(len=max_timedelta_str_len) :: delta1_string
  CHARACTER(len=max_eventname_str_len) :: event1_string
  LOGICAL :: l_active

  WRITE (0,*) "libmtime example."

  ! -----------------------------------------------
  ! part 1: set a date/time object and increment it
  ! -----------------------------------------------

  CALL setCalendar(PROLEPTIC_GREGORIAN)
  CALL calendarToString(calendar)
  WRITE (0,*) "calendar_type : ", TRIM(calendar)

  date1 => newDatetime("2010-01-01T00:00:00")
  CALL datetimeToString(date1, date1_string)
  WRITE (0,*) "initial date  : ", TRIM(date1_string)

  delta1 => newTimedelta("PT12H00M01S")
  CALL timedeltaToString(delta1, delta1_string)
  WRITE (0,*) "timedelta     : ", TRIM(delta1_string)

  date1 = date1 + delta1
  CALL datetimeToString(date1, date1_string)
  WRITE (0,*) "increm. date  : ", TRIM(date1_string)

  ! -----------------------
  ! part 2: define an event
  ! -----------------------

  CALL getPTStringFromMS(97100, date2_string)
  WRITE (0,*) "getPTStringFromMS = ", TRIM(date2_string)

  event1 => newEvent("regular output",            &  ! name
    &                "2010-01-01T00:00:00",       &  ! ref.date
    &                "2010-01-01T00:00:00",       &  ! first
    &                "2010-01-05T00:00:00",       &  ! last
    &                "P02DT00H00M01S")                  ! interval
  CALL eventToString(event1, event1_string)
  WRITE (0,*) "event : ", TRIM(event1_string)
  date2 =>  getTriggerNextEventAtDatetime(event1)
  CALL datetimeToString(date2, date2_string)
  WRITE (0,*) "next trigger date  : ", TRIM(date2_string)

  ! --------------------
  ! part 3: looping test
  ! --------------------

  EVENT_LOOP : DO
    date1 = getTriggerNextEventAtDatetime(event1)
    l_active = isCurrentEventActive(event1, date1)
    WRITE (0,*) "isCurrentEventActive = ", l_active
    IF (.NOT. l_active) EXIT EVENT_LOOP

    CALL datetimeToString(date1, date1_string)
    WRITE (0,*) "trigger date  : ", TRIM(date1_string)
  END DO EVENT_LOOP

  WRITE (0,*) "looping..."

  WRITE (0,*) "looping done."


  ! clean up
  CALL deallocateEvent(event1)
  CALL deallocateDatetime(date1)
  CALL deallocateTimedelta(delta1)
  CALL resetCalendar()

END PROGRAM main
