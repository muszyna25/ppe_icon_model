!> @file libmtime.f90
!!
!! @brief Providing the Fortran language bindings for libmtime
!!
!! @author  Luis Kornblueh, Max Planck Institute for Meteorology
!! @author  Rahul Sinha, Max Planck Institute for Meteorology
!!
!! @defgroup FortranBindings libmtime Fortran languange bindings
!! @{
!!
!___________________________________________________________________________________________________________
!>
!! @brief Provides the calendar to the users, abstracting the different calendars available.
!!
!! @details
!!
!! Three calendar types are provided:
!!
!!   - a proleptic Gregorian calendar
!!   - a calendar with 365 days per year without leap years
!!   - a calendar with 360 days per year and each month having 30 days
!!
!! To use a specific calendar a call to setCalendar with the respective selector must be
!! done. The implementation is based on a singleton concept meaning that only one calendar
!! can be active at a time. To release a calendar a call to resetCalendar has to be done.
!!
!___________________________________________________________________________________________________________
module mtime_calendar
  !
  use, intrinsic :: iso_c_binding, only: c_int, c_char, c_null_char
  !
  implicit none
  !
  private
  !   
#ifdef __SX__
  integer, parameter :: calendar_not_set    = 0 
  integer, parameter :: proleptic_gregorian = 1
  integer, parameter :: year_of_365_days    = 2
  integer, parameter :: year_of_360_days    = 3
#else
  enum, bind(c)
    enumerator :: calendar_not_set    = 0   ! calendar is not defined yet
    enumerator :: proleptic_gregorian = 1   ! proleptic Gregorian calendar
    enumerator :: year_of_365_days    = 2   ! 365 day year without leap years
    enumerator :: year_of_360_days    = 3   ! 360 day year with 30 day months
  end enum
#endif
  !
  public :: max_calendar_str_len
  public :: calendar_not_set, proleptic_gregorian, year_of_365_days, year_of_360_days
  public :: setCalendar
  public :: resetCalendar
  public :: calendarType
  public :: calendarToString
  !
  !> provides a string length for toString 
  integer, parameter :: max_calendar_str_len = 32
  !
  !> @cond DOXYGEN_IGNORE_THIS
  interface
    !
    subroutine setCalendar(ct) bind(c, name='initCalendar')
#ifdef __SX__
      use, intrinsic :: iso_c_binding, only: c_int  
#else
      import :: c_int
#endif
      integer(c_int), value :: ct
    end subroutine setCalendar
    !
    subroutine resetCalendar() bind(c, name='freeCalendar')
    end subroutine resetCalendar
    !
    function calendarType() bind(c, name='getCalendarType')
#ifdef __SX__
      use, intrinsic :: iso_c_binding, only: c_int  
#else
      import :: c_int
#endif
      integer(c_int) :: calendarType
    end function calendarType
    !
    subroutine my_calendartostring(calendar) bind(c, name='calendarToString') 
#ifdef __SX__
      use, intrinsic :: iso_c_binding, only: c_char  
#else
      import :: c_char
#endif
      character(c_char), dimension(*) :: calendar
    end subroutine my_calendartostring
    !
  end interface
  !> @endcond DOXYGEN_IGNORE_THIS
  !
contains
  !
#ifdef DOXYGEN_DOCUMENTATION_ONLY
  !> 
  !! @brief Initialize a new calendar.
  !!
  !! setCalendar is done at the very begining to select one of the
  !! provided calendar libraries. It intializes the calendar to one of:
  !!
  !! - proleptic_gregorian
  !! - year_of_365_days
  !! - year_of_360_days
  !! 
  !! The calendar type and hence it's behaviour (Calendar to Julian
  !! conversion and vice versa) is fixed for the lifetime of the
  !! selected calendar.  Attempts to change the calendar type on the
  !! fly is discouraged. The lib has built-in checks to reject
  !! change attempts at run time.  However, a calendar can be
  !! "re-initialized" after calling resetCalendar(), but this is not
  !! advised. 
  !!
  !! MANTRA: Know what you are doing before you do it and do it
  !! right the first time.
  !!
  !! @param[in]  ct      the calendar type
  subroutine setCalendar(ct) bind(c, name='initCalendar')
    integer(c_int), value :: ct
  end subroutine setCalendar
  !>
  !! @brief called to discard the selected calendar type 
  subroutine resetCalendar() bind(c, name='freeCalendar')
  end subroutine resetCalendar
  !>
  !! @brief to get an idea, which calendar type is selected
  !!
  !! @return an integer denoting the calendar in use
  function calendarType() bind(c, name='getCalendarType')
    integer(c_int) :: calendarType
  end function calendarType
  !
#endif
  !>
  !! @brief convert the calendar identifier into a human readable string
  !!
  !! @param[out]  string      the calendar type verbose
  !!
  subroutine calendarToString(string)
    character(len=max_calendar_str_len), intent(out) :: string
    integer :: i
    call my_calendartostring(string)
    char_loop: do i = 1 , len(string)
      if (string(i:i) == c_null_char) exit char_loop
    end do char_loop
    string(i:len(string)) = ' '
  end subroutine calendarToString
  !
end module mtime_calendar
!>
!! @brief Julian Day Calendar and some operations supported on julian dates.
!!
!! @details
!___________________________________________________________________________________________________________
module mtime_julianday
  !
  use, intrinsic :: iso_c_binding, only: c_int64_t, c_char, c_null_char, c_ptr, c_loc, c_f_pointer
  !
  implicit none
  !
  private
  !
  public :: max_julianday_str_len
  public :: julianday
  public :: newJulianday
  public :: deallocateJulianday
  public :: juliandayToString
  !
  !> provides a string length for toString 
  integer, parameter :: max_julianday_str_len = 32   
  !
  type, bind(c) :: julianday
    integer(c_int64_t) :: day  !> the actual Julian day
    integer(c_int64_t) :: ms   !> the milisecond on that particular day
  end type julianday
  !
  interface
     function my_newjulianday(day, ms) result(c_pointer) bind(c, name='newJulianDay')
#ifdef __SX__
      use, intrinsic :: iso_c_binding, only: c_int64_t, c_ptr
#else
      import :: c_int64_t, c_ptr
#endif
      type(c_ptr) :: c_pointer
      integer(c_int64_t), intent(in) :: day
      integer(c_int64_t), intent(in) :: ms
    end function my_newjulianday
    !
    subroutine my_deallocatejulianday(jd) bind(c,name='deallocateJulianDay')
#ifdef __SX__
      use, intrinsic :: iso_c_binding, only: c_ptr
#else
      import :: c_ptr
#endif
      type(c_ptr), value :: jd
    end subroutine my_deallocatejulianday
    !
    function my_juliandaytostring(my_julianday, string) result(string_ptr) bind(c, name='juliandayToString')
#ifdef __SX__
      use, intrinsic :: iso_c_binding, only: c_ptr, c_char
#else
      import :: c_ptr, c_char
#endif
      type(c_ptr) :: string_ptr
      type(c_ptr), value :: my_julianday
      character(c_char), dimension(*) :: string
    end function my_juliandaytostring
    !
  end interface
  !
contains
  !
  !>
  !! @brief construct a new Julian date 
  !!
  !! @param[in] day            the Julian day
  !! @param[in] ms             an integer denoting the actual milli seconds of a day 
  !! @return    ret_julianday  a pointer of type(julianday)
  function newJulianday(day, ms) result(ret_julianday)
    type(julianday), pointer :: ret_julianday
    integer(c_int64_t), intent(in) :: day
    integer(c_int64_t), intent(in) :: ms
    type(c_ptr) :: c_pointer
    c_pointer = my_newjulianday(day, ms)
    call c_f_pointer(c_pointer, ret_julianday)
  end function newJulianday
  !
  !>
  !! @brief destructor for a Julian date 
  !!
  !! @param     my_julianday   a pointer of type(julianday)
  subroutine deallocateJulianday(my_julianday)
    type(julianday), pointer :: my_julianday
    call my_deallocatejulianday(c_loc(my_julianday))
    my_julianday => null()
  end subroutine deallocateJulianday
  !
  !>
  !! @brief get Julian day as a string.
  !!
  !! @param[in]  my_julianday   a pointer to type(julianday). The Julian day to be converted to a string 
  !! @param[out] string         the Julian day verbose
  subroutine juliandayToString(my_julianday, string)
    type(julianday), pointer :: my_julianday
    character(len=max_julianday_str_len), intent(out) :: string
    type(c_ptr) :: dummy_ptr
    integer :: i
    dummy_ptr = my_juliandaytostring(c_loc(my_julianday), string)
    char_loop: do i = 1 , len(string)
      if (string(i:i) == c_null_char) exit char_loop
    end do char_loop
    string(i:len(string)) = ' '
  end subroutine juliandayToString
  !
end module mtime_julianday
!>
!! @brief Date and some operations supported on Date.
!!
!! @details
!!
!___________________________________________________________________________________________________________
module mtime_date
  !
  use, intrinsic :: iso_c_binding, only: c_int, c_int64_t, c_char, c_ptr, c_null_char, c_loc, c_f_pointer
  !
  implicit none
  !
  private
  !
  public :: max_date_str_len
  public :: date
  public :: newDate
  public :: deallocateDate
  public :: dateToString
  !
  !> provides a string length for toString 
  integer, parameter :: max_date_str_len = 32
  !
  type, bind(c) :: date
    integer(c_int64_t) :: year
    integer(c_int)     :: month
    integer(c_int)     :: day
  end type date
  !
  interface
    !
    function my_newdate(string) result(c_pointer) bind(c, name='newDate')
#ifdef __SX__
      use, intrinsic :: iso_c_binding, only: c_char, c_ptr
#else
      import :: c_char, c_ptr
#endif
      type(c_ptr) :: c_pointer
      character(c_char), dimension(*) :: string
    end function my_newdate
    !
    ! TODO: needs additional interface
    subroutine my_deallocatedate(d) bind(c, name='deallocateDate') 
#ifdef __SX__
      use, intrinsic :: iso_c_binding, only: c_ptr
#else
      import :: c_ptr
#endif
      type(c_ptr), value :: d
    end subroutine my_deallocatedate
    !
    function my_datetostring(my_date, string) result(string_ptr) bind(c, name='dateToString') 
#ifdef __SX__
      use, intrinsic :: iso_c_binding, only: c_ptr, c_char
#else
      import :: c_ptr, c_char
#endif
      type(c_ptr) :: string_ptr
      type(c_ptr), value :: my_date
      character(c_char), dimension(*) :: string
    end function my_datetostring
    !
  end interface
  !
contains
  !
  !>
  !! @brief construct a new date 
  !!
  !! @param[in] string         an ISO 8601 conforming date string
  !! @return    ret_date       a pointer of type(date)
  function newDate(string) result(ret_date)
    type(date), pointer :: ret_date
    character(len=*), intent(in) :: string
    type(c_ptr) :: c_pointer
    c_pointer = my_newdate(trim(string)//c_null_char)
    call c_f_pointer(c_pointer, ret_date)
  end function newDate
  !
  !! @brief destructor for a date 
  !!
  !! @param[in] my_date        a pointer of type(date)
  subroutine deallocateDate(my_date)
    type(date), pointer :: my_date
    call my_deallocatedate(c_loc(my_date))
    my_date => null()
  end subroutine deallocateDate
  !
  subroutine dateToString(my_date, string)
    type(date), pointer :: my_date
    character(len=max_date_str_len) :: string
    type(c_ptr) :: dummy_ptr
    integer :: i
    dummy_ptr = my_datetostring(c_loc(my_date), string)
    char_loop: do i = 1 , len(string)
      if (string(i:i) == c_null_char) exit char_loop
    end do char_loop
    string(i:len(string)) = ' '
  end subroutine dateToString
  !
end module mtime_date
!>
!! @brief Time and some operations supported on Time
!!
!! @details
!!
!! @author  Luis Kornblueh, Max Planck Institute for Meteorology
!! @author  Rahul Sinha, Max Planck Institute for Meteorology
!!
!___________________________________________________________________________________________________________
module mtime_time
  !
  use, intrinsic :: iso_c_binding, only: c_int, c_char, c_ptr, c_null_char, c_loc, c_f_pointer
  !
  implicit none
  !
  private
  !
  public :: max_time_str_len
  public :: time
  public :: newTime
  public :: deallocateTime
  public :: timeToString
  !
  !> provides a string length for toString 
  integer, parameter :: max_time_str_len = 32
  !
  type, bind(c) :: time
    integer(c_int) :: hour
    integer(c_int) :: minute
    integer(c_int) :: second
    integer(c_int) :: ms
  end type time
  !
  interface
    !
    function my_newtime(string) result(c_pointer) bind(c, name='newTime')
#ifdef __SX__
      use, intrinsic :: iso_c_binding, only: c_char, c_ptr
#else
      import :: c_char, c_ptr
#endif
      type(c_ptr) :: c_pointer
      character(c_char), dimension(*) :: string
    end function my_newtime
    !
    subroutine my_deallocatetime(t) bind(c,name='deallocateTime') 
#ifdef __SX__
      use, intrinsic :: iso_c_binding, only: c_ptr
#else
      import :: c_ptr
#endif
      type(c_ptr), value :: t
    end subroutine my_deallocatetime
    !
    function my_timetostring(my_time, string) result(string_ptr) bind(c, name='timeToString') 
#ifdef __SX__
      use, intrinsic :: iso_c_binding, only: c_ptr, c_char
#else
      import :: c_ptr, c_char
#endif
      type(c_ptr) :: string_ptr
      type(c_ptr), value :: my_time
      character(c_char), dimension(*) :: string
    end function my_timetostring
    !
  end interface
  !
contains
  !
  function newTime(string) result(ret_time)
    type(time), pointer :: ret_time
    character(len=*), intent(in) :: string
    type(c_ptr) :: c_pointer
    c_pointer = my_newtime(trim(string)//c_null_char)
    call c_f_pointer(c_pointer, ret_time)
  end function newTime
  !
  subroutine deallocateTime(my_time)
    type(time), pointer :: my_time
    call my_deallocatetime(c_loc(my_time))
    my_time => null()
  end subroutine deallocateTime
  !
  subroutine timeToString(my_time, string)
    type(time), pointer :: my_time
    character(len=max_time_str_len) :: string
    type(c_ptr) :: dummy_ptr
    integer :: i
    dummy_ptr = my_timetostring(c_loc(my_time), string)
    char_loop: do i = 1 , len(string)
      if (string(i:i) == c_null_char) exit char_loop
    end do char_loop
    string(i:len(string)) = ' '
  end subroutine timeToString
  !
end module mtime_time
!>
!! @brief DateTime and some operations supported on DateTime.
!!
!! @details
!!
!! @author  Luis Kornblueh, Max Planck Institute for Meteorology
!! @author  Rahul Sinha, Max Planck Institute for Meteorology
!!
!___________________________________________________________________________________________________________
module mtime_datetime
  !
  use, intrinsic :: iso_c_binding, only: c_int, c_int64_t, c_char, c_ptr, c_null_char, c_loc, c_f_pointer
  !
  use mtime_julianday
  use mtime_date
  use mtime_time
  !
  implicit none
  !
  private
  !
  public :: max_datetime_str_len
  public :: datetime
  public :: newDatetime
  public :: deallocateDatetime
  public :: datetimeToString
  public :: getJulianDayFromDatetime
  !
  public :: assignment(=)
  public :: operator(>)
  public :: operator(<)
  public :: operator(<=)
  public :: operator(>=)
  public :: operator(==)
  public :: operator(/=)
  !
  !> provides a string length for toString 
  integer, parameter :: max_datetime_str_len = 32
  !
  type, bind(c) :: datetime
    type(date) :: date
    type(time) :: time
  end type datetime
  !
  interface newDatetime
    module procedure newdatetimefromstring
    module procedure newdatetimefromraw
    module procedure newdatetimefromconstructandcopy
  end interface newDatetime
  !
  interface assignment (=)
    module procedure replacedatetime
  end interface assignment (=)
  !
  interface operator (>)
    module procedure datetime_gt
  end interface operator (>)
  !
  interface operator (<)
    module procedure datetime_lt
  end interface operator (<)
  !
  interface operator (<=)
    module procedure datetime_lt_or_eq
  end interface operator (<=)
  !
  interface operator (>=)
    module procedure datetime_gt_or_eq
  end interface operator (>=)
  !
  interface operator (==)
    module procedure datetime_eq
  end interface operator (==)
  !
  interface operator (/=)
    module procedure datetime_ne
  end interface operator (/=)
  !
  interface
    !
    function my_newdatetime(string) result(c_pointer) bind(c, name='newDateTime')
#ifdef __SX__
      use, intrinsic :: iso_c_binding, only: c_char, c_ptr
#else
      import :: c_char, c_ptr
#endif
      type(c_ptr) :: c_pointer
      character(c_char), dimension(*) :: string
    end function my_newdatetime
    !
    function my_newrawdatetime(year, month, day, hour, minute, second, ms) result(c_pointer) bind(c, name='newRawDateTime')
#ifdef __SX__
      use, intrinsic :: iso_c_binding, only: c_int64_t, c_int, c_ptr
#else
      import :: c_int64_t, c_int, c_ptr
#endif
      type(c_ptr) :: c_pointer
      integer(c_int64_t) :: year
      integer(c_int) :: month, day, hour, minute, second, ms
    end function my_newrawdatetime
    !
    function my_constructandcopydatetime(dt) result(c_pointer) bind(c,name='constructAndCopyDateTime')
#ifdef __SX__
      use, intrinsic :: iso_c_binding, only: c_ptr
#else
      import :: c_ptr
#endif
      type(c_ptr) :: c_pointer
      type(c_ptr), value :: dt
    end function my_constructandcopydatetime
    !
    subroutine my_deallocatedatetime(dt) bind(c,name='deallocateDateTime')
#ifdef __SX__
      use, intrinsic :: iso_c_binding, only: c_ptr  
#else
      import :: c_ptr
#endif
      type(c_ptr), value :: dt
    end subroutine my_deallocatedatetime
    !
    function my_datetimetostring(my_time, string) result(string_ptr) bind(c, name='datetimeToString')
#ifdef __SX__
      use, intrinsic :: iso_c_binding, only: c_ptr, c_char
#else
      import :: c_ptr, c_char
#endif
      type(c_ptr) :: string_ptr
      type(c_ptr), value :: my_time
      character(c_char), dimension(*) :: string
    end function my_datetimetostring
    !
    function my_replacedatetime(src, dest) result(ret_dest) bind(c, name='replaceDatetime')
#ifdef __SX__
      use, intrinsic :: iso_c_binding, only: c_ptr
#else
      import :: c_ptr
#endif
      type(c_ptr) :: ret_dest
      type(c_ptr), value :: src
      type(c_ptr), value :: dest
    end function my_replacedatetime
    !
    function my_comparedatetime(op1, op2) result(ret) bind(c, name='compareDatetime')
#ifdef __SX__
      use, intrinsic :: iso_c_binding, only: c_ptr, c_int  
#else
      import :: c_ptr, c_int
#endif
      integer(c_int) :: ret
      type(c_ptr), value :: op1
      type(c_ptr), value :: op2
    end function my_comparedatetime
    !
    function my_getnoofdaysinmonthdatetime(dt) bind(c, name='getNoOfDaysInMonthDateTime')
#ifdef __SX__
      use, intrinsic :: iso_c_binding, only: c_int, c_ptr  
#else
      import :: c_int, c_ptr
#endif
      integer(c_int) :: my_getnoofdaysinmonthdatetime
      type(c_ptr), value :: dt
    end function my_getnoofdaysinmonthdatetime
    !
    function my_getnoofdaysinyeardatetime(dt) bind(c, name='getNoOfDaysInYearDateTime')
#ifdef __SX__
      use, intrinsic :: iso_c_binding, only: c_int, c_ptr  
#else
      import :: c_int, c_ptr
#endif
      integer(c_int) :: my_getnoofdaysinyeardatetime
      type(c_ptr), value :: dt
    end function my_getnoofdaysinyeardatetime
    !
    function my_getdayofyearfromdatetime(dt) bind(c, name='getDayOfYearFromDateTime')
#ifdef __SX__
      use, intrinsic :: iso_c_binding, only: c_int, c_ptr  
#else
      import :: c_int, c_ptr
#endif
      integer(c_int) :: my_getdayofyearfromdatetime
      type(c_ptr), value :: dt
    end function my_getdayofyearfromdatetime
    !
    function my_getjuliandayfromdatetime(dt, jd) bind(c, name='getJulianDayFromDateTime')
#ifdef __SX__
      use, intrinsic :: iso_c_binding, only: c_ptr
#else
      import :: c_ptr
#endif
      type(c_ptr) :: my_getjuliandayfromdatetime
      type(c_ptr), value :: dt
      type(c_ptr), value :: jd
    end function my_getjuliandayfromdatetime
    !
  end interface
  !
contains
  !
  function newdatetimefromstring(string) result(ret_datetime)
    type(datetime), pointer :: ret_datetime
    character(len=*), intent(in) :: string
    type(c_ptr) :: c_pointer
    c_pointer = my_newdatetime(trim(string)//c_null_char)
    call c_f_pointer(c_pointer, ret_datetime)
  end function newdatetimefromstring
  !
  function newdatetimefromraw(year, month, day, hour, minute, second, ms) result(ret_datetime)
    type(datetime), pointer :: ret_datetime
    integer(c_int64_t), intent(in) :: year
    integer(c_int), intent(in) :: month
    integer(c_int), intent(in) :: day
    integer(c_int), intent(in) :: hour
    integer(c_int), intent(in) :: minute
    integer(c_int), intent(in) :: second
    integer(c_int), intent(in) :: ms
    type(c_ptr) :: c_pointer
    c_pointer = my_newrawdatetime(year, month, day, hour, minute, second, ms)
    call c_f_pointer(c_pointer, ret_datetime)
  end function newdatetimefromraw
  !
  function newdatetimefromconstructandcopy(src) result(dest)
    type(datetime), pointer :: dest
    type(datetime), pointer :: src
    type(c_ptr) :: c_pointer
    c_pointer = my_constructandcopydatetime(c_loc(src))
    call c_f_pointer(c_pointer, dest)
  end function newdatetimefromconstructandcopy
  !
  subroutine deallocateDatetime(my_datetime)
    type(datetime), pointer :: my_datetime
    call my_deallocatedatetime(c_loc(my_datetime))
    nullify(my_datetime)
  end subroutine deallocateDatetime
  !
  subroutine datetimeToString(my_datetime, string)
    type(datetime), pointer :: my_datetime
    character(len=max_datetime_str_len) :: string
    type(c_ptr) :: dummy_ptr
    integer :: i
    dummy_ptr = my_datetimetostring(c_loc(my_datetime), string)
    char_loop: do i = 1 , len(string)
      if (string(i:i) == c_null_char) exit char_loop
    end do char_loop
    string(i:len(string)) = ' '
  end subroutine datetimeToString
  !
  subroutine replacedatetime(dest, src)
    type(datetime), target, intent(inout) :: dest
    type(datetime), target, intent(in) :: src
    type(c_ptr) :: dummy_ptr
    dummy_ptr = my_replacedatetime(c_loc(src), c_loc(dest))
  end subroutine replacedatetime
  !
  function datetime_gt(op1, op2) result(gt)
    logical :: gt
    type(datetime), target, intent(in) :: op1
    type(datetime), target, intent(in) :: op2
    integer(c_int) :: ret
    ret = my_comparedatetime(c_loc(op1), c_loc(op2))
    if (ret > 0) then
      gt = .true.
    else
      gt = .false.
    endif
  end function datetime_gt
  !
  function datetime_lt(op1, op2) result(lt)
    logical :: lt
    type(datetime), target, intent(in) :: op1
    type(datetime), target, intent(in) :: op2
    integer(c_int) :: ret
    ret = my_comparedatetime(c_loc(op1), c_loc(op2))
    if (ret < 0) then
      lt = .true.
    else
      lt = .false.
    endif
  end function datetime_lt
  !
  function datetime_lt_or_eq(op1, op2) result(lt_or_eq)
    logical :: lt_or_eq
    type(datetime), target, intent(in) :: op1
    type(datetime), target, intent(in) :: op2
    integer(c_int) :: ret
    ret = my_comparedatetime(c_loc(op1), c_loc(op2))
    if (ret <= 0) then
      lt_or_eq = .true.
    else
      lt_or_eq = .false.
    endif
  end function datetime_lt_or_eq
  !
  function datetime_gt_or_eq(op1, op2) result(gt_or_eq)
    logical :: gt_or_eq
    type(datetime), target, intent(in) :: op1
    type(datetime), target, intent(in) :: op2
    integer(c_int) :: ret
    ret = my_comparedatetime(c_loc(op1), c_loc(op2))
    if (ret >= 0) then
      gt_or_eq = .true.
    else
      gt_or_eq = .false.
    endif
  end function datetime_gt_or_eq
  !
  function datetime_eq(op1, op2) result(eq)
    logical :: eq
    type(datetime), target, intent(in) :: op1
    type(datetime), target, intent(in) :: op2
    integer(c_int) :: ret
    ret = my_comparedatetime(c_loc(op1), c_loc(op2))
    if (ret == 0) then
      eq = .true.
    else
      eq = .false.
    endif
  end function datetime_eq
  !
  function datetime_ne(op1, op2) result(ne)
    logical :: ne
    type(datetime), target, intent(in) :: op1
    type(datetime), target, intent(in) :: op2
    integer(c_int) :: ret
    ret = my_comparedatetime(c_loc(op1), c_loc(op2))
    if (ret /= 0) then
      ne = .true.
    else
      ne = .false.
    endif
  end function datetime_ne
  !
  subroutine getJulianDayFromDatetime(dt, jd)
    type(datetime), target :: dt
    type(julianday), pointer :: jd
    type(c_ptr) :: dummy_ptr
    dummy_ptr = my_getjuliandayfromdatetime(c_loc(dt), c_loc(jd))
  end subroutine getJulianDayFromDatetime
  !
end module mtime_datetime
!>
!! @brief TimeDelta and some operations supported on TimeDelta.
!!
!! @details
!!
!! @author  Luis Kornblueh, Max Planck Institute for Meteorology
!! @author  Rahul Sinha, Max Planck Institute for Meteorology
!!
!___________________________________________________________________________________________________________
module mtime_timedelta
  !
  use, intrinsic :: iso_c_binding, only: c_int, c_int64_t, c_char, c_ptr, c_null_char, c_loc, c_f_pointer
  !
  use mtime_datetime
  !
  implicit none
  !
  private
  !
  public :: max_timedelta_str_len
  public :: timedelta
  public :: newTimedelta
  public :: deallocateTimedelta
  public :: timedeltaToString
  public :: operator(+)
  !
  !> provides a string length for toString 
  integer, parameter :: max_timedelta_str_len = 32
  !
  type, bind(c) :: timedelta
    character(c_char)  :: sign
    integer(c_int64_t) :: year
    integer(c_int)     :: month
    integer(c_int)     :: day
    integer(c_int)     :: hour
    integer(c_int)     :: minute
    integer(c_int)     :: second
    integer(c_int)     :: ms
  end type timedelta
  !
  interface operator (+)
    module procedure addtimedeltatodatetime
    module procedure adddatetimetotimedelta
  end interface operator (+)
  !
  interface
    !
    function my_newtimedelta(string) result(c_pointer) bind(c, name='newTimeDelta')
#ifdef __SX__
      use, intrinsic :: iso_c_binding, only: c_int  
#else
      import :: c_char, c_ptr
#endif
      type(c_ptr) :: c_pointer
      character(c_char), dimension(*) :: string
    end function my_newtimedelta
    !
    subroutine my_deallocatetimedelta(dt) bind(c, name='deallocateTimeDelta')
#ifdef __SX__
      use, intrinsic :: iso_c_binding, only: c_int  
#else
      import :: c_ptr
#endif
      type(c_ptr), value :: dt
    end subroutine my_deallocatetimedelta
    !
    function my_timedeltatostring(dt, tostring) result(string_ptr) bind(c, name='timedeltaToString')
#ifdef __SX__
      use, intrinsic :: iso_c_binding, only: c_int  
#else
      import :: c_ptr, c_char
#endif
      type(c_ptr) :: string_ptr
      type(c_ptr), value :: dt
      character(c_char), dimension(*) :: tostring
    end function my_timedeltatostring
    !
    function my_addtimedeltatodatetime(my_datetime, my_timedelta, ret_datetime) result(datetime_ptr) &
         &   bind(c, name='addTimeDeltaToDateTime')
#ifdef __SX__
      use, intrinsic :: iso_c_binding, only: c_int  
#else
      import :: c_ptr
#endif
      type(c_ptr), value :: my_datetime
      type(c_ptr), value :: my_timedelta
      type(c_ptr), value :: ret_datetime
      type(c_ptr) :: datetime_ptr
    end function my_addtimedeltatodatetime
    !
  end interface
  !
contains
  !
  function newTimedelta(string) result(ret_timedelta)
    type(timedelta), pointer :: ret_timedelta
    character(len=*), intent(in) :: string
    type(c_ptr) :: c_pointer
    c_pointer = my_newtimedelta(trim(string)//c_null_char)
    call c_f_pointer(c_pointer, ret_timedelta)
  end function newTimedelta
  !
  subroutine deallocateTimedelta(my_timedelta)
    type(timedelta), pointer :: my_timedelta
    call my_deallocatetimedelta(c_loc(my_timedelta))
    nullify(my_timedelta)
  end subroutine deallocateTimedelta
  !
  subroutine timedeltaToString(my_timedelta, string)
    type(timedelta), pointer :: my_timedelta
    character(len=max_timedelta_str_len) :: string
    type(c_ptr) :: dummy_ptr
    integer :: i
    dummy_ptr = my_timedeltatostring(c_loc(my_timedelta), string)
    char_loop: do i = 1 , len(string)
      if (string(i:i) == c_null_char) exit char_loop
    end do char_loop
    string(i:len(string)) = ' '
  end subroutine timedeltaToString
  !
  function addTimedeltaToDatetime(op1, op2) result(ret)
    type(datetime), target :: ret
    type(datetime), target, intent(in) :: op1
    type(timedelta), target, intent(in) :: op2
    type(c_ptr) :: dummy_ptr
    dummy_ptr = my_addtimedeltatodatetime(c_loc(op1), c_loc(op2), c_loc(ret))
  end function addTimedeltaToDatetime
  !
  function addDatetimeToTimedelta(op2, op1) result(ret)
    type(datetime), target :: ret
    type(datetime), target, intent(in) :: op1
    type(timedelta), target, intent(in) :: op2
    type(c_ptr) :: dummy_ptr
    dummy_ptr = my_addtimedeltatodatetime(c_loc(op1), c_loc(op2), c_loc(ret))
  end function addDatetimeToTimedelta
  !
end module mtime_timedelta
!>
!! @brief Definition of the basic event type and its methods.
!!
!! @details
!!
!___________________________________________________________________________________________________________
module mtime_events
  !
  use, intrinsic :: iso_c_binding, only: c_int64_t, c_char, c_null_char, c_bool, c_ptr, c_loc, c_f_pointer
  !
  implicit none
  !
  private
  !
  public :: max_eventname_str_len
  public :: event
  public :: newEvent
  public :: deallocateEvent
  public :: eventToString
  !
  !> provides a string length for toString 
  integer, parameter :: max_eventname_str_len = 132
  !
  type, bind(c) :: event
    type(c_ptr) :: nextEventInGroup
    integer(c_int64_t) :: eventId
    type(c_ptr) :: eventName
    type(c_ptr) :: eventReferenceDatetime
    type(c_ptr) :: eventFirstDatetime
    type(c_ptr) :: eventLastDatetime
    type(c_ptr) :: eventInterval
    logical(c_bool) :: triggerCurrentEvent
    logical(c_bool) :: nextEventIsFirst
    logical(c_bool) :: eventisFirstInDay
    logical(c_bool) :: eventisFirstInMonth
    logical(c_bool) :: eventisFirstInYear
    logical(c_bool) :: eventisLastInDay
    logical(c_bool) :: eventisLastInMonth
    logical(c_bool) :: eventisLastInYear
    type(c_ptr) :: triggerNextEventDateTime
    type(c_ptr) :: triggeredPreviousEventDateTime
  end type event
  !
  interface
    !
    function my_newevent(name, referenceDate, firstdate, lastDate, interval) result(c_pointer) bind(c, name='newEvent')
#ifdef __SX__
      use, intrinsic :: iso_c_binding, only: c_int  
#else
      import :: c_char, c_ptr
#endif
      type(c_ptr) :: c_pointer
      character(c_char), dimension(*) :: name
      character(c_char), dimension(*) :: referenceDate
      character(c_char), dimension(*) :: firstDate
      character(c_char), dimension(*) :: lastDate
      character(c_char), dimension(*) :: interval
    end function my_newevent
    !
    subroutine my_deallocateevent(ev) bind(c,name='deallocateEvent')
#ifdef __SX__
      use, intrinsic :: iso_c_binding, only: c_int  
#else
      import :: c_ptr
#endif
      type(c_ptr), value :: ev
    end subroutine my_deallocateevent
    !
    function my_eventtostring(my_event, string) result(string_ptr) bind(c, name='eventToString')
#ifdef __SX__
      use, intrinsic :: iso_c_binding, only: c_int  
#else
      import :: c_ptr, c_char
#endif
      type(c_ptr) :: string_ptr
      type(c_ptr), value :: my_event
      character(c_char), dimension(*) :: string
    end function my_eventtostring
    !
  end interface
  !
contains
  !
  function newEvent(name, referenceDate, firstdate, lastDate, interval) result(ret_event)
    type(event), pointer :: ret_event
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: referenceDate
    character(len=*), intent(in) :: firstDate
    character(len=*), intent(in) :: lastDate
    character(len=*), intent(in) :: interval
    type(c_ptr) :: c_pointer
    c_pointer = my_newevent(trim(name)//c_null_char,               &
         &                  trim(referenceDate)//c_null_char,      &
         &                  trim(firstDate)//c_null_char,          &
         &                  trim(lastDate)//c_null_char,           &
         &                  trim(interval)//c_null_char)
    call c_f_pointer(c_pointer, ret_event)
  end function newEvent
  !
  subroutine deallocateEvent(my_event)
    type(event), pointer :: my_event
    call my_deallocateevent(c_loc(my_event))
    my_event => null()
  end subroutine deallocateEvent
  !
  subroutine eventToString(my_event, string)
    type(event), pointer :: my_event
    character(len=max_eventname_str_len) :: string
    type(c_ptr) :: dummy_ptr
    integer :: i
    dummy_ptr = my_eventtostring(c_loc(my_event), string)
    char_loop: do i = 1 , len(string)
      if (string(i:i) == c_null_char) exit char_loop
    end do char_loop
    string(i:len(string)) = ' '
  end subroutine eventToString
  !
end module mtime_events
!>
!! @brief Event-groups which contains a list of events.
!!
!! @details
!!
!! @author  Luis Kornblueh, Max Planck Institute for Meteorology
!! @author  Rahul Sinha, Max Planck Institute for Meteorology
!!
!___________________________________________________________________________________________________________
module mtime_eventgroups
  !
  use, intrinsic :: iso_c_binding, only: c_int64_t, c_ptr, c_char, c_null_char, c_bool, c_loc, c_f_pointer, &
       &            c_associated
  !
  use mtime_events
  !
  implicit none
  !
  private
  !
  public :: max_groupname_str_len
  public :: eventgroup
  public :: newEventGroup
  public :: deallocateEventGroup
  public :: addEventToEventGroup
  public :: removeEventFromEventGroup
  public :: getFirstEventFromEventGroup
  public :: getNextEventFromEventGRoup
  !
  !> provides a string length for toString 
  integer, parameter :: max_groupname_str_len = 132
  !
  type, bind(c) :: eventgroup
    integer(c_int64_t) :: eventGroupId
    type(c_ptr) :: eventGroupName
    type(c_ptr) :: firstEventInGroup
  end type eventgroup
  !
  interface
    !
    function my_neweventgroup(name) result(c_pointer) bind(c, name='newEventGroup')
#ifdef __SX__
      use, intrinsic :: iso_c_binding, only: c_int  
#else
      import :: c_char, c_ptr
#endif
      type(c_ptr) :: c_pointer
      character(c_char), dimension(*) :: name
    end function my_neweventgroup
    !
    subroutine my_deallocateeventgroup(evgrp) bind(c,name='deallocateEventGroup')
#ifdef __SX__
      use, intrinsic :: iso_c_binding, only: c_int  
#else
      import :: c_ptr
#endif
      type(c_ptr), value :: evgrp
    end subroutine my_deallocateeventgroup
    !
    function my_addeventtoeventgroup(my_event, my_eventgroup) result(ret) bind(c, name='addNewEventToEventGroup')
#ifdef __SX__
      use, intrinsic :: iso_c_binding, only: c_int  
#else
      import :: c_bool, c_ptr
#endif
      logical(c_bool) :: ret
      type(c_ptr), value :: my_event
      type(c_ptr), value :: my_eventgroup
    end function my_addeventtoeventgroup
    !
    function my_removeeventfromeventgroup(evname, evgrp) result(ret) bind(c, name='removeEventFromEventGroup')
#ifdef __SX__
      use, intrinsic :: iso_c_binding, only: c_int  
#else
      import :: c_bool, c_char, c_ptr
#endif
      logical(c_bool) :: ret
      character(c_char), dimension(*) :: evname
      type(c_ptr), value :: evgrp
    end function my_removeeventfromeventgroup
    !
  end interface
  !
contains
  !
  function newEventGroup(name) result(ret_eventgroup)
    type(eventgroup), pointer :: ret_eventgroup
    character(len=*), intent(in) :: name
    type(c_ptr) :: c_pointer
    c_pointer = my_neweventgroup(trim(name)//c_null_char)
    call c_f_pointer(c_pointer, ret_eventgroup)
  end function newEventGroup
  !
  subroutine deallocateEventGroup(my_eventgroup)
    type(eventgroup), pointer :: my_eventgroup
    call my_deallocateeventgroup(c_loc(my_eventgroup))
    my_eventgroup => null()
  end subroutine deallocateEventGroup
  !
  function addEventToEventGroup(my_event, my_eventGroup) result(ret)
    logical :: ret
    type(event), pointer :: my_event
    type(eventgroup), pointer :: my_eventgroup
    ret = my_addeventtoeventgroup(c_loc(my_event), c_loc(my_eventgroup))
  end function addEventToEventGroup
  !
  function removeEventfromEventGroup(my_name, my_eventGroup) result(ret)
    logical :: ret
    character(len=*), intent(in) :: my_name
    type(eventgroup), pointer :: my_eventgroup
    ret = my_removeeventfromeventgroup(trim(my_name)//c_null_char, c_loc(my_eventgroup))
  end function removeEventFromEventGroup
  !
  function getFirstEventFromEventGroup(my_eventgroup) result(ret_event)
    type(event), pointer :: ret_event
    type(eventgroup), pointer :: my_eventgroup
    type(c_ptr) :: c_pointer
    c_pointer = my_eventgroup%firstEventInGroup
        if (c_associated(c_pointer)) then
      call c_f_pointer(c_pointer, ret_event)
    else
      ret_event => null()
    endif
  end function getFirstEventFromEventGroup
  !
  function getNextEventFromEventGroup(my_event) result(ret_event)
    type(event), pointer :: ret_event
    type(event), pointer :: my_event
    type(c_ptr) :: c_pointer
    c_pointer = my_event%nextEventInGroup
    if (c_associated(c_pointer)) then
      call c_f_pointer(c_pointer, ret_event)
    else
      ret_event => null()
    endif
  end function getNextEventFromEventGroup
  !
end module mtime_eventgroups
!>
!! @}
!>
!! @brief mtime is a compound module making all library components accessible via one module.
!!
!! @details
!!
!! @author  Luis Kornblueh, Max Planck Institute for Meteorology
!! @author  Rahul Sinha, Max Planck Institute for Meteorology
!!
!___________________________________________________________________________________________________________
module mtime
  !
  use, intrinsic :: iso_c_binding, only: c_int, c_int64_t, c_bool, c_ptr, c_char
  !
  use mtime_calendar
  use mtime_julianday
  use mtime_date
  use mtime_time
  use mtime_datetime
  use mtime_timedelta
  use mtime_events
  use mtime_eventgroups
  !
  implicit none
  !
  public
  !
end module mtime
