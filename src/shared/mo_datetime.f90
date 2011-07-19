!>
!! This modules contains a data structure for calendar date and timeinformation
!!
!! Time information generally depends on the type of calendar in use, and
!! is needed in the different forms. This modules provides the data type
!! "t_datetime" and methods to handle time information in the model.
!! 
!! Calendars: 
!! - Historic Julian/Gregorian based on Julian day (*)
!!   Julian calendar until   1582-10-04,
!!     with a leap year if MOD(year,4)=0
!!   Gregorian calendar from 1582-10-15,
!!     with a leap year if (MOD(year,4)=0 and MOD(year,100)/=0) or (MOD(year,400)=0)
!!   Dates from 1582-10-5 to 1582-10-14 are skipped
!!
!! - Proleptic Gregorian based on Julian day (*)
!!   Gregorian calendar used for all times,
!!     with a leap year if (MOD(year,4)=0 and MOD(year,100)/=0) or (MOD(year,400)=0)
!!
!! - 360 day/year calendar based on time since 0000-01-01, 00:00 Universal Time
!!     with all months equal 30 days
!!
!! (*) Julian day is the time elapsed since  (-4712)-01-01 12:00 Universal Time
!!                                         =(4713BC)-01-01 12:00 Universal Time
!!
!! The data structure contains information on:
!! - calendar type
!! - date in year, month and day
!! - time in hour, minute, second
!! - elapsed number of days and fraction of day since the reference point of time
!! - length of day, month and year
!! - time in the year, fraction of the year, day of the year
!! - time in the month, fraction of the month
!! - time in the day (=fraction of the day), seconds in the day
!!
!! Methods are defined to:
!! - check if the date and time information is valid
!! - convert (year,month,day,hour,minute,second) to (day,time) of the calendar
!! - convert (day,time) of the calendar to (year,month,day,hour,minute,second)
!! - compute auxiliary date and time information
!! - add a time increment (seconds,minutes,hours,days) (positive or negative)
!! - print date and time information
!!
!! References:
!!   Montenbruck, Oliver, "Practical Ephmeris Calculations", Ch. 2, pp 33-40.
!!   The 1992 Astronomical Almanac, page B4.
!!
!!
!! @author Marco Giorgetta, MPI-M Hamburg
!!
!!
!! @par Revision History
!! First version by Marco Giorgetta, MPI-M Hamburg (2010-11-19)
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
MODULE mo_datetime

  USE mo_kind,      ONLY: wp, i8
  USE mo_exception, ONLY: finish, message, message_text

  IMPLICIT NONE

  PRIVATE

  ! Date and time structure and methods
  !
  PUBLIC :: &
    ! types
    &   t_datetime          ,& ! type for date and time information
    ! parameters
    &   julian_gregorian    ,& ! historic Julian / Gregorian calendar
    &   proleptic_gregorian ,& ! proleptic Gregorian calendar
    &   cly360              ,& ! constant 30 dy/mo and 360 dy/yr calendar
    &   idaylen, rdaylen    ,& ! length of day in seconds, as integer and real
    ! subroutines
    &   check_date          ,& ! check if (yr,mo,dy,hr,mn,s) is valid
    &   date_to_time        ,& ! convert (yr,mo,dy,hr,mn,s) to time and compute auxiliary info
    &   time_to_date        ,& ! convert time to (yr,mo,dy,hr,mn,s) and compute auxiliary info
    &   add_time            ,& ! add a time increment in (dy,hr,mn,s) to a datetime variable
    &   print_datetime_all  ,& ! print a datetime variable, including: 
    &   print_calendar      ,& ! print calendar info
    &   print_datetime      ,& ! print date and time
    &   print_auxinfos      ,& ! print auxiliary date and time info
    &   iso8601                ! convert (yr,mo,dy,hr,mn,s) to a string of iso8601 format

 PUBLIC :: string_to_datetime, datetime_to_string,  date_len

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  ! Calendar types
  !
  INTEGER,  PARAMETER :: julian_gregorian    = 0 !< historic Julian / Gregorian
  INTEGER,  PARAMETER :: proleptic_gregorian = 1 !< proleptic Gregorian
  INTEGER,  PARAMETER :: cly360              = 2 !< constant 30 dy/mo and 360 dy/yr

  ! Day length
  !
  INTEGER,  PARAMETER :: idaylen=86400     ! [s]
  REAL(wp), PARAMETER :: rdaylen=86400._wp ! [s]

  INTEGER, PARAMETER :: date_len = 32
  !>
  !! Date and time structure
  !!
  !! Provides components for time information in different forms
  !!
  TYPE :: t_datetime
    !
    ! 1. calendar
    ! -----------
    INTEGER     :: calendar  ! 0 = historic Julian/Gregorian calendar
    !                        ! 1 = proleptic Gregorian calendar
    !                        ! 2 = constant 30 dy/mo and 360 dy/yr calendar
    !
    !------------------------------------------------
    !
    ! 2. time information
    ! -------------------
    !
    ! date and time
    !
    INTEGER     :: year      ! [yr]
    INTEGER     :: month     ! [mo]  1:12
    INTEGER     :: day       ! [dy]  1:{28,29,30,31}
    INTEGER     :: hour      ! [hr]  0:23
    INTEGER     :: minute    ! [mn]  0:59
    REAL(wp)    :: second    ! [s]   [0.,60.[
    !
    ! elapsed time since reference point of time
    !
    INTEGER(i8) :: calday    ! [dy] full days elapsed since the reference
    !                        !      point of time of the calendar
    REAL(wp)    :: caltime   ! [dy] remaining fraction of a day elapsed ...
    !
    !------------------------------------------------
    !
    ! 3. auxiliary information
    !
    ! lengths
    INTEGER     :: yealen  !< [dy]  length of year
    INTEGER     :: monlen  !< [dy]  length of month
    INTEGER     :: daylen  !< [s]   length of day
    !
    ! time in year
    REAL(wp)    :: yeatim  !< [dy]  time since year-01-01, 00:00 UT
    REAL(wp)    :: yeafrc  !< []    same, but as fraction of year length
    INTEGER     :: yeaday  !< [dy]  current day of the year
    !
    ! time in month
    REAL(wp)    :: montim  !< [dy]  time since the year-month-01, 00:00 UT
    REAL(wp)    :: monfrc  !< []    same, but as fraction of the month length
    !
    ! time in day
    REAL(wp)    :: daytim  !< [dy]  time since year-month-day, 00:00 UT
    REAL(wp)    :: daysec  !< [s]   same, but in [s]
    !
  END TYPE t_datetime

CONTAINS

  !>
  !! Subroutine to check if date and time information is valid.
  !!
  !! Check if:
  !! - the calendar index is known
  !! - month and day are valid in this calendar
  !! - hour, minute and second is valid
  !!
  !! @par Revision History
  !! First version by Marco Giorgetta, MPI-M Hamburg (2010-10-15)
  !!
  SUBROUTINE check_date(datetime)

    TYPE(t_datetime), INTENT(inout) :: datetime
    !
    ! input components
    ! - calendar
    ! - year, month, day, hour, minute, second
    !
    !
    ! output components
    ! - monlen
    !

    ! Length of month
    !
    SELECT CASE (datetime%calendar)
      !
    CASE (julian_gregorian, proleptic_gregorian)
      CALL julian_monlen(datetime)
      !
    CASE (cly360)
      datetime%monlen = 30
      !
    CASE DEFAULT
      CALL print_datetime(datetime)
      CALL finish ('mo_datetime/check_datetime','invalid calendar')
      !
    END SELECT

    ! Check for historic Julian Gregorian calendar if the date is
    ! in 5-14 October 1582 AD, which is not valid.
    IF ( datetime%calendar == julian_gregorian &
      &  .AND. datetime%year  == 1582          &
      &  .AND. datetime%month == 10            &
      &  .AND. datetime%day   >  4             &
      &  .AND. datetime%day   <  15 ) THEN
      CALL print_datetime(datetime)
      CALL finish ('mo_datetime/check_datetime','invalid date in 5-14 Oct 1582')
    END IF

    ! Check month
    IF ( datetime%month < 1 .OR. datetime%month > 12) THEN
      CALL print_datetime(datetime)
      CALL finish ('mo_datetime/check_datetime','invalid month')
    END IF

    ! Check day
    IF ( datetime%day < 1 .OR. datetime%day > datetime%monlen) THEN
      CALL print_datetime(datetime)
      CALL finish ('mo_datetime/check_datetime','invalid day')
    END IF

    ! Check time
    IF (     datetime%hour   < 0     .OR. datetime%hour   >  23       &
      & .OR. datetime%minute < 0     .OR. datetime%minute >  59       &
      & .OR. datetime%second < 0._wp .OR. datetime%second >= 60._wp ) THEN
      CALL print_datetime(datetime)
      CALL finish ('mo_datetime/check_datetime','invalid time')
    END IF

  END SUBROUTINE check_date

  !----------------------------------------------------------------------------

  !>
  !! Subroutine to convert a date given as (year,month,day,hour,minute,second)
  !! to elapsed time in days since the reference point of time.
  !!
  !! Dependent on the calendar specialzed subrotines are called to:
  !! - check the validity of (year,month,day,hour,minute,second)
  !! - to compute the time elapsed since the reference point of time
  !! - and to compute auxiliary time information.
  !!
  !! @par Revision History
  !! First version by Marco Giorgetta, MPI-M Hamburg (2010-10-14)
  !!
  SUBROUTINE date_to_time(datetime)

    TYPE(t_datetime), INTENT(inout) :: datetime
    !
    ! input components
    ! - calendar
    ! - year, month, day, hour, minute, second
    !
    ! output component
    ! - calday, caltime
    ! - yealen, monlen, daylen
    ! - yeatim, yeafrc, yeaday
    ! - montim, monfrc
    ! - daytim, daysec
    ! 

    ! Check date and time
    !
    CALL check_date(datetime)

    ! Compute time with respect to refrence point of time
    !
    SELECT CASE (datetime%calendar)
      !
    CASE (julian_gregorian, proleptic_gregorian)
      ! compute time based on Julian day 0 at (-4712)-01-01, 12:00 UT
      CALL date_to_julianday(datetime)
      !
    CASE (cly360)
      ! compute time based on 0 at 0000-01-01, 00:00 UT
      CALL date_to_cly360day(datetime)
      !
    END SELECT

    ! Compute auxiliary date and time information
    !
    CALL aux_datetime(datetime)

  END SUBROUTINE date_to_time

  !----------------------------------------------------------------------------

  !>
  !! Subroutine to convert time given as (year,month,day,hour,minute,second)
  !! to time based on Julian day 0 at (-4712)-01-01, 12:00.
  !!
  !! For the historic Julian/Gregorian and the proleptic Gregorian calendar.
  !!
  !! @par Revision History
  !! Based on mo_datetime.f90 by Marco Giorgetta, MPI-M Hamburg (2010-10-15)
  !!
  SUBROUTINE date_to_julianday(datetime)

    TYPE(t_datetime), INTENT(inout) :: datetime
    !
    ! input components
    ! - calendar
    ! - year, month, day, hour, minute, second
    !
    ! output component
    ! - calday, caltime
    !

    INTEGER  :: iy, im    !< modified year and month indices
    INTEGER  :: ib        !< leap year correction
    REAL(wp) :: time      !< calendar time in days

    ! Shift the calendar year by 2 months, i.e. with start on the 1st of March,
    ! so that a leap day is the last day of the shifted year
    IF (datetime%month <= 2) THEN
      iy = datetime%year  - 1
      im = datetime%month + 12
    ELSE
      iy = datetime%year
      im = datetime%month
    ENDIF
    !
    ! Correction for leap day
    !
    ! (1)  (proleptic) Gregorian calendar
    !
    IF ( iy < 0 ) THEN
      ib = INT((iy+1)/400)-INT((iy+1)/100)
    ELSE
      ib = INT(iy/400)-INT(iy/100)
    END IF
    !
    ! (2) Julian calendar until 1582-10-4
    !
    IF ( datetime%calendar == julian_gregorian &
      &  .AND. datetime%year <  1582  &
      &  .OR. (datetime%year == 1582 .AND. datetime%month <  10 ) &
      &  .OR. (datetime%year == 1582 .AND. datetime%month == 10 .AND. datetime%day <= 4 ) ) THEN
      !
      ib = -2
      !
    END IF
    !
    ! Compute Julian day
    !
    time = REAL(FLOOR(365.25_wp*REAL(iy,wp))+INT(30.6001_wp*REAL(im+1,wp))+ib+datetime%day ,wp) &
      & +(datetime%second+REAL(60*(datetime%minute+60*datetime%hour),wp))/rdaylen               &
      & +1720996.5_wp

    datetime%calday  = INT(time,i8)
    datetime%caltime = MOD(time,1._wp)

  END SUBROUTINE date_to_julianday

  !----------------------------------------------------------------------------

  !>
  !! Subroutine to convert time given as (year,month,day,hour,minute,second)
  !! to time since 0000-01-01, 00:00 for constant years of 360 days.
  !!
  !! This subroutine is used for the historic constant 360 days per year calendar.
  !!
  !! @par Revision History
  !! First version by Marco Giorgetta, MPI-M Hamburg (2010-10-14)
  !!
  SUBROUTINE date_to_cly360day(datetime)

    TYPE(t_datetime), INTENT(inout)  :: datetime
    !
    ! input components
    ! - calendar
    ! - year, month, day, hour, minute, second
    !
    ! output component
    ! - calday, caltime
    !

    REAL(wp) :: time     !< calendar time in days

    ! Time in days since 0 = 0000-01-01, 00:00 UT
    time = REAL(360*datetime%year + 30*(datetime%month-1) + (datetime%day-1),wp)  &
      &    +(datetime%second + REAL(60*(datetime%minute + 60*datetime%hour),wp))  &
      &    /rdaylen

    datetime%calday  = INT(time,i8)
    datetime%caltime = MOD(time,1._wp)

  END SUBROUTINE date_to_cly360day

  !----------------------------------------------------------------------------

  !>
  !! Subroutine to convert a time in days since the reference point of time
  !! to a date given as (year,month,day,hour,minute,second).
  !!
  !! Dependent on the calendar specialzed subrotines are called
  !! to compute the date information.
  !!
  !! @par Revision History
  !! First version by Marco Giorgetta, MPI-M Hamburg (2010-10-14)
  !!
  SUBROUTINE time_to_date(datetime)

    TYPE(t_datetime), INTENT(inout) :: datetime
    !
    ! input components
    ! - calendar
    ! - calday, caltime
    !
    ! output component
    ! - year, month, day, hour, minute, second
    ! - yealen, monlen, daylen
    ! - yeatim, yeafrc, yeaday
    ! - montim, monfrc
    ! - daytim, daysec
    !

    SELECT CASE (datetime%calendar)
      !
    CASE (julian_gregorian,proleptic_gregorian)
      ! based on Julian day, 0 at (-4712)-01-01, 12:00
      CALL julianday_to_date(datetime)
      !
    CASE (cly360)
      ! based on 0 at 0000-01-01, 00:00
      CALL cly360day_to_date(datetime)
      !
    END SELECT

    ! Check date and time
    !
    CALL check_date(datetime)

    ! Compute auxiliary date and time information
    !
    CALL aux_datetime(datetime)

  END SUBROUTINE time_to_date

  !----------------------------------------------------------------------------

  !>
  !! Subroutine to convert a time given as Julian day
  !! to a date given as (year,month,day,hour,minute,second)
  !!
  !! This subroutine is used for the historic Julian/Gregorian
  !! or the proleptic Gregorian calendar.
  !!
  !! @par Revision History
  !! Coded based on mo_time_base.f90 of echam6 by Marco Giorgetta, MPI-M Hamburg (2010-10-14)
  !!
  SUBROUTINE julianday_to_date(datetime)

    TYPE(t_datetime), INTENT(inout) :: datetime
    !
    ! input components
    ! - calendar
    ! - calday, caltime
    !
    ! output component
    ! - year, month, day, hour, minute, second

    REAL(wp) :: za, zb, zc, zd, ze, zf
    REAL(wp) :: time, daytim

    time = REAL(datetime%calday,wp)+datetime%caltime

    za = REAL(FLOOR(time+0.5_wp),wp)
    zb = REAL(FLOOR((za-1867216.25_wp)/36524.25_wp),wp)

    IF ( (datetime%calendar==julian_gregorian) .AND. (za < 2299161.0_wp) ) THEN
      ! Julian/Gregorian until 1582-10-04, 12:00 UT =  Julian day 2299161.0
      zc = za+1524.0_wp
    ELSE
      ! Julian/Gregorian after 1582-10-04, 12:00 UT and proleptic Gregorian
      zc = za+zb-REAL(FLOOR(zb/4._wp),wp)+1525._wp
    END IF

    zd = REAL(FLOOR((zc-122.1_wp)/365.25_wp),wp)
    ze = REAL(FLOOR(365.25_wp*zd),wp)
    zf = REAL(FLOOR((zc-ze)/30.6001_wp),wp)

    datetime%day    = INT(zc-ze- REAL(FLOOR(30.6001_wp*zf),wp))
    datetime%month  = INT(zf-REAL(1+12*FLOOR(zf/REAL(14,wp)),wp))
    datetime%year   = INT(zd-4715._wp-REAL((7+datetime%month)/10,wp))

    daytim          = MODULO(time-0.5_wp,1._wp)
    datetime%hour   = INT(daytim*24._wp)
    datetime%minute = INT(MODULO(daytim*1440._wp,60._wp))
    datetime%second = MODULO(MODULO(daytim*rdaylen,3600._wp),60._wp)

  END SUBROUTINE julianday_to_date

  !----------------------------------------------------------------------------

  !>
  !! Subroutine to convert a time given as day in the 360 day/year calendar
  !! to a date given as (year,month,day,hour,minute,second)
  !!
  !! This subroutine is used for the constant 30 day/month and 360 day/year calendar
  !!
  !! @par Revision History
  !! Coded based on mo_time_base.f90 of echam6 by Marco Giorgetta, MPI-M Hamburg (2010-10-14)
  !!
  SUBROUTINE cly360day_to_date(datetime)

    TYPE(t_datetime), INTENT(inout) :: datetime
    !
    ! input components
    ! - calendar
    ! - time
    !
    ! output component
    ! - year, month, day, hour, minute, second

    REAL(wp) :: time, daytim

    time               = REAL(datetime%calday,wp)+datetime%caltime

    datetime%year      = FLOOR(time/360._wp)
    datetime%month     = INT(MODULO(time,360._wp)/30._wp)+1
    datetime%day       = INT(MODULO(time,30._wp))+1

    daytim             = MODULO(time,1._wp)
    datetime%hour      = INT(daytim*24._wp)
    datetime%minute    = INT(MODULO(daytim*1440._wp,60._wp))
    datetime%second    = MODULO(MODULO(daytim*rdaylen,3600._wp),60._wp)

  END SUBROUTINE cly360day_to_date

  !----------------------------------------------------------------------------

  !>
  !! Subroutine to compute the day number in a year.
  !!
  !! @par Revision History
  !! Coded based on mo_time_base.f90 of echam6 by Marco Giorgetta, MPI-M Hamburg (2010-10-14)
  !!
  SUBROUTINE julian_yeaday(datetime)

    TYPE(t_datetime), INTENT(inout) :: datetime
    !
    ! input components
    ! - calendar
    ! - calday, caltime
    !
    ! output components
    ! - yeaday

    TYPE(t_datetime) :: year_01_01_00ut

    ! Set 
    year_01_01_00ut        = datetime
    year_01_01_00ut%month  = 1
    year_01_01_00ut%day    = 1
    year_01_01_00ut%hour   = 0
    year_01_01_00ut%minute = 0
    year_01_01_00ut%second = 0.0_wp

    CALL date_to_julianday(year_01_01_00ut)

    datetime%yeaday = INT( REAL(datetime%calday,wp)       + datetime%caltime           &
      &                   -REAL(year_01_01_00ut%calday,wp)- year_01_01_00ut%caltime )+1

  END SUBROUTINE julian_yeaday

  !----------------------------------------------------------------------------

  !>
  !! Subroutine to compute the length of the year in days
  !!
  !! This subroutine is used for the historic Julian/Gregorian
  !! or the proleptic Gregorian calendar.
  !! 
  !! @par Revision History
  !! Coded based on mo_time_base.f90 of echam6 by Marco Giorgetta, MPI-M Hamburg (2010-10-14)
  !!
  SUBROUTINE julian_yealen(datetime)

    TYPE(t_datetime), INTENT(inout) :: datetime
    !
    ! input components
    ! - calendar
    ! - year, month, day, hour, minute, second
    !
    ! output component
    ! - yealen
    ! 

    SELECT CASE (datetime%calendar)
      !
    CASE (julian_gregorian)
      IF (datetime%year == 1582) THEN
        datetime%yealen = 355
      ELSE IF ( (MOD(datetime%year,4)==0 .AND. MOD(datetime%year,100)/=0) &
        &       .OR. MOD(datetime%year,400)==0 ) THEN
        datetime%yealen = 366
      ELSE
        datetime%yealen = 365
      END IF
      !
    CASE (proleptic_gregorian)
      IF ( (MOD(datetime%year,4)==0 .AND. MOD(datetime%year,100)/=0) &
        &  .OR. MOD(datetime%year,400)==0 ) THEN
        datetime%yealen = 366
      ELSE
        datetime%yealen = 365
      END IF
      !
    END SELECT

  END SUBROUTINE julian_yealen

  !----------------------------------------------------------------------------

  !>
  !! Subroutine to compute the length of the month in days
  !!
  !! This subroutine is used for the historic Julian/Gregorian
  !! or the proleptic Gregorian calendar.
  !!
  !! @par Revision History
  !! Coded based on mo_time_base.f90 of echam6 by Marco Giorgetta, MPI-M Hamburg (2010-10-14)
  !!
  SUBROUTINE julian_monlen(datetime)

    TYPE(t_datetime), INTENT(inout) :: datetime
    !
    ! input components
    ! - calendar
    ! - year, month, day, hour, minute, second
    !
    ! output component
    ! - monlen
    ! 

    SELECT CASE(datetime%month)
      !
    CASE(1,3,5,7,8,12)
      datetime%monlen = 31
      !
    CASE(4,6,9,11)
      datetime%monlen = 30
      !
    CASE(10)
      IF (datetime%calendar == julian_gregorian .AND. datetime%year == 1582) THEN
        datetime%monlen = 21
      ELSE ! proleptic_gregorian
        datetime%monlen = 31
      END IF
      !
    CASE(2)
      ! Gregorian leap year
      IF ( (MOD(datetime%year,4)==0 .AND. MOD(datetime%year,100)/=0) &
        &  .OR. MOD(datetime%year,400)==0 ) THEN
        datetime%monlen = 29
      ELSE
        datetime%monlen = 28
      END IF
      ! Julian
      IF ( datetime%calendar == julian_gregorian .AND. datetime%year <=1582 ) THEN
        IF ( MOD(datetime%year,4) == 0 ) THEN
          datetime%monlen = 29
        ELSE
          datetime%monlen = 28
        END IF
      END IF
      !
    END SELECT

  END SUBROUTINE julian_monlen

  !----------------------------------------------------------------------------

  !>
  !! Subroutine to provide auxiliary time information depedent on the
  !! primary date and time information.
  !!
  !! Auxiliary components
  !! - length of year, month and day
  !! - 
  !! - month and day are valid in this calendar
  !! - hour, minute and second is valid
  !!
  !! @par Revision History
  !! First version by Marco Giorgetta, MPI-M Hamburg (2010-10-15)
  !!
  SUBROUTINE aux_datetime(datetime)

    TYPE(t_datetime), INTENT(inout) :: datetime
    !
    ! input components
    ! - calendar
    ! - year, month, day, hour, minute, second
    ! - calday, caltime
    !
    ! output components
    ! - yealen, monlen, daylen
    ! - yeaday
    ! - yeatim, yeafrc
    ! - montim, monfrc
    ! - daytim, daysec
    !

    ! Length of year, month and  day, and
    ! day in year
    !
    SELECT CASE (datetime%calendar)
      !
    CASE (julian_gregorian, proleptic_gregorian)
      CALL julian_yealen(datetime) ! [dy]
      CALL julian_monlen(datetime) ! [dy]
      datetime%daylen = idaylen    ! [s]
      !
      CALL julian_yeaday(datetime)

      !
    CASE (cly360)
      datetime%yealen =    360     ! [dy]
      datetime%monlen =     30     ! [dy]
      datetime%daylen = idaylen    ! [s]
      !
      datetime%yeaday = (datetime%month-1)*datetime%monlen + datetime%day
      !
    END SELECT

    ! Time in day and fraction of day
    datetime%daysec = (datetime%second + REAL(60*(datetime%minute+60*datetime%hour),wp))
    datetime%daytim =  datetime%daysec / rdaylen

    ! Time in month and fraction of month
    datetime%montim = REAL(datetime%day-1,wp) + datetime%daytim
    datetime%monfrc = datetime%montim / REAL(datetime%monlen,wp)

    ! Time in year and fraction of year
    datetime%yeatim = REAL(datetime%yeaday-1,wp) + datetime%daytim
    datetime%yeafrc = datetime%yeatim / REAL(datetime%yealen,wp)


  END SUBROUTINE aux_datetime

  !----------------------------------------------------------------------------

  !>
  !! Subroutine to add a time increment to a date
  !!
  !! @par Revision History
  !! First version by Marco Giorgetta, MPI-M Hamburg (2010-10-17)
  !!
  SUBROUTINE add_time(seconds,minutes,hours,days,datetime)

    REAL(wp)         , INTENT(in)    :: seconds
    INTEGER          , INTENT(in)    :: minutes,hours,days

    TYPE (t_datetime), INTENT(inout) :: datetime
    !
    ! input components
    ! - time
    !
    ! output components
    ! - year, month, day, hour, minute, second
    ! - time
    ! - yealen, monlen, daylen
    ! - yeatim, yeafrc, yeaday
    ! - montim, monfrc
    ! - daytim, daysec
    !

    REAL(wp) :: timesum

    timesum = datetime%caltime + (seconds+REAL(60*(minutes+60*hours),wp))/rdaylen

    datetime%caltime = MODULO(timesum,1._wp)
    datetime%calday  = datetime%calday + INT(days + FLOOR(timesum),i8)

    CALL time_to_date(datetime)
    CALL check_date  (datetime)
    CALL aux_datetime(datetime)

  END SUBROUTINE add_time

  !----------------------------------------------------------------------------

  !>
  !! Subroutine to print information on the calendar and time base.
  !!
  !! @par Revision History
  !! First version by Marco Giorgetta, MPI-M Hamburg (2010-10-17)
  !!
  SUBROUTINE print_datetime_all(datetime)

    TYPE (t_datetime), INTENT(in) :: datetime
    !
    ! input components
    ! - calendar
    ! - year, month, day, hour, minute, second
    ! - calday, caltime
    ! - yealen, monlen, daylen
    ! - yeatim, yeafrc, yeaday
    ! - montim, monfrc
    ! - daytim, daysec
    !

    CALL print_calendar(datetime)
    CALL print_datetime(datetime)
    CALL print_auxinfos(datetime)

  END SUBROUTINE print_datetime_all

  !----------------------------------------------------------------------------

  !>
  !! Subroutine to print information on the calendar and time base.
  !!
  !! @par Revision History
  !! First version by Marco Giorgetta, MPI-M Hamburg (2010-10-17)
  !!
  SUBROUTINE print_calendar(datetime)

    TYPE (t_datetime), INTENT(in) :: datetime
    !
    ! input components
    ! - calendar
    !

    WRITE(message_text,'(a,i1)') &
      &   'Calendar index : ',datetime%calendar
    CALL message('mo_datetime/print_calendar', message_text)
    !
    SELECT CASE (datetime%calendar)
    CASE (julian_gregorian)
      WRITE(message_text,'(a)')  &
        & 'Calendar       : Julian Gregorian'
      CALL message('mo_datetime/print_calendar', message_text)
      !
      WRITE(message_text,'(a)')  &
        & 'Time Base      : (-4712)-01-01, 12:00 UT'
      CALL message('mo_datetime/print_calendar', message_text)
      !
    CASE (proleptic_gregorian)
      WRITE(message_text,'(a)')  &
        & 'Calendar       : Proleptic Gregorian'
      CALL message('mo_datetime/print_calendar', message_text)
      !
      WRITE(message_text,'(a)')  &
        & 'Time Base      : (-4712)-01-01, 12:00 UT'
      CALL message('mo_datetime/print_calendar', message_text)
      !
    CASE (cly360)
      WRITE(message_text,'(a)')  &
        & 'Calendar       : Constant 30 dy/mo and 360 dy/yr'
      CALL message('mo_datetime/print_calendar', message_text)
      !
      WRITE(message_text,'(a)')  &
        & 'Time Base      : 0000-01-01, 00:00 UT'
      CALL message('mo_datetime/print_calendar', message_text)
      !
    END SELECT
    !
    WRITE(message_text,'(a)') ' '
    CALL message('mo_datetime/print_calendar', message_text)

  END SUBROUTINE print_calendar

  !----------------------------------------------------------------------------

  !>
  !! Subroutine to print information on date and time.
  !!
  !! @par Revision History
  !! First version by Marco Giorgetta, MPI-M Hamburg (2010-10-17)
  !!
  SUBROUTINE print_datetime(datetime)

    TYPE (t_datetime), INTENT(in) :: datetime
    !
    ! input components
    ! - year, month, day, hour, minute, second
    ! - time
    !
      
    WRITE(message_text,'(a,i6,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,f9.6,a)') &
      & 'Date/time      : ',datetime%year,'-',datetime%month,'-',datetime%day, &
      & '/',datetime%hour,':',datetime%minute,':',datetime%second,' UT'
    CALL message('mo_datetime/print_datetime', message_text)
    !
    WRITE(message_text,'(a,i9,a,f15.13)') &
      & 'Cal. day+dayfrc: ',datetime%calday,'+',datetime%caltime
    CALL message('mo_datetime/print_datetime', message_text)
    !
    WRITE(message_text,'(a)') ' '
    CALL message('mo_datetime/print_datetime', message_text)

  END SUBROUTINE print_datetime

  !----------------------------------------------------------------------------

  !>
  !! Subroutine to print auxiliary information on date and time
  !!
  !! @par Revision History
  !! First version by Marco Giorgetta, MPI-M Hamburg (2010-10-17)
  !!
  SUBROUTINE print_auxinfos(datetime)

    TYPE (t_datetime), INTENT(in) :: datetime
    !
    ! input components
    ! - year, month, day, hour, minute, second
    ! - time
    !
      
    WRITE(message_text,'(a,i5)')     'Length of year  [dy] : ',datetime%yealen
    CALL message('mo_datetime/print_auxinfos ', message_text)
    !
    WRITE(message_text,'(a,i5)')     'Length of month [dy] : ',datetime%monlen
    CALL message('mo_datetime/print_auxinfos ', message_text)
    !
    WRITE(message_text,'(a,i5)')     'Length of day   [s]  : ',datetime%daylen
    CALL message('mo_datetime/print_auxinfos ', message_text)
    !
    WRITE(message_text,'(a)')        ' '
    CALL message('mo_datetime/print_auxinfos ', message_text)


    WRITE(message_text,'(a,f16.10)') 'Time in year    [dy] : ',datetime%yeatim
    CALL message('mo_datetime/print_auxinfos ', message_text)
    !
    WRITE(message_text,'(a,f16.10)') 'Year fraction   []   : ',datetime%yeafrc
    CALL message('mo_datetime/print_auxinfos ', message_text)
    !
    WRITE(message_text,'(a,i5)')     'Day of year     [dy] : ',datetime%yeaday
    CALL message('mo_datetime/print_auxinfos ', message_text)
    !
    WRITE(message_text,'(a)')        ' '
    CALL message('mo_datetime/print_auxinfos ', message_text)


    WRITE(message_text,'(a,f16.10)') 'Time in month   [dy] : ',datetime%montim
    CALL message('mo_datetime/print_auxinfos ', message_text)
    !
    WRITE(message_text,'(a,f16.10)') 'Month fraction  []   : ',datetime%monfrc
    CALL message('mo_datetime/print_auxinfos ', message_text)
    !
    WRITE(message_text,'(a)')        ' '
    CALL message('mo_datetime/print_auxinfos ', message_text)


    WRITE(message_text,'(a,f16.10)') 'Time in day     [dy] : ',datetime%daytim
    CALL message('mo_datetime/print_auxinfos ', message_text)
    !
    WRITE(message_text,'(a,f16.10)') 'Seconds in day  [s]  : ',datetime%daysec
    CALL message('mo_datetime/print_auxinfos ', message_text)
    !
    WRITE(message_text,'(a)')        ' '
    CALL message('mo_datetime/print_auxinfos ', message_text)

  END SUBROUTINE print_auxinfos

  !---
  !>
  !! Conversion routines between ISO810 data/time string and t_datetime
  !!
  !! iso8601 used to create file names
  !!
  !! Input: year, month, day, hour, minute, second
  !! Output: string in the format of "yyyymmddThhmmssZ"
  !!
  FUNCTION iso8601( datetime )

    TYPE(t_datetime),INTENT(IN) :: datetime
    CHARACTER(len=16) :: iso8601 

    WRITE(iso8601,'(i4.4,2i2.2,a,3i2.2,a)')                            &
         & datetime%year, datetime%month,  datetime%day,          'T', &
         & datetime%hour, datetime%minute, NINT(datetime%second), 'Z'

  END FUNCTION iso8601


  !>
  !! Used to convert back and forth namelist input
  !!
  !! Subroutine to convert a string into a date/time struct
  !!
  !! Input: string in the format of "YEAR-mm-ddThh:mm:ssZ"
  !! Output: datetime struct year, month, day, hour, minute, second
  !!
  SUBROUTINE string_to_datetime ( datetime_str, datetime )

    IMPLICIT NONE

    CHARACTER(len=date_len), INTENT(in) :: datetime_str
    TYPE (t_datetime), INTENT(out)      :: datetime

    CHARACTER(len=32) :: datetime_format
    CHARACTER(len=16) :: year_fmt_str
    INTEGER           :: position_of_first_dash
    INTEGER           :: position_of_second_dash
    INTEGER           :: length
    INTEGER           :: year_format

    INTEGER           :: second

    ! Build the format
    ! -------------------------------------

    length = LEN_TRIM(datetime_str)

    position_of_first_dash = INDEX ( datetime_str, '-', .FALSE. )

    IF ( position_of_first_dash == 1 ) THEN

       ! we have negative years 

       position_of_second_dash = INDEX ( datetime_str(2:length),  '-', .FALSE. )

       year_format = position_of_second_dash - 1

       !TODO: IF ( year_format + 1 > 9 ) Abort!

       WRITE (year_fmt_str, '(A1,I1)' )  'I', year_format + 1

    ELSE

       year_format = position_of_first_dash - 1

       !TODO: IF ( year_format > 9 ) Abort!

       WRITE (year_fmt_str, '(A1,I1)' )  'I', year_format

    ENDIF

   !WRITE (datetime_format,'(A1,A2,A11)') '(',year_fmt_str,',5(X,I2),X)'
    WRITE (datetime_format,'(A1,A2,A13)') '(',year_fmt_str,',5(1X,I2),1X)'

    ! Read from string
    ! --------------------------------------

    READ(datetime_str, TRIM(datetime_format)) datetime%year,   &
         &                              datetime%month,  &
         &                              datetime%day,    &
         &                              datetime%hour,   &
         &                              datetime%minute, &
         &                              second

    datetime%second = REAL(second,wp)

  END SUBROUTINE string_to_datetime

  !>
  !! Subroutine to convert date/time struct into a string
  !!
  !! Input: year, month, day, hour, minute, second
  !! Output: string in the format of "YEAR-mm-ddThh:mm:ssZ"
  !!  where YEAR can have between 1 and 6 digits and can be negative.
  !!
  SUBROUTINE datetime_to_string ( datetime_str, datetime )

    IMPLICIT NONE

    CHARACTER(len=date_len), INTENT(OUT) :: datetime_str
    TYPE (t_datetime), INTENT(IN)        :: datetime

    CHARACTER(len=32) :: datetime_format
    CHARACTER(len=7)  :: year_fmt_str

    INTEGER           :: year_format
    INTEGER           :: second
    INTEGER           :: abs_year

    ! Build the format
    ! -------------------------------------

    ! ... first get the number of digits for the year

    abs_year = ABS(datetime%year)
    year_format = 0

    DO WHILE ( abs_year > 0 )
      abs_year = abs_year / 10
      year_format = year_format + 1
    ENDDO
 
    IF ( datetime%year < 0 ) year_format = year_format + 1

    !TODO: IF ( year_format > 9 ) Abort!

    WRITE (year_fmt_str, '(A1,I1)' )  'I', year_format

    WRITE (datetime_format,'(A1,A2,A15)') '(',year_fmt_str,',5(A1,I2.2),A1)'

    second = INT(datetime%second,i8)

    ! Write string
    ! --------------------------------------

    WRITE(datetime_str, datetime_format ) &
         &           datetime%year,  '-', &
         &           datetime%month, '-', &
         &           datetime%day,   'T', &
         &           datetime%hour,  ':', &
         &           datetime%minute,':', &
         &           second,         'Z'

  END SUBROUTINE datetime_to_string

END MODULE mo_datetime
