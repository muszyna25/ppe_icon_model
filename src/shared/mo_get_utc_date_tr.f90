!>
!! Computation of actual date of current timestep (SUBROUTINE get_utc_date_tr).
!!
!! @author Thorsten Reinhardt, AGeoBw, Offenbach (2010-10-18)
!!
!!
!! @par Revision History
!! Initial Release by Thorsten Reinhardt, AGeoBw, Offenbach (2010-10-18)
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
MODULE mo_get_utc_date_tr

  USE mo_kind,                 ONLY: wp, i8
  USE mo_time_nml,             ONLY: ini_datetime
  
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: get_utc_date_tr

CONTAINS

  !>
  !! This routine determines the actual date of this forecast step.
  !! Using the date of the forecast-start and the given forecast time,
  !! the actual date is calculated taking leap-years into consideration.
  !! The date is given in three different formats.
  !! Based on subroutine get_utc_date of the COSMO model version 4.6 (utilities.f90)
  !! (Modification: Use of simulation time p_sim_time instead of timestep).
  !!
  !! @par Revision History
  !! Initial release by Thorsten Reinhardt, AGeoBw (2010-10-18)
  !! Use initial date from mo_global_variables' ini_datetime (T.R., 2010-12-16)
  !! yactdate2 now optional (T.R., 2011-02-09)
  !!
  SUBROUTINE get_utc_date_tr (p_sim_time, itype_calendar, iyear, nactday, acthour)

    INTEGER   , INTENT(IN)   ::                           &
      itype_calendar    ! for specifying the calendar used

    REAL      (wp), INTENT(IN)      ::                           &
      p_sim_time         ! simulation time in seconds

!!$    CHARACTER (LEN=10), INTENT(OUT)                    ::        &
!!$      yactdate1  ! actual date in the form   yyyymmddhh
!!$
!!$    CHARACTER (LEN=22), INTENT(OUT), OPTIONAL          ::        &
!!$      yactdate2  ! actual date in the form   wd   dd.mm.yy  hh UTC


    INTEGER   , INTENT(OUT)  ::                           &
      iyear    ! year
    
    INTEGER   , INTENT(OUT)  ::                           &
      nactday    ! day of the year

    REAL      (wp), INTENT(OUT)     ::                           &
      acthour    ! actual hour of the day

    ! Local variables:
    INTEGER     ::                                       &
      month(12), monthsum(13), ileap, iy, m,                       &
      idd, imm, iyy, ihh, iday, imonth, ihour, immhours, iyyhours, &
      iyear_hours

    CHARACTER (LEN=3)            :: yweek(7)

    ! And for computing the amount of seconds of the whole forecast time,
    ! an 8-Byte INTEGER has to be used. Otherwise the computation fails after
    ! approx. 68 years!!

    REAL(wp)  :: zseconds

    !------------ End of header ---------------------------------------------------

    ! Begin subroutine get_utc_date

    DATA         month  / 31 ,  28 ,  31 ,  30 ,  31 ,  30 ,       &
      31 ,  31 ,  30 ,  31 ,  30 ,  31 /
    DATA         yweek  /'MON', 'TUE', 'WED', 'THU', 'FRI', 'SAT', 'SUN' /


    ! Statementfunction: ileap(yy) = 0:  no leap year,
    !                    ileap(yy) = 1:  leap year
    ileap (iy) = IABS( MOD(iy,4) - 4) / 4

    ! year, month, day, hours of initial date
    iyy = ini_datetime%year
    imm = ini_datetime%month
    idd = ini_datetime%day 
    ihh = ini_datetime%hour

    IF     (itype_calendar == 0) THEN
      month (2)    = 28 + ileap (iyy)
      monthsum(1) =  0
      DO m =  2 , 13
        monthsum(m) =  monthsum(m-1) + month(m-1)
      enddo
    ELSEIF (itype_calendar == 1) THEN
      monthsum(1) =  0
      DO m =  2 , 13
        monthsum(m) =  monthsum(m-1) + 30
      enddo
    ENDIF

    ! Determine how many hours have passed in this year
    iyyhours = (idd*24) + monthsum(imm)*24 + (ihh-24)
    iyyhours = iyyhours +                                           &
      INT (NINT (p_sim_time, i8)/3600_i8)

    ! Take turning of the year into account
    IF     (itype_calendar == 0) THEN
      iyear_hours = 8760 + ileap(iyy)*24
    ELSEIF (itype_calendar == 1) THEN
      iyear_hours = 8640
    ENDIF

    IF (iyyhours < 0) THEN
      iyear    = iyy-1
      IF     (itype_calendar == 0) THEN
        iyyhours = 8760 + ileap(iyear)*24 + iyyhours
      ELSEIF (itype_calendar == 1) THEN
        iyyhours = 8640                            + iyyhours
      ENDIF
    ELSE IF (iyyhours >= iyear_hours) THEN
      ! Take also into account if the run lasts for several years
      iyear    = iyy
      IF     (itype_calendar == 0) THEN
        iyear_hours = 8760 + ileap(iyear)*24
      ELSEIF (itype_calendar == 1) THEN
        iyear_hours = 8640
      ENDIF

      DO WHILE (iyyhours >= iyear_hours)
        iyyhours = iyyhours - iyear_hours
        iyear=iyear+1
        IF     (itype_calendar == 0) THEN
          iyear_hours = 8760 + ileap(iyear)*24
        ELSEIF (itype_calendar == 1) THEN
          iyear_hours = 8640
        ENDIF
      ENDDO
    ELSE
      iyear    =   iyy
    ENDIF

    ! calculate monthsum for actual year
    IF     (itype_calendar == 0) THEN
      month (2)    = 28 + ileap (iyear)
      monthsum(1) =  0
      DO m =  2 , 13
        monthsum(m) =  monthsum(m-1) + month(m-1)
      enddo
    ELSEIF (itype_calendar == 1) THEN
      monthsum(1) =  0
      DO m =  2 , 13
        monthsum(m) =  monthsum(m-1) + 30
      enddo
    ENDIF

    ! Determine the actual date from iyyhours
    m        = 1
    immhours = iyyhours
    DO WHILE (immhours >= 0)
      m        = m+1
      immhours = iyyhours - monthsum(m) * 24
    ENDDO
    imonth   = m-1

    immhours = iyyhours - monthsum(imonth)*24
    iday     = immhours/24 + 1
    ihour    = MOD(immhours,24)

    zseconds = p_sim_time / 3600.0_wp
    acthour  = REAL (ihour, wp) +                          &
      (zseconds - REAL(INT(zseconds, i8), wp))

!!$    ihour    = INT(acthour)
    nactday  = monthsum(imonth) + iday + INT(acthour/24._wp + 0.0001_wp)
!!$    iweek    = MOD(monthsum(imonth) + iday + (iyear-1901) + (iyear-1901)/4, 7)+1
!!$
!!$    WRITE ( yactdate1(1:4) , '(I4.4)' ) iyear
!!$    WRITE ( yactdate1(5:6) , '(I2.2)' ) imonth
!!$    WRITE ( yactdate1(7:8) , '(I2.2)' ) iday
!!$    WRITE ( yactdate1(9:10), '(I2.2)' ) ihour
!!$
!!$    IF ( PRESENT(yactdate2) ) THEN
!!$      IF     (itype_calendar == 0) THEN
!!$        yactdate2 = yweek(iweek)//' '//yactdate1(7:8)//'.'// yactdate1(5:6)//'.' &
!!$          //yactdate1(1:4)//'  '//yactdate1(9:10)//' UTC'
!!$      ELSEIF (itype_calendar == 1) THEN
!!$        yactdate2 = '    '//yactdate1(7:8)//'.'// yactdate1(5:6)//'.' &
!!$          //yactdate1(1:4)//'  '//yactdate1(9:10)//' UTC'
!!$      ENDIF
!!$    ENDIF

  END SUBROUTINE get_utc_date_tr

END MODULE mo_get_utc_date_tr
