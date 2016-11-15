!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
MODULE mo_time_interpolation
  ! mo_time_interpolation [module]
  !   routines for time interpolation of external data sets
  !
  ! Authors;
  !   J.S. Rast, MPI August 2013     base version
  !  
  !----------------------------------------------------------------------
  !
  ! time_weights_limm: calculates interpolation weights and indices for
  !   the interpolation of monthly means (valid at the middle of the month
  !   to any model time step. The suffix _limm is for linear interpolation
  !   of monthly means.
  !
  !----------------------------------------------------------------------
  
  USE mo_datetime,    ONLY: t_datetime, aux_datetime !, print_datetime_all
  USE mo_datetime,    ONLY: OPERATOR(==)
  USE mo_kind,        ONLY: wp
  USE mo_time_interpolation_weights, ONLY: t_wi_limm

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: time_weights_limm 

  CONTAINS

    SUBROUTINE time_weights_limm(event_date,wi)
    !
    TYPE(t_datetime), INTENT(in) :: event_date !any event date
    TYPE(t_wi_limm),  INTENT(out):: wi
    REAL(wp)                     :: zcmonfrc ! current month fraction
    REAL(wp)                     :: zevent_tim !time in month
    REAL(wp)                     :: zcmlen2, znmlen2 !half of current/nearest month length
    TYPE(t_datetime)             :: znevent_date
    
    INTEGER :: seconds_in_month, seconds_in_middle_of_month
    INTEGER :: seconds_in_middle_of_next_month, seconds_in_middle_of_previous_month

    IF (wi%time == event_date) RETURN

    wi%time=event_date !save event_date in wi
    zcmonfrc=event_date%monfrc
    zevent_tim=event_date%montim ! event time since first of month in frac. days
    zcmlen2=event_date%monlen*0.5_wp
    IF (zcmonfrc<=0.5_wp) THEN
      wi%inm1=event_date%month-1  !interpolate between value of previous and current month, 
                               !"nearest" is previous month
      wi%inm2=event_date%month
      znevent_date=event_date
      znevent_date%day=1
      IF (wi%inm1 == 0) THEN
        znevent_date%month=12
        znevent_date%year=event_date%year-1
      ELSE
        znevent_date%month=znevent_date%month-1        
      END IF
      CALL aux_datetime(znevent_date)
      seconds_in_middle_of_previous_month = 43200*znevent_date%monlen
      znmlen2=znevent_date%monlen*0.5_wp
      wi%wgt1=(zcmlen2-zevent_tim)/(zcmlen2+znmlen2)
      wi%wgt2=1._wp-wi%wgt1
    ELSE
      wi%inm1=event_date%month
      wi%inm2=event_date%month+1
      znevent_date=event_date
      znevent_date%day=1
      IF (wi%inm2 == 13) THEN
        znevent_date%month=1
        znevent_date%year=znevent_date%year+1
      ELSE
        znevent_date%month=znevent_date%month+1
      END IF
      CALL aux_datetime(znevent_date)
      seconds_in_middle_of_next_month = 43200*znevent_date%monlen
      znmlen2=znevent_date%monlen*0.5_wp
      wi%wgt2=(zevent_tim-zcmlen2)/(zcmlen2+znmlen2)
      wi%wgt1=1._wp-wi%wgt2
    END IF

    seconds_in_middle_of_month = 43200*event_date%monlen
    seconds_in_month = ((event_date%day-1)*86400+event_date%hour*3600+event_date%minute*60+event_date%second)

    IF (zcmonfrc<=0.5_wp) THEN
      wi%wgt1 = REAL(seconds_in_middle_of_month - seconds_in_month,wp) &
           &   /REAL(seconds_in_middle_of_month + seconds_in_middle_of_previous_month,wp)
      wi%wgt2 = 1.0_wp - wi%wgt1
    ELSE
      wi%wgt2 = REAL(seconds_in_month - seconds_in_middle_of_month,wp) &
           &   /REAL(seconds_in_middle_of_month + seconds_in_middle_of_next_month,wp)
      wi%wgt1 = 1.0_wp - wi%wgt2
    ENDIF
    
!!$  WRITE(0,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
!!$  CALL print_datetime_all(event_date)
!!$  WRITE(0,*) 'wi%inm1,wi%inm2,wi%wgt1,wi%wgt2= ',wi%inm1, wi%inm2, wi%wgt1, wi%wgt2
!!$  WRITE(0,*) '=============================================================='

  END SUBROUTINE time_weights_limm 

END MODULE mo_time_interpolation
