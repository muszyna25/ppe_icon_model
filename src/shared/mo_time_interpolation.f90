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
  
  USE mo_datetime,    ONLY: t_datetime, aux_datetime, print_datetime_all
  USE mo_kind,        ONLY: wp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: wgt1_limm, wgt2_limm, inm1_limm, inm2_limm, time_weights_limm 

  REAL(wp)                     :: wgt1_limm, wgt2_limm !interpolation weights for months with
  INTEGER                      :: inm1_limm, inm2_limm !inm1_limm and inm2_limm 
                                                       !as indices [0,..,13] respectively


  CONTAINS

    SUBROUTINE time_weights_limm(event_date)
    !
    TYPE(t_datetime), INTENT(in) :: event_date !any event date
    
    REAL(wp)                     :: zcmonfrc ! current month fraction
    REAL(wp)                     :: zevent_tim !time in month
    REAL(wp)                     :: zcmlen2, znmlen2 !half of current/nearest month length
    TYPE(t_datetime)             :: znevent_date

    zcmonfrc=event_date%monfrc
    zevent_tim=event_date%montim ! event time since first of month in frac. days
    zcmlen2=event_date%monlen*0.5_wp
    IF (zcmonfrc<=0.5_wp) THEN
      inm1_limm=event_date%month-1  !interpolate between value of previous and current month, 
                               !"nearest" is previous month
      inm2_limm=event_date%month
      znevent_date=event_date
      znevent_date%day=1
      IF (inm1_limm == 0) THEN
        znevent_date%month=12
        znevent_date%year=event_date%year-1
      END IF
      CALL aux_datetime(znevent_date)
      znmlen2=znevent_date%monlen*0.5_wp
      wgt1_limm=(zcmlen2-zevent_tim)/(zcmlen2+znmlen2)
      wgt2_limm=1._wp-wgt1_limm
    ELSE
      inm1_limm=event_date%month
      inm2_limm=event_date%month+1
      znevent_date=event_date
      znevent_date%day=1
      IF (inm2_limm == 13) THEN
        znevent_date%month=1
        znevent_date%year=znevent_date%year+1
      END IF
      CALL aux_datetime(znevent_date)
      znmlen2=znevent_date%monlen*0.5_wp
      wgt2_limm=(zevent_tim-zcmlen2)/(zcmlen2+znmlen2)
      wgt1_limm=1._wp-wgt2_limm
    END IF
!!$  WRITE(0,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
!!$  CALL print_datetime_all(event_date)
!!$  WRITE(0,*) 'inm1_limm,inm2_limm,wgt1_limm,wgt2_limm= ',inm1_limm, inm2_limm, wgt1_limm, wgt2_limm
!!$  WRITE(0,*) '=============================================================='
  END SUBROUTINE time_weights_limm 

END MODULE mo_time_interpolation
