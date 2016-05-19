!>
!! @file powadi.f90
!! @brief vertical diffusion with simultaneous dissolution,
!! implicit discretisation.
!!
!!
SUBROUTINE powadi ( j,  solrat, sedb1, sediso, bolven)

  USE mo_kind, ONLY       : wp

  USE mo_sedmnt, ONLY     : ks, sedict, seddzi, seddw, porwah, porwat

  USE mo_biomod, ONLY     : bolay


  IMPLICIT NONE

  !! Arguments

  INTEGER, INTENT(in)     :: j             !< zonal grid index

  REAL(wp), INTENT(in)    :: solrat(ks)    !< dissolution rate
  REAL(wp), INTENT(inout) :: sedb1(0:ks)   !< tracer at entry
  REAL(wp), INTENT(inout) :: sediso(0:ks)  !<
  REAL(wp), INTENT(in)    :: bolven        !< bottom layer ventilation rate

  !! Local variables

  INTEGER,  SAVE  :: k,l   
  REAL(wp), SAVE :: tredsy(0:ks,3)
  REAL(wp), SAVE :: asu,alo
  !
  !----------------------------------------------------------------------
  !
  DO k = 1, ks

     asu=sedict*seddzi(k)*porwah(k)
     alo = 0._wp

     IF (k < ks) alo = sedict*seddzi(k+1)*porwah(k+1)

        tredsy(k,1) = -asu
        tredsy(k,3) = -alo
        tredsy(k,2) = seddw(k)*porwat(k) - tredsy(k,1)       &
             &        - tredsy(k,3) + solrat(k)*porwat(k)*seddw(k)

  END DO

  k = 0

  asu = 0._wp
  alo = sedict*seddzi(1)*porwah(1)

     IF(bolay(j) > 0._wp) THEN
        tredsy(k,1) = -asu
        tredsy(k,3) = -alo
        tredsy(k,2) = bolven*bolay(j) - tredsy(k,1) - tredsy(k,3)
     ELSE
        tredsy(k,1) = 0._wp
        tredsy(k,3) = 0._wp
        tredsy(k,2) = 0._wp
     ENDIF


  DO k = 1, ks
        IF (bolay(j) > 0._wp) THEN
           tredsy(k-1,1) = tredsy(k,1) / tredsy(k-1,2)
           tredsy(k,2)   = tredsy(k,2)                       &
                &          - tredsy(k-1,3) * tredsy(k,1) / tredsy(k-1,2)
        ENDIF
  END DO

  DO k = 1, ks
        sedb1(k) = sedb1(k) - tredsy(k-1,1) * sedb1(k-1)
  END DO

  k = ks

     IF (bolay(j) > 0._wp)sediso(k) = sedb1(k) / tredsy(k,2)

  DO k = 1, ks
     l = ks-k
        IF (bolay(j) > 0._wp) sediso(l) =                           &
             &           ( sedb1(l) - tredsy(l,3) * sediso(l+1) )       &
             &           / tredsy(l,2)
  END DO

END SUBROUTINE 
