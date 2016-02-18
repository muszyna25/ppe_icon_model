!>
!! @file carchm.f90
!! @brief Inorganic carbon cycle.
!!
!! Computes carbonate chemisty in the water column, and sea-air gass exchange
!! for oxygen, DMS, O2, N2, and CO2. It computes dissolution of calcium carbonate.
!!
!!
!!

SUBROUTINE carchm ( start_idx, end_idx, klevs, pddpo, psao)

  USE mo_biomod, ONLY         : rrrcl, dremcalc
  USE mo_control_bgc, ONLY    : dtbgc, bgc_nproma, bgc_zlevs

  
  USE mo_param1_bgc, ONLY     : icalc, igasnit, ian2o, iatmco2,          &
       &                        iatmo2, iatmn2, ioxygen, ialkali,        &
       &                        isco212, kphosy

  USE mo_carbch, ONLY         : hi, aksp, akb3, akw3, ak13, ak23, co3, bgctra, bgctend
  USE mo_kind, ONLY           : wp
  
  IMPLICIT NONE

  !! Arguments

  INTEGER, INTENT(in) :: start_idx                
  INTEGER, INTENT(in) :: end_idx                
  INTEGER, INTENT(in) :: klevs(bgc_nproma)                  !< 3rd (vertical) REAL of model grid.

  REAL(wp),INTENT(in) :: pddpo(bgc_nproma,bgc_zlevs) !< size of scalar grid cell (3rd REAL) [m]
  REAL(wp),INTENT(in) :: psao(bgc_nproma,bgc_zlevs)  !< salinity

  !! Local variables

  INTEGER :: k, iter, j, kpke 

  INTEGER :: laumo1


  REAL(wp) :: supsat, undsa, dissol
  REAL(wp) :: dddhhh,dadh,a,h,c,alk,t1,t2
  REAL(wp) :: akbi,ak2i,ak1i

  REAL(wp) :: AK0,AK1,AK2,AKB,AKW,BT
  REAL(wp) :: supsatup,satdiff,depthdiff   ! needed to calculate depth of lysokline
  !
   !*********************************************************************
  !
  ! Dissolution of calcite, whole water column
  !
  !*********************************************************************


  !*        22. CHEMICAL CONSTANTS - water column
  !     -----------------------------------------------------------------

  DO iter = 1, 3
    DO j=start_idx,end_idx
     kpke=klevs(j)
     DO k = 1, kpke
              IF (pddpo(j,k) > 0.5_wp) THEN
                 h   = hi(j,k)
                 c   = bgctra(j,k,isco212)
                 t1  = h/ak13(j,k)
                 t2  = h/ak23(j,k)
                 ak2 = ak23(j,k)
                 akw = akw3(j,k)
                 bt  = rrrcl*psao(j,k)
                 akb = akb3(j,k)
                 alk = bgctra(j,k,ialkali)
                 ! Determine hydrogen ion HI so that ALK(DIC,BT,HI) matches given alk by Newton iteration
                 ! Actual mismatch
                 a = c * (2._wp + t2) / (1._wp+t2+t2*t1)+akw/h-h+bt/(1._wp+h/akb)-alk
                 ! Derivative
                 dadh = c * (1._wp / (ak2 * (1._wp + t2 + t2 * t1)) &
                      - (2._wp + t2) * (1._wp / ak2 + 2._wp * t1 / ak2)/ &
                      &          (1._wp + t2 + t2 * t1)**2)                          &
                      &          - akw / h**2 - 1._wp - (bt / akb) / (1._wp + h / akb)**2
                 dddhhh=a/dadh
                 h = MAX(h - dddhhh, 1.e-10_wp) ! Prevent overshooting to negative values at start of iteration
                 hi(j,k) = h
                 co3(j,k) = c/(1._wp+h*(1._wp+h/ak13(j,k))/ak23(j,k))
              ENDIF ! wet cell
        END DO
     END DO
  END DO

  DO j= start_idx, end_idx
       kpke=klevs(j)
        IF(pddpo(j,1) > 0.5_wp) THEN
          ! surface
          ! dissolution in surface layer
              supsat = co3(j,1)-97._wp*aksp(j,1)
              undsa  = MAX(0._wp, -supsat)
              dissol = MIN(undsa,dremcalc*bgctra(j,1,icalc))
              bgctra(j,1,icalc)   = bgctra(j,1,icalc)-dissol
              bgctra(j,1,ialkali) = bgctra(j,1,ialkali)+2._wp*dissol
              bgctra(j,1,isco212) = bgctra(j,1,isco212)+dissol
        ENDIF
        DO k = 2, kpke
           IF(pddpo(j,k) > 0.5_wp) THEN
   
           ! water column
           ! dissolution in subsurface layer and deeper layers 
              supsat = co3(j,k)-97._wp*aksp(j,k)   ! 97. = 1./1.03e-2 (MEAN TOTAL [CA++] IN SEAWATER [kmol/m3])
              undsa  = MAX(0._wp, -supsat)
             
              dissol = MIN(undsa,dremcalc*bgctra(j,k,icalc))
              bgctra(j,k,icalc)   = bgctra(j,k,icalc)-dissol
              !   bgctend(j,k,kphosy) = hi(j,k)
              bgctra(j,k,ialkali) = bgctra(j,k,ialkali)+2._wp*dissol

              bgctra(j,k,isco212) = bgctra(j,k,isco212)+dissol

            ENDIF   ! wet cell
         END DO
  END DO



END SUBROUTINE
