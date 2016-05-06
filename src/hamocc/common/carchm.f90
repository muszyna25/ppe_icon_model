!>
!! @file carchm.f90
!! @brief Inorganic carbon cycle.
!!
!! Computes dissolution of calcium carbonate.
!!
!!
!!
#include "hamocc_omp_definitions.inc"

SUBROUTINE carchm ( start_idx, end_idx, klevs, pddpo, psao)

  USE mo_biomod, ONLY         : rrrcl, dremcalc
  USE mo_control_bgc, ONLY    : bgc_nproma, bgc_zlevs

  
  USE mo_param1_bgc, ONLY     : icalc,            &
       &                        ialkali,        &
       &                        isco212

  USE mo_carbch, ONLY         : hi, aksp, akb3, akw3, ak13, ak23, co3, bgctra
  USE mo_kind, ONLY           : wp
  
  IMPLICIT NONE

  !! Arguments

  INTEGER, INTENT(in) :: start_idx                   !< start index for j loop (ICON cells, MPIOM lat dir)        
  INTEGER, INTENT(in) :: end_idx                     !< end index  for j loop  (ICON cells, MPIOM lat dir) 
  INTEGER, INTENT(in) :: klevs(bgc_nproma)           !< vertical levels

  REAL(wp),INTENT(in) :: pddpo(bgc_nproma,bgc_zlevs) !< size of scalar grid cell (3rd REAL) [m]
  REAL(wp),INTENT(in) :: psao(bgc_nproma,bgc_zlevs)  !< salinity

  !! Local variables

  INTEGER :: k, iter, j, kpke 

  REAL(wp) :: supsat, undsa, dissol
  REAL(wp) :: dddhhh,dadh,a,h,c,alk,t1,t2

  REAL(wp) :: AK2,AKB,AKW,BT
  !
   !*********************************************************************
  !
  ! Dissolution of calcite, whole water column
  !
  !*********************************************************************


  !*        22. CHEMICAL CONSTANTS - water column
  !     -----------------------------------------------------------------

!HAMOCC_OMP_PARALLEL
!HAMOCC_OMP_DO PRIVATE(iter,j,kpke,k,h,&
!HAMOCC_OMP           c,t1,t2,ak2,akw,bt,akb,alk,a,dadh,dddhhh) HAMOCC_OMP_DEFAULT_SCHEDULE
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
!HAMOCC_OMP_END_DO
!HAMOCC_OMP_END_PARALLEL

!HAMOCC_OMP_PARALLEL
!HAMOCC_OMP_DO PRIVATE(j,kpke,k,supsat,&
!HAMOCC_OMP            undsa,dissol) HAMOCC_OMP_DEFAULT_SCHEDULE 
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
!HAMOCC_OMP_END_DO
!HAMOCC_OMP_END_PARALLEL



END SUBROUTINE
