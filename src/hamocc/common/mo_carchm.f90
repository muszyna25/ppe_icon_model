#include "hamocc_omp_definitions.inc"
MODULE mo_carchm

!> @file mo_carchm.f90
!! @brief Inorganic carbon cycle.
!!        Calc dissolution, update of hydrogen ions
!!
!!
#include "hamocc_omp_definitions.inc"

USE mo_carbch, ONLY         : hi, aksp, akb3, akw3, ak13, ak23, co3, bgctra, bgctend
USE mo_kind, ONLY           : wp
USE mo_biomod, ONLY         : rrrcl, dremcalc
USE mo_control_bgc, ONLY    : dtbgc, bgc_nproma, bgc_zlevs

IMPLICIT NONE

PRIVATE

PUBLIC:: calc_dissol, update_hi

CONTAINS

SUBROUTINE calc_dissol ( start_idx, end_idx, klevs, pddpo, psao)

!! Computes calcium carbonate dissolution

  USE mo_param1_bgc, ONLY     : icalc, ialkali, isco212
  
  IMPLICIT NONE

  !! Arguments

  INTEGER, INTENT(in) :: start_idx             !< start index for j loop (ICON cells, MPIOM lat dir)    
  INTEGER, INTENT(in) :: end_idx               !< end index  for j loop  (ICON cells, MPIOM lat dir)    
  INTEGER, INTENT(in) :: klevs(bgc_nproma)     !<  vertical levels
           

  REAL(wp),INTENT(in) :: pddpo(bgc_nproma,bgc_zlevs) !< size of scalar grid cell (3rd REAL) [m]
  REAL(wp),INTENT(in) :: psao(bgc_nproma,bgc_zlevs)  !< salinity

  !! Local variables

  INTEGER :: k, j, kpke 

  REAL(wp) :: supsat, undsa, dissol
  !
   !*********************************************************************
  !
  ! Dissolution of calcite, whole water column
  !
  !*********************************************************************
!HAMOCC_OMP_PARALLEL
!HAMOCC_OMP_DO PRIVATE(k,supsat,undsa, dissol) HAMOCC_OMP_DEFAULT_SCHEDULE
  DO j= start_idx, end_idx

       kpke=klevs(j)
    
        DO k = 1, kpke

           IF(pddpo(j,k) > 0.5_wp) THEN

       
              hi(j,k) = update_hi(hi(j,k), bgctra(j,k,isco212), ak13(j,k) , &
          &          ak23(j,k), akw3(j,k),psao(j,k) , akb3(j,k), bgctra(j,k,ialkali) )

              co3(j,k) = bgctra(j,k,isco212)/(1._wp+hi(j,k)*(1._wp+hi(j,k)/ak13(j,k))/ak23(j,k))

              supsat = co3(j,k)-97._wp*aksp(j,k)   ! 97. = 1./1.03e-2 (MEAN TOTAL [CA++] IN SEAWATER [kmol/m3])
              undsa  = MAX(0._wp, -supsat)
             
              dissol = MIN(undsa,dremcalc*bgctra(j,k,icalc))
              bgctra(j,k,icalc)   = bgctra(j,k,icalc)-dissol
              bgctra(j,k,ialkali) = bgctra(j,k,ialkali)+2._wp*dissol

              bgctra(j,k,isco212) = bgctra(j,k,isco212)+dissol

            ENDIF   ! wet cell

         END DO
  END DO
!HAMOCC_OMP_END_DO
!HAMOCC_OMP_END_PARALLEL
END SUBROUTINE

FUNCTION update_hi(hi,c,ak1,ak2,akw,s,akb,alk) RESULT (h)
  
 REAL(wp) :: h
 REAL(wp), INTENT(in):: hi,c,ak1,ak2,akw,s,akb,alk
 
 ! LOCAL
 REAL(wp) :: bt, a, t2,t1, dadh, dddhhh
 INTEGER:: iter, niter

 niter=3

 h=hi
 bt  = rrrcl*s
 DO iter=1,niter

    t1  = h/ak1
    t2  = h/ak2
    ! Determine hydrogen ion HI so that ALK(DIC,BT,HI) matches given alk by Newton iteration
    ! Actual mismatch
    a = c * (2._wp + t2) / (1._wp+t2+t2*t1)+akw/h-h+bt/(1._wp+h/akb)-alk
    ! Derivative
    dadh = c * (1._wp / (ak2 * (1._wp + t2 + t2 * t1)) &
  &   - (2._wp + t2) * (1._wp / ak2 + 2._wp * t1 / ak2)/ &
  &          (1._wp + t2 + t2 * t1)**2)                          &
  &          - akw / h**2 - 1._wp - (bt / akb) / (1._wp + h / akb)**2
    dddhhh=a/dadh
    h = MAX(h - dddhhh, 1.e-10_wp) ! Prevent overshooting to negative values at start of iteration
 ENDDO

END FUNCTION

END MODULE
