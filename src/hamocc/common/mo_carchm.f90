#include "hamocc_omp_definitions.inc"
MODULE mo_carchm

!> @file mo_carchm.f90
!! @brief Inorganic carbon cycle.
!!        Calc dissolution, update of hydrogen ions
!!
!!
#include "hamocc_omp_definitions.inc"

USE mo_memory_bgc, ONLY     : hi, aksp, akb3, akw3, ak13, ak23, co3, bgctra,&
       &                      aks3,akf3,ak1p3,ak2p3,ak3p3,aksi3,rrrcl, dremcalc

USE mo_kind, ONLY           : wp
USE mo_control_bgc, ONLY    : bgc_nproma, bgc_zlevs

IMPLICIT NONE

PRIVATE

PUBLIC:: calc_dissol, update_hi

CONTAINS

SUBROUTINE calc_dissol ( start_idx, end_idx, klevs, pddpo, psao)

!! Computes calcium carbonate dissolution

  USE mo_param1_bgc, ONLY     : icalc, ialkali, isco212, isilica, iphosph
  
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
          &          ak23(j,k), akw3(j,k),aks3(j,k),akf3(j,k), aksi3(j,k),&
          &          ak1p3(j,k),ak2p3(j,k),ak3p3(j,k),psao(j,k) , akb3(j,k), &
          &          bgctra(j,k,isilica),bgctra(j,k,iphosph),bgctra(j,k,ialkali) )

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




FUNCTION update_hi(hi,c,ak1,ak2,akw,aks,akf,aksi,ak1p,ak2p,ak3p,s,akb,sit,pt,alk) RESULT (h)
 ! update hydrogen ion concentration
 
 REAL(wp) :: ah1, hi
 REAL(wp), INTENT(in):: ak1,ak2,akw,akb,aks,akf,aksi,c,ak1p,ak2p,&
&                       ak3p,sit,pt,alk 
 
 ! LOCAL
 REAL(wp) :: bt, sti,ft, hso4,hf,hsi,hpo4,ab,aw,ac,ah2o,ah2,erel,h,s
 INTEGER:: iter,jit


 ah1=hi
 bt  = rrrcl*s
! sulfate Morris & Riley (1966)
 sti   = 0.14_wp *  s*1.025_wp/1.80655_wp  / 96.062_wp
! fluoride Riley (1965)
 ft    = 0.000067_wp * s*1.025_wp/1.80655_wp / 18.9984_wp

 iter  = 0

 DO jit = 1,20

     hso4 = sti / ( 1._wp + aks / ( ah1 / ( 1._wp + sti / aks ) ) )
     hf   = 1._wp / ( 1._wp + akf / ah1 )
     hsi  = 1._wp/ ( 1._wp + ah1 / aksi )
     hpo4 = ( ak1p * ak2p * ( ah1 + 2._wp * ak3p ) - ah1**3 ) /    & 
           &        ( ah1**3 + ak1p * ah1**2 + ak1p * ak2p * ah1 + ak1p * ak2p*ak3p )
     ab   = bt / ( 1._wp + ah1 / akb )
     aw   = akw / ah1 - ah1 / ( 1._wp + sti / aks )
     ac   = alk + hso4 - sit * hsi - ab - aw + ft * hf - pt * hpo4
     ah2o = SQRT( ( c - ac )**2 + 4._wp * ( ac * ak2 / ak1 ) * ( 2._wp*c - ac ) )
     ah2  = 0.5_wp * ak1 / ac *( ( c - ac ) + ah2o )
     erel = ( ah2 - ah1 ) / ah2

     if (abs( erel ).ge.5.e-5_wp) then
         ah1 = ah2
         iter = iter + 1
       else
         ah1 = ah2
         exit
      endif
 ENDDO
 if(ah1.gt.0._wp)h=max(1.e-20_wp,ah1)

END FUNCTION

END MODULE
