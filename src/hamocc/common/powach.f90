
!>!! @file powach.f90
!!   @brief compute sediment chemistry, explicit discretisation
!!   call pore water diffusion
!!
#include "hamocc_omp_definitions.inc"

   SUBROUTINE POWACH(start_idx,end_idx,psao,pddpo)

! 
   USE mo_kind, ONLY        : wp

   USE mo_carbch, ONLY      : ak13, ak23, akb3, akw3, aksp,               &
        &                     aks3,akf3,ak1p3,ak2p3,ak3p3,aksi3          
   USE mo_sedmnt, ONLY      : sedlay, sedhpl, seddw, silpro,              &
        &                     powtra, prcaca, prorca, produs,             &
        &                     pown2bud, powh2obud, sred_sed,              &
        &                     porsol, pors2w, calcon,                     &
        &                     disso_op, disso_cal,  ks, o2thresh,         &
        &                     sedtend, isremins, isremino, isreminn,      &
        &                     silsat

   USE mo_biomod, ONLY      : bolay, kbo, ro2ut, rnit, nitdem, n2prod,    &
        &                     rrrcl, rcar
 
   USE mo_control_bgc, ONLY : dtbgc, bgc_nproma, bgc_zlevs
 
   USE mo_param1_bgc, ONLY  : ipowasi, issssil, ipowaox,     &
        &                     issso12, ipowaph, ipowno3, ipown2, &
        &                     ipowaal, ipowaic, isssc12, issster
 
   USE mo_hamocc_nml, ONLY  : disso_po
 
  IMPLICIT NONE

  !! Arguments

   INTEGER, INTENT(in)  :: start_idx     !< start index for j loop (ICON cells, MPIOM lat dir)          
   INTEGER, INTENT(in)  :: end_idx        !< end index  for j loop  (ICON cells, MPIOM lat dir)        

   REAL(wp), INTENT(in) :: psao(bgc_nproma,bgc_zlevs)   !< salinity [psu.].
   REAL(wp), INTENT(in) :: pddpo(bgc_nproma,bgc_zlevs)  !< size of scalar grid cell 
  
  !! Local variables

   INTEGER :: j,k, iter


   REAL(wp) :: orgsed                        ! sum of C12 and C13 in organic sed.
   REAL(wp) :: sssnew                        ! temporarily value of solid constituent after solution

   REAL(wp) :: powcar(ks)

   REAL(wp) :: o2lim          ! o2-limitation of processes in suboxic water
   REAL(wp) :: undsa, posol, pomax
   REAL(wp) :: umfa,denit,bt,alk,c
   REAL(wp) :: ak1,ak2,akb,akw
   REAL(wp) :: h,a,satlev
   REAL(wp) :: AKSI,AKF,AKS,AK1P,AK2P,AK3P
   REAL(wp) :: sti,ft,sit,pt,ah1       ! Chr. Heinze [H+]
   REAL(wp) :: hso4,hf,hsi,hpo4
   REAL(wp) :: ab,aw,ac,ah2o,ah2,erel,pco2_ch
   INTEGER  :: jit



!!! WE start with remineralisation of organic to estimate alkalinity changes first
!          
!HAMOCC_OMP_PARALLEL
!HAMOCC_OMP_DO PRIVATE(j,k,o2lim,pomax,posol,orgsed,sssnew,&
!HAMOCC_OMP           undsa,iter,bt,alk,c,ak1,ak2,akb,akw,h,&
!HAMOCC_OMP           a,satlev,powcar,sti,ft,sit,&
!HAMOCC_OMP           pt,ah1,aks,akf,ak1p,ak2p,ak3p,aksi,iter, jit,&
!HAMOCC_OMP           hso4,hf,hpo4,ab,aw,ac,ah2o,ah2,erel,ah2) HAMOCC_OMP_DEFAULT_SCHEDULE


    Do j = start_idx, end_idx

         IF(bolay(j) > 0._wp) THEN
            sedlay(j,1,issso12)                                     &
     &      = sedlay(j,1,issso12) + prorca(j)/(porsol(1)*seddw(1))
         !   prorca(j) = 0._wp
         ENDIF


! CALCULATE OXYGEN-POC CYCLE 
!*************************************************************
! This scheme is not based on undersaturation, but on O2 itself
! Calculate new solid sediment.
! Update pore water concentration.

      DO k=1,ks
         IF(bolay(j) > 0._wp) THEN
         IF (powtra(j, k, ipowaox) > 2.e-6_wp) THEN

         o2lim = powtra(j,k,ipowaox)/(o2thresh+powtra(j,k,ipowaox)) ! o2 limitation in oxic water

        ! maximal possible aerobe dissolution per time step 
        ! limited by available oxygen concentration, currently max 70 % of O2 
         pomax = disso_po * max(0._wp,sedlay(j,k,issso12))*powtra(j,k,ipowaox)
         posol = min(0.7_wp*powtra(j,k,ipowaox)/ro2ut,pomax*pors2w(k))

         
         sedlay(j,k,issso12) = sedlay(j,k,issso12) - posol

         powtra(j,k,ipowaic) = powtra(j,k,ipowaic) + posol*rcar*pors2w(k)
         powtra(j,k,ipowaph)=powtra(j,k,ipowaph)+posol*pors2w(k)
	 powtra(j,k,ipowno3)=powtra(j,k,ipowno3)+posol*rnit*pors2w(k)
         powtra(j,k,ipowaal)=powtra(j,k,ipowaal)-posol*rnit*pors2w(k)
         powtra(j,k,ipowaox)=powtra(j,k,ipowaox)-posol*pors2w(k)*ro2ut

         sedtend(j,k,isremino) = posol*pors2w(k)/dtbgc

         else
          sedtend(j,k,isremino) = 0._wp


         ENDIF ! oxygen > 2.
         ENDIF ! bolay > 0

      ENDDO


! CALCULATE NITRATE REDUCTION UNDER ANAEROBIC CONDITIONS EXPLICITELY
!*******************************************************************
! Denitrification rate constant of POP (disso) [1/sec]
! Store flux in array anaerob, for later computation of DIC and alkalinity.


      DO  k=1,ks
         IF (bolay( j) .GT. 0._wp) THEN
         IF (powtra( j, k, ipowaox) < 2.e-6_wp) THEN

           orgsed = max(0._wp,sedlay(j,k,issso12)) 

           posol = denit * MIN(0.5_wp * powtra(j, k, ipowno3)/nitdem, orgsed)
           sedlay(j,k,issso12)=sedlay(j,k,issso12)-posol
           powtra(j,k,ipowaph)=powtra(j,k,ipowaph)+posol*pors2w(k)
           powtra(j,k,ipowaic)=powtra(j,k,ipowaic)+rcar*posol*pors2w(k)

           powtra(j,k,ipowaal)=powtra(j,k,ipowaal)+posol*nitdem*pors2w(k)
           powtra(j,k,ipowno3)=powtra(j,k,ipowno3)-nitdem*posol*pors2w(k)
           powtra(j,k,ipown2)=powtra(j,k,ipown2)+n2prod*posol*pors2w(k)
           powh2obud(j,k)=powh2obud(j,k)+0.5_wp*n2prod*posol*pors2w(k)
           pown2bud(j,k) = pown2bud(j,k) + 2._wp*n2prod*posol*pors2w(k)

           sedtend(j,k,isreminn) = posol*pors2w(k)/dtbgc
         else
           sedtend(j,k,isreminn) = 0._wp

         ENDIF   ! oxygen <1.e-6
         endif   ! bolay
        ENDDO



!  sulphate reduction in sediments
! Denitrification rate constant of POP (disso) [1/sec]

      DO  k=1,ks
         IF (bolay( j) > 0._wp) THEN
         IF (powtra( j, k, ipowaox) < 1.e-6_wp) THEN
           orgsed=max(0._wp,sedlay(j,k,issso12))
         ! reduced by factor 100
           sssnew = orgsed/(1._wp + sred_sed)
           posol=orgsed - sssnew   ! change for org sed.

           sedlay(j,k,issso12)=sedlay(j,k,issso12)-posol
           powtra(j,k,ipowaic)=powtra(j,k,ipowaic)+posol*pors2w(k)*rcar
           powtra(j,k,ipowaph)=powtra(j,k,ipowaph)+posol*pors2w(k)
           powtra(j,k,ipowno3)=powtra(j,k,ipowno3)+posol*rnit*pors2w(k)
           powtra(j,k,ipowaal)=powtra(j,k,ipowaal)+posol*rnit*pors2w(k) ! alkalinity correction 2015
           powh2obud(j,k)=powh2obud(j,k)-posol*ro2ut*pors2w(k)
           pown2bud(j,k) = pown2bud(j,k) + 2._wp*rnit*posol*pors2w(k) 

           sedtend(j,k,isremins) = posol*pors2w(k)/dtbgc
         else
           sedtend(j,k,isremins) = 0._wp
         ENDIF   ! oxygen <1.e-6

         endif   ! bolay
       ENDDO
!!! End remineralization POC




! CALCULATE SILICATE-OPAL CYCLE
!******************************************************************
! Calculate updated degradation rate from updated undersaturation.
! Calculate new solid sediment.
! Update pore water concentration from new undersaturation.

! Add sediment flux of opal
         IF(bolay(j).GT.0._wp) THEN

             sedlay(j,1,issssil)=                                     &
      &      sedlay(j,1,issssil)+silpro(j)/(porsol(1)*seddw(1))
            ! silpro(j)=0._wp
          ENDIF


      DO  k=1,ks
         IF(bolay(j).GT.0._wp) THEN
            undsa=MAX(silsat-powtra(j,k,ipowasi),0._wp)
! explixit version     
! new implicit within layer
             sssnew = sedlay(j,k,issssil)/(1._wp+ disso_op*undsa)
             posol =  sedlay(j,k,issssil) - sssnew

             sedlay(j,k,issssil)=sedlay(j,k,issssil)-posol
             powtra(j,k,ipowasi)=powtra(j,k,ipowasi)+posol*pors2w(k)
         ENDIF
      ENDDO

!!! End dissolution opal

! CALCULATE CaCO3-CO3 CYCLE AND SIMULTANEOUS CO3-UNDERSATURATION DIFFUSION
!*************************************************************************
! COMPUTE NEW POWCAR=CARBONATE ION CONCENTRATION IN THE SEDIMENT
! FROM CHANGED ALKALINITY (NITRATE PRODUCTION DURING REMINERALISATION)
! AND DIC GAIN. ITERATE 5 TIMES. THIS CHANGES PH (SEDHPL) OF SEDIMENT.


! Add sediment flux of CaCO3
         IF(bolay(j).GT.0._wp) THEN
            sedlay(j,1,isssc12)=                                     &
     &      sedlay(j,1,isssc12)+prcaca(j)/(porsol(1)*seddw(1))
         !   prcaca(j)=0._wp
         ENDIF


      DO k = 1, ks

         IF((bolay(j).GT.0._wp).and.(pddpo(j,1)>0.5_wp)) THEN

               bt = rrrcl*psao(j,kbo(j))
             ! sulfate Morris & Riley (1966)
               sti   = 0.14_wp *  psao(j,kbo(j))*1.025_wp/1.80655_wp  / 96.062_wp
             ! fluoride Riley (1965)
               ft    = 0.000067_wp * psao(j,kbo(j))*1.025_wp/1.80655_wp / 18.9984_wp

               sit   = powtra(j,k,ipowasi) 
               pt    = powtra(j,k,ipowaph)

               alk  = powtra(j,k,ipowaal) 
               c    = powtra(j,k,ipowaic) 

               ak1  = ak13(j,kbo(j))
               ak2  = ak23(j,kbo(j))
               akb  = akb3(j,kbo(j))
               akw  = akw3(j,kbo(j))
               ah1  = sedhpl(j,k)
               aks  = aks3(j,kbo(j))
               akf  = akf3(j,kbo(j))
               ak1p = ak1p3(j,kbo(j))
               ak2p = ak2p3(j,kbo(j))
               ak3p = ak3p3(j,kbo(j))
               aksi = aksi3(j,kbo(j)) 

               iter  = 0
               DO jit = 1,20

                    hso4 = sti / ( 1._wp + aks / ( ah1 / ( 1._wp + sti / aks ) ))
                    hf   = 1._wp / ( 1._wp + akf / ah1 )
                    hsi  = 1._wp/ ( 1._wp + ah1 / aksi )
                    hpo4 = ( ak1p * ak2p * ( ah1 + 2._wp * ak3p ) - ah1**3 ) /& 
               &        ( ah1**3 + ak1p * ah1**2 + ak1p * ak2p * ah1 + ak1p *ak2p *ak3p )
                    ab   = bt / ( 1._wp + ah1 / akb )
                    aw   = akw / ah1 - ah1 / ( 1._wp + sti / aks )
                    ac   = alk + hso4 - sit * hsi - ab - aw + ft * hf - pt *hpo4
                    ah2o = SQRT( ( c - ac )**2 + 4._wp * ( ac * ak2 / ak1 ) * (2._wp *c - ac ) )
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

               if(ah1.gt.0._wp) then 
                 sedhpl(j,k)=max(1.e-20_wp,ah1)
               endif

                 powcar(k)  = c / (1._wp + sedhpl(j,k)/ak1 * (1._wp + sedhpl(j,k)/ak2))
              ENDIF
          END DO
      


! Evaluate boundary conditions for sediment-water column exchange.
! Current undersaturation of bottom water: sedb(i,0) and

! CO3 saturation concentration is aksp/calcon as in CARCHM
! (calcon defined in BELEG_BGC with 1.03e-2; 1/calcon =~ 97.)

! Calculate updated degradation rate from updated undersaturation.
! Calculate new solid sediment.
! No update of powcar pore water concentration from new undersaturation so far.
! Instead, only update DIC, and, of course, alkalinity.
! This also includes gains from aerobic and anaerobic decomposition.

      DO  k=1,ks
         IF((bolay(j).GT.0._wp).and.(pddpo(j,1)>0.5_wp)) THEN

            satlev=aksp(j,kbo(j))/calcon
!in oversaturated ( wrt calcite, powcar> satlev) water the "undersaturation" is negativ.
!there is evidence that in warm water like the Mediterranean spontaneous
!adsorption to existing calcareous shells happens. In most parts of the
!ocean this seems to be unlikely. Thus, we restrict on real undersaturation:

!!! emr used here only k=1 for undersaturation, seems to be wrong, 
!!! instead 1 use k
!!!            undsa = MAX(satlev - powcar(i, 1), 0._wp)
          undsa = MAX(satlev - powcar( k), 0._wp)
          sssnew = sedlay(j,k,isssc12)/(1._wp+ disso_cal*undsa)
          posol =  sedlay(j,k,isssc12) - sssnew

           sedlay(j,k,isssc12)=sedlay(j,k,isssc12)-posol
           powtra(j,k,ipowaic)=powtra(j,k,ipowaic)+posol*pors2w(k)
           powtra(j,k,ipowaal)=powtra(j,k,ipowaal)+2._wp*posol*pors2w(k)
         ENDIF
   ENDDO
 ENDDO
!HAMOCC_OMP_END_DO
!HAMOCC_OMP_END_PARALLEL

  CALL dipowa(start_idx,end_idx)

  DO j = start_idx, end_idx
        sedlay(j,1,issster) = sedlay(j,1,issster)                 &
             &                + produs(j)/(porsol(1)*seddw(1))
  ENDDO

      END SUBROUTINE POWACH
