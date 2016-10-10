!>
!! @file powach_impl.f90
!! @brief Computes sediment chemistry, implicit method, 
!!        calls powre water diffusion
!!
!!
!!
#include "hamocc_omp_definitions.inc"
SUBROUTINE powach_impl( start_idx, end_idx, psao )

  USE mo_kind, ONLY        : wp

  USE mo_carbch, ONLY      : bgctra, bgctend, satoxy,                     &
       &                     ak13, ak23, akb3, akw3, aksp, co3,           &
       &                     aks3,akf3,ak1p3,ak2p3,ak3p3,aksi3


  USE mo_sedmnt, ONLY      : sedlay, sedhpl, seddw, silpro, produs,      &
       &                     powtra, prcaca, prorca,              &
       &                     pown2bud, powh2obud,                        &
       &                     porsol, porwat, rno3, calcon,               &
       &                     sred_sed, ks, silsat,     &
       &                     seddenit,disso_op,disso_cal, &
       &                     sedtend, isremins, isremino, isreminn

  USE mo_biomod, ONLY      : bolay, kbo, ro2ut, rnit, nitdem, n2prod,    &
       &                     rrrcl, rcar, riron

  USE mo_control_bgc, ONLY : dtbgc, bgc_nproma, bgc_zlevs

  USE mo_param1_bgc, ONLY  : ipowasi, isilica, issssil, ipowaox, kaou,   &
       &                     ioxygen, issso12, ipowaph, ipowno3, ipown2, &
       &                     ipowaal, ipowaic, isssc12, ipowafe, issster

  USE mo_hamocc_nml, ONLY  : disso_po, denit_sed

  IMPLICIT NONE

  !! Arguments

  INTEGER, INTENT(in)  :: start_idx      !< start index for j loop (ICON cells, MPIOM lat dir)        
  INTEGER, INTENT(in)  :: end_idx        !< end index  for j loop  (ICON cells, MPIOM lat dir)         

  REAL(wp), INTENT(in) :: psao(bgc_nproma,bgc_zlevs)  !< salinity [psu.].

  !! Local variables

  INTEGER :: j,k, iter

  REAL(wp) :: sedb1(0:ks), sediso(0:ks)
  REAL(wp) :: solrat(ks), powcar(ks)
  REAL(wp) :: aerob(ks), anaerob(ks), ansulf(ks)

  REAL(wp) :: undsa, posol 
  REAL(wp) :: umfa, bt, alk, c
  REAL(wp) :: ak1, ak2, akb, akw
  REAL(wp) :: h, t1, t2, a, dadh, dddhhh, satlev
  REAL(wp) :: AKSI,AKF,AKS,AK1P,AK2P,AK3P
   REAL(wp) :: sti,ft,sit,pt,ah1       ! Chr. Heinze [H+]
   REAL(wp) :: hso4,hf,hsi,hpo4
   REAL(wp) :: ab,aw,ac,ah2o,ah2,erel,pco2_ch
   INTEGER  :: jit

  !
  REAL(wp) :: bolven
  !
  !----------------------------------------------------------------------
  !

  seddenit(:) = 0._wp
  solrat(:)   = 0._wp
  powcar(:)   = 0._wp
  anaerob(:)  = 0._wp
  aerob(:)    = 0._wp
  ansulf(:)   = 0._wp

  sedb1(:)    = 0._wp
  sediso(:)   = 0._wp

!!HAMOCC_OMP_PARALLEL
!!HAMOCC_OMP_DO PRIVATE(j,bolven,undsa,sedb1,solrat,posol,umfa,&
!!HAMOCC_OMP           aerob,anaerob,seddenit,ansulf,iter,k,bt,&
!!HAMMOC_OMP           alk,c,ak1,ak2,akb,akw,h,t1,t2,a,dadh,dddhhh,&
!!HAMOCC_OMP           satlev) HAMOCC_OMP_DEFAULT_SCHEDULE

  DO j = start_idx, end_idx
     ! calculate bottom ventilation rate for scaling of sediment-water exchange
      bolven = 1._wp


     ! CALCULATE SILICATE-OPAL CYCLE AND SIMULTANEOUS SILICATE DIFFUSION
     !******************************************************************

     ! Evaluate boundary conditions for sediment-water column exchange.
     ! Current undersaturation of bottom water: sedb(i,0) and
     ! Approximation for new solid sediment, as from sedimentation flux: solrat(1)

     IF (bolay(j) > 0._wp) THEN
         undsa=silsat-powtra(j,1,ipowasi)
         sedb1(0) = bolay(j)*(silsat-bgctra(j,kbo(j),isilica)) &
                &                 *bolven
         solrat(1)=                                                &
                &      (sedlay(j,1,issssil)+silpro(j)/(porsol(1)*seddw(1)))    &
                &      *disso_op/(1._wp + disso_op*undsa)*porsol(1)/porwat(1)
     ENDIF

     ! Evaluate sediment undersaturation and degradation.
     ! Current undersaturation in pore water: sedb(i,k) and
     ! Approximation for new solid sediment, as from degradation: solrat(k)

     DO k = 1, ks
           IF (bolay(j) > 0._wp) THEN
              undsa=silsat-powtra(j,k,ipowasi)
              sedb1(k)=seddw(k)*porwat(k)*(silsat-powtra(j,k,ipowasi))
              IF (k > 1) solrat(k) = sedlay(j,k,issssil)                 &
                   &                 * disso_op/(1._wp+disso_op*undsa)*porsol(k)/porwat(k)
           ENDIF
     END DO

     ! Solve for new undersaturation sediso, from current undersaturation sedb1,
     ! and first guess of new solid sediment solrat.

     CALL powadi(j,solrat(:),sedb1(:),sediso(:),bolven)

     ! Update water column silicate, and store the flux for budget.
     ! Add biogenic opal flux to top sediment layer.

     IF(bolay(j) > 0._wp) THEN

           bgctra(j,kbo(j),isilica) = silsat-sediso(0)
           sedlay(j,1,issssil)=                                    &
                &        sedlay(j,1,issssil)+silpro(j)/(porsol(1)*seddw(1))
           silpro(j)=0._wp
     ENDIF

     ! Calculate updated degradation rate from updated undersaturation.
     ! Calculate new solid sediment.
     ! Update pore water concentration from new undersaturation.

     DO k = 1, ks

           IF (bolay(j) > 0._wp) THEN
              solrat(k) = sedlay(j,k,issssil)                        &
                   &      * disso_op/(1._wp+disso_op*sediso(k))
              posol = sediso(k)*solrat(k)
              sedlay(j,k,issssil) = sedlay(j,k,issssil) - posol
              powtra(j,k,ipowasi) = silsat-sediso(k)
           ENDIF

     END DO


     ! CALCULATE OXYGEN-POC CYCLE AND SIMULTANEOUS OXYGEN DIFFUSION
     !*************************************************************

     ! This scheme is not based on undersaturation, but on O2 itself

     ! Evaluate boundary conditions for sediment-water column exchange.
     ! Current concentration of bottom water: sedb(i,0) and
     ! Approximation for new solid sediment, as from sedimentation flux: solrat(1)


        IF (bolay(j) > 0._wp) THEN
           undsa=powtra(j,1,ipowaox)
           sedb1(0)  = bolay(j)*bgctra(j,kbo(j),ioxygen)*bolven
           solrat(1) = (sedlay(j,1,issso12)+prorca(j) / (porsol(1)*seddw(1)))  &
                &      * ro2ut*disso_po/(1._wp+disso_po*undsa)*porsol(1)/porwat(1)
        ENDIF


     ! Evaluate sediment concentration and degradation.
     ! Current concentration in pore water: sedb(i,k) and
     ! Approximation for new solid sediment, as from degradation: solrat(k)

     DO k = 1, ks

           IF (bolay(j) > 0._wp) THEN
              undsa=powtra(j,k,ipowaox)
              sedb1(k)=seddw(k)*porwat(k)*powtra(j,k,ipowaox)
              IF (k > 1) solrat(k) = sedlay(j,k,issso12)               &
                   &                 * ro2ut*disso_po/(1._wp+disso_po*undsa) &
                   &                 * porsol(k)/porwat(k)
           ENDIF

     END DO

     ! Solve for new O2 concentration sediso, from current concentration sedb1,
     ! and first guess of new solid sediment solrat.

     CALL powadi(j,solrat(:),sedb1(:),sediso(:),bolven)

     ! Update water column oxygen, and store the flux for budget (opwflux). ! js: opwflux not in present model code
     ! Add organic carbon flux 'prorca' to top sediment layer.


        IF (bolay(j) > 0._wp) THEN
           bgctra(j,kbo(j),ioxygen)=sediso(0)
           bgctend(j,kbo(j),kaou) = satoxy(j,kbo(j)) - bgctra(j,kbo(j),ioxygen) ! update AOU
           sedlay(j,1,issso12)                                     &
                &      =sedlay(j,1,issso12)+prorca(j)/(porsol(1)*seddw(1))
        ENDIF


     ! Calculate updated degradation rate from updated concentration.
     ! Calculate new solid sediment.
     ! Update pore water concentration.
     ! Store flux in array aerob, for later computation of DIC and alkalinity.

     DO k = 1, ks
        umfa = porsol(k)/porwat(k)

           IF (bolay(j) > 0._wp) THEN
              solrat(k) = sedlay(j,k,issso12)                         &
                   &      * disso_po/(1._wp+disso_po*sediso(k))
              posol  = sediso(k)*solrat(k)
              aerob(k) = posol*umfa !this has P units: kmol P/m3 of pore water
              sedlay(j,k,issso12) = sedlay(j,k,issso12) - posol
              powtra(j,k,ipowaph) = powtra(j,k,ipowaph) + posol*umfa
              powtra(j,k,ipowno3) = powtra(j,k,ipowno3) + posol*rnit*umfa
              powtra(j,k,ipowaox) = sediso(k)

              sedtend(j,k,isremino) = posol*umfa/dtbgc
           ENDIF

     END DO

     ! CALCULATE NITRATE REDUCTION UNDER ANAEROBIC CONDITIONS EXPLICITELY
     !*******************************************************************

     ! Denitrification rate constant o
     ! Store flux in array anaerob, for later computation of DIC and alkalinity.

     DO k = 1, ks
        umfa = porsol(k)/porwat(k)

           IF (bolay( j) > 0._wp) THEN
              IF (powtra(j, k, ipowaox) < 1.e-6_wp) THEN

                 posol  = denit_sed * MIN(0.5_wp * powtra(j,k,ipowno3)/nitdem, &
                      &                       sedlay(j,k,issso12))
                 anaerob(k) = posol*umfa !this has P units: kmol P/m3 of pore water
                 seddenit(j) = seddenit(j) + nitdem*posol*umfa/dtbgc*seddw(k)
                 sedlay(j,k,issso12) = sedlay(j,k,issso12) - posol
                 powtra(j,k,ipowaph) = powtra(j,k,ipowaph) + posol*umfa

                 powtra(j,k,ipowno3) = powtra(j,k,ipowno3) - nitdem*posol*umfa
                 powtra(j,k,ipown2)  = powtra(j,k,ipown2)  + n2prod*posol*umfa
                 powh2obud(j,k)      = powh2obud(j,k)+0.5_wp*n2prod*posol*umfa

                 sedtend(j,k, isreminn) = posol*umfa/dtbgc
              ELSE
                 sedtend(j,k, isreminn) = 0._wp
              ENDIF   ! oxygen <1.e-6
           ENDIF   ! bolay

     END DO

     !    sulphate reduction in sediments
     DO k = 1, ks
        umfa = porsol(k)/porwat(k)

           IF (bolay(j) > 0._wp) THEN
              IF(powtra(j,k,ipowaox).LT.1.e-6_wp) THEN

                 posol  = sred_sed * sedlay(j,k,issso12)
                 ansulf(k) = posol*umfa !this has P units: kmol P/m3 of pore water
                 sedlay(j,k,issso12) = sedlay(j,k,issso12) - posol
                 powtra(j,k,ipowaph) = powtra(j,k,ipowaph) + posol*umfa
                 powtra(j,k,ipowno3) = powtra(j,k,ipowno3) + posol*umfa*rno3
                 powh2obud(j,k)      = powh2obud(j,k) - ro2ut*posol*umfa
                 !pown2bud update below
                 sedtend(j,k, isremins) = posol*umfa/dtbgc
              ELSE
                 sedtend(j,k, isremins) = 0._wp

              ENDIF
           ENDIF

     END DO


     ! CALCULATE CaCO3-CO3 CYCLE AND SIMULTANEOUS CO3-UNDERSATURATION DIFFUSION
     !*************************************************************************
     !
     ! COMPUTE NEW POWCAR=CARBONATE ION CONCENTRATION IN THE SEDIMENT
     ! FROM CHANGED ALKALINITY (NITRATE PRODUCTION DURING REMINERALISATION)
     ! AND DIC GAIN. THIS CHANGES PH (SEDHPL) OF SEDIMENT.

        DO k = 1, ks

         IF((bolay(j).GT.0._wp)) THEN

               bt = rrrcl*psao(j,kbo(j))
             ! sulfate Morris & Riley (1966)
               sti   = 0.14_wp *  psao(j,kbo(j))*1.025_wp/1.80655_wp  / 96.062_wp
             ! fluoride Riley (1965)
               ft    = 0.000067_wp * psao(j,kbo(j))*1.025_wp/1.80655_wp / 18.9984_wp

               sit   = powtra(j,k,ipowasi) 
               pt    = powtra(j,k,ipowaph)

               alk  = powtra(j,k,ipowaal) -( -ansulf(k)  +aerob(k))*rnit + nitdem*anaerob(k)
               c    = powtra(j,k,ipowaic) +(anaerob(k) +aerob(k) + ansulf(k))*rcar

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
     ! Approximation for new solid sediment, as from sedimentation flux: solrat(1)

     ! CO3 saturation concentration is aksp/calcon as in CARCHM
     ! (calcon defined in BELEG_BGC with 1.03e-2; 1/calcon =~ 97.)

        IF (bolay(j) > 0._wp) THEN
           satlev = aksp(j,kbo(j))/calcon+2.e-5_wp
           undsa  = MAX(satlev - powcar(1), 0._wp)
           sedb1(0) = bolay(j)*(satlev-co3(j,kbo(j)))  &
                &     * bolven
           solrat(1)= (sedlay(j,1,isssc12)+prcaca(j)/(porsol(1)*seddw(1)))  &
                &     * disso_cal/(1._wp + disso_cal*undsa)*porsol(1)/porwat(1)
        ENDIF

     ! Evaluate sediment undersaturation and degradation.
     ! Current undersaturation in pore water: sedb(i,k) and
     ! Approximation for new solid sediment, as from degradation: solrat(k)

     DO k = 1, ks
           IF (bolay(j) > 0._wp) THEN
              undsa = MAX(aksp( j, kbo( j)) / calcon - powcar(k), 0._wp)
              sedb1(k) = seddw(k)*porwat(k)*undsa
              IF (k > 1) solrat(k) = sedlay(j,k,isssc12)                 &
                   &                 * disso_cal/(1._wp+disso_cal*undsa)*porsol(k)/porwat(k)
              IF (undsa <= 0._wp) solrat(k) = 0._wp
           ENDIF
     END DO

     ! Solve for new undersaturation sediso, from current undersaturation sedb1,
     ! and first guess of new solid sediment solrat.

     CALL powadi(j,solrat(:),sedb1(:),sediso(:),bolven)

     ! There is no exchange between water and sediment with respect to co3 so far.
     ! Add calcite flux 'prcaca' to uppermost sediment layer.


    IF (bolay(j) > 0._wp) THEN
       sedlay(j,1,isssc12) =                                     &
               &      sedlay(j,1,isssc12) + prcaca(j)/(porsol(1)*seddw(1))
        prcaca(j) = 0._wp
     ENDIF


     ! Calculate updated degradation rate from updated undersaturation.
     ! Calculate new solid sediment.
     ! No update of powcar pore water concentration from new undersaturation so far.
     ! Instead, only update DIC, and, of course, alkalinity.
     ! This also includes gains from aerobic and anaerobic decomposition.

     DO k = 1, ks
        umfa = porsol(k)/porwat(k)


           IF (bolay(j) > 0._wp) THEN
              solrat(k) = sedlay(j,k,isssc12)                           &
                   &      * disso_cal/(1._wp+disso_cal*sediso(k))

              posol       = sediso(k)*solrat(k)

              sedlay(j,k,isssc12) = sedlay(j,k,isssc12)-posol

              powtra(j,k,ipowaic) = powtra(j,k,ipowaic)                 &
                   &                + posol*umfa+(aerob(k)                &
                   &                + anaerob(k) + ansulf(k))*rcar

              powtra(j,k,ipowaal) = powtra(j,k,ipowaal)                 &
                   &                + 2._wp*posol*umfa - rnit*(aerob(k)   &
                   &                - ansulf(k)) + nitdem*anaerob(k)

              pown2bud(j,k)       = pown2bud(j,k) +  2._wp*n2prod*anaerob(k) &
                                  + 2._wp*rnit*ansulf(k)

              powtra(j,k,ipowafe) = powtra(j,k,ipowafe)                 &
                   &                + (aerob(k)+anaerob(k)+ansulf(k))*riron
           ENDIF

     END DO

  END DO ! cells 
!!HAMOCC_OMP_END_DO
!!HAMOCC_OMP_END_PARALLEL

  CALL dipowa(start_idx,end_idx)

  DO j = start_idx, end_idx
        sedlay(j,1,issster) = sedlay(j,1,issster)                 &
             &                + produs(j)/(porsol(1)*seddw(1))
  ENDDO

!   DO j = start_idx, end_idx
!         silpro(j) = 0._wp
!         prorca(j) = 0._wp
!         prcaca(j) = 0._wp
!         produs(j) = 0._wp
!   END DO
! 
END SUBROUTINE powach_impl
