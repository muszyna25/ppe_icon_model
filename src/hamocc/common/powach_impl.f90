!>
!! @file powach.f90
!! @brief compute biological production, settling of debris, and related biogeochemistry
!!
!!
!! called by bgc
!!
!! @author Ernst Maier-Reimer, MPI-Met, HH
!!
!! @par Revision History
!!
!! First version by Ernst Maier-Reimer     (MPI-M)     Apr 10, 2001
!!
!!
SUBROUTINE powach_impl( start_idx, end_idx, klevs, psao )

  USE mo_kind, ONLY        : wp

  USE mo_carbch, ONLY      : bgctra, sedfluxo,                           &
       &                     ak13, ak23, akb3, akw3, aksp, co3

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

  USE mo_param1_bgc, ONLY  : ipowasi, isilica, issssil, ipowaox,     &
       &                     ioxygen, issso12, ipowaph, ipowno3, ipown2, &
       &                     ipowaal, ipowaic, isssc12, ipowafe, issster

  USE mo_hamocc_nml, ONLY  : disso_po, denit_sed

  IMPLICIT NONE

  !! Arguments

  INTEGER, INTENT(in)  :: start_idx                  !< 1st REAL of model grid.
  INTEGER, INTENT(in)  :: end_idx                  !< 2nd REAL of model grid.
  INTEGER, INTENT(in)  :: klevs(bgc_nproma)                  !< 3rd (vertical) REAL of model grid.

  REAL(wp), INTENT(in) :: psao(bgc_nproma,bgc_zlevs)  !< salinity [psu.].

  !! Local variables

  INTEGER :: j,k, iter, kpke

  REAL(wp) :: sedb1(0:ks), sediso(0:ks)
  REAL(wp) :: solrat(ks), powcar(ks)
  REAL(wp) :: aerob(ks), anaerob(ks), ansulf(ks)

  REAL(wp) :: undsa, posol 
  REAL(wp) :: umfa, bt, alk, c
  REAL(wp) :: ak1, ak2, akb, akw
  REAL(wp) :: h, t1, t2, a, dadh, dddhhh, satlev

  ! *****************************************************************
  ! accelerated sediment
  ! needed for boundary layer ventilation in fast sediment routine
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

  DO j = start_idx, end_idx
     prorca(j)=0._wp
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
     !      sedfluxo(j,ipowasi) = sedfluxo(j,ipowasi) +                &
     !           &  (silsat-sediso(0)-bgctra(j,kbo(j),isilica))*bolay(j)

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
           sedlay(j,1,issso12)                                     &
                &      =sedlay(j,1,issso12)+prorca(j)/(porsol(1)*seddw(1))
           prorca(j) = 0._wp
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

              sedtend(j,k,isremino) = posol/dtbgc
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

                 !tk27012010 changed no3 use and n2 production in denitrification
                 powtra(j,k,ipowno3) = powtra(j,k,ipowno3) - nitdem*posol*umfa
                 powtra(j,k,ipown2)  = powtra(j,k,ipown2)  + n2prod*posol*umfa
                 powh2obud(j,k)      = powh2obud(j,k)+0.5_wp*n2prod*posol*umfa

                 sedtend(j,k, isreminn) = posol/dtbgc
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
                ! seddiag(j,k,ksred_sed) = posol*umfa/dtbgc
                 sedlay(j,k,issso12) = sedlay(j,k,issso12) - posol
                 powtra(j,k,ipowaph) = powtra(j,k,ipowaph) + posol*umfa
                 powtra(j,k,ipowno3) = powtra(j,k,ipowno3) + posol*umfa*rno3
                 powh2obud(j,k)      = powh2obud(j,k) - ro2ut*posol*umfa
                 !pown2bud update below
                 sedtend(j,k, isremins) = posol/dtbgc
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
     ! AND DIC GAIN. ITERATE 5 TIMES. THIS CHANGES PH (SEDHPL) OF SEDIMENT.

     DO iter = 1, 5

        DO k = 1, ks

              IF (bolay(j) > 0._wp) THEN
                 bt = rrrcl*psao(j,kbo(j))
                 !alkalinity is increased during denitrification due to consumption of H+ via NO3 (see Wolf-Gladrow etal,2007)
                 alk  = powtra(j,k,ipowaal) -( -ansulf(k)  +aerob(k))*rnit + nitdem*anaerob(k)
                 c    = powtra(j,k,ipowaic) +(anaerob(k) +aerob(k) + ansulf(k))*rcar
                 ak1  = ak13(j,kbo(j))
                 ak2  = ak23(j,kbo(j))
                 akb  = akb3(j,kbo(j))
                 akw  = akw3(j,kbo(j))
                 h    = sedhpl(j,k)
                 t1   = h/ak1
                 t2   = h/ak2
                 a    = c*(2._wp+t2)/(1._wp+t2+t2*t1) + akw/h-h + bt/(1._wp + h/akb) - alk
                 dadh = c*( 1._wp/(ak2*(1._wp+t2+t2*t1))-(2._wp+t2)*(1._wp/ak2+2._wp*t1/ak2)/  &
                      &          (1._wp+t2+t2*t1)**2)                                        &
                      &          -akw/h**2-1._wp-(bt/akb)/(1._wp+h/akb)**2
                 dddhhh = a/dadh
                 sedhpl(j,k) = MAX(h - dddhhh, 1.e-11_wp)
                 powcar(k)  = c / (1._wp + t2 * (1._wp + t1))
              ENDIF

           END DO

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

  CALL dipowa(start_idx,end_idx)

  ! f(POC) [kg C] / f(total) [kg] = 0.05
  ! thus it is
  DO j = start_idx, end_idx
        sedlay(j,1,issster) = sedlay(j,1,issster)                 &
             &                + produs(j)/(porsol(1)*seddw(1))
  ENDDO

  DO j = start_idx, end_idx
        silpro(j) = 0._wp
        prorca(j) = 0._wp
        prcaca(j) = 0._wp
        produs(j) = 0._wp
  END DO

END SUBROUTINE powach_impl
