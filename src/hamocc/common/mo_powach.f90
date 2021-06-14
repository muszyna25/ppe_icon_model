!! @file mo_powach.f90
!! @brief sediment chemistry, implicit and explicit discretisation


#include "hamocc_omp_definitions.inc"
MODULE mo_powach

 USE mo_kind, ONLY        : wp
 USE mo_sedmnt_diffusion, ONLY: powadi, dipowa
 USE mo_hamocc_nml,       ONLY: ks,porwat
 USE mo_carchm,           ONLY: update_hi

 IMPLICIT NONE

  PRIVATE

  PUBLIC :: powach, &
            powach_impl


CONTAINS
   SUBROUTINE POWACH(start_idx,end_idx,psao,pddpo)
!!   @brief compute sediment chemistry, explicit discretisation
!!   call pore water diffusion
!!

! 

   USE mo_memory_bgc, ONLY  : ak13, ak23, akb3, akw3, aksp,               &
        &                     aks3,akf3,ak1p3,ak2p3,ak3p3,aksi3,          &          
        &                     bolay, kbo, ro2ut, rnit, nitdem, n2prod,    &
        &                     rrrcl, rcar, ralk, riron,                   &
        &                     nitrira, ro2ammo, anamoxra, bkno2, nitriox, &
        &                     no2denit, rnh4no2, rno2no3, rno3no2,rno3nh4,&
        &                     rno2n2, ro2nitri, alk_nrn2,                 &
        &                     o2thresh, o2den_lim

   USE mo_sedmnt, ONLY      : sedlay, sedhpl, seddw, silpro,              &
        &                     powtra, prcaca, prorca, produs,             &
        &                     pown2bud, powh2obud, sred_sed,              &
        &                     porsol, pors2w, calcon,                     &
        &                     disso_op, disso_cal,                        &
        &                     sedtend,      &
        &                     silsat

 
   USE mo_control_bgc, ONLY : dtbgc, bgc_nproma, bgc_zlevs
 
   USE mo_param1_bgc, ONLY  : ipowasi, issssil, ipowaox,ipowh2s,  &
        &                     issso12, ipowaph, ipowno3, ipown2, &
        &                     ipowaal, ipowaic, isssc12, issster, &
        &                     ipowafe, isremins, isremino, isreminn, &
        &                     ipownh4, ipowno2, ksammox, ksanam, &
        &                     ksdnra, ksdnrn, ksnrn2, ksnitox
 
   USE mo_hamocc_nml, ONLY  : disso_po,denit_sed, &
                              l_N_cycle, no3no2red, no3nh4red
 
  IMPLICIT NONE

  !! Arguments

   INTEGER, INTENT(in)  :: start_idx     !< start index for j loop (ICON cells, MPIOM lat dir)          
   INTEGER, INTENT(in)  :: end_idx        !< end index  for j loop  (ICON cells, MPIOM lat dir)        

   REAL(wp), INTENT(in) :: psao(bgc_nproma,bgc_zlevs)   !< salinity [psu.].
   REAL(wp), INTENT(in) :: pddpo(bgc_nproma,bgc_zlevs)  !< size of scalar grid cell 
  
  !! Local variables

   INTEGER :: j,k


   REAL(wp) :: orgsed                        ! sum of C12 and C13 in organic sed.
   REAL(wp) :: sssnew                        ! temporarily value of solid constituent after solution

   REAL(wp) :: powcar(ks)

   REAL(wp) :: undsa, posol, pomax
   REAL(wp) :: satlev
!   REAL(wp) :: o2lim

   !!!! extended N-cycle variables
   REAL(wp) :: posol_nit
   REAL(wp) :: dissot1
   REAL(wp) :: popot
   REAL(wp) :: o2lim
   REAL(wp) :: nitrif          ! local rate of NH4 conversion to NO3 (adapted from EMR)
   REAL(wp) :: ammox,nitox     ! oxidation of NH4 to NO2 and NO2 to NO3, light dep.
   REAL(wp) :: fammox,fnitox   ! local rate of ammox,nitox
   REAL(wp) :: detn,dnrn,dnra,nrn2,anamox ! processes of N-cycle
   REAL(wp) :: nh4a,no2a,no3a,nh4n       ! old NH4,NO2,NO3 after biological processes
   REAL(wp) :: no2c_act       ! absolut  change in no2 due to anamox and denitrification
   REAL(wp) :: no2c_max       ! max  change in no2 due to anamox and denitrification
   REAL(wp) :: no2rmax        ! max potential rate of change no2 due to anamox and denitrification
   REAL(wp) :: denpot         ! moa potential change in detn
   REAL(wp) :: fdnrn,fdnra,fnrn2,fdenit! fraction of total processes of N-cycle
   REAL(wp) :: annpot,anam     ! potential anamox in anaerobic water
   REAL(wp) :: oxpot      ! potential oxdiation of NH4 and NO2 in aerobic water
   REAL(wp) :: no3c_act       ! absolute change in NO3 due to DNRN/A and auto.denitrification  
   REAL(wp) :: no3c_max       ! max change in NO3 due to DNRN/A and auto.denitrification  
   REAL(wp) :: no3rmax        ! max potential rate of change in NO3 due to DNRN/A and auto.denitrification  
   REAL(wp) :: detc_act       ! absolute change in NO3 due to DNRN/A 
   REAL(wp) :: detc_max       ! max potential change in det due to DNRN/A
   REAL(wp) :: rdnrn,rdnra    ! local rates of DNRN/A
   REAL(wp) :: oxmax,oxact,r_nitox,r_ammox ! all for nitrification of NH4 and NO2

   REAL(wp) :: denlim_no3     ! saturation function for denitrification on low NO3
   REAL(wp) :: quadno3        ! for saturation function quadartic no3
   REAL(wp) :: newammo,newnitr

!!! WE start with remineralisation of organic to estimate alkalinity changes first
!          
!HAMOCC_OMP_PARALLEL
!HAMOCC_OMP_DO PRIVATE(j,k,pomax,posol,orgsed,sssnew,&
!HAMOCC_OMP           undsa,satlev,powcar) HAMOCC_OMP_DEFAULT_SCHEDULE


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

            !!!! N-cycle !!!!!!!!
            IF (l_N_cycle) THEN
               o2lim = powtra(j,k,ipowaox)/(o2thresh + powtra(j,k,ipowaox)) ! o2 limitation in oxic water
               pomax = o2lim*dissot1*powtra(j,k,ipowaox)
               sssnew = sedlay(j,k,issso12)/(1._wp + pomax)
               popot = sedlay(j,k,issso12) - sssnew   ! potential change for org sed.
               posol = min(0.9_wp*powtra(j,k,ipowaox)/(pors2w(k)*ro2ammo),popot)
            ELSE
               ! o2lim = powtra(j,k,ipowaox)/(o2thresh+powtra(j,k,ipowaox)) ! o2 limitation in oxic water

               ! maximal possible aerobe dissolution per time step 
               ! limited by available oxygen concentration, currently max 70 % of O2 
               pomax = disso_po * max(0._wp,sedlay(j,k,issso12))*powtra(j,k,ipowaox)
               posol = min(0.7_wp*powtra(j,k,ipowaox)/ro2ut,pomax*pors2w(k))
            ENDIF

         
            sedlay(j,k,issso12) = sedlay(j,k,issso12) - posol

            powtra(j,k,ipowaic) = powtra(j,k,ipowaic) + posol*rcar*pors2w(k)
            powtra(j,k,ipowaph)=powtra(j,k,ipowaph)+ posol*pors2w(k)
            powtra(j,k,ipowafe)=powtra(j,k,ipowafe)+ posol*riron*pors2w(k)
            !!!! N-cycle !!!!!!!!
            IF (l_N_cycle) THEN
               powtra(j,k,ipownh4)=powtra(j,k,ipownh4)+posol*rnit*pors2w(k)
               powtra(j,k,ipowaal)=powtra(j,k,ipowaal)+posol*(rnit - 1._wp)*pors2w(k)  !LR: watch out!
               powtra(j,k,ipowaox)=powtra(j,k,ipowaox)-posol*pors2w(k)*ro2ammo
            ELSE
               powtra(j,k,ipowno3)=powtra(j,k,ipowno3)+ posol*rnit*pors2w(k)
               powtra(j,k,ipowaal)=powtra(j,k,ipowaal)-posol*ralk*pors2w(k)
               powtra(j,k,ipowaox)=powtra(j,k,ipowaox)-posol*pors2w(k)*ro2ut
            ENDIF

            sedtend(j,k,isremino) = posol*pors2w(k)/dtbgc


            !!!! N-cycle !!!!!!!!
            !!!! Oxidation of NH4 and NO2
            IF (l_N_cycle) THEN
                fammox = nitrira          ! local rate of ammox per time step, no light dependence
                fnitox = nitriox          ! local rate of nitrite

                nh4a  =  max(0._wp,powtra(j,k,ipownh4))                         

                newammo = nh4a/(1._wp + fammox)
                ammox = nh4a - newammo
                no2a =  max(0._wp,powtra(j,k,ipowno2))  
                newnitr = no2a/ (1._wp + fnitox)                        ! change of nitrite

                nitox = no2a - newnitr

                !  ratio of change in o2 due to oxidation
                oxmax =  rno2no3*nitox +  rnh4no2*ammox   ! potential oxidation of NO2
                if  (oxmax > 0._wp) then
                   r_nitox =   rno2no3*nitox/ oxmax
                   r_ammox =   rnh4no2*ammox /oxmax
                else
                   r_nitox = 0._wp
                   r_ammox = 0._wp
                   oxmax = 0._wp
                endif
               !        oxact always > 0 because we are in o2 > 2
                oxact = min(powtra(j,k,ipowaox) - 0.5E-6_wp, oxmax)

                nitox = r_nitox*oxact/rno2no3
                ammox = r_ammox*oxact/rnh4no2
 
                powtra(j,k,ipownh4) = nh4a -ammox
                powtra(j,k,ipowno2) = powtra(j,k,ipowno2) + ammox-nitox ! change of nitrite
        
                powtra(j,k,ipowno3) = powtra(j,k,ipowno3) + nitox

                powtra(j,k,ipowaox) = powtra(j,k,ipowaox)- rno2no3*nitox &
                                       - rnh4no2*ammox ! O2 will be used during nitrification 

                powtra(j,k,ipowaal) = powtra(j,k,ipowaal) - 2._wp*ammox ! ocean with NH4 - alkalinity change
                                                                            ! according to Wolf-Gladrow etal Mar. Chem. (2007)

                sedtend(j,k,ksammox)=ammox/dtbgc
                sedtend(j,k,ksnitox)=nitox/dtbgc
            ENDIF

         ELSE

            sedtend(j,k,isremino) = 0._wp

            !!!! N-cycle !!!!!!!!
            IF (l_N_cycle) THEN
               sedtend(j,k,ksammox) = 0._wp
                sedtend(j,k,ksnitox) = 0._wp
            ENDIF

         ENDIF ! oxygen > 2.
         ENDIF ! bolay > 0

      ENDDO





! CALCULATE NITRATE REDUCTION UNDER ANAEROBIC CONDITIONS EXPLICITELY
!*******************************************************************
! Denitrification rate constant of POP (disso) [1/sec]
! Store flux in array anaerob, for later computation of DIC and alkalinity.

      DO  k=1,ks
         IF (bolay( j) .GT. 0._wp) THEN

         !!!! not N-cycle !!!!!!!!
         IF (.not. l_N_cycle) THEN
         IF (powtra( j, k, ipowaox) < 2.e-6_wp) THEN

           orgsed = max(0._wp,sedlay(j,k,issso12)) 

           posol = denit_sed * MIN(0.5_wp * powtra(j, k, ipowno3)/(nitdem-rnit), orgsed)
           sedlay(j,k,issso12)=sedlay(j,k,issso12)-posol
           powtra(j,k,ipowaph)=powtra(j,k,ipowaph)+posol*pors2w(k)
           powtra(j,k,ipowaic)=powtra(j,k,ipowaic)+rcar*posol*pors2w(k)

           powtra(j,k,ipowaal)=powtra(j,k,ipowaal)+posol*(nitdem-ralk)*pors2w(k)
           powtra(j,k,ipowno3)=powtra(j,k,ipowno3)-(2._wp*n2prod - rnit)*posol*pors2w(k)
           powtra(j,k,ipowafe)=powtra(j,k,ipowafe)+ posol*riron*pors2w(k)
           powtra(j,k,ipown2)=powtra(j,k,ipown2)+n2prod*posol*pors2w(k)
           powh2obud(j,k)=powh2obud(j,k)+0.5_wp*n2prod*posol*pors2w(k)
           pown2bud(j,k) = pown2bud(j,k) + 2._wp*n2prod*posol*pors2w(k)

           sedtend(j,k,isreminn) = posol*pors2w(k)/dtbgc
         else
           sedtend(j,k,isreminn) = 0._wp

         ENDIF   ! oxygen <1.e-6

         ELSE

         !!!! N-cycle !!!!!!!!
         IF (powtra( j, k, ipowaox) < o2thresh) THEN
            orgsed=max(0._wp,sedlay(j,k,issso12))

            ! o2 limitation identical for all suboxic processes 
            o2lim = 1._wp - max(0._wp,powtra(j,k, ipowaox)/o2thresh)

            ! convert detritus in P-units to N-units for nitrogen cycle changes  
            ! convert to from solid to water
            detn = max(0._wp,orgsed*rnit*pors2w(k))

            ! denitrification rate on NO3
            rdnrn = o2lim*no3no2red *detn/(powtra(j,k,ipowno3)+3.E-5_wp)
            rdnra = o2lim*no3nh4red *detn/(powtra(j,k,ipowno3)+3.E-5_wp)
 
            no3rmax = rdnrn + rdnra                       ! max pot loss of NO3

            IF(no3rmax > 0._wp) THEN 
               fdnrn = rdnrn/no3rmax                       ! fraction each process
               fdnra = rdnra/no3rmax

               !< implicit formulation to avoid neg. nitrate concentration
               no3a = powtra(j,k,ipowno3)/(1._wp +no3rmax)   ! max change in NO3  
               no3c_max = powtra(j,k,ipowno3) - no3a         ! corresponding max NO3 loss 
               detc_max=  no3c_max*(fdnrn/rno3no2+fdnra/rno3nh4)  ! corresponding max change in det in water part  
               detc_act = min ( detn ,detc_max) ! convert solid to water part

               dnrn = fdnrn*detc_act             ! in P units in water part
               dnra = fdnra*detc_act             ! in P untis in water part

               posol_nit = dnrn + dnra      ! change for DIC and PO4

               sedlay(j,k,issso12) = sedlay(j,k,issso12) -(dnrn + dnra)/pors2w(k)

               powtra(j,k,ipowno3)= powtra(j,k,ipowno3)    &    ! change in nitrate 
                            &       -rno3no2*dnrn          &    ! from DNRN
                            &       -rno3nh4*dnra               ! from DNRA

               powtra(j,k,ipownh4)= powtra(j,k,ipownh4)    &    ! change in ammonium 
                            &       +rnit*dnrn             &    ! from DNRN
                            &       +86._wp*dnra                ! from DNRA

               powtra(j,k,ipowno2) = powtra(j,k,ipowno2)   &    ! change in nitrite 
                             &       + rno3no2*dnrn             ! from DNRN

               powtra(j,k,ipowaal)=powtra(j,k,ipowaal)      &   ! change from DNRN and DNRA
                           &       + rnit*dnrn              &   ! from DNRN
                           &       + 86._wp*dnra            &   ! from DNRA
                           &       - posol_nit

               pown2bud(j, k) = pown2bud(j, k) - 70._wp*dnra


               powtra(j,k,ipowaph) = powtra(j,k,ipowaph) + posol_nit
               powtra(j,k,ipowaic) = powtra(j,k,ipowaic) + rcar*posol_nit
               powtra(j,k,ipowafe) = powtra(j,k,ipowafe) + riron*posol_nit

               sedtend(j,k,ksdnrn)=dnrn/dtbgc
               sedtend(j,k,ksdnra)=dnra/dtbgc

            ELSE
               sedtend(j,k,ksdnrn)=0._wp
               sedtend(j,k,ksdnra)=0._wp
               posol_nit = 0._wp
            ENDIF ! no3rmax > 0

         ENDIF   ! oxygen <1.e-6
         ENDIF   ! l_N_cycle
         ENDIF   ! bolay
        ENDDO





      !!!! N-cycle !!!!!!!!
      IF (l_N_cycle) THEN
      DO  k=1,ks
         IF (bolay( j) .GT. 0._wp) THEN

! DENITRIFICATION on NO2
         IF (powtra(j,k,ipowaox) < o2den_lim) THEN
            o2lim = 1._wp - max(0._wp,powtra(j,k, ipowaox)/o2thresh)
            orgsed=max(0._wp,sedlay(j,k,issso12))
            detn = orgsed*rnit*pors2w(k)   ! converted to water equiv. part
            no2rmax = o2lim*no2denit*detn/(powtra(j,k,ipowno2)+0.1E-6_wp)

            ! implicit formulation to avoid neg. nitrite concentration
            no2a = powtra(j,k,ipowno2)/(1._wp+no2rmax) 

            no2c_max = powtra(j,k,ipowno2) -no2a          ! maximal NO2 loss
            detc_max=  no2c_max/rno2n2                    ! corresponding max change in det;
                                                          ! rno2n2 conversion to P units 

            nrn2 = min (orgsed*pors2w(k),detc_max)        ! in P units
            ! changes in detritus in P-units 
            posol = nrn2

            powtra(j,k,ipowno2)= powtra(j,k,ipowno2) - rno2n2*nrn2  ! change in ammonium             
                                                                    ! from nitrite reduction to N2; NRN2

            powtra(j,k,ipownh4)= powtra(j,k,ipownh4) + rnit*nrn2    ! change in ammonium             
                                                                    ! from nitrite reduction to N2; NRN2

            powtra(j,k,ipown2)= powtra(j,k,ipown2) + rno2n2*nrn2/2._wp

            powtra(j,k,ipowaal)=powtra(j,k,ipowaal) + alk_nrn2*nrn2 - posol
            pown2bud(j,k) = pown2bud(j,k)  + (alk_nrn2-rnit)*nrn2

            powh2obud(j,k)= powh2obud(j,k) + rno2n2*nrn2*0.25_wp
            ! now change in DIC, PO4 and Det in Sediment 
            sedlay(j,k,issso12)=sedlay(j,k,issso12) -posol/pors2w(k)

            powtra(j,k,ipowaph)=powtra(j,k,ipowaph)+posol
            powtra(j,k,ipowaic)=powtra(j,k,ipowaic)+rcar*posol
            powtra(j,k,ipowafe)=powtra(j,k,ipowafe)+riron*posol

            sedtend(j,k,ksnrn2)=nrn2/(dtbgc*pors2w(k)) ! change in ssso12; P units
         ELSE
            sedtend(j,k,ksnrn2) = 0._wp
         ENDIF   ! oxygen <1.e-6

! ANAMMOX on NO2 and NH4
         IF (powtra(j,k,ipowaox) < o2thresh) THEN

            nh4a = max(0._wp,powtra(j,k,ipownh4))
            no2a = max(0._wp,powtra(j,k,ipowno2))

            o2lim = 1._wp - max(0._wp,powtra(j,k, ipowaox)/o2thresh)
            anam = o2lim*anamoxra*no2a/(no2a+bkno2)   

            nh4n= nh4a/(1._wp + anam)
            anamox = nh4a - nh4n 
            anamox = min(anamox,no2a/1.3_wp)

            powtra(j,k,ipownh4)= powtra(j,k,ipownh4) -anamox
         
            powtra(j,k,ipowno2) =powtra(j,k,ipowno2) - 1.3_wp*anamox     ! from anamox          
            powtra(j,k,ipown2)= powtra(j,k,ipown2) + anamox              ! from anamox 
            powtra(j,k,ipowno3)=powtra(j,k,ipowno3) + 0.3_wp*anamox
            powtra(j,k,ipowaal)=powtra(j,k,ipowaal) - 0.3_wp*anamox

            ! loss of Os from NO2 - gain NO3 - 0.5 from NH4 =1.3 - 0.3*1.5- 0.5 = +0.35   
            powh2obud(j,k)= powh2obud(j,k) + anamox*0.35_wp                         

            pown2bud(j,k) = pown2bud(j,k) + 1.7_wp*anamox ! alk is unchanged, but for alk mass we need anammox

            sedtend(j,k,ksanam)=2.0_wp*anamox / dtbgc  ! if given in TgN must be doubled 
         ELSE
            sedtend(j,k,ksanam)=0._wp
         ENDIF   ! oxygen <1.e-6

         ENDIF   ! bolay
        ENDDO
        ENDIF
        !!!! N-cycle !!!!!!!!



! sulphate reduction in sediments
! Denitrification rate constant of POP (disso) [1/sec]

      DO  k=1,ks
         IF (bolay( j) > 0._wp) THEN
         IF (.not. l_N_cycle .and. powtra(j,k,ipowaox)<1.e-6_wp .or. &
            & l_N_cycle .and. powtra(j,k,ipowaox)<o2den_lim .and. &
            & powtra(j,k,ipowno3) < 30.e-6_wp) THEN
         
           orgsed=max(0._wp,sedlay(j,k,issso12))
           ! reduced by factor 100
           sssnew = orgsed/(1._wp + sred_sed)
           posol=orgsed - sssnew   ! change for org sed.

           sedlay(j,k,issso12)=sedlay(j,k,issso12)-posol
           powtra(j,k,ipowaic)=powtra(j,k,ipowaic)+posol*pors2w(k)*rcar
           powtra(j,k,ipowaph)=powtra(j,k,ipowaph)+posol*pors2w(k)
           powtra(j,k,ipowafe)=powtra(j,k,ipowafe)+ posol*riron*pors2w(k)

           IF (l_N_cycle) THEN
              ! according to Thamdrup ammonium might be oxidized by sulfur to form n2
              ! no change in water  and a smaler (32 instead of 48) alk change 
              powtra(j,k,ipown2) = powtra(j,k,ipown2) + 0.5_wp*rnit*posol*pors2w(k)
              powtra(j,k,ipowaal) = powtra(j,k,ipowaal) + posol*(2._wp*rnit - 1._wp)*pors2w(k) ! alk change is +32
              pown2bud(j,k) = pown2bud(j, k) +3._wp*posol*rnit*pors2w(k) 
              powh2obud(j,k) = powh2obud(j,k)-posol*(ro2ammo+0.5_wp*rnit)*pors2w(k)
           ELSE
              powtra(j,k,ipowno3)=powtra(j,k,ipowno3)+posol*rnit*pors2w(k)
              powtra(j,k,ipowaal)=powtra(j,k,ipowaal)+posol*ralk*pors2w(k) 
              powtra(j,k,ipowh2s) = powtra(j,k,ipowh2s) + posol*pors2w(k)*ralk
              powh2obud(j,k)=powh2obud(j,k)-posol*ro2ut*pors2w(k)
              pown2bud(j,k) = pown2bud(j,k) + 2._wp*ralk*posol*pors2w(k)
           ENDIF

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
               sedhpl(j,k)= update_hi(sedhpl(j,k),powtra(j,k,ipowaic),ak13(j,kbo(j)),&
        &                             ak23(j,kbo(j)),akw3(j,kbo(j)),aks3(j,kbo(j)),&
        &                             akf3(j,kbo(j)),aksi3(j,kbo(j)),ak1p3(j,kbo(j)),&
        &                             ak2p3(j,kbo(j)),ak3p3(j,kbo(j)),psao(j,kbo(j)),&
        &                             akb3(j,kbo(j)),powtra(j,k,ipowasi),powtra(j,k,ipowaph),&
        &                             powtra(j,k,ipowaal))

                powcar(k)  = powtra(j,k,ipowaic) / (1._wp + sedhpl(j,k)/ak13(j,kbo(j)) &
        &                    * (1._wp + sedhpl(j,k)/ak23(j,kbo(j))))
         ENDIF
      END DO
      


! Evaluate boundary conditions for sediment-water column exchange.
! Current undersaturation of bottom water: sedb(i,0) and

! CO3 saturation concentration is aksp/calcon as in mo_carchm.f90
! (calcon defined in mo_ini_bgc.f90 with 1.03e-2; 1/calcon =~ 97.)

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
SUBROUTINE powach_impl( start_idx, end_idx, psao )
!>
!! @brief Computes sediment chemistry, implicit method, 
!!        calls powre water diffusion
!!
!!
!!


  USE mo_memory_bgc, ONLY   : bgctra, bgctend, satoxy,                   &
       &                     ak13, ak23, akb3, akw3, aksp, co3,          &
       &                     aks3,akf3,ak1p3,ak2p3,ak3p3,aksi3,          &
       &                     bolay, kbo, ro2ut, rnit, nitdem, n2prod,    &
       &                     rrrcl, rcar, riron, ralk


  USE mo_sedmnt, ONLY      : sedlay, sedhpl, seddw, silpro, produs,      &
       &                     powtra, prcaca, prorca,              &
       &                     pown2bud, powh2obud,                        &
       &                     porsol, rno3, calcon,               &
       &                     sred_sed, silsat,     &
       &                     seddenit,disso_op,disso_cal, &
       &                     sedtend


  USE mo_control_bgc, ONLY : dtbgc, bgc_nproma, bgc_zlevs

  USE mo_param1_bgc, ONLY  : ipowasi, isilica, issssil, ipowaox, kaou,   &
       &                     ioxygen, issso12, ipowaph, ipowno3, ipown2, &
       &                     ipowaal, ipowaic, isssc12, ipowafe, issster,&
       &                     ipowh2s, isremins, isremino, isreminn

  USE mo_hamocc_nml, ONLY  : disso_po, denit_sed

  IMPLICIT NONE

  !! Arguments

  INTEGER, INTENT(in)  :: start_idx      !< start index for j loop (ICON cells, MPIOM lat dir)        
  INTEGER, INTENT(in)  :: end_idx        !< end index  for j loop  (ICON cells, MPIOM lat dir)         

  REAL(wp), INTENT(in) :: psao(bgc_nproma,bgc_zlevs)  !< salinity [psu.].

  !! Local variables

  INTEGER :: j,k

  REAL(wp) :: sedb1(0:ks), sediso(0:ks)
  REAL(wp) :: solrat(ks), powcar(ks)
  REAL(wp) :: aerob(ks), anaerob(ks), ansulf(ks)

  REAL(wp) :: undsa, posol 
  REAL(wp) :: umfa, alk, c
  REAL(wp) :: satlev

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
!!HAMOCC_OMP           aerob,anaerob,seddenit,ansulf,k,&
!!HAMMOC_OMP           alk,c,t1,t2,a,satlev) HAMOCC_OMP_DEFAULT_SCHEDULE

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



        IF (bolay(j) > 0._wp) THEN
           bgctra(j,kbo(j),ioxygen)=sediso(0)
           bgctend(j,kbo(j),kaou) = satoxy(j,kbo(j)) - bgctra(j,kbo(j),ioxygen) ! update AOU
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

                 posol  = denit_sed * MIN(0.5_wp * powtra(j,k,ipowno3)/(nitdem-rnit), &
                      &                       sedlay(j,k,issso12))
                 anaerob(k) = posol*umfa !this has P units: kmol P/m3 of pore water
                 seddenit(j) = seddenit(j) + 2._wp*n2prod*posol*umfa/dtbgc*seddw(k)
                 sedlay(j,k,issso12) = sedlay(j,k,issso12) - posol
                 powtra(j,k,ipowaph) = powtra(j,k,ipowaph) + posol*umfa

                 powtra(j,k,ipowno3) = powtra(j,k,ipowno3) - (2._wp*n2prod - rnit)*posol*umfa
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
                 powtra(j,k,ipowh2s) = powtra(j,k,ipowh2s) + posol*umfa*ralk

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

               alk  = powtra(j,k,ipowaal) -( -ansulf(k)  +aerob(k)+anaerob(k))*ralk + nitdem*anaerob(k)
               c    = powtra(j,k,ipowaic) +(anaerob(k) +aerob(k) + ansulf(k))*rcar
                             
               sedhpl(j,k)= update_hi(sedhpl(j,k),c,ak13(j,kbo(j)),&
        &                             ak23(j,kbo(j)),akw3(j,kbo(j)),aks3(j,kbo(j)),&
        &                             akf3(j,kbo(j)),aksi3(j,kbo(j)),ak1p3(j,kbo(j)),&
        &                             ak2p3(j,kbo(j)),ak3p3(j,kbo(j)),psao(j,kbo(j)),&
        &                             akb3(j,kbo(j)),powtra(j,k,ipowasi),powtra(j,k,ipowaph),&
        &                             alk)

                powcar(k)  = c / (1._wp + sedhpl(j,k)/ak13(j,kbo(j)) &
        &                    * (1._wp + sedhpl(j,k)/ak23(j,kbo(j))))
          ENDIF
      END DO


     

     ! Evaluate boundary conditions for sediment-water column exchange.
     ! Current undersaturation of bottom water: sedb(i,0) and
     ! Approximation for new solid sediment, as from sedimentation flux: solrat(1)

     ! CO3 saturation concentration is aksp/calcon as in mo_carchm.f90
     ! (calcon defined in mo_ini_bgc.f90 with 1.03e-2; 1/calcon =~ 97.)

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
                   &                + 2._wp*posol*umfa - ralk*(aerob(k)   &
                   &                - ansulf(k) + anaerob(k)) + nitdem*anaerob(k)

              pown2bud(j,k)       = pown2bud(j,k) +  nitdem*anaerob(k) &
                                  + 2._wp*ralk*ansulf(k)

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
END MODULE mo_powach
