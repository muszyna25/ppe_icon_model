MODULE mo_bgc_surface
!! @file mo_bgc_surface.f90
!! @brief module contains gas exchange, weathering fluxes,
!!        dust & nitrogen deposition

  USE mo_kind, ONLY           : wp
  USE mo_control_bgc, ONLY    : dtbgc, bgc_nproma, bgc_zlevs

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: gasex, update_weathering, dust_deposition, nitrogen_deposition


contains

SUBROUTINE update_weathering ( start_idx,end_idx, pddpo, za)
! apply weathering rates

  USE mo_memory_bgc, ONLY     : bgctra, calcinp, orginp, silinp, bgcflux
  USE mo_param1_bgc, ONLY     : isco212, ialkali, idoc, isilica,  &
 &                              korginp, ksilinp, kcalinp

  ! Arguments
 
  INTEGER, INTENT(in)            :: start_idx              !< start index for j loop (ICON cells, MPIOM lat dir)  
  INTEGER, INTENT(in)            :: end_idx                !< end index  for j loop  (ICON cells, MPIOM lat dir) 

  REAL(wp), INTENT(in), TARGET   :: pddpo(bgc_nproma,bgc_zlevs)      !< size of scalar grid cell (3rd dimension) [m]
  REAL(wp), INTENT(in), TARGET   :: za(bgc_nproma)      !< surface height

  ! Local variables

  INTEGER :: jc

  DO jc = start_idx, end_idx

  if(pddpo(jc,1) > 0.5_wp) then

    bgctra(jc,1,idoc) = bgctra(jc,1,idoc) + orginp / (pddpo(jc,1) + za(jc))
    bgctra(jc,1,isco212) = bgctra(jc,1,isco212) + calcinp / (pddpo(jc,1) + za(jc))
    bgctra(jc,1,ialkali) = bgctra(jc,1,ialkali) + 2._wp * calcinp / (pddpo(jc,1) + za(jc))
    bgctra(jc,1,isilica) = bgctra(jc,1,isilica) + silinp / (pddpo(jc,1) + za(jc))
  
    bgcflux(jc,korginp) = orginp / (pddpo(jc,1) + za(jc))
    bgcflux(jc,ksilinp) = silinp / (pddpo(jc,1) + za(jc))
    bgcflux(jc,kcalinp) = calcinp / (pddpo(jc,1) + za(jc))

  endif

  ENDDO

END SUBROUTINE

SUBROUTINE nitrogen_deposition ( start_idx,end_idx, pddpo,za,nitinp)
! apply nitrogen deposition
  USE mo_memory_bgc, ONLY     : bgctra, bgctend
  USE mo_param1_bgc, ONLY     : iano3, ialkali, kn2b
  USE mo_bgc_constants, ONLY  : rmnit


  !Arguments

  INTEGER, INTENT(in)            :: start_idx              !< start index for j loop (ICON cells, MPIOM lat dir)  
  INTEGER, INTENT(in)            :: end_idx                !< end index  for j loop  (ICON cells, MPIOM lat dir) 

  REAL(wp),INTENT(in) :: nitinp(bgc_nproma )                         !< nitrogen input
  REAL(wp), INTENT(in), TARGET   :: pddpo(bgc_nproma,bgc_zlevs)      !< size of scalar grid cell (3rd dimension) [m]
  REAL(wp), INTENT(in), TARGET   :: za(bgc_nproma)                   !< surface height
  
  ! Local variables

  INTEGER :: jc
  REAL(wp) :: ninp 

  DO jc = start_idx, end_idx

  if(pddpo(jc,1) > 0.5_wp) then

      ! ndepo : CCMI wet+dry dep of NHx and NOy in kg (N) m-2 s-1
       ninp = nitinp(jc) / rmnit* dtbgc/(pddpo(jc,1)+za(jc)) ! kmol N m-3 time_step-1

       bgctra(jc,1,iano3) = bgctra(jc,1,iano3) + ninp
       bgctra(jc,1,ialkali) = bgctra(jc,1,ialkali) - ninp
       bgctend(jc,1,kn2b)   = bgctend(jc,1,kn2b) - ninp * (pddpo(jc,1) + za(jc)) 

  endif

  ENDDO


END SUBROUTINE
SUBROUTINE dust_deposition ( start_idx,end_idx, pddpo,za,dustinp)
! apply dust deposition
  USE mo_memory_bgc, ONLY      : bgctra,perc_diron 
  USE mo_param1_bgc, ONLY     : iiron, idust
  USE mo_control_bgc, ONLY    : dtb

  
  !Arguments

  INTEGER, INTENT(in)            :: start_idx              !< start index for j loop (ICON cells, MPIOM lat dir)  
  INTEGER, INTENT(in)            :: end_idx                !< end index  for j loop  (ICON cells, MPIOM lat dir) 

  REAL(wp),INTENT(in) :: dustinp(bgc_nproma )                        !< dust input
  REAL(wp), INTENT(in), TARGET   :: pddpo(bgc_nproma,bgc_zlevs)      !< size of scalar grid cell (3rd dimension) [m]
  REAL(wp), INTENT(in), TARGET   :: za(bgc_nproma)                   !< surface height
  
  ! Local variables

  INTEGER :: jc


  DO jc = start_idx, end_idx

  if(pddpo(jc,1) > 0.5_wp) then

   bgctra(jc,1,iiron) = bgctra(jc,1,iiron) + dustinp(jc)*dtb/365._wp/(pddpo(jc,1)+za(jc)) *perc_diron 

   bgctra(jc,1,idust) = bgctra(jc,1,idust) + dustinp(jc)*dtb/365._wp/(pddpo(jc,1)+za(jc))  

  endif

 ENDDO

END SUBROUTINE


SUBROUTINE gasex ( start_idx,end_idx, pddpo, za, psao, ptho,  &
     &              pfu10, psicomo )
!! @brief Computes sea-air gass exchange
!!         for oxygen, O2, N2, N2O, DMS, and CO2.
!!


  USE mo_param1_bgc, ONLY     : igasnit, ian2o,  iatmco2, iphosph,         &
       &                        ioxygen, isco212, isilica,       &
       &                        ialkali, kcflux, koflux, knflux,          &
       &                        kn2oflux, idms, kdmsflux 

  USE mo_memory_bgc, ONLY         : hi, &
       &                        solco2,satoxy,satn2,aksurf,    &
       &                        satn2o,            &
       &                        bgctra, atm, bgcflux

  USE mo_hamocc_nml, ONLY     : l_cpl_co2, atm_co2, atm_o2, atm_n2

  USE mo_bgc_constants, ONLY : cmh2ms

  USE mo_carchm,    ONLY: update_hi

  IMPLICIT NONE

  !! Arguments

  INTEGER, INTENT(in)  :: start_idx              !< start index for j loop (ICON cells, MPIOM lat dir)  
  INTEGER, INTENT(in)  :: end_idx                !< end index  for j loop  (ICON cells, MPIOM lat dir) 

  REAL(wp),INTENT(in) :: pddpo(bgc_nproma,bgc_zlevs) !< size of scalar grid cell (3rd REAL) [m]
  REAL(wp),INTENT(in) :: psao(bgc_nproma,bgc_zlevs)  !< salinity
  REAL(wp),INTENT(in) :: ptho(bgc_nproma,bgc_zlevs)  !< potential temperature
  REAL(wp),INTENT(in) :: pfu10(bgc_nproma)           !< forcing field wind speed
  REAL(wp),INTENT(in) :: psicomo(bgc_nproma)         !< sea ice concentration
  REAL(wp),INTENT(in) :: za(bgc_nproma)              !< sea surface height

  !! Local variables

  INTEGER :: j, k

  REAL(wp) :: fluxd,fluxu
  REAL(wp) :: kwco2,kwo2, kwdms, kwn2o
  REAL(wp) :: scco2,sco2, scdms, scn2o
  REAL(wp) :: oxflux,niflux,nlaughflux, dmsflux
  REAL(wp) :: ato2, atn2, atco2,pco2
  REAL(wp) :: thickness

  !
  !---------------------------------------------------------------------
  !


  k = 1      ! surface layer

  DO j = start_idx, end_idx


        IF (pddpo(j, 1) .GT. 0.5_wp) THEN

         
           !*********************************************************************
           !
           !  Compute the Schmidt number of CO2 in seawater and the transfer
           !  (piston) velocity using the formulation presented in
           !   Wanninkhof 2014 
           !*********************************************************************

           scco2 = 2116.8_wp - 136.25_wp*ptho(j,1) + 4.7353_wp*ptho(j,1)**2 &
                &   - 0.092307_wp*ptho(j,1)**3 + 0.0007555_wp*ptho(j,1)**4

           sco2 =  1920.4_wp - 135.6_wp*ptho(j,1)  + 5.2122_wp*ptho(j,1)**2 &
               & - 0.10939_wp*ptho(j,1)**3 + 0.00093777_wp*ptho(j,1)**4

           scn2o =  2356.2_wp - 166.38_wp*ptho(j,1)  + 6.3952_wp*ptho(j,1)**2 &
               & - 0.13422_wp*ptho(j,1)**3 + 0.0011506_wp*ptho(j,1)**4

           scdms =  2855.7_wp - 177.63_wp*ptho(j,1)  + 6.0438_wp*ptho(j,1)**2 &
               & - 0.11645_wp*ptho(j,1)**3 + 0.00094743_wp*ptho(j,1)**4
     
           !
           !  Compute the transfer (piston) velocity in m/s
           !  660 = Schmidt number of CO2 @ 20 degC in seawater

           kwco2 = (1._wp - psicomo( j)) * cmh2ms * pfu10( j)**2        &
                &           * (660._wp / scco2)**0.5_wp


           kwo2  = (1._wp - psicomo( j)) * cmh2ms * pfu10( j)**2        &
                &           * (660._wp / sco2)**0.5_wp

           kwdms = (1._wp - psicomo(j)) * cmh2ms * pfu10( j)**2        &
                &           * (660._wp / scdms)**0.5_wp

       
           kwn2o = (1._wp - psicomo(j)) * cmh2ms * pfu10( j)**2        &
                &           * (660._wp / scn2o)**0.5_wp

           if(l_cpl_co2)then 
            atco2 = atm(j,iatmco2)
           else
            atco2 = atm_co2
           endif
           ato2  = atm_o2
           atn2  = atm_n2

          !*********************************************************************
          !     
          ! Calculate air-sea exchange for O2, N2, N2O 
          !
          !*********************************************************************

           ! Surface flux of oxygen
           ! (Meiner-Reimer et. al, 2005, Eq. 74)

           oxflux = kwo2 * dtbgc * (bgctra(j,1,ioxygen)                &
                &  -satoxy(j,1) * (ato2 / 196800._wp)) ! *ppao(i,j)/101300. ! sea level pressure normalization


           bgcflux(j,koflux) = oxflux/dtbgc

           bgctra(j,1,ioxygen) = bgctra(j,1,ioxygen)                 &
                &                - oxflux/(pddpo(j,1)+za(j))

         
           ! Surface flux of gaseous nitrogen (same piston velocity as for O2)
           ! (Meiner-Reimer et. al, 2005, Eq. 75)

           niflux = kwo2 * dtbgc * (bgctra(j,1,igasnit)                &
                & -satn2(j)*(atn2/802000._wp)) ! *ppao(i,j)/101300.

           bgcflux(j,knflux) = niflux/dtbgc

           bgctra(j,1,igasnit) = bgctra(j,1,igasnit)                   &
                &                - niflux/(pddpo(j,1)+za(j))

          
           ! Surface flux of laughing gas (same piston velocity as for O2 and N2)
           ! (Meiner-Reimer et. al, 2005, Eq. 76)

           nlaughflux = kwn2o * dtbgc * (bgctra(j,1,ian2o)              &
                &     - satn2o(j))  

           bgctra(j,1,ian2o) = bgctra(j,1,ian2o)                     &
                &              - nlaughflux/(pddpo(j,1)+za(j))
           bgcflux(j,kn2oflux) = nlaughflux/dtbgc



           ! Surface flux of dms

           ! (Meiner-Reimer et. al, 2005, Eq. 77)

           dmsflux = kwdms*dtbgc*bgctra(j,1,idms)

           bgctra(j,1,idms) = bgctra(j,1,idms) - dmsflux/pddpo(j,1)

           bgcflux(j,kdmsflux) = dmsflux/dtbgc

        !*********************************************************************
        !     
        ! Calculate air sea exchange for CO2
        !
        !*********************************************************************

         ! Update hi

           hi(j,k) = update_hi(hi(j,k), bgctra(j,k,isco212), aksurf(j,1) , &
    &          aksurf(j,2),  aksurf(j,4), aksurf(j,7), aksurf(j,6), aksurf(j,5),&
    &          aksurf(j,8), aksurf(j,9),aksurf(j,10), psao(j,k) , aksurf(j,3), &
    &          bgctra(j,k,isilica), bgctra(j,k,iphosph),bgctra(j,k,ialkali) )


         !
         ! Calculate pCO2 [ppmv] from total dissolved inorganic carbon (DIC: SCO212)
         ! the calculation also includes solubility
         !
           pco2=  bgctra(j,k,isco212)  /((1._wp + aksurf(j,1) * (1._wp + aksurf(j,2)/hi(j,k))/hi(j,k)) * solco2(j))

           fluxd=atco2*kwco2*dtbgc*solco2(j) ! 
           fluxu=pco2 *kwco2*dtbgc*solco2(j) ! 

!         ! new concentrations ocean (kmol/m3 -->ppm)
           thickness = pddpo(j,1) + za(j)                             
           bgctra(j,1,isco212) = bgctra(j,1,isco212)+(fluxd-fluxu)/thickness
           bgcflux(j,kcflux) = (fluxu-fluxd)/dtbgc


        ENDIF ! wet cell
     END DO

END SUBROUTINE 
END MODULE mo_bgc_surface
