MODULE MO_BGC_SURFACE

  USE mo_kind, ONLY           : wp
  USE mo_control_bgc, ONLY    : dtbgc, bgc_nproma, bgc_zlevs

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: gasex, update_weathering, dust_deposition


contains

SUBROUTINE update_weathering ( start_idx,end_idx, pddpo, surface_height)

  USE mo_carbch, ONLY         : bgctra, calcinp, orginp, silinp, bgcflux
  USE mo_param1_bgc, ONLY     : isco212, icalc, ialkali, idoc, isilica, iiron, &
 &                              korginp, ksilinp, kcalinp

  INTEGER, INTENT(in) :: start_idx                  !< 1st REAL of model grid
  INTEGER, INTENT(in) :: end_idx                  !< 2nd REAL of model grid.
  REAL(wp),INTENT(in) :: pddpo(bgc_nproma,bgc_zlevs) !< size of scalar grid cell (3rd REAL) [m]
  REAL(wp),INTENT(in) :: surface_height(bgc_nproma) !< size of scalar grid cell (3rd REAL) [m]
  
  ! Local variables

  INTEGER :: jc

  DO jc = start_idx, end_idx

  if(pddpo(jc,1) > 0.5_wp) then
    bgctra(jc,1,idoc) = bgctra(jc,1,idoc) + orginp / (pddpo(jc,1) + surface_height(jc))
    bgctra(jc,1,isco212) = bgctra(jc,1,isco212) + calcinp / (pddpo(jc,1) + surface_height(jc))
    bgctra(jc,1,ialkali) = bgctra(jc,1,ialkali) + 2._wp * calcinp / (pddpo(jc,1) +surface_height(jc))
    bgctra(jc,1,isilica) = bgctra(jc,1,isilica) + silinp / (pddpo(jc,1) + surface_height(jc))
  
    bgcflux(jc,korginp) = orginp / (pddpo(jc,1) + surface_height(jc))
    bgcflux(jc,ksilinp) = silinp / (pddpo(jc,1) + surface_height(jc))
    bgcflux(jc,kcalinp) = calcinp / (pddpo(jc,1) + surface_height(jc))
  endif
  ENDDO

END SUBROUTINE

SUBROUTINE dust_deposition ( start_idx,end_idx, pddpo,za,dustinp)

  USE mo_carbch, ONLY         : bgctra
  USE mo_param1_bgc, ONLY     : iiron, idust
  USE mo_control_bgc, ONLY    : dtb
  USE mo_biomod,     ONLY     : perc_diron 

  INTEGER, INTENT(in) :: start_idx                  !< 1st REAL of model grid
  INTEGER, INTENT(in) :: end_idx                  !< 2nd REAL of model grid.
  REAL(wp),INTENT(in) :: dustinp(bgc_nproma ) !< size of scalar grid cell (3rd REAL) [m]
  REAL(wp),INTENT(in) :: pddpo(bgc_nproma,bgc_zlevs) !< size of scalar grid cell (3rd REAL) [m]
  REAL(wp),INTENT(in) :: za(bgc_nproma ) !< size of scalar grid cell (3rd REAL) [m]
  
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
!! Computes sea-air gass exchange
!! for oxygen, O2, N2, and CO2.
!! called by bgc
!!
!! @par Revision History
!!
!! First version by Ernst Maier-Reimer  (MPI-M)     Apr 04, 2001
!!

  USE mo_biomod, ONLY         : rrrcl


  USE mo_param1_bgc, ONLY     : icalc, igasnit, ian2o,  iatmco2,          &
       &                        iatmo2, iatmn2, ioxygen, isco212,         &
       &                        ialkali, kcflux, koflux, knflux,          &
       &                        kn2oflux, idms, kdmsflux 

  USE mo_carbch, ONLY         : hi, aksp, akb3, akw3, ak13, ak23, co3,         &
       &                        solco2,satoxy,satn2,aksurf, molw_co2,          &
       &                        globalmean_o2, globalmean_n2, globalmean_co2,  &
       &                        co2flux, co2trans, o2flux,            &
       &                        contppm, satn2o, n2oflux, n2flux,              &
       &                        co2flux_cpl, bgctra, atm, bgcflux

  USE mo_hamocc_nml, ONLY     : l_cpl_co2, l_diffat, atm_co2, atm_o2, atm_n2

  USE mo_bgc_constants, ONLY : cmh2ms

  USE mo_carchm,    ONLY: update_hi

  IMPLICIT NONE

  !! Arguments

  INTEGER, INTENT(in) :: start_idx                !< cell range start  
  INTEGER, INTENT(in) :: end_idx                  !< cell range end

  REAL(wp),INTENT(in) :: pddpo(bgc_nproma,bgc_zlevs) !< size of scalar grid cell (3rd REAL) [m]
  REAL(wp),INTENT(in) :: psao(bgc_nproma,bgc_zlevs)  !< salinity
  REAL(wp),INTENT(in) :: ptho(bgc_nproma,bgc_zlevs)  !< potential temperature
  REAL(wp),INTENT(in) :: pfu10(bgc_nproma)          !< forcing field wind speed
  REAL(wp),INTENT(in) :: psicomo(bgc_nproma)        !< sea ice concentration
  REAL(wp),INTENT(in) :: za(bgc_nproma)             !< sea surface height

  !! Local variables

  INTEGER :: j, k, kpke

  REAL(wp) :: fluxd,fluxu
  REAL(wp) :: kwco2,kwo2, kwdms
  REAL(wp) :: scco2,sco2, scdms
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
           !  Groeger&Mikolajewicz, Ocean Modeling,39,2011.yy
           !  Input is temperature in deg C.
           !  O2 Schmidt number after Keeling et al. (1998, Global Biogeochem.
           !                                                Cycles, 12, 141-163)
           !
           !DMS Schmidt number after Saltzmann et al. (1993, J. Geophys. Res. 98,
           !                                                      16,481-16,486)

           !*********************************************************************

           scco2 = 142.653_wp + 1955.9_wp *  EXP(-0.0663147_wp*ptho(j,1))

           sco2  = 1638.0_wp - 81.83_wp*ptho(j,1) + 1.483_wp*ptho(j,1)**2    &
                &       - 0.008004_wp*ptho(j,1)**3

           scdms = 186.560_wp + 2506.78_wp * EXP(-0.0618603_wp*ptho(j,1))
     
           !
           !  Compute the transfer (piston) velocity in m/s
           !  660 = Schmidt number of CO2 @ 20 degC in seawater

           kwco2 = (1._wp - psicomo( j)) * cmh2ms * pfu10( j)**2        &
                &           * (660._wp / scco2)**0.5_wp


           kwo2  = (1._wp - psicomo( j)) * cmh2ms * pfu10( j)**2        &
                &           * (660._wp / sco2)**0.5_wp

           kwdms = (1._wp - psicomo(j)) * cmh2ms * pfu10( j)**2        &
                &           * (660._wp / scdms)**0.5_wp

           atco2 = atm_co2
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

           nlaughflux = kwo2 * dtbgc * (bgctra(j,1,ian2o)              &
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
    &          aksurf(j,2),  aksurf(j,4),psao(j,k) , aksurf(j,3), bgctra(j,k,ialkali) )


         !
         ! Calculate pCO2 [ppmv] from total dissolved inorganic carbon (DIC: SCO212)
         ! the calculation also includes solubility
         !
           pco2=  bgctra(j,k,isco212)  /((1._wp + aksurf(j,1) * (1._wp + aksurf(j,2)/hi(j,k))/hi(j,k)) * solco2(j))

           fluxd=atm(j,iatmco2)*kwco2*dtbgc*solco2(j) ! *ppao(i,j)/101300.
           fluxu=pco2 *kwco2*dtbgc*solco2(j) ! *ppao(i,j)/101300.

!         ! new concentrations ocean (kmol/m3 -->ppm)
           thickness = pddpo(j,1) + za(j)                             
           bgctra(j,1,isco212) = bgctra(j,1,isco212)+(fluxd-fluxu)/thickness
           bgcflux(j,kcflux) = (fluxu-fluxd)/dtbgc


        ENDIF ! wet cell
     END DO

END SUBROUTINE 
END MODULE MO_BGC_SURFACE
