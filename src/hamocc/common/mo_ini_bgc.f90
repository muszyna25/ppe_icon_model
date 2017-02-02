#include "hamocc_omp_definitions.inc"
!>
!! @brief set start values for bgc variables
!!
!!
MODULE mo_ini_bgc

  USE mo_kind, ONLY        : wp
  USE mo_memory_bgc, ONLY   : hi, co3, totarea, bgctra, atdifv, atm, &
       &                     atmacon, atmacmol,    &
       &                     atcoa, ozkoa,    &
       &                     wpoc, calcinp,orginp,silinp, &
       &                     phytomi, grami, remido, dyphy, zinges,        &
       &                     epsher, grazra, spemor, gammap, gammaz, ecan, &
       &                     pi_alpha, fpar, bkphy, bkzoo, bkopal,         &
       &                     drempoc, dremdoc,            &
       &                     dremopal, dremn2o, sulfate_reduction,         &
       &                     dremcalc, n2_fixation, ro2ut, rcar, rnit,     &
       &                     rnoi, nitdem, n2prod, rcalc, ropal, calmax,   &
       &                     gutc, perc_diron, riron, fesoly, relaxfe,     &
       &                     denitrification, kbo, bolay, rn2,             &
       &                     dustd1, dustd2, dustsink, wdust, thresh_o2,   &
       &                     cycdec, pi_alpha_cya,cya_growth_max,          &
       &                     Topt_cya,T1_cya,T2_cya,bkcya_N, bkcya_P, bkcya_fe, &
       &                     remido_cya, dremdoc_cya, buoyancyspeed_cya, &
       &                     doccya_fac, thresh_aerob, thresh_sred, &
       &                     wopal, wcal, wcya, p2gtc, ro2bal, dmsp, prodn2o

  USE mo_sedmnt, ONLY      : powtra, sedlay, sedhpl,disso_op,disso_cal,&
       &                     o2ut, rno3, claydens, sred_sed, ks, silsat, &
                             o2thresh 

  USE mo_hamocc_nml, ONLY  : l_diffat, l_cpl_co2, l_cyadyn, l_diffat, i_settling, &
       &                     sinkspeed_poc, sinkspeed_opal, sinkspeed_calc, &
       &                     mc_fac, sinkspeed_martin_ez, mc_depth, denit_sed, disso_po, &
       &                     atm_co2, atm_o2, atm_n2, deltacalc, deltaorg, deltasil     


  USE mo_control_bgc, ONLY : ldtbgc, dtb, dtbgc, rmasko, rmasks, &
       &                     bgc_gin, bgc_arctic, bgc_lab, & 
       &                     bgc_natl, bgc_atl, bgc_tatl, &
       &                     bgc_tropac, &
       &                     bgc_land, bgc_ind, &
       &                     bgc_soce, bgc_npac, bgc_carb, &
       &                     bgc_nproma, bgc_zlevs

  USE mo_param1_bgc, ONLY  : ialkali, ian2o, iatmco2, iatmn2,     &
       &                     iano3, iatmo2, icalc, idet, idoc, igasnit,   &
       &                     iopal, ioxygen, idust,                       &
       &                     iphosph, iphy, ipowaal, ipowaic, ipowaox,    &
       &                     ipowaph, ipowasi, ipown2, ipowno3,           &
       &                     isco212, isilica, isssc12, issso12, issssil, &
       &                     izoo, ipowafe, issster, &
       &                     icya, iiron, idms
!  USE mo_planetary_constants, ONLY: g, rhoref_water
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: ini_aquatic_tracers,            &
       &    ini_pore_water_tracers,         &
       &    ini_atmospheric_concentrations, &
       &    ini_wpoc, bgc_param_conv_unit,  &
       &    ini_continental_carbon_input,   &
       &    set_parameters_bgc!,             &
!       &    level_ini





CONTAINS


  SUBROUTINE SET_PARAMETERS_BGC
    !
    ! Initialize overall time step counter.
    !
    ldtbgc = 0

    IF (l_diffat) THEN
       ! all concentrations will be calculated in carchm
    ELSE
       atm_co2 = 278._wp
       atcoa   = atm_co2  ! emr for use in maschk

       atm_o2  = 196800._wp
       atm_n2  = 802000._wp
    ENDIF

    ozkoa    = 3.e15_wp                      ! ozean kohlenstoff alt, start with dummy value in maschk
    atmacon  = 0.35e-3_wp * 5.1e14_wp*12._wp ! factor for atmospheric concentration of CO2 in g C -> ppm
    atmacmol = 0.35e-3_wp * 5.1e14_wp        ! factor for atmospheric concentration of CO2 in mol C -> ppm
    ! 5.1e14 = total surface of Earth in [m^2]

    !
    ! Biology
    !
    phytomi  = 1.e-11_wp   ! i.e. 1e-5 mmol P/m3 minimum concentration of phytoplankton
    grami    = 1.e-11_wp   ! i.e. 1e-5 mmol P/m3 minimum concentration of zooplankton | test e-11 ->e-8 very slow decay

    remido   = 0.008_wp     ! KS, JS, EMR 12/07/2007
    dyphy    = 0.008_wp     ! 1/d -mortality rate of phytoplankton
    zinges   = 0.6_wp            ! dimensionless fraction -assimilation efficiency of zooplankton
    epsher   = 0.8_wp            ! dimensionless fraction - (1-epsher)=fraction of grazing egested
    grazra   = 1.0_wp       ! 1/d -grazing rate [emr: 0.6-0.9]
    spemor   = 3.e6_wp      ! 1/d -mortality rate of zooplankton
    gammap   = 0.03_wp      ! 1/d -exudation rate
    gammaz   = 0.06_wp      ! 1/d -excretion rate
    ecan     = 0.95_wp           ! fraction of mortality as PO_4
    pi_alpha = 0.02_wp           ! initial slope of production vs irradiance
    fPAR     = 0.4_wp            ! fraction of Photosynthetic Active Radiation
    thresh_aerob = 5.e-8_wp      ! kmol m-3,  O2 threshold for aerob remineralization
    thresh_o2 = 10.e-6_wp      ! kmol m-3,  O2 threshold for aerob remineralization
    thresh_sred = 3.e-6_wp      ! kmol m-3,  O2 threshold for sulfate reduction
    prodn2o = 1.e-4_wp

    calmax = 0.15_wp ! maximum fraction (of "export") for calc production


    bkphy  = 1.e-8_wp 
    bkzoo  = 4.e-8_wp
    bkopal = 1.e-6_wp 


    ! water column remineralisation constants

    drempoc  = 0.025_wp     ! 1/d      ! 0.75/month. H3.1: 0.2/month k=1, 0.1/month k>1
    dremopal = 0.01_wp 

! ------ cyanobacteria
    buoyancyspeed_cya = 1._wp   ! daily buoyancy speed of cya  
    cycdec            = 0.1_wp 
    pi_alpha_cya      = 0.03_wp      ! m2 W-1 d-1
    cya_growth_max    = 0.2_wp      ! d-1
    Topt_cya          = 28._wp       ! deg C
    T1_cya            = 5.5_wp       ! deg C
    T2_cya            = 1._wp        ! deg C
    bkcya_P           = 1.e-8_wp     ! kmol/m3  
    bkcya_fe          = 90.e-8_wp     ! kmol/m3  
    bkcya_N           = 1.e-9_wp     ! kmol/m3  
    doccya_fac        = 0.1_wp
! ------

    dremn2o  = 0.01_wp      ! 1/d
    sulfate_reduction = 0.005_wp
    dremcalc = 0.075_wp     ! 

    ! nitrogen fixation 
    n2_fixation = 0.005_wp

    ! total denitrification rate is a fraction of aerob remineralisation rate drempoc
    denitrification = 0.05_wp   ! 1/d


    ! extended redfield ratio declaration
    ! Note: stoichiometric ratios are based on Takahashi etal. (1985)
    ! P:N:C:-O2 + 1:16:122:172

    ro2ut = 172._wp
    rcar  = 122._wp
    rnit  = 16._wp
    rnoi  = 1._wp/rnit
    rn2   = 2._wp
    ro2bal = ro2ut-rcar-2._wp-1.5_wp*rnit ! OM conversion for O2 mass balance
                                          ! OM { rcar C: rnit N : P }
                                          ! C -> CO2, N-> NO3, P -> PO4
                                          ! C[O2] = 1, N[O2] =1.5, P[O2]=2 
                                          ! total O2   = [O2]+[CO2]+2[PO4]+1.5[NO3] - ro2bal[OM]   

    p2gtc = rcar*12._wp*1.e-12_wp ! kmolP to GtC
    

    ! N consumption of denitrification corrected after Paulmier etal, 2009)
    n2prod = 68.8_wp       ! N2 production for 1 mol P remineralized in suboxic water
    nitdem = 2._wp*n2prod - rnit ! nitrate demand to remin. 1 mol P in suboxic water

    rcalc = 35._wp ! iris 40 !calcium carbonate to organic phosphorous production ratio
    IF (l_cpl_co2) THEN
       rcalc = 20._wp ! emr !calcium carbonate to organic phosphorous production ratio
    END IF

    ropal  = 25._wp ! opal to organic phosphorous production ratio

    ! for interaction with sediment module
    o2ut = 172._wp
    rno3 = 16._wp

    ! Sediment
    denit_sed = 0.01_wp      ! 1/d  sediment denitrification rate 
    sred_sed  = 0.0001_wp    ! sediment sulfate reduction rate
    disso_op  = 0.002592_wp  ! [m3/(kmol Si(OH)4)*1/d ]dissolution rate of opal 
    disso_cal = 0.00864_wp   ! [m3/(kmol CO3) 1/d] dissolution rate of CaCO3
    disso_po  = 0.01_wp      ! [m3/(kmol O2) 1/d] degradation rate organic C in sediment
    silsat    = 0.001_wp     ! [mol/m3] Silicate saturation constant  
    o2thresh  = 10.E-6_wp    ! set to 10 umol/l O2 threshold in sediment


    ! weight percent iron in dust deposition (0.035) times Fe solubility (0.01) /55.85 g--> mol
    perc_diron = 0.035_wp * 0.01_wp / 55.85_wp

    ! the three values below are from Johnson et al., 1997 Mar. Chemistry 57, 137-161
    ! riron   = 5.*rcar*1.e-6       ! 'Q'=5.*122.*1.e-6 = 6.1e-4   (iron to phosphate stoichiometric ratio * umol->mol)
    ! riron   = 2.*rcar*1.e-6       ! 14.2.06: 5.->2. (2.5umol P/nmol Fe) emr: 3.*rcar*1.e-6
    riron   = 3._wp * rcar*1.e-6_wp ! 06.3.06: 2.->3. coex90=12.2GTC/a before change

    fesoly  = 0.6_wp *1.e-9_wp      ! global mean/max (?) dissolved iron concentration in deep water (for relaxation) [mol/m3]
    ! 'fesoly' stands for 'apparent solubility value' of iron

    relaxfe = 0.05_wp/365._wp       ! relaxation time for iron to fesoly corresponds to 20 yrs
    !  relaxfe = 0.005_wp/365._wp*dtb       ! changed as suggested by emr
    ! (relaxfe is only used to remove iron, corresponding to Johnson's 200yr residence time)


    ! DMS
     dmsp(1) = 10._wp     ! temperature dependent release by phytoplankton
     dmsp(2) = 0.0075_wp  ! photolysis (uv-decay)
     dmsp(3) = 0.0096_wp  ! microbial consumption
     dmsp(4) = 1.25_wp * 0.107638502e0_wp  ! production with delcar(coccolithoporides) 
     dmsp(5) = 1.25_wp * 0.109784522e-1_wp  ! production with delsil(diatoms)
     dmsp(6) = 0.1e-07_wp            ! half saturation rate const. bacterial decomp

  END SUBROUTINE SET_PARAMETERS_BGC

  SUBROUTINE BGC_PARAM_CONV_UNIT

   
    remido   = remido * dtb    ! 
    dyphy    = dyphy * dtb    ! 1/d -mortality rate of phytoplankton
    grazra   = grazra * dtb      ! 1/d -grazing rate 
    spemor   = spemor * dtb     ! 1/d -mortality rate of zooplankton
    gammap   = gammap * dtb     ! 1/d -exudation rate
    gammaz   = gammaz * dtb     ! 1/d -excretion rate
    wopal = sinkspeed_opal *dtb
    wcal = sinkspeed_calc * dtb
    wdust = wdust * dtb
    wcya = buoyancyspeed_cya *dtb  !  buoyancy speed of cya  
    sinkspeed_martin_ez = sinkspeed_martin_ez * dtb
    sinkspeed_poc = sinkspeed_poc * dtb
    sulfate_reduction = sulfate_reduction * dtb
    drempoc  = drempoc  *dtb    ! 1/d      
    dremopal = dremopal * dtb  ! 1/d      
    remido_cya  = remido_cya * dtb
    dremn2o  = dremn2o * dtb      ! 1/d
    dremcalc = dremcalc *dtb    ! 
    denit_sed = denit_sed *dtb    ! sediment denitrification rate
    sred_sed = sred_sed *dtb    ! sediment sulfate reduction rate
    relaxfe = relaxfe *dtb       ! relaxation time for iron to fesoly 
    disso_op = disso_op * dtb
    disso_po = disso_po * dtb
    disso_cal = disso_cal * dtb

  END SUBROUTINE

  ! ---------------------------------------------------------------------
  SUBROUTINE ini_wpoc(ptiestw)
  ! initialize wpoc
  ! if lmartin==TRUE (mo_control_bgc, nml)
  ! wpoc increases linearly with depth below mc_depth (beleg, nml)
  ! otherwise the constant sinkspeed_poc (beleg, nml) is used
   REAL(wp),INTENT(in):: ptiestw(bgc_zlevs+1)

   INTEGER :: k 
   REAL(wp) :: at_mc_depth
   
   ! default case: constant sinking speed
   wpoc = sinkspeed_poc 

   IF(i_settling==1)then
   DO k = 1,bgc_zlevs
      ! Are we at z0 (below which the Martin curve starts)?
      at_mc_depth=merge(0._wp,1._wp,ptiestw(k+1)<=mc_depth)
      ! w=w0 + a*(z-z0)
      ! z0= mc_depth
      ! a=remin_rate/b  with F(z)=F(z0)(z/zo)**(-b) 
      wpoc(k) = sinkspeed_martin_ez + at_mc_depth * drempoc/mc_fac * (ptiestw(k+1) - mc_depth) 
   ENDDO
   ENDIF
  END SUBROUTINE ini_wpoc
  
  SUBROUTINE ini_aquatic_tracers (start_idx, end_idx , klevs, ibek )

    INTEGER, INTENT(in)  :: start_idx                  !< 1st REAL of model grid.
    INTEGER, INTENT(in)  :: end_idx                  !< 2nd REAL of model grid.
    INTEGER :: klevs(bgc_nproma)                  !< 3rd (vertical) REAL of model grid.
    INTEGER, INTENT(in)  :: ibek(bgc_nproma)                  !< 3rd (vertical) REAL of model grid.

    INTEGER :: j,k, m, kpke 

    REAL(wp) :: phosat, phosmed, phospac
    REAL(wp) :: oxyat, oxymed, oxypac
    REAL(wp) :: silat, silmed, silpac

    !  Initial values for aquatic (advected) ocean tracers (for new start)

    phosat  = 1.e-6_wp
    phospac = 2.5e-6_wp
    phosmed = 2.e-7_wp
    silat   = 2.e-5_wp
    silmed  = 5.e-6_wp
    silpac  = 1.e-4_wp
    oxypac  = 1.5e-4_wp
    oxyat   = 2.5e-4_wp
    oxymed  = 2.e-4_wp

!HAMOCC_OMP_PARALLEL
!HAMOCC_OMP_DO PRIVATE(j,k,kpke,m) HAMOCC_OMP_DEFAULT_SCHEDULE
   DO j = start_idx, end_idx
    kpke=klevs(j)
    m=ibek(j)
    DO k = 1, kpke

                bgctra(j,k,isco212) = 2.27e-3_wp     ! [kmol/m3]
                bgctra(j,k,ialkali) = 2.37e-3_wp
                bgctra(j,k,iphosph) = phosmed
                bgctra(j,k,ioxygen) = oxymed
                bgctra(j,k,isilica) = silmed
                bgctra(j,k,igasnit)= 1e-10_wp
                bgctra(j,k,idoc)   = 1.e-10_wp
                bgctra(j,k,iphy)   = 1.e-8_wp
                bgctra(j,k,izoo)   = 1.e-8_wp
                bgctra(j,k,idet)   = 1.e-8_wp
                bgctra(j,k,icalc)  = 0._wp
                bgctra(j,k,idms)   = 0._wp
                bgctra(j,k,iopal)  = 1.e-8_wp
                bgctra(j,k,idust)  = 0._wp
                bgctra(j,k,icya)   = 1.e-10_wp
                bgctra(j,k,ian2o)    = 1.e-10_wp
                hi(j,k)              = 3.e-9_wp
                !co3(j,k)             = 0._wp       ! this good for initialisation -> 2.e-4?
                co3(j,k)             = 2.e-4_wp       ! this good for initialisation -> 2.e-4?
                bgctra(j,k,iano3)  = rnit*bgctra(j,k,iphosph)
                bgctra(j,k,iiron)   = 0.6e-9_wp

                if((m.ge.bgc_soce).and.(m.le.bgc_npac))then
                   bgctra(j,k,iphosph) = phospac
                   bgctra(j,k,ioxygen) = oxypac
                   bgctra(j,k,isilica) = silpac
                   bgctra(j,k,iano3)  = rnit*bgctra(j,k,iphosph)
                elseif((m.ge.bgc_gin).and.(m.le.bgc_tatl))THEN
                   bgctra(j,k,iphosph) = phosat
                   bgctra(j,k,ioxygen) = oxyat
                   bgctra(j,k,isilica) = silat
                   bgctra(j,k,iano3)  = rnit*bgctra(j,k,iphosph)
                endif
                if(m.eq.bgc_land)then
                   bgctra(j,k,iphosph) = rmasko
                   bgctra(j,k,isilica) = rmasko
                   bgctra(j,k,ioxygen) = rmasko
                   bgctra(j,k,ialkali) = rmasko
                   bgctra(j,k,iiron)   = rmasko
                   bgctra(j,k,isco212) = rmasko
                   bgctra(j,k,iano3)   = rmasko
                   bgctra(j,k,igasnit) = rmasko
                   bgctra(j,k,idoc)    = rmasko
                   bgctra(j,k,iphy)    = rmasko
                   bgctra(j,k,izoo)    = rmasko
                   bgctra(j,k,idet)    = rmasko
                   bgctra(j,k,icalc)   = rmasko
                   bgctra(j,k,iopal)   = rmasko
                   bgctra(j,k,ian2o)   = rmasko
                   bgctra(j,k,iiron)   = rmasko
                   bgctra(j,k,icya)    = rmasko
                   bgctra(j,k,idms)    = rmasko
                   hi(j,k)             = rmasko
                   co3(j,k)            = rmasko
                ENDIF

       ENDDO
    ENDDO
!HAMOCC_OMP_END_DO
!HAMOCC_OMP_END_PARALLEL

  END SUBROUTINE ini_aquatic_tracers

  ! ---------------------------------------------------------------------

  SUBROUTINE ini_pore_water_tracers(start_idx,end_idx)

    INTEGER, INTENT(in) :: start_idx
    INTEGER, INTENT(in) :: end_idx

    INTEGER :: j, k
    !  Initial values for sediment pore water tracers. (solid components?)

!HAMOCC_OMP_PARALLEL
!HAMOCC_OMP_DO PRIVATE(j,k) HAMOCC_OMP_DEFAULT_SCHEDULE

   DO j = start_idx, end_idx
    DO k = 1, ks
             IF(bolay(j) > 0._wp) THEN
                powtra(j,k,ipowaic) = bgctra(j,kbo(j),isco212)
                powtra(j,k,ipowaal) = bgctra(j,kbo(j),ialkali)
                powtra(j,k,ipowaph) = bgctra(j,kbo(j),iphosph)
                powtra(j,k,ipowaox) = bgctra(j,kbo(j),ioxygen)
                powtra(j,k,ipown2)  = 0._wp
                powtra(j,k,ipowno3) = bgctra(j,kbo(j),iano3)
                powtra(j,k,ipowasi) = bgctra(j,kbo(j),isilica)
                powtra(j,k,ipowafe) = bgctra(j,kbo(j),iiron)
                sedlay(j,k,issso12) = 1.e-8_wp
                sedlay(j,k,isssc12) = 1.e-8_wp
                sedlay(j,k,issster) = 30._wp
                sedlay(j,k,issssil) = 0._wp
                sedhpl(j,k)         = hi(j,kbo(j))
             ELSE
                powtra(j,k,ipowno3) = rmasks   ! pore water
                powtra(j,k,ipown2)  = rmasks
                powtra(j,k,ipowaic) = rmasks
                powtra(j,k,ipowaal) = rmasks
                powtra(j,k,ipowaph) = rmasks
                powtra(j,k,ipowaox) = rmasks
                powtra(j,k,ipowasi) = rmasks
                powtra(j,k,ipowafe) = rmasks
                sedlay(j,k,issso12) = rmasks   ! solid sediment
                sedlay(j,k,isssc12) = rmasks
                sedlay(j,k,issssil) = rmasks
                sedlay(j,k,issster) = rmasks
                sedhpl(j,k)         = rmasks
             ENDIF

       ENDDO
    ENDDO
!HAMOCC_OMP_END_DO
!HAMOCC_OMP_END_PARALLEL

  END SUBROUTINE ini_pore_water_tracers

  ! ---------------------------------------------------------------------

  SUBROUTINE ini_continental_carbon_input(totarea)

    REAL(wp), INTENT(in):: totarea

    calcinp = deltacalc * dtbgc / totarea
    orginp  = deltaorg  * dtbgc / totarea
    silinp  = deltasil  * dtbgc / totarea

  END SUBROUTINE ini_continental_carbon_input

  SUBROUTINE ini_atmospheric_concentrations

    !  atmospheric concentrations (overwritten from restart file)

    atm(:,iatmco2) = atm_co2
    atm(:,iatmo2)  = atm_o2
    atm(:,iatmn2)  = atm_n2

  END SUBROUTINE ini_atmospheric_concentrations

END MODULE 
