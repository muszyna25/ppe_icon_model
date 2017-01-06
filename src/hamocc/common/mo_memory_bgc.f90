!>
!! @brief 
!!
!! Definition of variables and allocation of memory
!!
MODULE mo_memory_bgc

  USE mo_kind, ONLY   : wp
  USE mo_param1_bgc, ONLY: n_bgctra, natm, npowtra, nbgctend,nbgcflux
  USE mo_control_bgc, ONLY: rmasko, bgc_nproma, bgc_zlevs
  USE mo_hamocc_nml, ONLY: l_cpl_co2 
 

  IMPLICIT NONE

  PUBLIC

  INTEGER :: n90depth,n1000depth,n2000depth
  INTEGER, DIMENSION (:), ALLOCATABLE  :: kbo   !< k-index of bottom layer (2d)

  REAL(wp), ALLOCATABLE, TARGET :: bgctra(:,:,:)
  REAL(wp), ALLOCATABLE, TARGET :: bgctend(:,:,:)
  REAL(wp), ALLOCATABLE, TARGET :: bgcflux(:,:)
  REAL(wp), ALLOCATABLE         :: solco2(:)


  REAL(wp), ALLOCATABLE, TARGET :: co3 (:,:)
  REAL(wp), ALLOCATABLE, TARGET :: hi  (:,:)
  REAL(wp), ALLOCATABLE, TARGET :: aksp(:,:)

  REAL(wp), ALLOCATABLE :: swr_frac (:,:)
  REAL(wp), ALLOCATABLE :: meanswr(:,:)
  REAL(wp), ALLOCATABLE :: akw3(:,:)
  REAL(wp), ALLOCATABLE :: akb3(:,:)
  REAL(wp), ALLOCATABLE :: ak13(:,:)
  REAL(wp), ALLOCATABLE :: ak23(:,:)
  REAL(wp), ALLOCATABLE :: aks3(:,:)
  REAL(wp), ALLOCATABLE :: akf3(:,:)
  REAL(wp), ALLOCATABLE :: ak1p3(:,:)
  REAL(wp), ALLOCATABLE :: ak2p3(:,:)
  REAL(wp), ALLOCATABLE :: ak3p3(:,:)
  REAL(wp), ALLOCATABLE :: aksi3(:,:)
  REAL(wp), ALLOCATABLE :: satoxy(:,:)
  REAL(wp), ALLOCATABLE :: satn2(:)
  REAL(wp), ALLOCATABLE :: satn2o(:)
  REAL(wp), ALLOCATABLE :: atdifv(:)
  REAL(wp), ALLOCATABLE :: suppco2(:)
  REAL(wp), ALLOCATABLE :: sedfluxo(:,:)
  REAL(wp), ALLOCATABLE :: aksurf(:,:)    !> dissociation constants for DIC,bor and water 
                                            !> at the sea surface, no pressure dependency 
  REAL(wp), ALLOCATABLE :: atm(:,:)

  REAL(wp), ALLOCATABLE :: co2flux(:)     !> sea-air C-flux, ingetrated over whole simulation period
  REAL(wp), ALLOCATABLE :: o2flux(:)      !> sea-air O2-flux, ingetrated over whole simulation period
  REAL(wp), ALLOCATABLE :: n2flux(:)      !> sea-air N2-flux, ingetrated over whole simulation period
  REAL(wp), ALLOCATABLE :: n2oflux(:)     !> sea-air N2-flux, ingetrated over whole simulation period
  REAL(wp), ALLOCATABLE :: co2trans(:)    !> transfer coefficient for CO2 atmosph/ocean
  REAL(wp), ALLOCATABLE :: wpoc(:)        ! depth-dependent detritus settling speed
  REAL(wp), ALLOCATABLE :: co2conc(:)
  REAL(wp), ALLOCATABLE :: co2flux_cpl(:)

  REAL(wp), DIMENSION (:), ALLOCATABLE :: bolay !<  height of bottom cell
  REAL(wp), DIMENSION (:), ALLOCATABLE :: expoor
  REAL(wp), DIMENSION (:), ALLOCATABLE :: expoca
  REAL(wp), DIMENSION (:), ALLOCATABLE :: exposi
  REAL(wp), DIMENSION (:), ALLOCATABLE :: strahl
  REAL(wp), DIMENSION (:), ALLOCATABLE :: alar1max
  REAL(wp), DIMENSION (:), ALLOCATABLE :: TSFmax
  REAL(wp), DIMENSION (:), ALLOCATABLE :: TMFmax

  REAL(wp), DIMENSION(6):: dmsp

  REAL(wp) :: calcinp, orginp, silinp
  REAL(wp) :: roc13, atcoa
  REAL(wp) :: ozkoa, totalarea
  REAL(wp) :: globalmean_co2, globalmean_o2, globalmean_n2
  REAL(wp) :: atmacon, atmacmol
  REAL(wp) :: ems_per_step
  REAL(wp) :: totarea
  REAL(wp) :: phytomi, grami, grazra, rrrcl
  REAL(wp) :: remido, dyphy, zinges, epsher, spemor, gammap, gammaz, ecan
  REAL(wp) :: ro2ut, rcar, rnit, rnoi, rnit23, rnit13, rcalc, ropal, rn2, p2gtc
  REAL(wp) :: bkphy, bkzoo, bkopal, bifr13, bifr14, plafr13, plafr14
  REAL(wp) :: drempoc, dremdoc, dremcalc, dremn2o
  REAL(wp) :: remido_cya, dremdoc_cya
  REAL(wp) :: thresh_aerob ! O2 threshold for aerob POC remineralization
  REAL(wp) :: thresh_sred ! O2 threshold for sulfate reduction
  REAL(wp) :: thresh_o2, prodn2o 
  REAL(wp) :: mc_fac, mc_eu_speed, mc_depth
  REAL(wp) :: n2_fixation
  REAL(wp) :: sulfate_reduction
  REAL(wp) :: denitrification
  REAL(wp) :: wopal ! daily sinking speed of opal (namelist parameter)
  REAL(wp) :: wcal  ! daily sinking speed of cal (namelist parameter)
  REAL(wp) :: wcya  ! daily boyancy speed of cyanos
  REAL(wp) :: dremopal, calmax, gutc
  REAL(wp) :: psedi, csedi, ssedi
  REAL(wp) :: perc_diron, riron, fesoly, relaxfe, wdust,bolaymin
  REAL(wp) :: pi_alpha, fPAR
  REAL(wp) :: nitdem, n2prod, ro2nitri
  REAL(wp) :: dremn3o, ro2bal
  REAL(wp) :: dustd1, dustd2, dustd3, dustsink

  REAL(wp) :: cycdec, pi_alpha_cya,cya_growth_max,Topt_cya,T1_cya,T2_cya ! (namelist parameter)
  REAL(wp) :: bkcya_P, bkcya_fe, bkcya_N, doccya_fac           ! (namelist parameter)
  REAL(wp) :: buoyancyspeed_cya      



  REAL(wp), PARAMETER :: ppm2con=0.35e-3_wp     !> ppm2con: atmospheric weight: ~10000kg/m^2, avrg. ~29 g/mol
                                                !> --> 350 kmol/m^2 --> 1ppm ~ 0.35e-3 kmol/m^2
  REAL(wp), PARAMETER :: contppm=1._wp/ppm2con  !> inverse of conversion factor molC/m**2 to ppm (in AT)

  REAL(wp), PARAMETER :: molw_co2=44.011_wp     !> Total Molecular Weight of CO2 (44.011 g/mol)
 
  REAL(wp), PARAMETER :: molw_dry_air=28.970_wp !> mean pseudo-molecular weight of dry air



CONTAINS

  SUBROUTINE ALLOC_MEM_CARBCH

    ALLOCATE (wpoc(bgc_zlevs))

    ALLOCATE (bgctra(bgc_nproma,bgc_zlevs,n_bgctra))
    
    ALLOCATE (bgctend(bgc_nproma,bgc_zlevs,nbgctend))
    bgctend = 0._wp

    ALLOCATE (bgcflux(bgc_nproma,nbgcflux))
    bgcflux = 0._wp

    ALLOCATE (hi(bgc_nproma,bgc_zlevs))

    ALLOCATE (co3(bgc_nproma,bgc_zlevs))

    ALLOCATE (solco2(bgc_nproma))
    solco2(:) = rmasko

    ALLOCATE (satn2o(bgc_nproma))
    satn2o(:) = rmasko

    ALLOCATE (satoxy(bgc_nproma,bgc_zlevs))
    satoxy(:,:) = rmasko

    ALLOCATE (satn2(bgc_nproma))
    satn2(:) = rmasko

    ALLOCATE (sedfluxo(bgc_nproma,npowtra))
    sedfluxo(:,:) = 0._wp

    ALLOCATE (aksp(bgc_nproma,bgc_zlevs))
    aksp(:,:) = rmasko

    ALLOCATE (aks3(bgc_nproma,bgc_zlevs))
    aks3(:,:) = rmasko

    ALLOCATE (akf3(bgc_nproma,bgc_zlevs))
    akf3(:,:) = rmasko

    ALLOCATE (ak1p3(bgc_nproma,bgc_zlevs))
    ak1p3(:,:) = rmasko

    ALLOCATE (ak2p3(bgc_nproma,bgc_zlevs))
    ak2p3(:,:) = rmasko

    ALLOCATE (ak3p3(bgc_nproma,bgc_zlevs))
    ak3p3(:,:) = rmasko

    ALLOCATE (aksi3(bgc_nproma,bgc_zlevs))
    aksi3(:,:) = rmasko

    ALLOCATE (ak23(bgc_nproma,bgc_zlevs))

    ALLOCATE (ak13(bgc_nproma,bgc_zlevs))

    ALLOCATE (akb3(bgc_nproma,bgc_zlevs))

    ALLOCATE (akw3(bgc_nproma,bgc_zlevs))
    ALLOCATE (swr_frac(bgc_nproma,bgc_zlevs))
    swr_frac=1._wp
    ALLOCATE (meanswr(bgc_nproma,bgc_zlevs))
    meanswr=0._wp
    ALLOCATE (aksurf(bgc_nproma,10))
     aksurf(:,:) = rmasko
    ALLOCATE (atm(bgc_nproma,natm))
    ALLOCATE (atdifv(bgc_nproma))
    ALLOCATE (suppco2(bgc_nproma))
    suppco2(:) = 0.0_wp

    ALLOCATE (co2flux(bgc_nproma))
    co2flux(:) = 0.0_wp
    ALLOCATE (o2flux(bgc_nproma))
    o2flux(:) = 0.0_wp
    ALLOCATE (n2flux(bgc_nproma))
    n2flux(:) = 0.0_wp
    ALLOCATE (n2oflux(bgc_nproma))
    n2oflux(:) = 0.0_wp

    IF (l_cpl_co2) THEN
      ALLOCATE (co2trans(bgc_nproma))
      co2trans(:) = 0.0_wp
      ALLOCATE (co2conc(bgc_nproma))
      co2conc(:) = 0.0_wp
      ALLOCATE (co2flux_cpl(bgc_nproma))
      co2flux_cpl(:) = 0.0_wp
    ENDIF


  END SUBROUTINE ALLOC_MEM_CARBCH

  SUBROUTINE ALLOC_MEM_BIOMOD

    USE mo_control_bgc
    USE mo_param1_bgc


    ALLOCATE (expoor(bgc_nproma))
    expoor(:)=0._wp
    ALLOCATE (expoca(bgc_nproma))
    expoca(:)=0._wp
    ALLOCATE (exposi(bgc_nproma))
    exposi(:)=0._wp
    ALLOCATE (kbo(bgc_nproma))
    ALLOCATE (strahl(bgc_nproma))
    strahl=200._wp
    ALLOCATE (bolay(bgc_nproma))
    ALLOCATE (alar1max(bgc_nproma))
    ALLOCATE (TSFmax(bgc_nproma))
    ALLOCATE (TMFmax(bgc_nproma))

  END SUBROUTINE ALLOC_MEM_BIOMOD

END MODULE mo_memory_bgc
