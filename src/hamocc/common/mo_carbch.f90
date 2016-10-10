!>
!! @brief Variables for inorganic carbon cycle
!!
!! Definition of variables and allocation of memory
!!
!! @author S.Legutke (MPI-M)
!!
!! @par Revision History
!!
!! First version by S.Legutke            (MPI-M)    Oct 31, 2001
!!
!! 2002-2013 by MPI-M
!! This software is provided for non-commercial use only.
!! See the LICENSE and the WARRANTY conditions.
!!
MODULE mo_carbch

  USE mo_kind, ONLY   : wp
  USE mo_param1_bgc, ONLY: n_bgctra, natm, npowtra, nbgctend,nbgcflux
  USE mo_control_bgc, ONLY: rmasko, bgc_nproma, bgc_zlevs
  USE mo_hamocc_nml, ONLY: l_cpl_co2 
 
!  USE mo_commo1, ONLY : l_cpl_co2, l_cpl_dust

  IMPLICIT NONE

  PUBLIC

  REAL(wp), ALLOCATABLE, TARGET :: bgctra(:,:,:)
  REAL(wp), ALLOCATABLE, TARGET :: bgctend(:,:,:)
  REAL(wp), ALLOCATABLE, TARGET :: bgcflux(:,:)
  REAL(wp), ALLOCATABLE         :: solco2(:)

  REAL(wp), ALLOCATABLE :: wpoc(:)  ! depth-dependent detritus settling speed

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

  REAL(wp), ALLOCATABLE :: co2flux(:)     !> sea-air C-flux, ingetrated over whole simulation period
  REAL(wp), ALLOCATABLE :: o2flux(:)      !> sea-air O2-flux, ingetrated over whole simulation period
  REAL(wp), ALLOCATABLE :: n2flux(:)      !> sea-air N2-flux, ingetrated over whole simulation period
  REAL(wp), ALLOCATABLE :: n2oflux(:)     !> sea-air N2-flux, ingetrated over whole simulation period


  REAL(wp), ALLOCATABLE :: co2trans(:)    !> transfer coefficient for CO2 atmosph/ocean
                                            !> exchange, still to be multiplied by v**2
  REAL(wp), ALLOCATABLE :: co2conc(:)
  REAL(wp), ALLOCATABLE :: co2flux_cpl(:)
  REAL(wp), ALLOCATABLE :: atm(:,:)
  REAL(wp) :: calcinp, orginp, silinp


  REAL(wp) :: roc13, atcoa
  REAL(wp) :: ozkoa, totalarea

  !                              .35e-3 * 5.1e14*12. .35e-3 * 5.1e14
  REAL(wp) :: globalmean_co2, globalmean_o2, globalmean_n2
  REAL(wp) :: atmacon, atmacmol

  REAL(wp), PARAMETER :: ppm2con=0.35e-3_wp     !> ppm2con: atmospheric weight: ~10000kg/m^2, avrg. ~29 g/mol
                                                !> --> 350 kmol/m^2 --> 1ppm ~ 0.35e-3 kmol/m^2
  REAL(wp), PARAMETER :: contppm=1._wp/ppm2con  !> inverse of conversion factor molC/m**2 to ppm (in AT)

  REAL(wp), PARAMETER :: molw_co2=44.011_wp     !> Total Molecular Weight of CO2 (44.011 g/mol)
 
  REAL(wp), PARAMETER :: molw_dry_air=28.970_wp !> mean pseudo-molecular weight of dry air

  REAL(wp) :: ems_per_step

  REAL(wp) :: totarea




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
    !ALLOCATE (dusty(bgc_nproma))
    !dusty(:,:) = 0.0_wp


  END SUBROUTINE ALLOC_MEM_CARBCH

  


END MODULE mo_carbch
