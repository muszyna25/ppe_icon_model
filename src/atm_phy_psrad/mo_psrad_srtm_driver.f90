!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
!! @brief Module to provide interface to rrtmg shortwave radiation
!!
!! @remarks
!!   This module contains routines that provide the interface between ECHAM
!!   and the AER RRTMG radiation code.  Mostly it organizes and calculates the 
!!   information necessary to call the radiative transfer solvers.
!!
!! @author Bjorn Stevens, MPI-M, Hamburg (2010-08)
!!         Robert Pincus, U. Colorado, visiting MPI-M, Hamburg (2011-07)
!!
!! $ID: n/a$
!!
!! @par Origin
!!   Major segments of this code combines and rewrites (for the ICON standard) 
!!   code previously provided by AER and copyrighed by them.  The authors of the
!!   original AER code are: Eli J. Mlawer, Jennifer S. Delamere, Michael J. 
!!   Iacono and Shepard A. Clough with acknowledgments to Steven J. Taubman, 
!!   Karen Cady-Pereira, Patrick D. Brown, Ronald E. Farren, Luke Chen, Robert 
!!   Bergstrom. The rewrites were designed to better interface with the structure
!!   of the ICON family of models and elements of the ICON programming standard.
!!
!! @par Copyright
!!   The AER copyright
!!   on the original code is as follows: Copyright 2002-2009, Atmospheric and
!!   Environmental Research, Inc. (AER). This software may be used, copied, or
!!   redistributed as long as it is not sold and this copyright notice is
!!   reproduced on each copy made.  This model is provided as is without any
!!   express or implied warranties. (http://www.rtweb.aer.com/)               
!! 
!
MODULE mo_psrad_srtm_driver

  USE mo_psrad_general, ONLY : wp, nbndsw, ngptsw, ngas
  USE mo_psrad_radiation_parameters, ONLY : i_overlap, l_do_sep_clear_sky, &
    rad_undef
  USE mo_psrad_srtm_setup, ONLY : ngb, wavenum2, ssi_default, delwave
  USE mo_psrad_srtm_gas_optics, ONLY : gpt_taumol
  USE mo_psrad_rrtm_coeffs, ONLY : srtm_coeffs
  USE mo_psrad_srtm_solver, ONLY : srtm_solver, delta_scale, two_stream
  USE mo_psrad_cld_sampling, ONLY : sample_cld_state

  IMPLICIT NONE

  PRIVATE 
  PUBLIC :: srtm

  REAL (wp), PARAMETER :: &
    nir_vis_boundary   = 14500._wp, &
    par_lower_boundary = 14285.7143_wp, & ! equivalent to 700nm wavelength
    par_upper_boundary = 25000._wp, & ! equivalent to 400nm wavelength
    zepsec = 1.e-06_wp ! epsilon
  INTEGER :: i
  REAL (wp), PARAMETER :: &
    frc_par_array(nbndsw) = (/ (0.0_wp, i = 1, 8), 0.533725_wp, 1.0_wp, &
      0.550164_wp, (0.0_wp, i = 12, nbndsw) /), & 
    frc_vis(1:nbndsw) = MAX(0.0_wp, MIN(1.0_wp, & 
      (/ ((wavenum2(i) - nir_vis_boundary) / &
      delwave(i), i = 1, nbndsw) /) ))

CONTAINS

  SUBROUTINE srtm(kproma, kbdim, klev, play, tlay, wkl, coldry, &
    asdir, asdif, aldir, aldif, prmu0, daylght_frc, ssi_factor, psctm, &
    cld_frc, cld_tau_sw, cld_cg_sw, cld_piz_sw, aer_tau_sw, aer_cg_sw, &
    aer_piz_sw, rnseeds, flxd_sw, flxu_sw, flxd_sw_clr, flxu_sw_clr, &
    vis_dn_dir_sfc, par_dn_dir_sfc, nir_dn_dir_sfc, vis_dn_dff_sfc, &
    par_dn_dff_sfc, nir_dn_dff_sfc, vis_up_sfc, par_up_sfc, nir_up_sfc)

    INTEGER, INTENT(in) :: kproma, & ! Number of horizontal columns 
    kbdim, & ! Maximum number of columns as declared in calling (sub)prog.
    klev ! Number of model layers

    REAL(wp), INTENT(in) :: &
      psctm, & ! local solar constant
      play(KBDIM,klev), & ! Layer pressures [hPa, mb]
      tlay(KBDIM,klev), & ! Layer temperatures [K]
      prmu0(KBDIM), & ! Solar zenith angle
      daylght_frc(kbdim),&!< daylight fraction; with diurnal cycle 0 or 1, with zonal mean in [0,1]
      wkl(KBDIM,klev,ngas), & ! Gas volume mixing ratios [mol/frac]
      coldry(KBDIM,klev), &  ! Column dry amount
      asdir(KBDIM), & ! UV/vis direct sfc albedo
      aldir(KBDIM), & ! Near-IR direct sfc albedo
      asdif(KBDIM), & ! UV/vis diffuse sfc albedo 
      aldif(KBDIM), & ! Near-IR diffuse sfc albedo
      ssi_factor(nbndsw), & ! solar constant factor
      cld_frc(KBDIM,klev), &
      cld_tau_sw(KBDIM,klev,nbndsw), &
      cld_cg_sw (KBDIM,klev,nbndsw), &
      cld_piz_sw(KBDIM,klev,nbndsw), &
      aer_tau_sw(KBDIM,klev,nbndsw), &
      aer_cg_sw (KBDIM,klev,nbndsw), &
      aer_piz_sw(KBDIM,klev,nbndsw)

    INTEGER, INTENT(INOUT) :: &
      rnseeds(:,:) ! Seeds for random number generator (KBDIM, :) 

    REAL(WP), INTENT(OUT) :: &
      flxd_sw(KBDIM,klev+1), & ! downward flux total sky
      flxd_sw_clr(KBDIM,klev+1), & ! downward flux clear sky
      flxu_sw(KBDIM,klev+1), & ! upward flux total sky
      flxu_sw_clr(KBDIM,klev+1) ! upward flux clear sky
    REAL(WP), DIMENSION(KBDIM), INTENT(OUT) :: &
      ! Fluxes at surface: 
      ! vis_* => Visible (250-680) fraction of net surface radiation
      ! par_* => Photosynthetically Active Radiation 
      ! nir_* => Near infrared radiation
      ! dff/dir - diffuse/direct
      vis_dn_dir_sfc, par_dn_dir_sfc, nir_dn_dir_sfc, &
      vis_dn_dff_sfc, par_dn_dff_sfc, nir_dn_dff_sfc, &
      vis_up_sfc, par_up_sfc, nir_up_sfc

    ! BUG: PACK wrongly used (below) - but variables are unused...
    !LOGICAL :: sunUp(KBDIM) ! Mask for sunlit points
    !INTEGER :: idxSunUp(KBDIM), & ! Indicies of sunlit points
    INTEGER :: jl, &   ! column loop index
      jk, &  ! level index
      ib, &  ! band loop index (starting at 1)
      jb, &  ! band loop index (starting at first solar band)
      ig  ! cumulative g-point index

    REAL(wp) :: cossza(KBDIM), & ! Cosine of solar zenith angle
      adjflux(nbndsw), & ! adjustment for current Earth/Sun distance
      albdir(KBDIM,ngptsw), & ! surface albedo, direct          
      albdif(KBDIM,ngptsw) ! surface albedo, diffuse         

    ! Variables for gas optics calculations
    INTEGER :: laytrop(KBDIM), & ! 100 hPa layer index
      iabs(KBDIM,2,2,klev)
    INTEGER, DIMENSION(KBDIM,klev) :: jp, jp1, jt, jt1, indself, indfor

    ! column amounts of basic gases: Bjorn, colmol?!
    REAL(wp), DIMENSION(KBDIM,klev) :: colmol, selffac, selffrac, &
      forfac, forfrac
    REAL(wp) :: fac(KBDIM,2,2,klev), gases(KBDIM,klev,ngas), &
      bnd_wght(nbndsw)


    ! Cloud, aerosol, and total optical properties per sample
    REAL(wp), DIMENSION(KBDIM,klev,ngptsw) :: &
      ztauc, & ! cloud optical depth
      zasyc, & ! cloud asymmetry parameter 
      zomgc, & ! cloud single scattering albedo
      ztaua, & ! total aerosol optical depth
      zasya, & ! total aerosol asymmetry parameter 
      zomga, & ! total aerosol single scattering albedo
      ztaut, & ! total optical depth
      zasyt, & ! total asymmetry parameter 
      zomgt ! total single scattering albedo

    ! When computing separate clear-sky fluxes, variables need to 
    ! compute new optical properties for initally cloudy cells
    REAL(wp), DIMENSION(KBDIM,klev) :: &
      Rdir, Rdif, & ! Direct and diffuse reflectance, transmittance
      Tdir, Tdif ! from two-stream calculation

    ! Band-by-band fluxes, for computing spectrally-resolved surface fluxes
    REAL(wp), DIMENSION(KBDIM,nbndsw) :: &
      zbbfu, & ! All-sky   upward surface sw flux [w/m2]
      zbbfd, & ! All-sky downward surface sw flux [w/m2]
      zbbfddir! All-sky downward surface direct sw flux [w/m2]
    
    ! Narrow band (g-point, or gp) fluxes
    REAL(wp), DIMENSION(KBDIM,klev+1) :: &
      zgpfu, & ! Fullsky upward sw flux [w/m2]
      zgpfd, & ! Fullsky downward sw flux [w/m2]
      zgpcu, & ! Clearsky upward sw flux [w/m2]
      zgpcd ! Clearsky downward sw flux [w/m2]
    REAL(wp), DIMENSION(KBDIM) :: &
      ztdbt, & ! Fullsky downward direct sw flux [w/m2]
      ztdbtc ! Clearsky downward direct sw flux [w/m2]

    ! parameters from gas optics
    REAL(wp), DIMENSION(KBDIM,klev,ngptsw) :: ztaug, ztaur
    REAL(wp) :: zsflxzen(KBDIM,ngptsw)

    REAL(wp) :: zincflx(KBDIM)
    REAL(wp), DIMENSION(KBDIM,nbndsw) :: zfvis, zfnir, zfpar

    LOGICAL :: colcldMask(KBDIM, ngptsw), &
      cldMask(KBDIM, klev, ngptsw) ! cloud mask in each cell

    ! --- weight radiation within a band for the solar cycle ---
    ! psctm contains TSI (the "solar constant") scaled with the
    ! Sun-Earth distance. ssi_factor contains the relative contribution
    ! of each band to TSI. ssi_default is the (originally only
    ! implicitly defined) solar flux in the 14 bands.
    !
    bnd_wght(:) = psctm*ssi_factor(:) / ssi_default(:)
    DO jb = 1,nbndsw
      adjflux(jb) = bnd_wght(jb)
    ENDDO

    ! Which input points are sunlit? 
    !sunUp(1:kproma)    = prmu0(1:kproma) > 0._wp 
    !BUG: BROKEN - length of idxSunUp/result of pack is COUNT(sunUp)
    !idxSunUp(1:kproma) = PACK( (/ (jl, jl = 1, kproma) /), sunUp(1:kproma))

    ! Cloud optical depth is only saved for the band associated with 
    ! this g-point We sample clouds first because we may want to adjust 
    ! water vapor based on presence/absence of clouds
    CALL sample_cld_state(kproma, KBDIM, klev, ngptsw, rnseeds(:,:), &
      i_overlap, cld_frc, cldMask)

!IBM* ASSERT(NODEPS)
    DO ig = 1, ngptsw
      DO jl = 1, kproma
        ztauc(jl,:,ig) = &
          MERGE(cld_tau_sw(jl,1:klev,ngb(ig)), 0._wp, cldMask(jl,:,ig)) 
        zasyc(jl,:,ig) = &
          MERGE(cld_cg_sw (jl,1:klev,ngb(ig)), 0._wp, cldMask(jl,:,ig)) 
        zomgc(jl,:,ig) = &
          MERGE(cld_piz_sw(jl,1:klev,ngb(ig)), 1._wp, cldMask(jl,:,ig))
      END DO
    END DO ! Loop over samples - done with cloud optical depth calculations 


!IBM* ASSERT(NODEPS)
    DO ig = 1, ngptsw
      DO jl = 1, kproma  
        ib = ngb(ig)
        ! Aerosol optical properties
        ztaua(jl,1:klev,ig) = aer_tau_sw(jl,1:klev,ib)
        zasya(jl,1:klev,ig) = aer_cg_sw (jl,1:klev,ib)
        zomga(jl,1:klev,ig) = aer_piz_sw(jl,1:klev,ib)
        albdif(jl,ig) = asdif(jl)*frc_vis(ib) + &
          aldif(jl)*(1.0_wp - frc_vis(ib))
        albdir(jl,ig) = asdir(jl)*frc_vis(ib) + &
          aldir(jl)*(1.0_wp - frc_vis(ib))
      END DO
    END DO


    ! Cloud masks for sorting out clear skies - by cell and by column
    ! These calculations are only needed if l_do_sep_clear_sky is false
    IF(.not. l_do_sep_clear_sky) THEN 
      ! Are any layers cloudy? 
      colcldMask(1:kproma, :) = ANY(cldMask(1:kproma,:,:), DIM=2)
    END IF

    ! Calculate information needed by the radiative transfer routine
    ! that is specific to this atmosphere, especially some of the 
    ! coefficients and indices needed to compute the optical depths
    ! by interpolating data from stored reference atmospheres. 
    ! The coefficients are functions of temperature and pressure and 
    ! remain the same for all g-point samples.
    ! If gas concentrations, temperatures, or pressures vary with sample (ig) 
!thecoefficientsneedtobecalculatedinsidetheloopoversamples
    CALL srtm_coeffs(KBDIM, klev, play, tlay, coldry, wkl, laytrop, &
      jp, jp1, jt, jt1, iabs, gases, colmol, fac, selffac, selffrac, &
      indself, forfac, forfrac, indfor)

    ! Loop over g-points calculating gas optical properties. 
!IBM* ASSERT(NODEPS)
    CALL gpt_taumol(kproma, kbdim, klev, jp, jt, jt1, laytrop, &
      indself, indfor, gases, colmol, fac, selffac, &
      selffrac, forfac, forfrac, zsflxzen, ztaug, ztaur)

    ! Combine optical properties of aerosols, clouds, &
    ! gases (absorption + Rayleigh scattering)
    ztaut(1:kproma,:,:) = ztaur(1:kproma,:,:) + ztaua(1:kproma,:,:) + &
      ztaug(1:kproma,:,:) + ztauc(1:kproma,:,:)
    zomgt(1:kproma,:,:) = ztaur(1:kproma,:,:) * 1.0_wp + &
      ztaua(1:kproma,:,:) * zomga(1:kproma,:,:) + &
      ztauc(1:kproma,:,:) * zomgc(1:kproma,:,:) 
    zasyt(1:kproma,:,:) = (&
      zasya(1:kproma,:,:) * zomga(1:kproma,:,:) * ztaua(1:kproma,:,:) + &
      zasyc(1:kproma,:,:) * zomgc(1:kproma,:,:) *  ztauc(1:kproma,:,:)) / &
      zomgt(1:kproma,:,:)
    zomgt(1:kproma,:,:) = zomgt(1:kproma,:,:) / ztaut(1:kproma,:,:)

    ! Compute radiative transfer.
    flxu_sw(1:kproma,1:klev+1) = 0.0_wp
    flxd_sw(1:kproma,1:klev+1) = 0.0_wp
    flxu_sw_clr(1:kproma,1:klev+1) = 0.0_wp
    flxd_sw_clr(1:kproma,1:klev+1) = 0.0_wp
    zbbfu(:,:) = 0.0_wp
    zbbfd(:,:) = 0.0_wp
    zbbfddir(:,:) = 0.0_wp

    ! Solar illumination
    cossza(1:kproma) = MAX(0.01_wp, prmu0(1:kproma))
    ! The line above is required to compensate for improper input.
    ! Namely, ICON will pass all data columns further down to the
    ! SW radiation routines, regardless of there being sunlight.
    ! The input data is correct, with a negative angle corresponding to
    ! the sun below the horizon. Spuriously computed columns are simply
    ! discarded higher up in the call chain.

    ! Compute fluxes for each set of samples in turn
    DO ig = 1, ngptsw
      zincflx(1:kproma) = adjflux(ngb(ig)) * &
        zsflxzen(1:kproma,ig) * cossza(1:kproma) * daylght_frc(1:kproma)
      !
      ! All (or full) sky calculation 
      IF(l_do_sep_clear_sky) THEN
      ! Delta scale values of tau, g, w0 for solver. This scaling is done 
      ! inside the solver when l_do_sep_clear_sky, so these variables have 
      ! different values and the end of this if block depending on 
      ! l_do_sep_clear_sky (they aren't used again, though) 
        CALL delta_scale(kproma, KBDIM, klev, ztaut(:,:,ig), zasyt(:,:,ig), &
          zomgt(:,:,ig)) 
        CALL two_stream(kproma, KBDIM, klev, cossza, ztaut(:,:,ig), &
          zomgt(:,:,ig), zasyt(:,:,ig), Rdir, Rdif, Tdir, Tdif)
        CALL srtm_solver(kproma, KBDIM, klev, albdif(:,ig), albdir(:,ig), &
          cossza, ztaut(:,:,ig), Rdir, Rdif, Tdir, Tdif, zgpfd, zgpfu, ztdbt)
      ELSE
        ! Compute fluxes directly from optical properties
        CALL srtm_solver(kproma, KBDIM, klev, albdif(:,ig), albdir(:,ig), &
          cossza, ztaut(:,:,ig), zasyt(:,:,ig), zomgt(:,:,ig), zgpfd(:,:), &
          zgpfu, ztdbt)
      END IF
      ! Accumulate broadband (flx) clear-sky fluxes, and band-by-band (zbb) 
      ! clear sky fluxes
      flxu_sw(1:kproma,1:klev+1) = flxu_sw(1:kproma,1:klev+1) + &
        SPREAD(zincflx(1:kproma),DIM=2, NCOPIES=klev+1) * &
        zgpfu(1:kproma,:)
      flxd_sw(1:kproma,1:klev+1) = flxd_sw(1:kproma,1:klev+1) + &
        SPREAD(zincflx(1:kproma),DIM=2, NCOPIES=klev+1) * &
        zgpfd(1:kproma,:)

!IBM* ASSERT(NODEPS)
      DO jl = 1, kproma
        ib = ngb(ig)
        zbbfu(jl,ib) = zbbfu(jl,ib) + zincflx(jl) * &
          zgpfu(jl,klev+1)
        zbbfd(jl,ib) = zbbfd(jl,ib) + zincflx(jl) * &
          zgpfd(jl,klev+1)
        zbbfddir(jl,ib) = zbbfddir(jl,ib) + zincflx(jl) * &
          ztdbt(jl) 
      END DO

      IF(l_do_sep_clear_sky) THEN
        IF(ANY(cldMask(1:kproma,:,ig))) THEN 
          WHERE(cldMask(1:kproma,:,ig))
            ztaut(1:kproma,:,ig) = ztaur(1:kproma,:,ig) + &
              ztaua(1:kproma,:,ig)+ ztaug(1:kproma,:,ig)
            zomgt(1:kproma,:,ig) = ztaur(1:kproma,:,ig) * 1.0_wp + &
              ztaua(1:kproma,:,ig) * zomga(1:kproma,:,ig)
            zasyt(1:kproma,:,ig) =(zasya(1:kproma,:,ig) * &
              zomga(1:kproma,:,ig) * ztaua(1:kproma,:,ig)) / &
              zomgt(1:kproma,:,ig)
            zomgt(1:kproma,:,ig) = zomgt(1:kproma,:,ig) / ztaut(1:kproma,:,ig)
          END WHERE

          CALL delta_scale(kproma, KBDIM, klev, ztaut(:,:,ig), &
            zasyt(:,:,ig), zomgt(:,:,ig), update=cldMask(:,:,ig)) 
          CALL two_stream(kproma, KBDIM, klev, cossza, ztaut(:,:,ig), &
            zomgt(:,:,ig), zasyt(:,:,ig), Rdir, Rdif, Tdir, Tdif, &
            update = cldMask(:,:,ig)) 

          ! Compute fluxes using clear-sky optical properties
          CALL srtm_solver(kproma, KBDIM, klev, albdif(:,ig), albdir(:,ig), &
            cossza, ztaut(:,:,ig), Rdir, Rdif, Tdir, Tdif, zgpcd, zgpcu, &
            ztdbtc)

        ELSE 
          !
          ! Clear-sky and all-sky fluxes are the same
          !
          zgpcu(1:kproma,:) = zgpfu(1:kproma,:)
          zgpcd(1:kproma,:) = zgpfd(1:kproma,:)
        END IF
        !
        ! Accumulate broadband (flx) clear-sky fluxes
        !
        flxu_sw_clr(1:kproma,1:klev+1) = flxu_sw_clr(1:kproma,1:klev+1) + &
          SPREAD(zincflx(1:kproma),DIM=2, NCOPIES=klev+1) * &
          zgpcu(1:kproma,1:klev+1)
        flxd_sw_clr(1:kproma,1:klev+1) = flxd_sw_clr(1:kproma,1:klev+1) + &
          SPREAD(zincflx(1:kproma),DIM=2, NCOPIES=klev+1) * &
          zgpcd(1:kproma,1:klev+1)
      ELSE 
        ! Accumulate broadband (flx) clear-sky fluxes but here we exclude 
        ! cloudy subcolumns and weight to account for smaller sample size
!IBM* ASSERT(NODEPS)
        DO jk = 1, klev+1
          flxu_sw_clr(1:kproma,jk) = flxu_sw_clr(1:kproma,jk) + MERGE(0.0_wp,&
            zincflx(1:kproma) * zgpfu(1:kproma,jk),&
            colCldMask(1:kproma,ig)) 
          flxd_sw_clr(1:kproma,jk) = flxd_sw_clr(1:kproma,jk) + MERGE(0.0_wp,&
            zincflx(1:kproma) * zgpfd(1:kproma,jk),&
            colCldMask(1:kproma,ig)) 
        END DO
      END IF 
    END DO

    ! Derived calculations - diagnostics and heating rates
    ! Spectrally resolved fluxes of various kinds
    zfvis(1:kproma,1:nbndsw) = SPREAD(frc_vis(1:nbndsw), DIM=1, NCOPIES=kproma)
    zfnir(1:kproma,1:nbndsw) = SPREAD(1.0_wp - frc_vis(1:nbndsw), DIM=1, NCOPIES=kproma)
    zfpar(1:kproma,1:nbndsw) = SPREAD(frc_par_array(1:nbndsw), DIM=1, NCOPIES=kproma)

    vis_dn_dir_sfc(1:kproma) = SUM( zfvis(1:kproma,1:nbndsw) * zbbfddir(1:kproma,1:nbndsw), DIM = 2)
    par_dn_dir_sfc(1:kproma) = SUM( zfpar(1:kproma,1:nbndsw) * zbbfddir(1:kproma,1:nbndsw), DIM = 2)
    nir_dn_dir_sfc(1:kproma) = SUM( zfnir(1:kproma,1:nbndsw) * zbbfddir(1:kproma,1:nbndsw), DIM = 2)

    vis_dn_dff_sfc(1:kproma) = SUM( zfvis(1:kproma,1:nbndsw) * (zbbfd(1:kproma,1:nbndsw) - zbbfddir(1:kproma,1:nbndsw)), DIM = 2)
    par_dn_dff_sfc(1:kproma) = SUM( zfpar(1:kproma,1:nbndsw) * (zbbfd(1:kproma,1:nbndsw) - zbbfddir(1:kproma,1:nbndsw)), DIM = 2)
    nir_dn_dff_sfc(1:kproma) = SUM( zfnir(1:kproma,1:nbndsw) * (zbbfd(1:kproma,1:nbndsw) - zbbfddir(1:kproma,1:nbndsw)), DIM = 2)

    vis_up_sfc    (1:kproma) = SUM( zfvis(1:kproma,1:nbndsw) * zbbfu(1:kproma,1:nbndsw), DIM = 2)
    par_up_sfc    (1:kproma) = SUM( zfpar(1:kproma,1:nbndsw) * zbbfu(1:kproma,1:nbndsw), DIM = 2)
    nir_up_sfc    (1:kproma) = SUM( zfnir(1:kproma,1:nbndsw) * zbbfu(1:kproma,1:nbndsw), DIM = 2)


    ! If computing clear-sky fluxes from samples, flag any columns where 
    ! all samples were cloudy
    IF(.not. l_do_sep_clear_sky) THEN 
!IBM* ASSERT(NODEPS)
      DO jl = 1, kproma
        IF(ALL(colCldMask(jl,:))) THEN
          flxu_sw_clr(jl,1:klev+1) = rad_undef
          flxd_sw_clr(jl,1:klev+1) = rad_undef
        END IF
      END DO
    END IF
  END SUBROUTINE srtm

END MODULE mo_psrad_srtm_driver
