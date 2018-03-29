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

  USE mo_psrad_general, ONLY: wp, nbndsw
  USE mo_psrad_radiation_parameters, ONLY: i_overlap
  USE mo_psrad_radiation_parameters, ONLY: l_do_sep_clear_sky, rad_undef
  USE mo_psrad_srtm_kgs, ONLY: ngpt, wavenum2, delwave
  USE mo_psrad_solar_data, ONLY: ssi_default
  USE mo_psrad_srtm_gas_optics, ONLY: gas_optics_sw
  USE mo_psrad_srtm_solver, ONLY: srtm_solver_tr,  srtm_reftra_ec 
  USE mo_psrad_cld_sampling, ONLY: sample_cld_state
#ifdef PSRAD_TIMING
  USE mo_timer, ONLY: ltimer, timer_start, timer_stop, &
    timer_sample_cloud_sw, timer_gas_optics_sw
#endif
#ifdef _OPENMP
  USE omp_lib, ONLY: omp_get_thread_num
#endif
  IMPLICIT NONE

  PRIVATE 
  PUBLIC :: srtm, srtm_diags, srtm_diags_old

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

  SUBROUTINE srtm(kproma, kbdim, klev, &
    albdir_vis_UV, albdif_vis_UV, albdir_NIR, albdif_NIR, &
    cos_mu0, daylght_frc, ssi_factor, local_tsi, cld_frc, &
    tau_cloud, asymm_cloud, omega_cloud, &
    tau_aero_external, asymm_aero_external, omega_aero_external, rnseeds, &
    laytrop, jp, iabs, gases, colmol, fac, &
    h2o_factor, h2o_fraction, h2o_index, &
    flux_dn, flux_up, flux_dn_clr, flux_up_clr, &
    band_weight, per_band_flux)

    USE mo_psrad_general, ONLY: ngptsw, jSFC, jBELOW

    INTEGER, INTENT(IN) :: kproma, kbdim, klev
    REAL(wp), INTENT(IN) :: &
      local_tsi, & ! local solar constant
      cos_mu0(:), & ! Solar zenith angle
      daylght_frc(:),&!< daylight fraction; with diurnal cycle 0 or 1, with zonal mean in [0,1]
      albdir_vis_UV(:), & ! UV/vis direct sfc albedo
      albdir_NIR(:), & ! Near-IR direct sfc albedo
      albdif_vis_UV(:), & ! UV/vis diffuse sfc albedo 
      albdif_NIR(:), & ! Near-IR diffuse sfc albedo
      ssi_factor(:), & ! solar constant factor
      cld_frc(:,:), &
      ! The following are given per band, are spread over correspondingly
      ! named arrays with (:,:,ngptsw) without _per_band names
      tau_cloud(:,:,:), &
      asymm_cloud (:,:,:), &
      omega_cloud(:,:,:), &
      tau_aero_external(:,:,:), &
      asymm_aero_external (:,:,:), &
      omega_aero_external(:,:,:)

    INTEGER, INTENT(IN) :: laytrop(KBDIM), & ! tropopause layer index
      iabs(KBDIM,2,2,klev)
    INTEGER, DIMENSION(:,:), INTENT(IN) :: jp
    REAL(wp), DIMENSION(:,:,:), INTENT(IN) :: h2o_factor,h2o_fraction
    INTEGER, DIMENSION(:,:,:), INTENT(IN) :: h2o_index
    ! column amounts of basic gases: Bjorn, colmol?!
    REAL(wp), DIMENSION(:,:), INTENT(IN) :: colmol
    REAL(wp), INTENT(IN) :: fac(:,:,:,:)
    REAL(wp), TARGET :: gases(:,:,:)


    INTEGER, INTENT(INOUT) :: rnseeds(:,:)

    REAL(WP), INTENT(INOUT) :: &
      flux_dn(KBDIM,klev+1), & ! downward flux total sky
      flux_dn_clr(KBDIM,klev+1), & ! downward flux clear sky
      flux_up(KBDIM,klev+1), & ! upward flux total sky
      flux_up_clr(KBDIM,klev+1) ! upward flux clear sky
    ! All-sky per band fluxes [W/m2], 
    !for computing spectrally-resolved surface fluxes
    ! 1 - upward, 2 - downward, 3 - downward direct
    REAL(WP), INTENT(INOUT) :: per_band_flux(KBDIM,nbndsw,3), &
      band_weight(nbndsw) ! adjustment for current Earth/Sun distance
    
    INTEGER :: i, jk, &  ! level index
      ig, &  ! cumulative g-point index
      band, gpt_in

    REAL(wp), DIMENSION(KBDIM) :: &
      fixed_cos_mu0, & ! Cosine of solar zenith angle
      inv_cos_mu0, &
      albdir, & ! surface albedo, direct          
      albdif ! surface albedo, diffuse         

    ! Cloud, aerosol, and total optical properties per sample
    REAL(wp), DIMENSION(KBDIM,klev) :: &
      tau, & ! total optical depth
      asymm, & ! total asymmetry parameter 
      omega, & ! total single scattering albedo
      pre_asymm_base, &
      pre_omega_base, &
      pre_asymm_cloud, &
      pre_omega_cloud

    ! When computing separate clear-sky fluxes, variables need to 
    ! compute new optical properties for initally cloudy cells
    REAL(wp), DIMENSION(KBDIM,klev+1) :: &
      Rc, Rd, & ! Direct and diffuse reflectance, transmittance
      Tc, Td, & ! from two-stream calculation
      Tb
    ! Narrow band (g-point, or gp) "gptflux_"es in [w/m2]
    REAL(wp), DIMENSION(KBDIM,klev+1) :: &
      gptflux_up, &
      gptflux_dn
    ! "directflux_"es in [w/m2]
    REAL(wp), DIMENSION(KBDIM) :: &
      directflux_dn
    ! parameters from gas optics
    REAL(wp), DIMENSION(KBDIM,klev,ngptsw) :: tau_gas, tau_aer_internal
    REAL(wp) :: solar_flux(KBDIM,ngptsw)
    REAL(wp) :: incident_flux(KBDIM)
    LOGICAL :: cell_cloudy(KBDIM, klev, ngptsw), & ! cloud mask in each cell
      column_cloudy(KBDIM, ngptsw)
    LOGICAL :: true_mask(KBDIM,klev)
    REAL(wp) :: tmp

    INTEGER :: thread_id
#ifdef _OPENMP
    thread_id = omp_get_thread_num() + 1
#else
    thread_id = 1
#endif
    true_mask = .true.

    ! --- weight radiation within a band for the solar cycle ---
    ! local_tsi contains TSI (the "solar constant") scaled with the
    ! Sun-Earth distance. ssi_factor contains the relative contribution
    ! of each band to TSI. ssi_default is the (originally only
    ! implicitly defined) solar flux in the 14 bands.
    !
    band_weight(:) = local_tsi*ssi_factor(:) / ssi_default(:)
    !FUNNY: Lookup roundabout adjflux variable in original! :D

    ! Cloud optical depth is only saved for the band associated with 
    ! this g-point We sample clouds first because we may want to adjust 
    ! water vapor based on presence/absence of clouds
#ifdef PSRAD_TIMING
    IF (ltimer) CALL timer_start(timer_sample_cloud_sw)
#endif
    CALL sample_cld_state(kproma, KBDIM, klev, ngptsw, &
      rnseeds(:,:), i_overlap, cld_frc, cell_cloudy)
    IF (l_do_sep_clear_sky) THEN
      column_cloudy = .false.
      DO ig = 1, ngptsw
      DO jk = 1, klev
        WHERE (cell_cloudy(1:kproma,jk,ig))
          column_cloudy(1:kproma,ig) = column_cloudy(1:kproma,ig) .or. &
            cell_cloudy(1:kproma,jk, ig)
        ENDWHERE
      ENDDO
      ENDDO
    ENDIF

#ifdef PSRAD_TIMING
    IF (ltimer) CALL timer_stop(timer_sample_cloud_sw)
#endif

    ! Loop over g-points calculating gas optical properties. 
#ifdef PSRAD_TIMING
    IF (ltimer) CALL timer_start(timer_gas_optics_sw)
#endif

    CALL gas_optics_sw(kproma, KBDIM, klev, jp, fac, iabs, laytrop, &
      gases, h2o_factor, h2o_fraction, h2o_index, colmol, &
      solar_flux, tau_gas, tau_aer_internal)
#ifdef PSRAD_TIMING
    IF (ltimer) CALL timer_stop(timer_gas_optics_sw)
#endif

    ! Compute radiative transfer.
    flux_up(1:kproma,1:klev+1) = 0.0
    flux_dn(1:kproma,1:klev+1) = 0.0
    flux_up_clr(1:kproma,1:klev+1) = 0.0
    flux_dn_clr(1:kproma,1:klev+1) = 0.0
    per_band_flux = 0.0

    ! Solar illumination
    fixed_cos_mu0(1:kproma) = MAX(0.01_wp, cos_mu0(1:kproma)) !EVIL BUG?
    inv_cos_mu0(1:kproma) = 1._wp / fixed_cos_mu0(1:kproma)

    ! Compute fluxes for each set of samples in turn
    ig = 0
    DO band = 1, nbndsw
      pre_omega_base(1:kproma,:) = &
        tau_aero_external(1:kproma,:,band) * omega_aero_external(1:kproma,:,band)
      pre_asymm_base(1:kproma,:) = asymm_aero_external(1:kproma,:,band) * &
        pre_omega_base(1:kproma,:)
      pre_omega_cloud(1:kproma,:) = tau_cloud(1:kproma,:,band) * omega_cloud(1:kproma,:,band)
      pre_asymm_cloud(1:kproma,:) = asymm_cloud(1:kproma,:,band) * pre_omega_cloud(1:kproma,:)
      DO i = 1, kproma
        albdif(i) = albdif_vis_UV(i) * frc_vis(band) + &
          albdif_NIR(i) * (1.0_wp - frc_vis(band))
        albdir(i) = albdir_vis_UV(i) * frc_vis(band) + &
          albdir_NIR(i) * (1.0_wp - frc_vis(band))
      ENDDO
      DO gpt_in = 1, ngpt(band)
        ig = ig + 1
      ! Combine optical properties of aerosols, clouds, &
      ! gases (absorption + Rayleigh scattering)
        DO i = 1, kproma
          incident_flux(i) = band_weight(band) * &
            solar_flux(i,ig) * fixed_cos_mu0(i) * daylght_frc(i)
        ENDDO

        DO jk = 1, klev
        DO i = 1, kproma
          tau(i,jk) = tau_aer_internal(i,jk,ig) + &
            tau_aero_external(i,jk,band) + &
            tau_gas(i,jk,ig)
          tmp = tau_aer_internal(i,jk,ig) + pre_omega_base(i,jk)
          asymm(i,jk) = pre_asymm_base(i,jk) / tmp
          omega(i,jk) = tmp / tau(i,jk)
        ENDDO
        ENDDO

        IF (l_do_sep_clear_sky) THEN
          CALL srtm_reftra_ec(kproma, KBDIM, klev, .true., true_mask, &
            fixed_cos_mu0, inv_cos_mu0, tau, asymm, omega, &
            Rc, Rd, Tc, Td, Tb)
          ! Compute fluxes using clear-sky optical properties
          CALL srtm_solver_tr(kproma, KBDIM, klev, &
            albdif, albdir, &
            Rc, Rd, Tc, Td, Tb, &
            gptflux_dn, gptflux_up, &
            directflux_dn)
          DO jk = 1,klev+1
          DO i = 1, kproma
            flux_up_clr(i,jk) = flux_up_clr(i,jk) + &
              incident_flux(i) * gptflux_up(i,jk)
            flux_dn_clr(i,jk) = flux_dn_clr(i,jk) + &
              incident_flux(i) * gptflux_dn(i,jk)
          ENDDO
          ENDDO
        END IF
        DO jk = 1,klev
        DO i = 1, kproma
          IF (cell_cloudy(i,jk,ig)) THEN
            tau(i,jk) = tau_aer_internal(i,jk,ig) + &
              tau_aero_external(i,jk,band) + &
              tau_gas(i,jk,ig) + tau_cloud(i,jk,band)
            tmp = tau_aer_internal(i,jk,ig) + &
            pre_omega_base(i,jk) + pre_omega_cloud(i,jk)
            asymm(i,jk) = (pre_asymm_base(i,jk) + pre_asymm_cloud(i,jk)) /&
              tmp
            omega(i,jk) = tmp / tau(i,jk)
          ENDIF
        ENDDO
        ENDDO

        CALL srtm_reftra_ec(kproma, KBDIM, klev, &
          .not. l_do_sep_clear_sky, cell_cloudy(:,:,ig), &
          fixed_cos_mu0, inv_cos_mu0, tau, asymm, omega, &
          Rc, Rd, Tc, Td, Tb)
        CALL srtm_solver_tr(kproma, KBDIM, klev, &
          albdif, albdir, &
          Rc, Rd, Tc, Td, Tb, &
          gptflux_dn, gptflux_up, directflux_dn)

        DO jk = 1,klev+1
        DO i = 1, kproma
          flux_up(i,jk) = flux_up(i,jk) + &
            incident_flux(i) * gptflux_up(i,jk)
          flux_dn(i,jk) = flux_dn(i,jk) + &
            incident_flux(i) * gptflux_dn(i,jk)
        ENDDO
        ENDDO
        DO i = 1, kproma
          per_band_flux(i,band,1) = per_band_flux(i,band,1) + &
            incident_flux(i) * gptflux_up(i,jSFC+jBELOW)
          per_band_flux(i,band,2) = per_band_flux(i,band,2) + &
            incident_flux(i) * gptflux_dn(i,jSFC+jBELOW)
          per_band_flux(i,band,3) = per_band_flux(i,band,3) + &
            incident_flux(i) * directflux_dn(i)
        ENDDO

        IF(.not. l_do_sep_clear_sky) THEN
          ! Accumulate broadband (flx) clear-sky fluxes but here we exclude 
          ! cloudy subcolumns and weight to account for smaller sample size
          DO jk = 1, klev+1
          DO i = 1, kproma
            IF (.not. column_cloudy(i,ig)) THEN
              flux_up_clr(i,jk) = flux_up_clr(i,jk) + &
                  incident_flux(i) * gptflux_up(i,jk)
              flux_dn_clr(i,jk) = flux_dn_clr(i,jk) + &
                incident_flux(i) * gptflux_dn(i,jk)
            ENDIF
          END DO
          END DO
        END IF 
      END DO
    END DO

    ! If computing clear-sky fluxes from samples, flag any columns where 
    ! all samples were cloudy
    IF(.not. l_do_sep_clear_sky) THEN 
      DO i = 1, kproma
        IF(ALL(column_cloudy(i,:))) THEN
          flux_up_clr(i,1:klev+1) = rad_undef
          flux_dn_clr(i,1:klev+1) = rad_undef
        END IF
      END DO
    END IF
  END SUBROUTINE srtm

  ! Derived calculations - diagnostics and heating rates
  ! Spectrally resolved fluxes of various kinds
  SUBROUTINE srtm_diags(kproma, kbdim, per_band_flux, &
    vis_dn_dir_sfc, par_dn_dir_sfc, nir_dn_dir_sfc, vis_dn_dff_sfc, &
    par_dn_dff_sfc, nir_dn_dff_sfc, vis_up_sfc, par_up_sfc, nir_up_sfc)

    INTEGER, INTENT(IN) :: kproma, kbdim

    ! All-sky per band fluxes [W/m2], 
    !for computing spectrally-resolved surface fluxes
    ! 1 - upward, 2 - downward, 3 - downward direct
    REAL(WP), DIMENSION(KBDIM,nbndsw,3), INTENT(IN) :: per_band_flux

    REAL(WP), DIMENSION(KBDIM), INTENT(INOUT) :: &
      ! Fluxes at surface: 
      ! vis_* => Visible (250-680) fraction of net surface radiation
      ! par_* => Photosynthetically Active Radiation 
      ! nir_* => Near infrared radiation
      ! dff/dir - diffuse/direct
      vis_dn_dir_sfc, par_dn_dir_sfc, nir_dn_dir_sfc, &
      vis_dn_dff_sfc, par_dn_dff_sfc, nir_dn_dff_sfc, &
      vis_up_sfc, par_up_sfc, nir_up_sfc

    INTEGER :: band

    vis_dn_dir_sfc = 0.0_wp
    DO band = 1, nbndsw
      vis_dn_dir_sfc(1:kproma) = vis_dn_dir_sfc(1:kproma) + &
        frc_vis(band) * per_band_flux(1:kproma,band,3)
    ENDDO
    par_dn_dir_sfc = 0.0_wp
    DO band = 1, nbndsw
      par_dn_dir_sfc(1:kproma) = par_dn_dir_sfc(1:kproma) + &
        frc_par_array(band) * per_band_flux(1:kproma,band,3)
    ENDDO
    nir_dn_dir_sfc = 0.0_wp
    DO band = 1, nbndsw
      nir_dn_dir_sfc(1:kproma) = nir_dn_dir_sfc(1:kproma) + &
        (1.0_wp - frc_vis(band)) * per_band_flux(1:kproma,band,3)
    ENDDO

    vis_dn_dff_sfc = 0.0_wp
    DO band = 1, nbndsw
      vis_dn_dff_sfc(1:kproma) = vis_dn_dff_sfc(1:kproma) + &
        frc_vis(band) * ( &
        per_band_flux(1:kproma,band,2) - per_band_flux(1:kproma,band,3))
    ENDDO
    par_dn_dff_sfc = 0.0_wp
    DO band = 1, nbndsw
      par_dn_dff_sfc(1:kproma) = par_dn_dff_sfc(1:kproma) + &
        frc_par_array(band) * ( &
        per_band_flux(1:kproma,band,2) - per_band_flux(1:kproma,band,3))
    ENDDO
    nir_dn_dff_sfc = 0.0_wp
    DO band = 1, nbndsw
      nir_dn_dff_sfc(1:kproma) = nir_dn_dff_sfc(1:kproma) + &
        (1.0_wp - frc_vis(band)) * ( &
        per_band_flux(1:kproma,band,2) - per_band_flux(1:kproma,band,3))
    ENDDO


    vis_up_sfc = 0.0_wp
    DO band = 1, nbndsw
      vis_up_sfc(1:kproma) = vis_up_sfc(1:kproma) + &
        frc_vis(band) * per_band_flux(1:kproma,band,1)
    ENDDO
    par_up_sfc = 0.0_wp
    DO band = 1, nbndsw
      par_up_sfc(1:kproma) = par_up_sfc(1:kproma) + &
        frc_par_array(band) * per_band_flux(1:kproma,band,1)
    ENDDO
    nir_up_sfc = 0.0_wp
    DO band = 1, nbndsw
      nir_up_sfc(1:kproma) = nir_up_sfc(1:kproma) + &
        (1.0_wp - frc_vis(band)) * per_band_flux(1:kproma,band,1)
    ENDDO

  END SUBROUTINE srtm_diags

  SUBROUTINE srtm_diags_old(kproma, kbdim, klev, band_weight, &
    per_band_flux, flux_dn, flux_up, vis_frc_sfc, par_dn_sfc, nir_dff_frc, &
    vis_dff_frc, par_dff_frc)

    USE mo_psrad_general, ONLY: jSFC, jBELOW

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: kproma, kbdim, klev

    ! All-sky per band fluxes [W/m2], 
    !for computing spectrally-resolved surface fluxes
    ! 1 - upward, 2 - downward, 3 - downward direct
    REAL(WP), INTENT(IN) :: per_band_flux(KBDIM,nbndsw,3), &
      band_weight(nbndsw)
    REAL(WP), INTENT(IN) :: &
      flux_dn(KBDIM,klev+1), & ! downward flux total sky
      flux_up(KBDIM,klev+1) ! upward flux total sky
    REAL(WP), INTENT(INOUT) :: &
      vis_frc_sfc(KBDIM), & ! Visible (250-680) fraction of 
      !...net surface radiation
      par_dn_sfc(KBDIM), & ! Downward Photosynthetically Active Radiation 
      !...(PAR) at surface
      nir_dff_frc(KBDIM), & ! Diffuse fraction of downward surface 
      !...near-infrared radiation
      vis_dff_frc(KBDIM), & ! Diffuse fraction of downward surface 
      !...visible radiation 
      par_dff_frc(KBDIM)! Diffuse fraction of downward surface PAR

    REAL(wp), DIMENSION(KBDIM,nbndsw) :: zfvis, zfnir, zfpar

    zfvis = SPREAD(band_weight * frc_vis, DIM=1, NCOPIES=kbdim)
    zfnir = SPREAD(band_weight * (1.0_wp - frc_vis), DIM=1, NCOPIES=kbdim)
    zfpar = SPREAD(band_weight * frc_par_array, DIM=1, NCOPIES=kbdim)

    vis_frc_sfc(1:kproma) = SUM(zfvis(1:kproma,:) * &
      (per_band_flux(1:kproma,:,2) - per_band_flux(1:kproma,:,1)), DIM=2) / &
      (flux_dn(1:kproma,jSFC+jBELOW) - flux_up(1:kproma,jSFC+jBELOW) + zepsec)
    vis_frc_sfc(kproma+1:kbdim) = 0

    par_dn_sfc(1:kproma) = &
      SUM(zfpar(1:kproma,:)*(per_band_flux(1:kproma,:,2)), DIM=2) 
    par_dn_sfc(kproma+1:kbdim) = 0

    nir_dff_frc(1:kproma) = SUM(zfnir(1:kproma,:) * &
      (per_band_flux(1:kproma,:,2) - per_band_flux(1:kproma,:,3)), DIM = 2) / & 
      (SUM(zfnir(1:kproma,:) * per_band_flux(1:kproma,:,2), DIM=2) + zepsec)
    nir_dff_frc(kproma+1:kbdim) = 0

    vis_dff_frc(1:kproma) = SUM(zfvis(1:kproma,:) * &
      (per_band_flux(1:kproma,:,2) - per_band_flux(1:kproma,:,3)), DIM = 2) / &
      (SUM(zfvis(1:kproma,:) * per_band_flux(1:kproma,:,2), DIM=2) + zepsec)
    vis_dff_frc(kproma+1:kbdim) = 0

    par_dff_frc(1:kproma) = SUM(zfpar(1:kproma,:) * &
      (per_band_flux(1:kproma,:,2) - per_band_flux(1:kproma,:,3)), DIM = 2) / & 
      (SUM(zfpar(1:kproma,:) * per_band_flux(1:kproma,:,2), DIM=2) + zepsec) 
    par_dff_frc(kproma+1:kbdim) = 0

  END SUBROUTINE srtm_diags_old

END MODULE mo_psrad_srtm_driver
