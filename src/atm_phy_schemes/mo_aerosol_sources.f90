!>
!! Module to calculate source terms for different aerosol species.
!!
!! This module contains routines to calculate source terms for different aerosol species.
!! Implemented up to now:
!! - Mineral dust source term based Kok et al. (2014). This aerosol emission is then
!!   converted to an AOD source term for the 2D-aerosol scheme of ICON
!!
!! Literature references:
!! Fecan et al. (1998)   - Fécan, F., Marticorena, B., & Bergametti, G. (1998, December).
!!                         Parametrization of the increase of the aeolian erosion threshold wind
!!                         friction velocity due to soil moisture for arid and semi-arid areas.
!!                         In Annales Geophysicae (Vol. 17, No. 1, pp. 149-157). Springer-Verlag.
!! Kok et al. (2012)     - Kok, J. F., E. J. Parteli, T. I. Michaels, and D. B. Karam, 2012:
!!                         The physics of wind-blown sand and dust. Rep. prog. Phys., 75(10), 106901.
!! Kok et al. (2014)     - Kok, J., N. Mahowald, G. Fratini, J. Gillies, M. Ishizuka, J. Leys, M. Mikami,
!!                         M.-S. Park, S.-U. Park, R. Van Pelt, et al., 2014: An improved dust emission
!!                         model - part 1: Model description and comparison against measurements.
!!                         Atmos. Chem. Phys., 14(23), 13023-13041.
!! Kok et al. (2014b)    - Kok, J., S. Albani, N.M. Mahowald, and D.S. Ward, 2014: An improved dust 
!!                         emission model - Part 2: Evaluation in the Community Earth System Model, 
!!                         with implications for the use of dust source functions
!!                         Atmos. Chem. Phys., 14, 13043-13061
!! Raupach, M. R. (1993) - Dry deposition of gases and particles to vegetation.
!!                         Clean Air: Journal of the Clean Air Society of Australia and New Zealand,
!!                         27(4), 200.
!! Rieger et al. (2017)  - Rieger D., Steiner A., Bachmann V., Gasch P., Foerstner J., Deetz K., 
!!                         Vogel B., and Vogel H., 2017: Impact of the 4 April 2014 Saharan dust 
!!                         outbreak on the photovoltaic power generation in Germany. 
!!                         Atmos. Chem. and Phys. 17 (21), 13391
!! Shao, Y., & Lu, H. (2000) - A simple expression for wind erosion threshold friction velocity.
!!                         Journal of Geophysical Research: Atmospheres, 105(D17), 22437-22443.
!! Zender et al. (2003)  - Zender, C. S., Bian, H., & Newman, D. (2003).
!!                         Mineral Dust Entrainment and Deposition (DEAD) model:
!!                         Description and 1990s dust climatology.
!!                         Journal of Geophysical Research: Atmospheres, 108(D14).
!!
!! @author Daniel Rieger, DWD
!!
!!
!! @par Revision History
!! Initial Revision by Daniel Rieger, DWD (2021-06-29)
!!
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_aerosol_sources

  USE mo_kind,                          ONLY: wp
  USE mo_util_phys,                     ONLY: calc_ustar
  USE mo_physical_constants,            ONLY: grav
  USE mo_aerosol_sources_types,         ONLY: t_dust_source_const

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: aerosol_dust_aod_source

CONTAINS

  !>
  !! SUBROUTINE aerosol_dust_aod_source
  !!
  !! Calculates the source term for mineral dust optical depth based
  !! on external data
  !!
  !! @par Revision History
  !! Initial release by Daniel Rieger, DWD (2019-08-01)
  !!
  SUBROUTINE aerosol_dust_aod_source (this_source, dzsoil, w_so, h_snow, w_so_ice, &
    &                                 soiltyp, plcov, lc_class, rho_a, tcm, u, v, aod_flux, dust_flux)

    TYPE(t_dust_source_const), INTENT(in) :: &
      &  this_source           !< Constant information for dust emission
    REAL(wp), INTENT(in)  :: &
      &  dzsoil,             & !< Thickness of uppermost soil layer (m)
      &  w_so,               & !< Soil water content (m H2O) (liquid+ice)
      &  w_so_ice,           & !< Soil ice content (m H2O)
      &  h_snow,             & !< Snow height (m)
      &  plcov,              & !< Plant cover fraction (-)
      &  rho_a,              & !< Air density in the lowermost layer (kg m-3)
      &  tcm,                & !< Transfer coefficient for momentum
      &  u,                  & !< Meridional wind speed (m s-1)
      &  v                     !< Zonal wind speed (m s-1)
    INTEGER, INTENT(in)   :: &
      &  soiltyp,            & !< Soiltype index (0-9)
      &  lc_class              !< Land use class
    REAL(wp), INTENT(out) :: &
      &  aod_flux              !< Dust optical depth tendency (s-1)
    REAL(wp), OPTIONAL, INTENT(out) :: &
      &  dust_flux             !< Mineral dust emission flux (kg m-2 s-1)
    ! Local Variables
    REAL(wp)              :: &
      &  f_bare,             & !< Bare soil fraction, deducted from several land use classes (-)
      &  f_clay,             & !< Clay fraction (-)
      &  h_snow_fac,         & !< Snow factor, no emission with snow height > 5cm
      &  wgrav,              & !< Gravimetric soil moisture content (%)
      &  f_z0,               & !< Roughness correction factor (-)
      &  f_eta,              & !< Soil moisture correction factor (-)
      &  ustar,              & !< Friction wind speed (m s-1)
      &  ustart                !< Threshold friction wind speed (m s-1)

    ! Initializations
    dust_flux   = 0._wp

    ! Calculate gravimetric soil water content from soil moisture
    wgrav       = calc_wgrav(w_so, w_so_ice , dzsoil)
    ! Use f_clay and f_bare from look up tables
    f_clay      = this_source%f_clay(soiltyp)
    f_bare      = this_source%f_bare(lc_class)
    ! Calculate factors considering the impact of frozen soil and snow height
    h_snow_fac  = calc_hsnow_fac(h_snow)
    ! Calculate soil moisture correction
    f_eta       = calc_aerosol_dust_correta_fecan1998(wgrav, f_clay)
    ! Calculate roughness correction
    f_z0        = calc_aerosol_dust_corrz0_raupach1993(plcov)
    ! Calculate friction wind speed
    ustar       = calc_ustar(tcm, u, v)
    ! Calculate threshold friction wind speed
    ustart      = calc_aerosol_dust_ustart_shaolu2000 (f_eta, f_z0, rho_a)

    ! Calculate the emission flux
    dust_flux   = calc_aerosol_dust_mflux_kok2014 (ustar, ustart, f_bare, f_clay, rho_a)
    dust_flux   = dust_flux * h_snow_fac
    aod_flux    = calc_dust_aod (dust_flux)

  END SUBROUTINE aerosol_dust_aod_source

!-------------------------------------------------------------------------------------------------

  !>
  !! FUNCTION calc_aerosol_dust_correta_fecan1998
  !!
  !! Calculates a soil moisture correction factor for the threshold friction
  !! wind speed according to Fecan et al. (1998).
  !!
  !! @par Revision History
  !! Initial release by Daniel Rieger, DWD (2019-08-01)
  !!
  ELEMENTAL FUNCTION calc_aerosol_dust_correta_fecan1998 (wgrav, f_clay) RESULT (f_eta)

    REAL(wp), INTENT(in)    :: &
      &  wgrav,                & !< Input: Gravimetric soil moisture content       (kg kg-1)
      &  f_clay                  !< Input: Mass fraction of clay particles in soil (-)
    REAL(wp)                :: &
      &  f_eta                   !< Output: Roughness correction factor (-)
    ! Local Variables
    REAL(wp)                :: &
      &  wgravt,               & !< Gravimetric soil moisture content threshold
      &  atune                   !< Ad-hoc factor introduced by Zender et al. (2003)

    ! Range for atune highly model dependant, Kok et al. (2014b) use atune=1
    ! Zender et al. (2003) and ICON-ART (Rieger et al. 2017) use atune = 5
    atune  = 1._wp
    wgravt = atune * (14._wp * f_clay* f_clay + 17._wp * f_clay)

    ! Note that wgrav contains a soil ice content penalty
    f_eta  = SQRT( 1._wp + 1.21_wp * (MAX(0._wp,wgrav - wgravt))**0.68_wp)

  END FUNCTION calc_aerosol_dust_correta_fecan1998

!-------------------------------------------------------------------------------------------------

  !>
  !! FUNCTION calc_aerosol_dust_corrz0_raupach1993
  !!
  !! Calculates a roughness correction factor for the threshold friction
  !! wind speed according to Raupach 1993.
  !!
  !! @par Revision History
  !! Initial release by Daniel Rieger, DWD (2019-08-01)
  !!
  ELEMENTAL FUNCTION calc_aerosol_dust_corrz0_raupach1993 (plcov) RESULT(f_z0)

    REAL(wp), INTENT(in)    :: &
      &  plcov                   !< Input: Fraction of plant cover (-)
    REAL(wp)                :: &
      &  f_z0                    !< Output: Roughness correction factor
    ! Local Variables
    REAL(wp)                :: &
      &  lambda_p                !< plant cover dependant factor

    ! This function requires a check for plcov being smaller than 1 (<0.995)
    lambda_p = - 0.35_wp* LOG(1._wp - plcov)
    f_z0     = SQRT( (1._wp - 0.5_wp * lambda_p) * (1._wp + 45._wp * lambda_p) )

  END FUNCTION calc_aerosol_dust_corrz0_raupach1993

!-------------------------------------------------------------------------------------------------

  !>
  !! FUNCTION calc_aerosol_dust_ustart_shaolu2000
  !!
  !! Calculates the threshold friction wind speed according to Shao and Lu (2000):
  !! Only for friction wind speeds above this threshold, saltation occurs
  !!
  !! @par Revision History
  !! Initial release by Daniel Rieger, DWD (2019-08-01)
  !!
  ELEMENTAL FUNCTION calc_aerosol_dust_ustart_shaolu2000 (f_eta, f_z0, rho_a) RESULT(ustart)

    REAL(wp), INTENT(in)    :: &
      &  f_eta,                & !< Input: Soil moisture correction factor
      &  f_z0,                 & !< Input: Roughness correction factor
      &  rho_a                   !< Input: Air density (kg m-3)
    REAL(wp)                :: &
      &  ustart                  !< Output: Threshold friction wind speed (m s-1)
    ! Local Variables
    REAL(wp), PARAMETER     :: &
     &  D_0   = 75.e-6_wp,    & !< Saltation occurs for u_* > u_*t(D_0) (m)
     &  A_N   = 0.0123_wp,    & !< Coefficient (-)       (Shao and Lu, 2000)
     &  gamma = 3.e-4_wp,     & !< Coefficient (kg s-2)  (Shao and Lu, 2000)
     &  rho_s = 2650_wp,      & !< Bulk soil density (kg m-3)
     &  constval = A_N * ( rho_s * grav * D_0 + gamma / D_0)

    ! Use global minimum of threshold friction velocity d_min as in Rieger et al. (2017)
    !d_min = SQRT(gamma / rho_s / grav)
    ! Calculate threshold friction velocity after Shao and Lu (2000), eq. 24
    !ustart = f_eta * f_z0 * SQRT( A_N / rho_a * ( rho_s * grav * d_min + gamma / d_min) )
    ustart = f_eta * f_z0 * SQRT( constval / rho_a )

  END FUNCTION calc_aerosol_dust_ustart_shaolu2000

!-------------------------------------------------------------------------------------------------

  !>
  !! FUNCTION calc_aerosol_dust_mflux_kok2014
  !!
  !! Calculates the mineral dust emission flux according to Kok et al. (2014). 
  !!
  !! @par Revision History
  !! Initial release by Daniel Rieger, DWD (2019-07-29)
  !!
  ELEMENTAL FUNCTION calc_aerosol_dust_mflux_kok2014 (ustar, ustart, f_bare, f_clay, rho_a) RESULT(dust_flux)

    REAL(wp), INTENT(in)  :: &
      &  ustar,              & !< Input: Friction wind speed (m s-1)
      &  ustart,             & !< Input: Threshold friction wind speed (m s-1) 
      &  f_bare,             & !< Input: Bare soil fraction (-)
      &  f_clay,             & !< Input: Clay fraction (-)
      &  rho_a                 !< Input: Air density (kg m-3)
    REAL(wp)              :: &
      &  dust_flux             !< Output: Mineral dust emission flux (kg m-2 s-1)
    ! Local Variables
    REAL(wp), PARAMETER     :: &
      &  rho_a0   = 1.225_wp,  & !< Standard atmospheric density at sea level (kg m-3)
      &  C_d0     = 3.9e-5_wp, & !< Dimensionless coeff. (4.4 +/- 0.5)*10^(-5) (Kok et al. 2014)
      &  C_e      = 2.3_wp,    & !< Dimensionless coeff.  2.0 +/- 0.3          (Kok et al. 2014)
      &  C_a      = 1.7_wp,    & !< Dimensionless coeff.  2.7 +/- 1.0          (Kok et al. 2014)
      &  ustarst0 = 0.16_wp      !< From measurements (m s-1)                  (Kok et al. 2012)
    REAL(wp)              :: &
      &  C_d,                & !< Dimensionless coefficient
      &  frac_ust_ust0,      & !< (ustarst - ustarst0) / ustarst0
      &  ustarst               !< Standardized threshold friction wind speed (m s-1)

    ! Calculate standardized threshold friction wind speed as defined in Kok et al. (2014), eq. 6.
    ustarst = ustart * SQRT(rho_a / rho_a0)

    ! Calculate the following fraction only once
    frac_ust_ust0 = (ustarst-ustarst0) / (ustarst0)

    ! Calculate C_d according to Kok et al. (2014), eq. 18b
    C_d = C_d0 * exp( -C_e * frac_ust_ust0 )

    ! Calculate dust emission flux according to Kok et al. (2014), eq. 18a
    dust_flux = C_d * f_bare * f_clay * rho_a * (ustar**2 - ustart**2) / ustarst  &
      &        * ( ustar / ustart )**(C_a * frac_ust_ust0)

    dust_flux = MAX(0._wp, dust_flux)

  END FUNCTION calc_aerosol_dust_mflux_kok2014

!-------------------------------------------------------------------------------------------------

  !>
  !! SUBROUTINE calc_dust_aod
  !!
  !! Calculates the mineral dust optical depth.
  !! There is no information on the size distribution required at this point as the size distribution
  !! of emitted dust is assumed to be constant. Hence, the weighted contribution of differently sized
  !! particles can be considered already during the preprocessing step when deriving the factor k_etot.
  !!
  !! Here, a complex refractive index of m=1.55+0.0025j was used.
  !!
  !! @par Revision History
  !! Initial release by Daniel Rieger, DWD (2020-01-09)
  !!
  ELEMENTAL FUNCTION calc_dust_aod (dust_flux) RESULT(aod_flux)

    REAL(wp), INTENT(in)  :: &
      &  dust_flux             !< Mineral dust emission flux (kg m-2 s-1)
    REAL(wp)              :: &
      &  aod_flux              !< Output: Dust optical depth tendency (s-1)
    REAL(wp)              :: &
      &  k_etot                !< Sum of size-distribution weighted mass-specific extinction coefficient (m2 g-1)

    ! Value calculated for a refractive index of m=1.55+0.0025j
    ! Meaningful range:
    ! For size distribution up to 2.5 microns k_etot=0.0689_wp
    ! For size distribution up to 20  microns k_etot=0.2472_wp
    k_etot = 0.1_wp
    aod_flux = k_etot * dust_flux * 1.e+3_wp

  END FUNCTION calc_dust_aod

  !>
  !! FUNCTION calc_wgrav
  !!
  !! Calculates the gravimetric soil moisture content in percent
  !!
  !! @par Revision History
  !! Initial release by Daniel Rieger, DWD (2019-08-06)
  !!
  FUNCTION calc_wgrav(w_so, w_so_ice, dzsoil) RESULT (wgrav)

    REAL(wp), INTENT(in)  :: &
      &  w_so,               & !< Soil water+ice content (m H2O)
      &  w_so_ice,           & !< Soil ice content (m H2O)
      &  dzsoil                !< Thickness of uppermost soil layer (m)
    REAL(wp)              :: &
      &  wgrav                 !< Gravimetric soil moisture content (%)
!    REAL(wp)                 :: &
!      &  rho_s_bulk = 1.5e3_wp, & !< Bulk soil density (kg m-3)
!      &  rho_w= 1000._wp          !< Water density (kg m-3)

    !wgrav=100._wp*(w_so/dzsoil)*rho_s_bulk/rho_w
    ! Give w_so_ice a factor 10 penalty to efficiently inhibit
    ! mineral dust emission on frozen soil
    wgrav=150_wp*( (w_so + 9._wp * w_so_ice )/dzsoil)

  END FUNCTION calc_wgrav

  !>
  !! FUNCTION calc_hsnow_fac
  !!
  !! Calculate a factor to consider the impact of snow height. For a snow height 
  !! above 5 cm, the factor is 0.
  !! Use a simple linear fit in between f(h_snow=0) = 1 and f(h_snow=0.05) = 0
  !!
  !! @par Revision History
  !! Initial release by Daniel Rieger, DWD (2019-08-12)
  !!
  ELEMENTAL FUNCTION calc_hsnow_fac(h_snow) RESULT (h_snow_fac)

    REAL(wp), INTENT(in)  :: &
      &  h_snow                !< Snow height (m)
    REAL(wp)              :: &
      &  h_snow_fac            !< Snow height factor

    h_snow_fac = -20._wp * h_snow + 1._wp

    h_snow_fac = MAX(0._wp,h_snow_fac)

  END FUNCTION calc_hsnow_fac
  
END MODULE mo_aerosol_sources
