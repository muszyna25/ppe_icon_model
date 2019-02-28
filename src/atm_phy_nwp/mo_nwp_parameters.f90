!>
!!  Declares parameters computed during the initialization of the physics
!!  parameterizations that have to be domain-dependent
!!
!! @par Revision History
!!  Guenther Zaengl, DWD, 2011-12-08
!!  - Restructuring the namelists
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_nwp_parameters
!-------------------------------------------------------------------------
  USE mo_kind,               ONLY: wp


  IMPLICIT NONE

  PRIVATE


  TYPE t_phy_params
    !
    ! Parameters which are only computed if convection is switched on
    !
    ! Level parameters for convection scheme
    INTEGER  :: kcon1, kcon2
    ! resolution-dependent parameters for convection scheme
    REAL(wp) :: tau, mfcfl, tau0
    ! relative humidity below which sub-cloud evaporation of rain starts over land and water, respectively
    REAL(wp) :: rhebc_land, rhebc_ocean, rhebc_land_trop, rhebc_ocean_trop
    ! 'excess values' of temperature and QV used for convection triggering (test parcel ascent)
    REAL(wp) :: texc, qexc
    ! fractional area covered by convective precipitation
    REAL(wp) :: rcucov, rcucov_trop
    ! tuning coefficient for organized entrainment of deep convection
    REAL(wp) :: entrorg
    ! coefficient for conversion of cloud water into precipitation
    REAL(wp) :: rprcon
    ! maximum allowed depth of shallow convection (hPa)
    REAL(wp) :: rdepths
    ! switches for activation of shallow, midlevel and deep convection
    LOGICAL :: lmfscv, lmfmid, lmfpen
    ! switch for detrainment of rain and snow to gridscale scheme
    LOGICAL :: lmfdsnow
    !
    ! Parameters which are only computed if Gravity wave drag scheme is switched on
    !
    ! launch level for GWD scheme
    INTEGER  :: klaunch
    !
    ! Parameters which are only computed if Sub-grid Scale Orography (SSO) scheme is switched on
    !
    INTEGER  :: ngwdlim, ngwdtop, nktopg
    REAL(wp) :: gkwake, gkdrag, gfrcrit, grcrit
    !
    ! Parameters which are always computed
    !
    ! characteristic horizontal length scale (grid-scale) for 
    ! turbulence scheme and convection scheme
    REAL(wp) :: mean_charlen
    ! level index corresponding to the HAG of the 60hPa level (identical to kcon2, if computed)
    INTEGER :: k060
  END TYPE t_phy_params


  PUBLIC :: t_phy_params


END MODULE mo_nwp_parameters
