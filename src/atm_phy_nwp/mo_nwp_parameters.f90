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
    ! Level parameters for convection scheme
    INTEGER  :: kcon1, kcon2, kcon3, kcon4, kcon5
    ! resolution-dependent parameters for convection scheme
    REAL(wp) :: tau, mfcfl, tau0
    ! characteristic horizontal length scale (grid-scale) for 
    ! turbulence scheme and convection scheme
    REAL(wp) :: mean_charlen
    ! relative humidity below which sub-cloud evaporation of rain starts over land and water, respectively
    REAL(wp) :: rhebc_land, rhebc_ocean
    ! launch level for GWD scheme
    INTEGER  :: klaunch
  END TYPE t_phy_params


  PUBLIC :: t_phy_params


END MODULE mo_nwp_parameters
