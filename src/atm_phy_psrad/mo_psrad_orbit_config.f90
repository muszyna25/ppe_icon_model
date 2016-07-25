!>
!! Configuration of the orbit used in psrad radiation.
!!
!! @author S. Rast, MPI-M
!!
!! @par Revision History
!! First version by S. Rast, MPI (2016-01-28)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_psrad_orbit_config

  USE mo_kind,     ONLY: wp

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_psrad_orbit_config, psrad_orbit_config

  !>
  !! Derived type
  !!
  TYPE t_psrad_orbit_config

  REAL(wp) :: cecc           !< Eccentricity of Earth's Orbit
  REAL(wp) :: cobld          !< Obliquity of Earth [Deg]
  LOGICAL  :: l_orbvsop87    !< .TRUE. for VSOP87 orbit, 
                             !< .FALSE. for Kepler orbit
  LOGICAL  :: l_sph_symm_irr !< .TRUE. for globally averaged irradiation (RCE)
                             !< .FALSE. for lat (lon) dependent irradiation

  END TYPE t_psrad_orbit_config

  !>
  !! The configuration state (variable).
  !!
  TYPE(t_psrad_orbit_config) :: psrad_orbit_config

END MODULE mo_psrad_orbit_config
