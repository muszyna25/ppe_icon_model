!>
!! Constants used in the cumulus convection parameterization scheme 
!! of the ECHAM physics package
!!
!! @author Hui Wan, MPI-M
!!
!! @par Revision History
!! First version by Hui Wan, MPI (2010-07)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_echam_conv_constants

  USE mo_kind, ONLY: wp

  IMPLICIT NONE
  PUBLIC
  CHARACTER(len=*), PARAMETER, PRIVATE :: version = '$Id$'

  LOGICAL, PARAMETER :: lmfmid   = .true.    !< true when midlevel    convection is switched on
  LOGICAL, PARAMETER :: lmfdd    = .true.    !< true when cumulus downdraft      is switched on
  LOGICAL, PARAMETER :: lmfdudv  = .true.    !< true when cumulus friction       is switched on

  REAL(wp),PARAMETER :: entrmid  = 1.0E-4_wp !< average entrainment rate for midlevel convection
  REAL(wp),PARAMETER :: entrscv  = 3.0E-4_wp !< average entrainment rate for shallow convection
  REAL(wp),PARAMETER :: entrpen  = 1.0E-4_wp !< entrainment rate for penetrative convection
  REAL(wp),PARAMETER :: entrdd   = 2.0E-4_wp !< average entrainment rate for cumulus downdrafts

  REAL(wp),PARAMETER :: cprcon   = 2.E-4_wp  !< coefficient for determining conversion
                                             !< from cloud water to rain
  REAL(wp),PARAMETER :: cmfctop  = 0.21_wp   !< fractional convective mass flux across the top of cloud
  REAL(wp),PARAMETER :: cmfdeps  = 0.3_wp    !< fractional convective mass flux for downdrafts at lfs
  
  REAL(wp),PARAMETER :: cmfcmin  = 1.E-10_wp !< minimum massflux value (for safety)
  REAL(wp),PARAMETER :: cmfcmax  = 1.0_wp    !< maximum massflux value allowed for updrafts etc

  REAL(wp),PARAMETER :: cminbuoy = 0.1_wp    !< minimum excess buoyancy
  REAL(wp),PARAMETER :: cmaxbuoy = 1.0_wp    !< maximum excess buoyancy
  REAL(wp),PARAMETER :: cbfac    = 1.0_wp    !< factor for std dev of virtual pot temp
  REAL(wp),PARAMETER :: centrmax = 3.E-4_wp  !<

  REAL(wp),PARAMETER :: dlev     = 3.0E4_wp  !< "zdlev" in subroutine "cuasc".
                                             !< Critical thickness (unit: Pa) necessary for the
                                             !< onset of convective precipitation
  REAL(wp),PARAMETER :: cmftau   = 10800._wp !< characteristic adjustment time scale (s)
                                             !< (replaces "ztau" in "cumastr")

END MODULE mo_echam_conv_constants
