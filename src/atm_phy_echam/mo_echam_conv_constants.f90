!>
!! Constants used in the cumulus convection parameterization scheme 
!! of the ECHAM physics package
!!
!! @author Hui Wan, MPI-M
!!
!! @par Revision History
!! First version by Hui Wan, MPI (2010-07)
!!
!! @par Copyright
!! 2002-2011 by DWD and MPI-M
!! This software is provided for non-commercial use only.
!! See the LICENSE and the WARRANTY conditions.
!!
!! @par License
!! The use of ICON is hereby granted free of charge for an unlimited time,
!! provided the following rules are accepted and applied:
!! <ol>
!! <li> You may use or modify this code for your own non commercial and non
!!    violent purposes.
!! <li> The code may not be re-distributed without the consent of the authors.
!! <li> The copyright notice and statement of authorship must appear in all
!!    copies.
!! <li> You accept the warranty conditions (see WARRANTY).
!! <li> In case you intend to use the code commercially, we oblige you to sign
!!    an according license agreement with DWD and MPI-M.
!! </ol>
!!
!! @par Warranty
!! This code has been tested up to a certain level. Defects and weaknesses,
!! which may be included in the code, do not establish any warranties by the
!! authors.
!! The authors do not make any warranty, express or implied, or assume any
!! liability or responsibility for the use, acquisition or application of this
!! software.
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
