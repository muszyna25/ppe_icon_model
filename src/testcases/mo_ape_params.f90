!>
!! Physical constants and SST specification for the Aqua-Planet Experiments.
!! See <http://www.met.reading.ac.uk/~mike/APE/ape_spec.html>
!!
!! @author Hui Wan, MPI-M
!!
!! @par Revision History
!! First version by Hui Wan (2010-07-08)
!! added experimental setups 2-4 by Kristina Froehlich (2010-10-05s)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_ape_params

  USE mo_kind,               ONLY: wp

  USE mo_math_constants,     ONLY: pi
  USE mo_physical_constants, ONLY: tmelt
#ifndef __NO_ICON_ATMO__
  USE mo_nh_testcases_nml,   ONLY: ape_sst_val
#endif
  USE mo_impl_constants,     ONLY: max_char_length
  USE mo_exception,          ONLY: finish

  IMPLICIT NONE

  PUBLIC
  PRIVATE :: ape_sst1, ape_sst2, ape_sst3, ape_sst4, ape_sst_qobs, ape_sst_ice, ape_sst_const

  CHARACTER(len=*),PARAMETER,PRIVATE :: version = '$Id$'

  ! Requirements for APE

  REAL(wp),PARAMETER :: sst_min = tmelt        !< no sea ice in APE
  REAL(wp),PARAMETER :: sirrad  = 1365._wp     !< constant solar irradiance [W/m^2]
  REAL(wp),PARAMETER :: eccent  = 0._wp        !< Earth orbit param.
  REAL(wp),PARAMETER :: obliq   = 0._wp        !< Earth orbit param.
  REAL(wp),PARAMETER :: co2     = 348._wp      !< ppmv (as in AMIP II)

  ! Recommendations for APE

!   REAL(wp),PARAMETER :: omega = 7.292115e-5_wp !< Earth rotation rate [1/m]
!   REAL(wp),PARAMETER :: re    = 6.371e6_wp     !< mean Earth radius [m]
  REAL(wp),PARAMETER :: grav  = 9.79764_wp     !< mean surface gravity [m/s^2]
  REAL(wp),PARAMETER :: rd    = 287.04_wp      !< gas const. for dry air [J/kg/K]
  REAL(wp),PARAMETER :: cpd   = 1004.64_wp     !< specific heat capacity for dry air [J/kg/K]
  REAL(wp),PARAMETER :: rv    = 461.50_wp      !< water vapour gas constant [J/kg/K]
  REAL(wp),PARAMETER :: cpv   = 1846.0_wp      !< water vapour specific heat capacity [J/kg/K]
  REAL(wp),PARAMETER :: alv   = 2.501e6_wp     !< latent heat of vaporization at 0 degC [J/kg]
  REAL(wp),PARAMETER :: alf   = 3.337e5_wp     !< latent heat of fusion       at 0 degC [J/kg]
  REAL(wp),PARAMETER :: als   = 2.834e6_wp     !< latent heat of sublimation  at 0 degC [J/kg]

  REAL(wp),PARAMETER :: ch4   = 1650._wp       !< ppbv
  REAL(wp),PARAMETER :: n2o   = 306._wp        !< ppbv
  REAL(wp),PARAMETER :: haloc = 0.2_wp         !< radiative forcing from Halocarbon [W/m^2]

  REAL(wp),PARAMETER :: ps_dry_init = 101325._wp - 245._wp !< initial mean surface pressure [Pa]
  REAL(wp),PARAMETER :: q_tot_dry   = 25.006_wp            !< initial value of
  !< global moisture content [kg/m^2]
CONTAINS
  !>
  !! Zonally symmetric analytic distributions of sea surface temperature (SST)
  !! for the "control" case of the APE.
  !!
  ELEMENTAL FUNCTION ape_sst1( lat ) RESULT(sst)

    REAL(wp),INTENT(in)  :: lat
    REAL(wp)             :: sst

    REAL(wp)             :: lat0

    lat0     = pi/3._wp

    IF (ABS(lat)<lat0 ) THEN
      sst = tmelt + 27._wp*( 1._wp - SIN(1.5_wp*lat)**2 )
    ELSE
      sst = tmelt
    END IF

  END FUNCTION ape_sst1


  !> Zonally symmetric analytic distributions of sea surface temperature (SST)
  !! for the "peaked" case of the APE.
  !!
  ELEMENTAL FUNCTION ape_sst2( lat ) RESULT(sst)

    REAL(wp),INTENT(in)  :: lat
    REAL(wp)             :: sst

    REAL(wp)             :: lat0

    lat0     = pi/3._wp

    IF (ABS(lat) < lat0 ) THEN
      sst = tmelt + 27._wp*( 1._wp - (3._wp/pi*ABS(lat)) )
    ELSE
      sst = tmelt
    END IF

  END FUNCTION ape_sst2


  !> Zonally symmetric analytic distributions of sea surface temperature (SST)
  !! for the "flat" case of the APE.
  !!
  ELEMENTAL FUNCTION ape_sst3( lat ) RESULT(sst)

    REAL(wp),INTENT(in)  :: lat
    REAL(wp)             :: sst

    REAL(wp)             :: lat0

    lat0     = pi/3._wp

    IF (ABS(lat) < lat0 ) THEN
      sst = tmelt + 27._wp*( 1._wp - SIN(1.5_wp*lat)**4 )
    ELSE
      sst = tmelt
    END IF

  END FUNCTION ape_sst3


  !>
  !! Zonally symmetric analytic distributions of sea surface temperature (SST)
  !! with temperature peak at 5 N  of the APE.
  !!
  ELEMENTAL FUNCTION ape_sst4( lat ) RESULT(sst)

    REAL(wp),INTENT(in)  :: lat
    REAL(wp)             :: sst
    REAL(wp)             :: lat0, lat5

    lat0     = pi/3._wp
    lat5     = pi/36._wp

    IF (lat > lat5 .AND. lat < lat0 ) THEN
      sst = 27._wp*( 1._wp - SIN(1.6363_wp*(lat-lat5))**2 ) + tmelt
    ELSE IF (lat > -1._wp*lat0 .AND. lat < lat5 ) THEN
      sst = 27._wp*( 1._wp - SIN(1.3846_wp*(lat-lat5))**2 ) + tmelt
    ELSE
      sst = tmelt
    END IF

  END FUNCTION ape_sst4


  !>
  !! @brief APE SST, case "Qobs".
  !!
  ELEMENTAL FUNCTION ape_sst_qobs(lat) RESULT(sst)

    REAL(wp),INTENT(in)  :: lat
    REAL(wp)             :: sst

    REAL(wp)             :: zsst1, zsst3
    REAL(wp)             :: zlat0, zsinlat

    zlat0 = pi/3._wp

    IF (ABS(lat) < zlat0 ) THEN
      zsinlat = SIN(1.5_wp*lat)
      zsst1   = 27._wp*( 1._wp - zsinlat**2 )
      zsst3   = 27._wp*( 1._wp - zsinlat**4 )
      sst     = tmelt + 0.5_wp*( zsst1 + zsst3 )
    ELSE
      sst     = tmelt
    END IF

  END FUNCTION ape_sst_qobs

  !!KF attempt to generate freezing temperatures for the ocean
 !! Zonally symmetric analytic distributions of sea surface temperature (SST)
  !! for the "control" case of the APE.
  !!
  ELEMENTAL FUNCTION ape_sst_ice( lat ) RESULT(sst)

    REAL(wp),INTENT(in)  :: lat
    REAL(wp)             :: sst

    REAL(wp)             :: lat0

    lat0     = pi/2.8_wp

    IF (ABS(lat)<lat0 ) THEN
      sst = tmelt -1.9_wp + 27._wp*( 1._wp - SIN(1.5_wp*lat)**2 )
    ELSE
      sst = tmelt -1.9_wp
    END IF

  END FUNCTION ape_sst_ice


  !! Globally constant sea surface temperature (SST)
  !! for the DCMIP tropical cyclone case of the APE.
  !!
  FUNCTION ape_sst_const( ) RESULT(sst)

    REAL(wp) :: sst

#ifndef __NO_ICON_ATMO__
    sst = tmelt + ape_sst_val
#else
    sst = tmelt + 29.0_wp
#endif

  END FUNCTION ape_sst_const


  !-----------

  FUNCTION ape_sst(sst_case,lat) RESULT(sst)
    REAL(wp),INTENT(IN)                       :: lat
    CHARACTER(len=max_char_length),INTENT(IN) :: sst_case
    REAL(wp)                                  :: sst

    CHARACTER(len=max_char_length), PARAMETER :: FUNCTION = '(mo_ape_params) ape_sst:'

    SELECT CASE (sst_case)
    CASE ('sst1')
      sst=ape_sst1(lat)
    CASE ('sst2')
      sst=ape_sst2(lat)
    CASE ('sst3')
      sst=ape_sst3(lat)
    CASE ('sst4')
      sst=ape_sst4(lat)
    CASE ('sst_qobs')
      sst=ape_sst_qobs(lat)
    CASE ('sst_ice')
      sst=ape_sst_ice(lat)
    CASE ('sst_const')
      sst=ape_sst_const()
    CASE DEFAULT
      CALL finish( TRIM(FUNCTION),'wrong sst name, must be sst1, sst2, sst3, sst4 or sst_qobs')
    END SELECT
  END FUNCTION ape_sst

END MODULE mo_ape_params

