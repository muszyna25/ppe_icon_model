!>
!! Idealized q-flux corrections for the slab ocean model 
!! ...to run more complex cases, this will need further development
!!
!! @author Andrew IL Williams, Oxford
!!
!! @par Revision History
!! First version by Andrew IL Williams (2022-08-28)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_ml_qflux

  USE mo_kind,               ONLY: wp

  USE mo_math_constants,     ONLY: pi
#ifndef __NO_ICON_ATMO__
  USE mo_nh_testcases_nml,   ONLY: qflux_case, qflux_const
#endif
  USE mo_impl_constants,     ONLY: max_char_length
  USE mo_exception,          ONLY: finish

  IMPLICIT NONE

  PUBLIC
  PRIVATE :: simple_qflux_zero, simple_qflux_const, simple_qflux_mockWalker


CONTAINS
  
  !! ZERO QFLUX
  ELEMENTAL FUNCTION simple_qflux_zero( ) RESULT(qflux)

    REAL(wp)             :: qflux

    qflux = 0

  END FUNCTION simple_qflux_zero


  !! CONSTANT VALUE
  FUNCTION simple_qflux_const( ) RESULT(qflux)

    REAL(wp) :: qflux

    qflux = qflux_const

  END FUNCTION simple_qflux_const


  !! Linear tapering following Bretherton 2007
  ELEMENTAL FUNCTION simple_qflux_mockwalker( lon ) RESULT(qflux)

    REAL(wp),INTENT(in)  :: lon
    REAL(wp)             :: qflux

    ! n.b. ICON uses longitude in [-pi, pi] range!
    qflux = qflux_const + 50._wp * ( ABS(lon/pi) - 0.5_wp )

  END FUNCTION simple_qflux_mockwalker

  !-----------

  FUNCTION ml_qflux(qflux_case,lat,lon) RESULT(qflux)
    REAL(wp),INTENT(IN)                       :: lat
    REAL(wp),INTENT(IN)                       :: lon
    CHARACTER(len=max_char_length),INTENT(IN) :: qflux_case
    REAL(wp)                                  :: qflux

    CHARACTER(len=max_char_length), PARAMETER :: FUNCTION = '(mo_ml_qflux) ml_qflux:'

    SELECT CASE (qflux_case)
    CASE ('qflux_zero')
      qflux=simple_qflux_zero( )
    CASE ('qflux_const')
      qflux=simple_qflux_const( )
    CASE ('qflux_mockwalker')
      qflux=simple_qflux_mockwalker( lon )
    CASE DEFAULT
      CALL finish( TRIM(FUNCTION),'wrong qflux name, must be qflux_zero, qflux_const or qflux_mockwalker')
    END SELECT
  END FUNCTION ml_qflux

END MODULE mo_ml_qflux

