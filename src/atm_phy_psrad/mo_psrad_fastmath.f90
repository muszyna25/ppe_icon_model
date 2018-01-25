!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_psrad_fastmath
  USE mo_psrad_general, ONLY : wp, dp
#include "psrad_fastmath.inc"
  IMPLICIT NONE
  PRIVATE 
  PUBLIC :: transmit, tautrans
  
  REAL(wp), PARAMETER :: rec_6  = 1._wp/6._wp

CONTAINS

  FUNCTION transmit(x, n) 
    ! Compute transmittance 1 - EXP(-x) 
    INTEGER,  INTENT(in) :: n
    REAL(DP), INTENT(IN) :: x(n) 
    REAL(DP)             :: transmit(n) 

    ! MASS and MKL libraries have exp(x) - 1 functions; we could 
    !   use those here
    INV_EXPON(x(1:n),transmit(1:n))
    transmit(1:n) = 1._wp - transmit(1:n)
  END FUNCTION transmit

  FUNCTION tautrans(x, n) 
    ! Compute "tau transition" using linear-in-tau approximation
    INTEGER,  INTENT(in) :: n
    REAL(DP), INTENT(IN) :: x(n) 
    REAL(DP)             :: tautrans(n) 
    REAL(DP) :: y(n) 
    ! Default calculation is unstable (NaN) for the very lowest value 
    ! of tau (3.6e-4) 
    INV_EXPON(x(1:n),y(1:n))
    WHERE (x > 1.E-3_wp)
      tautrans = 1._wp - 2._wp*(1._wp/x - y/(1._wp-y))
    ELSEWHERE
      tautrans = x * rec_6
    END WHERE
  END FUNCTION tautrans

END MODULE mo_psrad_fastmath
