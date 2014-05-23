!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
MODULE mo_canopy

  USE mo_kind, ONLY: dp => wp

  IMPLICIT NONE

  ! Parameters for computation of unstressed canopy resistance using Eq. 3.3.2.12 in ECHAM3 manual
  REAL(dp), PARAMETER :: &
       k = 0.9_dp,      &
       a = 5000._dp,    &
       b = 10._dp,      &
       c = 100._dp

  PUBLIC :: unstressed_canopy_cond_par

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

CONTAINS
  
  SUBROUTINE unstressed_canopy_cond_par(lai, par, conductance)


    ! Computes the canopy conductance without water stress limitation according to ECHAM formalism, Eq. 3.3.2.12 in 
    ! ECHAM3 Manual

    REAL(dp), INTENT(in), DIMENSION(:) :: &
         lai,      &
         par
    REAL(dp), INTENT(inout) :: conductance(SIZE(lai,DIM=1)) ! out
  

    ! Local variables
    REAL(dp) :: d(SIZE(par))    !! Variable d in Eq. 3.3.2.12 in ECHAM3 manual
    REAL(dp) :: par_echam(SIZE(par,DIM=1))
    par_echam = max(1.e-10_dp, par)

    WHERE(lai > EPSILON(1._dp))
       d = (a + b*c) / (c * par_echam)
       conductance = (LOG((d * EXP(k*lai) + 1._dp) / &
                     (d+1._dp)) * b / (d*par_echam) - LOG((d + EXP(-k*lai)) / &
                     (d+1._dp)))/(k*c)
    ELSEWHERE
       conductance = 1.e-20_dp
    END WHERE
    
  END SUBROUTINE unstressed_canopy_cond_par

END MODULE mo_canopy
