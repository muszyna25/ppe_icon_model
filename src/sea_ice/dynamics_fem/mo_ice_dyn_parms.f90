!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
MODULE mo_ice_dyn_parms
!introduced by Ralph Timmermann, 17.9.2004
USE mo_kind,    ONLY: wp
implicit none
PUBLIC
REAL(wp), parameter :: cdwin = 2.25e-3_wp   ! drag coefficient atmosphere - ice FR284
REAL(wp), parameter :: cdwat = 5.00e-3_wp   ! drag coefficient ocean - ice      FR284
REAL(wp), parameter :: cdao  = 1.20e-3_wp   ! drag coefficient atmosphere - ocean !FR284

END MODULE mo_ice_dyn_parms
!==============================================================================
