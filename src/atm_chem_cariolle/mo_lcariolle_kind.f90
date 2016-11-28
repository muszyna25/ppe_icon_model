!>
!! @brief This module defines the kinds of the numeric types of the Cariolle
!! scheme. They have to be exactly the same as in the host model
!! documentation: cr2016_10_22_rjs
!!
!! @author Sebastian Rast, MPI-M
!!
!! @par Revision History
!!  Original version Sebastian Rast (2016)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_lcariolle_kind
  USE mo_kind, ONLY: dp, i4
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: wp, wi
!!$  INTEGER, PARAMETER :: pd =  12
!!$  INTEGER, PARAMETER :: rd = 307
!!$  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(pd,rd)
!!$  INTEGER, PARAMETER :: pi8= 8
!!$  INTEGER, PARAMETER :: i8 = SELECTED_INT_KIND(pi8)  
!!$  !
!!$  INTEGER, PARAMETER :: wp = dp
!!$  INTEGER, PARAMETER :: wi = i8

INTEGER, PARAMETER :: wp=dp
INTEGER, PARAMETER :: wi=i4
END MODULE mo_lcariolle_kind
