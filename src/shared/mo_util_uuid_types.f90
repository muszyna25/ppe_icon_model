!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.

MODULE mo_util_uuid_types
  
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_SIGNED_CHAR, C_INT
  
  IMPLICIT NONE 
  
  PRIVATE
  
  PUBLIC :: t_uuid
  PUBLIC :: UUID_STRING_LENGTH
  PUBLIC :: UUID_DATA_LENGTH
  
  INTEGER, PARAMETER :: UUID_STRING_LENGTH = 36
  INTEGER, PARAMETER :: UUID_DATA_LENGTH   = 16
  
  !> module name string
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_util_uuid_types'
  
  TYPE, BIND(C) :: t_uuid
    INTEGER(C_SIGNED_CHAR) :: DATA(16)
  END TYPE t_uuid

END MODULE mo_util_uuid_types
