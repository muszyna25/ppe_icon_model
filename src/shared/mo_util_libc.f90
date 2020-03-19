!> This module contains bindings to the functions of C standard library.
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!! -----------------------------------------------------------------------------------
MODULE mo_util_libc

  USE ISO_C_BINDING, ONLY: C_PTR, C_INT, C_SIZE_T

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: memset
  PUBLIC :: memcmp

CONTAINS

  SUBROUTINE memset(str, c, n)
    TYPE(C_PTR), VALUE :: str
    INTEGER(C_INT), VALUE :: c
    INTEGER(C_SIZE_T), VALUE :: n

    INTERFACE
      FUNCTION c_memset(str, c, n) BIND(C, NAME='memset') RESULT(c_result)
        IMPORT C_PTR, C_INT, C_SIZE_T
        TYPE(C_PTR), VALUE :: str
        INTEGER(C_INT), VALUE :: c
        INTEGER(C_SIZE_T), VALUE :: n
        TYPE(C_PTR) :: c_result
      END FUNCTION c_memset
    END INTERFACE

    str = c_memset(str, c, n)
  END SUBROUTINE memset

  FUNCTION memcmp(str1, str2, n) RESULT(f_result)
    TYPE(C_PTR), VALUE, INTENT(IN) :: str1, str2
    INTEGER(C_SIZE_T), VALUE :: n
    LOGICAL :: f_result

    INTERFACE
      FUNCTION c_memcmp(str1, str2, n) BIND(C, NAME='memcmp') RESULT(c_result)
        IMPORT C_PTR, C_INT, C_SIZE_T
        TYPE(C_PTR), VALUE, INTENT(IN) :: str1, str2
        INTEGER(C_SIZE_T), VALUE :: n
        INTEGER(C_INT) :: c_result
      END FUNCTION c_memcmp
    END INTERFACE

    f_result = (c_memcmp(str1, str2, n) /= 0)
  END FUNCTION memcmp

END MODULE mo_util_libc

