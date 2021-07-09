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

  PUBLIC :: memset, memcmp, memcpy

CONTAINS

  SUBROUTINE memset(str, c, n)
    TYPE(C_PTR), VALUE :: str
    INTEGER(C_INT), VALUE :: c
    INTEGER(C_SIZE_T), VALUE :: n

    INTERFACE
      TYPE(C_PTR) FUNCTION c_memset(str, c, n) BIND(C, NAME='memset')
        IMPORT C_PTR, C_INT, C_SIZE_T
        TYPE(C_PTR), VALUE :: str
        INTEGER(C_INT), VALUE :: c
        INTEGER(C_SIZE_T), VALUE :: n
      END FUNCTION c_memset
    END INTERFACE

    str = c_memset(str, c, n)
  END SUBROUTINE memset

  FUNCTION memcmp(str1, str2, n) RESULT(f_result)
    TYPE(C_PTR), VALUE, INTENT(IN) :: str1, str2
    INTEGER(C_SIZE_T), VALUE :: n
    LOGICAL :: f_result

    INTERFACE
      INTEGER(C_INT) FUNCTION c_memcmp(a, b, n) BIND(C, NAME='memcmp')
        IMPORT C_PTR, C_INT, C_SIZE_T
        TYPE(C_PTR), VALUE, INTENT(IN) :: a, b
        INTEGER(C_SIZE_T), VALUE :: n
      END FUNCTION c_memcmp
    END INTERFACE

    f_result = (c_memcmp(str1, str2, n) /= 0)
  END FUNCTION memcmp

  TYPE(C_PTR) FUNCTION memcpy(dest, src ,bsize) RESULT(ret)
    TYPE(C_PTR), VALUE :: dest
    TYPE(C_PTR), INTENT(IN), VALUE :: src
    INTEGER(C_SIZE_T), INTENT(IN) :: bsize

    INTERFACE
      TYPE(C_PTR) FUNCTION c_memcpy(a, b, s) BIND(C,NAME='memcpy')
        IMPORT C_SIZE_T, C_PTR
        TYPE(C_PTR), VALUE :: a
        TYPE(C_PTR), INTENT(IN), VALUE :: b
        INTEGER(C_SIZE_T), INTENT(IN), VALUE :: s
      END FUNCTION c_memcpy
    END INTERFACE

    ret = c_memcpy(dest, src, bsize)
  END FUNCTION memcpy

END MODULE mo_util_libc

