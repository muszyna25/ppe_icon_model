MODULE mo_util_system

  USE ISO_C_BINDING, ONLY: C_INT, C_CHAR, C_NULL_CHAR

  IMPLICIT NONE

  PRIVATE

#ifdef __XT3__
  PUBLIC :: util_base_iobuf
#endif
  PUBLIC :: util_exit
  PUBLIC :: util_abort
  PUBLIC :: util_system

  INTERFACE
#ifdef __XT3__
    SUBROUTINE util_base_iobuf() BIND(C)
    END SUBROUTINE util_base_iobuf
#endif
    SUBROUTINE util_exit(exit_no) BIND(C)
      IMPORT C_INT
      INTEGER(C_INT), VALUE :: exit_no
    END SUBROUTINE util_exit
    SUBROUTINE util_abort() BIND(C)
    END SUBROUTINE util_abort
  END INTERFACE

CONTAINS

  FUNCTION util_system(f_s) RESULT(f_result)
    INTEGER(C_INT) :: f_result
    CHARACTER(KIND = C_CHAR, LEN = *), INTENT(IN) :: f_s

    CHARACTER(KIND = C_CHAR) :: c_s(LEN(f_s) + 1)
    INTEGER :: i

    INTERFACE
      FUNCTION c_util_system(c_s) BIND(C, NAME = 'util_system') RESULT(c_result)
        IMPORT C_INT, C_CHAR
        INTEGER(C_INT) :: c_result
        CHARACTER(KIND = C_CHAR), INTENT(IN) :: c_s(*)
      END FUNCTION c_util_system
    END INTERFACE
    DO i = 1, LEN(f_s)
      c_s(i) = f_s(i:i)
    END DO
    c_s(LEN(f_s) + 1) = C_NULL_CHAR
    f_result = c_util_system(c_s)
  END FUNCTION util_system

END MODULE mo_util_system

