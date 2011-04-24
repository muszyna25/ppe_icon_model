MODULE mo_util_system

  USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT, C_CHAR, C_NULL_CHAR

  IMPLICIT NONE

  PRIVATE

  INTERFACE
    SUBROUTINE util_exit(exit_no)  BIND(C,NAME='util_exit')
      IMPORT :: C_INT
      INTEGER(C_INT), INTENT(in) :: exit_no
    END SUBROUTINE util_exit
  END INTERFACE

  INTERFACE
    SUBROUTINE util_abort() BIND(C,NAME='util_abort')
    END SUBROUTINE util_abort
  END INTERFACE

  INTERFACE
    FUNCTION private_util_system(cmd) RESULT(iret) BIND(C,NAME='util_system') 
      IMPORT :: C_INT, C_CHAR
      INTEGER(C_INT) :: iret
      CHARACTER(C_CHAR), DIMENSION(*), INTENT(in) :: cmd
    END FUNCTION private_util_system
  END INTERFACE

  PUBLIC :: util_exit
  PUBLIC :: util_abort
  PUBLIC :: util_system

CONTAINS

  FUNCTION util_system(cmd) RESULT(iret)
    CHARACTER(len=*), INTENT(in) :: cmd
    INTEGER :: iret
    iret = private_util_system(TRIM(cmd)//C_NULL_CHAR)
  END FUNCTION util_system

END MODULE mo_util_system
