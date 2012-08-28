MODULE mo_util_signal

  USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT

  IMPLICIT NONE

  INTERFACE
    FUNCTION signal_trap(create_dump, signals) RESULT(iret) BIND(C,NAME='signal_trap')
#if defined(__SX__) || defined (__SUNPRO_F95)
      USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT
#else
      IMPORT :: C_INT
#endif
      INTEGER(C_INT) :: iret
      INTEGER(C_INT), INTENT(in) :: create_dump
      INTEGER(C_INT), INTENT(in) :: signals(*) 
    END FUNCTION signal_trap
  END INTERFACE

END MODULE mo_util_signal
