MODULE mo_util_hash

  USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT, C_CHAR

  IMPLICIT NONE

  PRIVATE

  INTERFACE
    FUNCTION util_hashword(text, text_len, inithash) RESULT(hash) BIND(C,NAME='util_hashword')
#ifdef __SX__
      USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT, C_CHAR
#else
      IMPORT :: C_INT, C_CHAR
#endif
      INTEGER(C_INT) :: hash
      CHARACTER(C_CHAR), DIMENSION(*), INTENT(in) :: text
      INTEGER(C_INT), VALUE,           INTENT(in) :: text_len
      INTEGER(C_INT), VALUE,           INTENT(in) :: inithash
    END FUNCTION util_hashword
  END INTERFACE

  PUBLIC :: utiL_hashword

END MODULE mo_util_hash

