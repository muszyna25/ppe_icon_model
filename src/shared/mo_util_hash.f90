MODULE mo_util_hash

  USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT, C_CHAR

  IMPLICIT NONE

  PRIVATE

  INTERFACE
    FUNCTION hashword(text, text_len, inithash) RESULT(hash) &
         BIND(C,NAME='hashword')
      USE ISO_C_BINDING, ONLY: C_INT, C_CHAR
!DR      IMPORT :: C_INT, C_CHAR
      INTEGER(C_INT) :: hash
      CHARACTER(C_CHAR), DIMENSION(*), INTENT(in) :: text
      INTEGER(C_INT) :: text_len
      INTEGER(C_INT) :: inithash
    END FUNCTION hashword
  END INTERFACE

END MODULE mo_util_hash

