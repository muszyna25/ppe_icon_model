!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
MODULE mo_util_texthash

  USE ISO_C_BINDING, ONLY: C_INT, C_SIZE_T, C_CHAR, C_INT32_T
  USE mo_exception,  ONLY: finish

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: text_hash, text_hash_c, text_isEqual, sel_char

#ifdef __PGI
  TYPE, PUBLIC :: t_char_workaround
    CHARACTER(:), ALLOCATABLE :: c
  END TYPE t_char_workaround
#endif

  CHARACTER(*), PARAMETER :: modname = "mo_util_texthash"

CONTAINS

  FUNCTION sel_char(key, routine, err_msg) RESULT(ptr)
    CLASS(*), POINTER, INTENT(in) :: key
    CHARACTER(*), INTENT(IN) :: routine, err_msg
    CHARACTER(:), POINTER :: ptr

    SELECT TYPE(key)
#ifdef __PGI
    TYPE IS(t_char_workaround)
      ptr => key%c
#endif
    TYPE IS(CHARACTER(*))
      ptr => key
    CLASS DEFAULT
      CALL finish(routine, err_msg)
    END SELECT
  END FUNCTION sel_char

  INTEGER FUNCTION text_hash_c(key) RESULT(hash)
    CHARACTER(*), INTENT(IN), TARGET :: key
    CLASS(*), POINTER :: key_p

    key_p => key
    hash = INT(text_hash(key_p))
  END FUNCTION text_hash_c

  INTEGER(C_INT) FUNCTION text_hash(key) RESULT(hash)
    INTERFACE
      FUNCTION util_hashword(text, text_len, inithash) RESULT(hash) BIND(C,NAME='util_hashword')
#if defined(__SUNPRO_F95)
        USE ISO_C_BINDING, ONLY: C_INT32_T, C_CHAR, C_SIZE_T
#else
        IMPORT :: C_CHAR, C_SIZE_T, C_INT32_T
#endif
        INTEGER(C_INT32_T) :: hash
        CHARACTER(C_CHAR), DIMENSION(*), INTENT(IN) :: text
        INTEGER(kind=C_SIZE_T), VALUE,   INTENT(IN) :: text_len
        INTEGER(C_INT32_T), VALUE,       INTENT(IN) :: inithash
      END FUNCTION util_hashword
    END INTERFACE

    CLASS(*), POINTER, INTENT(IN) :: key
    CHARACTER(*), PARAMETER :: routine = modname//":text_hash_cs"
    CHARACTER(:), POINTER :: key_p

    key_p => sel_char(key, routine, "Unknown type for key.")
    hash = INT(util_hashword(key_p, INT(LEN(key_p), C_SIZE_T), 0_C_INT32_T), C_INT)
  END FUNCTION text_hash

  LOGICAL FUNCTION text_isEqual(keyA, keyB) RESULT(is_equal)
    CLASS(*), POINTER, INTENT(in) :: keyA, keyB
    CHARACTER(*), PARAMETER :: routine = modname//":text_isEqual_cs"
    CHARACTER(:), POINTER :: keyA_p, keyB_p

    keyA_p => sel_char(keyA, routine, "Unknown type for keyA.")
    keyB_p => sel_char(keyB, routine, "Unknown type for keyB.")
    is_equal = keyA_p == keyB_p
  END FUNCTION text_isEqual
END MODULE mo_util_texthash
