!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
MODULE mo_util_hash

  USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT, C_CHAR, C_SIZE_T

  IMPLICIT NONE

  PRIVATE

  INTERFACE
    FUNCTION util_hashword(text, text_len, inithash) RESULT(hash) BIND(C,NAME='util_hashword')
#if defined(__SUNPRO_F95)
      USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT, C_CHAR, C_SIZE_T
#else
      IMPORT :: C_INT, C_CHAR, C_SIZE_T
#endif
      INTEGER(C_INT) :: hash
      CHARACTER(C_CHAR), DIMENSION(*), INTENT(in) :: text
      INTEGER(kind=c_size_t), VALUE,   INTENT(in) :: text_len
      INTEGER(C_INT), VALUE,           INTENT(in) :: inithash
    END FUNCTION util_hashword
  END INTERFACE

  PUBLIC :: utiL_hashword

END MODULE mo_util_hash

