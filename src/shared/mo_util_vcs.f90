!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
MODULE mo_util_vcs

  USE, INTRINSIC :: iso_c_binding, ONLY: c_int, c_char, c_null_char

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: util_repository_url
  PUBLIC :: util_branch_name
  PUBLIC :: util_revision_key

  INTERFACE

    SUBROUTINE private_util_repository_url(name, actual_len) BIND(c,name='repository_url')
      IMPORT :: c_int, c_char
      CHARACTER(c_char), DIMENSION(*), INTENT(inout) :: name
      INTEGER(c_int), INTENT(inout) :: actual_len
    END SUBROUTINE private_util_repository_url

    SUBROUTINE private_util_branch_name(name, actual_len) BIND(c,name='branch_name')
      IMPORT :: c_int, c_char
      CHARACTER(c_char), DIMENSION(*), INTENT(inout) :: name
      INTEGER(c_int), INTENT(inout) :: actual_len
    END SUBROUTINE private_util_branch_name

    SUBROUTINE private_util_revision_key(name, actual_len) BIND(c,name='revision_key')
      IMPORT :: c_int, c_char
      CHARACTER(c_char), DIMENSION(*), INTENT(inout) :: name
      INTEGER(c_int), INTENT(inout) :: actual_len
    END SUBROUTINE private_util_revision_key

  END INTERFACE

CONTAINS

  SUBROUTINE util_repository_url(name, actual_len)
    CHARACTER(len=*), INTENT(out) :: name
    INTEGER, INTENT(inout) :: actual_len
    INTEGER :: i
    CALL private_util_repository_url(name, actual_len)
    char_loop: DO i = 1 , LEN(name)
      IF (name(i:i) == c_null_char) EXIT char_loop
    END DO char_loop
    name(i:LEN(name)) = ' '
  END SUBROUTINE util_repository_url

  SUBROUTINE util_branch_name(name, actual_len)
    CHARACTER(len=*), INTENT(out) :: name
    INTEGER, INTENT(inout) :: actual_len
    INTEGER :: i
    CALL private_util_branch_name(name, actual_len)
    char_loop: DO i = 1 , LEN(name)
      IF (name(i:i) == c_null_char) EXIT char_loop
    END DO char_loop
    name(i:LEN(name)) = ' '
  END SUBROUTINE util_branch_name

  SUBROUTINE util_revision_key(name, actual_len)
    CHARACTER(len=*), INTENT(out) :: name
    INTEGER, INTENT(inout) :: actual_len
    INTEGER :: i
    CALL private_util_revision_key(name, actual_len)
    char_loop: DO i = 1 , LEN(name)
      IF (name(i:i) == c_null_char) EXIT char_loop
    END DO char_loop
    name(i:LEN(name)) = ' '
  END SUBROUTINE util_revision_key

END MODULE mo_util_vcs
