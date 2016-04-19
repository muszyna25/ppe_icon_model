!> Wrapper module containing the Fortran-C-Interface for the
!  Fortran namelist scanner.
!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_util_nml

  USE, INTRINSIC ::  ISO_C_BINDING, ONLY: C_INT, C_CHAR, C_NULL_CHAR

  IMPLICIT NONE

  PRIVATE

  INTERFACE
    FUNCTION private_annotate_nml(in_filename, out_filename) RESULT(iret) &
      &      BIND(C,NAME='util_annotate_nml')
#if defined (__SUNPRO_F95)
      USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT, C_CHAR
#else
      IMPORT :: C_INT, C_CHAR
#endif
      INTEGER(C_INT) :: iret
#if defined (__SUNPRO_F95)
      CHARACTER(kind=C_CHAR,len=*),    INTENT(in) :: in_filename, out_filename
#else
      CHARACTER(C_CHAR), DIMENSION(*), INTENT(in) :: in_filename, out_filename
#endif
    END FUNCTION private_annotate_nml
  END INTERFACE


  PUBLIC :: util_annotate_nml

CONTAINS

  FUNCTION util_annotate_nml(in_filename, out_filename) RESULT(iret)
    INTEGER :: iret
    CHARACTER(len=*), INTENT(in) :: in_filename, out_filename
    iret = private_annotate_nml(TRIM(in_filename)//C_NULL_CHAR, &
      &                         TRIM(out_filename)//C_NULL_CHAR)
  END FUNCTION util_annotate_nml

END MODULE mo_util_nml
