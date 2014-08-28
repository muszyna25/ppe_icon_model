!>
!! Utility module:
!!
!! Fortran-C-Interface for string parsing routine.
!!
!! Initial revision: 08/2014 : F. Prill, DWD
!!
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!

MODULE mo_util_string_parse

  USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT, C_CHAR, C_NULL_CHAR

  IMPLICIT NONE

  PRIVATE

  INTERFACE
    SUBROUTINE private_do_parse_intlist(parse_line, nvalues, out_values, ierr) BIND(C,NAME='do_parse_intlist') 
#if defined(__SX__) || defined (__SUNPRO_F95) 
      USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT, C_CHAR
      CHARACTER(kind=C_CHAR,len=*), INTENT(IN) :: parse_line
#else
      IMPORT :: C_INT, C_CHAR
      CHARACTER(C_CHAR), DIMENSION(*), INTENT(IN) :: parse_line
#endif
      INTEGER(C_INT), VALUE, INTENT(IN)    :: nvalues
      INTEGER(C_INT), INTENT(INOUT) :: out_values(nvalues)
      INTEGER(C_INT), INTENT(INOUT) :: ierr
    END SUBROUTINE private_do_parse_intlist
  END INTERFACE

  PUBLIC :: util_do_parse_intlist

CONTAINS

  ! ---------------------------------------------------------------------
  ! Subroutine parsing the string parse_line containing integer numbers.  
  !
  !    Allowed is a comma- (or semicolon-) separated list of integers,
  !    and of integer ranges like "10...20".  One may also use the
  !    keyword "nlev" to denote the maximum integer (or, equivalently,
  !    "n" or "N").
  !
  !    Furthermore, arithmetic expressions like "(nlev - 2)" are
  !    possible.
  !
  !    Basic example:
  !       parse_line = "1,3,5...10,20...nlev"
  !    More complex example:
  !       parse_line = "1,2, 10 ...22;2;16-(3+11), N-2,16-(2+10);5"
  !
  ! 08/2014 : F. Prill, DWD
  !
  ! @param[in]  parse_line     string containing integer numbers
  ! @param[out] out_values     out_values[i] = 1 if "i" was in parse_line
  ! @param[out] ierr           error code != 0 if parser failed
  ! ---------------------------------------------------------------------
  SUBROUTINE util_do_parse_intlist(parse_line, nlev_value, out_values, ierr)
    CHARACTER(len=*), INTENT(IN)    :: parse_line
    INTEGER,          INTENT(IN)    :: nlev_value      !< number to substitute for "N"/"nlev"
    INTEGER,          INTENT(INOUT) :: out_values(0:)
    INTEGER,          INTENT(INOUT) :: ierr

    CALL private_do_parse_intlist(TRIM(parse_line)//C_NULL_CHAR, nlev_value, out_values, ierr)
  END SUBROUTINE util_do_parse_intlist

END MODULE mo_util_string_parse
