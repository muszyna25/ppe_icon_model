!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.

MODULE mo_util_uuid

  USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_CHAR, C_NULL_CHAR, C_DOUBLE, C_INT, C_SIGNED_CHAR

  IMPLICIT NONE 

  PRIVATE

  PUBLIC :: t_uuid

  PUBLIC :: uuid_generate
  PUBLIC :: compare_uuid

  PUBLIC :: uuid_parse
  PUBLIC :: uuid_unparse

  PUBLIC :: uuid2char
  PUBLIC :: char2uuid

  PUBLIC :: clear_uuid

  PUBLIC :: OPERATOR(==)

  PUBLIC :: UUID_STRING_LENGTH
  PUBLIC :: UUID_DATA_LENGTH
  PUBLIC :: UUID_EQUAL, UUID_EQUAL_LIMITED_ACCURACY, UUID_UNEQUAL

  INTEGER, PARAMETER :: UUID_STRING_LENGTH = 36
  INTEGER, PARAMETER :: UUID_DATA_LENGTH   = 16

  TYPE, BIND(C) :: t_uuid
    INTEGER(C_SIGNED_CHAR) :: data(16)
  END type t_uuid

  ENUM, BIND(c)
    ENUMERATOR :: UUID_EQUAL                  = 0
    ENUMERATOR :: UUID_EQUAL_LIMITED_ACCURACY = 1
    ENUMERATOR :: UUID_UNEQUAL                = 2
  END ENUM

  INTERFACE OPERATOR (==)
    MODULE PROCEDURE uuid_compare
  END INTERFACE OPERATOR (==)

  INTERFACE
    SUBROUTINE my_uuid_generate(val, nval, uuid) BIND(C,NAME='uuid_generate')
      IMPORT :: C_DOUBLE, C_INT, t_uuid
      REAL(c_double),   INTENT(IN)         :: val(*)
      INTEGER(c_int),   INTENT(IN), VALUE  :: nval
      TYPE(t_uuid),     INTENT(OUT)        :: uuid
    END SUBROUTINE my_uuid_generate
  END INTERFACE

  INTERFACE
    INTEGER(C_INT) FUNCTION my_compare_uuid(uuid_A, uuid_B, min_difference) BIND(C,NAME='compare_UUID')
      IMPORT :: t_uuid, C_INT, C_DOUBLE
      TYPE(t_uuid),     INTENT(IN), VALUE  :: uuid_A, uuid_B
      REAL(c_double),   INTENT(OUT)        :: min_difference
    END FUNCTION my_compare_uuid
  END INTERFACE

  INTERFACE
    SUBROUTINE my_uuid_unparse(uuid_string, uuid) BIND(C,NAME='uuid_unparse')
      IMPORT :: C_CHAR, t_uuid
      CHARACTER(C_CHAR), DIMENSION(*), INTENT(out) :: uuid_string
      TYPE(t_uuid),                    INTENT(in)  :: uuid
    END SUBROUTINE my_uuid_unparse
  END INTERFACE

  INTERFACE
    SUBROUTINE my_uuid_parse(uuid, uuid_string) BIND(C,NAME='uuid_parse')
      IMPORT :: C_CHAR, t_uuid
      TYPE(t_uuid),                    INTENT(out) :: uuid
      CHARACTER(C_CHAR), DIMENSION(*), INTENT(in)  :: uuid_string
    END SUBROUTINE my_uuid_parse
  END INTERFACE


  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_util_uuid'

CONTAINS
 
  SUBROUTINE uuid_parse(uuid_string, uuid)
    CHARACTER(len=*), INTENT(in)  :: uuid_string
    TYPE(t_uuid),     INTENT(out) :: uuid
    CALL my_uuid_parse(uuid, TRIM(uuid_string)//C_NULL_CHAR)
  END SUBROUTINE uuid_parse

  SUBROUTINE uuid_unparse(uuid, uuid_string)
    TYPE(t_uuid),     INTENT(in)  :: uuid
    CHARACTER(len=*), INTENT(out) :: uuid_string
    CALL my_uuid_unparse(uuid_string, uuid)
  END SUBROUTINE uuid_unparse

  FUNCTION uuid_compare(uuid1, uuid2)
    LOGICAL :: uuid_compare
    TYPE(t_uuid),     INTENT(in)  :: uuid1
    TYPE(t_uuid),     INTENT(in)  :: uuid2
    uuid_compare = .TRUE.
!CDIR NOVECTOR
    IF (ANY(uuid1%data /= uuid2%data)) THEN
      uuid_compare = .FALSE.
    ENDIF
  END FUNCTION uuid_compare
  
  SUBROUTINE uuid2char(uuid, string)
    TYPE(t_uuid), INTENT(in) :: uuid
    CHARACTER(len=1), INTENT(out) :: string(16)
    string = TRANSFER(uuid%data, string)
  END SUBROUTINE uuid2char

  SUBROUTINE char2uuid(string, uuid)
    CHARACTER(len=1), INTENT(in) :: string(16)
    TYPE(t_uuid), INTENT(out) :: uuid
    uuid%data  = TRANSFER(string, uuid%data)
  END SUBROUTINE char2uuid

  SUBROUTINE uuid_generate(val, nval, uuid)
    REAL(C_DOUBLE), INTENT(IN)  :: val(*)
    INTEGER(C_INT), INTENT(IN)  :: nval
    TYPE(t_uuid),   INTENT(out) :: uuid
    CALL my_uuid_generate(val, nval, uuid)
  END SUBROUTINE uuid_generate

  INTEGER(C_INT) FUNCTION compare_uuid(uuid_A, uuid_B) 
    TYPE(t_uuid),   INTENT(IN)  :: uuid_A, uuid_B    
    REAL(C_DOUBLE) :: min_difference
    compare_uuid = my_compare_uuid(uuid_A, uuid_B, min_difference) 
  END FUNCTION compare_uuid

  SUBROUTINE clear_uuid(uuid)
    TYPE(t_uuid), INTENT(inout) :: uuid
    uuid%data(:) = INT(0, C_SIGNED_CHAR)
  END SUBROUTINE clear_uuid

END MODULE mo_util_uuid
