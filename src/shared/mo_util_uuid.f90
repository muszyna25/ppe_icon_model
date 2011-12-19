#ifdef __SX__
MODULE my_mo_util_uuid_type

  USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_CHAR

  IMPLICIT NONE 

  PUBLIC

  TYPE, BIND(C) :: t_uuid
    CHARACTER(C_CHAR) :: data(16)
  END type t_uuid

END MODULE my_mo_util_uuid_type
#endif

MODULE mo_util_uuid

#ifdef __SX__
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_CHAR, C_NULL_CHAR
  USE my_mo_util_uuid_type, ONLY: t_uuid
#else
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_CHAR, C_SIGNED_CHAR, C_NULL_CHAR
#endif

  IMPLICIT NONE 

  PRIVATE

  PUBLIC :: t_uuid

  PUBLIC :: uuid_generate
  PUBLIC :: uuid_parse
  PUBLIC :: uuid_unparse

  PUBLIC :: OPERATOR(==)

  PUBLIC :: uuid_string_length

  INTEGER, PARAMETER :: uuid_string_length = 36

#ifndef __SX__
  TYPE, BIND(C) :: t_uuid
    INTEGER(C_SIGNED_CHAR) :: data(16)
  END type t_uuid
#endif

  INTERFACE OPERATOR (==)
    MODULE PROCEDURE uuid_compare
  END INTERFACE OPERATOR (==)

  INTERFACE
    SUBROUTINE my_uuid_get(uuid) BIND(C,NAME='uuid_get')
#ifdef __SX__
      USE my_mo_util_uuid_type, ONLY: t_uuid
#else
      IMPORT :: t_uuid
#endif
      TYPE(t_uuid),     INTENT(out) :: uuid
    END SUBROUTINE my_uuid_get
  END INTERFACE

  INTERFACE
    SUBROUTINE my_uuid_format(uuid_string, uuid) BIND(C,NAME='uuid_format')
#ifdef __SX__
      USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_CHAR
      USE my_mo_util_uuid_type, ONLY: t_uuid
      CHARACTER(kind=C_CHAR,len=*),    INTENT(out) :: uuid_string
#else
      IMPORT :: C_CHAR, t_uuid
      CHARACTER(C_CHAR), DIMENSION(*), INTENT(out) :: uuid_string
#endif
      TYPE(t_uuid),                    INTENT(in)  :: uuid
    END SUBROUTINE my_uuid_format
  END INTERFACE

  INTERFACE
    SUBROUTINE my_uuid_parse(uuid, uuid_string) BIND(C,NAME='uuid_parse')
#ifdef __SX__
      USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_CHAR
      USE my_mo_util_uuid_type, ONLY: t_uuid
#else
      IMPORT :: C_CHAR, t_uuid
#endif
      TYPE(t_uuid),                    INTENT(out) :: uuid
#ifdef __SX__
      CHARACTER(kind=C_CHAR, len=*),   INTENT(in)  :: uuid_string
#else
      CHARACTER(C_CHAR), DIMENSION(*), INTENT(in)  :: uuid_string
#endif
    END SUBROUTINE my_uuid_parse
  END INTERFACE

CONTAINS

  SUBROUTINE uuid_generate(uuid)
    TYPE(t_uuid), INTENT(out) :: uuid
    CALL my_uuid_get(uuid)
  END SUBROUTINE uuid_generate
  
  SUBROUTINE uuid_parse(uuid_string, uuid)
    CHARACTER(len=*), INTENT(in)  :: uuid_string
    TYPE(t_uuid),     INTENT(out) :: uuid
    CALL my_uuid_parse(uuid, TRIM(uuid_string)//C_NULL_CHAR)
  END SUBROUTINE uuid_parse

  SUBROUTINE uuid_unparse(uuid, uuid_string)
    TYPE(t_uuid),     INTENT(in)  :: uuid
    CHARACTER(len=*), INTENT(out) :: uuid_string
    CALL my_uuid_format(uuid_string, uuid)
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
  
END MODULE mo_util_uuid
