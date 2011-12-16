MODULE mo_util_uuid

  USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_CHAR, C_SIGNED_CHAR, C_NULL_CHAR

  IMPLICIT NONE 

  PRIVATE

  PUBLIC :: t_uuid

  PUBLIC :: uuid_generate
  PUBLIC :: uuid_parse
  PUBLIC :: uuid_unparse

  PUBLIC :: OPERATOR(==)

  PUBLIC :: uuid_string_length

  INTEGER, PARAMETER :: uuid_string_length = 36

  TYPE, BIND(C) :: t_uuid
    INTEGER(C_SIGNED_CHAR) :: data(16)
  END type t_uuid

  INTERFACE OPERATOR (==)
    MODULE PROCEDURE uuid_compare
  END INTERFACE OPERATOR (==)

  INTERFACE
    SUBROUTINE my_uuid_get(uuid) BIND(C,NAME='uuid_get')
      IMPORT :: t_uuid
      TYPE(t_uuid),     INTENT(out) :: uuid
    END SUBROUTINE my_uuid_get
  END INTERFACE

  INTERFACE
    SUBROUTINE my_uuid_format(uuid_string, uuid) BIND(C,NAME='uuid_format')
      IMPORT :: C_CHAR, t_uuid
      CHARACTER(C_CHAR), DIMENSION(*), INTENT(out) :: uuid_string
      TYPE(t_uuid),                    INTENT(in)  :: uuid
    END SUBROUTINE my_uuid_format
  END INTERFACE

  INTERFACE
    SUBROUTINE my_uuid_parse(uuid, uuid_string) BIND(C,NAME='uuid_parse')
      IMPORT :: C_CHAR, t_uuid
      TYPE(t_uuid),                    INTENT(out) :: uuid
      CHARACTER(C_CHAR), DIMENSION(*), INTENT(in)  :: uuid_string
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
    IF (ANY(uuid1%data /= uuid2%data)) THEN
      uuid_compare = .FALSE.
    ENDIF
  END FUNCTION uuid_compare
  
END MODULE mo_util_uuid
