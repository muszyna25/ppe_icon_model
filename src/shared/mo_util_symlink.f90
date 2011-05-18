MODULE mo_util_symlink

  USE, INTRINSIC ::  ISO_C_BINDING, ONLY: C_INT, C_CHAR, C_NULL_CHAR
  
  IMPLICIT NONE
  
  PRIVATE

  INTERFACE 
    FUNCTION private_symlink(file, link) RESULT(iret) BIND(C,NAME='symlink')
#ifdef __SX__
      USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT, C_CHAR
#else
      IMPORT :: C_INT, C_CHAR
#endif
      INTEGER(C_INT) :: iret
      CHARACTER(C_CHAR), DIMENSION(*), INTENT(in) :: file
      CHARACTER(C_CHAR), DIMENSION(*), INTENT(in) :: link
    END FUNCTION private_symlink
  END INTERFACE
  
  INTERFACE
    FUNCTION private_unlink(filename) RESULT(iret) BIND(C,NAME='unlink')
#ifdef __SX__
      USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT, C_CHAR
#else
      IMPORT :: C_INT, C_CHAR
#endif
      INTEGER(C_INT) :: iret
      CHARACTER(C_CHAR), DIMENSION(*), INTENT(in) :: filename
    END FUNCTION private_unlink
  END INTERFACE
    
  INTERFACE
    FUNCTION private_islink(filename) RESULT(iret) BIND(C,NAME='util_islink')
#ifdef __SX__
      USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT, C_CHAR
#else
      IMPORT :: C_INT, C_CHAR
#endif
      INTEGER(C_INT) :: iret
      CHARACTER(C_CHAR), DIMENSION(*), INTENT(in) :: filename
    END FUNCTION private_islink
  END INTERFACE

  INTERFACE 
    FUNCTION private_rename(old_filename, new_filename) RESULT(iret) BIND(C,NAME='rename')
#ifdef __SX__
      USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT, C_CHAR
#else
      IMPORT :: C_INT, C_CHAR
#endif
      INTEGER(C_INT) :: iret
      CHARACTER(C_CHAR), DIMENSION(*), INTENT(in) :: old_filename
      CHARACTER(C_CHAR), DIMENSION(*), INTENT(in) :: new_filename
    END FUNCTION private_rename
  END INTERFACE
  
  PUBLIC :: util_symlink
  PUBLIC :: util_unlink
  PUBLIC :: util_islink
  PUBLIC :: util_rename

CONTAINS

  FUNCTION util_symlink(file, link) RESULT(iret)
    INTEGER :: iret
    CHARACTER(len=*), INTENT(in) :: file
    CHARACTER(len=*), INTENT(in) :: link
    iret = private_symlink(TRIM(file)//C_NULL_CHAR, TRIM(link)//C_NULL_CHAR)
  END FUNCTION util_symlink

  FUNCTION util_unlink(filename) RESULT(iret)
    INTEGER :: iret
    CHARACTER(len=*), INTENT(in) :: filename
    iret = private_unlink(TRIM(filename)//C_NULL_CHAR)
  END FUNCTION util_unlink

  FUNCTION util_islink(filename) RESULT(islink)
    LOGICAL :: islink
    CHARACTER(len=*), INTENT(in) :: filename
    INTEGER :: iret
    iret = private_islink(TRIM(filename)//C_NULL_CHAR)
    islink = .FALSE.
    IF (iret == 1) islink = .TRUE.
  END FUNCTION util_islink

  FUNCTION util_rename(old_filename, new_filename) RESULT(iret)
    INTEGER :: iret
    CHARACTER(len=*), INTENT(in) :: old_filename
    CHARACTER(len=*), INTENT(in) :: new_filename
    iret = private_rename(TRIM(old_filename)//C_NULL_CHAR, TRIM(new_filename)//C_NULL_CHAR)
  END FUNCTION util_rename
    
END MODULE mo_util_symlink


