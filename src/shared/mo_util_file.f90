MODULE mo_util_file

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
#ifdef __SX__
      CHARACTER(kind=C_CHAR,len=*), INTENT(in) :: file
      CHARACTER(kind=C_CHAR,len=*), INTENT(in) :: link
#else
      CHARACTER(C_CHAR), DIMENSION(*), INTENT(in) :: file
      CHARACTER(C_CHAR), DIMENSION(*), INTENT(in) :: link
#endif
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
#ifdef __SX__
     CHARACTER(kind=C_CHAR,len=*), INTENT(in) :: filename
#else
      CHARACTER(C_CHAR), DIMENSION(*), INTENT(in) :: filename
#endif
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
#ifdef __SX__
      CHARACTER(kind=C_CHAR,len=*), INTENT(in) :: filename
#else
      CHARACTER(C_CHAR), DIMENSION(*), INTENT(in) :: filename
#endif
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
#ifdef __SX__
      CHARACTER(kind=C_CHAR,len=*), INTENT(in) :: old_filename
      CHARACTER(kind=C_CHAR,len=*), INTENT(in) :: new_filename
#else
      CHARACTER(C_CHAR), DIMENSION(*), INTENT(in) :: old_filename
      CHARACTER(C_CHAR), DIMENSION(*), INTENT(in) :: new_filename
#endif
    END FUNCTION private_rename
  END INTERFACE

  INTERFACE
    FUNCTION private_tmpnam_len() RESULT(maxlen) BIND(C,NAME='util_tmpnam_len')
#ifdef __SX__
      USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT
#else
      IMPORT :: C_INT
#endif
      INTEGER(C_INT) :: maxlen
    END FUNCTION private_tmpnam_len
  END INTERFACE

  INTERFACE
    FUNCTION private_tmpnam(filename) RESULT(flen) BIND(C,NAME='util_tmpnam')
#ifdef __SX__
      USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_CHAR, C_INT
#else
      IMPORT :: C_CHAR, C_INT
#endif
      INTEGER(C_INT) :: flen
      CHARACTER(C_CHAR), DIMENSION(*), INTENT(inout) :: filename
    END FUNCTION private_tmpnam
  END INTERFACE

  INTERFACE
    FUNCTION private_filesize(filename) RESULT(flen) BIND(C,NAME='util_filesize')
#ifdef __SX__
      USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT, C_CHAR
#else
      IMPORT :: C_INT, C_CHAR
#endif
      INTEGER(C_INT) :: flen
#ifdef __SX__
      CHARACTER(kind=C_CHAR,len=*), INTENT(in) :: filename
#else
      CHARACTER(C_CHAR), DIMENSION(*), INTENT(in) :: filename
#endif
    END FUNCTION private_filesize
  END INTERFACE

  PUBLIC :: util_symlink
  PUBLIC :: util_unlink
  PUBLIC :: util_islink
  PUBLIC :: util_rename
  PUBLIC :: util_tmpnam
  PUBLIC :: util_filesize

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
    
  FUNCTION util_tmpnam(filename, klen) RESULT(flen)
    INTEGER :: flen
#ifdef __SX__
    CHARACTER(len=*), INTENT(out) :: filename
#else
    CHARACTER, DIMENSION(*), INTENT(out) :: filename
#endif
    INTEGER,                 INTENT(in)  :: klen
#ifdef __SX__
    INTEGER :: i
#endif
    !
    CHARACTER(C_CHAR), ALLOCATABLE :: tf(:)    
    INTEGER :: maxlen
    !
    maxlen = private_tmpnam_len()
    ALLOCATE(tf(maxlen))
    flen = private_tmpnam(tf)
    IF (flen > klen) THEN
      flen = -1
    ELSE
#ifdef __SX__
      DO i = 1, flen
        filename(i:i) = tf(i)
      ENDDO
#else
      filename(1:flen) = tf(1:flen)
#endif
    ENDIF
    DEALLOCATE(tf)
  END FUNCTION util_tmpnam

  FUNCTION util_filesize(filename) RESULT(flen)
    INTEGER :: flen
    CHARACTER(len=*), INTENT(in) :: filename
    flen = private_filesize(TRIM(filename)//C_NULL_CHAR)
  END FUNCTION util_filesize

END MODULE mo_util_file


