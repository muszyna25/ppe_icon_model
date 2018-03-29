!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
MODULE mo_util_file

  USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT, C_CHAR, C_NULL_CHAR, C_LONG, C_SIZE_T
  USE mo_exception, ONLY: finish
  USE mo_impl_constants, ONLY: SUCCESS
  USE mo_kind, ONLY: i8

  IMPLICIT NONE

  PRIVATE

  INTERFACE 
    FUNCTION private_symlink(file, link) RESULT(iret) BIND(C,NAME='symlink')
#if defined (__SUNPRO_F95)
      USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT, C_CHAR
#else
      IMPORT :: C_INT, C_CHAR
#endif
      INTEGER(C_INT) :: iret
#if defined (__SUNPRO_F95)
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
#if defined (__SUNPRO_F95)
      USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT, C_CHAR
#else
      IMPORT :: C_INT, C_CHAR
#endif
      INTEGER(C_INT) :: iret
#if defined (__SUNPRO_F95)
     CHARACTER(kind=C_CHAR,len=*), INTENT(in) :: filename
#else
      CHARACTER(C_CHAR), DIMENSION(*), INTENT(in) :: filename
#endif
    END FUNCTION private_unlink
  END INTERFACE
    
  INTERFACE
    FUNCTION private_islink(filename) RESULT(iret) BIND(C,NAME='util_islink')
#if defined (__SUNPRO_F95)
      USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT, C_CHAR
#else
      IMPORT :: C_INT, C_CHAR
#endif
      INTEGER(C_INT) :: iret
#if defined (__SUNPRO_F95)
      CHARACTER(kind=C_CHAR,len=*), INTENT(in) :: filename
#else
      CHARACTER(C_CHAR), DIMENSION(*), INTENT(in) :: filename
#endif
    END FUNCTION private_islink
  END INTERFACE

  INTERFACE 
    FUNCTION private_rename(old_filename, new_filename) RESULT(iret) BIND(C,NAME='rename')
#if defined (__SUNPRO_F95)
      USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT, C_CHAR
#else
      IMPORT :: C_INT, C_CHAR
#endif
      INTEGER(C_INT) :: iret
#if defined (__SUNPRO_F95)
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
#if defined (__SUNPRO_F95)
      USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT
#else
      IMPORT :: C_INT
#endif
      INTEGER(C_INT) :: maxlen
    END FUNCTION private_tmpnam_len
  END INTERFACE

  INTERFACE
    FUNCTION private_tmpnam(filename) RESULT(flen) BIND(C,NAME='util_tmpnam')
#if defined (__SUNPRO_F95)
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
#if defined (__SUNPRO_F95)
      USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_LONG, C_CHAR
#else
      IMPORT :: C_LONG, C_CHAR
#endif
      INTEGER(C_LONG) :: flen
#if defined (__SUNPRO_F95)
      CHARACTER(kind=C_CHAR,len=*), INTENT(in) :: filename
#else
      CHARACTER(C_CHAR), DIMENSION(*), INTENT(in) :: filename
#endif
    END FUNCTION private_filesize
  END INTERFACE

  INTERFACE
    FUNCTION private_file_is_writable(filename) RESULT(iwritable) BIND(C,NAME='util_file_is_writable')
#if defined (__SUNPRO_F95)
      USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT, C_CHAR
#else
      IMPORT :: C_INT, C_CHAR
#endif
      INTEGER(C_INT) :: iwritable
#if defined (__SUNPRO_F95)
      CHARACTER(kind=C_CHAR,len=*), INTENT(in) :: filename
#else
      CHARACTER(C_CHAR), DIMENSION(*), INTENT(in) :: filename
#endif
    END FUNCTION private_file_is_writable
  END INTERFACE

  PUBLIC :: util_symlink
  PUBLIC :: util_unlink
  PUBLIC :: util_islink
  PUBLIC :: util_rename
  PUBLIC :: util_tmpnam
  PUBLIC :: util_filesize
  PUBLIC :: util_file_is_writable
  PUBLIC :: putFile
  PUBLIC :: createSymlink

  CHARACTER(*), PARAMETER :: modname = "mo_util_file"

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
    
  FUNCTION generate_tmpnam(filename, klen) RESULT(flen)
    INTEGER :: flen
    CHARACTER(len=*), INTENT(out) :: filename
    INTEGER,          INTENT(in)  :: klen
    ! local variables
    INTEGER :: i
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
      DO i = 1, flen
        filename(i:i) = tf(i)
      ENDDO
    ENDIF
    DEALLOCATE(tf)
  END FUNCTION generate_tmpnam


  FUNCTION util_tmpnam(filename, klen) RESULT(flen)
    INTEGER :: flen
    CHARACTER(len=*), INTENT(out) :: filename
    INTEGER,          INTENT(in)  :: klen
    ! local variables
    INTEGER, PARAMETER :: N_RETRIES = 50
    INTEGER              :: i
    LOGICAL              :: lexists
    CHARACTER (LEN=klen) :: new_filename

    ! Note: (At least) on the SX-9 it is not sufficient to generate a
    ! filename - the TMPDIR of the local file system is seldom tidied
    ! up. Therefore, we test N_RETRIES times, if the generated
    ! filename already exists.
    !
    ! try to find a file name for our temporary file that does not
    ! exist yet:
    TEST_LOOP : DO i=1,N_RETRIES
      flen   = generate_tmpnam(new_filename, klen)
      INQUIRE(file=new_filename(1:flen), exist=lexists)
      IF (.NOT. lexists) THEN
        filename(1:flen) = new_filename(1:flen)
        EXIT TEST_LOOP
      END IF
      IF (i == N_RETRIES) THEN
        WRITE (0,*) "mo_util_file::util_tmpnam : Failed to find a tmp filename!"
        STOP
      END IF
    END DO TEST_LOOP
  END FUNCTION util_tmpnam

  FUNCTION util_filesize(filename) RESULT(flen)
    INTEGER(KIND=i8) :: flen
    CHARACTER(len=*), INTENT(in) :: filename
    flen = private_filesize(TRIM(filename)//C_NULL_CHAR)
  END FUNCTION util_filesize

  FUNCTION util_file_is_writable(filename) RESULT(lwritable)
    LOGICAL :: lwritable
    CHARACTER(len=*), INTENT(in) :: filename
    lwritable = (private_file_is_writable(TRIM(filename)//C_NULL_CHAR) == 1)
  END FUNCTION util_file_is_writable

  INTEGER FUNCTION putFile(path, string, fileMode) RESULT(resultVar)
    CHARACTER(LEN = *), INTENT(IN) :: path, string
    INTEGER(C_INT), VALUE :: fileMode

    INTEGER :: error, i, pathLen, stringLen
    CHARACTER(KIND = C_CHAR) :: pathCopy(LEN(path) + 1) ! path IS passed zero terminated
    CHARACTER(KIND = C_CHAR), ALLOCATABLE :: stringCopy(:)
    CHARACTER(*), PARAMETER :: routine = modname//":putFile"

    INTERFACE
        INTEGER(C_INT) FUNCTION c_putFile(c_path, c_dataSize, c_data, c_fileMode) BIND(C, NAME = "putFile")
            IMPORT C_INT, C_CHAR, C_SIZE_T
            CHARACTER(KIND = C_CHAR) :: c_path(*), c_data(*)
            INTEGER(C_SIZE_T), VALUE :: c_dataSize
            INTEGER(C_INT), VALUE :: c_fileMode
        END FUNCTION c_putFile
    END INTERFACE

    pathLen = LEN(path)
    stringLen = LEN(string)

    ALLOCATE(stringCopy(stringLen), STAT = error) ! string IS passed with an explicit SIZE argument, so no need for a termination CHARACTER
    IF(error /= SUCCESS) CALL finish(routine, "memory allocation failure")

    DO i = 1, pathLen
        pathCopy(i) = path(i:i)
    END DO
    pathCopy(pathLen + 1) = C_NULL_CHAR

    DO i = 1, stringLen
        stringCopy(i) = string(i:i)
    END DO

    resultVar = c_putFile(pathCopy, INT(stringLen, C_SIZE_T), stringCopy, fileMode)
  END FUNCTION putFile

  INTEGER FUNCTION createSymlink(targetPath, linkName) RESULT(error)
    CHARACTER(*), INTENT(IN) :: targetPath, linkName

    INTEGER :: i
    CHARACTER(KIND = C_CHAR) :: linkNameCopy(LEN(linkName) + 1), targetPathCopy(LEN(targetPath) + 1)

    INTERFACE
        INTEGER(C_INT) FUNCTION c_createSymlink(c_targetPath, c_linkName) BIND(C, NAME = "createSymlink")
            IMPORT C_INT, C_CHAR
            CHARACTER(KIND = C_CHAR) :: c_targetPath(*), c_linkName(*)
        END FUNCTION c_createSymlink
    END INTERFACE

    DO i = 1, LEN(targetPath)
        targetPathCopy(i) = targetPath(i:i)
    END DO
    targetPathCopy(i) = C_NULL_CHAR

    DO i = 1, LEN(linkName)
        linkNameCopy(i) = linkName(i:i)
    END DO
    linkNameCopy(i) = C_NULL_CHAR

    error = c_createSymlink(targetPathCopy, linkNameCopy)
  END FUNCTION createSymlink

END MODULE mo_util_file


