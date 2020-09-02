!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
MODULE mo_restart_nml_and_att

  USE mo_exception,             ONLY: finish, message
  USE mo_impl_constants,        ONLY: SUCCESS
  USE mo_master_config,         ONLY: isRestart
  USE mo_mpi, ONLY: p_comm_rank
  USE mo_util_file, ONLY: util_tmpnam, util_filesize, util_unlink
  USE mo_util_string, ONLY: tocompact
  USE mo_io_units, ONLY: filename_max, find_next_free_unit
  USE mo_key_value_store, ONLY: t_key_value_store
  USE ISO_C_BINDING, ONLY: C_DOUBLE, C_INT
  USE mo_cdi, ONLY: DATATYPE_FLT64, DATATYPE_INT32, DATATYPE_TXT, CDI_GLOBAL, &
    & cdiInqNatts, cdiInqAtt, cdiInqAttFlt, cdiInqAttInt, cdiInqAttTxt, &
    & CDI_UNDEFID, CDI_MAX_NAME, cdidefAttInt, cdidefAttFlt, cdidefAttTxt
  USE mo_kind, ONLY: wp
  USE mo_hash_table, ONLY: t_HashIterator

  IMPLICIT NONE
  PRIVATE

  ! used throughout the icon code
  PUBLIC :: open_tmpfile, store_and_close_namelist
  PUBLIC :: open_and_restore_namelist, close_tmpfile

  PUBLIC :: restartAttributeList_make, restartAttributeList_write_to_cdi
  PUBLIC :: getAttributesForRestarting, restartAttributeList_read

  PUBLIC :: bcastNamelistStore

  TYPE(t_key_value_store), POINTER :: gAttributeStore => NULL()
  TYPE(t_key_value_store) :: gNamelistStore

  CHARACTER(LEN = *), PARAMETER :: modname = "mo_restart_nml_and_att"
  LOGICAL, SAVE, PUBLIC :: ocean_initFromRestart_OVERRIDE = .FALSE. 

CONTAINS

  INTEGER FUNCTION open_tmpfile() RESULT(funit)
    INTEGER :: flen
    CHARACTER(len=filename_max) :: filename

    flen = util_tmpnam(filename, filename_max)
    funit = find_next_free_unit(10,100)
    OPEN(UNIT=funit, FILE=filename(1:flen), ACTION='write', ACCESS='sequential', DELIM='apostrophe')
  END FUNCTION open_tmpfile

  SUBROUTINE store_and_close_namelist(funit, nmlname)
    INTEGER, INTENT(IN) :: funit
    CHARACTER(*), INTENT(in) :: nmlname
    CHARACTER(filename_max) :: filename
    INTEGER :: nmllen, ierr
    CHARACTER(:), ALLOCATABLE :: nmlbuf
    CHARACTER(LEN = *), PARAMETER :: routine = modname//":store_and_close_namelist"

    INQUIRE(funit, NAME=filename)
    CLOSE(funit)
    nmllen = INT(util_filesize(filename))
    IF (nmllen == 0) CALL message(routine, 'empty nml "'//TRIM(nmlname)//'" : not saved')
    ALLOCATE(CHARACTER(LEN = nmllen) :: nmlbuf, STAT=ierr)
    IF(ierr /= SUCCESS) CALL finish(routine, "memory allocation failed")
    nmlbuf(1:nmllen) = ''
    OPEN(UNIT=funit, FILE=TRIM(filename), ACTION='read', ACCESS='stream', FORM='unformatted')
    READ(funit) nmlbuf(1:nmllen)
    CLOSE(funit)
    CALL tocompact(nmlbuf)
    IF(.NOT.gNamelistStore%is_init) CALL gNamelistStore%init(.true.)
    CALL gNamelistStore%put('nml_'//nmlname, TRIM(nmlbuf))
    ierr = util_unlink(TRIM(filename))
  END SUBROUTINE store_and_close_namelist

  INTEGER FUNCTION open_and_restore_namelist(nmlname) RESULT(funit)
    CHARACTER(len=*), INTENT(in) :: nmlname
    INTEGER :: flen
    CHARACTER(len=filename_max) :: filename
    CHARACTER(:), ALLOCATABLE :: nmlbuf

    IF (.NOT.gNamelistStore%is_init) &
      & CALL finish(modname//':open_and_restore_namelist', "gNamelistStore not initialized")
    CALL gNamelistStore%get('nml_'//nmlname, nmlbuf)
    flen = util_tmpnam(filename, filename_max)
    funit = find_next_free_unit(10, 100)
    OPEN(UNIT=funit, FILE=filename(1:flen), ACTION='write', ACCESS='stream', FORM='unformatted')
    WRITE(funit) TRIM(nmlbuf)
    CLOSE(funit)
    OPEN(UNIT=funit, FILE=filename(1:flen), ACTION='read', ACCESS='sequential', RECL=65535, DELIM='apostrophe')
  END FUNCTION open_and_restore_namelist

  SUBROUTINE close_tmpfile(funit)
    INTEGER, INTENT(in) :: funit
    CHARACTER(len=filename_max) :: filename
    INTEGER :: ierr
    CHARACTER(*), PARAMETER :: routine = modname//":close_tmpfile"

    INQUIRE(funit,NAME=filename)
    CLOSE(funit)
    ierr = util_unlink(TRIM(filename))
    IF(ierr /= SUCCESS) CALL message(routine, "warning: error unlinking '"//TRIM(filename)//"'")
  END SUBROUTINE close_tmpfile

  SUBROUTINE getAttributesForRestarting(ptr)
    TYPE(t_key_value_store), POINTER, INTENT(OUT) :: ptr
    CHARACTER(*), PARAMETER :: routine = modname//":getAttributesForRestarting"

    IF (.NOT.ASSOCIATED(gAttributeStore) .AND. &
      & (isRestart() .OR. ocean_initFromRestart_OVERRIDE)) &
      & CALL finish(routine, "restart attributes not yet set up")
    ptr => gAttributeStore
  END SUBROUTINE getAttributesForRestarting

  SUBROUTINE restartAttributeList_read(vlistId, root_pe, comm)
    INTEGER, INTENT(IN) :: vlistId, root_pe, comm
    CHARACTER(*), PARAMETER :: routine = modname//":restartAttributeList_read"

    IF (.NOT.isRestart() .AND. .NOT.ocean_initFromRestart_OVERRIDE) &
      & CALL finish(routine, "not a restart run")
    IF (ASSOCIATED(gAttributeStore)) &
      & CALL finish(routine, "no second assignment of gAttributeStore allowed")
    ALLOCATE(gAttributeStore)
    IF (p_comm_rank(comm) == root_pe) CALL read_from_cdi(gAttributeStore)
    CALL gAttributeStore%bcast(root_pe, comm)
    IF (.NOT.ocean_initFromRestart_OVERRIDE) &
      CALL gAttributeStore%output(archive=gNamelistStore)
  CONTAINS

  SUBROUTINE read_from_cdi(kvs)
    TYPE(t_key_value_store), INTENT(INOUT) :: kvs
    INTEGER :: natt, i, ierr, attType, attLen, namLen
    CHARACTER(CDI_MAX_NAME) :: attName
    REAL(C_DOUBLE) :: oneDouble(1)
    INTEGER(C_INT) :: oneInt(1)
    CHARACTER(:), ALLOCATABLE :: attTxt
    CHARACTER(*), PARAMETER :: routine = "read_from_cdi"

    CALL kvs%init(.true.)
    IF (vlistID == CDI_UNDEFID) CALL finish(routine, "invalid vlistId")
    ALLOCATE(CHARACTER(LEN=4096) :: attTxt, STAT=ierr)
    IF(ierr /= SUCCESS) CALL finish(routine, "memory allocation error")
    ierr = cdiInqNatts(vlistID, CDI_GLOBAL, natt)
    IF (ierr /= SUCCESS) CALL finish(routine, "getting attribute count failed")
    DO i = 0, natt - 1
      ierr = cdiInqAtt(vlistId, CDI_GLOBAL, i, attName, attType, attLen)
      IF (ierr /= SUCCESS) CALL finish(routine, "failed reading attribute name")
      namLen = LEN_TRIM(attName)
      SELECT CASE(attType)
      CASE(DATATYPE_FLT64)
        ierr = cdiInqAttFlt(vlistID, CDI_GLOBAL, attName(1:namLen), 1, oneDouble)
        IF (ierr /= SUCCESS) CALL finish(routine, "failed reading attribute '"//attName(1:namLen)//"'")
        CALL kvs%put(attName(1:namLen), REAL(oneDouble(1), wp))
      CASE(DATATYPE_INT32)
        ierr = cdiInqAttInt(vlistID, CDI_GLOBAL, attName(1:namLen), 1, oneInt)
        IF (ierr /= SUCCESS) CALL finish(routine, "failed reading attribute '"//attName(1:namLen)//"'")
        IF ('bool_' == attName(1:MIN(5,namLen))) THEN
          CALL kvs%put(attName(6:namLen), INT(oneInt(1)) .NE. 0)
        ELSE
          CALL kvs%put(attName(1:namLen), INT(oneInt(1)))
        END IF
      CASE(DATATYPE_TXT)
        IF (attLen .GT. LEN(attTxt)) DEALLOCATE(attTxt)
        IF (.NOT.ALLOCATED(attTxt)) THEN
          ALLOCATE(CHARACTER(LEN=attLen) :: attTxt, STAT=ierr)
          IF (ierr /= SUCCESS) CALL finish(routine, "memory allocation error")
        END IF
        ierr = cdiInqAttTxt(vlistID, CDI_GLOBAL, attName(1:namLen), attLen, attTxt(1:attLen))
        IF (ierr /= SUCCESS) CALL finish(routine, "failed reading attribute '"//attName(1:namLen)//"'")
        CALL kvs%put(attName(1:namLen), attTxt(1:attLen))
      END SELECT
    ENDDO
  END SUBROUTINE read_from_cdi

  END SUBROUTINE restartAttributeList_read

 SUBROUTINE restartAttributeList_write_to_cdi(kvs, vlistID)
#ifdef __PGI
    USE mo_key_value_store, ONLY: t_char_workaround
#endif
    TYPE(t_key_value_store), INTENT(IN) :: kvs
    INTEGER, INTENT(IN) :: vlistId
    CHARACTER(*), PARAMETER :: routine = modname//"write_to_cdi"
    TYPE(t_HashIterator), POINTER :: iterator
    CLASS(*), POINTER :: curKey, curVal
    CHARACTER(:), POINTER :: ccKey, ccVal
    INTEGER :: ierr

    iterator => kvs%getIterator()
    DO WHILE (iterator%nextEntry(curKey, curVal))
      SELECT TYPE(curKey)
#ifdef __PGI
      TYPE IS(t_char_workaround)
        ccKey => curKey%c
#endif
      TYPE IS(CHARACTER(*))
        ccKey => curKey
      CLASS DEFAULT
        CALL finish(routine, "key: invalid type")
      END SELECT
      SELECT TYPE(curVal)
#ifdef __PGI
      TYPE IS(t_char_workaround)
        ierr = cdidefAttTxt(vlistId, CDI_GLOBAL, ccKey, LEN(curVal%c), curVal%c)
#else
      TYPE IS(CHARACTER(*))
        ierr = cdidefAttTxt(vlistId, CDI_GLOBAL, ccKey, LEN(curVal), curVal)
#endif
      TYPE IS(REAL(wp))
        ierr = cdidefAttFlt(vlistId, CDI_GLOBAL, ccKey, DATATYPE_FLT64, 1, [REAL(curVal, C_DOUBLE)])
      TYPE IS(INTEGER)
        ierr = cdidefAttInt(vlistId, CDI_GLOBAL, ccKey, DATATYPE_INT32, 1, [INT(curVal, C_INT)])
      TYPE IS(LOGICAL)
        ierr = cdidefAttInt(vlistId, CDI_GLOBAL, 'bool_'//ccKey, DATATYPE_INT32, 1, [INT(MERGE(1,0,curVal), C_INT)])
      CLASS DEFAULT
        CALL finish(routine, "val: invalid type")
      END SELECT
      IF(ierr /= SUCCESS) CALL finish(routine, "error writing attribute")
    END DO
    DEALLOCATE(iterator)
  END SUBROUTINE restartAttributeList_write_to_cdi

  SUBROUTINE restartAttributeList_make(ptr)
    TYPE(t_key_value_store), ALLOCATABLE, INTENT(OUT) :: ptr

    ALLOCATE(ptr)
    CALL gNamelistStore%output(copy=ptr)
  END SUBROUTINE restartAttributeList_make

  SUBROUTINE bcastNamelistStore(root, comm)
    INTEGER, INTENT(IN) :: root, comm
    
    CALL gNamelistStore%bcast(root, comm)
  END SUBROUTINE bcastNamelistStore

END MODULE mo_restart_nml_and_att
