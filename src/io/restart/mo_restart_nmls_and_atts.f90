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

  IMPLICIT NONE
  PRIVATE

  ! used throughout the icon code
  PUBLIC :: open_tmpfile, store_and_close_namelist
  PUBLIC :: open_and_restore_namelist, close_tmpfile

  PUBLIC :: restartAttributeList_make
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
    CHARACTER(*), PARAMETER :: routine = modname//":restartAttributeList_make"

    IF (.NOT.isRestart() .AND. .NOT.ocean_initFromRestart_OVERRIDE) &
      & CALL finish(routine, "not a restart run")
    IF (ASSOCIATED(gAttributeStore)) &
      & CALL finish(routine, "no second assignment of gAttributeStore allowed")
    ALLOCATE(gAttributeStore)
    IF (p_comm_rank(comm) == root_pe) CALL gAttributeStore%input(vlistId=vlistId)
    CALL gAttributeStore%bcast(root_pe, comm)
    IF (.NOT.ocean_initFromRestart_OVERRIDE) &
      CALL gAttributeStore%output(archive=gNamelistStore)
  END SUBROUTINE restartAttributeList_read

  SUBROUTINE restartAttributeList_make(ptr)
    TYPE(t_key_value_store), POINTER, INTENT(OUT) :: ptr

    ALLOCATE(ptr)
    CALL gNamelistStore%output(copy=ptr)
  END SUBROUTINE restartAttributeList_make

  SUBROUTINE bcastNamelistStore(root, comm)
    INTEGER, INTENT(IN) :: root, comm
    
    CALL gNamelistStore%bcast(root, comm)
  END SUBROUTINE bcastNamelistStore

END MODULE mo_restart_nml_and_att
