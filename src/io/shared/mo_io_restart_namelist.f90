!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
MODULE mo_io_restart_namelist

  USE ISO_C_BINDING, ONLY: C_CHAR
  USE mo_util_file,   ONLY: util_tmpnam, util_filesize, util_unlink
  USE mo_util_string, ONLY: tocompact, toCharacter, toCharArray
  USE mo_impl_constants, ONLY: SUCCESS
  USE mo_io_units,    ONLY: find_next_free_unit, filename_max
  USE mo_exception,   ONLY: message, finish
  USE mo_mpi,         ONLY: p_bcast
  USE mo_cdi,         ONLY: CDI_GLOBAL, vlistInqNatts, vlistInqAtt, vlistInqAttTxt

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: open_tmpfile
  PUBLIC :: store_and_close_namelist
  PUBLIC :: open_and_restore_namelist
  PUBLIC :: close_tmpfile
  PUBLIC :: read_and_bcast_restart_namelists
  PUBLIC :: get_restart_namelist
  PUBLIC :: delete_restart_namelists
  PUBLIC :: nmls
  PUBLIC :: restart_namelist
  PUBLIC :: set_restart_namelist

  PUBLIC :: t_att_namelist


  PUBLIC :: nmllen_max
!   INTEGER, PARAMETER :: nmllen_max = 4096
  INTEGER, PARAMETER :: nmllen_max = 65536

  TYPE t_att_namelist
    CHARACTER(len=64) :: name
    CHARACTER(KIND = C_CHAR), POINTER :: text(:)
  END TYPE t_att_namelist

  INTEGER, PARAMETER :: nmax_nmls = 64
  INTEGER, SAVE :: nmls = 0
  TYPE(t_att_namelist), ALLOCATABLE :: restart_namelist(:)

  INTERFACE get_restart_namelist
    MODULE PROCEDURE get_restart_namelist_by_name
    MODULE PROCEDURE get_restart_namelist_by_index
  END INTERFACE get_restart_namelist

  CHARACTER(LEN = *), PARAMETER :: modname = "mo_io_restart_namelist"

CONTAINS

  SUBROUTINE init_restart_namelists()
    INTEGER :: i, error
    CHARACTER(LEN = *), PARAMETER :: routine = modname//":init_restart_namelists"

    IF(ALLOCATED(restart_namelist)) CALL finish(routine, "assertion failed: init_restart_namelists() called several times")
    IF(nmls /= 0) CALL finish(routine, "assertion failed")

    ALLOCATE(restart_namelist(nmax_nmls), STAT = error)
    IF(error /= SUCCESS) CALL finish(routine, "memory allocation failed")
    DO i = 1, SIZE(restart_namelist, 1)
        restart_namelist(i)%text => NULL()
    END DO
  END SUBROUTINE init_restart_namelists

  SUBROUTINE delete_restart_namelists()
    INTEGER :: i

    IF (ALLOCATED(restart_namelist)) THEN
        DO i = 1, nmls
            IF(ASSOCIATED(restart_namelist(i)%text)) DEALLOCATE(restart_namelist(i)%text)
        END DO
        DEALLOCATE(restart_namelist)
    ENDIF
    nmls = 0
  END SUBROUTINE delete_restart_namelists

  SUBROUTINE set_restart_namelist(namelist_name, namelist_text)
    CHARACTER(len=*), INTENT(in) :: namelist_name
    CHARACTER(len=*), INTENT(in) :: namelist_text
    INTEGER :: i

    IF (.NOT. ALLOCATED(restart_namelist)) THEN
      CALL init_restart_namelists()
    ENDIF

    ! look, if entry has to be overwritten

    DO i = 1, nmls
      IF ('nml_'//TRIM(namelist_name) == TRIM(restart_namelist(i)%name)) THEN
        restart_namelist(i)%text => toCharArray(TRIM(namelist_text))
        RETURN
      ENDIF
    ENDDO

    ! ok, we have a new entry

    nmls = nmls+1
    IF (nmls > nmax_nmls) THEN
      CALL finish('set_restart_namelist', &
           &      'too many restart attributes for restart file')
    ELSE
      restart_namelist(nmls)%name = 'nml_'//TRIM(namelist_name)
      restart_namelist(nmls)%text => toCharArray(TRIM(namelist_text))
    ENDIF

  END SUBROUTINE set_restart_namelist

  SUBROUTINE get_restart_namelist_by_name(namelist_name, namelist_text)
    CHARACTER(len=*),              INTENT(in)  :: namelist_name
    CHARACTER(len=*),              INTENT(out) :: namelist_text
    INTEGER :: i
    CHARACTER(LEN = :), POINTER :: tempText

    namelist_text = ''

    DO i = 1, nmls
      IF ('nml_'//TRIM(namelist_name) == TRIM(restart_namelist(i)%name)) THEN
        tempText => toCharacter(restart_namelist(i)%text)
        namelist_text = tempText
        DEALLOCATE(tempText)
        RETURN
      ENDIF
    ENDDO

    CALL finish('','namelist '//TRIM(namelist_name)//' not available in restart file.')

  END SUBROUTINE get_restart_namelist_by_name

  SUBROUTINE get_restart_namelist_by_index(namelist_index, namelist_name)
    INTEGER,          INTENT(in)  :: namelist_index
    CHARACTER(len=*), INTENT(out) :: namelist_name

    IF (nmls > nmax_nmls) THEN
      CALL finish('get_restart_namelist','index out of range.')
    ELSE
      namelist_name = ''
      namelist_name = TRIM(restart_namelist(namelist_index)%name(5:))
    ENDIF

  END SUBROUTINE get_restart_namelist_by_index

  FUNCTION open_tmpfile() RESULT(funit)
    INTEGER :: funit

    INTEGER :: flen
    CHARACTER(len=filename_max) :: filename

    flen = util_tmpnam(filename, filename_max)
    funit = find_next_free_unit(10,100)
    OPEN(UNIT=funit, FILE=filename(1:flen), &
         ACTION='write',                    &
         ACCESS='sequential',               &
         DELIM='apostrophe')

  END FUNCTION open_tmpfile

  SUBROUTINE store_and_close_namelist(funit, name)
    INTEGER,          INTENT(in) :: funit
    CHARACTER(len=*), INTENT(in) :: name 

    CHARACTER(len=filename_max) :: filename
    INTEGER :: nmllen
    CHARACTER(len=nmllen_max) :: nmlbuf
    INTEGER :: iret

    INQUIRE(funit, NAME=filename)

    CLOSE(funit)

    nmllen = util_filesize(filename)
    IF (nmllen == 0) THEN
      CALL message('','namelist '//TRIM(name)//' is empty, saving in restart file fails.')
    ENDIF
    IF (nmllen > nmllen_max) THEN
      CALL message('', &
           'The problem could be solved by increasing nmllen_max in mo_io_restart_namelist.f90.')
      CALL finish('','namelist '//TRIM(name)//' is too long, saving in restart file fails.')
    ENDIF
    nmlbuf = ''

#ifdef __SX__
    ! requires in runscript (ksh/bash): export F_NORCW=65535
    ! (this SX environment variable specifies that a control record is
    !  not added/expected, s.t. the file content can be treated like a
    !  stream of characters)
    OPEN(UNIT=65535, FILE=TRIM(filename), ACTION='read', FORM='unformatted')
    READ(65535) nmlbuf(1:nmllen)
    CLOSE(65535)
#else
    OPEN(UNIT=funit, FILE=TRIM(filename), ACTION='read', ACCESS='stream', FORM='unformatted')
    READ(funit) nmlbuf(1:nmllen)
    CLOSE(funit)
#endif

    CALL tocompact(nmlbuf)
    CALL set_restart_namelist(TRIM(name), TRIM(nmlbuf))


    iret = util_unlink(TRIM(filename))

  END SUBROUTINE store_and_close_namelist

  FUNCTION open_and_restore_namelist(name) RESULT(funit)
    INTEGER :: funit
    CHARACTER(len=*), INTENT(in) :: name 

    INTEGER :: flen
    CHARACTER(len=filename_max) :: filename

    CHARACTER(len=nmllen_max) :: nmlbuf

    CALL get_restart_namelist(name, nmlbuf)

    funit = find_next_free_unit(10,100)
    flen = util_tmpnam(filename, filename_max)
#ifdef __SX__
    ! requires in runscript (ksh/bash): export F_NORCW=65535
    ! (this SX environment variable specifies that a control record is
    !  not added/expected, s.t. the file content can be treated like a
    !  stream of characters)
    OPEN(UNIT=65535, FILE=filename(1:flen), ACTION='write', FORM='unformatted')
    WRITE(65535) TRIM(nmlbuf)
    CLOSE(65535)
#else
    OPEN(UNIT=funit, FILE=filename(1:flen), ACTION='write', ACCESS='stream', FORM='unformatted')
    WRITE(funit) TRIM(nmlbuf)
    CLOSE(funit)
#endif


    OPEN(UNIT=funit, FILE=filename(1:flen), &
         ACTION='read',                     &
         ACCESS='sequential',               &
         RECL=65535,                        &
         DELIM='apostrophe')

  END FUNCTION open_and_restore_namelist

  SUBROUTINE close_tmpfile(funit)
    INTEGER, INTENT(in) :: funit
    CHARACTER(len=filename_max) :: filename
    INTEGER :: iret

    INQUIRE(funit,NAME=filename)

    CLOSE(funit)

    iret = util_unlink(TRIM(filename))

  END SUBROUTINE close_tmpfile

  SUBROUTINE read_and_bcast_restart_namelists(vlistID, lread_pe, root_pe, comm)
    INTEGER, INTENT(IN) :: vlistID      !< CDI vlist ID
    LOGICAL, INTENT(IN) :: lread_pe     !< .TRUE., if current PE has opened the file for reading
    INTEGER, INTENT(IN) :: root_pe      !< rank of broadcast root PE
    INTEGER, INTENT(IN) :: comm         !< MPI communicator

    ! local variables
    INTEGER             :: natts, att_type, att_len, nmllen, &
      &                    i, status
    CHARACTER(len=64)   :: att_name
    CHARACTER(len=nmllen_max) :: nmlbuf

    ! Reset the namelists.
    IF (ALLOCATED(restart_namelist)) THEN
      CALL delete_restart_namelists()
    ENDIF
    CALL init_restart_namelists()

    IF (lread_pe) THEN
      status = vlistInqNatts(vlistID, CDI_GLOBAL, natts)
    END IF
    CALL p_bcast(natts, root_pe, comm)

    DO i = 0, natts-1
      att_name = ''
      att_len  = 0
      IF (lread_pe) THEN
        status = vlistInqAtt(vlistID, CDI_GLOBAL, i, att_name, att_type, att_len)
      END IF
      CALL p_bcast(att_name, root_pe, comm)
      IF ( att_name(1:4) /= 'nml_') CYCLE ! skip this, it is not a namelist 

      CALL p_bcast(att_len, root_pe, comm)
      nmllen = att_len
      IF (nmllen > nmllen_max) THEN
        CALL message('', &
             &       'The problem could be solved by increasing nmllen_max '// &
             &       'in mo_io_restart_namelist.f90.')
        CALL finish('','namelist '//TRIM(att_name)// &
             &      ' is too long, reload from restart file fails.')
      ENDIF
      IF (lread_pe) THEN
        status = vlistInqAttTxt(vlistID, CDI_GLOBAL, TRIM(att_name), nmllen, nmlbuf)
      END IF
      CALL p_bcast(nmlbuf, root_pe, comm)

      nmls = nmls+1
      IF (nmls > nmax_nmls) THEN
        CALL finish('set_restart_namelist', &
             &      'too many restart attributes for restart file')
      ENDIF
      restart_namelist(nmls)%name = TRIM(att_name)
      restart_namelist(nmls)%text => toCharArray(nmlbuf(1:nmllen))
    ENDDO

  END SUBROUTINE read_and_bcast_restart_namelists

END MODULE mo_io_restart_namelist
