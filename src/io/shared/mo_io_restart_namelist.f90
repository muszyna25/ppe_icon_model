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
  USE mo_io_units,    ONLY: nerr, find_next_free_unit, filename_max
  USE mo_exception,   ONLY: message, finish
  USE mo_mpi,         ONLY: p_bcast, p_pe, p_comm_rank
  USE mo_cdi,         ONLY: CDI_GLOBAL, CDI_UNDEFID, vlistInqNatts, vlistInqAtt, vlistInqAttTxt

  IMPLICIT NONE

  PRIVATE

  ! used throughout the icon code
  PUBLIC :: open_tmpfile
  PUBLIC :: store_and_close_namelist
  PUBLIC :: open_and_restore_namelist
  PUBLIC :: close_tmpfile

  ! used ONLY by mo_io_restart AND mo_io_restart_async
  PUBLIC :: read_and_bcast_restart_namelists
  PUBLIC :: set_restart_namelist
  PUBLIC :: get_restart_namelist
  PUBLIC :: print_restart_name_lists
  PUBLIC :: delete_restart_namelists    ! also used by mo_atmo_model

  ! variables used by mo_io_restart AND mo_io_restart_async
  PUBLIC :: restart_namelist
  PUBLIC :: nmls


  PUBLIC :: t_att_namelist


  PUBLIC :: nmllen_max
!   INTEGER, PARAMETER :: nmllen_max = 4096
  INTEGER, PARAMETER :: nmllen_max = 65536

  TYPE t_att_namelist
    CHARACTER(len=64) :: name
    CHARACTER(KIND = C_CHAR), POINTER :: text(:)
  END TYPE t_att_namelist

  INTEGER, SAVE :: nmls = 0
  TYPE(t_att_namelist), POINTER, SAVE :: restart_namelist(:)
  LOGICAL, SAVE :: l_restart_namelist_initialized = .FALSE.

  INTERFACE get_restart_namelist
    MODULE PROCEDURE get_restart_namelist_by_name
    MODULE PROCEDURE get_restart_namelist_by_index
  END INTERFACE get_restart_namelist

  CHARACTER(LEN = *), PARAMETER :: modname = "mo_io_restart_namelist"

CONTAINS

  SUBROUTINE init_restart_namelists()
    INTEGER :: i, error
    CHARACTER(LEN = *), PARAMETER :: routine = modname//":init_restart_namelists"

    IF(l_restart_namelist_initialized) CALL finish(routine, "assertion failed: init_restart_namelists() called several times")
    IF(nmls /= 0) CALL finish(routine, "assertion failed")

    ALLOCATE(restart_namelist(8), STAT = error)
    IF(error /= SUCCESS) CALL finish(routine, "memory allocation failed")
    DO i = 1, SIZE(restart_namelist, 1)
        restart_namelist(i)%text => NULL()
    END DO
    l_restart_namelist_initialized = .TRUE.
  END SUBROUTINE init_restart_namelists

  FUNCTION find_namelist(namelist_name) RESULT(RESULT)
    CHARACTER(len=*), INTENT(in) :: namelist_name
    TYPE(t_att_namelist), POINTER :: RESULT

    CHARACTER(LEN = *), PARAMETER :: routine = modname//":find_namelist"
    INTEGER :: i, error
    TYPE(t_att_namelist), POINTER :: temp(:)

    ! first ensure that we actually have a list
    IF(.NOT. l_restart_namelist_initialized) CALL init_restart_namelists()

    ! check whether we already have an entry of that NAME
    DO i = 1, nmls
        IF('nml_'//TRIM(namelist_name) == TRIM(restart_namelist(i)%NAME)) THEN
            RESULT => restart_namelist(i)
            RETURN
        END IF
    END DO

    ! looks like we have to create an entry IF this point IS reached
    ! ensure that we have at least enough space for one more entry
    IF(nmls == SIZE(restart_namelist, 1)) THEN
        ALLOCATE(temp(2*nmls), STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation failed")
        DO i = 1, nmls
            temp(i)%NAME = restart_namelist(i)%NAME
            temp(i)%text => restart_namelist(i)%text
            restart_namelist(i)%text => NULL()
        END DO
        DO i = nmls + 1, SIZE(temp, 1)
            temp(i)%text => NULL()
        END DO
        DEALLOCATE(restart_namelist)
        restart_namelist => temp
    END IF

    ! add an entry
    nmls = nmls + 1
    RESULT => restart_namelist(nmls)
    RESULT%NAME = 'nml_'//TRIM(namelist_name)
  END FUNCTION find_namelist

  SUBROUTINE delete_restart_namelists()
    INTEGER :: i

    IF (l_restart_namelist_initialized) THEN
        DO i = 1, nmls
            IF(ASSOCIATED(restart_namelist(i)%text)) DEALLOCATE(restart_namelist(i)%text)
        END DO
        DEALLOCATE(restart_namelist)
        l_restart_namelist_initialized = .FALSE.
    ENDIF
    nmls = 0
  END SUBROUTINE delete_restart_namelists

  SUBROUTINE set_restart_namelist(namelist_name, namelist_text)
    CHARACTER(len=*), INTENT(in) :: namelist_name
    CHARACTER(len=*), INTENT(in) :: namelist_text
    TYPE(t_att_namelist), POINTER :: list_entry

    list_entry => find_namelist(namelist_name)
    IF(ASSOCIATED(list_entry%text)) DEALLOCATE(list_entry%text)
    list_entry%text => toCharArray(namelist_text)
  END SUBROUTINE set_restart_namelist

  SUBROUTINE get_restart_namelist_by_name(namelist_name, namelist_text)
    CHARACTER(len=*),              INTENT(in)  :: namelist_name
    CHARACTER(len=:), ALLOCATABLE, INTENT(out) :: namelist_text

    CHARACTER(LEN = *), PARAMETER :: routine = modname//":get_restart_namelist_by_name"
    TYPE(t_att_namelist), POINTER :: list_entry
    CHARACTER(LEN = :), POINTER :: tempText
    INTEGER :: textLength, error

    list_entry => find_namelist(namelist_name)
    IF(.NOT.ASSOCIATED(list_entry%text)) CALL finish(routine, 'namelist '//TRIM(namelist_name)//' not available in restart file.')
    tempText => toCharacter(list_entry%text)
    textLength = LEN(tempText)

    IF(ALLOCATED(namelist_text)) DEALLOCATE(namelist_text)
    ALLOCATE(CHARACTER(LEN = textLength) :: namelist_text, STAT = error)
    IF(error /= SUCCESS) CALL finish(routine, "memory allocation failed")

    namelist_text(1:textLength) = tempText(1:textLength)
    DEALLOCATE(tempText)
  END SUBROUTINE get_restart_namelist_by_name

  SUBROUTINE get_restart_namelist_by_index(namelist_index, namelist_name)
    INTEGER,          INTENT(in)  :: namelist_index
    CHARACTER(len=*), INTENT(out) :: namelist_name

    IF (namelist_index > nmls) CALL finish('get_restart_namelist_by_index','index out of range.')
    namelist_name = TRIM(restart_namelist(namelist_index)%name(5:))
  END SUBROUTINE get_restart_namelist_by_index

  SUBROUTINE print_restart_name_lists()
    INTEGER                           :: i
    CHARACTER(LEN=*), PARAMETER       :: routine = modname//':print_restart_name_list'

    WRITE(nerr, '(a,a,i3)') routine, ' p_pe=', p_pe
    PRINT *,'restart name lists count = ',nmls
    DO i = 1, nmls
        PRINT *, ' restart name list = ', TRIM(restart_namelist(i)%NAME), ' text = "', restart_namelist(i)%text, '"'
    ENDDO
  END SUBROUTINE print_restart_name_lists

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

    CHARACTER(LEN = :), ALLOCATABLE :: nmlbuf

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

  SUBROUTINE read_and_bcast_restart_namelists(vlistID, root_pe, comm)
    INTEGER, INTENT(IN) :: vlistID      !< CDI vlist ID
    INTEGER, INTENT(IN) :: root_pe      !< rank of the PE that reads the vlist AND broadcast the namelist information
    INTEGER, INTENT(IN) :: comm         !< MPI communicator

    ! local variables
    CHARACTER(LEN = *), PARAMETER :: routine = modname//":read_and_bcast_restart_namelists"
    INTEGER :: natts, att_type, att_len, i, status, my_rank
    CHARACTER(len=64) :: att_name
    CHARACTER(LEN = :), ALLOCATABLE :: tempText
    TYPE(t_att_namelist), POINTER :: list_entry

    ! reset the list
    CALL delete_restart_namelists()

    ! get the number of attributes so we can loop over them
    my_rank = p_comm_rank(comm)
    IF (my_rank == root_pe) THEN
      IF(vlistID == CDI_UNDEFID) CALL finish(routine, "assertion failed: root PE does not have a valid vlistID")
      status = vlistInqNatts(vlistID, CDI_GLOBAL, natts)
      IF(status /= SUCCESS) CALL finish(routine, "vlistInqNatts() returned an error")
    END IF
    CALL p_bcast(natts, root_pe, comm)

    DO i = 0, natts-1
        ! inquire the attribute NAME AND check whether it IS a namelist attribute
        att_name = ''
        att_len  = 0
        IF (my_rank == root_pe) THEN
            status = vlistInqAtt(vlistID, CDI_GLOBAL, i, att_name, att_type, att_len)
            IF(status /= SUCCESS) CALL finish(routine, "vlistInqAtt() returned an error")
        END IF
        CALL p_bcast(att_name, root_pe, comm)
        IF ( att_name(1:4) /= 'nml_') CYCLE ! skip this, it is not a namelist 

        ! it IS a namelist attribute, so we add it to the list
        list_entry => find_namelist(att_name(5:))
        CALL p_bcast(att_len, root_pe, comm)
        ALLOCATE(CHARACTER(len=att_len) :: tempText)
        IF (my_rank == root_pe) THEN
            status = vlistInqAttTxt(vlistID, CDI_GLOBAL, TRIM(att_name), att_len, tempText)
            IF(status /= SUCCESS) CALL finish(routine, "vlistInqAttTxt() returned an error")
        END IF
        CALL p_bcast(tempText, root_pe, comm)
        list_entry%text => toCharArray(tempText)
        DEALLOCATE(tempText)
    ENDDO

  END SUBROUTINE read_and_bcast_restart_namelists

END MODULE mo_io_restart_namelist
