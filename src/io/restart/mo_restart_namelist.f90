!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This module contains utility functions for storing the namelist parameters in the restart files.

MODULE mo_restart_namelist
    USE mo_util_file, ONLY: util_tmpnam, util_filesize, util_unlink
    USE mo_util_string, ONLY: tocompact
    USE mo_impl_constants, ONLY: SUCCESS
    USE mo_io_units, ONLY: nerr, find_next_free_unit, filename_max
    USE mo_exception, ONLY: message, finish
    USE mo_mpi, ONLY: p_bcast, p_pe, p_comm_rank, p_get_bcast_role
    USE mo_cdi, ONLY: CDI_GLOBAL, CDI_UNDEFID, CDI_MAX_NAME, vlistInqNatts, vlistInqAtt, vlistInqAttTxt, vlistDefAttTxt

    IMPLICIT NONE

    PRIVATE

    ! used throughout the icon code
    PUBLIC :: open_tmpfile
    PUBLIC :: store_and_close_namelist
    PUBLIC :: open_and_restore_namelist
    PUBLIC :: close_tmpfile

    ! used ONLY by mo_sync_restart AND mo_async_restart
    PUBLIC :: RestartNamelist_bcast
    PUBLIC :: read_and_bcast_restart_namelists
    PUBLIC :: RestartNamelist_writeToFile
    PUBLIC :: set_restart_namelist
    PUBLIC :: get_restart_namelist
    PUBLIC :: print_restart_name_lists
    PUBLIC :: delete_restart_namelists ! also used by mo_atmo_model

    PUBLIC :: t_att_namelist

    TYPE t_att_namelist
        CHARACTER(:), ALLOCATABLE :: name, text
    END TYPE t_att_namelist

    TYPE t_NamelistArchive
        INTEGER :: namelistCount
        TYPE(t_att_namelist), POINTER :: namelists(:)
    CONTAINS
        PROCEDURE :: setNamelist => namelistArchive_setNamelist
        PROCEDURE :: getNamelist => namelistArchive_getNamelist
        PROCEDURE :: print => namelistArchive_print
        PROCEDURE :: writeToFile => namelistArchive_writeToFile ! store the namelists as attributes to the given CDI vlistId

        PROCEDURE, PRIVATE :: construct => namelistArchive_construct
        PROCEDURE, PRIVATE :: find => namelistArchive_find  ! returns a new entry IF NONE exists yet
        PROCEDURE, PRIVATE :: destruct => namelistArchive_destruct
    END TYPE t_NamelistArchive

    TYPE(t_NamelistArchive), ALLOCATABLE, SAVE :: gArchive

    CHARACTER(LEN = *), PARAMETER :: modname = "mo_restart_namelist"

CONTAINS

    SUBROUTINE namelistArchive_construct(me)
        CLASS(t_NamelistArchive), INTENT(INOUT) :: me

        INTEGER :: error
        CHARACTER(*), PARAMETER :: routine = modname//":namelistArchive_construct"

        me%namelistCount = 0
        ALLOCATE(me%namelists(8), STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation failed")
    END SUBROUTINE namelistArchive_construct

    FUNCTION namelistArchive_find(me, namelist_name) RESULT(RESULT)
        CLASS(t_NamelistArchive), INTENT(INOUT) :: me
        CHARACTER(len=*), INTENT(in) :: namelist_name
        TYPE(t_att_namelist), POINTER :: RESULT

        INTEGER :: i, error
        TYPE(t_att_namelist), POINTER :: temp(:)
        CHARACTER(:), ALLOCATABLE :: fullName
        CHARACTER(LEN = *), PARAMETER :: routine = modname//":namelistArchive_find"

        ! check whether we already have an entry of that NAME
        fullName = 'nml_'//TRIM(namelist_name)
        DO i = 1, me%namelistCount
            IF(fullName == me%namelists(i)%NAME) THEN
                RESULT => me%namelists(i)
                DEALLOCATE(fullName)
                RETURN
            END IF
        END DO

        ! looks like we have to create an entry IF this point IS reached
        ! ensure that we have at least enough space for one more entry
        IF(me%namelistCount == SIZE(me%namelists, 1)) THEN
            ALLOCATE(temp(2*me%namelistCount), STAT = error)
            IF(error /= SUCCESS) CALL finish(routine, "memory allocation failed")
            DO i = 1, me%namelistCount
                CALL MOVE_ALLOC(me%namelists(i)%NAME, temp(i)%NAME)
                CALL MOVE_ALLOC(me%namelists(i)%text, temp(i)%text)
            END DO
            DEALLOCATE(me%namelists)
            me%namelists => temp
        END IF

        ! add an entry
        me%namelistCount = me%namelistCount + 1
        RESULT => me%namelists(me%namelistCount)
        CALL MOVE_ALLOC(fullName, RESULT%NAME)
    END FUNCTION namelistArchive_find

    SUBROUTINE namelistArchive_setNamelist(me, namelistName, namelistText)
        CLASS(t_NamelistArchive), INTENT(INOUT) :: me
        CHARACTER(*), INTENT(IN) :: namelistName, namelistText

        TYPE(t_att_namelist), POINTER :: list_entry
        INTEGER :: textLength, error
        CHARACTER(*), PARAMETER :: routine = modname//":namelistArchive_setNamelist"

        list_entry => me%find(namelistName)
        textLength = LEN_TRIM(namelistText)

        IF(ALLOCATED(list_entry%text)) DEALLOCATE(list_entry%text)
        ALLOCATE(CHARACTER(textLength) :: list_entry%text, STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation error")

        list_entry%text(1:textLength) = namelistText(1:textLength)
    END SUBROUTINE namelistArchive_setNamelist

    SUBROUTINE namelistArchive_getNamelist(me, namelistName, namelistText)
        CLASS(t_NamelistArchive), INTENT(INOUT) :: me
        CHARACTER(len=*), INTENT(in) :: namelistName
        CHARACTER(len=:), ALLOCATABLE, INTENT(out) :: namelistText

        TYPE(t_att_namelist), POINTER :: list_entry
        CHARACTER(LEN = *), PARAMETER :: routine = modname//":namelistArchive_getNamelist"

        list_entry => me%find(namelistName)
        IF(.NOT.ALLOCATED(list_entry%text)) CALL finish(routine, 'namelist '//TRIM(namelistName)//' not available in restart file.')
        namelistText = list_entry%text
    END SUBROUTINE namelistArchive_getNamelist

    SUBROUTINE namelistArchive_print(me)
        CLASS(t_NamelistArchive), INTENT(INOUT) :: me

        INTEGER :: i
        CHARACTER(LEN=*), PARAMETER :: routine = modname//':namelistArchive_print'

        WRITE(nerr, '(a,a,i3)') routine, ' p_pe=', p_pe
        PRINT *,'restart name lists count = ',me%namelistCount
        DO i = 1, me%namelistCount
            PRINT *, ' restart name list = "'//me%namelists(i)%NAME//'", text = "'//me%namelists(i)%text//'"'
        ENDDO
    END SUBROUTINE namelistArchive_print

    SUBROUTINE namelistArchive_writeToFile(me, cdiVlistId)
        CLASS(t_NamelistArchive), INTENT(INOUT) :: me
        INTEGER, VALUE :: cdiVlistId

        INTEGER :: i, error
        CHARACTER(LEN = *), PARAMETER :: routine = modname//":namelistArchive_writeToFile"

        DO i = 1, me%namelistCount
            error = vlistDefAttTxt(cdiVlistId, CDI_GLOBAL, me%namelists(i)%name, LEN(me%namelists(i)%text), &
                                  &me%namelists(i)%text)
            IF(error /= SUCCESS) CALL finish(routine, "error WHILE writing a namelist to a restart file")
        END DO
    END SUBROUTINE namelistArchive_writeToFile

    SUBROUTINE namelistArchive_destruct(me)
        CLASS(t_NamelistArchive), INTENT(INOUT) :: me

        INTEGER :: i

        DO i = 1, me%namelistCount
            IF(ALLOCATED(me%namelists(i)%NAME)) DEALLOCATE(me%namelists(i)%NAME)
            IF(ALLOCATED(me%namelists(i)%text)) DEALLOCATE(me%namelists(i)%text)
        END DO
        DEALLOCATE(me%namelists)
    END SUBROUTINE namelistArchive_destruct

    SUBROUTINE init_restart_namelists()
        INTEGER :: i, error
        CHARACTER(LEN = *), PARAMETER :: routine = modname//":init_restart_namelists"

        IF(ALLOCATED(gArchive)) CALL finish(routine, "assertion failed: init_restart_namelists() called several times")

        ALLOCATE(gArchive, STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation failed")
        CALL gArchive%construct()
    END SUBROUTINE init_restart_namelists

    SUBROUTINE delete_restart_namelists()
        IF (ALLOCATED(gArchive)) THEN
            CALL gArchive%destruct()
            DEALLOCATE(gArchive)
        ENDIF
    END SUBROUTINE delete_restart_namelists

    SUBROUTINE set_restart_namelist(namelist_name, namelist_text)
        CHARACTER(len=*), INTENT(in) :: namelist_name
        CHARACTER(len=*), INTENT(in) :: namelist_text

        IF(.NOT. ALLOCATED(gArchive)) CALL init_restart_namelists()
        CALL gArchive%setNamelist(namelist_name, namelist_text)
    END SUBROUTINE set_restart_namelist

    SUBROUTINE get_restart_namelist(namelist_name, namelist_text)
        CHARACTER(len=*), INTENT(in) :: namelist_name
        CHARACTER(len=:), ALLOCATABLE, INTENT(out) :: namelist_text

        IF(.NOT. ALLOCATED(gArchive)) CALL init_restart_namelists()
        CALL gArchive%getNamelist(namelist_name, namelist_text)
    END SUBROUTINE get_restart_namelist

    SUBROUTINE print_restart_name_lists()
        IF(.NOT.ALLOCATED(gArchive)) CALL init_restart_namelists()
        CALL gArchive%print()
    END SUBROUTINE print_restart_name_lists

    SUBROUTINE RestartNamelist_writeToFile(cdiVlistId)
        INTEGER, VALUE :: cdiVlistId

        IF(.NOT.ALLOCATED(gArchive)) CALL init_restart_namelists()
        CALL gArchive%writeToFile(cdiVlistId)
    END SUBROUTINE RestartNamelist_writeToFile

    FUNCTION open_tmpfile() RESULT(funit)
        INTEGER :: funit

        INTEGER :: flen
        CHARACTER(len=filename_max) :: filename

        flen = util_tmpnam(filename, filename_max)
        funit = find_next_free_unit(10,100)
        OPEN(UNIT=funit, FILE=filename(1:flen), &
             ACTION='write', &
             ACCESS='sequential', &
             DELIM='apostrophe')
    END FUNCTION open_tmpfile

    SUBROUTINE store_and_close_namelist(funit, name)
        INTEGER, INTENT(in) :: funit
        CHARACTER(len=*), INTENT(in) :: name

        CHARACTER(len=filename_max) :: filename
        INTEGER :: nmllen, error
        CHARACTER(len = :), ALLOCATABLE :: nmlbuf
        CHARACTER(LEN = *), PARAMETER :: routine = modname//":store_and_close_namelist"

        INQUIRE(funit, NAME=filename)

        CLOSE(funit)

        nmllen = util_filesize(filename)
        IF (nmllen == 0) THEN
            CALL message(routine, 'namelist '//TRIM(name)//' is empty, saving in restart file fails.')
        ENDIF
        ALLOCATE(CHARACTER(LEN = nmllen) :: nmlbuf, STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation failed")
        nmlbuf(1:nmllen) = ''

#ifdef __SX__
        ! requires in runscript (ksh/bash): export F_NORCW=65535
        ! (this SX environment variable specifies that a control record is
        ! not added/expected, s.t. the file content can be treated like a
        ! stream of characters)
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

        error = util_unlink(TRIM(filename))
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
        ! not added/expected, s.t. the file content can be treated like a
        ! stream of characters)
        OPEN(UNIT=65535, FILE=filename(1:flen), ACTION='write', FORM='unformatted')
        WRITE(65535) TRIM(nmlbuf)
        CLOSE(65535)
#else
        OPEN(UNIT=funit, FILE=filename(1:flen), ACTION='write', ACCESS='stream', FORM='unformatted')
        WRITE(funit) TRIM(nmlbuf)
        CLOSE(funit)
#endif

        OPEN(UNIT=funit, FILE=filename(1:flen), &
             ACTION='read', &
             ACCESS='sequential', &
             RECL=65535, &
             DELIM='apostrophe')
    END FUNCTION open_and_restore_namelist

    SUBROUTINE close_tmpfile(funit)
        INTEGER, INTENT(in) :: funit

        CHARACTER(len=filename_max) :: filename
        INTEGER :: error
        CHARACTER(*), PARAMETER :: routine = modname//":close_tmpfile"

        INQUIRE(funit,NAME=filename)
        CLOSE(funit)
        error = util_unlink(TRIM(filename))
        IF(error /= SUCCESS) CALL message(routine, "warning: error while unlinking file '"//TRIM(filename)//"'")
    END SUBROUTINE close_tmpfile

    SUBROUTINE read_and_bcast_restart_namelists(vlistID, root_pe, comm)
        INTEGER, INTENT(IN) :: vlistID, root_pe, comm

        INTEGER :: natts, att_type, att_len, i, status
        CHARACTER(LEN = CDI_MAX_NAME) :: att_name
        CHARACTER(LEN = :), ALLOCATABLE :: tempText
        CHARACTER(LEN = *), PARAMETER :: routine = modname//":read_and_bcast_restart_namelists"

        ! reset the list
        CALL delete_restart_namelists()

        ! the root pe reads the namelists from the file
        IF (p_comm_rank(comm) == root_pe) THEN
            IF(vlistID == CDI_UNDEFID) CALL finish(routine, "assertion failed: root PE does not have a valid vlistID")

            ! get the number of attributes so we can loop over them
            status = vlistInqNatts(vlistID, CDI_GLOBAL, natts)
            IF(status /= SUCCESS) CALL finish(routine, "vlistInqNatts() returned an error")

            DO i = 0, natts-1
                ! inquire the attribute NAME AND check whether it IS a namelist attribute
                att_name = ''
                att_len = 0
                status = vlistInqAtt(vlistID, CDI_GLOBAL, i, att_name, att_type, att_len)
                IF(status /= SUCCESS) CALL finish(routine, "vlistInqAtt() returned an error")

                IF(att_name(1:4) /= 'nml_') CYCLE ! skip this, it is not a namelist

                ! it IS a namelist attribute, so we add it to the list
                ALLOCATE(CHARACTER(len=att_len) :: tempText)
                status = vlistInqAttTxt(vlistID, CDI_GLOBAL, TRIM(att_name), att_len, tempText)
                IF(status /= SUCCESS) CALL finish(routine, "vlistInqAttTxt() returned an error")
                CALL set_restart_namelist(att_name(5:), tempText)
                DEALLOCATE(tempText)
            ENDDO
        END IF

        ! then we broadcast the RESULT
        CALL RestartNamelist_bcast(root_pe, comm)
    END SUBROUTINE read_and_bcast_restart_namelists

    !-------------------------------------------------------------------------------------------------
    !
    ! Transfers the restart name lists from the worker to the restart PEs.
    !
    ! Designed to work with both intra AND inter communicators.
    SUBROUTINE RestartNamelist_bcast(root, communicator)
        INTEGER, VALUE :: root, communicator

        LOGICAL :: lIsRoot, lIsReceiver
        INTEGER :: iv, nv, i, maxNameLength, maxAttributeLength, error
        CHARACTER(LEN = :), ALLOCATABLE :: list_name, list_text
        CHARACTER(LEN = *), PARAMETER :: routine = modname//":RestartNamelist_bcast"

        ! get our role IN the broadcast operation
        CALL p_get_bcast_role(root, communicator, lIsRoot, lIsReceiver)

        ! delete old name lists
        IF(lIsReceiver) CALL delete_restart_namelists

        ! get the number of name lists
        IF(lIsRoot) nv = gArchive%namelistCount
        CALL p_bcast(nv, root, communicator)

        maxNameLength = 1
        maxAttributeLength = 1
        IF(lIsRoot) THEN
            DO iv = 1, nv
                maxNameLength = MAX(maxNameLength, LEN(gArchive%namelists(iv)%NAME))
                maxAttributeLength = MAX(maxAttributeLength, LEN(gArchive%namelists(iv)%text))
            END DO
        END IF
        CALL p_bcast(maxNameLength, root, communicator)
        CALL p_bcast(maxAttributeLength, root, communicator)

        ALLOCATE(CHARACTER(LEN = maxNameLength) :: list_name, STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation failed")
        ALLOCATE(CHARACTER(LEN = maxAttributeLength) :: list_text, STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation failed")

        DO iv = 1, nv
            ! send name of the name list
            list_name(1:maxNameLength) = ''
            IF(lIsRoot) THEN
                ! The first four characters are skipped because that's the `nml_` which will be added again by set_restart_namelist().
                list_name(1:LEN(gArchive%namelists(iv)%NAME) - 4) = gArchive%namelists(iv)%NAME(5:)
            END IF
            CALL p_bcast(list_name(1:maxNameLength), root, communicator)

            ! send text of the name list
            list_text(1:maxAttributeLength) = ''
            IF(lIsRoot) THEN
                list_text(1:maxAttributeLength) = gArchive%namelists(iv)%text
            END IF
            CALL p_bcast(list_text(1:maxAttributeLength), root, communicator)

            ! store name list parameters
            IF(lIsReceiver) THEN
                CALL set_restart_namelist(list_name(1:maxNameLength), list_text(1:maxAttributeLength))
            ENDIF
        ENDDO

        DEALLOCATE(list_name)
        DEALLOCATE(list_text)
    END SUBROUTINE RestartNamelist_bcast

END MODULE mo_restart_namelist
