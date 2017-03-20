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
    USE mo_cdi, ONLY: CDI_GLOBAL, CDI_UNDEFID, CDI_MAX_NAME, cdiInqNatts, cdiInqAtt, cdiInqAttTxt, cdiDefAttTxt
    USE mo_exception, ONLY: message, finish
    USE mo_impl_constants, ONLY: SUCCESS
    USE mo_io_units, ONLY: nerr, find_next_free_unit, filename_max
    USE mo_mpi, ONLY: p_pe, p_get_bcast_role
    USE mo_packed_message, ONLY: t_PackedMessage, kPackOp, kUnpackOp
    USE mo_util_file, ONLY: util_tmpnam, util_filesize, util_unlink
    USE mo_util_string, ONLY: tocompact

    IMPLICIT NONE

    PRIVATE

    ! used throughout the icon code
    PUBLIC :: open_tmpfile
    PUBLIC :: store_and_close_namelist
    PUBLIC :: open_and_restore_namelist
    PUBLIC :: close_tmpfile

    ! used ONLY by mo_sync_restart AND mo_async_restart
    PUBLIC :: t_NamelistArchive
    PUBLIC :: namelistArchive

    TYPE t_Namelist
        CHARACTER(:), ALLOCATABLE :: name, text
    END TYPE t_Namelist

    TYPE t_NamelistArchive
        INTEGER :: namelistCount
        TYPE(t_Namelist), POINTER :: namelists(:)
    CONTAINS
        PROCEDURE :: setNamelist => namelistArchive_setNamelist
        PROCEDURE :: getNamelist => namelistArchive_getNamelist
        PROCEDURE :: print => namelistArchive_print

        ! store the namelists as attributes to the given CDI vlistId, noncollective:
        PROCEDURE :: writeToCdiVlist => namelistArchive_writeToCdiVlist

        PROCEDURE :: readFromFile => namelistArchive_readFromFile   ! noncollective!
        PROCEDURE :: packer => namelistArchive_packer
        PROCEDURE :: bcast => namelistArchive_bcast
        PROCEDURE :: reset => namelistArchive_reset

        ! singleton: namelistArchive() IS the ONLY point where a
        ! t_NamelistArchive IS constructed:
        PROCEDURE, PRIVATE :: construct => namelistArchive_construct

        PROCEDURE, PRIVATE :: find => namelistArchive_find  ! returns a new entry IF NONE exists yet
        ! no destructor since this IS a singleton that lives until the very END of the program run
    END TYPE t_NamelistArchive

    CHARACTER(LEN = *), PARAMETER :: modname = "mo_restart_namelist"

CONTAINS
    FUNCTION namelistArchive() RESULT(resultVar)
        TYPE(t_NamelistArchive), POINTER :: resultVar

        TYPE(t_NamelistArchive), ALLOCATABLE, TARGET, SAVE :: archive
        INTEGER :: error
        CHARACTER(*), PARAMETER :: routine = modname//":namelistArchive"

        IF(.NOT.ALLOCATED(archive)) THEN
            ALLOCATE(archive, STAT = error)
            IF(error /= SUCCESS) CALL finish(routine, "memory allocation error")
            CALL archive%construct()
        END IF
        resultVar => archive
    END FUNCTION namelistArchive

    SUBROUTINE namelistArchive_construct(me)
        CLASS(t_NamelistArchive), INTENT(INOUT) :: me

        INTEGER :: error
        CHARACTER(*), PARAMETER :: routine = modname//":namelistArchive_construct"

        me%namelistCount = 0
        ALLOCATE(me%namelists(8), STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation failed")
    END SUBROUTINE namelistArchive_construct

    FUNCTION namelistArchive_find(me, namelist_name) RESULT(resultVar)
        CLASS(t_NamelistArchive), INTENT(INOUT) :: me
        CHARACTER(len=*), INTENT(in) :: namelist_name
        TYPE(t_Namelist), POINTER :: resultVar

        INTEGER :: i, error
        TYPE(t_Namelist), POINTER :: temp(:)
        CHARACTER(:), ALLOCATABLE :: fullName
        CHARACTER(LEN = *), PARAMETER :: routine = modname//":namelistArchive_find"

        ! check whether we already have an entry of that NAME
        fullName = 'nml_'//TRIM(namelist_name)
        DO i = 1, me%namelistCount
            IF(fullName == me%namelists(i)%NAME) THEN
                resultVar => me%namelists(i)
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
                temp(i)%name = me%namelists(i)%name
                IF(temp(i)%name /= me%namelists(i)%name) CALL finish(routine, "assertion failed")
                temp(i)%text = me%namelists(i)%text
                IF(temp(i)%text /= me%namelists(i)%text) CALL finish(routine, "assertion failed")
            END DO
            DEALLOCATE(me%namelists)
            me%namelists => temp
        END IF

        ! add an entry
        me%namelistCount = me%namelistCount + 1
        resultVar => me%namelists(me%namelistCount)
        resultVar%NAME = fullName
        IF(resultVar%NAME /= fullName) CALL finish(routine, "assertion failed")
    END FUNCTION namelistArchive_find

    SUBROUTINE namelistArchive_setNamelist(me, namelistName, namelistText)
        CLASS(t_NamelistArchive), INTENT(INOUT) :: me
        CHARACTER(*), INTENT(IN) :: namelistName, namelistText

        TYPE(t_Namelist), POINTER :: list_entry
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

        TYPE(t_Namelist), POINTER :: list_entry
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

    SUBROUTINE namelistArchive_writeToCdiVlist(me, cdiVlistId)
        CLASS(t_NamelistArchive), INTENT(INOUT) :: me
        INTEGER, VALUE :: cdiVlistId

        INTEGER :: i, error
        CHARACTER(LEN = *), PARAMETER :: routine = modname//":namelistArchive_writeToCdiVlist"

        DO i = 1, me%namelistCount
            error = cdiDefAttTxt(cdiVlistId, CDI_GLOBAL, me%namelists(i)%name, LEN(me%namelists(i)%text), &
                                  &me%namelists(i)%text)
            IF(error /= SUCCESS) CALL finish(routine, "error WHILE writing a namelist to a restart file")
        END DO
    END SUBROUTINE namelistArchive_writeToCdiVlist

    SUBROUTINE namelistArchive_readFromFile(me, vlistID)
        CLASS(t_NamelistArchive), INTENT(INOUT) :: me
        INTEGER, VALUE :: vlistID

        INTEGER :: natts, att_type, att_len, i, status
        CHARACTER(LEN = CDI_MAX_NAME) :: att_name
        CHARACTER(LEN = :), ALLOCATABLE :: tempText
        CHARACTER(LEN = *), PARAMETER :: routine = modname//":namelistArchive_readFromFile"

        IF(vlistID == CDI_UNDEFID) CALL finish(routine, "assertion failed: invalid vlistId passed to "//routine//"()")

        ! reset the list
        CALL me%reset()

        ! get the number of attributes so we can loop over them
        status = cdiInqNatts(vlistID, CDI_GLOBAL, natts)
        IF(status /= SUCCESS) CALL finish(routine, "cdiInqNatts() returned an error")

        DO i = 0, natts-1
            ! inquire the attribute NAME AND check whether it IS a namelist attribute
            att_name = ''
            att_len = 0
            status = cdiInqAtt(vlistID, CDI_GLOBAL, i, att_name, att_type, att_len)
            IF(status /= SUCCESS) CALL finish(routine, "cdiInqAtt() returned an error")

            IF(att_name(1:4) /= 'nml_') CYCLE ! skip this, it is not a namelist

            ! it IS a namelist attribute, so we add it to the list
            ALLOCATE(CHARACTER(len=att_len) :: tempText)
            status = cdiInqAttTxt(vlistID, CDI_GLOBAL, TRIM(att_name), att_len, tempText)
            IF(status /= SUCCESS) CALL finish(routine, "cdiInqAttTxt() returned an error")
            CALL me%setNamelist(att_name(5:), tempText)
            DEALLOCATE(tempText)
        END DO
    END SUBROUTINE namelistArchive_readFromFile

    SUBROUTINE namelistArchive_packer(me, operation, packedMessage)
        CLASS(t_NamelistArchive), INTENT(INOUT) :: me
        INTEGER, VALUE :: operation
        TYPE(t_PackedMessage), INTENT(INOUT) :: packedMessage

        INTEGER :: allocSize, error, i, length
        CHARACTER(*), PARAMETER :: routine = modname//":namelistArchive_packer"

        IF(operation == kUnpackOp) CALL me%reset()  ! get rid of the old contents 

        CALL packedMessage%packer(operation, me%namelistCount)

        IF(operation == kUnpackOp) THEN
            ! ensure sufficient space for the contents of the message
            allocSize = MAX(8, me%namelistCount)
            IF(allocSize > SIZE(me%namelists)) THEN
                DEALLOCATE(me%namelists)
                ALLOCATE(me%namelists(allocSize), STAT = error)
                IF(error /= SUCCESS) CALL finish(routine, "memory allocation error")
            END IF
        END IF
        IF(me%namelistCount > SIZE(me%namelists)) CALL finish(routine, "assertion failed")

        ! (un)pack the payload contents
        DO i = 1, me%namelistCount
            length = 0
            IF(ALLOCATED(me%namelists(i)%NAME)) length = LEN(me%namelists(i)%NAME)
            CALL packedMessage%packer(operation, length)
            IF(operation == kUnpackOp) THEN
                IF(ALLOCATED(me%namelists(i)%NAME)) DEALLOCATE(me%namelists(i)%NAME)
                ALLOCATE(CHARACTER(LEN = length) :: me%namelists(i)%NAME, STAT = error)
                IF(error /= SUCCESS) CALL finish(routine, "memory allocation failure")
            END IF
            CALL packedMessage%packer(operation, me%namelists(i)%NAME)

            length = 0
            IF(ALLOCATED(me%namelists(i)%text)) length = LEN(me%namelists(i)%text)
            CALL packedMessage%packer(operation, length)
            IF(operation == kUnpackOp) THEN
                IF(ALLOCATED(me%namelists(i)%text)) DEALLOCATE(me%namelists(i)%text)
                ALLOCATE(CHARACTER(LEN = length) :: me%namelists(i)%text, STAT = error)
                IF(error /= SUCCESS) CALL finish(routine, "memory allocation failure")
            END IF
            CALL packedMessage%packer(operation, me%namelists(i)%text)
        END DO
    END SUBROUTINE namelistArchive_packer

    SUBROUTINE namelistArchive_bcast(me, root, communicator)
        CLASS(t_NamelistArchive), INTENT(INOUT) :: me
        INTEGER, VALUE :: root, communicator

        TYPE(t_PackedMessage) :: packedMessage
        LOGICAL :: lIsRoot, lIsReceiver

        CALL p_get_bcast_role(root, communicator, lIsRoot, lIsReceiver)

        CALL packedMessage%construct()
        IF(lIsRoot) CALL me%packer(kPackOp, packedMessage)
        CALL packedMessage%bcast(root, communicator)
        IF(lIsReceiver) CALL me%packer(kUnpackOp, packedMessage)
        CALL packedMessage%destruct()
    END SUBROUTINE namelistArchive_bcast

    SUBROUTINE namelistArchive_reset(me)
        CLASS(t_NamelistArchive), INTENT(INOUT) :: me

        INTEGER :: i

        DO i = 1, me%namelistCount
            IF(ALLOCATED(me%namelists(i)%NAME)) DEALLOCATE(me%namelists(i)%NAME)
            IF(ALLOCATED(me%namelists(i)%text)) DEALLOCATE(me%namelists(i)%text)
        END DO
    END SUBROUTINE namelistArchive_reset

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
        TYPE(t_NamelistArchive), POINTER :: archive
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

        OPEN(UNIT=funit, FILE=TRIM(filename), ACTION='read', ACCESS='stream', FORM='unformatted')
        READ(funit) nmlbuf(1:nmllen)
        CLOSE(funit)

        CALL tocompact(nmlbuf)
        archive => namelistArchive()
        CALL archive%setNamelist(TRIM(name), TRIM(nmlbuf))

        error = util_unlink(TRIM(filename))
    END SUBROUTINE store_and_close_namelist

    FUNCTION open_and_restore_namelist(name) RESULT(funit)
        INTEGER :: funit
        CHARACTER(len=*), INTENT(in) :: name

        TYPE(t_NamelistArchive), POINTER :: archive
        INTEGER :: flen
        CHARACTER(len=filename_max) :: filename
        CHARACTER(LEN = :), ALLOCATABLE :: nmlbuf

        archive => namelistArchive()
        CALL archive%getNamelist(name, nmlbuf)

        funit = find_next_free_unit(10,100)
        flen = util_tmpnam(filename, filename_max)
        OPEN(UNIT=funit, FILE=filename(1:flen), ACTION='write', ACCESS='stream', FORM='unformatted')
        WRITE(funit) TRIM(nmlbuf)
        CLOSE(funit)

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

END MODULE mo_restart_namelist
