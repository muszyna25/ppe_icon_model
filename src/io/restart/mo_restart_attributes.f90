!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
MODULE mo_restart_attributes

  USE ISO_C_BINDING,            ONLY: C_DOUBLE, C_INT, C_INT32_T, C_INT64_T
  USE mo_cdi,                   ONLY: DATATYPE_FLT64, DATATYPE_INT32, DATATYPE_TXT, CDI_GLOBAL, vlistInqNatts, vlistInqAtt, &
                                    & vlistInqAttFlt, vlistInqAttInt, vlistInqAttTxt, vlistDefAttTxt, vlistDefAttFlt, &
                                    & vlistDefAttInt
  USE mo_exception,             ONLY: finish
  USE mo_fortran_tools,         ONLY: t_Destructible
  USE mo_hash_table,            ONLY: t_HashTable, hashTable_make, t_HashIterator
  USE mo_impl_constants,        ONLY: SUCCESS
  USE mo_kind,                  ONLY: wp
  USE mo_master_config,         ONLY: isRestart
  USE mo_mpi,                   ONLY: p_bcast, p_comm_rank
  USE mo_util_string,           ONLY: real2string, int2string

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: t_RestartAttributeList, restartAttributeList_make
  PUBLIC :: setAttributesForRestarting, getAttributesForRestarting

  ! This IS basically a key/VALUE store that can handle text, REAL, INTEGER, AND LOGICAL values.
  ! The ONLY restart connected part about it IS, that it's also able to dump itself into the global attributes of a NetCDF file,
  ! AND to recreate a t_RestartAttributeList from that file, ignoring namelist attributes that may also be PRESENT IN the NetCDF file.
  TYPE, EXTENDS(t_Destructible) :: t_RestartAttributeList
    TYPE(t_HashTable), POINTER :: table

  CONTAINS
    PROCEDURE :: setText => RestartAttributeList_setText
    PROCEDURE :: setReal => RestartAttributeList_setReal
    PROCEDURE :: setInteger => RestartAttributeList_setInteger
    PROCEDURE :: setLogical => RestartAttributeList_setLogical
    PROCEDURE :: getText => RestartAttributeList_getText
    PROCEDURE :: getReal => RestartAttributeList_getReal
    PROCEDURE :: getInteger => RestartAttributeList_getInteger
    PROCEDURE :: getLogical => RestartAttributeList_getLogical

    PROCEDURE :: printAttributes => RestartAttributeList_printAttributes
    PROCEDURE :: writeToCdiVlist => RestartAttributeList_writeToCdiVlist

    PROCEDURE :: destruct => RestartAttributeList_destruct  ! override

    !-------------------------------------------------------------------------------------------------------------------------------

    PROCEDURE, PRIVATE :: setObject => RestartAttributeList_setObject
    PROCEDURE, PRIVATE :: set => RestartAttributeList_set
    PROCEDURE, PRIVATE :: getObject => RestartAttributeList_getObject
    PROCEDURE, PRIVATE :: get => RestartAttributeList_get

    PROCEDURE, PRIVATE :: readFromFile => RestartAttributeList_readFromFile
  END TYPE t_RestartAttributeList

  INTERFACE restartAttributeList_make
    MODULE PROCEDURE restartAttributeList_makeEmpty
    MODULE PROCEDURE restartAttributeList_makeFromFile
  END INTERFACE restartAttributeList_make

  !------------------------------------------------------------------------------------------------

  TYPE, EXTENDS(t_Destructible) :: t_Text
    CHARACTER(:), ALLOCATABLE :: text
  CONTAINS
    PROCEDURE :: destruct => text_destruct  !< override
  END TYPE t_Text

  TYPE, EXTENDS(t_Destructible) :: t_Real
    REAL(wp) :: VALUE
  CONTAINS
    PROCEDURE :: destruct => real_destruct  !< override
  END TYPE t_Real

  TYPE, EXTENDS(t_Destructible) :: t_Integer
    INTEGER(KIND = C_INT) :: VALUE
  CONTAINS
    PROCEDURE :: destruct => integer_destruct   !< override
  END TYPE t_Integer

  TYPE, EXTENDS(t_Destructible) :: t_Logical
    LOGICAL :: VALUE
  CONTAINS
    PROCEDURE :: destruct => logical_destruct   !< override
  END TYPE t_Logical

  TYPE(t_RestartAttributeList), SAVE, POINTER :: gRestartAttributes
  LOGICAL, SAVE :: gRestartAttributes_initialized = .FALSE.

  CHARACTER(LEN = *), PARAMETER :: modname = "mo_restart_attributes"

CONTAINS

  FUNCTION text_create(VALUE) RESULT(resultVar)
    CHARACTER(*), INTENT(IN) :: VALUE
    TYPE(t_Text), POINTER :: resultVar

    integer :: error
    CHARACTER(*), PARAMETER :: routine = modname//":text_create"

    ALLOCATE(resultVar, STAT = error)
    IF(error /= SUCCESS) CALL finish(routine, "memory allocation failed")
    resultVar%text = VALUE
  END FUNCTION text_create

  FUNCTION real_create(VALUE) RESULT(resultVar)
    REAL(wp), VALUE :: VALUE
    TYPE(t_Real), POINTER :: resultVar

    integer :: error
    CHARACTER(*), PARAMETER :: routine = modname//":real_create"

    ALLOCATE(resultVar, STAT = error)
    IF(error /= SUCCESS) CALL finish(routine, "memory allocation failed")
    resultVar%VALUE = VALUE
  END FUNCTION real_create

  FUNCTION integer_create(VALUE) RESULT(resultVar)
    INTEGER(KIND = C_INT), VALUE :: VALUE
    TYPE(t_Integer), POINTER :: resultVar

    integer :: error
    CHARACTER(*), PARAMETER :: routine = modname//":integer_create"

    ALLOCATE(resultVar, STAT = error)
    IF(error /= SUCCESS) CALL finish(routine, "memory allocation failed")
    resultVar%VALUE = VALUE
  END FUNCTION integer_create

  FUNCTION logical_create(VALUE) RESULT(resultVar)
    LOGICAL, VALUE :: VALUE
    TYPE(t_Logical), POINTER :: resultVar

    integer :: error
    CHARACTER(*), PARAMETER :: routine = modname//":logical_create"

    ALLOCATE(resultVar, STAT = error)
    IF(error /= SUCCESS) CALL finish(routine, "memory allocation failed")
    resultVar%VALUE = VALUE
  END FUNCTION logical_create

  ! Sets the restart attribute list that IS to be used for restarting. Must NOT be set when we are NOT restarting.
  SUBROUTINE setAttributesForRestarting(restartAttributes)
    TYPE(t_RestartAttributeList), POINTER, INTENT(INOUT) :: restartAttributes
    CHARACTER(LEN = *), PARAMETER :: routine = modname//":setAttributesForRestarting"

    IF(.NOT.isRestart()) CALL finish(routine, "setAttributesForRestarting() must only be called in a restarted run")
    IF(gRestartAttributes_initialized) CALL finish(routine, "setAttributesForRestarting() called several times")
    IF(.NOT.ASSOCIATED(restartAttributes)) THEN
        CALL finish(routine, "argument to setAttributesForRestarting() must be an associated pointer")
    END IF
    gRestartAttributes => restartAttributes
    gRestartAttributes_initialized = .TRUE.
  END SUBROUTINE setAttributesForRestarting

  ! Returns the restart attribute list that we are currently using to initialize the model. Returns NULL on non-restart runs.
  FUNCTION getAttributesForRestarting() RESULT(resultVar)
    TYPE(t_RestartAttributeList), POINTER :: resultVar

    CHARACTER(LEN = *), PARAMETER :: routine = modname//":getAttributesForRestarting"

    IF(isRestart()) THEN
        IF(.NOT.gRestartAttributes_initialized) THEN
            CALL finish(routine, "assertion failed: getAttributesForRestarting() called before restart attributes were read from &
                                 &file")
        END IF
        resultVar => gRestartAttributes
    ELSE
        resultVar => NULL()
    END IF
  END FUNCTION getAttributesForRestarting

  INTEGER(C_INT32_T) FUNCTION text_hash(me) RESULT(resultVar)
    CLASS(t_Destructible), POINTER, INTENT(IN) :: me

    ! Just some large primes to produce some good pseudorandom bits IN the hashes.
    INTEGER(C_INT64_T), PARAMETER :: prime1 = 729326603, prime2 = 941095657, &
                                   & mask = 2_C_INT64_T**31 - 1_C_INT64_T   ! a bitmask for the 31 low order bits
    CHARACTER(LEN = *), PARAMETER :: routine = modname//":text_hash"

    INTEGER :: i
    INTEGER(C_INT64_T) :: temp
    CLASS(*), POINTER :: meAlias    !just a workaround for a gcc compiler bug
    CHARACTER(LEN = :), POINTER :: textAlias    !just a workaround for a gcc compiler bug

    meAlias => me

    SELECT TYPE(meAlias)
        TYPE IS(t_Text)
            IF(.NOT.ALLOCATED(meAlias%text)) CALL finish(routine, "assertion failed: text object NOT initialized")
            textAlias => meAlias%text
            resultVar = 0
            DO i = 1, LEN(textAlias)
                temp = INT(IACHAR(textAlias(i:i)), C_INT64_T) * prime1 + i
                temp = IAND(temp, mask) * prime2 + INT(resultVar, C_INT64_T)
                resultVar = INT(IAND(temp, mask), C_INT32_T)
            END DO
        CLASS DEFAULT
            CALL finish(routine, "assertion failed: illegal argument TYPE")
    END SELECT
  END FUNCTION text_hash

  LOGICAL FUNCTION text_isEqual(me, other) RESULT(resultVar)
    CLASS(t_Destructible), POINTER, INTENT(IN) :: me, other

    CHARACTER(LEN = *), PARAMETER :: routine = modname//":text_isEqual"
    CLASS(*), POINTER :: meAlias, otherAlias    !just a workaround for a gcc compiler bug

    meAlias => me
    otherAlias => other

    SELECT TYPE(meAlias)
        TYPE IS(t_Text)
            SELECT TYPE(otherAlias)
                TYPE IS(t_Text)
                    IF(.NOT.ALLOCATED(meAlias%text)) CALL finish(routine, "assertion failed: text object NOT initialized")
                    IF(.NOT.ALLOCATED(otherAlias%text)) CALL finish(routine, "assertion failed: text object NOT initialized")
                    resultVar = meAlias%text == otherAlias%text
                CLASS DEFAULT
                    CALL finish(routine, "assertion failed: illegal argument TYPE")
            END SELECT
        CLASS DEFAULT
            CALL finish(routine, "assertion failed: illegal argument TYPE")
    END SELECT
  END FUNCTION text_isEqual

  SUBROUTINE text_destruct(me)
    CLASS(t_Text), INTENT(INOUT) :: me

    DEALLOCATE(me%text)
  END SUBROUTINE text_destruct

  SUBROUTINE real_destruct(me)
    CLASS(t_Real), INTENT(INOUT) :: me
  END SUBROUTINE real_destruct

  SUBROUTINE integer_destruct(me)
    CLASS(t_Integer), INTENT(INOUT) :: me
  END SUBROUTINE integer_destruct

  SUBROUTINE logical_destruct(me)
    CLASS(t_Logical), INTENT(INOUT) :: me
  END SUBROUTINE logical_destruct

  FUNCTION RestartAttributeList_makeEmpty() RESULT(resultVar)
    TYPE(t_RestartAttributeList), POINTER :: resultVar

    INTEGER :: error
    CHARACTER(LEN = *), PARAMETER :: routine = modname//':RestartAttributeList_makeEmpty'

    ALLOCATE(resultVar, STAT = error)
    IF(error /= SUCCESS) CALL finish(routine, "memory allocation failure")
    resultVar%table => hashTable_make(text_hash, text_isEqual)
  END FUNCTION RestartAttributeList_makeEmpty

  FUNCTION RestartAttributeList_makeFromFile(vlistId, root_pe, comm) RESULT(resultVar)
    INTEGER, VALUE :: vlistId, root_pe, comm
    TYPE(t_RestartAttributeList), POINTER :: resultVar

    resultVar => RestartAttributeList_makeEmpty()
    CALL resultVar%readFromFile(vlistId, root_pe, comm)
  END FUNCTION RestartAttributeList_makeFromFile

  ! This takes posession of the VALUE object!
  SUBROUTINE RestartAttributeList_setObject(me, key, VALUE)
    CLASS(t_RestartAttributeList), INTENT(INOUT) :: me
    CHARACTER(*), INTENT(IN) :: key
    CLASS(t_Destructible), POINTER, INTENT(INOUT) :: VALUE

    CLASS(t_Destructible), POINTER :: keyObject
    CHARACTER(*), PARAMETER :: routine = modname//":RestartAttributeList_setObject"

    ! Create a key object AND check whether this attribute has already been set.
    keyObject => text_create(key)
    IF(ASSOCIATED(me%table%getEntry(keyObject))) CALL finish(routine, "double definition of restart attribute '"//key//"'")

    ! Insert the key-VALUE pair into the hash table
    CALL me%table%setEntry(keyObject, VALUE)
  END SUBROUTINE RestartAttributeList_setObject

  SUBROUTINE RestartAttributeList_set(me, key, opt_text, opt_real, opt_integer, opt_logical)
    CLASS(t_RestartAttributeList), INTENT(INOUT) :: me
    CHARACTER(*), INTENT(IN) :: key
    CHARACTER(*), INTENT(IN), OPTIONAL :: opt_text
    REAL(KIND = wp), INTENT(IN), OPTIONAL :: opt_real
    INTEGER(KIND = C_INT), INTENT(IN), OPTIONAL :: opt_integer
    LOGICAL, INTENT(IN), OPTIONAL :: opt_logical

    CLASS(t_Destructible), POINTER :: valueObject
    CHARACTER(*), PARAMETER :: routine = modname//":RestartAttributeList_set"

    IF(PRESENT(opt_text)) THEN
        valueObject => text_create(opt_text)
    ELSE IF(PRESENT(opt_real)) THEN
        valueObject => real_create(opt_real)
    ELSE IF(PRESENT(opt_integer)) THEN
        valueObject => integer_create(opt_integer)
    ELSE IF(PRESENT(opt_logical)) THEN
        valueObject => logical_create(opt_logical)
    ELSE
        CALL finish(routine, "assertion failed: no VALUE passed to "//routine//"()")
    END IF
    CALL me%setObject(key, valueObject)
  END SUBROUTINE RestartAttributeList_set

  FUNCTION RestartAttributeList_getObject(me, key) RESULT(resultVar)
    CLASS(t_RestartAttributeList), INTENT(IN) :: me
    CHARACTER(*), INTENT(IN) :: key
    CLASS(t_Destructible), POINTER :: resultVar

    CLASS(t_Destructible), POINTER :: keyObject
    CHARACTER(*), PARAMETER :: routine = modname//":RestartAttributeList_getObject"

    resultVar => NULL()
    keyObject => text_create(key)
    resultVar => me%table%getEntry(keyObject)
    CALL keyObject%destruct()
    DEALLOCATE(keyObject)
  END FUNCTION RestartAttributeList_getObject

  LOGICAL FUNCTION RestartAttributeList_get(me, key, opt_text, opt_real, opt_integer, opt_logical) RESULT(resultVar)
    CLASS(t_RestartAttributeList), INTENT(IN) :: me
    CHARACTER(*), INTENT(IN) :: key
    CHARACTER(:), INTENT(OUT), ALLOCATABLE, OPTIONAL :: opt_text
    REAL(KIND = wp), INTENT(OUT), OPTIONAL :: opt_real
    INTEGER(KIND = C_INT), INTENT(OUT), OPTIONAL :: opt_integer
    LOGICAL, INTENT(OUT), OPTIONAL :: opt_logical

    CLASS(*), POINTER :: valueObject
    CHARACTER(*), PARAMETER :: routine = modname//"RestartAttributeList_get"

    resultVar = .FALSE.
    valueObject => me%getObject(key)
    IF(.NOT.ASSOCIATED(valueObject)) RETURN
    resultVar = .TRUE.

    SELECT TYPE(valueObject)
        TYPE IS(t_Text)
            IF(.NOT.PRESENT(opt_text)) CALL finish(routine, "type mismatch while reading restart attribute '"//key// &
                                                            "', which is of type TEXT")
            opt_text = valueObject%text
        TYPE IS(t_Real)
            IF(.NOT.PRESENT(opt_real)) CALL finish(routine, "type mismatch while reading restart attribute '"//key// &
                                                            "', which is of type REAL")
            opt_real = valueObject%VALUE
        TYPE IS(t_integer)
            IF(.NOT.PRESENT(opt_integer)) CALL finish(routine, "type mismatch while reading restart attribute '"//key// &
                                                               "', which is of type INTEGER")
            opt_integer = valueObject%VALUE
        TYPE IS(t_logical)
            IF(.NOT.PRESENT(opt_logical)) CALL finish(routine, "type mismatch while reading restart attribute '"//key// &
                                                               "', which is of type LOGICAL")
            opt_logical = valueObject%VALUE
        CLASS DEFAULT
            CALL finish(routine, "assertion failed")
    END SELECT
  END FUNCTION RestartAttributeList_get



  SUBROUTINE RestartAttributeList_setText(me, key, VALUE)
    CLASS(t_RestartAttributeList), INTENT(INOUT) :: me
    CHARACTER(*), INTENT(IN) :: key, VALUE
    CALL me%set(key, opt_text = VALUE)
  END SUBROUTINE RestartAttributeList_setText

  SUBROUTINE RestartAttributeList_setReal(me, key, VALUE)
    CLASS(t_RestartAttributeList), INTENT(INOUT) :: me
    CHARACTER(*), INTENT(IN) :: key
    REAL(KIND = wp), VALUE :: VALUE
    CALL me%set(key, opt_real = VALUE)
  END SUBROUTINE RestartAttributeList_setReal

  SUBROUTINE RestartAttributeList_setInteger(me, key, VALUE)
    CLASS(t_RestartAttributeList), INTENT(INOUT) :: me
    CHARACTER(*), INTENT(IN) :: key
    INTEGER(KIND = C_INT), VALUE :: VALUE
    CALL me%set(key, opt_integer = VALUE)
  END SUBROUTINE RestartAttributeList_setInteger

  SUBROUTINE RestartAttributeList_setLogical(me, key, VALUE)
    CLASS(t_RestartAttributeList), INTENT(INOUT) :: me
    CHARACTER(*), INTENT(IN) :: key
    LOGICAL, VALUE :: VALUE
    CALL me%set(key, opt_logical = VALUE)
  END SUBROUTINE RestartAttributeList_setLogical



  FUNCTION RestartAttributeList_getText(me, key) RESULT(resultVar)
    CLASS(t_RestartAttributeList), INTENT(IN) :: me
    CHARACTER(*), INTENT(IN) :: key
    CHARACTER(:), ALLOCATABLE :: resultVar
    CHARACTER(*), PARAMETER :: routine = ":RestartAttributeList_getText"
    IF(.NOT.me%get(key, opt_text = resultVar)) CALL finish(routine, "restart attribute '"//key//"' not found")
  END FUNCTION RestartAttributeList_getText

  FUNCTION RestartAttributeList_getReal(me, key) RESULT(resultVar)
    CLASS(t_RestartAttributeList), INTENT(IN) :: me
    CHARACTER(*), INTENT(IN) :: key
    REAL(KIND = wp) :: resultVar
    CHARACTER(*), PARAMETER :: routine = ":RestartAttributeList_getReal"
    IF(.NOT.me%get(key, opt_real = resultVar)) CALL finish(routine, "restart attribute '"//key//"' not found")
  END FUNCTION RestartAttributeList_getReal

  FUNCTION RestartAttributeList_getInteger(me, key, opt_default) RESULT(resultVar)
    CLASS(t_RestartAttributeList), INTENT(IN) :: me
    CHARACTER(*), INTENT(IN) :: key
    INTEGER, OPTIONAL, INTENT(IN) :: opt_default
    INTEGER(KIND = C_INT) :: resultVar
    CHARACTER(*), PARAMETER :: routine = ":RestartAttributeList_getInteger"
    IF(.NOT.me%get(key, opt_integer = resultVar)) THEN
        IF(.NOT.PRESENT(opt_default)) CALL finish(routine, "restart attribute '"//key//"' not found")
        resultVar = opt_default
    END IF
  END FUNCTION RestartAttributeList_getInteger

  FUNCTION RestartAttributeList_getLogical(me, key) RESULT(resultVar)
    CLASS(t_RestartAttributeList), INTENT(IN) :: me
    CHARACTER(*), INTENT(IN) :: key
    LOGICAL :: resultVar
    CHARACTER(*), PARAMETER :: routine = ":RestartAttributeList_getLogical"
    IF(.NOT.me%get(key, opt_logical = resultVar)) CALL finish(routine, "restart attribute '"//key//"' not found")
  END FUNCTION RestartAttributeList_getLogical



  ! Noncollective CALL, should ONLY be called by std_io process
  SUBROUTINE RestartAttributeList_printAttributes(me)
    CLASS(t_RestartAttributeList), INTENT(INOUT) :: me
    CHARACTER(LEN = *), PARAMETER :: routine = modname//":RestartAttributeList_printAttributes"
    TYPE(t_HashIterator) :: iterator
    CLASS(t_Destructible), POINTER :: curKey, curValue
    CLASS(*), POINTER :: curKeyAlias, curValueAlias   !XXX: workaround for gcc compiler bug

    WRITE(0, *) "list of restart attributes:"
    CALL iterator%init(me%table)
    DO
        IF(.NOT.iterator%nextEntry(curKey, curValue)) EXIT

        curKeyAlias => curKey
        curValueAlias => curValue

        SELECT TYPE(curKeyAlias)
            TYPE IS(t_Text)
                SELECT TYPE(curValueAlias)
                    TYPE IS(t_Text)
                        WRITE(0, *) "    '"//curKeyAlias%text//"' = '"//curValueAlias%text//"'"
                    TYPE IS(t_Real)
                        WRITE(0, *) "    '"//curKeyAlias%text//"' = "//TRIM(real2string(curValueAlias%VALUE))
                    TYPE IS(t_integer)
                        WRITE(0, *) "    '"//curKeyAlias%text//"' = INT("//TRIM(int2string(curValueAlias%VALUE))//")"
                    TYPE IS(t_logical)
                        IF(curValueAlias%VALUE) THEN
                            WRITE(0, *) "    '"//curKeyAlias%text//"' = .TRUE."
                        ELSE
                            WRITE(0, *) "    '"//curKeyAlias%text//"' = .FALSE."
                        END IF
                    CLASS DEFAULT
                        CALL finish(routine, "assertion failed")
                END SELECT
            CLASS DEFAULT
                CALL finish(routine, "assertion failed")
        END SELECT
    END DO
  END SUBROUTINE RestartAttributeList_printAttributes

  ! Noncollective CALL, should ONLY be called by the process that really does the writing.
  SUBROUTINE RestartAttributeList_writeToCdiVlist(me, vlistId)
    CLASS(t_RestartAttributeList), INTENT(INOUT) :: me
    INTEGER, VALUE :: vlistId

    CHARACTER(LEN = *), PARAMETER :: routine = modname//":RestartAttributeList_writeToCdiVlist"
    INTEGER :: error
    TYPE(t_HashIterator) :: iterator
    CLASS(t_Destructible), POINTER :: curKey, curValue
    CLASS(*), POINTER :: curKeyAlias, curValueAlias   !XXX: workaround for gcc compiler bug
    CHARACTER(LEN = :), POINTER :: textAlias    !XXX: workaround for gcc compiler bug

    CALL iterator%init(me%table)
    DO
        IF(.NOT.iterator%nextEntry(curKey, curValue)) EXIT

        curKeyAlias => curKey
        curValueAlias => curValue

        SELECT TYPE(curKeyAlias)
            TYPE IS(t_Text)
                SELECT TYPE(curValueAlias)
                    TYPE IS(t_Text)
                        textAlias => curValueAlias%text
                        error = vlistDefAttTxt(vlistId, CDI_GLOBAL, curKeyAlias%text, LEN(textAlias), curValueAlias%text)
                    TYPE IS(t_Real)
                        error = vlistDefAttFlt(vlistId, CDI_GLOBAL, curKeyAlias%text, DATATYPE_FLT64, 1, [curValueAlias%VALUE])
                    TYPE IS(t_integer)
                        error = vlistDefAttInt(vlistId, CDI_GLOBAL, curKeyAlias%text, DATATYPE_INT32, 1, [curValueAlias%VALUE])
                    TYPE IS(t_logical)
                        IF(curValueAlias%VALUE) THEN
                            error = vlistDefAttInt(vlistId, CDI_GLOBAL, 'bool_'//curKeyAlias%text, DATATYPE_INT32, 1, [1])
                        ELSE
                            error = vlistDefAttInt(vlistId, CDI_GLOBAL, 'bool_'//curKeyAlias%text, DATATYPE_INT32, 1, [0])
                        END IF
                    CLASS DEFAULT
                        CALL finish(routine, "assertion failed")
                END SELECT
                IF(error /= SUCCESS) CALL finish(routine, "error writing restart attribute to file")
            CLASS DEFAULT
                CALL finish(routine, "assertion failed")
        END SELECT
    END DO
  END SUBROUTINE RestartAttributeList_writeToCdiVlist

  ! Collective CALL: Reads AND broadcasts the restart attributes to all processes.
  ! TODO[NH]: Fuse into makeFromFile() to guarantee that we start with an empty AttributeList.
  SUBROUTINE RestartAttributeList_readFromFile(me, vlistId, root_pe, comm)
    CLASS(t_RestartAttributeList), INTENT(INOUT) :: me
    INTEGER, VALUE :: vlistId, root_pe, comm

    CHARACTER(LEN = *), PARAMETER :: routine = modname//":RestartAttributeList_readFromFile"
    CHARACTER(LEN = 256) :: attributeName, text
    INTEGER :: attributeCount, attributeType, attributeLength, error, i
    REAL(KIND = C_DOUBLE) :: oneDouble(1)
    INTEGER(KIND = C_INT) :: oneInt(1)
    LOGICAL :: lread_pe     !< .TRUE., if current PE has opened the file for reading

    lread_pe = p_comm_rank(comm) == root_pe
    IF (lread_pe) THEN
        error = vlistInqNatts(vlistID, CDI_GLOBAL, attributeCount)
        IF(error /= SUCCESS) CALL finish(routine, "couldn't inquire number of attributes from the restart file")
    END IF
    CALL p_bcast(attributeCount, root_pe, comm)

    DO i = 0, attributeCount - 1
        IF (lread_pe) THEN
            error = vlistInqAtt(vlistId, CDI_GLOBAL, i, attributeName, attributeType, attributeLength)
            IF(error /= SUCCESS) CALL finish(routine, "error while reading attributes from the restart file")
        END IF
        CALL p_bcast(attributeName, root_pe, comm)
        IF ( attributeName(1:4) == 'nml_') CYCLE ! skip this, it is a namelist

        CALL p_bcast(attributeType, root_pe, comm)
        SELECT CASE(attributeType)
            CASE(DATATYPE_FLT64)
                IF (lread_pe) THEN
                    error = vlistInqAttFlt(vlistID, CDI_GLOBAL, TRIM(attributeName), 1, oneDouble)
                    IF(error /= SUCCESS) CALL finish(routine, "error while reading restart attribute '"//TRIM(attributeName)//"'")
                END IF
                CALL p_bcast(oneDouble(1), root_pe, comm)
                CALL me%setReal(TRIM(attributeName), oneDouble(1))
            CASE(DATATYPE_INT32)
                IF (lread_pe) THEN
                    error = vlistInqAttInt(vlistID, CDI_GLOBAL, TRIM(attributeName), 1, oneInt)
                    IF(error /= SUCCESS) CALL finish(routine, "error while reading restart attribute '"//TRIM(attributeName)//"'")
                END IF
                CALL p_bcast(oneInt(1), root_pe, comm)
                IF (attributeName(1:5) == 'bool_') THEN
                    CALL me%setLogical(TRIM(attributeName(6:)), oneInt(1) /= 0)
                ELSE
                    CALL me%setInteger(TRIM(attributeName), oneInt(1))
                ENDIF
            CASE(DATATYPE_TXT)
                text = ''   !This is required because vlistInqAttTxt() seems not to output a properly zero terminated string in all cases.
                IF (lread_pe) THEN
                    error = vlistInqAttTxt(vlistID, CDI_GLOBAL, TRIM(attributeName), attributeLength, text)
                    IF(error /= SUCCESS) CALL finish(routine, "error while reading restart attribute '"//TRIM(attributeName)//"'")
                END IF
                CALL p_bcast(text, root_pe, comm)
                CALL me%setText(TRIM(attributeName), TRIM(text))
        END SELECT
    ENDDO
  END SUBROUTINE RestartAttributeList_readFromFile

  SUBROUTINE RestartAttributeList_destruct(me)
    CLASS(t_RestartAttributeList), INTENT(INOUT) :: me

    CALL me%table%destruct()
    DEALLOCATE(me%table)
  END SUBROUTINE RestartAttributeList_destruct

END MODULE mo_restart_attributes
