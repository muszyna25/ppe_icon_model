!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
MODULE mo_io_restart_attributes

  USE ISO_C_BINDING,            ONLY: C_DOUBLE, C_INT, C_INT32_T, C_INT64_T
  USE mo_cdi,                   ONLY: DATATYPE_FLT64, DATATYPE_INT32, DATATYPE_TXT, CDI_GLOBAL, vlistInqNatts, vlistInqAtt, &
                                    & vlistInqAttFlt, vlistInqAttInt, vlistInqAttTxt, vlistDefAttTxt, vlistDefAttFlt, &
                                    & vlistDefAttInt
  USE mo_exception,             ONLY: finish
  USE mo_fortran_tools,         ONLY: t_Destructible
  USE mo_hash_table,            ONLY: t_HashTable, hashTable_make, t_HashIterator
  USE mo_impl_constants,        ONLY: SUCCESS
  USE mo_kind,                  ONLY: wp, i8
  USE mo_mpi,                   ONLY: p_bcast, p_comm_rank, my_process_is_stdio
  USE mo_util_string,           ONLY: REAL2string, int2string

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: RestartAttributes_setText, RestartAttributes_getText
  PUBLIC :: RestartAttributes_setReal, RestartAttributes_getReal
  PUBLIC :: RestartAttributes_setInteger, RestartAttributes_getInteger
  PUBLIC :: RestartAttributes_setLogical, RestartAttributes_getLogical

  PUBLIC :: RestartAttributes_printAttributes, RestartAttributes_writeToFile, RestartAttributes_readFromFile

  PUBLIC :: RestartAttributes_reset

  !------------------------------------------------------------------------------------------------

  TYPE, EXTENDS(t_Destructible) :: t_Text
    CHARACTER(:), ALLOCATABLE :: text
  CONTAINS
    PROCEDURE :: destruct => Text_destruct  !< override
  END TYPE t_Text

  TYPE, EXTENDS(t_Destructible) :: t_Real
    REAL(wp) :: VALUE
  CONTAINS
    PROCEDURE :: destruct => Real_destruct  !< override
  END TYPE t_Real

  TYPE, EXTENDS(t_Destructible) :: t_Integer
    INTEGER(KIND = C_INT) :: VALUE
  CONTAINS
    PROCEDURE :: destruct => Integer_destruct   !< override
  END TYPE t_Integer

  TYPE, EXTENDS(t_Destructible) :: t_Logical
    LOGICAL :: VALUE
  CONTAINS
    PROCEDURE :: destruct => Logical_destruct   !< override
  END TYPE t_Logical

  TYPE(t_HashTable), SAVE, POINTER :: table
  LOGICAL, SAVE :: initialized = .FALSE.

  CHARACTER(LEN = *), PARAMETER :: modname = "mo_io_restart_attributes"

CONTAINS

  INTEGER(C_INT32_T) FUNCTION Text_hash(me) RESULT(RESULT)
    CLASS(t_Destructible), POINTER, INTENT(IN) :: me

    ! Just some large primes to produce some good pseudorandom bits IN the hashes.
    INTEGER(C_INT64_T), PARAMETER :: prime1 = 729326603, prime2 = 941095657, &
                                   & mask = 2_C_INT64_T**31 - 1_C_INT64_T   ! a bitmask for the 31 low order bits
    CHARACTER(LEN = *), PARAMETER :: routine = modname//":Text_hash"

    INTEGER :: i
    INTEGER(C_INT64_T) :: temp
    CLASS(*), POINTER :: meAlias    !just a workaround for a gcc compiler bug
    CHARACTER(LEN = :), POINTER :: textAlias    !just a workaround for a gcc compiler bug

    meAlias => me

    SELECT TYPE(meAlias)
        TYPE IS(t_Text)
            IF(.NOT.ALLOCATED(meAlias%text)) CALL finish(routine, "assertion failed: text object NOT initialized")
            textAlias => meAlias%text
            RESULT = 0
            DO i = 1, LEN(textAlias)
                temp = INT(IACHAR(textAlias(i:i)), C_INT64_T) * prime1 + i
                temp = IAND(temp, mask) * prime2 + INT(RESULT, C_INT64_T)
                RESULT = INT(IAND(temp, mask), C_INT32_T)
            END DO
        CLASS DEFAULT
            CALL finish(routine, "assertion failed: illegal argument TYPE")
    END SELECT
  END FUNCTION Text_hash

  LOGICAL FUNCTION Text_isEqual(me, other) RESULT(RESULT)
    CLASS(t_Destructible), POINTER, INTENT(IN) :: me, other

    CHARACTER(LEN = *), PARAMETER :: routine = modname//":Text_isEqual"
    CLASS(*), POINTER :: meAlias, otherAlias    !just a workaround for a gcc compiler bug

    meAlias => me
    otherAlias => other

    SELECT TYPE(meAlias)
        TYPE IS(t_Text)
            SELECT TYPE(otherAlias)
                TYPE IS(t_Text)
                    IF(.NOT.ALLOCATED(meAlias%text)) CALL finish(routine, "assertion failed: text object NOT initialized")
                    IF(.NOT.ALLOCATED(otherAlias%text)) CALL finish(routine, "assertion failed: text object NOT initialized")
                    RESULT = meAlias%text == otherAlias%text
                CLASS DEFAULT
                    CALL finish(routine, "assertion failed: illegal argument TYPE")
            END SELECT
        CLASS DEFAULT
            CALL finish(routine, "assertion failed: illegal argument TYPE")
    END SELECT
  END FUNCTION Text_isEqual

  SUBROUTINE Text_destruct(me)
    CLASS(t_Text), INTENT(INOUT) :: me

    DEALLOCATE(me%text)
  END SUBROUTINE Text_destruct

  SUBROUTINE Real_destruct(me)
    CLASS(t_Real), INTENT(INOUT) :: me
  END SUBROUTINE Real_destruct

  SUBROUTINE Integer_destruct(me)
    CLASS(t_Integer), INTENT(INOUT) :: me
  END SUBROUTINE Integer_destruct

  SUBROUTINE Logical_destruct(me)
    CLASS(t_Logical), INTENT(INOUT) :: me
  END SUBROUTINE Logical_destruct

  ! This takes posession of the VALUE object!
  SUBROUTINE RestartAttributes_setObject(key, VALUE)
    CHARACTER(*), INTENT(IN) :: key
    CLASS(t_Destructible), POINTER, INTENT(INOUT) :: VALUE

    INTEGER :: error
    CLASS(t_Destructible), POINTER :: keyObject
    CLASS(*), POINTER :: keyObjectAlias !XXX: workaround for gcc compiler bug
    CHARACTER(*), PARAMETER :: routine = modname//":RestartAttributes_setObject"

    ! Makes sure that there IS a hash table.
    IF(.NOT.initialized) THEN
        table => hashTable_make(Text_hash, Text_isEqual)
        initialized = .TRUE.
    END IF

    ! Create a key object AND check whether this attribute has already been set.
    ALLOCATE(t_Text :: keyObject, STAT = error)
    IF(error /= SUCCESS) CALL finish(routine, "memory allocation error")
    keyObjectAlias => keyObject
    SELECT TYPE(keyObjectAlias)
        TYPE IS(t_Text)
            keyObjectAlias%text = key
            keyObject => keyObjectAlias !Noop, but necessary to stop compilers from calling me%table%getEntry() before performing the text assignment
        CLASS DEFAULT
            CALL finish(routine, "assertion failed")
    END SELECT
    IF(ASSOCIATED(table%getEntry(keyObject))) CALL finish(routine, "double definition of restart attribute '"//key//"'")

    ! Insert the key-VALUE pair into the hash table
    CALL table%setEntry(keyObject, VALUE)
  END SUBROUTINE RestartAttributes_setObject

  SUBROUTINE RestartAttributes_reset()
    IF(initialized) THEN
        CALL table%destruct()
        DEALLOCATE(table)
        initialized = .FALSE.
    END IF
  END SUBROUTINE RestartAttributes_reset

  SUBROUTINE RestartAttributes_set(key, opt_text, opt_real, opt_integer, opt_logical)
    CHARACTER(*), INTENT(IN) :: key
    CHARACTER(*), INTENT(IN), OPTIONAL :: opt_text
    REAL(KIND = wp), INTENT(IN), OPTIONAL :: opt_real
    INTEGER(KIND = C_INT), INTENT(IN), OPTIONAL :: opt_integer
    LOGICAL, INTENT(IN), OPTIONAL :: opt_logical

    INTEGER :: error
    CLASS(t_Destructible), POINTER :: valueObject
    CLASS(*), POINTER :: valueObjectAlias   !XXX: workaround for gcc compiler bug
    CHARACTER(*), PARAMETER :: routine = modname//":RestartAttributes_set"

    IF(PRESENT(opt_text)) THEN
        ALLOCATE(t_Text :: valueObject, STAT = error)
    ELSE IF(PRESENT(opt_real)) THEN
        ALLOCATE(t_Real :: valueObject, STAT = error)
    ELSE IF(PRESENT(opt_integer)) THEN
        ALLOCATE(t_Integer :: valueObject, STAT = error)
    ELSE IF(PRESENT(opt_logical)) THEN
        ALLOCATE(t_Logical :: valueObject, STAT = error)
    ELSE
        CALL finish(routine, "assertion failed: no VALUE passed to "//routine//"()")
    END IF
    IF(error /= SUCCESS) CALL finish(routine, "memory allocation error")
    valueObjectAlias => valueObject

    SELECT TYPE(valueObjectAlias)
        TYPE IS(t_Text)
            valueObjectAlias%text = opt_text
        TYPE IS(t_Real)
            valueObjectAlias%VALUE = opt_real
        TYPE IS(t_integer)
            valueObjectAlias%VALUE = opt_integer
        TYPE IS(t_logical)
            valueObjectAlias%VALUE = opt_logical
        CLASS DEFAULT
            CALL finish(routine, "assertion failed")
    END SELECT

    SELECT TYPE(valueObjectAlias)
        CLASS IS(t_Destructible)
            valueObject => valueObjectAlias !Noop, but necessary to stop compilers from calling me%setObject() before performing the text assignment
    END SELECT

    CALL RestartAttributes_setObject(key, valueObject)
  END SUBROUTINE RestartAttributes_set

  FUNCTION RestartAttributes_getObject(key) RESULT(RESULT)
    CHARACTER(*), INTENT(IN) :: key
    CLASS(t_Destructible), POINTER :: RESULT

    INTEGER :: error
    CLASS(t_Destructible), POINTER :: keyObject
    CLASS(*), POINTER :: keyObjectAlias  !XXX: workaround for gcc compiler bug
    CHARACTER(*), PARAMETER :: routine = modname//":RestartAttributes_getObject"

    RESULT => NULL()
    IF(.NOT.initialized) RETURN

    ALLOCATE(t_Text :: keyObject, STAT = error)
    IF(error /= SUCCESS) CALL finish(routine, "memory allocation error")
    keyObjectAlias => keyObject

    SELECT TYPE(keyObjectAlias)
        TYPE IS(t_Text)
            keyObjectAlias%text = key
            keyObject => keyObjectAlias !Noop, but necessary to stop compilers from calling me%table%getEntry() before performing the text assignment
        CLASS DEFAULT
            CALL finish(routine, "assertion failed")
    END SELECT
    RESULT => table%getEntry(keyObject)
    CALL keyObject%destruct()
    DEALLOCATE(keyObject)
  END FUNCTION RestartAttributes_getObject

  LOGICAL FUNCTION RestartAttributes_get(key, opt_text, opt_real, opt_integer, opt_logical) RESULT(RESULT)
    CHARACTER(*), INTENT(IN) :: key
    CHARACTER(:), INTENT(OUT), ALLOCATABLE, OPTIONAL :: opt_text
    REAL(KIND = wp), INTENT(OUT), OPTIONAL :: opt_real
    INTEGER(KIND = C_INT), INTENT(OUT), OPTIONAL :: opt_integer
    LOGICAL, INTENT(OUT), OPTIONAL :: opt_logical

    CLASS(*), POINTER :: valueObject
    CHARACTER(*), PARAMETER :: routine = modname//"RestartAttributes_get"

    RESULT = .FALSE.
    valueObject => RestartAttributes_getObject(key)
    IF(.NOT.ASSOCIATED(valueObject)) RETURN
    RESULT = .TRUE.

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
  END FUNCTION RestartAttributes_get



  SUBROUTINE RestartAttributes_setText(key, VALUE)
    CHARACTER(*), INTENT(IN) :: key, VALUE
    CALL RestartAttributes_set(key, opt_text = VALUE)
  END SUBROUTINE RestartAttributes_setText

  SUBROUTINE RestartAttributes_setReal(key, VALUE)
    CHARACTER(*), INTENT(IN) :: key
    REAL(KIND = wp), VALUE :: VALUE
    CALL RestartAttributes_set(key, opt_real = VALUE)
  END SUBROUTINE RestartAttributes_setReal

  SUBROUTINE RestartAttributes_setInteger(key, VALUE)
    CHARACTER(*), INTENT(IN) :: key
    INTEGER(KIND = C_INT), VALUE :: VALUE
    CALL RestartAttributes_set(key, opt_integer = VALUE)
  END SUBROUTINE RestartAttributes_setInteger

  SUBROUTINE RestartAttributes_setLogical(key, VALUE)
    CHARACTER(*), INTENT(IN) :: key
    LOGICAL, VALUE :: VALUE
    CALL RestartAttributes_set(key, opt_logical = VALUE)
  END SUBROUTINE RestartAttributes_setLogical



  FUNCTION RestartAttributes_getText(key) RESULT(RESULT)
    CHARACTER(*), INTENT(IN) :: key
    CHARACTER(:), ALLOCATABLE :: RESULT
    CHARACTER(*), PARAMETER :: routine = ":RestartAttributes_getText"
    IF(.NOT.RestartAttributes_get(key, opt_text = RESULT)) CALL finish(routine, "restart attribute '"//key//"' not found")
  END FUNCTION RestartAttributes_getText

  FUNCTION RestartAttributes_getReal(key) RESULT(RESULT)
    CHARACTER(*), INTENT(IN) :: key
    REAL(KIND = wp) :: RESULT
    CHARACTER(*), PARAMETER :: routine = ":RestartAttributes_getReal"
    IF(.NOT.RestartAttributes_get(key, opt_real = RESULT)) CALL finish(routine, "restart attribute '"//key//"' not found")
  END FUNCTION RestartAttributes_getReal

  FUNCTION RestartAttributes_getInteger(key, opt_default) RESULT(RESULT)
    CHARACTER(*), INTENT(IN) :: key
    INTEGER, OPTIONAL, INTENT(IN) :: opt_default
    INTEGER(KIND = C_INT) :: RESULT
    CHARACTER(*), PARAMETER :: routine = ":RestartAttributes_getInteger"
    IF(.NOT.RestartAttributes_get(key, opt_integer = RESULT)) THEN
        IF(.NOT.PRESENT(opt_default)) CALL finish(routine, "restart attribute '"//key//"' not found")
        RESULT = opt_default
    END IF
  END FUNCTION RestartAttributes_getInteger

  FUNCTION RestartAttributes_getLogical(key) RESULT(RESULT)
    CHARACTER(*), INTENT(IN) :: key
    LOGICAL :: RESULT
    CHARACTER(*), PARAMETER :: routine = ":RestartAttributes_getLogical"
    IF(.NOT.RestartAttributes_get(key, opt_logical = RESULT)) CALL finish(routine, "restart attribute '"//key//"' not found")
  END FUNCTION RestartAttributes_getLogical



  ! Noncollective CALL, should ONLY be called by std_io process
  SUBROUTINE RestartAttributes_printAttributes()
    CHARACTER(LEN = *), PARAMETER :: routine = modname//":RestartAttributes_printAttributes"
    TYPE(t_HashIterator) :: iterator
    CLASS(t_Destructible), POINTER :: curKey, curValue
    CLASS(*), POINTER :: curKeyAlias, curValueAlias   !XXX: workaround for gcc compiler bug

    WRITE(0, *) "list of restart attributes:"
    CALL iterator%init(table)
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
  END SUBROUTINE RestartAttributes_printAttributes

  ! Noncollective CALL, should ONLY be called by the process that really does the writing.
  SUBROUTINE RestartAttributes_writeToFile(vlistId)
    INTEGER, VALUE :: vlistId

    CHARACTER(LEN = *), PARAMETER :: routine = modname//":RestartAttributes_writeToFile"
    INTEGER :: error
    TYPE(t_HashIterator) :: iterator
    CLASS(t_Destructible), POINTER :: curKey, curValue
    CLASS(*), POINTER :: curKeyAlias, curValueAlias   !XXX: workaround for gcc compiler bug
    CHARACTER(LEN = :), POINTER :: textAlias    !XXX: workaround for gcc compiler bug

    CALL iterator%init(table)
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
  END SUBROUTINE RestartAttributes_writeToFile

  ! Collective CALL: Reads AND broadcasts the restart attributes to all processes.
  SUBROUTINE RestartAttributes_readFromFile(vlistId, root_pe, comm)
    INTEGER, VALUE :: vlistId, root_pe, comm

    CHARACTER(LEN = *), PARAMETER :: routine = modname//":RestartAttributes_readFromFile"
    CHARACTER(LEN = 256) :: attributeName, text
    INTEGER :: attributeCount, attributeType, attributeLength, error, i
    REAL(KIND = C_DOUBLE) :: oneDouble(1)
    INTEGER(KIND = C_INT) :: oneInt(1)
    LOGICAL :: lread_pe     !< .TRUE., if current PE has opened the file for reading

    CALL RestartAttributes_reset()

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
                CALL RestartAttributes_setReal(TRIM(attributeName), oneDouble(1))
            CASE(DATATYPE_INT32)
                IF (lread_pe) THEN
                    error = vlistInqAttInt(vlistID, CDI_GLOBAL, TRIM(attributeName), 1, oneInt)
                    IF(error /= SUCCESS) CALL finish(routine, "error while reading restart attribute '"//TRIM(attributeName)//"'")
                END IF
                CALL p_bcast(oneInt(1), root_pe, comm)
                IF (attributeName(1:5) == 'bool_') THEN
                    CALL RestartAttributes_setLogical(TRIM(attributeName(6:)), oneInt(1) /= 0)
                ELSE
                    CALL RestartAttributes_setInteger(TRIM(attributeName), oneInt(1))
                ENDIF
            CASE(DATATYPE_TXT)
                text = ''   !This is required because vlistInqAttTxt() seems not to output a properly zero terminated string in all cases.
                IF (lread_pe) THEN
                    error = vlistInqAttTxt(vlistID, CDI_GLOBAL, TRIM(attributeName), attributeLength, text)
                    IF(error /= SUCCESS) CALL finish(routine, "error while reading restart attribute '"//TRIM(attributeName)//"'")
                END IF
                CALL p_bcast(text, root_pe, comm)
                CALL RestartAttributes_setText(TRIM(attributeName), TRIM(text))
        END SELECT
    ENDDO
  END SUBROUTINE RestartAttributes_readFromFile

END MODULE mo_io_restart_attributes
