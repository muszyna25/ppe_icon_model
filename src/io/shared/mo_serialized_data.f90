!> An altenative implementation of t_PackedMessage:
!! t_PackedMessage serializes the DATA via a CALL to `TRANSFER()`, which IS good for performance, however, it IS also platform dependent because of the endianess issue.
!! t_SerializedData produces a stream of ASCII coded, space separated decimal numbers, which makes the format platform independent, AND thus suitable for storing IN a file.
!!
!! Initial implementation: Nathanael Huebbe
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.

MODULE mo_serialized_data
    USE mo_exception, ONLY: finish
    USE mo_impl_constants, ONLY: SUCCESS
    USE mo_kind, ONLY: sp, dp
    USE mo_packed_message, ONLY: t_PackedMessage
    USE mo_util_string, ONLY: int2string

    IMPLICIT NONE

    PRIVATE

    PUBLIC :: t_SerializedData

    TYPE, EXTENDS(t_PackedMessage) :: t_SerializedData
    CONTAINS
        PROCEDURE :: packInt => SerializedData_packInt   ! override
        PROCEDURE :: packSingle => SerializedData_packSingle ! override
        PROCEDURE :: packDouble => SerializedData_packDouble ! override
        PROCEDURE :: packLogical => SerializedData_packLogical   ! override
        PROCEDURE :: packCharacter => SerializedData_packCharacter   ! override

        PROCEDURE :: unpackInt => SerializedData_unpackInt   ! override
        PROCEDURE :: unpackSingle => SerializedData_unpackSingle ! override
        PROCEDURE :: unpackDouble => SerializedData_unpackDouble ! override
        PROCEDURE :: unpackLogical => SerializedData_unpackLogical   ! override
        PROCEDURE :: unpackCharacter => SerializedData_unpackCharacter   ! override

        ! These provide access to the serialized DATA, allowing the caller to store it IN a file etc.
        PROCEDURE :: getData => SerializedData_getData
        PROCEDURE :: setData => SerializedData_setData

        PROCEDURE, PRIVATE :: unpackLine => SerializedData_unpackLine   ! used internally by the non-character unpackXXX() members to READ a line from the serialized DATA, which can THEN be interpreted as an INTEGER/REAL/LOGICAL
    END TYPE t_SerializedData

    CHARACTER(*), PARAMETER :: modname = "mo_serialized_data"

CONTAINS

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! packing routines !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SUBROUTINE SerializedData_packInt(me, val)
        CLASS(t_SerializedData), INTENT(INOUT) :: me
        INTEGER, VALUE :: val

        CHARACTER(*), PARAMETER :: routine = modname//":SerializedData_packInt"

        CALL me%packCharacter(TRIM(ADJUSTL(int2string(val))))
    END SUBROUTINE SerializedData_packInt

    SUBROUTINE SerializedData_packSingle(me, val)
        CLASS(t_SerializedData), INTENT(INOUT) :: me
        REAL(sp), VALUE :: val

        CHARACTER(14) :: decimalValue   ! 1 sign + 1 digit + 1 point + 7 digits + 'e' + 1 sign + 2 digits
        CHARACTER(*), PARAMETER :: routine = modname//":SerializedData_packSingle"

        WRITE(decimalValue, '(es14.7e2)') val
        CALL me%packCharacter(TRIM(ADJUSTL(decimalValue)))
    END SUBROUTINE SerializedData_packSingle

    SUBROUTINE SerializedData_packDouble(me, val)
        CLASS(t_SerializedData), INTENT(INOUT) :: me
        REAL(dp), VALUE :: val

        CHARACTER(23) :: decimalValue   ! 1 sign + 1 digit + 1 point + 15 digits + 'e' + 1 sign + 3 digits
        CHARACTER(*), PARAMETER :: routine = modname//":SerializedData_packDouble"

        WRITE(decimalValue, '(es23.15e3)') val
        CALL me%packCharacter(TRIM(ADJUSTL(decimalValue)))
    END SUBROUTINE SerializedData_packDouble

    SUBROUTINE SerializedData_packLogical(me, val)
        CLASS(t_SerializedData), INTENT(INOUT) :: me
        LOGICAL, VALUE :: val

        CHARACTER(*), PARAMETER :: routine = modname//":SerializedData_packLogical"

        IF(val) THEN
            CALL me%packCharacter("TRUE")
        ELSE
            CALL me%packCharacter("FALSE")
        END IF
    END SUBROUTINE SerializedData_packLogical

    ! This FUNCTION adds a newline CHARACTER after every string.
    ! In combination with the unpackLine() member, this allows storage of strings of flexible length,
    ! which IS used by the other packXXX() routines which provide trimmed strings.
    SUBROUTINE SerializedData_packCharacter(me, val)
        CLASS(t_SerializedData), INTENT(INOUT) :: me
        CHARACTER(*), INTENT(IN) :: val

        INTEGER :: i
        CHARACTER(*), PARAMETER :: routine = modname//":SerializedData_packCharacter"

        CALL me%ensureSpace(LEN(val) + 1)
        DO i = 1, LEN(val)
            me%messageBuffer(me%messageSize + i) = val(i:i)
        END DO
        me%messageBuffer(me%messageSize + LEN(val) + 1) = NEW_LINE(" ")  ! separate fields by a space
        me%messageSize = me%messageSize + LEN(val) + 1
    END SUBROUTINE SerializedData_packCharacter

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! unpacking routines !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SUBROUTINE SerializedData_unpackInt(me, val)
        CLASS(t_SerializedData), INTENT(INOUT) :: me
        INTEGER, INTENT(OUT) :: val

        INTEGER :: error
        CHARACTER(11) :: decimalValue   ! 1 sign + 10 digits
        CHARACTER(*), PARAMETER :: routine = modname//":SerializedData_unpackInt"

        CALL me%unpackLine(decimalValue)
        READ(decimalValue, '(i11)', iostat = error) val
        IF(error /= SUCCESS) CALL finish(routine, "error parsing the string '"//decimalValue//"' as an INTEGER")
    END SUBROUTINE SerializedData_unpackInt

    SUBROUTINE SerializedData_unpackSingle(me, val)
        CLASS(t_SerializedData), INTENT(INOUT) :: me
        REAL(sp), INTENT(OUT) :: val

        INTEGER :: error
        CHARACTER(14) :: decimalValue   ! 1 sign + 1 digit + 1 point + 7 digits + 'e' + 1 sign + 2 digits
        CHARACTER(*), PARAMETER :: routine = modname//":SerializedData_unpackSingle"

        CALL me%unpackLine(decimalValue)
        READ(decimalValue, '(f14.0)', iostat = error) val
        IF(error /= SUCCESS) CALL finish(routine, "error parsing the string '"//decimalValue//"' as a REAL")
    END SUBROUTINE SerializedData_unpackSingle

    SUBROUTINE SerializedData_unpackDouble(me, val)
        CLASS(t_SerializedData), INTENT(INOUT) :: me
        REAL(dp), INTENT(OUT) :: val

        INTEGER :: error
        CHARACTER(23) :: decimalValue   ! 1 sign + 1 digit + 1 point + 15 digits + 'e' + 1 sign + 3 digits
        CHARACTER(*), PARAMETER :: routine = modname//":SerializedData_unpackDouble"

        CALL me%unpackLine(decimalValue)
        READ(decimalValue, '(f23.0)', iostat = error) val
        IF(error /= SUCCESS) CALL finish(routine, "error parsing the string '"//decimalValue//"' as a REAL")
    END SUBROUTINE SerializedData_unpackDouble

    SUBROUTINE SerializedData_unpackLogical(me, val)
        CLASS(t_SerializedData), INTENT(INOUT) :: me
        LOGICAL, INTENT(OUT) :: val

        CHARACTER(5) :: stringValue ! LEN('FALSE')
        CHARACTER(*), PARAMETER :: routine = modname//":SerializedData_unpackLogical"

        CALL me%unpackLine(stringValue)
        IF(stringValue == "TRUE") THEN
            val = .TRUE.
        ELSE IF(stringValue == "FALSE") THEN
            val = .FALSE.
        ELSE
            CALL finish(routine, "fatal error: 'TRUE' OR 'FALSE' expected WHILE unpacking LOGICAL, &
                                 &'"//stringValue//"' found instead")
        END IF
    END SUBROUTINE SerializedData_unpackLogical

    ! The length of the CHARACTER argument must match exactly the length of the CHARACTER which was used for writing.
    ! To READ a string of flexible SIZE, unpackLine() must be used.
    SUBROUTINE SerializedData_unpackCharacter(me, val)
        CLASS(t_SerializedData), INTENT(INOUT) :: me
        CHARACTER(*), INTENT(INOUT) :: val

        INTEGER :: i
        CHARACTER(*), PARAMETER :: routine = modname//":SerializedData_unpackCharacter"

        ! buffer length check
        IF(me%readPosition + LEN(val) + 1 > me%messageSize) THEN
            CALL finish(routine, "attempt to read more data from a t_SerializedData than it contains (&
                                 &zero based read position = "//TRIM(int2string(me%readPosition))//", &
                                 &read size = "//TRIM(int2string(LEN(val) + 1))//", &
                                 &message size = "//TRIM(int2string(me%messageSize))//")")
        END IF

        ! copy the string
        DO i = 1, LEN(val)
            val(i:i) = me%messageBuffer(me%readPosition + i)
        END DO
        me%readPosition = me%readPosition + LEN(val)

        ! gobble up the trailing newline
        IF(me%messageBuffer(me%readPosition + 1) /= NEW_LINE(" ")) THEN
            CALL finish(routine, "fatal error: corrupted DATA IN a t_SerializedData detected. This could be the RESULT of reading &
                                 &a string with a different length than it was written with.")
        END IF
        me%readPosition = me%readPosition + 1
    END SUBROUTINE SerializedData_unpackCharacter

    ! Copy a string of non-newline characters to val. It IS a hard error to supply a string of insufficient SIZE.
    SUBROUTINE SerializedData_unpackLine(me, val)
        CLASS(t_SerializedData), INTENT(INOUT) :: me
        CHARACTER(*), INTENT(INOUT) :: val

        INTEGER :: i
        CHARACTER(*), PARAMETER :: routine = modname//":SerializedData_unpackLine"

        ! copy some non-newline characters to val
        DO i = 1, me%messageSize - me%readPosition
            IF(me%messageBuffer(me%readPosition + i) == NEW_LINE(" ")) EXIT
            IF(i > LEN(val)) CALL finish(routine, "fatal error: insufficient space in output string. This might be a buffer &
                                                    &size mismatch between a packXXX() and an unpackXXX() method of this class, &
                                                    &but could also be the result of corrupted data.")
            val(i:i) = me%messageBuffer(me%readPosition + i)
        END DO
        me%readPosition = me%readPosition + i - 1

        ! gobble up the trailing newline
        IF(me%messageBuffer(me%readPosition + 1) /= NEW_LINE(" ")) THEN
            CALL finish(routine, "fatal error: corrupted DATA IN a t_SerializedData detected. This could be the RESULT of reading &
                                 &a string with a different length than it was written with.")
        END IF
        me%readPosition = me%readPosition + 1
    END SUBROUTINE SerializedData_unpackLine

    FUNCTION SerializedData_getData(me) RESULT(resultVar)
        CLASS(t_SerializedData), INTENT(INOUT) :: me
        CHARACTER(:), ALLOCATABLE :: resultVar

        INTEGER :: error, i
        CHARACTER(*), PARAMETER :: routine = modname//":SerializedData_getData"

        IF(ALLOCATED(resultVar)) DEALLOCATE(resultVar)
        ALLOCATE(CHARACTER(me%messageSize) :: resultVar, STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation failure")

        DO i = 1, me%messageSize
            resultVar(i:i) = me%messageBuffer(i)
        END DO
    END FUNCTION SerializedData_getData

    SUBROUTINE SerializedData_setData(me, dataString)
        CLASS(t_SerializedData), INTENT(INOUT) :: me
        CHARACTER(*), INTENT(IN) :: dataString

        INTEGER :: error, i
        CHARACTER(*), PARAMETER :: routine = modname//":SerializedData_setData"

        me%messageSize = LEN(dataString)
        me%readPosition = 0

        ALLOCATE(me%messageBuffer(me%messageSize), STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation failure")

        DO i = 1, me%messageSize
            me%messageBuffer(i) = dataString(i:i)
        END DO
    END SUBROUTINE SerializedData_setData

END MODULE mo_serialized_data
