!> A CLASS to make sending/receiving of DATA packets easier.
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

MODULE mo_packed_message
    USE mo_exception, ONLY: finish
    USE mo_impl_constants, ONLY: SUCCESS
    USE mo_kind, ONLY: sp, dp
    USE mo_mpi, ONLY: p_int, p_comm_rank
    USE mo_util_string, ONLY: int2string

#ifndef NOMPI
    USE mo_mpi, ONLY: MPI_SUCCESS
    USE mpi
#endif

    IMPLICIT NONE
    PRIVATE

    PUBLIC :: t_PackedMessage
    INTEGER, PARAMETER, PUBLIC :: kPackOp = 0, kUnpackOp = 1

    TYPE :: t_PackedMessage
        INTEGER messageSize, readPosition
        CHARACTER, POINTER :: messageBuffer(:)
    CONTAINS
        PROCEDURE :: construct => PackedMessage_construct
        PROCEDURE :: reset => PackedMessage_reset   ! functionally equivalent to `CALL message%destruct(); CALL message%construct()`, but more efficient (buffer IS reused)

        ! pack routines
        PROCEDURE :: packInt => PackedMessage_packInt
        PROCEDURE :: packSingle => PackedMessage_packSingle
        PROCEDURE :: packDouble => PackedMessage_packDouble
        PROCEDURE :: packLogical => PackedMessage_packLogical

        PROCEDURE :: packIntArray => PackedMessage_packIntArray
        PROCEDURE :: packSingleArray => PackedMessage_packSingleArray
        PROCEDURE :: packDoubleArray => PackedMessage_packDoubleArray
        PROCEDURE :: packLogicalArray => PackedMessage_packLogicalArray

        GENERIC :: pack => packInt, packSingle, packDouble, packLogical, packIntArray, packSingleArray, packDoubleArray, &
                         & packLogicalArray

        ! unpack routines
        PROCEDURE :: unpackInt => PackedMessage_unpackInt
        PROCEDURE :: unpackSingle => PackedMessage_unpackSingle
        PROCEDURE :: unpackDouble => PackedMessage_unpackDouble
        PROCEDURE :: unpackLogical => PackedMessage_unpackLogical

        PROCEDURE :: unpackIntArray => PackedMessage_unpackIntArray
        PROCEDURE :: unpackSingleArray => PackedMessage_unpackSingleArray
        PROCEDURE :: unpackDoubleArray => PackedMessage_unpackDoubleArray
        PROCEDURE :: unpackLogicalArray => PackedMessage_unpackLogicalArray

        GENERIC :: unpack => unpackInt, unpackSingle, unpackDouble, unpackLogical, unpackIntArray, unpackSingleArray, &
                           & unpackDoubleArray, unpackLogicalArray

        ! routines to facilitate packing AND unpacking with the same code
        ! USE of these routine prohibits ANY errors by mismatches between packing AND unpacking code
        PROCEDURE :: executeInt => PackedMessage_executeInt
        PROCEDURE :: executeSingle => PackedMessage_executeSingle
        PROCEDURE :: executeDouble => PackedMessage_executeDouble
        PROCEDURE :: executeLogical => PackedMessage_executeLogical

        PROCEDURE :: executeIntArray => PackedMessage_executeIntArray
        PROCEDURE :: executeSingleArray => PackedMessage_executeSingleArray
        PROCEDURE :: executeDoubleArray => PackedMessage_executeDoubleArray
        PROCEDURE :: executeLogicalArray => PackedMessage_executeLogicalArray

        GENERIC :: execute => executeInt, executeSingle, executeDouble, executeLogical, executeIntArray, executeSingleArray, &
                            & executeDoubleArray, executeLogicalArray

        ! communication routines
        ! All of these will flush ANY contents of the reciever(s), replacing it with a copy of the sender's packet.
        PROCEDURE :: send => PackedMessage_send
        PROCEDURE :: recv => PackedMessage_recv
        PROCEDURE :: bcast => PackedMessage_bcast

        ! destructor
        PROCEDURE :: destruct => PackedMessage_destruct

        PROCEDURE, PRIVATE :: ensureSpace => PackedMessage_ensureSpace
        PROCEDURE, PRIVATE :: printState => PackedMessage_printState

        ! these two methods are ONLY ever invoked via the two preprocessor macros below
        PROCEDURE, PRIVATE :: packBytes => PackedMessage_packBytes
        PROCEDURE, PRIVATE :: unpackBytes => PackedMessage_unpackBytes  ! returns a POINTER to a subarray of the internal storage, so that POINTER IS invalidated by ANY (IN)direct CALL to ensureSpace()!
#define doPacking(me, VALUE) CALL me%packBytes(TRANSFER(VALUE, characterArrayMold))
#define doUnpacking(me, VALUE) VALUE = TRANSFER(me%unpackBytes(TRANSFER(VALUE, characterArrayMold)), VALUE)
    END TYPE t_PackedMessage

    CHARACTER :: characterMold, characterArrayMold(1)
    INTEGER :: integerMold
    REAL(sp) :: singleMold
    REAL(dp) :: doubleMold
    LOGICAL :: logicalMold

    CHARACTER(*), PARAMETER :: modname = "mo_packed_message"

CONTAINS

    SUBROUTINE PackedMessage_construct(me)
        CLASS(t_PackedMessage), INTENT(OUT) :: me

        INTEGER :: error
        CHARACTER(LEN = *), PARAMETER :: routine = modname//":PackedMessage_construct"

        me%messageSize = 0
        me%readPosition = 0
        ALLOCATE(me%messageBuffer(256), STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation failed")
    END SUBROUTINE PackedMessage_construct

    SUBROUTINE PackedMessage_reset(me)
        CLASS(t_PackedMessage), INTENT(INOUT) :: me

        me%messageSize = 0
        me%readPosition = 0
    END SUBROUTINE PackedMessage_reset

    SUBROUTINE PackedMessage_ensureSpace(me, requiredSpace)
        CLASS(t_PackedMessage), INTENT(INOUT) :: me
        INTEGER, VALUE :: requiredSpace

        INTEGER :: error, oldSize, newSize
        CHARACTER, POINTER :: newBuffer(:)
        CHARACTER(LEN = *), PARAMETER :: routine = modname//":PackedMessage_construct"

        !fast path IF the current buffer IS sufficient
        oldSize = SIZE(me%messageBuffer)
        IF(oldSize >= me%messageSize + requiredSpace) RETURN

        !ALLOCATE a new buffer that's at least twice as large as the old one, AND at least sufficient
        newSize = MAX(2*oldSize, me%messageSize + requiredSpace)
        ALLOCATE(newBuffer(newSize), STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation failed")

        !copy DATA over
        newBuffer(1:oldSize) = me%messageBuffer(1:oldSize)

        !get rid of the old buffer
        DEALLOCATE(me%messageBuffer)
        me%messageBuffer => newBuffer
    END SUBROUTINE PackedMessage_ensureSpace

    SUBROUTINE PackedMessage_packBytes(me, bytes)
        CLASS(t_PackedMessage), INTENT(INOUT) :: me
        CHARACTER, INTENT(IN) :: bytes(:)

        INTEGER :: i

        CALL me%ensureSpace(SIZE(bytes))
        DO i = 1, SIZE(bytes)
            me%messageBuffer(me%messageSize + i) = bytes(i)
        END DO
        me%messageSize = me%messageSize + SIZE(bytes)
    END SUBROUTINE PackedMessage_packBytes

    FUNCTION PackedMessage_unpackBytes(me, byteArrayMold) RESULT(RESULT)
        CLASS(t_PackedMessage), INTENT(INOUT) :: me
        CHARACTER, INTENT(IN) :: byteArrayMold(:)
        CHARACTER, POINTER :: RESULT(:)

        CHARACTER(LEN = *), PARAMETER :: routine = modname//":PackedMessage_unpackBytes"

        IF(me%readPosition + SIZE(byteArrayMold) > SIZE(me%messageBuffer)) THEN
            CALL finish(routine, "attempt to read more data from a t_PackedMessage than it contains (&
                                 &zero based read position = "//TRIM(int2string(me%readPosition))//", &
                                 &read size = "//TRIM(int2string(SIZE(byteArrayMold)))//", &
                                 &message size = "//TRIM(int2string(SIZE(me%messageBuffer)))//")")
        END IF
        RESULT => me%messageBuffer(me%readPosition + 1 : me%readPosition + SIZE(byteArrayMold))
        me%readPosition = me%readPosition + SIZE(byteArrayMold)
    END FUNCTION PackedMessage_unpackBytes

    SUBROUTINE PackedMessage_printState(me)
        CLASS(t_PackedMessage), INTENT(IN) :: me

        WRITE(0,*) "t_PackedMessage{ bufferSize = "//TRIM(int2string(SIZE(me%messageBuffer)))//"; &
                   &messageSize = "//TRIM(int2string(INT(me%messageSize)))//"; &
                   &readPosition = "//TRIM(int2string(INT(me%readPosition)))//"; }"
    END SUBROUTINE PackedMessage_printState

    ! pack routines for scalar values !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !XXX: We are using MPI_Pack*() with MPI_COMM_WORLD here, because we need to be able to forward a recieved packet via a different communicator.
    !     The correct way of implementing this would be to USE MPI_Pack*_external(), but there seems to be a bug IN the MPI implementation on the cray,
    !     which renders the MPI_Pack*_external() routines unusable.

    SUBROUTINE PackedMessage_packInt(me, VALUE)
        CLASS(t_PackedMessage), INTENT(INOUT) :: me
        INTEGER, VALUE :: VALUE

        doPacking(me, VALUE)
    END SUBROUTINE PackedMessage_packInt

    SUBROUTINE PackedMessage_packSingle(me, VALUE)
        CLASS(t_PackedMessage), INTENT(INOUT) :: me
        REAL(sp), VALUE :: VALUE

        doPacking(me, VALUE)
    END SUBROUTINE PackedMessage_packSingle

    SUBROUTINE PackedMessage_packDouble(me, VALUE)
        CLASS(t_PackedMessage), INTENT(INOUT) :: me
        REAL(dp), VALUE :: VALUE

        doPacking(me, VALUE)
    END SUBROUTINE PackedMessage_packDouble

    SUBROUTINE PackedMessage_packLogical(me, VALUE)
        CLASS(t_PackedMessage), INTENT(INOUT) :: me
        LOGICAL, VALUE :: VALUE

        doPacking(me, VALUE)
    END SUBROUTINE PackedMessage_packLogical

    ! pack routines for array values !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SUBROUTINE PackedMessage_packIntArray(me, array)
        CLASS(t_PackedMessage), INTENT(INOUT) :: me
        INTEGER, ALLOCATABLE :: array(:)

        INTEGER :: arraySize, i

        arraySize = 0
        IF(ALLOCATED(array)) arraySize = SIZE(array)
        CALL me%packInt(arraySize)
        DO i = 1, arraySize
            CALL me%packInt(array(i))
        END DO
    END SUBROUTINE PackedMessage_packIntArray

    SUBROUTINE PackedMessage_packSingleArray(me, array)
        CLASS(t_PackedMessage), INTENT(INOUT) :: me
        REAL(sp), ALLOCATABLE :: array(:)

        INTEGER :: arraySize, i

        arraySize = 0
        IF(ALLOCATED(array)) arraySize = SIZE(array)
        CALL me%packInt(arraySize)
        DO i = 1, arraySize
            CALL me%packSingle(array(i))
        END DO
    END SUBROUTINE PackedMessage_packSingleArray

    SUBROUTINE PackedMessage_packDoubleArray(me, array)
        CLASS(t_PackedMessage), INTENT(INOUT) :: me
        REAL(dp), ALLOCATABLE :: array(:)

        INTEGER :: arraySize, i

        arraySize = 0
        IF(ALLOCATED(array)) arraySize = SIZE(array)
        CALL me%packInt(arraySize)
        DO i = 1, arraySize
            CALL me%packDouble(array(i))
        END DO
    END SUBROUTINE PackedMessage_packDoubleArray

    SUBROUTINE PackedMessage_packLogicalArray(me, array)
        CLASS(t_PackedMessage), INTENT(INOUT) :: me
        LOGICAL, ALLOCATABLE :: array(:)

        INTEGER :: arraySize, i

        arraySize = 0
        IF(ALLOCATED(array)) arraySize = SIZE(array)
        CALL me%packInt(arraySize)
        DO i = 1, arraySize
            CALL me%packLogical(array(i))
        END DO
    END SUBROUTINE PackedMessage_packLogicalArray

    ! unpack routines for scalar values !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SUBROUTINE PackedMessage_unpackInt(me, VALUE)
        CLASS(t_PackedMessage), INTENT(INOUT) :: me
        INTEGER, INTENT(OUT) :: VALUE

        doUnpacking(me, VALUE)
    END SUBROUTINE PackedMessage_unpackInt

    SUBROUTINE PackedMessage_unpackSingle(me, VALUE)
        CLASS(t_PackedMessage), INTENT(INOUT) :: me
        REAL(sp), INTENT(OUT) :: VALUE

        doUnpacking(me, VALUE)
    END SUBROUTINE PackedMessage_unpackSingle

    SUBROUTINE PackedMessage_unpackDouble(me, VALUE)
        CLASS(t_PackedMessage), INTENT(INOUT) :: me
        REAL(dp), INTENT(OUT) :: VALUE

        doUnpacking(me, VALUE)
    END SUBROUTINE PackedMessage_unpackDouble

    SUBROUTINE PackedMessage_unpackLogical(me, VALUE)
        CLASS(t_PackedMessage), INTENT(INOUT) :: me
        LOGICAL, INTENT(OUT) :: VALUE

        doUnpacking(me, VALUE)
    END SUBROUTINE PackedMessage_unpackLogical

    ! unpack routines for array values !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SUBROUTINE PackedMessage_unpackIntArray(me, array)
        CLASS(t_PackedMessage), INTENT(INOUT) :: me
        INTEGER, ALLOCATABLE, INTENT(INOUT) :: array(:)

        INTEGER :: arraySize, i, error
        CHARACTER(LEN = *), PARAMETER :: routine = modname//":PackedMessage_unpackIntArray"

        CALL me%unpackInt(arraySize)
        IF(ALLOCATED(array)) THEN
            IF(arraySize == 0 .OR. arraySize /= SIZE(array)) DEALLOCATE(array)
        END IF
        IF(.NOT.ALLOCATED(array) .AND. arraySize /= 0) THEN
            ALLOCATE(array(arraySize), STAT = error)
            IF(error /= SUCCESS) CALL finish(routine, "memory allocation failed")
        END IF
        DO i = 1, arraySize
            CALL me%unpackInt(array(i))
        END DO
    END SUBROUTINE PackedMessage_unpackIntArray

    SUBROUTINE PackedMessage_unpackSingleArray(me, array)
        CLASS(t_PackedMessage), INTENT(INOUT) :: me
        REAL(sp), ALLOCATABLE, INTENT(INOUT) :: array(:)

        INTEGER :: arraySize, i, error
        CHARACTER(LEN = *), PARAMETER :: routine = modname//":PackedMessage_unpackSingleArray"

        CALL me%unpackInt(arraySize)
        IF(ALLOCATED(array)) THEN
            IF(arraySize == 0 .OR. arraySize /= SIZE(array)) DEALLOCATE(array)
        END IF
        IF(.NOT.ALLOCATED(array) .AND. arraySize /= 0) THEN
            ALLOCATE(array(arraySize), STAT = error)
            IF(error /= SUCCESS) CALL finish(routine, "memory allocation failed")
        END IF
        DO i = 1, arraySize
            CALL me%unpackSingle(array(i))
        END DO
    END SUBROUTINE PackedMessage_unpackSingleArray

    SUBROUTINE PackedMessage_unpackDoubleArray(me, array)
        CLASS(t_PackedMessage), INTENT(INOUT) :: me
        REAL(dp), ALLOCATABLE, INTENT(INOUT) :: array(:)

        INTEGER :: arraySize, i, error
        CHARACTER(LEN = *), PARAMETER :: routine = modname//":PackedMessage_unpackDoubleArray"

        CALL me%unpackInt(arraySize)
        IF(ALLOCATED(array)) THEN
            IF(arraySize == 0 .OR. arraySize /= SIZE(array)) DEALLOCATE(array)
        END IF
        IF(.NOT.ALLOCATED(array) .AND. arraySize /= 0) THEN
            ALLOCATE(array(arraySize), STAT = error)
            IF(error /= SUCCESS) CALL finish(routine, "memory allocation failed")
        END IF
        DO i = 1, arraySize
            CALL me%unpackDouble(array(i))
        END DO
    END SUBROUTINE PackedMessage_unpackDoubleArray

    SUBROUTINE PackedMessage_unpackLogicalArray(me, array)
        CLASS(t_PackedMessage), INTENT(INOUT) :: me
        LOGICAL, ALLOCATABLE, INTENT(INOUT) :: array(:)

        INTEGER :: arraySize, i, error
        CHARACTER(LEN = *), PARAMETER :: routine = modname//":PackedMessage_unpackLogicalArray"

        CALL me%unpackInt(arraySize)
        IF(ALLOCATED(array)) THEN
            IF(arraySize == 0 .OR. arraySize /= SIZE(array)) DEALLOCATE(array)
        END IF
        IF(.NOT.ALLOCATED(array) .AND. arraySize /= 0) THEN
            ALLOCATE(array(arraySize), STAT = error)
            IF(error /= SUCCESS) CALL finish(routine, "memory allocation failed")
        END IF
        DO i = 1, arraySize
            CALL me%unpackLogical(array(i))
        END DO
    END SUBROUTINE PackedMessage_unpackLogicalArray

    ! wrappers for unified (un)packing !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SUBROUTINE PackedMessage_executeInt(me, operation, value)
        CLASS(t_PackedMessage), INTENT(INOUT) :: me
        INTEGER, VALUE :: operation
        INTEGER, INTENT(INOUT) :: value
        CHARACTER(LEN = *), PARAMETER :: routine = modname//":PackedMessage_executeInt"

        SELECT CASE(operation)
            CASE(kPackOp)
                CALL me%pack(value)
            CASE(kUnpackOp)
                CALL me%unpack(value)
            CASE DEFAULT
                CALL finish(routine, "illegal operation")
        END SELECT
    END SUBROUTINE PackedMessage_executeInt

    SUBROUTINE PackedMessage_executeSingle(me, operation, value)
        CLASS(t_PackedMessage), INTENT(INOUT) :: me
        INTEGER, VALUE :: operation
        REAL(sp), INTENT(INOUT) :: value
        CHARACTER(LEN = *), PARAMETER :: routine = modname//":PackedMessage_executeSingle"

        SELECT CASE(operation)
            CASE(kPackOp)
                CALL me%pack(value)
            CASE(kUnpackOp)
                CALL me%unpack(value)
            CASE DEFAULT
                CALL finish(routine, "illegal operation")
        END SELECT
    END SUBROUTINE PackedMessage_executeSingle

    SUBROUTINE PackedMessage_executeDouble(me, operation, value)
        CLASS(t_PackedMessage), INTENT(INOUT) :: me
        INTEGER, VALUE :: operation
        REAL(dp), INTENT(INOUT) :: value
        CHARACTER(LEN = *), PARAMETER :: routine = modname//":PackedMessage_executeDouble"

        SELECT CASE(operation)
            CASE(kPackOp)
                CALL me%pack(value)
            CASE(kUnpackOp)
                CALL me%unpack(value)
            CASE DEFAULT
                CALL finish(routine, "illegal operation")
        END SELECT
    END SUBROUTINE PackedMessage_executeDouble

    SUBROUTINE PackedMessage_executeLogical(me, operation, value)
        CLASS(t_PackedMessage), INTENT(INOUT) :: me
        INTEGER, VALUE :: operation
        LOGICAL, INTENT(INOUT) :: value
        CHARACTER(LEN = *), PARAMETER :: routine = modname//":PackedMessage_executeLogical"

        SELECT CASE(operation)
            CASE(kPackOp)
                CALL me%pack(value)
            CASE(kUnpackOp)
                CALL me%unpack(value)
            CASE DEFAULT
                CALL finish(routine, "illegal operation")
        END SELECT
    END SUBROUTINE PackedMessage_executeLogical

    SUBROUTINE PackedMessage_executeIntArray(me, operation, value)
        CLASS(t_PackedMessage), INTENT(INOUT) :: me
        INTEGER, VALUE :: operation
        INTEGER, ALLOCATABLE, INTENT(INOUT) :: value(:)
        CHARACTER(LEN = *), PARAMETER :: routine = modname//":PackedMessage_executeIntArray"

        SELECT CASE(operation)
            CASE(kPackOp)
                CALL me%pack(value)
            CASE(kUnpackOp)
                CALL me%unpack(value)
            CASE DEFAULT
                CALL finish(routine, "illegal operation")
        END SELECT
    END SUBROUTINE PackedMessage_executeIntArray

    SUBROUTINE PackedMessage_executeSingleArray(me, operation, value)
        CLASS(t_PackedMessage), INTENT(INOUT) :: me
        INTEGER, VALUE :: operation
        REAL(sp), ALLOCATABLE, INTENT(INOUT) :: value(:)
        CHARACTER(LEN = *), PARAMETER :: routine = modname//":PackedMessage_executeSingleArray"

        SELECT CASE(operation)
            CASE(kPackOp)
                CALL me%pack(value)
            CASE(kUnpackOp)
                CALL me%unpack(value)
            CASE DEFAULT
                CALL finish(routine, "illegal operation")
        END SELECT
    END SUBROUTINE PackedMessage_executeSingleArray

    SUBROUTINE PackedMessage_executeDoubleArray(me, operation, value)
        CLASS(t_PackedMessage), INTENT(INOUT) :: me
        INTEGER, VALUE :: operation
        REAL(dp), ALLOCATABLE, INTENT(INOUT) :: value(:)
        CHARACTER(LEN = *), PARAMETER :: routine = modname//":PackedMessage_executeDoubleArray"

        SELECT CASE(operation)
            CASE(kPackOp)
                CALL me%pack(value)
            CASE(kUnpackOp)
                CALL me%unpack(value)
            CASE DEFAULT
                CALL finish(routine, "illegal operation")
        END SELECT
    END SUBROUTINE PackedMessage_executeDoubleArray

    SUBROUTINE PackedMessage_executeLogicalArray(me, operation, value)
        CLASS(t_PackedMessage), INTENT(INOUT) :: me
        INTEGER, VALUE :: operation
        LOGICAL, ALLOCATABLE, INTENT(INOUT) :: value(:)
        CHARACTER(LEN = *), PARAMETER :: routine = modname//":PackedMessage_executeLogicalArray"

        SELECT CASE(operation)
            CASE(kPackOp)
                CALL me%pack(value)
            CASE(kUnpackOp)
                CALL me%unpack(value)
            CASE DEFAULT
                CALL finish(routine, "illegal operation")
        END SELECT
    END SUBROUTINE PackedMessage_executeLogicalArray

    ! communication routines !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! These always send/receive two MPI messages: one with the SIZE of the packed message, AND one with the actual message.
    ! This allows us to ensure that the receiver(s) always have enough memory to store the entire incoming message.

#ifndef NOMPI
#define handleMpiError(error, routine) IF(error /= MPI_SUCCESS) CALL finish(routine, "MPI error: "//TRIM(int2string(error)))
#endif

    SUBROUTINE PackedMessage_send(me, destination, tag, communicator)
        CLASS(t_PackedMessage), INTENT(IN) :: me
        INTEGER, VALUE :: destination, tag, communicator

#ifndef NOMPI
        INTEGER :: error
        CHARACTER(LEN = *), PARAMETER :: routine = modname//":PackedMessage_send"

        CALL MPI_Send(me%messageSize, 1, p_int, destination, tag, communicator, error)
        handleMpiError(error, routine)
        CALL MPI_Send(me%messageBuffer, INT(me%messageSize), MPI_BYTE, destination, tag, communicator, error)
        handleMpiError(error, routine)
#endif
    END SUBROUTINE PackedMessage_send

    SUBROUTINE PackedMessage_recv(me, source, tag, communicator)
        CLASS(t_PackedMessage), INTENT(INOUT) :: me
        INTEGER, VALUE :: source, tag, communicator

#ifndef NOMPI
        INTEGER :: error, incomingSize, status(MPI_STATUS_SIZE)
        CHARACTER(LEN = *), PARAMETER :: routine = modname//":PackedMessage_recv"

        me%messageSize = 0
        me%readPosition = 0
        CALL MPI_Recv(incomingSize, 1, p_int, source, tag, communicator, status, error)
        handleMpiError(error, routine)
        CALL me%ensureSpace(incomingSize)
        CALL MPI_Recv(me%messageBuffer, incomingSize, MPI_PACKED, source, tag, communicator, status, error)
        handleMpiError(error, routine)
        me%messageSize = incomingSize
#endif
    END SUBROUTINE PackedMessage_recv

    SUBROUTINE PackedMessage_bcast(me, root, communicator)
        CLASS(t_PackedMessage), INTENT(INOUT) :: me
        INTEGER, VALUE :: root, communicator

#ifndef NOMPI
        INTEGER :: error, incomingSize
        LOGICAL :: isRoot
        CHARACTER(LEN = *), PARAMETER :: routine = modname//":PackedMessage_bcast"

        isRoot = root == p_comm_rank(communicator)
        IF(isRoot) THEN
            incomingSize = me%messageSize
        ELSE
            me%messageSize = 0
            me%readPosition = 0
        END IF
        CALL MPI_Bcast(incomingSize, 1, p_int, root, communicator, error)
        handleMpiError(error, routine)
        IF(.NOT.isRoot) CALL me%ensureSpace(incomingSize)
        CALL MPI_Bcast(me%messageBuffer, incomingSize, MPI_PACKED, root, communicator, error)
        handleMpiError(error, routine)
        me%messageSize = incomingSize
#endif
    END SUBROUTINE PackedMessage_bcast

    ! destructor !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SUBROUTINE PackedMessage_destruct(me)
        CLASS(t_PackedMessage), INTENT(INOUT) :: me

        DEALLOCATE(me%messageBuffer)
    END SUBROUTINE PackedMessage_destruct

END MODULE mo_packed_message
