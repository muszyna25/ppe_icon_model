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
    USE mo_fortran_tools, ONLY: t_Destructible
    USE mo_impl_constants, ONLY: SUCCESS
    USE mo_kind, ONLY: sp, dp, i8
    USE mo_util_string, ONLY: int2string

#ifndef NOMPI
    USE mo_mpi, ONLY: p_int, p_get_bcast_role, MPI_SUCCESS
    USE mpi
#endif

    IMPLICIT NONE
    PRIVATE

    PUBLIC :: t_PackedMessage
    INTEGER, PARAMETER, PUBLIC :: kPackOp = 0, kUnpackOp = 1

    ! A t_PackedMessage IS used to bundle a number of different values
    ! together into a single message, that can be communicated via a
    ! single CALL.  It IS possible to have ANY number of communication
    ! steps between the packing AND unpacking, including zero (a PE
    ! unpacks its own DATA), AND two (a PE recieves a packed message
    ! AND passes it on, possibly via a different communicator).
    !
    ! If NOMPI IS defined, the communication routines are simply
    ! noops, the packing/unpacking still works as expected.
    !
    ! As an added bonus, this provides packerXXX() routines IN
    ! addition to the packXXX() AND unpackXXX() routines, which allow
    ! folding the packing AND unpacking into the same code. Ie.,
    ! instead of writing a routine containing
    !
    !   message%pack(foo)
    !   message%pack(bar)
    !   message%pack(baz)
    !
    ! AND a second routine containing
    !
    !   message%unpack(foo)
    !   message%unpack(baz) !Error: messed up sequence!
    !   message%unpack(bar) !Error: messed up sequence!
    !
    ! you can WRITE a single routine containing
    !
    !   message%packer(operation, foo)
    !   message%packer(operation, baz)
    !   message%packer(operation, bar)
    !
    ! knowing that it will be simply impossible to mix up the sequence
    ! when setting operation to kUnpackOp to unpack the message.
    !
    ! XXX: This originated as a wrapper around MPI_Pack() AND friends
    ! that IS able to manage the buffer that's used to hold the packed
    ! message.  However, it turned OUT to be more sensible to DO the
    ! packing ourselves: We need to be able to pack/unpack even when
    ! NOMPI IS defined.

    TYPE, EXTENDS(t_Destructible) :: t_PackedMessage
        INTEGER messageSize, readPosition
        CHARACTER, POINTER :: messageBuffer(:)
      CONTAINS
        ! functionally equivalent to CALL message%destruct(); CALL
        ! message%construct(), but more efficient (buffer IS reused)
        PROCEDURE :: construct => PackedMessage_construct
        PROCEDURE :: reset => PackedMessage_reset   

        ! pack routines
        PROCEDURE :: packInt => PackedMessage_packInt
        PROCEDURE :: packLong => PackedMessage_packLong        
        PROCEDURE :: packSingle => PackedMessage_packSingle
        PROCEDURE :: packDouble => PackedMessage_packDouble
        PROCEDURE :: packLogical => PackedMessage_packLogical
        PROCEDURE :: packCharacter => PackedMessage_packCharacter

        PROCEDURE :: packAllocatableInt => PackedMessage_packAllocatableInt
        PROCEDURE :: packAllocatableLong => PackedMessage_packAllocatableLong        
        PROCEDURE :: packAllocatableSingle => PackedMessage_packAllocatableSingle
        PROCEDURE :: packAllocatableDouble => PackedMessage_packAllocatableDouble
        PROCEDURE :: packAllocatableLogical => PackedMessage_packAllocatableLogical
        !XXX: see PackedMessage_unpackAllocatableCharacter() for an explanation of the deactivation
        ! PROCEDURE :: packAllocatableCharacter => PackedMessage_packAllocatableCharacter

        PROCEDURE :: packIntArray => PackedMessage_packIntArray
        PROCEDURE :: packLongArray => PackedMessage_packLongArray        
        PROCEDURE :: packSingleArray => PackedMessage_packSingleArray
        PROCEDURE :: packDoubleArray => PackedMessage_packDoubleArray
        PROCEDURE :: packLogicalArray => PackedMessage_packLogicalArray

        GENERIC :: pack => packInt, packLong, &
             &             packSingle, packDouble, &
             &             packLogical, packCharacter, &
             &             packIntArray, packLongArray, &
             &             packSingleArray, packDoubleArray, &
             &             packLogicalArray
        GENERIC :: packAllocatable => packAllocatableInt, packAllocatableLong, &
             &                        packAllocatableSingle, packAllocatableDouble, &
             &                        packAllocatableLogical

        ! unpack routines
        PROCEDURE :: unpackInt => PackedMessage_unpackInt
        PROCEDURE :: unpackLong => PackedMessage_unpackLong        
        PROCEDURE :: unpackSingle => PackedMessage_unpackSingle
        PROCEDURE :: unpackDouble => PackedMessage_unpackDouble
        PROCEDURE :: unpackLogical => PackedMessage_unpackLogical
        PROCEDURE :: unpackCharacter => PackedMessage_unpackCharacter

        PROCEDURE :: unpackAllocatableInt => PackedMessage_unpackAllocatableInt
        PROCEDURE :: unpackAllocatableLong => PackedMessage_unpackAllocatableLong        
        PROCEDURE :: unpackAllocatableSingle => PackedMessage_unpackAllocatableSingle
        PROCEDURE :: unpackAllocatableDouble => PackedMessage_unpackAllocatableDouble
        PROCEDURE :: unpackAllocatableLogical => PackedMessage_unpackAllocatableLogical
        !XXX: see PackedMessage_unpackAllocatableCharacter() for an explanation of the deactivation
        ! PROCEDURE :: unpackAllocatableCharacter => PackedMessage_unpackAllocatableCharacter

        PROCEDURE :: unpackIntArray => PackedMessage_unpackIntArray
        PROCEDURE :: unpackLongArray => PackedMessage_unpackLongArray        
        PROCEDURE :: unpackSingleArray => PackedMessage_unpackSingleArray
        PROCEDURE :: unpackDoubleArray => PackedMessage_unpackDoubleArray
        PROCEDURE :: unpackLogicalArray => PackedMessage_unpackLogicalArray

        GENERIC :: unpack => unpackInt, unpackLong, &
             &               unpackSingle, unpackDouble, &
             &               unpackLogical, unpackCharacter, &
             &               unpackIntArray, unpackLongArray, &
             &               unpackSingleArray, unpackDoubleArray, &
             &               unpackLogicalArray
        GENERIC :: unpackAllocatable => unpackAllocatableInt, unpackAllocatableLong, &
             &                          unpackAllocatableSingle, unpackAllocatableDouble, &
             &                          unpackAllocatableLogical

        ! routines to facilitate packing AND unpacking with the same code
        ! USE of these routine prohibits ANY errors by mismatches between packing AND unpacking code
        PROCEDURE :: packerInt => PackedMessage_packerInt
        PROCEDURE :: packerLong => PackedMessage_packerLong        
        PROCEDURE :: packerSingle => PackedMessage_packerSingle
        PROCEDURE :: packerDouble => PackedMessage_packerDouble
        PROCEDURE :: packerLogical => PackedMessage_packerLogical
        PROCEDURE :: packerCharacter => PackedMessage_packerCharacter

        PROCEDURE :: packerAllocatableInt => PackedMessage_packerAllocatableInt
        PROCEDURE :: packerAllocatableLong => PackedMessage_packerAllocatableLong        
        PROCEDURE :: packerAllocatableSingle => PackedMessage_packerAllocatableSingle
        PROCEDURE :: packerAllocatableDouble => PackedMessage_packerAllocatableDouble
        PROCEDURE :: packerAllocatableLogical => PackedMessage_packerAllocatableLogical
        !XXX: see PackedMessage_unpackAllocatableCharacter() for an explanation of the deactivation
        ! PROCEDURE :: packerAllocatableCharacter => PackedMessage_packerAllocatableCharacter

        PROCEDURE :: packerIntArray => PackedMessage_packerIntArray
        PROCEDURE :: packerLongArray => PackedMessage_packerLongArray        
        PROCEDURE :: packerSingleArray => PackedMessage_packerSingleArray
        PROCEDURE :: packerDoubleArray => PackedMessage_packerDoubleArray
        PROCEDURE :: packerLogicalArray => PackedMessage_packerLogicalArray

        GENERIC :: packer => packerInt, packerLong, &
             &               packerSingle, packerDouble, &
             &               packerLogical, packerCharacter, &
             &               packerIntArray, packerLongArray, &
             &               packerSingleArray, packerDoubleArray, &
             &               packerLogicalArray
        GENERIC :: packerAllocatable => packerAllocatableInt, packerAllocatableLong, &
             &                          packerAllocatableSingle, packerAllocatableDouble, &
             &                          packerAllocatableLogical

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
        ! returns a POINTER to a subarray of the internal storage, so
        ! that POINTER IS invalidated by ANY (IN)direct CALL to
        ! ensureSpace()!
        PROCEDURE, PRIVATE :: unpackBytes => PackedMessage_unpackBytes  

#define doPacking(me, VALUE) CALL me%packBytes(TRANSFER(VALUE, characterArrayMold))
#define doUnpacking(me, VALUE) VALUE = TRANSFER(me%unpackBytes(TRANSFER(VALUE, characterArrayMold)), VALUE)

    END TYPE t_PackedMessage

    CHARACTER :: characterMold, characterArrayMold(1)
    INTEGER :: integerMold
    REAL(sp) :: singleMold
    REAL(dp) :: doubleMold
    LOGICAL :: logicalMold

    CHARACTER(len=*), PARAMETER :: modname = "mo_packed_message"

CONTAINS

    SUBROUTINE PackedMessage_construct(me)
        CLASS(t_PackedMessage), INTENT(OUT) :: me

        INTEGER :: error
        CHARACTER(len=*), PARAMETER :: routine = modname//":PackedMessage_construct"

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

    FUNCTION PackedMessage_unpackBytes(me, byteArrayMold) RESULT(resultVar)
        CLASS(t_PackedMessage), INTENT(INOUT) :: me
        CHARACTER, INTENT(IN) :: byteArrayMold(:)
        CHARACTER, POINTER :: resultVar(:)

        CHARACTER(LEN = *), PARAMETER :: routine = modname//":PackedMessage_unpackBytes"

        IF(me%readPosition + SIZE(byteArrayMold) > SIZE(me%messageBuffer)) THEN
            CALL finish(routine, "attempt to read more data from a t_PackedMessage than it contains (&
                                 &zero based read position = "//TRIM(int2string(me%readPosition))//", &
                                 &read size = "//TRIM(int2string(SIZE(byteArrayMold)))//", &
                                 &message size = "//TRIM(int2string(SIZE(me%messageBuffer)))//")")
        END IF
        resultVar => me%messageBuffer(me%readPosition + 1 : me%readPosition + SIZE(byteArrayMold))
        me%readPosition = me%readPosition + SIZE(byteArrayMold)
    END FUNCTION PackedMessage_unpackBytes

    SUBROUTINE PackedMessage_printState(me)
        CLASS(t_PackedMessage), INTENT(IN) :: me

        WRITE(0,*) "t_PackedMessage{ bufferSize = "//TRIM(int2string(SIZE(me%messageBuffer)))//"; &
                   &messageSize = "//TRIM(int2string(INT(me%messageSize)))//"; &
                   &readPosition = "//TRIM(int2string(INT(me%readPosition)))//"; }"
    END SUBROUTINE PackedMessage_printState

    ! pack routines for scalar values !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !XXX: Originally, this was designed to USE the MPI_Pack*() routines. However, this had some problems:
    !
    !       * The MPI_Pack*() routines require that the communicator IS passed already when packing, NOT just during communication.
    !         In theory, this would stop us from passing a message received via one communicator on via another communicator.
    !         While current MPI implementations don't seem to object such an abuse, it nevertheless invokes undefined behavior, AND may break with ANY MPI update.
    !         So, since we need to be able to pass on messages through different communicators, MPI_Pack*() routines are OUT of the game.
    !
    !       * The MPI_Pack*_external() routines don't have the problem above, however there IS a bug IN their implementation on the cray which renders them unusable.
    !
    !       * Since ICON can be built without MPI, AND since we nevertheless need to be able to pack/unpack a message as a local operation, we would need to provide an MPI-free fallback implementation.
    !
    !     Since the required fallback implementation IS sufficient to be also used IN the MPI CASE, I have scraped the USE of MPI_Pack*() AND MPI_Pack*_external() from this MODULE.

    SUBROUTINE PackedMessage_packInt(me, VALUE)
        CLASS(t_PackedMessage), INTENT(INOUT) :: me
        INTEGER, VALUE :: VALUE

        doPacking(me, VALUE)
    END SUBROUTINE PackedMessage_packInt

    SUBROUTINE PackedMessage_packLong(me, VALUE)
        CLASS(t_PackedMessage), INTENT(INOUT) :: me
        INTEGER(i8), VALUE :: VALUE

        doPacking(me, VALUE)
    END SUBROUTINE PackedMessage_packLong

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

    SUBROUTINE PackedMessage_packCharacter(me, VALUE)
        CLASS(t_PackedMessage), INTENT(INOUT) :: me
        CHARACTER(LEN = *), INTENT(IN) :: VALUE

        doPacking(me, VALUE)
    END SUBROUTINE PackedMessage_packCharacter

    ! pack routines for ALLOCATABLE scalars !

    SUBROUTINE PackedMessage_packAllocatableInt(me, VALUE)
        CLASS(t_PackedMessage), INTENT(INOUT) :: me
        INTEGER, ALLOCATABLE, INTENT(IN) :: VALUE

        CALL me%pack(ALLOCATED(VALUE))
        IF(ALLOCATED(VALUE)) CALL me%pack(VALUE)
    END SUBROUTINE PackedMessage_packAllocatableInt

    SUBROUTINE PackedMessage_packAllocatableLong(me, VALUE)
        CLASS(t_PackedMessage), INTENT(INOUT) :: me
        INTEGER(i8), ALLOCATABLE, INTENT(IN) :: VALUE

        CALL me%pack(ALLOCATED(VALUE))
        IF(ALLOCATED(VALUE)) CALL me%pack(VALUE)
    END SUBROUTINE PackedMessage_packAllocatableLong

    SUBROUTINE PackedMessage_packAllocatableSingle(me, VALUE)
        CLASS(t_PackedMessage), INTENT(INOUT) :: me
        REAL(sp), ALLOCATABLE, INTENT(IN) :: VALUE

        CALL me%pack(ALLOCATED(VALUE))
        IF(ALLOCATED(VALUE)) CALL me%pack(VALUE)
    END SUBROUTINE PackedMessage_packAllocatableSingle

    SUBROUTINE PackedMessage_packAllocatableDouble(me, VALUE)
        CLASS(t_PackedMessage), INTENT(INOUT) :: me
        REAL(dp), ALLOCATABLE, INTENT(IN) :: VALUE

        CALL me%pack(ALLOCATED(VALUE))
        IF(ALLOCATED(VALUE)) CALL me%pack(VALUE)
    END SUBROUTINE PackedMessage_packAllocatableDouble

    SUBROUTINE PackedMessage_packAllocatableLogical(me, VALUE)
        CLASS(t_PackedMessage), INTENT(INOUT) :: me
        LOGICAL, ALLOCATABLE, INTENT(IN) :: VALUE

        CALL me%pack(ALLOCATED(VALUE))
        IF(ALLOCATED(VALUE)) CALL me%pack(VALUE)
    END SUBROUTINE PackedMessage_packAllocatableLogical

!XXX: see PackedMessage_unpackAllocatableCharacter() for an explanation of the deactivation
!   SUBROUTINE PackedMessage_packAllocatableCharacter(me, VALUE)
!       CLASS(t_PackedMessage), INTENT(INOUT) :: me
!       CHARACTER(:), ALLOCATABLE, INTENT(IN) :: VALUE

!       CALL me%pack(ALLOCATED(VALUE))
!       IF(ALLOCATED(VALUE)) THEN
!           CALL me%pack(LEN(VALUE))
!           CALL me%pack(VALUE)
!       END IF
!   END SUBROUTINE PackedMessage_packAllocatableCharacter

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

    SUBROUTINE PackedMessage_packLongArray(me, array)
        CLASS(t_PackedMessage), INTENT(INOUT) :: me
        INTEGER(i8), ALLOCATABLE :: array(:)

        INTEGER :: arraySize, i

        arraySize = 0
        IF(ALLOCATED(array)) arraySize = SIZE(array)
        CALL me%packInt(arraySize)
        DO i = 1, arraySize
            CALL me%packLong(array(i))
        END DO
      END SUBROUTINE PackedMessage_packLongArray

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

    SUBROUTINE PackedMessage_unpackLong(me, VALUE)
        CLASS(t_PackedMessage), INTENT(INOUT) :: me
        INTEGER(i8), INTENT(OUT) :: VALUE

        doUnpacking(me, VALUE)
    END SUBROUTINE PackedMessage_unpackLong

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

    SUBROUTINE PackedMessage_unpackCharacter(me, VALUE)
        CLASS(t_PackedMessage), INTENT(INOUT) :: me
        CHARACTER(LEN = *), INTENT(INOUT) :: VALUE

        doUnpacking(me, VALUE)
    END SUBROUTINE PackedMessage_unpackCharacter

    ! unpack routines for ALLOCATABLE scalars !

    SUBROUTINE PackedMessage_unpackAllocatableInt(me, VALUE)
        CLASS(t_PackedMessage), INTENT(INOUT) :: me
        INTEGER, ALLOCATABLE, INTENT(INOUT) :: VALUE

        LOGICAL :: isAllocated
        INTEGER :: error
        CHARACTER(*), PARAMETER :: routine = modname//":PackedMessage_unpackAllocatableInt"

        CALL me%unpack(isAllocated)
        IF(isAllocated) THEN
            IF(.NOT.ALLOCATED(VALUE)) THEN
                ALLOCATE(VALUE, STAT = error)
                IF(error /= SUCCESS) CALL finish(routine, "memory allocation error")
            END IF
            CALL me%unpack(VALUE)
        ELSE
            IF(ALLOCATED(VALUE)) DEALLOCATE(VALUE)
        END IF
    END SUBROUTINE PackedMessage_unpackAllocatableInt

    SUBROUTINE PackedMessage_unpackAllocatableLong(me, VALUE)
        CLASS(t_PackedMessage), INTENT(INOUT) :: me
        INTEGER(i8), ALLOCATABLE, INTENT(INOUT) :: VALUE

        LOGICAL :: isAllocated
        INTEGER :: error
        CHARACTER(*), PARAMETER :: routine = modname//":PackedMessage_unpackAllocatableLong"

        CALL me%unpack(isAllocated)
        IF(isAllocated) THEN
            IF(.NOT.ALLOCATED(VALUE)) THEN
                ALLOCATE(VALUE, STAT = error)
                IF(error /= SUCCESS) CALL finish(routine, "memory allocation error")
            END IF
            CALL me%unpack(VALUE)
        ELSE
            IF(ALLOCATED(VALUE)) DEALLOCATE(VALUE)
        END IF
    END SUBROUTINE PackedMessage_unpackAllocatableLong

    SUBROUTINE PackedMessage_unpackAllocatableSingle(me, VALUE)
        CLASS(t_PackedMessage), INTENT(INOUT) :: me
        REAL(sp), ALLOCATABLE, INTENT(INOUT) :: VALUE

        LOGICAL :: isAllocated
        INTEGER :: error
        CHARACTER(*), PARAMETER :: routine = modname//":PackedMessage_unpackAllocatableSingle"

        CALL me%unpack(isAllocated)
        IF(isAllocated) THEN
            IF(.NOT.ALLOCATED(VALUE)) THEN
                ALLOCATE(VALUE, STAT = error)
                IF(error /= SUCCESS) CALL finish(routine, "memory allocation error")
            END IF
            CALL me%unpack(VALUE)
        ELSE
            IF(ALLOCATED(VALUE)) DEALLOCATE(VALUE)
        END IF
    END SUBROUTINE PackedMessage_unpackAllocatableSingle

    SUBROUTINE PackedMessage_unpackAllocatableDouble(me, VALUE)
        CLASS(t_PackedMessage), INTENT(INOUT) :: me
        REAL(dp), ALLOCATABLE, INTENT(INOUT) :: VALUE

        LOGICAL :: isAllocated
        INTEGER :: error
        CHARACTER(*), PARAMETER :: routine = modname//":PackedMessage_unpackAllocatableDouble"

        CALL me%unpack(isAllocated)
        IF(isAllocated) THEN
            IF(.NOT.ALLOCATED(VALUE)) THEN
                ALLOCATE(VALUE, STAT = error)
                IF(error /= SUCCESS) CALL finish(routine, "memory allocation error")
            END IF
            CALL me%unpack(VALUE)
        ELSE
            IF(ALLOCATED(VALUE)) DEALLOCATE(VALUE)
        END IF
    END SUBROUTINE PackedMessage_unpackAllocatableDouble

    SUBROUTINE PackedMessage_unpackAllocatableLogical(me, VALUE)
        CLASS(t_PackedMessage), INTENT(INOUT) :: me
        LOGICAL, ALLOCATABLE, INTENT(INOUT) :: VALUE

        LOGICAL :: isAllocated
        INTEGER :: error
        CHARACTER(*), PARAMETER :: routine = modname//":PackedMessage_unpackAllocatableLogical"

        CALL me%unpack(isAllocated)
        IF(isAllocated) THEN
            IF(.NOT.ALLOCATED(VALUE)) THEN
                ALLOCATE(VALUE, STAT = error)
                IF(error /= SUCCESS) CALL finish(routine, "memory allocation error")
            END IF
            CALL me%unpack(VALUE)
        ELSE
            IF(ALLOCATED(VALUE)) DEALLOCATE(VALUE)
        END IF
    END SUBROUTINE PackedMessage_unpackAllocatableLogical

!XXX: This is deactivated because it triggers a bug in gfortran that looses the contents of VALUE on return from the routine.
!     The only known workaround is to implement the functionality right where it is needed.
!
!   SUBROUTINE PackedMessage_unpackAllocatableCharacter(me, VALUE)
!       CLASS(t_PackedMessage), INTENT(INOUT) :: me
!       CHARACTER(:), ALLOCATABLE, INTENT(INOUT) :: VALUE

!       LOGICAL :: isAllocated
!       INTEGER :: error, length
!       CHARACTER(*), PARAMETER :: routine = modname//":PackedMessage_unpackAllocatableCharacter"

!       IF(ALLOCATED(VALUE)) DEALLOCATE(VALUE)  ! get VALUE into a known state (unallocated)

!       CALL me%unpack(isAllocated)
!       IF(isAllocated) THEN
!           CALL me%unpack(length)
!           ALLOCATE(CHARACTER(length) :: VALUE, STAT = error)
!           IF(error /= SUCCESS) CALL finish(routine, "memory allocation error")
!           CALL me%unpack(VALUE)
!       END IF
!   END SUBROUTINE PackedMessage_unpackAllocatableCharacter

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

    SUBROUTINE PackedMessage_unpackLongArray(me, array)
        CLASS(t_PackedMessage), INTENT(INOUT) :: me
        INTEGER(i8), ALLOCATABLE, INTENT(INOUT) :: array(:)

        INTEGER :: arraySize, i, error
        CHARACTER(LEN = *), PARAMETER :: routine = modname//":PackedMessage_unpackLongArray"

        CALL me%unpackInt(arraySize)
        IF(ALLOCATED(array)) THEN
            IF(arraySize == 0 .OR. arraySize /= SIZE(array)) DEALLOCATE(array)
        END IF
        IF(.NOT.ALLOCATED(array) .AND. arraySize /= 0) THEN
            ALLOCATE(array(arraySize), STAT = error)
            IF(error /= SUCCESS) CALL finish(routine, "memory allocation failed")
        END IF
        DO i = 1, arraySize
            CALL me%unpackLong(array(i))
        END DO
    END SUBROUTINE PackedMessage_unpackLongArray

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

    SUBROUTINE PackedMessage_packerInt(me, operation, value)
        CLASS(t_PackedMessage), INTENT(INOUT) :: me
        INTEGER, VALUE :: operation
        INTEGER, INTENT(INOUT) :: value
        CHARACTER(LEN = *), PARAMETER :: routine = modname//":PackedMessage_packerInt"

        SELECT CASE(operation)
            CASE(kPackOp)
                CALL me%pack(value)
            CASE(kUnpackOp)
                CALL me%unpack(value)
            CASE DEFAULT
                CALL finish(routine, "illegal operation")
        END SELECT
    END SUBROUTINE PackedMessage_packerInt

    SUBROUTINE PackedMessage_packerLong(me, operation, value)
        CLASS(t_PackedMessage), INTENT(INOUT) :: me
        INTEGER, VALUE :: operation
        INTEGER(i8), INTENT(INOUT) :: value
        CHARACTER(LEN = *), PARAMETER :: routine = modname//":PackedMessage_packerInt"

        SELECT CASE(operation)
            CASE(kPackOp)
                CALL me%pack(value)
            CASE(kUnpackOp)
                CALL me%unpack(value)
            CASE DEFAULT
                CALL finish(routine, "illegal operation")
        END SELECT
    END SUBROUTINE PackedMessage_packerLong

    SUBROUTINE PackedMessage_packerSingle(me, operation, value)
        CLASS(t_PackedMessage), INTENT(INOUT) :: me
        INTEGER, VALUE :: operation
        REAL(sp), INTENT(INOUT) :: value
        CHARACTER(LEN = *), PARAMETER :: routine = modname//":PackedMessage_packerSingle"

        SELECT CASE(operation)
            CASE(kPackOp)
                CALL me%pack(value)
            CASE(kUnpackOp)
                CALL me%unpack(value)
            CASE DEFAULT
                CALL finish(routine, "illegal operation")
        END SELECT
    END SUBROUTINE PackedMessage_packerSingle

    SUBROUTINE PackedMessage_packerDouble(me, operation, value)
        CLASS(t_PackedMessage), INTENT(INOUT) :: me
        INTEGER, VALUE :: operation
        REAL(dp), INTENT(INOUT) :: value
        CHARACTER(LEN = *), PARAMETER :: routine = modname//":PackedMessage_packerDouble"

        SELECT CASE(operation)
            CASE(kPackOp)
                CALL me%pack(value)
            CASE(kUnpackOp)
                CALL me%unpack(value)
            CASE DEFAULT
                CALL finish(routine, "illegal operation")
        END SELECT
    END SUBROUTINE PackedMessage_packerDouble

    SUBROUTINE PackedMessage_packerLogical(me, operation, value)
        CLASS(t_PackedMessage), INTENT(INOUT) :: me
        INTEGER, VALUE :: operation
        LOGICAL, INTENT(INOUT) :: value
        CHARACTER(LEN = *), PARAMETER :: routine = modname//":PackedMessage_packerLogical"

        SELECT CASE(operation)
            CASE(kPackOp)
                CALL me%pack(value)
            CASE(kUnpackOp)
                CALL me%unpack(value)
            CASE DEFAULT
                CALL finish(routine, "illegal operation")
        END SELECT
    END SUBROUTINE PackedMessage_packerLogical

    SUBROUTINE PackedMessage_packerCharacter(me, operation, value)
        CLASS(t_PackedMessage), INTENT(INOUT) :: me
        INTEGER, VALUE :: operation
        CHARACTER(LEN = *), INTENT(INOUT) :: VALUE
        CHARACTER(LEN = *), PARAMETER :: routine = modname//":PackedMessage_packerCharacter"

        SELECT CASE(operation)
            CASE(kPackOp)
                CALL me%pack(value)
            CASE(kUnpackOp)
                CALL me%unpack(value)
            CASE DEFAULT
                CALL finish(routine, "illegal operation")
        END SELECT
    END SUBROUTINE PackedMessage_packerCharacter

    SUBROUTINE PackedMessage_packerAllocatableInt(me, operation, VALUE)
        CLASS(t_PackedMessage), INTENT(INOUT) :: me
        INTEGER, VALUE :: operation
        INTEGER, ALLOCATABLE, INTENT(INOUT) :: VALUE

        CHARACTER(*), PARAMETER :: routine = modname//":PackedMessage_packerAllocatableInt"

        SELECT CASE(operation)
            CASE(kPackOp)
                CALL me%packAllocatable(VALUE)
            CASE(kUnpackOp)
                CALL me%unpackAllocatable(VALUE)
            CASE DEFAULT
                CALL finish(routine, "illegal operation")
        END SELECT
    END SUBROUTINE PackedMessage_packerAllocatableInt

    SUBROUTINE PackedMessage_packerAllocatableLong(me, operation, VALUE)
        CLASS(t_PackedMessage), INTENT(INOUT) :: me
        INTEGER, VALUE :: operation
        INTEGER(i8), ALLOCATABLE, INTENT(INOUT) :: VALUE

        CHARACTER(*), PARAMETER :: routine = modname//":PackedMessage_packerAllocatableLong"

        SELECT CASE(operation)
            CASE(kPackOp)
                CALL me%packAllocatable(VALUE)
            CASE(kUnpackOp)
                CALL me%unpackAllocatable(VALUE)
            CASE DEFAULT
                CALL finish(routine, "illegal operation")
        END SELECT
    END SUBROUTINE PackedMessage_packerAllocatableLong

    SUBROUTINE PackedMessage_packerAllocatableSingle(me, operation, VALUE)
        CLASS(t_PackedMessage), INTENT(INOUT) :: me
        INTEGER, VALUE :: operation
        REAL(sp), ALLOCATABLE, INTENT(INOUT) :: VALUE

        CHARACTER(*), PARAMETER :: routine = modname//":PackedMessage_packerAllocatableSingle"

        SELECT CASE(operation)
            CASE(kPackOp)
                CALL me%packAllocatable(VALUE)
            CASE(kUnpackOp)
                CALL me%unpackAllocatable(VALUE)
            CASE DEFAULT
                CALL finish(routine, "illegal operation")
        END SELECT
    END SUBROUTINE PackedMessage_packerAllocatableSingle

    SUBROUTINE PackedMessage_packerAllocatableDouble(me, operation, VALUE)
        CLASS(t_PackedMessage), INTENT(INOUT) :: me
        INTEGER, VALUE :: operation
        REAL(dp), ALLOCATABLE, INTENT(INOUT) :: VALUE

        CHARACTER(*), PARAMETER :: routine = modname//":PackedMessage_packerAllocatableDouble"

        SELECT CASE(operation)
            CASE(kPackOp)
                CALL me%packAllocatable(VALUE)
            CASE(kUnpackOp)
                CALL me%unpackAllocatable(VALUE)
            CASE DEFAULT
                CALL finish(routine, "illegal operation")
        END SELECT
    END SUBROUTINE PackedMessage_packerAllocatableDouble

    SUBROUTINE PackedMessage_packerAllocatableLogical(me, operation, VALUE)
        CLASS(t_PackedMessage), INTENT(INOUT) :: me
        INTEGER, VALUE :: operation
        LOGICAL, ALLOCATABLE, INTENT(INOUT) :: VALUE

        CHARACTER(*), PARAMETER :: routine = modname//":PackedMessage_packerAllocatableLogical"

        SELECT CASE(operation)
            CASE(kPackOp)
                CALL me%packAllocatable(VALUE)
            CASE(kUnpackOp)
                CALL me%unpackAllocatable(VALUE)
            CASE DEFAULT
                CALL finish(routine, "illegal operation")
        END SELECT
    END SUBROUTINE PackedMessage_packerAllocatableLogical

!XXX: see PackedMessage_unpackAllocatableCharacter() for an explanation of the deactivation
!   SUBROUTINE PackedMessage_packerAllocatableCharacter(me, operation, VALUE)
!       CLASS(t_PackedMessage), INTENT(INOUT) :: me
!       INTEGER, VALUE :: operation
!       CHARACTER(:), ALLOCATABLE, INTENT(INOUT) :: VALUE

!       CHARACTER(*), PARAMETER :: routine = modname//":PackedMessage_packerAllocatableCharacter"

!       SELECT CASE(operation)
!           CASE(kPackOp)
!               CALL me%packAllocatable(VALUE)
!           CASE(kUnpackOp)
!               CALL me%unpackAllocatable(VALUE)
!           CASE DEFAULT
!               CALL finish(routine, "illegal operation")
!       END SELECT
!   END SUBROUTINE PackedMessage_packerAllocatableCharacter

    SUBROUTINE PackedMessage_packerIntArray(me, operation, value)
        CLASS(t_PackedMessage), INTENT(INOUT) :: me
        INTEGER, VALUE :: operation
        INTEGER, ALLOCATABLE, INTENT(INOUT) :: value(:)
        CHARACTER(LEN = *), PARAMETER :: routine = modname//":PackedMessage_packerIntArray"

        SELECT CASE(operation)
            CASE(kPackOp)
                CALL me%pack(value)
            CASE(kUnpackOp)
                CALL me%unpack(value)
            CASE DEFAULT
                CALL finish(routine, "illegal operation")
        END SELECT
    END SUBROUTINE PackedMessage_packerIntArray
    
    SUBROUTINE PackedMessage_packerLongArray(me, operation, value)
        CLASS(t_PackedMessage), INTENT(INOUT) :: me
        INTEGER, VALUE :: operation
        INTEGER(i8), ALLOCATABLE, INTENT(INOUT) :: value(:)
        CHARACTER(LEN = *), PARAMETER :: routine = modname//":PackedMessage_packerLongArray"

        SELECT CASE(operation)
            CASE(kPackOp)
                CALL me%pack(value)
            CASE(kUnpackOp)
                CALL me%unpack(value)
            CASE DEFAULT
                CALL finish(routine, "illegal operation")
        END SELECT
    END SUBROUTINE PackedMessage_packerLongArray

    SUBROUTINE PackedMessage_packerSingleArray(me, operation, value)
        CLASS(t_PackedMessage), INTENT(INOUT) :: me
        INTEGER, VALUE :: operation
        REAL(sp), ALLOCATABLE, INTENT(INOUT) :: value(:)
        CHARACTER(LEN = *), PARAMETER :: routine = modname//":PackedMessage_packerSingleArray"

        SELECT CASE(operation)
            CASE(kPackOp)
                CALL me%pack(value)
            CASE(kUnpackOp)
                CALL me%unpack(value)
            CASE DEFAULT
                CALL finish(routine, "illegal operation")
        END SELECT
    END SUBROUTINE PackedMessage_packerSingleArray

    SUBROUTINE PackedMessage_packerDoubleArray(me, operation, value)
        CLASS(t_PackedMessage), INTENT(INOUT) :: me
        INTEGER, VALUE :: operation
        REAL(dp), ALLOCATABLE, INTENT(INOUT) :: value(:)
        CHARACTER(LEN = *), PARAMETER :: routine = modname//":PackedMessage_packerDoubleArray"

        SELECT CASE(operation)
            CASE(kPackOp)
                CALL me%pack(value)
            CASE(kUnpackOp)
                CALL me%unpack(value)
            CASE DEFAULT
                CALL finish(routine, "illegal operation")
        END SELECT
    END SUBROUTINE PackedMessage_packerDoubleArray

    SUBROUTINE PackedMessage_packerLogicalArray(me, operation, value)
        CLASS(t_PackedMessage), INTENT(INOUT) :: me
        INTEGER, VALUE :: operation
        LOGICAL, ALLOCATABLE, INTENT(INOUT) :: value(:)
        CHARACTER(LEN = *), PARAMETER :: routine = modname//":PackedMessage_packerLogicalArray"

        SELECT CASE(operation)
            CASE(kPackOp)
                CALL me%pack(value)
            CASE(kUnpackOp)
                CALL me%unpack(value)
            CASE DEFAULT
                CALL finish(routine, "illegal operation")
        END SELECT
    END SUBROUTINE PackedMessage_packerLogicalArray

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
        LOGICAL :: isRoot, isReceiver
        CHARACTER(LEN = *), PARAMETER :: routine = modname//":PackedMessage_bcast"

        ! reset the message IF this IS a receiving process
        CALL p_get_bcast_role(root, communicator, isRoot, isReceiver)
        IF(isReceiver) CALL me%reset()

        ! broadcast the message SIZE to ensure enough buffer space IN the receiving processes
        incomingSize = me%messageSize
        CALL MPI_Bcast(incomingSize, 1, p_int, root, communicator, error)
        handleMpiError(error, routine)
        IF(isReceiver) CALL me%ensureSpace(incomingSize)

        ! broadcast the message itself
        CALL MPI_Bcast(me%messageBuffer, incomingSize, MPI_PACKED, root, communicator, error)
        handleMpiError(error, routine)
        IF(isReceiver) me%messageSize = incomingSize
#endif
    END SUBROUTINE PackedMessage_bcast

    ! destructor !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SUBROUTINE PackedMessage_destruct(me)
        CLASS(t_PackedMessage), INTENT(INOUT) :: me

        DEALLOCATE(me%messageBuffer)
    END SUBROUTINE PackedMessage_destruct

END MODULE mo_packed_message
