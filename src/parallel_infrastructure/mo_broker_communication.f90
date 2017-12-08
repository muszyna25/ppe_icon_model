!> Module for reading multifile restart files
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.

MODULE mo_broker_communication
    USE mo_exception, ONLY: finish
    USE mo_util_sort, ONLY: t_Permutation
    USE mo_mpi, ONLY: p_comm_work, p_comm_size, p_comm_rank, p_alltoall, p_alltoallv
    USE mo_impl_constants, ONLY: SUCCESS

    IMPLICIT NONE

    PUBLIC :: t_BrokerCommunicationPattern

    ! This TYPE holds all the DATA that's necessary to communicate
    ! DATA from a given domain decomposition to a simple broker
    ! decomposition AND back.  Effectively, this behaves like a
    ! permutation that crosses process boundaries.  The communication
    ! IS performed via MPI_Alltoall() AND MPI_Alltoallv().
    TYPE t_BrokerCommunicationPattern
        INTEGER :: localPoints, brokerPoints
        TYPE(t_Permutation) :: localToBufferPerm, bufferToBrokerPerm
        INTEGER, ALLOCATABLE :: toBrokerCounts(:), fromProcessCounts(:)
        INTEGER, ALLOCATABLE :: toBrokerDisplacements(:), fromProcessDisplacements(:)
        INTEGER, ALLOCATABLE :: brokerBuffer(:), localBuffer(:)
    CONTAINS
        PROCEDURE :: construct => brokerCommunicationPattern_construct
        PROCEDURE :: communicateToBroker => brokerCommunicationPattern_communicateToBroker
        PROCEDURE :: communicateFromBroker => brokerCommunicationPattern_communicateFromBroker
        PROCEDURE :: destruct => brokerCommunicationPattern_destruct
    END TYPE t_BrokerCommunicationPattern

    CHARACTER(*), PARAMETER :: modname = "mo_broker_communication"

CONTAINS

    ! Determine how many DATA points belong into each bin.
    ! Bin i IS defined to contain all DATA points IN the interval binBounds(i) < x <= binBounds(i + 1).
    ! Thus, binBounds must contain exactly one more entry than counts.
    ! All DATA points must belong into one of the bins.
    SUBROUTINE countBinSizes(dataArray, binBounds, counts)
        INTEGER, INTENT(IN) :: dataArray(:), binBounds(:)
        INTEGER, INTENT(OUT) :: counts(:)

        INTEGER :: i, lower, middle, upper
        CHARACTER(*), PARAMETER :: routine = modname//":countBinSizes"

        ! check preconditions
        IF(SIZE(counts) + 1 /= SIZE(binBounds)) CALL finish(routine, "assertion failed: invalid arguments, array size mismatch")
        IF(ANY(dataArray <= binBounds(1) .OR. dataArray > binBounds(SIZE(binBounds)))) THEN
            CALL finish(routine, "assertion failed: OUT of bounds DATA detected")
        END IF

        counts(:) = 0
        DO i = 1, SIZE(dataArray)
            lower = 1
            upper = SIZE(binBounds)
            DO  ! binary search for the correct bin
                IF(upper == lower + 1) EXIT
                middle = (lower + upper)/2
                IF(dataArray(i) > binBounds(middle)) THEN
                    lower = middle
                ELSE
                    upper = middle
                END IF
            END DO
            counts(lower) = counts(lower) + 1
        END DO
    END SUBROUTINE countBinSizes

    ! Determine the displacements IN a buffer that's filled with consecutive chunks according to the given counts.
    ! displacements(1) will be set to zero, its SIZE IS equal to the SIZE of counts (it's intended to be used IN an MPI_Alltoallv() CALL)
    FUNCTION countsToDisplacements(counts) RESULT(displacements)
        INTEGER, ALLOCATABLE :: displacements(:)
        INTEGER, INTENT(IN) :: counts(:)

        INTEGER :: i, error
        CHARACTER(*), PARAMETER :: routine = modname//":countsToDisplacements"

        ! ALLOCATE memory
        ALLOCATE(displacements(SIZE(counts)), STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation failure")

        ! fill the array
        displacements(1) = 0
        DO i = 2, SIZE(counts)
            displacements(i) = displacements(i - 1) + counts(i - 1)
        END DO
    END FUNCTION countsToDisplacements

    ! Sets up a communication pattern so that calling
    ! `communicateToBroker(myGlobalIndices, resultBuffer)` results IN
    ! `resultBuffer == [(i, i = brokerBounds(myRank) + 1,
    ! brokerBounds(myRank + 1))]` on process `myRank`.
    !
    ! To this END, the myGlobalIndices arrays on the different
    ! processes must completely AND uniquely fill the range defined by
    ! the brokerBounds array.  brokerBounds must have one more entries
    ! than there are processes, AND must be identical on all
    ! processes.
    SUBROUTINE brokerCommunicationPattern_construct(me, myGlobalIndices, brokerBounds)
        CLASS(t_BrokerCommunicationPattern), INTENT(INOUT) :: me
        INTEGER, INTENT(IN) :: myGlobalIndices(:), brokerBounds(:)

        INTEGER :: procCount, myRank, error
        CHARACTER(*), PARAMETER :: routine = modname//":brokerCommunicationPattern_construct"

        procCount = p_comm_size(p_comm_work)
        myRank = p_comm_rank(p_comm_work)
        me%localPoints = SIZE(myGlobalIndices)
        me%brokerPoints = brokerBounds(myRank + 2) - brokerBounds(myRank + 1)

        ! check the precondition
        IF(SIZE(brokerBounds) /= procCount + 1) CALL finish(routine, "assertion failed: wrong size of argument array")

        ! ALLOCATE the memory that we don't ALLOCATE implicitly later
        ALLOCATE(me%toBrokerCounts(procCount), STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation failure")
        ALLOCATE(me%fromProcessCounts(procCount), STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation failure")
        ALLOCATE(me%localBuffer(me%localPoints), STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation failure")

        ! providers compute a permutation to sort their points into
        ! the buffer ready to be sent to the brokers AND fill the
        ! buffer with their provided points
        CALL me%localToBufferPerm%construct(myGlobalIndices)

        ! providers determine how many points they have for each PE IN the broker decomposition
        me%toBrokerCounts(:) = 0
        CALL countBinSizes(myGlobalIndices, brokerBounds, me%toBrokerCounts)

        ! providers send their DATA counts to the respective brokers
        CALL p_alltoall(me%toBrokerCounts, me%fromProcessCounts, p_comm_work)

        ! ALLOCATE the broker buffer
        ! this may be larger than the range expected from our broker
        ! range due to the possibility that processes request the same
        ! points several times (halo points!)
        ALLOCATE(me%brokerBuffer(SUM(me%fromProcessCounts)), STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation failure")

        ! compute the displacements for the MPI_Alltoallv() CALL
        me%toBrokerDisplacements = countsToDisplacements(me%toBrokerCounts)
        me%fromProcessDisplacements = countsToDisplacements(me%fromProcessCounts)

        ! providers send their provided indices to the respective brokers
        CALL me%localToBufferPerm%permute(myGlobalIndices, me%localBuffer)
        CALL p_alltoallv(me%localBuffer, me%toBrokerCounts, me%toBrokerDisplacements, &
                        &me%brokerBuffer, me%fromProcessCounts, me%fromProcessDisplacements, p_comm_work)

        ! brokers compute a permutation that sorts the points from the providers
        CALL me%bufferToBrokerPerm%construct(me%brokerBuffer)
    END SUBROUTINE brokerCommunicationPattern_construct

    SUBROUTINE brokerCommunicationPattern_communicateToBroker(me, input, output)
        CLASS(t_BrokerCommunicationPattern), INTENT(INOUT) :: me
        INTEGER, INTENT(IN) :: input(:)
        INTEGER, INTENT(OUT) :: output(:)

        CALL me%localToBufferPerm%permute(input, me%localBuffer)
        CALL p_alltoallv(me%localBuffer, me%toBrokerCounts, me%toBrokerDisplacements, &
                        &me%brokerBuffer, me%fromProcessCounts, me%fromProcessDisplacements, p_comm_work)
        CALL me%bufferToBrokerPerm%permute(me%brokerBuffer, output)
    END SUBROUTINE brokerCommunicationPattern_communicateToBroker

    SUBROUTINE brokerCommunicationPattern_communicateFromBroker(me, input, output)
        CLASS(t_BrokerCommunicationPattern), INTENT(INOUT) :: me
        INTEGER, INTENT(IN) :: input(:)
        INTEGER, INTENT(OUT) :: output(:)

        CALL me%bufferToBrokerPerm%reverse(input, me%brokerBuffer)
        CALL p_alltoallv(me%brokerBuffer, me%fromProcessCounts, me%fromProcessDisplacements, &
                        &me%localBuffer, me%toBrokerCounts, me%toBrokerDisplacements, p_comm_work)
        CALL me%localToBufferPerm%reverse(me%localBuffer, output)
    END SUBROUTINE brokerCommunicationPattern_communicateFromBroker

    SUBROUTINE brokerCommunicationPattern_destruct(me)
        CLASS(t_BrokerCommunicationPattern), INTENT(INOUT) :: me

        CALL me%localToBufferPerm%destruct()
        CALL me%bufferToBrokerPerm%destruct()
        DEALLOCATE(me%toBrokerCounts, me%fromProcessCounts, me%toBrokerDisplacements, me%fromProcessDisplacements, &
                  &me%brokerBuffer, me%localBuffer)
    END SUBROUTINE brokerCommunicationPattern_destruct

END MODULE mo_broker_communication
