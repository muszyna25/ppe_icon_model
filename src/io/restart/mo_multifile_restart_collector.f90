!> module for collecting the payload DATA onto the respective restart processes
!!
!! Initial implementation: Nathanael Hübbe
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.

MODULE mo_multifile_restart_collector
    USE mo_communication, ONLY: idx_no, blk_no
    USE mo_decomposition_tools, ONLY: t_grid_domain_decomp_info
    USE mo_exception, ONLY: finish
    USE mo_fortran_tools, ONLY: alloc, t_Destructible
    USE mo_impl_constants, ONLY: SUCCESS
    USE mo_kind, ONLY: dp, sp
    USE mo_mpi, ONLY: p_comm_work_restart, p_comm_rank, p_send, p_recv, my_process_is_work
    USE mo_util_string, ONLY: int2string

    IMPLICIT NONE

    PUBLIC :: t_MultifileRestartCollector, t_MultifileRestartCollector_ptr

    PRIVATE

    !This CLASS IS used to gather the restart payload DATA onto the restart writing processes.
    !
    !
    !Usage notes:
    !
    ! 1. The constructor IS a FUNCTION. It returns a field with the
    ! global indices of the points written by this PE.
    !
    ! 2. The roles of the different processes IS determined by the
    !    constructor arguments ONLY.  There IS no logic coded into
    !    this CLASS that requires a specific placement, selecting a
    !    good placement IS entirely up to the caller.
    !
    ! 3. The input of the collection IS a 2D slice (blk,idx), however, the output IS a 1D array of points.
    !
    !
    !Implementation notes:
    !
    !  * The collection IS split into two parts:
    !     1. Collection of the DATA to be sent into a 1D buffer.
    !        This IS a perfectly local operation.
    !     2. Collection of the buffer contents IN the RESULT 1D arrays on the writer PEs.
    !        This IS the communication part which IS implemented IN collectBuffer().
    !    The reason for this split IS that it allows the second part
    !    to be used already IN the constructor to produce the
    !    resulting array of global point indices.
    TYPE, EXTENDS(t_Destructible) :: t_MultifileRestartCollector
        PRIVATE
        !All ALLOCATABLE arrays are ALLOCATED on both work AND restart
        !processes, however, they may be empty on either dedicated
        !restart OR pure worker processes.
        INTEGER :: receivePointCount    !number of points received by this PE, zero on pure work PEs

        !SIZE of sourceProcs AND sourcePointCounts, number of PEs
        !sending DATA to this one, zero on pure work PEs:
        INTEGER :: sourceProcCount

        !the ranks of the source processes IN p_comm_work_restart,
        !empty on pure work procs:
        INTEGER, ALLOCATABLE :: sourceProcs(:)

        !the count of points sent by each source PE, empty on pure
        !work procs:
        INTEGER, ALLOCATABLE :: sourcePointCounts(:)    

        !SIZE of sendIdx, sendBlk, AND sendBuffer_x; zero on dedicated
        !restart procs:
        INTEGER :: sendPointCount

        !This defines how the points of this PE are collected into the
        !array that IS sent to the restart proc. ALLOCATED on both
        !restart AND work processes, but on dedicated restart procs,
        !there are no entries:
        INTEGER, ALLOCATABLE :: sendIdx(:), sendBlk(:)

        INTEGER :: destProc !the rank within p_comm_work_restart of the process writing our DATA
        REAL(dp), ALLOCATABLE :: sendBuffer_d(:)    !buffer for linerarizing the DATA before sending it to the writer process
        REAL(sp), ALLOCATABLE :: sendBuffer_s(:)    !buffer for linerarizing the DATA before sending it to the writer process
        REAL(sp), ALLOCATABLE :: sendBuffer_int(:)  !buffer for linerarizing the DATA before sending it to the writer process
    CONTAINS

        !collective across work AND restart (actually ONLY among some
        !subsets indicated by the arguments, but since the SUM of
        !these subsets IS the total set, there IS no REAL difference)
        PROCEDURE :: construct => multifileRestartCollector_construct

        !same collectivity as construct()
        PROCEDURE :: collectField_d   => multifileRestartCollector_collectField_d
        PROCEDURE :: collectField_s   => multifileRestartCollector_collectField_s
        PROCEDURE :: collectField_int => multifileRestartCollector_collectField_int
        GENERIC :: collectField => collectField_d, collectField_s, collectField_int

        PROCEDURE :: destruct => multifileRestartCollector_destruct    !override

        !same collectivity as construct()
        PROCEDURE, PRIVATE :: collectBuffer_d   => multifileRestartCollector_collectBuffer_d
        PROCEDURE, PRIVATE :: collectBuffer_s   => multifileRestartCollector_collectBuffer_s
        PROCEDURE, PRIVATE :: collectBuffer_int => multifileRestartCollector_collectBuffer_int
        GENERIC, PRIVATE :: collectBuffer => collectBuffer_d, collectBuffer_s, collectBuffer_int
    END TYPE t_MultifileRestartCollector

    TYPE t_MultifileRestartCollector_ptr
        TYPE(t_MultifileRestartCollector), POINTER :: p
    END TYPE t_MultifileRestartCollector_ptr

    CHARACTER(*), PARAMETER :: modname = "mo_multifile_restart_collector"

CONTAINS

    !On the restart writers, this returns an array with the global
    !indices of the points collected by this PE. The RETURN VALUE IS
    !NOT ALLOCATED on pure worker procs.
    FUNCTION multifileRestartCollector_construct(me, localPointCount, decompInfo, destProc, sourceProcs) RESULT(globalIndices)
        REAL(dp), POINTER :: globalIndices(:)   !should have been of TYPE INTEGER, but CDI does NOT support writing of integers
        CLASS(t_MultifileRestartCollector), INTENT(INOUT) :: me
        INTEGER, VALUE :: localPointCount   !the number of points PRESENT on this PE
        TYPE(t_grid_domain_decomp_info), INTENT(IN) :: decompInfo

        !rank within p_comm_work_restart to which this process should
        !send its DATA, must be the own rank on the restart writer
        !procs
        INTEGER, VALUE :: destProc

        !ranks within p_comm_work_restart from which DATA IS to be
        !collected, must include the own rank on the restart writer
        !procs; empty on pure work PEs
        INTEGER, INTENT(IN) :: sourceProcs(:)

        INTEGER :: myRank, error, i, j
        CHARACTER(*), PARAMETER :: routine = modname//":multifileRestartCollector_construct"

        myRank = p_comm_rank(p_comm_work_restart)

        !Copy the information about the processes IN our group, AND
        !ALLOCATE the sourceProcCount dependent arrays.
        me%sourceProcCount = SIZE(sourceProcs)
        CALL alloc(me%sourceProcs, me%sourceProcCount)
        CALL alloc(me%sourcePointCounts, me%sourceProcCount)
        me%sourceProcs(1:SIZE(sourceProcs)) = sourceProcs(:)
        me%destProc = destProc

        !Compute the number of points we want to send, AND ALLOCATE
        !the sendPointCount dependent arrays.
        me%sendPointCount = 0
        IF(my_process_is_work()) THEN
            DO i = 1, localPointCount
                IF(decompInfo%owner_mask(idx_no(i), blk_no(i))) me%sendPointCount = me%sendPointCount + 1
            END DO
        END IF
        CALL alloc(me%sendIdx, me%sendPointCount)
        CALL alloc(me%sendBlk, me%sendPointCount)
        CALL alloc(me%sendBuffer_d,   me%sendPointCount)
        CALL alloc(me%sendBuffer_s,   me%sendPointCount)
        CALL alloc(me%sendBuffer_int, me%sendPointCount)

        !Fill sendIdx, sendBlk, AND sendBuffer_d.
        IF(my_process_is_work()) THEN
            j = 1
            DO i = 1, localPointCount
                IF(decompInfo%owner_mask(idx_no(i), blk_no(i))) THEN
                    me%sendIdx(j) = idx_no(i)
                    me%sendBlk(j) = blk_no(i)
                    me%sendBuffer_d(j) = decompInfo%glb_index(i)
                    j = j + 1
                END IF
            END DO
        END IF

        !Collect the source point counts.
        DO i = 1, me%sourceProcCount
            IF(me%sourceProcs(i) == myRank) THEN
                me%sourcePointCounts(i) = me%sendPointCount
            ELSE
                CALL p_recv(me%sourcePointCounts(i), me%sourceProcs(i), 0, comm = p_comm_work_restart)
            END IF
        END DO
        IF(me%destProc /= myRank) CALL p_send(me%sendPointCount, me%destProc, 0, comm = p_comm_work_restart)

        !Compute the receivePointCount AND ALLOCATE the globalIndices array.
        me%receivePointCount = SUM(me%sourcePointCounts(1:me%sourceProcCount))
        ALLOCATE(globalIndices(me%receivePointCount), STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation failure")

        !Collect the global indices of the points written by this PE.
        CALL me%collectBuffer(globalIndices)
    END FUNCTION multifileRestartCollector_construct

    !Collect the DATA within the send buffers to the respective outputData arrays.
    SUBROUTINE multifileRestartCollector_collectBuffer_d(me, outputData)
        CLASS(t_MultifileRestartCollector), INTENT(INOUT) :: me
        REAL(dp), INTENT(OUT) :: outputData(:)

        INTEGER :: myRank, i, j
        CHARACTER(*), PARAMETER :: routine = modname//":multifileRestartCollector_collectBuffer_d"
        REAL(dp) :: dummy(1)

        !sanity check
        myRank = p_comm_rank(p_comm_work_restart)
        IF(myRank == me%destProc .AND. SIZE(outputData) /= me%receivePointCount) THEN
            CALL finish(routine, "assertion failed: wrong buffer size (expected "//TRIM(int2string(me%receivePointCount))//&
                               & ", got "//TRIM(int2string(SIZE(outputData)))//"), check the calling routine")
        END IF

        !Collect the DATA on the writer PEs.
        j = 1
        outputData(:) = 0._dp
        DO i = 1, me%sourceProcCount
            IF(me%sourceProcs(i) == myRank) THEN
                outputData(j:j-1+me%sourcePointCounts(i)) = me%sendBuffer_d(1:me%sourcePointCounts(i))
            ELSE
              IF (me%sourcePointCounts(i) > 0) THEN
                CALL p_recv(outputData(j:j-1+me%sourcePointCounts(i)), me%sourceProcs(i), 0, comm = p_comm_work_restart)
              ELSE
                CALL p_recv(dummy(1), me%sourceProcs(i), 0, comm = p_comm_work_restart)
              END IF
            END IF
            j = j + me%sourcePointCounts(i)
        END DO
        IF(me%destProc /= myRank) THEN
          IF (me%sendPointCount > 0) THEN
            CALL p_send(me%sendBuffer_d(1:me%sendPointCount), me%destProc, 0, comm = p_comm_work_restart)
          ELSE
            CALL p_send(dummy(1), me%destProc, 0, comm = p_comm_work_restart)
          END IF
        END IF

        !sanity check
        IF(j /= me%receivePointCount + 1) CALL finish(routine, "assertion failed")
    END SUBROUTINE multifileRestartCollector_collectBuffer_d

    !Collect the DATA within the send buffers to the respective outputData arrays.
    SUBROUTINE multifileRestartCollector_collectBuffer_s(me, outputData)
        CLASS(t_MultifileRestartCollector), INTENT(INOUT) :: me
        REAL(sp), INTENT(OUT) :: outputData(:)

        INTEGER :: myRank, i, j
        CHARACTER(*), PARAMETER :: routine = modname//":multifileRestartCollector_collectBuffer_s"
        REAL(sp) :: dummy(1)

        !sanity check
        myRank = p_comm_rank(p_comm_work_restart)
        IF(myRank == me%destProc .AND. SIZE(outputData) /= me%receivePointCount) THEN
            CALL finish(routine, "assertion failed, SIZE(outputData) = "//TRIM(int2string(SIZE(outputData)))// &
                               & ", me%receivePointCount = "//TRIM(int2string(me%receivePointCount)))
        END IF

        !Collect the DATA on the writer PEs.
        j = 1
        outputData(:) = 0._sp
        DO i = 1, me%sourceProcCount
            IF(me%sourceProcs(i) == myRank) THEN
                outputData(j:j-1+me%sourcePointCounts(i)) = me%sendBuffer_s(1:me%sourcePointCounts(i))
            ELSE
              IF (me%sourcePointCounts(i) > 0) THEN
                CALL p_recv(outputData(j:j-1+me%sourcePointCounts(i)), me%sourceProcs(i), 0, comm = p_comm_work_restart)
              ELSE
                CALL p_recv(dummy(1), me%sourceProcs(i), 0, comm = p_comm_work_restart)
              END IF
            END IF
            j = j + me%sourcePointCounts(i)
        END DO
        IF(me%destProc /= myRank) THEN
          IF (me%sendPointCount > 0) THEN
            CALL p_send(me%sendBuffer_s(1:me%sendPointCount), me%destProc, 0, comm = p_comm_work_restart)
          ELSE
            CALL p_send(dummy(1), me%destProc, 0, comm = p_comm_work_restart)
          END IF
        END IF

        !sanity check
        IF(j /= me%receivePointCount + 1) CALL finish(routine, "assertion failed")
    END SUBROUTINE multifileRestartCollector_collectBuffer_s

    !Collect the DATA within the send buffers to the respective outputData arrays.
    SUBROUTINE multifileRestartCollector_collectBuffer_int(me, outputData)
        CLASS(t_MultifileRestartCollector), INTENT(INOUT) :: me
        INTEGER, INTENT(OUT) :: outputData(:)

        INTEGER :: myRank, i, j
        CHARACTER(*), PARAMETER :: routine = modname//":multifileRestartCollector_collectBuffer_int"
        INTEGER :: dummy(1)

        !sanity check
        myRank = p_comm_rank(p_comm_work_restart)
        IF(myRank == me%destProc .AND. SIZE(outputData) /= me%receivePointCount) THEN
            CALL finish(routine, "assertion failed, SIZE(outputData) = "//TRIM(int2string(SIZE(outputData)))// &
                               & ", me%receivePointCount = "//TRIM(int2string(me%receivePointCount)))
        END IF

        !Collect the DATA on the writer PEs.
        j = 1
        outputData(:) = 0
        DO i = 1, me%sourceProcCount
            IF(me%sourceProcs(i) == myRank) THEN
                outputData(j:j-1+me%sourcePointCounts(i)) = me%sendBuffer_int(1:me%sourcePointCounts(i))
            ELSE
              IF (me%sourcePointCounts(i) > 0) THEN
                CALL p_recv(outputData(j:j-1+me%sourcePointCounts(i)), me%sourceProcs(i), 0, comm = p_comm_work_restart)
              ELSE
                CALL p_recv(dummy(1), me%sourceProcs(i), 0, comm = p_comm_work_restart)
              END IF
            END IF
            j = j + me%sourcePointCounts(i)
        END DO
        IF(me%destProc /= myRank) THEN
          IF (me%sendPointCount > 0) THEN
            CALL p_send(me%sendBuffer_int(1:me%sendPointCount), me%destProc, 0, comm = p_comm_work_restart)
          ELSE
            CALL p_send(dummy(1), me%destProc, 0, comm = p_comm_work_restart)
          END IF
        END IF

        !sanity check
        IF(j /= me%receivePointCount + 1) CALL finish(routine, "assertion failed")
      END SUBROUTINE multifileRestartCollector_collectBuffer_int

    SUBROUTINE multifileRestartCollector_collectField_d(me, inputData, outputData)
        CLASS(t_MultifileRestartCollector), INTENT(INOUT) :: me
        REAL(dp), INTENT(IN) :: inputData(:,:)
        REAL(dp), INTENT(OUT) :: outputData(:)

        INTEGER :: i

        !Fill the send buffer.
        DO i = 1, me%sendPointCount
            me%sendBuffer_d(i) = inputData(me%sendIdx(i), me%sendBlk(i))
        END DO

        !Actually collect the DATA.
        CALL me%collectBuffer(outputData)
    END SUBROUTINE multifileRestartCollector_collectField_d

    SUBROUTINE multifileRestartCollector_collectField_s(me, inputData, outputData)
        CLASS(t_MultifileRestartCollector), INTENT(INOUT) :: me
        REAL(sp), INTENT(IN) :: inputData(:,:)
        REAL(sp), INTENT(OUT) :: outputData(:)

        INTEGER :: i

        !Fill the send buffer.
        DO i = 1, me%sendPointCount
            me%sendBuffer_s(i) = inputData(me%sendIdx(i), me%sendBlk(i))
        END DO

        !Actually collect the DATA.
        CALL me%collectBuffer(outputData)
    END SUBROUTINE multifileRestartCollector_collectField_s

    SUBROUTINE multifileRestartCollector_collectField_int(me, inputData, outputData)
        CLASS(t_MultifileRestartCollector), INTENT(INOUT) :: me
        INTEGER, INTENT(IN)  :: inputData(:,:)
        INTEGER, INTENT(OUT) :: outputData(:)

        INTEGER :: i

        !Fill the send buffer.
        DO i = 1, me%sendPointCount
            me%sendBuffer_int(i) = inputData(me%sendIdx(i), me%sendBlk(i))
        END DO

        !Actually collect the DATA.
        CALL me%collectBuffer(outputData)
    END SUBROUTINE multifileRestartCollector_collectField_int

    SUBROUTINE multifileRestartCollector_destruct(me)
        CLASS(t_MultifileRestartCollector), INTENT(INOUT) :: me

        CHARACTER(*), PARAMETER :: routine = modname//":multifileRestartCollector_destruct"

        DEALLOCATE(me%sourceProcs, me%sourcePointCounts, me%sendIdx, me%sendBlk)
        DEALLOCATE(me%sendBuffer_d)
        DEALLOCATE(me%sendBuffer_s)
        DEALLOCATE(me%sendBuffer_int)
    END SUBROUTINE multifileRestartCollector_destruct

END MODULE mo_multifile_restart_collector
