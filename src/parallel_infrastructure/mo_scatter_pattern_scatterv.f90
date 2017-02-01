!>
!! An implementation of t_scatterPattern that uses MPI_Scatterv to distribute the data.
!!
!! @author N. Hübbe, DWD
!!
!!
!! @par Revision History
!! Initial hack: 2014-08-18 : N. Hübbe, DWD
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.

MODULE mo_scatter_pattern_scatterv
    USE mo_impl_constants, ONLY: SUCCESS
    USE mo_scatter_pattern_base
    USE mo_kind, ONLY: wp, dp, sp, i8
#ifndef NOMPI
    USE mpi
#endif
    USE mo_mpi, ONLY: p_io, my_process_is_stdio, &
    &                 p_real_dp, p_real_sp, p_int, &
    &                 p_comm_size, p_gather, p_gatherv
    USE mo_parallel_config, ONLY: blk_no, idx_no
    USE mo_exception, ONLY: finish

    IMPLICIT NONE

PUBLIC :: t_scatterPatternScatterV

    TYPE, EXTENDS(t_scatterPattern) :: t_scatterPatternScatterV
        !This global description is only created on p_io.
        INTEGER, POINTER :: pointCounts(:), displacements(:)    !The recvcounts/sendcounts & displs arrays for
                                                                !MPI_GATHERV/MPI_SCATTERV. Indexes into pointIndices or other
                                                                !arrays that need to be distributed.
        INTEGER, POINTER :: pointIndices(:)    !For each point requested by a process, this lists the global index of the point.
        INTEGER :: pointCount !size of pointIndices
    CONTAINS
        PROCEDURE :: construct       => constructScatterPatternScatterV !< override
        PROCEDURE :: distribute_dp   => distributeDataScatterV_dp       !< override
        PROCEDURE :: distribute_spdp => distributeDataScatterV_spdp     !< override
        PROCEDURE :: distribute_sp   => distributeDataScatterV_sp       !< override
        PROCEDURE :: distribute_int  => distributeDataScatterV_int      !< override
        PROCEDURE :: destruct        => destructScatterPatternScatterV  !< override
    END TYPE

PRIVATE

    CHARACTER(*), PARAMETER :: modname = "mo_grid_distribution_scatterv"
    LOGICAL, PARAMETER :: debugModule = .false.

CONTAINS

    !-------------------------------------------------------------------------------------------------------------------------------
    !> constructor
    !-------------------------------------------------------------------------------------------------------------------------------
    SUBROUTINE constructScatterPatternScatterV(me, jg, loc_arr_len, glb_index, communicator)
        CLASS(t_scatterPatternScatterV), TARGET, INTENT(OUT) :: me
        INTEGER, VALUE :: jg, loc_arr_len, communicator
        INTEGER, INTENT(IN) :: glb_index(:)

        CHARACTER(*), PARAMETER :: routine = modname//":constructScatterPatternScatterV"
        INTEGER :: procCount, i, ierr

        IF(debugModule .and. my_process_is_stdio()) WRITE(0,*) "entering ", routine

        CALL constructScatterPattern(me, jg, loc_arr_len, glb_index, communicator)
        IF(my_process_is_stdio()) THEN
            procCount = p_comm_size(communicator)
            ALLOCATE(me%pointCounts(procCount), stat = ierr)
            IF(ierr /= SUCCESS) CALL finish(routine, "error allocating memory")
            ALLOCATE(me%displacements(procCount), stat = ierr)
            IF(ierr /= SUCCESS) CALL finish(routine, "error allocating memory")
        ELSE
            ALLOCATE(me%pointCounts(1), stat = ierr)
            IF(ierr /= SUCCESS) CALL finish(routine, "error allocating memory")
            ALLOCATE(me%displacements(1), stat = ierr)
            IF(ierr /= SUCCESS) CALL finish(routine, "error allocating memory")
        END IF
        CALL p_gather(me%myPointCount, me%pointCounts, p_io, communicator)
        IF(my_process_is_stdio()) THEN
            me%pointCount = 0
            DO i = 1, procCount
                me%displacements(i) = me%pointCount
                me%pointCount = me%pointCount + me%pointCounts(i)
            END DO
            ALLOCATE(me%pointIndices(me%pointCount), stat = ierr)
            IF(ierr /= SUCCESS) CALL finish(routine, "error allocating memory")
        ELSE
            ALLOCATE(me%pointIndices(1), stat = ierr)
            IF(ierr /= SUCCESS) CALL finish(routine, "error allocating memory")
        END IF
        CALL p_gatherv(glb_index, loc_arr_len, me%pointIndices, me%pointCounts, me%displacements, p_io, communicator)

        IF(debugModule .and. my_process_is_stdio()) WRITE(0,*) "leaving ", routine
    END SUBROUTINE constructScatterPatternScatterV

    !-------------------------------------------------------------------------------------------------------------------------------
    !> implementation of t_scatterPattern::distribute_dp()
    !-------------------------------------------------------------------------------------------------------------------------------
    SUBROUTINE distributeDataScatterV_dp(me, globalArray, localArray, ladd_value)
        IMPLICIT NONE
        CLASS(t_scatterPatternScatterV), INTENT(INOUT) :: me
        REAL(dp), INTENT(INOUT) :: globalArray(:)
        REAL(wp), INTENT(INOUT) :: localArray(:,:)
        LOGICAL, INTENT(IN) :: ladd_value

        CHARACTER(*), PARAMETER :: routine = modname//":distributeDataScatterV_dp"
        REAL(dp), ALLOCATABLE :: sendArray(:), recvArray(:)
        INTEGER :: i, blk, idx, ierr

        IF(debugModule .and. my_process_is_stdio()) WRITE(0,*) "entering ", routine
        CALL me%startDistribution()

        IF(my_process_is_stdio()) THEN
            ALLOCATE(sendArray(me%pointCount), stat = ierr)
            IF(ierr /= SUCCESS) CALL finish(routine, "error allocating memory")
            DO i = 1, me%pointCount
                sendArray(i) = globalArray(me%pointIndices(i))
            END DO
        ELSE
            ALLOCATE(sendArray(1), stat = ierr)
            IF(ierr /= SUCCESS) CALL finish(routine, "error allocating memory")
        END IF
        ALLOCATE(recvArray(me%myPointCount), stat = ierr)
        IF(ierr /= SUCCESS) CALL finish(routine, "error allocating memory")

        !For some very weird reason, I always get a crash on my system if the following block is moved to its own subroutine.
!       CALL p_scatterv(sendArray, me%pointCounts, me%displacements, recvArray, me%myPointCount, p_io, me%communicator)
#ifndef NOMPI
        CALL MPI_Scatterv(sendArray, me%pointCounts, me%displacements, p_real_dp, recvArray, me%myPointCount, p_real_dp, p_io, &
        &                 me%communicator, ierr)
        IF (ierr /=  MPI_SUCCESS) CALL finish (routine, 'Error in MPI_Scatterv operation!')
#else
        recvArray(1:me%myPointCount) = sendArray((me%displacements(1)+1):(me%displacements(1)+me%myPointCount))
#endif

        IF(ladd_value) THEN
            DO i = 1, me%myPointCount
                blk = blk_no(i)
                idx = idx_no(i)
                localArray(idx, blk) = localArray(idx, blk) + recvArray(i)
            END DO
        ELSE
            DO i = 1, me%myPointCount
                localArray(idx_no(i), blk_no(i)) = recvArray(i)
            END DO
        END IF

        DEALLOCATE(recvArray)
        DEALLOCATE(sendArray)
        CALL me%endDistribution(INT(me%pointCount, i8) * 8_i8)
        IF(debugModule .and. my_process_is_stdio()) WRITE(0,*) "leaving ", routine
    END SUBROUTINE distributeDataScatterV_dp

    !-------------------------------------------------------------------------------------------------------------------------------
    !> implementation of t_scatterPattern::distribute_spdp()
    !-------------------------------------------------------------------------------------------------------------------------------
    SUBROUTINE distributeDataScatterV_spdp(me, globalArray, localArray, ladd_value)
        IMPLICIT NONE
        CLASS(t_scatterPatternScatterV), INTENT(INOUT) :: me
        REAL(sp), INTENT(INOUT) :: globalArray(:)
        REAL(wp), INTENT(INOUT) :: localArray(:,:)
        LOGICAL, INTENT(IN) :: ladd_value

        CHARACTER(*), PARAMETER :: routine = modname//":distributeDataScatterV_spdp"
        REAL(sp), ALLOCATABLE :: sendArray(:), recvArray(:)
        INTEGER :: i, blk, idx, ierr

        IF(debugModule .and. my_process_is_stdio()) WRITE(0,*) "entering ", routine
        CALL me%startDistribution()

        IF(my_process_is_stdio()) THEN
            ALLOCATE(sendArray(me%pointCount), stat = ierr)
            IF(ierr /= SUCCESS) CALL finish(routine, "error allocating memory")
            DO i = 1, me%pointCount
                sendArray(i) = globalArray(me%pointIndices(i))
            END DO
        ELSE
            ALLOCATE(sendArray(1), stat = ierr)
            IF(ierr /= SUCCESS) CALL finish(routine, "error allocating memory")
        END IF
        ALLOCATE(recvArray(me%myPointCount), stat = ierr)
        IF(ierr /= SUCCESS) CALL finish(routine, "error allocating memory")

        !For some very weird reason, I always get a crash on my system if the following block is moved to its own subroutine.
!       CALL p_scatterv(sendArray, me%pointCounts, me%displacements, recvArray, me%myPointCount, p_io, me%communicator)
#ifndef NOMPI
        CALL MPI_Scatterv(sendArray, me%pointCounts, me%displacements, p_real_sp, recvArray, me%myPointCount, p_real_sp, p_io, &
        &                 me%communicator, ierr)
        IF (ierr /=  MPI_SUCCESS) CALL finish (routine, 'Error in MPI_Scatterv operation!')
#else
        recvArray(1:me%myPointCount) = sendArray((me%displacements(1)+1):(me%displacements(1)+me%myPointCount))
#endif

        IF(ladd_value) THEN
            DO i = 1, me%myPointCount
                blk = blk_no(i)
                idx = idx_no(i)
                localArray(idx, blk) = localArray(idx, blk) + recvArray(i)
            END DO
        ELSE
            DO i = 1, me%myPointCount
                localArray(idx_no(i), blk_no(i)) = recvArray(i)
            END DO
        END IF

        DEALLOCATE(recvArray)
        DEALLOCATE(sendArray)
        CALL me%endDistribution(INT(me%pointCount, i8) * 4_i8)
        IF(debugModule .and. my_process_is_stdio()) WRITE(0,*) "leaving ", routine
    END SUBROUTINE distributeDataScatterV_spdp


    !-------------------------------------------------------------------------------------------------------------------------------
    !> implementation of t_scatterPattern::distribute_sp()
    !-------------------------------------------------------------------------------------------------------------------------------
    SUBROUTINE distributeDataScatterV_sp(me, globalArray, localArray, ladd_value)
        IMPLICIT NONE
        CLASS(t_scatterPatternScatterV), INTENT(INOUT) :: me
        REAL(sp), INTENT(INOUT) :: globalArray(:)
        REAL(sp), INTENT(INOUT) :: localArray(:,:)
        LOGICAL, INTENT(IN) :: ladd_value

        CHARACTER(*), PARAMETER :: routine = modname//":distributeDataScatterV_sp"
        REAL(sp), ALLOCATABLE :: sendArray(:), recvArray(:)
        INTEGER :: i, blk, idx, ierr

        IF(debugModule .and. my_process_is_stdio()) WRITE(0,*) "entering ", routine
        CALL me%startDistribution()

        IF(my_process_is_stdio()) THEN
            ALLOCATE(sendArray(me%pointCount), stat = ierr)
            IF(ierr /= SUCCESS) CALL finish(routine, "error allocating memory")
            DO i = 1, me%pointCount
                sendArray(i) = globalArray(me%pointIndices(i))
            END DO
        ELSE
            ALLOCATE(sendArray(1), stat = ierr)
            IF(ierr /= SUCCESS) CALL finish(routine, "error allocating memory")
        END IF
        ALLOCATE(recvArray(me%myPointCount), stat = ierr)
        IF(ierr /= SUCCESS) CALL finish(routine, "error allocating memory")

        !For some very weird reason, I always get a crash on my system if the following block is moved to its own subroutine.
!       CALL p_scatterv(sendArray, me%pointCounts, me%displacements, recvArray, me%myPointCount, p_io, me%communicator)
#ifndef NOMPI
        CALL MPI_Scatterv(sendArray, me%pointCounts, me%displacements, p_real_sp, recvArray, me%myPointCount, p_real_sp, p_io, &
        &                 me%communicator, ierr)
        IF (ierr /=  MPI_SUCCESS) CALL finish (routine, 'Error in MPI_Scatterv operation!')
#else
        recvArray(1:me%myPointCount) = sendArray((me%displacements(1)+1):(me%displacements(1)+me%myPointCount))
#endif

        IF(ladd_value) THEN
            DO i = 1, me%myPointCount
                blk = blk_no(i)
                idx = idx_no(i)
                localArray(idx, blk) = localArray(idx, blk) + recvArray(i)
            END DO
        ELSE
            DO i = 1, me%myPointCount
                localArray(idx_no(i), blk_no(i)) = recvArray(i)
            END DO
        END IF

        DEALLOCATE(recvArray)
        DEALLOCATE(sendArray)
        CALL me%endDistribution(INT(me%pointCount, i8) * 4_i8)
        IF(debugModule .and. my_process_is_stdio()) WRITE(0,*) "leaving ", routine
    END SUBROUTINE distributeDataScatterV_sp

    !-------------------------------------------------------------------------------------------------------------------------------
    !> implementation of t_scatterPattern::distribute_int()
    !-------------------------------------------------------------------------------------------------------------------------------
    SUBROUTINE distributeDataScatterV_int(me, globalArray, localArray, ladd_value)
        IMPLICIT NONE
        CLASS(t_scatterPatternScatterV), INTENT(INOUT) :: me
        INTEGER, INTENT(INOUT) :: globalArray(:)
        INTEGER, INTENT(INOUT) :: localArray(:,:)
        LOGICAL, INTENT(IN) :: ladd_value

        CHARACTER(*), PARAMETER :: routine = modname//":distributeDataScatterV_sp"
        INTEGER, ALLOCATABLE :: sendArray(:), recvArray(:)
        INTEGER :: i, blk, idx, ierr

        IF(debugModule .and. my_process_is_stdio()) WRITE(0,*) "entering ", routine
        CALL me%startDistribution()

        IF(my_process_is_stdio()) THEN
            ALLOCATE(sendArray(me%pointCount), stat = ierr)
            IF(ierr /= SUCCESS) CALL finish(routine, "error allocating memory")
            DO i = 1, me%pointCount
                sendArray(i) = globalArray(me%pointIndices(i))
            END DO
        ELSE
            ALLOCATE(sendArray(1), stat = ierr)
            IF(ierr /= SUCCESS) CALL finish(routine, "error allocating memory")
        END IF
        ALLOCATE(recvArray(me%myPointCount), stat = ierr)
        IF(ierr /= SUCCESS) CALL finish(routine, "error allocating memory")

        !For some very weird reason, I always get a crash on my system if the following block is moved to its own subroutine.
!       CALL p_scatterv(sendArray, me%pointCounts, me%displacements, recvArray, me%myPointCount, p_io, me%communicator)
#ifndef NOMPI
        CALL MPI_Scatterv(sendArray, me%pointCounts, me%displacements, p_int, recvArray, me%myPointCount, p_int, p_io, &
        &                 me%communicator, ierr)
        IF (ierr /=  MPI_SUCCESS) CALL finish (routine, 'Error in MPI_Scatterv operation!')
#else
        recvArray(1:me%myPointCount) = sendArray((me%displacements(1)+1):(me%displacements(1)+me%myPointCount))
#endif

        IF(ladd_value) THEN
            DO i = 1, me%myPointCount
                blk = blk_no(i)
                idx = idx_no(i)
                localArray(idx, blk) = localArray(idx, blk) + recvArray(i)
            END DO
        ELSE
            DO i = 1, me%myPointCount
                localArray(idx_no(i), blk_no(i)) = recvArray(i)
            END DO
        END IF

        DEALLOCATE(recvArray)
        DEALLOCATE(sendArray)
        CALL me%endDistribution(INT(me%pointCount, i8) * 4_i8)
        IF(debugModule .and. my_process_is_stdio()) WRITE(0,*) "leaving ", routine
    END SUBROUTINE distributeDataScatterV_int

    !-------------------------------------------------------------------------------------------------------------------------------
    !> destructor
    !-------------------------------------------------------------------------------------------------------------------------------
    SUBROUTINE destructScatterPatternScatterV(me)
        CLASS(t_scatterPatternScatterV), TARGET, INTENT(INOUT) :: me

        CHARACTER(*), PARAMETER :: routine = modname//":destructScatterPatternScatterV"

        IF(debugModule .and. my_process_is_stdio()) WRITE(0,*) "entering ", routine
        DEALLOCATE(me%pointCounts)
        DEALLOCATE(me%displacements)
        DEALLOCATE(me%pointIndices)
        CALL destructScatterPattern(me)
        IF(debugModule .and. my_process_is_stdio()) WRITE(0,*) "leaving ", routine
    END SUBROUTINE destructScatterPatternScatterV

END MODULE
