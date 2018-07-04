!> Module for collecting the payload data onto the respective restart
!> processes.
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

  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_ptr, C_F_POINTER, C_NULL_PTR

#ifndef NOMPI
  USE mpi
#endif
  USE mo_communication,       ONLY: idx_no, blk_no
  USE mo_decomposition_tools, ONLY: t_grid_domain_decomp_info
  USE mo_exception,           ONLY: finish, message_text
  USE mo_fortran_tools,       ONLY: alloc, ensureSize, no_copy
  USE mo_impl_constants,      ONLY: SUCCESS, SINGLE_T, REAL_T, INT_T
  USE mo_kind,                ONLY: dp, sp, i8
  USE mo_mpi,                 ONLY: p_comm_work_restart, p_comm_rank, p_send, p_recv, &
    &                               my_process_is_work, p_int, p_real_dp, p_real_sp
  USE mo_restart_util,        ONLY: my_process_is_restart_writer
  USE mo_timer,               ONLY: timer_start, timer_stop, timer_restart_collector_setup, &
    &                               timer_restart_indices_setup, timers_level
  
  IMPLICIT NONE
  
  PUBLIC :: t_MultifileRestartCollector,     &
    &       t_CollectorIndices,              &
    &       t_CollectorSendBuffer
  
  PRIVATE
  
  LOGICAL, SAVE :: have_MPI_RGet = .false.
  ! buffers for linearizing the data before sending it to the writer
  ! process

  TYPE :: ptr_arr_t
    REAL(KIND=dp), POINTER :: p(:)
  END TYPE ptr_arr_t

  TYPE :: t_CollectorSendBuffer
    REAL(dp), POINTER :: sendBuffer_d(:)
    REAL(sp), POINTER :: sendBuffer_s(:)
    INTEGER,  POINTER :: sendBuffer_int(:)
    INTEGER :: win_d, win_s, win_int
    INTEGER,           PRIVATE :: winComm, &
      &                          winClGroup, winSvGroup
    INTEGER(KIND=i8),  PRIVATE :: winSizes(3)
    TYPE(ptr_arr_t),   PRIVATE :: winPtr(3)
    LOGICAL,           PRIVATE :: allocd, handshaked, &
      &                          win_posted, win_started
  CONTAINS
    PROCEDURE :: construct    => t_CollectorSendBuffer_construct
    PROCEDURE, PRIVATE :: handshake  => t_CollectorSendBuffer_handshake
    PROCEDURE :: finalize     => t_CollectorSendBuffer_finalize
    PROCEDURE :: start_local_access  => t_CollectorSendBuffer_start_local_access
    PROCEDURE :: start_remote_access => t_CollectorSendBuffer_start_remote_access
  END TYPE t_CollectorSendBuffer

  TYPE :: t_CollectorIndices
    !All allocatable arrays are allocated on both work and restart
    !processes, however, they may be empty on either dedicated
    !restart or pure worker processes.
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
    ! the rank within p_comm_work_restart of the process writing our
    ! data
    INTEGER :: destProc
  CONTAINS
    !collective across work AND restart (actually ONLY among some
    !subsets indicated by the arguments, but since the SUM of
    !these subsets IS the total set, there IS no REAL difference)
    PROCEDURE :: construct => t_CollectorIndices_construct
    PROCEDURE :: finalize  => t_CollectorIndices_finalize
  END TYPE t_CollectorIndices

  ! This class is used to gather the restart payload DATA onto the
  ! restart writing processes.
  !
  !
  ! Usage notes:
  !
  ! 1. The constructor is a function. It returns a field with the
  ! global indices of the points written by this PE.
  !
  ! 2. The roles of the different processes IS determined by the
  !    constructor arguments ONLY.  There IS no logic coded into
  !    this CLASS that requires a specific placement, selecting a
  !    good placement IS entirely up to the caller.
  !
  ! 3. The input of the collection IS a 2D slice (blk,idx), however,
  !    the output IS a 1D array of points.
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
  !
  TYPE :: t_MultifileRestartCollector
    ! index arrays corresponding to this collector
    TYPE(t_CollectorIndices),     POINTER    :: idx
    ! pointer to global send buffer / MPI window
    TYPE(t_CollectorSendBuffer),  POINTER    :: glb_sendbuf
    ! send buffers (pointing into global send buffer)
    INTEGER(KIND=i8), ALLOCATABLE :: send_buffer_offset(:)
    ! on receiving PEs: start index in send buffer for every source
    ! PE:
    INTEGER, ALLOCATABLE :: sourcePointStart(:)
  CONTAINS
    PROCEDURE :: construct => t_multifileRestartCollector_construct
    PROCEDURE :: finalize  => multifileRestartCollector_finalize   
    PROCEDURE :: sendField_d   => multifileRestartCollector_sendField_d
    PROCEDURE :: sendField_s   => multifileRestartCollector_sendField_s
    PROCEDURE :: sendField_int => multifileRestartCollector_sendField_int
    PROCEDURE, PRIVATE :: sendField_generic => multifileRestartCollector_sendField_generic
    GENERIC :: sendField => sendField_d, sendField_s, sendField_int
    PROCEDURE :: receiveBuffer_d   => multifileRestartCollector_receiveBuffer_d
    PROCEDURE :: receiveBuffer_s   => multifileRestartCollector_receiveBuffer_s
    PROCEDURE :: receiveBuffer_int => multifileRestartCollector_receiveBuffer_int
    PROCEDURE, PRIVATE :: receiveBuffer_generic => multifileRestartCollector_receiveBuffer_generic
    GENERIC :: receiveBuffer => receiveBuffer_d, receiveBuffer_s, receiveBuffer_int
    PROCEDURE, PRIVATE :: checkArguments => multifileRestartCollector_checkArguments
  END TYPE t_MultifileRestartCollector

  CHARACTER(*), PARAMETER :: modname = "mo_multifile_restart_collector"

CONTAINS

  ! On the restart writers, this returns an array with the global
  ! indices of the points collected by this PE. The return value is
  ! not allocated on pure worker procs.
  !
  FUNCTION t_CollectorIndices_construct(me, localPointCount, decompInfo,        &
    &                                   destProc, sourceProcs, lthis_pe_active) &
    &   RESULT(globalIndices)
    !should have been of TYPE INTEGER, but CDI does NOT support
    !writing of integers:
    REAL(dp), POINTER :: globalIndices(:)
    CLASS(t_CollectorIndices),       INTENT(INOUT), TARGET :: me
    !the number of points PRESENT on this PE:
    INTEGER,                         INTENT(IN) :: localPointCount
    TYPE(t_grid_domain_decomp_info), INTENT(IN) :: decompInfo
    !rank within p_comm_work_restart to which this process should
    !send its DATA, must be the own rank on the restart writer
    !procs
    INTEGER,                         INTENT(IN) :: destProc
    !ranks within p_comm_work_restart from which DATA IS to be
    !collected, must include the own rank on the restart writer
    !procs; empty on pure work PEs
    INTEGER,                         INTENT(IN) :: sourceProcs(:)
    ! Flag. .FALSE. if this PE does not send any points (since it will
    ! write to another file).
    LOGICAL,                         INTENT(IN) :: lthis_pe_active
    CHARACTER(*), PARAMETER :: routine = modname//":t_CollectorIndices_construct"    
    INTEGER :: myRank, error, i, j
    REAL(dp), ALLOCATABLE :: sendBuffer_d(:)
    IF (timers_level >= 10)  CALL timer_start(timer_restart_indices_setup)
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
    IF (my_process_is_work() .AND. lthis_pe_active) THEN
      DO i = 1, localPointCount
        IF(decompInfo%owner_mask(idx_no(i), blk_no(i))) me%sendPointCount = me%sendPointCount + 1
      END DO
    END IF
    CALL alloc(me%sendIdx, me%sendPointCount)
    CALL alloc(me%sendBlk, me%sendPointCount)
    ALLOCATE(sendBuffer_d(me%sendPointCount), STAT=error)
    IF(error /= SUCCESS) CALL finish(routine, "memory allocation failure")
    !Fill sendIdx, sendBlk, and sendBuffer_d.
    IF (my_process_is_work() .AND. lthis_pe_active) THEN
      j = 1
      DO i = 1, localPointCount
        IF (decompInfo%owner_mask(idx_no(i), blk_no(i))) THEN
          me%sendIdx(j) = idx_no(i)
          me%sendBlk(j) = blk_no(i)
          sendBuffer_d(j) = decompInfo%glb_index(i)
          j = j + 1
        END IF
      END DO
    END IF
    !Collect the source point counts.
    IF (my_process_is_restart_writer()) THEN
      DO i = 1, me%sourceProcCount
        IF(me%sourceProcs(i) == myRank) THEN
          me%sourcePointCounts(i) = me%sendPointCount
        ELSE
          CALL p_recv(me%sourcePointCounts(i), me%sourceProcs(i), 0, comm = p_comm_work_restart)
        END IF
      END DO
    END IF
    IF (me%destProc /= myRank) THEN
      CALL p_send(me%sendPointCount, me%destProc, 0, comm = p_comm_work_restart)
    END IF
    !Compute the receivePointCount and allocate the globalIndices array.
    me%receivePointCount = SUM(me%sourcePointCounts(1:me%sourceProcCount))
    ALLOCATE(globalIndices(me%receivePointCount), STAT=error)
    IF(error /= SUCCESS) CALL finish(routine, "memory allocation failure")
    !Collect the global indices of the points written by this PE.
    CALL collectBuffer_blocking(me, sendBuffer_d, globalIndices)
    IF (timers_level >= 10)  CALL timer_stop(timer_restart_indices_setup)
  END FUNCTION t_CollectorIndices_construct

  SUBROUTINE t_CollectorIndices_finalize(me)
    CLASS(t_CollectorIndices), INTENT(INOUT) :: me
    CHARACTER(*), PARAMETER :: routine = modname//":t_CollectorIndices_finalize"
    INTEGER :: ierror

    DEALLOCATE(me%sourceProcs,       &
      &        me%sourcePointCounts, &
      &        me%sendIdx,           &
      &        me%sendBlk,  STAT=ierror)
    IF (ierror /= SUCCESS) CALL finish(routine, "memory deallocation failure")
  END SUBROUTINE t_CollectorIndices_finalize


  SUBROUTINE t_MultifileRestartCollector_construct(me, idx, glb_sendbuf, nlev, itype, &
    &                                              ioffset, sourceOffset)
    CLASS(t_MultifileRestartCollector),   INTENT(INOUT) :: me
    TYPE(t_CollectorIndices), TARGET,     INTENT(INOUT) :: idx
    TYPE(t_CollectorSendBuffer), TARGET,  INTENT(IN)    :: glb_sendbuf
    INTEGER,                              INTENT(IN)    :: nlev
    INTEGER,                              INTENT(IN)    :: itype
    INTEGER(i8),                          INTENT(INOUT) :: ioffset(:)
    !the start offset of a variable:
    INTEGER,                              INTENT(INOUT) :: sourceOffset(:,:)    
    CHARACTER(*), PARAMETER :: routine = modname//":t_multifileRestartCollector_construct"
    INTEGER     :: i, error, majver_mpi, minver_mpi

    IF (timers_level >= 10)  CALL timer_start(timer_restart_collector_setup)
    me%idx         => idx
    me%glb_sendbuf => glb_sendbuf
    CALL MPI_Get_version(majver_mpi, minver_mpi, error)
    IF (majver_mpi .GE. 3) have_MPI_RGet = .true.
    ! allocate the start offset:
    ALLOCATE(me%sourcePointStart(SIZE(sourceOffset,2)),  STAT=error)
    IF (error /= SUCCESS) CALL finish(routine, "memory allocation failure")
    IF (my_process_is_restart_writer()) THEN
      me%sourcePointStart(:) = sourceOffset(itype,:) 
      sourceOffset(itype,:)  = sourceOffset(itype,:) + nlev*idx%sourcePointCounts(:)
    END IF
    ! build send buffers for each variable (which point into the
    ! global send buffer):
    ALLOCATE(me%send_buffer_offset(nlev), STAT=error)
    IF(error /= SUCCESS) CALL finish(routine, "memory allocation failure")
    DO i = 1, nlev
      me%send_buffer_offset(i) = ioffset(itype) + (i - 1) * idx%sendPointCount
    END DO
    ioffset(itype) = ioffset(itype) + idx%sendPointCount * nlev
    IF (timers_level >= 10)  CALL timer_stop(timer_restart_collector_setup)
  END SUBROUTINE t_multifileRestartCollector_construct

  SUBROUTINE multifileRestartCollector_finalize(me)
    CLASS(t_MultifileRestartCollector), INTENT(INOUT) :: me
    CHARACTER(*), PARAMETER :: routine = modname//":multifileRestartCollector_finalize"
    INTEGER :: i, ierror

    IF (ALLOCATED(me%send_buffer_offset)) &
    &  DEALLOCATE(me%send_buffer_offset)
    NULLIFY(me%idx)
    DEALLOCATE(me%sourcePointStart, STAT=ierror)
    IF (ierror /= SUCCESS) CALL finish(routine, "memory deallocation failure")
  END SUBROUTINE multifileRestartCollector_finalize

  ! Collect the data within the send buffers to the respective
  ! outputData arrays.
  !
  ! Note: This blocking send/receive implementation is required only
  ! for setting up the index data structure.
  !
  SUBROUTINE collectBuffer_blocking(idx, sendBuffer_d, outputData)
    TYPE(t_CollectorIndices),  INTENT(IN)    :: idx
    REAL(dp),                  INTENT(INOUT) :: sendBuffer_d(:)
    REAL(dp),                  INTENT(INOUT) :: outputData(:)
    CHARACTER(*), PARAMETER :: routine = modname//":collectBuffer_blocking"
    INTEGER              :: myRank, i, j
    REAL(dp)             :: dummy_d(1)

    !sanity check
    myRank = p_comm_rank(p_comm_work_restart)
    IF(myRank == idx%destProc .AND. SIZE(outputData) /= idx%receivePointCount) THEN
      WRITE (message_text, *) "assertion failed: wrong buffer size (expected ",  &
        &                     idx%receivePointCount, ", got ", SIZE(outputData), &
        &                     "), check the calling routine"
      CALL finish(routine, message_text)
    END IF
    IF(idx%destProc /= myRank) THEN
      IF (idx%sendPointCount > 0) THEN
        CALL p_send(sendBuffer_d(1:idx%sendPointCount), idx%destProc, 0, &
          &         comm = p_comm_work_restart)
      ELSE
        CALL p_send(dummy_d(1), idx%destProc, 0, comm = p_comm_work_restart)
      END IF
    END IF
    !Collect the data on the writer PEs.
    j = 1
    outputData(:) = 0._dp
    DO i = 1, idx%sourceProcCount
      IF(idx%sourceProcs(i) == myRank) THEN
        outputData(j:j-1+idx%sourcePointCounts(i)) = sendBuffer_d(1:idx%sourcePointCounts(i))
      ELSE
        IF (idx%sourcePointCounts(i) > 0) THEN
          CALL p_recv(outputData(j:j-1+idx%sourcePointCounts(i)), idx%sourceProcs(i), 0, &
            &          comm = p_comm_work_restart)
        ELSE
          CALL p_recv(dummy_d(1), idx%sourceProcs(i), 0, comm = p_comm_work_restart)
        END IF
      END IF
      j = j + idx%sourcePointCounts(i)
    END DO
    !sanity check
    IF(j /= idx%receivePointCount + 1) CALL finish(routine, "assertion failed")
  END SUBROUTINE collectBuffer_blocking

  SUBROUTINE multifileRestartCollector_receiveBuffer_generic(me, levStart, used_size, levCount, outputData_s, outputData_d, outputData_i)
    CLASS(t_MultifileRestartCollector), INTENT(INOUT) :: me
    INTEGER,       INTENT(IN   )                      :: levStart, levCount
    INTEGER,       INTENT(INOUT)                      :: used_size
    REAL(KIND=sp), INTENT(INOUT), OPTIONAL, POINTER   :: outputData_s(:)
    REAL(KIND=dp), INTENT(INOUT), OPTIONAL, POINTER   :: outputData_d(:)
    INTEGER,       INTENT(INOUT), OPTIONAL, POINTER   :: outputData_i(:)
    CHARACTER(*), PARAMETER :: routine = modname//":multifileRestartCollector_receiveBuffer_generic"
#ifndef NOMPI
    CONTIGUOUS :: outputData_s, outputData_d, outputData_i
    INTEGER, PARAMETER :: n_openreqs_max = 32
    INTEGER(KIND=MPI_ADDRESS_KIND) :: target_disp
    INTEGER :: i, j, ierror, origin_addr, origin_count, origin_datatype, target_rank, &
      &        target_count, target_datatype, n_openreqs
    INTEGER :: get_reqs(n_openreqs_max), i_req, outType, stride, the_win, ddt

    stride = SUM(me%idx%sourcePointCounts(:))
    n_openreqs = 0
    IF (PRESENT(outputData_s)) THEN
      CALL ensureSize(outputData_s, stride*levCount, no_copy)
      outType = 1
      target_datatype = p_real_sp
      the_win = me%glb_sendbuf%win_s
!$OMP PARALLEL DO SCHEDULE(STATIC)
      DO i = 1, stride*levCount
        outputData_s(i) = 0._sp
      END DO
    END IF
    IF (PRESENT(outputData_d)) THEN
      CALL ensureSize(outputData_d, stride*levCount, no_copy)
      outType = 2
      target_datatype = p_real_dp
      the_win = me%glb_sendbuf%win_d
!$OMP PARALLEL DO SCHEDULE(STATIC)
      DO i = 1, stride*levCount
        outputData_d(i) = 0._dp
      END DO
    END IF
    IF (PRESENT(outputData_i)) THEN
      CALL ensureSize(outputData_i, stride*levCount, no_copy)
      outType = 3
      target_datatype = p_int
      the_win = me%glb_sendbuf%win_int
!$OMP PARALLEL DO SCHEDULE(STATIC)
      DO i = 1, stride*levCount
        outputData_i(i) = 0
      END DO
    END IF
    origin_addr  = 1
    origin_count = 1
    DO i = 1, me%idx%sourceProcCount
      ddt = MPI_DATATYPE_NULL
      IF (me%idx%sourcePointCounts(i) > 0) THEN
        CALL MPI_Type_vector(levCount, me%idx%sourcePointCounts(i), stride, &
                             target_datatype, ddt, ierror)
        IF (ierror /= MPI_SUCCESS) CALL finish(routine, "MPI error!")
        CALL MPI_Type_commit(ddt, ierror)
        IF (ierror /= MPI_SUCCESS) CALL finish(routine, "MPI error!")
        origin_datatype = ddt
        target_rank     = me%idx%sourceProcs(i)
        target_disp     = me%sourcePointStart(i) + (levStart-1)*me%idx%sourcePointCounts(i)
        target_count    = me%idx%sourcePointCounts(i) * levCount
        IF (have_MPI_RGet) THEN
          IF (n_openreqs .LT. n_openreqs_max) THEN
            n_openreqs = n_openreqs + 1
            i_req = n_openreqs
          ELSE
            CALL MPI_Waitany(n_openreqs, get_reqs, i_req, MPI_STATUS_IGNORE, ierror)
            get_reqs(i_req) = MPI_REQUEST_NULL
            IF (ierror /= MPI_SUCCESS) CALL finish(routine, "MPI error!")
          END IF
          SELECT CASE(outType)
            CASE(1)
              CALL MPI_RGET(outputData_s(origin_addr), origin_count, origin_datatype, &
                            target_rank, target_disp,  target_count, target_datatype, &
                            the_win, get_reqs(i_req), ierror)
            CASE(2)
              CALL MPI_RGET(outputData_d(origin_addr), origin_count, origin_datatype, &
                            target_rank, target_disp,  target_count, target_datatype, &
                            the_win, get_reqs(i_req), ierror)
            CASE(3)
              CALL MPI_RGET(outputData_i(origin_addr), origin_count, origin_datatype, &
                            target_rank, target_disp,  target_count, target_datatype, &
                            the_win, get_reqs(i_req), ierror)
          END SELECT
        ELSE
          CALL MPI_WIN_LOCK(MPI_LOCK_SHARED, target_rank, MPI_MODE_NOCHECK, &
            &               the_win, ierror)
          IF (ierror /= MPI_SUCCESS) CALL finish(routine, "MPI error!")
          SELECT CASE(outType)
            CASE(1)
              CALL MPI_GET(outputData_s(origin_addr), origin_count, origin_datatype, &
                           target_rank, target_disp,  target_count, target_datatype, &
                           the_win, ierror)
            CASE(2)
              CALL MPI_GET(outputData_d(origin_addr), origin_count, origin_datatype, &
                           target_rank, target_disp,  target_count, target_datatype, &
                           the_win, ierror)
            CASE(3)
              CALL MPI_GET(outputData_i(origin_addr), origin_count, origin_datatype, &
                           target_rank, target_disp,  target_count, target_datatype, &
                           the_win, ierror)
          END SELECT
        END IF
        IF (ierror /= MPI_SUCCESS) CALL finish(routine, "MPI error!")
        CALL MPI_Type_free(ddt, ierror)
        IF (ierror /= MPI_SUCCESS) CALL finish(routine, "MPI error!")
        origin_addr = origin_addr + me%idx%sourcePointCounts(i)
      END IF
    END DO
    IF (have_MPI_Rget) THEN
      CALL MPI_WAITALL(n_openreqs, get_reqs ,MPI_STATUSES_IGNORE, ierror)
      IF (ierror /= MPI_SUCCESS) CALL finish(routine, "MPI error!")
    ELSE
      DO i = 1, me%idx%sourceProcCount
        IF (me%idx%sourcePointCounts(i) > 0) &
          &  CALL MPI_WIN_UNLOCK(me%idx%sourceProcs(i), the_win, ierror)
        IF (ierror /= MPI_SUCCESS) CALL finish(routine, "MPI error!")
      ENDDO
    ENDIF
    IF (origin_addr .NE. stride + 1) CALL finish(routine, "Inconsistent!")
    used_size = levCount * stride
#else
    CALL finish(routine, "Not implemented!")
#endif
  END SUBROUTINE multifileRestartCollector_receiveBuffer_generic

  SUBROUTINE multifileRestartCollector_receiveBuffer_s(me, levStart, outputData, used_size, levCount)
    CLASS(t_MultifileRestartCollector), INTENT(INOUT) :: me
    INTEGER,                            INTENT(IN   ) :: levStart, levCount
    REAL(sp), CONTIGUOUS,      POINTER, INTENT(INOUT) :: outputData(:)
    INTEGER,                            INTENT(  OUT) :: used_size

    CALL me%receiveBuffer_generic(levStart, used_size, levCount, outputData_s=outputData)
  END SUBROUTINE multifileRestartCollector_receiveBuffer_s

  SUBROUTINE multifileRestartCollector_receiveBuffer_d(me, levStart, outputData, used_size, levCount)
    CLASS(t_MultifileRestartCollector), INTENT(INOUT) :: me
    INTEGER,                            INTENT(IN   ) :: levStart, levCount
    REAL(dp), CONTIGUOUS,      POINTER, INTENT(INOUT) :: outputData(:)
    INTEGER,                            INTENT(  OUT) :: used_size

    CALL me%receiveBuffer_generic(levStart, used_size, levCount, outputData_d=outputData)
  END SUBROUTINE multifileRestartCollector_receiveBuffer_d

  SUBROUTINE multifileRestartCollector_receiveBuffer_int(me, levStart, outputData, used_size, levCount)
    CLASS(t_MultifileRestartCollector), INTENT(INOUT) :: me
    INTEGER,                            INTENT(IN   ) :: levStart, levCount
    INTEGER, CONTIGUOUS,       POINTER, INTENT(INOUT) :: outputData(:)
    INTEGER,                            INTENT(  OUT) :: used_size

    CALL me%receiveBuffer_generic(levStart, used_size, levCount, outputData_i=outputData)
  END SUBROUTINE multifileRestartCollector_receiveBuffer_int

  SUBROUTINE multifileRestartCollector_sendField_generic(me, levStart, inputData_s, &
                                                         inputData_d, inputData_i,  &
                                                         levCount_in)
    CLASS(t_MultifileRestartCollector), INTENT(INOUT) :: me
    INTEGER, INTENT(IN)            :: levStart
    REAL(sp), INTENT(IN), OPTIONAL :: inputData_s(:,:)
    REAL(dp), INTENT(IN), OPTIONAL :: inputData_d(:,:)
    INTEGER, INTENT(IN),  OPTIONAL :: inputData_i(:,:)
    INTEGER, INTENT(IN), OPTIONAL  :: levCount_in
    CHARACTER(*), PARAMETER :: routine = modname//":multifileRestartCollector_sendField_d"
#ifndef NOMPI
    INTEGER     :: i, ierror, ilev, optinCount, inType, levCount
    INTEGER(i8) :: ioffset, glb_size
    CHARACTER(LEN=3) :: my_type_is
    
    CALL me%checkArguments(levStart)
    optinCount = 0
    IF (PRESENT(inputData_s)) THEN
      inType = 1
      optinCount = optinCount + 1
      my_type_is = "sp"
      glb_size = SIZE(me%glb_sendbuf%sendBuffer_s)
    END IF
    IF (PRESENT(inputData_d)) THEN
      inType = 2
      optinCount = optinCount + 1
      my_type_is = "dp"
      glb_size = SIZE(me%glb_sendbuf%sendBuffer_d)
    END IF
    IF (PRESENT(inputData_i)) THEN
      inType = 3
      optinCount = optinCount + 1
      my_type_is = "int"
      glb_size = SIZE(me%glb_sendbuf%sendBuffer_int)
    END IF
    IF(optinCount /= 1) CALL finish(routine, "assertion failed")
    IF (PRESENT(levCount_in)) THEN 
      levCount = levCount_in 
      CALL me%checkArguments(levStart - 1 + levCount)
    END IF
    IF (.NOT. ASSOCIATED(me%glb_sendbuf%sendBuffer_d)) THEN
      CALL finish(routine, "Unassociated send buffer!")
    END IF
    IF (me%idx%sendPointCount > 0) THEN
      DO ilev = levStart, levStart -1 + levCount
        !Fill the send buffer.
        ioffset = me%send_buffer_offset(ilev)
        IF (glb_size < ioffset + me%idx%sendPointCount) THEN
          WRITE (message_text,*) "Internal error: SIZE(sendBuffer_generic(",TRIM(my_type_is),")) = ", glb_size,    &
            &                    " < ioffset + me%idx%sendPointCount) = ", ioffset + me%idx%sendPointCount, &
            &                    "; sendPointCount = ", me%idx%sendPointCount, "; ilev=", ilev
          CALL finish(routine, message_text)
        END IF
        SELECT CASE(inType)
          CASE(1)
!$OMP PARALLEL DO SCHEDULE(STATIC)
            DO i = 1, me%idx%sendPointCount
              me%glb_sendbuf%sendBuffer_s(ioffset+i) = inputData_s(me%idx%sendIdx(i), me%idx%sendBlk(i))
            END DO
          CASE(2)
!$OMP PARALLEL DO SCHEDULE(STATIC)
            DO i = 1, me%idx%sendPointCount
              me%glb_sendbuf%sendBuffer_d(ioffset+i) = inputData_d(me%idx%sendIdx(i), me%idx%sendBlk(i))
            END DO
          CASE(3)
!$OMP PARALLEL DO SCHEDULE(STATIC)
            DO i = 1, me%idx%sendPointCount
              me%glb_sendbuf%sendBuffer_int(ioffset+i) = inputData_i(me%idx%sendIdx(i), me%idx%sendBlk(i))
            END DO
        END SELECT
      END DO
    END IF
#else
    CALL finish(routine, "Not implemented!")
#endif
  END SUBROUTINE multifileRestartCollector_sendField_generic

  SUBROUTINE multifileRestartCollector_sendField_s(me, inputData, levStart, levCount_in)
    CLASS(t_MultifileRestartCollector), INTENT(INOUT) :: me
    REAL(sp), INTENT(IN)    :: inputData(:,:)
    INTEGER, INTENT(IN)    :: levStart
    INTEGER, OPTIONAL, INTENT(IN) :: levCount_in
    INTEGER :: levCount

    levCount = 1
    IF (PRESENT(levCount_in)) levCount = levCount_in
    CALL me%sendField_generic(levStart, inputData_s=inputData, levCount_in=levCount)
  END SUBROUTINE multifileRestartCollector_sendField_s

  SUBROUTINE multifileRestartCollector_sendField_d(me, inputData, levStart, levCount_in)
    CLASS(t_MultifileRestartCollector), INTENT(INOUT) :: me
    REAL(dp), INTENT(IN)    :: inputData(:,:)
    INTEGER, INTENT(IN)    :: levStart
    INTEGER, OPTIONAL, INTENT(IN) :: levCount_in
    INTEGER :: levCount

    levCount = 1
    IF (PRESENT(levCount_in)) levCount = levCount_in
    CALL me%sendField_generic(levStart, inputData_d=inputData, levCount_in=levCount)
  END SUBROUTINE multifileRestartCollector_sendField_d

  SUBROUTINE multifileRestartCollector_sendField_int(me, inputData, levStart, levCount_in)
    CLASS(t_MultifileRestartCollector), INTENT(INOUT) :: me
    INTEGER, INTENT(IN)    :: inputData(:,:)
    INTEGER, INTENT(IN)    :: levStart
    INTEGER, OPTIONAL, INTENT(IN) :: levCount_in
    INTEGER :: levCount

    levCount = 1
    IF (PRESENT(levCount_in)) levCount = levCount_in
    CALL me%sendField_generic(levStart, inputData_i=inputData, levCount_in=levCount)
  END SUBROUTINE multifileRestartCollector_sendField_int

  SUBROUTINE multifileRestartCollector_checkArguments(me, ilev)
    CLASS(t_MultifileRestartCollector), INTENT(INOUT) :: me
    INTEGER, INTENT(IN)    :: ilev
    CHARACTER(*), PARAMETER :: routine = modname//"::multifileRestartCollector_checkArguments"

    IF (.NOT. ASSOCIATED(me%idx)) THEN
      CALL finish(routine, "Internal error!")
    END IF
    IF (.NOT. ALLOCATED(me%send_buffer_offset)) THEN
      CALL finish(routine, "Internal error: me%send_buffer_offset unallocated!")
    END IF
    IF ((ilev < 1) .OR. (ilev > SIZE(me%send_buffer_offset))) THEN
      WRITE (message_text, *) "Internal error: Invalid level index: ", ilev, "!"
      CALL finish(routine, message_text)
    END IF
  END SUBROUTINE multifileRestartCollector_checkArguments

  SUBROUTINE t_CollectorSendBuffer_construct(this, isize, idx)
    CLASS(t_CollectorSendBuffer), INTENT(INOUT) :: this
    INTEGER(i8),                  INTENT(IN)    :: isize(:)
    TYPE(t_CollectorIndices),     INTENT(INOUT) :: idx
    CHARACTER(*), PARAMETER :: routine = modname//":t_CollectorSendBuffer_add_domain"
#ifdef NOMPI
    CALL finish(routine, "Not implemented!")
#else
    INTEGER, PARAMETER :: typeId(3) = (/ REAL_T, SINGLE_T, INT_T /)
    INTEGER :: ierr, splitKey, wcGrp, wcGrp_size, i, j, ii, iStart, myRank
    INTEGER, ALLOCATABLE, DIMENSION(:) :: rank_list, rank_map, &
                                          tmpSrcCounts, tmpSrcProcs

    this%allocd      = .false.
    this%handshaked  = .false.
    this%win_posted  = .false.
    this%win_started = .false.
    DO i = 1, 3
      this%winSizes(i) = isize(typeId(i))
      NULLIFY(this%winPtr(i)%p)
    END DO
    myRank = p_comm_rank(p_comm_work_restart)
    splitKey = MERGE(1, myRank + 2, my_process_is_restart_writer())
    CALL MPI_Comm_split(p_comm_work_restart, idx%destProc, splitKey, this%winComm, ierr)
    IF (ierr /= MPI_SUCCESS) CALL finish(routine, "MPI error!")
    idx%destProc = 0
    CALL MPI_Comm_group(this%winComm, wcGrp, ierr)
    IF (ierr /= MPI_SUCCESS) CALL finish(routine, "MPI error!")
    CALL MPI_Group_size(wcGrp, wcGrp_size, ierr)
    IF (ierr /= MPI_SUCCESS) CALL finish(routine, "MPI error!")
    ALLOCATE(rank_list(wcGrp_size), rank_map(wcGrp_size))
    rank_list(1:wcGrp_size) = (/ (i, i = 0, wcGrp_size-1) /)
    CALL MPI_Group_incl(wcGrp, 1, rank_list(1:1), this%winSvGroup, ierr)
    IF (ierr /= MPI_SUCCESS) CALL finish(routine, "MPI error!")
    CALL MPI_Gather(myRank, 1, p_int, rank_map, 1, p_int, 0, this%winComm, ierr)
    IF (ierr /= MPI_SUCCESS) CALL finish(routine, "MPI error!")
    IF (my_process_is_restart_writer()) THEN
      iStart = MERGE(1, 2, my_process_is_work())
      CALL MPI_Group_incl(wcGrp, wcGrp_size-iStart+1, rank_list(iStart:wcGrp_size), &
        &                 this%winClGroup, ierr)
      IF (ierr /= MPI_SUCCESS) CALL finish(routine, "MPI error!")
      ALLOCATE(tmpSrcProcs (idx%sourceProcCount))
      tmpSrcProcs(:)  = idx%sourceProcs(:)
      DO i = 1, idx%sourceProcCount
        DO j = 1, wcGrp_size
          IF (rank_map(j) .EQ. tmpSrcProcs(i)) ii = j - 1
        END DO
        idx%sourceProcs(i) = ii
      END DO
      DEALLOCATE(tmpSrcProcs)
    ELSE
      this%winClGroup = MPI_GROUP_NULL
    END IF
    CALL MPI_Group_free(wcGrp, ierr)
    IF (ierr /= MPI_SUCCESS) CALL finish(routine, "MPI error!")
    DEALLOCATE(rank_list, rank_map)
#endif
  END SUBROUTINE t_CollectorSendBuffer_construct

  SUBROUTINE t_CollectorSendBuffer_handshake(this)
    CLASS(t_CollectorSendBuffer), INTENT(INOUT) :: this
    CHARACTER(*), PARAMETER :: routine = modname//":t_CollectorSendBuffer_construct"
#ifdef NOMPI
    CALL finish(routine, "Not implemented!")
#else
    INTEGER(KIND=MPI_ADDRESS_KIND) :: memSize(1), memBytes, typeBytes, typeLB, dpBytes
    TYPE(c_ptr) :: cMemPtr
    INTEGER :: ierr, iType, typeIdMPI(3)
    INTEGER, PARAMETER :: idummy(1) = (/ 1 /)

    typeIdMPI(1:3) = (/ p_real_dp, p_real_sp, p_int /)
    CALL MPI_Type_get_extent(p_real_dp, typeLB, dpBytes, ierr)
    DO iType = 1, 3
      CALL MPI_Type_get_extent(typeIdMPI(iType), typeLB, typeBytes, ierr)
      IF (ierr /= MPI_SUCCESS) CALL finish(routine, "MPI error!")
      memSize(1) = MAX(INT(this%winSizes(iType), MPI_ADDRESS_KIND), 1_MPI_ADDRESS_KIND)
      memBytes = MAX(memSize(1) * typeBytes, dpBytes)
! as of MPI3.0 standard the following MPI_Alloc_mem interface must be present, if
! the compiler provides ISO_C_BINDING
      CALL MPI_Alloc_mem(memBytes, MPI_INFO_NULL, cMemPtr, ierr)
      CALL C_F_POINTER(cMemPtr, this%winPtr(iType)%p, idummy)
      IF (ierr /= MPI_SUCCESS) CALL finish(routine, "MPI error!")
      SELECT CASE (iType)
      CASE (1)
        CALL C_F_POINTER(cMemPtr, this%sendBuffer_d, INT(memSize))
        CALL MPI_Win_create(this%sendBuffer_d, memBytes, INT(typeBytes), &
          &                  MPI_INFO_NULL, this%winComm, this%win_d, ierr)
      CASE (2)
        CALL C_F_POINTER(cMemPtr, this%sendBuffer_s, INT(memSize))
        CALL MPI_Win_create(this%sendBuffer_s, memBytes, INT(typeBytes), &
          &                 MPI_INFO_NULL, this%winComm, this%win_s, ierr)
      CASE (3)
        CALL C_F_POINTER(cMemPtr, this%sendBuffer_int, INT(memSize))
        CALL MPI_Win_create(this%sendBuffer_int, memBytes, INT(typeBytes), &
          &                 MPI_INFO_NULL, this%winComm, this%win_int, ierr)
      END SELECT
      IF (ierr /= MPI_SUCCESS) CALL finish(routine, "MPI error!")
    END DO
    CALL MPI_Comm_free(this%winComm, ierr)
    IF (ierr /= MPI_SUCCESS) CALL finish(routine, "MPI error!")
    this%winComm    = MPI_COMM_NULL
    this%handshaked = .true.
    this%allocd     = .true.
#endif
  END SUBROUTINE t_CollectorSendBuffer_handshake

  SUBROUTINE t_CollectorSendBuffer_finalize(this)
    CLASS(t_CollectorSendBuffer), INTENT(INOUT) :: this
    CHARACTER(*), PARAMETER :: routine = modname//":t_CollectorSendBuffer_finalize"
#ifndef NOMPI
    INTEGER :: ierr, i

    IF (this%handshaked) THEN
      IF(this%win_posted .OR. this%win_started) &
        CALL this%start_local_access()
      CALL MPI_Win_free(this%win_d, ierr)
      IF (ierr /= MPI_SUCCESS) CALL finish(routine, "MPI error!")
      CALL MPI_Win_free(this%win_s, ierr)
      IF (ierr /= MPI_SUCCESS) CALL finish(routine, "MPI error!")
      CALL MPI_Win_free(this%win_int, ierr)
      IF (ierr /= MPI_SUCCESS) CALL finish(routine, "MPI error!")
      this%handshaked  = .false.
    END IF
    IF (this%allocd) THEN
! cumbersome: MPI_Free_mem only accepts Fortran variables/references NOT TYPE(C_PTR)/MPI_ADDRESS_KIND as base-address, while MPI_Alloc_mem returns an TYPE(C_PTR)/MPI_ADDRESS_KIND !!
      DO i = 1, 3
        CALL MPI_Free_mem(this%winPtr(i)%p, ierr)
        IF (ierr /= MPI_SUCCESS) CALL finish(routine, "MPI error!")
      END DO
      NULLIFY(this%sendBuffer_d, this%sendBuffer_s, this%sendBuffer_int)
      this%allocd = .false.
    END IF
    IF (this%winClGroup .NE. MPI_GROUP_NULL) THEN
      CALL MPI_Group_free(this%winClGroup, ierr)
      IF (ierr /= MPI_SUCCESS) CALL finish(routine, "MPI error!")
    END IF
    CALL MPI_Group_free(this%winSvGroup, ierr)
    IF (ierr /= MPI_SUCCESS) CALL finish(routine, "MPI error!")
    IF (this%winComm .NE. MPI_COMM_NULL) THEN
      CALL MPI_Comm_free(this%winComm, ierr)
      IF (ierr /= MPI_SUCCESS) CALL finish(routine, "MPI error!")
      this%winComm    = MPI_COMM_NULL
    END IF
    this%win_posted  = .false.
    this%win_started = .false.
#endif
  END SUBROUTINE t_CollectorSendBuffer_finalize

  SUBROUTINE t_CollectorSendBuffer_start_local_access(this)
    CLASS(t_CollectorSendBuffer), INTENT(INOUT) :: this
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::t_CollectorSendBuffer_start_local_access"
#ifndef NOMPI
    INTEGER :: iassert, ierr

    IF (.NOT.this%handshaked) CALL this%handshake()
    IF (.NOT.this%win_started .OR. .NOT.this%win_posted) &
      CALL this%start_remote_access()
    IF (this%win_started) THEN
      CALL MPI_Win_complete(this%win_d, ierr)
      IF (ierr /= MPI_SUCCESS) CALL finish(routine, "MPI error!")
      CALL MPI_Win_complete(this%win_s, ierr)
      IF (ierr /= MPI_SUCCESS) CALL finish(routine, "MPI error!")
      CALL MPI_Win_complete(this%win_int, ierr)
      IF (ierr /= MPI_SUCCESS) CALL finish(routine, "MPI error!")
      this%win_started = .false.
    END IF
    IF (this%win_posted) THEN
      CALL MPI_Win_wait(this%win_d, ierr)
      IF (ierr /= MPI_SUCCESS) CALL finish(routine, "MPI error!")
      CALL MPI_Win_wait(this%win_s, ierr)
      IF (ierr /= MPI_SUCCESS) CALL finish(routine, "MPI error!")
      CALL MPI_Win_wait(this%win_int, ierr)
      IF (ierr /= MPI_SUCCESS) CALL finish(routine, "MPI error!")
      this%win_posted  = .false.
    ENDIF
#else
    CALL finish(routine, "Not implemented!")
#endif
  END SUBROUTINE t_CollectorSendBuffer_start_local_access

  ! Start MPI epoch where windows may not be accessed locally and
  ! not via MPI_PUT operations.
  !
  SUBROUTINE t_CollectorSendBuffer_start_remote_access(this)
    CLASS(t_CollectorSendBuffer), INTENT(INOUT) :: this
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::t_CollectorSendBuffer_start_remote_access"
#ifndef NOMPI
    INTEGER :: iassert, ierr

    IF (.NOT.this%handshaked) CALL finish(routine, "there is no buffer to get!")
    IF (this%win_started .OR. this%win_posted) RETURN ! nothing to do
    IF (my_process_is_work()) THEN
      iassert = MPI_MODE_NOPUT
      CALL MPI_Win_post(this%winSvGroup, iassert, this%win_d, ierr)
      IF (ierr /= MPI_SUCCESS) CALL finish(routine, "MPI error!")
      CALL MPI_Win_post(this%winSvGroup, iassert, this%win_s, ierr)
      IF (ierr /= MPI_SUCCESS) CALL finish(routine, "MPI error!")
      CALL MPI_Win_post(this%winSvGroup, iassert, this%win_int, ierr)
      IF (ierr /= MPI_SUCCESS) CALL finish(routine, "MPI error!")
      this%win_posted  = .true.
    END IF
    IF (my_process_is_restart_writer()) THEN
      iassert = 0
      CALL MPI_Win_start(this%winClGroup, iassert, this%win_d, ierr)
      IF (ierr /= MPI_SUCCESS) CALL finish(routine, "MPI error!")
      CALL MPI_Win_start(this%winClGroup, iassert, this%win_s, ierr)
      IF (ierr /= MPI_SUCCESS) CALL finish(routine, "MPI error!")
      CALL MPI_Win_start(this%winClGroup, iassert, this%win_int, ierr)
      IF (ierr /= MPI_SUCCESS) CALL finish(routine, "MPI error!")
      this%win_started  = .true.
    END IF
#else
    CALL finish(routine, "Not implemented!")
#endif
  END SUBROUTINE t_CollectorSendBuffer_start_remote_access

END MODULE mo_multifile_restart_collector
