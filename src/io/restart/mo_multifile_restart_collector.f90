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

  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_ptr, c_intptr_t, c_f_pointer

#ifndef NOMPI
  USE mpi
#endif
  USE mo_communication,       ONLY: idx_no, blk_no
  USE mo_decomposition_tools, ONLY: t_grid_domain_decomp_info
  USE mo_exception,           ONLY: finish, message_text
  USE mo_fortran_tools,       ONLY: alloc
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
  

  ! buffers for linearizing the data before sending it to the writer
  ! process
  TYPE :: t_CollectorSendBuffer
    
    REAL(dp), POINTER :: sendBuffer_d(:)
    REAL(sp), POINTER :: sendBuffer_s(:)
    INTEGER,  POINTER :: sendBuffer_int(:)  
    
    ! One-Sided MPI window handles
    INTEGER :: mpi_win_sendbuf_d, mpi_win_sendbuf_s, mpi_win_sendbuf_int  

  CONTAINS

    PROCEDURE :: construct    => t_CollectorSendBuffer_construct
    PROCEDURE :: finalize     => t_CollectorSendBuffer_finalize

    PROCEDURE :: start_local_access  => t_CollectorSendBuffer_start_local_access
    PROCEDURE :: start_remote_access => t_CollectorSendBuffer_start_remote_access

    PROCEDURE, PRIVATE :: allocate_mem => t_CollectorSendBuffer_allocate_mem
  END TYPE t_CollectorSendBuffer


  TYPE :: t_CollectorVarBuffer

    INTEGER(i8) :: send_buffer_offset
    INTEGER     :: send_buffer_size

  CONTAINS
  
    PROCEDURE :: construct => t_CollectorVarBuffer_construct
    PROCEDURE :: finalize  => t_CollectorVarBuffer_finalize

  END TYPE t_CollectorVarBuffer


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
    TYPE(t_CollectorVarBuffer), ALLOCATABLE  :: buf(:)

    ! on receiving PEs: start index in send buffer for every source
    ! PE:
    INTEGER, ALLOCATABLE :: sourcePointStart(:)
   
  CONTAINS
    PROCEDURE :: construct => t_MultifileRestartCollector_construct
    PROCEDURE :: finalize  => multifileRestartCollector_finalize   

    PROCEDURE :: sendField_d   => multifileRestartCollector_sendField_d
    PROCEDURE :: sendField_s   => multifileRestartCollector_sendField_s
    PROCEDURE :: sendField_int => multifileRestartCollector_sendField_int
    GENERIC :: sendField => sendField_d, sendField_s, sendField_int

    PROCEDURE :: receiveBuffer_d   => multifileRestartCollector_receiveBuffer_d
    PROCEDURE :: receiveBuffer_s   => multifileRestartCollector_receiveBuffer_s
    PROCEDURE :: receiveBuffer_int => multifileRestartCollector_receiveBuffer_int
    GENERIC :: receiveBuffer => receiveBuffer_d, receiveBuffer_s, receiveBuffer_int

    PROCEDURE, PRIVATE :: checkArguments => multifileRestartCollector_checkArguments
   
  END TYPE t_MultifileRestartCollector

 
 
  CHARACTER(*), PARAMETER :: modname = "mo_multifile_restart_collector"

CONTAINS

  SUBROUTINE t_CollectorVarBuffer_construct(me, ioffset, isize)
    CLASS(t_CollectorVarBuffer), INTENT(INOUT) :: me
    INTEGER(i8),                 INTENT(INOUT) :: ioffset
    INTEGER,                     INTENT(IN)    :: isize

    me%send_buffer_offset = ioffset
    me%send_buffer_size   = isize
    ioffset = ioffset + isize
  END SUBROUTINE t_CollectorVarBuffer_construct


  SUBROUTINE t_CollectorVarBuffer_finalize(me)
    CLASS(t_CollectorVarBuffer), INTENT(INOUT) :: me
    ! do nothing
  END SUBROUTINE t_CollectorVarBuffer_finalize


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
    DO i = 1, me%sourceProcCount
      IF(me%sourceProcs(i) == myRank) THEN
        me%sourcePointCounts(i) = me%sendPointCount
      ELSE
        CALL p_recv(me%sourcePointCounts(i), me%sourceProcs(i), 0, comm = p_comm_work_restart)
      END IF
    END DO
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


  SUBROUTINE t_multifileRestartCollector_construct(me, idx, glb_sendbuf, nlev, itype, &
    &                                              ioffset, sourceOffset)

    CLASS(t_multifileRestartCollector),   INTENT(INOUT) :: me
    TYPE(t_CollectorIndices), TARGET,     INTENT(INOUT) :: idx
    TYPE(t_CollectorSendBuffer), TARGET,  INTENT(IN)    :: glb_sendbuf
    INTEGER,                              INTENT(IN)    :: nlev
    INTEGER,                              INTENT(IN)    :: itype
    INTEGER(i8),                          INTENT(INOUT) :: ioffset(:)
    !the start offset of a variable:
    INTEGER,                              INTENT(INOUT) :: sourceOffset(:,:)    

    CHARACTER(*), PARAMETER :: routine = modname//":t_multifileRestartCollector_construct"
    INTEGER     :: i, error
    INTEGER(i8) :: itype_offset

    IF (timers_level >= 10)  CALL timer_start(timer_restart_collector_setup)

    me%idx         => idx
    me%glb_sendbuf => glb_sendbuf

    ! allocate the start offset:

    ALLOCATE(me%sourcePointStart(SIZE(sourceOffset,2)),  STAT=error)
    IF (error /= SUCCESS) CALL finish(routine, "memory allocation failure")

    IF (my_process_is_restart_writer()) THEN
      me%sourcePointStart(:) = sourceOffset(itype,:) 
      sourceOffset(itype,:)  = sourceOffset(itype,:) + nlev*idx%sourcePointCounts(:)
    END IF

    ! build send buffers for each variable (which point into the
    ! global send buffer):
      
    ALLOCATE(me%buf(nlev), STAT=error)
    IF(error /= SUCCESS) CALL finish(routine, "memory allocation failure")
      
    DO i=1,nlev
      ! note: the following call modifies the offset!
      itype_offset = ioffset(itype)
      CALL me%buf(i)%construct(itype_offset, idx%sendPointCount)
      ioffset(itype) = itype_offset
    END DO

    IF (timers_level >= 10)  CALL timer_stop(timer_restart_collector_setup)

  END SUBROUTINE t_multifileRestartCollector_construct


  SUBROUTINE multifileRestartCollector_finalize(me)
    CLASS(t_MultifileRestartCollector), INTENT(INOUT) :: me

    CHARACTER(*), PARAMETER :: routine = modname//":multifileRestartCollector_finalize"
    INTEGER :: i, ierror

    IF (ALLOCATED(me%buf)) THEN
      DO i=1,SIZE(me%buf)
        CALL me%buf(i)%finalize()
      END DO
      DEALLOCATE(me%buf, STAT=ierror)
      IF(ierror /= SUCCESS) CALL finish(routine, "memory deallocation failure")
    END IF
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


  SUBROUTINE multifileRestartCollector_receiveBuffer_d(me, ilev, outputData)
    CLASS(t_MultifileRestartCollector), INTENT(INOUT) :: me
    INTEGER,                    INTENT(IN)    :: ilev
    REAL(dp),                   INTENT(INOUT) :: outputData(:)

    CHARACTER(*), PARAMETER :: routine = modname//":multifileRestartCollector_receiveBuffer_d"
#ifndef NOMPI
    INTEGER(KIND=MPI_ADDRESS_KIND) :: target_disp
    INTEGER :: i, j, ierror, origin_addr, origin_count, origin_datatype, target_rank, &
      &        target_count, target_datatype

    !Collect the data on the writer PEs.
    j = 1
    outputData(:) = 0._dp
    DO i = 1, me%idx%sourceProcCount


      origin_addr     = j
      origin_count    = me%idx%sourcePointCounts(i)
      origin_datatype = p_real_dp
      target_rank     = me%idx%sourceProcs(i)
      target_disp     = me%sourcePointStart(i) + (ilev-1)*me%idx%sourcePointCounts(i)
      target_count    = origin_count
      target_datatype = origin_datatype

      IF (me%idx%sourcePointCounts(i) > 0) THEN
        CALL MPI_WIN_LOCK(MPI_LOCK_SHARED, me%idx%sourceProcs(i), MPI_MODE_NOCHECK, &
          &               me%glb_sendbuf%mpi_win_sendbuf_d,ierror)
        IF (ierror /= SUCCESS) CALL finish(routine, "MPI error!")

        CALL MPI_GET(outputData(origin_addr), origin_count, origin_datatype, target_rank, target_disp, &
          &          target_count, target_datatype, me%glb_sendbuf%mpi_win_sendbuf_d, ierror)
        IF (ierror /= SUCCESS) CALL finish(routine, "MPI error!")

        CALL MPI_WIN_UNLOCK(me%idx%sourceProcs(i), me%glb_sendbuf%mpi_win_sendbuf_d, ierror)
        IF (ierror /= SUCCESS) CALL finish(routine, "MPI error!")
      ENDIF

      j = j + me%idx%sourcePointCounts(i)
    END DO

    !sanity check
    IF(j /= me%idx%receivePointCount + 1) CALL finish(routine, "assertion failed")
#else
    CALL finish(routine, "Not implemented!")
#endif
  END SUBROUTINE multifileRestartCollector_receiveBuffer_d


  SUBROUTINE multifileRestartCollector_receiveBuffer_s(me, ilev, outputData)
    CLASS(t_MultifileRestartCollector), INTENT(INOUT) :: me
    INTEGER,                    INTENT(IN)    :: ilev
    REAL(sp),                   INTENT(INOUT) :: outputData(:)

    CHARACTER(*), PARAMETER :: routine = modname//":multifileRestartCollector_receiveBuffer_s"
#ifndef NOMPI
    INTEGER(KIND=MPI_ADDRESS_KIND) :: target_disp
    INTEGER :: i, j, ierror, origin_addr, origin_count, origin_datatype, target_rank, &
      &        target_count, target_datatype

    !Collect the data on the writer PEs.
    j = 1
    outputData(:) = 0._sp
    DO i = 1, me%idx%sourceProcCount
      

      origin_addr     = j
      origin_count    = me%idx%sourcePointCounts(i)
      origin_datatype = p_real_sp
      target_rank     = me%idx%sourceProcs(i)
      target_disp     = me%sourcePointStart(i) + (ilev-1)*me%idx%sourcePointCounts(i)
      target_count    = origin_count
      target_datatype = origin_datatype

      IF (me%idx%sourcePointCounts(i) > 0) THEN
        CALL MPI_WIN_LOCK(MPI_LOCK_SHARED, me%idx%sourceProcs(i), MPI_MODE_NOCHECK, &
          &               me%glb_sendbuf%mpi_win_sendbuf_s, ierror)
        IF (ierror /= SUCCESS) CALL finish(routine, "MPI error!")

        CALL MPI_GET(outputData(origin_addr), origin_count, origin_datatype, target_rank, target_disp, &
          &          target_count, target_datatype, me%glb_sendbuf%mpi_win_sendbuf_s, ierror)
        IF (ierror /= SUCCESS) CALL finish(routine, "MPI error!")

        CALL MPI_WIN_UNLOCK(me%idx%sourceProcs(i), me%glb_sendbuf%mpi_win_sendbuf_s, ierror)
        IF (ierror /= SUCCESS) CALL finish(routine, "MPI error!")
      ENDIF

      j = j + me%idx%sourcePointCounts(i)
    END DO

    !sanity check
    IF(j /= me%idx%receivePointCount + 1) CALL finish(routine, "assertion failed")
#else
    CALL finish(routine, "Not implemented!")
#endif
  END SUBROUTINE multifileRestartCollector_receiveBuffer_s


  SUBROUTINE multifileRestartCollector_receiveBuffer_int(me, ilev, outputData)
    CLASS(t_MultifileRestartCollector), INTENT(INOUT) :: me
    INTEGER,                    INTENT(IN)    :: ilev
    INTEGER,                    INTENT(INOUT) :: outputData(:)

    CHARACTER(*), PARAMETER :: routine = modname//":multifileRestartCollector_receiveBuffer_int"
#ifndef NOMPI
    INTEGER(KIND=MPI_ADDRESS_KIND) :: target_disp
    INTEGER :: i, j, ierror, origin_addr, origin_count, origin_datatype, target_rank, &
      &        target_count, target_datatype

    !Collect the data on the writer PEs.
    j = 1
    outputData(:) = 0
    DO i = 1, me%idx%sourceProcCount
      

      origin_addr     = j
      origin_count    = me%idx%sourcePointCounts(i)
      origin_datatype = p_int
      target_rank     = me%idx%sourceProcs(i)
      target_disp     = me%sourcePointStart(i) + (ilev-1)*me%idx%sourcePointCounts(i)
      target_count    = origin_count
      target_datatype = origin_datatype

      IF (me%idx%sourcePointCounts(i) > 0) THEN
        CALL MPI_WIN_LOCK(MPI_LOCK_SHARED, me%idx%sourceProcs(i), MPI_MODE_NOCHECK, &
          &               me%glb_sendbuf%mpi_win_sendbuf_int, ierror)
        IF (ierror /= SUCCESS) CALL finish(routine, "MPI error!")

        CALL MPI_GET(outputData(origin_addr), origin_count, origin_datatype, target_rank, target_disp, &
          &          target_count, target_datatype, me%glb_sendbuf%mpi_win_sendbuf_int, ierror)
        IF (ierror /= SUCCESS) CALL finish(routine, "MPI error!")

        CALL MPI_WIN_UNLOCK(me%idx%sourceProcs(i), me%glb_sendbuf%mpi_win_sendbuf_int, ierror)
        IF (ierror /= SUCCESS) CALL finish(routine, "MPI error!")
      ENDIF

      j = j + me%idx%sourcePointCounts(i)
    END DO

    !sanity check
    IF(j /= me%idx%receivePointCount + 1) CALL finish(routine, "assertion failed")
#else
    CALL finish(routine, "Not implemented!")
#endif
  END SUBROUTINE multifileRestartCollector_receiveBuffer_int


  SUBROUTINE multifileRestartCollector_sendField_d(me, inputData, ilev)
    CLASS(t_MultifileRestartCollector), INTENT(INOUT) :: me
    REAL(dp), INTENT(IN)    :: inputData(:,:)
    INTEGER,  INTENT(IN)    :: ilev

    CHARACTER(*), PARAMETER :: routine = modname//":multifileRestartCollector_sendField_d"
#ifndef NOMPI
    INTEGER     :: i, myRank, ierror
    INTEGER(i8) :: ioffset

    CALL me%checkArguments(ilev)
    IF (.NOT. ASSOCIATED(me%glb_sendbuf%sendBuffer_d)) THEN
      CALL finish(routine, "Unassociated send buffer!")
    END IF

    IF (me%idx%sendPointCount > 0) THEN
      myRank = p_comm_rank(p_comm_work_restart)
      CALL MPI_WIN_LOCK(MPI_LOCK_EXCLUSIVE, myRank, MPI_MODE_NOCHECK, me%glb_sendbuf%mpi_win_sendbuf_d, ierror)
      IF (ierror /= SUCCESS) CALL finish(routine, "MPI error!")

      !Fill the send buffer.
      ioffset = me%buf(ilev)%send_buffer_offset
      IF (SIZE(me%glb_sendbuf%sendBuffer_d) < ioffset + me%idx%sendPointCount) THEN
        WRITE (message_text,*) "Internal error: SIZE(sendBuffer_d) = ", SIZE(me%glb_sendbuf%sendBuffer_d),    &
          &                    " < ioffset + me%idx%sendPointCount) = ", ioffset + me%idx%sendPointCount, &
          &                    "; sendPointCount = ", me%idx%sendPointCount, "; ilev=", ilev
        CALL finish(routine, message_text)
      END IF

      DO i = 1, me%idx%sendPointCount
        me%glb_sendbuf%sendBuffer_d(ioffset+i) = inputData(me%idx%sendIdx(i), me%idx%sendBlk(i))
      END DO

      CALL MPI_WIN_UNLOCK(myRank, me%glb_sendbuf%mpi_win_sendbuf_d, ierror)
      IF (ierror /= SUCCESS) CALL finish(routine, "MPI error!")
    ENDIF
#else
    CALL finish(routine, "Not implemented!")
#endif
  END SUBROUTINE multifileRestartCollector_sendField_d


  SUBROUTINE multifileRestartCollector_sendField_s(me, inputData, ilev)
    CLASS(t_MultifileRestartCollector), INTENT(INOUT) :: me
    REAL(sp), INTENT(IN)    :: inputData(:,:)
    INTEGER,  INTENT(IN)    :: ilev

    CHARACTER(*), PARAMETER :: routine = modname//":multifileRestartCollector_sendField_s"
#ifndef NOMPI
    INTEGER     :: i, myRank, ierror
    INTEGER(i8) :: ioffset

    CALL me%checkArguments(ilev)
    IF (.NOT. ASSOCIATED(me%glb_sendbuf%sendBuffer_s)) THEN
      CALL finish(routine, "Unassociated send buffer!")
    END IF

    IF (me%idx%sendPointCount > 0) THEN
      myRank = p_comm_rank(p_comm_work_restart)
      CALL MPI_WIN_LOCK(MPI_LOCK_EXCLUSIVE, myRank, MPI_MODE_NOCHECK, me%glb_sendbuf%mpi_win_sendbuf_s, ierror)
      IF (ierror /= SUCCESS) CALL finish(routine, "MPI error!")

      !Fill the send buffer.
      ioffset = me%buf(ilev)%send_buffer_offset
      IF (SIZE(me%glb_sendbuf%sendBuffer_s) < ioffset + me%idx%sendPointCount) THEN
        WRITE (message_text,*) "Internal error: SIZE(sendBuffer_s) = ", SIZE(me%glb_sendbuf%sendBuffer_s),    &
          &                    " < ioffset + me%idx%sendPointCount) = ", ioffset + me%idx%sendPointCount, &
          &                    "; sendPointCount = ", me%idx%sendPointCount, "; ilev=", ilev
        CALL finish(routine, message_text)
      END IF

      DO i = 1, me%idx%sendPointCount
        me%glb_sendbuf%sendBuffer_s(ioffset+i) = inputData(me%idx%sendIdx(i), me%idx%sendBlk(i))
      END DO

      CALL MPI_WIN_UNLOCK(myRank, me%glb_sendbuf%mpi_win_sendbuf_s, ierror)
      IF (ierror /= SUCCESS) CALL finish(routine, "MPI error!")
    ENDIF
#else
    CALL finish(routine, "Not implemented!")
#endif
  END SUBROUTINE multifileRestartCollector_sendField_s


  SUBROUTINE multifileRestartCollector_sendField_int(me, inputData, ilev)
    CLASS(t_MultifileRestartCollector), INTENT(INOUT) :: me
    INTEGER, INTENT(IN)    :: inputData(:,:)
    INTEGER, INTENT(IN)    :: ilev

    CHARACTER(*), PARAMETER :: routine = modname//":multifileRestartCollector_sendField_int"
#ifndef NOMPI
    INTEGER     :: i, myRank, ierror
    INTEGER(i8) :: ioffset

    CALL me%checkArguments(ilev)
    IF (.NOT. ASSOCIATED(me%glb_sendbuf%sendBuffer_int)) THEN
      CALL finish(routine, "Unassociated send buffer!")
    END IF

    IF (me%idx%sendPointCount > 0) THEN
      myRank = p_comm_rank(p_comm_work_restart)
      CALL MPI_WIN_LOCK(MPI_LOCK_EXCLUSIVE, myRank, MPI_MODE_NOCHECK, me%glb_sendbuf%mpi_win_sendbuf_int, ierror)
      IF (ierror /= SUCCESS) CALL finish(routine, "MPI error!")

      !Fill the send buffer.
      ioffset = me%buf(ilev)%send_buffer_offset
      IF (SIZE(me%glb_sendbuf%sendBuffer_int) < ioffset + me%idx%sendPointCount) THEN
        WRITE (message_text,*) "Internal error: SIZE(sendBuffer_int) = ", SIZE(me%glb_sendbuf%sendBuffer_int), &
          &                    " < ioffset + me%idx%sendPointCount) = ", ioffset + me%idx%sendPointCount,  &
          &                    "; sendPointCount = ", me%idx%sendPointCount, "; ilev=", ilev
        CALL finish(routine, message_text)
      END IF

      DO i = 1, me%idx%sendPointCount
        me%glb_sendbuf%sendBuffer_int(ioffset+i) = inputData(me%idx%sendIdx(i), me%idx%sendBlk(i))
      END DO

      CALL MPI_WIN_UNLOCK(myRank, me%glb_sendbuf%mpi_win_sendbuf_int, ierror)
      IF (ierror /= SUCCESS) CALL finish(routine, "MPI error!")
    ENDIF
#else
    CALL finish(routine, "Not implemented!")
#endif
  END SUBROUTINE multifileRestartCollector_sendField_int


  SUBROUTINE multifileRestartCollector_checkArguments(me, ilev)
    CLASS(t_MultifileRestartCollector), INTENT(INOUT) :: me
    INTEGER, INTENT(IN)    :: ilev

    CHARACTER(*), PARAMETER :: routine = modname//"::multifileRestartCollector_checkArguments"

    IF (.NOT. ASSOCIATED(me%idx)) THEN
      CALL finish(routine, "Internal error!")
    END IF
    IF (.NOT. ALLOCATED(me%buf)) THEN
      CALL finish(routine, "Internal error: me%buf unallocated!")
    END IF
    IF ((ilev < 1) .OR. (ilev > SIZE(me%buf))) THEN
      WRITE (message_text, *) "Internal error: Invalid level index: ", ilev, "!"
      CALL finish(routine, message_text)
    END IF
    IF (me%idx%sendPointCount > me%buf(ilev)%send_buffer_size) THEN
      WRITE (message_text,*) "Internal error: send buffer has wrong size; sendPointCount=", &
        &                    me%idx%sendPointCount, "!"
      CALL finish(routine, message_text)
    END IF
  END SUBROUTINE multifileRestartCollector_checkArguments


  SUBROUTINE t_CollectorSendBuffer_construct(buf, isize)
    CLASS(t_CollectorSendBuffer), INTENT(INOUT) :: buf
    INTEGER(i8),                  INTENT(IN)    :: isize(:)

    CHARACTER(*), PARAMETER :: routine = modname//":t_CollectorSendBuffer_construct"

    CALL buf%allocate_mem(isize(REAL_T),   REAL_T,   buf%mpi_win_sendbuf_d)
    CALL buf%allocate_mem(isize(SINGLE_T), SINGLE_T, buf%mpi_win_sendbuf_s)
    CALL buf%allocate_mem(isize(INT_T),    INT_T,    buf%mpi_win_sendbuf_int)
  END SUBROUTINE t_CollectorSendBuffer_construct


  SUBROUTINE t_CollectorSendBuffer_finalize(buf)
    CLASS(t_CollectorSendBuffer), INTENT(INOUT) :: buf

    CHARACTER(*), PARAMETER :: routine = modname//":t_CollectorSendBuffer_finalize"
#ifndef NOMPI
    INTEGER :: mpierr

    CALL MPI_WIN_FREE( buf%mpi_win_sendbuf_d,   mpierr )
    IF (mpierr /= SUCCESS) CALL finish(routine, "MPI error!")

    CALL MPI_WIN_FREE( buf%mpi_win_sendbuf_s,   mpierr )
    IF (mpierr /= SUCCESS) CALL finish(routine, "MPI error!")

    CALL MPI_WIN_FREE( buf%mpi_win_sendbuf_int, mpierr )
    IF (mpierr /= SUCCESS) CALL finish(routine, "MPI error!")
#endif
  END SUBROUTINE t_CollectorSendBuffer_finalize


  SUBROUTINE t_CollectorSendBuffer_start_local_access(buf)
    CLASS(t_CollectorSendBuffer), INTENT(INOUT) :: buf

    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::t_CollectorSendBuffer_start_local_access"
#ifndef NOMPI
    INTEGER :: iassert, mpierr

    iassert = MPI_MODE_NOSTORE
    CALL MPI_WIN_FENCE(iassert, buf%mpi_win_sendbuf_d,    mpierr)
    IF (mpierr /= SUCCESS) CALL finish(routine, "MPI error!")
    CALL MPI_WIN_FENCE(iassert, buf%mpi_win_sendbuf_s,    mpierr)
    IF (mpierr /= SUCCESS) CALL finish(routine, "MPI error!")
    CALL MPI_WIN_FENCE(iassert, buf%mpi_win_sendbuf_int,  mpierr)
    IF (mpierr /= SUCCESS) CALL finish(routine, "MPI error!")
#else
    CALL finish(routine, "Not implemented!")
#endif

  END SUBROUTINE t_CollectorSendBuffer_start_local_access


  ! Start MPI epoch where windows may not be accessed locally and
  ! not via MPI_PUT operations.
  !
  SUBROUTINE t_CollectorSendBuffer_start_remote_access(buf)
    CLASS(t_CollectorSendBuffer), INTENT(INOUT) :: buf

    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::t_CollectorSendBuffer_start_remote_access"
#ifndef NOMPI
    INTEGER :: iassert, mpierr

    iassert = IOR(MPI_MODE_NOPUT, MPI_MODE_NOPRECEDE)
    CALL MPI_WIN_FENCE(iassert, buf%mpi_win_sendbuf_d,  mpierr)
    IF (mpierr /= SUCCESS) CALL finish(routine, "MPI error!")
    CALL MPI_WIN_FENCE(iassert, buf%mpi_win_sendbuf_s,  mpierr)
    IF (mpierr /= SUCCESS) CALL finish(routine, "MPI error!")
    CALL MPI_WIN_FENCE(iassert, buf%mpi_win_sendbuf_int,  mpierr)
    IF (mpierr /= SUCCESS) CALL finish(routine, "MPI error!")
#else
    CALL finish(routine, "Not implemented!")
#endif
  END SUBROUTINE t_CollectorSendBuffer_start_remote_access



  !------------------------------------------------------------------------------------------------
  !> allocate amount of memory needed with MPI_Alloc_mem
  !
  !  @note Implementation for non-Cray pointers
  !
  SUBROUTINE t_CollectorSendBuffer_allocate_mem(buf, isize, itype, mpi_win)
    CLASS(t_CollectorSendBuffer), INTENT(INOUT) :: buf
    INTEGER(i8),              INTENT(IN)    :: isize
    INTEGER,                  INTENT(IN)    :: itype
    INTEGER,                  INTENT(INOUT) :: mpi_win
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::t_CollectorSendBuffer_allocate_mem"

#ifndef NOMPI
    INTEGER (KIND=MPI_ADDRESS_KIND) :: mem_size, mem_bytes
    TYPE(c_ptr)                     :: c_mem_ptr
    INTEGER                         :: mpierr, nbytes

    ! Get the amount of bytes per REAL*8 or REAL*4 variable (as used in MPI
    ! communication)
    SELECT CASE (itype)
    CASE (REAL_T)
      CALL MPI_Type_extent(p_real_dp, nbytes, mpierr)
    CASE (SINGLE_T)
      CALL MPI_Type_extent(p_real_sp, nbytes, mpierr)
    CASE (INT_T)
      CALL MPI_Type_extent(p_int,     nbytes, mpierr)
    CASE DEFAULT
      CALL finish(routine, "Internal error!")
    END SELECT
    IF (mpierr /= SUCCESS) CALL finish(routine, "MPI error!")

    ! For the IO PEs the amount of memory needed is 0 - allocate at least 1 word there:
    mem_size  = isize
    mem_bytes = MAX(mem_size,1_i8)*INT(nbytes,i8)

    ! TYPE(c_ptr) and INTEGER(KIND=MPI_ADDRESS_KIND) do NOT necessarily have the same size!!!
    ! So check if at least c_intptr_t and MPI_ADDRESS_KIND are the same, else we may get
    ! into deep, deep troubles!
    ! There is still a slight probability that TYPE(c_ptr) does not have the size indicated
    ! by c_intptr_t since the standard only requires c_intptr_t is big enough to hold pointers
    ! (so it may be bigger than a pointer), but I hope no vendor screws up its ISO_C_BINDING
    ! in such a way!!!
    ! If c_intptr_t<=0, this type is not defined and we can't do this check, of course.

    IF(c_intptr_t > 0 .AND. c_intptr_t /= MPI_ADDRESS_KIND) &
     & CALL finish(routine,'c_intptr_t /= MPI_ADDRESS_KIND, too dangerous to proceed!')

    CALL MPI_ALLOC_MEM(mem_bytes, MPI_INFO_NULL, c_mem_ptr, mpierr)
    IF (mpierr /= SUCCESS) CALL finish(routine, "MPI error!")

    ! Create memory window for communication
    SELECT CASE (itype)
    CASE (REAL_T)
      NULLIFY(buf%sendBuffer_d)
      CALL C_F_POINTER(c_mem_ptr, buf%sendBuffer_d, (/ mem_size /) )
      CALL MPI_WIN_CREATE( buf%sendBuffer_d, mem_bytes, nbytes, MPI_INFO_NULL,&
        &                  p_comm_work_restart, mpi_win, mpierr )
      IF (mpierr /= SUCCESS) CALL finish(routine, "MPI error!")
      !
    CASE (SINGLE_T)
      NULLIFY(buf%sendBuffer_s)
      CALL C_F_POINTER(c_mem_ptr, buf%sendBuffer_s, (/ mem_size /) )
      CALL MPI_WIN_CREATE( buf%sendBuffer_s, mem_bytes, nbytes, MPI_INFO_NULL,&
        &                  p_comm_work_restart, mpi_win, mpierr )
      IF (mpierr /= SUCCESS) CALL finish(routine, "MPI error!")
      !
    CASE (INT_T)
      NULLIFY(buf%sendBuffer_int)
      CALL C_F_POINTER(c_mem_ptr, buf%sendBuffer_int, (/ mem_size /) )
      CALL MPI_WIN_CREATE( buf%sendBuffer_int, mem_bytes, nbytes, MPI_INFO_NULL,&
        &                  p_comm_work_restart, mpi_win, mpierr )
      IF (mpierr /= SUCCESS) CALL finish(routine, "MPI error!")
      !
    CASE DEFAULT
      CALL finish(routine, "Internal error!")
    END SELECT

    ! No local operations prior to this epoch, so give an assertion
    CALL MPI_WIN_FENCE(MPI_MODE_NOPRECEDE, mpi_win, mpierr)
    IF (mpierr /= SUCCESS) CALL finish(routine, "MPI error!")

#else
    CALL finish(routine, "Not implemented!")
#endif
  END SUBROUTINE t_CollectorSendBuffer_allocate_mem

END MODULE mo_multifile_restart_collector
