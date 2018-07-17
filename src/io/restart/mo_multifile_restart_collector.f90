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

  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_ptr, C_F_POINTER, C_NULL_PTR, C_LOC

#ifndef NOMPI
  USE mpi
#else
#define MPI_ADDRESS_KIND i8
#endif
  USE mo_communication,       ONLY: idx_no, blk_no
  USE mo_decomposition_tools, ONLY: t_grid_domain_decomp_info
  USE mo_exception,           ONLY: finish, message_text
  USE mo_fortran_tools,       ONLY: alloc, ensureSize, no_copy, &
    &                               t_ptr_2d, t_ptr_2d_sp, t_ptr_2d_int
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
    &       t_CollectorSendBuffer, commonBuf_t, dataPtrs_t
  
  PRIVATE
  
  TYPE :: commonBuf_t
    REAL(dp), POINTER :: d(:)
    REAL(sp), POINTER :: s(:)
    INTEGER,  POINTER :: i(:)
  END TYPE commonBuf_t

  TYPE :: dataPtrs_t
    TYPE(t_ptr_2d),     ALLOCATABLE :: d(:)
    TYPE(t_ptr_2d_sp),  ALLOCATABLE :: s(:)
    TYPE(t_ptr_2d_int), ALLOCATABLE :: i(:)
  END TYPE dataPtrs_t

  TYPE :: ptr_arr_t
    REAL(KIND=dp), POINTER :: p(:)
  END TYPE ptr_arr_t

  TYPE :: t_CollectorSendBuffer
    TYPE(commonBuf_t) :: sendBuffer
    INTEGER :: win
    INTEGER(KIND=MPI_ADDRESS_KIND), ALLOCATABLE :: tOffSv(:)
    INTEGER(KIND=MPI_ADDRESS_KIND) :: facDpSp, facIntSp
    INTEGER(KIND=MPI_ADDRESS_KIND), PRIVATE, ALLOCATABLE :: tOffCl(:)
    INTEGER,          PRIVATE :: wComm, wClGrp, wSvGrp
    INTEGER(KIND=i8), PRIVATE :: wSizes(3)
    TYPE(ptr_arr_t),  PRIVATE :: wPtr
    LOGICAL,          PRIVATE :: allocd, handshakd, &
      &                          wPosted, wStarted
  CONTAINS
    PROCEDURE :: construct           => CollectorSendBuffer_construct
    PROCEDURE, PRIVATE :: handshake  => CollectorSendBuffer_handshake
    PROCEDURE :: finalize            => CollectorSendBuffer_finalize
    PROCEDURE :: start_local_access  => CollectorSendBuffer_start_local_access
    PROCEDURE :: start_remote_access => CollectorSendBuffer_start_remote_access
  END TYPE t_CollectorSendBuffer

  TYPE :: t_CollectorIndices
    INTEGER :: recvPntCnt, srcProcCnt, sendPntCnt, destProc
    INTEGER, ALLOCATABLE :: srcProc(:), srcPntCnts(:), sendIdx(:), sendBlk(:)
  CONTAINS
    PROCEDURE :: construct => CollectorIndices_construct
    PROCEDURE :: finalize  => CollectorIndices_finalize
  END TYPE t_CollectorIndices

  TYPE :: t_CollectorIndices_ptr
    TYPE(t_CollectorIndices), POINTER :: p
  END TYPE t_CollectorIndices_ptr

  TYPE :: t_MultifileRestartCollector
    PRIVATE
    TYPE(t_CollectorIndices_ptr) :: idx(3)
    TYPE(t_CollectorSendBuffer), POINTER :: glb_sendbuf
    INTEGER, ALLOCATABLE :: vGrid(:), vType(:)
  CONTAINS
    PROCEDURE :: construct => multifileRestartCollector_construct
    PROCEDURE :: defVar    => multifileRestartCollector_defVar
    PROCEDURE :: finalize  => multifileRestartCollector_finalize
    PROCEDURE :: sendField => multifileRestartCollector_sendField
    PROCEDURE :: fetch     => multifileRestartCollector_fetch
  END TYPE t_MultifileRestartCollector

  CHARACTER(*), PARAMETER :: modname = "mo_multifile_restart_collector"

CONTAINS

  ! On the restart writers, this returns an array with the global
  ! indices of the points collected by this PE. The return value is
  ! not allocated on pure worker procs.
  FUNCTION CollectorIndices_construct(me, localPointCount, decompInfo,        &
    &                                   destProc, sourceProcs, lthis_pe_active) &
    &   RESULT(globalIndices)
    REAL(dp), POINTER :: globalIndices(:)
    CLASS(t_CollectorIndices),       INTENT(INOUT), TARGET :: me
    INTEGER,                         INTENT(IN) :: localPointCount
    TYPE(t_grid_domain_decomp_info), INTENT(IN) :: decompInfo
    INTEGER,                         INTENT(IN) :: destProc
    INTEGER,                         INTENT(IN) :: sourceProcs(:)
    LOGICAL,                         INTENT(IN) :: lthis_pe_active
    CHARACTER(*), PARAMETER :: routine = modname//":CollectorIndices_construct"    
    INTEGER :: myRank, error, i, j
    REAL(dp), ALLOCATABLE :: sendBuffer_d(:)
    IF (timers_level >= 10)  CALL timer_start(timer_restart_indices_setup)
    myRank = p_comm_rank(p_comm_work_restart)
    !Copy the information about the processes IN our group, AND
    !ALLOCATE the sourceProcCount dependent arrays.
    me%srcProcCnt = SIZE(sourceProcs)
    CALL alloc(me%srcProc, me%srcProcCnt)
    CALL alloc(me%srcPntCnts, me%srcProcCnt)
    me%srcProc(1:SIZE(sourceProcs)) = sourceProcs(:)
    me%destProc = destProc
    !Compute the number of points we want to send, AND ALLOCATE
    !the sendPointCount dependent arrays.
    me%sendPntCnt = 0
    IF (my_process_is_work() .AND. lthis_pe_active) THEN
      DO i = 1, localPointCount
        IF(decompInfo%owner_mask(idx_no(i), blk_no(i))) &
          & me%sendPntCnt = me%sendPntCnt + 1
      END DO
    END IF
    CALL alloc(me%sendIdx, me%sendPntCnt)
    CALL alloc(me%sendBlk, me%sendPntCnt)
    ALLOCATE(sendBuffer_d(me%sendPntCnt), STAT=error)
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
      DO i = 1, me%srcProcCnt
        IF(me%srcProc(i) == myRank) THEN
          me%srcPntCnts(i) = me%sendPntCnt
        ELSE
          CALL p_recv(me%srcPntCnts(i), me%srcProc(i), 0, comm = p_comm_work_restart)
        END IF
      END DO
    END IF
    IF (me%destProc /= myRank .AND. lthis_pe_active) THEN
      CALL p_send(me%sendPntCnt, me%destProc, 0, comm = p_comm_work_restart)
    END IF
    !Compute the receivePointCount and allocate the globalIndices array.
    me%recvPntCnt = SUM(me%srcPntCnts(1:me%srcProcCnt))
    ALLOCATE(globalIndices(me%recvPntCnt), STAT=error)
    !Collect the global indices of the points written by this PE.
    IF (lthis_pe_active) THEN
      CALL collectBuffer_blocking(me, sendBuffer_d, globalIndices)
    END IF
    IF (timers_level >= 10)  CALL timer_stop(timer_restart_indices_setup)
  END FUNCTION CollectorIndices_construct

  SUBROUTINE CollectorIndices_finalize(me)
    CLASS(t_CollectorIndices), INTENT(INOUT) :: me
    CHARACTER(*), PARAMETER :: routine = modname//":CollectorIndices_finalize"
    INTEGER :: ierror

    DEALLOCATE(me%srcProc, me%srcPntCnts, me%sendIdx, me%sendBlk,  STAT=ierror)
    IF (ierror /= SUCCESS) CALL finish(routine, "memory deallocation failure")
  END SUBROUTINE CollectorIndices_finalize

  SUBROUTINE MultifileRestartCollector_construct(this, idx_cell, idx_edge, &
      &                                            idx_vert, glb_sendbuf, nVar)
    CLASS(t_MultifileRestartCollector),  INTENT(INOUT) :: this
    TYPE(t_CollectorIndices), TARGET,    INTENT(INOUT) :: idx_cell, idx_edge, idx_vert
    TYPE(t_CollectorSendBuffer), TARGET, INTENT(IN)    :: glb_sendbuf
    INTEGER,                             INTENT(IN)    :: nVar
    CHARACTER(*), PARAMETER :: routine = modname//":multifileRestartCollector_construct"
    INTEGER     :: ierr

    IF (timers_level >= 10)  CALL timer_start(timer_restart_collector_setup)
    this%idx(1)%p => idx_cell
    this%idx(2)%p => idx_edge
    this%idx(3)%p => idx_vert
    this%glb_sendbuf => glb_sendbuf
    ALLOCATE(this%vGrid(nVar), this%vType(nVar), STAT=ierr)
    IF (timers_level >= 10)  CALL timer_stop(timer_restart_collector_setup)
  END SUBROUTINE multifileRestartCollector_construct

  SUBROUTINE multifileRestartCollector_defVar(this, iVar, nLevs, iType, iGrid, iOffset)
  USE mo_cdi_constants,               ONLY: GRID_UNSTRUCTURED_CELL, &
    &                                       GRID_UNSTRUCTURED_EDGE, &
    &                                       GRID_UNSTRUCTURED_VERT
    CLASS(t_MultifileRestartCollector), INTENT(INOUT) :: this
    INTEGER,                            INTENT(IN   ) :: iVar, nLevs, iType, iGrid
    INTEGER(KIND=i8),                   INTENT(INOUT) :: iOffset(:)
    CHARACTER(*), PARAMETER :: routine = modname//":multifileRestartCollector_defVar"

    IF (timers_level >= 10)  CALL timer_start(timer_restart_collector_setup)
    SELECT CASE(iGrid)
    CASE(GRID_UNSTRUCTURED_CELL)
      this%vGrid(iVar) = 1
    CASE(GRID_UNSTRUCTURED_EDGE)
      this%vGrid(iVar) = 2
    CASE(GRID_UNSTRUCTURED_VERT)
      this%vGrid(iVar) = 3
    END SELECT
    SELECT CASE(iType)
    CASE(REAL_T)
      this%vType(iVar) = 1
    CASE(SINGLE_T)
      this%vType(iVar) = 2
    CASE(INT_T)
      this%vType(iVar) = 3
    END SELECT
    ioffset(itype) = ioffset(itype) + this%idx(this%vGrid(iVar))%p%sendPntCnt * nLevs
    IF (timers_level >= 10)  CALL timer_stop(timer_restart_collector_setup)
  END SUBROUTINE multifileRestartCollector_defVar

  SUBROUTINE multifileRestartCollector_finalize(this)
    CLASS(t_MultifileRestartCollector), INTENT(INOUT) :: this
    CHARACTER(*), PARAMETER :: routine = modname//":multifileRestartCollector_finalize"

    NULLIFY(this%idx(1)%p, this%idx(2)%p, this%idx(3)%p)
    IF (ALLOCATED(this%vGrid)) DEALLOCATE(this%vGrid, this%vType)
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

    myRank = p_comm_rank(p_comm_work_restart)
    IF(myRank == idx%destProc .AND. SIZE(outputData) /= idx%recvPntCnt) THEN
      WRITE (message_text, *) "assertion failed: wrong buffer size (expected ",  &
        &                     idx%recvPntCnt, ", got ", SIZE(outputData), &
        &                     "), check the calling routine"
      CALL finish(routine, message_text)
    END IF
    IF(idx%destProc /= myRank) THEN
      IF (idx%sendPntCnt > 0) THEN
        CALL p_send(sendBuffer_d(1:idx%sendPntCnt), idx%destProc, 0, &
          &         comm = p_comm_work_restart)
      ELSE
        CALL p_send(dummy_d(1), idx%destProc, 0, comm = p_comm_work_restart)
      END IF
    END IF
    !Collect the data on the writer PEs.
    j = 1
    outputData(:) = 0._dp
    DO i = 1, idx%srcProcCnt
      IF(idx%srcProc(i) == myRank) THEN
        outputData(j:j-1+idx%srcPntCnts(i)) = sendBuffer_d(1:idx%srcPntCnts(i))
      ELSE
        IF (idx%srcPntCnts(i) > 0) THEN
          CALL p_recv(outputData(j:j-1+idx%srcPntCnts(i)), idx%srcProc(i), 0, &
            &          comm = p_comm_work_restart)
        ELSE
          CALL p_recv(dummy_d(1), idx%srcProc(i), 0, comm = p_comm_work_restart)
        END IF
      END IF
      j = j + idx%srcPntCnts(i)
    END DO
    IF(j /= idx%recvPntCnt + 1) CALL finish(routine, "assertion failed")
  END SUBROUTINE collectBuffer_blocking

  SUBROUTINE multifileRestartCollector_fetch(me, vStart, vCnt, vSkip, lCnt, &
                                             used_size, output, srcOffset)
    CLASS(t_MultifileRestartCollector), INTENT(INOUT) :: me
    INTEGER,           INTENT(IN   ) :: vStart, vCnt, lCnt(vCnt)
    LOGICAL,           INTENT(IN   ) :: vSkip(vCnt)
    INTEGER,           INTENT(  OUT) :: used_size(vCnt)
    TYPE(commonBuf_t), INTENT(INOUT) :: output
    INTEGER(KIND=i8), ALLOCATABLE, INTENT(INOUT) :: srcOffset(:,:)
    CHARACTER(*), PARAMETER :: routine = &
      & modname//":multifileRestartCollector_receiveBuffer_generic"
#ifndef NOMPI
    INTEGER, PARAMETER :: nOpenReqMax = 16, addr_k = MPI_ADDRESS_KIND
    INTEGER(KIND=addr_k) :: t_dsp, fDsp, tBy, tLB, rDsp(vCnt), rD(vCnt)
    INTEGER :: i, ierr, o_dType, dType, nOpenReq, iReq, bSize, &
      &        getReq(nOpenReqMax), iType, iV, iV_o, iV_s
    INTEGER, DIMENSION(vCnt) :: ones, stride, tmp_ddt
    INTEGER, ALLOCATABLE :: pct(:,:), t_cnt(:)
    TYPE(t_CollectorIndices), POINTER :: idx

    idx => me%idx(me%vGrid(vStart))%p
    iType = me%vType(vStart)
    IF (.NOT.ALLOCATED(srcOffset)) THEN
      ALLOCATE(srcOffset(idx%srcProcCnt, 3))
      srcOffset(:,:) = 0_i8
    END IF
    bSize = 0
    ones(:) = 1
    stride(:) = 0
    ALLOCATE(t_cnt(idx%srcProcCnt), pct(idx%srcProcCnt, vCnt))
    t_cnt(:) = 0
    rDsp(:) = 0_addr_k
    DO iV_o = 1, vCnt
      IF (vSkip(iV_o)) CYCLE
      iV = vStart - 1 + iV_o
      pct(:,iV_o) = me%idx(me%vGrid(iV))%p%srcPntCnts(:)
      t_cnt(:) = t_cnt(:) + pct(:,iV_o) * lCnt(iV_o)
      stride(iV_o) = SUM(pct(:,iV_o))
      rDsp(iV_o) = bSize
      bSize = bSize + stride(iV_o) * lCnt(iV_o)
    END DO
    nOpenReq = 0
    SELECT CASE(iType)
    CASE(1)
      CALL ensureSize(output%d, bSize, no_copy)
      fDsp = me%glb_sendbuf%facDpSp
      dType = p_real_dp
    CASE(2)
      CALL ensureSize(output%s, bSize, no_copy)
      fDsp = 1_addr_k
      dType = p_real_sp
    CASE(3)
      CALL ensureSize(output%i, bSize, no_copy)
      fDsp = me%glb_sendbuf%facIntSp
      dType = p_int
    END SELECT
    CALL MPI_Type_get_extent(dType, tLB, tBy, ierr)
    DO i = 1, idx%srcProcCnt
      o_dtype = MPI_DATATYPE_NULL
      IF (t_cnt(i) > 0) THEN
        iV_s = 0
        DO iV_o = 1, vCnt
          IF (vSkip(iV_o)) CYCLE
          iV = vStart - 1 + iV_o
          iV_s = iV_s + 1
          CALL MPI_Type_vector(lCnt(iV_o), pct(i,iV_o), stride(iV_o), dType, tmp_ddt(iV_s), ierr)
          rD(iV_s) = rDsp(iV_o) * tBy
        END DO
        CALL MPI_Type_create_struct(iV_s, ones, rD, tmp_ddt, o_dType, ierr)
        DO iV_o = 1, iV_s
          CALL MPI_Type_free(tmp_ddt(iV_o), ierr)
        END DO
        CALL MPI_Type_commit(o_dType, ierr)
        t_dsp  = srcOffset(i, iType) * fDsp + me%glb_sendbuf%tOffSv(4*(i-1)+iType)
        IF (nOpenReq .LT. nOpenReqMax) THEN
          nOpenReq = nOpenReq + 1
          iReq = nOpenReq
        ELSE
          CALL MPI_Waitany(nOpenReq, getReq, iReq, MPI_STATUS_IGNORE, ierr)
          getReq(iReq) = MPI_REQUEST_NULL
          IF (ierr /= MPI_SUCCESS) CALL finish(routine, "MPI error!")
        END IF
        SELECT CASE(iType)
        CASE(1)
          CALL MPI_RGET(output%d(1), 1, o_dType, idx%srcProc(i), t_dsp, &
            &           t_cnt(i), dType, me%glb_sendbuf%win, getReq(iReq), ierr)
        CASE(2)
          CALL MPI_RGET(output%s(1), 1, o_dType, idx%srcProc(i), t_dsp, &
            &           t_cnt(i), dType, me%glb_sendbuf%win, getReq(iReq), ierr)
        CASE(3)
          CALL MPI_RGET(output%i(1), 1, o_dType, idx%srcProc(i), t_dsp, &
            &           t_cnt(i), dType, me%glb_sendbuf%win, getReq(iReq), ierr)
        END SELECT
        IF (ierr /= MPI_SUCCESS) CALL finish(routine, "MPI error!")
        CALL MPI_Type_free(o_dType, ierr)
        rDsp(:) = rDsp(:) + pct(i,:)
      END IF
    END DO
    CALL MPI_Waitall(nOpenReq, getReq ,MPI_STATUSES_IGNORE, ierr)
    IF (ierr /= MPI_SUCCESS) CALL finish(routine, "MPI error!")
    used_size(:) = stride(:) * lCnt(:)
    srcOffset(:, iType) = srcOffset(:, iType) + INT(t_cnt(:), i8)
    DEALLOCATE(t_cnt, pct)
#else
    CALL finish(routine, "Not implemented!")
#endif
  END SUBROUTINE multifileRestartCollector_fetch

  SUBROUTINE multifileRestartCollector_sendField(me, lStart, lCnt, iV, input, ioffset)
    CLASS(t_MultifileRestartCollector), INTENT(INOUT) :: me
    INTEGER,          INTENT(IN   ) :: lStart, lCnt, iV
    TYPE(dataPtrs_t), INTENT(INOUT) :: input
    INTEGER(i8),      INTENT(INOUT) :: ioffset(3)
    CHARACTER(*), PARAMETER :: routine = modname//":multifileRestartCollector_sendField"
#ifndef NOMPI
    INTEGER     :: i, iLev, inType
    INTEGER(i8) :: offset
    TYPE(t_CollectorIndices), POINTER :: idx
    
    idx => me%idx(me%vGrid(iV))%p
    inType = me%vType(iV)
    IF (.NOT. ASSOCIATED(me%glb_sendbuf%sendBuffer%d)) THEN
      CALL finish(routine, "Unassociated send buffer!")
    END IF
    offset = ioffset(inType)
    IF (idx%sendPntCnt > 0) THEN
      DO iLev = lStart, lStart -1 + lCnt
        SELECT CASE(inType)
        CASE(1)
!$OMP PARALLEL DO SCHEDULE(STATIC)
          DO i = 1, idx%sendPntCnt
            me%glb_sendbuf%sendBuffer%d(offset+i) = &
              &  input%d(iLev)%p(idx%sendIdx(i), idx%sendBlk(i))
          END DO
        CASE(2)
!$OMP PARALLEL DO SCHEDULE(STATIC)
          DO i = 1, idx%sendPntCnt
            me%glb_sendbuf%sendBuffer%s(offset+i) = &
              &  input%s(iLev)%p(idx%sendIdx(i), idx%sendBlk(i))
          END DO
        CASE(3)
!$OMP PARALLEL DO SCHEDULE(STATIC)
          DO i = 1, idx%sendPntCnt
            me%glb_sendbuf%sendBuffer%i(offset+i) = &
              &  input%i(iLev)%p(idx%sendIdx(i), idx%sendBlk(i))
          END DO
        END SELECT
        offset = offset + INT(idx%sendPntCnt, i8)
      END DO
      SELECT CASE(inType)
      CASE(1)
        DEALLOCATE(input%d)
      CASE(2)
        DEALLOCATE(input%s)
      CASE(3)
        DEALLOCATE(input%i)
      END SELECT
      ioffset(inType) = offset
    END IF
#else
    CALL finish(routine, "Not implemented!")
#endif
  END SUBROUTINE multifileRestartCollector_sendField

  SUBROUTINE collectorSendBuffer_construct(this, isize, idx_c, idx_e, idx_v)
    CLASS(t_CollectorSendBuffer), INTENT(INOUT) :: this
    INTEGER(i8),                  INTENT(IN)    :: isize(:)
    TYPE(t_CollectorIndices),     INTENT(INOUT) :: idx_c, idx_e, idx_v
    CHARACTER(*), PARAMETER :: routine = modname//":CollectorSendBuffer_add_domain"
#ifdef NOMPI
    CALL finish(routine, "Not implemented!")
#else
    INTEGER, PARAMETER :: typeId(3) = (/ REAL_T, SINGLE_T, INT_T /)
    INTEGER :: ierr, splitKey, wcGrp, wcGrp_size, i, j, ii, iStart, myRank, &
      &        typeIdMPI(3), p_addr, addrBytes
    INTEGER, ALLOCATABLE, DIMENSION(:) :: rank_list, rank_map, tmpSrcProcs
    INTEGER, PARAMETER :: addr = MPI_ADDRESS_KIND
    INTEGER(KIND=addr), PARAMETER :: one = 1_addr, zero = 0_addr
    INTEGER(KIND=addr) :: typeBytes(3), typeLB
    INTEGER(KIND=addr), ALLOCATABLE :: tmpOffSv(:)

    typeIdMPI(1:3) = (/ p_real_dp, p_real_sp, p_int /)
    this%allocd      = .false.
    this%handshakd  = .false.
    this%wPosted  = .false.
    this%wStarted = .false.
    NULLIFY(this%wPtr%p)
    DO i = 1, 3
      this%wSizes(i) = isize(typeId(i))
      CALL MPI_Type_get_extent(typeIdMPI(i), typeLB, typeBytes(i), ierr)
      IF (ierr /= MPI_SUCCESS) CALL finish(routine, "MPI error!")
    END DO
    this%facDpSp = typeBytes(1) / typeBytes(2)
    IF (this%facDpSp .NE. (typeBytes(1)+typeBytes(2)-1)/typeBytes(2)) &
      CALL finish(routine, "sizeof(DP) not an integer multiple of sizeof(SP)!")
    this%facIntSp = typeBytes(3) / typeBytes(2)
    IF (this%facIntSp .NE. (typeBytes(2)+typeBytes(3)-1)/typeBytes(3)) &
      CALL finish(routine, "sizeof(DP) not an integer multiple of sizeof(INT)!")
    ALLOCATE(this%tOffCl(4))
    this%tOffCl(:) = zero
    IF (my_process_is_work()) THEN
      this%tOffCl(2) = this%tOffCl(1) + INT(this%wSizes(1), addr) * this%facDpSp
      this%tOffCl(3) = this%tOffCl(2) + INT(this%wSizes(2), addr)
      this%tOffCl(4) = this%tOffCl(3) + INT(this%wSizes(3), addr) * this%facIntSp - one
    END IF
    myRank = p_comm_rank(p_comm_work_restart)
    splitKey = MERGE(1, myRank + 2, my_process_is_restart_writer())
    CALL MPI_Comm_split(p_comm_work_restart, idx_c%destProc, splitKey, this%wComm, ierr)
    IF (ierr /= MPI_SUCCESS) CALL finish(routine, "MPI error!")
    idx_c%destProc = 0
    idx_e%destProc = 0
    idx_v%destProc = 0
    CALL MPI_Comm_group(this%wComm, wcGrp, ierr)
    CALL MPI_Group_size(wcGrp, wcGrp_size, ierr)
    ALLOCATE(rank_list(wcGrp_size))
    rank_list(1:wcGrp_size) = (/ (i, i = 0, wcGrp_size-1) /)
    CALL MPI_Group_incl(wcGrp, 1, rank_list(1:1), this%wSvGrp, ierr)
    IF (my_process_is_restart_writer()) THEN
      ALLOCATE(tmpOffSv(4*wcGrp_size), rank_map(wcGrp_size))
    ELSE
       ALLOCATE(tmpOffSv(1), rank_map(1))
    END IF
    CALL MPI_Gather(myRank, 1, p_int, rank_map, 1, p_int, 0, this%wComm, ierr)
    IF (ierr /= MPI_SUCCESS) CALL finish(routine, "MPI error!")
    CALL MPI_Sizeof(one, addrBytes, ierr)
    CALL MPI_Type_match_size(MPI_TYPECLASS_INTEGER, addrBytes, p_addr, ierr)
    CALL MPI_Gather(this%tOffCl, 4, p_addr, tmpOffSv, 4, p_addr, 0, this%wComm, ierr)
    IF (ierr /= MPI_SUCCESS) CALL finish(routine, "MPI error!")
    IF (my_process_is_restart_writer()) THEN
      iStart = MERGE(1, 2, my_process_is_work())
      CALL MPI_Group_incl(wcGrp, wcGrp_size-iStart+1, rank_list(iStart:wcGrp_size), &
        &                 this%wClGrp, ierr)
      ALLOCATE(tmpSrcProcs (idx_c%srcProcCnt), &
               this%tOffSv(4*idx_c%srcProcCnt))
      tmpSrcProcs(:)  = idx_c%srcProc(:)
      DO i = 1, idx_c%srcProcCnt
        DO j = 1, wcGrp_size
          IF (rank_map(j) .EQ. tmpSrcProcs(i)) ii = j - 1
        END DO
        idx_c%srcProc(i) = ii
        idx_e%srcProc(i) = ii
        idx_v%srcProc(i) = ii
        this%tOffSv((i-1)*4+1:i*4) = tmpOffSv(ii*4+1:(ii+1)*4)
      END DO
      DEALLOCATE(tmpSrcProcs)
    ELSE
      this%wClGrp = MPI_GROUP_NULL
    END IF
    CALL MPI_Group_free(wcGrp, ierr)
    DEALLOCATE(rank_list, rank_map, tmpOffSv)
    CALL this%handshake()
#endif
  END SUBROUTINE collectorSendBuffer_construct

  SUBROUTINE collectorSendBuffer_handshake(this)
    CLASS(t_CollectorSendBuffer), INTENT(INOUT) :: this
    CHARACTER(*), PARAMETER :: routine = modname//":CollectorSendBuffer_construct"
#ifdef NOMPI
    CALL finish(routine, "Not implemented!")
#else
    INTEGER, PARAMETER :: idummy(1) = (/ 1 /), addr = MPI_ADDRESS_KIND
    INTEGER(KIND=addr) :: memSize(1), memBytes, typeLB, spBytes
    INTEGER(KIND=addr), PARAMETER :: one = 1_addr
    TYPE(c_ptr) :: cMemPtr
    INTEGER :: ierr
    REAL(KIND=sp), POINTER :: tmp_sp(:)

    CALL MPI_Type_get_extent(p_real_sp, typeLB, spBytes, ierr)
    memSize(1) = this%tOffCl(4)
    memBytes = MAX(memSize(1) * spBytes, spBytes)
! as of MPI3.0 standard the following MPI_Alloc_mem interface must be 
! present, if the compiler provides ISO_C_BINDING
    CALL MPI_Alloc_mem(memBytes, MPI_INFO_NULL, cMemPtr, ierr)
    IF (ierr /= MPI_SUCCESS) CALL finish(routine, "MPI error!")
    CALL C_F_POINTER(cMemPtr, this%wPtr%p, idummy)
    CALL C_F_POINTER(cMemPtr, tmp_sp, INT(memSize))
    CALL MPI_Win_create(tmp_sp, memBytes, INT(spBytes), MPI_INFO_NULL, &
      &                 this%wComm, this%win, ierr)
    IF (ierr /= MPI_SUCCESS) CALL finish(routine, "MPI error!")
    memSize(1) = MAX(one, this%tOffCl(4))
    CALL C_F_POINTER(cMemPtr, tmp_sp, INT(memSize))
    cMemPtr = C_LOC(tmp_sp(this%tOffCl(1) + one))
    memSize(1) = INT(this%wSizes(1), addr)
    CALL C_F_POINTER(cMemPtr, this%sendBuffer%d,   INT(memSize))
    cMemPtr = C_LOC(tmp_sp(this%tOffCl(2) + one))
    memSize(1) = INT(this%wSizes(2), addr)
    CALL C_F_POINTER(cMemPtr, this%sendBuffer%s,   INT(memSize))
    cMemPtr = C_LOC(tmp_sp(this%tOffCl(3)+ one))
    memSize(1) = INT(this%wSizes(3), addr)
    CALL C_F_POINTER(cMemPtr, this%sendBuffer%i, INT(memSize))
    DEALLOCATE(this%tOffCl)
    CALL MPI_Comm_free(this%wComm, ierr)
    IF (ierr /= MPI_SUCCESS) CALL finish(routine, "MPI error!")
    this%wComm    = MPI_COMM_NULL
    this%handshakd = .true.
    this%allocd     = .true.
#endif
  END SUBROUTINE collectorSendBuffer_handshake

  SUBROUTINE CollectorSendBuffer_finalize(this)
    CLASS(t_CollectorSendBuffer), INTENT(INOUT) :: this
    CHARACTER(*), PARAMETER :: routine = modname//":CollectorSendBuffer_finalize"
#ifndef NOMPI
    INTEGER :: ierr

    IF (this%handshakd) THEN
      IF(this%wPosted .OR. this%wStarted) &
        CALL this%start_local_access()
      CALL MPI_Win_free(this%win, ierr)
      IF (ierr /= MPI_SUCCESS) CALL finish(routine, "MPI error!")
      this%handshakd  = .false.
    END IF
    IF (this%allocd) THEN
      CALL MPI_Free_mem(this%wPtr%p, ierr)
      IF (ierr /= MPI_SUCCESS) CALL finish(routine, "MPI error!")
      NULLIFY(this%wPtr%p, this%sendBuffer%d, this%sendBuffer%s, this%sendBuffer%i)
      this%wSizes(:) = 0_i8
      this%allocd = .false.
    END IF
    IF (this%wClGrp .NE. MPI_GROUP_NULL) THEN
      CALL MPI_Group_free(this%wClGrp, ierr)
    END IF
    CALL MPI_Group_free(this%wSvGrp, ierr)
    IF (this%wComm .NE. MPI_COMM_NULL) THEN
      CALL MPI_Comm_free(this%wComm, ierr)
      IF (ierr /= MPI_SUCCESS) CALL finish(routine, "MPI error!")
      this%wComm    = MPI_COMM_NULL
    END IF
    IF (ALLOCATED(this%tOffSv)) DEALLOCATE(this%tOffSv)
    IF (ALLOCATED(this%tOffCl)) DEALLOCATE(this%tOffCl)
    this%wPosted  = .false.
    this%wStarted = .false.
#endif
  END SUBROUTINE collectorSendBuffer_finalize

  SUBROUTINE collectorSendBuffer_start_local_access(this)
    CLASS(t_CollectorSendBuffer), INTENT(INOUT) :: this
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::CollectorSendBuffer_start_local_access"
#ifndef NOMPI
    INTEGER :: ierr

    IF (.NOT.this%allocd) CALL finish(routine, "there is no buffer allocd to fill!")
    IF (this%wStarted) THEN
      CALL MPI_Win_complete(this%win, ierr)
      IF (ierr /= MPI_SUCCESS) CALL finish(routine, "MPI error!")
      this%wStarted = .false.
    END IF
    IF (this%wPosted) THEN
      CALL MPI_Win_wait(this%win, ierr)
      IF (ierr /= MPI_SUCCESS) CALL finish(routine, "MPI error!")
      this%wPosted  = .false.
    ENDIF
#else
    CALL finish(routine, "Not implemented!")
#endif
  END SUBROUTINE collectorSendBuffer_start_local_access

  SUBROUTINE collectorSendBuffer_start_remote_access(this)
    CLASS(t_CollectorSendBuffer), INTENT(INOUT) :: this
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::CollectorSendBuffer_start_remote_access"
#ifndef NOMPI
    INTEGER :: iassert, ierr

    IF (.NOT.this%handshakd) CALL finish(routine, "there is no window to expose!")
    IF (my_process_is_work() .AND. .NOT.this%wPosted) THEN
      iassert = MPI_MODE_NOPUT
      CALL MPI_Win_post(this%wSvGrp, iassert, this%win, ierr)
      IF (ierr /= MPI_SUCCESS) CALL finish(routine, "MPI error!")
      this%wPosted  = .true.
    END IF
    IF (my_process_is_restart_writer() .AND. .NOT.this%wStarted) THEN
      iassert = 0
      CALL MPI_Win_start(this%wClGrp, iassert, this%win, ierr)
      IF (ierr /= MPI_SUCCESS) CALL finish(routine, "MPI error!")
      this%wStarted  = .true.
    END IF
#else
    CALL finish(routine, "Not implemented!")
#endif
  END SUBROUTINE collectorSendBuffer_start_remote_access

END MODULE mo_multifile_restart_collector
