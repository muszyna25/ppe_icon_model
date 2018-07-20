!> The multifile version of t_RestartPatchData.
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

MODULE mo_multifile_restart_patch_data
  USE mo_cdi,                         ONLY: FILETYPE_NC4, DATATYPE_FLT64, streamWriteVar, streamWriteVarSlice, &
    &                                       streamWriteVarSliceF
  USE mo_zaxis_type,                  ONLY: ZA_SURFACE
  USE mo_cdi_ids,                     ONLY: t_CdiIds
  USE mtime,                          ONLY: datetime
  USE mo_decomposition_tools,         ONLY: t_grid_domain_decomp_info
  USE mo_exception,                   ONLY: finish
  USE mo_fortran_tools,               ONLY: t_ptr_2d, t_ptr_2d_sp, t_ptr_2d_int, ensureSize,                   &
   &                                        no_copy, DO_PTR_DEALLOCATE
  USE mo_impl_constants,              ONLY: SUCCESS, SINGLE_T, REAL_T, INT_T
  USE mo_cdi_constants,               ONLY: GRID_UNSTRUCTURED_CELL, GRID_UNSTRUCTURED_EDGE,                    &
    &                                       GRID_UNSTRUCTURED_VERT
  USE mo_kind,                        ONLY: dp, sp, i8
  USE mo_model_domain,                ONLY: p_patch
  USE mo_mpi,                         ONLY: my_process_is_work
  USE mo_multifile_restart_collector, ONLY: t_MultifileRestartCollector, t_CollectorIndices, t_CollectorSendBuffer
  USE mo_multifile_restart_util,      ONLY: multifilePayloadPath, isAsync, iAmRestartWriter, commonBuf_t,      &
                                            dataPtrs_t, vNames_glbIdx
  USE mo_restart_attributes,          ONLY: t_RestartAttributeList
  USE mo_restart_descriptor,          ONLY: t_RestartPatchData
  USE mo_restart_patch_description,   ONLY: t_restart_patch_description
  USE mo_restart_util,                ONLY: t_restart_args 
  USE mo_restart_var_data,            ONLY: has_valid_time_level, getLevelPointers
  USE mo_timer,                       ONLY: timer_start, timer_stop, timer_write_restart_io,                   &
    &                                       timer_write_restart_communication, timers_level,                   &
    &                                       timer_write_restart_setup, timer_write_restart_wait
  USE mo_var_metadata_types,          ONLY: t_var_metadata
  USE mo_parallel_config,             ONLY: restart_chunk_size

  IMPLICIT NONE
  
  PUBLIC :: t_MultifilePatchData
  PUBLIC :: toMultifilePatchData

  PRIVATE

  TYPE :: ptr_arr_t
    REAL(KIND=dp), POINTER :: p(:)
  END TYPE ptr_arr_t

  TYPE, EXTENDS(t_RestartPatchData) :: t_MultifilePatchData
    PRIVATE
    TYPE(commonBuf_t) :: RecvBuffer
    TYPE(ptr_arr_t) :: glb_idx(3)
    INTEGER,  POINTER :: varReordered(:)
    TYPE (t_CollectorSendBuffer) :: glb_sendbuf
    TYPE(t_CollectorIndices) :: cellColl, edgeColl, vertColl
    TYPE(t_MultifileRestartCollector) :: collector
  CONTAINS
    PROCEDURE :: construct => multifilePatchData_construct
    PROCEDURE :: createCollectors => multifilePatchData_createCollectors
    PROCEDURE :: openPayloadFile => multifilePatchData_openPayloadFile  ! restart writers only
    PROCEDURE :: start_local_access  => multifilePatchData_start_local_access
    PROCEDURE :: start_remote_access => multifilePatchData_start_remote_access
    PROCEDURE :: exposeData  => multifilePatchData_exposeData    
    PROCEDURE :: collectData => multifilePatchData_collectData
    PROCEDURE :: writeFile => multifilePatchData_writeFile  ! override
    PROCEDURE :: destruct => multifilePatchData_destruct    ! override
  END TYPE t_MultifilePatchData
  
  CHARACTER(*), PARAMETER :: modname = "mo_multifile_restart_patch_data"
  
CONTAINS

  FUNCTION toMultifilePatchData(me) RESULT(resultVar)
    CLASS(t_RestartPatchData), TARGET, INTENT(INOUT) :: me
    TYPE(t_MultifilePatchData), POINTER :: resultVar
    CHARACTER(*), PARAMETER :: routine = modname//":toMultifilePatchData"

    SELECT TYPE(me)
    TYPE IS(t_MultifilePatchData)
      resultVar => me
    CLASS DEFAULT
      CALL finish(routine, "assertion failed: wrong dynamic type")
    END SELECT
  END FUNCTION toMultifilePatchData


  SUBROUTINE multifilePatchData_construct(me, modelType, domain)
    CLASS(t_MultifilePatchData), INTENT(INOUT) :: me
    CHARACTER(*),                INTENT(IN)    :: modelType
    INTEGER,                     INTENT(IN)    :: domain
    CHARACTER(*), PARAMETER :: routine = modname//":multifilePatchData_construct"

    NULLIFY(me%RecvBuffer%d, me%RecvBuffer%s, me%RecvBuffer%i, &
     &      me%glb_idx(1)%p, me%glb_idx(2)%p, me%glb_idx(3)%p)
    CALL me%t_RestartPatchData%construct(modelType, domain)
    IF (isAsync()) CALL me%transferToRestart()
  END SUBROUTINE multifilePatchData_construct


  ! Note: This should have been part of the constructor, however,
  ! Fortran does not allow changing the signature while overriding a
  ! method, nor does it allow overloading the overridden method, so we
  ! need to use a different name for this
  SUBROUTINE multifilePatchData_createCollectors(me, writerRank, sourceRanks, lthis_pe_active)
    CLASS(t_MultifilePatchData), INTENT(INOUT), TARGET :: me
    INTEGER,                     INTENT(IN) :: writerRank, sourceRanks(:)
    LOGICAL,                     INTENT(IN) :: lthis_pe_active
    CHARACTER(*), PARAMETER :: routine = modname//":multifilePatchData_createCollectors"
    INTEGER, PARAMETER :: TYPE_LEN = MAXVAL([REAL_T, SINGLE_T, INT_T])
    INTEGER                         :: curVar, jg, nLevs, iType, nVar, i
    INTEGER(i8)                     :: iOffset(TYPE_LEN)
    TYPE(t_var_metadata), POINTER   :: curInfo
    TYPE(t_grid_domain_decomp_info) :: dummyInfo

    IF (timers_level >= 7) CALL timer_start(timer_write_restart_setup)
    jg = me%description%id
    ! create a global receive buffer (for joint proc mode this is the
    ! only receive buffer which is actually allocated):
    ! dedicated restart PEs MUST NOT access p_patch! (this obviously was NEVER used!)
    IF (my_process_is_work()) THEN
      me%glb_idx(1)%p => me%cellColl%construct(p_patch(jg)%n_patch_cells, &
        &                                      p_patch(jg)%cells%decomp_info, &
        &                                      writerRank, sourceRanks, lthis_pe_active)
      me%glb_idx(2)%p => me%edgeColl%construct(p_patch(jg)%n_patch_edges, &
        &                                      p_patch(jg)%edges%decomp_info, &
        &                                      writerRank, sourceRanks, lthis_pe_active)
      me%glb_idx(3)%p => me%vertColl%construct(p_patch(jg)%n_patch_verts, &
        &                                      p_patch(jg)%verts%decomp_info, &
        &                                      writerRank, sourceRanks, lthis_pe_active)
    ELSE
      me%glb_idx(1)%p => me%cellColl%construct(0, dummyInfo, writerRank, &
        &                                      sourceRanks, lthis_pe_active)
      me%glb_idx(2)%p => me%edgeColl%construct(0, dummyInfo, writerRank, &
        &                                      sourceRanks, lthis_pe_active)
      me%glb_idx(3)%p => me%vertColl%construct(0, dummyInfo, writerRank, &
        &                                      sourceRanks, lthis_pe_active)
    END IF
    nVar = SIZE(me%varData)
    ioffset(:) = 0_i8
    ALLOCATE(me%varReordered(nVar))
    i = 0
    DO iType = 1, TYPE_LEN   
      DO curVar = 1, nVar
        curInfo => me%varData(curVar)%info
        IF (curInfo%data_type .EQ. iType) THEN
          i = i + 1
          me%varReordered(i) = curVar
        END IF
      END DO
    END DO
    IF (i .NE. nVar) CALL finish(routine, "inconsistentcy!!!")
    CALL me%collector%construct(me%cellColl, me%edgeColl, me%vertColl, me%glb_sendbuf, nVar)
    DO curVar = 1, nVar
      curInfo => me%varData(me%varReordered(curVar))%info
      iType = curInfo%data_type
      !get the effective level count
      nLevs = 1
      IF (curInfo%ndims == 3) nLevs = curInfo%used_dimensions(2)
      IF (curInfo%ndims > 3) CALL finish(routine, "ndims > 3 is not supported")
      ! get the correct collector and buffer for the first level
      CALL me%collector%defVar(curVar, nLevs, iType, curInfo%hgrid, iOffset) !, sourceOffset)
    END DO
    ! finally, allocate the global send buffers:
    CALL me%glb_sendbuf%construct(iOffset, me%cellColl, me%edgeColl, me%vertColl)
    IF (timers_level >= 7) CALL timer_stop(timer_write_restart_setup)
  END SUBROUTINE multifilePatchData_createCollectors

  INTEGER FUNCTION defineIndexVar(cdiIds, hgrid, name_str) RESULT(varId)
    TYPE(t_CdiIds), INTENT(INOUT) :: cdiIds
    INTEGER,        INTENT(IN)    :: hgrid
    CHARACTER(*),   INTENT(IN)    :: name_str

    varId = cdiIds%defineVariableSimple(hgrid, ZA_SURFACE, DATATYPE_FLT64, name_str)
  END FUNCTION defineIndexVar

  FUNCTION multifilePatchData_openPayloadFile(me, filename, ifile, this_datetime, bytesWritten) RESULT(cdiIds)
    TYPE(t_CdiIds)                                     :: cdiIds
    CLASS(t_MultifilePatchData), TARGET, INTENT(INOUT) :: me
    CHARACTER(*),                        INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ifile
    TYPE(datetime), POINTER,             INTENT(IN)    :: this_datetime
    INTEGER(i8),                         INTENT(INOUT) :: bytesWritten
    CHARACTER(*), PARAMETER                            :: routine = modname//":multifilePatchData_openPayloadFile"
    CHARACTER(:), ALLOCATABLE                          :: effectiveFilename
    CLASS(t_restart_patch_description), POINTER        :: desc
    INTEGER                                            :: iVarId(3), curVar, iG
    REAL(dp)                                           :: dummy_d(1)

    dummy_d(1) = 0._dp
    IF (timers_level >= 7) CALL timer_start(timer_write_restart_io)
    desc => me%description
    CALL multifilePayloadPath(filename, desc%id, ifile, effectiveFilename)
    CALL cdiIds%init()
    IF (ALLOCATED(desc%opt_pvct)) THEN
      CALL cdiIds%openRestartAndCreateIds(effectiveFilename, FILETYPE_NC4,              &
        &                                 SIZE(me%glb_idx(1)%p), SIZE(me%glb_idx(3)%p), &
        &                                 SIZE(me%glb_idx(2)%p), desc%cell_type,        &
        &                                 desc%v_grid_defs(1:desc%v_grid_count),        &
        &                                 desc%opt_pvct)
    ELSE
      CALL cdiIds%openRestartAndCreateIds(effectiveFilename, FILETYPE_NC4,              &
        &                                 SIZE(me%glb_idx(1)%p), SIZE(me%glb_idx(3)%p), &
        &                                 SIZE(me%glb_idx(2)%p), desc%cell_type,        &
        &                                 desc%v_grid_defs(1:desc%v_grid_count))
    END IF
    iVarId(1) = defineIndexVar(cdiIds, GRID_UNSTRUCTURED_CELL, vNames_glbIdx(1))
    iVarId(2) = defineIndexVar(cdiIds, GRID_UNSTRUCTURED_EDGE, vNames_glbIdx(2))
    iVarId(3) = defineIndexVar(cdiIds, GRID_UNSTRUCTURED_VERT, vNames_glbIdx(3))
    ! define all restart variables with a valid time level
    DO curVar = 1, SIZE(me%varData)
      IF (has_valid_time_level(me%varData(curVar)%info, desc%id, desc%nnew, &
        &                     desc%nnew_rcf)) THEN
        CALL cdiIds%defineVariable(me%varData(curVar)%info)
      END IF
    END DO
    CALL cdiIds%finalizeVlist(this_datetime)
    ! write the arrays with the global indices of the points
    DO iG = 1, 3
      IF (SIZE(me%glb_idx(iG)%p) > 0) THEN
        CALL streamWriteVar(cdiIds%fHndl, iVarId(iG), me%glb_idx(iG)%p, 0)
        bytesWritten = bytesWritten + INT(SIZE(me%glb_idx(iG)%p), i8)*8_i8
      ELSE
        ! even zero-sized payload file must exist...
        CALL streamWriteVar(cdiIds%fHndl, iVarId(iG), dummy_d, 0)
      END IF
    END DO
    IF (timers_level >= 7) CALL timer_stop(timer_write_restart_io)
  END FUNCTION multifilePatchData_openPayloadFile

  SUBROUTINE multifilePatchData_start_local_access(me)
    CLASS(t_MultifilePatchData), INTENT(INOUT) :: me

    IF (timers_level >= 7) CALL timer_start(timer_write_restart_wait)
    CALL me%glb_sendbuf%start_local_access()
    IF (timers_level >= 7) CALL timer_stop(timer_write_restart_wait)
  END SUBROUTINE multifilePatchData_start_local_access

  ! Start MPI epoch where windows may not be accessed locally and
  ! not via MPI_PUT operations.
  SUBROUTINE multifilePatchData_start_remote_access(me)
    CLASS(t_MultifilePatchData), INTENT(INOUT) :: me

    IF (timers_level >= 7) CALL timer_start(timer_write_restart_wait)
    CALL me%glb_sendbuf%start_remote_access()
    IF (timers_level >= 7) CALL timer_stop(timer_write_restart_wait)
  END SUBROUTINE multifilePatchData_start_remote_access

  ! Collect the data, optionally writing it immediately.
  !
  SUBROUTINE multifilePatchData_exposeData(me)
    CLASS(t_MultifilePatchData), INTENT(INOUT) :: me
    CHARACTER(*), PARAMETER :: routine = modname//":multifilePatchData_exposeData"
    TYPE(t_var_metadata), POINTER   :: curInfo
    TYPE(dataPtrs_t)                :: dataPointers
    INTEGER                         :: startLevel, nLevel, curVar, curVar_o
    INTEGER(KIND=i8)                :: ioffset(3)

    IF (.NOT. my_process_is_work()) THEN
      CALL finish(routine, "assertion failed.")
    END IF
    ioffset(:) = 0_i8
    IF (timers_level >= 7) CALL timer_start(timer_write_restart_io)
    DO curVar_o = 1, SIZE(me%varData)
      curVar = me%varReordered(curVar_o)
      curInfo => me%varData(curVar)%info
      IF (.NOT.has_valid_time_level(curInfo, me%description%id, me%description%nnew, &
        &                          me%description%nnew_rcf)) CYCLE
      startLevel = 1
      nLevel = 1
      IF (curInfo%ndims == 3) nLevel = curInfo%used_dimensions(2)
      CALL getLevelPointers(curInfo, me%varData(curVar), dataPointers)
      CALL me%collector%sendField(startLevel, nLevel, curVar_o, dataPointers, ioffset)
    END DO
    CALL dataPointers%free()
    IF (timers_level >= 7) CALL timer_stop(timer_write_restart_io)
  END SUBROUTINE multifilePatchData_exposeData

  ! Collect the data, optionally writing it immediately.
  !
  SUBROUTINE multifilePatchData_collectData(me, fileId, bytesWritten)
    CLASS(t_MultifilePatchData), INTENT(INOUT) :: me
    INTEGER,                     INTENT(IN)    :: fileId
    INTEGER(i8),                 INTENT(INOUT) :: bytesWritten
    CHARACTER(*), PARAMETER :: routine = modname//":multifilePatchData_collectData"
    TYPE(t_var_metadata), POINTER   :: curInfo
    REAL(KIND=dp), POINTER, DIMENSION(:) :: ptr2level_dp
    REAL(KIND=sp), POINTER, DIMENSION(:) :: ptr2level_sp
    INTEGER, ALLOCATABLE, DIMENSION(:) :: used_size, lCnt, lStart
    LOGICAL, ALLOCATABLE, DIMENSION(:) :: vSkip
    INTEGER :: lFree, vNext, lNext, vCnt, vFrst, lOff, vOff, lInc, &
               lFrst, tFrst, i, iV, iL, lN, iG
    INTEGER(KIND=i8), ALLOCATABLE :: srcOff(:,:)

    DO iG = 1, 3
      IF (.NOT.ASSOCIATED(me%glb_idx(iG)%p)) THEN
        RETURN
      ELSE
        IF (SIZE(me%glb_idx(iG)%p) .LE. 0) RETURN
      END IF
    END DO
    NULLIFY(ptr2level_dp, ptr2level_sp)
    vNext = 1
    lNext = 1
    DO WHILE (vNext .LE. SIZE(me%varData))
      lFree = restart_chunk_size
      vCnt = 0
      vFrst = vNext
      lFrst = lNext
      DO iV = vFrst, SIZE(me%varData)
        IF (lFree .EQ. 0) CYCLE
        curInfo => me%varData(me%varReordered(iV))%info
        IF (.NOT.has_valid_time_level(curInfo, me%description%id, &
          &        me%description%nnew, me%description%nnew_rcf)) THEN
          IF (vFrst .EQ. iV) THEN
            vFrst = iV + 1
          ELSE
            vCnt = vCnt + 1
          END IF
          vNext = iV + 1
          CYCLE
        END IF
        IF (vFrst .EQ. iV) tFrst = curInfo%data_type
        IF (tFrst .NE. curInfo%data_type) lFree = 0
        vCnt = vCnt + 1
        lN = 1
        IF (curInfo%ndims == 3) lN = curInfo%used_dimensions(2)
        IF (lFree .GE. lN - lNext + 1) THEN
          lFree = lFree - (lN - lNext + 1)
          lNext = 1
          vNext = iV + 1
        ELSE
          lNext = lNext + lFree
          lFree = 0
        END IF
      END DO
      IF (vCnt .LE. 0) CYCLE
      ALLOCATE(used_size(vCnt), lCnt(vCnt), lStart(vCnt), vSkip(vCnt))
      used_size(:) = 0
      lStart(1) = lFrst
      lCnt(:) = 1
      IF (vCnt .GT. 1) lStart(2:vCnt) = 1
      DO iV = 1, vCnt
        curInfo => me%varData(me%varReordered(iV + vFrst - 1))%info
        vSkip(iV) = .NOT.has_valid_time_level(curInfo, me%description%id, &
          &                me%description%nnew, me%description%nnew_rcf)
        IF (curInfo%ndims == 3) lCnt(iV) = curInfo%used_dimensions(2)
        IF (iV + vFrst - 1 .EQ. vNext) lCnt(iV) = lNext - 1
        lCnt(iV) = lCnt(iV) - lStart(iV) + 1
      END DO
      IF (timers_level >= 7) CALL timer_start(timer_write_restart_communication)
      CALL me%collector%fetch(vFrst, vCnt, vSkip, lCnt, used_size, me%RecvBuffer, srcOff)
      IF (timers_level >= 7) CALL timer_stop(timer_write_restart_communication)
      IF (timers_level >= 7) CALL timer_start(timer_write_restart_io)
      vOff = 0
      DO iV = 1, vCnt
        IF (vSkip(iV)) CYCLE
        curInfo => me%varData(me%varReordered(iV + vFrst - 1))%info
        lOff = 0
        lInc = used_size(iV)/MAX(lCnt(iV), 1)
        DO iL = 1, lCnt(iV)
          SELECT CASE(tFrst)
          CASE(REAL_T)
            ptr2level_dp => me%recvBuffer%d(1+vOff+lOff:vOff+lOff+lInc)
            CALL streamWriteVarSlice(fileId, curInfo%cdiVarId, iL+lStart(iV)-2, ptr2level_dp, 0)
            bytesWritten = bytesWritten + INT(lInc, i8) * 8_i8
          CASE(SINGLE_T)
            ptr2level_sp => me%recvBuffer%s(1+vOff+lOff:vOff+lOff+lInc)
            CALL streamWriteVarSliceF(fileId, curInfo%cdiVarId, iL+lStart(iV)-2, ptr2level_sp, 0)
            bytesWritten = bytesWritten + INT(lInc, i8) * 4_i8
          CASE(INT_T)
            CALL ensureSize(me%RecvBuffer%d, lInc, no_copy)
            !$OMP PARALLEL DO SCHEDULE(STATIC)
            DO i = 1, lInc
              me%RecvBuffer%d(i) = REAL(me%RecvBuffer%i(vOff+lOff+i), dp)
            END DO
            ptr2level_dp => me%recvBuffer%d(1:lInc)
            CALL streamWriteVarSlice(fileId, curInfo%cdiVarId, iL+lStart(iV)-2, ptr2level_dp, 0)
            bytesWritten = bytesWritten + INT(lInc, i8) * 8_i8
          CASE DEFAULT
            CALL finish(routine, "Internal error! Variable "//TRIM(curInfo%name))
          END SELECT
          lOff = lOff + lInc
        END DO
        vOff = vOff + used_size(iV)
      END DO
      IF (timers_level >= 7) CALL timer_stop(timer_write_restart_io)
      DEALLOCATE(used_size, lCnt, lStart, vSkip)
    END DO
    IF(ALLOCATED(srcOff)) DEALLOCATE(srcOff)
  END SUBROUTINE multifilePatchData_collectData

  SUBROUTINE multifilePatchData_writeFile(me, restartAttributes, restartArgs, lIsWriteProcess)
    CLASS(t_MultifilePatchData),  INTENT(INOUT) :: me
    TYPE(t_RestartAttributeList), INTENT(INOUT) :: restartAttributes
    TYPE(t_restart_args),         INTENT(IN)    :: restartArgs
    LOGICAL, VALUE :: lIsWriteProcess
    CHARACTER(*), PARAMETER :: routine = modname//":multifilePatchData_writeFile"

    CALL finish(routine, "assertion failed: this implementation must not be called")
  END SUBROUTINE multifilePatchData_writeFile

  SUBROUTINE multifilePatchData_destruct(me)
    CLASS(t_MultifilePatchData), INTENT(INOUT) :: me
    CHARACTER(*), PARAMETER :: routine = modname//":multifilePatchData_destruct"

    CALL me%cellColl%finalize()
    CALL me%edgeColl%finalize()
    CALL me%vertColl%finalize()
    CALL me%collector%finalize()
    CALL me%t_RestartPatchData%destruct()
    CALL me%glb_sendbuf%finalize()
    CALL DO_PTR_DEALLOCATE(me%RecvBuffer%d)
    CALL DO_PTR_DEALLOCATE(me%RecvBuffer%s)
    CALL DO_PTR_DEALLOCATE(me%RecvBuffer%i)
    CALL DO_PTR_DEALLOCATE(me%glb_idx(1)%p)
    CALL DO_PTR_DEALLOCATE(me%glb_idx(2)%p)
    CALL DO_PTR_DEALLOCATE(me%glb_idx(3)%p)
    CALL DO_PTR_DEALLOCATE(me%varReordered)
  END SUBROUTINE multifilePatchData_destruct

END MODULE mo_multifile_restart_patch_data
