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
  USE mo_fortran_tools,               ONLY: t_ptr_2d, t_ptr_2d_sp, t_ptr_2d_int, ensureSize, no_copy
  USE mo_impl_constants,              ONLY: SUCCESS, SINGLE_T, REAL_T, INT_T
  USE mo_cdi_constants,               ONLY: GRID_UNSTRUCTURED_CELL, GRID_UNSTRUCTURED_EDGE,                    &
    &                                       GRID_UNSTRUCTURED_VERT
  USE mo_kind,                        ONLY: dp, sp, i8
  USE mo_model_domain,                ONLY: p_patch
  USE mo_mpi,                         ONLY: my_process_is_work
  USE mo_multifile_restart_collector, ONLY: t_MultifileRestartCollector, commonBuf_t, dataPtrs_t,              &
    &                                       t_CollectorIndices, t_CollectorSendBuffer
  USE mo_multifile_restart_util,      ONLY: multifilePayloadPath, kVarName_globalCellIndex,                    &
    &                                       kVarName_globalEdgeIndex, kVarName_globalVertIndex
  USE mo_restart_attributes,          ONLY: t_RestartAttributeList
  USE mo_restart_descriptor,          ONLY: t_RestartPatchData
  USE mo_restart_patch_description,   ONLY: t_restart_patch_description
  USE mo_restart_util,                ONLY: my_process_is_restart_writer,                                      &
    &                                       isDedicatedProcMode, t_restart_args 
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

  TYPE, EXTENDS(t_RestartPatchData) :: t_MultifilePatchData
    PRIVATE
    TYPE(commonBuf_t) :: RecvBuffer
    REAL(dp), POINTER :: glb_indices_cell(:)
    REAL(dp), POINTER :: glb_indices_edge(:)
    REAL(dp), POINTER :: glb_indices_vert(:)
    INTEGER,  POINTER :: varReordered(:)
    TYPE (t_CollectorSendBuffer) :: glb_sendbuf
    TYPE(t_CollectorIndices) :: cellCollector, edgeCollector, vertCollector
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
     &      me%glb_indices_cell, me%glb_indices_edge, me%glb_indices_vert)
    
    CALL me%t_RestartPatchData%construct(modelType, domain)
    IF (isDedicatedProcMode()) CALL me%transferToRestart()
  END SUBROUTINE multifilePatchData_construct


  ! Note: This should have been part of the constructor, however,
  ! Fortran does not allow changing the signature while overriding a
  ! method, nor does it allow overloading the overridden method, so we
  ! need to use a different name for this
  SUBROUTINE multifilePatchData_createCollectors(me, writerRank, sourceRanks, lthis_pe_active)
    CLASS(t_MultifilePatchData), INTENT(INOUT), TARGET :: me
    INTEGER,                     INTENT(IN) :: writerRank
    INTEGER,                     INTENT(IN) :: sourceRanks(:)
    ! Flag. .FALSE. if this PE does not send any points (since it will
    ! write to another file).
    LOGICAL,                     INTENT(IN) :: lthis_pe_active
    CHARACTER(*), PARAMETER :: routine = modname//":multifilePatchData_createCollectors"
    INTEGER, PARAMETER :: TYPE_LEN = MAXVAL([REAL_T, SINGLE_T, INT_T])
    INTEGER                         :: curVar, ierr, jg, nLevs, iType, nVar, i
    INTEGER(i8)                     :: iOffset(TYPE_LEN)
    TYPE(t_var_metadata), POINTER   :: curInfo
    INTEGER, ALLOCATABLE            :: sourceOffset(:,:)
    TYPE(t_grid_domain_decomp_info) :: dummyInfo

    IF (timers_level >= 7) CALL timer_start(timer_write_restart_setup)
    jg = me%description%id
    ! create a global receive buffer (for joint proc mode this is the
    ! only receive buffer which is actually allocated):
    ! dedicated restart PEs MUST NOT access p_patch! (this obviously was NEVER used!)
    IF (my_process_is_work()) THEN
      me%glb_indices_cell => me%cellCollector%construct(p_patch(jg)%n_patch_cells, &
        &                                     p_patch(jg)%cells%decomp_info, &
        &                                     writerRank, sourceRanks, lthis_pe_active)
      me%glb_indices_edge => me%edgeCollector%construct(p_patch(jg)%n_patch_edges, &
        &                                     p_patch(jg)%edges%decomp_info, &
        &                                     writerRank, sourceRanks, lthis_pe_active)
      me%glb_indices_vert => me%vertCollector%construct(p_patch(jg)%n_patch_verts, &
        &                                     p_patch(jg)%verts%decomp_info, &
        &                                     writerRank, sourceRanks, lthis_pe_active)
    ELSE
      me%glb_indices_cell => me%cellCollector%construct(0, dummyInfo, writerRank, &
        &                                     sourceRanks, lthis_pe_active)
      me%glb_indices_edge => me%edgeCollector%construct(0, dummyInfo, writerRank, &
        &                                     sourceRanks, lthis_pe_active)
      me%glb_indices_vert => me%vertCollector%construct(0, dummyInfo, writerRank, &
        &                                     sourceRanks, lthis_pe_active)
    END IF
    nVar = SIZE(me%varData)
    IF (my_process_is_restart_writer()) THEN
      ALLOCATE(sourceOffset(TYPE_LEN, me%cellCollector%sourceProcCount), STAT=ierr)
      IF ( &
        & me%cellCollector%sourceProcCount .NE. me%edgeCollector%sourceProcCount .OR. &
        & me%cellCollector%sourceProcCount .NE. me%vertCollector%sourceProcCount .OR. &
        & me%edgeCollector%sourceProcCount .NE. me%vertCollector%sourceProcCount) &
        & CALL finish(routine, "inconsistent!!!")
    ELSE
      ALLOCATE(sourceOffset(TYPE_LEN, 1), STAT=ierr)
    END IF
    IF (ierr /= SUCCESS) CALL finish(routine, "memory allocation failure")
    sourceOffset(:,:) = 0
    ioffset(:) = 0
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
    CALL me%collector%construct(me%cellCollector, me%edgeCollector, me%vertCollector, &
                                me%glb_sendbuf, nVar)
    DO curVar = 1, nVar
      curInfo => me%varData(me%varReordered(curVar))%info
      iType = curInfo%data_type
      !get the effective level count
      nLevs = 1
      IF (curInfo%ndims == 3) nLevs = curInfo%used_dimensions(2)
      IF (curInfo%ndims > 3) CALL finish(routine, "ndims > 3 is not supported")
      ! get the correct collector and buffer for the first level
      CALL me%collector%defVar(curVar, nLevs, iType, curInfo%hgrid, iOffset, sourceOffset)
    END DO
    ! finally, allocate the global send buffers:
    CALL me%glb_sendbuf%construct(ioffset, me%cellCollector, me%edgeCollector, me%vertCollector)
    ! clean up
    DEALLOCATE(sourceOffset, STAT=ierr)
    IF (ierr /= SUCCESS) CALL finish(routine, "memory deallocation failure")
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
    CLASS(t_restart_patch_description), POINTER        :: description
    INTEGER                                            :: cellIndexVarId, edgeIndexVarId, vertIndexVarId, curVar
    REAL(dp)                                           :: dummy_d(1)

    dummy_d(1) = 0._dp
    IF (timers_level >= 7) CALL timer_start(timer_write_restart_io)
    description => me%description
    CALL multifilePayloadPath(filename, description%id, ifile, effectiveFilename)
    CALL cdiIds%init()
    IF (ALLOCATED(description%opt_pvct)) THEN
      CALL cdiIds%openRestartAndCreateIds(effectiveFilename, FILETYPE_NC4,                     &
        &                                 SIZE(me%glb_indices_cell),                           &
        &                                 SIZE(me%glb_indices_vert),                           &
        &                                 SIZE(me%glb_indices_edge),                           &
        &                                 description%cell_type,                               &
        &                                 description%v_grid_defs(1:description%v_grid_count), &
        &                                 description%opt_pvct)
    ELSE
      CALL cdiIds%openRestartAndCreateIds(effectiveFilename, FILETYPE_NC4,                     &
        &                                 SIZE(me%glb_indices_cell),                           &
        &                                 SIZE(me%glb_indices_vert),                           &
        &                                 SIZE(me%glb_indices_edge),                           &
        &                                 description%cell_type,                               &
        &                                 description%v_grid_defs(1:description%v_grid_count))
    END IF
    cellIndexVarId = defineIndexVar(cdiIds, GRID_UNSTRUCTURED_CELL, kVarName_globalCellIndex)
    edgeIndexVarId = defineIndexVar(cdiIds, GRID_UNSTRUCTURED_EDGE, kVarName_globalEdgeIndex)
    vertIndexVarId = defineIndexVar(cdiIds, GRID_UNSTRUCTURED_VERT, kVarName_globalVertIndex)
    ! define all restart variables with a valid time level
    DO curVar = 1, SIZE(me%varData)
      IF (has_valid_time_level(me%varData(curVar)%info, me%description%id, description%nnew, &
        &                     description%nnew_rcf)) THEN
        CALL cdiIds%defineVariable(me%varData(curVar)%info)
      END IF
    END DO
    CALL cdiIds%finalizeVlist(this_datetime)
    ! write the arrays with the global indices of the points
    IF (SIZE(me%glb_indices_cell) > 0) THEN
      CALL streamWriteVar(cdiIds%file, cellIndexVarId, me%glb_indices_cell, 0)
      bytesWritten = bytesWritten + INT(SIZE(me%glb_indices_cell), i8)*8_i8
    ELSE
      ! Note that we may have a cell count of 0 for this payload file;
      ! this may happen in the case of processor splitting. Still, we
      ! have to write *something*, otherwise the CDI will not store
      ! NetCDF attributes, which, in turn, is necessary to detect
      ! size-0 payload files!
      CALL streamWriteVar(cdiIds%file, cellIndexVarId, dummy_d, 0)
    END IF
    IF (SIZE(me%glb_indices_edge) > 0) THEN
      CALL streamWriteVar(cdiIds%file, edgeIndexVarId, me%glb_indices_edge, 0)
      bytesWritten = bytesWritten + INT(SIZE(me%glb_indices_edge), i8)*8_i8
    ELSE
      CALL streamWriteVar(cdiIds%file, edgeIndexVarId, dummy_d, 0)
    END IF
    IF (SIZE(me%glb_indices_vert) > 0) THEN
      CALL streamWriteVar(cdiIds%file, vertIndexVarId, me%glb_indices_vert, 0)
      bytesWritten = bytesWritten + INT(SIZE(me%glb_indices_vert), i8)*8_i8
    ELSE
      CALL streamWriteVar(cdiIds%file, vertIndexVarId, dummy_d, 0)
    END IF
    IF (timers_level >= 7) CALL timer_stop(timer_write_restart_io)
  END FUNCTION multifilePatchData_openPayloadFile

  SUBROUTINE multifilePatchData_start_local_access(me)
    CLASS(t_MultifilePatchData), INTENT(INOUT) :: me

    IF (timers_level >= 7) CALL timer_start(timer_write_restart_communication)
    CALL me%glb_sendbuf%start_local_access()
    IF (timers_level >= 7) CALL timer_stop(timer_write_restart_communication)
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

    IF (.NOT. my_process_is_work()) THEN
      CALL finish(routine, "assertion failed.")
    END IF
    IF (timers_level >= 7) CALL timer_start(timer_write_restart_io)
    LOOP_VAR1 : DO curVar_o = 1, SIZE(me%varData)
      curVar = me%varReordered(curVar_o)
      curInfo => me%varData(curVar)%info
      ! no valid time level -> no output:
      IF (.NOT.has_valid_time_level(curInfo, me%description%id, me%description%nnew, &
        &                          me%description%nnew_rcf)) CYCLE LOOP_VAR1
      startLevel = 1
      nLevel   = 1
      IF (curInfo%ndims == 3) nLevel = curInfo%used_dimensions(2)
      SELECT CASE(curInfo%data_type)
      CASE(REAL_T)
        CALL getLevelPointers(curInfo, me%varData(curVar)%r_ptr, dataPointers%d)
      CASE(SINGLE_T)
        CALL getLevelPointers(curInfo, me%varData(curVar)%s_ptr, dataPointers%s)
      CASE(INT_T)
        CALL getLevelPointers(curInfo, me%varData(curVar)%i_ptr, dataPointers%i)
      CASE DEFAULT
        CALL finish(routine, "Internal error! Variable "//TRIM(curInfo%name))
      END SELECT
      CALL me%collector%sendField(startLevel, nLevel, curVar_o, dataPointers)
    END DO LOOP_VAR1
    IF (timers_level >= 7) CALL timer_stop(timer_write_restart_io)
  END SUBROUTINE multifilePatchData_exposeData

  ! Collect the data, optionally writing it immediately.
  !
  SUBROUTINE multifilePatchData_collectData(me, fileId, bytesWritten)
    CLASS(t_MultifilePatchData), INTENT(INOUT) :: me
    INTEGER,                     INTENT(IN)    :: fileId
    INTEGER(i8),                 INTENT(INOUT) :: bytesWritten
    CHARACTER(*), PARAMETER :: routine = modname//":multifilePatchData_collectData"
    INTEGER                         :: curVar
    TYPE(t_var_metadata), POINTER   :: curInfo
    INTEGER                         :: startLevel, endLevel, used_size, ilevel
    INTEGER                         :: nchunks, ichunk, nlevels, ilevel_start, i
    REAL(KIND=dp), POINTER, DIMENSION(:) :: ptr2level_dp
    REAL(KIND=sp), POINTER, DIMENSION(:) :: ptr2level_sp

    IF (.NOT.ASSOCIATED(me%glb_indices_cell)) THEN !checking one of cell/edge/vert should suffice
      RETURN
    ELSE
      IF (SIZE(me%glb_indices_cell) .LE. 0) &
        RETURN
    END IF
    NULLIFY(ptr2level_dp, ptr2level_sp)
    LOOP_VAR2 : DO curVar = 1, SIZE(me%varData)
      curInfo => me%varData(me%varReordered(curVar))%info
      startLevel = 1
      endLevel   = 1
      IF (curInfo%ndims == 3) endLevel = curInfo%used_dimensions(2)
      nchunks = (endLevel - startLevel + 1 + restart_chunk_size -1) &
        &       / restart_chunk_size
      ! no valid time level -> no output:
      IF (.NOT.has_valid_time_level(curInfo, me%description%id, me%description%nnew, &
        &                          me%description%nnew_rcf)) CYCLE LOOP_VAR2
      LOOP_IO: DO ichunk = 1, nchunks
        ilevel_start = restart_chunk_size*(ichunk - 1) + 1
        nlevels      = MIN(restart_chunk_size * ichunk, endLevel) - ilevel_start + 1
        IF (timers_level >= 7) CALL timer_start(timer_write_restart_communication)
        CALL me%collector%fetch(curVar, ilevel_start, used_size, nlevels, me%RecvBuffer)
        IF (timers_level >= 7) CALL timer_stop(timer_write_restart_communication)
        IF (timers_level >= 7) CALL timer_start(timer_write_restart_io)
        SELECT CASE(curInfo%data_type)
        CASE(REAL_T)
          DO ilevel = 1, nlevels
            ptr2level_dp => me%recvBuffer%d(1+(ilevel-1)*used_size/nlevels:ilevel*used_size/nlevels)
            CALL streamWriteVarSlice(fileId, curInfo%cdiVarId,  ilevel+ilevel_start-2, ptr2level_dp, 0)
          END DO
          bytesWritten = bytesWritten + INT(used_size, i8) * 8_i8
        CASE(SINGLE_T)
          DO ilevel = 1, nlevels
            ptr2level_sp => me%RecvBuffer%s(1+(ilevel-1)*used_size/nlevels:ilevel*used_size/nlevels)
            CALL streamWriteVarSliceF(fileId, curInfo%cdiVarId, ilevel-2+ilevel_start, ptr2level_sp, 0)
          END DO
          bytesWritten = bytesWritten + INT(used_size, i8) * 4_i8
        CASE(INT_T)
          ! integer data type: copy to double precision field
          CALL ensureSize(me%RecvBuffer%d, used_size, no_copy)
!$OMP PARALLEL DO SCHEDULE(STATIC)
          DO i = 1, used_size
            me%RecvBuffer%d(i) = REAL(me%RecvBuffer%i(i), dp)
          END DO
          DO ilevel = 1, nlevels
            ptr2level_dp => me%RecvBuffer%d(1+(ilevel-1)*used_size/nlevels:ilevel*used_size/nlevels)
            CALL streamWriteVarSlice(fileId, curInfo%cdiVarId, ilevel-2+ilevel_start, ptr2level_dp, 0)
          END DO
          bytesWritten = bytesWritten + INT(used_size, i8) * 8_i8
        CASE DEFAULT
          CALL finish(routine, "Internal error! Variable "//TRIM(curInfo%name))
        END SELECT
        IF (timers_level >= 7) CALL timer_stop(timer_write_restart_io)
      END DO LOOP_IO
    END DO LOOP_VAR2
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
    INTEGER :: i, ierr

    CALL me%cellCollector%finalize()
    CALL me%edgeCollector%finalize()
    CALL me%vertCollector%finalize()
    CALL me%collector%finalize()
    CALL me%t_RestartPatchData%destruct()
    CALL me%glb_sendbuf%finalize()
    CALL dealloc_pointerarr_dp(me%RecvBuffer%d)
    IF (ASSOCIATED(me%RecvBuffer%s)) THEN
      DEALLOCATE(me%RecvBuffer%s, STAT=ierr)
      IF (ierr /= SUCCESS) &
        CALL finish(routine, "memory deallocation failure")
    END IF
    NULLIFY(me%RecvBuffer%s)
    CALL dealloc_pointerarr_int(me%RecvBuffer%i)
    CALL dealloc_pointerarr_dp(me%glb_indices_cell)
    CALL dealloc_pointerarr_dp(me%glb_indices_edge)
    CALL dealloc_pointerarr_dp(me%glb_indices_vert)
    CALL dealloc_pointerarr_int(me%varReordered)
  CONTAINS

    SUBROUTINE dealloc_pointerarr_dp(ptr)
      REAL(KIND=dp), POINTER, INTENT(INOUT) :: ptr(:)
      INTEGER :: ierr

      IF (ASSOCIATED(ptr)) THEN
        DEALLOCATE(ptr, STAT=ierr)
        IF (ierr /= SUCCESS) &
         CALL finish(routine, "memory deallocation failure")
      END IF
      NULLIFY(ptr)
    END SUBROUTINE dealloc_pointerarr_dp

    SUBROUTINE dealloc_pointerarr_int(ptr)
      INTEGER, POINTER, INTENT(INOUT) :: ptr(:)
      INTEGER :: ierr

      IF (ASSOCIATED(ptr)) THEN
        DEALLOCATE(ptr, STAT=ierr)
        IF (ierr /= SUCCESS) &
         CALL finish(routine, "memory deallocation failure")
      END IF
      NULLIFY(ptr)
    END SUBROUTINE dealloc_pointerarr_int
  END SUBROUTINE multifilePatchData_destruct

END MODULE mo_multifile_restart_patch_data
