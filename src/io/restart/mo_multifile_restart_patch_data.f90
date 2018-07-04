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
  USE mo_multifile_restart_collector, ONLY: t_MultifileRestartCollector,                                       &
    &                                       t_CollectorIndices, t_CollectorSendBuffer
  USE mo_multifile_restart_util,      ONLY: multifilePayloadPath, kVarName_globalCellIndex,                    &
    &                                       kVarName_globalEdgeIndex, kVarName_globalVertIndex,                &
    &                                       t_ptr_1d_generic
  USE mo_restart_attributes,          ONLY: t_RestartAttributeList
  USE mo_restart_descriptor,          ONLY: t_RestartPatchData
  USE mo_restart_patch_description,   ONLY: t_restart_patch_description
  USE mo_restart_util,                ONLY: my_process_is_restart_writer,                                      &
    &                                       isDedicatedProcMode, t_restart_args 
  USE mo_restart_var_data,            ONLY: has_valid_time_level, getLevelPointers
  USE mo_timer,                       ONLY: timer_start, timer_stop, timer_write_restart_io,                   &
    &                                       timer_write_restart_communication, timers_level,                   &
    &                                       timer_write_restart_setup
  USE mo_var_metadata_types,          ONLY: t_var_metadata
  USE mo_parallel_config,             ONLY: restart_chunk_size

  IMPLICIT NONE
  
  PUBLIC :: t_MultifilePatchData
  PUBLIC :: toMultifilePatchData

  PRIVATE
  
  TYPE, EXTENDS(t_RestartPatchData) :: t_MultifilePatchData
    ! global receive buffers
    REAL(dp), POINTER :: commonRecvBuf_dp(:)
    REAL(sp), POINTER :: commonRecvBuf_sp(:)
    INTEGER,  POINTER :: commonRecvBuf_int(:)
    REAL(dp), POINTER :: glb_indices_cell(:)
    REAL(dp), POINTER :: glb_indices_edge(:)
    REAL(dp), POINTER :: glb_indices_vert(:)
    ! --- send buffer
    ! global sender buffers attached to single MPI memory windows; all
    ! collectors point into these data structures.
    TYPE (t_CollectorSendBuffer) :: glb_sendbuf_cell, &
      &                             glb_sendbuf_edge, &
      &                             glb_sendbuf_vert
    ! --- "collector" objects     
    ! three index sets which we need to aggregate the data on the
    ! restart writers
    TYPE(t_CollectorIndices) :: cellCollector, edgeCollector, vertCollector
    ! for each variable: "collector" aggregating the data, including
    ! the send buffer for each variable
    TYPE(t_MultifileRestartCollector), ALLOCATABLE :: collectors(:)
  CONTAINS
    ! override, collective across both work and restart procs, must
    ! be followed by a CALL to createCollectors()
    PROCEDURE :: construct => multifilePatchData_construct
    ! no communication, but entered by all procs
    PROCEDURE :: createCollectors => multifilePatchData_createCollectors
    PROCEDURE :: openPayloadFile => multifilePatchData_openPayloadFile  ! restart writers only
    ! collective across both work and restart procs:
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

    NULLIFY(me%commonRecvBuf_dp, me%commonRecvBuf_sp, me%commonRecvBuf_int, &
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
    INTEGER                         :: curVar, error, jg, nlevs, itype
    INTEGER(i8)                     :: ioffset_cell(TYPE_LEN), ioffset_edge(TYPE_LEN), &
      &                                ioffset_vert(TYPE_LEN)
    TYPE(t_var_metadata), POINTER   :: curInfo
    INTEGER, ALLOCATABLE            :: sourceOffset_cell(:,:), sourceOffset_Edge(:,:), &
      &                                sourceOffset_Vert(:,:)
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
    ! allocate variable indexed arrays
    ALLOCATE(me%collectors(SIZE(me%varData)), STAT = error)
    IF (error /= SUCCESS) CALL finish(routine, "memory allocation failure")

    IF (my_process_is_restart_writer()) THEN
      ALLOCATE(sourceOffset_Cell(TYPE_LEN, me%cellCollector%sourceProcCount), &
        &      sourceOffset_Edge(TYPE_LEN, me%edgeCollector%sourceProcCount), &
        &      sourceOffset_Vert(TYPE_LEN, me%vertCollector%sourceProcCount), STAT=error)
    ELSE
      ALLOCATE(sourceOffset_Cell(TYPE_LEN, 1), &
        &      sourceOffset_Edge(TYPE_LEN, 1), &
        &      sourceOffset_Vert(TYPE_LEN, 1), STAT=error)
    END IF
    IF (error /= SUCCESS) CALL finish(routine, "memory allocation failure")
    sourceOffset_Cell(:,:) = 0
    sourceOffset_Edge(:,:) = 0
    sourceOffset_Vert(:,:) = 0
    ioffset_cell(:) = 0
    ioffset_edge(:) = 0
    ioffset_vert(:) = 0
    DO curVar = 1, SIZE(me%varData)
      curInfo => me%varData(curVar)%info
      itype = curInfo%data_type
      !get the effective level count
      nlevs = 1
      IF (curInfo%ndims == 3) nlevs = curInfo%used_dimensions(2)
      IF (curInfo%ndims > 3) CALL finish(routine, "ndims > 3 is not supported")
      ! get the correct collector and buffer for the first level
      SELECT CASE(curInfo%hgrid)
      CASE(GRID_UNSTRUCTURED_CELL)
        ! note: the send buffers are allocated for all levels.
        CALL me%collectors(curVar)%construct(me%cellCollector, me%glb_sendbuf_cell, nlevs, itype, &
          &                                  ioffset_cell, sourceOffset_Cell)
      CASE(GRID_UNSTRUCTURED_EDGE)
        CALL me%collectors(curVar)%construct(me%edgeCollector, me%glb_sendbuf_edge, nlevs, itype, &
          &                                  ioffset_edge, sourceOffset_Edge)
      CASE(GRID_UNSTRUCTURED_VERT)
        CALL me%collectors(curVar)%construct(me%vertCollector, me%glb_sendbuf_vert, nlevs, itype, &
          &                                  ioffset_vert, sourceOffset_Vert)
      CASE DEFAULT
        CALL finish(routine, "Internal error!")
      END SELECT
    END DO
    ! finally, allocate the global send buffers:
    CALL me%glb_sendbuf_cell%construct(ioffset_cell, me%cellCollector)
    CALL me%glb_sendbuf_edge%construct(ioffset_edge, me%edgeCollector)
    CALL me%glb_sendbuf_vert%construct(ioffset_vert, me%vertCollector)
    ! clean up
    DEALLOCATE(sourceOffset_cell, sourceOffset_Edge, sourceOffset_Vert, STAT=error)
    IF (error /= SUCCESS) CALL finish(routine, "memory deallocation failure")
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
    CALL me%glb_sendbuf_cell%start_local_access()
    CALL me%glb_sendbuf_edge%start_local_access()
    CALL me%glb_sendbuf_vert%start_local_access()
    IF (timers_level >= 7) CALL timer_stop(timer_write_restart_communication)
  END SUBROUTINE multifilePatchData_start_local_access

  ! Start MPI epoch where windows may not be accessed locally and
  ! not via MPI_PUT operations.
  SUBROUTINE multifilePatchData_start_remote_access(me)
    CLASS(t_MultifilePatchData), INTENT(INOUT) :: me

    IF (timers_level >= 7) CALL timer_start(timer_write_restart_communication)
    CALL me%glb_sendbuf_cell%start_remote_access()
    CALL me%glb_sendbuf_edge%start_remote_access()
    CALL me%glb_sendbuf_vert%start_remote_access()
    IF (timers_level >= 7) CALL timer_stop(timer_write_restart_communication)
  END SUBROUTINE multifilePatchData_start_remote_access

  ! Collect the data, optionally writing it immediately.
  !
  SUBROUTINE multifilePatchData_exposeData(me)
    CLASS(t_MultifilePatchData), INTENT(INOUT) :: me
    CHARACTER(*), PARAMETER :: routine = modname//":multifilePatchData_exposeData"
    INTEGER                         :: curVar, curLevel
    TYPE(t_var_metadata), POINTER   :: curInfo
    TYPE(t_ptr_2d),     ALLOCATABLE :: dataPointers_d(:)
    TYPE(t_ptr_2d_sp),  ALLOCATABLE :: dataPointers_s(:)
    TYPE(t_ptr_2d_int), ALLOCATABLE :: dataPointers_int(:)
    INTEGER                         :: startLevel, endLevel, error

    IF (.NOT. my_process_is_work()) THEN
      CALL finish(routine, "assertion failed.")
    END IF
    IF (timers_level >= 7) CALL timer_start(timer_write_restart_communication)
    LOOP_VAR1 : DO curVar = 1, SIZE(me%varData)
      curInfo => me%varData(curVar)%info
      ! no valid time level -> no output:
      IF (.NOT.has_valid_time_level(curInfo, me%description%id, me%description%nnew, &
        &                          me%description%nnew_rcf)) CYCLE LOOP_VAR1
      ! get the actual data pointers
      SELECT CASE(curInfo%data_type)
      CASE(REAL_T)
        CALL getLevelPointers(curInfo, me%varData(curVar)%r_ptr, dataPointers_d)
      CASE(SINGLE_T)
        CALL getLevelPointers(curInfo, me%varData(curVar)%s_ptr, dataPointers_s)
      CASE(INT_T)
        CALL getLevelPointers(curInfo, me%varData(curVar)%i_ptr, dataPointers_int)
      CASE DEFAULT
        CALL finish(routine, "Internal error! Variable "//TRIM(curInfo%name))
      END SELECT
      startLevel = 1
      endLevel   = 1
      IF (curInfo%ndims == 3) endLevel = curInfo%used_dimensions(2)

      LOOP_COLLECT: DO curLevel = startLevel, endLevel
        SELECT CASE(curInfo%data_type)
        CASE(REAL_T)
          CALL me%collectors(curVar)%sendField(dataPointers_d(curLevel)%p,   curLevel)
        CASE(SINGLE_T)
          CALL me%collectors(curVar)%sendField(dataPointers_s(curLevel)%p,   curLevel)
        CASE(INT_T)
          CALL me%collectors(curVar)%sendField(dataPointers_int(curLevel)%p, curLevel)
        CASE DEFAULT
          CALL finish(routine, "Internal error! Variable "//TRIM(curInfo%name))
        END SELECT
      END DO LOOP_COLLECT
    END DO LOOP_VAR1
    ! clean-up
    IF (ALLOCATED(dataPointers_d))    DEALLOCATE(dataPointers_d)
    IF (ALLOCATED(dataPointers_s))    DEALLOCATE(dataPointers_s)
    IF (ALLOCATED(dataPointers_int))  DEALLOCATE(dataPointers_int)
    IF (timers_level >= 7) CALL timer_stop(timer_write_restart_communication)
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
      curInfo => me%varData(curVar)%info
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
        SELECT CASE(curInfo%data_type)
        CASE(REAL_T)
          CALL me%collectors(curVar)%receiveBuffer(ilevel_start, me%commonRecvBuf_dp, used_size, nlevels)
          IF (timers_level >= 7) CALL timer_stop(timer_write_restart_communication)
          IF (timers_level >= 7) CALL timer_start(timer_write_restart_io)
          DO ilevel = 1, nlevels
            ptr2level_dp => me%commonRecvBuf_dp(1+(ilevel-1)*used_size/nlevels:ilevel*used_size/nlevels)
            CALL streamWriteVarSlice(fileId, curInfo%cdiVarId,  ilevel+ilevel_start-2, ptr2level_dp, 0)
          END DO
          bytesWritten = bytesWritten + INT(used_size, i8) * 8_i8
          !
        CASE(SINGLE_T)
          CALL me%collectors(curVar)%receiveBuffer(ilevel_start, me%commonRecvBuf_sp, used_size, nlevels)
          IF (timers_level >= 7) CALL timer_stop(timer_write_restart_communication)
          IF (timers_level >= 7) CALL timer_start(timer_write_restart_io)
          DO ilevel = 1, nlevels
            ptr2level_sp => me%commonRecvBuf_sp(1+(ilevel-1)*used_size/nlevels:ilevel*used_size/nlevels)
            CALL streamWriteVarSliceF(fileId, curInfo%cdiVarId, ilevel-2+ilevel_start, ptr2level_sp, 0)
          END DO
          bytesWritten = bytesWritten + INT(used_size, i8) * 4_i8
          !               
        CASE(INT_T)
          ! integer data type: copy to double precision field
          CALL me%collectors(curVar)%receiveBuffer(ilevel_start, me%commonRecvBuf_int, used_size, nlevels)
          IF (timers_level >= 7) CALL timer_stop(timer_write_restart_communication)
          IF (timers_level >= 7) CALL timer_start(timer_write_restart_io)
          CALL ensureSize(me%commonRecvBuf_dp, used_size, no_copy)
!$OMP PARALLEL DO SCHEDULE(STATIC)
          DO i = 1, used_size
            me%commonRecvBuf_dp(i) = REAL(me%commonRecvBuf_int(i), dp)
          END DO
          DO ilevel = 1, nlevels
            ptr2level_dp => me%commonRecvBuf_dp(1+(ilevel-1)*used_size/nlevels:ilevel*used_size/nlevels)
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
    INTEGER :: i, ierrstat

    CALL me%cellCollector%finalize()
    CALL me%edgeCollector%finalize()
    CALL me%vertCollector%finalize()
    DO i=1,SIZE(me%collectors)
      CALL me%collectors(i)%finalize()
    END DO
    DEALLOCATE(me%collectors, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "memory deallocation failure")
    CALL me%t_RestartPatchData%destruct()
    CALL me%glb_sendbuf_cell%finalize()
    CALL me%glb_sendbuf_edge%finalize()
    CALL me%glb_sendbuf_vert%finalize()
    CALL dealloc_pointerarr_dp(me%commonRecvBuf_dp)
    IF (ASSOCIATED(me%commonRecvBuf_sp)) THEN
      DEALLOCATE(me%commonRecvBuf_sp, STAT=ierrstat)
      IF (ierrstat /= SUCCESS) &
        CALL finish(routine, "memory deallocation failure")
    END IF
    IF (ASSOCIATED(me%commonRecvBuf_int)) THEN
      DEALLOCATE(me%commonRecvBuf_int, STAT=ierrstat)
      IF (ierrstat /= SUCCESS) &
        CALL finish(routine, "memory deallocation failure")
    END IF
    NULLIFY(me%commonRecvBuf_sp, me%commonRecvBuf_int)
    CALL dealloc_pointerarr_dp(me%glb_indices_cell)
    CALL dealloc_pointerarr_dp(me%glb_indices_edge)
    CALL dealloc_pointerarr_dp(me%glb_indices_vert)
  CONTAINS

    SUBROUTINE dealloc_pointerarr_dp(ptr)
      REAL(KIND=dp), POINTER, INTENT(INOUT) :: ptr(:)
      INTEGER :: ierrstat

      IF (ASSOCIATED(ptr)) THEN
        DEALLOCATE(ptr, STAT=ierrstat)
        IF (ierrstat /= SUCCESS) &
         CALL finish(routine, "memory deallocation failure")
      END IF
      NULLIFY(ptr)
    END SUBROUTINE dealloc_pointerarr_dp

  END SUBROUTINE multifilePatchData_destruct

END MODULE mo_multifile_restart_patch_data
