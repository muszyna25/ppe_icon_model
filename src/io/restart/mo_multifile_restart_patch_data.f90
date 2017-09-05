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
  USE mo_cdi_constants,               ONLY: GRID_UNSTRUCTURED_CELL, GRID_UNSTRUCTURED_EDGE,                    &
    &                                       GRID_UNSTRUCTURED_VERT, ZA_SURFACE
  USE mo_cdi_ids,                     ONLY: t_CdiIds
  USE mtime,                          ONLY: datetime
  USE mo_decomposition_tools,         ONLY: t_grid_domain_decomp_info
  USE mo_exception,                   ONLY: finish
  USE mo_fortran_tools,               ONLY: t_ptr_2d, t_ptr_2d_sp, t_ptr_2d_int
  USE mo_impl_constants,              ONLY: SUCCESS, SINGLE_T, REAL_T, INT_T
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
  
  IMPLICIT NONE
  
  PUBLIC :: t_MultifilePatchData
  PUBLIC :: toMultifilePatchData

  PRIVATE
  
  TYPE t_CollectorRecvBuffer
    
    !each array buffers the complete data of a level of a variable
    !that is written by this PE
    REAL(dp), POINTER :: data_d(:)
    REAL(sp), POINTER :: data_s(:)
    INTEGER,  POINTER :: data_int(:)

    ! --- global index array
    !
    ! The global indices are transferred during the initial
    ! construction of the collectors. They are added to every payload
    ! file and must therefore be stored here.

    REAL(dp), ALLOCATABLE :: glb_indices(:)
        
  CONTAINS
    
    PROCEDURE :: construct => t_CollectorRecvBuffer_construct
    PROCEDURE :: finalize  => t_CollectorRecvBuffer_finalize
    
  END TYPE t_CollectorRecvBuffer

  
  TYPE, EXTENDS(t_RestartPatchData) :: t_MultifilePatchData
  
    ! --- receive buffer
  
    ! global receive buffers
    TYPE (t_CollectorRecvBuffer) :: glb_recvbuf_cell, &
      &                             glb_recvbuf_edge, &
      &                             glb_recvbuf_vert
    
    ! pointers to the three target buffers above, reusing the
    ! allocations:
    TYPE(t_ptr_1d_generic), ALLOCATABLE :: recvbuf(:)
    
    
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

  SUBROUTINE t_CollectorRecvBuffer_construct(me, writerRank, sourceRanks,   &
    &                                        isize, decomp_info, collector, &
    &                                        lthis_pe_active)

    CLASS(t_CollectorRecvBuffer),     INTENT(INOUT) :: me
    INTEGER,                          INTENT(IN)    :: writerRank
    INTEGER,                          INTENT(IN)    :: sourceRanks(:)
    INTEGER,                          INTENT(IN)    :: isize
    TYPE(t_grid_domain_decomp_info),  INTENT(IN)    :: decomp_info
    TYPE(t_CollectorIndices), TARGET, INTENT(INOUT) :: collector

    ! Flag. .FALSE. if this PE does not send any points (since it will
    ! write to another file).
    LOGICAL,                          INTENT(IN)    :: lthis_pe_active

    CHARACTER(*), PARAMETER :: routine = modname//"::t_CollectorRecvBuffer_construct"
    INTEGER                         :: ierrstat
    TYPE(t_grid_domain_decomp_info) :: dummyInfo

    IF (my_process_is_work()) THEN
      me%data_d => collector%construct(isize, decomp_info, writerRank, sourceRanks, lthis_pe_active)
    ELSE
      !must NOT access p_patch on dedicated restart processes
      me%data_d => collector%construct(0, dummyInfo, writerRank, sourceRanks, lthis_pe_active)
    END IF

    ALLOCATE(me%glb_indices(SIZE(me%data_d)), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "memory allocation failure")
    me%glb_indices(:) = me%data_d(:)

    ALLOCATE(me%data_s(SIZE(me%data_d)), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "memory allocation failure")
    ALLOCATE(me%data_int(SIZE(me%data_d)), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "memory allocation failure")
  END SUBROUTINE t_CollectorRecvBuffer_construct


  SUBROUTINE t_CollectorRecvBuffer_finalize(me)
    CLASS(t_CollectorRecvBuffer), INTENT(INOUT) :: me
    CHARACTER(*), PARAMETER :: routine = modname//":t_CollectorRecvBuffer_finalize"
    INTEGER :: ierrstat

    DEALLOCATE(me%data_d, me%data_s, me%data_int, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "memory allocation failure")

    IF (ALLOCATED(me%glb_indices)) THEN
      DEALLOCATE(me%glb_indices, STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, "memory deallocation failure")
    END IF
  END SUBROUTINE t_CollectorRecvBuffer_finalize


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

    IF (timers_level >= 7) CALL timer_start(timer_write_restart_setup)

    jg = me%description%id

    ! create a global receive buffer (for joint proc mode this is the
    ! only receive buffer which is actually allocated):
    CALL me%glb_recvbuf_cell%construct(writerRank, sourceRanks, p_patch(jg)%n_patch_cells, &
      &                                p_patch(jg)%cells%decomp_info, me%cellCollector, lthis_pe_active)
    CALL me%glb_recvbuf_edge%construct(writerRank, sourceRanks, p_patch(jg)%n_patch_edges, &
      &                                p_patch(jg)%edges%decomp_info, me%edgeCollector, lthis_pe_active)
    CALL me%glb_recvbuf_vert%construct(writerRank, sourceRanks, p_patch(jg)%n_patch_verts, &
      &                                p_patch(jg)%verts%decomp_info, me%vertCollector, lthis_pe_active)

    ! allocate variable indexed arrays
    ALLOCATE(me%recvbuf(SIZE(me%varData)),    &
      &      me%collectors(SIZE(me%varData)), STAT = error)
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
        me%recvbuf(curVar)%dp  => me%glb_recvbuf_cell%data_d
        me%recvbuf(curVar)%sp  => me%glb_recvbuf_cell%data_s
        me%recvbuf(curVar)%int => me%glb_recvbuf_cell%data_int
        ! note: the send buffers are allocated for all levels.
        CALL me%collectors(curVar)%construct(me%cellCollector, me%glb_sendbuf_cell, nlevs, itype, &
          &                                  ioffset_cell, sourceOffset_Cell)
      CASE(GRID_UNSTRUCTURED_EDGE)
        me%recvbuf(curVar)%dp  => me%glb_recvbuf_edge%data_d
        me%recvbuf(curVar)%sp  => me%glb_recvbuf_edge%data_s
        me%recvbuf(curVar)%int => me%glb_recvbuf_edge%data_int
        CALL me%collectors(curVar)%construct(me%edgeCollector, me%glb_sendbuf_edge, nlevs, itype, &
          &                                  ioffset_edge, sourceOffset_Edge)
      CASE(GRID_UNSTRUCTURED_VERT)
        me%recvbuf(curVar)%dp  => me%glb_recvbuf_vert%data_d
        me%recvbuf(curVar)%sp  => me%glb_recvbuf_vert%data_s
        me%recvbuf(curVar)%int => me%glb_recvbuf_vert%data_int
        CALL me%collectors(curVar)%construct(me%vertCollector, me%glb_sendbuf_vert, nlevs, itype, &
          &                                  ioffset_vert, sourceOffset_Vert)
      CASE DEFAULT
        CALL finish(routine, "Internal error!")
      END SELECT

    END DO

    ! finally, allocate the global send buffers:

    CALL me%glb_sendbuf_cell%construct(ioffset_cell)
    CALL me%glb_sendbuf_edge%construct(ioffset_edge)
    CALL me%glb_sendbuf_vert%construct(ioffset_vert)

    ! clean up
    DEALLOCATE(sourceOffset_cell, sourceOffset_Edge, sourceOffset_Vert, STAT=error)
    IF (error /= SUCCESS) CALL finish(routine, "memory deallocation failure")

    IF (timers_level >= 7) CALL timer_stop(timer_write_restart_setup)
  END SUBROUTINE multifilePatchData_createCollectors


  INTEGER FUNCTION defineIndexVar(cdiIds, hgrid, name) RESULT(varId)
    TYPE(t_CdiIds), INTENT(INOUT) :: cdiIds
    INTEGER,        INTENT(IN)    :: hgrid
    CHARACTER(*),   INTENT(IN)    :: name

    varId = cdiIds%defineVariableSimple(hgrid, ZA_SURFACE, DATATYPE_FLT64, name)
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
    IF (.NOT.ALLOCATED(me%recvbuf)) THEN
      CALL finish(routine, "assertion failed, call to createCollectors() missing")
    END IF

    IF (timers_level >= 7) CALL timer_start(timer_write_restart_io)

    description => me%description
    effectiveFilename = multifilePayloadPath(filename, description%id, ifile)

    CALL cdiIds%init()
    IF (ALLOCATED(description%opt_pvct)) THEN
      CALL cdiIds%openRestartAndCreateIds(effectiveFilename, FILETYPE_NC4,                               &
        &                                 SIZE(me%glb_recvbuf_cell%data_d),                              &
        &                                 SIZE(me%glb_recvbuf_vert%data_d),                              &
        &                                 SIZE(me%glb_recvbuf_edge%data_d),                              &
        &                                 description%cell_type,                                         &
        &                                 description%v_grid_defs(1:description%v_grid_count),           &
        &                                 description%opt_pvct)
    ELSE
      CALL cdiIds%openRestartAndCreateIds(effectiveFilename, FILETYPE_NC4,                               &
        &                                 SIZE(me%glb_recvbuf_cell%data_d),                              &
        &                                 SIZE(me%glb_recvbuf_vert%data_d),                              &
        &                                 SIZE(me%glb_recvbuf_edge%data_d),                              &
        &                                 description%cell_type,                                         &
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
    IF (SIZE(me%glb_recvbuf_cell%glb_indices) > 0) THEN
      CALL streamWriteVar(cdiIds%file, cellIndexVarId, me%glb_recvbuf_cell%glb_indices, 0)
      bytesWritten = bytesWritten + INT(SIZE(me%glb_recvbuf_cell%glb_indices), i8)*8_i8
    ELSE
      ! Note that we may have a cell count of 0 for this payload file;
      ! this may happen in the case of processor splitting. Still, we
      ! have to write *something*, otherwise the CDI will not store
      ! NetCDF attributes, which, in turn, is necessary to detect
      ! size-0 payload files!
      CALL streamWriteVar(cdiIds%file, cellIndexVarId, dummy_d, 0)
    END IF
    IF (SIZE(me%glb_recvbuf_edge%glb_indices) > 0) THEN
      CALL streamWriteVar(cdiIds%file, edgeIndexVarId, me%glb_recvbuf_edge%glb_indices, 0)
      bytesWritten = bytesWritten + INT(SIZE(me%glb_recvbuf_edge%glb_indices), i8)*8_i8
    ELSE
      CALL streamWriteVar(cdiIds%file, edgeIndexVarId, dummy_d, 0)
    END IF
    IF (SIZE(me%glb_recvbuf_vert%glb_indices) > 0) THEN
      CALL streamWriteVar(cdiIds%file, vertIndexVarId, me%glb_recvbuf_vert%glb_indices, 0)
      bytesWritten = bytesWritten + INT(SIZE(me%glb_recvbuf_vert%glb_indices), i8)*8_i8
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
    INTEGER, ALLOCATABLE            :: var_nlev(:)
    INTEGER                         :: startLevel, endLevel, error

    IF (.NOT.ALLOCATED(me%recvbuf)) THEN
      CALL finish(routine, "assertion failed, call to createCollectors() missing")
    END IF
    IF (.NOT. my_process_is_work()) THEN
      CALL finish(routine, "assertion failed.")
    END IF

    ! store for each variable the number of levels
    ALLOCATE(var_nlev(SIZE(me%varData)), STAT=error)
    IF (error /= SUCCESS) CALL finish(routine, "memory deallocation failure")
    DO curVar = 1, SIZE(me%varData)
      curInfo => me%varData(curVar)%info
      var_nlev(curVar) = 1
      IF (curInfo%ndims == 3) var_nlev(curVar) = curInfo%used_dimensions(2)
    END DO

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
      endLevel   = var_nlev(curVar)

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

    ! clean-up
    DEALLOCATE(var_nlev, STAT=error)
    IF (error /= SUCCESS) CALL finish(routine, "memory deallocation failure")

  END SUBROUTINE multifilePatchData_exposeData


  ! Collect the data, optionally writing it immediately.
  !
  SUBROUTINE multifilePatchData_collectData(me, fileId, bytesWritten)
    CLASS(t_MultifilePatchData), INTENT(INOUT) :: me
    INTEGER,                     INTENT(IN)    :: fileId
    INTEGER(i8),                 INTENT(INOUT) :: bytesWritten

    CHARACTER(*), PARAMETER :: routine = modname//":multifilePatchData_collectData"
    INTEGER                         :: curVar, curLevel
    TYPE(t_var_metadata), POINTER   :: curInfo
    INTEGER, ALLOCATABLE            :: var_nlev(:)
    REAL(dp), ALLOCATABLE           :: buffer_dp(:)
    INTEGER                         :: startLevel, endLevel, error

    IF (.NOT.ALLOCATED(me%recvbuf)) THEN
      CALL finish(routine, "assertion failed, call to createCollectors() missing")
    END IF

    ! store for each variable the number of levels
    ALLOCATE(var_nlev(SIZE(me%varData)), STAT=error)
    IF (error /= SUCCESS) CALL finish(routine, "memory deallocation failure")
    DO curVar = 1, SIZE(me%varData)
      curInfo => me%varData(curVar)%info
      var_nlev(curVar) = 1
      IF (curInfo%ndims == 3) var_nlev(curVar) = curInfo%used_dimensions(2)
    END DO
   
    LOOP_VAR2 : DO curVar = 1, SIZE(me%varData)

      curInfo => me%varData(curVar)%info

      startLevel = 1
      endLevel   = var_nlev(curVar)

      ! no valid time level -> no output:
      IF (.NOT.has_valid_time_level(curInfo, me%description%id, me%description%nnew, &
        &                          me%description%nnew_rcf)) CYCLE LOOP_VAR2

      LOOP_IO: DO curLevel=startLevel, endLevel

        SELECT CASE(curInfo%data_type)
        CASE(REAL_T)
          IF (SIZE(me%recvbuf(curVar)%dp) > 0) THEN
            IF (timers_level >= 7) CALL timer_start(timer_write_restart_communication)
            CALL me%collectors(curVar)%receiveBuffer(curLevel, me%recvbuf(curVar)%dp)
            IF (timers_level >= 7) CALL timer_stop(timer_write_restart_communication)
            IF (timers_level >= 7) CALL timer_start(timer_write_restart_io)
            CALL streamWriteVarSlice(fileId, curInfo%cdiVarId, curLevel-1, &
              &                      me%recvbuf(curVar)%dp, 0)
            IF (timers_level >= 7) CALL timer_stop(timer_write_restart_io)
            bytesWritten = bytesWritten + INT(SIZE(me%recvbuf(curVar)%dp), i8)*8_i8
          END IF
          !
        CASE(SINGLE_T)
          IF (SIZE(me%recvbuf(curVar)%sp) > 0) THEN
            IF (timers_level >= 7) CALL timer_start(timer_write_restart_communication)
            CALL me%collectors(curVar)%receiveBuffer(curLevel, me%recvbuf(curVar)%sp)
            IF (timers_level >= 7) CALL timer_stop(timer_write_restart_communication)
            IF (timers_level >= 7) CALL timer_start(timer_write_restart_io)
            CALL streamWriteVarSliceF(fileId, curInfo%cdiVarId, curLevel-1, &
              &                       me%recvbuf(curVar)%sp, 0)
            IF (timers_level >= 7) CALL timer_stop(timer_write_restart_io)
            bytesWritten = bytesWritten + INT(SIZE(me%recvbuf(curVar)%sp), i8)*4_i8
          END IF
          !               
        CASE(INT_T)
          ! integer data type: copy to double precision field
          IF (SIZE(me%recvbuf(curVar)%int) > 0) THEN
            IF (timers_level >= 7) CALL timer_start(timer_write_restart_communication)
            CALL me%collectors(curVar)%receiveBuffer(curLevel, me%recvbuf(curVar)%int)
            IF (timers_level >= 7) CALL timer_stop(timer_write_restart_communication)

            ALLOCATE(buffer_dp(SIZE(me%recvbuf(curVar)%int)), STAT = error)
            IF (error /= SUCCESS) CALL finish(routine, "memory allocation failure")
            buffer_dp(:) = REAL(me%recvbuf(curVar)%int(:), dp)
            IF (timers_level >= 7) CALL timer_start(timer_write_restart_io)
            CALL streamWriteVarSlice(fileId, curInfo%cdiVarId, curLevel-1, &
              &                      buffer_dp, 0)
            IF (timers_level >= 7) CALL timer_stop(timer_write_restart_io)
            bytesWritten = bytesWritten + INT(SIZE(me%recvbuf(curVar)%int), i8)*8_i8
            DEALLOCATE(buffer_dp, STAT = error)
            IF (error /= SUCCESS) CALL finish(routine, "memory deallocation failure")
          END IF
          !
        CASE DEFAULT
          CALL finish(routine, "Internal error! Variable "//TRIM(curInfo%name))
        END SELECT

      END DO LOOP_IO

    END DO LOOP_VAR2

    ! clean-up
    DEALLOCATE(var_nlev, STAT=error)
    IF (error /= SUCCESS) CALL finish(routine, "memory deallocation failure")

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

    IF (.NOT.ALLOCATED(me%recvbuf)) THEN
      CALL finish(routine, "assertion failed, call to createCollectors() missing")
    END IF

    CALL me%glb_recvbuf_cell%finalize()
    CALL me%glb_recvbuf_edge%finalize()
    CALL me%glb_recvbuf_vert%finalize()

    DEALLOCATE(me%recvbuf, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "memory deallocation failure")

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

  END SUBROUTINE multifilePatchData_destruct

END MODULE mo_multifile_restart_patch_data
