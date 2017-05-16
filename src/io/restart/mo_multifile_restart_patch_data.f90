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
    USE mo_fortran_tools,               ONLY: t_ptr_1d_generic_ptr_1d, t_ptr_2d, t_ptr_2d_sp, t_ptr_2d_int
    USE mo_impl_constants,              ONLY: SUCCESS, SINGLE_T, REAL_T, INT_T
    USE mo_kind,                        ONLY: dp, sp, i8
    USE mo_model_domain,                ONLY: p_patch
    USE mo_mpi,                         ONLY: my_process_is_work, my_process_is_restart
    USE mo_multifile_restart_collector, ONLY: t_MultifileRestartCollector, t_MultifileRestartCollector_ptr
    USE mo_multifile_restart_util,      ONLY: multifilePayloadPath, kVarName_globalCellIndex,                    &
      &                                       kVarName_globalEdgeIndex, kVarName_globalVertIndex
    USE mo_restart_attributes,          ONLY: t_RestartAttributeList
    USE mo_restart_descriptor,          ONLY: t_RestartPatchData
    USE mo_restart_patch_description,   ONLY: t_restart_patch_description
    USE mo_restart_util,                ONLY: restartWorkProcId, my_process_is_restart_writer,                   &
      &                                       isDedicatedProcMode, t_restart_args
    USE mo_restart_var_data,            ONLY: has_valid_time_level, getLevelPointers
    USE mo_timer,                       ONLY: timer_start, timer_stop, timer_write_restart_io,                   &
      &                                       timer_write_restart_communication, timers_level
    USE mo_var_metadata_types,          ONLY: t_var_metadata

    IMPLICIT NONE

    PUBLIC :: t_MultifilePatchData
    PUBLIC :: toMultifilePatchData

    PRIVATE

    TYPE, EXTENDS(t_RestartPatchData) :: t_MultifilePatchData
      !these arrays are used on the restart writers to buffer the
      !data between collecting it from the (other) workers AND
      !writing it to the file
      
      !each array buffers the complete data of a level of a variable
      !that is written by this PE
      REAL(dp), POINTER :: cellData_d(:), edgeData_d(:), vertData_d(:)

      !each array buffers the complete data of a level of a variable
      !that is written by this PE
      REAL(sp), POINTER :: cellData_s(:), edgeData_s(:), vertData_s(:)   

      !each array buffers the complete data of a level of a variable
      !that is written by this PE
      INTEGER,  POINTER :: cellData_int(:), edgeData_int(:), vertData_int(:)   

      !on dedicated restart writers, these have separate allocations;
      !on joint proc writers, they are pointers to the three TARGET
      !buffers above, reusing the allocations
      TYPE(t_ptr_1d_generic_ptr_1d), ALLOCATABLE :: buffers(:)

      !the three collectors we need to aggregate the data on the restart writers
      TYPE(t_MultifileRestartCollector), POINTER :: cellCollector, edgeCollector, vertCollector

      !POINTERs to the three targets above, same index as varData(:)
      TYPE(t_MultifileRestartCollector_ptr), ALLOCATABLE :: collectors(:)
    CONTAINS

      ! override, collective across both work AND restart procs, must
      ! be followed by a CALL to createCollectors()
      PROCEDURE :: construct => multifilePatchData_construct

      ! no communication, but entered by all procs
      PROCEDURE :: createCollectors => multifilePatchData_createCollectors

      PROCEDURE :: openPayloadFile => multifilePatchData_openPayloadFile  ! restart writers only

      ! collective across both work AND restart procs:
      PROCEDURE :: collectData => multifilePatchData_collectData

      ! dedicated restart writers only, no communication
      PROCEDURE :: writeFromBuffer => multifilePatchData_writeFromBuffer

      PROCEDURE :: writeFile => multifilePatchData_writeFile  ! override
      PROCEDURE :: destruct => multifilePatchData_destruct  ! override
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
        CHARACTER(*), INTENT(IN) :: modelType
        INTEGER, INTENT(IN) :: domain

        CHARACTER(*), PARAMETER :: routine = modname//":multifilePatchData_construct"

        CALL me%t_RestartPatchData%construct(modelType, domain)

        IF(isDedicatedProcMode()) CALL me%transferToRestart()
    END SUBROUTINE multifilePatchData_construct

    !this should have been part of the constructor, however, fortran
    !does NOT allow changing the signature WHILE overriding a method,
    !nor does it allow overloading the overridden method, so we need
    !to USE a different NAME for this
    SUBROUTINE multifilePatchData_createCollectors(me, writerRank, sourceRanks)
        CLASS(t_MultifilePatchData), INTENT(INOUT) :: me
        INTEGER, INTENT(IN) :: writerRank
        INTEGER, INTENT(IN) :: sourceRanks(:)

        INTEGER :: curVar, levelCount, curLevel, levelSize, error, jg
        TYPE(t_var_metadata), POINTER :: curInfo
        TYPE(t_grid_domain_decomp_info) :: dummyInfo
        CHARACTER(*), PARAMETER :: routine = modname//":multifilePatchData_createCollectors"

        IF(timers_level >= 7) CALL timer_start(timer_write_restart_communication)

        !ALLOCATE the collectors
        ALLOCATE(me%cellCollector, STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation failure")
        ALLOCATE(me%edgeCollector, STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation failure")
        ALLOCATE(me%vertCollector, STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation failure")

        !initialize them
        jg = me%description%id
        IF(my_process_is_work()) THEN
            me%cellData_d => me%cellCollector%construct(p_patch(jg)%n_patch_cells, p_patch(jg)%cells%decomp_info, &
                                                     &writerRank, sourceRanks)
            me%edgeData_d => me%edgeCollector%construct(p_patch(jg)%n_patch_edges, p_patch(jg)%edges%decomp_info, &
                                                     &writerRank, sourceRanks)
            me%vertData_d => me%vertCollector%construct(p_patch(jg)%n_patch_verts, p_patch(jg)%verts%decomp_info, &
                                                     &writerRank, sourceRanks)
        ELSE
            !must NOT access p_patch on dedicated restart processes
            me%cellData_d => me%cellCollector%construct(0, dummyInfo, writerRank, sourceRanks)
            me%edgeData_d => me%edgeCollector%construct(0, dummyInfo, writerRank, sourceRanks)
            me%vertData_d => me%vertCollector%construct(0, dummyInfo, writerRank, sourceRanks)
        END IF
        ALLOCATE(me%cellData_s(SIZE(me%cellData_d)), STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation failure")
        ALLOCATE(me%edgeData_s(SIZE(me%edgeData_d)), STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation failure")
        ALLOCATE(me%vertData_s(SIZE(me%vertData_d)), STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation failure")

        ALLOCATE(me%cellData_int(SIZE(me%cellData_d)), STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation failure")
        ALLOCATE(me%edgeData_int(SIZE(me%edgeData_d)), STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation failure")
        ALLOCATE(me%vertData_int(SIZE(me%vertData_d)), STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation failure")

        !ALLOCATE the two variable indexed arrays
        ALLOCATE(me%buffers(SIZE(me%varData)), STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation failure")
        ALLOCATE(me%collectors(SIZE(me%varData)), STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation failure")

        !fill the two arrays
        DO curVar = 1, SIZE(me%varData)
            me%buffers(curVar)%p => NULL()
            me%collectors(curVar)%p => NULL()

            curInfo => me%varData(curVar)%info

            !get the effective level count
            levelCount = 1
            IF(curInfo%ndims == 3) levelCount = curInfo%used_dimensions(2)
            IF(curInfo%ndims > 3) CALL finish(routine, "ndims > 3 is not supported")

            !ALLOCATE the level array
            ALLOCATE(me%buffers(curVar)%p(levelCount), STAT = error)
            IF(error /= SUCCESS) CALL finish(routine, "memory allocation failure")

            !get the correct collector AND buffer for the first level
            SELECT CASE(curInfo%hgrid)
                CASE(GRID_UNSTRUCTURED_CELL)
                    me%buffers(curVar)%p(1)%dp  => me%cellData_d
                    me%buffers(curVar)%p(1)%sp  => me%cellData_s
                    me%buffers(curVar)%p(1)%int => me%cellData_int
                    me%collectors(curVar)%p => me%cellCollector
                CASE(GRID_UNSTRUCTURED_EDGE)
                    me%buffers(curVar)%p(1)%dp  => me%edgeData_d
                    me%buffers(curVar)%p(1)%sp  => me%edgeData_s
                    me%buffers(curVar)%p(1)%int => me%edgeData_int
                    me%collectors(curVar)%p => me%edgeCollector
                CASE(GRID_UNSTRUCTURED_VERT)
                    me%buffers(curVar)%p(1)%dp  => me%vertData_d
                    me%buffers(curVar)%p(1)%sp  => me%vertData_s
                    me%buffers(curVar)%p(1)%int => me%vertData_int
                    me%collectors(curVar)%p => me%vertCollector
            END SELECT

            !init the rest of the level pointers, either by aliasing
            !the buffer for the first level, OR by creating separate
            !allocations (on dedicated restart writers)
            levelSize = SIZE(me%buffers(curVar)%p(1)%dp)
            DO curLevel = 1, levelCount
                IF(my_process_is_restart()) THEN
                  !dedicated proc mode: each level has its own
                  !allocation, so that we can first fill all the
                  !me%buffers, AND THEN write their contents
                  !asynchronously WHILE the work processes continue
                  !the integration
                  
                  me%buffers(curVar)%p(curLevel)%dp  => NULL()                    
                  me%buffers(curVar)%p(curLevel)%sp  => NULL()                    
                  me%buffers(curVar)%p(curLevel)%int => NULL()                    
                  
                  SELECT CASE(curInfo%data_type)
                  CASE(REAL_T)
                    ALLOCATE(me%buffers(curVar)%p(curLevel)%dp(levelSize), STAT = error)
                    !
                  CASE(SINGLE_T)
                    ALLOCATE(me%buffers(curVar)%p(curLevel)%sp(levelSize), STAT = error)                    
                    !
                  CASE(INT_T)
                    ALLOCATE(me%buffers(curVar)%p(curLevel)%int(levelSize), STAT = error)
                    !
                  CASE DEFAULT
                    CALL finish(routine, "Internal error! Variable "//TRIM(curInfo%name))
                  END SELECT
                  
                  IF(error /= SUCCESS) CALL finish(routine, "memory allocation failure")
                ELSE
                    !joint proc mode/pure workers: just alias the
                    !me%buffers with the first level, which is aliased
                    !with the other variables of the same SIZE.
                    !Separate allocations must be avoided IN joint
                    !proc mode to keep our memory footprint down.
                    me%buffers(curVar)%p(curLevel)%dp  => me%buffers(curVar)%p(1)%dp
                    me%buffers(curVar)%p(curLevel)%sp  => me%buffers(curVar)%p(1)%sp
                    me%buffers(curVar)%p(curLevel)%int => me%buffers(curVar)%p(1)%int
                END IF
            END DO
        END DO

        IF(timers_level >= 7) CALL timer_stop(timer_write_restart_communication)
    END SUBROUTINE multifilePatchData_createCollectors

    INTEGER FUNCTION defineIndexVar(cdiIds, hgrid, NAME) RESULT(varId)
        TYPE(t_CdiIds), INTENT(INOUT) :: cdiIds
        INTEGER, VALUE :: hgrid
        CHARACTER(*), INTENT(IN) :: NAME

        varId = cdiIds%defineVariableSimple(hgrid, ZA_SURFACE, DATATYPE_FLT64, NAME)
    END FUNCTION defineIndexVar

    FUNCTION multifilePatchData_openPayloadFile(me, filename, this_datetime, bytesWritten) RESULT(cdiIds)
        TYPE(t_CdiIds) :: cdiIds
        CLASS(t_MultifilePatchData), TARGET, INTENT(INOUT) :: me
        CHARACTER(*), INTENT(IN) :: filename
        TYPE(datetime), POINTER, INTENT(IN) :: this_datetime
        INTEGER(i8), INTENT(INOUT) :: bytesWritten

        CHARACTER(:), ALLOCATABLE :: effectiveFilename
        CLASS(t_restart_patch_description), POINTER :: description
        INTEGER :: cellIndexVarId, edgeIndexVarId, vertIndexVarId, curVar
        CHARACTER(*), PARAMETER :: routine = modname//":multifilePatchData_openPayloadFile"
        REAL(dp) :: dummy_d(1)

        dummy_d(1) = 0._dp
        IF(.NOT.ALLOCATED(me%buffers)) CALL finish(routine, "assertion failed, call to createCollectors() missing")

        IF(timers_level >= 7) CALL timer_start(timer_write_restart_io)

        effectiveFilename = multifilePayloadPath(filename, me%description%id, restartWorkProcId())
        description => me%description

        CALL cdiIds%init()
        IF(ALLOCATED(description%opt_pvct)) THEN
            CALL cdiIds%openRestartAndCreateIds(effectiveFilename, FILETYPE_NC4,                               &
              &                                 SIZE(me%cellData_d), SIZE(me%vertData_d), SIZE(me%edgeData_d), &
              &                                 description%cell_type,                                         &
              &                                 description%v_grid_defs(1:description%v_grid_count),           &
              &                                 description%opt_pvct)
        ELSE
            CALL cdiIds%openRestartAndCreateIds(effectiveFilename, FILETYPE_NC4,                               &
              &                                 SIZE(me%cellData_d), SIZE(me%vertData_d), SIZE(me%edgeData_d), &
              &                                 description%cell_type,                                         &
              &                                 description%v_grid_defs(1:description%v_grid_count))
        END IF
        cellIndexVarId = defineIndexVar(cdiIds, GRID_UNSTRUCTURED_CELL, kVarName_globalCellIndex)
        edgeIndexVarId = defineIndexVar(cdiIds, GRID_UNSTRUCTURED_EDGE, kVarName_globalEdgeIndex)
        vertIndexVarId = defineIndexVar(cdiIds, GRID_UNSTRUCTURED_VERT, kVarName_globalVertIndex)

        !define all restart variables with a valid time level
        DO curVar = 1, SIZE(me%varData)
            IF(has_valid_time_level(me%varData(curVar)%info, me%description%id, description%nnew, &
              &                     description%nnew_rcf)) THEN
                CALL cdiIds%defineVariable(me%varData(curVar)%info)
            END IF
        END DO

        CALL cdiIds%finalizeVlist(this_datetime)

        !write the arrays with the global indices of the points
        IF (SIZE(me%cellData_d) > 0) THEN
          CALL streamWriteVar(cdiIds%file, cellIndexVarId, me%cellData_d, 0)
          bytesWritten = bytesWritten + INT(SIZE(me%cellData_d), i8)*8_i8
        ELSE
          ! Note that we may have a cell count of 0 for this payload
          ! file; this may happen in the case of processor
          ! splitting. Still, we have to write *something*, otherwise
          ! the CDI will not store NetCDF attributes, which, in turn,
          ! is necessary to detect size-0 payload files!
          CALL streamWriteVar(cdiIds%file, cellIndexVarId, dummy_d, 0)
        END IF
        IF (SIZE(me%edgeData_d) > 0) THEN
          CALL streamWriteVar(cdiIds%file, edgeIndexVarId, me%edgeData_d, 0)
          bytesWritten = bytesWritten + INT(SIZE(me%edgeData_d), i8)*8_i8
        ELSE
          CALL streamWriteVar(cdiIds%file, edgeIndexVarId, dummy_d, 0)
        END IF
        IF (SIZE(me%vertData_d) > 0) THEN
          CALL streamWriteVar(cdiIds%file, vertIndexVarId, me%vertData_d, 0)
          bytesWritten = bytesWritten + INT(SIZE(me%vertData_d), i8)*8_i8
        ELSE
          CALL streamWriteVar(cdiIds%file, vertIndexVarId, dummy_d, 0)
        END IF

        IF(timers_level >= 7) CALL timer_stop(timer_write_restart_io)
    END FUNCTION multifilePatchData_openPayloadFile

    !collect the data, optionally writing it OUT immediately
    SUBROUTINE multifilePatchData_collectData(me, opt_fileId, opt_bytesWritten)
        CLASS(t_MultifilePatchData), INTENT(INOUT) :: me
        INTEGER, OPTIONAL, INTENT(IN) :: opt_fileId
        INTEGER(i8), OPTIONAL, INTENT(INOUT) :: opt_bytesWritten

        INTEGER :: curVar, curLevel
        TYPE(t_var_metadata), POINTER :: curInfo
        TYPE(t_ptr_2d), ALLOCATABLE :: dataPointers_d(:)
        TYPE(t_ptr_2d_sp), ALLOCATABLE :: dataPointers_s(:)
        TYPE(t_ptr_2d_int), ALLOCATABLE :: dataPointers_int(:)
        REAL(dp) :: dataDummy_d(1,1)
        REAL(sp) :: dataDummy_s(1,1)
        INTEGER  :: dataDummy_int(1,1)
        REAL(dp), ALLOCATABLE :: buffer_dp(:)
        CHARACTER(*), PARAMETER :: routine = modname//":multifilePatchData_collectData"

        IF(.NOT.ALLOCATED(me%buffers)) CALL finish(routine, "assertion failed, call to createCollectors() missing")

        DO curVar = 1, SIZE(me%varData)
            curInfo => me%varData(curVar)%info
            IF(.NOT.has_valid_time_level(curInfo, me%description%id, me%description%nnew, &
              &                          me%description%nnew_rcf)) CYCLE !no valid time level -> no output

            !on work processes, get the actual data pointers
            IF(my_process_is_work()) THEN
              
              SELECT CASE(curInfo%data_type)
              CASE(REAL_T)
                CALL getLevelPointers(curInfo, me%varData(curVar)%r_ptr, dataPointers_d)
                IF(SIZE(dataPointers_d) /= SIZE(me%buffers(curVar)%p)) CALL finish(routine, "assertion failed")
                !
              CASE(SINGLE_T)
                CALL getLevelPointers(curInfo, me%varData(curVar)%s_ptr, dataPointers_s)
                IF(SIZE(dataPointers_s) /= SIZE(me%buffers(curVar)%p)) CALL finish(routine, "assertion failed")
                !
              CASE(INT_T)
                CALL getLevelPointers(curInfo, me%varData(curVar)%i_ptr, dataPointers_int)
                IF(SIZE(dataPointers_int) /= SIZE(me%buffers(curVar)%p)) CALL finish(routine, "assertion failed")
              CASE DEFAULT
                CALL finish(routine, "Internal error! Variable "//TRIM(curInfo%name))
              END SELECT
            END IF

            DO curLevel = 1, SIZE(me%buffers(curVar)%p)
              SELECT CASE(curInfo%data_type)
              CASE(REAL_T)
                IF(timers_level >= 7) CALL timer_start(timer_write_restart_communication)
                IF(my_process_is_work()) THEN
                  CALL me%collectors(curVar)%p%collectField(dataPointers_d(curLevel)%p, &
                    &                                       me%buffers(curVar)%p(curLevel)%dp)
                ELSE
                  CALL me%collectors(curVar)%p%collectField(dataDummy_d, &
                    &                                       me%buffers(curVar)%p(curLevel)%dp)
                END IF
                IF(timers_level >= 7) CALL timer_stop(timer_write_restart_communication)
                IF(PRESENT(opt_fileId)) THEN
                  IF(timers_level >= 7) CALL timer_start(timer_write_restart_io)
                  IF (SIZE(me%buffers(curVar)%p(curLevel)%dp) > 0) THEN
                    CALL streamWriteVarSlice(opt_fileId, curInfo%cdiVarId, curLevel - 1, me%buffers(curVar)%p(curLevel)%dp, 0)
                  END IF
                  IF(timers_level >= 7) CALL timer_stop(timer_write_restart_io)
                  IF(PRESENT(opt_bytesWritten)) THEN
                    opt_bytesWritten = opt_bytesWritten + INT(SIZE(me%buffers(curVar)%p(curLevel)%dp), i8)*8_i8
                  END IF
                END IF
                !
              CASE(SINGLE_T)
                IF(timers_level >= 7) CALL timer_start(timer_write_restart_communication)
                IF(my_process_is_work()) THEN
                  CALL me%collectors(curVar)%p%collectField(dataPointers_s(curLevel)%p, me%buffers(curVar)%p(curLevel)%sp)
                ELSE
                  CALL me%collectors(curVar)%p%collectField(dataDummy_s, me%buffers(curVar)%p(curLevel)%sp)
                END IF
                IF(timers_level >= 7) CALL timer_stop(timer_write_restart_communication)
                IF(PRESENT(opt_fileId)) THEN
                  IF(timers_level >= 7) CALL timer_start(timer_write_restart_io)
                  IF (SIZE(me%buffers(curVar)%p(curLevel)%sp) > 0) THEN
                    CALL streamWriteVarSliceF(opt_fileId, curInfo%cdiVarId, curLevel - 1, me%buffers(curVar)%p(curLevel)%sp, 0)
                  END IF
                  IF(timers_level >= 7) CALL timer_stop(timer_write_restart_io)
                  IF(PRESENT(opt_bytesWritten)) THEN
                    opt_bytesWritten = opt_bytesWritten + INT(SIZE(me%buffers(curVar)%p(curLevel)%sp), i8)*4_i8
                  END IF
                END IF
                !               
              CASE(INT_T)
                IF(timers_level >= 7) CALL timer_start(timer_write_restart_communication)
                IF(my_process_is_work()) THEN
                  CALL me%collectors(curVar)%p%collectField(dataPointers_int(curLevel)%p, me%buffers(curVar)%p(curLevel)%int)
                ELSE
                  CALL me%collectors(curVar)%p%collectField(dataDummy_int, me%buffers(curVar)%p(curLevel)%int)
                END IF
                IF(timers_level >= 7) CALL timer_stop(timer_write_restart_communication)
                IF(PRESENT(opt_fileId)) THEN
                  IF(timers_level >= 7) CALL timer_start(timer_write_restart_io)
                  ! copy to double precision field
                  IF (SIZE(me%buffers(curVar)%p(curLevel)%int) > 0) THEN
                    ALLOCATE(buffer_dp(SIZE(me%buffers(curVar)%p(curLevel)%int)))
                    buffer_dp(:) = REAL(me%buffers(curVar)%p(curLevel)%int(:), dp)
                    CALL streamWriteVarSlice(opt_fileId, curInfo%cdiVarId, curLevel - 1, buffer_dp, 0)
                    DEALLOCATE(buffer_dp)
                  END IF
                  IF(timers_level >= 7) CALL timer_stop(timer_write_restart_io)
                  IF(PRESENT(opt_bytesWritten)) THEN
                    opt_bytesWritten = opt_bytesWritten + INT(SIZE(me%buffers(curVar)%p(curLevel)%int), i8)*8_i8
                  END IF
                END IF
                !
              CASE DEFAULT
                CALL finish(routine, "Internal error! Variable "//TRIM(curInfo%name))
              END SELECT
            END DO
        END DO
    END SUBROUTINE multifilePatchData_collectData

    !write the data from the buffers to the file.
    !
    !This is for dedicated proc mode only, joint proc mode must USE
    !collectData() to write the data immediately after it is collected
    !because is that case there is only three buffers for all
    !variables AND levels.
    SUBROUTINE multifilePatchData_writeFromBuffer(me, fileId, bytesWritten)
        CLASS(t_MultifilePatchData), INTENT(INOUT) :: me
        INTEGER, VALUE :: fileId
        INTEGER(i8), INTENT(INOUT) :: bytesWritten

        INTEGER :: curVar, curLevel, isize
        TYPE(t_var_metadata), POINTER :: curInfo
        REAL(dp), ALLOCATABLE :: buffer_dp(:)
        CHARACTER(*), PARAMETER :: routine = modname//":multifilePatchData_writeFromBuffer"

        IF(.NOT.ALLOCATED(me%buffers)) THEN
          CALL finish(routine, "assertion failed, call to createCollectors() missing")
        END IF

        IF(.NOT.my_process_is_restart_writer()) THEN
          CALL finish(routine, "assertion failed: must only be called by restart writers");
        END IF
        IF(.NOT.isDedicatedProcMode()) THEN
          CALL finish(routine, "assertion failed: must not be called in joint proc mode");
        END IF

        DO curVar = 1, SIZE(me%varData)
            curInfo => me%varData(curVar)%info

            !no valid time level -> no output:
            IF(.NOT.has_valid_time_level(curInfo, me%description%id, &
              &                          me%description%nnew, me%description%nnew_rcf)) CYCLE

            SELECT CASE(curInfo%data_type)
            CASE(REAL_T)
              DO curLevel = 1, SIZE(me%buffers(curVar)%p)
                IF(timers_level >= 7) CALL timer_start(timer_write_restart_io)
                IF (SIZE(me%buffers(curVar)%p(curLevel)%dp) > 0) THEN
                  CALL streamWriteVarSlice(fileId, curInfo%cdiVarId, curLevel - 1, me%buffers(curVar)%p(curLevel)%dp, 0)
                END IF
                IF(timers_level >= 7) CALL timer_stop(timer_write_restart_io)
                bytesWritten = bytesWritten + INT(SIZE(me%buffers(curVar)%p(curLevel)%dp), i8)*8_i8
              END DO
              !
            CASE(SINGLE_T)
              DO curLevel = 1, SIZE(me%buffers(curVar)%p)
                IF(timers_level >= 7) CALL timer_start(timer_write_restart_io)
                IF (SIZE(me%buffers(curVar)%p(curLevel)%sp) > 0) THEN
                  CALL streamWriteVarSliceF(fileId, curInfo%cdiVarId, curLevel - 1, me%buffers(curVar)%p(curLevel)%sp, 0)
                END IF
                IF(timers_level >= 7) CALL timer_stop(timer_write_restart_io)
                bytesWritten = bytesWritten + INT(SIZE(me%buffers(curVar)%p(curLevel)%sp), i8)*4_i8
              END DO
              !
            CASE(INT_T)
              IF (SIZE(me%buffers(curVar)%p(1)%int) > 0) THEN
                isize = SIZE(me%buffers(curVar)%p(1)%int)
                ALLOCATE(buffer_dp(isize))
                DO curLevel = 1, SIZE(me%buffers(curVar)%p)
                  IF(timers_level >= 7) CALL timer_start(timer_write_restart_io)
                  buffer_dp(:) = REAL(me%buffers(curVar)%p(curLevel)%int(:), dp)
                  CALL streamWriteVarSlice(fileId, curInfo%cdiVarId, curLevel - 1, buffer_dp, 0)
                  IF(timers_level >= 7) CALL timer_stop(timer_write_restart_io)
                  bytesWritten = bytesWritten + INT(SIZE(me%buffers(curVar)%p(curLevel)%dp), i8)*8_i8
                END DO
                DEALLOCATE(buffer_dp)
              END IF
              !
            CASE DEFAULT
              CALL finish(routine, "Internal error! Variable "//TRIM(curInfo%name))
            END SELECT
        END DO
    END SUBROUTINE multifilePatchData_writeFromBuffer

    SUBROUTINE multifilePatchData_writeFile(me, restartAttributes, restartArgs, lIsWriteProcess)
        CLASS(t_MultifilePatchData), INTENT(INOUT) :: me
        TYPE(t_RestartAttributeList), INTENT(INOUT) :: restartAttributes
        TYPE(t_restart_args), INTENT(IN) :: restartArgs
        LOGICAL, VALUE :: lIsWriteProcess

        CHARACTER(*), PARAMETER :: routine = modname//":multifilePatchData_writeFile"

        CALL finish(routine, "assertion failed: this implementation must not be called")
    END SUBROUTINE multifilePatchData_writeFile

    SUBROUTINE multifilePatchData_destruct(me)
        CLASS(t_MultifilePatchData), INTENT(INOUT) :: me

        INTEGER :: i, j
        CHARACTER(*), PARAMETER :: routine = modname//":multifilePatchData_destruct"

        IF(.NOT.ALLOCATED(me%buffers)) CALL finish(routine, "assertion failed, call to createCollectors() missing")

        DEALLOCATE(me%cellData_d)
        DEALLOCATE(me%cellData_s)
        DEALLOCATE(me%edgeData_d)
        DEALLOCATE(me%edgeData_s)
        DEALLOCATE(me%vertData_d)
        DEALLOCATE(me%vertData_s)

        DO i = 1, SIZE(me%buffers)
            IF(my_process_is_restart()) THEN    ! we only have separate allocations on dedicated restart procs
                DO j = 1, SIZE(me%buffers(i)%p)
                    IF(ASSOCIATED(me%buffers(i)%p(j)%dp)) DEALLOCATE(me%buffers(i)%p(j)%dp)
                    IF(ASSOCIATED(me%buffers(i)%p(j)%sp)) DEALLOCATE(me%buffers(i)%p(j)%sp)
                    IF(ASSOCIATED(me%buffers(i)%p(j)%int)) DEALLOCATE(me%buffers(i)%p(j)%int)
                END DO
            END IF
            DEALLOCATE(me%buffers(i)%p)
        END DO
        DEALLOCATE(me%buffers)

        CALL me%cellCollector%destruct()
        CALL me%edgeCollector%destruct()
        CALL me%vertCollector%destruct()

        DEALLOCATE(me%collectors)

        CALL me%t_RestartPatchData%destruct()
    END SUBROUTINE multifilePatchData_destruct

END MODULE mo_multifile_restart_patch_data
