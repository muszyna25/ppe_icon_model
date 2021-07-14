!> The multifile version of t_RestartPatchData.
!!
!! Initial implementation: Nathanael Hübbe
!! 2018-08: Major revision / revamp / refactoring : Harald Braun (Atos SE)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
#include "omp_definitions.inc"
MODULE mo_multifile_restart_patch_data
  USE mo_cdi,                         ONLY: FILETYPE_NC4, DATATYPE_FLT64, streamWriteVar, streamWriteVarF,     &
    &                                       streamWriteVarSlice, streamWriteVarSliceF
  USE mo_zaxis_type,                  ONLY: ZA_SURFACE
  USE mo_cdi_ids,                     ONLY: t_CdiIds
  USE mtime,                          ONLY: datetime
  USE mo_exception,                   ONLY: finish
  USE mo_impl_constants,              ONLY: SINGLE_T, REAL_T, INT_T
  USE mo_cdi_constants, ONLY: GRID_UNSTRUCTURED_CELL, GRID_UNSTRUCTURED_EDGE, GRID_UNSTRUCTURED_VERT
  USE mo_kind,                        ONLY: dp, sp, i8
  USE mo_mpi,                         ONLY: my_process_is_work
  USE mo_multifile_restart_collector, ONLY: t_MultifileRestartCollector
  USE mo_multifile_restart_util,      ONLY: multifilePayloadPath, isAsync, vNames_glbIdx, commonBuf_t, &
    & typeMax, typeID, typeMap, typeByte
  USE mo_restart_patch_data,          ONLY: t_RestartPatchData
  USE mo_restart_patch_description,   ONLY: t_restart_patch_description
  USE mo_restart_var_data,            ONLY: has_valid_time_level, get_var_3d_ptr
  USE mo_timer,                       ONLY: timer_start, timer_stop, timer_write_restart_io,                   &
    &                                       timer_write_restart_communication, timers_level,                   &
    &                                       timer_write_restart_setup, timer_write_restart_wait
  USE mo_var_metadata_types,          ONLY: t_var_metadata
  USE mo_var,                         ONLY: t_var_ptr
  USE mo_parallel_config,             ONLY: restart_chunk_size
  USE mo_var_list_register_utils,     ONLY: vlr_select_restart_vars

  IMPLICIT NONE
  PRIVATE
  
  PUBLIC :: t_MultifilePatchData

  TYPE, EXTENDS(t_RestartPatchData) :: t_MultifilePatchData
    INTEGER :: cnkLvs
    LOGICAL :: shortcut
    TYPE(t_MultifileRestartCollector) :: coll
  CONTAINS
    PROCEDURE :: construct => multifilePatchData_construct
    PROCEDURE :: createCollectors => multifilePatchData_createCollectors
    PROCEDURE :: openPayloadFile => multifilePatchData_openPayloadFile  ! restart writers only
    PROCEDURE :: start_local_access  => multifilePatchData_start_local_access
    PROCEDURE :: start_remote_access => multifilePatchData_start_remote_access
    PROCEDURE :: exposeData  => multifilePatchData_exposeData    
    PROCEDURE :: collectData => multifilePatchData_collectData
    PROCEDURE :: destruct => multifilePatchData_destruct    ! override
    PROCEDURE :: select_vars => multifilePatchData_select_vars
    PROCEDURE :: writeData => multifilePatchData_writeData
  END TYPE t_MultifilePatchData
  
  CHARACTER(*), PARAMETER :: modname = "mo_multifile_restart_patch_data"
  REAL(dp), TARGET :: dummy_d(0)

CONTAINS

  SUBROUTINE multifilePatchData_writeData(me, file_handle)
    CLASS(t_MultifilePatchData), INTENT(INOUT) :: me
    INTEGER,                     INTENT(IN)    :: file_handle

    CALL finish(modname//"writeData", "not implemented!")
  END SUBROUTINE multifilePatchData_writeData

  SUBROUTINE multifilePatchData_construct(me, modelType, jg)
    CLASS(t_MultifilePatchData), INTENT(INOUT) :: me
    CHARACTER(*),                INTENT(IN)    :: modelType
    INTEGER,                     INTENT(IN)    :: jg

    CALL me%description%init(jg)
    CALL vlr_select_restart_vars(me%varData, jg, modelType, me%restartType)
    IF (isAsync()) CALL me%description%transferToRestart()
  END SUBROUTINE multifilePatchData_construct

  SUBROUTINE multifilePatchData_createCollectors(me, wRnk, srcRnks, lactive)
    CLASS(t_MultifilePatchData), INTENT(INOUT), TARGET :: me
    INTEGER,                     INTENT(IN) :: wRnk, srcRnks(:)
    LOGICAL,                     INTENT(IN) :: lactive
    CHARACTER(*), PARAMETER :: routine = modname//":multifilePatchData_createCollectors"
    INTEGER :: curVar, jg, nLev, iType, nVar, i
    INTEGER(i8)                     :: iOffset(typeMax)
    TYPE(t_var_metadata), POINTER   :: curInfo
    TYPE(t_var_ptr), ALLOCATABLE :: varReordered(:)

    IF (timers_level >= 7) CALL timer_start(timer_write_restart_setup)
    jg = me%description%id
    nVar = SIZE(me%varData)
    CALL me%coll%construct(jg, nVar, wRnk, srcRnks, lactive)
    ioffset(:) = 0_i8
    ALLOCATE(varReordered(nVar))
    i = 0
    DO iType = 1, typeMax
      DO curVar = 1, nVar
        curInfo => me%varData(curVar)%p%info
        IF (curInfo%data_type .EQ. iType) THEN
          i = i + 1
          varReordered(i) = me%varData(curVar)
        END IF
      END DO
    END DO
    CALL MOVE_ALLOC(varReordered, me%varData)
    IF (i .NE. nVar) CALL finish(routine, "inconsistency!!!")
    me%shortcut = .false.
    IF (SIZE(srcRnks) .EQ. 1) &
      me%shortcut = .NOT.isAsync() .AND. wRnk .EQ. srcRnks(1)
    me%cnkLvs = 1
    DO curVar = 1, nVar
      curInfo => me%varData(curVar)%p%info
      iType = curInfo%data_type
      IF (curInfo%ndims > 3) CALL finish(routine, "ndims > 3 is not supported")
      nlev = MERGE(curInfo%used_dimensions(2), 1, curInfo%ndims == 3)
      CALL me%coll%defVar(curVar, nLev, iType, curInfo%hgrid, iOffset)
      me%cnkLvs = MAX(me%cnkLvs, nLev)
    END DO
    me%cnkLvs = MERGE(restart_chunk_size, me%cnkLvs, restart_chunk_size .GT. 0)
    CALL me%coll%init_win(iOffset, me%shortcut)
    IF (timers_level >= 7) CALL timer_stop(timer_write_restart_setup)
  END SUBROUTINE multifilePatchData_createCollectors

  FUNCTION multifilePatchData_openPayloadFile(me, filename, ifile, this_datetime, bytesWritten) RESULT(cdiIds)
    TYPE(t_CdiIds)                                     :: cdiIds
    CLASS(t_MultifilePatchData), TARGET, INTENT(INOUT) :: me
    CHARACTER(*),                        INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ifile
    TYPE(datetime), POINTER,             INTENT(IN)    :: this_datetime
    INTEGER(i8),                         INTENT(INOUT) :: bytesWritten
    CHARACTER(:), ALLOCATABLE                          :: effectiveFilename
    CLASS(t_restart_patch_description), POINTER        :: desc
    INTEGER :: iVarId(3), curVar, iG, sub(3), lCnt, iV
    REAL(dp), POINTER :: buf_d(:)
    LOGICAL, DIMENSION(SIZE(me%varData)) :: vWrDone
    TYPE(t_var_metadata), POINTER   :: curInfo
    INTEGER, ALLOCATABLE :: vWrNow(:)

    IF (timers_level >= 7) CALL timer_start(timer_write_restart_io)
    desc => me%description
    CALL multifilePayloadPath(filename, desc%id, ifile, effectiveFilename)
    CALL cdiIds%init()
    IF (ALLOCATED(desc%opt_pvct)) THEN
      CALL cdiIds%openRestartAndCreateIds(effectiveFilename, FILETYPE_NC4,  &
        & me%coll%idx(1)%nRecv, me%coll%idx(2)%nRecv, me%coll%idx(3)%nRecv, &
        & desc%cell_type, desc%v_grid_defs(1:desc%v_grid_count), desc%opt_pvct)
    ELSE
      CALL cdiIds%openRestartAndCreateIds(effectiveFilename, FILETYPE_NC4,  &
        & me%coll%idx(1)%nRecv, me%coll%idx(2)%nRecv, me%coll%idx(3)%nRecv, &
        & desc%cell_type, desc%v_grid_defs(1:desc%v_grid_count))
    END IF
    iVarId(1) = cdiIds%defineVariableSimple(GRID_UNSTRUCTURED_CELL, ZA_SURFACE, DATATYPE_FLT64, vNames_glbIdx(1))
    iVarId(2) = cdiIds%defineVariableSimple(GRID_UNSTRUCTURED_VERT, ZA_SURFACE, DATATYPE_FLT64, vNames_glbIdx(2))
    iVarId(3) = cdiIds%defineVariableSimple(GRID_UNSTRUCTURED_EDGE, ZA_SURFACE, DATATYPE_FLT64, vNames_glbIdx(3))
    vWrDone(:) = .false.
    sub(:) = 0
    DO WHILE (ANY(.NOT.vWrDone(:)))
      CALL me%select_vars(vWrDone, vWrNow, sub)
      DO iV = 1, SIZE(vWrNow)
        curVar = vWrNow(iV)
        curInfo => me%varData(curVar)%p%info
        lCnt = MERGE(curInfo%used_dimensions(2), 1, curInfo%ndims == 3)
        IF (sub(3) .EQ. curVar) THEN
          IF (sub(2) .EQ. lCnt) sub(:) = 0
          IF (sub(3) .NE. 0) CYCLE
        END IF
        CALL cdiIds%defineVariable(me%varData(curVar)%p%info)
      END DO
    END DO
    IF (ALLOCATED(vWrNow)) DEALLOCATE(vWrNow)
    CALL cdiIds%finalizeVlist(this_datetime)
    DO iG = 1, 3
      buf_d => dummy_d
      IF (me%coll%idx(iG)%nRecv .GT. 0) buf_d => me%coll%glb_idx(iG)%p
      CALL streamWriteVar(cdiIds%fHndl, iVarId(iG), buf_d, 0)
      bytesWritten = bytesWritten + INT(SIZE(buf_d), i8) * 8_i8
    END DO
    IF (timers_level >= 7) CALL timer_stop(timer_write_restart_io)
  END FUNCTION multifilePatchData_openPayloadFile

  SUBROUTINE multifilePatchData_start_local_access(me)
    CLASS(t_MultifilePatchData), INTENT(INOUT) :: me

    IF (timers_level >= 7) CALL timer_start(timer_write_restart_wait)
    CALL me%coll%local_access()
    IF (timers_level >= 7) CALL timer_stop(timer_write_restart_wait)
  END SUBROUTINE multifilePatchData_start_local_access

  SUBROUTINE multifilePatchData_start_remote_access(me)
    CLASS(t_MultifilePatchData), INTENT(INOUT) :: me

    IF (timers_level >= 7) CALL timer_start(timer_write_restart_wait)
    CALL me%coll%remote_access()
    IF (timers_level >= 7) CALL timer_stop(timer_write_restart_wait)
  END SUBROUTINE multifilePatchData_start_remote_access

  SUBROUTINE multifilePatchData_select_vars(me, vWrDone, vWrNow, sub)
    CLASS(t_MultifilePatchData), INTENT(IN) :: me
     CHARACTER(*), PARAMETER :: routine = modname//":multifilePatchData_select_vars"
    LOGICAL, INTENT(INOUT) :: vWrDone(SIZE(me%varData))
    INTEGER, ALLOCATABLE, INTENT(OUT) :: vWrNow(:)
    INTEGER, INTENT(INOUT) :: sub(3)
    TYPE(t_var_metadata), POINTER   :: curInfo
    INTEGER :: lFree, iV, lN, tFrst, nW
    LOGICAL :: l_vWrNow(SIZE(me%varData))

    l_vWrNow(:) = .false.
    lFree = me%cnkLvs
    tFrst = -1
    nW = 0
    DO iV = 1, SIZE(me%varData)
      IF (lFree .EQ. 0) EXIT
      IF (vWrDone(iV)) CYCLE
      curInfo => me%varData(iV)%p%info
      IF (.NOT.has_valid_time_level(curInfo, me%description%id, &
        &        me%description%nnew, me%description%nnew_rcf)) THEN
        vWrDone(iV) = .true. ! nothing to do
        CYCLE
      END IF
      IF (tFrst .EQ. -1) tFrst = curInfo%data_type
      IF (.NOT.ANY(tFrst .EQ. typeID(:))) CALL finish(routine, &
        & "data type not recognized! Variable" // TRIM(curInfo%name))
      IF (tFrst .NE. curInfo%data_type) EXIT
      lN = MERGE(curInfo%used_dimensions(2), 1, curInfo%ndims .EQ. 3)
      IF (lFree .GE. lN .OR. me%shortcut) THEN
        IF (.NOT.me%shortcut) lFree = lFree - lN
        vWrDone(iV) = .true.
        IF (lN .NE. 0) l_vWrNow(iV) = .true.
      ELSE IF (me%cnkLvs .LT. lN .AND. (sub(3) .EQ. iV .OR. sub(3) .EQ. 0)) THEN
        l_vWrNow(iV) = .true.
        sub(1) = sub(2)
        sub(2) = MIN(sub(2) + me%cnkLvs, lN)
        sub(3) = iV
        lFree = lFree - sub(2) + sub(1)
        IF (sub(2) .EQ. lN) vWrDone(iV) = .true.
      END IF
      IF (l_vWrNow(iV)) nW = nW + 1
    END DO
    IF (ALLOCATED(vWrNow)) DEALLOCATE(vWrNow)
    ALLOCATE(vWrNow(nW))
    IF (nW .GT. 0) vWrNow(:) = PACK([(iV, iV = 1, SIZE(me%varData))], l_vWrNow(:))
  END SUBROUTINE multifilePatchData_select_vars

  SUBROUTINE multifilePatchData_exposeData(me)
    CLASS(t_MultifilePatchData), INTENT(INOUT) :: me
    CHARACTER(*), PARAMETER :: routine = modname//":multifilePatchData_exposeData"
    TYPE(t_var_metadata), POINTER   :: curInfo
    INTEGER :: iV, sub(3), lCnt, curVar
    INTEGER(i8) :: offset(3)
    REAL(KIND=dp), POINTER :: r_ptr_3d(:,:,:)
    REAL(KIND=sp), POINTER :: s_ptr_3d(:,:,:)
    INTEGER, POINTER :: i_ptr_3d(:,:,:)
    LOGICAL, DIMENSION(SIZE(me%varData)) :: vWrDone
    INTEGER, ALLOCATABLE :: vWrNow(:)

    offset(:) = 0_i8
    IF (.NOT. my_process_is_work()) CALL finish(routine, "only workers allowed")
    IF (timers_level >= 7) CALL timer_start(timer_write_restart_io)
    vWrDone(:) = .false.
    sub(:) = 0
    DO WHILE (ANY(.NOT.vWrDone(:)))
      CALL me%select_vars(vWrDone, vWrNow, sub)
      DO iV = 1, SIZE(vWrNow)
        curVar = vWrNow(iV)
        curInfo => me%varData(curVar)%p%info
        lCnt = MERGE(curInfo%used_dimensions(2), 1, curInfo%ndims == 3)
        IF (sub(3) .EQ. curVar) THEN
          IF (sub(2) .EQ. lCnt) sub(:) = 0
          IF (sub(3) .NE. 0) CYCLE
        END IF
        SELECT CASE (curInfo%data_type)
        CASE(REAL_T)
          CALL get_var_3d_ptr(me%varData(curVar)%p, r_ptr_3d)
          CALL me%coll%sendField(curVar, r_ptr_3d, offset(1))
        CASE(SINGLE_T)
          CALL get_var_3d_ptr(me%varData(curVar)%p, s_ptr_3d)
          CALL me%coll%sendField(curVar, s_ptr_3d, offset(2))
        CASE(INT_T)
          CALL get_var_3d_ptr(me%varData(curVar)%p, i_ptr_3d)
          CALL me%coll%sendField(curVar, i_ptr_3d, offset(3))
        END SELECT
      END DO
    END DO
    IF (ALLOCATED(vWrNow)) DEALLOCATE(vWrNow)
    IF (timers_level >= 7) CALL timer_stop(timer_write_restart_io)
  END SUBROUTINE multifilePatchData_exposeData

  SUBROUTINE multifilePatchData_collectData(me, fId, bytesWritten)
    CLASS(t_MultifilePatchData), INTENT(INOUT) :: me
    INTEGER,                     INTENT(IN)    :: fId
    INTEGER(i8),                 INTENT(INOUT) :: bytesWritten
    TYPE(t_var_metadata), POINTER   :: curInfo
    REAL(KIND=dp), POINTER :: buf_d(:)
    INTEGER :: lOff, vOff, lInc, i, iV, iL, lCnt, sub(3), curVar
    TYPE(commonBuf_t) :: rBuf
    INTEGER, ALLOCATABLE :: vSize(:), vWrNow(:)
    LOGICAL, DIMENSION(SIZE(me%varData)) :: vWrDone
    INTEGER(KIND=i8), ALLOCATABLE :: srcOff(:,:)
 
    DO i = 1, 3
      IF (.NOT.ASSOCIATED(me%coll%glb_idx(i)%p)) RETURN
      IF (SIZE(me%coll%glb_idx(i)%p) .LE. 0) RETURN
    END DO
    buf_d => dummy_d
    vWrDone(:) = .false.
    sub(:) = 0
    DO WHILE (ANY(.NOT.vWrDone(:)))
      CALL me%select_vars(vWrDone, vWrNow, sub)
      IF (timers_level >= 7) CALL timer_start(timer_write_restart_communication)
      CALL me%coll%fetch(vWrNow, sub, vSize, rBuf, srcOff)
      IF (timers_level >= 7) CALL timer_stop(timer_write_restart_communication)
      IF (timers_level >= 7) CALL timer_start(timer_write_restart_io)
      vOff = 0
      DO iV = 1, SIZE(vWrNow)
        curVar = vWrNow(iV)
        curInfo => me%varData(curVar)%p%info
        lCnt = MERGE(curInfo%used_dimensions(2), 1, curInfo%ndims == 3)
        IF (sub(3) .EQ. curVar) THEN
          lOff = vOff
          lInc = vSize(iV) / (sub(2) - sub(1))
          DO iL = sub(1), sub(2) - 1
            SELECT CASE(curInfo%data_type)
            CASE(REAL_T)
              CALL streamWriteVarSlice(fId, curInfo%cdiVarId, iL, rBuf%d(1+lOff:lOff+lInc), 0)
            CASE(SINGLE_T)
              CALL streamWriteVarSliceF(fId, curInfo%cdiVarId, iL, rBuf%s(1+lOff:lOff+lInc), 0)
            CASE(INT_T)
              IF (lInc .GT. SIZE(buf_d)) THEN
                IF (SIZE(buf_d) .GT. 0) DEALLOCATE(buf_d)
                ALLOCATE(buf_d(lInc))
              END IF
              !ICON_OMP PARALLEL DO SCHEDULE(STATIC)
              DO i = 1, lInc
                buf_d(i) = REAL(rBuf%i(lOff+i), dp)
              END DO
              CALL streamWriteVarSlice(fId, curInfo%cdiVarId, iL, buf_d(1:lInc), 0)
            END SELECT
            lOff = lOff + lInc
          END DO
          IF (sub(2) .EQ. lCnt) sub(:) = 0
        ELSE
          SELECT CASE(curInfo%data_type)
          CASE(REAL_T)
            CALL streamWriteVar(fId, curInfo%cdiVarId, rBuf%d(1+vOff:vOff+vSize(iV)), 0)
          CASE(SINGLE_T)
            CALL streamWriteVarF(fId, curInfo%cdiVarId, rBuf%s(1+vOff:vOff+vSize(iV)), 0)
          CASE(INT_T)
            IF (vSize(iV) .GT. SIZE(buf_d)) THEN
              IF (SIZE(buf_d) .GT. 0) DEALLOCATE(buf_d)
              ALLOCATE(buf_d(vSize(iV)))
            END IF
            !ICON_OMP PARALLEL DO SCHEDULE(STATIC)
            DO i = 1, vSize(iV)
              buf_d(i) = REAL(rBuf%i(vOff+i), dp)
            END DO
            CALL streamWriteVar(fId, curInfo%cdiVarId, buf_d(1:vSize(iV)), 0)
          END SELECT
        END IF
        vOff = vOff + vSize(iV)
        bytesWritten = bytesWritten + INT(vSize(iV), i8) * typeByte(typeMap(curInfo%data_type))
      END DO
      IF (timers_level >= 7) CALL timer_stop(timer_write_restart_io)
    END DO
    IF(ALLOCATED(srcOff)) DEALLOCATE(srcOff)
    IF(SIZE(buf_d) .GT. 0) DEALLOCATE(buf_d)
    IF(ALLOCATED(vSize)) DEALLOCATE(vSize, vWrNow)
  END SUBROUTINE multifilePatchData_collectData

  SUBROUTINE multifilePatchData_destruct(me)
    CLASS(t_MultifilePatchData), INTENT(INOUT) :: me

    CALL me%coll%finalize()
    IF (ALLOCATED(me%varData)) DEALLOCATE(me%varData)
  END SUBROUTINE multifilePatchData_destruct

END MODULE mo_multifile_restart_patch_data
