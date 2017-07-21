!> Module for reading single file restart files
!!
!! Note: The asynchronous implementation of the restart output can be
!!       found in the module "mo_async_restart"
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.

MODULE mo_load_singlefile_restart
    USE mo_cdi,                ONLY: vlistInqTaxis, taxisInqVdate, taxisInqVtime, streamOpenRead, &
      &                              streamInqVlist, vlistNvars, vlistInqVarName, streamClose,    &
      &                              streamReadVarSlice, streamReadVarSliceF, vlistInqVarGrid, gridInqSize
    USE mo_cdi_constants,      ONLY: GRID_UNSTRUCTURED_CELL, GRID_UNSTRUCTURED_VERT, GRID_UNSTRUCTURED_EDGE
    USE mo_communication,      ONLY: t_ScatterPattern
    USE mo_exception,          ONLY: message, warning, finish
    USE mo_fortran_tools,      ONLY: t_ptr_2d, t_ptr_2d_sp, t_ptr_2d_int, alloc
    USE mo_impl_constants,     ONLY: SUCCESS, SINGLE_T, REAL_T, INT_T
    USE mo_kind,               ONLY: dp, sp
    USE mo_model_domain,       ONLY: t_patch
    USE mo_mpi,                ONLY: my_process_is_mpi_workroot, p_bcast, p_comm_work
    USE mo_restart_attributes, ONLY: t_RestartAttributeList, getAttributesForRestarting
    USE mo_restart_util,       ONLY: restartSymlinkName
    USE mo_restart_var_data,   ONLY: t_restartVarData, getLevelPointers
    USE mo_timer,              ONLY: timer_start, timer_stop, timer_load_restart_io, &
      &                              timer_load_restart_communication, timers_level
    USE mo_util_string,        ONLY: int2string
    USE mo_var_metadata_types, ONLY: t_var_metadata

    IMPLICIT NONE

    PUBLIC :: singlefileCheckRestartFiles, singlefileReadPatch

    CHARACTER(*), PARAMETER :: modname = "mo_load_singlefile_restart"

CONTAINS

  ! This IS just an existence check for the different patch files, a
  ! warning IS generated each domain for which no corresponding
  ! restart file IS found.
  !
  ! XXX: The code that I found READ the restart attributes from all
  ! the files, *appending* them to the list of attributes. Since all
  ! restart files contain the same attributes, this ONLY served to
  ! duplicate them. The current code just ignores the attributes IN
  ! all but the first file, just like the original code ignored the
  ! namelists IN all but the first file.  However, it might be a good
  ! idea to add some consistency checking on the attributes/namelists
  ! of the other files.
  SUBROUTINE singlefileCheckRestartFiles(modelType)
    CHARACTER(*), INTENT(IN) :: modelType

    INTEGER :: jg, n_dom
    LOGICAL :: lexists
    CHARACTER(:), ALLOCATABLE :: filename
    TYPE(t_RestartAttributeList), POINTER :: restartAttributes
    CHARACTER(*), PARAMETER :: routine = modname//":singlefileCheckRestartFiles"

    ! since we do not know about the total number of domains yet,
    ! we have to ask the restart file for this information:
    restartAttributes => getAttributesForRestarting()
    n_dom = restartAttributes%getInteger('n_dom')

    ! check whether we have all the restart files we expect
    DO jg = 1, n_dom
        filename = restartSymlinkName(modelType, jg, n_dom)
        INQUIRE(file = filename, exist = lexists)
        IF (lexists) THEN
            CALL message(routine, "found restart file at '"//filename//"'")
        ELSE
            CALL warning(routine, 'domain '//TRIM(int2string(jg))//' is not active at restart time')
        END IF
    END DO
  END SUBROUTINE singlefileCheckRestartFiles

  FUNCTION vtimeString(vlistId) RESULT(resultVar)
    CHARACTER(:), ALLOCATABLE :: resultVar
    INTEGER, VALUE :: vlistId

    INTEGER :: taxisId

    taxisId = vlistInqTaxis(vlistId)
    resultVar = TRIM(int2string(taxisInqVdate(taxisId)))//"T"//&
              & TRIM(int2string(taxisInqVtime(taxisId)))//"Z"
  END FUNCTION vtimeString

  FUNCTION getScatterPattern(p_patch, info) RESULT(resultVar)
    CLASS(t_ScatterPattern), POINTER :: resultVar
    TYPE(t_patch), TARGET, INTENT(IN) :: p_patch
    TYPE(t_var_metadata), INTENT(IN) :: info

    CHARACTER(*), PARAMETER :: routine = modname//":getScatterPattern"

    SELECT CASE(info%hgrid)
        CASE (GRID_UNSTRUCTURED_CELL)
            resultVar => p_patch%comm_pat_scatter_c
        CASE (GRID_UNSTRUCTURED_VERT)
            resultVar => p_patch%comm_pat_scatter_v
        CASE (GRID_UNSTRUCTURED_EDGE)
            resultVar => p_patch%comm_pat_scatter_e
        CASE default
            CALL finish(routine, 'unknown grid type')
    END SELECT
  END FUNCTION getScatterPattern

  SUBROUTINE singlefileReadPatch(varData, modelType, p_patch, opt_ndom)
    TYPE(t_restartVarData), INTENT(INOUT) :: varData(:)
    CHARACTER(*), INTENT(IN) :: modelType
    TYPE(t_patch), TARGET, INTENT(IN) :: p_patch
    INTEGER, OPTIONAL, INTENT(IN) :: opt_ndom

    TYPE(t_ptr_2d), ALLOCATABLE :: levelPointers_d(:)
    TYPE(t_ptr_2d_sp), ALLOCATABLE :: levelPointers_s(:)
    TYPE(t_ptr_2d_int), ALLOCATABLE :: levelPointers_i(:)
    CHARACTER(len=80) :: NAME
    CHARACTER(:), ALLOCATABLE :: restart_filename
    INTEGER :: fileID, vlistID, gridID, varID, nvars, level, trash, varIndex, levelCount
    REAL(dp), ALLOCATABLE :: r1d_d(:)
    REAL(sp), ALLOCATABLE :: r1d_s(:)
    INTEGER, ALLOCATABLE :: i1d(:)
    CLASS(t_ScatterPattern), POINTER :: scatter_pattern
    CHARACTER(*), PARAMETER :: routine = modname//":singlefileReadPatch"

    CALL alloc(r1d_d, 1)  !dummy allocation so that there IS no process without a valid allocation
    CALL alloc(r1d_s, 1)  !dummy allocation so that there IS no process without a valid allocation
    CALL alloc(i1d,   1)  !dummy allocation so that there IS no process without a valid allocation

    restart_filename = restartSymlinkName(modelType, p_patch%id, opt_ndom)

    IF(my_process_is_mpi_workroot()) THEN
        WRITE(0,*) "streamOpenRead ", TRIM(restart_filename)

        fileID  = streamOpenRead(restart_filename)
        IF(fileID < 0) CALL finish(routine, "could not open file '"//TRIM(restart_filename)//"'")
        vlistID = streamInqVlist(fileID)

        WRITE(0,*) routine//': Read restart for: '//vtimeString(vlistId)//' from '//TRIM(restart_filename)

        nvars = vlistNvars(vlistID)
    END IF
    CALL p_bcast(nvars, 0, comm=p_comm_work)

    DO varID = 0, nvars - 1
        IF(my_process_is_mpi_workroot()) CALL vlistInqVarName(vlistID, varID, NAME)
        CALL p_bcast(NAME, 0, comm=p_comm_work)

        !lookup the varData corresponding to this variable AND get the corresponding DATA pointers
        DO varIndex = 1, SIZE(varData)
            IF(NAME == varData(varIndex)%info%NAME) EXIT
        END DO
        IF(varIndex > SIZE(varData)) THEN
            CALL warning(routine, "variable '"//TRIM(NAME)//"' from restart file not found in the list of restart variables")
            CYCLE
        END IF

        SELECT CASE(varData(varIndex)%getDatatype())
        CASE(REAL_T)
          CALL getLevelPointers(varData(varIndex)%info, varData(varIndex)%r_ptr, levelPointers_d)
          levelCount = SIZE(levelPointers_d)
          !
        CASE(SINGLE_T)
          CALL getLevelPointers(varData(varIndex)%info, varData(varIndex)%s_ptr, levelPointers_s)
          levelCount = SIZE(levelPointers_s)
          !
        CASE(INT_T)
          ! Read-in of INTEGER fields: we read them as
          ! REAL-valued arrays and transform them using NINT.
          CALL getLevelPointers(varData(varIndex)%info, varData(varIndex)%i_ptr, levelPointers_i)
          levelCount = SIZE(levelPointers_i)
        CASE DEFAULT
          CALL finish(routine, "Internal error! Variable "//TRIM(varData(varIndex)%info%name))
        END SELECT

        !lookup the correct scatter pattern
        scatter_pattern => getScatterPattern(p_patch, varData(varIndex)%info)

        !READ IN the DATA AND distribute it to the processes
        DO level = 1, levelCount
          IF(my_process_is_mpi_workroot()) THEN
            gridID  = vlistInqVarGrid(vlistID, varID)
            SELECT CASE(varData(varIndex)%getDatatype())
            CASE(REAL_T)
              CALL alloc(r1d_s, 1)
              CALL alloc(r1d_d, gridInqSize(gridID))
              CALL alloc(i1d, 1)
              IF(timers_level >= 7) CALL timer_start(timer_load_restart_io)
              CALL streamReadVarSlice(fileId, varId, level - 1, r1d_d, trash)
              IF(timers_level >= 7) CALL timer_stop(timer_load_restart_io)
              !
            CASE(SINGLE_T)
              CALL alloc(r1d_d, 1)
              CALL alloc(r1d_s, gridInqSize(gridID))
              CALL alloc(i1d, 1)
              IF(timers_level >= 7) CALL timer_start(timer_load_restart_io)
              CALL streamReadVarSliceF(fileId, varId, level - 1, r1d_s, trash)
              IF(timers_level >= 7) CALL timer_stop(timer_load_restart_io)
              !
            CASE(INT_T)
              CALL alloc(r1d_s, 1)
              CALL alloc(r1d_d, gridInqSize(gridID))
              CALL alloc(i1d, gridInqSize(gridID))
              IF(timers_level >= 7) CALL timer_start(timer_load_restart_io)
              CALL streamReadVarSlice(fileId, varId, level - 1, r1d_d, trash)
              i1d(:) = INT(r1d_d)
              IF(timers_level >= 7) CALL timer_stop(timer_load_restart_io)
              !
            CASE DEFAULT
              CALL finish(routine, "Internal error! Variable "//TRIM(varData(varIndex)%info%name))
            END SELECT
          END IF
          IF(timers_level >= 7) CALL timer_start(timer_load_restart_communication)

          SELECT CASE(varData(varIndex)%getDatatype())
          CASE(REAL_T)
            CALL scatter_pattern%distribute(r1d_d, levelPointers_d(level)%p, .FALSE.)
            !
          CASE(SINGLE_T)
            CALL scatter_pattern%distribute(r1d_s, levelPointers_s(level)%p, .FALSE.)
            !
          CASE(INT_T)
            CALL scatter_pattern%distribute(i1d, levelPointers_i(level)%p, .FALSE.)
            !
          CASE DEFAULT
            CALL finish(routine, "Internal error! Variable "//TRIM(varData(varIndex)%info%name))
          END SELECT

          IF(timers_level >= 7) CALL timer_stop(timer_load_restart_communication)
        END DO

        IF (my_process_is_mpi_workroot()) WRITE (0,*) ' ... read '//TRIM(NAME)
    END DO

    IF (my_process_is_mpi_workroot())  CALL streamClose(fileID)

    DEALLOCATE (r1d_d)
    DEALLOCATE (r1d_s)
    DEALLOCATE (i1d)
  END SUBROUTINE singlefileReadPatch

END MODULE mo_load_singlefile_restart
