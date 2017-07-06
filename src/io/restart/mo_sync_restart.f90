!> Module for writing restart files (synchronously)
!!
!! Note: The asynchronous implementation of the restart output can be
!!       found in the module "mo_async_restart"
!!
!! Initial implementation: L. Kornblueh
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!! --------------------------------------------------------------------------------
!!
!! Generated files:
!! ================
!!
!! 1. For each domain 1, ..., n_dom, and for each restart output time step:
!!
!!    Restart data file
!!    -----------------
!!      "<gridfile>_restart_<modeltype>_<timestamp>.nc",     e.g.
!!      "iconR2B06_DOM01_restart_atm_20110101T001200Z.nc"    (NetCDF format)
!!
!!    This filename can be customized using the namelist parameter
!!    "mo_run_nml/restart_filename".
!!
!!    This file contains
!!    -  data
!!    -  namelists
!!    -  several attributes
!!
!!    Note:
!!    -  We read the namelists only once and assume that these
!!       are identical for all domains.
!!    -  Since we do not know about the total number of domains at startup,
!!       we have to ask the current restart file for the attribute "n_dom"
!!
!! 2. For each domain 1, ..., n_dom, and for the LAST restart output time step:
!!
!!    Symbolic link to data file
!!    --------------------------
!!      "restart_<modeltype>_DOMxx.nc"
!!
!!    Note:
!!    -  The domain-dependent suffix "...DOMxx" is also required for non-nested setups.
!!
!! --------------------------------------------------------------------------------
!!
!OPTION! -pvctl conflict
MODULE mo_sync_restart

  USE mo_communication,             ONLY: t_comm_gather_pattern, exchange_data
  USE mtime,                        ONLY: datetime
  USE mo_dynamics_config,           ONLY: nold, nnow, nnew, nnew_rcf, nnow_rcf
  USE mo_exception,                 ONLY: finish, message
  USE mo_fortran_tools,             ONLY: t_ptr_2d, t_ptr_2d_sp, t_ptr_2d_int
  USE mo_grid_config,               ONLY: n_dom
  USE mo_impl_constants,            ONLY: SUCCESS, SINGLE_T, REAL_T, INT_T
  USE mo_kind,                      ONLY: dp, sp
  USE mo_mpi,                       ONLY: my_process_is_mpi_workroot, my_process_is_mpi_test
  USE mo_restart_attributes,        ONLY: t_RestartAttributeList, RestartAttributeList_make
  USE mo_restart_descriptor,        ONLY: t_RestartDescriptor, t_RestartPatchData
  USE mo_restart_file,              ONLY: t_RestartFile
  USE mo_restart_patch_description, ONLY: t_restart_patch_description
  USE mo_restart_util,              ONLY: setDynamicPatchRestartAttributes, setGeneralRestartAttributes, &
    &                                     setPhysicsRestartAttributes, t_restart_args
  USE mo_restart_var_data,          ONLY: getLevelPointers, has_valid_time_level
  USE mo_timer,                     ONLY: timer_start, timer_stop, timer_write_restart, timer_write_restart_io, &
                                        & timer_write_restart_communication, timers_level
  USE mo_util_string,               ONLY: int2string
  USE mo_var_metadata_types,        ONLY: t_var_metadata

  IMPLICIT NONE

  PRIVATE

  INCLUDE 'netcdf.inc'

  PUBLIC :: t_SyncRestartDescriptor

  TYPE, EXTENDS(t_RestartPatchData) :: t_SyncPatchData
  CONTAINS
    PROCEDURE :: writeData => syncPatchData_writeData  ! override
  END TYPE t_SyncPatchData

  TYPE, EXTENDS(t_RestartDescriptor) :: t_SyncRestartDescriptor
  CONTAINS
    PROCEDURE :: construct => syncRestartDescriptor_construct   ! override
    PROCEDURE :: writeRestart => syncRestartDescriptor_writeRestart ! override
    PROCEDURE :: destruct => syncRestartDescriptor_destruct ! override
  END TYPE t_SyncRestartDescriptor

  !> module name string
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_sync_restart'

  !------------------------------------------------------------------------------------------------
CONTAINS

  SUBROUTINE syncRestartDescriptor_construct(me, modelType)
    CLASS(t_SyncRestartDescriptor), INTENT(INOUT) :: me
    INTEGER :: jg, error
    CHARACTER(LEN = *), INTENT(IN) :: modelType
    CHARACTER(LEN=*), PARAMETER :: routine = modname//':restartDescriptor_construct'

    IF(timers_level >= 5) CALL timer_start(timer_write_restart)

    ! ensure that our subtimers are consistently seen as subtimers
    IF(timers_level >= 7) THEN
        CALL timer_start(timer_write_restart_io)
        CALL timer_stop(timer_write_restart_io)
        CALL timer_start(timer_write_restart_communication)
        CALL timer_stop(timer_write_restart_communication)
    END IF

    CALL me%restartDescriptor_construct(modelType)

    ! allocate patch data structure
    ALLOCATE(t_SyncPatchData :: me%patchData(n_dom), STAT = error)
    IF(error /= SUCCESS) CALL finish(routine, "memory allocation failed")

    ! initialize the patch data structures
    DO jg = 1, n_dom
        CALL me%patchData(jg)%construct(modelType, jg)
    END DO
    IF(timers_level >= 5) CALL timer_stop(timer_write_restart)
  END SUBROUTINE syncRestartDescriptor_construct

  SUBROUTINE syncRestartDescriptor_defineRestartAttributes(me, restartAttributes, this_datetime, jstep, opt_output_jfile)
    CLASS(t_SyncRestartDescriptor), TARGET, INTENT(IN) :: me
    TYPE(t_RestartAttributeList), INTENT(INOUT) :: restartAttributes
    TYPE(datetime), POINTER, INTENT(IN) :: this_datetime
    INTEGER, VALUE :: jstep
    INTEGER, INTENT(IN), OPTIONAL :: opt_output_jfile(:)

    INTEGER :: jg, effectiveDomainCount
    CHARACTER(LEN = 2) :: jgString
    TYPE(t_restart_patch_description), POINTER :: curDescription
    CHARACTER(LEN = *), PARAMETER :: routine = modname//":syncRestartDescriptor_defineRestartAttributes"

    ! first the attributes that are independent of the domain
    effectiveDomainCount = 1
    IF(ALLOCATED(me%patchData(1)%description%opt_ndom)) effectiveDomainCount = me%patchData(1)%description%opt_ndom
    CALL setGeneralRestartAttributes(restartAttributes, this_datetime, effectiveDomainCount, jstep, opt_output_jfile)

    ! now the stuff that depends on the domain
    DO jg = 1, n_dom
        curDescription => me%patchData(jg)%description
        jgString = TRIM(int2string(jg, "(i2.2)"))
        CALL setDynamicPatchRestartAttributes(restartAttributes, jg, nold(jg), nnow(jg), nnew(jg), nnow_rcf(jg), nnew_rcf(jg))

        !----------------
        ! additional restart-output for nonhydrostatic model, not available anymore LK
        !IF(ALLOCATED(curDescription%opt_sim_time)) CALL restartAttributes%setReal('sim_time_DOM'//jgString, &
        !                                                                         &curDescription%opt_sim_time )

        !-------------------------------------------------------------
        ! DR
        ! WORKAROUND FOR FIELDS WHICH NEED TO GO INTO THE RESTART FILE,
        ! BUT SO FAR CANNOT BE HANDELED CORRECTLY BY ADD_VAR OR
        ! SET_RESTART_ATTRIBUTE
        !-------------------------------------------------------------

        IF(ALLOCATED(curDescription%opt_ndyn_substeps)) THEN
            CALL restartAttributes%setInteger('ndyn_substeps_DOM'//jgString, curDescription%opt_ndyn_substeps)
        END IF
        IF(ALLOCATED(curDescription%opt_jstep_adv_marchuk_order)) THEN
            CALL restartAttributes%setInteger('jstep_adv_marchuk_order_DOM'//jgString, curDescription%opt_jstep_adv_marchuk_order)
        END IF

        CALL setPhysicsRestartAttributes(restartAttributes, jg, curDescription%opt_t_elapsed_phy(:))
    END DO
  END SUBROUTINE syncRestartDescriptor_defineRestartAttributes

  SUBROUTINE syncRestartDescriptor_writeRestart(me, this_datetime, jstep, opt_output_jfile)
    CLASS(t_SyncRestartDescriptor), INTENT(INOUT) :: me
    TYPE(datetime), POINTER, INTENT(IN) :: this_datetime
    INTEGER, INTENT(IN) :: jstep
    INTEGER, INTENT(IN), OPTIONAL :: opt_output_jfile(:)

    INTEGER :: jg
    TYPE(t_RestartAttributeList), POINTER :: restartAttributes
    TYPE(t_restart_args) :: restartArgs
    CHARACTER(*), PARAMETER :: routine = modname//":syncRestartDescriptor_writeRestart"

    IF(timers_level >= 5) CALL timer_start(timer_write_restart)

    DO jg = 1, n_dom
        CALL me%patchData(jg)%description%setTimeLevels() !update the time levels
    END DO

    restartAttributes => RestartAttributeList_make()
    CALL restartArgs%construct(this_datetime, jstep, me%modelType, opt_output_jfile)
    CALL me%defineRestartAttributes(restartAttributes, restartArgs)

    DO jg = 1, n_dom
        CALL me%patchData(jg)%writeFile(restartAttributes, restartArgs, my_process_is_mpi_workroot())
    END DO

    CALL restartArgs%destruct()
    CALL restartAttributes%destruct()
    DEALLOCATE(restartAttributes)

    IF(timers_level >= 5) CALL timer_stop(timer_write_restart)
  END SUBROUTINE syncRestartDescriptor_writeRestart

  SUBROUTINE syncRestartDescriptor_destruct(me)
    CLASS(t_SyncRestartDescriptor), INTENT(INOUT) :: me

    INTEGER :: i

    IF(timers_level >= 5) CALL timer_start(timer_write_restart)

    IF(ASSOCIATED(me%patchData)) THEN
        DO i = 1, SIZE(me%patchData, 1)
            CALL me%patchData(i)%destruct()
        END DO
        DEALLOCATE(me%patchData)
    END IF
    IF(timers_level >= 5) CALL timer_stop(timer_write_restart)
  END SUBROUTINE syncRestartDescriptor_destruct

  ! loop over all var_lists for restart
  SUBROUTINE syncPatchData_writeData(me, file)
    CLASS(t_SyncPatchData), INTENT(INOUT) :: me
    TYPE(t_RestartFile), INTENT(INOUT) :: file

    INTEGER                              :: domain, i, gridSize, error, level
    TYPE(t_var_metadata), POINTER        :: info
    TYPE(t_comm_gather_pattern), POINTER :: gatherPattern
    REAL(dp), ALLOCATABLE                :: gatherBuffer_dp(:)
    REAL(sp), ALLOCATABLE                :: gatherBuffer_sp(:)
    INTEGER, ALLOCATABLE                 :: gatherBuffer_int(:)
    TYPE(t_ptr_2d), ALLOCATABLE          :: levelPointers_dp(:)
    TYPE(t_ptr_2d_sp), ALLOCATABLE       :: levelPointers_sp(:)
    TYPE(t_ptr_2d_int), ALLOCATABLE      :: levelPointers_int(:)
    CHARACTER(*), PARAMETER              :: routine = modname//":syncPatchData_writeData"

    IF(my_process_is_mpi_test()) RETURN

    domain = me%description%id
    DO i = 1, SIZE(me%varData)
        info => me%varData(i)%info
        IF(.NOT.has_valid_time_level(info, domain, nnew(domain), nnew_rcf(domain))) CYCLE

        ! we are committed to writing now
        IF(my_process_is_mpi_workroot()) write (0,*)' ... write '//TRIM(info%name)

        ! ALLOCATE the global array to gather the DATA on the master process
        gridSize = me%description%getGlobalGridSize(info%hgrid)
        gatherPattern => me%description%getGatherPattern(info%hgrid)

        ! get pointers to the local DATA
        SELECT CASE(info%data_type)
        CASE(REAL_T)
          !
          ALLOCATE(gatherBuffer_dp(MERGE(gridSize, 0, my_process_is_mpi_workroot())), STAT = error)
          IF(error /= SUCCESS) CALL finish(routine, "memory allocation failed")

          CALL getLevelPointers(info, me%varData(i)%r_ptr, levelPointers_dp)

          ! gather the data in the master process and write it to disk
          DO level = 1, SIZE(levelPointers_dp)
            IF(timers_level >= 7) CALL timer_start(timer_write_restart_communication)
            CALL exchange_data(in_array = levelPointers_dp(level)%p, out_array = gatherBuffer_dp, &
              &                gather_pattern = gatherPattern)
            IF(timers_level >= 7) CALL timer_stop(timer_write_restart_communication)
            IF(my_process_is_mpi_workroot()) THEN
              CALL file%writeLevel(info%cdiVarID, level - 1, gatherBuffer_dp)
            END IF
          END DO
          ! deallocate temporary global arrays
          DEALLOCATE(gatherBuffer_dp)
          !
        CASE(SINGLE_T)
          !
          ALLOCATE(gatherBuffer_sp(MERGE(gridSize, 0, my_process_is_mpi_workroot())), STAT = error)
          IF(error /= SUCCESS) CALL finish(routine, "memory allocation failed")

          CALL getLevelPointers(info, me%varData(i)%s_ptr, levelPointers_sp)

          ! gather the data in the master process and write it to disk
          DO level = 1, SIZE(levelPointers_sp)
            IF(timers_level >= 7) CALL timer_start(timer_write_restart_communication)
            CALL exchange_data(in_array = levelPointers_sp(level)%p, out_array = gatherBuffer_sp, &
              &                gather_pattern = gatherPattern)
            IF(timers_level >= 7) CALL timer_stop(timer_write_restart_communication)
            IF(my_process_is_mpi_workroot()) THEN
              CALL file%writeLevel(info%cdiVarID, level - 1, gatherBuffer_sp)
            END IF
          END DO
          ! deallocate temporary global arrays
          DEALLOCATE(gatherBuffer_sp)
          !
        CASE(INT_T)
          !
          ALLOCATE(gatherBuffer_dp(MERGE(gridSize, 0, my_process_is_mpi_workroot())), STAT = error)
          IF(error /= SUCCESS) CALL finish(routine, "memory allocation failed")
          ALLOCATE(gatherBuffer_int(SIZE(gatherBuffer_dp)), STAT = error)
          IF(error /= SUCCESS) CALL finish(routine, "memory allocation failed")

          CALL getLevelPointers(info, me%varData(i)%i_ptr, levelPointers_int)

          ! gather the data in the master process and write it to disk
          DO level = 1, SIZE(levelPointers_int)
            CALL exchange_data(in_array = levelPointers_int(level)%p, out_array = gatherBuffer_int, &
              &                gather_pattern = gatherPattern)
            IF(my_process_is_mpi_workroot()) THEN
              gatherBuffer_dp(:) = REAL(gatherBuffer_int(:), dp)
              CALL file%writeLevel(info%cdiVarID, level - 1, gatherBuffer_dp)
            END IF
          END DO
          ! deallocate temporary global arrays
          DEALLOCATE(gatherBuffer_int, gatherBuffer_dp)
          !
        CASE DEFAULT
          CALL finish(routine, "Internal error! Variable "//TRIM(info%name))
        END SELECT

        ! no deallocation of levelPointers so that the next invocation
        ! of getLevelPointers() may reuse the allocation
    END DO

    CALL message('','Finished writing restart')
  END SUBROUTINE syncPatchData_writeData

END MODULE mo_sync_restart
