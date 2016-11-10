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
  USE mo_communication, ONLY: t_comm_gather_pattern, exchange_data
  USE mtime, ONLY: datetime
  USE mo_dynamics_config, ONLY: nold, nnow, nnew, nnew_rcf, nnow_rcf
  USE mo_exception, ONLY: finish, message
  USE mo_fortran_tools, ONLY: t_ptr_2d
  USE mo_grid_config, ONLY: n_dom
  USE mo_impl_constants, ONLY: SUCCESS
  USE mo_kind, ONLY: wp
  USE mo_model_domain, ONLY: t_patch, p_patch
  USE mo_mpi, ONLY: my_process_is_mpi_workroot, my_process_is_mpi_test, my_process_is_work
  USE mo_restart_attributes, ONLY: t_RestartAttributeList, RestartAttributeList_make
  USE mo_restart_descriptor, ONLY: t_RestartDescriptor, t_RestartPatchData
  USE mo_restart_file, ONLY: t_RestartFile
  USE mo_restart_patch_description, ONLY: t_restart_patch_description
  USE mo_restart_util, ONLY: setGeneralRestartAttributes, setDynamicPatchRestartAttributes, setPhysicsRestartAttributes, &
                           & create_restart_file_link, t_restart_args
  USE mo_restart_var_data, ONLY: getLevelPointers, has_valid_time_level, createRestartVarData
  USE mo_run_config, ONLY: ltimer
  USE mo_timer, ONLY: timer_start, timer_stop, timer_write_restart_file
  USE mo_util_string, ONLY: int2string
  USE mo_var_metadata_types, ONLY: t_var_metadata

  IMPLICIT NONE

  PRIVATE

  INCLUDE 'netcdf.inc'

  PUBLIC :: t_SyncRestartDescriptor

  TYPE, EXTENDS(t_RestartPatchData) :: t_SyncPatchData
  CONTAINS
    PROCEDURE :: writeData => syncPatchData_writeData  ! override
  END TYPE t_SyncPatchData

  TYPE, EXTENDS(t_RestartDescriptor) :: t_SyncRestartDescriptor
    CHARACTER(:), ALLOCATABLE :: modelType
  CONTAINS
    PROCEDURE :: construct => syncRestartDescriptor_construct   ! override
    PROCEDURE :: writeRestart => syncRestartDescriptor_writeRestart ! override
    PROCEDURE :: destruct => syncRestartDescriptor_destruct ! override

    PROCEDURE, PRIVATE :: defineRestartAttributes => syncRestartDescriptor_defineRestartAttributes
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

    me%modelType = modelType

    ! allocate patch data structure
    ALLOCATE(t_SyncPatchData :: me%patchData(n_dom), STAT = error)
    IF(error /= SUCCESS) CALL finish(routine, "memory allocation failed")

    ! initialize the patch data structures
    DO jg = 1, n_dom
        CALL me%patchData(jg)%construct(modelType, jg)
    END DO
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

        IF(ALLOCATED(curDescription%opt_t_elapsed_phy) .AND. ALLOCATED(curDescription%opt_lcall_phy)) THEN
            CALL setPhysicsRestartAttributes(restartAttributes, jg, curDescription%opt_t_elapsed_phy(:), &
                                                                  & curDescription%opt_lcall_phy(:))
        END IF
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

    IF(ltimer) CALL timer_start(timer_write_restart_file)

    restartAttributes => RestartAttributeList_make()
    CALL me%defineRestartAttributes(restartAttributes, this_datetime, jstep, opt_output_jfile)
    CALL restartArgs%construct(this_datetime, jstep, me%modelType, opt_output_jfile)

    DO jg = 1, n_dom
        CALL me%patchData(jg)%description%setTimeLevels() !update the time levels
        CALL me%patchData(jg)%writeFile(restartAttributes, restartArgs, 0, my_process_is_mpi_workroot())
    END DO

    CALL restartArgs%destruct()
    CALL restartAttributes%destruct()
    DEALLOCATE(restartAttributes)

    IF(ltimer) CALL timer_stop(timer_write_restart_file)
  END SUBROUTINE syncRestartDescriptor_writeRestart

  SUBROUTINE syncRestartDescriptor_destruct(me)
    CLASS(t_SyncRestartDescriptor), INTENT(INOUT) :: me

    DEALLOCATE(me%patchData)
  END SUBROUTINE syncRestartDescriptor_destruct

  SUBROUTINE defineRestartAttributes(restartAttributes, this_datetime, jstep, opt_ndom, opt_ndyn_substeps, &
                                    &opt_jstep_adv_marchuk_order, opt_output_jfile, opt_sim_time, opt_t_elapsed_phy, opt_lcall_phy)
    TYPE(t_RestartAttributeList), INTENT(INOUT) :: restartAttributes
    TYPE(datetime), POINTER, INTENT(IN) :: this_datetime
    INTEGER, VALUE :: jstep
    INTEGER, INTENT(IN), OPTIONAL :: opt_ndom, opt_ndyn_substeps, opt_jstep_adv_marchuk_order, opt_output_jfile(:)
    REAL(wp), INTENT(IN), OPTIONAL :: opt_sim_time, opt_t_elapsed_phy(:,:)
    LOGICAL , INTENT(IN), OPTIONAL :: opt_lcall_phy(:,:)

    INTEGER :: jg, effectiveDomainCount
    CHARACTER(LEN = 2) :: jgString
    CHARACTER(LEN = *), PARAMETER :: routine = modname//":defineRestartAttributes"

    ! first the attributes that are independent of the domain
    effectiveDomainCount = 1
    IF(PRESENT(opt_ndom)) effectiveDomainCount = opt_ndom
    CALL setGeneralRestartAttributes(restartAttributes, this_datetime, effectiveDomainCount, jstep, opt_output_jfile)

    ! now the stuff that depends on the domain
    DO jg = 1, n_dom
        jgString = TRIM(int2string(jg, "(i2.2)"))
        CALL setDynamicPatchRestartAttributes(restartAttributes, jg, nold(jg), nnow(jg), nnew(jg), nnow_rcf(jg), nnew_rcf(jg))

        !----------------
        ! additional restart-output for nonhydrostatic model
        IF (PRESENT(opt_sim_time)) CALL restartAttributes%setReal('sim_time_DOM'//jgString, opt_sim_time )

        !-------------------------------------------------------------
        ! DR
        ! WORKAROUND FOR FIELDS WHICH NEED TO GO INTO THE RESTART FILE,
        ! BUT SO FAR CANNOT BE HANDELED CORRECTLY BY ADD_VAR OR
        ! SET_RESTART_ATTRIBUTE
        !-------------------------------------------------------------

        IF (PRESENT(opt_ndyn_substeps)) CALL restartAttributes%setInteger('ndyn_substeps_DOM'//jgString, opt_ndyn_substeps)
        IF (PRESENT(opt_jstep_adv_marchuk_order)) CALL restartAttributes%setInteger('jstep_adv_marchuk_order_DOM'//jgString, &
                                                                                   &opt_jstep_adv_marchuk_order)

        IF (PRESENT(opt_t_elapsed_phy) .AND. PRESENT(opt_lcall_phy)) THEN
            CALL setPhysicsRestartAttributes(restartAttributes, jg, opt_t_elapsed_phy(jg,:), opt_lcall_phy(jg,:))
        ENDIF
    END DO
  END SUBROUTINE defineRestartAttributes

  ! loop over all var_lists for restart
  SUBROUTINE syncPatchData_writeData(me, file)
    CLASS(t_SyncPatchData), INTENT(INOUT) :: me
    TYPE(t_RestartFile), INTENT(INOUT) :: file

    INTEGER :: domain, i, gridSize, error, level
    TYPE(t_var_metadata), POINTER :: info
    TYPE(t_comm_gather_pattern), POINTER :: gatherPattern
    REAL(wp), ALLOCATABLE :: gatherBuffer(:)
    TYPE(t_ptr_2d), ALLOCATABLE :: levelPointers(:)
    CHARACTER(*), PARAMETER :: routine = modname//":syncPatchData_writeData"

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
        ALLOCATE(gatherBuffer(MERGE(gridSize, 0, my_process_is_mpi_workroot())), STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation failed")

        ! get pointers to the local DATA
        CALL getLevelPointers(info, me%varData(i)%r_ptr, levelPointers)

        ! gather the DATA IN the master process AND WRITE it to disk
        DO level = 1, SIZE(levelPointers)
            CALL exchange_data(in_array = levelPointers(level)%p, out_array = gatherBuffer, gather_pattern = gatherPattern)
            IF(my_process_is_mpi_workroot()) THEN
                CALL file%writeLevel(info%cdiVarID, level - 1, gatherBuffer)
            END IF
        END DO

        ! deallocate temporary global arrays
        DEALLOCATE(gatherBuffer)
        ! no deallocation of levelPointers so that the next invocation of getLevelPointers() may reuse the allocation
    END DO

    CALL message('','Finished writing restart')
  END SUBROUTINE syncPatchData_writeData

END MODULE mo_sync_restart
