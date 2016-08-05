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
  USE mo_datetime, ONLY: t_datetime
  USE mo_dynamics_config, ONLY: nold, nnow, nnew, nnew_rcf, nnow_rcf
  USE mo_exception, ONLY: finish, message
  USE mo_fortran_tools, ONLY: t_ptr_2d
  USE mo_grid_config, ONLY: n_dom
  USE mo_impl_constants, ONLY: SUCCESS
  USE mo_kind, ONLY: wp
  USE mo_model_domain, ONLY: t_patch, p_patch
  USE mo_mpi, ONLY: my_process_is_mpi_workroot, my_process_is_mpi_test
  USE mo_restart_attributes, ONLY: t_RestartAttributeList, RestartAttributeList_make
  USE mo_restart_descriptor, ONLY: t_RestartDescriptor
  USE mo_restart_file, ONLY: t_RestartFile
  USE mo_restart_patch_description, ONLY: t_restart_patch_description
  USE mo_restart_util, ONLY: setGeneralRestartAttributes, setDynamicPatchRestartAttributes, setPhysicsRestartAttributes, &
                           & create_restart_file_link, t_restart_args
  USE mo_restart_var_data, ONLY: getLevelPointers, has_valid_time_level
  USE mo_restart_var_data, ONLY: t_RestartVarData, createRestartVarData
  USE mo_run_config, ONLY: ltimer
  USE mo_timer, ONLY: timer_start, timer_stop, timer_write_restart_file
  USE mo_util_string, ONLY: int2string
  USE mo_var_metadata_types, ONLY: t_var_metadata

  IMPLICIT NONE

  PRIVATE

  INCLUDE 'netcdf.inc'

  PUBLIC :: t_SyncRestartDescriptor

  ! this type stores all the information that we need to know about a patch and its variables
  TYPE t_PatchData
    TYPE(t_restart_patch_description) :: description
    TYPE(t_RestartVarData), POINTER :: varData(:)
    INTEGER :: restartType
  CONTAINS
    PROCEDURE :: writeFile => patchData_writeFile

    PROCEDURE, PRIVATE :: writeData => patchData_writeData  ! implementation detail of writeFile()
  END TYPE t_PatchData

  TYPE, EXTENDS(t_RestartDescriptor) :: t_SyncRestartDescriptor
    TYPE(t_PatchData), ALLOCATABLE :: patchData(:)
    CHARACTER(:), ALLOCATABLE :: modelType
  CONTAINS
    PROCEDURE :: construct => syncRestartDescriptor_construct   ! override
    PROCEDURE :: updatePatch => syncRestartDescriptor_updatePatch   ! override
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
    ALLOCATE(me%patchData(n_dom), STAT = error)
    IF(error /= SUCCESS) CALL finish(routine, "memory allocation failed")

    ! initialize the patch data structures
    DO jg = 1, n_dom
        ! construct the subobjects
        CALL me%patchData(jg)%description%init(p_patch(jg))
        me%patchData(jg)%varData => createRestartVarData(jg, modelType, opt_out_restartType = me%patchData(jg)%restartType)
    END DO
  END SUBROUTINE syncRestartDescriptor_construct

  SUBROUTINE syncRestartDescriptor_updatePatch(me, patch, opt_pvct, opt_t_elapsed_phy, opt_lcall_phy, opt_sim_time, &
                                        &opt_ndyn_substeps, opt_jstep_adv_marchuk_order, opt_depth_lnd, &
                                        &opt_nlev_snow, opt_nice_class, opt_ndom, opt_ocean_zlevels, &
                                        &opt_ocean_zheight_cellMiddle, opt_ocean_zheight_cellInterfaces)
    CLASS(t_SyncRestartDescriptor), INTENT(INOUT) :: me
    TYPE(t_patch), INTENT(IN) :: patch
    INTEGER, INTENT(IN), OPTIONAL :: opt_depth_lnd, opt_ndyn_substeps, opt_jstep_adv_marchuk_order, &
                                   & opt_nlev_snow, opt_nice_class, opt_ndom, opt_ocean_zlevels
    REAL(wp), INTENT(IN), OPTIONAL :: opt_sim_time, opt_pvct(:), opt_t_elapsed_phy(:), opt_ocean_zheight_cellMiddle(:), &
                                    & opt_ocean_zheight_cellInterfaces(:)
    LOGICAL, INTENT(IN), OPTIONAL :: opt_lcall_phy(:)

    INTEGER :: jg
    CHARACTER(LEN = *), PARAMETER :: routine = modname//":restartDescriptor_updatePatch"

    jg = patch%id
    IF(jg < 1 .OR. jg > SIZE(me%patchData)) CALL finish(routine, "assertion failed: patch id is out of range")
    IF(me%patchData(jg)%description%id /= jg) CALL finish(routine, "assertion failed: patch id doesn't match its array index")
    CALL me%patchData(jg)%description%update(patch, opt_pvct, opt_t_elapsed_phy, opt_lcall_phy, opt_sim_time, &
                                             &opt_ndyn_substeps, opt_jstep_adv_marchuk_order, opt_depth_lnd, &
                                             &opt_nlev_snow, opt_nice_class, opt_ndom, opt_ocean_zlevels, &
                                             &opt_ocean_zheight_cellMiddle, opt_ocean_zheight_cellInterfaces)
  END SUBROUTINE syncRestartDescriptor_updatePatch

  SUBROUTINE syncRestartDescriptor_defineRestartAttributes(me, restartAttributes, datetime, jstep, opt_output_jfile)
    CLASS(t_SyncRestartDescriptor), TARGET, INTENT(IN) :: me
    TYPE(t_RestartAttributeList), INTENT(INOUT) :: restartAttributes
    TYPE(t_datetime), INTENT(IN) :: datetime
    INTEGER, VALUE :: jstep
    INTEGER, INTENT(IN), OPTIONAL :: opt_output_jfile(:)

    INTEGER :: jg, effectiveDomainCount
    CHARACTER(LEN = 2) :: jgString
    TYPE(t_restart_patch_description), POINTER :: curDescription
    CHARACTER(LEN = *), PARAMETER :: routine = modname//":syncRestartDescriptor_defineRestartAttributes"

    ! first the attributes that are independent of the domain
    effectiveDomainCount = 1
    IF(me%patchData(1)%description%l_opt_ndom) effectiveDomainCount = me%patchData(1)%description%opt_ndom
    CALL setGeneralRestartAttributes(restartAttributes, datetime, effectiveDomainCount, jstep, opt_output_jfile)

    ! now the stuff that depends on the domain
    DO jg = 1, n_dom
        curDescription => me%patchData(jg)%description
        jgString = TRIM(int2string(jg, "(i2.2)"))
        CALL setDynamicPatchRestartAttributes(restartAttributes, jg, nold(jg), nnow(jg), nnew(jg), nnow_rcf(jg), nnew_rcf(jg))

        !----------------
        ! additional restart-output for nonhydrostatic model
        IF(curDescription%l_opt_sim_time) CALL restartAttributes%setReal('sim_time_DOM'//jgString, curDescription%opt_sim_time )

        !-------------------------------------------------------------
        ! DR
        ! WORKAROUND FOR FIELDS WHICH NEED TO GO INTO THE RESTART FILE,
        ! BUT SO FAR CANNOT BE HANDELED CORRECTLY BY ADD_VAR OR
        ! SET_RESTART_ATTRIBUTE
        !-------------------------------------------------------------

        IF(curDescription%l_opt_ndyn_substeps) THEN
            CALL restartAttributes%setInteger('ndyn_substeps_DOM'//jgString, curDescription%opt_ndyn_substeps)
        END IF
        IF(curDescription%l_opt_jstep_adv_marchuk_order) THEN
            CALL restartAttributes%setInteger('jstep_adv_marchuk_order_DOM'//jgString, curDescription%opt_jstep_adv_marchuk_order)
        END IF

        IF(ALLOCATED(curDescription%opt_t_elapsed_phy) .AND. ALLOCATED(curDescription%opt_lcall_phy)) THEN
            CALL setPhysicsRestartAttributes(restartAttributes, jg, curDescription%opt_t_elapsed_phy(:), &
                                                                  & curDescription%opt_lcall_phy(:))
        END IF
    END DO
  END SUBROUTINE syncRestartDescriptor_defineRestartAttributes

  SUBROUTINE patchData_writeFile(me, restartAttributes, restartArgs)
    CLASS(t_PatchData), INTENT(INOUT) :: me
    TYPE(t_RestartAttributeList) :: restartAttributes
    TYPE(t_restart_args), INTENT(IN) :: restartArgs

    TYPE(t_RestartFile) :: file
    CHARACTER(LEN = *), PARAMETER :: routine = modname//":patchData_writeFile"

    IF(ALLOCATED(me%description%opt_ocean_zheight_cellMiddle)) THEN
      IF(.NOT. ALLOCATED(me%description%opt_ocean_Zheight_CellInterfaces) .OR. .NOT. me%description%l_opt_ocean_Zlevels) THEN
          CALL finish('patchData_writeFile','Ocean level parameteres not complete')
      END IF
    END IF

    CALL me%description%defineVGrids()
    CALL me%description%setTimeLevels() !update the time levels

    IF(ASSOCIATED(me%varData)) THEN ! no restart variables => no restart file
        IF(my_process_is_mpi_workroot()) CALL file%open(me%description, me%varData, restartArgs, restartAttributes)
        CALL me%writeData(file)
        IF(my_process_is_mpi_workroot()) THEN
            IF(me%description%l_opt_ndom) THEN
                CALL create_restart_file_link(TRIM(file%filename), TRIM(restartArgs%modelType), 0, me%description%id, &
                                             &opt_ndom = me%description%opt_ndom)
            ELSE
                CALL create_restart_file_link(TRIM(file%filename), TRIM(restartArgs%modelType), 0, me%description%id)
            END IF
            CALL file%close()
        END IF
    END IF
  END SUBROUTINE patchData_writeFile

  SUBROUTINE syncRestartDescriptor_writeRestart(me, datetime, jstep, opt_output_jfile)
    CLASS(t_SyncRestartDescriptor), INTENT(INOUT) :: me
    TYPE(t_datetime), INTENT(IN) :: datetime
    INTEGER, INTENT(IN) :: jstep
    INTEGER, INTENT(IN), OPTIONAL :: opt_output_jfile(:)

    INTEGER :: jg
    TYPE(t_RestartAttributeList), POINTER :: restartAttributes
    TYPE(t_restart_args) :: restartArgs

    IF(ltimer) CALL timer_start(timer_write_restart_file)

    restartAttributes => RestartAttributeList_make()
    CALL me%defineRestartAttributes(restartAttributes, datetime, jstep, opt_output_jfile)
    CALL restartArgs%construct(datetime, jstep, me%modelType, opt_output_jfile)

    DO jg = 1, n_dom
        CALL me%patchData(jg)%writeFile(restartAttributes, restartArgs)
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

  SUBROUTINE defineRestartAttributes(restartAttributes, datetime, jstep, opt_ndom, opt_ndyn_substeps, &
                                    &opt_jstep_adv_marchuk_order, opt_output_jfile, opt_sim_time, opt_t_elapsed_phy, opt_lcall_phy)
    TYPE(t_RestartAttributeList), INTENT(INOUT) :: restartAttributes
    TYPE(t_datetime), INTENT(IN) :: datetime
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
    CALL setGeneralRestartAttributes(restartAttributes, datetime, effectiveDomainCount, jstep, opt_output_jfile)

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
  SUBROUTINE patchData_writeData(me, file)
    CLASS(t_PatchData), INTENT(IN) :: me
    TYPE(t_RestartFile) :: file

    INTEGER :: domain, i, gridSize, error, level
    TYPE(t_var_metadata), POINTER :: info
    TYPE(t_comm_gather_pattern), POINTER :: gatherPattern
    REAL(wp), ALLOCATABLE :: gatherBuffer(:)
    TYPE(t_ptr_2d), ALLOCATABLE :: levelPointers(:)
    CHARACTER(*), PARAMETER :: routine = modname//":patchData_writeData"

    IF(my_process_is_mpi_test()) RETURN

    domain = me%description%id
    DO i = 1, SIZE(me%varData)
        info => me%varData(i)%info
        IF(.NOT.has_valid_time_level(info, domain, nnew(domain), nnew_rcf(domain))) CYCLE

        ! we are committed to writing now
        IF(my_process_is_mpi_workroot()) write (0,*)' ... write ',info%name

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
  END SUBROUTINE patchData_writeData

END MODULE mo_sync_restart
