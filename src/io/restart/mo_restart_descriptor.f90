!>
!! The base CLASS for the PUBLIC INTERFACE for restart writing.
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.

MODULE mo_restart_descriptor
    USE mo_datetime, ONLY: t_datetime
    USE mo_exception, ONLY: finish
    USE mo_fortran_tools, ONLY: t_Destructible
    USE mo_kind, ONLY: wp
    USE mo_model_domain, ONLY: t_patch
    USE mo_mpi, ONLY: my_process_is_work
    USE mo_restart_attributes, ONLY: t_RestartAttributeList
    USE mo_restart_file, ONLY: t_RestartFile
    USE mo_restart_patch_description, ONLY: t_restart_patch_description
    USE mo_restart_util, ONLY: t_restart_args, create_restart_file_link
    USE mo_restart_var_data, ONLY: t_RestartVarData, createRestartVarData

    IMPLICIT NONE

    PUBLIC :: t_RestartDescriptor
    PUBLIC :: t_RestartPatchData

    PRIVATE

    ! this type stores all the information that we need to know about a patch and its variables
    TYPE, EXTENDS(t_Destructible) :: t_RestartPatchData
        TYPE(t_restart_patch_description) :: description
        TYPE(t_RestartVarData), POINTER :: varData(:)
        INTEGER :: restartType
    CONTAINS
        PROCEDURE :: construct => restartPatchData_construct
        PROCEDURE :: writeData => restartPatchData_writeDataDummy    ! >>> Must be overriden by derived class. I would have made this deferred, had it not been for a bug in the NAG compiler that makes calling of an inherited function from an abstract base impossible.
        PROCEDURE :: writeFile => restartPatchData_writeFile
        PROCEDURE :: destruct => restartPatchData_destruct
    END TYPE t_RestartPatchData

    ! This IS the actual INTERFACE to the restart writing code (apart from the restart_main_proc PROCEDURE). Its USE IS as follows:
    !
    ! First, AND ONLY once during a run, a t_RestartDescriptor IS created.
    ! This IS done via the factory FUNCTION createRestartDescriptor(), which IS found IN mo_restart.
    !
    ! Then, for each restart that IS to be written, the updatePatch() method IS used to set the current time dependend information for each patch.
    ! Once all patches are updated, a single CALL to writeRestart() triggers the actual restart writing.
    ! The updatePatch() - writeRestart() sequence can be repeated ANY number of times.
    !
    ! Finally, destruct() must be called for cleanup. This IS especially important IN the CASE of asynchronous restart writing,
    ! because the destruct() CALL will signal the restart PEs to finish their work, AND wait for them to stop.
    TYPE, ABSTRACT, EXTENDS(t_Destructible) :: t_RestartDescriptor
        !XXX: Using ALLOCATABLE instead of POINTER here seems to trigger a bug IN the cray compiler that leads to later allocations to fail.
        CLASS(t_RestartPatchData), POINTER :: patchData(:)   ! must be ALLOCATED IN the subclass constructor
    CONTAINS
        PROCEDURE(restartDescriptor_construct), DEFERRED :: construct
        PROCEDURE :: updatePatch => restartDescriptor_updatePatch
        PROCEDURE(restartDescriptor_writeRestart), DEFERRED :: writeRestart
    END TYPE t_RestartDescriptor

    ABSTRACT INTERFACE

        ! Constructor. Not called directly from user code, USE the factory FUNCTION createRestartDescriptor() instead, which IS found IN mo_restart.
        SUBROUTINE restartDescriptor_construct(me, modelType)
            IMPORT t_RestartDescriptor
            CLASS(t_RestartDescriptor), INTENT(INOUT) :: me
            CHARACTER(LEN = *), INTENT(IN) :: modelType
        END SUBROUTINE restartDescriptor_construct

        ! Actually WRITE a restart.
        SUBROUTINE restartDescriptor_writeRestart(me, datetime, jstep, opt_output_jfile)
            IMPORT t_RestartDescriptor, t_datetime
            CLASS(t_RestartDescriptor), INTENT(INOUT) :: me
            TYPE(t_datetime), INTENT(IN) :: datetime
            INTEGER, INTENT(IN) :: jstep
            INTEGER, INTENT(IN), OPTIONAL :: opt_output_jfile(:)
        END SUBROUTINE restartDescriptor_writeRestart

    END INTERFACE

    CHARACTER(*), PARAMETER :: modname = "mo_restart_descriptor"

CONTAINS

    ! Update the internal description of the given patch. This should be called once for every patch before every CALL to writeRestart().
    SUBROUTINE restartDescriptor_updatePatch(me, patch, opt_pvct, opt_t_elapsed_phy, opt_lcall_phy, opt_sim_time, &
                                            &opt_ndyn_substeps, opt_jstep_adv_marchuk_order, opt_depth_lnd, &
                                            &opt_nlev_snow, opt_nice_class, opt_ndom, opt_ocean_zlevels, &
                                            &opt_ocean_zheight_cellMiddle, opt_ocean_zheight_cellInterfaces)
        CLASS(t_RestartDescriptor), INTENT(INOUT) :: me
        TYPE(t_patch), INTENT(IN) :: patch
        INTEGER, INTENT(IN), OPTIONAL :: opt_depth_lnd, opt_ndyn_substeps, opt_jstep_adv_marchuk_order, &
                                       & opt_nlev_snow, opt_nice_class, opt_ndom, opt_ocean_zlevels
        REAL(wp), INTENT(IN), OPTIONAL :: opt_sim_time, opt_pvct(:), opt_t_elapsed_phy(:), opt_ocean_zheight_cellMiddle(:), &
                                       & opt_ocean_zheight_cellInterfaces(:)
        LOGICAL, INTENT(IN), OPTIONAL :: opt_lcall_phy(:)

        INTEGER :: jg
        CHARACTER(LEN = *), PARAMETER :: routine = modname//":restartDescriptor_updatePatch"

        IF(.NOT.my_process_is_work()) CALL finish(routine, "assertion failed")
        jg = patch%id
        IF(jg < 1 .OR. jg > SIZE(me%patchData)) CALL finish(routine, "assertion failed: patch id IS OUT of range")
        IF(me%patchData(jg)%description%id /= jg) CALL finish(routine, "assertion failed: patch id doesn't match its array index")
        CALL me%patchData(jg)%description%update(patch, opt_pvct, opt_t_elapsed_phy, opt_lcall_phy, opt_sim_time, &
                                                 &opt_ndyn_substeps, opt_jstep_adv_marchuk_order, opt_depth_lnd, &
                                                 &opt_nlev_snow, opt_nice_class, opt_ndom, opt_ocean_zlevels, &
                                                 &opt_ocean_zheight_cellMiddle, opt_ocean_zheight_cellInterfaces)
    END SUBROUTINE restartDescriptor_updatePatch

    SUBROUTINE restartPatchData_construct(me, modelType, domain)
        CLASS(t_RestartPatchData), INTENT(INOUT) :: me
        CHARACTER(*), INTENT(IN) :: modelType
        INTEGER, INTENT(IN) :: domain

        CALL me%description%init(domain)
        me%varData => createRestartVarData(domain, modelType, me%restartType)
    END SUBROUTINE restartPatchData_construct

    ! Write the payload DATA to an opened file.
    ! Depending on the implementation, this CALL may OR may NOT be collective (it IS collective for sync, AND non-collective for async restart writing).
    !
    ! Pure dummy implementation to avoid making t_RestartPatchData abstract.
    SUBROUTINE restartPatchData_writeDataDummy(me, file)
        CLASS(t_RestartPatchData), INTENT(INOUT) :: me
        TYPE(t_RestartFile), INTENT(INOUT) :: file

        CHARACTER(*), PARAMETER :: routine = modname//":restartPatchData_writeDataDummy"

        CALL finish(routine, "this dummy implementation must not be called")
    END SUBROUTINE restartPatchData_writeDataDummy

    ! Write a restart file if varData IS set. lIsWriteProcess should ONLY be set on one process.
    SUBROUTINE restartPatchData_writeFile(me, restartAttributes, restartArgs, procId, lIsWriteProcess)
        CLASS(t_RestartPatchData), INTENT(INOUT) :: me
        TYPE(t_RestartAttributeList), INTENT(INOUT) :: restartAttributes
        TYPE(t_restart_args), INTENT(IN) :: restartArgs
        INTEGER, VALUE :: procId
        LOGICAL, VALUE :: lIsWriteProcess

        TYPE(t_RestartFile) :: file

        IF(ASSOCIATED(me%varData)) THEN ! no restart variables => no restart file
            CALL me%description%defineVGrids()

            IF(lIsWriteProcess) CALL file%open(me%description, me%varData, restartArgs, restartAttributes, me%restartType)

            CALL me%writeData(file)

            IF(lIsWriteProcess) THEN
                IF(ALLOCATED(me%description%opt_ndom)) THEN
                    CALL create_restart_file_link(TRIM(file%filename), TRIM(restartArgs%modelType), procId, me%description%id, &
                                                 &opt_ndom = me%description%opt_ndom)
                ELSE
                    CALL create_restart_file_link(TRIM(file%filename), TRIM(restartArgs%modelType), procId, me%description%id)
                END IF
                CALL file%close()
            END IF
        END IF
    END SUBROUTINE restartPatchData_writeFile

    SUBROUTINE restartPatchData_destruct(me)
        CLASS(t_RestartPatchData), INTENT(INOUT) :: me

        IF(ASSOCIATED(me%varData)) DEALLOCATE(me%varData)
    END SUBROUTINE restartPatchData_destruct

END MODULE mo_restart_descriptor
