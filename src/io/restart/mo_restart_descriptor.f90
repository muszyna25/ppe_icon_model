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
    USE mo_fortran_tools, ONLY: t_destructible
    USE mo_kind, ONLY: wp
    USE mo_model_domain, ONLY: t_patch

    IMPLICIT NONE

    PUBLIC :: t_RestartDescriptor

    PRIVATE

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
    TYPE, ABSTRACT, EXTENDS(t_destructible) :: t_RestartDescriptor
    CONTAINS
        PROCEDURE(restartDescriptor_construct), DEFERRED :: construct
        PROCEDURE(restartDescriptor_updatePatch), DEFERRED :: updatePatch
        PROCEDURE(restartDescriptor_writeRestart), DEFERRED :: writeRestart
    END TYPE t_RestartDescriptor

    ABSTRACT INTERFACE
        ! Constructor. Not called directly from user code, USE the factory FUNCTION createRestartDescriptor() instead, which IS found IN mo_restart.
        SUBROUTINE restartDescriptor_construct(me)
            IMPORT t_RestartDescriptor
            CLASS(t_RestartDescriptor), INTENT(INOUT) :: me
        END SUBROUTINE restartDescriptor_construct

        ! Update the internal description of the given patch. This should be called once for every patch before every CALL to writeRestart().
        SUBROUTINE restartDescriptor_updatePatch(me, patch, opt_pvct, opt_t_elapsed_phy, opt_lcall_phy, opt_sim_time, &
                                                &opt_ndyn_substeps, opt_jstep_adv_marchuk_order, opt_depth, opt_depth_lnd, &
                                                &opt_nlev_snow, opt_nice_class, opt_ndom, opt_ocean_zlevels, &
                                                &opt_ocean_zheight_cellMiddle, opt_ocean_zheight_cellInterfaces)
            IMPORT t_RestartDescriptor, t_patch, wp
            CLASS(t_RestartDescriptor), INTENT(INOUT) :: me
            TYPE(t_patch), INTENT(IN) :: patch
            INTEGER, INTENT(IN), OPTIONAL :: opt_depth, opt_depth_lnd, opt_ndyn_substeps, opt_jstep_adv_marchuk_order, &
                                           & opt_nlev_snow, opt_nice_class, opt_ndom, opt_ocean_zlevels
            REAL(wp), INTENT(IN), OPTIONAL :: opt_sim_time, opt_pvct(:), opt_t_elapsed_phy(:), opt_ocean_zheight_cellMiddle(:), &
                                           & opt_ocean_zheight_cellInterfaces(:)
            LOGICAL, INTENT(IN), OPTIONAL :: opt_lcall_phy(:)
        END SUBROUTINE restartDescriptor_updatePatch

        ! Actually WRITE a restart.
        SUBROUTINE restartDescriptor_writeRestart(me, datetime, jstep, modelType, opt_output_jfile)
            IMPORT t_RestartDescriptor, t_datetime
            CLASS(t_RestartDescriptor), INTENT(INOUT) :: me
            TYPE(t_datetime), INTENT(IN) :: datetime
            INTEGER, INTENT(IN) :: jstep
            CHARACTER(LEN = *), INTENT(IN) :: modelType
            INTEGER, INTENT(IN), OPTIONAL :: opt_output_jfile(:)
        END SUBROUTINE restartDescriptor_writeRestart
    END INTERFACE

END MODULE mo_restart_descriptor
