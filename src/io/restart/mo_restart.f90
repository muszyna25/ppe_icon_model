!>
!! This is the high-level module that is used to write a restart.
!! Basically it reexports t_RestartDescriptor and provides a factory function to create restart descriptors.
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.

MODULE mo_restart
#ifndef NOMPI
    USE mo_async_restart, ONLY: t_AsyncRestartDescriptor, asyncRestart_mainLoop
#endif
    USE mo_exception, ONLY: finish, message
    USE mo_impl_constants, ONLY: SUCCESS
    USE mo_io_config, ONLY: restartWritingParameters, kSyncRestartModule, kAsyncRestartModule, kMultifileRestartModule
    USE mo_mpi, ONLY: my_process_is_restart, process_mpi_restart_size
    USE mo_multifile_restart, ONLY: t_MultifileRestartDescriptor
#ifndef NOMPI
    USE mo_multifile_restart, ONLY: multifileRestart_mainLoop
#endif
    USE mo_restart_descriptor, ONLY: t_RestartDescriptor
    USE mo_restart_util, ONLY: becomeDedicatedRestartProc, shutdownRestartProc
    USE mo_sync_restart, ONLY: t_SyncRestartDescriptor
#ifdef YAC_coupling
  USE mo_coupling_config, ONLY: is_coupled_run
  USE mo_io_coupling, ONLY: construct_io_coupler, destruct_io_coupler
#endif

    IMPLICIT NONE

    ! documentation for t_RestartDescriptor IS found IN mo_restart_descriptor
    PUBLIC :: t_RestartDescriptor, createRestartDescriptor, deleteRestartDescriptor, detachRestartProcs

    PRIVATE

    CHARACTER(LEN = *), PARAMETER :: modname = "mo_restart"

CONTAINS

    ! Factory FUNCTION to create the appropriate restart descriptor.
    FUNCTION createRestartDescriptor(modelType) RESULT(resultVar)
        CHARACTER(*), INTENT(IN) :: modelType
        CLASS(t_RestartDescriptor), POINTER :: resultVar

        INTEGER :: error, restartModule
        LOGICAL :: lDedicatedProcMode
        CHARACTER(LEN = *), PARAMETER :: routine = modname//":createRestartDescriptor"

        CALL restartWritingParameters(opt_restartModule      = restartModule, &
          &                           opt_lDedicatedProcMode = lDedicatedProcMode)
        SELECT CASE(restartModule)
            CASE(kSyncRestartModule)
                CALL message('','synchronous restart writing selected.')
                ALLOCATE(t_SyncRestartDescriptor :: resultVar, STAT = error)
            CASE(kAsyncRestartModule)
                CALL message('','asynchronous restart writing selected.')
#ifdef NOMPI
                CALL finish(routine, "this executable was compiled without MPI support, hence async restart writing is not &
                                     &available")
#else
                ALLOCATE(t_AsyncRestartDescriptor :: resultVar, STAT = error)
#endif
            CASE(kMultifileRestartModule)
              IF (lDedicatedProcMode) THEN
                CALL message('','multifile restart writing selected, with dedicated procs.')
              ELSE
                CALL message('','multifile restart writing selected, joint proc mode.')
              END IF
              ALLOCATE(t_MultifileRestartDescriptor :: resultVar, STAT = error)
        END SELECT
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation failure")
        CALL resultVar%construct(modelType)
    END FUNCTION createRestartDescriptor

    ! Convenience FUNCTION for destroying a restart descriptor.
    SUBROUTINE deleteRestartDescriptor(descriptor)
        CLASS(t_RestartDescriptor), POINTER, INTENT(INOUT) :: descriptor

        CALL descriptor%destruct()
        DEALLOCATE(descriptor)
    END SUBROUTINE deleteRestartDescriptor

    ! Enter the restart main proc IF this IS a pure restart PE AND set use_async_restart_output accordingly.
    !
    ! This routine does NOT RETURN on dedicated restart processes.
    SUBROUTINE detachRestartProcs()
        INTEGER :: restartModule
        CHARACTER(*), PARAMETER :: routine = modname//":detachRestartProcs"

        IF(process_mpi_restart_size <= 0) RETURN    ! no dedicated restart processes configured -> noop
        IF(.NOT.my_process_is_restart()) RETURN ! this IS NOT a dedicated restart process -> noop

#ifdef NOMPI
        CALL finish('', 'this executable was compiled without MPI support, hence restart writing with dedicated restart processes &
                        &is not available')
#else
        ! Actually detach the restart processes.
#ifdef YAC_coupling
    ! The initialisation of YAC needs to be called by all (!) MPI processes
    ! in MPI_COMM_WORLD. Thus we do it here for the restart processes.

    IF ( is_coupled_run() ) CALL construct_io_coupler ( "restart_io" )
#endif

        CALL becomeDedicatedRestartProc()

        CALL restartWritingParameters(opt_restartModule = restartModule)
        SELECT CASE(restartModule)
            CASE(kSyncRestartModule)
                CALL finish(routine, "assertion failed: value of process_mpi_restart_size is inconsistent with result of &
                                     &restartWritingParameters()")

            CASE(kAsyncRestartModule)
                CALL asyncRestart_mainLoop()

            CASE(kMultifileRestartModule)
                CALL multifileRestart_mainLoop()
        END SELECT
#ifdef YAC_coupling
        IF ( is_coupled_run() ) CALL destruct_io_coupler ( "restart_io" )
#endif
        ! This IS a FUNCTION of no RETURN!
        CALL shutdownRestartProc()
#endif
    END SUBROUTINE detachRestartProcs

END MODULE mo_restart
