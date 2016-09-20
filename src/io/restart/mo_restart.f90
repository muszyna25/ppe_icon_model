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
    USE mo_async_restart, ONLY: t_AsyncRestartDescriptor
#endif
    USE mo_exception, ONLY: finish
    USE mo_impl_constants, ONLY: SUCCESS
    USE mo_parallel_config, ONLY: use_async_restart_output
    USE mo_restart_descriptor, ONLY: t_RestartDescriptor
    USE mo_sync_restart, ONLY: t_SyncRestartDescriptor

    IMPLICIT NONE

    ! documentation for t_RestartDescriptor IS found IN mo_restart_descriptor
    PUBLIC :: t_RestartDescriptor, createRestartDescriptor, deleteRestartDescriptor

    PRIVATE

    CHARACTER(LEN = *), PARAMETER :: modname = "mo_restart"

CONTAINS

    ! Factory FUNCTION to create the appropriate restart descriptor.
    FUNCTION createRestartDescriptor(modelType) RESULT(RESULT)
        CHARACTER(*), INTENT(IN) :: modelType
        CLASS(t_RestartDescriptor), POINTER :: RESULT

        INTEGER :: error
        CHARACTER(LEN = *), PARAMETER :: routine = modname//":createRestartDescriptor"

        IF(use_async_restart_output) THEN
#ifndef NOMPI
            ALLOCATE(t_AsyncRestartDescriptor :: RESULT, STAT = error)
#else
            CALL finish(routine, "this executable was compiled without MPI support, hence async restart writing is not available")
#endif
        ELSE
            ALLOCATE(t_SyncRestartDescriptor :: RESULT, STAT = error)
        END IF
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation failure")
        CALL RESULT%construct(modelType)
    END FUNCTION createRestartDescriptor

    ! Convenience FUNCTION for destroying a restart descriptor.
    SUBROUTINE deleteRestartDescriptor(descriptor)
        CLASS(t_RestartDescriptor), POINTER, INTENT(INOUT) :: descriptor

        CALL descriptor%destruct()
        DEALLOCATE(descriptor)
    END SUBROUTINE deleteRestartDescriptor

END MODULE mo_restart
