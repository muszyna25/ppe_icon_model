!>
!! Contains routines for asynchronous restart Output
!! --------------------------------------------------------
!!
!! Note: The synchronous implementation of the restart output can be
!!       found in the module "mo_sync_restart". See module header for
!!       more details on generated files.
!!
!! @par Revision History
!! Initial implementation by Joerg Benkenstein (2013-01-15)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!!

MODULE mo_async_restart
! There is no point in pretending this module is usable if NOMPI is defined.
#ifndef NOMPI

  USE mo_async_restart_comm_data,   ONLY: t_AsyncRestartCommData
  USE mo_exception,                 ONLY: finish, message
  USE mo_fortran_tools,             ONLY: t_ptr_2d, t_ptr_2d_sp, t_ptr_2d_int
  USE mo_kind,                      ONLY: wp, i8, dp, sp
  USE mtime,                        ONLY: datetime
  USE mo_io_units,                  ONLY: nerr
  USE mo_restart_attributes,        ONLY: t_RestartAttributeList, RestartAttributeList_make
  USE mo_impl_constants,            ONLY: SUCCESS, SINGLE_T, REAL_T, INT_T
  USE mo_cdi_constants,             ONLY: GRID_UNSTRUCTURED_CELL, GRID_UNSTRUCTURED_VERT, GRID_UNSTRUCTURED_EDGE
  USE mo_restart_file,              ONLY: t_RestartFile
#ifdef DEBUG
  USE mo_restart_namelist,          ONLY: t_NamelistArchive, namelistArchive
#endif
  USE mo_parallel_config,           ONLY: restart_chunk_size
  USE mo_grid_config,               ONLY: n_dom
  USE mo_run_config,                ONLY: msg_level
  USE mo_packed_message,            ONLY: t_PackedMessage, kPackOp, kUnpackOp
  USE mo_restart_descriptor,        ONLY: t_RestartDescriptor, t_RestartPatchData
  USE mo_restart_util,              ONLY: t_restart_args
  USE mo_restart_var_data,          ONLY: t_RestartVarData, getLevelPointers, has_valid_time_level
  USE mo_timer,                     ONLY: timer_start, timer_stop, timer_write_restart, timer_write_restart_io, &
                                        & timer_write_restart_communication, timers_level
  USE mo_mpi,                       ONLY: p_pe, p_pe_work, p_restart_pe0, p_comm_work, p_work_pe0,       &
    &                                     num_work_procs, MPI_SUCCESS, p_send, p_recv,         &
    &                                     p_barrier, my_process_is_restart, my_process_is_work, &
    &                                     process_mpi_restart_size, p_mpi_wtime,  &
    &                                     process_mpi_all_comm
#ifdef __SUNPRO_F95
  INCLUDE "mpif.h"
#else
  USE mpi,                          ONLY: MPI_ADDRESS_KIND
#endif

  IMPLICIT NONE

  ! public routines
  PUBLIC :: t_AsyncRestartDescriptor
  PUBLIC :: asyncRestart_mainLoop

  PRIVATE

  ! byte sizes for double and single precision floats
  INTEGER, PARAMETER :: BYTES_DP                  = 8
  INTEGER, PARAMETER :: BYTES_SP                  = 4

  ! tags for communication between compute PEs and restart PEs
  INTEGER, PARAMETER :: MSG_RESTART_START         = 111111
  INTEGER, PARAMETER :: MSG_RESTART_DONE          = 222222
  INTEGER, PARAMETER :: MSG_RESTART_SHUTDOWN      = 999999

  ! common constant strings
  CHARACTER(LEN=*), PARAMETER :: modname                  = 'mo_async_restart'

  CHARACTER(LEN=*), PARAMETER :: NET_CDF_ERROR_FORMAT     = '(a,i5,a)'
  CHARACTER(LEN=*), PARAMETER :: FORMAT_VALS3             = '(a,a,i3)'
  CHARACTER(LEN=*), PARAMETER :: FORMAT_VALS5             = '(a,a,i3,a,i3)'
  CHARACTER(LEN=*), PARAMETER :: FORMAT_VALS5I            = '(a,a,i3,a,a)'
  CHARACTER(LEN=*), PARAMETER :: FORMAT_VALS7             = '(a,a,i3,a,i6,a,i3)'
  CHARACTER(LEN=*), PARAMETER :: FORMAT_VALS7I            = '(a,a,i3,a,a,a,i8)'

  ! combine the DATA that describes a patch for restart purposes with
  ! the infos required for the asynchronous fetching of the DATA from
  ! the compute PEs
  TYPE, EXTENDS(t_RestartPatchData) :: t_AsyncPatchData
    TYPE(t_AsyncRestartCommData) :: commData
  CONTAINS
    PROCEDURE          :: construct         => asyncPatchData_construct   ! override
    PROCEDURE          :: writeData         => asyncPatchData_writeData   ! override
    PROCEDURE          :: destruct          => asyncPatchData_destruct    ! override
  END TYPE t_AsyncPatchData

  ! This IS the actual INTERFACE to the restart writing code (apart
  ! from the asyncRestart_mainLoop PROCEDURE). Its USE IS as follows:
  !
  ! First, AND ONLY once during a run, a t_AsyncRestartDescriptor IS constructed.
  !
  ! Then, for each restart that IS to be written, the updatePatch()
  ! method IS used to set the current time dependend information for
  ! each patch.
  ! Once all patches are updated, a single CALL to writeRestart() triggers the actual restart writing.
  ! The updatePatch() - writeRestart() sequence can be repeated ANY number of times.
  !
  ! Finally, destruct() must be called to signal the restart PEs to
  ! finish their work, AND to wait for them to stop.
  TYPE, EXTENDS(t_RestartDescriptor) :: t_AsyncRestartDescriptor
  CONTAINS
    PROCEDURE :: construct => asyncRestartDescriptor_construct  ! override
    PROCEDURE :: writeRestart => asyncRestartDescriptor_writeRestart    ! override
    PROCEDURE :: destruct => asyncRestartDescriptor_destruct    ! override

    ! methods called ONLY by the restart processes
    PROCEDURE, PRIVATE :: restartWriteAsyncRestart => asyncRestartDescriptor_restartWriteAsyncRestart
  END TYPE t_AsyncRestartDescriptor

CONTAINS

  ! collective across restart AND worker PEs
  SUBROUTINE asyncPatchData_construct(me, modelType, domain)
    CLASS(t_AsyncPatchData), INTENT(INOUT) :: me
    CHARACTER(*), INTENT(IN) :: modelType
    INTEGER, INTENT(IN) :: domain

    CHARACTER(*), PARAMETER :: routine = modname//":asyncPatchData_construct"

    CALL me%t_RestartPatchData%construct(modelType, domain)
    CALL me%commData%construct(domain, me%varData)

    CALL me%transferToRestart()

    ! consistency checks
    IF(me%description%n_patch_cells_g /= me%commData%cells%n_glb) THEN
        CALL finish(routine, "assertion failed: mismatch of global cell count")
    END IF
    IF(me%description%n_patch_verts_g /= me%commData%verts%n_glb) THEN
        CALL finish(routine, "assertion failed: mismatch of global vert count")
    END IF
    IF(me%description%n_patch_edges_g /= me%commData%edges%n_glb) THEN
        CALL finish(routine, "assertion failed: mismatch of global edge count")
    END IF
  END SUBROUTINE asyncPatchData_construct

  SUBROUTINE asyncPatchData_destruct(me)
    CLASS(t_AsyncPatchData), INTENT(INOUT) :: me

    CALL me%commData%destruct()
    CALL me%t_RestartPatchData%destruct()
  END SUBROUTINE asyncPatchData_destruct

  FUNCTION toAsyncPatchData(patchData) RESULT(resultVar)
    CLASS(t_RestartPatchData), TARGET :: patchData
    TYPE(t_AsyncPatchData), POINTER :: resultVar

    resultVar => NULL()
    SELECT TYPE(patchData)
        TYPE IS(t_AsyncPatchData)
            resultVar => patchData
    END SELECT
  END FUNCTION toAsyncPatchData

  !------------------------------------------------------------------------------------------------
  !
  !> Prepare the asynchronous restart (collective call).
  !
  SUBROUTINE asyncRestartDescriptor_construct(me, modelType)
    CLASS(t_AsyncRestartDescriptor), INTENT(INOUT) :: me
    CHARACTER(*), INTENT(IN) :: modelType

    INTEGER :: jg, error
    CHARACTER(LEN=*), PARAMETER :: routine = modname//':asyncRestartDescriptor_construct'

    IF(.NOT. (my_process_is_work() .OR. my_process_is_restart())) RETURN

    IF(timers_level >= 5) CALL timer_start(timer_write_restart)

    ! ensure that our subtimers are consistently seen as subtimers
    IF(timers_level >= 7) THEN
        CALL timer_start(timer_write_restart_io)
        CALL timer_stop(timer_write_restart_io)
        CALL timer_start(timer_write_restart_communication)
        CALL timer_stop(timer_write_restart_communication)
    END IF

    CALL me%restartDescriptor_construct(modelType)
    CALL me%transferGlobalParameters()

    ! allocate patch data structure
    ALLOCATE(t_AsyncPatchData :: me%patchData(n_dom), STAT = error)
    IF(error /= SUCCESS) CALL finish(routine, "memory allocation failure")

    ! initialize the patch data structures
    DO jg = 1, n_dom
        CALL me%patchData(jg)%construct(me%modelType, jg)
    END DO

    IF(timers_level >= 5) CALL timer_stop(timer_write_restart)
  END SUBROUTINE asyncRestartDescriptor_construct

  !------------------------------------------------------------------------------------------------
  !
  !> Set patch-dependent dynamic data for asynchronous restart.
  !
  !------------------------------------------------------------------------------------------------
  !
  !> Writes all restart data into one or more files (one file per patch, collective across work processes).
  !
  SUBROUTINE asyncRestartDescriptor_writeRestart(me, this_datetime, jstep, opt_output_jfile)
    CLASS(t_AsyncRestartDescriptor), INTENT(INOUT) :: me
    TYPE(datetime), POINTER, INTENT(IN) :: this_datetime
    INTEGER, INTENT(IN) :: jstep
    INTEGER, INTENT(IN), OPTIONAL :: opt_output_jfile(:)

    INTEGER :: idx
    TYPE(t_restart_args) :: restart_args
    TYPE(t_PackedMessage) :: packedMessage
    CHARACTER(LEN=*), PARAMETER :: routine = modname//':asyncRestartDescriptor_writeRestart'

    IF(timers_level >= 5) CALL timer_start(timer_write_restart)

#ifdef DEBUG
    WRITE (nerr,FORMAT_VALS3)routine,' p_pe=',p_pe
#endif

    ! check kind of process
    IF(.NOT. my_process_is_work()) RETURN

    CALL compute_wait_for_restart()
    CALL restart_args%construct(this_datetime, jstep, me%modelType, opt_output_jfile)
    CALL compute_start_restart(restart_args, me%patchData, packedMessage)
    CALL restart_args%destruct()

    ! do the restart output
    DO idx = 1, SIZE(me%patchData)
      ! collective call to write the restart variables
      IF(me%patchData(idx)%description%l_dom_active) CALL compute_write_var_list(me%patchData(idx))
    END DO
    ! make sure all buffers are ready
    CALL p_barrier(comm=p_comm_work)
    IF (p_pe_work == 0) THEN
      CALL packedMessage%send(p_restart_pe0, 0, process_mpi_all_comm)
      CALL packedMessage%destruct()
    END IF

    IF(timers_level >= 5) CALL timer_stop(timer_write_restart)
  END SUBROUTINE asyncRestartDescriptor_writeRestart

  SUBROUTINE restart_write_patch(restart_args, patchData, restartAttributes, modelType)
    TYPE(t_restart_args), INTENT(IN) :: restart_args
    CLASS(t_RestartPatchData), INTENT(INOUT) :: patchData
    TYPE(t_RestartAttributeList), INTENT(INOUT) :: restartAttributes
    CHARACTER(*), INTENT(IN) :: modelType

    CHARACTER(*), PARAMETER :: routine = modname//":restart_write_patch"
#ifdef DEBUG
    TYPE(t_NamelistArchive), POINTER :: namelists
    namelists => namelistArchive()
#endif

    ! check if the patch is active
    IF(.NOT. patchData%description%l_dom_active) RETURN

    ! consider the right restart process
    IF(p_pe - p_restart_pe0 == MOD(patchData%description%id-1, process_mpi_restart_size)) THEN
#ifdef DEBUG
        CALL print_restart_arguments(restart_args, patchData)
        CALL restartAttributes%printAttributes()
        CALL namelists%print()
#endif
        CALL patchData%writeFile(restartAttributes, restart_args, .TRUE.)
    END IF
  END SUBROUTINE restart_write_patch

  !> Writes all restart data into one or more files (one file per patch, collective across restart processes).
  SUBROUTINE asyncRestartDescriptor_restartWriteAsyncRestart(me, restart_args)
    CLASS(t_AsyncRestartDescriptor), INTENT(INOUT) :: me
    TYPE(t_restart_args), INTENT(IN) :: restart_args

    INTEGER :: idx
    TYPE(t_RestartAttributeList), POINTER :: restartAttributes
    CHARACTER(LEN=*), PARAMETER :: routine = modname//':asyncRestartDescriptor_restartWriteAsyncRestart'

#ifdef DEBUG
    WRITE (nerr,FORMAT_VALS3)routine,' p_pe=',p_pe
#endif

    ! check kind of process
    IF(.NOT. my_process_is_restart()) CALL finish(routine, "assertion failed: "//routine//"() called by non-restart process")

    restartAttributes => RestartAttributeList_make()
    CALL me%defineRestartAttributes(restartAttributes, restart_args)

    ! do the restart output
    DO idx = 1, SIZE(me%patchData)
        CALL restart_write_patch(restart_args, me%patchData(idx), restartAttributes, me%modelType)
    END DO

    CALL restartAttributes%destruct()
    DEALLOCATE(restartAttributes)
  END SUBROUTINE asyncRestartDescriptor_restartWriteAsyncRestart

  !------------------------------------------------------------------------------------------------
  !
  !> Closes asynchronous restart (collective call).
  !
  SUBROUTINE asyncRestartDescriptor_destruct(me)
    CLASS(t_AsyncRestartDescriptor), INTENT(INOUT) :: me

    INTEGER :: i
    CHARACTER(LEN=*), PARAMETER :: routine = modname//':asyncRestartDescriptor_destruct'

    IF(timers_level >= 5) CALL timer_start(timer_write_restart)

    ! check kind of process
    IF (.NOT. my_process_is_work() .AND. .NOT. my_process_is_restart()) RETURN

    IF (my_process_is_work()) THEN

      CALL compute_wait_for_restart()
      CALL compute_shutdown_restart

    ENDIF

    DO i = 1, SIZE(me%patchData, 1)
        CALL me%patchData(i)%destruct()
    END DO
    DEALLOCATE(me%patchData)

    IF(timers_level >= 5) CALL timer_stop(timer_write_restart)
  END SUBROUTINE asyncRestartDescriptor_destruct

  !-------------------------------------------------------------------------------------------------
  !>
  !! Main routine for restart PEs.
  SUBROUTINE asyncRestart_mainLoop
    LOGICAL :: done
    TYPE(t_AsyncRestartDescriptor) :: restartDescriptor
    TYPE(t_restart_args) :: restart_args
    CHARACTER(LEN=*), PARAMETER :: routine = modname//':asyncRestart_mainLoop'

#ifdef DEBUG
    WRITE (nerr,FORMAT_VALS3)routine,' p_pe=',p_pe
#endif

    ! check kind of process
    IF (.NOT. my_process_is_restart()) RETURN

    ! prepare restart (collective call)
    CALL restartDescriptor%construct('')    ! the model TYPE will be supplied by the work processes

    ! tell the compute PEs that we are ready to work
    CALL restart_send_ready

    ! enter restart loop
    done = .FALSE.
    DO
      ! wait for a message from the compute PEs to start
      CALL restart_wait_for_start(restart_args, restartDescriptor, done)   ! this constructs the restart_args

      IF(done) EXIT ! leave loop, we are done

      ! read and write restart variable lists (collective call)
      IF(timers_level >= 5) CALL timer_start(timer_write_restart)
      CALL restartDescriptor%restartWriteAsyncRestart(restart_args)
      CALL restart_args%destruct()
      IF(timers_level >= 5) CALL timer_stop(timer_write_restart)

      ! inform compute PEs that the restart is done
      CALL restart_send_ready
    ENDDO

    ! finalization sequence (collective call)
    CALL restartDescriptor%destruct()
  END SUBROUTINE asyncRestart_mainLoop

  !-------------------------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------------------------
  ! Flow control routines between compute and restart procs ...
  !-------------------------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------------------------
  ! ... called on restart procs:
  !-------------------------------------------------------------------------------------------------
  !>
  !! restart_send_ready: Send a message to the compute PEs that the restart is ready.
  !! The counterpart on the compute side is compute_wait_for_restart.
  !
  SUBROUTINE restart_send_ready

    REAL(wp) :: msg

#ifdef DEBUG
    CHARACTER(LEN=*), PARAMETER :: routine = modname//':restart_send_ready'

    WRITE (nerr,FORMAT_VALS5)routine,' p_pe=',p_pe, &
      & ' call p_barrier with communicator=',p_comm_work
#endif
    ! make sure all are done
    CALL p_barrier(comm=p_comm_work)

    ! simply send a message from restart PE 0 to compute PE 0
    IF (p_pe_work == 0) THEN
      msg = REAL(MSG_RESTART_DONE, wp)
#ifdef DEBUG
      WRITE (nerr,FORMAT_VALS7)routine,' p_pe=',p_pe, &
        & ' send message=',INT(msg),' to pe=',p_work_pe0
#endif
      CALL p_send(msg, p_work_pe0, 0)
    ENDIF

  END SUBROUTINE restart_send_ready

  SUBROUTINE restartMetadataPacker(operation, restart_args, patchData, packedMessage)
    INTEGER, VALUE :: operation
    TYPE(t_restart_args), INTENT(INOUT) :: restart_args
    CLASS(t_RestartPatchData), INTENT(INOUT) :: patchData(:)
    TYPE(t_PackedMessage), INTENT(INOUT) :: packedMessage

    INTEGER :: i
    CHARACTER(LEN = *), PARAMETER :: routine = modname//":restartMetadataPacker"

    ! (un)pack patch independent arguments
    CALL restart_args%packer(operation, packedMessage)

    ! (un)pack the patch descriptions
    DO i = 1, SIZE(patchData)
        CALL patchData(i)%description%packer(operation, packedMessage)
    END DO
  END SUBROUTINE restartMetadataPacker

  !-------------------------------------------------------------------------------------------------
  !>
  !! restart_wait_for_start: Wait for a message from compute PEs that we should start restart or finish.
  !! The counterpart on the compute side is compute_start_restart/compute_shutdown_restart.
  !
  SUBROUTINE restart_wait_for_start(restart_args, restartDescriptor, done)
    TYPE(t_restart_args), INTENT(INOUT) :: restart_args
    CLASS(t_AsyncRestartDescriptor), INTENT(INOUT) :: restartDescriptor
    LOGICAL, INTENT(OUT)           :: done ! flag if we should shut down

    INTEGER :: iheader
    TYPE(t_PackedMessage) :: packedMessage
    CHARACTER(LEN=*), PARAMETER :: routine = modname//':restart_wait_for_start'

#ifdef DEBUG
    WRITE (nerr,FORMAT_VALS3)routine,' is called, p_pe=',p_pe
#endif

    IF(timers_level >= 7) CALL timer_start(timer_write_restart_communication)

    CALL packedMessage%construct()

    ! receive message that we may start restart (or should finish)
    IF(p_pe_work == 0) CALL packedMessage%recv(p_work_pe0, 0, process_mpi_all_comm)
    CALL packedMessage%bcast(0, p_comm_work)

    ! unpack AND interpret the message
    CALL packedMessage%unpack(iheader)
    SELECT CASE(iheader)
      CASE(MSG_RESTART_START)
        CALL restartMetadataPacker(kUnpackOp, restart_args, restartDescriptor%patchData, packedMessage)
        done = .FALSE.

      CASE(MSG_RESTART_SHUTDOWN)
        done = .TRUE.

      CASE DEFAULT
        ! anything else is an error
        CALL finish(routine,'restart PE: Got illegal restart tag')

    END SELECT

    CALL packedMessage%destruct() ! cleanup

    IF(timers_level >= 7) CALL timer_stop(timer_write_restart_communication)
  END SUBROUTINE restart_wait_for_start

  !-------------------------------------------------------------------------------------------------
  ! ... called on compute procs:
  !-------------------------------------------------------------------------------------------------
  !>
  !! compute_wait_for_restart: Wait for a message that the restart is ready.
  !! The counterpart on the restart side is restart_send_ready.
  !
  SUBROUTINE compute_wait_for_restart

    REAL(wp) :: msg

    CHARACTER(LEN=*), PARAMETER :: routine = modname//':compute_wait_for_restart'

#ifdef DEBUG
    WRITE (nerr,FORMAT_VALS3)routine,' is called, p_pe=',p_pe
#endif

    IF(timers_level >= 7) CALL timer_start(timer_write_restart_communication)

    ! first compute PE receives message from restart leader
    IF(p_pe_work == 0) THEN
      CALL message(routine, "Waiting for restart to finish (might take a while)")
      CALL p_recv(msg, p_restart_pe0, 0)
#ifdef DEBUG
      WRITE (nerr,FORMAT_VALS7)routine,' p_pe=',p_pe, &
        & ' p_recv got msg=',INT(msg),' from pe=',p_restart_pe0
#endif
      ! just for safety: Check if we got the correct tag
      IF(INT(msg) /= MSG_RESTART_DONE) THEN
        CALL finish(routine,'Compute PE: Got illegal restart tag')
      ENDIF
    ENDIF

#ifdef DEBUG
    WRITE (nerr,FORMAT_VALS5)routine,' p_pe=',p_pe, &
      & ' call p_barrier with communicator=',p_comm_work
#endif
    ! wait in barrier until message is here
    CALL p_barrier(comm=p_comm_work)

    IF(timers_level >= 7) CALL timer_stop(timer_write_restart_communication)
  END SUBROUTINE compute_wait_for_restart

  !-------------------------------------------------------------------------------------------------
  !>
  !! compute_start_restart: Send a message to restart PEs that they should start restart.
  !! The counterpart on the restart side is restart_wait_for_start.
  !
  SUBROUTINE compute_start_restart(restart_args, patchData, packedMessage)
    TYPE(t_restart_args), INTENT(INOUT) :: restart_args
    CLASS(t_RestartPatchData), INTENT(INOUT) :: patchData(:)
    TYPE(t_PackedMessage), INTENT(inout) :: packedMessage

    INTEGER :: i, trash
    CHARACTER(LEN=*), PARAMETER :: routine = modname//':compute_start_restart'

#ifdef DEBUG
    WRITE (nerr,FORMAT_VALS5)routine,' p_pe=',p_pe, ' call p_barrier with communicator=',p_comm_work
#endif

    IF(timers_level >= 7) CALL timer_start(timer_write_restart_communication)

    ! if processor splitting is applied, the time-dependent data need to be transferred
    ! from the subset master PE to PE0, from where they are communicated to the output PE(s)
    DO i = 1, SIZE(patchData)
      CALL patchData(i)%description%updateOnMaster()
    END DO

    CALL packedMessage%construct()
    IF(p_pe_work == 0) THEN
      ! send the DATA to the restart master
      CALL packedMessage%pack(MSG_RESTART_START)  ! set command id
      CALL restartMetadataPacker(kPackOp, restart_args, patchData, packedMessage)    ! all the other DATA
    ENDIF
    ! broadcast a copy among the compute processes, so that all
    ! processes have a consistent view of the restart patch
    ! descriptions
    CALL packedMessage%bcast(0, p_comm_work)
    CALL packedMessage%unpack(trash)  ! ignore the command id
    CALL restartMetadataPacker(kUnpackOp, restart_args, patchData, packedMessage)

    IF (p_pe_work /= 0) CALL packedMessage%destruct()

    IF(timers_level >= 7) CALL timer_stop(timer_write_restart_communication)
  END SUBROUTINE compute_start_restart

  !-------------------------------------------------------------------------------------------------
  !>
  !! compute_shutdown_restart: Send a message to restart PEs that they should shut down.
  !! The counterpart on the restart side is restart_wait_for_start.
  !
  SUBROUTINE compute_shutdown_restart
    TYPE(t_PackedMessage) :: packedMessage
    CHARACTER(LEN=*), PARAMETER :: routine = modname//':compute_shutdown_restart'

#ifdef DEBUG
    WRITE (nerr,FORMAT_VALS5)routine,' p_pe=',p_pe, &
      & ' call p_barrier with communicator=',p_comm_work
#endif

    IF(timers_level >= 7) CALL timer_start(timer_write_restart_communication)

    ! make sure all are here
    CALL p_barrier(comm=p_comm_work)

    IF(p_pe_work == 0) THEN
      CALL packedMessage%construct()  ! create message array

      ! sent the shutdown message
      CALL packedMessage%pack(MSG_RESTART_SHUTDOWN)
      CALL packedMessage%send(p_restart_pe0, 0, process_mpi_all_comm)

      CALL packedMessage%destruct()   ! cleanup
    ENDIF

    IF(timers_level >= 7) CALL timer_stop(timer_write_restart_communication)
  END SUBROUTINE compute_shutdown_restart

  !------------------------------------------------------------------------------------------------
  !
  ! Common helper routines for error handling and releasing of resources.
  !
  !------------------------------------------------------------------------------------------------
  !
  !  A simple helper routine to check the result of the last MPI call.
  !
  SUBROUTINE check_mpi_error(routine, mpi_call, mpi_error, l_finish)

    CHARACTER(LEN=*), INTENT(IN)    :: routine, mpi_call
    INTEGER, INTENT(IN)             :: mpi_error
    LOGICAL, INTENT(IN)             :: l_finish

    CHARACTER(LEN = 256) :: error_message

    IF (mpi_error /= MPI_SUCCESS) THEN
      IF (l_finish) THEN
        WRITE (error_message, '(2a,i5)')TRIM(mpi_call), &
          &                    ' returned with error=',mpi_error
        CALL finish(routine, TRIM(error_message))
      ELSE
        WRITE (error_message, '(4a,i5)')TRIM(routine), ".", TRIM(mpi_call), &
          &                    ' returned with error=',mpi_error
        WRITE (nerr, TRIM(error_message))
      ENDIF
    ENDIF

  END SUBROUTINE check_mpi_error

  !------------------------------------------------------------------------------------------------
  !
  ! Common helper routines to processing lists.
  !

  !------------------------------------------------------------------------------------------------
  !
  !  Print restart arguments.
  !
  SUBROUTINE print_restart_arguments(restart_args, patchData)
    TYPE(t_restart_args), INTENT(IN) :: restart_args
    CLASS(t_RestartPatchData), INTENT(IN) :: patchData
    CHARACTER(LEN=*), PARAMETER   :: routine = modname//':print_restart_arguments'

    WRITE (nerr,FORMAT_VALS3)routine,' is called for p_pe=',p_pe

    !LK removed because arguments do not exist anymore:    CALL restart_args%print(routine//": ")

    ! patch informations
    PRINT *,routine, ' SIZE(patchData%description%opt_pvct) = ', SIZE(patchData%description%opt_pvct)
    PRINT *,routine, ' SIZE(patchData%description%opt_t_elapsed_phy) =', SIZE(patchData%description%opt_t_elapsed_phy)
  END SUBROUTINE print_restart_arguments

  !------------------------------------------------------------------------------------------------
  !
  ! Common helper routines to preparing the restart.
  !
  !-------------------------------------------------------------------------------------------------

  !------------------------------------------------------------------------------------------------
  !
  ! Write restart variable list for a restart PE.
  !
  SUBROUTINE asyncPatchData_writeData(me, file)
    CLASS(t_AsyncPatchData), INTENT(INOUT) :: me
    TYPE(t_RestartFile), INTENT(INOUT) :: file

    INTEGER                         :: iv, nval, ierrstat, nlevs, ilev, pointCount
    INTEGER(KIND=MPI_ADDRESS_KIND)  :: ioff(0:num_work_procs-1)
    REAL(dp), ALLOCATABLE           :: buffer_dp(:,:)
    REAL(sp), ALLOCATABLE           :: buffer_sp(:,:)
    INTEGER                         :: ichunk, nchunks, chunk_start, chunk_end

    ! For timing
    REAL(dp) :: t_get, t_write
    INTEGER(i8) :: bytesGet, bytesWrite
    LOGICAL :: flag_dp

    CHARACTER(LEN=*), PARAMETER     :: routine = modname//':asyncPatchData_writeData'

#ifdef DEBUG
    WRITE (nerr,FORMAT_VALS3)routine,' p_pe=',p_pe
#endif

    ! check process
    IF (.NOT. my_process_is_restart()) CALL finish(routine, 'Must be called on a restart PE!')

    t_get   = 0.d0
    t_write = 0.d0
    bytesGet = 0_i8
    bytesWrite = 0_i8

    ! check the contained array of restart variables
    IF (.NOT. ASSOCIATED(me%varData)) RETURN

    ! get maximum number of data points in a slice and allocate tmp. variables
    nval = me%commData%maxLevelSize()

    ! allocate RMA memory
    !
    ! TODO[FP] avoid allocating two of these buffers , e.g. by looping
    ! twice over the variables, first over double precision, then over
    ! single precision.
    ! 
    ALLOCATE(buffer_dp(nval,restart_chunk_size), &
      &      buffer_sp(nval,restart_chunk_size), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, "memory allocation failure")

    ioff(:) = 0

    ! go over the all restart variables in the associated array
    VAR_LOOP : DO iv = 1, SIZE(me%varData)

#ifdef DEBUG
      WRITE (nerr,FORMAT_VALS5I) routine,' p_pe=',p_pe,' restart pe processes field=', &
        &                        TRIM(me%varData(iv)%info%name)
#endif

      ! check time level of the field
      IF (.NOT. has_valid_time_level(me%varData(iv)%info, me%description%id, &
        &                            me%description%nnew, me%description%nnew_rcf)) CYCLE

      ! get current level
      IF(me%varData(iv)%info%ndims == 2) THEN
        nlevs = 1
      ELSE
        nlevs = me%varData(iv)%info%used_dimensions(2)
      ENDIF

      ! get pointer to reorder data
      IF(me%varData(iv)%info%hgrid == GRID_UNSTRUCTURED_CELL) pointCount = me%description%n_patch_cells_g
      IF(me%varData(iv)%info%hgrid == GRID_UNSTRUCTURED_VERT) pointCount = me%description%n_patch_verts_g
      IF(me%varData(iv)%info%hgrid == GRID_UNSTRUCTURED_EDGE) pointCount = me%description%n_patch_edges_g

      ! check if this is single or double precision:
      SELECT CASE(me%varData(iv)%getDatatype())
      CASE(REAL_T, INT_T)
        ! INTEGER fields: we write them as REAL-valued arrays
        flag_dp = .TRUE.
      CASE(SINGLE_T)
        flag_dp = .FALSE.
      CASE DEFAULT
        CALL finish(routine, "Internal error! Variable "//TRIM(me%varData(iv)%info%name))
      END SELECT

      ! no. of chunks of levels (each of size "restart_chunk_size"):
      nchunks = (nlevs-1)/restart_chunk_size + 1
      ! loop over all chunks (of levels)
      LEVELS : DO ichunk=1,nchunks
        chunk_start = (ichunk-1)*restart_chunk_size + 1
        chunk_end = MIN(chunk_start+restart_chunk_size-1, nlevs)
        IF (flag_dp) THEN
          CALL me%commData%collectData(me%varData(iv)%info%hgrid, chunk_end - chunk_start + 1, &
            &                          buffer_dp, ioff, t_get, bytesGet)
        ELSE
          CALL me%commData%collectData(me%varData(iv)%info%hgrid, chunk_end - chunk_start + 1, &
            &                          buffer_sp, ioff, t_get, bytesGet)
        END IF

        ! write field content into a file
        t_write = t_write - p_mpi_wtime()
        IF (flag_dp) THEN
          DO ilev=chunk_start, chunk_end
            CALL file%writeLevel(me%varData(iv)%info%cdiVarID, (ilev-1), buffer_dp(:, ilev - chunk_start + 1))
            bytesWrite = bytesWrite + pointCount*BYTES_DP
          END DO
        ELSE
          ! single precision data
          DO ilev=chunk_start, chunk_end
            CALL file%writeLevel(me%varData(iv)%info%cdiVarID, (ilev-1), buffer_sp(:, ilev - chunk_start + 1))
            bytesWrite = bytesWrite + pointCount*BYTES_SP
          END DO
        END IF
        t_write = t_write + p_mpi_wtime()

      ENDDO LEVELS

#ifdef DEBUG
      WRITE(nerr, FORMAT_VALS7I) routine, ' p_pe=', p_pe, ' restart pe writes field=', &
        &                        TRIM(me%varData(iv)%info%name), ' data=', pointCount*nlevs
#endif
    ENDDO VAR_LOOP

    IF (ALLOCATED(buffer_dp)) THEN
      DEALLOCATE(buffer_dp, STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed!')
    END IF
    IF (ALLOCATED(buffer_sp)) THEN
      DEALLOCATE(buffer_sp, STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed!')
    END IF

    IF (msg_level >= 12) THEN
      WRITE (0,'(10(a,f10.3))') ' Restart: Got ', REAL(bytesGet, dp)*1.d-6, ' MB, time get: ', t_get, ' s [', &
           & REAL(bytesGet, dp)*1.d-6/MAX(1.e-6_wp, t_get), ' MB/s], time write: ', t_write, ' s [', &
           & REAL(bytesWrite, dp)*1.d-6/MAX(1.e-6_wp,t_write), ' MB/s]'
    ENDIF
  END SUBROUTINE asyncPatchData_writeData


  !------------------------------------------------------------------------------------------------
  !
  ! Write restart variable lists for a compute PE.
  !
  SUBROUTINE compute_write_var_list(patchData)
    CLASS(t_RestartPatchData), TARGET, INTENT(INOUT) :: patchData

    TYPE(t_AsyncPatchData), POINTER :: asyncPatchData
    TYPE(t_RestartVarData), POINTER :: p_vars(:)
    TYPE(t_ptr_2d),     ALLOCATABLE :: dataPointers_dp(:)
    TYPE(t_ptr_2d_sp) , ALLOCATABLE :: dataPointers_sp(:)
    TYPE(t_ptr_2d_int), ALLOCATABLE :: dataPointers_int(:)
    INTEGER                         :: iv
    INTEGER(i8)                     :: offset
    CHARACTER(LEN=*), PARAMETER     :: routine = modname//':compute_write_var_list'

#ifdef DEBUG
    WRITE (nerr,FORMAT_VALS3)routine,' p_pe=',p_pe
#endif

    ! dynamic downcast of patchData argument
    asyncPatchData => toAsyncPatchData(patchData)
    IF(.NOT.ASSOCIATED(asyncPatchData)) CALL finish(routine, "assertion failed: wrong type of patchData")

    ! check the array of restart variables
    p_vars => asyncPatchData%varData
    IF (.NOT. ASSOCIATED(p_vars)) RETURN

    ! offset in RMA window for async restart
    offset = 0

    ! go over the all restart variables in the associated array
    DO iv = 1, SIZE(p_vars)
#ifdef DEBUG
        WRITE (nerr,FORMAT_VALS5I)routine,' p_pe=',p_pe,' compute pe processes field=',TRIM(p_vars(iv)%info%name)
#endif

        ! check time level of the field
        IF (has_valid_time_level(p_vars(iv)%info, asyncPatchData%description%id, &
          &                      asyncPatchData%description%nnew,                &
          &                      asyncPatchData%description%nnew_rcf)) THEN

          SELECT CASE(p_vars(iv)%getDatatype())
          CASE(REAL_T)
            CALL getLevelPointers(p_vars(iv)%info, p_vars(iv)%r_ptr, dataPointers_dp)
            CALL asyncPatchData%commData%postData(p_vars(iv)%info%hgrid, dataPointers_dp, offset)
          CASE(SINGLE_T)
            CALL getLevelPointers(p_vars(iv)%info, p_vars(iv)%s_ptr, dataPointers_sp)
            CALL asyncPatchData%commData%postData(p_vars(iv)%info%hgrid, dataPointers_sp, offset)
          CASE(INT_T)
            CALL getLevelPointers(p_vars(iv)%info, p_vars(iv)%i_ptr, dataPointers_int)
            CALL asyncPatchData%commData%postData(p_vars(iv)%info%hgrid, dataPointers_int, offset)
          CASE DEFAULT
            CALL finish(routine, "Internal error! Variable "//TRIM(p_vars(iv)%info%name))
          END SELECT

          ! no deallocation of dataPointers, so that the next
          ! invocation of getLevelPointers() may reuse the last
          ! allocation
        END IF
    END DO
  END SUBROUTINE compute_write_var_list

#endif
END MODULE mo_async_restart
