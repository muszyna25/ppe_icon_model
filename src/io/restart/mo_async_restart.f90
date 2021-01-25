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
MODULE mo_async_restart
#ifndef NOMPI
  USE mo_async_restart_comm_data,   ONLY: t_AsyncRestartCommData
  USE mo_exception,                 ONLY: finish, message, message_text
  USE mo_kind,                      ONLY: dp, sp
  USE mtime,                        ONLY: datetime
  USE mo_impl_constants,            ONLY: SUCCESS, SINGLE_T, REAL_T, INT_T
  USE mo_grid_config,               ONLY: n_dom
  USE mo_packed_message,            ONLY: t_PackedMessage, kPackOp, kUnpackOp
  USE mo_restart_descriptor,        ONLY: t_RestartDescriptor
  USE mo_restart_patch_data,        ONLY: t_RestartPatchData
  USE mo_async_restart_patch_data,  ONLY: t_asyncPatchData, toAsyncPatchData
  USE mo_restart_patch_description, ONLY: t_restart_patch_description
  USE mo_restart_util,              ONLY: t_restart_args
  USE mo_restart_var_data,          ONLY: get_var_3d_ptr, has_valid_time_level
  USE mo_var_list_element,          ONLY: t_p_var_list_element
  USE mo_timer,                     ONLY: timer_start, timer_stop, timer_write_restart, timer_write_restart_io, &
                                        & timer_write_restart_communication, timers_level
  USE mo_mpi,                       ONLY: p_pe_work, p_comm_work, &
    &                                     p_send, p_recv, p_comm_work_2_restart, &
    &                                     p_barrier, my_process_is_restart, my_process_is_work, &
    &                                     process_mpi_restart_size, p_mpi_wtime

  IMPLICIT NONE
  PRIVATE
  ! public routines
  PUBLIC :: t_AsyncRestartDescriptor
  PUBLIC :: asyncRestart_mainLoop

  INTEGER, PARAMETER :: MSG_RESTART_START    = 111111
  INTEGER, PARAMETER :: MSG_RESTART_DONE     = 222222
  INTEGER, PARAMETER :: MSG_RESTART_SHUTDOWN = 999999
  CHARACTER(*), PARAMETER :: modname = 'mo_async_restart'

  TYPE, EXTENDS(t_RestartDescriptor) :: t_AsyncRestartDescriptor
  CONTAINS
    PROCEDURE :: construct => asyncRestartDescriptor_construct  ! override
    PROCEDURE :: writeRestart => asyncRestartDescriptor_writeRestart    ! override
    PROCEDURE :: destruct => asyncRestartDescriptor_destruct    ! override
  END TYPE t_AsyncRestartDescriptor

CONTAINS
  !
  !> Prepare the asynchronous restart (collective call).
  SUBROUTINE asyncRestartDescriptor_construct(me, modelType)
    CLASS(t_AsyncRestartDescriptor), INTENT(INOUT), TARGET :: me
    CHARACTER(*), INTENT(IN) :: modelType
    INTEGER :: jg, error
    CHARACTER(LEN=*), PARAMETER :: routine = modname//':asyncRestartDescriptor_construct'

    IF(.NOT. (my_process_is_work() .OR. my_process_is_restart())) RETURN
    IF(timers_level >= 5) CALL timer_start(timer_write_restart)
    IF(timers_level >= 7) THEN
      CALL timer_start(timer_write_restart_io)
      CALL timer_stop(timer_write_restart_io)
      CALL timer_start(timer_write_restart_communication)
      CALL timer_stop(timer_write_restart_communication)
    END IF
    me%modelType = modelType
    CALL me%transferGlobalParameters()
    ALLOCATE(t_AsyncPatchData :: me%patchData(n_dom), STAT = error)
    IF(error /= SUCCESS) CALL finish(routine, "memory allocation failure")
    DO jg = 1, n_dom
      CALL me%patchData(jg)%construct(me%modelType, jg)
    END DO
    IF(timers_level >= 5) CALL timer_stop(timer_write_restart)
  END SUBROUTINE asyncRestartDescriptor_construct

  !> Writes all restart data into one or more files (one file per patch, collective across work processes).
  SUBROUTINE asyncRestartDescriptor_writeRestart(me, this_datetime, jstep, opt_output_jfile)
    CLASS(t_AsyncRestartDescriptor), INTENT(INOUT), TARGET :: me
    TYPE(datetime), POINTER, INTENT(IN) :: this_datetime
    INTEGER, INTENT(IN) :: jstep
    INTEGER, INTENT(IN), OPTIONAL :: opt_output_jfile(:)
    INTEGER :: idx
    TYPE(t_restart_args) :: restart_args
    TYPE(t_PackedMessage) :: packedMessage
    TYPE(t_AsyncPatchData), POINTER :: patchData
    CHARACTER(LEN=*), PARAMETER :: routine = modname//':asyncRestartDescriptor_writeRestart'

    IF(timers_level >= 5) CALL timer_start(timer_write_restart)
    IF(.NOT. my_process_is_work()) RETURN
    CALL compute_wait_for_restart()
    CALL restart_args%construct(this_datetime, jstep, me%modelType, opt_output_jfile)
    CALL compute_prepare_restart()
    CALL restart_args%destruct()
    DO idx = 1, SIZE(me%patchData)
      patchData => toAsyncPatchData(me%patchData(idx))
      IF (.NOT. ASSOCIATED(patchData)) CALL finish(routine, "wrong type of patchData")
      IF (patchData%description%l_dom_active .AND. SIZE(patchData%varData) > 0) &
        & CALL compute_write_var_list(patchData%commData, patchData%description, &
        &                           patchData%varData)
    END DO
    CALL p_barrier(comm=p_comm_work)
    IF (p_pe_work == 0) & ! actual trigger for restart writing
      & CALL packedMessage%send(0, 0, p_comm_work_2_restart)
    IF(timers_level >= 5) CALL timer_stop(timer_write_restart)
  CONTAINS

  SUBROUTINE compute_wait_for_restart
    INTEGER :: msg
    CHARACTER(LEN=*), PARAMETER :: routine = modname//':compute_wait_for_restart'
    REAL(dp) :: t_wait

    IF(timers_level >= 7) CALL timer_start(timer_write_restart_communication)
    IF(p_pe_work == 0) THEN
      t_wait = p_mpi_wtime()
      CALL message(routine, "Waiting for restart servers to accomplish")
      CALL p_recv(msg, 0, 0, comm=p_comm_work_2_restart)
      t_wait = p_mpi_wtime() - t_wait
      WRITE(message_text, '(a,e10.3,a)') 'restart servers ready... waited ', &
           t_wait, ' sec'
      CALL message(routine, message_text)
      IF (msg /= MSG_RESTART_DONE) CALL finish(routine,'message content does not match')
    ENDIF
    CALL p_barrier(comm=p_comm_work)
    IF(timers_level >= 7) CALL timer_stop(timer_write_restart_communication)
  END SUBROUTINE compute_wait_for_restart

  SUBROUTINE compute_prepare_restart()
    INTEGER :: i, trash

    IF(timers_level >= 7) CALL timer_start(timer_write_restart_communication)
    DO i = 1, SIZE(me%patchData)
      CALL me%patchData(i)%description%updateOnMaster()
    END DO
    IF (p_pe_work == 0) THEN
      CALL packedMessage%pack(MSG_RESTART_START)  ! set command id
      CALL restartMetadataPacker(kPackOp, restart_args, me%patchData, packedMessage)    ! all the other DATA
    END IF
    CALL packedMessage%bcast(0, p_comm_work)
    CALL packedMessage%unpack(trash)  ! ignore the command id
    CALL restartMetadataPacker(kUnpackOp, restart_args, me%patchData, packedMessage)
    IF(timers_level >= 7) CALL timer_stop(timer_write_restart_communication)
  END SUBROUTINE compute_prepare_restart

  SUBROUTINE compute_write_var_list(commData, description, vars)
    TYPE(t_AsyncRestartCommData), INTENT(inout) :: commData
    TYPE(t_restart_patch_description), INTENT(in) :: description
    TYPE(t_p_var_list_element), INTENT(in) :: vars(:)
    REAL(dp), POINTER               :: r_ptr_3d(:,:,:)
    REAL(sp), POINTER               :: s_ptr_3d(:,:,:)
    INTEGER, POINTER                :: i_ptr_3d(:,:,:)
    INTEGER                         :: iv, offset
    CHARACTER(LEN=*), PARAMETER     :: routine = modname//':compute_write_var_list'

    offset = 0
    CALL commData%sync(.FALSE.)
    DO iv = 1, SIZE(vars)
      IF (has_valid_time_level(vars(iv)%p%info, description%id, &
        &                      description%nnew, description%nnew_rcf)) THEN
        SELECT CASE(vars(iv)%p%info%data_type)
        CASE(REAL_T)
          CALL get_var_3d_ptr(vars(iv)%p, r_ptr_3d)
          CALL commData%postData(vars(iv)%p%info%hgrid, r_ptr_3d, offset)
        CASE(SINGLE_T)
          CALL get_var_3d_ptr(vars(iv)%p, s_ptr_3d)
          CALL commData%postData(vars(iv)%p%info%hgrid, s_ptr_3d, offset)
        CASE(INT_T)
          CALL get_var_3d_ptr(vars(iv)%p, i_ptr_3d)
          CALL commData%postData(vars(iv)%p%info%hgrid, i_ptr_3d, offset)
        CASE DEFAULT
          CALL finish(routine, "Internal error! Variable "//TRIM(vars(iv)%p%info%name))
        END SELECT
      END IF
    END DO
    CALL commData%sync(.TRUE.)
  END SUBROUTINE compute_write_var_list

  END SUBROUTINE asyncRestartDescriptor_writeRestart

  !------------------------------------------------------------------------------------------------
  !
  !> Closes asynchronous restart (workers only).
  !
  SUBROUTINE asyncRestartDescriptor_destruct(me)
    CLASS(t_AsyncRestartDescriptor), INTENT(INOUT) :: me
    TYPE(t_PackedMessage) :: packedMessage
    INTEGER :: i, msg
    REAL(dp) :: t_wait
    CHARACTER(LEN=*), PARAMETER :: routine = modname//':asyncRestartDescriptor_destruct'

    IF (my_process_is_restart()) RETURN
    IF(timers_level >= 5) CALL timer_start(timer_write_restart)
    IF(timers_level >= 7) CALL timer_start(timer_write_restart_communication)
    IF(p_pe_work == 0) THEN
      t_wait = p_mpi_wtime()
      CALL message(routine, "Waiting for restart to accomplish")
      CALL p_recv(msg, 0, 0, comm=p_comm_work_2_restart)
      t_wait = p_mpi_wtime() - t_wait
      WRITE(message_text, '(a,e10.3,a)') 'restart servers ready for shutdown ', t_wait, ' sec'
      CALL message(routine, message_text)
      IF (msg /= MSG_RESTART_DONE) CALL finish(routine,'message content does not match')
    ENDIF
    CALL p_barrier(comm=p_comm_work)
    IF(p_pe_work == 0) THEN
      CALL packedMessage%pack(MSG_RESTART_SHUTDOWN)
      CALL packedMessage%send(0, 0, p_comm_work_2_restart)
      CALL packedMessage%reset()
    ENDIF
    IF(timers_level >= 7) CALL timer_stop(timer_write_restart_communication)
    DO i = 1, SIZE(me%patchData, 1)
      CALL me%patchData(i)%destruct()
    END DO
    DEALLOCATE(me%patchData)
    IF(timers_level >= 5) CALL timer_stop(timer_write_restart)
  END SUBROUTINE asyncRestartDescriptor_destruct

  !-------------------------------------------------------------------------------------------------
  !>
  !! Main routine for restart PEs.
  SUBROUTINE asyncRestart_mainLoop()
    LOGICAL :: at_service
    TYPE(t_AsyncRestartDescriptor), TARGET :: desc
    TYPE(t_restart_args) :: restart_args
    INTEGER :: j

    IF (.NOT. my_process_is_restart()) RETURN
    CALL desc%construct('')
    IF (p_pe_work == 0) CALL p_send(MSG_RESTART_DONE, 0, 0, comm=p_comm_work_2_restart)
    at_service = .true.
    DO WHILE (at_service)
      CALL restart_wait_for_trigger()   ! this constructs the restart_args
      IF (.NOT.at_service) CYCLE
      IF(timers_level >= 5) CALL timer_start(timer_write_restart)
      CALL restart_sync(.FALSE.)
      CALL desc%writeFiles(restart_args, isSync=.FALSE.)
      CALL restart_sync(.TRUE.)
      CALL restart_args%destruct()
      CALL p_barrier(comm=p_comm_work)
      IF (p_pe_work == 0) &
        & CALL p_send(MSG_RESTART_DONE, 0, 0, comm=p_comm_work_2_restart)
      IF(timers_level >= 5) CALL timer_stop(timer_write_restart)
    ENDDO
    DO j = 1, SIZE(desc%patchData, 1)
      CALL desc%patchData(j)%destruct()
    END DO
    DEALLOCATE(desc%patchData)
  CONTAINS

  SUBROUTINE restart_wait_for_trigger()
    TYPE(t_PackedMessage) :: packedMessage
    CHARACTER(LEN=*), PARAMETER :: routine = modname//':restart_wait_for_start'
    INTEGER :: steer

    IF(timers_level >= 7) CALL timer_start(timer_write_restart_communication)
    IF (p_pe_work == 0) CALL packedMessage%recv(0, 0, p_comm_work_2_restart)
    CALL packedMessage%bcast(0, p_comm_work)
    CALL packedMessage%unpack(steer)
    SELECT CASE(steer)
    CASE(MSG_RESTART_START)
      CALL restartMetadataPacker(kUnpackOp, restart_args, desc%patchData, packedMessage)
    CASE(MSG_RESTART_SHUTDOWN)
      at_service = .false.
    CASE DEFAULT
      CALL finish(routine,'command code not recognized')
    END SELECT
    IF(timers_level >= 7) CALL timer_stop(timer_write_restart_communication)
  END SUBROUTINE restart_wait_for_trigger

  SUBROUTINE restart_sync(mode)
    LOGICAL, INTENT(IN) :: mode
    TYPE(t_asyncPatchData), POINTER :: apData
    INTEGER :: i

    DO i = 1, SIZE(desc%patchData)
      apData => toAsyncPatchData(desc%patchData(i))
      IF(p_pe_work == MOD(apData%description%id-1, process_mpi_restart_size)) THEN
        CALL apData%commData%sync(mode)
      ELSE IF (.NOT.mode) THEN
        CALL apData%commData%sync(.FALSE.)
        CALL apData%commData%sync(.TRUE.)
      END IF
    END DO
  END SUBROUTINE restart_sync

  END SUBROUTINE asyncRestart_mainLoop

  SUBROUTINE restartMetadataPacker(operation, restart_args, patchData, packedMessage)
    INTEGER, INTENT(IN) :: operation
    TYPE(t_restart_args), INTENT(INOUT) :: restart_args
    CLASS(t_RestartPatchData), INTENT(INOUT) :: patchData(:)
    TYPE(t_PackedMessage), INTENT(INOUT) :: packedMessage
    INTEGER :: i

    CALL restart_args%packer(operation, packedMessage)
    DO i = 1, SIZE(patchData)
      CALL patchData(i)%description%packer(operation, packedMessage)
    END DO
  END SUBROUTINE restartMetadataPacker

#endif
END MODULE mo_async_restart
