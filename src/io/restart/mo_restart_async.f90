!>
!! Contains routines for asynchronous restart Output
!! --------------------------------------------------------
!!
!! Note: The synchronous implementation of the restart output can be
!!       found in the module "mo_restart". See module header for
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

MODULE mo_restart_async

  USE mo_decomposition_tools,     ONLY: t_grid_domain_decomp_info
  USE mo_exception,               ONLY: finish
  USE mo_fortran_tools,           ONLY: t_ptr_2d, ensureSize
  USE mo_kind,                    ONLY: wp, i8, dp
  USE mo_datetime,                ONLY: t_datetime
  USE mo_io_units,                ONLY: nerr, filename_max
  USE mo_var_list,                ONLY: nvar_lists, var_lists, new_var_list, delete_var_lists
  USE mo_linked_list,             ONLY: t_list_element, t_var_list
  USE mo_restart_attributes,   ONLY: t_RestartAttributeList, RestartAttributeList_make
  USE mo_dynamics_config,         ONLY: iequations
  USE mo_grid_config,             ONLY: l_limited_area
  USE mo_impl_constants,          ONLY: IHS_ATM_TEMP, IHS_ATM_THETA, ISHALLOW_WATER, INH_ATMOSPHERE, LEAPFROG_EXPL, LEAPFROG_SI, &
                                      & SUCCESS, TLEV_NNOW, TLEV_NNOW_RCF
  USE mo_var_metadata_types,      ONLY: t_var_metadata
  USE mo_restart_namelist,     ONLY: print_restart_name_lists, RestartNamelist_bcast
  USE mo_communication,           ONLY: idx_no, blk_no
  USE mo_parallel_config,         ONLY: restart_chunk_size
  USE mo_grid_config,             ONLY: n_dom
  USE mo_run_config,              ONLY: msg_level
  USE mo_ha_dyn_config,           ONLY: ha_dyn_config
  USE mo_model_domain,            ONLY: p_patch, t_patch
  USE mo_cdi,                     ONLY: CDI_UNDEFID, FILETYPE_NC2, FILETYPE_NC4, streamWriteVarSlice
  USE mo_cdi_constants,           ONLY: GRID_UNSTRUCTURED_EDGE, GRID_UNSTRUCTURED_VERT, GRID_UNSTRUCTURED_CELL
  USE mo_packed_message,          ONLY: t_PackedMessage, kPackOp, kUnpackOp
  USE mo_restart_patch_description, ONLY: t_restart_patch_description
  USE mo_util_restart,            ONLY: t_restart_cdi_ids, setGeneralRestartAttributes, create_restart_file_link, t_var_data, &
                                      & getRestartFilename, t_restart_args, getLevelPointers

#ifndef NOMPI
  USE mo_mpi,                     ONLY: p_pe, p_pe_work, p_restart_pe0, p_comm_work, p_work_pe0, num_work_procs, MPI_SUCCESS, &
                                      & stop_mpi, p_send, p_recv, p_barrier, p_bcast, my_process_is_restart, my_process_is_work, &
                                      & p_comm_work_2_restart, p_n_work, p_int, process_mpi_restart_size, p_int_i8, p_real_dp, &
                                      & p_comm_work_restart, p_mpi_wtime, process_mpi_all_comm, p_get_bcast_role

  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_ptr, c_intptr_t, c_f_pointer

#ifdef __SUNPRO_F95
  INCLUDE "mpif.h"
#else
  USE mpi,                        ONLY: MPI_ADDRESS_KIND, MPI_INFO_NULL, MPI_ROOT, MPI_LOCK_SHARED, MPI_MODE_NOCHECK, &
                                      & MPI_PROC_NULL, MPI_WIN_NULL, MPI_LOCK_EXCLUSIVE
#endif
#endif

  IMPLICIT NONE

  INCLUDE 'netcdf.inc'

  ! public routines
  PUBLIC :: t_restart_descriptor
  PUBLIC :: restart_main_proc

  PRIVATE

  ! tags for communication between compute PEs and restart PEs
  INTEGER, PARAMETER :: MSG_RESTART_START         = 111111
  INTEGER, PARAMETER :: MSG_RESTART_DONE          = 222222
  INTEGER, PARAMETER :: MSG_RESTART_SHUTDOWN      = 999999

  ! maximum text lengths in this module
  INTEGER, PARAMETER :: MAX_NAME_LENGTH           = 128
  INTEGER, PARAMETER :: MAX_ERROR_LENGTH          = 256

  ! common constant strings
  CHARACTER(LEN=*), PARAMETER :: modname                  = 'mo_restart_async'
  CHARACTER(LEN=*), PARAMETER :: ALLOCATE_FAILED          = 'ALLOCATE failed!'
  CHARACTER(LEN=*), PARAMETER :: DEALLOCATE_FAILED        = 'DEALLOCATE failed!'
  CHARACTER(LEN=*), PARAMETER :: UNKNOWN_GRID_TYPE        = 'Unknown grid type!'
  CHARACTER(LEN=*), PARAMETER :: UNKNOWN_FILE_FORMAT      = 'Unknown file format for restart file.'
  CHARACTER(LEN=*), PARAMETER :: ASYNC_RESTART_REQ_MPI    = 'asynchronous restart can only run with MPI!'
  CHARACTER(LEN=*), PARAMETER :: NO_COMPUTE_PE            = 'Must be called on a compute PE!'
  CHARACTER(LEN=*), PARAMETER :: NO_RESTART_PE            = 'Must be called on a restart PE!'
  CHARACTER(LEN=*), PARAMETER :: NET_CDF_ERROR_FORMAT     = '(a,i5,a)'

  CHARACTER(LEN=*), PARAMETER :: FORMAT_VALS3             = '(a,a,i3)'
  CHARACTER(LEN=*), PARAMETER :: FORMAT_VALS5             = '(a,a,i3,a,i3)'
  CHARACTER(LEN=*), PARAMETER :: FORMAT_VALS5I            = '(a,a,i3,a,a)'
  CHARACTER(LEN=*), PARAMETER :: FORMAT_VALS7             = '(a,a,i3,a,i6,a,i3)'
  CHARACTER(LEN=*), PARAMETER :: FORMAT_VALS7I            = '(a,a,i3,a,a,a,i8)'

  ! TYPE t_restart_file
  TYPE t_restart_file
    ! the following data can be set before opening the restart file
    TYPE(t_var_data), POINTER   :: var_data(:)

    ! the following members are set during open
    CHARACTER(LEN=filename_max) :: filename
    CHARACTER(LEN=32)           :: model_type
    TYPE(t_restart_cdi_ids)     :: cdiIds
  CONTAINS
    PROCEDURE, PRIVATE :: construct => restartFile_construct
  END TYPE t_restart_file

  !------------------------------------------------------------------------------------------------
  ! TYPE t_reorder_data describes how local cells/edges/verts
  ! have to be reordered to get the global array.
  ! Below, "points" refers to either cells, edges or verts.
  !
  TYPE t_reorder_data
    INTEGER :: n_glb  ! Global number of points per logical patch
    INTEGER :: n_own  ! Number of own points (without halo, only belonging to logical patch). Only set on compute PEs, set to 0 on restart PEs.
    INTEGER, ALLOCATABLE :: own_idx(:), own_blk(:)  ! idx and blk for own points, only set on compute PEs
    INTEGER, ALLOCATABLE :: pe_own(:)   ! n_own, gathered for all compute PEs (set on all PEs)

    ! Index how to reorder the contributions of all compute PEs into the global array (set on all PEs)
    INTEGER, ALLOCATABLE :: inverse_reorder_index(:)    ! given a gathered array index, this gives the index within the global array
  CONTAINS
    PROCEDURE, PRIVATE :: construct => reorderData_construct    ! Collective across work AND restart. This IS a two step process utilizing the two methods below.
    PROCEDURE, PRIVATE :: constructCompute => reorderData_constructCompute  ! The reorder DATA IS constructed on the work processes AND THEN transferred to the restart processes.
    PROCEDURE, PRIVATE :: transferToRestart => reorderData_transferToRestart    ! Constructs the reorder DATA on the restart processes.

    PROCEDURE, PRIVATE :: packLevel => reorderData_packLevel    ! pack the contents of a single level/variable into our memory window
    PROCEDURE, PRIVATE :: reorderLevelFromPe => reorderData_reorderLevelFromPe  ! sort the DATA of a single level from a single PE into a global array.

    PROCEDURE, PRIVATE :: destruct => reorderData_destruct
  END TYPE t_reorder_data

  ! this combines all the DATA that's relevant for transfering the payload DATA from the worker PEs to the restart PEs
  TYPE t_restart_comm_data
    ! DATA for remote memory access
    INTEGER :: mpiWindow
    REAL(dp), POINTER :: windowPtr(:)

    ! reorder data
    TYPE(t_reorder_data) :: cells
    TYPE(t_reorder_data) :: edges
    TYPE(t_reorder_data) :: verts
  CONTAINS
    PROCEDURE, PRIVATE :: construct => restartCommData_construct
    PROCEDURE, PRIVATE :: postData => restartCommData_postData  ! called by the compute processes to WRITE their DATA to their memory window
    PROCEDURE, PRIVATE :: collectData => restartCommData_collectData    ! called by the restart processes to fetch the DATA from the compute processes
    PROCEDURE, PRIVATE :: destruct => restartCommData_destruct
  END TYPE t_restart_comm_data

  ! combine the DATA that describes a patch for restart purposes with the infos required for the asynchronous fetching of the DATA from the compute PEs
  TYPE t_patch_data
    TYPE(t_restart_patch_description) :: description
    TYPE(t_restart_comm_data) :: commData
    TYPE(t_restart_file) :: restart_file
  END TYPE t_patch_data

  ! This IS the actual INTERFACE to the restart writing code (apart from the restart_main_proc PROCEDURE). Its USE IS as follows:
  !
  ! First, AND ONLY once during a run, a t_restart_descriptor IS constructed.
  !
  ! Then, for each restart that IS to be written, the updatePatch() method IS used to set the current time dependend information for each patch.
  ! Once all patches are updated, a single CALL to writeRestart() triggers the actual restart writing.
  ! The updatePatch() - writeRestart() sequence can be repeated ANY number of times.
  !
  ! Finally, destruct() must be called to signal the restart PEs to finish their work, AND to wait for them to stop.
  TYPE t_restart_descriptor
    TYPE(t_patch_data), ALLOCATABLE :: patch_data(:)
  CONTAINS
    PROCEDURE :: construct => restartDescriptor_construct
    PROCEDURE :: updatePatch => restartDescriptor_updatePatch
    PROCEDURE :: writeRestart => restartDescriptor_writeRestart
    PROCEDURE :: destruct => restartDescriptor_destruct

    ! methods called ONLY by the restart processes
    PROCEDURE, PRIVATE :: restartWriteAsyncRestart => restartDescriptor_restartWriteAsyncRestart
  END TYPE t_restart_descriptor

CONTAINS

#ifndef NOMPI
  ! Broadcast root for intercommunicator broadcasts from compute PEs to restart PEs using p_comm_work_2_restart.
  INTEGER FUNCTION bcastRoot() RESULT(RESULT)
    IF(my_process_is_restart()) THEN
        ! root is proc 0 on the compute PEs
        RESULT = 0
    ELSE
        ! Special root setting for intercommunicators:
        ! The PE really sending must use MPI_ROOT, the others MPI_PROC_NULL.
        IF(p_pe_work == 0) THEN
            RESULT = MPI_ROOT
        ELSE
            RESULT = MPI_PROC_NULL
        END IF
    END IF
  END FUNCTION bcastRoot
#endif

  !------------------------------------------------------------------------------------------------
  !
  ! public routines
  !
  !------------------------------------------------------------------------------------------------
  !
  !> Prepare the asynchronous restart (collective call).
  !
  SUBROUTINE restartDescriptor_construct(me)
    CLASS(t_restart_descriptor), INTENT(INOUT) :: me

    INTEGER :: jg, error
    CHARACTER(LEN=*), PARAMETER :: routine = modname//':restartDescriptor_construct'

#ifdef NOMPI
    CALL finish(routine, ASYNC_RESTART_REQ_MPI)
#else

    IF(.NOT. (my_process_is_work() .OR. my_process_is_restart())) RETURN

    ! TRANSFER some global DATA to the restart processes
    CALL p_bcast(n_dom, bcastRoot(), p_comm_work_2_restart)
    CALL bcastRestartVarlists(bcastRoot(), p_comm_work_2_restart)
    CALL RestartNamelist_bcast(bcastRoot(), p_comm_work_2_restart)

    ! allocate patch data structure
    ALLOCATE(me%patch_data(n_dom), STAT = error)
    IF(error /= SUCCESS) CALL finish(routine, ALLOCATE_FAILED)

    ! initialize the patch data structures
    DO jg = 1, n_dom
        CALL create_patch_description(me%patch_data(jg)%description, jg)
        CALL me%patch_data(jg)%restart_file%construct(jg)
        CALL me%patch_data(jg)%commData%construct(jg, me%patch_data(jg)%restart_file%var_data)
    END DO
#endif
  END SUBROUTINE restartDescriptor_construct

  !------------------------------------------------------------------------------------------------
  !
  !> Set patch-dependent dynamic data for asynchronous restart.
  !
  SUBROUTINE restartDescriptor_updatePatch(me, patch, opt_pvct, opt_t_elapsed_phy, opt_lcall_phy, opt_sim_time, &
                                          &opt_ndyn_substeps, opt_jstep_adv_marchuk_order, opt_depth, opt_depth_lnd, &
                                          &opt_nlev_snow, opt_nice_class, opt_ndom)
    CLASS(t_restart_descriptor), INTENT(INOUT) :: me
    TYPE(t_patch), INTENT(IN) :: patch
    INTEGER, INTENT(IN), OPTIONAL :: opt_depth, opt_depth_lnd, opt_ndyn_substeps, opt_jstep_adv_marchuk_order, &
                                   & opt_nlev_snow, opt_nice_class, opt_ndom
    REAL(wp), INTENT(IN), OPTIONAL :: opt_sim_time, opt_pvct(:), opt_t_elapsed_phy(:)
    LOGICAL, INTENT(IN), OPTIONAL :: opt_lcall_phy(:)

    INTEGER :: jg
    CHARACTER(LEN = *), PARAMETER :: routine = modname//":restartDescriptor_updatePatch"

#ifdef NOMPI
    CALL finish(routine, ASYNC_RESTART_REQ_MPI)
#else
    IF(my_process_is_work()) THEN
        jg = patch%id
        IF(jg < 1 .OR. jg > SIZE(me%patch_data)) CALL finish(routine, "assertion failed: patch id IS OUT of range")
        IF(me%patch_data(jg)%description%id /= jg) CALL finish(routine, "assertion failed: patch id doesn't match its array index")
        CALL me%patch_data(jg)%description%update(patch, opt_pvct, opt_t_elapsed_phy, opt_lcall_phy, opt_sim_time, &
                                                 &opt_ndyn_substeps, opt_jstep_adv_marchuk_order, opt_depth, opt_depth_lnd, &
                                                 &opt_nlev_snow, opt_nice_class, opt_ndom)
    END IF
#endif

  END SUBROUTINE restartDescriptor_updatePatch

  !------------------------------------------------------------------------------------------------
  !
  !> Writes all restart data into one or more files (one file per patch, collective across work processes).
  !
  SUBROUTINE restartDescriptor_writeRestart(me, datetime, jstep, opt_output_jfile)
    CLASS(t_restart_descriptor), INTENT(INOUT) :: me
    TYPE(t_datetime), INTENT(IN) :: datetime
    INTEGER, INTENT(IN) :: jstep
    INTEGER, INTENT(IN), OPTIONAL :: opt_output_jfile(:)

    INTEGER :: idx
    TYPE(t_restart_args) :: restart_args
    CHARACTER(LEN=*), PARAMETER :: routine = modname//':restartDescriptor_writeRestart'

#ifdef NOMPI
    CALL finish (routine, ASYNC_RESTART_REQ_MPI)
#else

#ifdef DEBUG
    WRITE (nerr,FORMAT_VALS3)routine,' p_pe=',p_pe
#endif

    ! check kind of process
    IF(.NOT. my_process_is_work()) RETURN

    CALL compute_wait_for_restart()
    CALL restart_args%construct(datetime, jstep, opt_output_jfile)
    CALL compute_start_restart(restart_args, me%patch_data)
    CALL restart_args%destruct()

    ! do the restart output
    DO idx = 1, SIZE(me%patch_data)
      ! collective call to write the restart variables
      IF(me%patch_data(idx)%description%l_dom_active) CALL compute_write_var_list(me%patch_data(idx))
    END DO
#endif
  END SUBROUTINE restartDescriptor_writeRestart

#ifndef NOMPI
  SUBROUTINE restart_write_patch(restart_args, patch_data, restartAttributes)
    TYPE(t_restart_args), INTENT(IN) :: restart_args
    TYPE(t_patch_data), INTENT(INOUT) :: patch_data
    TYPE(t_RestartAttributeList), INTENT(INOUT) :: restartAttributes

    ! check if the patch is actice
    IF(.NOT. patch_data%description%l_dom_active) RETURN

    ! consider the right restart process
    IF(p_pe == patch_data%description%restart_proc_id) THEN
        ! set global restart attributes/lists
        CALL patch_data%description%defineVGrids()

#ifdef DEBUG
        CALL print_restart_arguments()
        CALL restartAttributes%printAttributes()
        CALL print_restart_name_lists()
#endif
        CALL open_restart_file(patch_data, restart_args, restartAttributes)

        ! collective call to write the restart variables
        CALL restart_write_var_list(patch_data)
        IF(patch_data%description%l_opt_ndom) THEN
            CALL create_restart_file_link(TRIM(patch_data%restart_file%filename), TRIM(patch_data%restart_file%model_type), &
                                         &patch_data%description%restart_proc_id - p_restart_pe0, patch_data%description%id, &
                                         &opt_ndom = patch_data%description%opt_ndom)
        ELSE
            CALL create_restart_file_link(TRIM(patch_data%restart_file%filename), TRIM(patch_data%restart_file%model_type), &
                                         &patch_data%description%restart_proc_id - p_restart_pe0, patch_data%description%id)
        END IF
        CALL close_restart_file(patch_data%restart_file)
    ENDIF
  END SUBROUTINE restart_write_patch
#endif

  !> Writes all restart data into one or more files (one file per patch, collective across restart processes).
  SUBROUTINE restartDescriptor_restartWriteAsyncRestart(me, restart_args)
    CLASS(t_restart_descriptor), INTENT(INOUT) :: me
    TYPE(t_restart_args), INTENT(IN) :: restart_args

    INTEGER :: idx
    TYPE(t_RestartAttributeList), POINTER :: restartAttributes
    CHARACTER(LEN=*), PARAMETER :: routine = modname//':restartDescriptor_restartWriteAsyncRestart'

#ifdef NOMPI
    CALL finish (routine, ASYNC_RESTART_REQ_MPI)
#else

#ifdef DEBUG
    WRITE (nerr,FORMAT_VALS3)routine,' p_pe=',p_pe
#endif

    ! check kind of process
    IF(.NOT. my_process_is_restart()) CALL finish(routine, "assertion failed: "//routine//"() called by non-restart process")

    restartAttributes => RestartAttributeList_make()
    CALL set_restart_attributes(restartAttributes, restart_args, me%patch_data)

    ! do the restart output
    DO idx = 1, SIZE(me%patch_data)
        CALL restart_write_patch(restart_args, me%patch_data(idx), restartAttributes)
    END DO

    CALL restartAttributes%destruct()
    DEALLOCATE(restartAttributes)
#endif
  END SUBROUTINE restartDescriptor_restartWriteAsyncRestart

  !------------------------------------------------------------------------------------------------
  !
  !> Closes asynchronous restart (collective call).
  !
  SUBROUTINE restartDescriptor_destruct(me)
    CLASS(t_restart_descriptor), INTENT(INOUT) :: me

    INTEGER :: i
    CHARACTER(LEN=*), PARAMETER :: routine = modname//':restartDescriptor_destruct'

#ifdef NOMPI
    CALL finish(routine, ASYNC_RESTART_REQ_MPI)
#else
    ! check kind of process
    IF (.NOT. my_process_is_work() .AND. .NOT. my_process_is_restart()) RETURN

    IF (my_process_is_work()) THEN

      CALL compute_wait_for_restart()
      CALL compute_shutdown_restart

    ENDIF

    DO i = 1, SIZE(me%patch_data, 1)
        CALL me%patch_data(i)%commData%destruct()
    END DO
    DEALLOCATE(me%patch_data)
#endif

  END SUBROUTINE restartDescriptor_destruct

  !-------------------------------------------------------------------------------------------------
  !>
  !! Main routine for restart PEs.
  !! Please note that this routine never returns.
  SUBROUTINE restart_main_proc
    LOGICAL :: done
    TYPE(t_restart_descriptor) :: restartDescriptor
    TYPE(t_restart_args) :: restart_args
    CHARACTER(LEN=*), PARAMETER :: routine = modname//':restart_main_proc'

#ifdef NOMPI
    CALL finish(routine, ASYNC_RESTART_REQ_MPI)
#else

#ifdef DEBUG
    WRITE (nerr,FORMAT_VALS3)routine,' p_pe=',p_pe
#endif

    ! check kind of process
    IF (.NOT. my_process_is_restart()) RETURN

    ! prepare restart (collective call)
    CALL restartDescriptor%construct()

    ! tell the compute PEs that we are ready to work
    CALL restart_send_ready

    ! enter restart loop
    done = .FALSE.
    DO
      ! wait for a message from the compute PEs to start
      CALL restart_wait_for_start(restart_args, restartDescriptor, done)   ! this constructs the restart_args

      IF(done) EXIT ! leave loop, we are done

      ! read and write restart variable lists (collective call)
      CALL restartDescriptor%restartWriteAsyncRestart(restart_args)
      CALL restart_args%destruct()

      ! inform compute PEs that the restart is done
      CALL restart_send_ready
    ENDDO

    ! finalization sequence (collective call)
    CALL restartDescriptor%destruct()

    ! shut down MPI
    CALL stop_mpi

    STOP

#endif

  END SUBROUTINE restart_main_proc

!---------------------------------------------------------------------------------------------------
!
! All other routines are needed for MPI compilation.
!
#ifndef NOMPI
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

  SUBROUTINE restartMetadataPacker(operation, restart_args, patch_data, message)
    INTEGER, VALUE :: operation
    TYPE(t_restart_args), INTENT(INOUT) :: restart_args
    TYPE(t_patch_data), INTENT(INOUT) :: patch_data(:)
    TYPE(t_PackedMessage), INTENT(INOUT) :: message

    INTEGER :: i, calday, patchId
    CHARACTER(LEN = *), PARAMETER :: routine = modname//":restartMetadataPacker"

    ! (un)pack patch independent arguments
    CALL restart_args%packer(operation, message)

    ! (un)pack the patch descriptions
    DO i = 1, SIZE(patch_data)
        CALL patch_data(i)%description%packer(operation, message)
    END DO
  END SUBROUTINE restartMetadataPacker

  !-------------------------------------------------------------------------------------------------
  !>
  !! restart_wait_for_start: Wait for a message from compute PEs that we should start restart or finish.
  !! The counterpart on the compute side is compute_start_restart/compute_shutdown_restart.
  !
  SUBROUTINE restart_wait_for_start(restart_args, restartDescriptor, done)
    TYPE(t_restart_args), INTENT(INOUT) :: restart_args
    CLASS(t_restart_descriptor), INTENT(INOUT) :: restartDescriptor
    LOGICAL, INTENT(OUT)           :: done ! flag if we should shut down

    INTEGER :: iheader
    TYPE(t_PackedMessage) :: message
    CHARACTER(LEN=*), PARAMETER :: routine = modname//':restart_wait_for_start'

#ifdef DEBUG
    WRITE (nerr,FORMAT_VALS3)routine,' is called, p_pe=',p_pe
#endif

    CALL message%construct() ! create message array

    ! receive message that we may start restart (or should finish)
    IF(p_pe_work == 0) CALL message%recv(p_work_pe0, 0, process_mpi_all_comm)
    CALL message%bcast(0, p_comm_work)

    ! unpack AND interpret the message
    CALL message%unpack(iheader)
    SELECT CASE(iheader)
      CASE(MSG_RESTART_START)
        CALL restartMetadataPacker(kUnpackOp, restart_args, restartDescriptor%patch_data, message)
        done = .FALSE.

      CASE(MSG_RESTART_SHUTDOWN)
        done = .TRUE.

      CASE DEFAULT
        ! anything else is an error
        CALL finish(routine,'restart PE: Got illegal restart tag')

    END SELECT

    CALL message%destruct() ! cleanup
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

    ! first compute PE receives message from restart leader
    IF(p_pe_work == 0) THEN
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

  END SUBROUTINE compute_wait_for_restart

  ! the patch's master process sends its patch description to the work master
  ! for all other processes, this IS a noop
  SUBROUTINE sendDescriptionToMaster(patchDescription)
    TYPE(t_restart_patch_description), INTENT(INOUT) :: patchDescription

    TYPE(t_PackedMessage) :: message

    IF(patchDescription%work_pe0_id == 0) RETURN   ! nothing to communicate IF PE0 IS already the subset master

    CALL message%construct()

    IF (p_pe_work == 0) THEN
        ! recieve the package for this patch
        CALL message%recv(patchDescription%work_pe0_id, 0, process_mpi_all_comm)
        CALL patchDescription%packer(kUnpackOp, message)
    ELSE IF (p_pe_work == patchDescription%work_pe0_id) THEN
        ! send the time dependent DATA to process 0
        CALL patchDescription%packer(kPackOp, message)
        CALL message%send(0, 0, process_mpi_all_comm)
    END IF

    CALL message%destruct()
  END SUBROUTINE sendDescriptionToMaster
  !-------------------------------------------------------------------------------------------------
  !>
  !! compute_start_restart: Send a message to restart PEs that they should start restart.
  !! The counterpart on the restart side is restart_wait_for_start.
  !
  SUBROUTINE compute_start_restart(restart_args, patch_data)
    TYPE(t_restart_args), INTENT(INOUT) :: restart_args
    TYPE(t_patch_data), INTENT(INOUT) :: patch_data(:)

    INTEGER :: i, trash
    TYPE(t_PackedMessage) :: message
    CHARACTER(LEN=*), PARAMETER :: routine = modname//':compute_start_restart'

#ifdef DEBUG
    WRITE (nerr,FORMAT_VALS5)routine,' p_pe=',p_pe, ' call p_barrier with communicator=',p_comm_work
#endif

    ! make sure all are here
    CALL p_barrier(comm=p_comm_work)

    ! if processor splitting is applied, the time-dependent data need to be transferred
    ! from the subset master PE to PE0, from where they are communicated to the output PE(s)
    CALL message%construct()
    DO i = 1, SIZE(patch_data)
      CALL patch_data(i)%description%setTimeLevels() ! copy the global variables (nold, ...) to the patch description
      CALL sendDescriptionToMaster(patch_data(i)%description)
    END DO

    CALL message%reset()
    IF(p_pe_work == 0) THEN
      ! send the DATA to the restart master
      CALL message%pack(MSG_RESTART_START)  ! set command id
      CALL restartMetadataPacker(kPackOp, restart_args, patch_data, message)    ! all the other DATA
      CALL message%send(p_restart_pe0, 0, process_mpi_all_comm)
    ENDIF
    ! broadcast a copy among the compute processes, so that all processes have a consistent view of the restart patch descriptions
    CALL message%bcast(0, p_comm_work)
    CALL message%unpack(trash)  ! ignore the command id
    CALL restartMetadataPacker(kUnpackOp, restart_args, patch_data, message)

    CALL message%destruct()
  END SUBROUTINE compute_start_restart

  !-------------------------------------------------------------------------------------------------
  !>
  !! compute_shutdown_restart: Send a message to restart PEs that they should shut down.
  !! The counterpart on the restart side is restart_wait_for_start.
  !
  SUBROUTINE compute_shutdown_restart
    TYPE(t_PackedMessage) :: message
    CHARACTER(LEN=*), PARAMETER :: routine = modname//':compute_shutdown_restart'

#ifdef DEBUG
    WRITE (nerr,FORMAT_VALS5)routine,' p_pe=',p_pe, &
      & ' call p_barrier with communicator=',p_comm_work
#endif

    ! make sure all are here
    CALL p_barrier(comm=p_comm_work)

    IF(p_pe_work == 0) THEN
      CALL message%construct()  ! create message array

      ! sent the shutdown message
      CALL message%pack(MSG_RESTART_SHUTDOWN)
      CALL message%send(p_restart_pe0, 0, process_mpi_all_comm)

      CALL message%destruct()   ! cleanup
    ENDIF

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

    CHARACTER (LEN=MAX_ERROR_LENGTH):: error_message

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
  SUBROUTINE print_restart_arguments(restart_args, patch_data)
    TYPE(t_restart_args), INTENT(IN) :: restart_args
    TYPE(t_patch_data), INTENT(IN) :: patch_data(:)
    CHARACTER(LEN=*), PARAMETER   :: routine = modname//':print_restart_arguments'

    WRITE (nerr,FORMAT_VALS3)routine,' is called for p_pe=',p_pe

    CALL restart_args%print(routine//": ")

    ! patch informations
    PRINT *,routine, ' size of patches=', SIZE(patch_data)
    PRINT *,routine, ' SIZE(patch_data(1)%description%opt_pvct) = ', SIZE(patch_data(1)%description%opt_pvct)
    PRINT *,routine, ' SIZE(patch_data(1)%description%opt_lcall_phy) =', SIZE(patch_data(1)%description%opt_lcall_phy)
    PRINT *,routine, ' SIZE(patch_data(1)%description%opt_t_elapsed_phy) =', SIZE(patch_data(1)%description%opt_t_elapsed_phy)
  END SUBROUTINE print_restart_arguments

  !------------------------------------------------------------------------------------------------
  !
  ! Gets the number of  restart variables for the given logical patch ident.
  !
  SUBROUTINE get_var_list_number(all_fld_cnt, patch_id)

    INTEGER,    INTENT(INOUT)     :: all_fld_cnt
    INTEGER,    INTENT(IN)        :: patch_id

    INTEGER                       :: i, fld_cnt, list_cnt
    TYPE(t_list_element), POINTER :: element

#ifdef DEBUG
    TYPE(t_list_element), POINTER :: element_list
    CHARACTER(LEN=*), PARAMETER   :: routine = modname//':get_var_list_number'

    WRITE (nerr,FORMAT_VALS5)routine,' p_pe=',p_pe,' patch_id=',patch_id
#endif

    all_fld_cnt = 0
    list_cnt = 0

    DO i = 1, nvar_lists
      ! skip, if var_list is not required for restart
      IF (.NOT. var_lists(i)%p%lrestart) CYCLE

      ! check the given logical patch id
      IF (var_lists(i)%p%patch_id /= patch_id) CYCLE

      ! check, if the list has valid restart fields
      fld_cnt = 0
      element => var_lists(i)%p%first_list_element
      DO
        IF(.NOT. ASSOCIATED(element)) EXIT
        IF (element%field%info%lrestart) THEN
          fld_cnt = fld_cnt + 1
        ENDIF
        element => element%next_list_element
      ENDDO
      all_fld_cnt = all_fld_cnt + fld_cnt

#ifdef DEBUG
      IF (my_process_is_restart()) THEN
        IF (fld_cnt > 0) THEN
          list_cnt = list_cnt + 1
          WRITE(nerr,'(i4,3a,i4,2a,i3)') &
            & list_cnt,'. restart var_list ',TRIM(var_lists(i)%p%name), &
            & '(',fld_cnt,')', &
            & ' Patch: ',var_lists(i)%p%patch_id
          element_list => var_lists(i)%p%first_list_element
          DO
            IF(.NOT. ASSOCIATED(element_list)) EXIT
            IF (element_list%field%info%lrestart) THEN
              WRITE (nerr,'(4a)') &
                  &     '    ',TRIM(element_list%field%info%name), &
                  &       '  ',TRIM(element_list%field%info%cf%long_name)
            ENDIF
            element_list => element_list%next_list_element
          ENDDO
        ENDIF
      ENDIF
#endif

    ENDDO

  END SUBROUTINE get_var_list_number

  !------------------------------------------------------------------------------------------------
  !
  ! Common helper routines to preparing the restart.
  !
  !-------------------------------------------------------------------------------------------------
  SUBROUTINE restartVarlistPacker(operation, message)
    INTEGER, VALUE :: operation
    TYPE(t_PackedMessage), INTENT(INOUT) :: message

    INTEGER :: info_size, iv, nv, nelems, patch_id, restart_type, vlevel_type, n, ierrstat
    INTEGER, ALLOCATABLE            :: info_storage(:)
    TYPE(t_list_element), POINTER   :: element, newElement
    TYPE(t_var_metadata)            :: info
    TYPE(t_var_list)                :: p_var_list
    CHARACTER(LEN=MAX_NAME_LENGTH)  :: var_list_name
    CHARACTER(LEN=32)               :: model_type
    LOGICAL                         :: lrestart

    CHARACTER(LEN=*), PARAMETER :: routine = modname//':restartVarlistPacker'

    ! delete old var lists
    IF(operation == kUnpackOp) CALL delete_var_lists

    ! get the size - in default INTEGER words - which is needed to
    ! hold the contents of TYPE(t_var_metadata)
    info_size = SIZE(TRANSFER(info, (/ 0 /)))
    ALLOCATE(info_storage(info_size), STAT=ierrstat)
    IF(ierrstat /= SUCCESS) CALL finish(routine, ALLOCATE_FAILED)

    ! get the number of var lists
    nv = nvar_lists
    CALL message%packer(operation, nv)

    ! for each var list, get its components
    DO iv = 1, nv
        IF(operation == kPackOp) THEN
            ! copy the values needed for the new_var_list() CALL to local variables
            lrestart = var_lists(iv)%p%lrestart
            var_list_name = var_lists(iv)%p%name
            model_type = var_lists(iv)%p%model_type
            patch_id = var_lists(iv)%p%patch_id
            restart_type = var_lists(iv)%p%restart_type
            vlevel_type = var_lists(iv)%p%vlevel_type

            ! count the number of variable restart entries
            element => var_lists(iv)%p%first_list_element
            nelems = 0
            DO
                IF(.NOT.ASSOCIATED(element)) EXIT
                IF(element%field%info%lrestart) nelems = nelems+1
                element => element%next_list_element
            END DO
        END IF
        CALL message%packer(operation, lrestart)
        CALL message%packer(operation, var_list_name)
        CALL message%packer(operation, model_type)
        CALL message%packer(operation, patch_id)
        CALL message%packer(operation, restart_type)
        CALL message%packer(operation, vlevel_type)
        CALL message%packer(operation, nelems)

        IF(.NOT. lrestart) CYCLE  ! transfer only a restart var_list
        IF(nelems == 0) CYCLE ! check if there are valid restart fields

        IF(operation == kPackOp) THEN
            element => var_lists(iv)%p%first_list_element
            DO
                IF(.NOT. ASSOCIATED(element)) EXIT
                IF(element%field%info%lrestart) THEN
                    info_storage = TRANSFER(element%field%info, (/ 0 /))
                    CALL message%packer(operation, info_storage)
                END IF
                element => element%next_list_element
            END DO
        END IF

        IF(operation == kUnpackOp) THEN
            ! create var list
            CALL new_var_list(p_var_list, var_list_name, patch_id=patch_id, restart_type=restart_type, vlevel_type=vlevel_type, &
                             &lrestart=.TRUE.)
            p_var_list%p%model_type = TRIM(model_type)

            ! insert elements into var list
            DO n = 1, nelems
                ! ALLOCATE a new element
                ALLOCATE(newElement, STAT=ierrstat)
                IF(ierrstat /= SUCCESS) CALL finish(routine, ALLOCATE_FAILED)
                IF(n == 1) THEN   ! the first element pointer needs to be stored IN a different variable than the later pointers (there are no double pointers IN FORTRAN...)
                    p_var_list%p%first_list_element => newElement
                ELSE
                    element%next_list_element => newElement
                END IF
                element => newElement
                element%next_list_element => NULL()

                ! these pointers don't make sense on the restart PEs, NULLIFY them
                NULLIFY(element%field%r_ptr, element%field%i_ptr, element%field%l_ptr)
                element%field%var_base_size = 0 ! Unknown here

                ! set info structure from binary representation in info_storage
                CALL message%packer(operation, info_storage)
                element%field%info = TRANSFER(info_storage, info)
            END DO
        END IF

        DEALLOCATE(info_storage)
    END DO
  END SUBROUTINE restartVarlistPacker

  !
  ! Transfers the restart var lists from the worker to the restart PEs.
  !
  SUBROUTINE bcastRestartVarlists(root, communicator)
    INTEGER, VALUE :: root, communicator

    TYPE(t_PackedMessage) :: message
    LOGICAL :: lIsSender, lIsReceiver
    CHARACTER(LEN=*), PARAMETER :: routine = modname//':bcastRestartVarlists'

    CALL p_get_bcast_role(root, communicator, lIsSender, lIsReceiver)

    CALL message%construct()
    IF(lIsSender) CALL restartVarlistPacker(kPackOp, message)
    CALL message%bcast(root, communicator)
    IF(lIsReceiver) CALL restartVarlistPacker(kUnpackOp, message)
    CALL message%destruct()
  END SUBROUTINE bcastRestartVarlists

  ! collective across restart AND worker PEs
  SUBROUTINE create_patch_description(description, domain)
    TYPE(t_restart_patch_description), INTENT(INOUT) :: description
    INTEGER, VALUE :: domain

    TYPE(t_PackedMessage) :: message

    ! initialize on work PEs
    CALL message%construct()
    IF(my_process_is_work()) THEN
        CALL description%init(p_patch(domain))
        CALL description%packer(kPackOp, message)
    END IF

    ! transfer data to restart PEs
    CALL message%bcast(bcastRoot(), p_comm_work_2_restart)
    CALL description%packer(kUnpackOp, message)

    ! initialize the fields that we DO NOT communicate from the worker PEs to the restart PEs
    description%restart_proc_id = MOD(description%id-1, process_mpi_restart_size) + p_restart_pe0
    description%v_grid_count = 0

    CALL message%destruct() ! cleanup
  END SUBROUTINE create_patch_description

  ! collective across restart AND worker PEs
  ! returns the memory window offset for the next t_restart_comm_data object
  SUBROUTINE restartCommData_construct(me, jg, var_data)
    CLASS(t_restart_comm_data), TARGET, INTENT(INOUT) :: me
    INTEGER, VALUE :: jg
    TYPE(t_var_data), POINTER, INTENT(IN) :: var_data(:)

    INTEGER :: error, iv, nlevs
    TYPE(t_grid_domain_decomp_info) :: dummyInfo
    INTEGER(KIND = MPI_ADDRESS_KIND) :: memWindowSize
    TYPE(t_reorder_data), POINTER :: reorderData
    CHARACTER(LEN = *), PARAMETER :: routine = modname//":restartCommData_construct"

    ! initialize the reorder info
    IF(my_process_is_work()) THEN
        CALL me%cells%construct(p_patch(jg)%n_patch_cells_g, p_patch(jg)%n_patch_cells, p_patch(jg)%cells%decomp_info)
        CALL me%edges%construct(p_patch(jg)%n_patch_edges_g, p_patch(jg)%n_patch_edges, p_patch(jg)%edges%decomp_info)
        CALL me%verts%construct(p_patch(jg)%n_patch_verts_g, p_patch(jg)%n_patch_verts, p_patch(jg)%verts%decomp_info)
    ELSE
        !must not access p_patch on the restart processes, nevertheless the restart processes need to take part in the collective calls
        CALL me%cells%construct(0, 0, dummyInfo)
        CALL me%edges%construct(0, 0, dummyInfo)
        CALL me%verts%construct(0, 0, dummyInfo)
    END IF

    ! initialize the memory window offsets
    memWindowSize = 0
    IF(ASSOCIATED(var_data)) THEN
        ! compute the SIZE of the memory window
        DO iv = 1, SIZE(var_data)
            nlevs = 1
            IF(var_data(iv)%info%ndims /= 2) nlevs = var_data(iv)%info%used_dimensions(2)

            reorderData => get_reorder_ptr(me, var_data(iv)%info, routine)
            memWindowSize = memWindowSize + INT(nlevs*reorderData%n_own, i8)
        END DO

        ! actually open the memory window
        CALL openMpiWindow(memWindowSize, p_comm_work_restart, me%windowPtr, me%mpiWindow)
    END IF
  END SUBROUTINE restartCommData_construct

  SUBROUTINE restartCommData_postData(me, varInfo, dataPointers, offset)
    CLASS(t_restart_comm_data), TARGET, INTENT(INOUT) :: me
    TYPE(t_var_metadata), INTENT(IN) :: varInfo
    TYPE(t_ptr_2d), ALLOCATABLE, INTENT(IN) :: dataPointers(:)
    !TODO[NH]: Replace this by an index within the t_restart_comm_data structure.
    INTEGER(i8), INTENT(INOUT) :: offset

    TYPE(t_reorder_data), POINTER :: reorderData
    INTEGER :: mpi_error, i
    CHARACTER(LEN = *), PARAMETER :: routine = modname//":restartCommData_postData"

    ! get pointer to reorder data
    reorderData => get_reorder_ptr(me, varInfo, routine)

#ifdef DEBUG
    WRITE(nerr,FORMAT_VALS7I) routine,' p_pe=',p_pe,' compute pe writes field=', TRIM(varInfo%name),' data=', &
                            & SIZE(dataPointers)*reorderData%n_own
#endif

    ! lock own window before writing to it
    CALL MPI_Win_lock(MPI_LOCK_EXCLUSIVE, p_pe_work, MPI_MODE_NOCHECK, me%mpiWindow, mpi_error)
    CALL check_mpi_error(routine, 'MPI_Win_lock', mpi_error, .TRUE.)

    ! just copy the OWN data points to the memory window
    DO i = 1, SIZE(dataPointers)
        CALL reorderData%packLevel(dataPointers(i)%p, me%windowPtr, offset)
    END DO

    ! unlock RMA window
    CALL MPI_Win_unlock(p_pe_work, me%mpiWindow, mpi_error)
    CALL check_mpi_error(routine, 'MPI_Win_unlock', mpi_error, .TRUE.)
  END SUBROUTINE restartCommData_postData

  SUBROUTINE restartCommData_collectData(me, varInfo, levelCount, dest, offsets, elapsedTime, bytesFetched)
    CLASS(t_restart_comm_data), TARGET, INTENT(INOUT) :: me
    TYPE(t_var_metadata), POINTER :: varInfo
    INTEGER, VALUE :: levelCount
    REAL(dp), INTENT(INOUT) :: dest(:,:)
    !TODO[NH]: Replace this by an index within the t_restart_comm_data structure.
    INTEGER(KIND=MPI_ADDRESS_KIND), INTENT(INOUT) :: offsets(0:num_work_procs-1)   ! This remembers the current offsets within all the different processes RMA windows. Pass a zero initialized array to the first collectData() CALL, THEN keep passing the array without external modifications.
    REAL(dp), INTENT(INOUT) :: elapsedTime
    INTEGER(i8), INTENT(INOUT) :: bytesFetched

    INTEGER :: sourceProc, doubleCount, reorderOffset, curLevel, mpi_error
    TYPE(t_reorder_data), POINTER :: reorderData
    REAL(dp), POINTER, SAVE :: buffer(:) => NULL()
    CHARACTER(LEN = *), PARAMETER :: routine = modname//":restartCommData_collectData"

    reorderData => get_reorder_ptr(me, varInfo, routine)

    ! retrieve part of variable from every worker PE using MPI_Get
    DO sourceProc = 0, num_work_procs-1
        IF(reorderData%pe_own(sourceProc) == 0) CYCLE
        CALL ensureSize(buffer, reorderData%pe_own(sourceProc))

        ! number of words to transfer
        doubleCount = reorderData%pe_own(sourceProc) * levelCount
        elapsedTime = elapsedTime -  p_mpi_wtime()
        CALL MPI_Win_lock(MPI_LOCK_SHARED, sourceProc, MPI_MODE_NOCHECK, me%mpiWindow, mpi_error)
        CALL check_mpi_error(routine, 'MPI_Win_lock', mpi_error, .TRUE.)

        CALL MPI_Get(buffer(1), doubleCount, p_real_dp, sourceProc, offsets(sourceProc), &
                              & doubleCount, p_real_dp, me%mpiWindow, mpi_error)
        CALL check_mpi_error(routine, 'MPI_Get', mpi_error, .TRUE.)

        CALL MPI_Win_unlock(sourceProc, me%mpiWindow, mpi_error)
        CALL check_mpi_error(routine, 'MPI_Win_unlock', mpi_error, .TRUE.)

        elapsedTime  = elapsedTime + p_mpi_wtime()
        bytesFetched = bytesFetched + INT(doubleCount, i8)*8

        ! update the offset in the RMA window on compute PEs
        offsets(sourceProc) = offsets(sourceProc) + INT(doubleCount, i8)

        ! separate the levels received from PE "sourceProc":
        reorderOffset = 1
        DO curLevel = 1, levelCount
            CALL reorderData%reorderLevelFromPe(curLevel, sourceProc, buffer(reorderOffset:), dest)
            reorderOffset = reorderOffset + reorderData%pe_own(sourceProc)
        END DO
    END DO
  END SUBROUTINE restartCommData_collectData

  SUBROUTINE restartCommData_destruct(me)
    CLASS(t_restart_comm_data), INTENT(INOUT) :: me

    CALL closeMpiWindow(me%windowPtr, me%mpiWindow)

    CALL me%cells%destruct()
    CALL me%verts%destruct()
    CALL me%edges%destruct()
  END SUBROUTINE restartCommData_destruct

  !------------------------------------------------------------------------------------------------
  !
  ! Sets the restart file data with the given logical patch ident.
  !
  SUBROUTINE restartFile_construct(me, patch_id)
    CLASS(t_restart_file), INTENT (INOUT) :: me
    INTEGER, VALUE :: patch_id

    INTEGER                               :: ierrstat, i, i2, num_vars
    TYPE (t_list_element), POINTER        :: element
    CHARACTER(LEN=*), PARAMETER           :: routine = modname//':restartFile_construct'

#ifdef DEBUG
    WRITE (nerr,FORMAT_VALS5)routine,' is called for p_pe=',p_pe,' patch_id=',patch_id
#endif

    ! init. main variables
    me%var_data => NULL()
    me%filename = ''

    CALL me%cdiIds%init()

    ! counts number of restart variables for this file (logical patch ident)
    num_vars = 0
    CALL get_var_list_number(num_vars, patch_id)

#ifdef DEBUG
    WRITE (nerr,FORMAT_VALS3)routine,' numvars=',num_vars
#endif

    IF (num_vars <= 0) RETURN
    ! allocate the array of restart variables
    ALLOCATE (me%var_data(num_vars), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, ALLOCATE_FAILED)

    ! fill the array of restart variables
    i2 = 0
    DO i = 1, nvar_lists
      ! skip, if var_list is not required for restart
      IF (.NOT. var_lists(i)%p%lrestart) CYCLE

      ! check the given logical patch id
      IF (var_lists(i)%p%patch_id /= patch_id) CYCLE

      ! check, if the list has valid restart fields
      element => var_lists(i)%p%first_list_element
      DO
        IF(.NOT. ASSOCIATED(element)) EXIT
        IF (element%field%info%lrestart) THEN
          i2 = i2 + 1
          me%var_data(i2)%info = element%field%info
          IF (my_process_is_work() .AND. ASSOCIATED(element%field%r_ptr)) THEN
            me%var_data(i2)%r_ptr => element%field%r_ptr
          ELSE
            me%var_data(i2)%r_ptr => NULL()
          ENDIF
        ENDIF
        element => element%next_list_element
      ENDDO
    ENDDO

  END SUBROUTINE restartFile_construct

  !------------------------------------------------------------------------------------------------
  !
  ! Sets the reorder data for cells/edges/verts.
  !
  SUBROUTINE reorderData_constructCompute(me, n_points_g, n_points, decomp_info)
    CLASS(t_reorder_data), INTENT(INOUT) :: me ! Result: reorder data
    INTEGER, INTENT(IN) :: n_points_g      ! Global number of cells/edges/verts in logical patch
    INTEGER, INTENT(IN) :: n_points        ! Local number of cells/edges/verts in logical patch
    TYPE(t_grid_domain_decomp_info), INTENT(IN) :: decomp_info

    INTEGER :: i, n, il, ib, mpi_error, ierrstat
    LOGICAL, ALLOCATABLE :: owner_mask_1d(:) ! non-blocked owner mask for (logical) patch
    INTEGER, ALLOCATABLE :: glbidx_own(:), glbidx_glb(:), reorder_index_log_dom(:), pe_off(:), reorder_index(:)
    CHARACTER (LEN=MAX_ERROR_LENGTH) :: error_message

    CHARACTER(LEN=*), PARAMETER :: routine = modname//':reorderData_constructCompute'

#ifdef DEBUG
    WRITE (nerr,FORMAT_VALS3)routine,' is called for p_pe=',p_pe
#endif

    ! just for safety
    IF(my_process_is_restart()) CALL finish(routine, NO_COMPUTE_PE)

    ! set the non-blocked patch owner mask
    ALLOCATE(owner_mask_1d(n_points), STAT = ierrstat)
    IF(ierrstat /= SUCCESS) CALL finish(routine, "allocation of owner_mask_1d failed")
    DO i = 1, n_points
      il = idx_no(i)
      ib = blk_no(i)
      owner_mask_1d(i) = decomp_info%owner_mask(il,ib)
    ENDDO

    ! get number of owned cells/edges/verts (without halos)
    me%n_own = COUNT(owner_mask_1d(:))

    ! set index arrays to own cells/edges/verts
    ALLOCATE(me%own_idx(me%n_own), STAT=ierrstat)
    IF(ierrstat /= SUCCESS) CALL finish(routine, "allocation of me%own_idx failed")
    ALLOCATE(me%own_blk(me%n_own), STAT=ierrstat)
    IF(ierrstat /= SUCCESS) CALL finish(routine, "allocation of me%own_blk failed")

    ! global index of my own points
    ALLOCATE(glbidx_own(me%n_own), STAT=ierrstat)
    IF(ierrstat /= SUCCESS) CALL finish(routine, "allocation of glbidx_own failed")

    n = 0
    DO i = 1, n_points
      IF(owner_mask_1d(i)) THEN
        n = n+1
        me%own_idx(n) = idx_no(i)
        me%own_blk(n) = blk_no(i)
        glbidx_own(n)  = decomp_info%glb_index(i)
      ENDIF
    ENDDO
    DEALLOCATE(owner_mask_1d)

    ! gather the number of own points for every PE into me%pe_own
    ALLOCATE(me%pe_own(0:p_n_work-1), STAT=ierrstat)
    IF(ierrstat /= SUCCESS) CALL finish(routine, "allocation of me%pe_own failed")
    ALLOCATE(pe_off(0:p_n_work-1), STAT=ierrstat)
    IF(ierrstat /= SUCCESS) CALL finish(routine, "allocation of pe_off failed")

    CALL MPI_Allgather(me%n_own,  1, p_int, &
                       me%pe_own, 1, p_int, &
                       p_comm_work, mpi_error)
    CALL check_mpi_error(routine, 'MPI_Allgather', mpi_error, .TRUE.)

    ! get offset within result array
    pe_off(0) = 0
    DO i = 1, p_n_work-1
      pe_off(i) = pe_off(i-1) + me%pe_own(i-1)
    ENDDO

    ! get global number of points for current patch
    me%n_glb = SUM(me%pe_own(:))

    ! Get the global index numbers of the data when it is gathered on PE 0
    ! exactly in the same order as it is retrieved later during restart.
    ALLOCATE(glbidx_glb(me%n_glb), STAT=ierrstat)
    IF(ierrstat /= SUCCESS) CALL finish(routine, "allocation of glbidx_glb failed")

    CALL MPI_Allgatherv(glbidx_own, me%n_own, p_int, &
                        glbidx_glb, me%pe_own, pe_off, p_int, &
                        p_comm_work, mpi_error)
    CALL check_mpi_error(routine, 'MPI_Allgatherv', mpi_error, .TRUE.)
    DEALLOCATE(glbidx_own)
    DEALLOCATE(pe_off)

    ! spans the complete logical domain
    ALLOCATE(reorder_index_log_dom(n_points_g), STAT=ierrstat)
    IF(ierrstat /= SUCCESS) CALL finish(routine, "allocation of reorder_index_log_dom failed")
    reorder_index_log_dom(:) = 0

    DO i = 1, me%n_glb
      ! reorder_index_log_dom stores where a global point in logical domain comes from.
      ! It is nonzero only at the patch locations
      reorder_index_log_dom(glbidx_glb(i)) = i
    ENDDO
    DEALLOCATE(glbidx_glb)

    ! remove the zero entries from reorder_index_log_dom to form the reorder_index
    ALLOCATE(reorder_index(me%n_glb), STAT=ierrstat)
    IF(ierrstat /= SUCCESS) CALL finish(routine, "allocation of reorder_index failed")

    n = 0
    DO i = 1, n_points_g
      IF(reorder_index_log_dom(i)>0) THEN
        n = n+1
        reorder_index(n) = reorder_index_log_dom(i)
      ENDIF
    ENDDO
    DEALLOCATE(reorder_index_log_dom)
    IF(n /= me%n_glb) CALL finish(routine, "assertion failed: reorder_index is not surjective")

    ! compute the inverse of the reorder_index
    ALLOCATE(me%inverse_reorder_index(me%n_glb), STAT=ierrstat)
    IF(ierrstat /= SUCCESS) CALL finish(routine, "allocation of me%inverse_reorder_index failed")

    me%inverse_reorder_index(:) = 0
    DO i = 1, me%n_glb
        IF(me%inverse_reorder_index(reorder_index(i)) /= 0) THEN
            CALL finish(routine, "assertion failed: reorder_index is not injective")
        END IF
        me%inverse_reorder_index(reorder_index(i)) = i
    END DO
    DEALLOCATE(reorder_index)
  END SUBROUTINE reorderData_constructCompute

  !------------------------------------------------------------------------------------------------
  !
  ! Transfers reorder data to restart PEs.
  !
  SUBROUTINE reorderData_transferToRestart(me)
    CLASS(t_reorder_data), INTENT(INOUT) :: me
    INTEGER                             :: ierrstat

    CHARACTER(LEN=*), PARAMETER :: routine = modname//':reorderData_transferToRestart'

    ! transfer the global number of points, this is not yet known on restart PEs
    CALL p_bcast(me%n_glb,  bcastRoot(), p_comm_work_2_restart)

    IF(my_process_is_restart()) THEN
      me%n_own = 0  ! on restart PEs: n_own = 0, own_idx and own_blk are not allocated

      ! pe_own must be allocated for num_work_procs, not for p_n_work
      ALLOCATE(me%pe_own(0:num_work_procs-1), STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, ALLOCATE_FAILED)

      ALLOCATE(me%inverse_reorder_index(me%n_glb), STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, ALLOCATE_FAILED)
    ENDIF

    CALL p_bcast(me%pe_own, bcastRoot(), p_comm_work_2_restart)
    CALL p_bcast(me%inverse_reorder_index, bcastRoot(), p_comm_work_2_restart)
  END SUBROUTINE reorderData_transferToRestart

  SUBROUTINE reorderData_construct(me, n_points_g, n_points, decomp_info)
    CLASS(t_reorder_data), INTENT(INOUT) :: me
    INTEGER, INTENT(IN) :: n_points_g      ! Global number of cells/edges/verts in logical patch
    INTEGER, INTENT(IN) :: n_points        ! Local number of cells/edges/verts in logical patch
    TYPE(t_grid_domain_decomp_info), INTENT(IN) :: decomp_info

    IF(my_process_is_work()) THEN
        CALL me%constructCompute(n_points_g, n_points, decomp_info)
    END IF
    CALL me%transferToRestart()
  END SUBROUTINE reorderData_construct

  SUBROUTINE reorderData_packLevel(me, source, dest, offset)
    CLASS(t_reorder_data), INTENT(IN) :: me
    REAL(wp), INTENT(IN) :: source(:,:)
    REAL(dp), INTENT(INOUT) :: dest(:)
    !TODO[NH]: Replace this by an INTENT(IN) level count
    INTEGER(i8), INTENT(INOUT) :: offset

    INTEGER :: i

    DO i = 1, me%n_own
        offset = offset + 1
        dest(offset) = REAL(source(me%own_idx(i), me%own_blk(i)), dp)
    END DO
  END SUBROUTINE reorderData_packLevel

  SUBROUTINE reorderData_reorderLevelFromPe(me, level, pe, dataIn, globalArray)
    CLASS(t_reorder_data), INTENT(IN) :: me
    INTEGER, VALUE :: level, pe
    REAL(dp), INTENT(IN) :: dataIn(:)
    REAL(dp), INTENT(INOUT) :: globalArray(:,:)

    INTEGER :: offset, i

    offset = SUM(me%pe_own(0:pe - 1))
    DO i = 1, me%pe_own(pe)
        globalArray(me%inverse_reorder_index(offset + i), level) = dataIn(i)
    END DO
  END SUBROUTINE reorderData_reorderLevelFromPe

  SUBROUTINE reorderData_destruct(me)
    CLASS(t_reorder_data), INTENT(INOUT) :: me

    IF(ALLOCATED(me%own_idx)) DEALLOCATE(me%own_idx)
    IF(ALLOCATED(me%own_blk)) DEALLOCATE(me%own_blk)
    DEALLOCATE(me%pe_own)
    DEALLOCATE(me%inverse_reorder_index)
  END SUBROUTINE reorderData_destruct

  !------------------------------------------------------------------------------------------------
  !
  ! Opens an MPI memory window for the given amount of double values, returning both the ALLOCATED buffer AND the MPI window handle.
  !
  SUBROUTINE openMpiWindow(doubleCount, communicator, mem_ptr_dp, mpi_win)
    INTEGER(KIND=MPI_ADDRESS_KIND), VALUE :: doubleCount   ! the requested memory window SIZE IN doubles
    INTEGER, VALUE :: communicator
    REAL(dp), POINTER, INTENT(OUT) :: mem_ptr_dp(:) ! this returns a fortran POINTER to the memory buffer
    INTEGER, INTENT(OUT) :: mpi_win ! this returns the handle to the MPI window

    INTEGER :: nbytes_real, mpi_error, rma_cache_hint
    INTEGER (KIND=MPI_ADDRESS_KIND) :: mem_bytes
    TYPE(c_ptr) :: c_mem_ptr
    CHARACTER(LEN=*), PARAMETER :: routine = modname//':openMpiWindow'

#ifdef DEBUG
    WRITE (nerr,FORMAT_VALS3)routine,' is called for p_pe=',p_pe
#endif

    mpi_win = MPI_WIN_NULL

    ! doubleCount is calculated as number of variables above, get number of bytes
    ! get the amount of bytes per REAL*8 variable (as used in MPI communication)
    CALL MPI_Type_extent(p_real_dp, nbytes_real, mpi_error)
    CALL check_mpi_error(routine, 'MPI_Type_extent', mpi_error, .TRUE.)

    ! for the restart PEs the amount of memory needed is 0 - allocate at least 1 word there:
    mem_bytes = MAX(doubleCount, 1_i8)*INT(nbytes_real, i8)

    ! allocate amount of memory needed with MPI_Alloc_mem
    ! 
    ! Depending on wether the Fortran 2003 C interoperability features
    ! are available, one needs to use non-standard language extensions
    ! for calls from Fortran, namely Cray Pointers, since
    ! MPI_Alloc_mem wants a C pointer argument.
    !
    ! see, for example: http://www.lrz.de/services/software/parallel/mpi/onesided/
    ! TYPE(c_ptr) and INTEGER(KIND=MPI_ADDRESS_KIND) do NOT necessarily have the same size!!!
    ! So check if at least c_intptr_t and MPI_ADDRESS_KIND are the same, else we may get
    ! into deep, deep troubles!
    ! There is still a slight probability that TYPE(c_ptr) does not have the size indicated
    ! by c_intptr_t since the standard only requires c_intptr_t is big enough to hold pointers
    ! (so it may be bigger than a pointer), but I hope no vendor screws up its ISO_C_BINDING
    ! in such a way!!!
    ! If c_intptr_t<=0, this type is not defined and we can't do this check, of course.
    IF(c_intptr_t > 0 .AND. c_intptr_t /= MPI_ADDRESS_KIND) &
     & CALL finish(routine,'c_intptr_t /= MPI_ADDRESS_KIND, too dangerous to proceed!')

    CALL MPI_Alloc_mem(mem_bytes, MPI_INFO_NULL, c_mem_ptr, mpi_error)
    CALL check_mpi_error(routine, 'MPI_Alloc_mem', mpi_error, .TRUE.)

    ! The NEC requires a standard INTEGER array as 3rd argument for c_f_pointer,
    ! although it would make more sense to have it of size MPI_ADDRESS_KIND.
    NULLIFY(mem_ptr_dp)
#ifdef __SX__
    CALL C_F_POINTER(c_mem_ptr, mem_ptr_dp, [INT(doubleCount)] )
#else
    CALL C_F_POINTER(c_mem_ptr, mem_ptr_dp, [doubleCount] )
#endif

    rma_cache_hint = MPI_INFO_NULL
#ifdef __xlC__
    ! IBM specific RMA hint, that we don't want window caching
    CALL MPI_Info_create(rma_cache_hint, mpi_error);
    CALL check_mpi_error(routine, 'MPI_Info_create', mpi_error, .TRUE.)
    CALL MPI_Info_set(rma_cache_hint, "IBM_win_cache","0", mpi_error)
    CALL check_mpi_error(routine, 'MPI_Info_set', mpi_error, .TRUE.)
#endif

    ! create memory window for communication
    mem_ptr_dp(:) = 0._dp
    CALL MPI_Win_create(mem_ptr_dp, mem_bytes, nbytes_real, MPI_INFO_NULL, communicator, mpi_win, mpi_error)
    CALL check_mpi_error(routine, 'MPI_Win_create', mpi_error, .TRUE.)

#ifdef __xlC__
    CALL MPI_Info_free(rma_cache_hint, mpi_error);
    CALL check_mpi_error(routine, 'MPI_Info_free', mpi_error, .TRUE.)
#endif

  END SUBROUTINE openMpiWindow

  SUBROUTINE closeMpiWindow(windowPtr, mpiWindow)
    REAL(dp), POINTER, INTENT(INOUT) :: windowPtr(:)
    INTEGER, INTENT(INOUT) :: mpiWindow

    INTEGER :: mpi_error
    CHARACTER(LEN = *), PARAMETER :: routine = modname//":closeMpiWindow"

#ifdef NOMPI
    CALL finish(routine, "assertion failed: "//routine//"() must NOT be called without MPI")
#else
    ! release RMA window
    CALL MPI_Win_fence(0, mpiWindow, mpi_error)
    CALL check_mpi_error(routine, 'MPI_Win_fence', mpi_error, .FALSE.)
    CALL MPI_Win_free(mpiWindow, mpi_error)
    CALL check_mpi_error(routine, 'MPI_Win_free', mpi_error, .FALSE.)

    ! release RMA memory
    CALL MPI_Free_mem(windowPtr, mpi_error)
    CALL check_mpi_error(routine, 'MPI_Free_mem', mpi_error, .FALSE.)

    ! safety
    mpiWindow = MPI_WIN_NULL
    windowPtr => NULL()
#endif
  END SUBROUTINE closeMpiWindow

  !------------------------------------------------------------------------------------------------
  !
  !  Set global restart attributes.
  !
  SUBROUTINE set_restart_attributes (restartAttributes, restart_args, patch_data)
    TYPE(t_RestartAttributeList), INTENT(INOUT) :: restartAttributes
    TYPE(t_restart_args), INTENT(IN) :: restart_args
    TYPE(t_patch_data), INTENT(IN) :: patch_data(:)

    INTEGER :: i, effectiveDomainCount
    CHARACTER(LEN=*), PARAMETER :: routine = modname//':set_restart_attributes'

#ifdef DEBUG
    WRITE (nerr,FORMAT_VALS3)routine,' is called for p_pe=',p_pe
#endif

    effectiveDomainCount = 1
    IF(patch_data(1)%description%l_opt_ndom) effectiveDomainCount = patch_data(1)%description%opt_ndom
    IF(ALLOCATED(restart_args%output_jfile)) THEN
        CALL setGeneralRestartAttributes(restartAttributes, restart_args%datetime, effectiveDomainCount, restart_args%jstep, &
                                        &restart_args%output_jfile)
    ELSE
        CALL setGeneralRestartAttributes(restartAttributes, restart_args%datetime, effectiveDomainCount, restart_args%jstep)
    END IF

    ! set the domain dependent attributes
    DO i = 1, SIZE(patch_data)
        CALL patch_data(i)%description%setRestartAttributes(restartAttributes)
    END DO
  END SUBROUTINE set_restart_attributes

  !------------------------------------------------------------------------------------------------
  !
  ! Common helper routines to write a restart file.
  !
  !------------------------------------------------------------------------------------------------
  !
  ! Returns true, if the time level of the given field is valid, else false.
  !
  FUNCTION has_valid_time_level(p_info, patchDescription)
    TYPE(t_var_metadata), INTENT(IN) :: p_info
    TYPE(t_restart_patch_description), INTENT(IN) :: patchDescription

    LOGICAL                       :: has_valid_time_level

    INTEGER                       :: idx, time_level
    LOGICAL                       :: lskip_timelev, lskip_extra_timelevs

#ifdef DEBUG
    CHARACTER(LEN=*), PARAMETER   :: routine = modname//':has_valid_time_level'
#endif

    has_valid_time_level = .FALSE.
    IF (.NOT. p_info%lrestart) RETURN
    has_valid_time_level = .TRUE.

    lskip_timelev = .FALSE.
    IF (iequations == INH_ATMOSPHERE .AND. .NOT. (l_limited_area .AND. patchDescription%id == 1)) THEN
      lskip_extra_timelevs = .TRUE.
    ELSE
      lskip_extra_timelevs = .FALSE.
    ENDIF

    ! get time index of the given field
    idx = INDEX(p_info%name,'.TL')
    IF (idx == 0) THEN
      time_level = -1
    ELSE
      time_level = ICHAR(p_info%name(idx+3:idx+3)) - ICHAR('0')

      ! get information about time level to be skipped for current field
      IF (p_info%tlev_source == TLEV_NNOW) THEN
        IF (time_level == patchDescription%nnew) lskip_timelev = .TRUE.
        ! this is needed to skip the extra time levels allocated for nesting
        IF (lskip_extra_timelevs .AND. time_level > 2) lskip_timelev = .TRUE.
      ELSE IF (p_info%tlev_source == TLEV_NNOW_RCF) THEN
        IF (time_level == patchDescription%nnew_rcf) lskip_timelev = .TRUE.
      ENDIF
    ENDIF

    SELECT CASE (iequations)
      CASE(IHS_ATM_TEMP, IHS_ATM_THETA, ISHALLOW_WATER)

        IF ( lskip_timelev                        &
          & .AND. ha_dyn_config%itime_scheme/=LEAPFROG_EXPL &
          & .AND. ha_dyn_config%itime_scheme/=LEAPFROG_SI   ) &
          & has_valid_time_level = .FALSE.
      CASE default
        IF ( lskip_timelev ) has_valid_time_level = .FALSE.
    END SELECT

#ifdef DEBUG
    IF (.NOT. has_valid_time_level) THEN
      WRITE (nerr,'(2a,i3,a,i4,3a)')routine,' p_pe=',p_pe,' time level=',time_level, &
        &                           ' of field=',TRIM(p_info%name),' is invalid'
    ENDIF
#endif

  END FUNCTION has_valid_time_level

  !------------------------------------------------------------------------------------------------
  !
  ! Returns the pointer of the reorder data for the given field.
  !
  FUNCTION get_reorder_ptr(commData, p_info,routine)
    TYPE(t_restart_comm_data), TARGET, INTENT(IN) :: commData
    TYPE(t_var_metadata), INTENT(IN) :: p_info
    CHARACTER(LEN=*), INTENT(IN) :: routine

    TYPE(t_reorder_data), POINTER :: get_reorder_ptr

    get_reorder_ptr => NULL()

    SELECT CASE (p_info%hgrid)
      CASE (GRID_UNSTRUCTURED_CELL)
        get_reorder_ptr => commData%cells
      CASE (GRID_UNSTRUCTURED_EDGE)
        get_reorder_ptr => commData%edges
      CASE (GRID_UNSTRUCTURED_VERT)
        get_reorder_ptr => commData%verts
      CASE default
        CALL finish(routine, UNKNOWN_GRID_TYPE)
    END SELECT

  END FUNCTION get_reorder_ptr

  !------------------------------------------------------------------------------------------------
  !
  ! Check the status of the last netCDF file operation.
  !
  SUBROUTINE check_netcdf_status(status, routine)

    INTEGER, INTENT(IN)           :: status
    CHARACTER(LEN=*), INTENT(IN)  :: routine

    CHARACTER(LEN=128)            :: error_text


    IF (status /= nf_noerr) THEN
      WRITE (error_text, NET_CDF_ERROR_FORMAT)'netCDF error ', status , &
        &                                     ' - ' // nf_strerror(status)
      CALL finish(routine, error_text)
    ENDIF

  END SUBROUTINE check_netcdf_status

  !------------------------------------------------------------------------------------------------
  !
  ! Write restart variable list for a restart PE.
  !
  SUBROUTINE restart_write_var_list(p_pd)
    TYPE(t_patch_data), TARGET, INTENT(INOUT) :: p_pd

    TYPE(t_restart_file), POINTER   :: p_rf
    TYPE(t_var_metadata), POINTER   :: p_info
    TYPE(t_reorder_data), POINTER   :: p_ri
    TYPE(t_var_data), POINTER       :: p_vars(:)

    INTEGER                         :: iv, nval, ierrstat, nlevs, ilev
    INTEGER(KIND=MPI_ADDRESS_KIND)  :: ioff(0:num_work_procs-1)
    REAL(dp), ALLOCATABLE           :: buffer(:,:)
    INTEGER                         :: ichunk, nchunks, chunk_start, chunk_end

    ! For timing
    REAL(dp) :: t_get, t_write
    INTEGER(i8) :: bytesGet, bytesWrite

    CHARACTER(LEN=*), PARAMETER     :: routine = modname//':restart_write_var_list'

#ifdef DEBUG
    WRITE (nerr,FORMAT_VALS3)routine,' p_pe=',p_pe
#endif

    ! check process
    IF (.NOT. my_process_is_restart()) CALL finish(routine, NO_RESTART_PE)

    t_get   = 0.d0
    t_write = 0.d0
    bytesGet = 0_i8
    bytesWrite = 0_i8

    p_rf => p_pd%restart_file

    ! check the contained array of restart variables
    p_vars => p_rf%var_data
    IF (.NOT. ASSOCIATED(p_vars)) RETURN

    ! get maximum number of data points in a slice and allocate tmp. variables
    nval = MAX(p_pd%commData%cells%n_glb, p_pd%commData%edges%n_glb, p_pd%commData%verts%n_glb)

    ! allocate RMA memory
    ALLOCATE(buffer(nval,restart_chunk_size), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, ALLOCATE_FAILED)

    ioff(:) = 0

    ! go over the all restart variables in the associated array
    VAR_LOOP : DO iv = 1, SIZE(p_vars)

      ! get pointer to metadata
      p_info => p_vars(iv)%info

#ifdef DEBUG
      WRITE (nerr,FORMAT_VALS5I)routine,' p_pe=',p_pe,' restart pe processes field=',TRIM(p_info%name)
#endif

      ! check time level of the field
      IF (.NOT. has_valid_time_level(p_info, p_pd%description)) CYCLE

      ! get current level
      IF(p_info%ndims == 2) THEN
        nlevs = 1
      ELSE
        nlevs = p_info%used_dimensions(2)
      ENDIF

      ! get pointer to reorder data
      p_ri => get_reorder_ptr(p_pd%commData, p_info, routine)

      ! no. of chunks of levels (each of size "restart_chunk_size"):
      nchunks = (nlevs-1)/restart_chunk_size + 1
      ! loop over all chunks (of levels)
      LEVELS : DO ichunk=1,nchunks
        chunk_start = (ichunk-1)*restart_chunk_size + 1
        chunk_end = MIN(chunk_start+restart_chunk_size-1, nlevs)
        CALL p_pd%commData%collectData(p_info, chunk_end - chunk_start + 1, buffer, ioff, t_get, bytesGet)

        ! write field content into a file
        t_write = t_write - p_mpi_wtime()
        DO ilev=chunk_start, chunk_end
          CALL streamWriteVarSlice(p_rf%cdiIds%file, p_info%cdiVarID, (ilev-1), buffer(1:p_ri%n_glb, ilev - chunk_start + 1), 0)
          bytesWrite = bytesWrite + p_ri%n_glb*8
        END DO
        t_write = t_write + p_mpi_wtime()

      ENDDO LEVELS

#ifdef DEBUG
      WRITE (nerr,FORMAT_VALS7I)routine,' p_pe=',p_pe,' restart pe writes field=', &
        &                               TRIM(p_info%name),' data=',p_ri%n_glb*nlevs
#endif
    ENDDO VAR_LOOP

    DEALLOCATE(buffer, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, DEALLOCATE_FAILED)

    IF (msg_level >= 12) THEN
      WRITE (0,'(10(a,f10.3))') ' Restart: Got ', REAL(bytesGet, dp)*1.d-6, ' MB, time get: ', t_get, ' s [', &
           & REAL(bytesGet, dp)*1.d-6/MAX(1.e-6_wp, t_get), ' MB/s], time write: ', t_write, ' s [', &
           & REAL(bytesWrite, dp)*1.d-6/MAX(1.e-6_wp,t_write), ' MB/s]'
    ENDIF
  END SUBROUTINE restart_write_var_list

  !------------------------------------------------------------------------------------------------
  !
  ! Write restart variable lists for a compute PE.
  !
  SUBROUTINE compute_write_var_list(p_pd)
    TYPE(t_patch_data), TARGET, INTENT(INOUT) :: p_pd

    TYPE(t_reorder_data), POINTER   :: p_ri
    TYPE(t_var_data), POINTER       :: p_vars(:)
    TYPE(t_ptr_2d), ALLOCATABLE     :: dataPointers(:)
    INTEGER                         :: iv
    INTEGER(i8)                     :: offset
    CHARACTER(LEN=*), PARAMETER     :: routine = modname//':compute_write_var_list'

#ifdef DEBUG
    WRITE (nerr,FORMAT_VALS3)routine,' p_pe=',p_pe
#endif

    ! check process
    IF (.NOT. my_process_is_work()) CALL finish(routine, NO_COMPUTE_PE)

    ! check the array of restart variables
    p_vars => p_pd%restart_file%var_data
    IF (.NOT. ASSOCIATED(p_vars)) RETURN

    ! offset in RMA window for async restart
    offset = 0

    ! go over the all restart variables in the associated array
    DO iv = 1, SIZE(p_vars)
#ifdef DEBUG
      WRITE (nerr,FORMAT_VALS5I)routine,' p_pe=',p_pe,' compute pe processes field=',TRIM(p_vars(iv)%info%name)
#endif

      ! check time level of the field
      IF (.NOT. has_valid_time_level(p_vars(iv)%info, p_pd%description)) CYCLE

      CALL getLevelPointers(p_vars(iv)%info, p_vars(iv)%r_ptr, dataPointers)
      CALL p_pd%commData%postData(p_vars(iv)%info, dataPointers, offset)
      ! no deallocation of dataPointers, so that the next invocation of getLevelPointers() may reuse the last allocation
    END DO
  END SUBROUTINE compute_write_var_list

  !------------------------------------------------------------------------------------------------
  !
  ! Initialize the variable list of the given restart file.
  !
  SUBROUTINE init_restart_variables(p_rf, patchDescription)
    TYPE(t_restart_file), TARGET, INTENT(IN) :: p_rf
    TYPE(t_restart_patch_description), INTENT(IN) :: patchDescription

    TYPE(t_var_data), POINTER     :: p_vars(:)
    TYPE(t_var_metadata), POINTER :: p_info
    CHARACTER(LEN=*), PARAMETER   :: routine = modname//':init_restart_variables'
    INTEGER                       :: iv

#ifdef DEBUG
    WRITE (nerr,FORMAT_VALS3)routine,' p_pe=',p_pe
#endif

    ! check the contained array of restart variables
    p_vars => p_rf%var_data
    IF (.NOT. ASSOCIATED(p_vars)) RETURN

    ! go over the all restart variables in the associated array AND define those that have a valid time level
    DO iv = 1, SIZE(p_vars)
      p_info => p_vars(iv)%info
      IF(has_valid_time_level(p_info, patchDescription)) CALL p_rf%cdiIds%defineVariable(p_info)
    ENDDO
  END SUBROUTINE init_restart_variables

  !------------------------------------------------------------------------------------------------
  !
  ! Opens the restart file from the given parameters.
  !
  SUBROUTINE open_restart_file(p_pd, restart_args, restartAttributes)
    TYPE(t_patch_data), TARGET, INTENT(INOUT) :: p_pd
    TYPE(t_restart_args), INTENT(IN) :: restart_args
    TYPE(t_RestartAttributeList), INTENT(INOUT) :: restartAttributes

    TYPE(t_restart_patch_description), POINTER :: description
    TYPE(t_restart_comm_data), POINTER :: commData
    TYPE(t_restart_file), POINTER :: p_rf
    TYPE(t_var_list), POINTER     :: p_re_list
    INTEGER                       :: restart_type, i
    CHARACTER(LEN=*), PARAMETER   :: routine = modname//':open_restart_file'

#ifdef DEBUG
    WRITE (nerr,FORMAT_VALS3)routine,' p_pe=',p_pe
#endif

    ! just for safety
    IF(.NOT. my_process_is_restart()) CALL finish(routine, NO_RESTART_PE)

    description => p_pd%description
    p_rf => p_pd%restart_file

    ! check the contained array of restart variables
    IF (.NOT. ASSOCIATED(p_rf%var_data)) RETURN

    ! get the first restart list from the global lists
    p_re_list => NULL()
    DO i = 1, nvar_lists
      ! skip, if var_list is not required for restart
      IF (.NOT. var_lists(i)%p%lrestart) CYCLE

      ! save pointer
      p_re_list => var_lists(i)
      EXIT
    ENDDO
    ! assume all restart variables have the same model name
    p_rf%model_type = p_re_list%p%model_type

    ! assume all restart variables uses the same file format
    SELECT CASE (p_re_list%p%restart_type)
      CASE (FILETYPE_NC2, FILETYPE_NC4)
        restart_type = p_re_list%p%restart_type
      CASE default
        CALL finish(routine, UNKNOWN_FILE_FORMAT)
    END SELECT

    p_rf%filename = getRestartFilename(description%base_filename, description%id, restart_args%datetime, p_rf%model_type)

    commData => p_pd%commData
    IF(ALLOCATED(description%opt_pvct)) THEN
        CALL p_rf%cdiIds%openRestartAndCreateIds(TRIM(p_rf%filename), restart_type, restartAttributes, commData%cells%n_glb, &
                                                &commData%verts%n_glb, commData%edges%n_glb, description%cell_type, &
                                                &description%v_grid_defs(1:description%v_grid_count), &
                                                &description%opt_pvct)
    ELSE
        CALL p_rf%cdiIds%openRestartAndCreateIds(TRIM(p_rf%filename), restart_type, restartAttributes, commData%cells%n_glb, &
                                                &commData%verts%n_glb, commData%edges%n_glb, description%cell_type, &
                                                &description%v_grid_defs(1:description%v_grid_count))
    END IF

#ifdef DEBUG
    WRITE (nerr, FORMAT_VALS5)routine,' p_pe=',p_pe,' open netCDF file with ID=',p_rf%cdiIds%file
#endif

    ! init list of restart variables
    CALL init_restart_variables(p_rf, p_pd%description)

    CALL p_rf%cdiIds%finalizeVlist(restart_args%datetime)
  END SUBROUTINE open_restart_file

  !------------------------------------------------------------------------------------------------
  !
  ! Closes the given restart file.
  !
  SUBROUTINE close_restart_file(rf)

    TYPE (t_restart_file), INTENT(INOUT)  :: rf
    CHARACTER(LEN=*), PARAMETER           :: routine = modname//':close_restart_file'

#ifdef DEBUG
    WRITE (nerr,FORMAT_VALS3)routine,' p_pe=',p_pe
#endif

    ! just for safety
    IF(.NOT. my_process_is_restart()) CALL finish(routine, NO_RESTART_PE)

    IF (rf%cdiIds%file /= CDI_UNDEFID) THEN
#ifdef DEBUG
      WRITE (nerr,'(3a)')routine,' try to close restart file=',TRIM(rf%filename)
      WRITE (nerr, FORMAT_VALS5)routine,' p_pe=',p_pe,' close netCDF file with ID=',rf%cdiIds%file
#endif
    ENDIF

    CALL rf%cdiIds%closeAndDestroyIds()
    rf%filename = ''

  END SUBROUTINE close_restart_file

#endif

END MODULE mo_restart_async


