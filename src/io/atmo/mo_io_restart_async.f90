!>
!! Contains routines for asynchronous restart Output
!! --------------------------------------------------------
!!
!! Note: The synchronous implementation of the restart output can be
!!       found in the module "mo_io_restart". See module header for
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

MODULE mo_io_restart_async

  USE mo_decomposition_tools,     ONLY: t_grid_domain_decomp_info
  USE mo_exception,               ONLY: finish, message, message_text, get_filename_noext
  USE mo_fortran_tools,           ONLY: assign_if_present, assign_if_present_allocatable
  USE mo_kind,                    ONLY: wp, i8, dp
  USE mo_datetime,                ONLY: t_datetime, iso8601, iso8601extended
  USE mo_io_units,                ONLY: nerr, filename_max
  USE mo_var_list,                ONLY: nvar_lists, var_lists, new_var_list, delete_var_lists
  USE mo_linked_list,             ONLY: t_list_element, t_var_list
  USE mo_io_restart_attributes,   ONLY: t_RestartAttributeList, RestartAttributeList_make
  USE mo_dynamics_config,         ONLY: iequations
  USE mo_grid_config,             ONLY: l_limited_area
  USE mo_impl_constants,          ONLY: IHS_ATM_TEMP, IHS_ATM_THETA, ISHALLOW_WATER, INH_ATMOSPHERE, &
    &                                   LEAPFROG_EXPL, LEAPFROG_SI, SUCCESS, MAX_CHAR_LENGTH,        &
    &                                   TLEV_NNOW, TLEV_NNOW_RCF
  USE mo_var_metadata_types,      ONLY: t_var_metadata
  USE mo_io_restart_namelist,     ONLY: RestartNamelist_writeToFile, delete_restart_namelists, set_restart_namelist, &
                                      & get_restart_namelist, print_restart_name_lists, RestartNamelist_bcast
#ifdef USE_CRAY_POINTER
  USE mo_name_list_output_init,   ONLY: set_mem_ptr_dp
#endif
  USE mo_communication,           ONLY: idx_no, blk_no
  USE mo_parallel_config,         ONLY: nproma, restart_chunk_size
  USE mo_grid_config,             ONLY: n_dom
  USE mo_run_config,              ONLY: msg_level, restart_filename
  USE mo_ha_dyn_config,           ONLY: ha_dyn_config
  USE mo_model_domain,            ONLY: p_patch, t_patch
  USE mo_cdi,                     ONLY: CDI_UNDEFID, FILETYPE_NC2, FILETYPE_NC4, CDI_GLOBAL, DATATYPE_FLT64, &
                                      & TAXIS_ABSOLUTE, ZAXIS_DEPTH_BELOW_SEA, ZAXIS_GENERIC, ZAXIS_HEIGHT, ZAXIS_HYBRID, &
                                      & ZAXIS_HYBRID_HALF, ZAXIS_LAKE_BOTTOM, ZAXIS_MIX_LAYER, ZAXIS_SEDIMENT_BOTTOM_TW, &
                                      & ZAXIS_SURFACE, ZAXIS_TOA, TIME_VARIABLE, ZAXIS_DEPTH_BELOW_LAND, GRID_UNSTRUCTURED, &
                                      & vlistDefVar, cdiEncodeDate, cdiEncodeTime, streamDefTimestep, gridDestroy, &
                                      & streamWriteVarSlice, streamDefVlist, vlistDefVarDatatype, vlistDefVarName, &
                                      & vlistDefVarLongname, vlistDefVarUnits, vlistDefVarMissval, taxisDefVdate, taxisDefVtime
  USE mo_util_cdi,                ONLY: cdiGetStringError
  USE mo_cdi_constants,           ONLY: GRID_UNSTRUCTURED_EDGE, GRID_UNSTRUCTURED_VERT, GRID_UNSTRUCTURED_CELL
  USE mo_cf_convention
  USE mo_packed_message,          ONLY: t_PackedMessage, kPackOp, kUnpackOp
  USE mo_util_string,             ONLY: t_keyword_list, associate_keyword, with_keywords, &
    &                                   int2string, toCharacter
  USE mo_util_restart,            ONLY: t_v_grid, t_restart_cdi_ids, setGeneralRestartAttributes, &
                                      & setDynamicPatchRestartAttributes, setPhysicsRestartAttributes, create_restart_file_link, &
                                      & t_restart_patch_description, t_var_data

#ifndef NOMPI
  USE mo_mpi,                     ONLY: p_pe, p_pe_work, p_restart_pe0, p_comm_work, p_work_pe0, num_work_procs, MPI_SUCCESS, &
                                      & stop_mpi, p_send, p_recv, p_barrier, p_bcast, my_process_is_restart, my_process_is_work, &
                                      & p_comm_work_2_restart, p_n_work, p_int, process_mpi_restart_size, p_int_i8, p_real_dp, &
                                      & p_comm_work_restart, p_mpi_wtime, p_int_byte, get_my_mpi_work_id, p_comm_rank, &
                                      & process_mpi_all_comm

#ifndef USE_CRAY_POINTER
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_ptr, c_intptr_t, c_f_pointer
#endif

#ifdef __SUNPRO_F95
  INCLUDE "mpif.h"
#else
  USE mpi,                        ONLY: MPI_ADDRESS_KIND, MPI_INFO_NULL, MPI_ROOT, &
    &                                   MPI_LOCK_SHARED, MPI_MODE_NOCHECK, &
    &                                   MPI_PROC_NULL, MPI_WIN_NULL, MPI_LOCK_EXCLUSIVE
#endif
#endif

  IMPLICIT NONE

  INCLUDE 'netcdf.inc'

  ! public routines
  PUBLIC :: set_data_async_restart
  PUBLIC :: prepare_async_restart
  PUBLIC :: write_async_restart
  PUBLIC :: close_async_restart
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
  CHARACTER(LEN=*), PARAMETER :: modname                  = 'shared/mo_io_restart_async/'
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
    INTEGER                     :: cdiTimeIndex
  END TYPE t_restart_file

  !------------------------------------------------------------------------------------------------
  ! TYPE t_reorder_data describes how local cells/edges/verts
  ! have to be reordered to get the global array.
  ! Below, "points" refers to either cells, edges or verts.
  !
  TYPE t_reorder_data
    INTEGER :: n_glb  ! Global number of points per logical patch
    INTEGER :: n_own  ! Number of own points (without halo, only belonging to logical patch)
                      ! Only set on compute PEs, set to 0 on restart PEs
    INTEGER, ALLOCATABLE :: own_idx(:), own_blk(:)
                      ! idx and blk for own points, only set on compute PEs
    INTEGER, ALLOCATABLE :: pe_own(:)
                      ! n_own, gathered for all compute PEs (set on all PEs)
    INTEGER, ALLOCATABLE :: pe_off(:)
                      ! offset of contributions of PEs (set on all PEs)
    INTEGER, ALLOCATABLE :: reorder_index(:)
                      ! Index how to reorder the contributions of all compute PEs
                      ! into the global array (set on all PEs)
  END TYPE t_reorder_data

  ! this combines all the DATA that's relevant for transfering the payload DATA from the worker PEs to the restart PEs
  TYPE t_restart_comm_data
    ! DATA for remote memory access
    INTEGER(i8) :: my_mem_win_off
    INTEGER(i8), ALLOCATABLE :: mem_win_off(:)

    ! reorder data
    TYPE(t_reorder_data) :: cells
    TYPE(t_reorder_data) :: edges
    TYPE(t_reorder_data) :: verts
  END TYPE t_restart_comm_data

  ! combine the DATA that describes a patch for restart purposes with the infos required for the asynchronous fetching of the DATA from the compute PEs
  TYPE t_patch_data
    TYPE(t_restart_patch_description) :: description
    TYPE(t_restart_comm_data) :: commData
    TYPE(t_restart_file) :: restart_file
  END TYPE t_patch_data

  TYPE(t_patch_data), ALLOCATABLE, TARGET :: patch_data(:)

  !------------------------------------------------------------------------------------------------
  ! patch independent arguments
  !
  TYPE t_restart_args
    TYPE(t_datetime)  :: datetime
    INTEGER           :: jstep
    INTEGER, ALLOCATABLE  :: output_jfile(:)
  END TYPE t_restart_args

#ifndef NOMPI
  !------------------------------------------------------------------------------------------------
  ! Currently, we use only 1 MPI window for all output files
  INTEGER mpi_win

  ! MPI memory pointer
#ifdef USE_CRAY_POINTER
  INTEGER (KIND=MPI_ADDRESS_KIND) :: iptr
#else
  TYPE(c_ptr) :: c_mem_ptr
#endif
! Fortran pointer to memory window (REAL*8)
  REAL(dp), POINTER :: mem_ptr_dp(:)

  !------------------------------------------------------------------------------------------------
  ! Broadcast root for intercommunicator broadcasts from compute PEs to restart PEs using
  ! p_comm_work_2_restart.
  INTEGER :: bcast_root

#endif

CONTAINS

  !------------------------------------------------------------------------------------------------
  !
  ! public routines
  !
  !------------------------------------------------------------------------------------------------
  !
  !> Prepare the asynchronous restart (collective call).
  !
  SUBROUTINE prepare_async_restart ()
    CHARACTER(LEN=*), PARAMETER :: routine = modname//'prepare_async_restart'

#ifdef NOMPI
    CALL finish(routine, ASYNC_RESTART_REQ_MPI)
#else

#ifdef DEBUG
    WRITE (nerr,FORMAT_VALS3)routine,' p_pe=',p_pe
#endif

    IF(.NOT. (my_process_is_work() .OR. my_process_is_restart())) RETURN

    ! set broadcast root for intercommunicator broadcasts
    IF(my_process_is_restart()) THEN
      ! root is proc 0 on the compute PEs
      bcast_root = 0
    ELSE
      ! Special root setting for intercommunicators:
      ! The PE really sending must use MPI_ROOT, the others MPI_PROC_NULL.
      IF(p_pe_work == 0) THEN
        bcast_root = MPI_ROOT
      ELSE
        bcast_root = MPI_PROC_NULL
      ENDIF
    ENDIF

    ! transfer restart varlists
    CALL transfer_restart_var_lists

    ! transfer restart namelists
!    CALL RestartNamelist_bcast(0, p_comm_work)
    CALL RestartNamelist_bcast(bcast_root, p_comm_work_2_restart)

    ! create and transfer patch data
    CAll create_and_transfer_patch_data()

    ! init. remote memory access
    CALL init_remote_memory_access

#endif

  END SUBROUTINE prepare_async_restart

  !------------------------------------------------------------------------------------------------
  !
  !> Set patch-dependent dynamic data for asynchronous restart.
  !
  SUBROUTINE set_data_async_restart(patch_id, l_dom_active,       &
                                &   opt_pvct,                     &
                                &   opt_t_elapsed_phy,            &
                                &   opt_lcall_phy,                &
                                &   opt_sim_time,                 &
                                &   opt_ndyn_substeps,            &
                                &   opt_jstep_adv_marchuk_order,  &
                                &   opt_depth,                    &
                                &   opt_depth_lnd,                &
                                &   opt_nlev_snow,                &
                                &   opt_nice_class,               &
                                &   opt_ndom)

    INTEGER,              INTENT(IN)           :: patch_id
    LOGICAL,              INTENT(IN)           :: l_dom_active
                          
    INTEGER,              INTENT(IN), OPTIONAL :: opt_depth
    INTEGER,              INTENT(IN), OPTIONAL :: opt_depth_lnd
    INTEGER,              INTENT(IN), OPTIONAL :: opt_ndyn_substeps
    INTEGER,              INTENT(IN), OPTIONAL :: opt_jstep_adv_marchuk_order
    INTEGER,              INTENT(IN), OPTIONAL :: opt_nlev_snow
    INTEGER,              INTENT(IN), OPTIONAL :: opt_nice_class
    REAL(wp),             INTENT(IN), OPTIONAL :: opt_sim_time
    LOGICAL ,             INTENT(IN), OPTIONAL :: opt_lcall_phy(:)
    REAL(wp),             INTENT(IN), OPTIONAL :: opt_pvct(:)
    REAL(wp),             INTENT(IN), OPTIONAL :: opt_t_elapsed_phy(:)
    INTEGER,              INTENT(IN), OPTIONAL :: opt_ndom            !< no. of domains (appended to symlink name)

    TYPE(t_restart_patch_description), POINTER :: p_pd
    INTEGER                        :: ierrstat, i
    CHARACTER(LEN=*), PARAMETER    :: routine = modname//'set_data_async_restart'

#ifdef NOMPI
    CALL finish(routine, ASYNC_RESTART_REQ_MPI)
#else

#ifdef DEBUG
    WRITE (nerr,FORMAT_VALS3)routine,' p_pe=',p_pe
#endif

    IF(.NOT. (my_process_is_work())) RETURN

    ! find patch
    p_pd => find_patch_description(patch_id, routine)

    ! set activity flag - this needs to be done on all compute PEs because the restart
    ! file may be incomplete otherwise when a nest is started during runtime
    p_pd%l_dom_active = l_dom_active

    ! otherwise, only the first compute PE needs the dynamic restart arguments
    ! in the case of processor splitting, the first PE of the split subset needs them as well
    IF (p_pe_work == 0 .OR. p_pe_work == p_pd%work_pe0_id) THEN

      ! Patch-dependent attributes
      CALL assign_if_present_allocatable(p_pd%opt_pvct, opt_pvct)
      CALL assign_if_present_allocatable(p_pd%opt_t_elapsed_phy, opt_t_elapsed_phy)
      CALL assign_if_present_allocatable(p_pd%opt_lcall_phy, opt_lcall_phy)

      CALL assign_if_present(p_pd%opt_ndyn_substeps, opt_ndyn_substeps)
      p_pd%l_opt_ndyn_substeps = PRESENT(opt_ndyn_substeps)

      CALL assign_if_present(p_pd%opt_jstep_adv_marchuk_order, opt_jstep_adv_marchuk_order)
      p_pd%l_opt_jstep_adv_marchuk_order = PRESENT(opt_jstep_adv_marchuk_order)

      CALL assign_if_present(p_pd%opt_depth, opt_depth)
      p_pd%l_opt_depth = PRESENT(opt_depth)

      CALL assign_if_present(p_pd%opt_depth_lnd, opt_depth_lnd)
      p_pd%l_opt_depth_lnd = PRESENT(opt_depth_lnd)

      CALL assign_if_present(p_pd%opt_nlev_snow, opt_nlev_snow)
      p_pd%l_opt_nlev_snow = PRESENT(opt_nlev_snow)

      CALL assign_if_present(p_pd%opt_nice_class, opt_nice_class)
      p_pd%l_opt_nice_class = PRESENT(opt_nice_class)

      CALL assign_if_present(p_pd%opt_sim_time, opt_sim_time)
      p_pd%l_opt_sim_time = PRESENT(opt_sim_time)

      CALL assign_if_present(p_pd%opt_ndom, opt_ndom)
      p_pd%l_opt_ndom = PRESENT(opt_ndom)
    ENDIF ! (pe_work == 0)
#endif

  END SUBROUTINE set_data_async_restart

  !------------------------------------------------------------------------------------------------
  !
  !> Writes all restart data into one or more files (one file per patch, collective across work processes).
  !
  SUBROUTINE write_async_restart(datetime, jstep, opt_output_jfile)
    TYPE(t_datetime), INTENT(IN) :: datetime
    INTEGER, INTENT(IN) :: jstep
    INTEGER, INTENT(IN), OPTIONAL :: opt_output_jfile(:)

    INTEGER :: idx
    TYPE(t_restart_args) :: restart_args
    CHARACTER(LEN=*), PARAMETER :: routine = modname//'write_async_restart'

#ifdef NOMPI
    CALL finish (routine, ASYNC_RESTART_REQ_MPI)
#else

#ifdef DEBUG
    WRITE (nerr,FORMAT_VALS3)routine,' p_pe=',p_pe
#endif

    ! check kind of process
    IF(.NOT. my_process_is_work()) RETURN

    CALL compute_wait_for_restart
    restart_args%datetime = datetime
    restart_args%jstep = jstep
    ! otherwise, only the first compute PE needs the dynamic restart arguments
    ! in the case of processor splitting, the first PE of the split subset needs them as well
    CALL assign_if_present_allocatable(restart_args%output_jfile, opt_output_jfile)
    CALL compute_start_restart(restart_args)

    ! do the restart output
    DO idx = 1, SIZE(patch_data)
      ! collective call to write the restart variables
      IF(patch_data(idx)%description%l_dom_active) CALL compute_write_var_list(patch_data(idx))
    END DO
#endif
  END SUBROUTINE write_async_restart

  !> Writes all restart data into one or more files (one file per patch, collective across restart processes).
  SUBROUTINE restart_write_async_restart(restart_args)
    TYPE(t_restart_args), INTENT(IN) :: restart_args

    TYPE(t_patch_data), POINTER :: p_pd
    TYPE(t_restart_patch_description), POINTER :: description
    INTEGER :: idx
    TYPE(t_RestartAttributeList), POINTER :: restartAttributes
    CHARACTER(LEN=*), PARAMETER :: routine = modname//'restart_write_async_restart'

#ifdef NOMPI
    CALL finish (routine, ASYNC_RESTART_REQ_MPI)
#else

#ifdef DEBUG
    WRITE (nerr,FORMAT_VALS3)routine,' p_pe=',p_pe
#endif

    ! check kind of process
    IF(.NOT. my_process_is_restart()) RETURN

    ! do the restart output
    DO idx = 1, SIZE(patch_data)
      p_pd => patch_data(idx)
      description => p_pd%description

      ! check if the patch is actice
      IF (.NOT. description%l_dom_active) CYCLE

      ! consider the right restart process
      IF (p_pe == description%restart_proc_id) THEN
        ! set global restart attributes/lists
        restartAttributes => RestartAttributeList_make()
        CALL set_restart_attributes(restartAttributes, restart_args)
        CALL description%defineVGrids()

#ifdef DEBUG
        CALL print_restart_arguments()
        CALL restartAttributes%printAttributes()
        CALL print_restart_name_lists()
#endif
        CALL open_restart_file(p_pd, restart_args, restartAttributes)

        CALL restartAttributes%destruct()
        DEALLOCATE(restartAttributes)

        ! collective call to write the restart variables
        CALL restart_write_var_list(p_pd, restart_args)
        IF(description%l_opt_ndom) THEN
            CALL create_restart_file_link(TRIM(p_pd%restart_file%filename), TRIM(p_pd%restart_file%model_type), &
                                         &description%restart_proc_id - p_restart_pe0, description%id, &
                                         &opt_ndom = description%opt_ndom)
        ELSE
            CALL create_restart_file_link(TRIM(p_pd%restart_file%filename), TRIM(p_pd%restart_file%model_type), &
                                         &description%restart_proc_id - p_restart_pe0, description%id)
        END IF
        CALL close_restart_file(p_pd%restart_file)
      ENDIF
    ENDDO
#endif
  END SUBROUTINE restart_write_async_restart

  !------------------------------------------------------------------------------------------------
  !
  !> Closes asynchronous restart (collective call).
  !
  SUBROUTINE close_async_restart

#ifdef NOMPI
    CHARACTER(LEN=*), PARAMETER :: routine = modname//'close_async_restart'

    CALL finish(routine, ASYNC_RESTART_REQ_MPI)
#else
    ! check kind of process
    IF (.NOT. my_process_is_work() .AND. .NOT. my_process_is_restart()) RETURN

    IF (my_process_is_work()) THEN

      CALL compute_wait_for_restart
      CALL compute_shutdown_restart

    ENDIF

    CALL release_resources
#endif

  END SUBROUTINE close_async_restart

  !-------------------------------------------------------------------------------------------------
  !>
  !! Main routine for restart PEs.
  !! Please note that this routine never returns.
  SUBROUTINE restart_main_proc
    LOGICAL :: done
    TYPE(t_restart_args) :: restart_args
    CHARACTER(LEN=*), PARAMETER :: routine = modname//'restart_main_proc'

#ifdef NOMPI
    CALL finish(routine, ASYNC_RESTART_REQ_MPI)
#else

#ifdef DEBUG
    WRITE (nerr,FORMAT_VALS3)routine,' p_pe=',p_pe
#endif

    ! check kind of process
    IF (.NOT. my_process_is_restart()) RETURN

    ! prepare restart (collective call)
    CALL prepare_async_restart()

    ! tell the compute PEs that we are ready to work
    CALL restart_send_ready

    ! enter restart loop
    done = .FALSE.
    DO
      ! wait for a message from the compute PEs to start
      CALL restart_wait_for_start(restart_args, done)

      IF(done) EXIT ! leave loop, we are done

      ! read and write restart variable lists (collective call)
      CALL restart_write_async_restart(restart_args)

      ! inform compute PEs that the restart is done
      CALL restart_send_ready
    ENDDO

    ! finalization sequence (collective call)
    CALL close_async_restart

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
    CHARACTER(LEN=*), PARAMETER :: routine = modname//'restart_send_ready'

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

  SUBROUTINE restartMetadataPacker(operation, restart_args, message)
    INTEGER, VALUE :: operation
    TYPE(t_restart_args), INTENT(INOUT) :: restart_args
    TYPE(t_PackedMessage), INTENT(INOUT) :: message

    INTEGER :: i, calday, patchId
    CHARACTER(LEN = *), PARAMETER :: routine = modname//":restartMetadataPacker"

    ! set patch independent arguments
    CALL message%execute(operation, restart_args%datetime%year)
    CALL message%execute(operation, restart_args%datetime%month)
    CALL message%execute(operation, restart_args%datetime%day)
    CALL message%execute(operation, restart_args%datetime%hour)
    CALL message%execute(operation, restart_args%datetime%minute)
    CALL message%execute(operation, restart_args%datetime%second)
    CALL message%execute(operation, restart_args%datetime%caltime)
    calday = INT(restart_args%datetime%calday)
    CALL message%execute(operation, calday)
    restart_args%datetime%calday = INT(calday,i8)
    CALL message%execute(operation, restart_args%datetime%daysec)
    CALL message%execute(operation, restart_args%jstep)
    CALL message%execute(operation, restart_args%output_jfile)

    ! set data of all patches
    DO i = 1, SIZE(patch_data)
        CALL patch_data(i)%description%packer(operation, message)
    END DO
  END SUBROUTINE restartMetadataPacker

  !-------------------------------------------------------------------------------------------------
  !>
  !! restart_wait_for_start: Wait for a message from compute PEs that we should start restart or finish.
  !! The counterpart on the compute side is compute_start_restart/compute_shutdown_restart.
  !
  SUBROUTINE restart_wait_for_start(restart_args, done)
    TYPE(t_restart_args), INTENT(INOUT) :: restart_args
    LOGICAL, INTENT(OUT)           :: done ! flag if we should shut down

    INTEGER :: iheader
    TYPE(t_PackedMessage) :: message
    CHARACTER(LEN=*), PARAMETER :: routine = modname//'restart_wait_for_start'

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
        CALL restartMetadataPacker(kUnpackOp, restart_args, message)
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

    CHARACTER(LEN=*), PARAMETER :: routine = modname//'compute_wait_for_restart'

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

  !-------------------------------------------------------------------------------------------------
  !>
  !! compute_start_restart: Send a message to restart PEs that they should start restart.
  !! The counterpart on the restart side is restart_wait_for_start.
  !
  SUBROUTINE compute_start_restart(restart_args)
    TYPE(t_restart_args), INTENT(INOUT) :: restart_args

    TYPE(t_restart_patch_description), POINTER :: curPatch
    INTEGER :: i, trash
    TYPE(t_PackedMessage) :: message
    CHARACTER(LEN=*), PARAMETER :: routine = modname//'compute_start_restart'

#ifdef DEBUG
    WRITE (nerr,FORMAT_VALS5)routine,' p_pe=',p_pe, ' call p_barrier with communicator=',p_comm_work
#endif

    ! make sure all are here
    CALL p_barrier(comm=p_comm_work)

    ! if processor splitting is applied, the time-dependent data need to be transferred
    ! from the subset master PE to PE0, from where they are communicated to the output PE(s)
    CALL message%construct()
    DO i = 1, SIZE(patch_data)
      curPatch => patch_data(i)%description
      CALL message%reset()

      CALL curPatch%setTimeLevels() ! copy the global variables (nold, ...) to the patch description

      IF (curPatch%work_pe0_id /= 0) THEN   ! nothing to communicate IF PE0 IS already the subset master
        IF (p_pe_work == 0) THEN
          ! recieve the package for this patch
          CALL message%recv(curPatch%work_pe0_id, 0, process_mpi_all_comm)
          CALL curPatch%packer(kUnpackOp, message)
        ELSE IF (p_pe_work == curPatch%work_pe0_id) THEN
          ! send the time dependent DATA to process 0
          CALL curPatch%packer(kPackOp, message)
          CALL message%send(0, 0, process_mpi_all_comm)
        END IF
      END IF
    END DO

    CALL message%reset()
    IF(p_pe_work == 0) THEN
      ! send the DATA to the restart master
      CALL message%pack(MSG_RESTART_START)  ! set command id
      CALL restartMetadataPacker(kPackOp, restart_args, message)    ! all the other DATA
      CALL message%send(p_restart_pe0, 0, process_mpi_all_comm)
    ENDIF
    ! broadcast a copy among the compute processes, so that all processes have a consistent view of the restart patch descriptions
    CALL message%bcast(0, p_comm_work)
    CALL message%unpack(trash)  ! ignore the command id
    CALL restartMetadataPacker(kUnpackOp, restart_args, message)

    CALL message%destruct()
  END SUBROUTINE compute_start_restart

  !-------------------------------------------------------------------------------------------------
  !>
  !! compute_shutdown_restart: Send a message to restart PEs that they should shut down.
  !! The counterpart on the restart side is restart_wait_for_start.
  !
  SUBROUTINE compute_shutdown_restart
    TYPE(t_PackedMessage) :: message
    CHARACTER(LEN=*), PARAMETER :: routine = modname//'compute_shutdown_restart'

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
  !  Release dynamically created memory of reorder data.
  !
  SUBROUTINE release_reorder_data (rd)

    TYPE(t_reorder_data), INTENT(INOUT) :: rd

    IF (ALLOCATED(rd%own_idx))       DEALLOCATE(rd%own_idx)
    IF (ALLOCATED(rd%own_blk))       DEALLOCATE(rd%own_blk)
    IF (ALLOCATED(rd%pe_own))        DEALLOCATE(rd%pe_own)
    IF (ALLOCATED(rd%pe_off))        DEALLOCATE(rd%pe_off)
    IF (ALLOCATED(rd%reorder_index)) DEALLOCATE(rd%reorder_index)

  END SUBROUTINE release_reorder_data

  !------------------------------------------------------------------------------------------------
  !
  !  Release resources of a restart file.
  !
  SUBROUTINE release_restart_file (rf)
    TYPE(t_restart_file), INTENT(INOUT) :: rf

    CHARACTER(LEN=*), PARAMETER   :: routine = modname//'release_restart_file'

#ifdef DEBUG
    WRITE (nerr,FORMAT_VALS3)routine,' is called for p_pe=',p_pe
#endif

    IF (ASSOCIATED(rf%var_data))   DEALLOCATE(rf%var_data)
    IF (my_process_is_restart()) rf%cdiTimeIndex = CDI_UNDEFID
  END SUBROUTINE release_restart_file

  !------------------------------------------------------------------------------------------------
  !
  !  Release all resource of the restart process.
  !
  SUBROUTINE release_resources

    TYPE(t_patch_data), POINTER   :: p_pd
    INTEGER                       :: idx, mpi_error

    CHARACTER(LEN=*), PARAMETER   :: routine = modname//'release_resources'

#ifdef DEBUG
    WRITE (nerr,FORMAT_VALS3)routine,' is called for p_pe=',p_pe
#endif

    ! release patch data
    IF (ALLOCATED(patch_data)) THEN
      DO idx = 1, SIZE(patch_data)
        p_pd => patch_data(idx)

        ! release patch arrays
        IF (ASSOCIATED(p_pd%description%v_grid_defs))      DEALLOCATE(p_pd%description%v_grid_defs)
        p_pd%description%v_grid_count = 0
        IF (ALLOCATED(p_pd%description%opt_pvct))          DEALLOCATE(p_pd%description%opt_pvct)
        IF (ALLOCATED(p_pd%description%opt_lcall_phy))     DEALLOCATE(p_pd%description%opt_lcall_phy)
        IF (ALLOCATED(p_pd%description%opt_t_elapsed_phy)) DEALLOCATE(p_pd%description%opt_t_elapsed_phy)

        ! release restart file data
        CALL release_restart_file(p_pd%restart_file)

        ! release communication data
        IF (ALLOCATED(p_pd%commData%mem_win_off)) DEALLOCATE(p_pd%commData%mem_win_off)
        CALL release_reorder_data(p_pd%commData%cells)
        CALL release_reorder_data(p_pd%commData%verts)
        CALL release_reorder_data(p_pd%commData%edges)
      END DO

      DEALLOCATE(patch_data)
    END IF

    ! release RMA window
    IF (mpi_win /= MPI_WIN_NULL) THEN
      CALL MPI_Win_fence(0, mpi_win, mpi_error)
      CALL check_mpi_error(routine, 'MPI_Win_fence', mpi_error, .FALSE.)
      CALL MPI_Win_free(mpi_win, mpi_error)
      CALL check_mpi_error(routine, 'MPI_Win_free', mpi_error, .FALSE.)
      mpi_win = MPI_WIN_NULL
    END IF

    ! release RMA memory
    CALL MPI_Free_mem(mem_ptr_dp, mpi_error)
    CALL check_mpi_error(routine, 'MPI_Free_mem', mpi_error, .FALSE.)

  END SUBROUTINE release_resources

  !------------------------------------------------------------------------------------------------
  !
  ! Common helper routines to processing lists.
  !

  !------------------------------------------------------------------------------------------------
  !
  !  Print restart arguments.
  !
  SUBROUTINE print_restart_arguments(restart_args)
    TYPE(t_restart_args), INTENT(IN) :: restart_args
    CHARACTER(LEN=*), PARAMETER   :: routine = modname//'print_restart_arguments'

    WRITE (nerr,FORMAT_VALS3)routine,' is called for p_pe=',p_pe

    PRINT *,routine, ' current_caltime=', restart_args%datetime%caltime
    PRINT *,routine, ' current_calday=',  restart_args%datetime%calday
    PRINT *,routine, ' current_daysec=',  restart_args%datetime%daysec

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
    CHARACTER(LEN=*), PARAMETER   :: routine = modname//'get_var_list_number'

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

    CHARACTER(LEN=*), PARAMETER :: routine = modname//'restartVarlistPacker'

    ! delete old var lists
    IF(operation == kUnpackOp) CALL delete_var_lists

    ! get the size - in default INTEGER words - which is needed to
    ! hold the contents of TYPE(t_var_metadata)
    info_size = SIZE(TRANSFER(info, (/ 0 /)))
    ALLOCATE(info_storage(info_size), STAT=ierrstat)
    IF(ierrstat /= SUCCESS) CALL finish(routine, ALLOCATE_FAILED)

    ! get the number of var lists
    nv = nvar_lists
    CALL message%execute(operation, nv)

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
        CALL message%execute(operation, lrestart)
        CALL message%execute(operation, var_list_name)
        CALL message%execute(operation, model_type)
        CALL message%execute(operation, patch_id)
        CALL message%execute(operation, restart_type)
        CALL message%execute(operation, vlevel_type)
        CALL message%execute(operation, nelems)

        IF(.NOT. lrestart) CYCLE  ! transfer only a restart var_list
        IF(nelems == 0) CYCLE ! check if there are valid restart fields

        IF(operation == kPackOp) THEN
            element => var_lists(iv)%p%first_list_element
            DO
                IF(.NOT. ASSOCIATED(element)) EXIT
                IF(element%field%info%lrestart) THEN
                    info_storage = TRANSFER(element%field%info, (/ 0 /))
                    CALL message%execute(operation, info_storage)
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
                CALL message%execute(operation, info_storage)
                element%field%info = TRANSFER(info_storage, info)
            END DO
        END IF

        DEALLOCATE(info_storage)
    END DO
  END SUBROUTINE restartVarlistPacker

  !
  ! Transfers the restart var lists from the worker to the restart PEs.
  !
  SUBROUTINE transfer_restart_var_lists()
    TYPE(t_PackedMessage) :: message
    CHARACTER(LEN=*), PARAMETER :: routine = modname//'transfer_restart_var_lists'

    CALL message%construct()
    IF(.NOT.my_process_is_restart()) CALL restartVarlistPacker(kPackOp, message)
    CALL message%bcast(bcast_root, p_comm_work_2_restart)
    IF(my_process_is_restart()) CALL restartVarlistPacker(kUnpackOp, message)
    CALL message%destruct()
  END SUBROUTINE transfer_restart_var_lists

  ! collective across restart AND worker PEs
  SUBROUTINE create_patch_description(description, domain)
    TYPE(t_restart_patch_description), INTENT(INOUT) :: description
    INTEGER, VALUE :: domain

    TYPE(t_PackedMessage) :: message

    ! initialize on work PEs
    CALL message%construct()
    IF(my_process_is_work()) THEN
        CALL description%setPatch(p_patch(domain))
        CALL description%packer(kPackOp, message)
    END IF

    ! transfer data to restart PEs
    CALL message%bcast(bcast_root, p_comm_work_2_restart)
    CALL description%packer(kUnpackOp, message)

    ! initialize the fields that we DO NOT communicate from the worker PEs to the restart PEs
    description%restart_proc_id = MOD(description%id-1, process_mpi_restart_size) + p_restart_pe0
    description%v_grid_defs => NULL()
    description%v_grid_count = 0

    CALL message%destruct() ! cleanup
  END SUBROUTINE create_patch_description

  ! collective across restart AND worker PEs
  SUBROUTINE createCommData(commData, domain)
    TYPE(t_restart_comm_data), INTENT(INOUT) :: commData
    INTEGER, VALUE :: domain

    ! initialize on work PEs
    IF(my_process_is_work()) THEN
        CALL set_reorder_data(p_patch(domain)%n_patch_cells_g, p_patch(domain)%n_patch_cells, p_patch(domain)%cells%decomp_info, &
                             &commData%cells)
        CALL set_reorder_data(p_patch(domain)%n_patch_edges_g, p_patch(domain)%n_patch_edges, p_patch(domain)%edges%decomp_info, &
                             &commData%edges)
        CALL set_reorder_data(p_patch(domain)%n_patch_verts_g, p_patch(domain)%n_patch_verts, p_patch(domain)%verts%decomp_info, &
                             &commData%verts)
    END IF

    ! transfer data to restart PEs
    CALL transfer_reorder_data(commData%cells)
    CALL transfer_reorder_data(commData%edges)
    CALL transfer_reorder_data(commData%verts)

    commData%my_mem_win_off = 0_i8
  END SUBROUTINE createCommData

  !-------------------------------------------------------------------------------------------------
  !
  ! Create patch data and transfers this data from the worker to the restart PEs.
  !
  SUBROUTINE create_and_transfer_patch_data()
    INTEGER                        :: jg, ierrstat
    CHARACTER(LEN=*), PARAMETER :: routine = modname//'create_and_transfer_patch_data'

#ifdef DEBUG
    WRITE (nerr,FORMAT_VALS3)routine,' p_pe=',p_pe
#endif

    ! replicate domain setup
    CALL p_bcast(n_dom, bcast_root, p_comm_work_2_restart)

    ! allocate patch data structure
    ALLOCATE(patch_data(n_dom), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, ALLOCATE_FAILED)

    ! initialize the patch_data structures
    DO jg = 1, n_dom
        CALL create_patch_description(patch_data(jg)%description, jg)
        CALL createCommData(patch_data(jg)%commData, jg)
        CALL set_restart_file_data(patch_data(jg)%restart_file, jg)
    END DO
  END SUBROUTINE create_and_transfer_patch_data

  !------------------------------------------------------------------------------------------------
  !
  ! Sets the restart file data with the given logical patch ident.
  !
  SUBROUTINE set_restart_file_data (rf, patch_id)

    TYPE (t_restart_file),  INTENT (INOUT) :: rf
    INTEGER,               INTENT (IN)    :: patch_id

    INTEGER                               :: ierrstat, i, i2, num_vars
    TYPE (t_list_element), POINTER        :: element

    CHARACTER(LEN=*), PARAMETER           :: routine = modname//'set_restart_file_data'

#ifdef DEBUG
    WRITE (nerr,FORMAT_VALS5)routine,' is called for p_pe=',p_pe,' patch_id=',patch_id
#endif

    ! init. main variables
    rf%var_data                   => NULL()
    rf%filename                   = ''

    CALL rf%cdiIds%init()
    rf%cdiTimeIndex               = CDI_UNDEFID

    ! counts number of restart variables for this file (logical patch ident)
    num_vars = 0
    CALL get_var_list_number(num_vars, patch_id)

#ifdef DEBUG
    WRITE (nerr,FORMAT_VALS3)routine,' numvars=',num_vars
#endif

    IF (num_vars <= 0) RETURN
    ! allocate the array of restart variables
    ALLOCATE (rf%var_data(num_vars), STAT=ierrstat)
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
          rf%var_data(i2)%info = element%field%info
          IF (my_process_is_work() .AND. ASSOCIATED(element%field%r_ptr)) THEN
            rf%var_data(i2)%r_ptr => element%field%r_ptr
          ELSE
            rf%var_data(i2)%r_ptr => NULL()
          ENDIF
        ENDIF
        element => element%next_list_element
      ENDDO
    ENDDO

  END SUBROUTINE set_restart_file_data

  !------------------------------------------------------------------------------------------------
  !
  ! Sets the reorder data for cells/edges/verts.
  !
  SUBROUTINE set_reorder_data(n_points_g, n_points, decomp_info, reo)
    INTEGER, INTENT(IN) :: n_points_g      ! Global number of cells/edges/verts in logical patch
    INTEGER, INTENT(IN) :: n_points        ! Local number of cells/edges/verts in logical patch
    TYPE(t_grid_domain_decomp_info), INTENT(IN) :: decomp_info

    TYPE(t_reorder_data), INTENT(INOUT) :: reo ! Result: reorder data

    INTEGER :: i, n, il, ib, mpi_error, ierrstat
    LOGICAL, ALLOCATABLE :: owner_mask_1d(:) ! non-blocked owner mask for (logical) patch
    INTEGER, ALLOCATABLE :: glbidx_own(:), glbidx_glb(:), reorder_index_log_dom(:)
    CHARACTER (LEN=MAX_ERROR_LENGTH) :: error_message

    CHARACTER(LEN=*), PARAMETER :: routine = modname//'set_reorder_data'

#ifdef DEBUG
    WRITE (nerr,FORMAT_VALS3)routine,' is called for p_pe=',p_pe
#endif

    ! just for safety
    IF(my_process_is_restart()) CALL finish(routine, NO_COMPUTE_PE)

    ! set the non-blocked patch owner mask
    ALLOCATE(owner_mask_1d(n_points))
    DO i = 1, n_points
      il = idx_no(i)
      ib = blk_no(i)
      owner_mask_1d(i) = decomp_info%owner_mask(il,ib)
    ENDDO

    ! get number of owned cells/edges/verts (without halos)
    reo%n_own = COUNT(owner_mask_1d(:))

    ! set index arrays to own cells/edges/verts
    ALLOCATE(reo%own_idx(reo%n_own), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, ALLOCATE_FAILED)
    ALLOCATE(reo%own_blk(reo%n_own), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, ALLOCATE_FAILED)

    ! global index of my own points
    ALLOCATE(glbidx_own(reo%n_own), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, ALLOCATE_FAILED)

    n = 0
    DO i = 1, n_points
      IF(owner_mask_1d(i)) THEN
        n = n+1
        reo%own_idx(n) = idx_no(i)
        reo%own_blk(n) = blk_no(i)
        glbidx_own(n)  = decomp_info%glb_index(i)
      ENDIF
    ENDDO

    ! gather the number of own points for every PE into reo%pe_own
    ALLOCATE(reo%pe_own(0:p_n_work-1), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, ALLOCATE_FAILED)
    ALLOCATE(reo%pe_off(0:p_n_work-1), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, ALLOCATE_FAILED)

    CALL MPI_Allgather(reo%n_own,  1, p_int, &
                       reo%pe_own, 1, p_int, &
                       p_comm_work, mpi_error)
    CALL check_mpi_error(routine, 'MPI_Allgather', mpi_error, .TRUE.)

    ! get offset within result array
    reo%pe_off(0) = 0
    DO i = 1, p_n_work-1
      reo%pe_off(i) = reo%pe_off(i-1) + reo%pe_own(i-1)
    ENDDO

    ! get global number of points for current patch
    reo%n_glb = SUM(reo%pe_own(:))

    ! Get the global index numbers of the data when it is gathered on PE 0
    ! exactly in the same order as it is retrieved later during restart.
    ALLOCATE(glbidx_glb(reo%n_glb), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, ALLOCATE_FAILED)

    CALL MPI_Allgatherv(glbidx_own, reo%n_own, p_int, &
                        glbidx_glb, reo%pe_own, reo%pe_off, p_int, &
                        p_comm_work, mpi_error)
    CALL check_mpi_error(routine, 'MPI_Allgatherv', mpi_error, .TRUE.)

    ! get reorder_index
    ALLOCATE(reo%reorder_index(reo%n_glb), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, ALLOCATE_FAILED)

    ! spans the complete logical domain
    ALLOCATE(reorder_index_log_dom(n_points_g), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, ALLOCATE_FAILED)
    reorder_index_log_dom(:) = 0

    DO i = 1, reo%n_glb
      ! reorder_index_log_dom stores where a global point in logical domain comes from.
      ! It is nonzero only at the patch locations
      reorder_index_log_dom(glbidx_glb(i)) = i
    ENDDO

    ! gather the reorder index
    n = 0
    DO i = 1, n_points_g
      IF(reorder_index_log_dom(i)>0) THEN
        n = n+1
        reo%reorder_index(n) = reorder_index_log_dom(i)
      ENDIF
    ENDDO

#ifdef DEBUG
  WRITE (nerr, '(a,i8,a,i8)') 'Reordering number: n=',n, &
         & ', reo%n_glb=',reo%n_glb
#endif

    ! safety check
    IF(n/=reo%n_glb) THEN
      WRITE (error_message, '(a,i8,a,i8)') 'Reordering failed: n=',n, &
                                         & ' /= reo%n_glb=',reo%n_glb
      CALL finish(routine,TRIM(error_message))
    ENDIF

    DEALLOCATE(owner_mask_1d)
    DEALLOCATE(glbidx_own)
    DEALLOCATE(glbidx_glb)
    DEALLOCATE(reorder_index_log_dom)

  END SUBROUTINE set_reorder_data

  !------------------------------------------------------------------------------------------------
  !
  ! Transfers reorder data to restart PEs.
  !
  SUBROUTINE transfer_reorder_data(reo)

    TYPE(t_reorder_data), INTENT(INOUT) :: reo
    INTEGER                             :: ierrstat

    CHARACTER(LEN=*), PARAMETER :: routine = modname//'transfer_reorder_data'

    ! transfer the global number of points, this is not yet known on restart PEs
    CALL p_bcast(reo%n_glb,  bcast_root, p_comm_work_2_restart)

    IF(my_process_is_restart()) THEN

      ! on restart PEs: n_own = 0, own_idx and own_blk are not allocated
      reo%n_own = 0

      ! pe_own must be allocated for num_work_procs, not for p_n_work
      ALLOCATE(reo%pe_own(0:num_work_procs-1), STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, ALLOCATE_FAILED)

      ALLOCATE(reo%pe_off(0:num_work_procs-1), STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, ALLOCATE_FAILED)

      ALLOCATE(reo%reorder_index(reo%n_glb), STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, ALLOCATE_FAILED)
    ENDIF

    CALL p_bcast(reo%pe_own, bcast_root, p_comm_work_2_restart)
    CALL p_bcast(reo%pe_off, bcast_root, p_comm_work_2_restart)
    CALL p_bcast(reo%reorder_index, bcast_root, p_comm_work_2_restart)

  END SUBROUTINE transfer_reorder_data

  !------------------------------------------------------------------------------------------------
  !
  ! Initializes the remote memory access for asynchronous restart.
  !
  SUBROUTINE init_remote_memory_access
    INTEGER :: i, iv, nlevs, ierrstat
    INTEGER :: nbytes_real, mpi_error, rma_cache_hint
    INTEGER (KIND=MPI_ADDRESS_KIND) :: mem_size, mem_bytes
#ifdef USE_CRAY_POINTER
    REAL(dp) :: tmp_dp
    POINTER(tmp_ptr_dp,tmp_dp(*))
#endif
    TYPE(t_var_data), POINTER :: p_vars(:)
    TYPE(t_restart_comm_data), POINTER :: commData
    CHARACTER(LEN=*), PARAMETER :: routine = modname//'init_remote_memory_access'

#ifdef DEBUG
    WRITE (nerr,FORMAT_VALS3)routine,' is called for p_pe=',p_pe
#endif

    mpi_win = MPI_WIN_NULL

    ! get size and offset of the data for every restart file
    mem_size = 0_i8

    ! go over all patches
    DO i = 1, SIZE(patch_data)
      commData => patch_data(i)%commData
      commData%my_mem_win_off = mem_size

      ! go over all restart variables for this restart file
      p_vars => patch_data(i)%restart_file%var_data
      IF (.NOT. ASSOCIATED(p_vars)) CYCLE
      DO iv = 1, SIZE(p_vars)

        IF (p_vars(iv)%info%ndims == 2) THEN
          nlevs = 1
        ELSE
          nlevs = p_vars(iv)%info%used_dimensions(2)
        ENDIF

        SELECT CASE (p_vars(iv)%info%hgrid)
          CASE (GRID_UNSTRUCTURED_CELL)
            mem_size = mem_size + INT(nlevs*commData%cells%n_own,i8)
          CASE (GRID_UNSTRUCTURED_EDGE)
            mem_size = mem_size + INT(nlevs*commData%edges%n_own,i8)
          CASE (GRID_UNSTRUCTURED_VERT)
            mem_size = mem_size + INT(nlevs*commData%verts%n_own,i8)
          CASE DEFAULT
            CALL finish(routine,UNKNOWN_GRID_TYPE)
        END SELECT

      ENDDO

      ! get the offset on all PEs
      ALLOCATE(commData%mem_win_off(0:num_work_procs-1), STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, ALLOCATE_FAILED)
      IF(.NOT.my_process_is_restart()) THEN
        CALL MPI_Allgather(commData%my_mem_win_off, 1, p_int_i8, &
                           commData%mem_win_off, 1, p_int_i8,    &
                           p_comm_work, mpi_error)
        CALL check_mpi_error(routine, 'MPI_Allgather', mpi_error, .TRUE.)
      ENDIF

      CALL p_bcast(commData%mem_win_off, bcast_root, p_comm_work_2_restart)

    ENDDO

    ! mem_size is calculated as number of variables above, get number of bytes
    ! get the amount of bytes per REAL*8 variable (as used in MPI communication)
    CALL MPI_Type_extent(p_real_dp, nbytes_real, mpi_error)
    CALL check_mpi_error(routine, 'MPI_Type_extent', mpi_error, .TRUE.)

    ! for the restart PEs the amount of memory needed is 0 - allocate at least 1 word there:
    mem_bytes = MAX(mem_size,1_i8)*INT(nbytes_real,i8)

    ! allocate amount of memory needed with MPI_Alloc_mem
    ! 
    ! Depending on wether the Fortran 2003 C interoperability features
    ! are available, one needs to use non-standard language extensions
    ! for calls from Fortran, namely Cray Pointers, since
    ! MPI_Alloc_mem wants a C pointer argument.
    !
    ! see, for example: http://www.lrz.de/services/software/parallel/mpi/onesided/
#ifdef USE_CRAY_POINTER
    CALL MPI_Alloc_mem(mem_bytes, MPI_INFO_NULL, iptr, mpi_error)
    CALL check_mpi_error(routine, 'MPI_Alloc_mem', mpi_error, .TRUE.)

    tmp_ptr_dp = iptr
    CALL set_mem_ptr_dp(tmp_dp, INT(mem_size))
#else
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
    CALL C_F_POINTER(c_mem_ptr, mem_ptr_dp, (/ INT(mem_size) /) )
#else
    CALL C_F_POINTER(c_mem_ptr, mem_ptr_dp, (/ mem_size /) )
#endif
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
    CALL MPI_Win_create(mem_ptr_dp,mem_bytes,nbytes_real,MPI_INFO_NULL,&
      &                 p_comm_work_restart,mpi_win,mpi_error )
    CALL check_mpi_error(routine, 'MPI_Win_create', mpi_error, .TRUE.)

#ifdef __xlC__
    CALL MPI_Info_free(rma_cache_hint, mpi_error);
    CALL check_mpi_error(routine, 'MPI_Info_free', mpi_error, .TRUE.)
#endif

  END SUBROUTINE init_remote_memory_access

  !------------------------------------------------------------------------------------------------
  !
  ! Common helper routines to init. the restart.
  !
  !------------------------------------------------------------------------------------------------
  !
  !  Find the patch of the given id.
  !
  FUNCTION find_patch_description(id, routine) RESULT(RESULT)
    INTEGER, INTENT(IN) :: id
    CHARACTER(LEN = *), INTENT(IN) :: routine
    TYPE(t_restart_patch_description), POINTER :: RESULT

    IF(id < 1 .OR. id > SIZE(patch_data)) CALL finish(routine, "assertion failed: patch id IS OUT of range")
    RESULT => patch_data(id)%description
    IF(RESULT%id /= id) CALL finish(routine, "assertion failed: patch id does NOT match its array index")
  END FUNCTION find_patch_description

  !------------------------------------------------------------------------------------------------
  !
  !  Set global restart attributes.
  !
  SUBROUTINE set_restart_attributes (restartAttributes, restart_args)
    TYPE(t_RestartAttributeList), INTENT(INOUT) :: restartAttributes
    TYPE(t_restart_args), INTENT(IN) :: restart_args

    TYPE(t_restart_patch_description), POINTER :: p_pd
    CHARACTER(LEN=MAX_NAME_LENGTH) :: attrib_name
    INTEGER                        :: jp, jp_end, jg, i, current_jfile, effectiveDomainCount

    CHARACTER(LEN=*), PARAMETER    :: routine = modname//'set_restart_attributes'
    CHARACTER(LEN=*), PARAMETER    :: attrib_format_int  = '(a,i2.2)'
    CHARACTER(LEN=*), PARAMETER    :: attrib_format_int2 = '(a,i2.2,a,i2.2)'

    CHARACTER(len=MAX_CHAR_LENGTH) :: attname   ! attribute name

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
        p_pd => patch_data(i)%description

        ! set time levels
        jg = p_pd%id
        CALL setDynamicPatchRestartAttributes(restartAttributes, jg, p_pd%nold, p_pd%nnow, p_pd%nnew, p_pd%nnow_rcf, p_pd%nnew_rcf)

        ! additional restart-output for nonhydrostatic model
        IF (p_pd%l_opt_sim_time) THEN
            WRITE(attrib_name, attrib_format_int) 'sim_time_DOM', jg
            CALL restartAttributes%setReal (TRIM(attrib_name), p_pd%opt_sim_time)
        END IF

        !-------------------------------------------------------------
        ! DR
        ! WORKAROUND FOR FIELDS WHICH NEED TO GO INTO THE RESTART FILE,
        ! BUT SO FAR CANNOT BE HANDELED CORRECTLY BY ADD_VAR OR
        ! SET_RESTART_ATTRIBUTE
        !-------------------------------------------------------------
        IF (p_pd%l_opt_ndyn_substeps) THEN
            WRITE(attrib_name, attrib_format_int) 'ndyn_substeps_DOM', jg
            CALL restartAttributes%setInteger(TRIM(attrib_name), p_pd%opt_ndyn_substeps)
        END IF

        IF (p_pd%l_opt_jstep_adv_marchuk_order) THEN
            WRITE(attrib_name, attrib_format_int) 'jstep_adv_marchuk_order_DOM', jg
            CALL restartAttributes%setInteger(TRIM(attrib_name), p_pd%opt_jstep_adv_marchuk_order)
        END IF

        IF (ALLOCATED(p_pd%opt_t_elapsed_phy) .AND. ALLOCATED(p_pd%opt_lcall_phy)) THEN
            CALL setPhysicsRestartAttributes(restartAttributes, jg, p_pd%opt_t_elapsed_phy, p_pd%opt_lcall_phy)
        END IF
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
    CHARACTER(LEN=*), PARAMETER   :: routine = modname//'has_valid_time_level'
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
  SUBROUTINE restart_write_var_list(p_pd, restart_args)
    TYPE(t_patch_data), TARGET, INTENT(IN) :: p_pd
    TYPE(t_restart_args), INTENT(IN) :: restart_args

    TYPE(t_restart_file), POINTER   :: p_rf
    TYPE(t_var_metadata), POINTER   :: p_info
    TYPE(t_reorder_data), POINTER   :: p_ri
    TYPE(t_var_data), POINTER       :: p_vars(:)

    INTEGER                         :: iv, nval, ierrstat, nlevs, nv_off, &
      &                                np, mpi_error, i, idate, itime, status, ilev
    INTEGER(KIND=MPI_ADDRESS_KIND)  :: ioff(0:num_work_procs-1)
    REAL(dp), ALLOCATABLE           :: var1_dp(:), var2_dp(:,:), var3_dp(:)
    INTEGER                         :: ichunk, nchunks, chunk_start, chunk_end,     &
      &                                this_chunk_nlevs, ioff2

    CHARACTER(LEN=*), PARAMETER     :: routine = modname//'restart_write_var_list'
    ! For timing
    REAL(dp)                        :: t_get, t_write, t_0, mb_get, mb_wr

#ifdef DEBUG
    WRITE (nerr,FORMAT_VALS3)routine,' p_pe=',p_pe
#endif

    ! check process
    IF (.NOT. my_process_is_restart()) CALL finish(routine, NO_RESTART_PE)

    t_get   = 0.d0
    t_write = 0.d0
    mb_get  = 0.d0
    mb_wr   = 0.d0

    ! write restart time
    idate = cdiEncodeDate(restart_args%datetime%year, restart_args%datetime%month, restart_args%datetime%day)
    itime = cdiEncodeTime(restart_args%datetime%hour, restart_args%datetime%minute, NINT(restart_args%datetime%second))

    p_rf => p_pd%restart_file
    CALL taxisDefVdate(p_rf%cdiIds%taxis, idate)
    CALL taxisDefVtime(p_rf%cdiIds%taxis, itime)
    status = streamDefTimestep(p_rf%cdiIds%file, p_rf%cdiTimeIndex)

    p_rf%cdiTimeIndex = p_rf%cdiTimeIndex + 1

    ! check the contained array of restart variables
    p_vars => p_rf%var_data
    IF (.NOT. ASSOCIATED(p_vars)) RETURN

    ! get maximum number of data points in a slice and allocate tmp. variables
    nval = MAX(p_pd%commData%cells%n_glb, p_pd%commData%edges%n_glb, p_pd%commData%verts%n_glb)

    ! allocate RMA memory
    ALLOCATE(var1_dp(nval*restart_chunk_size), var2_dp(nval,restart_chunk_size), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, ALLOCATE_FAILED)

    ioff(:) = p_pd%commData%mem_win_off(:)

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
      
      ! var1 is stored in the order in which the variable was stored on compute PEs,
      ! get it back into the global storage order
      ALLOCATE(var3_dp(p_ri%n_glb), STAT=ierrstat) ! Must be allocated to exact size
      IF (ierrstat /= SUCCESS) CALL finish (routine, ALLOCATE_FAILED)

      ! no. of chunks of levels (each of size "restart_chunk_size"):
      nchunks = (nlevs-1)/restart_chunk_size + 1
      ! loop over all chunks (of levels)
      LEVELS : DO ichunk=1,nchunks
        chunk_start       = (ichunk-1)*restart_chunk_size + 1
        chunk_end         = MIN(chunk_start+restart_chunk_size-1, nlevs)
        this_chunk_nlevs  = (chunk_end - chunk_start + 1)

        ! retrieve part of variable from every worker PE using MPI_Get
        nv_off  = 0
        DO np = 0, num_work_procs-1
          IF(p_ri%pe_own(np) == 0) CYCLE
          
          ! number of words to transfer
          nval = p_ri%pe_own(np) * this_chunk_nlevs
          t_0 = p_mpi_wtime()
          CALL MPI_Win_lock(MPI_LOCK_SHARED, np, MPI_MODE_NOCHECK, mpi_win, mpi_error)
   !       CALL check_mpi_error(routine, 'MPI_Win_lock', mpi_error, .TRUE.)
          
          CALL MPI_Get(var1_dp(1), nval, p_real_dp, np, ioff(np), &
            &          nval, p_real_dp, mpi_win, mpi_error)
  !        CALL check_mpi_error(routine, 'MPI_Get', mpi_error, .TRUE.)
          
          CALL MPI_Win_unlock(np, mpi_win, mpi_error)
  !        CALL check_mpi_error(routine, 'MPI_Win_unlock', mpi_error, .TRUE.)
          
          t_get  = t_get  + p_mpi_wtime() - t_0
          mb_get = mb_get + nval

          ! update the offset in the RMA window on compute PEs
          ioff(np) = ioff(np) + INT(nval,i8)

          ! separate the levels received from PE "np":
          ioff2 = 0
          DO ilev=1,this_chunk_nlevs
            DO i=1,p_ri%pe_own(np)
              var2_dp(i+nv_off,ilev) = var1_dp(ioff2 + i)
            END DO
            ioff2 = ioff2 + p_ri%pe_own(np)
          END DO
          ! update the offset in var2
          nv_off = nv_off + p_ri%pe_own(np)
        END DO
        t_0 = p_mpi_wtime()

        ! write field content into a file
        DO ilev=chunk_start, chunk_end
!$OMP PARALLEL DO
          DO i = 1, p_ri%n_glb
            var3_dp(i) = var2_dp(p_ri%reorder_index(i),(ilev-chunk_start+1))
          ENDDO
!$OMP END PARALLEL DO
          CALL streamWriteVarSlice(p_rf%cdiIds%file, p_info%cdiVarID, (ilev-1), var3_dp(:), 0)
          mb_wr = mb_wr + REAL(SIZE(var3_dp), wp)
        END DO
        t_write = t_write + p_mpi_wtime() - t_0

      ENDDO LEVELS

#ifdef DEBUG
      WRITE (nerr,FORMAT_VALS7I)routine,' p_pe=',p_pe,' restart pe writes field=', &
        &                               TRIM(p_info%name),' data=',p_ri%n_glb*nlevs
#endif

      DEALLOCATE(var3_dp, STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish (routine, DEALLOCATE_FAILED)

    ENDDO VAR_LOOP

    DEALLOCATE(var1_dp, var2_dp, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, DEALLOCATE_FAILED)
    mb_get = mb_get*8*1.d-6
    mb_wr  = mb_wr*8*1.d-6

    IF (msg_level >= 12) THEN
      WRITE (0,'(10(a,f10.3))') &
           & ' Restart: Got ',mb_get,' MB, time get: ',t_get,' s [',mb_get/MAX(1.e-6_wp,t_get), &
           & ' MB/s], time write: ',t_write,' s [',mb_wr/MAX(1.e-6_wp,t_write),        &
           & ' MB/s]'
    ENDIF

  END SUBROUTINE restart_write_var_list

  !------------------------------------------------------------------------------------------------
  !
  ! Write restart variable lists for a compute PE.
  !
  SUBROUTINE compute_write_var_list(p_pd)
    TYPE(t_patch_data), TARGET, INTENT(IN) :: p_pd

    TYPE(t_var_metadata), POINTER   :: p_info
    TYPE(t_reorder_data), POINTER   :: p_ri
    TYPE(t_var_data), POINTER       :: p_vars(:)
    REAL(wp), POINTER               :: r_ptr(:,:,:)
    INTEGER                         :: iv, mpi_error, nindex, ierrstat, nlevs, i, jk, &
      &                                var_ref_pos
    INTEGER(i8)                     :: ioff
    CHARACTER(LEN=*), PARAMETER     :: routine = modname//'compute_write_var_list'

#ifdef DEBUG
    WRITE (nerr,FORMAT_VALS3)routine,' p_pe=',p_pe
#endif

    ! check process
    IF (.NOT. my_process_is_work()) CALL finish(routine, NO_COMPUTE_PE)

    ! check the array of restart variables
    p_vars => p_pd%restart_file%var_data
    IF (.NOT. ASSOCIATED(p_vars)) RETURN

    ! offset in RMA window for async restart
    ioff = p_pd%commData%my_mem_win_off

    ! in case of async restart: Lock own window before writing to it
    CALL MPI_Win_lock(MPI_LOCK_EXCLUSIVE, p_pe_work, MPI_MODE_NOCHECK, mpi_win, mpi_error)
    CALL check_mpi_error(routine, 'MPI_Win_lock', mpi_error, .TRUE.)

    ! go over the all restart variables in the associated array
    DO iv = 1, SIZE(p_vars)

      ! get pointer to metadata
      p_info => p_vars(iv)%info

#ifdef DEBUG
      WRITE (nerr,FORMAT_VALS5I)routine,' p_pe=',p_pe,' compute pe processes field=',TRIM(p_info%name)
#endif

      ! check time level of the field
      IF (.NOT. has_valid_time_level(p_info, p_pd%description)) CYCLE

      ! Check if first dimension of array is nproma.
      ! Otherwise we got an array which is not suitable for this output scheme.
      IF (p_info%used_dimensions(1) /= nproma) &
        CALL finish(routine,'1st dim is not nproma: '//TRIM(p_info%name))

      ! init. data pointer
      r_ptr => NULL()

      ! get data index
      IF (p_info%lcontained) THEN
        nindex = p_info%ncontained
      ELSE
        nindex = 1
      ENDIF

      ! get data pointer
      SELECT CASE (p_info%ndims)
        CASE (1)
          CALL message(routine, p_info%name)
          CALL finish(routine,'1d arrays not handled yet.')
        CASE (2)
          ! make a 3D copy of the array
          ALLOCATE(r_ptr(p_info%used_dimensions(1),1,p_info%used_dimensions(2)), &
            &      STAT=ierrstat)
          IF (ierrstat /= SUCCESS) CALL finish(routine, ALLOCATE_FAILED)
          var_ref_pos = 3
          IF (p_info%lcontained)  var_ref_pos = p_info%var_ref_pos
          SELECT CASE(var_ref_pos)
          CASE (1)
            r_ptr(:,1,:) = p_vars(iv)%r_ptr(nindex,:,:,1,1)
          CASE (2)
            r_ptr(:,1,:) = p_vars(iv)%r_ptr(:,nindex,:,1,1)
          CASE (3)
            r_ptr(:,1,:) = p_vars(iv)%r_ptr(:,:,nindex,1,1)
          CASE default
            CALL finish(routine, "internal error!")
          END SELECT
        CASE (3)
          ! copy the pointer
          var_ref_pos = 4
          IF (p_info%lcontained)  var_ref_pos = p_info%var_ref_pos
          SELECT CASE(var_ref_pos)
          CASE (1)
            r_ptr => p_vars(iv)%r_ptr(nindex,:,:,:,1)
          CASE (2)
            r_ptr => p_vars(iv)%r_ptr(:,nindex,:,:,1)
          CASE (3)
            r_ptr => p_vars(iv)%r_ptr(:,:,nindex,:,1)
          CASE (4)
            r_ptr => p_vars(iv)%r_ptr(:,:,:,nindex,1)
          CASE default
            CALL finish(routine, "internal error!")
          END SELECT
        CASE (4)
          CALL message(routine, p_info%name)
          CALL finish(routine,'4d arrays not handled yet.')
        CASE (5)
          CALL message(routine, p_info%name)
          CALL finish(routine,'5d arrays not handled yet.')
        CASE DEFAULT
          CALL message(routine, p_info%name)
          CALL finish(routine,'dimension not set.')
      END SELECT

      ! get number of data levels
      IF(p_info%ndims == 2) THEN
        nlevs = 1
      ELSE
        nlevs = p_info%used_dimensions(2)
      ENDIF

      ! get pointer to reorder data
      p_ri => get_reorder_ptr(p_pd%commData, p_info, routine)

#ifdef DEBUG
      WRITE (nerr,FORMAT_VALS7I)routine,' p_pe=',p_pe,' compute pe writes field=', &
        &                               TRIM(p_info%name),' data=',nlevs*p_ri%n_own
#endif

      ! just copy the OWN data points to the memory window
      DO jk = 1, nlevs
        DO i = 1, p_ri%n_own
          mem_ptr_dp(ioff+INT(i,i8)) = REAL(r_ptr(p_ri%own_idx(i),jk,p_ri%own_blk(i)),dp)
        ENDDO
        ioff = ioff + INT(p_ri%n_own,i8)
      END DO

      ! deallocate temp. 2D array
      IF(p_info%ndims == 2) DEALLOCATE(r_ptr)

    ENDDO

    ! unlock RMA window
    CALL MPI_Win_unlock(p_pe_work, mpi_win, mpi_error)
    CALL check_mpi_error(routine, 'MPI_Win_unlock', mpi_error, .TRUE.)

  END SUBROUTINE compute_write_var_list

  !------------------------------------------------------------------------------------------------
  !
  ! Initialize the variable list of the given restart file.
  !
  SUBROUTINE init_restart_variables(p_rf,patch_id)
    TYPE(t_restart_file), TARGET, INTENT(IN) :: p_rf
    INTEGER, INTENT(IN)           :: patch_id

    TYPE(t_var_data), POINTER     :: p_vars(:)
    TYPE(t_var_metadata), POINTER :: p_info
    CHARACTER(LEN=*), PARAMETER   :: routine = modname//'init_restart_variables'
    INTEGER                       :: iv

#ifdef DEBUG
    WRITE (nerr,FORMAT_VALS3)routine,' p_pe=',p_pe
#endif

    ! check the contained array of restart variables
    p_vars => p_rf%var_data
    IF (.NOT. ASSOCIATED(p_vars)) RETURN

    ! go over the all restart variables in the associated array
    DO iv = 1, SIZE(p_vars)

      ! get pointer to metadata
      p_info => p_vars(iv)%info

      ! check time level of the field
      IF (.NOT. has_valid_time_level(p_info, find_patch_description(patch_id, routine))) CYCLE

      ! define the CDI variable
      !XXX: The code I found here simply assumed that all variables are of TYPE REAL. I have NOT changed this behavior, however, the .FALSE. constants should be replaced by something more sensible IN the future.
      CALL p_rf%cdiIds%defineVariable(p_info, .FALSE., .FALSE.)
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
    CHARACTER(LEN=32)             :: datetime
    INTEGER                       :: restart_type, i
    CHARACTER(LEN=*), PARAMETER   :: routine = modname//'open_restart_file'
    TYPE (t_keyword_list), POINTER :: keywords => NULL()

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

    datetime = iso8601(restart_args%datetime)

    ! build the file name
    CALL associate_keyword("<gridfile>",   TRIM(get_filename_noext(description%base_filename)),   keywords)
    CALL associate_keyword("<idom>",       TRIM(int2string(description%id, "(i2.2)")),            keywords)
    CALL associate_keyword("<rsttime>",    TRIM(datetime),                                 keywords)
    CALL associate_keyword("<mtype>",      TRIM(p_rf%model_type),                          keywords)
    ! replace keywords in file name
    p_rf%filename = TRIM(with_keywords(keywords, TRIM(restart_filename)))

    commData => p_pd%commData
    IF(ALLOCATED(description%opt_pvct)) THEN
        CALL p_rf%cdiIds%openRestartAndCreateIds(TRIM(p_rf%filename), restart_type, restartAttributes, commData%cells%n_glb, &
                                                &commData%verts%n_glb, commData%edges%n_glb, description%cell_type, &
                                                &description%v_grid_defs(1:description%v_grid_count), description%opt_pvct)
    ELSE
        CALL p_rf%cdiIds%openRestartAndCreateIds(TRIM(p_rf%filename), restart_type, restartAttributes, commData%cells%n_glb, &
                                                &commData%verts%n_glb, commData%edges%n_glb, description%cell_type, &
                                                &description%v_grid_defs(1:description%v_grid_count))
    END IF

#ifdef DEBUG
    WRITE (nerr, FORMAT_VALS5)routine,' p_pe=',p_pe,' open netCDF file with ID=',p_rf%cdiIds%file
#endif

    ! set cdi internal time index to 0 for writing time slices in netCDF
    p_rf%cdiTimeIndex = 0

    ! init list of restart variables
    CALL init_restart_variables(p_rf, p_pd%description%id)

    CALL streamDefVlist(p_rf%cdiIds%file, p_rf%cdiIds%vlist)

  END SUBROUTINE open_restart_file

  !------------------------------------------------------------------------------------------------
  !
  ! Closes the given restart file.
  !
  SUBROUTINE close_restart_file(rf)

    TYPE (t_restart_file), INTENT(INOUT)  :: rf
    CHARACTER(LEN=*), PARAMETER           :: routine = modname//'close_restart_file'

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

END MODULE mo_io_restart_async


