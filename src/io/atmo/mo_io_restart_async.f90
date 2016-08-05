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

  USE mo_util_file,               ONLY: util_symlink, util_unlink, util_islink
  USE mo_exception,               ONLY: finish, message, message_text, get_filename_noext
  USE mo_fortran_tools,           ONLY: assign_if_present, assign_if_present_allocatable
  USE mo_kind,                    ONLY: wp, i8, dp
  USE mo_datetime,                ONLY: t_datetime, iso8601, iso8601extended
  USE mo_io_units,                ONLY: nerr, filename_max
  USE mo_var_list,                ONLY: nvar_lists, var_lists, new_var_list, delete_var_lists
  USE mo_linked_list,             ONLY: t_list_element, t_var_list
  USE mo_io_restart_attributes,   ONLY: t_RestartAttributeList, RestartAttributeList_make
  USE mo_dynamics_config,         ONLY: nold, nnow, nnew, nnew_rcf, nnow_rcf, iequations
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
  USE mo_model_domain,            ONLY: p_patch
  USE mo_cdi,                     ONLY: CDI_UNDEFID, FILETYPE_NC2, FILETYPE_NC4, CDI_GLOBAL, DATATYPE_FLT64, &
                                      & TAXIS_ABSOLUTE, ZAXIS_DEPTH_BELOW_SEA, ZAXIS_GENERIC, ZAXIS_HEIGHT, ZAXIS_HYBRID, &
                                      & ZAXIS_HYBRID_HALF, ZAXIS_LAKE_BOTTOM, ZAXIS_MIX_LAYER, ZAXIS_SEDIMENT_BOTTOM_TW, &
                                      & ZAXIS_SURFACE, ZAXIS_TOA, TIME_VARIABLE, ZAXIS_DEPTH_BELOW_LAND, GRID_UNSTRUCTURED, &
                                      & vlistDefVar, cdiEncodeDate, cdiEncodeTime, streamDefTimestep, gridDestroy, &
                                      & streamWriteVarSlice, streamDefVlist, vlistDefVarDatatype, vlistDefVarName, &
                                      & vlistDefVarLongname, vlistDefVarUnits, vlistDefVarMissval, taxisDefVdate, taxisDefVtime
  USE mo_util_cdi,                ONLY: cdiGetStringError
  USE mo_cdi_constants,           ONLY: ZA_SURFACE, ZA_HYBRID, ZA_HYBRID_HALF, ZA_DEPTH_BELOW_LAND, ZA_DEPTH_BELOW_LAND_P1, &
                                      & ZA_SNOW, ZA_SNOW_HALF, ZA_HEIGHT_2M, ZA_HEIGHT_10M, ZA_TOA, ZA_LAKE_BOTTOM, ZA_MIX_LAYER, &
                                      & ZA_LAKE_BOTTOM_HALF, ZA_SEDIMENT_BOTTOM_TW_HALF, ZA_DEPTH_BELOW_SEA, &
                                      & ZA_DEPTH_BELOW_SEA_HALF, ZA_GENERIC_ICE, ZA_DEPTH_RUNOFF_S, ZA_DEPTH_RUNOFF_G, ZA_COUNT, &
                                      & GRID_UNSTRUCTURED_EDGE, GRID_UNSTRUCTURED_VERT, GRID_UNSTRUCTURED_CELL, cdi_zaxis_types, &
                                      & GRID_UNSTRUCTURED_COUNT
  USE mo_cf_convention
  USE mo_packed_message,          ONLY: t_PackedMessage
  USE mo_util_string,             ONLY: t_keyword_list, associate_keyword, with_keywords, &
    &                                   int2string, toCharacter
  USE mo_util_restart,            ONLY: t_v_grid, t_restart_cdi_ids, set_vertical_grid, closeAndDestroyIds, &
                                      & openRestartAndCreateIds, defineVariable, setGeneralRestartAttributes, &
                                      & setDynamicPatchRestartAttributes, setPhysicsRestartAttributes

#ifndef NOMPI
  USE mo_mpi,                     ONLY: p_pe, p_pe_work, p_restart_pe0, p_comm_work, p_work_pe0, num_work_procs, MPI_SUCCESS, &
                                      & stop_mpi, p_send, p_recv, p_barrier, p_bcast, my_process_is_restart, my_process_is_work, &
                                      & p_comm_work_2_restart, p_n_work, p_int, process_mpi_restart_size, p_int_i8, p_real_dp, &
                                      & p_comm_work_restart, p_mpi_wtime, p_send_packed, p_recv_packed, p_bcast_packed, &
                                      & p_pack_int, p_pack_bool, p_pack_real, p_unpack_int, p_unpack_bool, p_unpack_real, &
                                      & p_int_byte, get_my_mpi_work_id, p_comm_rank, p_unpack_allocatable_real, &
                                      & p_unpack_allocatable_logical, p_pack_allocatable_real, p_pack_allocatable_logical, &
                                      & p_pack_allocatable_int, p_unpack_allocatable_int, process_mpi_all_comm

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

  ! maximum no. of output_nml output files
  INTEGER, parameter :: MAX_NML_OUTPUT_FILES      = 100

  ! minimum number of dynamic restart arguments
  INTEGER, PARAMETER :: MIN_DYN_RESTART_ARGS      = 16 + MAX_NML_OUTPUT_FILES

  ! minimum number of dynamic restart patch data
  ! id, l_dom_active, time levels and optional attributes
  INTEGER, PARAMETER :: MIN_DYN_RESTART_PDATA     = 23

  ! maximumm number of verticale axes
  INTEGER, PARAMETER :: MAX_VERTICAL_AXES         = 19

  ! common constant strings
  CHARACTER(LEN=*), PARAMETER :: modname                  = 'shared/mo_io_restart_async/'
  CHARACTER(LEN=*), PARAMETER :: ALLOCATE_FAILED          = 'ALLOCATE failed!'
  CHARACTER(LEN=*), PARAMETER :: DEALLOCATE_FAILED        = 'DEALLOCATE failed!'
  CHARACTER(LEN=*), PARAMETER :: UNKNOWN_GRID_TYPE        = 'Unknown grid type!'
  CHARACTER(LEN=*), PARAMETER :: UNKNOWN_FILE_FORMAT      = 'Unknown file format for restart file.'
  CHARACTER(LEN=*), PARAMETER :: UNKNOWN_VERT_GRID_DESCR  = 'Vertical grid description not found.'
  CHARACTER(LEN=*), PARAMETER :: ASYNC_RESTART_REQ_MPI    = 'asynchronous restart can only run with MPI!'
  CHARACTER(LEN=*), PARAMETER :: NO_COMPUTE_PE            = 'Must be called on a compute PE!'
  CHARACTER(LEN=*), PARAMETER :: NO_RESTART_PE            = 'Must be called on a restart PE!'
  CHARACTER(LEN=*), PARAMETER :: NET_CDF_ERROR_FORMAT     = '(a,i5,a)'
  CHARACTER(LEN=*), PARAMETER :: WRONG_ARRAY_SIZE         =  &
    & 'No or wrong array size set in prepare_async_restart() for='

  CHARACTER(LEN=*), PARAMETER :: FORMAT_VALS3             = '(a,a,i3)'
  CHARACTER(LEN=*), PARAMETER :: FORMAT_VALS5             = '(a,a,i3,a,i3)'
  CHARACTER(LEN=*), PARAMETER :: FORMAT_VALS5I            = '(a,a,i3,a,a)'
  CHARACTER(LEN=*), PARAMETER :: FORMAT_VALS7             = '(a,a,i3,a,i6,a,i3)'
  CHARACTER(LEN=*), PARAMETER :: FORMAT_VALS7I            = '(a,a,i3,a,a,a,i8)'
  CHARACTER(LEN=*), PARAMETER :: FORMAT_VALS9             = '(a,a,i3,a,i6,a,i3,a,i3)'

  !------------------------------------------------------------------------------------------------
  ! TYPE t_var_data (restart variable)
  !
  TYPE t_var_data
    REAL(wp), POINTER     :: r_ptr(:,:,:,:,:)
    TYPE(t_var_metadata)  :: info
  END TYPE t_var_data

  !------------------------------------------------------------------------------------------------
  ! TYPE t_restart_file
  !
  TYPE t_restart_file
    ! the following data can be set before opening the restart file
    TYPE(t_var_data), POINTER   :: var_data(:)

    !-----------------------------------
    ! used for remote memory access
    INTEGER(i8)                 :: my_mem_win_off
    INTEGER(i8), ALLOCATABLE    :: mem_win_off(:)
    !-----------------------------------

    ! the following members are set during open
    CHARACTER(LEN=filename_max) :: filename
    CHARACTER(LEN=32)           :: model_type
    CHARACTER(len=64)           :: linkname
    CHARACTER(len=10)           :: linkprefix
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

  !------------------------------------------------------------------------------------------------
  ! TYPE t_patch_data contains the reordering data for cells, edges and verts
  !
  TYPE t_patch_data
    ! reorder data
    TYPE(t_reorder_data) :: cells
    TYPE(t_reorder_data) :: edges
    TYPE(t_reorder_data) :: verts

    ! vertical grid definitions
    TYPE(t_v_grid), POINTER :: v_grid_defs(:)
    INTEGER :: v_grid_count

    ! restart file data
    TYPE(t_restart_file) :: restart_file

    ! logical patch id
    INTEGER :: id

    ! current model domain activity flag
    LOGICAL :: l_dom_active

    ! number of full levels
    INTEGER :: nlev

    ! cell type
    INTEGER :: cell_type

    ! total # of cells, # of vertices per cell
    INTEGER :: n_patch_cells_g
    ! total # of cells, shape of control volume for edge
    INTEGER :: n_patch_edges_g
    ! total # of vertices, # of vertices per dual cell
    INTEGER :: n_patch_verts_g

    ! global number of blocks
    INTEGER :: nblks_glb_c, nblks_glb_v, nblks_glb_e

    ! process id
    INTEGER :: restart_proc_id

    ! id of PE0 of working group (/= 0 in case of processor splitting)
    INTEGER :: work_pe0_id

    ! base file name contains already logical patch ident
    CHARACTER(LEN=filename_max) :: base_filename

    ! dynamic patch arguments (mandatory)
    INTEGER :: nold,nnow,nnew,nnew_rcf,nnow_rcf

    ! dynamic patch arguments (optionally)
    LOGICAL               :: l_opt_depth
    INTEGER               :: opt_depth
    LOGICAL               :: l_opt_depth_lnd
    INTEGER               :: opt_depth_lnd
    LOGICAL               :: l_opt_nlev_snow
    INTEGER               :: opt_nlev_snow
    LOGICAL               :: l_opt_nice_class
    INTEGER               :: opt_nice_class
    LOGICAL               :: l_opt_ndyn_substeps
    INTEGER               :: opt_ndyn_substeps
    LOGICAL               :: l_opt_jstep_adv_marchuk_order
    INTEGER               :: opt_jstep_adv_marchuk_order
    LOGICAL               :: l_opt_sim_time
    REAL(wp)              :: opt_sim_time
    LOGICAL               :: l_opt_ndom
    INTEGER               :: opt_ndom
    !
    REAL(wp), ALLOCATABLE :: opt_pvct(:)
    LOGICAL, ALLOCATABLE  :: opt_lcall_phy(:)
    REAL(wp), ALLOCATABLE :: opt_t_elapsed_phy(:)
  END TYPE t_patch_data
  TYPE(t_patch_data), ALLOCATABLE, TARGET :: patch_data (:)

  !------------------------------------------------------------------------------------------------
  ! patch independent arguments
  !
  TYPE t_restart_args
    TYPE(t_datetime)  :: datetime
    INTEGER           :: jstep

    INTEGER, ALLOCATABLE  :: opt_output_jfile(:)
  END TYPE t_restart_args
  TYPE(t_restart_args), TARGET :: restart_args

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
  !> Set dynamically data for asynchronous restart.
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
                                &   opt_ndom,                     &
                                &   opt_output_jfile )

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
    INTEGER,              INTENT(IN), OPTIONAL :: opt_output_jfile(:)
    INTEGER,              INTENT(IN), OPTIONAL :: opt_ndom            !< no. of domains (appended to symlink name)

    TYPE(t_patch_data),   POINTER  :: p_pd
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
    p_pd => find_patch(patch_id, routine)

    ! set activity flag - this needs to be done on all compute PEs because the restart
    ! file may be incomplete otherwise when a nest is started during runtime
    p_pd%l_dom_active = l_dom_active

    ! otherwise, only the first compute PE needs the dynamic restart arguments
    ! in the case of processor splitting, the first PE of the split subset needs them as well
    IF (p_pe_work == 0 .OR. p_pe_work == p_pd%work_pe0_id) THEN

      ! Patch-independent attributes (only communicated through patch 1)
      IF (patch_id == 1) CALL assign_if_present_allocatable(restart_args%opt_output_jfile, opt_output_jfile)

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
  !> Writes all restart data into one or more files (one file per patch, collective call).
  !
  SUBROUTINE write_async_restart (datetime, jstep)

    TYPE(t_datetime), INTENT(IN)    :: datetime
    INTEGER,          INTENT(IN)    :: jstep

    TYPE(t_patch_data), POINTER     :: p_pd
    INTEGER                         :: idx
    TYPE(t_RestartAttributeList), POINTER :: restartAttributes

    CHARACTER(LEN=*), PARAMETER :: routine = modname//'write_async_restart'

#ifdef NOMPI
    CALL finish (routine, ASYNC_RESTART_REQ_MPI)
#else

#ifdef DEBUG
    WRITE (nerr,FORMAT_VALS3)routine,' p_pe=',p_pe
#endif

    ! check kind of process
    IF (.NOT. my_process_is_work() .AND. .NOT. my_process_is_restart()) RETURN

    IF (my_process_is_work()) THEN
      CALL compute_wait_for_restart
      CALL compute_start_restart (datetime, jstep)
    END IF

    ! do the restart output
    DO idx = 1, SIZE(patch_data)

      p_pd => patch_data(idx)

      ! check if the patch is actice
      IF (.NOT. p_pd%l_dom_active) CYCLE

      ! write the variable restart lists
      IF (my_process_is_restart()) THEN

        ! consider the right restart process
        IF (p_pe == p_pd%restart_proc_id) THEN

          ! set global restart attributes/lists
          restartAttributes => RestartAttributeList_make()
          CALL set_restart_attributes(restartAttributes)
          CALL defineVerticalGrids(p_pd)

#ifdef DEBUG
          CALL print_restart_arguments()
          CALL restartAttributes%printAttributes()
          CALL print_restart_name_lists()
#endif
          CALL open_restart_file(p_pd, restartAttributes)

          CALL restartAttributes%destruct()
          DEALLOCATE(restartAttributes)

          ! collective call to write the restart variables
          CALL restart_write_var_list(p_pd)
          CALL create_restart_file_link(p_pd%restart_file, p_pd%restart_proc_id, p_pd%id, &
            &                           p_pd%l_opt_ndom, p_pd%opt_ndom )
          CALL close_restart_file(p_pd%restart_file)

        ENDIF

      ELSE
        ! collective call to write the restart variables
        CALL compute_write_var_list(p_pd)
      ENDIF

    ENDDO

#endif

  END SUBROUTINE write_async_restart

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

    LOGICAL                       :: done

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
      CALL restart_wait_for_start(done)

      IF(done) EXIT ! leave loop, we are done

      ! read and write restart variable lists (collective call)
      CALL write_async_restart (restart_args%datetime,    &
                              & restart_args%jstep)

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

  SUBROUTINE packRestartMetadata(datetime, jstep, message)
    TYPE(t_datetime), INTENT(IN) :: datetime
    INTEGER, VALUE :: jstep
    TYPE(t_PackedMessage), INTENT(INOUT) :: message

    INTEGER :: i
    TYPE(t_patch_data), POINTER :: curPatch

    ! set patch independent arguments
    CALL message%pack(datetime%year)
    CALL message%pack(datetime%month)
    CALL message%pack(datetime%day)
    CALL message%pack(datetime%hour)
    CALL message%pack(datetime%minute)
    CALL message%pack(datetime%second)
    CALL message%pack(datetime%caltime)
    CALL message%pack(INT(datetime%calday))
    CALL message%pack(datetime%daysec)
    CALL message%pack(jstep)
    CALL message%pack(restart_args%opt_output_jfile)

    ! set data of all patches
    DO i = 1, SIZE(patch_data)
        curPatch => patch_data(i)

        ! patch id
        CALL message%pack(curPatch%id)

        ! activity flag
        CALL message%pack(curPatch%l_dom_active)

        ! time levels
        CALL message%pack(nold(curPatch%id))
        CALL message%pack(nnow(curPatch%id))
        CALL message%pack(nnew(curPatch%id))
        CALL message%pack(nnew_rcf(curPatch%id))
        CALL message%pack(nnow_rcf(curPatch%id))

        ! optional parameter values
        CALL message%pack(curPatch%l_opt_depth)
        CALL message%pack(curPatch%opt_depth)
        CALL message%pack(curPatch%l_opt_depth_lnd)
        CALL message%pack(curPatch%opt_depth_lnd)
        CALL message%pack(curPatch%l_opt_nlev_snow)
        CALL message%pack(curPatch%opt_nlev_snow)
        CALL message%pack(curPatch%l_opt_nice_class)
        CALL message%pack(curPatch%opt_nice_class)
        CALL message%pack(curPatch%l_opt_ndyn_substeps)
        CALL message%pack(curPatch%opt_ndyn_substeps)
        CALL message%pack(curPatch%l_opt_jstep_adv_marchuk_order)
        CALL message%pack(curPatch%opt_jstep_adv_marchuk_order)
        CALL message%pack(curPatch%l_opt_sim_time)
        CALL message%pack(curPatch%opt_sim_time)
        CALL message%pack(curPatch%l_opt_ndom)
        CALL message%pack(curPatch%opt_ndom)

        ! optional parameter arrays
        CALL message%pack(curPatch%opt_pvct)
        CALL message%pack(curPatch%opt_lcall_phy)
        CALL message%pack(curPatch%opt_t_elapsed_phy)
    END DO
  END SUBROUTINE packRestartMetadata

  SUBROUTINE unpackRestartMetadata(datetime, jstep, message)
    TYPE(t_datetime), INTENT(INOUT) :: datetime
    INTEGER, VALUE :: jstep
    TYPE(t_PackedMessage), INTENT(INOUT) :: message

    INTEGER :: i, calday, this_patch
    TYPE(t_patch_data), POINTER :: curPatch
    CHARACTER(LEN = *), PARAMETER :: routine = modname//":unpackRestartMetadata"

    ! get patch independent arguments
    CALL message%unpack(datetime%year)
    CALL message%unpack(datetime%month)
    CALL message%unpack(datetime%day)
    CALL message%unpack(datetime%hour)
    CALL message%unpack(datetime%minute)
    CALL message%unpack(datetime%second)
    CALL message%unpack(datetime%caltime)
    CALL message%unpack(calday)
    datetime%calday = INT(calday,i8)
    CALL message%unpack(datetime%daysec)
    CALL message%unpack(jstep)
    CALL message%unpack(restart_args%opt_output_jfile)

    ! get patch dependent arguments
    DO i = 1, SIZE(patch_data)
        ! find the patch of the current patch id
        CALL message%unpack(this_patch)
        curPatch => find_patch(this_patch, routine)

        ! activity flag
        CALL message%unpack(curPatch%l_dom_active)

        ! time levels
        CALL message%unpack(nold(this_patch))
        CALL message%unpack(nnow(this_patch))
        CALL message%unpack(nnew(this_patch))
        CALL message%unpack(nnew_rcf(this_patch))
        CALL message%unpack(nnow_rcf(this_patch))

        ! optional parameter values
        CALL message%unpack(curPatch%l_opt_depth)
        CALL message%unpack(curPatch%opt_depth)
        CALL message%unpack(curPatch%l_opt_depth_lnd)
        CALL message%unpack(curPatch%opt_depth_lnd)
        CALL message%unpack(curPatch%l_opt_nlev_snow)
        CALL message%unpack(curPatch%opt_nlev_snow)
        CALL message%unpack(curPatch%l_opt_nice_class)
        CALL message%unpack(curPatch%opt_nice_class)
        CALL message%unpack(curPatch%l_opt_ndyn_substeps)
        CALL message%unpack(curPatch%opt_ndyn_substeps)
        CALL message%unpack(curPatch%l_opt_jstep_adv_marchuk_order)
        CALL message%unpack(curPatch%opt_jstep_adv_marchuk_order)
        CALL message%unpack(curPatch%l_opt_sim_time)
        CALL message%unpack(curPatch%opt_sim_time)
        CALL message%unpack(curPatch%l_opt_ndom)
        CALL message%unpack(curPatch%opt_ndom)

        ! optional parameter arrays
        CALL message%unpack(curPatch%opt_pvct)
        CALL message%unpack(curPatch%opt_lcall_phy)
        CALL message%unpack(curPatch%opt_t_elapsed_phy)
    END DO
  END SUBROUTINE unpackRestartMetadata

  !-------------------------------------------------------------------------------------------------
  !>
  !! restart_wait_for_start: Wait for a message from compute PEs that we should start restart or finish.
  !! The counterpart on the compute side is compute_start_restart/compute_shutdown_restart.
  !
  SUBROUTINE restart_wait_for_start(done)

    LOGICAL, INTENT(OUT)           :: done ! flag if we should shut down

    INTEGER :: i, iheader
    TYPE(t_patch_data), POINTER :: curPatch
    TYPE(t_PackedMessage) :: message
    CHARACTER(LEN=*), PARAMETER :: routine = modname//'restart_wait_for_start'

#ifdef DEBUG
    WRITE (nerr,FORMAT_VALS3)routine,' is called, p_pe=',p_pe
#endif

    CALL message%construct() ! create message array

    ! receive message that we may start restart (or should finish)
    IF(p_pe_work == 0) CALL message%recv(p_work_pe0, 0, process_mpi_all_comm)
    CALL message%bcast(0, p_comm_work)

    ! set output parameter to default value
    done = .FALSE.

    ! unpack AND interpret the message
    CALL message%unpack(iheader)
    SELECT CASE(iheader)
      CASE(MSG_RESTART_START)
        CALL unpackRestartMetadata(restart_args%datetime, restart_args%jstep, message)

        ! update the patch_data
        DO i = 1, SIZE(patch_data)
            curPatch => patch_data(i)
            curPatch%nold = nold(curPatch%id)
            curPatch%nnow = nnow(curPatch%id)
            curPatch%nnew = nnew(curPatch%id)
            curPatch%nnow_rcf = nnow_rcf(curPatch%id)
            curPatch%nnew_rcf = nnew_rcf(curPatch%id)
        END DO

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
  SUBROUTINE compute_start_restart(datetime, jstep)

    TYPE(t_datetime), INTENT(IN)  :: datetime
    INTEGER,          INTENT(IN)  :: jstep

    TYPE(t_patch_data),   POINTER  :: p_pd
    CHARACTER, POINTER             :: p_msg(:)
    INTEGER                        :: i, position, messageSize, error
    TYPE(t_PackedMessage) :: message
    CHARACTER(LEN=*), PARAMETER :: routine = modname//'compute_start_restart'

#ifdef DEBUG
    WRITE (nerr,FORMAT_VALS5)routine,' p_pe=',p_pe, ' call p_barrier with communicator=',p_comm_work
#endif

    ! make sure all are here
    CALL p_barrier(comm=p_comm_work)

    ! if processor splitting is applied, the time-dependent data need to be transferred
    ! from the subset master PE to PE0, from where they are communicated to the output PE(s)
    DO i = 2, SIZE(patch_data)
      p_pd => patch_data(i)
      IF (p_pd%work_pe0_id /= 0) THEN
        IF (p_pe_work == 0) THEN
          ! recieve the package for this patch
          CALL p_recv(messageSize, p_pd%work_pe0_id, 0, comm = process_mpi_all_comm)
          ALLOCATE(p_msg(messageSize), STAT = error)
          IF(error /= SUCCESS) CALL finish(routine, "memory allocation failed")
          CALL p_recv_packed(p_msg, p_pd%work_pe0_id, 0, messageSize, comm = process_mpi_all_comm)
          position = 0

          CALL p_unpack_bool               (p_msg, position, p_pd%l_dom_active,                process_mpi_all_comm)
          CALL p_unpack_int                (p_msg, position, nnow(p_pd%id),                    process_mpi_all_comm)
          CALL p_unpack_int                (p_msg, position, nnew(p_pd%id),                    process_mpi_all_comm)
          CALL p_unpack_int                (p_msg, position, nnow_rcf(p_pd%id),                process_mpi_all_comm)
          CALL p_unpack_int                (p_msg, position, nnew_rcf(p_pd%id),                process_mpi_all_comm)

          CALL p_unpack_int                (p_msg, position, p_pd%opt_ndyn_substeps,           process_mpi_all_comm)
          CALL p_unpack_int                (p_msg, position, p_pd%opt_jstep_adv_marchuk_order, process_mpi_all_comm)
          CALL p_unpack_real               (p_msg, position, p_pd%opt_sim_time,                process_mpi_all_comm)

          CALL p_unpack_allocatable_logical(p_msg, position, p_pd%opt_lcall_phy,               process_mpi_all_comm)
          CALL p_unpack_allocatable_real   (p_msg, position, p_pd%opt_t_elapsed_phy,           process_mpi_all_comm)

          DEALLOCATE(p_msg)
        ELSE IF (p_pe_work == p_pd%work_pe0_id) THEN
          !create package of the DATA that we need to send to process 0
          p_msg => get_message_array(routine)
          position     = 0

          CALL p_pack_bool               (p_pd%l_dom_active,                p_msg, position, process_mpi_all_comm)
          CALL p_pack_int                (nnow(p_pd%id),                    p_msg, position, process_mpi_all_comm)
          CALL p_pack_int                (nnew(p_pd%id),                    p_msg, position, process_mpi_all_comm)
          CALL p_pack_int                (nnow_rcf(p_pd%id),                p_msg, position, process_mpi_all_comm)
          CALL p_pack_int                (nnew_rcf(p_pd%id),                p_msg, position, process_mpi_all_comm)

          CALL p_pack_int                (p_pd%opt_ndyn_substeps,           p_msg, position, process_mpi_all_comm)
          CALL p_pack_int                (p_pd%opt_jstep_adv_marchuk_order, p_msg, position, process_mpi_all_comm)
          CALL p_pack_real               (p_pd%opt_sim_time,                p_msg, position, process_mpi_all_comm)

          CALL p_pack_allocatable_logical(p_pd%opt_lcall_phy,               p_msg, position, process_mpi_all_comm)
          CALL p_pack_allocatable_real   (p_pd%opt_t_elapsed_phy,           p_msg, position, process_mpi_all_comm)

          !send the package
          CALL p_send(position, 0, 0, comm = process_mpi_all_comm)
          CALL p_send_packed(p_msg, 0, 0, position, comm = process_mpi_all_comm)

          DEALLOCATE(p_msg)
        END IF
      END IF
    END DO


    IF(p_pe_work == 0) THEN
      CALL message%construct()   ! create message array

      CALL message%pack(MSG_RESTART_START)  ! set command id
      CALL packRestartMetadata(datetime, jstep, message)    ! all the other DATA

      CALL message%send(p_restart_pe0, 0, process_mpi_all_comm)

      CALL message%destruct()
    ENDIF

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
  !  Get a message array to transfer all dynamical restart arguments
  !  between compute and restart PEs.
  !
  FUNCTION get_message_array (routine) 
    CHARACTER, POINTER            :: get_message_array(:)
    CHARACTER(LEN=*), INTENT(in)  :: routine

    INTEGER                       :: ierrstat, n_msg
    TYPE(t_patch_data), POINTER   :: p_pd

    ! set minimum size to transfer patch data
    n_msg = MIN_DYN_RESTART_PDATA

    ! considerate dynamic attributes
    p_pd => patch_data(1)
    n_msg = n_msg + SIZE(p_pd%opt_pvct)
    n_msg = n_msg + SIZE(p_pd%opt_lcall_phy)
    n_msg = n_msg + SIZE(p_pd%opt_t_elapsed_phy)

    ! calculate summary of all data
    n_msg = MIN_DYN_RESTART_ARGS + (SIZE(patch_data) * n_msg)

#ifdef DEBUG
    WRITE (nerr,FORMAT_VALS5)routine,' p_pe=',p_pe, &
      & ' calculated message size=',n_msg
#endif

    ! allocate memory (with 2x safety margin ;) )
    ALLOCATE (get_message_array(2*n_msg*p_int_byte), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, ALLOCATE_FAILED)

  END FUNCTION get_message_array


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

    INTEGER                             :: i

    CHARACTER(LEN=*), PARAMETER   :: routine = modname//'release_restart_file'

#ifdef DEBUG
    WRITE (nerr,FORMAT_VALS3)routine,' is called for p_pe=',p_pe
#endif

    IF (ALLOCATED(rf%mem_win_off)) DEALLOCATE(rf%mem_win_off)
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
        IF (ASSOCIATED(p_pd%v_grid_defs))      DEALLOCATE(p_pd%v_grid_defs)
        p_pd%v_grid_count = 0
        IF (ALLOCATED(p_pd%opt_pvct))          DEALLOCATE(p_pd%opt_pvct)
        IF (ALLOCATED(p_pd%opt_lcall_phy))     DEALLOCATE(p_pd%opt_lcall_phy)
        IF (ALLOCATED(p_pd%opt_t_elapsed_phy)) DEALLOCATE(p_pd%opt_t_elapsed_phy)

        ! release reorder data
        CALL release_reorder_data(p_pd%cells)
        CALL release_reorder_data(p_pd%verts)
        CALL release_reorder_data(p_pd%edges)

        ! release restart file data
        CALL release_restart_file(p_pd%restart_file)

      ENDDO

      DEALLOCATE(patch_data)

    ENDIF

    ! release RMA window
    IF (mpi_win /= MPI_WIN_NULL) THEN
      CALL MPI_Win_fence(0, mpi_win, mpi_error)
      CALL check_mpi_error(routine, 'MPI_Win_fence', mpi_error, .FALSE.)
      CALL MPI_Win_free(mpi_win, mpi_error)
      CALL check_mpi_error(routine, 'MPI_Win_free', mpi_error, .FALSE.)
      mpi_win = MPI_WIN_NULL
    ENDIF

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
  SUBROUTINE print_restart_arguments()
    CHARACTER(LEN=*), PARAMETER   :: routine = modname//'print_restart_arguments'
    TYPE(t_patch_data), POINTER   :: p_pd

    WRITE (nerr,FORMAT_VALS3)routine,' is called for p_pe=',p_pe

    PRINT *,routine, ' current_caltime=', restart_args%datetime%caltime
    PRINT *,routine, ' current_calday=',  restart_args%datetime%calday
    PRINT *,routine, ' current_daysec=',  restart_args%datetime%daysec

    ! patch informations
    PRINT *,routine, ' size of patches=', SIZE(patch_data)
    p_pd => patch_data(1)
    PRINT *,routine, ' SIZE(p_pd%opt_pvct) = ', SIZE(p_pd%opt_pvct)
    PRINT *,routine, ' SIZE(p_pd%opt_lcall_phy) =', SIZE(p_pd%opt_lcall_phy)
    PRINT *,routine, ' SIZE(p_pd%opt_t_elapsed_phy) =', SIZE(p_pd%opt_t_elapsed_phy)
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
  !
  ! Transfers the restart var lists from the worker to the restart PEs.
  !
  SUBROUTINE transfer_restart_var_lists

    INTEGER :: info_size, iv, nv, nelems, n, list_info(4), ierrstat

    INTEGER, ALLOCATABLE            :: info_storage(:,:)
    TYPE(t_list_element), POINTER   :: element
    TYPE(t_var_metadata)            :: info
    TYPE(t_var_list)                :: p_var_list
    CHARACTER(LEN=MAX_NAME_LENGTH)  :: var_list_name
    CHARACTER(LEN=32)               :: model_type
    LOGICAL                         :: lrestart

    CHARACTER(LEN=*), PARAMETER :: routine = modname//'transfer_restart_var_lists'

    ! delete old var lists
    IF (my_process_is_restart()) CALL delete_var_lists

    ! get the size - in default INTEGER words - which is needed to
    ! hold the contents of TYPE(t_var_metadata)
    info_size = SIZE(TRANSFER(info, (/ 0 /)))

    ! get the number of var lists
    IF(.NOT. my_process_is_restart()) nv = nvar_lists
    CALL p_bcast(nv, bcast_root, p_comm_work_2_restart)

    ! for each var list, get its components
    DO iv = 1, nv
      ! transfer only a restart var_list
      IF (.NOT. my_process_is_restart()) lrestart = var_lists(iv)%p%lrestart
      CALL p_bcast(lrestart, bcast_root, p_comm_work_2_restart)
      IF (.NOT. lrestart) CYCLE

      ! send name of the var list
      var_list_name = ' '
      IF (.NOT. my_process_is_restart()) var_list_name = var_lists(iv)%p%name
      CALL p_bcast(var_list_name, bcast_root, p_comm_work_2_restart)

      ! send model type
      model_type = ' '
      IF (.NOT. my_process_is_restart()) model_type = var_lists(iv)%p%model_type
      CALL p_bcast(model_type, bcast_root, p_comm_work_2_restart)

      IF (.NOT. my_process_is_restart()) THEN
        ! count the number of variable restart entries
        element => var_lists(iv)%p%first_list_element
        nelems = 0
        DO
          IF(.NOT.ASSOCIATED(element)) EXIT
          IF (element%field%info%lrestart) nelems = nelems+1
          element => element%next_list_element
        ENDDO

        ! gather the components needed for name list restart and send them.
        list_info(1) = nelems
        list_info(2) = var_lists(iv)%p%patch_id
        list_info(3) = var_lists(iv)%p%restart_type
        list_info(4) = var_lists(iv)%p%vlevel_type
      ENDIF

      ! send basic info
      CALL p_bcast(list_info, bcast_root, p_comm_work_2_restart)
      ! sheck if there are valid restart fields
      IF (list_info(1) == 0) CYCLE

      IF(my_process_is_restart()) THEN
        nelems = list_info(1)
        ! create var list
        CALL new_var_list( p_var_list, var_list_name, patch_id=list_info(2), &
          &                restart_type=list_info(3), vlevel_type=list_info(4), &
          &                lrestart=.TRUE.)
        p_var_list%p%model_type = TRIM(model_type)
      ENDIF

      ALLOCATE(info_storage(info_size, nelems), STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, ALLOCATE_FAILED)

      IF(.NOT. my_process_is_restart()) THEN
        element => var_lists(iv)%p%first_list_element
        nelems = 0
        DO
          IF(.NOT. ASSOCIATED(element)) EXIT
          IF (element%field%info%lrestart) THEN
            nelems = nelems+1
            info_storage(:,nelems) = TRANSFER(element%field%info, (/ 0 /))
          ENDIF
          element => element%next_list_element
        ENDDO
      ENDIF

      ! send binary representation of all info members
      CALL p_bcast(info_storage, bcast_root, p_comm_work_2_restart)

      IF(my_process_is_restart()) THEN
        ! insert elements into var list
        p_var_list%p%first_list_element => NULL()
        element => NULL()

        DO n = 1, nelems
          IF(.NOT. ASSOCIATED(p_var_list%p%first_list_element)) THEN
            ALLOCATE(p_var_list%p%first_list_element, STAT=ierrstat)
            IF (ierrstat /= SUCCESS) CALL finish(routine, ALLOCATE_FAILED)

            element => p_var_list%p%first_list_element
          ELSE
            ALLOCATE(element%next_list_element, STAT=ierrstat)
            IF (ierrstat /= SUCCESS) CALL finish(routine, ALLOCATE_FAILED)
            element => element%next_list_element
          ENDIF

          element%next_list_element => NULL()

          ! nullify all pointers in element%field, they don't make sense on the restart PEs
          element%field%r_ptr => NULL()
          element%field%i_ptr => NULL()
          element%field%l_ptr => NULL()
          element%field%var_base_size = 0 ! Unknown here

          ! set info structure from binary representation in info_storage
          element%field%info = TRANSFER(info_storage(:, n), info)
        ENDDO
      ENDIF
      DEALLOCATE(info_storage)
    ENDDO

  END SUBROUTINE transfer_restart_var_lists

  !-------------------------------------------------------------------------------------------------
  !
  ! Create patch data and transfers this data from the worker to the restart PEs.
  !
  SUBROUTINE create_and_transfer_patch_data()
    INTEGER                        :: jg, jl, ierrstat
    TYPE(t_patch_data), POINTER    :: p_pd
    CHARACTER(LEN=*), PARAMETER :: routine = modname//'create_and_transfer_patch_data'

#ifdef DEBUG
    WRITE (nerr,FORMAT_VALS3)routine,' p_pe=',p_pe
#endif

    ! replicate domain setup
    CALL p_bcast(n_dom, bcast_root, p_comm_work_2_restart)

    ! allocate patch data structure
    ALLOCATE(patch_data(n_dom), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, ALLOCATE_FAILED)

    ! set number of global cells/edges/verts and patch ID
    DO jg = 1, n_dom

      p_pd => patch_data(jg)
      IF (my_process_is_work()) THEN
        p_pd%id              = p_patch(jg)%id
        p_pd%work_pe0_id     = p_patch(jg)%proc0
        p_pd%nlev            = p_patch(jg)%nlev
        p_pd%cell_type       = p_patch(jg)%geometry_info%cell_type
        p_pd%nblks_glb_c     = (p_patch(jg)%n_patch_cells-1)/nproma + 1
        p_pd%nblks_glb_e     = (p_patch(jg)%n_patch_edges-1)/nproma + 1
        p_pd%nblks_glb_v     = (p_patch(jg)%n_patch_verts-1)/nproma + 1
        p_pd%base_filename   = TRIM(p_patch(jg)%grid_filename)
        p_pd%n_patch_cells_g = p_patch(jg)%n_patch_cells_g
        p_pd%n_patch_verts_g = p_patch(jg)%n_patch_verts_g
        p_pd%n_patch_edges_g = p_patch(jg)%n_patch_edges_g
      END IF

      ! transfer data to restart PEs
      CALL p_bcast(p_pd%id,              bcast_root, p_comm_work_2_restart)
      CALL p_bcast(p_pd%work_pe0_id,     bcast_root, p_comm_work_2_restart)
      CALL p_bcast(p_pd%nlev,            bcast_root, p_comm_work_2_restart)
      CALL p_bcast(p_pd%cell_type,       bcast_root, p_comm_work_2_restart)
      CALL p_bcast(p_pd%nblks_glb_c,     bcast_root, p_comm_work_2_restart)
      CALL p_bcast(p_pd%nblks_glb_e,     bcast_root, p_comm_work_2_restart)
      CALL p_bcast(p_pd%nblks_glb_v,     bcast_root, p_comm_work_2_restart)
      CALL p_bcast(p_pd%base_filename,   bcast_root, p_comm_work_2_restart)
      CALL p_bcast(p_pd%n_patch_cells_g, bcast_root, p_comm_work_2_restart)
      CALL p_bcast(p_pd%n_patch_verts_g, bcast_root, p_comm_work_2_restart)
      CALL p_bcast(p_pd%n_patch_edges_g, bcast_root, p_comm_work_2_restart)

    ENDDO

    ! set reorder data
    DO jg = 1, n_dom
      p_pd => patch_data(jg)
      jl = p_pd%id

      IF(my_process_is_work()) THEN

        ! set reorder data on work PE
        CALL set_reorder_data(jg, p_patch(jl)%n_patch_cells_g, p_patch(jl)%n_patch_cells, &
                              p_patch(jl)%cells%decomp_info%owner_mask,                   &
                              p_patch(jl)%cells%decomp_info%glb_index, p_pd%cells)

        CALL set_reorder_data(jg, p_patch(jl)%n_patch_edges_g, p_patch(jl)%n_patch_edges, &
                              p_patch(jl)%edges%decomp_info%owner_mask,                   &
                              p_patch(jl)%edges%decomp_info%glb_index, p_pd%edges)

        CALL set_reorder_data(jg, p_patch(jl)%n_patch_verts_g, p_patch(jl)%n_patch_verts, &
                              p_patch(jl)%verts%decomp_info%owner_mask,                   &
                              p_patch(jl)%verts%decomp_info%glb_index, p_pd%verts)
      ENDIF

      ! transfer reorder data to restart PEs
      CALL transfer_reorder_data(p_pd%cells)
      CALL transfer_reorder_data(p_pd%edges)
      CALL transfer_reorder_data(p_pd%verts)

      ! set restart process ident
      p_pd%restart_proc_id = MOD(jg-1, process_mpi_restart_size) + p_restart_pe0

      ! reset pointer
      p_pd%v_grid_defs => NULL()
      p_pd%v_grid_count = 0

      ! set restart file data
      CALL set_restart_file_data(p_pd%restart_file, p_pd%id)

    ENDDO

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
    rf%my_mem_win_off             = 0_i8
    rf%var_data                   => NULL()
    rf%filename                   = ''
    rf%linkname                   = ''
    rf%linkprefix                 = ''

    rf%cdiIds%file                  = CDI_UNDEFID
    rf%cdiIds%vlist                 = CDI_UNDEFID
    rf%cdiIds%taxis                 = CDI_UNDEFID

    rf%cdiIds%hgrids(:)             = CDI_UNDEFID
    rf%cdiIds%vgrids(:)             = CDI_UNDEFID
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
  SUBROUTINE set_reorder_data(patch_id, n_points_g, n_points, owner_mask, glb_index, reo)

    INTEGER, INTENT(IN) :: patch_id        ! Logical patch ID
    INTEGER, INTENT(IN) :: n_points_g      ! Global number of cells/edges/verts in logical patch
    INTEGER, INTENT(IN) :: n_points        ! Local number of cells/edges/verts in logical patch
    LOGICAL, INTENT(IN) :: owner_mask(:,:) ! owner_mask for logical patch
    INTEGER, INTENT(IN) :: glb_index(:)    ! glb_index for logical patch

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
      owner_mask_1d(i) = owner_mask(il,ib)
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
        glbidx_own(n)  = glb_index(i)
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

    CHARACTER(LEN=*), PARAMETER :: routine = modname//'init_remote_memory_access'

#ifdef DEBUG
    WRITE (nerr,FORMAT_VALS3)routine,' is called for p_pe=',p_pe
#endif

    mpi_win = MPI_WIN_NULL

    ! get size and offset of the data for every restart file
    mem_size = 0_i8

    ! go over all patches
    DO i = 1, SIZE(patch_data)

      patch_data(i)%restart_file%my_mem_win_off = mem_size

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
            mem_size = mem_size + INT(nlevs*patch_data(i)%cells%n_own,i8)
          CASE (GRID_UNSTRUCTURED_EDGE)
            mem_size = mem_size + INT(nlevs*patch_data(i)%edges%n_own,i8)
          CASE (GRID_UNSTRUCTURED_VERT)
            mem_size = mem_size + INT(nlevs*patch_data(i)%verts%n_own,i8)
          CASE DEFAULT
            CALL finish(routine,UNKNOWN_GRID_TYPE)
        END SELECT

      ENDDO

      ! get the offset on all PEs
      ALLOCATE(patch_data(i)%restart_file%mem_win_off(0:num_work_procs-1), STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, ALLOCATE_FAILED)
      IF(.NOT.my_process_is_restart()) THEN
        CALL MPI_Allgather(patch_data(i)%restart_file%my_mem_win_off, 1, p_int_i8, &
                           patch_data(i)%restart_file%mem_win_off, 1, p_int_i8,    &
                           p_comm_work, mpi_error)
        CALL check_mpi_error(routine, 'MPI_Allgather', mpi_error, .TRUE.)
      ENDIF

      CALL p_bcast(patch_data(i)%restart_file%mem_win_off, bcast_root, p_comm_work_2_restart)

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
  FUNCTION find_patch(id, routine)

    INTEGER, INTENT(IN)             :: id
    CHARACTER(LEN=*), INTENT(IN)    :: routine

    INTEGER                         :: i
    TYPE(t_patch_data), POINTER     :: find_patch
    CHARACTER(LEN=MAX_ERROR_LENGTH) :: err_message

    ! try to find the patch in the modul array
    find_patch => NULL()
    DO i = 1, SIZE(patch_data)
      IF (patch_data(i)%id == id) THEN
        find_patch => patch_data(i)
        EXIT
      ENDIF
    ENDDO
    IF (.NOT. ASSOCIATED(find_patch)) THEN
      WRITE (err_message, '(a,i5)') ' patch data not found for id=',id
      CALL finish(routine, err_message)
    ENDIF
  END FUNCTION find_patch

  !------------------------------------------------------------------------------------------------
  !
  !  Set vertical grid definition.
  !
  SUBROUTINE defineVerticalGrids(patchData)
    TYPE(t_patch_data), TARGET, INTENT(INOUT) :: patchData

    INTEGER :: nlev_soil, nlev_snow, nlev_ocean, nice_class, ierrstat
    CHARACTER(LEN = *), PARAMETER :: routine = modname//":defineVerticalGrids"

    ! DEFAULT values for the level counts
    nlev_soil = 0
    nlev_snow = 0
    nlev_ocean = 0
    nice_class = 1

    ! replace DEFAULT values by the overrides provided IN the patchData
    IF (patchData%l_opt_depth_lnd) nlev_soil = patchData%opt_depth_lnd
    IF (patchData%l_opt_nlev_snow) nlev_snow = patchData%opt_nlev_snow
    IF (patchData%l_opt_depth) nlev_ocean = patchData%opt_depth
    IF (patchData%l_opt_nice_class) nice_class = patchData%opt_nice_class

    ! set vertical grid definitions
    ALLOCATE(patchData%v_grid_defs(MAX_VERTICAL_AXES), STAT=ierrstat)
    patchData%v_grid_count = 0
    IF (ierrstat /= SUCCESS) CALL finish(routine, ALLOCATE_FAILED)

    CALL set_vertical_grid(patchData%v_grid_defs, patchData%v_grid_count, ZA_SURFACE, 0._wp)
    CALL set_vertical_grid(patchData%v_grid_defs, patchData%v_grid_count, ZA_HYBRID, patchData%nlev)
    CALL set_vertical_grid(patchData%v_grid_defs, patchData%v_grid_count, ZA_HYBRID_HALF, patchData%nlev+1)
    CALL set_vertical_grid(patchData%v_grid_defs, patchData%v_grid_count, ZA_HEIGHT_2M, 2._wp)
    CALL set_vertical_grid(patchData%v_grid_defs, patchData%v_grid_count, ZA_HEIGHT_10M, 10._wp)
    CALL set_vertical_grid(patchData%v_grid_defs, patchData%v_grid_count, ZA_TOA, 1._wp)
    CALL set_vertical_grid(patchData%v_grid_defs, patchData%v_grid_count, ZA_LAKE_BOTTOM, 1._wp)
    CALL set_vertical_grid(patchData%v_grid_defs, patchData%v_grid_count, ZA_MIX_LAYER, 1._wp)
    CALL set_vertical_grid(patchData%v_grid_defs, patchData%v_grid_count, ZA_LAKE_BOTTOM_HALF, 1._wp)
    CALL set_vertical_grid(patchData%v_grid_defs, patchData%v_grid_count, ZA_SEDIMENT_BOTTOM_TW_HALF, 0._wp)
    CALL set_vertical_grid(patchData%v_grid_defs, patchData%v_grid_count, ZA_GENERIC_ICE, 1._wp)
    CALL set_vertical_grid(patchData%v_grid_defs, patchData%v_grid_count, ZA_DEPTH_RUNOFF_S, 1)
    CALL set_vertical_grid(patchData%v_grid_defs, patchData%v_grid_count, ZA_DEPTH_RUNOFF_G, 1)
    IF(patchData%l_opt_depth_lnd) CALL set_vertical_grid(patchData%v_grid_defs, patchData%v_grid_count, ZA_DEPTH_BELOW_LAND, &
                                                        &nlev_soil)
    IF(patchData%l_opt_depth_lnd) CALL set_vertical_grid(patchData%v_grid_defs, patchData%v_grid_count, ZA_DEPTH_BELOW_LAND_P1, &
                                                        &nlev_soil+1)
    IF(patchData%l_opt_nlev_snow) CALL set_vertical_grid(patchData%v_grid_defs, patchData%v_grid_count, ZA_SNOW, nlev_snow)
    IF(patchData%l_opt_nlev_snow) CALL set_vertical_grid(patchData%v_grid_defs, patchData%v_grid_count, ZA_SNOW_HALF, nlev_snow+1)
    IF(patchData%l_opt_depth) CALL set_vertical_grid(patchData%v_grid_defs, patchData%v_grid_count, ZA_DEPTH_BELOW_SEA, &
                                                        &nlev_ocean)
    IF(patchData%l_opt_depth) CALL set_vertical_grid(patchData%v_grid_defs, patchData%v_grid_count, ZA_DEPTH_BELOW_SEA_HALF, &
                                                        &nlev_ocean+1)
  END SUBROUTINE defineVerticalGrids

  !------------------------------------------------------------------------------------------------
  !
  !  Set global restart attributes.
  !
  SUBROUTINE set_restart_attributes (restartAttributes)
    TYPE(t_RestartAttributeList), POINTER, INTENT(INOUT) :: restartAttributes

    TYPE(t_patch_data), POINTER :: p_pd
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
    IF(patch_data(1)%l_opt_ndom) effectiveDomainCount = patch_data(1)%opt_ndom
    IF(ALLOCATED(restart_args%opt_output_jfile)) THEN
        CALL setGeneralRestartAttributes(restartAttributes, restart_args%datetime, effectiveDomainCount, restart_args%jstep, &
                                        &restart_args%opt_output_jfile)
    ELSE
        CALL setGeneralRestartAttributes(restartAttributes, restart_args%datetime, effectiveDomainCount, restart_args%jstep)
    END IF

    ! set the domain dependent attributes
    DO i = 1, SIZE(patch_data)
        p_pd => patch_data(i)

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
  FUNCTION has_valid_time_level (p_info,id)

    TYPE(t_var_metadata), POINTER, INTENT(IN) :: p_info
    INTEGER, INTENT(IN)           :: id

    LOGICAL                       :: has_valid_time_level

    INTEGER                       :: idx, time_level
    LOGICAL                       :: lskip_timelev, lskip_extra_timelevs

#ifdef DEBUG
    CHARACTER(LEN=*), PARAMETER   :: routine = modname//'has_valid_time_level'
#endif

    has_valid_time_level = .TRUE.
    IF (.NOT. p_info%lrestart) THEN
      has_valid_time_level = .FALSE.
      RETURN
    ENDIF

    lskip_timelev = .FALSE.
    IF (iequations == INH_ATMOSPHERE .AND. .NOT. (l_limited_area .AND. id == 1)) THEN
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
        IF (time_level == nnew(id))                    lskip_timelev = .TRUE.
        ! this is needed to skip the extra time levels allocated for nesting
        IF (lskip_extra_timelevs .AND. time_level > 2) lskip_timelev = .TRUE.
      ELSE IF (p_info%tlev_source == TLEV_NNOW_RCF) THEN
        IF (time_level == nnew_rcf(id)) lskip_timelev = .TRUE.
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
  FUNCTION get_reorder_ptr (p_pd, p_info,routine)

    TYPE(t_patch_data), POINTER, INTENT(IN)   :: p_pd
    TYPE(t_var_metadata), POINTER, INTENT(IN) :: p_info
    CHARACTER(LEN=*), INTENT(IN)              :: routine

    TYPE(t_reorder_data), POINTER             :: get_reorder_ptr

    get_reorder_ptr => NULL()

    SELECT CASE (p_info%hgrid)
      CASE (GRID_UNSTRUCTURED_CELL)
        get_reorder_ptr => p_pd%cells
      CASE (GRID_UNSTRUCTURED_EDGE)
        get_reorder_ptr => p_pd%edges
      CASE (GRID_UNSTRUCTURED_VERT)
        get_reorder_ptr => p_pd%verts
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

    TYPE(t_patch_data), POINTER, INTENT(IN) :: p_pd

    TYPE(t_restart_file), POINTER   :: p_rf
    TYPE(t_var_metadata), POINTER   :: p_info
    TYPE(t_reorder_data), POINTER   :: p_ri
    TYPE(t_datetime), POINTER       :: dt
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
    dt => restart_args%datetime
    idate = cdiEncodeDate(dt%year, dt%month, dt%day)
    itime = cdiEncodeTime(dt%hour, dt%minute, NINT(dt%second))

    p_rf => p_pd%restart_file
    CALL taxisDefVdate(p_rf%cdiIds%taxis, idate)
    CALL taxisDefVtime(p_rf%cdiIds%taxis, itime)
    status = streamDefTimestep(p_rf%cdiIds%file, p_rf%cdiTimeIndex)

    p_rf%cdiTimeIndex = p_rf%cdiTimeIndex + 1

    ! check the contained array of restart variables
    p_vars => p_pd%restart_file%var_data
    IF (.NOT. ASSOCIATED(p_vars)) RETURN

    ! get maximum number of data points in a slice and allocate tmp. variables
    nval = MAX(p_pd%cells%n_glb, p_pd%edges%n_glb, p_pd%verts%n_glb)

    ! allocate RMA memory
    ALLOCATE(var1_dp(nval*restart_chunk_size), var2_dp(nval,restart_chunk_size), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, ALLOCATE_FAILED)

    ioff(:) = p_rf%mem_win_off(:)

    ! go over the all restart variables in the associated array
    VAR_LOOP : DO iv = 1, SIZE(p_vars)

      ! get pointer to metadata
      p_info => p_vars(iv)%info

#ifdef DEBUG
      WRITE (nerr,FORMAT_VALS5I)routine,' p_pe=',p_pe,' restart pe processes field=',TRIM(p_info%name)
#endif

      ! check time level of the field
      IF (.NOT. has_valid_time_level(p_info,p_pd%id)) CYCLE

      ! get current level
      IF(p_info%ndims == 2) THEN
        nlevs = 1
      ELSE
        nlevs = p_info%used_dimensions(2)
      ENDIF

      ! get pointer to reorder data
      p_ri => get_reorder_ptr (p_pd, p_info, routine)
      
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

    TYPE(t_patch_data), POINTER, INTENT(IN) :: p_pd

    TYPE(t_restart_file), POINTER   :: p_rf
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
    p_rf => p_pd%restart_file
    ioff = p_rf%my_mem_win_off

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
      IF (.NOT. has_valid_time_level(p_info,p_pd%id)) CYCLE

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
      p_ri => get_reorder_ptr(p_pd, p_info, routine)

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
    TYPE(t_restart_file), POINTER, INTENT(IN) :: p_rf
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
      IF (.NOT. has_valid_time_level(p_info,patch_id)) CYCLE

      ! define the CDI variable
      !XXX: The code I found here simply assumed that all variables are of TYPE REAL. I have NOT changed this behavior, however, the .FALSE. constants should be replaced by something more sensible IN the future.
      CALL defineVariable(p_rf%cdiIds, p_info, .FALSE., .FALSE.)
    ENDDO
  END SUBROUTINE init_restart_variables

  !------------------------------------------------------------------------------------------------
  !
  ! Opens the restart file from the given parameters.
  !
  SUBROUTINE open_restart_file(p_pd, restartAttributes)
    TYPE(t_patch_data), POINTER, INTENT(IN) :: p_pd
    TYPE(t_RestartAttributeList), POINTER, INTENT(INOUT) :: restartAttributes

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
    CALL associate_keyword("<gridfile>",   TRIM(get_filename_noext(p_pd%base_filename)),   keywords)
    CALL associate_keyword("<idom>",       TRIM(int2string(p_pd%id, "(i2.2)")),            keywords)
    CALL associate_keyword("<rsttime>",    TRIM(datetime),                                 keywords)
    CALL associate_keyword("<mtype>",      TRIM(p_rf%model_type),                          keywords)
    ! replace keywords in file name
    p_rf%filename = TRIM(with_keywords(keywords, TRIM(restart_filename)))

    IF(ALLOCATED(p_pd%opt_pvct)) THEN
        CALL openRestartAndCreateIds(p_rf%cdiIds, TRIM(p_rf%filename), restart_type, restartAttributes, p_pd%cells%n_glb, &
                                    &p_pd%verts%n_glb, p_pd%edges%n_glb, p_pd%cell_type, p_pd%v_grid_defs(1:p_pd%v_grid_count), &
                                    &p_pd%opt_pvct)
    ELSE
        CALL openRestartAndCreateIds(p_rf%cdiIds, TRIM(p_rf%filename), restart_type, restartAttributes, p_pd%cells%n_glb, &
                                    &p_pd%verts%n_glb, p_pd%edges%n_glb, p_pd%cell_type, p_pd%v_grid_defs(1:p_pd%v_grid_count))
    END IF

#ifdef DEBUG
    WRITE (nerr, FORMAT_VALS5)routine,' p_pe=',p_pe,' open netCDF file with ID=',p_rf%cdiIds%file
#endif

    ! set cdi internal time index to 0 for writing time slices in netCDF
    p_rf%cdiTimeIndex = 0

    ! init list of restart variables
    CALL init_restart_variables (p_rf,p_pd%id)

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

    CALL closeAndDestroyIds(rf%cdiIds)
    rf%filename = ''
    rf%linkname = ''
    rf%linkprefix = ''

  END SUBROUTINE close_restart_file

  !------------------------------------------------------------------------------------------------
  !
  ! Creates a symbolic link from the given restart file.
  !
  SUBROUTINE create_restart_file_link (rf, proc_id, jg, l_opt_ndom, opt_ndom)

    TYPE (t_restart_file), INTENT(INOUT)  :: rf
    INTEGER,               INTENT(IN)     :: proc_id
    INTEGER,               INTENT(IN)     :: jg                   !< patch ID
    LOGICAL                               :: l_opt_ndom
    INTEGER                               :: opt_ndom

    INTEGER                               :: iret, id
    CHARACTER(LEN=5)                      :: str_id
    CHARACTER(LEN=*), PARAMETER           :: routine = modname//'create_restart_file_link'

    ! build link name
    id = proc_id - p_restart_pe0
    IF (id == 0) THEN
      str_id = ' '
    ELSE
      WRITE(str_id, '(I5)')id
    ENDIF
    rf%linkprefix = 'restart'//TRIM(str_id)
    IF (l_opt_ndom .AND. (opt_ndom > 1)) THEN
      rf%linkname = TRIM(rf%linkprefix)//'_'//TRIM(rf%model_type)//"_DOM"//TRIM(int2string(jg, "(i2.2)"))//'.nc'
    ELSE
      rf%linkname = TRIM(rf%linkprefix)//'_'//TRIM(rf%model_type)//'_DOM01.nc'
    END IF

    ! delete old symbolic link, if exists
    IF (util_islink(TRIM(rf%linkname))) THEN
      iret = util_unlink(TRIM(rf%linkname))
      IF (iret /= SUCCESS) THEN
          WRITE (nerr,'(3a)')routine,' cannot unlink ',TRIM(rf%linkname)
      ENDIF
    ENDIF

    ! create a new symbolic link
    iret = util_symlink(TRIM(rf%filename),TRIM(rf%linkname))
    IF (iret /= SUCCESS) THEN
      WRITE (nerr,'(5a)')routine,' cannot create symbolic link ', &
        & TRIM(rf%linkname),' for ', TRIM(rf%filename)
    ENDIF

  END SUBROUTINE create_restart_file_link

#endif

END MODULE mo_io_restart_async


