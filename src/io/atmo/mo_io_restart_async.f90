!>
!! Contains routines for asynchronous restart Output
!! --------------------------------------------------------
!!
!!
!! @par Revision History
!! Initial implementation by Joerg Benkenstein (2013-01-15)
!!
!! @par Copyright
!! 2002-2013 by DWD and MPI-M
!! This software is provided for non-commercial use only.
!! See the LICENSE and the WARRANTY conditions.
!!
!! @par License
!! The use of ICON is hereby granted free of charge for an unlimited time,
!! provided the following rules are accepted and applied:
!! <ol>
!! <li> You may use or modify this code for your own non commercial and non
!!    violent purposes.
!! <li> The code may not be re-distributed without the consent of the authors.
!! <li> The copyright notice and statement of authorship must appear in all
!!    copies.
!! <li> You accept the warranty conditions (see WARRANTY).
!! <li> In case you intend to use the code commercially, we oblige you to sign
!!    an according license agreement with DWD and MPI-M.
!! </ol>
!!
!! @par Warranty
!! This code has been tested up to a certain level. Defects and weaknesses,
!! which may be included in the code, do not establish any warranties by the
!! authors.
!! The authors do not make any warranty, express or implied, or assume any
!! liability or responsibility for the use, acquisition or application of this
!! software.
!!
!!

#if ! (defined (__GNUC__) || defined(__SX__) || defined(__SUNPRO_F95) || defined(__INTEL_COMPILER) || defined (__PGI))
#define HAVE_F2003
#endif
MODULE mo_io_restart_async

  USE mo_util_file,               ONLY: util_symlink, util_unlink, util_islink
  USE mo_exception,               ONLY: finish, message, message_text, get_filename_noext
  USE mo_kind,                    ONLY: wp, i8, dp
  USE mo_datetime,                ONLY: t_datetime, iso8601
  USE mo_io_config,               ONLY: out_expname
  USE mo_io_units,                ONLY: nerr, filename_max, find_next_free_unit
  USE mo_var_list,                ONLY: nvar_lists, var_lists, new_var_list, delete_var_lists
  USE mo_linked_list,             ONLY: t_list_element, t_var_list
  USE mo_io_restart_attributes,   ONLY: set_restart_attribute, delete_attributes, get_restart_attribute, &
    &                                   restart_attributes_count_int, restart_attributes_count_text, &
    &                                   restart_attributes_count_real, restart_attributes_count_bool
  USE mo_dynamics_config,         ONLY: nold, nnow, nnew, nnew_rcf, nnow_rcf, iequations
  USE mo_impl_constants,          ONLY: IHS_ATM_TEMP, IHS_ATM_THETA, ISHALLOW_WATER, &
    &                                   LEAPFROG_EXPL, LEAPFROG_SI, SUCCESS, MAX_CHAR_LENGTH
  ! USE mo_oce_state,               ONLY: set_zlev
  USE mo_var_metadata,            ONLY: t_var_metadata
  USE mo_io_restart_namelist,     ONLY: nmls, restart_namelist, delete_restart_namelists, &
    &                                   set_restart_namelist, get_restart_namelist
  USE mo_name_list_output_init,   ONLY: output_file
  USE mo_name_list_output_types,  ONLY: max_z_axes
#ifdef USE_CRAY_POINTER
  USE mo_name_list_output_init,   ONLY: set_mem_ptr_dp
#endif
#ifndef HAVE_F2003
  USE mo_io_restart_namelist,     ONLY: nmllen_max
#endif
  USE mo_communication,           ONLY: idx_no, blk_no
  USE mo_parallel_config,         ONLY: nproma
  USE mo_grid_config,             ONLY: n_dom
  USE mo_ha_dyn_config,           ONLY: ha_dyn_config
  USE mo_model_domain,            ONLY: p_patch
  USE mo_util_sysinfo,            ONLY: util_user_name, util_os_system, util_node_name
  USE mo_cdi_constants

#ifndef NOMPI
  USE mo_mpi,                     ONLY: p_pe, p_pe_work, p_restart_pe0, p_comm_work, &
    &                                   p_work_pe0, num_work_procs, MPI_SUCCESS, &
    &                                   p_stop, p_send, p_recv, p_barrier, p_bcast, &
    &                                   my_process_is_restart, my_process_is_work, &
    &                                   p_comm_work_2_restart, p_n_work, p_int, &
    &                                   process_mpi_restart_size, p_int_i8, p_real_dp, &
    &                                   p_comm_work_restart

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
#ifndef HAVE_F2003
  INTEGER, PARAMETER :: MAX_ATTRIB_TLENGTH        = nmllen_max
#endif

  ! maximum no. of output_nml output files
  INTEGER, parameter :: MAX_NML_OUTPUT_FILES      = 100

  ! minimum number of dynamic restart arguments
  INTEGER, PARAMETER :: MIN_DYN_RESTART_ARGS      = 16 + MAX_NML_OUTPUT_FILES

  ! minimum number of dynamic restart patch data
  ! id, l_dom_active, time levels and optional attributes
  INTEGER, PARAMETER :: MIN_DYN_RESTART_PDATA     = 21

  ! maximumm number of verticale axes
  INTEGER, PARAMETER :: MAX_VERTICAL_AXES         = 19

  ! common constant strings
  CHARACTER(LEN=*), PARAMETER :: MODUL_NAME               = 'shared/mo_io_restart_async/'
  CHARACTER(LEN=*), PARAMETER :: MODEL_TITLE              = 'ICON simulation'
  CHARACTER(LEN=*), PARAMETER :: MODEL_INSTITUTION        = &
    &                            'Max Planck Institute for Meteorology/Deutscher Wetterdienst'
  CHARACTER(LEN=*), PARAMETER :: MODEL_VERSION            = '1.2.2'
  CHARACTER(LEN=*), PARAMETER :: MODEL_REFERENCES         = 'see MPIM/DWD publications'
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

#ifdef DEBUG
  CHARACTER(LEN=*), PARAMETER :: FORMAT_VALS3             = '(a,a,i3)'
  CHARACTER(LEN=*), PARAMETER :: FORMAT_VALS5             = '(a,a,i3,a,i3)'
  CHARACTER(LEN=*), PARAMETER :: FORMAT_VALS5I            = '(a,a,i3,a,a)'
  CHARACTER(LEN=*), PARAMETER :: FORMAT_VALS7             = '(a,a,i3,a,i6,a,i3)'
  CHARACTER(LEN=*), PARAMETER :: FORMAT_VALS7I            = '(a,a,i3,a,a,a,i8)'
  CHARACTER(LEN=*), PARAMETER :: FORMAT_VALS9             = '(a,a,i3,a,i6,a,i3,a,i3)'
#endif

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
    INTEGER                     :: cdiFileID
    INTEGER                     :: cdiVlistID
    INTEGER                     :: cdiCellGridID
    INTEGER                     :: cdiVertGridID
    INTEGER                     :: cdiEdgeGridID
    INTEGER                     :: cdiTaxisID
    INTEGER                     :: cdiZaxisIDs(max_z_axes)
    INTEGER                     :: cdiTimeIndex

  END TYPE t_restart_file

  !------------------------------------------------------------------------------------------------
  ! TYPE t_reorder_data describes how local cells/edges/verts
  ! have to be reordered to get the global array.
  ! Below, "points" refers to either cells, edges or verts.
  !
  TYPE t_reorder_data
    INTEGER :: n_glb  ! Global number of points per physical patch
    INTEGER :: n_own  ! Number of own points (without halo, only belonging to phyiscal patch)
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
  ! TYPE t_v_grid contains the data of a vertical grid definition.
  !
  TYPE t_v_grid
    INTEGER :: type
    INTEGER :: nlevels
  END type t_v_grid

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

    ! base file name contains already logical patch ident
    CHARACTER(LEN=filename_max) :: base_filename

    ! dynamic patch arguments (mandotary)
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
    LOGICAL               :: l_opt_jstep_adv_ntstep
    INTEGER               :: opt_jstep_adv_ntstep
    LOGICAL               :: l_opt_jstep_adv_marchuk_order
    INTEGER               :: opt_jstep_adv_marchuk_order
    LOGICAL               :: l_opt_sim_time
    REAL(wp)              :: opt_sim_time
    INTEGER               :: n_opt_pvct
    REAL(wp), ALLOCATABLE :: opt_pvct(:)
    INTEGER               :: n_opt_lcall_phy
    LOGICAL, ALLOCATABLE  :: opt_lcall_phy(:)
    INTEGER               :: n_opt_t_elapsed_phy
    REAL(wp), ALLOCATABLE :: opt_t_elapsed_phy(:)
  END TYPE t_patch_data
  TYPE(t_patch_data), ALLOCATABLE, TARGET :: patch_data (:)

  !------------------------------------------------------------------------------------------------
  ! patch independent arguments
  !
  TYPE t_restart_args
    TYPE(t_datetime)  :: datetime
    INTEGER           :: jstep
#ifdef HAVE_F2003
    INTEGER           :: n_max_attrib_tlen
#endif
    INTEGER           :: noutput_files
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
  SUBROUTINE prepare_async_restart (opt_lcall_phy_size, opt_t_elapsed_phy_size, &
     &                              opt_pvct_size)

    INTEGER,  INTENT(IN), OPTIONAL :: opt_lcall_phy_size, opt_t_elapsed_phy_size, &
      &                               opt_pvct_size

    CHARACTER(LEN=*), PARAMETER :: subname = MODUL_NAME//'prepare_async_restart'

#ifdef NOMPI
    CALL finish(subname, ASYNC_RESTART_REQ_MPI)
#else

#ifdef DEBUG
    WRITE (nerr,FORMAT_VALS3)subname,' p_pe=',p_pe
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
    CALL transfer_restart_name_lists

    ! create and transfer patch data
    CAll create_and_transfer_patch_data(opt_lcall_phy_size, &
      &                                 opt_t_elapsed_phy_size, &
      &                                 opt_pvct_size)

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
                                &   opt_jstep_adv_ntstep,         &
                                &   opt_jstep_adv_marchuk_order,  &
                                &   opt_depth,                    &
                                &   opt_depth_lnd,                &
                                &   opt_nlev_snow,                &
                                &   opt_nice_class)

    INTEGER,  INTENT(IN)           :: patch_id
    LOGICAL,  INTENT(IN)           :: l_dom_active

    INTEGER,  INTENT(IN), OPTIONAL :: opt_depth
    INTEGER,  INTENT(IN), OPTIONAL :: opt_depth_lnd
    INTEGER,  INTENT(IN), OPTIONAL :: opt_jstep_adv_ntstep
    INTEGER,  INTENT(IN), OPTIONAL :: opt_jstep_adv_marchuk_order
    INTEGER,  INTENT(IN), OPTIONAL :: opt_nlev_snow
    INTEGER,  INTENT(IN), OPTIONAL :: opt_nice_class
    REAL(wp), INTENT(IN), OPTIONAL :: opt_sim_time
    LOGICAL , INTENT(IN), OPTIONAL :: opt_lcall_phy(:)
    REAL(wp), INTENT(IN), OPTIONAL :: opt_pvct(:)
    REAL(wp), INTENT(IN), OPTIONAL :: opt_t_elapsed_phy(:)

    TYPE(t_patch_data), POINTER    :: p_pd
    INTEGER                        :: ierrstat
    CHARACTER(LEN=*), PARAMETER    :: subname = MODUL_NAME//'set_data_async_restart'

#ifdef NOMPI
    CALL finish(subname, ASYNC_RESTART_REQ_MPI)
#else

#ifdef DEBUG
    WRITE (nerr,FORMAT_VALS3)subname,' p_pe=',p_pe
#endif

    IF(.NOT. (my_process_is_work())) RETURN

    ! only the first compute PE needs the dynamic restart arguments
    IF(p_pe_work == 0) THEN

      ! find patch
      p_pd => find_patch(patch_id, subname)

      ! set activity flag
      p_pd%l_dom_active = l_dom_active

      ! copy optional array parameter
      IF (PRESENT(opt_pvct)) THEN
        IF (SIZE(opt_pvct) /= p_pd%n_opt_pvct) THEN
          CALL finish(subname, WRONG_ARRAY_SIZE//'opt_pvct')
        ENDIF
        IF (.NOT. ALLOCATED (p_pd%opt_pvct)) THEN
          ALLOCATE(p_pd%opt_pvct(p_pd%n_opt_pvct), STAT=ierrstat)
          IF (ierrstat /= SUCCESS) CALL finish(subname, ALLOCATE_FAILED)
        ENDIF
        p_pd%opt_pvct = opt_pvct
      ENDIF

      IF (PRESENT(opt_t_elapsed_phy)) THEN
        IF (SIZE(opt_t_elapsed_phy) /= p_pd%n_opt_t_elapsed_phy) THEN
          CALL finish(subname, WRONG_ARRAY_SIZE//'opt_t_elapsed_phy')
        ENDIF
        IF (.NOT. ALLOCATED(p_pd%opt_t_elapsed_phy)) THEN
          ALLOCATE(p_pd%opt_t_elapsed_phy(p_pd%n_opt_t_elapsed_phy), STAT=ierrstat)
          IF (ierrstat /= SUCCESS) CALL finish(subname, ALLOCATE_FAILED)
        ENDIF
        p_pd%opt_t_elapsed_phy = opt_t_elapsed_phy
      ENDIF
!
      IF (PRESENT(opt_lcall_phy)) THEN
        IF (SIZE(opt_lcall_phy) /= p_pd%n_opt_lcall_phy) THEN
          CALL finish(subname, WRONG_ARRAY_SIZE//'opt_lcall_phy')
        ENDIF
        IF (.NOT. ALLOCATED(p_pd%opt_lcall_phy)) THEN
          ALLOCATE(p_pd%opt_lcall_phy(p_pd%n_opt_lcall_phy), STAT=ierrstat)
          IF (ierrstat /= SUCCESS) CALL finish(subname, ALLOCATE_FAILED)
        ENDIF
        p_pd%opt_lcall_phy = opt_lcall_phy
      ENDIF

      ! copy optional value parameter
      IF (PRESENT(opt_jstep_adv_ntstep)) THEN
        p_pd%opt_jstep_adv_ntstep = opt_jstep_adv_ntstep
        p_pd%l_opt_jstep_adv_ntstep = .TRUE.
      ELSE
        p_pd%l_opt_jstep_adv_ntstep = .FALSE.
      ENDIF

      IF (PRESENT(opt_jstep_adv_marchuk_order)) THEN
        p_pd%opt_jstep_adv_marchuk_order = opt_jstep_adv_marchuk_order
        p_pd%l_opt_jstep_adv_marchuk_order = .TRUE.
      ELSE
        p_pd%l_opt_jstep_adv_marchuk_order = .FALSE.
      ENDIF

      IF (PRESENT(opt_depth)) THEN
        p_pd%opt_depth = opt_depth
        p_pd%l_opt_depth = .TRUE.
      ELSE
        p_pd%l_opt_depth = .FALSE.
      ENDIF

      IF (PRESENT(opt_depth_lnd)) THEN
        p_pd%opt_depth_lnd = opt_depth_lnd
        p_pd%l_opt_depth_lnd = .TRUE.
      ELSE
        p_pd%l_opt_depth_lnd = .FALSE.
      ENDIF

      IF (PRESENT(opt_nlev_snow)) THEN
        p_pd%opt_nlev_snow = opt_nlev_snow
        p_pd%l_opt_nlev_snow = .TRUE.
      ELSE
        p_pd%l_opt_nlev_snow = .FALSE.
      ENDIF

      IF (PRESENT(opt_nice_class)) THEN
        p_pd%opt_nice_class = opt_nice_class
        p_pd%l_opt_nice_class = .TRUE.
      ELSE
        p_pd%l_opt_nice_class = .FALSE.
      ENDIF

      IF (PRESENT(opt_sim_time)) THEN
        p_pd%opt_sim_time = opt_sim_time
        p_pd%l_opt_sim_time = .TRUE.
      ELSE
        p_pd%l_opt_sim_time = .FALSE.
      ENDIF
  ENDIF
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
    INTEGER                         :: idx, noutput_files, j

    CHARACTER(LEN=*), PARAMETER :: subname = MODUL_NAME//'write_async_restart'

#ifdef NOMPI
    CALL finish (subname, ASYNC_RESTART_REQ_MPI)
#else

#ifdef DEBUG
    WRITE (nerr,FORMAT_VALS3)subname,' p_pe=',p_pe
#endif

    ! check kind of process
    IF (.NOT. my_process_is_work() .AND. .NOT. my_process_is_restart()) RETURN

    IF (my_process_is_work()) THEN
      CALL compute_wait_for_restart
      noutput_files = SIZE(output_file,1)
      CALL compute_start_restart (datetime, jstep, noutput_files)
    END IF

    ! do the restart output
    DO idx = 1, SIZE(patch_data)

      p_pd => patch_data(idx)

      ! check, if the patch is actice
      IF (.NOT. p_pd%l_dom_active) CYCLE

      ! write the variable restart lists
      IF (my_process_is_restart()) THEN

        ! considerate the right restart process
        IF (p_pe == p_pd%restart_proc_id) THEN

          ! set global restart attributes/lists
          CALL set_restart_attributes(p_pd)

#ifdef DEBUG
          CALL print_restart_arguments
          CALL print_restart_attributes
          CALL print_restart_name_lists
#endif
          CALL open_restart_file(p_pd)

          ! collective call to write the restart variables
          CALL restart_write_var_list(p_pd)

          CALL create_restart_file_link(p_pd%restart_file, p_pd%restart_proc_id)
          CALL create_restart_info_file(p_pd%restart_file)
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
    CHARACTER(LEN=*), PARAMETER :: subname = MODUL_NAME//'close_async_restart'

    CALL finish(subname, ASYNC_RESTART_REQ_MPI)
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

    TYPE(t_restart_args), POINTER :: p_ra
    LOGICAL                       :: done

    CHARACTER(LEN=*), PARAMETER :: subname = MODUL_NAME//'restart_main_proc'

#ifdef NOMPI
    CALL finish(subname, ASYNC_RESTART_REQ_MPI)
#else

#ifdef DEBUG
    WRITE (nerr,FORMAT_VALS3)subname,' p_pe=',p_pe
#endif

    ! check kind of process
    IF (.NOT. my_process_is_restart()) RETURN

    ! prepare restart (collective call)
    CALL prepare_async_restart()

    ! tell the compute PEs that we are ready to work
    CALL restart_send_ready

    ! enter restart loop
    done = .FALSE.
    p_ra => restart_args
    DO
      ! wait for a message from the compute PEs to start
      CALL restart_wait_for_start(done)

      IF(done) EXIT ! leave loop, we are done

      ! read and write restart variable lists (collective call)
      CALL write_async_restart (p_ra%datetime,    &
                              & p_ra%jstep)

      ! inform compute PEs that the restart is done
      CALL restart_send_ready
    ENDDO

    ! finalization sequence (collective call)
    CALL close_async_restart

    ! shut down MPI
    CALL p_stop

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
    CHARACTER(LEN=*), PARAMETER :: subname = MODUL_NAME//'restart_send_ready'

    WRITE (nerr,FORMAT_VALS5)subname,' p_pe=',p_pe, &
      & ' call p_barrier with communicator=',p_comm_work
#endif
    ! make sure all are done
    CALL p_barrier(comm=p_comm_work)

    ! simply send a message from restart PE 0 to compute PE 0
    IF (p_pe_work == 0) THEN
      msg = REAL(MSG_RESTART_DONE, wp)
#ifdef DEBUG
      WRITE (nerr,FORMAT_VALS7)subname,' p_pe=',p_pe, &
        & ' send message=',INT(msg),' to pe=',p_work_pe0
#endif
      CALL p_send(msg, p_work_pe0, 0)
    ENDIF

  END SUBROUTINE restart_send_ready

  !-------------------------------------------------------------------------------------------------
  !>
  !! restart_wait_for_start: Wait for a message from compute PEs that we should start restart or finish.
  !! The counterpart on the compute side is compute_start_restart/compute_shutdown_restart.
  !
  SUBROUTINE restart_wait_for_start(done)

    LOGICAL, INTENT(OUT)           :: done ! flag if we should shut down

    TYPE(t_patch_data), POINTER    :: p_pd
    TYPE(t_restart_args), POINTER  :: p_ra
    INTEGER                        :: i, j, k, ierrstat
    REAL(wp), POINTER              :: p_msg(:)
    CHARACTER(LEN=*), PARAMETER    :: subname = MODUL_NAME//'restart_wait_for_start'

    ! set output parameter to default value
    done = .FALSE.

#ifdef DEBUG
    WRITE (nerr,FORMAT_VALS3)subname,' is called, p_pe=',p_pe
#endif

    ! create message array
    p_msg => get_message_array(subname)

    ! receive message that we may start restart (or should finish)
    IF(p_pe_work == 0) THEN
      CALL p_recv(p_msg, p_work_pe0, 0)
#ifdef DEBUG
      WRITE (nerr,FORMAT_VALS9)subname,' p_pe=',p_pe, &
        & ' p_recv got msg=',INT(p_msg(1)),' len=',SIZE(p_msg), &
        & ' from pe=',p_work_pe0
#endif
    ENDIF

#ifdef DEBUG
    WRITE (nerr,FORMAT_VALS5)subname,' p_pe=',p_pe, &
      & ' call p_bcast with communicator=',p_comm_work
#endif
    CALL p_bcast(p_msg, 0, comm=p_comm_work)

    SELECT CASE(INT(p_msg(1)))

      CASE(MSG_RESTART_START)

#ifdef DEBUG
        DO i = 1, SIZE(p_msg)
          PRINT *,subname, ' p_msg(i)=',i,': ',p_msg(i)
        ENDDO
#endif

        done = .FALSE.

        ! get patch independent arguments
        p_ra => restart_args
        p_ra%datetime%year    = INT (p_msg(4))
        p_ra%datetime%month   = INT (p_msg(5))
        p_ra%datetime%day     = INT (p_msg(6))
        p_ra%datetime%hour    = INT (p_msg(7))
        p_ra%datetime%minute  = INT (p_msg(8))
        p_ra%datetime%second  = REAL(p_msg(9),  wp)
        p_ra%datetime%caltime = REAL(p_msg(10), wp)
        p_ra%datetime%calday  = INT (p_msg(11), i8)
        p_ra%datetime%daysec  = REAL(p_msg(12), wp)
        nnew(1)               = INT (p_msg(13))
        nnew_rcf(1)           = INT (p_msg(14))

        p_ra%noutput_files    = INT (p_msg(15))
        p_ra%jstep            = INT (p_msg(16))

        ! get patch dependent arguments
        i = MIN_DYN_RESTART_ARGS

        DO j = 1, SIZE(patch_data)

          ! find the patch of the current patch id
          p_pd => find_patch(INT(p_msg(incr(i))), subname)

          ! activity flag
          p_pd%l_dom_active = get_flag(p_msg(incr(i)))

          ! time levels
          p_pd%nold     = INT(p_msg(incr(i)))
          p_pd%nnow     = INT(p_msg(incr(i)))
          p_pd%nnew     = INT(p_msg(incr(i)))
          p_pd%nnew_rcf = INT(p_msg(incr(i)))
          p_pd%nnow_rcf = INT(p_msg(incr(i)))

          ! optional parameter values
          p_pd%l_opt_depth                    = get_flag(p_msg(incr(i)))
          p_pd%opt_depth                      = p_msg(incr(i))
          p_pd%l_opt_depth_lnd                = get_flag(p_msg(incr(i)))
          p_pd%opt_depth_lnd                  = p_msg(incr(i))
          p_pd%l_opt_nlev_snow                = get_flag(p_msg(incr(i)))
          p_pd%opt_nlev_snow                  = p_msg(incr(i))
          p_pd%l_opt_nice_class               = get_flag(p_msg(incr(i)))
          p_pd%opt_nice_class                 = p_msg(incr(i))
          p_pd%l_opt_jstep_adv_ntstep         = get_flag(p_msg(incr(i)))
          p_pd%opt_jstep_adv_ntstep           = p_msg(incr(i))
          p_pd%l_opt_jstep_adv_marchuk_order  = get_flag(p_msg(incr(i)))
          p_pd%opt_jstep_adv_marchuk_order    = p_msg(incr(i))
          p_pd%l_opt_sim_time                 = get_flag(p_msg(incr(i)))
          p_pd%opt_sim_time                   = p_msg(incr(i))

          ! optional parameter arrays
          IF (p_pd%n_opt_pvct > 0) THEN
            IF (.NOT. ALLOCATED(p_pd%opt_pvct)) THEN
              ALLOCATE(p_pd%opt_pvct(p_pd%n_opt_pvct), STAT=ierrstat)
              IF (ierrstat /= SUCCESS) CALL finish(subname, ALLOCATE_FAILED)
            ENDIF
            DO k = 1, SIZE(p_pd%opt_pvct)
              p_pd%opt_pvct(k) = p_msg(incr(i))
            ENDDO
          ENDIF

          IF (p_pd%n_opt_lcall_phy > 0) THEN
            IF (.NOT. ALLOCATED(p_pd%opt_lcall_phy)) THEN
              ALLOCATE(p_pd%opt_lcall_phy(p_pd%n_opt_lcall_phy), STAT=ierrstat)
              IF (ierrstat /= SUCCESS) CALL finish(subname, ALLOCATE_FAILED)
            ENDIF
            DO k = 1, SIZE(p_pd%opt_lcall_phy)
               p_pd%opt_lcall_phy(k) = get_flag(p_msg(incr(i)))
            ENDDO
          ENDIF

          IF (p_pd%n_opt_t_elapsed_phy > 0) THEN
            IF (.NOT. ALLOCATED(p_pd%opt_t_elapsed_phy)) THEN
              ALLOCATE(p_pd%opt_t_elapsed_phy(p_pd%n_opt_t_elapsed_phy), STAT=ierrstat)
              IF (ierrstat /= SUCCESS) CALL finish(subname, ALLOCATE_FAILED)
            ENDIF
            DO k = 1, SIZE(p_pd%opt_t_elapsed_phy)
               p_pd%opt_t_elapsed_phy(k) = p_msg(incr(i))
            ENDDO
          ENDIF
        ENDDO

      CASE(MSG_RESTART_SHUTDOWN)
        done = .TRUE.

      CASE DEFAULT
        ! anything else is an error
        CALL finish(subname,'restart PE: Got illegal restart tag')

    END SELECT

    IF (ASSOCIATED(p_msg)) DEALLOCATE(p_msg)

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

    CHARACTER(LEN=*), PARAMETER :: subname = MODUL_NAME//'compute_wait_for_restart'

#ifdef DEBUG
    WRITE (nerr,FORMAT_VALS3)subname,' is called, p_pe=',p_pe
#endif

    ! first compute PE receives message from restart leader
    IF(p_pe_work == 0) THEN
      CALL p_recv(msg, p_restart_pe0, 0)
#ifdef DEBUG
      WRITE (nerr,FORMAT_VALS7)subname,' p_pe=',p_pe, &
        & ' p_recv got msg=',INT(msg),' from pe=',p_restart_pe0
#endif
      ! just for safety: Check if we got the correct tag
      IF(INT(msg) /= MSG_RESTART_DONE) THEN
        CALL finish(subname,'Compute PE: Got illegal restart tag')
      ENDIF
    ENDIF

#ifdef DEBUG
    WRITE (nerr,FORMAT_VALS5)subname,' p_pe=',p_pe, &
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
  SUBROUTINE compute_start_restart(datetime, jstep, noutput_files)

    TYPE(t_datetime), INTENT(IN)  :: datetime
    INTEGER,          INTENT(IN)  :: jstep
    INTEGER,          INTENT(IN)  :: noutput_files

    TYPE(t_patch_data), POINTER   :: p_pd
    REAL(wp), POINTER             :: p_msg(:)
    INTEGER                       :: i, j, k
    CHARACTER(LEN=*), PARAMETER :: subname = MODUL_NAME//'compute_start_restart'

#ifdef DEBUG
    WRITE (nerr,FORMAT_VALS5)subname,' p_pe=',p_pe, &
      & ' call p_barrier with communicator=',p_comm_work
#endif

    ! make sure all are here
    CALL p_barrier(comm=p_comm_work)

    IF(p_pe_work == 0) THEN

      ! create message array
      p_msg => get_message_array(subname)

      ! set command id
      p_msg(1) = REAL(MSG_RESTART_START,  wp)

      ! set patch independent arguments
      p_msg(4)  = REAL(datetime%year,     wp)
      p_msg(5)  = REAL(datetime%month,    wp)
      p_msg(6)  = REAL(datetime%day,      wp)
      p_msg(7)  = REAL(datetime%hour,     wp)
      p_msg(8)  = REAL(datetime%minute,   wp)
      p_msg(9)  = REAL(datetime%second,   wp)
      p_msg(10) = REAL(datetime%caltime,  wp)
      p_msg(11) = REAL(datetime%calday,   wp)
      p_msg(12) = REAL(datetime%daysec,   wp)
      p_msg(13) = REAL(nnew(1),           wp)
      p_msg(14) = REAL(nnew_rcf(1),       wp)

      p_msg(15) = REAL(noutput_files,     wp)
      p_msg(16) = REAL(jstep,             wp)

      ! set data of all patches
      i = MIN_DYN_RESTART_ARGS

      DO j = 1, SIZE(patch_data)

        p_pd => patch_data(j)

        ! patch id
        p_msg(incr(i)) = p_pd%id

        ! activity flag
        p_msg(incr(i)) = MERGE(1._wp, 0._wp, p_pd%l_dom_active)

        ! time levels
        p_msg(incr(i)) = REAL(nold(p_pd%id),    wp)
        p_msg(incr(i)) = REAL(nnow(p_pd%id),    wp)
        p_msg(incr(i)) = REAL(nnew(p_pd%id),    wp)
        p_msg(incr(i)) = REAL(nnew_rcf(p_pd%id),wp)
        p_msg(incr(i)) = REAL(nnow_rcf(p_pd%id),wp)

        ! optional parameter values
        p_msg(incr(i)) = MERGE(1._wp, 0._wp, p_pd%l_opt_depth)
        p_msg(incr(i)) = REAL(p_pd%opt_depth, wp)
        p_msg(incr(i)) = MERGE(1._wp, 0._wp, p_pd%l_opt_depth_lnd)
        p_msg(incr(i)) = REAL(p_pd%opt_depth_lnd, wp)
        p_msg(incr(i)) = MERGE(1._wp, 0._wp, p_pd%l_opt_nlev_snow)
        p_msg(incr(i)) = REAL(p_pd%opt_nlev_snow, wp)
        p_msg(incr(i)) = MERGE(1._wp, 0._wp, p_pd%l_opt_nice_class)
        p_msg(incr(i)) = REAL(p_pd%opt_nice_class, wp)
        p_msg(incr(i)) = MERGE(1._wp, 0._wp, p_pd%l_opt_jstep_adv_ntstep)
        p_msg(incr(i)) = REAL(p_pd%opt_jstep_adv_ntstep, wp)
        p_msg(incr(i)) = MERGE(1._wp, 0._wp, p_pd%l_opt_jstep_adv_marchuk_order)
        p_msg(incr(i)) = REAL(p_pd%opt_jstep_adv_marchuk_order, wp)
        p_msg(incr(i)) = MERGE(1._wp, 0._wp, p_pd%l_opt_sim_time)
        p_msg(incr(i)) = p_pd%opt_sim_time

        ! optional parameter arrays
        IF (ALLOCATED(p_pd%opt_pvct)) THEN
          DO k = 1, SIZE(p_pd%opt_pvct)
             p_msg(incr(i)) = REAL(p_pd%opt_pvct(k), wp)
          ENDDO
        ENDIF
        IF (ALLOCATED(p_pd%opt_lcall_phy)) THEN
          DO k = 1, SIZE(p_pd%opt_lcall_phy)
             p_msg(incr(i)) = MERGE(1._wp, 0._wp, p_pd%opt_lcall_phy(k))
          ENDDO
        ENDIF
        IF (ALLOCATED(p_pd%opt_t_elapsed_phy)) THEN
          DO k = 1, SIZE(p_pd%opt_t_elapsed_phy)
             p_msg(incr(i)) = REAL(p_pd%opt_t_elapsed_phy(k), wp)
          ENDDO
        ENDIF
      ENDDO

#ifdef DEBUG
      DO i = 1, SIZE(p_msg)
        PRINT *,subname, ' p_msg(i)=',i,': ',p_msg(i)
      ENDDO

      WRITE (nerr,FORMAT_VALS9)subname,' p_pe=',p_pe, &
        & ' send message=',INT(p_msg(1)),' len=',SIZE(p_msg), &
        & ' to pe=',p_restart_pe0
#endif
      CALL p_send(p_msg, p_restart_pe0, 0)

      IF (ASSOCIATED(p_msg)) DEALLOCATE(p_msg)
    ENDIF

  END SUBROUTINE compute_start_restart

  !-------------------------------------------------------------------------------------------------
  !>
  !! compute_shutdown_restart: Send a message to restart PEs that they should shut down.
  !! The counterpart on the restart side is restart_wait_for_start.
  !
  SUBROUTINE compute_shutdown_restart

    REAL(wp), POINTER           :: p_msg(:)
    CHARACTER(LEN=*), PARAMETER :: subname = MODUL_NAME//'compute_shutdown_restart'

#ifdef DEBUG
    WRITE (nerr,FORMAT_VALS5)subname,' p_pe=',p_pe, &
      & ' call p_barrier with communicator=',p_comm_work
#endif

    ! make sure all are here
    CALL p_barrier(comm=p_comm_work)

    IF(p_pe_work == 0) THEN

      ! create message array
      p_msg => get_message_array(subname)
      p_msg(1) = REAL(MSG_RESTART_SHUTDOWN, wp)

#ifdef DEBUG
      WRITE (nerr,FORMAT_VALS9)subname,' p_pe=',p_pe, &
        & ' send message=',INT(p_msg(1)),' len=',SIZE(p_msg), &
        & ' to pe=',p_restart_pe0
#endif
      CALL p_send(p_msg, p_restart_pe0, 0)

      IF (ASSOCIATED(p_msg)) DEALLOCATE(p_msg)
    ENDIF

  END SUBROUTINE compute_shutdown_restart

  !------------------------------------------------------------------------------------------------
  !
  !  Get a message array to transfer all dynamical restart arguments
  !  between compute and restart PEs.
  !
  FUNCTION get_message_array (subname) 
    REAL(wp), POINTER             :: get_message_array(:)
    CHARACTER(LEN=*), INTENT(in)  :: subname

    INTEGER                       :: ierrstat, n_msg
    TYPE(t_patch_data), POINTER   :: p_pd

    ! set minimum size to transfer patch data
    n_msg = MIN_DYN_RESTART_PDATA

    ! considerate dynamic attributes
    p_pd => patch_data(1)
    n_msg = n_msg + p_pd%n_opt_pvct
    n_msg = n_msg + p_pd%n_opt_lcall_phy
    n_msg = n_msg + p_pd%n_opt_t_elapsed_phy

    ! calculate summary of all data
    n_msg = MIN_DYN_RESTART_ARGS + (SIZE(patch_data) * n_msg)

#ifdef DEBUG
    WRITE (nerr,FORMAT_VALS5)subname,' p_pe=',p_pe, &
      & ' calculated message size=',n_msg
#endif

    ! allocate memory
    ALLOCATE (get_message_array(n_msg), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (subname, ALLOCATE_FAILED)
    get_message_array = 0._wp ! initialize array to avoid Cray compiler warning

  END FUNCTION get_message_array

  !------------------------------------------------------------------------------------------------
  !
  !  Helper function to increment the given integer value.
  !
  FUNCTION incr(ival)

    INTEGER, INTENT(INOUT) :: ival
    INTEGER                :: incr

    ival = ival + 1
    incr = ival

  END FUNCTION incr

  !------------------------------------------------------------------------------------------------
  !
  !  Helper function to convert a given real value into a flag.
  !
  FUNCTION get_flag(rval)

    REAL(wp), INTENT(IN) :: rval
    LOGICAL              :: get_flag

    IF (rval == 1._wp) THEN
      get_flag = .TRUE.
    ELSE
      get_flag = .FALSE.
    ENDIF

  END FUNCTION get_flag

  !------------------------------------------------------------------------------------------------
  !
  ! Common helper routines for error handling and releasing of resources.
  !
  !------------------------------------------------------------------------------------------------
  !
  !  A simple helper routine to check the result of the last MPI call.
  !
  SUBROUTINE check_mpi_error(subname, mpi_call, mpi_error, l_finish)

    CHARACTER(LEN=*), INTENT(IN)    :: subname, mpi_call
    INTEGER, INTENT(IN)             :: mpi_error
    LOGICAL, INTENT(IN)             :: l_finish

    CHARACTER (LEN=MAX_ERROR_LENGTH):: error_message

    IF (mpi_error /= MPI_SUCCESS) THEN
      IF (l_finish) THEN
        WRITE (error_message, '(2a,i5)')TRIM(mpi_call), &
          &                    ' returned with error=',mpi_error
        CALL finish(subname, TRIM(error_message))
      ELSE
        WRITE (error_message, '(4a,i5)')TRIM(subname), ".", TRIM(mpi_call), &
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
  !  Destroy the given CDI GRID handle.
  !
  SUBROUTINE destroy_cdi_grid (iID)

    INTEGER, INTENT(INOUT) :: iID

    IF (iID /= CDI_UNDEFID) THEN
      CALL gridDestroy(iID)
      iID = CDI_UNDEFID
    ENDIF

  END SUBROUTINE destroy_cdi_grid

  !------------------------------------------------------------------------------------------------
  !
  !  Release resources of a restart file.
  !
  SUBROUTINE release_restart_file (rf)

    TYPE(t_restart_file), INTENT(INOUT) :: rf

    INTEGER                             :: i

    CHARACTER(LEN=*), PARAMETER   :: subname = MODUL_NAME//'release_restart_file'

#ifdef DEBUG
    WRITE (nerr,FORMAT_VALS3)subname,' is called for p_pe=',p_pe
#endif

    IF (ALLOCATED(rf%mem_win_off)) DEALLOCATE(rf%mem_win_off)
    IF (ASSOCIATED(rf%var_data))   DEALLOCATE(rf%var_data)

    IF (my_process_is_restart()) THEN

      CALL close_restart_file(rf)

      CALL destroy_cdi_grid(rf%cdiCellGridID)
      CALL destroy_cdi_grid(rf%cdiVertGridID)
      CALL destroy_cdi_grid(rf%cdiEdgeGridID)

      IF (rf%cdiTaxisID /= CDI_UNDEFID) THEN
        CALL taxisDestroy(rf%cdiTaxisID)
        rf%cdiTaxisID = CDI_UNDEFID
      ENDIF

      DO i = 1, SIZE(rf%cdiZaxisIDs)
        IF (rf%cdiZaxisIDs(i) /= CDI_UNDEFID) THEN
          CALL zaxisDestroy(rf%cdiZaxisIDs(i))
          rf%cdiZaxisIDs(i) = CDI_UNDEFID
        ENDIF
      ENDDO

      IF (rf%cdiVlistID /= CDI_UNDEFID) THEN
        CALL vlistDestroy(rf%cdiVlistID)
        rf%cdiVlistID = CDI_UNDEFID
      ENDIF

      rf%cdiTimeIndex = CDI_UNDEFID

    ENDIF

  END SUBROUTINE release_restart_file

  !------------------------------------------------------------------------------------------------
  !
  !  Release all resource of the restart process.
  !
  SUBROUTINE release_resources

    TYPE(t_patch_data), POINTER   :: p_pd
    INTEGER                       :: idx, mpi_error

    CHARACTER(LEN=*), PARAMETER   :: subname = MODUL_NAME//'release_resources'

#ifdef DEBUG
    WRITE (nerr,FORMAT_VALS3)subname,' is called for p_pe=',p_pe
#endif

    ! release patch data
    IF (ALLOCATED(patch_data)) THEN
      DO idx = 1, SIZE(patch_data)

        p_pd => patch_data(idx)

        ! release patch arrays
        IF (ASSOCIATED(p_pd%v_grid_defs))      DEALLOCATE(p_pd%v_grid_defs)
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
      CALL check_mpi_error(subname, 'MPI_Win_fence', mpi_error, .FALSE.)
      CALL MPI_Win_free(mpi_win, mpi_error)
      CALL check_mpi_error(subname, 'MPI_Win_free', mpi_error, .FALSE.)
      mpi_win = MPI_WIN_NULL
    ENDIF

    ! release RMA memory
    CALL MPI_Free_mem(mem_ptr_dp, mpi_error)
    CALL check_mpi_error(subname, 'MPI_Free_mem', mpi_error, .FALSE.)

  END SUBROUTINE release_resources

  !------------------------------------------------------------------------------------------------
  !
  ! Common helper routines to processing lists.
  !

  !------------------------------------------------------------------------------------------------
  !
  ! Prints the restart name lists.
  !
  SUBROUTINE print_restart_name_lists

#ifdef DEBUG
    INTEGER                           :: iv
    CHARACTER(LEN=MAX_NAME_LENGTH)    :: list_name
#ifdef HAVE_F2003
    CHARACTER(LEN=:), ALLOCATABLE     :: list_text
#else
    CHARACTER(LEN=MAX_ATTRIB_TLENGTH) :: list_text
#endif
    CHARACTER(LEN=*), PARAMETER       :: subname = MODUL_NAME//'print_restart_name_list'

    WRITE (nerr,FORMAT_VALS3)subname,' p_pe=',p_pe

    PRINT *,'restart name lists count=',nmls
    DO iv = 1, nmls
      CALL get_restart_namelist(iv, list_name)
      CALL get_restart_namelist(list_name, list_text)

      PRINT *,' restart name list=',TRIM(list_name),' text=',TRIM(list_text)
#ifdef HAVE_F2003
      IF (ALLOCATED(list_text)) DEALLOCATE(list_text)
#endif
    ENDDO
#endif

    END SUBROUTINE print_restart_name_lists

  !------------------------------------------------------------------------------------------------
  !
  !  Print restart arguments.
  !
  SUBROUTINE print_restart_arguments

#ifdef DEBUG
    CHARACTER(LEN=*), PARAMETER   :: subname = MODUL_NAME//'print_restart_arguments'
    TYPE(t_restart_args), POINTER :: p_ra
    TYPE(t_patch_data), POINTER   :: p_pd

    WRITE (nerr,FORMAT_VALS3)subname,' is called for p_pe=',p_pe

    p_ra => restart_args
    PRINT *,subname, ' current_caltime=', p_ra%datetime%caltime
    PRINT *,subname, ' current_calday=',  p_ra%datetime%calday
    PRINT *,subname, ' current_daysec=',  p_ra%datetime%daysec

#ifdef HAVE_F2003
    PRINT *,subname, ' max. attr. text length=', p_ra%n_max_attrib_tlen
#endif

    ! patch informations
    PRINT *,subname, ' size of patches=',        SIZE(patch_data)
    p_pd => patch_data(1)
    PRINT *,subname, ' pd%n_opt_pvct=',          p_pd%n_opt_pvct
    PRINT *,subname, ' pd%n_opt_lcall_phy=',     p_pd%n_opt_lcall_phy
    PRINT *,subname, ' pd%n_opt_t_elapsed_phy=', p_pd%n_opt_t_elapsed_phy
#endif

  END SUBROUTINE print_restart_arguments

  !------------------------------------------------------------------------------------------------
  !
  !  Print restart attributes.
  !
  SUBROUTINE print_restart_attributes

#ifdef DEBUG
    CHARACTER(LEN=MAX_NAME_LENGTH)    :: attrib_name
#ifdef HAVE_F2003
    TYPE(t_restart_args), POINTER     :: p_ra
    INTEGER                           :: ierrstat
    CHARACTER(LEN=:), ALLOCATABLE     :: attrib_txt
#else
    CHARACTER(LEN=MAX_ATTRIB_TLENGTH) :: attrib_txt
#endif

    INTEGER                           :: attrib_int, i
    REAL(wp)                          :: attrib_real
    LOGICAL                           :: attrib_bool
    CHARACTER(LEN=*), PARAMETER       :: subname = MODUL_NAME//'print_restart_attributes'

    WRITE (nerr,FORMAT_VALS3)subname,' is called for p_pe=',p_pe

#ifdef HAVE_F2003
    p_ra => restart_args
    ALLOCATE(CHARACTER(LEN=p_ra%n_max_attrib_tlen) :: attrib_txt, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(subname, ALLOCATE_FAILED)
#endif

    ! check text attributes
    DO i = 1, restart_attributes_count_text()
      CALL get_restart_attribute(i, attrib_name, attrib_txt)
      PRINT *,'restart text attribute: ', TRIM(attrib_name),'=',TRIM(attrib_txt)
    ENDDO

#ifdef HAVE_F2003
    DEALLOCATE(attrib_txt)
#endif

    ! check integer attributes
    DO i= 1, restart_attributes_count_int()
      CALL get_restart_attribute(i, attrib_name, attrib_int)
      PRINT *,'restart integer attribute: ', TRIM(attrib_name),'=',attrib_int
    ENDDO

    ! check real attributes
    DO i = 1, restart_attributes_count_real()
      CALL get_restart_attribute(i, attrib_name, attrib_real)
      PRINT *,'restart real attribute: ', TRIM(attrib_name),'=',attrib_real
    ENDDO

    ! check boolean attributes
    DO i = 1, restart_attributes_count_bool()
      CALL get_restart_attribute(i, attrib_name, attrib_bool)
      PRINT *,'restart booloean attribute: ', TRIM(attrib_name),'=',attrib_bool
    ENDDO
#endif

  END SUBROUTINE print_restart_attributes

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
    CHARACTER(LEN=*), PARAMETER   :: subname = MODUL_NAME//'get_var_list_number'

    WRITE (nerr,FORMAT_VALS5)subname,' p_pe=',p_pe,' patch_id=',patch_id
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

    CHARACTER(LEN=*), PARAMETER :: subname = MODUL_NAME//'transfer_restart_var_lists'

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
      IF (ierrstat /= SUCCESS) CALL finish(subname, ALLOCATE_FAILED)

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
            IF (ierrstat /= SUCCESS) CALL finish(subname, ALLOCATE_FAILED)

            element => p_var_list%p%first_list_element
          ELSE
            ALLOCATE(element%next_list_element, STAT=ierrstat)
            IF (ierrstat /= SUCCESS) CALL finish(subname, ALLOCATE_FAILED)
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
  ! Transfers the restart name lists from the worker to the restart PEs.
  !
  SUBROUTINE transfer_restart_name_lists
    INTEGER                           :: iv, nv
    CHARACTER(LEN=MAX_NAME_LENGTH)    :: list_name
#ifndef HAVE_F2003
    CHARACTER(LEN=MAX_ATTRIB_TLENGTH) :: list_text
#else
    TYPE(t_restart_args),             POINTER :: p_ra
    INTEGER                           :: ierrstat
    CHARACTER(LEN=:), ALLOCATABLE     :: list_text
    CHARACTER(LEN=*), PARAMETER       :: subname = MODUL_NAME//'transfer_restart_name_lists'

#ifdef DEBUG
    WRITE (nerr,FORMAT_VALS3)subname,' p_pe=',p_pe
#endif
#endif

    ! delete old name lists
    IF (my_process_is_restart()) CALL delete_restart_namelists

    ! get the number of name lists
    IF(.NOT. my_process_is_restart()) nv = nmls
    CALL p_bcast(nv, bcast_root, p_comm_work_2_restart)

#ifdef HAVE_F2003
    ! get the maximum text length
    p_ra => restart_args
    p_ra%n_max_attrib_tlen = 0
    IF (my_process_is_work()) THEN
      DO iv = 1, nv
        p_ra%n_max_attrib_tlen = MAX(p_ra%n_max_attrib_tlen, &
          &                          LEN_TRIM(restart_namelist(iv)%text))
      ENDDO

#ifdef DEBUG
      PRINT *,subname,' maximum attribute text length=',p_ra%n_max_attrib_tlen
#endif
    ENDIF

    ! allocate text buffer
    CALL p_bcast(p_ra%n_max_attrib_tlen, bcast_root, p_comm_work_2_restart)
    ALLOCATE(CHARACTER(LEN=p_ra%n_max_attrib_tlen) :: list_text, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(subname, ALLOCATE_FAILED)
#endif

    DO iv = 1, nv
      ! send name of the name list
      list_name = ''
      IF (my_process_is_work()) CALL get_restart_namelist(iv, list_name)
      CALL p_bcast(list_name, bcast_root, p_comm_work_2_restart)

      ! send text of the name list
      list_text = ''
      IF (my_process_is_work()) list_text = TRIM(restart_namelist(iv)%text)
      CALL p_bcast(list_text, bcast_root, p_comm_work_2_restart)

      ! store name list parameters
      IF (my_process_is_restart()) THEN
        CALL set_restart_namelist(list_name, list_text)
      ENDIF
    ENDDO

#ifdef HAVE_F2003
    DEALLOCATE(list_text)
#endif

  END SUBROUTINE transfer_restart_name_lists

  !-------------------------------------------------------------------------------------------------
  !
  ! Create patch data and transfers this data from the worker to the restart PEs.
  !
  SUBROUTINE create_and_transfer_patch_data(opt_lcall_phy_size, &
    &                                       opt_t_elapsed_phy_size, &
    &                                       opt_pvct_size)

    INTEGER,  INTENT(IN), OPTIONAL :: opt_lcall_phy_size, opt_t_elapsed_phy_size, &
      &                               opt_pvct_size

    INTEGER                        :: jg, jl, ierrstat
    TYPE(t_patch_data), POINTER    :: p_pd

    CHARACTER(LEN=*), PARAMETER :: subname = MODUL_NAME//'create_and_transfer_patch_data'

#ifdef DEBUG
    WRITE (nerr,FORMAT_VALS3)subname,' p_pe=',p_pe
#endif

    ! replicate domain setup
    CALL p_bcast(n_dom, bcast_root, p_comm_work_2_restart)

    ! allocate patch data structure
    ALLOCATE(patch_data(n_dom), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(subname, ALLOCATE_FAILED)

    ! set number of global cells/edges/verts and patch ID
    DO jg = 1, n_dom

      p_pd => patch_data(jg)
      IF (my_process_is_work()) THEN
        p_pd%id              = p_patch(jg)%id
        p_pd%l_dom_active    = p_patch(jg)%ldom_active
        p_pd%nlev            = p_patch(jg)%nlev
        p_pd%cell_type       = p_patch(jg)%cell_type
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
      CALL p_bcast(p_pd%l_dom_active,    bcast_root, p_comm_work_2_restart)
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
                              p_patch(jl)%cells%decomp_info%owner_mask, p_patch(jl)%cells%phys_id,    &
                              p_patch(jl)%cells%decomp_info%glb_index, p_pd%cells)

        CALL set_reorder_data(jg, p_patch(jl)%n_patch_edges_g, p_patch(jl)%n_patch_edges, &
                              p_patch(jl)%edges%decomp_info%owner_mask, p_patch(jl)%edges%phys_id,    &
                              p_patch(jl)%edges%decomp_info%glb_index, p_pd%edges)

        CALL set_reorder_data(jg, p_patch(jl)%n_patch_verts_g, p_patch(jl)%n_patch_verts, &
                              p_patch(jl)%verts%decomp_info%owner_mask, p_patch(jl)%verts%phys_id,    &
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

      ! init. optional parameter arrays
      p_pd%n_opt_lcall_phy = 0
      IF (my_process_is_work() .AND. PRESENT(opt_lcall_phy_size)) THEN
        p_pd%n_opt_lcall_phy = opt_lcall_phy_size
      ENDIF
      CALL p_bcast(p_pd%n_opt_lcall_phy, bcast_root, p_comm_work_2_restart)

      p_pd%n_opt_t_elapsed_phy = 0
      IF (my_process_is_work() .AND. PRESENT(opt_t_elapsed_phy_size)) THEN
        p_pd%n_opt_t_elapsed_phy = opt_t_elapsed_phy_size
      ENDIF
      CALL p_bcast(p_pd%n_opt_t_elapsed_phy, bcast_root, p_comm_work_2_restart)

      p_pd%n_opt_pvct = 0
      IF (my_process_is_work() .AND. PRESENT(opt_pvct_size)) THEN
        p_pd%n_opt_pvct = opt_pvct_size
      ENDIF
      CALL p_bcast(p_pd%n_opt_pvct, bcast_root, p_comm_work_2_restart)

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

    CHARACTER(LEN=*), PARAMETER           :: subname = MODUL_NAME//'set_restart_file_data'

#ifdef DEBUG
    WRITE (nerr,FORMAT_VALS5)subname,' is called for p_pe=',p_pe,' patch_id=',patch_id
#endif

    ! init. main variables
    rf%my_mem_win_off             = 0_i8
    rf%var_data                   => NULL()
    rf%filename                   = ''
    rf%linkname                   = ''
    rf%linkprefix                 = ''

    rf%cdiFileID                  = CDI_UNDEFID
    rf%cdiVlistID                 = CDI_UNDEFID
    rf%cdiTaxisID                 = CDI_UNDEFID

    rf%cdiCellGridID              = CDI_UNDEFID
    rf%cdiVertGridID              = CDI_UNDEFID
    rf%cdiEdgeGridID              = CDI_UNDEFID
    rf%cdiZaxisIDs(:)             = CDI_UNDEFID
    rf%cdiTimeIndex               = CDI_UNDEFID

    ! counts number of restart variables for this file (logical patch ident)
    num_vars = 0
    CALL get_var_list_number(num_vars, patch_id)

#ifdef DEBUG
    WRITE (nerr,FORMAT_VALS3)subname,' numvars=',num_vars
#endif

    IF (num_vars <= 0) RETURN
    ! allocate the array of restart variables
    ALLOCATE (rf%var_data(num_vars), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(subname, ALLOCATE_FAILED)

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
  SUBROUTINE set_reorder_data(phys_patch_id, n_points_g, n_points, owner_mask, phys_id, &
                              glb_index, reo)

    INTEGER, INTENT(IN) :: phys_patch_id   ! Physical patch ID
    INTEGER, INTENT(IN) :: n_points_g      ! Global number of cells/edges/verts in logical patch
    INTEGER, INTENT(IN) :: n_points        ! Local number of cells/edges/verts in logical patch
    LOGICAL, INTENT(IN) :: owner_mask(:,:) ! owner_mask for logical patch
    INTEGER, INTENT(IN) :: phys_id(:,:)    ! phys_id for logical patch
    INTEGER, INTENT(IN) :: glb_index(:)    ! glb_index for logical patch

    TYPE(t_reorder_data), INTENT(INOUT) :: reo ! Result: reorder data

    INTEGER :: i, n, il, ib, mpi_error, ierrstat
    LOGICAL, ALLOCATABLE :: phys_owner_mask(:) ! owner mask for physical patch
    INTEGER, ALLOCATABLE :: glbidx_own(:), glbidx_glb(:), reorder_index_log_dom(:)
    CHARACTER (LEN=MAX_ERROR_LENGTH) :: error_message

    CHARACTER(LEN=*), PARAMETER :: subname = MODUL_NAME//'set_reorder_data'

#ifdef DEBUG
    WRITE (nerr,FORMAT_VALS3)subname,' is called for p_pe=',p_pe
#endif

    ! just for safety
    IF(my_process_is_restart()) CALL finish(subname, NO_COMPUTE_PE)

    ! set the physical patch owner mask
    ALLOCATE(phys_owner_mask(n_points))
    DO i = 1, n_points
      il = idx_no(i)
      ib = blk_no(i)
      phys_owner_mask(i) = owner_mask(il,ib)
      phys_owner_mask(i) = phys_owner_mask(i) .AND. (phys_id(il,ib) == phys_patch_id)
    ENDDO

    ! get number of owned cells/edges/verts (without halos, physical patch only)
    reo%n_own = COUNT(phys_owner_mask(:))

    ! set index arrays to own cells/edges/verts
    ALLOCATE(reo%own_idx(reo%n_own), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(subname, ALLOCATE_FAILED)
    ALLOCATE(reo%own_blk(reo%n_own), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(subname, ALLOCATE_FAILED)

    ! global index of my own points
    ALLOCATE(glbidx_own(reo%n_own), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(subname, ALLOCATE_FAILED)

    n = 0
    DO i = 1, n_points
      IF(phys_owner_mask(i)) THEN
        n = n+1
        reo%own_idx(n) = idx_no(i)
        reo%own_blk(n) = blk_no(i)
        glbidx_own(n)  = glb_index(i)
      ENDIF
    ENDDO

    ! gather the number of own points for every PE into reo%pe_own
    ALLOCATE(reo%pe_own(0:p_n_work-1), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(subname, ALLOCATE_FAILED)
    ALLOCATE(reo%pe_off(0:p_n_work-1), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(subname, ALLOCATE_FAILED)

    CALL MPI_Allgather(reo%n_own,  1, p_int, &
                       reo%pe_own, 1, p_int, &
                       p_comm_work, mpi_error)
    CALL check_mpi_error(subname, 'MPI_Allgather', mpi_error, .TRUE.)

    ! get offset within result array
    reo%pe_off(0) = 0
    DO i = 1, p_n_work-1
      reo%pe_off(i) = reo%pe_off(i-1) + reo%pe_own(i-1)
    ENDDO

    ! get global number of points for current (physical!) patch
    reo%n_glb = SUM(reo%pe_own(:))

    ! Get the global index numbers of the data when it is gathered on PE 0
    ! exactly in the same order as it is retrieved later during restart.
    ALLOCATE(glbidx_glb(reo%n_glb), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(subname, ALLOCATE_FAILED)

    CALL MPI_Allgatherv(glbidx_own, reo%n_own, p_int, &
                        glbidx_glb, reo%pe_own, reo%pe_off, p_int, &
                        p_comm_work, mpi_error)
    CALL check_mpi_error(subname, 'MPI_Allgatherv', mpi_error, .TRUE.)

    ! get reorder_index
    ALLOCATE(reo%reorder_index(reo%n_glb), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(subname, ALLOCATE_FAILED)

    ! spans the complete logical domain
    ALLOCATE(reorder_index_log_dom(n_points_g), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(subname, ALLOCATE_FAILED)
    reorder_index_log_dom(:) = 0

    DO i = 1, reo%n_glb
      ! reorder_index_log_dom stores where a global point in logical domain comes from.
      ! It is nonzero only at the physical patch locations
      reorder_index_log_dom(glbidx_glb(i)) = i
    ENDDO

    ! gather the reorder index for the physical domain
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
      CALL finish(subname,TRIM(error_message))
    ENDIF

    DEALLOCATE(phys_owner_mask)
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

    CHARACTER(LEN=*), PARAMETER :: subname = MODUL_NAME//'transfer_reorder_data'

    ! transfer the global number of points, this is not yet known on restart PEs
    CALL p_bcast(reo%n_glb,  bcast_root, p_comm_work_2_restart)

    IF(my_process_is_restart()) THEN

      ! on restart PEs: n_own = 0, own_idx and own_blk are not allocated
      reo%n_own = 0

      ! pe_own must be allocated for num_work_procs, not for p_n_work
      ALLOCATE(reo%pe_own(0:num_work_procs-1), STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(subname, ALLOCATE_FAILED)

      ALLOCATE(reo%pe_off(0:num_work_procs-1), STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(subname, ALLOCATE_FAILED)

      ALLOCATE(reo%reorder_index(reo%n_glb), STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(subname, ALLOCATE_FAILED)
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

    CHARACTER(LEN=*), PARAMETER :: subname = MODUL_NAME//'init_remote_memory_access'

#ifdef DEBUG
    WRITE (nerr,FORMAT_VALS3)subname,' is called for p_pe=',p_pe
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
            CALL finish(subname,UNKNOWN_GRID_TYPE)
        END SELECT

      ENDDO

      ! get the offset on all PEs
      ALLOCATE(patch_data(i)%restart_file%mem_win_off(0:num_work_procs-1), STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(subname, ALLOCATE_FAILED)
      IF(.NOT.my_process_is_restart()) THEN
        CALL MPI_Allgather(patch_data(i)%restart_file%my_mem_win_off, 1, p_int_i8, &
                           patch_data(i)%restart_file%mem_win_off, 1, p_int_i8,    &
                           p_comm_work, mpi_error)
        CALL check_mpi_error(subname, 'MPI_Allgather', mpi_error, .TRUE.)
      ENDIF

      CALL p_bcast(patch_data(i)%restart_file%mem_win_off, bcast_root, p_comm_work_2_restart)

    ENDDO

    ! mem_size is calculated as number of variables above, get number of bytes
    ! get the amount of bytes per REAL*8 variable (as used in MPI communication)
    CALL MPI_Type_extent(p_real_dp, nbytes_real, mpi_error)
    CALL check_mpi_error(subname, 'MPI_Type_extent', mpi_error, .TRUE.)

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
    CALL check_mpi_error(subname, 'MPI_Alloc_mem', mpi_error, .TRUE.)

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
     & CALL finish(subname,'c_intptr_t /= MPI_ADDRESS_KIND, too dangerous to proceed!')

    CALL MPI_Alloc_mem(mem_bytes, MPI_INFO_NULL, c_mem_ptr, mpi_error)
    CALL check_mpi_error(subname, 'MPI_Alloc_mem', mpi_error, .TRUE.)

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
    CALL check_mpi_error(subname, 'MPI_Info_create', mpi_error, .TRUE.)
    CALL MPI_Info_set(rma_cache_hint, "IBM_win_cache","0", mpi_error)
    CALL check_mpi_error(subname, 'MPI_Info_set', mpi_error, .TRUE.)
#endif

    ! create memory window for communication
    mem_ptr_dp(:) = 0._dp
    CALL MPI_Win_create(mem_ptr_dp,mem_bytes,nbytes_real,MPI_INFO_NULL,&
      &                 p_comm_work_restart,mpi_win,mpi_error )
    CALL check_mpi_error(subname, 'MPI_Win_create', mpi_error, .TRUE.)

#ifdef __xlC__
    CALL MPI_Info_free(rma_cache_hint, mpi_error);
    CALL check_mpi_error(subname, 'MPI_Info_free', mpi_error, .TRUE.)
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
  FUNCTION find_patch(id, subname)

    INTEGER, INTENT(IN)             :: id
    CHARACTER(LEN=*), INTENT(IN)    :: subname

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
      CALL finish(subname, err_message)
    ENDIF
  END FUNCTION find_patch

  !------------------------------------------------------------------------------------------------
  !
  !  Set vertical grid definition.
  !
  SUBROUTINE set_vertical_grid_def(vg, ntype, nlevels)

    TYPE(t_v_grid), INTENT(INOUT) :: vg
    INTEGER, INTENT(IN)           :: ntype, nlevels

    vg%type    = ntype
    vg%nlevels = nlevels

  END SUBROUTINE set_vertical_grid_def

  !------------------------------------------------------------------------------------------------
  !
  !  Set global restart attributes.
  !
  SUBROUTINE set_restart_attributes (p_pd)

    TYPE(t_patch_data), POINTER, INTENT(INOUT) :: p_pd

    TYPE(t_restart_args),  POINTER :: p_ra
    CHARACTER(LEN=MAX_NAME_LENGTH) :: attrib_name
    INTEGER                        :: jp, jp_end, jg, nlev_soil, &
      &                               nlev_snow, nlev_ocean, nice_class, ierrstat, &
      &                               nlena, nlenb, nlenc, nlend

    CHARACTER(LEN=*), PARAMETER    :: subname = MODUL_NAME//'set_restart_attributes'
    CHARACTER(LEN=*), PARAMETER    :: attrib_format_int  = '(a,i2.2)'
    CHARACTER(LEN=*), PARAMETER    :: attrib_format_int2 = '(a,i2.2,a,i2.2)'

    CHARACTER(LEN=256) :: executable, user_name, os_name, host_name, tmp_string
    CHARACTER(LEN=  8) :: date_string
    CHARACTER(LEN= 10) :: time_string

#ifdef DEBUG
    WRITE (nerr,FORMAT_VALS3)subname,' is called for p_pe=',p_pe
#endif

    ! delete old attributes
    CALL delete_attributes()

    ! get environment attributes
    CALL get_command_argument(0, executable, nlend)
    CALL date_and_time(date_string, time_string)

    tmp_string = ''
    CALL util_os_system (tmp_string, nlena)
    os_name = tmp_string(1:nlena)

    tmp_string = ''
    CALL util_user_name (tmp_string, nlenb)
    user_name = tmp_string(1:nlenb)

    tmp_string = ''
    CALL util_node_name (tmp_string, nlenc)
    host_name = tmp_string(1:nlenc)

    ! set CD-Convention required restart attributes
    CALL set_restart_attribute('title',       &
         MODEL_TITLE)
    CALL set_restart_attribute('institution', &
         MODEL_INSTITUTION)
    CALL set_restart_attribute('source',      &
         TRIM(out_expname)//'-'//MODEL_VERSION)
    CALL set_restart_attribute('history',     &
         executable(1:nlend)//' at '//date_string(1:8)//' '//time_string(1:6))
    CALL set_restart_attribute('references',  &
         MODEL_REFERENCES)
    CALL set_restart_attribute('comment',     &
         TRIM(user_name)//' on '//TRIM(host_name)//' ('//TRIM(os_name)//')')

    ! set restart time
    p_ra => restart_args
    CALL set_restart_attribute ('current_caltime', p_ra%datetime%caltime)
    CALL set_restart_attribute ('current_calday' , p_ra%datetime%calday)
    CALL set_restart_attribute ('current_daysec' , p_ra%datetime%daysec)

    ! set simulation step
    CALL set_restart_attribute( 'jstep', p_ra%jstep )

    ! set time levels
    jg = p_pd%id
    CALL set_restart_attribute( 'nold'    , p_pd%nold)
    CALL set_restart_attribute( 'nnow'    , p_pd%nnow)
    CALL set_restart_attribute( 'nnew'    , p_pd%nnew)
    CALL set_restart_attribute( 'nnow_rcf', p_pd%nnow_rcf)
    CALL set_restart_attribute( 'nnew_rcf', p_pd%nnew_rcf)

    ! additional restart-output for nonhydrostatic model
    IF (p_pd%l_opt_sim_time) THEN
      WRITE(attrib_name, attrib_format_int) 'sim_time_DOM', jg
      CALL set_restart_attribute (TRIM(attrib_name), &
        &                         p_pd%opt_sim_time)
    ENDIF

    !-------------------------------------------------------------
    ! DR
    ! WORKAROUND FOR FIELDS WHICH NEED TO GO INTO THE RESTART FILE,
    ! BUT SO FAR CANNOT BE HANDELED CORRECTLY BY ADD_VAR OR
    ! SET_RESTART_ATTRIBUTE
    !-------------------------------------------------------------
    IF (p_pd%l_opt_jstep_adv_ntstep) THEN
        WRITE(attrib_name, attrib_format_int) 'jstep_adv_ntsteps_DOM', jg
        CALL set_restart_attribute (TRIM(attrib_name), &
          &                         p_pd%opt_jstep_adv_ntstep)
    ENDIF

    IF (p_pd%l_opt_jstep_adv_marchuk_order) THEN
        WRITE(attrib_name, attrib_format_int) 'jstep_adv_marchuk_order_DOM', jg
        CALL set_restart_attribute (TRIM(attrib_name), &
          &                         p_pd%opt_jstep_adv_marchuk_order)
    ENDIF

    IF (ALLOCATED(p_pd%opt_t_elapsed_phy) .AND. &
      & ALLOCATED(p_pd%opt_lcall_phy)) THEN
      jp_end = SIZE(p_pd%opt_t_elapsed_phy)
      DO jp = 1, jp_end
        WRITE(attrib_name, attrib_format_int2) 't_elapsed_phy_DOM',jg,'_PHY',jp
        CALL set_restart_attribute (TRIM(attrib_name), p_pd%opt_t_elapsed_phy(jp))
      ENDDO
      jp_end = SIZE(p_pd%opt_lcall_phy)
      DO jp = 1, jp_end
        WRITE(attrib_name, attrib_format_int2) 'lcall_phy_DOM',jg,'_PHY', jp
        CALL set_restart_attribute (TRIM(attrib_name), p_pd%opt_lcall_phy(jp) )
      ENDDO
    ENDIF

    ! geometrical depth for land module
    IF (p_pd%l_opt_depth_lnd) THEN
      nlev_soil = p_pd%opt_depth_lnd
    ELSE
      nlev_soil = 0
    ENDIF

    ! number of snow levels (multi layer snow model)
    IF (p_pd%l_opt_nlev_snow) THEN
      nlev_snow = p_pd%opt_nlev_snow
    ELSE
      nlev_snow = 0
    ENDIF

    ! ocean depth
    IF (p_pd%l_opt_depth) THEN
      nlev_ocean = p_pd%opt_depth
    ELSE
      nlev_ocean = 0
    END IF

    IF (p_pd%l_opt_nice_class) THEN
      nice_class = p_pd%opt_nice_class
    ELSE
      nice_class = 1
    END IF

    ! set vertical grid definitions
    ALLOCATE(p_pd%v_grid_defs(MAX_VERTICAL_AXES), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(subname, ALLOCATE_FAILED)

    CALL set_vertical_grid_def(p_pd%v_grid_defs(1),  ZA_SURFACE             , 1            )
    CALL set_vertical_grid_def(p_pd%v_grid_defs(2),  ZA_HYBRID              , p_pd%nlev    )
    CALL set_vertical_grid_def(p_pd%v_grid_defs(3),  ZA_HYBRID_HALF         , p_pd%nlev+1  )
    CALL set_vertical_grid_def(p_pd%v_grid_defs(4),  ZA_DEPTH_BELOW_LAND    , nlev_soil    )
    CALL set_vertical_grid_def(p_pd%v_grid_defs(5),  ZA_DEPTH_BELOW_LAND_P1 , nlev_soil+1  )
    CALL set_vertical_grid_def(p_pd%v_grid_defs(6),  ZA_SNOW                , nlev_snow    )
    CALL set_vertical_grid_def(p_pd%v_grid_defs(7),  ZA_SNOW_HALF           , nlev_snow+1  )
    CALL set_vertical_grid_def(p_pd%v_grid_defs(8),  ZA_HEIGHT_2M           , 1            )
    CALL set_vertical_grid_def(p_pd%v_grid_defs(9),  ZA_HEIGHT_10M          , 1            )
    CALL set_vertical_grid_def(p_pd%v_grid_defs(10), ZA_TOA                 , 1            )
    CALL set_vertical_grid_def(p_pd%v_grid_defs(11), ZA_LAKE_BOTTOM         , 1            )
    CALL set_vertical_grid_def(p_pd%v_grid_defs(12), ZA_MIX_LAYER           , 1            )
    CALL set_vertical_grid_def(p_pd%v_grid_defs(13), ZA_LAKE_BOTTOM_HALF    , 1            )
    CALL set_vertical_grid_def(p_pd%v_grid_defs(14), ZA_SEDIMENT_BOTTOM_TW_HALF, 1         )
    CALL set_vertical_grid_def(p_pd%v_grid_defs(15), ZA_DEPTH_BELOW_SEA     , nlev_ocean   )
    CALL set_vertical_grid_def(p_pd%v_grid_defs(16), ZA_DEPTH_BELOW_SEA_HALF, nlev_ocean+1 )
    CALL set_vertical_grid_def(p_pd%v_grid_defs(17), ZA_GENERIC_ICE         , nice_class   )
    CALL set_vertical_grid_def(p_pd%v_grid_defs(18), ZA_DEPTH_RUNOFF_S      , 1            )
    CALL set_vertical_grid_def(p_pd%v_grid_defs(19), ZA_DEPTH_RUNOFF_G      , 1            )

  END SUBROUTINE set_restart_attributes

  !------------------------------------------------------------------------------------------------
  !
  ! Common helper routines to write a restart file.
  !
  !------------------------------------------------------------------------------------------------
  !
  ! Returns true, if the time level of the given field is valid, else false.
  !
  FUNCTION has_valid_time_level (p_info)

    TYPE(t_var_metadata), POINTER, INTENT(IN) :: p_info
    LOGICAL                       :: has_valid_time_level

    INTEGER                       :: idx, time_level, tlev_skip

#ifdef DEBUG
    CHARACTER(LEN=*), PARAMETER   :: subname = MODUL_NAME//'has_valid_time_level'
#endif

    has_valid_time_level = .TRUE.
    IF (.NOT. p_info%lrestart) THEN
      has_valid_time_level = .FALSE.
      RETURN
    ENDIF

    ! get time index of the given field
    idx = INDEX(p_info%name,'.TL')
    IF (idx == 0) THEN
      time_level = -1
    ELSE
      time_level = ICHAR(p_info%name(idx+3:idx+3)) - ICHAR('0')
    ENDIF

    ! get information about time level to be skipped for current field
    ! for the time being this will work with the global patch only
    IF (p_info%tlev_source == 0) THEN
      tlev_skip = nnew(1)          ! ATTENTION: 1 (global patch) hardcoded
    ELSE IF (p_info%tlev_source == 1) THEN
      tlev_skip = nnew_rcf(1)      ! ATTENTION: 1 (global patch) hardcoded
    ELSE
      tlev_skip = -99
    ENDIF

    SELECT CASE (iequations)
      CASE(IHS_ATM_TEMP, IHS_ATM_THETA, ISHALLOW_WATER)

        IF ( time_level == tlev_skip                        &
          & .AND. ha_dyn_config%itime_scheme/=LEAPFROG_EXPL &
          & .AND. ha_dyn_config%itime_scheme/=LEAPFROG_SI   ) &
          & has_valid_time_level = .FALSE.
      CASE default
        IF ( time_level == tlev_skip ) has_valid_time_level = .FALSE.
    END SELECT

#ifdef DEBUG
    IF (.NOT. has_valid_time_level) THEN
      WRITE (nerr,'(2a,i3,a,i4,3a)')subname,' p_pe=',p_pe,' time level=',time_level, &
        &                           ' of field=',TRIM(p_info%name),' is invalid'
    ENDIF
#endif

  END FUNCTION has_valid_time_level

  !------------------------------------------------------------------------------------------------
  !
  ! Returns the pointer of the reorder data for the given field.
  !
  FUNCTION get_reorder_ptr (p_pd, p_info,subname)

    TYPE(t_patch_data), POINTER, INTENT(IN)   :: p_pd
    TYPE(t_var_metadata), POINTER, INTENT(IN) :: p_info
    CHARACTER(LEN=*), INTENT(IN)              :: subname

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
        CALL finish(subname, UNKNOWN_GRID_TYPE)
    END SELECT

  END FUNCTION get_reorder_ptr

  !------------------------------------------------------------------------------------------------
  !
  ! Check the status of the last netCDF file operation.
  !
  SUBROUTINE check_netcdf_status(status, subname)

    INTEGER, INTENT(IN)           :: status
    CHARACTER(LEN=*), INTENT(IN)  :: subname

    CHARACTER(LEN=128)            :: error_text


    IF (status /= nf_noerr) THEN
      WRITE (error_text, NET_CDF_ERROR_FORMAT)'netCDF error ', status , &
        &                                     ' - ' // nf_strerror(status)
      CALL finish(subname, error_text)
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

    INTEGER                         :: iv, nval, nlev_max, ierrstat, nlevs, nv_off, &
      &                                np, mpi_error, jk, i, idate, itime, status
    INTEGER(KIND=MPI_ADDRESS_KIND)  :: ioff(0:num_work_procs-1)
    INTEGER                         :: voff(0:num_work_procs-1)
    REAL(dp), ALLOCATABLE           :: var1_dp(:), var2_dp(:), var3_dp(:,:)

    CHARACTER(LEN=*), PARAMETER     :: subname = MODUL_NAME//'restart_write_var_list'

#ifdef DEBUG
    WRITE (nerr,FORMAT_VALS3)subname,' p_pe=',p_pe
#endif

    ! check process
    IF (.NOT. my_process_is_restart()) CALL finish(subname, NO_RESTART_PE)

    ! write restart time
    dt => restart_args%datetime
    idate = cdiEncodeDate(dt%year, dt%month, dt%day)
    itime = cdiEncodeTime(dt%hour, dt%minute, NINT(dt%second))

    p_rf => p_pd%restart_file
    CALL taxisDefVdate(p_rf%cdiTaxisID, idate)
    CALL taxisDefVtime(p_rf%cdiTaxisID, itime)
    status = streamDefTimestep(p_rf%cdiFileID, p_rf%cdiTimeIndex)

    p_rf%cdiTimeIndex = p_rf%cdiTimeIndex + 1

    ! check the contained array of restart variables
    p_vars => p_pd%restart_file%var_data
    IF (.NOT. ASSOCIATED(p_vars)) RETURN

    ! get maximum number of data points in a slice and allocate tmp. variables
    nval = MAX(p_pd%cells%n_glb, p_pd%edges%n_glb, p_pd%verts%n_glb)

    ! get max. level
    nlev_max = 1
    DO iv = 1, SIZE(p_vars)
      p_info => p_vars(iv)%info
      IF(p_info%ndims == 3) nlev_max = MAX(nlev_max, p_info%used_dimensions(2))
    ENDDO

    ! allocate RMA memory
    ALLOCATE(var1_dp(nval*nlev_max), var2_dp(-1:nval), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (subname, ALLOCATE_FAILED)

    ioff(:) = p_rf%mem_win_off(:)

    ! go over the all restart variables in the associated array
    DO iv = 1, SIZE(p_vars)

      ! get pointer to metadata
      p_info => p_vars(iv)%info

#ifdef DEBUG
      WRITE (nerr,FORMAT_VALS5I)subname,' p_pe=',p_pe,' restart pe processes field=',TRIM(p_info%name)
#endif

      ! check time level of the field
      IF (.NOT. has_valid_time_level(p_info)) CYCLE

      ! get current level
      IF(p_info%ndims == 2) THEN
        nlevs = 1
      ELSE
        nlevs = p_info%used_dimensions(2)
      ENDIF

      ! get pointer to reorder data
      p_ri => get_reorder_ptr (p_pd, p_info, subname)

      ! retrieve part of variable from every worker PE using MPI_Get
      nv_off = 0
      DO np = 0, num_work_procs-1

        IF(p_ri%pe_own(np) == 0) CYCLE

        ! number of words to transfer
        nval = p_ri%pe_own(np)*nlevs

        CALL MPI_Win_lock(MPI_LOCK_SHARED, np, MPI_MODE_NOCHECK, mpi_win, mpi_error)
        CALL check_mpi_error(subname, 'MPI_Win_lock', mpi_error, .TRUE.)

        CALL MPI_Get(var1_dp(nv_off+1), nval, p_real_dp, np, ioff(np), &
          &          nval, p_real_dp, mpi_win, mpi_error)
        CALL check_mpi_error(subname, 'MPI_Get', mpi_error, .TRUE.)

        CALL MPI_Win_unlock(np, mpi_win, mpi_error)
        CALL check_mpi_error(subname, 'MPI_Win_unlock', mpi_error, .TRUE.)

        ! update the offset in var1
        nv_off = nv_off + nval

        ! update the offset in the RMA window on compute PEs
        ioff(np) = ioff(np) + INT(nval,i8)

      ENDDO

      ! compute the total offset for each PE
      nv_off = 0
      DO np = 0, num_work_procs-1
        voff(np) = nv_off
        nval     = p_ri%pe_own(np)*nlevs
        nv_off   = nv_off + nval
      END DO

      ! var1 is stored in the order in which the variable was stored on compute PEs,
      ! get it back into the global storage order
      ALLOCATE(var3_dp(p_ri%n_glb, nlevs), STAT=ierrstat) ! Must be allocated to exact size
      IF (ierrstat /= SUCCESS) CALL finish (subname, ALLOCATE_FAILED)

      ! go over all levels
      DO jk = 1, nlevs
        nv_off = 0
        DO np = 0, num_work_procs-1
          var2_dp(-1) = 0._dp ! special value for lon-lat areas overlapping local patches
          var2_dp(nv_off+1:nv_off+p_ri%pe_own(np)) = var1_dp(voff(np)+1:voff(np)+p_ri%pe_own(np))
          nv_off = nv_off+p_ri%pe_own(np)
          voff(np) = voff(np)+p_ri%pe_own(np)
        ENDDO
        DO i = 1, p_ri%n_glb
          var3_dp(i,jk) = var2_dp(p_ri%reorder_index(i))
        ENDDO
      ENDDO

#ifdef DEBUG
      WRITE (nerr,FORMAT_VALS7I)subname,' p_pe=',p_pe,' restart pe writes field=', &
        &                               TRIM(p_info%name),' data=',p_ri%n_glb*nlevs
#endif

      ! write field content into a file
      CALL streamWriteVar(p_rf%cdiFileID, p_info%cdiVarID, var3_dp, 0)

      DEALLOCATE(var3_dp, STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish (subname, DEALLOCATE_FAILED)

    ENDDO

    DEALLOCATE(var1_dp, var2_dp, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (subname, DEALLOCATE_FAILED)

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
    INTEGER                         :: iv, mpi_error, nindex, ierrstat, nlevs, i, jk
    INTEGER(i8)                     :: ioff
    CHARACTER(LEN=*), PARAMETER     :: subname = MODUL_NAME//'compute_write_var_list'

#ifdef DEBUG
    WRITE (nerr,FORMAT_VALS3)subname,' p_pe=',p_pe
#endif

    ! check process
    IF (.NOT. my_process_is_work()) CALL finish(subname, NO_COMPUTE_PE)

    ! check the array of restart variables
    p_vars => p_pd%restart_file%var_data
    IF (.NOT. ASSOCIATED(p_vars)) RETURN

    ! offset in RMA window for async restart
    p_rf => p_pd%restart_file
    ioff = p_rf%my_mem_win_off

    ! in case of async restart: Lock own window before writing to it
    CALL MPI_Win_lock(MPI_LOCK_EXCLUSIVE, p_pe_work, MPI_MODE_NOCHECK, mpi_win, mpi_error)
    CALL check_mpi_error(subname, 'MPI_Win_lock', mpi_error, .TRUE.)

    ! go over the all restart variables in the associated array
    DO iv = 1, SIZE(p_vars)

      ! get pointer to metadata
      p_info => p_vars(iv)%info

#ifdef DEBUG
      WRITE (nerr,FORMAT_VALS5I)subname,' p_pe=',p_pe,' compute pe processes field=',TRIM(p_info%name)
#endif

      ! check time level of the field
      IF (.NOT. has_valid_time_level(p_info)) CYCLE

      ! Check if first dimension of array is nproma.
      ! Otherwise we got an array which is not suitable for this output scheme.
      IF (p_info%used_dimensions(1) /= nproma) &
        CALL finish(subname,'1st dim is not nproma: '//TRIM(p_info%name))

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
          CALL message(subname, p_info%name)
          CALL finish(subname,'1d arrays not handled yet.')
        CASE (2)
          ! make a 3D copy of the array
          ALLOCATE(r_ptr(p_info%used_dimensions(1),1,p_info%used_dimensions(2)), &
            &      STAT=ierrstat)
          IF (ierrstat /= SUCCESS) CALL finish(subname, ALLOCATE_FAILED)
          r_ptr(:,1,:) = p_vars(iv)%r_ptr(:,:,nindex,1,1)
        CASE (3)
          ! copy the pointer
          r_ptr => p_vars(iv)%r_ptr(:,:,:,nindex,1)
        CASE (4)
          CALL message(subname, p_info%name)
          CALL finish(subname,'4d arrays not handled yet.')
        CASE (5)
          CALL message(subname, p_info%name)
          CALL finish(subname,'5d arrays not handled yet.')
        CASE DEFAULT
          CALL message(subname, p_info%name)
          CALL finish(subname,'dimension not set.')
      END SELECT

      ! get number of data levels
      IF(p_info%ndims == 2) THEN
        nlevs = 1
      ELSE
        nlevs = p_info%used_dimensions(2)
      ENDIF

      ! get pointer to reorder data
      p_ri => get_reorder_ptr(p_pd, p_info, subname)

#ifdef DEBUG
      WRITE (nerr,FORMAT_VALS7I)subname,' p_pe=',p_pe,' compute pe writes field=', &
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
    CALL check_mpi_error(subname, 'MPI_Win_unlock', mpi_error, .TRUE.)

  END SUBROUTINE compute_write_var_list

  !------------------------------------------------------------------------------------------------
  !
  ! Creates a horizontal grid definition from the given parameters.
  !
  SUBROUTINE create_cdi_hgrid_def(iID, iCnt, iNVert, cNameX, cLNameX, cUnitsX, &
    &                                                cNameY, cLNameY, cUnitsY)

    INTEGER, INTENT(INOUT)        :: iID
    INTEGER, INTENT(IN)           :: iCnt, iNVert
    CHARACTER(LEN=*), INTENT(IN)  :: cNameX, cLNameX, cUnitsX, cNameY, cLNameY, cUnitsY

    iID = gridCreate      (GRID_UNSTRUCTURED, iCnt)
    CALL gridDefNvertex   (iID, iNVert)

    CALL gridDefXname     (iID, TRIM(cNameX))
    CALL gridDefXlongname (iID, TRIM(cLNameX))
    CALL gridDefXunits    (iID, TRIM(cUnitsX))

    CALL gridDefYname     (iID, TRIM(cNameY))
    CALL gridDefYlongname (iID, TRIM(cLNameY))
    CALL gridDefYunits    (iID, TRIM(cUnitsY))

  END SUBROUTINE create_cdi_hgrid_def

  !------------------------------------------------------------------------------------------------
  !
  ! Creates a Z-axis grid definition from the given parameters.
  !
  SUBROUTINE create_cdi_zaxis(iID, iGridID, iLevels, iDefLevels, rDefLevelVal, lOcean)

    INTEGER,  INTENT(INOUT)        :: iID
    INTEGER,  INTENT(IN)           :: iGridID, iLevels
    INTEGER,  INTENT(IN), OPTIONAL :: iDefLevels
    REAL(wp), INTENT(IN), OPTIONAL :: rDefLevelVal
    LOGICAL,  INTENT(IN), OPTIONAL :: lOcean

    REAL(wp), ALLOCATABLE          :: rDefLevelVec(:), rDefLevelVecH(:)
    INTEGER                        :: ierrstat, i, iUsedDefLevels

    CHARACTER(LEN=*), PARAMETER    :: subname = MODUL_NAME//'create_cdi_zaxis'

    ! create cdi handle
    iID = zaxisCreate(iGridID, iLevels)

    ! create default vector
    IF (PRESENT(iDefLevels)) THEN
      iUsedDefLevels = iDefLevels
    ELSE
      iUsedDefLevels = iLevels
    ENDIF

!    ! considerate ocean
!    IF (PRESENT(lOcean)) THEN
!      ! check if half
!      IF (.NOT. lOcean) THEN
!        iUsedDefLevels = iUsedDefLevels - 1
!      ENDIF
!      ALLOCATE(rDefLevelVec(iUsedDefLevels), STAT=ierrstat)
!      IF (ierrstat /= SUCCESS) CALL finish(subname, ALLOCATE_FAILED)
!      ALLOCATE(rDefLevelVecH(iUsedDefLevels+1), STAT=ierrstat)
!      IF (ierrstat /= SUCCESS) CALL finish(subname, ALLOCATE_FAILED)
!      CALL set_zlev(rDefLevelVecH, rDefLevelVec)
!      IF (lOcean) THEN
!        CALL zaxisDefLevels(iID, rDefLevelVec)
!      ELSE
!        CALL zaxisDefLevels(iID, rDefLevelVecH)
!      ENDIF
!      DEALLOCATE(rDefLevelVec)
!      DEALLOCATE(rDefLevelVecH)
!    ELSE
      ALLOCATE(rDefLevelVec(iUsedDefLevels), STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(subname, ALLOCATE_FAILED)
      DO i = 1, SIZE(rDefLevelVec)
        IF (PRESENT(rDefLevelVal)) THEN
          rDefLevelVec(i) = rDefLevelVal
        ELSE
          rDefLevelVec(i) = REAL(i, wp)
        ENDIF
      END DO
      CALL zaxisDefLevels(iID, rDefLevelVec)
      DEALLOCATE(rDefLevelVec)
!    ENDIF

  END SUBROUTINE create_cdi_zaxis

  !------------------------------------------------------------------------------------------------
  !
  ! Initialize the variable list of the given restart file.
  !
  SUBROUTINE init_restart_variables(p_rf)

    TYPE(t_restart_file), POINTER, INTENT(IN) :: p_rf

    TYPE(t_var_data), POINTER     :: p_vars(:)
    TYPE(t_var_metadata), POINTER :: p_info
    CHARACTER(LEN=*), PARAMETER   :: subname = MODUL_NAME//'init_restart_variables'
    INTEGER                       :: gridID, zaxisID, varID, vlistID, iv

#ifdef DEBUG
    WRITE (nerr,FORMAT_VALS3)subname,' p_pe=',p_pe
#endif

    ! check the contained array of restart variables
    p_vars => p_rf%var_data
    IF (.NOT. ASSOCIATED(p_vars)) RETURN

    vlistID = p_rf%cdiVlistID

    ! go over the all restart variables in the associated array
    DO iv = 1, SIZE(p_vars)

      ! get pointer to metadata
      p_info => p_vars(iv)%info

      ! check time level of the field
      IF (.NOT. has_valid_time_level(p_info)) CYCLE

      ! set grid ID
      gridID = CDI_UNDEFID
      SELECT CASE (p_info%hgrid)
        CASE(GRID_UNSTRUCTURED_CELL)
          p_info%cdiGridID = p_rf%cdiCellGridID
        CASE(GRID_UNSTRUCTURED_VERT)
          p_info%cdiGridID = p_rf%cdiVertGridID
        CASE(GRID_UNSTRUCTURED_EDGE)
          p_info%cdiGridID = p_rf%cdiEdgeGridID
      END SELECT

      gridID = p_info%cdiGridID
      IF (gridID == CDI_UNDEFID) THEN
        CALL finish(subname, 'Grid type not defined for field '//TRIM(p_info%name))
      ENDIF

      ! set z axis ID
      zaxisID = p_rf%cdiZaxisIDs(p_info%vgrid)
      IF (zaxisID /= CDI_UNDEFID) THEN
        p_info%cdiZaxisID = zaxisID
      ELSE
        CALL finish(subname, 'Z axis not defined for field '//TRIM(p_info%name))
      ENDIF

      ! define the CDI variable
      varID = vlistDefVar(vlistID, gridID, zaxisID, TIME_VARIABLE)
      p_info%cdiVarID = varID
      CALL vlistDefVarDatatype(vlistID, varID, DATATYPE_FLT64)

      ! set optional parameters
      CALL vlistDefVarName(vlistID, varID, TRIM(p_info%name))
      IF (LEN_TRIM(p_info%cf%long_name) > 0) THEN
        CALL vlistDefVarLongname(vlistID, varID, p_info%cf%long_name)
      ENDIF
      IF (LEN_TRIM(p_info%cf%units) > 0) THEN
        CALL vlistDefVarUnits(vlistID, varID, p_info%cf%units)
      ENDIF

      ! currently only real valued variables are allowed, so we can always use info%missval%rval
      IF (p_info%lmiss) CALL vlistDefVarMissval(vlistID, varID, p_info%missval%rval)

#ifdef DEBUG
      IF (varID == CDI_UNDEFID) THEN
        CALL finish(subname,'CDI variable could not be defined='//TRIM(p_info%name))
      ENDIF
#endif

    ENDDO

  END SUBROUTINE init_restart_variables

  !------------------------------------------------------------------------------------------------
  !
  ! Initialize the CDI variable list of the given restart file.
  !
  SUBROUTINE init_restart_vlist(p_pd)

    TYPE(t_patch_data), POINTER, INTENT(IN) :: p_pd

    TYPE(t_restart_file), POINTER     :: p_rf
    TYPE(t_v_grid), POINTER           :: p_vgd
    REAL(wp)                          :: real_attribute
    INTEGER                           :: i, status, int_attribute, nlevp1
    LOGICAL                           :: bool_attribute
    CHARACTER(LEN=MAX_NAME_LENGTH)    :: attribute_name, text_attribute
    CHARACTER(LEN=*), PARAMETER       :: subname = MODUL_NAME//'init_restart_vlist'

    p_rf => p_pd%restart_file

#ifdef DEBUG
    WRITE (nerr,FORMAT_VALS3)subname,' p_pe=',p_pe
#endif

    p_rf%cdiVlistID = vlistCreate ()

    ! 2. add global attributes
    ! 2.1. namelists as text attributes
    DO i = 1, nmls
      status = vlistDefAttTxt(p_rf%cdiVlistID, CDI_GLOBAL,        &
           &                  TRIM(restart_namelist(i)%name),     &
           &                  LEN_TRIM(restart_namelist(i)%text), &
           &                  TRIM(restart_namelist(i)%text))
#ifdef DEBUG
      CALL check_netcdf_status(status, 'vlistDefAttTxt='// &
        &                      TRIM(attribute_name))
#endif
    ENDDO

    ! 2.2. text attributes
    DO i = 1, restart_attributes_count_text()
      CALL get_restart_attribute(i, attribute_name, text_attribute)
      status = vlistDefAttTxt(p_rf%cdiVlistID, CDI_GLOBAL,           &
           &                  TRIM(attribute_name),                  &
           &                  LEN_TRIM(text_attribute),              &
           &                  TRIM(text_attribute))
#ifdef DEBUG
      CALL check_netcdf_status(status, 'vlistDefAttTxt='// &
        &                      TRIM(attribute_name))
#endif
    END DO

    ! 2.3. real attributes
    DO i = 1, restart_attributes_count_real()
      CALL get_restart_attribute(i, attribute_name, real_attribute)
      status = vlistDefAttFlt(p_rf%cdiVlistID, CDI_GLOBAL,           &
           &                  TRIM(attribute_name),                  &
           &                  DATATYPE_FLT64,                        &
           &                  1,                                     &
           &                  real_attribute)
#ifdef DEBUG
      CALL check_netcdf_status(status, 'vlistDefAttFlt='// &
        &                      TRIM(attribute_name))
#endif
    ENDDO

    ! 2.4. integer attributes
    DO i = 1, restart_attributes_count_int()
      CALL get_restart_attribute(i, attribute_name, int_attribute)
      status = vlistDefAttInt(p_rf%cdiVlistID, CDI_GLOBAL,           &
           &                  TRIM(attribute_name),                  &
           &                  DATATYPE_INT32,                        &
           &                  1,                                     &
           &                  int_attribute)
#ifdef DEBUG
      CALL check_netcdf_status(status, 'vlistDefAttInt='// &
        &                      TRIM(attribute_name))
#endif
    ENDDO

    ! 2.5. logical attributes
    DO i = 1, restart_attributes_count_bool()
      CALL get_restart_attribute(i, attribute_name, bool_attribute)
      IF (bool_attribute) THEN
        int_attribute = 1
      ELSE
        int_attribute = 0
      ENDIF
      status = vlistDefAttInt(p_rf%cdiVlistID, CDI_GLOBAL,           &
           &                  TRIM(attribute_name),                  &
           &                  DATATYPE_INT32,                        &
           &                  1,                                     &
           &                  int_attribute)
#ifdef DEBUG
      CALL check_netcdf_status(status, 'vlistDefAttInt='// &
        &                      TRIM(attribute_name))
#endif
    ENDDO

    ! 3. add horizontal grid descriptions
    ! 3.1. cells
    CALL create_cdi_hgrid_def(p_rf%cdiCellGridID, p_pd%cells%n_glb, p_pd%cell_type,  &
      &                      'clon', 'center longitude', 'radian',                   &
      &                      'clat', 'center latitude',  'radian')

    ! 3.2. certs
    CALL create_cdi_hgrid_def(p_rf%cdiVertGridID, p_pd%verts%n_glb, 9-p_pd%cell_type,  &
      &                      'vlon', 'vertex longitude', 'radian',                     &
      &                      'vlat', 'vertex latitude',  'radian')


    ! 3.3. edges
    CALL create_cdi_hgrid_def(p_rf%cdiEdgeGridID, p_pd%edges%n_glb, 4,     &
      &                      'elon', 'edge midpoint longitude', 'radian',  &
      &                      'elat', 'edge midpoint latitude',  'radian')

    ! 4. add vertical grid descriptions
    DO i = 1, SIZE(p_pd%v_grid_defs)

      p_vgd => p_pd%v_grid_defs(i)

      SELECT CASE (p_vgd%type)

        CASE (ZA_SURFACE)
          CALL create_cdi_zaxis(p_rf%cdiZaxisIDs(ZA_SURFACE), ZAXIS_SURFACE, &
            &                   p_vgd%nlevels, 1, 0.0_wp)
        CASE (ZA_HYBRID)
          CALL create_cdi_zaxis(p_rf%cdiZaxisIDs(ZA_HYBRID), ZAXIS_HYBRID, p_vgd%nlevels)
          IF (ALLOCATED(p_pd%opt_pvct)) THEN
            nlevp1 = p_vgd%nlevels+1
            CALL zaxisDefVct(p_rf%cdiZaxisIDs(ZA_HYBRID), 2*nlevp1, &
              &              p_pd%opt_pvct(1:2*nlevp1))
          ENDIF
        CASE (ZA_HYBRID_HALF)
          CALL create_cdi_zaxis(p_rf%cdiZaxisIDs(ZA_HYBRID_HALF), ZAXIS_HYBRID_HALF, &
            &                   p_vgd%nlevels)
          IF (ALLOCATED(p_pd%opt_pvct)) THEN
            nlevp1 = p_vgd%nlevels
            CALL zaxisDefVct(p_rf%cdiZaxisIDs(ZA_HYBRID_HALF), 2*nlevp1, &
              &                               p_pd%opt_pvct(1:2*nlevp1))
          ENDIF
        CASE (ZA_DEPTH_BELOW_LAND)
          IF (p_pd%l_opt_depth_lnd) THEN
            CALL create_cdi_zaxis(p_rf%cdiZaxisIDs(ZA_DEPTH_BELOW_LAND), ZAXIS_DEPTH_BELOW_LAND, &
            &                     p_vgd%nlevels, p_pd%opt_depth_lnd)
          ENDIF
        CASE (ZA_DEPTH_BELOW_LAND_P1)
          IF (p_pd%l_opt_depth_lnd) THEN
            CALL create_cdi_zaxis(p_rf%cdiZaxisIDs(ZA_DEPTH_BELOW_LAND_P1), ZAXIS_DEPTH_BELOW_LAND, &
            &                     p_vgd%nlevels, p_pd%opt_depth_lnd+1)
          ENDIF
        CASE (ZA_DEPTH_RUNOFF_S)
          CALL create_cdi_zaxis(p_rf%cdiZaxisIDs(ZA_DEPTH_RUNOFF_S), ZAXIS_DEPTH_BELOW_LAND, &
            &                   p_vgd%nlevels)
        CASE (ZA_DEPTH_RUNOFF_G)
          CALL create_cdi_zaxis(p_rf%cdiZaxisIDs(ZA_DEPTH_RUNOFF_G), ZAXIS_DEPTH_BELOW_LAND, &
            &                   p_vgd%nlevels)
        CASE (ZA_SNOW)
          IF (p_pd%l_opt_nlev_snow) THEN
            CALL create_cdi_zaxis(p_rf%cdiZaxisIDs(ZA_SNOW), ZAXIS_GENERIC, &
            &                     p_vgd%nlevels, p_pd%opt_nlev_snow)
          ENDIF
        CASE (ZA_SNOW_HALF)
          IF (p_pd%l_opt_nlev_snow) THEN
            CALL create_cdi_zaxis(p_rf%cdiZaxisIDs(ZA_SNOW_HALF), ZAXIS_GENERIC, &
            &                     p_vgd%nlevels, p_pd%opt_nlev_snow+1)
          ENDIF
        CASE (ZA_TOA)
          CALL create_cdi_zaxis(p_rf%cdiZaxisIDs(ZA_TOA), ZAXIS_TOA, &
            &                   p_vgd%nlevels, 1, 1.0_wp)
        CASE (ZA_LAKE_BOTTOM)
          CALL create_cdi_zaxis(p_rf%cdiZaxisIDs(ZA_LAKE_BOTTOM), ZAXIS_LAKE_BOTTOM, &
            &                   p_vgd%nlevels, 1, 1.0_wp)
        CASE (ZA_MIX_LAYER)
          CALL create_cdi_zaxis(p_rf%cdiZaxisIDs(ZA_MIX_LAYER), ZAXIS_MIX_LAYER, &
            &                   p_vgd%nlevels, 1, 1.0_wp)
        CASE (ZA_LAKE_BOTTOM_HALF)
          CALL create_cdi_zaxis(p_rf%cdiZaxisIDs(ZA_LAKE_BOTTOM_HALF), ZAXIS_LAKE_BOTTOM, &
            &                   p_vgd%nlevels, 1, 1.0_wp)
        CASE (ZA_SEDIMENT_BOTTOM_TW_HALF)
          CALL create_cdi_zaxis(p_rf%cdiZaxisIDs(ZA_SEDIMENT_BOTTOM_TW_HALF), ZAXIS_SEDIMENT_BOTTOM_TW, &
            &                   p_vgd%nlevels, 1, 0.0_wp)
        CASE (ZA_HEIGHT_2M)
          CALL create_cdi_zaxis(p_rf%cdiZaxisIDs(ZA_HEIGHT_2M), ZAXIS_HEIGHT, &
            &                   p_vgd%nlevels, 1, 2.0_wp)
        CASE (ZA_HEIGHT_10M)
          CALL create_cdi_zaxis(p_rf%cdiZaxisIDs(ZA_HEIGHT_10M), ZAXIS_HEIGHT, &
            &                   p_vgd%nlevels, 1, 10.0_wp)
        CASE (ZA_DEPTH_BELOW_SEA)
          IF (p_pd%l_opt_depth) THEN
            CALL create_cdi_zaxis(p_rf%cdiZaxisIDs(ZA_DEPTH_BELOW_SEA), ZAXIS_DEPTH_BELOW_SEA, &
            &                     p_vgd%nlevels, p_pd%opt_depth, lOcean=.TRUE.)
          ENDIF
        CASE (ZA_DEPTH_BELOW_SEA_HALF)
          IF (p_pd%l_opt_depth) THEN
            CALL create_cdi_zaxis(p_rf%cdiZaxisIDs(ZA_DEPTH_BELOW_SEA_HALF), ZAXIS_DEPTH_BELOW_SEA, &
            &                     p_vgd%nlevels, p_pd%opt_depth+1, lOcean=.FALSE.)
          ENDIF
        CASE (ZA_GENERIC_ICE)
          !!!!!!!!! ATTENTION: !!!!!!!!!!!
          ! As soon as i_no_ice_thick_class is set to i_no_ice_thick_class>1 this no longer works
          CALL create_cdi_zaxis(p_rf%cdiZaxisIDs(ZA_GENERIC_ICE), ZAXIS_GENERIC, &
            &                   p_vgd%nlevels, 1, 1.0_wp)
        CASE DEFAULT
          CALL finish(subname, UNKNOWN_VERT_GRID_DESCR)
        END SELECT
    ENDDO

    ! 5. restart does contain absolute time
    p_rf%cdiTaxisID = taxisCreate(TAXIS_ABSOLUTE)
    CALL vlistDefTaxis(p_rf%cdiVlistID, p_rf%cdiTaxisID)

    ! set cdi internal time index to 0 for writing time slices in netCDF
    p_rf%cdiTimeIndex = 0

    ! init list of restart variables
    CALL init_restart_variables (p_rf)

  END SUBROUTINE init_restart_vlist

  !------------------------------------------------------------------------------------------------
  !
  ! Opens the restart file from the given parameters.
  !
  SUBROUTINE open_restart_file(p_pd)

    TYPE(t_patch_data), POINTER, INTENT(IN) :: p_pd

    TYPE(t_restart_file), POINTER :: p_rf
    TYPE(t_var_list), POINTER     :: p_re_list
    CHARACTER(LEN=32)             :: datetime
    INTEGER                       :: restart_type, i
    CHARACTER(LEN=*), PARAMETER   :: subname = MODUL_NAME//'open_restart_file'

#ifdef DEBUG
    WRITE (nerr,FORMAT_VALS3)subname,' p_pe=',p_pe
#endif

    ! just for safety
    IF(.NOT. my_process_is_restart()) CALL finish(subname, NO_RESTART_PE)

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
      CASE (FILETYPE_NC2)
        restart_type = FILETYPE_NC2
      CASE (FILETYPE_NC4)
        restart_type = FILETYPE_NC4
      CASE default
        CALL finish(subname, UNKNOWN_FILE_FORMAT)
    END SELECT

    datetime = iso8601(restart_args%datetime)

    ! build the file name
    p_rf%filename = 'restart.'// &
      &             TRIM(get_filename_noext(p_pd%base_filename))//'_'// &
      &             TRIM(datetime)//'_'// &
      &             TRIM(p_rf%model_type)//'.nc'

    p_rf%cdiFileID = streamOpenWrite(p_rf%filename, restart_type)

    IF (p_rf%cdiFileID < 0) THEN
      WRITE(message_text,'(a)') cdiStringError(p_rf%cdiFileID)
      CALL message('', message_text, all_print=.TRUE.)
      CALL finish (subname, 'open failed on '//TRIM(p_rf%filename))
    ELSE
      CALL message (subname, 'opened '//TRIM(p_rf%filename), all_print=.TRUE.)
    END IF

#ifdef DEBUG
    WRITE (nerr, FORMAT_VALS5)subname,' p_pe=',p_pe,' open netCDF file with ID=',p_rf%cdiFileID
#endif

    CALL init_restart_vlist(p_pd)

    CALL streamDefVlist(p_rf%cdiFileID, p_rf%cdiVlistID)

  END SUBROUTINE open_restart_file

  !------------------------------------------------------------------------------------------------
  !
  ! Closes the given restart file.
  !
  SUBROUTINE close_restart_file(rf)

    TYPE (t_restart_file), INTENT(INOUT)  :: rf
    CHARACTER(LEN=*), PARAMETER           :: subname = MODUL_NAME//'close_restart_file'

#ifdef DEBUG
    WRITE (nerr,FORMAT_VALS3)subname,' p_pe=',p_pe
#endif

    ! just for safety
    IF(.NOT. my_process_is_restart()) CALL finish(subname, NO_RESTART_PE)

    IF (rf%cdiFileID /= CDI_UNDEFID) THEN

#ifdef DEBUG
      WRITE (nerr,'(3a)')subname,' try to close restart file=',TRIM(rf%filename)
#endif
      CALL streamClose(rf%cdiFileID)

#ifdef DEBUG
      WRITE (nerr, FORMAT_VALS5)subname,' p_pe=',p_pe,' close netCDF file with ID=',rf%cdiFileID
#endif

      rf%cdiFileID  = CDI_UNDEFID
      rf%filename   = ''
      rf%linkname   = ''
      rf%linkprefix = ''

    ENDIF

  END SUBROUTINE close_restart_file

  !------------------------------------------------------------------------------------------------
  !
  ! Creates a symbolic link from the given restart file.
  !
  SUBROUTINE create_restart_file_link (rf, proc_id)

    TYPE (t_restart_file), INTENT(INOUT)  :: rf
    INTEGER, INTENT(IN)                   :: proc_id

    INTEGER                               :: iret, id
    CHARACTER(LEN=5)                      :: str_id
    CHARACTER(LEN=*), PARAMETER           :: subname = MODUL_NAME//'create_restart_file_link'

    ! build link name
    id = proc_id - p_restart_pe0
    IF (id == 0) THEN
      str_id = ' '
    ELSE
      WRITE(str_id, '(I5)')id
    ENDIF
    rf%linkprefix = 'restart'//TRIM(str_id)
    rf%linkname = TRIM(rf%linkprefix)//'_'//TRIM(rf%model_type)//'.nc'

    ! delete old symbolic link, if exists
    IF (util_islink(TRIM(rf%linkname))) THEN
      iret = util_unlink(TRIM(rf%linkname))
      IF (iret /= SUCCESS) THEN
          WRITE (nerr,'(3a)')subname,' cannot unlink ',TRIM(rf%linkname)
      ENDIF
    ENDIF

    ! create a new symbolic link
    iret = util_symlink(TRIM(rf%filename),TRIM(rf%linkname))
    IF (iret /= SUCCESS) THEN
      WRITE (nerr,'(5a)')subname,' cannot create symbolic link ', &
        & TRIM(rf%linkname),' for ', TRIM(rf%filename)
    ENDIF

  END SUBROUTINE create_restart_file_link

  !------------------------------------------------------------------------------------------------
  !
  ! Creates an information file from the given restart file.
  !
  SUBROUTINE create_restart_info_file(rf)

    TYPE (t_restart_file), INTENT(IN)  :: rf
    INTEGER :: nrf
    CHARACTER(LEN=20)                  :: info_file_name
    CHARACTER(LEN=*), PARAMETER        :: subname = MODUL_NAME//'create_restart_info_file'

    nrf = find_next_free_unit(10,100)
    info_file_name = TRIM(rf%linkprefix)//'.info'
    OPEN(nrf, file=info_file_name)
    WRITE(nrf, '(4a)')'gridspec: ',TRIM(rf%linkname),' ! ',TRIM(rf%filename)
    CLOSE(nrf)
    CALL message('',TRIM(info_file_name)//' written')

  END SUBROUTINE create_restart_info_file

#endif

END MODULE mo_io_restart_async

