!----------------------------------
#include "dsl_definitions.inc"
#include "omp_definitions.inc"
!>
!! A collection of MPI communication tools
!!
!! @par Revision History
!! First version by Leonidas Linardakis,  MPI-M, November 2011.
!!
!! @par
!!
!! @par Copyright
!! 2002-2010 by DWD and MPI-M
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
!! $Id: n/a$
!!
MODULE mo_icon_comm_lib

  USE mo_kind,            ONLY: wp
  USE mo_io_units,        ONLY: filename_max
  USE mo_exception,       ONLY: message_text, message, finish, warning
  USE mo_parallel_config, ONLY: nproma, icon_comm_debug, max_send_recv_buffer_size, &
    & icon_comm_method, icon_comm_openmp

  USE mo_communication,   ONLY: blk_no, idx_no
  USE mo_model_domain,    ONLY: t_patch
  USE mo_mpi,             ONLY: p_send, p_recv, p_irecv, p_wait, p_isend, &
     & p_real_dp, p_int, p_bool, my_process_is_mpi_seq,   &
     & process_mpi_all_comm, work_mpi_barrier, p_stop, &
     & get_my_mpi_work_communicator, get_my_mpi_work_comm_size, &
     & get_my_mpi_work_id
  USE mo_timer,           ONLY: ltimer, timer_start, timer_stop, timer_icon_comm_sync, &
    & activate_sync_timers, timer_icon_comm_fillrecv, timer_icon_comm_wait, &
    & timer_icon_comm_ircv, timer_icon_comm_fillsend, timer_icon_comm_fillandsend, &
    & timer_icon_comm_isend

  USE mo_master_control,  ONLY: get_my_process_name
#ifndef NOMPI
  USE mpi
#endif
  
  USE mo_impl_constants,     ONLY: HALO_LEVELS_CEILING
  
!   USE mo_impl_constants,  ONLY: &
!     & min_rlcell, max_rlcell, min_rlcell_int, &
!     & min_rlvert, max_rlvert,                 & ! min_rlvert_int,
!     & min_rledge, max_rledge, min_rledge_int

#ifdef _OPENMP
  USE omp_lib, ONLY: omp_in_parallel
#endif

  IMPLICIT NONE

  PRIVATE

  ! public constants
  PUBLIC :: cells_not_owned, cells_not_in_domain, cells_one_edge_in_domain
  PUBLIC :: edges_not_owned, edges_not_in_domain
  PUBLIC :: verts_not_owned, verts_not_in_domain
  
  PUBLIC :: is_ready, until_sync
   
  ! public methods
  PUBLIC :: init_icon_comm_lib       ! the first call to the icon_comm_lib 
  PUBLIC :: init_icon_comm_patterns  ! initilize comm_patterns for a patch
                                     ! must be called before any communication variable is
                                     ! created for this patch
  PUBLIC :: destruct_icon_comm_lib
  
  PUBLIC :: new_icon_comm_variable
  PUBLIC :: delete_icon_comm_variable
  
  PUBLIC :: icon_comm_var_is_ready  ! The comm_variable can be communicated
  PUBLIC :: icon_comm_sync
  PUBLIC :: icon_comm_sync_all

!   PUBLIC :: icon_comm_show

  CHARACTER(LEN=*), PARAMETER :: version = '$Id$'
  !--------------------------------------------------------------
  
  !--------------------------------------------------------------
  ! internal parameters
  INTEGER, PARAMETER ::  max_no_of_comm_variables = 20
  INTEGER, PARAMETER ::  max_no_of_comm_processes = 32
  INTEGER, PARAMETER ::  max_no_of_patches = 1
!  INTEGER, PARAMETER ::  max_send_recv_buffer_size = 131072
  REAL(wp), PARAMETER::  header_separator = 99999999.0_wp
  INTEGER, PARAMETER ::  checksum_size = 2
  ! TAG parameters
  INTEGER, PARAMETER ::  halo_tag = 1
  INTEGER, PARAMETER ::  internal_tag = 2
  !--------------------------------------------------------------
  ! status / requests
  INTEGER, PARAMETER ::  not_active = 0
  INTEGER, PARAMETER ::  active = 1
  INTEGER, PARAMETER ::  is_ready = 10    ! essentially the same as communicate
  INTEGER, PARAMETER ::  communicate = 11

  !--------------------------------------------------------------
  ! grid locations.
  ! Note that these are used also as indexes for the communication patterns

  ! Halo comm patterns: these will eventually be specialized
  INTEGER ::  cells_not_in_domain = -1
  INTEGER ::  cells_not_owned = -1          ! = cells_not_in_domain
  INTEGER ::  cells_one_edge_in_domain = -1
  INTEGER ::  edges_not_owned = -1
  INTEGER ::  edges_not_in_domain = -1
  INTEGER ::  verts_not_owned = -1
  INTEGER ::  verts_not_in_domain = -1

  ! non halo comm patterns, these should eventually go to a different universe
  INTEGER :: radiation_repartition

  !--------------------------------------------------------------
  INTEGER ::  max_comm_patterns = 0
  INTEGER, PARAMETER ::  allocated_comm_patterns = 12

  !--------------------------------------------------------------
  ! grid dimensions
  INTEGER, PARAMETER ::  grid_2D = 2
  INTEGER, PARAMETER ::  grid_3D = 3
  !--------------------------------------------------------------
  ! scope parameters
  INTEGER, PARAMETER ::  global = 0
  INTEGER, PARAMETER ::  until_sync = 1

  !--------------------------------------------------------------
  ! TYPE definitions

  !--------------------------------------------------------------
  !> Holds the buffer pointers to procs to be sent and recieved
  TYPE t_comm_process_buffer
  
     INTEGER :: pid             ! the process to communicate to
     
     INTEGER :: start_index      ! the start point in the whole bufffer
     INTEGER :: end_index        ! the end   point in the whole bufffer
     INTEGER :: current_index    ! the current point in the whole bufffer
     INTEGER :: buffer_size      ! the size of this comm to process bufffer
       
  END TYPE t_comm_process_buffer
  !--------------------------------------------------------------

  !--------------------------------------------------------------
  !> Holds the communcation pattern for each of the process to which we communicate
  TYPE t_process_comm_pattern

    INTEGER :: pid           ! the process to which we communicate
    INTEGER :: no_of_points  ! the (horizontal) number of values to be communicated
    INTEGER, ALLOCATABLE :: global_index(:) ! the global index of the variable
    INTEGER, ALLOCATABLE :: block_no(:)  ! the local block of the variable
    INTEGER, ALLOCATABLE :: index_no(:)  ! the local index of the variable
  
    INTEGER :: buffer_index ! the index of the send/recv buffer handler for this process
       
  END TYPE t_process_comm_pattern
  !--------------------------------------------------------------
  
  !--------------------------------------------------------------
  !> Holds the communication patterns to all procs to be communicated.
  !  One structure represents one overall pattern for cells, edges, etc.
  TYPE t_grid_comm_pattern
  
    INTEGER :: id   ! the id of this, it is the index in the array it lives
  
    ! the send communication pattern for each process to sent to
    INTEGER :: no_of_send_procs
    TYPE(t_process_comm_pattern), POINTER :: send(:)
    ! the recv communication pattern for each process to receive from
    INTEGER :: no_of_recv_procs
    TYPE(t_process_comm_pattern), POINTER :: recv(:)
       
    CHARACTER(len=32) :: name
    
    INTEGER :: status 
       
  END TYPE t_grid_comm_pattern
  !--------------------------------------------------------------
  
  !--------------------------------------------------------------
  !> Holds the variables to be communicated
  TYPE t_comm_variable_real
    INTEGER :: status                          ! not_active, active, ...
    INTEGER :: request                         ! not_active, ...
    INTEGER :: comm_status                     ! not_active, ...
    INTEGER :: comm_buffer                     ! for future use
    INTEGER :: scope                           ! until_sync = temporary, will be deleted after sync
                                               ! global = will be kept until deleted
     
    INTEGER :: comm_pattern_index                   ! location where variable lives:
                                                 ! edges, vertices, cells
    INTEGER :: grid_dim                        ! 1D, 2D 3D
    TYPE(t_patch), POINTER :: p_patch          ! the patch the varibale exists
    TYPE(t_grid_comm_pattern), POINTER :: grid_comm_pattern    ! the communication pattern

    INTEGER :: no_of_variables                 ! how many variables stored in the array
                                                 ! if =0 then only one
    
!     INTEGER :: dim_1                         ! should be nproma
    INTEGER :: vertical_layers                 ! if 3D: nlevels
!     INTEGER :: dim_3                         !
    INTEGER :: dim_4                           ! =no_of_variables
        
    REAL(wp), POINTER :: values_2d(:,:)     ! nproma, nblocks
    REAL(wp), GENERAL_3D, POINTER :: values_3d   ! if 3D: nproma, vertical layers, nblocks
    REAL(wp), POINTER :: values_4d(:,:,:,:) ! nproma, vertical layers, nblocks

    CHARACTER(len=32) :: name
    
  END TYPE t_comm_variable_real

  !--------------------------------------------------------------

  
  !> The array of the grid objects.
  TYPE(t_comm_variable_real), TARGET :: comm_variable(max_no_of_comm_variables)
  
  TYPE(t_comm_process_buffer), TARGET :: send_procs_buffer(max_no_of_comm_processes)
  
  TYPE(t_comm_process_buffer), TARGET :: recv_procs_buffer(max_no_of_comm_processes)
  
  TYPE(t_grid_comm_pattern), TARGET ::  &
    & grid_comm_pattern_list(allocated_comm_patterns, max_no_of_patches)

  ! At the moment we only have on sett of communication buffers
  ! This means that only one communation bulk can be executed each time
  REAL(wp), ALLOCATABLE :: send_buffer(:)
  REAL(wp), ALLOCATABLE :: recv_buffer(:)
  INTEGER  :: buffer_comm_status
 
  !> The number of actual active comm_variables
  INTEGER :: active_comm_variables
  !> The max id of the active comm_variables
  INTEGER :: max_active_comm_variables
  
  !> The number of actual active send_buffers
  INTEGER :: active_send_buffers
  
  !> The number of actual active recv_buffers
  INTEGER :: active_recv_buffers
   
   
  LOGICAL :: comm_lib_is_initialized

  INTEGER :: log_file_id

  !-------------------------------------------------------------------------
  INTEGER :: my_work_communicator
  INTEGER :: my_work_comm_size
  INTEGER :: my_mpi_work_id
  LOGICAL :: this_is_mpi_sequential
  
  !-------------------------------------------------------------------------
  INTERFACE new_icon_comm_variable
    MODULE PROCEDURE new_comm_variable_r2d
    MODULE PROCEDURE new_comm_variable_r3d
!     MODULE PROCEDURE new_comm_variable_r3d_target
    MODULE PROCEDURE new_comm_variable_r4d
  END INTERFACE
  !-------------------------------------------------------------------------
  INTERFACE icon_comm_sync
    MODULE PROCEDURE icon_comm_sync_2D_1
    MODULE PROCEDURE icon_comm_sync_3D_1
    MODULE PROCEDURE icon_comm_sync_3D_2
    MODULE PROCEDURE icon_comm_sync_3D_3
  END INTERFACE
  
  !-------------------------------------------------------------------------
CONTAINS

  !-----------------------------------------------------------------------
  !> NOTE: Not completed yet!
  SUBROUTINE destruct_icon_comm_lib()
    CHARACTER(*), PARAMETER :: method_name = "destruct_icon_comm_lib"

    IF (this_is_mpi_sequential) RETURN
#ifdef _OPENMP
    IF (omp_in_parallel()) &
      CALL finish(method_name, 'cannot be called from openmp parallel')
#endif
    
    DEALLOCATE(send_buffer, recv_buffer)
    CLOSE(log_file_id)
    
  END SUBROUTINE destruct_icon_comm_lib
  !-----------------------------------------------------------------------
  
  !-----------------------------------------------------------------------
  !>
  SUBROUTINE init_icon_comm_lib()

    INTEGER :: i,k, return_status
    LOGICAL :: unit_is_occupied
    
    CHARACTER(*), PARAMETER :: method_name = "init_icon_comm_lib"

    this_is_mpi_sequential = my_process_is_mpi_seq()
!     write(0,*) method_name, " this_is_mpi_sequential=", this_is_mpi_sequential
    
    IF(this_is_mpi_sequential) RETURN
    
#ifdef _OPENMP
    IF (omp_in_parallel()) &
      CALL finish(method_name, 'cannot be called from openmp parallel')
#endif

    IF ( comm_lib_is_initialized ) THEN
      CALL message(method_name, 'cannot be called more than once')
      RETURN
    END IF

    ALLOCATE(send_buffer(max_send_recv_buffer_size), &
      & recv_buffer(max_send_recv_buffer_size),stat=return_status)
    IF (return_status > 0) &
      CALL finish (method_name, 'ALLOCATE send,recv buffers failed')

    active_comm_variables = 0
    max_active_comm_variables = 0
    active_send_buffers = 0
    active_recv_buffers = 0
    comm_lib_is_initialized = .TRUE.
    buffer_comm_status = not_active
    
    DO i=1,max_no_of_comm_variables
      CALL clear_comm_variable(i)
    ENDDO

    DO i=1,max_no_of_patches
      DO k=1,max_comm_patterns
        grid_comm_pattern_list(k,i)%status = not_active
      ENDDO
    ENDDO

    ! initialize the communicators
    my_work_communicator = get_my_mpi_work_communicator()
    my_work_comm_size    = get_my_mpi_work_comm_size()
    my_mpi_work_id       = get_my_mpi_work_id()
    
!    IF (icon_comm_debug) THEN
    DO log_file_id = 500, 5000
      INQUIRE (UNIT=log_file_id, OPENED=unit_is_occupied)
      IF ( .NOT. unit_is_occupied ) EXIT
    ENDDO
    IF (unit_is_occupied) &
      CALL finish(method_name, "Cannot find avaliable file unit")
    WRITE(message_text,'(a,a,a,i4.4)') 'log.', TRIM(get_my_process_name()), &
      & ".icon_comm.", my_mpi_work_id
    OPEN (log_file_id, FILE=TRIM(message_text))
!    ENDIF
    
    RETURN

  END SUBROUTINE init_icon_comm_lib
  !-----------------------------------------------------------------------
  
  !-----------------------------------------------------------------------
  !>
  SUBROUTINE init_icon_comm_patterns(p_patch)
    TYPE(t_patch), INTENT(IN) :: p_patch

    INTEGER :: i
    CHARACTER(*), PARAMETER :: method_name = "init_icon_comm_patterns"
    
    max_comm_patterns = 0
    IF(this_is_mpi_sequential) RETURN
    ! set id of grid_comm_pattern_list to identity
    DO i = 1, max_comm_patterns
    ENDDO

#ifdef _OPENMP
    IF (omp_in_parallel()) &
      CALL finish(method_name, 'cannot be called from openmp parallel')
#endif
    
    IF ( p_patch%id > max_no_of_patches) &
      CALL finish(method_name, "p_patch%id > max_no_of_patches")

    ! halo cells comm_pattern
!     CALL work_mpi_barrier()
!     write(0,*) my_mpi_work_id, method_name, "setup_grid_comm_pattern cells_not_in_domain..."
    cells_not_in_domain = new_halo_comm_pattern(p_patch%id, &
      & p_patch%n_patch_cells,   p_patch%cells%owner_local, &
      & p_patch%cells%glb_index, p_patch%cells%loc_index,   &
      & name="cells_not_in_domain" )
    cells_not_owned = cells_not_in_domain
               
    cells_one_edge_in_domain = new_halo_comm_pattern(p_patch%id,&
      & p_patch%n_patch_cells,   p_patch%cells%owner_local, &
      & p_patch%cells%glb_index, p_patch%cells%loc_index,   &
      & halo_level=p_patch%cells%halo_level, level_start=1, level_end=1,&
      & name="cells_not_in_domain" )
            
    ! halo edges comm_pattern
!     CALL work_mpi_barrier()
!     write(0,*) my_mpi_work_id, method_name, "setup_grid_comm_pattern edges_not_owned..."
    edges_not_owned = new_halo_comm_pattern(p_patch%id, &
      & p_patch%n_patch_edges,   p_patch%edges%owner_local, &
      & p_patch%edges%glb_index, p_patch%edges%loc_index,   &
      & name="edges_not_owned")
    
    edges_not_in_domain = new_halo_comm_pattern(p_patch%id,&
      & p_patch%n_patch_edges,   p_patch%edges%owner_local, &
      & p_patch%edges%glb_index, p_patch%edges%loc_index,   &
      & halo_level=p_patch%edges%halo_level, level_start=2, level_end=HALO_LEVELS_CEILING,&
      & name="edges_not_in_domain")
    
    ! halo verts comm_pattern
!     CALL work_mpi_barrier()
!     write(0,*) my_mpi_work_id, method_name, "setup_grid_comm_pattern verts_not_owned..."
    verts_not_owned = new_halo_comm_pattern(p_patch%id, &
      & p_patch%n_patch_verts,   p_patch%verts%owner_local, &
      & p_patch%verts%glb_index, p_patch%verts%loc_index,   &
      & name="verts_not_owned" )
        
    verts_not_in_domain = new_halo_comm_pattern(p_patch%id, &
      & p_patch%n_patch_verts,   p_patch%verts%owner_local, &
      & p_patch%verts%glb_index, p_patch%verts%loc_index,   &
      & halo_level=p_patch%verts%halo_level, level_start=2, level_end=HALO_LEVELS_CEILING,&
      & name="verts_not_in_domain" )
        
    CALL print_grid_comm_stats(grid_comm_pattern_list(cells_not_in_domain, p_patch%id))
    CALL print_grid_comm_stats(grid_comm_pattern_list(edges_not_owned, p_patch%id))
    CALL print_grid_comm_stats(grid_comm_pattern_list(verts_not_owned, p_patch%id))
    IF ( icon_comm_debug) THEN
      CALL print_grid_comm_pattern(grid_comm_pattern_list(cells_not_in_domain, p_patch%id))
      CALL print_grid_comm_pattern(grid_comm_pattern_list(edges_not_owned, p_patch%id))
      CALL print_grid_comm_pattern(grid_comm_pattern_list(verts_not_owned, p_patch%id))
    ENDIF    
        
 !   CALL finish("init_icon_comm_patterns","ends")

  END SUBROUTINE init_icon_comm_patterns
  !-----------------------------------------------------------------------
  
  !-----------------------------------------------------------------------
  !>
  INTEGER FUNCTION new_halo_comm_pattern(p_patch_id, total_no_of_points, &
    & owner, global_index, local_index, halo_level, level_start, level_end, name)

    INTEGER, INTENT(IN) :: p_patch_id
    INTEGER, INTENT(in) :: total_no_of_points
    INTEGER, INTENT(in) :: owner(:), global_index(:), local_index(:)
    INTEGER, INTENT(in), OPTIONAL :: halo_level(:,:), level_start, level_end
    CHARACTER(*), INTENT(in) :: name

    max_comm_patterns = max_comm_patterns + 1
    IF (max_comm_patterns > allocated_comm_patterns) THEN
      CALL finish("new_halo_comm_pattern","max_comm_patterns > allocated_comm_patterns")
    ENDIF
       
    grid_comm_pattern_list(max_comm_patterns,p_patch_id)%id = max_comm_patterns
    
    CALL setup_grid_comm_pattern(grid_comm_pattern_list(max_comm_patterns, p_patch_id), &
      & total_no_of_points, owner, global_index, local_index, &
      & halo_level, level_start, level_end, name)

    new_halo_comm_pattern = max_comm_patterns

  END FUNCTION new_halo_comm_pattern
  !-----------------------------------------------------------------------
  
  !-----------------------------------------------------------------------
  !>
  SUBROUTINE setup_grid_comm_pattern(grid_comm_pattern, total_no_of_points, &
    & owner, global_index, local_index, halo_level, level_start, level_end, name)
!    & owner, global_index, local_index, name)

    
    TYPE(t_grid_comm_pattern), INTENT(inout) :: grid_comm_pattern
    INTEGER, INTENT(in) :: total_no_of_points
    INTEGER, INTENT(in) :: owner(:), global_index(:), local_index(:)
    INTEGER, INTENT(in), OPTIONAL :: halo_level(:,:), level_start, level_end
    CHARACTER(*), INTENT(in) :: name

    TYPE(t_process_comm_pattern), POINTER :: p_comm_pattern
    
    INTEGER :: comm_points(max_no_of_comm_processes)
    INTEGER :: procs_id(max_no_of_comm_processes)
    INTEGER :: buffer_id(max_no_of_comm_processes)
    INTEGER :: comm_of_buffer_id(max_no_of_comm_processes)

    INTEGER, ALLOCATABLE :: recv_requests(:), total_requests(:)
    INTEGER, ALLOCATABLE :: recv_global_indexes(:,:), send_global_indexes(:,:)
    
    INTEGER :: i, point_idx, bfid, no_comm_procs, no_of_points
    INTEGER :: owner_id, max_comm_points, max_buffer_size, no_of_recv_requests
    INTEGER :: local_idx, return_status
    
    CHARACTER(*), PARAMETER :: method_name = "setup_grid_comm_pattern"

#ifndef NOMPI

!     IF (my_mpi_work_id==1) CALL work_mpi_barrier()
!     write(0,*) my_mpi_work_id, "global_index:", global_index
!     write(0,*) my_mpi_work_id, "owner:", owner
!     IF (my_mpi_work_id==0) CALL work_mpi_barrier()
    
    ! first count how many recv procs we have
    ! and how many recv points per procs
    comm_of_buffer_id(:) = -1    
    no_comm_procs=0
    comm_points(:) = 0
    DO point_idx = 1, total_no_of_points
    
      owner_id = owner(point_idx)
      IF(owner_id < 0) CYCLE
      IF(owner_id == my_mpi_work_id) CYCLE
      IF(PRESENT(halo_level)) THEN
        IF (halo_level(idx_no(point_idx), blk_no(point_idx)) < level_start .OR. &
          & halo_level(idx_no(point_idx), blk_no(point_idx)) > level_end)  CYCLE
      ENDIF
      
      bfid = get_recvbuffer_id_of_pid(owner_id)
      IF (comm_points(bfid) == 0) THEN
        no_comm_procs = no_comm_procs + 1
        procs_id(no_comm_procs)  = owner_id
        buffer_id(no_comm_procs) = bfid
        comm_of_buffer_id(bfid)  = no_comm_procs
      ENDIF
      comm_points(bfid) = comm_points(bfid) + 1
      
    ENDDO

    ! allocate the recv patterns
    grid_comm_pattern%no_of_recv_procs = no_comm_procs
    ALLOCATE(grid_comm_pattern%recv(no_comm_procs),stat=return_status)
    IF (return_status > 0) &
      CALL finish (method_name, 'ALLOCATE(grid_comm_pattern%recv(no_comm_procs)')
    
    DO i=1,no_comm_procs
      p_comm_pattern => grid_comm_pattern%recv(i)
      p_comm_pattern%pid          = procs_id(i)
      p_comm_pattern%buffer_index = buffer_id(i)
      no_of_points                = comm_points(buffer_id(i))
      p_comm_pattern%no_of_points = 0 ! this will be used as a counter for filling
         ! the pattern indexes. At the end will be checked against
         ! comm_points(procs_comm_pattern%buffer_index) for consistency
      ALLOCATE(p_comm_pattern%global_index(no_of_points), &
        & p_comm_pattern%block_no(no_of_points), &
        & p_comm_pattern%index_no(no_of_points), &
        & stat=return_status)
      IF (return_status > 0) &
        CALL finish (method_name, 'ALLOCATE recv patterns')
                       
    ENDDO
  
    ! fill the indexes of the recv patterns
    DO point_idx = 1, total_no_of_points
    
      IF(owner(point_idx) < 0) CYCLE
      IF(owner(point_idx) == my_mpi_work_id) CYCLE
      IF(PRESENT(halo_level)) THEN
        IF (halo_level(idx_no(point_idx), blk_no(point_idx)) < level_start .OR. &
          & halo_level(idx_no(point_idx), blk_no(point_idx)) > level_end)  CYCLE
      ENDIF
    
      bfid = get_recvbuffer_id_of_pid(owner(point_idx))
      p_comm_pattern => grid_comm_pattern%recv(comm_of_buffer_id(bfid))
      p_comm_pattern%no_of_points = p_comm_pattern%no_of_points + 1
      p_comm_pattern%global_index(p_comm_pattern%no_of_points) = &
        & global_index(point_idx)
      p_comm_pattern%block_no(p_comm_pattern%no_of_points) = &
        & blk_no(point_idx)
      p_comm_pattern%index_no(p_comm_pattern%no_of_points) = &
        & idx_no(point_idx)
      
    ENDDO

    ! Check if the p_comm_pattern%recv(i)%no_of_points
    ! are consistent with comm_points(bfid)
    ! Compute the max_comm_points
    max_comm_points = 0
    DO i=1,no_comm_procs
      IF (grid_comm_pattern%recv(i)%no_of_points /= &
        & comm_points(buffer_id(i))) THEN
        CALL finish (method_name, 'incosistent recv comm_points')
      ENDIF
      max_comm_points = MAX(max_comm_points, grid_comm_pattern%recv(i)%no_of_points)
    ENDDO

    ! The receive patterns have been constructed
    ! Next construct the send patterns from information
    ! received from the requesting processes.

    ! We use MPI_ALLREDUCE sum to calculate the total number
    ! of requests for send
    ! Not the best approach, but not the worst neither
    ALLOCATE(recv_requests(0:my_work_comm_size-1),&
      & total_requests(0:my_work_comm_size-1), stat=return_status)
    IF (return_status > 0) &
      CALL finish (method_name, 'ALLOCATE(recv_requests)')
    recv_requests(:)=0
    DO i=1,no_comm_procs
      recv_requests(grid_comm_pattern%recv(i)%pid) = 1
    ENDDO

    CALL MPI_ALLREDUCE(recv_requests, total_requests, my_work_comm_size, &
      & MPI_INTEGER, MPI_SUM, my_work_communicator, return_status)
    IF (return_status /= MPI_SUCCESS) THEN
      CALL finish (method_name, 'MPI_ALLREDUCE failed')
    ENDIF

    ! we have total_requests(my_mpi_work_id)
    ! for this setup this should be equal to no_comm_procs
    no_of_recv_requests = total_requests(my_mpi_work_id)
    IF (no_of_recv_requests /= no_comm_procs) THEN
      write(0,*) name, " no_of_recv_requests=", no_of_recv_requests,&
        & "no_comm_procs", no_comm_procs
      CALL warning (method_name, 'total_requests /= no_comm_procs')
    ENDIF
    DEALLOCATE(recv_requests, total_requests)

    ! Allocate the buffers for sending and receiving the global indexes
    max_buffer_size =  (max_comm_points * 2) + 16 
    ALLOCATE(recv_global_indexes(max_buffer_size, no_of_recv_requests),&
      & send_global_indexes(max_buffer_size, no_comm_procs), stat=return_status)
    IF (return_status > 0) &
      CALL finish (method_name, 'ALLOCATE(recv_global_indexes)')

    ! start receiving global indexes
    DO i=1, no_of_recv_requests
      CALL p_irecv(recv_global_indexes(:,i), MPI_ANY_SOURCE, &
          & p_tag=internal_tag, p_count=max_buffer_size, comm=my_work_communicator)
    ENDDO
      
    ! Fill and sent global indexes
    DO i=1, no_comm_procs
      p_comm_pattern => grid_comm_pattern%recv(i)
      no_of_points = p_comm_pattern%no_of_points
      send_global_indexes(1,i) = my_mpi_work_id
      send_global_indexes(2,i) = no_of_points
      send_global_indexes(3:no_of_points+2,i) = p_comm_pattern%global_index(1:no_of_points)
      CALL p_isend(send_global_indexes(:,i), p_comm_pattern%pid, &
          & p_tag=internal_tag, p_count=no_of_points+2, comm=my_work_communicator)
    ENDDO
    
    ! Wait for all requests to finish
    CALL p_wait

    ! Make some room
    DEALLOCATE(send_global_indexes)
    
    ! Allocate the send patterns
    grid_comm_pattern%no_of_send_procs = no_of_recv_requests
    ALLOCATE(grid_comm_pattern%send(no_of_recv_requests),stat=return_status)
    IF (return_status > 0) &
      CALL finish (method_name, 'ALLOCATE(grid_comm_pattern%end(no_of_recv_requests))')
    
    ! Fill the send patterns from the recv_global_indexes
    DO i=1,no_of_recv_requests
      p_comm_pattern => grid_comm_pattern%send(i)
      p_comm_pattern%pid          = recv_global_indexes(1,i)
      p_comm_pattern%no_of_points = recv_global_indexes(2,i)      
      p_comm_pattern%buffer_index = get_sendbuffer_id_of_pid(p_comm_pattern%pid)
      ALLOCATE(p_comm_pattern%global_index(p_comm_pattern%no_of_points), &
        & p_comm_pattern%block_no(p_comm_pattern%no_of_points), &
        & p_comm_pattern%index_no(p_comm_pattern%no_of_points), &
        & stat=return_status)
      IF (return_status > 0) &
        CALL finish (method_name, 'ALLOCATE send patterns')
      ! fill the global_index
      p_comm_pattern%global_index(1:p_comm_pattern%no_of_points) = &
        & recv_global_indexes(3:p_comm_pattern%no_of_points+2,i)                       
    ENDDO
  
    ! calculate the local indexes
    DO i=1,no_of_recv_requests
      p_comm_pattern => grid_comm_pattern%send(i)
      DO point_idx = 1, p_comm_pattern%no_of_points
        local_idx = local_index(p_comm_pattern%global_index(point_idx))
        IF ( local_idx < 0 ) &
          & CALL finish(method_name,'Wrong local index')
        p_comm_pattern%block_no(point_idx) = blk_no(local_idx)
        p_comm_pattern%index_no(point_idx) = idx_no(local_idx)
      ENDDO
    ENDDO

   
    DEALLOCATE(recv_global_indexes)
    
    grid_comm_pattern%status = active
    grid_comm_pattern%name   = TRIM(name)
#endif
    
    CALL work_mpi_barrier()
 
  END SUBROUTINE setup_grid_comm_pattern
  !-----------------------------------------------------------------------
  
  !-----------------------------------------------------------------------
  !>
  SUBROUTINE print_grid_comm_stats(grid_comm_pattern)
    TYPE(t_grid_comm_pattern), INTENT(in) :: grid_comm_pattern

    INTEGER :: i, min_points, max_points, tot_points
           
    write(log_file_id,*) " === Communication stats for ", TRIM(grid_comm_pattern%name), &
      & " id=", grid_comm_pattern%id
    write(log_file_id,*) " --- recv: no_of_procs= ", grid_comm_pattern%no_of_recv_procs
    min_points = 99999999
    max_points = 0
    tot_points = 0
    DO i=1,grid_comm_pattern%no_of_recv_procs
      min_points = MIN(min_points, grid_comm_pattern%recv(i)%no_of_points)
      max_points = MAX(max_points, grid_comm_pattern%recv(i)%no_of_points)
      tot_points = tot_points + grid_comm_pattern%recv(i)%no_of_points
    ENDDO
    write(log_file_id,*) " - min halos=",min_points, " max halos=",max_points, &
      & " total halos=", tot_points
       
    write(log_file_id,*) " --- send: no_of_procs= ", grid_comm_pattern%no_of_send_procs
    min_points = 99999999
    max_points = 0
    tot_points = 0
    DO i=1,grid_comm_pattern%no_of_send_procs
      min_points = MIN(min_points, grid_comm_pattern%send(i)%no_of_points)
      max_points = MAX(max_points, grid_comm_pattern%send(i)%no_of_points)
      tot_points = tot_points + grid_comm_pattern%send(i)%no_of_points
    ENDDO
    write(log_file_id,*) " - min halos=",min_points, " max halos=",max_points, &
      & " total halos=", tot_points
      
        
    write(log_file_id,*) " === End Communication stats for ", TRIM(grid_comm_pattern%name)
    
  END SUBROUTINE print_grid_comm_stats
  !-----------------------------------------------------------------------

    
  !-----------------------------------------------------------------------
  !>
  SUBROUTINE print_grid_comm_pattern(grid_comm_pattern)
    TYPE(t_grid_comm_pattern), INTENT(in) :: grid_comm_pattern

    INTEGER :: i
        
    IF ( .NOT. icon_comm_debug) RETURN
    
    write(log_file_id,*) " === Communication pattern info for ", TRIM(grid_comm_pattern%name), &
      & " id=", grid_comm_pattern%id
    write(log_file_id,*) " - no_of_recv_procs= ", grid_comm_pattern%no_of_recv_procs
    write(log_file_id,*) " - no_of_send_procs= ", grid_comm_pattern%no_of_send_procs
    
    DO i=1,grid_comm_pattern%no_of_recv_procs
      write(log_file_id,*) " - recv from ", grid_comm_pattern%recv(i)%pid, " size=", &
        & grid_comm_pattern%recv(i)%no_of_points
      write(log_file_id,*) " - global indexes=", grid_comm_pattern%recv(i)%global_index(:)
    ENDDO
    
    DO i=1,grid_comm_pattern%no_of_send_procs
      write(log_file_id,*) " - send to ", grid_comm_pattern%send(i)%pid, " size=", &
        & grid_comm_pattern%send(i)%no_of_points
      write(log_file_id,*) " - global indexes=", grid_comm_pattern%send(i)%global_index(:)
    ENDDO
        
    write(log_file_id,*) " === End Communication pattern info for ", TRIM(grid_comm_pattern%name)
    
  END SUBROUTINE print_grid_comm_pattern
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  INTEGER FUNCTION get_recvbuffer_id_of_pid( pid )
    INTEGER, INTENT(IN)       :: pid

    INTEGER :: bfid
    CHARACTER(*), PARAMETER :: method_name = "get_recvbuffer_id_of_pid"
    
    DO bfid = 1, active_recv_buffers
      IF ( recv_procs_buffer(bfid)%pid == pid ) THEN
        get_recvbuffer_id_of_pid = bfid
        RETURN
      ENDIF
    ENDDO

    ! add a new recv_procs_buffer
    IF (active_recv_buffers >= max_no_of_comm_processes ) &
      CALL finish(method_name, "active_recv_buffers >= max_no_of_comm_processes")
    active_recv_buffers = active_recv_buffers + 1
    recv_procs_buffer(active_recv_buffers)%pid = pid
    get_recvbuffer_id_of_pid = active_recv_buffers
    
  END FUNCTION get_recvbuffer_id_of_pid
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  INTEGER FUNCTION get_sendbuffer_id_of_pid( pid )
    INTEGER, INTENT(IN)       :: pid

    INTEGER :: bfid
    CHARACTER(*), PARAMETER :: method_name = "get_sendbuffer_id_of_pid"
    
    DO bfid = 1, active_send_buffers
      IF ( send_procs_buffer(bfid)%pid == pid ) THEN
        get_sendbuffer_id_of_pid = bfid
        RETURN
      ENDIF
    ENDDO

    ! add a new send_procs_buffer
    IF (active_send_buffers >= max_no_of_comm_processes ) &
      CALL finish(method_name, "active_send_buffers >= max_no_of_comm_processes")
    active_send_buffers = active_send_buffers + 1
    send_procs_buffer(active_send_buffers)%pid = pid
    get_sendbuffer_id_of_pid = active_send_buffers
    
  END FUNCTION get_sendbuffer_id_of_pid
  !-----------------------------------------------------------------------
  
  !-----------------------------------------------------------------------
  !>
  !! Creates a new comm_variable and returns its id.
  INTEGER FUNCTION new_comm_variable_r4d(var,  comm_pattern_index, p_patch, vertical_layers, &
    & no_of_variables, status, scope, name )
    INTEGER, INTENT(IN)       :: comm_pattern_index
    TYPE(t_patch), INTENT(IN) :: p_patch
!     REAL(wp), POINTER, INTENT(INOUT)   :: var(:,:,:,:)
    REAL(wp), POINTER  :: var(:,:,:,:)
    
     INTEGER, INTENT(IN), OPTIONAL :: no_of_variables
!     INTEGER, INTENT(IN), OPTIONAL :: var_dim
    INTEGER, INTENT(IN), OPTIONAL :: vertical_layers
!     TYPE(t_comm_pattern), INTENT(IN), POINTER, OPTIONAL :: comm_pattern
    INTEGER, INTENT(IN), OPTIONAL :: status
    INTEGER, INTENT(IN), OPTIONAL :: scope
    CHARACTER(*), INTENT(IN), OPTIONAL :: name
    
    REAL(wp), POINTER  :: var_3d(:,:,:)
    
    CHARACTER(*), PARAMETER :: method_name = "new_comm_variable_r4d"

    
    IF  (this_is_mpi_sequential) THEN
      new_comm_variable_r4d = 0
      RETURN
    ENDIF

    var_3d => var(:,:,:,1)
    IF (PRESENT(vertical_layers)) THEN
      new_comm_variable_r4d = &
        & new_comm_variable_r3d(var_3d,  comm_pattern_index, p_patch, vertical_layers)
    ELSE
      new_comm_variable_r4d = &
        & new_comm_variable_r3d(var_3d,  comm_pattern_index, p_patch)
    ENDIF
    
    comm_variable(new_comm_variable_r4d)%values_4d => var
    NULLIFY(comm_variable(new_comm_variable_r4d)%values_3d)

    IF (PRESENT(no_of_variables)) THEN   
      comm_variable(new_comm_variable_r4d)%no_of_variables = no_of_variables
    ELSE
      comm_variable(new_comm_variable_r4d)%no_of_variables = SIZE(var,4)
    ENDIF

    comm_variable(new_comm_variable_r4d)%dim_4 = &
      comm_variable(new_comm_variable_r4d)%no_of_variables

    IF (PRESENT(status)) THEN
      IF (status == is_ready) CALL icon_comm_var_is_ready(new_comm_variable_r4d)
    ENDIF
    
    IF (PRESENT(scope)) THEN
      comm_variable(new_comm_variable_r4d)%scope = scope
    ELSE
      comm_variable(new_comm_variable_r4d)%scope = global   
    ENDIF
    
    IF (PRESENT(name)) THEN
      comm_variable(new_comm_variable_r4d)%name = TRIM(name)
    ELSE
      comm_variable(new_comm_variable_r4d)%name = ""
    ENDIF
        
  END FUNCTION new_comm_variable_r4d
  !-----------------------------------------------------------------------
  
  !-----------------------------------------------------------------------
  !>
  !! Creates a new comm_variable and returns its id.
!   INTEGER FUNCTION new_comm_variable_r3d_target(var,  comm_pattern_index, p_patch, &
!     & vertical_layers, status, scope, name)
!     
!     INTEGER, INTENT(IN)       :: comm_pattern_index
!     TYPE(t_patch), INTENT(IN) :: p_patch
! !     REAL(wp), POINTER, INTENT(INOUT)   :: var(:,:,:)
!     REAL(wp), TARGET, INTENT(inout) :: var(:,:,:)
!     
! !     INTEGER, INTENT(IN), OPTIONAL :: no_of_variables
! !     INTEGER, INTENT(IN), OPTIONAL :: var_dim
!     INTEGER, INTENT(IN), OPTIONAL :: vertical_layers
! !     TYPE(t_comm_pattern), INTENT(IN), POINTER, OPTIONAL :: comm_pattern
!     INTEGER, INTENT(IN), OPTIONAL :: status    
!     INTEGER, INTENT(IN), OPTIONAL :: scope
!     CHARACTER(*), INTENT(IN), OPTIONAL :: name
!     
!     REAL(wp), POINTER :: p_var(:,:,:)
! 
!     p_var => var
!     new_comm_variable_r3d_target = &
!       & new_comm_variable_r3d(p_var, comm_pattern_index, p_patch, &
!       & vertical_layers, status, scope, name)
!       
!   END FUNCTION new_comm_variable_r3d_target
  !-----------------------------------------------------------------------
  

  !-----------------------------------------------------------------------
  !>
  !! Creates a new comm_variable and returns its id.
  INTEGER FUNCTION new_comm_variable_r3d(var,  comm_pattern_index, p_patch, &
    & vertical_layers, status, scope, name) !, &
   ! & var_dim, no_of_variables, vertical_layers)
    
    INTEGER, INTENT(IN)       :: comm_pattern_index
    TYPE(t_patch), INTENT(IN) :: p_patch
!     REAL(wp), POINTER, INTENT(INOUT)   :: var(:,:,:)
    REAL(wp), POINTER   :: var(:,:,:)
    
!     INTEGER, INTENT(IN), OPTIONAL :: no_of_variables
!     INTEGER, INTENT(IN), OPTIONAL :: var_dim
    INTEGER, INTENT(IN), OPTIONAL :: vertical_layers
!     TYPE(t_comm_pattern), INTENT(IN), POINTER, OPTIONAL :: comm_pattern
    INTEGER, INTENT(IN), OPTIONAL :: status    
    INTEGER, INTENT(IN), OPTIONAL :: scope
    CHARACTER(*), INTENT(IN), OPTIONAL :: name
    
    CHARACTER(*), PARAMETER :: method_name = "new_comm_variable_r3d"
        
    IF(this_is_mpi_sequential) THEN
      new_comm_variable_r3d = 0
      RETURN
    ENDIF

    new_comm_variable_r3d = get_new_comm_variable()

    ! where the variable lives
!     SELECT CASE (comm_pattern_index)
!     CASE ( cells_not_in_domain, cells_one_edge_in_domain, edges_not_owned, verts_not_owned)
!       comm_variable(new_comm_variable_r3d)%comm_pattern_index = comm_pattern_index
!    CASE default
!       CALL finish(method_name, "Unrecoginzed comm_pattern_index")
!     END SELECT
    
    comm_variable(new_comm_variable_r3d)%comm_pattern_index = comm_pattern_index

    ! check if comm_pattern is initialized
    IF ( p_patch%id > max_no_of_patches ) &
      CALL finish(method_name, "p_patch%id > max_no_of_patches")
    IF ( grid_comm_pattern_list(comm_pattern_index, p_patch%id)%status == not_active ) &
      CALL finish(method_name, "grid_comm_pattern status = not_active")
    
    comm_variable(new_comm_variable_r3d)%grid_comm_pattern => &
        & grid_comm_pattern_list(comm_pattern_index, p_patch%id)
    
    ! this is for a 3D variable
    comm_variable(new_comm_variable_r3d)%grid_dim = grid_3D

    ! check the dim_1
!     IF( SIZE(var,1) /= nproma ) THEN
!        CALL finish(method_name, 'SIZE(var,1) /= nproma')
!     ENDIF
!     comm_variable(new_comm_variable_r3d)%dim_1 = SIZE(var,1)
    
    ! check the vertical_layers
    comm_variable(new_comm_variable_r3d)%vertical_layers =  SIZE(var,LEVELS_POSITION)
    IF ( PRESENT(vertical_layers) ) THEN
      IF ( vertical_layers <=  comm_variable(new_comm_variable_r3d)%vertical_layers) THEN
         comm_variable(new_comm_variable_r3d)%vertical_layers = vertical_layers
      ELSE
         CALL finish(method_name, "vertical_layers are greater than SIZE(var,2)")
      END IF
    END IF    

    comm_variable(new_comm_variable_r3d)%values_3d => var
    comm_variable(new_comm_variable_r3d)%no_of_variables = 1
    
    IF (PRESENT(status)) THEN
      IF (status == is_ready) CALL icon_comm_var_is_ready(new_comm_variable_r3d )
    ENDIF
    
    IF (PRESENT(scope)) THEN
      comm_variable(new_comm_variable_r3d)%scope = scope
    ELSE
      comm_variable(new_comm_variable_r3d)%scope = global
    ENDIF
    
    IF (PRESENT(name)) THEN
      comm_variable(new_comm_variable_r3d)%name = TRIM(name)
    ELSE
      comm_variable(new_comm_variable_r3d)%name = ""      
    ENDIF

  END FUNCTION new_comm_variable_r3d
  !-----------------------------------------------------------------------
  
  !-----------------------------------------------------------------------
  !>
  !! Creates a new comm_variable and returns its id.
  INTEGER FUNCTION new_comm_variable_r2d(var,  comm_pattern_index, p_patch, &
    & vertical_layers, status, scope, name) !, &
   ! & var_dim, no_of_variables, vertical_layers)
    
    INTEGER, INTENT(IN)       :: comm_pattern_index
    TYPE(t_patch), INTENT(IN) :: p_patch
!     REAL(wp), POINTER, INTENT(INOUT)   :: var(:,:)
    REAL(wp), POINTER   :: var(:,:)
    
!     INTEGER, INTENT(IN), OPTIONAL :: no_of_variables
!     INTEGER, INTENT(IN), OPTIONAL :: var_dim
     INTEGER, INTENT(IN), OPTIONAL :: vertical_layers
!     TYPE(t_comm_pattern), INTENT(IN), POINTER, OPTIONAL :: comm_pattern
    INTEGER, INTENT(IN), OPTIONAL :: status    
    INTEGER, INTENT(IN), OPTIONAL :: scope
    CHARACTER(*), INTENT(IN), OPTIONAL :: name
    
    CHARACTER(*), PARAMETER :: method_name = "new_comm_variable_r2d"
        
    IF(this_is_mpi_sequential) THEN
      new_comm_variable_r2d = 0
      RETURN
    ENDIF

    new_comm_variable_r2d = get_new_comm_variable()

    ! where the variable lives
!     SELECT CASE (comm_pattern_index)
!     CASE ( cells_not_in_domain, cells_one_edge_in_domain, edges_not_owned, verts_not_owned)
!       comm_variable(new_comm_variable_r2d)%comm_pattern_index = comm_pattern_index
!     CASE default
!       CALL finish(method_name, "Unrecoginzed comm_pattern_index")
!     END SELECT
    comm_variable(new_comm_variable_r2d)%comm_pattern_index = comm_pattern_index

    ! check if comm_pattern is initialized
    IF ( p_patch%id > max_no_of_patches ) &
      CALL finish(method_name, "p_patch%id > max_no_of_patches")
    IF ( grid_comm_pattern_list(comm_pattern_index, p_patch%id)%status == not_active ) &
      CALL finish(method_name, "grid_comm_pattern status = not_active")
    
    comm_variable(new_comm_variable_r2d)%grid_comm_pattern => &
        & grid_comm_pattern_list(comm_pattern_index, p_patch%id)
    
    ! this is for a 3D variable
    comm_variable(new_comm_variable_r2d)%grid_dim = grid_2D

    ! check the dim_1
!     IF( SIZE(var,1) /= nproma ) THEN
!        CALL finish(method_name, 'SIZE(var,1) /= nproma')
!     ENDIF
!     comm_variable(new_comm_variable_r2d)%dim_1 = nproma
    
    ! check the vertical_layers
    comm_variable(new_comm_variable_r2d)%vertical_layers = 1

    comm_variable(new_comm_variable_r2d)%values_2d => var
    comm_variable(new_comm_variable_r2d)%no_of_variables = 1
    
    IF (PRESENT(status)) THEN
      IF (status == is_ready) CALL icon_comm_var_is_ready(new_comm_variable_r2d )
    ENDIF
    
    IF (PRESENT(scope)) THEN
      comm_variable(new_comm_variable_r2d)%scope = scope
    ELSE
      comm_variable(new_comm_variable_r2d)%scope = global
    ENDIF
    
    IF (PRESENT(name)) THEN
      comm_variable(new_comm_variable_r2d)%name = TRIM(name)
    ELSE
      comm_variable(new_comm_variable_r2d)%name = ""
    ENDIF

  END FUNCTION new_comm_variable_r2d
  !-----------------------------------------------------------------------
  

  !-----------------------------------------------------------------------
  !>
  !! Returns the id of a new comm_variable
  INTEGER FUNCTION get_new_comm_variable()

    INTEGER :: i

    CHARACTER(*), PARAMETER :: method_name = "get_new_comm_variable"
    
#ifdef _OPENMP
    IF (omp_in_parallel()) &
      CALL finish(method_name, 'cannot be called from openmp parallel');
#endif
    
    IF ( .NOT. comm_lib_is_initialized ) &
      CALL finish(method_name, 'comm_lib is not initialized')

    IF (max_active_comm_variables /= active_comm_variables) THEN
      ! we have a hole (inactive grid) in the list of max_active_grids
      DO i = 1, max_active_comm_variables
        IF ( comm_variable(i)%status == not_active ) THEN
          get_new_comm_variable = i
          EXIT
        ENDIF
      ENDDO
    ELSE
      ! check if we reached the end 
      IF (active_comm_variables >= max_no_of_comm_variables) THEN
        CALL finish(method_name, 'exceeded max_no_of_comm_variables');
      ENDIF
      max_active_comm_variables = max_active_comm_variables + 1
      get_new_comm_variable = max_active_comm_variables
    ENDIF

    active_comm_variables = active_comm_variables + 1
    
    comm_variable(get_new_comm_variable)%status = active
    
  END FUNCTION get_new_comm_variable
  !-----------------------------------------------------------------------
  
  !-----------------------------------------------------------------------
  !>
  !! Deletes the comm_variable
  SUBROUTINE delete_icon_comm_variable(comm_variable_id)
    INTEGER, INTENT(in) :: comm_variable_id
    
    CHARACTER(*), PARAMETER :: method_name = "delete_icon_comm_variable"

!     CALL check_active_comm_variable(comm_variable_id)
#ifdef _OPENMP
    IF (omp_in_parallel()) &
      CALL finish(method_name, 'cannot be called from openmp parallel');
#endif
    IF(this_is_mpi_sequential) RETURN

    comm_variable(comm_variable_id)%status = not_active
    active_comm_variables = active_comm_variables - 1
    IF ( comm_variable_id == max_active_comm_variables ) &
      max_active_comm_variables = max_active_comm_variables - 1

    CALL clear_comm_variable(comm_variable_id)
    
  END SUBROUTINE delete_icon_comm_variable
  !-----------------------------------------------------------------------
  
  !-----------------------------------------------------------------------
  SUBROUTINE clear_comm_variable(id)
    INTEGER, INTENT(in) :: id
       
    comm_variable(id)%status = not_active
    comm_variable(id)%request = not_active
    comm_variable(id)%comm_status = not_active
    comm_variable(id)%comm_buffer = 0
    
    NULLIFY(comm_variable(id)%p_patch)  
          
    comm_variable(id)%no_of_variables = 0
!     comm_variable(id)%dim_1 = 0
    comm_variable(id)%vertical_layers = 0
!     comm_variable(id)%dim_3 = 0
    comm_variable(id)%dim_4 = 0
          
    NULLIFY(comm_variable(id)%values_2d)
    NULLIFY(comm_variable(id)%values_3d)
    NULLIFY(comm_variable(id)%values_4d)
  
  END SUBROUTINE clear_comm_variable
  !-----------------------------------------------------------------------
  
  
  !-----------------------------------------------------------------------
  !>
  !! The comm_variable can be communicated
  SUBROUTINE icon_comm_var_is_ready(comm_variable_id)
    INTEGER, INTENT(in) :: comm_variable_id
    
     IF(this_is_mpi_sequential) RETURN

!     CALL check_active_comm_variable(comm_variable_id)
!$OMP SINGLE
      comm_variable(comm_variable_id)%request = communicate
!$OMP END SINGLE
    
  END SUBROUTINE icon_comm_var_is_ready
  !-----------------------------------------------------------------------
  
  !-----------------------------------------------------------------------
  !>
  SUBROUTINE icon_comm_sync_2D_1(var,  comm_pattern_index, patch)
    INTEGER, INTENT(IN)       :: comm_pattern_index
    TYPE(t_patch), INTENT(IN) :: patch
!    REAL(wp), POINTER, INTENT(INOUT)   :: var(:,:,:)
    REAL(wp), POINTER   :: var(:,:)
    
    INTEGER :: comm_var

    IF(this_is_mpi_sequential) RETURN

    comm_var = new_icon_comm_variable(var,  comm_pattern_index, patch, &
      & status=is_ready, scope=until_sync )
    CALL icon_comm_sync_all
 
  END SUBROUTINE icon_comm_sync_2D_1
  !-----------------------------------------------------------------------
  
  !-----------------------------------------------------------------------
  !>
  SUBROUTINE icon_comm_sync_3D_1(var,  comm_pattern_index, patch)
    INTEGER, INTENT(IN)       :: comm_pattern_index
    TYPE(t_patch), INTENT(IN) :: patch
!    REAL(wp), POINTER, INTENT(INOUT)   :: var(:,:,:)
    REAL(wp), POINTER   :: var(:,:,:)
    
    INTEGER :: comm_var

    IF(this_is_mpi_sequential) RETURN

    comm_var = new_icon_comm_variable(var,  comm_pattern_index, patch, &
      & status=is_ready, scope=until_sync )
    CALL icon_comm_sync_all
 
  END SUBROUTINE icon_comm_sync_3D_1
  !-----------------------------------------------------------------------
  
  !-----------------------------------------------------------------------
  !>
  SUBROUTINE icon_comm_sync_3D_2(var1,  var2, comm_pattern_index, patch)
    INTEGER, INTENT(IN)       :: comm_pattern_index
    TYPE(t_patch), INTENT(IN) :: patch
!     REAL(wp), POINTER, INTENT(INOUT)   :: var1(:,:,:)
!     REAL(wp), POINTER, INTENT(INOUT)   :: var2(:,:,:)
    REAL(wp), POINTER   :: var1(:,:,:)
    REAL(wp), POINTER   :: var2(:,:,:)
    
    INTEGER :: comm_var_1, comm_var_2

    IF(this_is_mpi_sequential) RETURN

    comm_var_1 = new_icon_comm_variable(var1,  comm_pattern_index, patch, &
      & status=is_ready, scope=until_sync)
    comm_var_2 = new_icon_comm_variable(var2,  comm_pattern_index, patch,&
      & status=is_ready, scope=until_sync)
    CALL icon_comm_sync_all
 
  END SUBROUTINE icon_comm_sync_3D_2
  !-----------------------------------------------------------------------
          
  !-----------------------------------------------------------------------
  !>
  SUBROUTINE icon_comm_sync_3D_3(var1,  var2, var3, comm_pattern_index, patch)
    INTEGER, INTENT(IN)       :: comm_pattern_index
    TYPE(t_patch), INTENT(IN) :: patch
!     REAL(wp), POINTER, INTENT(INOUT)   :: var1(:,:,:)
!     REAL(wp), POINTER, INTENT(INOUT)   :: var2(:,:,:)
    REAL(wp), POINTER   :: var1(:,:,:)
    REAL(wp), POINTER   :: var2(:,:,:)
    REAL(wp), POINTER   :: var3(:,:,:)
    
    INTEGER :: comm_var_1, comm_var_2, comm_var_3

    IF(this_is_mpi_sequential) RETURN

    comm_var_1 = new_icon_comm_variable(var1,  comm_pattern_index, patch, &
      & status=is_ready, scope=until_sync)
    comm_var_2 = new_icon_comm_variable(var2,  comm_pattern_index, patch,&
      & status=is_ready, scope=until_sync)
    comm_var_3 = new_icon_comm_variable(var3,  comm_pattern_index, patch,&
      & status=is_ready, scope=until_sync)
    CALL icon_comm_sync_all
 
  END SUBROUTINE icon_comm_sync_3D_3
  !-----------------------------------------------------------------------
        
  
  !-----------------------------------------------------------------------
  SUBROUTINE icon_comm_sync_all()

    LOGICAL :: exist_communication_var
    INTEGER :: comm_var

    CHARACTER(*), PARAMETER :: method_name = "icon_comm_sync_all"
    
    
#ifdef _OPENMP
    IF (omp_in_parallel()) &
      CALL finish(method_name, 'cannot be called from openmp parallel');
#endif
   
    IF(this_is_mpi_sequential) RETURN

    ! check if we have any variables to communicate
    exist_communication_var = .false.
    DO comm_var = 1, max_active_comm_variables     
      IF ( comm_variable(comm_var)%request == communicate ) THEN
        exist_communication_var = .true.
        EXIT
      ENDIF
    ENDDO
    
    IF (.NOT. exist_communication_var) RETURN
    
    IF (ltimer) CALL timer_start(timer_icon_comm_sync)
    
    IF (buffer_comm_status == not_active) THEN
      ! no communication steps have been taken
      ! go through the whole process
      SELECT CASE(icon_comm_method)
      
      CASE(1,2)
        CALL start_recv_active_buffers() 
        CALL fill_send_buffers()
        CALL sent_active_buffers()
        CALL finalize_recv_active_buffers()
        
      CASE(3,4)
        CALL start_recv_active_buffers() 
        CALL fill_and_send_buffers
        CALL finalize_recv_active_buffers()

      CASE DEFAULT
        CALL finish( method_name,'unknown icon_comm_method.')
      END SELECT
      
      CALL clear_comm()

    ELSE
      
      CALL finish(method_name, 'buffer_comm_status /= not_active');

    ENDIF
    
    IF (ltimer) CALL timer_stop(timer_icon_comm_sync)

!     CALL p_barrier(process_mpi_all_comm)

!    CALL p_stop
    
  END SUBROUTINE icon_comm_sync_all
  !-----------------------------------------------------------------------
  
  !-----------------------------------------------------------------------
  !>
  SUBROUTINE clear_comm()
    
    INTEGER :: comm_var! ,  bfid
    
    buffer_comm_status = not_active
    
    IF(this_is_mpi_sequential) RETURN

    ! clear buffer sizes
!     DO bfid = 1, active_send_buffers
!       send_procs_buffer(bfid)%buffer_size = 0
!     ENDDO

    ! go through all variables and clean comm requests
    DO comm_var = 1, max_active_comm_variables
      IF (comm_variable(comm_var)%scope == until_sync) THEN
         CALL delete_icon_comm_variable(comm_var)
      ELSE
        comm_variable(comm_var)%request     = not_active
        comm_variable(comm_var)%comm_status = not_active
      ENDIF
    ENDDO
      
    CALL adjust_max_active_variables()

  END SUBROUTINE clear_comm
  !-----------------------------------------------------------------------
  
  !-----------------------------------------------------------------------
  !>
  SUBROUTINE adjust_max_active_variables()
    INTEGER :: comm_var, new_max_active_variables
    
    new_max_active_variables=0
    DO comm_var = 1, max_active_comm_variables
      IF ( comm_variable(comm_var)%status /= not_active ) &
        new_max_active_variables = comm_var
    ENDDO
    max_active_comm_variables = new_max_active_variables
    
  END SUBROUTINE adjust_max_active_variables
  !-----------------------------------------------------------------------
  
  
  !-----------------------------------------------------------------------
  SUBROUTINE fill_send_buffers()
        
    INTEGER :: comm_var, var_no    
!     CHARACTER(*), PARAMETER :: method_name = "compute_send_buffer_sizes"

    IF (activate_sync_timers) CALL timer_start(timer_icon_comm_fillsend)
    
    CALL compute_send_buffer_sizes()

    DO comm_var = 1, max_active_comm_variables
      ! 
      IF ( comm_variable(comm_var)%request /= communicate ) CYCLE

        DO var_no = 1, comm_variable(comm_var)%no_of_variables
          CALL fill_send_buffers_var(comm_var, var_no)
        END DO  
      
    END DO  
    IF (activate_sync_timers) CALL timer_stop(timer_icon_comm_fillsend)
       
  END SUBROUTINE fill_send_buffers
  !-----------------------------------------------------------------------
  
  !-----------------------------------------------------------------------  
  SUBROUTINE fill_and_send_buffers()
        
    TYPE(t_grid_comm_pattern), POINTER :: grid_comm_pattern
    INTEGER :: comm_var, var_no, bfid, np    
!     CHARACTER(*), PARAMETER :: method_name = "compute_send_buffer_sizes"

    IF (activate_sync_timers) CALL timer_start(timer_icon_comm_fillandsend)
    
    CALL compute_send_buffer_sizes()    

!ICON_OMP_PARALLEL IF(icon_comm_openmp)
!ICON_OMP_DO PRIVATE(bfid, comm_var, grid_comm_pattern, np, var_no)
    DO bfid = 1, active_send_buffers
      IF ( send_procs_buffer(bfid)%buffer_size == 0) CYCLE

      DO comm_var = 1, max_active_comm_variables
        IF ( comm_variable(comm_var)%request /= communicate ) CYCLE
        grid_comm_pattern => comm_variable(comm_var)%grid_comm_pattern        
!         vertical_layers = comm_variable(comm_var)%vertical_layers
        
        ! go through the requested send pids
        ! and choose the requested
        DO np = 1, grid_comm_pattern%no_of_send_procs ! loop over PEs where to send the data
          IF (bfid /= grid_comm_pattern%send(np)%buffer_index) CYCLE

          DO var_no = 1, comm_variable(comm_var)%no_of_variables
            CALL fill_onepid_send_buffer_var(bfid, comm_var, np, var_no)
           END DO
      
        ENDDO 
      
      END DO ! comm_var = 1, max_active_comm_variables
      
      CALL p_isend(send_buffer(send_procs_buffer(bfid)%start_index:), &
        & send_procs_buffer(bfid)%pid, p_tag=halo_tag,                &
        &  p_count=send_procs_buffer(bfid)%buffer_size,               &
        & comm=my_work_communicator)
            
    ENDDO !bfid = 1, active_send_buffers
!ICON_OMP_END_DO_NOWAIT
!ICON_OMP_END_PARALLEL
    IF (activate_sync_timers) CALL timer_stop(timer_icon_comm_fillandsend)
           
  END SUBROUTINE fill_and_send_buffers
  !-----------------------------------------------------------------------
  
  !-----------------------------------------------------------------------
  !>
  !! Computes the sizes and the indexes for buffers to be sent
  SUBROUTINE compute_send_buffer_sizes()
    
    TYPE(t_grid_comm_pattern), POINTER :: grid_comm_pattern
    
    INTEGER :: comm_var, vertical_layers, np, bfid, total_send_size, buffer_start
    INTEGER :: no_of_variables
    
    CHARACTER(*), PARAMETER :: method_name = "compute_send_buffer_sizes"

    ! clear buffer sizes
    DO bfid = 1, active_send_buffers
      send_procs_buffer(bfid)%buffer_size = 0
    ENDDO

    ! go through all variables and compute the requested send sizes
    DO comm_var = 1, max_active_comm_variables
      ! 
      IF ( comm_variable(comm_var)%request /= communicate ) CYCLE

      grid_comm_pattern => comm_variable(comm_var)%grid_comm_pattern
      vertical_layers = comm_variable(comm_var)%vertical_layers
      no_of_variables = comm_variable(comm_var)%no_of_variables
      
      ! go through the requested send pids
      DO np = 1, grid_comm_pattern%no_of_send_procs ! loop over PEs where to send the data
        
        total_send_size = &
          & ((grid_comm_pattern%send(np)%no_of_points * vertical_layers) &
          & + checksum_size) * no_of_variables
          
        bfid = grid_comm_pattern%send(np)%buffer_index
        send_procs_buffer(bfid)%buffer_size = send_procs_buffer(bfid)%buffer_size &
          &  + total_send_size
        
      ENDDO
      
    END DO  
                  
    ! now compute the buffers' start, current and end indexes
    buffer_start = 1
    DO bfid = 1, active_send_buffers
      send_procs_buffer(bfid)%start_index = buffer_start
      buffer_start = buffer_start + send_procs_buffer(bfid)%buffer_size
      send_procs_buffer(bfid)%end_index = buffer_start - 1
      ! set the current_index to start
      send_procs_buffer(bfid)%current_index = send_procs_buffer(bfid)%start_index
    ENDDO
    
    ! check if we went over the max_send_recv_buffer_size
    IF (buffer_start >= max_send_recv_buffer_size) &
      CALL finish(method_name, "buffer_start >= max_send_recv_buffer_size")
 
  END SUBROUTINE compute_send_buffer_sizes
  !-----------------------------------------------------------------------
  
  !-----------------------------------------------------------------------
  SUBROUTINE fill_send_buffers_var(comm_var, var_no)
    INTEGER, INTENT(in) :: comm_var, var_no
    
    TYPE(t_grid_comm_pattern), POINTER :: grid_comm_pattern
    TYPE(t_patch), POINTER :: p_patch
    REAL(wp), POINTER :: send_var_2d(:,:)
    REAL(wp), GENERAL_3D, POINTER :: send_var_3d
    
    INTEGER :: vertical_layers, np, bfid, var_send_size, current_buffer_index
    INTEGER :: i, k
    
!     CHARACTER(*), PARAMETER :: method_name = "fill_send_buffers_var"

    IF ( comm_variable(comm_var)%request /= communicate ) RETURN

    p_patch => comm_variable(comm_var)%p_patch
    grid_comm_pattern => comm_variable(comm_var)%grid_comm_pattern
    vertical_layers = comm_variable(comm_var)%vertical_layers
    IF (comm_variable(comm_var)%grid_dim == grid_2D ) THEN
      send_var_2d => comm_variable(comm_var)%values_2d(:,:)
    ELSEIF ( comm_variable(comm_var)%dim_4 > 0 ) THEN
      send_var_3d => comm_variable(comm_var)%values_4d(:,:,:,var_no)
    ELSE    
      send_var_3d => comm_variable(comm_var)%values_3d
    ENDIF

    ! go through the requested send pids
    DO np = 1, grid_comm_pattern%no_of_send_procs ! loop over PEs where to send the data

      bfid = grid_comm_pattern%send(np)%buffer_index
      current_buffer_index = send_procs_buffer(bfid)%current_index
      var_send_size = grid_comm_pattern%send(np)%no_of_points * vertical_layers
      ! fill the header
      ! this is the id of the variabkle and the size we sent
      send_buffer(current_buffer_index) = REAL(comm_var,wp)
      current_buffer_index = current_buffer_index + 1
      send_buffer(current_buffer_index) = REAL(var_send_size,wp)
      current_buffer_index = current_buffer_index + 1
      
      ! fill the buffer
      IF (comm_variable(comm_var)%grid_dim == grid_2D ) THEN
        ! fill 2d
        DO i = 1, grid_comm_pattern%send(np)%no_of_points
          send_buffer(current_buffer_index) = send_var_2d &
            & ( grid_comm_pattern%send(np)%index_no(i),   &
                grid_comm_pattern%send(np)%block_no(i) )

            current_buffer_index = current_buffer_index + 1
        ENDDO
        
        IF (icon_comm_debug) THEN
          DO i = 1, grid_comm_pattern%send(np)%no_of_points
            write(log_file_id,*) TRIM(comm_variable(comm_var)%name), " sent to ", &
              & send_procs_buffer(bfid)%pid, ":", &
              & i, send_var_2d( grid_comm_pattern%send(np)%index_no(i), &
                    grid_comm_pattern%send(np)%block_no(i) )
          ENDDO
        ENDIF

      ELSE
      
        ! fill 3d
!         DO k = 1, vertical_layers
        DO i = 1, grid_comm_pattern%send(np)%no_of_points
            send_buffer(current_buffer_index : (current_buffer_index + vertical_layers - 1)) = &
               send_var_3d ( grid_comm_pattern%send(np)%index_no(i), 1:vertical_layers, &
                  grid_comm_pattern%send(np)%block_no(i) )

            current_buffer_index = current_buffer_index + vertical_layers
!             current_buffer_index = current_buffer_index + 1            
        ENDDO
!         ENDDO

        IF (icon_comm_debug) THEN
          k=1
          DO i = 1, grid_comm_pattern%send(np)%no_of_points
            write(log_file_id,*) TRIM(comm_variable(comm_var)%name), " sent to ", &
              & send_procs_buffer(bfid)%pid, ":", &
              & i,k, send_var_3d( grid_comm_pattern%send(np)%index_no(i), k, &
                    grid_comm_pattern%send(np)%block_no(i) )
          ENDDO                
        ENDIF
       
      ENDIF
        
      send_procs_buffer(bfid)%current_index = current_buffer_index

    ENDDO
      
  END SUBROUTINE fill_send_buffers_var
  !-----------------------------------------------------------------------
 
  !-----------------------------------------------------------------------
  SUBROUTINE fill_onepid_send_buffer_var(bfid, comm_var, send_proc_index, var_no)
    INTEGER, INTENT(in) :: bfid, comm_var, send_proc_index, var_no
    
    TYPE(t_grid_comm_pattern), POINTER :: grid_comm_pattern
    TYPE(t_patch), POINTER :: p_patch
    REAL(wp), POINTER :: send_var_2d(:,:)
    REAL(wp), GENERAL_3D, POINTER :: send_var_3d
    
    INTEGER :: vertical_layers, np, var_send_size, current_buffer_index
    INTEGER :: i, k
    
!     CHARACTER(*), PARAMETER :: method_name = "fill_send_buffers_var"

    p_patch => comm_variable(comm_var)%p_patch
    grid_comm_pattern => comm_variable(comm_var)%grid_comm_pattern
    vertical_layers = comm_variable(comm_var)%vertical_layers
    IF (comm_variable(comm_var)%grid_dim == grid_2D ) THEN
      send_var_2d => comm_variable(comm_var)%values_2d(:,:)
    ELSEIF ( comm_variable(comm_var)%dim_4 > 0 ) THEN
      send_var_3d => comm_variable(comm_var)%values_4d(:,:,:,var_no)
    ELSE    
      send_var_3d => comm_variable(comm_var)%values_3d
    ENDIF

    ! go through the requested send pids
    np = send_proc_index
      
    current_buffer_index = send_procs_buffer(bfid)%current_index
    var_send_size = grid_comm_pattern%send(np)%no_of_points * vertical_layers
    ! fill the header
    ! this is the id of the variabkle and the size we sent
    send_buffer(current_buffer_index) = REAL(comm_var,wp)
    send_buffer(current_buffer_index + 1) = REAL(var_send_size,wp)
    current_buffer_index = current_buffer_index + 2

    ! fill the buffer
    IF (comm_variable(comm_var)%grid_dim == grid_2D ) THEN
      ! fill 2d
      DO i = 1, grid_comm_pattern%send(np)%no_of_points
        send_buffer(current_buffer_index) = send_var_2d &
          & ( grid_comm_pattern%send(np)%index_no(i),   &
              grid_comm_pattern%send(np)%block_no(i) )

          current_buffer_index = current_buffer_index + 1
      ENDDO

      IF (icon_comm_debug) THEN
        DO i = 1, grid_comm_pattern%send(np)%no_of_points
          write(log_file_id,*) TRIM(comm_variable(comm_var)%name), " sent to ", &
            & send_procs_buffer(bfid)%pid, ":", &
            & i, send_var_2d( grid_comm_pattern%send(np)%index_no(i), &
                  grid_comm_pattern%send(np)%block_no(i) )
        ENDDO
      ENDIF

    ELSE

      ! fill 3d
      DO i = 1, grid_comm_pattern%send(np)%no_of_points
          send_buffer(current_buffer_index : (current_buffer_index + vertical_layers - 1)) = &
            send_var_3d ( grid_comm_pattern%send(np)%index_no(i), 1:vertical_layers, &
                grid_comm_pattern%send(np)%block_no(i) )

          current_buffer_index = current_buffer_index + vertical_layers
      ENDDO

      IF (icon_comm_debug) THEN
        k=1
        DO i = 1, grid_comm_pattern%send(np)%no_of_points
          write(log_file_id,*) TRIM(comm_variable(comm_var)%name), " sent to ", &
            & send_procs_buffer(bfid)%pid, ":", &
            & i,k, send_var_3d( grid_comm_pattern%send(np)%index_no(i), k, &
                  grid_comm_pattern%send(np)%block_no(i) )
        ENDDO
      ENDIF

    ENDIF ! (comm_variable(comm_var)%grid_dim == grid_2D )

    send_procs_buffer(bfid)%current_index = current_buffer_index

      
  END SUBROUTINE fill_onepid_send_buffer_var
  !-----------------------------------------------------------------------
 
  !-----------------------------------------------------------------------
  !>
  SUBROUTINE sent_active_buffers(  )

    INTEGER :: bfid, buffer_start, buffer_size
!     CHARACTER(*), PARAMETER :: method_name = "sent_all_buffers"

    IF (activate_sync_timers) CALL timer_start(timer_icon_comm_isend)
    
    DO bfid = 1, active_send_buffers

      buffer_start = send_procs_buffer(bfid)%start_index
      buffer_size  = send_procs_buffer(bfid)%buffer_size
      
      IF ( buffer_size > 0 ) THEN
        CALL p_isend(send_buffer(buffer_start:), send_procs_buffer(bfid)%pid, &
          & p_tag=halo_tag, p_count=buffer_size, comm=my_work_communicator)
      ENDIF

    ENDDO
    
    IF (activate_sync_timers) CALL timer_stop(timer_icon_comm_isend)
          
  END SUBROUTINE sent_active_buffers
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  SUBROUTINE start_recv_active_buffers(  )

    INTEGER :: bfid, buffer_start, buffer_size

    CALL compute_recv_buffer_sizes()
    
    IF (activate_sync_timers) CALL timer_start(timer_icon_comm_ircv)
     
    ! start recieve
    DO bfid = 1, active_recv_buffers

      buffer_start = recv_procs_buffer(bfid)%start_index
      buffer_size  = recv_procs_buffer(bfid)%buffer_size
      
      IF ( buffer_size > 0 ) THEN
        CALL p_irecv(recv_buffer(buffer_start:), recv_procs_buffer(bfid)%pid, &
          & p_tag=halo_tag, p_count=buffer_size, comm=my_work_communicator)
      ENDIF

    ENDDO
    IF (activate_sync_timers) CALL timer_stop(timer_icon_comm_ircv)
          
  END SUBROUTINE start_recv_active_buffers
  !-----------------------------------------------------------------------
  
  !-----------------------------------------------------------------------
  !>
  !! set the current recv index to start
  SUBROUTINE set_recv_current_index()
    INTEGER :: bfid
    DO bfid = 1, active_recv_buffers
      recv_procs_buffer(bfid)%current_index = recv_procs_buffer(bfid)%start_index
    ENDDO
  END SUBROUTINE set_recv_current_index
  !-----------------------------------------------------------------------

  
  !-----------------------------------------------------------------------
  !>
  SUBROUTINE finalize_recv_active_buffers()

    INTEGER :: comm_var, var_no, bfid
        
    
    IF (activate_sync_timers) CALL timer_start(timer_icon_comm_wait)
    ! Wait for all outstanding requests to finish
    CALL p_wait
    IF (activate_sync_timers) THEN
      CALL timer_stop(timer_icon_comm_wait)
      CALL timer_start(timer_icon_comm_fillrecv)
    ENDIF
    
    ! set the current recv index to start
    DO bfid = 1, active_recv_buffers
      recv_procs_buffer(bfid)%current_index = recv_procs_buffer(bfid)%start_index
    ENDDO
    
    ! go through all active variables and fill from buffer
    DO comm_var = 1, max_active_comm_variables
      IF ( comm_variable(comm_var)%request /= communicate ) CYCLE
    
        DO var_no = 1, comm_variable(comm_var)%no_of_variables
          CALL fill_var_from_recv_buffers(comm_var, var_no)
        END DO

    ENDDO
    IF (activate_sync_timers) CALL timer_stop(timer_icon_comm_fillrecv)
  
  END SUBROUTINE finalize_recv_active_buffers
  !-----------------------------------------------------------------------
  
  !-----------------------------------------------------------------------
  SUBROUTINE fill_var_from_recv_buffers(comm_var, var_no)
    INTEGER, INTENT(in) :: comm_var, var_no
    
    TYPE(t_grid_comm_pattern), POINTER :: grid_comm_pattern
    TYPE(t_patch), POINTER :: p_patch
    REAL(wp), POINTER :: recv_var_2d(:,:)
    REAL(wp), GENERAL_3D, POINTER :: recv_var_3d
    
    INTEGER :: vertical_layers, np, bfid, var_recv_size, current_buffer_index
    INTEGER :: recv_var, recv_size
    INTEGER :: i, k
    
    CHARACTER(*), PARAMETER :: method_name = "fill_var_from_recv_buffers"

    IF ( comm_variable(comm_var)%request /= communicate ) RETURN

    p_patch => comm_variable(comm_var)%p_patch
    grid_comm_pattern => comm_variable(comm_var)%grid_comm_pattern
    vertical_layers = comm_variable(comm_var)%vertical_layers
    IF (comm_variable(comm_var)%grid_dim == grid_2D ) THEN
      recv_var_2d => comm_variable(comm_var)%values_2d(:,:)
    ELSEIF ( comm_variable(comm_var)%dim_4 > 0 ) THEN
      recv_var_3d => comm_variable(comm_var)%values_4d(:,:,:,var_no)
    ELSE
      recv_var_3d => comm_variable(comm_var)%values_3d
    ENDIF

    ! go through the requested recieve pids
    DO np = 1, grid_comm_pattern%no_of_recv_procs ! loop over PEs from where to recieve the data

      bfid = grid_comm_pattern%recv(np)%buffer_index
      current_buffer_index = recv_procs_buffer(bfid)%current_index
      var_recv_size = grid_comm_pattern%recv(np)%no_of_points * vertical_layers
      ! fill the header
      ! this is the id of the variable and the size we sent
      recv_var = INT(recv_buffer(current_buffer_index))
      current_buffer_index = current_buffer_index + 1
      recv_size = INT(recv_buffer(current_buffer_index))
      current_buffer_index = current_buffer_index + 1
      IF ( recv_var /= comm_var ) &
        CALL finish(method_name, 'recv_var /= comm_var')
      IF ( recv_size /= var_recv_size ) &
        CALL finish(method_name, 'recv_size /= var_recv_size')
      
      ! fill from the buffer
      IF (comm_variable(comm_var)%grid_dim == grid_2D ) THEN
         ! fill 2d
        DO i = 1, grid_comm_pattern%recv(np)%no_of_points
          recv_var_2d &
            & ( grid_comm_pattern%recv(np)%index_no(i),    &
                grid_comm_pattern%recv(np)%block_no(i) ) = &
          recv_buffer(current_buffer_index)

          current_buffer_index = current_buffer_index + 1
        ENDDO
        
        IF (icon_comm_debug) THEN
          DO i = 1, grid_comm_pattern%recv(np)%no_of_points
            write(log_file_id,*) TRIM(comm_variable(comm_var)%name), " recv from ", &
              & recv_procs_buffer(bfid)%pid, ":", &
              & i, recv_var_2d( grid_comm_pattern%recv(np)%index_no(i), &
                    grid_comm_pattern%recv(np)%block_no(i) )
          ENDDO
        ENDIF
        
      
      ELSE
        ! fill 3d
!         DO k = 1, vertical_layers
        DO i = 1, grid_comm_pattern%recv(np)%no_of_points
            recv_var_3d( grid_comm_pattern%recv(np)%index_no(i), 1:vertical_layers, &
                  grid_comm_pattern%recv(np)%block_no(i) ) = &
            recv_buffer(current_buffer_index : (current_buffer_index + vertical_layers - 1))
            
            current_buffer_index = current_buffer_index + vertical_layers
!             current_buffer_index = current_buffer_index + 1
!           ENDDO
        ENDDO
        
        IF (icon_comm_debug) THEN
          k=1
          DO i = 1, grid_comm_pattern%recv(np)%no_of_points
            write(log_file_id,*) TRIM(comm_variable(comm_var)%name), " recv from ", &
              & recv_procs_buffer(bfid)%pid, ":", &
              & i,k, recv_var_3d( grid_comm_pattern%recv(np)%index_no(i), k, &
                    grid_comm_pattern%recv(np)%block_no(i) )
          ENDDO
        ENDIF
        
      ENDIF
        
      recv_procs_buffer(bfid)%current_index = current_buffer_index

    ENDDO
      
  END SUBROUTINE fill_var_from_recv_buffers
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  !! Computes the sizes and the indexes for buffers to be recieved
  SUBROUTINE compute_recv_buffer_sizes()
    
    TYPE(t_grid_comm_pattern), POINTER :: grid_comm_pattern
    
    INTEGER :: comm_var, vertical_layers, np, bfid, total_recv_size, buffer_start
    INTEGER :: no_of_variables
    
    CHARACTER(*), PARAMETER :: method_name = "compute_recv_buffer_sizes"

    ! clear buffer sizes
    DO bfid = 1, active_recv_buffers
      recv_procs_buffer(bfid)%buffer_size = 0
    ENDDO

    ! go through all variables and compute the requested sizes
    DO comm_var = 1, max_active_comm_variables
      
      IF ( comm_variable(comm_var)%request /= communicate ) CYCLE

      grid_comm_pattern => comm_variable(comm_var)%grid_comm_pattern
      vertical_layers = comm_variable(comm_var)%vertical_layers
      no_of_variables = comm_variable(comm_var)%no_of_variables

      ! go through the requested recieve pids
      DO np = 1, grid_comm_pattern%no_of_recv_procs ! loop over PEs from where to recieve the data
        
        total_recv_size = &
          & ((grid_comm_pattern%recv(np)%no_of_points * vertical_layers) + checksum_size) &
          & * no_of_variables
        bfid = grid_comm_pattern%recv(np)%buffer_index
        recv_procs_buffer(bfid)%buffer_size = recv_procs_buffer(bfid)%buffer_size &
          &  + total_recv_size
        
      ENDDO
      
    END DO  
                  
    ! now compute the buffers' start, current and end indexes
    buffer_start = 1
    DO bfid = 1, active_recv_buffers
      recv_procs_buffer(bfid)%start_index = buffer_start
      buffer_start = buffer_start + recv_procs_buffer(bfid)%buffer_size
      recv_procs_buffer(bfid)%end_index = buffer_start - 1
      ! set the current_index to start
      recv_procs_buffer(bfid)%current_index = recv_procs_buffer(bfid)%start_index
    ENDDO
    
    ! check if we went over the max_send_recv_buffer_size
    IF (buffer_start >= max_send_recv_buffer_size) &
      CALL finish(method_name, "buffer_start >= max_send_recv_buffer_size")
 
  END SUBROUTINE compute_recv_buffer_sizes
  !-----------------------------------------------------------------------
  
  !-----------------------------------------------------------------------
!   SUBROUTINE icon_comm_show(var,  comm_pattern)
!     TYPE(t_comm_pattern), INTENT(IN) :: comm_pattern
!     REAL(wp), POINTER, INTENT(IN)   :: var(:,:,:)
! 
!     
!     INTEGER :: dim_2, np, current_buffer_index
!     INTEGER :: i, k, endidx
!     
!    IF(this_is_mpi_sequential) RETURN
!   
!     dim_2=SIZE(var,2)
!     dim_2=1
!     
!     DO np = 1, comm_pattern%np_send
!      
!        endidx = comm_pattern%send_startidx(np) + &
!          & comm_pattern%send_count(np) - 1
!                 
!         DO i = comm_pattern%send_startidx(np), endidx
!           DO k = 1, dim_2
!              
!             write(file_id(my_mpi_work_id),*) " is to send to ", comm_pattern%pelist_send(np), ":", &
!               i,k, var( comm_pattern%send_src_idx(i), k, &
!                      comm_pattern%send_src_blk(i) )
!           ENDDO
!         ENDDO
!         
!     ENDDO
!     
!     DO np = 1, comm_pattern%np_recv
!            
!         endidx = comm_pattern%recv_startidx(np) + &
!           & comm_pattern%recv_count(np) - 1
!                  
!         DO i = comm_pattern%recv_startidx(np), endidx
!           DO k = 1, dim_2
!              
!             write(file_id(my_mpi_work_id),*) " got from ", comm_pattern%pelist_recv(np), ":", &
!               i,k, var( comm_pattern%recv_dst_idx(i), k, &
!                      comm_pattern%recv_dst_blk(i) )
!           ENDDO
!         ENDDO
!         
!     ENDDO
! 
!     DO i = 1, comm_pattern%n_send
!         DO k = 1, dim_2    
!           write(file_id(my_mpi_work_id),*) " send out :", i, k, &
!              var(comm_pattern%send_src_idx(i),k,comm_pattern%send_src_blk(i))
!        ENDDO
!     ENDDO
!     
!     DO i = 1, comm_pattern%n_pnts
!         DO k = 1, dim_2    
!           write(file_id(my_mpi_work_id),*) " received :", i, k, &
!            var(comm_pattern%recv_dst_idx(i),k,comm_pattern%recv_dst_blk(i))
!        ENDDO
!    ENDDO
!     
!       
!   END SUBROUTINE icon_comm_show
  !-----------------------------------------------------------------------

END MODULE mo_icon_comm_lib
