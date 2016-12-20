MODULE mo_reorder_info
  USE mo_mpi, ONLY: p_bcast, p_comm_remote_size, p_allgather, p_allgatherv, &
       p_comm_size
  USE mo_exception, ONLY: finish, message_text
  USE mo_communication,             ONLY: idx_no, blk_no
  IMPLICIT NONE
  PRIVATE
  !------------------------------------------------------------------------------------------------
  ! TYPE t_reorder_info describes how local cells/edges/verts
  ! have to be reordered to get the global array.
  ! Below, "points" refers to either cells, edges or verts.
  !
  ! TODO[FP] Note that the "reorder_info" contains fields of *global*
  !          size (reorder_index). On the compute PEs these fields
  !          could be deallocated after the call to
  !          "transfer_reorder_info" in the setup phase!

  TYPE t_reorder_info
    INTEGER                    :: n_glb  ! Global number of points per physical patch
    INTEGER                    :: n_own  ! Number of own points (without halo, only belonging to phyiscal patch)
    ! Only set on compute PEs, set to 0 on IO PEs
    INTEGER, ALLOCATABLE       :: own_idx(:), own_blk(:)
    ! dest idx and blk for own points, only set on sequential/test PEs
    INTEGER, ALLOCATABLE       :: pe_own(:)
    ! n_own, gathered for all compute PEs (set on all PEs)
    INTEGER, ALLOCATABLE       :: pe_off(:)
    ! offset of contributions of PEs (set on all PEs)
    INTEGER, ALLOCATABLE       :: reorder_index(:)
    ! Index how to reorder the contributions of all compute PEs
    ! into the global array (set on all PEs)

    ! grid information: geographical locations of cells, edges, and
    ! vertices which is first collected on working PE 0 - from where
    ! it will be broadcasted to the pure I/O PEs.
  CONTAINS
    PROCEDURE :: finalize => t_reorder_info_finalize
  END TYPE t_reorder_info
  PUBLIC :: t_reorder_info
  PUBLIC :: transfer_reorder_info
  PUBLIC :: mask2reorder_info
  CHARACTER(len=*), PARAMETER :: modname = 'mo_reorder_info'
CONTAINS

  SUBROUTINE t_reorder_info_finalize(reorder_data)
    CLASS(t_reorder_info), INTENT(INOUT) :: reorder_data

    !CALL message("", 't_reorder_data_finalize')

    IF (ALLOCATED(reorder_data%reorder_index))  DEALLOCATE(reorder_data%reorder_index)
    IF (ALLOCATED(reorder_data%own_idx))        DEALLOCATE(reorder_data%own_idx)
    IF (ALLOCATED(reorder_data%own_blk))        DEALLOCATE(reorder_data%own_blk)
    IF (ALLOCATED(reorder_data%pe_own))         DEALLOCATE(reorder_data%pe_own)
    IF (ALLOCATED(reorder_data%pe_off))         DEALLOCATE(reorder_data%pe_off)
  END SUBROUTINE t_reorder_info_finalize

  SUBROUTINE mask2reorder_info(ri, mask, n_points_g, glb_index, group_comm, &
       retained_reorder_index_log_dom)
    TYPE(t_reorder_info), INTENT(inout) :: ri
    LOGICAL, INTENT(in) :: mask(:)
    INTEGER, INTENT(in) :: n_points_g, glb_index(:), group_comm
    INTEGER, ALLOCATABLE, OPTIONAL, INTENT(out) :: &
         retained_reorder_index_log_dom(:)

    INTEGER :: n_points, i, il, n, group_comm_size
    INTEGER, ALLOCATABLE :: glbidx_own(:), glbidx_glb(:), &
         reorder_index_log_dom(:)
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::set_reorder_data"

    n_points = SIZE(mask)
    n = COUNT(mask)
    ! Get number of owned cells/edges/verts (without halos, physical patch only)
    ri%n_own = n

    ! Set index arrays to own cells/edges/verts
    ! Global index of my own points
    ALLOCATE(ri%own_idx(n), ri%own_blk(n), glbidx_own(n))

    n = 0
    DO i = 1, n_points
      IF (mask(i)) THEN
        n = n+1
        ri%own_idx(n) = idx_no(i)
        ri%own_blk(n) = blk_no(i)
        glbidx_own(n) = glb_index(i)
      ENDIF
    ENDDO

    ! Gather the number of own points for every PE into ri%pe_own
    group_comm_size = p_comm_size(group_comm)
    ALLOCATE(ri%pe_own(0:group_comm_size-1), ri%pe_off(0:group_comm_size-1))
    CALL p_allgather(ri%n_own, ri%pe_own, group_comm)

    ! Get offset within result array
    il = 0
    DO i = 0, group_comm_size-1
      ri%pe_off(i) = il
      il = il + ri%pe_own(i)
    ENDDO

    ! Get global number of points for current (physical!) patch
    ri%n_glb = il

    ! Get the global index numbers of the data when it is gathered on PE 0
    ! exactly in the same order as it is retrieved later during I/O

    ALLOCATE(glbidx_glb(ri%n_glb))
    CALL p_allgatherv(glbidx_own(1:ri%n_own), glbidx_glb, &
      &               ri%pe_own, ri%pe_off, group_comm)

    ! Get reorder_index
    ALLOCATE(ri%reorder_index(ri%n_glb), &
      ! spans the complete logical domain
      &      reorder_index_log_dom(n_points_g))
    reorder_index_log_dom(:) = 0

    DO i = 1, ri%n_glb
      ! reorder_index_log_dom stores where a global point in logical domain comes from.
      ! It is nonzero only at the physical patch locations
      reorder_index_log_dom(glbidx_glb(i)) = i
    ENDDO

    ! Gather the reorder index for the physical domain
    n = 0
    DO i = 1, n_points_g
      IF(reorder_index_log_dom(i)>0) THEN
        n = n+1
        ri%reorder_index(reorder_index_log_dom(i)) = n
      ENDIF
    ENDDO

    ! Safety check
    IF(n/=ri%n_glb) THEN
      WRITE (message_text, '(2(a,i0))') 'Reordering failed: n=',n, &
           & ' /= ri%n_glb=',ri%n_glb
      CALL finish(routine,TRIM(message_text))
    ENDIF

    DEALLOCATE(glbidx_own, glbidx_glb)
    IF (PRESENT(retained_reorder_index_log_dom)) THEN
      CALL MOVE_ALLOC(reorder_index_log_dom, retained_reorder_index_log_dom)
    ELSE
      DEALLOCATE(reorder_index_log_dom)
    END IF
  END SUBROUTINE mask2reorder_info

  !------------------------------------------------------------------------------------------------
  !> Transfers reorder_info from clients to servers for IO/async latbc
  !
  SUBROUTINE transfer_reorder_info(ri, is_server, bcast_root, &
       client_server_intercomm)
    TYPE(t_reorder_info), INTENT(INOUT) :: ri ! Result: reorder info
    LOGICAL, INTENT(in) :: is_server
    INTEGER, INTENT(in) :: bcast_root, client_server_intercomm

    INTEGER :: client_group_size
    ! Transfer the global number of points, this is not yet known on IO PEs
    CALL p_bcast(ri%n_glb,  bcast_root, client_server_intercomm)

    IF (is_server) THEN

      ! On IO PEs: n_own = 0, own_idx and own_blk are not allocated
      ri%n_own = 0
      client_group_size = p_comm_remote_size(client_server_intercomm)
      ! pe_own/pe_off must be allocated for client_group_size, not for p_n_work
      ALLOCATE(ri%pe_own(0:client_group_size-1))
      ALLOCATE(ri%pe_off(0:client_group_size-1))

      ALLOCATE(ri%reorder_index(ri%n_glb))
    ENDIF

    CALL p_bcast(ri%pe_own, bcast_root, client_server_intercomm)
    CALL p_bcast(ri%pe_off, bcast_root, client_server_intercomm)

    CALL p_bcast(ri%reorder_index, bcast_root, client_server_intercomm)

  END SUBROUTINE transfer_reorder_info

END MODULE mo_reorder_info
