MODULE mo_reorder_info
  USE mo_mpi, ONLY: p_bcast, p_comm_remote_size
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
  END TYPE t_reorder_info
  PUBLIC :: t_reorder_info
  PUBLIC :: transfer_reorder_info
CONTAINS
  !------------------------------------------------------------------------------------------------
  !> Transfers reorder_info from clients to servers for IO/async latbc
  !
  SUBROUTINE transfer_reorder_info(p_ri, is_server, bcast_root, &
       client_server_intercomm)
    TYPE(t_reorder_info), INTENT(INOUT) :: p_ri ! Result: reorder info
    LOGICAL, INTENT(in) :: is_server
    INTEGER, INTENT(in) :: bcast_root, client_server_intercomm

    INTEGER :: client_group_size
    ! Transfer the global number of points, this is not yet known on IO PEs
    CALL p_bcast(p_ri%n_glb,  bcast_root, client_server_intercomm)

    IF (is_server) THEN

      ! On IO PEs: n_own = 0, own_idx and own_blk are not allocated
      p_ri%n_own = 0
      client_group_size = p_comm_remote_size(client_server_intercomm)
      ! pe_own/pe_off must be allocated for client_group_size, not for p_n_work
      ALLOCATE(p_ri%pe_own(0:client_group_size-1))
      ALLOCATE(p_ri%pe_off(0:client_group_size-1))

      ALLOCATE(p_ri%reorder_index(p_ri%n_glb))
    ENDIF

    CALL p_bcast(p_ri%pe_own, bcast_root, client_server_intercomm)
    CALL p_bcast(p_ri%pe_off, bcast_root, client_server_intercomm)

    CALL p_bcast(p_ri%reorder_index, bcast_root, client_server_intercomm)

  END SUBROUTINE transfer_reorder_info

END MODULE mo_reorder_info
