MODULE mo_reorder_info
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
END MODULE mo_reorder_info
