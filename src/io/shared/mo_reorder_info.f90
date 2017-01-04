MODULE mo_reorder_info
  USE mo_kind, ONLY: i4, i8, dp, sp
  USE mo_mpi, ONLY: p_bcast, p_comm_remote_size, p_allgather, p_allgatherv, &
       p_comm_size, p_int_i8
  USE mo_exception, ONLY: finish, message_text
  USE mo_communication,             ONLY: idx_no, blk_no
#ifndef NOMPI
  USE mpi
#endif
  USE mo_util_sort, ONLY: quicksort
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
    !> Only set on compute PEs, set to 0 on IO PEs
    INTEGER, ALLOCATABLE       :: own_idx(:), own_blk(:)
    !> dest idx and blk for own points, only set on sequential/test PEs
    INTEGER, ALLOCATABLE       :: pe_own(:)
    !> n_own, gathered for all compute PEs (set on all PEs)
    INTEGER, ALLOCATABLE       :: pe_off(:)
    !> this ranks reordering indices, i.e. place each packed element goes into
    !! in final array
    INTEGER, ALLOCATABLE       :: reorder_index_own(:)
    !> offset of contributions of all ranks concatenated
    !! (should ONLY be set on I/O servers)
    INTEGER, ALLOCATABLE       :: reorder_index(:)
    ! Index how to reorder the contributions of all compute PEs
    ! into the global array (set on all PEs)

    ! grid information: geographical locations of cells, edges, and
    ! vertices which is first collected on working PE 0 - from where
    ! it will be broadcasted to the pure I/O PEs.
  CONTAINS
    PROCEDURE :: finalize => t_reorder_info_finalize
  END TYPE t_reorder_info
  INTERFACE ri_cpy_part2whole
    MODULE PROCEDURE ri_part2whole_1d_dp_dp, ri_part2whole_1d_sp_sp, &
         ri_part2whole_1d_sp_dp, ri_part2whole_1d_dp_sp, ri_part2whole_1d_i4_i4
    MODULE PROCEDURE ri_part2whole_2d_dp_dp, ri_part2whole_2d_sp_sp, &
         ri_part2whole_2d_sp_dp, ri_part2whole_2d_dp_sp, ri_part2whole_2d_i4_i4
  END INTERFACE ri_cpy_part2whole
  INTERFACE ri_cpy_blk2part
    MODULE PROCEDURE ri_blk2part_2d_dp_dp, ri_blk2part_2d_sp_sp, &
         ri_blk2part_2d_sp_dp, ri_blk2part_2d_i4_i4, ri_blk2part_2d_i4_dp
  END INTERFACE ri_cpy_blk2part
  PUBLIC :: t_reorder_info
  PUBLIC :: transfer_reorder_info
  PUBLIC :: mask2reorder_info
  PUBLIC :: ri_cpy_part2whole
  PUBLIC :: ri_cpy_blk2part
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
       retained_occupation_mask)
    TYPE(t_reorder_info), INTENT(inout) :: ri
    LOGICAL, INTENT(in) :: mask(:)
    INTEGER, INTENT(in) :: n_points_g, glb_index(:), group_comm
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
    CONTIGUOUS :: glb_index, mask
#endif
    INTEGER(i8), ALLOCATABLE, OPTIONAL, INTENT(out) :: &
         retained_occupation_mask(:)

    INTEGER :: n_points, i, il, n, group_comm_size
    INTEGER :: nocc, ierror
    INTEGER(i8) :: pos, apos, bpos
    INTEGER, ALLOCATABLE :: glbidx_own(:), permutation(:), occ_pfxsum(:), &
         buf(:)
    INTEGER(i8), ALLOCATABLE :: occupation_mask(:)
    INTEGER(i8) :: occ_temp, occ_accum
    INTEGER(i8), PARAMETER :: nbits_i8 = BIT_SIZE(occ_temp)
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::set_reorder_data"

    n_points = SIZE(mask)
    n = COUNT(mask)
    ! Get number of owned cells/edges/verts (without halos, physical patch only)
    ri%n_own = n

    ! Set index arrays to own cells/edges/verts
    ! Global index of my own points
    ALLOCATE(ri%own_idx(n), ri%own_blk(n), ri%reorder_index_own(ri%n_own))
    group_comm_size = p_comm_size(group_comm)
    ALLOCATE(ri%pe_own(0:group_comm_size-1), ri%pe_off(0:group_comm_size-1))
    ALLOCATE(glbidx_own(n), permutation(n), buf(n))
    n = 0
    DO i = 1, n_points
      IF (mask(i)) THEN
        n = n+1
        ri%own_idx(n) = idx_no(i)
        ri%own_blk(n) = blk_no(i)
        glbidx_own(n) = glb_index(i)
      ENDIF
    ENDDO

    ! sort used global indices so that target-side re-ordering has lower overhead
    DO i = 1, n
      permutation(i) = i
    END DO
    CALL quicksort(glbidx_own, permutation)
    CALL permute(ri%own_idx, permutation, buf)
    CALL permute(ri%own_blk, permutation, buf)
    DEALLOCATE(buf, permutation)
    ! Gather the number of own points for every PE into ri%pe_own
    CALL p_allgather(ri%n_own, ri%pe_own, group_comm)

    ! Get offset within result array
    il = 0
    DO i = 0, group_comm_size-1
      ri%pe_off(i) = il
      il = il + ri%pe_own(i)
    ENDDO

    ! Get global number of points for current (physical!) patch
    ri%n_glb = il

    ! Compute the global ordering of the local data when it is gathered
    ! on PE 0 or I/O servers exactly as it is retrieved later during I/O
    ! 1. set bit for every global index used
    ! todo: occupation_mask is still a global size array but only uses one bit
    ! per grid cell, i.e. decomposition makes sense but takes additional coding
    ! effort
    nocc = (INT(n_points_g, i8) + nbits_i8 - 1_i8) / nbits_i8
    ALLOCATE(occupation_mask(0:INT(nocc-1_i8)))
    occupation_mask = 0_i8
    DO i = 1, ri%n_own
      pos = INT(glbidx_own(i), i8)
      apos = (pos - 1_i8)/nbits_i8
      bpos = MOD(pos - 1_i8, nbits_i8)
      occupation_mask(apos) = IBSET(occupation_mask(apos), bpos)
    END DO
#ifndef NOMPI
    CALL mpi_allreduce(mpi_in_place, occupation_mask, INT(nocc), p_int_i8, &
      &                mpi_bor, group_comm, ierror)
    IF (ierror /= MPI_SUCCESS) CALL finish(routine, 'mpi_allreduce failed')
#endif

    ! now compute number of bits set in preceding entries of occupation_mask.
    ! After the following loop, occ_pfxsum(i) equals the number of set bits
    ! (i.e. number of global indices used in output) contained in
    ! occupation_mask(0:i-1)
    ALLOCATE(occ_pfxsum(0:INT(nocc-1)))
    occ_accum = 0
    DO i = 0, INT(nocc-1)
      occ_pfxsum(i) = INT(occ_accum)
      occ_temp = POPCNT(occupation_mask(i))
      occ_accum = occ_accum + occ_temp
    END DO
    IF (occ_accum /= ri%n_glb) THEN
      WRITE (message_text, '(2(a,i0))') 'Bit-counting failed: n=', occ_accum, &
           & ' /= ri%n_glb=', ri%n_glb
      CALL finish(routine,TRIM(message_text))
    ENDIF
    ! given the above two arrays, one can now compute for each global index its
    ! position in the output array in O(1)
    DO i = 1, ri%n_own
      pos = INT(glbidx_own(i), i8)
      apos = (pos - 1_i8)/nbits_i8
      bpos = MOD(pos - 1_i8, nbits_i8)
      occ_accum = IAND(ISHFT(1_i8, bpos) - 1_i8, occupation_mask(apos))
      ri%reorder_index_own(i) = occ_pfxsum(apos) + POPCNT(occ_accum) + 1
    END DO

    DEALLOCATE(glbidx_own)
    IF (PRESENT(retained_occupation_mask)) THEN
      CALL MOVE_ALLOC(occupation_mask, retained_occupation_mask)
    ELSE
      DEALLOCATE(occupation_mask)
    END IF
  END SUBROUTINE mask2reorder_info

  SUBROUTINE permute(a, permutation, a_perm)
    INTEGER, ALLOCATABLE, INTENT(inout) :: a(:)
    INTEGER, ALLOCATABLE, INTENT(inout) :: a_perm(:)
    INTEGER, INTENT(in) :: permutation(:)
    CHARACTER(len=*), PARAMETER :: routine=modname//"::permute"
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
    CONTIGUOUS :: permutation
#endif
    INTEGER, ALLOCATABLE :: temp(:)
    INTEGER :: i, n
    n = SIZE(a)
    IF (SIZE(permutation) /= n) THEN
      CALL finish(routine, "invalid permutation")
    END IF
    IF (ALLOCATED(a_perm)) THEN
      IF (SIZE(a_perm) /= n) THEN
        DEALLOCATE(a_perm)
        ALLOCATE(a_perm(n))
      END IF
    ELSE
      ALLOCATE(a_perm(n))
    END IF
    DO i = 1, n
      a_perm(i) = a(permutation(i))
    END DO
    CALL MOVE_ALLOC(a, temp)
    CALL MOVE_ALLOC(a_perm, a)
    CALL MOVE_ALLOC(temp, a_perm)
  END SUBROUTINE permute
  !------------------------------------------------------------------------------------------------
  !> Transfers reorder_info from clients to servers for IO/async latbc
  !
  SUBROUTINE transfer_reorder_info(ri, is_server, bcast_root, &
       client_server_intercomm)
    TYPE(t_reorder_info), INTENT(INOUT) :: ri ! Result: reorder info
    LOGICAL, INTENT(in) :: is_server
    INTEGER, INTENT(in) :: bcast_root, client_server_intercomm

    INTEGER :: other_group_size, dummy(1)
    INTEGER, ALLOCATABLE :: zcounts(:)
    ! Transfer the global number of points, this is not yet known on IO PEs
    CALL p_bcast(ri%n_glb,  bcast_root, client_server_intercomm)

    other_group_size = p_comm_remote_size(client_server_intercomm)
    IF (is_server) THEN

      ! On IO PEs: n_own = 0, own_idx and own_blk are not allocated
      ri%n_own = 0
      ! pe_own/pe_off must be allocated for client_group_size, not for p_n_work
      ALLOCATE(ri%pe_own(0:other_group_size-1))
      ALLOCATE(ri%pe_off(0:other_group_size-1))

      ALLOCATE(ri%reorder_index(ri%n_glb))
    ENDIF

    CALL p_bcast(ri%pe_own, bcast_root, client_server_intercomm)
    CALL p_bcast(ri%pe_off, bcast_root, client_server_intercomm)

    IF (is_server) THEN
      CALL p_allgatherv(dummy(1:0), ri%reorder_index, &
        &               ri%pe_own, ri%pe_off, client_server_intercomm)
    ELSE
      ALLOCATE(zcounts(0:other_group_size-1))
      zcounts = 0
      CALL p_allgatherv(ri%reorder_index_own, dummy(1:0), &
        &               zcounts, zcounts, client_server_intercomm)
    END IF

  END SUBROUTINE transfer_reorder_info

  SUBROUTINE ri_part2whole_1d_dp_dp(ri, part_idx, part_data, whole_data)
    TYPE(t_reorder_info), INTENT(IN) :: ri
    INTEGER, INTENT(in) :: part_idx
    REAL(dp), INTENT(in) :: part_data(:)
    REAL(dp), INTENT(inout) :: whole_data(:)
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
    CONTIGUOUS :: part_data, whole_data
#endif
    INTEGER :: i, n, ofs

    ofs = ri%pe_off(part_idx)
    n = ri%pe_own(part_idx)
    DO i = 1, n
      whole_data(ri%reorder_index(ofs + i)) = part_data(i)
    END DO
  END SUBROUTINE ri_part2whole_1d_dp_dp

  SUBROUTINE ri_part2whole_1d_sp_sp(ri, part_idx, part_data, whole_data)
    TYPE(t_reorder_info), INTENT(IN) :: ri
    INTEGER, INTENT(in) :: part_idx
    REAL(sp), INTENT(in) :: part_data(:)
    REAL(sp), INTENT(inout) :: whole_data(:)
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
    CONTIGUOUS :: part_data, whole_data
#endif
    INTEGER :: i, n, ofs

    ofs = ri%pe_off(part_idx)
    n = ri%pe_own(part_idx)
    DO i = 1, n
      whole_data(ri%reorder_index(ofs + i)) = part_data(i)
    END DO
  END SUBROUTINE ri_part2whole_1d_sp_sp

  SUBROUTINE ri_part2whole_1d_sp_dp(ri, part_idx, part_data, whole_data)
    TYPE(t_reorder_info), INTENT(IN) :: ri
    INTEGER, INTENT(in) :: part_idx
    REAL(sp), INTENT(in) :: part_data(:)
    REAL(dp), INTENT(inout) :: whole_data(:)
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
    CONTIGUOUS :: part_data, whole_data
#endif
    INTEGER :: i, n, ofs

    ofs = ri%pe_off(part_idx)
    n = ri%pe_own(part_idx)
    DO i = 1, n
      whole_data(ri%reorder_index(ofs + i)) = part_data(i)
    END DO
  END SUBROUTINE ri_part2whole_1d_sp_dp

  SUBROUTINE ri_part2whole_1d_dp_sp(ri, part_idx, part_data, whole_data)
    TYPE(t_reorder_info), INTENT(IN) :: ri
    INTEGER, INTENT(in) :: part_idx
    REAL(dp), INTENT(in) :: part_data(:)
    REAL(sp), INTENT(inout) :: whole_data(:)
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
    CONTIGUOUS :: part_data, whole_data
#endif
    INTEGER :: i, n, ofs

    ofs = ri%pe_off(part_idx)
    n = ri%pe_own(part_idx)
    DO i = 1, n
      whole_data(ri%reorder_index(ofs + i)) = REAL(part_data(i), sp)
    END DO
  END SUBROUTINE ri_part2whole_1d_dp_sp

  SUBROUTINE ri_part2whole_1d_i4_i4(ri, part_idx, part_data, whole_data)
    TYPE(t_reorder_info), INTENT(IN) :: ri
    INTEGER, INTENT(in) :: part_idx
    INTEGER(i4), INTENT(in) :: part_data(:)
    INTEGER(i4), INTENT(inout) :: whole_data(:)
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
    CONTIGUOUS :: part_data, whole_data
#endif
    INTEGER :: i, n, ofs

    ofs = ri%pe_off(part_idx)
    n = ri%pe_own(part_idx)
    DO i = 1, n
      whole_data(ri%reorder_index(ofs + i)) = part_data(i)
    END DO
  END SUBROUTINE ri_part2whole_1d_i4_i4

  SUBROUTINE ri_part2whole_2d_dp_dp(ri, part_idx, part_data, whole_data)
    TYPE(t_reorder_info), INTENT(IN) :: ri
    INTEGER, INTENT(in) :: part_idx
    REAL(dp), INTENT(in) :: part_data(:,:)
    REAL(dp), INTENT(inout) :: whole_data(:,:)
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
    CONTIGUOUS :: part_data, whole_data
#endif
    INTEGER :: i, k, n, nlev, ofs

    ofs = ri%pe_off(part_idx)
    n = ri%pe_own(part_idx)
    nlev = SIZE(part_data, 2)
    DO k = 1, nlev
      DO i = 1, n
        whole_data(ri%reorder_index(ofs + i), k) = part_data(i, k)
      END DO
    END DO
  END SUBROUTINE ri_part2whole_2d_dp_dp

  SUBROUTINE ri_part2whole_2d_sp_sp(ri, part_idx, part_data, whole_data)
    TYPE(t_reorder_info), INTENT(IN) :: ri
    INTEGER, INTENT(in) :: part_idx
    REAL(sp), INTENT(in) :: part_data(:,:)
    REAL(sp), INTENT(inout) :: whole_data(:,:)
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
    CONTIGUOUS :: part_data, whole_data
#endif
    INTEGER :: i, k, n, nlev, ofs

    ofs = ri%pe_off(part_idx)
    n = ri%pe_own(part_idx)
    nlev = SIZE(part_data, 2)
    DO k = 1, nlev
      DO i = 1, n
        whole_data(ri%reorder_index(ofs + i), k) = part_data(i, k)
      END DO
    END DO
  END SUBROUTINE ri_part2whole_2d_sp_sp

  SUBROUTINE ri_part2whole_2d_sp_dp(ri, part_idx, part_data, whole_data)
    TYPE(t_reorder_info), INTENT(IN) :: ri
    INTEGER, INTENT(in) :: part_idx
    REAL(sp), INTENT(in) :: part_data(:,:)
    REAL(dp), INTENT(inout) :: whole_data(:,:)
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
    CONTIGUOUS :: part_data, whole_data
#endif
    INTEGER :: i, k, n, nlev, ofs

    ofs = ri%pe_off(part_idx)
    n = ri%pe_own(part_idx)
    nlev = SIZE(part_data, 2)
    DO k = 1, nlev
      DO i = 1, n
        whole_data(ri%reorder_index(ofs + i), k) = part_data(i, k)
      END DO
    END DO
  END SUBROUTINE ri_part2whole_2d_sp_dp

  SUBROUTINE ri_part2whole_2d_dp_sp(ri, part_idx, part_data, whole_data)
    TYPE(t_reorder_info), INTENT(IN) :: ri
    INTEGER, INTENT(in) :: part_idx
    REAL(dp), INTENT(in) :: part_data(:,:)
    REAL(sp), INTENT(inout) :: whole_data(:,:)
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
    CONTIGUOUS :: part_data, whole_data
#endif
    INTEGER :: i, k, n, nlev, ofs

    ofs = ri%pe_off(part_idx)
    n = ri%pe_own(part_idx)
    nlev = SIZE(part_data, 2)
    DO k = 1, nlev
      DO i = 1, n
        whole_data(ri%reorder_index(ofs + i), k) = REAL(part_data(i, k), sp)
      END DO
    END DO
  END SUBROUTINE ri_part2whole_2d_dp_sp

  SUBROUTINE ri_part2whole_2d_i4_i4(ri, part_idx, part_data, whole_data)
    TYPE(t_reorder_info), INTENT(IN) :: ri
    INTEGER, INTENT(in) :: part_idx
    INTEGER(i4), INTENT(in) :: part_data(:,:)
    INTEGER(i4), INTENT(inout) :: whole_data(:,:)
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
    CONTIGUOUS :: part_data, whole_data
#endif
    INTEGER :: i, k, n, nlev, ofs

    ofs = ri%pe_off(part_idx)
    n = ri%pe_own(part_idx)
    nlev = SIZE(part_data, 2)
    DO k = 1, nlev
      DO i = 1, n
        whole_data(ri%reorder_index(ofs + i), k) = part_data(i, k)
      END DO
    END DO
  END SUBROUTINE ri_part2whole_2d_i4_i4

  ! in the following routines, part_data has the POINTER attribute because
  ! ifort is being a dick about copying in the array otherwise
  SUBROUTINE ri_blk2part_2d_dp_dp(ri, blk_data, part_data, offset)
    TYPE(t_reorder_info), INTENT(IN) :: ri
    REAL(dp), INTENT(in) :: blk_data(:,:)
    REAL(dp), INTENT(inout) :: part_data(:)
    INTEGER, INTENT(inout) :: offset
    INTEGER :: i, n

    n = ri%n_own
    DO i = 1, n
      part_data(offset + i) = blk_data(ri%own_idx(i), ri%own_blk(i))
    END DO
    offset = offset + n
  END SUBROUTINE ri_blk2part_2d_dp_dp

  SUBROUTINE ri_blk2part_2d_sp_sp(ri, blk_data, part_data, offset)
    TYPE(t_reorder_info), INTENT(IN) :: ri
    REAL(sp), INTENT(in) :: blk_data(:,:)
    REAL(sp), INTENT(inout) :: part_data(:)
    INTEGER, INTENT(inout) :: offset
    INTEGER :: i, n

    n = ri%n_own
    DO i = 1, n
      part_data(offset + i) = blk_data(ri%own_idx(i), ri%own_blk(i))
    END DO
    offset = offset + n
  END SUBROUTINE ri_blk2part_2d_sp_sp

  SUBROUTINE ri_blk2part_2d_sp_dp(ri, blk_data, part_data, offset)
    TYPE(t_reorder_info), INTENT(IN) :: ri
    REAL(sp), INTENT(in) :: blk_data(:,:)
    REAL(dp), INTENT(inout) :: part_data(:)
    INTEGER, INTENT(inout) :: offset
    INTEGER :: i, n

    n = ri%n_own
    DO i = 1, n
      part_data(offset + i) = REAL(blk_data(ri%own_idx(i), ri%own_blk(i)), dp)
    END DO
    offset = offset + n
  END SUBROUTINE ri_blk2part_2d_sp_dp

  SUBROUTINE ri_blk2part_2d_i4_i4(ri, blk_data, part_data, offset)
    TYPE(t_reorder_info), INTENT(IN) :: ri
    INTEGER(i4), INTENT(in) :: blk_data(:,:)
    INTEGER(i4), INTENT(inout) :: part_data(:)
    INTEGER, INTENT(inout) :: offset
    INTEGER :: i, n

    n = ri%n_own
    DO i = 1, n
      part_data(offset + i) = blk_data(ri%own_idx(i), ri%own_blk(i))
    END DO
    offset = offset + n
  END SUBROUTINE ri_blk2part_2d_i4_i4

  SUBROUTINE ri_blk2part_2d_i4_dp(ri, blk_data, part_data, offset)
    TYPE(t_reorder_info), INTENT(IN) :: ri
    INTEGER(i4), INTENT(in) :: blk_data(:,:)
    REAL(dp), INTENT(inout) :: part_data(:)
    INTEGER, INTENT(inout) :: offset
    INTEGER :: i, n

    n = ri%n_own
    DO i = 1, n
      part_data(offset + i) = REAL(blk_data(ri%own_idx(i), ri%own_blk(i)), dp)
    END DO
    offset = offset + n
  END SUBROUTINE ri_blk2part_2d_i4_dp

END MODULE mo_reorder_info
