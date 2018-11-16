!#define HAVE_YAXT
!>
!!               This module provides the yaxt based communication routines.
!!
!!               This module provides the yaxt based communication routines
!! for parallel runs
!!
!! @par Revision History
!! Initial version by Moritz Hanke, April 2016
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!

!----------------------------
#include "icon_definitions.inc"
#include "crayftn_ptr_fail.inc"
!----------------------------
MODULE mo_communication_yaxt
#ifdef HAVE_YAXT

USE mpi
USE mo_kind,               ONLY: dp, sp
USE mo_exception,          ONLY: finish, message
USE mo_mpi,                ONLY: p_comm_work, p_barrier, &
  &                              p_real_dp, p_real_sp, p_int, p_address_kind, &
  &                              p_alltoall, p_alltoallv, p_bcast, &
  &                              p_comm_size, p_comm_rank
USE mo_parallel_config, ONLY: nproma, itype_exch_barrier
USE mo_timer,           ONLY: timer_start, timer_stop, activate_sync_timers, &
  &                           timer_exch_data, timer_barrier
USE mo_decomposition_tools, ONLY: t_glb2loc_index_lookup, get_local_index
USE mo_parallel_config, ONLY: blk_no, idx_no, idx_1d
USE yaxt, ONLY: xt_initialized, xt_initialize, xt_idxlist, &
  &             xt_int_kind, xt_idxvec_new, &
  &             xt_xmap, xt_xmap_delete, xt_xmap_intersection_new, &
  &             xt_idxlist_delete, xt_redist, xt_redist_p2p_off_new, &
  &             xt_redist_s_exchange1, xt_redist_delete, xt_redist_p2p_new, &
  &             xt_redist_collection_new, xt_redist_s_exchange, &
  &             xt_xmap_get_num_sources, xt_xmap_get_num_destinations, &
  &             xt_xmap_get_source_ranks, xt_com_list, xt_redist_repeat_new, &
  &             xt_mpi_comm_mark_exclusive, xt_int_kind, xt_slice_c_loc
USE mo_fortran_tools,        ONLY: t_ptr_3d, t_ptr_3d_sp
USE iso_c_binding, ONLY: c_int, c_loc, c_ptr
USE mo_communication_types, ONLY: t_comm_pattern, t_p_comm_pattern, &
  & t_comm_pattern_collection, xfer_list


IMPLICIT NONE

PRIVATE

!modules interface-------------------------------------------
!subroutines
PUBLIC :: t_comm_pattern_yaxt
PUBLIC :: t_comm_pattern_collection_yaxt
!
!variables

!--------------------------------------------------------------------------------------------------
!
TYPE t_comm_pattern_redist

  PRIVATE

  INTEGER :: nfields
  INTEGER, ALLOCATABLE :: dst_nlev(:)
  INTEGER, ALLOCATABLE :: src_nlev(:)
  INTEGER, ALLOCATABLE :: nshift(:)
  INTEGER :: mpi_type
  TYPE(xt_redist) :: redist

END TYPE t_comm_pattern_redist

TYPE t_comm_pattern_contiguous_data_type

  PRIVATE

  INTEGER :: base_type
  INTEGER :: count
  INTEGER :: type

END TYPE t_comm_pattern_contiguous_data_type

TYPE, EXTENDS(t_comm_pattern) :: t_comm_pattern_yaxt

  PRIVATE

    ! if anything is changed here, please ensure that the deep copy in
    ! setup_comm_pattern_collection still works correctly
    INTEGER :: src_n_points, dst_n_points, comm = mpi_comm_null
    LOGICAL, ALLOCATABLE :: dst_mask(:)
    TYPE(xt_xmap) :: xmap
    TYPE(t_comm_pattern_redist), ALLOCATABLE :: redists(:)
    TYPE(t_comm_pattern_contiguous_data_type), &
         ALLOCATABLE :: contiguous_data_types(:)
    LOGICAL :: inplace

  CONTAINS

    PROCEDURE :: setup => setup_comm_pattern
    PROCEDURE :: setup2 => setup_comm_pattern2
    PROCEDURE :: delete => delete_comm_pattern
    PROCEDURE :: exchange_data_r3d => exchange_data_r3d
    PROCEDURE :: exchange_data_s3d => exchange_data_s3d
    PROCEDURE :: exchange_data_i3d => exchange_data_i3d
    PROCEDURE :: exchange_data_r2d => exchange_data_r2d
    PROCEDURE :: exchange_data_s2d => exchange_data_s2d
    PROCEDURE :: exchange_data_i2d => exchange_data_i2d
    PROCEDURE :: exchange_data_l2d => exchange_data_l2d
    PROCEDURE :: exchange_data_l3d => exchange_data_l3d
    PROCEDURE :: exchange_data_mult => exchange_data_mult_dp
    PROCEDURE :: exchange_data_mult_mixprec => exchange_data_mult_mixprec
    PROCEDURE :: exchange_data_4de1 => exchange_data_4de1
    PROCEDURE :: get_np_recv => get_np_recv
    PROCEDURE :: get_np_send => get_np_send
    PROCEDURE :: get_pelist_recv => get_pelist_recv

END TYPE t_comm_pattern_yaxt

TYPE t_comm_pattern_coll_redist

  PRIVATE

  INTEGER :: nfields
  INTEGER, ALLOCATABLE :: dst_nlev(:)
  INTEGER, ALLOCATABLE :: src_nlev(:)
  TYPE(xt_redist) :: redist

END TYPE t_comm_pattern_coll_redist

TYPE t_p_comm_pattern_yaxt
  TYPE(t_comm_pattern_yaxt), POINTER :: p
END TYPE t_p_comm_pattern_yaxt

  TYPE, EXTENDS(t_comm_pattern_collection) :: t_comm_pattern_collection_yaxt
    PRIVATE

    TYPE(t_comm_pattern_coll_redist), ALLOCATABLE :: redists(:)
    TYPE(t_p_comm_pattern_yaxt), ALLOCATABLE :: patterns(:)

  CONTAINS

    PROCEDURE :: setup => setup_comm_pattern_collection
    PROCEDURE :: delete => delete_comm_pattern_collection
    PROCEDURE :: exchange_data_grf => exchange_data_grf

END TYPE t_comm_pattern_collection_yaxt

#ifndef HAVE_IS_CONTIGUOUS
INTERFACE IS_CONTIGUOUS
  MODULE PROCEDURE MY_IS_CONTIGUOUS_DP_3D
  MODULE PROCEDURE MY_IS_CONTIGUOUS_SP_3D
  MODULE PROCEDURE MY_IS_CONTIGUOUS_DP_4D
  MODULE PROCEDURE MY_IS_CONTIGUOUS_SP_4D
END INTERFACE
#endif

!--------------------------------------------------------------------------------------------------
!

CHARACTER(*), PARAMETER :: modname = "mo_communication_yaxt"

!-------------------------------------------------------------------------

CONTAINS

#ifndef HAVE_IS_CONTIGUOUS
FUNCTION MY_IS_CONTIGUOUS_DP_3D(array)
  REAL(dp), INTENT(IN), TARGET :: array(:,:,:)
  LOGICAL :: MY_IS_CONTIGUOUS_DP_3D

  REAL(dp) :: dummy(2)
  INTEGER(KIND=p_address_kind) :: addr(2), type_size
  INTEGER :: ierror

  CALL MPI_Get_address(dummy(1), addr(1), ierror)
  CALL MPI_Get_address(dummy(2), addr(2), ierror)
  type_size = addr(2) - addr(1)
  CALL MPI_Get_address(array(LBOUND(array,1),LBOUND(array,2),LBOUND(array,3)), &
    &                  addr(1), ierror)
  CALL MPI_Get_address(array(UBOUND(array,1),UBOUND(array,2),UBOUND(array,3)), &
    &                  addr(2), ierror)
  MY_IS_CONTIGUOUS_DP_3D = ((addr(2) - addr(1)) / type_size + 1) == SIZE(array)
END FUNCTION MY_IS_CONTIGUOUS_DP_3D

FUNCTION MY_IS_CONTIGUOUS_SP_3D(array)
  REAL(sp), INTENT(IN), TARGET :: array(:,:,:)
  LOGICAL :: MY_IS_CONTIGUOUS_SP_3D

  REAL(sp) :: dummy(2)
  INTEGER(KIND=p_address_kind) :: addr(2), type_size
  INTEGER :: ierror

  CALL MPI_Get_address(dummy(1), addr(1), ierror)
  CALL MPI_Get_address(dummy(2), addr(2), ierror)
  type_size = addr(2) - addr(1)
  CALL MPI_Get_address(array(LBOUND(array,1),LBOUND(array,2),LBOUND(array,3)), &
    &                  addr(1), ierror)
  CALL MPI_Get_address(array(UBOUND(array,1),UBOUND(array,2),UBOUND(array,3)), &
    &                  addr(2), ierror)
  MY_IS_CONTIGUOUS_SP_3D = ((addr(2) - addr(1)) / type_size + 1) == SIZE(array)
END FUNCTION MY_IS_CONTIGUOUS_SP_3D

FUNCTION MY_IS_CONTIGUOUS_DP_4D(array)
  REAL(dp), INTENT(IN), TARGET :: array(:,:,:,:)
  LOGICAL :: MY_IS_CONTIGUOUS_DP_4D

  REAL(dp) :: dummy(2)
  INTEGER(KIND=p_address_kind) :: addr(2), type_size
  INTEGER :: ierror

  CALL MPI_Get_address(dummy(1), addr(1), ierror)
  CALL MPI_Get_address(dummy(2), addr(2), ierror)
  type_size = addr(2) - addr(1)
  CALL MPI_Get_address(array(LBOUND(array,1),LBOUND(array,2),LBOUND(array,3), &
                             LBOUND(array,4)), addr(1), ierror)
  CALL MPI_Get_address(array(UBOUND(array,1),UBOUND(array,2),UBOUND(array,3), &
                             UBOUND(array,4)), addr(2), ierror)
  MY_IS_CONTIGUOUS_DP_4D = ((addr(2) - addr(1)) / type_size + 1) == SIZE(array)
END FUNCTION MY_IS_CONTIGUOUS_DP_4D

FUNCTION MY_IS_CONTIGUOUS_SP_4D(array)
  REAL(sp), INTENT(IN), TARGET :: array(:,:,:,:)
  LOGICAL :: MY_IS_CONTIGUOUS_SP_4D

  REAL(sp) :: dummy(2)
  INTEGER(KIND=p_address_kind) :: addr(2), type_size
  INTEGER :: ierror

  CALL MPI_Get_address(dummy(1), addr(1), ierror)
  CALL MPI_Get_address(dummy(2), addr(2), ierror)
  type_size = addr(2) - addr(1)
  CALL MPI_Get_address(array(LBOUND(array,1),LBOUND(array,2),LBOUND(array,3), &
                       LBOUND(array,4)), addr(1), ierror)
  CALL MPI_Get_address(array(UBOUND(array,1),UBOUND(array,2),UBOUND(array,3), &
                       UBOUND(array,4)), addr(2), ierror)
  MY_IS_CONTIGUOUS_SP_4D = ((addr(2) - addr(1)) / type_size + 1) == SIZE(array)
END FUNCTION MY_IS_CONTIGUOUS_SP_4D
#endif


!-------------------------------------------------------------------------
!
!

!>
!! Sets up a communication pattern for exchanging data.
!!
!! Note: This setup routine works only for the trivial communication
!!       patterns in sequential runs.
!!
!! dst_n_points     Total number of points in the RECEIVER array,
!!                  not every point is necessarily set during exchange
!!                  (see owner!)
!!
!! dst_owner        Owner PE number of every point in the RECEIVER array,
!!                  if owner(.) == -1, this point will not be set during exchange.
!!                  If owner(.) == p_pe, this point will be exchanged,
!!                  this is necessary if sender and receiver arrays are
!!                  different (e.g. feedback, gather, scatter)
!!
!! dst_global_index Global index of of every point in the RECEIVER array
!!                  There may be more than 1 point with the same global index,
!!                  in this case the point is exchanged only once and
!!                  locally distributed.
!!                  - If this argument is not present, we assume global_index=1,2.3,...
!! inplace          In-place data exchanges are allowed (source and destination
!!                  arrays can be identically for the exchange)
!!                  - if inplace == true, the user guarantees that
!!                    (src_n_points == dst_n_points) and that points, which will
!!                    have to be sent to other processes, are disjunct from the
!!                    points that will have to be received
!!                  - in case the user only provides the receive array to an
!!                    exchange call and not send array, the exchange will be
!!                    faster if inplace == true
!!
!! send_decomp_info domain decomposition information for the SENDER array
!!
!! @par Revision History
!! Initial version by Rainer Johanni, Nov 2009
!! yaxt version by Moritz Hanke, April 2016
!!
SUBROUTINE setup_comm_pattern(p_pat, dst_n_points, dst_owner, &
                              dst_global_index, send_glb2loc_index, &
                              src_n_points, src_owner, src_global_index, &
                              inplace, comm)

   CLASS(t_comm_pattern_yaxt), INTENT(OUT) :: p_pat
   INTEGER, INTENT(IN)           :: dst_n_points        ! Total number of points
   INTEGER, INTENT(IN)           :: dst_owner(:)        ! Owner of every point
   INTEGER, INTENT(IN)           :: dst_global_index(:) ! Global index of every point
   TYPE(t_glb2loc_index_lookup), INTENT(IN) :: send_glb2loc_index
                                                        ! global to local index
                                                        ! lookup information
                                                        ! of the SENDER array
   INTEGER, INTENT(IN)           :: src_n_points        ! Total number of points
   INTEGER, INTENT(IN)           :: src_owner(:)        ! Owner of every point
   INTEGER, INTENT(IN)           :: src_global_index(:) ! Global index of every point

   ! FIXME: We might want to do this more elegant. And we might want to check
   ! whether xt_int_kind is compatible with the size of integer.
   INTEGER(xt_int_kind), ALLOCATABLE :: src_global_index_cpy(:)
   INTEGER(xt_int_kind), ALLOCATABLE :: dst_global_index_cpy(:)

   LOGICAL, OPTIONAL, INTENT(IN) :: inplace
   INTEGER, OPTIONAL, INTENT(IN) :: comm

   TYPE(xt_idxlist) :: src_idxlist, dst_idxlist
   TYPE(xt_com_list), ALLOCATABLE :: src_com(:), dst_com(:)
   INTEGER :: src_com_size, dst_com_size

   INTEGER :: dst_count, src_count, accum_dst, n_dst_com, n_src_com
   INTEGER(xt_int_kind), ALLOCATABLE :: receive_indices(:), send_indices(:)
   INTEGER, ALLOCATABLE :: dst_indices_displ(:), src_indices_displ(:), &
        src_count_per_rank(:), dst_count_per_rank(:)

   INTEGER :: i, ierror, dst_owner_i
   CHARACTER(len=*), PARAMETER :: routine = modname//'::setup_comm_pattern'
   INTEGER :: pcomm, comm_size, comm_rank
   LOGICAL :: any_dst_owner_lt_0

!-----------------------------------------------------------------------

   IF (PRESENT(comm)) THEN
     pcomm = comm
   ELSE
     pcomm = p_comm_work
   END IF
   CALL mpi_comm_dup(pcomm, p_pat%comm, ierror)
   comm_size = p_comm_size(p_pat%comm)
   comm_rank = p_comm_rank(p_pat%comm)
   ALLOCATE(src_indices_displ(0:comm_size), dst_indices_displ(0:comm_size), &
        src_count_per_rank(0:comm_size-1), dst_count_per_rank(0:comm_size-1))
   IF (.NOT. xt_initialized()) CALL xt_initialize(p_pat%comm)

   IF (ierror /= mpi_success) &
        CALL finish(routine, 'failed to duplicate communicator')
   CALL xt_mpi_comm_mark_exclusive(p_pat%comm)

   dst_count_per_rank = 0

   DO i = 1, dst_n_points
     dst_owner_i = dst_owner(i)
     IF (dst_owner_i >= 0) &
       dst_count_per_rank(dst_owner_i) = dst_count_per_rank(dst_owner_i) + 1
   END DO

   CALL p_alltoall(dst_count_per_rank, src_count_per_rank, p_pat%comm)

   dst_count = 0
   n_dst_com = 0
   src_count = 0
   n_src_com = 0
   DO i = 0, comm_size-1
     dst_indices_displ(i) = dst_count
     dst_count = dst_count + dst_count_per_rank(i)
     n_dst_com = n_dst_com + MERGE(1, 0, dst_count_per_rank(i) > 0)
     src_indices_displ(i) = src_count
     src_count = src_count + src_count_per_rank(i)
     n_src_com = n_src_com + MERGE(1, 0, src_count_per_rank(i) > 0)
   END DO
   src_indices_displ(comm_size) = src_count

   ALLOCATE(receive_indices(dst_count), send_indices(src_count))

   any_dst_owner_lt_0 = .FALSE.
   DO i = 1, dst_n_points
     dst_owner_i = dst_owner(i)
     any_dst_owner_lt_0 = any_dst_owner_lt_0 .OR. dst_owner_i < 0
     IF (dst_owner_i >= 0) THEN
       dst_indices_displ(dst_owner_i) = dst_indices_displ(dst_owner_i) + 1
       receive_indices(dst_indices_displ(dst_owner_i)) &
         = INT(dst_global_index(i), xt_int_kind)
     END IF
   END DO

   IF (any_dst_owner_lt_0) THEN
     ALLOCATE(p_pat%dst_mask(dst_n_points))
     p_pat%dst_mask = dst_owner >= 0
   END IF

   accum_dst = 0
   DO i = 0, comm_size-1
     dst_indices_displ(i) = accum_dst
     accum_dst = accum_dst + dst_count_per_rank(i)
   END DO
   dst_indices_displ(comm_size) = accum_dst

   CALL p_alltoallv(receive_indices, dst_count_per_rank, dst_indices_displ, &
      &             send_indices, src_count_per_rank, src_indices_displ, &
      &             p_pat%comm)

   DEALLOCATE(src_count_per_rank, dst_count_per_rank)
   ALLOCATE(src_com(n_src_com), dst_com(n_dst_com))

   src_com_size = 0
   dst_com_size = 0

   DO i = 0, comm_size-1

     src_count = src_indices_displ(i+1) - src_indices_displ(i)
     IF (src_count > 0) THEN

       src_com_size = src_com_size + 1
       src_com(src_com_size)%rank = i
       src_com(src_com_size)%list = &
            xt_idxvec_new(send_indices(src_indices_displ(i)+1: &
            &                          src_indices_displ(i+1)))
      END IF

      dst_count = dst_indices_displ(i+1) - dst_indices_displ(i)
      IF (dst_count > 0) THEN

         dst_com_size = dst_com_size + 1
         dst_com(dst_com_size)%rank = i
         dst_com(dst_com_size)%list = &
              xt_idxvec_new(receive_indices(dst_indices_displ(i)+1: &
              &                             dst_indices_displ(i+1)))
      END IF
   END DO

   p_pat%src_n_points = src_n_points
   p_pat%dst_n_points = dst_n_points

   ALLOCATE(src_global_index_cpy(src_n_points))
   src_global_index_cpy = MERGE(src_global_index(1:src_n_points), -1,&
                                  src_owner(1:src_n_points) == comm_rank)
   src_idxlist = xt_idxvec_new(src_global_index_cpy)

   ALLOCATE(dst_global_index_cpy(dst_n_points))
   dst_global_index_cpy = dst_global_index(1:dst_n_points)
   IF (ALLOCATED(p_pat%dst_mask)) THEN
     dst_idxlist = xt_idxvec_new(PACK(dst_global_index_cpy, p_pat%dst_mask))
   ELSE
     dst_idxlist = xt_idxvec_new(dst_global_index_cpy)
   END IF

   p_pat%xmap = xt_xmap_intersection_new(src_com_size, src_com, dst_com_size, &
      &                                  dst_com, src_idxlist, dst_idxlist, &
      &                                  p_pat%comm)
   CALL xt_idxlist_delete(dst_idxlist)
   CALL xt_idxlist_delete(src_idxlist)
   DO i = 1, src_com_size
      CALL xt_idxlist_delete(src_com(i)%list)
   END DO
   DO i = 1, dst_com_size
      CALL xt_idxlist_delete(dst_com(i)%list)
   END DO
   IF (PRESENT(inplace)) THEN
     p_pat%inplace = inplace
   ELSE
     p_pat%inplace = .false.
   END IF

END SUBROUTINE setup_comm_pattern


  SUBROUTINE setup_comm_pattern2(p_pat, comm, recv_msg, send_msg, &
       glb2loc_index_recv, glb2loc_index_send, inplace)
    CLASS(t_comm_pattern_yaxt), INTENT(out) :: p_pat
    INTEGER, INTENT(in) :: comm
    TYPE(xfer_list), INTENT(in) :: recv_msg(:), send_msg(:)
    TYPE(t_glb2loc_index_lookup), INTENT(IN) :: glb2loc_index_recv, &
         glb2loc_index_send
    LOGICAL, OPTIONAL, INTENT(in) :: inplace
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
    CONTIGUOUS :: recv_msg, send_msg
#endif
    INTEGER :: i, np_recv, np_send, nlocal
    TYPE(xt_com_list), ALLOCATABLE :: src_com(:), dst_com(:)
    TYPE(xt_idxlist) :: src_idxlist, dst_idxlist
    CHARACTER(len=*), PARAMETER :: routine &
         = 'mo_communication_yaxt::setup_comm_pattern2'

    np_send = SIZE(send_msg)
    np_recv = SIZE(recv_msg)
    ALLOCATE(src_com(np_send), dst_com(np_recv))
    nlocal = SIZE(glb2loc_index_send%outer_glb_index) &
         + SIZE(glb2loc_index_send%inner_glb_index)
    p_pat%src_n_points = nlocal
    CALL compose_lists(src_com, src_idxlist, nlocal, glb2loc_index_send, &
         send_msg)
    nlocal = SIZE(glb2loc_index_recv%outer_glb_index) &
         + SIZE(glb2loc_index_recv%inner_glb_index)
    p_pat%dst_n_points = nlocal
    CALL compose_lists(dst_com, dst_idxlist, nlocal, glb2loc_index_recv, &
         recv_msg, p_pat%dst_mask)
    p_pat%xmap = xt_xmap_intersection_new(src_com, dst_com, src_idxlist, &
         dst_idxlist, comm)
    CALL xt_idxlist_delete(dst_idxlist)
    CALL xt_idxlist_delete(src_idxlist)
    DO i = 1, np_send
      CALL xt_idxlist_delete(src_com(i)%list)
    END DO
    DO i = 1, np_recv
      CALL xt_idxlist_delete(dst_com(i)%list)
    END DO
    IF (PRESENT(inplace)) THEN
      p_pat%inplace = inplace
    ELSE
      p_pat%inplace = .FALSE.
    END IF
  CONTAINS
    SUBROUTINE compose_lists(com, list, nlocal, glb2loc_index, msg, mask)
      TYPE(xt_com_list), INTENT(out) :: com(:)
      TYPE(xt_idxlist), INTENT(out) :: list
      INTEGER, INTENT(in) :: nlocal
      TYPE(xfer_list), INTENT(in) :: msg(:)
      type(t_glb2loc_index_lookup), INTENT(IN) :: glb2loc_index
      LOGICAL, OPTIONAL, ALLOCATABLE, INTENT(out) :: mask(:)
      INTEGER(xt_int_kind) :: indices(nlocal)
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
      CONTIGUOUS :: com, msg
#endif
      INTEGER :: np, nidx, i, j, jl, glbidx
      indices = -1_xt_int_kind
      np = SIZE(msg)
      DO i = 1, np
        com(i)%rank = INT(msg(i)%rank, c_int)
        com(i)%list = xt_idxvec_new(int(msg(i)%glob_idx, xt_int_kind))
        nidx = SIZE(msg(i)%glob_idx)
        DO j = 1, nidx
          glbidx = msg(i)%glob_idx(j)
          jl = get_local_index(glb2loc_index, glbidx)
          indices(jl) = INT(glbidx, xt_int_kind)
        END DO
      END DO
      IF (PRESENT(mask)) THEN
        ALLOCATE(mask(nlocal))
        mask = indices >= 0
        j = 1
        i = 1
        DO WHILE (i < nlocal .AND. .NOT. mask(i))
          i = i + 1
        END DO
        IF (.NOT. mask(i)) THEN
          nidx = 0
        ELSE
          DO i = i, nlocal
            IF (mask(i)) THEN
              indices(j) = indices(i)
              j = j + 1
            END IF
          END DO
          nidx = j - 1
        END IF
      ELSE
        nidx = nlocal
      END IF
      list = xt_idxvec_new(indices(1:nidx))
    END SUBROUTINE compose_lists

  END SUBROUTINE setup_comm_pattern2


FUNCTION generate_single_field_redist(p_pat, dst_nlev, src_nlev, nshift, &
  &                                   mpi_type, src_is_blocked, dst_is_blocked)

  CLASS(t_comm_pattern_yaxt), INTENT(IN) :: p_pat
  INTEGER, INTENT(IN) :: dst_nlev
  INTEGER, INTENT(IN) :: src_nlev
  INTEGER, INTENT(IN) :: nshift
  INTEGER, INTENT(IN) :: mpi_type
  LOGICAL, OPTIONAL, INTENT(IN) :: src_is_blocked, dst_is_blocked

  TYPE(xt_redist) :: generate_single_field_redist

  INTEGER :: i
  INTEGER, ALLOCATABLE :: src_offsets(:)
  INTEGER, ALLOCATABLE, TARGET :: dst_offsets(:)
  INTEGER, POINTER :: p_dst_offsets(:)
  INTEGER(p_address_kind) :: src_extent, dst_extent
  TYPE(xt_redist) :: redist

  IF (dst_nlev > src_nlev) THEN
     CALL finish('generate_single_field_redist', &
                 'ERROR: receive array has more levels than send array.')
  ENDIF

  IF ((dst_nlev == 1) .AND. (src_nlev == 1) .AND. &
      (.NOT. ALLOCATED(p_pat%dst_mask))) THEN

    redist = xt_redist_p2p_new(p_pat%xmap, mpi_type)

  ELSE
    ALLOCATE(src_offsets(p_pat%src_n_points), dst_offsets(p_pat%dst_n_points))

    CALL generate_offs_displs(src_offsets, src_extent, p_pat%src_n_points, &
      &                       src_nlev, nshift, mpi_type, src_is_blocked)
    CALL generate_offs_displs(dst_offsets, dst_extent, p_pat%dst_n_points, &
      &                       dst_nlev, nshift, mpi_type, dst_is_blocked)

    IF (ALLOCATED(p_pat%dst_mask)) THEN
      ALLOCATE(p_dst_offsets(COUNT(p_pat%dst_mask)))
      p_dst_offsets = PACK(dst_offsets, p_pat%dst_mask)
    ELSE
      p_dst_offsets => dst_offsets
    END IF

    redist = xt_redist_p2p_off_new(p_pat%xmap, src_offsets, p_dst_offsets, &
      &                            mpi_type)

    IF (ALLOCATED(p_pat%dst_mask)) DEALLOCATE(p_dst_offsets)

    DEALLOCATE(src_offsets, dst_offsets)
  END IF

  IF ((dst_nlev == 1) .AND. (src_nlev == 1)) THEN

    generate_single_field_redist = redist

  ELSE

    generate_single_field_redist = &
      xt_redist_repeat_new(redist, src_extent, dst_extent, dst_nlev - nshift, &
      &                    (/(i, i = nshift, dst_nlev - 1)/))

    CALL xt_redist_delete(redist)
  END IF

CONTAINS

  SUBROUTINE generate_offs_displs(offsets, div_disp, n, nlev, nshift, &
    &                             mpi_type, is_blocked)

    INTEGER, INTENT(IN) :: mpi_type
    INTEGER, INTENT(IN) :: n, nlev, nshift
    LOGICAL, OPTIONAL, INTENT(IN) :: is_blocked
    INTEGER, INTENT(OUT) :: offsets(n)
    INTEGER(p_address_kind), INTENT(OUT) :: div_disp

    LOGICAL :: is_blk
    INTEGER :: i, ierror
    INTEGER(p_address_kind) :: lb, extent


    IF (PRESENT(is_blocked)) THEN
      is_blk = is_blocked
    ELSE
      is_blk = .TRUE.
    END IF

    CALL mpi_type_get_extent(mpi_type, lb, extent, ierror)

    IF (is_blk) THEN
      DO i = 1, n
        offsets(i) = idx_no(i) + (blk_no(i) - 1) * (nproma * nlev) - 1
      END DO
      div_disp = nproma * extent
    ELSE
      DO i = 1, n
        offsets(i) = (i - 1) * nlev
      END DO
      div_disp = extent
    END IF
  END SUBROUTINE generate_offs_displs

END FUNCTION generate_single_field_redist

FUNCTION comm_pattern_get_redist(p_pat, nfields, nlev, mpi_type, &
  &                              nshift)

  CLASS(t_comm_pattern_yaxt), INTENT(INOUT) :: p_pat
  INTEGER, INTENT(IN) :: nfields
  ! nlev(i, 1) is number of levels for dst array i
  ! nlev(i, 2) is number of levels for src array i
  INTEGER, INTENT(IN) :: nlev(nfields, 2)
  INTEGER, INTENT(IN) :: mpi_type
  INTEGER, OPTIONAL, INTENT(IN) :: nshift(nfields)

  INTEGER :: kshift(nfields)
  TYPE(xt_redist) :: comm_pattern_get_redist
  TYPE(t_comm_pattern_redist), ALLOCATABLE :: tmp_redists(:)

  INTEGER :: i, n

  TYPE(xt_redist) :: redists(nfields)

  IF (PRESENT(nshift)) THEN
    kshift = nshift
  ELSE
    kshift = 0
  END IF

  IF (ALLOCATED(p_pat%redists)) THEN

    n = SIZE(p_pat%redists)
    DO i = 1, n

      IF ((p_pat%redists(i)%nfields == nfields) .AND. &
        & (p_pat%redists(i)%mpi_type == mpi_type)) THEN

        IF (ALL(p_pat%redists(i)%dst_nlev(:) == nlev(:,1)) .AND. &
            ALL(p_pat%redists(i)%src_nlev(:) == nlev(:,2)) .AND. &
            ALL(p_pat%redists(i)%nshift(:) == kshift(:))) THEN

          comm_pattern_get_redist = p_pat%redists(i)%redist
          RETURN
        END IF
      END IF

    END DO

    ALLOCATE(tmp_redists(n + 1))
    tmp_redists(1:n) = p_pat%redists
    CALL MOVE_ALLOC(tmp_redists, p_pat%redists)
    n = n + 1

  ELSE
    n = 1
    ALLOCATE(p_pat%redists(n))
  END IF

  p_pat%redists(n)%nfields = nfields
  ALLOCATE(p_pat%redists(n)%dst_nlev(nfields), &
    &      p_pat%redists(n)%src_nlev(nfields))
  p_pat%redists(n)%dst_nlev = nlev(:,1)
  p_pat%redists(n)%src_nlev = nlev(:,2)
  ALLOCATE(p_pat%redists(n)%nshift(nfields))
  p_pat%redists(n)%nshift = kshift
  p_pat%redists(n)%mpi_type = mpi_type

  IF (nfields == 1) THEN

    p_pat%redists(n)%redist = &
      generate_single_field_redist(p_pat, nlev(1, 1), nlev(1, 2), kshift(1), &
        &                          mpi_type)

  ELSE

    DO i = 1, nfields
      redists(i) = generate_single_field_redist(p_pat, nlev(i, 1), &
        &                                       nlev(i, 2), kshift(i), &
        &                                       mpi_type)
    END DO

    p_pat%redists(n)%redist = &
      xt_redist_collection_new(redists, nfields, -1, p_pat%comm)

    CALL xt_redist_delete(redists)
  END IF

  comm_pattern_get_redist = p_pat%redists(n)%redist

END FUNCTION comm_pattern_get_redist

FUNCTION comm_pattern_get_contiguous_data_type(p_pat, base_type, count)

  CLASS(t_comm_pattern_yaxt), INTENT(INOUT) :: p_pat
  INTEGER, INTENT(IN) :: base_type
  INTEGER, INTENT(IN) :: count
  INTEGER :: comm_pattern_get_contiguous_data_type

  TYPE(t_comm_pattern_contiguous_data_type), ALLOCATABLE :: &
    tmp_contiguous_data_types(:)

  INTEGER :: i, n, ierror

  IF (ALLOCATED(p_pat%contiguous_data_types)) THEN

    n = SIZE(p_pat%contiguous_data_types)
    DO i = 1, n

      IF ((p_pat%contiguous_data_types(i)%base_type == base_type) .AND. &
        & (p_pat%contiguous_data_types(i)%count == count)) THEN

        comm_pattern_get_contiguous_data_type = &
          p_pat%contiguous_data_types(i)%type
        RETURN
      END IF

    END DO

    ALLOCATE(tmp_contiguous_data_types(n + 1))
    tmp_contiguous_data_types(1:n) = p_pat%contiguous_data_types
    CALL MOVE_ALLOC(tmp_contiguous_data_types, p_pat%contiguous_data_types)
    n = n + 1

  ELSE
    n = 1
    ALLOCATE(p_pat%contiguous_data_types(n))
  END IF

  p_pat%contiguous_data_types(n)%base_type = base_type
  p_pat%contiguous_data_types(n)%count = count

  CALL MPI_Type_contiguous(count, base_type, p_pat%contiguous_data_types(n)%type, &
    &                      ierror)
  CALL MPI_Type_commit(p_pat%contiguous_data_types(n)%type, ierror)

  comm_pattern_get_contiguous_data_type = p_pat%contiguous_data_types(n)%type

END FUNCTION comm_pattern_get_contiguous_data_type

!-------------------------------------------------------------------------

FUNCTION comm_pattern_collection_get_redist(p_pat_coll, nfields, dst_nlev, &
  &                                         src_nlev)

  CLASS(t_comm_pattern_collection_yaxt), INTENT(INOUT) :: p_pat_coll
  INTEGER, INTENT(IN) :: nfields
  INTEGER, INTENT(IN) :: dst_nlev(nfields)
  INTEGER, INTENT(IN) :: src_nlev(nfields)

  TYPE(xt_redist) :: comm_pattern_collection_get_redist
  TYPE(t_comm_pattern_coll_redist), ALLOCATABLE :: tmp_redists(:)

  INTEGER :: i, j, n

  TYPE(xt_redist) :: redists(nfields * SIZE(p_pat_coll%patterns))

  IF (ALLOCATED(p_pat_coll%redists)) THEN

    n = SIZE(p_pat_coll%redists)
    DO i = 1, n
      IF ((p_pat_coll%redists(i)%nfields == nfields)) THEN

        IF (ALL(p_pat_coll%redists(i)%dst_nlev(:) == dst_nlev(:)) .AND. &
            ALL(p_pat_coll%redists(i)%src_nlev(:) == src_nlev(:))) THEN
          comm_pattern_collection_get_redist = p_pat_coll%redists(i)%redist
          RETURN
        END IF
      END IF
    END DO

    ALLOCATE(tmp_redists(n + 1))
    tmp_redists(1:n) = p_pat_coll%redists
    CALL MOVE_ALLOC(tmp_redists, p_pat_coll%redists)
    n = n + 1
  ELSE
    n = 1
    ALLOCATE(p_pat_coll%redists(n))
  END IF

  p_pat_coll%redists(n)%nfields = nfields
  ALLOCATE(p_pat_coll%redists(n)%dst_nlev(nfields), &
    &      p_pat_coll%redists(n)%src_nlev(nfields))
  p_pat_coll%redists(n)%dst_nlev = dst_nlev
  p_pat_coll%redists(n)%src_nlev = src_nlev

  DO i = 1, nfields
    DO j = 1, SIZE(p_pat_coll%patterns)
      redists(j + (i-1)*SIZE(p_pat_coll%patterns)) = &
        generate_single_field_redist(p_pat_coll%patterns(j)%p, dst_nlev(i), &
          &                          src_nlev(i), 0, p_real_dp, &
          &                          src_is_blocked=.FALSE.)
    END DO
  END DO

  p_pat_coll%redists(n)%redist = &
    xt_redist_collection_new(redists, -1, p_pat_coll%patterns(1)%p%comm)

  CALL xt_redist_delete(redists)

  comm_pattern_collection_get_redist = p_pat_coll%redists(n)%redist

END FUNCTION comm_pattern_collection_get_redist

SUBROUTINE setup_comm_pattern_collection(pattern_collection, patterns)

  CLASS(t_comm_pattern_collection_yaxt), INTENT(OUT) :: pattern_collection
  TYPE(t_p_comm_pattern), INTENT(IN) :: patterns(:)

  CHARACTER(len=*), PARAMETER :: &
       routine = modname//'::setup_comm_pattern_collection'
  INTEGER :: i, n

  ALLOCATE(pattern_collection%patterns(SIZE(patterns)))

  DO i = 1, SIZE(patterns)
    SELECT TYPE (pattern_yaxt => patterns(i)%p)
      TYPE IS (t_comm_pattern_yaxt)
        pattern_collection%patterns(i)%p => pattern_yaxt
      CLASS DEFAULT
        CALL finish(routine, "wrong t_comm_pattern type")
    END SELECT
  END DO

END SUBROUTINE setup_comm_pattern_collection

!-------------------------------------------------------------------------
!
!>
!! Deletes a communication pattern
!!
!! @par Revision History
!! Initial version by Rainer Johanni, Oct 2011
!! yaxt version by Moritz Hanke, April 2016
!!
!
SUBROUTINE delete_comm_pattern(p_pat)

   CLASS(t_comm_pattern_yaxt), INTENT(INOUT) :: p_pat

   CHARACTER(len=*), PARAMETER :: routine = modname//'::delete_comm_pattern'
   INTEGER :: i, n, ierror

   ! deallocate arrays
   IF (ALLOCATED(p_pat%dst_mask)) DEALLOCATE(p_pat%dst_mask)
   IF (ALLOCATED(p_pat%contiguous_data_types)) THEN
     n = SIZE(p_pat%contiguous_data_types)
     DO i = 1, n
      CALL MPI_Type_free(p_pat%contiguous_data_types(i)%type, ierror)
     END DO
     DEALLOCATE(p_pat%contiguous_data_types)
   END IF
   IF (ALLOCATED(p_pat%redists)) THEN
     n = SIZE(p_pat%redists)
     DO i = 1, n
       DEALLOCATE(p_pat%redists(i)%dst_nlev)
       DEALLOCATE(p_pat%redists(i)%src_nlev)
       DEALLOCATE(p_pat%redists(i)%nshift)
       CALL xt_redist_delete(p_pat%redists(i)%redist)
     END DO
     DEALLOCATE(p_pat%redists)
   END IF
   CALL xt_xmap_delete(p_pat%xmap)
   IF (p_pat%comm /= mpi_comm_null) THEN
     CALL mpi_comm_free(p_pat%comm, ierror)
     IF (ierror /= mpi_success) &
       CALL finish(routine, 'failed to duplicate communicator')
   END IF
END SUBROUTINE delete_comm_pattern

!-------------------------------------------------------------------------

SUBROUTINE delete_comm_pattern_collection(pattern_collection)

   CLASS(t_comm_pattern_collection_yaxt), INTENT(INOUT) :: pattern_collection

   INTEGER :: i, n

   IF (ALLOCATED(pattern_collection%redists)) THEN
     n = SIZE(pattern_collection%redists)
     DO i = 1, n
       DEALLOCATE(pattern_collection%redists(i)%dst_nlev)
       DEALLOCATE(pattern_collection%redists(i)%src_nlev)
       CALL xt_redist_delete(pattern_collection%redists(i)%redist)
     END DO
     DEALLOCATE(pattern_collection%redists)
   END IF

   DO i = 1, SIZE(pattern_collection%patterns)
     CALL pattern_collection%patterns(i)%p%delete()
     DEALLOCATE(pattern_collection%patterns(i)%p)
   END DO
   DEALLOCATE(pattern_collection%patterns)

END SUBROUTINE delete_comm_pattern_collection

!-------------------------------------------------------------------------
!
!
!>
!! Does data exchange according to a communication pattern (in p_pat).
!!
!!
!! @par Revision History
!! Initial version by Rainer Johanni, Nov 2009
!! Modified by Guenther Zaengl for vectorization
!! yaxt version by Moritz Hanke, April 2016
!!
!================================================================================================
! REAL SECTION ----------------------------------------------------------------------------------
!
SUBROUTINE exchange_data_r3d(p_pat, recv, send, add)

   CLASS(t_comm_pattern_yaxt), INTENT(INOUT) :: p_pat
   REAL(dp), INTENT(INOUT), TARGET :: recv(:,:,:)
   REAL(dp), INTENT(IN), OPTIONAL, TARGET :: send(:,:,:)
   REAL(dp), INTENT(IN), OPTIONAL, TARGET :: add (:,:,:)

   TYPE(xt_redist) :: redist

   INTEGER :: i, j, k, nlev(1,2), m, n, o

   IF(SIZE(recv,1) /= nproma) THEN
     CALL finish('exchange_data_r3d','Illegal first dimension of data array')
   ENDIF

   IF (itype_exch_barrier == 1 .OR. itype_exch_barrier == 3) THEN
     IF (activate_sync_timers) CALL timer_start(timer_barrier)
     CALL p_barrier(p_pat%comm)
     IF (activate_sync_timers) CALL timer_stop(timer_barrier)
   ENDIF

   start_sync_timer(timer_exch_data)

   m = SIZE(recv, 1)
   n = SIZE(recv, 2)
   o = SIZE(recv, 3)

   nlev(1,1) = n
   IF (PRESENT(send)) THEN
     nlev(1,2) = SIZE(send,2)
   ELSE
     nlev(1,2) = n
   END IF
   redist = comm_pattern_get_redist(p_pat, 1, nlev, p_real_dp)

   IF (PRESENT(send)) THEN
     CALL xt_redist_s_exchange1_contiguous(redist, send, recv)
   ELSE
     IF (p_pat%inplace) THEN
       CALL xt_redist_s_exchange1_contiguous_inplace(redist, recv)
     ELSE
       ! make copy of recv
       CALL xt_redist_s_exchange1_contiguous_copy(redist, recv, SIZE(recv))
     END IF
   END IF

   IF (itype_exch_barrier == 2 .OR. itype_exch_barrier == 3) THEN
     IF (activate_sync_timers) CALL timer_start(timer_barrier)
     CALL p_barrier(p_pat%comm)
     IF (activate_sync_timers) CALL timer_stop(timer_barrier)
   ENDIF

   IF (PRESENT(add)) THEN
     IF (ALLOCATED(p_pat%dst_mask)) THEN
       DO k = 1, o
         DO j = 1, n
           DO i = 1, m
             IF (p_pat%dst_mask(idx_1d(i,j))) THEN
               recv(i,j,k) = recv(i,j,k) + add(i,j,k)
             END IF
           END DO
         END DO
       END DO
     ELSE
       DO k = 1, o
         DO j = 1, n
           DO i = 1, m
             recv(i,j,k) = recv(i,j,k) + add(i,j,k)
           END DO
         END DO
       END DO
     END IF
   END IF

   stop_sync_timer(timer_exch_data)

CONTAINS

  ! this wrapper is needed because we have to ensure that the input array are in
  ! contiguous memory (the keyword CONTIGUOUS is Fortran2008, which is not
  ! required by ICON)
  SUBROUTINE xt_redist_s_exchange1_contiguous(redist, send, recv)

    TYPE(xt_redist), INTENT(IN) :: redist
    REAL(dp), INTENT(INOUT), TARGET :: recv(*)
    REAL(dp), INTENT(IN), TARGET :: send(*)

    CALL xt_redist_s_exchange1(redist, c_loc(send), c_loc(recv))

  END SUBROUTINE

  SUBROUTINE xt_redist_s_exchange1_contiguous_inplace(redist, recv)

    TYPE(xt_redist), INTENT(IN) :: redist
    REAL(dp), INTENT(INOUT), TARGET :: recv(*)

    CALL xt_redist_s_exchange1(redist, c_loc(recv), c_loc(recv))

  END SUBROUTINE

    SUBROUTINE xt_redist_s_exchange1_contiguous_copy(redist, recv, n)

      TYPE(xt_redist), INTENT(IN) :: redist
      INTEGER, INTENT(in) :: n
      REAL(dp), INTENT(INOUT), TARGET :: recv(n)
      REAL(dp), TARGET :: send(n)

      send = recv
      CALL xt_redist_s_exchange1(redist, c_loc(send), c_loc(recv))

    END SUBROUTINE xt_redist_s_exchange1_contiguous_copy


END SUBROUTINE exchange_data_r3d

SUBROUTINE exchange_data_s3d(p_pat, recv, send, add)

   CLASS(t_comm_pattern_yaxt), INTENT(INOUT) :: p_pat
   REAL(sp), INTENT(INOUT), TARGET :: recv(:,:,:)
   REAL(sp), INTENT(IN), OPTIONAL, TARGET :: send(:,:,:)
   REAL(sp), INTENT(IN), OPTIONAL, TARGET :: add (:,:,:)

   TYPE(xt_redist) :: redist

   INTEGER :: i, j, k, nlev(1, 2), m, n, o

   IF(SIZE(recv,1) /= nproma) THEN
     CALL finish('exchange_data_s3d','Illegal first dimension of data array')
   ENDIF

   IF (itype_exch_barrier == 1 .OR. itype_exch_barrier == 3) THEN
     IF (activate_sync_timers) CALL timer_start(timer_barrier)
     CALL p_barrier(p_pat%comm)
     IF (activate_sync_timers) CALL timer_stop(timer_barrier)
   ENDIF

   start_sync_timer(timer_exch_data)

   m = SIZE(recv, 1)
   n = SIZE(recv, 2)
   o = SIZE(recv, 3)

   nlev(1,1) = n
   IF (PRESENT(send)) THEN
     nlev(1, 2) = SIZE(send,2)
   ELSE
     nlev(1, 2) = n
   END IF
   redist = comm_pattern_get_redist(p_pat, 1, nlev, p_real_sp)

   IF (PRESENT(send)) THEN
     CALL xt_redist_s_exchange1_contiguous(redist, send, recv)
   ELSE
     IF (p_pat%inplace) THEN
       CALL xt_redist_s_exchange1_contiguous_inplace(redist, recv)
     ELSE
       ! make copy of recv
       CALL xt_redist_s_exchange1_contiguous_copy(redist, recv, SIZE(recv))
     END IF
   END IF

   IF (itype_exch_barrier == 2 .OR. itype_exch_barrier == 3) THEN
     IF (activate_sync_timers) CALL timer_start(timer_barrier)
     CALL p_barrier(p_pat%comm)
     IF (activate_sync_timers) CALL timer_stop(timer_barrier)
   ENDIF

   IF (PRESENT(add)) THEN
     IF (ALLOCATED(p_pat%dst_mask)) THEN
       DO k = 1, o
         DO j = 1, n
           DO i = 1, m
             IF (p_pat%dst_mask(idx_1d(i,j))) THEN
               recv(i,j,k) = recv(i,j,k) + add(i,j,k)
             END IF
           END DO
         END DO
       END DO
     ELSE
       DO k = 1, o
         DO j = 1, n
           DO i = 1, m
             recv(i,j,k) = recv(i,j,k) + add(i,j,k)
           END DO
         END DO
       END DO
     END IF
   END IF

   stop_sync_timer(timer_exch_data)

CONTAINS

  ! this wrapper is needed because we have to ensure that the input array are in
  ! contiguous memory (the keyword CONTIGUOUS is Fortran2008, which is not
  ! required by ICON)
  SUBROUTINE xt_redist_s_exchange1_contiguous(redist, send, recv)

    TYPE(xt_redist), INTENT(IN) :: redist
    REAL(sp), INTENT(INOUT), TARGET :: recv(*)
    REAL(sp), INTENT(IN), TARGET :: send(*)

    CALL xt_redist_s_exchange1(redist, c_loc(send), c_loc(recv))

  END SUBROUTINE

  SUBROUTINE xt_redist_s_exchange1_contiguous_inplace(redist, recv)

    TYPE(xt_redist), INTENT(IN) :: redist
    REAL(sp), INTENT(INOUT), TARGET :: recv(*)

    CALL xt_redist_s_exchange1(redist, c_loc(recv), c_loc(recv))

  END SUBROUTINE

    SUBROUTINE xt_redist_s_exchange1_contiguous_copy(redist, recv, n)

      TYPE(xt_redist), INTENT(IN) :: redist
      INTEGER, INTENT(in) :: n
      REAL(sp), INTENT(INOUT), TARGET :: recv(n)
      REAL(sp), TARGET :: send(n)

      send = recv
      CALL xt_redist_s_exchange1(redist, c_loc(send), c_loc(recv))

    END SUBROUTINE xt_redist_s_exchange1_contiguous_copy

END SUBROUTINE exchange_data_s3d

!================================================================================================
! INTEGER SECTION -------------------------------------------------------------------------------
!
SUBROUTINE exchange_data_i3d(p_pat, recv, send, add)

   CLASS(t_comm_pattern_yaxt), INTENT(INOUT) :: p_pat
   INTEGER, INTENT(INOUT), TARGET :: recv(:,:,:)
   INTEGER, INTENT(IN), OPTIONAL, TARGET :: send(:,:,:)
   INTEGER, INTENT(IN), OPTIONAL, TARGET :: add (:,:,:)

   TYPE(xt_redist) :: redist

   INTEGER :: i, j, k, nlev(1, 2), m, n, o

   IF(SIZE(recv,1) /= nproma) THEN
     CALL finish('exchange_data_i3d','Illegal first dimension of data array')
   ENDIF

   IF (itype_exch_barrier == 1 .OR. itype_exch_barrier == 3) THEN
     IF (activate_sync_timers) CALL timer_start(timer_barrier)
     CALL p_barrier(p_pat%comm)
     IF (activate_sync_timers) CALL timer_stop(timer_barrier)
   ENDIF

   start_sync_timer(timer_exch_data)

   m = SIZE(recv, 1)
   n = SIZE(recv, 2)
   o = SIZE(recv, 3)

   nlev(1, 1) = n
   IF (PRESENT(send)) THEN
     nlev(1, 2) = SIZE(send,2)
   ELSE
     nlev(1, 2) = n
   END IF
   redist = comm_pattern_get_redist(p_pat, 1, nlev, p_int)

   IF (PRESENT(send)) THEN
     CALL xt_redist_s_exchange1_contiguous(redist, send, recv)
   ELSE
     IF (p_pat%inplace) THEN
       CALL xt_redist_s_exchange1_contiguous_inplace(redist, recv)
     ELSE
       ! make copy of recv
       CALL xt_redist_s_exchange1_contiguous_copy(redist, recv, SIZE(recv))
     END IF
   END IF

   IF (itype_exch_barrier == 2 .OR. itype_exch_barrier == 3) THEN
     IF (activate_sync_timers) CALL timer_start(timer_barrier)
     CALL p_barrier(p_pat%comm)
     IF (activate_sync_timers) CALL timer_stop(timer_barrier)
   ENDIF

   IF (PRESENT(add)) THEN
     IF (ALLOCATED(p_pat%dst_mask)) THEN
       DO k = 1, o
         DO j = 1, n
           DO i = 1, m
             IF (p_pat%dst_mask(idx_1d(i,j))) THEN
               recv(i,j,k) = recv(i,j,k) + add(i,j,k)
             END IF
           END DO
         END DO
       END DO
     ELSE
       DO k = 1, o
         DO j = 1, n
           DO i = 1, m
             recv(i,j,k) = recv(i,j,k) + add(i,j,k)
           END DO
         END DO
       END DO
     END IF
   END IF

   stop_sync_timer(timer_exch_data)

CONTAINS

  ! this wrapper is needed because we have to ensure that the input array are in
  ! contiguous memory (the keyword CONTIGUOUS is Fortran2008, which is not
  ! required by ICON)
  SUBROUTINE xt_redist_s_exchange1_contiguous(redist, send, recv)

    TYPE(xt_redist), INTENT(IN) :: redist
    INTEGER, INTENT(INOUT), TARGET :: recv(*)
    INTEGER, INTENT(IN), TARGET :: send(*)

    CALL xt_redist_s_exchange1(redist, c_loc(send), c_loc(recv))

  END SUBROUTINE

  SUBROUTINE xt_redist_s_exchange1_contiguous_inplace(redist, recv)

    TYPE(xt_redist), INTENT(IN) :: redist
    INTEGER, INTENT(INOUT), TARGET :: recv(*)

    CALL xt_redist_s_exchange1(redist, c_loc(recv), c_loc(recv))

  END SUBROUTINE

    SUBROUTINE xt_redist_s_exchange1_contiguous_copy(redist, recv, n)

      TYPE(xt_redist), INTENT(IN) :: redist
      INTEGER, INTENT(in) :: n
      INTEGER, INTENT(INOUT), TARGET :: recv(n)
      INTEGER, TARGET :: send(n)

      send = recv
      CALL xt_redist_s_exchange1(redist, c_loc(send), c_loc(recv))

    END SUBROUTINE xt_redist_s_exchange1_contiguous_copy

END SUBROUTINE exchange_data_i3d

  SUBROUTINE exchange_data_l3d(p_pat, recv, send)
    CLASS(t_comm_pattern_yaxt), INTENT(INOUT) :: p_pat
    LOGICAL, INTENT(INOUT), TARGET :: recv(:,:,:)
    LOGICAL, INTENT(IN), OPTIONAL, TARGET :: send(:,:,:)

    TYPE(xt_redist) :: redist

    INTEGER :: nlev(1, 2), n

    IF(SIZE(recv,1) /= nproma) THEN
      CALL finish('exchange_data_i3d','Illegal first dimension of data array')
    ENDIF

    IF (itype_exch_barrier == 1 .OR. itype_exch_barrier == 3) THEN
      IF (activate_sync_timers) CALL timer_start(timer_barrier)
      CALL p_barrier(p_pat%comm)
      IF (activate_sync_timers) CALL timer_stop(timer_barrier)
    ENDIF

    start_sync_timer(timer_exch_data)

    n = SIZE(recv, 2)

    nlev(1, 1) = n
    IF (PRESENT(send)) THEN
      nlev(1, 2) = SIZE(send,2)
    ELSE
      nlev(1, 2) = n
    END IF
    redist = comm_pattern_get_redist(p_pat, 1, nlev, p_int)

    IF (PRESENT(send)) THEN
      CALL xt_redist_s_exchange1_contiguous(redist, send, recv)
    ELSE
      IF (p_pat%inplace) THEN
        CALL xt_redist_s_exchange1_contiguous_inplace(redist, recv)
      ELSE
        ! make copy of recv
        CALL xt_redist_s_exchange1_contiguous_copy(redist, recv, SIZE(recv))
      END IF
    END IF

    IF (itype_exch_barrier == 2 .OR. itype_exch_barrier == 3) THEN
      IF (activate_sync_timers) CALL timer_start(timer_barrier)
      CALL p_barrier(p_pat%comm)
      IF (activate_sync_timers) CALL timer_stop(timer_barrier)
    ENDIF

    stop_sync_timer(timer_exch_data)

  CONTAINS

    ! this wrapper is needed because we have to ensure that the input array are in
    ! contiguous memory (the keyword CONTIGUOUS is Fortran2008, which is not
    ! required by ICON)
    SUBROUTINE xt_redist_s_exchange1_contiguous(redist, send, recv)

      TYPE(xt_redist), INTENT(IN) :: redist
      LOGICAL, INTENT(INOUT), TARGET :: recv(*)
      LOGICAL, INTENT(IN), TARGET :: send(*)
      TYPE(c_ptr) :: send_ptr, recv_ptr

      CALL xt_slice_c_loc(send, send_ptr)
      CALL xt_slice_c_loc(recv, recv_ptr)
      CALL xt_redist_s_exchange1(redist, send_ptr, recv_ptr)

    END SUBROUTINE xt_redist_s_exchange1_contiguous

    SUBROUTINE xt_redist_s_exchange1_contiguous_inplace(redist, recv)

      TYPE(xt_redist), INTENT(IN) :: redist
      LOGICAL, INTENT(INOUT), TARGET :: recv(*)
      TYPE(c_ptr) :: recv_ptr

      CALL xt_slice_c_loc(recv, recv_ptr)
      CALL xt_redist_s_exchange1(redist, recv_ptr, recv_ptr)

    END SUBROUTINE xt_redist_s_exchange1_contiguous_inplace

    SUBROUTINE xt_redist_s_exchange1_contiguous_copy(redist, recv, n)

      TYPE(xt_redist), INTENT(IN) :: redist
      INTEGER, INTENT(in) :: n
      LOGICAL, INTENT(INOUT), TARGET :: recv(n)
      LOGICAL, TARGET :: send(n)
      TYPE(c_ptr) :: send_ptr, recv_ptr

      send = recv
      CALL xt_slice_c_loc(send, send_ptr)
      CALL xt_slice_c_loc(recv, recv_ptr)
      CALL xt_redist_s_exchange1(redist, send_ptr, recv_ptr)

    END SUBROUTINE xt_redist_s_exchange1_contiguous_copy

  END SUBROUTINE exchange_data_l3d

!>
!! Does data exchange according to a communication pattern (in p_pat).
!!
!!
!! @par Revision History
!! Initial version by Rainer Johanni, Nov 2009
!! Optimized version by Guenther Zaengl to process 4D fields or up to seven 3D fields
!! in one step
!! yaxt version by Moritz Hanke, April 2016
!!
SUBROUTINE exchange_data_mult_dp(p_pat, ndim2tot, &
   recv, send, nshift)

   CLASS(t_comm_pattern_yaxt), INTENT(INOUT) :: p_pat

   TYPE(t_ptr_3d), PTR_INTENT(in) :: recv(:)
   TYPE(t_ptr_3d), OPTIONAL, PTR_INTENT(in) :: send(:)

   INTEGER, INTENT(IN)           :: ndim2tot
   INTEGER, OPTIONAL, INTENT(IN) :: nshift

   INTEGER :: nlev(SIZE(recv), 2)
   LOGICAL :: iscont(SIZE(recv), 2)
   INTEGER :: nproma1, nblk, nl, cpy_shape(3), cpy_size
   INTEGER :: i, nfields
   LOGICAL :: lsend, nproma_mismatch_found, cpy_recv

   CHARACTER(len=*), PARAMETER :: &
     routine = "mo_communication::exchange_data_mult_dp"

!-----------------------------------------------------------------------

   IF (itype_exch_barrier == 1 .OR. itype_exch_barrier == 3) THEN
     IF (activate_sync_timers) CALL timer_start(timer_barrier)
     CALL p_barrier(p_pat%comm)
     IF (activate_sync_timers) CALL timer_stop(timer_barrier)
   ENDIF

   start_sync_timer(timer_exch_data)

   lsend = PRESENT(send)
   cpy_recv = .NOT. (lsend .OR. p_pat%inplace)

   nfields = SIZE(recv)
   nproma1 = SIZE(recv(1)%p, 1)
   nproma_mismatch_found = .FALSE.
   cpy_size = 0
   DO i = 1, nfields
     iscont(i, 1) = IS_CONTIGUOUS(recv(i)%p)
     cpy_shape = SHAPE(recv(i)%p)
     nproma_mismatch_found = nproma_mismatch_found &
       &          .OR. (i > 1 .AND. cpy_shape(1) /= nproma1)
     nl = cpy_shape(2)
     nblk = cpy_shape(3)
     nlev(i, 1) = nl
     cpy_size =   cpy_size &
       &        + nproma1 * nl * nblk &
       &             * (MERGE(0, 1, iscont(i, 1)) + MERGE(1, 0, cpy_recv))
   END DO
   IF (lsend) THEN
     DO i = 1, nfields
       iscont(i, 2) = IS_CONTIGUOUS(send(i)%p)
       cpy_shape = SHAPE(send(i)%p)
       nproma_mismatch_found = nproma_mismatch_found &
         &        .OR. cpy_shape(1) /= nproma1
       nl = cpy_shape(2)
       nblk = cpy_shape(3)
       cpy_size =   cpy_size &
         &        + nproma1 * MERGE(0, nl, iscont(i, 2)) * nblk
       nlev(i, 2) = nl
     END DO
   ELSE
     iscont(:, 2) = .TRUE.
     nlev(:, 2) = nlev(:, 1)
   END IF

   IF (nproma_mismatch_found) &
     CALL finish(routine, "inconsistent array shapes detected")

   CALL exchange_data_mult_dp_bottom(p_pat, cpy_size, nlev, iscont, &
        recv, send, nshift)

   stop_sync_timer(timer_exch_data)

END SUBROUTINE exchange_data_mult_dp

  !>
  !! Does data exchange according to a communication pattern (in p_pat).
  !!
  !!
  !! @par Revision History
  !! Initial version by Rainer Johanni, Nov 2009
  !! Optimized version by Guenther Zaengl to process 4D fields or up to seven 3D fields
  !! in one step
  !! yaxt version by Moritz Hanke, April 2016
  !!
  SUBROUTINE exchange_data_mult_dp_bottom(p_pat, cpy_size, nlev, iscont, &
    recv, send, nshift)

    CLASS(t_comm_pattern_yaxt), INTENT(INOUT) :: p_pat

    TYPE(t_ptr_3d), PTR_INTENT(in) :: recv(:)
    TYPE(t_ptr_3d), OPTIONAL, PTR_INTENT(in) :: send(:)

    INTEGER, INTENT(IN) :: cpy_size, nlev(SIZE(recv), 2)
    LOGICAL, INTENT(IN)           :: iscont(SIZE(recv), 2)
    INTEGER, OPTIONAL, INTENT(IN) :: nshift

    INTEGER :: i, nfields, cpy_psum, ofs, last, nproma, nblk, incr, nl, &
         kshift(SIZE(recv))
    LOGICAL :: lsend, cpy_recv
    REAL(dp), TARGET :: cpy_buf(cpy_size)
    REAL(dp), POINTER :: cpy(:, :, :)
    TYPE(c_ptr) :: src_data_cptr(SIZE(recv)), dst_data_cptr(SIZE(recv))
    TYPE(xt_redist) :: redist_coll

    lsend = PRESENT(send)
    cpy_recv = .NOT. (lsend .OR. p_pat%inplace)

    nfields = SIZE(recv)
    cpy_psum = 0
    nproma = SIZE(recv(1)%p, 1)
    ! set up C pointers
    DO i = 1, nfields
      IF (iscont(i, 1)) THEN
        dst_data_cptr(i) = C_LOC(recv(i)%p(1,1,1))
      ELSE
        nblk = SIZE(recv(i)%p, 3)
        ofs = cpy_psum + 1
        nl = nlev(i, 1)
        incr = nproma * nl * nblk
        last = cpy_psum + incr
        cpy(1:nproma, 1:nl, 1:nblk) => cpy_buf(ofs:last)
        cpy_psum = cpy_psum + incr
        cpy(:, :, :) = recv(i)%p
        dst_data_cptr(i) = C_LOC(cpy(1,1,1))
      END IF
    END DO
    IF (lsend) THEN
      DO i = 1, nfields
        IF (iscont(i, 2)) THEN
          src_data_cptr(i) = C_LOC(send(i)%p(1,1,1))
        ELSE
          nblk = SIZE(send(i)%p, 3)
          ofs = cpy_psum + 1
          nl = nlev(i, 2)
          incr = nproma * nl * nblk
          last = cpy_psum + incr
          cpy(1:nproma, 1:nl, 1:nblk) => cpy_buf(ofs:last)
          cpy_psum = cpy_psum + incr
          cpy(:, :, :) = send(i)%p
          src_data_cptr(i) = C_LOC(cpy(1,1,1))
        END IF
      END DO
    ELSE IF (cpy_recv) THEN
      DO i = 1, nfields
        nblk = SIZE(recv(i)%p, 3)
        ofs = cpy_psum + 1
        nl = nlev(i, 1)
        incr = nproma * nl * nblk
        last = cpy_psum + incr
        cpy(1:nproma, 1:nl, 1:nblk) => cpy_buf(ofs:last)
        cpy_psum = cpy_psum + incr
        cpy(:, :, :) = recv(i)%p
        src_data_cptr(i) = C_LOC(cpy(1,1,1))
      END DO
    ELSE
      src_data_cptr = dst_data_cptr
    END IF

    ! Reset kshift to 0 if 2D fields are passed together with 3D fields
    IF (PRESENT(nshift)) THEN
      DO i = 1, nfields
        kshift(i) = MERGE(0, nshift, nlev(i, 1) == 1)
      ENDDO
    ELSE
      kshift = 0
    ENDIF

    redist_coll = comm_pattern_get_redist(p_pat, nfields, &
         nlev, p_real_dp, kshift)

    CALL xt_redist_s_exchange(redist_coll, src_data_cptr, dst_data_cptr)

    cpy_psum = 0
    DO i = 1, nfields
      IF (.NOT. iscont(i, 1)) THEN
        nblk = SIZE(recv(i)%p, 3)
        ofs = cpy_psum + 1
        nl = nlev(i, 1)
        incr = nproma * nl * nblk
        last = cpy_psum + incr
        cpy(1:nproma, 1:nl, 1:nblk) => cpy_buf(ofs:last)
        cpy_psum = cpy_psum + incr
        recv(i)%p(:, :, :) = cpy
      END IF
    END DO

  END SUBROUTINE exchange_data_mult_dp_bottom

  !>
!! Does data exchange according to a communication pattern (in p_pat).
!!
!!
!! @par Revision History
!! Initial version by Rainer Johanni, Nov 2009
!! Optimized version by Guenther Zaengl to process 4D fields or up to seven 3D fields
!! in one step
!! yaxt version by Moritz Hanke, April 2016
!!
SUBROUTINE exchange_data_mult_sp(p_pat, ndim2tot, &
   recv, send, nshift)

   CLASS(t_comm_pattern_yaxt), INTENT(INOUT) :: p_pat

   TYPE(t_ptr_3d_sp), PTR_INTENT(in) :: recv(:)
   TYPE(t_ptr_3d_sp), OPTIONAL, PTR_INTENT(in) :: send(:)

   INTEGER, INTENT(IN)           :: ndim2tot
   INTEGER, OPTIONAL, INTENT(IN) :: nshift

   INTEGER :: nlev(SIZE(recv), 2)
   LOGICAL :: iscont(SIZE(recv), 2)
   INTEGER :: nproma1, nblk, nl, cpy_shape(3), cpy_size
   INTEGER :: i, nfields
   LOGICAL :: lsend, nproma_mismatch_found, cpy_recv

   CHARACTER(len=*), PARAMETER :: &
     routine = "mo_communication::exchange_data_mult_sp"

!-----------------------------------------------------------------------

   IF (itype_exch_barrier == 1 .OR. itype_exch_barrier == 3) THEN
     IF (activate_sync_timers) CALL timer_start(timer_barrier)
     CALL p_barrier(p_pat%comm)
     IF (activate_sync_timers) CALL timer_stop(timer_barrier)
   ENDIF

   start_sync_timer(timer_exch_data)

   lsend = PRESENT(send)
   cpy_recv = .NOT. (lsend .OR. p_pat%inplace)

   nfields = SIZE(recv)
   nproma1 = SIZE(recv(1)%p, 1)
   nproma_mismatch_found = .FALSE.
   cpy_size = 0
   DO i = 1, nfields
     iscont(i, 1) = IS_CONTIGUOUS(recv(i)%p)
     cpy_shape = SHAPE(recv(i)%p)
     nproma_mismatch_found = nproma_mismatch_found &
       &          .OR. (i > 1 .AND. cpy_shape(1) /= nproma1)
     nl = cpy_shape(2)
     nblk = cpy_shape(3)
     nlev(i, 1) = nl
     cpy_size =   cpy_size &
       &        + nproma1 * nl * nblk &
       &             * (MERGE(0, 1, iscont(i, 1)) + MERGE(1, 0, cpy_recv))
   END DO
   IF (lsend) THEN
     DO i = 1, nfields
       iscont(i, 2) = IS_CONTIGUOUS(send(i)%p)
       cpy_shape = SHAPE(send(i)%p)
       nproma_mismatch_found = nproma_mismatch_found &
         &        .OR. cpy_shape(1) /= nproma1
       nl = cpy_shape(2)
       nblk = cpy_shape(3)
       cpy_size =   cpy_size &
         &        + nproma1 * MERGE(0, nl, iscont(i, 2)) * nblk
       nlev(i, 2) = nl
     END DO
   ELSE
     iscont(:, 2) = .TRUE.
     nlev(:, 2) = nlev(:, 1)
   END IF

   IF (nproma_mismatch_found) &
     CALL finish(routine, "inconsistent array shapes detected")

   CALL exchange_data_mult_sp_bottom(p_pat, cpy_size, nlev, iscont, &
        recv, send, nshift)

   stop_sync_timer(timer_exch_data)

END SUBROUTINE exchange_data_mult_sp

  !>
  !! Does data exchange according to a communication pattern (in p_pat).
  !!
  !!
  !! @par Revision History
  !! Initial version by Rainer Johanni, Nov 2009
  !! Optimized version by Guenther Zaengl to process 4D fields or up to seven 3D fields
  !! in one step
  !! yaxt version by Moritz Hanke, April 2016
  !!
  SUBROUTINE exchange_data_mult_sp_bottom(p_pat, cpy_size, nlev, iscont, &
    recv, send, nshift)

    CLASS(t_comm_pattern_yaxt), INTENT(INOUT) :: p_pat

    TYPE(t_ptr_3d_sp), INTENT(in) :: recv(:)
    TYPE(t_ptr_3d_sp), OPTIONAL, INTENT(in) :: send(:)

    INTEGER, INTENT(IN) :: cpy_size, nlev(SIZE(recv), 2)
    LOGICAL, INTENT(IN)           :: iscont(SIZE(recv), 2)
    INTEGER, OPTIONAL, INTENT(IN) :: nshift

    INTEGER :: i, nfields, cpy_psum, ofs, last, nproma, nblk, incr, nl, &
         kshift(SIZE(recv))
    LOGICAL :: lsend, cpy_recv
    REAL(sp), TARGET :: cpy_buf(cpy_size)
    REAL(sp), POINTER :: cpy(:, :, :)
    TYPE(c_ptr) :: src_data_cptr(SIZE(recv)), dst_data_cptr(SIZE(recv))
    TYPE(xt_redist) :: redist_coll

    lsend = PRESENT(send)
    cpy_recv = .NOT. (lsend .OR. p_pat%inplace)

    nfields = SIZE(recv)
    cpy_psum = 0
    nproma = SIZE(recv(1)%p, 1)
    ! set up C pointers
    DO i = 1, nfields
      IF (iscont(i, 1)) THEN
        dst_data_cptr(i) = C_LOC(recv(i)%p(1,1,1))
      ELSE
        nblk = SIZE(recv(i)%p, 3)
        ofs = cpy_psum + 1
        nl = nlev(i, 1)
        incr = nproma * nl * nblk
        last = cpy_psum + incr
        cpy(1:nproma, 1:nl, 1:nblk) => cpy_buf(ofs:last)
        cpy_psum = cpy_psum + incr
        cpy(:, :, :) = recv(i)%p
        dst_data_cptr(i) = C_LOC(cpy(1,1,1))
      END IF
    END DO
    IF (lsend) THEN
      DO i = 1, nfields
        IF (iscont(i, 2)) THEN
          src_data_cptr(i) = C_LOC(send(i)%p(1,1,1))
        ELSE
          nblk = SIZE(send(i)%p, 3)
          ofs = cpy_psum + 1
          nl = nlev(i, 2)
          incr = nproma * nl * nblk
          last = cpy_psum + incr
          cpy(1:nproma, 1:nl, 1:nblk) => cpy_buf(ofs:last)
          cpy_psum = cpy_psum + incr
          cpy(:, :, :) = send(i)%p
          src_data_cptr(i) = C_LOC(cpy(1,1,1))
        END IF
      END DO
    ELSE IF (cpy_recv) THEN
      DO i = 1, nfields
        nblk = SIZE(recv(i)%p, 3)
        ofs = cpy_psum + 1
        nl = nlev(i, 1)
        incr = nproma * nl * nblk
        last = cpy_psum + incr
        cpy(1:nproma, 1:nl, 1:nblk) => cpy_buf(ofs:last)
        cpy_psum = cpy_psum + incr
        cpy(:, :, :) = recv(i)%p
        src_data_cptr(i) = C_LOC(cpy(1,1,1))
      END DO
    ELSE
      src_data_cptr = dst_data_cptr
    END IF

    ! Reset kshift to 0 if 2D fields are passed together with 3D fields
    IF (PRESENT(nshift)) THEN
      DO i = 1, nfields
        kshift(i) = MERGE(0, nshift, nlev(i, 1) == 1)
      ENDDO
    ELSE
      kshift = 0
    ENDIF

    redist_coll = comm_pattern_get_redist(p_pat, nfields, &
         nlev, p_real_sp, kshift)

    CALL xt_redist_s_exchange(redist_coll, src_data_cptr, dst_data_cptr)

    cpy_psum = 0
    DO i = 1, nfields
      IF (.NOT. iscont(i, 1)) THEN
        nblk = SIZE(recv(i)%p, 3)
        ofs = cpy_psum + 1
        nl = nlev(i, 1)
        incr = nproma * nl * nblk
        last = cpy_psum + incr
        cpy(1:nproma, 1:nl, 1:nblk) => cpy_buf(ofs:last)
        cpy_psum = cpy_psum + incr
        recv(i)%p(:, :, :) = cpy
      END IF
    END DO

  END SUBROUTINE exchange_data_mult_sp_bottom

!>
!! Does data exchange according to a communication pattern (in p_pat).
!!
!!
!! @par Revision History
!! Initial version by Rainer Johanni, Nov 2009
!! Optimized version by Guenther Zaengl to process 4D fields and 3D fields with either single
!! precision or double precision
!! yaxt version by Moritz Hanke, April 2016
!!
SUBROUTINE exchange_data_mult_mixprec(p_pat, nfields_dp, ndim2tot_dp, &
     nfields_sp, ndim2tot_sp, recv_dp, send_dp, recv_sp, send_sp, nshift)

  CLASS(t_comm_pattern_yaxt), INTENT(INOUT) :: p_pat

    TYPE(t_ptr_3d), PTR_INTENT(in), OPTIONAL :: recv_dp(:)
    TYPE(t_ptr_3d), PTR_INTENT(in), OPTIONAL :: send_dp(:)
    TYPE(t_ptr_3d_sp), PTR_INTENT(in), OPTIONAL :: recv_sp(:)
    TYPE(t_ptr_3d_sp), PTR_INTENT(in), OPTIONAL :: send_sp(:)

    INTEGER, INTENT(IN)           :: nfields_dp, ndim2tot_dp, &
      nfields_sp, ndim2tot_sp
  INTEGER, OPTIONAL, INTENT(IN) :: nshift

  IF (nfields_dp > 0) THEN
    CALL exchange_data_mult_dp(p_pat=p_pat, &
      ndim2tot=ndim2tot_dp, recv=recv_dp, send=send_dp, nshift=nshift)
  END IF

  IF (nfields_sp > 0) THEN
    CALL exchange_data_mult_sp(p_pat=p_pat, &
      ndim2tot=ndim2tot_sp, recv=recv_sp, send=send_sp, nshift=nshift)
  END IF

END SUBROUTINE exchange_data_mult_mixprec

!>
!! Does data exchange according to a communication pattern (in p_pat).
!!
!!
!! @par Revision History
!! Initial version by Rainer Johanni, Nov 2009
!! Optimized version by Guenther Zaengl to process a 4D field whose extra dimension
!! is on the first index
!!
SUBROUTINE exchange_data_4de1(p_pat, nfields, ndim2tot, recv, send)

   CLASS(t_comm_pattern_yaxt), INTENT(INOUT) :: p_pat

   REAL(dp), INTENT(INOUT)           :: recv(:,:,:,:)
   REAL(dp), INTENT(IN   ), OPTIONAL :: send(:,:,:,:)

   INTEGER, INTENT(IN)           :: nfields, ndim2tot

   INTEGER :: data_type
   TYPE(xt_redist) :: redist

   INTEGER :: nlev(1, 2), i, comm_size
   LOGICAL, SAVE :: first_call = .TRUE., mpi_is_buggy

   CHARACTER(len=*), PARAMETER :: &
     routine = modname//"::exchange_data_4de1"

!-----------------------------------------------------------------------

   comm_size = p_comm_size(p_pat%comm)
   IF (comm_size > 1 .AND. first_call) THEN

     first_call = .FALSE.
     mpi_is_buggy = test_mpich_bug(p_pat%comm)
     IF (mpi_is_buggy) THEN
       CALL message(routine, &
            "WARNING: your MPI contains a serious bug, in order to &
            &avoid problems a workaround has been activated in &
            &exchange_data_4de1 that significatly decreases its performance")
     END IF
   END IF

   IF (mpi_is_buggy) THEN
     ! MPICH seems to have a problem with that datatype generated by
     ! exchange_data_4de1. This affects intelmpi and mvapich.
     ! This issue is currently under investigation. Until it is fixed, we need
     ! the following workaround.
     IF (PRESENT(send)) THEN
       DO i = 1, SIZE(recv, 1)
           CALL exchange_data_r3d(p_pat, recv(i,:,:,:), send(i,:,:,:))
       END DO
     ELSE
       DO i = 1, SIZE(recv, 1)
           CALL exchange_data_r3d(p_pat, recv(i,:,:,:))
       END DO
     END IF
     RETURN
   END IF

   IF (itype_exch_barrier == 1 .OR. itype_exch_barrier == 3) THEN
     IF (activate_sync_timers) CALL timer_start(timer_barrier)
     CALL p_barrier(p_pat%comm)
     IF (activate_sync_timers) CALL timer_stop(timer_barrier)
   ENDIF

   start_sync_timer(timer_exch_data)

   stop_sync_timer(timer_exch_data)

   data_type = comm_pattern_get_contiguous_data_type(p_pat, p_real_dp, nfields)

   nlev(1, 1) = SIZE(recv,3)
   IF (PRESENT(send)) THEN
     nlev(1, 2) = SIZE(send,3)
   ELSE
     nlev(1, 2) = nlev(1, 1)
   END IF
   redist = comm_pattern_get_redist(p_pat, 1, nlev, data_type)

   IF (PRESENT(send)) THEN
     CALL xt_redist_s_exchange1_contiguous(redist, send, recv)
   ELSE
     IF (p_pat%inplace) THEN
       CALL xt_redist_s_exchange1_contiguous_inplace(redist, recv)
     ELSE
       CALL xt_redist_s_exchange1_contiguous_copy(redist, recv, SIZE(recv))
     END IF
   ENDIF

CONTAINS

  ! this wrapper is needed because we have to ensure that the input array are in
  ! contiguous memory (the keyword CONTIGUOUS is Fortran2008, which is not
  ! required by ICON)

  SUBROUTINE xt_redist_s_exchange1_contiguous(redist, send, recv)

    TYPE(xt_redist), INTENT(IN), TARGET :: redist
    REAL(dp), INTENT(INOUT), TARGET :: recv(*)
    REAL(dp), INTENT(IN), TARGET :: send(*)

    CALL xt_redist_s_exchange1(redist, c_loc(send), c_loc(recv))

  END SUBROUTINE

  SUBROUTINE xt_redist_s_exchange1_contiguous_inplace(redist, recv)

    TYPE(xt_redist), INTENT(IN), TARGET :: redist
    REAL(dp), INTENT(INOUT), TARGET :: recv(*)

    CALL xt_redist_s_exchange1(redist, c_loc(recv), c_loc(recv))

  END SUBROUTINE

    SUBROUTINE xt_redist_s_exchange1_contiguous_copy(redist, recv, n)

      TYPE(xt_redist), INTENT(IN) :: redist
      INTEGER, INTENT(in) :: n
      REAL(dp), INTENT(INOUT), TARGET :: recv(n)
      REAL(dp), TARGET :: send(n)

      send = recv
      CALL xt_redist_s_exchange1(redist, c_loc(send), c_loc(recv))

    END SUBROUTINE xt_redist_s_exchange1_contiguous_copy

END SUBROUTINE exchange_data_4de1

!>
!! Does data exchange according to a communication pattern (in p_pat).
!!
!!
!! @par Revision History
!! Initial version by Rainer Johanni, Nov 2009
!! Optimized version by Guenther Zaengl to process up to two 4D fields or up to six 3D fields
!! for an array-sized communication pattern (as needed for boundary interpolation) in one step
!! yaxt version by Moritz Hanke, April 2016
!!
SUBROUTINE exchange_data_grf(p_pat_coll, nfields, ndim2tot, recv1, send1, &
                             recv2, send2, recv3, send3, recv4, send4, &
                             recv5, send5, recv6, send6, recv4d1, send4d1, &
                             recv4d2, send4d2)

   CLASS(t_comm_pattern_collection_yaxt), TARGET, INTENT(INOUT) :: p_pat_coll

   REAL(dp), INTENT(INOUT), TARGET, OPTIONAL ::  &
     recv1(:,:,:), recv2(:,:,:), recv3(:,:,:), recv4d1(:,:,:,:), &
     recv4(:,:,:), recv5(:,:,:), recv6(:,:,:), recv4d2(:,:,:,:)
   ! Note: the last index of the send fields corresponds to the dimension of p_pat
   ! On the other hand, they are not blocked and have the vertical index first
   REAL(dp), INTENT(IN   ), TARGET, OPTIONAL ::  &
     send1(:,:,:), send2(:,:,:), send3(:,:,:), send4d1(:,:,:,:), &
     send4(:,:,:), send5(:,:,:), send6(:,:,:), send4d2(:,:,:,:)

   INTEGER, INTENT(IN)           :: nfields  ! total number of input fields
   INTEGER, INTENT(IN)           :: ndim2tot ! sum of vertical levels of input fields

   CHARACTER(len=*), PARAMETER :: routine = "mo_communication::exchange_data_grf"

   INTEGER :: npats, dst_nlev(nfields), src_nlev(nfields), &
    &         src_fsize4d(nfields), dst_fsize4d(nfields)
   TYPE(xt_redist) :: redist_coll
   INTEGER, TARGET, SAVE :: dummy

!-----------------------------------------------------------------------

   IF (itype_exch_barrier == 1 .OR. itype_exch_barrier == 3) THEN
     IF (activate_sync_timers) CALL timer_start(timer_barrier)
     CALL p_barrier(p_pat_coll%patterns(1)%p%comm)
     IF (activate_sync_timers) CALL timer_stop(timer_barrier)
   ENDIF

   start_sync_timer(timer_exch_data)

!-----------------------------------------------------------------------

   ! Set pointers to input fields
   IF (PRESENT(recv4d1) .AND. .NOT. PRESENT(recv4d2)) THEN
     dst_nlev = SIZE(recv4d1, 2)
     src_nlev = SIZE(send4d1, 1)
   ELSE IF (PRESENT(recv4d1) .AND. PRESENT(recv4d2)) THEN
    dst_nlev(1:nfields/2) = SIZE(recv4d1, 2)
    dst_nlev(1+nfields/2:) = SIZE(recv4d2, 2)
    src_nlev(1:nfields/2) = SIZE(send4d1, 1)
    src_nlev(1+nfields/2:) = SIZE(send4d2, 1)
   ELSE
     dst_nlev = 0
     src_nlev = 0
     src_fsize4d = 0
     dst_fsize4d = 0
     IF (PRESENT(recv1)) THEN
       dst_nlev(1) = SIZE(recv1, 2)
       src_nlev(1) = SIZE(send1, 1)
       src_fsize4d(1) = SIZE(send1, 1) * SIZE(send1, 2)
       dst_fsize4d(1) = SIZE(recv1)
       IF (PRESENT(recv2)) THEN
         dst_nlev(2) = SIZE(recv2, 2)
         src_nlev(2) = SIZE(send2, 1)
         src_fsize4d(2) = SIZE(send2, 1) * SIZE(send2, 2)
         dst_fsize4d(2) = SIZE(recv2)
         IF (PRESENT(recv3)) THEN
           dst_nlev(3) = SIZE(recv3, 2)
           src_nlev(3) = SIZE(send3, 1)
           src_fsize4d(3) = SIZE(send3, 1) * SIZE(send3, 2)
           dst_fsize4d(3) = SIZE(recv3)
           IF (PRESENT(recv4)) THEN
             dst_nlev(4) = SIZE(recv4, 2)
             src_nlev(4) = SIZE(send4, 1)
             src_fsize4d(4) = SIZE(send4, 1) * SIZE(send4, 2)
             dst_fsize4d(4) = SIZE(recv4)
             IF (PRESENT(recv5)) THEN
               dst_nlev(5) = SIZE(recv5, 2)
               src_nlev(5) = SIZE(send5, 1)
               src_fsize4d(5) = SIZE(send5, 1) * SIZE(send5, 2)
               dst_fsize4d(5) = SIZE(recv5)
               IF (PRESENT(recv6)) THEN
                 dst_nlev(6) = SIZE(recv6, 2)
                 src_nlev(6) = SIZE(send6, 1)
                 src_fsize4d(6) = SIZE(send6, 1) * SIZE(send6, 2)
                 dst_fsize4d(6) = SIZE(recv6)
               ENDIF
             ENDIF
           ENDIF
         ENDIF
       ENDIF
     ENDIF
   ENDIF

   stop_sync_timer(timer_exch_data)

   redist_coll = comm_pattern_collection_get_redist(p_pat_coll, nfields, &
    &                                               dst_nlev, src_nlev)

   ! Number of communication patterns provided on input
   npats = SIZE(p_pat_coll%patterns)

   IF (PRESENT(recv4d1) .AND. .NOT. PRESENT(recv4d2)) THEN

     CALL xt_redist_s_exchange_grf_4d( &
       redist_coll, nfields, npats, SIZE(send4d1,1)*SIZE(send4d1,2), &
       SIZE(recv4d1,1)*SIZE(recv4d1,2)*SIZE(recv4d1,3), recv4d1, send4d1)

   ELSE IF (PRESENT(recv4d1) .AND. PRESENT(recv4d2)) THEN

     CALL xt_redist_s_exchange_grf_2x4d( &
       redist_coll, nfields, npats, &
       (/SIZE(send4d1,1)*SIZE(send4d1,2), SIZE(send4d2,1)*SIZE(send4d2,2)/), &
       (/SIZE(recv4d1,1)*SIZE(recv4d1,2)*SIZE(recv4d1,3), &
         SIZE(recv4d2,1)*SIZE(recv4d2,2)*SIZE(recv4d2,3)/), &
       recv4d1, send4d1, recv4d2, send4d2)

   ELSE

     SELECT CASE (nfields)
       CASE (1)
         CALL xt_redist_s_exchange_grf_1x3d(redist_coll, npats, src_fsize4d(1), &
          dst_fsize4d(1), recv1, send1)
       CASE (2)
         CALL xt_redist_s_exchange_grf_2x3d(redist_coll, npats, src_fsize4d, &
          dst_fsize4d, recv1, send1, recv2, send2)
       CASE (3)
         CALL xt_redist_s_exchange_grf_3x3d(redist_coll, npats, src_fsize4d, &
          dst_fsize4d, recv1, send1, recv2, send2, recv3, send3)
       CASE (4)
         CALL xt_redist_s_exchange_grf_4x3d(redist_coll, npats, src_fsize4d, &
          dst_fsize4d, recv1, send1, recv2, send2, recv3, send3, recv4, send4)
       CASE (5)
         CALL xt_redist_s_exchange_grf_5x3d(redist_coll, npats, src_fsize4d, &
          dst_fsize4d, recv1, send1, recv2, send2, recv3, send3, recv4, send4, &
          recv5, send5)
       CASE (6)
         CALL xt_redist_s_exchange_grf_6x3d(redist_coll, npats, src_fsize4d, &
          dst_fsize4d, recv1, send1, recv2, send2, recv3, send3, recv4, send4, &
          recv5, send5, recv6, send6)
     END SELECT
   ENDIF

CONTAINS

  !> grf data is often empty (not present) on tasks not involved.
  !! Since empty slices must not be indexed into, wrap the decision
  !! here and return dummy pointers in case the slice is empty.
  SUBROUTINE get_grf_src_data_cptr(data_cptr, var, npats, fsize4d)
    INTEGER, INTENT(in) :: npats, fsize4d
    TYPE(c_ptr), INTENT(out) :: data_cptr(npats)
    REAL(dp), INTENT(IN), TARGET ::  var(*)
    INTEGER :: i
    IF (fsize4d > 0) THEN
      DO i = 1, npats
        data_cptr(i) = C_LOC(var((i-1) * fsize4d + 1))
      END DO
    ELSE
      data_cptr = C_LOC(dummy)
    END IF
  END SUBROUTINE get_grf_src_data_cptr
  SUBROUTINE get_grf_dst_data_cptr(data_cptr, var, npats, nfields, fsize4d)
    INTEGER, INTENT(in) :: npats, fsize4d, nfields
    TYPE(c_ptr), INTENT(out) :: data_cptr(npats * nfields)
    REAL(dp), INTENT(IN), TARGET ::  var(*)
    INTEGER :: i
    IF (fsize4d > 0) THEN
      DO i = 1, nfields
        data_cptr((i-1)*npats+1:i*npats) = &
          c_loc(var((i-1) * fsize4d + 1))
      END DO
    ELSE
      data_cptr = C_LOC(dummy)
    END IF
  END SUBROUTINE get_grf_dst_data_cptr


  ! the following wrappers are needed because we have to ensure that the input
  ! arrays are in contiguous memory (the keyword CONTIGUOUS is Fortran2008,
  ! which is not required by ICON)

  SUBROUTINE xt_redist_s_exchange_grf_2x4d(redist_coll, nfields, npats, &
    src_fsize4d, dst_fsize4d, recv4d1, send4d1, recv4d2, send4d2)

    TYPE(xt_redist) :: redist_coll
    INTEGER, INTENT(IN) :: nfields, npats, src_fsize4d(2), dst_fsize4d(2)
    REAL(dp), INTENT(INOUT), TARGET :: recv4d1(*), recv4d2(*)
    REAL(dp), INTENT(IN), TARGET ::  send4d1(*), send4d2(*)

    INTEGER :: nfields2
    TYPE(c_ptr) :: src_data_cptr(nfields*npats), dst_data_cptr(nfields*npats)

    nfields2 = ((nfields)/2)
    CALL get_grf_src_data_cptr(src_data_cptr(1:nfields2*npats), send4d1, &
         nfields2*npats, src_fsize4d(1))
    CALL get_grf_src_data_cptr(src_data_cptr(nfields2*npats+1:nfields*npats), &
         send4d2, nfields2*npats, src_fsize4d(2))
    CALL get_grf_dst_data_cptr(dst_data_cptr(1:nfields2*npats), recv4d1, &
      &                        npats, nfields2, dst_fsize4d(1))
    CALL get_grf_dst_data_cptr(dst_data_cptr(nfields2*npats+1:), recv4d2, &
      &                        npats, nfields2, dst_fsize4d(2))

    CALL xt_redist_s_exchange(redist_coll, src_data_cptr, dst_data_cptr)

  END SUBROUTINE xt_redist_s_exchange_grf_2x4d

  SUBROUTINE xt_redist_s_exchange_grf_4d(redist_coll, nfields, npats, &
    src_fsize4d, dst_fsize4d, recv4d, send4d)

    TYPE(xt_redist) :: redist_coll
    INTEGER, INTENT(IN) :: nfields, npats, src_fsize4d, dst_fsize4d
    REAL(dp), INTENT(INOUT), TARGET :: recv4d(*)
    REAL(dp), INTENT(IN), TARGET ::  send4d(*)

    TYPE(c_ptr) :: src_data_cptr(nfields*npats), dst_data_cptr(nfields*npats)

    CALL get_grf_src_data_cptr(src_data_cptr, send4d, nfields*npats, &
         src_fsize4d)
    CALL get_grf_dst_data_cptr(dst_data_cptr, recv4d, npats, nfields, &
      &                        dst_fsize4d)

    CALL xt_redist_s_exchange(redist_coll, src_data_cptr, dst_data_cptr)

  END SUBROUTINE xt_redist_s_exchange_grf_4d

  SUBROUTINE xt_redist_s_exchange_grf_6x3d(redist_coll, npats, &
    src_fsize4d, dst_fsize4d, recv1, send1, recv2, send2, recv3, send3, &
    recv4, send4, recv5, send5, recv6, send6)

    TYPE(xt_redist) :: redist_coll
    INTEGER, INTENT(IN) :: npats, src_fsize4d(6), dst_fsize4d(6)
    REAL(dp), INTENT(INOUT), TARGET :: recv1(*), recv2(*), recv3(*), &
      &                                recv4(*), recv5(*), recv6(*)
    REAL(dp), INTENT(IN), TARGET ::  send1(*), send2(*), send3(*), &
      &                              send4(*), send5(*), send6(*)

    TYPE(c_ptr) :: src_data_cptr(6*npats), dst_data_cptr(6*npats)

    CALL get_grf_src_data_cptr(src_data_cptr(0*npats+1:1*npats), send1, npats, &
         src_fsize4d(1))
    CALL get_grf_src_data_cptr(src_data_cptr(1*npats+1:2*npats), send2, npats, &
         src_fsize4d(2))
    CALL get_grf_src_data_cptr(src_data_cptr(2*npats+1:3*npats), send3, npats, &
         src_fsize4d(3))
    CALL get_grf_src_data_cptr(src_data_cptr(3*npats+1:4*npats), send4, npats, &
         src_fsize4d(4))
    CALL get_grf_src_data_cptr(src_data_cptr(4*npats+1:5*npats), send5, npats, &
         src_fsize4d(5))
    CALL get_grf_src_data_cptr(src_data_cptr(5*npats+1:6*npats), send6, npats, &
         src_fsize4d(6))
    CALL get_grf_dst_data_cptr(dst_data_cptr(0*npats+1:1*npats), recv1, npats, &
      &                        1, dst_fsize4d(1))
    CALL get_grf_dst_data_cptr(dst_data_cptr(1*npats+1:2*npats), recv2, npats, &
      &                        1, dst_fsize4d(2))
    CALL get_grf_dst_data_cptr(dst_data_cptr(2*npats+1:3*npats), recv3, npats, &
      &                        1, dst_fsize4d(3))
    CALL get_grf_dst_data_cptr(dst_data_cptr(3*npats+1:4*npats), recv4, npats, &
      &                        1, dst_fsize4d(4))
    CALL get_grf_dst_data_cptr(dst_data_cptr(4*npats+1:5*npats), recv5, npats, &
      &                        1, dst_fsize4d(5))
    CALL get_grf_dst_data_cptr(dst_data_cptr(5*npats+1:6*npats), recv6, npats, &
      &                        1, dst_fsize4d(6))

    CALL xt_redist_s_exchange(redist_coll, src_data_cptr, dst_data_cptr)

  END SUBROUTINE xt_redist_s_exchange_grf_6x3d

  SUBROUTINE xt_redist_s_exchange_grf_5x3d(redist_coll, npats, &
    src_fsize4d, dst_fsize4d, recv1, send1, recv2, send2, recv3, send3, &
    recv4, send4, recv5, send5)

    TYPE(xt_redist) :: redist_coll
    INTEGER, INTENT(IN) :: npats, src_fsize4d(5), dst_fsize4d(5)
    REAL(dp), INTENT(INOUT), TARGET :: recv1(*), recv2(*), recv3(*), &
      &                                recv4(*), recv5(*)
    REAL(dp), INTENT(IN), TARGET ::  send1(*), send2(*), send3(*), &
      &                              send4(*), send5(*)

    TYPE(c_ptr) :: src_data_cptr(5*npats), dst_data_cptr(5*npats)

    CALL get_grf_src_data_cptr(src_data_cptr(0*npats+1:1*npats), send1, npats, &
         src_fsize4d(1))
    CALL get_grf_src_data_cptr(src_data_cptr(1*npats+1:2*npats), send2, npats, &
         src_fsize4d(2))
    CALL get_grf_src_data_cptr(src_data_cptr(2*npats+1:3*npats), send3, npats, &
         src_fsize4d(3))
    CALL get_grf_src_data_cptr(src_data_cptr(3*npats+1:4*npats), send4, npats, &
         src_fsize4d(4))
    CALL get_grf_src_data_cptr(src_data_cptr(4*npats+1:5*npats), send5, npats, &
         src_fsize4d(5))
    CALL get_grf_dst_data_cptr(dst_data_cptr(0*npats+1:1*npats), recv1, npats, &
      &                        1, dst_fsize4d(1))
    CALL get_grf_dst_data_cptr(dst_data_cptr(1*npats+1:2*npats), recv2, npats, &
      &                        1, dst_fsize4d(2))
    CALL get_grf_dst_data_cptr(dst_data_cptr(2*npats+1:3*npats), recv3, npats, &
      &                        1, dst_fsize4d(3))
    CALL get_grf_dst_data_cptr(dst_data_cptr(3*npats+1:4*npats), recv4, npats, &
      &                        1, dst_fsize4d(4))
    CALL get_grf_dst_data_cptr(dst_data_cptr(4*npats+1:5*npats), recv5, npats, &
      &                        1, dst_fsize4d(5))

    CALL xt_redist_s_exchange(redist_coll, src_data_cptr, dst_data_cptr)

  END SUBROUTINE xt_redist_s_exchange_grf_5x3d

  SUBROUTINE xt_redist_s_exchange_grf_4x3d(redist_coll, npats, &
    src_fsize4d, dst_fsize4d, recv1, send1, recv2, send2, recv3, send3, &
    recv4, send4)

    TYPE(xt_redist) :: redist_coll
    INTEGER, INTENT(IN) :: npats, src_fsize4d(4), dst_fsize4d(4)
    REAL(dp), INTENT(INOUT), TARGET :: recv1(*), recv2(*), recv3(*), recv4(*)
    REAL(dp), INTENT(IN), TARGET ::  send1(*), send2(*), send3(*), send4(*)

    TYPE(c_ptr) :: src_data_cptr(4*npats), dst_data_cptr(4*npats)

    CALL get_grf_src_data_cptr(src_data_cptr(1:npats), send1, npats, src_fsize4d(1))
    CALL get_grf_src_data_cptr(src_data_cptr(1*npats+1:2*npats), send2, npats, &
         src_fsize4d(2))
    CALL get_grf_src_data_cptr(src_data_cptr(2*npats+1:3*npats), send3, npats, &
         src_fsize4d(3))
    CALL get_grf_src_data_cptr(src_data_cptr(3*npats+1:4*npats), send4, npats, &
         src_fsize4d(4))
    CALL get_grf_dst_data_cptr(dst_data_cptr(0*npats+1:1*npats), recv1, npats, &
      &                        1, dst_fsize4d(1))
    CALL get_grf_dst_data_cptr(dst_data_cptr(1*npats+1:2*npats), recv2, npats, &
      &                        1, dst_fsize4d(2))
    CALL get_grf_dst_data_cptr(dst_data_cptr(2*npats+1:3*npats), recv3, npats, &
      &                        1, dst_fsize4d(3))
    CALL get_grf_dst_data_cptr(dst_data_cptr(3*npats+1:4*npats), recv4, npats, &
      &                        1, dst_fsize4d(4))

    CALL xt_redist_s_exchange(redist_coll, src_data_cptr, dst_data_cptr)

  END SUBROUTINE xt_redist_s_exchange_grf_4x3d

  SUBROUTINE xt_redist_s_exchange_grf_3x3d(redist_coll, npats, &
    src_fsize4d, dst_fsize4d, recv1, send1, recv2, send2, recv3, send3)

    TYPE(xt_redist) :: redist_coll
    INTEGER, INTENT(IN) :: npats, src_fsize4d(3), dst_fsize4d(3)
    REAL(dp), INTENT(INOUT), TARGET :: recv1(*), recv2(*), recv3(*)
    REAL(dp), INTENT(IN), TARGET ::  send1(*), send2(*), send3(*)

    TYPE(c_ptr) :: src_data_cptr(3*npats), dst_data_cptr(3*npats)

    CALL get_grf_src_data_cptr(src_data_cptr(1:npats), send1, npats, src_fsize4d(1))
    CALL get_grf_src_data_cptr(src_data_cptr(1*npats+1:2*npats), send2, npats, &
         src_fsize4d(2))
    CALL get_grf_src_data_cptr(src_data_cptr(2*npats+1:3*npats), send3, npats, &
         src_fsize4d(3))
    CALL get_grf_dst_data_cptr(dst_data_cptr(0*npats+1:1*npats), recv1, npats, &
      &                        1, dst_fsize4d(1))
    CALL get_grf_dst_data_cptr(dst_data_cptr(1*npats+1:2*npats), recv2, npats, &
      &                        1, dst_fsize4d(2))
    CALL get_grf_dst_data_cptr(dst_data_cptr(2*npats+1:3*npats), recv3, npats, &
      &                        1, dst_fsize4d(3))

    CALL xt_redist_s_exchange(redist_coll, src_data_cptr, dst_data_cptr)

  END SUBROUTINE xt_redist_s_exchange_grf_3x3d

  SUBROUTINE xt_redist_s_exchange_grf_2x3d(redist_coll, npats, &
    src_fsize4d, dst_fsize4d, recv1, send1, recv2, send2)

    TYPE(xt_redist) :: redist_coll
    INTEGER, INTENT(IN) :: npats, src_fsize4d(2), dst_fsize4d(2)
    REAL(dp), INTENT(INOUT), TARGET :: recv1(*), recv2(*)
    REAL(dp), INTENT(IN), TARGET ::  send1(*), send2(*)

    TYPE(c_ptr) :: src_data_cptr(2*npats), dst_data_cptr(2*npats)

    CALL get_grf_src_data_cptr(src_data_cptr(1:npats), send1, npats, src_fsize4d(1))
    CALL get_grf_src_data_cptr(src_data_cptr(1*npats+1:2*npats), send2, npats, &
         src_fsize4d(2))
    CALL get_grf_dst_data_cptr(dst_data_cptr(0*npats+1:1*npats), recv1, npats, &
      &                        1, dst_fsize4d(1))
    CALL get_grf_dst_data_cptr(dst_data_cptr(1*npats+1:2*npats), recv2, npats, &
      &                        1, dst_fsize4d(2))

    CALL xt_redist_s_exchange(redist_coll, src_data_cptr, dst_data_cptr)

  END SUBROUTINE xt_redist_s_exchange_grf_2x3d

  SUBROUTINE xt_redist_s_exchange_grf_1x3d(redist_coll, npats, &
    src_fsize4d, dst_fsize4d, recv1, send1)

    TYPE(xt_redist) :: redist_coll
    INTEGER, INTENT(IN) :: npats, src_fsize4d, dst_fsize4d
    REAL(dp), INTENT(INOUT), TARGET :: recv1(*)
    REAL(dp), INTENT(IN), TARGET ::  send1(*)

    TYPE(c_ptr) :: src_data_cptr(npats), dst_data_cptr(npats)

    CALL get_grf_src_data_cptr(src_data_cptr(1:npats), send1, npats, src_fsize4d)
    CALL get_grf_dst_data_cptr(dst_data_cptr(0*npats+1:1*npats), recv1, npats, &
      &                        1, dst_fsize4d)

    CALL xt_redist_s_exchange(redist_coll, src_data_cptr, dst_data_cptr)

  END SUBROUTINE xt_redist_s_exchange_grf_1x3d

END SUBROUTINE exchange_data_grf

!-------------------------------------------------------------------------
!
!

!>
!! Interface for 2D arrays for exchange_data.
!!
!! Just reshapes the arrays and calls exchange_data.
!!
!! @par Revision History
!! Initial version by Rainer Johanni, Nov 2009
!! yaxt version by Moritz Hanke, April 2016
!!
!================================================================================================
! REAL SECTION ----------------------------------------------------------------------------------
!
SUBROUTINE exchange_data_r2d(p_pat, recv, send, add, l_recv_exists)
   !
   CLASS(t_comm_pattern_yaxt), INTENT(INOUT), TARGET :: p_pat
   REAL(dp), INTENT(INOUT), TARGET        :: recv(:,:)
   REAL(dp), INTENT(IN), OPTIONAL, TARGET :: send(:,:)
   REAL(dp), INTENT(IN), OPTIONAL, TARGET :: add (:,:)
   LOGICAL, OPTIONAL :: l_recv_exists

   CHARACTER(len=*), PARAMETER :: routine = "mo_communication::exchange_data_r2d"
   REAL(dp) :: tmp_recv(SIZE(recv,1),1,SIZE(recv,2))

   !-----------------------------------------------------------------------

   IF (PRESENT(send) .AND. .NOT. PRESENT(l_recv_exists)) THEN
     tmp_recv(:,1,:) = 0._dp
   ELSE
     tmp_recv(:,1,:) = recv(:,:)
   ENDIF

   IF (PRESENT(send)) THEN
      IF (PRESENT(add)) THEN
         CALL exchange_data_r3d(p_pat, tmp_recv, &
           send=RESHAPE(send,(/SIZE(send,1),1,SIZE(send,2)/)), &
           add =RESHAPE(add, (/SIZE(add ,1),1,SIZE(add ,2)/)))
      ELSE
         CALL exchange_data_r3d(p_pat, tmp_recv, &
           send=RESHAPE(send,(/SIZE(send,1),1,SIZE(send,2)/)))
      ENDIF
   ELSE
      IF (PRESENT(add)) THEN
         CALL exchange_data_r3d(p_pat, tmp_recv, &
           add =RESHAPE(add, (/SIZE(add ,1),1,SIZE(add ,2)/)))
      ELSE
         CALL exchange_data_r3d(p_pat, tmp_recv)
      ENDIF
   ENDIF

   recv(:,:) = tmp_recv(:,1,:)

END SUBROUTINE exchange_data_r2d

SUBROUTINE exchange_data_s2d(p_pat, recv, send, add, l_recv_exists)
   !
   CLASS(t_comm_pattern_yaxt), INTENT(INOUT), TARGET :: p_pat
   REAL(sp), INTENT(INOUT), TARGET        :: recv(:,:)
   REAL(sp), INTENT(IN), OPTIONAL, TARGET :: send(:,:)
   REAL(sp), INTENT(IN), OPTIONAL, TARGET :: add (:,:)
   LOGICAL, OPTIONAL :: l_recv_exists

   CHARACTER(len=*), PARAMETER :: routine = "mo_communication::exchange_data_s2d"
   REAL(sp) :: tmp_recv(SIZE(recv,1),1,SIZE(recv,2))

   !-----------------------------------------------------------------------

   IF (PRESENT(send) .AND. .NOT. PRESENT(l_recv_exists)) THEN
     tmp_recv(:,1,:) = 0._sp
   ELSE
     tmp_recv(:,1,:) = recv(:,:)
   ENDIF

   IF (PRESENT(send)) THEN
      IF (PRESENT(add)) THEN
         CALL exchange_data_s3d(p_pat, tmp_recv, &
           send=RESHAPE(send,(/SIZE(send,1),1,SIZE(send,2)/)), &
           add =RESHAPE(add, (/SIZE(add ,1),1,SIZE(add ,2)/)))
      ELSE
         CALL exchange_data_s3d(p_pat, tmp_recv, &
           send=RESHAPE(send,(/SIZE(send,1),1,SIZE(send,2)/)))
      ENDIF
   ELSE
      IF (PRESENT(add)) THEN
         CALL exchange_data_s3d(p_pat, tmp_recv, &
           add =RESHAPE(add, (/SIZE(add ,1),1,SIZE(add ,2)/)))
      ELSE
         CALL exchange_data_s3d(p_pat, tmp_recv)
      ENDIF
   ENDIF

   recv(:,:) = tmp_recv(:,1,:)

END SUBROUTINE exchange_data_s2d

!================================================================================================
! INTEGER SECTION -------------------------------------------------------------------------------
!
SUBROUTINE exchange_data_i2d(p_pat, recv, send, add, l_recv_exists)
   !
   CLASS(t_comm_pattern_yaxt), INTENT(INOUT), TARGET :: p_pat
   INTEGER, INTENT(INOUT), TARGET        :: recv(:,:)
   INTEGER, INTENT(IN), OPTIONAL, TARGET :: send(:,:)
   INTEGER, INTENT(IN), OPTIONAL, TARGET :: add (:,:)
   LOGICAL, OPTIONAL :: l_recv_exists

   CHARACTER(len=*), PARAMETER :: routine = "mo_communication::exchange_data_i2d"
   INTEGER :: tmp_recv(SIZE(recv,1),1,SIZE(recv,2))

   !-----------------------------------------------------------------------

   IF (PRESENT(send) .AND. .NOT. PRESENT(l_recv_exists)) THEN
     tmp_recv(:,1,:) = 0
   ELSE
     tmp_recv(:,1,:) = recv(:,:)
   ENDIF

   IF (PRESENT(send)) THEN
      IF (PRESENT(add)) THEN
         CALL exchange_data_i3d(p_pat, tmp_recv, &
           send=RESHAPE(send,(/SIZE(send,1),1,SIZE(send,2)/)), &
           add =RESHAPE(add, (/SIZE(add ,1),1,SIZE(add ,2)/)))
      ELSE
         CALL exchange_data_i3d(p_pat, tmp_recv, &
           send=RESHAPE(send,(/SIZE(send,1),1,SIZE(send,2)/)))
      ENDIF
   ELSE
      IF (PRESENT(add)) THEN
         CALL exchange_data_i3d(p_pat, tmp_recv, &
           add =RESHAPE(add, (/SIZE(add ,1),1,SIZE(add ,2)/)))
      ELSE
         CALL exchange_data_i3d(p_pat, tmp_recv)
      ENDIF
   ENDIF

   recv(:,:) = tmp_recv(:,1,:)

END SUBROUTINE exchange_data_i2d

  SUBROUTINE exchange_data_l2d(p_pat, recv, send, l_recv_exists)
    !
    CLASS(t_comm_pattern_yaxt), INTENT(INOUT), TARGET :: p_pat
    LOGICAL, INTENT(INOUT), TARGET        :: recv(:,:)
    LOGICAL, INTENT(IN), OPTIONAL, TARGET :: send(:,:)
    LOGICAL, OPTIONAL :: l_recv_exists

    CHARACTER(len=*), PARAMETER :: routine = "mo_communication::exchange_data_i2d"
    LOGICAL :: tmp_recv(SIZE(recv,1),1,SIZE(recv,2))

    !-----------------------------------------------------------------------

    IF (PRESENT(send) .AND. .NOT. PRESENT(l_recv_exists)) THEN
      tmp_recv(:,1,:) = .FALSE.
    ELSE
      tmp_recv(:,1,:) = recv(:,:)
    ENDIF

    IF (PRESENT(send)) THEN
      CALL exchange_data_l3d(p_pat, tmp_recv, &
           send=RESHAPE(send,(/SIZE(send,1),1,SIZE(send,2)/)))
    ELSE
      CALL exchange_data_l3d(p_pat, tmp_recv)
    ENDIF

    recv(:,:) = tmp_recv(:,1,:)

  END SUBROUTINE exchange_data_l2d

FUNCTION get_np_recv(comm_pat)
  CLASS(t_comm_pattern_yaxt), INTENT(IN) :: comm_pat
  INTEGER :: get_np_recv

  get_np_recv = xt_xmap_get_num_sources(comm_pat%xmap)
END FUNCTION get_np_recv

FUNCTION get_np_send(comm_pat)
  CLASS(t_comm_pattern_yaxt), INTENT(IN) :: comm_pat
  INTEGER :: get_np_send

  get_np_send = xt_xmap_get_num_destinations(comm_pat%xmap)
END FUNCTION get_np_send

SUBROUTINE get_pelist_recv(comm_pat, pelist_recv)
  CLASS(t_comm_pattern_yaxt), INTENT(IN) :: comm_pat
  INTEGER, INTENT(OUT) :: pelist_recv(:)

  IF (SIZE(pelist_recv) /= xt_xmap_get_num_sources(comm_pat%xmap)) &
    CALL finish("get_pelist_recv", &
      &         "pelist_recv size is different from np_recv")

  CALL xt_xmap_get_source_ranks(comm_pat%xmap, pelist_recv)
END SUBROUTINE get_pelist_recv

! This tests checks for a bug that can occur when using MPI datatypes.
! The following MPI implementations have been tested:
! buggy MPIs:
! intelmpi/2017.0.098
! intelmpi/4.1.3.049
! intelmpi/5.0.3.048
! intelmpi/5.1.0.038
! intelmpi/5.1.1.109
! intelmpi/5.1.2.150
! intelmpi/5.1.3.223
! mvapich2/1.9b
! mvapich2/2.1
!
! working MPIs:
! openmpi/1.10.1
! openmpi/1.10.2
! openmpi/1.8.4
!
! Description of the underlying problem:
! On rank 0 we have a 1D data array of size [NBLK * NIDX * VECTOR_SIZE * NLEV].
! If stored in a 3D array, the layout would be [NBLK][NLEV][NIDX][VECTOR_SIZE].
! The variables count, displs, lenghts describe points accessing the vectors of
! single level. Values for each of these points for each level is sent from
! rank 1 to rank 0.
! Each point consists of VECTOR_SIZE elements. The data sent by rank 1 is
! level-wise; data from different levels does not overlap in the send buffer.
!
FUNCTION test_mpich_bug(comm)
  INTEGER, INTENT(in) :: comm
  LOGICAL :: test_mpich_bug

  INTEGER, PARAMETER :: VECTOR_SIZE = 2, &
                        NLEV = 90, &
                        NIDX = 8, &
                        NBLK = 7, &
                        TAG = 0
  INTEGER :: displs(9)
  INTEGER :: lenghts(9)

  INTEGER :: comm_rank
  INTEGER :: num_exchange_indices, ierror, i, j, lev, vector
  INTEGER :: datatype_vector, datatype_1lev, datatype_with_extent, datatype_nlev
  INTEGER(KIND=MPI_ADDRESS_KIND) :: lb, extent, dummy
  REAL(dp) :: v
  REAL(dp), ALLOCATABLE :: ref_data(:), recv_data(:), send_data(:)

  displs = (/   1, & ! idx 1 blk 0
              720, & ! idx 0 blk 1
             1440, & ! idx 0 blk 2
             1447, & ! idx 7 blk 2
             2160, & ! idx 0 blk 3
             2166, & ! idx 6 blk 3
             2880, & ! idx 0 blk 4
             3600, & ! idx 0 blk 5
             4320/)  ! idx 0 blk 6
  lenghts = (/7, 8, 1, 1, 5, 2, 8, 6, 8/)

  ! compute the number of points that will be exchanged
  ! (not including VECTOR_SIZE and NLEV)
  num_exchange_indices = SUM(lenghts)

  comm_rank = p_comm_rank(comm)
  IF (comm_rank == 0) THEN

    ! create contiguous datatype for the vector
    CALL MPI_Type_contiguous(VECTOR_SIZE, MPI_DOUBLE_PRECISION, &
                             datatype_vector, ierror)

    ! generate datatype for all points of a single level
    CALL MPI_Type_indexed(SIZE(lenghts), lenghts, displs, datatype_vector, &
                          datatype_1lev, ierror)

    ! create a datatype for receiving all levels at once (each level is
    ! VECTOR_SIZE * NIDX elements appart)
    CALL MPI_Type_get_extent(datatype_1lev, lb, dummy, ierror)
    CALL MPI_Type_get_extent(MPI_DOUBLE_PRECISION, dummy, extent, ierror)
    extent = extent * VECTOR_SIZE * NIDX
    CALL MPI_Type_create_resized(datatype_1lev, lb, extent, &
                                 datatype_with_extent, ierror)
    CALL MPI_Type_contiguous(NLEV, datatype_with_extent, datatype_nlev, ierror)
    CALL MPI_Type_commit(datatype_nlev, ierror);
    ! generate reference data
    ALLOCATE(ref_data(NIDX * NBLK * NLEV * VECTOR_SIZE))
    ref_data = -1.0_dp
    v = 1.0
    DO lev = 0, NLEV - 1
      DO i = 1, 9
        DO j = 0, lenghts(i) - 1
          DO vector = 1, VECTOR_SIZE
            ref_data(vector + VECTOR_SIZE * (displs(i) + j) + &
              &      lev * VECTOR_SIZE * NIDX) = v
            v = v + 1.0_dp
          END DO
        END DO
      END DO
    END DO

    ! allocate and initialize receive data
    ALLOCATE(recv_data(NIDX * NBLK * NLEV * VECTOR_SIZE))
    recv_data = -1.0_dp

    ! receive data from rank 1
    CALL MPI_Recv(recv_data, 1, datatype_nlev, 1, TAG, comm, &
                  MPI_STATUS_IGNORE, ierror);

    ! check data
    test_mpich_bug = ANY(ref_data /= recv_data)
    !DO i = 1,  NIDX * NBLK * NLEV * VECTOR_SIZE
    !  IF (ref_data(i) /= recv_data(i)) &
    !    print *, "i", i, "ref_data", ref_data(i), "recv_data", recv_data(i)
    !END DO

    DEALLOCATE(ref_data, recv_data)
    ! free generated datatypes
    CALL MPI_Type_free(datatype_nlev, ierror)
    CALL MPI_Type_free(datatype_with_extent, ierror)
    CALL MPI_Type_free(datatype_1lev, ierror)
    CALL MPI_Type_free(datatype_vector, ierror)

  ELSE IF (comm_rank == 1) THEN

    ! send requested data
    ALLOCATE(send_data(VECTOR_SIZE * NLEV * num_exchange_indices))
    DO i = 1, VECTOR_SIZE * NLEV * num_exchange_indices
      send_data(i) = REAL(i, dp)
    END DO
    CALL MPI_Send(send_data, VECTOR_SIZE * NLEV * num_exchange_indices, &
                  p_real_dp, 0, TAG, comm, ierror)
    DEALLOCATE(send_data)
  END IF

  CALL p_bcast(test_mpich_bug, 0, comm)

END FUNCTION test_mpich_bug

#endif
! HAVE_YAXT
END MODULE mo_communication_yaxt
