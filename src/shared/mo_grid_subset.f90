!>
!! This module contains utilty types
!!
!! @par Revision History
!!  Created by Leonidas Linardakis, MPIM (2012-03-06)
!!
!! @par Copyright
!! 2002-2007 by DWD and MPI-M
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
!! This code has been tested up to a certain mask. Defects and weaknesses,
!! which may be included in the code, do not establish any warranties by the
!! authors.
!! The authors do not make any warranty, express or implied, or assume any
!! liability or responsibility for the use, acquisition or application of this
!! software.
!!
!!
MODULE mo_grid_subset

  USE mo_kind,           ONLY: wp
  USE mo_exception,      ONLY: warning, finish
  USE mo_model_domain,   ONLY: t_patch, t_subset_range, t_subset_range_index, t_subset_indexed
  USE mo_mpi,            ONLY: get_my_mpi_work_id
  USE mo_impl_constants, ONLY: on_cells, on_edges, on_vertices
  USE mo_parallel_config, ONLY: nproma

  IMPLICIT NONE
  INCLUDE 'netcdf.inc'

  CHARACTER(len=*), PARAMETER, PRIVATE :: version = '$Id$'

  PUBLIC :: t_subset_range, t_subset_range_index, t_subset_indexed

  PUBLIC :: fill_subset, get_index_range
  PUBLIC :: fill_subset_from_global_index
  PUBLIC :: read_subset, write_subset
  PUBLIC :: block_no, index_no
  PUBLIC :: get_oriented_edges_from_global_vertices

CONTAINS

  !----------------------------------------------------
  !>
  ! Fills the subset_range with the indexes which statisfy
  !     start_mask <= mask <= end_mask
  !
  ! Assumes that invalid places in mask have value outside the
  ! start_mask - end_mask
  !
  ! Assumes that mask is of the shape (1:,1:)
  SUBROUTINE fill_subset(subset, patch, mask, start_mask, end_mask, subset_name, located)
    TYPE(t_subset_range), INTENT(inout) :: subset
    TYPE(t_patch), TARGET, INTENT(in) :: patch  ! nag does not return the values in subset
    INTEGER, OPTIONAL, INTENT(in) :: located
    CHARACTER(len=32), OPTIONAL :: subset_name
                                                ! unless the patch is declared INTENT(in)!
    INTEGER, INTENT(in) :: mask(:,:), start_mask, end_mask

    INTEGER :: masks_size(2)
    INTEGER :: block, index_in_block, start_index, end_index

    CHARACTER(*), PARAMETER :: method_name = "mo_grid_subset:fill_subset"

    masks_size = SHAPE(mask)

    subset%start_block = -1
    subset%start_index = -1
    subset%end_block   = -2
    subset%end_index   = -2
    subset%block_size  = masks_size(1)
    subset%no_of_holes = 0
    subset%size        = 0
    subset%entity_location  = 0
    subset%patch       => patch

    IF (PRESENT(located)) subset%entity_location = located

    DO block=1,masks_size(2)
      DO index_in_block=1,subset%block_size

        IF (mask(index_in_block, block) >= start_mask .AND. &
            mask(index_in_block, block) <= end_mask) THEN
          ! we found an elemant in range
          IF (subset%start_block < 0) THEN
            ! this is the first element
            subset%start_block = block
            subset%start_index = index_in_block
          ENDIF
          ! this is the up-to-now last element
          subset%end_block = block
          subset%end_index = index_in_block
        ENDIF

      ENDDO
    ENDDO

!     IF (subset%start_block == -1) THEN
!       CALL warning(method_name, "Empty range subset")
!     ENDIF

    IF (subset%start_block > -1) THEN
      ! count the holes
      DO block = subset%start_block, subset%end_block

        start_index = 1
        end_index = subset%block_size
        IF (block == subset%start_block) start_index = subset%start_index
        IF (block == subset%end_block)   end_index = subset%end_index

        DO index_in_block = start_index, end_index

          IF (mask(index_in_block, block) < start_mask .OR. &
              mask(index_in_block, block) > end_mask) THEN
            ! element is a hole in the range
            subset%no_of_holes = subset%no_of_holes + 1
          ENDIF

        ENDDO
      ENDDO
    ENDIF

    ! compute size
    subset%size = 0
    IF (subset%end_block > 0) THEN
      IF ((subset%end_block - subset%start_block) > 1) &
        subset%size = (subset%end_block - subset%start_block -1) * subset%block_size
      subset%size = subset%size + subset%end_index + (subset%block_size - subset%start_index + 1)
    ENDIF
    IF (PRESENT(subset_name)) &
      & subset%name = TRIM(subset_name)

!     IF (subset%no_of_holes > 0) THEN
!       CALL warning(method_name, "We have holes in the range subset")
!     ENDIF

  END SUBROUTINE fill_subset
  !----------------------------------------------------

  !----------------------------------------------------
  !>
  ! Fills the subset_indexed with the indexes defined by the global_index_array
  ! The global_index_array end is define either by its size or the first non-positive integer
  SUBROUTINE fill_subset_from_global_index(subset, patch, global_index_array, subset_name, located)
    TYPE(t_subset_indexed), INTENT(inout) :: subset
    TYPE(t_patch), TARGET, INTENT(in) :: patch  ! nag does not return the values in subset
    INTEGER :: global_index_array(:)   ! intent in
    INTEGER, INTENT(in) :: located     ! =on_cells, on_edges, or on_vertices
    CHARACTER(len=32), OPTIONAL :: subset_name

    INTEGER :: my_proc_id, i, max_allocation_size, owned_indexes, local_idx
    INTEGER, ALLOCATABLE :: tmp_local_index_array(:)
    INTEGER, POINTER :: local_index_array(:), owner_local(:)

    CHARACTER(*), PARAMETER :: method_name = "mo_grid_subset:fill_subset_from_global_index"
    subset%size               = 0
    subset%recommended_stride = 0
    subset%entity_location    = 0
    subset%patch              => patch

    my_proc_id = get_my_mpi_work_id()

    subset%entity_location = located
    SELECT CASE( subset%entity_location )
      CASE( on_cells )
        local_index_array => patch%cells%loc_index
        owner_local       => patch%cells%owner_local
      CASE( on_edges )
        local_index_array => patch%edges%loc_index
        owner_local       => patch%edges%owner_local
      CASE( on_vertices )
        local_index_array => patch%verts%loc_index
        owner_local       => patch%verts%owner_local
      CASE default
        CALL finish(method_name, "Unknown subset%entity_location")
    END SELECT

    ! temporary array for keeping track of what's local
    max_allocation_size = SIZE(global_index_array)
    DO i=1, max_allocation_size
      IF ( global_index_array(i) <= 0 ) EXIT
    ENDDO
    IF ( global_index_array(i) <= 0) &
       max_allocation_size = i - 1

    ALLOCATE(tmp_local_index_array(max_allocation_size))

    owned_indexes = 0
    DO i=1, max_allocation_size
      local_idx = local_index_array(global_index_array(i))
      IF (local_idx > 0) THEN
        IF (owner_local(local_idx) == my_proc_id) THEN
          owned_indexes = owned_indexes + 1
          tmp_local_index_array(owned_indexes) = local_idx
        ENDIF
     ENDIF
    ENDDO

    !now fill the subset if not empty
    IF (owned_indexes > 0) THEN

      subset%size  = owned_indexes
      subset%recommended_stride = 1  ! needs to be calculated
      ALLOCATE(subset%block(owned_indexes), subset%idx(owned_indexes))

      DO i=1, owned_indexes
        subset%block(i) = block_no(tmp_local_index_array(i))
        subset%idx(i)   = index_no(tmp_local_index_array(i))
      ENDDO

    ENDIF

    DEALLOCATE(tmp_local_index_array)

  END SUBROUTINE fill_subset_from_global_index
  !----------------------------------------------------

  !----------------------------------------------------
  !>
  ! Fills the subset_indexed with the indexes defined by the global_index_array
  ! The global_index_array end is define either by its size or the first non-positive integer
  SUBROUTINE get_oriented_edges_from_global_vertices(edge_subset, orientation, patch, global_vertex_array, subset_name)
    TYPE(t_subset_indexed), INTENT(inout) :: edge_subset
    REAL(wp), ALLOCATABLE :: orientation(:)
    TYPE(t_patch), TARGET, INTENT(in) :: patch  ! nag does not return the values in subset
    INTEGER :: global_vertex_array(:)   ! intent in
    CHARACTER(len=32), OPTIONAL :: subset_name

    INTEGER :: my_proc_id, i, max_allocation_size, owned_edges
    INTEGER :: edge_index, edge_block, edge_orientation
    INTEGER :: local_vertex_1d_index(2), vertex_index(2), vertex_block(2)
    INTEGER, ALLOCATABLE :: tmp_edge_block_array(:), tmp_edge_index_array(:)
    REAL(wp), ALLOCATABLE :: tmp_orientation(:)
    INTEGER, POINTER :: local_vertex_array(:), owner_edge_local(:)

    CHARACTER(*), PARAMETER :: method_name = "mo_grid_subset:get_oriented_edges_from_global_vertices"
    edge_subset%size               = 0
    edge_subset%recommended_stride = 0
    edge_subset%entity_location    = on_edges
    edge_subset%patch              => patch
    local_vertex_array             => patch%verts%loc_index
    owner_edge_local               => patch%edges%owner_local

    my_proc_id = get_my_mpi_work_id()

    ! temporary array for keeping track of what's local
    max_allocation_size = SIZE(global_vertex_array)
    DO i=1, max_allocation_size
      IF ( global_vertex_array(i) <= 0 ) EXIT
    ENDDO
    IF ( global_vertex_array(i) <= 0) &
       max_allocation_size = i - 1

    ALLOCATE(tmp_edge_block_array(max_allocation_size), &
      & tmp_edge_index_array(max_allocation_size),      &
      & tmp_orientation(max_allocation_size))
    owned_edges = 0
    IF (max_allocation_size > 1) &
      & local_vertex_1d_index(1) = local_vertex_array(global_vertex_array(1))
    DO i=2, max_allocation_size

      ! get a pair of vertices
      local_vertex_1d_index(2) = local_vertex_array(global_vertex_array(i))

      IF (local_vertex_1d_index(1) > 0 .AND. local_vertex_1d_index(2) > 0 ) THEN
        ! find the edge and orientation of these vertices
        vertex_block(1) = block_no(local_vertex_1d_index(1))
        vertex_index(1) = index_no(local_vertex_1d_index(1))
        vertex_block(2) = block_no(local_vertex_1d_index(2))
        vertex_index(2) = index_no(local_vertex_1d_index(2))

        ! CALL find_oriented_edge_from_vertices(vertex_block, vertex_index, edge_block, edge_index, edge_orientation)

        IF (owner_edge_local(index_1d(idx=edge_index, block=edge_block)) == my_proc_id) THEN
          owned_edges = owned_edges + 1
          tmp_edge_block_array(owned_edges) = edge_block
          tmp_edge_index_array(owned_edges) = edge_index
          tmp_orientation(owned_edges)      = edge_orientation
        ENDIF

     ENDIF

     ! get next pair
     local_vertex_1d_index(1) = local_vertex_1d_index(2)

    ENDDO

    !now fill the edge_subset if not empty
    IF (owned_edges > 0) THEN

      edge_subset%size = owned_edges
      edge_subset%recommended_stride = 1  ! needs to be calculated
      ALLOCATE(edge_subset%block(owned_edges), &
        &      edge_subset%idx(owned_edges),   &
        &      orientation(owned_edges))

      DO i=1, owned_edges
        edge_subset%block(i) = tmp_edge_block_array(i)
        edge_subset%idx(i)   = tmp_edge_index_array(i)
        orientation(i)       = tmp_orientation(i)
      ENDDO

    ENDIF

    DEALLOCATE(tmp_edge_block_array, tmp_edge_index_array, tmp_orientation)

  END SUBROUTINE get_oriented_edges_from_global_vertices
  !----------------------------------------------------

  !----------------------------------------------------
  SUBROUTINE get_index_range(subset_range, current_block, start_index, end_index)
    TYPE(t_subset_range), INTENT(in) :: subset_range
    INTEGER, INTENT(in) :: current_block
    INTEGER, INTENT(out) :: start_index, end_index

    start_index = 1
    end_index = subset_range%block_size

    IF (current_block == subset_range%start_block) &
      start_index = subset_range%start_index
    IF (current_block == subset_range%end_block) &
      end_index = subset_range%end_index

  END SUBROUTINE get_index_range
  !----------------------------------------------------

  !-------------------------------------------------------------------------
  ! The following functions are for conversion of 1D to 2D indices and vice versa
  !
  ! Treatment of 0 (important for empty patches) and negative numbers:
  !
  ! Converting 1D => 2D:
  !
  ! 0 always is mapped to blk_no = 0, idx_no = 0
  ! negative numbers: Convert usings ABS(j) and negate idx_no
  !
  ! Thus: blk_no >= 0 always!
  !       idx_no > 0  for j > 0
  !       idx_no = 0  for j = 0
  !       idx_no < 0  for j < 0
  !
  ! This mimics mostly the behaviour of reshape_idx in mo_model_domimp_patches
  ! with a difference for nproma=1 and j=0 (where reshape_idx returns blk_no=0, idx_no=1)
  !
  ! The consisten treatment of 0 in the above way is very important for empty patches
  ! where start_index=1, end_index=0
  !
  !-------------------------------------------------------------------------
  ELEMENTAL INTEGER FUNCTION block_no(j)
    INTEGER, INTENT(in) :: j

    IF ( j == 0 ) THEN
      block_no = 0
    ELSE
      block_no = (ABS(j)-1)/nproma + 1
    ENDIF

  END FUNCTION block_no
  !-------------------------------------------------------------------------
  ELEMENTAL INTEGER FUNCTION index_no(j)
    INTEGER, INTENT(in) :: j

    IF( j == 0 ) THEN
      index_no = 0
    ELSE
      index_no = SIGN(MOD(ABS(j)-1,nproma)+1, j)
    ENDIF

  END FUNCTION index_no
  !-------------------------------------------------------------------------
  ELEMENTAL INTEGER FUNCTION index_1d(idx,block)
    INTEGER, INTENT(in) :: idx, block

    IF( block <= 0 ) THEN
      index_1d = 0 ! This covers the special case nproma==1,jb=0,jl=1
                 ! All other cases are invalid and get also a 0 returned
    ELSE
      index_1d = SIGN((block-1)*nproma + ABS(idx), idx)
    ENDIF

  END FUNCTION index_1d
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE read_subset(ncid, subset_range, patch)
    INTEGER, INTENT(in) :: ncid
    TYPE(t_subset_range), INTENT(inout) :: subset_range
    TYPE(t_patch), TARGET, INTENT(in) :: patch  ! nag does not return the values in subset_range

    INTEGER :: netcd_status
    CHARACTER(*), PARAMETER :: method_name = "read_subset_range"

    netcd_status = nf_get_att_int(ncid, nf_global,TRIM(subset_range%name)//'.start_block', subset_range%start_block)
    IF (netcd_status /= nf_noerr) &
      & CALL finish(method_name, "Could not read start_block")
    netcd_status = nf_get_att_int(ncid, nf_global,TRIM(subset_range%name)//'.start_index', subset_range%start_index)
    IF (netcd_status /= nf_noerr) &
      & CALL finish(method_name, "Could not read start_index")

    netcd_status = nf_get_att_int(ncid, nf_global,TRIM(subset_range%name)//'.end_block', subset_range%end_block)
    IF (netcd_status /= nf_noerr) &
      & CALL finish(method_name, "Could not read end_block")
    netcd_status = nf_get_att_int(ncid, nf_global,TRIM(subset_range%name)//'.end_index', subset_range%end_index)
    IF (netcd_status /= nf_noerr) &
      & CALL finish(method_name, "Could not read end_index")

    netcd_status = nf_get_att_int(ncid, nf_global,TRIM(subset_range%name)//'.block_size', subset_range%block_size)
    IF (netcd_status /= nf_noerr) &
      & CALL finish(method_name, "Could not read block_size")

    netcd_status = nf_get_att_int(ncid, nf_global,TRIM(subset_range%name)//'.entity_location', subset_range%entity_location)
    IF (netcd_status /= nf_noerr) &
      & CALL finish(method_name, "Could not read entity_type")

    netcd_status = nf_get_att_int(ncid, nf_global,TRIM(subset_range%name)//'.no_of_holes', subset_range%no_of_holes)
    IF (netcd_status /= nf_noerr) &
      & CALL finish(method_name, "Could not read no_of_holes")

    netcd_status = nf_get_att_int(ncid, nf_global,TRIM(subset_range%name)//'.is_in_domain', subset_range%is_in_domain)
    IF (netcd_status /= nf_noerr) &
      & CALL finish(method_name, "Could not read is_in_domain")

    subset_range%patch => patch

  END SUBROUTINE read_subset
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE write_subset(ncid, subset_range)
    INTEGER, INTENT(in) :: ncid
    TYPE(t_subset_range), INTENT(inout) :: subset_range

    INTEGER :: netcd_status
    CHARACTER(*), PARAMETER :: method_name = "write_subset_range"

    CALL nf(nf_put_att_int(ncid, nf_global,TRIM(subset_range%name)//'.start_block', nf_int, 1,     &
      &  subset_range%start_block))
    CALL nf(nf_put_att_int(ncid, nf_global,TRIM(subset_range%name)//'.start_index', nf_int, 1,     &
      &  subset_range%start_index))

    CALL nf(nf_put_att_int(ncid, nf_global,TRIM(subset_range%name)//'.end_block',   nf_int, 1,     &
      & subset_range%end_block))
    CALL nf(nf_put_att_int(ncid, nf_global,TRIM(subset_range%name)//'.end_index',   nf_int, 1,     &
      & subset_range%end_index))

    CALL nf(nf_put_att_int(ncid, nf_global,TRIM(subset_range%name)//'.block_size',  nf_int, 1,     &
      & subset_range%block_size))

    CALL nf(nf_put_att_int(ncid, nf_global,TRIM(subset_range%name)//'.entity_location', nf_int, 1,     &
      & subset_range%entity_location))

    CALL nf(nf_put_att_int(ncid, nf_global,TRIM(subset_range%name)//'.no_of_holes', nf_int, 1,     &
      & subset_range%no_of_holes))

    CALL nf(nf_put_att_int(ncid, nf_global,TRIM(subset_range%name)//'.is_in_domain',nf_int, 1,     &
      & subset_range%is_in_domain))

  END SUBROUTINE write_subset
  !-------------------------------------------------------------------------

  !--------------------------------------------------------------------
  SUBROUTINE nf(return_status)
    INTEGER, INTENT(in) :: return_status

    IF (return_status /= nf_noerr) THEN
      CALL finish('mo_io_grid netCDF error', nf_strerror(return_status))
    ENDIF

  END SUBROUTINE nf
  !-------------------------------------------------------------------------

END MODULE mo_grid_subset

