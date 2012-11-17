!-------------------------------------------------------------------------------------
!>
!! mo_local_grid_hierarchy first implementation
!! Provides basic tools for creating the grid hierarchy
!
!! $Id: n/a$
!!
!! @par Revision History
!! Previous implementations in
!! mo_gridrefinement.f90, mo_geometry.f90 mo_topolgy.f90 by
!!
!! Current implementation by
!!   Leonidas Linardakis, MPI-M, 2010-01-12
!!
!! @par Copyright
!! 2002-2009 by DWD and MPI-M
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
MODULE mo_local_grid_hierarchy
#include "grid_definitions.inc"
  !-------------------------------------------------------------------------
!  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: message_text, message, finish
  USE mo_local_grid
  USE mo_io_local_grid,      ONLY: read_netcdf_grid, write_netcdf_grid
!  USE mo_base_geometry,      ONLY: gvec2cvec
  USE mo_impl_constants,     ONLY: min_rledge, max_rledge, min_rlvert, &
    & max_rlvert, min_rledge_int, min_rlvert_int

  IMPLICIT NONE


  ! !VERSION CONTROL:
  CHARACTER(LEN=*), PARAMETER :: version = '$Id$'

  PRIVATE

  PUBLIC :: create_grid_hierarchy
  PUBLIC :: set_bdy_indexing_depth

  INTEGER :: max_boundary_depth
  INTEGER :: inner_boundary_depth  ! = -MAX_BOUNDARY_LEVEL  ?
  INTEGER :: outer_boundary_depth  ! = 2*MAX_BOUNDARY_LEVEL ?
  INTEGER :: inner_edge_boundary_depth, outer_edge_boundary_depth
  INTEGER :: bdy_indexing_depth = 0 ! if > outer_boundary_depth, use this
  INTEGER :: min_edge_allocate, min_allocate ! min index maybe smaller than inner_boundary_depth
  !-------------------------------------------------------------------------

!   INTEGER,  ALLOCATABLE:: inner_boundary_start(:,:),outer_boundary_start(:)
!   INTEGER,  ALLOCATABLE:: inner_boundary_end(:,:),  outer_boundary_end(:)

#define parent_is_triangle(parent_type) (parent_type >= parenttype_triangle)


CONTAINS


  !-------------------------------------------------------------------------
  SUBROUTINE set_bdy_indexing_depth(depth)
    INTEGER, INTENT(in) :: depth

    bdy_indexing_depth = depth

  END SUBROUTINE set_bdy_indexing_depth
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE define_boundary_depths()

    max_boundary_depth        = max_rlvert
    inner_boundary_depth      = min_rlvert_int
    outer_boundary_depth      = max_boundary_depth
    inner_edge_boundary_depth = min_rledge_int
    outer_edge_boundary_depth = max_rledge
    min_edge_allocate         = min_rledge
    min_allocate              = min_rlvert
    
!     write(0,*) 'boundary_depths:', inner_boundary_depth,&
!       outer_boundary_depth
!     write(0,*) 'edge_boundary_depths:', inner_edge_boundary_depth,&
!       outer_edge_boundary_depth
    ! it should be:
    ! INNER_EDGE_BOUNDARY_DEPTH = 2*INNER_BOUNDARY_DEPTH
    ! OUTER_EDGE_BOUNDARY_DEPTH = 2*OUTER_BOUNDARY_DEPTH

  END SUBROUTINE define_boundary_depths
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE create_grid_hierarchy(grid_id)
    INTEGER, INTENT(in) :: grid_id

    CALL define_boundary_depths()
    CALL create_rec_grid_hierarchy(grid_id)

  END SUBROUTINE create_grid_hierarchy
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !   RECURSIVE SUBROUTINE create_grid_hierarchy(grid_id)
  !>
  !! Main method for creating the grid hierarchy.
  !! Reads a (nested) set of grids, assigns the child pointers,
  !! computes the outer and inner boundary levels and reindexes
  !! the grid entities according to the levels.
  !! Private
  RECURSIVE SUBROUTINE create_rec_grid_hierarchy(grid_id)
    INTEGER, INTENT(in) :: grid_id

    INTEGER :: current_grid_id, child_grid_id, no_of_children
    INTEGER :: i

    current_grid_id = grid_id
    no_of_children = get_grid_no_of_children(current_grid_id)

    WRITE(message_text,'(i3)') current_grid_id
    CALL message ('create_grid_hierarchy for grid ', TRIM(message_text))

    ! first process the children
    DO i=1,no_of_children
      CALL create_rec_grid_hierarchy(get_grid_child_id(current_grid_id,i))
    ENDDO

    ! read grid
    IF (.not. grid_is_filled(current_grid_id)) &
      CALL read_netcdf_grid(current_grid_id)

    CALL message ('  get the children pointers ', '')
    ! the children info is cleaned by the construction
    ! get the children pointers from the children
    DO i=1,no_of_children
      CALL get_child_pointers(current_grid_id, get_grid_child_id(current_grid_id,i))
    ENDDO

    ! get the boundary levels, ie compute the refine control parameter
    CALL message ('  compute_boundary_levels ', '')
    CALL compute_boundary_levels(current_grid_id)

    !re-index
    CALL message ('  re-index ', '')
    CALL re_index_grid(current_grid_id)

    ! write and delete the children
    CALL message ('  write and delete the children ', '')
    DO i=1,no_of_children
      child_grid_id = get_grid_child_id(current_grid_id,i)
      ! CALL check_parent_child_grid(child_grid_id)
      CALL write_netcdf_grid(child_grid_id)
      CALL delete_grid(child_grid_id)
    ENDDO

    WRITE(message_text,'(i3,a)') current_grid_id, ' is done.'
    CALL message ('create_grid_hierarchy for grid ', TRIM(message_text))
    RETURN

  END SUBROUTINE create_rec_grid_hierarchy
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  ! NOTE: verts%parent_index with edge parents are not updated
  ! for the moment we do not use verts%parent_index
  SUBROUTINE re_index_grid(grid_id)
    INTEGER, INTENT(in) :: grid_id

    TYPE(t_grid),          POINTER :: current_grid,child_grid
    TYPE(t_grid_vertices), POINTER :: verts
    TYPE(t_grid_edges),    POINTER :: edges
    TYPE(t_grid_cells),    POINTER :: cells
    INTEGER :: no_of_verts, no_of_edges, no_of_cells

    INTEGER, POINTER :: vertex_new_index(:) ! (0:no_of_verts)
    INTEGER, POINTER :: edge_new_index(:)   ! (0:no_of_edges)
    INTEGER, POINTER :: cell_new_index(:)   ! (0:no_of_cells)
    INTEGER, POINTER :: verts_inner_boundary_start(:,:), verts_inner_boundary_end(:,:)
    INTEGER, POINTER :: verts_outer_boundary_start(:),   verts_outer_boundary_end(:)
    INTEGER, POINTER :: edges_inner_boundary_start(:,:), edges_inner_boundary_end(:,:)
    INTEGER, POINTER :: edges_outer_boundary_start(:),   edges_outer_boundary_end(:)
    INTEGER, POINTER :: cells_inner_boundary_start(:,:), cells_inner_boundary_end(:,:)
    INTEGER, POINTER :: cells_outer_boundary_start(:),   cells_outer_boundary_end(:)

    ! INTEGER :: child_id, child_id2, no_of_children
    INTEGER :: child_grid_id, no_of_children, new_grid_id
    INTEGER :: i, k,i_status

    current_grid  => get_grid(grid_id)
    no_of_children =  current_grid%no_of_children
    verts         => current_grid%verts
    edges         => current_grid%edges
    cells         => current_grid%cells
    no_of_verts    = verts%no_of_existvertices
    no_of_edges    = edges%no_of_existedges
    no_of_cells    = cells%no_of_existcells

    !--------------------------------
!     i = MIN(inner_boundary_depth,inner_edge_boundary_depth)
!     ALLOCATE(inner_boundary_size(i:-1,1:max_no_of_grid_objects), &
!       & inner_boundary_start(i:-1,1:max_no_of_grid_objects), &
!       & inner_boundary_end(i:-1,1:max_no_of_grid_objects),stat=i_status)
!     IF (i_status > 0) &
!       & CALL finish ('re_index_grid', 'ALLOCATE(boundarySize)')
!     i = MAX(outer_boundary_depth,outer_edge_boundary_depth)
!     ALLOCATE(outer_boundary_size(0:i), &
!       & outer_boundary_start(0:i), &
!       & outer_boundary_end(0:i),stat=i_status)
!     IF (i_status > 0) &
!       & CALL finish ('re_index_grid', 'ALLOCATE(boundaryStart)')
    !--------------------------------

    !--------------------------------
    ! get new indeces
    CALL index_level_entities(grid_id, no_of_verts, verts%refin_ctrl, verts%child_id, &
      & inner_boundary_depth,outer_boundary_depth, vertex_new_index, &
      & verts_inner_boundary_start, verts_inner_boundary_end,        &
      & verts_outer_boundary_start, verts_outer_boundary_end)

    CALL index_level_entities(grid_id, no_of_cells, cells%refin_ctrl, cells%child_id, &
      & inner_boundary_depth,outer_boundary_depth, cell_new_index,   &
      & cells_inner_boundary_start, cells_inner_boundary_end,        &
      & cells_outer_boundary_start, cells_outer_boundary_end)

    CALL index_level_entities(grid_id, no_of_edges, edges%refin_ctrl, edges%child_id, &
      & inner_edge_boundary_depth,outer_edge_boundary_depth, edge_new_index, &
      & edges_inner_boundary_start, edges_inner_boundary_end,        &
      & edges_outer_boundary_start, edges_outer_boundary_end)

    new_grid_id = -1
    CALL copy_grid(grid_id, new_grid_id, &
      & vertex_new_index, edge_new_index, cell_new_index)
    CALL replace_grid(grid_id, new_grid_id)
    ! update pointers
    current_grid  => get_grid(grid_id)
    verts         => current_grid%verts
    edges         => current_grid%edges
    cells         => current_grid%cells

    ! reallocate start_idx, end_idx
    IF (no_of_children > 1) THEN
      DEALLOCATE (verts%start_idx,verts%end_idx)
      DEALLOCATE (edges%start_idx,edges%end_idx)
      DEALLOCATE (cells%start_idx,cells%end_idx)
      ALLOCATE(verts%start_idx(min_allocate:outer_boundary_depth, no_of_children),&
        & verts%end_idx(min_allocate:outer_boundary_depth, no_of_children),&
        & cells%start_idx(min_allocate:outer_boundary_depth, no_of_children),&
        & cells%end_idx(min_allocate:outer_boundary_depth, no_of_children),&
        & edges%start_idx(min_edge_allocate:outer_edge_boundary_depth, no_of_children),&
        & edges%end_idx(min_edge_allocate:outer_edge_boundary_depth, no_of_children), &
        & stat=i_status)
      IF (i_status > 0) &
        & CALL finish ('re_index_grid', 'REALLOCATE(start_idx, end_idx)')
     ENDIF

    ! reindex children's parent indexes
    DO i=1,no_of_children
      child_grid_id = current_grid%child_grid_ids(i)
      child_grid => get_grid(child_grid_id)

      CALL re_index_list(child_grid%verts%parent_index, vertex_new_index, &
        & child_grid%verts%no_of_existvertices)

      CALL re_index_list(child_grid%edges%parent_index, edge_new_index, &
        & child_grid%edges%no_of_existedges)

      CALL re_index_list(child_grid%cells%parent_index, cell_new_index, &
        & child_grid%cells%no_of_existcells)

      ! reindex verts%parent_index when parent is an edge
      ! for the moment we do not use verts%parent_index
      child_grid%verts%parent_index(:) = -child_grid%verts%parent_index(:)
      CALL re_index_list(child_grid%verts%parent_index, edge_new_index, &
        & child_grid%verts%no_of_existvertices)
      child_grid%verts%parent_index(:) = -child_grid%verts%parent_index(:)
      child_grid%verts%parent_index(:) = -child_grid%verts%parent_index(:)

      DO k=inner_boundary_depth,-1
        verts%start_idx(k,i) = verts_inner_boundary_start(k,child_grid_id)
        verts%end_idx(k,i)   = verts_inner_boundary_end(k,child_grid_id) - 1

        cells%start_idx(k,i) = cells_inner_boundary_start(k,child_grid_id)
        cells%end_idx(k,i)   = cells_inner_boundary_end(k,child_grid_id) - 1
      ENDDO
      ! from min_rlvert to min_rlvert_int nullify indexes
      DO k=min_allocate, inner_boundary_depth-1
        verts%start_idx(k,i) = no_of_verts+1
        verts%end_idx(k,i)   = no_of_verts

        cells%start_idx(k,i) = no_of_cells+1
        cells%end_idx(k,i)   = no_of_cells
      ENDDO

      child_grid%verts%parent_index(:) = -child_grid%verts%parent_index(:)
      DO k=inner_edge_boundary_depth,-1
        edges%start_idx(k,i) = edges_inner_boundary_start(k,child_grid_id)
        edges%end_idx(k,i)   = edges_inner_boundary_end(k,child_grid_id) - 1
      ENDDO
      ! from min_rledge to min_rledge_int nullify indexes
      DO k=min_edge_allocate,inner_edge_boundary_depth-1
        edges%start_idx(k,i) = no_of_edges+1
        edges%end_idx(k,i)   = no_of_edges
      ENDDO

    ENDDO

    DO k=0,outer_boundary_depth
      verts%start_idx(k,:) = verts_outer_boundary_start(k)
      verts%end_idx(k,:)   = verts_outer_boundary_end(k) - 1
    ENDDO
    DO k=0,outer_boundary_depth
      cells%start_idx(k,:) = cells_outer_boundary_start(k)
      cells%end_idx(k,:)   = cells_outer_boundary_end(k) - 1
    ENDDO
    DO k=0,outer_edge_boundary_depth
      edges%start_idx(k,:) = edges_outer_boundary_start(k)
      edges%end_idx(k,:)   = edges_outer_boundary_end(k) - 1
    ENDDO

    DEALLOCATE(vertex_new_index)
    DEALLOCATE(edge_new_index)
    DEALLOCATE(cell_new_index)

    DEALLOCATE(verts_inner_boundary_start, verts_outer_boundary_start)
    DEALLOCATE(verts_inner_boundary_end,   verts_outer_boundary_end)
    DEALLOCATE(edges_inner_boundary_start, edges_outer_boundary_start)
    DEALLOCATE(edges_inner_boundary_end,   edges_outer_boundary_end)
    DEALLOCATE(cells_inner_boundary_start, cells_outer_boundary_start)
    DEALLOCATE(cells_inner_boundary_end,   cells_outer_boundary_end)

  END SUBROUTINE re_index_grid
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  PURE SUBROUTINE re_index_list(the_list, new_index, list_size)
    INTEGER, INTENT(inout) :: the_list(:)
    INTEGER, INTENT(in)    :: new_index(0:), list_size

    INTEGER :: i

    DO i=1,list_size
      IF (the_list(i) > 0) &
        ! WRITE(*,*) 'reIndexList:',theList(i),newIndex(theList(i))
        the_list(i) = new_index(the_list(i))
    ENDDO

  END SUBROUTINE re_index_list
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE index_level_entities(grid_id, no_of_entities, refin_ctrl_list, child_id_list, &
      & inner_boundary_depth,outer_boundary_depth, new_index_list,   &
      & inner_boundary_start, inner_boundary_end,        &
      & outer_boundary_start, outer_boundary_end)

    INTEGER, INTENT(in) :: grid_id, no_of_entities,inner_boundary_depth,outer_boundary_depth
    INTEGER, INTENT(in) :: refin_ctrl_list(no_of_entities), child_id_list(no_of_entities)
    INTEGER, POINTER :: new_index_list(:) ! (0:no_of_entities)
    INTEGER, POINTER :: inner_boundary_start(:,:), inner_boundary_end(:,:)
    INTEGER, POINTER :: outer_boundary_start(:), outer_boundary_end(:)

    TYPE(t_grid), POINTER ::  grid_obj

    INTEGER,  ALLOCATABLE:: inner_boundary_size(:,:), outer_boundary_size(:)
    INTEGER :: current_level, start_index, no_of_children, child_grid_id
    INTEGER :: i, i_status

   !--------------------------------
    i = MIN(inner_boundary_depth,inner_edge_boundary_depth)
    ALLOCATE(inner_boundary_size(i:-1,1:max_no_of_grid_objects), &
      & inner_boundary_start(i:-1,1:max_no_of_grid_objects), &
      & inner_boundary_end(i:-1,1:max_no_of_grid_objects),stat=i_status)
    IF (i_status > 0) &
      & CALL finish ('re_index_grid', 'ALLOCATE(boundarySize)')
    i = MAX(outer_boundary_depth,outer_edge_boundary_depth)
    ALLOCATE(outer_boundary_size(0:i), &
      & outer_boundary_start(0:i), &
      & outer_boundary_end(0:i),stat=i_status)
    IF (i_status > 0) &
      & CALL finish ('re_index_grid', 'ALLOCATE(boundaryStart)')
    !--------------------------------

    !WRITE(*,*) " indexLevelEntities..."
    ! CALL flush(6)
    grid_obj => get_grid(grid_id)
    no_of_children =  grid_obj%no_of_children


    ALLOCATE(new_index_list(0:no_of_entities),stat=i_status)
    IF (i_status > 0) &
      & CALL finish ('indexLevelEntities', 'ALLOCATE(newIndexList)')
    new_index_list(0) = 0

    ! count the number of entities for each level/subdomain
    inner_boundary_size(:,:) = 0
    outer_boundary_size(:)   = 0

    DO i = 1, no_of_entities
      current_level = refin_ctrl_list(i)
      child_grid_id = child_id_list(i)
      !  WRITE(*,*) i,'currentLevel=', currentLevel
      IF (current_level < 0) THEN
        IF (current_level < inner_boundary_depth) &
          & CALL finish ('indexLevelEntities', 'refin_ctrl < innerBoundaryDepth')
        IF (child_grid_id < 1) THEN
          ! WRITE(*,*) i, currentLevel, " childID = ",childID
          ! CALL FLUSH(6)
          CALL finish ('indexLevelEntities', 'childID problem')
        ENDIF
        inner_boundary_size(current_level,child_grid_id) = &
          & inner_boundary_size(current_level,child_grid_id) + 1
      ELSE
        IF (current_level > outer_boundary_depth) &
          & current_level = 0
!         IF (current_level > outer_boundary_depth) &
!           & CALL finish ('indexLevelEntities', 'refin_ctrl > outerBoundaryDepth')
        outer_boundary_size(current_level) = outer_boundary_size(current_level) + 1
      ENDIF
    ENDDO

    ! compute startLevel
    start_index = 1
    DO current_level=1,outer_boundary_depth
      outer_boundary_start(current_level) = start_index
      start_index = start_index + outer_boundary_size(current_level)
      ! WRITE(*,*) currentLevel, ' outerBoundaryStart=', outerBoundaryStart(currentLevel), &
      !      & outerBoundarySize(currentLevel)
    ENDDO
    current_level=0
    outer_boundary_start(current_level) = start_index
    start_index = start_index + outer_boundary_size(current_level)
    ! WRITE(*,*) currentLevel, ' outerBoundaryStart=', outerBoundaryStart(currentLevel), &
    !      & outerBoundarySize(currentLevel)
    DO i=1,no_of_children
      child_grid_id = grid_obj%child_grid_ids(i)
      DO current_level=-1,inner_boundary_depth,-1
        inner_boundary_start(current_level,child_grid_id) = start_index
        start_index = start_index + inner_boundary_size(current_level, child_grid_id)
        ! WRITE(*,*) childID,currentLevel,' innerBoundaryStart=', &
        !     innerBoundaryStart(currentLevel,childID),innerBoundarySize(currentLevel,childID)
      ENDDO
    ENDDO

    ! compute new indexes
    inner_boundary_end(:,:) = inner_boundary_start(:,:)
    outer_boundary_end(:)   = outer_boundary_start(:)
    DO i = 1, no_of_entities
      current_level = refin_ctrl_list(i)
      child_grid_id = child_id_list(i)
      !  WRITE(*,*) i,'currentLevel=', currentLevel
      IF (current_level < 0) THEN
        new_index_list(i) = inner_boundary_end(current_level,child_grid_id)
        inner_boundary_end(current_level,child_grid_id) = &
          & inner_boundary_end(current_level,child_grid_id) + 1
      ELSE
        IF (current_level > outer_boundary_depth) &
          & current_level = 0
        new_index_list(i) = outer_boundary_end(current_level)
        outer_boundary_end(current_level) = outer_boundary_end(current_level) + 1
      ENDIF
      !    WRITE(*,*) 'n  new index:',i,newIndexList(i)
    ENDDO

    DEALLOCATE(inner_boundary_size,outer_boundary_size)

  END SUBROUTINE index_level_entities
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE compute_boundary_levels(grid_id)
    INTEGER, INTENT(in) :: grid_id

    TYPE(t_grid),       POINTER :: current_grid
    TYPE(t_grid_vertices), POINTER :: verts
    TYPE(t_grid_edges), POINTER :: edges
    TYPE(t_grid_cells), POINTER :: cells
    INTEGER :: no_of_verts, no_of_edges, no_of_cells

    TYPE(t_integer_list) :: boundary_edge_index
    INTEGER :: no_of_boundary_edges

    INTEGER :: edge_index, cell1, cell2, max_cell_vertices
    INTEGER :: child_id1, child_id2, no_of_children
    INTEGER :: cell1_level, cell2_level, level_diff
    INTEGER :: i,j,k, i_status

    current_grid    => get_grid(grid_id)
    no_of_children  =  current_grid%no_of_children
    verts           => current_grid%verts
    edges           => current_grid%edges
    cells           => current_grid%cells
    no_of_verts    = verts%no_of_existvertices
    no_of_edges    = edges%no_of_existedges
    no_of_cells    = cells%no_of_existcells
    max_cell_vertices =  cells%max_no_of_vertices

    ! clean the boundary levels
    verts%refin_ctrl(:) = 0
    edges%refin_ctrl(:) = 0
    cells%refin_ctrl(:) = 0

    ALLOCATE(boundary_edge_index%value(no_of_edges),stat=i_status)
    IF (i_status > 0) &
      & CALL finish ('computeBoundaryLevels', 'ALLOCATE(boundaryEdgeIndex(0:noOfEdges))')
    boundary_edge_index%value(:) = 0
    no_of_boundary_edges    = 0

    ! mark boundary edges
    DO edge_index = 1, no_of_edges
      ! first check if it is outer boundary
      cell1 = edges%get_cell_index(edge_index,1)
      cell2 = edges%get_cell_index(edge_index,2)
      IF (cell1 == 0 .or. cell2 == 0 ) THEN
        ! add edge to boundary list
        no_of_boundary_edges = no_of_boundary_edges + 1
        boundary_edge_index%value(no_of_boundary_edges) = edge_index
        ! print *, 'mark boundary edges: edges%subgrid_id=', edges%subgrid_id(edge_index)
        IF (edges%child_id(edge_index) /= 0) THEN
           CALL print_grid_edge(grid_id,edge_index)
           CALL print_grid_edge(edges%child_id(edge_index),edges%child_index(edge_index,1))
           CALL finish ('computeBoundaryLevels', 'outer boundary overlapping with subdomain')
        ENDIF
        CYCLE ! no need to check further
      ENDIF

      ! check if edge is an inner boundary
      ! ie one cell is in a subdomain, the other is not
      child_id1 = cells%child_id(cell1)
      child_id2 = cells%child_id(cell2)
      IF (child_id1 /= child_id2) THEN
        ! add edge to boundary list
        no_of_boundary_edges = no_of_boundary_edges + 1
        boundary_edge_index%value(no_of_boundary_edges) = edge_index
        IF (child_id1 /= 0 .and. child_id2 /= 0) &
          & CALL finish ('computeBoundaryLevels', 'overlaping subdomains')
      ENDIF ! (childID1 /= childID2)
    ENDDO !edgeIndex = 1, noOfEdges

    ! WRITE(*,*) "computeBoundaryLevels: Boundary edges = ", noOfBoundaryEdges
    boundary_edge_index%list_size = no_of_boundary_edges

    ! boundary edges are marked
    ! calculate the boundary levels (ie refine control) for the outer boundary
    ! and for each of the boundaries with subdomains (childs)

    ! the outer boundary refine control will get values from 1 to  max_rlvert (for vertices)
    child_id1 = 0
    CALL compute_boundary_entities(grid_id, boundary_edge_index, child_id1, &
      & 1, max(outer_boundary_depth,bdy_indexing_depth))
    ! process the internal (nested) boundaries (along the children)
    DO i=1,no_of_children
      CALL compute_boundary_entities(grid_id, boundary_edge_index,current_grid%child_grid_ids(i),&
        & -1, inner_boundary_depth)
    ENDDO

    DEALLOCATE(boundary_edge_index%value)

    ! mark all remaining cells/verts with childs as INNER_BOUNDARY_DEPTH
    DO i=1,no_of_cells
      IF (cells%child_id(i) /= 0 .and. cells%refin_ctrl(i) == 0) THEN
        cells%refin_ctrl(i) = inner_boundary_depth
        DO j=1,max_cell_vertices
          k = cells%get_vertex_index(i,j)
          IF (k > 0) THEN
            IF (verts%child_id(k) /= cells%child_id(i)) THEN
              CALL finish('computeBoundaryLevels','childID of cell and vertex do not match')
            ENDIF
            IF (verts%refin_ctrl(k) == 0) &
              & verts%refin_ctrl(k) = inner_boundary_depth
          ENDIF
        ENDDO
      ENDIF
    ENDDO

    ! finally compute boundary level (refin_ctrl) for the edges
    DO edge_index = 1, no_of_edges
      cell1 = edges%get_cell_index(edge_index,1)
      cell2 = edges%get_cell_index(edge_index,2)

      IF (cell1 == 0 .or. cell2 == 0 ) THEN
        ! we have an outer boundary edge
        edges%refin_ctrl(edge_index) = 1
        !        WRITE(*,*) edgeIndex," edges%refin_ctrl=",edges%refin_ctrl(edgeIndex)
        CYCLE ! no need to check further
      ENDIF

      cell1_level = cells%refin_ctrl(cell1)
      cell2_level = cells%refin_ctrl(cell2)
      level_diff   = cell2_level - cell1_level
!      edges%refin_ctrl(edge_index) = 2 * cell1_level + level_diff
      edges%refin_ctrl(edge_index) =  cell1_level +  cell2_level
      IF ((edges%refin_ctrl(edge_index)) == 1) THEN
         WRITE(*,*) edge_index," edges%refin_ctrl=",edges%refin_ctrl(edge_index) , &
         &   cell1,cell1_level,cell2,cell2_level,level_diff
         CALL print_grid_cell(grid_id, cell1)
         CALL print_grid_cell(grid_id, cell2)
         CALL finish('compute_boundary_levels','refin_ctrl(edgeIndex)) == 1')
      ENDIF
      ! if the edge is between a boundary level and an inner cell will get 0
      IF (ABS(level_diff) > 1) THEN
        IF (cell1_level /= 0 .and. cell2_level /= 0) THEN
          WRITE(message_text,'(i9,i9,a,i3,i3)') cell1,cell2," cell level difference > 1:",&
            & cell1_level, cell2_level
          CALL message('compute_boundary_levels. Warning', message_text)
          ! CALL finish ('compute_boundary_levels', message_text)
        ENDIF
        edges%refin_ctrl(edge_index) = 0
      ENDIF

      IF (edges%refin_ctrl(edge_index) < inner_edge_boundary_depth .or. &
        & edges%refin_ctrl(edge_index) > outer_edge_boundary_depth) THEN
        !        WRITE(*,*) "Out of limits edges%refin_ctrl=",edgeIndex, &
        !             &     edges%refin_ctrl(edgeIndex) ,cell1_level,cell2_level
        !        call flush(6)
        ! this is caused by inconsistence of the definition of OUTER_EDGE_BOUNDARY_DEPTH
        ! for the moment just zero it
        edges%refin_ctrl(edge_index) = 0
        !       CALL finish ('computeBoundaryLevels', 'edges%refin_ctrl out of limits')
      ENDIF


    ENDDO ! edgeIndex = 1, noOfEdges

    RETURN

  END SUBROUTINE compute_boundary_levels
  !-------------------------------------------------------------------------


  ! !-------------------------------------------------------------------------
  ! SUBROUTINE draw_grid_boundary_levels(current_grid)
  !   TYPE(grid),  POINTER  :: current_grid
  !
  !   TYPE(grid_cells), POINTER  :: cells
  !   INTEGER :: no_of_cells
  !
  !   TYPE(grid)  :: boundary_grid
  !   INTEGER, ALLOCATABLE :: boundary_grid_cells(:)
  !   INTEGER :: no_ofboundary_grid_cells
  !   CHARACTER(len=filename_max) :: boundary_file_name
  !
  !  INTEGER :: i, i_status
  !
  !   cells        => current_grid%cells
  !   no_of_cells    =  current_grid%ncells
  !
  !   ! draw external boundary elements
  !   ALLOCATE(boundary_grid_cells(0:no_of_cells),stat=i_status)
  !   IF (i_status > 0) &
  !      CALL finish ('drawGridBoundaryLevels', 'ALLOCATE(boundaryGridCells(0:noOfCells))')
  !   boundary_grid_cells(:)  = 0
  !   no_ofboundary_grid_cells = 0
  !   DO i=1,no_of_cells
  !     IF (cells%refin_ctrl(i) > 0) THEN
  !        no_ofboundary_grid_cells = no_ofboundary_grid_cells + 1
  !        boundary_grid_cells(no_ofboundary_grid_cells) = i
  !     ENDIF
  !   ENDDO
  !   IF (no_ofboundary_grid_cells > 0) THEN
  !     boundary_grid_cells(0) = no_ofboundary_grid_cells
  !     CALL create_gridfrom_cell_list(current_grid, boundary_grid, boundary_grid_cells)
  !     WRITE(boundary_file_name,'(i2.2,a)' ) current_grid%grid_id, ".nc"
  !     WRITE(*,*) "writitng boundary to:",trim(boundary_file_name)
  !     CALL write_grid(boundary_grid, boundary_file_name)
  !     CALL destruct_grid(patch(root_id)%grid_id)
  !   ENDIF
  !   DEALLOCATE(boundary_grid_cells)
  !
  ! END SUBROUTINE draw_grid_boundary_levels
  ! !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE compute_boundary_entities(grid_id, boundary_edge_index, child_id, &
      & start_level, end_level)
    ! computes the boundary levels (refin_ctrl) for cells and vertices
    INTEGER, INTENT(in) :: grid_id, child_id, start_level, end_level
    TYPE(t_integer_list), INTENT(in) :: boundary_edge_index

    TYPE(t_grid),       POINTER :: current_grid
    TYPE(t_grid_edges), POINTER :: edges
    TYPE(t_grid_cells), POINTER :: cells
    TYPE(t_grid_vertices), POINTER :: verts
    INTEGER :: no_of_edges, no_of_boundary_edges, no_of_cells, no_of_verts
    INTEGER :: max_vertex_cells, max_cell_vertices

    INTEGER, ALLOCATABLE :: init_vertex_stack(:), next_vertex_stack(:)
    INTEGER :: stack_index, init_vertex_stack_size, next_vertex_stack_size

    INTEGER :: vertex_index, cell_vertex_index, cell_index, edge_index
    INTEGER :: i,j,level,step_level,i_status
    
    CHARACTER(*), PARAMETER :: method_name = "compute_boundary_entities"

    current_grid => get_grid(grid_id)
    edges       => current_grid%edges
    cells       => current_grid%cells
    verts       => current_grid%verts

    no_of_verts       = verts%no_of_existvertices
    no_of_edges       = edges%no_of_existedges
    no_of_cells       = cells%no_of_existcells
    max_vertex_cells  = verts%max_connectivity
    max_cell_vertices = cells%max_no_of_vertices

    no_of_boundary_edges = boundary_edge_index%list_size

    WRITE(*,*) method_name, ": nodeID, childID=", grid_id,child_id
    ! WRITE(*,*) "computeBoundaryEntities, nodeID, childID, noOfBoundaryEdges=", &
    !      & nodeID, childID, noOfBoundaryEdges
    ! fill initial stack of boundary cells
    ALLOCATE(init_vertex_stack(no_of_verts),next_vertex_stack(no_of_verts),stat=i_status)
    IF (i_status > 0) &
      & CALL finish (method_name, 'ALLOCATE(init_vertex_stack...)')
    init_vertex_stack_size = 0
    DO i=1, no_of_boundary_edges
      edge_index = boundary_edge_index%value(i)
      IF (edges%child_id(edge_index) == child_id) THEN
        ! we have a hit
        ! add the related vertices to the stack
        DO j=1,2
          vertex_index = edges%get_vertex_index(edge_index,j)
          IF (verts%refin_ctrl(vertex_index) == 0) THEN
            ! WRITE(*,*) "adding init edge,vertex ",edgeIndex,vertexIndex,startLevel
            ! add to stack and mark
            init_vertex_stack_size = init_vertex_stack_size + 1
            init_vertex_stack(init_vertex_stack_size) = vertex_index
            verts%refin_ctrl(vertex_index)        = start_level
            IF (verts%child_id(vertex_index) /= child_id) THEN
              WRITE(message_text,'(a,i9,i3,i3)') " Boundary vertex, child_id mismatch",&
                & vertex_index, verts%child_id(vertex_index),child_id
              CALL finish (method_name, message_text)
            ENDIF
            ! verts%child_id(vertex_index)          = child_id
          ENDIF
        ENDDO
      ENDIF
    ENDDO

    ! we have the initial stack of points
    ! go through all boundary levels and mark the grid boundary entities
    step_level = 1
    IF (start_level > end_level) step_level = -1
    DO level = start_level, end_level, step_level
      ! go through the initial stack and
      !   1. mark the cells
      !   2. add the next level of points to a new stack and mark them
      ! move the new stack to the previous
      next_vertex_stack_size = 0
      ! WRITE(*,*) level, " level,  StackSize:",initVertexStackSize
      DO stack_index=1, init_vertex_stack_size
        vertex_index = init_vertex_stack(stack_index)
        DO i = 1,max_vertex_cells
          cell_index = verts%get_cell_index(vertex_index,i)
          !         WRITE(*,*) " got cell ", cellIndex
          IF (cell_index /= 0) THEN
!               IF (cell_index == 1069375) &
!                 & print*,'checking cell=',cell_index
            ! if the cell is not marked and matches the current child ID,
            ! mark it and add the vertices to the next stack
            !  WRITE(*,*) " checking cell ", cellIndex, cells%refin_ctrl(cellIndex), &
            !       &     cells%child_id(cellIndex)
            IF (cells%refin_ctrl(cell_index) == 0 &
              & .and. cells%child_id(cell_index) == child_id) THEN
              ! mark the cell
              cells%refin_ctrl(cell_index) = level
              ! if this is not the last level then
              ! add the cell vertices to the new stack and mark them (if not alread marked)
              IF (level /= end_level) THEN
                !-------------------------------------
                DO j=1,max_cell_vertices
                  cell_vertex_index = cells%get_vertex_index(cell_index, j)
                  IF (cell_vertex_index /= 0) THEN
!                     IF (cell_index == 1069375) &
!                      & print*,'cell=',cell_index, ' checking vertex ',&
!                      cell_vertex_index,verts%refin_ctrl(cell_vertex_index)
                    IF (verts%refin_ctrl(cell_vertex_index) == 0) THEN
                      ! add newVertex to stack and mark
!                       IF (cell_index == 1069375) &
!                         & print*,'cell=',cell_index,' level=',level,&
!                         & ' adding vertex:',cell_vertex_index
                      next_vertex_stack_size = next_vertex_stack_size + 1
                      next_vertex_stack(next_vertex_stack_size) = cell_vertex_index
                      verts%refin_ctrl(cell_vertex_index)        = level + step_level
                      IF (verts%child_id(vertex_index) /= child_id) THEN
                        WRITE(message_text,'(i2,a,i9,i3,i3)') level, &
                          & " level. Vertex, child_id mismatch", &
                          & vertex_index, verts%child_id(vertex_index),child_id
                        CALL finish (method_name, message_text)
                      ENDIF
                      ! verts%child_id(cell_vertex_index)          = child_id
                      ! WRITE(*,*) "  adding vertex to new stack ",nextVertexStackSize, &
                      !      &     cellVertexIndex
                    ENDIF ! (verts%refin_ctrl(cellVertexIndex) == 0)
                  ENDIF ! (cellVertexIndex /= 0)
                ENDDO ! DO j=1,maxCellVertices
                !-------------------------------------
              ENDIF !(level /= endLevel)

            ENDIF ! (cells%refin_ctrl(cellIndex) == 0 .AND. cells%child_id(cellIndex) == childID)
          ENDIF ! (cellIndex /= 0)
        ENDDO ! i = 1,maxVertexCells

      ENDDO ! stackIndex=1, initVertexStackSize

      ! WRITE(*,*) level, " Level is done"
      ! CALL flush(6)

      ! move the new stack to the initStack
      DO i=1,next_vertex_stack_size
        init_vertex_stack(i) = next_vertex_stack(i)
      ENDDO
      init_vertex_stack_size = next_vertex_stack_size

      ! get it to the next level
    ENDDO ! level = startLevel, endLevel, stepLevel


    DEALLOCATE(init_vertex_stack,next_vertex_stack)


    RETURN
  END SUBROUTINE compute_boundary_entities
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  SUBROUTINE get_child_pointers(parent_grid_id, child_grid_id)
    INTEGER, INTENT(in) :: parent_grid_id, child_grid_id

    TYPE(t_grid),   POINTER :: parent_grid, child_grid
    TYPE(t_grid_cells), POINTER :: parent_cells, child_cells
    TYPE(t_grid_edges), POINTER :: parent_edges, child_edges
    TYPE(t_grid_vertices), POINTER :: parent_verts, child_verts
    INTEGER :: no_of_parent_cells, no_of_parent_edges, no_of_parent_verts
    INTEGER :: no_of_child_cells, no_of_child_edges, no_of_child_verts
    INTEGER :: parent_type, index_in_parent

    INTEGER, ALLOCATABLE :: no_of_children(:)
    INTEGER :: j, parent_index, child_index, i_status
!     REAL(wp) :: child_edge_vec(3), parent_edge_vec(3), orientation

    parent_grid     => get_grid(parent_grid_id)
    parent_cells    => parent_grid%cells
    parent_edges    => parent_grid%edges
    parent_verts    => parent_grid%verts
    no_of_parent_cells = parent_cells%no_of_existcells
    no_of_parent_edges = parent_edges%no_of_existedges
    no_of_parent_verts = parent_verts%no_of_existvertices

    child_grid     => get_grid(child_grid_id)
    child_cells    => child_grid%cells
    child_edges    => child_grid%edges
    child_verts    => child_grid%verts
    no_of_child_cells = child_cells%no_of_existcells
    no_of_child_edges = child_edges%no_of_existedges
    no_of_child_verts = child_verts%no_of_existvertices

    IF (child_grid%parent_grid_id /= parent_grid_id) &
      & CALL finish ('get_child_pointers', 'childGrid%parentGridID /= parentID')

    !-------------------------------
    ! fill child cells
    ALLOCATE(no_of_children(no_of_parent_cells),stat=i_status)
    IF (i_status > 0) &
      & CALL finish ('get_child_pointers', 'ALLOCATE(noOfChilds(noOfParentCells))')
    no_of_children(:) = 0
    DO child_index=1, no_of_child_cells
      parent_index = child_cells%parent_index(child_index)
      !     WRITE(*,*) "Cell, child,parent=",childIndex,parentIndex
      IF (parent_index < 1 .or. parent_index > no_of_parent_cells) THEN
        WRITE(0,*) 'parent_index:', parent_index, ' > ', no_of_parent_cells
        CALL finish ('get_child_pointers', 'cell parentIndex out of limits')
      ENDIF

      ! this is only for checking
      no_of_children(parent_index) = no_of_children(parent_index) + 1
      IF (no_of_children(parent_index) > 4) &
        & CALL finish ('get_child_pointers', 'cells: more than 4 children')

      ! the inner child triangle has parent_child_type = parenttype_triangle
      ! the other three have parent_child_type = parenttype_triangle + 1,2,3
      ! By convention of calculating the refine interpolation coeffs
      ! the inner triangle has to be in the 3rd place of the child_index
      index_in_parent = child_cells%parent_child_type(child_index) - parenttype_triangle
      SELECT CASE (index_in_parent)
        CASE (3)
           index_in_parent = 4
        CASE (0)
           index_in_parent = 3
      END SELECT
            
      parent_cells%child_index(parent_index, index_in_parent) = child_index
      IF (parent_cells%child_id(parent_index) /= 0 &
        & .and. parent_cells%child_id(parent_index) /= child_grid_id) &
        & CALL finish ('get_child_pointers', 'more than 1 child for the same cell')

      parent_cells%child_id(parent_index) = child_grid_id

    ENDDO
    DEALLOCATE(no_of_children)

    !-------------------------------
    ! fill child edges
    !   ALLOCATE(noOfChildren(noOfParentEdges),STAT=i_status)
    !   IF (i_status > 0) &
    !      CALL finish ('initTree', 'ALLOCATE(noOfChilds(noOfParentCells))')
    !   noOfChildren(:) = 0
    DO child_index=1, no_of_child_edges
      parent_index = child_edges%parent_index(child_index)
      parent_type  = child_edges%parent_child_type(child_index)
      IF (parent_index < 1 .or. parent_index > no_of_parent_edges) &
        & CALL finish ('get_child_pointers', 'edge parentIndex out of limits')
      ! get the second index j for the parentEdges%child_index
      IF (parent_is_triangle(parent_type)) THEN
        IF (parent_edges%child_index(parent_index, 3) == 0) THEN
          j = 3
        ELSE
          j = 4
        ENDIF
      ELSE
        IF (parent_edges%child_index(parent_index, 1) == 0) THEN
          j = 1
        ELSE
          j = 2
        ENDIF
      ENDIF
      ! compute the orientation relation between normals of the parent/child edge
      ! positive parentEdges%child_index means the normals have the same orientation
      ! negative means they have opposite orientation
      ! NOT used any more
!       CALL gvec2cvec(child_edges%primal_normal(child_index)%v1,   &
!         & child_edges%primal_normal(child_index)%v2,   &
!         & child_edges%center(child_index)%lon,         &
!         & child_edges%center(child_index)%lat,         &
!         & child_edge_vec(1), child_edge_vec(2),        &
!         & child_edge_vec(3))
!       CALL gvec2cvec(parent_edges%primal_normal(parent_index)%v1, &
!         & parent_edges%primal_normal(parent_index)%v2, &
!         & parent_edges%center(parent_index)%lon,       &
!         & parent_edges%center(parent_index)%lat,       &
!         & parent_edge_vec(1), parent_edge_vec(2),      &
!         & parent_edge_vec(3))
!       orientation = DOT_PRODUCT(child_edge_vec,parent_edge_vec)
!
!       IF (orientation >= 0.0_wp) THEN
        parent_edges%child_index(parent_index, j) = child_index
!       ELSE
!         parent_edges%child_index(parent_index, j) = -child_index
!       ENDIF
      IF (parent_edges%child_id(parent_index) /= 0 .and. &
        & parent_edges%child_id(parent_index) /= child_grid_id) &
        & CALL finish ('get_child_pointers', 'more than 1 child for the same edge')
      parent_edges%child_id(parent_index) = child_grid_id

      ! WRITE(*,*) "Edges, child,parent,j,orientation=",parentEdges%child_index(parentIndex, j), &
      !      &     parentIndex,j,orientation
      !     noOfChildren(parentIndex) = noOfChildren(parentIndex) + 1
      !     IF (noOfChildren(parentIndex) > 4) &
      !        CALL finish ('get_child_pointers', 'edges: more than 4 children')

    ENDDO
    !  DEALLOCATE(noOfChildren)

    !-------------------------------
    ! fill child_id for verts
    ! the vertex child is not used (only for consistency at the moment)
    DO child_index=1, no_of_child_verts
      parent_index = child_verts%parent_index(child_index)
      ! WRITE(*,*) parentID, childID, " vertex childIndex,parentIndex=",childIndex,parentIndex
      IF (parent_index > no_of_parent_verts) THEN
        WRITE(*,*) "childIndex,parentIndex=",child_index,parent_index
        CALL finish ('get_child_pointers', 'vertex parentIndex out of limits')
      ENDIF
      IF (parent_index > 0) THEN ! if parentIndex < 0 then the parent is an edge, not a vertex
        parent_verts%child_id(parent_index) = child_grid_id
      ENDIF

    ENDDO

    RETURN

  END SUBROUTINE get_child_pointers
  !-------------------------------------------------------------------------

END MODULE mo_local_grid_hierarchy

