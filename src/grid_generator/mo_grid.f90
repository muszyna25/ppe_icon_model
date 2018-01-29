!>
!!  The grid data structure that is actually used by the numerical subroutines.
!!
!!  The grid data structure that is actually used by the numerical subroutines
!!  is defined in this module. The topologic and geometric information are
!!  computed by the subroutines in <i>m_topology</i> and <i>m_geometry</i>, respectively,
!!   using the internal data representation of the grid generator. The computed quantities
!!  are then loaded onto this <i>grid</i> data structure by the subroutine <i>init_grid.</i>
!!
!! @par Revision History
!! Initial version  by:
!! @par
!! Luis Kornblueh,  MPI-M, Hamburg, January 2004
!! @par
!! Luis Kornblueh, MPI-M, Hamburg, February 2004
!!          - change to compiling and running version
!! Luca Bonaventura,  MPI-M, Hamburg, October 2004
!!          - include documentation and Protex headers
!!          - introduce first version of vertex edge structure
!! @par
!! Luca Bonaventura,  MPI-M, Hamburg, January 2005
!! @par
!!          - include storage of edge structures, modification of data structures
!! @par
!! Luca Bonaventura,  MPI-M, Hamburg, August 2005
!! @par
!!          - different computation of orientations
!! @par
!! Modifications by Th.Heinze, DWD (2006-10-25):
!! - changed index to idx in TYPE declarations of edge_info, vertex_info,
!!   triangle_info, grid_cells, grid_edges and grid_vertices
!! Modification by Thomas Heinze, DWD (2006-10-31):
!! - removed start_level in PUBLIC as it is never used elsewhere and was not
!!   determined here
!! Modification by P. Ripodas, DWD (2007-02):
!! - New array "system_orientation" added to TYPE grid_edges (system
!!   orientation of vectors (v2-v1), (c2-c1), v for vertice and c for cell)
!! - Where needed, memory allocation and deallocation is done for the
!!   new "edges%system_orientation" array. Also "init_grid" is modified
!!   to store the system_orientation values
!! Modifications by Almut Gassmann, MPI-M (2007-03)
!! - reorganized optimization variables
!! Modifications by Almut Gassmann, MPI-M (2007-04)
!! - reorganized destruction of cells, edges and vertices
!! - removed unused edges%neighbor_index
!! - put computation of area_edge to mo_geometry
!! Modifications by Guenther Zaengl, DWD (2008-10-10)
!! - added fields needed for mesh refinement control
!! Modifications by Almut Gassmann, MPI-M (2009-01-15)
!! - cleanings, adding planar option, adding dummy values of grid refinement
!!   so that invoking patch generation could be avoided if not necessary
!! Modifications by Almut Gassmann, MPI-M (2009-01-16)
!! - edge_vert_length and edge_cell_length as new length measures
!! Modifications by Almut Gassmann, MPI-M (2009-03-02)
!! - added beta_spring variable parameter
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
MODULE mo_grid
  !
  !
  !
  USE mo_kind,            ONLY: wp

  USE mo_io_units,        ONLY: filename_max

  USE mo_exception,       ONLY: message_text, message, finish

  USE mo_math_types,      ONLY: t_geographical_coordinates,  &
    & t_tangent_vectors

  USE mo_impl_constants,  ONLY: min_rlcell, max_rlcell, &
    & min_rlvert, max_rlvert, &
    & min_rledge, max_rledge

  IMPLICIT NONE

  PRIVATE


  !--------------------------------------------------------------
  ! Public types
  PUBLIC :: t_grid, t_grid_cells, t_grid_edges, t_grid_vertices

  ! Public subroutines
  PUBLIC :: construct_grid, destruct_grid, construct_grid_ext, &
    & construct_cells, construct_edges, construct_verts, &
    & destruct_cells, destruct_edges, destruct_verts



  ! Public constants
  PUBLIC :: cell_direction_c1_c2,                  &
    & undefined, cut_off_grid, refined_bisection_grid,                &
    & parent_child_identical, parenttype_edge, parenttype_triangle,   &
    & root_node, leaf_node, inner_node,                               &
    & equal_to, greater_than, less_than ,                             &
    & is_land, is_sea
  !--------------------------------------------------------------
  ! definitions of parameters (constants)
  ! Note: 0 is always reserved for undefined value
  INTEGER, PARAMETER ::  undefined = 0

  ! -----------------------------
  ! conditional parameters
  INTEGER, PARAMETER ::  equal_to     = 1
  INTEGER, PARAMETER ::  greater_than = 2
  INTEGER, PARAMETER ::  less_than    = -1

  ! -----------------------------
  ! Sea land parameters
  INTEGER, PARAMETER ::  is_land = 1
  INTEGER, PARAMETER ::  is_sea = -1

  ! -----------------------------
  ! Grid nesting hierarchy parameters
  INTEGER, PARAMETER ::  root_node = 1
  INTEGER, PARAMETER ::  leaf_node = 3
  INTEGER, PARAMETER ::  inner_node = 2
  ! type of parent-child
  INTEGER, PARAMETER ::  parent_child_identical = 1


  INTEGER, PARAMETER ::  parenttype_edge   = 100
  !   For edges:
  !   PARENTTYPE_EDGE+j, j=1,2, signifies the child edge,
  !   that lies beween the j and middle of the parent edge

  INTEGER, PARAMETER ::  parenttype_triangle   = 200
  !   For cells:
  !   PARENTTYPE_TRIANGLE, signifies the child triangle in the middle
  !   that lies beween the j and j+1 edge of the parent  triangle
  !   PARENTTYPE_TRIANGLE+j, signifies the peripheral child triangle,
  !   that lies beween the j and j+1 edge of the parent  triangle

  ! -----------------------------
  ! types of grid creation
  INTEGER, PARAMETER ::  cut_off_grid = 1
  INTEGER, PARAMETER ::  refined_bisection_grid = 2

  ! -----------------------------
  ! types of grid optimization


  ! -----------------------------
  ! parameters for orientation

  ! Edge normal is positive from cell 1 to cell 2, -cell_direction_c1_c2 otherwise
  INTEGER, PARAMETER ::  cell_direction_c1_c2 = 1

  !  END of parameters
  !--------------------------------------------------------------


  !--------------------------------------------------------------
  ! TYPE definitions

  TYPE t_grid_vertices
    INTEGER :: no_of_vertices    ! the number of vertices in this structure
    INTEGER, POINTER :: idx(:)              ! the index of the vertex

    ! geometry part
    TYPE(t_geographical_coordinates), POINTER :: vertex(:) ! the vertex coordinates
    REAL(wp), POINTER :: dual_area(:)                    ! area of dual cell

    ! connectivity part
    INTEGER :: max_connectivity  ! the max number of neigbors
    ! connectivity number = no of neighboring vertices, cells, edges
    INTEGER, POINTER :: no_of_neigbors(:)
    INTEGER, POINTER :: neighbor_index(:,:) ! neigbor idx, 1 to noOfNeigbors
    INTEGER, POINTER :: cell_index(:,:)     ! cell idx, 1 to noOfNeigbors
    INTEGER, POINTER :: edge_index(:,:)     ! edge idx, 1 to noOfNeigbors
    ! +1: from this to neigbor vertex is the same orientation
    ! as the edge unit tangent vector
    ! -1: opposite orientation
    INTEGER, POINTER :: edge_orientation(:,:)


    ! Grid nesting hierarchy part
    INTEGER, POINTER :: parent_index(:)      ! negative values mean that parent is an edge
    INTEGER, POINTER :: parent_child_type(:) !
    INTEGER, POINTER :: refin_ctrl(:)        !
    INTEGER, POINTER :: indlist(:,:)         ! only used by mo_gridrefinement.f90
    INTEGER, POINTER :: start_idx(:,:)       !
    INTEGER, POINTER :: end_idx(:,:)         !
    INTEGER, POINTER :: child_id(:)          !
    INTEGER, POINTER :: phys_id(:)           !
  END TYPE t_grid_vertices
  !--------------------------------------------------------------

  !--------------------------------------------------------------
  TYPE t_grid_edges
    !-----------------------
    !   CONVENTION NOTE: positive normal vector is from cell_index(:,1) to cell_index(:,2)
    !-----------------------
    INTEGER :: no_of_edges
    INTEGER, POINTER :: idx(:)                ! the index of the edge

    ! geometry part

    ! center coordinates. Note: these are the coordinates of
    ! the intersection of the prime and dual edge,
    ! NOT necessarily the middle of the edge
    TYPE(t_geographical_coordinates), POINTER :: center(:)

    TYPE(t_tangent_vectors), POINTER :: primal_normal(:)   ! unit normal vector
    TYPE(t_tangent_vectors), POINTER :: dual_normal(:)     ! dual normal vector (=tangent unit)
    INTEGER, POINTER :: system_orientation(:) !
    !
    REAL(wp), POINTER :: primal_edge_length(:)           ! length of primal edge
    REAL(wp), POINTER :: dual_edge_length(:)             ! length of dual edge
    ! length from edge center to vertices, from 1 to 2
    REAL(wp), POINTER :: edge_vert_length(:,:)
    ! length from edge center to cell centers, from 1 to 2
    REAL(wp), POINTER :: edge_cell_length(:,:)

    ! connectivity part
    INTEGER, POINTER :: cell_index(:,:)       ! indexes of neigboring cells, from 1 to 2
    INTEGER, POINTER :: vertex_index(:,:)     ! indexes of connecting vertices, from 1 to 2
    !
    ! Grid nesting hierarchy part
    INTEGER, POINTER :: refin_ctrl(:)            !
    INTEGER, POINTER :: indlist(:,:)                    ! only used by mo_gridrefinement.f90
    INTEGER, POINTER :: start_idx(:,:)        !
    INTEGER, POINTER :: end_idx(:,:)        !
    INTEGER, POINTER :: parent_index(:)          !
    INTEGER, POINTER :: parent_child_type(:)  !
    INTEGER, POINTER :: child_index(:,:)         !
    INTEGER, POINTER :: child_id(:)            !
    INTEGER, POINTER :: phys_id(:)          !
  END TYPE t_grid_edges
  !--------------------------------------------------------------

  !--------------------------------------------------------------
  TYPE t_grid_cells
    INTEGER :: no_of_cells
    INTEGER, POINTER ::  idx(:)           ! the index of the edge

    ! geometry part
    TYPE(t_geographical_coordinates), POINTER :: center(:) ! the geometric center
    REAL(wp), POINTER :: area(:)             ! cell area

    ! NOTE : negative elevations are below sea level (sea depth)
    REAL(wp), POINTER :: elevation(:)        !
    INTEGER,  POINTER :: sea_land_mask(:)
    REAL(wp)          :: min_sea_depth


    ! connectivity part
    INTEGER :: max_no_of_vertices   ! max number of vertices in one cell
    ! no of cell vertices = no of edges = no of neibgoring cells
    INTEGER, POINTER :: no_of_vertices(:)
    INTEGER, POINTER :: neighbor_index(:,:)   ! neibgoring cells indeces, from 1 to noOfVertices
    INTEGER, POINTER :: vertex_index(:,:)     ! vertex indeces, from 1 to noOfVertices
    INTEGER, POINTER :: edge_index(:,:)       ! edge indeces, from 1 to noOfVertices

    ! edge_orientation defined according to Gauss formula
    ! +1 : the normal to the edge is outwards, -1 : inwards
    ! for LEFT-HAND coord system the same is true for tangent vectors and Stokes form
    ! for RIGHT-HAND coord system it is reveresed for Stokes form: -edge_normal_direction()
    INTEGER, POINTER :: edge_orientation(:,:)

    ! Grid nesting hierarchy part
    INTEGER, POINTER :: parent_index(:)  !
    INTEGER, POINTER :: parent_child_type(:)  !
    INTEGER, POINTER :: child_index(:,:) !
    INTEGER, POINTER :: child_id(:)      !
    INTEGER, POINTER :: phys_id(:)          !
    INTEGER, POINTER :: curr_id(:)       ! ID of current domain (needed in
    ! mo_gridrefinement only)
    INTEGER, POINTER :: refin_ctrl(:)       !
    INTEGER, POINTER :: indlist(:,:)        ! only used by mo_gridrefinement.f90
    INTEGER, POINTER :: start_idx(:,:)      !
    INTEGER, POINTER :: end_idx(:,:)        !
  END TYPE t_grid_cells
  !--------------------------------------------------------------


  !--------------------------------------------------------------
  TYPE t_grid

    INTEGER :: grid_id
    CHARACTER(LEN=filename_max) ::  file_name
    INTEGER :: level    ! grid level ( = number of bisections)

    INTEGER :: ncells   ! no of total cells in the grid
    INTEGER :: nedges   ! no of total edges in the grid
    INTEGER :: nverts   ! no of total vertices in the grid

    ! grid entities
    TYPE(t_grid_cells)    :: cells
    TYPE(t_grid_edges)    :: edges
    TYPE(t_grid_vertices) :: verts


    ! grid nesting hierarchy info
    INTEGER :: patch_id  ! 0=no grid hierarchy
    INTEGER :: root_grid_id
    INTEGER :: parent_grid_id
    INTEGER :: no_of_children
    INTEGER, POINTER :: child_grid_ids(:)

    ! grid status info
    LOGICAL :: is_allocated

  END TYPE t_grid
  !--------------------------------------------------------------

  !-------------------------------------------------------------------------
CONTAINS
  !-------------------------------------------------------------------------



  !>
  !! interface to construct_grid_ext().
  !!
  !!
  !! @par Revision History
  !! Developed  by Leonidas  (2009).
  !!
  SUBROUTINE construct_grid (g,no_of_cells,k_ed,k_ve)
    !-----------------------------------------------------------------------

    INTEGER, INTENT(in) :: no_of_cells,k_ed,k_ve
    TYPE(t_grid), INTENT(inout) :: g

    CALL construct_grid_ext(g,no_of_cells,k_ed,k_ve,3,6)

  END SUBROUTINE construct_grid

  !-------------------------------------------------------------------------
  !
  !
  !>
  !!  Allocates arrays for single grid object.
  !!
  !!
  !! @par Revision History
  !! Developed  by Luca Bonaventura  (2005).
  !! Added maxCellVertices, maxVertexConnect, Leonidas, May 2009
  !!
  SUBROUTINE construct_grid_ext (g,no_of_cells,k_ed,k_ve,max_cell_vertices,max_vertex_connect)
    INTEGER, INTENT(in) :: no_of_cells,k_ed,k_ve,max_cell_vertices,max_vertex_connect
    TYPE(t_grid), INTENT(inout) :: g

    g%cells%max_no_of_vertices = max_cell_vertices
    g%verts%max_connectivity   = max_vertex_connect
    g%ncells  = no_of_cells
    g%nedges  = k_ed
    g%nverts  = k_ve

    g%verts%no_of_vertices          = g%nverts
    g%cells%no_of_cells             = g%ncells
    g%edges%no_of_edges             = g%nedges

    CALL construct_verts(g%verts)
    CALL construct_edges(g%edges)
    CALL construct_cells(g%cells)


  END SUBROUTINE construct_grid_ext
  !--------------------------------------------------------------

  !--------------------------------------------------------------
  SUBROUTINE construct_cells(cells)
    TYPE(t_grid_cells), INTENT(inout) :: cells

    INTEGER :: no_of_cells, max_cell_vertices
    INTEGER :: istat, ist

    no_of_cells       = cells%no_of_cells
    max_cell_vertices = cells%max_no_of_vertices

    ist = 0

    ALLOCATE(cells%idx(no_of_cells),stat=istat)
    ist=ist+istat
    cells%idx(:)=0

    ALLOCATE(cells%parent_index(no_of_cells),stat=istat)
    ist=ist+istat
    cells%parent_index(:)=0
    ALLOCATE(cells%parent_child_type(no_of_cells),stat=istat)
    ist=ist+istat
    cells%parent_child_type(:)=0

    ALLOCATE(cells%child_index(no_of_cells,4),stat=istat)
    ist=ist+istat
    cells%child_index(:,:)=0

    ALLOCATE(cells%child_id(no_of_cells),stat=istat)
    ist=ist+istat
    cells%child_id(:)=0

    ALLOCATE(cells%phys_id(no_of_cells),stat=istat)
    ist=ist+istat
    cells%phys_id(:)=0

    ALLOCATE(cells%curr_id(no_of_cells),stat=istat)
    ist=ist+istat
    cells%curr_id(:)=0

    ALLOCATE(cells%neighbor_index(no_of_cells,max_cell_vertices) ,stat=istat)
    ist=ist+istat
    cells%neighbor_index(:,:)=0

    ALLOCATE(cells%edge_index (no_of_cells,max_cell_vertices),stat=istat)
    ist=ist+istat
    cells%edge_index (:,:)=0

    ALLOCATE(cells%vertex_index(no_of_cells,max_cell_vertices),stat=istat)
    ist=ist+istat
    cells%vertex_index(:,:)=0

    ALLOCATE(cells%edge_orientation(no_of_cells,max_cell_vertices),stat=istat)
    ist=ist+istat
    cells%edge_orientation(:,:)=0

    ALLOCATE(cells%center(no_of_cells),stat=istat)
    ist=ist+istat
    cells%center(:)%lon=0.0_wp
    cells%center(:)%lat=0.0_wp

    ALLOCATE(cells%area(no_of_cells),stat=istat)
    ist=ist+istat
    cells%area(:)=0.0_wp

    ALLOCATE(cells%refin_ctrl(no_of_cells),stat=istat)
    ist=ist+istat
    cells%refin_ctrl(:)=0

    ALLOCATE(cells%indlist(min_rlcell:max_rlcell,2),stat=istat)
    ist=ist+istat
    cells%indlist(:,:)=0

    ALLOCATE(cells%start_idx(min_rlcell:max_rlcell,1),stat=istat)
    ist=ist+istat
    cells%start_idx(:,:)=0

    ALLOCATE(cells%end_idx(min_rlcell:max_rlcell,1),stat=istat)
    ist=ist+istat
    cells%end_idx(:,:)=0

    ALLOCATE(cells%no_of_vertices(no_of_cells),stat=istat)
    ist=ist+istat
    cells%no_of_vertices(:)=0

    ALLOCATE(cells%elevation(no_of_cells),stat=istat)
    ist=ist+istat
    cells%elevation(1:) = 0.0_wp

    ALLOCATE(cells%sea_land_mask(no_of_cells),stat=istat)
    ist=ist+istat
    cells%sea_land_mask(1:) = undefined

    IF (ist>0) THEN
      WRITE (message_text, '(a,i8,a)') &
        & 'construct_cells with ', no_of_cells, ' cells.'
      CALL finish ('construct_cells', TRIM(message_text))
    ENDIF

  END SUBROUTINE construct_cells
  !--------------------------------------------------------------

  !--------------------------------------------------------------
  SUBROUTINE construct_edges(edges)
    TYPE(t_grid_edges), INTENT(inout) :: edges

    INTEGER :: k_ed
    INTEGER :: istat, ist

    k_ed = edges%no_of_edges

    ist = 0

    ALLOCATE(edges%idx(k_ed),stat=istat)
    ist=ist+istat
    edges%idx(:)=0

    ALLOCATE(edges%cell_index (k_ed,2),stat=istat)
    ist=ist+istat
    edges%cell_index(:,:)=0

    ALLOCATE(edges%vertex_index(k_ed,2),stat=istat)
    ist=ist+istat
    edges%vertex_index(:,:)=0

    ALLOCATE(edges%system_orientation(k_ed),stat=istat)
    ist=ist+istat
    edges%system_orientation(:)=0

    ALLOCATE(edges%center(k_ed),stat=istat)
    ist=ist+istat
    edges%center(:)%lat=0.0_wp
    edges%center(:)%lon=0.0_wp

    ALLOCATE(edges%primal_normal(k_ed),stat=istat)
    ist=ist+istat
    edges%primal_normal(:)%v1=0.0_wp
    edges%primal_normal(:)%v2=0.0_wp

    ALLOCATE(edges%dual_normal(k_ed),stat=istat)
    ist=ist+istat
    edges%dual_normal(:)%v1=0.0_wp
    edges%dual_normal(:)%v2=0.0_wp

    ALLOCATE(edges%primal_edge_length(k_ed),stat=istat)
    ist=ist+istat
    edges%primal_edge_length(:)=0.0_wp

    ALLOCATE(edges%dual_edge_length(k_ed),stat=istat)
    ist=ist+istat
    edges%dual_edge_length(:)=0.0_wp

    ALLOCATE(edges%edge_vert_length(k_ed,2),stat=istat)
    ist=ist+istat
    edges%edge_vert_length(:,:)=0.0_wp

    ALLOCATE(edges%edge_cell_length(k_ed,2),stat=istat)
    ist=ist+istat
    edges%edge_cell_length(:,:)=0.0_wp

    ALLOCATE(edges%refin_ctrl(k_ed),stat=istat)
    ist=ist+istat
    edges%refin_ctrl(:)=0

    ALLOCATE(edges%indlist(min_rledge:max_rledge,2),stat=istat)
    ist=ist+istat
    edges%indlist(:,:)=0

    ALLOCATE(edges%start_idx(min_rledge:max_rledge,1),stat=istat)
    ist=ist+istat
    edges%start_idx(:,:)=0

    ALLOCATE(edges%end_idx(min_rledge:max_rledge,1),stat=istat)
    ist=ist+istat
    edges%end_idx(:,:)=0

    ALLOCATE(edges%parent_index(k_ed),stat=istat)
    ist=ist+istat
    edges%parent_index(:)=0
    ALLOCATE(edges%parent_child_type(k_ed),stat=istat)
    ist=ist+istat
    edges%parent_child_type(:)=0

    ALLOCATE(edges%child_index(k_ed,4),stat=istat)
    ist=ist+istat
    edges%child_index(:,:)=0

    ALLOCATE(edges%child_id(k_ed),stat=istat)
    ist=ist+istat
    edges%child_id(:)=0

    ALLOCATE(edges%phys_id(k_ed),stat=istat)
    ist=ist+istat
    edges%phys_id(:)=0


    IF (ist>0) THEN
      WRITE (message_text, '(a,i8,a)') &
        & 'construct_edges with ', k_ed, ' edges.'
      CALL finish ('construct_edges', TRIM(message_text))
    ENDIF

  END SUBROUTINE construct_edges
  !--------------------------------------------------------------

  !--------------------------------------------------------------
  SUBROUTINE construct_verts(verts)
    TYPE(t_grid_vertices), INTENT(inout) :: verts

    INTEGER :: k_ve,max_vertex_connect
    INTEGER :: istat, ist

    k_ve               = verts%no_of_vertices
    max_vertex_connect = verts%max_connectivity

    ist = 0

    ALLOCATE(verts%idx(k_ve),stat=istat)
    ist=ist+istat
    verts%idx(:)=0

    ALLOCATE(verts%parent_index(k_ve),stat=istat)
    ist=ist+istat
    verts%parent_index(:)=0
    ALLOCATE(verts%parent_child_type(k_ve),stat=istat)
    ist=ist+istat
    verts%parent_child_type(:)=0

    ALLOCATE(verts%no_of_neigbors(k_ve),stat=istat)
    ist=ist+istat
    verts%no_of_neigbors(:)=0

    ALLOCATE(verts%neighbor_index(k_ve,max_vertex_connect) ,stat=istat)
    ist=ist+istat
    verts%neighbor_index(:,:)=0

    ALLOCATE(verts%cell_index (k_ve,max_vertex_connect),stat=istat)
    ist=ist+istat
    verts%cell_index (:,:)=0

    ALLOCATE(verts%edge_index (k_ve,max_vertex_connect),stat=istat)
    ist=ist+istat
    verts%edge_index (:,:)=0

    ALLOCATE(verts%edge_orientation(k_ve,max_vertex_connect),stat=istat)
    ist=ist+istat
    verts%edge_orientation(:,:)=0

    ALLOCATE(verts%dual_area(k_ve),stat=istat)
    ist=ist+istat
    verts%dual_area(:)=0.0_wp

    ALLOCATE(verts%vertex(k_ve),stat=istat)
    ist=ist+istat
    verts%vertex(:)%lon=0.0_wp
    verts%vertex(:)%lat=0.0_wp

    ALLOCATE(verts%refin_ctrl(k_ve),stat=istat)
    ist=ist+istat
    verts%refin_ctrl(:)=0

    ALLOCATE(verts%indlist(min_rlvert:max_rlvert,2),stat=istat)
    ist=ist+istat
    verts%indlist(:,:)=0

    ALLOCATE(verts%start_idx(min_rlvert:max_rlvert,1),stat=istat)
    ist=ist+istat
    verts%start_idx(:,:)=0

    ALLOCATE(verts%end_idx(min_rlvert:max_rlvert,1),stat=istat)
    ist=ist+istat
    verts%end_idx(:,:)=0

    ALLOCATE(verts%child_id(k_ve),stat=istat)
    ist=ist+istat
    verts%child_id(:)=0

    ALLOCATE(verts%phys_id(k_ve),stat=istat)
    ist=ist+istat
    verts%phys_id(:)=0

    IF (ist>0) THEN
      WRITE (message_text, '(a,i8,a)') &
        & 'construct_verts with ', k_ve, ' vertices.'
      CALL finish ('construct_verts', TRIM(message_text))
    ENDIF

  END SUBROUTINE construct_verts

  !-------------------------------------------------------------------------
  !>
  !!  Deallocates arrays for single grid object.
  !!
  !! @par Revision History
  !! Developed  by Luca Bonaventura  (2005).
  !!
  SUBROUTINE destruct_grid(g)
    TYPE(t_grid), INTENT(inout) :: g

    CALL destruct_cells(g%cells)
    CALL destruct_edges(g%edges)
    CALL destruct_verts(g%verts)

  END SUBROUTINE destruct_grid
  !-----------------------------------------------------------------------


  !-----------------------------------------------------------------------
  SUBROUTINE destruct_cells(cells)
    TYPE(t_grid_cells), INTENT(inout) :: cells

    INTEGER :: istat, ist

    ist=0
    DEALLOCATE(cells%idx,stat=istat)
    ist=ist+istat
    DEALLOCATE(cells%parent_index,stat=istat)
    ist=ist+istat
    DEALLOCATE(cells%parent_child_type,stat=istat)
    ist=ist+istat
    DEALLOCATE(cells%child_index,stat=istat)
    ist=ist+istat
    DEALLOCATE(cells%child_id,stat=istat)
    ist=ist+istat
    DEALLOCATE(cells%phys_id,stat=istat)
    ist=ist+istat
    DEALLOCATE(cells%curr_id,stat=istat)
    ist=ist+istat
    DEALLOCATE(cells%no_of_vertices,stat=istat)
    ist=ist+istat
    DEALLOCATE(cells%neighbor_index ,stat=istat)
    ist=ist+istat
    DEALLOCATE(cells%edge_index,stat=istat)
    ist=ist+istat
    DEALLOCATE(cells%vertex_index,stat=istat)
    ist=ist+istat
    DEALLOCATE(cells%edge_orientation,stat=istat)
    ist=ist+istat
    DEALLOCATE(cells%center,stat=istat)
    ist=ist+istat
    DEALLOCATE(cells%area,stat=istat)
    ist=ist+istat
    DEALLOCATE(cells%refin_ctrl,stat=istat)
    ist=ist+istat
    DEALLOCATE(cells%indlist,stat=istat)
    ist=ist+istat
    DEALLOCATE(cells%start_idx,stat=istat)
    ist=ist+istat
    DEALLOCATE(cells%end_idx,stat=istat)
    ist=ist+istat
    DEALLOCATE(cells%elevation,stat=istat)
    ist=ist+istat
    DEALLOCATE(cells%sea_land_mask,stat=istat)
    ist=ist+istat

    IF (ist>0) THEN
      WRITE (message_text, '(a)') 'Deallocate grid.'
      CALL message ('destruct_cells', 'Deallocate grid.')
    ENDIF

  END SUBROUTINE destruct_cells
  !-----------------------------------------------------------------------


  !-----------------------------------------------------------------------
  SUBROUTINE destruct_edges(edges)
    TYPE(t_grid_edges), INTENT(inout) :: edges

    INTEGER :: istat, ist

    ist=0

    DEALLOCATE(edges%idx,stat=istat)
    ist=ist+istat
    DEALLOCATE(edges%cell_index ,stat=istat)
    ist=ist+istat
    DEALLOCATE(edges%vertex_index,stat=istat)
    ist=ist+istat
    DEALLOCATE(edges%system_orientation,stat=istat)
    ist=ist+istat
    DEALLOCATE(edges%center,stat=istat)
    ist=ist+istat
    DEALLOCATE(edges%primal_normal,stat=istat)
    ist=ist+istat
    DEALLOCATE(edges%dual_normal,stat=istat)
    ist=ist+istat
    DEALLOCATE(edges%primal_edge_length,stat=istat)
    ist=ist+istat
    DEALLOCATE(edges%dual_edge_length,stat=istat)
    ist=ist+istat
    DEALLOCATE(edges%edge_vert_length,stat=istat)
    ist=ist+istat
    DEALLOCATE(edges%edge_cell_length,stat=istat)
    ist=ist+istat
    DEALLOCATE(edges%refin_ctrl,stat=istat)
    ist=ist+istat
    DEALLOCATE(edges%indlist,stat=istat)
    ist=ist+istat
    DEALLOCATE(edges%start_idx,stat=istat)
    ist=ist+istat
    DEALLOCATE(edges%end_idx,stat=istat)
    ist=ist+istat
    DEALLOCATE(edges%parent_index,stat=istat)
    ist=ist+istat
    DEALLOCATE(edges%parent_child_type,stat=istat)
    ist=ist+istat
    DEALLOCATE(edges%child_index,stat=istat)
    ist=ist+istat
    DEALLOCATE(edges%child_id,stat=istat)
    ist=ist+istat
    DEALLOCATE(edges%phys_id,stat=istat)
    ist=ist+istat

    IF (ist>0) THEN
      WRITE (message_text, '(a)') 'Deallocate grid.'
      CALL message ('destruct_edges', 'Deallocate grid.')
    ENDIF

  END SUBROUTINE destruct_edges
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  SUBROUTINE destruct_verts(verts)
    TYPE(t_grid_vertices), INTENT(inout) :: verts

    INTEGER :: istat, ist

    ist=0

    DEALLOCATE(verts%idx,stat=istat)
    ist=ist+istat
    DEALLOCATE(verts%parent_index,stat=istat)
    ist=ist+istat
    DEALLOCATE(verts%parent_child_type,stat=istat)
    ist=ist+istat
    DEALLOCATE(verts%no_of_neigbors,stat=istat)
    ist=ist+istat
    DEALLOCATE(verts%neighbor_index,stat=istat)
    ist=ist+istat
    DEALLOCATE(verts%cell_index,stat=istat)
    ist=ist+istat
    DEALLOCATE(verts%edge_index,stat=istat)
    ist=ist+istat
    DEALLOCATE(verts%edge_orientation,stat=istat)
    ist=ist+istat
    DEALLOCATE(verts%dual_area,stat=istat)
    ist=ist+istat
    DEALLOCATE(verts%vertex,stat=istat)
    ist=ist+istat
    DEALLOCATE(verts%refin_ctrl,stat=istat)
    ist=ist+istat
    DEALLOCATE(verts%indlist,stat=istat)
    ist=ist+istat
    DEALLOCATE(verts%start_idx,stat=istat)
    ist=ist+istat
    DEALLOCATE(verts%end_idx,stat=istat)
    ist=ist+istat
    DEALLOCATE(verts%child_id,stat=istat)
    ist=ist+istat
    DEALLOCATE(verts%phys_id,stat=istat)
    ist=ist+istat

    IF (ist>0) THEN
      WRITE (message_text, '(a)') 'Deallocate grid.'
      CALL message ('destruct_verts', 'Deallocate grid.')
    ENDIF

  END SUBROUTINE destruct_verts
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------

END MODULE mo_grid
