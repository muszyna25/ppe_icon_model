!>
!! Defines the grid data structure used by the local grid refinement.
!! Basic "object styled" methods are provided through this module.
!!
!!
!! @par Revision History
!! Current version by Leonidas Linardakis,  MPI-M, January 2010.
!! @par
!! Based on the previous mo_grid.f90 developed by:
!! Luis Kornblueh, MPI-M.
!! Luca Bonaventura, MPI-M.
!! Th.Heinze, DWD.
!! P. Ripodas, DWD.
!! Almut Gassmann, MPI-M.
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
MODULE mo_local_grid
#include "grid_definitions.inc"

  USE mo_kind,            ONLY: wp
  USE mo_io_units,        ONLY: filename_max
  USE mo_exception,       ONLY: message_text, message, finish
  USE mo_base_geometry,   ONLY: t_geographical_coordinates, t_cartesian_coordinates, &
    & t_tangent_vectors
  USE mo_impl_constants,  ONLY: &
    & min_rlcell, max_rlcell, min_rlcell_int, &
    & min_rlvert, max_rlvert, min_rlvert_int,  &
    & min_rledge, max_rledge, min_rledge_int

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: version = '$Id$'

  !--------------------------------------------------------------
  ! Public types
  PUBLIC :: t_grid, t_grid_cells, t_grid_edges, t_grid_vertices, t_vertical_ocean_structure, &
    & t_integer_list, t_float_list

  ! Public subroutines
  PUBLIC :: construct_cells, construct_edges, construct_verts, &
    & destruct_cells, destruct_edges, destruct_verts,          &
    & allocate_vertical_columns, deallocate_vertical,          &
    & allocate_vertical_zlayers, allocate_vertical,            &
    & set_nest_defaultindexes, get_next_cell_edge_vertex

  ! Public object-oriented methods
  PUBLIC :: detsruct_grid_objects, new_grid, delete_grid, get_grid,  &
    & allocate_grid_object, set_grid_filename, set_grid_creation,    &
    & get_cells, get_vertices, replace_grid,  copy_grid,             &
    & grid_set_exist_eq_allocated, grid_set_allocated_eq_exist,      &
    & get_edges, set_grid_vertical_structure,grid_get_parent_pointers,&
    & get_grid_child_id, get_grid_no_of_children,grid_is_filled,     &
    & grid_set_parents_from, set_grid_parent_id,                     &
    & set_no_of_subgrids, set_start_subgrids, grid_set_sea_depth,    &
    & grid_set_min_sea_depth, set_grid_level, get_grid_level,        &
    & set_grid_netcdf_flags, get_number_of_vertices

  PUBLIC :: print_grid_cell, print_grid_edge, print_grid_vertex

  ! Public parameters
  PUBLIC :: max_no_of_grid_objects, cell_direction_c1_c2,             &
    & undefined, cut_off_grid, refined_bisection_grid,                &
    & parent_child_identical, parenttype_edge, parenttype_triangle,   &
    & root_node, leaf_node, inner_node,                               &
    & equal_to, greater_than, less_than,                              &
    & land_inner, land_boundary, sea_inner, sea_boundary,             &
    & linear_interpolation, min_interpolation, max_interpolation,     &
    & parents_from_idpointers, parents_from_parentpointers,           &
    & netcdf_CF_1_1_convention, sphere_geometry, torus_geometry
  
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
  ! interpolation methods parameters
  INTEGER, PARAMETER :: linear_interpolation = 1
  INTEGER, PARAMETER :: min_interpolation    = 2
  INTEGER, PARAMETER :: max_interpolation    = 3

  ! -----------------------------
  ! Sea land parameters
  INTEGER, PARAMETER ::  land_inner    = 2
  INTEGER, PARAMETER ::  land_boundary = 1
  INTEGER, PARAMETER ::  sea_boundary  = -1
  INTEGER, PARAMETER ::  sea_inner     = -2

  ! -----------------------------
  ! Grid nesting hierarchy parameters
  INTEGER, PARAMETER ::  root_node = 1
  INTEGER, PARAMETER ::  leaf_node = 3
  INTEGER, PARAMETER ::  inner_node = 2

  INTEGER, PARAMETER ::  parents_from_idpointers = 1
  INTEGER, PARAMETER ::  parents_from_parentpointers = 2

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
  ! types of grid geometries
  INTEGER, PARAMETER ::  sphere_geometry = 1
  INTEGER, PARAMETER ::  torus_geometry  = 2

  ! -----------------------------
  ! types of grid optimization


  ! -----------------------------
  ! parameters for orientation

  ! Edge normal is positive from cell 1 to cell 2, -cell_direction_c1_c2 otherwise
  INTEGER, PARAMETER ::  cell_direction_c1_c2 = 1

  !--------------------------------------------------------------
  ! other parameters
  INTEGER, PARAMETER ::  netcdf_CF_1_1_convention = 1


  !  END of parameters
  !--------------------------------------------------------------


  !--------------------------------------------------------------
  ! TYPE definitions
  !> Holds a list of initegers
  TYPE t_integer_list
    INTEGER :: list_size
    INTEGER, POINTER :: value(:)
  END TYPE t_integer_list

  !> Holds a list of reals
  TYPE t_float_list
    INTEGER :: list_size
    REAL(wp), POINTER :: value(:)
  END TYPE t_float_list
  !--------------------------------------------------------------

  !--------------------------------------------------------------
  !> The vertical structure for ocean grid
  ! Includes the layers and columns information
  TYPE t_vertical_ocean_structure
    ! the z layers info
    INTEGER :: no_of_zlayers
    REAL(wp), POINTER :: layer_thicknes(:)
    REAL(wp), POINTER :: layer_bed(:)
    REAL(wp), POINTER :: layer_middle(:)
    ! the grid column info
    INTEGER :: no_of_columns
    INTEGER, POINTER :: column_cell_index(:,:) ! (no_of_zlayers,no_of_columns)
    INTEGER, POINTER :: column_size(:) ! (no_of_columns)
  END TYPE t_vertical_ocean_structure
  !--------------------------------------------------------------

  !--------------------------------------------------------------
  !> The vertices structure of the grid
  TYPE t_grid_vertices
    !> The number of allocated vertices in this structure
    INTEGER :: no_of_allocatedvertices
    !> The number of actual vertices in this structure
    INTEGER :: no_of_existvertices
    !> The index of this vertex
    INTEGER, POINTER :: idx(:)

    ! connectivity part
    !> The maximum number of neigbors (vertices, edges)
    INTEGER :: max_connectivity
    ! connectivity number = no of neighboring vertices, edges.
    !> The number of neigbors (vertices, edges)
    INTEGER, POINTER :: no_of_neigbors(:)
    !> The neigboring vertices indexes
    INTEGER, POINTER :: a_neighbor_index(:,:)
    !> The neigboring cells indexes
    INTEGER, POINTER :: a_cell_index(:,:)
     !> The neigboring edges indexes
    INTEGER, POINTER :: a_edge_index(:,:)
    !> The orientation of the tangent vector to the neigboring edges.
    !! +1: from this to the neigbor vertex has the same orientation as the tangent vector.
    !! -1: opposite orientation
    INTEGER, POINTER :: a_edge_orientation(:,:)

    ! geometry part
    !> The vertex geographical coordinates
    TYPE(t_geographical_coordinates), POINTER :: vertex(:)
    !> The vertex cartesian (3D) coordinates on the unit sphere
    TYPE(t_cartesian_coordinates), POINTER ::  cartesian(:)
    !> The area of dual cell
    REAL(wp), POINTER :: dual_area(:)


    ! Grid nesting hierarchy part
    INTEGER, POINTER :: parent_index(:)      ! negative values mean that parent is an edge
    INTEGER, POINTER :: parent_child_type(:) !
    INTEGER, POINTER :: refin_ctrl(:)        !
    ! INTEGER, POINTER :: indlist(:,:)         ! only used by mo_gridrefinement.f90
    INTEGER, POINTER :: start_idx(:,:)       !
    INTEGER, POINTER :: end_idx(:,:)         !
    INTEGER, POINTER :: child_id(:)          !
  END TYPE t_grid_vertices
  !--------------------------------------------------------------

  !--------------------------------------------------------------
  !> The edges structure of the grid
  TYPE t_grid_edges
    !-----------------------
    !   CONVENTION NOTE: positive normal vector is from cell_index(:,1) to cell_index(:,2)
    !-----------------------
    !> The number of allocated edges
    INTEGER :: no_of_allocatededges
    !> The number of actual edges
    INTEGER :: no_of_existedges
    !> The index of this edge
    INTEGER, POINTER :: idx(:)    ! the index of this edge

    ! connectivity part
    INTEGER, POINTER :: a_cell_index(:,:)       ! indexes of neigboring cells, from 1 to 2
    INTEGER, POINTER :: a_vertex_index(:,:)     ! indexes of connecting vertices, from 1 to 2
    ! length from edge center to vertices, from 1 to 2
    REAL(wp), POINTER :: a_edge_vert_length(:,:)
    ! length from edge center to cell centers, from 1 to 2
    REAL(wp), POINTER :: a_edge_cell_length(:,:)

    ! geometry part
    !> The geographical coordinates of the edge center.
    !! These are the coordinates of the intersection of the prime and dual edge,
    !! It is not necessarily the middle of the edge.
    TYPE(t_geographical_coordinates), POINTER :: center(:)
    !> The center in cartesian (3D) coordinates on the unit sphere
    TYPE(t_cartesian_coordinates), POINTER ::  cartesian_center(:)
    !> The normal to the edge unit vector in local coordinates
    !! (on the tangent plane to the sphere).
    !! It always points from vertex_index(,1) to vertex_index(,2).
    TYPE(t_tangent_vectors), POINTER :: primal_normal(:)
    !> The normal in cartesian coordinates
    TYPE(t_cartesian_coordinates), POINTER :: cartesian_primal_normal(:)
    !> The tangent to the edge unit vector in local coordinates
    !! The (primal_normal,dual_normal) forms a left-handed system
    TYPE(t_tangent_vectors), POINTER :: dual_normal(:)
    !> The dual_normal in cartesian coordinates
    TYPE(t_cartesian_coordinates), POINTER :: cartesian_dual_normal(:)
    !> Is equal to +1 if the vector from vertex_index(,1) vertex_index(,2)
    !!  has the same orientation as the tangent vector.
    !! Is equal to -1 otherwise.
    INTEGER, POINTER :: system_orientation(:) !
    !> The length of the primal edge
    REAL(wp), POINTER :: primal_edge_length(:)
    !> The length of the dual edge
    REAL(wp), POINTER :: dual_edge_length(:)
    !> area of the quadrilateral formed by two cells adjacent to the edge
    REAL(wp), ALLOCATABLE :: quad_area(:)

    ! Ocean part
    !> The edge topography. Negative values signify sea depth.
    REAL(wp), POINTER :: elevation(:)
    !> The number of active vertical layers for each edge
    INTEGER,  POINTER :: no_of_zlayers(:)
    !> Holds the SEA,LAND information (see the parameters in this module).
    INTEGER,  POINTER :: sea_land_mask(:)

    ! Implicit blocking part
    INTEGER :: no_of_blocks
    INTEGER,  POINTER :: block_end_idx(:)

    ! Grid nesting hierarchy part
    INTEGER,  POINTER :: subgrid_id(:)
    INTEGER, POINTER :: refin_ctrl(:)            !
    ! INTEGER, POINTER :: indlist(:,:)                    ! only used by mo_gridrefinement.f90
    INTEGER, POINTER :: start_idx(:,:)        !
    INTEGER, POINTER :: end_idx(:,:)        !
    INTEGER, POINTER :: parent_index(:)          !
    INTEGER, POINTER :: parent_child_type(:)  !
    INTEGER, POINTER :: child_index(:,:)         !
    INTEGER, POINTER :: child_id(:)            !
  END TYPE t_grid_edges
  !--------------------------------------------------------------

  !--------------------------------------------------------------
  !> The cells structure of the grid
  TYPE t_grid_cells
    !> The number of allocated cells
    INTEGER :: no_of_allocatedcells
    !> The number of actual cells
    INTEGER :: no_of_existcells
    !> The index of this cell
    INTEGER, POINTER ::  idx(:)     ! the index of this cell

    ! connectivity part
    !> The max number of vertices in one cell
    INTEGER :: max_no_of_vertices
    ! no of cell vertices = no of edges = no of neibgoring cells
    !> The actual number of vertices for each cell
    INTEGER, POINTER :: no_of_vertices(:)
    !> The indexes of the cell vetrices.
    INTEGER, POINTER :: a_vertex_index(:,:)
    !> The indexes of the cell edges.
    INTEGER, POINTER :: a_edge_index(:,:)
    !> The indexes of the neighboring cells.
    INTEGER, POINTER :: a_neighbor_index(:,:)
    !> The edge_orientation defined according to Gauss formula.
    !! Is equal to +1 if the normal to the edge is outwards.
    !! Is -1 if it's inwards.
    !! For the left-hand coordinate system (which is in use),
    !! the same is true for tangent vectors and Stokes formula.
    INTEGER, POINTER :: a_edge_orientation(:,:)

    ! geometry part
    !> The geographical coordinates of the geometric center (cirmucenter) of the cell.
    TYPE(t_geographical_coordinates), POINTER :: center(:) ! the geometric center
    !> The center in cartesian (3D) coordinates on the unit sphere
    TYPE(t_cartesian_coordinates), POINTER ::  cartesian_center(:)
    !> The area of the cell.
    REAL(wp), POINTER :: area(:)             ! cell area

    ! Ocean part
    ! NOTE : negative elevations are below sea level (sea depth)
    !> The depth below which the cells are considered to be sea.
    REAL(wp)          :: min_sea_depth
    !> The cell topography. Negative values signify sea depth.
    REAL(wp), POINTER :: elevation(:)
    !> The number of active vertical layers for each cell
    INTEGER,  POINTER :: no_of_zlayers(:)
    !> Holds the SEA,LAND information (see the parameters in this module).
    INTEGER,  POINTER :: sea_land_mask(:)

    ! Implicit blocking part
    INTEGER :: no_of_blocks
    INTEGER,  POINTER :: block_end_idx(:)

    ! Grid nesting hierarchy part
    INTEGER, POINTER :: subgrid_id(:)
    INTEGER, POINTER :: parent_index(:)  !
    INTEGER, POINTER :: parent_child_type(:)  !
    INTEGER, POINTER :: child_index(:,:) !
    INTEGER, POINTER :: child_id(:)      !
    INTEGER, POINTER :: refin_ctrl(:)       !
    ! INTEGER, POINTER :: indlist(:,:)        ! only used by mo_gridrefinement.f90
    INTEGER, POINTER :: start_idx(:,:)      !
    INTEGER, POINTER :: end_idx(:,:)        !
  END TYPE t_grid_cells
  !--------------------------------------------------------------


  !--------------------------------------------------------------
  !> The grid structure
  TYPE t_grid
    !grid characteristics
    !> Filename of this grid
    CHARACTER(LEN=filename_max) ::  file_name
    !> Output filename for this grid
    CHARACTER(LEN=filename_max) ::  out_file_name
    !> flags for the netcdf files
    INTEGER :: netcdf_flags
    !> The creation process of the grid (cut_off, refined, etc, see parameters).
    INTEGER :: grid_creation_process
    !> The grid optimization procces
    INTEGER :: grid_optimization_process
    !> The grid domain geometry
    INTEGER :: grid_geometry

    ! grid entities
    !> The number of cells in the grid
    INTEGER :: ncells
    !> The number of edges in the grid
    INTEGER :: nedges
    !> The number of vertices in the grid
    INTEGER :: nverts

    !> The cells in the grid
    TYPE(t_grid_cells)    :: cells
    !> The edges in the grid
    TYPE(t_grid_edges)    :: edges
    !> The vertices in the grid
    TYPE(t_grid_vertices) :: verts

    INTEGER :: no_of_subgrids
    INTEGER :: start_subgrid_id

    ! grid vertical parameters
    !> The ocean vertical structure (layers and columns) of the grid
    TYPE(t_vertical_ocean_structure), POINTER :: vertical_structure

    ! grid object parameters
    !> The ID of this grid object
    INTEGER :: grid_obj_id
    !> True if the grid entities (cells,edges,vertices) are allocated.
    LOGICAL :: is_allocated
    !> True if the grid entities (cells,edges,vertices) are filled.
    LOGICAL :: is_filled

    ! grid nesting hierarchy info
    INTEGER :: patch_id  ! less than 0= no grid hierarchy
    INTEGER :: hierarchy_node_type  ! 0=no grid hierarchy, 1= root,  2=inter, 3=leaf
    INTEGER :: parent_grid_id
    !> The grid optimization procces
    INTEGER :: parents_from
    INTEGER :: no_of_children
    INTEGER, POINTER :: child_grid_ids(:)

    ! other
    INTEGER :: level    ! grid level ( = number of bisections)

  END TYPE t_grid
  !--------------------------------------------------------------

  !--------------------------------------------------------------
  ! The objects data
  !--------------------------------------------------------------
  !> The maximum number of grid objects.
  INTEGER, PARAMETER ::  max_no_of_grid_objects = 50
  !> The array of the grid objects.
  TYPE(t_grid), ALLOCATABLE, TARGET :: grid_object_list(:)
  TYPE t_grid_pointer_type
    TYPE(t_grid), POINTER ::  pnt
  END TYPE t_grid_pointer_type
  TYPE(t_grid_pointer_type), ALLOCATABLE :: grid_pnt_array(:)
  !> The number of allocated grid objects.
  INTEGER :: no_of_allocated_grids = 0        ! the size of the grid objects array
  !> The number of actual active grid objects.
  INTEGER :: active_grids
  !> The maximum id of the active grids
  INTEGER :: max_active_grids
  !> True if the grid object is active.
  LOGICAL, ALLOCATABLE :: grid_is_active(:)
  !> The index where the grid actually exists.
#define get_grid_object(id)  grid_pnt_array(id)%pnt
#define get_verts_object(id) grid_pnt_array(id)%pnt%verts
#define get_edges_object(id) grid_pnt_array(id)%pnt%edges
#define get_cell_object(id)  grid_pnt_array(id)%pnt%cells
  !-------------------------------------------------------------------------
CONTAINS
  !-------------------------------------------------------------------------
  ! Public Methods
  !-------------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  !! Creates a new grid object and returns its id.
  !! The grid arrays are not allocated.
  INTEGER FUNCTION new_grid()

    INTEGER :: i

    ! check if the grid objects have been construced
    IF (no_of_allocated_grids == 0) THEN
      CALL construct_grid_objects()
    ENDIF

    IF (max_active_grids /= active_grids) THEN
      ! we have a hole (inactive grid) in the list of max_active_grids
      DO i = 1, max_active_grids
        IF (.not. grid_is_active(i)) THEN
          new_grid = i
          EXIT
        ENDIF
      ENDDO
    ELSE
      ! add a new grid object
      IF (max_active_grids >= no_of_allocated_grids) THEN
        CALL finish('new_grid', 'exceeded no_of_allocated_grids');
      ENDIF
      max_active_grids = max_active_grids + 1
      new_grid = max_active_grids
    ENDIF

    active_grids = active_grids + 1
    CALL init_grid_object(new_grid)

  END FUNCTION new_grid
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  !! Deletes the grid object.
  !! The grid arrays are deallocated.
  !! The grid_id is invalid after the delete call.
  !! It maybe associated with another grid in the future.
  SUBROUTINE delete_grid(grid_id)
    INTEGER, INTENT(in) :: grid_id

    CALL check_active_grid_id(grid_id)

    CALL deallocate_grid_object(grid_id)

    grid_is_active(grid_id) = .false.
    active_grids = active_grids - 1

  END SUBROUTINE delete_grid
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  !! Returns a pointer to the grid structure associated with the grid_id
  FUNCTION get_grid(grid_id) result(grid_obj)
    INTEGER, INTENT(in) :: grid_id

    TYPE(t_grid), POINTER :: grid_obj

    CALL check_active_grid_id(grid_id)

    grid_obj => get_grid_object(grid_id)

  END FUNCTION get_grid
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  !! Returns a pointer to the cells structure associated with the grid_id
  FUNCTION get_cells(grid_id) result(cell_obj)
    INTEGER, INTENT(in) :: grid_id

    TYPE(t_grid_cells), POINTER :: cell_obj

    CALL check_active_grid_id(grid_id)

    cell_obj => get_grid_object(grid_id)%cells

  END FUNCTION get_cells
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  !! Returns a pointer to the edges structure associated with the grid_id
  FUNCTION get_edges(grid_id) result(edge_obj)
    INTEGER, INTENT(in) :: grid_id

    TYPE(t_grid_edges), POINTER :: edge_obj

    CALL check_active_grid_id(grid_id)

    edge_obj => get_grid_object(grid_id)%edges

  END FUNCTION get_edges
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  !! Returns a pointer to the vertices structure associated with the grid_id
  FUNCTION get_vertices(grid_id) result(vert_obj)
    INTEGER, INTENT(in) :: grid_id

    TYPE(t_grid_vertices), POINTER :: vert_obj

    CALL check_active_grid_id(grid_id)

    vert_obj => get_grid_object(grid_id)%verts

  END FUNCTION get_vertices
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  !! Returns a pointer to the vertices structure associated with the grid_id
  INTEGER FUNCTION get_number_of_vertices(grid_id)
    INTEGER, INTENT(in) :: grid_id

    CALL check_active_grid_id(grid_id)

    get_number_of_vertices = get_grid_object(grid_id)%verts%no_of_existvertices

  END FUNCTION get_number_of_vertices
  !-----------------------------------------------------------------------


  !-----------------------------------------------------------------------
  LOGICAL FUNCTION grid_is_filled(grid_id)
    INTEGER, INTENT(in) ::grid_id

    grid_is_filled = &
      & get_grid_object(grid_id)%is_filled

  END FUNCTION grid_is_filled
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  INTEGER FUNCTION get_grid_no_of_children(grid_id)
    INTEGER, INTENT(in) ::grid_id

    get_grid_no_of_children = &
      & get_grid_object(grid_id)%no_of_children

  END FUNCTION get_grid_no_of_children
  !-----------------------------------------------------------------------


  !-----------------------------------------------------------------------
  INTEGER FUNCTION get_grid_child_id(grid_id,child_no)
    INTEGER, INTENT(in) ::grid_id,child_no

    get_grid_child_id = &
      & get_grid_object(grid_id)%child_grid_ids(child_no)

  END FUNCTION get_grid_child_id
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  !! Allocates the gird arrays associated with the grid_id.
  !! The array sizes are determined by the variables of the grid structure:
  !! grid\%ncells,grid\%nedges,grid\%nverts
  !! grid\%cells\%max_no_of_vertices, grid\%verts\%max_connectivity
  !! The above fields in the grid object must have been filled
  !! before calling allocate_grid_object
  SUBROUTINE allocate_grid_object(grid_id)
    INTEGER, INTENT(in) :: grid_id

    TYPE(t_grid), POINTER :: grid_obj

    CALL check_active_grid_id(grid_id)
    grid_obj => get_grid_object(grid_id)

    IF (grid_obj%is_allocated) &
      & CALL deallocate_grid_object(grid_id)

    grid_obj%cells%no_of_allocatedcells = grid_obj%ncells
    grid_obj%edges%no_of_allocatededges = grid_obj%nedges
    grid_obj%verts%no_of_allocatedvertices = grid_obj%nverts

    CALL construct_verts(grid_obj%verts)
    CALL construct_edges(grid_obj%edges)
    CALL construct_cells(grid_obj%cells)

    grid_obj%is_allocated  = .true.
    CALL id_grid_idxs(grid_id)

  END SUBROUTINE allocate_grid_object
  !-----------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Sets the filename associated with the grid_id
  SUBROUTINE set_grid_filename(grid_id, file_name)
    INTEGER, INTENT(in) :: grid_id
    CHARACTER(LEN=*), INTENT(in) :: file_name

    CALL check_active_grid_id(grid_id)

    get_grid_object(grid_id)%file_name = file_name

  END SUBROUTINE set_grid_filename
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Sets the grid_creation_process flag associated with the grid_id
  SUBROUTINE set_grid_creation(grid_id, creation_id)
    INTEGER, INTENT(in) :: grid_id, creation_id

    CALL check_active_grid_id(grid_id)

    get_grid_object(grid_id)%grid_creation_process = creation_id

  END SUBROUTINE set_grid_creation
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Sets the grid_creation_process flag associated with the grid_id
  SUBROUTINE set_grid_netcdf_flags(grid_id, netcdf_flags)
    INTEGER, INTENT(in) :: grid_id, netcdf_flags

    CALL check_active_grid_id(grid_id)

    get_grid_object(grid_id)%netcdf_flags = netcdf_flags

  END SUBROUTINE set_grid_netcdf_flags
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Sets the start_subgrid_id
  SUBROUTINE set_start_subgrids(grid_id, start_subgrids)
    INTEGER, INTENT(in) :: grid_id, start_subgrids

    CALL check_active_grid_id(grid_id)

    get_grid_object(grid_id)%start_subgrid_id = start_subgrids

  END SUBROUTINE set_start_subgrids
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  !! Sets the grid level
  SUBROUTINE set_grid_level(grid_id, new_level)
    INTEGER, INTENT(in) :: grid_id, new_level

    CALL check_active_grid_id(grid_id)

    get_grid_object(grid_id)%level = new_level

  END SUBROUTINE set_grid_level
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Gets the grid level
  INTEGER FUNCTION get_grid_level(grid_id)
    INTEGER, INTENT(in) :: grid_id

    CALL check_active_grid_id(grid_id)

    get_grid_level = get_grid_object(grid_id)%level

  END FUNCTION get_grid_level
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Sets the no_of_subgrids
  SUBROUTINE set_no_of_subgrids(grid_id, no_of_subgrids)
    INTEGER, INTENT(in) :: grid_id, no_of_subgrids

    CALL check_active_grid_id(grid_id)

    get_grid_object(grid_id)%no_of_subgrids = no_of_subgrids

  END SUBROUTINE set_no_of_subgrids
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Sets what pointers to use for the parents
  SUBROUTINE grid_set_sea_depth(grid_id, depth, cell_list)
    INTEGER, INTENT(in) :: grid_id
    REAL(wp), INTENT(in) :: depth
    TYPE(t_integer_list),INTENT(in), OPTIONAL :: cell_list

    TYPE(t_grid_cells), POINTER  :: cells
    INTEGER :: i
    
    CALL check_active_grid_id(grid_id)
    
    IF (PRESENT(cell_list) .AND. ASSOCIATED(cell_list%value)) THEN
      cells => get_grid_object(grid_id)%cells
      DO i=1,cell_list%list_size
        cells%elevation(cell_list%value(i)) = depth
      ENDDO
      RETURN
    ENDIF
      
    get_grid_object(grid_id)%cells%elevation(:) = depth
    get_grid_object(grid_id)%edges%elevation(:) = depth

  END SUBROUTINE grid_set_sea_depth
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Sets what pointers to use for the parents
  SUBROUTINE grid_set_min_sea_depth(grid_id, depth, cell_list)
    INTEGER, INTENT(in) :: grid_id
    REAL(wp), INTENT(in) :: depth
    TYPE(t_integer_list),INTENT(in), OPTIONAL :: cell_list

    REAL(wp), POINTER   :: elevation(:)
    INTEGER :: i

    CALL check_active_grid_id(grid_id)
    elevation =>  get_grid_object(grid_id)%cells%elevation

    IF (PRESENT(cell_list) .AND. ASSOCIATED(cell_list%value)) THEN
      DO i=1,cell_list%list_size
        IF (elevation(cell_list%value(i)) > depth) &
          elevation(cell_list%value(i)) = depth
      ENDDO
      RETURN
    ENDIF
    
    ! check all
    WHERE (elevation > depth )
      elevation = depth
    END WHERE

  END SUBROUTINE grid_set_min_sea_depth
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  !! Sets what pointers to use for the parents
  SUBROUTINE set_grid_parent_id(grid_id, parent_id)
    INTEGER, INTENT(in) :: grid_id, parent_id

    CALL check_active_grid_id(grid_id)

    get_grid_object(grid_id)%parent_grid_id = parent_id

  END SUBROUTINE set_grid_parent_id
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  !! Sets what pointers to use for the prents
  SUBROUTINE grid_set_parents_from(grid_id, parents_from)
    INTEGER, INTENT(in) :: grid_id, parents_from

    CALL check_active_grid_id(grid_id)
!     IF (PRESENT(redirect_id)) THEN
!       CALL finish('grid_set_parents_from','OPTIONAL redirect_id is not accepted')
!     ENDIF

    get_grid_object(grid_id)%parents_from = parents_from

  END SUBROUTINE grid_set_parents_from
  !-------------------------------------------------------------------------

 !-----------------------------------------------------------------------
  SUBROUTINE replace_grid(grid_id, new_grid_id)
    INTEGER, INTENT(in)      :: grid_id
    INTEGER, INTENT(inout)   :: new_grid_id

    TYPE(t_grid), POINTER :: switch_grid

    CALL check_active_grid_id(new_grid_id)
    switch_grid => get_grid_object(new_grid_id)
    get_grid_object(new_grid_id) => get_grid_object(grid_id)
    get_grid_object(grid_id) => switch_grid

    CALL delete_grid(new_grid_id)
    grid_is_active(grid_id) = .true.
    get_grid_object(grid_id)%grid_obj_id = grid_id

  END SUBROUTINE replace_grid
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !> copies from_grid_id to_grid_id
  !! An optional reindexing array maybe given.
  !! If the target grid is created, then the copy appends it.
  SUBROUTINE copy_grid(from_grid_id, to_grid_id, &
    & opt_new_vertindex, opt_new_edgeindex, opt_new_cellindex)
    INTEGER, INTENT(in)      :: from_grid_id
    INTEGER, INTENT(inout)   :: to_grid_id
    INTEGER, OPTIONAL, INTENT(in) :: opt_new_vertindex(0:)
    INTEGER, OPTIONAL, INTENT(in) :: opt_new_edgeindex(0:)
    INTEGER, OPTIONAL, INTENT(in) :: opt_new_cellindex(0:)

    TYPE(t_grid), POINTER          :: from_grid, to_grid
    TYPE(t_grid_vertices), POINTER :: from_verts, to_verts
    TYPE(t_grid_edges), POINTER    :: from_edges, to_edges
    TYPE(t_grid_cells), POINTER    :: from_cells, to_cells
    INTEGER, POINTER :: parent_cell_pnt(:), parent_edge_pnt(:), parent_vertex_pnt(:)

    LOGICAL :: reindex
    INTEGER :: max_vertex_connectivity, max_no_of_cell_vertices
    INTEGER :: start_no_of_verts,start_no_of_edges,start_no_of_cells
    INTEGER :: no_of_verts,no_of_edges,no_of_cells
    INTEGER :: no_of_copy_verts,no_of_copy_edges,no_of_copy_cells
    INTEGER :: no_of_out_verts,no_of_out_edges,no_of_out_cells
    INTEGER :: i,j,to_index
    INTEGER :: add_subgrid_id

    reindex = PRESENT(opt_new_vertindex)
    IF (reindex .AND. .NOT. &
      & (PRESENT(opt_new_edgeindex) .AND. PRESENT(opt_new_cellindex) )) &
        & CALL finish('mo_local_grid: copy grid entities','missing optional index')


    CALL check_active_grid_id(from_grid_id)
    from_grid  => get_grid_object(from_grid_id)
    from_verts => from_grid%verts
    from_edges => from_grid%edges
    from_cells => from_grid%cells
    no_of_verts = from_verts%no_of_existvertices
    no_of_edges = from_edges%no_of_existedges
    no_of_cells = from_cells%no_of_existcells
    max_vertex_connectivity =  from_verts%max_connectivity
    max_no_of_cell_vertices   = from_cells%max_no_of_vertices
    CALL grid_get_parent_pointers(from_grid_id, parent_cell_pnt, &
      & parent_edge_pnt, parent_vertex_pnt)

    IF (reindex) THEN
      no_of_copy_verts = MAXVAL(opt_new_vertindex)
      no_of_copy_edges = MAXVAL(opt_new_edgeindex)
      no_of_copy_cells = MAXVAL(opt_new_cellindex)
    ELSE
      no_of_copy_verts = no_of_verts
      no_of_copy_edges = no_of_edges
      no_of_copy_cells = no_of_cells
    ENDIF

    IF (to_grid_id < 1) THEN
      ! create the output grid
      to_grid_id = new_grid()
      to_grid    => get_grid_object(to_grid_id)
      start_no_of_verts = 0
      start_no_of_edges = 0
      start_no_of_cells = 0
      to_grid%nverts = no_of_copy_verts
      to_grid%nedges = no_of_copy_edges
      to_grid%ncells = no_of_copy_cells
      no_of_out_verts = no_of_copy_verts + start_no_of_verts
      no_of_out_edges = no_of_copy_edges + start_no_of_edges
      no_of_out_cells = no_of_copy_cells + start_no_of_cells
      to_grid%cells%max_no_of_vertices = max_no_of_cell_vertices
      to_grid%verts%max_connectivity   = max_vertex_connectivity
      to_grid%start_subgrid_id = from_grid%start_subgrid_id
      CALL allocate_grid_object(to_grid_id)
    ELSE
      ! check the to_grid_id
      to_grid   => get_grid(to_grid_id)
      start_no_of_verts = to_grid%verts%no_of_existvertices
      start_no_of_edges = to_grid%edges%no_of_existedges
      start_no_of_cells = to_grid%cells%no_of_existcells
      no_of_out_verts = no_of_copy_verts + start_no_of_verts
      no_of_out_edges = no_of_copy_edges + start_no_of_edges
      no_of_out_cells = no_of_copy_cells + start_no_of_cells
    ENDIF

    to_verts   => to_grid%verts
    to_edges   => to_grid%edges
    to_cells   => to_grid%cells
    ! check if we have enough allocated space
    IF (no_of_out_verts > to_verts%no_of_allocatedvertices) &
       & CALL finish('grid_copy_vertices','no_of_out_verts > no_of_allocatedvertices')
    IF (no_of_out_edges > to_edges%no_of_allocatededges) &
      & CALL finish('grid_copy_edges','no_of_out_edges > no_of_allocatededges')
    IF (no_of_out_cells > to_cells%no_of_allocatedcells) &
      & CALL finish('grid_copy_cells','no_of_out_cells > no_of_allocatedcells')
    to_verts%no_of_existvertices = no_of_out_verts
    to_edges%no_of_existedges = no_of_out_edges
    to_cells%no_of_existcells = no_of_out_cells

    to_grid%file_name                 = from_grid%file_name
    to_grid%out_file_name             = from_grid%out_file_name
    to_grid%netcdf_flags              = from_grid%netcdf_flags
    
    to_grid%grid_creation_process     = from_grid%grid_creation_process
    to_grid%grid_optimization_process = from_grid%grid_optimization_process
    to_grid%grid_geometry             = from_grid%grid_geometry
    to_grid%level                     = from_grid%level
    to_grid%vertical_structure        => from_grid%vertical_structure
    to_grid%patch_id                  = from_grid%patch_id
    to_grid%hierarchy_node_type       = from_grid%hierarchy_node_type
    ! to_grid%parents_from              = from_grid%parents_from
    to_grid%parent_grid_id            = from_grid%parent_grid_id
    to_grid%no_of_children            = from_grid%no_of_children
    to_grid%is_allocated              = .true.
    to_grid%is_filled                 = .true.


    add_subgrid_id                    = to_grid%start_subgrid_id + to_grid%no_of_subgrids &
      & - from_grid%start_subgrid_id
    to_grid%no_of_subgrids            = from_grid%no_of_subgrids + to_grid%no_of_subgrids

!     print *, 'from_grid%start_subgrid_id:',from_grid%start_subgrid_id
!     print *, 'to_grid%start_subgrid_id:',to_grid%start_subgrid_id
!     print *, 'from_grid%no_of_subgrids :', from_grid%no_of_subgrids
!     print *, 'to_grid%no_of_subgrids :', to_grid%no_of_subgrids
!     print *, 'add_subgrid_id :', add_subgrid_id


   ! print *, 'add_subgrid_id:',add_subgrid_id, ' to_grid%no_of_subgrids :',to_grid%no_of_subgrids

    IF (to_grid%no_of_children > 0) THEN
      ALLOCATE(to_grid%child_grid_ids(to_grid%no_of_children),stat=i)
      IF (i > 0) THEN
        WRITE (message_text, '(a,i8,a)') &
          & 'ALLOCATE grid_obj%child_grid_ids with ', &
          & to_grid%no_of_children, '  children.'
        CALL finish ('copy_grid', TRIM(message_text))
      ENDIF
      to_grid%child_grid_ids(1:to_grid%no_of_children) = &
        & from_grid%child_grid_ids(1:to_grid%no_of_children)
    ENDIF
    !-----------------------------------------------------------------------
    ! Note: this is unsafe
    to_verts%start_idx(:,:) = from_verts%start_idx(:,:)
    to_verts%end_idx(:,:)   = from_verts%end_idx(:,:)

    to_edges%start_idx(:,:) = from_edges%start_idx(:,:)
    to_edges%end_idx(:,:)   = from_edges%end_idx(:,:)

    to_cells%min_sea_depth  = from_cells%min_sea_depth
    to_cells%max_no_of_vertices = from_cells%max_no_of_vertices
    to_cells%start_idx(:,:) = from_cells%start_idx(:,:)
    to_cells%end_idx(:,:)   = from_cells%end_idx(:,:)
    !-----------------------------------------------------------------------

!$OMP PARALLEL
    !-----------------------------------------------------------------------
    ! copy vertices
!$OMP DO PRIVATE(i,to_index,j)
    DO i=1,no_of_verts

      IF (reindex) THEN
        to_index = opt_new_vertindex(i)
        IF (to_index == 0) CYCLE

        to_index = to_index + start_no_of_verts

        DO j=1,max_vertex_connectivity
          to_verts%get_neighbor_index(to_index,j) = &
            & opt_new_vertindex(from_verts%get_neighbor_index(i,j))
          to_verts%get_edge_index    (to_index,j) = &
            & opt_new_edgeindex(from_verts%get_edge_index    (i,j))
          to_verts%get_cell_index    (to_index,j) =  &
            & opt_new_cellindex(from_verts%get_cell_index    (i,j))
        ENDDO

      ELSE

        to_index = i + start_no_of_verts

        DO j=1,max_vertex_connectivity
          to_verts%get_neighbor_index(to_index,j) = &
            & from_verts%get_neighbor_index(i,j)
          to_verts%get_edge_index    (to_index,j) = &
            & from_verts%get_edge_index    (i,j)
          to_verts%get_cell_index    (to_index,j) =  &
            & from_verts%get_cell_index    (i,j)
        ENDDO

      ENDIF !(reindex)

      ! print *,  'Reindex verts:',i,'->',to_index
      to_verts%idx              (to_index)   = to_index
      to_verts%get_edge_orient  (to_index,:) = from_verts%get_edge_orient (i,:)
      to_verts%dual_area        (to_index)   = from_verts%dual_area        (i)
      to_verts%vertex           (to_index)   = from_verts%vertex           (i)
      to_verts%cartesian        (to_index)   = from_verts%cartesian        (i)
      to_verts%parent_index     (to_index)   = parent_vertex_pnt           (i)
      to_verts%parent_child_type(to_index)   = from_verts%parent_child_type(i)
      to_verts%refin_ctrl       (to_index)   = from_verts%refin_ctrl       (i)
      to_verts%child_id         (to_index)   = from_verts%child_id         (i)
      to_verts%no_of_neigbors   (to_index)   = from_verts%no_of_neigbors   (i)

    ENDDO ! i=1,no_of_verts
!$OMP END DO

    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    ! copy edges
!$OMP DO PRIVATE(i,to_index)
    DO i=1,no_of_edges
      IF (reindex) THEN
        to_index = opt_new_edgeindex(i)
        IF (to_index == 0) CYCLE

        to_index = to_index + start_no_of_edges

        to_edges%get_vertex_index(to_index,1) = &
          & opt_new_vertindex(from_edges%get_vertex_index(i,1))
        to_edges%get_vertex_index(to_index,2) = &
          & opt_new_vertindex(from_edges%get_vertex_index(i,2))
        to_edges%get_cell_index  (to_index,1) = &
          & opt_new_cellindex(from_edges%get_cell_index(i,1))
        to_edges%get_cell_index  (to_index,2) = &
          & opt_new_cellindex(from_edges%get_cell_index(i,2))

      ELSE
        to_index = i + start_no_of_edges

        to_edges%get_vertex_index(to_index,:) = &
          & from_edges%get_vertex_index(i,:)
        to_edges%get_cell_index  (to_index,:) = &
          & from_edges%get_cell_index(i,:)

      ENDIF

      ! print *, 'Reindex edges:',i,'->',to_index
      to_edges%idx               (to_index)   = to_index
      to_edges%system_orientation(to_index)   = from_edges%system_orientation(i)
      to_edges%center            (to_index)   = from_edges%center            (i)
      to_edges%cartesian_center  (to_index)   = from_edges%cartesian_center  (i)
      to_edges%primal_normal     (to_index)   = from_edges%primal_normal     (i)
      to_edges%cartesian_primal_normal(to_index)= from_edges%cartesian_primal_normal(i)
      to_edges%dual_normal       (to_index)   = from_edges%dual_normal       (i)
      to_edges%cartesian_dual_normal(to_index)= from_edges%cartesian_dual_normal(i)
      to_edges%primal_edge_length(to_index)   = from_edges%primal_edge_length(i)
      to_edges%dual_edge_length  (to_index)   = from_edges%dual_edge_length  (i)
      to_edges%quad_area         (to_index)   = from_edges%quad_area         (i)
      to_edges%get_edge_vert_length(to_index,:) = from_edges%get_edge_vert_length(i,:)
      to_edges%get_edge_cell_length(to_index,:) = from_edges%get_edge_cell_length(i,:)

      to_edges%elevation         (to_index)   = from_edges%elevation         (i)
      to_edges%sea_land_mask     (to_index)   = from_edges%sea_land_mask     (i)
      to_edges%no_of_zlayers     (to_index)   = from_edges%no_of_zlayers     (i)

      to_edges%refin_ctrl        (to_index)   = from_edges%refin_ctrl        (i)
      to_edges%subgrid_id        (to_index)   = from_edges%subgrid_id        (i)
      to_edges%parent_index      (to_index)   = parent_edge_pnt              (i)
      to_edges%parent_child_type (to_index)   = from_edges%parent_child_type (i)
      to_edges%child_index       (to_index,:) = from_edges%child_index       (i,:)
      to_edges%child_id          (to_index)   = from_edges%child_id          (i)
      to_edges%subgrid_id        (to_index)   = from_edges%subgrid_id        (i) + add_subgrid_id

    ENDDO ! i=1,noOfEdges
!$OMP END DO
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    ! copy cells
!$OMP DO PRIVATE(i,to_index,j)
    DO i=1,no_of_cells

      IF (reindex) THEN
        to_index = opt_new_cellindex(i)
        IF (to_index == 0) CYCLE

        to_index = to_index + start_no_of_cells

        DO j=1,max_no_of_cell_vertices

          to_cells%get_neighbor_index(to_index,j) = &
            & opt_new_cellindex(from_cells%get_neighbor_index(i,j))
          to_cells%get_vertex_index  (to_index,j) = &
            & opt_new_vertindex(from_cells%get_vertex_index  (i,j))
          to_cells%get_edge_index    (to_index,j) = &
            & opt_new_edgeindex(from_cells%get_edge_index    (i,j))

        ENDDO

      ELSE

        to_index = i + start_no_of_cells

        to_cells%get_neighbor_index(to_index,:) = &
          & from_cells%get_neighbor_index(i,:)
        to_cells%get_vertex_index  (to_index,:) = &
          & from_cells%get_vertex_index  (i,:)
        to_cells%get_edge_index    (to_index,:) = &
          & from_cells%get_edge_index    (i,:)

      ENDIF

      ! print *, 'Reindex cells:',i,'->',to_index
      to_cells%idx          (to_index) = to_index
      to_cells%center       (to_index) = from_cells%center       (i)
      to_cells%cartesian_center(to_index) = from_cells%cartesian_center(i)
      to_cells%area         (to_index) = from_cells%area         (i)
      to_cells%elevation    (to_index) = from_cells%elevation    (i)
      to_cells%sea_land_mask(to_index) = from_cells%sea_land_mask(i)
      to_cells%no_of_zlayers(to_index) = from_cells%no_of_zlayers(i)

      to_cells%no_of_vertices    (to_index)   = from_cells%no_of_vertices    (i)
      to_cells%get_edge_orient   (to_index,:) = from_cells%get_edge_orient   (i,:)

      to_cells%subgrid_id      (to_index)   = from_cells%subgrid_id   (i)
      to_cells%refin_ctrl      (to_index)   = from_cells%refin_ctrl   (i)
      to_cells%parent_index    (to_index)   = parent_cell_pnt             (i)
      to_cells%parent_child_type(to_index)  = from_cells%parent_child_type(i)
      to_cells%child_index     (to_index,:) = from_cells%child_index     (i,:)
      to_cells%child_id        (to_index)   = from_cells%child_id        (i)
      to_cells%subgrid_id      (to_index)   = from_cells%subgrid_id      (i) + add_subgrid_id

    ENDDO
!$OMP END DO

    !-----------------------------------------------------------------------
    !   add starting points to the connectivity parts
    IF (start_no_of_verts > 0) THEN
!$OMP DO PRIVATE(i)
      DO i=start_no_of_verts + 1, no_of_out_verts
        WHERE (to_verts%get_neighbor_index(i,:) /= 0)
               to_verts%get_neighbor_index(i,:) = &
            &  to_verts%get_neighbor_index(i,:) + start_no_of_verts
        END WHERE
        WHERE (to_verts%get_edge_index(i,:) /= 0)
               to_verts%get_edge_index(i,:) = &
            &  to_verts%get_edge_index(i,:) + start_no_of_edges
        END WHERE
        WHERE (to_verts%get_cell_index(i,:) /= 0)
               to_verts%get_cell_index(i,:) = &
            &  to_verts%get_cell_index(i,:) + start_no_of_cells
        END WHERE
      ENDDO
!$OMP END DO NOWAIT
!$OMP DO PRIVATE(i)
      DO i=start_no_of_edges + 1, no_of_out_edges
        WHERE (to_edges%get_vertex_index(i,:) /= 0)
               to_edges%get_vertex_index(i,:) = &
          &    to_edges%get_vertex_index(i,:) + start_no_of_verts
        ENDWHERE
        WHERE (to_edges%get_cell_index  (i,:) /= 0)
               to_edges%get_cell_index  (i,:) = &
          &    to_edges%get_cell_index  (i,:) + start_no_of_cells
        ENDWHERE
      ENDDO
!$OMP END DO NOWAIT
!$OMP DO PRIVATE(i)
      DO i=start_no_of_cells + 1, no_of_out_cells
        WHERE (to_cells%get_neighbor_index(i,:) /= 0)
               to_cells%get_neighbor_index(i,:) = &
          &    to_cells%get_neighbor_index(i,:) + start_no_of_cells
        ENDWHERE
        WHERE (to_cells%get_vertex_index  (i,:) /= 0)
               to_cells%get_vertex_index  (i,:) = &
          &    to_cells%get_vertex_index  (i,:) + start_no_of_verts
        ENDWHERE
        WHERE (to_cells%get_edge_index    (i,:) /= 0)
               to_cells%get_edge_index    (i,:) = &
          &    to_cells%get_edge_index    (i,:) + start_no_of_edges
        ENDWHERE
      ENDDO
!$OMP END DO NOWAIT

     ENDIF !(start_no_of_verts > 0)
!$OMP END PARALLEL


  END SUBROUTINE copy_grid
  !-----------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !   grid_get_parent_pointers(parent_grid_id, parent_cell, parent_edge, parent_vertex)
  !>
  !! Returns the actual pointers of parent entities for parent_grid_id
  SUBROUTINE grid_get_parent_pointers(grid_id, parent_cell, parent_edge, parent_vertex)
    INTEGER, INTENT(in) :: grid_id

    TYPE(t_grid), POINTER :: parent_grid
    INTEGER,    POINTER :: parent_cell(:), parent_edge(:), parent_vertex(:)
    INTEGER :: parents_from

    parent_grid => get_grid(grid_id)
    parents_from = parent_grid%parents_from
    ! parent indexes will climb-up when a CUT_OFF_GRID
    IF (parents_from == undefined) THEN
!       print *, 'parents_from == undefined'
!       print *, 'grid_creation_process:',parent_grid%grid_creation_process
!       print *, 'parent_grid_id:',parent_grid%parent_grid_id
      IF (parent_grid%grid_creation_process == cut_off_grid &
        & .OR. parent_grid%parent_grid_id /= undefined) THEN
        parents_from = parents_from_parentpointers
      ELSE
        parents_from = parents_from_idpointers
      ENDIF
    ENDIF

    IF (parents_from == parents_from_idpointers) THEN
      parent_cell   => parent_grid%cells%idx
      parent_edge   => parent_grid%edges%idx
      parent_vertex => parent_grid%verts%idx
!      print *, 'grid_get_parent_pointers form idx'
    ELSE
      parent_cell   => parent_grid%cells%parent_index
      parent_edge   => parent_grid%edges%parent_index
      parent_vertex => parent_grid%verts%parent_index
!      print *, 'grid_get_parent_pointers form parent_index'
    ENDIF

  END SUBROUTINE grid_get_parent_pointers
  !-------------------------------------------------------------------------


  !-----------------------------------------------------------------------
  !>
  !! Sets the actual number of entities equal to the allocated
  SUBROUTINE grid_set_exist_eq_allocated(grid_id)
    INTEGER, INTENT(in) :: grid_id

    TYPE(t_grid), POINTER :: grid_obj

    CALL check_active_grid_id(grid_id)
    grid_obj => get_grid_object(grid_id)
      grid_obj%verts%no_of_existvertices = &
        & grid_obj%verts%no_of_allocatedvertices
      grid_obj%edges%no_of_existedges    = &
        & grid_obj%edges%no_of_allocatededges
      grid_obj%cells%no_of_existcells    = &
        & grid_obj%cells%no_of_allocatedcells
  END SUBROUTINE grid_set_exist_eq_allocated
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  !! Sets the allocated number of entities equal to actual
  SUBROUTINE grid_set_allocated_eq_exist(grid_id)
    INTEGER, INTENT(in) :: grid_id

    TYPE(t_grid), POINTER :: grid_obj

    CALL check_active_grid_id(grid_id)

    grid_obj => get_grid_object(grid_id)
    grid_obj%ncells = grid_obj%cells%no_of_existcells
    grid_obj%nedges = grid_obj%edges%no_of_existedges
    grid_obj%nverts = grid_obj%verts%no_of_existvertices

    grid_obj%verts%no_of_allocatedvertices = grid_obj%nverts
    grid_obj%edges%no_of_allocatededges    = grid_obj%nedges
    grid_obj%cells%no_of_allocatedcells    = grid_obj%ncells

  END SUBROUTINE grid_set_allocated_eq_exist
  !-----------------------------------------------------------------------


  !-------------------------------------------------------------------------
  SUBROUTINE print_grid_edge(in_grid_id, edge_no)
    INTEGER, INTENT(in) :: in_grid_id, edge_no

    TYPE(t_grid), POINTER :: in_grid
    TYPE(t_grid_edges), POINTER :: edges
    TYPE(t_grid_vertices), POINTER :: verts

    INTEGER :: vertex_no

    in_grid => get_grid(in_grid_id)
    edges => in_grid%edges
    verts  => in_grid%verts
    print *, 'Grid:',in_grid_id, " edge no:", edge_no
    print *, "  cells:",edges%get_cell_index(edge_no,:)
    print *, "  center:",edges%center(edge_no)
    print *, "  normal:",edges%primal_normal(edge_no)
    vertex_no = edges%get_vertex_index(edge_no, 1)
    print *, "  vertex 1:",vertex_no, &
      & ' Lon,Lat:', verts%vertex(vertex_no)%lon, verts%vertex(vertex_no)%lat
    vertex_no = edges%get_vertex_index(edge_no, 2)
    print *, "  vertex 2:",vertex_no, &
      & ' Lon,Lat:', verts%vertex(vertex_no)%lon, verts%vertex(vertex_no)%lat
    print *, "  parent_id:",in_grid%parent_grid_id
    print *, "  parent_index:", edges%parent_index(edge_no)
    print *, "  child_id:", edges%child_id(edge_no)
    print *, "  child_index:", edges%child_index(edge_no,:)

    !-----------------------------------------------------------

    RETURN

  END SUBROUTINE print_grid_edge
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE print_grid_cell(in_grid_id, cell_no, opt_text)
    INTEGER, INTENT(in) :: in_grid_id, cell_no
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: opt_text

    TYPE(t_grid), POINTER :: in_grid
    TYPE(t_grid_cells), POINTER :: cell
    TYPE(t_grid_edges), POINTER :: edge
    TYPE(t_grid_vertices), POINTER :: vert

    INTEGER :: j,edge_no,vert_no

    in_grid => get_grid(in_grid_id)
    cell => in_grid%cells
    edge => in_grid%edges
    vert => in_grid%verts

    print*, '======================================================'
    print*, 'Grid ID:', in_grid_id, ' filename=', TRIM(in_grid%file_name), &
      & ' Cell No=', cell_no
    IF (PRESENT(opt_text)) print *, TRIM(opt_text)
    print*, '  subgrid_id :', cell%subgrid_id(cell_no)
    print*, '  refin_ctrl :', cell%refin_ctrl(cell_no)
    print*, '  child_id :', cell%child_id(cell_no)
    print*, '  Cell Edges:', cell%get_edge_index(cell_no,:)
    print*, '  Cell edge orientation:', cell%get_edge_orient(cell_no,:)
    print*, '  Edges:'
    DO j=1,cell%max_no_of_vertices
      edge_no = cell%get_edge_index(cell_no,j)
      print*, '   edge:',j, edge_no, edge%refin_ctrl(edge_no),edge%child_id(edge_no)
      print*, '      edge orientation:',edge%system_orientation(edge_no)
      print*, '      edge cells:',edge%get_cell_index(edge_no,:)
      print*, '      edge verts:',edge%get_vertex_index(edge_no,:)
    ENDDO
    print*, '  Vertices:'
    DO j=1,cell%max_no_of_vertices
      vert_no = cell%get_vertex_index(cell_no,j)
      print*, '   vertex:', vert_no, vert%refin_ctrl(vert_no),vert%child_id(vert_no)
      print*, '    vertex cells:', vert%get_cell_index(vert_no,:)
      print*, '    vertex edges:', vert%get_edge_index(vert_no,:)
      print*, '    vertex lon,lat:', vert%vertex(vert_no)
      print*, '    vertex cartesian:', vert%cartesian(vert_no)
    ENDDO
    print*, 'END print Cell ', cell_no
    print*, '======================================================'

  END SUBROUTINE print_grid_cell
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE print_grid_vertex(in_grid_id, vertex_no)
    INTEGER, INTENT(in) :: in_grid_id, vertex_no

    TYPE(t_grid_vertices), POINTER :: vertices

    vertices => get_vertices(in_grid_id)

    print*, 'Grid ID:', in_grid_id, ' Vertex No=', vertex_no
    print*, "   LonLat:", vertices%vertex(vertex_no)
    print*, "   Cartesian:", vertices%cartesian(vertex_no)
    print*, 'END print Vertex ', vertex_no

  END SUBROUTINE print_grid_vertex
  !-------------------------------------------------------------------------


!--------------------------------------------------------------
  !>
  !! Allocates the cells arrays.
  !! The array sizes are defined by:
  !! cells\%no_of_allocatedcells, cells\%max_no_of_vertices
  SUBROUTINE construct_cells(cells)
    TYPE(t_grid_cells), INTENT(inout) :: cells

    INTEGER :: no_of_cells, max_cell_vertices
    INTEGER :: istat, ist

    no_of_cells       = cells%no_of_allocatedcells
    max_cell_vertices = cells%max_no_of_vertices
    cells%no_of_existcells = 0

    ist = 0

    ALLOCATE(cells%idx(no_of_cells),stat=istat)
    ist=ist+istat
    cells%idx(:)=0

    ALLOCATE(cells%get_neighbor_index(no_of_cells,max_cell_vertices) ,stat=istat)
    ist=ist+istat
    cells%get_neighbor_index(:,:)=0

    ALLOCATE(cells%get_edge_index (no_of_cells,max_cell_vertices),stat=istat)
    ist=ist+istat
    cells%get_edge_index (:,:)=0

    ALLOCATE(cells%get_vertex_index(no_of_cells,max_cell_vertices),stat=istat)
    ist=ist+istat
    cells%get_vertex_index(:,:)=0

    ALLOCATE(cells%get_edge_orient(no_of_cells,max_cell_vertices),stat=istat)
    ist=ist+istat
    cells%get_edge_orient(:,:)=0

    ALLOCATE(cells%center(no_of_cells),stat=istat)
    ist=ist+istat
    cells%center(:)%lon=0.0_wp
    cells%center(:)%lat=0.0_wp

    ALLOCATE(cells%cartesian_center(no_of_cells),stat=istat)
    ist=ist+istat
    cells%cartesian_center(:)%x(1) =0.0_wp
    cells%cartesian_center(:)%x(2) =0.0_wp
    cells%cartesian_center(:)%x(3) =0.0_wp

    ALLOCATE(cells%area(no_of_cells),stat=istat)
    ist=ist+istat
    cells%area(:)=0.0_wp

    ALLOCATE(cells%subgrid_id(no_of_cells),stat=istat)
    ist=ist+istat
    cells%subgrid_id(:)=undefined

    ALLOCATE(cells%refin_ctrl(no_of_cells),stat=istat)
    ist=ist+istat
    cells%refin_ctrl(:)=0

    !ALLOCATE(cells%indlist(min_rlcell:max_rlcell,2),stat=istat)
    !ist=ist+istat
    !cells%indlist(:,:)=0

    ALLOCATE(cells%start_idx(min_rlcell:max_rlcell,1),stat=istat)
    ist=ist+istat
    cells%start_idx(:,:)=0

    ALLOCATE(cells%end_idx(min_rlcell:max_rlcell,1),stat=istat)
    ist=ist+istat
    cells%end_idx(:,:)=0

    ALLOCATE(cells%no_of_vertices(no_of_cells),stat=istat)
    ist=ist+istat
    cells%no_of_vertices(:)=0

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

    ALLOCATE(cells%elevation(no_of_cells),stat=istat)
    ist=ist+istat
    cells%elevation(1:) = 0.0_wp

    ALLOCATE(cells%sea_land_mask(no_of_cells),stat=istat)
    ist=ist+istat
    cells%sea_land_mask(1:) = undefined

    ALLOCATE(cells%no_of_zlayers(no_of_cells),stat=istat)
    ist=ist+istat
    cells%no_of_zlayers(1:) = undefined

    IF (ist>0) THEN
      WRITE (message_text, '(a,i8,a)') &
        & 'construct_cells with ', no_of_cells, ' cells.'
      CALL finish ('construct_cells', TRIM(message_text))
    ENDIF

  END SUBROUTINE construct_cells
  !--------------------------------------------------------------

  !--------------------------------------------------------------
  !>
  !! Allocates the edges arrays.
  !! The array sizes are defined by:
  !! edges\%no_of_allocatededge
  SUBROUTINE construct_edges(edges)
    TYPE(t_grid_edges), INTENT(inout) :: edges

    INTEGER :: no_of_edges
    INTEGER :: istat, ist

    no_of_edges = edges%no_of_allocatededges
    edges%no_of_existedges = 0

    ist = 0

    ALLOCATE(edges%idx(no_of_edges),stat=istat)
    ist=ist+istat
    edges%idx(:)=0

    ALLOCATE(edges%get_cell_index (no_of_edges,2),stat=istat)
    ist=ist+istat
    edges%get_cell_index(:,:)=0

    ALLOCATE(edges%get_vertex_index(no_of_edges,2),stat=istat)
    ist=ist+istat
    edges%get_vertex_index(:,:)=0

    ALLOCATE(edges%system_orientation(no_of_edges),stat=istat)
    ist=ist+istat
    edges%system_orientation(:)=0

    ALLOCATE(edges%center(no_of_edges),stat=istat)
    ist=ist+istat
    edges%center(:)%lat=0.0_wp
    edges%center(:)%lon=0.0_wp

    ALLOCATE(edges%cartesian_center(no_of_edges),stat=istat)
    ist=ist+istat
    edges%cartesian_center(:)%x(1) =0.0_wp
    edges%cartesian_center(:)%x(2) =0.0_wp
    edges%cartesian_center(:)%x(3) =0.0_wp

    ALLOCATE(edges%primal_normal(no_of_edges),stat=istat)
    ist=ist+istat
    edges%primal_normal(:)%v1=0.0_wp
    edges%primal_normal(:)%v2=0.0_wp

    ALLOCATE(edges%cartesian_primal_normal(no_of_edges),stat=istat)
    ist=ist+istat
    edges%cartesian_primal_normal(:)%x(1)=0.0_wp
    edges%cartesian_primal_normal(:)%x(2)=0.0_wp
    edges%cartesian_primal_normal(:)%x(3)=0.0_wp

    ALLOCATE(edges%dual_normal(no_of_edges),stat=istat)
    ist=ist+istat
    edges%dual_normal(:)%v1=0.0_wp
    edges%dual_normal(:)%v2=0.0_wp

    ALLOCATE(edges%cartesian_dual_normal(no_of_edges),stat=istat)
    ist=ist+istat
    edges%cartesian_dual_normal(:)%x(1) = 0.0_wp
    edges%cartesian_dual_normal(:)%x(2) = 0.0_wp
    edges%cartesian_dual_normal(:)%x(3) = 0.0_wp

    ALLOCATE(edges%primal_edge_length(no_of_edges),stat=istat)
    ist=ist+istat
    edges%primal_edge_length(:)=0.0_wp

    ALLOCATE(edges%dual_edge_length(no_of_edges),stat=istat)
    ist=ist+istat
    edges%dual_edge_length(:)=0.0_wp

    ALLOCATE(edges%quad_area(no_of_edges),stat=istat)
    ist=ist+istat
    edges%quad_area(:)=0.0_wp

    ALLOCATE(edges%get_edge_vert_length(no_of_edges,2),stat=istat)
    ist=ist+istat
    edges%get_edge_vert_length(:,:)=0.0_wp

    ALLOCATE(edges%get_edge_cell_length(no_of_edges,2),stat=istat)
    ist=ist+istat
    edges%get_edge_cell_length(:,:)=0.0_wp

    ALLOCATE(edges%elevation(no_of_edges),stat=istat)
    ist=ist+istat
    edges%elevation(:)=0.0_wp

    ALLOCATE(edges%sea_land_mask(no_of_edges),stat=istat)
    ist=ist+istat
    edges%sea_land_mask(:)=undefined

    ALLOCATE(edges%no_of_zlayers(no_of_edges),stat=istat)
    ist=ist+istat
    edges%no_of_zlayers(:)=-1

    ALLOCATE(edges%subgrid_id(no_of_edges),stat=istat)
    ist=ist+istat
    edges%subgrid_id(:)=undefined

    ALLOCATE(edges%refin_ctrl(no_of_edges),stat=istat)
    ist=ist+istat
    edges%refin_ctrl(:)=0

    !ALLOCATE(edges%indlist(min_rledge:max_rledge,2),stat=istat)
    !ist=ist+istat
    !edges%indlist(:,:)=0

    ALLOCATE(edges%start_idx(min_rledge:max_rledge,1),stat=istat)
    ist=ist+istat
    edges%start_idx(:,:)=0

    ALLOCATE(edges%end_idx(min_rledge:max_rledge,1),stat=istat)
    ist=ist+istat
    edges%end_idx(:,:)=0

    ALLOCATE(edges%parent_index(no_of_edges),stat=istat)
    ist=ist+istat
    edges%parent_index(:)=0

    ALLOCATE(edges%parent_child_type(no_of_edges),stat=istat)
    ist=ist+istat
    edges%parent_child_type(:)=0

    ALLOCATE(edges%child_index(no_of_edges,4),stat=istat)
    ist=ist+istat
    edges%child_index(:,:)=0

    ALLOCATE(edges%child_id(no_of_edges),stat=istat)
    ist=ist+istat
    edges%child_id(:)=0


    IF (ist>0) THEN
      WRITE (message_text, '(a,i8,a)') &
        & 'construct_edges with ', no_of_edges, ' edges.'
      CALL finish ('construct_edges', TRIM(message_text))
    ENDIF

  END SUBROUTINE construct_edges
  !--------------------------------------------------------------

  !--------------------------------------------------------------
  !>
  !! Allocates the vertices arrays.
  !! The array sizes are defined by:
  !! verts\%no_of_allocatedvertices,verts\%max_connectivity
  SUBROUTINE construct_verts(verts)
    TYPE(t_grid_vertices), INTENT(inout) :: verts

    INTEGER :: no_of_verts,max_vertex_connect
    INTEGER :: istat, ist

    no_of_verts        = verts%no_of_allocatedvertices
    max_vertex_connect = verts%max_connectivity
    verts%no_of_existvertices = 0

    ist = 0

    ALLOCATE(verts%idx(no_of_verts),stat=istat)
    ist=ist+istat
    verts%idx(:)=0

    ALLOCATE(verts%parent_index(no_of_verts),stat=istat)
    ist=ist+istat
    verts%parent_index(:)=0
    ALLOCATE(verts%parent_child_type(no_of_verts),stat=istat)
    ist=ist+istat
    verts%parent_child_type(:)=0

    ALLOCATE(verts%no_of_neigbors(no_of_verts),stat=istat)
    ist=ist+istat
    verts%no_of_neigbors(:)=0

    ALLOCATE(verts%get_neighbor_index(no_of_verts,max_vertex_connect) ,stat=istat)
    ist=ist+istat
    verts%get_neighbor_index(:,:)=0

    ALLOCATE(verts%get_cell_index (no_of_verts,max_vertex_connect),stat=istat)
    ist=ist+istat
    verts%get_cell_index (:,:)=0

    ALLOCATE(verts%get_edge_index (no_of_verts,max_vertex_connect),stat=istat)
    ist=ist+istat
    verts%get_edge_index (:,:)=0

    ALLOCATE(verts%get_edge_orient(no_of_verts,max_vertex_connect),stat=istat)
    ist=ist+istat
    verts%get_edge_orient(:,:)=0

    ALLOCATE(verts%dual_area(no_of_verts),stat=istat)
    ist=ist+istat
    verts%dual_area(:)=0.0_wp

    ALLOCATE(verts%vertex(no_of_verts),stat=istat)
    ist=ist+istat
    verts%vertex(:)%lon=0.0_wp
    verts%vertex(:)%lat=0.0_wp

    ALLOCATE(verts%cartesian(no_of_verts),stat=istat)
    ist=ist+istat
    verts%cartesian(:)%x(1)=0.0_wp
    verts%cartesian(:)%x(2)=0.0_wp
    verts%cartesian(:)%x(3)=0.0_wp

    ALLOCATE(verts%refin_ctrl(no_of_verts),stat=istat)
    ist=ist+istat
    verts%refin_ctrl(:)=0

!     ALLOCATE(verts%indlist(min_rlvert:max_rlvert,2),stat=istat)
!     ist=ist+istat
!     verts%indlist(:,:)=0

    ALLOCATE(verts%start_idx(min_rlvert:max_rlvert,1),stat=istat)
    ist=ist+istat
    verts%start_idx(:,:)=0

    ALLOCATE(verts%end_idx(min_rlvert:max_rlvert,1),stat=istat)
    ist=ist+istat
    verts%end_idx(:,:)=0

    ALLOCATE(verts%child_id(no_of_verts),stat=istat)
    ist=ist+istat
    verts%child_id(:)=0

    IF (ist>0) THEN
      WRITE (message_text, '(a,i8,a)') &
        & 'construct_verts with ', no_of_verts, ' vertices.'
      CALL finish ('construct_verts', TRIM(message_text))
    ENDIF

  END SUBROUTINE construct_verts


  !-----------------------------------------------------------------------
  !>
  !! Deletes all grid objects
  !! Note: Should be called only at the stop of a program.
  SUBROUTINE detsruct_grid_objects ()

    INTEGER :: i

    IF (no_of_allocated_grids == 0) THEN
      CALL finish('detsruct_grid_objects', 'grid_objects have not been constructed');
    ENDIF

    DO i=1,max_active_grids
      IF (grid_is_active(i)) THEN
        CALL delete_grid(i)
      ENDIF
    ENDDO

    DEALLOCATE(grid_object_list,grid_is_active,grid_pnt_array)

    no_of_allocated_grids  = 0
    active_grids     = 0
    max_active_grids = 0

  END SUBROUTINE detsruct_grid_objects
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  !! Deallocates the cells arrays.
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
    DEALLOCATE(cells%no_of_vertices,stat=istat)
    ist=ist+istat
    DEALLOCATE(cells%a_neighbor_index ,stat=istat)
    ist=ist+istat
    DEALLOCATE(cells%a_edge_index,stat=istat)
    ist=ist+istat
    DEALLOCATE(cells%a_vertex_index,stat=istat)
    ist=ist+istat
    DEALLOCATE(cells%a_edge_orientation,stat=istat)
    ist=ist+istat
    DEALLOCATE(cells%center,stat=istat)
    ist=ist+istat
    DEALLOCATE(cells%cartesian_center,stat=istat)
    ist=ist+istat
    DEALLOCATE(cells%area,stat=istat)
    ist=ist+istat
    DEALLOCATE(cells%subgrid_id,stat=istat)
    ist=ist+istat
    DEALLOCATE(cells%refin_ctrl,stat=istat)
    ist=ist+istat
!     DEALLOCATE(cells%indlist,stat=istat)
!     ist=ist+istat
    DEALLOCATE(cells%start_idx,stat=istat)
    ist=ist+istat
    DEALLOCATE(cells%end_idx,stat=istat)
    ist=ist+istat
    DEALLOCATE(cells%elevation,stat=istat)
    ist=ist+istat
    DEALLOCATE(cells%sea_land_mask,stat=istat)
    ist=ist+istat
    DEALLOCATE(cells%no_of_zlayers,stat=istat)
    ist=ist+istat

    IF (ist>0) THEN
      WRITE (message_text, '(a)') 'Deallocate grid.'
      CALL message ('destruct_cells', 'Deallocate grid.')
    ENDIF

  END SUBROUTINE destruct_cells
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  !! Deallocates the edges arrays.
  SUBROUTINE destruct_edges(edges)
    TYPE(t_grid_edges), INTENT(inout) :: edges

    INTEGER :: istat, ist

    ist=0

    DEALLOCATE(edges%idx,stat=istat)
    ist=ist+istat
    DEALLOCATE(edges%a_cell_index ,stat=istat)
    ist=ist+istat
    DEALLOCATE(edges%a_vertex_index,stat=istat)
    ist=ist+istat
    DEALLOCATE(edges%system_orientation,stat=istat)
    ist=ist+istat
    DEALLOCATE(edges%center,stat=istat)
    ist=ist+istat
    DEALLOCATE(edges%cartesian_center,stat=istat)
    ist=ist+istat
    DEALLOCATE(edges%primal_normal,stat=istat)
    ist=ist+istat
    DEALLOCATE(edges%cartesian_primal_normal,stat=istat)
    ist=ist+istat
    DEALLOCATE(edges%dual_normal,stat=istat)
    ist=ist+istat
    DEALLOCATE(edges%cartesian_dual_normal,stat=istat)
    ist=ist+istat
    DEALLOCATE(edges%primal_edge_length,stat=istat)
    ist=ist+istat
    DEALLOCATE(edges%dual_edge_length,stat=istat)
    ist=ist+istat
    DEALLOCATE(edges%quad_area,stat=istat)
    ist=ist+istat
    DEALLOCATE(edges%a_edge_vert_length,stat=istat)
    ist=ist+istat
    DEALLOCATE(edges%a_edge_cell_length,stat=istat)
    ist=ist+istat
    DEALLOCATE(edges%elevation,stat=istat)
    ist=ist+istat
    DEALLOCATE(edges%sea_land_mask,stat=istat)
    ist=ist+istat
    DEALLOCATE(edges%no_of_zlayers,stat=istat)
    ist=ist+istat
    DEALLOCATE(edges%subgrid_id,stat=istat)
    ist=ist+istat
    DEALLOCATE(edges%refin_ctrl,stat=istat)
    ist=ist+istat
    ! DEALLOCATE(edges%indlist,stat=istat)
    ! ist=ist+istat
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

    IF (ist>0) THEN
      WRITE (message_text, '(a)') 'Deallocate grid.'
      CALL message ('destruct_edges', 'Deallocate grid.')
    ENDIF

  END SUBROUTINE destruct_edges
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  !! Deallocates the vertices arrays.
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
    DEALLOCATE(verts%a_neighbor_index,stat=istat)
    ist=ist+istat
    DEALLOCATE(verts%a_cell_index,stat=istat)
    ist=ist+istat
    DEALLOCATE(verts%a_edge_index,stat=istat)
    ist=ist+istat
    DEALLOCATE(verts%a_edge_orientation,stat=istat)
    ist=ist+istat
    DEALLOCATE(verts%dual_area,stat=istat)
    ist=ist+istat
    DEALLOCATE(verts%vertex,stat=istat)
    ist=ist+istat
    DEALLOCATE(verts%cartesian,stat=istat)
    ist=ist+istat
    DEALLOCATE(verts%refin_ctrl,stat=istat)
    ist=ist+istat
    DEALLOCATE(verts%child_id,stat=istat)
    ist=ist+istat
    DEALLOCATE(verts%start_idx,stat=istat)
    ist=ist+istat
    DEALLOCATE(verts%end_idx,stat=istat)
    ist=ist+istat

    IF (ist>0) THEN
      WRITE (message_text, '(a)') 'Deallocate grid.'
      CALL message ('destruct_verts', 'Deallocate grid.')
    ENDIF

  END SUBROUTINE destruct_verts
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !   SUBROUTINE set_nest_defaultindexes(of_grid_id)
  !>
  !! Sets the grid nest interpolation indexes to default
  SUBROUTINE set_nest_defaultindexes(of_grid_id)
    INTEGER, INTENT(in):: of_grid_id

    TYPE(t_grid), POINTER :: of_grid

    INTEGER :: no_of_input_cells, no_of_input_edges, no_of_input_verts

    !--------------------------------------------------------------
    ! set interpolation indexes to default
    !--------------------------------------------------------------
    of_grid => get_grid_object(of_grid_id)
    no_of_input_cells = of_grid%cells%no_of_existcells
    no_of_input_edges = of_grid%edges%no_of_existedges
    no_of_input_verts = of_grid%verts%no_of_existvertices

    of_grid%verts%start_idx(1:max_rlvert ,1) = 1
    of_grid%verts%start_idx(min_rlvert:0, 1) = no_of_input_verts + 1
    of_grid%verts%end_idx(1:max_rlvert ,  1)   = 0
    of_grid%verts%end_idx(min_rlvert:0 ,  1)   = no_of_input_verts

    of_grid%edges%start_idx(1:max_rledge ,1) = 1
    of_grid%edges%start_idx(min_rledge:0, 1) = no_of_input_edges + 1
    of_grid%edges%end_idx(1:max_rledge ,  1) = 0
    of_grid%edges%end_idx(min_rledge:0 ,  1) = no_of_input_edges
    of_grid%edges%refin_ctrl        (:)      = min_rledge_int


    of_grid%cells%start_idx(1:max_rlcell ,1) = 1
    of_grid%cells%start_idx(min_rlcell:0, 1) = no_of_input_cells + 1
    of_grid%cells%end_idx(1:max_rlcell ,  1) = 0
    of_grid%cells%end_idx(min_rlcell:0 ,  1) = no_of_input_cells
     of_grid%cells%refin_ctrl          (:  ) = min_rlcell_int

    

  END SUBROUTINE set_nest_defaultindexes
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !    INTEGER FUNCTION get_next_cell_edge(in_grid_id, cell_no, edge_no, vertex_no)
  !>
  !! Returns the next edge of edge_no in the cell, defined by vertex_no
  !! This means the get_next_cell_edge is the edge in the cell that shares vertex_no
  !! with the edge_no
  SUBROUTINE get_next_cell_edge_vertex(in_grid_id, cell_no, edge_no, vertex_no,&
    & next_edge, next_vertex, edge_in_cell)
    INTEGER, INTENT(in)  :: in_grid_id, cell_no, edge_no, vertex_no
    INTEGER, INTENT(out)  :: next_edge, next_vertex,edge_in_cell

    TYPE(t_grid),       POINTER :: in_grid
    TYPE(t_grid_cells), POINTER :: cells
    TYPE(t_grid_edges), POINTER :: edges
    INTEGER :: max_no_of_vertices

    in_grid => get_grid_object(in_grid_id)
    cells=>in_grid%cells
    edges=>in_grid%edges
    max_no_of_vertices = cells%max_no_of_vertices
    ! WRITE(0,*) "cell_no, edge_no, vertex_no:", cell_no, edge_no, vertex_no

    DO edge_in_cell=1,max_no_of_vertices
      next_edge = cells%get_edge_index(cell_no,edge_in_cell)

      IF (next_edge /= edge_no) THEN
        IF (edges%get_vertex_index(next_edge,1) == vertex_no) THEN
          next_vertex =  edges%get_vertex_index(next_edge,2)
          RETURN
        ENDIF
        IF (edges%get_vertex_index(next_edge,2) == vertex_no) THEN
          next_vertex =  edges%get_vertex_index(next_edge,1)
          RETURN
        ENDIF
      ENDIF
    ENDDO

    WRITE(0,*) "cell_no, edge_no, vertex_no:", cell_no, edge_no, vertex_no
    WRITE(0,*) "cell_edges:", cells%get_edge_index(cell_no,:)
    WRITE(0,*) "cell_verts:", cells%get_vertex_index(cell_no,:)
    CALL finish ('get_next_cell_edge_vertex', 'Cannot find next edge in cell')

  END SUBROUTINE get_next_cell_edge_vertex
  !-------------------------------------------------------------------------

  !-----------------------------------------------------------------------
  ! Vertical Structure Methods
  !-----------------------------------------------------------------------

  !--------------------------------------------------------------
  !>
  !! Allocates the arrays of the vertical structure
  !! The following field define the sizes:
  !! vertical_srtc\%no_of_columns
  !! vertical_srtc\%no_of_zlayers
  SUBROUTINE allocate_vertical(vertical_srtc)
    TYPE(t_vertical_ocean_structure), INTENT(inout) :: vertical_srtc

    CALL allocate_vertical_zlayers(vertical_srtc)
    CALL allocate_vertical_columns(vertical_srtc)

  END SUBROUTINE allocate_vertical
  !--------------------------------------------------------------

  !--------------------------------------------------------------
  !>
  !! Deallocates the arrays of the vertical structure
  SUBROUTINE deallocate_vertical(vertical_srtc)
    TYPE(t_vertical_ocean_structure), INTENT(inout) :: vertical_srtc

    CALL deallocate_vertical_zlayers(vertical_srtc)
    CALL deallocate_vertical_columns(vertical_srtc)

  END SUBROUTINE deallocate_vertical
  !--------------------------------------------------------------

  !--------------------------------------------------------------
  !>
  !! Allocates the column arrays of the vertical structure
  !! The following field define the sizes:
  !! vertical_srtc\%no_of_columns
  !! vertical_srtc\%no_of_zlayers
  SUBROUTINE allocate_vertical_columns(vertical_srtc, opt_no_of_columns)
    TYPE(t_vertical_ocean_structure), INTENT(inout) :: vertical_srtc
    INTEGER, OPTIONAL, INTENT(in) :: opt_no_of_columns

    INTEGER :: i

    IF (PRESENT(opt_no_of_columns)) &
      &  vertical_srtc%no_of_columns = opt_no_of_columns

    ALLOCATE(vertical_srtc%column_cell_index &
      & (vertical_srtc%no_of_zlayers,vertical_srtc%no_of_columns),&
      &  vertical_srtc%column_size(vertical_srtc%no_of_columns),stat=i)
    IF (i>0) THEN
      CALL finish ('allocate_vertical_columns', 'Failed')
    ENDIF
  END SUBROUTINE allocate_vertical_columns
  !--------------------------------------------------------------

  !--------------------------------------------------------------
  !>
  !! Deallocates the column arrays of the vertical structure
  SUBROUTINE deallocate_vertical_columns(vertical_srtc)
    TYPE(t_vertical_ocean_structure), INTENT(inout) :: vertical_srtc

    DEALLOCATE(vertical_srtc%column_cell_index,&
      &  vertical_srtc%column_size)
  END SUBROUTINE deallocate_vertical_columns
  !--------------------------------------------------------------

  !--------------------------------------------------------------
  !>
  !! Allocates the layer arrays of the vertical structure
  !! The following field define the sizes:
  !! vertical_srtc\%no_of_zlayers
  SUBROUTINE allocate_vertical_zlayers(vertical_srtc, opt_no_of_zlayers)
    TYPE(t_vertical_ocean_structure), INTENT(inout) :: vertical_srtc
    INTEGER, OPTIONAL, INTENT(in) :: opt_no_of_zlayers

    INTEGER :: i,nlayers

    IF (PRESENT(opt_no_of_zlayers)) &
      &  vertical_srtc%no_of_zlayers = opt_no_of_zlayers

    nlayers = vertical_srtc%no_of_zlayers
    ALLOCATE(vertical_srtc%layer_thicknes(nlayers), &
      & vertical_srtc%layer_bed(nlayers),           &
      & vertical_srtc%layer_middle(nlayers),stat=i)
    IF (i>0) THEN
      CALL finish ('allocate_vertical_columns', 'Failed')
    ENDIF
  END SUBROUTINE allocate_vertical_zlayers
  !--------------------------------------------------------------

  !--------------------------------------------------------------
  !>
  !! Deallocates the layer arrays of the vertical structure
  SUBROUTINE deallocate_vertical_zlayers(vertical_srtc)
    TYPE(t_vertical_ocean_structure), INTENT(inout) :: vertical_srtc

    DEALLOCATE(vertical_srtc%layer_thicknes, &
      & vertical_srtc%layer_bed,           &
      & vertical_srtc%layer_middle)
  END SUBROUTINE deallocate_vertical_zlayers
  !--------------------------------------------------------------

  !--------------------------------------------------------------
  !>
  !! Sets the vertical structure pointer of the grid to the given pointer
  SUBROUTINE set_grid_vertical_structure(grid_id, vertical_srtc)
    INTEGER, INTENT(in) :: grid_id
    TYPE(t_vertical_ocean_structure), TARGET, INTENT(in) :: vertical_srtc

    CALL check_active_grid_id(grid_id)
    get_grid_object(grid_id)%vertical_structure => vertical_srtc

  END SUBROUTINE set_grid_vertical_structure
  !--------------------------------------------------------------

  !-------------------------------------------------------------------------
  !-------------------------------------------------------------------------
  ! Private Methods
  !-------------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  !! "Class construction" method.
  !! Called once to initiallize the "class".
  SUBROUTINE construct_grid_objects ()

    INTEGER :: total_status,STATUS,i

    IF (no_of_allocated_grids /= 0) THEN
      ! do some consistency checks, return if everything ok
      IF (no_of_allocated_grids /= max_no_of_grid_objects) THEN
        CALL finish('construct_grid_objects', &
          & 'no_of_allocated_grids/= max_no_of_grid_objects');
      ENDIF
      IF (.not. ALLOCATED(grid_object_list)) THEN
        CALL finish('construct_grid_objects', '.NOT. ALLOCATED(grid_object_list)');
      ENDIF
      IF (SIZE(grid_object_list) /= no_of_allocated_grids) THEN
        CALL finish('construct_grid_objects', &
          & 'SIZE(grid_object_list) /= no_of_allocated_grids');
      ENDIF

      RETURN
    ENDIF

    IF (ALLOCATED(grid_object_list)) THEN
      CALL finish('construct_grid_objects', 'ALLOCATED(grid_object_list)');
    ENDIF

    no_of_allocated_grids  = max_no_of_grid_objects
    active_grids     = 0
    max_active_grids = 0

    total_status = 0
    ALLOCATE(grid_object_list(no_of_allocated_grids),stat=STATUS)
    total_status=total_status + STATUS
    ALLOCATE(grid_pnt_array(no_of_allocated_grids),stat=STATUS)
    total_status=total_status + STATUS
    ALLOCATE(grid_is_active(no_of_allocated_grids),stat=STATUS)
    total_status=total_status + STATUS
    IF (total_status /= 0) THEN
      CALL finish('construct_grid_objects', 'failed to ALLOCATE(grid_object_list)');
    ENDIF
    DO i=1,no_of_allocated_grids
      grid_pnt_array(i)%pnt => grid_object_list(i)
    ENDDO
    grid_is_active(:) = .false.

  END SUBROUTINE construct_grid_objects
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  !! Initializes a new grid object
  SUBROUTINE init_grid_object(grid_id)
    INTEGER, INTENT(in) :: grid_id

    TYPE(t_grid), POINTER :: grid_obj

    grid_obj => get_grid_object(grid_id)
    grid_obj%grid_obj_id              = grid_id
    grid_obj%file_name                = ''
    grid_obj%out_file_name    = ''
    grid_obj%grid_creation_process    = undefined
    grid_obj%grid_optimization_process= undefined
    grid_obj%grid_geometry            = undefined
    grid_obj%netcdf_flags             = undefined
    grid_obj%netcdf_flags             = netcdf_CF_1_1_convention

    grid_obj%ncells   = 0
    grid_obj%nedges   = 0
    grid_obj%nverts   = 0
    grid_obj%verts%no_of_allocatedvertices = 0
    grid_obj%edges%no_of_allocatededges = 0
    grid_obj%cells%no_of_allocatedcells = 0
    grid_obj%verts%no_of_existvertices = 0
    grid_obj%edges%no_of_existedges = 0
    grid_obj%cells%no_of_existcells = 0
    grid_obj%level                    = 0
    NULLIFY(grid_obj%vertical_structure)
    grid_obj%cells%min_sea_depth      = 0.0_wp

    ! zero the grid nesting info
!     grid_obj%patch_id            = undefined
    grid_obj%patch_id            = -1
    grid_obj%hierarchy_node_type = undefined
    grid_obj%parent_grid_id      = undefined
    grid_obj%parents_from        = undefined
    grid_obj%no_of_children      = 0
    grid_obj%no_of_subgrids      = 0
    grid_obj%start_subgrid_id    = 0
    NULLIFY(grid_obj%child_grid_ids)

    grid_obj%is_allocated  = .false.
    grid_obj%is_filled  = .false.
    grid_is_active(grid_id) = .true.

  END SUBROUTINE init_grid_object
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  !! Deallocates a grid object
  SUBROUTINE deallocate_grid_object(grid_id)
    INTEGER, INTENT(in) :: grid_id

    TYPE(t_grid), POINTER :: grid_obj

    CALL check_active_grid_id(grid_id)
    grid_obj => get_grid_object(grid_id)

    IF (.not. grid_obj%is_allocated) &
      & RETURN

    ! CALL message('destruct_cells','...')
    CALL destruct_cells(grid_obj%cells)
    ! CALL message('destruct_edges','...')
    CALL destruct_edges(grid_obj%edges)
    ! CALL message('destruct_verts','...')
    CALL destruct_verts(grid_obj%verts)

    IF (ASSOCIATED(grid_obj%child_grid_ids)) THEN
      ! CALL message('child_grid_ids','...')
      DEALLOCATE(grid_obj%child_grid_ids)
    ENDIF

    grid_obj%is_allocated  = .false.
    ! CALL message('deallocate_grid_object','is done')

  END SUBROUTINE deallocate_grid_object
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  !! Checks if a the grid object associated with the grid_id is active.
  !! If the grid object is not active the program stops with an error.
  SUBROUTINE check_active_grid_id (grid_id)
    INTEGER :: grid_id

    IF (grid_id > max_active_grids) THEN
      CALL finish('get_grid', 'grid_id > max_active_grids');
    ENDIF
    IF (.not. grid_is_active(grid_id)) THEN
      CALL finish('get_grid', 'grid is not active');
    ENDIF
  END SUBROUTINE check_active_grid_id
  !-----------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !   SUBROUTINE id_grid_idxs(in_grid_id)
  !>
  !! Forces grid entities idx to identity, idx(i) = i
  SUBROUTINE id_grid_idxs(in_grid_id)
    INTEGER, INTENT(in)  :: in_grid_id

    TYPE(t_grid), POINTER :: in_grid
    TYPE(t_grid_cells), POINTER :: cells
    TYPE(t_grid_edges), POINTER :: edges
    TYPE(t_grid_vertices), POINTER :: verts

    INTEGER :: no_of_input_cells, no_of_input_edges, no_of_input_verts
    INTEGER :: i

    in_grid => get_grid_object(in_grid_id)
    no_of_input_cells = in_grid%ncells
    no_of_input_edges = in_grid%nedges
    no_of_input_verts = in_grid%nverts
    verts=>in_grid%verts
    edges=>in_grid%edges
    cells=>in_grid%cells

    DO i=1,no_of_input_cells
      cells%idx(i) = i
    ENDDO
    DO i=1,no_of_input_edges
      edges%idx(i) = i
    ENDDO
    DO i=1,no_of_input_verts
      verts%idx(i) = i
    ENDDO

    RETURN

  END SUBROUTINE id_grid_idxs
  !-------------------------------------------------------------------------

END MODULE mo_local_grid
