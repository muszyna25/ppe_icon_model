!----------------------------------
!>
!!               The module contains the definition of the.
!!
!!               The module contains the definition of the
!!  grid and patch data structure that is
!!          actually used by the model.
!! It contains abstract data types that
!! are copies of those in
!! m_grid, m_hierarchy and m_base_geometry.
!!
!! @par Revision History
!! Initial version  by: Peter Korn,  MPI-M, Hamburg, June 2005
!! Cartesian coordinates type added by
!! Luca Bonaventura,  Polimi, Milano, October 2005
!!  Modification by Thomas Heinze (2006-02-21):
!!  - renamed m_modules to mo_modules
!! Component gitter of type 'patch' changed to grid, P. Korn (2006).
!! Modifications by Th.Heinze, DWD (2006-10-25):
!! - changed index to idx in TYPE declarations of grid_cells, grid_edges
!!    and grid_vertices
!! Modification by Peter Korn, MPI-M, (2006-11-23):
!!  - replacements in TYPE patch: ic by l2g_c, ie by l2g_e, iv by l2g_v,
!!    iic by g2l_c, iie by g2l_e, iiv by g2l_v
!!  - replaced edge_index by edge_idx
!!  - replaced vertex_index by vertex_idx
!!  - replaced cell_index by cell_idx
!!  - replaced neighbor_index by neighbor_idx
!!  - replaced child_index by child_idx
!! Modification by P. Ripodas (2007-01-31)
!!  - added system_orientation to TYPE grid_edges
!! Modification by Hui Wan, MPI-M, (2007-02-22):
!!  - type cartesian_coordinates and type geographical_coordinates
!!    moved to mo_math_utilities
!! Modification by A. Gassmann, MPI-M, (2007-04-03)
!!  - added list type
!!  - reorganized patch structure
!! Modification by Almut Gassman, MPI-M, (2007-04-13)
!!  - remove grid type
!!  - store external data (topography and land_sea_mask) in TYPE external data
!!  - k_list as variable name for integer lists
!! Modification by Jochen Foesrtner, DWD, (2008-07-16)
!!  - new fields in the derived type for the edges:
!!    grid_edges%primal_cart_normal (Cartesian normal to edge),
!!    grid_edges%quad_idx, grid_edges%quad_area and grid_edges%quad_orientation
!!    (indices of edges and area of the quadrilateral formed by two adjacent cells)
!!    up to now these new fields are initialized in mo_model_domain_import.f90
!!    rather than read from a grid/patch file.
!! Modification by Almut Gassmann, MPI-M, (2008-10-30)
!!  - Coriolis parameter is part of the external data
!! Modification by Almut Gassmann, MPI-M, (2009-03-06)
!!  - added 3d metrics part
!! Modification by Stephan Lorenz, MPI-M, (2010-02-02)
!!  - added type patch_ocean including 3-dim vertical topography (bathymetry)
!! Modification by Guenther Zaengl, DWD (2010-05-05)
!! - Moved metrics fields into NH state to avoid memory problems with MPI parallelization
!! Modification by Peter Korn, MPI-M, (2010-05-31)
!!  - added new data types patch_ocean for reconstruction process
!! Modification by Daniel Reinert, DWD (2010-07-21)
!!  - completely removed atmospheric external data from patch
!! Modification by Stephan Lorenz, MPI-M, (2010-10-26)
!!  - added data structures in patch_ocean for new "PtP" reconstruction process
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
MODULE mo_model_domain
  !-------------------------------------------------------------------------
  !
  !    ProTeX FORTRAN source: Style 2
  !    modified for ICON project, DWD/MPI-M 2006
  !
  !-------------------------------------------------------------------------
  !
  !
  !
  USE mo_kind
  USE mo_math_utilities, ONLY: t_geographical_coordinates, t_cartesian_coordinates
  USE mo_impl_constants, ONLY: max_dom, max_phys_dom
  USE mo_communication,  ONLY: t_comm_pattern, t_comm_gather_pattern
  USE mo_io_units,       ONLY: filename_max
  USE mo_util_uuid,      ONLY: t_uuid
  USE mo_grid_geometry_info, ONLY: t_grid_geometry_info
  USE mo_decomposition_tools, ONLY: t_grid_domain_decomp_info

  IMPLICIT NONE

  PRIVATE

  ! ! the following abstract data types are taken from mo_grid,

  PUBLIC :: t_patch, t_pre_patch
  PUBLIC :: t_grid_cells
  PUBLIC :: t_grid_edges
  PUBLIC :: t_grid_vertices
  PUBLIC :: t_phys_patch
  PUBLIC :: t_subset_range, t_subset_range_index, t_subset_indexed

  !PUBLIC :: t_patch_ocean
  PUBLIC :: t_patch_3D
  PUBLIC :: t_patch_vert

  PUBLIC :: t_tangent_vectors


  !----------------------------------------------------
  ! TODO Abstraction of the two types of subsets
  !----------------------------------------------------
  !> Defines a subset explicitly using an array of indexes
  TYPE :: t_subset_indexed

    INTEGER, ALLOCATABLE :: block(:)
    INTEGER, ALLOCATABLE :: idx(:)

    INTEGER :: size
    INTEGER :: recommended_stride  ! in case of parallelization/vectorization

    TYPE(t_patch), POINTER :: patch
    TYPE(t_patch_3D), POINTER :: patch_3D
    INTEGER :: entity_location ! on_cells, on_edges, on_verts

    INTEGER, POINTER :: vertical_levels(:,:)  ! if not null, points to the number of vertical levels array

    LOGICAL :: is_in_domain

    CHARACTER(len=32) :: name

  END TYPE t_subset_indexed
  !----------------------------------------------------

  !----------------------------------------------------
  !> Defines a subset in a range (in terms of blocks)
  TYPE :: t_subset_range

    INTEGER :: start_block
    INTEGER :: start_index
    INTEGER :: end_block
    INTEGER :: end_index
    INTEGER :: block_size

    INTEGER :: size

    TYPE(t_patch), POINTER :: patch
    INTEGER :: entity_location ! on_cells, on_edges, on_verts

    INTEGER :: max_vertical_levels
    INTEGER, POINTER :: vertical_levels(:,:)  ! if not null, points to the number of verticall levels array

    INTEGER :: no_of_holes ! the number of holes in the subset
    LOGICAL :: is_in_domain

    CHARACTER(len=32) :: name

  END TYPE t_subset_range
  !----------------------------------------------------


  !----------------------------------------------------
  !> Defines an index for a subset_range
  TYPE :: t_subset_range_index

    INTEGER, POINTER :: block ! the current block in the subset
    INTEGER, POINTER :: index ! the current index in the subset

    INTEGER :: start_index ! the current start index within the current block,
    INTEGER :: end_index   ! the current end index within the current block,

    TYPE(t_subset_range), POINTER :: subset_range

  END TYPE t_subset_range_index
  !----------------------------------------------------

  ! tangent vector class
  TYPE t_tangent_vectors
    REAL(wp) :: v1
    REAL(wp) :: v2
  END TYPE t_tangent_vectors

  ! !grid_cell class - corresponds to triangles

  TYPE t_grid_cells

    INTEGER :: max_connectivity
    ! number of edges connected to cell
    ! index1=1,nproma, index2=1,nblks_c
    INTEGER, ALLOCATABLE :: num_edges(:,:)

    ! line index of parent triangle:
    ! index1=1,nproma, index2=1,nblks_c
    INTEGER, ALLOCATABLE :: parent_idx(:,:)
    ! block index of parent triangle:
    ! index1=1,nproma, index2=1,nblks_c
    INTEGER, ALLOCATABLE :: parent_blk(:,:)
    ! parent child index, number of current cell in parent's child_idx/child_blk:
    ! index1=1,nproma, index2=1,nblks_c
    INTEGER, ALLOCATABLE :: pc_idx(:,:)

    ! line indices of child triangles:
    ! index1=1,nproma, index2=1,nblks_c, index3=1,4
    INTEGER, ALLOCATABLE :: child_idx(:,:,:)
    ! block indices of child triangles:
    ! index1=1,nproma, index2=1,nblks_c, index3=1,4
    INTEGER, ALLOCATABLE :: child_blk(:,:,:)
    ! domain ID of child triangles:
    ! index1=1,nproma, index2=1,nblks_c
    INTEGER, ALLOCATABLE :: child_id(:,:)
    ! physical domain ID of triangles
    ! (may differ from the "normal" domain ID in case of domain merging):
    ! index1=1,nproma, index2=1,nblks_c
    INTEGER, ALLOCATABLE :: phys_id(:,:)

    ! line indices of triangles next to each cell:
    ! index1=1,nproma, index2=1,nblks_c, index3=1,3
    INTEGER, ALLOCATABLE :: neighbor_idx(:,:,:)
    ! block indices of triangles next to each cell:
    ! index1=1,nproma, index2=1,nblks_c, index3=1,3
    INTEGER, ALLOCATABLE :: neighbor_blk(:,:,:)

    ! line indices of edges of triangle:
    ! index1=1,nproma, index2=1,nblks_c, index3=1,3
    INTEGER, ALLOCATABLE :: edge_idx(:,:,:)
    ! block indices of edges of triangle:
    ! index1=1,nproma, index2=1,nblks_c, index3=1,3
    INTEGER, ALLOCATABLE :: edge_blk(:,:,:)

    ! line indices of verts of triangle:
    ! index1=1,nproma, index2=1,nblks_c, index3=1,3
    INTEGER, ALLOCATABLE :: vertex_idx(:,:,:)
    ! block indices of verts of triangle:
    ! index1=1,nproma, index2=1,nblks_c, index3=1,3
    INTEGER, ALLOCATABLE :: vertex_blk(:,:,:)

    ! orientation of vectors normal to cell edges:
    ! index1=nproma, index2=1,nblks_c, index3=1,3
    REAL(wp), ALLOCATABLE :: edge_orientation(:,:,:)

    ! cell geometry

    ! longitude & latitude of centers of triangular cells
    ! index1=nproma, index2=1,nblks_c
    TYPE(t_geographical_coordinates), ALLOCATABLE ::  &
      & center(:,:)

    ! area of triangle
    ! index1=nproma, index2=1,nblks_c
    REAL(wp), POINTER :: area(:,:)

    ! Coriolis parameter at cell centers
    ! index1=1,nproma, index2=1,nblks_c
    REAL(wp), ALLOCATABLE :: f_c(:,:)

    !----------------------------------
    ! cell geometry auxiliary variables
    ! the cartesian coordinates of the cell centers on the unit sphere
    TYPE(t_cartesian_coordinates), POINTER ::  &
      & cartesian_center(:,:)
    !----------------------------------

    ! refinement control flag
    ! index1=1,nproma, index2=1,nblks_c
    INTEGER, ALLOCATABLE :: refin_ctrl(:,:)

    ! list of start indices for each refin_ctrl level
    ! index1=min_rlcell,max_rlcell (defined in mo_impl_constants), index2=n_childdom
    INTEGER, ALLOCATABLE :: start_idx(:,:) ! to be removed soon
    INTEGER, ALLOCATABLE :: start_index(:) ! revised implementation

    ! list of end indices for each refin_ctrl level
    ! index1=min_rlcell,max_rlcell, index2=n_childdom
    INTEGER, ALLOCATABLE :: end_idx(:,:) ! to be removed soon
    INTEGER, ALLOCATABLE :: end_index(:) ! revised implementation

    ! list of start block for each refin_ctrl level
    ! index1=min_rlcell,max_rlcell, index2=n_childdom
    INTEGER, ALLOCATABLE :: start_blk(:,:) ! to be removed soon
    INTEGER, ALLOCATABLE :: start_block(:) ! revised implementation

    ! list of end block for each refin_ctrl level
    ! index1=min_rlcell,max_rlcell, index2=n_childdom
    INTEGER, ALLOCATABLE :: end_blk(:,:) ! to be removed soon
    INTEGER, ALLOCATABLE :: end_block(:) ! revised implementation

    ! information on domain decomposition
    TYPE(t_grid_domain_decomp_info) :: decomp_info

    ! The owner when running the radiation
    ! only used with redistributed radiation
    INTEGER, POINTER :: radiation_owner(:)

    ! Please note that the following array is only needed on local parent patches
    ! for storing the corresponding variable from nh_metrics.
    ! It is not allocated/deallocated with the regular patch (de)allocation routines
    ! nor does it need to be dumped/restored or subdivided in patch subdivision.
    REAL(wp), ALLOCATABLE :: ddqz_z_full(:,:,:)

    ! define basic subsets
    TYPE(t_subset_range) :: ALL          ! these are the all valid entities, including all valid halos
    TYPE(t_subset_range) :: owned         ! these are the owned entities
    TYPE(t_subset_range) :: in_domain     ! these are the computation domain entities, for cells = owned
    TYPE(t_subset_range) :: not_owned     ! these are all the halo entities
    TYPE(t_subset_range) :: not_in_domain ! for cells = not_owned
    TYPE(t_subset_range) :: one_edge_in_domain ! cells with exactly one edge in domain. these are always halo cells

    INTEGER :: dummy_cell_block, dummy_cell_index     ! in case of existing a dummy cell for closing the
      ! local domain, these point to the dummy cell. This cell does not actually "exist"!
      ! =0 if no dummy cell

  END TYPE t_grid_cells

  ! !grid_edge class

  TYPE t_grid_edges

    ! line index of parent edge:
    ! index1=1,nproma, index2=1,nblks_e
    INTEGER, ALLOCATABLE :: parent_idx(:,:)
    ! block index of parent edge:
    ! index1=1,nproma, index2=1,nblks_e
    INTEGER, ALLOCATABLE :: parent_blk(:,:)
    ! parent child index, number of current edge in parent's child_idx/child_blk:
    ! index1=1,nproma, index2=1,nblks_e
    INTEGER, ALLOCATABLE :: pc_idx(:,:)

    ! line indices of child edges:
    ! index1=1,nproma, index2=1,nblks_e, index3=1,4
    INTEGER, ALLOCATABLE :: child_idx(:,:,:)
    ! block indices of child edges:
    ! index1=1,nproma, index2=1,nblks_e, index3=1,4
    INTEGER, ALLOCATABLE :: child_blk(:,:,:)
    ! domain ID of child edges:
    ! index1=1,nproma, index2=1,nblks_e
    INTEGER, ALLOCATABLE :: child_id(:,:)
    ! physical domain ID of edges
    ! (may differ from the "normal" domain ID in case of domain merging):
    ! index1=1,nproma, index2=1,nblks_e
    INTEGER, ALLOCATABLE :: phys_id(:,:)

    ! line indices of adjacent cells:
    ! index1=1,nproma, index2=1,nblks_e, index3=1,2
    INTEGER, ALLOCATABLE :: cell_idx(:,:,:)
    ! block indices of adjacent cells:
    ! index1=1,nproma, index2=1,nblks_e, index3=1,2
    INTEGER, ALLOCATABLE :: cell_blk(:,:,:)

    ! line indices of edge vertices:
    ! vertex indices 3 and 4 are the non-edge-aligned vertices of cells 1 and 2
    ! index1=1,nproma, index2=1,nblks_e, index3=1,4
    INTEGER, ALLOCATABLE :: vertex_idx(:,:,:)
    ! block indices of edge vertices:
    ! index1=1,nproma, index2=1,nblks_e, index3=1,4
    INTEGER, ALLOCATABLE :: vertex_blk(:,:,:)

    ! =1 if vector product of vector from vertex1 to vertex 2 (v2-v1) by vector
    ! from cell c1 to cell c2 (c2-c1) goes outside the sphere
    ! =-1 if vector product ...       goes inside  the sphere
    ! index=1,nproma, index2=1,nblks_e
    REAL(wp), ALLOCATABLE :: system_orientation(:,:)

    ! line indices of the  of the quadrilateral formed by two adjacent cells:
    ! index1=1,nproma, index2=1,nblks_e, index3=1,4
    INTEGER, ALLOCATABLE :: quad_idx(:,:,:)
    ! block indices of the  of the quadrilateral formed by two adjacent cells:
    ! index1=1,nproma, index2=1,nblks_e, index3=1,4
    INTEGER, ALLOCATABLE :: quad_blk(:,:,:)

    ! orientation of vectors normal to quad cell edges:
    ! index1=nproma, index2=1,nblks_e, index3=1,4
    REAL(wp), ALLOCATABLE :: quad_orientation(:,:,:)

    ! line indices of cells which share the two vertices
    ! of a given edge. These cells are subdivided into two classes:
    ! Cells which are neighbors of edge neighbor 1 and cells which are
    ! neighbors of edge neighbor 2. These cells are then numbered
    ! according to the number of the edge-vertex they share.
    ! index1=nproma, index2=1,nblks_e, index3=1,2 (cell), index4=1,2 (vert)
    INTEGER, ALLOCATABLE :: butterfly_idx(:,:,:,:)
    ! block indices of cells which share the two vertices
    ! of a given edge. These cells are subdivided into two classes:
    ! Cells which are neighbors of edge neighbor 1 and cells which are
    ! neighbors of edge neighbor 2. These cells are then numbered
    ! according to the number of the edge-vertex they share.
    ! index1=nproma, index2=1,nblks_e, index3=1,2 (cell), index4=1,2 (vert)
    INTEGER, ALLOCATABLE :: butterfly_blk(:,:,:,:)


    !-------------------------------------------------
    ! edges geometry

    ! longitude & latitude of edge midpoint
    ! index=1,nproma, index2=1,nblks_e
    TYPE(t_geographical_coordinates), ALLOCATABLE ::  &
      & center(:,:)

    ! normal to triangle edge
    ! index=1,nproma, index2=1,nblks_e
    TYPE(t_tangent_vectors), ALLOCATABLE ::  &
      & primal_normal(:,:)

    ! Cartesian normal to triangle edge
    ! index=1,nproma, index2=1,nblks_e
    TYPE(t_cartesian_coordinates), POINTER ::  &
      & primal_cart_normal(:,:)

    ! Cartesian dual to triangle edge
    ! index=1,nproma, index2=1,nblks_e
    TYPE(t_cartesian_coordinates), POINTER ::  &
      & dual_cart_normal(:,:)

    ! normal to hexagon/pentagon edge
    ! index=1,nproma, index2=1,nblks_e
    TYPE(t_tangent_vectors), ALLOCATABLE ::  &
      & dual_normal(:,:)

    ! normal to triangle edge, projected to the location of the neighbor cells
    ! index=1,nproma, index2=1,nblks_e, index3=1,2
    TYPE(t_tangent_vectors), ALLOCATABLE ::  &
      & primal_normal_cell(:,:,:)

    ! tangent to triangle edge, projected to the location of the neighbor cells
    ! index=1,nproma, index2=1,nblks_e, index3=1,2
    TYPE(t_tangent_vectors), ALLOCATABLE ::  &
      & dual_normal_cell(:,:,:)

    ! normal to triangle edge, projected to the location of the vertices
    ! index=1,nproma, index2=1,nblks_e, index3=1,4
    TYPE(t_tangent_vectors), ALLOCATABLE ::  &
      & primal_normal_vert(:,:,:)

    ! tangent to triangle edge, projected to the location of the vertices
    ! index=1,nproma, index2=1,nblks_e, index3=1,4
    TYPE(t_tangent_vectors), ALLOCATABLE ::  &
      & dual_normal_vert(:,:,:)

    ! length of triangle edge
    ! index=1,nproma, index2=1,nblks_e
    REAL(wp), POINTER :: primal_edge_length(:,:)

    ! inverse length of triangle edge
    ! index=1,nproma, index2=1,nblks_e
    REAL(wp), ALLOCATABLE :: inv_primal_edge_length(:,:)

    ! length of hexagon/pentagon edge
    ! index=1,nproma, index2=1,nblks_e
    REAL(wp), ALLOCATABLE :: dual_edge_length(:,:)

    ! inverse length of hexagon/pentagon edge
    ! index=1,nproma, index2=1,nblks_e
    REAL(wp), ALLOCATABLE :: inv_dual_edge_length(:,:)

    ! length of edge midpoint to vertex
    ! index=1,nproma, index2=1,nblks_e,2
    REAL(wp), ALLOCATABLE :: edge_vert_length(:,:,:)

    ! length of edge midpoint to cell center
    ! index=1,nproma, index2=1,nblks_e,2
    REAL(wp), ALLOCATABLE :: edge_cell_length(:,:,:)

    ! inverse distance between outer vertices of adjacent cells
    ! index=1,nproma, index2=1,nblks_e
    REAL(wp), ALLOCATABLE :: inv_vert_vert_length(:,:)

    ! area of the quadrilateral given by the adjacent vertices and centers
    ! index=1,nproma, index2=1,nblks_e
    REAL(wp), ALLOCATABLE :: area_edge(:,:)

    ! area of the quadrilateral formed by two cells adjacent to the edge
    ! index=1,nproma, index2=1,nblks_e
    REAL(wp), ALLOCATABLE :: quad_area(:,:)

    ! edge geometry auxiliary variables
    TYPE(t_cartesian_coordinates), POINTER ::  &
      & cartesian_center(:,:)

    TYPE(t_cartesian_coordinates), POINTER ::  &
      & cartesian_dual_middle(:,:)

    ! Coriolis parameter at cell edges
    ! index1=1,nproma, index2=1,nblks_e
    REAL(wp), ALLOCATABLE :: f_e(:,:)

    ! refinement control flag
    ! index1=1,nproma, index2=1,nblks_e
    INTEGER, ALLOCATABLE :: refin_ctrl(:,:)

    ! list of start indices for each refin_ctrl level
    ! index1=min_rledge,max_rledge (defined in mo_impl_constants), index2=n_childdom
    INTEGER, ALLOCATABLE :: start_idx(:,:) ! to be removed soon
    INTEGER, ALLOCATABLE :: start_index(:) ! revised implementation

    ! list of end indices for each refin_ctrl level
    ! index1=min_rledge,max_rledge, index2=n_childdom
    INTEGER, ALLOCATABLE :: end_idx(:,:) ! to be removed soon
    INTEGER, ALLOCATABLE :: end_index(:) ! revised implementation

    ! list of start block for each refin_ctrl level
    ! index1=min_rledge,max_rledge, index2=n_childdom
    INTEGER, ALLOCATABLE :: start_blk(:,:) ! to be removed soon
    INTEGER, ALLOCATABLE :: start_block(:) ! revised implementation

    ! list of end block for each refin_ctrl level
    ! index1=min_rledge,max_rledge, index2=n_childdom
    INTEGER, ALLOCATABLE :: end_blk(:,:) ! to be removed soon
    INTEGER, ALLOCATABLE :: end_block(:) ! revised implementation

    ! information on domain decomposition
    TYPE(t_grid_domain_decomp_info) :: decomp_info

    ! define basic subsets
    TYPE(t_subset_range) :: ALL          ! these are the all valid entities, including all valid halos
    TYPE(t_subset_range) :: owned         ! these are the owned entities
    TYPE(t_subset_range) :: in_domain     ! these are the computation domain entities,
    ! this includes all edges of the domain cells
    TYPE(t_subset_range) :: not_owned     ! these are all the halo entities
    TYPE(t_subset_range) :: not_in_domain ! all - in_domain

  END TYPE t_grid_edges

  !  !grid_vertices class

  TYPE t_grid_vertices

    INTEGER :: max_connectivity
    ! physical domain ID of verts
    ! (may differ from the "normal" domain ID in case of domain merging):
    ! index1=1,nproma, index2=1,nblks_v
    INTEGER, ALLOCATABLE :: phys_id(:,:)

    ! line indices of neighbor vertices:
    ! index1=1,nproma, index2=1,nblks_v, index3=1,6
    INTEGER, ALLOCATABLE :: neighbor_idx(:,:,:)
    ! block indices of neighbor vertices:
    ! index1=1,nproma, index2=1,nblks_v, index3=1,6
    INTEGER, ALLOCATABLE :: neighbor_blk(:,:,:)

    ! line indices of cells around each vertex:
    ! index1=1,nproma, index2=1,nblks_v, index3=1,6
    INTEGER, ALLOCATABLE :: cell_idx(:,:,:)
    ! block indices of cells around each vertex:
    ! index1=1,nproma, index2=1,nblks_v, index3=1,6
    INTEGER, ALLOCATABLE :: cell_blk(:,:,:)

    ! line indices of edges around a vertex:
    ! index1=1,nproma, index2=1,nblks_v, index3=1,6
    INTEGER, ALLOCATABLE :: edge_idx(:,:,:)
    ! block indices of edges around a vertex:
    ! index1=1,nproma, index2=1,nblks_v, index3=1,6
    INTEGER, ALLOCATABLE :: edge_blk(:,:,:)

    ! (xx)
    ! index1=1,nproma, index2=1,nblks_v, index3=1,6
    REAL(wp), ALLOCATABLE :: edge_orientation(:,:,:)

    ! number of edges connected to vertex
    ! index1=1,nproma, index2=1,nblks_v
    INTEGER, ALLOCATABLE :: num_edges(:,:)

    ! longitude & latitude of vertex:
    ! index1=1,nproma, index2=1,nblks_v
    TYPE(t_geographical_coordinates), ALLOCATABLE :: vertex(:,:)

    ! area of hexagon/pentagon of which vertex is center:
    ! index1=1,nproma, index2=1,nblks_v
    REAL(wp), POINTER :: dual_area(:,:)

    ! Coriolis parameter at cell vertices
    ! index1=1,nproma, index2=1,nblks_v
    REAL(wp), ALLOCATABLE :: f_v(:,:)

    ! vertex geometry auxiliary variables
    TYPE(t_cartesian_coordinates), POINTER ::  &
      & cartesian(:,:)

    ! refinement control flag
    ! index1=1,nproma, index2=1,nblks_v
    INTEGER, ALLOCATABLE :: refin_ctrl(:,:)

    ! list of start indices for each refin_ctrl level
    ! index1=min_rlvert,max_rlvert (defined in mo_impl_constants), index2=n_childdom
    INTEGER, ALLOCATABLE :: start_idx(:,:) ! to be removed soon
    INTEGER, ALLOCATABLE :: start_index(:) ! revised implementation

    ! list of end indices for each refin_ctrl level
    ! index1=min_rlvert,max_rlvert, index2=n_childdom
    INTEGER, ALLOCATABLE :: end_idx(:,:) ! to be removed soon
    INTEGER, ALLOCATABLE :: end_index(:) ! revised implementation

    ! list of start block for each refin_ctrl level
    ! index1=min_rlvert,max_rlvert, index2=n_childdom
    INTEGER, ALLOCATABLE :: start_blk(:,:) ! to be removed soon
    INTEGER, ALLOCATABLE :: start_block(:) ! revised implementation

    ! list of end block for each refin_ctrl level
    ! index1=min_rlvert,max_rlvert, index2=n_childdom
    INTEGER, ALLOCATABLE :: end_blk(:,:) ! to be removed soon
    INTEGER, ALLOCATABLE :: end_block(:) ! revised implementation

    ! information on domain decomposition
    TYPE(t_grid_domain_decomp_info) :: decomp_info

    ! define basic subsets
    TYPE(t_subset_range) :: ALL          ! these are the all valid entities, including all valid halos
    TYPE(t_subset_range) :: owned         ! these are the owned entities
    TYPE(t_subset_range) :: in_domain     ! these are the computation domain entities,
    ! this includes all edges of the domain cells
    TYPE(t_subset_range) :: not_owned     ! these are all the halo entities
    TYPE(t_subset_range) :: not_in_domain ! all - in_domain

  END TYPE t_grid_vertices


  ! !patch class

  TYPE t_patch

    !
    ! !  level in grid hierarchy on which patch lives
    !
    CHARACTER(LEN=filename_max) :: grid_filename
    !
    ! uuid of grid
    TYPE(t_uuid) :: grid_uuid
    !
    ! grid level
    INTEGER :: level
    !
    ! domain ID of current domain
    INTEGER :: id
    !
    ! indicator if current model domain is active
    LOGICAL :: ldom_active

    !-------------------------------------
    !> The grid domain geometry parameters
    ! cell type =3 or 6
    INTEGER :: cell_type

    INTEGER :: geometry_type

    TYPE(t_grid_geometry_info) :: geometry_info
    !-------------------------------------
    INTEGER :: boundary_depth_index  ! when is limited area grid, this is the  number of boundary levels based
                                     ! on the edge-connected cells (the fisrt level are cells that have
                                     ! at least one boundary edge, the next level are cells shering an edge with
                                     ! level 1 cells, etc
    
    !
    ! domain ID of parent domain
    INTEGER :: parent_id
    !
    ! child domain index of current domain as seen from parent domain
    ! In other words: I am the nth child of my parents (n=parent_child_index)
    INTEGER :: parent_child_index
    !
    ! list of child domain ID's
    INTEGER :: child_id(max_dom)
    !
    ! actual number of child domains
    INTEGER :: n_childdom
    !
    ! total number of child domains in the calling tree (over all nest levels)
    INTEGER :: n_chd_total
    !
    ! corresponding list of child domain ID's
    INTEGER :: child_id_list(max_dom)
    !
    ! maximum number of child domains
    INTEGER :: max_childdom
    !

    ! total number of locally allocated cells, edges and vertices
    INTEGER :: n_patch_cells
    INTEGER :: n_patch_edges
    INTEGER :: n_patch_verts

    ! total number of exist cells
!    INTEGER :: n_exist_cells
!    INTEGER :: n_cell_blocks
!    INTEGER :: last_cell_block_size

    ! in case of a dummy cell, we keep the blk and index of it
!    INTEGER :: dummy_cell_blk
!    INTEGER :: dummy_cell_idx


    !
    ! ! number of cells, edges and vertices in the global patch
    INTEGER :: n_patch_cells_g
    INTEGER :: n_patch_edges_g
    INTEGER :: n_patch_verts_g
    !
    ! ! values for the blocking
    !
    INTEGER :: alloc_cell_blocks  ! number of allocated cell blocks
    ! number of blocks and chunk length in last block
    ! ... for the cells
    INTEGER :: nblks_c
    INTEGER :: npromz_c
    ! ... for the edges
    INTEGER :: nblks_e
    INTEGER :: npromz_e
    ! ... for the vertices
    INTEGER :: nblks_v
    INTEGER :: npromz_v
    !
    ! ! vertical full and half levels
    !
    ! number of full and half levels
    INTEGER :: nlev
    INTEGER :: nlevp1
    !
    ! half level of parent domain (jg-1) that coincides
    ! with the upper margin of the current domain jg
    INTEGER :: nshift
    !
    ! total shift of model top with respect to global domain
    INTEGER :: nshift_total
    !
    ! the same information seen from the parent level (duplication needed to simplify flow control)
    INTEGER :: nshift_child


    !
    ! ! Mask and bathymetry for ocean patch
    !  TYPE(t_patch_ocean) ::  &
    !    &  patch_oce

    !
    ! ! grid information on the patch
    !
    TYPE(t_grid_cells) ::  &
      & cells
    TYPE(t_grid_edges) ::  &
      & edges
    TYPE(t_grid_vertices) ::  &
      & verts

    !
    ! communication patterns for parallelization
    !
    ! Boundary exchange within patches, defined on regular patches and local parents
    TYPE(t_comm_pattern) :: comm_pat_c
    TYPE(t_comm_pattern) :: comm_pat_c1 ! reduced communication pattern, only level-1 halo cells
    TYPE(t_comm_pattern) :: comm_pat_e
    TYPE(t_comm_pattern) :: comm_pat_v

    ! Interpolation for grid refinement, defined only on regular patches
    TYPE(t_comm_pattern) :: comm_pat_interpolation_c
    TYPE(t_comm_pattern) :: comm_pat_interpol_vec_grf(4)
    TYPE(t_comm_pattern) :: comm_pat_interpol_scal_grf(4)
    TYPE(t_comm_pattern) :: comm_pat_interpol_vec_ubc(4)
    TYPE(t_comm_pattern) :: comm_pat_interpol_scal_ubc(4)

    ! Gather complete patch to proc 0
    ! Useful only for regular patches (defined but unused on local parents)
    TYPE(t_comm_gather_pattern) :: comm_pat_gather_c
    TYPE(t_comm_gather_pattern) :: comm_pat_gather_e
    TYPE(t_comm_gather_pattern) :: comm_pat_gather_v

    ! Communication between local parent and its global counterpart,
    ! defined only on local parents.
    ! Please note that these communicate between
    ! p_patch_local_parent(jg) and p_patch(p_patch(jg)%parent_id)
    TYPE(t_comm_pattern) :: comm_pat_glb_to_loc_c
    TYPE(t_comm_pattern) :: comm_pat_glb_to_loc_e
    TYPE(t_comm_pattern) :: comm_pat_loc_to_glb_c_fbk
    TYPE(t_comm_pattern) :: comm_pat_loc_to_glb_e_fbk

    ! Halo comm patterns for the icon_commm_lib
    INTEGER ::  sync_cells_not_in_domain
    INTEGER ::  sync_cells_not_owned       ! = cells_not_in_domain
    INTEGER ::  sync_cells_one_edge_in_domain
    INTEGER ::  sync_edges_not_owned
    INTEGER ::  sync_edges_not_in_domain
    INTEGER ::  sync_verts_not_owned
    INTEGER ::  sync_verts_not_in_domain

    ! basic parallelization info for this patch
    LOGICAL :: compute_is_parallel  ! if true more than 1 processes are used
    LOGICAL :: is_in_parallel_test  ! if true there's one process that runs seq and compares the results to the parallel runs
    LOGICAL :: is_test_parallel_process  ! if true this process runs seq and compares the results to the parallel runs

    ! basic communicators concerning this patch
    !> the work universe of this patch, halo exchange, global sum,
    ! etc, takes place in this universe


    INTEGER :: work_communicator
    !> the work universe plus the sequential test process,
    ! in case it exists (otherwise the same as work_communicator)
    INTEGER :: parallel_test_communicator

    ! MPI communicator for this patch
    INTEGER :: comm
    !
    ! my index among the processors on this patch
    INTEGER :: rank
    !
    ! total number of processors on this patch
    INTEGER :: n_proc
    !
    ! global id of processor with rank 0 (within the working set p_comm_work)
    INTEGER :: proc0

  END TYPE t_patch
  !-----------------------------------------------------------------------------------
  ! !prepare patch class

  TYPE t_pre_patch

    !
    ! !  level in grid hierarchy on which patch lives
    !
    CHARACTER(LEN=filename_max) :: grid_filename
    !
    ! uuid of grid
    TYPE(t_uuid) :: grid_uuid
    !
    ! grid level
    INTEGER :: level
    !
    ! domain ID of current domain
    INTEGER :: id
    !
    ! indicator if current model domain is active
    LOGICAL :: ldom_active

    !-------------------------------------
    !> The grid domain geometry parameters
    ! cell type =3 or 6
    INTEGER :: cell_type

    INTEGER :: geometry_type

    TYPE(t_grid_geometry_info) :: geometry_info
    !-------------------------------------
    INTEGER :: boundary_depth_index  ! when is limited area grid, this is the  number of boundary levels based
                                     ! on the edge-connected cells (the fisrt level are cells that have
                                     ! at least one boundary edge, the next level are cells shering an edge with
                                     ! level 1 cells, etc
    
    !
    ! domain ID of parent domain
    INTEGER :: parent_id
    !
    ! child domain index of current domain as seen from parent domain
    ! In other words: I am the nth child of my parents (n=parent_child_index)
    INTEGER :: parent_child_index
    !
    ! list of child domain ID's
    INTEGER :: child_id(max_dom)
    !
    ! actual number of child domains
    INTEGER :: n_childdom
    !
    ! total number of child domains in the calling tree (over all nest levels)
    INTEGER :: n_chd_total
    !
    ! corresponding list of child domain ID's
    INTEGER :: child_id_list(max_dom)
    !
    ! maximum number of child domains
    INTEGER :: max_childdom
    !

    ! total number of locally allocated cells, edges and vertices
    INTEGER :: n_patch_cells
    INTEGER :: n_patch_edges
    INTEGER :: n_patch_verts

    ! total number of exist cells
!    INTEGER :: n_exist_cells
!    INTEGER :: n_cell_blocks
!    INTEGER :: last_cell_block_size

    ! in case of a dummy cell, we keep the blk and index of it
!    INTEGER :: dummy_cell_blk
!    INTEGER :: dummy_cell_idx


    !
    ! ! number of cells, edges and vertices in the global patch
    INTEGER :: n_patch_cells_g
    INTEGER :: n_patch_edges_g
    INTEGER :: n_patch_verts_g
    !
    ! ! values for the blocking
    !
    INTEGER :: alloc_cell_blocks  ! number of allocated cell blocks
    ! number of blocks and chunk length in last block
    ! ... for the cells
    INTEGER :: nblks_c
    INTEGER :: npromz_c
    ! ... for the edges
    INTEGER :: nblks_e
    INTEGER :: npromz_e
    ! ... for the vertices
    INTEGER :: nblks_v
    INTEGER :: npromz_v
    !
    ! ! vertical full and half levels
    !
    ! number of full and half levels
    INTEGER :: nlev
    INTEGER :: nlevp1
    !
    ! half level of parent domain (jg-1) that coincides
    ! with the upper margin of the current domain jg
    INTEGER :: nshift
    !
    ! total shift of model top with respect to global domain
    INTEGER :: nshift_total
    !
    ! the same information seen from the parent level (duplication needed to simplify flow control)
    INTEGER :: nshift_child


    !
    ! ! Mask and bathymetry for ocean patch
    !  TYPE(t_patch_ocean) ::  &
    !    &  patch_oce

    !
    ! ! grid information on the patch
    !
    TYPE(t_grid_cells) ::  &
      & cells
    TYPE(t_grid_edges) ::  &
      & edges
    TYPE(t_grid_vertices) ::  &
      & verts

    !
    ! communication patterns for parallelization
    !
    ! Boundary exchange within patches, defined on regular patches and local parents
    TYPE(t_comm_pattern) :: comm_pat_c
    TYPE(t_comm_pattern) :: comm_pat_c1 ! reduced communication pattern, only level-1 halo cells
    TYPE(t_comm_pattern) :: comm_pat_e
    TYPE(t_comm_pattern) :: comm_pat_v

    ! Interpolation for grid refinement, defined only on regular patches
    TYPE(t_comm_pattern) :: comm_pat_interpolation_c
    TYPE(t_comm_pattern) :: comm_pat_interpol_vec_grf(4)
    TYPE(t_comm_pattern) :: comm_pat_interpol_scal_grf(4)
    TYPE(t_comm_pattern) :: comm_pat_interpol_vec_ubc(4)
    TYPE(t_comm_pattern) :: comm_pat_interpol_scal_ubc(4)

    ! Gather complete patch to proc 0
    ! Useful only for regular patches (defined but unused on local parents)
    TYPE(t_comm_gather_pattern) :: comm_pat_gather_c
    TYPE(t_comm_gather_pattern) :: comm_pat_gather_e
    TYPE(t_comm_gather_pattern) :: comm_pat_gather_v

    ! Communication between local parent and its global counterpart,
    ! defined only on local parents.
    ! Please note that these communicate between
    ! p_patch_local_parent(jg) and p_patch(p_patch(jg)%parent_id)
    TYPE(t_comm_pattern) :: comm_pat_glb_to_loc_c
    TYPE(t_comm_pattern) :: comm_pat_glb_to_loc_e
    TYPE(t_comm_pattern) :: comm_pat_loc_to_glb_c_fbk
    TYPE(t_comm_pattern) :: comm_pat_loc_to_glb_e_fbk

    ! Halo comm patterns for the icon_commm_lib
    INTEGER ::  sync_cells_not_in_domain
    INTEGER ::  sync_cells_not_owned       ! = cells_not_in_domain
    INTEGER ::  sync_cells_one_edge_in_domain
    INTEGER ::  sync_edges_not_owned
    INTEGER ::  sync_edges_not_in_domain
    INTEGER ::  sync_verts_not_owned
    INTEGER ::  sync_verts_not_in_domain

    ! basic parallelization info for this patch
    LOGICAL :: compute_is_parallel  ! if true more than 1 processes are used
    LOGICAL :: is_in_parallel_test  ! if true there's one process that runs seq and compares the results to the parallel runs
    LOGICAL :: is_test_parallel_process  ! if true this process runs seq and compares the results to the parallel runs

    ! basic communicators concerning this patch
    !> the work universe of this patch, halo exchange, global sum,
    ! etc, takes place in this universe


    INTEGER :: work_communicator
    !> the work universe plus the sequential test process,
    ! in case it exists (otherwise the same as work_communicator)
    INTEGER :: parallel_test_communicator

    ! MPI communicator for this patch
    INTEGER :: comm
    !
    ! my index among the processors on this patch
    INTEGER :: rank
    !
    ! total number of processors on this patch
    INTEGER :: n_proc
    !
    ! global id of processor with rank 0 (within the working set p_comm_work)
    INTEGER :: proc0

  END TYPE t_pre_patch

  !-----------------------------------------------------------------------------------
  ! Description of physical patches
  TYPE t_phys_patch
    !
    ! domain ID of associated logical domain
    INTEGER :: logical_id
    !
    ! total number of cells, edges and vertices (always global)
    INTEGER :: n_patch_cells
    INTEGER :: n_patch_edges
    INTEGER :: n_patch_verts

    ! Gather the physical patch to proc 0
    TYPE(t_comm_gather_pattern) :: comm_pat_gather_c
    TYPE(t_comm_gather_pattern) :: comm_pat_gather_e
    TYPE(t_comm_gather_pattern) :: comm_pat_gather_v

  END TYPE t_phys_patch

  TYPE(t_patch), PUBLIC, TARGET, ALLOCATABLE :: p_patch(:)


  ! Definition of local parent patches
  ! For any given patch p_patch(jg) and jgp = p_patch(jg)%parent_id,
  ! p_patch_local_parent(jg) has the same resolution as p_patch(jgp)
  ! but it covers only the area of p_patch(jgp) which is covered by its child p_patch(jg)
  ! and it is divided in the same manner as p_patch(jg).
  ! Please note that p_patch_local_parent(1) is undefined if n_dom_start = 1

  TYPE(t_patch), PUBLIC, TARGET, ALLOCATABLE :: p_patch_local_parent(:)


  ! Please note: There is currently no means of determining the number
  ! of physical patches until they are actually assembled
  ! since this number is missing in the input file.
  ! Therefore p_phys_patch is declared with its maximum size.
  ! This shouldn't hurt since since TYPE t_phys_patch is small.

  TYPE(t_phys_patch), PUBLIC, TARGET :: p_phys_patch(max_phys_dom)

  !--------------------------------------------------------------------

  TYPE t_patch_vert

    !! The ocean uses z-coordinates in meters in the vertical.
    !! The following data are required:
    !!
    !! n_zlev: number of z-coordinate surfaces
    !! n_zlvp: number of intermediate levels (+1)
    !! n_zlvm: number of z-coordinate distances (-1)
    INTEGER :: n_zlev, n_zlvp, n_zlvm

    !! del_zlev_m: thickness (height) of elemental prism, defined as the
    !!             distance between top and bottom of elemental prism,
    !!             i.e. the distance between two intermediate z-coordinate
    !!             surfaces. These data are provided by the user, all other
    !!             vertical information is calculated from this array of
    !!             thicknesses.  Dimension: n_zlev
    REAL(wp), ALLOCATABLE :: del_zlev_m(:)

    !! the inverse of the above
    REAL(wp), ALLOCATABLE :: inv_del_zlev_m(:)

    !! zlev_m    : position of the vertical cell centers, i.e. below zero surface;
    !!             Numbering starts from surface and increases downwards to bottom.
    !!             Dimension: n_zlev
    !!             At these surfaces the horizontal velocities, vorticity, divergence
    !!             and scalar variables are evaluated.
    REAL(wp), ALLOCATABLE :: zlev_m(:)

    !! zlev_i    : vertical position of the UPPER BORDER of the vertical cell
    !!             i.e. the position of the top of elemental prisms.
    !!             Position of first surface is 0.
    !!             Dimension: n_zlvp = n_zlev + 1
    !!             The vertical velocities are evaluated at such surfaces.
    REAL(wp), ALLOCATABLE :: zlev_i(:)

    !! del_zlev_i: distance between two z-coordinate surfaces. The first is
    !!             the distance from the ocean surface = zlev_m(1)
    !!             Dimension: n_zlev
    REAL(wp), ALLOCATABLE :: del_zlev_i(:)

    ! To simplify the acess to the required information within these loops
    ! we store an cell and edge based version of the deepest ocean layer
    ! in column. dolic_e(edge1) and dolic_c(cell1) are identical if 'edge1'
    ! is one of the edges of 'cell1'.
    ! If the ocean bottom is flat dolic_c and dolic_e are identical and equal
    ! to the number of z-coodinate surfaces.
    !
    INTEGER, POINTER :: dolic_c(:,:)    ! index1=1,nproma, index2=1,nblks_c
    INTEGER, POINTER :: dolic_e(:,:)    ! index1=1,nproma, index2=1,nblks_e

    REAL(wp), POINTER :: ocean_area   (:)  ! global ocean area for each vertical level
    REAL(wp), POINTER :: ocean_volume (:)  ! global ocean volume for each vertical level

    REAL(wp), POINTER ::                &
      &  prism_thick_c(:,:,:),          & ! individual prism thickness at cells. Unit [m].
                                          ! This array includes the free surface. dimension: (nproma, n_zlev, nblks_c)
      &  prism_volume(:,:,:),           & ! individual prism volume (cells). Unit [m^3].
                                          ! This array includes the free surface. dimension: (nproma, n_zlev, nblks_c)
      &  prism_thick_e(:,:,:),          & ! individual prism thickness at edges. Unit [m].
                                          ! This array includes the free surface. dimension: (nproma, n_zlev, nblks_e)
      &  prism_thick_flat_sfc_c(:,:,:) ,& ! individual fluid column thickness at cells. Unit [m].
                                          ! This array assumes a flat surface. dimension: (nproma, n_zlev, nblks_c)
      &  prism_thick_flat_sfc_e(:,:,:) ,& ! individual fluid column thickness at edges. Unit [m].
                                          ! This array assumes a flat surface. dimension: (nproma, n_zlev, nblks_c)
      &  prism_center_dist_c(:,:,:),    & ! distance between prism centers at cells. Unit [m].
                                          ! dimension: (nproma,n_zlev, nblks_c)
      &  inv_prism_thick_c(:,:,:) ,     & ! inverse individual fluid column thickness at cells. Unit [m].
                                          ! This array assumes a flat surface dimension: (nproma, n_zlev, nblks_c)
      &  inv_prism_thick_e(:,:,:) ,     & ! inverse individual fluid column thickness at edges. Unit [m].
                                          ! This array assumes a flat surface. dimension: (nproma, n_zlev, nblks_e)
      &  inv_prism_center_dist_c(:,:,:),& ! inverse vertical distance between prism centers at cells. Unit [m].
                                          ! This array assumes a flat surface  dimension: (nproma, n_zlev, nblks_c)
      &  inv_prism_center_dist_e(:,:,:),& ! inverse vertical distance between prism centers at edges. Unit [m].
                                          ! This array assumes a flat surface dimension: (nproma, n_zlev, nblks_e)
      &  depth_CellMiddle(:,:,:),      & ! depth of the middle of the prism, update according to the current h
      &  depth_CellInterface(:,:,:)      ! depths at the interface (size is levels + 1)

  END TYPE t_patch_vert



  TYPE t_patch_3D

    TYPE(t_patch),     POINTER :: p_patch_2D(:)
    TYPE(t_patch_vert),POINTER :: p_patch_1D(:)

    ! land-sea-mask for ocean has 3 dimensions (the 2nd is the number of
    ! vertical levels)
    ! sea=-2, sea_boundary=-1, boundary (edges only)=0, land_boundary=1, land=2
    !
    ! land-sea-mask for cell centers
    ! index1=1,nproma, index2=1,n_zlev, index3=1,nblks_c
    INTEGER, POINTER :: lsm_c(:,:,:)

    ! land-sea-mask for edges
    ! index1=1,nproma, index2=1,n_zlev, index3=1,nblks_e
    INTEGER, POINTER :: lsm_e(:,:,:)

    ! surface land-sea-mask for cells, edges, vertices
    INTEGER, POINTER :: surface_cell_sea_land_mask(:,:)
    INTEGER, POINTER :: surface_edge_sea_land_mask(:,:)
    INTEGER, POINTER :: surface_vertex_sea_land_mask(:,:)

    ! To simply set land points to zero we store additional 3-dim wet points
    ! dimensions as in lsm_oce:
    REAL(wp), POINTER :: wet_c(:,:,:)  ! cells
    REAL(wp), POINTER :: wet_e(:,:,:)  ! edges
    ! For calculation of global sum and area including lsm the halo must be set to zero:
    REAL(wp), POINTER :: wet_halo_zero_c(:,:,:)  !  cells
    REAL(wp), POINTER :: wet_halo_zero_e(:,:,:)  !  edges

    ! For diagnosis like stream functions and area calculations we add surface arrays
    ! index1=1,nproma, index2=1,nblks_c
    INTEGER,  POINTER :: basin_c(:,:)  ! basin information Atlantic/Indian/Pacific
    INTEGER,  POINTER :: regio_c(:,:)  ! area information like tropical Atlantic etc.
    REAL(wp), POINTER :: bottom_thick_c(:,:)  ! individual bottom prism thickness at cells. Unit [m].
    REAL(wp), POINTER :: bottom_thick_e(:,:)  ! individual bottom prism thickness at edges. Unit [m].
    REAL(wp), POINTER :: column_thick_c(:,:)  ! individual column thickness at cells, no elevation. Unit [m].
    REAL(wp), POINTER :: column_thick_e(:,:)  ! individual column thickness at edges, no elevation. Unit [m].

  END TYPE t_patch_3D
  !--------------------------------------------------------------------

END MODULE mo_model_domain

