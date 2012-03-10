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
!! This code has been tested up to a certain level. Defects and weaknesses,
!! which may be included in the code, do not establish any warranties by the
!! authors.
!! The authors do not make any warranty, express or implied, or assume any
!! liability or responsibility for the use, acquisition or application of this
!! software.
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
USE mo_communication,  ONLY: t_comm_pattern
USE mo_io_units,       ONLY: filename_max
USE mo_util_uuid,      ONLY: t_uuid
USE mo_util_subset,    ONLY: t_subset_range

IMPLICIT NONE

PRIVATE

CHARACTER(len=*), PARAMETER :: version = '$Id$'

! ! the following abstract data types are taken from mo_grid,
! ! mo_hierarchy and mo_base_geometry
!

PUBLIC :: t_patch
PUBLIC :: t_grid_cells
PUBLIC :: t_grid_edges
PUBLIC :: t_grid_vertices
PUBLIC :: t_phys_patch
  
!PUBLIC :: t_patch_ocean

PUBLIC :: t_tangent_vectors

! tangent vector class
TYPE t_tangent_vectors
  REAL(wp) :: v1
  REAL(wp) :: v2
END TYPE t_tangent_vectors

! !grid_cell class - corresponds to triangles

TYPE t_grid_cells
  
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
    &  center(:,:)

  ! area of triangle
  ! index1=nproma, index2=1,nblks_c
  REAL(wp), POINTER :: area(:,:)
  
  ! Coriolis parameter at cell centers
  ! index1=1,nproma, index2=1,nblks_c
  REAL(wp), ALLOCATABLE :: f_c(:,:)

  !----------------------------------
  ! cell geometry auxiliary variables
  ! the cartesian coordinates of the cell centers on the unit sphere
  TYPE(t_cartesian_coordinates), ALLOCATABLE ::  &
    &  cartesian_center(:,:)
  !----------------------------------
  
  ! refinement control flag
  ! index1=1,nproma, index2=1,nblks_c
  INTEGER, ALLOCATABLE :: refin_ctrl(:,:)

  ! list of start indices for each refin_ctrl level
  ! index1=min_rlcell,max_rlcell (defined in mo_impl_constants), index2=n_childdom
  INTEGER, ALLOCATABLE :: start_idx(:,:)

  ! list of end indices for each refin_ctrl level
  ! index1=min_rlcell,max_rlcell, index2=n_childdom
  INTEGER, ALLOCATABLE :: end_idx(:,:)

  ! list of start block for each refin_ctrl level
  ! index1=min_rlcell,max_rlcell, index2=n_childdom
  INTEGER, ALLOCATABLE :: start_blk(:,:)

  ! list of end block for each refin_ctrl level
  ! index1=min_rlcell,max_rlcell, index2=n_childdom
  INTEGER, ALLOCATABLE :: end_blk(:,:)

  ! Owner mask:
  ! For cells this is the same as decomp_domain(:,:)==0
  ! index1=nproma, index2=1,nblks_c
  LOGICAL, ALLOCATABLE :: owner_mask(:,:)

  ! The following is only used internally for the coupler
  ! and the icon_comm_lib
  INTEGER, ALLOCATABLE :: owner_local(:)
  
  ! The following is only used internally for domain decomposition

  INTEGER, ALLOCATABLE :: glb_index(:)
  INTEGER, ALLOCATABLE :: loc_index(:)

  ! Global array of owners
  INTEGER, ALLOCATABLE :: owner_g(:)

  ! Please note that the following array is only needed on local parent patches
  ! for storing the corresponding variable from nh_metrics.
  ! It is not allocated/deallocated with the regular patch (de)allocation routines
  ! nor does it need to be dumped/restored or subdivided in patch subdivision.
  REAL(wp), ALLOCATABLE :: ddqz_z_full(:,:,:)

  ! Domain decomposition flag:
  ! decomp_domain==0: inner domain, decomp_domain>0: boundary, decomp_domain<0: undefined
  ! 0=owned, 1=shared edge with owned, 2=shared vertex with ownded
  ! index1=nproma, index2=1,nblks_c
  INTEGER, POINTER :: decomp_domain(:,:)
  
  INTEGER, POINTER :: halo_level(:,:)! just points to the decomp_domain as a more accurate name
  
  ! define basic subsets
  TYPE(t_subset_range) :: owned         ! these are the owned entities
  TYPE(t_subset_range) :: in_domain     ! these are the computation domain entities, for cells = owned
  TYPE(t_subset_range) :: not_owned     ! these are all the halo entities
  TYPE(t_subset_range) :: not_in_domain ! for cells = not_owned  
  TYPE(t_subset_range) :: halos_with_edge_in_domain ! halo cells that at least one edge is in the
                                        ! computation domain
  
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

  !-------------------------------------------------
  ! edges geometry
  
  ! longitude & latitude of edge midpoint
  ! index=1,nproma, index2=1,nblks_e
  TYPE(t_geographical_coordinates), ALLOCATABLE ::  &
    &  center(:,:)

  ! normal to triangle edge
  ! index=1,nproma, index2=1,nblks_e
  TYPE(t_tangent_vectors), ALLOCATABLE ::  &
    &  primal_normal(:,:)

  ! Cartesian normal to triangle edge
  ! index=1,nproma, index2=1,nblks_e
  TYPE(t_cartesian_coordinates), ALLOCATABLE ::  &
    &  primal_cart_normal(:,:)

  ! Cartesian dual to triangle edge
  ! index=1,nproma, index2=1,nblks_e
  TYPE(t_cartesian_coordinates), ALLOCATABLE ::  &
    &  dual_cart_normal(:,:)

  ! normal to hexagon/pentagon edge
  ! index=1,nproma, index2=1,nblks_e
  TYPE(t_tangent_vectors), ALLOCATABLE ::  &
    &  dual_normal(:,:)

  ! normal to triangle edge, projected to the location of the neighbor cells
  ! index=1,nproma, index2=1,nblks_e, index3=1,2
  TYPE(t_tangent_vectors), ALLOCATABLE ::  &
    &  primal_normal_cell(:,:,:)

  ! tangent to triangle edge, projected to the location of the neighbor cells
  ! index=1,nproma, index2=1,nblks_e, index3=1,2
  TYPE(t_tangent_vectors), ALLOCATABLE ::  &
    &  dual_normal_cell(:,:,:)

  ! normal to triangle edge, projected to the location of the vertices
  ! index=1,nproma, index2=1,nblks_e, index3=1,4
  TYPE(t_tangent_vectors), ALLOCATABLE ::  &
    &  primal_normal_vert(:,:,:)

  ! tangent to triangle edge, projected to the location of the vertices
  ! index=1,nproma, index2=1,nblks_e, index3=1,4
  TYPE(t_tangent_vectors), ALLOCATABLE ::  &
    &  dual_normal_vert(:,:,:)

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
  TYPE(t_cartesian_coordinates), ALLOCATABLE ::  &
    &  cartesian_center(:,:)
  
  TYPE(t_cartesian_coordinates), ALLOCATABLE ::  &
    &  cartesian_dual_middle(:,:)
  
  ! Coriolis parameter at cell edges
  ! index1=1,nproma, index2=1,nblks_e
  REAL(wp), ALLOCATABLE :: f_e(:,:)

  ! refinement control flag
  ! index1=1,nproma, index2=1,nblks_e
  INTEGER, ALLOCATABLE :: refin_ctrl(:,:)

  ! list of start indices for each refin_ctrl level
  ! index1=min_rledge,max_rledge (defined in mo_impl_constants), index2=n_childdom
  INTEGER, ALLOCATABLE :: start_idx(:,:)

  ! list of end indices for each refin_ctrl level
  ! index1=min_rledge,max_rledge, index2=n_childdom
  INTEGER, ALLOCATABLE :: end_idx(:,:)

  ! list of start block for each refin_ctrl level
  ! index1=min_rledge,max_rledge, index2=n_childdom
  INTEGER, ALLOCATABLE :: start_blk(:,:)

  ! list of end block for each refin_ctrl level
  ! index1=min_rledge,max_rledge, index2=n_childdom
  INTEGER, ALLOCATABLE :: end_blk(:,:)

  ! Owner mask:
  ! For edges, this can not be derived from decomp_domain:
  ! edges at the border are assigned the PE with the bigger number
  ! index1=nproma, index2=1,nblks_e
  LOGICAL, ALLOCATABLE :: owner_mask(:,:)

  ! The following is only used internally for domain decomposition
  ! and communication
  INTEGER, ALLOCATABLE :: glb_index(:)
  INTEGER, ALLOCATABLE :: loc_index(:)

  ! Global array of owners
  INTEGER, ALLOCATABLE :: owner_g(:)
  
  ! The following is only used internally for the coupler
  ! and the icon_comm_lib
  INTEGER, ALLOCATABLE :: owner_local(:)

  ! Domain decomposition flag:
  ! decomp_domain==0: inner domain, decomp_domain>0: boundary, decomp_domain<0: undefined
  ! index1=nproma, index2=1,nblks_e
  ! 0=owned, 1=on owned cell=in domain, 2=exaclty one shared vertex with owned cells
  INTEGER, POINTER :: decomp_domain(:,:)  
  
  INTEGER, POINTER :: halo_level(:,:)! just points to the decomp_domain as a more accurate name
  
  ! define basic subsets
  TYPE(t_subset_range) :: owned         ! these are the owned entities
  TYPE(t_subset_range) :: in_domain     ! these are the computation domain entities,
                                        ! this includes all edges of the domain cells
  TYPE(t_subset_range) :: not_owned     ! these are all the halo entities
  TYPE(t_subset_range) :: not_in_domain ! all - in_domain 

END TYPE t_grid_edges

!  !grid_vertices class

TYPE t_grid_vertices

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
  TYPE(t_cartesian_coordinates), ALLOCATABLE ::  &
    &  cartesian(:,:)
  
  ! refinement control flag
  ! index1=1,nproma, index2=1,nblks_v
  INTEGER, ALLOCATABLE :: refin_ctrl(:,:)

  ! list of start indices for each refin_ctrl level
  ! index1=min_rlvert,max_rlvert (defined in mo_impl_constants), index2=n_childdom
  INTEGER, ALLOCATABLE :: start_idx(:,:)

  ! list of end indices for each refin_ctrl level
  ! index1=min_rlvert,max_rlvert, index2=n_childdom
  INTEGER, ALLOCATABLE :: end_idx(:,:)

  ! list of start block for each refin_ctrl level
  ! index1=min_rlvert,max_rlvert, index2=n_childdom
  INTEGER, ALLOCATABLE :: start_blk(:,:)

  ! list of end block for each refin_ctrl level
  ! index1=min_rlvert,max_rlvert, index2=n_childdom
  INTEGER, ALLOCATABLE :: end_blk(:,:)

  ! Owner mask:
  ! For verts, this can not be derived from decomp_domain:
  ! verts at the border are assigned the PE with the bigger number
  ! index1=nproma, index2=1,nblks_v
  LOGICAL, ALLOCATABLE :: owner_mask(:,:)

  ! The following is only used internally for domain decomposition

  INTEGER, ALLOCATABLE :: glb_index(:)
  INTEGER, ALLOCATABLE :: loc_index(:)

  ! Global array of owners
  INTEGER, ALLOCATABLE :: owner_g(:)
  
  ! The following is only used internally for the coupler
  ! and the icon_comm_lib
  INTEGER, ALLOCATABLE :: owner_local(:)
  
  ! Domain decomposition flag:
  ! decomp_domain==0: inner domain, decomp_domain>0: boundary, decomp_domain<0: undefined
  ! 0=owned, 1=on owned cell=in domain, 2=on level 1 cells
  ! index1=nproma, index2=1,nblks_v
  INTEGER, POINTER :: decomp_domain(:,:)

  INTEGER, POINTER :: halo_level(:,:) ! just points to the decomp_domain as a more accurate name
  
  ! define basic subsets
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

  ! cell type =3 or 6
  INTEGER :: cell_type
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

  !
  ! ! total number of cells, edges and vertices
  !
  INTEGER :: n_patch_cells
  INTEGER :: n_patch_edges
  INTEGER :: n_patch_verts
  !
  ! ! number of cells, edges and vertices in the global patch
  ! ! used only for debugging purposes
  !
  INTEGER :: n_patch_cells_g
  INTEGER :: n_patch_edges_g
  INTEGER :: n_patch_verts_g
  !
  ! ! values for the blocking
  !
  ! total...
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

  ! internal...
  ! LL: these numbers are the same as nblks_c, etc.
  !     See mo_subdivision, line 1100
  !     They are confusing and should be removed.
  ! number of blocks and chunk length in last block
  ! ... for the cells
  INTEGER :: nblks_int_c
  INTEGER :: npromz_int_c
  ! ... for the edges
  INTEGER :: nblks_int_e
  INTEGER :: npromz_int_e
  ! ... for the vertices
  INTEGER :: nblks_int_v
  INTEGER :: npromz_int_v

  !
  ! ! Mask and bathymetry for ocean patch
!  TYPE(t_patch_ocean) ::  &
!    &  patch_oce

  !
  ! ! grid information on the patch
  !
  TYPE(t_grid_cells) ::  &
    &  cells
  TYPE(t_grid_edges) ::  &
    &  edges
  TYPE(t_grid_vertices) ::  &
    &  verts

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
  TYPE(t_comm_pattern) :: comm_pat_gather_c
  TYPE(t_comm_pattern) :: comm_pat_gather_e
  TYPE(t_comm_pattern) :: comm_pat_gather_v

  ! Communication between local parent and its global counterpart,
  ! defined only on local parents.
  ! Please note that these communicate between
  ! p_patch_local_parent(jg) and p_patch(p_patch(jg)%parent_id)
  TYPE(t_comm_pattern) :: comm_pat_glb_to_loc_c
  TYPE(t_comm_pattern) :: comm_pat_glb_to_loc_e
  TYPE(t_comm_pattern) :: comm_pat_loc_to_glb_c_fbk
  TYPE(t_comm_pattern) :: comm_pat_loc_to_glb_e_fbk

END TYPE t_patch

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
  TYPE(t_comm_pattern) :: comm_pat_gather_c
  TYPE(t_comm_pattern) :: comm_pat_gather_e
  TYPE(t_comm_pattern) :: comm_pat_gather_v

END TYPE t_phys_patch

TYPE(t_patch), PUBLIC, TARGET, ALLOCATABLE :: p_patch(:), p_patch_local_parent(:)

! Please note: There is currently no means of determining the number
! of physical patches until they are actually assembled
! since this number is missing in the input file.
! Therefore p_phys_patch is declared with its maximum size.
! This shouldn't hurt since since TYPE t_phys_patch is small.

TYPE(t_phys_patch), PUBLIC, TARGET :: p_phys_patch(max_phys_dom)

!--------------------------------------------------------------------

END MODULE mo_model_domain

