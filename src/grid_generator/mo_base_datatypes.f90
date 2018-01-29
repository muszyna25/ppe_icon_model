!>
!!  Contains the definition of basic used defined data types.
!!
!!  Contains the definition of basic used defined data types
!!  to store information about the grid in the INTERNAL grid representation
!!  used by the grid generator.  In particular,
!!  types for triangles, edges and vertices are defined.
!!  This representation is not convenient for
!!  use by the numerical model since edge and vertex entities are replicated
!!  for each grid cell.
!!
!! @par Revision History
!!  Initial version  by Luis Kornblueh (2004)
!!  Modifications and additions by Luca Bonaventura (2004-5).
!!  Modifications by Th.Heinze, DWD (2006-10-25):
!!  - changed index to idx in TYPE declarations of edge_info, vertex_info,
!!    triangle_info, grid_cells, grid_edges and grid_vertices
!! Modification by P. Ripodas (2007-02):
!!  - A new variable is added to TYPE edge (to store the system orientation
!!  of vectors  (v2-v1), (c2-c1), where v is for vertice and c is for cell)
!!  Modifications by Almut Gassmann (2007-03)
!!  - added some specifications for equal area subdivision in TYPE edge
!!  Modifications by Almut Gassmann (2007-04)
!!  - added area_edge to TYPE edge
!!  Modifications by Th.Heinze, DWD (2007-05-10):
!!  - moved definitions of dummy_e, dummy_c and dummy_v from mo_topology.f90 to
!!    this module
!!  Modification by Hui Wan (MPI-M 2007-08-06)
!!  - The attribute of dummy_e, dummy_c, and dummy_v
!!    was changed from TARGET,SAVE to POINTER.
!!  Modification by Almut Gassman, MPI-M (2009-01-14)
!!  - remove area_edge and add mid_edgex_len and edgex_orientation for cells
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_base_datatypes
  !-------------------------------------------------------------------------
  !
  !    ProTeX FORTRAN source: Style 2
  !    modified for ICON project, DWD/MPI-M 2004
  !
  !-------------------------------------------------------------------------
  !
  !
  !

  USE mo_kind,          ONLY: wp
  USE mo_math_types,    ONLY: t_cartesian_coordinates, t_tangent_vectors

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: t_edge_info, edge
  PUBLIC :: vertex, t_vertex_info
  PUBLIC :: t_triangle, t_triangle_info, t_triangle_reference
  PUBLIC :: dummy_e, dummy_v, dummy_c

  PUBLIC :: t_spheres

  ! Definition of several basic data types:
  !
  ! The definitions consits of several basic components. First of all
  ! is the definition of the edges and vertices. They are build by two major
  ! components. One is the informational part. The other part is the associated
  ! geometry. the whole triangle object does contain as well the topology and graph
  ! informations to its neighbors in all topological directions.
  !
  ! Note: in general it would be nice to have the pure topological information in
  !       mo_topology and the data types extended in mo_geometry, but that will
  !       only be available within the next Fortran standard.

  TYPE t_edge_info
    INTEGER :: idx
  END TYPE t_edge_info

  TYPE edge
    TYPE(t_edge_info) :: info
    TYPE(edge), POINTER :: parent => NULL() ! if edge is obtained from
    ! bisection of coarse grid edge,
    ! it has a parent, otherwise,
    ! it is parent of itself
    ! can be useful in the computation of fluxes
    ! at the coarse/fine interfaces and
    ! for the construction
    ! of topographies at different resolutions
    TYPE(t_cartesian_coordinates) :: edge_center

    TYPE(t_cartesian_coordinates) :: sc_pole !/ these both concern small circle
    REAL(wp) :: cos_rel_co_lat             !\ subdivision

    INTEGER ::  system_orientation   ! orientation of system
    ! vector (v2-v1), vector (c2-c1) at edge
    ! where v is for vertice and c for cell
    !
    !   MODIFIED to be directly of tangent_vectors type,
    !   in order to simplify re-implementation of their
    !   computation along the lines of the old generator
    !
    TYPE(t_tangent_vectors) :: edge_primal_normal
    TYPE(t_tangent_vectors) :: edge_dual_normal

    REAL(wp) :: edge_primal_arc  ! length of arc between
    ! the respective vertices
    REAL(wp) :: edge_dual_arc    ! length of arc between the
    ! respective triangle circum centers
    REAL(wp) :: edge_vert0_arc   ! arc length: edge midpoint to vertex0
    REAL(wp) :: edge_vert1_arc   ! arc length: edge midpoint to vertex1
    REAL(wp) :: edge_cell0_arc   ! arc length: edge midpoint to cell0
    REAL(wp) :: edge_cell1_arc   ! arc length: edge midpoint to cell1

    TYPE(t_triangle), POINTER :: triangle0 => NULL()     ! adjacent triangles
    TYPE(t_triangle), POINTER :: triangle1 => NULL()     !

    TYPE(vertex), POINTER :: vertex0 => NULL()    ! adjacent vertices
    TYPE(vertex), POINTER :: vertex1 => NULL()    !
  END TYPE edge

  TYPE t_vertex_info
    INTEGER :: idx
  END TYPE t_vertex_info

  TYPE vertex
    TYPE(t_vertex_info) :: info
    TYPE(t_cartesian_coordinates) :: vertex
    TYPE(t_cartesian_coordinates) :: vert_vel !velocity of vertex for spring dynamics
    REAL(wp) :: area_dual_cell   ! area of dual cell centered at vertex

    ! adjacent vertices (faces of dual cells)

    TYPE(vertex), POINTER :: neighbor0 => NULL()
    TYPE(vertex), POINTER :: neighbor1 => NULL()
    TYPE(vertex), POINTER :: neighbor2 => NULL()
    TYPE(vertex), POINTER :: neighbor3 => NULL()
    TYPE(vertex), POINTER :: neighbor4 => NULL()
    TYPE(vertex), POINTER :: neighbor5 => NULL()

    ! adjacent triangles (triangles that own that vertex)

    TYPE(t_triangle), POINTER :: triangle0 => NULL()     !
    TYPE(t_triangle), POINTER :: triangle1 => NULL()     !
    TYPE(t_triangle), POINTER :: triangle2 => NULL()     !
    TYPE(t_triangle), POINTER :: triangle3 => NULL()     !
    TYPE(t_triangle), POINTER :: triangle4 => NULL()     !
    TYPE(t_triangle), POINTER :: triangle5 => NULL()     !

    ! adjacent edges

    TYPE(edge), POINTER :: edge0 => NULL()     !
    TYPE(edge), POINTER :: edge1 => NULL()     !
    TYPE(edge), POINTER :: edge2 => NULL()     !
    TYPE(edge), POINTER :: edge3 => NULL()     !
    TYPE(edge), POINTER :: edge4 => NULL()     !
    TYPE(edge), POINTER :: edge5 => NULL()

  END TYPE vertex

  TYPE t_triangle_info
    INTEGER :: idx                ! absolute number of triangle on respective level
    INTEGER :: triangle_number    ! for checks with respect to development
    INTEGER :: parent_number = 0
    INTEGER :: orientation
    INTEGER :: configuration = 0  ! decimal number pattern containing parents
    INTEGER :: level = -1
  END TYPE t_triangle_info

  ! A very simplistic triangle type; pointer arrays are not allowed by
  ! Fortran95, so we get some few components numbered by their variable name.

  TYPE t_triangle
    TYPE(t_triangle_info) :: info
    REAL(wp) :: triangle_area
    TYPE(t_cartesian_coordinates) :: triangle_center

    TYPE(t_triangle), POINTER :: parent => NULL()

    ! The triangles are organized the following way:
    !
    !   0
    ! 1 2 3   for row-wise view
    !
    ! and so on for larger kroot (more subdivision)

    TYPE(t_triangle), POINTER :: sub_triangle0 => NULL()
    TYPE(t_triangle), POINTER :: sub_triangle1 => NULL()
    TYPE(t_triangle), POINTER :: sub_triangle2 => NULL()
    TYPE(t_triangle), POINTER :: sub_triangle3 => NULL()

    ! The orientation is with respect to the North Pole
    !
    ! UP:            DOWN:  base
    !         /\            _____
    ! left   /  \   right   \   /  left
    !       ------           \ /
    !        base
    !
    ! downward oriented triangles are rotated so that the same
    ! numbering scheme can be used!

    ! Note: neighbors should be counter-clockwise as well!

    TYPE(t_triangle), POINTER :: neighbor0 => NULL()    ! right
    TYPE(t_triangle), POINTER :: neighbor1 => NULL()    ! left
    TYPE(t_triangle), POINTER :: neighbor2 => NULL()    ! base

    ! Note: edgeX == neigborX->edge? (to be done ...); i.e., an edge
    ! is shared by two neighboring triangles

    TYPE(edge), POINTER :: edge0 => NULL()
    TYPE(edge), POINTER :: edge1 => NULL()
    TYPE(edge), POINTER :: edge2 => NULL()

    INTEGER :: edge0_orientation ! is 1 if the normal points outwards
    INTEGER :: edge1_orientation ! is -1 if the normal points inwards
    INTEGER :: edge2_orientation ! needed for divergence operator

    ! Note: vertexX == not_neigborX->vertex? (to be done ...); i.e., a vertex
    ! is shared by two not neighboring triangles

    TYPE(vertex), POINTER :: vertex0 => NULL()
    TYPE(vertex), POINTER :: vertex1 => NULL()
    TYPE(vertex), POINTER :: vertex2 => NULL()

  END TYPE t_triangle

  ! To go around restrictions in Fortran95 OO features, this constructs
  ! allows for generating arrays of pointers

  TYPE t_triangle_reference
    TYPE(t_triangle), POINTER :: t => NULL()
  END TYPE t_triangle_reference

  ! This defines a vector not a linked list of triangles on a single level.
  ! This would support as well refined patches - change to linked list?

  TYPE t_spheres

    INTEGER :: level
    INTEGER :: idx           = 0     ! required for recursive tree build
    INTEGER :: no_triangles  = 0
    INTEGER :: no_vertices   = 0
    INTEGER :: no_edges      = 0
    TYPE(t_triangle), POINTER :: ts(:) ! array of triangles on each level
    TYPE(edge)    , POINTER :: es(:) ! array of edges    on each level
    TYPE(vertex)  , POINTER :: vs(:) ! array of vertices on each level


  END TYPE t_spheres

  ! dummy instances of the above defined data types

  TYPE(t_triangle), POINTER :: dummy_c
  TYPE(vertex), POINTER :: dummy_v
  TYPE(edge), POINTER :: dummy_e

  !--------------------------------------------------------------------

END MODULE mo_base_datatypes

