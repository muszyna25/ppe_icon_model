!>
!!  Contains the routines which build the tree representation.
!!
!!  Contains the routines which build the tree representation
!!  of the grid based on refinement of the regular icosahedron.
!!  This tree representation is used by the grid generator and
!!  the data is stored in arrays indexed automatically according
!!  to an ordering resulting from this tree representation.
!!  Also the global indexing of all grid entities is builts
!!  along with the respective neighbourhood relationships and grid
!!  connectivities.
!!  The grid level -1 represents the original icosahedron with
!!  20 faces. It is stored in the array of triangle pointers
!!  <i>icosahedron</i>.
!!  Grid level 0 is the level where the initial subdivision of
!!  the later grid can be selected via the value of the
!!  parameter <i>kroot</i>. The grid is stored in the array of
!!  triangle pointers  <i>root</i>. The division can
!!  be selected arbitrarily, but useful for later use are
!!  divisions by 2, 3, and 5 only. The grid levels  from 2 onwards are
!!  divided always by two only (bisection). The achievable grid
!!  resolutions should allow for flexible enough selections.
!!  Numerical solvers are intended to work on level 1 and higher. Level 0
!!  does allow for flexible, fine-grained selections of grid resolutions.
!!
!!  The numbering scheme for vertices and edges is defined by the
!!  following rule:
!!
!!  Vertex 0 is the vertex at the tip of a triangle which would
!!  constitute the subtriangle 0 if subdivided. The order of the
!!  other vertices is counter clock-wise. Edge 0 is between the
!!  vertex 0 and 1 and so on as well ordered counter clock-wise.
!!  Care has to be taken with the up- and downward pointing of
!!  triangles, where upward  and downward refers to the orientation
!!  with respect to the North Pole.
!!
!! @par Revision History
!! Initial version  by:
!! @par
!!  Peter Sanders, MPI-I, Saarbruecken, January 2004
!! @par
!!  Luca Bonaventura, Luis Kornblueh, Uwe Schulzweida  MPI-M, Hamburg,
!!  January 2004
!! @par
!!  Luis Kornblueh, MPI-M, Hamburg, February 2004
!!          - change to compiling and running version
!!          - create level 0 subdivision, connect triangle pointers
!! Luis Kornblueh, MPI-M, Hamburg, March 2004
!!          - include parents
!!          - full working version of tree and graph representation
!!  Luca Bonaventura,  MPI-M, Hamburg, October 2004
!!          - include documentation and Protex headers
!!          - changed names to some routines for consistency
!!  Luca Bonaventura,  MPI-M, Hamburg, April 2005
!!          - complete neighborhood relationships and connectivities
!!  Modifications by Th.Heinze, DWD (2006-10-25):
!!  - changed index to idx in TYPE declarations of edge_info, vertex_info,
!!    triangle_info, grid_cells, grid_edges and grid_vertices
!!  Modifications by Th.Heinze, DWD (2007-05-07):
!!  - adapted init_topology for grid and graph generator:
!!  - created init_topology_graph and init_topology_grid, deleted init_topology
!!  Modifications by Th.Heinze, DWD (2007-05-09):
!!  - created read_graph
!!  Modifications by Th.Heinze, DWD (2007-05-10):
!!  - created write_graph
!!  - moved definitions of dummy_e, dummy_c and dummy_v to mo_base_datatypes.f90
!!  Modification by Hui Wan (MPI-M 2007-08-06)
!!  - the attribute of dummy_e, dummy_c, and dummy_v in mo_base_datatypes
!!    was changed from TARGET,SAVE to POINTER. Accordingly, allocation and
!!    deallocation of these items were added to subroutine build_full_graph.
!!  - deallocation of spheres_on_levels(base_level)%ts(j)%vertex0/1/2 removed.
!!  Modifications by Th.Heinze, DWD (2007-08-07):
!!  - allocated dummy_e, dummy_c and dummy_v according to the POINTER attribute
!!  - moved the allocation from build_full_graph to init_topology_grid resp.
!!    init_topology_graph
!!  - deallocated of spheres_on_levels(base_level)%ts(j)%vertex0/1/2 again.
!!    (seems necessary on DWD's Linux machine using g95 compiler)
!!  Modifications by G. Zaengl, DWD (2008-12-08):
!!  - removal of unnecessary computations for efficiency improvement
!!  Modifications by Almut Gassmann, MPI-M (2009-01-15)
!!  - double periodic planar grid option
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
MODULE mo_topology
  !
  !
  !
  !
  !
  ! Note: The grid level -1 is dedicated to the original icosahedron with
  !       20 faces.
  !       Grid level 0 is the level where the initial subdivision of
  !       the later grid can be selected via kroot. The division can
  !       be selected arbitrarily, but useful for later use are
  !       divisions by 2, 3, and 5 only. The following grid levels are
  !       divided always by two only (bisection). The achievable grid
  !       resolutions should allow for flexible enough selections.
  !
  !       The numbering scheme for vertices and edges is defined by the
  !       following rule:
  !
  !                            vertex 0
  !                               /\
  !                neighbor 0    /  \
  !                     edge 0  /    \ edge 2
  !                            /      \  neighbor 2
  !                           /        \
  !                 vertex 1 ------------ vertex 2
  !                             edge 1
  !                            neighbor 1
  !
  !                            neighbor 1
  !                             edge 1
  !                 vertex 2 ------------ vertex 1
  !                           \        /
  !                     edge 2 \      / edge 0
  !                  neighbor 2 \    /   neighbor 0
  !                              \  /
  !                               \/
  !                            vertex 0
  !
  !       vertex 0 is the vertex at the tip of a triangle which would
  !       constitute the subtriangle 0 if subdivided. The order of the
  !       other vertices is counter clock-wise. edge 0 is between the
  !       vertex 0 and 1 and so on as well ordered counter clock-wise.
  !       Care has to be taken with the up- and downward pointing of
  !       triangles, where upward  and downward refers to the orientation
  !       with respect to the North Pole.


  USE mo_base_datatypes, ONLY: t_triangle, edge, vertex, t_triangle_reference, &
    & t_spheres, dummy_e, dummy_c, dummy_v

  USE mo_exception,      ONLY: message, message_text, finish
  USE mo_io_units,       ONLY: nerr
  USE mo_io_graph,       ONLY: input_graph,output_graph

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: version = '$Id$'

  PUBLIC :: init_topology_grid, init_topology_graph
  PUBLIC :: generate_tree, generate_graph, read_graph, write_graph
  PUBLIC :: destruct_topology_grid, destruct_topology_graph
  PUBLIC :: icosahedron, root, spheres_on_levels
  PUBLIC :: base_level, number_of_levels, kroot, up_down

  ! kroot is defining the initial sub division of
  ! the icosahedron. Values to be used are 2, 3, or 5.
  ! The algorithms implemented allow nevertheless for
  ! arbitrary subdivision.

  INTEGER :: kroot

  ! Bit pattern usage for comparison of special cases (start indexing by 1
  ! as required by Fortran95 intrinsics):

  INTEGER, PARAMETER :: up     = 1
  INTEGER, PARAMETER :: down   = 2
  INTEGER, PARAMETER :: pole   = 3
  INTEGER, PARAMETER :: main   = 4
  INTEGER, PARAMETER :: border = 5
  INTEGER, PARAMETER :: left   = 6
  INTEGER, PARAMETER :: right  = 7

  INTEGER, PARAMETER :: up_down = IOR(IBSET(0,up),IBSET(0,down))

  ! This is to define level -1, contains special handling for setting up
  ! that is the reason for the reference usage here and for the following
  ! root.

  TYPE(t_triangle_reference), ALLOCATABLE :: icosahedron(:)

  ! Reference of the root triangles on level 0 with respect to the triangles
  ! of the original icosahedron, first index represents the level 0 local index
  ! within a main triangle of the original icosahedron, the second
  ! index the originating triangle of the icosahedron

  TYPE(t_triangle_reference), ALLOCATABLE, TARGET :: root(:,:)

  ! Array of triangles arrays on each refinement level (index represent level)

  TYPE(t_spheres), ALLOCATABLE, TARGET :: spheres_on_levels(:)

  ! Refinement level of original icosahedron

  INTEGER, PARAMETER :: base_level = -1

  ! The number of levels is input parameter for the data structure generation.
  ! Initiliazed with the base_level.

  INTEGER :: number_of_levels = base_level

CONTAINS

  !-------------------------------------------------------------------------
  !
  !-------------------------------------------------------------------------
  !
  !

  !---------------------------------------------------------------------------
  !
  ! generate tree representation and allocate container array(s)

  !>
  !!               Allocation of arrays needed to store the original icosahedron,.
  !!
  !!               Allocation of arrays needed to store the original icosahedron,
  !!               the root grid obtained with the first refinements
  !!               and all the refined grids. The variable <i>levels</i>
  !!               denotes the number of refinements to be performed below
  !!               level 0, which corresponds to the root  grid. The original
  !!               icosahedron corresponds to level -1.
  !!               The variable <i>root_subdivision</i> gives the number of times
  !!               the edge of a main triangle in the original icosahedron is
  !!               subdivided in order to generate the root grid. The indices of
  !!               triangles on the root grid are stored in the array <i>root</i>.
  !!               These triangles are identified by a local index
  !!               (first index of <i>root</i> ) and by index of the corresponding
  !!               main triangle belonging to the icosahedron (second index of
  !!               <i>root</i>).
  !!               \\newline
  !!               If we want to generate a plane a pseudo icosahedron consists of
  !!               8 triangles in the form
  !!                    __  __
  !!                  /\\0 /\\1 /
  !!                 /2_\\/3_\\/
  !!                /\\ 4/\\ 5/
  !!               /6_\\/7_\\/
  !!
  !! @par Revision History
  !!  Luis Kornblueh, MPI-M, Hamburg, (2004-03)
  !!  Modification by Thomas Heinze, DWD, (2007-05-07):
  !!  - renamed init_topology to init_topology_graph and adapted structure to
  !!    graph generator
  !!  Modification by Almut Gassmann, MPI-M, (2009-01-08):
  !!  - added facility to generate a planar grid
  !!
  SUBROUTINE init_topology_graph (k_levels, k_root_subdivision, lplane)
    !

    ! number of levels to be generated below the
    ! root grid (level 0)

    INTEGER, INTENT(in) :: k_levels

    ! number of times the edge of
    ! a main triangle of the original icosahedron (level -1)
    ! is subdivided to obtain the root grid (level 0)

    INTEGER, INTENT(in) :: k_root_subdivision

    LOGICAL, INTENT(in) :: lplane !needed for generating a plane

    INTEGER :: i_not, i_nov, i_noe
    INTEGER :: ji, istat, n_faces_m1

    !-------------------------------------------------

    ! The number of faces minus 1 of the icosahedron
    IF (.not. lplane) THEN
      n_faces_m1 = 19  ! icosahedron
    ELSE
      n_faces_m1 = 7   ! pseudo icosahedron on a plane
    ENDIF

    ! allocations added by Hui Wan (MPI-M, 2007-08-06)
    ! moved from build_full_graph here by Thomas Heinze, DWD, 2007-08-07

    ALLOCATE(dummy_c, stat=istat)
    IF( istat/=0 ) THEN
      CALL finish('init_topology_graph','allocation of dummy_c failed')
    ENDIF

    ALLOCATE(dummy_e, stat=istat)
    IF( istat/=0 ) THEN
      CALL finish('init_topology_graph','allocation of dummy_e failed')
    ENDIF

    ALLOCATE(dummy_v, stat=istat)
    IF( istat/=0 ) THEN
      CALL finish('init_topology_graph','allocation of dummy_v failed')
    ENDIF

    ! now set the final numbers of levels and kroot

    number_of_levels = k_levels
    kroot = k_root_subdivision

    ! The icosahedron itself
    ALLOCATE (icosahedron(0:n_faces_m1), stat=istat)

    IF (istat >0) THEN
      WRITE (message_text,'(a)')    'Problem in allocating icosahedron '
      CALL finish ('init_topology', message_text)
    ENDIF

    ! root of 'tree', starting on level 0 - not on the original icosahedron
    !
    ! it does the allocation of the root matrix, containing the topology of the
    ! root grid; triangles on the root grid are identified by a local index
    ! (first index of root) and by index of the corresponding main triangle
    ! belonging to the icosahedron

    ALLOCATE (root(0:kroot*kroot-1,0:n_faces_m1), stat=istat)

    IF (istat >0) THEN
      WRITE (message_text,'(a)')    'Problem in allocating root'
      CALL finish ('init_topology', message_text)
    ENDIF

    ! Container for all spherical grids on all levels

    ALLOCATE (spheres_on_levels(base_level:number_of_levels), stat=istat)

    IF (istat >0) THEN
      WRITE (message_text,'(a)')    'Problem in allocating spheres_on_levels'
      CALL finish ('init_topology', message_text)
    ENDIF

    !--------------------------------------------------------------------------
    ! level base_level does contain the original icosahedron information

    IF (.not. lplane) THEN
      i_not = 20 ! triangles (faces)
      i_nov = 12 ! vertices
      i_noe = 30 ! edges
    ELSE ! lplane ==.TRUE.
      i_not = 8  ! triangles (faces)
      i_nov = 4  ! vertices
      i_noe = 12 ! edges
    ENDIF

    spheres_on_levels(base_level)%level = base_level

    ! required for the recursive generation of topology information
    spheres_on_levels(base_level)%idx         = 0

    spheres_on_levels(base_level)%no_triangles  = i_not
    spheres_on_levels(base_level)%no_vertices   = i_nov
    spheres_on_levels(base_level)%no_edges      = i_noe

    ALLOCATE(spheres_on_levels(base_level)%ts(0:i_not-1), stat=istat)

    IF (istat >0) THEN
      WRITE (message_text,'(a)')'Problem in allocating spheres_on_levels%ts, base_level'
      CALL finish ('init_topology', message_text)
    ENDIF

    ALLOCATE(spheres_on_levels(base_level)%es(0:i_noe-1), stat=istat)

    IF (istat >0) THEN
      WRITE (message_text,'(a)')'Problem in allocating spheres_on_levels%es, base_level'
      CALL finish ('init_topology', message_text)
    ENDIF

    ALLOCATE(spheres_on_levels(base_level)%vs(0:i_nov-1), stat=istat)

    IF (istat >0) THEN
      WRITE (message_text,'(a)')'Problem in allocating spheres_on_levels%vs, base_level'
      CALL finish ('init_topology', message_text)
    ENDIF

    WRITE (message_text,'(a,i3,a,i9)') &
      & 'Number of triangular (primal) cells on grid level ', base_level, &
      & ', is ', spheres_on_levels(base_level)%no_triangles
    CALL message ('',TRIM(message_text))

    WRITE (message_text,'(a,i3,a,i9)') &
      & 'Number of hexagonal  (dual)   cells on grid level ', base_level, &
      & ', is ', spheres_on_levels(base_level)%no_vertices
    CALL message ('',TRIM(message_text))

    !   !--------------------------------------------------------------------------
    !   ! setup usable grid levels - allocate those arrays
    !   ! these allocated triangles get connected in build_tree

    DO ji = 0, number_of_levels

      IF (.not. lplane) THEN
        i_not = (20*kroot*kroot)*4**ji
        i_nov = (10*kroot*kroot)*4**ji + 2
        i_noe = (30*kroot*kroot)*4**ji
      ELSE
        i_not = ( 8*kroot*kroot)*4**ji
        i_nov = ( 4*kroot*kroot)*4**ji
        i_noe = (12*kroot*kroot)*4**ji
      ENDIF

      spheres_on_levels(ji)%level = ji

      ! required for the recursive generation of topology information

      spheres_on_levels(ji)%idx         = 0

      spheres_on_levels(ji)%no_triangles  = i_not
      spheres_on_levels(ji)%no_vertices   = i_nov
      spheres_on_levels(ji)%no_edges      = i_noe

      ALLOCATE(spheres_on_levels(ji)%ts(0:i_not-1), stat=istat)
      IF (istat >0) THEN
        WRITE (message_text,'(a,i2)')'Problem in allocating spheres_on_levels%ts, level ',ji
        CALL finish ('init_topology', message_text)
      ENDIF

      ALLOCATE(spheres_on_levels(ji)%es(0:i_noe-1), stat=istat)
      IF (istat >0) THEN
        WRITE (message_text,'(a,i2)')'Problem in allocating spheres_on_levels%es, level ',ji
        CALL finish ('init_topology', message_text)
      ENDIF

      ALLOCATE(spheres_on_levels(ji)%vs(0:i_nov-1), stat=istat)
      IF (istat >0) THEN
        WRITE (message_text,'(a,i2)')'Problem in allocating spheres_on_levels%ts, level ',ji
        CALL finish ('init_topology', message_text)
      ENDIF

      WRITE (message_text,'(a,i3,a,i9)')                                     &
        & 'Number of triangular (primal) cells on grid level ', ji,       &
        & ', is ', spheres_on_levels(ji)%no_triangles
      CALL message ('',TRIM(message_text))

      WRITE (message_text,'(a,i3,a,i9)')                                     &
        & 'Number of hexagonal  (dual)   cells on grid level ', ji,      &
        & ', is ', spheres_on_levels(ji)%no_vertices
      CALL message ('',TRIM(message_text))

    END DO

  END SUBROUTINE init_topology_graph

  !-------------------------------------------------------------------------
  !
  !

  !---------------------------------------------------------------------------
  !
  ! generate tree representation and allocate container array(s)

  !>
  !!  Allocation of arrays needed to store the original icosahedron, the.
  !!
  !!  Allocation of arrays needed to store the original icosahedron, the
  !!  root grid obtained with the first refinements
  !!   and all the refined grids. The variable   <i>levels</i>
  !!  denotes the number of refinements to be performed below
  !!  level 0, which corresponds to the root  grid. The original icosahedron
  !!  corresponds to level -1.
  !!  The variable   <i>root_subdivision</i> gives the number of times
  !!  the edge of a main triangle in the original icosahedron is subdivided
  !!  in order to generate the root grid. The indices of triangles on the root grid
  !!  are stored in the array   <i>root</i>. These triangles are identified by a local index
  !!   (first index of  <i>root</i> ) and by index of the corresponding main triangle
  !!   belonging to the icosahedron (second index of   <i>root</i>).
  !!
  !! @par Revision History
  !!  Luis Kornblueh, MPI-M, Hamburg, (2004-03)
  !!  Modification by Thomas Heinze, DWD, (2007-05-07):
  !!  - renamed init_topology to init_topology_graph and adapted structure to
  !!    graph generator
  !!
  SUBROUTINE init_topology_grid (k_levels, k_root_subdivision, lplane)
    !

    ! number of levels to be generated below the
    ! root grid (level 0)

    INTEGER, INTENT(in) :: k_levels

    ! number of times the edge of
    ! a main triangle of the original icosahedron (level -1)
    ! is subdivided to obtain the root grid (level 0)

    INTEGER, INTENT(in) :: k_root_subdivision

    LOGICAL, INTENT(in) :: lplane

    INTEGER :: i_not, i_nov, i_noe
    INTEGER :: ji, istat, n_faces_m1

    !-------------------------------------------------

    ! allocations added by Hui Wan (MPI-M, 2007-08-06)
    ! moved from build_full_graph here by Thomas Heinze, DWD, 2007-08-07

    ALLOCATE(dummy_c, stat=istat)
    IF( istat/=0 ) THEN
      CALL finish('init_topology_grid','allocation of dummy_c failed')
    ENDIF

    ALLOCATE(dummy_e, stat=istat)
    IF( istat/=0 ) THEN
      CALL finish('init_topology_grid','allocation of dummy_e failed')
    ENDIF

    ALLOCATE(dummy_v, stat=istat)
    IF( istat/=0 ) THEN
      CALL finish('init_topology_grid','allocation of dummy_v failed')
    ENDIF

    ! now set the final numbers of levels and kroot

    number_of_levels = k_levels
    kroot = k_root_subdivision

    IF (.not.lplane) THEN
      n_faces_m1 = 19
    ELSE
      n_faces_m1 = 7
    ENDIF

    ! The icosahedron itself

    ALLOCATE (icosahedron(0:n_faces_m1), stat=istat)

    IF (istat >0) THEN
      WRITE (message_text,'(a)')    'Problem in allocating icosahedron '
      CALL message ('',TRIM(message_text))
      CALL finish ('init_topology', message_text)
    ENDIF

    ! root of 'tree', starting on level 0 - not on the original icosahedron
    !
    ! it does the allocation of the root matrix, containing the topology of the
    ! root grid; triangles on the root grid are identified by a local index
    ! (first index of root) and by index of the corresponding main triangle
    ! belonging to the icosahedron

    ALLOCATE (root(0:kroot*kroot-1,0:n_faces_m1), stat=istat)

    IF (istat >0) THEN
      WRITE (message_text,'(a)')    'Problem in allocating root'
      CALL message ('',TRIM(message_text))
      CALL finish ('init_topology', message_text)
    ENDIF

    ! Container for all spherical grids on all levels

    ALLOCATE (spheres_on_levels(base_level:number_of_levels), stat=istat)

    IF (istat >0) THEN
      WRITE (message_text,'(a)')    'Problem in allocating spheres_on_levels'
      CALL message ('',TRIM(message_text))
      CALL finish ('init_topology', message_text)
    ENDIF

    !--------------------------------------------------------------------------
    ! level base_level does contain the original icosahedron information

    IF (.not. lplane) THEN
      i_not = 20 ! triangles (faces)
      i_nov = 12 ! vertices
      i_noe = 30 ! edges
    ELSE
      i_not =  8 ! triangles (faces)
      i_nov =  4 ! vertices
      i_noe = 12 ! edges
    ENDIF

    spheres_on_levels(base_level)%level = base_level

    ! required for the recursive generation of topology information
    spheres_on_levels(base_level)%idx         = 0

    spheres_on_levels(base_level)%no_triangles  = i_not
    spheres_on_levels(base_level)%no_vertices   = i_nov
    spheres_on_levels(base_level)%no_edges      = i_noe

    ALLOCATE(spheres_on_levels(base_level)%ts(0:i_not-1), stat=istat)

    IF (istat >0) THEN
      WRITE (message_text,'(a)')    'Problem in allocating spheres_on_levels%ts, base_level'
      CALL message ('',TRIM(message_text))
      CALL finish ('init_topology', message_text)
    ENDIF

    ALLOCATE(spheres_on_levels(base_level)%es(0:i_noe-1), stat=istat)

    IF (istat >0) THEN
      WRITE (message_text,'(a)')    'Problem in allocating spheres_on_levels%es, base_level'
      CALL message ('',TRIM(message_text))
      CALL finish ('init_topology', message_text)
    ENDIF

    ALLOCATE(spheres_on_levels(base_level)%vs(0:i_nov-1), stat=istat)

    IF (istat >0) THEN
      WRITE (message_text,'(a)')    'Problem in allocating spheres_on_levels%vs, base_level'
      CALL message ('',TRIM(message_text))
      CALL finish ('init_topology', message_text)
    ENDIF

    WRITE (message_text,'(a,i3,a,i9)')                                         &
      & 'Number of triangular (primal) cells on grid level ', base_level,    &
      & ', is ', spheres_on_levels(base_level)%no_triangles
    CALL message ('',TRIM(message_text))

    WRITE (message_text,'(a,i3,a,i9)')                                         &
      & 'Number of hexagonal  (dual)   cells on grid level ', base_level,    &
      & ', is ', spheres_on_levels(base_level)%no_vertices
    CALL message ('',TRIM(message_text))

    !   !--------------------------------------------------------------------------
    !   ! setup usable grid levels - allocate those arrays
    !   ! these allocated triangles get connected in build_tree

    DO ji = 0, number_of_levels

      IF (.not.lplane) THEN
        i_not = (20*kroot*kroot)*4**ji
        i_nov = (10*kroot*kroot)*4**ji + 2
        i_noe = (30*kroot*kroot)*4**ji
      ELSE
        i_not = ( 8*kroot*kroot)*4**ji
        i_nov = ( 4*kroot*kroot)*4**ji
        i_noe = (12*kroot*kroot)*4**ji
      ENDIF

      spheres_on_levels(ji)%level = ji

      ! required for the recursive generation of topology information

      spheres_on_levels(ji)%idx         = 0

      spheres_on_levels(ji)%no_triangles  = i_not
      spheres_on_levels(ji)%no_vertices   = i_nov
      spheres_on_levels(ji)%no_edges      = i_noe

      ALLOCATE(spheres_on_levels(ji)%ts(0:i_not-1), stat=istat)

      IF (istat >0) THEN
        WRITE (message_text,'(a,i2)')    'Problem in allocating spheres_on_levels%ts, level ',ji
        CALL message ('',TRIM(message_text))
        CALL finish ('init_topology', message_text)
      ENDIF

      ALLOCATE(spheres_on_levels(ji)%es(0:i_noe-1), stat=istat)

      IF (istat >0) THEN
        WRITE (message_text,'(a,i2)')    'Problem in allocating spheres_on_levels%es, level ',ji
        CALL message ('',TRIM(message_text))
        CALL finish ('init_topology', message_text)
      ENDIF

      ALLOCATE(spheres_on_levels(ji)%vs(0:i_nov-1), stat=istat)

      IF (istat >0) THEN
        WRITE (message_text,'(a,i2)')    'Problem in allocating spheres_on_levels%ts, level ',ji
        CALL message ('',TRIM(message_text))
        CALL finish ('init_topology', message_text)
      ENDIF

      WRITE (message_text,'(a,i3,a,i9)')                                     &
        & 'Number of triangular (primal) cells on grid level ', ji,       &
        & ', is ', spheres_on_levels(ji)%no_triangles
      CALL message ('',TRIM(message_text))

      WRITE (message_text,'(a,i3,a,i9)')                                     &
        & 'Number of hexagonal  (dual)   cells on grid level ', ji,      &
        & ', is ', spheres_on_levels(ji)%no_vertices
      CALL message ('',TRIM(message_text))

    END DO

  END SUBROUTINE init_topology_grid

  !-------------------------------------------------------------------------

  !

  !>
  !! Deallocates everything belonging to topology.
  !!
  !!
  !! @par Revision History
  !!  Thomas Heinze, DWD, Offenbach, 2006-10-25
  !!  Modification by Thomas Heinze, DWD, Offenbach, (2007-05-07):
  !!  - renamed destruct_topology to destruct_topology_graph
  !!  - adapted structure to graph generator
  !!
  SUBROUTINE destruct_topology_graph (k_levels)
    !

    ! number of levels to be generated below the
    ! root grid (level 0)

    INTEGER, INTENT(in) :: k_levels

    INTEGER :: number_of_levels
    INTEGER :: ji, istat

    !-------------------------------------------------

    ! now set the final numbers of levels and kroot
    DEALLOCATE(dummy_c, stat=istat)
    IF (istat >0) THEN
      WRITE (message_text,'(a)') 'Problem in deallocating dummy_c'
      CALL finish ('destruct_topology_graph', TRIM(message_text))
    ENDIF

    DEALLOCATE(dummy_e, stat=istat)
    IF (istat >0) THEN
      WRITE (message_text,'(a)') 'Problem in deallocating dummy_e'
      CALL finish ('destruct_topology_graph', TRIM(message_text))
    ENDIF

    DEALLOCATE(dummy_v, stat=istat)
    IF (istat >0) THEN
      WRITE (message_text,'(a)') 'Problem in deallocating dummy_v'
      CALL finish ('destruct_topology_graph', TRIM(message_text))
    ENDIF

    number_of_levels = k_levels

    ! nag fial to deallocate the top level
    ! Temporarly disable it
    DO ji = number_of_levels-1, base_level, -1
      istat=0
      IF (ASSOCIATED(spheres_on_levels(ji)%ts)) &
        & DEALLOCATE(spheres_on_levels(ji)%ts, stat=istat)
      IF (istat >0) THEN
        WRITE (message_text,'(a,i2)')    'Problem in deallocating spheres_on_levels%ts, level ',ji
        CALL message ('',TRIM(message_text))
        CALL finish ('destruct_topology', message_text)
      ENDIF

      IF (ASSOCIATED(spheres_on_levels(ji)%es)) &
         & DEALLOCATE(spheres_on_levels(ji)%es, stat=istat)
      IF (istat >0) THEN
        WRITE (message_text,'(a,i2)')    'Problem in deallocating spheres_on_levels%es, level ',ji
        CALL message ('',TRIM(message_text))
        CALL finish ('destruct_topology', message_text)
      ENDIF

      IF (ASSOCIATED(spheres_on_levels(ji)%es)) &
        & DEALLOCATE(spheres_on_levels(ji)%es, stat=istat)
      IF (istat >0) THEN
        WRITE (message_text,'(a,i2)')    'Problem in deallocating spheres_on_levels%ts, level ',ji
        CALL message ('',TRIM(message_text))
        CALL finish ('destruct_topology', message_text)
      ENDIF

    END DO

    ! root of 'tree', starting on level 0 - not on the original icosahedron
    !
    ! it does the allocation of the root matrix, containing the topology of the
    ! root grid; triangles on the root grid are identified by a local index
    ! (first index of root) and by index of the corresponding main triangle
    ! belonging to the icosahedron

    DEALLOCATE (root, stat=istat)

    IF (istat >0) THEN
      WRITE (message_text,'(a)')    'Problem in deallocating root'
      CALL message ('',TRIM(message_text))
      CALL finish ('destruct_topology', message_text)
    ENDIF

    ! Container for all spherical grids on all levels

    DEALLOCATE (spheres_on_levels, stat=istat)

    IF (istat >0) THEN
      WRITE (message_text,'(a)')    'Problem in deallocating spheres_on_levels'
      CALL message ('',TRIM(message_text))
      CALL finish ('destruct_topology', message_text)
    ENDIF

    ! The icosahedron itself

    DEALLOCATE (icosahedron, stat=istat)

    IF (istat >0) THEN
      WRITE (message_text,'(a)')    'Problem in dedeallocating icosahedron '
      CALL message ('',TRIM(message_text))
      CALL finish ('destruct_topology', message_text)
    ENDIF

  END SUBROUTINE destruct_topology_graph

  !-------------------------------------------------------------------------
  !
  !

  !>
  !! Deallocates everything belonging to topology.
  !!
  !!
  !! @par Revision History
  !!  Thomas Heinze, DWD, Offenbach, 2006-10-25
  !!  Modification by Thomas Heinze, DWD, Offenbach, (2007-05-07):
  !!  - renamed destruct_topology to destruct_topology_graph
  !!  - adapted structure to graph generator
  !!
  SUBROUTINE destruct_topology_grid (k_levels)
    !

    ! number of levels to be generated below the
    ! root grid (level 0)

    INTEGER, INTENT(in) :: k_levels

    INTEGER :: number_of_levels
    INTEGER :: ji, istat

    !-------------------------------------------------

    DEALLOCATE(dummy_c, stat=istat)
    IF (istat >0) THEN
      WRITE (message_text,'(a)') 'Problem in deallocating dummy_c'
      CALL finish ('destruct_topology_grid', TRIM(message_text))
    ENDIF

    DEALLOCATE(dummy_e, stat=istat)
    IF (istat >0) THEN
      WRITE (message_text,'(a)') 'Problem in deallocating dummy_e'
      CALL finish ('destruct_topology_grid', TRIM(message_text))
    ENDIF

    DEALLOCATE(dummy_v, stat=istat)
    IF (istat >0) THEN
      WRITE (message_text,'(a)') 'Problem in deallocating dummy_v'
      CALL finish ('destruct_topology_grid', TRIM(message_text))
    ENDIF

    ! now deallocate spheres_on_levels structures

    number_of_levels = k_levels

    ! nag fails to deallocate the top level, probably deallocated elsewhere
    ! temporarly disable to top level deallocation
    DO ji = number_of_levels-1, base_level, -1
      istat=0
      IF (ASSOCIATED(spheres_on_levels(ji)%ts)) &
        & DEALLOCATE(spheres_on_levels(ji)%ts, stat=istat)

      IF (istat >0) THEN
        WRITE (message_text,'(a,i2)')    'Problem in deallocating spheres_on_levels%ts, level ',ji
        CALL message ('',TRIM(message_text))
        CALL finish ('destruct_topology', message_text)
      ENDIF

      IF (ASSOCIATED(spheres_on_levels(ji)%es)) &
        & DEALLOCATE(spheres_on_levels(ji)%es, stat=istat)

      IF (istat >0) THEN
        WRITE (message_text,'(a,i2)')    'Problem in deallocating spheres_on_levels%es, level ',ji
        CALL message ('',TRIM(message_text))
        CALL finish ('destruct_topology', message_text)
      ENDIF

        IF (ASSOCIATED(spheres_on_levels(ji)%vs)) &
          & DEALLOCATE(spheres_on_levels(ji)%vs, stat=istat)

        IF (istat >0) THEN
          WRITE (message_text,'(a,i2)')   &
             'Problem in deallocating spheres_on_levels%ts, level ',ji
          CALL message ('',TRIM(message_text))
          CALL finish ('destruct_topology', message_text)
        ENDIF

    END DO

    ! root of 'tree', starting on level 0 - not on the original icosahedron
    !
    ! it does the allocation of the root matrix, containing the topology of the
    ! root grid; triangles on the root grid are identified by a local index
    ! (first index of root) and by index of the corresponding main triangle
    ! belonging to the icosahedron

    DEALLOCATE (root, stat=istat)

    IF (istat >0) THEN
      WRITE (message_text,'(a)')    'Problem in deallocating root'
      CALL message ('',TRIM(message_text))
      CALL finish ('destruct_topology', message_text)
    ENDIF

    ! Container for all spherical grids on all levels

    DEALLOCATE (spheres_on_levels, stat=istat)

    IF (istat >0) THEN
      WRITE (message_text,'(a)')    'Problem in deallocating spheres_on_levels'
      CALL message ('',TRIM(message_text))
      CALL finish ('destruct_topology', message_text)
    ENDIF

    ! The icosahedron itself

    DEALLOCATE (icosahedron, stat=istat)

    IF (istat >0) THEN
      WRITE (message_text,'(a)')    'Problem in deallocating icosahedron '
      CALL message ('',TRIM(message_text))
      CALL finish ('destruct_topology', message_text)
    ENDIF

  END SUBROUTINE destruct_topology_grid

  !-------------------------------------------------------------------------
  !
  !

  !>
  !!  This subroutine is a driver for the call to the recursive.
  !!
  !!  This subroutine is a driver for the call to the recursive
  !!  function <i>build_tree, </i> which sets up the tree structure
  !!  underlying each triangle of the original icosahedron.
  !!
  !! @par Revision History
  !!  Luis Kornblueh, MPI-M, Hamburg (2004-03)
  !!  Modification by Almut Gassmann, MPI-M (2009-01-09)
  !!  - adapt for planar graph generation
  !!
  SUBROUTINE generate_tree (lplane)

    LOGICAL, INTENT(in) :: lplane
    INTEGER :: ji, j, n_faces_m1
    !-------------------------------------------------
    ! now fill data structures

    IF(.not. lplane) THEN
      n_faces_m1 = 19
    ELSE
      n_faces_m1 = 7
    ENDIF

    DO j = 0, n_faces_m1
      DO ji = 0, kroot*kroot-1
        root(ji,j)%t => build_tree (0)
      END DO
    END DO

  END SUBROUTINE generate_tree
  !EOC
  !-------------------------------------------------------------------------
  !>
  !!   Recursive function that sets up the tree structure.
  !!
  !!   Recursive function that sets up the tree structure
  !!   underlying each triangle of the original icosahedron.
  !!   The pointers to parent and offspring triangles are initialized.
  !!
  !! @par Revision History
  !!  Luis Kornblueh, MPI-M, Hamburg, (2004-03)
  !!
  RECURSIVE FUNCTION build_tree (k_lvl) result(current)
    TYPE(t_triangle), POINTER :: current

    INTEGER, INTENT(in) :: k_lvl

    INTEGER :: i_cidx
    !-------------------------------------------------
    !BOC
    i_cidx = spheres_on_levels(k_lvl)%idx
    current => spheres_on_levels(k_lvl)%ts(i_cidx)

    current%info%level = k_lvl

    IF (k_lvl == number_of_levels) THEN
      current%sub_triangle0 => NULL()
      current%sub_triangle1 => NULL()
      current%sub_triangle2 => NULL()
      current%sub_triangle3 => NULL()
    ELSE
      current%sub_triangle0 => build_tree(k_lvl+1)
      current%sub_triangle1 => build_tree(k_lvl+1)
      current%sub_triangle2 => build_tree(k_lvl+1)
      current%sub_triangle3 => build_tree(k_lvl+1)

      ! what about parents? - set them here for convenience

      current%sub_triangle0%parent => current
      current%sub_triangle1%parent => current
      current%sub_triangle2%parent => current
      current%sub_triangle3%parent => current
    END IF

    spheres_on_levels(k_lvl)%idx = i_cidx+1

  END FUNCTION build_tree

  !-------------------------------------------------------------------------
  !>
  !!  Generates the graph that is used internally by the grid generator.
  !!
  !!  Generates the graph that is used internally by the grid generator
  !!  to represent the neighborhood relationships among cells
  !!  on the icosahedral grids. First, the neighborhood
  !!  relationships for the  original icosahedron are
  !!  initialized by <i>init_icosahedron.</i> Then,  neighborhood
  !!  relationships for the root grid are  initialized by <i>init_root.</i>
  !!  Finally,  neighborhood  relationships for the finer grid levels  are  initialized
  !!  recursively by <i>build_cell_graph.</i>
  !!
  !! @par Revision History
  !!  Luis Kornblueh, MPI-M, Hamburg, (2004-03)
  !!  Modification by Almut Gassmann, MPI-M (2009-01-09)
  !!
  SUBROUTINE generate_graph(k_lvl,lplane)
    INTEGER, INTENT(in) :: k_lvl
    LOGICAL, INTENT(in) :: lplane
    INTEGER :: ji, j, n_faces_m1
    !-------------------------------------------------
    !BOC
    CALL init_icosahedron(lplane)
    CALL init_root(lplane)

    ! now fill data structures:
    !
    ! first build graph of
    ! neighbourhood relationships among triangular cells
    ! (part of it)
    IF(.not. lplane) THEN
      n_faces_m1 = 19
    ELSE
      n_faces_m1 = 7
    ENDIF

    DO j = 0, n_faces_m1
      DO ji = 0, kroot*kroot-1
        CALL build_cell_graph (root(ji,j)%t)
      END DO
    END DO

    !
    ! then complete cell graph and build graphs for edges and vertices
    !
    DO ji = 0,k_lvl
      CALL build_full_graph(spheres_on_levels(ji))
    ENDDO

  END SUBROUTINE generate_graph

  !-------------------------------------------------------------------------
  !
  !
  !>
  !!  Generates the graph that is used internally by the grid generator.
  !!
  !!  Generates the graph that is used internally by the grid generator
  !!  to represent the neighborhood relationships among cells
  !!  on the icosahedral grids. First, the neighborhood
  !!  relationships for the  original icosahedron are
  !!  initialized by <i>init_icosahedron.</i> Then,  neighborhood
  !!  relationships for the root grid are  initialized by <i>init_root.</i>
  !!  Finally,  neighborhood  relationships for the finer grid levels
  !!  are read from files created earlier by graph generator.
  !!
  !! @par Revision History
  !!  Original version by Thomas Heinze, DWD, (2007-05-03) derived from
  !!  generate_graph.
  !!  Modification by Almut Gassmann, MPI-M (2009-01-09)
  !!  - adaption for planar grid
  !!
  SUBROUTINE read_graph(k_lvl,lplane)
    !
    INTEGER, INTENT(in) :: k_lvl     ! finest level
    LOGICAL, INTENT(in) :: lplane

    INTEGER :: ji, j     ! loop indices
    INTEGER :: n_faces_m1
    !-------------------------------------------------

    CALL init_icosahedron (lplane)
    CALL init_root(lplane)

    ! now fill data structures:
    !
    ! first build graph of
    ! neighbourhood relationships among triangular cells
    ! (part of it)
    !
    IF(.not. lplane) THEN
      n_faces_m1 = 19
    ELSE
      n_faces_m1 = 7
    ENDIF

    DO j = 0, n_faces_m1
      DO ji = 0, kroot*kroot-1
        CALL build_cell_graph (root(ji,j)%t)
      END DO
    END DO

    !
    ! then complete cell graph and build graphs for edges and vertices
    !

    DO ji = 0, k_lvl
      CALL input_graph(kroot, spheres_on_levels(ji),lplane)
    ENDDO

  END SUBROUTINE read_graph

  !-------------------------------------------------------------------------
  !
  !
  !>
  !!  Generates the graph that is used internally by the grid generator.
  !!
  !!  Generates the graph that is used internally by the grid generator
  !!  to represent the neighborhood relationships among cells
  !!  on the icosahedral grids. First, the neighborhood
  !!  relationships for the  original icosahedron are
  !!  initialized by <i>init_icosahedron.</i> Then,  neighborhood
  !!  relationships for the root grid are  initialized by <i>init_root.</i>
  !!  Finally,  neighborhood  relationships for the finer grid levels
  !!  are read from files created earlier by graph generator.
  !!
  !! @par Revision History
  !!  Original version by Thomas Heinze, DWD, (2007-05-03)
  !!
  SUBROUTINE write_graph(k_root, k_lvl, lplane)
    !
    INTEGER, INTENT(in) :: k_root, k_lvl     ! root and finest level
    LOGICAL, INTENT(in) :: lplane

    INTEGER :: j         ! loop index

    !-------------------------------------------------

    ! write full graphs to GRAPHMAP files

    DO j = 0, k_lvl
      CALL output_graph(k_root, spheres_on_levels(j),lplane)
    ENDDO

  END SUBROUTINE write_graph

  !-------------------------------------------------------------------------

  !>
  !!  Defines the addressing, the orientation and the neighborhood relationships.
  !!
  !!  Defines the addressing, the orientation and the neighborhood relationships
  !!  for triangles on the original icosahedron. The original icosahedron
  !!  has 20 triangular faces, 12 vertices (which are faces of the dual grid)
  !!  30 edges. The triangles are characterised by orientations <i>up</i> if one
  !!  of their vertices is pointing towards the North Pole or <i>down</i> if one
  !!  of their vertices is pointing towards the South Pole.
  !!  The 20 faces are numbered as follows:
  !! \\vskip 1.cm
  !! @f{tabular}{
  !! [c]{|c|c|c|c|c|c|}  \hline
  !! 0    & 1  & 2 &  3 & 4 & up, vertex at NP \\\
  !! 5    & 6  & 7  & 8 & 9 &  down, across Equator\\\
  !! 10    & 11  & 12  &  13 & 14 & up, across Equator \\\
  !! 15    & 16  & 17  &  18 & 19 & down, vertex at SP \\\
  !! \hline
  !! @f}
  !!
  !! @par Revision History
  !!  Luis Kornblueh, MPI-M, Hamburg, (2004-03)
  !!  Modification by Almut Gassmann, MPI-M (2009-01-09)
  !!
  SUBROUTINE init_icosahedron(lplane)

    LOGICAL, INTENT(in) :: lplane
    INTEGER :: j, n_faces_m1
    !-------------------------------------------------

    ! The number of faces minus 1 of the icosahedron
    IF (.not. lplane) THEN
      n_faces_m1 = 19  ! icosahedron
    ELSE
      n_faces_m1 = 7   ! pseudo icosahedron on a plane
    ENDIF

    DO j = 0, n_faces_m1
      icosahedron(j)%t => spheres_on_levels(base_level)%ts(j)
    ENDDO

    DO j = 0, n_faces_m1
      ! set just for identification the triangle number and configuration
      icosahedron(j)%t%info%triangle_number = j
      icosahedron(j)%t%info%configuration   = j
      ! nullify parent pointer if we are on the base level
      icosahedron(j)%t%parent => NULL()
    END DO

    ! have to set neighbors on this level by hand

    IF (.not. lplane) THEN

      icosahedron(0)%t%info%orientation = IOR(IBSET(0,pole),IBSET(0,up))
      icosahedron(0)%t%neighbor0 => icosahedron(4)%t
      icosahedron(0)%t%neighbor1 => icosahedron(5)%t
      icosahedron(0)%t%neighbor2 => icosahedron(1)%t

      icosahedron(1)%t%info%orientation = IOR(IBSET(0,pole),IBSET(0,up))
      icosahedron(1)%t%neighbor0 => icosahedron(0)%t
      icosahedron(1)%t%neighbor1 => icosahedron(6)%t
      icosahedron(1)%t%neighbor2 => icosahedron(2)%t

      icosahedron(2)%t%info%orientation = IOR(IBSET(0,pole),IBSET(0,up))
      icosahedron(2)%t%neighbor0 => icosahedron(1)%t
      icosahedron(2)%t%neighbor1 => icosahedron(7)%t
      icosahedron(2)%t%neighbor2 => icosahedron(3)%t

      icosahedron(3)%t%info%orientation = IOR(IBSET(0,pole),IBSET(0,up))
      icosahedron(3)%t%neighbor0 => icosahedron(2)%t
      icosahedron(3)%t%neighbor1 => icosahedron(8)%t
      icosahedron(3)%t%neighbor2 => icosahedron(4)%t

      icosahedron(4)%t%info%orientation = IOR(IBSET(0,pole),IBSET(0,up))
      icosahedron(4)%t%neighbor0 => icosahedron(3)%t
      icosahedron(4)%t%neighbor1 => icosahedron(9)%t
      icosahedron(4)%t%neighbor2 => icosahedron(0)%t

      icosahedron(5)%t%info%orientation = IOR(IBSET(0,main),IBSET(0,down))
      icosahedron(5)%t%neighbor2 => icosahedron(14)%t
      icosahedron(5)%t%neighbor1 => icosahedron(0)%t
      icosahedron(5)%t%neighbor0 => icosahedron(10)%t

      icosahedron(6)%t%info%orientation = IOR(IBSET(0,main),IBSET(0,down))
      icosahedron(6)%t%neighbor2 => icosahedron(10)%t
      icosahedron(6)%t%neighbor1 => icosahedron(1)%t
      icosahedron(6)%t%neighbor0 => icosahedron(11)%t

      icosahedron(7)%t%info%orientation = IOR(IBSET(0,main),IBSET(0,down))
      icosahedron(7)%t%neighbor2 => icosahedron(11)%t
      icosahedron(7)%t%neighbor1 => icosahedron(2)%t
      icosahedron(7)%t%neighbor0 => icosahedron(12)%t

      icosahedron(8)%t%info%orientation = IOR(IBSET(0,main),IBSET(0,down))
      icosahedron(8)%t%neighbor2 => icosahedron(12)%t
      icosahedron(8)%t%neighbor1 => icosahedron(3)%t
      icosahedron(8)%t%neighbor0 => icosahedron(13)%t

      icosahedron(9)%t%info%orientation = IOR(IBSET(0,main),IBSET(0,down))
      icosahedron(9)%t%neighbor2 => icosahedron(13)%t
      icosahedron(9)%t%neighbor1 => icosahedron(4)%t
      icosahedron(9)%t%neighbor0 => icosahedron(14)%t

      icosahedron(10)%t%info%orientation = IOR(IBSET(0,main),IBSET(0,up))
      icosahedron(10)%t%neighbor0 => icosahedron(5)%t
      icosahedron(10)%t%neighbor1 => icosahedron(15)%t
      icosahedron(10)%t%neighbor2 => icosahedron(6)%t

      icosahedron(11)%t%info%orientation = IOR(IBSET(0,main),IBSET(0,up))
      icosahedron(11)%t%neighbor0 => icosahedron(6)%t
      icosahedron(11)%t%neighbor1 => icosahedron(16)%t
      icosahedron(11)%t%neighbor2 => icosahedron(7)%t

      icosahedron(12)%t%info%orientation = IOR(IBSET(0,main),IBSET(0,up))
      icosahedron(12)%t%neighbor0 => icosahedron(7)%t
      icosahedron(12)%t%neighbor1 => icosahedron(17)%t
      icosahedron(12)%t%neighbor2 => icosahedron(8)%t

      icosahedron(13)%t%info%orientation = IOR(IBSET(0,main),IBSET(0,up))
      icosahedron(13)%t%neighbor0 => icosahedron(8)%t
      icosahedron(13)%t%neighbor1 => icosahedron(18)%t
      icosahedron(13)%t%neighbor2 => icosahedron(9)%t

      icosahedron(14)%t%info%orientation = IOR(IBSET(0,main),IBSET(0,up))
      icosahedron(14)%t%neighbor0 => icosahedron(9)%t
      icosahedron(14)%t%neighbor1 => icosahedron(19)%t
      icosahedron(14)%t%neighbor2 => icosahedron(5)%t

      icosahedron(15)%t%info%orientation = IOR(IBSET(0,pole),IBSET(0,down))
      icosahedron(15)%t%neighbor2 => icosahedron(19)%t
      icosahedron(15)%t%neighbor1 => icosahedron(10)%t
      icosahedron(15)%t%neighbor0 => icosahedron(16)%t

      icosahedron(16)%t%info%orientation = IOR(IBSET(0,pole),IBSET(0,down))
      icosahedron(16)%t%neighbor2 => icosahedron(15)%t
      icosahedron(16)%t%neighbor1 => icosahedron(11)%t
      icosahedron(16)%t%neighbor0 => icosahedron(17)%t

      icosahedron(17)%t%info%orientation = IOR(IBSET(0,pole),IBSET(0,down))
      icosahedron(17)%t%neighbor2 => icosahedron(16)%t
      icosahedron(17)%t%neighbor1 => icosahedron(12)%t
      icosahedron(17)%t%neighbor0 => icosahedron(18)%t

      icosahedron(18)%t%info%orientation = IOR(IBSET(0,pole),IBSET(0,down))
      icosahedron(18)%t%neighbor2 => icosahedron(17)%t
      icosahedron(18)%t%neighbor1 => icosahedron(13)%t
      icosahedron(18)%t%neighbor0 => icosahedron(19)%t

      icosahedron(19)%t%info%orientation = IOR(IBSET(0,pole),IBSET(0,down))
      icosahedron(19)%t%neighbor2 => icosahedron(18)%t
      icosahedron(19)%t%neighbor1 => icosahedron(14)%t
      icosahedron(19)%t%neighbor0 => icosahedron(15)%t

    ELSE ! lplane ==.TRUE.

      icosahedron(0)%t%info%orientation = IBSET(0,down)
      icosahedron(0)%t%neighbor0 => icosahedron(3)%t
      icosahedron(0)%t%neighbor1 => icosahedron(7)%t
      icosahedron(0)%t%neighbor2 => icosahedron(2)%t

      icosahedron(1)%t%info%orientation = IBSET(0,down)
      icosahedron(1)%t%neighbor0 => icosahedron(2)%t
      icosahedron(1)%t%neighbor1 => icosahedron(6)%t
      icosahedron(1)%t%neighbor2 => icosahedron(3)%t

      icosahedron(2)%t%info%orientation = IBSET(0,up)
      icosahedron(2)%t%neighbor0 => icosahedron(1)%t
      icosahedron(2)%t%neighbor1 => icosahedron(4)%t
      icosahedron(2)%t%neighbor2 => icosahedron(0)%t

      icosahedron(3)%t%info%orientation = IBSET(0,up)
      icosahedron(3)%t%neighbor0 => icosahedron(0)%t
      icosahedron(3)%t%neighbor1 => icosahedron(5)%t
      icosahedron(3)%t%neighbor2 => icosahedron(1)%t

      icosahedron(4)%t%info%orientation = IBSET(0,down)
      icosahedron(4)%t%neighbor0 => icosahedron(7)%t
      icosahedron(4)%t%neighbor1 => icosahedron(2)%t
      icosahedron(4)%t%neighbor2 => icosahedron(6)%t

      icosahedron(5)%t%info%orientation = IBSET(0,down)
      icosahedron(5)%t%neighbor0 => icosahedron(6)%t
      icosahedron(5)%t%neighbor1 => icosahedron(3)%t
      icosahedron(5)%t%neighbor2 => icosahedron(7)%t

      icosahedron(6)%t%info%orientation = IBSET(0,up)
      icosahedron(6)%t%neighbor0 => icosahedron(5)%t
      icosahedron(6)%t%neighbor1 => icosahedron(1)%t
      icosahedron(6)%t%neighbor2 => icosahedron(4)%t

      icosahedron(7)%t%info%orientation = IBSET(0,up)
      icosahedron(7)%t%neighbor0 => icosahedron(4)%t
      icosahedron(7)%t%neighbor1 => icosahedron(0)%t
      icosahedron(7)%t%neighbor2 => icosahedron(5)%t

    ENDIF

  END SUBROUTINE init_icosahedron
  !-------------------------------------------------------------------------
  !>
  !!  This is the driver routine for the call to.
  !!
  !!  This is the driver routine for the call to
  !!   <i>init_root_graph_neighbors</i>, which defines the addressing,
  !!   the orientation and the neighborhood relationships
  !!   for triangles on the root level grid.
  !!
  !! @par Revision History
  !!  Luis Kornblueh, MPI-M, Hamburg, (2004-03)
  !!  Modification by Almut Gassmann, MPI-M (2009-01-09)
  !!  - adaption for planar grid
  !!
  SUBROUTINE init_root (lplane)

    LOGICAL, INTENT(in) :: lplane
    TYPE(t_triangle), POINTER :: current
    INTEGER :: j, n_faces_m1
    !-------------------------------------------------

    ! The number of faces minus 1 of the icosahedron
    IF (.not. lplane) THEN
      n_faces_m1 = 19  ! icosahedron
    ELSE
      n_faces_m1 = 7   ! pseudo icosahedron on a plane
    ENDIF

    DO j = 0, n_faces_m1

      current => icosahedron(j)%t
      CALL init_root_graph_neighbors(current)

    END DO

  END SUBROUTINE init_root
  !-------------------------------------------------------------------------

  !>
  !!   This subroutine defines the addressing,.
  !!
  !!   This subroutine defines the addressing,
  !!   the orientation and the neighborhood relationships
  !!   for triangles on the root level grid, which is obtained
  !!   from the original icosahedron by subdivision of the edges
  !!   of the 20 initial triangles in <i>kroot -1 </i> parts.
  !!   For example, with <i>kroot=1,</i> the original triangles are used as root grid.
  !!   With <i>kroot=2,</i> the root grid coincides with the first
  !!   possible dyadic refinement of the original icosahedron; one original triangle
  !!   is subdivided into 4 triangles by halving the triangle edges.
  !!   With <i>kroot=3,</i> the root grid edges are trisected; one original triangle
  !!   is subdivided into 9 triangles by dividing the edge into three parts.
  !!   In principle any value of  <i>kroot </i> is admissible, but reccommended
  !!   for practical purposes are 2, 3 or 5. With   <i>kroot=2,</i> the most usual
  !!   icosahedral grid hierarchy is recovered.
  !!
  !! @par Revision History
  !!  Luis Kornblueh, MPI-M, Hamburg, (2004-03)
  !!
  SUBROUTINE init_root_graph_neighbors (root_triangle)


    ! triangle for which neighbors have
    ! to be given (number, orientation)
    TYPE(t_triangle), POINTER :: root_triangle    ! INTENT(in)

    ! UP:            DOWN:  base
    !         /\            _____
    ! left   /  \   right   \   /  left
    !       ------           \ /
    !        base
    !
    ! downward oriented triangles are rotated so that the same
    ! numbering scheme can be used!
    INTEGER :: ji, ii, j, ict, ils, ile, iils, iile
    INTEGER :: icase, iorientation
    INTEGER :: icn, iln, ibn, irn
    INTEGER :: icrn, icrp, ictb, ikbn
    LOGICAL :: leven_line, leven, llastrow

    TYPE(t_triangle), POINTER :: current


    !-------------------------------------------------
    !BOC

    current => root_triangle

    icn = current%info%triangle_number;

    iln = current%neighbor0%info%triangle_number;
    ibn = current%neighbor1%info%triangle_number;
    irn = current%neighbor2%info%triangle_number;

    icase = current%info%orientation

    IF (BTEST(icase,up)) THEN
      iorientation = IBSET(0,up)
    ELSE
      iorientation = IBSET(0,down)
    END IF

    ! setup the root structure

    ! loop over sub-divisions of original icosahedron triangle faces (rows)

    DO ji = 0, kroot-1

      ! This are the indices for the upward looking triangles
      ! with respect to north.

      ils = ji*ji           ! number of first triangle in row
      ile = ji*ji+2*ji       ! number of last triangle in row

      ! This are the indices for the downward looking triangles
      ! with respect to north.
      ! Indexing is rotated for all triangles, so the indexing for the
      ! downward looking trianlges must only be inverted (top-down).

      ii   = kroot-1-ji    ! invert direction for interior triangle neighbors
      iils = ii*ii        ! number of last triangle in row
      iile = ii*ii+2*ii   ! number of first triangle in row

      icrn = (ji+1)*(ji+1)  ! number of triangle in the row below the current one
      icrp = (ji-1)*(ji-1)  ! number of triangle in the row above the current one

      ! line type

      IF (MOD(ji,2) == 0) THEN
        leven_line = .true.
      ELSE
        leven_line = .false.
      END IF

      llastrow = .false.
      IF (ji == kroot-1) llastrow = .true.

      DO j = 0, 2*ji      ! loop over triangles in each row

        IF (MOD(j,2) == 0) THEN
          leven = .true.
        ELSE
          leven = .false.
        END IF

        ! index number of current triangle

        ict = ji*ji+j

        root(ict,icn)%t%info%triangle_number = ict

        root(ict,icn)%t%info%configuration = 100*icn+ict

        ! connect level 0 triangle to its icosahedron parent

        root(ict,icn)%t%info%parent_number = icn
        root(ict,icn)%t%parent => current

        IF (BTEST(icase,pole) .and. BTEST(icase,up)) THEN

          ! first case: pole and upward looking, northern hemispheric pole

          ! select indices for base neighbor inside a triangle

          ikbn = icn
          IF (leven_line) THEN
            IF (MOD(ict,2) == 0) THEN
              ictb = icrn+1+j
              IF (llastrow) THEN
                ikbn  = ibn
                ictb = ile-j
              END IF
            ELSE
              ictb = icrp-1+j
            END IF
          ELSE
            IF (MOD(ict,2) == 0) THEN
              ictb = icrp-1+j
            ELSE
              ictb = icrn+1+j
              IF (llastrow) THEN
                ikbn  = ibn
                ictb = ile-j
              END IF
            END IF
          END IF

          IF (ict == ils .and. ict == ile) THEN  ! topmost triangle

            root(ict,icn)%t%neighbor0 => root(ils,iln)%t
            root(ict,icn)%t%neighbor1 => root(ictb,ikbn)%t
            root(ict,icn)%t%neighbor2 => root(ile,irn)%t

            root(ict,icn)%t%info%orientation = IOR(IBSET(0,pole),IBSET(0,up))

          ELSE IF (ict == ils) THEN              ! triangle on the left bound

            root(ict,icn)%t%neighbor0 => root(ile,iln)%t
            root(ict,icn)%t%neighbor2 => root(ict+1,icn)%t
            root(ict,icn)%t%neighbor1 => root(ictb,ikbn)%t

            root(ict,icn)%t%info%orientation = &
              & IOR(IOR(IBSET(0,border),IBSET(0,up)),IBSET(0,left))

          ELSE IF (ict == ile) THEN              ! triangle on the right bound

            root(ict,icn)%t%neighbor0 => root(ict-1,icn)%t
            root(ict,icn)%t%neighbor2 => root(ils,irn)%t
            root(ict,icn)%t%neighbor1 => root(ictb,ikbn)%t

            root(ict,icn)%t%info%orientation = &
              & IOR(IOR(IBSET(0,border),IBSET(0,up)),IBSET(0,right))

          ELSE   ! (ils < ict .AND. ict < ile)   ! triangle in the interiour

            root(ict,icn)%t%neighbor1 => root(ictb,ikbn)%t

            IF (leven) THEN
              root(ict,icn)%t%neighbor0 => root(ict-1,icn)%t
              root(ict,icn)%t%neighbor2 => root(ict+1,icn)%t
              root(ict,icn)%t%info%orientation = IBSET(0,up)
            ELSE
              root(ict,icn)%t%neighbor0 => root(ict+1,icn)%t
              root(ict,icn)%t%neighbor2 => root(ict-1,icn)%t
              root(ict,icn)%t%info%orientation = IBSET(0,down)
            END IF

          END IF

        ELSE IF (BTEST(icase,pole) .and. BTEST(icase,down)) THEN

          ! second case: pole and downward looking, southern hemispheric pole

          ! select indices for base neighbor inside a triangle
          ikbn = icn
          IF (leven_line) THEN
            IF (MOD(ict,2) == 0) THEN
              ictb = icrn+1+j
              IF (llastrow) THEN
                ikbn  = ibn
                ictb = ile-j
              END IF
            ELSE
              ictb = icrp-1+j
            END IF
          ELSE
            IF (MOD(ict,2) == 0) THEN
              ictb = icrp-1+j
            ELSE
              ictb = icrn+1+j
              IF (llastrow) THEN
                ikbn  = ibn
                ictb = ile-j
              END IF
            END IF
          END IF

          IF (ict == ils .and. ict == ile) THEN  ! lowest triangle

            root(ict,icn)%t%neighbor2 => root(ils,irn)%t
            root(ict,icn)%t%neighbor1 => root(ictb,ikbn)%t
            root(ict,icn)%t%neighbor0 => root(ile,iln)%t

            root(ict,icn)%t%info%orientation = IOR(IBSET(0,pole),IBSET(0,down))

          ELSE IF (ict == ils) THEN              ! triangle on the left bound

            root(ict,icn)%t%neighbor2 => root(ict+1,icn)%t
            root(ict,icn)%t%neighbor1 => root(ictb,ikbn)%t
            root(ict,icn)%t%neighbor0 => root(ile,iln)%t

            root(ict,icn)%t%info%orientation = &
              & IOR(IOR(IBSET(0,border),IBSET(0,down)),IBSET(0,left))

          ELSE IF (ict == ile) THEN              ! triangle on the right bound

            root(ict,icn)%t%neighbor2 => root(ils,irn)%t
            root(ict,icn)%t%neighbor1 => root(ictb,ikbn)%t
            root(ict,icn)%t%neighbor0 => root(ict-1,icn)%t

            root(ict,icn)%t%info%orientation = &
              & IOR(IOR(IBSET(0,border),IBSET(0,down)),IBSET(0,right))

          ELSE   ! (ils < ict .AND. ict < ile)   ! triangle in the interior

            root(ict,icn)%t%neighbor1 => root(ictb,ikbn)%t

            IF (leven) THEN
              root(ict,icn)%t%neighbor0 => root(ict-1,icn)%t
              root(ict,icn)%t%neighbor2 => root(ict+1,icn)%t
              root(ict,icn)%t%info%orientation = IBSET(0,down)
            ELSE
              root(ict,icn)%t%neighbor0 => root(ict+1,icn)%t
              root(ict,icn)%t%neighbor2 => root(ict-1,icn)%t
              root(ict,icn)%t%info%orientation = IBSET(0,up)
            END IF

          END IF

        ELSE   ! interior of icosahedron

          ! third case

          ! select indices for base neighbor inside a triangle
          ikbn = icn
          IF (leven_line) THEN
            IF (MOD(ict,2) == 0) THEN
              ictb = icrn+1+j
              IF (llastrow) THEN
                ikbn  = ibn
                ictb = ile-j
              END IF
            ELSE
              ictb = icrp-1+j
            END IF
          ELSE
            IF (MOD(ict,2) == 0) THEN
              ictb = icrp-1+j
            ELSE
              ictb = icrn+1+j
              IF (llastrow) THEN
                ikbn  = ibn
                ictb = ile-j
              END IF
            END IF
          END IF

          IF (leven) THEN
            root(ict,icn)%t%info%orientation = iorientation
          ELSE
            root(ict,icn)%t%info%orientation = IEOR(up_down,iorientation)
          END IF

          IF (ict == ils .and. ict == ile) THEN  ! lowest/topmost triangle

            root(ict,icn)%t%neighbor0 => root(iils,iln)%t
            root(ict,icn)%t%neighbor1 => root(ictb,ikbn)%t
            root(ict,icn)%t%neighbor2 => root(iile,irn)%t

          ELSE IF (ict == ils) THEN              ! triangle on the left bound

            root(ict,icn)%t%neighbor0 => root(iils,iln)%t
            root(ict,icn)%t%neighbor1 => root(ictb,ikbn)%t
            root(ict,icn)%t%neighbor2 => root(ict+1,icn)%t

          ELSE IF (ict == ile) THEN              ! triangle on the right bound

            root(ict,icn)%t%neighbor0 => root(ict-1,icn)%t
            root(ict,icn)%t%neighbor1 => root(ictb,ikbn)%t
            root(ict,icn)%t%neighbor2 => root(iile,irn)%t

          ELSE   ! (ils < ict .AND. ict < ile)   ! triangle in the interior

            root(ict,icn)%t%neighbor1 => root(ictb,ikbn)%t

            IF (leven) THEN
              root(ict,icn)%t%neighbor0 => root(ict-1,icn)%t
              root(ict,icn)%t%neighbor2 => root(ict+1,icn)%t
            ELSE
              root(ict,icn)%t%neighbor0 => root(ict+1,icn)%t
              root(ict,icn)%t%neighbor2 => root(ict-1,icn)%t
            END IF

          END IF
        END IF
      ENDDO
    ENDDO

  END SUBROUTINE init_root_graph_neighbors
  !EOC
  !-------------------------------------------------------------------------
  !>
  !!  This recursive subroutine defines the addressing,.
  !!
  !!  This recursive subroutine defines the addressing,
  !!   the orientation and the neighborhood relationships
  !!   for triangles on the higher level grids, which are obtained
  !!   by bisection from the root grid (level 0).
  !!   The addressing is inherited from the tree structure
  !!   defined in the <i>build_tree </i> subroutine and should allow
  !!   to preserve data locality. The data used by the numerical model
  !!   inherit again the same ordering for the triangular grid cells
  !!   and orderings that are derived naturally from it for the edges
  !!   and vertices. In order to set up the  neighborhood relationships,
  !!   this subroutine uses the fact that
  !! <ul>
  !! <li> [1]        both neighbors on the left and right are of the same orientation,
  !!                 which are the pole triangles;
  !! <li> [2]        both neighbors on the left and right are of the opposite
  !!                 orientation, which are the triangles in the interior
  !! <li> [3]        one neighbor and the current have the same orientation and
  !!                 the neighbor on the other side is oriented opposite, which
  !!                 are the triangles on the edge off the original level -1 pole
  !!                 triangles.
  !! </ul>
  !!
  !! @par Revision History
  !!  Luis Kornblueh, MPI-M, Hamburg, (2004-03)
  !!
  RECURSIVE SUBROUTINE build_cell_graph (cn)


    TYPE(t_triangle), POINTER :: cn    ! INTENT(in)

    INTEGER :: iorientation
    !-------------------------------------------------
    !BOC
    IF (cn%info%level < number_of_levels) THEN

      ! Compute the right parameters for the next level of recursion:
      !
      ! - Number the sub triangles and set the orientations


      cn%sub_triangle0%info%triangle_number = 0
      cn%sub_triangle0%info%configuration   = 10*cn%info%configuration+0

      cn%sub_triangle0%info%parent_number   =  cn%info%idx


      cn%sub_triangle1%info%triangle_number = 1
      cn%sub_triangle1%info%configuration   = 10*cn%info%configuration+1

      cn%sub_triangle1%info%parent_number   =  cn%info%idx

      cn%sub_triangle2%info%triangle_number = 2
      cn%sub_triangle2%info%configuration   = 10*cn%info%configuration+2

      cn%sub_triangle2%info%parent_number   =  cn%info%idx

      cn%sub_triangle3%info%triangle_number = 3
      cn%sub_triangle3%info%configuration   = 10*cn%info%configuration+3

      cn%sub_triangle3%info%parent_number   =  cn%info%idx

      iorientation     = cn%info%orientation

      IF (BTEST(iorientation,pole)) THEN

        cn%sub_triangle0%info%orientation = iorientation
        cn%sub_triangle1%info%orientation = &
          & IBCLR(iorientation,pole)+IBSET(0,border)+IBSET(0,left)
        cn%sub_triangle2%info%orientation = &
          & IBCLR(IBCLR(IBCLR(IEOR(up_down,iorientation),pole),left),right)
        cn%sub_triangle3%info%orientation = &
          & IBCLR(iorientation,pole)+IBSET(0,border)+IBSET(0,right)

      ELSE IF (BTEST(iorientation,border)) THEN

        IF (BTEST(iorientation,left)) THEN
          cn%sub_triangle0%info%orientation = iorientation
          cn%sub_triangle1%info%orientation = iorientation
          cn%sub_triangle2%info%orientation = &
            & IBCLR(IBCLR(IBCLR(IEOR(up_down,iorientation),border),left),right)
          cn%sub_triangle3%info%orientation = &
            & IAND(IBCLR(iorientation,border),IBCLR(iorientation,left))
        ELSEIF (BTEST(iorientation,right)) THEN
          cn%sub_triangle0%info%orientation = iorientation
          cn%sub_triangle1%info%orientation = &
            & IAND(IBCLR(iorientation,border),IBCLR(iorientation,right))
          cn%sub_triangle2%info%orientation = &
            & IBCLR(IBCLR(IBCLR(IEOR(up_down,iorientation),border),left),right)
          cn%sub_triangle3%info%orientation = iorientation
        ELSE
          cn%sub_triangle0%info%orientation = IBCLR(iorientation,border)
          cn%sub_triangle1%info%orientation = IBCLR(iorientation,border)
          cn%sub_triangle2%info%orientation = IBCLR(IEOR(up_down,iorientation),border)
          cn%sub_triangle3%info%orientation = IBCLR(iorientation,border)
        END IF

      ELSE

        cn%sub_triangle0%info%orientation = iorientation
        cn%sub_triangle1%info%orientation = iorientation
        cn%sub_triangle2%info%orientation = IEOR(up_down,iorientation)
        cn%sub_triangle3%info%orientation = iorientation

      END IF

      ! - Now the terrible part - set the neighbors. I think there are
      !   three cases:
      !   1)  both neighbors on the left and right are of the same orientation,
      !       which are the pole triangles.
      !   2)  both neighbors on the left and right are of the opposite
      !       orientation, which are the triangles in the interior
      !   3)  one neighbor and the current have the same orientation and
      !       the neighbor on the other side is oriented opposite, which
      !       are the triangles on the edge off the original level -1 pole
      !       triangles.

      cn%sub_triangle2%neighbor0 => cn%sub_triangle3
      cn%sub_triangle2%neighbor1 => cn%sub_triangle0
      cn%sub_triangle2%neighbor2 => cn%sub_triangle1

      IF (BTEST(cn%sub_triangle0%info%orientation,pole)) THEN

        cn%sub_triangle0%neighbor0 => cn%neighbor0%sub_triangle0
        cn%sub_triangle0%neighbor1 => cn%sub_triangle2
        cn%sub_triangle0%neighbor2 => cn%neighbor2%sub_triangle0

        cn%sub_triangle1%neighbor0 => cn%neighbor0%sub_triangle3
        cn%sub_triangle1%neighbor1 => cn%neighbor1%sub_triangle3
        cn%sub_triangle1%neighbor2 => cn%sub_triangle2

        cn%sub_triangle3%neighbor0 => cn%sub_triangle2
        cn%sub_triangle3%neighbor1 => cn%neighbor1%sub_triangle1
        cn%sub_triangle3%neighbor2 => cn%neighbor2%sub_triangle1

      ELSE IF (BTEST(cn%sub_triangle0%info%orientation,border)) THEN

        IF (BTEST(iorientation,left)) THEN

          cn%sub_triangle0%neighbor0 => cn%neighbor0%sub_triangle0
          cn%sub_triangle0%neighbor1 => cn%sub_triangle2
          cn%sub_triangle0%neighbor2 => cn%neighbor2%sub_triangle3

          cn%sub_triangle1%neighbor0 => cn%neighbor0%sub_triangle3
          cn%sub_triangle1%neighbor1 => cn%neighbor1%sub_triangle3
          cn%sub_triangle1%neighbor2 => cn%sub_triangle2

          cn%sub_triangle3%neighbor0 => cn%sub_triangle2
          cn%sub_triangle3%neighbor1 => cn%neighbor1%sub_triangle1
          cn%sub_triangle3%neighbor2 => cn%neighbor2%sub_triangle0

        ELSEIF (BTEST(iorientation,right)) THEN

          cn%sub_triangle0%neighbor0 => cn%neighbor0%sub_triangle1
          cn%sub_triangle0%neighbor1 => cn%sub_triangle2
          cn%sub_triangle0%neighbor2 => cn%neighbor2%sub_triangle0

          cn%sub_triangle1%neighbor0 => cn%neighbor0%sub_triangle0
          cn%sub_triangle1%neighbor1 => cn%neighbor1%sub_triangle3
          cn%sub_triangle1%neighbor2 => cn%sub_triangle2

          cn%sub_triangle3%neighbor0 => cn%sub_triangle2
          cn%sub_triangle3%neighbor1 => cn%neighbor1%sub_triangle1
          cn%sub_triangle3%neighbor2 => cn%neighbor2%sub_triangle1

        ELSE

          cn%sub_triangle0%neighbor0 => cn%neighbor0%sub_triangle1
          cn%sub_triangle0%neighbor1 => cn%sub_triangle2
          cn%sub_triangle0%neighbor2 => cn%neighbor2%sub_triangle3

          cn%sub_triangle1%neighbor0 => cn%neighbor0%sub_triangle0
          cn%sub_triangle1%neighbor1 => cn%neighbor1%sub_triangle3
          cn%sub_triangle1%neighbor2 => cn%sub_triangle2

          cn%sub_triangle3%neighbor0 => cn%sub_triangle2
          cn%sub_triangle3%neighbor1 => cn%neighbor1%sub_triangle1
          cn%sub_triangle3%neighbor2 => cn%neighbor2%sub_triangle0

        END IF

      ELSE

        cn%sub_triangle0%neighbor0 => cn%neighbor0%sub_triangle1
        cn%sub_triangle0%neighbor1 => cn%sub_triangle2
        cn%sub_triangle0%neighbor2 => cn%neighbor2%sub_triangle3

        cn%sub_triangle1%neighbor0 => cn%neighbor0%sub_triangle0
        cn%sub_triangle1%neighbor1 => cn%neighbor1%sub_triangle3
        cn%sub_triangle1%neighbor2 => cn%sub_triangle2

        cn%sub_triangle3%neighbor0 => cn%sub_triangle2
        cn%sub_triangle3%neighbor1 => cn%neighbor1%sub_triangle1
        cn%sub_triangle3%neighbor2 => cn%neighbor2%sub_triangle0

      END IF

      CALL build_cell_graph(cn%sub_triangle0)
      CALL build_cell_graph(cn%sub_triangle1)
      CALL build_cell_graph(cn%sub_triangle2)
      CALL build_cell_graph(cn%sub_triangle3)

    ENDIF

  END SUBROUTINE build_cell_graph

  !
  !EOC
  !-------------------------------------------------------------------------
  !
  !

  !>
  !! This subroutine builds the complete topological structure.
  !!
  !! This subroutine builds the complete topological structure
  !! for the grid entities other than cells (edges, vertices),
  !! including the construction of their neighborhood
  !! relationships and the grid connectivities relating e.g. edges to
  !! cells, cells to vertices. The main burden of the work is related
  !! to the construction of the vertex connectivities, which are
  !! built starting with the <i>cchain</i> subroutine.
  !! This is also the slowest part of the code so far and might
  !! have to be improved in the future.
  !!
  !! @par Revision History
  !!  Luca Bonaventura, MPI-M, Hamburg, January 2005
  !!  Modification by Thomas Heinze, DWD, (2007-04-27):
  !!  - replaced double loop to set vertices at ends of edge by a single
  !!    loop over all vertices
  !!
  SUBROUTINE build_full_graph (sph)
    !


    TYPE(t_spheres) :: sph
    INTEGER :: i_not, i_nov, i_noe
    INTEGER :: j,jk,i0,i1,i2,ie_count,iv_count,i_count,j1,j2,jj
    INTEGER :: itemp(1:3),jtemp(1:3)
    INTEGER :: ieknown1,ieknown2,ieknown3
    INTEGER :: i_edge_idx
    INTEGER :: i_kk, istat
    LOGICAL :: lchecks


    INTEGER, ALLOCATABLE :: cell_chains(:,:)

    TYPE(t_triangle), POINTER :: ct(:) ! local array of triangles
    TYPE(edge)    , POINTER :: ce(:) ! local array of edges
    TYPE(vertex)  , POINTER :: cv(:) ! local array of vertices


    !-------------------------------------------------

    ! set dummy edge, vertex and cells(triangles)

    dummy_e%info%idx=-1
    dummy_c%info%idx=-1
    dummy_v%info%idx=-1

    dummy_e%vertex0=>dummy_v
    dummy_e%vertex1=>dummy_v
    dummy_v%triangle5=>dummy_c
    dummy_v%neighbor5=>dummy_v
    dummy_v%edge5=>dummy_e

    ! set number of edges, vertices and triangles(cells)

    i_not=sph%no_triangles
    i_nov=sph%no_vertices
    i_noe=sph%no_edges

    ! point local arrays to global one

    ct=> sph%ts
    ce=> sph%es
    cv=> sph%vs


    ALLOCATE(cell_chains(0:i_not-1,6), stat=istat)

    IF (istat >0) THEN
      WRITE (message_text,'(a)')    'Problem in allocating cell_chains'
      CALL message ('',TRIM(message_text))
      CALL finish ('build_full_graph', message_text)
    ENDIF

    cell_chains = -1

    ! good place to set absolute cell index on level

    DO j = 0, i_not-1
      ct(j)%info%idx = j
    ENDDO
    DO j = 0, i_noe-1
      ce(j)%info%idx = -1
    ENDDO
    DO j = 0, i_nov-1
      cv(j)%info%idx = -1
      cv(j)%edge0=>dummy_e
      cv(j)%edge1=>dummy_e
      cv(j)%edge2=>dummy_e
      cv(j)%edge3=>dummy_e
      cv(j)%edge4=>dummy_e
      cv(j)%edge5=>dummy_e
      cv(j)%triangle0=>dummy_c
      cv(j)%triangle1=>dummy_c
      cv(j)%triangle2=>dummy_c
      cv(j)%triangle3=>dummy_c
      cv(j)%triangle4=>dummy_c
      cv(j)%triangle5=>dummy_c
    ENDDO
    !
    ! first:  set edge indices
    !

    ie_count=0
    !   loop over triangular cells
    DO j=0,i_not-1
      !
      !        first check that edge has not been counted already
      !
      ieknown1=-1
      ieknown2=-1
      ieknown3=-1

      i0=ct(j)%neighbor0%info%idx
      i1=ct(j)%neighbor1%info%idx
      i2=ct(j)%neighbor2%info%idx

      !       DO jk=0,j-1  ! Unnecessary nested do loop, restructured for efficiency improvement

      !          IF (jk==i0) THEN
      IF ((i0>=0).and.(i0<=j-1)) THEN
        jk = i0

        IF(j==ct(jk)%neighbor0%info%idx) THEN

          ieknown1=ct(jk)%edge0%info%idx

        ELSEIF(j==ct(jk)%neighbor1%info%idx) THEN

          ieknown1=ct(jk)%edge1%info%idx

        ELSEIF(j==ct(jk)%neighbor2%info%idx) THEN

          ieknown1=ct(jk)%edge2%info%idx

        ENDIF
      ENDIF
      !         ELSEIF (jk==i1) THEN
      IF ((i1>=0).and.(i1<=j-1)) THEN
        jk = i1

        IF(j==ct(jk)%neighbor0%info%idx) THEN

          ieknown2=ct(jk)%edge0%info%idx

        ELSEIF(j==ct(jk)%neighbor1%info%idx) THEN

          ieknown2=ct(jk)%edge1%info%idx

        ELSEIF(j==ct(jk)%neighbor2%info%idx) THEN
          ieknown2=ct(jk)%edge2%info%idx
        ENDIF
      ENDIF
      !          ELSEIF (jk==i2) THEN
      IF ((i2>=0).and.(i2<=j-1)) THEN
        jk = i2
        IF(j==ct(jk)%neighbor0%info%idx) THEN
          ieknown3=ct(jk)%edge0%info%idx
        ELSEIF(j==ct(jk)%neighbor1%info%idx) THEN
          ieknown3=ct(jk)%edge1%info%idx
        ELSEIF(j==ct(jk)%neighbor2%info%idx) THEN
          ieknown3=ct(jk)%edge2%info%idx
        ENDIF
      ENDIF

      !       ENDDO
      !
      ! then add edge to the table
      !

      IF((ieknown1<0).and.(ieknown2<0).and.(ieknown3<0)) THEN

        ce(ie_count)%info%idx=ie_count
        ce(ie_count)%triangle0=>ct(j)
        ce(ie_count)%triangle1=>ct(j)%neighbor0
        ct(j)%edge0=>ce(ie_count)
        ie_count=ie_count+1

        ce(ie_count)%info%idx=ie_count
        ce(ie_count)%triangle0=>ct(j)
        ce(ie_count)%triangle1=>ct(j)%neighbor1
        ct(j)%edge1=>ce(ie_count)
        ie_count=ie_count+1

        ce(ie_count)%info%idx=ie_count
        ce(ie_count)%triangle0=>ct(j)
        ce(ie_count)%triangle1=>ct(j)%neighbor2
        ct(j)%edge2=>ce(ie_count)
        ie_count=ie_count+1

      ELSE IF((ieknown1>=0).and.(ieknown2<0).and.(ieknown3<0)) THEN

        ct(j)%edge0=>ce(ieknown1)

        ce(ie_count)%info%idx=ie_count
        ce(ie_count)%triangle0=>ct(j)
        ce(ie_count)%triangle1=>ct(j)%neighbor1
        ct(j)%edge1=>ce(ie_count)
        ie_count=ie_count+1

        ce(ie_count)%info%idx=ie_count
        ce(ie_count)%triangle0=>ct(j)
        ce(ie_count)%triangle1=>ct(j)%neighbor2
        ct(j)%edge2=>ce(ie_count)
        ie_count=ie_count+1

      ELSE IF((ieknown2 >=0).and.(ieknown1 <0).and.(ieknown3<0)) THEN
        ct(j)%edge1=>ce(ieknown2)

        ce(ie_count)%info%idx=ie_count
        ce(ie_count)%triangle0=>ct(j)
        ce(ie_count)%triangle1=>ct(j)%neighbor0

        ct(j)%edge0=>ce(ie_count)
        ie_count=ie_count+1

        ce(ie_count)%info%idx=ie_count
        ce(ie_count)%triangle0=>ct(j)
        ce(ie_count)%triangle1=>ct(j)%neighbor2
        ct(j)%edge2=>ce(ie_count)
        ie_count=ie_count+1

      ELSE IF((ieknown3>=0).and.(ieknown1<0).and.(ieknown2<0)) THEN
        ct(j)%edge2=>ce(ieknown3)
        ce(ie_count)%info%idx=ie_count
        ce(ie_count)%triangle0=>ct(j)
        ce(ie_count)%triangle1=>ct(j)%neighbor0

        ct(j)%edge0=>ce(ie_count)
        ie_count=ie_count+1

        ce(ie_count)%info%idx=ie_count
        ce(ie_count)%triangle0=>ct(j)
        ce(ie_count)%triangle1=>ct(j)%neighbor1
        ct(j)%edge1=>ce(ie_count)
        ie_count=ie_count+1


      ELSE IF((ieknown1>=0).and.(ieknown2>=0).and.(ieknown3<0)) THEN

        ct(j)%edge0=>ce(ieknown1)
        ct(j)%edge1=>ce(ieknown2)

        ce(ie_count)%info%idx=ie_count
        ce(ie_count)%triangle0=>ct(j)
        ce(ie_count)%triangle1=>ct(j)%neighbor2
        ct(j)%edge2=>ce(ie_count)
        ie_count=ie_count+1

      ELSE IF((ieknown2>=0).and.(ieknown3>=0).and.(ieknown1<0)) THEN
        ct(j)%edge1=>ce(ieknown2)
        ct(j)%edge2=>ce(ieknown3)

        ce(ie_count)%info%idx=ie_count
        ce(ie_count)%triangle0=>ct(j)
        ce(ie_count)%triangle1=>ct(j)%neighbor0
        ct(j)%edge0=>ce(ie_count)
        ie_count=ie_count+1


      ELSE IF(( ieknown1 >=0).and.(ieknown3 >=0).and.(ieknown2<0)) THEN


        ct(j)%edge0=>ce(ieknown1)
        ct(j)%edge2=>ce(ieknown3)

        ce(ie_count)%info%idx=ie_count
        ce(ie_count)%triangle0=>ct(j)
        ce(ie_count)%triangle1=>ct(j)%neighbor1
        ct(j)%edge1=>ce(ie_count)
        ie_count=ie_count+1

      ELSE IF((ieknown1>=0).and.(ieknown2>=0).and.(ieknown3>=0)) THEN

        ct(j)%edge0=>ce(ieknown1)

        ct(j)%edge1=>ce(ieknown2)

        ct(j)%edge2=>ce(ieknown3)

        CONTINUE

      ENDIF


    ENDDO
    !
    ! second:  set vertex indices. This is the most complex part
    !

    !
    ! determine chains of cells around vertices
    !
    DO j=0,i_not-1
      CALL cchain(j, cell_chains(j,:))
    ENDDO

    iv_count=0

    DO j=0,i_not-1

      !
      !  check if cell chain corresponds to already counted vertex
      !
      lchecks=.false.
      ! The only points we have to check for double counts are those
      ! belonging to the current cell chain; a cell chain not originating
      ! at one cell of the current cell chain cannot be identical to it!
      ! Compared to the previous implementation, this reduces the computing time
      ! by several orders of magnitude
      DO jj=1,6
        jk = cell_chains(j,jj)
        IF ((jk < 0) .or. (jk >= j)) CYCLE
        i_count=0
        DO j1=1,6
          DO j2=1,6
            IF (cell_chains(j,j1)==cell_chains(jk,j2)) THEN
              i_count=i_count+1
              EXIT
            ENDIF
          ENDDO
          IF (i_count < j1) EXIT ! No chance to get 6 counts if one index does not occur
          ! in the other sample.
        ENDDO
        IF (i_count==6) THEN
          lchecks=.true.
          EXIT
        ENDIF
      ENDDO

      IF (lchecks) THEN
        CONTINUE
      ELSE
        !
        !  add vertex if it had not been counted before
        !
        cv(iv_count)%info%idx=iv_count

        cv(iv_count)%triangle0=>ct(cell_chains(j,1))
        cv(iv_count)%triangle1=>ct(cell_chains(j,2))
        cv(iv_count)%triangle2=>ct(cell_chains(j,3))
        cv(iv_count)%triangle3=>ct(cell_chains(j,4))
        cv(iv_count)%triangle4=>ct(cell_chains(j,5))



        IF(cell_chains(j,6)>=0) THEN
          cv(iv_count)%triangle5=>ct(cell_chains(j,6))
        ELSE
          !
          !   special case of pentagonal dual cell
          !
          !         IF((levels==0).AND.(j==0)) write(*,*) 'here'

          cv(iv_count)%triangle5=>dummy_c
          cv(iv_count)%neighbor5=>dummy_v
          cv(iv_count)%edge5=>dummy_e
          !       write(*,*) iv_count,ASSOCIATED(sph%vs(iv_count)%neighbor5)

        ENDIF
        !   write(*,*) ASSOCIATED(sph%vs(0)%neighbor5)
        !

        !     set edges of dual cell
        !

        itemp(1)=ct(cell_chains(j,1))%edge0%info%idx
        itemp(2)=ct(cell_chains(j,1))%edge1%info%idx
        itemp(3)=ct(cell_chains(j,1))%edge2%info%idx

        jtemp(1)=ct(cell_chains(j,2))%edge0%info%idx
        jtemp(2)=ct(cell_chains(j,2))%edge1%info%idx
        jtemp(3)=ct(cell_chains(j,2))%edge2%info%idx
        i_kk=find_edge(itemp(:),jtemp(:))
        cv(iv_count)%edge0=>ce(i_kk)



        itemp(1)=ct(cell_chains(j,2))%edge0%info%idx
        itemp(2)=ct(cell_chains(j,2))%edge1%info%idx
        itemp(3)=ct(cell_chains(j,2))%edge2%info%idx


        jtemp(1)=ct(cell_chains(j,3))%edge0%info%idx
        jtemp(2)=ct(cell_chains(j,3))%edge1%info%idx
        jtemp(3)=ct(cell_chains(j,3))%edge2%info%idx

        i_kk=find_edge(itemp(:),jtemp(:))
        cv(iv_count)%edge1=>ce(i_kk)


        itemp(1)=ct(cell_chains(j,3))%edge0%info%idx
        itemp(2)=ct(cell_chains(j,3))%edge1%info%idx
        itemp(3)=ct(cell_chains(j,3))%edge2%info%idx

        jtemp(1)=ct(cell_chains(j,4))%edge0%info%idx
        jtemp(2)=ct(cell_chains(j,4))%edge1%info%idx
        jtemp(3)=ct(cell_chains(j,4))%edge2%info%idx

        i_kk= find_edge(itemp,jtemp)

        cv(iv_count)%edge2=>ce(i_kk)


        itemp(1)=ct(cell_chains(j,4))%edge0%info%idx
        itemp(2)=ct(cell_chains(j,4))%edge1%info%idx
        itemp(3)=ct(cell_chains(j,4))%edge2%info%idx

        jtemp(1)=ct(cell_chains(j,5))%edge0%info%idx
        jtemp(2)=ct(cell_chains(j,5))%edge1%info%idx
        jtemp(3)=ct(cell_chains(j,5))%edge2%info%idx

        i_kk= find_edge(itemp,jtemp)

        cv(iv_count)%edge3=>ce(i_kk)



        IF(cell_chains(j,6)<0) THEN
          !
          !   special case of pentagonal dual cell
          !

          itemp(1)=ct(cell_chains(j,5))%edge0%info%idx
          itemp(2)=ct(cell_chains(j,5))%edge1%info%idx
          itemp(3)=ct(cell_chains(j,5))%edge2%info%idx

          jtemp(1)=ct(cell_chains(j,1))%edge0%info%idx
          jtemp(2)=ct(cell_chains(j,1))%edge1%info%idx
          jtemp(3)=ct(cell_chains(j,1))%edge2%info%idx

          i_kk=find_edge(itemp,jtemp)

          cv(iv_count)%edge4=>ce(i_kk)
          cv(iv_count)%edge5=>dummy_e


        ELSE
          !
          !    case of hexagonal dual cell
          !

          itemp(1)=ct(cell_chains(j,5))%edge0%info%idx
          itemp(2)=ct(cell_chains(j,5))%edge1%info%idx
          itemp(3)=ct(cell_chains(j,5))%edge2%info%idx

          jtemp(1)=ct(cell_chains(j,6))%edge0%info%idx
          jtemp(2)=ct(cell_chains(j,6))%edge1%info%idx
          jtemp(3)=ct(cell_chains(j,6))%edge2%info%idx

          i_kk= find_edge(itemp,jtemp)

          cv(iv_count)%edge4=>ce(i_kk)

          itemp(1)=ct(cell_chains(j,6))%edge0%info%idx
          itemp(2)=ct(cell_chains(j,6))%edge1%info%idx
          itemp(3)=ct(cell_chains(j,6))%edge2%info%idx

          jtemp(1)=ct(cell_chains(j,1))%edge0%info%idx
          jtemp(2)=ct(cell_chains(j,1))%edge1%info%idx
          jtemp(3)=ct(cell_chains(j,1))%edge2%info%idx

          i_kk=find_edge(itemp,jtemp)

          cv(iv_count)%edge5=>ce(i_kk)


        ENDIF


        iv_count=iv_count+1

      ENDIF

    ENDDO

    !
    ! set vertices at ends of edge
    !

    DO j=0,i_noe-1
      ce(j)%vertex0 => dummy_v
      ce(j)%vertex1 => dummy_v
    ENDDO

    ! loop over all vertices (they already now the adjecent edges)

    DO j=0,i_nov-1

      ! edge0

      i_edge_idx = cv(j)%edge0%info%idx

      ! check if edge is asociated to another vertex

      IF(ce(i_edge_idx)%vertex0%info%idx < 0) THEN
        ce(i_edge_idx)%vertex0 => cv(j)
      ELSE   ! edge already associated to a vertex
        ce(i_edge_idx)%vertex1 => cv(j)
      ENDIF

      ! edge1

      i_edge_idx = cv(j)%edge1%info%idx

      ! check if edge is asociated to another vertex

      IF(ce(i_edge_idx)%vertex0%info%idx < 0) THEN
        ce(i_edge_idx)%vertex0 => cv(j)
      ELSE   ! edge already associated to a vertex
        ce(i_edge_idx)%vertex1 => cv(j)
      ENDIF

      ! edge2

      i_edge_idx = cv(j)%edge2%info%idx

      ! check if edge is asociated to another vertex

      IF(ce(i_edge_idx)%vertex0%info%idx < 0) THEN
        ce(i_edge_idx)%vertex0 => cv(j)
      ELSE   ! edge already associated to a vertex
        ce(i_edge_idx)%vertex1 => cv(j)
      ENDIF

      ! edge3

      i_edge_idx = cv(j)%edge3%info%idx

      ! check if edge is asociated to another vertex

      IF(ce(i_edge_idx)%vertex0%info%idx < 0) THEN
        ce(i_edge_idx)%vertex0 => cv(j)
      ELSE   ! edge already associated to a vertex
        ce(i_edge_idx)%vertex1 => cv(j)
      ENDIF

      ! edge4

      i_edge_idx = cv(j)%edge4%info%idx

      ! check if edge is asociated to another vertex

      IF(ce(i_edge_idx)%vertex0%info%idx < 0) THEN
        ce(i_edge_idx)%vertex0 => cv(j)
      ELSE   ! edge already associated to a vertex
        ce(i_edge_idx)%vertex1 => cv(j)
      ENDIF

      ! edge5, special case (might not exist)

      i_edge_idx = cv(j)%edge5%info%idx

      IF (i_edge_idx > -1) THEN

        ! check if edge is asociated to another vertex

        IF(ce(i_edge_idx)%vertex0%info%idx < 0) THEN
          ce(i_edge_idx)%vertex0 => cv(j)
        ELSE   ! edge already associated to a vertex
          ce(i_edge_idx)%vertex1 => cv(j)
        ENDIF
      ENDIF

    ENDDO

    !
    !  set neighbors of vertex (only necessary for two time level schemes!!)
    !

    DO j=0,i_nov-1
      IF(cv(j)%edge0%vertex0%info%idx==j) THEN
        cv(j)%neighbor0=>cv(j)%edge0%vertex1
      ELSE
        cv(j)%neighbor0=>cv(j)%edge0%vertex0
      ENDIF
      IF(cv(j)%edge1%vertex0%info%idx==j) THEN
        cv(j)%neighbor1=>cv(j)%edge1%vertex1
      ELSE
        cv(j)%neighbor1=>cv(j)%edge1%vertex0
      ENDIF
      IF(cv(j)%edge2%vertex0%info%idx==j) THEN
        cv(j)%neighbor2=>cv(j)%edge2%vertex1
      ELSE
        cv(j)%neighbor2=>cv(j)%edge2%vertex0
      ENDIF
      IF(cv(j)%edge3%vertex0%info%idx==j) THEN
        cv(j)%neighbor3=>cv(j)%edge3%vertex1
      ELSE
        cv(j)%neighbor3=>cv(j)%edge3%vertex0
      ENDIF
      IF(cv(j)%edge4%vertex0%info%idx==j) THEN
        cv(j)%neighbor4=>cv(j)%edge4%vertex1
      ELSE
        cv(j)%neighbor4=>cv(j)%edge4%vertex0
      ENDIF
      IF((cv(j)%edge5%vertex0%info%idx==j).and.(cv(j)%edge5%info%idx>=0)) THEN
        cv(j)%neighbor5=>cv(j)%edge5%vertex1
      ELSE IF((cv(j)%edge5%vertex1%info%idx==j).and.(cv(j)%edge5%info%idx>=0)) THEN
        cv(j)%neighbor5=>cv(j)%edge5%vertex0
      ELSE
        cv(j)%neighbor5=>dummy_v
      ENDIF
    ENDDO

    !
    !
    !

    DO j=0,i_not-1

      IF(ct(j)%edge0%vertex0%info%idx==ct(j)%edge2%vertex1%info%idx) THEN
        ct(j)%vertex0=>ct(j)%edge0%vertex0
        ct(j)%vertex1=>ct(j)%edge0%vertex1
        ct(j)%vertex2=>ct(j)%edge2%vertex0
      ELSE IF(ct(j)%edge0%vertex0%info%idx==ct(j)%edge2%vertex0%info%idx) THEN
        ct(j)%vertex0=>ct(j)%edge0%vertex0
        ct(j)%vertex1=>ct(j)%edge0%vertex1
        ct(j)%vertex2=>ct(j)%edge2%vertex1
      ELSE IF(ct(j)%edge0%vertex1%info%idx==ct(j)%edge2%vertex0%info%idx) THEN
        ct(j)%vertex0=>ct(j)%edge0%vertex1
        ct(j)%vertex1=>ct(j)%edge0%vertex0
        ct(j)%vertex2=>ct(j)%edge2%vertex1
      ELSE IF(ct(j)%edge0%vertex1%info%idx==ct(j)%edge2%vertex1%info%idx) THEN
        ct(j)%vertex0=>ct(j)%edge0%vertex1
        ct(j)%vertex1=>ct(j)%edge0%vertex0
        ct(j)%vertex2=>ct(j)%edge2%vertex0
      ELSE
        WRITE (message_text,'(a,i3)') &
          & 'ERROR IN SETTING VERTEX NEIGHBOURS FOR TRIANGLE ', j
        CALL message ('',TRIM(message_text))
      ENDIF

    ENDDO

    DEALLOCATE(cell_chains)


  CONTAINS

    !-------------------------------------------------------------------------

    !>
    !!               This subroutine builds a chain of grid cells.
    !!
    !!               This subroutine builds a chain of grid cells
    !!               around each vertex, which then allows to build
    !!               vertex neighbourhood relationships and
    !!               vertex connectitvities with the other grid entities.
    !!               It uses ONLY topological information, which is neat but leads
    !!               to slow performance compared to the previous grid generator.
    !!               This might require improvement in the future.
    !!
    !! @par Revision History
    !!  Luca Bonaventura, MPI-M, Hamburg, (2005-03)
    !!
    SUBROUTINE cchain(j,k_ii)

      IMPLICIT NONE

      !INTEGER, INTENT(OUT) :: k_ii(:) ! indexes of cells around a vertex
      INTEGER, INTENT(inout) :: k_ii(:) ! indexes of cells around a vertex
      INTEGER, INTENT(in) ::  j       ! vertex index
      INTEGER :: itemp0,itemp1,itemp2
      INTEGER :: i0,i1,i00,i01,i10,i11,i000,i001,i010,i011
      LOGICAL:: lchain

      !-------------------------------------------------------------------------
      !BOC

      ! initialization

      i0   = 0
      i1   = 0
      i00  = 0
      i01  = 0
      i10  = 0
      i11  = 0
      i000 = 0
      i001 = 0
      i010 = 0
      i011 = 0

      k_ii=-1

      ! first cells of lchain

      lchain=.false.

      i0=ct(j)%neighbor0%info%idx
      i1=ct(j)%neighbor2%info%idx

      k_ii(2)=j
      k_ii(3)=i0

      itemp0=ct(i0)%neighbor0%info%idx
      itemp1=ct(i0)%neighbor1%info%idx
      itemp2=ct(i0)%neighbor2%info%idx
      !
      !    neighbours of i0
      !
      IF(itemp0==j) THEN
        i00=itemp1
        i01=itemp2
      ELSE IF (itemp1==j) THEN
        i00=itemp0
        i01=itemp2
      ELSE IF (itemp2==j) THEN
        i00=itemp0
        i01=itemp1
      ELSE
        WRITE(*,*) 'error'
      ENDIF
      !
      !    neighbours of i00
      !
      itemp0=ct(i00)%neighbor0%info%idx
      itemp1=ct(i00)%neighbor1%info%idx
      itemp2=ct(i00)%neighbor2%info%idx

      IF(itemp0==i0) THEN
        i000=itemp1
        i001=itemp2
      ELSE IF (itemp1==i0) THEN
        i000=itemp0
        i001=itemp2
      ELSE IF (itemp2==i0) THEN
        i000=itemp0
        i001=itemp1
      ELSE
        WRITE(*,*) 'error'
      ENDIF
      !
      !    neighbours of i01
      !

      itemp0=ct(i01)%neighbor0%info%idx
      itemp1=ct(i01)%neighbor1%info%idx
      itemp2=ct(i01)%neighbor2%info%idx

      IF(itemp0==i0) THEN
        i010=itemp1
        i011=itemp2
      ELSE IF (itemp1==i0) THEN
        i010=itemp0
        i011=itemp2
      ELSE IF (itemp2==i0) THEN
        i010=itemp0
        i011=itemp1
      ELSE
        WRITE(*,*) 'error'
      ENDIF
      !
      !    neighbours of i000: first check if chain is closed
      !
      itemp0=ct(i000)%neighbor0%info%idx
      itemp1=ct(i000)%neighbor1%info%idx
      itemp2=ct(i000)%neighbor2%info%idx

      IF((itemp0==i1).or.(itemp1==i1).or.(itemp2==i1)) THEN
        k_ii(1)=i1
        k_ii(5)=i000
        k_ii(4)=i00
        lchain=.true.
      ENDIF
      !
      !    neighbours of i001: first check if chain is closed
      !

      itemp0=ct(i001)%neighbor0%info%idx
      itemp1=ct(i001)%neighbor1%info%idx
      itemp2=ct(i001)%neighbor2%info%idx

      IF((itemp0==i1).or.(itemp1==i1).or.(itemp2==i1)) THEN
        k_ii(1)=i1
        k_ii(5)=i001
        k_ii(4)=i00
        lchain=.true.
      ENDIF
      !
      !    neighbours of i010: first check if chain is closed
      !

      itemp0=ct(i010)%neighbor0%info%idx
      itemp1=ct(i010)%neighbor1%info%idx
      itemp2=ct(i010)%neighbor2%info%idx

      IF((itemp0==i1).or.(itemp1==i1).or.(itemp2==i1)) THEN
        k_ii(1)=i1
        k_ii(5)=i010
        k_ii(4)=i01
        lchain=.true.
      ENDIF
      !
      !    neighbours of i010: first check if chain is closed
      !
      itemp0=ct(i011)%neighbor0%info%idx
      itemp1=ct(i011)%neighbor1%info%idx
      itemp2=ct(i011)%neighbor2%info%idx

      IF((itemp0==i1).or.(itemp1==i1).or.(itemp2==i1)) THEN
        k_ii(1)=i1
        k_ii(5)=i011
        k_ii(4)=i01
        lchain=.true.
      ENDIF

      IF (.not.(lchain)) THEN


        !
        !    if chain is not closed look around i1
        !

        itemp0=ct(i1)%neighbor0%info%idx
        itemp1=ct(i1)%neighbor1%info%idx
        itemp2=ct(i1)%neighbor2%info%idx

        IF(itemp0==j) THEN
          i10=itemp1
          i11=itemp2
        ELSE IF (itemp1==j) THEN
          i10=itemp0
          i11=itemp2
        ELSE IF (itemp2==j) THEN
          i10=itemp0
          i11=itemp1
        ELSE
          WRITE(*,*) 'error'
        ENDIF




        !
        !  now second round of checks
        !

        !
        !    neighbours of i1
        !
        itemp0=ct(i10)%neighbor0%info%idx
        itemp1=ct(i10)%neighbor1%info%idx
        itemp2=ct(i10)%neighbor2%info%idx

        IF((itemp0==i000).or.(itemp1==i000).or.(itemp2==i000)) THEN
          k_ii(6)=i10
          k_ii(1)=i1
          k_ii(5)=i000
          k_ii(4)=i00
          lchain=.true.
        ENDIF

        IF((itemp0==i001).or.(itemp1==i001).or.(itemp2==i001)) THEN
          k_ii(6)=i10
          k_ii(1)=i1
          k_ii(5)=i001
          k_ii(4)=i00
          lchain=.true.
        ENDIF

        IF((itemp0==i010).or.(itemp1==i010).or.(itemp2==i010)) THEN
          k_ii(6)=i10
          k_ii(1)=i1
          k_ii(5)=i010
          k_ii(4)=i00
          lchain=.true.
        ENDIF

        IF((itemp0==i011).or.(itemp1==i011).or.(itemp2==i011)) THEN
          k_ii(6)=i10
          k_ii(1)=i1
          k_ii(5)=i011
          k_ii(4)=i00
          lchain=.true.
        ENDIF

        itemp0=ct(i11)%neighbor0%info%idx
        itemp1=ct(i11)%neighbor1%info%idx
        itemp2=ct(i11)%neighbor2%info%idx

        IF((itemp0==i000).or.(itemp1==i000).or.(itemp2==i000)) THEN
          k_ii(6)=i11
          k_ii(1)=i1
          k_ii(5)=i000
          k_ii(4)=i00
          lchain=.true.
        ENDIF

        IF((itemp0==i001).or.(itemp1==i001).or.(itemp2==i001)) THEN
          k_ii(6)=i11
          k_ii(1)=i1
          k_ii(5)=i001
          k_ii(4)=i00
          lchain=.true.
        ENDIF

        IF((itemp0==i010).or.(itemp1==i010).or.(itemp2==i010)) THEN
          k_ii(6)=i11
          k_ii(1)=i1
          k_ii(5)=i010
          k_ii(4)=i00
          lchain=.true.
        ENDIF

        IF((itemp0==i011).or.(itemp1==i011).or.(itemp2==i011)) THEN
          k_ii(6)=i11
          k_ii(1)=i1
          k_ii(5)=i011
          k_ii(4)=i00
          lchain=.true.
        ENDIF


      ENDIF

    END SUBROUTINE cchain

    !-------------------------------------------------------------------------

    !>
    !!               This subroutine checks if two chains of grid cells.
    !!
    !!               This subroutine checks if two chains of grid cells
    !! around vertices coincide
    !!
    !! @par Revision History
    !!  Luca Bonaventura, MPI-M, Hamburg, (2005-03)
    !!
    SUBROUTINE check_chain(k_a,k_b,lchecked)



      IMPLICIT NONE
      INTEGER ::ji,j,i_count
      INTEGER,INTENT(in):: k_a(:),k_b(:)
      LOGICAL, INTENT(inout) :: lchecked

      !-------------------------------------------------



      i_count=0

      DO ji=1,SIZE(k_a)

        DO j=1,SIZE(k_b)

          IF (k_a(ji)==k_b(j)) THEN
            i_count=i_count+1
            EXIT
          ENDIF

        ENDDO
        IF (i_count < ji) RETURN ! No chance to get 6 counts if one index
        ! does not occur in the other sample
      ENDDO
      IF(i_count==SIZE(k_a)) lchecked=.true.

    END SUBROUTINE check_chain

    !-------------------------------------------------------------------------

    !>
    !!               This function finds the edge that is.
    !!
    !!               This function finds the edge that is
    !!               common to two adjacent grid cells.
    !!
    !! @par Revision History
    !!  Luca Bonaventura, MPI-M, Hamburg, (2005-03)
    !!
    FUNCTION find_edge(k_ite,k_jte) result(k2)




      IMPLICIT NONE
      INTEGER :: k2
      INTEGER ::   ji1,jj1
      INTEGER, INTENT(in) :: k_ite(3),k_jte(3)
      LOGICAL ::lequal
      !-------------------------------------------------
      k2 = 0
      lequal=.false.
      DO ji1=1,3
        DO jj1=1,3
          IF (k_ite(ji1)==k_jte(jj1)) THEN
            k2=k_jte(jj1)
            lequal=.true.
          ENDIF
          IF(lequal) EXIT
        ENDDO
        IF(lequal) EXIT
      ENDDO


    END FUNCTION find_edge

  END SUBROUTINE build_full_graph


  !-------------------------------------------------------------------------
  !>
  !!               Prints information on triangles of the <i>root</i> grid.
  !!
  !!               Produces a huge amount of data, to be called only for debugging
  !!               purposes.
  !!
  !! @par Revision History
  !!  Luis Kornblueh, MPI-M, Hamburg, (2004-03)
  !!
  SUBROUTINE print_root

    CHARACTER(LEN=4) :: conf_self, conf_left, conf_right, conf_base
    INTEGER :: ji, j
    !-------------------------------------------------
    !BOC
    WRITE (nerr,*) '----------------------------------------------------------'
    WRITE (nerr,*) ' root triangles neighborhood relations'
    DO j = 0, 19
      DO ji = 0, kroot*kroot-1
        CALL config_string (root(ji,j)%t%info%orientation,           conf_self)
        CALL config_string (root(ji,j)%t%neighbor0%info%orientation, conf_left)
        CALL config_string (root(ji,j)%t%neighbor2%info%orientation, conf_right)
        CALL config_string (root(ji,j)%t%neighbor1%info%orientation, conf_base)
        WRITE (nerr,*) '------------------------------------------------------'
        WRITE (nerr,'(3(a,i4.4,1x,a4))')                                               &
          & ' left  = ', root(ji,j)%t%neighbor0%info%configuration, conf_left,       &
          & ' self  = ', root(ji,j)%t%info%configuration , conf_self,                &
          & ' right = ', root(ji,j)%t%neighbor2%info%configuration, conf_right
        WRITE (nerr,'(a,i4.4,1x,a4)') '                   base  = ',                   &
          & root(ji,j)%t%neighbor1%info%configuration, conf_base
      END DO
    END DO

  END SUBROUTINE print_root
  !EOC
  !-------------------------------------------------------------------------
  !>
  !!               Defines for each.
  !!
  !!               Defines for each
  !!               triangle on the level 1 <i>root</i> grid
  !!               a string that summarizes the main triangle
  !!               characteristics.
  !!
  !! @par Revision History
  !!  Luis Kornblueh, MPI-M, Hamburg, (2004-03)
  !!
  SUBROUTINE config_string (k_iorientation, string)
    INTEGER,          INTENT(in)  :: k_iorientation
    CHARACTER(LEN=*), INTENT(out) :: string
    !-------------------------------------------------
    string = '----'

    IF (BTEST(k_iorientation,up))     string(4:4) = 'u'
    IF (BTEST(k_iorientation,down))   string(4:4) = 'd'
    IF (BTEST(k_iorientation,pole))   string(3:3) = 'p'
    IF (BTEST(k_iorientation,main))   string(3:3) = 'm'
    IF (BTEST(k_iorientation,border)) string(2:2) = 'b'
    IF (BTEST(k_iorientation,right))  string(1:1) = 'r'
    IF (BTEST(k_iorientation,left))   string(1:1) = 'l'

  END SUBROUTINE config_string

END MODULE mo_topology

!------------------------------------













