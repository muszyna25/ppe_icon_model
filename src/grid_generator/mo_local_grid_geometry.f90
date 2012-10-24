!-------------------------------------------------------------------------------------
! mo_local_grid_geometry first implementation
!>
!! Computes the geometry of a triangular grid on the sphere
!!
!! @author Leonidas Linardakis, MPI-M
!!
!! @par Revision History
!!   First implementation by Leonidas Linardakis, MPI-M, 2010-02-01
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
!-------------------------------------------------------------------------------------
MODULE mo_local_grid_geometry
#include "grid_definitions.inc"

  !-------------------------------------------------------------------------
  USE mo_kind,               ONLY: wp
  USE mo_math_constants,     ONLY: rad2deg
  USE mo_exception,          ONLY: finish! , message
  USE mo_local_grid
  USE mo_base_geometry,  ONLY: t_cartesian_coordinates, vector_product, &
    & circum_center, cc2gc, arc_length,  triangle_area, &
    & inter_section, t_geographical_coordinates, angle_of_vectors,  &
    & sphere_cartesian_midpoint, cartesian_to_geographical,         &
    & geographical_to_cartesian
  USE mo_io_units,       ONLY: nnml, filename_max
  USE mo_namelist,       ONLY: position_nml, open_nml, positioned
  USE mo_timer,          ONLY: new_timer, timer_start, timer_stop, print_timer, delete_timer
  USE mo_io_local_grid,  ONLY: read_new_netcdf_grid, write_netcdf_grid

  IMPLICIT NONE

  PRIVATE

  ! !VERSION CONTROL:
  CHARACTER(LEN=*), PARAMETER :: version = '$Id$'

  PUBLIC :: compute_sphere_geometry
  PUBLIC :: compute_sphere_grid_geometry
  PUBLIC :: order_cell_connectivity     ! Reorders the cell vertices, edges, neigbors
  PUBLIC :: get_cell_barycenters
  PUBLIC :: get_triangle_circumcenters
    ! in counterclock way
  PUBLIC :: get_common_edge_vertex
  PUBLIC :: edges_cell_angle, edges_normal_angle, edges_3D_angle, print_cell_angles
  PUBLIC :: get_celledge_vertices
  PUBLIC :: use_barycenters, use_cartesian_centers
  PUBLIC :: NO_ANGLE
  !---------------------------------------------------------------------------

  ! !DEFINED PARAMETERS:
  !
  INTEGER, PARAMETER ::  use_barycenters = 1
  INTEGER, PARAMETER ::  use_cartesian_centers = 2
  INTEGER, PARAMETER ::  lonlat_rad_coordinates = 1
  INTEGER, PARAMETER ::  cartesian_meter_coordinates = 2
  REAL(wp), PARAMETER  :: NO_ANGLE = 1000.0_wp


CONTAINS


  !-------------------------------------------------------------------------
  !>
  !!Computes the sphere geometry for a netcdf grid
  !-------------------------------------------------------------------------
  SUBROUTINE compute_sphere_geometry(param_file_name)
    CHARACTER(LEN=*), INTENT(in) :: param_file_name
    
    CHARACTER(LEN=filename_max) :: in_file, out_file
    INTEGER :: grid_id
    INTEGER :: i_status
    
    NAMELIST /file_names/ in_file, out_file

    in_file=""
    out_file=""
    CALL open_nml(param_file_name)
    CALL position_nml('file_names',STATUS=i_status)
    IF (i_status == positioned) THEN
      READ (nnml,file_names)
    ELSE
      CALL finish("compute_sphere_geometry","file_names section not found")
    ENDIF
    
    grid_id = read_new_netcdf_grid(in_file, read_grid_ids=.true.)
    CALL set_default_geometry_parameters(to_grid_id=grid_id, param_file_name=param_file_name, &
      & from_grid_id=grid_id)
    CALL compute_sphere_grid_geometry(grid_id)
    
    CALL write_netcdf_grid(grid_id, out_file)
    CALL delete_grid(grid_id)
    RETURN

  END SUBROUTINE compute_sphere_geometry
  !-------------------------------------------------------------------------

  !---------------------------------------------------------------------------
  SUBROUTINE compute_sphere_grid_geometry(in_grid_id)
    INTEGER, INTENT(in) :: in_grid_id

    TYPE(t_grid), POINTER :: compute_grid
    TYPE(t_grid_cells), POINTER :: cells
    TYPE(t_grid_edges), POINTER :: edges
    TYPE(t_grid_vertices), POINTER :: verts

    TYPE(t_cartesian_coordinates) :: cartesian_v(3)
    TYPE(t_cartesian_coordinates) :: cartesian_c(2)
    TYPE(t_cartesian_coordinates) :: cartesian_center,edge_vector, edge_normal_vector, x,y
    TYPE(t_cartesian_coordinates) :: circumcenters_vector,tmp_vector

    REAL(wp) :: real_tmp,sphere_radius,sphere_radius_squared,lon,lat

    INTEGER :: no_of_cells, no_of_edges, no_of_verts
    INTEGER :: i,j,cell_index,edge_index,vertex_index
    INTEGER :: cell_1, cell_2, vertex_1, vertex_2
    ! INTEGER :: v1,v2, k
! !$  INTEGER OMP_GET_NUM_THREADS

    INTEGER :: timer_set_sphere_geom_grid

    timer_set_sphere_geom_grid = new_timer("compute_sphere_grid_geometry")
    CALL timer_start(timer_set_sphere_geom_grid)


    compute_grid => get_grid(in_grid_id)
    verts=>compute_grid%verts
    edges=>compute_grid%edges
    cells=>compute_grid%cells
    no_of_cells = cells%no_of_existcells
    no_of_edges = edges%no_of_existedges
    no_of_verts = verts%no_of_existvertices
    compute_grid%geometry_type = sphere_geometry

    IF (cells%max_no_of_vertices /= 3) THEN
      CALL finish('compute_sphere_grid_geometry','Not a triangular grid')
    ENDIF

!$OMP PARALLEL
! !$  PRINT*,  "OMP_GET_NUM_THREADS=", OMP_GET_NUM_THREADS()
!$OMP DO PRIVATE(vertex_index)
    ! get the vertices geocoordinates from the cartesian
    DO vertex_index=1,no_of_verts
        verts%vertex(vertex_index) = cc2gc(verts%cartesian(vertex_index))
    ENDDO
!$OMP END DO

!$OMP DO PRIVATE(cell_index, cartesian_v)
    DO cell_index = 1, no_of_cells
      cartesian_v(1) = verts%cartesian(cells%get_vertex_index(cell_index,1))
      cartesian_v(2) = verts%cartesian(cells%get_vertex_index(cell_index,2))
      cartesian_v(3) = verts%cartesian(cells%get_vertex_index(cell_index,3))
      ! Compute triangle centers
      cells%cartesian_center(cell_index) = &
        & circum_center(cartesian_v(1), cartesian_v(2), cartesian_v(3))

      ! Compute triangle areas
      cells%area(cell_index) = triangle_area(cartesian_v(1),cartesian_v(2),cartesian_v(3))
      cells%center(cell_index) = cc2gc(cells%cartesian_center(cell_index))
    ENDDO ! cellIndex = 1, noOfInputCells
!$OMP END DO

!$OMP DO PRIVATE(cell_index,edge_index,i,j)
    DO cell_index = 1, no_of_cells
      ! compute the cell  orientation
      DO i=1,cells%max_no_of_vertices
        edge_index=cells%get_edge_index(cell_index,i)

        IF (cell_index == edges%get_cell_index(edge_index,1)) THEN
          j = 1
          ! compute edge orientation according to Gauss formula (positive points outwards)
          cells%get_edge_orient(cell_index,i) = 1
        ELSE
          j = 2
          cells%get_edge_orient(cell_index,i) = -1
        ENDIF

      ENDDO

    ENDDO
!$OMP END DO

    ! the following cell geometric properties have been calculated
    !    gridCells%center
    !    gridCells%area
    !    gridCells%edge_orient

    !   write(*,*) "grid_sphericalTriangleGridGeometry: compute edges..."
    !   call flush(6)

    ! compute edges geometric properties
!$OMP DO PRIVATE(edge_index,cell_1, cell_2, cartesian_v,cartesian_c,cartesian_center,&
!$OMP circumcenters_vector,x,y,edge_vector,edge_normal_vector,real_tmp,tmp_vector,&
!$OMP lon,lat, vertex_1, vertex_2)
    DO edge_index = 1, no_of_edges
      !------------------------------------------
      vertex_1 = edges%get_vertex_index(edge_index,1)
      vertex_2 = edges%get_vertex_index(edge_index,2)
      cell_1 = edges%get_cell_index(edge_index,1)
      cell_2 = edges%get_cell_index(edge_index,2)
      
      cartesian_v(1) = verts%cartesian(vertex_1)
      cartesian_v(2) = verts%cartesian(vertex_2)
      !------------------------------------------
      IF (cell_1 > 0) THEN
        cartesian_c(1) = cells%cartesian_center(cell_1)
      ENDIF
      IF (cell_2 > 0) THEN
        cartesian_c(2) = cells%cartesian_center(cell_2)
      ENDIF
      !------------------------------------------
      IF ( cell_1 < 1 .and.  cell_2 < 1 ) &
        & CALL finish("grid_sphericalTriangleGridGeometry",".not. associated(cartesian_c both)");
      !------------------------------------------
      ! edge center
      IF (cell_1 > 0 .and. cell_2 > 0 ) THEN
        cartesian_center = inter_section ( &
          & cartesian_c(1), cartesian_c(2), &
          & cartesian_v(1), cartesian_v(2))        
      ELSE
        cartesian_center = sphere_cartesian_midpoint(cartesian_v(1), cartesian_v(2))
        IF (cell_1 <=0) cartesian_c(1) = cartesian_center
        IF (cell_2 <=0) cartesian_c(2) = cartesian_center
      ENDIF
      edges%cartesian_center(edge_index) = cartesian_center
      edges%center(edge_index)     = cc2gc(cartesian_center)
      edges%dual_cartesian_center(edge_index) = &
        & sphere_cartesian_midpoint(cartesian_c(1), cartesian_c(2))
      !------------------------------------------
      ! compute lengths
      edges%primal_edge_length(edge_index) = &
        & arc_length (cartesian_v(1), cartesian_v(2))

      edges%get_edge_vert_length(edge_index,1) = 0.5_wp*edges%primal_edge_length(edge_index)
      edges%get_edge_vert_length(edge_index,2) = 0.5_wp*edges%primal_edge_length(edge_index)

      !------------------------------------------
      IF (cell_1 > 0) THEN
        edges%get_edge_cell_length(edge_index,1) = &
          & arc_length (cartesian_center, cartesian_c(1))
      ELSE
        cartesian_c(1) = cartesian_center 
        edges%get_edge_cell_length(edge_index,1) = 0.0_wp
      ENDIF
      IF ( cell_2 > 0) THEN
        edges%get_edge_cell_length(edge_index,2) = &
          & arc_length (cartesian_center, cartesian_c(2))
      ELSE
        cartesian_c(2) = cartesian_center
        edges%get_edge_cell_length(edge_index,2) = 0.0_wp
      ENDIF
      
      edges%dual_edge_length(edge_index) = &
          & edges%get_edge_cell_length(edge_index,1) + edges%get_edge_cell_length(edge_index,2)

      !------------------------------------------
      circumcenters_vector%x = cartesian_c(2)%x - cartesian_c(1)%x
      circumcenters_vector%x = circumcenters_vector%x /&
       & d_norma_3d(circumcenters_vector)
      ! update quad area (just an approximation, not in use)
      edges%quad_area(edge_index) = 0.5_wp * edges%primal_edge_length(edge_index)  * &
          & edges%dual_edge_length(edge_index)
      !------------------------------------------
      ! compute normals
      ! define the coordinate system tangent on the sphere at the edge center
      ! CALL sphere_tanget_coordinates(edges%center(edge_index),x,y)
      !-------------------------------------------------------------------------
      lon = edges%center(edge_index)%lon
      lat = edges%center(edge_index)%lat
      x%x(1) = -SIN(lon)
      x%x(2) =  COS(lon)
      x%x(3) =  0._wp
      y%x(1) = -COS(lon) * SIN(lat)
      y%x(2) = -SIN(lon) * SIN(lat)
      y%x(3) =  COS(lat)
      !-------------------------------------------------------------------------


      ! the vector between the two vertices
      edge_vector%x = cartesian_v(2)%x - cartesian_v(1)%x
      ! the normal to the edge, as a product of the normal to the sphere
      ! times the edge vector
      edge_normal_vector   = vector_product(edge_vector, cartesian_center)
      edge_normal_vector%x = edge_normal_vector%x / &
       & d_norma_3d(edge_normal_vector)

      !      write(*,*) " : edge_normal_vector=",edge_normal_vector%x, d_norma_3d(edge_normal_vector)
      !      call flush(6)

      ! set the correct orientation from cartesian_c1 to cartesian_c2
      ! ie. from the center of the first element to the center of the second
      real_tmp = DOT_PRODUCT(edge_normal_vector%x,circumcenters_vector%x)
      IF (real_tmp < 0.0_wp) edge_normal_vector%x = -1.0_wp * edge_normal_vector%x

      ! get the components of normal in the local coordinates x,y
      edges%primal_normal(edge_index)%v1 = DOT_PRODUCT(edge_normal_vector%x,x%x)
      edges%primal_normal(edge_index)%v2 = DOT_PRODUCT(edge_normal_vector%x,y%x)

      !      write(*,*) " : primal_normal=",gridEdges%primal_normal(edgeIndex)
      !      call flush(6)
      ! get the dual normal (tangent to the edge)
      ! this is a right handed for dual-primal system (left handed for primal-dual)
      edges%dual_normal(edge_index)%v1 = edges%primal_normal(edge_index)%v2
      edges%dual_normal(edge_index)%v2 = -edges%primal_normal(edge_index)%v1

      !      write(*,*) " : dual_normal=",gridEdges%dual_normal(edgeIndex)
      !      call flush(6)

      ! the product edgeVector * circumcentersVector
      tmp_vector = vector_product(edge_vector,circumcenters_vector)
      edges%system_orientation(edge_index)= &
        & INT(SIGN(1._wp,DOT_PRODUCT(tmp_vector%x,cartesian_center%x)))

      ! this should be the same as the sign of the edge_vector*dual_normal
      tmp_vector = vector_product(circumcenters_vector,cartesian_center)
      IF (edges%system_orientation(edge_index) /= &
        & INT(SIGN(1._wp,DOT_PRODUCT(edge_vector%x,tmp_vector%x)))) THEN
        CALL finish ('grid_sphericalTriangleGridGeometry',&
          & 'gridEdges%system_orientation weird')
      ENDIF

    ENDDO
!$OMP END DO
    ! the following edge geometric properties have been calculated
    !    gridEdges%center
    !    gridEdges%primal_edge_length
    !    gridEdges%dual_edge_length
    !    gridEdges%edge_vert_length
    !    gridEdges%edge_cell_length
    !    gridEdges%primal_normal
    !    gridEdges%dual_normal
    !    gridEdges%system_orientation
    !------------------------------------------
!$OMP END PARALLEL
        
     CALL order_cell_connectivity(in_grid_id)
    
!$OMP PARALLEL
    !------------------------------------------
    ! compute verts%dual_area
!$OMP DO PRIVATE(vertex_index,edge_index,i,cell_1,cell_2,cartesian_c)
    DO vertex_index = 1, no_of_verts
      verts%dual_area(vertex_index) = 0.0_wp 
      DO i=1,verts%max_connectivity
        edge_index = verts%get_edge_index(vertex_index,i)
        IF (edge_index > 0) THEN
          ! add the dual area corresponding to the edge
          cell_1 = edges%get_cell_index(edge_index,1)
          cell_2 = edges%get_cell_index(edge_index,2)
          !------------------------------------------
          IF (cell_1 > 0) THEN
            cartesian_c(1) = cells%cartesian_center(cell_1)
          ELSE
            cartesian_c(1) = edges%cartesian_center(edge_index)            
          ENDIF
          IF (cell_2 > 0) THEN
            cartesian_c(2) = cells%cartesian_center(cell_2)
          ELSE
            cartesian_c(2) = edges%cartesian_center(edge_index)
          ENDIF
          
          verts%dual_area(vertex_index) = verts%dual_area(vertex_index) + &
            triangle_area(cartesian_c(1), verts%cartesian(vertex_index), cartesian_c(2))
       ENDIF
     ENDDO

    ENDDO
!$OMP END DO
    !------------------------------------------
    
    !------------------------------------------
    ! calculate verts%edge_orientation,
    ! this is done in order_cell_connectivity
! !$OMP DO PRIVATE(vertex_index,edge_index,i)
!     DO vertex_index=1,no_of_verts
!       DO i=1,verts%max_connectivity
!         edge_index = verts%get_edge_index(vertex_index,i)
!         IF (edge_index > 0) THEN
!           ! compute orientation
!           IF (edges%get_vertex_index(edge_index,1) == vertex_index) THEN
!             verts%get_edge_orient(vertex_index,i) = edges%system_orientation(edge_index)
!           ELSE
!             verts%get_edge_orient(vertex_index,i) = -edges%system_orientation(edge_index)
!           ENDIF
!        ENDIF
!      ENDDO
!    ENDDO
! !$OMP END DO


    !----------------------------------------------------------
    ! Finally,rescale distances by radius of the sphere
    sphere_radius         = compute_grid%sphere_radius
    sphere_radius_squared = sphere_radius*sphere_radius
!$OMP DO PRIVATE(vertex_index)
    DO vertex_index=1,no_of_verts
      verts%dual_area(vertex_index) = &
        & sphere_radius_squared * verts%dual_area(vertex_index)
    ENDDO
!$OMP END DO
!$OMP DO PRIVATE(cell_index)
    DO cell_index=1,no_of_cells
      cells%area(cell_index) = sphere_radius_squared * cells%area(cell_index)
    ENDDO
!$OMP END DO
!$OMP DO PRIVATE(edge_index,i)
     DO edge_index=1, no_of_edges
      edges%primal_edge_length(edge_index) = sphere_radius * edges%primal_edge_length(edge_index)
      edges%dual_edge_length(edge_index)   = sphere_radius * edges%dual_edge_length(edge_index)
      DO i=1,2
        edges%get_edge_vert_length(edge_index, i) = &
          & sphere_radius * edges%get_edge_vert_length(edge_index, i)
        edges%get_edge_cell_length(edge_index, i) = &
          & sphere_radius * edges%get_edge_cell_length(edge_index, i)
      ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

    !----------------------------------
    ! calculate remaining properties
    CALL get_cell_barycenters(in_grid_id, output_centers=cells%cartesian_barycenter)
    CALL cartesian_to_geographical(cells%cartesian_barycenter, no_of_cells, cells%barycenter)
    !----------------------------------

    !----------------------------------
    CALL timer_stop(timer_set_sphere_geom_grid)
    CALL print_timer(timer_set_sphere_geom_grid)
    CALL delete_timer(timer_set_sphere_geom_grid)
    !----------------------------------

  END SUBROUTINE compute_sphere_grid_geometry
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !   SUBROUTINE reorder_cell_connectivity(in_grid_id)
  !   Reorders the cell vertices edges, neigbors as follows:
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
  !>
  !! Reorders the cell vertices, edges, neigbors in counterclock way.
  !! The geometry of the gird has to be calculated beofre this method is called
  SUBROUTINE order_cell_connectivity(in_grid_id)
    INTEGER, INTENT(in) :: in_grid_id

    TYPE(t_grid),         POINTER :: in_grid
    TYPE(t_grid_cells),   POINTER :: cells
    TYPE(t_grid_edges),   POINTER :: edges
    TYPE(t_grid_vertices),   POINTER :: verts
    INTEGER :: no_of_cells, max_no_of_vertices
    INTEGER :: no_of_verts, max_vert_connectivity

    INTEGER :: vertex_no, cell_idx, vert_idx, no_of_neigbors
    
    INTEGER :: cell_no, edge_1, edge_2,  vertex_1, vertex_2
!     INTEGER :: bnd_edge_1, bnd_edge_2, bnd_cell_1, bnd_cell_2
    INTEGER :: current_edges, ncell, edge_in_cell
    INTEGER, ALLOCATABLE :: cell_edges(:), cell_verts(:), cell_neigbors(:), cell_orientation(:)

    INTEGER :: no_of_first_cells, no_of_vertex_cells, i
    INTEGER, ALLOCATABLE :: first_cell(:),first_edge(:)
    
!     LOGICAL :: polygone_is_around_pole
!     REAL(wp) :: lon_1, lon_2
    
    CHARACTER(len=*), PARAMETER :: method_name = 'order_cell_connectivity'


    !-------------------------------------------------------------------------
!    CALL message(method_name, 'starts ...')
!    write(0,*) 'in_grid_id=',in_grid_id
    in_grid => get_grid(in_grid_id)
    cells   => in_grid%cells
    edges   => in_grid%edges
    verts   => in_grid%verts
    
    no_of_cells = cells%no_of_existcells
    max_no_of_vertices = cells%max_no_of_vertices
    no_of_verts=verts%no_of_existvertices
    max_vert_connectivity=verts%max_connectivity
    
!    write(0,*) 'max_no_of_vertices=',max_no_of_vertices

!$OMP PARALLEL PRIVATE(cell_edges, cell_verts, cell_neigbors, cell_orientation, i, &
!$OMP                  first_cell, first_edge)
    ALLOCATE (cell_edges(max_no_of_vertices),cell_verts(max_no_of_vertices),&
      & cell_neigbors(max_no_of_vertices), cell_orientation(max_no_of_vertices), stat=i)
    IF (i >0) THEN
      CALL finish (method_name, 'Problem in allocating auxiliary arrays')
    ENDIF
      
!     print *, "======================================="

!------------------------------------------------------
!   Reorder cell connectivity
!------------------------------------------------------

!$OMP DO PRIVATE(cell_no, edge_1, edge_2, edge_in_cell, current_edges, &
!$OMP   vertex_1, vertex_2, ncell, i)
    DO cell_no=1,no_of_cells
!      write(0,*) 'cell_no=',cell_no

      cell_edges(:) = 0
      cell_verts(:) = 0
      cell_neigbors(:) = 0
      ! find counterclock-wise edge order
      edge_1 = cells%get_edge_index(cell_no,1)
      vertex_1 = edges%get_vertex_index(edge_1,1)
      vertex_2 = edges%get_vertex_index(edge_1,2)
!       write(0,*) cell_no, cells%get_edge_orient(cell_no,1), edges%system_orientation(edge_1)
      IF (cells%get_edge_orient(cell_no,1) * edges%system_orientation(edge_1) > 0) THEN
        !reverse order
        i = vertex_1
        vertex_1 = vertex_2
        vertex_2 = i
      ENDIF

      ! store the first edge, vert,neigbor
!      write(0,*) 'store the first edge, vert,neigbor=', edge_1,vertex_1
      
      cell_edges(1) = edge_1
      cell_orientation(1) = cells%get_edge_orient(cell_no,1)
      cell_verts(1) = vertex_1
      IF (edges%get_cell_index(edge_1,1) == cell_no) THEN
        cell_neigbors(1) = edges%get_cell_index(edge_1,2)
      ELSE
        cell_neigbors(1) = edges%get_cell_index(edge_1,1)
      ENDIF
      current_edges = 1
      ! store the rest of edges, verts,neigbors
      DO WHILE (vertex_2 /= cell_verts(1))
        CALL get_next_cell_edge_vertex(in_grid_id, cell_no, edge_1, vertex_2, &
          & edge_2, vertex_1, edge_in_cell)

        current_edges = current_edges + 1
        cell_edges(current_edges)       = edge_2
        cell_orientation(current_edges) = cells%get_edge_orient(cell_no,edge_in_cell)
        cell_verts(current_edges)       = vertex_2
        IF (edges%get_cell_index(edge_2,1) == cell_no) THEN
          cell_neigbors(current_edges) = edges%get_cell_index(edge_2,2)
        ELSE
          cell_neigbors(current_edges) = edges%get_cell_index(edge_2,1)
        ENDIF

        ! get the next one
        edge_1   = edge_2
        vertex_2 = vertex_1
      ENDDO
      cells%get_edge_index(cell_no,1:max_no_of_vertices)   = cell_edges(1:max_no_of_vertices)
      cells%get_edge_orient(cell_no,1:max_no_of_vertices)  = &
        & cell_orientation(1:max_no_of_vertices)
      cells%get_vertex_index(cell_no,1:max_no_of_vertices)  = cell_verts(1:max_no_of_vertices)
      cells%get_neighbor_index(cell_no,1:max_no_of_vertices)= &
        & cell_neigbors(1:max_no_of_vertices)

    ENDDO !cell_no=1,no_of_cells
!$OMP END DO

    DEALLOCATE (cell_edges, cell_verts, cell_neigbors, cell_orientation)
! !$OMP BARRIER
    !check edges
!$OMP DO PRIVATE(cell_no, edge_1, edge_2, edge_in_cell, current_edges, &
!$OMP               vertex_1, vertex_2, ncell, i)
    DO cell_no=1,no_of_cells
      DO current_edges=1,max_no_of_vertices
        ncell = cells%get_neighbor_index(cell_no,current_edges)

        IF (ncell < cell_no .AND. ncell /= 0) THEN
          edge_1 = cells%get_edge_index(cell_no,current_edges)
          vertex_1 = cells%get_vertex_index(cell_no,current_edges)
          ! find edge in the neighbor
          DO i=1,max_no_of_vertices
            IF (cells%get_edge_index(ncell,i) == edge_1) &
              EXIT
          ENDDO
          IF (cells%get_vertex_index(ncell,i) == vertex_1) THEN
            CALL print_cell_angles(in_grid_id, cell_no, 'Wrong orientation')
            CALL print_cell_angles(in_grid_id, ncell, 'Wrong orientation')
            CALL finish(method_name, 'Wrong orientation')
          ENDIF
        ENDIF

      ENDDO!check edges

    ENDDO !cell_no=1,no_of_cells
!$OMP END DO

!------------------------------------------------------
!   Reorder vertex connectivity
!------------------------------------------------------

   ALLOCATE (first_cell(max_vert_connectivity),first_edge(max_vert_connectivity),&
      & stat=i)
    IF (i >0) THEN
      CALL finish (method_name, 'Problem in allocating auxiliary arrays')
    ENDIF
! !$OMP BARRIER


!$OMP DO PRIVATE(vertex_no,no_of_first_cells,cell_idx, cell_no,vert_idx,edge_1,&
!$OMP  no_of_neigbors, no_of_vertex_cells,ncell)
    DO vertex_no=1,no_of_verts
      ! first get the start cells
      ! start cells are the ones that have no neigbor in the clock-wise side
      ! if all cells have neigbor, then get the first non-zero cell
      first_edge(:) = 0
      no_of_first_cells=0
      DO cell_idx=1,max_vert_connectivity
        cell_no  = verts%get_cell_index(vertex_no,cell_idx)
        IF (cell_no == 0) CYCLE
        
        ! get the triangle edge containing this vertex counter clockwise
        ! works only for triangles
        vert_idx = 1        
        DO WHILE ( cells%get_vertex_index(cell_no,vert_idx) /= vertex_no)
          vert_idx = vert_idx + 1
        ENDDO
        edge_1 = cells%get_edge_index(cell_no,vert_idx)
        IF (MINVAL(edges%get_cell_index(edge_1,:)) == 0) THEN
          ! add this cell to the first_cell list
          no_of_first_cells = no_of_first_cells + 1
          first_cell(no_of_first_cells) = cell_no
          first_edge(no_of_first_cells) = edge_1
        ENDIF
      ENDDO

      IF (no_of_first_cells == 0) THEN
        cell_idx = 1
        DO WHILE ( verts%get_cell_index(vertex_no,cell_idx) == 0 )
          cell_idx = cell_idx + 1
        ENDDO
        no_of_first_cells=1
        first_cell(no_of_first_cells) = verts%get_cell_index(vertex_no,cell_idx)
      ENDIF


      no_of_neigbors = 0
      no_of_vertex_cells = 0
     ! go through the list of first cells
      DO i=1,no_of_first_cells
        cell_no=first_cell(i)

        IF ( first_edge(i) /= 0 ) THEN
          edge_1 = first_edge(i)
          no_of_neigbors = no_of_neigbors + 1
          ! fill the neigbors
          verts%get_edge_index(vertex_no,no_of_neigbors)     = &
            &  first_edge(i)
          
          IF (edges%get_vertex_index(edge_1,1) == vertex_no) THEN
            verts%get_neighbor_index(vertex_no,no_of_neigbors) = &
              & edges%get_vertex_index(edge_1,2)
          ELSE
            verts%get_neighbor_index(vertex_no,no_of_neigbors) = &
              & edges%get_vertex_index(edge_1,1)
          ENDIF            
        ENDIF !  ( first_edge(i) /= 0 )

        DO WHILE(.true.)
          ! add this cell/edge/vertex
          no_of_neigbors = no_of_neigbors + 1
          no_of_vertex_cells = no_of_vertex_cells + 1
!           write(0,*) vertex_no,no_of_neigbors,no_of_vertex_cells,cell_no

          verts%get_cell_index(vertex_no,no_of_vertex_cells) = cell_no
          
          ! get the vertex before this one
          vert_idx = 1
          DO WHILE ( cells%get_vertex_index(cell_no,vert_idx) /= vertex_no)
            vert_idx = vert_idx + 1
          ENDDO
          IF (vert_idx == 1) THEN
            vert_idx = max_no_of_vertices
          ELSE
            vert_idx = vert_idx-1
          ENDIF

          ! fill the neigbors
          edge_1=cells%get_edge_index(cell_no,vert_idx)
          verts%get_neighbor_index(vertex_no,no_of_neigbors) = &
            & cells%get_vertex_index(cell_no,vert_idx)
          verts%get_edge_index(vertex_no,no_of_neigbors)     = &
            & edge_1

          ! get the next cell/edge/vertex
          ncell = edges%get_cell_index(edge_1, 2)
          IF (ncell == cell_no) &
            & ncell = edges%get_cell_index(edge_1, 1)
          IF (ncell == first_cell(i)) EXIT ! we are done
          IF (ncell == 0) EXIT ! we are done

          cell_no  = ncell

        ENDDO !WHILE(.true.)

      ENDDO ! i=1,no_of_first_cells

      !check if polygone_is_around_pole
!       write(0,*) "z:", verts%cartesian(vertex_no)%x(3)
!       polygone_is_around_pole = .true.
!       lon_1 = cells%center(verts%get_cell_index(vertex_no,1))%lon
!       DO vert_idx = 2, no_of_neigbors
!         lon_2 = cells%center(verts%get_cell_index(vertex_no,vert_idx))%lon
!         IF (lon_1 * lon_2 < 0.0_wp) THEN
!            IF (lon_1 < 0.0_wp) THEN
!              write(0,*) "lon_1,lon_2 < 0", lon_1,lon_2
!              polygone_is_around_pole = .false.
!              EXIT
!            ENDIF
!         ELSE
!           IF (lon_1 >= lon_2) THEN
!             write(0,*) "lon_1 >= lon_2 ", lon_1,lon_2
!             polygone_is_around_pole = .false.
!             EXIT
!           ENDIF
!         ENDIF
!         lon_1 = lon_2  
!       ENDDO
!       polygone_is_around_pole = .false.
!       IF (polygone_is_around_pole) THEN
!         write(0,*) "polygone_is_around_pole:", verts%cartesian(vertex_no)%x(3)
!         write(0,*) "cells:", verts%get_cell_index(vertex_no,  :)
!         ! reverse order
!         DO vertex_1=1,no_of_neigbors / 2
!           
!           vertex_2 = no_of_neigbors - vertex_1 + 1
!           print *, "exhange:",vertex_1,  vertex_2
!           
!           vert_idx = verts%get_neighbor_index(vertex_no,vertex_1)
!           edge_1   = verts%get_edge_index(vertex_no,    vertex_1)
!           cell_no  = verts%get_cell_index(vertex_no,    vertex_1)
! 
!           verts%get_neighbor_index(vertex_no,vertex_1) = &
!             verts%get_neighbor_index(vertex_no,vertex_2)
!           verts%get_edge_index(vertex_no,    vertex_1) = &
!             verts%get_edge_index(vertex_no,    vertex_2) 
!           verts%get_cell_index(vertex_no,    vertex_1) = &
!             verts%get_cell_index(vertex_no,    vertex_2)
!             
!           verts%get_neighbor_index(vertex_no,vertex_2) = vert_idx
!           verts%get_edge_index(vertex_no,    vertex_2) = edge_1
!           verts%get_cell_index(vertex_no,    vertex_2) = cell_no
!         ENDDO
!          write(0,*) "cells:", verts%get_cell_index(vertex_no,  :)
!        ENDIF

      ! clean the rest of the indexes
      verts%no_of_neigbors(vertex_no) = no_of_neigbors
      DO vert_idx = no_of_neigbors + 1, max_vert_connectivity
        verts%get_edge_index(vertex_no,vert_idx)     = 0
        verts%get_neighbor_index(vertex_no,vert_idx) = 0
      ENDDO
      DO vert_idx = no_of_vertex_cells + 1, max_vert_connectivity
        verts%get_cell_index(vertex_no,vert_idx)     = 0
      ENDDO


    ! calculate verts%edge_orientation
      DO i=1,no_of_neigbors
        edge_1 = verts%get_edge_index(vertex_no,i)
        ! compute orientation
        IF (edges%get_vertex_index(edge_1,1) == vertex_no) THEN
          verts%get_edge_orient(vertex_no,i) = edges%system_orientation(edge_1)
        ELSE
          verts%get_edge_orient(vertex_no,i) = -edges%system_orientation(edge_1)
        ENDIF
      ENDDO

    ENDDO ! vertex_no=1,no_of_verts
!$OMP END DO


!$OMP END PARALLEL
!    CALL message(method_name, 'ended')
!     print *, "======================================="
    RETURN

  END SUBROUTINE order_cell_connectivity
  !-------------------------------------------------------------------------
    

  !-------------------------------------------------------------------------
  SUBROUTINE get_celledge_vertices(in_grid_id, edge1, edge2, vertex1, common_vertex, vertex2)
    INTEGER, INTENT(in) :: in_grid_id, edge1, edge2
    INTEGER, INTENT(out):: vertex1, common_vertex, vertex2

    TYPE(t_grid_edges), POINTER :: edges
    INTEGER :: p1, p2

    edges => get_edges(in_grid_id)

    common_vertex = get_common_edge_vertex(edges, edge1, edge2, p1, p2)

    vertex1 = edges%get_vertex_index(edge1,1)
    IF (vertex1 == common_vertex) &
      & vertex1 =  edges%get_vertex_index(edge1,2)

    vertex2 = edges%get_vertex_index(edge2,1)
    IF (vertex2 == common_vertex) &
      & vertex2 =  edges%get_vertex_index(edge2,2)

    RETURN
  END SUBROUTINE get_celledge_vertices
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !   INTEGER FUNCTION get_common_edge_vertex(edges, edge_index1, edge_index2 p1, p2)
  !>
  !! Returns the index of the common vertex of the two edges. Private
  INTEGER FUNCTION get_common_edge_vertex(edges, edge_index1, edge_index2, p1, p2)
    !  TYPE(grid_cells), POINTER, INTENT(IN)  :: cells
    TYPE(t_grid_edges), POINTER :: edges
    INTEGER, INTENT(in) :: edge_index1, edge_index2
    INTEGER, INTENT(out) :: p1 ,p2

    p1 = 1
    p2 = 1
    get_common_edge_vertex = edges%get_vertex_index(edge_index1,1)
    IF (get_common_edge_vertex == edges%get_vertex_index(edge_index2,1)) RETURN
    p2 = 2
    IF (get_common_edge_vertex == edges%get_vertex_index(edge_index2,2)) RETURN
    p1 = 2
    get_common_edge_vertex = edges%get_vertex_index(edge_index1,2)
    IF (get_common_edge_vertex == edges%get_vertex_index(edge_index2,2)) RETURN
    p2 = 1
    IF (get_common_edge_vertex == edges%get_vertex_index(edge_index2,1)) RETURN

    get_common_edge_vertex = -1
    CALL finish ('getCommonEdgeVertex', 'No common vertex was found!')

    RETURN
  END FUNCTION get_common_edge_vertex
  !-------------------------------------------------------------------------

 !-------------------------------------------------------------------------
  FUNCTION edges_normal_angle(in_grid_id_1, edge_1, in_grid_id_2, edge_2) RESULT(angle)
    INTEGER, INTENT(in) :: in_grid_id_1, edge_1, in_grid_id_2, edge_2
    REAL(wp)  :: angle

    TYPE(t_grid_edges), POINTER :: edges_1, edges_2
!    TYPE(cartesian_coordinates) :: n1, n2
    REAL(wp)  :: prod

    edges_1 => get_edges(in_grid_id_1)
    edges_2 => get_edges(in_grid_id_2)

    prod = edges_1%primal_normal(edge_1)%v1 * edges_2%primal_normal(edge_2)%v1 + &
      &    edges_1%primal_normal(edge_1)%v2 * edges_2%primal_normal(edge_2)%v2
     angle = ACOS(prod) * rad2deg
    RETURN
    !-----------------------------------------------------------
    ! check if the angle between parent,child is ~ 90
!     CALL gvec2cvec(edges_1%primal_normal(edge_1)%v1,&
!                     edges_1%primal_normal(edge_1)%v2, &
!                     edges_1%center(edge_1)%lon, &
!                     edges_1%center(edge_1)%lat, &
!                     n1%x(1),n1%x(2),n1%x(3))
!     CALL gvec2cvec(edges_2%primal_normal(edge_2)%v1,&
!                     edges_2%primal_normal(edge_2)%v2, &
!                     edges_2%center(edge_2)%lon, &
!                     edges_2%center(edge_2)%lat, &
!                     n2%x(1),n2%x(2),n2%x(3))
!
!     angle = angle_of_vectors(n1,n2) * rad2deg
!
!     RETURN

  END FUNCTION edges_normal_angle
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  FUNCTION edges_3D_angle(in_grid_id_1, edge_1, in_grid_id_2, edge_2) RESULT(angle)
    INTEGER, INTENT(in) :: in_grid_id_1, edge_1, in_grid_id_2, edge_2
    REAL(wp)  :: angle

    TYPE(t_grid), POINTER :: in_grid
    TYPE(t_grid_edges), POINTER :: edges
    TYPE(t_grid_vertices), POINTER :: verts
    TYPE(t_cartesian_coordinates) :: v1, v2

    in_grid => get_grid(in_grid_id_1)
    edges => in_grid%edges
    verts => in_grid%verts
    v1%x = verts%cartesian(edges%get_vertex_index(edge_1,2))%x - &
         &  verts%cartesian(edges%get_vertex_index(edge_1,1))%x

    in_grid => get_grid(in_grid_id_2)
    edges => in_grid%edges
    verts => in_grid%verts
    v2%x = verts%cartesian(edges%get_vertex_index(edge_2,2))%x - &
         & verts%cartesian(edges%get_vertex_index(edge_2,1))%x

    angle = angle_of_vectors(v1,v2)

    RETURN

  END FUNCTION edges_3D_angle
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  FUNCTION edges_cell_angle(in_grid_id, cell_no, edge1_in_cell) RESULT(angle)
    INTEGER, INTENT(in) :: in_grid_id, cell_no, edge1_in_cell
    REAL(wp)  :: angle

    TYPE(t_grid), POINTER :: in_grid
    TYPE(t_grid_cells), POINTER :: cells
    ! TYPE(grid_edges), POINTER :: edges
    TYPE(t_grid_vertices), POINTER :: verts
    TYPE(t_cartesian_coordinates) :: v1, v2
    INTEGER :: edge2_in_cell, edge1, edge2
    INTEGER :: vertex1, vertex2, common_vertex

    angle = NO_ANGLE
    in_grid => get_grid(in_grid_id)
    cells => in_grid%cells
    verts => in_grid%verts

    IF (edge1_in_cell == cells%max_no_of_vertices) THEN
     edge2_in_cell = 1
    ELSE
      edge2_in_cell = edge1_in_cell + 1
    ENDIF
    edge1 = cells%get_edge_index(cell_no,edge1_in_cell)
    edge2 = cells%get_edge_index(cell_no,edge2_in_cell)
    IF (edge1 == 0 .OR. edge2 == 0) RETURN

    CALL get_celledge_vertices(in_grid_id, edge1, edge2, vertex1, common_vertex, vertex2)

    v1%x = verts%cartesian(vertex1)%x - &
         &  verts%cartesian(common_vertex)%x
    v2%x = verts%cartesian(vertex2)%x - &
         &  verts%cartesian(common_vertex)%x

    angle = angle_of_vectors(v1,v2)

    RETURN

  END FUNCTION edges_cell_angle
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE print_cell_angles(in_grid_id, cell_no, opt_text)
    INTEGER, INTENT(in) :: in_grid_id, cell_no
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: opt_text

    TYPE(t_grid), POINTER :: in_grid
    TYPE(t_grid_cells), POINTER :: cell

    INTEGER :: j
    REAL(wp):: angle

    CALL print_grid_cell(in_grid_id, cell_no, opt_text)

    in_grid => get_grid(in_grid_id)
    cell => in_grid%cells

    print*, '======================================================'
    print*, 'angles of cell ', cell_no, ":"
    DO j=1,cell%max_no_of_vertices
      angle = edges_cell_angle(in_grid_id, cell_no, j)
      IF (angle /= NO_ANGLE) THEN
        print*, '     ', j, ' angle:', angle*rad2deg
      ELSE
        print*, '     ', j,' angle: error'
      ENDIF
    ENDDO
    print*, '======================================================'

  END SUBROUTINE print_cell_angles
  !-------------------------------------------------------------------------

  !---------------------------------------------------------------------------
  SUBROUTINE get_triangle_circumcenters(in_grid_id, output_centers)
    INTEGER, INTENT(in) :: in_grid_id
    TYPE(t_cartesian_coordinates), POINTER, OPTIONAL :: output_centers(:)

    TYPE(t_grid), POINTER :: in_grid
    TYPE(t_grid_cells), POINTER :: cells
    TYPE(t_grid_vertices), POINTER :: verts
    TYPE(t_cartesian_coordinates), POINTER :: centers(:)

    INTEGER :: no_of_cells,max_no_of_vertices,cell_index,i

    in_grid => get_grid(in_grid_id)
    cells => in_grid%cells
    verts => in_grid%verts
    no_of_cells = cells%no_of_existcells
    max_no_of_vertices = cells%max_no_of_vertices

    IF (PRESENT(output_centers)) THEN
      ALLOCATE(output_centers(no_of_cells), stat=i)
      IF (i > 0) &
        CALL finish ('set_triangle_circumcenter', 'ALLOCATE')
      centers => output_centers
    ELSE
      centers => cells%cartesian_center
    ENDIF

    IF (max_no_of_vertices /= 3) &
      & CALL finish('set_triangle_circumcenter', 'max_no_of_vertices /= 3')

!$OMP PARALLEL 
!$OMP DO PRIVATE(cell_index)
    DO cell_index = 1,no_of_cells
      centers(cell_index) = &
        & circum_center(verts%cartesian(cells%get_vertex_index(cell_index,1)), &
        &               verts%cartesian(cells%get_vertex_index(cell_index,2)), &
        &               verts%cartesian(cells%get_vertex_index(cell_index,3)) )
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE get_triangle_circumcenters
  !---------------------------------------------------------------------------

  !---------------------------------------------------------------------------
  SUBROUTINE get_cell_barycenters(in_grid_id, center_type, output_centers)
    INTEGER, INTENT(in) :: in_grid_id
    INTEGER, INTENT(in), OPTIONAL  :: center_type
    TYPE(t_cartesian_coordinates), POINTER, OPTIONAL :: output_centers(:)

    INTEGER :: use_center_type
    TYPE(t_grid), POINTER :: in_grid
    TYPE(t_grid_cells), POINTER :: cells
    TYPE(t_grid_edges), POINTER :: edges
    TYPE(t_grid_vertices), POINTER :: verts

    TYPE(t_cartesian_coordinates) :: v(32),c_center,t_center
    TYPE(t_cartesian_coordinates), POINTER :: centers(:)
    INTEGER :: no_of_cells,max_no_of_vertices,cell_no
    INTEGER :: k,l,cell_edges,edge
    REAL(wp):: area

    in_grid => get_grid(in_grid_id)
    cells => in_grid%cells
    edges => in_grid%edges
    verts => in_grid%verts
    no_of_cells = cells%no_of_existcells
    max_no_of_vertices = cells%max_no_of_vertices

    IF (PRESENT(center_type)) THEN
       use_center_type = center_type
    ELSE
       use_center_type =  use_barycenters
    ENDIF

    IF (PRESENT(output_centers)) THEN
      IF (.NOT. ASSOCIATED(output_centers)) THEN
        ALLOCATE(output_centers(no_of_cells), stat=k)
        IF (k > 0) &
          CALL finish ('set_cell_barycenters', 'ALLOCATE')
        ENDIF
      centers => output_centers
    ELSE
      centers => cells%cartesian_center
    ENDIF

    IF (max_no_of_vertices == 3) THEN
      ! triangle barycenters are easy to compute
!$OMP PARALLEL 
!$OMP DO PRIVATE(cell_no,v)
      DO cell_no = 1,no_of_cells
        v(1)%x = verts%cartesian(cells%get_vertex_index(cell_no,1))%x
        v(2)%x = verts%cartesian(cells%get_vertex_index(cell_no,2))%x
        v(3)%x = verts%cartesian(cells%get_vertex_index(cell_no,3))%x
        !centers(cell_no)%x = (v(1)%x + v(2)%x + v(3)%x)
        !d_normalize(centers(cell_no))
        centers(cell_no) = triangle_norm_barycenter(v(1),v(2),v(3))
      ENDDO
!$OMP END DO
!$OMP END PARALLEL
      RETURN
    ENDIF

    ! general case
!$OMP PARALLEL 
!$OMP DO PRIVATE(cell_no, cell_edges, l, k, edge, v, c_center, t_center, area)
    DO cell_no = 1,no_of_cells
      ! get the edge vertices
      cell_edges = 0
      l = 1
      DO k=1,max_no_of_vertices
        edge = cells%get_edge_index(cell_no,k)
        IF (edge > 0) THEN
          ! add the vertices to the list
          cell_edges = cell_edges+1
          v(l)%x   = verts%cartesian(edges%get_vertex_index(edge,1))%x
          v(l+1)%x = verts%cartesian(edges%get_vertex_index(edge,2))%x
!           WRITE(0,*) edge, edges%get_vertex_index(edge,1), l, 'v(l)=', v(l)
!           WRITE(0,*) edge, edges%get_vertex_index(edge,2), l+1, 'v(l+1)=',v(l+1)
          l = l+2
        ENDIF
      ENDDO

      ! Get the cartesian center of the cell
      c_center = cartesian_norm_center(v,cell_edges*2)
      IF (use_center_type ==  use_cartesian_centers) THEN
        centers(cell_no) = c_center
        CYCLE
      ENDIF

      ! calculate the weighted center of the triangles between edges and c_center
      l = 1
!      total_area = 0.0_wp
      DO k=1,cell_edges
        !t_center = triangle_norm_barycenter(v(l), v(l+1), c_center)
        t_center%x = (v(l)%x + v(l+1)%x + c_center%x)
        area     = triangle_area(v(l), v(l+1), c_center)
        v(k)%x = t_center%x * area
        ! total_area = total_area + area
        l = l+2
      ENDDO
      ! finally get barycenter
      centers(cell_no) = cartesian_norm_center(v,cell_edges)
      ! cells%cartesian_center(cell_no)%x = cells%cartesian_center(cell_no)%x / total_area

      ! print *, cells%cartesian_center(cell_no)%x

    ENDDO ! cell_no = 1,no_of_cells
!$OMP END DO
!$OMP END PARALLEL


  END SUBROUTINE get_cell_barycenters
  !---------------------------------------------------------------------------

  !---------------------------------------------------------------------------
  ELEMENTAL FUNCTION triangle_norm_barycenter(v1,v2,v3) result(center)

    TYPE(t_cartesian_coordinates), INTENT(in) :: v1,v2,v3
    TYPE(t_cartesian_coordinates) :: center

    center%x = (v1%x + v2%x + v3%x)
    d_normalize(center)
!     write(0,*) 'triangle_norm_barycenter 1:', v1
!     write(0,*) 'triangle_norm_barycenter 2:', v2
!     write(0,*) 'triangle_norm_barycenter 3:', v3
!     write(0,*) 'triangle_norm_barycenter center:', center
  END FUNCTION triangle_norm_barycenter
  !-------------------------------------------------------------------------

  !---------------------------------------------------------------------------
  PURE FUNCTION cartesian_norm_center(v,v_size) result(center)

    TYPE(t_cartesian_coordinates), INTENT(in) :: v(:)
    INTEGER, INTENT(in) :: v_size
    TYPE(t_cartesian_coordinates) :: center

    INTEGER :: i

    center%x = 0.0_wp
    DO i=1,v_size
      center%x = center%x + v(i)%x
    ENDDO
    !center%x = center%x / REAL(v_size,wp)
    d_normalize(center)

  END FUNCTION cartesian_norm_center
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
END MODULE mo_local_grid_geometry

