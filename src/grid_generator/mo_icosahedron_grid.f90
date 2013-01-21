! ! #ifdef __xlC__
! ! @PROCESS NOOPTIMIZE
! ! #endif
!>
!! Create spherical grids based on the icosahedron
!!
!! @par Revision History
!! Initial Release by Leonidas Linardakis (2011-01-10)
!!
!! @par Copyright
!! 2002-2008 by DWD and MPI-M
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
MODULE mo_icosahedron_grid
#include "grid_definitions.inc"

  !-------------------------------------------------------------------------
  USE mo_kind,           ONLY: wp
  USE mo_exception,      ONLY: message_text, message, finish
  USE mo_io_units,       ONLY: nnml, filename_max
  USE mo_namelist,       ONLY: position_nml, open_nml, positioned
!   USE mo_physical_constants, ONLY: re

  USE mo_math_utilities, ONLY: t_cartesian_coordinates!, cc2gc
  USE mo_math_constants, ONLY: pi_5
  USE mo_io_local_grid,  ONLY: read_new_netcdf_grid, write_netcdf_grid
  USE mo_grid_toolbox,   ONLY: get_basic_dual_grid
  USE mo_local_grid_geometry,     ONLY: compute_sphere_grid_geometry, get_cell_barycenters!, &
!    & use_cartesian_centers
  USE mo_local_grid_refinement,   ONLY: refine_grid_edgebisection, refine_grid_insert_centers, &
    & complete_grid_connectivity
  USE mo_local_grid_optimization, ONLY: read_grid_optimization_param, optimize_grid
  USE mo_local_grid
!  USE mo_grid_checktools,         ONLY: check_grid
  USE mo_grid_decomposition, ONLY: decompose_all_cells, grid_cluster_subdomains, &
    & get_no_of_domains, get_max_subdomain_cells, write_ascii_grid_decomposition

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: create_icon_grid

  CHARACTER(LEN=*), PARAMETER :: version = '$Id$'

  INTEGER :: no_of_levels, start_level, refinement_method, start_optimize, end_optimize
  INTEGER :: dual_decomposition_domains = -2
  LOGICAL :: use_clustered_decompositions
!   INTEGER :: decompose_roundrobin_at_level = -2
  
  REAL(wp) :: icon_norm_edge, icon_norm_dual_edge
  CHARACTER(LEN=filename_max) :: input_file, output_file, decomposition_ascii_ext
  CHARACTER(LEN=filename_max) :: optimization_extension
  
  INTEGER, PARAMETER :: edge_bisection = 1    ! for refinement_method
  INTEGER, PARAMETER :: dual_cell_centers = 2 ! for refinement_method

CONTAINS

  !-------------------------------------------------------------------------
  !! SUBROUTINE read_grid_optimization_param
  !>
  !! Reads the parameters for local grid optimization
  SUBROUTINE read_icosahedron_grid_param(param_file_name)

    CHARACTER(LEN=*), INTENT(in) :: param_file_name

    INTEGER :: i_status

    NAMELIST /icosahedron_grid/ no_of_levels, start_level, &
      & refinement_method, input_file, output_file, optimization_extension, &
      & start_optimize, end_optimize, use_clustered_decompositions, &
      & dual_decomposition_domains, decomposition_ascii_ext
!       & decompose_roundrobin_at_level

    ! set default values
    no_of_levels = 1
    refinement_method = edge_bisection
    start_level = 0
    output_file = 'test'
    input_file = 'NULL'
    decomposition_ascii_ext = "_cell_domain_ids"
    start_optimize = 1
    end_optimize = -1
    dual_decomposition_domains = -2
    use_clustered_decompositions = .true.
    
    ! read namelist
    CALL open_nml(param_file_name)
    CALL position_nml('icosahedron_grid',STATUS=i_status)
    IF (i_status == positioned) THEN
      READ (nnml,icosahedron_grid)
    ELSE
      WRITE(message_text,'(a)') " File", param_file_name, " not POSITIONED"
      CALL finish ('read_icosahedron_grid_param', message_text)
    ENDIF
    CLOSE(nnml)

    WRITE(message_text,'(a)') "===================================="
    CALL message ('', TRIM(message_text))
    CALL message ('', '   icosahedron_grid')
    WRITE(message_text,'(a,i9)')   'no_of_levels=', no_of_levels
    CALL message ('', TRIM(message_text))
    SELECT CASE(refinement_method)
      CASE(edge_bisection)
        CALL message ('', 'refinement_method = edge_bisection')
      CASE(dual_cell_centers)
        CALL message ('', 'refinement_method = dual_cell_centers insertion')
      CASE default
        CALL finish('read_icosahedron_grid_param','Unkown refinement_method')
    END SELECT
    WRITE(message_text,'(a,a)')    'input_file=', TRIM(input_file)
    CALL message ('', TRIM(message_text))
    WRITE(message_text,'(a,a)')    'output_file=', TRIM(output_file)
    CALL message ('', TRIM(message_text))

    IF (end_optimize < 0) end_optimize=no_of_levels

    CALL read_grid_optimization_param(param_file_name)

  END SUBROUTINE read_icosahedron_grid_param
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE create_icon_grid(param_file_name)
    CHARACTER(LEN=*), INTENT(in) :: param_file_name

    INTEGER :: base_grid_id, next_grid_id, level, end_level
    INTEGER :: no_of_domains
    LOGICAL :: is_dual_decomposed
    CHARACTER(LEN=filename_max) :: file_name

    is_dual_decomposed = .false.
    
    CALL read_icosahedron_grid_param(param_file_name)

    IF (input_file == 'NULL') THEN
      !level = -1, iccosahedron
      base_grid_id = get_icosahedron_grid()
      CALL set_default_geometry_parameters(to_grid_id=base_grid_id, &
        & param_file_name=param_file_name)
      CALL compute_sphere_grid_geometry(base_grid_id)
      CALL set_grid_level(base_grid_id, -1)
      CALL set_grid_parent_id(base_grid_id, 0)
!       IF (-1 == decompose_cells_at_level) THEN
!         CALL decompose_all_cells(base_grid_id, 1)
!       ENDIF
      WRITE(file_name,'(a,a)')  TRIM(output_file), '_icosahedron.nc'
      CALL write_netcdf_grid(base_grid_id, file_name)
      !  print *, ' icon_norm_edge=', icon_norm_edge, ' icon_norm_dual_edge=', icon_norm_dual_edge
      start_level = 0
      end_level = no_of_levels
    ELSE
      base_grid_id = read_new_netcdf_grid(input_file)
      start_level = get_grid_level(base_grid_id) + 1
      end_level = start_level + no_of_levels - 1
   ENDIF

   DO level=start_level,end_level

      !get next level by edge bisection
      WRITE(message_text,'(a,i2.2)') "level=",level
      CALL message ('create_icon_grid', TRIM(message_text))
      
      SELECT CASE(refinement_method)

        CASE(edge_bisection)

          IF (get_number_of_vertices(base_grid_id) == dual_decomposition_domains) THEN
            is_dual_decomposed = .true.
            next_grid_id = get_basic_dual_grid(base_grid_id)
            CALL delete_grid(base_grid_id)
            CALL decompose_all_cells(next_grid_id, 1)
!             CALL get_cell_barycenters(next_grid_id)
            base_grid_id = refine_grid_insert_centers(next_grid_id)
            CALL delete_grid(next_grid_id)
            
          ELSE
            next_grid_id = refine_grid_edgebisection(base_grid_id)
            CALL delete_grid(base_grid_id)
            base_grid_id = next_grid_id
          ENDIF

        CASE(dual_cell_centers)
!          CALL get_cell_barycenters(base_grid_id)
          next_grid_id = get_basic_dual_grid(base_grid_id)
          CALL delete_grid(base_grid_id)
          ! CALL get_cell_barycenters(next_grid_id, use_cartesian_centers)
!          CALL get_cell_barycenters(next_grid_id)
          base_grid_id = refine_grid_insert_centers(next_grid_id)
          CALL delete_grid(next_grid_id)

        CASE default
          CALL finish('create_icon_grid','Unkown refinement_method')
      END SELECT


      !icon_norm_edge = icon_norm_edge * 0.5114_wp
      !icon_norm_dual_edge = icon_norm_dual_edge * 0.5_wp
      !print *, ' icon_norm_edge=', icon_norm_edge, ' icon_norm_dual_edge=', icon_norm_dual_edge
      ! CALL optimize_grid_set_lengths(icon_norm_edge, icon_norm_dual_edge)

      ! apply otimizations
      IF (level >= start_optimize .AND. level <= end_optimize) THEN
        CALL optimize_grid(base_grid_id)
      ENDIF

      CALL compute_sphere_grid_geometry(base_grid_id)
      CALL set_grid_parent_id(base_grid_id, 0)
            
      IF (is_dual_decomposed) THEN
        IF ( use_clustered_decompositions) &
          CALL grid_cluster_subdomains(base_grid_id,1, 1)
       
        no_of_domains = get_no_of_domains(base_grid_id, 1)
!         WRITE(file_name,'(a,i2.2,a,a,i4.4,a)')  TRIM(output_file), level,  &
!           TRIM(optimization_extension), "_hex_dd_", no_of_domains, &
!           TRIM(decomposition_ascii_ext)
        WRITE(message_text,'(i6,a)') get_max_subdomain_cells(base_grid_id, 1), &
          TRIM(decomposition_ascii_ext)
        message_text = ADJUSTL(message_text)
        WRITE(file_name,'(a,i2.2,a,".",a)')  TRIM(output_file), level,  &
          TRIM(optimization_extension), TRIM(message_text)
        CALL write_ascii_grid_decomposition(base_grid_id, 1, file_name)
      ENDIF
      
      WRITE(file_name,'(a,i2.2,a,a)')  TRIM(output_file), level,  &
        TRIM(optimization_extension), ".nc"
      CALL write_netcdf_grid(base_grid_id, file_name)

      ! if decomposition takes place write the ascci decomposition file
      !  if this level >= decompose level
      
!       IF (decompose_cells_at_level > -2 .AND. &
!         & level >= decompose_cells_at_level) THEN        
!         no_of_domains = get_no_of_domains(base_grid_id, 1)
!         WRITE(file_name,'(a,i2.2,a,a,i4.4,a)')  TRIM(output_file), level,  &
!           TRIM(optimization_extension), "_tri_dd_", no_of_domains, &
!           TRIM(decomposition_ascii_ext)
!         CALL write_ascii_grid_decomposition(base_grid_id, 1, file_name)
!       ENDIF
      
    ENDDO ! level=1,no_of_levels

    CALL delete_grid(base_grid_id)

    IF (dual_decomposition_domains > 0 .AND. .NOT. is_dual_decomposed) THEN
      CALL finish("create_icon_grid","grid is NOT dual_decomposed")
    ENDIF

  END SUBROUTINE create_icon_grid
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  INTEGER FUNCTION get_icosahedron_grid() result(out_grid_id)
    ! !DESCRIPTION:
    ! Calculates the Cartesian coordinates of the gridpoints of the
    ! icosahedral grid on the unit sphere. The vertices are saved
    ! in an array from 1 till 12 for comparison with other published
    ! solutions. Computes the  angles associated with the base icosahedron.
    ! Calculate the latitudes of the two times 5 vertices
    ! (the other two important vertices are south and north pole).
    !
    ! !REVISION HISTORY:
    !  Initial Version by Luis Kornblueh, MPI-M, Hamburg, March 2004
    !  Modified by L. Linardakis, MPI-M, Hamburg, August 2010

    TYPE(t_cartesian_coordinates) :: base_vertex(12)

    REAL(wp) :: z_w, z_cosw, z_sinw, z_rlon
    INTEGER :: i_mdist(2:11) ! distribute vertices on latitude ring
    INTEGER :: i_msgn, j


    TYPE(t_grid), POINTER :: out_grid
    TYPE(t_grid_cells), POINTER :: new_cells
    TYPE(t_grid_edges), POINTER :: new_edges
    TYPE(t_grid_vertices), POINTER :: new_verts
    INTEGER :: no_of_output_cells, no_of_output_edges, no_of_output_verts
    INTEGER :: max_cell_vertices, max_vertex_connect
    INTEGER :: vertex, edge, cell, vertex1, vertex2, max_edge
    INTEGER :: vertex_edge(12,12)


    CALL message('==========================','')
    CALL message('Creating icosahedron grid','...')
    !-------------------------------------------------------
    ! first define the vertices of the icosahedron
    z_w      = 2.0_wp*acos(1.0_wp/(2.0_wp*sin(pi_5)))
    z_cosw   = COS(z_w)
    z_sinw   = SIN(z_w)

    ! Northern hemisphere are the first 6 elements of base_vertex(1:6)
    ! Southern hemisphere are the other 6 elements of base_vertex(7:12)
    !
    ! set poles first - it is simple
    base_vertex( 1) = t_cartesian_coordinates( (/ 0.0_wp, 0.0_wp,  1.0_wp /) )
    base_vertex(12) = t_cartesian_coordinates( (/ 0.0_wp, 0.0_wp, -1.0_wp /) )

    ! now set the vertices on the two latitude rings
    DO j = 1, 10
      IF (MOD(j,2) == 0) THEN
        i_mdist(j/2+6)     = -1+(j-1)-10*((j-1)/7)
      ELSE
        i_mdist((j+1)/2+1) = -1+(j-1)-10*((j-1)/7)
      END IF
    END DO

    DO j = 2,11
      ! toggle the hemisphere
      IF (j > 6) THEN
        i_msgn = -1    ! southern
      ELSE
        i_msgn =  1    ! northern
      ENDIF

      ! compute the meridian angle for the base vertex.
      z_rlon =  (1.0_wp+REAL(i_mdist(j),wp))*pi_5

      ! now initialize the coordinates
      base_vertex(j) = t_cartesian_coordinates( &
        & (/ z_sinw*cos(z_rlon), z_sinw*sin(z_rlon), z_cosw*REAL(i_msgn,wp) /) &
        & )
    END DO

    ! the vertices of the icosahedron are defined
    !-------------------------------------------------------
    ! create the grid
    CALL message('Creating icosahedron grid','create the grid...')
    out_grid_id = new_grid()
    out_grid => get_grid(out_grid_id)
    no_of_output_cells = 20
    no_of_output_edges = 30
    no_of_output_verts = 12
    max_cell_vertices  = 3
    max_vertex_connect = 6

    !  allocate grid
    out_grid%ncells = no_of_output_cells
    out_grid%nedges = no_of_output_edges
    out_grid%nverts = no_of_output_verts
    out_grid%cells%max_no_of_vertices = max_cell_vertices
    out_grid%verts%max_connectivity   = max_vertex_connect
    CALL allocate_grid_object(out_grid_id)
    CALL grid_set_exist_eq_allocated(out_grid_id)

    ! fill out_grid
    out_grid%level = 0
    new_cells=>out_grid%cells
    new_edges=>out_grid%edges
    new_verts=>out_grid%verts

    !----------------------------------------------
    ! fill connectivity

    CALL message('Creating icosahedron grid','fill cell vertices...')
    ! fill cell vertices
    new_cells%get_vertex_index(1,1) = 1
    new_cells%get_vertex_index(1,2) = 2
    new_cells%get_vertex_index(1,3) = 3

    new_cells%get_vertex_index(2,1) = 1
    new_cells%get_vertex_index(2,2) = 3
    new_cells%get_vertex_index(2,3) = 4

    new_cells%get_vertex_index(3,1) = 1
    new_cells%get_vertex_index(3,2) = 4
    new_cells%get_vertex_index(3,3) = 5

    new_cells%get_vertex_index(4,1) = 1
    new_cells%get_vertex_index(4,2) = 5
    new_cells%get_vertex_index(4,3) = 6

    new_cells%get_vertex_index(5,1) = 1
    new_cells%get_vertex_index(5,2) = 6
    new_cells%get_vertex_index(5,3) = 2

    new_cells%get_vertex_index(6,1) = 7
    new_cells%get_vertex_index(6,2) = 3
    new_cells%get_vertex_index(6,3) = 2

    new_cells%get_vertex_index(7,1) = 8
    new_cells%get_vertex_index(7,2) = 4
    new_cells%get_vertex_index(7,3) = 3

    new_cells%get_vertex_index(8,1) = 9
    new_cells%get_vertex_index(8,2) = 5
    new_cells%get_vertex_index(8,3) = 4

    new_cells%get_vertex_index(9,1) = 10
    new_cells%get_vertex_index(9,2) = 6
    new_cells%get_vertex_index(9,3) = 5

    new_cells%get_vertex_index(10,1) = 11
    new_cells%get_vertex_index(10,2) = 2
    new_cells%get_vertex_index(10,3) = 6

    new_cells%get_vertex_index(11,1) = 3
    new_cells%get_vertex_index(11,2) = 7
    new_cells%get_vertex_index(11,3) = 8

    new_cells%get_vertex_index(12,1) = 4
    new_cells%get_vertex_index(12,2) = 8
    new_cells%get_vertex_index(12,3) = 9

    new_cells%get_vertex_index(13,1) = 5
    new_cells%get_vertex_index(13,2) = 9
    new_cells%get_vertex_index(13,3) = 10

    new_cells%get_vertex_index(14,1) = 6
    new_cells%get_vertex_index(14,2) = 10
    new_cells%get_vertex_index(14,3) = 11

    new_cells%get_vertex_index(15,1) = 2
    new_cells%get_vertex_index(15,2) = 11
    new_cells%get_vertex_index(15,3) = 7

    new_cells%get_vertex_index(16,1) = 12
    new_cells%get_vertex_index(16,2) = 8
    new_cells%get_vertex_index(16,3) = 7

    new_cells%get_vertex_index(17,1) = 12
    new_cells%get_vertex_index(17,2) = 9
    new_cells%get_vertex_index(17,3) = 8

    new_cells%get_vertex_index(18,1) = 12
    new_cells%get_vertex_index(18,2) = 10
    new_cells%get_vertex_index(18,3) = 9

    new_cells%get_vertex_index(19,1) = 12
    new_cells%get_vertex_index(19,2) = 11
    new_cells%get_vertex_index(19,3) = 10

    new_cells%get_vertex_index(20,1) = 12
    new_cells%get_vertex_index(20,2) = 7
    new_cells%get_vertex_index(20,3) = 11

    ! fill edges
    CALL message('Creating icosahedron grid','fill edges...')
    vertex_edge(:,:) = 0
    max_edge = 0
    DO cell=1,no_of_output_cells
      vertex1 = new_cells%get_vertex_index(cell,max_cell_vertices)
      DO j=1,max_cell_vertices
        vertex2 = new_cells%get_vertex_index(cell,j)
        IF (vertex1 > no_of_output_verts .OR. vertex1 < 0) &
          CALL finish('get_icosahedron_grid','vertex1 > no_of_output_verts .OR. vertex1 < 0')
        IF (vertex2 > no_of_output_verts .OR. vertex2 < 0) &
          CALL finish('get_icosahedron_grid','vertex2 > no_of_output_verts .OR. vertex2 < 0')
        IF (vertex_edge(vertex1,vertex2) == 0) THEN
          ! add new edge
          max_edge = max_edge + 1
          edge = max_edge
          vertex_edge(vertex1,vertex2) = edge
          vertex_edge(vertex2,vertex1) = edge
          new_edges%get_vertex_index(edge,1) = vertex1
          new_edges%get_vertex_index(edge,2) = vertex2
          new_edges%get_cell_index(edge,1) = cell

        ELSE
          edge = vertex_edge(vertex1,vertex2)
          new_edges%get_cell_index(edge,2) = cell
        ENDIF

        new_cells%get_edge_index(cell,j) = edge

        ! add triangle to vertex2
        new_verts%no_of_neigbors(vertex2) = new_verts%no_of_neigbors(vertex2) + 1
        IF (new_verts%no_of_neigbors(vertex2) > max_vertex_connect) &
          & CALL finish('get_icosahedron_grid', &
          & 'new_verts%no_of_neigbors(vertex2) > max_vertex_connect')
        new_verts%get_cell_index(vertex2,   new_verts%no_of_neigbors(vertex2)) = cell

        ! shift vertices
        vertex1 = vertex2

      ENDDO ! j=1,max_cell_vertices
    ENDDO ! cell=1,no_of_output_cells

    ! check if we got the right numbers
    IF (max_edge /= no_of_output_edges) &
      CALL finish('get_icosahedron_grid','max_edge /= no_of_output_edges')

     !--------------------------------------------------------------
   ! fill vertices
    CALL message('Creating icosahedron grid','fill vertex coordinates...')
    DO vertex = 1,no_of_output_verts
      write(0,*) vertex, " base vertex:", base_vertex(vertex)
      new_verts%cartesian(vertex)%x = base_vertex(vertex)%x
!      new_verts%vertex(vertex) = cc2gc(new_verts%cartesian(vertex))
    ENDDO

    !--------------------------------------------------------------
    CALL message('Creating icosahedron grid','fill connectivity...')
    CALL complete_grid_connectivity(out_grid_id)
    !--------------------------------------------------------------
    ! connectivity is created
    ! temporarly set  parameters to default
    CALL set_nest_defaultindexes(out_grid_id)
    CALL set_no_of_subgrids(out_grid_id, 1)
    CALL set_start_subgrids(out_grid_id, 0)
    out_grid%is_filled = .true.
    !--------------------------------------------------------------
    ! create geometry
    ! this will be created by the caller
!     CALL compute_sphere_grid_geometry(out_grid_id)
!     icon_norm_edge = new_edges%primal_edge_length(1) / re
!     icon_norm_dual_edge = new_edges%dual_edge_length(1) / re

    RETURN
    
!     new_verts%no_of_neigbors(:) = 0
!     DO edge=1,no_of_output_edges
!       ! fill the vertex connectivity
!       vertex1 = new_edges%get_vertex_index(edge,1)
!       vertex2 = new_edges%get_vertex_index(edge,2)
!       IF (vertex1 > no_of_output_verts .OR. vertex1 < 0) &
!         CALL finish('get_icosahedron_grid','vertex1 > no_of_output_verts .OR. vertex1 < 0')
!       IF (vertex2 > no_of_output_verts .OR. vertex2 < 0) &
!         CALL finish('get_icosahedron_grid','vertex2 > no_of_output_verts .OR. vertex2 < 0')
!       new_verts%no_of_neigbors(vertex1) = new_verts%no_of_neigbors(vertex1) + 1
!       new_verts%no_of_neigbors(vertex2) = new_verts%no_of_neigbors(vertex2) + 1
!       IF (new_verts%no_of_neigbors(vertex1) > max_vertex_connect) &
!         & CALL finish('get_icosahedron_grid', &
!         & 'new_verts%no_of_neigbors(vertex1) > max_vertex_connect')
!       IF (new_verts%no_of_neigbors(vertex2) > max_vertex_connect) &
!         & CALL finish('get_icosahedron_grid', &
!         & 'new_verts%no_of_neigbors(vertex2) > max_vertex_connect')
! 
!       ! add edge to the vertices
!       new_verts%get_edge_index(vertex1,   new_verts%no_of_neigbors(vertex1)) = edge
!       new_verts%get_edge_index(vertex2,   new_verts%no_of_neigbors(vertex2)) = edge
! 
!       ! add the vertices to the neigbors list
!       new_verts%get_neighbor_index(vertex1,   new_verts%no_of_neigbors(vertex1)) = vertex2
!       new_verts%get_neighbor_index(vertex2,   new_verts%no_of_neigbors(vertex2)) = vertex1
! 
!       !fill cell neigbors
!       cell1 = new_edges%get_cell_index(edge,1)
!       cell2 = new_edges%get_cell_index(edge,2)
!       IF (cell1 > no_of_output_cells .OR. cell1 < 0) &
!         CALL finish('get_icosahedron_grid','cell1 > no_of_output_cells .OR. cell1 < 0')
!       IF (cell2 > no_of_output_cells .OR. cell2 < 0) &
!         CALL finish('get_icosahedron_grid','cell2 > no_of_output_cells .OR. cell2 < 0')
!       new_cells%no_of_vertices(cell1) =   new_cells%no_of_vertices(cell1) + 1
!       new_cells%get_neighbor_index(cell1, new_cells%no_of_vertices(cell1)) = cell2
!       new_cells%no_of_vertices(cell2) =   new_cells%no_of_vertices(cell2) + 1
!       new_cells%get_neighbor_index(cell2, new_cells%no_of_vertices(cell2)) = cell1
!     ENDDO ! edge=1,no_of_output_edges
! 
!     !--------------------------------------------------------------
!     ! connectivity is created
!     ! temporarly set  parameters to default
!     CALL set_nest_defaultindexes(out_grid_id)
!     CALL set_no_of_subgrids(out_grid_id, 1)
!     CALL set_start_subgrids(out_grid_id, 0)
!     out_grid%geometry_type   = sphere_geometry
!     out_grid%is_filled = .true.
!     !--------------------------------------------------------------
!     ! create geometry
!     CALL compute_sphere_grid_geometry(out_grid_id)
!     icon_norm_edge = new_edges%primal_edge_length(1) / re
!     icon_norm_dual_edge = new_edges%dual_edge_length(1) / re
! 
!     RETURN

  END FUNCTION get_icosahedron_grid
  !---------------------------------------------------------------------------------------------

END MODULE mo_icosahedron_grid


