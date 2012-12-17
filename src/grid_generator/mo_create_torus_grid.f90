!>
!! mo_create_torus provides the functionality for creating a torus grid
!!
!!
! The torus grid has the following form:
!
!      x=0       1       2       3       0
!      y=0____c_______d______e_______f____
!         \      /\      /\      /\      /\
!        a \    /  \    /  \    /  \    /  \ a
!           \  /    \  /    \  /    \  /    \
!      y=1   \/______\/______\/______\/______\
!            /\      /\      /\      /\      /
!         b /  \    /  \    /  \    /  \    / b
!          /    \  /    \  /    \  /    \  /
!         /___c__\/__d___\/__e___\/__f___\/
!      y=0
!      x=0       1       2       3       0
!
! where the edges with ths same letters are identical.
! The result is a double wrapped geometry.
! Wrapping along the y axis creates two surfaces (front and back)
! of a strip, while wrapping along the x axis creates a "cylinder"
! of the strip. Finally, the geometry is taken to lay on a plain.
!
! All triangles  are planar and equilateral
!
! The number of rows and columns may vary and are determined by
! the input parameters
!
! Creation Method
! The creation is based on a discrete modular coordinate system
! as shown above, representing the vertices of the grid.
! Each vertex, edge, and cell receives a unique index, based on
! this coordinate system. The connectivity is also created
! using these vertex coordinates, throught the indexing functions,
! edgeIndex_leftDiagonal(x,y), cellIndex_topRight(x,y), etc
!
! Main Routines
! SUBROUTINE create_torus_grid() :  creates the torus grid
! SUBROUTINE read_torus_grid_parameters() : reads the parameters of the grid
!       from the file 'create_torus.in'
!!
! Notes on orientation values
!
! . Normal primal vectors are stored on primal edges (edges%primal_normal).
!   They always point from edges%cell_index(edge_index,1) to edges%cell_index(edge_index,2)
!
! . Tangent primal vectors are stored on primal edges (edges%dual_normal).
!   They are on the right of the primal normal vectors, so that normal-tangent
!   form a left-handed system.
!
!   Orientation information on the primal edges
! . edges%system_orientation(edge_index) =
!        +1 if the tangent vector along the edge is from
!           edges%vertex_index(edge_index,1) to edges%vertex_index(edge_index,2)
!        otherwise is -1
!
!   Orientation information on the primal vertices
!   verts%edge_orient(this_vertex,neigbor_vertex)=
!        +1 if from this_vertex to neighbor_index(this_vertex,neigbor_vertex)
!           is the same orientation as the edge tangent vector
!        -1 if the opposite direction
!
!   Orientation information on the primal cells
!   cells%edge_orient(cell_index,cell_edge) =+1 if the primal normal on the edge defined by
!                                                    cells%edge_index(cell_index,cell_edge)
!                                                    points outwards form this cell
!                           -1 if the primal normal points inwards to this cell
!
!   Orientation with respect to the Gauss and Stokes formulae
!
! . For the Gauss formula on the primal grid:
!      Use cells%edge_orient: +1= positive (points outwards from the primal cell),
!                                  -1=negative (points inwards the primal cell)
!
! . For the Stokes formula on the primal grid. Positive orientation is counterclock-wise.
!      Use cells%edge_orient: -1=
!          positive orientation, +1=negative orienation (for the left-handed system)
!
! . For the Gauss formula on the dual grid:
!      Use verts%edge_orient: +1=positive (points outwards),
!                                  -1=negative (points inwards)
! . For the Stokes formula on the dual grid:
!      Use verts%edge_orient: +1=positive (is counterclock-wise),
!                                  -1=negative
!
!!
!! @par Revision History
!! First version  by Leonidas Linardakis, MPI-M, (2009-08-10)
!! Corrections    by Almut Grassman,      MPI-M, (2009-08-15)
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
MODULE mo_create_torus_grid
#include "grid_definitions.inc"

  USE mo_kind,            ONLY: wp
  USE mo_io_units,        ONLY: nnml, filename_max
  USE mo_namelist,        ONLY: position_nml, open_nml, positioned
  USE mo_exception,       ONLY: message, message_text, finish
  USE mo_grid_geometry_info,  ONLY: planar_torus_geometry
  USE mo_local_grid,      ONLY: t_grid, &
    & new_grid, delete_grid,get_grid, allocate_grid_object,        &
    & undefined, set_grid_creation, grid_set_exist_eq_allocated,   &
    & set_nest_defaultindexes
  !  & grid_cells, grid_edges, grid_vertices
  USE mo_io_local_grid,   ONLY: write_netcdf_grid
  USE mo_math_constants,  ONLY: pi
  USE mo_local_grid_geometry,  ONLY: order_cell_connectivity
  USE mo_math_utilities,  ONLY: plane_torus_distance, plane_torus_closest_coordinates

  IMPLICIT NONE

  PRIVATE

  !--------------------------------------------------------------
  !  Public subroutines
  PUBLIC create_torus_grid

  !--------------------------------------------------------------
  !  Grid Parameters
  INTEGER  :: y_no_of_rows=4
  INTEGER  :: x_no_of_columns=8
  REAL(wp) :: edge_length=1000.0_wp
  REAL(wp) :: x_center=0.0_wp
  REAL(wp) :: y_center=0.0_wp

  ! Grid structures
  TYPE(t_grid), POINTER :: torus_grid     ! the torus grid
  TYPE(t_grid), POINTER :: unfolded_grid  ! the unfolded grid for displaying the results
  INTEGER :: torus_grid_id, unfolded_grid_id

  INTEGER :: max_cell_vertices, max_vertex_connect
  INTEGER :: no_of_cells,no_of_edges,no_of_vertices
  INTEGER :: no_of_unfolded_cells,no_of_unfolded_edges,no_of_unfolded_vertices
  INTEGER :: unfoldno_hor_edges, unfoldno_diag_edges

  CHARACTER(LEN=filename_max) :: out_file_name='torus_triangle_grid.nc'
  CHARACTER(LEN=filename_max) :: unfolded_torus_file_name='unfoldedTorus_triangle_grid.nc'
  CHARACTER(LEN=filename_max) :: ascii_filename='unfoldedTorus_triangle_grid.ascii'


CONTAINS

  !---------------------------------------------------------
  !>
  !! SUBROUTINE read_torus_grid_parameters()
  !! Reads the parameters for creating the torus grid
  !-------------------------------------------------------------------------
  SUBROUTINE read_torus_grid_parameters(param_file_name)

    CHARACTER(LEN=*), INTENT(in) :: param_file_name

    INTEGER :: i_status
    CHARACTER(*), PARAMETER :: method_name = "read_torus_grid_parameters"

    NAMELIST /torus_grid_parameters/ y_no_of_rows, x_no_of_columns, &
      & edge_length, x_center, y_center, &
      & out_file_name, unfolded_torus_file_name, ascii_filename

    ! set default values
    y_no_of_rows = 4
    x_no_of_columns = 120
    edge_length = 1.0_wp
    x_center = 0.0_wp
    y_center = 0.0_wp
    out_file_name = "torus_grid.nc"
    unfolded_torus_file_name = "unfoldedTorus_grid.nc"
    ascii_filename = "unfoldedTorus_grid.ascii"

    ! read namelist
    CALL open_nml(param_file_name)
    CALL position_nml('torus_grid_parameters',STATUS=i_status)
    IF (i_status == positioned) THEN
      READ (nnml,torus_grid_parameters)
    ELSE
      WRITE(message_text,'(a)') " File", param_file_name, " not POSITIONED"
      CALL finish ('read_torus_grid_parameters', message_text)
    ENDIF
    CLOSE(nnml)

    WRITE(message_text,'(a)') "===================================="
    CALL message ('', TRIM(message_text))
    WRITE(message_text,'(a,i9)')   'y_noOfRows=', y_no_of_rows
    CALL message ('', TRIM(message_text))
    WRITE(message_text,'(a,i9)')   'x_noOfColumns=', x_no_of_columns
    CALL message ('', TRIM(message_text))
    WRITE(message_text,'(a,f8.2)') 'edge_length=', edge_length
    CALL message ('', TRIM(message_text))
    WRITE(message_text,'(a,f8.2)') 'x_center=', x_center
    CALL message ('', TRIM(message_text))
    WRITE(message_text,'(a,f8.2)') 'y_center=', y_center
    CALL message ('', TRIM(message_text))
    WRITE(message_text,'(a,a)')    'out_file_name=', TRIM(out_file_name)
    CALL message ('', TRIM(message_text))
    WRITE(message_text,'(a,a)')    'unfoldedTorusFileName=', TRIM(unfolded_torus_file_name)
    CALL message ('', TRIM(message_text))
    WRITE(message_text,'(a,a)')    'ascii_filename=', TRIM(ascii_filename)
    CALL message ('', TRIM(message_text))

    IF (y_no_of_rows < 2) &
      CALL finish(method_name, "y_no_of_rows must be at least 2")
    IF (x_no_of_columns < 2) &
      CALL finish(method_name, "x_no_of_columns must be at least 2")
    IF (MOD(y_no_of_rows, 2) /= 0) &
      CALL finish(method_name, "y_no_of_rows must be even")
    IF (MOD(x_no_of_columns, 2) /= 0) &
      CALL finish(method_name, "x_no_of_columns must be even")

  END SUBROUTINE read_torus_grid_parameters
  !-------------------------------------------------------------------------


  !---------------------------------------------------------
  !   SUBROUTINE create_torus_grid()
  !>
  !! The main subroutine for creating the torus grid
  !!
  SUBROUTINE create_torus_grid(param_file_name)
    CHARACTER(LEN=*), INTENT(in) :: param_file_name

    !--------------------------------------------------------------
    CALL read_torus_grid_parameters(param_file_name)
    WRITE(message_text,*)  'Create Torus Triangular Grid ...'
    CALL message ('', TRIM(message_text))

    max_cell_vertices  = 3
    max_vertex_connect = 6

    no_of_cells     = (y_no_of_rows * x_no_of_columns) * 2
    no_of_edges     = (y_no_of_rows * x_no_of_columns) * 3
    no_of_vertices  = (y_no_of_rows * x_no_of_columns)

    torus_grid_id = new_grid()
    torus_grid => get_grid(torus_grid_id)

    torus_grid%ncells = no_of_cells
    torus_grid%nedges = no_of_edges
    torus_grid%nverts = no_of_vertices
    torus_grid%cells%max_no_of_vertices = max_cell_vertices
    torus_grid%verts%max_connectivity   = max_vertex_connect
    
    torus_grid%geometry_type = planar_torus_geometry

    CALL allocate_grid_object(torus_grid_id)
    CALL grid_set_exist_eq_allocated(torus_grid_id)
    CALL set_grid_creation(torus_grid_id, undefined)
    !--------------------------------------------------------------

    CALL create_torus_topology()
    CALL create_planar_torus_geometry()
    CALL order_cell_connectivity(torus_grid_id)
    CALL set_nest_defaultindexes(torus_grid_id)

    WRITE(message_text,*)  torus_grid%ncells, ' triangles were created for the torus grid.'
    CALL message ('', TRIM(message_text))
    !--------------------------------------------------------------
    !  write netcdf file
!    CALL set_grid_creation(torus_grid_id, undefined)
    CALL write_netcdf_grid(torus_grid_id, out_file_name)

    !  print to standard output
    CALL print_torus_grid()
    CALL delete_grid(torus_grid_id)
    !--------------------------------------------------------------

    !--------------------------------------------------------------
!almut    ! create unfolded grid
!almut    CALL create_unfolded_torus()
!almut    !  write netcdf file
!almut    CALL write_netcdf_grid(unfolded_grid_id, unfolded_torus_file_name)
!almut    !  write unfolded_grid to ascii file
!     CALL write_grid_ascii(torus_grid_id, ascii_filename)
    
    !--------------------------------------------------------------
  END SUBROUTINE create_torus_grid


  !--------------------------------------------------------------
  !  SUBROUTINE create_torus_topology()
  !>
  !! Creates the torus topology
  !!
  !! See the module prologue for a general description of the method
  SUBROUTINE create_torus_topology()

    INTEGER :: x, y, index_no, index_pnt, vertex_1,vertex_2

    !--------------------------------------------------------------
    !   get vertices
    DO x=0, x_no_of_columns-1
      DO y=0, y_no_of_rows-1
        index_no=vertex_index(x,y)
        torus_grid%verts%idx(index_no) = index_no
        ! we will fill the noOfNeigbors when we get the edges
        torus_grid%verts%no_of_neigbors(index_no) = 0
      ENDDO
    ENDDO

    !--------------------------------------------------------------
    !   get edges
    !   the vertex index is defined so that the edges%system_orientation is always positive
    DO x=0, x_no_of_columns-1
      DO y=0, y_no_of_rows-1
        ! left diagonal edges
        index_no=edge_index_left_diagonal(x,y)
        torus_grid%edges%idx(index_no) = index_no
        vertex_1 = vertex_index_left_diagonal(x,y)
        vertex_2 = vertex_index(x,y)
        CALL update_torus_edge(index_no,vertex_1,vertex_2, 1)
        ! must be done after updateEdge
        torus_grid%edges%get_cell_index(index_no,1) = cell_index_top_right(x-1,y)
        torus_grid%edges%get_cell_index(index_no,2) = cell_index_top(x,y)

        ! right_diagonal edges
        index_no=edge_index_right_diagonal(x,y)
        torus_grid%edges%idx(index_no) = index_no
        vertex_1 = vertex_index_right_diagonal(x,y)
        vertex_2 = vertex_index(x,y)
        CALL update_torus_edge(index_no,vertex_1,vertex_2, 1)
        ! must be done after updateEdge
        torus_grid%edges%get_cell_index(index_no,1) = cell_index_top(x,y)
        torus_grid%edges%get_cell_index(index_no,2) = cell_index_top_right(x,y)

        ! horizontal edges
        index_no=edge_index_horizontal(x,y)
        torus_grid%edges%idx(index_no) = index_no
        vertex_1 = vertex_index(x,y)
        vertex_2 = vertex_index(x+1,y)
        CALL update_torus_edge(index_no,vertex_1,vertex_2, 1)
        torus_grid%edges%get_cell_index(index_no,1) = cell_index_bottom_right(x,y)
        torus_grid%edges%get_cell_index(index_no,2) = cell_index_top_right(x,y)

      ENDDO
    ENDDO

    !--------------------------------------------------------------
    !   zero vertices for cell updating
    DO x=0, x_no_of_columns-1
      DO y=0, y_no_of_rows-1
        index_no=vertex_index(x,y)
        torus_grid%verts%no_of_neigbors(index_no) = 0
      ENDDO
    ENDDO
    !--------------------------------------------------------------

    !--------------------------------------------------------------
    !   get cells
    DO x=0, x_no_of_columns-1
      DO y=0, y_no_of_rows-1

        !--------------------------------------------------------------
        index_no=cell_index_top(x,y)

        torus_grid%cells%idx(index_no) = index_no
        torus_grid%cells%no_of_vertices(index_no) = 3

        ! update cell vertices
        CALL update_torus_cell_vertex(index_no,x,y,1)
        CALL update_torus_cell_vertex(index_no,x-MOD(y,2),y+1,2)
        CALL update_torus_cell_vertex(index_no,x+1-MOD(y,2),y+1,3)

        ! update cell edges and neigbor cells
        index_pnt=edge_index_left_diagonal(x,y)
        torus_grid%cells%get_edge_index(index_no,1) = index_pnt
        index_pnt=cell_index_top_right(x-1,y)
        torus_grid%cells%get_neighbor_index(index_no,1) = index_pnt
        torus_grid%cells%get_edge_orient(index_no,1) = -1  ! normal is inwards

        index_pnt=edge_index_horizontal(x-MOD(y,2),y+1)
        torus_grid%cells%get_edge_index(index_no,2) = index_pnt
        index_pnt=cell_index_top_right(x-MOD(y,2),y+1)
        torus_grid%cells%get_neighbor_index(index_no,2) = index_pnt
        torus_grid%cells%get_edge_orient(index_no,2) = 1  ! normal is outwards

        index_pnt=edge_index_right_diagonal(x,y)
        torus_grid%cells%get_edge_index(index_no,3) = index_pnt
        index_pnt=cell_index_top_right(x,y)
        torus_grid%cells%get_neighbor_index(index_no,3) = index_pnt
        torus_grid%cells%get_edge_orient(index_no,3) = 1  ! normal is outwards
        !--------------------------------------------------------------
        index_no=cell_index_top_right(x,y)

        torus_grid%cells%idx(index_no) = index_no
        torus_grid%cells%no_of_vertices(index_no) = 3

        ! update cell vertices
        CALL update_torus_cell_vertex(index_no,x,y,1)
        CALL update_torus_cell_vertex(index_no,x+1-MOD(y,2),y+1,2)
        CALL update_torus_cell_vertex(index_no,x+1,y,3)

        ! update cell edges and neigbor cells
        index_pnt=edge_index_right_diagonal(x,y)
        torus_grid%cells%get_edge_index(index_no,1) = index_pnt
        index_pnt=cell_index_top(x,y)
        torus_grid%cells%get_neighbor_index(index_no,1) = index_pnt
        torus_grid%cells%get_edge_orient(index_no,1) = -1  ! normal is inwards

        index_pnt=edge_index_left_diagonal(x+1,y)
        torus_grid%cells%get_edge_index(index_no,2) = index_pnt
        index_pnt=cell_index_top(x+1,y)
        torus_grid%cells%get_neighbor_index(index_no,2) = index_pnt
        torus_grid%cells%get_edge_orient(index_no,2) = 1  ! normal is outwards

        index_pnt=edge_index_horizontal(x,y)
        torus_grid%cells%get_edge_index(index_no,3) = index_pnt
        index_pnt=cell_index_bottom_right(x,y)
        torus_grid%cells%get_neighbor_index(index_no,3) = index_pnt
        torus_grid%cells%get_edge_orient(index_no,3) = -1  ! normal is inwards
        !--------------------------------------------------------------

      ENDDO
    ENDDO

  END SUBROUTINE create_torus_topology
  !--------------------------------------------------------------

  !--------------------------------------------------------------
  !  SUBROUTINE create_planar_torus_geometry()
  !>
  !! Calculates the torus geometry (distances and areas)
  !!
  !! The geometry is planar (although it is not "real", it does not create conflicts)
  !! All triangles  are on the plain and equilateral
  !
  ! NOTE: the x,y coordinates are for display purposes
  ! they DO NOT wrap around and DO NOT represent the geometry correctly
  !
  ! The normals, lengths and the areas though are correct.
  !--------------------------------------------------------------
  SUBROUTINE create_planar_torus_geometry()

    INTEGER :: x, y, index_no
    ! the geometry (distances, areas) fields
    REAL(wp) :: x_lon_start,x_lon_row_starts,y_lat_start,x_lon_ref,y_lat_ref,x_lon_step,y_lat_step
    REAL(wp) :: x_start, x_row_starts,y_start,x_ref,y_ref,x_step,y_step
    REAL(wp) :: dual_edge_length,triangle_area,hexagon_area
    REAL(wp) :: sin60

    REAL(wp), PARAMETER :: max_lat = pi / 18.0_wp

    sin60 = SQRT(0.75_wp) 
    !--------------------------------------------------------------
    ! GEOMETRY PART
    dual_edge_length  = edge_length * (1.0_wp / SQRT(3.0_wp))
    triangle_area     = edge_length * edge_length * SQRT(0.1875_wp)
    hexagon_area      = dual_edge_length * edge_length * 1.5_wp
    ! find x_lon_step, y_lat_step, x_lon_start, y_lat_start
    ! Almut: scale to the length of 2pi to be similar to Radians (as in the Model)
    ! it is centered at lon,lat=0,0
    x_lon_step  = (2.0_wp * pi)/  REAL(x_no_of_columns,wp)
    x_lon_start = - (REAL(x_no_of_columns,wp) * 0.5_wp) * x_lon_step
    y_lat_step  = (max_lat * 2.0_wp) / REAL(y_no_of_rows,wp)
    y_lat_start = - (REAL(y_no_of_rows,wp) * 0.5_wp) * y_lat_step
    
    x_step  = edge_length
    y_step  = edge_length * sin60
    x_start = x_center - (REAL(x_no_of_columns,wp) * 0.5_wp) * x_step
    y_start = y_center - (REAL(y_no_of_rows,wp) * 0.5_wp) * y_step
    
    !--------------------------------------------------------------
    ! write planar torus geometry properties
    torus_grid%planar_torus_info%center%x(1)  = x_center
    torus_grid%planar_torus_info%center%x(2)  = y_center
    torus_grid%planar_torus_info%center%x(3)  = 0.0_wp
    torus_grid%planar_torus_info%cell_edge_length  = edge_length
    torus_grid%planar_torus_info%length  = x_step * REAL(x_no_of_columns,wp)
    torus_grid%planar_torus_info%height  = y_step * REAL(y_no_of_rows,wp)
    
   !--------------------------------------------------------------
    !  get  coordinates
    DO y=0, y_no_of_rows-1
      x_lon_row_starts = x_lon_start - x_lon_step * 0.5_wp * REAL(MOD(y,2),wp)
      y_lat_ref        = y_lat_start + REAL(y,wp) * y_lat_step
      x_row_starts     = x_start - x_step * 0.5_wp * REAL(MOD(y,2),wp)
      y_ref            = y_start + REAL(y,wp) * y_step
    
      DO x=0, x_no_of_columns-1
        x_lon_ref        = x_lon_row_starts + REAL(x,wp) * x_lon_step        
        x_ref            = x_row_starts + REAL(x,wp) * x_step

        ! for vertices
        index_no=vertex_index(x,y)
        torus_grid%verts%vertex(index_no)%lon     = x_lon_ref
        torus_grid%verts%vertex(index_no)%lat     = y_lat_ref
        torus_grid%verts%cartesian(index_no)%x(1) = x_ref
        torus_grid%verts%cartesian(index_no)%x(2) = y_ref
        torus_grid%verts%cartesian(index_no)%x(3) = 0.0_wp
        torus_grid%verts%dual_area(index_no)      = hexagon_area

        ! for cells (top)
        index_no=cell_index_top(x,y)
        torus_grid%cells%center(index_no)%lon            = x_lon_ref
        torus_grid%cells%center(index_no)%lat            = y_lat_ref + y_lat_step * 2.0_wp / 3.0_wp
        torus_grid%cells%cartesian_center(index_no)%x(1) = x_ref
        torus_grid%cells%cartesian_center(index_no)%x(2) = y_ref     + y_step * 2.0_wp / 3.0_wp
        torus_grid%cells%cartesian_center(index_no)%x(3) = 0.0_wp

        ! for cells (right)
        index_no=cell_index_top_right(x,y)
        torus_grid%cells%center(index_no)%lon            = x_lon_ref + 0.5_wp * x_lon_step
        torus_grid%cells%center(index_no)%lat            = y_lat_ref + y_lat_step / 3.0_wp
        torus_grid%cells%cartesian_center(index_no)%x(1) = x_ref     + 0.5_wp * x_step
        torus_grid%cells%cartesian_center(index_no)%x(2) = y_ref     + y_step / 3.0_wp
        torus_grid%cells%cartesian_center(index_no)%x(3) = 0.0_wp

      ENDDO
    ENDDO

    !--------------------------------------------------------------
    !  get edge geometry
    ! the dual_normal is on the right of the primal_normal (left-handed system)
    DO y=0, y_no_of_rows-1
      x_lon_row_starts = x_lon_start - x_lon_step * 0.5_wp * REAL(MOD(y,2),wp)
      y_lat_ref        =  y_lat_start + REAL(y,wp) * y_lat_step
      x_row_starts     = x_start - x_step * 0.5_wp * REAL(MOD(y,2),wp)
      y_ref            =  y_start + REAL(y,wp) * y_step
      
      DO x=0, x_no_of_columns-1

        x_lon_ref = x_lon_row_starts + REAL(x,wp) * x_lon_step
        x_ref     = x_row_starts + REAL(x,wp) * x_step

        ! left_diagonal edges
        index_no = edge_index_left_diagonal(x,y)
        torus_grid%edges%center(index_no)%lon            = x_lon_ref - x_lon_step * 0.25_wp
        torus_grid%edges%center(index_no)%lat            = y_lat_ref + (y_lat_step * 0.5_wp)
        torus_grid%edges%cartesian_center(index_no)%x(1) = x_ref     - x_step * 0.25_wp
        torus_grid%edges%cartesian_center(index_no)%x(2) = y_ref     + (y_step * 0.5_wp)
        torus_grid%edges%cartesian_center(index_no)%x(3) = 0.0_wp        
        torus_grid%edges%primal_normal(index_no)%v1 = sin60
        torus_grid%edges%primal_normal(index_no)%v2 = 0.5_wp
        torus_grid%edges%dual_normal(index_no)%v1   = 0.5_wp
        torus_grid%edges%dual_normal(index_no)%v2   = -sin60

        ! right_diagonal edges
        index_no=edge_index_right_diagonal(x,y)
        torus_grid%edges%center(index_no)%lon            = x_lon_ref + x_lon_step * 0.25_wp
        torus_grid%edges%center(index_no)%lat            = y_lat_ref + (y_lat_step * 0.5_wp)
        torus_grid%edges%cartesian_center(index_no)%x(1) = x_ref     + x_step * 0.25_wp
        torus_grid%edges%cartesian_center(index_no)%x(2) = y_ref     + (y_step * 0.5_wp)
        torus_grid%edges%cartesian_center(index_no)%x(3) = 0.0_wp        
        torus_grid%edges%primal_normal(index_no)%v1 = sin60
        torus_grid%edges%primal_normal(index_no)%v2 = -0.5_wp
        torus_grid%edges%dual_normal(index_no)%v1   = -0.5_wp
        torus_grid%edges%dual_normal(index_no)%v2   = -sin60

        ! horizontal edges
        index_no=edge_index_horizontal(x,y)
        torus_grid%edges%center(index_no)%lon            = x_lon_ref + x_lon_step * 0.5_wp
        torus_grid%edges%center(index_no)%lat            = y_lat_ref
        torus_grid%edges%cartesian_center(index_no)%x(1) = x_ref     + x_step * 0.5_wp
        torus_grid%edges%cartesian_center(index_no)%x(2) = y_ref
        torus_grid%edges%cartesian_center(index_no)%x(3) = 0.0_wp        
        torus_grid%edges%primal_normal(index_no)%v1 = 0.0_wp
        torus_grid%edges%primal_normal(index_no)%v2 = 1.0_wp
        torus_grid%edges%dual_normal(index_no)%v1   = 1.0_wp
        torus_grid%edges%dual_normal(index_no)%v2   = 0.0_wp
      ENDDO
    ENDDO

    DO index_no=1, no_of_edges
      torus_grid%edges%primal_edge_length(index_no) = edge_length
      torus_grid%edges%get_edge_vert_length(index_no,1) = edge_length * 0.5_wp
      torus_grid%edges%get_edge_vert_length(index_no,2) = edge_length * 0.5_wp

      torus_grid%edges%dual_edge_length(index_no)   = dual_edge_length
      torus_grid%edges%get_edge_cell_length(index_no,1) = dual_edge_length * 0.5_wp
      torus_grid%edges%get_edge_cell_length(index_no,2) = dual_edge_length * 0.5_wp
    ENDDO

    !--------------------------------------------------------------
    !  get cell geometry
    DO index_no=1, no_of_cells
      torus_grid%cells%area(index_no) = triangle_area
    ENDDO
    !--------------------------------------------------------------
  END SUBROUTINE create_planar_torus_geometry
  !--------------------------------------------------------------


  !--------------------------------------------------------------
  !  SUBROUTINE print_torus_grid()
  !>
  !! Prints the  torus grid to standard output
  !!
  SUBROUTINE print_torus_grid()

    INTEGER :: i, x, y, index_no, counter
    INTEGER :: vertex_no, next_vertex_no

    !--------------------------------------------------------------
    WRITE(message_text,*) '-------------------------------'
    WRITE(*,*) TRIM(message_text)
    WRITE(message_text,*) '    ==  TORUS GRID  ==='
    WRITE(*,*) TRIM(message_text)
    write(*,*) " Center:", torus_grid%planar_torus_info%center%x 
    WRITE(message_text,*) '-------------------------------'
    WRITE(*,*) TRIM(message_text)
    ! print vertices
    DO x=0, x_no_of_columns-1
      DO y=0, y_no_of_rows-1
        index_no=vertex_index(x,y)
        counter=torus_grid%verts%no_of_neigbors(index_no)
        WRITE(message_text,*) 'vertex:',index_no, '(',x,y,')'
        WRITE(*,*) TRIM(message_text)
        WRITE(message_text,*) 'lon  ', torus_grid%verts%vertex(index_no)%lon
        WRITE(*,*) TRIM(message_text)
        WRITE(message_text,*) 'lat  ', torus_grid%verts%vertex(index_no)%lat
        WRITE(*,*) TRIM(message_text)
        WRITE(*,*) " cartesian coord:", torus_grid%verts%cartesian(index_no)%x
        WRITE(message_text,*) 'neigs', &
          (torus_grid%verts%get_neighbor_index(index_no,i),i=1,counter)
        WRITE(*,*) TRIM(message_text)
        WRITE(message_text,*) 'cells', &
          (torus_grid%verts%get_cell_index(index_no,i),i=1,counter)
        WRITE(*,*) TRIM(message_text)
        WRITE(message_text,*) 'edges', &
          (torus_grid%verts%get_edge_index(index_no,i),i=1,counter)
        WRITE(*,*) TRIM(message_text)
        WRITE(message_text,*) 'edori', &
          (torus_grid%verts%get_edge_orient(index_no,i),i=1,counter)
        WRITE(*,*) TRIM(message_text)
      ENDDO
    ENDDO
    ! print edges
    DO x=0, x_no_of_columns-1
      DO y=0, y_no_of_rows-1
        index_no=edge_index_left_diagonal(x,y)
        WRITE(message_text,*) 'vertical edge:',index_no, '(',x,y,')'
        write(*,*) TRIM(message_text)
        WRITE(*,*) TRIM(message_text)
        WRITE(message_text,*) 'sysor', torus_grid%edges%system_orientation(index_no)
        WRITE(*,*) TRIM(message_text)
        WRITE(message_text,*) 'lon  ', torus_grid%edges%center(index_no)%lon
        WRITE(*,*) TRIM(message_text)
        WRITE(message_text,*) 'lat  ', torus_grid%edges%center(index_no)%lat
        WRITE(*,*) TRIM(message_text)
        WRITE(*,*) "edge center cartesian coord:", torus_grid%edges%cartesian_center(index_no)%x
        WRITE(message_text,*) 'verts', &
          (torus_grid%edges%get_vertex_index(index_no,i),i=1,2)
        WRITE(*,*) TRIM(message_text)
        WRITE(message_text,*) 'cells', &
          (torus_grid%edges%get_cell_index(index_no,i),i=1,2)
        WRITE(*,*) TRIM(message_text)

        index_no=edge_index_right_diagonal(x,y)
        WRITE(message_text,*) 'diagonal edge:',index_no, '(',x,y,')'
        WRITE(*,*) TRIM(message_text)
        WRITE(message_text,*) 'sysor', &
          torus_grid%edges%system_orientation(index_no)
        WRITE(*,*) TRIM(message_text)
        WRITE(message_text,*) 'lon  ', &
          torus_grid%edges%center(index_no)%lon
        WRITE(*,*) TRIM(message_text)
        WRITE(message_text,*) 'lat  ', torus_grid%edges%center(index_no)%lat
        WRITE(*,*) TRIM(message_text)
        WRITE(*,*) "edge center cartesian coord:", torus_grid%edges%cartesian_center(index_no)%x
        WRITE(message_text,*) 'verts', &
          (torus_grid%edges%get_vertex_index(index_no,i),i=1,2)
        WRITE(*,*) TRIM(message_text)
        WRITE(message_text,*) 'cells', &
          (torus_grid%edges%get_cell_index(index_no,i),i=1,2)
        WRITE(*,*) TRIM(message_text)

        index_no=edge_index_horizontal(x,y)
        WRITE(message_text,*) 'horizontal edge:',index_no, '(',x,y,')'
        WRITE(*,*) TRIM(message_text)
        WRITE(message_text,*) 'sysor', &
          torus_grid%edges%system_orientation(index_no)
        WRITE(*,*) TRIM(message_text)
        WRITE(message_text,*) 'lon  ', torus_grid%edges%center(index_no)%lon
        WRITE(*,*) TRIM(message_text)
        WRITE(message_text,*) 'lat  ', torus_grid%edges%center(index_no)%lat
        WRITE(*,*) TRIM(message_text)
        WRITE(*,*) "edge center cartesian coord:", torus_grid%edges%cartesian_center(index_no)%x
        WRITE(message_text,*) 'verts', &
          (torus_grid%edges%get_vertex_index(index_no,i),i=1,2)
        WRITE(*,*) TRIM(message_text)
        WRITE(message_text,*) 'cells', &
          (torus_grid%edges%get_cell_index(index_no,i),i=1,2)
        WRITE(*,*) TRIM(message_text)

      ENDDO
    ENDDO 
    ! print cells
    DO x=0, x_no_of_columns-1
      DO y=0, y_no_of_rows-1
        index_no=cell_index_top(x,y)
        WRITE(message_text,*) 'top triangle:',index_no, '(',x,y,')'
        WRITE(*,*) TRIM(message_text)
        WRITE(message_text,*) 'lon  ', torus_grid%cells%center(index_no)%lon
        WRITE(*,*) TRIM(message_text)
        WRITE(message_text,*) 'lat  ', torus_grid%cells%center(index_no)%lat
        WRITE(*,*) TRIM(message_text)
        WRITE(*,*) "cell center cartesian coord:", torus_grid%cells%cartesian_center(index_no)%x
        WRITE(message_text,*) 'neigs', &
          (torus_grid%cells%get_neighbor_index(index_no,i),i=1,3)
        WRITE(*,*) TRIM(message_text)
        WRITE(message_text,*) 'verts', &
          (torus_grid%cells%get_vertex_index(index_no,i),i=1,3)
        WRITE(*,*) TRIM(message_text)
        WRITE(message_text,*) 'edges', &
          (torus_grid%cells%get_edge_index(index_no,i),i=1,3)
        WRITE(*,*) TRIM(message_text)
        WRITE(message_text,*) 'edori', &
          (torus_grid%cells%get_edge_orient(index_no,i),i=1,3)
        WRITE(*,*) TRIM(message_text)

        index_no=cell_index_top_right(x,y)
        WRITE(message_text,*) 'right triangle:',index_no, '(',x,y,')'
        WRITE(*,*) TRIM(message_text)
        WRITE(message_text,*) 'lon  ', torus_grid%cells%center(index_no)%lon
        WRITE(*,*) TRIM(message_text)
        WRITE(message_text,*) 'lat  ', torus_grid%cells%center(index_no)%lat
        WRITE(*,*) TRIM(message_text)
        WRITE(message_text,*) 'neigs', &
          (torus_grid%cells%get_neighbor_index(index_no,i),i=1,3)
        WRITE(*,*) TRIM(message_text)
        WRITE(message_text,*) 'verts', &
          (torus_grid%cells%get_vertex_index(index_no,i),i=1,3)
        WRITE(*,*) TRIM(message_text)
        WRITE(message_text,*) 'edges', &
          (torus_grid%cells%get_edge_index(index_no,i),i=1,3)
        WRITE(*,*) TRIM(message_text)
        WRITE(message_text,*) 'edori', &
          (torus_grid%cells%get_edge_orient(index_no,i),i=1,3)
        WRITE(*,*) TRIM(message_text)
      ENDDO
    ENDDO

    ! check the plane_torus_distance, plane_torus_closest_coordinates routines
    DO vertex_no = 1, torus_grid%verts%no_of_existvertices
      DO next_vertex_no = 1, torus_grid%verts%no_of_existvertices
        write(*,*) "---------- check plane_torus_distance, plane_torus_closest_coordinates ------------"
        write(*,*) " Coordinates 0:", torus_grid%verts%cartesian(vertex_no)%x
        write(*,*) " Coordinates 1:", torus_grid%verts%cartesian(next_vertex_no)%x
        write(*,*) " New coordinates 1:", plane_torus_closest_coordinates(&
          & torus_grid%verts%cartesian(vertex_no),      &
          & torus_grid%verts%cartesian(next_vertex_no), &
          & torus_grid%planar_torus_info%length, torus_grid%planar_torus_info%height)
        write(*,*) " distance:", plane_torus_distance(  &
          & torus_grid%verts%cartesian(vertex_no),      &
          & torus_grid%verts%cartesian(next_vertex_no), &
          & torus_grid%planar_torus_info%length, torus_grid%planar_torus_info%height)
      ENDDO
    ENDDO
    
    
  END SUBROUTINE print_torus_grid





  !--------------------------------------------------------------
  !  SUBROUTINE write_grid_ascii()
  !>
  !! Writes a grid in ascii format, for displaying
  !!
  ! SUBROUTINE write_grid_ascii(write_grid, outfile)
  ! !
  !   USE mo_kind
  !   USE mo_io_units,  ONLY: filename_max
  !   USE mo_exception, ONLY: message, message_text
  !   USE mo_grid
  !
  ! !--------------------------------------------------------------
  !   IMPLICIT NONE
  !
  !   TYPE(grid), INTENT(in) :: write_grid
  !   CHARACTER(len=filename_max), INTENT(in) :: outfile
  !
  !   INTEGER :: no_of_cells,no_of_edges,no_of_vertices,no_of_cell_verts
  !   INTEGER :: i, k
  ! !--------------------------------------------------------------
  !
  !
  !     OPEN(500, FILE=trim(outfile))
  !     WRITE(500,*) 'Grid Ascii File'
  !
  !     no_of_cells = write_grid%ncells
  !     no_of_edges = write_grid%nedges
  !     no_of_vertices = write_grid%nverts
  !
  !     WRITE(500,'(i8," ",i8)') no_of_cells, no_of_vertices
  !     DO i=1,no_of_vertices
  !        WRITE(500,'(i8,3(" ",f13.6))') i, write_grid%verts%vertex(i)%lon, &
  !                 write_grid%verts%vertex(i)%lat, 0.0_wp
  !     ENDDO
  !     DO i=1,no_of_cells
  !        no_of_cell_verts = write_grid%cells%no_of_vertices(i)
  !        WRITE(500,'(i8,12(" ",i8))') i, no_of_cell_verts, &
  !           (write_grid%cells% get vertex_index(i,k), k=1,no_of_cell_verts)
  !     ENDDO
  !     WRITE(500,*)'0 \n0 \n0 \n0 '
  !     CLOSE(500)
  !
  !
  ! END SUBROUTINE write_grid_ascii
  !--------------------------------------------------------------

  !--------------------------------------------------------------
  !  SUBROUTINE updateEdge(torus_grid,indexNo,vertex_1,vertex_2, orientation)
  !>
  !! Adds the vertices of the edge(indexNo)
  !!
  !! vertex_1,vertex_2 are the two edge vertices
  !! For the orientation, see notes in the module prologue
  SUBROUTINE update_torus_edge(edge_index_no,vertex_1,vertex_2, orientation)
    INTEGER, INTENT(in) :: edge_index_no,vertex_1,vertex_2, orientation

    INTEGER :: counter
    !--------------------------------------------------------------

    torus_grid%edges%get_vertex_index(edge_index_no,1) = vertex_1
    torus_grid%edges%get_vertex_index(edge_index_no,2) = vertex_2
    torus_grid%edges%system_orientation(edge_index_no) = orientation
    torus_grid%edges%get_cell_index(edge_index_no,1) = -1
    torus_grid%edges%get_cell_index(edge_index_no,2) = -1

    ! update the two vertices
    counter=torus_grid%verts%no_of_neigbors(vertex_1)+1
    torus_grid%verts%get_neighbor_index(vertex_1,counter)= vertex_2
    torus_grid%verts%get_edge_index(vertex_1,counter)= edge_index_no
    torus_grid%verts%get_edge_orient(vertex_1,counter) = orientation
    torus_grid%verts%no_of_neigbors(vertex_1) = counter

    counter=torus_grid%verts%no_of_neigbors(vertex_2)+1
    torus_grid%verts%get_neighbor_index(vertex_2,counter)= vertex_1
    torus_grid%verts%get_edge_index(vertex_2,counter)= edge_index_no
    torus_grid%verts%get_edge_orient(vertex_2,counter) = -orientation
    torus_grid%verts%no_of_neigbors(vertex_2) = counter

  END SUBROUTINE update_torus_edge
  !--------------------------------------------------------------

  !--------------------------------------------------------------
  !  SUBROUTINE updateCellEdge(torus_grid,indexNo,indexPnt,k)
  !>
  !! Adds a new edge (edgeIndex) to cell(cellIndexNo) at position k
  SUBROUTINE update_torus_cell_edge(cell_index_no,edge_index,k)
    INTEGER, INTENT(in) :: cell_index_no,edge_index,k
    !--------------------------------------------------------------
    torus_grid%cells%get_edge_index(cell_index_no,k) = edge_index
  END SUBROUTINE update_torus_cell_edge

  !--------------------------------------------------------------
  !  SUBROUTINE updateCellVertex(torus_grid, indexNo,x,y,k)
  !>
  !! Adds a new vertex at (x,y) to the cell cellIndexNo, at position k
  SUBROUTINE update_torus_cell_vertex(cell_index_no,x,y,k)
    INTEGER, INTENT(in) :: cell_index_no, x, y,k

    INTEGER :: vertex_index_pnt,counter
    !--------------------------------------------------------------
    vertex_index_pnt=vertex_index(x,y)
    torus_grid%cells%get_vertex_index(cell_index_no,k) = vertex_index_pnt
    counter = torus_grid%verts%no_of_neigbors(vertex_index_pnt)+1
    torus_grid%verts%get_cell_index(vertex_index_pnt,counter)=cell_index_no
    torus_grid%verts%no_of_neigbors(vertex_index_pnt) = counter
  END SUBROUTINE update_torus_cell_vertex
  !--------------------------------------------------------------

  !--------------------------------------------------------------
  !FUNCTION cellIndex_topRight(x,y) RESULT(i)
  !>
  !! returns the cell index at the top-right of (x,y) vertex
  ELEMENTAL FUNCTION cell_index_top_right(x,y) result(i)
    INTEGER, INTENT(in) :: x, y
    INTEGER :: i,x_index,y_index
    !--------------------------------------------------------------
    x_index=MOD(x,x_no_of_columns)
    IF (x_index < 0) x_index = x_index + x_no_of_columns
    y_index=MOD(y,y_no_of_rows)
    IF (y_index < 0) y_index = y_index + y_no_of_rows

    i=(x_index * y_no_of_rows + y_index)*2 + 2

  END FUNCTION cell_index_top_right
  !--------------------------------------------------------------

  !--------------------------------------------------------------
  !FUNCTION cellIndex_bottomRight(x,y) RESULT(i)
  !>
  !! returns the cell index at the bottom-right of (x,y) vertex
  ELEMENTAL FUNCTION cell_index_bottom_right(x,y) result(i)
    INTEGER, INTENT(in) :: x, y
    INTEGER :: i,x_index,y_index
    !--------------------------------------------------------------
    y_index=y-1
    x_index=x+1-MOD(y,2)

    i=cell_index_top(x_index,y_index)

  END FUNCTION cell_index_bottom_right
  !--------------------------------------------------------------

  !--------------------------------------------------------------
  !FUNCTION cellIndex_top(x,y) RESULT(i)
  !>
  !! returns the cell index at the top of (x,y) vertex
  ELEMENTAL FUNCTION cell_index_top(x,y) result(i)
    INTEGER, INTENT(in) :: x, y
    INTEGER :: i,x_index,y_index
    !--------------------------------------------------------------
    x_index=MOD(x,x_no_of_columns)
    IF (x_index < 0) x_index = x_index + x_no_of_columns
    y_index=MOD(y,y_no_of_rows)
    IF (y_index < 0) y_index = y_index + y_no_of_rows

    i=(x_index * y_no_of_rows + y_index)*2 + 1

  END FUNCTION cell_index_top
  !--------------------------------------------------------------


  !--------------------------------------------------------------
  !FUNCTION vertexIndex_leftDiagonal(x,y) RESULT(i)
  !>
  !! returns the vertex index to which the left diagonal edge from (x,y) ends
  ELEMENTAL FUNCTION vertex_index_left_diagonal(x,y) result(i)
    INTEGER, INTENT(in) :: x, y
    INTEGER :: i,x_index,y_index
    !--------------------------------------------------------------
    y_index=y+1
    x_index=x-MOD(y,2)
    i=vertex_index(x_index,y_index)

  END FUNCTION vertex_index_left_diagonal
  !--------------------------------------------------------------

  !--------------------------------------------------------------
  !FUNCTION vertexIndex_rightDiagonal(x,y) RESULT(i)
  !>
  !! returns the vertex index to which the right diagonal edge from (x,y) ends
  ELEMENTAL FUNCTION vertex_index_right_diagonal(x,y) result(i)
    INTEGER, INTENT(in) :: x, y
    INTEGER :: i,x_index,y_index
    !--------------------------------------------------------------
    y_index=y+1
    x_index=x+1-MOD(y,2)
    i=vertex_index(x_index,y_index)

  END FUNCTION vertex_index_right_diagonal
  !--------------------------------------------------------------


  !--------------------------------------------------------------
  !FUNCTION vertexIndex(x,y) RESULT(i)
  !>
  !! returns the vertex index at (x,y)
  ELEMENTAL FUNCTION vertex_index(x,y) result(i)
    INTEGER, INTENT(in) :: x, y
    INTEGER :: i,x_index,y_index
    !--------------------------------------------------------------
    x_index=MOD(x,x_no_of_columns)
    IF (x_index < 0) x_index = x_index + x_no_of_columns
    y_index=MOD(y,y_no_of_rows)
    IF (y_index < 0) y_index = y_index + y_no_of_rows

    i=x_index * y_no_of_rows + y_index + 1

  END FUNCTION vertex_index
  !--------------------------------------------------------------


  !--------------------------------------------------------------
  !FUNCTION edgeIndex_leftDiagonal(x,y) RESULT(i)
  !>
  !! returns the index of the left-diagonal edge on the top of (x,y) point
  ELEMENTAL FUNCTION edge_index_left_diagonal(x,y) result(i)
    INTEGER, INTENT(in) :: x, y
    INTEGER :: i,x_index,y_index
    !--------------------------------------------------------------
    x_index=MOD(x,x_no_of_columns)
    IF (x_index < 0) x_index = x_index + x_no_of_columns
    y_index=MOD(y,y_no_of_rows)
    IF (y_index < 0) y_index = y_index + y_no_of_rows

    i=(x_index * y_no_of_rows + y_index) * 3 + 1

  END FUNCTION edge_index_left_diagonal
  !--------------------------------------------------------------

  !--------------------------------------------------------------
  !FUNCTION edgeIndex_horizontal(x,y) RESULT(i)
  !>
  !! returns the index of the horizonal edge on the right of the (x,y) point
  ELEMENTAL FUNCTION edge_index_horizontal(x,y) result(i)
    INTEGER, INTENT(in) :: x, y
    INTEGER :: i,x_index,y_index
    !--------------------------------------------------------------
    x_index=MOD(x,x_no_of_columns)
    IF (x_index < 0) x_index = x_index + x_no_of_columns
    y_index=MOD(y,y_no_of_rows)
    IF (y_index < 0) y_index = y_index + y_no_of_rows

    i=(x_index * y_no_of_rows + y_index) * 3 + 3

  END FUNCTION edge_index_horizontal
  !--------------------------------------------------------------

  !--------------------------------------------------------------
  !FUNCTION eedgeIndex_rightDiagonal(x,y) RESULT(i)
  !>
  !! returns the index of the right_diagonal edge on the top-right of the (x,y) point
  ELEMENTAL FUNCTION edge_index_right_diagonal(x,y) result(i)
    INTEGER, INTENT(in) :: x, y
    INTEGER :: i,x_index,y_index
    !--------------------------------------------------------------
    x_index=MOD(x,x_no_of_columns)
    IF (x_index < 0) x_index = x_index + x_no_of_columns
    y_index=MOD(y,y_no_of_rows)
    IF (y_index < 0) y_index = y_index + y_no_of_rows

    i=(x_index * y_no_of_rows + y_index) * 3 + 2

  END FUNCTION edge_index_right_diagonal
  !--------------------------------------------------------------

  !--------------------------------------------------------------
  ! the calculation for the unfolded indexes
  ELEMENTAL INTEGER FUNCTION unfolded_vertex_index(x,y)
    INTEGER, INTENT(in) :: x,y
    unfolded_vertex_index = ((x) * (y_no_of_rows+1) + (y) + 1)
  END FUNCTION


  ELEMENTAL INTEGER FUNCTION unfoldleftright_diag_edgeindex(x,y)
    INTEGER, INTENT(in) :: x,y
    unfoldleftright_diag_edgeindex=((x)*y_no_of_rows+(y)+unfoldno_hor_edges+unfoldno_diag_edges+1)
  END FUNCTION

  ELEMENTAL INTEGER FUNCTION unfoldright_diag_edgeindex(x,y)
    INTEGER, INTENT(in) :: x,y
    unfoldright_diag_edgeindex = ((x) * y_no_of_rows + (y)  + unfoldno_hor_edges + 1)
  END FUNCTION

  ELEMENTAL INTEGER FUNCTION unfolded_horizontal_edge_index(x,y)
    INTEGER, INTENT(in) :: x,y
    unfolded_horizontal_edge_index = ((x) * (y_no_of_rows + 1) + (y) + 1)
  END FUNCTION

  INTEGER FUNCTION unfolded_cell_index_top(x,y)
    INTEGER, INTENT(in) :: x,y
    unfolded_cell_index_top = cell_index_top(x,y)
  END FUNCTION

  INTEGER FUNCTION unfolded_cell_index_right(x,y)
    INTEGER, INTENT(in) :: x,y
    unfolded_cell_index_right = cell_index_top_right(x,y)
  END FUNCTION



  !--------------------------------------------------------------
  !  SUBROUTINE create_unfolded_torus()
  !>
  !! Creates the unfolded version of the torus grid for displaying purposes
  SUBROUTINE create_unfolded_torus()
    INTEGER :: x, y, index_no,  unfolded_index_no
    ! the geometry (distances, areas) fields
    REAL(wp) :: x_lon_start,y_lat_start,x_lon_ref,y_lat_ref,x_lon_step,y_lat_step


    !--------------------------------------------------------------
    ! create unfolded grid for visualization only
    ! the idx points to the torus entities, where we will get the values for visualization
    !--------------------------------------------------------------
    max_cell_vertices  = 3
    max_vertex_connect = 6
    unfoldno_hor_edges = ((y_no_of_rows+1) * x_no_of_columns)
    unfoldno_diag_edges = (y_no_of_rows * x_no_of_columns)

    no_of_unfolded_cells     = (y_no_of_rows*x_no_of_columns) * 2 ! the same as the torus
    no_of_unfolded_edges     = (y_no_of_rows*x_no_of_columns) * 3 &
      & + y_no_of_rows + x_no_of_columns     ! torus +  the cut edges
    no_of_unfolded_vertices  = (y_no_of_rows*x_no_of_columns)     &
      & + y_no_of_rows + x_no_of_columns + 1 ! torus +  the cut vertices
    !--------------------------------------------------------------
    x_lon_step  = (2.0_wp * pi)/ (REAL(x_no_of_columns,wp)*1.0_wp)
    y_lat_step  = x_lon_step * SQRT(0.75_wp)
    x_lon_start = x_center - (REAL(x_no_of_columns,wp) * 0.5_wp) * x_lon_step
    y_lat_start = y_center - (REAL(y_no_of_rows,wp) * 0.5_wp) * y_lat_step
    !--------------------------------------------------------------

    !--------------------------------------------------------------
    unfolded_grid_id = new_grid()
    unfolded_grid => get_grid(unfolded_grid_id)

    unfolded_grid%ncells = no_of_cells
    unfolded_grid%nedges = no_of_edges
    unfolded_grid%nverts = no_of_vertices
    unfolded_grid%cells%max_no_of_vertices = max_cell_vertices
    unfolded_grid%verts%max_connectivity   = max_vertex_connect

    CALL allocate_grid_object(unfolded_grid_id)
    CALL grid_set_exist_eq_allocated(unfolded_grid_id)

    !--------------------------------------------------------------
    !   get unfolded vertices-edges-cells
    DO x=0, x_no_of_columns-1
      DO y=0, y_no_of_rows-1
        ! vertices
        index_no=vertex_index(x,y)
        unfolded_index_no = unfolded_vertex_index(x,y)
        unfolded_grid%verts%idx(unfolded_index_no) = index_no
        unfolded_grid%verts%vertex(unfolded_index_no)%lon = torus_grid%verts%vertex(index_no)%lon
        unfolded_grid%verts%vertex(unfolded_index_no)%lat = torus_grid%verts%vertex(index_no)%lat
        ! edges
        index_no=edge_index_left_diagonal(x,y)
        unfolded_index_no = unfoldleftright_diag_edgeindex(x,y)
        unfolded_grid%edges%idx(unfolded_index_no) = index_no
        unfolded_grid%edges%center(unfolded_index_no)%lon = torus_grid%edges%center(index_no)%lon
        unfolded_grid%edges%center(unfolded_index_no)%lat = torus_grid%edges%center(index_no)%lat
        unfolded_grid%edges%get_vertex_index(unfolded_index_no,1) = &
          & unfolded_vertex_index(x,y)
        unfolded_grid%edges%get_vertex_index(unfolded_index_no,2) = &
          & unfolded_vertex_index(x,y+1)

        index_no=edge_index_right_diagonal(x,y)
        unfolded_index_no = unfoldright_diag_edgeindex(x,y)
        unfolded_grid%edges%idx(unfolded_index_no) = index_no
        unfolded_grid%edges%center(unfolded_index_no)%lon = torus_grid%edges%center(index_no)%lon
        unfolded_grid%edges%center(unfolded_index_no)%lat = torus_grid%edges%center(index_no)%lat
        unfolded_grid%edges%get_vertex_index(unfolded_index_no,1) = &
          & unfolded_vertex_index(x,y)
        unfolded_grid%edges%get_vertex_index(unfolded_index_no,2) = &
          & unfolded_vertex_index(x+1,y+1)

        index_no=edge_index_horizontal(x,y)
        unfolded_index_no = unfolded_horizontal_edge_index(x,y)
        unfolded_grid%edges%idx(unfolded_index_no) = index_no
        unfolded_grid%edges%center(unfolded_index_no)%lon = torus_grid%edges%center(index_no)%lon
        unfolded_grid%edges%center(unfolded_index_no)%lat = torus_grid%edges%center(index_no)%lat
        unfolded_grid%edges%get_vertex_index(unfolded_index_no,1) = &
          & unfolded_vertex_index(x,y+1)
        unfolded_grid%edges%get_vertex_index(unfolded_index_no,2) = &
          & unfolded_vertex_index(x,y)

        ! cells
        index_no=cell_index_top(x,y)
        unfolded_index_no=unfolded_cell_index_top(x,y)
        unfolded_grid%cells%idx(unfolded_index_no) = index_no
        unfolded_grid%cells%no_of_vertices(unfolded_index_no) = 3
        unfolded_grid%cells%get_vertex_index(unfolded_index_no,1) = unfolded_vertex_index(x,y)
        unfolded_grid%cells%get_vertex_index(unfolded_index_no,2) = unfolded_vertex_index(x,y+1)
        unfolded_grid%cells%get_vertex_index(unfolded_index_no,3) = unfolded_vertex_index(x+1,y+1)

        index_no=cell_index_top_right(x,y)
        unfolded_index_no=unfolded_cell_index_right(x,y)
        unfolded_grid%cells%idx(unfolded_index_no) = index_no
        unfolded_grid%cells%no_of_vertices(unfolded_index_no) = 3
        unfolded_grid%cells%get_vertex_index(unfolded_index_no,1) = unfolded_vertex_index(x,y)
        unfolded_grid%cells%get_vertex_index(unfolded_index_no,2) = unfolded_vertex_index(x+1,y+1)
        unfolded_grid%cells%get_vertex_index(unfolded_index_no,3) = unfolded_vertex_index(x+1,y)

      ENDDO
    ENDDO
    !--------------------------------------------------------------
    !   get cuts
    !--------------------------------------------------------------
    !   left_diagonal cut on the right
    x=x_no_of_columns
    DO y=0, y_no_of_rows-1
      ! vertices
      index_no=vertex_index(x,y)
      unfolded_index_no = unfolded_vertex_index(x,y)
      unfolded_grid%verts%idx(unfolded_index_no) = index_no
      x_lon_ref =  x_lon_start + REAL(x,wp) * x_lon_step - x_lon_step * 0.5_wp * REAL(y,wp)
      y_lat_ref=  y_lat_start + REAL(y,wp) * y_lat_step
      unfolded_grid%verts%vertex(unfolded_index_no)%lon = x_lon_ref
      unfolded_grid%verts%vertex(unfolded_index_no)%lat = y_lat_ref
      ! left_diagonal edges
      index_no=edge_index_left_diagonal(x,y)
      unfolded_index_no = unfoldleftright_diag_edgeindex(x,y)
      unfolded_grid%edges%idx(unfolded_index_no) = index_no
      unfolded_grid%edges%get_vertex_index(unfolded_index_no,1) = unfolded_vertex_index(x,y)
      unfolded_grid%edges%get_vertex_index(unfolded_index_no,2) = unfolded_vertex_index(x,y+1)
      !       WRITE(message_text,*) 'left_diagonal edge cut at:', x,y,unfoldedIndexNo,'->',indexNo
    ENDDO
    !   horizontal  cut on the top
    y=y_no_of_rows
    DO x=0, x_no_of_columns-1
      ! vertices
      index_no=vertex_index(x,y)
      unfolded_index_no = unfolded_vertex_index(x,y)
      unfolded_grid%verts%idx(unfolded_index_no) = index_no
      x_lon_ref =  x_lon_start + REAL(x,wp) * x_lon_step - x_lon_step * 0.5_wp * REAL(y,wp)
      y_lat_ref=  y_lat_start + REAL(y,wp) * y_lat_step
      unfolded_grid%verts%vertex(unfolded_index_no)%lon = x_lon_ref
      unfolded_grid%verts%vertex(unfolded_index_no)%lat = y_lat_ref
      ! horizontal edges
      index_no=edge_index_horizontal(x,y)
      unfolded_index_no = unfolded_horizontal_edge_index(x,y)
      unfolded_grid%edges%idx(unfolded_index_no) = index_no
      unfolded_grid%edges%get_vertex_index(unfolded_index_no,1) = unfolded_vertex_index(x+1,y)
      unfolded_grid%edges%get_vertex_index(unfolded_index_no,2) = unfolded_vertex_index(x,y)
      !        WRITE(message_text,*) 'Horizontal edge cut at:', x,y,unfoldedIndexNo,'->',indexNo
    ENDDO
    ! get top right vertex
    y=y_no_of_rows
    x=x_no_of_columns
    index_no=vertex_index(x,y)
    unfolded_index_no = unfolded_vertex_index(x,y)
    unfolded_grid%verts%idx(unfolded_index_no) = index_no
    x_lon_ref =  x_lon_start + REAL(x,wp) * x_lon_step - x_lon_step * 0.5_wp * REAL(y,wp)
    y_lat_ref=  y_lat_start + REAL(y,wp) * y_lat_step
    unfolded_grid%verts%vertex(unfolded_index_no)%lon = x_lon_ref
    unfolded_grid%verts%vertex(unfolded_index_no)%lat = y_lat_ref
    !  WRITE(message_text,*) 'Top Right Vertex:',unfoldedIndexNo,noOfUnfoldedVertices,'->',indexNo
    !--------------------------------------------------------------

    !--------------------------------------------------------------
  END SUBROUTINE create_unfolded_torus

  !--------------------------------------------------------------
  !  SUBROUTINE print_unfolded_torus()
  !>
  !! Prints the  unfolded torus grid to standard output
  SUBROUTINE print_unfolded_torus()

    INTEGER :: i, x, y, index_no, counter
    !!$  INTEGER :: noOfUnfoldedEdges
    !--------------------------------------------------------------
    ! print unfolded_grid to standard output
    !--------------------------------------------------------------
    ! print vertices
    WRITE(message_text,*) '-------------------------------'
    WRITE(message_text,*) '    ==  UNFOLDED GRID  ==='
    WRITE(message_text,*) '-------------------------------'
    DO x=0, x_no_of_columns
      DO y=0, y_no_of_rows
        index_no=unfolded_vertex_index(x,y)
        counter=unfolded_grid%verts%no_of_neigbors(index_no)
        WRITE(message_text,*) 'vertex:',index_no, '(',x,y,') ->',unfolded_grid%verts%idx(index_no)
      ENDDO
    ENDDO
    ! print left_diagonal edges
    DO x=0, x_no_of_columns
      DO y=0, y_no_of_rows-1
        index_no=unfoldleftright_diag_edgeindex(x,y)
        WRITE(message_text,*) &
          & 'vertical edge:',index_no, '(',x,y,') ->',unfolded_grid%edges%idx(index_no)
        WRITE(message_text,*) 'vertices', (unfolded_grid%edges%get_vertex_index(index_no,i),i=1,2)
      ENDDO
    ENDDO
    ! print right_diagonal edges
    DO x=0, x_no_of_columns-1
      DO y=0, y_no_of_rows-1
        index_no=unfoldright_diag_edgeindex(x,y)
        WRITE(message_text,*) &
          & 'diagonal edge:',index_no, '(',x,y,') ->',unfolded_grid%edges%idx(index_no)
        WRITE(message_text,*) 'vertices', (unfolded_grid%edges%get_vertex_index(index_no,i),i=1,2)
      ENDDO
    ENDDO
    ! print horizontal edges
    DO x=0, x_no_of_columns-1
      DO y=0, y_no_of_rows
        index_no=unfolded_horizontal_edge_index(x,y)
        WRITE(message_text,*) &
          & 'horizontal edge:',index_no, '(',x,y,') ->',unfolded_grid%edges%idx(index_no)
        WRITE(message_text,*) 'vertices', (unfolded_grid%edges%get_vertex_index(index_no,i),i=1,2)
      ENDDO
    ENDDO
    ! print cells
    DO x=0, x_no_of_columns-1
      DO y=0, y_no_of_rows-1
        index_no=unfolded_cell_index_top(x,y)
        WRITE(message_text,*) 'top triangle:',index_no, '(',x,y,')'
        WRITE(message_text,*) 'vertices', (unfolded_grid%cells%get_vertex_index(index_no,i),i=1,3)

        index_no=unfolded_cell_index_right(x,y)
        WRITE(message_text,*) 'right triangle:',index_no, '(',x,y,')'
        WRITE(message_text,*) 'vertices', (unfolded_grid%cells%get_vertex_index(index_no,i),i=1,3)
      ENDDO
    ENDDO
    ! print all edges
    DO index_no=1, no_of_unfolded_edges
      WRITE(message_text,*) 'unfolded edge:',index_no, ' ->',unfolded_grid%edges%idx(index_no)
      WRITE(message_text,*) 'vertices', (unfolded_grid%edges%get_vertex_index(index_no,i),i=1,2)
    ENDDO

    !--------------------------------------------------------------
  END SUBROUTINE print_unfolded_torus



END MODULE mo_create_torus_grid
