!-------------------------------------------------------------------------------------
! mo_grid_conditions
!>
!! Set od geometric conditions for defining a sub-grid
!!
!! @author Leonidas Linardakis, MPI-M
!!
!! @par Revision History
!!   First implementation by Leonidas Linardakis, MPI-M, 2010-03-01
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
MODULE mo_grid_conditions
#include "grid_definitions.inc"
  !-------------------------------------------------------------------------
  USE mo_kind,               ONLY: wp
  USE mo_math_constants,     ONLY: rad2deg,pi
  USE mo_exception,          ONLY: message_text, message, finish
  USE mo_io_units,           ONLY: nnml, filename_max
  USE mo_namelist,           ONLY: position_nml, open_nml, positioned
  USE mo_local_grid
  USE mo_base_geometry,      ONLY: t_cartesian_coordinates, t_geographical_coordinates, gc2cc, &
    & arc_length
!  USE mo_local_grid_geometry,ONLY: geographical_to_cartesian
  USE mo_grid_toolbox,       ONLY: smooth_boundaryfrom_cell_list, &
    & get_grid_from_cell_list
  USE mo_io_local_grid,      ONLY: read_new_netcdf_grid, write_netcdf_grid
  USE mo_timer

  IMPLICIT NONE

  PRIVATE

  ! !VERSION CONTROL:
  CHARACTER(LEN=*), PARAMETER :: version = '$Id$'

  PUBLIC :: cut_local_grid, cut_local_grid_ascii
  PUBLIC :: cut_conditional_grid, read_grid_conditions
  PUBLIC :: get_conditional_cells
  PUBLIC :: get_boundary_vertices, get_inner_vertices
  !----------------------------------------

  ! !DEFINE PARAMETERS:
  INTEGER, PARAMETER ::  max_no_of_conditions = 40
  INTEGER, PARAMETER ::  rectangle_shape = 1
  INTEGER, PARAMETER ::  circle_shape = 2

  !-------------------------------------------------------------------------
  INTEGER :: no_of_conditions = 0
  REAL(wp) :: patch_center_x(max_no_of_conditions), patch_center_y(max_no_of_conditions)
  TYPE(t_geographical_coordinates) :: patch_center_geocoord(max_no_of_conditions)
  TYPE(t_cartesian_coordinates)    :: patch_center_cartesian(max_no_of_conditions)
  REAL(wp) :: rectangle_xradious(max_no_of_conditions), rectangle_yradious(max_no_of_conditions)
  REAL(wp) :: circle_radious(max_no_of_conditions)
  INTEGER :: patch_shape(max_no_of_conditions)
  INTEGER :: smooth_boundary_iterations

  CHARACTER(LEN=filename_max) :: input_file, output_file
  !-------------------------------------------------------------------------


CONTAINS

  !-------------------------------------------------------------------------
  !   FUNCTION read_grid_conditions(param_file_name)
  !>
  !! reads the parameters for the cutting or refining a grid. Private
  FUNCTION read_grid_conditions(param_file_name) result(no_of_read_conditions)

    CHARACTER(LEN=*), INTENT(in) :: param_file_name
    INTEGER :: no_of_read_conditions

    INTEGER :: i_status

    NAMELIST /grid_geometry_conditions/ input_file, output_file, &
      & no_of_conditions, patch_shape,   &
      & patch_center_x, patch_center_y,         &
      & rectangle_xradious, rectangle_yradious, &
      & circle_radious, smooth_boundary_iterations

    ! set default values
    input_file=''
    output_file = ''
    no_of_conditions = 0
    patch_shape=0
    rectangle_xradious = 0.0_WP
    rectangle_yradious = 0.0_WP
    circle_radious = 0.0_WP
    patch_center_x = 0.0_WP
    patch_center_y = 0.0_WP
    no_of_conditions = 0
    no_of_read_conditions = 0
    smooth_boundary_iterations = 1
    
    ! read namelist
    CALL open_nml(param_file_name)
    CALL position_nml('grid_geometry_conditions',STATUS=i_status)
    IF (i_status == positioned) THEN
      READ (nnml,grid_geometry_conditions)
    ELSE
      RETURN
!       WRITE(message_text,'(a,a)') " File", param_file_name, " not POSITIONED"
!       CALL finish ('read_grid_conditions', message_text)
    ENDIF
    CLOSE(nnml)

    WRITE(message_text,'(a)') "===================================="
    CALL message ('', TRIM(message_text))
    WRITE(message_text,'(a)') "----- grid_geometry_conditions -----"
    CALL message ('', TRIM(message_text))
    IF (input_file /= '') THEN
      WRITE(message_text,'(a,a)') "inputFile=",  TRIM(input_file)
      CALL message ('', TRIM(message_text))
      WRITE(message_text,'(a,a)') "outputFile=", TRIM(output_file)
      CALL message ('', TRIM(message_text))
    ENDIF
    WRITE(message_text,'(a,i2)') "no_of_conditions=", no_of_conditions
    CALL message ('', TRIM(message_text))
    WRITE(message_text,'(a,40i2)') "patch_shapes=", patch_shape(1:no_of_conditions)
    CALL message ('', TRIM(message_text))
    WRITE(message_text,'(a,40F6.2)') "patch_center_x=", patch_center_x(1:no_of_conditions)
    CALL message ('', TRIM(message_text))
    WRITE(message_text,'(a,40F6.2)') "patch_center_y=", patch_center_y(1:no_of_conditions)
    CALL message ('', TRIM(message_text))
    WRITE(message_text,'(a,40F6.2)') "rectangle_xradious=", &
      & rectangle_xradious(1:no_of_conditions)
    CALL message ('', TRIM(message_text))
    WRITE(message_text,'(a,40F6.2)') "rectangle_yradious=", &
      & rectangle_yradious(1:no_of_conditions)
    CALL message ('', TRIM(message_text))
    WRITE(message_text,'(a,40F16.3)') "circle_radious=", &
      & circle_radious(1:no_of_conditions)
    CALL message ('', TRIM(message_text))
    WRITE(message_text,'(a)') "===================================="
    CALL message ('', TRIM(message_text))

    CALL ref_set_param()
    no_of_read_conditions = no_of_conditions
    
  END FUNCTION read_grid_conditions
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !   SUBROUTINE cut_local_grid(param_file_name)
  !>
  !! Main routine for cutting-off a sub-grid. Public
  SUBROUTINE cut_local_grid(param_file_name)

    CHARACTER(LEN=*), INTENT(in) :: param_file_name

    INTEGER :: in_grid_id, out_grid_id, tmp

    tmp=read_grid_conditions(param_file_name)

    in_grid_id = read_new_netcdf_grid(input_file)

    out_grid_id = cut_conditional_grid(in_grid_id)
    IF (out_grid_id == in_grid_id .OR. out_grid_id < 1)&
      & CALL finish ('grid_cutLocalGrid', 'Failed to cut grid')

    CALL write_netcdf_grid(out_grid_id, output_file)
    CALL delete_grid(in_grid_id)
    CALL delete_grid(out_grid_id)

  END SUBROUTINE cut_local_grid
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !   SUBROUTINE cut_local_grid(param_file_name)
  !>
  !! Main routine for cutting-off a sub-grid. Public
  SUBROUTINE cut_local_grid_ascii(param_file_name)

   CHARACTER(LEN=*), INTENT(in) :: param_file_name

    TYPE(t_grid_cells), POINTER :: cells

    TYPE(t_integer_list)  :: cut_cell_list, smooth_cell_list
    INTEGER :: in_grid_id, i, k

    i = read_grid_conditions(param_file_name)

    in_grid_id = read_new_netcdf_grid(input_file)

    IF (no_of_conditions  < 1) RETURN

    CALL get_conditional_cells(in_grid_id, cut_cell_list)
    
    CALL smooth_boundaryfrom_cell_list(in_grid_id, cut_cell_list, smooth_cell_list, &
      & opt_iterations = smooth_boundary_iterations )

    cells => get_cells(in_grid_id)
    
    DO i = 1, smooth_cell_list%list_size
    
      k = smooth_cell_list%value(i)
      WRITE(0,*) "cell=", k, cells%center(k)%lon, cells%center(k)%lat
      
    ENDDO


    DEALLOCATE(cut_cell_list%value)
    DEALLOCATE(smooth_cell_list%value)

    CALL delete_grid(in_grid_id)

  END SUBROUTINE cut_local_grid_ascii
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Returns a list of boundary vertices with vertex_depth <= depth
  FUNCTION get_boundary_vertices(grid_id, depth) result(boundary_verts_list)
    INTEGER, INTENT(in) :: grid_id, depth
    TYPE(t_integer_list) :: boundary_verts_list

    TYPE(t_integer_list)  :: verts_depth_list
    INTEGER :: no_of_verts, vertex_index
    INTEGER :: error_status

    verts_depth_list = get_vertex_depths(grid_id,depth)

    no_of_verts = get_number_of_vertices(grid_id)
    ALLOCATE(boundary_verts_list%value(no_of_verts),stat=error_status)
    IF (error_status > 0) &
      & CALL finish ('get_boudary_edges', 'ALLOCATE(boundary_edges_list)')    
   
    boundary_verts_list%list_size = 0
    DO vertex_index=1,no_of_verts
      IF (verts_depth_list%value(vertex_index) >= 0 .AND. &
        & verts_depth_list%value(vertex_index) <= depth) THEN
        boundary_verts_list%list_size = boundary_verts_list%list_size +1
        boundary_verts_list%value(boundary_verts_list%list_size) = vertex_index
      ENDIF
    ENDDO
    
    DEALLOCATE(verts_depth_list%value)
    
  END FUNCTION get_boundary_vertices
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Returns a list of inner vertices with vertex_depth > depth
  FUNCTION get_inner_vertices(grid_id, depth) result(inner_verts_list)
    INTEGER, INTENT(in) :: grid_id, depth
    TYPE(t_integer_list) :: inner_verts_list

    TYPE(t_integer_list)  :: verts_depth_list
    INTEGER :: no_of_verts, vertex_index
    INTEGER :: error_status

    verts_depth_list = get_vertex_depths(grid_id,depth)

    no_of_verts = get_number_of_vertices(grid_id)
    ALLOCATE(inner_verts_list%value(no_of_verts),stat=error_status)
    IF (error_status > 0) &
      & CALL finish ('get_boudary_edges', 'ALLOCATE(boundary_edges_list)')

    inner_verts_list%list_size = 0
    DO vertex_index=1,no_of_verts
      IF (verts_depth_list%value(vertex_index) < 0 .OR. &
        & verts_depth_list%value(vertex_index) > depth) THEN
        inner_verts_list%list_size = inner_verts_list%list_size +1
        inner_verts_list%value(inner_verts_list%list_size) = vertex_index
      ENDIF
    ENDDO

    DEALLOCATE(verts_depth_list%value)

  END FUNCTION get_inner_vertices
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  !! Returns the depth of ecah vertex from the boundary
  !! 0=boundary
  !! -1 = vertex_depth is > depth
  !! if there are no boundaries then all vertexes get -1
  FUNCTION get_vertex_depths(grid_id, depth) result(verts_depth_list)
    INTEGER, INTENT(in) :: grid_id, depth
    TYPE(t_integer_list) :: verts_depth_list

    TYPE(t_grid),       POINTER :: current_grid
    TYPE(t_grid_edges), POINTER :: edges
    TYPE(t_grid_vertices), POINTER :: verts
    TYPE(t_integer_list)  :: boundary_edges_list
    INTEGER, POINTER :: vertex_stack(:)
    INTEGER :: stack_size
    
    INTEGER :: no_of_verts, vertex_index, stack_vertex_index, neigbor_idx
    INTEGER :: edge_index, edge_list_index, stack_vertex_depth
    INTEGER :: error_status

    IF (depth < 0) &
      & CALL finish ('get_vertex_depths', 'depth < 0')

    current_grid  => get_grid(grid_id)
    verts         => current_grid%verts
    edges         => current_grid%edges
    no_of_verts    = verts%no_of_existvertices

    ALLOCATE(verts_depth_list%value(no_of_verts),stat=error_status)
    IF (error_status > 0) &
      & CALL finish ('get_vertex_depths', 'ALLOCATE(verts_depth_list)')

    boundary_edges_list = get_boundary_edges(grid_id)

    !mark all vertices as inner
    verts_depth_list%value(1:no_of_verts) = -1

    IF (boundary_edges_list%list_size == 0) THEN
      ! all vertices are internal, retuen
      DEALLOCATE(boundary_edges_list%value)
      RETURN
    ENDIF

    ! make a stack
    ALLOCATE(vertex_stack(no_of_verts),stat=error_status)
    IF (error_status > 0) &
      & CALL finish ('get_vertex_depths', 'ALLOCATE(vertex_stack)')
    stack_size=0
    
    ! mark boundary vertices as 0 and add them to the stack
    DO edge_list_index = 1, boundary_edges_list%list_size
      edge_index = boundary_edges_list%value(edge_list_index)
      
      vertex_index =  edges%get_vertex_index(edge_index,1)
      IF (verts_depth_list%value(vertex_index) < 0) THEN
        stack_size=stack_size+1
        vertex_stack(stack_size) = vertex_index
        verts_depth_list%value(vertex_index) = 0
      ENDIF
      
      vertex_index =  edges%get_vertex_index(edge_index,2)
      IF (verts_depth_list%value(vertex_index) < 0) THEN
        stack_size=stack_size+1
        vertex_stack(stack_size) = vertex_index
        verts_depth_list%value(vertex_index) = 0
      ENDIF
   ENDDO

   ! while the stack is not empty
   !   pop a vertex
   !   if the vertex is < depth, process its neigbors
   DO WHILE(stack_size > 0)
     stack_vertex_index =  vertex_stack(stack_size)
     stack_size = stack_size - 1
     stack_vertex_depth = verts_depth_list%value(stack_vertex_index)
     IF (verts_depth_list%value(stack_vertex_index) < depth) THEN
        ! check its neigbors
        DO neigbor_idx=1, verts%max_connectivity
          vertex_index=verts%get_neighbor_index(stack_vertex_index,neigbor_idx)
          IF (vertex_index > 0) THEN
          
            IF (verts_depth_list%value(vertex_index) < 0) THEN
              stack_size=stack_size+1
              vertex_stack(stack_size) = vertex_index
              verts_depth_list%value(vertex_index) = stack_vertex_depth + 1
            ENDIF
            
          ENDIF          
        ENDDO 
     ENDIF

   ENDDO

   DEALLOCATE(vertex_stack)
    
  END FUNCTION get_vertex_depths
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !!
  FUNCTION get_boundary_edges(grid_id) result(boundary_edges_list)
    INTEGER, INTENT(in) :: grid_id
    TYPE(t_integer_list)  :: boundary_edges_list

    TYPE(t_grid),       POINTER :: current_grid
    TYPE(t_grid_edges), POINTER :: edges
    INTEGER :: no_of_edges, no_of_boundary_edges, edge_index
    INTEGER :: error_status

    current_grid  => get_grid(grid_id)
    edges         => current_grid%edges
    no_of_edges   = edges%no_of_existedges

    ALLOCATE(boundary_edges_list%value(no_of_edges),stat=error_status)
    IF (error_status > 0) &
      & CALL finish ('get_boudary_edges', 'ALLOCATE(boundary_edges_list)')
!     boundary_edges_list%value(:) = 0

    no_of_boundary_edges    = 0
    DO edge_index = 1, no_of_edges
      ! check if it is outer boundary
      IF (edges%get_cell_index(edge_index,1) == 0 .OR. &
        & edges%get_cell_index(edge_index,2) == 0) THEN
        ! add edge to boundary list
        no_of_boundary_edges = no_of_boundary_edges + 1
        boundary_edges_list%value(no_of_boundary_edges) = edge_index
      ENDIF
    ENDDO
    
    boundary_edges_list%list_size = no_of_boundary_edges

  END FUNCTION get_boundary_edges
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  !! Cuts-off the cut_grid_id from in_grid_id. Private
  INTEGER FUNCTION cut_conditional_grid(in_grid_id, param_file) result(cut_grid_id)
    INTEGER, INTENT(in) :: in_grid_id
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: param_file
!     INTEGER, INTENT(out):: cut_grid_id
!     INTEGER, INTENT(out):: cut_status

    TYPE(t_integer_list)  :: cut_cell_list, smooth_cell_list
    INTEGER :: tmp

!    cut_status = -
    IF (PRESENT(param_file)) tmp=read_grid_conditions(param_file)
    cut_grid_id = in_grid_id

    ! IF (patch_shape == undefined) RETURN
    IF (no_of_conditions  < 1) RETURN

    CALL get_conditional_cells(in_grid_id, cut_cell_list)
    CALL smooth_boundaryfrom_cell_list(in_grid_id, cut_cell_list, smooth_cell_list, &
      & opt_iterations = smooth_boundary_iterations )
    cut_grid_id = get_grid_from_cell_list(in_grid_id, smooth_cell_list)
    CALL set_grid_creation(cut_grid_id, cut_off_grid, from_grid_id=in_grid_id)

    DEALLOCATE(cut_cell_list%value)
    DEALLOCATE(smooth_cell_list%value)

    ! cut_status = 0

  END FUNCTION cut_conditional_grid
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !   SUBROUTINE ref_set_param
  !>
  !! Set the actual parameters for cutting and refining a grid. Private
  !! Changes degrees to rads, km to unit-sphere length
  SUBROUTINE ref_set_param()

    INTEGER :: i

    patch_center_x(1:no_of_conditions) = patch_center_x(1:no_of_conditions) / rad2deg
    patch_center_y(1:no_of_conditions) = patch_center_y(1:no_of_conditions) / rad2deg
    patch_center_geocoord(1:no_of_conditions)%lon = patch_center_x(1:no_of_conditions)
    patch_center_geocoord(1:no_of_conditions)%lat = patch_center_y(1:no_of_conditions)
    DO i=1,no_of_conditions
      patch_center_cartesian(i)    = gc2cc(patch_center_geocoord(i))
    ENDDO

    rectangle_xradious(1:no_of_conditions) = rectangle_xradious(1:no_of_conditions) / rad2deg
    rectangle_yradious(1:no_of_conditions) = rectangle_yradious(1:no_of_conditions) / rad2deg

    circle_radious(1:no_of_conditions) = circle_radious(1:no_of_conditions) / rad2deg

    WRITE(message_text,'(a,40F9.6)') "circle_radious=", circle_radious(1:no_of_conditions)
    CALL message ('', TRIM(message_text))

  END SUBROUTINE ref_set_param
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !   SUBROUTINE get_conditional_cells(in_grid_id, cell_list)
  !>
  !! Returns the cells in cell_list that are contained into an area. Private
  SUBROUTINE get_conditional_cells(in_grid_id, t_cell_list)

    INTEGER, INTENT(in):: in_grid_id
    TYPE(t_integer_list) :: t_cell_list

    TYPE(t_grid_cells), POINTER :: cells
    INTEGER :: no_of_input_cells,no_of_list_cells, cell_index
    ! TYPE(cartesian_coordinates), POINTER :: cartesian_cellcenter(:)
    TYPE(t_cartesian_coordinates) :: c_cellcenter
    REAL(wp) :: cell_lon,cell_lat,x_dist,y_dist

    INTEGER ::  istat, condition_index
    !-------------------------------------------------------------------------

    cells => get_cells(in_grid_id)
    no_of_input_cells = cells%no_of_allocatedcells

    ! cartesian_cellcenter => cells%cartesian_center
    !CALL geographical_to_cartesian(cells%center,no_of_input_cells, cartesian_cellcenter)

    ALLOCATE (t_cell_list%value(no_of_input_cells),stat=istat)
    IF (istat >0) THEN
      CALL finish ('get_conditional_cells', 'Problem in allocating cell_list')
    ENDIF

    no_of_list_cells=0
    DO cell_index = 1, no_of_input_cells
      cell_lon = cells%center(cell_index)%lon
      cell_lat = cells%center(cell_index)%lat
      c_cellcenter = cells%cartesian_center(cell_index)

      DO condition_index = 1, no_of_conditions

        SELECT CASE(patch_shape(condition_index))

        CASE(circle_shape)
          IF (arc_length(c_cellcenter, patch_center_cartesian(condition_index)) &
            & <= circle_radious(condition_index)) THEN
            ! cell center in circle
            ! add cell and exit loop
            no_of_list_cells=no_of_list_cells+1
            t_cell_list%value(no_of_list_cells) = cell_index
            EXIT
          ENDIF

        CASE(rectangle_shape)
          y_dist =  ABS(cell_lat - patch_center_y(condition_index))
          IF ( y_dist <= rectangle_yradious(condition_index)) THEN
            !  check the x distance
            x_dist = cell_lon - patch_center_x(condition_index)
            IF ( x_dist < -pi) THEN
              x_dist = x_dist + 2._wp*pi
            ELSE IF ( x_dist > pi) THEN
              x_dist = x_dist - 2._wp
            ENDIF
            IF (ABS(x_dist) <= rectangle_xradious(condition_index)) THEN
              ! cell center in rectangle
              ! add cell and exit loop
              no_of_list_cells=no_of_list_cells+1
              t_cell_list%value(no_of_list_cells) = cell_index
              EXIT
            ENDIF
          ENDIF

        CASE default
          CALL finish ('get_conditional_cells', 'Unkown patch_shape')

        END SELECT

      ENDDO !condition_index = 1, no_of_conditions
    ENDDO ! cell_index = 1, no_of_input_cells

    !  allocate cellList adjusted to the no_of_list_cells
    ! IF (ALLOCATED(cell_list)) DEALLOCATE (cell_list,stat=istat)

    ! fill cellList
    t_cell_list%list_size = no_of_list_cells
    ! cell_list%value(:) = tmp_cell_list%value(1:no_of_list_cells)

    ! clean-up
    ! DEALLOCATE (tmp_cell_list%value,stat=istat)
    RETURN

  END SUBROUTINE get_conditional_cells
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !   SUBROUTINE cells_center_in_rectangle(cells, tagged_cell_list)
  !>
  !! Returns in tagged_cell_list the cells whose center is in a rectangle. Private
  !   SUBROUTINE cells_center_in_rectangle(cells_id, tagged_cell_list)
  !     INTEGER, INTENT(in):: cells_id
  !     TYPE(integer_list) :: tagged_cell_list
  !
  !     TYPE(grid_cells), POINTER :: cells
  !     INTEGER :: no_of_input_cells,no_of_list_cells
  !     INTEGER :: cell_index,istat
  !     REAL(wp) :: cell_lon,cell_lat,x_dist,y_dist
  !
  !     cells => get_cells(cells_id)
  !     no_of_input_cells = cells%no_of_cells
  !
  !     WRITE(message_text,'(a)') "===================================="
  !     CALL message ('', TRIM(message_text))
  !     WRITE(message_text,'(a)')  "get cells in rectangle"
  !     CALL message ('', TRIM(message_text))
  !
  !     ALLOCATE (tagged_cell_list%value(no_of_input_cells),stat=istat)
  !     IF (istat >0) THEN
  !       CALL finish ('grid_getConditionalCellIndex', 'Problem in allocating tmpCellList')
  !     ENDIF
  !
  !     !  fill tagged_cell_list
  !     no_of_list_cells=0
  !     DO cell_index = 1, no_of_input_cells
  !       cell_lon = cells%center(cell_index)%lon
  !       cell_lat = cells%center(cell_index)%lat
  !
  !       y_dist =  ABS(cell_lat - patch_center_y)
  !       IF ( y_dist <= rectangle_yradious) THEN
  !         !  check the x distance
  !         x_dist = cell_lon - patch_center_x
  !         IF ( x_dist < -pi) THEN
  !           x_dist = x_dist + 2._wp*pi
  !         ELSE IF ( x_dist > pi) THEN
  !           x_dist = x_dist - 2._wp
  !         ENDIF
  !         IF (ABS(x_dist) <= rectangle_xradious) THEN
  !           ! cell center in rectangle
  !           no_of_list_cells=no_of_list_cells+1
  !           tagged_cell_list%value(no_of_list_cells) = cell_index
  !         ENDIF
  !       ENDIF
  !     ENDDO
  !
  !     tagged_cell_list%list_size = no_of_list_cells
  !     WRITE(message_text,'(a, i9, i9)') "input/output Cells:", no_of_input_cells, no_of_list_cells
  !     CALL message ('', TRIM(message_text))
  !
  !     RETURN
  !
  !   END SUBROUTINE cells_center_in_rectangle
  !   !-------------------------------------------------------------------------
  !
  !   !-------------------------------------------------------------------------
  !   !   SUBROUTINE cells_center_in_circle(cells, tagged_cell_list)
  !   !>
  !   !! Returns in tagged_cell_list the cells whose center is in a rectangle. Private
  !   SUBROUTINE cells_center_in_circle(cells_id, tagged_cell_list)
  !     INTEGER, INTENT(in):: cells_id
  !     TYPE(integer_list) :: tagged_cell_list
  !
  !     TYPE(grid_cells), POINTER :: cells
  !     INTEGER :: no_of_input_cells,no_of_list_cells
  !     TYPE(cartesian_coordinates), POINTER :: cartesian_cellcenter(:)
  !
  !     INTEGER :: cell_index,istat
  !
  !     cells => get_cells(cells_id)
  !     no_of_input_cells = cells%no_of_cells
  !
  !     WRITE(message_text,'(a)') "===================================="
  !     CALL message ('', TRIM(message_text))
  !     WRITE(message_text,'(a)')  "get cells in circle"
  !     CALL message ('', TRIM(message_text))
  !
  !     CALL geographical_to_cartesian(cells%center,no_of_input_cells, cartesian_cellcenter)
  !
  !     ALLOCATE (tagged_cell_list%value(no_of_input_cells),stat=istat)
  !     IF (istat >0) THEN
  !       CALL finish ('grid_getConditionalCellIndex', 'Problem in allocating tmpCellList')
  !     ENDIF
  !
  !     !  fill tagged_cell_list
  !     no_of_list_cells=0
  !     DO cell_index = 1, no_of_input_cells
  !       IF (arc_length(cartesian_cellcenter(cell_index), patch_center_cartesian) &
  !         & <= circle_radious) THEN
  !         ! cell center in circle
  !         no_of_list_cells=no_of_list_cells+1
  !         tagged_cell_list%value(no_of_list_cells) = cell_index
  !       ENDIF
  !     ENDDO
  !
  !     tagged_cell_list%list_size = no_of_list_cells
  !     WRITE(message_text,'(a, i9, i9)') "input/output Cells:", no_of_input_cells, no_of_list_cells
  !     CALL message ('', TRIM(message_text))
  !
  !     DEALLOCATE(cartesian_cellcenter)
  !
  !     RETURN
  !
  !   END SUBROUTINE cells_center_in_circle
  !-------------------------------------------------------------------------



END MODULE mo_grid_conditions

