!-------------------------------------------------------------------------------------
! mo_grid_toolbox first implementation
!>
!! A collection of basic methods for the ICON grid
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
MODULE mo_grid_toolbox
#include "grid_definitions.inc"
  !-------------------------------------------------------------------------
  USE mo_kind,           ONLY: wp
  USE mo_io_units,       ONLY: filename_max
  USE mo_exception,      ONLY: message_text, message, finish
!   USE mo_impl_constants, ONLY: min_rlcell, max_rlcell,  &
!     & min_rledge, max_rledge,  &
!     & min_rlvert, max_rlvert
  USE mo_local_grid
  USE mo_io_local_grid,  ONLY: read_new_netcdf_grid,write_netcdf_grid
!  USE mo_local_grid_geometry,ONLY: geographical_to_cartesian

  IMPLICIT NONE
  !-------------------------------------------------------------------------
  PRIVATE
  ! !VERSION CONTROL:
  CHARACTER(LEN=*), PARAMETER :: version = '$Id$'


  PUBLIC :: until_convergence             ! constant that define max number od iteration
    ! for boundary smoothing
  PUBLIC :: smooth_boundaryfrom_cell_list ! smooths the boundary of a cell list
  PUBLIC :: get_grid_from_cell_list       ! creates a nea grid from a cell list
  PUBLIC :: get_conditional_list          ! gets a list of entities that satisfy
                                          ! a simple condition (>,=,<)
  PUBLIC :: get_boundary_edges            ! returns the boundary edges in a list
  PUBLIC :: get_grid_conditional_cells    ! Grid from a subset of grid cells, satisfying
    ! the given condition for drawing purposes.
  PUBLIC :: set_edge_elev_fromcells       !Calculates the edge elevation by linear
    ! interpolation of two cells elevation
  PUBLIC :: write_conditional_cells
  PUBLIC :: concatenate_grids
  PUBLIC :: concatenate_grid_files
  PUBLIC :: create_dual
  PUBLIC :: get_dual_grid, get_basic_dual_grid         !  returns the dual grid
  PUBLIC :: inverse_connectivity_verts
  PUBLIC :: shift_grid_ids
  PUBLIC :: add_to_list_if_not_exist
   
  INTEGER, PARAMETER :: until_convergence = 40


  INTERFACE get_conditional_list
    MODULE PROCEDURE get_float_conditional_list
    MODULE PROCEDURE get_int_conditional_list
  END INTERFACE
  INTERFACE get_grid_conditional_cells
    MODULE PROCEDURE get_grid_float_cond_cells
    MODULE PROCEDURE get_grid_int_cond_cells
  END INTERFACE
  !-------------------------------------------------------------------------

CONTAINS


  !-------------------------------------------------------------------------
  ! INTEGER FUNCTION get_grid_conditional_cells(grid_id, cells_check, condition, refer_value)
  !>
  !! Returns a grid containing a subset of grid cells, satisfying the given condition.
  INTEGER FUNCTION get_grid_int_cond_cells(grid_id, cells_check, condition, refer_value,&
    &  option_out_grid_id)
    INTEGER, INTENT(in) :: grid_id
    INTEGER, POINTER :: cells_check(:)
    INTEGER, INTENT(in) :: condition
    INTEGER, INTENT(in) :: refer_value
    INTEGER, OPTIONAL, INTENT(in) :: option_out_grid_id
    !-------------------------------------------------------------------------

    TYPE(t_grid),         POINTER :: in_grid
    TYPE(t_integer_list) :: cell_checklist, out_list

    !-------------------------------------------------------------------------
    in_grid => get_grid(grid_id)
    cell_checklist%list_size = in_grid%cells%no_of_existcells
    cell_checklist%value => cells_check
    CALL get_int_conditional_list(cell_checklist, condition, refer_value, -1, &
      & out_list)

    IF (PRESENT(option_out_grid_id)) THEN
      get_grid_int_cond_cells = get_grid_from_cell_list(grid_id, out_list,option_out_grid_id)
    ELSE
      get_grid_int_cond_cells = get_grid_from_cell_list(grid_id, out_list)
    ENDIF

  END FUNCTION get_grid_int_cond_cells
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  ! INTEGER FUNCTION get_grid_conditional_cells(grid_id, cells_check, condition, refer_value)
  !>
  !! Returns a grid containing a subset of grid cells, satisfying the given condition.
  INTEGER FUNCTION get_grid_float_cond_cells(grid_id, cells_check, condition, refer_value,&
    &  option_out_grid_id)
    INTEGER, INTENT(in) :: grid_id
    REAL(wp), POINTER :: cells_check(:)
    INTEGER, INTENT(in) :: condition
    REAL(wp), INTENT(in) :: refer_value
    INTEGER, OPTIONAL, INTENT(in) :: option_out_grid_id
    !-------------------------------------------------------------------------

    TYPE(t_grid),         POINTER :: in_grid
    TYPE(t_integer_list) :: out_list
    TYPE(t_float_list)   :: cell_checklist

    !-------------------------------------------------------------------------
    in_grid => get_grid(grid_id)
    cell_checklist%list_size = in_grid%cells%no_of_existcells
    cell_checklist%value => cells_check
    CALL get_float_conditional_list(cell_checklist, condition, refer_value, -1, &
      & out_list)

    IF (PRESENT(option_out_grid_id)) THEN
      get_grid_float_cond_cells = get_grid_from_cell_list(grid_id, out_list,option_out_grid_id)
    ELSE
      get_grid_float_cond_cells = get_grid_from_cell_list(grid_id, out_list)
    ENDIF

  END FUNCTION get_grid_float_cond_cells
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Returns in out_cell_list a smoothed subset of in_cell_list
  !! Smoothed cells have at most one boundary edge
  SUBROUTINE smooth_boundaryfrom_cell_list(in_grid_id, in_cell_list, out_cell_list, opt_iterations)
    INTEGER, INTENT(in) :: in_grid_id
    TYPE(t_integer_list), INTENT(in) :: in_cell_list
    TYPE(t_integer_list) :: out_cell_list
    INTEGER, OPTIONAL :: opt_iterations   ! not of iterations to be executed

    INTEGER :: iteration, max_iterations,in_cells, out_cells
    TYPE(t_integer_list) :: wrk_cell_list


    max_iterations = 1
    IF (PRESENT(opt_iterations)) &
      & max_iterations = opt_iterations

    ! smooth the boundary cells until converegnce
    in_cells = in_cell_list%list_size
    
    IF (in_cells < 1) THEN
      out_cell_list%list_size = 0
      NULLIFY(out_cell_list%value)
      RETURN
    ENDIF


    wrk_cell_list%list_size = in_cells
    wrk_cell_list%value => in_cell_list%value
  
    
    out_cells = 0
    iteration = 0
    DO WHILE (in_cells /= out_cells .AND. iteration <  max_iterations)
      NULLIFY(out_cell_list%value)
      iteration = iteration + 1
      CALL smooth_boundary_cell_list_prv(in_grid_id, wrk_cell_list, out_cell_list)
      in_cells  = wrk_cell_list%list_size
      out_cells = out_cell_list%list_size
      IF (iteration > 1) THEN
        DEALLOCATE(wrk_cell_list%value)
      ENDIF
      wrk_cell_list%list_size = out_cell_list%list_size
      wrk_cell_list%value => out_cell_list%value
    ENDDO

  END SUBROUTINE smooth_boundaryfrom_cell_list

  !-------------------------------------------------------------------------
  !   SUBROUTINE smooth_boundaryfrom_cell_list(in_grid_id, in_cell_list, out_cell_list)
  !>
  !! Returns in out_cell_list a smoothed subset of in_cell_list
  !! Smoothed cells have at most one boundary edge
  SUBROUTINE smooth_boundary_cell_list_prv(in_grid_id, in_cell_list, out_cell_list)
    INTEGER, INTENT(in) :: in_grid_id
    TYPE(t_integer_list), INTENT(in) :: in_cell_list
    TYPE(t_integer_list)    :: out_cell_list
    !-------------------------------------------------------------------------

    TYPE(t_grid),         POINTER :: in_grid
    TYPE(t_grid_cells),   POINTER :: in_cells
    TYPE(t_grid_edges),   POINTER :: in_edges
    INTEGER :: no_of_input_cells, no_of_input_edges
    INTEGER :: no_of_cell_edges, no_of_list_cells, no_of_tagged_edges, no_of_tagged_cells

    INTEGER, POINTER :: edge_tag(:)
    TYPE(t_integer_list) :: cell_tag

    INTEGER :: istat, i,j, edge_index, cell_index

    !-------------------------------------------------------------------------
    in_grid => get_grid(in_grid_id)
    in_cells=>in_grid%cells
    in_edges=>in_grid%edges
    no_of_input_cells = in_cells%no_of_existcells
    no_of_input_edges = in_edges%no_of_existedges

    WRITE(message_text,'(a)') "===================================="
    CALL message ('', TRIM(message_text))
    WRITE(message_text,'(a)') "smooth_boundaryfrom_cell_list"
    CALL message ('', TRIM(message_text))

    ! mark edges in the cell list
    ALLOCATE (edge_tag(0:no_of_input_edges), stat=istat)
    IF (istat >0) THEN
      CALL finish ('smooth_boundaryfrom_cell_list', 'Problem in allocating edgeTag')
    ENDIF
    edge_tag(:) = 0
    no_of_list_cells = in_cell_list%list_size
    no_of_cell_edges = in_cells%max_no_of_vertices
    DO i=1,no_of_list_cells
      cell_index = in_cell_list%value(i)
      DO j=1,no_of_cell_edges
        edge_tag(in_cells%get_edge_index(cell_index, j)) = &
          & edge_tag(in_cells%get_edge_index(cell_index, j)) + 1
      ENDDO
    ENDDO

    ! mark cells for smoothing
    ALLOCATE (cell_tag%value(no_of_input_cells), stat=istat)
    IF (istat >0) THEN
      CALL finish ('smooth_boundaryfrom_cell_list', 'Problem in allocating cellTag')
    ENDIF
    cell_tag%value(:) = 0
    cell_tag%list_size = no_of_input_cells

    no_of_tagged_cells = 0
    DO i=1,no_of_input_cells
      no_of_tagged_edges = 0
      DO j=1,no_of_cell_edges
        edge_index = in_cells%get_edge_index(i, j)
        IF (edge_index > 0 .and. edge_tag(edge_index) > 1) &
          & no_of_tagged_edges = no_of_tagged_edges+1
      ENDDO
      IF (no_of_tagged_edges > 1) THEN
        cell_tag%value(i) = 1
        no_of_tagged_cells = no_of_tagged_cells + 1
      ENDIF
    ENDDO

    DEALLOCATE(edge_tag)

    ! get the out_cell_list from the cell_tag
    CALL get_int_conditional_list(cell_tag, equal_to, 1, no_of_tagged_cells, out_cell_list)
    DEALLOCATE(cell_tag%value)

    WRITE(message_text,'(a, i9, i9)') "input/output Cells:", &
      & in_cell_list%list_size, out_cell_list%list_size
    CALL message ('', TRIM(message_text))

  END SUBROUTINE smooth_boundary_cell_list_prv
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !   SUBROUTINE get_int_conditional_list(tag_list, condition, tag, no_of_tagged_objects, out_list)
  !>
  !! Returns in out_list the subset of in_list that satisfies  the condition (>,=,<) to refer_value
  SUBROUTINE get_int_conditional_list(in_list, condition, refer_value, no_of_tagged_objects, &
      & out_list)
    TYPE(t_integer_list)  :: in_list
    INTEGER, INTENT(in) :: condition
    INTEGER, INTENT(in) :: refer_value
    INTEGER, INTENT(in) :: no_of_tagged_objects
    TYPE(t_integer_list)  :: out_list

    INTEGER :: in_list_size, out_list_size, allocate_size
    INTEGER :: i

    in_list_size = in_list%list_size
    IF (no_of_tagged_objects < 0) THEN
      allocate_size = in_list_size
    ELSE
      allocate_size = no_of_tagged_objects
    ENDIF

    ALLOCATE (out_list%value(allocate_size),stat=i)
    IF (i >0) THEN
      CALL finish ('get_tagged_list', 'Problem in allocating out_list')
    ENDIF

    out_list_size = 0
    SELECT CASE (condition)
    CASE (greater_than)
      DO i=1,in_list_size
        IF (in_list%value(i) > refer_value) THEN
          out_list_size = out_list_size + 1
          out_list%value(out_list_size) = i
        ENDIF
      ENDDO
    CASE (less_than)
      DO i=1,in_list_size
        IF (in_list%value(i) < refer_value) THEN
          out_list_size = out_list_size + 1
          out_list%value(out_list_size) = i
        ENDIF
      ENDDO
    CASE (equal_to)
      DO i=1,in_list_size
        IF (in_list%value(i) == refer_value) THEN
          out_list_size = out_list_size + 1
          out_list%value(out_list_size) = i
        ENDIF
      ENDDO
    CASE default
      CALL finish ('get_tagged_list', 'unrecognized condition')
    END SELECT

    IF (no_of_tagged_objects >= 0 .and. &
      & no_of_tagged_objects /= out_list_size) THEN
      ! WRITE(*,*) no_of_tagged_objects, out_list_size, in_list_size
      CALL finish ('get_tagged_list', 'no_of_tagged_objects /= out_list_size')
    ENDIF

    out_list%list_size = out_list_size

  END SUBROUTINE get_int_conditional_list
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !  SUBROUTINE get_float_conditional_list(in_list, condition, tag, no_of_tagged_objects, out_list)
  !>
  !! Returns in out_list the subset of in_list that satisfies  the condition (>,=,<) to refer_value
  !-------------------------------------------------------------------------
  SUBROUTINE get_float_conditional_list(in_list, condition, refer_value, no_of_tagged_objects, &
      & out_list)
    TYPE(t_float_list)  :: in_list
    INTEGER, INTENT(in) :: condition
    REAL(wp), INTENT(in) :: refer_value
    INTEGER, INTENT(in) :: no_of_tagged_objects
    TYPE(t_integer_list)  :: out_list

    INTEGER :: in_list_size, out_list_size, allocate_size
    INTEGER :: i

    in_list_size = in_list%list_size
    IF (no_of_tagged_objects < 0) THEN
      allocate_size = in_list_size
    ELSE
      allocate_size = no_of_tagged_objects
    ENDIF

    ALLOCATE (out_list%value(allocate_size),stat=i)
    IF (i >0) THEN
      CALL finish ('get_tagged_list', 'Problem in allocating out_list')
    ENDIF

    out_list_size = 0
    SELECT CASE (condition)
    CASE (greater_than)
      DO i=1,in_list_size
        IF (in_list%value(i) > refer_value) THEN
          out_list_size = out_list_size + 1
          out_list%value(out_list_size) = i
        ENDIF
      ENDDO
    CASE (less_than)
      DO i=1,in_list_size
        IF (in_list%value(i) < refer_value) THEN
          out_list_size = out_list_size + 1
          out_list%value(out_list_size) = i
        ENDIF
      ENDDO
    CASE (equal_to)
      DO i=1,in_list_size
        IF (in_list%value(i) == refer_value) THEN
          out_list_size = out_list_size + 1
          out_list%value(out_list_size) = i
        ENDIF
      ENDDO
    CASE default
      CALL finish ('get_tagged_list', 'unrecognized condition')
    END SELECT

    IF (no_of_tagged_objects >= 0 .and. &
      & no_of_tagged_objects /= out_list_size) THEN
      ! WRITE(*,*) no_of_tagged_objects, out_list_size, in_list_size
      CALL finish ('get_tagged_list', 'no_of_tagged_objects /= out_list_size')
    ENDIF

    out_list%list_size = out_list_size

  END SUBROUTINE get_float_conditional_list
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  !! Adds the value to the list, if it does not already exist
  INTEGER FUNCTION add_to_list_if_not_exist(inout_list, value)
    TYPE(t_integer_list), INTENT(inout)  :: inout_list
    INTEGER, INTENT(in) :: value

    INTEGER :: i
    
    CHARACTER(*), PARAMETER :: method_name = "add_to_list_if_not_exist"

    DO i=1,inout_list%list_size
      IF (inout_list%value(i) == value) THEN
        add_to_list_if_not_exist = i
        RETURN
      ENDIF
    ENDDO
    
    IF (inout_list%list_size >= inout_list%allocated_size) THEN
      CALL finish(method_name, "list_size >= allocated_size")
    ENDIF

    ! add the value to the list
    inout_list%list_size = inout_list%list_size + 1
    inout_list%value(inout_list%list_size) = value
    add_to_list_if_not_exist = inout_list%list_size
    
  END FUNCTION add_to_list_if_not_exist
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !  SUBROUTINE get_boundary_edges(grid_id, boundary_edge_list)
  !>
  ! Returns the boundary edges in boundary_edge_list
  SUBROUTINE get_boundary_edges(grid_id, boundary_edge_list)
    INTEGER, INTENT(in) :: grid_id
    TYPE(t_integer_list)  :: boundary_edge_list

    TYPE(t_grid),       POINTER :: current_grid
    TYPE(t_grid_edges), POINTER :: edges

    INTEGER :: no_of_edges, edge_index, no_of_boundary_edges
    INTEGER :: i_status


    current_grid    => get_grid(grid_id)
    edges           => current_grid%edges
    no_of_edges     =  edges%no_of_existedges

    ALLOCATE(boundary_edge_list%value(no_of_edges),stat=i_status)
    IF (i_status > 0) &
      & CALL finish ('get_boundary_edges', 'ALLOCATE(boundary_edge_list(noOfEdges))')
    no_of_boundary_edges    = 0

    ! mark boundary edges
    DO edge_index = 1, no_of_edges
      IF (edges%get_cell_index(edge_index,1) == 0 .or. edges%get_cell_index(edge_index,2) == 0 ) THEN
        ! add edge to boundary list
        no_of_boundary_edges = no_of_boundary_edges + 1
        boundary_edge_list%value(no_of_boundary_edges) = edge_index
      ENDIF
    ENDDO

    boundary_edge_list%list_size = no_of_boundary_edges
    WRITE(message_text,'(a,i9,a,i9)') "Total no of edges=",no_of_edges,&
      & " Boundary edges=", no_of_boundary_edges
    CALL message ('', TRIM(message_text))


  END SUBROUTINE get_boundary_edges
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !  SUBROUTINE set_edge_elev_fromcells(grid_id, interpolation_method)
  !>
  ! Calculates the edge elevation by interpolation_method
  ! based on the two neigboring cells elevation
  SUBROUTINE set_edge_elev_fromcells(grid_id, interpolation_method)
    INTEGER, INTENT(in) :: grid_id
    INTEGER, INTENT(in) :: interpolation_method

    TYPE(t_grid),       POINTER :: current_grid
    TYPE(t_grid_edges), POINTER :: edges
    TYPE(t_grid_cells), POINTER :: cells

    INTEGER :: no_of_edges, edge_index, cell_index1, cell_index2
    REAL(wp) :: lng1, lng2


    current_grid    => get_grid(grid_id)
    edges           => current_grid%edges
    no_of_edges     =  edges%no_of_existedges
    cells           => current_grid%cells

    ! mark boundary edges
    DO edge_index = 1, no_of_edges
      cell_index1 = edges%get_cell_index(edge_index,1)
      cell_index2 = edges%get_cell_index(edge_index,2)

      IF (cell_index1 == 0 .or. cell_index2 == 0 ) THEN
        edges%elevation(edge_index) = 0.0_wp
        CYCLE
      ENDIF

      SELECT CASE(interpolation_method)

      CASE(linear_interpolation)
        lng1 = edges%get_edge_cell_length(edge_index, 1)
        lng2 = edges%get_edge_cell_length(edge_index, 2)
        edges%elevation(edge_index) = &
           (lng2 * cells%elevation(cell_index1) + lng1 * cells%elevation(cell_index2)) &
           / (lng1 + lng2)

      CASE(min_interpolation)
        edges%elevation(edge_index) = &
        & MIN(cells%elevation(cell_index1), cells%elevation(cell_index2))

      CASE(max_interpolation)
        edges%elevation(edge_index) = &
        & MAX(cells%elevation(cell_index1), cells%elevation(cell_index2))

      CASE DEFAULT
        CALL finish("set_edge_elev_fromcells","unkown interpolation method")
      END SELECT

    ENDDO

  END SUBROUTINE set_edge_elev_fromcells
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !    INTEGER FUNCTION get_grid_from_cell_list(in_grid_id, cell_list)
  !>
  !! Returns the sub-grid of in_grid_id defined by the cell_list
  INTEGER FUNCTION get_grid_from_cell_list(in_grid_id, t_cell_list, option_out_grid_id)
    INTEGER, INTENT(in)  :: in_grid_id
    TYPE(t_integer_list), OPTIONAL   :: t_cell_list
    INTEGER, OPTIONAL, INTENT(in) :: option_out_grid_id

    INTEGER :: out_grid_id
    TYPE(t_grid),       POINTER :: in_grid
    TYPE(t_grid_cells), POINTER :: in_cells
    TYPE(t_grid_edges), POINTER :: in_edges
    TYPE(t_grid_vertices), POINTER :: in_verts
    INTEGER, ALLOCATABLE :: vertices_old_to_new_pointer(:), edges_old_to_new_pointer(:), &
      & cells_old_to_new_pointer(:)

    INTEGER :: no_of_list_cells, no_of_input_cells, no_of_input_edges, no_of_input_verts
    INTEGER :: no_of_output_cells, no_of_output_edges, no_of_output_verts
    INTEGER :: max_cell_vertices, max_vertex_connect
    INTEGER :: start_outcellno, start_outvertno, start_outedgeno
    INTEGER :: istat, i,j, cell_index, edge_index, vertex_index

    !-------------------------------------------------------------------------
    in_grid => get_grid(in_grid_id)

    in_cells=>in_grid%cells
    in_edges=>in_grid%edges
    in_verts=>in_grid%verts
    no_of_input_cells = in_cells%no_of_existcells
    no_of_input_edges = in_edges%no_of_existedges
    no_of_input_verts = in_verts%no_of_existvertices
    max_cell_vertices  = in_cells%max_no_of_vertices
    max_vertex_connect = in_verts%max_connectivity

    ! allocate and clear the oldToNewPointer arrays
    ALLOCATE (vertices_old_to_new_pointer(0:no_of_input_verts),&
      & edges_old_to_new_pointer(0:no_of_input_edges),&
      & cells_old_to_new_pointer(0:no_of_input_cells), stat=istat)
    IF (istat >0) THEN
      CALL finish ('get_grid_from_cell_list', 'Problem in allocating oldToNewPointer')
    ENDIF

    ! check if we are given the out_grid
!     IF (PRESENT(option_out_grid_id)) THEN
!       ! add to an existing grid
!       out_grid_id = option_out_grid_id
!       out_grid    => get_grid(out_grid_id)
!       start_outcellno = out_grid%cells%no_of_existcells
!       start_outedgeno = out_grid%edges%no_of_existedges
!       start_outvertno = out_grid%verts%no_of_existvertices
!     ELSE
!       start_outcellno = 0
!       start_outedgeno = 0
!       start_outvertno = 0
!     ENDIF

    start_outcellno = 0
    start_outedgeno = 0
    start_outvertno = 0
    vertices_old_to_new_pointer(:) = 0
    edges_old_to_new_pointer(:)    = 0
    cells_old_to_new_pointer(:)    = 0
    ! fill old to new pointers
    IF (PRESENT(t_cell_list)) THEN
      !--------------------------------------------------------------
      ! take cells from cell_list
      no_of_list_cells   = t_cell_list%list_size
      no_of_output_cells = start_outcellno + no_of_list_cells
      no_of_output_edges = start_outedgeno
      no_of_output_verts = start_outvertno
      DO i = 1, no_of_list_cells
        cell_index = t_cell_list%value(i)
        cells_old_to_new_pointer(cell_index) = i + start_outcellno
  !       IF (PRESENT(option_out_grid_id)) &
  !         print *, 'cells=',i,cell_index,'->', cells_old_to_new_pointer(cell_index)
        DO j = 1, max_cell_vertices

          edge_index = in_cells%get_edge_index(cell_index, j)
          IF ( edge_index > 0 .AND. edges_old_to_new_pointer(edge_index) == 0) THEN
            no_of_output_edges                  = no_of_output_edges + 1
            edges_old_to_new_pointer(edge_index) = no_of_output_edges
          ENDIF

          vertex_index = in_cells%get_vertex_index(cell_index, j)
          IF ( vertex_index > 0 .AND. vertices_old_to_new_pointer(vertex_index) == 0) THEN
            no_of_output_verts = no_of_output_verts + 1
            vertices_old_to_new_pointer(vertex_index) = no_of_output_verts
          ENDIF

        ENDDO ! j = 1, noOfCellVertices
      ENDDO ! i = 1, noOfOutputCells
      !--------------------------------------------------------------
    ELSE
      !--------------------------------------------------------------
      ! take all cells
      no_of_output_cells = start_outcellno + no_of_input_cells
      no_of_output_edges = start_outedgeno + no_of_input_edges
      no_of_output_verts = start_outvertno + no_of_input_verts
      DO i = 1, no_of_input_cells
        cells_old_to_new_pointer(i) = i + start_outcellno
      ENDDO
      DO i = 1, no_of_input_edges
        edges_old_to_new_pointer(i) = i + start_outedgeno
      ENDDO
      DO i = 1, no_of_input_verts
        vertices_old_to_new_pointer(i) = i + start_outvertno
      ENDDO
      !--------------------------------------------------------------

    ENDIF

    !--------------------------------------------------------------
    IF (.NOT.PRESENT(option_out_grid_id)) THEN
      out_grid_id = -1
      !  create out_grid
!       out_grid_id = new_grid()
!       out_grid    => get_grid(out_grid_id)
!       out_grid%ncells = no_of_output_cells
!       out_grid%nedges = no_of_output_edges
!       out_grid%nverts = no_of_output_verts
!       out_grid%cells%max_no_of_vertices = max_cell_vertices
!       out_grid%verts%max_connectivity   = max_vertex_connect
!       CALL allocate_grid_object(out_grid_id)
    ELSE
      out_grid_id = option_out_grid_id
!       ! check if we have enough space
!       IF (no_of_output_cells > out_grid%cells%no_of_allocatedcells .OR.&
!         & no_of_output_edges > out_grid%edges%no_of_allocatededges .OR. &
!         & no_of_output_verts > out_grid%verts%no_of_allocatedvertices) THEN
!         CALL finish("get_grid_from_cell_list", "Allocated grid memory is not enough")
!       ENDIF
!       ! update the number of entities
!       out_grid%cells%no_of_existcells    = no_of_output_cells
!       out_grid%edges%no_of_existedges    = no_of_output_edges
!       out_grid%verts%no_of_existvertices = no_of_output_verts
    ENDIF
    !--------------------------------------------------------------

    !--------------------------------------------------------------
    ! fill out_grid
    CALL copy_grid(in_grid_id, out_grid_id, &
      & vertices_old_to_new_pointer, edges_old_to_new_pointer, cells_old_to_new_pointer)

    !--------------------------------------------------------------
    CALL set_nest_defaultindexes(out_grid_id)
    CALL set_grid_creation(out_grid_id, cut_off_grid)
    CALL set_grid_parent_id(out_grid_id, in_grid_id)
    !--------------------------------------------------------------
    DEALLOCATE (vertices_old_to_new_pointer,&
      & edges_old_to_new_pointer, cells_old_to_new_pointer)

    get_grid_from_cell_list = out_grid_id
    !--------------------------------------------------------------
    WRITE(message_text,'(a)') "===================================="
    CALL message ('', TRIM(message_text))
    WRITE(message_text,'(a)') "get_grid_from_cell_list"
    CALL message ('', TRIM(message_text))
    WRITE(message_text,'(a, i9, i9,a, i9)') "input/output Cells:", no_of_input_cells,&
        no_of_output_cells, ' starting at ', start_outcellno

    CALL message ('', TRIM(message_text))
    WRITE(message_text,'(a)') "===================================="
    CALL message ('', TRIM(message_text))
    !--------------------------------------------------------------

  END FUNCTION get_grid_from_cell_list
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !    INTEGER FUNCTION concatenate_grids(grid_1_id, grid_2_id)
  !>
  !! Returns the concatenate grid of the two grids
  INTEGER FUNCTION concatenate_grids(in_grid_1_id, in_grid_2_id) result(out_grid_id)
    INTEGER, INTENT(in) :: in_grid_1_id, in_grid_2_id

    TYPE(t_grid), POINTER :: in_grid, out_grid
    INTEGER :: no_of_output_cells, no_of_output_edges, no_of_output_verts
    INTEGER :: max_cell_vertices, max_vertex_connect

    in_grid => get_grid(in_grid_1_id)
    no_of_output_cells = in_grid%cells%no_of_existcells
    no_of_output_edges = in_grid%edges%no_of_existedges
    no_of_output_verts = in_grid%verts%no_of_existvertices
    max_cell_vertices  = in_grid%cells%max_no_of_vertices
    max_vertex_connect = in_grid%verts%max_connectivity

    in_grid => get_grid(in_grid_2_id)
    no_of_output_cells = no_of_output_cells + in_grid%cells%no_of_existcells
    no_of_output_edges = no_of_output_edges + in_grid%edges%no_of_existedges
    no_of_output_verts = no_of_output_verts + in_grid%verts%no_of_existvertices
    max_cell_vertices  = max(max_cell_vertices,in_grid%cells%max_no_of_vertices)
    max_vertex_connect = max(max_vertex_connect,in_grid%verts%max_connectivity)


    out_grid_id = new_grid()
    out_grid    => get_grid(out_grid_id)
    out_grid%ncells = no_of_output_cells
    out_grid%nedges = no_of_output_edges
    out_grid%nverts = no_of_output_verts
    out_grid%cells%max_no_of_vertices = max_cell_vertices
    out_grid%verts%max_connectivity   = max_vertex_connect
    CALL allocate_grid_object(out_grid_id)

    CALL grid_set_parents_from(in_grid_1_id, parents_from_parentpointers)
    CALL grid_set_parents_from(in_grid_2_id, parents_from_parentpointers)
    CALL grid_set_parents_from(out_grid_id, parents_from_parentpointers)

    CALL copy_grid(in_grid_1_id, out_grid_id)
    CALL copy_grid(in_grid_2_id, out_grid_id)

END FUNCTION concatenate_grids
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE concatenate_grid_files(in_file_1, in_file_2, out_file)
    CHARACTER(LEN=filename_max), INTENT(in) :: in_file_1,in_file_2,out_file

    INTEGER :: in_grid_id_1, in_grid_id_2, out_grid_id

    in_grid_id_1 = read_new_netcdf_grid(in_file_1)
    in_grid_id_2 = read_new_netcdf_grid(in_file_2)

    out_grid_id = concatenate_grids(in_grid_id_1,in_grid_id_2)
    CALL write_netcdf_grid(out_grid_id, out_file)

    CALL delete_grid(in_grid_id_1)
    CALL delete_grid(in_grid_id_2)
    CALL delete_grid(out_grid_id)


  END SUBROUTINE concatenate_grid_files

 !-------------------------------------------------------------------------
  !   SUBROUTINE write_conditional_cells(grid_id, cells_tag, condition, refer_value, file_name)
  !>
  !! Writes a subset of grid cells, satisfying the given condition, for drawing purposes.
  SUBROUTINE write_conditional_cells(grid_id, cells_tag, condition, refer_value, file_name)
    INTEGER, INTENT(in) :: grid_id
    INTEGER, POINTER :: cells_tag(:)
    INTEGER, INTENT(in) :: condition
    INTEGER, INTENT(in) :: refer_value
    CHARACTER(LEN=filename_max), INTENT(in) :: file_name
    !-------------------------------------------------------------------------
    INTEGER :: draw_grid_id

    !-------------------------------------------------------------------------
    draw_grid_id = get_grid_conditional_cells(grid_id, cells_tag, condition, refer_value)
    CALL write_netcdf_grid(draw_grid_id, file_name)
    CALL delete_grid(draw_grid_id)

  END SUBROUTINE write_conditional_cells
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE create_dual(in_file, out_file)
    CHARACTER(LEN=filename_max), INTENT(in) :: in_file,out_file

    INTEGER :: in_grid_id, out_grid_id

    in_grid_id = read_new_netcdf_grid(in_file)
    out_grid_id = get_dual_grid(in_grid_id)
    CALL write_netcdf_grid(out_grid_id, out_file)

    CALL delete_grid(in_grid_id)
    CALL delete_grid(out_grid_id)

  END SUBROUTINE create_dual
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Constructs the basic dual grid of the input_grid.
  !! The constructed dual grid has only the basic info, so that the dual of the dual can be created.
  !! The info includes all connectivity, cartesian position of vertices and cell centers.
  !! Both grids should be 2 dimensional.
  !!
  !! Constructs the dual grid of the input_grid Both grids shoud be 2 dimentional
  !! The algebraic (topological) dual procedure maps:
  !!      vertices -> dual cells
  !!      edges -> dual edges
  !!      cells -> dual vertices
  !!
  !! TODO:
  !! Further work for accomodating the vertex indexing conventions is required.
  !!
  !! @par Revision History
  !! Developed  by  Leonidas Linardakis, MPI-M  (2009).
  !!
  INTEGER FUNCTION get_basic_dual_grid(input_grid_id) result(dual_grid_id)
    INTEGER, INTENT(in):: input_grid_id

    TYPE(t_grid), POINTER :: input_grid, dual_grid

    ! Local Variables
    INTEGER :: no_of_dual_cells,no_of_dual_edges,no_of_dual_vertices
    INTEGER :: dual_max_cell_vertices, dual_max_vertex_connect
    !-----------------------------------------------------------------------

    !--------------------------------------------------------------
    input_grid => get_grid(input_grid_id)
    no_of_dual_cells     = input_grid%nverts
    no_of_dual_edges     = input_grid%nedges
    no_of_dual_vertices  = input_grid%ncells
    dual_max_cell_vertices  = input_grid%verts%max_connectivity
    dual_max_vertex_connect = input_grid%cells%max_no_of_vertices
    !--------------------------------------------------------------

    dual_grid_id = new_grid()
    dual_grid    => get_grid(dual_grid_id)
    dual_grid%ncells = no_of_dual_cells
    dual_grid%nedges = no_of_dual_edges
    dual_grid%nverts = no_of_dual_vertices
    dual_grid%cells%max_no_of_vertices = dual_max_cell_vertices
    dual_grid%verts%max_connectivity   = dual_max_vertex_connect
    CALL allocate_grid_object(dual_grid_id)
    CALL grid_set_exist_eq_allocated(dual_grid_id)
    !--------------------------------------------------------------
    dual_grid%level = input_grid%level
    !--------------------------------------------------------------
    ! fill the dual edges
    dual_grid%edges%idx(1:no_of_dual_edges)                = &
      & input_grid%edges%idx(1:no_of_dual_edges)
    dual_grid%edges%get_cell_index(1:no_of_dual_edges,1)   = &
      & input_grid%edges%get_vertex_index(1:no_of_dual_edges,1)
    dual_grid%edges%get_cell_index(1:no_of_dual_edges,2)   = &
      & input_grid%edges%get_vertex_index(1:no_of_dual_edges,2)
    dual_grid%edges%get_vertex_index(1:no_of_dual_edges,1) = &
      & input_grid%edges%get_cell_index(1:no_of_dual_edges,1)
    dual_grid%edges%get_vertex_index(1:no_of_dual_edges,2) = &
      & input_grid%edges%get_cell_index(1:no_of_dual_edges,2)

    !--------------------------------------------------------------
    ! fill the dual vertices
    dual_grid%verts%idx(1:no_of_dual_vertices)            = &
      & input_grid%cells%idx(1:no_of_dual_vertices)
    dual_grid%verts%cartesian(1:no_of_dual_vertices)%x(1) = &
      & input_grid%cells%cartesian_center(1:no_of_dual_vertices)%x(1)
    dual_grid%verts%cartesian(1:no_of_dual_vertices)%x(2) = &
      & input_grid%cells%cartesian_center(1:no_of_dual_vertices)%x(2)
    dual_grid%verts%cartesian(1:no_of_dual_vertices)%x(3) = &
      & input_grid%cells%cartesian_center(1:no_of_dual_vertices)%x(3)
    dual_grid%verts%no_of_neigbors(1:no_of_dual_vertices) = &
      & input_grid%cells%no_of_vertices(1:no_of_dual_vertices)
    dual_grid%verts%get_neighbor_index(1:no_of_dual_vertices,:) = &
      & input_grid%cells%get_neighbor_index(1:no_of_dual_vertices,:)
    dual_grid%verts%get_cell_index(1:no_of_dual_vertices,:) = &
      & input_grid%cells%get_vertex_index(1:no_of_dual_vertices,:)
    dual_grid%verts%get_edge_index(1:no_of_dual_vertices,:) = &
      & input_grid%cells%get_edge_index(1:no_of_dual_vertices,:)

    !--------------------------------------------------------------
    ! fill the dual cells
!     call message('fill the dual cells','')
    dual_grid%cells%idx(1:no_of_dual_cells)                   = &
      & input_grid%verts%idx(1:no_of_dual_cells)
    dual_grid%cells%cartesian_center(1:no_of_dual_cells)%x(1) = &
      & input_grid%verts%cartesian(1:no_of_dual_cells)%x(1)
    dual_grid%cells%cartesian_center(1:no_of_dual_cells)%x(2) = &
      & input_grid%verts%cartesian(1:no_of_dual_cells)%x(2)
    dual_grid%cells%cartesian_center(1:no_of_dual_cells)%x(3) = &
      & input_grid%verts%cartesian(1:no_of_dual_cells)%x(3)
    dual_grid%cells%no_of_vertices(1:no_of_dual_cells)        = &
      & input_grid%verts%no_of_neigbors(1:no_of_dual_cells)
    dual_grid%cells%get_neighbor_index(1:no_of_dual_cells,:)    = &
      & input_grid%verts%get_neighbor_index(1:no_of_dual_cells,:)
    dual_grid%cells%get_vertex_index(1:no_of_dual_cells,:)      = &
      & input_grid%verts%get_cell_index(1:no_of_dual_cells,:)
    dual_grid%cells%get_edge_index(1:no_of_dual_cells,:)        = &
      & input_grid%verts%get_edge_index(1:no_of_dual_cells,:)

    CALL set_nest_defaultindexes(dual_grid_id)


  END FUNCTION get_basic_dual_grid
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Constructs the dual grid of the input_grid.
  !! Both grids should be 2 dimensional.
  !!
  !! Constructs the dual grid of the input_grid Both grids shoud be 2 dimentional
  !! The algebraic (topological) dual procedure maps:
  !!      vertices -> dual cells
  !!      edges -> dual edges
  !!      cells -> dual vertices
  !!
  !! Properties:
  !! Dual of a Voronoi diagram is a Delaunay triangulation
  !! Dual of a Delaunay triangulation is a Voronoi diagram
  !! NOTE: The center of Deleunay edges is the middle, but not for Voronoi edges
  !! Dual form: Edge centers are equidistant from cell centers in Voronoi diagrams,
  !!            not true for Delaunay triangles
  !!
  !! TODO:
  !! Further work for accomodating the vertex indexing conventions is required.
  !!
  !! @par Revision History
  !! Developed  by  Leonidas Linardakis, MPI-M  (2009).
  !!
  INTEGER FUNCTION get_dual_grid(input_grid_id) result(dual_grid_id)
    INTEGER, INTENT(in):: input_grid_id

    TYPE(t_grid), POINTER :: input_grid, dual_grid

    ! Local Variables
    INTEGER :: i,j,edge_index,direction
    INTEGER :: no_of_dual_cells,no_of_dual_edges,no_of_dual_vertices
    INTEGER :: dual_max_cell_vertices, dual_max_vertex_connect
    !-----------------------------------------------------------------------

    !--------------------------------------------------------------
    input_grid => get_grid(input_grid_id)
    !--------------------------------------------------------------

    dual_grid_id = get_basic_dual_grid(input_grid_id)
    dual_grid    => get_grid(dual_grid_id)
    no_of_dual_cells    = dual_grid%cells%no_of_existcells
    no_of_dual_edges    = dual_grid%edges%no_of_existedges
    no_of_dual_vertices = dual_grid%verts%no_of_existvertices
    dual_max_cell_vertices = dual_grid%cells%max_no_of_vertices
    dual_max_vertex_connect = dual_grid%verts%max_connectivity
    !--------------------------------------------------------------
    ! The basic info of the dual is filled by get_basic_dual_grid
    ! Fill remaining fields
    !--------------------------------------------------------------

    !--------------------------------------------------------------
    ! fill the dual edges
    dual_grid%edges%center(1:no_of_dual_edges)%lon         = &
      & input_grid%edges%center(1:no_of_dual_edges)%lon
    dual_grid%edges%center(1:no_of_dual_edges)%lat         = &
      & input_grid%edges%center(1:no_of_dual_edges)%lat
    
    ! NOTE: additional transport for the cartesian variables is required!
    ! This is not saved in file, we ignore it for the time
    dual_grid%edges%cartesian_center(:)%x(1) = input_grid%edges%cartesian_center(:)%x(1)
    dual_grid%edges%cartesian_center(:)%x(2) = input_grid%edges%cartesian_center(:)%x(2)
    dual_grid%edges%cartesian_center(:)%x(3) = input_grid%edges%cartesian_center(:)%x(3)
    
    dual_grid%edges%primal_edge_length(1:no_of_dual_edges) = &
      & input_grid%edges%dual_edge_length(1:no_of_dual_edges)
    dual_grid%edges%dual_edge_length(1:no_of_dual_edges)   = &
      & input_grid%edges%primal_edge_length(1:no_of_dual_edges)
    dual_grid%edges%get_edge_vert_length(1:no_of_dual_edges,:) = &
      & input_grid%edges%get_edge_cell_length(1:no_of_dual_edges,:)
    dual_grid%edges%get_edge_cell_length(1:no_of_dual_edges,:) = &
      & input_grid%edges%get_edge_vert_length(1:no_of_dual_edges,:)
    dual_grid%edges%dual_normal(1:no_of_dual_edges)%v1     = &
      & -input_grid%edges%primal_normal(1:no_of_dual_edges)%v1
    dual_grid%edges%dual_normal(1:no_of_dual_edges)%v2     = &
      & -input_grid%edges%primal_normal(1:no_of_dual_edges)%v2

    DO i = 1, no_of_dual_edges
      IF (input_grid%edges%system_orientation(i) == 1) THEN
        !         we keep the dual normal, inverse the prime normal to keep orientation
        dual_grid%edges%primal_normal(i)%v1   = input_grid%edges%dual_normal(i)%v1
        dual_grid%edges%primal_normal(i)%v2   = input_grid%edges%dual_normal(i)%v2
        dual_grid%edges%dual_normal(i)%v1     = -input_grid%edges%primal_normal(i)%v1
        dual_grid%edges%dual_normal(i)%v2     = -input_grid%edges%primal_normal(i)%v2
        dual_grid%edges%system_orientation(i)  = -1
      ELSE
        !         inverse the dual normal, keep the prime normal to keep orientation
        dual_grid%edges%primal_normal(i)%v1   = -input_grid%edges%dual_normal(i)%v1
        dual_grid%edges%primal_normal(i)%v2   = -input_grid%edges%dual_normal(i)%v2
        dual_grid%edges%dual_normal(i)%v1     = input_grid%edges%primal_normal(i)%v1
        dual_grid%edges%dual_normal(i)%v2     = input_grid%edges%primal_normal(i)%v2
        dual_grid%edges%system_orientation(i)= 1
      ENDIF

    ENDDO

    !--------------------------------------------------------------
    ! fill the dual vertices
    dual_grid%verts%dual_area(:)  = input_grid%cells%area(:)
    dual_grid%verts%vertex(:)%lon = input_grid%cells%center(:)%lon
    dual_grid%verts%vertex(:)%lat = input_grid%cells%center(:)%lat
    ! we do not store the cells%cartesian_center in file, so we have to recompute
    ! This is filled in get_basic_dual_grid
    !CALL (dual_grid%verts%vertex, no_of_dual_vertices,&
    !  & dual_grid%verts%cartesian)

    DO i = 1, no_of_dual_vertices
!      DO j = 1, dual_grid%verts%no_of_neigbors(i)
      DO j = 1, dual_max_vertex_connect
        edge_index                              = dual_grid%verts%get_edge_index(i,j)
        IF (edge_index > 0) THEN
          IF (i == dual_grid%edges%get_vertex_index(edge_index,1)) THEN
            direction = 1
          ELSE
            direction = -1
          ENDIF
          dual_grid%verts%get_edge_orient(i,j) = direction &
            & * dual_grid%edges%system_orientation(edge_index)
        ENDIF
      ENDDO
    ENDDO

    !--------------------------------------------------------------
    ! fill the dual cells
    dual_grid%cells%area(1:no_of_dual_cells)       = input_grid%verts%dual_area(:)
    dual_grid%cells%center(:)%lon = input_grid%verts%vertex(:)%lon
    dual_grid%cells%center(:)%lat = input_grid%verts%vertex(:)%lat

    DO i = 1, no_of_dual_cells
!      DO j = 1, dual_grid%cells%no_of_vertices(i)
      DO j = 1, dual_max_cell_vertices
        edge_index = dual_grid%cells%get_edge_index(i,j)
        IF (edge_index > 0) THEN
          IF (i == dual_grid%edges%get_cell_index(edge_index,1)) THEN
            direction = cell_direction_c1_c2
          ELSE
            direction = -cell_direction_c1_c2
          ENDIF
          dual_grid%cells%get_edge_orient(i,j) = direction ! +1 if outwards
        ENDIF
      ENDDO
    ENDDO

  END FUNCTION get_dual_grid
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  SUBROUTINE inverse_connectivity_verts(in_grid_id)
    INTEGER, INTENT(in) :: in_grid_id

    TYPE(t_grid), POINTER :: in_grid
!    TYPE(t_grid_cells), POINTER :: cells
    TYPE(t_grid_vertices), POINTER :: verts
    TYPE(t_grid_edges), POINTER :: edges

    INTEGER :: max_connectivity, vert_no, edge, i, j, no_of_conn_edges
    
    INTEGER, POINTER :: cell_index(:), edge_index(:), edge_orientation(:)


    CALL message("",'Inversing connectivity for vertices...')
    in_grid => get_grid(in_grid_id)
    verts => in_grid%verts
    edges => in_grid%edges
    max_connectivity = verts%max_connectivity
      
    ALLOCATE(cell_index(max_connectivity),&
      & edge_index(max_connectivity), edge_orientation(max_connectivity), stat=i)
    IF (i > 0) THEN
      CALL finish ('inverse_connectivity_verts', "ALLOCATE failed")
    ENDIF

    DO vert_no=1, verts%no_of_existvertices
      no_of_conn_edges=0
      DO i=1,max_connectivity      
        edge = verts%get_edge_index(vert_no, i)
        IF (edge /= 0) THEN
          no_of_conn_edges=no_of_conn_edges+1
          edge_index(no_of_conn_edges) = edge
          edge_orientation(no_of_conn_edges) = verts%get_edge_orient(vert_no,i)
          cell_index(no_of_conn_edges) = verts%get_cell_index(vert_no,i)
        ENDIF
!         cell = verts%get_cell_index(vert_no, i)                
      ENDDO !i=1,max_connectivity

      j=no_of_conn_edges
      DO i=1,no_of_conn_edges               
        verts%get_edge_index(vert_no, i) = edge_index(j)
        verts%get_neighbor_index(vert_no, i) = &
          & MERGE(edges%get_vertex_index(edge_index(j), 2),&
                & edges%get_vertex_index(edge_index(j), 1),&
                & edges%get_vertex_index(edge_index(j), 1) == vert_no)        
        verts%get_edge_orient(vert_no,i) = edge_orientation(j)
        
!         verts%get_cell_index(vert_no,i) = cell_index(i)
        
        j=j-1
      ENDDO !i=1,no_of_conn_edges
      DO i=no_of_conn_edges+1,max_connectivity
        verts%get_edge_index(vert_no, i) = 0
        verts%get_neighbor_index(vert_no, i) = 0       
        verts%get_edge_orient(vert_no,i) = 0
      ENDDO !i=no_of_conn_edges+1,max_connectivity
      
      verts%no_of_neigbors(vert_no) = no_of_conn_edges
      
    ENDDO !vert_no=1, verts%no_of_existvertices 

    DEALLOCATE(cell_index, edge_index, edge_orientation)
  
  END SUBROUTINE inverse_connectivity_verts
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Shifts the grid ids for this gird, the parent and children
  !! by the given value
  SUBROUTINE shift_grid_ids(file_name, shift_grid_id_value)
    CHARACTER(LEN=*), INTENT(in) :: file_name
    INTEGER, INTENT(in) :: shift_grid_id_value

    TYPE(t_grid), POINTER :: grid_obj
    TYPE(t_grid_cells), POINTER :: cells
    TYPE(t_grid_edges), POINTER :: edges
    TYPE(t_grid_vertices), POINTER :: verts
    INTEGER :: grid_id
    INTEGER :: no_of_cells, no_of_edges, no_of_verts
    INTEGER :: i

    grid_id = read_new_netcdf_grid(file_name, read_grid_ids=.true.)

    grid_obj => get_grid(grid_id)
    cells => grid_obj%cells
    edges => grid_obj%edges
    verts => grid_obj%verts
    no_of_cells = cells%no_of_existcells
    no_of_edges = edges%no_of_existedges
    no_of_verts = verts%no_of_existvertices

    grid_obj%patch_id = grid_obj%patch_id + shift_grid_id_value
    grid_obj%parent_grid_id = grid_obj%parent_grid_id + shift_grid_id_value
    ! The grid_obj%child_grid_ids(:) are not actually stored,
    ! skip it
    ! The verts%child_id(:) is not actually used.
    ! skip it
    DO i=1,no_of_edges    
      edges%child_id(i) = edges%child_id(i) + shift_grid_id_value
      edges%subgrid_id(i) = edges%subgrid_id(i) + shift_grid_id_value
    ENDDO
    DO i=1,no_of_cells
      cells%child_id(i) = cells%child_id(i) + shift_grid_id_value
      cells%subgrid_id(i) = cells%subgrid_id(i) + shift_grid_id_value
    ENDDO
    
    CALL write_netcdf_grid(grid_id, file_name)
    CALL delete_grid(grid_id)
    RETURN

  END SUBROUTINE shift_grid_ids
  !-------------------------------------------------------------------------
  
END MODULE mo_grid_toolbox

