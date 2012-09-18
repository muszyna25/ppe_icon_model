!-------------------------------------------------------------------------------------
! mo_localGridRefinement first implementation
!>
!! Refining of a triangular grid on the sphere by edge bisecting
!!
!! @author Leonidas Linardakis, MPI-M
!!
!! @par Revision History
!!   First implementation by Leonidas Linardakis, MPI-M, 2009-12-4
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
MODULE mo_local_grid_refinement
#include "grid_definitions.inc"
  !-------------------------------------------------------------------------
  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: message_text, message, finish
  USE mo_io_units,           ONLY: nnml, filename_max
  USE mo_namelist,           ONLY: position_nml, open_nml, positioned
  USE mo_timer,              ONLY: new_timer, timer_start, timer_stop, print_timer, delete_timer

  USE mo_base_geometry,      ONLY: sphere_cartesian_midpoint

  USE mo_local_grid
!   USE mo_local_grid,         ONLY: t_grid, t_grid_cells, t_grid_edges, t_grid_vertices,           &
!     & new_grid, delete_grid,get_grid, allocate_grid_object, refined_bisection_grid,       &
!     & parent_child_identical, parenttype_edge, parenttype_triangle, set_grid_creation,    &
!     & grid_set_exist_eq_allocated, set_nest_defaultindexes, grid_get_parent_pointers,     &
!     & set_grid_parent_id, set_no_of_subgrids, set_start_subgrids

  USE mo_grid_conditions,    ONLY: cut_conditional_grid
  USE mo_local_grid_geometry,ONLY: compute_sphere_grid_geometry, get_common_edge_vertex
!   USE mo_local_grid_hierarchy, ONLY: create_grid_hierarchy
  USE mo_local_grid_optimization, ONLY: read_grid_optimization_param, optimize_grid

  USE mo_io_local_grid,      ONLY: read_new_netcdf_grid, write_netcdf_grid

  IMPLICIT NONE

  PRIVATE

  ! !VERSION CONTROL:
  CHARACTER(LEN=*), PARAMETER :: version = '$Id$'

  PUBLIC :: grid_refine
  PUBLIC :: refine_grid_edgebisection
  PUBLIC :: refine_grid_insert_centers
  PUBLIC :: complete_grid_connectivity
  PUBLIC :: coarsen_grid_file
  !----------------------------------------

  INTEGER :: max_cell_vertices, max_vertex_connect

  !-------------------------------------------------------------------------
  ! namelist parameters
  CHARACTER(LEN=filename_max) :: input_file, output_file
  INTEGER :: refine_depth
  LOGICAL :: create_hierarchy
  !-------------------------------------------------------------------------


CONTAINS

  !-------------------------------------------------------------------------
  !   SUBROUTINE read_param(param_file_name)
  !>
  !! reads the parameters for the refining a grid. Private
  SUBROUTINE read_param(param_file_name)

    CHARACTER(LEN=*), INTENT(in) :: param_file_name

    INTEGER :: i_status

    NAMELIST /local_grid_refine/ input_file, output_file, &
      & refine_depth, create_hierarchy

    ! default values
    refine_depth=1
    create_hierarchy=.false.
    
    ! read namelist
    CALL open_nml(param_file_name)
    CALL position_nml('local_grid_refine',STATUS=i_status)
    IF (i_status == positioned) THEN
      READ (nnml,local_grid_refine)
    ELSE
      WRITE(message_text,'(a,a)') " File", param_file_name, " not POSITIONED"
      CALL finish ('read_param', message_text)
    ENDIF
    CLOSE(nnml)

    WRITE(message_text,'(a)') "===================================="
    CALL message ('', TRIM(message_text))
    WRITE(message_text,'(a)') "------ local grid refinement -------"
    CALL message ('', TRIM(message_text))
    WRITE(message_text,'(a,a)') "inputFile=",  TRIM(input_file)
    CALL message ('', TRIM(message_text))
    WRITE(message_text,'(a,a)') "outputFile=", TRIM(output_file)
    CALL message ('', TRIM(message_text))
    WRITE(message_text,'(a,i2.2)') "refine_depth=", refine_depth
    CALL message ('', TRIM(message_text))
    WRITE(message_text,'(a,L2)') "create_hierarchy=", create_hierarchy
    CALL message ('', TRIM(message_text))
    WRITE(message_text,'(a)') "===================================="
    CALL message ('', TRIM(message_text))

  END SUBROUTINE read_param
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !   SUBROUTINE cut_local_grid(param_file_name)
  !>
  !! Main routine for refining a grid by edge-bisection. Public
  SUBROUTINE grid_refine(param_file_name)

    CHARACTER(LEN=*), INTENT(in) :: param_file_name

    CHARACTER(LEN=filename_max) :: file_name
    INTEGER :: in_grid_id, out_grid_id, cutoff_grid_id
    INTEGER :: timer_total_grid_refine

    CALL read_param(param_file_name)

    in_grid_id = read_new_netcdf_grid(input_file)

!     CALL message ('grid_refine', 'start timer')
    timer_total_grid_refine = new_timer("total grid_refine")
    CALL timer_start(timer_total_grid_refine)

!     CALL message ('grid_refine', 'cut_conditional_grid...')
    cutoff_grid_id = cut_conditional_grid(in_grid_id, param_file_name)
    IF (cutoff_grid_id /= in_grid_id) THEN
      CALL delete_grid(in_grid_id)
    ENDIF

!     CALL message ('grid_refine', 'refine_grid_edgebisection...')
    out_grid_id = refine_grid_edgebisection(cutoff_grid_id)
    CALL delete_grid(cutoff_grid_id)
!     CALL create_grid_hierarchy(out_grid_id)

    CALL read_grid_optimization_param(param_file_name)
    CALL optimize_grid(out_grid_id)
    CALL compute_sphere_grid_geometry(out_grid_id)

    IF (refine_depth < 2) THEN
      WRITE(file_name,'(a)')  TRIM(output_file)
    ELSE
      WRITE(file_name,'(a,i2.2,a)')  TRIM(output_file), 1,  '.nc'
    ENDIF
    CALL write_netcdf_grid(out_grid_id, output_file)

    CALL timer_stop(timer_total_grid_refine)
    CALL print_timer(timer_total_grid_refine)
    CALL delete_timer(timer_total_grid_refine)

    CALL delete_grid(out_grid_id)

  END SUBROUTINE grid_refine
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !   SUBROUTINE cut_local_grid(param_file_name)
  !>
  !! Main routine for refining a grid by edge-bisection. Public
  SUBROUTINE coarsen_grid_file(child_filename, parent_filename, out_filename)

    CHARACTER(LEN=*), INTENT(in) :: child_filename, parent_filename, out_filename

    INTEGER :: child_grid_id, parent_grid_id

    child_grid_id  = read_new_netcdf_grid(child_filename)
    parent_grid_id = read_new_netcdf_grid(parent_filename)

!     CALL message ('grid_refine', 'start timer')
    
    CALL coarsen_child_parent_grid(child_grid_id, parent_grid_id)

    CALL write_netcdf_grid(parent_grid_id, out_filename)

    CALL delete_grid(child_grid_id)
    CALL delete_grid(parent_grid_id)

  END SUBROUTINE coarsen_grid_file
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Read-adjust the parent vertxes coordinates to the child's
  !! Effectivel carsens the child grid
  SUBROUTINE coarsen_child_parent_grid(child_grid_id, parent_grid_id)
    INTEGER, INTENT(in)  :: child_grid_id, parent_grid_id
    !-------------------------------------------------------------------------
!     TYPE(t_grid), POINTER :: parent_grid, child_grid

    TYPE(t_grid_vertices), POINTER :: parent_verts, child_verts

    INTEGER :: no_of_child_verts, no_of_parent_verts
    INTEGER :: child_vertex,parent_vertex
    INTEGER :: timer_coarsen

    timer_coarsen = new_timer("coarsen_child_parent_grid")
    CALL timer_start(timer_coarsen)

    parent_verts => get_vertices(parent_grid_id)
    child_verts  => get_vertices(child_grid_id)
    no_of_child_verts = child_verts%no_of_existvertices
    no_of_parent_verts = parent_verts%no_of_existvertices
!$OMP PARALLEL
!$OMP DO PRIVATE(child_vertex,parent_vertex)
    DO child_vertex = 1, no_of_child_verts
      parent_vertex = child_verts%parent_index(child_vertex)
      IF (parent_vertex > no_of_parent_verts)  &
         & CALL finish('coarsen_child_parent_grid',&
         & 'parent_vertex > no_of_parent_verts')
      IF (parent_vertex > 0) THEN
         parent_verts%cartesian(parent_vertex) = child_verts%cartesian(child_vertex)
      ENDIF
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

    CALL timer_stop(timer_coarsen)
    CALL print_timer(timer_coarsen)
    CALL delete_timer(timer_coarsen)

    CALL compute_sphere_grid_geometry(parent_grid_id)

  END SUBROUTINE coarsen_child_parent_grid
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
!   ELEMENTAL INTEGER FUNCTION newcellindex_cellvertex(cell, vertex_no_in_cell)
!     INTEGER, INTENT(in) :: cell, vertex_no_in_cell
!     newcellindex_cellvertex  = (((cell)-1) * 4 + (vertex_no_in_cell))
!   END FUNCTION
!
!   ELEMENTAL INTEGER FUNCTION newcellindex_cell(cell)
!     INTEGER, INTENT(in) :: cell
!     newcellindex_cell = ((cell) * 4)
!   END FUNCTION
!
!   ELEMENTAL INTEGER FUNCTION newvertexindex_vertex(vertex)
!     INTEGER, INTENT(in) :: vertex
!     newvertexindex_vertex = (vertex)
!   END FUNCTION
!
!   ELEMENTAL INTEGER FUNCTION newvertexindex_edgecenter(edge, no_of_input_verts)
!     INTEGER, INTENT(in) :: edge, no_of_input_verts
!     newvertexindex_edgecenter = (no_of_input_verts + (edge))
!   END FUNCTION
!
!   ELEMENTAL INTEGER FUNCTION newedgeindex_edgevertex(edge, vertex_in_edge)
!     INTEGER, INTENT(in) :: edge, vertex_in_edge
!     newedgeindex_edgevertex = (((edge)-1) * 2 + (vertex_in_edge))
!   END FUNCTION
!
!   ELEMENTAL INTEGER FUNCTION newedgeindex_cellvertex(cell, vertex_no_in_cell, no_of_input_edges)
!     INTEGER, INTENT(in) :: cell, vertex_no_in_cell, no_of_input_edges
!     newedgeindex_cellvertex = ((2 * no_of_input_edges) + ((cell)-1) * 3 + (vertex_no_in_cell))
!   END FUNCTION
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !   SUBROUTINE refine_grid_edgebisection(in_grid_id, out_grid_id)
  !>
  !! Refines the in_grid_id into out_grid_id by edge bisection. Public
  !! Note: It ony creates the topology, not the gerometry.
  !! The geometry should be created according to the underlying 2D surface (plane, sphere, etc)
  INTEGER FUNCTION refine_grid_edgebisection(in_grid_id) result(out_grid_id)
    INTEGER, INTENT(in):: in_grid_id

    !-------------------------------------------------------------------------
    TYPE(t_grid), POINTER :: in_grid, out_grid

    TYPE(t_grid_cells), POINTER :: in_cells, out_cells
    TYPE(t_grid_edges), POINTER :: in_edges, out_edges
    TYPE(t_grid_vertices), POINTER :: in_verts, out_verts
    INTEGER, POINTER :: parent_cell_pnt(:), parent_edge_pnt(:), parent_vertex_pnt(:)

    INTEGER :: no_of_input_cells, no_of_input_edges, no_of_input_verts
    INTEGER :: no_of_output_cells, no_of_output_edges, no_of_output_verts

    INTEGER :: out_cell_index, in_cell_index, out_cell1
    INTEGER :: out_edge_index, in_edge_index, in_edge1, in_edge2, in_edge3, opposite_edge
    INTEGER :: out_vertex_index, in_vertex_index, mid_vertex_index, vertex_no_in_cell
    INTEGER :: out_vertex1, out_vertex2

    INTEGER :: j,k, p1, p2

    REAL(wp) :: in_elevation
    INTEGER :: in_sea_land_mask

    ! reindexing arrays
    INTEGER, ALLOCATABLE :: newcellindex_cellvertex(:,:)
    INTEGER, ALLOCATABLE :: newcellindex_cell(:)
    INTEGER, ALLOCATABLE :: newvertexindex_vertex(:)
    INTEGER, ALLOCATABLE :: newvertexindex_edgecenter(:)
    INTEGER, ALLOCATABLE :: newedgeindex_edgevertex(:,:)
    INTEGER, ALLOCATABLE :: newedgeindex_cellvertex(:,:)


    CHARACTER(*), PARAMETER :: method_name="refine_grid_edgebisection"

    INTEGER :: timer_grid_refine

    !-------------------------------------------------------------------------
    timer_grid_refine = new_timer("refine_grid_edgebisection")
    CALL timer_start(timer_grid_refine)
    !-------------------------------------------------------------------------
    out_grid_id = new_grid()
    in_grid  => get_grid(in_grid_id)
    out_grid => get_grid(out_grid_id)
    !-------------------------------------------------------------------------
    !WARNING  these normally are read from the in_grid
    max_cell_vertices  = 3
    max_vertex_connect = 6

    !-------------------------------------------------------------------------
    in_verts=>in_grid%verts
    in_edges=>in_grid%edges
    in_cells=>in_grid%cells
    no_of_input_cells = in_cells%no_of_existcells
    no_of_input_edges = in_edges%no_of_existedges
    no_of_input_verts = in_verts%no_of_existvertices
    CALL grid_get_parent_pointers(in_grid_id, parent_cell_pnt, parent_edge_pnt, parent_vertex_pnt)

    !------------------------------------------------
    no_of_output_cells = 4 * no_of_input_cells
    no_of_output_edges = 2 * no_of_input_edges + 3 * no_of_input_cells
    no_of_output_verts = no_of_input_verts + no_of_input_edges

    !------------------------------------------------
    ! Allocate and fill the new indexes
    ALLOCATE(newcellindex_cellvertex(no_of_input_cells, max_cell_vertices),&
      & newcellindex_cell(no_of_input_cells),          &
      !
      & newvertexindex_vertex(no_of_input_verts),      &
      & newvertexindex_edgecenter(no_of_input_edges),  &
      !
      & newedgeindex_edgevertex(no_of_input_edges, 2), &
      & newedgeindex_cellvertex(no_of_input_cells, max_cell_vertices), &
      & stat=j)
    IF (j > 0) THEN
      CALL finish (method_name, 'ALLOCATE new_index')
    ENDIF

    ! zero new indexes
    newcellindex_cellvertex(:,:) = 0
    newcellindex_cell(:)         = 0
    newvertexindex_vertex(:)     = 0
    newvertexindex_edgecenter(:) = 0
    newedgeindex_edgevertex(:,:) = 0
    newedgeindex_cellvertex(:,:) = 0
    !fill the new  indexes
    out_vertex_index = 0
    out_edge_index   = 0
    out_cell_index   = 0
    in_edge3 = 0
    opposite_edge = 0
    DO in_cell_index=1, no_of_input_cells
      ! fill the center indexes
      out_cell_index = out_cell_index + 1
      newcellindex_cell(in_cell_index) = out_cell_index
      DO j=1, max_cell_vertices

        in_edge_index = in_cells%get_edge_index(in_cell_index, j)
        IF (in_edge_index < 1) &
          CALL finish(method_name,'in_edge_index < 1')
        ! fill middle vertices index
        IF (newvertexindex_edgecenter(in_edge_index) == 0) THEN
          out_vertex_index = out_vertex_index + 1
          newvertexindex_edgecenter(in_edge_index) =  out_vertex_index
        ENDIF

        !fill center edge index
        out_edge_index = out_edge_index + 1
        newedgeindex_cellvertex(in_cell_index,j) = out_edge_index
      ENDDO !j=1, max_cell_vertices

      ! fill the peripheral indexes
      DO j=1, max_cell_vertices
        !fill the peripheral cell indexes
        out_cell_index = out_cell_index + 1
        newcellindex_cellvertex(in_cell_index,j) = out_cell_index

        ! fill the peripheral vertex indexes
        in_vertex_index = in_cells%get_vertex_index(in_cell_index, j)
        IF (in_vertex_index < 1) &
          CALL finish(method_name,'in_vertex_index < 1')
        IF (newvertexindex_vertex(in_vertex_index) == 0) THEN
          out_vertex_index = out_vertex_index + 1
          newvertexindex_vertex(in_vertex_index) = out_vertex_index
        ENDIF

        ! go through the edges and find the two edges that have this vertex
        DO k=1,max_cell_vertices
          in_edge_index = in_cells%get_edge_index(in_cell_index, k)
          IF (in_edges%get_vertex_index(in_edge_index,1) == in_vertex_index) THEN
            IF (newedgeindex_edgevertex(in_edge_index, 1) == 0) THEN
              out_edge_index = out_edge_index + 1
              newedgeindex_edgevertex(in_edge_index, 1) = out_edge_index
            ENDIF
          ELSEIF (in_edges%get_vertex_index(in_edge_index,2) == in_vertex_index) THEN
            IF (newedgeindex_edgevertex(in_edge_index, 2) == 0) THEN
              out_edge_index = out_edge_index + 1
              newedgeindex_edgevertex(in_edge_index, 2) = out_edge_index
            ENDIF
          ENDIF
        ENDDO !k=1, max_cell_vertices

      ENDDO !j=1, max_cell_vertices

    ENDDO !in_cell_index=1, no_of_input_cells

    ! indexes are filled, do some checks
    IF (out_vertex_index /= no_of_output_verts) &
      CALL finish(method_name,'out_vertex_index /= no_of_output_verts')
    IF (out_edge_index /= no_of_output_edges) &
      CALL finish(method_name,'out_edge_index /= no_of_output_edges')
    IF (out_cell_index /= no_of_output_cells) &
      CALL finish(method_name,'out_cell_index /= no_of_output_cells')

    !------------------------------------------------
    WRITE(message_text,'(a)') "===================================="
    CALL message ('', TRIM(message_text))
    CALL message ('', TRIM(method_name))
    WRITE(message_text,'(a,i9,i9)') "no Of Input/Output Cells = ", &
      & no_of_input_cells, no_of_output_cells
    CALL message ('', TRIM(message_text))
    WRITE(message_text,'(a,i9,i9)') "                   Edges = ", &
      & no_of_input_edges, no_of_output_edges
    CALL message ('', TRIM(message_text))
    WRITE(message_text,'(a,i9,i9)') "                   Verts = ", &
      & no_of_input_verts, no_of_output_verts
    CALL message ('', TRIM(message_text))
    !------------------------------------------------


    ! #define newcellindex_cellvertex(cell, vertex_no_in_cell) (((cell)-1) * 4 + (vertex_no_in_cell))
    ! #define newcellindex_cell(cell) ((cell) * 4)
    ! #define newvertexindex_vertex(vertex) (vertex)
    ! #define newvertexindex_edgecenter(edge) (no_of_input_verts + (edge))
    ! #define newedgeindex_edgevertex(edge, vertex_in_edge) (((edge)-1) * 2 + (vertex_in_edge))
    ! #define newedgeindex_cellvertex(cell, vertex_no_in_cell) ((2 * no_of_input_edges) + ((cell)-1) * 3 + (vertex_no_in_cell))

    !--------------------------------------------------------------
    !  create out_grid
    out_grid%ncells = no_of_output_cells
    out_grid%nedges = no_of_output_edges
    out_grid%nverts = no_of_output_verts
    out_grid%cells%max_no_of_vertices = max_cell_vertices
    out_grid%verts%max_connectivity   = max_vertex_connect
    CALL allocate_grid_object(out_grid_id)
    CALL grid_set_exist_eq_allocated(out_grid_id)

    !--------------------------------------------------------------
    ! fill out_grid
    out_grid%level = in_grid%level+1
    out_cells=>out_grid%cells
    out_edges=>out_grid%edges
    out_verts=>out_grid%verts

    !--------------------------------------------------------------
    ! construct new vertices
    out_verts%no_of_neigbors(:) = 0
    DO in_vertex_index = 1, no_of_input_verts
      out_vertex_index = newvertexindex_vertex(in_vertex_index)
      out_verts%idx(out_vertex_index) = out_vertex_index
      out_verts%parent_index(out_vertex_index) = parent_vertex_pnt(in_vertex_index)
      out_verts%parent_child_type(out_vertex_index) = parent_child_identical
!      out_verts%vertex(out_vertex_index) = in_verts%vertex(in_vertex_index)
      out_verts%cartesian(out_vertex_index) = in_verts%cartesian(in_vertex_index)
    ENDDO
    DO in_edge_index = 1, no_of_input_edges
      out_vertex_index = newvertexindex_edgecenter(in_edge_index)
      out_verts%idx(out_vertex_index) = out_vertex_index

      ! negative values mean that parent is an edge
      out_verts%parent_index(out_vertex_index) = -parent_edge_pnt(in_edge_index)
      out_verts%parent_child_type(out_vertex_index) = parenttype_edge

      out_verts%cartesian(out_vertex_index) = &
        & sphere_cartesian_midpoint(&
          & in_verts%cartesian(in_edges%get_vertex_index(in_edge_index,1)), &
          & in_verts%cartesian(in_edges%get_vertex_index(in_edge_index,2)) )
    ENDDO

    !--------------------------------------------------------------
    ! construct new cells
    out_cells%no_of_domains(:) = in_cells%no_of_domains(:)
    DO in_cell_index = 1, no_of_input_cells
      in_elevation = in_cells%elevation(in_cell_index)
      in_sea_land_mask = in_cells%sea_land_mask(in_cell_index)

      !this is the new triangle in the center of the old
      out_cell1 = newcellindex_cell(in_cell_index)

      out_cells%idx(out_cell1) = out_cell1
      out_cells%parent_index(out_cell1) = parent_cell_pnt(in_cell_index)
!        print *, ' in_cell_index:', in_cell_index, &
!          & ' parent:', parent_cell_pnt(in_cell_index)
      out_cells%parent_child_type(out_cell1) = parenttype_triangle

      !inherit the values from the parent
      out_cells%elevation(out_cell1)     = in_elevation
      out_cells%sea_land_mask(out_cell1) = in_sea_land_mask
      out_cells%get_domain_id(:,out_cell1)= in_cells%get_domain_id(:,in_cell_index)

      ! get the vertices and the edges
      DO j=1,max_cell_vertices
        in_edge1 = in_cells%get_edge_index(in_cell_index,j)
        out_cells%get_vertex_index(out_cell1,j) = &
          & newvertexindex_edgecenter(in_edge1)
        out_cells%get_edge_index(out_cell1,j) = &
          & newedgeindex_cellvertex(in_cell_index, j)
      ENDDO


      ! construct outer cells
      DO j=1,max_cell_vertices

        in_edge1 = in_cells%get_edge_index(in_cell_index,j)
        SELECT CASE (j)
          CASE(1)
            in_edge2 = in_cells%get_edge_index(in_cell_index,2)
            in_edge3 = in_cells%get_edge_index(in_cell_index,3)
            opposite_edge = 3
          CASE(2)
            in_edge2 = in_cells%get_edge_index(in_cell_index,3)
            in_edge3 = in_cells%get_edge_index(in_cell_index,1)
            opposite_edge = 1
          CASE(3)
            in_edge2 = in_cells%get_edge_index(in_cell_index,1)
            in_edge3 = in_cells%get_edge_index(in_cell_index,2)
            opposite_edge = 2
        END SELECT
        in_vertex_index = get_common_edge_vertex(in_edges,in_edge1,in_edge2, p1, p2)

        vertex_no_in_cell = get_vertex_index_in_cell(in_cells, in_cell_index, in_vertex_index)
        out_cell_index = newcellindex_cellvertex(in_cell_index, vertex_no_in_cell)

        out_cells%idx(out_cell_index)              = out_cell_index
        out_cells%parent_index(out_cell_index)     = parent_cell_pnt(in_cell_index)
        out_cells%parent_child_type(out_cell_index)= parenttype_triangle + j

        !inherit the values from the parent
        out_cells%elevation(out_cell_index)        = in_elevation
        out_cells%sea_land_mask(out_cell_index)    = in_sea_land_mask
      
        out_cells%get_domain_id(:,out_cell_index)  = in_cells%get_domain_id(:,in_cell_index)

        ! fill the three vertices of the new triangle
        out_vertex1 = newvertexindex_edgecenter(in_edge1)
        out_vertex2 = newvertexindex_edgecenter(in_edge2)
        out_cells%get_vertex_index(out_cell_index,1) = out_vertex1
        out_cells%get_vertex_index(out_cell_index,2) = &
          & newvertexindex_vertex(in_vertex_index)
        out_cells%get_vertex_index(out_cell_index,3) = out_vertex2

        ! fill the three edges of the new triangle
        out_cells%get_edge_index(out_cell_index,1) = &
          & newedgeindex_edgevertex(in_edge1,p1)
        out_cells%get_edge_index(out_cell_index,2) = &
          & newedgeindex_edgevertex(in_edge2,p2)
        out_edge_index = newedgeindex_cellvertex(in_cell_index, vertex_no_in_cell)
        out_cells%get_edge_index(out_cell_index,3) = out_edge_index

        ! construct the inner edge
        out_edges%idx(out_edge_index)               = out_edge_index
        !  the parent will be the opposite (parallel) edge of the triangle
        out_edges%parent_index(out_edge_index)      = parent_edge_pnt(in_edge3)
        out_edges%parent_child_type(out_edge_index) = parenttype_triangle + opposite_edge
        out_edges%get_vertex_index(out_edge_index,1)= out_vertex1
        out_edges%get_vertex_index(out_edge_index,2)= out_vertex2


        ! fill the cell_index of the vertices
!         DO k=1,3
!           vertex_index = out_cells%get_vertex_index(out_cell_index,k)
!           no_of_neigbors = out_verts%no_of_neigbors(vertex_index) + 1
!           out_verts%no_of_neigbors(vertex_index) = no_of_neigbors
!           out_verts%get_cell_index(vertex_index,no_of_neigbors) = out_cell_index
!         ENDDO


      ENDDO !j=1,oldCells%noOfVertices(oldCellIndex)

      ! fill the cell_index of the vertices
!       DO k=1,3
!         vertex_index = out_cells%get_vertex_index(out_cell1,k)
!         no_of_neigbors = out_verts%no_of_neigbors(vertex_index) + 1
!         out_verts%no_of_neigbors(vertex_index) = no_of_neigbors
!         out_verts%get_cell_index(vertex_index,no_of_neigbors) = out_cell1
!       ENDDO

    ENDDO !oldCellIndex = 1, noOfInputCells

    !--------------------------------------------------------------
    ! construct new edges by splitting the existing edges
    DO in_edge_index = 1, no_of_input_edges
      mid_vertex_index = newvertexindex_edgecenter(in_edge_index)
      ! in_cell_index_list  = in_edges%get_cell_index(in_edge_index,:)

      DO j=1,2  ! construct the 2 new child edges
        in_vertex_index = in_edges%get_vertex_index(in_edge_index,j)
        out_edge_index = newedgeindex_edgevertex(in_edge_index,j)
        out_edges%idx(out_edge_index)              = out_edge_index
        out_edges%parent_index(out_edge_index)     = parent_edge_pnt(in_edge_index)
        out_edges%parent_child_type(out_edge_index) = parenttype_edge + j

        ! fill the 2 new vertices
        out_vertex_index = newvertexindex_vertex(in_vertex_index)
        IF (j == 1) THEN
          out_edges%get_vertex_index(out_edge_index,1) = out_vertex_index
          out_edges%get_vertex_index(out_edge_index,2) = mid_vertex_index
        ELSE
          out_edges%get_vertex_index(out_edge_index,2) = out_vertex_index
          out_edges%get_vertex_index(out_edge_index,1) = mid_vertex_index
        ENDIF
        ! fill the neighbor info in the vertices
!         no_of_neigbors = out_verts%no_of_neigbors(out_vertex_index) + 1
!         out_verts%no_of_neigbors(out_vertex_index) = no_of_neigbors
!         out_verts%get_edge_index(out_vertex_index,no_of_neigbors) = out_edge_index
!         out_verts%get_neighbor_index(out_vertex_index,no_of_neigbors) = mid_vertex_index
!
!         no_of_neigbors = out_verts%no_of_neigbors(mid_vertex_index) + 1
!         out_verts%no_of_neigbors(mid_vertex_index) = no_of_neigbors
!         out_verts%get_edge_index(mid_vertex_index,no_of_neigbors) = out_edge_index
!         out_verts%get_neighbor_index(mid_vertex_index,no_of_neigbors) = out_vertex_index

        ! fill the 2 new neighbor cells of the edge
!         DO k=1,2
!           IF (in_cell_index_list(k) > 0) THEN
!             vertex_no_in_cell  = &
!               & get_vertex_index_in_cell(in_cells, in_cell_index_list(k), in_vertex_index)
!             out_edges%get_cell_index(out_edge_index,k) = &
!               & newcellindex_cellvertex(in_cell_index_list(k), vertex_no_in_cell)
!           ELSE
!             out_edges%get_cell_index(out_edge_index,k) = 0
!           ENDIF
!         ENDDO ! k=1,2

        ! fill the edge and neighbor info of the cells
!         out_cell_index_list(:) = out_edges%get_cell_index(out_edge_index,:)
!         DO k=1,2
!           out_cell_index = out_cell_index_list(k)
!           IF (out_cell_index > 0) THEN
!             no_of_neigbors = out_cells%no_of_vertices(out_cell_index) + 1
!             out_cells%no_of_vertices(out_cell_index) = no_of_neigbors
!             out_cells%get_edge_index(out_cell_index, no_of_neigbors) = out_edge_index
!             IF (k == 1) THEN
!               out_cells%get_neighbor_index(out_cell_index, no_of_neigbors) = out_cell_index_list(2)
!             ELSE
!               out_cells%get_neighbor_index(out_cell_index, no_of_neigbors) = out_cell_index_list(1)
!             ENDIF
!           ENDIF
!        ENDDO ! k=1,2

      ENDDO ! j=1,2  ! construct the 2 new child edges
    ENDDO ! oldEdgeIndex = 1, noOfInputEdges

    ! second, create the edges inside the old triangles
!     DO in_cell_index = 1, no_of_input_cells
!       DO j=1,max_cell_vertices
!
!         in_edge1 = in_cells%get_edge_index(in_cell_index,j)
!         IF (j == max_cell_vertices) THEN
!           in_edge2 = in_cells%get_edge_index(in_cell_index,1)
!         ELSE
!           in_edge2 = in_cells%get_edge_index(in_cell_index,j+1)
!         ENDIF
!         IF (j == 1) THEN ! this works only for triangles
!           opposite_edge = 3
!         ELSE
!           opposite_edge = j-1
!         ENDIF
!         in_edge3 = in_cells%get_edge_index(in_cell_index,opposite_edge)
!         in_vertex_index = get_common_edge_vertex(in_edges,in_edge1,in_edge2)
!         vertex_no_in_cell = get_vertex_index_in_cell(in_cells, in_cell_index, in_vertex_index)
!         out_edge_index   = newedgeindex_cellvertex(in_cell_index, vertex_no_in_cell)
!
!         out_vertex1 = newvertexindex_edgecenter(in_edge1)
!         out_vertex2 = newvertexindex_edgecenter(in_edge2)
!
!         out_edges%idx(out_edge_index)              = out_edge_index
!         !  the parent will be the opposite (parallel) edge of the triangle
!         out_edges%parent_index(out_edge_index)     = parent_edge_pnt(in_edge3)
!         out_edges%parent_child_type(out_edge_index) = parenttype_triangle + opposite_edge
!         out_edges%get_vertex_index(out_edge_index,1)   = out_vertex1
!         out_edges%get_vertex_index(out_edge_index,2)   = out_vertex2
!
!         ! fill the neighbor info in the vertices
!         no_of_neigbors = out_verts%no_of_neigbors(out_vertex1) + 1
!         out_verts%no_of_neigbors(out_vertex1) = no_of_neigbors
!         out_verts%get_edge_index(out_vertex1,no_of_neigbors) = out_edge_index
!         out_verts%get_neighbor_index(out_vertex1,no_of_neigbors) = out_vertex2
!         no_of_neigbors = out_verts%no_of_neigbors(out_vertex2) + 1
!         out_verts%no_of_neigbors(out_vertex2) = no_of_neigbors
!         out_verts%get_edge_index(out_vertex2,no_of_neigbors) = out_edge_index
!         out_verts%get_neighbor_index(out_vertex2,no_of_neigbors) = out_vertex1
!
!         ! fill the 2 new neighbor cells of the edge
!         out_cell1 = newcellindex_cellvertex(in_cell_index, vertex_no_in_cell)
!         out_edges%get_cell_index(out_edge_index,1) = out_cell1
!         out_cell2 = newcellindex_cell(in_cell_index)
!         out_edges%get_cell_index(out_edge_index,2) = out_cell2
!
!         ! fill the edge and neighbor info of the cells
!         no_of_neigbors = out_cells%no_of_vertices(out_cell1) + 1
!         out_cells%no_of_vertices(out_cell1) = no_of_neigbors
!         out_cells%get_edge_index(out_cell1, no_of_neigbors) = out_edge_index
!         out_cells%get_neighbor_index(out_cell1, no_of_neigbors) = out_cell2
!
!         no_of_neigbors = out_cells%no_of_vertices(out_cell2) + 1
!         out_cells%no_of_vertices(out_cell2) = no_of_neigbors
!         out_cells%get_edge_index(out_cell2, no_of_neigbors) = out_edge_index
!         out_cells%get_neighbor_index(out_cell2, no_of_neigbors) = out_cell1
!
!       ENDDO ! j=1,oldCells%noOfVertices(oldCellIndex)
!     ENDDO ! oldCellIndex = 1, noOfInputCells

    DEALLOCATE(newcellindex_cellvertex,&
      & newcellindex_cell,          &
      & newvertexindex_vertex,      &
      & newvertexindex_edgecenter,  &
      & newedgeindex_edgevertex,    &
      & newedgeindex_cellvertex)

    !--------------------------------------------------------------
    CALL complete_grid_connectivity(out_grid_id)
    !--------------------------------------------------------------
    ! temporarly set  parameters to default
    CALL set_nest_defaultindexes(out_grid_id)
    CALL set_grid_parent_id(out_grid_id, in_grid_id)
    CALL set_grid_creation(out_grid_id, refined_bisection_grid, from_grid_id=in_grid_id)
    CALL set_no_of_subgrids(out_grid_id, 1)
    CALL set_start_subgrids(out_grid_id, 0)
    out_grid%cells%min_sea_depth = in_grid%cells%min_sea_depth
    out_grid%is_filled = .true.

    !--------------------------------------------------------------
    CALL timer_stop(timer_grid_refine)
    CALL print_timer(timer_grid_refine)
    CALL delete_timer(timer_grid_refine)
    !--------------------------------------------------------------
    RETURN

  END FUNCTION refine_grid_edgebisection
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !   INTEGER FUNCTION refine_grid_insert_centers(in_grid_id) result(out_grid_id)
  !>
  !! Refines the in_grid_id into a triangular grid by inserting the cell cenetrs
  !! and connecting them to the cell vertices. Public
  !! Note: It ony creates the topology, not the gerometry.
  !! The geometry should be created according to the underlying 2D surface (plane, sphere, etc)
  INTEGER FUNCTION refine_grid_insert_centers(in_grid_id) result(out_grid_id)
    INTEGER, INTENT(in):: in_grid_id

    !-------------------------------------------------------------------------
    TYPE(t_grid), POINTER :: in_grid, out_grid

    TYPE(t_grid_cells), POINTER :: in_cells, out_cells
    TYPE(t_grid_edges), POINTER :: in_edges, out_edges
    TYPE(t_grid_vertices), POINTER :: in_verts, out_verts

    INTEGER :: no_of_input_cells, no_of_input_edges, no_of_input_verts, max_in_cell_vertices
    INTEGER :: no_of_output_cells, no_of_output_edges, no_of_output_verts

    INTEGER :: out_cell_index, in_cell_index
    INTEGER :: out_edge_index, in_edge_index, out_edge1
    INTEGER :: out_vertex_index
    INTEGER :: in_vertex1, in_vertex2
    INTEGER :: out_vertex1, out_vertex2, cell_center
    INTEGER :: vertex1_in_cell, vertex2_in_cell
    INTEGER, ALLOCATABLE ::new_vertex_index(:), new_edge_index(:), new_cell_edge_index(:)

    INTEGER :: j, i
    INTEGER :: cell_creation_order(6) ! this will identify the order with which the cells
                                      ! will be created for each vertex (dual cell)

    !-------------------------------------------------------------------------
    in_grid  => get_grid(in_grid_id)
    out_grid_id = new_grid()
    out_grid => get_grid(out_grid_id)
    !-------------------------------------------------------------------------

    !-------------------------------------------------------------------------
    in_verts=>in_grid%verts
    in_edges=>in_grid%edges
    in_cells=>in_grid%cells
    no_of_input_cells = in_cells%no_of_existcells
    no_of_input_edges = in_edges%no_of_existedges
    no_of_input_verts = in_verts%no_of_existvertices
    max_in_cell_vertices  = in_cells%max_no_of_vertices

    ALLOCATE(new_vertex_index(no_of_input_verts),&
      & new_edge_index(no_of_input_edges),       &
      & new_cell_edge_index(max_in_cell_vertices),stat=j)
    IF (j > 0) THEN
      CALL finish ('refine_grid_insert_centers', 'ALLOCATE new_index')
    ENDIF
    new_vertex_index(:) = 0
    new_edge_index(:) = 0

    !------------------------------------------------
    ! note that in each iteration the number of cells are x3
    no_of_output_cells = 2 * no_of_input_edges
    no_of_output_edges = 3 * no_of_input_edges
    no_of_output_verts = no_of_input_cells + no_of_input_verts

    !------------------------------------------------
    WRITE(message_text,'(a)') "===================================="
    CALL message ('', TRIM(message_text))
    WRITE(message_text,'(a)')  "refine_grid_insert_centers"
    CALL message ('', TRIM(message_text))
    WRITE(message_text,'(a,i9,i9)') "no Of Input/Output Cells = ", &
      & no_of_input_cells, no_of_output_cells
    CALL message ('', TRIM(message_text))
    WRITE(message_text,'(a,i9,i9)') "                   Edges = ", &
      & no_of_input_edges, no_of_output_edges
    CALL message ('', TRIM(message_text))
    WRITE(message_text,'(a,i9,i9)') "                   Verts = ", &
      & no_of_input_verts, no_of_output_verts
    CALL message ('', TRIM(message_text))
    !------------------------------------------------

    !--------------------------------------------------------------
    !  create out_grid
    out_grid%ncells = no_of_output_cells
    out_grid%nedges = no_of_output_edges
    out_grid%nverts = no_of_output_verts
    out_grid%cells%max_no_of_vertices = 3
    out_grid%verts%max_connectivity   = MAX(2*in_grid%verts%max_connectivity, max_in_cell_vertices)
    CALL allocate_grid_object(out_grid_id)
    CALL grid_set_exist_eq_allocated(out_grid_id)
    !--------------------------------------------------------------
    ! fill out_grid
    out_grid%level = in_grid%level+1
    out_cells=>out_grid%cells
    out_edges=>out_grid%edges
    out_verts=>out_grid%verts
    out_cells%no_of_vertices(:) = 3

    !--------------------------------------------------------------
    ! construct new vertices,edges,cells
    out_vertex_index = 0
    out_edge_index = 0
    out_cell_index = 0
    out_cells%no_of_domains(:) = in_cells%no_of_domains(:)
    cell_creation_order(:) = (/ 2, 1, 3, 5, 4, 6/)
    DO in_cell_index=1,no_of_input_cells
      ! insert cell center
      out_vertex_index = out_vertex_index + 1
      out_verts%idx(out_vertex_index) = out_vertex_index
      out_verts%cartesian(out_vertex_index) = in_cells%cartesian_center(in_cell_index)
      cell_center = out_vertex_index

      new_cell_edge_index(:) = 0
      ! go through the cell vertices and create a new vertex and edge for each vertex
      ! use the connectivity matrix
      DO i = 1, max_in_cell_vertices
        j = cell_creation_order(i)
        in_vertex1 = in_cells%get_vertex_index(in_cell_index, j)
        IF (in_vertex1 /= 0) THEN
          out_vertex1 = new_vertex_index(in_vertex1)
          IF (out_vertex1 == 0) THEN
            ! insert this vertex
            out_vertex_index = out_vertex_index + 1
            out_verts%idx(out_vertex_index) = out_vertex_index
            out_verts%cartesian(out_vertex_index) = in_verts%cartesian(in_vertex1)
            new_vertex_index(in_vertex1) = out_vertex_index
            out_vertex1 = out_vertex_index
          ENDIF

          out_edge_index = out_edge_index + 1
          out_edges%idx(out_edge_index) = out_edge_index
          out_edges%get_vertex_index(out_edge_index,1) = cell_center
          out_edges%get_vertex_index(out_edge_index,2) = out_vertex1
          new_cell_edge_index(j) = out_edge_index
        ENDIF
      ENDDO

      ! go through the cell edges and create a new triagle for each edge
      DO j=1,max_in_cell_vertices
        in_edge_index = in_cells%get_edge_index(in_cell_index, j)
        IF (in_edge_index /= 0) THEN
          !get the two new vertices
          in_vertex1 = in_edges%get_vertex_index(in_edge_index,1)
          in_vertex2 = in_edges%get_vertex_index(in_edge_index,2)
          out_vertex1 = new_vertex_index(in_vertex1)
          out_vertex2 = new_vertex_index(in_vertex2)

          !get the new edge corresponding to in_edge_index
          out_edge1 = new_edge_index(in_edge_index)
          IF (out_edge1 == 0) THEN
            ! insert this edge
            out_edge_index = out_edge_index + 1
            out_edges%idx(out_edge_index) = out_edge_index
            out_edges%get_vertex_index(out_edge_index,1) = out_vertex1
            out_edges%get_vertex_index(out_edge_index,2) = out_vertex2
            out_edge1 = out_edge_index
            new_edge_index(in_edge_index) = out_edge_index
          ENDIF

          ! construct the two new edges
          ! first find the location of the two vertices in the old cell
          DO vertex1_in_cell=1,max_in_cell_vertices
            IF (in_vertex1 == in_cells%get_vertex_index(in_cell_index, vertex1_in_cell)) EXIT
          ENDDO
          DO vertex2_in_cell=1,max_in_cell_vertices
            IF (in_vertex2 == in_cells%get_vertex_index(in_cell_index, vertex2_in_cell)) EXIT
          ENDDO
          IF (in_cells%get_vertex_index(in_cell_index, vertex1_in_cell) /=  in_vertex1 .OR. &
              in_cells%get_vertex_index(in_cell_index, vertex2_in_cell) /=  in_vertex2) &
                CALL finish('refine_grid_insert_centers', 'vertex not found in cell')

          ! create a new triangle and add the three edges
          out_cell_index = out_cell_index + 1
          out_cells%idx(out_cell_index) = out_cell_index
          out_cells%get_edge_index(out_cell_index,1) = out_edge1
          out_cells%get_edge_index(out_cell_index,2) = new_cell_edge_index(vertex2_in_cell)
          out_cells%get_edge_index(out_cell_index,3) = new_cell_edge_index(vertex1_in_cell)
          out_cells%get_vertex_index(out_cell_index,1) = out_vertex2
          out_cells%get_vertex_index(out_cell_index,2) = cell_center
          out_cells%get_vertex_index(out_cell_index,3) = out_vertex1

          out_cells%get_domain_id(:,out_cell_index) = in_cells%get_domain_id(:,in_cell_index)

          ! WRITE(0,*) 'cell edges:', out_cells%get_edge_index(out_cell_index,:)
        ENDIF
      ENDDO
    ENDDO

    ! some checks
    IF (out_vertex_index /= no_of_output_verts) THEN
      CALL finish('refine_grid_insert_centers','out_vertex_index /= no_of_output_verts')
    ENDIF
    IF (out_edge_index /= no_of_output_edges) THEN
      WRITE(0,*) 'out_edge_index=',out_edge_index
      CALL finish('refine_grid_insert_centers','out_edge_index /= no_of_output_edges')
    ENDIF
    IF (out_cell_index /= no_of_output_cells) THEN
      CALL finish('refine_grid_insert_centers','out_cell_index /= no_of_output_cells')
    ENDIF

    DEALLOCATE(new_vertex_index, new_edge_index, new_cell_edge_index)
    ! complete the connectivity
    CALL complete_grid_connectivity(out_grid_id)

    CALL set_nest_defaultindexes(out_grid_id)
    !CALL set_grid_parent_id(out_grid_id, in_grid_id)
    CALL set_grid_creation(out_grid_id, dualy_refined_grid, from_grid_id=in_grid_id)
    CALL set_no_of_subgrids(out_grid_id, 1)
    CALL set_start_subgrids(out_grid_id, 0)
    out_grid%cells%min_sea_depth = in_grid%cells%min_sea_depth
    out_grid%is_filled = .true.
    !--------------------------------------------------------------

  END FUNCTION refine_grid_insert_centers
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !   SUBROUTINE complete_grid_connectivity(in_grid_id)
  !>
  !! Completes the connectivity fields of a grid. Private
  !! Assumes that the edges have the vertex indices, and cells have the vertex and edge indexes.
  !! It fills the rest of the connectivity fields.
  SUBROUTINE complete_grid_connectivity(in_grid_id)
    INTEGER, INTENT(in):: in_grid_id

    !-------------------------------------------------------------------------
    TYPE(t_grid), POINTER :: in_grid

    TYPE(t_grid_cells), POINTER :: in_cells
    TYPE(t_grid_edges), POINTER :: in_edges
    TYPE(t_grid_vertices), POINTER :: in_verts

    INTEGER :: no_of_input_cells, no_of_input_edges, no_of_input_verts, max_in_cell_vertices

    INTEGER :: in_cell_index, in_edge_index, in_vertex_index
    INTEGER :: in_vertex1, in_vertex2

    INTEGER :: j

    !CALL message ('complete_grid_connectivity','starts...')
    !-------------------------------------------------------------------------
    in_grid  => get_grid(in_grid_id)
    in_verts=>in_grid%verts
    in_edges=>in_grid%edges
    in_cells=>in_grid%cells
    no_of_input_cells = in_cells%no_of_existcells
    no_of_input_edges = in_edges%no_of_existedges
    no_of_input_verts = in_verts%no_of_existvertices
    max_in_cell_vertices  = in_cells%max_no_of_vertices

    !clean the connectivity fileds to be filled
    in_verts%get_neighbor_index(:,:) = 0
    in_verts%get_edge_index(:,:) = 0
    in_verts%get_cell_index(:,:) = 0

    in_edges%get_cell_index(:,:) = 0

    in_cells%get_neighbor_index(:,:) = 0

    ! complete edge-vertex connectivity
    in_verts%no_of_neigbors(:) = 0
    DO in_edge_index=1,no_of_input_edges
      ! complete vertex connectivity
      in_vertex1 = in_edges%get_vertex_index(in_edge_index,1)
      in_vertex2 = in_edges%get_vertex_index(in_edge_index,2)
      in_verts%no_of_neigbors(in_vertex1) = in_verts%no_of_neigbors(in_vertex1) + 1
      in_verts%no_of_neigbors(in_vertex2) = in_verts%no_of_neigbors(in_vertex2) + 1
      in_verts%get_edge_index(in_vertex1, in_verts%no_of_neigbors(in_vertex1)) = in_edge_index
      in_verts%get_edge_index(in_vertex2, in_verts%no_of_neigbors(in_vertex2)) = in_edge_index
      in_verts%get_neighbor_index(in_vertex1, in_verts%no_of_neigbors(in_vertex1)) = in_vertex2
      in_verts%get_neighbor_index(in_vertex2, in_verts%no_of_neigbors(in_vertex2)) = in_vertex1
    ENDDO


    ! complete cell-vertex/edge connectivity
    !CALL message ('complete_grid_connectivity','complete cell-vertex/edge connectivity')
    in_verts%no_of_neigbors(:) = 0
    DO in_cell_index=1,no_of_input_cells
      ! fill the edge info
      DO j=1,max_in_cell_vertices
        in_edge_index = in_cells%get_edge_index(in_cell_index,j)
        IF (in_edge_index /= 0) THEN
          ! add this cell to the edge
          IF (in_edges%get_cell_index(in_edge_index,1) == 0) THEN
            in_edges%get_cell_index(in_edge_index,1) = in_cell_index
          ELSE
            in_edges%get_cell_index(in_edge_index,2) = in_cell_index
          ENDIF
        ENDIF

        in_vertex_index = in_cells%get_vertex_index(in_cell_index,j)
        IF (in_vertex_index /= 0) THEN
          ! add this cell to the vertex
          in_verts%no_of_neigbors(in_vertex_index) = &
            & in_verts%no_of_neigbors(in_vertex_index) + 1
          in_verts%get_cell_index(in_vertex_index, in_verts%no_of_neigbors(in_vertex_index)) = &
            & in_cell_index
        ENDIF

      ENDDO
    ENDDO

    ! finally fill the cell-cell neigboring connectivity
    !CALL message ('complete_grid_connectivity','fill the cell-cell neigboring connectivity')
    !WRITE(0,*) "--------- verts in cells -----------------"
    DO in_cell_index=1,no_of_input_cells
      DO j=1,max_in_cell_vertices
        in_edge_index = in_cells%get_edge_index(in_cell_index,j)
        IF (in_edges%get_cell_index(in_edge_index,1) == in_cell_index) THEN
          in_cells%get_neighbor_index(in_cell_index,j) = &
            & in_edges%get_cell_index(in_edge_index,2)
        ELSE
          in_cells%get_neighbor_index(in_cell_index,j) = &
            & in_edges%get_cell_index(in_edge_index,1)
        ENDIF
      ENDDO
    !  WRITE(0,*) in_cells%get_vertex_index(in_cell_index, :)
    ENDDO
    !WRITE(0,*) "--------- cells in verts -----------------"
    !DO  in_vertex_index=1,no_of_input_verts
    !  WRITE(0,*) in_verts%get_cell_index(in_vertex_index, :)
    !ENDDO

    !CALL message ('complete_grid_connectivity','ends')
  END SUBROUTINE complete_grid_connectivity
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !   INTEGER FUNCTION get_vertex_index_in_cell(cells, cell_index, vertex_index)
  !>
  !! Returns the index of vertex_index of the cells%vertex_index. Private
  INTEGER FUNCTION get_vertex_index_in_cell(cells, cell_index, vertex_index)
    !  TYPE(grid_cells), POINTER, INTENT(IN)  :: cells
    TYPE(t_grid_cells), POINTER :: cells
    INTEGER, INTENT(in) :: cell_index, vertex_index

    DO get_vertex_index_in_cell = 1, max_cell_vertices
      IF (cells%get_vertex_index(cell_index,get_vertex_index_in_cell) == vertex_index) return;
    ENDDO

    get_vertex_index_in_cell = -1
    WRITE(message_text,'(a,i9)') 'maxCellVertices=', max_cell_vertices
    CALL message ('', TRIM(message_text))
    WRITE(message_text,'(a,i9)') 'cellIndex=', cell_index
    CALL message ('', TRIM(message_text))
    WRITE(message_text,'(a,i9)') 'vertexIndex=', vertex_index
    CALL message ('', TRIM(message_text))

    CALL finish ('get_vertex_index_in_cell', 'Vertex was not found!')

  END FUNCTION get_vertex_index_in_cell
  !-------------------------------------------------------------------------


END MODULE mo_local_grid_refinement

