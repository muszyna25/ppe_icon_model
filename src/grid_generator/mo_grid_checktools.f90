!-------------------------------------------------------------------------------------
! mo_grid_checktools
!>
!! A collection of basic methods for checking the connectivity
!! and geometry of the the ICON grid
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
MODULE mo_grid_checktools
#include "grid_definitions.inc"
  !-------------------------------------------------------------------------
  USE mo_kind,           ONLY: wp
  USE mo_io_units,       ONLY: find_next_free_unit, filename_max
  USE mo_exception,      ONLY: message, finish
  USE mo_local_grid
  USE mo_io_local_grid
  USE mo_base_geometry,  ONLY: t_cartesian_coordinates, t_geographical_coordinates, &
    & gc2cc, arc_length, norma, sin_cc, gvec2cvec, triangle_area!, sphere_tanget_coordinates
  USE mo_math_constants, ONLY: rad2deg, deg2rad
  USE mo_local_grid_geometry,  ONLY: edges_cell_angle, edges_normal_angle, no_angle, &
    & get_cell_barycenters
    !, get_triangle_circumcenters, geographical_to_cartesian
  USE mo_grid_toolbox , ONLY :  inverse_connectivity_verts! get_basic_dual_grid
  USE mo_statistics_tools
  USE mo_physical_constants, ONLY: earth_radius
  
  IMPLICIT NONE

  !-------------------------------------------------------------------------
  PRIVATE

  ! !VERSION CONTROL:
  CHARACTER(LEN=*), PARAMETER :: version = '$Id$'

  PUBLIC :: check_read_write_grid
  PUBLIC :: check_grid_connectivity_edges
  PUBLIC :: check_grid_connectivity_cells
  PUBLIC :: check_grid_connectivity
  PUBLIC :: check_grid_geometry
  PUBLIC :: check_grid_geometry_verts
  PUBLIC :: check_grid
  PUBLIC :: check_parent_child_grid
  PUBLIC :: check_grid_file, grid_statistics_file
  PUBLIC :: check_inverse_connect_verts
  PUBLIC :: calculate_triangle_properties
  
  !----------------------------------------
  REAL(wp), PARAMETER  :: EDGE_CENTER_LIMIT = 1e-3_wp
  REAL(wp), PARAMETER  :: ANGLE_INTRIANGLE_LIMIT = 8.0e-0_wp
  REAL(wp), PARAMETER  :: PARENT_CHILD_ANGLE = 20.0e-0_wp
  INTEGER :: latex_file_v = 500
  INTEGER :: latex_file_e = 501
  INTEGER :: latex_file_c = 502

CONTAINS

  !-------------------------------------------------------------------------
  SUBROUTINE calculate_triangle_properties(lon1, lat1, lon2, lat2, lon3, lat3 )
    REAL(wp), INTENT(in) :: lon1, lat1, lon2, lat2, lon3, lat3

    TYPE(t_cartesian_coordinates) :: x0, x1, x2
    TYPE(t_geographical_coordinates) :: v
    REAL(wp) :: area


    v%lon = lon1 * deg2rad
    v%lat = lat1 * deg2rad    
    x0 = gc2cc(v)
    
    v%lon = lon2 * deg2rad
    v%lat = lat2 * deg2rad    
    x1 = gc2cc(v)

    v%lon = lon3 * deg2rad
    v%lat = lat3 * deg2rad    
    x2 = gc2cc(v)

    area = triangle_area(x0, x1, x2)
    write(0,*) " Normed triangle area=", area
    write(0,*) " Earth triangle area (km) =", area * earth_radius * earth_radius * 1e-6_wp
    
  END SUBROUTINE calculate_triangle_properties
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  SUBROUTINE check_grid_file(in_file)
    CHARACTER(LEN=filename_max), INTENT(in) :: in_file

    INTEGER :: grid_id

    CALL message('===============','===============')
    CALL message('check_grid_file',TRIM(in_file))
    grid_id = read_new_netcdf_grid(in_file)
    CALL check_grid_connectivity(grid_id)
    CALL delete_grid(grid_id)

  END SUBROUTINE check_grid_file
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE grid_statistics_file(in_file, in_latex_file_name)
    CHARACTER(LEN=filename_max), INTENT(in) :: in_file
    CHARACTER(LEN=filename_max), INTENT(in) :: in_latex_file_name

    CHARACTER(LEN=filename_max) :: latex_file_name
    INTEGER :: grid_id, error_status

    CALL message('===============','===============')
    CALL message('check_grid_file',TRIM(in_file))
    grid_id = read_new_netcdf_grid(in_file)

    ! latex files
    ! for vertices
    WRITE(latex_file_name,'(a,a)') TRIM(in_latex_file_name), "_verts.tex"
    latex_file_v = find_next_free_unit(100,1000)
    OPEN (latex_file_v, FILE=TRIM(latex_file_name),IOSTAT=error_status, &
      & POSITION='APPEND')
    IF (error_status /= 0) &
      & CALL finish("unable to open",TRIM(latex_file_name))
    WRITE(latex_file_v, '("\verb+",a,"+ ")', advance='no') TRIM(in_file)
    ! for edges
    WRITE(latex_file_name,'(a,a)') TRIM(in_latex_file_name), "_edges.tex"
    latex_file_e = find_next_free_unit(100,1000)
    OPEN (latex_file_e, FILE=TRIM(latex_file_name),IOSTAT=error_status, &
      & POSITION='APPEND')
    IF (error_status /= 0) &
      & CALL finish("unable to open",TRIM(latex_file_name))
    WRITE(latex_file_e, '("\verb+",a,"+ ")', advance='no') TRIM(in_file)
    ! for cells
    WRITE(latex_file_name,'(a,a)') TRIM(in_latex_file_name), "_cells.tex"
    latex_file_c = find_next_free_unit(100,1000)
    OPEN (latex_file_c, FILE=TRIM(latex_file_name),IOSTAT=error_status, &
      & POSITION='APPEND')
    IF (error_status /= 0) &
      & CALL finish("unable to open",TRIM(latex_file_name))
    WRITE(latex_file_c, '("\verb+",a,"+ ")', advance='no') TRIM(in_file)
    
    CALL check_grid_geometry(grid_id)
    
    CALL delete_grid(grid_id)
    
    WRITE(latex_file_v, '("\\ ")')
    CLOSE(latex_file_v)
    WRITE(latex_file_e, '("\\ ")')
    CLOSE(latex_file_e)
    WRITE(latex_file_c, '("\\ ")')
    CLOSE(latex_file_c)

  END SUBROUTINE grid_statistics_file
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE check_read_write_grid(in_file, out_file)
    CHARACTER(LEN=filename_max), INTENT(in) :: in_file, out_file

    INTEGER :: grid_id

    grid_id = read_new_netcdf_grid(in_file)
    CALL write_netcdf_grid(grid_id, out_file)
    CALL delete_grid(grid_id)

  END SUBROUTINE check_read_write_grid
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  SUBROUTINE check_inverse_connect_verts(in_file, out_file)
    CHARACTER(LEN=filename_max), INTENT(in) :: in_file, out_file

    INTEGER :: grid_id

    grid_id = read_new_netcdf_grid(in_file)
    CALL inverse_connectivity_verts(grid_id)
    CALL write_netcdf_grid(grid_id, out_file)
    CALL delete_grid(grid_id)

  END SUBROUTINE check_inverse_connect_verts
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE check_grid(in_grid_id)
    INTEGER, INTENT(in) :: in_grid_id

    CALL check_grid_connectivity(in_grid_id)
    CALL check_grid_geometry(in_grid_id)

  END SUBROUTINE check_grid
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE check_parent_child_grid(child_grid_id)
    INTEGER, INTENT(in) :: child_grid_id

    INTEGER :: parent_grid_id
    TYPE(t_grid)         , POINTER :: child_grid , parent_grid
!    TYPE(t_grid_cells)   , POINTER :: child_cells, parent_cells
    TYPE(t_grid_edges)   , POINTER :: child_edges, parent_edges
!    TYPE(t_grid_vertices), POINTER :: child_verts, parent_verts

    REAL(wp)  :: angle
    INTEGER child_edge,parent_edge


    CALL message('check_parent_child_grid...','')
    child_grid  => get_grid(child_grid_id)
!    child_cells => child_grid%cells
    child_edges => child_grid%edges
!    child_verts => child_grid%verts

    parent_grid_id  = child_grid%parent_grid_id

    IF (parent_grid_id == 0) THEN
      CALL message('','parent_grid_id == 0, nothing to do')
      RETURN
    ENDIF

    parent_grid     => get_grid(parent_grid_id)
!    parent_cells    => parent_grid%cells
    parent_edges    => parent_grid%edges
!    parent_verts    => parent_grid%verts

    ! check edges
    DO child_edge=1,child_edges%no_of_existedges
      parent_edge = child_edges%parent_index(child_edge)

      !-----------------------------------------------------------
      ! check if the angle between parent,child is ~ 90
      angle = edges_normal_angle(parent_grid_id, parent_edge, child_grid_id, child_edge)
      IF (angle > 100.0_wp) angle = angle - 180.0_wp

!       IF (ABS(angle) > PARENT_CHILD_ANGLE) &
! !       .OR. child_edge == 1001 .OR. child_edge == 5689 &
! !       .OR. child_edge == 5852 .OR. child_edge == 6020 &
!        & ) THEN
!          print *, parent_grid_id, "parent edge:", parent_edge
!          print *, "  center:",parent_edges%center(parent_edge)
!          print *, "  normal:",parent_edges%primal_normal(parent_edge)
!          print *, child_grid_id, "child edge:",child_edge
!          print *, "  center:",child_edges%center(child_edge)
!          print *, "  normal:",child_edges%primal_normal(child_edge)
!          print *, "angle:", angle
!       ENDIF
      IF (ABS(angle) > PARENT_CHILD_ANGLE ) THEN
         print *, parent_grid_id, "parent edge:", parent_edge
         print *, "  center:",parent_edges%center(parent_edge)
         print *, "  normal:",parent_edges%primal_normal(parent_edge)
         print *, child_grid_id, "child edge:",child_edge
         print *, "  center:",child_edges%center(child_edge)
         print *, "  normal:",child_edges%primal_normal(child_edge)
         print *, "angle:", angle
         CALL finish('check_parent_child_grid','parent,child edge not oriented')
      ENDIF
      !-----------------------------------------------------------

    ENDDO

  END SUBROUTINE check_parent_child_grid
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE check_grid_geometry(in_grid_id)
    INTEGER, INTENT(in) :: in_grid_id

     CALL check_grid_geometry_verts(in_grid_id)
     CALL check_grid_geometry_edges(in_grid_id)
     CALL check_grid_geometry_cells(in_grid_id)

  END SUBROUTINE check_grid_geometry
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE check_grid_geometry_verts(in_grid_id)
    INTEGER, INTENT(in) :: in_grid_id

    TYPE(t_grid), POINTER :: in_grid
!    TYPE(t_grid_cells), POINTER :: cells
    TYPE(t_grid_edges), POINTER :: edges
    TYPE(t_grid_vertices), POINTER :: verts

    TYPE(t_cartesian_coordinates) :: x,d
    INTEGER :: max_connectivity, vert_no, edge, i
    REAL(wp)  :: norm
    REAL(wp):: max_dual_area, min_dual_area
    REAL(wp):: max_dual_edge, min_dual_edge, max_dual_edge_ratio
    REAL(wp):: max_edge, min_edge, max_edge_ratio

!     INTEGER :: dual_grid_id
!     TYPE(t_grid), POINTER :: dual_grid
!     TYPE(t_cartesian_coordinates), POINTER :: dual_barycenters(:)
    REAL(wp):: max_barycenter_distance

    INTEGER :: dual_area_stats, dual_edge_ratio_stats, edge_ratio_stats

    in_grid => get_grid(in_grid_id)
!    cells => in_grid%cells
    edges => in_grid%edges
    verts => in_grid%verts
    max_connectivity = verts%max_connectivity

!     dual_grid_id = get_basic_dual_grid(in_grid_id)
!     CALL get_cell_barycenters(dual_grid_id)
!     dual_grid => get_grid(dual_grid_id)
!     dual_barycenters => dual_grid%cells%cartesian_center

    WRITE(0,*) "-----------------------------"
    WRITE(0,*) " Checking geometry of ", verts%no_of_existvertices, " verts."

    dual_area_stats       = new_statistic()
    dual_edge_ratio_stats = new_statistic()
    edge_ratio_stats      = new_statistic()
    
    max_barycenter_distance = 0.0_wp
    DO vert_no=1, verts%no_of_existvertices
      x = gc2cc(verts%vertex(vert_no))
      d%x = verts%cartesian(vert_no)%x - x%x
      norm = norma(d)
      IF (norm > 1.0E-8_wp) THEN
        CALL print_grid_vertex(in_grid_id, vert_no)
        WRITE(0,*) 'cartesian - lonlat =',d
        WRITE(0,*) 'norma =',norm
        CALL finish('check_grid_geometry_verts',&
          & 'Not acceptable cartesian - lonlat')
      ENDIF

      CALL add_statistic_to(dual_area_stats, verts%dual_area(vert_no))
      !       max_barycenter_distance = MAX(max_barycenter_distance, &
!        & arc_length(verts%cartesian(vert_no), dual_barycenters(vert_no)))

      max_dual_edge = 0.0_wp
      min_dual_edge = 1E100_wp
      max_edge = 0.0_wp
      min_edge = 1E100_wp
      DO i=1,max_connectivity
        edge = verts%get_edge_index(vert_no, i)
        IF (edge > 0) THEN
          max_dual_edge = MAX(max_dual_edge, edges%dual_edge_length(edge))
          min_dual_edge = MIN(min_dual_edge, edges%dual_edge_length(edge))
          max_edge = MAX(max_edge, edges%primal_edge_length(edge))
          min_edge = MIN(min_edge, edges%primal_edge_length(edge))
        ENDIF
      ENDDO !i=1,max_connectivity
      
      CALL add_statistic_to(dual_edge_ratio_stats, max_dual_edge/min_dual_edge)
      CALL add_statistic_to(edge_ratio_stats, max_edge/min_edge)
      
    ENDDO

    min_dual_area       = min_statistic_of(dual_area_stats)
    max_dual_area       = max_statistic_of(dual_area_stats)
    max_dual_edge_ratio = max_statistic_of(dual_edge_ratio_stats)
    max_edge_ratio      = max_statistic_of(edge_ratio_stats)

    CALL delete_statistic(dual_area_stats)
    CALL delete_statistic(dual_edge_ratio_stats)
    
    WRITE(0,*) 'Min/Max dual area:', min_dual_area, max_dual_area, max_dual_area/min_dual_area
    WRITE(0,*) 'Max dual cell dual edge ratio:', max_dual_edge_ratio
    WRITE(0,*) 'Max dual cell edge ratio:', max_edge_ratio
    WRITE(0,*) 'Max vertex-dual baryceneter distance:', max_barycenter_distance * &
      & in_grid%sphere_radius

    WRITE(latex_file_v, '(" & ", f8.2, " & ", f6.3," & ",f6.3," & ",f6.3)', advance='no')  &
      & max_barycenter_distance *  in_grid%sphere_radius &
      & / 1000.0_wp, max_dual_area/min_dual_area, &
      & max_dual_edge_ratio, max_edge_ratio

!     CALL delete_grid(dual_grid_id)

  END SUBROUTINE check_grid_geometry_verts
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE check_grid_geometry_edges(in_grid_id)
    INTEGER, INTENT(in) :: in_grid_id

    TYPE(t_grid), POINTER :: in_grid
    TYPE(t_grid_cells), POINTER :: cells
    TYPE(t_grid_edges), POINTER :: edges
    TYPE(t_grid_vertices), POINTER :: verts
    INTEGER :: edge_no
    INTEGER :: vert1,vert2, cell1, cell2
    REAL(wp):: min_prime_edge, min_dual_edge, max_prime_edge, max_dual_edge
    REAL(wp):: area_ratio, max_min_cell_area, max_min_dual_area, max_edge_cell_ratio
    TYPE(t_cartesian_coordinates) :: x1,x2,x3,v1,v2,v3,v4,v5,vn
    TYPE(t_cartesian_coordinates) :: x_unit,y_unit,z_unit,xy
    REAL(wp)  :: dproduct,sin1,sin2

    INTEGER :: edge_cell_offcenter_stat, edge_vert_offcenter_stat, dual_edge_stat
    REAL(wp)::  max_edge_vert_offcenter, max_edge_cell_offcenter

      
    CHARACTER(*), PARAMETER :: method_name = "check_grid_geometry_edges"
   
    edge_cell_offcenter_stat   = new_statistic(mode=ADD_MAX_RATIO)
    edge_vert_offcenter_stat   = new_statistic(mode=ADD_MAX_RATIO)
    dual_edge_stat             = new_statistic()

    x_unit%x = (/ 1._wp,0._wp,0._wp /)
    y_unit%x = (/ 0._wp,1._wp,0._wp /)
    z_unit%x = (/ 0._wp,0._wp,1._wp /)
    xy%x = (/ 1._wp,1._wp,0._wp /)
!     write(0,*), '0=',sin_cc(x_unit,x_unit)
!     write(0,*), '0=',sin_cc(y_unit,y_unit)
!     write(0,*), '0=',sin_cc(z_unit,z_unit)
!     write(0,*), '1=',sin_cc(x_unit,y_unit)
!     write(0,*), '1=',sin_cc(x_unit,z_unit)
!     write(0,*), '1=',sin_cc(y_unit,z_unit)
!     write(0,*), SQRT(2.0_wp)*0.5_wp,"=",sin_cc(xy,x_unit)
!     write(0,*), SQRT(2.0_wp)*0.5_wp,"=",sin_cc(xy,y_unit)
!     write(0,*), "1=",sin_cc(xy,z_unit)


    in_grid => get_grid(in_grid_id)
    cells => in_grid%cells
    edges => in_grid%edges
    verts => in_grid%verts

    WRITE(0,*) "-----------------------------"
    WRITE(0,*) " Checking geometry of ", edges%no_of_existedges, " edges."
    min_prime_edge = 1E100_wp
    min_dual_edge  = 1E100_wp
    max_prime_edge = 0.0_wp
    max_dual_edge  = 0.0_wp
    max_min_cell_area = 0.0_wp
    max_min_dual_area = 0.0_wp
    max_edge_cell_ratio = 0.0_wp
    DO edge_no=1,edges%no_of_existedges
      min_prime_edge = MIN(min_prime_edge, edges%primal_edge_length(edge_no))
      max_prime_edge = MAX(max_prime_edge, edges%primal_edge_length(edge_no))

      CALL add_statistic_to(dual_edge_stat, edges%dual_edge_length(edge_no))

      CALL add_statistic_to(edge_vert_offcenter_stat, &
        & edges%get_edge_vert_length(edge_no,1), edges%get_edge_vert_length(edge_no,2))

      vert1 = edges%get_vertex_index(edge_no,1)
      vert2 = edges%get_vertex_index(edge_no,2)
      x1 = verts%cartesian(vert1)
      x2 = verts%cartesian(vert2)
      x3 = gc2cc(edges%center(edge_no))
      ! project x3 on the line x1x2
      v1%x = x2%x - x1%x
      sin1 = sin_cc(x3,v1)
      sin2 = sin_cc(x1,v1)
      v2%x = x3%x * (sin2/sin1)

      ! check if center is on x2-x1
      v3%x = v2%x - x1%x
!      IF (MINVAL(ABS(v3%x)) > 1e-8) THEN
      IF (MINVAL(ABS(v3%x)) < 0.0_wp) THEN
        v4%x = v1%x / v3%x
        v5%x = ABS(v4%x - (SUM(v4%x) / 3.0_wp))
        IF (MAXVAL(v5%x) > EDGE_CENTER_LIMIT) THEN
          write(0,*) 'edges verts:'
          write(0,*) x1
          write(0,*) x2
          write(0,*) 'edges center:'
          write(0,*) x3
          write(0,*) 'projection'
          write(0,*) v2
          write(0,*) 'x1x2:'
          write(0,*) v1
          write(0,*) 'x1x3:'
          write(0,*) v3
          write(0,*) 'ratio:'
          write(0,*) v4
          write(0,*) 'diff'
          write(0,*) v5
          CALL finish(method_name,'edges center not on edges arc')
        ENDIF
      ENDIF

      ! check if normal is perpenticular to the edges
      ! get local coordinates
     CALL gvec2cvec(edges%primal_normal(edge_no)%v1,&
                    edges%primal_normal(edge_no)%v2, &
                     edges%center(edge_no)%lon, &
                     edges%center(edge_no)%lat, &
                     vn%x(1),vn%x(2),vn%x(3))
       !see if they are orthogonal
       dproduct = DOT_PRODUCT(v1%x,vn%x)

!       CALL sphere_tanget_coordinates(edges%center(edge_no),x,y)
!       ! get x1x2 on local coordinates
!       v3%x(1) = DOT_PRODUCT(v1%x,x%x)
!       v3%x(2) = DOT_PRODUCT(v1%x,y%x)
!       v3%x(3) = 0.0_wp
!       ! get primal normal
!       v4%x(1) = edges%primal_normal(edge_no)%v1
!       v4%x(2) = edges%primal_normal(edge_no)%v2
!       v4%x(3) = 0.0_wp
!       !see if they are orthogonal
!       dproduct = DOT_PRODUCT(v3%x,v4%x)
      IF (ABS(dproduct) > 5e-3_wp) THEN
        write(0,*) "edges :"
        write(0,*) v1
        write(0,*) "normal :"
        write(0,*) vn
        write(0,*) "dproduct:", dproduct
        CALL finish(method_name,'normal not orthogonal to edges')
      ENDIF

      cell1 = edges%get_cell_index(edge_no, 1)
      cell2 = edges%get_cell_index(edge_no, 2)
      IF (cell1 < 1 .AND. cell2 < 1) THEN
        CALL finish(method_name,'cell1 < 1 .AND. cell2 < 1')
      ENDIF
      IF (cell1 > 0 .AND. cell2 > 0) THEN
        
        CALL add_statistic_to(edge_cell_offcenter_stat, &
          & edges%get_edge_cell_length(edge_no,1), edges%get_edge_cell_length(edge_no,2))
      
        IF (cells%area(cell1) > cells%area(cell2)) THEN
          area_ratio = cells%area(cell1) / cells%area(cell2)
        ELSE
          area_ratio = cells%area(cell2) / cells%area(cell1)
        ENDIF
        max_min_cell_area = MAX(max_min_cell_area, area_ratio)
      ENDIF

      IF (verts%dual_area(vert1) > verts%dual_area(vert2)) THEN
        area_ratio = verts%dual_area(vert1) / verts%dual_area(vert2)
      ELSE
        area_ratio = verts%dual_area(vert2) / verts%dual_area(vert1)
      ENDIF
      max_min_dual_area = MAX(max_min_dual_area, area_ratio)
      ! check angles between edges
!       DO k=1,2
!         cell_idx = edges%get_cell_index(edge_no,1)
!         IF (cell_idx == 0) CYCLE
!         DO j=1,3
!           cell_edge = cells%get_edge_index(cell_idx,j)
!           IF (cell_edge == edge_no) CYCLE
!           angle = edges_normal_angle(in_grid_id, edge_no, in_grid_id, cell_edge)
!           IF (angle > 100.0_wp) angle = angle - 60.0_wp
!           IF (ABS(angle -60.0_wp) > ANGLE_INTRIANGLE_LIMIT) THEN
!             CALL print_grid_cell(in_grid_id, cell_idx)
!             CALL print_grid_edge(in_grid_id, edge_no)
!             CALL print_grid_edge(in_grid_id, cell_edge)
!             print *, 'angle = ', angle
!             CALL finish('check_grid_geometry_edges',&
!              'angle between cell edges are not ~ 60')
!           ENDIF

          ! check the quads
!           q_cell = edges% get _cell_index(cell_edge,1)
!           IF (q_cell == cell_idx) q_cell = edges% get _cell_index(cell_edge,2)
!           IF (q_cell == 0) CYCLE
!           DO q=1,3
!             q_edge = cells% get _edge_index(q_cell,q)
!             IF (q_edge == cell_edge) CYCLE
!
!             IF ( edges_common_vertex(in_grid_id,edge_no,q_edge) < 0) THEN
!               angle = edges_normal_angle(in_grid_id, edge_no, in_grid_id, cell_edge)
!               IF (angle > 110.0_wp) angle = angle - 90.0_wp
!               IF (ABS(angle) > 35.0E0_wp) THEN
!                 CALL print_grid_cell(in_grid_id,q_cell)
!                 CALL print_grid_edge(in_grid_id, edge_no)
!                 CALL print_grid_edge(in_grid_id, q_edge)
!                 print *, 'angle in quad = ', angle
!                 CALL finish('check_grid_geometry_edges',&
!                 'angle between quad edges is not ~ 0')
!               ENDIF
!             ENDIF
!
!           ENDDO


!        ENDDO
!      ENDDO

    ENDDO ! edge_no=1,edges%no_of_existedges

    min_dual_edge  = min_statistic_of(dual_edge_stat)
    max_dual_edge  = max_statistic_of(dual_edge_stat)
    
    WRITE(0,*) "Adjacent dual max area ratio:", max_min_dual_area
    WRITE(0,*) "Adjacent cell max area ratio:", max_min_cell_area
    WRITE(0,*) "global min max prime edges length:", min_prime_edge, max_prime_edge, &
      & max_prime_edge/min_prime_edge
    WRITE(0,*) "global min max dual edges length:",  min_dual_edge, max_dual_edge, &
      & max_dual_edge/min_dual_edge
    WRITE(0,*) "mean dual edges length:",  mean_statistic_of(dual_edge_stat)

    max_edge_vert_offcenter = max_statistic_of(edge_vert_offcenter_stat)
    max_edge_cell_offcenter = max_statistic_of(edge_cell_offcenter_stat)
    
    WRITE(0,*) "max_edge_vert_offcenter:",  max_edge_vert_offcenter
    WRITE(0,*) "max_edge_cell_offcenter:",  max_edge_cell_offcenter
    
    CALL delete_statistic(edge_vert_offcenter_stat)
    CALL delete_statistic(edge_cell_offcenter_stat)
    CALL delete_statistic(dual_edge_stat)
   
    WRITE(latex_file_e, '( " & ", f6.3, " & ", f6.3," & ",f6.3, " & ",f6.3, " & ",f6.3)',&
      advance='no')  &
      & max_min_cell_area, max_min_dual_area, &
      & max_prime_edge/min_prime_edge, max_dual_edge/min_dual_edge, &
      & max_edge_cell_offcenter

  END SUBROUTINE check_grid_geometry_edges
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE check_grid_geometry_cells(in_grid_id)
    INTEGER, INTENT(in) :: in_grid_id

    TYPE(t_grid), POINTER :: in_grid
    TYPE(t_grid_cells), POINTER :: cells
    TYPE(t_grid_edges), POINTER :: edges
!    TYPE(t_grid_vertices), POINTER :: verts
    INTEGER :: no_of_cells
    INTEGER :: cell_no,edge,k,max_no_of_vertices

    REAL(wp):: angle,min_cell_angle, max_cell_angle, ratio
    REAL(wp):: max_edge_ratio, max_prime_edge_cell, min_prime_edge_cell
    REAL(wp):: max_dualedge_ratio, max_dual_edge_cell, min_dual_edge_cell
    REAL(wp):: max_prime_edge, min_prime_edge, max_dual_edge, min_dual_edge
    INTEGER :: min_cell_angle_idx, max_cell_angle_idx

    TYPE(t_cartesian_coordinates), POINTER :: barycenters(:)
    REAL(wp):: max_centers_diff, centers_diff
    REAL(wp):: max_area, min_area
    REAL(wp):: edge_cell_length, max_cell_edge_cell_length, min_cell_edge_cell_length
    REAL(wp):: max_edge_cell_length_ratio, max_centers_diff_ratio

    in_grid => get_grid(in_grid_id)
    cells => in_grid%cells
    edges => in_grid%edges
!    verts => in_grid%verts
    no_of_cells = cells%no_of_existcells
    max_no_of_vertices =cells%max_no_of_vertices

    WRITE(0,*) "-----------------------------"
    WRITE(0,*) " Checking geometry of ", no_of_cells, " cells."

    min_cell_angle = 10.0_wp
    max_cell_angle = -10.0_wp
    max_edge_ratio = 0.0_wp
    max_dualedge_ratio = 0.0_wp
    max_centers_diff = 0.0_wp
    max_area = 0.0_wp
    min_area = 1E100_wp
    max_edge_cell_length_ratio = 0.0_wp
    max_centers_diff_ratio = 0.0_wp

    ! get barycenters for testing
    NULLIFY(barycenters)
    CALL get_cell_barycenters(in_grid_id, output_centers=barycenters)
    ! CALL geographical_to_cartesian(cells%center, no_of_cells,cells%cartesian_center)
    ! CALL get_triangle_circumcenters(in_grid_id)

    DO cell_no=1,no_of_cells
      max_prime_edge = 0.0_wp
      min_prime_edge = 1E100_wp
      max_dual_edge = 0.0_wp
      min_dual_edge = 1E100_wp
      max_cell_edge_cell_length =  0.0_wp
      min_cell_edge_cell_length =  1E100_wp
      DO k=1,max_no_of_vertices

        ! check angles
        angle = edges_cell_angle(in_grid_id, cell_no, k)
        IF (angle == no_angle) CYCLE
        IF (angle < min_cell_angle) THEN
          min_cell_angle = angle
          min_cell_angle_idx = cell_no
        ENDIF
        IF (angle > max_cell_angle) THEN
          max_cell_angle = angle
          max_cell_angle_idx = cell_no
        ENDIF

        ! check edges
        edge = cells%get_edge_index(cell_no, k)
        IF (edge /= 0) THEN
          min_prime_edge = MIN(min_prime_edge, edges%primal_edge_length(edge))
          max_prime_edge = MAX(max_prime_edge, edges%primal_edge_length(edge))
          min_dual_edge =  MIN(min_dual_edge,  edges%dual_edge_length(edge))
          max_dual_edge =  MAX(max_dual_edge,  edges%dual_edge_length(edge))
          IF (edges%get_cell_index(edge,1) == cell_no) THEN
            edge_cell_length = edges%get_edge_cell_length(edge,1)
          ELSE
            edge_cell_length = edges%get_edge_cell_length(edge,2)
          ENDIF
          max_cell_edge_cell_length = MAX(max_cell_edge_cell_length,edge_cell_length)
          min_cell_edge_cell_length = MIN(min_cell_edge_cell_length,edge_cell_length)
          
     !     print *, "max/min dual", max_dual_edge, min_dual_edge
        ENDIF

      ENDDO  ! k=1,max_no_of_vertices

      max_edge_cell_length_ratio = MAX(max_edge_cell_length_ratio, &
        & max_cell_edge_cell_length/min_cell_edge_cell_length)
      
      ratio = max_prime_edge / min_prime_edge
      IF (ratio > max_edge_ratio) THEN
        max_edge_ratio = ratio
        min_prime_edge_cell =  min_prime_edge
        max_prime_edge_cell =  max_prime_edge
      ENDIF
      ratio = max_dual_edge / min_dual_edge
      IF (ratio > max_dualedge_ratio) THEN
        max_dualedge_ratio = ratio
        min_dual_edge_cell =  min_dual_edge
        max_dual_edge_cell =  max_dual_edge
      ENDIF

      centers_diff = arc_length(cells%cartesian_center(cell_no), barycenters(cell_no)) &
        & *  in_grid%sphere_radius
      max_centers_diff = MAX(max_centers_diff, centers_diff)
      max_centers_diff_ratio = MAX(max_centers_diff_ratio,&
        & centers_diff/min_cell_edge_cell_length)
!       write(0,*) "centers_diff:", centers_diff, min_cell_edge_cell_length, &
!         max_centers_diff_ratio

      max_area = MAX(max_area, cells%area(cell_no))
      min_area = MIN(min_area, cells%area(cell_no))

      !print *, "center:", cells%cartesian_center(cell_no)
      !print *, "barycenter:", barycenters(cell_no)
      !print *, "max_centers_diff=",max_centers_diff

    ENDDO ! cell_no=1,no_of_cells

    write(0,*) 'Max barycenter-center distance:', max_centers_diff *  in_grid%sphere_radius
    write(0,*) 'Min/Max cell area:', min_area, max_area, max_area/min_area
    write(0,*) 'Min/Max cell angle:', min_cell_angle * rad2deg, max_cell_angle * rad2deg
    write(0,*) 'cell max_edge_ratio:', min_prime_edge_cell, max_prime_edge_cell, max_edge_ratio
    write(0,*) 'cell max_dualedge_ratio:', &
      & min_dual_edge_cell, max_dual_edge_cell, max_dualedge_ratio

    DEALLOCATE(barycenters)

    WRITE(latex_file_c, '( " & ", f6.3, " & ", f6.3," & ",&
      & f6.3 ," & ",f6.3, " & ",f8.2," & ",f6.3)', &
      & advance='no')  &
      & max_area/min_area, max_edge_ratio, max_dualedge_ratio, max_edge_cell_length_ratio, &
      & max_centers_diff /1000.0_wp, max_centers_diff_ratio
   
   ! CALL print_cell_angles(in_grid_id, min_cell_angle_idx, 'Cell with min edge angle')
   ! CALL print_cell_angles(in_grid_id, max_cell_angle_idx, 'Cell with max edge angle')

  END SUBROUTINE check_grid_geometry_cells
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE check_grid_connectivity(in_grid_id)
    INTEGER, INTENT(in) :: in_grid_id

    CALL check_grid_connectivity_verts(in_grid_id)
    CALL check_grid_connectivity_edges(in_grid_id)
    CALL check_grid_connectivity_cells(in_grid_id)

  END SUBROUTINE check_grid_connectivity
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE check_grid_connectivity_cells(in_grid_id)
    INTEGER, INTENT(in) :: in_grid_id

    TYPE(t_grid), POINTER :: in_grid
    TYPE(t_grid_cells), POINTER :: cells
    TYPE(t_grid_edges), POINTER :: edges
    TYPE(t_grid_vertices), POINTER :: verts
    INTEGER :: no_of_cells
    INTEGER :: cell,ncell,j,k,vertex
    LOGICAL :: found

    in_grid => get_grid(in_grid_id)
    cells => in_grid%cells
    edges => in_grid%edges
    verts => in_grid%verts
    no_of_cells = cells%no_of_existcells

    CALL message("",'Checking connectivity for cells...')
    !print*, 'Level 1 indexes=',cell%start_idx(1,1),cell%end_idx(1,1)
    DO cell=1,no_of_cells
      DO j=1,cells%max_no_of_vertices
        !----------------------  
        !check neigbors
        found = .false.
        ncell = cells%get_neighbor_index(cell,j) 
        IF (ncell /= 0) THEN
          DO k=1,cells%max_no_of_vertices
            IF (cells%get_neighbor_index(ncell,k) == cell) THEN
              found=.true.
              EXIT
            ENDIF
          ENDDO
          IF (.not. found) THEN
            write(0,*) cell, ncell, ' Neigbors dont match.'
            CALL finish('check_grid_connectivity_cells', 'Wrong cell connectivity')
          ENDIF
        ENDIF
        !----------------------
        !check vertices
        found = .false.
        vertex = cells%get_vertex_index(cell,j)
        IF (vertex /= 0) THEN
          DO k=1,verts%max_connectivity
            IF (verts%get_cell_index(vertex,k) == cell) THEN
              found=.true.
              EXIT
            ENDIF
          ENDDO
          IF (.not. found) THEN
            write(0,*) cell, vertex, ' cell was not found in vertex.'
            CALL finish('check_grid_connectivity_cells', 'Wrong cell connectivity')
          ENDIF
        ENDIF

        
      ENDDO
    ENDDO
   

!     DO i=cell%start_idx(1,1),cell%end_idx(1,1)
!       neib_cells = 0
!       IF (cell%refin_ctrl(i) /= 1) THEN
!         print*, 'cell,rfn=',i,cell%refin_ctrl(i)
!         CALL finish('','refin_ctrl(i) /= 1')
!       ENDIF
! 
!       DO j=1,cell%max_no_of_vertices
!         edge_no = cell%get_edge_index(i,j)
! !         IF (edge_no == 0) THEN
! !           print*, 'cell=',i
! !           CALL finish('edge_no == 0','')
! !         ENDIF
!         IF (cell%get_neighbor_index(i,j) /= 0) &
!            neib_cells = neib_cells + 1
!       ENDDO
! 
! !       IF (neib_cells >= cell%max_no_of_vertices) THEN
! !         print*, 'cell, neighbors=',i, neib_cells
! !         CALL finish('neighbors >= max_no_of_vertice','')
! !       ENDIF
! 
!     ENDDO
! 
!     ! print*, 'Level 2 indexes=',cell%start_idx(2,1),cell%no_of_existcells
!     DO i=cell%start_idx(2,1),cell%no_of_existcells
!       neib_cells = 0
!       IF (cell%refin_ctrl(i) == 1) THEN
!         print*, 'cell,rfn=',i,cell%refin_ctrl(i)
!         CALL finish('','refin_ctrl(i) == 1')
!       ENDIF
! 
!       DO j=1,cell%max_no_of_vertices
!         edge_no = cell%get_edge_index(i,j)
! !         IF (edge_no == 0) THEN
! !           print*, 'cell=',i
! !           CALL finish('edge_no == 0','')
! !         ENDIF
!         IF (cell%get_neighbor_index(i,j) /= 0) &
!            neib_cells = neib_cells + 1
!       ENDDO
! 
! !       IF (neib_cells /= cell%max_no_of_vertices) THEN
! !         print*, 'cell, neighbors=',i, neib_cells
! !         CALL finish('neighbors >= max_no_of_vertice','')
! !       ENDIF
! 
!     ENDDO

  END SUBROUTINE check_grid_connectivity_cells
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE check_grid_connectivity_edges(in_grid_id)
    INTEGER, INTENT(in) :: in_grid_id

    TYPE(t_grid), POINTER :: in_grid
!    TYPE(t_grid_cells), POINTER :: cells
    TYPE(t_grid_edges), POINTER :: edges
    TYPE(t_grid_vertices), POINTER :: verts
    INTEGER :: no_of_edges
    INTEGER :: cell1,cell2
    INTEGER :: i,j,k,edge,vertex
    LOGICAL :: found

    CALL message("",'Checking connectivity for edges...')
    in_grid => get_grid(in_grid_id)
!    cells => in_grid%cells
    edges => in_grid%edges
    verts => in_grid%verts
    no_of_edges = edges%no_of_existedges

    DO edge=1,no_of_edges
      DO j=1,2
        found=.false.
        vertex=edges%get_vertex_index(edge,j)
        DO k=1,verts%max_connectivity
          IF (verts%get_edge_index(vertex,k) == edge) THEN
            found=.true.
            EXIT
          ENDIF
        ENDDO
        IF (.not. found) THEN
          write(0,*) edge, vertex, ' edge was not found in vertex.'
          CALL finish('check_grid_connectivity_edges', 'Wrong edge connectivity')
        ENDIF        
      ENDDO !j=1,2
    ENDDO ! edge=1,no_of_edges
    
    !print*, 'Level 1 indexes=',edge%start_idx(1,1),edge%end_idx(1,1)
    DO i=edges%start_idx(1,1),edges%end_idx(1,1)
      cell1 = edges%get_cell_index(i,1)
      cell2 = edges%get_cell_index(i,2)
      IF (cell1 /= 0 .and.cell2 /= 0) THEN

        CALL check_edge_incell(in_grid_id,i, cell1)
        CALL check_edge_incell(in_grid_id,i, cell2)
        ! CALL finish('','cell1 /= 0 .and.cell2 /= 0')
      ENDIF
      IF (cell1 == 0 .and.cell2 == 0) THEN
        print*, 'edge,cells=',i,cell1,cell2
        CALL finish('','cell1 == 0 .and.cell2 == 0')
      ENDIF
      IF (cell1 /= 0) CALL check_edge_incell(in_grid_id,i, cell1)
      IF (cell2 /= 0) CALL check_edge_incell(in_grid_id,i, cell2)
      IF (edges%refin_ctrl(i) /= 1) THEN
        print*, 'edge,cells,rfn=',i,cell1,cell2,edges%refin_ctrl(i)
        CALL finish('','refin_ctrl(i) /= 1')
      ENDIF

    ENDDO

    ! print*, 'Level 2 indexes=',edge%start_idx(2,1),edge%no_of_existedges
    DO i=edges%start_idx(2,1),edges%no_of_existedges
      cell1 = edges%get_cell_index(i,1)
      cell2 = edges%get_cell_index(i,2)
     ! print*, 'edge,cells=',i,cell1,cell2
      IF (cell1 == 0 .or.cell2 == 0) THEN
        print*, 'edge,cells=',i,cell1,cell2
        CALL message('','cell1 /= 0 .or.cell2 /= 0')
      ENDIF
      IF (cell1 /=0) CALL check_edge_incell(in_grid_id,i, cell1)
      IF (cell2 /=0) CALL check_edge_incell(in_grid_id,i, cell2)

    ENDDO

  END SUBROUTINE check_grid_connectivity_edges
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE check_grid_connectivity_verts(in_grid_id)
    INTEGER, INTENT(in) :: in_grid_id

    TYPE(t_grid), POINTER :: in_grid
!    TYPE(t_grid_cells), POINTER :: cells
    TYPE(t_grid_edges), POINTER :: edges
    TYPE(t_grid_vertices), POINTER :: verts
    TYPE(t_grid_cells), POINTER :: cells

    INTEGER :: max_connectivity, vert_no, edge, cell, i, j
    LOGICAL :: found


    CALL message("",'Checking connectivity for vertices...')
    in_grid => get_grid(in_grid_id)
    cells => in_grid%cells
    edges => in_grid%edges
    verts => in_grid%verts
    max_connectivity = verts%max_connectivity


!     WRITE(0,*) "-----------------------------"
!     WRITE(0,*) " Checking geometry of ", verts%no_of_existvertices, " verts."

    DO vert_no=1, verts%no_of_existvertices
      DO i=1,max_connectivity
      
        ! check edges
        edge = verts%get_edge_index(vert_no, i)
        IF (edge > 0) THEN
          IF (edges%get_vertex_index(edge,1) /= vert_no .AND. &
            & edges%get_vertex_index(edge,2) /= vert_no) THEN
            WRITE(0,*) "vertex:", vert_no, " edge:", edge
            WRITE(0,*) "vertex edges:", verts%get_edge_index(vert_no, :)
            WRITE(0,*) "edge verts:",   edges%get_vertex_index(edge,:)
            CALL finish("check_grid_connectivity_verts", "vertex not found in edges")
          ENDIF
        ENDIF
        
        ! check cells
        cell = verts%get_cell_index(vert_no, i)
        IF (cell > 0) THEN
          found = .false.
          DO j=1,cells%max_no_of_vertices
            IF (cells%get_vertex_index(cell,j) == vert_no) THEN
              found = .true.
              EXIT
            ENDIF
          ENDDO
          IF (.not. found) THEN
            write(0,*) vert_no, cell, ' vertex was not found in cell.'
            CALL finish('check_grid_connectivity_verts', 'Wrong vertex connectivity')
          ENDIF        
        ENDIF
        
        
      ENDDO !i=1,max_connectivity
    ENDDO

  END SUBROUTINE check_grid_connectivity_verts
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE  check_edge_incell(in_grid_id,edge_no,cell_no)
    INTEGER, INTENT(in) :: in_grid_id,edge_no,cell_no

    TYPE(t_grid_cells), POINTER :: cells
    INTEGER :: i

    cells => get_cells(in_grid_id)
    DO i=1,cells%max_no_of_vertices
      IF (cells%get_edge_index(cell_no,i) == edge_no) &
        RETURN
    ENDDO
    print*, 'cell=',cell_no
    print*, 'with edges=',cells%get_edge_index(cell_no,:)
    print*, ' does not contain edge_no=', edge_no
    CALL finish('check_edge_incell failed','')
  END SUBROUTINE check_edge_incell
  !-------------------------------------------------------------------------

END MODULE mo_grid_checktools

