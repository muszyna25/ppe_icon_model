!>
!! Perform grid optimization using a combination of optimizations:
!!
!! A. Edge spring dynamics (Tomita,H. et al., 2001: Shallow Water Model on a Modified
!!   Icosahedral Geodesic Grid by Using Spring Dynamics, JCP 174, 579-613
!!   Tomita,H. et al., 2002: An Optimization of the Icosahedral Grid Modified by
!!   Spring Dynamics
!!   (use_edge_springs)
!!
!! B. Adaptive/local edge spring reference lentgh (improves smoothness)
!!    (use_adaptive_spring_length, use_local_reference_length)
!!
!! C. Cell-centers to vertices springs. Gives a "dual" force.
!!    (dual_use_spring_cellcenters, use_prime_spring_cellcenters)
!!    Keeps cell areas ratio small and redistributes the local edge force
!!
!! D. Cell shape edge-spring reference length correction. 
!!    (use_centers_spring_correction) 
!!    The calculation of the edge-spring reference length is adapted according to the shape
!!    (isotropy) of the adjacent triangles.
!!
!! E. Barycenter forcing. (use_barycenter_force)
!!    Vertexes are moved close to the barycenter of the dual.
!!
!! F. Isotropy forcing. (use_isotropy_force)
!!    Aims to improve shape quality and smootheness
!!
!! @par Revision History
!! Initial Release by Almut Gassmann (2008-08-28)
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
MODULE mo_local_grid_optimization
#include "grid_definitions.inc"

  !-------------------------------------------------------------------------
  USE mo_kind,           ONLY: wp
  USE mo_exception,      ONLY: message_text, message, finish
  USE mo_io_units,       ONLY: nnml, filename_max
  USE mo_namelist,       ONLY: position_nml, open_nml, positioned
  USE mo_math_constants, ONLY: pi
  USE mo_local_grid_geometry,ONLY: get_triangle_circumcenters, get_cell_barycenters,&
    & use_cartesian_centers
  USE mo_io_local_grid,  ONLY: read_new_netcdf_grid, write_netcdf_grid
  USE mo_grid_toolbox,   ONLY: get_basic_dual_grid
  USE mo_timer,          ONLY: new_timer, timer_start, timer_stop, print_timer, delete_timer
  USE mo_local_grid
  USE mo_base_geometry
  USE mo_grid_conditions, ONLY: get_inner_vertices

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: version = '$Id$'

  PUBLIC :: optimize_grid, optimize_grid_file, read_grid_optimization_param
  ! PUBLIC :: optimize_grid_set_lengths

  INTEGER, PARAMETER :: max_min_condition_reached = 4
  INTEGER :: max_iterations

  REAL(wp) :: prime_ref_length_coeff, centers_ref_length_coeff
  REAL(wp) :: minedge_ref_length_weight

  REAL(wp) :: edge_ratio_condition
  REAL(wp) :: p_total_force_condition, p_max_force_condition, max_min_condition
  REAL(wp) :: centers_vertex_condition

  REAL(wp) :: prime_hard_spring_stiffness, prime_soft_spring_stiffness
  REAL(wp) :: spring_friction
  REAL(wp) :: centers_spring_stiffness
  
  REAL(wp) :: local_reference_length_coeff
  REAL(wp) :: centers_springcorrection_coeff
  REAL(wp) :: no_of_cell_edges
  
  REAL(wp) :: spring_dt
  REAL(wp) :: max_dt_distance_ratio
  REAL(wp) :: barycenter_force_coeff
  REAL(wp) :: isotropy_rotation_coeff, isotropy_stretch_coeff
  
  LOGICAL  :: use_optimization
  LOGICAL  :: use_prime_spring_cellcenters, use_dual_spring_cellcenters
  LOGICAL  :: use_adaptive_dt, use_adaptive_spring_length, use_local_reference_length
  LOGICAL  :: use_isotropy_force, use_barycenter_force
  LOGICAL  :: use_edge_springs, use_centers_spring_correction

  INTEGER :: optimize_vertex_depth

  ! R refinement parameters
  INTEGER, PARAMETER :: R_refine_none=0
  INTEGER, PARAMETER :: R_refine_springedges_linear=1
  
  ! R refinement namelist
  INTEGER  :: R_refine_method
  REAL(wp) :: R_refine_center_lon, R_refine_center_lat
  REAL(wp) :: R_refine_ratio, R_refine_flat_radius
  
  ! R refinement variables
  TYPE(t_cartesian_coordinates) :: R_refine_center
  


  CHARACTER(LEN=filename_max) :: input_file, output_file

CONTAINS

  !-------------------------------------------------------------------------
! #define compute_edge_length arc_length
! #define compute_edge_length norma
! #define compute_edge_length arc_length_normalsphere
#define compute_edge_length d_norma_3d
#define compute_arc_length(v) d_arc_of_hord_normalsphere(d_norma_3d(v))
#define compute_arc_of_length d_arc_of_hord_normalsphere
        
  !-------------------------------------------------------------------------
  FUNCTION get_tangent_unit(in_vector, sphere_point) result(tan_unit_vector)
    TYPE(t_cartesian_coordinates), INTENT(in) :: in_vector, sphere_point
    TYPE(t_cartesian_coordinates) :: tan_unit_vector

    tan_unit_vector%x =  in_vector%x - &
    & DOT_PRODUCT( in_vector%x, sphere_point%x) * sphere_point%x
    d_normalize(tan_unit_vector)
        
  END FUNCTION get_tangent_unit
  
  FUNCTION get_tangent_unit_x(in_vector, sphere_point) result(x)
    TYPE(t_cartesian_coordinates), INTENT(in) :: in_vector, sphere_point
    REAL(wp) :: x(3)

    x =  in_vector%x - &
    & DOT_PRODUCT( in_vector%x, sphere_point%x) * sphere_point%x
    d_normalize_x(x)
    
        
  END FUNCTION get_tangent_unit_x
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Spring dynamics adapted from NICAM.
  !!
  !! Abortion and exit criteria are slightly
  !! changed compared to that. Here, it is measured whether the potential energy
  !! might create new maxima, which means that the method is unstable and the
  !! 'beta' must be decreased. If the method is stable, the exit criterion looks
  !! for that state of the system, when kinetic energy has dropped to one tenth
  !! of the maximum value.
  !! Here: d= 2pi/ sqrt(number of triangles) instead of in NICAM: d=2 pi/10*2^(l-1)
  !!
  !! @par Revision History
  !! Initial Release by Almut Gassmann (2008-09-28)
  !! Modification by Almut Gassmann (2009-02-02)
  !! - spring_ref_length_coeff comes via the namelist input
  !!

  !---------------------------------------------------------
  !! SUBROUTINE read_grid_optimization_param
  !>
  !! Reads the parameters for local grid optimization
  !-------------------------------------------------------------------------
  SUBROUTINE read_grid_optimization_param(param_file_name)

    CHARACTER(LEN=*), INTENT(in) :: param_file_name

    TYPE(t_geographical_coordinates) :: geocoord
    INTEGER :: i_status

    NAMELIST /grid_optimization/ use_optimization, &
      & spring_dt, p_total_force_condition, p_max_force_condition, max_min_condition, &
      & spring_friction, input_file,  minedge_ref_length_weight,  &
      & output_file, use_prime_spring_cellcenters, use_dual_spring_cellcenters, &
      & centers_spring_stiffness, &
      & prime_hard_spring_stiffness, prime_soft_spring_stiffness, use_adaptive_dt, &
      & max_dt_distance_ratio, use_adaptive_spring_length,           &
      & prime_ref_length_coeff, centers_ref_length_coeff,                &
      & use_local_reference_length, &
      & centers_vertex_condition, use_isotropy_force,  &
      & local_reference_length_coeff, use_barycenter_force, barycenter_force_coeff, &
      & isotropy_rotation_coeff, isotropy_stretch_coeff, use_edge_springs,    &
      & centers_springcorrection_coeff, use_centers_spring_correction, &
      & R_refine_method, R_refine_center_lon, R_refine_center_lat, R_refine_ratio, &
      & R_refine_flat_radius, optimize_vertex_depth, max_iterations


    ! set default values
    max_iterations = 10000
    use_optimization = .false.
    prime_ref_length_coeff = 1.0_wp
    minedge_ref_length_weight = 0.25_wp
    centers_ref_length_coeff = 1.0_wp

    prime_hard_spring_stiffness = 1.0_wp
    prime_soft_spring_stiffness = 1.0_wp
    centers_spring_stiffness = 0.9_wp
    local_reference_length_coeff = 0.0_wp
    centers_springcorrection_coeff = 0.0_wp
    barycenter_force_coeff = 0.0_wp
    isotropy_rotation_coeff = 0.0_wp
    isotropy_stretch_coeff = 0.0_wp
    spring_friction = 1.0_wp

    spring_dt   = 1.6e-2_wp
    max_dt_distance_ratio = 0.75_wp

    p_total_force_condition = 0.5_wp
    p_max_force_condition = 0.25_wp
    centers_vertex_condition = 1.0_wp
    max_min_condition = 1.0_wp

    use_prime_spring_cellcenters = .false.
    use_dual_spring_cellcenters = .false.
    use_adaptive_dt = .false.
    use_adaptive_spring_length = .false.
    use_local_reference_length = .false.
    use_isotropy_force = .false.
    use_barycenter_force = .false.
    use_edge_springs = .true.

    input_file = ''
    output_file = ''
    optimize_vertex_depth = 1

    R_refine_method=R_refine_none
    R_refine_ratio=0.5_wp
    R_refine_flat_radius=0.05_wp
    
    
    ! read namelist
    CALL open_nml(param_file_name)
    CALL position_nml('grid_optimization',STATUS=i_status)
    IF (i_status == positioned) THEN
      READ (nnml,grid_optimization)
    ELSE
      WRITE(message_text,'(a)') " File", param_file_name, " not POSITIONED"
      CALL finish ('read_grid_optimization_param', message_text)
    ENDIF
    CLOSE(nnml)

    IF (R_refine_method /= R_refine_none) THEN
      geocoord%lon = R_refine_center_lon
      geocoord%lat = R_refine_center_lat
      R_refine_center = gc2cc(geocoord)
    ENDIF

    IF (optimize_vertex_depth < 1) &
      CALL finish('grid_optimization',&
        & 'optimize_vertex_depth must be > 0')
    
    WRITE(message_text,'(a)') "===================================="
    CALL message ('', TRIM(message_text))
    CALL message ('', '   Grid Optimization')
    WRITE(message_text,*) 'use_prime_spring_cellcenters=', use_prime_spring_cellcenters
    CALL message ('', TRIM(message_text))
    WRITE(message_text,*) 'use_adaptive_spring_length=', use_adaptive_spring_length
    CALL message ('', TRIM(message_text))
    WRITE(message_text,*) 'use_isotropy_force=', use_isotropy_force
    CALL message ('', TRIM(message_text))
    WRITE(message_text,*) 'use_adaptive_dt=', use_adaptive_dt
    CALL message ('', TRIM(message_text))
    WRITE(message_text,*) 'use_barycenter_force=', use_barycenter_force
    CALL message ('', TRIM(message_text))
    WRITE(message_text,*) 'barycenter_force_coeff=', barycenter_force_coeff
    CALL message ('', TRIM(message_text))
    WRITE(message_text,*) 'isotropy_rotation_coeff=', isotropy_rotation_coeff
    CALL message ('', TRIM(message_text))
    WRITE(message_text,'(a,f8.4)') 'max_dt_distance_ratio=', max_dt_distance_ratio
    CALL message ('', TRIM(message_text))
    WRITE(message_text,'(a,f8.4)') 'prime_ref_length_coeff=', prime_ref_length_coeff
    CALL message ('', TRIM(message_text))
    WRITE(message_text,'(a,f8.4)') 'centers_ref_length_coeff=', &
      & centers_ref_length_coeff
    CALL message ('', TRIM(message_text))
    WRITE(message_text,'(a,f8.4)') 'prime_soft_spring_stiffness=', prime_soft_spring_stiffness
    CALL message ('', TRIM(message_text))
    WRITE(message_text,'(a,f8.4)') 'prime_hard_spring_stiffness=', prime_hard_spring_stiffness
    CALL message ('', TRIM(message_text))
    WRITE(message_text,'(a,f8.4)') 'centers_spring_stiffness=', &
      & centers_spring_stiffness
    CALL message ('', TRIM(message_text))
    WRITE(message_text,'(a,f8.4)') 'local_reference_length_coeff=', &
      & local_reference_length_coeff
    CALL message ('', TRIM(message_text))
    WRITE(message_text,*) 'use_centers_spring_correction=', &
      & use_centers_spring_correction
    CALL message ('', TRIM(message_text))
    WRITE(message_text,'(a,f8.4)') 'centers_springcorrection_coeff=', &
      & centers_springcorrection_coeff
    CALL message ('', TRIM(message_text))
    WRITE(message_text,'(a,f8.4)') 'spring_friction=', spring_friction
    CALL message ('', TRIM(message_text))
     WRITE(message_text,'(a,f8.4)') 'spring_dt=', spring_dt
    CALL message ('', TRIM(message_text))
    WRITE(message_text,'(a,f8.4)') 'p_total_force_condition=', p_total_force_condition
    CALL message ('', TRIM(message_text))
    WRITE(message_text,'(a,f8.4)') 'p_max_force_condition=', p_max_force_condition
    CALL message ('', TRIM(message_text))
    WRITE(message_text,'(a,f8.4)') 'max_min_condition=', max_min_condition
    CALL message ('', TRIM(message_text))
    WRITE(message_text,'(a,f8.4)') 'centers_vertex_condition=', centers_vertex_condition
    CALL message ('', TRIM(message_text))
    WRITE(message_text,'(a,a)')    'input_file=', TRIM(input_file)
    CALL message ('', TRIM(message_text))
    WRITE(message_text,'(a,a)')    'output_file=', TRIM(output_file)
    CALL message ('', TRIM(message_text))

  END SUBROUTINE read_grid_optimization_param
  !-------------------------------------------------------------------------

  !---------------------------------------------------------

  !---------------------------------------------------------
  SUBROUTINE optimize_grid_file(param_file_name)

    CHARACTER(LEN=*), INTENT(in) :: param_file_name

    INTEGER :: grid_id

    CALL read_grid_optimization_param(param_file_name)

    grid_id = read_new_netcdf_grid(input_file)

    CALL optimize_grid(grid_id)

    CALL write_netcdf_grid(grid_id, output_file)
    CALL delete_grid(grid_id)

  END SUBROUTINE optimize_grid_file
  !---------------------------------------------------------

  !---------------------------------------------------------
  SUBROUTINE optimize_grid(grid_id, depth_level)

    INTEGER, INTENT(inout) :: grid_id
    INTEGER, INTENT(in), OPTIONAL :: depth_level

    INTEGER :: opt_result
    INTEGER :: timer_optimize_grid

    !--------------------------------------------------------------
    IF (.NOT. use_optimization) THEN
      CALL message('optimize_grid','No optimization')
      RETURN
    ENDIF
    
    timer_optimize_grid = new_timer("optimize_grid")
    CALL timer_start(timer_optimize_grid)
      
    opt_result = optimize_grid_methods(grid_id)
    CALL timer_stop(timer_optimize_grid)
    CALL print_timer(timer_optimize_grid)
    CALL delete_timer(timer_optimize_grid)
      
    RETURN
          
  END SUBROUTINE optimize_grid
  !---------------------------------------------------------

  !-------------------------------------------------------------------------
  INTEGER FUNCTION optimize_grid_methods(grid_id) result(opt_result)

    INTEGER, INTENT(inout) :: grid_id
!     INTEGER, INTENT(in), OPTIONAL :: depth_level
    
    TYPE(t_grid), POINTER :: in_grid
    TYPE(t_grid_cells), POINTER :: cells
    TYPE(t_grid_edges), POINTER :: edges
    TYPE(t_grid_vertices), POINTER :: verts

    TYPE(t_integer_list)  :: inner_verts_list

    INTEGER :: no_of_cells, no_of_edges, no_of_vertices
    INTEGER :: max_cell_vertices, max_vertex_connect
    INTEGER :: edge, cell, vertex, start_vertex, end_vertex
!     INTEGER :: min_edge_no, max_edge_no
   
    INTEGER :: iteration, j

    REAL (wp) :: dt, new_dt, old_dt
    REAL (wp) :: friction
    REAL (wp) :: ekin, total_force, oldtotal_force, maxekin, maxtotal_force
      
    REAL(wp) :: orientation, force, max_force, max_max_force, old_force
    REAL(wp) :: max_force_1, ekin_1, total_force_1
    REAL(wp) :: max_velocity, old_velocity

    REAL (wp) :: spring_ref_length, global_ref_length
    
    REAL(wp) :: length, average_length, min_length, max_length, centers_ref_length
    REAL(wp) :: triangle_area
    
    REAL(wp) :: max_dt_distance
    
    REAL(wp) :: spring_stiffness

    REAL(wp), ALLOCATABLE :: edge_force(:)

    REAL(wp) :: spring_dt_force_coeff, spring_dt_velocity_coeff
    REAL(wp) :: force_ceoff, velocity_coeff
    REAL(wp) :: spring_ref_length_coeff, half_ref_length_coeff

    REAL(wp) :: soft_spring_stiffness, hard_spring_stiffness
    REAL(wp) :: s_max_force_condition
    REAL(wp) :: s_total_force_condition
    ! REAL(wp), PARAMETER :: inv_sqrt_3 = 1.0_wp/SQRT(3.0_wp)
    ! REAL(wp), PARAMETER :: sqrt_2 = SQRT(2.0_wp)
    TYPE(t_cartesian_coordinates), ALLOCATABLE :: vertex_velocity(:), vertex_force(:)
    TYPE(t_cartesian_coordinates), ALLOCATABLE :: cell_isotropy_vector(:)
!     TYPE(t_cartesian_coordinates), ALLOCATABLE :: cell_vertex_centers_vector(:,:)
    TYPE(t_cartesian_coordinates), POINTER :: cell_center
!     TYPE(t_cartesian_coordinates), ALLOCATABLE :: cell_centers_vector(:)
!    TYPE(cartesian_coordinates) :: cell_center, vertex_vector, spring_vector
    TYPE(t_cartesian_coordinates) :: vertex_centers_vector(16)  ! should be >= max_vertex_connect
    TYPE(t_cartesian_coordinates) :: vertex_cell_force(16)
    REAL(wp) :: vertex_centers_length(16)

    TYPE(t_cartesian_coordinates) :: edge_vector, vertex_vector, force_vector
    REAL(wp),  ALLOCATABLE :: vertex_edge_ref_length(:)
    REAL(wp),  ALLOCATABLE :: cell_center_ref_length(:)
    REAL(wp),  ALLOCATABLE :: edge_ref_length(:)
    REAL(wp),  ALLOCATABLE :: cell_edge_ref(:)
!    REAL(wp) :: ref_length
    REAL(wp) :: diff_length

    INTEGER :: vertex_edges, vertex_list_idx

    INTEGER :: dual_grid_id, istat
    TYPE(t_grid), POINTER :: dual_grid
    TYPE(t_cartesian_coordinates), POINTER :: dual_barycenters(:)

    REAL(wp), POINTER :: vertex_ref_length(:)
    REAL(wp) :: R_refine_ratio_coeff

    LOGICAL :: is_triangle_grid!, use_cell_center_correction
    LOGICAL :: use_R_refine, use_vertex_ref_length, use_cell_centers
    LOGICAL :: use_edge_ref_length
        
    CHARACTER(*), PARAMETER :: method_name = "optimize_grid_methods"
    ! INTEGER :: dual_grid_id

    opt_result = 0
    in_grid  => get_grid(grid_id)
    verts=>in_grid%verts
    edges=>in_grid%edges
    cells=>in_grid%cells
    no_of_vertices = verts%no_of_existvertices
    no_of_edges = edges%no_of_existedges
    no_of_cells = cells%no_of_existcells
    max_cell_vertices  = cells%max_no_of_vertices
    max_vertex_connect = verts%max_connectivity
!     IF (PRESENT(depth_level)) THEN
!       start_vertex = verts%start_idx(depth_level,1)
!       end_vertex   = verts%end_idx(depth_level,1)
!     ELSE
    start_vertex = 1
    end_vertex   = no_of_vertices
!     ENDIF

    inner_verts_list=get_inner_vertices(grid_id,optimize_vertex_depth)
!     write(*,*) 'no_of_vertices:', no_of_vertices, inner_verts_list%list_size

    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    ! Calculate some conditions    
    is_triangle_grid = max_cell_vertices == 3
    use_R_refine = (R_refine_method == R_refine_none)
    !---------------
    IF (use_dual_spring_cellcenters .AND. use_prime_spring_cellcenters) THEN
      CALL finish(method_name, &
        & "Cannot use both use_dual_spring_cellcenters .AND. use_prime_spring_cellcenters")
    ENDIF
    !---------------
    use_edge_ref_length = use_prime_spring_cellcenters .OR. use_edge_springs
    use_vertex_ref_length = use_local_reference_length .OR. use_R_refine
    use_cell_centers = use_isotropy_force .OR. use_prime_spring_cellcenters
!     use_cell_center_correction = use_isotropy_force .OR. use_centers_spring_correction
    !-----------------------------------------------------------------------
    
    ! allocate auxiliary vectors
    ALLOCATE(edge_force(no_of_edges), &
      & vertex_force(no_of_vertices), &
      & vertex_velocity(no_of_vertices), stat=istat)
    IF (istat > 0) THEN
        CALL finish (method_name, 'ALLOCATE(vertex_force,vertex_velocity')
    ENDIF

    IF (use_edge_ref_length) THEN
      ALLOCATE(edge_ref_length(no_of_edges), stat=istat)
      IF (istat > 0) THEN
        CALL finish (method_name, 'ALLOCATE(edge_ref_length')
      ENDIF
    ENDIF
    

    IF (use_local_reference_length) THEN
      ALLOCATE(vertex_edge_ref_length(no_of_vertices), stat=istat)
      IF (istat > 0) &
        CALL finish (method_name, &
          & 'ALLOCATE(vertex_edge_ref_length, ...')
    ENDIF
    
    IF (use_prime_spring_cellcenters) THEN
      ALLOCATE(cell_center_ref_length(no_of_cells), &
        cell_edge_ref(max_cell_vertices), stat=istat)
      IF (istat > 0) CALL finish (method_name, &
          & 'ALLOCATE(cell_center_ref_length')
    ENDIF

    IF (use_isotropy_force) THEN
      ALLOCATE(cell_isotropy_vector(no_of_cells), stat=istat)
      IF (istat > 0) &
        CALL finish (method_name, &
          & 'ALLOCATE(cell_isotropy_vector)')
    ENDIF
    
!     IF (use_cell_center_correction) THEN
!       ALLOCATE(cell_vertex_centers_vector(max_vertex_connect,no_of_cells), stat=istat)
!       IF (istat > 0)&
!         CALL finish (method_name, &
!           & 'ALLOCATE(cell_vertex_centers_vector')
!     ENDIF
    
    IF (use_vertex_ref_length) THEN
      ALLOCATE(vertex_ref_length(no_of_vertices), stat=istat)
      IF (istat > 0)  CALL finish (method_name, &
          & 'ALLOCATE(vertex_ref_length)')
    ENDIF

    
     !-----------------------------------------------------------------------
    !lambda   = 2.0_wp*pi/(SQRT(REAL(no_of_cells,wp)*1.0_wp))
    !spring_ref_length     = lambda
    friction = spring_friction
    dt       = spring_dt
    !-----------------------------------------------------------------------
    ! Compute R_refine constants
    IF (R_refine_method == R_refine_springedges_linear) THEN
      IF (R_refine_flat_radius > 0.9_wp ) &
        CALL finish(method_name, 'R_refine_flat_radius > 0.9_wp')
      R_refine_ratio_coeff =  &
        (1.0_wp - R_refine_ratio)/ (1.0_wp - R_refine_flat_radius)
    ENDIF    
    !-----------------------------------------------------------------------
    s_max_force_condition   = p_max_force_condition
    s_total_force_condition = p_total_force_condition
    spring_ref_length_coeff = prime_ref_length_coeff
    soft_spring_stiffness   = prime_soft_spring_stiffness
    hard_spring_stiffness   = prime_hard_spring_stiffness

    !-----------------------------------------------------------------------
    ! compute spring_ref_length
    average_length = 0.0_wp
    min_length = 1E100_wp
    max_length = 0.0_wp
    DO edge = 1, no_of_edges

      edge_vector%x = &
        & verts%cartesian(edges%get_vertex_index(edge,2))%x - &
        & verts%cartesian(edges%get_vertex_index(edge,1))%x

      length = &
        compute_edge_length(edge_vector)
      edges%primal_edge_length(edge) = compute_arc_of_length(length)
      edges%cartesian_dual_normal(edge)%x = edge_vector%x / length

    ENDDO ! edge = 1, no_of_edges
    min_length = MINVAL(edges%primal_edge_length(1:no_of_edges))
    max_length = MAXVAL(edges%primal_edge_length(1:no_of_edges))
    average_length = SUM(edges%primal_edge_length(1:no_of_edges))/ REAL(no_of_edges,wp)

    half_ref_length_coeff = spring_ref_length_coeff * 0.5_wp
    global_ref_length = (minedge_ref_length_weight * min_length) + &
      & ((1.0_wp-minedge_ref_length_weight) * average_length)
    spring_ref_length = global_ref_length * spring_ref_length_coeff
!    spring_ref_length  = min_length
    centers_ref_length = average_length * centers_ref_length_coeff

    !-----------------------------------------------------------------------
    ! initial velocities of vertices is zero
    DO vertex = 1, no_of_vertices
      vertex_velocity(vertex)%x = 0.0_wp
    ENDDO

    !-----------------------------------------------------------------------
    ! initial values: force_ceoff = dt/2, changes after step 2
    force_ceoff = dt * 0.5_wp
    velocity_coeff = 0.0_wp

    !-----------------------------------------------------------------------
    ! time stepping values
    spring_dt_force_coeff = EXP(- friction * dt * 0.5_wp)
    spring_dt_velocity_coeff = spring_dt_force_coeff * spring_dt_force_coeff
    spring_dt_force_coeff = spring_dt_force_coeff * dt

    !-----------------------------------------------------------------------
    maxekin  = 0.0_wp
    maxtotal_force  = 0.0_wp
    max_max_force =  0.0_wp

    oldtotal_force  = 1E100_wp
    old_force = 1E100_wp
    old_velocity = 1E100_wp

!     write(*,*) 'Start iteration'
    !-----------------------------------------------------------------------
    DO iteration=1, max_iterations
!       write(*,*) "Iteration:",iteration
      !-----------------------------------------------------------------------
      ! get cell barycenters
      !CALL get_cell_barycenters(grid_id, center_type=use_cartesian_centers)
      IF (use_cell_centers) THEN 
        IF (is_triangle_grid) THEN
          CALL get_triangle_circumcenters(grid_id)
        ELSE
          CALL get_cell_barycenters(grid_id)
        ENDIF
      ENDIF
      !-----------------------------------------------------------------------


      !--------------------------------------------------------------------
      ! compute edge length and normal vectors
      average_length = 0.0_wp
      min_length = 1E100_wp
      max_length = 0.0_wp
      DO edge = 1, no_of_edges

        edge_vector%x = &
          & verts%cartesian(edges%get_vertex_index(edge,2))%x - &
          & verts%cartesian(edges%get_vertex_index(edge,1))%x

        length = &
          compute_edge_length(edge_vector)
!         edges%primal_edge_length(edge) = length
        edges%primal_edge_length(edge) = compute_arc_of_length(length)
        edges%cartesian_dual_normal(edge)%x = edge_vector%x / length

      ENDDO ! edge = 1, no_of_edges
      min_length = MINVAL(edges%primal_edge_length(1:no_of_edges))
      max_length = MAXVAL(edges%primal_edge_length(1:no_of_edges))
      average_length = SUM(edges%primal_edge_length(1:no_of_edges))/ REAL(no_of_edges,wp)

      !--------------------------------------------------------------------
      ! adaptive spring_ref_length
      max_dt_distance = (max_length - min_length) * max_dt_distance_ratio
      IF (use_adaptive_spring_length) THEN
        global_ref_length = (minedge_ref_length_weight * min_length) + &
          & ((1.0_wp-minedge_ref_length_weight) * average_length)
        spring_ref_length = global_ref_length * spring_ref_length_coeff
        centers_ref_length = average_length * centers_ref_length_coeff
      ENDIF
!       IF (max_length / min_length <= max_min_condition) THEN
!         opt_result = max_min_condition_reached
!         write(*,*) TRIM(method_name), " iterated ", iteration, " times."
!         CALL message(method_name, &
!           & 'max_min_condition reached. Exit')
!         EXIT
!       ENDIF
      
      !--------------------------------------------------------------------
!$OMP PARALLEL PRIVATE(vertex_centers_length, istat)
      ! compute edge_ref_length
      SELECT CASE (R_refine_method)
        CASE (R_refine_none)
        
          IF (use_vertex_ref_length) THEN
!$OMP DO PRIVATE(vertex)
            DO vertex = 1, no_of_vertices
              vertex_ref_length(vertex) = global_ref_length
            ENDDO
!$OMP ENDDO
          ENDIF
        CASE (R_refine_springedges_linear)
!$OMP DO PRIVATE(vertex)
          DO vertex=1,no_of_vertices
            ! vertex_ref_length gets a value from
            ! R_refine_ratio  to 1.0, x global_ref_length,
            ! as a linear function
            ! of the distance from the R_refine_center
            vertex_ref_length(vertex) = &
              (( MAX((arc_length_normalsphere(verts%cartesian(vertex), R_refine_center) / pi) &
               - R_refine_flat_radius, 0.0_wp) * R_refine_ratio_coeff) + R_refine_ratio )     &
               * global_ref_length
          
          ENDDO
!$OMP ENDDO
!           write(*,*) 'min,max edge_ref_length:', tmp_min, tmp_max,  tmp_min/tmp_max
        CASE default
          CALL finish(TRIM(method_name), "Uknown R_refine_method")          
      END SELECT
      
      
      !--------------------------------------------------------------------
      IF (use_local_reference_length) THEN

        ! compute edge reference length for each vertex
!$OMP DO PRIVATE(vertex, vertex_edges, j, edge)
        DO vertex=1,no_of_vertices
          vertex_edge_ref_length(vertex) = 0.0_wp
          vertex_edges = 0
          DO j=1,max_vertex_connect
            edge=verts%get_edge_index(vertex,j)
            IF (edge /= 0) THEN
              vertex_edges = vertex_edges + 1
              vertex_edge_ref_length(vertex) = vertex_edge_ref_length(vertex) &
                & + edges%primal_edge_length(edge)
            ENDIF              
          ENDDO
          vertex_edge_ref_length(vertex) = vertex_edge_ref_length(vertex) &
            & / REAL(vertex_edges,wp)
          ! update edge_reference
          vertex_ref_length(vertex) = &
            & (1.0_wp-local_reference_length_coeff) * vertex_ref_length(vertex) + &
            & local_reference_length_coeff * vertex_edge_ref_length(vertex)
        ENDDO
!$OMP ENDDO

      ENDIF !(use_local_reference_length)

        
      !-----------------------------------------------------------------------
      ! compute edge_ref_length
      IF (use_edge_ref_length) THEN
        IF (use_vertex_ref_length) THEN
!$OMP DO PRIVATE(edge)
          DO edge = 1, no_of_edges
            edge_ref_length(edge) = half_ref_length_coeff * &
              (vertex_ref_length(edges%get_vertex_index(edge,1)) + &
                vertex_ref_length(edges%get_vertex_index(edge,2)))
          ENDDO
!$OMP ENDDO
        ELSE
!$OMP DO PRIVATE(edge)
          DO edge = 1, no_of_edges
            edge_ref_length(edge) =  global_ref_length         
          ENDDO
!$OMP ENDDO
        ENDIF
      ENDIF ! use_edge_ref_length
      
      !-----------------------------------------------------------------------
!       IF (use_isotropy_force .OR. use_centers_spring_correction) THEN
!       ! compute cell_centers_vector (isotropy)
! !$OMP DO PRIVATE(cell, cell_center, j, vertex)
!         DO cell=1,no_of_cells
!           cell_center => cells%cartesian_center(cell)
!           cell_centers_vector(cell)%x = 0.0_wp
!           ! get cell_vertex_centers_vector, cell_centers_vector
!           DO j=1,max_cell_vertices
!             vertex=cells%get_vertex_index(cell,j)
!             IF (vertex > 0) THEN
!               cell_centers_vector(cell)%x = cell_centers_vector(cell)%x + &
!                 verts%cartesian(vertex)%x - cell_center%x             
!               cell_vertex_centers_vector(j,cell) = &
!                 normalize(verts%cartesian(vertex) - cell_center)
!             ENDIF
!           ENDDO !j=1,max_cell_vertices
!         ENDDO !cell=1,no_of_cells
! !$OMP ENDDO
!       ENDIF !(use_isotropy_force .OR. use_centers_spring_correction)
      

      !-----------------------------------------------------------------------
      ! If use_centers_spring_correction adjust the edge_ref_length
!       IF (use_centers_spring_correction) THEN
! !$OMP DO PRIVATE(edge, j, cell)
!         DO edge=1,no_of_edges
!           DO j=1,2
!             cell=edges%get_cell_index(edge,j)
!             IF (cell > 0) THEN
!                edge_ref_length(edge) = edge_ref_length(edge) - &
!                  norma_of_vector_product(cell_centers_vector(cell), &
!                  edges%cartesian_dual_normal(edge)) * centers_springcorrection_coeff
!             ENDIF
!           ENDDO
!         ENDDO
! !$OMP ENDDO
!       ENDIF

      !========================================================================
      ! FORCING
      !========================================================================
      
      ! zero vertex_force
!$OMP DO PRIVATE(vertex)
      DO vertex = start_vertex, end_vertex
        vertex_force(vertex)%x = 0.0_wp
      ENDDO
!$OMP ENDDO
      ! zero edge_force
!$OMP DO PRIVATE(edge)
      DO edge=1,no_of_edges
        edge_force(edge) = 0.0_wp
      ENDDO
!$OMP ENDDO
      
      !--------------------------------------------------------------------
      ! compute edge force from springs
      IF (use_edge_springs) THEN 
!         tmp_max=0._wp
!         tmp_min=99999._wp
!$OMP DO PRIVATE(edge, diff_length, spring_stiffness)
        DO edge = 1, no_of_edges
          diff_length =  edges%primal_edge_length(edge) - edge_ref_length(edge)
          spring_stiffness = MERGE(soft_spring_stiffness, hard_spring_stiffness, &
            & diff_length > 0.0_wp)
          edge_force(edge) = diff_length * spring_stiffness
          
!           IF (edge_ref_length(edge) /=spring_ref_length) THEN
!             write(*,*) edge_ref_length(edge), spring_ref_length
!             CALL finish('compute edge force from springs',&
!               'edge_ref_length(edge) /=spring_ref_length')
!           ENDIF
!           tmp_max=MAX(tmp_max,edge_ref_length(edge))
!           tmp_min=MIN(tmp_min,edge_ref_length(edge))
        ENDDO
!         write(*,*) 'min,max edge_ref_length:', tmp_min, tmp_max,  tmp_min/tmp_max
!$OMP ENDDO
      ENDIF
            
      !--------------------------------------------------------------------
      ! if requested use forcing from the centers of the cells
      IF (use_prime_spring_cellcenters) THEN
!         max_vertex_edge_ratio = 0.0_wp
!           write(*,*) " ======= use_prime_spring_cellcenters ========"

!$OMP DO PRIVATE(cell,no_of_cell_edges,j,edge, cell_edge_ref, length, triangle_area)
        DO cell=1,no_of_cells
          !--------------------------------------------
          ! this works only for triangles.
          cell_edge_ref(1) = edge_ref_length(cells%get_edge_index(cell,1))
          cell_edge_ref(2) = edge_ref_length(cells%get_edge_index(cell,2))
          cell_edge_ref(3) = edge_ref_length(cells%get_edge_index(cell,3))
          length = (cell_edge_ref(1) + cell_edge_ref(2) + cell_edge_ref(3)) * 0.5_wp
          triangle_area = SQRT(length * (length - cell_edge_ref(1)) * &
            & (length - cell_edge_ref(2)) * (length - cell_edge_ref(3)))
          
          cell_center_ref_length(cell) = &
            & (cell_edge_ref(1) * cell_edge_ref(2) * cell_edge_ref(3) * 0.25_wp &
            & * centers_ref_length_coeff) / triangle_area
          !--------------------------------------------
        ENDDO
! !$OMP ENDDO
        
        
!$OMP DO PRIVATE(vertex,vertex_vector,j, cell,vertex_cell_force,vertex_centers_vector, &
!$OMP  vertex_centers_length)
! vertex_edges,&
        DO vertex = start_vertex, end_vertex
!           max_vertex_edge = 0.0_wp
!           min_vertex_edge = 1E100_wp
!           vertex_edges = 0
          vertex_vector%x = verts%cartesian(vertex)%x
!           write(*,*) vertex, " verts%cartesian:", vertex_vector%x

          ! 1. The center-vertex length is measured on the sphere (ie. in great circles)
          ! 2. Each center-vertx force is moved (not projected) onto the tangent to
          !    the vertex plane
          DO j=1,max_vertex_connect
            cell = verts%get_cell_index(vertex,j)
            vertex_cell_force(j)%x = 0.0_wp
            
            IF (cell /= 0) THEN
                        
              vertex_centers_vector(j)%x = &
                & cells%cartesian_center(cell)%x - vertex_vector%x
! Unfortunately this line is taken as > 99 
              vertex_centers_length(j) = &
& compute_arc_length(vertex_centers_vector(j))
           
              vertex_cell_force(j)%x = &
                & get_tangent_unit_x(vertex_centers_vector(j), vertex_vector) * &
                & (vertex_centers_length(j) - cell_center_ref_length(cell)) *    &
                & centers_spring_stiffness
              
            ENDIF ! (cell /= 0) 
          ENDDO ! j=1,max_vertex_connect
          
          DO j=1,max_vertex_connect
            vertex_force(vertex)%x = vertex_force(vertex)%x + vertex_cell_force(j)%x
          ENDDO ! j=1,max_vertex_connect
          
        ENDDO! vertex = start_vertex, end_vertex                
!$OMP ENDDO
      ENDIF ! use_prime_spring_cellcenters
      !-----------------------------------------------------------------------
      
      !-----------------------------------------------------------------------
      ! If use_isotropy_force apply rotation and dual force
      IF (use_isotropy_force) THEN
!$OMP DO PRIVATE(cell, cell_center, j, vertex)
        DO cell=1,no_of_cells
          ! compute vertex_centers_vector and isotropy cell_centers_vector
          
          cell_center => cells%cartesian_center(cell)
          cell_isotropy_vector(cell)%x = 0.0_wp
          
          DO j=1,max_cell_vertices
            vertex=cells%get_vertex_index(cell,j)
            IF (vertex > 0) THEN
              cell_isotropy_vector(cell)%x =  cell_isotropy_vector(cell)%x + &
                & verts%cartesian(vertex)%x - cell_center%x  
            ENDIF
          ENDDO
        ENDDO ! cell=1,no_of_cells
!$OMP ENDDO


!$OMP DO PRIVATE(vertex,vertex_vector, j, cell, vertex_centers_vector, vertex_cell_force)
        ! compute rotation force for each vertwex
        DO vertex=start_vertex,end_vertex
          vertex_vector%x = verts%cartesian(vertex)%x
          
          DO j=1,max_vertex_connect
            cell=verts%get_cell_index(vertex,j)
!             vertex_cell_force(j)%x = 0.0_wp
            
            IF (cell /= 0) THEN
                            
              vertex_centers_vector(j) = &
                normalize(vertex_vector - cells%cartesian_center(cell))
                            
              vertex_cell_force(j)%x = vertex_centers_vector(j)%x * &
                (isotropy_rotation_coeff + isotropy_stretch_coeff) * &
                DOT_PRODUCT(vertex_centers_vector(j)%x,cell_isotropy_vector(cell)%x) -&
                cell_isotropy_vector(cell)%x  * isotropy_rotation_coeff
                
              ! add force tangent to the sphere at this vertex              
              vertex_force(vertex)%x = vertex_force(vertex)%x + &
                & (vertex_cell_force(j)%x - &
                & DOT_PRODUCT( vertex_cell_force(j)%x, vertex_vector%x) &
                & * vertex_vector%x)
            
            ENDIF ! (cell /= 0)
            
          ENDDO !j=1,max_cell_vertices
          
          ! compute stretching force
!           DO j=1,max_cell_vertices
!             edge=cells%get_edge_index(cell,j)
!             IF (edge > 0) THEN
!               IF ((edges%primal_edge_length(edge) - ref_length) &
!                    & * isotropy_stretch_coeff > 0.0_wp) &
!                 write(*,*) "Stretching:", &
!                 & (edges%primal_edge_length(edge) - ref_length) * isotropy_stretch_coeff
!               edge_force(edge) = edge_force(edge) + &
!                 & (edges%primal_edge_length(edge) - ref_length) * isotropy_stretch_coeff
!             ENDIF
!           ENDDO !j=1,max_cell_vertices
          
        ENDDO ! vertex=start_vertex,end_vertex
!$OMP ENDDO          
        
      ENDIF ! use_isotropy_force
!$OMP END PARALLEL
      !--------------------------------------------------------------------

      !--------------------------------------------------------------------
      
      !--------------------------------------------------------------------
      ! if requested also use the barycenter force
      IF (use_barycenter_force) THEN
        dual_grid_id = get_basic_dual_grid(grid_id)
        CALL get_cell_barycenters(dual_grid_id)
        dual_grid => get_grid(dual_grid_id)
        dual_barycenters => dual_grid%cells%cartesian_center

        DO vertex = start_vertex, end_vertex

          vertex_force(vertex)%x = vertex_force(vertex)%x + &
            (dual_barycenters(vertex)%x - verts%cartesian(vertex)%x) * barycenter_force_coeff

        ENDDO! vertex = start_vertex, end_vertex

        CALL delete_grid(dual_grid_id)
      ENDIF 

      !--------------------------------------------------------------------
      ! compute vertex force from grid edges
!       write(*,*) " ======= edge_spring ========"
      DO vertex = start_vertex, end_vertex
        DO j=1,max_vertex_connect
          edge = verts%get_edge_index(vertex,j)
          
          IF (edge /= 0) THEN
            orientation = MERGE(1.0_wp, -1.0_wp, edges%get_vertex_index(edge,1) == vertex)
        
            ! make force tangent to the sphere
            ! this is not a projection, the norm is invarient
!             force_vector = &
!               & get_tangent_unit(edges%cartesian_dual_normal(edge), verts%cartesian(vertex))
            
            vertex_force(vertex)%x = vertex_force(vertex)%x + &
              & edge_force(edge) *  orientation * &
              get_tangent_unit_x(edges%cartesian_dual_normal(edge), verts%cartesian(vertex))
            
          

          ENDIF ! (edge /= 0)

        ENDDO !(j=1,max_vertex_connect)
      ENDDO !vertex = start_vertex, end_vertex
!       write(*,*) " ======= end edge_spring ========"
      
      !--------------------------------------------------------------------
      !--------------------------------------------------------------------
      ! total vertex force
      max_force = 0.0_wp
      total_force = 0.0_wp
      DO vertex_list_idx=1,inner_verts_list%list_size
        vertex = inner_verts_list%value(vertex_list_idx)
        ! remove normal component from force
        ! force should not have a normal componenet from the construction
!         vertex_force(vertex)%x = vertex_force(vertex)%x - &
!           & DOT_PRODUCT(vertex_force(vertex)%x, verts%cartesian(vertex)%x) &
!           & * verts%cartesian(vertex)%x

        ! get max, total force
        ! force = SQRT(DOT_PRODUCT(vertex_force(vertex)%x, vertex_force(vertex)%x))
        force = DOT_PRODUCT(vertex_force(vertex)%x, vertex_force(vertex)%x)
        max_force = MAX(max_force, force)
        total_force = total_force + force
      ENDDO      
      total_force = SQRT(total_force) ! an approximation for checking conidtions
      max_force = SQRT(max_force)

      !--------------------------------------------------------------------
      ! compute adaptive dt if requested
      IF (use_adaptive_dt) THEN
         old_dt = dt
         new_dt = SQRT(2.0_wp * max_dt_distance / max_force)
         dt = MIN(spring_dt, new_dt)

!          print *, 'max/min distance:', max_length, min_length, max_length-min_length
!          print *, 'max_dt_distance:', max_dt_distance, max_dt_distance_ratio
!          print *, 'max_force:', max_force, max_dt_distance/max_force
!          print *, 'dt:', dt
!          print *, 'a. force_ceoff, velocity_coeff:', force_ceoff, velocity_coeff
!          print *,"dt, new dt:", dt, new_dt
        ! initial values
        force_ceoff = dt * 0.5_wp
        velocity_coeff = 0.0_wp
        IF (iteration > 1) THEN
          ! not initial values
          IF (dt == new_dt) THEN
!             print *,"use new dt:", dt
             force_ceoff = EXP(- friction * dt * 0.5_wp)
             velocity_coeff = force_ceoff * force_ceoff
             force_ceoff = force_ceoff * dt
          ELSE
            force_ceoff = spring_dt_force_coeff
            velocity_coeff = spring_dt_velocity_coeff
          ENDIF ! dt == new_dt
        ENDIF
!         print *, 'b. force_ceoff, velocity_coeff:', force_ceoff, velocity_coeff

      ENDIF !use_adaptive_dt

      !--------------------------------------------------------------------
      ! compute vertex velocities and new positions
      max_velocity = 0.0_wp
      ekin = 0.0_wp
      DO vertex_list_idx=1,inner_verts_list%list_size
        vertex = inner_verts_list%value(vertex_list_idx)

        ! solve for spring equation
        ! semi-implicit Stroermer-Verlet scheme
        vertex_velocity(vertex)%x = &
          vertex_velocity(vertex)%x * velocity_coeff + &
          vertex_force(vertex)%x * force_ceoff

        ! - > Horizontalize
        !vertex_velocity(vertex)%x = vertex_velocity(vertex)%x - &
        !  & DOT_PRODUCT(vertex_velocity(vertex)%x, verts%cartesian(vertex)%x) &
        !  & * verts%cartesian(vertex)%x

        ! solve the position equation and normalize (get vertex on the sphere)
        verts%cartesian(vertex)%x = verts%cartesian(vertex)%x + &
          & dt * vertex_velocity(vertex)%x
        verts%cartesian(vertex) = normalize(verts%cartesian(vertex))

        max_velocity = MAX(max_velocity,&
          & DOT_PRODUCT(vertex_velocity(vertex)%x, vertex_velocity(vertex)%x))

        ! kinetic energy for testing
        !ekin = 0.5_wp * DOT_PRODUCT(vertex_velocity(vertex)%x,&
        !  & vertex_velocity(vertex)%x) + ekin

        force_ceoff = spring_dt_force_coeff
        velocity_coeff = spring_dt_velocity_coeff

      ENDDO ! vertex_list_idx=1,inner_verts_list%list_size
      !--------------------------------------------------------------------

      IF (iteration==1) THEN
        ekin_1 = ekin
        maxekin = ekin
        total_force_1 = total_force
        maxtotal_force = total_force
        max_force_1 = max_force
        max_max_force = max_force
      ELSE
        maxekin = MAX(ekin,maxekin)
        maxtotal_force = MAX(total_force,maxtotal_force)
        max_max_force = MAX(max_max_force, max_force)
      ENDIF
 !     write(0,*) "1-current total,max force:",total_force_1,total_force,max_force_1,max_force

      ! print*, max_force, max_max_force

      ! print *, iteration,' potential :', total_force,maxtotal_force, ' kinetic:', ekin, maxekin
      ! print *, iteration,' max_force :', max_force, ' max_velocity:', max_velocity

      IF (iteration > 5 .AND. &
        max_force <= max_force_1 * s_max_force_condition .AND. &
        total_force <= total_force_1 * s_total_force_condition) THEN
          write(0,*) method_name," iterated ", iteration, " times."
          CALL message(TRIM(method_name), &
            & 'Potential conditions reached. Exit')
          EXIT
      ENDIF

!      IF (iteration > 5 .and. total_force == maxtotal_force) THEN
      IF (iteration > 5 .and. total_force > oldtotal_force) THEN
!        CALL finish('optimize_grid_methods', &
!          & 'potential energy increases. Please adjust spring_ref_length_coeff in the range 0.9..1.11')
        write(0,*) "Iteration ", iteration, ". Potential energy increases!"
!        EXIT
      ENDIF
      ! kinet energy should decrease to 1/10 of the maxium
!      IF (iteration > 5 .and. ekin < 0.1_wp*maxekin ) EXIT
      IF (iteration == max_iterations) THEN
        CALL finish(method_name,'no convergence in grid iteration')
      ENDIF

      oldtotal_force = total_force
      old_force = max_force
      old_velocity = max_velocity


    ENDDO ! iteration=1, max_iterations

    DEALLOCATE(inner_verts_list%value)
    DEALLOCATE(edge_force,vertex_force,vertex_velocity)

    IF (use_edge_ref_length) THEN
      DEALLOCATE(edge_ref_length)
    ENDIF
    
    IF (use_local_reference_length) THEN
      DEALLOCATE(vertex_edge_ref_length)
    ENDIF

    IF (use_isotropy_force) THEN
      DEALLOCATE(cell_isotropy_vector)
    ENDIF
    
!     IF (use_cell_centers) THEN
!       DEALLOCATE(cell_centers_vector)
!     ENDIF
    
    IF (use_prime_spring_cellcenters) THEN
      DEALLOCATE(cell_center_ref_length, cell_edge_ref)
    ENDIF
    
!     IF (use_cell_center_correction) THEN
!       DEALLOCATE(cell_vertex_centers_vector)
!     ENDIF
    
    IF (use_vertex_ref_length) THEN
      DEALLOCATE(vertex_ref_length)
    ENDIF

    RETURN
    
  END FUNCTION optimize_grid_methods


END MODULE mo_local_grid_optimization


