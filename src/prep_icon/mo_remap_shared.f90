MODULE mo_remap_shared

  USE mo_kind,               ONLY: wp
  USE mo_parallel_config,    ONLY: nproma
  USE mo_exception,          ONLY: finish
  USE mo_impl_constants,     ONLY: SUCCESS
  USE mo_communication,      ONLY: idx_1d, blk_no, idx_no
  USE mo_math_constants,     ONLY: pi_4, pi_2, pi2, pi_180
  USE mo_math_utilities,     ONLY: t_geographical_coordinates
  USE mo_model_domain,       ONLY: t_patch
  USE mo_gnat_gridsearch,    ONLY: gnat_tree, gnat_recursive_query,  &
    &                              gnat_std_radius, MAX_RANGE
  USE mo_remap_config,       ONLY: dbg_level, MAX_NAME_LENGTH

  IMPLICIT NONE

  PRIVATE

  ! subroutines:
  PUBLIC :: dist_deg, dist_cc
  PUBLIC :: compute_length
  PUBLIC :: edge2line
  PUBLIC :: cell2poly
  PUBLIC :: is_valid
  PUBLIC :: ccw
  PUBLIC :: intersect
  PUBLIC :: compute_intersection
  PUBLIC :: inside
  PUBLIC :: transform_lambert_azimuthal
  PUBLIC :: backtransform_lambert_azimuthal
  PUBLIC :: latitude_crossing
  PUBLIC :: compute_index_lists
  PUBLIC :: finalize_index_lists
  PUBLIC :: compute_coordinate_transform
  PUBLIC :: copy_vertex_coords
  PUBLIC :: clear_lookup_tables
  PUBLIC :: deallocate_lookup_tables
  PUBLIC :: allocate_lookup_tables
  PUBLIC :: finalize_vertex_coords
  PUBLIC :: get_containing_cell_generic
  PUBLIC :: get_containing_cell_gnat
  PUBLIC :: ccw_orientation
  PUBLIC :: gnat_recursive_query
  PUBLIC :: normalized_coord
  ! types and variables:
  PUBLIC :: t_line, t_poly
  PUBLIC :: t_regular_grid
  PUBLIC :: t_lookup_tbl
  PUBLIC :: t_index_list
  PUBLIC :: t_grid
  PUBLIC :: GRID_TYPE_REGULAR
  PUBLIC :: GRID_TYPE_ICON    
  PUBLIC :: LIST_NPOLE, LIST_SPOLE, LIST_DEFAULT, LIST_NAME
  PUBLIC :: pole_thresh
  PUBLIC :: npole, spole

  ! grid structure types:
  INTEGER, PARAMETER :: GRID_TYPE_REGULAR   = 1
  INTEGER, PARAMETER :: GRID_TYPE_ICON      = 2

  ! index lists (for better readibility):
  INTEGER, PARAMETER :: LIST_NPOLE          = 1
  INTEGER, PARAMETER :: LIST_SPOLE          = 2
  INTEGER, PARAMETER :: LIST_DEFAULT        = 3
  CHARACTER(LEN=7), PARAMETER :: LIST_NAME(3) = &
    & (/ "NPOLE  ", "SPOLE  ", "DEFAULT" /) 

  ! pole factor (constant used for Lambert transform):
  REAL(wp), PARAMETER :: r_ns(3) = (/ 1.0_wp, -1.0_wp, 0.0_wp /)

  !> north pole
  TYPE(t_geographical_coordinates), PARAMETER :: &
    &    npole = t_geographical_coordinates(0._wp,  90._wp)

  !> south pole
  TYPE(t_geographical_coordinates), PARAMETER :: &
    &    spole = t_geographical_coordinates(0._wp, -90._wp)

  !> threshold: cells within "pole_thresh*grid%char_length" to one of
  !  the poles are treated with Lambert transform:
  REAL(wp), PARAMETER :: pole_thresh = 2_wp

  TYPE t_line
    TYPE (t_geographical_coordinates) :: p(2)
  END TYPE t_line

  TYPE t_poly
    TYPE (t_geographical_coordinates) :: p(4)
  END TYPE t_poly

  !> Data structure for regular grid topologies
  TYPE t_regular_grid
    ! no. of points of the regular grid:
    INTEGER :: nxpoints, nypoints
    ! points (= "cell centers") of the regular grid:
    REAL(wp), ALLOCATABLE :: xvals1D(:), yvals1D(:)
    ! edge positions of the regular grid
    REAL(wp), ALLOCATABLE :: edge_lon(:), edge_lat(:)
    ! permutation array (needed for reordering)
    INTEGER, ALLOCATABLE  :: xperm1D(:), yperm1D(:)
  END TYPE t_regular_grid

  !> index lists for edges with at least one vertex within
  !  "edge_thresh" to north/south pole plus the complementary set
  !
  TYPE t_index_list
    LOGICAL :: l_initialized

    INTEGER               :: cnlist(3)
    INTEGER,  ALLOCATABLE :: clist_idx(:,:,:), clist_blk(:,:,:)
  END TYPE t_index_list

  !> Internal, temporary data (lookup tables)
  !
  TYPE t_lookup_tbl
    LOGICAL :: l_initialized
    ! index/block of containing cell (in this grid) for a given vertex
    ! (of a different grid): (idx,blk,coord_transform,thread)
    INTEGER, ALLOCATABLE :: vertex_c_idx(:,:,:,:), vertex_c_blk(:,:,:,:)
  END TYPE t_lookup_tbl


  !> Data structure for a general grid topology
  ! 
  !  Essentially, this derived type simply contains an instance of
  !  type "t_patch" for unstructured 2D meshes. However, the data
  !  structure applies also to structured rectangular grids where many
  !  performance optimizations are possible. Additional information is
  !  therefore required.
  !
  TYPE t_grid
    CHARACTER(LEN=MAX_NAME_LENGTH) :: name !< name string (for screen messages)

    INTEGER              :: structure
    TYPE(t_patch)        :: p_patch
    
    ! index lists for edges near north/south pole
    TYPE (t_index_list) :: index_list

    ! vertex coordinate list: This array contains a copy of
    ! "p_patch%verts%vertex" and a Lambert transform
    TYPE(t_geographical_coordinates), ALLOCATABLE :: vertex(:,:,:)

    ! regular grid data (not necessarily defined):
    TYPE(t_regular_grid) :: regular_grid

    ! Internal, temporary data (lookup tables)
    TYPE(t_lookup_tbl)   :: lookup_tbl

    ! characteristic grid size (deg)
    REAL(wp)             :: char_length

    ! index/block of vertex-neighbors (cells)
    INTEGER, ALLOCATABLE :: vertex_nb_idx(:,:,:), vertex_nb_blk(:,:,:)

    ! size of vertex-neighbor cell stencil
    INTEGER, ALLOCATABLE :: vertex_nb_stencil(:,:)

    ! For distributed coefficient computation, we create "patch
    ! coverings", i.e.  local grid partitions used as source grids,
    ! which are larger than the local grid partition.

    ! for "patch coverings": cell owner PE
    INTEGER, ALLOCATABLE :: owner_c(:)
    ! for "patch coverings": local cell index ( /= local index in covering!)
    INTEGER, ALLOCATABLE :: local_c(:)
    ! for grids having "patch coverings": local index in covering
    INTEGER, ALLOCATABLE :: cov_c(:)
  END TYPE t_grid


CONTAINS

  PURE FUNCTION dist_deg(p1, p2)
    REAL(wp) :: dist_deg
    TYPE (t_geographical_coordinates), INTENT(IN)  :: p1,p2
    ! local variables
    REAL(wp) :: v
    ! spherical distance:
    v = SIN(p1%lat*pi_180)*SIN(p2%lat*pi_180) &
      &  +  COS(p1%lat*pi_180)*COS(p2%lat*pi_180)*COS((p1%lon-p2%lon)*pi_180)
    dist_deg = ACOS( MIN(MAX(-1._wp, v), 1._wp) )
  END FUNCTION dist_deg


  PURE FUNCTION dist_cc(p1, p2)
    REAL(wp) :: dist_cc
    TYPE (t_geographical_coordinates), INTENT(IN)  :: p1,p2
    dist_cc = SQRT((p2%lon-p1%lon)*(p2%lon-p1%lon) + &
      &            (p2%lat-p1%lat)*(p2%lat-p1%lat))
  END FUNCTION dist_cc

  ! Compute length of a line segment.
  !
  ! This "length computation" may be just an approximation and must
  ! not necessarily identical to the "dist" function.
  !
  ELEMENTAL FUNCTION compute_length(l1) RESULT(length)
    REAL(wp) :: length
    TYPE (t_line), INTENT(IN) :: l1
    length = dist_deg(l1%p(1), l1%p(2))
  END FUNCTION compute_length

  !> Converts a grid edge to a "line" object in Euclidean space.
  !
  ! @note For treatments of lon-lat coords in Euclidean space, we have
  !       to make sure that we do not cross the +-180 degrees.
  !
  FUNCTION edge2line(grid, jc,jb, coord_transform) RESULT(line)
    TYPE(t_line) :: line
    TYPE(t_grid),   INTENT(IN) :: grid
    INTEGER,        INTENT(IN) :: jc,jb
    INTEGER,        INTENT(IN) :: coord_transform
    ! local variables
    INTEGER :: edge_vertex_idx(2), edge_vertex_blk(2)

    edge_vertex_idx(1:2) = grid%p_patch%edges%vertex_idx(jc, jb, 1:2)
    edge_vertex_blk(1:2) = grid%p_patch%edges%vertex_blk(jc, jb, 1:2)
    ! note: we use the (possibly transformed) vertex coordinates!
    line%p(1) = grid%vertex(edge_vertex_idx(1), edge_vertex_blk(1), coord_transform)
    line%p(2) = grid%vertex(edge_vertex_idx(2), edge_vertex_blk(2), coord_transform)

    IF ((MINVAL(line%p(:)%lon) < -90._wp) .AND. (MAXVAL(line%p(:)%lon) > 90._wp)) THEN
      WHERE (line%p(1:2)%lon < -90._wp)
        line%p(1:2)%lon = line%p(1:2)%lon + 360._wp
      END WHERE
    END IF
  END FUNCTION edge2line


  !> Converts a grid cell to a "polygon" object in Euclidean space.
  !
  ! @note For treatments of lon-lat coords in euclidean space, we have
  !       to make sure that we do not cross the +-180 degrees.
  !
  SUBROUTINE cell2poly(grid, jc,jb, coord_transform, poly, ne)
    TYPE(t_grid),   INTENT(IN)    :: grid
    INTEGER,        INTENT(IN)    :: jc,jb
    INTEGER,        INTENT(IN)    :: coord_transform
    TYPE(t_poly),   INTENT(OUT), TARGET :: poly
    INTEGER,        INTENT(OUT)   :: ne
    ! local variables
    REAL(wp), PARAMETER :: TOL = 1e-6
    INTEGER  :: vertex_idx(4), vertex_blk(4), i
    TYPE (t_geographical_coordinates) :: t
    REAL(wp) :: min_lon, max_lon, lon
    TYPE (t_geographical_coordinates), POINTER :: pp_ne

    ne = grid%p_patch%cell_type

    vertex_idx(1:ne) = grid%p_patch%cells%vertex_idx(jc,jb,1:ne)
    vertex_blk(1:ne) = grid%p_patch%cells%vertex_blk(jc,jb,1:ne)
    ! note: we use the (possibly transformed) vertex coordinates!
    ne = 1
    pp_ne => poly%p(ne)
    t = grid%vertex(vertex_idx(1), vertex_blk(1),coord_transform)
    pp_ne = t
    min_lon = t%lon
    max_lon = t%lon
    DO i=2,grid%p_patch%cell_type
      t = grid%vertex(vertex_idx(i), vertex_blk(i),coord_transform)
      ! skip degenerate edges
      IF ((ABS(pp_ne%lon - t%lon) > TOL) .OR.  &
        & (ABS(pp_ne%lat - t%lat) > TOL)) THEN
        ne = ne + 1
        pp_ne => poly%p(ne) 
        pp_ne = t

        min_lon = MIN(min_lon, t%lon)
        max_lon = MAX(max_lon, t%lon)
      END IF
    END DO

    IF ((min_lon < -90._wp) .AND. (max_lon > 90._wp)) THEN
      DO i=1,ne
        lon = poly%p(i)%lon
        IF (lon < -90._wp) THEN
          poly%p(i)%lon = lon + 360._wp
        END IF
      END DO
    END IF
  END SUBROUTINE cell2poly


  !> @return .TRUE. if this is a valid edge (with non-zero length)
  !
  !  @note   This routine must be *very* cheap. Therefore, one
  !          should probably resort to precomputed flags!
  !
  FUNCTION is_valid(line)
    LOGICAL :: is_valid
    TYPE (t_line), INTENT(IN) :: line
    ! local variables:
    REAL(wp), PARAMETER :: TOL = 1e-10
    ! skip degenerate edges
    is_valid = ((ABS(line%p(1)%lon - line%p(2)%lon) > TOL) .OR.  &
      &         (ABS(line%p(1)%lat - line%p(2)%lat) > TOL))
  END FUNCTION is_valid


  !>  test for counter-clockwise direction when travelling 
  !   from point p0 to p1 to p2.
  !
  !   @return -1 if clockwise direction or         
  !           if pts are collinear and p0 is between p1 and p2
  !   @return  0 if pts are collinear and p2 is between p0 and p1
  !   @return +1 if counter-clockwise direction or 
  !            if pts are collinear and p1 is between p0 and p2
  ELEMENTAL FUNCTION ccw(p0, p1, p2)
    INTEGER :: ccw
    TYPE (t_geographical_coordinates), INTENT(IN) :: p0,p1,p2
    ! local variables
    REAL(wp) :: dx1, dx2, dy1, dy2

    dx1 = p1%lon - p0%lon
    dy1 = p1%lat - p0%lat
    dx2 = p2%lon - p0%lon
    dy2 = p2%lat - p0%lat
    IF (dx1*dy2 > dy1*dx2) THEN
      ccw = +1
      RETURN
    END IF
    IF (dx1*dy2 < dy1*dx2) THEN
      ccw = -1
      RETURN
    END IF
    IF ((dx1*dx2 < 0) .OR. (dy1*dy2 < 0)) THEN
      ccw = -1
      RETURN
    END IF
    IF ((dx1*dx1 + dy1*dy1) < (dx2*dx2 + dy2*dy2)) THEN
      ccw = +1
      RETURN
    END IF
    ccw = 0;
  END FUNCTION ccw


  !> See, e.g., the chapter on elementary geometric methods in
  !  Sedgewick, Robert: Algorithms in C. Reading, MA: Addison Wesley Longman.
  !
  FUNCTION intersect(l1_in, l2_in)
    LOGICAL :: intersect
    TYPE (t_line), INTENT(INOUT) :: l1_in
    TYPE (t_line), INTENT(INOUT) :: l2_in
    ! local variables
    TYPE (t_line) :: l1, l2
    REAL(wp) :: lon1(2), lon2(2)
    INTEGER :: i

    l1 = l1_in
    l2 = l2_in

    lon1 = l1%p(1:2)%lon
    lon2 = l2%p(1:2)%lon
    IF (((MAXVAL(lon1) >  90._wp) .AND. (MINVAL(lon2) < -90._wp)) .OR. &
      & ((MAXVAL(lon2) >  90._wp) .AND. (MINVAL(lon1) < -90._wp))) THEN
      DO i=1,2
        IF (lon1(i) < -90._wp)  l1%p(i)%lon = lon1(i) + 360._wp
        IF (lon2(i) < -90._wp)  l2%p(i)%lon = lon2(i) + 360._wp
      END DO
    END IF

    intersect =  ((ccw(l1%p(1), l1%p(2), l2%p(1)) * ccw(l1%p(1), l1%p(2), l2%p(2))) <= 0)  .AND.   &
      &          ((ccw(l2%p(1), l2%p(2), l1%p(1)) * ccw(l2%p(1), l2%p(2), l1%p(2))) <= 0)

  END FUNCTION intersect


  SUBROUTINE compute_intersection(l1_in, l2_in, p) 
    TYPE (t_line), INTENT(IN) :: l1_in, l2_in
    TYPE (t_geographical_coordinates), INTENT(OUT) :: p
    ! local variables
    REAL(wp), PARAMETER :: TOL = 1e-10
    REAL(wp) :: dx1, dy1, dx2, dy2, dxa, det2, det11
    TYPE (t_line) :: l1, l2

    l1 = l1_in
    l2 = l2_in

    IF (((MAXVAL(l1%p(1:2)%lon) >  90._wp) .AND. (MINVAL(l2%p(1:2)%lon) < -90._wp)) .OR. &
      & ((MAXVAL(l2%p(1:2)%lon) >  90._wp) .AND. (MINVAL(l1%p(1:2)%lon) < -90._wp))) THEN
      WHERE (l1%p(1:2)%lon < -90._wp)
        l1%p(1:2)%lon = l1%p(1:2)%lon + 360._wp
      END WHERE
      WHERE (l2%p(1:2)%lon < -90._wp)
        l2%p(1:2)%lon = l2%p(1:2)%lon + 360._wp
      END WHERE
    END IF

    dx1=l1%p(2)%lon-l1%p(1)%lon
    dx2=l2%p(2)%lon-l2%p(1)%lon
    dxa=l2%p(1)%lon-l1%p(1)%lon
    IF (dx1 >  180.0_wp) dx1 = dx1 - 360._wp
    IF (dx1 < -180.0_wp) dx1 = dx1 + 360._wp
    IF (dx2 >  180.0_wp) dx2 = dx2 - 360._wp
    IF (dx2 < -180.0_wp) dx2 = dx2 + 360._wp
    IF (dxa >  180.0_wp) dxa = dxa - 360._wp
    IF (dxa < -180.0_wp) dxa = dxa + 360._wp

    dy1=l1%p(2)%lat-l1%p(1)%lat
    dy2=l2%p(2)%lat-l2%p(1)%lat

    ! line 1 parametrized as (x,y) = l1%p(1) + s1*(dx1,dx2)
    ! intersection: s1=det11/det2
    det2 = dx1*dy2 - dx2*dy1
  
    IF (ABS(det2) < TOL) THEN
      p = l1%p(1)
    ELSE
      det11 = dxa*dy2 - dx2*(l2%p(1)%lat-l1%p(1)%lat);
      p = l1%p(1)
      p%lon = p%lon + det11/det2*dx1
      p%lat = p%lat + det11/det2*dy1
    END IF
  END SUBROUTINE compute_intersection


  !> See, e.g., the chapter on elementary geometric methods in
  !  Sedgewick, Robert: Algorithms in C. Reading, MA: Addison Wesley Longman.
  !
  !  @todo enhance performance with special implementations for
  !        triangles and quads.
  !
  FUNCTION inside(t, npts, poly)
    LOGICAL :: inside
    INTEGER,                           INTENT(IN) :: npts
    TYPE (t_geographical_coordinates), INTENT(IN) :: t
    TYPE (t_poly),                     INTENT(IN) :: poly
    ! local variables:
    INTEGER :: p2,i,v, ones, m_ones

    ones   = 0
    m_ones = 0
    DO i=1,npts
      IF (i==npts) THEN
        p2=1 
      ELSE
        p2=i+1
      END IF
      v = ccw(t,poly%p(i),poly%p(p2))
      IF (v == +1)   ones =   ones + 1
      IF (v == -1) m_ones = m_ones + 1
    END DO
    inside = (ones*m_ones==0)
  END FUNCTION inside


  !> computes a Lambert equivalent azimuthal projection from
  !  spherical coordinates to the open disc with radius 2.
  !
  !  @note For technical reasons, the result data type is
  !        "t_geographical_coordinates", though it contains cartesian
  !        coordinate values.
  !
  ELEMENTAL FUNCTION transform_lambert_azimuthal(p_in, coord_transform) RESULT(x)
    TYPE (t_geographical_coordinates)             :: x    !< result in cartesian coords
    TYPE (t_geographical_coordinates), INTENT(IN) :: p_in    !< point in spherical coordinates (deg)
    INTEGER,                           INTENT(IN) :: coord_transform !< Flag, +1.: north pole, -1.: south pole
    ! local variables
    TYPE (t_geographical_coordinates) :: p
    REAL(wp)                          :: t, this_r_ns

    p = p_in
    this_r_ns = r_ns(coord_transform)
    t = 2._wp * SIN(this_r_ns*pi_4 -pi_180*p%lat/2._wp)
    x%lon = this_r_ns*t*COS(p%lon*pi_180)
    x%lat =           t*SIN(p%lon*pi_180)

  END FUNCTION transform_lambert_azimuthal


  !> computes a Lambert equivalent azimuthal projection from
  !  spherical coordinates to the open disc with radius 2.
  !
  !  @note For technical reasons, the result data type is
  !        "t_geographical_coordinates", though it contains Cartesian
  !        coordinate values.
  !
  FUNCTION backtransform_lambert_azimuthal(x, coord_transform) RESULT(p)
    TYPE (t_geographical_coordinates)             :: p    !< result in spherical coordinates (deg)
    TYPE (t_geographical_coordinates), INTENT(IN) :: x    !< point in cartesian coords
    INTEGER,                           INTENT(IN) :: coord_transform !< Flag, +1.: north pole, -1.: south pole
    ! local variables
    REAL(wp), PARAMETER :: TOL = 1.e-5_wp
    REAL(wp)            :: this_r_ns

    this_r_ns = r_ns(coord_transform)
    IF (ABS(x%lon) > TOL) THEN
      p%lon = this_r_ns*ATAN2(x%lat, x%lon)
      IF (p%lon < 0._wp) p%lon=p%lon+pi2
      p%lat = this_r_ns*pi_2 - 2._wp*ASIN(this_r_ns*x%lon/(2._wp*COS(p%lon)))
    ELSE IF (ABS(x%lat) > TOL) THEN
      p%lon = this_r_ns*ATAN2(x%lat, x%lon)
      IF (p%lon < 0._wp) p%lon=p%lon+pi2
      p%lat = this_r_ns*pi_2 - 2._wp*ASIN(x%lat/(2._wp*SIN(p%lon)))
    ELSE
      p%lat = this_r_ns*pi_2
    END IF
    p%lon = p%lon/pi_180
    p%lat = p%lat/pi_180
    p = normalized_coord(p) 
  END FUNCTION backtransform_lambert_azimuthal


  FUNCTION latitude_crossing(x1_deg,x2_deg,rlat, coord_transform)
    TYPE(t_geographical_coordinates) :: latitude_crossing
    TYPE(t_geographical_coordinates), INTENT(IN) :: x1_deg,x2_deg
    REAL(wp),                         INTENT(IN) :: rlat
    INTEGER,                          INTENT(IN) :: coord_transform
    ! local variables
    REAL(wp) :: p,q,s, n, this_r_ns, t
    TYPE(t_geographical_coordinates) :: d, x1, x2, res

    this_r_ns = r_ns(coord_transform)    
    t = 2._wp * SIN(this_r_ns*pi_4 -pi_180*rlat/2._wp)

    x1 = transform_lambert_azimuthal(x1_deg,coord_transform)
    x2 = transform_lambert_azimuthal(x2_deg,coord_transform)

    d%lon = x2%lon - x1%lon
    d%lat = x2%lat - x1%lat
    IF (d%lon >  180.0_wp) d%lon = d%lon - 360._wp
    IF (d%lon < -180.0_wp) d%lon = d%lon + 360._wp

    n = d%lon*d%lon + d%lat*d%lat
    p = (2*x1%lon*d%lon + 2*x1%lat*d%lat)/n
    q = (x1%lon*x1%lon + x1%lat*x1%lat - t*t)/n
    s = -1.0*p/2 + SQRT(p*p/4 - q)
    IF ((s<0._wp) .OR. (s>1._wp)) &
      & s = -1.0*p/2 - SQRT(p*p/4 - q)
    res%lon = x1%lon + s*d%lon
    res%lat = x1%lat + s*d%lat
    latitude_crossing = backtransform_lambert_azimuthal(res,coord_transform)
  END FUNCTION latitude_crossing


  !> Create index lists for edges with at least one vertex within
  !  "edge_thresh" to north/south pole plus the complementary set
  SUBROUTINE compute_index_lists(grid, thresh)
    TYPE (t_grid),        INTENT(INOUT) :: grid
    REAL(wp),             INTENT(IN)    :: thresh
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = TRIM('mo_remap_shared::compute_index_lists')
    INTEGER,          PARAMETER :: NLISTS  = 3
    INTEGER                           ::  ierrstat, start_idx, end_idx, &
      &                                   start_blk, end_blk, jb, jc,   &
      &                                   nblks_e, nblks_v, npromz_e,   &
      &                                   npromz_v,                     &
      &                                   new_jc, new_jb, ilist,        &
      &                                   npromz_c, nblks_c
    REAL(wp)                          ::  rlat, distn, dists

    ! allocate data structures:
    nblks_c  = grid%p_patch%nblks_c
    nblks_e  = grid%p_patch%nblks_e
    nblks_v  = grid%p_patch%nblks_v
    npromz_c = grid%p_patch%npromz_c
    npromz_e = grid%p_patch%npromz_e
    npromz_v = grid%p_patch%npromz_v

    grid%index_list%cnlist(1:NLISTS) = 0
    ALLOCATE(grid%index_list%clist_idx(nproma, nblks_c,NLISTS),      &
      &      grid%index_list%clist_blk(nproma, nblks_c,NLISTS),      &
      &      STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
    grid%index_list%l_initialized = .TRUE.

    ! first, divide cells into categories NORTH, SOUTH, DEFAULT:
    IF (dbg_level >= 5) WRITE (0,*) "# Create cell lists"
    start_blk = 1
    end_blk   = nblks_c
    DO jb=start_blk,end_blk
      start_idx = 1
      end_idx   = nproma
      if (jb == end_blk) end_idx = npromz_c
      DO jc=start_idx,end_idx
        ! cell center "distance" to the pole; not necessarily a norm but
        ! merely a function of the latitude value:
        rlat = grid%p_patch%cells%center(jc,jb)%lat
        distn = ABS(npole%lat - rlat)
        dists = ABS(spole%lat - rlat)
        IF (distn < thresh) THEN
          ilist = LIST_NPOLE
        ELSE IF (dists < thresh) THEN
          ilist = LIST_SPOLE
        ELSE
          ilist = LIST_DEFAULT
        END IF

        grid%index_list%cnlist(ilist) = grid%index_list%cnlist(ilist) + 1
        new_jc = idx_no(grid%index_list%cnlist(ilist))
        new_jb = blk_no(grid%index_list%cnlist(ilist))
        grid%index_list%clist_idx(new_jc, new_jb, ilist) = jc
        grid%index_list%clist_blk(new_jc, new_jb, ilist) = jb
      END DO ! jc
    END DO !jb
  END SUBROUTINE compute_index_lists


  !> Clean up data structures.
  !
  SUBROUTINE finalize_index_lists(grid)
    TYPE (t_grid),        INTENT(INOUT) :: grid
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = TRIM('mo_remap_shared::finalize_index_lists')
    INTEGER :: ierrstat

    ! clear index list (if allocated)
    IF (grid%index_list%l_initialized) THEN
      DEALLOCATE(grid%index_list%clist_idx, grid%index_list%clist_blk,     &
        &        STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")
      grid%index_list%l_initialized = .FALSE.
    END IF
  END SUBROUTINE finalize_index_lists


  !> Creates a copy of the vertex coordinates of a given grid which is
  !  then transformed using the Lambert equivalent azimuthal
  !  projection.
  !
  ! @note For simplicity we transform *all* vertices though we need
  !       only the poly region. The problem is that we do not know a priori,
  !       where this region exactly ends.
  !
  ! @todo Use at least a threshold to skip vertices which are clearly
  !       outside the pole region.
  !
  SUBROUTINE compute_coordinate_transform(grid, ilist)
    TYPE (t_grid), INTENT(INOUT) :: grid  !< grid data structure
    INTEGER,       INTENT(IN)    :: ilist
    ! local variables
    INTEGER                      ::  start_idx, end_idx, &
      &                              start_blk, end_blk, jb, jc
    TYPE (t_geographical_coordinates) :: p

    ! loop over vertices:
    IF (dbg_level >= 5) &
      &   WRITE (0,*) "# compute ", TRIM(LIST_NAME(ilist)), &
      &               " coordinate transform for ", TRIM(grid%name)
    start_blk = 1
    end_blk   = grid%p_patch%nblks_v
!$OMP PARALLEL DO &
!$OMP PRIVATE(jb,start_idx,end_idx,jc,p)
    DO jb=start_blk,end_blk
      start_idx = 1
      end_idx   = nproma
      if (jb == end_blk) end_idx = grid%p_patch%npromz_v
      DO jc=start_idx,end_idx
        p = grid%p_patch%verts%vertex(jc,jb)
        grid%vertex(jc,jb,ilist) = transform_lambert_azimuthal(p, ilist)

        ! *slightly* disturb those vertices which would otherwise land
        ! *exactly in the origin:
        IF ((ABS(grid%vertex(jc,jb,ilist)%lon) < 1e-10) .AND.  &
          & (ABS(grid%vertex(jc,jb,ilist)%lat) < 1e-10) .AND. &
          & (grid%p_patch%cell_type == 3)) THEN
          IF (ilist == LIST_NPOLE) p%lat = p%lat - 1e-2
          IF (ilist == LIST_SPOLE) p%lat = p%lat + 1e-2
          grid%vertex(jc,jb,ilist) = transform_lambert_azimuthal(p, ilist)
          grid%vertex(jc,jb,3)             = p
          grid%p_patch%verts%vertex(jc,jb) = p
        END IF
      END DO ! jc
    END DO !jb
!$OMP END PARALLEL DO

  END SUBROUTINE compute_coordinate_transform


  !> Copy the vertex coordinate list from "p_patch%vertes%vertex"
  !
  SUBROUTINE copy_vertex_coords(grid)
    TYPE (t_grid), INTENT(INOUT) :: grid
    CHARACTER(LEN=*), PARAMETER :: routine = TRIM('mo_remap_shared::copy_vertex_coords')
    INTEGER :: ierrstat

    IF (.NOT. ALLOCATED(grid%vertex)) THEN
      ALLOCATE(grid%vertex(nproma, grid%p_patch%nblks_v,3), &
        &      STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
    END IF
    ! copy data structure:
    grid%vertex(:,:,LIST_DEFAULT) = grid%p_patch%verts%vertex(:,:)
  END SUBROUTINE copy_vertex_coords


  !> Allocated temporary internal data structures for remapping
  !  algorithm.
  !
  SUBROUTINE allocate_lookup_tables(grid1, grid2, nthreads)
    TYPE (t_grid), INTENT(INOUT) :: grid1, grid2
    INTEGER,       INTENT(IN)    :: nthreads
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = TRIM('mo_remap_shared::allocate_lookup_tables')
    INTEGER :: ierrstat

  IF (.NOT. grid1%lookup_tbl%l_initialized) THEN
    ALLOCATE(grid1%lookup_tbl%vertex_c_idx(nproma, grid2%p_patch%nblks_v,3,nthreads), &
      &      grid1%lookup_tbl%vertex_c_blk(nproma, grid2%p_patch%nblks_v,3,nthreads), &
      &      STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
    grid1%lookup_tbl%l_initialized = .TRUE.
  END IF
  IF (.NOT. grid2%lookup_tbl%l_initialized) THEN
    ALLOCATE(grid2%lookup_tbl%vertex_c_idx(nproma, grid1%p_patch%nblks_v,3,nthreads), &
      &      grid2%lookup_tbl%vertex_c_blk(nproma, grid1%p_patch%nblks_v,3,nthreads), &
      &      STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
    grid2%lookup_tbl%l_initialized = .TRUE.
  END IF

  END SUBROUTINE allocate_lookup_tables


  !> Clear temporary internal data structures.
  !
  SUBROUTINE clear_lookup_tables(grid1, grid2)
    TYPE (t_grid), INTENT(INOUT) :: grid1, grid2

    IF (grid1%lookup_tbl%l_initialized) THEN
      grid1%lookup_tbl%vertex_c_idx(:,:,:,:) = -1
      grid1%lookup_tbl%vertex_c_blk(:,:,:,:) = -1
    END IF
    IF (grid2%lookup_tbl%l_initialized) THEN
      grid2%lookup_tbl%vertex_c_idx(:,:,:,:) = -1
      grid2%lookup_tbl%vertex_c_blk(:,:,:,:) = -1
    END IF
  END SUBROUTINE clear_lookup_tables


  !> Deallocate internal data structures.
  !
  SUBROUTINE deallocate_lookup_tables(grid1, grid2)
    TYPE (t_grid), INTENT(INOUT) :: grid1, grid2
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = TRIM('mo_remap_shared::deallocate_lookup_tables')
    INTEGER :: ierrstat

    IF (grid1%lookup_tbl%l_initialized) THEN
      DEALLOCATE(grid1%lookup_tbl%vertex_c_idx, grid1%lookup_tbl%vertex_c_blk, &
        &        STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")
      grid1%lookup_tbl%l_initialized = .FALSE.
    END IF
    IF (grid2%lookup_tbl%l_initialized) THEN
      DEALLOCATE(grid2%lookup_tbl%vertex_c_idx, grid2%lookup_tbl%vertex_c_blk, &
        &        STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")
      grid2%lookup_tbl%l_initialized = .FALSE.
    END IF
  END SUBROUTINE deallocate_lookup_tables


  !> clean up data structure.
  !
  SUBROUTINE finalize_vertex_coords(grid)
    TYPE (t_grid), INTENT(INOUT) :: grid
    CHARACTER(LEN=*), PARAMETER :: routine = TRIM('mo_remap_shared::finalize_vertex_coords')
    INTEGER :: ierrstat

    IF (ALLOCATED(grid%vertex)) THEN
      DEALLOCATE(grid%vertex, STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
    END IF
  END SUBROUTINE finalize_vertex_coords


  ! --------------------------------------------------------------------
  !> @return global index of containing cell.
  !          generic implementation for unstructured grids.
  !
  FUNCTION get_containing_cell_generic(grid, p_in, coord_transform) RESULT(global_idx)
    INTEGER :: global_idx 
    TYPE (t_grid),                     INTENT(IN) :: grid
    TYPE (t_geographical_coordinates), INTENT(IN) :: p_in
    INTEGER,                           INTENT(IN) :: coord_transform
    ! local variables
    INTEGER :: jc, jb, start_idx, start_blk, end_idx, end_blk, ne
    TYPE (t_geographical_coordinates) :: p,q
    TYPE (t_poly) :: poly

    IF (coord_transform == LIST_DEFAULT) THEN
      p = normalized_coord(p_in)
    ELSE
      p = p_in
    END IF

    ! we simply loop over all cells and check...
    ne = 0
    start_blk = 1
    end_blk   = grid%p_patch%nblks_c
    global_idx = -1
    DO jb=start_blk,end_blk
      start_idx = 1
      end_idx   = nproma
      if (jb == end_blk) end_idx = grid%p_patch%npromz_c

      LOOP : DO jc=start_idx,end_idx
        ! make sure that we do not mix up pole regions when dealing
        ! with transformed coordinates:
        IF (coord_transform /= LIST_DEFAULT) THEN
          IF ((grid%p_patch%cells%center(jc,jb)%lat * r_ns(coord_transform)) < 0._wp)  CYCLE LOOP
        END IF
        CALL cell2poly(grid, jc,jb, coord_transform, poly, ne)
        q = p
        IF ( (q%lon < -90._wp) .AND. (MAXVAL(poly%p(1:ne)%lon)) > 180._wp) THEN
          q%lon = q%lon + 360._wp
        END IF
        IF (inside(q, ne, poly)) THEN
          global_idx = idx_1d(jc,jb)
          IF (dbg_level >= 12)  THEN
            WRITE (0,*) "# found for ", p, ": ", poly%p(1:ne)
            WRITE (0,*) "#   ", jc, jb
          END IF
        END IF
      END DO LOOP
    END DO
  END FUNCTION get_containing_cell_generic


  !> Find containing cell using a GNAT data structure for fast search.
  !
  FUNCTION get_containing_cell_gnat(grid, p_in, coord_transform) RESULT(global_idx)
    INTEGER :: global_idx 
    TYPE (t_grid),                     INTENT(INOUT) :: grid
    TYPE (t_geographical_coordinates), INTENT(IN)    :: p_in
    INTEGER,                           INTENT(IN)    :: coord_transform
    ! local variables
    REAL (wp) :: r, min_dist, v(2)
    INTEGER   :: min_node_idx(3), count_dist, ne, icell, jb, jc
    TYPE (t_geographical_coordinates) :: p,q
    TYPE (t_poly) :: poly

    IF (coord_transform == LIST_DEFAULT) THEN
      p = normalized_coord(p_in)
      q = p
    ELSE
      p = p_in
      q = backtransform_lambert_azimuthal(p_in, coord_transform)
    END IF

    global_idx = -1
    count_dist =  0
    min_node_idx(1:2) = 0
    min_dist          = MAX_RANGE

    v = (/ q%lon * pi_180, q%lat * pi_180 /)  ! search point
    r = gnat_std_radius(gnat_tree)                  ! search radius
    CALL gnat_recursive_query(gnat_tree, v, r, min_dist, min_node_idx, count_dist)

    ! For local grids: abort if there is no "nearest cell"
    IF (min_node_idx(1) == 0) RETURN
    ! loop over vertex neighbors
    LOOP : DO icell=1,grid%vertex_nb_stencil(min_node_idx(1), min_node_idx(2))
      jc = grid%vertex_nb_idx(min_node_idx(1), min_node_idx(2), icell)
      jb = grid%vertex_nb_blk(min_node_idx(1), min_node_idx(2), icell)
      IF (.NOT. (jc > 0)) CYCLE
      ! make sure that we do not mix up pole regions when dealing
      ! with transformed coordinates:
      IF (coord_transform /= LIST_DEFAULT) THEN
        IF ((grid%p_patch%cells%center(jc,jb)%lat * r_ns(coord_transform)) < 0._wp)  CYCLE LOOP
      END IF
      ! check if cell found
      CALL cell2poly(grid, jc,jb, coord_transform, poly, ne)
      IF (coord_transform == LIST_DEFAULT) THEN
        q = p
        IF ( (q%lon < -90._wp) .AND. (MAXVAL(poly%p(1:ne)%lon)) > 180._wp) THEN
          q%lon = q%lon + 360._wp
        END IF
      ELSE
        q = p_in
      END IF
      IF (inside(q, ne, poly)) THEN
        global_idx = idx_1d(jc,jb)
        EXIT LOOP
      END IF
    END DO LOOP

    ! last resort: sequential search
    !    IF (global_idx == -1) THEN
    !      
    !      global_idx = get_containing_cell_generic(grid,p_in,coord_transform)
    !      IF (dbg_level >= 5)  WRITE (0,*) "# p=",p_in,": resort to ", global_idx
    !    END IF
  END FUNCTION get_containing_cell_gnat

  
  !> @return +1.0, if edge is oriented counter-clockwise wrt. a given cell,
  !          otherwise return -1.0
  !
  !  @todo Increase performance of this routine by generating a
  !        look-up table in advance!
  FUNCTION ccw_orientation(edge_idx, edge_blk, jc, jb, grid)
    REAL(wp)                  :: ccw_orientation
    INTEGER,       INTENT(IN) :: edge_idx, edge_blk, jc, jb
    TYPE (t_grid), INTENT(IN) :: grid    
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = TRIM('mo_remap_shared::ccw_orientation')
    INTEGER :: ne, i, local_idx, cedge_idx(4), cedge_blk(4), &
      &        cvertex_idx(4), cvertex_blk(4), evertex_idx1, evertex_blk1

    ne = grid%p_patch%cell_type
    ! loop over the cell's edges (which are assumed to be ordered ccw):
    local_idx = -1
    cedge_idx(1:ne) = grid%p_patch%cells%edge_idx(jc,jb,1:ne)
    cedge_blk(1:ne) = grid%p_patch%cells%edge_blk(jc,jb,1:ne)
    LOOP : DO i=1,ne
      IF ((cedge_idx(i) == edge_idx) .AND. (cedge_blk(i) == edge_blk)) THEN
        local_idx = i
        EXIT LOOP
      END IF
    END DO LOOP
    IF (local_idx == -1) CALL finish(routine, "Edge not found!")

    ! now compare the corresponding edge vertices with the cell's
    ! vertices (which are assumed to be ordered ccw):
    cvertex_idx(1:ne) = grid%p_patch%cells%vertex_idx(jc, jb, 1:ne)
    cvertex_blk(1:ne) = grid%p_patch%cells%vertex_blk(jc, jb, 1:ne)
    evertex_idx1 = grid%p_patch%edges%vertex_idx(edge_idx, edge_blk, 1)
    evertex_blk1 = grid%p_patch%edges%vertex_blk(edge_idx, edge_blk, 1)
    IF ((cvertex_idx(local_idx) == evertex_idx1) .AND. &
      & (cvertex_blk(local_idx) == evertex_blk1)) THEN
      ccw_orientation =  1.0_wp
    ELSE
      ccw_orientation = -1.0_wp
    END IF
  END FUNCTION ccw_orientation


  ! --------------------------------------------------------------------
  !> normalizes geographical coordinates (in degrees) to the 
  !  interval [-180., 180]/[-90.,90]
  !
  ELEMENTAL FUNCTION normalized_coord(p)
    TYPE (t_geographical_coordinates) :: normalized_coord
    TYPE (t_geographical_coordinates), INTENT(IN) :: p

    normalized_coord = p
    IF (p%lon > 180._wp) normalized_coord%lon = p%lon - 360._wp
    IF (normalized_coord%lon > 180._wp) normalized_coord%lon = normalized_coord%lon - 360._wp
    normalized_coord%lat = MIN(MAX(-90._wp, normalized_coord%lat), 90._wp)
  END FUNCTION normalized_coord

END MODULE mo_remap_shared
