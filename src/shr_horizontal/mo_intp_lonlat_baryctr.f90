#ifdef __xlC__
  @PROCESS smp=noopt
  @PROCESS noopt
#endif
#ifdef __PGI
  !pgi$g opt=1
#endif

!>
!! Contains the implementation of interpolation onto regular grids
!! with barycentric interpolation.
!!
!! @par Revision History
!! Moved from mo_intp_lonlat : 2015-05-29, F. Prill (DWD)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
  MODULE mo_intp_lonlat_baryctr
    !-------------------------------------------------------------------------
    !
    !    ProTeX FORTRAN source: Style 2
    !    modified for ICON project, DWD/MPI-M 2006
    !
    !-------------------------------------------------------------------------
    !
!$  USE OMP_LIB
    USE mo_kind,                ONLY: wp
    USE mo_exception,           ONLY: message, finish
    USE mo_impl_constants,      ONLY: SUCCESS, min_rlcell_int, min_rlcell
    USE mo_model_domain,        ONLY: t_patch
    USE mo_math_utilities,      ONLY: gc2cc, cc2gc, t_cartesian_coordinates,                &
      &                               t_geographical_coordinates, cc_dot_product
    USE mo_math_constants,      ONLY: pi, pi_2
    USE mo_parallel_config,     ONLY: nproma
    USE mo_loopindices,         ONLY: get_indices_c
    USE mo_intp_data_strc,      ONLY: t_lon_lat_intp
    USE mo_mpi,                 ONLY: get_my_mpi_work_id,                                   &
      &                               my_process_is_stdio, p_n_work
    USE mo_communication,       ONLY: idx_1d, blk_no, idx_no
    USE mo_delaunay_types,      ONLY: t_point_list, point_list, point, t_spherical_cap,     &
      &                               spherical_cap, OPERATOR(/), t_triangulation,          &
      &                               ccw_spherical, t_point, OPERATOR(+), triangulation,   &
      &                               t_sphcap_list
    USE mo_delaunay,            ONLY: point_cloud_diam, triangulate,write_triangulation_vtk,&
      &                               triangulate_mthreaded, create_thin_covering
    USE mo_util_string,         ONLY: int2string
    USE mo_octree,              ONLY: t_range_octree, octree_init, octree_count_point,      &
      &                               octree_query_point, octree_finalize

    IMPLICIT NONE

    !> level of output verbosity (for debugging purposes)
    INTEGER, PARAMETER  :: dbg_level = 0

    CHARACTER(LEN=*), PARAMETER :: modname = 'mo_intp_lonlat_baryctr'

    PRIVATE

    PUBLIC :: baryctr_interpol_lonlat
    PUBLIC :: setup_barycentric_intp_lonlat
    PUBLIC :: setup_barycentric_intp_lonlat_repartition
    PUBLIC :: compute_auxiliary_triangulation

    INTERFACE compute_auxiliary_triangulation
      MODULE PROCEDURE compute_triangulation_local_partition
      MODULE PROCEDURE compute_triangulation_repartition
    END INTERFACE

  CONTAINS

    !> Generate a sequence of points that are evenly distributed
    !> (approximately) on the sphere.
    !
    !  This spiral approximation formula originates from
    !
    !  E.A Rakhmanov, E.B Saff, Y.M Zhou: Electrons on the
    !  sphere. Series in Approximations and Decompositions, 1994
    !
    SUBROUTINE compute_point_distribution(n, points)
      INTEGER,             INTENT(IN)    :: n         !< number of points to distribute
      TYPE (t_point_list), INTENT(INOUT) :: points
      ! local variables
      TYPE(t_geographical_coordinates) :: points_gc(n)
      TYPE(t_cartesian_coordinates)    :: cc
      TYPE(t_point)                    :: pt
      INTEGER  :: k
      REAL(wp) :: h

      CALL points%initialize()
      DO k=1,n
        h = MIN(MAX(-1._wp + 2._wp*(k-1)/(N-1), -1._wp), 1._wp)
        points_gc(k)%lat = ACOS(h) 
        IF ((k==1) .OR. (k==n)) THEN
          points_gc(k)%lon = 0._wp
        ELSE
          points_gc(k)%lon = MOD( (points_gc(k-1)%lon + 3.6_wp/SQRT(n*(1._wp-h*h))) , 2._wp*pi ) 
        END IF
        points_gc(k)%lon = points_gc(k)%lon - pi
        points_gc(k)%lat = points_gc(k)%lat - pi_2
        cc = gc2cc(points_gc(k))
        pt = point(cc%x(1),cc%x(2),cc%x(3), -1, -1._wp)
        CALL points%push_back(pt)
      END DO
    END SUBROUTINE compute_point_distribution


    !-------------------------------------------------------------------------

    !> compute barycentric coordinates u(1...3) for the point pt
    !> inside the triangle pidx(1...3). The vertex positions are
    !> provided through the array @p vertex.
    !
    SUBROUTINE compute_barycentric_coords(pt, p1,p2,p3, u)
      TYPE(t_geographical_coordinates), INTENT(IN)  :: pt                  !< query point (longitude/latitude)  
      REAL(wp),                         INTENT(IN)  :: p1(3), p2(3), p3(3) !< triangle vertices
      REAL(wp),                         INTENT(OUT) :: u(3)                !< barycentric coordinates (dim: 3)
      ! local variables
      TYPE(t_cartesian_coordinates) :: a,b,c,p
      REAL(wp) :: B1(2), B2(2), r(2), det
      
      ! Cartesian coordinates of triangle vertices
      a%x(:) = p1
      b%x(:) = p2
      c%x(:) = p3
      ! query point
      p = gc2cc(pt)
      ! solve linear system for the barycentric coordinates
      !
      ! define matrix and right hand side
      a%x(1:3) = a%x(1:3) - c%x(1:3)
      b%x(1:3) = b%x(1:3) - c%x(1:3)
      p%x(1:3) = p%x(1:3) - c%x(1:3)
      B1(1:2) = (/ cc_dot_product(a,a) , cc_dot_product(a,b) /)   ! (a-c)^2     ,  (a-c)*(b-c)
      B2(1:2) = (/ B1(2)               , cc_dot_product(b,b) /)   ! (a-c)*(b-c) ,  (b-c)^2
      r(1:2)  = (/ cc_dot_product(p,a) , cc_dot_product(p,b) /)   ! [ (p-c)*(a-c), (p-c)*(b-c) ]
      ! solve using Cramer's rule:
      det     = B2(2)*B1(1) - B1(2)*B2(1) 
      u(1)    = ( r(1)*B2(2) - r(2)*B1(2) )/det
      u(2)    = ( r(2)*B1(1) - r(1)*B2(1) )/det
      u(3)    = 1._wp - (u(1) + u(2))
    END SUBROUTINE compute_barycentric_coords


    !-------------------------------------------------------------------------
    !> Simple test if a point is inside a triangle on the unit sphere.
    !
    FUNCTION inside_triangle(v, v1,v2,v3)
      LOGICAL :: inside_triangle
      REAL(wp),       INTENT(IN)     :: v(3)          !< query point
      TYPE(t_point),  INTENT(IN)     :: v1,v2,v3      !< vertex longitudes/latitudes
      ! local variables
      LOGICAL       :: c1,c2,c3
      TYPE(t_point) :: p

      p%x = v(1)
      p%y = v(2)
      p%z = v(3)

      c1  = ccw_spherical(v1,v2,p)
      c2  = ccw_spherical(v2,v3,p)
      c3  = ccw_spherical(v3,v1,p)
      inside_triangle = ((      c1) .AND. (      c2) .AND. (      c3)) .OR. &
        &               ((.NOT. c1) .AND. (.NOT. c2) .AND. (.NOT. c3))
    END FUNCTION inside_triangle


    !-------------------------------------------------------------------------
    !> Build a Delaunay triangulation connecting the cell circumcenters
    !  of "ptr_patch".
    !
    SUBROUTINE compute_triangulation_repartition (ptr_patch, tri_global, p_global, g2l_index)
      TYPE(t_patch), TARGET,  INTENT(IN)           :: ptr_patch          !< data structure containing grid info:
      TYPE (t_triangulation), INTENT(INOUT)        :: tri_global         !< resulting triangulation
      TYPE (t_point_list),    INTENT(INOUT)        :: p_global           !< resulting point set
      INTEGER, ALLOCATABLE,   INTENT(INOUT)        :: g2l_index(:)       !< point index mapping: global->local
      ! Local parameters:
      CHARACTER(*), PARAMETER :: routine = modname//"::compute_triangulation_repartition"
      INTEGER, ALLOCATABLE              :: permutation(:)
      TYPE(t_point)                     :: centroid
      TYPE (t_point_list), TARGET       :: p_local
      TYPE(t_cartesian_coordinates)     :: p_x
      TYPE (t_point_list)               :: pivot_points
      TYPE (t_sphcap_list)              :: subset_list
      INTEGER                           :: dim, errstat,                     &
        &                                  i, idx, jb, jc, nthreads,         &
        &                                  i_startblk, i_endblk,             &
        &                                  rl_start, rl_end, i_startidx,     &
        &                                  i_endidx, i_nchdom
      TYPE (t_spherical_cap)            :: subset
      TYPE(t_cartesian_coordinates)     :: cc          ! Cart. coordinates
      TYPE(t_geographical_coordinates)  :: gc          ! geo. coordinates
      REAL(wp)                          :: perturbation
!$  DOUBLE PRECISION                    :: time_s, toc

      ! --- create an array-like data structure containing the local
      ! --- mass points
      IF (dbg_level > 10) THEN
        WRITE (0,*) "# create an array-like data structure containing the local mass points"
      END IF
      CALL p_local%initialize()
      CALL p_local%reserve(ptr_patch%n_patch_cells)

      rl_start   = 1
      rl_end     = min_rlcell_int
      i_nchdom   = MAX(1,ptr_patch%n_childdom)
      i_startblk = ptr_patch%cells%start_blk(rl_start,1)
      i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)
      i = 0
      DO jb = i_startblk, i_endblk
        CALL get_indices_c(ptr_patch, jb,        &
          &                i_startblk, i_endblk, &
          &                i_startidx, i_endidx, &
          &                rl_start, rl_end)
        DO jc = i_startidx, i_endidx
          IF(.NOT. ptr_patch%cells%decomp_info%owner_mask(jc,jb)) CYCLE
          p_x = gc2cc(ptr_patch%cells%center(jc,jb))
          i = i + 1
          idx = idx_1d(jc,jb)
          CALL p_local%push_back(point(p_x%x(1),p_x%x(2),p_x%x(3), &
            &                          iindex=ptr_patch%cells%decomp_info%glb_index(idx)))
        END DO
      END DO
      IF (dbg_level > 10) THEN
        WRITE (0,*) "# min ICON grid coord: ", MINVAL(p_local%a(0:(p_local%nentries-1))%z)
        WRITE (0,*) "# max ICON grid coord: ", MAXVAL(p_local%a(0:(p_local%nentries-1))%z)
      END IF

      ! --- create a global copy of all points

      IF (dbg_level > 10)  WRITE (0,*) "# create a global copy of all points"

      p_global = point_list(p_local)
      CALL p_local%destructor()
      CALL p_global%sync()

      IF (dbg_level > 10) THEN
        WRITE (0,*) "# total no. of points to triangulate: ", p_global%nentries
      END IF

      CALL p_global%quicksort()

      ! slightly disturb symmetric coordinates; this should make the
      ! Delaunay triangulation unique, cf. [Lawson1984]
      !
      ! we disturb the triangulation by adding
      !  f(lon,lat) = 1.e-10 * cos(lon) * sin(lat)

      IF (dbg_level > 10)  WRITE (0,*) "# slightly disturb symmetric coordinates"

      dim = 0
!$OMP PARALLEL DO PRIVATE(perturbation,gc,cc)
      DO i=0,(p_global%nentries-1)
        cc%x(1:3)       = (/ p_global%a(i)%x, p_global%a(i)%y, p_global%a(i)%z /)
        gc              = cc2gc(cc)
        perturbation    = COS(gc%lon) * SIN(gc%lat) * 1.e-10_wp
        gc%lon          = gc%lon + perturbation
        cc              = gc2cc(gc)
        p_global%a(i)%x = cc%x(1)
        p_global%a(i)%y = cc%x(2)
        p_global%a(i)%z = cc%x(3)
      END DO
!$OMP END PARALLEL DO

      ALLOCATE(permutation(0:(p_global%nentries-1)), STAT=errstat)
      IF (errstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed')

      DO i=0,(p_global%nentries-1)
        IF (p_global%a(i)%gindex /= -1) THEN
          permutation(i) = g2l_index(p_global%a(i)%gindex)
        ELSE
          permutation(i) = -1
        END IF
        p_global%a(i)%gindex = i
      END DO
      DEALLOCATE(g2l_index, STAT=errstat)
      IF (errstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed')
      CALL MOVE_ALLOC(from=permutation, to=g2l_index)

      ! --- create an auxiliary triangulation, where the vertices are
      ! --- the mass points
      IF (dbg_level > 10) THEN
        WRITE (0,*) "# create an auxiliary triangulation, where the vertices are the mass points"
      END IF
      CALL tri_global%initialize()

      nthreads = 1
      !$  nthreads = omp_get_max_threads()
      IF (p_n_work < 3) THEN
        ! note: for regular grids the position of the spherical cap
        ! should be one of the poles (meridional convergence of grid
        ! lines)
        centroid = point(0._wp, 0._wp, 1._wp)
        subset = spherical_cap(centroid, -1._wp)
        CALL triangulate_mthreaded(p_global, tri_global, subset, nthreads, ignore_completeness=.FALSE.)
      ELSE
!$    time_s = omp_get_wtime()

        ! generate list of spherical cap centers (one for each MPI task)
        IF (dbg_level >= 10) THEN
          WRITE (0,*) "# triangulate pivot points"
        END IF
        CALL compute_point_distribution(p_n_work, pivot_points)
        pivot_points%a(get_my_mpi_work_id())%gindex = get_my_mpi_work_id()
        CALL subset_list%initialize()
        CALL create_thin_covering(pivot_points, subset_list, get_my_mpi_work_id())
        IF (subset_list%nentries > 1)  CALL finish(routine, "Internal error!")
        IF (dbg_level >= 10) THEN
          WRITE (0,*) "# triangulate subset with radius ", subset_list%a(0)%radius
        END IF
        CALL triangulate_mthreaded(p_global, tri_global, subset_list%a(0), nthreads, ignore_completeness=.FALSE.)

        IF (dbg_level >= 10) THEN
          WRITE (0,'(a,i0,a,a,i0)') "# done. triangulation: ", tri_global%nentries, " triangles.", &
            &                       "; pts = ", p_global%nentries
        END IF
!$    toc = omp_get_wtime() - time_s
!$    IF (dbg_level > 10) THEN
!$      WRITE (0,*) get_my_mpi_work_id()," :: elapsed time: ", toc
!$    END IF

!$    time_s = omp_get_wtime()
        CALL tri_global%sync()
!$    toc = omp_get_wtime() - time_s
!$    IF (dbg_level > 10) THEN
!$      WRITE (0,*) get_my_mpi_work_id()," :: triangulation sync, elapsed time: ", toc
!$    END IF
        ! clean up
        CALL pivot_points%destructor()
        CALL subset_list%destructor()
      END IF
      IF (dbg_level > 10) THEN
        WRITE (0,*) "# no. of cells in auxiliary triangulation: ", tri_global%nentries
      END IF

      ! --- plotting for debugging purposes:
      !
      IF (my_process_is_stdio() .AND. (dbg_level > 20)) THEN
        CALL write_triangulation_vtk("tri_global.vtk", p_global, tri_global)
      END IF

    END SUBROUTINE compute_triangulation_repartition


    !-------------------------------------------------------------------------
    !> Compute Delaunay triangulation of mass points.
    !
    ! @par Revision History
    !      Initial implementation  by  F. Prill, DWD (2015-04)
    !
    SUBROUTINE compute_triangulation_local_partition (ptr_patch, tri, p_global)
      ! data structure containing grid info:
      TYPE(t_patch), TARGET,  INTENT(IN)           :: ptr_patch
      TYPE (t_triangulation), INTENT(INOUT)        :: tri
      TYPE (t_point_list),    INTENT(INOUT)        :: p_global
      ! Local Parameters:
      CHARACTER(*), PARAMETER :: routine = modname//"::compute_triangulation_local_partition"
      ! enlarge the local triangulation area by this factor
      REAL(wp),     PARAMETER :: RADIUS_FACTOR = 1.1_wp

      INTEGER                         :: jb, jc, i_nchdom,                            &
        &                                i_startidx, i_endidx, i_startblk, i_endblk,  &
        &                                rl_start, rl_end, i, errstat,                &
        &                                idx, nthreads, dim, ierrstat
      TYPE (t_point_list)             :: p_local
      TYPE (t_spherical_cap)          :: subset
!$    DOUBLE PRECISION                :: time_s, toc
      TYPE(t_point)                   :: centroid
      TYPE(t_cartesian_coordinates)   :: p_x
      INTEGER, ALLOCATABLE            :: g2l_index(:)
      TYPE (t_triangulation)          :: tri_global
      INTEGER, ALLOCATABLE            :: permutation(:) ! point index permutation: sorted -> ICON ordering

      !-----------------------------------------------------------------------

      CALL message(routine, '')

      ! --- create an array-like data structure containing the local
      ! --- mass points

      CALL p_local%initialize()
      CALL p_local%reserve(ptr_patch%n_patch_cells)
      rl_start   = 1
      rl_end     = min_rlcell_int
      i_nchdom   = MAX(1,ptr_patch%n_childdom)
      i_startblk = ptr_patch%cells%start_blk(rl_start,1)
      i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)
      i = 0
      DO jb = i_startblk, i_endblk
        CALL get_indices_c(ptr_patch, jb,        &
          &                i_startblk, i_endblk, &
          &                i_startidx, i_endidx, &
          &                rl_start, rl_end)
        DO jc = i_startidx, i_endidx
          IF(.NOT. ptr_patch%cells%decomp_info%owner_mask(jc,jb)) CYCLE
          p_x = gc2cc(ptr_patch%cells%center(jc,jb))
          i = i + 1
          idx = idx_1d(jc,jb)
          CALL p_local%push_back(point(p_x%x(1),p_x%x(2),p_x%x(3), &
            &                          iindex=ptr_patch%cells%decomp_info%glb_index(idx)))
        END DO
      END DO

      IF (dbg_level > 1) THEN
        WRITE (0,*) "min ICON grid coord: ", MINVAL(p_local%a(0:(p_local%nentries-1))%z)
        WRITE (0,*) "max ICON grid coord: ", MAXVAL(p_local%a(0:(p_local%nentries-1))%z)
      END IF

      ! --- create a translation table global 1D index -> local index
      ALLOCATE(g2l_index(ptr_patch%n_patch_cells_g), STAT=errstat)
      IF (errstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed')
      g2l_index(:) = -1

      rl_start   = 2
      rl_end     = min_rlcell ! note that we include halo cells!
      i_nchdom   = MAX(1,ptr_patch%n_childdom)
      i_startblk = ptr_patch%cells%start_blk(rl_start,1)
      i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)
      i = 0
!$OMP PARALLEL DO PRIVATE(jb,jc,i_startidx,i_endidx,idx)
      DO jb = i_startblk, i_endblk
        CALL get_indices_c(ptr_patch, jb,        &
          &                i_startblk, i_endblk, &
          &                i_startidx, i_endidx, &
          &                rl_start, rl_end)
        DO jc = i_startidx, i_endidx
          idx = idx_1d(jc,jb)
          g2l_index(ptr_patch%cells%decomp_info%glb_index(idx)) = idx
        END DO
      END DO
!$OMP END PARALLEL DO

      ! --- create a global copy of all points
      
      p_global = point_list(p_local)
      CALL p_global%sync()

      ALLOCATE(permutation(0:(p_global%nentries-1)), STAT=errstat)
      IF (errstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed')

      CALL p_global%quicksort()

      ! slightly disturb symmetric coordinates; this should make the
      ! Delaunay triangulation unique, cf. [Lawson1984]
      dim = 0
      DO i=0,(p_global%nentries-1)
        SELECT CASE(dim)
        CASE (0)
          p_global%a(i)%x = p_global%a(i)%x + 1.e-10_wp
        CASE (1)
          p_global%a(i)%y = p_global%a(i)%y + 1.e-10_wp
        CASE (2)
          p_global%a(i)%z = p_global%a(i)%z + 1.e-10_wp
        END SELECT
        dim = MOD(dim+1,3)
      END DO

      DO i=0,(p_global%nentries-1)
        IF (p_global%a(i)%gindex /= -1) THEN
          permutation(i) = g2l_index(p_global%a(i)%gindex)
        ELSE
          permutation(i) = -1
        END IF
        p_global%a(i)%gindex = i
      END DO
      DEALLOCATE(g2l_index, STAT=errstat)
      IF (errstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed')

      ! --- create an auxiliary triangulation, where the vertices are
      ! --- the mass points

      ! create a spherical cap around the centroid of local mass
      ! points:
      centroid = point(0._wp, 0._wp, 0._wp)
      DO i=0,(p_local%nentries-1)
        centroid = centroid + p_local%a(i)
      END DO
      centroid = centroid/REAL(p_local%nentries,wp)
      subset = spherical_cap(centroid, MIN(RADIUS_FACTOR*point_cloud_diam(p_local, centroid), pi - 1.e-12_wp))      
      IF (dbg_level > 1) THEN
        WRITE (0,*) "spherical cap around ", p_local%a(0)%x, p_local%a(0)%y, p_local%a(0)%z, "; radius ", subset%radius
      END IF
      CALL p_local%destructor()

      CALL tri%initialize()
!$    time_s = omp_get_wtime()
      nthreads = 1
!$    nthreads = omp_get_max_threads()

      ! for local domains we do not force complete Delaunay
      ! triangulations, since these domains contain pathological
      ! triangles near the boundary which would lead to a
      ! time-consuming triangulation process.
      IF (nthreads > 1) THEN
        CALL triangulate_mthreaded(p_global, tri, subset, nthreads, &
          &                        ignore_completeness = (ptr_patch%id > 1))
      ELSE
        CALL triangulate(p_global, tri, subset, ignore_completeness = (ptr_patch%id > 1))
      END IF
!$    toc = omp_get_wtime() - time_s
      IF (dbg_level > 1) THEN
!$      WRITE (0,*) get_my_mpi_work_id()," :: elapsed time: ", toc, " (radius was ", subset%radius, ")"
        WRITE (0,*) "no. of cells in auxiliary triangulation: ", tri%nentries
      END IF

      DO i=0,(p_global%nentries-1)
        p_global%a(i)%gindex = permutation(i)
      END DO
      DEALLOCATE(permutation, STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')


      ! --- plotting for debugging purposes:
      !
      ! --- write a plot of the local triangulation
      IF (dbg_level > 2) THEN
        WRITE (0,'(a,i0,a)') "# formed ", tri%nentries, " triangles."
        IF (dbg_level > 10) THEN
          CALL write_triangulation_vtk("test"//TRIM(int2string(get_my_mpi_work_id()))//".vtk", p_global, tri)
        END IF
      END IF

      ! --- write a plot of the global triangulation
      CALL tri%quicksort() 
      tri_global=triangulation(tri)
      CALL tri_global%sync()
      IF (my_process_is_stdio() .AND. (dbg_level > 10)) THEN
        CALL write_triangulation_vtk("tri_global.vtk", p_global, tri_global)
      END IF

    END SUBROUTINE compute_triangulation_local_partition


    !-------------------------------------------------------------------------
    !> Build an octree out of triangles bounding boxes. Triangles
    !  outside a given range are ommitted.
    !
    SUBROUTINE compute_triangle_bboxes(p_global, tri_global, minrange, maxrange, octree)
      TYPE (t_point_list),    INTENT(IN)    :: p_global             !< set of global points
      TYPE (t_triangulation), INTENT(IN)    :: tri_global           !< global auxiliary triangulation
      REAL(wp),               INTENT(IN)    :: minrange(0:2)        !< corner of minimum coords (Cart. coords)
      REAL(wp),               INTENT(IN)    :: maxrange(0:2)        !< corner of maximum coords (Cart. coords)
      TYPE (t_range_octree),  INTENT(INOUT) :: octree               !< octree data structure

      ! local parameters
      CHARACTER(*), PARAMETER :: routine = modname//"::compute_triangle_bboxes"
      ! enlarge the triangle bounding boxes to prevent empty queries
      REAL(wp),     PARAMETER :: BBOX_MARGIN = 1.e-4_wp
      ! local variables
      INTEGER                         :: nlocal_triangles, l, i, j, k, errstat
      REAL(wp)                        :: pmin0(0:2), pmax0(0:2)
      REAL(wp), ALLOCATABLE           :: pmin(:,:), pmax(:,:)
      INTEGER,  ALLOCATABLE           :: glb_index_tri(:)
      REAL(wp)                        :: pp(0:2)
      REAL(wp)                        :: brange(2,3)          !< box range (min/max, dim=1,2,3)
      !$  DOUBLE PRECISION                :: time_s, toc

      ! --- count the no. of triangles that are not far-off
      IF (dbg_level > 10)  WRITE (0,*) "# count the no. of triangles that are not far-off"

      !$  time_s = omp_get_wtime()

      nlocal_triangles = 0 
      ! TODO: OpenMP parallelization
      DO l=0,(tri_global%nentries-1)
        pmin0(:) =  99._wp
        pmax0(:) = -99._wp
        DO j=0,2
          pp(0) = p_global%a(tri_global%a(l)%p(j))%x
          pp(1) = p_global%a(tri_global%a(l)%p(j))%y
          pp(2) = p_global%a(tri_global%a(l)%p(j))%z
          DO k=0,2
            pmin0(k) = MIN(pmin0(k), pp(k))
            pmax0(k) = MAX(pmax0(k), pp(k))
          END DO
        END DO
        ! [FP] enlarge the triangle bounding boxes to prevent empty queries
        pmin0(:) = pmin0(:) - BBOX_MARGIN
        pmax0(:) = pmax0(:) + BBOX_MARGIN

        IF (.NOT. ((pmax0(0) < minrange(0)) .OR. (pmin0(0) > maxrange(0)) .OR.  &
          &        (pmax0(1) < minrange(1)) .OR. (pmin0(1) > maxrange(1)) .OR.  &
          &        (pmax0(2) < minrange(2)) .OR. (pmin0(2) > maxrange(2)))) THEN
          nlocal_triangles = nlocal_triangles + 1
        END IF
      END DO
      IF (dbg_level > 10) THEN
        WRITE (0,*) "# ", nlocal_triangles, " triangles are local to this PE."
      END IF

      !$  toc = omp_get_wtime() - time_s
      !$  IF (dbg_level > 10) THEN
      !$    WRITE (0,*) get_my_mpi_work_id()," :: count the no. of triangles that are not far-off; elapsed time: ", toc
      !$  END IF

      ! --- build a list of triangle bounding boxes s.t. we can find
      ! --- the triangles containing our lon-lat points

      IF (dbg_level > 10)  WRITE (0,*) "# build a list of triangle bounding boxes"

      !$  time_s = omp_get_wtime()

      ALLOCATE(pmin(nlocal_triangles,3), pmax(nlocal_triangles,3), &
        &      glb_index_tri(nlocal_triangles), STAT=errstat)
      IF (errstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed')
      i = 1
      DO l=0,(tri_global%nentries-1)
        pmin0(:) =  99._wp
        pmax0(:) = -99._wp
        DO j=0,2
          pp(0) = p_global%a(tri_global%a(l)%p(j))%x
          pp(1) = p_global%a(tri_global%a(l)%p(j))%y
          pp(2) = p_global%a(tri_global%a(l)%p(j))%z
          DO k=0,2
            pmin0(k) = MIN(pmin0(k), pp(k))
            pmax0(k) = MAX(pmax0(k), pp(k))
          END DO
        END DO
        ! [FP] enlarge the triangle bounding boxes to prevent empty queries
        pmin0(:) = pmin0(:) - BBOX_MARGIN
        pmax0(:) = pmax0(:) + BBOX_MARGIN

        IF (.NOT. ((pmax0(0) < minrange(0)) .OR. (pmin0(0) > maxrange(0)) .OR.  &
          &        (pmax0(1) < minrange(1)) .OR. (pmin0(1) > maxrange(1)) .OR.  &
          &        (pmax0(2) < minrange(2)) .OR. (pmin0(2) > maxrange(2)))) THEN
          ! create a translation table "local -> global"
          glb_index_tri(i) = l
          pmin(i,1:3) = pmin0(0:2)
          pmax(i,1:3) = pmax0(0:2)
          i = i + 1
        END IF
      END DO

      ! --- insert local triangles into a tree-like data structure

      IF (dbg_level > 10)  WRITE (0,*) "# insert local triangles into a tree-like data structure"
      brange(1,:) = (/ -1._wp, -1._wp, -1._wp /)
      brange(2,:) = (/  1._wp,  1._wp,  1._wp /)
      CALL octree_init(octree, brange, pmin, pmax, opt_index=glb_index_tri)
      DEALLOCATE(pmin, pmax, glb_index_tri, STAT=errstat)
      IF (errstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed')

      !$  toc = omp_get_wtime() - time_s
      !$  IF (dbg_level > 10) THEN
      !$    WRITE (0,*) get_my_mpi_work_id()," :: build a list of triangle bounding boxes; elapsed time: ", toc
      !$  END IF
    END SUBROUTINE compute_triangle_bboxes

    
    !-------------------------------------------------------------------------
    !> Compute barycentric coordinates for a set of points, based on a
    !  given auxiliary triangulation.
    !
    SUBROUTINE compute_barycentric_coordinates(ptr_patch,p_global, tri_global, octree, &
      &                                        g2l_index, ptr_int_lonlat)
      TYPE(t_patch), TARGET,  INTENT(IN)    :: ptr_patch          !< data structure containing grid info:
      TYPE (t_point_list),    INTENT(IN)    :: p_global             !< set of global points
      TYPE (t_triangulation), INTENT(IN)    :: tri_global           !< global auxiliary triangulation
      TYPE (t_range_octree),  INTENT(IN)    :: octree               !< octree data structure
      INTEGER,                INTENT(IN)    :: g2l_index(:)         !< point index mapping: global->local
      ! Indices of source points and interpolation coefficients
      TYPE (t_lon_lat_intp), TARGET, INTENT(INOUT) :: ptr_int_lonlat
      ! local parameters
      CHARACTER(*), PARAMETER :: routine = modname//"::compute_barycentric_coordinates"
      ! we use the barycentric coords for the "point in triangle
      ! test"; this is the threshold for this test
      REAL(wp),     PARAMETER :: INSIDETEST_TOL = 1.e-6
      ! max. no. of triangles (bounding boxes) containing a single lat-lon point.
      INTEGER,      PARAMETER :: NMAX_HITS = 99
      ! local variables
      INTEGER                               :: jb, jc, start_idx, end_idx, nobjects, &
        &                                      i, j, k, idx0, idx1(3), nblks_lonlat, &
        &                                      npromz_lonlat
      !$  DOUBLE PRECISION                  :: time_s, toc
      INTEGER                               :: obj_list(NMAX_HITS)  !< query result (triangle search)
      TYPE(t_cartesian_coordinates)         :: ll_point_c           !< cartes. coordinates of lon-lat points
      REAL(wp)                              :: v(0:2,3)
      LOGICAL                               :: inside_test
      TYPE(t_cartesian_coordinates)         :: p_x
      INTEGER                               :: last_idx1(3)

      IF (dbg_level > 10)  WRITE (0,*) "# compute barycentric coordinates"

      ! make sure that the interpolation data structure for the
      ! barycentric interpolation has been allocated:
      IF (.NOT. ALLOCATED(ptr_int_lonlat%baryctr_coeff)) THEN
        CALL finish(routine, "Data structure for the barycentric interpolation not allocated!")
      END IF

      ! set local values for "nblks" and "npromz"
      nblks_lonlat  = blk_no(ptr_int_lonlat%nthis_local_pts)
      npromz_lonlat = ptr_int_lonlat%nthis_local_pts - (nblks_lonlat-1)*nproma

      !$  time_s = omp_get_wtime()

!$OMP PARALLEL DO PRIVATE(jb,jc,start_idx,end_idx,ll_point_c,nobjects,obj_list,      &
!$OMP                     idx0, idx1, v,i,j,k, inside_test, p_x, last_idx1)
      DO jb=1,nblks_lonlat
        start_idx = 1
        end_idx   = nproma
        IF (jb == nblks_lonlat) end_idx = npromz_lonlat

        DO jc=start_idx,end_idx

          ptr_int_lonlat%baryctr_idx(1:3,jc,jb)   = 1
          ptr_int_lonlat%baryctr_blk(1:3,jc,jb)   = 1
          ptr_int_lonlat%baryctr_coeff(1:3,jc,jb) = (/ 1._wp,  0._wp, 0._wp /)

          ! --- determine the triangle in our auxiliary triangulation
          ! --- which contains the lon-lat grid point:

          ! lon-lat point in Cartesian coordinates
          ll_point_c = gc2cc( ptr_int_lonlat%ll_coord(jc,jb) )
          ! query triangles whose bounding boxes contain this point:
          nobjects = octree_count_point(octree, ll_point_c%x(1:3))
          IF (nobjects > NMAX_HITS) THEN
            WRITE (0,*) "point ", ll_point_c%x, " hits ", nobjects, " objects."
            CALL finish(routine, "Internal error!")
          ELSE
            CALL octree_query_point(octree, ll_point_c%x(1:3), obj_list)
          END IF

          ! now test which of the triangles in "obj_list" actually
          ! contains "ll_point_c":
          idx0         = -1
          last_idx1(:) = -1
          LOOP: DO i=1,nobjects
            j = obj_list(i)

            DO k=0,2
              v(k,:) = (/ p_global%a(tri_global%a(j)%p(k))%x, &
                &         p_global%a(tri_global%a(j)%p(k))%y, &
                &         p_global%a(tri_global%a(j)%p(k))%z /)
            END DO

            ! --- compute the barycentric interpolation weights for
            ! --- this triangle

            CALL compute_barycentric_coords(ptr_int_lonlat%ll_coord(jc,jb),     &
              &                             v(0,:),v(1,:),v(2,:),               &
              &                             ptr_int_lonlat%baryctr_coeff(1:3,jc,jb))

            ! test if either the barycentric interpolation weights
            ! indicate that "ll_point_c" lies inside the triangle or
            ! if the test by dot-product succeeds:
            inside_test = ( ALL((ptr_int_lonlat%baryctr_coeff(1:3,jc,jb)) >= -1._wp*INSIDETEST_TOL)  .AND. &
              &             ALL(ptr_int_lonlat%baryctr_coeff(1:3,jc,jb)   <=  1._wp+INSIDETEST_TOL))
            IF (.NOT. inside_test) THEN
              inside_test = inside_triangle(ll_point_c%x, &
                &                           p_global%a(tri_global%a(j)%p(0)), &
                &                           p_global%a(tri_global%a(j)%p(1)), &
                &                           p_global%a(tri_global%a(j)%p(2)))
            END IF

            IF (inside_test) THEN
              idx0    = j
              ! get global index of cell circumcenters:
              idx1(:) = p_global%a(tri_global%a(idx0)%p(0:2))%gindex+1
              ! a query point may lie exactly on the edge between two
              ! triangles. We need to make sure that the result is
              ! processor-independent by choosing the triangles with
              ! larger indices.
              IF  (.NOT.  ((idx1(1) >  last_idx1(1))                                    .OR. &
                &         ((idx1(1) == last_idx1(1)) .AND. (idx1(2) >  last_idx1(2)))   .OR. &
                &         ((idx1(1) == last_idx1(1)) .AND. (idx1(2) == last_idx1(2)) &
                &                                    .AND. (idx1(3) >  last_idx1(3))))) THEN
                CYCLE LOOP
              END IF
              ! convert global to local indices
              idx1(:) = g2l_index(idx1(:))
              last_idx1(:) = idx1(:)

              IF (ANY(idx1(:) == -1)) THEN
                ! the containing triangle is not local for this PE?
                WRITE (0,*) "indices: ", idx1
                CALL finish(routine, "Internal error: The containing triangle is not local for this PE!")
              END IF
              
            END IF
          END DO LOOP
          
          ! no triangle was matching?
          IF (last_idx1(1) == -1) THEN
            CALL finish(routine, "Internal error: no triangle was matching!")
          END IF

          ! since the point set for the auxiliary triangulation
          ! had been distorted, we need to re-compute the
          ! barycentric coordinates based on the original
          ! triangle circumcenters:
          DO k=0,2
            p_x = gc2cc(ptr_patch%cells%center(idx_no(last_idx1(k+1)),blk_no(last_idx1(k+1))))
            v(k,:) = p_x%x(1:3)
          END DO
          CALL compute_barycentric_coords(ptr_int_lonlat%ll_coord(jc,jb),     &
            &                             v(0,:),v(1,:),v(2,:),               &
            &                             ptr_int_lonlat%baryctr_coeff(1:3,jc,jb))
            
          ! get indices of the containing triangle
          ptr_int_lonlat%baryctr_idx(1:3,jc,jb) = idx_no(last_idx1(1:3))
          ptr_int_lonlat%baryctr_blk(1:3,jc,jb) = blk_no(last_idx1(1:3))
                
          IF (ANY(ptr_int_lonlat%baryctr_idx(1:3,jc,jb) <= 0)) THEN
            WRITE (0,*) "g2l_index(tri_global%a(idx0)%p(0:2)) = ", last_idx1(:)
            CALL finish(routine, "Internal error!")
          END IF
          
        END DO
      END DO
!$OMP END PARALLEL DO

      !$  toc = omp_get_wtime() - time_s
      !$  IF (dbg_level > 10) THEN
      !$    WRITE (0,*) get_my_mpi_work_id()," :: compute barycentric coordinates; elapsed time: ", toc
      !$  END IF
    END SUBROUTINE compute_barycentric_coordinates


    !-------------------------------------------------------------------------
    !> Setup routine for barycentric interpolation at lon-lat grid
    !  points for an arbitrary grid.
    !
    ! @par Revision History
    !      Initial implementation  by  F. Prill, DWD (2015-03)
    !
    SUBROUTINE setup_barycentric_intp_lonlat_repartition(ptr_patch, ptr_int_lonlat)
      ! data structure containing grid info:
      TYPE(t_patch), TARGET,  INTENT(IN)           :: ptr_patch
      ! Indices of source points and interpolation coefficients
      TYPE (t_lon_lat_intp), TARGET, INTENT(INOUT) :: ptr_int_lonlat

      ! Local Parameters:
      CHARACTER(*), PARAMETER :: routine = modname//"::setup_barycentric_intp_lonlat_repartition"

      INTEGER                           :: jb, jc, start_idx, end_idx,       &
        &                                  i, errstat, idx,                  &
        &                                  i_startblk, i_endblk,             &
        &                                  rl_start, rl_end, i_startidx,     &
        &                                  i_endidx, i_nchdom, nblks_lonlat, &
        &                                  npromz_lonlat 
      TYPE (t_point_list)               :: p_global
      INTEGER, ALLOCATABLE              :: g2l_index(:)
      TYPE (t_range_octree)             :: octree               !< octree data structure
      TYPE(t_cartesian_coordinates)     :: ll_point_c           !< cartes. coordinates of lon-lat points
      TYPE (t_triangulation)            :: tri_global
      REAL(wp)                          :: minrange(0:2), maxrange(0:2)

      CALL message(routine, '')

      ! --- determine bounding box of query points

      ! set local values for "nblks" and "npromz"
      nblks_lonlat  = blk_no(ptr_int_lonlat%nthis_local_pts)
      npromz_lonlat = ptr_int_lonlat%nthis_local_pts - (nblks_lonlat-1)*nproma

      ! TODO: OpenMP parallelization
      minrange(:) =  99._wp
      maxrange(:) = -99._wp
      DO jb=1,nblks_lonlat
        start_idx = 1
        end_idx   = nproma
        IF (jb == nblks_lonlat) end_idx = npromz_lonlat

        DO jc=start_idx,end_idx
          ! lon-lat point in Cartesian coordinates
          ll_point_c = gc2cc( ptr_int_lonlat%ll_coord(jc,jb) )

          DO i=0,2
            minrange(i) = MIN(minrange(i), ll_point_c%x(i+1))
            maxrange(i) = MAX(maxrange(i), ll_point_c%x(i+1))
          END DO
        END DO
      END DO

      ! --- create a translation table global 1D index -> local index

      IF (dbg_level > 10)  WRITE (0,*) "# create a translation table global 1D index -> local index"
      ALLOCATE(g2l_index(ptr_patch%n_patch_cells_g), STAT=errstat)
      IF (errstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed')
      g2l_index(:) = -1

      rl_start   = 2
      rl_end     = min_rlcell ! note that we include halo cells
      i_nchdom   = MAX(1,ptr_patch%n_childdom)
      i_startblk = ptr_patch%cells%start_blk(rl_start,1)
      i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)
      i = 0
      DO jb = i_startblk, i_endblk
        CALL get_indices_c(ptr_patch, jb,        &
          &                i_startblk, i_endblk, &
          &                i_startidx, i_endidx, &
          &                rl_start, rl_end)
        DO jc = i_startidx, i_endidx
          i = i + 1
          idx = idx_1d(jc,jb)
          g2l_index(ptr_patch%cells%decomp_info%glb_index(idx)) = idx
        END DO
      END DO

      ! --- compute a triangulation of the mass points; this modifies
      !     the "g2l_index" array since the point ordering is altered.
      CALL compute_auxiliary_triangulation(ptr_patch, tri_global, p_global, g2l_index)

      ! --- determine bounding box of triangles
      CALL compute_triangle_bboxes(p_global, tri_global, minrange, maxrange, octree)

      ! --- compute barycentric coordinates
      CALL compute_barycentric_coordinates(ptr_patch,p_global, tri_global, octree,     &
        &                                  g2l_index, ptr_int_lonlat)

      ! clean up
      CALL p_global%destructor()
      CALL tri_global%destructor()
      CALL octree_finalize(octree)
      DEALLOCATE(g2l_index, STAT=errstat)
      IF (errstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed')

    END SUBROUTINE setup_barycentric_intp_lonlat_repartition


    !-------------------------------------------------------------------------
    !> Setup routine for barycentric interpolation at lon-lat grid
    !  points for an arbitrary grid.
    !
    ! @par Revision History
    !      Initial implementation  by  F. Prill, DWD (2015-01)
    !
    SUBROUTINE setup_barycentric_intp_lonlat(tri, p_global, ptr_int_lonlat)
      ! triangulation of mass points.
      TYPE (t_triangulation),        INTENT(IN)    :: tri
      TYPE (t_point_list),           INTENT(IN)    :: p_global
      ! Indices of source points and interpolation coefficients
      TYPE (t_lon_lat_intp), TARGET, INTENT(INOUT) :: ptr_int_lonlat
      ! Local Parameters:
      CHARACTER(*), PARAMETER :: routine = modname//"::setup_barycentric_intp_lonlat"
      ! max. no. of triangle (bounding boxes) containing a single lat-lon point.
      INTEGER,      PARAMETER :: NMAX_HITS = 99
      ! enlarge the triangle bounding boxes to prevent empty queries
      REAL(wp),     PARAMETER :: BBOX_MARGIN = 1.e-3_wp
      ! we use the barycentric coords for the "point in triangle
      ! test"; this is the threshold for this test
      REAL(wp),     PARAMETER :: INSIDETEST_TOL = 1.e-6

      INTEGER                         :: nblks_lonlat, npromz_lonlat, jb, jc,    &
        &                                i_startidx, i_endidx, i, j, k, errstat, &
        &                                nobjects, idx0
      REAL(wp)                        :: pp(3),v1(3),v2(3),v3(3)
      TYPE (t_range_octree)           :: octree               !< octree data structure
      REAL(wp)                        :: brange(2,3)          !< box range (min/max, dim=1,2,3)
      REAL(wp), ALLOCATABLE           :: pmin(:,:), pmax(:,:)
      INTEGER                         :: obj_list(NMAX_HITS)  !< query result (triangle search)
      TYPE(t_cartesian_coordinates)   :: ll_point_c           !< cartes. coordinates of lon-lat points
      LOGICAL                         :: inside_test1, inside_test2

      !-----------------------------------------------------------------------

      CALL message(routine, '')

      ! make sure that the interpolation data structure for the
      ! barycentric interpolation has been allocated:
      IF (.NOT. ALLOCATED(ptr_int_lonlat%baryctr_coeff)) THEN
        CALL finish(routine, "Data structure for the barycentric interpolation not allocated!")
      END IF

      ALLOCATE(pmin(tri%nentries,3), pmax(tri%nentries,3), STAT=errstat)
      IF (errstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed')

      ! --- build a list of triangle bounding boxes s.t. we can find
      ! --- the triangles containing our lon-lat points

      DO i=1,tri%nentries
        pmin(i,:) =  99._wp
        pmax(i,:) = -99._wp
        DO j=0,2
          pp(1) = p_global%a(tri%a(i-1)%p(j))%x
          pp(2) = p_global%a(tri%a(i-1)%p(j))%y
          pp(3) = p_global%a(tri%a(i-1)%p(j))%z
          DO k=1,3
            pmin(i,k) = MIN(pmin(i,k), pp(k))
            pmax(i,k) = MAX(pmax(i,k), pp(k))
          END DO
        END DO

        ! [FP] enlarge the triangle bounding boxes to prevent empty queries
        pmin(i,:) = pmin(i,:) - BBOX_MARGIN
        pmax(i,:) = pmax(i,:) + BBOX_MARGIN
      END DO

      ! --- insert local triangles into a tree-like data structure

      brange(1,:) = (/ -1._wp, -1._wp, -1._wp /)
      brange(2,:) = (/  1._wp,  1._wp,  1._wp /)
      CALL octree_init(octree, brange, pmin, pmax)

      ! --- compute barycentric coordinates

      ! set local values for "nblks" and "npromz"
      nblks_lonlat  = blk_no(ptr_int_lonlat%nthis_local_pts)
      npromz_lonlat = ptr_int_lonlat%nthis_local_pts - (nblks_lonlat-1)*nproma

!$OMP PARALLEL DO PRIVATE(jb,jc,i_startidx,i_endidx,ll_point_c,nobjects,obj_list, &
!$OMP                     idx0, v1,v2,v3,i,j,inside_test1,inside_test2 )
      DO jb=1,nblks_lonlat
        i_startidx = 1
        i_endidx   = nproma
        IF (jb == nblks_lonlat) i_endidx = npromz_lonlat

        DO jc=i_startidx,i_endidx

          ptr_int_lonlat%baryctr_idx(1:3,jc,jb)   = 1
          ptr_int_lonlat%baryctr_blk(1:3,jc,jb)   = 1
          ptr_int_lonlat%baryctr_coeff(1:3,jc,jb) = (/ 1._wp,  0._wp, 0._wp /)

          ! --- determine the triangle in our auxiliary triangulation
          ! --- which contains the lon-lat grid point:

          ! lon-lat point in Cartesian coordinates
          ll_point_c = gc2cc( ptr_int_lonlat%ll_coord(jc,jb) )
          ! query triangles whose bounding boxes contain this point:
          nobjects = octree_count_point(octree, ll_point_c%x(1:3))
          IF (nobjects > NMAX_HITS) THEN
            WRITE (0,*) "point ", ll_point_c%x, " hits ", nobjects, " objects."
            CALL finish(routine, "Internal error!")
          ELSE
            CALL octree_query_point(octree, ll_point_c%x(1:3), obj_list)
          END IF

          ! now test which of the triangles in "obj_list" actually
          ! contains "ll_point_c":
          idx0 = -1
          LOOP: DO i=1,nobjects
            j = obj_list(i) - 1
            v1(:) = (/ p_global%a(tri%a(j)%p(0))%x, p_global%a(tri%a(j)%p(0))%y, p_global%a(tri%a(j)%p(0))%z /)
            v2(:) = (/ p_global%a(tri%a(j)%p(1))%x, p_global%a(tri%a(j)%p(1))%y, p_global%a(tri%a(j)%p(1))%z /)
            v3(:) = (/ p_global%a(tri%a(j)%p(2))%x, p_global%a(tri%a(j)%p(2))%y, p_global%a(tri%a(j)%p(2))%z /)

            ! --- compute the barycentric interpolation weights for
            ! --- this triangle
            CALL compute_barycentric_coords(ptr_int_lonlat%ll_coord(jc,jb),       &
              &                             v1,v2,v3,                             &
              &                             ptr_int_lonlat%baryctr_coeff(1:3,jc,jb))

            ! test if either the barycentric interpolation weights
            ! indicate that "ll_point_c" lies inside the triangle or
            ! if the test by dot-product succeeds:
            inside_test1 = ( ALL((ptr_int_lonlat%baryctr_coeff(1:3,jc,jb)) >= -1._wp*INSIDETEST_TOL)  .AND. &
              &              ALL(ptr_int_lonlat%baryctr_coeff(1:3,jc,jb)   <=  1._wp+INSIDETEST_TOL))
            inside_test2 = inside_triangle(ll_point_c%x, p_global%a(tri%a(j)%p(0)), p_global%a(tri%a(j)%p(1)), &
              &                            p_global%a(tri%a(j)%p(2)))

            IF (inside_test1 .OR. inside_test2) THEN
              idx0 = j

              IF (ALL(p_global%a(tri%a(idx0)%p(0:2))%gindex /= -1)) THEN
                ! get indices of the containing triangle
                ptr_int_lonlat%baryctr_idx(1:3,jc,jb) = idx_no(p_global%a(tri%a(idx0)%p(0:2))%gindex)
                ptr_int_lonlat%baryctr_blk(1:3,jc,jb) = blk_no(p_global%a(tri%a(idx0)%p(0:2))%gindex)

                IF (ANY(ptr_int_lonlat%baryctr_idx(1:3,jc,jb) <= 0)) THEN
                  WRITE (0,*) "permutation(tri%a(idx0)%p(0:2)) = ", p_global%a(tri%a(idx0)%p(0:2))%gindex
                  CALL finish(routine, "Internal error!")
                END IF
              ELSE
                ! the containing triangle is not local for this PE;
                ! this may happen for nested regions; we therefore do not stop.
              END IF

              EXIT LOOP
            END IF
          END DO LOOP

        END DO
      END DO
!$OMP END PARALLEL DO

      ! clean up
      CALL octree_finalize(octree)
      DEALLOCATE(pmin, pmax, STAT=errstat)
      IF (errstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed')

    END SUBROUTINE setup_barycentric_intp_lonlat


    !-------------------------------------------------------------------------
    !> Performs barycentric interpolation
    !
    ! @par Revision History
    !      Initial implementation  by  F. Prill, DWD (2015-01)
    !
    SUBROUTINE baryctr_interpol_lonlat( p_cell_in, ptr_int,                   &
      &                                 p_out, nblks_lonlat, npromz_lonlat,   &
      &                                 opt_slev, opt_elev)
      ! INPUT PARAMETERS
      !
      ! input cell-based variable for which gradient at cell center is computed
      REAL(wp),                      INTENT(IN)    :: p_cell_in(:,:,:) ! dim: (nproma,nlev,nblks_c)
      ! Indices of source points and interpolation coefficients
      TYPE (t_lon_lat_intp), TARGET, INTENT(IN)    :: ptr_int
      ! reconstructed scalar value at lon-lat point
      REAL(wp),                      INTENT(INOUT) :: p_out(:,:,:)     ! dim: (nproma,nlev,nblks_lonlat)
      ! lon-lat grid blocking info
      INTEGER,                       INTENT(IN)    :: nblks_lonlat, npromz_lonlat
      ! optional vertical start/end level
      INTEGER,                       INTENT(IN), OPTIONAL :: opt_slev, opt_elev

      ! LOCAL VARIABLES
      INTEGER :: slev, elev,                 & ! vertical start and end level
        &        i_startidx, i_endidx,       & ! start/end index
        &        jc, jb, jk                    ! integer over lon-lat points, levels

      slev = 1
      elev = UBOUND(p_cell_in,2)
      ! check optional arguments
      IF ( PRESENT(opt_slev) ) slev = opt_slev
      IF ( PRESENT(opt_elev) ) elev = opt_elev

      ! initialize output array with zeros
      p_out(:,slev:elev,:) = 0._wp

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,jc), SCHEDULE(runtime)
      DO jb = 1,nblks_lonlat
        i_startidx = 1
        i_endidx   = nproma
        IF (jb == nblks_lonlat) i_endidx = npromz_lonlat
#ifdef __LOOP_EXCHANGE
          DO jc = i_startidx, i_endidx
            DO jk = slev, elev
#else
!CDIR UNROLL=3
          DO jk = slev, elev
            DO jc = i_startidx, i_endidx
#endif

              p_out(jc,jk,jb) = &
                &    ptr_int%baryctr_coeff(1,jc,jb)*                                            &
                &    p_cell_in(ptr_int%baryctr_idx(1,jc,jb), jk, ptr_int%baryctr_blk(1,jc,jb))  &
                &  + ptr_int%baryctr_coeff(2,jc,jb)*                                            &
                &    p_cell_in(ptr_int%baryctr_idx(2,jc,jb), jk, ptr_int%baryctr_blk(2,jc,jb))  &
                &  + ptr_int%baryctr_coeff(3,jc,jb)*                                            &
                &    p_cell_in(ptr_int%baryctr_idx(3,jc,jb), jk, ptr_int%baryctr_blk(3,jc,jb))

            ENDDO
          ENDDO
      END DO
!$OMP END DO
!$OMP END PARALLEL
    END SUBROUTINE baryctr_interpol_lonlat
      
  END MODULE mo_intp_lonlat_baryctr
