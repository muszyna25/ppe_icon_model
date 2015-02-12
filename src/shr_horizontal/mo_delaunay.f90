!> Parallel Bowyer/Watson algorithm for the construction of a Delaunay
!! triangulation on the sphere.
!!
!! @author F. Prill, DWD
!!
!! This algorithm borrows from
!!
!! (1) Paul Bourke: "Efficient Triangulation Algorithm Suitable for
!!      Terrain Modelling or An Algorithm for Interpolating
!!      Irregularly-Spaced Data with Applications in Terrain
!!      Modelling." Pan Pacific Computer Conference, Beijing, China,
!!      1989.
!!
!! (2) Jacobsen, D. W.; Gunzburger, M.; Ringler, T.; Burkardt, J.;
!!      Peterson, J.: "Parallel algorithms for planar and spherical
!!      Delaunay construction with an application to centroidal
!!      Voronoi tessellations."  Geoscientific Model Development;2013,
!!      Vol. 6 Issue 4, p.1353 (global Delaunay condition for
!!      parallelization)
!!
!! (3) Renka, R. J. "Interpolation of Data on the Surface of a
!!      Sphere." ACM Trans. Math. Softw., ACM, 1984, 10, 417-436
!!
!! (4) J. Shewchuck: "Lecture Notes on Delaunay Triangulation."
!!      Univ. of California, Berkeley. 2012. (ghost point instead of
!!      "super-triangle", as described in Section 3.4 of)
!!
!! with modifications: 
!! - sorting of points wrt. spherical cap (point subset for parallelization)
!! - multi-threading and indexing shortcuts
!! - global Delaunay condition is used as abort criterion for triangulation
!!
!!
!! @par Revision History
!! Initial implementation,            F. Prill, DWD (2014-11-12)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.


MODULE mo_delaunay

  !$  USE OMP_LIB
  USE mo_exception,         ONLY: finish
  USE mo_impl_constants,    ONLY: SUCCESS
  USE mo_kind,              ONLY: wp
  USE mo_delaunay_types, ONLY: t_edge, t_point, t_triangle, t_point_list, t_triangulation,     &
    &                          t_spherical_cap, t_sphcap_list,                                 &
    &                          point_list, triangle, point, spherical_cap,                     &
    &                          OPERATOR(<), OPERATOR(/), OPERATOR(*), OPERATOR(==),            &
    &                          ccw_spherical, circum_circle_spherical
  IMPLICIT NONE
  
  PRIVATE

  PUBLIC :: triangulate
  PUBLIC :: triangulate_mthreaded
  PUBLIC :: triangulate_parallel
  PUBLIC :: write_triangulation_vtk
  PUBLIC :: point_cloud_diam

  CHARACTER(LEN=*), PARAMETER :: modname   = 'mo_delaunay'
  INTEGER,          PARAMETER :: dbg_level = 4

  ! guarantee the exact evaluation of determinants in the
  ! algorithm by rounding coordinates to 2^-17 (losing, however,
  ! accuracy), cf. [Lawson1984]
  INTEGER,      PARAMETER :: ACCURACY = wp

CONTAINS

  ! --------------------------------------------------------------------
  !> Triangulation subroutine
  SUBROUTINE triangulate(points, tri, subset)

    TYPE(t_point_list),    INTENT(IN)            :: points
    TYPE(t_triangulation), INTENT(INOUT), TARGET :: tri
    TYPE(t_spherical_cap), INTENT(IN)            :: subset
    ! local variables

    !> triangle count that triggers clean-up of "tri" data structure
    INTEGER, PARAMETER :: MAX_EDGES      = 1000   !< size of edge buffer

    ! For the triangles that are created from the edges 0,1,2: compute
    ! the edge opposite to "ghost point" (if "ghost point" part of this
    ! triangle tri[j]). This lookup table is best figured out on a
    ! sheet of paper...
    INTEGER, PARAMETER :: next_oedge(0:2,-1:2) = &
      &  RESHAPE( (/ -1, -1, -1, &
      &              -1,  2,  1, &
      &               1, -1,  2, &
      &               2,  1, -1 /), (/3,4/) )
    TYPE(t_point_list)        :: pxyz
    TYPE(t_point)             :: ipoint
    TYPE(t_edge)              :: bedge
    TYPE(t_triangle), POINTER :: jtri
    INTEGER                   :: i,j,k, npts, gp, j0, ndiscard, ipt, ne, ntri, &
      &                          jmin, oedge, jmin0, this_edge, endtri
    LOGICAL                   :: inside
    REAL(wp)                  :: cos_radius, rlimit

    ! edges buffer for re-Delaunay step (together with
    ! "edge-index-opposite-to-ghost-point")
    TYPE(t_edge)              :: edges(0:(MAX_EDGES-1))
    INTEGER                   :: edges_oedge(0:(MAX_EDGES-1))
    LOGICAL                   :: edges_valid(0:(MAX_EDGES-1))
                              
    ! copy all points since they will be re-ordered:
    npts = points%nentries
    CALL pxyz%initialize()
    pxyz = point_list(points)
    DO i=0,(npts-1)
      pxyz%a(i)        = pxyz%a(i)/pxyz%a(i)%norm2()
      ! --- guarantee the exact evaluation of determinants in the
      ! --- algorithm (losing, however, accuracy):
      pxyz%a(i)%x      = REAL(REAL(pxyz%a(i)%x, ACCURACY), wp)
      pxyz%a(i)%y      = REAL(REAL(pxyz%a(i)%y, ACCURACY), wp)
      pxyz%a(i)%z      = REAL(REAL(pxyz%a(i)%z, ACCURACY), wp)
      pxyz%a(i)%ps     = pxyz%a(i) * subset%sorting_direction
    END DO
    CALL pxyz%quicksort()

    ! assume triangulation covers entire sphere (data structure grows
    ! dynamically if not)
    CALL tri%reserve(2*npts-4)

    ! set up the initial triangle: insert the first three points
    IF (ccw_spherical(pxyz%a(0), pxyz%a(1), pxyz%a(2)) > 0._wp) THEN
      CALL tri%push_back(triangle(0, 2, 1)) ! counter-clockwise starting triangle
    ELSE 
      CALL tri%push_back(triangle(0, 1, 2)) ! clockwise starting triangle
    END IF

    ! connect triangle edges to a "ghost point"
    gp = npts
    CALL pxyz%push_back(point(-1._wp,-1._wp,-1._wp))
    ! add three initial "ghost triangles"
    CALL tri%push_back(triangle(tri%a(0)%p(1),tri%a(0)%p(0),gp, 0))
    CALL tri%push_back(triangle(tri%a(0)%p(2),tri%a(0)%p(1),gp, 0))
    CALL tri%push_back(triangle(tri%a(0)%p(0),tri%a(0)%p(2),gp, 0))
    CALL tri%a(0)%compute_circumcenter(pxyz, subset)

    ipt      = 3           ! current point index    
    j0       = 0           ! first triangle index that is not "complete"
    ndiscard = 0           ! no. of triangles that violate global Delaunay condition
    ipoint   = pxyz%a(ipt) ! current point
    IF (subset%radius < 0._wp) THEN
      cos_radius = 99._wp
    ELSE
      cos_radius = -1._wp*COS(subset%radius)
    END IF

    ! Include each point one at a time into the existing mesh. Abort if
    ! all points have been processed or the requested radius has been
    ! reached without violating the global Delaunay condition.
    LOOP : DO ipt=3,(npts-1)
      IF (.NOT. ((ipoint%ps <= cos_radius) .OR. (ndiscard>0)))  EXIT LOOP

      ipoint       = pxyz%a(ipt)
      ntri         = tri%nentries
      endtri       = ntri - 1
      ne           =  0 ! no. of entries in edge buffer
      jmin         = -1

      j = j0
      ! find the next "complete" triangle from the back
      DO
        IF ((tri%a(endtri)%complete==1) .OR. (endtri==j))  EXIT
        endtri = endtri - 1
      END DO 

      ! Set up the edge buffer.
      !
      ! If the point pxyz[i] lies inside the circumcircle then the
      ! three edges of that triangle are added to the edge buffer and
      ! that triangle is removed.
      LOOP2 : DO 
        IF (j>=ntri)  EXIT LOOP2

        ! test if triangle remains on "discard" list
        jtri => tri%a(j)
        IF (ipoint%ps > jtri%rdiscard) THEN
          ndiscard      = ndiscard-1
          jtri%rdiscard = 1._wp
        END IF

        IF (jtri%complete == 0) THEN
          oedge = jtri%oedge
          jmin0 = jmin
          IF (jmin == -1)  jmin = j
          IF (oedge /= -1) THEN
            ! boundary edge
            bedge  = jtri%edge(oedge)
            inside = (ccw_spherical(pxyz%a(bedge%p1), pxyz%a(bedge%p2), ipoint) <= 0._wp)
          ELSE
            ! interior edge
            inside = circum_circle_spherical(ipoint, pxyz, jtri%p)
            IF (.NOT. inside .AND. ((ipoint%ps - jtri%cc%ps) > jtri%r)) THEN
              jtri%complete =     1
              jmin          = jmin0
            END IF
          END IF
          
          IF (inside) THEN
            ! push triangle's edges onto edge list
            this_edge = ne
            ne = ne + 3
            DO k=0,2
              edges(this_edge)       = jtri%edge(k)
              edges_oedge(this_edge) = next_oedge(k,oedge)
              edges_valid(this_edge) = .TRUE.
              this_edge = this_edge + 1
            END DO

            IF (jtri%rdiscard < 1._wp)  ndiscard=ndiscard-1            
            ntri = ntri - 1
            IF (endtri > j) THEN
              tri%a(j)      = tri%a(endtri)
              tri%a(endtri) = tri%a(ntri)
              ! find the next "complete" triangle from the back
              LOOP5 : DO
                IF ((tri%a(endtri)%complete==1) .OR. (endtri==j))  EXIT LOOP5
                endtri = endtri - 1
              END DO LOOP5
            ELSE
              tri%a(j) = tri%a(ntri)
            END IF
            CYCLE LOOP2
          END IF
        END IF
        j = j + 1
      END DO LOOP2
      CALL tri%resize(ntri)

      ! "j0" denotes the smallest triangle index which is not yet
      ! "complete" (significant speedup!)
      IF (jmin >= j0)  j0 = jmin
      
      ! remove multiple edges
      DO j=0,(ne-1)
        IF (edges_valid(j)) THEN
          DO k=j+1,(ne-1)
            IF (edges(j) == edges(k)) THEN
              edges_valid(j) = .FALSE.
              edges_valid(k) = .FALSE.
            END IF
          END DO
        END IF
      END DO

      ! form new triangles for the current point. All edges are
      ! arranged in clockwise order.
      DO j=0,(ne-1)
        IF (edges_valid(j))  CALL tri%push_back( triangle(edges(j)%p1, edges(j)%p2, ipt, edges_oedge(j)) )
      END DO
      DO j=ntri,(tri%nentries-1)
          IF (tri%a(j)%oedge == -1) THEN  
            ! triangle without ghost point
            CALL tri%a(j)%compute_circumcenter(pxyz, subset)
            ! mark triangle if it violates the global Delaunay condition:
            IF ((cos_radius > ipoint%ps) .AND. (tri%a(j)%cap_distance < subset%radius)) THEN 
              ndiscard          = ndiscard + 1
              tri%a(j)%rdiscard = -1._wp*COS(tri%a(j)%cap_distance+2._wp*tri%a(j)%r)
            END IF
          END IF
      END DO
    END DO LOOP

    ! compute limit radius (last point inserted into triangulation)
    IF (subset%radius < 0._wp) THEN
      rlimit = 99._wp
    ELSE
      rlimit = subset%point%spherical_dist(ipoint)
    END IF
    ! discard triangles:
    ntri = tri%nentries
    j = 0
    LOOP6 : DO
      IF (j>=ntri)  EXIT LOOP6
      ! remove all ghost point triangles or triangles that will not
      ! hold the global Delaunay condition [Jacobson2013]:
      IF ((tri%a(j)%oedge /= -1) .OR. &
        & ((tri%a(j)%cap_distance+2._wp*tri%a(j)%r) > rlimit)) THEN
        ntri     = ntri - 1
        tri%a(j) = tri%a(ntri)
      ELSE
        CALL tri%a(j)%normalize_indices(pxyz)
        j = j + 1
      END IF
    END DO LOOP6
    CALL tri%resize(ntri)
    CALL pxyz%destructor()
  END SUBROUTINE triangulate


  ! --------------------------------------------------------------------
  !> Triangulation subroutine - multi-threaded version.
  !
  SUBROUTINE triangulate_mthreaded(points, tri, subset, nthreads)

    TYPE(t_point_list),    INTENT(IN)            :: points
    TYPE(t_triangulation), INTENT(INOUT), TARGET :: tri
    TYPE(t_spherical_cap), INTENT(IN)            :: subset
    INTEGER,               INTENT(IN)            :: nthreads
    ! local variables

    !> triangle count that triggers clean-up of "tri" data structure
    INTEGER, PARAMETER :: cleanup_limit  = 50000  
    INTEGER, PARAMETER :: MAX_EDGES      = 1000   !< size of edge buffer

    ! For the triangles that are created from the edges 0,1,2: compute
    ! the edge opposite to "ghost point" (if "ghost point" part of this
    ! triangle tri[j]). This lookup table is best figured out on a
    ! sheet of paper...
    INTEGER, PARAMETER :: next_oedge(0:2,-1:2) = &
      &  RESHAPE( (/ -1, -1, -1, &
      &              -1,  2,  1, &
      &               1, -1,  2, &
      &               2,  1, -1 /), (/3,4/) )
    TYPE(t_point_list)        :: pxyz
    TYPE(t_point)             :: ipoint
    TYPE(t_edge)              :: bedge
    TYPE(t_triangle), POINTER :: jtri
    INTEGER                   :: i,j, npts, gp, ninvalid, j0, ndiscard, ipt, ne, ntri, &
      &                          new_ninvalid, new_ndiscard, ithread, jmin,            &
      &                          oedge, jmin0, this_edge,k, endtri
    LOGICAL                   :: inside
    REAL(wp)                  :: cos_radius, rlimit

    ! edges buffer for re-Delaunay step (together with
    ! "edge-index-opposite-to-ghost-point")
    TYPE(t_edge)              :: edges(0:(MAX_EDGES-1))
    INTEGER                   :: edges_oedge(0:(MAX_EDGES-1))
    LOGICAL                   :: edges_valid(0:(MAX_EDGES-1))
                              
    INTEGER                   :: jmin_t

    ! copy all points since they will be re-ordered:
    npts = points%nentries
    CALL pxyz%initialize()
    pxyz = point_list(points)
    DO i=0,(npts-1)
      pxyz%a(i)        = pxyz%a(i)/pxyz%a(i)%norm2()
      ! --- guarantee the exact evaluation of determinants in the
      ! --- algorithm (losing, however, accuracy):
      pxyz%a(i)%x      = REAL(REAL(pxyz%a(i)%x, ACCURACY), wp)
      pxyz%a(i)%y      = REAL(REAL(pxyz%a(i)%y, ACCURACY), wp)
      pxyz%a(i)%z      = REAL(REAL(pxyz%a(i)%z, ACCURACY), wp)
      pxyz%a(i)%ps     = pxyz%a(i) * subset%sorting_direction
    END DO
    CALL pxyz%quicksort()

    ! assume triangulation covers entire sphere (data structure grows
    ! dynamically if not)
    CALL tri%reserve(2*npts-4 + cleanup_limit)
    tri%nentries = 0

    ! set up the initial triangle: insert the first three points
    IF (ccw_spherical(pxyz%a(0), pxyz%a(1), pxyz%a(2)) > 0._wp) THEN
      CALL tri%push_back(triangle(0, 2, 1)) ! counter-clockwise starting triangle
    ELSE 
      CALL tri%push_back(triangle(0, 1, 2)) ! clockwise starting triangle
    END IF

    ! connect triangle edges to a "ghost point"
    gp = npts
    CALL pxyz%push_back(point(-1._wp,-1._wp,-1._wp))
    ! add three initial "ghost triangles"
    CALL tri%push_back(triangle(tri%a(0)%p(1),tri%a(0)%p(0),gp, 0))
    CALL tri%push_back(triangle(tri%a(0)%p(2),tri%a(0)%p(1),gp, 0))
    CALL tri%push_back(triangle(tri%a(0)%p(0),tri%a(0)%p(2),gp, 0))
    CALL tri%a(0)%compute_circumcenter(pxyz, subset)
    
    ninvalid = 0       ! no. of invalid triangles in "tri" DATA structure
    j0       = 0       ! first triangle index that is not "complete"
    ndiscard = 0       ! no. of triangles that violate global Delaunay condition

    ipt      = 3           ! current point index
    ipoint   = pxyz%a(ipt) ! current point
    IF (subset%radius < 0._wp) THEN
      cos_radius = 99._wp
    ELSE
      cos_radius = -1._wp*COS(subset%radius)
    END IF

    ! Include each point one at a time into the existing mesh. Abort if
    ! all points have been processed or the requested radius has been
    ! reached without violating the global Delaunay condition.
    LOOP : DO ipt=3,(npts-1)
      IF (.NOT. ((ipoint%ps <= cos_radius) .OR. (ndiscard>0)))  EXIT LOOP

      ipoint = pxyz%a(ipt)
      ! Set up the edge buffer.
      !
      ! If the point pxyz[i] lies inside the circumcircle then the
      ! three edges of that triangle are added to the edge buffer and
      ! that triangle is removed.

      ne           = 0 ! no. of entries in edge buffer
      ntri         = tri%nentries
      jmin_t       = ntri
      new_ninvalid = 0
      new_ndiscard = 0
!$omp parallel private(ithread,jtri,jmin,j,inside, oedge, jmin0, bedge, this_edge,k) &
!$omp          reduction(min:jmin_t)
      ithread = 0
!$    ithread = omp_get_thread_num()
      jmin    = ntri
!$omp do reduction(+:new_ninvalid,new_ndiscard) 
      DO j=j0,(ntri-1)
        jtri => tri%a(j)
        ! test if triangle remains on "discard" list
        IF (ipoint%ps > jtri%rdiscard) THEN
          new_ndiscard  = new_ndiscard-1
          jtri%rdiscard = 1._wp
        END IF

        IF (jtri%complete == 0) THEN
          oedge = jtri%oedge
          jmin0 = jmin
          IF (jmin == ntri)  jmin = j
          IF (oedge /= -1) THEN
            ! boundary edge
            bedge  = jtri%edge(oedge)
            inside = (ccw_spherical(pxyz%a(bedge%p1), pxyz%a(bedge%p2), ipoint) <= 0._wp)
          ELSE
            ! interior edge
            inside = circum_circle_spherical(ipoint, pxyz, jtri%p)
            IF (.NOT. inside .AND. ((ipoint%ps - jtri%cc%ps) > jtri%r)) THEN
              jtri%complete =     1
              jmin          = jmin0
            END IF
          END IF
          
          IF (inside) THEN
            ! push triangle's edges onto edge list
!$omp atomic capture
            this_edge = ne
            ne = ne + 3
!$omp end atomic
            DO k=0,2
              edges(this_edge)       = jtri%edge(k)
              edges_oedge(this_edge) = next_oedge(k,oedge)
              edges_valid(this_edge) = .TRUE.
              this_edge = this_edge + 1
            END DO
            
            jtri%complete =       2 ! mark triangle for removal
            new_ninvalid = new_ninvalid + 1
            IF (jtri%rdiscard < 1._wp)  new_ndiscard=new_ndiscard-1
          END IF
        END IF
      END DO
!$omp end do
      IF (jmin >= j0) THEN
        jmin_t = jmin
      ELSE
        jmin_t = ntri
      END IF
      
      ! remove multiple edges
!$omp do
      DO j=0,(ne-1)
        IF (edges_valid(j)) THEN
          DO k=j+1,(ne-1)
            IF (edges(j) == edges(k)) THEN
              edges_valid(j) = .FALSE.
              edges_valid(k) = .FALSE.
            END IF
          END DO
        END IF
      END DO
!$omp end do
!$omp end parallel

      ninvalid = ninvalid + new_ninvalid
      IF (jmin_t /= ntri)  j0 = jmin_t

      ! form new triangles for the current point. All edges are
      ! arranged in clockwise order.
      DO j=0,(ne-1)
        IF (edges_valid(j))  CALL tri%push_back( triangle(edges(j)%p1, edges(j)%p2, ipt, edges_oedge(j)) )
      END DO
      ndiscard     = ndiscard + new_ndiscard
      new_ndiscard = 0

!$omp parallel do reduction(+:new_ndiscard)
      DO j=ntri,(tri%nentries-1)
          IF (tri%a(j)%oedge == -1) THEN  
            ! triangle without ghost point
            CALL tri%a(j)%compute_circumcenter(pxyz, subset)
            ! mark triangle if it violates the global Delaunay condition:
            IF ((cos_radius > ipoint%ps) .AND. (tri%a(j)%cap_distance < subset%radius)) THEN 
              new_ndiscard  = new_ndiscard + 1
              tri%a(j)%rdiscard = -1._wp*COS(tri%a(j)%cap_distance+2._wp*tri%a(j)%r)
            END IF
          END IF
      END DO
!$omp end parallel do
      ndiscard = ndiscard + new_ndiscard
      
      ! remove "invalid" triangles at regular intervals:
      IF (ninvalid > cleanup_limit) THEN
        ! first, replace invalid triangles by "complete" ones
        ntri   = tri%nentries
        endtri = ntri - 1
        LOOP3 : DO j=0,(ntri-1)
          IF (tri%a(j)%complete == 2) THEN
	    ! find the next "complete" triangle from the back
            LOOP4 : DO
              IF ((tri%a(endtri)%complete==1) .OR. (endtri==j))  EXIT LOOP4
              endtri = endtri - 1
            END DO LOOP4
            IF (endtri /= j) THEN
              tri%a(j)               = tri%a(endtri)
              tri%a(endtri)%complete = 2
            ELSE
              IF (j<j0)  j0=j
              EXIT LOOP3 ! only incomplete triangles follow
            END IF
          END IF
        END DO LOOP3
        ! if there are no "complete" triangles left at the list end,
        ! simply remove the invalid triangles:
        endtri = tri%nentries - 1
        LOOP5 : DO 
          IF (j==endtri)  EXIT LOOP5
          IF (tri%a(j)%complete == 2) THEN
            tri%a(j) = tri%a(endtri)
            j      = j      - 1            
            endtri = endtri - 1
            ntri   = ntri   - 1
          END IF
          j = j + 1
        END DO LOOP5
        CALL tri%resize(ntri)
        ninvalid = 0
      END IF
    END DO LOOP

    ! compute limit radius (last point inserted into triangulation)
    IF (subset%radius < 0._wp) THEN
      rlimit = 99._wp
    ELSE
      rlimit = subset%point%spherical_dist(ipoint)
    END IF
    ! discard triangles:
    ntri = tri%nentries

    j = 0
    LOOP6 : DO
      IF (j>=ntri)  EXIT LOOP6
      ! remove all ghost point triangles or triangles that will not
      ! hold the global Delaunay condition [Jacobson2013]:
      IF ((tri%a(j)%complete==2) .OR. &
        & (tri%a(j)%oedge /= -1) .OR. &
        & ((tri%a(j)%cap_distance+2._wp*tri%a(j)%r) > rlimit)) THEN
        ntri     = ntri - 1
        tri%a(j) = tri%a(ntri)
      ELSE
        CALL tri%a(j)%normalize_indices(pxyz)
        j = j + 1
      END IF
    END DO LOOP6
    CALL tri%resize(ntri)
    CALL pxyz%destructor()
  END SUBROUTINE triangulate_mthreaded


  !> perform triangulation in parallel
  SUBROUTINE triangulate_parallel(p_local, this_rank, nparallel, nsubsets_local, tri)
    TYPE(t_point_list),    INTENT(IN)             :: p_local
    INTEGER,               INTENT(IN)             :: this_rank       !< local MPI rank no.
    INTEGER,               INTENT(IN)             :: nparallel
    INTEGER,               INTENT(IN)             :: nsubsets_local  !< triangulation tasks/rank
    TYPE(t_triangulation), INTENT(INOUT)          :: tri
    ! local variables
    CHARACTER(*), PARAMETER :: routine = modname//":triangulate_parallel"
    INTEGER                             :: itri, i,j, npivot, isubset, n1,n2
    TYPE (t_triangulation)              :: tri_list(0:(nparallel-1))
    TYPE (t_point_list)                 :: p_global, pivot_points, pts0
    TYPE (t_sphcap_list)                :: subset
    TYPE (t_point)                      :: default_direction
    
    IF (dbg_level >= 2)  WRITE (0,'(a,i0)') "# no. of local points: ", p_local%nentries

    ! --- create a copy of all points (global)

    p_global = point_list(p_local)
    DO i=0,(p_global%nentries-1)
      p_global%a(i) = p_global%a(i)/p_global%a(i)%norm2()
    END DO
    CALL p_global%sync()
    IF (dbg_level >= 2)  WRITE (0,'(a,i0)') "# no. of global points: ", p_global%nentries

    ! --- select "pivot points" out of the local point set
    !     each point generates an independent Delaunay triangulation.

    npivot = nsubsets_local
    CALL pivot_points%initialize()
    default_direction = point(-1._wp, 0._wp, 0._wp)
    j = 0
    DO i=0,(npivot-1)
      CALL pivot_points%push_back(p_local%a(j))
      pivot_points%a(i)%gindex = this_rank
      pivot_points%a(i)%ps     = pivot_points%a(i) * default_direction
      j = j + MAX(INT(p_local%nentries)/npivot,1)
    END DO

    ! --- merge pivot points into a global set

    pts0 = point_list(pivot_points)
    IF (dbg_level >= 3)  WRITE (0,'(a,i0)') "# total no. of local pivot points: ", pivot_points%nentries
    CALL pts0%sync()
    IF (dbg_level >= 1)  WRITE (0,'(a,i0)') "# total no. of pivot points: ", pts0%nentries

    ! --- build a coarse auxiliary triangulation 
    !     based on the pivot points
    IF (dbg_level >= 2)  WRITE (0,'(a)') "# build a small auxiliary triangulation"
    CALL subset%initialize()
    CALL create_thin_covering( pts0, subset, this_rank )

    ! --- parallel triangulation phase
    IF (dbg_level >= 2)  WRITE (0,'(a)') "# parallel triangulation phase"

    DO i=0,(nparallel-1) 
      CALL tri_list(i)%initialize()
    END DO

    isubset = 0
    DO 
      IF (isubset >= subset%nentries)  EXIT

      DO i=0,(nparallel-1)
        tri_list(i)%nentries = 0
      END DO
!$omp parallel do schedule(runtime)
      DO itri=isubset,(isubset + nparallel - 1)
        IF (itri >= subset%nentries)  CYCLE
        
        IF (dbg_level >= 3)  WRITE (0,'(a,f10.3)') "# subset radius: ", subset%a(itri)%radius
        CALL triangulate(p_global, tri_list(itri-isubset), subset%a(itri))
        CALL tri_list(itri-isubset)%quicksort() ! re-order triangles for global merge
        IF (dbg_level >= 1)  WRITE (0,'(a,i0,a)') "# formed ", tri_list(itri-isubset)%nentries, " triangles."
      END DO
!$omp end parallel do
      isubset = (isubset + nparallel)

      ! --- merge process: translate point indices into global ones:
      DO itri=0,(nparallel-1)
        IF (tri_list(itri)%nentries == 0) CYCLE

        n1 = tri%nentries
        n2 = tri_list(itri)%nentries
        CALL tri%merge(tri_list(itri))
      
        IF ((dbg_level >= 1) .AND. (n1 > 0)) THEN
          WRITE (0,'(a,i0,a,i0,a,i0,a)')                                                     &
            &           "# merge with triangulation ", itri,                                 &
            &           "; resulting triangulation of size ", tri%nentries,                  &
            &           " (overlap: ", (n1+n2-tri%nentries), ")"
        END IF
      END DO
    END DO
    DO i=0,(nparallel-1)
      CALL tri_list(i)%destructor()
    END DO
    CALL tri%sync()
    IF (dbg_level >= 1)  WRITE (0,'(a,i0,a)') "Resulting triangulation: ", tri%nentries, " triangles."

    ! clean up
    CALL pivot_points%destructor()
    CALL subset%destructor()
  END SUBROUTINE triangulate_parallel

  
  !> Upper bound for point cloud diameter
  REAL(wp)  FUNCTION point_cloud_diam(pts, p0)
    TYPE(t_point_list), INTENT(IN) :: pts
    TYPE(t_point),      INTENT(IN) :: p0
    ! local variables
    INTEGER       :: i
    REAL(wp)      :: dist

    dist = 0._wp
    DO i=1,(pts%nentries-1)
      dist = MAX(dist, p0%spherical_dist(pts%a(i)))
    END DO
    point_cloud_diam = dist ! upper bound
  END FUNCTION point_cloud_diam


  !> For a given list of points on the sphere: Compute radii, such
  !  that the spherical cap at these points cover the whole sphere.
  !
  SUBROUTINE create_thin_covering(point_set, subset, this_rank)
    TYPE(t_point_list),  INTENT(IN)    :: point_set
    TYPE(t_sphcap_list), INTENT(INOUT) :: subset
    INTEGER,             INTENT(IN)    :: this_rank !< local MPI rank no.
    ! local variables
    CHARACTER(*), PARAMETER :: routine = modname//":create_thin_covering"
    INTEGER               :: i, j, j1, j2, p1, p1_local, ierrstat
    INTEGER, ALLOCATABLE  :: local_idx(:)
    TYPE(t_triangulation) :: tri0
    TYPE(t_point_list)    :: pts0

    ! triangulate the partition points and use this triangulation to
    ! define the subset radii:
    ALLOCATE(local_idx(0:(point_set%nentries-1)), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")    

    pts0 = point_list(point_set)
    j = 0
    DO i=0,(point_set%nentries-1)
      IF (point_set%a(i)%gindex == this_rank) THEN
        CALL subset%push_back(spherical_cap(point_set%a(i), 0._wp))
        local_idx(i) = j
        j = j + 1
      END IF
      pts0%a(i)%gindex = i
    END DO
    CALL tri0%initialize()
    CALL triangulate(pts0, tri0, spherical_cap(point(-1._wp, 0._wp, 0._wp), -1._wp))

    IF ((dbg_level >= 3) .AND. (this_rank == 0)) THEN
      ! debugging output to file
      CALL write_triangulation_vtk("tri0.vtk", pts0, tri0)
    END IF
    
    ! radii are set to max. length of coarse Delaunay triangulation:
    DO i=0,(tri0%nentries-1)
      DO j1=0,2
        p1 = tri0%a(i)%p(j1)
        IF (point_set%a(p1)%gindex == this_rank) THEN
          p1_local = local_idx(p1)
          DO j2=0,2
            IF (j1 /= j2) THEN
              subset%a(p1_local)%radius =         &
                &  MAX(subset%a(p1_local)%radius, &
                &      point_set%a(p1)%spherical_dist(point_set%a(tri0%a(i)%p(j2))))
            END IF
          END DO
        END IF
      END DO
    END DO

    ! clean up
    DEALLOCATE(local_idx, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")    
  END SUBROUTINE create_thin_covering


  !> Write out triangulation in VTK ASCII format.
  !
  SUBROUTINE write_triangulation_vtk(filename, p, tri)
    CHARACTER(LEN=*)      :: filename
    TYPE(t_point_list)    :: p
    TYPE(t_triangulation) :: tri
    ! local variables
    INTEGER :: out_unit, i
    
    ! write triangles in vtk format
    out_unit=20
    OPEN (unit=out_unit,file=TRIM(filename),action="write",status="replace")
    WRITE (out_unit,'(a)') "# vtk DataFile Version 2.0"
    WRITE (out_unit,'(a)') "# delaunay example"
    WRITE (out_unit,'(a)') "ASCII"
    WRITE (out_unit,'(a)') "DATASET UNSTRUCTURED_GRID"
    WRITE (out_unit,'(a,i0,a)') "POINTS ", p%nentries, " float"
    DO i=0,(p%nentries-1)
      WRITE (out_unit,'(3f12.8)') p%a(i)%x, p%a(i)%y, p%a(i)%z
    END DO
    WRITE (out_unit,'(a,i0,a,i0)') "CELLS ", tri%nentries," ", 4*tri%nentries
    DO i=0,(tri%nentries-1)
      WRITE (out_unit,'(3(a,i0))') "3 ", tri%a(i)%p(0)," ",tri%a(i)%p(1)," ",tri%a(i)%p(2)
    END DO
    WRITE (out_unit,'(a,i0)') "CELL_TYPES ", tri%nentries
    DO i=0,(tri%nentries-1)
      WRITE (out_unit,'(a)') "5"
    END DO
    CLOSE(out_unit)
  END SUBROUTINE write_triangulation_vtk
  
END MODULE mo_delaunay
