!>
!!               The subroutines in this module compute.
!!
!!               The subroutines in this module compute
!!               all the geometric information for the grid hierarchy
!!               in the internal data representation used by the grid generator.
!!
!! @par Revision History
!! Initial version  by:
!! @par
!! Luis Kornblueh,  MPI-M, Hamburg, January 2004
!! @par
!! Luis Kornblueh, MPI-M, Hamburg, February 2004
!!  - change to compiling and running version
!!  - create level 0 subdivision, connect triangle pointers
!! Luis Kornblueh, MPI-M, Hamburg, March 2004
!!  - include parents
!!  - full working version of tree and graph representation
!! Luca Bonaventura,  MPI-M, Hamburg, October 2004
!!  - include documentation and Protex headers
!! Luca Bonaventura,  MPI-M, Hamburg, April 2005
!!  - revision to account for patches
!! Luca Bonaventura,  MPI-M, Hamburg, August 2005
!!  - revision to include grid optimization and
!!    correct computation of normal vectors to edges
!! Thomas Heinze,  DWD, Offenbach, 2005-09
!!  - cleaning up code according to latest programming guide
!!  - debugging for grid optimization: changed index of optimize_heikes
!!  Modifications by Th.Heinze, DWD (2006-10-25):
!!  - changed index to idx in TYPE declarations of edge_info, vertex_info,
!!    triangle_info, grid_cells, grid_edges and grid_vertices
!!  Modification by P. Ripodas, DWD, (2007-02):
!!  - The edge%system_orientation is calculated (system orientation
!!  of (v2-v1), (c2-c1)
!! Almut Gassmannn, MPI-M, (2007-03)
!!  - Equal area subdivision and small circle subdivision parts added
!! Luis Kornblueh, MPI-M, March 2008
!!  - removed unnecessary vertex0/1/2 references
!! Guenther Zaengl, DWD, 2008-12-08
!!  - add option for stopping grid optimization at a certain level
!!    (needed for combining spring dynamics with nesting)
!! Almut Gassmann, MPI-M (2009-01-15)
!!  - adding a planar option
!!  - cleaning and restructuring
!!  - Earth dimensions used already here for avoiding the patch generator
!!    in every case
!! Almut Gassmann, MPI-M (2009-03-02)
!!  - tuning parameter beta_spring for spring dynamics comes via namelist input
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!!
MODULE mo_geometry
  !-------------------------------------------------------------------------
  !
  !    ProTeX FORTRAN source: Style 2
  !    modified for ICON project, DWD/MPI-M 2005
  !
  !-------------------------------------------------------------------------
  !
  !
  !
  USE mo_kind,           ONLY: wp
  USE mo_exception,      ONLY: message_text, message, finish
  USE mo_math_constants, ONLY: pi

  USE mo_math_utilities, ONLY: t_cartesian_coordinates, vector_product, &
    & circum_center, cc2gc, arc_length,      &
    & spherical_intersection2, spherical_intersection2_v, t_geographical_coordinates

  USE mo_physical_constants, ONLY: earth_radius

  USE mo_base_datatypes, ONLY: t_triangle, t_spheres, dummy_e

  USE mo_topology,       ONLY: root, kroot, up_down, number_of_levels, &
    & icosahedron, spheres_on_levels

  USE mo_icosahedron_geometry, ONLY: init_icosahedron_vertices,     &
    & get_icosahedron_vertex, init_planar_vertices

  ! USE mo_grid,           ONLY: grid

  USE mo_grid_levels,    ONLY: init_grid,   &
    & itype_optimize, &
     ! = 0 : natural grid
     ! = 1 : Heikes optimization
     ! = 2 : equal area subdivision
     ! = 3 : c-grid small circle constraint
     ! = 4 : spring dynamics
     ! = 5 : spring dynamics with convergence acceleration
    & l_c_grid, &
     ! = T : last subdivision step uses
     !       C-grid small circle constraint
    & maxlev_optim, &
     ! highest level for which grid optimization
     ! is executed (for itype_optimize=1 only)
    & tria_arc_km, & ! arc length in km for finest grid
    & beta_spring ! tuning parameter for the spring
  ! dynamics optimization (refer to
  ! the corresponding subroutine)

  USE mo_optimize,       ONLY: optimize_heikes

  USE mo_small_circle_c, ONLY: scc_inner, scc_outer, scc_inner_root, &
    & scc_prepare

  USE mo_equal_area,     ONLY: equal_area, equal_area_nine, equal_area_twen5

  USE mo_spring_dynamics, ONLY: spring_dynamics

!   USE mo_impl_constants,  ONLY: min_rlcell, max_rlcell, &
!     & min_rlvert, max_rlvert, &
!     & min_rledge, max_rledge

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: generate_geometry, init_grid

  ! LOGICAL :: lgrid_initialized = .false.

CONTAINS

  !-------------------------------------------------------------------------
  !
  !
  !>
  !!  This is the main driver routine for the generation of the geometric.
  !!
  !!  This is the main driver routine for the generation of the geometric
  !!  information. The gridpoints of the original icosahedron (level -1) are
  !!  computed first by <i>init_icosahedron_coordinates</i>. The gridpoints of
  !!  the root grid (level 0) are computed by  <i>init_root_coordinates</i>.
  !!  The vertices of the grid cells are then computed recursively by
  !!  <i>build_vertices.</i> The remaining geometric quantities are computed
  !!  for each level by  <i>calculate_geometry</i>.
  !!
  !! @par Revision History
  !!  Luis Kornblueh, MPI-M, Hamburg, 2004-03
  !!  Thomas Heinze,  DWD, Offenbach, 2005-09
  !!  - cleaning up code according to latest programming guide
  !!  - debugging for grid optimization: changed index of optimize_heikes
  !!  Almut Gassmann, MPI-M, Hamburg, 2007-03
  !!  - added parts for small circle subdivision
  !!  Almut Gassmann, MPI-M, Hamburg, (2009-01-09)
  !!  - added part for planar grid construction
  !!
  SUBROUTINE generate_geometry(lplane)
    !

    LOGICAL, INTENT(in) :: lplane

    TYPE(t_spheres), POINTER :: ptr_tl
    TYPE(t_triangle), POINTER :: ptr_tr

    INTEGER :: ji, j, jl
    INTEGER :: i_nots, i_c_grid
    LOGICAL :: lopt

    !-------------------------------------------------------------------------
    !BOC

    CALL init_icosahedron_coordinates(lplane)

    IF (.not. lplane) THEN

      ! Decide whether the last level should be treated with the C-grid constraint
      ! The number i_c_grid is then subtracted from number_of_levels.

      IF (l_c_grid) THEN
        i_c_grid=1
      ELSE
        i_c_grid=0
      ENDIF

      ! now fill data structures

      SELECT CASE (itype_optimize)

      CASE (0)
        !
        !  build unoptimized grid beginning with root level
        !  up to number_of_levels
        !

        CALL init_root_coordinates

        DO j = 0, 19
          DO ji = 0, kroot*kroot-1
            CALL build_vertices (root(ji,j)%t,number_of_levels-i_c_grid,lplane)
          ENDDO
        ENDDO

      CASE (1)
        !
        !  Heikes grid optimization is called here
        !
        !  first build unoptimized grid on root level
        !  to simplify things, this level will not be
        !  optimized (it is never used for numerical computations)
        !
        !  note that build_vertices has been
        !  changed, it has a second argument
        !  that tells how deep the recursion is allowed
        !  to go in the grid hierarchy
        !
        CALL init_root_coordinates
        DO j = 0, 19
          DO ji = 0, kroot*kroot-1
            CALL build_vertices (root(ji,j)%t,0,lplane)
          ENDDO
        ENDDO
        !
        !  now go through all levels from 1 to number_of_levels
        !
        DO jl = 1, number_of_levels-i_c_grid

          ptr_tl => spheres_on_levels(jl-1)
          i_nots = ptr_tl%no_triangles
          !
          ! 1. for each level build the unoptimized grid on that level
          !
          DO j = 0, i_nots-1
            ptr_tr=>ptr_tl%ts(j)
            CALL build_vertices (ptr_tr, jl,lplane)
          ENDDO
          !
          ! 2. perform the optimization on that level
          ! thh:
          ! changed index from jl-1 to jl (no optimization on level 0,
          ! but optimization on finest level
          !
          IF (jl <= maxlev_optim) THEN
            CALL optimize_heikes(spheres_on_levels(jl))
          ENDIF

        ENDDO

      CASE (2)
        !
        ! Equal area subdivision
        !
        ! start with the root level
        !
        CALL init_root_coordinates

        DO j=0, 19
        SELECT CASE (kroot)
          CASE (2)
          CALL equal_area(icosahedron(j)%t)
          CASE (3)
            CALL equal_area_nine(j)
          CASE (5)
            CALL equal_area_twen5(j)
          CASE default
            WRITE (message_text,'(a,i2)') &
              & 'No equal area subdivision implemented for kroot= ', kroot
            CALL finish('',TRIM(message_text))
          END SELECT
        ENDDO
        !
        ! do subdivision on refined triangles
        !
        DO jl = 1, number_of_levels-i_c_grid

          ptr_tl => spheres_on_levels(jl-1)
          i_nots = ptr_tl%no_triangles

          DO j = 0, i_nots-1
            ptr_tr=>ptr_tl%ts(j)
            CALL equal_area(ptr_tr)
          ENDDO

        ENDDO

        CASE (3)
        !
        ! C-grid constraint for subdivision
        !
        CALL init_root_coordinates

        DO j=0, 19
          CALL scc_inner_root(j)
        ENDDO
        !
        DO jl = 1, number_of_levels
          !
          ptr_tl => spheres_on_levels(jl-1)
          i_nots = ptr_tl%no_triangles
          !
          ! subdivide all triangles with small circle subdivision
          !
          DO j = 0, i_nots-1
            ptr_tr=>ptr_tl%ts(j)
            CALL scc_inner(ptr_tr)
          ENDDO
          !
          ! adjust the outer boundary of the triangles, too
          !
          CALL scc_outer(ptr_tl)
          !
        ENDDO

        CASE (4, 5)
          !
          ! Spring dynamics
          !
          ! 1) build unoptimized grid beginning with root level
          !    up to number_of_levels
          !

          CALL init_root_coordinates

          DO j = 0, 19
            DO ji = 0, kroot*kroot-1
              CALL build_vertices (root(ji,j)%t,0,lplane)
            ENDDO
          ENDDO

          IF (itype_optimize == 5) THEN
            lopt = .TRUE.  ! use convergence acceleration
          ELSE
            lopt = .FALSE.
          ENDIF
          !
          ! 2) Call spring dynamics
          !
          !  now go through all levels from 1 to number_of_levels
          !
          DO jl = 1, number_of_levels-i_c_grid

            ptr_tl => spheres_on_levels(jl-1)
            i_nots = ptr_tl%no_triangles
            !
            ! 1. for each level build the unoptimized grid on that level
            !
            DO j = 0, i_nots-1
              ptr_tr=>ptr_tl%ts(j)
              CALL build_vertices (ptr_tr, jl,lplane)
            ENDDO
            !
            ! 2. perform the optimization on that level
            !
            IF (jl <= maxlev_optim) THEN
              CALL spring_dynamics(spheres_on_levels(jl),beta_spring,lopt)
            ENDIF

          ENDDO


        CASE default
          !
          ! give an error message
          !
          WRITE (message_text,'(a,i2)') &
            & 'No grid construction available for itype_optimize= ', itype_optimize
          CALL finish('',TRIM(message_text))

        END SELECT

        !
        ! If the last level should be treated with the C-grid constraint...
        !
        IF (l_c_grid .and. itype_optimize/=3 ) THEN
          !
          ptr_tl => spheres_on_levels(number_of_levels-1)
          i_nots = ptr_tl%no_triangles
          !
          ! Prepare values needed for c-grid subdivision,
          ! if they are not yet available
          !
          IF (itype_optimize <=1 .or. itype_optimize>=4) THEN
            CALL scc_prepare(ptr_tl)
          ENDIF
          !
          ! subdivide all triangles with small circle subdivision
          !
          DO j = 0, i_nots-1
            ptr_tr=>ptr_tl%ts(j)
            CALL scc_inner(ptr_tr)
          ENDDO
          !
          ! adjust the outer boundary of the triangles, too
          !
          CALL scc_outer(ptr_tl)
          !
        ENDIF

        !
        !  for all levels calculate geometric information and print some
        !  in GMT output format
        !

        DO j = 0, number_of_levels

          WRITE(message_text,'(a,i0)') 'Compute geometry on level ', j
          CALL message('', TRIM(message_text))

          CALL compute_geometry(j)

        ENDDO

      ELSEIF(lplane) THEN

        ! The planar grid. Treat differently because the coordinates give not
        ! really a connected surface in coordinate space.
        DO j = 0, 7
          CALL init_root_coordinates_plane(j)
          DO ji = 0, kroot*kroot-1
            CALL build_vertices (root(ji,j)%t,number_of_levels,lplane)
          ENDDO
        ENDDO

        DO j = 0, number_of_levels

          WRITE(message_text,'(a,i0)') 'Compute geometry on level ', j
          CALL message('', TRIM(message_text))

          CALL compute_planar_geometry(j)

        ENDDO

      ENDIF

    END SUBROUTINE generate_geometry

    !-------------------------------------------------------------------------
    !
    !
    !>
    !! Initializes Cartesian coordinates of the vertices of the original icosahedron.
    !!
    !! Initializes Cartesian coordinates of the vertices of the original icosahedron.
    !!
    !! @par Revision History
    !!  Luis Kornblueh, MPI-M, Hamburg, March 2004
    !!  Almut Gassmann, MPI-M, Hamburg, March 2007
    !!   - added parts for small circle subdisvision
    !!  Almut Gassmann, MPI-M, Hamburg, (2008-09)
    !!   - dummy edge as a target
    !!
    SUBROUTINE init_icosahedron_coordinates(lplane)
      !

      LOGICAL, INTENT(in) :: lplane
      INTEGER :: j

      ! initialize Cartesian coordinates of the icosahedron vertices
      !-------------------------------------------------

      IF (.not. lplane) THEN
        CALL init_icosahedron_vertices
      ELSE
        CALL init_planar_vertices
      ENDIF

      ! store in basic triangle data structure

      IF (itype_optimize ==2 .or. itype_optimize == 3) THEN

        dummy_e%cos_rel_co_lat=0.0_wp
        dummy_e%sc_pole%x(:)=0.0_wp
        dummy_e%edge_primal_arc=2.0_wp*acos(1.0_wp/(2.0_wp*sin(0.2_wp*pi)))

        DO j = 0, 19

          icosahedron(j)%t%triangle_area=0.2_wp*pi
          icosahedron(j)%t%edge0 => dummy_e
          icosahedron(j)%t%edge1 => dummy_e
          icosahedron(j)%t%edge2 => dummy_e

        END DO

      ENDIF

    END SUBROUTINE init_icosahedron_coordinates

    !-------------------------------------------------------------------------

    !>
    !!               Initializes Cartesian coordinates of vertices of.
    !!
    !!               Initializes Cartesian coordinates of vertices of
    !!              the triangles on the root grid.
    !!
    !! @par Revision History
    !!  Luis Kornblueh, MPI-M, Hamburg, March 2004
    !!
    SUBROUTINE init_root_coordinates
      !

      !-------------------------------------------------------------------------

      INTEGER :: ji, j, jcn, ict, iu, id, iorientation

      TYPE(t_cartesian_coordinates) :: v(0:2), va(0:kroot,0:kroot)
      REAL(wp) :: zchord   ! Cartesian distance between x1 and x2
      REAL(wp) :: ztheta   ! great circle angle between x1 and x2
      REAL(wp) :: zalpha, zbeta, zgamma, za, zb, zg

      ! fraction of distance between two vertices as part of subtriangle

      DO jcn = 0, 19

        v(0) = get_icosahedron_vertex(0,jcn)
        v(1) = get_icosahedron_vertex(1,jcn)
        v(2) = get_icosahedron_vertex(2,jcn)

        va(0,0)         = v(0)
        va(kroot,0)     = v(1)
        va(kroot,kroot) = v(2)

        DO ji = 0, kroot

          ! Calculate the weighting factors which follow from the condition
          ! that x is a point on the unit-sphere, too.

          zchord = SQRT(SUM((v(0)%x-v(1)%x)**2))
          ztheta = 2._wp*asin (0.5_wp*zchord)
          zgamma = REAL(ji,wp)/REAL(kroot,wp)
          zbeta  = SIN (zgamma*ztheta)/SIN (ztheta)
          zalpha = SIN ((1.0_wp-zgamma)*ztheta)/SIN (ztheta)
          va(ji,0)%x = zalpha*v(0)%x+zbeta*v(1)%x

          zchord = SQRT(SUM((v(0)%x-v(2)%x)**2))
          ztheta = 2.0_wp*asin (0.5_wp*zchord)
          zgamma = REAL(ji,wp)/REAL(kroot,wp)
          zbeta  = SIN (zgamma*ztheta)/SIN (ztheta)
          zalpha = SIN ((1.0_wp-zgamma)*ztheta)/SIN (ztheta)
          va(ji,ji)%x = zalpha*v(0)%x+zbeta*v(2)%x

          IF (ji < 2) CYCLE

          DO j = 1, ji-1
            zchord = SQRT(SUM((va(ji,0)%x-va(ji,ji)%x)**2))
            ztheta = 2._wp*asin (0.5_wp*zchord)
            zg = REAL(j,wp)/REAL(ji,wp)
            zb = SIN (zg*ztheta)/SIN (ztheta)
            za = SIN ((1.0_wp-zg)*ztheta)/SIN (ztheta)
            va(ji,j)%x = za*va(ji,0)%x+zb*va(ji,ji)%x
          ENDDO

        ENDDO

        ! now restore in root structure vertex 0 looking sub triangles

        iorientation = IAND(root(0,jcn)%t%info%orientation,up_down)

        DO ji = 0, kroot-1
          iu = 0           ! counter of vertices in a line
          id = 1           ! counter of vertices in a line
          DO j = 0, 2*ji

            ict = ji*ji+j

            IF (IAND(root(ict,jcn)%t%info%orientation,up_down) == iorientation) THEN
              root(ict,jcn)%t%vertex0%vertex = va(ji  ,iu  )
              root(ict,jcn)%t%vertex1%vertex = va(ji+1,iu  )
              root(ict,jcn)%t%vertex2%vertex = va(ji+1,iu+1)
              iu = iu+1
            ELSE
              root(ict,jcn)%t%vertex0%vertex = va(ji+1,id  )
              root(ict,jcn)%t%vertex1%vertex = va(ji  ,id  )
              root(ict,jcn)%t%vertex2%vertex = va(ji  ,id-1)
              id = id+1
            END IF

          ENDDO
        ENDDO

      ENDDO

    END SUBROUTINE init_root_coordinates
    !-------------------------------------------------------------------------

    !>
    !!               Initializes Cartesian coordinates of vertices of.
    !!
    !!               Initializes Cartesian coordinates of vertices of
    !!              the triangles on the root grid. Planar case.
    !!
    !! @par Revision History
    !! Almut Gassmann, MPI-M, (2009-01-12)
    !! -adaption for the plane
    !!
    SUBROUTINE init_root_coordinates_plane(jcn)
      !

      !-------------------------------------------------------------------------
      INTEGER, INTENT(in) :: jcn
      INTEGER :: ji, j, ict, iu, id, iorientation

      TYPE(t_cartesian_coordinates) :: v(0:2), va(0:kroot,0:kroot)
      REAL(wp) :: zalpha, zbeta, zgamma, za, zb, zg

      ! fraction of distance between two vertices as part of subtriangle

      v(0) = get_icosahedron_vertex(0,jcn)
      v(1) = get_icosahedron_vertex(1,jcn)
      v(2) = get_icosahedron_vertex(2,jcn)

      va(0,0)         = v(0)
      va(kroot,0)     = v(1)
      va(kroot,kroot) = v(2)

      DO ji = 0, kroot

        ! Calculate the weighting factors which follow from the condition
        ! that x is a point on the unit-sphere, too.

        zgamma = REAL(ji,wp)/REAL(kroot,wp)
        zbeta  = zgamma
        zalpha = 1.0_wp-zgamma
        va(ji,0)%x = zalpha*v(0)%x+zbeta*v(1)%x

        zgamma = REAL(ji,wp)/REAL(kroot,wp)
        zbeta  = zgamma
        zalpha = 1.0_wp-zgamma
        va(ji,ji)%x = zalpha*v(0)%x+zbeta*v(2)%x

        IF (ji < 2) CYCLE

        DO j = 1, ji-1
          zg = REAL(j,wp)/REAL(ji,wp)
          zb = zg
          za = 1.0_wp-zg
          va(ji,j)%x = za*va(ji,0)%x+zb*va(ji,ji)%x
        ENDDO

      ENDDO

      ! now restore in root structure vertex 0 looking sub triangles

      iorientation = IAND(root(0,jcn)%t%info%orientation,up_down)

      DO ji = 0, kroot-1
        iu = 0           ! counter of vertices in a line
        id = 1           ! counter of vertices in a line
        DO j = 0, 2*ji

          ict = ji*ji+j

          IF (IAND(root(ict,jcn)%t%info%orientation,up_down) == iorientation) THEN
            root(ict,jcn)%t%vertex0%vertex = va(ji  ,iu  )
            root(ict,jcn)%t%vertex1%vertex = va(ji+1,iu  )
            root(ict,jcn)%t%vertex2%vertex = va(ji+1,iu+1)
            iu = iu+1
          ELSE
            root(ict,jcn)%t%vertex0%vertex = va(ji+1,id  )
            root(ict,jcn)%t%vertex1%vertex = va(ji  ,id  )
            root(ict,jcn)%t%vertex2%vertex = va(ji  ,id-1)
            id = id+1
          END IF

        ENDDO
      ENDDO

    END SUBROUTINE init_root_coordinates_plane
    !-------------------------------------------------------------------------
    !>
    !!               Initializes Cartesian coordinates of vertices of.
    !!
    !!               Initializes Cartesian coordinates of vertices of
    !!              the triangles on higher grid levels (1 and above).
    !!              The computation of all the other geometrical properties is done
    !!              in a subsequent step because of
    !!              the recursive nature of this routine.
    !!
    !! @par Revision History
    !!
    RECURSIVE SUBROUTINE build_vertices (cn,klev,lplane)
      !

      ! Original version by Luis Kornblueh, MPI-M, Hamburg, March 2004.
      ! Modified by Luca Bonaventura, MPi-M, Hamburg, August 2005
      ! with inclusion of second argument klev to stop recursion at
      ! a given level (is used in the implementation of the grid
      ! optimization.
      !
      !
      !


      INTEGER, INTENT(in) :: klev

      TYPE(t_triangle), POINTER :: cn

      LOGICAL, INTENT(in) :: lplane

      TYPE(t_cartesian_coordinates) :: v(0:2), vn(0:2), w


      REAL(wp) :: ztheta, zchord
      REAL(wp) :: zalpha, zbeta, zgamma

      !-------------------------------------------------------------------------
      !BOC

      IF (lplane) THEN
        CALL plane_edge_cell_coords(cn)
      ENDIF

      IF (cn%info%level < klev) THEN

        v(0) = cn%vertex0%vertex
        v(1) = cn%vertex1%vertex
        v(2) = cn%vertex2%vertex

        w%x(1) = v(1)%x(1)-v(2)%x(1)
        w%x(2) = v(1)%x(2)-v(2)%x(2)
        w%x(3) = v(1)%x(3)-v(2)%x(3)
        zchord = SQRT(DOT_PRODUCT(w%x,w%x))
        IF (.not.lplane) THEN
          ztheta = 2._wp*asin (0.5_wp*zchord)
          zgamma = 0.5_wp
          zbeta  = SIN (zgamma*ztheta)/SIN (ztheta)
          zalpha = SIN ((1.0_wp-zgamma)*ztheta)/SIN (ztheta)
        ELSE
          zgamma = 0.5_wp
          zbeta  = zgamma
          zalpha = 1.0_wp-zgamma
        ENDIF
        vn(1)%x(1) = zalpha*v(1)%x(1)+zbeta*v(2)%x(1)
        vn(1)%x(2) = zalpha*v(1)%x(2)+zbeta*v(2)%x(2)
        vn(1)%x(3) = zalpha*v(1)%x(3)+zbeta*v(2)%x(3)

        w%x(1) = v(0)%x(1)-v(2)%x(1)
        w%x(2) = v(0)%x(2)-v(2)%x(2)
        w%x(3) = v(0)%x(3)-v(2)%x(3)
        zchord = SQRT(DOT_PRODUCT(w%x,w%x))
        IF (.not.lplane) THEN
          ztheta = 2._wp*asin (0.5_wp*zchord)
          zgamma = 0.5_wp
          zbeta  = SIN (zgamma*ztheta)/SIN (ztheta)
          zalpha = SIN ((1.0_wp-zgamma)*ztheta)/SIN (ztheta)
        ELSE
          zgamma = 0.5_wp
          zbeta  = zgamma
          zalpha = 1.0_wp-zgamma
        ENDIF
        vn(2)%x(1) = zalpha*v(0)%x(1)+zbeta*v(2)%x(1)
        vn(2)%x(2) = zalpha*v(0)%x(2)+zbeta*v(2)%x(2)
        vn(2)%x(3) = zalpha*v(0)%x(3)+zbeta*v(2)%x(3)

        w%x(1) = v(0)%x(1)-v(1)%x(1)
        w%x(2) = v(0)%x(2)-v(1)%x(2)
        w%x(3) = v(0)%x(3)-v(1)%x(3)
        zchord = SQRT(DOT_PRODUCT(w%x,w%x))
        IF (.not.lplane) THEN
          ztheta = 2._wp*asin (0.5_wp*zchord)
          zgamma = 0.5_wp
          zbeta  = SIN (zgamma*ztheta)/SIN (ztheta)
          zalpha = SIN ((1.0_wp-zgamma)*ztheta)/SIN (ztheta)
        ELSE
          zgamma = 0.5_wp
          zbeta  = zgamma
          zalpha = 1.0_wp-zgamma
        ENDIF
        vn(0)%x(1) = zalpha*v(0)%x(1)+zbeta*v(1)%x(1)
        vn(0)%x(2) = zalpha*v(0)%x(2)+zbeta*v(1)%x(2)
        vn(0)%x(3) = zalpha*v(0)%x(3)+zbeta*v(1)%x(3)

        cn%sub_triangle0%vertex0%vertex = v(0)
        cn%sub_triangle0%vertex1%vertex = vn(0)
        cn%sub_triangle0%vertex2%vertex = vn(2)

        cn%sub_triangle1%vertex0%vertex = vn(0)
        cn%sub_triangle1%vertex1%vertex = v(1)
        cn%sub_triangle1%vertex2%vertex = vn(1)

        cn%sub_triangle2%vertex0%vertex = vn(1)
        cn%sub_triangle2%vertex1%vertex = vn(2)
        cn%sub_triangle2%vertex2%vertex = vn(0)

        cn%sub_triangle3%vertex0%vertex = vn(2)
        cn%sub_triangle3%vertex1%vertex = vn(1)
        cn%sub_triangle3%vertex2%vertex = v(2)


        CALL build_vertices (cn%sub_triangle0,klev,lplane)
        CALL build_vertices (cn%sub_triangle1,klev,lplane)
        CALL build_vertices (cn%sub_triangle2,klev,lplane)
        CALL build_vertices (cn%sub_triangle3,klev,lplane)

      ENDIF

    END SUBROUTINE build_vertices
    !-------------------------------------------------------------------------
    !
    !>
    !!               Computes coordinates other than vertices for the planar grid.
    !!
    !!
    !! @par Revision History
    !! Almut Gassmann, MPI-M (2009-01-12)
    !!
    SUBROUTINE plane_edge_cell_coords(cn)
      !
      TYPE(t_triangle), POINTER :: cn

      INTEGER :: j2, j1, j0, j, i_orientation
      TYPE (t_cartesian_coordinates) :: vec_d, vec_t
      REAL (wp) :: z_tmp
      !
      !-------------------------------------------------------------------------

      !edge coords
      cn%edge0%edge_center%x = 0.5_wp*(cn%vertex0%vertex%x+cn%vertex1%vertex%x)
      cn%edge1%edge_center%x = 0.5_wp*(cn%vertex1%vertex%x+cn%vertex2%vertex%x)
      cn%edge2%edge_center%x = 0.5_wp*(cn%vertex2%vertex%x+cn%vertex0%vertex%x)

      !centers
      cn%triangle_center%x = (cn%vertex0%vertex%x+cn%vertex1%vertex%x+&
        & cn%vertex2%vertex%x)/3.0_wp

      !edge_orientations for cells
      j =cn%info%idx
      j0=cn%neighbor0%info%idx
      j1=cn%neighbor1%info%idx
      j2=cn%neighbor2%info%idx
      cn%edge0_orientation = SIGN(1,j0-j)
      cn%edge1_orientation = SIGN(1,j1-j)
      cn%edge2_orientation = SIGN(1,j2-j)

      !primal and dual normals
      vec_d%x = (cn%edge0%edge_center%x-cn%triangle_center%x)* &
        &  REAL(cn%edge0_orientation,wp)
      z_tmp = SQRT(DOT_PRODUCT(vec_d%x,vec_d%x))
      vec_d%x = vec_d%x/z_tmp
      cn%edge0%edge_primal_normal%v1 = vec_d%x(1)
      cn%edge0%edge_primal_normal%v2 = vec_d%x(2)
      cn%edge0%edge_dual_normal%v1   =  cn%edge0%edge_primal_normal%v2
      cn%edge0%edge_dual_normal%v2   = -cn%edge0%edge_primal_normal%v1

      vec_d%x = (cn%edge1%edge_center%x-cn%triangle_center%x)* &
        &  REAL(cn%edge1_orientation,wp)
      z_tmp = SQRT(DOT_PRODUCT(vec_d%x,vec_d%x))
      vec_d%x = vec_d%x/z_tmp
      cn%edge1%edge_primal_normal%v1 = vec_d%x(1)
      cn%edge1%edge_primal_normal%v2 = vec_d%x(2)
      cn%edge1%edge_dual_normal%v1   =  cn%edge1%edge_primal_normal%v2
      cn%edge1%edge_dual_normal%v2   = -cn%edge1%edge_primal_normal%v1

      vec_d%x = (cn%edge2%edge_center%x-cn%triangle_center%x)* &
        & REAL(cn%edge2_orientation,wp)
      z_tmp = SQRT(DOT_PRODUCT(vec_d%x,vec_d%x))
      vec_d%x = vec_d%x/z_tmp
      cn%edge2%edge_primal_normal%v1 = vec_d%x(1)
      cn%edge2%edge_primal_normal%v2 = vec_d%x(2)
      cn%edge2%edge_dual_normal%v1   =  cn%edge2%edge_primal_normal%v2
      cn%edge2%edge_dual_normal%v2   = -cn%edge2%edge_primal_normal%v1

      ! system orientation
      j0=cn%vertex0%info%idx
      j1=cn%vertex1%info%idx
      i_orientation = SIGN(1,j1-j0)
      vec_d%x    = (cn%vertex1%vertex%x-cn%vertex0%vertex%x)* &
        & REAL(i_orientation,wp)
      vec_t%x(1) = cn%edge0%edge_dual_normal%v1
      vec_t%x(2) = cn%edge0%edge_dual_normal%v2
      vec_t%x(3) = 0.0_wp
      z_tmp      = DOT_PRODUCT(vec_t%x,vec_d%x)
      cn%edge0%system_orientation = NINT(SIGN(1.0_wp,z_tmp))

      j1=cn%vertex1%info%idx
      j2=cn%vertex2%info%idx
      i_orientation = SIGN(1,j2-j1)
      vec_d%x    = (cn%vertex2%vertex%x-cn%vertex1%vertex%x)* &
        & REAL(i_orientation,wp)
      vec_t%x(1) = cn%edge1%edge_dual_normal%v1
      vec_t%x(2) = cn%edge1%edge_dual_normal%v2
      vec_t%x(3) = 0.0_wp
      z_tmp      = DOT_PRODUCT(vec_t%x,vec_d%x)
      cn%edge1%system_orientation = NINT(SIGN(1.0_wp,z_tmp))

      j2=cn%vertex2%info%idx
      j0=cn%vertex0%info%idx
      i_orientation = SIGN(1,j0-j2)
      vec_d%x    = (cn%vertex0%vertex%x-cn%vertex2%vertex%x)* &
        & REAL(i_orientation,wp)
      vec_t%x(1) = cn%edge2%edge_dual_normal%v1
      vec_t%x(2) = cn%edge2%edge_dual_normal%v2
      vec_t%x(3) = 0.0_wp
      z_tmp      = DOT_PRODUCT(vec_t%x,vec_d%x)
      cn%edge2%system_orientation = NINT(SIGN(1.0_wp,z_tmp))

    END SUBROUTINE plane_edge_cell_coords
    !-------------------------------------------------------------------------
    !
    !
    !>
    !!               Computes geometry of the planar grid.
    !!
    !!
    !! @par Revision History
    !! Almut Gassmann, MPI_M (2009-01-12)
    !!
    SUBROUTINE compute_planar_geometry (klev)
      !

      INTEGER, INTENT(in) :: klev

      INTEGER :: nots, noes, novs, j
      REAL(wp) :: z_primal_length
      TYPE(t_spheres), POINTER :: ptr_tl
      !-------------------------------------------------------------------------

      ptr_tl => spheres_on_levels(klev)
      nots = ptr_tl%no_triangles
      noes = ptr_tl%no_edges
      novs = ptr_tl%no_vertices

      ! 1) go through all coordinates and adjust them to a quadrilateral plane
      !-----------------------------------------------------------------------

      ! x-direction
      DO j = 0, nots-1
        IF( ptr_tl%ts(j)%triangle_center%x(1) >= 1.0_wp ) THEN
          ptr_tl%ts(j)%triangle_center%x(1) = &
            & ptr_tl%ts(j)%triangle_center%x(1) - 2.0_wp
        ENDIF
        IF( ptr_tl%ts(j)%triangle_center%x(1) < -1.0_wp ) THEN
          ptr_tl%ts(j)%triangle_center%x(1) = &
            & ptr_tl%ts(j)%triangle_center%x(1) + 2.0_wp
        ENDIF
      ENDDO
      DO j = 0, noes-1
        IF( ptr_tl%es(j)%edge_center%x(1) >= 1.0_wp ) THEN
          ptr_tl%es(j)%edge_center%x(1) = &
            & ptr_tl%es(j)%edge_center%x(1) - 2.0_wp
        ENDIF
        IF( ptr_tl%es(j)%edge_center%x(1) < -1.0_wp ) THEN
          ptr_tl%es(j)%edge_center%x(1) = &
            & ptr_tl%es(j)%edge_center%x(1) + 2.0_wp
        ENDIF
      ENDDO
      DO j = 0, novs-1
        IF( ptr_tl%vs(j)%vertex%x(1) >= 1.0_wp ) THEN
          ptr_tl%vs(j)%vertex%x(1) = &
            & ptr_tl%vs(j)%vertex%x(1) - 2.0_wp
        ENDIF
        IF( ptr_tl%vs(j)%vertex%x(1) < -1.0_wp ) THEN
          ptr_tl%vs(j)%vertex%x(1) = &
            & ptr_tl%vs(j)%vertex%x(1) + 2.0_wp
        ENDIF
      ENDDO

      ! y-direction
      DO j = 0, nots-1
        IF( ptr_tl%ts(j)%triangle_center%x(2) >= SQRT(3.0_wp)/2.0_wp ) THEN
          ptr_tl%ts(j)%triangle_center%x(2) = -SQRT(3.0_wp)/2.0_wp
        ENDIF
        ptr_tl%ts(j)%triangle_center%x(2) = MAX(-SQRT(3.0_wp)/2.0_wp , &
          & ptr_tl%ts(j)%triangle_center%x(2))
      ENDDO
      DO j = 0, noes-1
        IF( ptr_tl%es(j)%edge_center%x(2) >= SQRT(3.0_wp)/2.0_wp ) THEN
          ptr_tl%es(j)%edge_center%x(2) = -SQRT(3.0_wp)/2.0_wp
        ENDIF
        ptr_tl%es(j)%edge_center%x(2) = MAX(-SQRT(3.0_wp)/2.0_wp , &
          & ptr_tl%es(j)%edge_center%x(2))
      ENDDO
      DO j = 0, novs-1
        IF( ptr_tl%vs(j)%vertex%x(2) >= SQRT(3.0_wp)/2.0_wp ) THEN
          ptr_tl%vs(j)%vertex%x(2) = -SQRT(3.0_wp)/2.0_wp
        ENDIF
        ptr_tl%vs(j)%vertex%x(2) = MAX(-SQRT(3.0_wp)/2.0_wp , &
          & ptr_tl%vs(j)%vertex%x(2))
      ENDDO

      !2) All lengths and areas are the same
      !-------------------------------------

      !  z_primal_length = 1.0_wp/(kroot*(2**klev)) !:dimensionless
      z_primal_length = 1000.0_wp*tria_arc_km*(2.0_wp**(number_of_levels-klev))
      DO j = 0, noes-1
        ptr_tl%es(j)%edge_primal_arc = z_primal_length
        ptr_tl%es(j)%edge_dual_arc   = z_primal_length/SQRT(3.0_wp)
        ptr_tl%es(j)%edge_vert0_arc  = 0.5_wp*z_primal_length
        ptr_tl%es(j)%edge_vert1_arc  = 0.5_wp*z_primal_length
        ptr_tl%es(j)%edge_cell0_arc  = 0.5_wp*z_primal_length/SQRT(3.0_wp)
        ptr_tl%es(j)%edge_cell1_arc  = 0.5_wp*z_primal_length/SQRT(3.0_wp)
      ENDDO
      DO j = 0, nots-1
        ptr_tl%ts(j)%triangle_area   = z_primal_length*z_primal_length &
          & *sqrt(3.0_wp)*0.25_wp
      ENDDO
      DO j = 0, novs-1
        ptr_tl%vs(j)%area_dual_cell  = z_primal_length*z_primal_length &
          & *sqrt(3.0_wp)*0.5_wp
      ENDDO

    END SUBROUTINE compute_planar_geometry
    !-------------------------------------------------------------------------
    !
    !
    !>
    !!               Computes geometric quantities on grid level klev.
    !!
    !!
    !! @par Revision History
    !!  Luis Kornblueh, MPI-M, Hamburg, March 2004
    !!  Changes and corrections by Luca Bonaventura, MPI-M, August 2005
    !!  Cleanup by Luca Bonaventura, Polimi, February 2006
    !!  Modification by Th. Heinze, DWD, (2006-09-22):
    !!  - changed calculation of edge_primal_normal
    !!  Modification by P. Ripodas, DWD, (2007-02):
    !!  - The edge%system_orientation is calculated
    !!  Modification by Almut Gassmann, MPI-M, (2007-03-23)
    !!  - change for optional small circle subdivision
    !!  Modification by Almut Gassmann, MPI-M, (2007-04-17)
    !!  - added computation of area_edge and bug correction in computation of small
    !!    circle related triangle area
    !!  Modification fy Almut Gassmann, MPI-M, (2008-04-23)
    !!  - areas are no longer spherical by local "r**2 cos-phi-dlambda * dphi"
    !!    patches for the edge areas. Primal and dual areas are derived from them.
    !!
    SUBROUTINE compute_geometry (klev)
      !

      INTEGER, INTENT(in) :: klev

      INTEGER :: j,j1,j2, nots, noeds, novs, ie

      REAL(wp) :: z_tmp

      TYPE(t_cartesian_coordinates)   :: xt,yt,ctmp, c_ec, c_vd, c_nv, c_vp
      TYPE(t_geographical_coordinates):: geog

      TYPE(t_spheres), POINTER :: ptr_tl

#if (defined (_CRAYFTN) || defined(__INTEL_COMPILER) || defined(__SX__))
      TYPE(t_cartesian_coordinates), DIMENSION (0:spheres_on_levels(klev)%no_triangles-1) :: &
        &  tc, n0, n1, n2, v0, v1, v2, e_out
#endif

      !-------------------------------------------------------------------------

      ! As a reminder: the numbering scheme for vertices and edges
      ! is defined by the following rule:
      !
      !          vertex 0
      !       /\
      !    neighbor 0    /  \
      !   edge 0  /    \ edge 2
      !          /      \   neighbor 2
      !         /        \
      !     vertex 1 ------------ vertex 2
      !     edge 1
      !      neighbor 1
      !

      ptr_tl => spheres_on_levels(klev)
      nots = ptr_tl%no_triangles
      noeds= ptr_tl%no_edges
      novs = ptr_tl%no_vertices

      IF ((itype_optimize <=1 .or. itype_optimize >=4).and. &
        & (.not. (l_c_grid .and. klev==number_of_levels))) THEN

        !
        ! The following values are already known if the small circle computations
        ! were done.
        !

        ! triangle_center, edge_primal_arc
        !----------------
!$OMP PARALLEL

!$OMP DO
        DO j = 0, nots-1
          ptr_tl%ts(j)%triangle_center =  circum_center( &
            & ptr_tl%ts(j)%vertex0%vertex, &
            & ptr_tl%ts(j)%vertex1%vertex, &
            & ptr_tl%ts(j)%vertex2%vertex)
        ENDDO
!$OMP END DO

!$OMP DO
        DO j = 0, nots-1

          ptr_tl%ts(j)%edge0%edge_primal_arc = arc_length ( &
            & ptr_tl%ts(j)%vertex0%vertex, &
            & ptr_tl%ts(j)%vertex1%vertex)

          ptr_tl%ts(j)%edge1%edge_primal_arc = arc_length ( &
            & ptr_tl%ts(j)%vertex1%vertex, &
            & ptr_tl%ts(j)%vertex2%vertex)

          ptr_tl%ts(j)%edge2%edge_primal_arc = arc_length ( &
            & ptr_tl%ts(j)%vertex2%vertex, &
            & ptr_tl%ts(j)%vertex0%vertex)

        ENDDO
!$OMP END DO


        ! edge_center
        !--------------------------------------------

#if ( defined (_CRAYFTN) || defined(__SX__) || defined(__INTEL_COMPILER) )
! Remark (GZ): the unvectorized version causes a segfault with the Cray compiler
!$OMP DO
        DO j = 0, nots-1
          tc(j) = ptr_tl%ts(j)%triangle_center
          n0(j) = ptr_tl%ts(j)%neighbor0%triangle_center
          n1(j) = ptr_tl%ts(j)%neighbor1%triangle_center
          n2(j) = ptr_tl%ts(j)%neighbor2%triangle_center
          v0(j) = ptr_tl%ts(j)%vertex0%vertex
          v1(j) = ptr_tl%ts(j)%vertex1%vertex
          v2(j) = ptr_tl%ts(j)%vertex2%vertex
        ENDDO
!$OMP END DO
!$OMP END PARALLEL

        CALL spherical_intersection2_v (tc,n0,v0,v1,e_out,nots)
        DO j = 0, nots-1
          ptr_tl%ts(j)%edge0%edge_center = e_out(j)
        ENDDO

        CALL spherical_intersection2_v (tc,n1,v1,v2,e_out,nots)
        DO j = 0, nots-1
          ptr_tl%ts(j)%edge1%edge_center = e_out(j)
        ENDDO

        CALL spherical_intersection2_v (tc,n2,v2,v0,e_out,nots)
        DO j = 0, nots-1
          ptr_tl%ts(j)%edge2%edge_center = e_out(j)
        ENDDO

#else

!$OMP END PARALLEL

!GZ: Due to the indirectly addressed target field, OpenMP parallelization of this loop would cause
!    a race condition. Perhaps we should use the vectorized version above for all platfoms?
        DO j = 0, nots-1

          ptr_tl%ts(j)%edge0%edge_center     = spherical_intersection2 ( &
            & ptr_tl%ts(j)%triangle_center, &
            & ptr_tl%ts(j)%neighbor0%triangle_center, &
            & ptr_tl%ts(j)%vertex0%vertex, &
            & ptr_tl%ts(j)%vertex1%vertex)

          ptr_tl%ts(j)%edge1%edge_center     = spherical_intersection2 ( &
            & ptr_tl%ts(j)%triangle_center, &
            & ptr_tl%ts(j)%neighbor1%triangle_center, &
            & ptr_tl%ts(j)%vertex1%vertex,  &
            & ptr_tl%ts(j)%vertex2%vertex)

          ptr_tl%ts(j)%edge2%edge_center     = spherical_intersection2 ( &
            & ptr_tl%ts(j)%triangle_center, &
            & ptr_tl%ts(j)%neighbor2%triangle_center, &
            & ptr_tl%ts(j)%vertex2%vertex, &
            & ptr_tl%ts(j)%vertex0%vertex)

        ENDDO
#endif
      ENDIF

      !
      ! The dual arc lengths are unknown so far for all itype_optimize
      !
!$OMP PARALLEL               
!$OMP DO
      DO j = 0, nots-1

        ptr_tl%ts(j)%edge0%edge_dual_arc   = arc_length ( &
          & ptr_tl%ts(j)%triangle_center, &
          & ptr_tl%ts(j)%neighbor0%triangle_center)

        ptr_tl%ts(j)%edge1%edge_dual_arc   = arc_length ( &
          & ptr_tl%ts(j)%triangle_center, &
          & ptr_tl%ts(j)%neighbor1%triangle_center)

        ptr_tl%ts(j)%edge2%edge_dual_arc   = arc_length ( &
          & ptr_tl%ts(j)%triangle_center, &
          & ptr_tl%ts(j)%neighbor2%triangle_center)

      ENDDO
!$OMP END DO

      !
      ! distances to edge midpoints
      !
!$OMP DO
      DO j = 0, noeds-1

        ptr_tl%es(j)%edge_vert0_arc = 0.5_wp*ptr_tl%es(j)%edge_primal_arc
        ptr_tl%es(j)%edge_vert1_arc = 0.5_wp*ptr_tl%es(j)%edge_primal_arc

        ptr_tl%es(j)%edge_cell0_arc = arc_length ( &
          & ptr_tl%es(j)%edge_center,&
          & ptr_tl%es(j)%triangle0%triangle_center)
        ptr_tl%es(j)%edge_cell1_arc = arc_length ( &
          & ptr_tl%es(j)%edge_center,&
          & ptr_tl%es(j)%triangle1%triangle_center)

      ENDDO
!$OMP END DO

      ! This part  copies exactly what had been done in the previous grid generator.
      ! The normal vectors to primal and dual cell edges are oriented in a way
      ! consistent with the values stored in the index arrays e.g., the normal vector
      ! to the primal cell edge points from edge%triangle0 to edge%triangle1
      ! (and later, in the 'grid' data structure, from cell%edge_index(.,1) to
      ! cell%edge_index(.,2). The previous computation did not guarantee this.
      ! What is stored in   ptr_tl%es(j)%edge_primal_normal%v1 and
      ! ptr_tl%es(j)%edge_primal_normal%v2 are the scalar products of the normal
      ! vectors with the i,j unit vectors on the tangent plane (e.g. i,j in the
      ! notation of Holton's book or Williamson's shallow water test paper. These
      ! quantities are used to compute the normal velocity field from the u,v
      ! components at each cell edge

      ! changed by Th.Heinze, now normal vector is tangential to the sphere in the
      ! midpoint of the edges, was not in previous version.

!$OMP DO PRIVATE(j1,j2,ctmp,z_tmp,c_vd,c_ec,c_nv,geog,xt,yt,c_vp)
      DO j = 0, noeds-1

        ! determine vector between mass points (needed for orientation)

        j1 = ptr_tl%es(j)%triangle0%info%idx
        j2 = ptr_tl%es(j)%triangle1%info%idx

        ctmp%x(:) = ptr_tl%ts(j2)%triangle_center%x(:) &
          & - ptr_tl%ts(j1)%triangle_center%x(:)

        z_tmp = SQRT(DOT_PRODUCT(ctmp%x,ctmp%x))

        ctmp%x = ctmp%x/z_tmp

        ! determine vector between vorticity points

        c_vd%x = ptr_tl%es(j)%vertex1%vertex%x &
          & - ptr_tl%es(j)%vertex0%vertex%x

        ! determine position vector in center of edge (outward unit vector to sphere)

        c_ec  = ptr_tl%es(j)%edge_center

        ! determine normal vector in center of edge as vector product of c_vd and c_ec

        c_nv = vector_product(c_vd, c_ec)

        ! normalize normal vector c_nv

        z_tmp  = SQRT(DOT_PRODUCT(c_nv%x,c_nv%x))

        c_nv%x = c_nv%x/z_tmp

        ! determine orientation of normal vector

        z_tmp = DOT_PRODUCT(c_nv%x,ctmp%x)

        IF (z_tmp <0._wp) c_nv%x = -1._wp * c_nv%x

        ! calculate geographical coordinates of edge_center

        geog = cc2gc(ptr_tl%es(j)%edge_center)

        ! changed calculation by thh: removed factor COS(geog%lat) in first two comp.

        xt%x(1) = -SIN(geog%lon)
        xt%x(2) =  COS(geog%lon)
        xt%x(3) =  0._wp

        ptr_tl%es(j)%edge_primal_normal%v1 = DOT_PRODUCT(c_nv%x,xt%x)

        yt%x(1) = -COS(geog%lon)*sin(geog%lat)
        yt%x(2) = -SIN(geog%lon)*sin(geog%lat)
        yt%x(3) =  COS(geog%lat)

        ptr_tl%es(j)%edge_primal_normal%v2=DOT_PRODUCT(c_nv%x,yt%x)

        ptr_tl%es(j)%edge_dual_normal%v1 =  ptr_tl%es(j)%edge_primal_normal%v2
        ptr_tl%es(j)%edge_dual_normal%v2 = -ptr_tl%es(j)%edge_primal_normal%v1

        c_vp = vector_product(c_vd,ctmp)
        ptr_tl%es(j)%system_orientation = INT(SIGN(1._wp,DOT_PRODUCT(c_vp%x,c_ec%x)))

      ENDDO
!$OMP END DO

      ! Compute triangle_area and area_dual_cell. These are local
      ! computations using the edge lengths and thus no longer spherical
      ! approximations.

!$OMP DO PRIVATE(j1)
      DO j=0, nots-1

        ptr_tl%ts(j)%triangle_area=0.0_wp

        j1 = ptr_tl%ts(j)%edge0%triangle0%info%idx
        IF(j1==j)THEN
          ptr_tl%ts(j)%triangle_area =  ptr_tl%ts(j)%triangle_area +&
            & 0.5_wp*ptr_tl%ts(j)%edge0%edge_primal_arc* &
            & ptr_tl%ts(j)%edge0%edge_cell0_arc
        ELSE
          ptr_tl%ts(j)%triangle_area =  ptr_tl%ts(j)%triangle_area +&
            & 0.5_wp*ptr_tl%ts(j)%edge0%edge_primal_arc* &
            & ptr_tl%ts(j)%edge0%edge_cell1_arc
        ENDIF
        j1 = ptr_tl%ts(j)%edge1%triangle0%info%idx
        IF(j1==j)THEN
          ptr_tl%ts(j)%triangle_area =  ptr_tl%ts(j)%triangle_area +&
            & 0.5_wp*ptr_tl%ts(j)%edge1%edge_primal_arc* &
            & ptr_tl%ts(j)%edge1%edge_cell0_arc
        ELSE
          ptr_tl%ts(j)%triangle_area =  ptr_tl%ts(j)%triangle_area +&
            & 0.5_wp*ptr_tl%ts(j)%edge1%edge_primal_arc* &
            & ptr_tl%ts(j)%edge1%edge_cell1_arc
        ENDIF
        j1 = ptr_tl%ts(j)%edge2%triangle0%info%idx
        IF(j1==j)THEN
          ptr_tl%ts(j)%triangle_area =  ptr_tl%ts(j)%triangle_area +&
            & 0.5_wp*ptr_tl%ts(j)%edge2%edge_primal_arc* &
            & ptr_tl%ts(j)%edge2%edge_cell0_arc
        ELSE
          ptr_tl%ts(j)%triangle_area =  ptr_tl%ts(j)%triangle_area +&
            & 0.5_wp*ptr_tl%ts(j)%edge2%edge_primal_arc* &
            & ptr_tl%ts(j)%edge2%edge_cell1_arc
        ENDIF

      ENDDO
!$OMP END DO

!$OMP DO
      DO j=0, novs-1

        ptr_tl%vs(j)%area_dual_cell = 0.25_wp*(&
          & ptr_tl%vs(j)%edge0%edge_primal_arc*ptr_tl%vs(j)%edge0%edge_dual_arc + &
          & ptr_tl%vs(j)%edge1%edge_primal_arc*ptr_tl%vs(j)%edge1%edge_dual_arc + &
          & ptr_tl%vs(j)%edge2%edge_primal_arc*ptr_tl%vs(j)%edge2%edge_dual_arc + &
          & ptr_tl%vs(j)%edge3%edge_primal_arc*ptr_tl%vs(j)%edge3%edge_dual_arc + &
          & ptr_tl%vs(j)%edge4%edge_primal_arc*ptr_tl%vs(j)%edge4%edge_dual_arc )

        IF (ptr_tl%vs(j)%edge5%info%idx /= -1) THEN

          ptr_tl%vs(j)%area_dual_cell = ptr_tl%vs(j)%area_dual_cell+ 0.25_wp*&
            & ptr_tl%vs(j)%edge5%edge_primal_arc*ptr_tl%vs(j)%edge5%edge_dual_arc

        ENDIF

      ENDDO
!$OMP END DO

      !ALMUT: moved from mo_grid (2009-01-13)
      !   compute primal cell edge orientations:
      !   computation is done by an algebraic formula
      !   (see Casulli and Walters, IJNMF, 1998)
      !   which gives the correct sign
      ! IF ONE ASSUMES THAT THE POSITIVE DIRECTION
      ! FOR THE NORMAL VELOCITY COMPONENT
      ! IS THAT GOING FROM cg%cells%edge_index(j,2)
      ! TO cg%cells%edge_index(j,1), as it is
      !  guaranteed by the previous definitions
      !  (it was so in the old generator)
      ! ALMUT: this explanation seems to mix up the direction, but the code is ok

!$OMP DO PRIVATE(ie)
      DO j=0,nots-1
        ie= ptr_tl%ts(j)%edge0%info%idx
        ptr_tl%ts(j)%edge0_orientation=&
          & (ptr_tl%es(ie)%triangle1%info%idx+ptr_tl%es(ie)%triangle0%info%idx-2*j) &
          & /(ptr_tl%es(ie)%triangle1%info%idx-ptr_tl%es(ie)%triangle0%info%idx)
        ie= ptr_tl%ts(j)%edge1%info%idx
        ptr_tl%ts(j)%edge1_orientation=&
          & (ptr_tl%es(ie)%triangle1%info%idx+ptr_tl%es(ie)%triangle0%info%idx-2*j) &
          & /(ptr_tl%es(ie)%triangle1%info%idx-ptr_tl%es(ie)%triangle0%info%idx)
        ie= ptr_tl%ts(j)%edge2%info%idx
        ptr_tl%ts(j)%edge2_orientation=&
          & (ptr_tl%es(ie)%triangle1%info%idx+ptr_tl%es(ie)%triangle0%info%idx-2*j) &
          & /(ptr_tl%es(ie)%triangle1%info%idx-ptr_tl%es(ie)%triangle0%info%idx)
      ENDDO
!$OMP END DO

      ! Almut (2009-01-13): moved from input_grid in mo_io_grid at this position
      !   rescale geometric quantities by radius of the Earth
!$OMP DO
      DO j=0, novs-1
        ptr_tl%vs(j)%area_dual_cell   = earth_radius*earth_radius*ptr_tl%vs(j)%area_dual_cell
      ENDDO
!$OMP END DO
!$OMP DO
      DO j=0,nots-1
        ptr_tl%ts(j)%triangle_area    = earth_radius*earth_radius*ptr_tl%ts(j)%triangle_area
      ENDDO
!$OMP END DO
!$OMP DO
      DO j=0, noeds-1
        ptr_tl%es(j)%edge_dual_arc    = earth_radius*   ptr_tl%es(j)%edge_dual_arc
        ptr_tl%es(j)%edge_primal_arc  = earth_radius*   ptr_tl%es(j)%edge_primal_arc
        ptr_tl%es(j)%edge_vert0_arc   = earth_radius*   ptr_tl%es(j)%edge_vert0_arc
        ptr_tl%es(j)%edge_vert1_arc   = earth_radius*   ptr_tl%es(j)%edge_vert1_arc
        ptr_tl%es(j)%edge_cell0_arc   = earth_radius*   ptr_tl%es(j)%edge_cell0_arc
        ptr_tl%es(j)%edge_cell1_arc   = earth_radius*   ptr_tl%es(j)%edge_cell1_arc
      ENDDO
!$OMP END DO NOWAIT

!$OMP END PARALLEL

    END SUBROUTINE compute_geometry


    !------------------------------------------------------------------------------

  END MODULE mo_geometry

