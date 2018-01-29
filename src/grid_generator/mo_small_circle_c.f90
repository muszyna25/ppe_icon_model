!>
!! This module provides all subroutines necessary for the small circle.
!!
!! This module provides all subroutines necessary for the small circle
!! subdivision. Small circle subdivision is an alternative to grid optimization
!! and defines small circles as edges of the triangles such that the distances
!! between the edge midpoint and the adjacent triangle centers are equal. Thus,
!! second order accuracy may be achieved for the gradient operator.
!! The boundary of dual grid cells are formed by great circles. Consequently, the
!! Voronoi-Delaunay properties are kept when using this type of subdivision.
!!
!! This module contains two main subroutines:
!! <ul>
!! <li>  {\\tt scc_inner:}
!!       Computes the subdivision vertices and inner small circle
!!       properties (co-latitude and pole of the small circle).
!!       The small circle properties of the outer edges of the
!!       coarse triangle are copied to the edges of the refined
!!       triangles that form the coarse boundary. These values
!!       may be overwritten by the subroutine {\\tt scc_outer}.
!! <li>  {\\tt scc_outer:}
!!       Computes the small circle properties of the edges of the
!!       refined triangles forming the boundary of the coarse triangle.
!! </ul>
!!
!! @par Revision History
!! Initial release by Almut Gassmann, MPI-M (2007-03-23)
!! Modification by Almut Gassmann, MPI-M (2007-04-17)
!! - imported mo_exception
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
MODULE mo_small_circle_c
  !-------------------------------------------------------------------------
  !
  !    ProTeX FORTRAN source: Style 2
  !    modified for ICON project, DWD/MPI-M 2006
  !
  !-------------------------------------------------------------------------
  !
  !
  !
  USE mo_kind,          ONLY: wp
  USE mo_math_types,    ONLY: t_cartesian_coordinates
  USE mo_math_utilities, ONLY: circum_center, arc_length
  USE mo_base_datatypes,ONLY: t_triangle, t_spheres
  USE mo_topology,      ONLY: root, kroot, icosahedron
  USE mo_equal_area,    ONLY: edge_midpoint
  USE mo_exception,     ONLY: finish

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: scc_inner, scc_outer, scc_prepare, scc_inner_root

CONTAINS

  !-------------------------------------------------------------------------
  !-------------------------------------------------------------------------
  !
  !
  !>
  !! This subroutine is thought for c-grid optimization of the inner parts.
  !!
  !! This subroutine is thought for c-grid optimization of the inner parts
  !! of a four fold triangle.
  !!
  !! @par Revision History
  !! Initial release by Almut Gassmann (2007-03)
  !!
  SUBROUTINE scc_inner(cn)
    !
    TYPE(t_triangle), POINTER :: cn
    !
    TYPE(t_triangle), POINTER :: s0, s1, s2, s3
    REAL(wp) :: z_crcl, z_scal
    TYPE(t_cartesian_coordinates) :: x_m0,  x_m1,  x_m2,       &
      & x_mm0, x_mm1, x_mm2,      &
      & x_c0,  x_c1,  x_c2, x_c3, &
      & x_pole
    !
    !-------------------------------------------------------------------------
    !-----------------------
    ! Section 0: preparation
    !-----------------------

    s0  => cn%sub_triangle0
    s1  => cn%sub_triangle1
    s2  => cn%sub_triangle2
    s3  => cn%sub_triangle3
    x_m0 = cn%edge0%edge_center
    x_m1 = cn%edge1%edge_center
    x_m2 = cn%edge2%edge_center

    !---------------------------------------------
    ! Section 1: compute new vertices on the edges
    !---------------------------------------------

    s0%vertex1%vertex = x_m0
    s1%vertex0%vertex = x_m0
    s2%vertex2%vertex = x_m0

    s1%vertex2%vertex = x_m1
    s2%vertex0%vertex = x_m1
    s3%vertex1%vertex = x_m1

    s0%vertex2%vertex = x_m2
    s2%vertex1%vertex = x_m2
    s3%vertex0%vertex = x_m2

    s0%vertex0%vertex = cn%vertex0%vertex
    s1%vertex1%vertex = cn%vertex1%vertex
    s3%vertex2%vertex = cn%vertex2%vertex

    !-------------------------------------------
    ! Section 2: compute centers of subtriangles
    !-------------------------------------------

    x_c0 = circum_center(cn%vertex0%vertex,x_m0,x_m2)
    x_c1 = circum_center(x_m0,cn%vertex1%vertex,x_m1)
    x_c2 = circum_center(x_m1,x_m2,x_m0)
    x_c3 = circum_center(x_m2,x_m1,cn%vertex2%vertex)

    s0%triangle_center = x_c0
    s1%triangle_center = x_c1
    s2%triangle_center = x_c2
    s3%triangle_center = x_c3

    !----------------------------------------------------------------------
    ! Section 3: compute edge midpoints of subtriangles as midpoints of the
    !            cell centers on great arcs
    !----------------------------------------------------------------------

    ! "small circle" properties of a great circle
    x_pole%x(:)=0.0_wp
    z_crcl     =0.0_wp
    z_scal     =0.0_wp ! just set it for it is defined at the call of the subroutine

    x_mm0 = edge_midpoint(x_c3,x_c2,z_crcl,x_pole)
    x_mm1 = edge_midpoint(x_c0,x_c2,z_crcl,x_pole)
    x_mm2 = edge_midpoint(x_c1,x_c2,z_crcl,x_pole)

    s0%edge1%edge_center = x_mm1
    s1%edge2%edge_center = x_mm2
    s3%edge0%edge_center = x_mm0

    s2%edge0%edge_center = x_mm0
    s2%edge1%edge_center = x_mm1
    s2%edge2%edge_center = x_mm2

    !--------------------------------------------------------------
    ! Section 4: compute pole and cos_rel_co_lat of the inner edges
    !--------------------------------------------------------------

    CALL small_circle_features(x_m0,x_m2,x_mm1,z_crcl,x_pole,z_scal)
    s0%edge1%cos_rel_co_lat = z_crcl
    s0%edge1%sc_pole        = x_pole
    s0%edge1%edge_primal_arc= z_scal
    s2%edge1%cos_rel_co_lat = z_crcl
    s2%edge1%sc_pole        = x_pole
    s2%edge1%edge_primal_arc= z_scal

    CALL small_circle_features(x_m0,x_m1,x_mm2,z_crcl,x_pole,z_scal)
    s1%edge2%cos_rel_co_lat = z_crcl
    s1%edge2%sc_pole        = x_pole
    s1%edge2%edge_primal_arc= z_scal
    s2%edge2%cos_rel_co_lat = z_crcl
    s2%edge2%sc_pole        = x_pole
    s2%edge2%edge_primal_arc= z_scal


    CALL small_circle_features(x_m1,x_m2,x_mm0,z_crcl,x_pole,z_scal)
    s3%edge0%cos_rel_co_lat = z_crcl
    s3%edge0%sc_pole        = x_pole
    s3%edge0%edge_primal_arc= z_scal
    s2%edge0%cos_rel_co_lat = z_crcl
    s2%edge0%sc_pole        = x_pole
    s2%edge0%edge_primal_arc= z_scal

    !-------------------------------------------------------------------
    ! Section 5: copy the small circle properties from the outer edge to
    !            the subtriangle's outer edges
    !-------------------------------------------------------------------

    s0%edge0%cos_rel_co_lat = cn%edge0%cos_rel_co_lat
    s0%edge2%cos_rel_co_lat = cn%edge2%cos_rel_co_lat
    s1%edge0%cos_rel_co_lat = cn%edge0%cos_rel_co_lat
    s1%edge1%cos_rel_co_lat = cn%edge1%cos_rel_co_lat
    s3%edge1%cos_rel_co_lat = cn%edge1%cos_rel_co_lat
    s3%edge2%cos_rel_co_lat = cn%edge2%cos_rel_co_lat

    s0%edge0%sc_pole = cn%edge0%sc_pole
    s0%edge2%sc_pole = cn%edge2%sc_pole
    s1%edge0%sc_pole = cn%edge0%sc_pole
    s1%edge1%sc_pole = cn%edge1%sc_pole
    s3%edge1%sc_pole = cn%edge1%sc_pole
    s3%edge2%sc_pole = cn%edge2%sc_pole

    s0%edge0%edge_primal_arc = 0.5_wp * cn%edge0%edge_primal_arc
    s0%edge2%edge_primal_arc = 0.5_wp * cn%edge2%edge_primal_arc
    s1%edge0%edge_primal_arc = 0.5_wp * cn%edge0%edge_primal_arc
    s1%edge1%edge_primal_arc = 0.5_wp * cn%edge1%edge_primal_arc
    s3%edge1%edge_primal_arc = 0.5_wp * cn%edge1%edge_primal_arc
    s3%edge2%edge_primal_arc = 0.5_wp * cn%edge2%edge_primal_arc

    !-------------------------------------------------------------------------
    ! Section 6: compute the edge midpoints of outer boundary for completeness
    !-------------------------------------------------------------------------

    x_pole = cn%edge0%sc_pole
    z_crcl = cn%edge0%cos_rel_co_lat
    s0%edge0%edge_center = edge_midpoint(cn%vertex0%vertex,x_m0,z_crcl,x_pole)
    s1%edge0%edge_center = edge_midpoint(cn%vertex1%vertex,x_m0,z_crcl,x_pole)

    x_pole = cn%edge1%sc_pole
    z_crcl = cn%edge1%cos_rel_co_lat
    s1%edge1%edge_center = edge_midpoint(cn%vertex1%vertex,x_m1,z_crcl,x_pole)
    s3%edge1%edge_center = edge_midpoint(cn%vertex2%vertex,x_m1,z_crcl,x_pole)

    x_pole = cn%edge2%sc_pole
    z_crcl = cn%edge2%cos_rel_co_lat
    s3%edge2%edge_center = edge_midpoint(cn%vertex2%vertex,x_m2,z_crcl,x_pole)
    s0%edge2%edge_center = edge_midpoint(cn%vertex0%vertex,x_m2,z_crcl,x_pole)

  END SUBROUTINE scc_inner
  !
  !-------------------------------------------------------------------------
  !-------------------------------------------------------------------------
  !
  !
  !>
  !! This subroutine computes the default small circle features on the root level.
  !!
  !!
  !! @par Revision History
  !! Initial release by Almut Gassmann (2007-03)
  !!
  SUBROUTINE scc_inner_root(jdx)
    !
    INTEGER, INTENT(in) :: jdx
    !
    INTEGER :: ji, j, ict
    REAL(wp) :: z_crcl
    TYPE(t_cartesian_coordinates) :: x_pole
    !
    !-------------------------------------------------------------------------

    DO ji=0,kroot-1
      DO j=0,2*ji
        ict = ji*ji+j
        root(ict,jdx)%t%triangle_center=&
          & circum_center(root(ict,jdx)%t%vertex0%vertex, &
          & root(ict,jdx)%t%vertex1%vertex, &
          & root(ict,jdx)%t%vertex2%vertex)
      ENDDO
    ENDDO

    DO ji=0,kroot-1
      DO j=0,2*ji
        ict = ji*ji+j
        IF (MODULO(j+1,2)==0) CYCLE
        IF (ji/=kroot-1) THEN
          ! lower edge of inner triangles
          x_pole%x(:) = 0.0_wp
          z_crcl      = 0.0_wp
          root(ict,jdx)%t%edge1%edge_center=edge_midpoint(&
            & root(ict,jdx)%t%triangle_center, &
            & root(ict,jdx)%t%neighbor1%triangle_center,z_crcl,x_pole)
          CALL small_circle_features(root(ict,jdx)%t%vertex1%vertex,&
            & root(ict,jdx)%t%vertex2%vertex,&
            & root(ict,jdx)%t%edge1%edge_center,&
            & root(ict,jdx)%t%edge1%cos_rel_co_lat,&
            & root(ict,jdx)%t%edge1%sc_pole,&
            & root(ict,jdx)%t%edge1%edge_primal_arc)
        ELSE ! lower boundary, these edges stay great circles
          x_pole%x(:) = 0.0_wp
          z_crcl      = 0.0_wp
          root(ict,jdx)%t%edge1%edge_center    = edge_midpoint(&
            & root(ict,jdx)%t%vertex1%vertex, &
            & root(ict,jdx)%t%vertex2%vertex,z_crcl,x_pole)
          root(ict,jdx)%t%edge1%cos_rel_co_lat = 0.0_wp
          root(ict,jdx)%t%edge1%sc_pole%x(:)   = 0.0_wp
          root(ict,jdx)%t%edge1%edge_primal_arc=&
            & icosahedron(jdx)%t%edge1%edge_primal_arc/(REAL(kroot,wp))
        ENDIF
        IF (j/=0) THEN
          ! left edge of inner triangles
          x_pole%x(:) = 0.0_wp
          z_crcl      = 0.0_wp
          root(ict,jdx)%t%edge0%edge_center=edge_midpoint(&
            & root(ict,jdx)%t%triangle_center, &
            & root(ict,jdx)%t%neighbor0%triangle_center,z_crcl,x_pole)
          CALL small_circle_features(root(ict,jdx)%t%vertex0%vertex,&
            & root(ict,jdx)%t%vertex1%vertex,&
            & root(ict,jdx)%t%edge0%edge_center,&
            & root(ict,jdx)%t%edge0%cos_rel_co_lat,&
            & root(ict,jdx)%t%edge0%sc_pole,&
            & root(ict,jdx)%t%edge0%edge_primal_arc)
        ELSE ! left boundary, these edges stay great circles
          x_pole%x(:) = 0.0_wp
          z_crcl      = 0.0_wp
          root(ict,jdx)%t%edge0%edge_center    = edge_midpoint(&
            & root(ict,jdx)%t%vertex0%vertex, &
            & root(ict,jdx)%t%vertex1%vertex,z_crcl,x_pole)
          root(ict,jdx)%t%edge0%cos_rel_co_lat = 0.0_wp
          root(ict,jdx)%t%edge0%sc_pole%x(:)   = 0.0_wp
          root(ict,jdx)%t%edge0%edge_primal_arc=&
            & icosahedron(jdx)%t%edge0%edge_primal_arc/(REAL(kroot,wp))
        ENDIF
        IF (j/=2*ji) THEN
          ! right edge of inner triangles
          x_pole%x(:) = 0.0_wp
          z_crcl      = 0.0_wp
          root(ict,jdx)%t%edge2%edge_center=edge_midpoint(&
            & root(ict,jdx)%t%triangle_center, &
            & root(ict,jdx)%t%neighbor2%triangle_center,z_crcl,x_pole)
          CALL small_circle_features(root(ict,jdx)%t%vertex0%vertex,&
            & root(ict,jdx)%t%vertex2%vertex,&
            & root(ict,jdx)%t%edge2%edge_center,&
            & root(ict,jdx)%t%edge2%cos_rel_co_lat,&
            & root(ict,jdx)%t%edge2%sc_pole,&
            & root(ict,jdx)%t%edge2%edge_primal_arc)
        ELSE ! right boundary, these edges stay great circles
          x_pole%x(:) = 0.0_wp
          z_crcl      = 0.0_wp
          root(ict,jdx)%t%edge2%edge_center    = edge_midpoint(&
            & root(ict,jdx)%t%vertex2%vertex, &
            & root(ict,jdx)%t%vertex0%vertex,z_crcl,x_pole)
          root(ict,jdx)%t%edge2%cos_rel_co_lat = 0.0_wp
          root(ict,jdx)%t%edge2%sc_pole%x(:)   = 0.0_wp
          root(ict,jdx)%t%edge2%edge_primal_arc=&
            & icosahedron(jdx)%t%edge2%edge_primal_arc/(REAL(kroot,wp))
        ENDIF
      ENDDO
    ENDDO

  END SUBROUTINE scc_inner_root
  !
  !-------------------------------------------------------------------------
  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Provides all small circle features on the before last level if the option.
  !!
  !! Provides all small circle features on the before last level if the option
  !! l_c_grid=.TRUE. is chosen on the unoptimized or Heikes optimized grid.
  !!
  !! @par Revision History
  !! Initial release by Almut Gassmann (2007-03)
  !!
  SUBROUTINE scc_prepare (sph)
    !
    TYPE(t_spheres) :: sph
    !
    INTEGER :: i_nots, i_noes, j
    !
    !-------------------------------------------------------------------------
    !
    i_nots = sph%no_triangles
    i_noes = sph%no_edges
    DO j = 0, i_nots-1
      sph%ts(j)%triangle_center = circum_center(sph%ts(j)%vertex0%vertex,&
        & sph%ts(j)%vertex1%vertex,&
        & sph%ts(j)%vertex2%vertex)
    ENDDO

    DO j = 0, i_noes-1
      sph%es(j)%cos_rel_co_lat  = 0.0_wp
      sph%es(j)%sc_pole%x       = 0.0_wp
      sph%es(j)%edge_center     = edge_midpoint(sph%es(j)%vertex0%vertex,&
        & sph%es(j)%vertex1%vertex,&
        & sph%es(j)%cos_rel_co_lat,&
        & sph%es(j)%sc_pole )
      sph%es(j)%edge_primal_arc = arc_length(sph%es(j)%vertex0%vertex,&
        & sph%es(j)%vertex1%vertex)
    ENDDO

  END SUBROUTINE scc_prepare
  !
  !-------------------------------------------------------------------------
  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Here, the small circle features of the outer boundary edges of a four.
  !!
  !! Here, the small circle features of the outer boundary edges of a four
  !! fold triangle are computed.
  !!
  !! @par Revision History
  !! Initial release by Almut Gassmann (2007-03)
  !!
  SUBROUTINE scc_outer (sph_coarse)
    !
    TYPE (t_spheres) :: sph_coarse
    !
    INTEGER :: i_nots, j
    TYPE (t_cartesian_coordinates) :: v0, v1, c0, c1, xpole, em
    REAL (wp) :: z_crcl, z_scal
    !
    !-------------------------------------------------------------------------

    i_nots = sph_coarse%no_triangles

    !-------------------------------
    ! Loop over all coarse triangles
    !-------------------------------
    DO j = 0, i_nots-1

      !---------------------------------------------------
      ! compute both fine edges that form the coarse edge0
      !---------------------------------------------------

      IF(sph_coarse%ts(j)%neighbor0%info%idx > j) THEN

        ! first edge
        !-----------
        v0 = sph_coarse%ts(j)%vertex0%vertex
        v1 = sph_coarse%ts(j)%edge0%edge_center
        c0 = sph_coarse%ts(j)%sub_triangle0%triangle_center
        c1 = sph_coarse%ts(j)%sub_triangle0%neighbor0%triangle_center

        xpole%x(:)=0.0_wp
        z_crcl=0.0_wp
        z_scal=0.0_wp ! just for definition at the subroutine call
        em = edge_midpoint(c0,c1,z_crcl,xpole)
        sph_coarse%ts(j)%sub_triangle0%edge0%edge_center = em

        CALL small_circle_features(v0,v1,em,z_crcl,xpole,z_scal)
        sph_coarse%ts(j)%sub_triangle0%edge0%cos_rel_co_lat  = z_crcl
        sph_coarse%ts(j)%sub_triangle0%edge0%sc_pole         = xpole
        sph_coarse%ts(j)%sub_triangle0%edge0%edge_primal_arc = z_scal

        ! second edge
        !------------

        v0 = sph_coarse%ts(j)%edge0%edge_center
        v1 = sph_coarse%ts(j)%vertex1%vertex
        c0 = sph_coarse%ts(j)%sub_triangle1%triangle_center
        c1 = sph_coarse%ts(j)%sub_triangle1%neighbor0%triangle_center

        xpole%x(:)=0.0_wp
        z_crcl=0.0_wp
        z_scal=0.0_wp ! just for definition at the subroutine call
        em = edge_midpoint(c0,c1,z_crcl,xpole)
        sph_coarse%ts(j)%sub_triangle1%edge0%edge_center = em

        CALL small_circle_features(v0,v1,em,z_crcl,xpole,z_scal)
        sph_coarse%ts(j)%sub_triangle1%edge0%cos_rel_co_lat  = z_crcl
        sph_coarse%ts(j)%sub_triangle1%edge0%sc_pole         = xpole
        sph_coarse%ts(j)%sub_triangle1%edge0%edge_primal_arc = z_scal

      ENDIF

      !---------------------------------------------------
      ! compute both fine edges that form the coarse edge1
      !---------------------------------------------------

      IF(sph_coarse%ts(j)%neighbor1%info%idx > j) THEN

        ! first edge
        !-----------

        v0 = sph_coarse%ts(j)%vertex1%vertex
        v1 = sph_coarse%ts(j)%edge1%edge_center
        c0 = sph_coarse%ts(j)%sub_triangle1%triangle_center
        c1 = sph_coarse%ts(j)%sub_triangle1%neighbor1%triangle_center

        xpole%x(:)=0.0_wp
        z_crcl=0.0_wp
        z_scal=0.0_wp ! just for definition at the subroutine call
        em = edge_midpoint(c0,c1,z_crcl,xpole)
        sph_coarse%ts(j)%sub_triangle1%edge1%edge_center = em

        CALL small_circle_features(v0,v1,em,z_crcl,xpole,z_scal)
        sph_coarse%ts(j)%sub_triangle1%edge1%cos_rel_co_lat  = z_crcl
        sph_coarse%ts(j)%sub_triangle1%edge1%sc_pole         = xpole
        sph_coarse%ts(j)%sub_triangle1%edge1%edge_primal_arc = z_scal

        ! second edge
        !------------

        v0 = sph_coarse%ts(j)%edge1%edge_center
        v1 = sph_coarse%ts(j)%vertex2%vertex
        c0 = sph_coarse%ts(j)%sub_triangle3%triangle_center
        c1 = sph_coarse%ts(j)%sub_triangle3%neighbor1%triangle_center

        xpole%x(:)=0.0_wp
        z_crcl=0.0_wp
        z_scal=0.0_wp ! just for definition at the subroutine call
        em = edge_midpoint(c0,c1,z_crcl,xpole)
        sph_coarse%ts(j)%sub_triangle3%edge1%edge_center = em

        CALL small_circle_features(v0,v1,em,z_crcl,xpole,z_scal)
        sph_coarse%ts(j)%sub_triangle3%edge1%cos_rel_co_lat  = z_crcl
        sph_coarse%ts(j)%sub_triangle3%edge1%sc_pole         = xpole
        sph_coarse%ts(j)%sub_triangle3%edge1%edge_primal_arc = z_scal

      ENDIF

      !---------------------------------------------------
      ! compute both fine edges that form the coarse edge2
      !---------------------------------------------------

      IF(sph_coarse%ts(j)%neighbor2%info%idx > j) THEN

        ! first edge
        !-----------

        v0 = sph_coarse%ts(j)%vertex2%vertex
        v1 = sph_coarse%ts(j)%edge2%edge_center
        c0 = sph_coarse%ts(j)%sub_triangle3%triangle_center
        c1 = sph_coarse%ts(j)%sub_triangle3%neighbor2%triangle_center

        xpole%x(:)=0.0_wp
        z_crcl=0.0_wp
        z_scal=0.0_wp ! just for definition at the subroutine call
        em = edge_midpoint(c0,c1,z_crcl,xpole)
        sph_coarse%ts(j)%sub_triangle3%edge2%edge_center = em

        CALL small_circle_features(v0,v1,em,z_crcl,xpole,z_scal)
        sph_coarse%ts(j)%sub_triangle3%edge2%cos_rel_co_lat  = z_crcl
        sph_coarse%ts(j)%sub_triangle3%edge2%sc_pole         = xpole
        sph_coarse%ts(j)%sub_triangle3%edge2%edge_primal_arc = z_scal

        ! second edge
        !------------

        v0 = sph_coarse%ts(j)%edge2%edge_center
        v1 = sph_coarse%ts(j)%vertex0%vertex
        c0 = sph_coarse%ts(j)%sub_triangle0%triangle_center
        c1 = sph_coarse%ts(j)%sub_triangle0%neighbor2%triangle_center

        xpole%x(:)=0.0_wp
        z_crcl=0.0_wp
        z_scal=0.0_wp ! just for definition at the subroutine call
        em = edge_midpoint(c0,c1,z_crcl,xpole)
        sph_coarse%ts(j)%sub_triangle0%edge2%edge_center = em

        CALL small_circle_features(v0,v1,em,z_crcl,xpole,z_scal)
        sph_coarse%ts(j)%sub_triangle0%edge2%cos_rel_co_lat  = z_crcl
        sph_coarse%ts(j)%sub_triangle0%edge2%sc_pole         = xpole
        sph_coarse%ts(j)%sub_triangle0%edge2%edge_primal_arc = z_scal

      ENDIF

    ENDDO

  END SUBROUTINE scc_outer
  !
  !-------------------------------------------------------------------------
  !-------------------------------------------------------------------------
  !
  !
  !>
  !! This subroutine computes the cos_rel_co_lat and pole of a small circle.
  !!
  !! This subroutine computes the cos_rel_co_lat and pole of a small circle
  !! given by the vertices and a midpoint.
  !! The procedure to obtain those values uses the point @f$\vec{O'}@f$, which is
  !! defined as the center of the small circle (for reference see figure 6 of
  !! Kimerling et al., 2002) and lies somewhere on the line between the sphere
  !! midpoint and the pole @f$\vec{O'}=\vec{P}\cos\alpha@f$. Here, @f$\alpha@f$ is the
  !! co-latitude of the small circle with pole @f$\vec{P}@f$. If we know the edge
  !! midpoint @f$\vec{E}@f$ and the vertices @f$\vec{V_1}, \vec{V_2}@f$ of the small
  !! circle arc segment, the point @f$\vec{O'}@f$ can be obtained by simple geometrical
  !! thoughts from
  !! @f{equation}{
  !! \vec{O'}= \frac{1}{2}\left((\vec{M}+\vec{E})+(\vec{M}-\vec{E})
  !! \frac{|\vec{D}|^2}{|\vec{M}-\vec{E}|^2}\right)
  !! @f}
  !! with @f$\vec{D}=(\vec{V_1}-\vec{V_2})/2@f$ and @f$\vec{M}=(\vec{V_1}+\vec{V_2})/2@f$.
  !! Because we know that @f$|\vec{P}|=1@f$, we obtain
  !! @f{equation}{
  !! \cos\alpha=|\vec{O'}|.
  !! @f}
  !! The length of the small circle arc is
  !! @f{equation}{
  !! l_{sc}=r\varphi = |\vec{E}-\vec{O'}|\varphi=|\vec{E}-\vec{O'}|
  !! \frac{(\vec{V_1}-\vec{O'})\cdot (\vec{V_2}-\vec{O'})}
  !!  {|\vec{V_1}-\vec{O'}||\vec{V_2}-\vec{O'}|}.
  !! @f}
  !!
  !! @par Revision History
  !! Initial release by Almut Gassmann (2007-03)
  !! Modification by Almut Gassmann (2007-04)
  !! - finish if to the arc length is more than a semi circle
  !!
  SUBROUTINE small_circle_features(v0,v1,em,crcl,p,scal)
    !
    TYPE (t_cartesian_coordinates), INTENT(in) :: v0, v1, em
    !
    TYPE (t_cartesian_coordinates), INTENT(out) :: p
    REAL (wp), INTENT(out) :: crcl,scal
    !
    TYPE (t_cartesian_coordinates) :: z_mean, z_dist, z_diff, z_add, z_o, x0, x1
    REAL (wp) :: z_x2, z_z2
    !
    !-------------------------------------------------------------------------
    !
    z_mean%x = 0.5_wp * (v0%x + v1%x)
    z_dist%x = 0.5_wp * (v0%x - v1%x)
    z_diff%x = z_mean%x - em%x
    z_add%x  = z_mean%x + em%x
    z_x2     = DOT_PRODUCT(z_dist%x,z_dist%x)
    z_z2     = DOT_PRODUCT(z_diff%x,z_diff%x)

    z_o%x    = 0.5_wp * ( z_add%x + z_diff%x * z_x2/z_z2)

    crcl     = SQRT(DOT_PRODUCT(z_o%x,z_o%x))

    x0%x = v0%x-z_o%x
    x1%x = v1%x-z_o%x

    ! Be sure that no arcs greater than pi occur
    IF(DOT_PRODUCT(x1%x,x1%x)<dot_product(z_diff%x,z_diff%x))THEN
      CALL finish ('small_circle_features',&
        & 'Error and abort: small circle arc would be greater than pi')
    ENDIF

    IF ( crcl < 1.0e-8_wp ) THEN
      crcl=0.0_wp
      p%x =0.0_wp
      scal = arc_length (x0, x1)
    ELSE
      p%x = z_o%x/crcl
      scal = arc_length (x0, x1) * SQRT(DOT_PRODUCT(em%x-z_o%x,em%x-z_o%x))
    ENDIF

  END SUBROUTINE small_circle_features

  !-------------------------------------------------------------------------
  !-------------------------------------------------------------------------

END MODULE mo_small_circle_c

