!>
!!  Contains routines used by the grid optimization procedure.
!!
!!  More or less faithful translation of what was done by Thuburn.
!!
!! @par Revision History
!! Initial version  by Ross Heikes, Colorado State University, 1995
!! Present implementation originally written in f77 by John Thuburn,
!! Reading University (approx 1997).
!! @par
!! Recoded in f90 and adapted to the ICON data structure by
!!  Luca Bonaventura,  MPI-M, Hamburg, August 2005
!!  Modified by Thomas Heinze,  DWD, Offenbach, 2005-09
!!  - cleaning up code according to latest programming guide
!!  - debugging for grid optimization
!!  Modifications by Th.Heinze, DWD (2006-10-25):
!!  - changed index to idx in TYPE declarations of edge_info, vertex_info,
!!    triangle_info, grid_cells, grid_edges and grid_vertices
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
MODULE mo_optimize
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
  USE mo_math_constants
  USE mo_base_datatypes, ONLY: t_triangle, vertex, edge, t_spheres
  USE mo_exception,      ONLY: message_text, message, finish
  USE mo_math_types,     ONLY: t_cartesian_coordinates, t_geographical_coordinates
  USE mo_math_utilities, ONLY: cc2gc, gc2cc, circum_center
  !thh
  USE mo_topology,       ONLY: spheres_on_levels
  !thh

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: optimize_heikes

CONTAINS

  !--------------------------------------------------------------------
  !-------------------------------------------------------------------------
  !
  !
  !>
  !!  Executes the grid optimization on.
  !!
  !!  Executes the grid optimization on
  !!  the grid stored in the ptr_sph data structure.
  !!
  !! @par Revision History
  !! Initial version  by Ross Heikes, Colorado State University, 1995
  !! Present implementation originally written in f77 by John Thuburn,
  !! Reading University (approx 1997).
  !! Recoded in f90 and adapted to the ICON data structure by
  !! Luca Bonaventura,  MPI-M, Hamburg, (2005-08)
  !! Modified by Thomas Heinze,  DWD, Offenbach, 2005-09
  !! - cleaning up code according to latest programming guide
  !! - debugging for grid optimization
  !!
  SUBROUTINE optimize_heikes(ptr_sph)
    !
    REAL(wp), PARAMETER:: ftol=1.e-5_wp  ! tolerance for iterative minimization
    ! process
    INTEGER, PARAMETER :: maxiter =10    ! maximum number of iterations:
    ! usually more than 10 is just a waste
    ! !INPUT/OUTPUT  PARAMETERS:
    TYPE(t_spheres)               :: ptr_sph

    TYPE(t_spheres), POINTER :: ptr_parentsph            ! thh
    TYPE(t_triangle), POINTER :: ptr_tt, ptr_t1, ptr_t2
    TYPE(vertex), POINTER :: ptr_vv
    TYPE(t_cartesian_coordinates) :: ve1, ve2, ce1, ce2

    REAL(wp) :: z_cost, z_totcst, z_totcst0, z_ototcst

    INTEGER :: vertices(3)
    INTEGER :: i_lev
    INTEGER :: j_ie1, j_iter

    !--------------------------------------------------------------------

    i_lev=ptr_sph%level

    !thh
    ! the optimization operates on interior subtriangle
    !
    ptr_parentsph=>spheres_on_levels(i_lev-1)
    !
    !thh

    !  compute initial value of cost (penalty) function
    !  by summing contributions for each edge
    !
    z_totcst0 = 0.0_wp
    z_cost = 0.0_wp

    DO j_ie1 = 0, ptr_sph%no_edges-1

      ve1=ptr_sph%es(j_ie1)%vertex0%vertex
      ve2=ptr_sph%es(j_ie1)%vertex1%vertex

      !thh
      !
      ! set triangle centers
      ! first triangle
      !
      ptr_tt=>ptr_sph%es(j_ie1)%triangle0
      ptr_tt%triangle_center = circum_center( ptr_tt%vertex0%vertex,            &
        & ptr_tt%vertex1%vertex,            &
        & ptr_tt%vertex2%vertex )
      ce1=ptr_tt%triangle_center
      !
      ! set triangle centers
      ! second triangle
      !
      ptr_tt=>ptr_sph%es(j_ie1)%triangle1
      ptr_tt%triangle_center = circum_center( ptr_tt%vertex0%vertex,            &
        & ptr_tt%vertex1%vertex,            &
        & ptr_tt%vertex2%vertex )
      ce2=ptr_tt%triangle_center
      !
      ! triangles_center have not been set before!
      !       ce1=ptr_sph%es(j_ie1)%triangle0%triangle_center
      !       ce2=ptr_sph%es(j_ie1)%triangle1%triangle_center
      !thh

      z_cost=penalty_edge(ve1,ve2,ce1,ce2)

      z_totcst0 = z_totcst0 + z_cost

    ENDDO

    WRITE (message_text,'(a,e23.17)')  &
      & 'Total cost function before minimization: J = ', z_totcst0
    CALL message ('', TRIM(message_text))

    z_ototcst=0._wp
    z_totcst=z_totcst0

    main_opt_loop:  DO j_iter = 1, maxiter

      z_ototcst=z_totcst
      !
      ! Loop over all edges of parent sphere:
      ! The inner subtriangle of both adjacent triangles have one common vertex. This
      ! vertex is the vertex on the edge in mind.
      !
      DO j_ie1= 0, ptr_parentsph%no_edges-1
        !
        ! point onto inner subtriangle of both adjacent triangles
        !
        ptr_t1=>ptr_parentsph%es(j_ie1)%triangle0%sub_triangle2
        ptr_t2=>ptr_parentsph%es(j_ie1)%triangle1%sub_triangle2
        !
        ! store vertices of triangle1
        !
        vertices(1)=ptr_t2%vertex0%info%idx
        vertices(2)=ptr_t2%vertex1%info%idx
        vertices(3)=ptr_t2%vertex2%info%idx
        !
        ! compare them with vertices of triangle0
        !
        IF (ANY(vertices == ptr_t1%vertex0%info%idx)) THEN
          ptr_vv=>ptr_t1%vertex0
        ELSE IF (ANY(vertices == ptr_t1%vertex1%info%idx)) THEN
          ptr_vv=>ptr_t1%vertex1
        ELSE
          ptr_vv=>ptr_t1%vertex2
        ENDIF
        !
        ! do optimization
        !
        CALL powell(ptr_vv, ftol, ptr_sph)

      ENDDO
      !
      ! recompute value of cost function
      !
      z_totcst = 0.0_wp

      DO j_ie1 = 0, ptr_sph%no_edges-1

        ve1=ptr_sph%es(j_ie1)%vertex0%vertex
        ve2=ptr_sph%es(j_ie1)%vertex1%vertex
        ce1=ptr_sph%es(j_ie1)%triangle0%triangle_center
        ce2=ptr_sph%es(j_ie1)%triangle1%triangle_center

        z_cost=penalty_edge(ve1,ve2,ce1,ce2)

        z_totcst = z_totcst + z_cost

      ENDDO

      WRITE (message_text,'(a,i2,a,e23.17)')                                    &
        & 'Total cost function at iteration ', j_iter,                     &
        & ' of minimization: J = ', z_totcst
      CALL message ('', TRIM(message_text))
      !
      !   check with tolerance
      !
      IF(ABS(z_totcst-z_ototcst)<=ftol*z_totcst0) EXIT

    ENDDO main_opt_loop

    WRITE (message_text,'(a,e23.17)') 'Final cost function J = ', z_totcst
    CALL message ('', TRIM(message_text))

    WRITE (message_text,'(a,i2)') 'Done grid level ', i_lev
    CALL message ('', TRIM(message_text))

  END SUBROUTINE optimize_heikes
  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Find the contribution to the Heikes+Randall penalty function.
  !!
  !! Find the contribution to the Heikes+Randall penalty function
  !! from a given edge. It was a subroutine in Thuburn's
  !! implementation, it has been turned into a function
  !! that returns the value <i>p_cost</i>.<br>
  !! x1,x2 are vertices of the edge on the dual grid
  !! (i.e., triangle centers), x3,x4 are vertices
  !! of the edge on the primal grid
  !!
  !! @par Revision History
  !!  Originally developed  by R.Heikes (1995).
  !!  Present version obtained from J.Thuburn (2002).
  !!  Adapted to fortran90  and recoded  by Luca Bonaventura(2002-5).
  !!  Modified by Thomas Heinze (2005-09):
  !!  - cleaning up code according to latest programming guide
  !!  - exchanged points of primal and dual grid due to description
  !!
  !FUNCTION penalty_edge(x1,x2,x3,x4) RESULT(p_cost)
  FUNCTION penalty_edge(x3,x4,x1,x2) result(p_cost)
    !
    TYPE(t_cartesian_coordinates),INTENT(in) :: x1,x2,x3,x4

    REAL (wp) :: p_cost

    TYPE(t_cartesian_coordinates) :: x, xm, dx
    REAL(wp)                    :: l2, r2, mag

    !-----------------------------------------------------------------------
    !
    ! Calculate length and midpoint of dual edge
    ! whose vertices are x1,x2 (centers of adjacent primal cells)
    !
    dx%x(:) = x1%x(:) - x2%x(:)
    xm%x(:) = 0.5_wp*(x1%x(:)+x2%x(:))
    mag = 1.0_wp/SQRT(DOT_PRODUCT(xm%x(:),xm%x(:)))
    xm%x(:) = xm%x(:)*mag
    !
    ! Calculate normalization l2
    !
    l2 = DOT_PRODUCT(dx%x(:),dx%x(:))
    !
    ! Find midpoint of primal edge whose vertices are x3,x4
    !
    x%x(:) = 0.5_wp*(x3%x(:)+x4%x(:))
    mag = 1.0_wp/SQRT(DOT_PRODUCT(x%x(:),x%x(:)))
    x%x(:) = x%x(:)*mag
    !
    ! Contribution to penalty function
    !
    dx%x(:) = x%x(:) - xm%x(:)
    r2 = DOT_PRODUCT(dx%x(:),dx%x(:))/l2
    p_cost = r2*r2

  END FUNCTION penalty_edge
  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Minimizes the (grid dependent) function <i>penalty_vertex</i>.
  !!
  !! Minimizes the (grid dependent) function <i>penalty_vertex</i>
  !! using Powell's method, following its  'Numerical Recipes'
  !! implementation. The basic minimization call is that
  !! to the subroutine <i>brent</i>, which calls in turn
  !! <i>penalty_vertex</i>.
  !!
  !! @par Revision History
  !! Initial version  by Ross Heikes, Colorado State University, 1995
  !! Present implementation originally written in f77 by John Thuburn,
  !! Reading University (approx 1997).
  !! Recoded in f90 and adapted to the ICON data structure by
  !! Luca Bonaventura,  MPI-M, Hamburg, August 2005
  !! Modified by Thomas Heinze, DWD, Offenbach, according to latest programming
  !! guide (2005-09)
  !!
  SUBROUTINE powell(ptr_vv,ftol,ptr_sph)  !(if1,igrid,ftol,ngchck)
    !
!!$    REAL(wp), PARAMETER :: cgold = 0.3819660_wp   ! golden ratio
    REAL(wp), PARAMETER :: tol = 1.0e-5_wp        ! tolerance

    INTEGER, PARAMETER :: n = 2, itmax = 100

    TYPE(vertex), POINTER :: ptr_vv                 ! vertex

    REAL (wp)                      :: ftol                   ! given tolerance

    TYPE(t_spheres)                  :: ptr_sph                ! current sphere

    TYPE(t_triangle), POINTER :: ptr_tt                 ! triangle
    TYPE(t_geographical_coordinates) :: position

    REAL(wp) :: p(2), xi(2,2), pt(2), ptt(2), xit(2)
    REAL(wp) :: fp, fptt, ax, bx, cx,  del, fret, t, fac

    INTEGER :: j, j_iter, ibig, j_i, i_lev

    !-----------------------------------------------------------------------
    !
    i_lev=ptr_sph%level                              ! current level
    !
    ! Brackets for 1D search
    !
    fac = 2.0_wp**(-(i_lev+2))  ! igrid
    del = 1.18_wp*fac
    ax = -del
    bx = 0.0_wp
    !
    ! original :   BX=(2.0*CGOLD-1.0)*DEL
    !
    cx = del

    position=cc2gc(ptr_vv%vertex)

    p(1) = position%lon
    p(2) = position%lat

    fret = penalty_vertex(ptr_vv,ptr_sph)
    !
    ! Initialize search directions
    !
    DO j_i = 1, n
      DO j = 1, n
        xi(j,j_i) = 0.0_wp
      ENDDO
    ENDDO

    xi(1,1) = 1.0_wp/COS(p(2))
    xi(2,2) = 1.0_wp
    !
    ! Save initial point
    !
    DO j = 1, n
      pt(j) = p(j)
    ENDDO
    !
    ! MAIN minimization LOOP
    !
    j_iter = 0

    main_loop: DO
      j_iter = j_iter + 1
      fp = fret
      ibig = 0
      del = 0.0_wp
      !
      ! Loop over all directions
      !
      DO j_i = 1, n
        !
        ! Copy the direction
        !
        DO j = 1, n
          xit(j) = xi(j,j_i)
        ENDDO
        fptt = fret
        !
        ! Line minimization along direction XIT from P
        !
        CALL brent(ax,bx,cx,tol,n,p,xit,fret,ptr_vv,ptr_sph)

        IF (ABS(fptt-fret) > del) THEN
          del = ABS(fptt-fret)
          ibig = j_i
        ENDIF

      ENDDO
      !
      ! Are we finished?
      !
      IF (2.0_wp*abs(fp-fret) <= ftol*(ABS(fp)+ABS(fret))) EXIT
      !
      ! Construct extrapolated point and average direction moved
      !
      DO j = 1, n
        ptt(j) = 2.0_wp*p(j) - pt(j)
        xit(j) = p(j) - pt(j)
        pt(j) = p(j)
      ENDDO
      !
      ! Function value at extrapolated point
      !
      position%lon = ptt(1)
      position%lat = ptt(2)

      ptr_vv%vertex=gc2cc(position)

      fptt = penalty_vertex(ptr_vv,ptr_sph)

      IF (fptt < fp) THEN

        t = 2.0_wp*(fp-2.0_wp*fret+fptt)*sqrt(MAX(fp-fret-del,0.0_wp)) &
          & - del*sqrt(fp-fptt)

        IF (t < 0.0_wp) THEN

          CALL brent(ax,bx,cx,tol,n,p,xit,fret,ptr_vv,ptr_sph)
          DO j = 1, n
            xi(j,ibig) = xi(j,n)
            xi(j,n) = xit(j)
          ENDDO

        ENDIF

      ENDIF

      IF (j_iter < itmax) CYCLE

      CALL message ('powell minimization', 'Too many iterations in POWELL')
      CALL finish  ('powell minimization', 'Too many iterations in POWELL')

    ENDDO main_loop
    !
    ! Reset long and lat to best value
    !
    position%lon = p(1)
    position%lat = p(2)

    ptr_vv%vertex=gc2cc(position)
    !
    ! Reset positions of affected vertices
    !
    ptr_tt=>ptr_vv%triangle0

    ptr_tt%triangle_center = circum_center(   ptr_tt%vertex0%vertex, &
      & ptr_tt%vertex1%vertex, &
      & ptr_tt%vertex2%vertex )

    ptr_tt=>ptr_vv%triangle1

    ptr_tt%triangle_center = circum_center(   ptr_tt%vertex0%vertex, &
      & ptr_tt%vertex1%vertex, &
      & ptr_tt%vertex2%vertex )

    ptr_tt=>ptr_vv%triangle2

    ptr_tt%triangle_center = circum_center(   ptr_tt%vertex0%vertex, &
      & ptr_tt%vertex1%vertex, &
      & ptr_tt%vertex2%vertex )
    ptr_tt=>ptr_vv%triangle3

    ptr_tt%triangle_center = circum_center(   ptr_tt%vertex0%vertex, &
      & ptr_tt%vertex1%vertex, &
      & ptr_tt%vertex2%vertex )

    ptr_tt=>ptr_vv%triangle4

    ptr_tt%triangle_center = circum_center(   ptr_tt%vertex0%vertex, &
      & ptr_tt%vertex1%vertex, &
      & ptr_tt%vertex2%vertex )

    ptr_tt=>ptr_vv%triangle5

    ptr_tt%triangle_center = circum_center(   ptr_tt%vertex0%vertex, &
      & ptr_tt%vertex1%vertex, &
      & ptr_tt%vertex2%vertex )

  END SUBROUTINE powell


  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Find the contribution to the Heikes+Randall penalty function.
  !!
  !! Find the contribution to the Heikes+Randall penalty function
  !! from a given dual cell face, identified by the corresponding
  !! vertex of a primal triangle.
  !!
  !! @par Revision History
  !!  Originally developed  by R.Heikes (1995).
  !!  Present version obtained from J.Thuburn (2002).
  !!  Adapted to fortran90  and recoded by Luca Bonaventura(2002-5).
  !!  Modified by Thomas Heinze, DWD, Offenbach, (2005-09)
  !!  - according to latest programming guide
  !!  - vertex checking
  !!
  FUNCTION penalty_vertex(v1,ptr_sph) result (totcst)
    !

    TYPE(vertex), INTENT(in)     :: v1

    TYPE(t_spheres), INTENT(inout) :: ptr_sph

    REAL (wp)                    :: totcst

    TYPE(t_triangle), POINTER :: ptr_tt
    TYPE(edge), POINTER :: ptr_ee
    TYPE(t_cartesian_coordinates)  :: cc(12,2)

    REAL(wp)                     :: z_cost

    INTEGER :: ielist(12)
    INTEGER :: iv1, ie1, je1

    !-----------------------------------------------------------------------
    !
    ! Find triangles  that are affected  by changes
    ! in vertex v1 and update their coordinates
    ! (vertex v1  will never be the center of a pentagon,
    ! since only the grid levels higher than
    ! the root level are optimized,  so only consider hexagons).
    !
    ! List edges that are affected by vertex v1 in the ielist(:)
    ! array; there are a total of 12 affected edges
    !
    je1 = 0
    !
    ! first triangle
    !
    ptr_tt=>v1%triangle0
    iv1=ptr_tt%info%idx

    ptr_ee=>ptr_tt%edge0
    ie1 =ptr_ee%info%idx
    IF (.not. (ANY(ielist(1:je1) == ie1))) THEN
      je1 = je1+1
      ielist(je1) = ie1
    ENDIF

    ptr_ee=>ptr_tt%edge1
    ie1 =ptr_ee%info%idx
    IF (.not. (ANY(ielist(1:je1) == ie1))) THEN
      je1 = je1+1
      ielist(je1) = ie1
    ENDIF

    ptr_ee=>ptr_tt%edge2
    ie1 =ptr_ee%info%idx
    IF (.not. (ANY(ielist(1:je1) == ie1))) THEN
      je1 = je1+1
      ielist(je1) = ie1
    ENDIF
    !
    ! second triangle
    !
    ptr_tt=>v1%triangle1
    iv1=ptr_tt%info%idx

    ptr_ee=>ptr_tt%edge0
    ie1 =ptr_ee%info%idx
    IF (.not. (ANY(ielist(1:je1) == ie1))) THEN
      je1 = je1+1
      ielist(je1) = ie1
    ENDIF

    ptr_ee=>ptr_tt%edge1
    ie1 =ptr_ee%info%idx
    IF (.not. (ANY(ielist(1:je1) == ie1))) THEN
      je1 = je1+1
      ielist(je1) = ie1
    ENDIF

    ptr_ee=>ptr_tt%edge2
    ie1 =ptr_ee%info%idx
    IF (.not. (ANY(ielist(1:je1) == ie1))) THEN
      je1 = je1+1
      ielist(je1) = ie1
    ENDIF
    !
    ! third triangle
    !
    ptr_tt=>v1%triangle2
    iv1=ptr_tt%info%idx

    ptr_ee=>ptr_tt%edge0
    ie1 =ptr_ee%info%idx
    IF (.not. (ANY(ielist(1:je1) == ie1))) THEN
      je1 = je1+1
      ielist(je1) = ie1
    ENDIF

    ptr_ee=>ptr_tt%edge1
    ie1 =ptr_ee%info%idx
    IF (.not. (ANY(ielist(1:je1) == ie1))) THEN
      je1 = je1+1
      ielist(je1) = ie1
    ENDIF

    ptr_ee=>ptr_tt%edge2
    ie1 =ptr_ee%info%idx
    IF (.not. (ANY(ielist(1:je1) == ie1))) THEN
      je1 = je1+1
      ielist(je1) = ie1
    ENDIF
    !
    !  fourth triangle
    !
    ptr_tt=>v1%triangle3
    iv1=ptr_tt%info%idx

    ptr_ee=>ptr_tt%edge0
    ie1 =ptr_ee%info%idx
    IF (.not. (ANY(ielist(1:je1) == ie1))) THEN
      je1 = je1+1
      ielist(je1) = ie1
    ENDIF

    ptr_ee=>ptr_tt%edge1
    ie1 =ptr_ee%info%idx
    IF (.not. (ANY(ielist(1:je1) == ie1))) THEN
      je1 = je1+1
      ielist(je1) = ie1
    ENDIF

    ptr_ee=>ptr_tt%edge2
    ie1 =ptr_ee%info%idx
    IF (.not. (ANY(ielist(1:je1) == ie1))) THEN
      je1 = je1+1
      ielist(je1) = ie1
    ENDIF
    !
    !  fifth triangle
    !
    ptr_tt=>v1%triangle4
    iv1=ptr_tt%info%idx

    ptr_ee=>ptr_tt%edge0
    ie1 =ptr_ee%info%idx
    IF (.not. (ANY(ielist(1:je1) == ie1))) THEN
      je1 = je1+1
      ielist(je1) = ie1
    ENDIF

    ptr_ee=>ptr_tt%edge1
    ie1 =ptr_ee%info%idx
    IF (.not. (ANY(ielist(1:je1) == ie1))) THEN
      je1 = je1+1
      ielist(je1) = ie1
    ENDIF

    ptr_ee=>ptr_tt%edge2
    ie1 =ptr_ee%info%idx
    IF (.not. (ANY(ielist(1:je1) == ie1))) THEN
      je1 = je1+1
      ielist(je1) = ie1
    ENDIF
    !
    ! sixth triangle
    !
    ptr_tt=>v1%triangle5
    iv1=ptr_tt%info%idx

    ptr_ee=>ptr_tt%edge0
    ie1 =ptr_ee%info%idx
    IF (.not. (ANY(ielist(1:je1) == ie1))) THEN
      je1 = je1+1
      ielist(je1) = ie1
    ENDIF

    ptr_ee=>ptr_tt%edge1
    ie1 =ptr_ee%info%idx
    IF (.not. (ANY(ielist(1:je1) == ie1))) THEN
      je1 = je1+1
      ielist(je1) = ie1
    ENDIF

    ptr_ee=>ptr_tt%edge2
    ie1 =ptr_ee%info%idx

    IF (.not. (ANY(ielist(1:je1) == ie1))) THEN
      je1 = je1+1
      ielist(je1) = ie1
    ENDIF
    !
    ! Calculate new triangle centers of all effected triangles
    !
    DO je1=1,12

      ie1 = ielist(je1)

      ptr_tt=>ptr_sph%es(ie1)%triangle0
      !
      !thh
      ! Check if triangle out of bounce
      !
      IF (ptr_tt%info%idx>(ptr_sph%no_triangles -1 )) THEN

        WRITE (message_text,'(a,i4,a,i2,a,i4,a)')                            &
          & 'Index problem: Triangle no.',ptr_tt%info%idx,                   &
          & ' does not belong to level',ptr_sph%level,' (max index:',          &
          & ptr_sph%no_triangles-1,')'
        CALL message ('',TRIM(message_text))
        CALL finish ('mo_optimize', message_text)
      ENDIF
      !thh

      cc(je1,1)= circum_center(      ptr_tt%vertex0%vertex,                     &
        & ptr_tt%vertex1%vertex,                     &
        & ptr_tt%vertex2%vertex)

      ptr_tt=>ptr_sph%es(ie1)%triangle1

      IF (ptr_tt%info%idx>(ptr_sph%no_triangles -1 )) THEN

        WRITE (message_text,'(a,i4,a,i2,a,i4,a)')                            &
          & 'Index problem: Triangle no.',ptr_tt%info%idx,                   &
          & ' does not belong to level',ptr_sph%level,' (max index:',          &
          & ptr_sph%no_triangles-1,')'
        CALL message ('',TRIM(message_text))
        CALL finish ('mo_optimize', message_text)
      ENDIF

      cc(je1,2)= circum_center(      ptr_tt%vertex0%vertex,                     &
        & ptr_tt%vertex1%vertex,                     &
        & ptr_tt%vertex2%vertex)

    ENDDO
    !
    ! Now add up contributions to penalty function
    !
    totcst = 0.0_wp

    DO je1 = 1, 12

      ie1 = ielist(je1)

      z_cost=penalty_edge(ptr_sph%es(ie1)%vertex0%vertex,                        &
        & ptr_sph%es(ie1)%vertex1%vertex, cc(je1,1), cc(je1,2))

      totcst = totcst + z_cost

    ENDDO


  END FUNCTION penalty_vertex

  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Perform line minimization of the function <i>penalty_vertex</i>.
  !!
  !! Perform line minimization of the function <i>penalty_vertex</i>
  !! using Brent's method.
  !! Based on Numerical Recipes
  !! AX,BX,CX is assumed to bracket the minimum, i.e. F(BX)<F(AX),F(CX)
  !!
  !! @par Revision History
  !!  Originally developed  by R.Heikes (1995).
  !!  Present version obtained from J.Thuburn (2002)
  !!  Adapted to fortran90 by Luca Bonaventura  (2002).
  !!  Modified by Thomas Heinze, DWD, Offenbach, (2005-09)
  !!  - slight changes in prolog according to programming guide
  !!
  SUBROUTINE brent(ax, bx, cx, tol, n, p, xit, fmin, ptr_vv,ptr_sph)
    !
    INTEGER, PARAMETER :: itmax = 100             ! Maximum no. of iterations
    REAL (wp), PARAMETER :: cgold = 0.3819660_wp    ! Golden ratio
    REAL (wp), PARAMETER :: zeps  = EPSILON(0.0_wp) ! Divide by zero protection

    INTEGER, INTENT(in)            :: n
    REAL(wp), INTENT(in)           :: ax, bx, cx, tol

    TYPE(t_spheres)                  :: ptr_sph

    REAL(wp), INTENT(inout)        :: xit(n), p(n)

    REAL(wp), INTENT(out)          :: fmin


    TYPE(vertex)                   :: ptr_vv
    TYPE(t_geographical_coordinates) :: position

    REAL (wp) :: xmin, a, b, x, u, v, w, fx, fu, fv, fw, tol1, tol2
    REAL (wp) :: e, d, pp, q, r, etemp, xm

    INTEGER :: j_iter

    !-----------------------------------------------------------------------
    !
    ! Sort A, B, into ascending order
    !
    IF (ax < bx) THEN
      a = ax
      b = cx
    ELSE
      a = cx
      b = ax
    ENDIF
    !
    ! Initialize search points and function values
    !
    x = bx
    w = bx
    v = bx

    position%lon = p(1) + x*xit(1)
    position%lat = p(2) + x*xit(2)

    ptr_vv%vertex=gc2cc(position)

    fx = penalty_vertex(ptr_vv,ptr_sph)

    fw = fx
    fv = fx

    e = 0.0_wp
    d = 0.0_wp
    !
    ! main iteration loop
    !
    j_iter = 0
    main_loop: DO

      j_iter = j_iter + 1
      xm = 0.5_wp*(a+b)
      !
      ! Check for convergence:   TOL1=TOL*ABS(X)+ZEPS
      ! In the case where the min happens to fall at x=0 but P is
      ! non-zero, better to put a typical P value (1.0) instead of x
      !
      tol1 = tol + zeps
      tol2 = 2.0_wp*tol1

      IF (ABS(x-xm) <= tol2-0.5_wp*(b-a)) THEN
        xmin = x
        fmin = fx
        EXIT
      ENDIF
      !
      ! Construct a trial parabolic fit
      !
      IF (ABS(e) > tol1) THEN
        r = (x-w)*(fx-fv)
        q = (x-v)*(fx-fw)
        pp = (x-v)*q - (x-w)*r
        q = 2.0_wp*(q-r)

        IF (q > 0.0_wp) pp = -pp

        q = ABS(q)
        etemp = e
        e = d

        IF (ABS(pp) >= ABS(0.5_wp*q*etemp)  &
          & .or. pp <= q*(a-x)         &
          & .or. pp >= q*(b-x)) THEN
          IF (x >= xm) THEN
            e = a - x
          ELSE
            e = b - x
          ENDIF
          d = cgold*e
        ELSE
          d = pp/q
          u = x + d
          IF (u-a < tol2 .or. b-u < tol2) d = SIGN(tol1,xm-x)
        ENDIF
      ELSE
        IF (x >= xm) THEN
          e = a - x
        ELSE
          e = b - x
        ENDIF
        d = cgold*e
      ENDIF

      IF (ABS(d) >= tol1) THEN
        u = x + d
      ELSE
        u = x + SIGN(tol1,d)
      ENDIF

      position%lon = p(1) + u*xit(1)
      position%lat = p(2) + u*xit(2)

      ptr_vv%vertex=gc2cc(position)

      fu = penalty_vertex(ptr_vv,ptr_sph)

      IF (fu <= fx) THEN

        IF (u >= x) THEN
          a = x
        ELSE
          b = x
        ENDIF

        v = w
        w = x
        x = u
        fv = fw
        fw = fx
        fx = fu

      ELSE

        IF (u < x) THEN
          a = u
        ELSE
          b = u
        ENDIF

        IF (fu <= fw .or. w == x) THEN
          v = w
          w = u
          fv = fw
          fw = fu
        ELSEIF (fu <= fv .or. v == x .or. v == w) THEN
          v = u
          fv = fu
        ENDIF

      ENDIF

      IF (j_iter < itmax) CYCLE

      CALL message ('subroutine BRENT', 'Maximum iterations exceeded in BRENT')
      CALL finish ('subroutine BRENT', 'Maximum iterations exceeded in BRENT')

    ENDDO main_loop
    !
    ! Save vector displacement and new position of min
    !
    xit(:) = xmin*xit(:)
    p(:) = p(:)+xit(:)

  END SUBROUTINE brent

END MODULE mo_optimize

