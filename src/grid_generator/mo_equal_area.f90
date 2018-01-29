!>
!! Equal area subdivision aimes at constructing successively refined triangles.
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
!! Equal area subdivision aimes at constructing successively refined triangles
!! of equal size per level by using small circles as triangle edges.
!! The idea is presented in Song et al.&nbsp;(2002). We adapt it here for the grid
!! generator. We shortly repeat the steps outlined in detail in the mentioned
!! article.
!!
!! <ol>
!! <li>  Knowing the following quantities of a parent triangle
!! <ul>
!! <li>  vertices,
!! <li>  edge centers,
!! <li>  poles of the small circle forming the edge
!! <li>  cosines of the relative co-latitude of the small circle
!! <li>  triangle area
!! </ul>
!! we are able to set already some values to the unknown children of the parent.
!! We will enumerate the children with numbers (0,1,2,3) such that the child
!! number 2 is surrounded by the child triangles 0 (top), 1 (left) and 3 (right).
!! As in the grid generator, the edges and vertices are counted counterclockwise
!! starting at the top of each triangle. Now we find for instance the edge 0
!! between the vertices 0 and 1.
!!
!! Considering our generated child-triangles we know already the vertices and the
!! edges of the small triangles that coincinde with edges of the parent triangle.
!! The area of the small triangles is known as one fourth of the parent triangle
!! area.
!!
!! <li>  Now we consider each of the vertex triangles (0,1,3) indepentently.
!! The small circle adjacent to the inner triangle 2 is unknown and is searched
!! to be defined such that the prescribed triangle area is achieved.
!! We can compute triangle areas for triangles formed by great circles by
!! @f{equation}{
!! A_{gc}=\angle A +\angle B +\angle C - \pi
!! @f}
!! and the wedge areas of polar triangles build from a small circle and two
!! equal great circles
!! @f{equation}{
!! W_{sc}=\angle P ( 1-\cos\alpha )=:\angle P ( 1-q )
!! @f}
!! where @f$\alpha@f$ is the co-latitude of the small circle and @f$P@f$ the pole of the
!! small circle.
!!
!! An arbitrary triangle area with small circle edges  is  now given by a great
!! circle area and the sum of the semi-lunes formed by the great circle edge
!! and the small circle edge between two vertices
!! @f{equation}{
!! A_{sc} = A_{gc}+ S_0+S_1+S_2.
!! \label {areasc}
!! @f}
!! The semi-lune is computed from the difference between the polar wedges formed
!! by a small circle or a great circle, respectively
!! @f{equation}{
!! |S| = W_{sc}-W_{gc}
!! @f}
!! the great circle wedge ist given by
!! @f{equation}{
!! W_{gc}=A_{gc}^p=\angle P + 2 \angle D -\pi.
!! @f}
!! The semi-lune in that definition is always positive, since the great circle
!! always bends nearer to the pole. The the sign of the semi-lunes in
!! (\\ref{areasc}) has to be carefully defined: if the small circle edge is
!! outside the great circle edge viewed from the center of the triangle, the
!! semi-lune is given positve sign, and vice versa.
!!
!! With the given fomulas, the semi-lunes of the known edges of the vertex
!! triangles  may be computed. With the aid of the target area @f$A_t@f$ (one
!! fourth of the parent triangle area) we find for the unknown semi-lune adjacent
!! to triangle 2 (here for the vertex triangle 1)
!! @f{equation}{
!! S_2=A_t-A_{gc}-S_0-S_1
!! @f}
!!
!! <li>  Now the small circle properties of the small circle edge adjacent to the
!! triangle 2 left to be determined. We have several unknowns: The cosine of the
!! relative co-latitude, the polar angle and the angles forming the polar great
!! circle triangle. Regarding the way a great circle polar triangle area is
!! computed, we find
!! @f{equation}{
!! A_{gc}^p=\arccos(z_0) + \arccos (z_1) + \arccos(z_1) -\pi
!! @f}
!! with
!! @f{eqnarray*}{
!! z_0&=&\frac{-1}{|P\times V_1||V_2\times P|}\left((P\cdot V_2)(V_1\cdot P)-
!!                                                  (V_1\cdot V_2)\right)\\\
!! z_1&=&\frac{-1}{|V_1\times V_2||P\times V_1|}\left((V_1\cdot P)(V_2\cdot V_1)-
!!                                                  (P\cdot V_2)\right)\\\
!! z_2&=&\frac{-1}{|V_2\times P||V_1\times V_2|}\left((V_2\cdot V_1)(P\cdot V_2)-
!!                                                  (P\cdot V_1)\right)
!! @f}
!! The vector and scalar products of the vertices @f$V_1@f$ and @f$V_2@f$ are known to
!! @f{equation}{
!! V_1\cdot V_2 = \cos(l_{gc})=:r,\quad\quad\quad |V_1\times V_2| =
!! \sin(l_{gc})=\sqrt{1-\cos^2(l_{gc})}=\sqrt{1-r^2}
!! @f}
!! with @f$l_{gc}@f$ the length of the great circle arc. Furthermore we know that
!! @f{eqnarray*}{
!! P\cdot V_1=P\cdot V_2 & =  &\cos\alpha=:q\\\
!! |P\times V_1|=|P\times V_2| &=&\sin\alpha= \sqrt{1-\cos^2\alpha}=:\sqrt{1-q^2}
!! @f}
!!
!! Inserting all the found relationsships into the equation for the area of the
!! semi lune we find
!! @f{eqnarray*}{
!! S_2&=&\angle P(1-q)-(\angle P +2\angle D -\pi)\\\
!! S_2&=&\pi-q\angle P-2\angle D\\\
!! S_2&=&\pi-q\arccos(z_0)-2\arccos(z_1)\\\
!! S_2&=&\pi-q\arccos(\frac{r-q^2}{1-q^2})-2\arccos(\frac{q u}{\sqrt{1-q^2}})
!! @f}
!! with
!! @f{equation}{
!! u=(1-r)/\sqrt{1-r^2}=\sqrt{\frac{1-r}{1+r}}
!! @f}
!!
!! This transcendent equation may now be solved for @f$q@f$ by the secant or Newton
!! method in iteration steps.
!!
!! <li>  We now have found @f$q@f$ and for completeness we need the pole of the
!! small circle edge. We have to solve the system
!! @f{equation}{
!! P^2=1, \quad\quad\quad V_1\cdot P=q, \quad\quad\quad V_2\cdot P=q
!! @f}
!! This gives a quadratic equation for the pole coordinates. Ensuring that no
!! division by a small number occurs, we chose the equation to solve according
!! to the largest amount of the determinant @f$d_i@f$, given by
!! @f$d=V_1\times V_2@f$. The quadratic equation to solve is then
!! @f{equation}{
!! P_i^2+P_i2q\frac{f_{i+1}d_{i-1}-f_{i-1}d_{i+1}}{\sum_id_i^2}+
!! \frac{q^2(f_{i+1}^2+f_{i-1}^2)-d_i^2}{\sum_id_i^2}=0
!! @f}
!! where @f$f=V_1-V_2@f$.
!! The other components of the pole recover from
!! @f{eqnarray*}{
!! P_{i+1}&=&(-qf_{i-1}+P_i d_{i+1})/d_i\\\
!! P_{i-1}&=&(qf_{i+1}+P_i d_{i-1})/d_i
!! @f}
!! That pole is chosen as the significant one which is nearer to the triangle
!! center if the semi lune is positive, or which is further to the triangle
!! center if the semi lune is negative.
!!
!! <li>  Last, we need to find the edge midpoint on the small scale edges of all
!! child triangles so that the vertices may be recovered in the next level. This
!! is done exactly like described in Song et al. (2002) using their equations
!! (2)-(6).
!!
!! <b>Reference</b><br>
!! Song, L., Kimerling, A.&nbsp;J., and Sahr,K., 2002: Developing an equal area global
!! grid by small circle subdivision. In M.&nbsp;Goodchild ans A.&nbsp;J.&nbsp;Kimerling (Eds.),
!! Discrete Global Grids. Santa Barbara, CA, USA: National Center for Geographic
!! Information and Analysis.
!!
!! </ol>
!!
!! @par Revision History
!! Initial release by Almut Gassmann (2007-03-22)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_equal_area
  !-------------------------------------------------------------------------
  !
  !    ProTeX FORTRAN source: Style 2
  !    modified for ICON project, DWD/MPI-M 2006
  !
  !-------------------------------------------------------------------------
  !
  !

  USE mo_kind,          ONLY: wp
  USE mo_math_types,    ONLY: t_cartesian_coordinates
  USE mo_math_utilities, ONLY: &
    & triangle_area,         &
    & circum_center,         &
    & arc_length,            &
    & cos_arc_length,        &
    & vector_product
  USE mo_math_constants,ONLY: pi, eps
  USE mo_base_datatypes,ONLY: t_triangle, edge
  USE mo_topology,      ONLY: root, icosahedron
  USE mo_exception,     ONLY: finish, message_text
  USE mo_icosahedron_geometry, ONLY: get_icosahedron_vertex

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: equal_area,       &
    & equal_area_nine,  &
    & equal_area_twen5, &
    & edge_midpoint,    &
    & sc_triangle_area

CONTAINS

  !-------------------------------------------------------------------------
  !-------------------------------------------------------------------------
  !
  !

  !>
  !! This subroutine is thought for four fold sudivision.
  !!
  !!
  !! @par Revision History
  !! Initial release by Almut Gassmann (2007-03-22)
  !! Treatment of icosahedron level by Almut Gassmann (2008-09-02)
  !!
  SUBROUTINE equal_area(cn)



    TYPE(t_triangle), POINTER :: cn

    TYPE(t_triangle), POINTER :: s0, s1, s2, s3
    REAL (wp) :: z_sum_sl, z_sl, z_cgal
    INTEGER :: idx
    !-------------------------------------------------------------------------
    !-----------------------
    ! Section 0: preparation
    !-----------------------

    IF (cn%info%level == -1) THEN
      idx = cn%info%triangle_number
      s0 => root(0,idx)%t
      s1 => root(1,idx)%t
      s2 => root(2,idx)%t
      s3 => root(3,idx)%t
    ELSE
      s0 => cn%sub_triangle0
      s1 => cn%sub_triangle1
      s2 => cn%sub_triangle2
      s3 => cn%sub_triangle3
    ENDIF
    !---------------------------------------------
    ! Section 1: compute new vertices on the edges
    !---------------------------------------------
    IF (cn%info%level == -1) THEN
      ! Icosahedron
      s2%vertex0%vertex = edge_midpoint( &
        & get_icosahedron_vertex(1,idx),&
        & get_icosahedron_vertex(2,idx),&
        & 0.0_wp,cn%edge1%sc_pole)
      s2%vertex1%vertex = edge_midpoint( &
        & get_icosahedron_vertex(2,idx),&
        & get_icosahedron_vertex(0,idx),&
        & 0.0_wp,cn%edge2%sc_pole)
      s2%vertex2%vertex = edge_midpoint( &
        & get_icosahedron_vertex(0,idx),&
        & get_icosahedron_vertex(1,idx),&
        & 0.0_wp,cn%edge0%sc_pole)
      s0%vertex0%vertex = get_icosahedron_vertex(0,idx)
      s1%vertex1%vertex = get_icosahedron_vertex(1,idx)
      s3%vertex2%vertex = get_icosahedron_vertex(2,idx)
    ELSE
      s2%vertex0%vertex%x(:) = cn%edge1%edge_center%x(:)
      s2%vertex1%vertex%x(:) = cn%edge2%edge_center%x(:)
      s2%vertex2%vertex%x(:) = cn%edge0%edge_center%x(:)
      s0%vertex0%vertex%x(:) = cn%vertex0%vertex%x(:)
      s1%vertex1%vertex%x(:) = cn%vertex1%vertex%x(:)
      s3%vertex2%vertex%x(:) = cn%vertex2%vertex%x(:)
    ENDIF

    !--------------------------------------------------
    ! Section 2: set already known small scale features
    !--------------------------------------------------


    s0%edge0%sc_pole%x(:)    = cn%edge0%sc_pole%x(:)
    s0%edge0%cos_rel_co_lat  = cn%edge0%cos_rel_co_lat
    s0%edge0%edge_primal_arc = cn%edge0%edge_primal_arc*0.5_wp
    s0%edge0%edge_center     = edge_midpoint(s0%vertex0%vertex,&
      & s0%vertex1%vertex,&
      & s0%edge0%cos_rel_co_lat,&
      & s0%edge0%sc_pole)
    s0%edge2%sc_pole%x(:)    = cn%edge2%sc_pole%x(:)
    s0%edge2%cos_rel_co_lat  = cn%edge2%cos_rel_co_lat
    s0%edge2%edge_primal_arc = cn%edge2%edge_primal_arc*0.5_wp
    s0%edge2%edge_center     = edge_midpoint(s0%vertex0%vertex,&
      & s0%vertex2%vertex,&
      & s0%edge2%cos_rel_co_lat,&
      & s0%edge2%sc_pole)


    s1%edge0%sc_pole%x(:)    = cn%edge0%sc_pole%x(:)
    s1%edge0%cos_rel_co_lat  = cn%edge0%cos_rel_co_lat
    s1%edge0%edge_primal_arc = cn%edge0%edge_primal_arc*0.5_wp
    s1%edge0%edge_center     = edge_midpoint(s1%vertex0%vertex,&
      & s1%vertex1%vertex,&
      & s1%edge0%cos_rel_co_lat,&
      & s1%edge0%sc_pole)
    s1%edge1%sc_pole%x(:)    = cn%edge1%sc_pole%x(:)
    s1%edge1%cos_rel_co_lat  = cn%edge1%cos_rel_co_lat
    s1%edge1%edge_primal_arc = cn%edge1%edge_primal_arc*0.5_wp
    s1%edge1%edge_center     = edge_midpoint(s1%vertex1%vertex,&
      & s1%vertex2%vertex,&
      & s1%edge1%cos_rel_co_lat,&
      & s1%edge1%sc_pole)


    s3%edge2%sc_pole%x(:)    = cn%edge2%sc_pole%x(:)
    s3%edge2%cos_rel_co_lat  = cn%edge2%cos_rel_co_lat
    s3%edge2%edge_primal_arc = cn%edge2%edge_primal_arc*0.5_wp
    s3%edge2%edge_center     = edge_midpoint(s3%vertex0%vertex,&
      & s3%vertex2%vertex,&
      & s3%edge2%cos_rel_co_lat,&
      & s3%edge2%sc_pole)
    s3%edge1%sc_pole%x(:)    = cn%edge1%sc_pole%x(:)
    s3%edge1%cos_rel_co_lat  = cn%edge1%cos_rel_co_lat
    s3%edge1%edge_primal_arc = cn%edge1%edge_primal_arc*0.5_wp
    s3%edge1%edge_center     = edge_midpoint(s3%vertex1%vertex,&
      & s3%vertex2%vertex,&
      & s3%edge1%cos_rel_co_lat,&
      & s3%edge1%sc_pole)

    !-------------------------------------------
    ! Section 2: compute centers of subtriangles
    !-------------------------------------------

    s0%triangle_center = circum_center(s0%vertex0%vertex,&
      & s0%vertex1%vertex,&
      & s0%vertex2%vertex)
    s1%triangle_center = circum_center(s1%vertex0%vertex,&
      & s1%vertex1%vertex,&
      & s1%vertex2%vertex)
    s2%triangle_center = circum_center(s2%vertex0%vertex,&
      & s2%vertex1%vertex,&
      & s2%vertex2%vertex)
    s3%triangle_center = circum_center(s3%vertex0%vertex,&
      & s3%vertex1%vertex,&
      & s3%vertex2%vertex)

    s0%triangle_area = 0.25_wp*cn%triangle_area
    s1%triangle_area = 0.25_wp*cn%triangle_area
    s2%triangle_area = 0.25_wp*cn%triangle_area
    s3%triangle_area = 0.25_wp*cn%triangle_area

    !---------------------------------------------------------------------
    ! Section 3: solve for the small circles
    !            a) compute the sum of all semi_lunes for a triangle
    !            b) compute the areas of the unknown semi_lunes
    !            c) we need the cosine of the great arc length, too
    !            d) we need the solver range (note: qmin is always 0)
    !            e) call solver for the cosine of the relative co-latitude
    !            f) recover the poles, small arc lengths, and edge centers
    !---------------------------------------------------------------------

    ! Triangle s0
    !------------
    z_cgal   = cos_arc_length(s0%vertex1%vertex,s0%vertex2%vertex)
    z_sum_sl = triangle_area(s0%vertex0%vertex,s0%vertex1%vertex,s0%vertex2%vertex)
    z_sum_sl = s0%triangle_area - z_sum_sl
    z_sl     = semi_lune(z_sum_sl,s0%edge0,s0%edge2,s0%triangle_center)
    CALL solve_small_circle(z_sl,z_cgal,s0%edge1%cos_rel_co_lat)
    CALL compute_edge(z_sl,z_cgal,s0%edge1,s0%triangle_center)

    ! Triangle s1
    !------------
    z_cgal   = cos_arc_length(s1%vertex0%vertex,s1%vertex2%vertex)
    z_sum_sl = triangle_area(s1%vertex0%vertex,s1%vertex1%vertex,s1%vertex2%vertex)
    z_sum_sl = s1%triangle_area - z_sum_sl
    z_sl     = semi_lune(z_sum_sl,s1%edge0,s1%edge1,s1%triangle_center)
    CALL solve_small_circle(z_sl,z_cgal,s1%edge2%cos_rel_co_lat)
    CALL compute_edge(z_sl,z_cgal,s1%edge2,s1%triangle_center)

    ! Triangle s3
    !------------
    z_cgal   = cos_arc_length(s3%vertex0%vertex,s3%vertex1%vertex)
    z_sum_sl = triangle_area(s3%vertex0%vertex,s3%vertex1%vertex,s3%vertex2%vertex)
    z_sum_sl = s3%triangle_area - z_sum_sl
    z_sl     = semi_lune(z_sum_sl,s3%edge1,s3%edge2,s3%triangle_center)
    CALL solve_small_circle(z_sl,z_cgal,s3%edge0%cos_rel_co_lat)
    CALL compute_edge(z_sl,z_cgal,s3%edge0,s3%triangle_center)


  END SUBROUTINE equal_area

  !-------------------------------------------------------------------------
  !-------------------------------------------------------------------------
  !
  !
  !>
  !! This function computes the semi lune of an unknown triangle edge if the.
  !!
  !! This function computes the semi lune of an unknown triangle edge if the
  !! sum of sum of all semi lunes of the triangle is given. Then, the semi lunes
  !! of the known edges are computed and from the residual the unknown semi lune
  !! is obtained.
  !!
  FUNCTION semi_lune(sum_sl, e0, e1, tc) result(sl)
    !
    REAL (wp), INTENT(in) :: sum_sl ! area sum of the semi lunes
    TYPE (edge), INTENT(in) :: e0, e1  ! edges with with known features
    TYPE (t_cartesian_coordinates), INTENT(in) :: tc ! triangle center
    !
    REAL (wp) :: sl ! area of the semi lune with unknown features
    !
    REAL (wp) :: z_p_sc_area0, z_p_sc_area1, z_p_gc_area0, z_p_gc_area1, &
      & z_dists, z_distg, z_sign0, z_sign1
    TYPE (t_cartesian_coordinates) :: z_dummy, z_help
    !
    !-------------------------------------------------------------------------

    z_dummy%x = 0.0_wp

    ! areas of polar wedges with small circle arcs
    z_p_sc_area0 = sc_triangle_area( e0 )
    z_p_sc_area1 = sc_triangle_area( e1 )

    ! areas of polar wedges with great circle arcs
    IF ( e0%cos_rel_co_lat  /= 0.0_wp) THEN
      z_p_gc_area0 = triangle_area(e0%sc_pole,e0%vertex0%vertex,e0%vertex1%vertex)

      z_help = edge_midpoint(e0%vertex0%vertex,e0%vertex1%vertex,0.0_wp,z_dummy)
      z_distg = arc_length(z_help,tc)
      z_dists = arc_length(e0%edge_center,tc)

      z_sign0 = SIGN(1.0_wp,z_dists-z_distg)
    ELSE
      z_p_gc_area0 = z_p_sc_area0
      z_sign0 = 0.0_wp
    ENDIF
    IF ( e1%cos_rel_co_lat  /= 0.0_wp) THEN
      z_p_gc_area1 = triangle_area(e1%sc_pole,e1%vertex0%vertex,e1%vertex1%vertex)

      z_help = edge_midpoint(e1%vertex0%vertex,e1%vertex1%vertex,0.0_wp,z_dummy)
      z_distg = arc_length(z_help,tc)
      z_dists = arc_length(e1%edge_center,tc)

      z_sign1 = SIGN(1.0_wp,z_dists-z_distg)
    ELSE
      z_p_gc_area1 = z_p_sc_area1
      z_sign1 = 0.0_wp
    ENDIF

    sl = sum_sl - z_sign0*(z_p_sc_area0-z_p_gc_area0)&
      & - z_sign1*(z_p_sc_area1-z_p_gc_area1)

  END FUNCTION semi_lune
  !
  !-------------------------------------------------------------------------
  !-------------------------------------------------------------------------
  !
  !
  !>
  !! This function computes the edge midpoint if the adjacent vertices and small.
  !!
  !! This function computes the edge midpoint if the adjacent vertices and small
  !! circle properties are given.
  !!
  !! @par Revision History
  !! Initial release by Almut Gassmann (2007-03-22)
  !!
  FUNCTION edge_midpoint(v0,v1,crcl,p) result (em)
    !
    TYPE (t_cartesian_coordinates), INTENT(in) :: v0, v1, p
    REAL (wp), INTENT(in) :: crcl
    !
    TYPE (t_cartesian_coordinates) :: em
    !
    TYPE (t_cartesian_coordinates) :: z_o
    REAL (wp) :: z_help
    !
    !--------------------------------------------------------------------------

    em%x(:)  = 0.5_wp*(v0%x+v1%x)
    z_o%x(:) = p%x(:) * crcl
    z_help   = SQRT((1.0_wp-crcl*crcl)/&
      & DOT_PRODUCT(em%x(:)-z_o%x(:),em%x(:)-z_o%x(:)))
    em%x(:)  = z_o%x(:)+(em%x(:)-z_o%x(:))*z_help
    ! Ensure modulus 1 for midpoints
    em%x(:)  = em%x(:)/SQRT(DOT_PRODUCT(em%x(:),em%x(:)))

  END FUNCTION edge_midpoint
  !-------------------------------------------------------------------------
  !-------------------------------------------------------------------------
  !
  !
  !>
  !! This function computes the polar wedge area with a small circle as edge.
  !!
  !!
  !! @par Revision History
  !! Initial release by Almut Gassmann (2007-03-22)
  !!
  FUNCTION sc_triangle_area( ed ) result (scta)
    !
    TYPE(edge), INTENT(in) :: ed
    !
    REAL (wp) :: scta
    !
    REAL (wp) :: z_r, z_q
    !
    !-------------------------------------------------------------------------
    z_r  = cos_arc_length(ed%vertex0%vertex, ed%vertex1%vertex)
    z_q  = ed%cos_rel_co_lat
    scta = (1.0_wp-z_q) * ACOS( (z_r-z_q*z_q)/(1.0_wp-z_q*z_q) )

  END FUNCTION sc_triangle_area
  !-------------------------------------------------------------------------
  !-------------------------------------------------------------------------
  !
  !
  !>
  !! This subroutine computes the cosine of the co-latitude of an edge given.
  !!
  !! This subroutine computes the cosine of the co-latitude of an edge given
  !! its semi-lune and the great circle distance to the adjacent vertices.
  !! This is achieved by computing the root of the transcendent equation
  !! arising from system (22) of Song et al. (2002). A combination of bisection
  !! and the Newton-Raphson method is applied for that purpose (see Press et al.)
  !!
  !! @par Revision History
  !! Initial release by Almut Gassmann (2007-03-22)
  !!
  SUBROUTINE solve_small_circle(sl,r,q)
    !
    REAL (wp), INTENT(in) :: sl, &! Area of semi lune
      & r    ! cosine of great arc length
    !
    REAL (wp), INTENT(out):: q    ! cosine of relative co latitude
    !
    REAL (wp) :: z_qmax
    !
    !-------------------------------------------------------------------------
    !
    ! The right bound of the function has an infinite gradient ->
    ! avoid division by zero by setting the qmax-value a bit to
    ! the left. Note: chosing zero as left value gives always a
    ! positive function value because f=(1-cos(sl)) is always
    ! positive. The right function value will always be negative.
    z_qmax  = SQRT(0.5_wp*(r+1.0_wp))
    z_qmax  = z_qmax-eps

    q = rtsafe(0.0_wp,z_qmax,eps)

  CONTAINS

    !---------------------------------------------------------------
    ! adapted from "rtsafe" of "Numerical Recipes" (Press et al.)

    FUNCTION rtsafe(x1,x2,xacc) result (rtsafe_val)

      IMPLICIT NONE

      REAL(wp), INTENT(in) :: x1, x2, xacc
      REAL(wp) :: rtsafe_val

      INTEGER, PARAMETER :: maxit = 100
      INTEGER :: j
      REAL(wp) :: df, dx, dxold, f, fh, fl, temp, xh, xl

      CALL funcd(x1, fl, df)
      CALL funcd(x2, fh, df)

      IF ((fl > 0.0_wp .and. fh > 0.0_wp) .or. (fl < 0.0_wp .and. fh < 0.0_wp)) THEN
        WRITE(message_text,'(a,6e12.6)') &
          & 'Root must be bracketed ', x1, fl, x2, fh, r, sl
        CALL finish('rtsafe', message_text)
      ENDIF
      IF(fl == 0.0_wp) THEN
        rtsafe_val = x1
        RETURN
      ELSE IF(fh == 0.0_wp) THEN
        rtsafe_val = x2
        RETURN
      ELSE IF(fl < 0.0_wp) THEN
        xl = x1
        xh = x2
      ELSE
        xh = x1
        xl = x2
      ENDIF

      rtsafe_val = 0.5_wp*(x1+x2)
      dxold = ABS(x2-x1)
      dx = dxold

      CALL funcd(rtsafe_val, f, df)

      DO j = 1, maxit
        IF(((rtsafe_val-xh)*df-f)*((rtsafe_val-xl)*df-f) >= 0.0_wp .or. &
          & ABS(2.0_wp*f) > ABS(dxold*df) ) THEN
          dxold = dx
          dx = 0.5_wp*(xh-xl)
          rtsafe_val = xl+dx
          IF(xl == rtsafe_val) RETURN
        ELSE
          dxold = dx
          dx = f/df
          temp = rtsafe_val
          rtsafe_val = rtsafe_val-dx
          IF(temp == rtsafe_val) RETURN
        ENDIF
        IF(ABS(dx) < xacc) RETURN
        CALL funcd(rtsafe_val, f, df)
        IF(f < 0.0_wp) THEN
          xl = rtsafe_val
        ELSE
          xh = rtsafe_val
        ENDIF
      END DO

      CALL finish('rtsafe','exceeding maximum iterations')

    END FUNCTION rtsafe

    !-------------------------------------------

    SUBROUTINE funcd(x,f,fd)
      REAL (wp), INTENT(in) :: x
      REAL (wp), INTENT(out) :: f, fd  ! function and first derivative
      REAL (wp) :: z_h1, z_h2, z_h3, z_h4

      z_h1 = pi-ABS(sl)
      z_h2 = (1.0_wp-r)/(1.0_wp+r)
      z_h3 = 1.0_wp-x*x
      z_h4 = (r-x*x)/z_h3

      f  = z_h1-x*acos(z_h4)-2.0_wp*acos(x*sqrt(z_h2/z_h3))
      fd = -ACOS(z_h4)+2.0_wp/z_h3*(SQRT(z_h2/(1.0_wp-x*x*(1.0_wp+z_h2)))-&
        & x*x*sqrt((1.0_wp-r)/(1.0_wp+r-2.0_wp*x*x)))

    END SUBROUTINE funcd

    !--------------------------------------------

  END SUBROUTINE solve_small_circle
  !-------------------------------------------------------------------------
  !-------------------------------------------------------------------------
  !
  !
  !>
  !! This subroutine computes all missing edge properties if the cosine of the.
  !!
  !! This subroutine computes all missing edge properties if the cosine of the
  !! relative co-latitude is found. The pole of the small circle arc is computed
  !! by a quadratic equation (see scientific documentation). Afterwards the
  !! arc length and the edge center are obtained.
  !!
  !! @par Revision History
  !! Initial release by Almut Gassmann (2007-03-22)
  !! Modification by Almut Gassmann (2007-04-17)
  !! - finish if the small circle arc bends to much
  !!
  SUBROUTINE compute_edge(sl,cgal,ed,tc)
    !
    REAL(wp), INTENT(in)  :: sl,     &! area of the semi lune
      & cgal     ! cosine great arc length
    TYPE(t_cartesian_coordinates), INTENT(in) :: tc  ! triangle center
    !
    TYPE(edge), INTENT(inout) :: ed   ! edge
    !
    TYPE(t_cartesian_coordinates) :: x_pole(2), x_f, x_d, x_a, x_b
    REAL(wp) :: z_det, z_helpa, z_helpb, z_dist1, z_dist2, z_crcl
    INTEGER :: i, im, ip
    !
    !-------------------------------------------------------------------------

    z_crcl = ed%cos_rel_co_lat
    x_a    = ed%vertex0%vertex
    x_b    = ed%vertex1%vertex

    IF (z_crcl<eps) THEN
      ed%sc_pole%x(:)    = 0.0_wp
      ed%cos_rel_co_lat  = 0.0_wp
      ed%edge_primal_arc = arc_length(x_a,x_b)
      ed%edge_center     = edge_midpoint(x_a,x_b,ed%cos_rel_co_lat,ed%sc_pole)
      RETURN
    ENDIF

    ! solve for pole coords (eqns:23,24)
    x_f%x(:) = (x_a%x(:)-x_b%x(:))*z_crcl
    x_d      = vector_product (x_a,x_b)
    z_det    = DOT_PRODUCT(x_d%x,x_d%x)

    i  = MAXLOC(ABS(x_d%x),1)
    im = i-1
    ip = i+1
    IF (im == 0 ) im = 3
    IF (ip == 4 ) ip = 1

    z_helpa         = ( x_f%x(ip)*x_d%x(im)-x_f%x(im)*x_d%x(ip) ) / z_det
    z_helpb         = ( x_f%x(ip)**2 + x_f%x(im)**2-x_d%x(i)**2 ) / z_det
    x_pole(1)%x(i)  = - z_helpa + SQRT(z_helpa**2 - z_helpb)
    x_pole(2)%x(i)  = - z_helpa - SQRT(z_helpa**2 - z_helpb)
    x_pole(1)%x(ip) = (-x_f%x(im)+x_pole(1)%x(i)*x_d%x(ip))/x_d%x(i)
    x_pole(2)%x(ip) = (-x_f%x(im)+x_pole(2)%x(i)*x_d%x(ip))/x_d%x(i)
    x_pole(1)%x(im) = ( x_f%x(ip)+x_pole(1)%x(i)*x_d%x(im))/x_d%x(i)
    x_pole(2)%x(im) = ( x_f%x(ip)+x_pole(2)%x(i)*x_d%x(im))/x_d%x(i)

    x_pole(1)%x(:)=x_pole(1)%x(:)/SQRT(DOT_PRODUCT(x_pole(1)%x(:),x_pole(1)%x(:)))
    x_pole(2)%x(:)=x_pole(2)%x(:)/SQRT(DOT_PRODUCT(x_pole(2)%x(:),x_pole(2)%x(:)))

    z_dist1 = arc_length(x_pole(1),tc)
    z_dist2 = arc_length(x_pole(2),tc)

    ! cartes. coords of pole (eqn:25)
    IF ((sl > 0.0_wp .and. z_dist2 > z_dist1) .or. &
      & (sl < 0.0_wp .and. z_dist2 < z_dist1)) THEN
      ed%sc_pole = x_pole(1)
    ELSE
      ed%sc_pole = x_pole(2)
    ENDIF

    ! Test if the small circle is not larger than a semi circle
    IF( triangle_area(ed%sc_pole,x_a,x_b) < ABS(sl)) THEN
      CALL finish( 'compute_edge',&
        & 'Error and abort: small circle arc would be larger than pi')
    ENDIF

    ! arc length of small circle
    ed%edge_primal_arc = SQRT(1.0_wp-z_crcl*z_crcl)*&
      & ACOS( (cgal-z_crcl*z_crcl)/(1.0_wp-z_crcl*z_crcl))

    ! edge center
    ed%edge_center = edge_midpoint(x_a,x_b,z_crcl,ed%sc_pole)

  END SUBROUTINE compute_edge
  !
  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  !
  !
  !>
  !! Nine fold subdivision on root level.
  !!
  !! The symmetrical icosahedron face allows us to solve only for 2 kinds of
  !! semi lunes. One semi lune is in the corner of the icosahedron and 2 equal
  !! semi lunes connect the midpoint with the vertices of the first semi lune.
  !! All other features may then be reconstructed easily.
  !!
  !! @par Revision History
  !! Initial release by Almut Gassmann (2007-03-22)
  !!
  SUBROUTINE equal_area_nine(idx)
    !
    REAL (wp), PARAMETER :: z_ratio = 1.0_wp/3.0_wp
    !
    INTEGER, INTENT(in) :: idx  ! icosahedron index
    !
    TYPE (t_triangle), POINTER :: s0, s1, ico
    REAL (wp) :: z_sl, z_sl0, z_sl1, z_cgal, z_cgal0, z_cgal1
    INTEGER :: j
    !
    !------------------------------------------------------------------------------
    !------------------------
    ! Section 0: Preparations
    !------------------------
    s0    => root(0,idx)%t
    s1    => root(2,idx)%t
    ico   => icosahedron(idx)%t

    !------------------------------
    ! Section 1: Determine vertices
    !------------------------------

    root(0,idx)%t%vertex1%vertex = divide_icos(ico%vertex0%vertex,&
      & ico%vertex1%vertex,z_ratio)
    root(0,idx)%t%vertex2%vertex = divide_icos(ico%vertex0%vertex,&
      & ico%vertex2%vertex,z_ratio)
    root(0,idx)%t%vertex0%vertex = ico%vertex0%vertex

    root(4,idx)%t%vertex0%vertex = divide_icos(ico%vertex1%vertex,&
      & ico%vertex0%vertex,z_ratio)
    root(4,idx)%t%vertex2%vertex = divide_icos(ico%vertex1%vertex,&
      & ico%vertex2%vertex,z_ratio)
    root(4,idx)%t%vertex1%vertex = ico%vertex1%vertex

    root(8,idx)%t%vertex0%vertex = divide_icos(ico%vertex2%vertex,&
      & ico%vertex0%vertex,z_ratio)
    root(8,idx)%t%vertex1%vertex = divide_icos(ico%vertex2%vertex,&
      & ico%vertex1%vertex,z_ratio)
    root(8,idx)%t%vertex2%vertex = ico%vertex2%vertex

    root(2,idx)%t%vertex0%vertex = circum_center(ico%vertex0%vertex,&
      & ico%vertex1%vertex,&
      & ico%vertex2%vertex)

    !--------------------------------------
    ! Section 2: Triangle_centers and areas
    !--------------------------------------
    DO j = 0, 8
      root(j,idx)%t%triangle_area = ico%triangle_area/9.0_wp
      root(j,idx)%t%triangle_center = circum_center(root(j,idx)%t%vertex0%vertex,&
        & root(j,idx)%t%vertex1%vertex,&
        & root(j,idx)%t%vertex2%vertex)
    ENDDO

    !------------------------------------------------------------
    ! Section 3: For simplicity give all edges default properties
    !------------------------------------------------------------
    z_sl=0.0_wp
    z_cgal=0.0_wp
    DO j = 0,8
      root(j,idx)%t%edge0%cos_rel_co_lat = 0.0_wp
      root(j,idx)%t%edge1%cos_rel_co_lat = 0.0_wp
      root(j,idx)%t%edge2%cos_rel_co_lat = 0.0_wp
      CALL compute_edge(z_sl,z_cgal,root(j,idx)%t%edge0,&
        & root(j,idx)%t%triangle_center)
      CALL compute_edge(z_sl,z_cgal,root(j,idx)%t%edge1,&
        & root(j,idx)%t%triangle_center)
      CALL compute_edge(z_sl,z_cgal,root(j,idx)%t%edge2,&
        & root(j,idx)%t%triangle_center)
    ENDDO

    !-------------------------------------------------------------------------------
    ! Section 4: Prototype triangles s0 and s1:
    !            solve for cosines of the relative co-latitude
    !-------------------------------------------------------------------------------
    z_sl0   = s0%triangle_area-triangle_area(s0%vertex0%vertex,&
      & s0%vertex1%vertex,&
      & s0%vertex2%vertex)
    z_cgal0 = cos_arc_length(s0%vertex1%vertex,s0%vertex2%vertex)
    CALL solve_small_circle(z_sl0,z_cgal0,s0%edge1%cos_rel_co_lat)

    z_sl1   = -z_sl0 ! the semi-lune has now opposite sign
    z_sl1   = 0.5_wp*(s1%triangle_area-z_sl1-&
      & triangle_area(s1%vertex0%vertex,&
      & s1%vertex1%vertex,&
      & s1%vertex2%vertex))
    z_cgal1 = cos_arc_length(s1%vertex0%vertex,s1%vertex1%vertex)
    CALL solve_small_circle(z_sl1,z_cgal1,s1%edge0%cos_rel_co_lat)


    !----------------------------------------------------------------------
    ! Section 5: All relevant values are found, compute all missing values.
    !----------------------------------------------------------------------

    root(4,idx)%t%edge2%cos_rel_co_lat = s0%edge1%cos_rel_co_lat
    root(8,idx)%t%edge0%cos_rel_co_lat = s0%edge1%cos_rel_co_lat
    CALL compute_edge(z_sl0,z_cgal0,root(0,idx)%t%edge1,&
      & root(0,idx)%t%triangle_center)
    CALL compute_edge(z_sl0,z_cgal0,root(4,idx)%t%edge2,&
      & root(4,idx)%t%triangle_center)
    CALL compute_edge(z_sl0,z_cgal0,root(8,idx)%t%edge0,&
      & root(8,idx)%t%triangle_center)

    root(2,idx)%t%edge2%cos_rel_co_lat = s1%edge0%cos_rel_co_lat
    root(5,idx)%t%edge0%cos_rel_co_lat = s1%edge0%cos_rel_co_lat
    root(5,idx)%t%edge1%cos_rel_co_lat = s1%edge0%cos_rel_co_lat
    root(7,idx)%t%edge1%cos_rel_co_lat = s1%edge0%cos_rel_co_lat
    root(7,idx)%t%edge2%cos_rel_co_lat = s1%edge0%cos_rel_co_lat
    CALL compute_edge(z_sl1,z_cgal1,root(2,idx)%t%edge0,&
      & root(2,idx)%t%triangle_center)
    CALL compute_edge(z_sl1,z_cgal1,root(2,idx)%t%edge2,&
      & root(2,idx)%t%triangle_center)
    CALL compute_edge(z_sl1,z_cgal1,root(5,idx)%t%edge0,&
      & root(5,idx)%t%triangle_center)
    CALL compute_edge(z_sl1,z_cgal1,root(5,idx)%t%edge1,&
      & root(5,idx)%t%triangle_center)
    CALL compute_edge(z_sl1,z_cgal1,root(7,idx)%t%edge1,&
      & root(7,idx)%t%triangle_center)
    CALL compute_edge(z_sl1,z_cgal1,root(7,idx)%t%edge2,&
      & root(7,idx)%t%triangle_center)

  END SUBROUTINE equal_area_nine
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  !
  !
  !>
  !! 25-fold subdivision on root level.
  !!
  !! To achieve symmetry in the subdivison, the equal area subdivsion is done
  !! in two steps. In the figure besides, double edged fourfold triangles are
  !! subdivided first. Afterwards, the small triangle edges are treated.
  !!
  !! @verbatim
  !!                  //\\
  !!                 //0_\\
  !!                //\ 2/\\
  !!               //1_\/_3\\
  !!              /==========\
  !!             / \ 5//\\ 7/ \
  !!            /_4_\//_6\\/_8_\
  !!          //\\10//\12/\\14//\\
  !!         //9 \\//11\/13\\//15\\
  !!        //----\/========\/----\\
  !!       //\ 17 /\\19/\21//\ 23 /\\
  !!      //16\  /18\\/20\//22\  /24\\
  !!     //==========\\--//===========
  !! @endverbatim
  !!
  !! @par Revision History
  !! Initial release by Almut Gassmann (2007-03-22)
  !!
  SUBROUTINE equal_area_twen5(idx)
    !
    INTEGER, INTENT(in) :: idx  ! icosahedron index
    !
    TYPE (t_triangle), POINTER :: ico
    INTEGER :: j, ij, istat, iastat
    REAL (wp) :: z_ratio, z_triangle_area, &
      & z_sl0, z_sl1, z_sl2, z_sl3, &
      & z_cgal0, z_cgal1, z_cgal2, z_cgal3, &
      & z_crcl0, z_crcl1, z_crcl2, z_crcl3
    TYPE (t_cartesian_coordinates) :: z_v0, z_v1, z_v2, z_center
    TYPE (edge) :: z_ed0, z_ed1, z_ed2, z_ed3, z_ed4, z_ed5
    !
    !------------------------------------------------------------------------------
    !-----------------------------------------------------------------------
    ! Section 0: Allocate dummy edge components and abbreviate icosahedron%t
    !-----------------------------------------------------------------------

    ico => icosahedron(idx)%t

    iastat=0
    ALLOCATE(z_ed0%vertex0,stat=istat)
    iastat=ABS(istat)+iastat
    ALLOCATE(z_ed0%vertex1,stat=istat)
    iastat=ABS(istat)+iastat
    ALLOCATE(z_ed1%vertex0,stat=istat)
    iastat=ABS(istat)+iastat
    ALLOCATE(z_ed1%vertex1,stat=istat)
    iastat=ABS(istat)+iastat
    ALLOCATE(z_ed2%vertex0,stat=istat)
    iastat=ABS(istat)+iastat
    ALLOCATE(z_ed2%vertex1,stat=istat)
    iastat=ABS(istat)+iastat
    ALLOCATE(z_ed3%vertex0,stat=istat)
    iastat=ABS(istat)+iastat
    ALLOCATE(z_ed3%vertex1,stat=istat)
    iastat=ABS(istat)+iastat
    ALLOCATE(z_ed4%vertex0,stat=istat)
    iastat=ABS(istat)+iastat
    ALLOCATE(z_ed4%vertex1,stat=istat)
    iastat=ABS(istat)+iastat
    ALLOCATE(z_ed5%vertex0,stat=istat)
    iastat=ABS(istat)+iastat
    ALLOCATE(z_ed5%vertex1,stat=istat)
    iastat=ABS(istat)+iastat
    IF( iastat /= 0 ) THEN
      WRITE(message_text,'(a)') 'Problem in allocating dummy vertices'
      CALL finish('equal_area_twen5',message_text)
    ENDIF


    !--------------------------------------------------------------
    ! Section 1: Solve for each 4 triangles in the 3 vertex regions
    !--------------------------------------------------------------

    z_triangle_area = ico%triangle_area*4.0_wp/25.0_wp
    z_ratio = 0.4_wp

    z_v1 = divide_icos(ico%vertex0%vertex,ico%vertex1%vertex,z_ratio)
    z_v2 = divide_icos(ico%vertex0%vertex,ico%vertex2%vertex,z_ratio)
    z_v0 = ico%vertex0%vertex
    z_center = circum_center(z_v0,z_v1,z_v2)
    z_sl0 = z_triangle_area-triangle_area(z_v0,z_v1,z_v2)
    z_cgal0 = cos_arc_length(z_v1,z_v2)
    z_crcl0 = 0.0_wp
    CALL solve_small_circle(z_sl0,z_cgal0,z_crcl0)
    z_ed0%cos_rel_co_lat = z_crcl0
    z_ed0%vertex0%vertex = z_v1
    z_ed0%vertex1%vertex = z_v2
    CALL compute_edge(z_sl0,z_cgal0,z_ed0,z_center)

    z_v0 = divide_icos(ico%vertex1%vertex,ico%vertex0%vertex,z_ratio)
    z_v2 = divide_icos(ico%vertex1%vertex,ico%vertex2%vertex,z_ratio)
    z_v1 = ico%vertex1%vertex
    z_center = circum_center(z_v0,z_v1,z_v2)
    z_ed1%vertex0%vertex = z_v0
    z_ed1%vertex1%vertex = z_v2
    z_ed1%cos_rel_co_lat = z_crcl0
    CALL compute_edge(z_sl0,z_cgal0,z_ed1,z_center)

    z_v0 = divide_icos(ico%vertex2%vertex,ico%vertex0%vertex,z_ratio)
    z_v1 = divide_icos(ico%vertex2%vertex,ico%vertex1%vertex,z_ratio)
    z_v2 = ico%vertex2%vertex
    z_center = circum_center(z_v0,z_v1,z_v2)
    z_ed2%vertex0%vertex = z_v0
    z_ed2%vertex1%vertex = z_v1
    z_ed2%cos_rel_co_lat = z_crcl0
    CALL compute_edge(z_sl0,z_cgal0,z_ed2,z_center)


    !------------------------------------------
    ! Section2: Solve for the inner 4 triangles
    !------------------------------------------

    z_v0 = z_ed0%edge_center
    z_v1 = z_ed1%edge_center
    z_v2 = z_ed2%edge_center
    z_center = circum_center(z_v0,z_v1,z_v2)
    z_sl1 = (z_triangle_area-triangle_area(z_v0,z_v1,z_v2))/3.0_wp
    z_cgal1 = cos_arc_length(z_v0,z_v1)
    z_crcl1 = 0.0_wp
    CALL solve_small_circle(z_sl1,z_cgal1,z_crcl1)
    z_ed3%vertex0%vertex = z_v0
    z_ed3%vertex1%vertex = z_v1
    z_ed3%cos_rel_co_lat = z_crcl1
    CALL compute_edge(z_sl1,z_cgal1,z_ed3,z_center)
    z_ed4%vertex0%vertex = z_v1
    z_ed4%vertex1%vertex = z_v2
    z_ed4%cos_rel_co_lat = z_crcl1
    CALL compute_edge(z_sl1,z_cgal1,z_ed4,z_center)
    z_ed5%vertex0%vertex = z_v2
    z_ed5%vertex1%vertex = z_v0
    z_ed5%cos_rel_co_lat = z_crcl1
    CALL compute_edge(z_sl1,z_cgal1,z_ed5,z_center)

    !-----------------------------------------------------------
    ! Section 3: Determine vertices of all triangles explicitely
    !-----------------------------------------------------------

    z_ratio=0.2_wp

    root(0,idx)%t%vertex0%vertex = ico%vertex0%vertex
    root(0,idx)%t%vertex1%vertex = divide_icos(ico%vertex0%vertex,&
      & ico%vertex1%vertex,z_ratio)
    root(0,idx)%t%vertex2%vertex = divide_icos(ico%vertex0%vertex,&
      & ico%vertex2%vertex,z_ratio)
    root(5,idx)%t%vertex0%vertex = z_ed3%edge_center
    root(5,idx)%t%vertex1%vertex = z_ed0%edge_center
    root(5,idx)%t%vertex2%vertex = z_ed0%vertex0%vertex
    root(8,idx)%t%vertex0%vertex = z_ed0%vertex1%vertex
    root(8,idx)%t%vertex1%vertex = z_ed5%edge_center
    root(8,idx)%t%vertex2%vertex = z_ed2%vertex0%vertex
    root(9,idx)%t%vertex0%vertex = z_ed1%vertex0%vertex
    root(9,idx)%t%vertex1%vertex = divide_icos(ico%vertex1%vertex,&
      & ico%vertex0%vertex,z_ratio)
    root(9,idx)%t%vertex2%vertex = z_ed1%edge_center
    root(20,idx)%t%vertex0%vertex = z_ed4%edge_center
    root(20,idx)%t%vertex1%vertex = z_ed1%vertex1%vertex
    root(20,idx)%t%vertex2%vertex = z_ed2%vertex1%vertex
    root(23,idx)%t%vertex0%vertex = divide_icos(ico%vertex2%vertex,&
      & ico%vertex1%vertex,z_ratio)
    root(23,idx)%t%vertex1%vertex = divide_icos(ico%vertex2%vertex,&
      & ico%vertex0%vertex,z_ratio)
    root(23,idx)%t%vertex2%vertex = z_ed2%edge_center
    root(16,idx)%t%vertex1%vertex = ico%vertex1%vertex
    root(16,idx)%t%vertex2%vertex = divide_icos(ico%vertex1%vertex,&
      & ico%vertex2%vertex,z_ratio)
    root(24,idx)%t%vertex2%vertex = ico%vertex2%vertex

    !------------------------------------------------
    ! Section 4: Determine triangle centers and areas
    !------------------------------------------------

    DO j=0,24
      root(j,idx)%t%triangle_area=0.04_wp*ico%triangle_area
      root(j,idx)%t%triangle_center=circum_center(root(j,idx)%t%vertex0%vertex,&
        & root(j,idx)%t%vertex1%vertex,&
        & root(j,idx)%t%vertex2%vertex)
    ENDDO

    !--------------------------------------------------------------
    ! Section 5: Small circle features of prorated icosahedron arcs
    !--------------------------------------------------------------

    DO j = 0,4
      !left edge of icosahedron
      ij=j*j
      root(ij,idx)%t%edge0%cos_rel_co_lat=0.0_wp
      root(ij,idx)%t%edge0%sc_pole%x=0.0_wp
      root(ij,idx)%t%edge0%edge_center=&
        & edge_midpoint(root(ij,idx)%t%vertex0%vertex,&
        & root(ij,idx)%t%vertex1%vertex,&
        & root(ij,idx)%t%edge0%cos_rel_co_lat,&
        & root(ij,idx)%t%edge0%sc_pole)
      root(ij,idx)%t%edge0%edge_primal_arc=&
        & arc_length(root(ij,idx)%t%vertex0%vertex,&
        & root(ij,idx)%t%vertex1%vertex)
      !right edge of icosahedron
      ij=(j+1)*(j+1)-1
      root(ij,idx)%t%edge2%cos_rel_co_lat=0.0_wp
      root(ij,idx)%t%edge2%sc_pole%x=0.0_wp
      root(ij,idx)%t%edge2%edge_center=&
        & edge_midpoint(root(ij,idx)%t%vertex0%vertex,&
        & root(ij,idx)%t%vertex2%vertex,&
        & root(ij,idx)%t%edge2%cos_rel_co_lat,&
        & root(ij,idx)%t%edge2%sc_pole)
      root(ij,idx)%t%edge2%edge_primal_arc=&
        & arc_length(root(ij,idx)%t%vertex0%vertex,&
        & root(ij,idx)%t%vertex2%vertex)
      !lower edge of icosahedron
      ij=16+2*j
      root(ij,idx)%t%edge1%cos_rel_co_lat=0.0_wp
      root(ij,idx)%t%edge1%sc_pole%x=0.0_wp
      root(ij,idx)%t%edge1%edge_center=&
        & edge_midpoint(root(ij,idx)%t%vertex1%vertex,&
        & root(ij,idx)%t%vertex2%vertex,&
        & root(ij,idx)%t%edge1%cos_rel_co_lat,&
        & root(ij,idx)%t%edge1%sc_pole)
      root(ij,idx)%t%edge1%edge_primal_arc=&
        & arc_length(root(ij,idx)%t%vertex1%vertex,&
        & root(ij,idx)%t%vertex2%vertex)
    ENDDO

    !-----------------------------------------------------------------------
    ! Section 6: Already known small circle features from fourfold triangles
    !-----------------------------------------------------------------------

    ! Refering to edge 3
    root(10,idx)%t%edge0%cos_rel_co_lat  = z_ed3%cos_rel_co_lat
    root(10,idx)%t%edge0%sc_pole         = z_ed3%sc_pole
    root(10,idx)%t%edge0%edge_primal_arc = 0.5_wp*z_ed3%edge_primal_arc
    root(10,idx)%t%edge0%edge_center     = &
      & edge_midpoint(root(10,idx)%t%vertex0%vertex,&
      & root(10,idx)%t%vertex1%vertex,&
      & root(10,idx)%t%edge0%cos_rel_co_lat,&
      & root(10,idx)%t%edge0%sc_pole)
    root( 5,idx)%t%edge0%cos_rel_co_lat  = z_ed3%cos_rel_co_lat
    root( 5,idx)%t%edge0%sc_pole         = z_ed3%sc_pole
    root( 5,idx)%t%edge0%edge_primal_arc = 0.5_wp*z_ed3%edge_primal_arc
    root( 5,idx)%t%edge0%edge_center     = &
      & edge_midpoint(root( 5,idx)%t%vertex0%vertex,&
      & root( 5,idx)%t%vertex1%vertex,&
      & root( 5,idx)%t%edge0%cos_rel_co_lat,&
      & root( 5,idx)%t%edge0%sc_pole)


    ! Refering to edge 0
    root( 5,idx)%t%edge1%cos_rel_co_lat  = z_ed0%cos_rel_co_lat
    root( 5,idx)%t%edge1%sc_pole         = z_ed0%sc_pole
    root( 5,idx)%t%edge1%edge_primal_arc = 0.5_wp*z_ed0%edge_primal_arc
    root( 5,idx)%t%edge1%edge_center     = &
      & edge_midpoint(root( 5,idx)%t%vertex1%vertex,&
      & root( 5,idx)%t%vertex2%vertex,&
      & root( 5,idx)%t%edge1%cos_rel_co_lat,&
      & root( 5,idx)%t%edge1%sc_pole)
    root( 7,idx)%t%edge1%cos_rel_co_lat  = z_ed0%cos_rel_co_lat
    root( 7,idx)%t%edge1%sc_pole         = z_ed0%sc_pole
    root( 7,idx)%t%edge1%edge_primal_arc = 0.5_wp*z_ed0%edge_primal_arc
    root( 7,idx)%t%edge1%edge_center     = &
      & edge_midpoint(root( 7,idx)%t%vertex1%vertex,&
      & root( 7,idx)%t%vertex2%vertex,&
      & root( 7,idx)%t%edge1%cos_rel_co_lat,&
      & root( 7,idx)%t%edge1%sc_pole)

    ! Refering to edge 5
    root( 7,idx)%t%edge2%cos_rel_co_lat  = z_ed5%cos_rel_co_lat
    root( 7,idx)%t%edge2%sc_pole         = z_ed5%sc_pole
    root( 7,idx)%t%edge2%edge_primal_arc = 0.5_wp*z_ed5%edge_primal_arc
    root( 7,idx)%t%edge2%edge_center     = &
      & edge_midpoint(root( 7,idx)%t%vertex2%vertex,&
      & root( 7,idx)%t%vertex0%vertex,&
      & root( 7,idx)%t%edge2%cos_rel_co_lat,&
      & root( 7,idx)%t%edge2%sc_pole)
    root(14,idx)%t%edge2%cos_rel_co_lat  = z_ed5%cos_rel_co_lat
    root(14,idx)%t%edge2%sc_pole         = z_ed5%sc_pole
    root(14,idx)%t%edge2%edge_primal_arc = 0.5_wp*z_ed5%edge_primal_arc
    root(14,idx)%t%edge2%edge_center     = &
      & edge_midpoint(root(14,idx)%t%vertex2%vertex,&
      & root(14,idx)%t%vertex0%vertex,&
      & root(14,idx)%t%edge2%cos_rel_co_lat,&
      & root(14,idx)%t%edge2%sc_pole)

    ! Refering to edge 2
    root(14,idx)%t%edge0%cos_rel_co_lat  = z_ed2%cos_rel_co_lat
    root(14,idx)%t%edge0%sc_pole         = z_ed2%sc_pole
    root(14,idx)%t%edge0%edge_primal_arc = 0.5_wp*z_ed2%edge_primal_arc
    root(14,idx)%t%edge0%edge_center     = &
      & edge_midpoint(root(14,idx)%t%vertex0%vertex,&
      & root(14,idx)%t%vertex1%vertex,&
      & root(14,idx)%t%edge0%cos_rel_co_lat,&
      & root(14,idx)%t%edge0%sc_pole)
    root(21,idx)%t%edge0%cos_rel_co_lat  = z_ed2%cos_rel_co_lat
    root(21,idx)%t%edge0%sc_pole         = z_ed2%sc_pole
    root(21,idx)%t%edge0%edge_primal_arc = 0.5_wp*z_ed2%edge_primal_arc
    root(21,idx)%t%edge0%edge_center     = &
      & edge_midpoint(root(21,idx)%t%vertex0%vertex,&
      & root(21,idx)%t%vertex1%vertex,&
      & root(21,idx)%t%edge0%cos_rel_co_lat,&
      & root(21,idx)%t%edge0%sc_pole)

    ! Refering to edge 4
    root(21,idx)%t%edge1%cos_rel_co_lat  = z_ed4%cos_rel_co_lat
    root(21,idx)%t%edge1%sc_pole         = z_ed4%sc_pole
    root(21,idx)%t%edge1%edge_primal_arc = 0.5_wp*z_ed4%edge_primal_arc
    root(21,idx)%t%edge1%edge_center     = &
      & edge_midpoint(root(21,idx)%t%vertex1%vertex,&
      & root(21,idx)%t%vertex2%vertex,&
      & root(21,idx)%t%edge1%cos_rel_co_lat,&
      & root(21,idx)%t%edge1%sc_pole)
    root(19,idx)%t%edge1%cos_rel_co_lat  = z_ed4%cos_rel_co_lat
    root(19,idx)%t%edge1%sc_pole         = z_ed4%sc_pole
    root(19,idx)%t%edge1%edge_primal_arc = 0.5_wp*z_ed4%edge_primal_arc
    root(19,idx)%t%edge1%edge_center     = &
      & edge_midpoint(root(19,idx)%t%vertex1%vertex,&
      & root(19,idx)%t%vertex2%vertex,&
      & root(19,idx)%t%edge1%cos_rel_co_lat,&
      & root(19,idx)%t%edge1%sc_pole)

    ! Refering to edge 1
    root(19,idx)%t%edge2%cos_rel_co_lat  = z_ed1%cos_rel_co_lat
    root(19,idx)%t%edge2%sc_pole         = z_ed1%sc_pole
    root(19,idx)%t%edge2%edge_primal_arc = 0.5_wp*z_ed1%edge_primal_arc
    root(19,idx)%t%edge2%edge_center     = &
      & edge_midpoint(root(19,idx)%t%vertex2%vertex,&
      & root(19,idx)%t%vertex0%vertex,&
      & root(19,idx)%t%edge2%cos_rel_co_lat,&
      & root(19,idx)%t%edge2%sc_pole)
    root(10,idx)%t%edge2%cos_rel_co_lat  = z_ed1%cos_rel_co_lat
    root(10,idx)%t%edge2%sc_pole         = z_ed1%sc_pole
    root(10,idx)%t%edge2%edge_primal_arc = 0.5_wp*z_ed1%edge_primal_arc
    root(10,idx)%t%edge2%edge_center     = &
      & edge_midpoint(root(10,idx)%t%vertex2%vertex,&
      & root(10,idx)%t%vertex0%vertex,&
      & root(10,idx)%t%edge2%cos_rel_co_lat,&
      & root(10,idx)%t%edge2%sc_pole)

    !------------------------------------------------
    ! Section 7: Solve for remaining  prototype edges
    !------------------------------------------------

    z_sl0 = root(0,idx)%t%triangle_area-triangle_area(&
      & root(0,idx)%t%vertex0%vertex,&
      & root(0,idx)%t%vertex1%vertex,&
      & root(0,idx)%t%vertex2%vertex)
    z_cgal0 = cos_arc_length(root(0,idx)%t%vertex1%vertex,&
      & root(0,idx)%t%vertex2%vertex)
    CALL solve_small_circle(z_sl0,z_cgal0,z_crcl0)
    z_ed0%vertex0%vertex = root(0,idx)%t%vertex1%vertex
    z_ed0%vertex1%vertex = root(0,idx)%t%vertex2%vertex
    z_ed0%cos_rel_co_lat = z_crcl0
    CALL compute_edge(z_sl0,z_cgal0,z_ed0,root(0,idx)%t%triangle_center)

    z_sl1 =(root(2,idx)%t%triangle_area-triangle_area(&
      & root(2,idx)%t%vertex0%vertex,&
      & root(2,idx)%t%vertex1%vertex,&
      & root(2,idx)%t%vertex2%vertex) - (-z_sl0))*0.5_wp
    z_cgal1 = cos_arc_length(root(2,idx)%t%vertex0%vertex,&
      & root(2,idx)%t%vertex1%vertex)
    CALL solve_small_circle(z_sl1,z_cgal1,z_crcl1)
    z_ed1%vertex0%vertex = root(0,idx)%t%vertex0%vertex
    z_ed1%vertex1%vertex = root(0,idx)%t%vertex1%vertex
    z_ed1%cos_rel_co_lat = z_crcl1
    CALL compute_edge(z_sl1,z_cgal1,z_ed1,root(2,idx)%t%triangle_center)

    z_sl2 =(root(4,idx)%t%triangle_area-triangle_area(&
      & root(4,idx)%t%vertex0%vertex,&
      & root(4,idx)%t%vertex1%vertex,&
      & root(4,idx)%t%vertex2%vertex))*0.5_wp
    z_cgal2 = cos_arc_length(root(4,idx)%t%vertex0%vertex,&
      & root(4,idx)%t%vertex2%vertex)
    z_crcl2=0.0_wp
    CALL solve_small_circle(z_sl2,z_cgal2,z_crcl2)
    z_ed2%vertex0%vertex = root(4,idx)%t%vertex0%vertex
    z_ed2%vertex1%vertex = root(4,idx)%t%vertex2%vertex
    z_ed2%cos_rel_co_lat = z_crcl2
    CALL compute_edge(z_sl2,z_cgal2,z_ed2,root(4,idx)%t%triangle_center)

    z_sl3 =(root(12,idx)%t%triangle_area-triangle_area(&
      & root(12,idx)%t%vertex0%vertex,&
      & root(12,idx)%t%vertex1%vertex,&
      & root(12,idx)%t%vertex2%vertex))/3.0_wp
    z_cgal3 = cos_arc_length(root(12,idx)%t%vertex0%vertex,&
      & root(12,idx)%t%vertex2%vertex)
    z_crcl3=0.0_wp
    CALL solve_small_circle(z_sl3,z_cgal3,z_crcl3)
    z_ed3%vertex0%vertex = root(12,idx)%t%vertex0%vertex
    z_ed3%vertex1%vertex = root(12,idx)%t%vertex2%vertex
    z_ed3%cos_rel_co_lat = z_crcl3
    CALL compute_edge(z_sl3,z_cgal3,z_ed3,root(12,idx)%t%triangle_center)


    !----------------------------------------------------------------
    ! Section 8: Assign actual values of cos_rel_co_lats to the edges
    !----------------------------------------------------------------

    root(12,idx)%t%edge0%cos_rel_co_lat = z_ed3%cos_rel_co_lat
    root(12,idx)%t%edge1%cos_rel_co_lat = z_ed3%cos_rel_co_lat
    root(12,idx)%t%edge2%cos_rel_co_lat = z_ed3%cos_rel_co_lat

    root( 2,idx)%t%edge1%cos_rel_co_lat = z_ed0%cos_rel_co_lat
    root( 2,idx)%t%edge0%cos_rel_co_lat = z_ed1%cos_rel_co_lat
    root( 2,idx)%t%edge2%cos_rel_co_lat = z_ed1%cos_rel_co_lat

    root(17,idx)%t%edge2%cos_rel_co_lat = z_ed0%cos_rel_co_lat
    root(17,idx)%t%edge0%cos_rel_co_lat = z_ed1%cos_rel_co_lat
    root(17,idx)%t%edge1%cos_rel_co_lat = z_ed1%cos_rel_co_lat

    root(23,idx)%t%edge0%cos_rel_co_lat = z_ed0%cos_rel_co_lat
    root(23,idx)%t%edge2%cos_rel_co_lat = z_ed1%cos_rel_co_lat
    root(23,idx)%t%edge1%cos_rel_co_lat = z_ed1%cos_rel_co_lat

    root( 4,idx)%t%edge1%cos_rel_co_lat = z_ed2%cos_rel_co_lat
    root( 4,idx)%t%edge2%cos_rel_co_lat = z_ed2%cos_rel_co_lat

    root( 8,idx)%t%edge1%cos_rel_co_lat = z_ed2%cos_rel_co_lat
    root( 8,idx)%t%edge0%cos_rel_co_lat = z_ed2%cos_rel_co_lat

    root(20,idx)%t%edge2%cos_rel_co_lat = z_ed2%cos_rel_co_lat
    root(20,idx)%t%edge0%cos_rel_co_lat = z_ed2%cos_rel_co_lat

    !-------------------------------------------------------------------------
    ! Section 9: Compute the edge properties (Be careful with the signs of the
    !            semi-lunes.)
    !-------------------------------------------------------------------------

    CALL compute_edge(z_sl3,z_cgal3,root(12,idx)%t%edge0,&
      & root(12,idx)%t%triangle_center)
    CALL compute_edge(z_sl3,z_cgal3,root(12,idx)%t%edge1,&
      & root(12,idx)%t%triangle_center)
    CALL compute_edge(z_sl3,z_cgal3,root(12,idx)%t%edge2,&
      & root(12,idx)%t%triangle_center)

    CALL compute_edge(-z_sl0,z_cgal0,root( 2,idx)%t%edge1,&
      & root( 2,idx)%t%triangle_center)
    CALL compute_edge(z_sl1,z_cgal1,root( 2,idx)%t%edge0,&
      & root( 2,idx)%t%triangle_center)
    CALL compute_edge(z_sl1,z_cgal1,root( 2,idx)%t%edge2,&
      & root( 2,idx)%t%triangle_center)

    CALL compute_edge(-z_sl0,z_cgal0,root(17,idx)%t%edge2,&
      & root(17,idx)%t%triangle_center)
    CALL compute_edge(z_sl1,z_cgal1,root(17,idx)%t%edge0,&
      & root(17,idx)%t%triangle_center)
    CALL compute_edge(z_sl1,z_cgal1,root(17,idx)%t%edge1,&
      & root(17,idx)%t%triangle_center)


    CALL compute_edge(-z_sl0,z_cgal0,root(23,idx)%t%edge0,&
      & root(23,idx)%t%triangle_center)
    CALL compute_edge(z_sl1,z_cgal1,root(23,idx)%t%edge1,&
      & root(23,idx)%t%triangle_center)
    CALL compute_edge(z_sl1,z_cgal1,root(23,idx)%t%edge2,&
      & root(23,idx)%t%triangle_center)

    CALL compute_edge(z_sl2,z_cgal2,root( 4,idx)%t%edge1,&
      & root( 4,idx)%t%triangle_center)
    CALL compute_edge(z_sl2,z_cgal2,root( 4,idx)%t%edge2,&
      & root( 4,idx)%t%triangle_center)

    CALL compute_edge(z_sl2,z_cgal2,root( 8,idx)%t%edge1,&
      & root( 8,idx)%t%triangle_center)
    CALL compute_edge(z_sl2,z_cgal2,root( 8,idx)%t%edge0,&
      & root( 8,idx)%t%triangle_center)

    CALL compute_edge(z_sl2,z_cgal2,root(20,idx)%t%edge2,&
      & root(20,idx)%t%triangle_center)
    CALL compute_edge(z_sl2,z_cgal2,root(20,idx)%t%edge0,&
      & root(20,idx)%t%triangle_center)

    ! Section 10: Deallocate dummy edge vertices
    iastat=0
    DEALLOCATE(z_ed0%vertex0,stat=istat)
    iastat=ABS(istat)+iastat
    DEALLOCATE(z_ed0%vertex1,stat=istat)
    iastat=ABS(istat)+iastat
    DEALLOCATE(z_ed1%vertex0,stat=istat)
    iastat=ABS(istat)+iastat
    DEALLOCATE(z_ed1%vertex1,stat=istat)
    iastat=ABS(istat)+iastat
    DEALLOCATE(z_ed2%vertex0,stat=istat)
    iastat=ABS(istat)+iastat
    DEALLOCATE(z_ed2%vertex1,stat=istat)
    iastat=ABS(istat)+iastat
    DEALLOCATE(z_ed3%vertex0,stat=istat)
    iastat=ABS(istat)+iastat
    DEALLOCATE(z_ed3%vertex1,stat=istat)
    iastat=ABS(istat)+iastat
    DEALLOCATE(z_ed4%vertex0,stat=istat)
    iastat=ABS(istat)+iastat
    DEALLOCATE(z_ed4%vertex1,stat=istat)
    iastat=ABS(istat)+iastat
    DEALLOCATE(z_ed5%vertex0,stat=istat)
    iastat=ABS(istat)+iastat
    DEALLOCATE(z_ed5%vertex1,stat=istat)
    IF( iastat /= 0 ) THEN
      WRITE(message_text,'(a)') 'Problem in deallocating dummy vertices'
      CALL finish('equal_area_twen5',message_text)
    ENDIF

  END SUBROUTINE equal_area_twen5
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  !
  !
  !>
  !! Gives the point on an icosahedron arc that divides the arc by "ratio".
  !!
  !! Gives the point on an icosahedron arc that divides the arc by "ratio"
  !! looking from v0 to v1.
  !!
  !! @par Revision History
  !! Initial release by Almut Gassmann (2007-03-22)
  !!
  FUNCTION divide_icos(v0,v1,ratio) result(vd)
    !
    TYPE (t_cartesian_coordinates), INTENT(in) :: v0, v1 ! vertices of icosahedron
    REAL (wp), INTENT(in) :: ratio  ! dividing ratio
    !
    TYPE (t_cartesian_coordinates) :: vd  ! point on the icos. arc
    !
    TYPE (t_cartesian_coordinates) :: z_help1, z_help2
    REAL (wp) :: z_ico_alen, z_tan_arlen
    !
    !-----------------------------------------------------------------------------

    z_help1%x    = v0%x+v1%x
    z_help2%x    = v0%x-v1%x
    z_ico_alen   = 2.0_wp*acos(1.0_wp/(2.0_wp*sin(0.2_wp*pi)))
    z_tan_arlen  = TAN((0.5_wp-ratio)*z_ico_alen)
    vd%x = z_help1%x/SQRT(DOT_PRODUCT(z_help1%x,z_help1%x))+&
      & z_help2%x/SQRT(DOT_PRODUCT(z_help2%x,z_help2%x))*z_tan_arlen
    vd%x = vd%x / SQRT(DOT_PRODUCT(vd%x,vd%x))

  END FUNCTION divide_icos

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------

END MODULE mo_equal_area
!
