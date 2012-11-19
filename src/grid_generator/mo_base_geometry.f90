!>
!!  Contains the definition of basic used defined data types.
!!
!!  Contains the definition of basic used defined data types
!!  to store the gridpoint coordinates either in cartesian or
!!  spherical coordinates. It also contains a number of (mostly
!!  ELEMENTAL) functions used to compute vectors and geometric quantities
!!  needed by the grid generator.
!!
!! @par Revision History
!!  Initial version  by Luis Kornblueh (2004)
!!  Modified to include tangent vectors and Protex headers
!!  by Luca Bonaventura (2005)
!! @par
!!  Guenther Zaengl, DWD, 2008-10-10:
!!  Add subroutine gvec2cvec (copied from mo_math_utilities)
!!
!! @par Copyright
!! 2002-2007 by DWD and MPI-M
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
MODULE mo_base_geometry

  USE mo_kind, ONLY: wp
  USE mo_math_constants, ONLY: pi, pi_2, dbl_eps, eps
  USE mo_exception, ONLY: finish

  IMPLICIT NONE
#include "grid_definitions.inc"

  PRIVATE

  PUBLIC :: t_cartesian_coordinates
  PUBLIC :: t_coordinates_reference
  PUBLIC :: t_geographical_coordinates
  PUBLIC :: t_tangent_vectors
  PUBLIC :: sphere_cartesian_midpoint
  PUBLIC :: norma, normalize, angle_of_vectors, sin_cc, sphere_tanget_coordinates
  PUBLIC :: angle_of_points, arc_length_normalsphere
  PUBLIC :: cc2gc, gc2cc, cc2tv, gvec2cvec
  PUBLIC :: triangle_area, circum_center, arc_length, cos_arc_length, bary_center
  PUBLIC :: normal_vector, inter_section, inter_section2, vector_product
  PUBLIC :: rotate_z,rotate_x,rotate_y
  PUBLIC :: norma_of_vector_product
  PUBLIC :: x_rot_angle,y_rot_angle,z_rot_angle
  PUBLIC :: middle
  PUBLIC :: cartesian_to_geographical, geographical_to_cartesian
  PUBLIC :: project_point
  PUBLIC :: operator(+), operator(-)

  TYPE t_cartesian_coordinates
    REAL(wp) :: x(3)
  END TYPE t_cartesian_coordinates

  TYPE t_coordinates_reference
    TYPE(t_cartesian_coordinates), POINTER :: px
  END TYPE t_coordinates_reference

  TYPE t_geographical_coordinates
    REAL(wp) :: lon
    REAL(wp) :: lat
  END TYPE t_geographical_coordinates

  TYPE t_tangent_vectors
    REAL(wp) :: v1
    REAL(wp) :: v2
  END TYPE t_tangent_vectors

  REAL(wp) :: x_rot_angle,y_rot_angle,z_rot_angle

  INTERFACE OPERATOR(+)
    MODULE PROCEDURE add_cartesian
  END INTERFACE
  INTERFACE OPERATOR(-)
    MODULE PROCEDURE subtract_cartesian
  END INTERFACE

CONTAINS


  !--------------------------------------------------------------------
  ELEMENTAL FUNCTION add_cartesian(v1, v2) result(v)
    TYPE(t_cartesian_coordinates), INTENT(in) :: v1, v2
    TYPE(t_cartesian_coordinates)  :: v
    v%x = v1%x + v2%x
  END FUNCTION add_cartesian

  ELEMENTAL FUNCTION subtract_cartesian(v1, v2) result(v)
    TYPE(t_cartesian_coordinates), INTENT(in) :: v1, v2
    TYPE(t_cartesian_coordinates)  :: v
    v%x = v1%x - v2%x
  END FUNCTION subtract_cartesian

  !--------------------------------------------------------------------
  !>
  ELEMENTAL FUNCTION middle(p1,p2)

    TYPE(t_cartesian_coordinates), INTENT(in) :: p1,p2
    TYPE(t_cartesian_coordinates) :: middle

    middle%x = (p1%x + p2%x) / 2.0_wp
  
  END FUNCTION middle
  !--------------------------------------------------------------------

  !--------------------------------------------------------------------
  !>
  !! Converts cartesian coordinates to geographical.
  !!
  !! @par Revision History
  !! Developed  by Luis Kornblueh  (2004).
  !! Completely new version by Thomas Heinze (2006-07-20)
  !!
  ELEMENTAL FUNCTION cc2gc(x) result (position)

    TYPE(t_cartesian_coordinates), INTENT(in) :: x
    TYPE(t_geographical_coordinates)          :: position

    REAL(wp)                                :: z_x, z_y, z_z, z_r
    !-----------------------------------------------------------------------

    z_x = x%x(1)
    z_y = x%x(2)
    z_z = x%x(3)

    z_r = z_x * z_x + z_y * z_y
    z_r = SQRT(z_r)

    IF (ABS(z_r) < dbl_eps) THEN    ! one of the poles

      IF (z_z > 0.0_wp) THEN
        position%lat = pi_2
      ELSE
        position%lat = -1._wp * pi_2
      END IF
      position%lon = 0._wp

    ELSE

      position%lat = ATAN2 ( z_z, z_r)

      IF (ABS(z_x) < dbl_eps) THEN    ! z_x == 0 ?

        IF (z_y >= 0.0_wp) THEN
          position%lon = pi_2
        ELSE
          position%lon = -1._wp * pi_2
        END IF

      ELSE
        position%lon = ATAN2( z_y, z_x)
      END IF

    END IF

  END FUNCTION cc2gc
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  SUBROUTINE cartesian_to_geographical(cartesian, no_of_points, geo_coordinates, &
    & start_point, end_point )
    TYPE(t_cartesian_coordinates), POINTER :: cartesian(:)
    TYPE(t_geographical_coordinates), POINTER :: geo_coordinates(:)
    INTEGER, INTENT(in) :: no_of_points
    INTEGER, INTENT(in), OPTIONAL :: start_point, end_point

    INTEGER :: i, start_p, end_p

    IF (PRESENT(start_point)) THEN
      start_p = start_point
      end_p   = end_point
    ELSE
      start_p = 1
      end_p   = no_of_points
    ENDIF
    IF (.NOT. ASSOCIATED(geo_coordinates)) THEN
      ALLOCATE (geo_coordinates(start_p:end_p), stat=i)
      IF (i >0) THEN
        CALL finish ('cartesian_to_geographical', 'Problem in allocating local arrays')
      ENDIF
    ENDIF
    
!$OMP PARALLEL 
!$OMP DO PRIVATE(i)
    DO i=start_p, end_p
        geo_coordinates(i) = cc2gc(cartesian(i))
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
    
  END SUBROUTINE cartesian_to_geographical
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE geographical_to_cartesian(geo_coordinates, no_of_points, cartesian)
    TYPE(t_geographical_coordinates), POINTER :: geo_coordinates(:)
    TYPE(t_cartesian_coordinates), POINTER :: cartesian(:)
    INTEGER, INTENT(in) :: no_of_points

    INTEGER :: i
    REAL(wp) :: cos_lat
    REAL(wp) :: check

    IF (.NOT. ASSOCIATED(cartesian)) THEN
      ALLOCATE (cartesian(no_of_points), stat=i)
      IF (i >0) THEN
        CALL finish ('geographical_to_cartesian', 'Problem in allocating local arrays')
      ENDIF
    ENDIF


!$OMP PARALLEL 
!$OMP DO PRIVATE(i, cos_lat, check)
    DO i=1,no_of_points
      cos_lat           = COS(geo_coordinates(i)%lat)
      cartesian(i)%x(1) = COS(geo_coordinates(i)%lon) * cos_lat
      cartesian(i)%x(2) = SIN(geo_coordinates(i)%lon) * cos_lat
      cartesian(i)%x(3) = SIN(geo_coordinates(i)%lat)

      check = d_norma_3d(cartesian(i))
      IF (ABS(check-1.0_wp) > 1.0e-8_wp) &
        & CALL finish ('geographicalToCartesian', 'Error in calculation')

      !    WRITE(*,*) i, " LonLat=", geoCoordinates(i)
      !    WRITE(*,*) i, " Cart=", cartesian(i)

    ENDDO ! i=1,noOfPoints
!$OMP END DO
!$OMP END PARALLEL

    RETURN

  END SUBROUTINE geographical_to_cartesian
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Converts longitude and latitude to cartesian coordinates.
  !!
  !! @par Revision History
  !! Developed  by Luis Kornblueh  (2004).
  ELEMENTAL FUNCTION gc2cc(position) result(x)
    !-----------------------------------------------------------------------

    TYPE(t_geographical_coordinates), INTENT(in) :: position
    TYPE(t_cartesian_coordinates) :: x

    REAL (wp) :: z_cln, z_sln, z_clt, z_slt

    z_sln = SIN(position%lon)
    z_cln = COS(position%lon)
    z_slt = SIN(position%lat)
    z_clt = COS(position%lat)

    x%x(1) = z_cln*z_clt
    x%x(2) = z_sln*z_clt
    x%x(3) = z_slt

  END FUNCTION gc2cc
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Converts  vectors (in)  cartesian coordinate representation
  !! to tangent vectors.
  !!
  !! @par Revision History
  !! Developed  by Luca Bonaventura  (2005).
  FUNCTION cc2tv(xx,position) result (tt)
    !-----------------------------------------------------------------------


    TYPE(t_cartesian_coordinates), INTENT(in) :: xx
    TYPE(t_geographical_coordinates),INTENT(in)  :: position
    TYPE(t_tangent_vectors):: tt

    REAL(wp) :: z_sinlo,z_coslo,z_sinlacoslo,z_sinlasinlo

    z_sinlo=SIN(position%lon)
    z_coslo=COS(position%lon)
    z_sinlacoslo=SIN(position%lat)*z_coslo
    z_sinlasinlo=SIN(position%lat)*z_sinlo
    tt%v1= -z_sinlo*xx%x(1)+z_coslo*xx%x(2)
    tt%v2= -z_sinlacoslo*xx%x(1)-z_sinlasinlo*xx%x(2) + COS(position%lat)*xx%x(3)

  END FUNCTION cc2tv
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  !! Computes area of triangular cell.
  !!
  !! @par Revision History
  !! Developed  by Luis Kornblueh  (2004).
  ELEMENTAL FUNCTION triangle_area (x0, x1, x2) result(area)

    TYPE(t_cartesian_coordinates), INTENT(in) :: x0, x1, x2

    REAL(wp) :: area
    REAL(wp) :: z_s12, z_s23, z_s31, z_ca1, z_ca2, z_ca3, z_a1, z_a2, z_a3

    TYPE(t_cartesian_coordinates) :: u12, u23, u31

    ! This variant to calculate the area of a spherical triangle
    ! is more precise.

    !  Compute cross products Uij = Vi x Vj.
    u12 = vector_product (x0, x1)
    u23 = vector_product (x1, x2)
    u31 = vector_product (x2, x0)

    !  Normalize Uij to unit vectors.
    z_s12 = DOT_PRODUCT ( u12%x(1:3), u12%x(1:3) )
    z_s23 = DOT_PRODUCT ( u23%x(1:3), u23%x(1:3) )
    z_s31 = DOT_PRODUCT ( u31%x(1:3), u31%x(1:3) )

    !  Test for a degenerate triangle associated with collinear vertices.
    IF ( z_s12 == 0.0_wp .or. z_s23 == 0.0_wp  .or. z_s31 == 0.0_wp ) THEN
      area = 0.0_wp
      RETURN
    END IF

    z_s12 = SQRT(z_s12)
    z_s23 = SQRT(z_s23)
    z_s31 = SQRT(z_s31)

    u12%x(1:3) = u12%x(1:3)/z_s12
    u23%x(1:3) = u23%x(1:3)/z_s23
    u31%x(1:3) = u31%x(1:3)/z_s31

    !  Compute interior angles Ai as the dihedral angles between planes:
    !  CA1 = cos(A1) = -<U12,U31>
    !  CA2 = cos(A2) = -<U23,U12>
    !  CA3 = cos(A3) = -<U31,U23>
    z_ca1 = -u12%x(1)*u31%x(1)-u12%x(2)*u31%x(2)-u12%x(3)*u31%x(3)
    z_ca2 = -u23%x(1)*u12%x(1)-u23%x(2)*u12%x(2)-u23%x(3)*u12%x(3)
    z_ca3 = -u31%x(1)*u23%x(1)-u31%x(2)*u23%x(2)-u31%x(3)*u23%x(3)

    IF (z_ca1 < -1.0_wp) z_ca1 = -1.0_wp
    IF (z_ca1 >  1.0_wp) z_ca1 =  1.0_wp
    IF (z_ca2 < -1.0_wp) z_ca2 = -1.0_wp
    IF (z_ca2 >  1.0_wp) z_ca2 =  1.0_wp
    IF (z_ca3 < -1.0_wp) z_ca3 = -1.0_wp
    IF (z_ca3 >  1.0_wp) z_ca3 =  1.0_wp

    z_a1 = ACOS(z_ca1)
    z_a2 = ACOS(z_ca2)
    z_a3 = ACOS(z_ca3)

    !  Compute areas = z_a1 + z_a2 + z_a3 - pi.
    area = z_a1+z_a2+z_a3-pi

    IF ( area < 0.0_wp ) area = 0.0_wp

  END FUNCTION triangle_area
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Computes the vector product of x0,x1.
  !!
  !! @par Revision History
  !! Developed  by Luis Kornblueh  (2004).
  ELEMENTAL FUNCTION vector_product (x0, x1) result(x2)

    TYPE(t_cartesian_coordinates), INTENT(in) :: x0, x1
    TYPE(t_cartesian_coordinates) :: x2

    x2%x(1) = x0%x(2)*x1%x(3) - x0%x(3)*x1%x(2)
    x2%x(2) = x0%x(3)*x1%x(1) - x0%x(1)*x1%x(3)
    x2%x(3) = x0%x(1)*x1%x(2) - x0%x(2)*x1%x(1)

  END FUNCTION vector_product
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Computes the norma of the vector product of x0,x1.
  !!
  !! @par Revision History
  !! Developed  by Leonidas Linardakis (2020).
  ELEMENTAL FUNCTION norma_of_vector_product (x0, x1) result(r)

    TYPE(t_cartesian_coordinates), INTENT(in) :: x0, x1
    REAL(wp) :: r

    REAL(wp) :: y1, y2, y3

    y1 = x0%x(2)*x1%x(3) - x0%x(3)*x1%x(2)
    y2 = x0%x(3)*x1%x(1) - x0%x(1)*x1%x(3)
    y3 = x0%x(1)*x1%x(2) - x0%x(2)*x1%x(1)

    r = SQRT(y1*y1 + y2*y2 + y3*y3)

  END FUNCTION norma_of_vector_product
  !-------------------------------------------------------------------------



  !-------------------------------------------------------------------------
  !>
  !! Computes the vector product of of 2 unit vectors.
  ELEMENTAL FUNCTION vecpro(geog1,geog2) result (xx)

    TYPE(t_geographical_coordinates), INTENT(in) :: geog2, geog1

    TYPE(t_cartesian_coordinates)  :: xx
    REAL (wp)  :: long1, lat1,long2,lat2
    REAL (wp) :: sn1,sn2,sn3,sn4,cs1,cs2

    long1=geog1%lon
    lat1=geog1%lat
    long2=geog2%lon
    lat2=geog2%lat

    sn1=SIN(lat1-lat2)
    sn2=SIN(lat1+lat2)
    sn3=SIN(0.5_wp*(long1-long2))
    sn4=SIN(0.5_wp*(long1+long2))
    cs1=COS(0.5_wp*(long1-long2))
    cs2=COS(0.5_wp*(long1+long2))

    xx%x(1)=sn1*sn4*cs1-sn2*cs2*sn3
    xx%x(2)=sn1*cs2*cs1+sn2*sn4*sn3
    xx%x(3)=COS(lat1)*cos(lat2)*sin(long1-long2)

  END FUNCTION vecpro
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Rotate counterclockwise (seen from positive axis direction)
  !! vector x1 around axis z for phirot radians
  !!
  !! @par Revision History
  !! Developed  by Luca Bonaventura  (2005).
  ELEMENTAL FUNCTION rotate_z (p_phirot,x1) result(x2)

    REAL(wp), INTENT(in) :: p_phirot
    TYPE(t_cartesian_coordinates), INTENT(in) :: x1
    TYPE(t_cartesian_coordinates) :: x2

    x2%x(1) = COS(p_phirot)*x1%x(1) - SIN(p_phirot)*x1%x(2)
    x2%x(2) = SIN(p_phirot)*x1%x(1) + COS(p_phirot)*x1%x(2)
    x2%x(3) = x1%x(3)

  END FUNCTION rotate_z
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Rotate counterclockwise (seen from positive axis direction)
  !! vector x1 around axis y for phirot radians
  !!
  !! @par Revision History
  !! Developed  by Luca Bonaventura  (2005).
  ELEMENTAL FUNCTION rotate_y (p_phirot,x1) result(x2)

    !-----------------------------------------------------------------------
    REAL(wp), INTENT(in) :: p_phirot
    TYPE(t_cartesian_coordinates), INTENT(in) :: x1
    TYPE(t_cartesian_coordinates) :: x2

    x2%x(1) =  COS(p_phirot)*x1%x(1) + SIN(p_phirot)*x1%x(3)
    x2%x(3) = -SIN(p_phirot)*x1%x(1) + COS(p_phirot)*x1%x(3)
    x2%x(2) =  x1%x(2)

  END FUNCTION rotate_y
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Rotate counterclockwise (seen from positive axis direction)
  !! vector x1 around axis x for phirot radians
  !!
  !! @par Revision History
  !! Developed  by Luca Bonaventura  (2005).
  ELEMENTAL FUNCTION rotate_x (p_phirot,x1) result(x2)

    !-----------------------------------------------------------------------
    REAL(wp), INTENT(in) :: p_phirot
    TYPE(t_cartesian_coordinates), INTENT(in) :: x1
    TYPE(t_cartesian_coordinates) :: x2

    x2%x(2) =  COS(p_phirot)*x1%x(2) + SIN(p_phirot)*x1%x(3)
    x2%x(3) = -SIN(p_phirot)*x1%x(2) + COS(p_phirot)*x1%x(3)
    x2%x(1) =  x1%x(1)

  END FUNCTION rotate_x
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  !! Computes normal vector to plane determined by x0,x1.
  !!
  !! @par Revision History
  !! Developed  by Luis Kornblueh  (2004).
  ELEMENTAL FUNCTION normal_vector (x0, x1) result(x2)

    TYPE(t_cartesian_coordinates), INTENT(in) :: x0, x1
    TYPE(t_cartesian_coordinates) :: x2

    x2%x(1) = x0%x(2)*x1%x(3) - x0%x(3)*x1%x(2)
    x2%x(2) = x0%x(3)*x1%x(1) - x0%x(1)*x1%x(3)
    x2%x(3) = x0%x(1)*x1%x(2) - x0%x(2)*x1%x(1)

!    x2%x = x2%x/ SQRT(x2%x(1)**2+x2%x(2)**2+x2%x(3)**2)
    x2%x = x2%x/ d_norma_3d(x2)

  END FUNCTION normal_vector
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  ! Projects a point to plane defined by three points:
  ! to_v0, to_v1, to_v2
  ELEMENTAL FUNCTION project_point (in_point, to_v0, to_v1, to_v2) result(projected_point)
    TYPE(t_cartesian_coordinates), INTENT(in) :: in_point, to_v0, to_v1, to_v2
    TYPE(t_cartesian_coordinates) :: projected_point

    TYPE(t_cartesian_coordinates) :: in_p, v1, v2, norm_to_plane

    ! translate everybody to the center
    v1%x = to_v1%x - to_v0%x
    v2%x = to_v2%x - to_v0%x
    in_p%x = in_point%x - to_v0%x

    norm_to_plane = vector_product (v1, v2)
    d_normalize(norm_to_plane)
    norm_to_plane%x = norm_to_plane%x * &
      & DOT_PRODUCT(norm_to_plane%x, (v1%x - in_p%x))

    projected_point%x = in_p%x + norm_to_plane%x + to_v0%x

  END FUNCTION project_point
  !-------------------------------------------------------------------------
    
    
  

  !-------------------------------------------------------------------------
  !>
  !! Determines the bary_center of triangle with vertices v0,v1,v2.
  !!
  !! @par Revision History
  !! original version by Thomas Heinze (2006-08-17).
  !!
  !! @par Remarks
  !! currently not used in the code, just for testing purposes
  ELEMENTAL FUNCTION bary_center (v0, v1, v2) result(center)
    !> the coordinates of the three triangle vertices (unit vectors)
    !! in counter clockwise order.
    TYPE(t_cartesian_coordinates), INTENT(in) :: v0, v1, v2
    !> the coordinates of the bary_center, unless co-linear
    TYPE(t_cartesian_coordinates) :: center

    !  Local variables:
    TYPE(t_cartesian_coordinates) :: bc ! barycenter of the underlying planar
    ! triangle:
    REAL(wp) :: z_cnorm               ! norm of e1

    !-----------------------------------------------------------------------

    ! get barycenter of planar triangle

    bc%x = v0%x + v1%x + v2%x
    bc%x = bc%x/3._wp

    ! norm of barycenter
    z_cnorm = d_norma_3d(bc)

    ! new center
    center%x = bc%x/z_cnorm

  END FUNCTION bary_center
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Determines the circum_center of triangle with vertices v0,v1,v2.
  !!
  !! @par Revision History
  !! Developed  by Luis Kornblueh  (2004).
  ELEMENTAL FUNCTION circum_center (v0, v1, v2) result(center)
    !> the coordinates of the three triangle vertices (unit vectors)
    !! in counter clockwise order.
    TYPE(t_cartesian_coordinates), INTENT(in) :: v0, v1, v2
    !> the coordinates of the circumcenter unless co-linear
    TYPE(t_cartesian_coordinates) :: center

    !  Local variables:
    TYPE(t_cartesian_coordinates) :: e1, e2 ! edges of the underlying planar triangle:
    ! v1-v0 ands v2-v0, respectively
    TYPE(t_cartesian_coordinates) :: cu     ! dot product of center:  e1 x e2
    REAL(wp) :: z_cnorm                     ! norm of cu

    !-----------------------------------------------------------------------

    e1%x = v1%x - v0%x
    e2%x = v2%x - v0%x

    ! compute cu = e1 x e2 and cnorm**2.
    cu = vector_product (e1, e2)

    IF (DOT_PRODUCT(cu%x, v0%x) < 0.0_wp) THEN
      cu%x = -cu%x
    END IF

    z_cnorm = d_norma_3d(cu)
    center%x = cu%x/z_cnorm

  END FUNCTION circum_center
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Computes length of geodesic arc with endpoints x0,x1.
  !!
  !! @par Revision History
  !! Developed by Th.Heinze (2006-09-19).
  !! Previous version by Luis Kornblueh (2004) discarded.
  ELEMENTAL FUNCTION arc_length (x0, x1) result(arc)

    TYPE(t_cartesian_coordinates), INTENT(in) :: x0, x1
    REAL(wp) :: arc
    REAL(wp) :: z_d0,  z_d1,  z_cc

    z_d0 = d_norma_3d(x0)
    z_d1 = d_norma_3d(x1)

    z_cc = DOT_PRODUCT(x0%x, x1%x)/(z_d0*z_d1)

    ! in case we get numerically incorrect solutions
    IF (z_cc > 1._wp )  z_cc =  1._wp
    IF (z_cc < -1._wp ) z_cc = -1._wp

    arc = ACOS(z_cc)

  END FUNCTION arc_length
  !-----------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Computes length of geodesic arc with endpoints x0,x1.
  !!
  !! @par Revision History
  !! Developed by Th.Heinze (2006-09-19).
  !! Previous version by Luis Kornblueh (2004) discarded.
  ELEMENTAL FUNCTION arc_length_normalsphere (x0, x1) result(arc)

    TYPE(t_cartesian_coordinates), INTENT(in) :: x0, x1
    REAL(wp) :: arc
    REAL(wp) :: z_cc


    z_cc = DOT_PRODUCT(x0%x, x1%x)

    ! in case we get numerically incorrect solutions
    IF (z_cc > 1._wp )  z_cc =  1._wp
    IF (z_cc < -1._wp ) z_cc = -1._wp

    arc = ACOS(z_cc)

  END FUNCTION arc_length_normalsphere
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  !! Computes the cosine of the length of geodesic arc with endpoints x0,x1.
  !!
  !! @par Revision History
  !! Developed by Almut Gassmann (2007-03-13).
  ELEMENTAL FUNCTION cos_arc_length (x0, x1) result(carc)
   TYPE(t_cartesian_coordinates), INTENT(in) :: x0, x1
    REAL(wp) :: carc
    REAL(wp) :: z_d0,  z_d1,  z_cc

    z_d0 = d_norma_3d(x0)
    z_d1 = d_norma_3d(x1)

    z_cc = DOT_PRODUCT(x0%x, x1%x)/(z_d0*z_d1)

    ! in case we get numerically incorrect solutions
    IF (z_cc > 1._wp )  z_cc =  1._wp
    IF (z_cc < -1._wp ) z_cc = -1._wp

    carc = z_cc

  END FUNCTION cos_arc_length
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  !! Computes the sin between two vectors x0, x1
  !!
  !! @par Revision History
  !! Developed by Leonidas Linardakis (2010-05-13).
  ELEMENTAL FUNCTION sin_cc (x0, x1)
    TYPE(t_cartesian_coordinates), INTENT(in) :: x0, x1
    REAL(wp) :: sin_cc

    TYPE(t_cartesian_coordinates) :: r

    r = vector_product(x0, x1)
    sin_cc = d_norma_3d(r) / &
      & (d_norma_3d(x0) * d_norma_3d(x1))

  END FUNCTION sin_cc
  !-----------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Finds the intersection of two great circles.
  !!
  !! Cannot be made elemental, as long as the
  !! subroutine calls to message and finish are maintained.
  !!
  !! @par Revision History
  !! Developed  by Luis Kornblueh  (2004).
  !!
  ELEMENTAL FUNCTION inter_section (p0, p1, v0, v1) result(p)
    TYPE(t_cartesian_coordinates), INTENT(in) :: p0, p1
    TYPE(t_cartesian_coordinates), INTENT(in) :: v0, v1
    TYPE(t_cartesian_coordinates) :: p

    TYPE(t_cartesian_coordinates) :: en, cn

    REAL(wp) :: z_pn

    ! normal to plane where edge lies

    en = vector_product (v0, v1)

    ! normal tp plane where triangle centers lies

    cn = vector_product (p0, p1)

    p = vector_product (en, cn)

    ! rescale on sphere surface

    z_pn = d_norma_3d(p)

    p%x = p%x/z_pn

    IF (DOT_PRODUCT(p0%x, p%x) < 0.0_wp) THEN
      p%x = -p%x
    END IF

  END FUNCTION inter_section
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Finds the intersection of two great circles.
  !!
  !! Cannot be made elemental, as long as the
  !! subroutine calls to message and finish are maintained.
  !!
  !! @par Revision History
  !! Developed  by Luis Kornblueh  (2004).
  !!
#ifdef __SX__
  SUBROUTINE inter_section2 (p0, p1, v0, v1, p, nlen)
    TYPE(t_cartesian_coordinates), INTENT(in) :: p0(0:nlen-1), p1(0:nlen-1)
    TYPE(t_cartesian_coordinates), INTENT(in) :: v0(0:nlen-1), v1(0:nlen-1)
    TYPE(t_cartesian_coordinates), INTENT(out):: p(0:nlen-1)
    INTEGER, INTENT(in)                       :: nlen

    TYPE(t_geographical_coordinates) :: g1,g2
    TYPE(t_cartesian_coordinates) :: en, cn

    REAL(wp) :: z_pn
    INTEGER  :: ji

    DO ji = 0, nlen-1
      g1=cc2gc(v0(ji))
      g2=cc2gc(v1(ji))

      en = vecpro (g1, g2)

      g1=cc2gc(p0(ji))
      g2=cc2gc(p1(ji))

      cn = vecpro (g1, g2)

      g1=cc2gc(en)
      g2=cc2gc(cn)

      !  p = vector_product (en, cn)
      p(ji) = vecpro (g1, g2)

      ! rescale on sphere surface
      z_pn = d_norma_3d(p(ji))

      p(ji)%x = p(ji)%x/z_pn

      IF (DOT_PRODUCT(p0(ji)%x, p(ji)%x) < 0.0_wp) THEN
        p(ji)%x = -p(ji)%x
      END IF
    ENDDO

  END SUBROUTINE inter_section2
#else
  ELEMENTAL FUNCTION inter_section2 (p0, p1, v0, v1) result(p)
    TYPE(t_cartesian_coordinates), INTENT(in) :: p0, p1
    TYPE(t_cartesian_coordinates), INTENT(in) :: v0, v1
    TYPE(t_cartesian_coordinates) :: p

    TYPE(t_geographical_coordinates) :: g1,g2
    TYPE(t_cartesian_coordinates) :: en, cn

    !!$    REAL(wp), PARAMETER :: voronoi_tolerance = 1.0e-9_wp

    REAL(wp) :: z_pn

    ! normal to plane where edge lies
    !  en = vector_product (v0, v1)

    g1=cc2gc(v0)
    g2=cc2gc(v1)

    en = vecpro (g1, g2)

    ! normal tp plane where triangle centers lies
    ! cn = vector_product (p0, p1)


    g1=cc2gc(p0)
    g2=cc2gc(p1)

    cn = vecpro (g1, g2)

    ! check Voronoi/Delaunay properties
    !!$    IF (ABS(DOT_PRODUCT(en%x, cn%x)) > voronoi_tolerance) THEN
    !!$      WRITE (message_text,'(a,e12.6)') &
    !!$           'Orthogonality = ',                &
    !!$           ABS(DOT_PRODUCT(en%x, cn%x))
    !!$      CALL message ('inter_sect', TRIM(message_text))
    !!$      CALL finish ('inter_sect','Voronoi orthogonality violated.')
    !!$    END IF

    g1=cc2gc(en)
    g2=cc2gc(cn)

    !  p = vector_product (en, cn)
    p = vecpro (g1, g2)

    ! rescale on sphere surface
    z_pn = d_norma_3d(p)

    p%x = p%x/z_pn

    IF (DOT_PRODUCT(p0%x, p%x) < 0.0_wp) THEN
      p%x = -p%x
    END IF

  END FUNCTION inter_section2
#endif
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Converts zonal @f$p\_gu@f$ and meridional vector components @f$p\_gv@f$ into cartesian.
  !!
  !! Converts zonal @f$p\_gu@f$ and meridional vector components @f$p\_gv@f$ into cartesian
  !! ones @f$(p\_cu, p\_cv, p\_cw)@f$ using the conversion
  !! @f{align*}{
  !! \begin{pmatrix} p\_cu \\ p\_cv \\ p\_cw \end{pmatrix} =
  !! \begin{pmatrix}
  !! -\sin p\_long & -\sin p\_lat \cdot \cos p\_long \\\
  !! \cos p\_long & -\sin p\_lat \cdot \sin p\_long \\\
  !! 0 & -\cos p\_lat
  !! \end{pmatrix} \cdot
  !! \begin{pmatrix} p\_gu \\ p\_gv \end{pmatrix}
  !! @f}
  !!
  !! @par Revision History
  !! Original version by Tobias Ruppert and Thomas Heinze, DWD (2006-11-14)
  PURE SUBROUTINE gvec2cvec (p_gu, p_gv, p_long, p_lat, p_cu, p_cv, p_cw)
    !
    REAL(wp), INTENT(in)  :: p_gu, p_gv     ! zonal and meridional vec. component
    REAL(wp), INTENT(in)  :: p_long, p_lat  ! geo. coord. of data point

    REAL(wp), INTENT(out) :: p_cu, p_cv, p_cw            ! Cart. vector

    REAL(wp)              :: z_cln, z_sln, z_clt, z_slt  ! sin and cos of
    ! p_long and p_lat

    !-------------------------------------------------------------------------

    z_sln = SIN(p_long)
    z_cln = COS(p_long)
    z_slt = SIN(p_lat)
    z_clt = COS(p_lat)

    p_cu = z_sln * p_gu + z_slt * z_cln * p_gv
    p_cu = -1._wp * p_cu
    p_cv = z_cln * p_gu - z_slt * z_sln * p_gv
    p_cw = z_clt * p_gv

  END SUBROUTINE gvec2cvec
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Returns a coordinate system (x,y),
  !! for the tanget plane at the spehere on the given point (lon,lat)
  PURE SUBROUTINE sphere_tanget_coordinates (geo_point,x,y)
    TYPE(t_geographical_coordinates), INTENT(in)  :: geo_point
    TYPE(t_cartesian_coordinates), INTENT(out) :: x,y

    REAL(wp) :: lon,lat
    !-------------------------------------------------------------------------
     lon = geo_point%lon
     lat = geo_point%lat
     x%x(1) = -SIN(lon)
     x%x(2) =  COS(lon)
     x%x(3) =  0._wp
     y%x(1) = -COS(lon) * SIN(lat)
     y%x(2) = -SIN(lon) * SIN(lat)
     y%x(3) =  COS(lat)
    !-------------------------------------------------------------------------

  END SUBROUTINE sphere_tanget_coordinates
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  ELEMENTAL FUNCTION norma (v) result(length)
    TYPE(t_cartesian_coordinates), INTENT(in) :: v
    REAL(wp) :: length

    ! length = SQRT(v%x(1)**2+v%x(2)**2+v%x(3)**2)
    length = d_norma_3d(v)

  END FUNCTION norma
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  PURE FUNCTION normalize(v) result(n)
    TYPE(t_cartesian_coordinates), INTENT(in) :: v
    TYPE(t_cartesian_coordinates) :: n

 !   v%x = v%x / (SQRT(v%x(1)**2+v%x(2)**2+v%x(3)**2))
    n%x = v%x / d_norma_3d(v)

  END FUNCTION normalize
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  ELEMENTAL FUNCTION sphere_cartesian_midpoint (p1, p2) result(mid_point)
!  FUNCTION sphere_cartesian_midpoint (p1, p2) result(mid_point)
    TYPE(t_cartesian_coordinates), INTENT(in) :: p1, p2
    TYPE(t_cartesian_coordinates) :: mid_point

!    TYPE(cartesian_coordinates) :: w
!    REAL(wp) :: zchord,ztheta,zgamma,zbeta,zalpha

      mid_point%x = (p1%x + p2%x)
      mid_point%x = mid_point%x / &
        & d_norma_3d(mid_point)

!     w%x = p1%x - p2%x
!     zchord = norma(w)
!     ztheta = 2._wp*asin (0.5_wp*zchord)
!     zgamma = 0.5_wp
!     zbeta  = SIN (zgamma*ztheta)/SIN (ztheta)
!     zalpha = SIN ((1.0_wp-zgamma)*ztheta)/SIN (ztheta)
!     mid_point%x = zalpha * p1%x + zbeta * p2%x

!     print *, 'zchord,ztheta,zgamma,zbeta,zalpha',&
!       & zchord,ztheta,zgamma,zbeta,zalpha
!     print *, 'p1:',p1
!     print *, 'p2:',p2
!     print *, 'mid_point:',mid_point

  END FUNCTION sphere_cartesian_midpoint
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  ELEMENTAL FUNCTION angle_of_vectors(v1,v2) result(angle)
    TYPE(t_cartesian_coordinates), INTENT(in) :: v1,v2
    REAL(wp) :: angle

    REAL(wp) :: prod
   ! TYPE(cartesian_coordinates) :: v3

    prod = DOT_PRODUCT(v1%x,v2%x)/ &
      & (d_norma_3d(v1) * d_norma_3d(v2))
    angle = ATAN2(SQRT(1.0_wp - prod*prod), prod)

!    angle = ACOS(prod)

  END FUNCTION angle_of_vectors
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  ELEMENTAL FUNCTION angle_of_points(p1, p2, p3) result(angle)
    TYPE(t_cartesian_coordinates), INTENT(in) :: p1,p2,p3

    REAL(wp) :: prod,angle
    TYPE(t_cartesian_coordinates) :: v1,v2

    v1%x = p1%x - p2%x
    v2%x = p3%x - p2%x
    prod = DOT_PRODUCT(v1%x,v2%x)/ (norma(v1) * norma(v2))
    angle = ACOS(prod)

  END FUNCTION angle_of_points
  !-------------------------------------------------------------------------

END MODULE mo_base_geometry
!----------------------------------------------------------------------------









