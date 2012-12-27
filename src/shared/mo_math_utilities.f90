
! LL: xlc has trouble optimizing routines with implicit shaped parameters
! #ifdef __xlC__
! @process nohot
! ! @PROCESS NOOPTIMIZE
! #endif
!>
!!   Contains the implementation of various mathematical algorithms.
!!
!!   Contains the implementation of various mathematical algorithms
!!   used by the shallow water model.
!!
!! @par Revision History
!!  Developed  by Luca Bonaventura (2002-4).
!!  Modified to ProTeX-style by  Luca Bonaventura and Thomas Heinze (2004).
!!  Adapted to new data structure by Thomas Heinze,
!!  Peter Korn and Luca Bonaventura (2005).
!!  Modification by Thomas Heinze (2006-02-21):
!!  - renamed m_modules to mo_modules
!!  Modification by Peter Korn and Luca Bonaventura (2006-07-21):
!!  - moved here linear algebra subroutines and other auxiliary functions
!!  - added some Protex documentation and great cleanup
!!  Modification by Thomas Heinze (2006-09-07):
!!  - added functions barycenter and bary_center
!!  Modification by Thomas Heinze (2006-10-19):
!!  - added functions vector_product and integral_over_triangle
!!  Modification by Tobias Ruppert and Thomas Heinze (2006-11-14):
!!  - added functions solve_chol and choldec
!!  - renamed function solve to solve_lu
!!  Modification by Thomas Heinze (2006-11-16):
!!  - added functions cvec2gvec and gvec2cvec
!!  Modification by Peter Korn, MPI-M, (2006-11-23):
!!  - replaced vertex_index by vertex_idx
!!  Modification by Hui Wan (2007-02-22):
!!  - functions barycenter and bary_center removed.
!!    (These two functions used TYPE grid, thus should not be put
!!     in to shr_general.)
!!  - function ll2xyz removed because it had been replaced by
!!    gc2cc and was not used anymore.
!!  - type cartesian_coordinates and type geographical_coordinates
!!    moved to this module from mo_model_domain
!!  - functions func_f and dxg moved from mo_functions to here
!!  Problem related to the range of longitude in the spherical coordnate
!!  identified by Th. Heinze and P. Ripodas (2007-02-26). Correction in 'cc2gc'
!!  made by H. Wan (2007-02-27).
!!  Modification by Thomas Heinze, DWD, (2007-07-26):
!!  - including all the improvements of Tobias Ruppert's diploma thesis
!!  - several changes according to the programming guide
!!  - functions func_f and dxg moved from here to mo_interpolation
!!  Modification by Daniel Reinert, DWD, (2009-10-30)
!!  - added subroutine gnomonic_proj which uses a gnomonic projection in
!!    order to project a point (lat_1,lon_1) onto a tangent plane with
!!    the origin at (lat_0,lon_0). The results are the local cartesian
!!    coordinates (x_1,y_1)
!!  - added subroutine orthogr_proj. performs orthographic projection of
!!    point (lat_1,lon_1) onto a tangent plane with the origin at (lat_0,lon_0).
!!  Modification by Daniel Reinert, DWD, (2012-04-04)
!!  - added function which can be used to check whether to line segments
!!    intersect (in 2D cartesian system)
!!  Modification by Daniel Reinert, DWD, (2012-04-05)
!!  - added function which computes the intersection point between 2 lines
!!    (2D cartesian)
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
!!
MODULE mo_math_utilities
  !-------------------------------------------------------------------------
  !
  !    ProTeX FORTRAN source: Style 2
  !    modified for ICON project, DWD/MPI-M 2006
  !
  !-------------------------------------------------------------------------
  !
  USE mo_kind,                ONLY: wp
  USE mo_math_constants,      ONLY: pi, pi_2, dbl_eps
  USE mo_exception,           ONLY: message, finish
  USE mo_parallel_config,     ONLY: nproma
  USE mo_lonlat_grid,         ONLY: t_lon_lat_grid
  USE mo_grid_config,         ONLY: grid_sphere_radius
  USE mo_grid_geometry_info,  ONLY: t_grid_geometry_info
  USE mo_math_types
  IMPLICIT NONE
  
  PRIVATE
  
  CHARACTER(LEN=*), PARAMETER :: version = '$Id$'
  
  PUBLIC :: t_cartesian_coordinates
  PUBLIC :: t_geographical_coordinates
  PUBLIC :: t_line
  PUBLIC :: t_tangent_vectors
  
  PUBLIC :: cc2gc
  PUBLIC :: gc2cc
  PUBLIC :: cc2tv
  PUBLIC :: gvec2cvec
  PUBLIC :: cvec2gvec
  PUBLIC :: arc_length
  PUBLIC :: arc_length_v
  PUBLIC :: cartesian_to_geographical, geographical_to_cartesian
  PUBLIC :: arc_length_on_unitsphere
  PUBLIC :: cos_arc_length
  PUBLIC :: plane_torus_distance, plane_torus_closest_coordinates
  PUBLIC :: plane_angle_of_three_points
  PUBLIC :: sphere_cartesian_midpoint
  PUBLIC :: norma, normalize
  PUBLIC :: angle_of_vectors, sin_cc
  PUBLIC :: triangle_area
  PUBLIC :: circum_center
  PUBLIC :: bary_center
  PUBLIC :: middle
  PUBLIC :: sphere_tangent_coordinates
  PUBLIC :: normal_vector
  PUBLIC :: spherical_intersection, spherical_intersection2
  PUBLIC :: rotate_z,rotate_x,rotate_y
  PUBLIC :: norma_of_vector_product
  PUBLIC :: project_point_to_plane
  
  PUBLIC :: disp_new
  PUBLIC :: disp_new_vect
  
  PUBLIC :: cc_dot_product
  PUBLIC :: cc_norm
  PUBLIC :: vector_product
  PUBLIC :: integral_over_triangle
  PUBLIC :: rotate_latlon
  PUBLIC :: rotate_latlon_vec
  PUBLIC :: rotate_latlon_grid
  PUBLIC :: latlon_compute_area_weights
  
  PUBLIC :: gnomonic_proj
  PUBLIC :: orthogr_proj
  PUBLIC :: az_eqdist_proj
  PUBLIC :: gamma_fct
!   PUBLIC :: sphere_cell_mean_char_length
  PUBLIC :: ccw
  PUBLIC :: line_intersect
  PUBLIC :: lintersect
  
  PUBLIC :: OPERATOR(+)
  PUBLIC :: OPERATOR(-)
  PUBLIC :: OPERATOR(*)    
  
  INTERFACE OPERATOR(+)
    MODULE PROCEDURE cartesian_coordinates_plus
  END INTERFACE
  INTERFACE OPERATOR(-)
    MODULE PROCEDURE cartesian_coordinates_minus
  END INTERFACE
  INTERFACE OPERATOR(*)
    MODULE PROCEDURE cartesian_coordinates_prod
  END INTERFACE

!-----------------------------------------------------------------------
! Basic geometry functions definitions
#define d_norma_3d(v) SQRT(DOT_PRODUCT(v%x,v%x))
#define d_norma_3d_x(x) SQRT(DOT_PRODUCT(x,x))
#define d_normalize(v) v%x=v%x/d_norma_3d(v)
#define d_normalize_x(x) x=x/d_norma_3d_x(x)
#define d_arc_of_hord_normalsphere(l) (2.0_wp*ASIN((l)*0.5_wp))
#define d_normalize_f(v) v%x/d_norma_3d(v)
#define d_sqrdistance_3d(v1,v2) DOT_PRODUCT((v1%x-v2%x),(v1%x-v2%x))
  
CONTAINS
  
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
  !!
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
  !! Converts cartesian velocity vector @f$(p\_cu, p\_cv, p\_cw)@f$.
  !!
  !! Converts cartesian velocity vector @f$(p\_cu, p\_cv, p\_cw)@f$
  !! into zonal @f$p\_gu@f$ and meridional @f$g\_v@f$ velocity components
  !! using the conversion
  !! @f{align*}{
  !! \begin{pmatrix} p\_gu \\ p\_gv \end{pmatrix} =
  !! \begin{pmatrix}
  !! -\sin p\_long & \cos p\_long & 0 \\\
  !! -\sin p\_lat \cdot \cos p\_long & -\sin p\_lat \cdot \sin p\_long & \cos p\_lat
  !! \end{pmatrix} \cdot
  !! \begin{pmatrix} p\_cu \\ p\_cv \\ p\_cw \end{pmatrix}
  !! @f}
  !!
  !! @par Revision History
  !! Original version by Thomas Heinze, DWD (2006-11-16)
  !!
  PURE SUBROUTINE cvec2gvec (p_cu, p_cv, p_cw, p_long, p_lat, p_gu, p_gv)
    !
    REAL(wp), INTENT(in)  :: p_cu, p_cv, p_cw  ! Cart. vector
    REAL(wp), INTENT(in)  :: p_long, p_lat     ! geo. coord. of data point
    
    REAL(wp), INTENT(out) :: p_gu, p_gv        ! zonal and meridional vec. comp.
    
    REAL(wp)              :: z_cln, z_clt, z_sln, z_slt  ! sin and cos of
    ! p_long and p_lat
    
    !-------------------------------------------------------------------------    
    z_sln = SIN(p_long)
    z_cln = COS(p_long)
    z_slt = SIN(p_lat)
    z_clt = COS(p_lat)
    
    p_gu = z_cln * p_cv - z_sln * p_cu
    p_gv = z_cln * p_cu + z_sln * p_cv
    p_gv = z_slt * p_gv
    p_gv = z_clt * p_cw - p_gv
    
  END SUBROUTINE cvec2gvec
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !>
  !! Converts  vectors (in)  cartesian coordinate representation
  !! to tangent vectors.
  !!
  !! @par Revision History
  !! Developed  by Luca Bonaventura  (2005).
  FUNCTION cc2tv(xx,position) result (tt)

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
  ELEMENTAL FUNCTION norma (v) result(length)
    TYPE(t_cartesian_coordinates), INTENT(in) :: v
    REAL(wp) :: length

    length = d_norma_3d(v)

  END FUNCTION norma
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  PURE FUNCTION normalize(v) result(n)
    TYPE(t_cartesian_coordinates), INTENT(in) :: v
    TYPE(t_cartesian_coordinates) :: n

    n%x = v%x / d_norma_3d(v)

  END FUNCTION normalize
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

  END FUNCTION angle_of_vectors
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !>
  !! Returns a coordinate system (x,y),
  !! for the tanget plane at the spehere on the given point (lon,lat)
  PURE SUBROUTINE sphere_tangent_coordinates (geo_point,x,y)
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

  END SUBROUTINE sphere_tangent_coordinates
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  !! Computes normal vector to plane determined by the vectors x0,x1.
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

  !--------------------------------------------------------------------
  !>
  ELEMENTAL FUNCTION middle(p1,p2)

    TYPE(t_cartesian_coordinates), INTENT(in) :: p1,p2
    TYPE(t_cartesian_coordinates) :: middle

    middle%x = (p1%x + p2%x) / 2.0_wp
  
  END FUNCTION middle
  !--------------------------------------------------------------------
  
  
  !-------------------------------------------------------------------------
  !>
  !! Computes length of geodesic arc with endpoints @f$p\_x, p\_y@f$.
  !!
  !! @par Revision History
  !! Developed by Th.Heinze (2006-09-19).
  !! Previous version by Luis Kornblueh (2004) discarded.
  !!
  ELEMENTAL FUNCTION arc_length (p_x, p_y)  result (p_arc)
  
    TYPE(t_cartesian_coordinates), INTENT(in) :: p_x, p_y  ! endpoints
    
    REAL(wp)            :: p_arc          ! length of geodesic arc
    
    REAL(wp)            :: z_lx,  z_ly    ! length of vector p_x and p_y
    REAL(wp)            :: z_cc           ! cos of angle between endpoints
    
    !-----------------------------------------------------------------------
    
    z_lx = d_norma_3d(p_x)
    z_ly = d_norma_3d(p_y)
    
    z_cc = DOT_PRODUCT(p_x%x, p_y%x)/(z_lx*z_ly)
    
    ! in case we get numerically incorrect solutions    
    IF (z_cc > 1._wp )  z_cc =  1._wp
    IF (z_cc < -1._wp ) z_cc = -1._wp
    
    p_arc = ACOS(z_cc)
    
  END FUNCTION arc_length
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !>
  !! Computes length of geodesic arc with endpoints x0,x1.
  !!
  !! @par Revision History
  !! Developed by Th.Heinze (2006-09-19).
  !! Previous version by Luis Kornblueh (2004) discarded.
  ELEMENTAL FUNCTION arc_length_on_unitsphere (x0, x1) result(arc)

    TYPE(t_cartesian_coordinates), INTENT(in) :: x0, x1
    REAL(wp) :: arc
    REAL(wp) :: z_cc


    z_cc = DOT_PRODUCT(x0%x, x1%x)

    ! in case we get numerically incorrect solutions
    IF (z_cc > 1._wp )  z_cc =  1._wp
    IF (z_cc < -1._wp ) z_cc = -1._wp

    arc = ACOS(z_cc)

  END FUNCTION arc_length_on_unitsphere
  !-----------------------------------------------------------------------
    
  !-------------------------------------------------------------------------
  !>
  !! Computes length of geodesic arc with endpoints @f$p\_x, p\_y@f$.
  !!
  !! Vectorizable version
  !!
  !! @par Revision History
  !! Developed by Th.Heinze (2006-09-19).
  !! Previous version by Luis Kornblueh (2004) discarded.
  !! Vectorizable version developed by Guenther Zaengl, DWD (2009-04-20)
  !!
  PURE FUNCTION arc_length_v (p_x, p_y)  result (p_arc)
    REAL(wp), INTENT(in) :: p_x(3), p_y(3)  ! endpoints
    
    REAL(wp)            :: p_arc          ! length of geodesic arc
    
    REAL(wp)            :: z_lx,  z_ly    ! length of vector p_x and p_y
    REAL(wp)            :: z_cc           ! cos of angle between endpoints
    
    !-----------------------------------------------------------------------
    
    z_lx = SQRT(DOT_PRODUCT(p_x,p_x))
    z_ly = SQRT(DOT_PRODUCT(p_y,p_y))
    
    z_cc = DOT_PRODUCT(p_x, p_y)/(z_lx*z_ly)
    
    ! in case we get numerically incorrect solutions
    
    IF (z_cc > 1._wp )  z_cc =  1._wp
    IF (z_cc < -1._wp ) z_cc = -1._wp
    
    p_arc = ACOS(z_cc)
    
  END FUNCTION arc_length_v
  !-------------------------------------------------------------------------
  
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
  ! returns the torus modulo coordinates of v1 that are
  ! closest to v0
  FUNCTION plane_torus_closest_coordinates(v0, v1, &
    & geometry_info) result(new_v1_coord)
    
    REAL(wp), INTENT(in) :: v0(3), v1(3)
    TYPE(t_grid_geometry_info), INTENT(in) :: geometry_info
    
    REAL(wp) :: length_of_torus, height_of_torus
    TYPE(t_cartesian_coordinates) :: new_v1_coord

    length_of_torus = geometry_info%domain_length
    height_of_torus = geometry_info%domain_height
    ! check the x coordinate
    IF ( ABS(v0(1) - v1(1)) >  length_of_torus * 0.5_wp) THEN
      ! we will wrap around + or -  length_of_torus       
      IF (v0(1) > v1(1)) THEN
        new_v1_coord%x(1) = v1(1) + length_of_torus
      ELSE
        new_v1_coord%x(1) = v1(1) - length_of_torus
      ENDIF
    ELSE
      ! keep the same 
      new_v1_coord%x(1) = v1(1)
    ENDIF

    ! check the y coordinate
    IF ( ABS(v0(2) - v1(2)) >  height_of_torus * 0.5_wp) THEN
      ! we will wrap around + or -  length_of_torus       
      IF (v0(2) > v1(2)) THEN
        new_v1_coord%x(2) = v1(2) + height_of_torus
      ELSE
        new_v1_coord%x(2) = v1(2) - height_of_torus
      ENDIF
    ELSE
      ! keep the same 
      new_v1_coord%x(2) = v1(2)
    ENDIF

    ! this should be zero
    new_v1_coord%x(3) = v1(3)
    
  END FUNCTION plane_torus_closest_coordinates
  !--------------------------------------------------------------------

  !--------------------------------------------------------------------
  REAL(wp) FUNCTION plane_torus_distance(v0, v1, &
    & geometry_info)
    REAL(wp), INTENT(in) :: v0(3), v1(3)
    TYPE(t_grid_geometry_info), INTENT(in) :: geometry_info

    TYPE(t_cartesian_coordinates) :: dv

    dv = plane_torus_closest_coordinates(v0, v1, geometry_info)
    dv%x = dv%x - v0
    plane_torus_distance = d_norma_3d(dv)
    
  END FUNCTION plane_torus_distance
  !--------------------------------------------------------------------

  !-------------------------------------------------------------------------
  ELEMENTAL FUNCTION sphere_cartesian_midpoint (p1, p2) result(mid_point)
    TYPE(t_cartesian_coordinates), INTENT(in) :: p1, p2
    TYPE(t_cartesian_coordinates) :: mid_point

      mid_point%x = (p1%x + p2%x)
      mid_point%x = mid_point%x / &
        & d_norma_3d(mid_point)

  END FUNCTION sphere_cartesian_midpoint
  !-------------------------------------------------------------------------
   
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
  !! Computes spherical area of a triangle
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
  !! Finds the intersection of two great circles.
  !!
  !! Cannot be made elemental, as long as the
  !! subroutine calls to message and finish are maintained.
  !!
  !! @par Revision History
  !! Developed  by Luis Kornblueh  (2004).
  !!
  ELEMENTAL FUNCTION spherical_intersection (p0, p1, v0, v1) result(p)
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

  END FUNCTION spherical_intersection
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  ! Projects a point to a plane defined by three points:
  ! to_v0, to_v1, to_v2
  ELEMENTAL FUNCTION project_point_to_plane (in_point, to_v0, to_v1, to_v2) result(projected_point)
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

  END FUNCTION project_point_to_plane
  !-------------------------------------------------------------------------
    
  !-------------------------------------------------------------------------
  ELEMENTAL FUNCTION plane_angle_of_three_points(p1, p2, p3) result(angle)
    TYPE(t_cartesian_coordinates), INTENT(in) :: p1,p2,p3

    REAL(wp) :: prod,angle
    TYPE(t_cartesian_coordinates) :: v1,v2

    v1%x = p1%x - p2%x
    v2%x = p3%x - p2%x
    prod = DOT_PRODUCT(v1%x,v2%x)/ (norma(v1) * norma(v2))
    angle = ACOS(prod)

  END FUNCTION plane_angle_of_three_points
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
  SUBROUTINE spherical_intersection2 (p0, p1, v0, v1, p, nlen)
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

  END SUBROUTINE spherical_intersection2
#else
  ELEMENTAL FUNCTION spherical_intersection2 (p0, p1, v0, v1) result(p)
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

  END FUNCTION spherical_intersection2
#endif
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Computes the vector product of of 2 unit vectors in lon,lat
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
  !!  Calculates coordinates of a grid point (longin, latin) in relation to
  !!  a rotated coordinate system with origin in (reflongin, reflatin)
  !!  with rotation matrix and Euler angles. <br>
  !!  The rotation is done by two basic rotations. First around unit vector in
  !!  noth pole direction (z-direction), second around new (y-direction).
  !!  The basic roatation matrices are:
  !!  @f{align*}{
  !!  R_z (\lambda) &=
  !!  \begin{pmatrix}
  !!   \cos\lambda & \sin\lambda & 0 \\\
  !!   -\sin\lambda & \cos\lambda & 0 \\\
  !!   0 & 0 & 1
  !!  \end{pmatrix}\\\
  !!  R_y (\varphi) &=
  !!  \begin{pmatrix}
  !!   \cos\varphi & 0 & \sin\varphi \\\
  !!   0 & 1 & 0 \\\
  !!   -\sin\varphi & 0 &\cos\varphi
  !!  \end{pmatrix}
  !!  @f}
  !!  Thus, the rotation matrix @f$R@f$ is
  !!  @f{align*}{
  !!  R(\lambda,\varphi) ~=~  R_y (\varphi) \cdot R_z (\lambda) &=
  !!  \begin{pmatrix}
  !!   \cos\lambda \cos\varphi & \sin\lambda \cos\varphi & \sin\varphi\\\
  !!   -\sin\lambda & \cos\lambda & 0\\\
  !!   -\cos\lambda \sin\varphi & -\sin\lambda \sin\varphi & \cos\varphi
  !!  \end{pmatrix}
  !!  @f}
  !!  Given the Cartesian coordinates of point
  !!  @f$(\lambda, \varphi) \mapsto (x,y,z)^T@f$
  !!  and applying the rotation matrix
  !!  @f{align*}{
  !!  R(\lambda,\varphi) \cdot (x,y,z)^T &= (1,0,0)^T
  !!  @f}
  !!  this point is the origin of the rotated coordinate system. All other points
  !!  are then rotated to this new system in the similar way by applying the
  !!  rotation matrix to their Cartesian coordinates.<br>
  !!  Algorithm:
  !!  <ol>
  !!  <li>  Given the geographical coordinates of a
  !!   point @f$(\lambda_{in}, \varphi_{in})@f$ and the origin of the rotated system
  !!  @f$(\lambda_{in}^{ref}, \varphi_{in}^{ref})@f$
  !!  <li>  Get Cartesian coordinates of @f$(\lambda_{in}, \varphi_{in})@f$ by
  !!  <i>gc2cc</i>
  !!  <li>  Apply roatation matrix @f$R(\lambda_{in}^{ref}, \varphi_{in}^{ref})@f$
  !!  <li>  Get rotated geographical coordinates @f$(\lambda_{dis}, \varphi_{dis})@f$
  !!  of new Cartesian ones by <i>cc2gc</i>
  !!  </ol>
  !!
  !! @par Revision History
  !! Developed and tested by Th. Heinze (2006-07-20)
  !!
  PURE SUBROUTINE disp_new(longin, latin, reflongin, reflatin, dislong, dislat)
    
    REAL(wp), INTENT(in)   :: longin, latin, reflongin, reflatin
    REAL(wp), INTENT(out)  :: dislong, dislat
    
    TYPE(t_geographical_coordinates) :: gpoint, gpoint_rot
    TYPE(t_cartesian_coordinates)    :: cpoint, row, cpoint_rot
    
    REAL(wp)               :: z_long, z_lat, z_reflong, z_reflat
    REAL(wp)               :: z_cos_a, z_cos_b, z_sin_a, z_sin_b
    REAL(wp)               :: z_twopi
    
    z_twopi = 2._wp * pi
    
    ! check if all angles are in the wright range    
    z_long = longin
    IF (z_long < 0._wp) z_long = z_long + z_twopi
    
    z_reflong = reflongin
    IF (z_reflong < 0._wp) z_reflong = z_reflong + z_twopi
    
    z_lat    = latin
    z_reflat = reflatin
    
    ! get cartesian coordinates of (z_long, z_lat)
    
    gpoint%lon = z_long
    gpoint%lat = z_lat
    
    cpoint =  gc2cc(gpoint)
    
    ! calculate rotation matrix    
    z_cos_a = COS(z_reflong)
    z_sin_a = SIN(z_reflong)
    z_cos_b = COS(z_reflat)
    z_sin_b = SIN(z_reflat)
    
    ! first row stored in a cartesian vector (just for using cc_dot_product)    
    row%x(1) = z_cos_a * z_cos_b
    row%x(2) = z_sin_a * z_cos_b
    row%x(3) = z_sin_b
    
    ! first component of rotated cpoint    
    cpoint_rot%x(1) = cc_dot_product (row, cpoint)
    
    ! second row stored in a cartesian vector (just for using cc_dot_product)    
    row%x(1) = -1._wp * z_sin_a
    row%x(2) = z_cos_a
    row%x(3) = 0._wp
    
    ! second component of rotated cpoint    
    cpoint_rot%x(2) = cc_dot_product (row, cpoint)
    
    ! third row stored in a cartesian vector (just for using cc_dot_product)    
    row%x(1) = -1._wp * z_cos_a * z_sin_b
    row%x(2) = -1._wp * z_sin_a * z_sin_b
    row%x(3) = z_cos_b
    
    ! third component of rotated cpoint    
    cpoint_rot%x(3) = cc_dot_product (row, cpoint)
    
    ! get geographic coordinates of cpoint_rot
    gpoint_rot = cc2gc(cpoint_rot)
    
    dislong = gpoint_rot%lon
    dislat  = gpoint_rot%lat
    
  END SUBROUTINE disp_new
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !>
  !! Compute the zonal and meridional components of a vector with respect.
  !!
  !! Compute the zonal and meridional components of a vector with respect
  !! to a rotated coordinate system. This function uses the same
  !! conventions and formalism as "disp_new". The algortithm is as
  !! follows:
  !! 1) convert the vector into Cartesian coordinates
  !! 2) apply the rotation to the vector in cartesian components
  !! 3) transform the rotated vector in geographical coordinates, using
  !! the new values of latitude and longitude.
  !!
  !! @par Revision History
  !! Developed by Marco Restelli (2008-03-04)
  !!
  PURE SUBROUTINE disp_new_vect(p_gu,p_gv,lon,lat,barlon,barlat, &
    & new_lon,new_lat,new_p_gu,new_p_gv)
    
    REAL(wp), INTENT(in) :: &
      & p_gu, p_gv,       & ! original lon-lat components
      & lon, lat,         & ! original lon and lat
      & barlon, barlat,   & ! new origin
      & new_lon, new_lat    ! point (lon,lat) in the new system (this
    ! value should be computed with disp_new)
    !
    REAL(wp), INTENT(out) :: &
      & new_p_gu, new_p_gv  ! new geographical components
    !
    REAL(wp) :: p_c(3), r1(3,3), r2(3,3), rot_p_c(3)
    !
    !--------------------------------------------------------------------
    
    ! 1) convert to Cartesian coordinates
    CALL gvec2cvec(p_gu,p_gv,lon,lat,p_c(1),p_c(2),p_c(3))
    
    ! 2) apply the rotation
    r1(1,:) = (/  COS(barlon) , SIN(barlon) , 0.0_wp /)
    r1(2,:) = (/ -SIN(barlon) , COS(barlon) , 0.0_wp /)
    r1(3,:) = (/     0.0_wp   ,    0.0_wp   , 1.0_wp /)
    
    r2(1,:) = (/  COS(barlat) , 0.0_wp , SIN(barlat) /)
    r2(2,:) = (/     0.0_wp   , 1.0_wp ,    0.0_wp   /)
    r2(3,:) = (/ -SIN(barlat) , 0.0_wp , COS(barlat) /)
    
    rot_p_c = MATMUL(r2,MATMUL(r1,p_c))
    
    ! 3) back to geographical coordinates
    CALL cvec2gvec(rot_p_c(1),rot_p_c(2),rot_p_c(3), &
      & new_lon,new_lat,new_p_gu,new_p_gv)
    
  END SUBROUTINE disp_new_vect
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !>
  !! Converts cartesian coordinates to geographical.
  !!
  !!
  !! @par Revision History
  !! Developed  by Luis Kornblueh  (2004).
  !! Completely new version by Thomas Heinze (2006-07-20)
  !!
  ELEMENTAL FUNCTION cc2gc(p_x) result (p_pos)
    
    TYPE(t_cartesian_coordinates), INTENT(in) :: p_x          ! Cart. coordinates    
    TYPE(t_geographical_coordinates)          :: p_pos        ! geo. coordinates
    
    REAL(wp)                                :: z_x, z_y, z_z, z_r    
    !-----------------------------------------------------------------------
    
    z_x = p_x%x(1)
    z_y = p_x%x(2)
    z_z = p_x%x(3)
    
    z_r = z_x * z_x + z_y * z_y
    z_r = SQRT(z_r)
    
    IF (ABS(z_r) < dbl_eps) THEN    ! one of the poles
          
      IF (z_z > 0.0_wp) THEN
        p_pos%lat = pi_2
      ELSE
        p_pos%lat = -1._wp * pi_2
      END IF
      p_pos%lon = 0._wp
      
    ELSE
      
      p_pos%lat = ATAN2 ( z_z, z_r)
          
      IF (ABS(z_x) < dbl_eps) THEN    ! z_x == 0 ?
              
        IF (z_y >= 0.0_wp) THEN
          p_pos%lon = pi_2
        ELSE
          p_pos%lon = -1._wp * pi_2
        END IF
        
      ELSE
        p_pos%lon = ATAN2( z_y, z_x)
      END IF
      
    END IF
    
  END FUNCTION cc2gc
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !>
  !! Converts geographical to cartesian coordinates.
  !!
  !!
  !! @par Revision History
  !! Developed  by Luis Kornblueh  (2004).
  !!
  ELEMENTAL FUNCTION gc2cc (p_pos)  result(p_x)
    
    TYPE(t_geographical_coordinates), INTENT(in) :: p_pos     ! geo. coordinates
    
    TYPE(t_cartesian_coordinates)                :: p_x       ! Cart. coordinates
    
    REAL (wp)                                  :: z_cln, z_sln, z_clt, z_slt
    
    !-----------------------------------------------------------------------    
    z_sln = SIN(p_pos%lon)
    z_cln = COS(p_pos%lon)
    z_slt = SIN(p_pos%lat)
    z_clt = COS(p_pos%lat)
    
    p_x%x(1) = z_cln*z_clt
    p_x%x(2) = z_sln*z_clt
    p_x%x(3) = z_slt
    
  END FUNCTION gc2cc
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !>
  !!  Calculates dot product of to cartesian coordinates.
  !!
  !! @par Revision History
  !! Developed by Thomas Heinze (2006-07-05).
  !!
  ELEMENTAL FUNCTION cc_dot_product(cc_x1, cc_x2)  result (p_prod)
    !
    TYPE(t_cartesian_coordinates), INTENT(in):: cc_x1, cc_x2 ! cart. coordinates
    
    REAL(wp) :: p_prod    ! scalar product of cart.coordinates
    
    p_prod = DOT_PRODUCT (cc_x1%x, cc_x2%x)
    
  END FUNCTION cc_dot_product
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !>
  !!  Calculates 2 norm for vectors in cartesian coordinates.
  !!
  !!
  !! @par Revision History
  !! Developed by Marco Restelli (2007-11-22)
  !!
  ELEMENTAL FUNCTION cc_norm(cc_x1)  result (norm)
    !
    TYPE(t_cartesian_coordinates), INTENT(in) :: cc_x1 ! cart. coordinates
    
    REAL(wp) :: norm   
    
    norm = SQRT(cc_dot_product(cc_x1,cc_x1))
    
  END FUNCTION cc_norm
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !>
  !! Computes vector product of x0,x1.
  !!
  !!
  !! @par Revision History
  !! Developed  by Luis Kornblueh  (2004).
  !!
  ELEMENTAL FUNCTION vector_product (x0, x1) result(x2)
    
    TYPE(t_cartesian_coordinates), INTENT(in) :: x0, x1
    
    TYPE(t_cartesian_coordinates) :: x2
    
    !-----------------------------------------------------------------------    
    x2%x(1) = x0%x(2)*x1%x(3) - x0%x(3)*x1%x(2)
    x2%x(2) = x0%x(3)*x1%x(1) - x0%x(1)*x1%x(3)
    x2%x(3) = x0%x(1)*x1%x(2) - x0%x(2)*x1%x(1)
    
  END FUNCTION vector_product
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------  
  PURE FUNCTION cartesian_coordinates_plus(x,y) result(z)
    
    TYPE(t_cartesian_coordinates) :: z
    TYPE(t_cartesian_coordinates), INTENT(in) :: x,y
    
    z%x = x%x + y%x
    
  END FUNCTION cartesian_coordinates_plus
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  PURE FUNCTION cartesian_coordinates_minus(x,y) result(z)
    
    TYPE(t_cartesian_coordinates) :: z
    TYPE(t_cartesian_coordinates), INTENT(in) :: x,y
    
    z%x = x%x - y%x
    
  END FUNCTION cartesian_coordinates_minus
  !-------------------------------------------------------------------------
  
  !-----------------------------------------------------------------------  
  PURE FUNCTION cartesian_coordinates_prod(a,x) result(y)
    
    TYPE(t_cartesian_coordinates) :: y
    REAL(wp), INTENT(in) :: a
    TYPE(t_cartesian_coordinates), INTENT(in) :: x
    
    y%x = a * x%x
    
  END FUNCTION cartesian_coordinates_prod
  !-----------------------------------------------------------------------
  
  !-----------------------------------------------------------------------
  !>
  !! Calculates integral @f$\int_{\tau} \Psi@f$ over spherical triangle @f$\tau@f$.
  !!
  !! Calculates integral @f$\int_{\tau} \Psi@f$ over spherical triangle @f$\tau@f$
  !! using second order quadrature rule
  !! @f{align*}{
  !!  I(\Psi, \tau) &= \frac{1}{6} \sum_{i=1}^3 \omega_i \Psi_i
  !! @f}
  !! See Boal and Sayas (2004), page 66 for details.
  !!
  !! @par Revision History
  !! Original version by Thomas Heinze (2006-10-19).
  !!
  FUNCTION integral_over_triangle (weight, values) result(p_integ)
    
    REAL(wp), INTENT(in) :: weight(3), values(3)
    
    REAL(wp)             :: p_integ
    
    !-----------------------------------------------------------------------
    
    p_integ = DOT_PRODUCT( weight, values)
    p_integ = p_integ / 6._wp
    
  END FUNCTION integral_over_triangle
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
  !! Rotates latitude and longitude for more accurate computation.
  !!
  !! Rotates latitude and longitude for more accurate computation
  !! of bilinear interpolation
  !!
  !! @par Revision History
  !!  developed by Guenther Zaengl, 2009-03-13
  !!
  SUBROUTINE rotate_latlon( lat, lon, pollat, pollon )
    
    REAL(wp), INTENT(inout) :: lat, lon
    REAL(wp), INTENT(in)    :: pollat, pollon
        
    REAL(wp) :: rotlat, rotlon
    
    !-----------------------------------------------------------------------    
    rotlat = ASIN(SIN(lat)*SIN(pollat) + COS(lat)*COS(pollat)*COS(lon-pollon))
    rotlon = ATAN2( COS(lat)*SIN(lon-pollon) , &
      & (COS(lat)*SIN(pollat)*COS(lon-pollon)- SIN(lat)*COS(pollat)) )
    
    lat = rotlat
    lon = rotlon
    
  END SUBROUTINE rotate_latlon
  !-------------------------------------------------------------------------
  
  
  !-------------------------------------------------------------------------
  !> Rotates lon-lat grid
  !!
  !! Rotates latitude and longitude for all grid points to standard grid
  !! coordinates.
  !!
  !! @par Revision History
  !!  Initial revision                   : F. Prill,  2011-08-04
  !!  floating point exception handling  : G. Zaengl, 2012-11-20
  !!
  SUBROUTINE rotate_latlon_grid( lon_lat_grid, rotated_pts )
    
    TYPE (t_lon_lat_grid), INTENT(in)    :: lon_lat_grid
    REAL(wp),              INTENT(inout) :: rotated_pts(:,:,:)
    
    ! Local parameters
    REAL(wp) :: sincos_pole(2,2) ! (lon/lat, sin/cos)
    REAL(wp) :: sincos_lon(lon_lat_grid%lon_dim,2), &
      & sincos_lat(lon_lat_grid%lat_dim,2)
    INTEGER :: k, j
    REAL(wp) :: rlon_lat, pi_180, npole_rad(2), arg1, arg2
    
    !-----------------------------------------------------------------------
    
    ! convert north pole: degree -> rad:
    pi_180 = ATAN(1._wp)/45._wp
    npole_rad(:) = lon_lat_grid%north_pole(:)*pi_180
    
    sincos_pole(:,1) = SIN(npole_rad(:))
    sincos_pole(:,2) = COS(npole_rad(:))
    
    DO k=1,lon_lat_grid%lon_dim
      rlon_lat        = lon_lat_grid%start_corner(1) + REAL(k-1,wp)*lon_lat_grid%delta(1)
      sincos_lon(k,:) = (/ SIN(rlon_lat), COS(rlon_lat) /)
    END DO
    DO k=1,lon_lat_grid%lat_dim
      rlon_lat        = lon_lat_grid%start_corner(2) + REAL(k-1,wp)*lon_lat_grid%delta(2)
      sincos_lat(k,:) = (/ SIN(rlon_lat), COS(rlon_lat) /)
    END DO
    
    
    DO j = 1, lon_lat_grid%lat_dim
      DO k = 1, lon_lat_grid%lon_dim
        
        ! ATAN2(COS(phi)*SIN(lambda), SIN(poleY)*COS(phi)*COS(lambda) - SIN(phi)*COS(poleY)) + poleX
        arg1 = sincos_lat(j,2)*sincos_lon(k,1)
        arg2 = sincos_pole(2,1)*sincos_lat(j,2)*sincos_lon(k,2) - sincos_lat(j,1)*sincos_pole(2,2)
        
        IF (arg1 /= 0.0_wp .OR. arg2 /= 0.0_wp) THEN
          rotated_pts(k,j,1) = ATAN2( arg1, arg2 ) + npole_rad(1)
        ELSE
          rotated_pts(k,j,1) = 0.0_wp ! ATAN2(0,0) is undefined, so we just have to set something
        ENDIF
        
        ! ASIN( SIN(phi)*SIN(poleY) + COS(phi)*COS(lambda)*COS(poleY) )
        rotated_pts(k,j,2) = &
          & ASIN( sincos_lat(j,1)*sincos_pole(2,1) + &
          & sincos_lat(j,2)*sincos_lon(k,2)*sincos_pole(2,2) )
        
      ENDDO
    ENDDO
    
  END SUBROUTINE rotate_latlon_grid
  !-------------------------------------------------------------------------
  
  
  !-------------------------------------------------------------------------
  !> Compute normalized area weights for lon-lat grid
  !!
  !! @return 1D array (index: latitude)
  !!
  !! @par Revision History
  !!  developed by F. Prill, 2012-05-24
  !!
  SUBROUTINE latlon_compute_area_weights( grid, area )
    
    TYPE (t_lon_lat_grid), INTENT(in)    :: grid
    REAL(wp),              INTENT(inout) :: area(:)
    ! local variables
    REAL(wp) :: start_lat, delta_lon, delta_lat, delta_lat_2,  &
      & pi_180, radius, rr_dlon, tot_area
    REAL(wp) :: latitude(grid%lat_dim)
    INTEGER :: k, pole1, pole2
    
    radius = grid_sphere_radius ! earth's radius (average)
    pi_180 = ATAN(1._wp)/45._wp
    start_lat   = grid%reg_lat_def(1) * pi_180
    delta_lon   = grid%reg_lon_def(2) * pi_180
    delta_lat   = grid%reg_lat_def(2) * pi_180
    delta_lat_2 = delta_lat / 2._wp
    rr_dlon     = radius*radius * delta_lon
    
    ! for each latitude, compute area of a grid box with a lon-lat
    ! point at its center
    DO k = 1, grid%lat_dim
      latitude(k)  = start_lat + REAL(k-1,wp)*delta_lat
      area(k)      = 2._wp*rr_dlon * SIN(delta_lat_2) * COS(latitude(k))
    END DO
    ! treat the special case of the poles (compute area of "triangle"
    ! with one vertex at the pole and the opposite side with constant
    ! latitude)
    pole1 = INT(( 90._wp - grid%reg_lat_def(1))/grid%reg_lat_def(2)) + 1
    pole2 = INT((-90._wp - grid%reg_lat_def(1))/grid%reg_lat_def(2)) + 1
    IF ((pole1 > 0) .AND. (pole1 <= grid%lat_dim)) THEN
      area(pole1) = rr_dlon*(1._wp-SIN(latitude(pole1)-delta_lat_2))
    END IF
    IF ((pole2 > 0) .AND. (pole2 <= grid%lat_dim)) THEN
      area(pole2) = rr_dlon*(1._wp-ABS(SIN(latitude(pole2)+delta_lat_2)))
    END IF
    ! normalize
    tot_area = REAL(grid%lon_dim,wp) * SUM(area)
    DO k = 1, grid%lat_dim
      area(k) = area(k)/tot_area
    END DO
    
  END SUBROUTINE latlon_compute_area_weights
  !-----------------------------------------------------------------------
  
  !-----------------------------------------------------------------------
  !>
  !! Provides rotation angle between coordinate systems
  !!
  !! Provides sin(delta) and cos(delta), the entries of the rotation matrix
  !! for vectors from the geographical to the rotated coordinate system,
  !! of which the pole is known.
  !!
  !! @par Revision History
  !!  developed by Almut Gassmann, MPI-M, 2010-01-20
  !!
  SUBROUTINE rotate_latlon_vec( lon, lat, pollon, pollat, sin_d, cos_d )
    !
    !
    REAL(wp), INTENT(in) :: lat, lon
    REAL(wp), INTENT(in) :: pollat, pollon
    REAL(wp), INTENT(out):: sin_d, cos_d
    
    REAL(wp) :: z_lamdiff, z_a, z_b, z_sq
    !-----------------------------------------------------------------------
    
    z_lamdiff = pollon - lon
    z_a       = COS(pollat)*SIN(z_lamdiff)
    z_b       = COS(lat)*SIN(pollat)-SIN(lat)*COS(pollat)*COS(z_lamdiff)
    z_sq      = SQRT(z_a*z_a + z_b*z_b)
    sin_d     = z_a/z_sq
    cos_d     = z_b/z_sq
    
  END SUBROUTINE rotate_latlon_vec
  !-------------------------------------------------------------------------
    
  !-------------------------------------------------------------------------
  !>
  !! gnomonic projection for unit sphere.
  !!
  !! Projects a point
  !! (lat,lon) onto a tangent plane with the origin at (lat_c,lon_c).
  !! The results are the local cartesian coordinates (x,y).
  !!
  !!
  !! Features:
  !! - azimuthal
  !! - neither conformal nor equal-area
  !! - no distortion at the center only
  !! - directions from the center are true
  !! - radial scale factor increases as distances increase from the center
  !! - scale factor in a direction perpendicular to a line radiating from
  !!   the center increases as distances increase from the center
  !! - all great circles are shown as straight lines
  !!
  !! @par Revision History
  !! developed by Daniel Reinert, DWD, 2009-10-30
  !!
  SUBROUTINE gnomonic_proj( lon_c, lat_c, lon, lat, x, y )
    !
    ! !LITERATURE:
    ! Map Projections: A Working Manual, Snyder, 1987, p. 165
    !
    !
    REAL(wp), INTENT(in) :: lat_c, lon_c  ! center on tangent plane
    REAL(wp), INTENT(in) :: lat, lon      ! point to be projected
    REAL(wp), INTENT(out):: x, y          ! coordinates of projected point
    
    REAL(wp) :: zk   ! scale factor perpendicular to the radius from the
    ! center of the map
    REAL(wp) :: cosc ! cosine of the angular distance of the given point
    ! (lat,lon) from the center of projection
    
    !-----------------------------------------------------------------------
    
    cosc = SIN(lat_c)*SIN(lat) + COS(lat_c)*COS(lat)*COS(lon-lon_c)
    zk   = 1._wp/cosc
    
    x    = zk * COS(lat)*SIN(lon-lon_c)
    y    = zk * ( COS(lat_c)*SIN(lat) - SIN(lat_c)*COS(lat)*COS(lon-lon_c) )
    
  END SUBROUTINE gnomonic_proj
  !-------------------------------------------------------------------------
  
  
  !-------------------------------------------------------------------------
  !>
  !! azimuthal equidistant projection for unit sphere.
  !!
  !! Projects a point
  !! (lat,lon) onto a tangent plane with the origin at (lat_c,lon_c).
  !! The results are the local cartesian coordinates (x,y).
  !!
  !!
  !! Features:
  !! - azimuthal
  !! - neither conformal nor equal-area
  !! - no distortion at the center only
  !! - directions from the center are true
  !! - radial scale is true
  !! - scale factor in a direction perpendicular to a line radiating from
  !!   the center linearly increases as distances increase from the center
  !!
  !!
  !! @par Revision History
  !! developed by Daniel Reinert, DWD, 2011-04-27
  !!
  SUBROUTINE az_eqdist_proj( lon_c, lat_c, lon, lat, x, y )
    !
    ! !LITERATURE:
    ! Map Projections: A Working Manual, Snyder, 1987, p. 200
    !
    !
    REAL(wp), INTENT(in) :: lat_c, lon_c  ! center on tangent plane
    REAL(wp), INTENT(in) :: lat, lon      ! point to be projected
    REAL(wp), INTENT(out):: x, y          ! coordinates of projected point
    
    REAL(wp) :: zk   ! scale factor perpendicular to the radius from the
    ! center of the map
    REAL(wp) :: c    ! angular distance of the given point (lat,lon)
    ! from the center of projection
    
    !-----------------------------------------------------------------------
    
    c  = ACOS(SIN(lat_c)*SIN(lat) + COS(lat_c)*COS(lat)*COS(lon-lon_c))
    zk = c/SIN(c)
    
    x  = zk * COS(lat)*SIN(lon-lon_c)
    y  = zk * ( COS(lat_c)*SIN(lat) - SIN(lat_c)*COS(lat)*COS(lon-lon_c) )
    
  END SUBROUTINE az_eqdist_proj
  !-------------------------------------------------------------------------
  
  
  !-------------------------------------------------------------------------
  !>
  !! orthographic projection for unit sphere.
  !!
  !! Projects a point
  !! (lat,lon) onto a tangent plane with the origin at (lat_c,lon_c).
  !! The results are the local cartesian coordinates (x,y).
  !!
  !! Features:
  !! - azimuthal
  !! - neither conformal nor equal-area
  !! - no distortion at the center only
  !! - directions from the center are true
  !! - radial scale factor decreases as distances increase from the center
  !! - scale in a direction perpendicular to a line radiating from the center is true
  !!
  !! @par Revision History
  !! developed by Daniel Reinert, DWD, 2009-10-30
  !!
  SUBROUTINE orthogr_proj( lon_c, lat_c, lon, lat, x, y )
    !
    ! !LITERATURE:
    ! Map Projections: A Working Manual, Snyder, 1987,  p. 149
    !
    !
    REAL(wp), INTENT(in) :: lat_c, lon_c  ! center on tangent plane
    REAL(wp), INTENT(in) :: lat, lon      ! point to be projected
    REAL(wp), INTENT(out):: x, y          ! coordinates of projected point
    
    !-----------------------------------------------------------------------
    
    x = COS(lat)*SIN(lon-lon_c)
    y = COS(lat_c)*SIN(lat) - SIN(lat_c)*COS(lat)*COS(lon-lon_c)
    
  END SUBROUTINE orthogr_proj
  !-------------------------------------------------------------------------
  
  
  !!-----------------------------------------------------------------------
  !!>
  !! The following functions are needed for the ECHAM physics
  !!
  FUNCTION betacf(p,q,x)
    
    !! Description:
    !!
    !!  used by betai: evaluates continued fraction for incomplete
    !!  beta function by modi ed lentz's method ( x 5.2).
    !!  first step of lentz's method.
    !!
    !! Method:
    !!   See Numerical Recipes (Fortran)
    !!
    !! @author   A. Tompkins, MPI, 2000
    !!    
    USE mo_kind,      ONLY: wp
    USE mo_exception, ONLY: finish
    
    IMPLICIT NONE
    
    REAL(wp)             :: betacf
    REAL(wp), INTENT(in) :: p,q, & ! beta shape parameters
      & x      ! integration limit
    
    INTEGER :: maxit = 100, m, m2
    REAL(wp) :: zeps = 3.e-7_wp, fpmin = 1.e-30_wp, &
      & aa, c, d, del, h, qab, qam, qap
    
    qab = p+q
    
    !! these q's will be used in factors that occur in the coe cients (6.4.6).
    
    qap = p+1.0_wp
    qam = p-1.0_wp
    c = 1.0_wp
    d = 1.0_wp-qab*x/qap
    IF (ABS(d) < fpmin) d = fpmin
    d = 1.0_wp/d
    h = d
    m = 1
    del = 2.0_wp
    DO WHILE (ABS(del-1.0_wp) > zeps)
      m2 = 2*m
      aa = REAL(m,wp)*(q-REAL(m,wp))*x/((qam+REAL(m2,wp))*(p+REAL(m2,wp)))
      d = 1.0_wp+aa*d  ! one step (the even one) of the recurrence.
      IF (ABS(d) < fpmin) d = fpmin
      c = 1.0_wp+aa/c
      IF (ABS(c) < fpmin) c = fpmin
      d = 1.0_wp/d
      h = h*d*c
      aa = -(p+REAL(m,wp))*(qab+REAL(m,wp))*x/((p+REAL(m2,wp))*(qap+REAL(m2,wp)))
      d = 1.0_wp+aa*d ! next step of the recurrence (the odd one).
      IF (ABS(d) < fpmin) d = fpmin
      c = 1.0_wp+aa/c
      IF (ABS(c) < fpmin) c = fpmin
      d = 1.0_wp/d
      del = d*c
      h = h*del
      m = m+1            ! AMT
      IF (m > maxit) THEN
        CALL finish ('betacf','a or b too big, or maxit too small in betacf')
        del = 1.0_wp
      ENDIF
    ENDDO
    betacf = h
    
  END FUNCTION betacf
  !-----------------------------------------------------------------------
  
  
  !-----------------------------------------------------------------------  
  FUNCTION gammln(xx)
    
    !! Description:
    !!
    !! Gamma function calculation
    !! returns the value ln[g(xx)] for xx > 0.
    !!
    !! Method:
    !!   See Numerical Recipes
    !!
    !! @author  A. Tompkins, MPI, 2000
    !    
    USE mo_kind, ONLY: wp
    
    IMPLICIT NONE
    
    REAL(wp)             :: gammln
    REAL(wp), INTENT(in) :: xx
    
    INTEGER :: j
    
    REAL(wp) :: ser, tmp, x, y
    REAL(wp), PARAMETER :: cof(6) = (/ &
      & 76.18009172947146_wp, -86.50532032941677_wp, &
      & 24.01409824083091_wp, -1.231739572450155_wp, &
      & 0.1208650973866179e-2_wp, -0.5395239384953e-5_wp /)
    REAL(wp), PARAMETER :: stp = 2.5066282746310005_wp
    
    x = xx
    y = x
    tmp = x+5.5_wp
    tmp = (x+0.5_wp)*LOG(tmp)-tmp
    ser = 1.000000000190015_wp
    DO j =1, 6
      y = y+1.0_wp
      ser = ser+cof(j)/y
    ENDDO
    gammln = tmp+LOG(stp*ser/x)
    
  END FUNCTION gammln
  !-----------------------------------------------------------------------  
  
  
  !-----------------------------------------------------------------------  
  FUNCTION betai(p,q,x)
    
    !! Description:
    !!
    !! Uses betacf, gammln; returns the incomplete beta function I x (a; b).
    !!
    !! Method:
    !!   See Numerical Recipes (Fortran)
    !!
    !! @author   A. Tompkins, MPI, 2000
    !!    
    USE mo_kind, ONLY: wp
    
    IMPLICIT NONE
    
    REAL(wp)             :: betai
    REAL(wp), INTENT(in) :: p,q, & ! beta shape parameters
      & x      ! integration limit
    !  local scalars:
    REAL(wp) :: bt
    !  REAL FUNCTION (wp):: betacf,gammln
    
    IF (x > 0.0_wp .AND. x < 1.0_wp ) THEN  ! factors in front of the continued fraction.
      bt = EXP(gammln(p+q)-gammln(p)-gammln(q)+p*LOG(x)+q*LOG(1.0_wp-x))
    ELSE
      bt = 0.0_wp
    ENDIF
    IF (x < (p+1.0_wp)/(p+q+2.0_wp)) THEN ! use continued fraction directly.
      betai = bt*betacf(p,q,x)/p
    ELSE     ! use continued fraction after making the symmetry transformation.
      betai = 1.0_wp-bt*betacf(q,p,1.0_wp-x)/q !
    ENDIF
    
  END FUNCTION betai
  !-----------------------------------------------------------------------  
  
  !-----------------------------------------------------------------------
  !>
  !! Description:
  !!  Gamma-function from Numerical Recipes (F77)
  !! Method:
  !!
  FUNCTION gamma_fct(x)
    
    USE mo_kind, ONLY: wp
    
    IMPLICIT NONE
    
    REAL (wp):: gamma_fct
    
    REAL (wp):: cof(6) = (/76.18009173_wp, -86.50532033_wp, &
      & 24.01409822_wp, -1.231739516_wp, &
      & 0.120858003E-2_wp, -0.536382E-5_wp/)
    REAL (wp):: stp=2.50662827465_wp,                           &
      & x, xx, tmp, ser, gamma
    
    INTEGER ::  j
    
    xx  = x  - 1.0_wp
    tmp = xx + 5.5_wp
    tmp = (xx + 0.5_wp) * LOG(tmp) - tmp
    ser = 1.0_wp
    DO j = 1, 6
      xx  = xx  + 1.0_wp
      ser = ser + cof(j) / xx
    ENDDO
    gamma = tmp + LOG(stp*ser)
    gamma = EXP(gamma)
    
    gamma_fct = gamma
    
  END FUNCTION gamma_fct
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !>
  !! Calculates the domain mean of the characteristical
  !! length scale of the area. Needed for physics.
  !!
  !! @par Revision History
  !! Implemented by Kristina Froehlich, DWD (2010-10-29).
  !! moved to a more general place, Kristina Froehlich, MPI-M (2011-10-06)
  !! Changed to use total number odf cells insted of root/level by LL, MPI-M (2012-12)
  !! Not used since the charecteristic lenntgh os now part of the grid_geometry_info. LL, MPI-M (2012-12)
!   SUBROUTINE sphere_cell_mean_char_length( total_number_of_cells, mean_charlen ) ! output
!     
!     INTEGER , INTENT(in)  :: total_number_of_cells 
!     REAL(wp), INTENT(out) :: mean_charlen
!     
!     mean_charlen = SQRT (4._wp*pi*grid_sphere_radius**2 /REAL(total_number_of_cells,wp))
!     
!   END SUBROUTINE sphere_cell_mean_char_length
  !-------------------------------------------------------------------------
  
  
  !-------------------------------------------------------------------------
  !>
  !! Checking turn when travelling along three points
  !!
  !! Given three points in a cartesian system, we want to know
  !! whether, in travelling from the first to the second to the third
  !! we turn counterclockwise or clockwise.
  !! Can be used to check whether two line segments intersect, or not.
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert  (2012-04-03)
  !!
  !! @par LITERATURE
  !! Sedgewick, R. (1988): Algorithms, 2nd edition, pp. 350
  !!
  ELEMENTAL FUNCTION ccw( p0, p1, p2 )
    !
    
    IMPLICIT NONE
    
    TYPE(t_geographical_coordinates), INTENT(in) :: p0
    TYPE(t_geographical_coordinates), INTENT(in) :: p1
    TYPE(t_geographical_coordinates), INTENT(in) :: p2
    
    INTEGER :: ccw
    
    REAL(wp) :: dx1, dx2, dy1, dy2  ! segment lengths in x and y direction
    
    REAL(wp) :: dx1dy2, dy1dx2
    
    LOGICAL :: lccw
    
    !-----------------------------------------------------------------------    
    ! segment lengths in x and y direction P0-->P1
    dx1 = p1%lon - p0%lon
    dy1 = p1%lat - p0%lat
    
    ! segment lengths in x and y direction P0-->P2
    dx2 = p2%lon - p0%lon
    dy2 = p2%lat - p0%lat
    
    dx1dy2 = dx1 * dy2
    dy1dx2 = dy1 * dx2
    
    
    ! in this case we turn counterclockwise
    ! dy2/dx2 > dy1/dx1
    lccw = dx1dy2 > dy1dx2
    
    ! we set ccw=1 when turning counterclockwise and -1 when
    ! turning clockwise.
    ! The case dy2/dx2 = dy1/dx1 is neglected. In this case
    ! we set ccw = -1
    ccw = MERGE(1, -1, lccw)
    
  END FUNCTION ccw
  !-------------------------------------------------------------------------
  
  
  !-------------------------------------------------------------------------
  !>
  !! Checks, whether two lines intersect
  !!
  !! Checks whether two lines intersect (2D cartesian geometry)
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert  (2012-04-05)
  !!
  !! @par LITERATURE
  !! Sedgewick, R. (1988): Algorithms, 2nd edition, pp. 351
  !!
  ELEMENTAL FUNCTION lintersect( line1, line2 )
        
    TYPE(t_line), INTENT(in) :: line1
    TYPE(t_line), INTENT(in) :: line2
    
    INTEGER :: intersect1, intersect2
    
    LOGICAL :: lintersect
    
    !-----------------------------------------------------------------------    
    intersect1 = ccw(line1%p1,line1%p2,line2%p1)      &
      & * ccw(line1%p1,line1%p2,line2%p2)
    
    intersect2 = ccw(line2%p1,line2%p2,line1%p1)      &
      & * ccw(line2%p1,line2%p2,line1%p2)
    
    lintersect = (intersect1 + intersect2) == -2
    
  END FUNCTION lintersect
  !-------------------------------------------------------------------------
    
  
  !-------------------------------------------------------------------------
  !>
  !! Computes intersection point of two lines in 2D
  !!
  !! Computes intersection point of two lines in 2D (cartesian geometry).
  !! Note that this function does not check, whether the two lines
  !! do intersect at all. So it is up to the user to make sure that
  !! the two lines indeed do intersect (e.g. by making use of the
  !! ccw-function)
  !!
  !! The two lines are given in the from:
  !! f1(x) = y_a1 + m1*(x-x_a1)
  !! f2(x) = y_b1 + m2*(x-x_b1)
  !! with
  !! m = (y_2-y_1)/(x_2-x_1)
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert  (2012-04-03)
  !!
  FUNCTION line_intersect( line1, line2 ) result(intersect)
    
    TYPE(t_line), INTENT(in) :: line1
    TYPE(t_line), INTENT(in) :: line2
    
    REAL(wp) :: m1, m2          !< slopes    
    REAL(wp) :: intersect(2)    ! coordinates of intersection point 
    
    !-----------------------------------------------------------------------
    
    ! determine slopes of the two lines
    m1 = (line1%p2%lat - line1%p1%lat)/(line1%p2%lon - line1%p1%lon)
    m2 = (line2%p2%lat - line2%p1%lat)/(line2%p2%lon - line2%p1%lon)
    
    intersect(1) = (line2%p1%lat - line1%p1%lat + m1*line1%p1%lon - m2*line2%p1%lon) &
      & / (m1 - m2)
    
    intersect(2) = line1%p1%lat + m1*(intersect(1) - line1%p1%lon)
    
  END FUNCTION line_intersect
  !-------------------------------------------------------------------------
  
END MODULE mo_math_utilities

