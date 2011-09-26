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
!!  Modification by Daniel Reinert, DWD, (2009-07-20)
!!  - added subroutine qrdec
!!  Modification by Daniel Reinert, DWD, (2009-10-30)
!!  - added subroutine gnomonic_proj which uses a gnomonic projection in
!!    order to project a point (lat_1,lon_1) onto a tangent plane with
!!    the origin at (lat_0,lon_0). The results are the local cartesian
!!    coordinates (x_1,y_1)
!!  - added subroutine orthogr_proj. performs orthographic projection of
!!    point (lat_1,lon_1) onto a tangent plane with the origin at (lat_0,lon_0).
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
!

!

  USE mo_kind,                ONLY: wp
  USE mo_math_constants
  USE mo_physical_constants
  USE mo_exception,           ONLY: message, finish
  USE mo_parallel_config,  ONLY: nproma
  
  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  PUBLIC :: t_cartesian_coordinates
  PUBLIC :: t_geographical_coordinates
  PUBLIC :: t_lon_lat_grid
  PUBLIC :: gvec2cvec
  PUBLIC :: cvec2gvec
  PUBLIC :: arc_length
  PUBLIC :: arc_length_v
  PUBLIC :: solve_lu
  PUBLIC :: solve_chol
  PUBLIC :: solve_chol_v
  PUBLIC :: ludec
  PUBLIC :: choldec
  PUBLIC :: choldec_v
  PUBLIC :: inv_mat
  PUBLIC :: disp_new
  PUBLIC :: disp_new_vect
  PUBLIC :: cc2gc
  PUBLIC :: gc2cc
  PUBLIC :: cc_dot_product
  PUBLIC :: cc_norm
  PUBLIC :: vector_product
  PUBLIC :: integral_over_triangle
  PUBLIC :: operator(+)
  PUBLIC :: operator(*)
  PUBLIC :: rotate_latlon
  PUBLIC :: rotate_latlon_vec
  PUBLIC :: rotate_latlon_grid
  PUBLIC :: qrdec
  PUBLIC :: gnomonic_proj
  PUBLIC :: orthogr_proj
  PUBLIC :: az_eqdist_proj
  PUBLIC :: gamma_fct

! ! cartesian coordinate class

  TYPE t_cartesian_coordinates
    REAL(wp) :: x(3)
  END TYPE t_cartesian_coordinates

! ! geographical coordinate class

  TYPE t_geographical_coordinates
    REAL(wp) :: lon
    REAL(wp) :: lat
  END TYPE t_geographical_coordinates

! ! specification of (rotated) lon-lat grid

  TYPE t_lon_lat_grid
    ! grid points: start_corner(lon,lat) + delta(lon,lat)*[0,1,2, ..., dim(lon,lat)-1]
    REAL(wp) ::                     &
      &  delta       (2),           &     ! lon-lat grid resolution,                unit:rad
      &  start_corner(2),           &     ! south western corner of area (lon/lat), unit:rad
      &  poleN       (2)                  ! position of north pole (lon,lat),       unit:rad
    INTEGER  :: dimen(2)                  ! grid dimensions

    ! computed from above values:
    INTEGER  :: total_dim              ! total number of grid points 
    INTEGER  :: nblks, npromz          ! blocking info

  END TYPE t_lon_lat_grid

  INTERFACE OPERATOR(+)
    MODULE PROCEDURE cartesian_coordinates_plus
  END INTERFACE
  INTERFACE OPERATOR(*)
    MODULE PROCEDURE cartesian_coordinates_prod
  END INTERFACE

  CONTAINS

!-------------------------------------------------------------------------
!
!-------------------------------------------------------------------------
!
!
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
!
!
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
!
!
  !>
  !! Computes length of geodesic arc with endpoints @f$p\_x, p\_y@f$.
  !!
  !!
  !! @par Revision History
  !! Developed by Th.Heinze (2006-09-19).
  !! Previous version by Luis Kornblueh (2004) discarded.
  !!
  ELEMENTAL FUNCTION arc_length (p_x, p_y)  RESULT (p_arc)
  TYPE(t_cartesian_coordinates), INTENT(IN) :: p_x, p_y  ! endpoints

  REAL(wp)            :: p_arc          ! length of geodesic arc

  REAL(wp)            :: z_lx,  z_ly    ! length of vector p_x and p_y
  REAL(wp)            :: z_cc           ! cos of angle between endpoints

!-----------------------------------------------------------------------

  z_lx = SQRT(DOT_PRODUCT(p_x%x,p_x%x))
  z_ly = SQRT(DOT_PRODUCT(p_y%x,p_y%x))

  z_cc = DOT_PRODUCT(p_x%x, p_y%x)/(z_lx*z_ly)

! in case we get numerically incorrect solutions

  IF (z_cc > 1._wp )  z_cc =  1._wp
  IF (z_cc < -1._wp ) z_cc = -1._wp

  p_arc = ACOS(z_cc)

  END FUNCTION arc_length

!-------------------------------------------------------------------------
!
!
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
  PURE FUNCTION arc_length_v (p_x, p_y)  RESULT (p_arc)
  REAL(wp), INTENT(IN) :: p_x(3), p_y(3)  ! endpoints

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
!
!
  !>
  !!  Home made routine for LU decomposition.
  !!
  !!  Home made routine for LU decomposition
  !!  to be substituted by call from linear algebra library;
  !!  provides input matrix for routine <i>solve</i>;<br>
  !!  @f$k\_dim@f$  : matrix dimension <br>
  !!  @f$p\_a@f$    : input matrix, overwritten by lu decomposition <br>
  !!  @f$k\_pivot@f$: integer index array recording row pivot line exchanges
  !!
  !! @par Revision History
  !! Developed  by P.Korn (2006).
  !! Revised by Th.Heinze, DWD, (2006-11-01):
  !! - changed stop criterion from 1.e-08 to 1.e-12
  !!
  SUBROUTINE ludec (k_dim, p_a, k_pivot)
!
  INTEGER, INTENT(in)     :: k_dim              ! matrix dimension

  REAL(wp), INTENT(inout) :: p_a(k_dim,k_dim)   ! input matrix

  INTEGER, INTENT(out)    :: k_pivot(k_dim)     ! pivoting index

  REAL(wp)                :: z_rmax(k_dim)      ! reciprocal of max value for
                                                ! each row
  REAL(wp)                :: z_temp             ! temporary variable
  REAL(wp)                :: z_max              ! max value
  REAL(wp)                :: z_sum, z_dum       ! variables generated by matrix
                                                ! entries
  INTEGER                 :: ji, jj, jk         ! integer over dimension
  INTEGER                 :: imax

!-----------------------------------------------------------------------

  imax = 0

  DO ji = 1, k_dim
    z_max = 0._wp
    DO jj = 1, k_dim
      z_temp = ABS(p_a(ji,jj))
      IF(z_temp >= z_max) THEN
        z_max = z_temp
      ENDIF
    ENDDO

    IF (z_max <= 1.e-12_wp) THEN
      CALL message ('mo_math_utilities:ludec', &
                   & 'error in matrix inversion, nearly singular matrix')
      CALL finish  ('mo_math_utilities:ludec', &
                   &'error in matrix inversion, nearly singular matrix')
    ENDIF

    z_rmax(ji) = 1._wp / z_max
  ENDDO

  DO jj = 1, k_dim
    DO ji = 1, jj-1
      z_sum = p_a(ji,jj)
      DO jk = 1, ji-1
        z_sum = z_sum - p_a(ji,jk) * p_a(jk,jj)
      ENDDO
      p_a(ji,jj) = z_sum
    ENDDO

    z_max = 0._wp

    DO ji = jj, k_dim
      z_sum = p_a(ji,jj)
      DO jk = 1, jj-1
        z_sum = z_sum - p_a(ji,jk) * p_a(jk,jj)
      ENDDO
        p_a(ji,jj) = z_sum
        z_dum = z_rmax(ji) * ABS(z_sum)
        IF (z_dum >= z_max) THEN
          z_max = z_dum
          imax  = ji
        ENDIF
    ENDDO

    IF (jj /= imax) THEN

      DO jk=1,k_dim
        z_dum        = p_a(imax,jk)
        p_a(imax,jk) = p_a(jj,jk)
        p_a(jj,jk)   = z_dum
      ENDDO
      z_rmax(imax) = z_rmax(jj)

    ENDIF

    k_pivot(jj) = imax

    IF (ABS(p_a(jj,jj)) <= 1.e-12_wp) THEN
      CALL message ('mo_math_utilities:ludec', &
                   & 'error  in matrix inversion, nearly singular matrix 2')
      CALL finish  ('mo_math_utilities:ludec', &
                   & 'error  in matrix inversion, nearly singular matrix 2')
    ENDIF

    IF (jj /= k_dim) THEN
      z_dum = 1._wp / p_a(jj,jj)
      DO ji = jj + 1, k_dim
        p_a(ji,jj) = p_a(ji,jj) * z_dum
      ENDDO
    ENDIF
  ENDDO

  END SUBROUTINE ludec


!-------------------------------------------------------------------------
!
!
  !>
  !!  QR decomposition of (m x n) matrix  Applies modified Gram-Schmidt algorithm.
  !!
  !!  QR decomposition of (m x n) matrix  Applies modified Gram-Schmidt algorithm
  !!  @f$p\_inmat@f$ : input matrix (m x n)<br>
  !!  @f$p\_qmat@f$  : Q matrix (orthogonal) (m x n) <br>
  !!  @f$p\_rmat@f$  : R matrix (upper triangular) (n x n)
  !!
  !! @par Revision History
  !! Developed and tested by Daniel Reinert, DWD (2009-11-26)
  !!
  SUBROUTINE qrdec (nrow, ncol, i_startidx, i_endidx, p_inmat, p_qmat, p_rmat)
!
! !LITERATURE
! Ruenger et al. (2005): Comparison of different parallel modified
!                        Gram-Schmidt Algorithms. Lecture Notes in
!                        Computer Science, 3648, 826-836
!
  INTEGER, INTENT(in)  :: nrow          ! number of rows (m) of matrix
                                        ! p\_inmat

  INTEGER, INTENT(in)  :: ncol          ! number of columns (n) of matrix
                                        ! p\_inmat

  INTEGER, INTENT(in)  :: i_startidx,  & ! start and and indices for
                          i_endidx       ! loop over cells

  REAL(wp), INTENT(in) :: p_inmat(:,:,:)   ! matrix for which QR-
                                         ! factorization is calculated

  REAL(wp), INTENT(out)  :: p_qmat(nproma,nrow,ncol)  ! Q matrix
  REAL(wp), INTENT(out)  :: p_rmat(nproma,ncol,ncol)  ! R matrix

  INTEGER  :: i, k, jc                   ! loop indices
  REAL(wp) :: z_q(nproma,nrow)          ! column of input matrix
  REAL(wp) :: d_product(nproma)

!-----------------------------------------------------------------------

  !
  ! Modified Gram-Schmidt algorithm
  !
  p_qmat(i_startidx:i_endidx,:,:) = p_inmat(i_startidx:i_endidx,:,:)

  DO i=1, ncol

    DO jc = i_startidx, i_endidx
      z_q(jc,1:nrow)      = p_qmat(jc,1:nrow,i)                          ! ith column
      p_rmat(jc,i,i) = SQRT(DOT_PRODUCT(z_q(jc,1:nrow),z_q(jc,1:nrow)))  ! vector length
      z_q(jc,1:nrow)      = z_q(jc,1:nrow)/p_rmat(jc,i,i)                ! normalization
      p_qmat(jc,1:nrow,i) = z_q(jc,1:nrow)

      ! loop over all old columns of the input array
      DO k=i+1,ncol
        d_product(jc) = dot_product(p_qmat(jc,1:nrow,k),z_q(jc,1:nrow))
        p_qmat(jc,1:nrow,k) = p_qmat(jc,1:nrow,k) - d_product(jc) * z_q(jc,1:nrow)
        p_rmat(jc,i,k) = d_product(jc)
      ENDDO
    ENDDO

  ENDDO

  END SUBROUTINE qrdec


!-------------------------------------------------------------------------
!
!
  !>
  !! Given a positive-definite symmetric matrix @f$p\_a(1:k\_dim,1:k\_dim)@f$, with.
  !!
  !! Given a positive-definite symmetric matrix @f$p\_a(1:k\_dim,1:k\_dim)@f$, with
  !! physical dimension @f$k\_dim@f$, this routine constructs its Cholesky
  !! decomposition, @f$p\_a = LL^T@f$. On input, only the upper triangle of @f$p\_a@f$
  !! need be given; it is not modified. The Cholesky factor @f$L@f$ is returned
  !! in the lower triangle of @f$p\_a@f$, except for its diagonal elements which are
  !! returned in @f$p\_diag(1:k\_dim)@f$.
  !!
  !! @par Revision History
  !!  Original version by Tobias Ruppert and Thomas Heinze, DWD (2006-11-14)
  !!
  SUBROUTINE choldec (k_dim, p_a, p_diag)
!

  INTEGER, INTENT(in)       :: k_dim                ! matrix dimension

  REAL(wp), INTENT(inout)   :: p_a(:,:)             ! input matrix

  REAL(wp), INTENT(out)     :: p_diag(:)            ! diagonal of output matrix

  REAL(wp)                  :: z_sum                ! variable generated by
                                                    ! matrix entries
  INTEGER                   :: ji, jj, jk           ! integer over dimension

!-----------------------------------------------------------------------

  DO ji = 1, k_dim
    DO jj = ji, k_dim
      z_sum = p_a(ji,jj)

      DO jk = ji-1, 1, -1
        z_sum = z_sum - p_a(ji,jk) * p_a(jj,jk)
      ENDDO

      IF (ji == jj) THEN
        IF (z_sum <= 1.e-12_wp) THEN
          CALL message ('mo_math_utilities:choldec',                           &
                      & 'error in matrix inversion, nearly singular matrix')
          CALL finish  ('mo_math_utilities:choldec',                           &
                      & 'error in matrix inversion, nearly singular matrix')
        ELSE
          p_diag(ji) = SQRT(z_sum)
        ENDIF
      ELSE
        p_a(jj,ji) = z_sum / p_diag(ji)
      ENDIF
    ENDDO
  ENDDO

  END SUBROUTINE choldec
!-------------------------------------------------------------------------
!
!
!>
!! Vectorized version of choldec.
!!
!! Vectorized version of choldec
!! Given a positive-definite symmetric matrix @f$p\_a(1:k\_dim,1:k\_dim)@f$, with
!! physical dimension @f$k\_dim@f$, this routine constructs its Cholesky
!! decomposition, @f$p\_a = LL^T@f$. On input, only the upper triangle of @f$p\_a@f$
!! need be given; it is not modified. The Cholesky factor @f$L@f$ is returned
!! in the lower triangle of @f$p\_a@f$, except for its diagonal elements which are
!! returned in @f$p\_diag(1:k\_dim)@f$.
!!
!! @par Revision History
!!  Original version by Tobias Ruppert and Thomas Heinze, DWD (2006-11-14)
!!  Vectorized version by Guenther Zaengl, DWD (2009-04-20)
!!
#ifdef __SX__
SUBROUTINE choldec_v (istart, iend, k_dim, maxdim, p_a, p_diag)
#else
SUBROUTINE choldec_v (istart, iend, k_dim, p_a, p_diag)
#endif
!

  INTEGER, INTENT(in)       :: istart, iend         ! start and end index

  INTEGER, INTENT(in)       :: k_dim(:)                ! matrix dimension

#ifdef __SX__
  INTEGER, INTENT(in)       :: maxdim                  ! maximum matrix dimension
#endif
  REAL(wp), INTENT(inout)   :: p_a(:,:,:)             ! input matrix

  REAL(wp), INTENT(out)     :: p_diag(:,:)            ! diagonal of output matrix

#ifdef __SX__
  REAL(wp)                  :: z_sum(nproma)        ! variable generated by
                                                    ! matrix entries
#else
  REAL(wp)                  :: z_sum                ! variable generated by
                                                    ! matrix entries
#endif

  INTEGER                   :: jc, ji, jj, jk           ! integer over dimension

!-----------------------------------------------------------------------

#ifdef __SX__
DO ji = 1, maxdim
  DO jj = ji, maxdim
    DO jc = istart, iend
      IF ((ji > k_dim(jc)) .OR. (jj > k_dim(jc))) CYCLE
      z_sum(jc) = p_a(jc,ji,jj)
    ENDDO
    DO jk = ji-1, 1, -1
      DO jc = istart, iend
        IF ((ji > k_dim(jc)) .OR. (jj > k_dim(jc))) CYCLE
        z_sum(jc) = z_sum(jc) - p_a(jc,ji,jk) * p_a(jc,jj,jk)
      ENDDO
    ENDDO

    IF (ji == jj) THEN
      DO jc = istart, iend
        IF (ji > k_dim(jc)) CYCLE
        p_diag(jc,ji) = SQRT(z_sum(jc))
      ENDDO
    ELSE
      DO jc = istart, iend
        IF ((ji > k_dim(jc)) .OR. (jj > k_dim(jc))) CYCLE
        p_a(jc,jj,ji) = z_sum(jc) / p_diag(jc,ji)
      ENDDO
    ENDIF
  ENDDO
ENDDO
#else
DO jc = istart, iend
  DO ji = 1, k_dim(jc)
    DO jj = ji, k_dim(jc)
      z_sum = p_a(jc,ji,jj)

      DO jk = ji-1, 1, -1
        z_sum = z_sum - p_a(jc,ji,jk) * p_a(jc,jj,jk)
      ENDDO

      IF (ji == jj) THEN
        IF (z_sum <= 1.e-12_wp) THEN
          WRITE (*,*) p_a
          CALL message ('mo_math_utilities:choldec',                           &
                      & 'error in matrix inversion, nearly singular matrix')
          CALL finish  ('mo_math_utilities:choldec',                           &
                      & 'error in matrix inversion, nearly singular matrix')
        ELSE
          p_diag(jc,ji) = SQRT(z_sum)
        ENDIF
      ELSE
        p_a(jc,jj,ji) = z_sum / p_diag(jc,ji)
      ENDIF
    ENDDO
  ENDDO
ENDDO
#endif

END SUBROUTINE choldec_v

!-------------------------------------------------------------------------
!
!
  !>
  !! Compute the inverse of matrix a (use only for VERY small matrices!).
  !!
  !!
  !! @par Revision History
  !! Original version by Marco Restelli (2007-11-22)
  !!
  SUBROUTINE inv_mat (a)
!

  REAL(wp), INTENT(inout)   :: a(:,:)     ! input matrix

  INTEGER :: k_pivot(size(a,1)), j
  REAL(wp) :: b(size(a,1),size(a,1)), e(size(a,1))
!-----------------------------------------------------------------------

   b = a
   CALL ludec(SIZE(a,1),b,k_pivot)

   DO j=1,size(a,1)
     e = 0.0_wp
     e(j) = 1.0_wp
     CALL solve_lu(SIZE(a,1),b,e,a(:,j),k_pivot)
   ENDDO

  END SUBROUTINE inv_mat

!-------------------------------------------------------------------------
!
!
  !>
  !! Solves the set of @f$k\_dim@f$ linear equations @f$p\_a \cdot p\_x = p\_b@f$, where.
  !!
  !! Solves the set of @f$k\_dim@f$ linear equations @f$p\_a \cdot p\_x = p\_b@f$, where
  !! @f$p\_a@f$ is a positive-definite symmetric matrix with physical dimension @f$k\_dim@f$.
  !! @f$p\_a@f$ and @f$p\_diag@f$ are input as the output of the routine {\\tt choldec}.
  !! Only the lower triangle of @f$p\_a@f$ is accessed. @f$p\_b(1:k\_dim)@f$ is input as
  !! the right-hand side vector. The solution vector is returned in @f$p\_x(1:k\_dim)@f$.
  !! @f$p\_a@f$ and @f$p\_diag@f$ are not modified and can be left in place for successive
  !! calls with different right-hand sides @f$p\_b@f$. @f$p\_b@f$ is not modified unless
  !! you identify @f$p\_b@f$ and @f$p\_x@f$ in the calling sequence, which is allowed.
  !!
  !! @par Revision History
  !!  Original version by Tobias Ruppert and Thomas Heinze, DWD (2006-11-14)
  !!
  SUBROUTINE solve_chol (k_dim, p_a, p_diag, p_b, p_x)
!
  INTEGER,  INTENT(in)   :: k_dim         ! matrix dimension
  REAL(wp), INTENT(in)   :: p_a(:,:)      ! (k_dim,k_dim) matrix
  REAL(wp), INTENT(in)   :: p_diag(:)     ! (k_dim) diagonal of matrix

  REAL(wp), INTENT(inout):: p_b(:)        ! (k_dim) input vector

  REAL(wp), INTENT(inout):: p_x(:)        ! (k_dim) output vector

  REAL(wp)               :: z_sum         ! variable generated by
                                          ! matrix entries
  INTEGER                :: ji, jj        ! integer over dimensions

!-----------------------------------------------------------------------

! Solve L*y = b, storing y in x.

    DO ji = 1, k_dim
      z_sum = p_b(ji)
      DO jj = ji - 1, 1, -1
        z_sum = z_sum - p_a(ji,jj) * p_x(jj)
      ENDDO
      p_x(ji) = z_sum / p_diag(ji)
    ENDDO

! Solve L^T * x = y

    DO ji = k_dim, 1, -1
      z_sum = p_x(ji)
      DO jj = ji+1, k_dim
        z_sum = z_sum - p_a(jj,ji) * p_x(jj)
      ENDDO
      p_x(ji) = z_sum / p_diag(ji)
    ENDDO

  END SUBROUTINE solve_chol

!-------------------------------------------------------------------------
!
!
!>
!! Solves the set of @f$k\_dim@f$ linear equations @f$p\_a \cdot p\_x = p\_b@f$, where.
!!
!! Solves the set of @f$k\_dim@f$ linear equations @f$p\_a \cdot p\_x = p\_b@f$, where
!! @f$p\_a@f$ is a positive-definite symmetric matrix with physical dimension @f$k\_dim@f$.
!! @f$p\_a@f$ and @f$p\_diag@f$ are input as the output of the routine {\\tt choldec}.
!! Only the lower triangle of @f$p\_a@f$ is accessed. @f$p\_b(1:k\_dim)@f$ is input as
!! the right-hand side vector. The solution vector is returned in @f$p\_x(1:k\_dim)@f$.
!! @f$p\_a@f$ and @f$p\_diag@f$ are not modified and can be left in place for successive
!! calls with different right-hand sides @f$p\_b@f$. @f$p\_b@f$ is not modified unless
!! you identify @f$p\_b@f$ and @f$p\_x@f$ in the calling sequence, which is allowed.
!!
!! @par Revision History
!!  Original version by Tobias Ruppert and Thomas Heinze, DWD (2006-11-14)
!!
#ifdef __SX__
SUBROUTINE solve_chol_v (istart, iend, k_dim, maxdim, p_a, p_diag, p_b, p_x)
#else
SUBROUTINE solve_chol_v (istart, iend, k_dim, p_a, p_diag, p_b, p_x)
#endif
!
  INTEGER,  INTENT(in)   :: istart, iend    ! start and end indices
  INTEGER,  INTENT(in)   :: k_dim(:)        ! (nproma) matrix dimension
#ifdef __SX__
  INTEGER,  INTENT(in)   :: maxdim          ! maximum matrix dimension
#endif
  REAL(wp), INTENT(in)   :: p_a(:,:,:)      ! (nproma,k_dim,k_dim) matrix
  REAL(wp), INTENT(in)   :: p_diag(:,:)     ! (nproma,k_dim) diagonal of matrix

  REAL(wp), INTENT(inout):: p_b(:,:)        ! (nproma,k_dim) input vector

  REAL(wp), INTENT(inout):: p_x(:,:)        ! (k_dim,nproma) output vector

#ifdef __SX__
  REAL(wp)               :: z_sum(nproma) ! variable generated by
                                          ! matrix entries
#else
  REAL(wp)               :: z_sum         ! variable generated by
                                          ! matrix entries
#endif
  INTEGER                :: jc, ji, jj    ! integer over dimensions

!-----------------------------------------------------------------------

#ifdef __SX__
! Solve L*y = b, storing y in x.

  DO ji = 1,maxdim
    DO jc = istart, iend
      IF (ji > k_dim(jc)) CYCLE
      z_sum(jc) = p_b(jc,ji)
    ENDDO
    DO jj = ji - 1, 1, -1
      DO jc = istart, iend
        IF (ji > k_dim(jc)) CYCLE
        z_sum(jc) = z_sum(jc) - p_a(jc,ji,jj) * p_x(jj,jc)
      ENDDO
    ENDDO
    DO jc = istart, iend
      IF (ji > k_dim(jc)) CYCLE
      p_x(ji,jc) = z_sum(jc) / p_diag(jc,ji)
    ENDDO
  ENDDO

! Solve L^T * x = y

  DO ji = maxdim, 1, -1
    DO jc = istart, iend
      IF (ji > k_dim(jc)) CYCLE
      z_sum(jc) = p_x(ji,jc)
    ENDDO
    DO jj = ji+1, maxdim
      DO jc = istart, iend
        IF ((ji > k_dim(jc)) .OR. (jj > k_dim(jc))) CYCLE
        z_sum(jc) = z_sum(jc) - p_a(jc,jj,ji) * p_x(jj,jc)
      ENDDO
    ENDDO
    DO jc = istart, iend
      IF (ji > k_dim(jc)) CYCLE
      p_x(ji,jc) = z_sum(jc) / p_diag(jc,ji)
    ENDDO
  ENDDO
#else
! Solve L*y = b, storing y in x.

  DO jc = istart, iend
    DO ji = 1, k_dim(jc)
      z_sum = p_b(jc,ji)
      DO jj = ji - 1, 1, -1
        z_sum = z_sum - p_a(jc,ji,jj) * p_x(jj,jc)
      ENDDO
      p_x(ji,jc) = z_sum / p_diag(jc,ji)
    ENDDO

! Solve L^T * x = y

    DO ji = k_dim(jc), 1, -1
      z_sum = p_x(ji,jc)
      DO jj = ji+1, k_dim(jc)
        z_sum = z_sum - p_a(jc,jj,ji) * p_x(jj,jc)
      ENDDO
      p_x(ji,jc) = z_sum / p_diag(jc,ji)
    ENDDO
  ENDDO
#endif

END SUBROUTINE solve_chol_v

!-------------------------------------------------------------------------
!
!
  !>
  !!  home made routine for backsubstitution.
  !!
  !!  home made routine for backsubstitution
  !!  to be substituted by call from linear algebra library
  !!  dim: matrix dimension
  !!  a: input matrix  lu decomposed by ludec (otherwise it does not work)
  !!  b: right hand side, row pivot  exchanges to be performed using ind
  !!  ind: integer index array recording row pivot  exchanges
  !!
  SUBROUTINE solve_lu (k_dim, p_a, p_b, p_x, k_pivot)
!
  INTEGER,  INTENT(in)    :: k_dim                   ! matrix dimension
  INTEGER,  INTENT(in)    :: k_pivot(k_dim)          ! pivoting index
  REAL(wp), INTENT(in)    :: p_a(k_dim,k_dim)        ! matrix

  REAL(wp), INTENT(inout) :: p_b(k_dim)              ! input vector

  REAL(wp), INTENT(out)   :: p_x(k_dim)              ! output vector

  INTEGER                 :: ji, jj, ip              ! integer over dimensions
  REAL(wp)                :: z_sum                   ! sum

!-----------------------------------------------------------------------

! forward

  DO ji = 1, k_dim
    ip      = k_pivot(ji)
    z_sum   = p_b(ip)
    p_b(ip) = p_b(ji)
    DO jj = 1, ji-1
      z_sum = z_sum - p_a(ji,jj) * p_x(jj)
    ENDDO
    p_x(ji) = z_sum
  ENDDO

! backward

  DO ji = k_dim, 1, -1
    z_sum = p_x(ji)
    DO jj = ji + 1, k_dim
      z_sum = z_sum - p_a(ji,jj) * p_x(jj)
    ENDDO
      p_x(ji) = z_sum / p_a(ji,ji)
  ENDDO

  END SUBROUTINE solve_lu

!-------------------------------------------------------------------------
!
!

  !>
  !!  Calculates coordinates of a grid point (longin, latin) in relation to.
  !!
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
!

!--------------------------------------------------------------------

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
!
!

  !>
  !! Compute the zonal and meridional components of a vector with respect.
  !!
  !! Compute the zonal and meridional components of a vector with respect
  !! to a rotated coordinate system. This function uses the same
  !! conventions and formalism as "disp_new". The algortithm is as
  !! follwos:
  !! 1) convert the vector into Cartesian coordinates
  !! 2) apply the rotation to the vector in cartesian components
  !! 3) transform the rotated vector in geographical coordinates, using
  !! the new values of latitude and longitude.
  !!
  !! @par Revision History
  !! Developed by Marco Restelli (2008-03-04)
  !!
  PURE SUBROUTINE disp_new_vect(p_gu,p_gv,lon,lat,barlon,barlat, &
                                new_lon,new_lat,new_p_gu,new_p_gv)
!
  REAL(wp), INTENT(in) :: &
    p_gu, p_gv,       & ! original lon-lat components
    lon, lat,         & ! original lon and lat
    barlon, barlat,   & ! new origin
    new_lon, new_lat    ! point (lon,lat) in the new system (this
                        ! value should be computed with disp_new)
!
  REAL(wp), INTENT(out) :: &
    new_p_gu, new_p_gv  ! new geographical components
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
                 new_lon,new_lat,new_p_gu,new_p_gv)

  END SUBROUTINE disp_new_vect

!-------------------------------------------------------------------------
!
!
!
  !>
  !! Converts cartesian coordinates to geographical.
  !!
  !!
  !! @par Revision History
  !! Developed  by Luis Kornblueh  (2004).
  !! Completely new version by Thomas Heinze (2006-07-20)
  !!
  ELEMENTAL FUNCTION cc2gc(p_x) RESULT (p_pos)
!

  TYPE(t_cartesian_coordinates), INTENT(IN) :: p_x          ! Cart. coordinates

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

  END FUNCTION  cc2gc

!-------------------------------------------------------------------------
!
!

  !>
  !! Converts geographical to cartesian coordinates.
  !!
  !!
  !! @par Revision History
  !! Developed  by Luis Kornblueh  (2004).
  !!
  ELEMENTAL FUNCTION gc2cc (p_pos)  RESULT(p_x)

  TYPE(t_geographical_coordinates), INTENT(IN) :: p_pos     ! geo. coordinates

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
!
!
!
  !>
  !!  Calculates dot product of to cartesian coordinates.
  !!
  !!
  !! @par Revision History
  !! Developed by Thomas Heinze (2006-07-05).
  !!
  ELEMENTAL FUNCTION cc_dot_product(cc_x1, cc_x2)  RESULT (p_prod)
!
  TYPE(t_cartesian_coordinates), INTENT(IN) :: cc_x1, cc_x2 ! cart. coordinates

  REAL(wp)                                :: p_prod    ! scalar product of cart.
                                                       ! coordinates

!-----------------------------------------------------------------------

  p_prod = DOT_PRODUCT (cc_x1%x, cc_x2%x)

  END FUNCTION  cc_dot_product

!-------------------------------------------------------------------------
!
!
!
  !>
  !!  Calculates 2 norm for vectors in cartesian coordinates.
  !!
  !!
  !! @par Revision History
  !! Developed by Marco Restelli (2007-11-22)
  !!
  ELEMENTAL FUNCTION cc_norm(cc_x1)  RESULT (norm)
!
  TYPE(t_cartesian_coordinates), INTENT(IN) :: cc_x1 ! cart. coordinates

  REAL(wp)                                :: norm    ! scalar product of cart.
                                                     ! coordinates

!-----------------------------------------------------------------------

  norm = SQRT(cc_dot_product(cc_x1,cc_x1))

  END FUNCTION  cc_norm

!-------------------------------------------------------------------------
!
!

  !>
  !! Computes vector product of x0,x1.
  !!
  !!
  !! @par Revision History
  !! Developed  by Luis Kornblueh  (2004).
  !!
  ELEMENTAL FUNCTION vector_product (x0, x1) RESULT(x2)

    TYPE(t_cartesian_coordinates), INTENT(IN) :: x0, x1

    TYPE(t_cartesian_coordinates) :: x2

!-----------------------------------------------------------------------

    x2%x(1) = x0%x(2)*x1%x(3) - x0%x(3)*x1%x(2)
    x2%x(2) = x0%x(3)*x1%x(1) - x0%x(1)*x1%x(3)
    x2%x(3) = x0%x(1)*x1%x(2) - x0%x(2)*x1%x(1)

  END FUNCTION vector_product

!-------------------------------------------------------------------------

  PURE FUNCTION cartesian_coordinates_plus(x,y) RESULT(z)

  TYPE(t_cartesian_coordinates) :: z
  TYPE(t_cartesian_coordinates), INTENT(in) :: x,y

  z%x = x%x + y%x

  END FUNCTION cartesian_coordinates_plus

!-----------------------------------------------------------------------

  PURE FUNCTION cartesian_coordinates_prod(a,x) RESULT(y)

  TYPE(t_cartesian_coordinates) :: y
  REAL(wp), INTENT(in) :: a
  TYPE(t_cartesian_coordinates), INTENT(in) :: x

  y%x = a * x%x

  END FUNCTION cartesian_coordinates_prod

!-----------------------------------------------------------------------
!
!

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
  FUNCTION integral_over_triangle (weight, values) RESULT(p_integ)

    REAL(wp), INTENT(IN) :: weight(3), values(3)

    REAL(wp)             :: p_integ

!-----------------------------------------------------------------------

    p_integ = DOT_PRODUCT( weight, values)
    p_integ = p_integ / 6._wp

  END FUNCTION integral_over_triangle

!-------------------------------------------------------------------------
!
!
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
!

!
!  patch on which computation is performed
!
REAL(wp), INTENT(inout) :: lat, lon
REAL(wp), INTENT(in)    :: pollat, pollon


REAL(wp) :: rotlat, rotlon

!-----------------------------------------------------------------------

rotlat = ASIN(SIN(lat)*SIN(pollat) + COS(lat)*COS(pollat)*COS(lon-pollon))
rotlon = ATAN2( COS(lat)*SIN(lon-pollon) , &
               (COS(lat)*SIN(pollat)*COS(lon-pollon)- SIN(lat)*COS(pollat)) )

lat = rotlat
lon = rotlon

END SUBROUTINE rotate_latlon


!-------------------------------------------------------------------------
!> Rotates lon-lat grid
!!
!! Rotates latitude and longitude for all grid points to standard grid
!! coordinates.
!!
!! @par Revision History
!!  developed by F. Prill, 2011-08-04
!!
SUBROUTINE rotate_latlon_grid( lon_lat_grid, rotated_pts )

  TYPE (t_lon_lat_grid), INTENT(IN)    :: lon_lat_grid
  REAL(wp),              INTENT(INOUT) :: rotated_pts(:,:,:)
  
  ! Local parameters
  REAL(wp) :: sincos_pole(2,2) ! (lon/lat, sin/cos)
  REAL(wp) :: sincos_lon(lon_lat_grid%dimen(1),2), &
    &         sincos_lat(lon_lat_grid%dimen(2),2)
  INTEGER  :: k
  REAL(wp) :: rlon_lat

!-----------------------------------------------------------------------

  sincos_pole(:,1) = SIN(lon_lat_grid%poleN(:))
  sincos_pole(:,2) = COS(lon_lat_grid%poleN(:))

  DO k=1,lon_lat_grid%dimen(1)
    rlon_lat        = lon_lat_grid%start_corner(1) + (k-1)*lon_lat_grid%delta(1)
    sincos_lon(k,:) = (/ SIN(rlon_lat), COS(rlon_lat) /)
  END DO
  DO k=1,lon_lat_grid%dimen(2)
    rlon_lat        = lon_lat_grid%start_corner(2) + (k-1)*lon_lat_grid%delta(2)
    sincos_lat(k,:) = (/ SIN(rlon_lat), COS(rlon_lat) /)
  END DO

  FORALL (k=1:lon_lat_grid%dimen(1))

    ! ASIN( SIN(phi)*SIN(poleY) + COS(phi)*COS(lambda)*COS(poleY) )
    rotated_pts(k,:,2) = &
      ASIN( sincos_lat(:,1)*sincos_pole(2,1) + sincos_lat(:,2)*sincos_lon(k,2)*sincos_pole(2,2) )

    ! ATAN2(COS(phi)*SIN(lambda), SIN(poleY)*COS(phi)*COS(lambda) - SIN(phi)*COS(poleY)) + poleX
    rotated_pts(k,:,1) = &
      ATAN2( sincos_lat(:,2)*sincos_lon(k,1), &
      &      sincos_pole(2,1)*sincos_lat(:,2)*sincos_lon(k,2) - sincos_lat(:,1)*sincos_pole(2,2)) &
      &      + lon_lat_grid%poleN(1)

  END FORALL

END SUBROUTINE rotate_latlon_grid

!-----------------------------------------------------------------------
!
!
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
!
!
!>
!! gnomonic projection for unit sphere.
!!
!! Projects a point
!! (lat,lon) onto a tangent plane with the origin at (lat_c,lon_c).
!! The results are the local cartesian coordinates (x,y).
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

  cosc = sin(lat_c)*sin(lat) + cos(lat_c)*cos(lat)*cos(lon-lon_c)
  zk   = 1._wp/cosc

  x    = zk * cos(lat)*sin(lon-lon_c)
  y    = zk * ( cos(lat_c)*sin(lat) - sin(lat_c)*cos(lat)*cos(lon-lon_c) )

END SUBROUTINE gnomonic_proj


!-------------------------------------------------------------------------
!
!
!>
!! azimuthal equidistant projection for unit sphere.
!!
!! Projects a point
!! (lat,lon) onto a tangent plane with the origin at (lat_c,lon_c).
!! The results are the local cartesian coordinates (x,y).
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

  c  = acos(sin(lat_c)*sin(lat) + cos(lat_c)*cos(lat)*cos(lon-lon_c))
  zk = c/sin(c)

  x  = zk * cos(lat)*sin(lon-lon_c)
  y  = zk * ( cos(lat_c)*sin(lat) - sin(lat_c)*cos(lat)*cos(lon-lon_c) )

END SUBROUTINE az_eqdist_proj


!-------------------------------------------------------------------------
!
!
!>
!! orthographic projection for unit sphere.
!!
!! Projects a point
!! (lat,lon) onto a tangent plane with the origin at (lat_c,lon_c).
!! The results are the local cartesian coordinates (x,y).
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

  x = cos(lat)*sin(lon-lon_c)
  y = cos(lat_c)*sin(lat) - sin(lat_c)*cos(lat)*cos(lon-lon_c)

END SUBROUTINE orthogr_proj


!!-----------------------------------------------------------------------
!!>
!! The following functions are needed for the ECHAM physics
!!
!-----------------------------------------------------------------------
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
                          x      ! integration limit

  INTEGER :: maxit = 100, m, m2
  REAL(wp) :: zeps = 3.e-7_wp, fpmin = 1.e-30_wp, &
              aa, c, d, del, h, qab, qam, qap

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
       76.18009172947146_wp, -86.50532032941677_wp, &
       24.01409824083091_wp, -1.231739572450155_wp, &
       0.1208650973866179e-2_wp, -0.5395239384953e-5_wp /)
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
!<
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
                          x      ! integration limit
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
FUNCTION gamma_fct(x)
!!------------------------------------------------------------------------------
!!
!! Description:
!!  Gamma-function from Numerical Recipes (F77)
!! Method:
!!
!!------------------------------------------------------------------------------

USE mo_kind, ONLY: wp

  IMPLICIT NONE

REAL   (wp):: gamma_fct

REAL   (wp):: cof(6) = (/76.18009173_wp, -86.50532033_wp, &
                                   24.01409822_wp, -1.231739516_wp, &
                                   0.120858003E-2_wp, -0.536382E-5_wp/)
REAL   (wp):: stp=2.50662827465_wp,                           &
                        x, xx, tmp, ser, gamma

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


END MODULE mo_math_utilities

