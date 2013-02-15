!>
!! Type definition and utility routines of lon-lat grids used in the model.
!!
!! @author F. Prill, DWD
!!
!! @par Revision History
!! Initial implementation  by  F. Prill, DWD (2012-03-20)
!!
!! @par Copyright
!! 2002-2010 by DWD and MPI-M
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
!! -----------------------------------------------------------------------------------
MODULE mo_lonlat_grid

  USE mo_kind,                ONLY: wp
  USE mo_physical_constants,  ONLY: earth_radius

  IMPLICIT NONE

  ! subroutines + functions:
  PUBLIC :: rotate_latlon_grid
  PUBLIC :: latlon_compute_area_weights
  PUBLIC :: compute_lonlat_specs
  PUBLIC :: compute_lonlat_blocking
  ! types and variables:
  PUBLIC :: t_lon_lat_grid


  !> specification of (rotated) lon-lat grid
  !
  !  @note Variables of this type contain redundant information. The
  !        user must take care if all entries have been consistently
  !        initialized.
  !
  TYPE t_lon_lat_grid

    ! description of regular grid (in degrees) as in namelist
    REAL(wp) :: reg_lon_def(3)            ! start, increment, end longitude in degrees
    REAL(wp) :: reg_lat_def(3)            ! start, increment, end latitude in degrees
    REAL(wp) :: north_pole (2)            ! position of north pole (lon,lat) in degrees
                                           
    INTEGER  :: lon_dim                   ! Number of points in lon direction
    INTEGER  :: lat_dim                   ! Number of points in lat direction

    ! Alternative description of regular grid in rad,
    !   start_corner(lon,lat) + delta(lon,lat)*[0,1,2, ..., dim(lon,lat)-1]
    REAL(wp) ::                     &
      &  delta       (2),           &     ! lon-lat grid resolution,                unit:rad
      &  start_corner(2)                  ! south western corner of area (lon/lat), unit:rad

    ! Blocking information, computed from above values:
    INTEGER  :: total_dim                 ! total number of grid points 
    INTEGER  :: nblks, npromz             ! blocking info
    
  END TYPE t_lon_lat_grid

  
CONTAINS

  !---------------------------------------------------------------
  !> Computes rad-unit entries in lon-lat-grid definition
  !
  SUBROUTINE compute_lonlat_specs(lonlat_grid)
    TYPE (t_lon_lat_grid), INTENT(INOUT) :: lonlat_grid
    REAL(wp) :: pi_180

    pi_180 = ATAN(1._wp)/45._wp
    lonlat_grid%delta(1)        = lonlat_grid%reg_lon_def(2) * pi_180
    lonlat_grid%delta(2)        = lonlat_grid%reg_lat_def(2) * pi_180
    lonlat_grid%start_corner(1) = lonlat_grid%reg_lon_def(1) * pi_180
    lonlat_grid%start_corner(2) = lonlat_grid%reg_lat_def(1) * pi_180
  END SUBROUTINE compute_lonlat_specs


  !---------------------------------------------------------------
  !> Setup of nproma-blocking for lon-lat-grid
  !
  SUBROUTINE compute_lonlat_blocking(lonlat_grid, inproma)
    TYPE (t_lon_lat_grid), INTENT(INOUT) :: lonlat_grid
    INTEGER,               INTENT(IN)    :: inproma

    lonlat_grid%total_dim = lonlat_grid%lon_dim * lonlat_grid%lat_dim
    lonlat_grid%nblks     = (lonlat_grid%total_dim - 1)/inproma + 1
    lonlat_grid%npromz    = lonlat_grid%total_dim - (lonlat_grid%nblks-1)*inproma
  END SUBROUTINE compute_lonlat_blocking


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
    REAL(wp), PARAMETER :: ZERO_TOL = 1.e-15

    REAL(wp) :: sincos_pole(2,2),                   &  ! (lon/lat, sin/cos)
      &         sincos_lon(lon_lat_grid%lon_dim,2), &
      &         sincos_lat(lon_lat_grid%lat_dim,2), &
      &         pi_180, npole_rad(2), arg1, arg2,   &
      &         rlon(lon_lat_grid%lon_dim),         &
      &         rlat(lon_lat_grid%lat_dim)
    INTEGER  :: k, j
    LOGICAL  :: ltrivial_rotation
    
    !-----------------------------------------------------------------------
    
    ltrivial_rotation = ((ABS(90._wp - lon_lat_grid%north_pole(2)) < ZERO_TOL) .AND.  &
      &                   ABS( 0._wp - lon_lat_grid%north_pole(1)) < ZERO_TOL)
    DO k=1,lon_lat_grid%lon_dim
      rlon(k)         = lon_lat_grid%start_corner(1) + REAL(k-1,wp)*lon_lat_grid%delta(1)
      sincos_lon(k,:) = (/ SIN(rlon(k)), COS(rlon(k)) /)
    END DO
    DO k=1,lon_lat_grid%lat_dim
      rlat(k)         = lon_lat_grid%start_corner(2) + REAL(k-1,wp)*lon_lat_grid%delta(2)
      sincos_lat(k,:) = (/ SIN(rlat(k)), COS(rlat(k)) /)
    END DO

    ! special treatment for "trivial rotation" (no.pole at +90,0):
    IF (ltrivial_rotation) THEN

      DO j = 1, lon_lat_grid%lat_dim
        DO k = 1, lon_lat_grid%lon_dim
          rotated_pts(k,j,1) = rlon(k)
          rotated_pts(k,j,2) = rlat(j)
        END DO
      END DO

    ELSE

      ! convert north pole: degree -> rad:
      pi_180 = ATAN(1._wp)/45._wp
      npole_rad(:) = lon_lat_grid%north_pole(:)*pi_180

      sincos_pole(:,1) = SIN(npole_rad(:))
      sincos_pole(:,2) = COS(npole_rad(:))  

      DO j = 1, lon_lat_grid%lat_dim
        DO k = 1, lon_lat_grid%lon_dim

          ! ATAN2(COS(phi)*SIN(lambda), SIN(poleY)*COS(phi)*COS(lambda) - SIN(phi)*COS(poleY)) + poleX
          arg1 = sincos_lat(j,2)*sincos_lon(k,1)
          arg2 = sincos_pole(2,1)*sincos_lat(j,2)*sincos_lon(k,2) - sincos_lat(j,1)*sincos_pole(2,2)

          IF ((ABS(arg1) > ZERO_TOL) .OR. (ABS(arg2) > ZERO_TOL)) THEN
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
    END IF
  END SUBROUTINE rotate_latlon_grid


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
    
    radius = earth_radius ! earth's radius (average)
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

END MODULE mo_lonlat_grid
