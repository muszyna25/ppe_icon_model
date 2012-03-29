!>
!! Type definition and global list of lon-lat grids used in the model.
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

  USE mo_kind,  ONLY: wp

  IMPLICIT NONE

  PUBLIC :: compute_lonlat_specs
  PUBLIC :: compute_lonlat_blocking
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


END MODULE mo_lonlat_grid
