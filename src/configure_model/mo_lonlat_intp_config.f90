!>
!! Setup of lon-lat config based on namelist parameters.
!!
!! @par Revision History
!! initial implementation by F. Prill, DWD (2011-09-15)
!!
!! @par Copyright
!! 2002-2011 by DWD and MPI-M
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
MODULE mo_lonlat_intp_config

  USE mo_impl_constants, ONLY: MAX_CHAR_LENGTH, max_dom
  USE mo_exception,      ONLY: message, message_text, finish
  USE mo_math_utilities, ONLY: t_lon_lat_grid
  USE mo_util_string,    ONLY: split_string
  USE mo_kind,           ONLY: wp

  IMPLICIT NONE

  INTEGER, PARAMETER :: DIM_UNDEFINED = -1

  !--------------------------------------------------------------------------
  ! Derived type 
  !--------------------------------------------------------------------------

  ! lon-lat grid for (optional) interpolation of output variables:
  TYPE t_lonlat_intp_config
    ! Flag. True, if interpolation onto lon-lat grid is enabled
    LOGICAL                  :: l_enabled
    ! Flag. True, if standard variable is not written for lon-lat vars.
    LOGICAL              :: l_supersede
    ! lon-lat grid specification for each patch
    ! (we don't know about "n_dom_start" or "ndom" here)
    TYPE (t_lon_lat_grid)    :: lonlat_grid

    ! string with a list of variable names due for lon-lat
    ! interpolation (or "all")
    ! Note: This list is identical on all patches.
    CHARACTER (len=MAX_CHAR_LENGTH)  :: zlist

    ! list, split into its components:
    INTEGER                  :: n_list
    INTEGER                  :: pos_list(50), ilen_list(50)
  END TYPE t_lonlat_intp_config

  ! parameters specifying interpolation to (rotated) lon-lat grids
  TYPE(t_lonlat_intp_config), TARGET :: lonlat_intp_config(1:max_dom)

  ! (optional:) second corner of lon-lat area
  REAL(wp) :: lonlat_corner2(2, 1:max_dom)

  PUBLIC :: configure_lonlat_intp
  PUBLIC :: t_lonlat_intp_config
  PUBLIC :: lonlat_intp_config
  PUBLIC :: lonlat_corner2
  PUBLIC :: DIM_UNDEFINED

CONTAINS

  !>
  !! Setup of lon-lat config based on namelist parameters.
  !!
  !! @par Revision History
  !! initial implementation by F. Prill, DWD (2011-09-15)
  SUBROUTINE configure_lonlat_intp(maxdom, nproma, num_io_procs)

    INTEGER, INTENT(IN) :: maxdom, nproma, num_io_procs
    ! local variables
    CHARACTER(len=*), PARAMETER    :: routine = &
      &  'mo_lonlat_intp_config:configure_lonlat_intp'
    TYPE (t_lon_lat_grid), POINTER :: grid ! lon-lat grid
    INTEGER :: idom

    DO idom=1,maxdom
      IF (lonlat_intp_config(idom)%l_enabled) THEN
        ! string with a list of variable names due for lon-lat
        ! interpolation
        CALL split_string(lonlat_intp_config(idom)%zlist, &
          & lonlat_intp_config(idom)%n_list,     &
          & lonlat_intp_config(idom)%pos_list,   &
          & lonlat_intp_config(idom)%ilen_list)
        
        ! if the user has specified a second corner, override the
        ! grid point numbers by the necessary values:
        grid => lonlat_intp_config(idom)%lonlat_grid
        IF (ALL(lonlat_corner2(:,idom) /= DIM_UNDEFINED)) THEN
          grid%dimen(:) = NINT((lonlat_corner2(:,idom) - grid%start_corner(:))/grid%delta(:)) + 1
        END IF

        ! values for the blocking
        grid%total_dim      = PRODUCT(grid%dimen(:))
        grid%nblks          =  (grid%total_dim - 1)/nproma + 1
        grid%npromz         =  grid%total_dim - (grid%nblks-1)*nproma

        IF (grid%total_dim > 0) THEN
          WRITE(message_text,*) "read lon-lat grid definition, total size = ", grid%total_dim
          CALL message(routine, TRIM(message_text))
        END IF
      END IF
    END DO

    ! Throw an error message for asynchronous IO mode, as this has not
    ! yet been implemented
    IF (ANY(lonlat_intp_config(:)%l_enabled) .AND. &
      & (num_io_procs > 0)) THEN
      CALL finish(routine, &
        &    "Lon-lat interpolation not yet implemented for asynchronous IO mode.")
    END IF

  END SUBROUTINE configure_lonlat_intp

END MODULE mo_lonlat_intp_config
