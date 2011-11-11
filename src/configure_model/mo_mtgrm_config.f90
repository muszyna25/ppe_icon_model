!>
!! Setup of meteogram output.
!!
!! @par Revision History
!! initial implementation by F. Prill, DWD (2011-09-29)
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
MODULE mo_meteogram_config

  USE mo_impl_constants, ONLY: max_dom
  USE mo_math_utilities, ONLY: t_geographical_coordinates
  USE mo_exception,      ONLY: message

  IMPLICIT NONE

  INTEGER, PARAMETER :: MAX_NAME_LENGTH      =   32  !<  max. name string length   
  INTEGER, PARAMETER :: MAX_NUM_STATIONS     =   40  !<  max. number of meteogram locations (global)
  INTEGER, PARAMETER :: FTYPE_NETCDF         =    1


  !--------------------------------------------------------------------------
  ! Derived type 
  !--------------------------------------------------------------------------

  !>
  !! Input data specifying a meteogram station.
  !!
  TYPE t_station_list
    CHARACTER(len=MAX_NAME_LENGTH)   :: zname    !< station name
    TYPE(t_geographical_coordinates) :: location !< geographical position
  END TYPE t_station_list

  !>
  !! TODO[FP] : Make some fixed-size arrays ALLOCATABLE?
  !!
  TYPE t_meteogram_output_config
    LOGICAL                        :: lenabled     !< Flag. True for output.
    CHARACTER(len=MAX_NAME_LENGTH) :: zprefix      !< file prefix string    
    INTEGER                        :: ftype        !< file type (NetCDF, ...)
    LOGICAL                        :: ldistributed !< Flag. Separate files for each PE

    ! Initial step and step interval for meteogram output
    INTEGER                        :: n0_mtgrm, ninc_mtgrm

    ! Specification of meteogram stations.
    ! Note: This info is the same for all patches.
    INTEGER                           :: nstations
    TYPE(t_station_list), POINTER     :: station_list(:,:) !< (idx, block)
    INTEGER                           :: nblks, npromz
  END TYPE t_meteogram_output_config


  ! parameters specifying meteogram output
  TYPE(t_meteogram_output_config), TARGET   :: meteogram_output_config(1:max_dom)

  PUBLIC :: t_meteogram_output_config, t_station_list
  PUBLIC :: meteogram_output_config
  PUBLIC :: FTYPE_NETCDF, MAX_NAME_LENGTH, MAX_NUM_STATIONS
  PUBLIC :: check_meteogram_configuration

  !-------------------------------------------------------------------------

CONTAINS

  SUBROUTINE check_meteogram_configuration(num_io_procs)
    INTEGER, INTENT(IN) :: num_io_procs
    ! local variables
    CHARACTER(*), PARAMETER :: routine = TRIM("mo_mtgrm_config:check_meteogram_configuration")
    INTEGER :: idom

    ! Asynchronous output does not work with distributed meteogram
    ! file output!
    IF (num_io_procs > 0) THEN
      DO idom=1,max_dom
        IF ( meteogram_output_config(idom)%lenabled     .AND. &
          &  meteogram_output_config(idom)%ldistributed ) THEN
          CALL message(routine, 'Flag "ldistributed" collides with async IO. Resetting.')
          meteogram_output_config(idom)%ldistributed = .FALSE.
        END IF
      END DO
    END IF

  END SUBROUTINE check_meteogram_configuration

END MODULE mo_meteogram_config
