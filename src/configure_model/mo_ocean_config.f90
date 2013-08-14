!>
!!        Contains the variables to set up the ocean configuration
!!
!! @par Revision History
!!  Leonidas Linardakis, MPI-M, 2011/7/7
!!
!! @par Copyright
!! 2002-2006 by DWD and MPI-M
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
MODULE mo_ocean_config
!-------------------------------------------------------------------------
  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: message_text, finish

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: version = '$Id$'

  PRIVATE

  PUBLIC :: ignore_land_points
  
  !------------------------------------------------------------------------
  LOGICAL  :: ignore_land_points = .false.

CONTAINS

  !-------------------------------------------------------------------------
  SUBROUTINE init_ocean_configuration

  END SUBROUTINE init_ocean_configuration
  !-------------------------------------------------------------------------


END MODULE mo_ocean_config
