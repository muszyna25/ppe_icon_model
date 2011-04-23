!>
!! @brief Contains surface type (tile) indices used by the 
!! turbulent mixing parameterization.
!!
!! @remarks
!! In the turbulent mixing parameterization, some of the prognostic
!! variables (currently u, v, T, q, TKE) are subject to turbulent tranfer
!! between the Earth's surface and the lowest model level.
!! Boundary condition of this transfer is formulated separately
!! for different surface types. This module contains parameters
!! specifying how many different types are considered in a particular
!! simulation (variable "nsfc_type"), and their corresponding indices
!! (e.g., iwtr for water, ilnd for land, etc.). In addition, the
!! aggregated grid-box mean value is given the index "igbm".
!!
!! @author Hui Wan, MPI-M
!!
!! @par Revision History
!! First version by Hui Wan, MPI-M (2010-09-21)
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
MODULE mo_icoham_sfc_indices

  USE mo_exception, ONLY: message, message_text

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: nsfc_type, iwtr, iice, ilnd, igbm   !< variables
  PUBLIC :: init_sfc_indices                    !< subroutine

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  INTEGER :: nsfc_type   !< total number of surface types
  INTEGER :: iwtr = 1    !< index for water-covered surface
  INTEGER :: iice = 2    !< index for ice-covered   surface
  INTEGER :: ilnd = 3    !< index for land          surface
  INTEGER :: igbm        !< index for grid-box mean

CONTAINS
  !>
  !! Set surface indices according to the simulation setup
  !! (e.g., dynamical core test, aqua-planet, or real
  !! climate simulation).
  !!
  SUBROUTINE init_sfc_indices( ltestcase, ctest_name )

    LOGICAL,         INTENT(IN) :: ltestcase
    CHARACTER(len=*),INTENT(IN) :: ctest_name

    IF (ltestcase) THEN
      SELECT CASE(TRIM(ctest_name))
      CASE('APE')
      ! Aqua-planet simulation, no land, no ice;
      ! No needed to distinguish the aggregated grid-box mean
      ! and the value on different types of surface

        igbm      = 1
        iwtr      = 1
        nsfc_type = 1
        iice      = 999
        ilnd      = 999

      CASE('JWw-Moist','LDF-Moist')
      ! Baroclinic wave test, no land, no ice.

        igbm      = 1
        iwtr      = 1
        nsfc_type = 1
        iice      = 999
        ilnd      = 999

      ! maybe worth trying later:
      ! igbm      = 0
      ! iwtr      = 1
      ! iice      = 2
      ! nsfc_type = 2
      ! ilnd      = 999

      ! A wild idea: a dry or completely frozen planet
      ! igbm      = 0
      ! iice      = 1
      ! ilnd      = 2
      ! nsfc_type = 2
      ! iwtr      = 999

      END SELECT
    ELSE
    ! Standard setup for real climate simulation.
    ! Three surface types are considered.

      igbm      = 0
      iwtr      = 1
      iice      = 2
      ilnd      = 3
      nsfc_type = 3
    ENDIF

    WRITE(message_text,*) " "
    CALL message("mo_icoham_sfc_indices/init_sfc_indices",TRIM(message_text))

    WRITE(message_text,'(i3,a)')    &
      & nsfc_type, " surface type(s) activated."
    CALL message("mo_icoham_sfc_indices/init_sfc_indices",TRIM(message_text))

    WRITE(message_text,'(a,4i4,a)') &
      & "Indices for grid-box mean, water, ice, and land are ", &
      & igbm, iwtr, iice, ilnd, ", respectively."
    CALL message("mo_icoham_sfc_indices/init_sfc_indices",TRIM(message_text))

    WRITE(message_text,*) " "
    CALL message("mo_icoham_sfc_indices/init_sfc_indices",TRIM(message_text))

  END SUBROUTINE init_sfc_indices
  !-------------

END MODULE mo_icoham_sfc_indices
