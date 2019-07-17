!>
!! @brief Contains surface type (tile) indices used by the 
!! turbulent mixing parameterization.
!!
!! @remarks
!! In the turbulent mixing parameterization some of the prognostic
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
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_echam_sfc_indices

  USE mo_exception, ONLY: message, message_text

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: nsfc_type, iwtr, iice, ilnd, igbm   !< index variables
  PUBLIC :: csfc                                !< sfc names
  PUBLIC :: init_sfc_indices                    !< subroutine

  INTEGER :: nsfc_type   !< total number of surface types
  INTEGER :: iwtr = 1    !< index for water-covered surface
  INTEGER :: iice = 2    !< index for ice-covered   surface
  INTEGER :: ilnd = 3    !< index for land          surface
  INTEGER :: igbm        !< index for grid-box mean

  CHARACTER(LEN=3) :: csfc(3) = (/'wtr','ice','lnd'/)

CONTAINS
  !>
  !! Set surface indices according to the simulation setup
  !! (e.g., dynamical core test, aqua-planet, or real
  !! climate simulation).
  !!
  SUBROUTINE init_sfc_indices( ctest_name )

    CHARACTER(len=*),INTENT(IN) :: ctest_name

    SELECT CASE(TRIM(ctest_name))
    CASE('APE','APE_echam','RCE','RCE_glb','RCE_Tconst','RCEhydro')
      ! Aqua-planet simulation, no land, no ice;
      ! No needed to distinguish the aggregated grid-box mean
      ! and the value on different types of surface

      iwtr      = 1
      nsfc_type = 1
      igbm      = 0
      iice      = 999
      ilnd      = 999

    CASE('APEi','APEc','APEc_nh')
      ! Aqua-planet simulation with ice, but no land;

      iwtr      = 1
      iice      = 2
      nsfc_type = 2
      igbm      = 0
      ilnd      = 999

    CASE('JWw-Moist','LDF-Moist','jabw_m')
      ! Baroclinic wave test, no land, no ice.

      iwtr      = 1
      nsfc_type = 1
      igbm      = 0
      iice      = 999
      ilnd      = 999

      ! maybe worth trying later:
      ! iwtr      = 1
      ! iice      = 2
      ! nsfc_type = 2
      ! igbm      = 0
      ! ilnd      = 999

      ! A wild idea: a dry or completely frozen planet
      ! iice      = 1
      ! ilnd      = 2
      ! nsfc_type = 2
      ! igbm      = 0
      ! iwtr      = 999

    CASE('TPEo','TPEc')
      ! Terra-planet simulation
      ! no ocean, no ice but lakes and ice on lakes ... therefore have to use iice and iwtr
      ! No need to distinguish the aggregated grid-box mean
      ! and the value on different types of surface

      nsfc_type = 3
      iwtr      = 1
      iice      = 2
      ilnd      = 3
      igbm      = 0

    CASE DEFAULT
      ! Standard setup for real-world climate simulation.
      ! Three surface types are considered.

      iwtr      = 1
      iice      = 2
      ilnd      = 3
      nsfc_type = 3
      igbm      = 0

    END SELECT


    WRITE(message_text,*) " "
    CALL message("mo_echam_sfc_indices/init_sfc_indices",TRIM(message_text))

    WRITE(message_text,'(i3,a)')    &
      & nsfc_type, " surface type(s) activated."
    CALL message("mo_echam_sfc_indices/init_sfc_indices",TRIM(message_text))

    WRITE(message_text,'(a,4i4,a)') &
      & "Indices for water, ice, land, and grid-box mean are ", &
      & iwtr, iice, ilnd, igbm, ", respectively."
    CALL message("mo_echam_sfc_indices/init_sfc_indices",TRIM(message_text))

    WRITE(message_text,*) " "
    CALL message("mo_echam_sfc_indices/init_sfc_indices",TRIM(message_text))

  END SUBROUTINE init_sfc_indices
  !-------------

END MODULE mo_echam_sfc_indices
