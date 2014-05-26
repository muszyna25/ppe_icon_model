!>
!! @brief configuration setup for JSBACH land scheme
!!
!! configuration setup for JSBACH land scheme
!! <Details of procedures are documented below with their definitions.>
!! <Include any applicable external references inline as module::procedure,>
!! <external_procedure(), or by using @see.>
!! <Don't forget references to literature.>
!!
!! @author Reiner Schnur, MPI-M (2012-04-16)
!! @author <name, affiliation>
!!
!!
!! @par Revision History
!! <Description of activity> by <name, affiliation> (<YYYY-MM-DD>)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_lnd_jsbach_config

  USE mo_grid_config,    ONLY: n_dom
  USE mo_kind,           ONLY: wp
  USE mo_exception,      ONLY: finish

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: t_lnd_jsbach_config, lnd_jsbach_config  !< derived type and variable for configuration
  PUBLIC :: configure_lnd_jsbach                    !< subroutine

  !>
  !! Derived type containing main switches for JSBACH land surface scheme
  !!
  TYPE t_lnd_jsbach_config                       !< Derived type for configuration state

    INTEGER :: ntiles           !< Number of tiles
    INTEGER :: nsoil            !< Number of soil layers
    INTEGER :: ntsoil           !< Number of soil layers for soil temperature
    REAL(wp), POINTER :: zlev_soil (:)  !< Soil layers [m]
    REAL(wp), POINTER :: ztlev_soil(:)  !< Soil layers for soil temperature [m]

  END TYPE t_lnd_jsbach_config

  TYPE(t_lnd_jsbach_config), ALLOCATABLE :: lnd_jsbach_config(:) !< Configuration state variable

CONTAINS
  !>
  !!
  SUBROUTINE configure_lnd_jsbach

!!$    CHARACTER(LEN=*),PARAMETER  :: &
!!$             & routine ='mo_lnd_jsbach_config:configure_lnd_jsbach'

    ALLOCATE(lnd_jsbach_config(n_dom))

  END SUBROUTINE configure_lnd_jsbach

END MODULE mo_lnd_jsbach_config
