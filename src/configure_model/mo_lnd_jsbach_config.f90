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
!! @par Copyright
!! 2002-2012 by DWD and MPI-M
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
MODULE mo_lnd_jsbach_config

  USE mo_kind,           ONLY: wp
  USE mo_exception,      ONLY: finish

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: t_lnd_jsbach_config, lnd_jsbach_config  !< derived type and variable for configuration
  PUBLIC :: configure_lnd_jsbach                    !< subroutine
  PUBLIC :: get_ntiles_jsbach, get_nsoil_jsbach     !< functions
  PUBLIC :: get_ntsoil_jsbach                       !< functions

  CHARACTER(len=*), PARAMETER, PRIVATE :: version = '$Id$'

  !>
  !! Derived type containing main switches for JSBACH land surface scheme
  !!
  TYPE t_lnd_jsbach_config                       !< Derived type for configuration state

    INTEGER :: ntiles           !< Number of tiles
    INTEGER :: nsoil            !< Number of soil layers
    INTEGER :: ntsoil           !< Number of soil layers for soil temperature

  END TYPE t_lnd_jsbach_config

  TYPE(t_lnd_jsbach_config) :: lnd_jsbach_config !< Configuration state variable

CONTAINS
  !>
  !!
  SUBROUTINE configure_lnd_jsbach(ltestcase, ctest_name)

    LOGICAL,         INTENT(IN) :: ltestcase
    CHARACTER(LEN=*),INTENT(IN) :: ctest_name
    CHARACTER(LEN=*),PARAMETER  :: &
             & routine ='mo_lnd_jsbach_config:configure_lnd_jsbach'

  END SUBROUTINE configure_lnd_jsbach

  INTEGER FUNCTION get_ntiles_jsbach()
    get_ntiles_jsbach = lnd_jsbach_config%ntiles
  END FUNCTION get_ntiles_jsbach

  INTEGER FUNCTION get_nsoil_jsbach()
    get_nsoil_jsbach = lnd_jsbach_config%nsoil
  END FUNCTION get_nsoil_jsbach

  INTEGER FUNCTION get_ntsoil_jsbach()
    get_ntsoil_jsbach = lnd_jsbach_config%ntsoil
  END FUNCTION get_ntsoil_jsbach

END MODULE mo_lnd_jsbach_config
