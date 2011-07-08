!>
!! <Short description of module for listings and indices>
!!
!! <Describe the concepts of the procedures and algorithms used in the module.>
!! <Details of procedures are documented below with their definitions.>
!! <Include any applicable external references inline as module::procedure,>
!! <external_procedure(), or by using @see.>
!! <Don't forget references to literature.>
!!
!! @author <name, affiliation>
!! @author <name, affiliation>
!!
!!
!! @par Revision History
!! <Description of activity> by <name, affiliation> (<YYYY-MM-DD>)
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
MODULE mo_ha_dyn_config

  USE mo_kind, ONLY: wp

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_ha_dyn_config, ha_dyn_config

  CHARACTER(len=*),PARAMETER,PRIVATE :: version = '$Id$'

  !>
  !!--------------------------------------------------------------------------
  !! Derived type containing control variables specific to the hydrostatic 
  !! atm dynamical core
  !!--------------------------------------------------------------------------
  TYPE :: t_ha_dyn_config

    INTEGER  :: ileapfrog_startup  !<
    REAL(wp) :: asselin_coeff      !<

    INTEGER  :: si_expl_scheme     !< Scheme for the explicit part of the
                                   !< 2-time-level semi-implicit time integration.
                                   !< See mo_atm_constants for the options.
    REAL(wp) :: si_2tls            !< = 0 : explicit scheme
    REAL(wp) :: si_coeff           !< = 0 : explicit scheme
                                   !< = 1 : semi implicit scheme.
                                   !< in (0,1): a weighted scheme
    REAL(wp) :: si_offctr          !< Weighting parameter used in calculating the
                                   !< second temporal derivatives in the semi-implicit
                                   !< correction scheme. The value read from namelist are
                                   !< assumed to be the offcentering (i.e. between 0 and 1).
    REAL(wp) :: si_cmin            !< Min. phase speed of the decomposed modes to be
                                   !< solved by the semi-implicit correction scheme
    REAL(wp) :: si_rtol            !< Relative tolerance
    LOGICAL  :: lsi_3d             !< If .true., solve the 3D equation

    LOGICAL :: ldry_dycore !< If .TRUE., ignore the effact of water vapor,
                           !< cloud liquid and cloud ice on virtual temperature.
    LOGICAL :: lref_temp   !< If .TRUE., involve the reference temperature profile
                           !< in the calculation of pressure gradient force.

  END TYPE t_ha_dyn_config 

  !>
  !!
  TYPE(t_ha_dyn_config) :: ha_dyn_config

END MODULE mo_ha_dyn_config

