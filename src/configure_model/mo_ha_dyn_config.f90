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

  USE mo_kind,           ONLY: wp
  USE mo_exception,      ONLY: message, print_value
  USE mo_impl_constants, ONLY: LEAPFROG_EXPL, LEAPFROG_SI, TWO_TL_SI 

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_ha_dyn_config, ha_dyn_config, configure_ha_dyn

  CHARACTER(len=*),PARAMETER,PRIVATE :: version = '$Id$'

  !>
  !!--------------------------------------------------------------------------
  !! Derived type containing control variables specific to the hydrostatic 
  !! atm dynamical core
  !!--------------------------------------------------------------------------
  TYPE :: t_ha_dyn_config

    INTEGER  :: itime_scheme       !< Time stepping scheme. See mo_impl_constants.f90 for 
                                   !< the available choices.
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

    LOGICAL :: ltheta_dyn  !< If .TRUE., use delta-p*potential-temperature as thermodynamic
                           !< prognostic variable

  END TYPE t_ha_dyn_config 

  !>
  !!
  TYPE(t_ha_dyn_config) :: ha_dyn_config

CONTAINS
  !>
  !! Currently contains only printing
  !!
  SUBROUTINE configure_ha_dyn


    CALL message('','')
    CALL message('','------ Hydrostatic atm dynamical core ------')

    SELECT CASE(ha_dyn_config% itime_scheme)
    CASE (LEAPFROG_EXPL)

      CALL message('', '--- itime_scheme      : LEAPFROG_EXPL')
      CALL print_value('ileapfrog_startup', ha_dyn_config% ileapfrog_startup)
      CALL print_value('asseline_coeff   ', ha_dyn_config% asselin_coeff)

    CASE (LEAPFROG_SI)

      CALL message('', '--- itime_scheme      : LEAPFROG_SI')
      CALL print_value('ileapfrog_startup', ha_dyn_config% ileapfrog_startup)
      CALL print_value('asseline_coeff   ', ha_dyn_config% asselin_coeff    )
      CALL print_value('si_coeff         ', ha_dyn_config% si_coeff         )
      CALL print_value('si_offctr        ', ha_dyn_config% si_offctr        )
      CALL print_value('si_cmin          ', ha_dyn_config% si_cmin          )
      CALL print_value('si_rtol          ', ha_dyn_config% si_rtol          )
      CALL print_value('lsi_3d           ', ha_dyn_config% lsi_3d           )

    CASE (TWO_TL_SI)

      CALL message('', '--- itime_scheme      : TWO_TL_SI')
      CALL print_value('si_expl_scheme   ', ha_dyn_config% si_expl_scheme)
      CALL print_value('si_2tls          ', ha_dyn_config% si_2tls       )
      CALL print_value('si_rtol          ', ha_dyn_config% si_rtol       )

    END SELECT

    CALL print_value('ldry_dycore      ', ha_dyn_config% ldry_dycore)
    CALL print_value('lref_temp        ', ha_dyn_config% lref_temp  )
    CALL print_value('ltheta_dyn       ', ha_dyn_config% ltheta_dyn )

    CALL message('','--------------------------')
    CALL message('','')

  END SUBROUTINE configure_ha_dyn
  !-----------
END MODULE mo_ha_dyn_config

