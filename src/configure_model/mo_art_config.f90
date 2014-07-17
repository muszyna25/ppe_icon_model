!>
!! @brief configuration setup for ART-package
!!
!! configuration setup for ART-package
!! <Details of procedures are documented below with their definitions.>
!! <Include any applicable external references inline as module::procedure,>
!! <external_procedure(), or by using @see.>
!! <Don't forget references to literature.>
!!
!! @author Daniel Reinert, DWD
!!
!!
!! @par Revision History
!! Initial revision by Daniel Reinert, DWD (2011-12-08)
!! Modifications by Kristina Lundgren, KIT (2012-07-03)
!! Modifications by Daniel Rieger,     KIT (2014-17-06)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_art_config
  USE mo_kind,                 ONLY: wp
  USE mo_impl_constants,       ONLY: max_dom
  USE mo_math_utilities,       ONLY: t_geographical_coordinates  
  IMPLICIT NONE


  PUBLIC 


  !!--------------------------------------------------------------------------
  !! Basic configuration setup for ICON-ART
  !!--------------------------------------------------------------------------
   INTEGER, PARAMETER  :: max_volc_input  = 20 !Maximum number of volcanoes in input namelist art_volclist_tot
   INTEGER             :: nart_tendphy  = 0    !Maximum number of tracers that are effected by deep convective transport 

  TYPE t_volc_list
    CHARACTER(len=20)                :: zname    ! < name of volcanoe or location
    TYPE(t_geographical_coordinates) :: location !< geographical position
  END TYPE t_volc_list

  TYPE t_art_config ! Namelist variables for ART

    ! Namelist variables

    ! General variables (Details: cf. Tab. 2.1 ICON-ART User Guide)
    CHARACTER(LEN=120) :: cart_folder  !< Absolute Path to ART source code
    INTEGER :: iart_ntracer            !< number of transported ART tracers
    
    ! Atmospheric Chemistry (Details: cf. Tab. 2.2 ICON-ART User Guide)
    LOGICAL :: lart_chem               !< Main switch to enable chemistry
    INTEGER :: iart_chem_mechanism     !< Selects the chemical mechanism
    
    ! Atmospheric Aerosol (Details: cf. Tab. 2.3 ICON-ART User Guide)
    LOGICAL :: lart_aerosol            !< Main switch for the treatment of atmospheric aerosol
    INTEGER :: iart_seasalt            !< Treatment of sea salt aerosol
    INTEGER :: iart_dust               !< Treatment of mineral dust aerosol
    INTEGER :: iart_anthro             !< Treatment of anthropogenic aerosol
    INTEGER :: iart_fire               !< Treatment of wildfire aerosol
    INTEGER :: iart_volcano            !< Treatment of volcanic ash aerosol
    CHARACTER(LEN=120) :: cart_volcano_file  !< Absolute path + filename of input file for volcanoes
    INTEGER :: iart_radioact           !< Treatment of radioactive particles
    CHARACTER(LEN=120) :: cart_radioact_file !< Absolute path + filename of input file for radioactive emissions
    INTEGER :: iart_pollen             !< Treatment of pollen
    
    ! Feedback processes (Details: cf. Tab. 2.4 ICON-ART User Guide)
    INTEGER :: iart_aci_warm           !< Nucleation of aerosol to cloud droplets
    INTEGER :: iart_aci_cold           !< Nucleation of aerosol to cloud ice
    INTEGER :: iart_ari                !< Direct interaction of aerosol with radiation
    
    ! Fast Physics Processes (Details: cf. Tab. 2.5 ICON-ART User Guide)
    LOGICAL :: lart_conv               !< Convection of aerosol (TRUE/FALSE)
    INTEGER :: nconv_tracer            !< number of tracers in convection 
    LOGICAL :: lart_turb               !< Turbulent diffusion of aerosol (TRUE/FALSE)
    INTEGER :: nturb_tracer            !< number of tracers in turbulence

  END TYPE t_art_config

  TYPE(t_art_config), TARGET :: art_config(0:max_dom)


CONTAINS

  !>
  !! setup components of ICON-ART depending on this namelist
  !!
  !! Setup of additional ICON-ART control variables depending on the 
  !! art-NAMELIST and potentially other namelists. This routine is 
  !! called, after all namelists have been read and a synoptic consistency 
  !! check has been done.
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2011-12-08)
  !! Modification by Kristina Lundgren, KIT (2012-11-27)
  !
  SUBROUTINE configure_art(jg)
    INTEGER, INTENT(IN) :: jg          !< patch 

    art_config(jg)%nconv_tracer=0
    art_config(jg)%nturb_tracer=0
 
  END SUBROUTINE configure_art


END MODULE mo_art_config
