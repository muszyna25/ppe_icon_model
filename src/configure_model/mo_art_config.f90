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
  
  PRIVATE

  PUBLIC :: nart_tendphy
  PUBLIC :: t_art_config
  PUBLIC :: art_config
  PUBLIC :: configure_art
  ! Tracer indices
  PUBLIC :: iash1,iash2,iash3,iash4,iash5,iash6                          !Running index for Volcanic Ash in ICON-ART 
  PUBLIC :: iasha, iashb, iashc, iasha0, iashb0, iashc0                  !Running index for Volcanic Ash in ICON-ART
  PUBLIC :: iCS137,iI131,iTE132,iZR95,iXE133,iI131g,iI131o,iBA140,iRU103 !Running index for radioactive nuclides  in ICON-ART
  PUBLIC :: iseasa,iseasb,iseasc,iseasa0,iseasb0,iseasc0                 !Running index for sea salt in ICON-ART
  PUBLIC :: idusta,idustb,idustc,idusta0,idustb0,idustc0                 !Running index for mineral dust in ICON-ART
  PUBLIC :: iTRCHBR3,iTRCH2BR2,iTRBRy                                    !Running index for chemical tracer in ICON-ART - VSLS-BRy
  PUBLIC :: iTRCH4,iTRCO2,iTRCO,iTRH2O,iTRO3                             !Running index for chemical tracer in ICON-ART - CH4-CO-CO2-H2O-O3
  PUBLIC :: iTRCH3COCH3,iTRC2H6,iTRSF6,iTRN2O                            !Running index for chemical tracer in ICON-ART - CH3COCH3,C2H6,SF6,N2O
  PUBLIC :: iTR1,iTR2,iTR3,iTR4,iTR5                                     !Running index for chemical tracer in ICON-ART - artificial tracer
  
  
  
  !!--------------------------------------------------------------------------
  !! Tracer indices of ICON-ART species
  !!--------------------------------------------------------------------------
  INTEGER :: & !< Volcanic ash tracer indicies (bulk scheme)
    &  iash1, iash2, iash3, iash4, iash5, iash6
  INTEGER :: & !< Volcanic ash tracer indicies (modal scheme)
    &  iasha, iashb, iashc, iasha0, iashb0, iashc0
  INTEGER :: & !< Radioactive nuclides (bulk scheme)
    &  iCS137, iI131, iTE132, iZR95, iXE133, &
    &  iI131g, iI131o, iBA140, iRU103
  INTEGER :: & !< seasalt aerosol tracer indicies (modal scheme)
    &  iseasa, iseasb, iseasc, iseasa0, iseasb0, iseasc0
  INTEGER :: & !< mineral dust aerosol tracer indicies (modal scheme)
    &  idusta, idustb, idustc, idusta0, idustb0, idustc0
  INTEGER :: & !< Chemical tracers
    &  iTRCHBR3, iTRCH2BR2, iTRBRy, &
    &  iTRCH4, iTRCO2, iTRCO,       &
    &  iTRCH3COCH3, iTRC2H6, iTRH2O,&
    &  iTRO3, iTRSF6, iTRN2O,       &
    &  iTR1, iTR2, iTR3, iTR4, iTR5
  
  !!--------------------------------------------------------------------------
  !! Basic configuration setup for ICON-ART
  !!--------------------------------------------------------------------------
  INTEGER, PARAMETER  :: npreslay      = 7  !Number of pressure layers for diagnostic output of maximum concentration

  INTEGER             :: nart_tendphy  = 0  !Maximum number of tracers that are effected by deep convective transport 
  
  
  TYPE t_art_config ! Namelist variables for ART
    
    ! Namelist variables
    
    ! General control variables (Details: cf. Tab. 2.2 ICON-ART User Guide)
    CHARACTER(LEN=120) :: cart_folder  !< Absolute Path to ART source code
    INTEGER :: iart_ntracer            !< number of transported ART tracers
    INTEGER :: iart_init_aero          !< Initialization of aerosol species
    INTEGER :: iart_init_gas           !< Initialization of gaseous species
    LOGICAL :: lart_diag_out           !< Enable output of diagnostic fields
    
    ! Atmospheric Chemistry (Details: cf. Tab. 2.3 ICON-ART User Guide)
    LOGICAL :: lart_chem               !< Main switch to enable chemistry
    INTEGER :: iart_chem_mechanism     !< Selects the chemical mechanism
    CHARACTER(LEN=120) :: cart_emiss_table_path  !< path of tex-files with meta data of emissions. ! MiW
    CHARACTER(LEN=120) :: cart_emiss_table_file  !< file names of tex-files with meta data of emissions without "_DOM??.tex" at the end. ! MiW
 
    ! Atmospheric Aerosol (Details: cf. Tab. 2.4 ICON-ART User Guide)
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
    
    ! Feedback processes (Details: cf. Tab. 2.5 ICON-ART User Guide)
    INTEGER :: iart_aci_warm           !< Nucleation of aerosol to cloud droplets
    INTEGER :: iart_aci_cold           !< Nucleation of aerosol to cloud ice
    INTEGER :: iart_ari                !< Direct interaction of aerosol with radiation
    
    ! Fast Physics Processes (Details: cf. Tab. 2.6 ICON-ART User Guide)
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
