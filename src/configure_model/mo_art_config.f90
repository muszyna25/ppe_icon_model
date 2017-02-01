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

  USE mo_impl_constants,       ONLY: max_dom
  USE mo_impl_constants,       ONLY: MAX_CHAR_LENGTH
  
  IMPLICIT NONE
  
  PRIVATE

  PUBLIC :: nart_tendphy
  PUBLIC :: npreslay
  PUBLIC :: t_art_config
  PUBLIC :: art_config
  PUBLIC :: configure_art
  PUBLIC :: IART_PATH_LEN
  ! Running tracer indices in ICON-ART
  PUBLIC :: iash1,iash2,iash3,iash4,iash5,iash6                          !Volcanic Ash 
  PUBLIC :: iasha, iashb, iashc, iasha0, iashb0, iashc0                  !Volcanic Ash
  PUBLIC :: iCS137,iI131,iTE132,iZR95,iXE133,iI131g,iI131o,iBA140,iRU103 !radioactive nuclides
  PUBLIC :: iseasa,iseasb,iseasc,iseasa0,iseasb0,iseasc0                 !sea salt 
  PUBLIC :: idusta,idustb,idustc,idusta0,idustb0,idustc0                 !mineral dust
  PUBLIC :: iTRCHBR3,iTRCH2BR2                                           !chemical tracer - VSLS
  PUBLIC :: iTRCH4,iTRC2H6,iTRC3H8,iTRC5H8,iTRCH3COCH3,iTRCO,iTRCO2      !chemical tracer - CH4-C2H6-C3H5-C5H8-CH3COCH3-CO-CO2
  PUBLIC :: iTRH2O,iTRO3,iTRN2O,iTRNH3,iTRSO2,iTRH2SO4,iTRHNO3,iTRAGE    !chemical tracer - others
  PUBLIC :: iTR_vortex,iTR_stn,iTR_stt,iTR_sts                           !artificial tracer
  PUBLIC :: iTR_trn,iTR_trt,iTR_trs,iTR_tiln,iTR_tils                    !artificial tracer
  PUBLIC :: iTR_nh,iTR_sh                                                !artificial tracer
  PUBLIC :: iTR_nin,iTR_sin,iTR_ech,iTR_sea,iTR_sib,iTR_eur              !artificial tracer 
  PUBLIC :: iTR_med,iTR_naf,iTR_saf,iTR_mdg,iTR_aus,iTR_nam              !artificial tracer
  PUBLIC :: iTR_sam,iTR_tpo,iTR_tao,iTR_tio,iTR_bgn,iTR_bgs              !artificial tracer
  PUBLIC :: iTR_art                                                      !artificial tracer
  
  PUBLIC :: ctracer_art                                                  !generic tracer list, contains name of tracer
  
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
    &  iTRCHBR3,iTRCH2BR2,                                        &
    &  iTRCH4,iTRC2H6,iTRC3H8,iTRC5H8,iTRCH3COCH3,iTRCO,iTRCO2,   &
    &  iTRH2O,iTRO3,iTRN2O,iTRNH3,iTRSO2,iTRH2SO4,iTRHNO3,iTRAGE, &
    &  iTR_vortex,iTR_stn,iTR_stt,iTR_sts,                        &
    &  iTR_trn,iTR_trt,iTR_trs,iTR_tiln,iTR_tils,                 &
    &  iTR_nh,iTR_sh,                                             &
    &  iTR_nin,iTR_sin,iTR_ech,iTR_sea,iTR_sib,iTR_eur,           &
    &  iTR_med,iTR_naf,iTR_saf,iTR_mdg,iTR_aus,iTR_nam,           &
    &  iTR_sam,iTR_tpo,iTR_tao,iTR_tio,iTR_bgn,iTR_bgs,           &
    &  iTR_art
  
  !!--------------------------------------------------------------------------
  !! Basic configuration setup for ICON-ART
  !!--------------------------------------------------------------------------
  INTEGER, PARAMETER  :: npreslay      = 7   !Number of pressure layers for diagnostic output of maximum concentration
  INTEGER, PARAMETER  :: IART_PATH_LEN = 200 !Maximum length of file- and pathnames

  INTEGER             :: nart_tendphy  = 0  !Maximum number of tracers that are effected by deep convective transport 
  
  
  TYPE t_art_config ! Namelist variables for ART
    
    ! Namelist variables
    
    ! General control variables (Details: cf. Tab. 2.2 ICON-ART User Guide)
    CHARACTER(LEN=IART_PATH_LEN) :: &
      &  cart_input_folder               !< Absolute Path to ART input files 
    INTEGER :: iart_ntracer              !< number of transported ART tracers
    INTEGER :: iart_init_aero            !< Initialization of aerosol species
    INTEGER :: iart_init_gas             !< Initialization of gaseous species
    INTEGER :: iart_init_passive         !< Initialization of passive species
    LOGICAL :: lart_diag_out             !< Enable output of diagnostic fields
    LOGICAL :: lart_pntSrc               !< Enables point sources
    CHARACTER(LEN=20) :: cart_io_suffix  !< user given suffix instead of automatically generated grid number 
                                         !  in ICON-ART input filename convention: 
                                         !  ART_iconR<n>B<kk>-grid-<yyyy-mm-dd-hh>_<grid_suffix>.nc
    
    ! Atmospheric Chemistry (Details: cf. Tab. 2.3 ICON-ART User Guide)
    LOGICAL :: lart_chem               !< Main switch to enable chemistry
    LOGICAL :: lart_passive            !< Main switch to enable chemistry
    INTEGER :: iart_chem_mechanism     !< Selects the chemical mechanism
    INTEGER :: iart_psc                !< integer which indicates how to compute PSCs (0: no PSCs, 1: climatology, 2: online)
    CHARACTER(LEN=IART_PATH_LEN) :: &
      &  cart_emiss_xml_file          !< Path and file name of the xml file containing meta information of the emissions.
    CHARACTER(LEN=IART_PATH_LEN) :: &
      &  cart_vortex_init_date        !< Date of vortex initialization
    CHARACTER(LEN=IART_PATH_LEN) :: &
      &  cart_cheminit_file           !< Path to chemical initialization file
    CHARACTER(LEN=IART_PATH_LEN) :: &
      &  cart_cheminit_coord          !< Path to chemical initialization coordinate file
    CHARACTER(LEN=IART_PATH_LEN) :: &
      &  cart_cheminit_type           !< Type of  chemical initialization coordinate file
    ! Paths and filenames of XML configuration
    CHARACTER(LEN=IART_PATH_LEN) :: &
      &  cart_chemistry_xml           !< Path to XML file for chemical tracers
    CHARACTER(LEN=IART_PATH_LEN) :: &
      &  cart_aerosol_xml             !< Path to XML file for aerosol tracers
    CHARACTER(LEN=IART_PATH_LEN) :: &
      &  cart_passive_xml             !< Path to XML file for passive tracers
    CHARACTER(LEN=IART_PATH_LEN) :: &
      &  cart_modes_xml               !< Path to XML file for modes
    CHARACTER(LEN=IART_PATH_LEN) :: &
      &  cart_pntSrc_xml              !< Path to XML file for point sources

    ! Atmospheric Aerosol (Details: cf. Tab. 2.4 ICON-ART User Guide)
    LOGICAL :: lart_aerosol            !< Main switch for the treatment of atmospheric aerosol
    INTEGER :: iart_seasalt            !< Treatment of sea salt aerosol
    INTEGER :: iart_dust               !< Treatment of mineral dust aerosol
    INTEGER :: iart_anthro             !< Treatment of anthropogenic aerosol
    INTEGER :: iart_fire               !< Treatment of wildfire aerosol
    INTEGER :: iart_volcano            !< Treatment of volcanic ash aerosol
    CHARACTER(LEN=IART_PATH_LEN) :: &
      &  cart_volcano_file             !< Absolute path + filename of input file for volcanoes
    INTEGER :: iart_radioact           !< Treatment of radioactive particles
    CHARACTER(LEN=IART_PATH_LEN) :: &
      &  cart_radioact_file            !< Absolute path + filename of input file for radioactive emissions
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


    !------------- generic tracer list for NWP and ECHAM physics
  CHARACTER(len=MAX_CHAR_LENGTH), ALLOCATABLE :: ctracer_art(:) !< tracer acronyms

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
