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

  PUBLIC :: ctracer_art                                                  !generic tracer list, contains name of tracer

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
