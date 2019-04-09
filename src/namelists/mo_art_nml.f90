!>
!! Namelist for ART-package
!!
!! Subroutine is called by read_atmo_namelists for setting up the ART-package
!!
!! @author Daniel Reinert, DWD
!! @author Kristina Lundgren, KIT
!!
!! @par Revision History
!! Initial revision by Daniel Reinert, DWD (2011-12-08)
!! Modification of namelist parameters, Kristina Lundgren (2012-03-21)
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_art_nml
 
  USE mo_exception,           ONLY: message, finish
  USE mo_run_config,          ONLY: lart
  USE mo_io_units,            ONLY: nnml, nnml_output
  USE mo_impl_constants,      ONLY: max_dom
  USE mo_namelist,            ONLY: position_nml, POSITIONED, open_nml, close_nml
  USE mo_master_control,      ONLY: use_restart_namelists
  USE mo_mpi,                 ONLY: my_process_is_stdio
  USE mo_restart_namelist,    ONLY: open_tmpfile, store_and_close_namelist,     &
    &                               open_and_restore_namelist, close_tmpfile
  USE mo_art_config,          ONLY: art_config, IART_PATH_LEN
  USE mo_nml_annotate,        ONLY: temp_defaults, temp_settings

  USE mo_art_init_interface,  ONLY: art_calc_number_of_art_tracers_xml
  
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: read_art_namelist

  !----------------------------------!
  ! art_nml namelist variables       !
  !----------------------------------!

  ! General variables (Details: cf. Tab. 2.1 ICON-ART User Guide)
  CHARACTER(LEN=IART_PATH_LEN)  :: &
    &  cart_input_folder             !< Absolute Path to ART source code
  INTEGER :: iart_ntracer            !< number transported ART tracers
  INTEGER :: iart_init_aero          !< Initialization of aerosol species
  INTEGER :: iart_init_passive       !< Initialization of passive species
  INTEGER :: iart_init_gas           !< Initialization of gaseous species
  LOGICAL :: lart_diag_out           !< Enable output of diagnostic fields
  LOGICAL :: lart_pntSrc             !< Enables point sources
  LOGICAL :: lart_emiss_turbdiff     !< Switch if emissions should be included as surface flux condition
  CHARACTER(LEN=20) :: & 
   &  cart_io_suffix(1:max_dom)      !< user given suffix instead of automatically generated grid number 
                                     !  in ICON-ART input filename convention: 
                                     !  ART_iconR<n>B<kk>-grid-<yyyy-mm-dd-hh>_<grid_suffix>.nc

  ! Atmospheric Chemistry (Details: cf. Tab. 2.2 ICON-ART User Guide)
  LOGICAL :: lart_chem               !< Main switch to enable chemistry
  LOGICAL :: lart_passive            !< Main switch to enable passive tracers
  LOGICAL :: lart_psc                !< switch for computation of PSCs 
  INTEGER :: iart_chem_mechanism     !< Selects the chemical mechanism
  CHARACTER(LEN=IART_PATH_LEN)  :: &
    &  cart_vortex_init_date         !< Date of vortex initialization
  CHARACTER(LEN=IART_PATH_LEN)  :: &
    &  cart_cheminit_file(max_dom)   !< Path to chemical initialization file
  CHARACTER(LEN=IART_PATH_LEN)  :: &
    &  cart_cheminit_coord           !< Path to chemical initialization coordinate file
  CHARACTER(LEN=IART_PATH_LEN)  :: &
      &  cart_cheminit_type          !< type of chemical initialization coordinate file
  ! Paths and filenames of XML configuration
  CHARACTER(LEN=IART_PATH_LEN)  :: &
    &  cart_chemistry_xml            !< Path to XML file for chemical tracers
  CHARACTER(LEN=IART_PATH_LEN)  :: &
    &  cart_aerosol_xml              !< Path to XML file for aerosol tracers
  CHARACTER(LEN=IART_PATH_LEN)  :: &
    &  cart_passive_xml              !< Path to XML file for passive tracers
  CHARACTER(LEN=IART_PATH_LEN)  :: &
    &  cart_modes_xml                !< Path to XML file for modes
  CHARACTER(LEN=IART_PATH_LEN)  :: &
    &  cart_pntSrc_xml               !< Path to XML file for point sources
  CHARACTER(LEN=IART_PATH_LEN)  :: &
    &  cart_diagnostics_xml          !< Path to XML file for aerosol diagnostics (GRIB2 meta data)
  CHARACTER(LEN=IART_PATH_LEN)  :: &
    &  cart_emiss_xml_file           !< path and file name of the xml files for emission metadata
  CHARACTER(LEN=IART_PATH_LEN)  :: &
    &  cart_ext_data_xml             !< Path to XML file for metadata of datasets 
                                     !  that can prescribe tracers
  ! Atmospheric Aerosol (Details: cf. Tab. 2.3 ICON-ART User Guide)
  LOGICAL :: lart_aerosol            !< Main switch for the treatment of atmospheric aerosol
  INTEGER :: iart_seasalt            !< Treatment of sea salt aerosol
  INTEGER :: iart_dust               !< Treatment of mineral dust aerosol
  INTEGER :: iart_anthro             !< Treatment of anthropogenic aerosol
  INTEGER :: iart_fire               !< Treatment of wildfire aerosol
  INTEGER :: iart_volcano            !< Treatment of volcanic ash aerosol
  INTEGER :: iart_nonsph             !< Treatment of nonspherical particles
  CHARACTER(LEN=IART_PATH_LEN)  :: &
    &  cart_volcano_file             !< Absolute path + filename of input file for volcanoes
  INTEGER :: iart_radioact           !< Treatment of radioactive particles
  CHARACTER(LEN=IART_PATH_LEN)  :: &
    &  cart_radioact_file            !< Absolute path + filename of input file for radioactive emissions
  INTEGER :: iart_pollen             !< Treatment of pollen

  ! Feedback processes (Details: cf. Tab. 2.4 ICON-ART User Guide)
  INTEGER :: iart_aci_warm           !< Nucleation of aerosol to cloud droplets
  INTEGER :: iart_aci_cold           !< Nucleation of aerosol to cloud ice
  INTEGER :: iart_ari                !< Direct interaction of aerosol with radiation

  ! Fast Physics Processes (Details: cf. Tab. 2.5 ICON-ART User Guide)
  LOGICAL :: lart_conv               !< Convection of aerosol (TRUE/FALSE)
  LOGICAL :: lart_turb               !< Turbulent diffusion of aerosol (TRUE/FALSE)

  INTEGER :: iart_echam_ghg          !< integer for number tracers of hard coded echam greenhouse gases
  



  NAMELIST/art_nml/ cart_input_folder, lart_chem, lart_passive,                        &
   &                iart_chem_mechanism, cart_io_suffix, lart_pntSrc,                  &
   &                lart_aerosol, iart_seasalt, iart_dust, iart_anthro, iart_fire,     &
   &                iart_volcano, cart_volcano_file, iart_radioact,                    &
   &                cart_radioact_file, iart_pollen, iart_nonsph,                      &
   &                iart_aci_warm, iart_aci_cold, iart_ari,                            &
   &                lart_conv, iart_ntracer, lart_turb, iart_init_aero, iart_init_gas, &
   &                lart_diag_out, cart_emiss_xml_file, cart_ext_data_xml,             &
   &                cart_vortex_init_date , cart_cheminit_file, cart_cheminit_coord,   &
   &                cart_cheminit_type,                                                &
   &                lart_emiss_turbdiff,                                               &
   &                cart_chemistry_xml, cart_aerosol_xml, cart_passive_xml,            &
   &                cart_modes_xml, cart_pntSrc_xml, cart_diagnostics_xml,             &
   &                iart_init_passive, lart_psc, iart_echam_ghg

CONTAINS
  !-------------------------------------------------------------------------
  !>
  !! Read Namelist for ART-package.
  !!
  !! This subroutine
  !! - reads the Namelist for the ART-package
  !! - sets default values
  !! - potentially overwrites the defaults by values used in a
  !!   previous integration (if this is a resumed run)
  !! - reads the user's (new) specifications
  !! - stores the Namelist for restart
  !! - fills the configuration state (partly)
  !!
  !! @par Revision History
  !!  by Daniel Reinert, DWD (2011-12-08)
  !!  by Daniel Rieger,  KIT (2014-17-06)
  !!
  SUBROUTINE read_art_namelist( filename )

    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER :: istat, funit
    INTEGER :: jg          !< patch loop index
    INTEGER :: auto_ntracer      !< automatically computed number of tracers
    INTEGER :: auto_ntracer_xml  !< art ntracer from one xml file
    LOGICAL :: l_exist     !< variable for inquiring if the xml file 
                                        !   and emission base path exist.
    CHARACTER(len=*), PARAMETER ::  &
      &  routine = 'mo_art_nml: read_art_nml'
    INTEGER :: iunit

    !-----------------------
    ! 1. default settings
    !-----------------------

    ! General variables (Details: cf. Tab. 2.1 ICON-ART User Guide)
    cart_input_folder          = ''
    iart_ntracer               = -1    !< default value if it is not given
    iart_init_aero             = 0
    iart_init_passive          = 0
    iart_init_gas              = 0
    lart_diag_out              = .FALSE.
    lart_pntSrc                = .FALSE.
    lart_emiss_turbdiff        = .FALSE.
    cart_io_suffix(1:max_dom)  = 'grid-number'

    ! Atmospheric Chemistry (Details: cf. Tab. 2.2 ICON-ART User Guide)
    lart_chem             = .FALSE.
    lart_passive          = .FALSE.
    lart_psc              = .FALSE.
    iart_chem_mechanism   = 0
    cart_vortex_init_date = ''
    cart_cheminit_file(:) = ''
    cart_cheminit_coord   = ''
    cart_cheminit_type    = ''

    ! Paths and filenames of XML configuration
    cart_chemistry_xml    = ''
    cart_aerosol_xml      = ''
    cart_passive_xml      = ''
    cart_modes_xml        = ''
    cart_pntSrc_xml       = ''
    cart_diagnostics_xml  = ''
    cart_emiss_xml_file   = ''
    cart_ext_data_xml     = ''

    ! Atmospheric Aerosol (Details: cf. Tab. 2.3 ICON-ART User Guide)
    lart_aerosol        = .FALSE.
    iart_seasalt        = 0
    iart_dust           = 0
    iart_anthro         = 0
    iart_fire           = 0
    iart_volcano        = 0
    cart_volcano_file   = ''
    iart_radioact       = 0
    cart_radioact_file  = ''
    iart_pollen         = 0
    iart_nonsph         = 0

    ! Feedback processes (Details: cf. Tab. 2.4 ICON-ART User Guide)
    iart_aci_warm       = 0
    iart_aci_cold       = 0
    iart_ari            = 0

    ! Fast Physics Processes (Details: cf. Tab. 2.5 ICON-ART User Guide)
    lart_conv           = .TRUE.
    lart_turb           = .TRUE.

    iart_echam_ghg      = 0

    !------------------------------------------------------------------
    ! 2. If this is a resumed integration, overwrite the defaults above
    !    by values used in the previous integration.
    !------------------------------------------------------------------
    IF (use_restart_namelists()) THEN
      funit = open_and_restore_namelist('art_nml')
      READ(funit,NML=art_nml)
      CALL close_tmpfile(funit)
    END IF

    !--------------------------------------------------------------------
    ! 3. Read user's (new) specifications (Done so far by all MPI processes)
    !--------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('art_nml', STATUS=istat)
    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, art_nml)    ! write defaults to temporary text file
    END IF
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, art_nml)                                        ! overwrite default settings
      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, art_nml)    ! write settings to temporary text file
      END IF
    END SELECT
    CALL close_nml



    !----------------------------------------------------
    ! 4. Sanity check (only if lart is true)
    !----------------------------------------------------

    IF (lart) THEN

    
      IF (iart_aci_cold == 6 .AND. iart_dust == 0) THEN
        CALL finish('mo_art_nml:read_art_namelist',  &
          &         'Invalid combination: iart_aci_cold = 6 and iart_dust = 0')
      ENDIF
      IF (iart_aci_cold == 7 .AND. iart_dust == 0) THEN
        CALL finish('mo_art_nml:read_art_namelist',  &
          &         'Invalid combination: iart_aci_cold = 7 and iart_dust = 0')
      ENDIF
  
      ! Emission paths and file
      IF (TRIM(cart_emiss_xml_file) /= '') THEN
        INQUIRE(file = TRIM(cart_emiss_xml_file), EXIST = l_exist)
      
        IF (.NOT. l_exist) THEN
          CALL finish('mo_art_nml:read_art_namelist',  &
                      TRIM(cart_emiss_xml_file)//  &
                      & ' could not be found. Check cart_emiss_xml_file.')
        END IF
      END IF
  
      auto_ntracer = 0
      ! chemistry xml file
      IF (lart_chem) THEN
        IF (TRIM(cart_chemistry_xml) == '') THEN
          CALL finish('mo_art_nml:read_art_namelist','namelist parameter cart_chemistry_xml' &
                    //' has to be given for lart_chem == .TRUE.')
        ELSE
          INQUIRE(file = TRIM(cart_chemistry_xml), EXIST = l_exist)
  
          IF (l_exist) THEN
            CALL art_calc_number_of_art_tracers_xml(TRIM(cart_chemistry_xml),  &
                                   &                auto_ntracer_xml)
            auto_ntracer = auto_ntracer + auto_ntracer_xml
          ELSE
            CALL finish('mo_art_nml:read_art_namelist',  &
                        TRIM(cart_chemistry_xml)//  &
                        & ' could not be found. Check cart_chemistry_xml.')
          END IF
        END IF
      END IF
  
      ! aerosol xml file
      IF (lart_aerosol) THEN
        IF (TRIM(cart_aerosol_xml) == '') THEN
          CALL finish('mo_art_nml:read_art_namelist','namelist parameter cart_aerosol_xml' &
                    //' has to be given for lart_aerosol == .TRUE.')
        ELSE
          INQUIRE(file = TRIM(cart_aerosol_xml), EXIST = l_exist)
  
          IF (l_exist) THEN
            CALL art_calc_number_of_art_tracers_xml(TRIM(cart_aerosol_xml),  &
                                   &                auto_ntracer_xml)
            auto_ntracer = auto_ntracer + auto_ntracer_xml
          ELSE
            CALL finish('mo_art_nml:read_art_namelist',  &
                        TRIM(cart_aerosol_xml)//  &
                        & ' could not be found. Check cart_aerosol_xml.')
          END IF
        END IF
      END IF
  
      ! passive xml file
      IF (lart_passive) THEN
        IF (TRIM(cart_passive_xml) == '') THEN
          CALL finish('mo_art_nml:read_art_namelist','namelist parameter cart_passive_xml' &
                    //' has to be given for lart_passive == .TRUE.')
        ELSE
          INQUIRE(file = TRIM(cart_passive_xml), EXIST = l_exist)
  
          IF (l_exist) THEN
            CALL art_calc_number_of_art_tracers_xml(TRIM(cart_passive_xml),  &
                                   &                auto_ntracer_xml)
            auto_ntracer = auto_ntracer + auto_ntracer_xml
          ELSE
            CALL finish('mo_art_nml:read_art_namelist',  &
                        TRIM(cart_passive_xml)//  &
                        & ' could not be found. Check cart_passive_xml.')
          END IF
        END IF
      END IF

      IF (iart_ntracer > -1) THEN
        CALL message('WARNING',  &
          &          'Namelist parameter iart_ntracer of art_nml is obsolete '  &
          &        //'and will be removed soon.')
      END IF


      IF ((iart_ntracer > -1)  &
         &  .AND. (auto_ntracer /= iart_ntracer)) THEN
        CALL finish('mo_art_nml:read_art_namelist',                              &
              &     'The given namelist parameter iart_ntracer is not equal to ' &
              &   //'the automatically computed one. This namelist parameter '   &
              &   //'is obsolete so just remove it from your art_nml.')
      END IF

            

    END IF  ! lart


    ! Diagnostics paths and file
    IF (TRIM(cart_diagnostics_xml) /= '') THEN
      INQUIRE(file = TRIM(cart_diagnostics_xml), EXIST = l_exist)

      IF (.NOT. l_exist) THEN
        CALL finish('mo_art_nml:read_art_namelist',  &
                    TRIM(cart_diagnostics_xml)//  &
                    & ' could not be found. Check cart_diagnostics_xml.')
      END IF
    END IF

    !----------------------------------------------------
    ! 5. Fill the configuration state
    !----------------------------------------------------

    DO jg= 1,max_dom !< Do not take into account reduced radiation grid
      ! General variables (Details: cf. Tab. 2.1 ICON-ART User Guide)
      art_config(jg)%cart_input_folder   = TRIM(cart_input_folder)
      art_config(jg)%iart_init_aero      = iart_init_aero
      art_config(jg)%iart_init_gas       = iart_init_gas
      art_config(jg)%iart_init_passive   = iart_init_passive
      art_config(jg)%lart_diag_out       = lart_diag_out
      art_config(jg)%lart_pntSrc         = lart_pntSrc
      art_config(jg)%lart_emiss_turbdiff = lart_emiss_turbdiff
      art_config(jg)%cart_io_suffix      = TRIM(cart_io_suffix(jg))

      ! Atmospheric Chemistry (Details: cf. Tab. 2.2 ICON-ART User Guide)
      art_config(jg)%lart_chem             = lart_chem
      art_config(jg)%lart_passive          = lart_passive
      art_config(jg)%iart_chem_mechanism   = iart_chem_mechanism
      art_config(jg)%lart_psc              = lart_psc
      art_config(jg)%cart_vortex_init_date = TRIM(cart_vortex_init_date)
      art_config(jg)%cart_cheminit_file    = TRIM(cart_cheminit_file(jg))
      art_config(jg)%cart_cheminit_coord   = TRIM(cart_cheminit_coord)
      art_config(jg)%cart_cheminit_type    = TRIM(cart_cheminit_type)

      ! Paths and filenames of XML configuration
      art_config(jg)%cart_chemistry_xml    = TRIM(cart_chemistry_xml)
      art_config(jg)%cart_aerosol_xml      = TRIM(cart_aerosol_xml)
      art_config(jg)%cart_passive_xml      = TRIM(cart_passive_xml)
      art_config(jg)%cart_modes_xml        = TRIM(cart_modes_xml)
      art_config(jg)%cart_pntSrc_xml       = TRIM(cart_pntSrc_xml)
      art_config(jg)%cart_diagnostics_xml  = TRIM(cart_diagnostics_xml)
      art_config(jg)%cart_emiss_xml_file   = TRIM(cart_emiss_xml_file)
      art_config(jg)%cart_ext_data_xml     = TRIM(cart_ext_data_xml)


      ! Atmospheric Aerosol (Details: cf. Tab. 2.3 ICON-ART User Guide)
      art_config(jg)%lart_aerosol        = lart_aerosol
      art_config(jg)%iart_seasalt        = iart_seasalt
      art_config(jg)%iart_dust           = iart_dust
      art_config(jg)%iart_anthro         = iart_anthro
      art_config(jg)%iart_fire           = iart_fire
      art_config(jg)%iart_volcano        = iart_volcano
      art_config(jg)%iart_nonsph         = iart_nonsph
      art_config(jg)%cart_volcano_file   = TRIM(cart_volcano_file)
      art_config(jg)%iart_radioact       = iart_radioact
      art_config(jg)%cart_radioact_file  = TRIM(cart_radioact_file)
      art_config(jg)%iart_pollen         = iart_pollen

     ! Feedback processes (Details: cf. Tab. 2.4 ICON-ART User Guide)
      art_config(jg)%iart_aci_warm       = iart_aci_warm
      art_config(jg)%iart_aci_cold       = iart_aci_cold
      art_config(jg)%iart_ari            = iart_ari

      ! Fast Physics Processes (Details: cf. Tab. 2.5 ICON-ART User Guide)
      art_config(jg)%lart_conv           = lart_conv
      art_config(jg)%lart_turb           = lart_turb

      ! art number of tracers
      art_config(jg)%iart_ntracer        = auto_ntracer 


      art_config(jg)%iart_echam_ghg      = iart_echam_ghg 
    ENDDO !jg

    !-----------------------------------------------------
    ! 6. Store the namelist for restart
    !-----------------------------------------------------
    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=art_nml)
      CALL store_and_close_namelist(funit, 'art_nml')
    ENDIF

    ! 7. write the contents of the namelist to an ASCII file
    !
    IF(my_process_is_stdio()) WRITE(nnml_output,nml=art_nml)


  END SUBROUTINE read_art_namelist

END MODULE mo_art_nml
