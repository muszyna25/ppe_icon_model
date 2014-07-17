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

  USE mo_exception,           ONLY: message, finish, message_text
  USE mo_kind                ,ONLY: wp
  USE mo_parallel_config     ,ONLY: nproma
  USE mo_exception           ,ONLY: finish
  USE mo_io_units            ,ONLY: nnml, nnml_output
  USE mo_master_control      ,ONLY: is_restart_run
  USE mo_impl_constants      ,ONLY: max_dom
  USE mo_namelist            ,ONLY: position_nml, POSITIONED, open_nml, close_nml
  USE mo_mpi                 ,ONLY: my_process_is_stdio
  USE mo_io_restart_namelist ,ONLY: open_tmpfile, store_and_close_namelist,     &
    &                               open_and_restore_namelist, close_tmpfile
  USE mo_art_config          ,ONLY: art_config ,t_volc_list,max_volc_input
  USE mo_nml_annotate,   ONLY: temp_defaults, temp_settings

  
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: read_art_namelist

  !----------------------------------!
  ! art_nml namelist variables       !
  !----------------------------------!
  
  ! General variables (Details: cf. Tab. 2.1 ICON-ART User Guide)
  CHARACTER(LEN=120) :: cart_folder  !< Absolute Path to ART source code
  INTEGER :: iart_ntracer            !< number transported ART tracers
    
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
  LOGICAL :: lart_turb               !< Turbulent diffusion of aerosol (TRUE/FALSE)

  NAMELIST/art_nml/ cart_folder, lart_chem, iart_chem_mechanism,                       &
   &                lart_aerosol, iart_seasalt, iart_dust, iart_anthro, iart_fire,     &
   &                iart_volcano, cart_volcano_file, iart_radioact,                    &
   &                cart_radioact_file, iart_pollen,                                   &
   &                iart_aci_warm, iart_aci_cold, iart_ari,                            &
   &                lart_conv, lart_turb, iart_ntracer

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
    INTEGER :: jb, jc, nblks, npromz, ivolc   
    CHARACTER(len=*), PARAMETER ::  &
      &  routine = 'mo_art_nml: read_art_nml'
    INTEGER :: iunit

    !-----------------------
    ! 1. default settings   
    !-----------------------
      
    ! General variables (Details: cf. Tab. 2.1 ICON-ART User Guide)
    cart_folder         = './art/'
    iart_ntracer        = 0
      
    ! Atmospheric Chemistry (Details: cf. Tab. 2.2 ICON-ART User Guide)
    lart_chem           = .FALSE.
    iart_chem_mechanism = 0
      
    ! Atmospheric Aerosol (Details: cf. Tab. 2.3 ICON-ART User Guide)
    lart_aerosol        = .FALSE.
    iart_seasalt        = 0
    iart_dust           = 0
    iart_anthro         = 0
    iart_fire           = 0
    iart_volcano        = 0
    cart_volcano_file   = './volcanofile'
    iart_radioact       = 0
    cart_radioact_file  = './radioactfile'
    iart_pollen         = 0
      
    ! Feedback processes (Details: cf. Tab. 2.4 ICON-ART User Guide)
    iart_aci_warm       = 0
    iart_aci_cold       = 0
    iart_ari            = 0
      
    ! Fast Physics Processes (Details: cf. Tab. 2.5 ICON-ART User Guide)
    lart_conv           = .TRUE.
    lart_turb           = .TRUE.
    
    !------------------------------------------------------------------
    ! 2. If this is a resumed integration, overwrite the defaults above 
    !    by values used in the previous integration.
    !------------------------------------------------------------------
    IF (is_restart_run()) THEN
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
    ! 4. Sanity check
    !----------------------------------------------------



    !----------------------------------------------------
    ! 5. Fill the configuration state
    !----------------------------------------------------
    
    DO jg= 0,max_dom
      ! General variables (Details: cf. Tab. 2.1 ICON-ART User Guide)
      art_config(jg)%cart_folder         = TRIM(cart_folder)
      art_config(jg)%iart_ntracer        = iart_ntracer
    
      ! Atmospheric Chemistry (Details: cf. Tab. 2.2 ICON-ART User Guide)
      art_config(jg)%lart_chem           = lart_chem
      art_config(jg)%iart_chem_mechanism = iart_chem_mechanism
      
      ! Atmospheric Aerosol (Details: cf. Tab. 2.3 ICON-ART User Guide)
      art_config(jg)%lart_aerosol        = lart_aerosol
      art_config(jg)%iart_seasalt        = iart_seasalt
      art_config(jg)%iart_dust           = iart_dust
      art_config(jg)%iart_anthro         = iart_anthro
      art_config(jg)%iart_fire           = iart_fire
      art_config(jg)%iart_volcano        = iart_volcano
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
