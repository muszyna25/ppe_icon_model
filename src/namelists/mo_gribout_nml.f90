!>
!! Namelist for Grib output
!!
!! These subroutines are called by  read_atmo_namelists and do the transport 
!! setup.
!!
!! @author Daniel Reinert, DWD
!!
!!
!! @par Revision History
!! Initial revision by Daniel Reinert, DWD (2011-04-20)
!! - moved here from mo_advection_utils
!! Modification by Daniel Reinert, DWD (2011-04-20)
!! - some updates on the structure
!!
!! @par Copyright
!! 2002-2010 by DWD and MPI-M
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
MODULE mo_gribout_nml

  USE mo_io_units,            ONLY: nnml, nnml_output
  USE mo_impl_constants,      ONLY: max_dom
  USE mo_master_control,      ONLY: is_restart_run
  USE mo_namelist,            ONLY: position_nml, POSITIONED, open_nml, close_nml
  USE mo_mpi,                 ONLY: my_process_is_stdio
  USE mo_io_restart_namelist, ONLY: open_tmpfile, store_and_close_namelist,     &
    &                               open_and_restore_namelist, close_tmpfile
  USE mo_gribout_config,      ONLY: gribout_config 
  USE mo_nml_annotate,        ONLY: temp_defaults, temp_settings
  USE mo_util_string,         ONLY: int2string
  USE mo_exception,           ONLY: finish, message

  
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: read_gribout_namelist


  !> module name
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_output_event_handler'

  ! constant for better readability
  INTEGER, PARAMETER :: UNDEFINED = -1

  !----------------------------------!
  ! gribout_nml namelist variables   !
  !----------------------------------!

  ! Namelist parameter "preset": main switch (string), possible options are
  !  - "none"
  !  - "deterministic"
  !  - "ensemble"
  !
  ! Setting this different to "none" enables a couple of defaults for
  ! the other gribout_nml namelist parameters. If, additionally, the
  ! user tries to set any of these other parameters to a conflicting
  ! value, an error message is thrown.
  CHARACTER(LEN=32) :: preset


  INTEGER :: &                          ! Table 1.2
    & significanceOfReferenceTime       ! 0: Analysis
                                        ! 1: Start of forecast
                                        ! 2: Verifying time of forecast
                                        ! 4: ... 

  INTEGER :: &                          ! Table 1.3
    & productionStatusOfProcessedData   ! 0: Oper. products
                                        ! 1: Oper. test products
                                        ! 2: Research products
                                        ! 3: ...

  INTEGER :: &                          ! Table 1.4
    & typeOfProcessedData               ! 0: Analysis products
                                        ! 1: Forecast products
                                        ! 2: Analysis and forecast products
                                        ! 3: ...

  INTEGER :: &                          ! Table 4.3
    & typeOfGeneratingProcess           ! 0: Analysis
                                        ! 1: Initialization
                                        ! 2: Forecast
                                        ! 3: ...

  INTEGER :: &                          ! Table: backgroundProcess
    & backgroundProcess                 ! 0: main run
                                        ! 1: pre-assimilation
                                        ! 2: assimilation
                                        ! 3: ... 

  INTEGER :: &                          ! Table: generatingProcessIdentifier
    & generatingProcessIdentifier(0:max_dom) ! 1: icogl
                                        ! 2: icrgl
                                        ! 3: icoeu
                                        ! 4: ...

  INTEGER :: &                          ! Table: local.78.254.def
    & localDefinitionNumber             ! 252: Ensemble system incl. postprocessing
                                        ! 253: Ensemble system 
                                        ! 254: Deterministic system

  INTEGER :: &                          ! Table: local.78.254.def
    & localNumberOfExperiment           !  


  INTEGER :: &                          ! Output generating center
    & generatingCenter                  !  


  INTEGER :: &                          ! Output generating subcenter
    & generatingSubcenter               !  


  LOGICAL :: ldate_grib_act             ! add Creation date to GRIB file
                                        ! .TRUE. : activated
                                        ! .FALSE.: deactivated (use dummy date/time) 

  LOGICAL :: lgribout_24bit             ! write thermodynamic fields rho, theta_v, T, p
                                        ! with 24bit precision

  ! Local definition for ensemble products
  INTEGER :: productDefinitionTemplateNumber

  INTEGER :: typeOfEnsembleForecast,        &
    &        localTypeOfEnsembleForecast,   &
    &        numberOfForecastsInEnsemble,   &
    &        perturbationNumber


  NAMELIST/gribout_nml/  &
    &                    preset,                          &
    &                    significanceOfReferenceTime,     &
    &                    productionStatusOfProcessedData, &
    &                    typeOfProcessedData,             &
    &                    typeOfGeneratingProcess,         &
    &                    backgroundProcess,               &
    &                    generatingProcessIdentifier,     &
    &                    localDefinitionNumber,           &
    &                    localNumberOfExperiment,         &
    &                    generatingCenter,                &
    &                    generatingSubcenter,             &
    &                    ldate_grib_act,                  &
    &                    productDefinitionTemplateNumber, &
    &                    typeOfEnsembleForecast,          &
    &                    localTypeOfEnsembleForecast,     &
    &                    numberOfForecastsInEnsemble,     &
    &                    perturbationNumber,              &
    &                    lgribout_24bit


CONTAINS


  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Read Namelist for gribout. 
  !!
  !! This subroutine 
  !! - reads the Namelist for gribout
  !! - sets default values
  !! - potentially overwrites the defaults by values used in a 
  !!   previous integration (if this is a resumed run)
  !! - reads the user's (new) specifications
  !! - stores the Namelist for restart
  !! - fills the configuration state (partly)   
  !!
  !! @par Revision History
  !!  by Daniel Reinert, DWD (2013-01-29)
  !!
  SUBROUTINE read_gribout_namelist( filename )

    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER :: istat, funit
    INTEGER :: jg          !< patch loop index

    CHARACTER(len=*), PARAMETER ::  &
      &  routine = 'mo_gribout_nml: read_gribout_nml'

    !-----------------------
    ! 1. default settings   
    !-----------------------
    preset                               = "deterministic"

    significanceOfReferenceTime          = 1   ! 1: Start of forecast
    productionStatusOfProcessedData      = 1   ! 1: Oper. test products
    backgroundProcess                    = 0   ! 0: main run
    generatingProcessIdentifier(:)       = 1   ! 1: icogl
    localNumberOfExperiment              = 1
    ldate_grib_act                       = .TRUE.
    lgribout_24bit                       = .FALSE.  ! use 16bit precision for all fields

    typeOfProcessedData                  = UNDEFINED
    typeOfGeneratingProcess              = UNDEFINED
    localDefinitionNumber                = UNDEFINED
    generatingCenter                     = UNDEFINED  ! output generating center
    generatingSubcenter                  = UNDEFINED  ! output generating subcenter 
    productDefinitionTemplateNumber      = UNDEFINED  ! (undefined, will not be set if unchanged)
    typeOfEnsembleForecast               = UNDEFINED  ! (undefined, will not be set if unchanged)
    localTypeOfEnsembleForecast          = UNDEFINED  ! (undefined, will not be set if unchanged)
    numberOfForecastsInEnsemble          = UNDEFINED  ! (undefined, will not be set if unchanged)
    perturbationNumber                   = UNDEFINED  ! (undefined, will not be set if unchanged)

    !------------------------------------------------------------------
    ! 2. If this is a resumed integration, overwrite the defaults above 
    !    by values used in the previous integration.
    !------------------------------------------------------------------
    IF (is_restart_run()) THEN
      funit = open_and_restore_namelist('gribout_nml')
      READ(funit,NML=gribout_nml)
      CALL close_tmpfile(funit)
    END IF

    !--------------------------------------------------------------------
    ! 3. Read user's (new) specifications (Done so far by all MPI processes)
    !--------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('gribout_nml', STATUS=istat)
    IF (my_process_is_stdio()) WRITE(temp_defaults(), gribout_nml)  ! write defaults to temporary text file
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, gribout_nml, iostat=istat)                          ! overwrite default settings
      ! Preset values, when main switch is provided.
      CALL preset_namelist()
      IF (my_process_is_stdio()) WRITE(temp_settings(), gribout_nml)  ! write settings to temporary text file
    CASE DEFAULT
      CALL preset_namelist()
    END SELECT
    CALL close_nml

    !----------------------------------------------------
    ! 5. Fill the configuration state
    !----------------------------------------------------

    DO jg= 0,max_dom
      gribout_config(jg)%significanceOfReferenceTime       = &
        &                significanceOfReferenceTime       
      gribout_config(jg)%productionStatusOfProcessedData   = &
        &                productionStatusOfProcessedData   
      gribout_config(jg)%typeOfProcessedData               = &
        &                typeOfProcessedData               
      gribout_config(jg)%typeOfGeneratingProcess           = &
        &                typeOfGeneratingProcess           
      gribout_config(jg)%backgroundProcess                 = &
        &                backgroundProcess                 
      gribout_config(jg)%generatingProcessIdentifier       = &
        &                generatingProcessIdentifier(jg)   
      gribout_config(jg)%localDefinitionNumber             = &
        &                localDefinitionNumber             
      gribout_config(jg)%localNumberOfExperiment           = &
        &                localNumberOfExperiment           
      gribout_config(jg)%generatingCenter                  = &
        &                generatingCenter                  
      gribout_config(jg)%generatingSubcenter               = &
        &                generatingSubcenter               
      gribout_config(jg)%ldate_grib_act                    = &
        &                ldate_grib_act                    
      gribout_config(jg)%productDefinitionTemplateNumber   = &
        &                productDefinitionTemplateNumber       
      gribout_config(jg)%typeOfEnsembleForecast            = &
        &                typeOfEnsembleForecast       
      gribout_config(jg)%localTypeOfEnsembleForecast       = &
        &                localTypeOfEnsembleForecast
      gribout_config(jg)%numberOfForecastsInEnsemble       = &
        &                numberOfForecastsInEnsemble
      gribout_config(jg)%perturbationNumber                = &
        &                perturbationNumber
      gribout_config(jg)%lgribout_24bit                    = &
        &                lgribout_24bit
    ENDDO



    !-----------------------------------------------------
    ! 6. Store the namelist for restart
    !-----------------------------------------------------
    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=gribout_nml)                    
      CALL store_and_close_namelist(funit, 'gribout_nml')             
    ENDIF

    ! 7. write the contents of the namelist to an ASCII file
    !
    IF(my_process_is_stdio()) WRITE(nnml_output,nml=gribout_nml)


  END SUBROUTINE read_gribout_namelist


  !> Preset values, when main switch is provided.
  !
  SUBROUTINE preset_namelist()
    IF (TRIM(preset) == "deterministic") THEN
      !
      ! deterministic forecast
      !
      ! 2: Forecast
      CALL preset_value("typeOfGeneratingProcess",typeOfGeneratingProcess,             2 , quiet=.TRUE.)
      ! 254: Deterministic system
      CALL preset_value("localDefinitionNumber",  localDefinitionNumber,             254 , quiet=.TRUE. )
      ! 1: Forecast products
      CALL preset_value("typeOfProcessedData",    typeOfProcessedData,                 1 , quiet=.TRUE. )

    ELSE IF (TRIM(preset) == "ensemble") THEN
      !
      ! ensemble forecast
      !
      ! values are preset according to
      !  GME   : pp_makepdt.f90
      !  COSMO : io_metadata.f90 
      !
      ! The user is expected to set only
      !
      ! - perturbationNumber          (COSMO: "iepsmem")
      ! - numberOfForecastsInEnsemble (COSMO: "iepstot")
      ! - localTypeOfEnsembleForecast (COSMO: "iepstyp")
      !
      CALL preset_value("typeOfGeneratingProcess",         typeOfGeneratingProcess,             4 , quiet=.FALSE. )
      ! 253: ensemble
      CALL preset_value("localDefinitionNumber",           localDefinitionNumber,             253 , quiet=.FALSE. )
      ! 5  : "control and perturbed forecast products"
      CALL preset_value("typeOfProcessedData",             typeOfProcessedData,                 5 , quiet=.FALSE. )
      ! 1  : individual ensemble forecast
      CALL preset_value("productDefinitionTemplateNumber", productDefinitionTemplateNumber,     1 , quiet=.FALSE. )
      ! Note: atmospheric chemical constituents -> 41
      !       statistically processed data      -> 11
      CALL preset_value("typeOfEnsembleForecast",          typeOfEnsembleForecast,            192 , quiet=.FALSE. )
    END IF
  END SUBROUTINE preset_namelist


  !> Sets a variable to a given value, but only if current value was
  !> UNDEFINED
  !
  SUBROUTINE preset_value(name, ival, preset_val, quiet)
    CHARACTER(LEN=*), INTENT(IN)    :: name              !< name of this parameter (for screen output)
    INTEGER,          INTENT(INOUT) :: ival              !< value to be altered
    INTEGER,          INTENT(IN)    :: preset_val        !< value that shall be preset
    LOGICAL,          INTENT(IN)    :: quiet             !< LOGICAL: if .FALSE. we print some screen output
    ! local variables
    CHARACTER(len=*), PARAMETER :: routine = modname//"::preset_value"

    IF (ival == UNDEFINED) THEN
      ival = preset_val
      IF (.NOT. quiet) THEN
        CALL message(routine, "presetting namelist parameter '"//TRIM(name)//"' as "//TRIM(int2string(preset_val)))
      END IF
    ELSE
      IF (ival /= preset_val) THEN
        ! obviously, the user tried to set both: the main switch for
        ! presetting values and the local value
        CALL finish(routine, "Namelist setting of '"//TRIM(name)//"' contradicts to preset values!")
      END IF
    END IF
    
  END SUBROUTINE preset_value

END MODULE mo_gribout_nml
