!! -------------------------------------------------------------------------
!>
!! Module containing utility routines for setting GRIB2 keys.
!!
!! @author Daniel Reinert, DWD
!!
!! @par Revision History
!! Initial implementation  by  D. Reinert (2014)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
MODULE mo_grib2_util

  USE mo_impl_constants,     ONLY: MAX_CHAR_LENGTH
  USE mo_exception,          ONLY: finish, message
  USE mo_cdi,                ONLY: streamInqVlist, vlistInqVarTypeOfGeneratingProcess, vlistInqVarTsteptype, vlistInqTaxis, &
                                 & taxisInqTunit, TSTEP_CONSTANT, TSTEP_AVG, TSTEP_ACCUM, TSTEP_MAX, TSTEP_MIN, TUNIT_SECOND, &
                                 & TUNIT_MINUTE, TUNIT_HOUR, vlistDefVarIntKey, vlistDefVarProductDefinitionTemplate, &
                                 & vlistDefVarTypeOfGeneratingProcess
  USE mo_gribout_config,     ONLY: t_gribout_config
  USE mo_var_metadata_types, ONLY: t_var_metadata, CLASS_TILE, CLASS_SYNSAT, &
    &                              CLASS_TILE_LAND
  USE mo_action,             ONLY: ACTION_RESET, getActiveAction
  USE mo_util_string,        ONLY: one_of
#ifndef __NO_ICON_ATMO__
  USE mo_lnd_nwp_config,     ONLY: getNumberOfTiles, select_tile, t_tile
#endif
  ! calendar operations
  USE mtime,                 ONLY: timedelta, newTimedelta,                 &
    &                              datetime, newDatetime,                   &
    &                              deallocateTimedelta, deallocateDatetime, &
    &                              MAX_DATETIME_STR_LEN,                    &
    &                              PROLEPTIC_GREGORIAN, setCalendar,        &
    &                              MAX_DATETIME_STR_LEN,                    &
    &                              OPERATOR(-)

  IMPLICIT NONE

  PRIVATE

  ! Subroutines/Functions
  PUBLIC :: set_GRIB2_additional_keys
  PUBLIC :: set_GRIB2_ensemble_keys
  PUBLIC :: set_GRIB2_synsat_keys
  PUBLIC :: set_GRIB2_local_keys
  PUBLIC :: set_GRIB2_tile_keys
  PUBLIC :: set_GRIB2_art_keys
  PUBLIC :: set_GRIB2_timedep_keys
  PUBLIC :: set_GRIB2_timedep_local_keys

CONTAINS

  !------------------------------------------------------------------------------------------------
  !> Set additional GRIB2 keys
  !
  !
  ! Set additional GRIB2 keys. These are added to each single variable, even though
  ! adding it to the vertical or horizontal grid description may be more elegant.
  !
  SUBROUTINE set_GRIB2_additional_keys(vlistID, varID, grib_conf)
    INTEGER,                INTENT(IN) :: vlistID, varID
    TYPE(t_gribout_config), INTENT(IN) :: grib_conf
    !
    ! Local
    INTEGER :: steptype


    ! inquire steptype (needed below)
    steptype = vlistInqVarTsteptype(vlistID, varID)


    ! SECTION 1: Identification Section
    !
    ! Load correct tables
    !
    ! set tablesVersion=11
    CALL vlistDefVarIntKey(vlistID, varID, "tablesVersion", 11)
    !
    CALL vlistDefVarIntKey(vlistID, varID, "significanceOfReferenceTime",     &
      &                    grib_conf%significanceOfReferenceTime)
    CALL vlistDefVarIntKey(vlistID, varID, "productionStatusOfProcessedData", &
      &                    grib_conf%productionStatusOfProcessedData)
    CALL vlistDefVarIntKey(vlistID, varID, "typeOfProcessedData",             &
      &                    grib_conf%typeOfProcessedData)


    ! SECTION 3: Grid definition Section
    CALL vlistDefVarIntKey(vlistID, varID, "shapeOfTheEarth", 6)



    ! SECTION 4: Product definition
    CALL vlistDefVarIntKey(vlistID, varID, "backgroundProcess",               &
      &                    grib_conf%backgroundProcess)
    ! in case of lon-lat output, "1" has to be added to generatingProcessIdentifier
    CALL vlistDefVarIntKey(vlistID, varID, "generatingProcessIdentifier",     &
      &                    grib_conf%generatingProcessIdentifier)


    IF (ANY((/TSTEP_AVG,TSTEP_ACCUM,TSTEP_MAX,TSTEP_MIN/) == steptype)) THEN
      ! Always set
      !   typeOfTimeIncrement = 2 
      !   "Successive times processed have same start time of forecast, 
      !    forecast time is incremented"
      ! since this is the only type of time processing available in ICON
      CALL vlistDefVarIntKey(vlistID, varID, "typeOfTimeIncrement", 2)
    ENDIF

    IF (grib_conf%lspecialdate_invar) THEN
      ! Use special date for invariant and climatological fields
      !
      IF ( steptype == TSTEP_CONSTANT ) THEN
        ! invariant data 
        CALL vlistDefVarTypeOfGeneratingProcess(vlistID, varID, 196)
      ELSE
        CALL vlistDefVarTypeOfGeneratingProcess(vlistID, varID, grib_conf%typeOfGeneratingProcess)
      ENDIF
    ELSE
      ! no special treatment of invariant and climatological fields
      !
      CALL vlistDefVarTypeOfGeneratingProcess(vlistID, varID, grib_conf%typeOfGeneratingProcess)
    ENDIF

  END SUBROUTINE set_GRIB2_additional_keys



  !------------------------------------------------------------------------------------------------

  !>
  !! Set local GRIB2 keys (SECTION 2)
  !!
  !! Set DWD specific local keys (SECTION 2)
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2014-01-20)
  !! Modification by Daniel Reinert, DWD (2015-01-20)
  !! - extract local settings from routine set_additional_GRIB2_keys 
  !!
  SUBROUTINE set_GRIB2_local_keys(vlistID, varID, grib_conf)
    INTEGER,                INTENT(IN) :: vlistID, varID
    TYPE(t_gribout_config), INTENT(IN) :: grib_conf

  !----------------------------------------------------------------

    ! SECTION 2: Initialize local use section
    CALL vlistDefVarIntKey(vlistID, varID, "grib2LocalSectionPresent", 1)


    ! SECTION 2: DWD specific settings (local use)
    !
    ! DWD :: center    = 78
    !        subcenter = 255
    IF ((grib_conf%generatingCenter == 78) .AND. (grib_conf%generatingSubcenter == 255)) THEN

      CALL vlistDefVarIntKey(vlistID, varID, "localDefinitionNumber"  ,         &
        &                    grib_conf%localDefinitionNumber)

      CALL vlistDefVarIntKey(vlistID, varID, "localNumberOfExperiment",         &
        &                    grib_conf%localNumberOfExperiment)



      IF (grib_conf%localDefinitionNumber == 254) THEN
        !
        ! -------------------------------------------
        ! Local definition for deterministic forecast
        ! -------------------------------------------

      ELSE IF (grib_conf%localDefinitionNumber == 253) THEN
        !
        ! --------------------------------------
        ! Local definition for ensemble products
        ! --------------------------------------

        IF (grib_conf%localTypeOfEnsembleForecast /= -1)                                  &
          &   CALL vlistDefVarIntKey(vlistID, varID, "localTypeOfEnsembleForecast" ,      &
          &                          grib_conf%localTypeOfEnsembleForecast)

      END IF ! localDefinitionNumber
    END IF

  END SUBROUTINE set_GRIB2_local_keys


  !------------------------------------------------------------------------------------------------

  !>
  !! Set GRIB2 ensemble keys
  !!
  !! Set GRIB2 ensemble keys (SECTION 4)
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2014-01-20)
  !! Modification by Daniel Reinert, DWD (2015-01-20)
  !! - extract ensemble settings from routine set_additional_GRIB2_keys
  !!
  !! ATTENTION: To be called AFTER set_GRIB2_additional_keys
  !!            due to its dependency on typeOfGeneratingProcess that may 
  !!            be changed for invariant fields in set_GRIB2_additional_keys
  !!
  SUBROUTINE set_GRIB2_ensemble_keys(vlistID, varID, grib_conf)
    INTEGER,                INTENT(IN) :: vlistID, varID
    TYPE(t_gribout_config), INTENT(IN) :: grib_conf

    ! Local
    INTEGER  :: typeOfGeneratingProcess 
  !----------------------------------------------------------------

    ! get typeOfGeneratingProcess
    ! We do not make use of grib_conf%typeOfGeneratingProcess, since 
    ! typeOfGeneratingProcess is modified for invariant fields.
    typeOfGeneratingProcess = vlistInqVarTypeOfGeneratingProcess(vlistID, varID)


    IF (typeOfGeneratingProcess == 4) THEN  ! Ensemble forecast

      ! SECTION 4: Product definition Section

      IF (grib_conf%typeOfEnsembleForecast /= -1)                                       &
        &   CALL vlistDefVarIntKey(vlistID, varID, "typeOfEnsembleForecast" ,           &
        &                          grib_conf%typeOfEnsembleForecast)
      IF (grib_conf%numberOfForecastsInEnsemble /= -1)                                  &
        &   CALL vlistDefVarIntKey(vlistID, varID, "numberOfForecastsInEnsemble" ,      &
        &                          grib_conf%numberOfForecastsInEnsemble)
      IF (grib_conf%perturbationNumber /= -1)                                           &
        &   CALL vlistDefVarIntKey(vlistID, varID, "perturbationNumber" ,               &
        &                          grib_conf%perturbationNumber)
    END IF ! typeOfGeneratingProcess

  END SUBROUTINE set_GRIB2_ensemble_keys


  !>
  !! Set synsat-specific keys
  !!
  !! Set GRIB2 keys which are specific to synthetic satellite products.
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2015-06-10)
  !!
  SUBROUTINE set_GRIB2_synsat_keys (vlistID, varID, info)
    
    INTEGER,                INTENT(IN) :: vlistID, varID
    TYPE (t_var_metadata),  INTENT(IN) :: info
    
    ! ----------------------------------------------------------------
    
    ! Skip inapplicable fields
    IF ( info%var_class /= CLASS_SYNSAT ) RETURN
    
    ! change product definition template
    CALL vlistDefVarProductDefinitionTemplate(vlistID, varID, 32)
    
  END SUBROUTINE set_GRIB2_synsat_keys


  !>
  !! Set tile-specific keys
  !!
  !! Set GRIB2 keys which are specific to tile-variables.
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2014-01-20)
  !!
  SUBROUTINE set_GRIB2_tile_keys (vlistID, varID, info, i_lctype)

    INTEGER,                INTENT(IN) :: vlistID, varID
    TYPE (t_var_metadata),  INTENT(IN) :: info
    INTEGER,                INTENT(IN) :: i_lctype  !< Tile classification

#ifndef __NO_ICON_ATMO__
    ! local
    INTEGER      :: typeOfGeneratingProcess
    INTEGER      :: productDefinitionTemplate        ! Tile template number 
    TYPE(t_tile) :: tile

  !----------------------------------------------------------------

    ! Skip inapplicable fields
    IF ( ALL((/CLASS_TILE,CLASS_TILE_LAND/) /= info%var_class) ) RETURN

    typeOfGeneratingProcess = vlistInqVarTypeOfGeneratingProcess(vlistID, varID)

    ! change product definition template
    !
    IF (typeOfGeneratingProcess == 4) THEN
      ! ensemble
      productDefinitionTemplate = 40456     ! FIXME (temporary)
    ELSE
      ! deterministic
      productDefinitionTemplate = 40455     ! FIXME (temporary)
    ENDIF
    CALL vlistDefVarProductDefinitionTemplate(vlistID, varID, productDefinitionTemplate)

    ! Set tile classification
    CALL vlistDefVarIntKey(vlistID, varID, "tileClassification" , i_lctype)

    ! Set total number of tile/attribute pairs
    CALL vlistDefVarIntKey(vlistID, varID, "totalNumberOfTileAttributePairs" , info%maxcontained)

    ! Set number of used tiles
    CALL vlistDefVarIntKey(vlistID, varID, "numberOfTiles" , &
      &                    getNumberOfTiles(class=info%var_class))

    ! get the following attributes:
    ! - tileIndex
    ! - numberOfTileAttributes
    !
    ! Select tile info orresponding to index info%ncontained
    tile = select_tile(info%ncontained)

    ! Set tile index
    CALL vlistDefVarIntKey(vlistID, varID, "tileIndex" , tile%GRIB2_tile%tileIndex)

    ! Set total number of tile attributes for given tile
    CALL vlistDefVarIntKey(vlistID, varID, "numberOfTileAttributes" , tile%GRIB2_tile%numberOfTileAttributes)

    ! Set tile attribute
    CALL vlistDefVarIntKey(vlistID, varID, "tileAttribute" , tile%GRIB2_att%tileAttribute)
#endif

  END SUBROUTINE set_GRIB2_tile_keys


  !>
  !! Set ART-specific keys
  !!
  !! Set GRIB2 keys which are specific to ART-variables.
  !!
  !! @par Revision History
  !! Initial revision by Jochen Foerstner, DWD (2015-09-07)
  !!
  SUBROUTINE set_GRIB2_art_keys(vlistID, varID, info)

    INTEGER,                INTENT(IN) :: vlistID, varID
    TYPE (t_var_metadata),  INTENT(IN) :: info

#ifdef __ICON_ART

    ! local
    INTEGER      :: productDefinitionTemplate
    INTEGER      :: numberOfDistributionFunctionParameter
    INTEGER, ALLOCATABLE ::  &
      &  scaledValueOfDistributionFunctionParameter(:),  &
      &  scaleFactorOfDistributionFunctionParameter(:)

  !----------------------------------------------------------------
    CHARACTER(LEN=100) :: message_text
    WRITE(message_text,'(a,i4,a,i4)') 'vlistID = ', vlistID, '  varID = ', varID
    CALL message(' ==> set_GRIB2_art_keys :',TRIM(message_text),0,5,.TRUE.)
    CALL message(' ==> set_GRIB2_art_keys :','tracer_class = '//TRIM(info%tracer%tracer_class),0,5,.TRUE.)
 
    ! change product definition template
    SELECT CASE(TRIM(info%tracer%tracer_class))
    CASE ('volcash')
      CALL message(' ==> set_GRIB2_art_keys :','volcash --> PDT=57',0,5,.TRUE.)
      productDefinitionTemplate = 57
      numberOfDistributionFunctionParameter = 1
    CASE ('volcash_diag_mc', 'volcash_diag_mc_max', 'volcash_diag_mc_vi', 'volcash_diag_hml')
      CALL message(' ==> set_GRIB2_art_keys :','volcash_diag_... --> PDT=40',0,5,.TRUE.)
      productDefinitionTemplate = 40
      numberOfDistributionFunctionParameter = 0
    CASE ('dust', 'dust_number')
      CALL message(' ==> set_GRIB2_art_keys :','dust[_number] --> PDT=57',0,5,.TRUE.)
      productDefinitionTemplate = 57
      numberOfDistributionFunctionParameter = 2
    CASE ('dust_diag_tau')
      CALL message(' ==> set_GRIB2_art_keys :','dust_diag_tau --> PDT=48',0,5,.TRUE.)
      productDefinitionTemplate = 48
      numberOfDistributionFunctionParameter = 0
    CASE DEFAULT
      ! skip inapplicable fields
      RETURN
    END SELECT

    ! set product definition template
    CALL vlistDefVarProductDefinitionTemplate(vlistID, varID, productDefinitionTemplate)

    IF ( numberOfDistributionFunctionParameter /= 0 ) THEN
      ALLOCATE( scaledValueOfDistributionFunctionParameter(numberOfDistributionFunctionParameter),    &
        &       scaleFactorOfDistributionFunctionParameter(numberOfDistributionFunctionParameter) )
    END IF
    
    SELECT CASE(TRIM(info%tracer%tracer_class))

    CASE ('volcash')
      CALL vlistDefVarIntKey(vlistID, varID, "typeOfDistributionFunction", 1)
      CALL vlistDefVarIntKey(vlistID, varID, "numberOfModeOfDistribution", 6)
      CALL vlistDefVarIntKey(vlistID, varID, "numberOfDistributionFunctionParameter",  &
        &                    numberOfDistributionFunctionParameter)
      CALL vlistDefVarIntKey(vlistID, varID, "constituentType", 62025)
      CALL vlistDefVarIntKey(vlistID, varID, "modeNumber", info%tracer%mode_number)
      scaledValueOfDistributionFunctionParameter(1) = info%tracer%diameter
      scaleFactorOfDistributionFunctionParameter(1) = 6
      CALL vlistDefVarIntArrayKey(vlistID, varID, "scaledValueOfDistributionFunctionParameter",     &
        &   numberOfDistributionFunctionParameter, scaledValueOfDistributionFunctionParameter(:))
      CALL vlistDefVarIntArrayKey(vlistID, varID, "scaleFactorOfDistributionFunctionParameter",    &
        &   numberOfDistributionFunctionParameter, scaleFactorOfDistributionFunctionParameter(:))

    CASE ('volcash_diag_mc', 'volcash_diag_mc_max')
      CALL vlistDefVarIntKey(vlistID, varID, "constituentType", 62025)
      CALL vlistDefVarIntKey(vlistID, varID, "decimalScaleFactor", 9)
      
    CASE ('volcash_diag_mc_vi')
      CALL vlistDefVarIntKey(vlistID, varID, "constituentType", 62025)
      CALL vlistDefVarIntKey(vlistID, varID, "decimalScaleFactor", 3)
      
    CASE ('volcash_diag_hml')
      CALL vlistDefVarIntKey(vlistID, varID, "constituentType", 62025)
! JF:       CALL vlistDefVarCharKey(vlistID, varID, "typeOfLevel", "entireAtmosphere")
      
    CASE ('dust', 'dust_number')
      CALL vlistDefVarIntKey(vlistID, varID, "typeOfDistributionFunction", 7)
      CALL vlistDefVarIntKey(vlistID, varID, "numberOfModeOfDistribution", 3)
      CALL vlistDefVarIntKey(vlistID, varID, "numberOfDistributionFunctionParameter",  &
        &                    numberOfDistributionFunctionParameter)
      CALL vlistDefVarIntKey(vlistID, varID, "constituentType", 62001)
      CALL vlistDefVarIntKey(vlistID, varID, "modeNumber", info%tracer%mode_number)
      scaledValueOfDistributionFunctionParameter(1) = info%tracer%variance
      scaleFactorOfDistributionFunctionParameter(1) = 1
      scaledValueOfDistributionFunctionParameter(2) = 2650
      scaleFactorOfDistributionFunctionParameter(2) = 0
      CALL vlistDefVarIntArrayKey(vlistID, varID, "scaledValueOfDistributionFunctionParameter",     &
        &   numberOfDistributionFunctionParameter, scaledValueOfDistributionFunctionParameter(:))
      CALL vlistDefVarIntArrayKey(vlistID, varID, "scaleFactorOfDistributionFunctionParameter",    &
        &   numberOfDistributionFunctionParameter, scaleFactorOfDistributionFunctionParameter(:))
      IF (TRIM(info%tracer%tracer_class) == 'dust') THEN
        CALL vlistDefVarIntKey(vlistID, varID, "decimalScaleFactor", 9)
      END IF
      
    CASE ('dust_diag_tau')
      CALL vlistDefVarIntKey(vlistID, varID, "aerosolType", 62001)
      CALL vlistDefVarIntKey(vlistID, varID, "typeOfSizeInterval", 192)
      CALL vlistDefVarIntKey(vlistID, varID, "typeOfWavelengthInterval", 11)
      CALL vlistDefVarIntKey(vlistID, varID, "scaledValueOfFirstWavelength", info%tracer%tau_wavelength)
      CALL vlistDefVarIntKey(vlistID, varID, "scaleFactorOfFirstWavelength", 9)
      
    END SELECT

    IF ( numberOfDistributionFunctionParameter /= 0 ) THEN
      DEALLOCATE( scaledValueOfDistributionFunctionParameter,    &
        &         scaleFactorOfDistributionFunctionParameter )
    END IF
    
#endif

  END SUBROUTINE set_GRIB2_art_keys


  !------------------------------------------------------------------------------------------------
  !> Set additional, time-dependent GRIB2 keys
  !!  
  !!  This subroutine sets all GRIB2 keys that may change during the
  !!  simulation. Currently this is the case for the start and end time 
  !!  of the statistical process interval.
  !!
  !!  Description of how the GRIB2 keys 'forecastTime' and 'lengthOfTimeRange' are computed.
  !!
  !!  |===================*======================|================================> time axis
  !!  ^start_time         |                      ^current_time (output)
  !!                      |                      |
  !!                      ^EventLastTriggerTime  |
  !!                      |                      |
  !!                      |                      |
  !!  <-------------------><--------------------->
  !!     forecastTime         lengthOfTimeRange
  !!
  !!  @author D. Reinert, F. Prill (DWD)
  !!
  !!  CAVEATs
  !!
  !!  - we implicitly assume that actions are ordered according to increasing forecast time.
  !!  - we implicitly assume that the statistical process time range is always smaller than one month
  !!
  SUBROUTINE set_GRIB2_timedep_keys(streamID, varID, info, start_date, cur_date)
    INTEGER                             ,INTENT(IN) :: streamID, varID
    TYPE (t_var_metadata)               ,INTENT(IN) :: info
    CHARACTER(LEN=MAX_DATETIME_STR_LEN) ,INTENT(IN) :: start_date, cur_date
    ! local variables
    TYPE(timedelta), POINTER      :: mtime_lengthOfTimeRange       ! length of time range (mtime format)
    INTEGER                       :: ilengthOfTimeRange_secs       ! length of time range in seconds
    INTEGER                       :: ilengthOfTimeRange            ! length of time range in appropriate time units
    INTEGER                       :: taxis_tunit                   ! time axis units
    INTEGER                       :: vlistID, taxisID
    TYPE(timedelta), POINTER      :: forecast_delta                ! forecast time (mtime format)
    INTEGER                       :: forecast_secs                 ! forecast time in seconds
    INTEGER                       :: forecast_time                 ! forecast time in appropriate time units
    TYPE(datetime),  POINTER      :: mtime_start                   ! model start (initialization) time
    TYPE(datetime),  POINTER      :: mtime_cur                     ! model current time (rounded)
    TYPE(datetime)                :: statProc_startDateTime        ! statistical process starting DateTime
    INTEGER                       :: var_actionId                  ! action from which metainfo is used
    CHARACTER(len=*), PARAMETER   :: routine = 'set_GRIB2_timedep_keys'

    ! special fields for which time-dependent metainfos should be set even though they are not of 
    ! steptype TSTEP_MAX or TSTEP_MIN. These fields are special in the sense that averaging is not 
    ! performed over the entire model run but over only some intervals.
    CHARACTER(LEN=8) :: ana_avg_vars(5) = (/"u_avg   ", "v_avg   ", "pres_avg", "temp_avg", "qv_avg  "/)


    CALL setCalendar(PROLEPTIC_GREGORIAN)

    !---------------------------------------------------------
    ! Set time-dependent metainfo
    !---------------------------------------------------------
    !
    ! Skip inapplicable fields
    ! Currently all TSTEP_AVG and TSTEP_ACC fields are skipped, except for special ones 
    ! listed in ana_avg_vars
    IF ((ALL((/TSTEP_MAX, TSTEP_MIN/) /= info%isteptype)) .AND. &
      & (one_of(TRIM(info%name),ana_avg_vars) == -1) ) RETURN


    ! get vlistID. Note that the stream-internal vlistID must be used. 
    ! It is obtained via streamInqVlist(streamID)
    vlistID     = streamInqVlist(streamID)
    taxisID     = vlistInqTaxis(vlistID)
    taxis_tunit = taxisInqTunit(taxisID)

    ! get current DateTime (rounded)
    mtime_cur   => newDatetime(TRIM(cur_date))
    ! get model start date
    mtime_start => newDatetime(TRIM(start_date))

    IF (info%action_list%n_actions > 0) THEN

      ! more than one RESET action may be defined for a single variable.
      ! get ID of currently active RESET action
      var_actionId = getActiveAction(info, ACTION_RESET, mtime_cur)

      IF (var_actionId == -1) THEN
        write(0,*) 'set_timedependent_GRIB2_keys: no active action of type ACTION_RESET found. '//&
          &             'lengthOfTimeRange may not be set correctly'
        CALL finish (routine, 'Illegal actionId')
      ENDIF

      ! get latest (intended) triggering time, which is equivalent to 
      ! the statistical process starting time
      statProc_startDateTime = info%action_list%action(var_actionId)%EventLastTriggerDate


      ! get time interval, over which statistical process has been performed
      ! It is the time difference between the current time (rounded) mtime_cur and 
      ! the last time the nullify-action took place (rounded) statProc_startDateTime.
      mtime_lengthOfTimeRange  => newTimedelta("P01D")  ! init
      !
      ! mtime_lengthOfTimeRange = mtime_cur - statProc_startDateTime
      mtime_lengthOfTimeRange = mtime_cur - statProc_startDateTime

      ! time interval over which statistical process has been performed (in secs)    
      ilengthOfTimeRange_secs = 86400 *INT(mtime_lengthOfTimeRange%day)    &
           &                  + 3600  *INT(mtime_lengthOfTimeRange%hour)   &
           &                  + 60    *INT(mtime_lengthOfTimeRange%minute) &
           &                  +        INT(mtime_lengthOfTimeRange%second) 
           

      ! cleanup
      CALL deallocateTimedelta(mtime_lengthOfTimeRange)
    ELSE
      ilengthOfTimeRange_secs = 0
    END IF


    ! get forecast_time: forecast_time = statProc_startDateTime - model_startDateTime
    ! Note that for statistical quantities, the forecast time is the time elapsed between the 
    ! model start time and the start time of the statistical process
    forecast_delta => newTimedelta("P01D")
    forecast_delta = statProc_startDateTime - mtime_start

    ! forecast time in seconds
    forecast_secs =    forecast_delta%second    +   &
      &             60*(forecast_delta%minute   +   & 
      &                 60*(forecast_delta%hour +   &
      &                       24*forecast_delta%day))


    SELECT CASE (taxis_tunit)
    CASE (TUNIT_SECOND)
       ilengthOfTimeRange          = ilengthOfTimeRange_secs
       forecast_time               = forecast_secs
    CASE (TUNIT_MINUTE)
       ilengthOfTimeRange          = INT(ilengthOfTimeRange_secs/60)
       forecast_time               = INT(forecast_secs/60)
    CASE (TUNIT_HOUR)
       ilengthOfTimeRange          = INT(ilengthOfTimeRange_secs/3600)
       forecast_time               = INT(forecast_secs/3600)
    CASE DEFAULT
    END SELECT

    ! set forecast time: statProc_startDateTime - model_startDateTime
    CALL vlistDefVarIntKey(vlistID, varID, "forecastTime", forecast_time) 

    !
    ! set length of time range: current time - statProc_startDateTime
    CALL vlistDefVarIntKey(vlistID, varID, "lengthOfTimeRange", ilengthOfTimeRange)
    ! Note that if one of the statistics templates 4.8 or 4.11 is selected, the time unit 
    ! (GRIB2 key "indicatorOfUnitForTimeRange") is set automatically by CDI.
    ! It is always set identical to "indicatorOfUnitOFTimeRange"


    ! cleanup
    CALL deallocateDatetime(mtime_start)
    CALL deallocateDatetime(mtime_cur)
    CALL deallocateTimedelta(forecast_delta)

  END SUBROUTINE set_GRIB2_timedep_keys


  !------------------------------------------------------------------------------------------------
  !>
  !! Set time-dependent local GRIB2 keys (SECTION 2)
  !!
  !! Set DWD specific time-dependent local keys (SECTION 2)
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2015-02-06)
  !!
  SUBROUTINE set_GRIB2_timedep_local_keys(streamID, varID, grib_conf)
    INTEGER,                INTENT(IN) :: streamID, varID
    TYPE(t_gribout_config), INTENT(IN) :: grib_conf

    ! Local
    CHARACTER(len=MAX_CHAR_LENGTH)     :: ydate, ytime
    INTEGER :: cent, year, month, day    ! date
    INTEGER :: hour, minute, second      ! time
    INTEGER :: vlistID
  !----------------------------------------------------------------


    ! get vlistID. Note that the stream-internal vlistID must be used. 
    ! It is obtained via streamInqVlist(streamID)
    vlistID = streamInqVlist(streamID)

    ! SECTION 2: DWD specific settings (local use)
    !
    ! DWD :: center    = 78
    !        subcenter = 255
    IF ((grib_conf%generatingCenter == 78) .AND. (grib_conf%generatingSubcenter == 255)) THEN

      IF (grib_conf%ldate_grib_act) THEN
        ! get date and time
        ! ydate : ccyymmdd, ytime : hhmmss.sss
        CALL date_and_time(ydate,ytime)
        READ(ydate,'(4i2)') cent, year, month, day
        READ(ytime,'(3i2)') hour, minute, second
      ELSE ! set date to "01010101" (for better comparability of GRIB files)
        cent  = 1
        year  = 1
        month = 1
        hour  = 1
        minute= 1
        second= 1
      ENDIF


      CALL vlistDefVarIntKey(vlistID, varID, "localCreationDateYear"  , 100*cent+year)
      CALL vlistDefVarIntKey(vlistID, varID, "localCreationDateMonth" , month)
      CALL vlistDefVarIntKey(vlistID, varID, "localCreationDateDay"   , day)
      CALL vlistDefVarIntKey(vlistID, varID, "localCreationDateHour"  , hour)
      CALL vlistDefVarIntKey(vlistID, varID, "localCreationDateMinute", minute)
      CALL vlistDefVarIntKey(vlistID, varID, "localCreationDateSecond", second)
      ! CALL vlistDefVarIntKey(vlistID, varID, "localValidityDateYear"  , 2013)

    END IF

  END SUBROUTINE set_GRIB2_timedep_local_keys


END MODULE mo_grib2_util
