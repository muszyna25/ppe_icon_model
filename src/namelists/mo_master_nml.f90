!>
!! @brief Master namelist.
!!        
!! @par Revision History
!! Created by Rene Redler (2011-03-22)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_master_nml

  USE mo_exception,      ONLY: finish
  USE mo_io_units,       ONLY: filename_max, nnml
  USE mo_namelist,       ONLY: open_nml, position_nml, POSITIONED
  USE mo_util_string,    ONLY: t_keyword_list, associate_keyword, with_keywords, tolower
  USE mo_nml_annotate,   ONLY: temp_defaults, temp_settings
  USE mo_mpi,            ONLY: my_process_is_stdio, get_my_global_mpi_communicator
  USE mtime,             ONLY: max_calendar_str_len, max_datetime_str_len, max_timedelta_str_len, &
       &                       proleptic_gregorian, year_of_365_days,  year_of_360_days,          &
       &                       setCalendar
  USE mo_master_config,  ONLY: master_component_models, addModel, noOfModels, maxNoOfModels, &
       &                       setRestart, setModelBaseDir,                                  &
       &                       setExpRefdate,                                                &
       &                       setExpStartdate, setExpStopdate,                              &
       &                       setStartdate, setStopdate,                                    &
       &                       setCheckpointTimeInterval,  setRestartTimeInterval


  
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC :: read_master_namelist
    
CONTAINS

  !>
  !! Initialization of variables that contain general information
  !! about the coupled model run. The configuration is read from
  !! namelist 'master_nml'.
  !!
  !! @par Revision History
  !!
  INTEGER FUNCTION read_master_namelist(namelist_filename)
    
    CHARACTER(len=*), INTENT(in) :: namelist_filename

    ! Local variables
    
    ! Namelist variables

    LOGICAL :: lRestart
    CHARACTER(len=filename_max) :: modelBaseDir
    
    CHARACTER(len=132)          :: modelName
    CHARACTER(len=filename_max) :: modelNamelistFilename
    
    INTEGER :: modelType
    INTEGER :: modelMinRank
    INTEGER :: modelMaxRank
    INTEGER :: modelIncRank
    
    CHARACTER(len=max_calendar_str_len) :: calendar
    
    CHARACTER(len=max_datetime_str_len) :: experimentReferenceDate   
    CHARACTER(len=max_datetime_str_len) :: experimentStartDate
    CHARACTER(len=max_datetime_str_len) :: experimentStopDate
    CHARACTER(len=max_datetime_str_len) :: startdate
    CHARACTER(len=max_datetime_str_len) :: stopdate
    
    CHARACTER(len=max_timedelta_str_len) :: checkpointTimeIntval
    CHARACTER(len=max_timedelta_str_len) :: restartTimeIntval
    
    NAMELIST /master_nml/              &
         &    lRestart,                &
         &    modelBaseDir
    
    NAMELIST /master_time_control_nml/ &
         &    calendar,                &
         &    experimentReferenceDate, &   
         &    experimentStartDate,     &
         &    experimentStopDate,      &
         &    startdate,               &
         &    stopdate,                &
         &    checkpointTimeIntval,    &
         &    restartTimeIntval        
    
    NAMELIST /master_model_nml/        &
         &    modelName,               &
         &    modelNamelistFilename,   &
         &    modelType,               &
         &    modelMinRank,            &
         &    modelMaxRank,            &
         &    modelIncRank               

    INTEGER :: icalendar
    INTEGER :: istat
    INTEGER :: iunit
    LOGICAL :: lrewind
    
    CHARACTER(len=*), PARAMETER :: routine = 'mo_master_nml:read_master_namelist'
    
    TYPE (t_keyword_list), POINTER :: keywords         => NULL()
!    TYPE (t_keyword_list), POINTER :: keywords_restart => NULL()
    
    ! Read  master_nml (done so far by all MPI processes)
    
    lRestart     = .FALSE.
    modelBaseDir = ''
    
    istat = 0
    CALL open_nml(namelist_filename, lwarn=.TRUE., istat=istat)
    IF (istat /= 0) THEN
      read_master_namelist = -1
      RETURN
    ENDIF
    
    CALL position_nml('master_nml', STATUS=istat)
    IF (istat == POSITIONED) THEN
      READ (nnml, master_nml)
    ENDIF

    ! save namelist variables in configuration

    CALL setRestart(lRestart)
    CALL setModelBaseDir(modelBaseDir)
    
    CALL position_nml('master_time_control_nml', STATUS=istat)
    IF (istat == POSITIONED) THEN
      READ (nnml, master_time_control_nml)
    ENDIF

    SELECT CASE (toLower(calendar))
    CASE ('proleptic gregorian')
      icalendar  = proleptic_gregorian
    CASE ('365 day year')  
      icalendar = year_of_365_days
    CASE ('360 day year')  
      icalendar = year_of_360_days
    CASE default
      CALL finish(routine,'Selected calendar unknown.')
    END SELECT
    CALL setCalendar(icalendar)
    CALL setExpRefdate(experimentReferenceDate)
    CALL setExpStartdate(experimentStartDate)
    CALL setExpStopdate(experimentStopDate)
    CALL setStartdate(startdate)
    CALL setStopdate(stopdate)
    CALL setCheckpointTimeInterval(checkpointTimeIntval)
    CALL setRestartTimeInterval(restartTimeIntval)

    ! Read  master namelist (done so far by all MPI processes)
    
    ! for positioning to the first entry of the namelist 
    lrewind = .TRUE.

    DO
      CALL position_nml('master_model_nml', lrewind=lrewind, status=istat)
      IF ( istat /= POSITIONED ) EXIT

      IF (noOfModels() >= maxNoOfModels) THEN
        CALL finish(routine, 'no of models >= max no of models')
      ENDIF
      
      ! change to be able to allow fetching the next namelist entry
      lrewind = .FALSE.
      
      ! default values

      modelName             = ''
      modelNamelistFilename = ''
      modelType             = -1
      modelMinRank          = 0
      modelMaxRank          = -1 
      modelIncRank          = 1
      
      IF (my_process_is_stdio()) THEN
        iunit = temp_defaults()
        WRITE(iunit, master_model_nml)  ! write defaults to temporary text file
      END IF
      
      READ (nnml, master_model_nml)     ! overwrite default settings

      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, master_model_nml)  ! write settings to temporary text file
      END IF

      CALL addModel()
      master_component_models(noOfModels())%model_name = modelName

      CALL associate_keyword("<path>", TRIM(modelBaseDir), keywords)

      master_component_models(noOfModels())%model_namelist_filename = TRIM(with_keywords(keywords, modelNamelistFilename))

      master_component_models(noOfModels())%model_type = modelType
      master_component_models(noOfModels())%model_min_rank = modelMinRank
      master_component_models(noOfModels())%model_max_rank = modelMaxRank
      master_component_models(noOfModels())%model_inc_rank = modelIncRank

    ENDDO
      
    CLOSE (nnml, IOSTAT=istat)

    read_master_namelist = 0

  END FUNCTION read_master_namelist

END MODULE mo_master_nml
