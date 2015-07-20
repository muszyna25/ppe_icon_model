!>
!! @brief Master namelist.
!!        
!! @par Revision History
!! Created by Rene Redler (2011-03-22)
!! Major revision by Luis Kornblueh (2015-04-15)
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

  USE mo_exception,      ONLY: finish, message, message_text
  USE mo_io_units,       ONLY: filename_max, nnml
  USE mo_namelist,       ONLY: open_nml, position_nml, POSITIONED
  USE mo_util_string,    ONLY: t_keyword_list, associate_keyword, with_keywords, tolower
  USE mo_nml_annotate,   ONLY: temp_defaults, temp_settings
  USE mo_mpi,            ONLY: my_process_is_stdio, get_my_global_mpi_communicator
  USE mtime,             ONLY: max_calendar_str_len, setCalendar,  calendarToString,         &
       &                       proleptic_gregorian, year_of_365_days,  year_of_360_days,     &
       &                       max_datetime_str_len, datetimeToString,                       &
       &                       datetime, newDatetime, deallocateDatetime,                    &
       &                       max_timedelta_str_len, timedeltaToString,                     &
       &                       OPERATOR(+), OPERATOR(<)
  USE mo_master_config,  ONLY: master_component_models, addModel, noOfModels, maxNoOfModels, &
       &                       setRestart, isREstart, setModelBaseDir,                       &
       &                       setExpRefdate,                                                &
       &                       setExpStartdate, setExpStopdate,                              &
       &                       setCheckpointTimeInterval,  setRestartTimeInterval,           &
       &                       tc_exp_refdate, tc_exp_startdate, tc_exp_stopdate,            &
       &                       tc_dt_checkpoint, tc_dt_restart

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

    LOGICAL :: lRestart = .FALSE.
    CHARACTER(len=filename_max) :: modelBaseDir = ''
    
    CHARACTER(len=132)          :: modelName = ''
    CHARACTER(len=filename_max) :: modelNamelistFilename = ''
    
    INTEGER :: modelType 
    INTEGER :: modelMinRank
    INTEGER :: modelMaxRank
    INTEGER :: modelIncRank
    
    CHARACTER(len=max_calendar_str_len) :: calendar = ''
    
    CHARACTER(len=max_datetime_str_len) :: experimentReferenceDate = ''   
    CHARACTER(len=max_datetime_str_len) :: experimentStartDate = ''
    CHARACTER(len=max_datetime_str_len) :: experimentStopDate = ''
    
    CHARACTER(len=max_timedelta_str_len) :: checkpointTimeIntval = ''
    CHARACTER(len=max_timedelta_str_len) :: restartTimeIntval = ''
    
    NAMELIST /master_nml/              &
         &    lRestart,                &
         &    modelBaseDir
    
    NAMELIST /master_time_control_nml/ &
         &    calendar,                &
         &    experimentReferenceDate, &   
         &    experimentStartDate,     &
         &    experimentStopDate,      &
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

    CHARACTER(len=max_datetime_str_len) :: dstring
    
    CHARACTER(len=*), PARAMETER :: routine = 'mo_master_nml:read_master_namelist'
    
    TYPE (t_keyword_list), POINTER :: keywords         => NULL()
    
    ! Read  master_nml (done so far by all MPI processes)
    
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

    ! set calendar (singleton, so do not change later!)
    
    SELECT CASE (toLower(calendar))
    CASE ('proleptic gregorian')
      icalendar  = proleptic_gregorian
    CASE ('365 day year')  
      icalendar = year_of_365_days
    CASE ('360 day year')  
      icalendar = year_of_360_days
    CASE default
      icalendar  = proleptic_gregorian
      CALL message('','No calendar selected! Use default proleptic Gregorian.')
    END SELECT

    CALL setCalendar(icalendar)
      
    ! take care that the reference date and start date strings are set 
    
    IF (experimentReferenceDate /= '') THEN
      CALL setExpRefdate(experimentReferenceDate)
      IF (experimentStartDate /= '') THEN
        CALL setExpStartdate(experimentStartDate)
      ELSE
        CALL setExpStartdate(experimentReferenceDate)
      ENDIF
    ELSE
      IF (experimentStartDate /= '') THEN
        CALL setExpRefdate(experimentStartDate)
        CALL setExpStartdate(experimentStartDate)
      ENDIF
    ENDIF

    IF (experimentStopDate /= '') THEN
      CALL setExpStopdate(experimentStopDate)
    ENDIF

    ! next check for given restart, checkpoint interval

    IF (checkpointTimeIntval /= '') THEN
      CALL setCheckpointTimeInterval(checkpointTimeIntval)
      IF (restartTimeIntval /= '') THEN
        CALL setRestartTimeInterval(restartTimeIntval)
      ELSE
        CALL setRestartTimeInterval(checkpointTimeIntval)
      ENDIF
    ELSE
      IF (restartTimeIntval /= '') THEN
        CALL setRestartTimeInterval(restartTimeIntval)
        CALL setCheckpointTimeInterval(restartTimeIntval)        
      ELSE
        CALL message('','No restart and checkpoint time interval given: depend on restart file.')
      ENDIF
    ENDIF



    ! inform about time setup

    call message('','')
    
    CALL calendarToString(dstring)
    CALL message('','Calendar: '//TRIM(dstring))

    call message('','')
    
    IF (ASSOCIATED(tc_exp_refdate)) THEN
      CALL datetimeToString(tc_exp_refdate, dstring)
      WRITE(message_text,'(a,a)') 'Experiment reference date: ', dstring
      CALL message('',message_text)
    ENDIF

    IF (ASSOCIATED(tc_exp_startdate)) THEN    
      CALL datetimeToString(tc_exp_startdate, dstring)
      WRITE(message_text,'(a,a)') 'Experiment start date    : ', dstring
      CALL message('',message_text)
    ENDIF

    IF (ASSOCIATED(tc_exp_stopdate)) THEN
      CALL datetimeToString(tc_exp_stopdate, dstring)
      WRITE(message_text,'(a,a)') 'Experiment stop date     : ', dstring
      CALL message('',message_text)
    ENDIF

    call message('','')
    
    IF (ASSOCIATED(tc_dt_restart)) THEN
      CALL timedeltaToString(tc_dt_restart, dstring)
      WRITE(message_text,'(a,a)') 'Restart interval         : ', dstring
      CALL message('',message_text)
    ENDIF

    IF (ASSOCIATED(tc_dt_checkpoint)) THEN
      CALL timedeltaToString(tc_dt_checkpoint, dstring)
      WRITE(message_text,'(a,a)') 'Checkpointing interval   : ', dstring
      CALL message('',message_text)
    ENDIF

    call message('','')
    
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
