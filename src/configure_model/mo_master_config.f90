!>
!! @brief Master configuration.
!!        
!! @par Revision History
!! Created by Luis Kornblueh (2015-04-08)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_master_config

  USE mtime,       ONLY: datetime, timedelta, newDatetime, newTimedelta, &
    &                    MAX_DATETIME_STR_LEN, MAX_CALENDAR_STR_LEN,     &
    &                    MAX_TIMEDELTA_STR_LEN
  USE mo_io_units, ONLY: filename_max
  
  IMPLICIT NONE
  
  PRIVATE

  PUBLIC :: addModel, noOfModels, maxNoOfModels
  PUBLIC :: setInstitution
  PUBLIC :: setModelBaseDir, getModelBaseDir
  PUBLIC :: setRestart, setRestartWriteLast, isRestart
  PUBLIC :: setExpRefdate
  PUBLIC :: setExpStartdate, setExpStopdate
  PUBLIC :: setStartdate, setStopdate
  PUBLIC :: setCheckpointTimeInterval
  PUBLIC :: setRestartTimeInterval
  PUBLIC :: setCurrentdate
  PUBLIC :: master_component_models
  PUBLIC :: tc_exp_refdate
  PUBLIC :: tc_exp_startdate
  PUBLIC :: tc_exp_stopdate
  PUBLIC :: tc_startdate
  PUBLIC :: tc_stopdate
  PUBLIC :: tc_current_date
  PUBLIC :: tc_dt_checkpoint
  PUBLIC :: tc_dt_restart
  PUBLIC :: lrestart_write_last
  PUBLIC :: calendar
  PUBLIC :: experimentReferenceDate, experimentStartDate, experimentStopDate
  PUBLIC :: checkpointTimeIntval, restartTimeIntval
  
  ! component model configuration
  !_______________________________________________________________________________________________
  !
  !> Holds name, type, and position of processes in the global communicator
  !  of a component model
  
  TYPE t_master_component_model_config
    CHARACTER(len=132)          :: model_name
    CHARACTER(len=filename_max) :: model_namelist_filename
    INTEGER :: model_type
    INTEGER :: model_min_rank
    INTEGER :: model_max_rank
    INTEGER :: model_inc_rank
  END TYPE t_master_component_model_config

  INTEGER, PARAMETER :: maxNoOfModels = 16
  INTEGER :: no_of_models = 0
  TYPE(t_master_component_model_config) :: master_component_models(maxNoOfModels)

  ! defaults to DWD to be on the safe side for operational needs  
  CHARACTER(len=256), PROTECTED :: institution = '' 
  
  ! whole experiment and sinlge run time information
  !_______________________________________________________________________________________________
  !
  ! experiment

  TYPE(datetime), POINTER, PROTECTED :: tc_exp_refdate => NULL()
                                                          
  TYPE(datetime), POINTER, PROTECTED :: tc_exp_startdate => NULL()
  TYPE(datetime), POINTER, PROTECTED :: tc_exp_stopdate => NULL()

  ! single run 
  
  TYPE(datetime), POINTER, PROTECTED :: tc_startdate => NULL()
  TYPE(datetime), POINTER, PROTECTED :: tc_stopdate => NULL()

  ! current model date

  TYPE(datetime), POINTER, PROTECTED :: tc_current_date => NULL()
  
  ! checkpoint and restart time interval (needs to be equal for all component models) 
  
  TYPE(timedelta), POINTER, PROTECTED :: tc_dt_checkpoint => NULL()
  TYPE(timedelta), POINTER, PROTECTED :: tc_dt_restart => NULL()
  
  ! restart flag, required for consistent handling in all component models
  
  LOGICAL, PROTECTED :: lrestart = .false.

  !> Flag: True, if model run should create restart at experiment end.
  !  This is independent from the settings of the restart interval.
  LOGICAL, PROTECTED ::  lrestart_write_last = .FALSE.

  CHARACTER(len=filename_max), PROTECTED :: model_base_dir = ''

  !> namelist parameters (as raw character strings):

  CHARACTER(len=MAX_CALENDAR_STR_LEN)  :: calendar                 = ''
  CHARACTER(len=MAX_DATETIME_STR_LEN)  :: experimentReferenceDate  = ''
  CHARACTER(len=MAX_DATETIME_STR_LEN)  :: experimentStartDate      = ''
  CHARACTER(len=MAX_DATETIME_STR_LEN)  :: experimentStopDate       = ''
  CHARACTER(len=MAX_TIMEDELTA_STR_LEN) :: checkpointTimeIntval     = ''
  CHARACTER(len=MAX_TIMEDELTA_STR_LEN) :: restartTimeIntval        = ''
  
CONTAINS

  SUBROUTINE setInstitution(update_institute)
    CHARACTER(len=*) :: update_institute
    institution = ''  ! needed?
    institution = update_institute
  END SUBROUTINE setInstitution
  
  SUBROUTINE setRestart(lr)
    LOGICAL, INTENT(in) :: lr
    lrestart = lr
  END SUBROUTINE setRestart

  SUBROUTINE setRestartWriteLast(lr)
    LOGICAL, INTENT(in) :: lr
    lrestart_write_last = lr
  END SUBROUTINE setRestartWriteLast

  LOGICAL FUNCTION isRestart()
    isRestart = lrestart
  END FUNCTION isRestart

  SUBROUTINE setModelBaseDir(mbd)
    CHARACTER(len=*), INTENT(in) :: mbd
    model_base_dir = mbd
  END SUBROUTINE setModelBaseDir
  
  CHARACTER(len=filename_max) FUNCTION getModelBaseDir()
    getModelBaseDir = model_base_dir 
  END FUNCTION getModelBaseDir

  SUBROUTINE addModel()
    no_of_models = no_of_models+1
  END SUBROUTINE addModel

  INTEGER FUNCTION noOfModels()
    noOfModels = no_of_models
  END FUNCTION noOfModels
  
  SUBROUTINE setExpRefdate(experimentReferenceDate)   
    CHARACTER(len=*), INTENT(in) :: experimentReferenceDate   
    tc_exp_refdate => newDatetime(experimentReferenceDate)
  END SUBROUTINE setExpRefdate

  SUBROUTINE setExpStartdate(experimentStartDate)   
    CHARACTER(len=*), INTENT(in) :: experimentStartDate   
    tc_exp_startdate => newDatetime(experimentStartDate)   
  END SUBROUTINE setExpStartdate

  SUBROUTINE setExpStopdate(experimentStopDate)
    CHARACTER(len=*), INTENT(in) :: experimentStopDate
    tc_exp_stopdate => newDatetime(experimentStopDate)
  END SUBROUTINE setExpStopdate

  SUBROUTINE setStartdate(startdate)
    CHARACTER(len=*), INTENT(in) :: startdate
    tc_startdate => newDatetime(startdate)
  END SUBROUTINE setStartdate

  SUBROUTINE setStopdate(stopdate)
    CHARACTER(len=*), INTENT(in) :: stopdate
    tc_stopdate => newDatetime(stopdate)
  END SUBROUTINE setStopdate

  SUBROUTINE setCheckpointTimeInterval(checkpointTimeIntval)
    CHARACTER(len=*), INTENT(in) :: checkpointTimeIntval
    tc_dt_checkpoint => newTimedelta(checkpointTimeIntval)
  END SUBROUTINE setCheckpointTimeInterval
  
  SUBROUTINE setRestartTimeInterval(restartTimeIntval)
    CHARACTER(len=*), INTENT(in) :: restartTimeIntval   
    tc_dt_restart => newTimedelta(restartTimeIntval)
  END SUBROUTINE setRestartTimeInterval

  SUBROUTINE setCurrentdate(current_date)
    CHARACTER(len=*), INTENT(in) :: current_date
    tc_current_date => newDatetime(current_date)
  END SUBROUTINE setCurrentdate


END MODULE mo_master_config
