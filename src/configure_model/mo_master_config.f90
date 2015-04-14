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

  USE mtime,       ONLY: datetime, timedelta, newDatetime, newTimedelta
  USE mo_io_units, ONLY: filename_max
  
  IMPLICIT NONE
  
  PRIVATE

  PUBLIC :: addModel, noOfModels, maxNoOfModels
  PUBLIC :: setModelBaseDir, getModelBaseDir
  PUBLIC :: setRestart, isRestart
  PUBLIC :: setExpRefdate
  PUBLIC :: setExpStartdate, setExpStopdate
  PUBLIC :: setStartdate, setStopdate
  PUBLIC :: setCheckpointTimeInterval
  PUBLIC :: setRestartTimeInterval
  PUBLIC :: master_component_models
  
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
  
  ! checkpoint and restart time interval (needs to be equal for all component models) 
  
  TYPE(timedelta), POINTER, PROTECTED :: tc_dt_checkpoint => NULL()
  TYPE(timedelta), POINTER, PROTECTED :: tc_dt_restart => NULL()
  
  ! restart flag, required for consistent handling in all component models
  
  LOGICAL, PROTECTED :: lrestart = .false.

  CHARACTER(len=filename_max), PROTECTED :: model_base_dir = ''
  
CONTAINS

  SUBROUTINE setRestart(lr)
    LOGICAL, INTENT(in) :: lr
    lrestart = lr
  END SUBROUTINE setRestart

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

END MODULE mo_master_config
