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

  USE mtime,       ONLY: MAX_CALENDAR_STR_LEN, MAX_DATETIME_STR_LEN, &
    &                    MAX_TIMEDELTA_STR_LEN
  USE mo_io_units, ONLY: filename_max
  
  IMPLICIT NONE
  
  PRIVATE

  PUBLIC :: calendar_str
  PUBLIC :: addModel, noOfModels, maxNoOfModels
  PUBLIC :: setInstitution
  PUBLIC :: setModelBaseDir, getModelBaseDir
  PUBLIC :: master_component_models
  PUBLIC :: lrestart_write_last
  PUBLIC :: experimentReferenceDate, experimentStartDate, experimentStopDate
  PUBLIC :: checkpointTimeIntval, restartTimeIntval
  PUBLIC :: setRestart, setRestartWriteLast
  PUBLIC :: isRestart
  
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
  
 
  ! restart flag, required for consistent handling in all component models
  
  LOGICAL, PROTECTED :: lrestart = .false.

  !> Flag: True, if model run should create restart at experiment end.
  !  This is independent from the settings of the restart interval.
  LOGICAL, PROTECTED ::  lrestart_write_last = .TRUE.

  CHARACTER(len=filename_max), PROTECTED :: model_base_dir = ''

  !> namelist parameters (as raw character strings):

  CHARACTER(len=MAX_CALENDAR_STR_LEN)  :: calendar_str             = ''
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
  
END MODULE mo_master_config
