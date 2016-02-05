!>
!! <Short description of module for listings and indices>
!!
!! <Describe the concepts of the procedures and algorithms used in the module.>
!! <Details of procedures are documented below with their definitions.>
!! <Include any applicable external references inline as module::procedure,>
!! <external_procedure(), or by using @see.>
!! <Don't forget references to literature.>
!!
!! @author Hui Wan, MPI-M (2011-07-15)
!! @author Kristina Froehlich, MPI-M (2011-07-15)
!!
!!
!! @par Revision History
!! <Description of activity> by <name, affiliation> (<YYYY-MM-DD>)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_time_config

  USE mo_kind,                  ONLY: wp
  USE mtime,                    ONLY: datetime, timedelta, newDatetime, newTimedelta, &
    &                                 deallocateDatetime, MAX_CALENDAR_STR_LEN
  USE mo_impl_constants,        ONLY: proleptic_gregorian,                            &
                                    & julian_gregorian, cly360
  USE mo_util_string,           ONLY: tolower
 
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: calendar_index2string
  PUBLIC :: ini_datetime_string, end_datetime_string, icalendar
  PUBLIC :: restart_ini_datetime_string, restart_end_datetime_string, restart_calendar
  PUBLIC :: dt_restart, is_relative_time
  PUBLIC :: t_time_config, time_config
  PUBLIC :: setExpRefdate
  PUBLIC :: setExpStartdate, setExpStopdate
  PUBLIC :: setStartdate, setStopdate
  PUBLIC :: setCheckpointTimeInterval
  PUBLIC :: setRestartTimeInterval
  PUBLIC :: setCurrentdate
  PUBLIC :: setIsRelativeTime
  PUBLIC :: setTimeConfigCalendar
  PUBLIC :: setModelTimeStep

  !> namelist parameters (as raw character strings):
  !
  ! these are the namelist settings originating from the restart file:
  CHARACTER(len=32)                   :: restart_ini_datetime_string
  CHARACTER(len=32)                   :: restart_end_datetime_string
  INTEGER                             :: restart_calendar
  !
  !  these are namelist setting which may originate from the restart
  !  file, but with user modifications in the current run:
  INTEGER                             :: icalendar
  CHARACTER(len=32)                   :: ini_datetime_string
  CHARACTER(len=32)                   :: end_datetime_string
  REAL(wp)                            :: dt_restart          !< Length of restart cycle in seconds
  LOGICAL                             :: is_relative_time

  !>
  !! Derived type containing information for time control. 
  !!
  TYPE t_time_config

    ! from namelist 

    REAL(wp)         :: dt_restart         !< Length of restart cycle in seconds
    INTEGER          :: calendar           !< calendar type

    ! not directly from namelist  

    !> LOGICAL is_relative_time: .TRUE., if time loop shall start with
    !> step 0 regardless whether we are in a standard run or in a
    !> restarted run (which means re-initialized run):
    LOGICAL          ::  is_relative_time

    ! whole experiment and single run time information
    ! -----------------------------------------------------------------
    !
    ! experiment
    
    TYPE(datetime),  POINTER :: tc_exp_refdate   => NULL()

    TYPE(datetime),  POINTER :: tc_exp_startdate => NULL()
    TYPE(datetime),  POINTER :: tc_exp_stopdate  => NULL()

    ! single run 

    TYPE(datetime),  POINTER :: tc_startdate     => NULL()
    TYPE(datetime),  POINTER :: tc_stopdate      => NULL()

    ! current model date

    TYPE(datetime),  POINTER :: tc_current_date  => NULL()

    ! check point and restart time interval (needs to be equal for all
    ! component models)

    TYPE(timedelta), POINTER :: tc_dt_checkpoint => NULL()
    TYPE(timedelta), POINTER :: tc_dt_restart    => NULL()

    TYPE(timedelta), POINTER :: tc_dt_model      => NULL()
 
  END TYPE t_time_config
  !>
  !! 
  !! The actual variable
  !!
  TYPE(t_time_config), PROTECTED, SAVE :: time_config

  CHARACTER(LEN = *), PARAMETER :: modname = "mo_time_config"

CONTAINS

    !> Convert the calendar setting (which is an integer value for this
  !  namelist) into a string. The naming scheme is then compatible
  !  with concurrent namelist settings of the calendar (mtime).
  !
  FUNCTION calendar_index2string(icalendar) RESULT(ret)
    CHARACTER(LEN=MAX_CALENDAR_STR_LEN) :: ret
    INTEGER, INTENT(IN) :: icalendar
    ! local variables
    CHARACTER(len=*), PARAMETER ::  routine = modname//'::calendar_index2string'

    ret = ""
    SELECT CASE(icalendar)
    CASE(julian_gregorian)
      ret = 'julian gregorian'
    CASE(proleptic_gregorian)
      ret = 'proleptic gregorian'
    CASE(cly360)
      ret = '360 day year'
    END SELECT
  END FUNCTION calendar_index2string


  !> Convert the calendar setting (which is an integer value for this
  !  namelist) into a string. The naming scheme is then compatible
  !  with concurrent namelist settings of the calendar (mtime).
  !
  FUNCTION calendar_string2index(cal_str) RESULT(ret)
    INTEGER :: ret
    CHARACTER(LEN=*), INTENT(IN) :: cal_str
    ! local variables
    CHARACTER(len=*), PARAMETER ::  routine = modname//'::calendar_string2index'

    ret = -1
    IF (TRIM(tolower(cal_str)) == 'julian gregorian') THEN
      ret = julian_gregorian
    ELSE IF (TRIM(tolower(cal_str)) == 'proleptic gregorian') THEN
      ret = proleptic_gregorian
    ELSE IF (TRIM(tolower(cal_str)) == '360 day year') THEN
      ret = cly360
    END IF
  END FUNCTION calendar_string2index



  SUBROUTINE setExpRefdate(experimentReferenceDate)   
    CHARACTER(len=*), INTENT(in) :: experimentReferenceDate   
    time_config%tc_exp_refdate => newDatetime(experimentReferenceDate)
  END SUBROUTINE setExpRefdate

  SUBROUTINE setExpStartdate(experimentStartDate)   
    CHARACTER(len=*), INTENT(in) :: experimentStartDate   
    time_config%tc_exp_startdate => newDatetime(experimentStartDate)   
  END SUBROUTINE setExpStartdate

  SUBROUTINE setExpStopdate(experimentStopDate)
    CHARACTER(len=*), INTENT(in) :: experimentStopDate
    time_config%tc_exp_stopdate => newDatetime(experimentStopDate)
  END SUBROUTINE setExpStopdate

  SUBROUTINE setStartdate(startdate)
    CHARACTER(len=*), INTENT(in) :: startdate
    time_config%tc_startdate => newDatetime(startdate)
  END SUBROUTINE setStartdate

  SUBROUTINE setStopdate(stopdate)
    CHARACTER(len=*), INTENT(in) :: stopdate
    time_config%tc_stopdate => newDatetime(stopdate)
  END SUBROUTINE setStopdate

  SUBROUTINE setCheckpointTimeInterval(checkpointTimeIntval)
    CHARACTER(len=*), INTENT(in) :: checkpointTimeIntval
    time_config%tc_dt_checkpoint => newTimedelta(checkpointTimeIntval)
  END SUBROUTINE setCheckpointTimeInterval
  
  SUBROUTINE setRestartTimeInterval(restartTimeIntval)
    CHARACTER(len=*), INTENT(in) :: restartTimeIntval   
    time_config%tc_dt_restart => newTimedelta(restartTimeIntval)
  END SUBROUTINE setRestartTimeInterval

  SUBROUTINE setCurrentdate(current_date)
    CHARACTER(len=*), INTENT(in) :: current_date
    IF (ASSOCIATED(time_config%tc_current_date)) THEN
      CALL deallocateDatetime(time_config%tc_current_date)
    END IF
    time_config%tc_current_date => newDatetime(current_date)
  END SUBROUTINE setCurrentdate

  SUBROUTINE setTimeConfigCalendar(icalendar)
    INTEGER, INTENT(IN) :: icalendar
    time_config%calendar = icalendar
  END SUBROUTINE setTimeConfigCalendar

  SUBROUTINE setIsRelativeTime(lvalue)
    LOGICAL, INTENT(IN) :: lvalue
    time_config%is_relative_time = lvalue
  END SUBROUTINE setIsRelativeTime

  SUBROUTINE setModelTimeStep(modelTimeStep)
    CHARACTER(len=*), INTENT(in) :: modelTimeStep
    time_config%tc_dt_model => newTimedelta(modelTimeStep)
  END SUBROUTINE setModelTimeStep
 
END MODULE mo_time_config

