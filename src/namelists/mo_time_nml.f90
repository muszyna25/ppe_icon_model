!>
!!        
!! @par Revision History
!!       Kristina Froehlich, MPI-M 2011-07-05
!! first implementation for all time control variables
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!!
MODULE mo_time_nml

  USE mo_kind,                  ONLY: wp
  USE mo_datetime,              ONLY: proleptic_gregorian,                            &
                                    & julian_gregorian, cly360
  USE mtime,                    ONLY: MAX_CALENDAR_STR_LEN
  USE mo_time_config,           ONLY: cfg_dt_restart          => dt_restart,          &
    &                                 cfg_calendar            => calendar,            &
    &                                 cfg_is_relative_time    => is_relative_time,    &
    &                                 cfg_ini_datetime_string => ini_datetime_string, &
    &                                 cfg_end_datetime_string => end_datetime_string, &
    &                                 restart_ini_datetime_string,                    &
    &                                 restart_end_datetime_string,                    &
    &                                 restart_calendar
  USE mo_io_units,              ONLY: nnml, nnml_output
  USE mo_master_config,         ONLY: isRestart
  USE mo_namelist,              ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_mpi,                   ONLY: my_process_is_stdio 

  USE mo_io_restart_namelist,   ONLY: open_and_restore_namelist, close_tmpfile, &
                                    & open_tmpfile, store_and_close_namelist
  USE mo_nml_annotate,          ONLY: temp_defaults, temp_settings

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: read_time_namelist

  !--------------------
  ! Namelist variables
  !--------------------

  CHARACTER(LEN = *), PARAMETER :: modname = "mo_time_nml"

  INTEGER            ::  calendar         ! calendar type
  REAL(wp)           ::  dt_restart       ! [s] length of a restart cycle 

  CHARACTER(len=32)  ::  ini_datetime_string, end_datetime_string

  !> LOGICAL is_relative_time: .TRUE., if time loop shall start with
  !> step 0 regardless whether we are in a standard run or in a
  !> restarted run (which means re-initialized run):
  LOGICAL            ::  is_relative_time

  NAMELIST /time_nml/ calendar,            &
    &                 ini_datetime_string, &
    &                 end_datetime_string, &
    &                 dt_restart,          &
    &                 is_relative_time

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


  !-------------------------------------------------------------------------
  !>
  !!  Initialization of variables that contain general information.
  !!
  !!  Initialization of variables that contain general information
  !!  about the model run. The configuration is read from
  !!  namelist 'time_nml'.
  !!
  !! @par Revision History
  !!  Reading of the namelist and checking of the validity were
  !!  in some other modules in the earlier version of the shallow water model.
  !!  Moved to this module by Hui Wan, MPI-M (2007-02-23)
  !!  The character parameter <i>routine</i> was introduced and used
  !!  for error information by Hui Wan, MPI-M (2007-02-23).
  !!  Modified by Almut Gassmann, MPI-M (2008-09-23)
  !!  - introduced  lidealized and lshallow_water
  !!  Modified by Marco Giorgetta, MPI-M (2009-02-23)
  !!  - lidealized replaced by ltestcase
  !!  Modification by Constantin Junk, MPI-M (2010-02-22)
  !!  - changes to consistency checks
  !!  Modification by Daniel Reinert, DWD (2011-11-10)
  !!  - in order for the restart to produce bitwise identical results, 
  !!    it is necessary to read caltime/calday directly from restart. So 
  !!    far, caltime/calday had been re-computed from the date. However, 
  !!    this led to significant roundoff errors.
  !!
  SUBROUTINE read_time_namelist( filename )

   CHARACTER(LEN=*), INTENT(IN)        :: filename
   ! local variables
   CHARACTER(len=*), PARAMETER         :: routine = modname//'::read_time_namelist'

   INTEGER                             :: istat, funit
   INTEGER                             :: iunit

   !------------------------------------------------------------------------
   ! Default values
   !------------------------------------------------------------------------

   ! Note: The default needs to be empty, since there exist
   ! concurrent namelist parameters to specify these values:
   calendar  = -1   ! unspecified

   ! Initialize start and end date by empty strings, which means
   ! "undefined". The same information may be set through other
   ! namelist parameters from "master_time_control_nml".
   ini_datetime_string = ""     ! initial date and time
   end_datetime_string = ""     ! end date and time
   dt_restart          = 0._wp  ! length of restart cycle (unspecified)

    is_relative_time = .FALSE.

    IF (isRestart()) THEN
 
      ! 2.1 Overwrite the defaults above by values in the restart file

      funit = open_and_restore_namelist('time_nml')
      READ(funit,NML=time_nml)
      CALL close_tmpfile(funit) 

      ! store the namelist settings originating from the restart file:
      restart_calendar            = calendar
      restart_ini_datetime_string = ini_datetime_string
      restart_end_datetime_string = end_datetime_string

    END IF

   !------------------------------------------------------------------------
   !  Read user's (new) specifications. (Done so far by all MPI processes)
   !------------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml('time_nml', STATUS=istat)
    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, time_nml)  ! write defaults to temporary text file
    END IF
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, time_nml)                                      ! overwrite default settings
      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, time_nml)  ! write settings to temporary text file
      END IF
    END SELECT
    CALL close_nml

    !----------------------------------------------------
    ! 4. Fill the configuration state
    !----------------------------------------------------

    ! Convert the calendar setting (which is an integer value for this
    ! namelist) into a string. The naming scheme is then compatible
    ! with concurrent namelist settings of the calendar (mtime):
    cfg_calendar = calendar_index2string(calendar)

    cfg_ini_datetime_string = ini_datetime_string
    cfg_end_datetime_string = end_datetime_string

    cfg_dt_restart          = dt_restart
    cfg_is_relative_time    = is_relative_time

    !-----------------------------------------------------
    ! Store the namelist for restart
    !-----------------------------------------------------
    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=time_nml)
      CALL store_and_close_namelist(funit, 'time_nml')
    ENDIF
    ! write the contents of the namelist to an ASCII file
    IF(my_process_is_stdio()) WRITE(nnml_output,nml=time_nml)

 END SUBROUTINE read_time_namelist

END MODULE mo_time_nml
