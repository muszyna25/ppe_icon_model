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

  USE mo_kind,                  ONLY: wp, i8
  USE mo_datetime,              ONLY: proleptic_gregorian, time_to_date, &
                                    & date_to_time, string_to_datetime
  USE mo_time_config,           ONLY: time_config
  USE mo_io_units,              ONLY: nnml, nnml_output
  USE mo_master_control,        ONLY: use_restart_namelists, isRestart
  USE mo_namelist,              ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_mpi,                   ONLY: my_process_is_stdio 

  USE mo_io_restart_attributes, ONLY: t_RestartAttributeList, getRestartAttributes
  USE mo_io_restart_namelist,   ONLY: open_and_restore_namelist, close_tmpfile, &
                                    & open_tmpfile, store_and_close_namelist
  USE mo_nml_annotate,          ONLY: temp_defaults, temp_settings

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: read_time_namelist

  !--------------------
  ! Namelist variables
  !--------------------

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

   CHARACTER(LEN=*), INTENT(IN) :: filename
   INTEGER  :: istat, funit,calendar_old
   CHARACTER(len=32) :: ini_datetime_string_old

   INTEGER(i8) :: restart_calday
   REAL(wp)    :: restart_caltime
   REAL(wp)    :: restart_daysec
   INTEGER     :: iunit
   TYPE(t_RestartAttributeList), POINTER :: restartAttributes


   !0!CHARACTER(len=*), PARAMETER ::  routine = 'mo_time_nml:read_time_namelist'

   !------------------------------------------------------------------------
   ! Default values
   !------------------------------------------------------------------------

   ! initial date and time
    calendar            = proleptic_gregorian
    ini_datetime_string = "2008-09-01T00:00:00Z"
   !
   ! end date and time
    end_datetime_string = "2008-09-01T01:40:00Z"
   !
   ! length of integration = (number of timesteps)*(length of timestep)
   ! - If nsteps is set to a non-zero positive value, then the end date is computed
   !   from the initial date and time, the time step dtime, and nsteps.
   ! - Else if run_day, run_hour, run_minute or run_second is set to a non-zero,
   !   positive value, then the initial date and time and the run_... variables are
   !   used to compute the end date and time and, using dtime, nsteps.
   !   Else nsteps is computed from the initial and end date and time and dtime.
   !
   ! length of restart cycle
    dt_restart     = 86400._wp*30._wp   ! = 30 days

    is_relative_time = .FALSE.


    restartAttributes => getRestartAttributes()
    IF (ASSOCIATED(restartAttributes)) THEN
 
      ! 2.1 Overwrite the defaults above by values in the restart file
      funit = open_and_restore_namelist('time_nml')
      READ(funit,NML=time_nml)
      CALL close_tmpfile(funit) 

      calendar_old            = calendar
      ini_datetime_string_old = ini_datetime_string

      ! 2.2 Inquire the date/time at which the previous run stopped

      restart_calday = restartAttributes%getInteger('current_calday')
      restart_caltime = restartAttributes%getReal('current_caltime')
      restart_daysec = restartAttributes%getReal('current_daysec')

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

    time_config%calendar = calendar

    CALL string_to_datetime( ini_datetime_string, time_config%ini_datetime )

    time_config%ini_datetime%calendar = time_config%calendar
    CALL date_to_time( time_config%ini_datetime )

    IF (isRestart()) THEN
      ! In a resumed integration, if the calendar or initial date/time 
      ! is different from those in the restart file,
      ! we regard this integration as a new one with its own calendar. 
      ! Model time at which the previous run stopped is thus not relevant. 
      ! Simulation will start from the user-specified initial date/time,
      ! which is also the current model date/time.

      IF (time_config%calendar /= calendar_old .OR.     &
          ini_datetime_string_old /= ini_datetime_string) THEN

        time_config%cur_datetime = time_config%ini_datetime

      ELSE
        time_config%cur_datetime%calendar = time_config%calendar
        time_config%cur_datetime%caltime  = restart_caltime
        time_config%cur_datetime%calday   = restart_calday
        time_config%cur_datetime%daysec   = restart_daysec

        CALL time_to_date(time_config%cur_datetime) ! fill date time structure

      END IF

    ELSE
      ! In an initial run, current date/time is, naturally, the initial date/time
      time_config%cur_datetime = time_config%ini_datetime

    END IF !isRestart()

    CALL string_to_datetime( end_datetime_string, time_config%end_datetime ) 

    time_config%end_datetime%calendar = time_config%calendar
    CALL date_to_time( time_config%end_datetime )

    time_config%dt_restart       = dt_restart
    time_config%is_relative_time = is_relative_time

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
