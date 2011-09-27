!>
!!        
!! @par Revision History
!!       Kristina Froehlich, MPI-M 2011-07-05
!! first implementation for all time control variables
!!
!! @par Copyright
!! 2002-2011 by DWD and MPI-M
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
!!
MODULE mo_time_nml

  USE mo_kind,                  ONLY: wp
  USE mo_datetime,              ONLY: proleptic_gregorian, &
                                    & date_to_time, string_to_datetime
  USE mo_time_config,           ONLY: time_config
  USE mo_io_units,              ONLY: nnml, nnml_output
  USE mo_master_control,     ONLY: is_restart_run
  USE mo_namelist,              ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_mpi,                   ONLY: my_process_is_stdio 

  USE mo_io_restart_attributes, ONLY: get_restart_attribute
  USE mo_io_restart_namelist,   ONLY: open_and_restore_namelist, close_tmpfile, &
                                    & open_tmpfile, store_and_close_namelist

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: read_time_namelist

  CHARACTER(len=*), PARAMETER, PRIVATE :: version = '$Id$'

  !--------------------
  ! Namelist variables
  !--------------------

  INTEGER            ::  calendar         ! calendar type
  REAL(wp)           ::  dt_restart       ! [s] length of a restart cycle 

  CHARACTER(len=32)  ::  ini_datetime_string, end_datetime_string

  NAMELIST /time_nml/ calendar,            &
    &                 ini_datetime_string, &
    &                 end_datetime_string, &
    &                 dt_restart

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
  !!
  SUBROUTINE read_time_namelist( filename )

   CHARACTER(LEN=*), INTENT(IN) :: filename
   INTEGER  :: istat, funit,calendar_old
   CHARACTER(len=32) :: ini_datetime_string_old
   INTEGER  :: restart_year, restart_month,  restart_day
   INTEGER  :: restart_hour, restart_minute
   REAL(wp) :: restart_second

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

  IF (is_restart_run()) THEN
 
   ! 2.1 Overwrite the defaults above by values in the restart file

      funit = open_and_restore_namelist('time_nml')
      READ(funit,NML=time_nml)
      CALL close_tmpfile(funit) 

      calendar_old            = calendar
      ini_datetime_string_old = ini_datetime_string

      ! 2.2 Inquire the date/time at which the previous run stopped

      CALL get_restart_attribute( 'current_year'  , restart_year   )
      CALL get_restart_attribute( 'current_month' , restart_month  )
      CALL get_restart_attribute( 'current_day'   , restart_day    )
      CALL get_restart_attribute( 'current_hour'  , restart_hour   )
      CALL get_restart_attribute( 'current_minute', restart_minute )
      CALL get_restart_attribute( 'current_second', restart_second )

  END IF

   !------------------------------------------------------------------------
   !  Read user's (new) specifications. (Done so far by all MPI processes)
   !------------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml('time_nml', STATUS=istat)
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, time_nml)
    END SELECT
    CALL close_nml

    !----------------------------------------------------
    ! 4. Fill the configuration state
    !----------------------------------------------------

    time_config%calendar = calendar

    CALL string_to_datetime( ini_datetime_string, time_config%ini_datetime )

    time_config%ini_datetime%calendar = time_config%calendar
    CALL date_to_time( time_config%ini_datetime )

    IF (is_restart_run()) THEN
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
         time_config%cur_datetime%year     = restart_year
         time_config%cur_datetime%month    = restart_month
         time_config%cur_datetime%day      = restart_day
         time_config%cur_datetime%hour     = restart_hour
         time_config%cur_datetime%minute   = restart_minute
         time_config%cur_datetime%second   = restart_second

        CALL date_to_time(time_config%cur_datetime) ! fill date time structure
      END IF

    ELSE
      ! In an initial run, current date/time is, naturally, the initial date/time
      time_config%cur_datetime = time_config%ini_datetime

    END IF !is_restart_run()

    CALL string_to_datetime( end_datetime_string, time_config%end_datetime ) 

    time_config%end_datetime%calendar = time_config%calendar
    CALL date_to_time( time_config%end_datetime )

    time_config%dt_restart = dt_restart

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
