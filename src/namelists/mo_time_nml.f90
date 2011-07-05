!>
!!        
!! @par Revision History
!!       Kristina Froehlich, MPI-M 2011-07-05
!! first implementation for all time control variables
!!
!! @par Copyright
!! 2002-2006 by DWD and MPI-M
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

  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: message, message_text, finish
  USE mo_impl_constants,     ONLY: max_char_length, max_dom, max_ntracer
  USE mo_physical_constants, ONLY: grav
  USE mo_datetime,           ONLY: t_datetime, proleptic_gregorian,          &
                                 & date_to_time, add_time, print_datetime_all
  USE mo_io_units,           ONLY: nnml, nnml_output
  USE mo_namelist,           ONLY: position_nml, positioned
  USE mo_mpi,                ONLY: p_pe, p_io
  USE mo_master_nml,         ONLY: lrestart
  USE mo_io_restart_attributes, ONLY: get_restart_attribute
  USE mo_io_restart_namelist,ONLY: open_tmpfile, store_and_close_namelist,   &
                                 & open_and_restore_namelist, close_tmpfile

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: version = '$Id$'

  PUBLIC
  PRIVATE :: time_ctl, calendar,                                              &
    &        ini_year, ini_month, ini_day, ini_hour, ini_minute, ini_second, &
    &        end_year, end_month, end_day, end_hour, end_minute, end_second, &
    &                             run_day, run_hour, run_minute, run_second


  ! ------------------------------------------------------------------------
  ! 1.0 Namelist variables and auxiliary parameters
  ! ------------------------------------------------------------------------

  ! initialization
  ! --------------
  INTEGER            :: iinit               ! model initialization:

  ! time information
  ! ----------------
  !
  ! calendar type
  INTEGER            :: calendar
  !
  ! initial date and time
  ! - namelist variables
  INTEGER            :: ini_year, ini_month, ini_day
  INTEGER            :: ini_hour, ini_minute
  REAL(wp)           :: ini_second
  ! - data and time structure
  TYPE(t_datetime)   :: ini_datetime
  !
  ! current model time, not a namelist variable
  TYPE(t_datetime)   :: current_datetime
  !
  ! end date and time
  ! - namelist variables
  INTEGER            :: end_year, end_month, end_day
  INTEGER            :: end_hour, end_minute
  REAL(wp)           :: end_second
  ! - data and time structure
  TYPE(t_datetime)   :: end_datetime
  !
  ! run length
  ! - in day,hr,min,sec
  INTEGER            :: run_day
  INTEGER            :: run_hour, run_minute
  REAL(wp)           :: run_second
  ! - in time steps
  INTEGER            :: nsteps              ! number of time steps
  REAL(wp)           :: dtime               ! [s] length of a time step

  ! restart interval
  ! ----------------
  REAL(wp)           :: dt_restart          ! [s] length of a restart cycle 

  ! timer
  ! -----
  LOGICAL  :: ltimer     ! if .TRUE.,  the timer is switched on
  INTEGER  :: timers_level = 1  ! what level of timers to run

  NAMELIST /time_ctl/calendar,        &
    &                ini_year, ini_month,  ini_day,        &
    &                ini_hour, ini_minute, ini_second,     &
    &                end_year, end_month,  end_day,        &
    &                end_hour, end_minute, end_second,     &
    &                run_day,                              &
    &                run_hour, run_minute, run_second,     &
    &                dt_restart, dtime


CONTAINS
  !-------------------------------------------------------------------------
  !>
  !!  Initialization of variables that contain general information.
  !!
  !!  Initialization of variables that contain general information
  !!  about the model run. The configuration is read from
  !!  namelist 'time_ctl'.
  !!
  !! @par Revision History
  !!  Reading of the namelist and checking of the validity were
  !!  in some other modules in the earlier version of the shallow water model.
  !!  Moved to this module by Hui Wan, MPI-M (2007-02-23)
  !!  The character parameter <i>routine</i> was introduced and used
  !!  for error information by Hui Wan, MPI-M (2007-02-23).
  !!  Modified by Almut Gassmann, MPI-M (2008-09-23)
  !!  - introduced i_cell_type, lidealized and lshallow_water
  !!  Modified by Marco Giorgetta, MPI-M (2009-02-23)
  !!  - lidealized replaced by ltestcase
  !!  Modification by Constantin Junk, MPI-M (2010-02-22)
  !!  - changes to consistency checks
  !!
  SUBROUTINE time_nml_setup
                                               
   INTEGER  :: istat, funit, calendar_old
   INTEGER  :: ini_year_old, ini_month_old, ini_day_old, ini_hour_old, ini_minute_old
   INTEGER  :: restart_year, restart_month, restart_day, restart_hour, restart_minute
   REAL(wp) :: ini_second_old
   REAL(wp) :: restart_second
   REAL(wp) :: cur_datetime_calsec, end_datetime_calsec, length_sec

   CHARACTER(len=max_char_length), PARAMETER ::   &
            &  routine = 'mo_time_nml/time_nml_setup'

   ! initial date and time
   calendar       = proleptic_gregorian
   ini_year       = 2008
   ini_month      = 9
   ini_day        = 1
   ini_hour       = 0
   ini_minute     = 0
   ini_second     = 0.0_wp
   !
   ! end date and time
   end_year       = 2008
   end_month      = 9
   end_day        = 1
   end_hour       = 1
   end_minute     = 40
   end_second     = 0.0_wp
   !
   ! length of integration = (number of timesteps)*(length of timestep)
   ! - If nsteps is set to a non-zero positive value, then the end date is computed
   !   from the initial date and time, the time step dtime, and nsteps.
   ! - Else if run_day, run_hour, run_minute or run_second is set to a non-zero,
   !   positive value, then the initial date and time and the run_... variables are
   !   used to compute the end date and time and, using dtime, nsteps.
   !   Else nsteps is computed from the initial and end date and time and dtime.
   !
   ! initialize run_... variables with zero
   run_day        = 0
   run_hour       = 0
   run_minute     = 0
   run_second     = 0.0_wp
   !
   ! initialize nsteps with zero
   nsteps         = 0
   !
   ! length of restart cycle
   dt_restart     = 86400._wp*30._wp   ! = 30 days
   !
   ! time step
   dtime          = 600._wp   ! [s] for R2B04 + semi-implicit time steppping

    !------------------------------------------------------------------------                  
    ! 2. If this is a resumed integration...
    !------------------------------------------------------------------------                  
    IF (lrestart) THEN                                                                 

      ! 2.1 Overwrite the defaults above by values in the restart file

      funit = open_and_restore_namelist('time_ctl')
      READ(funit,NML=time_ctl)
      CALL close_tmpfile(funit)

      ! 2.2 Save the calendar and initial date/time of the old run

      calendar_old    = calendar
      ini_year_old    = ini_year
      ini_month_old   = ini_month
      ini_day_old     = ini_day
      ini_hour_old    = ini_hour
      ini_minute_old  = ini_minute
      ini_second_old  = ini_second

      ! 2.2 Inquire the date/time at which the previous run stopped

      CALL get_restart_attribute( 'current_year'  , restart_year   )
      CALL get_restart_attribute( 'current_month' , restart_month  )
      CALL get_restart_attribute( 'current_day'   , restart_day    )
      CALL get_restart_attribute( 'current_hour'  , restart_hour   )
      CALL get_restart_attribute( 'current_minute', restart_minute )
      CALL get_restart_attribute( 'current_second', restart_second )

    END IF

    !------------------------------------------------------------------------
    ! 3. Read user's (new) specifications. (Done so far by all MPI processes)
    !------------------------------------------------------------------------
    CALL position_nml('time_ctl', STATUS=istat)
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, time_ctl)
    END SELECT

    ! time step
    IF (dtime  <= 0._wp) CALL finish(routine,'"dtime" must be positive')

    !---------------------------------------------------------------
    ! 5. Set up model time
    !---------------------------------------------------------------
    ! 5.1 Initial date and time

    ini_datetime%calendar = calendar
    ini_datetime%year     = ini_year
    ini_datetime%month    = ini_month
    ini_datetime%day      = ini_day
    ini_datetime%hour     = ini_hour
    ini_datetime%minute   = ini_minute
    ini_datetime%second   = ini_second

    CALL date_to_time(ini_datetime) ! fill date time structure
    CALL message(' ',' ')
    CALL message(routine,'Initial date and time')
    CALL message(routine,'---------------------')
    CALL print_datetime_all(ini_datetime)  ! print all date and time components

    ! 5.2 Current date and time:

    IF (lrestart) THEN
      ! In a resumed integration, if the calendar or initial date/time 
      ! is different from those in the restart file,
      ! we regard this integration as a new one with its own calendar. 
      ! Model time at which the previous run stopped is thus not relevant. 
      ! Simulation will start from the user-specified initial date/time,
      ! which is also the current model date/time.

      IF (calendar  /=calendar_old   .OR.                                 &
          ini_year  /=ini_year_old   .OR. ini_month  /=ini_month_old .OR. &
          ini_day   /=ini_day_old    .OR. ini_hour   /=ini_hour_old  .OR. &
          ini_minute/=ini_minute_old .oR. ini_second /=ini_second_old     ) THEN

        current_datetime = ini_datetime

      ELSE
      ! Otherwise we start from the point when the previous integration stopped.

        current_datetime%calendar = calendar
        current_datetime%year     = restart_year
        current_datetime%month    = restart_month
        current_datetime%day      = restart_day
        current_datetime%hour     = restart_hour
        current_datetime%minute   = restart_minute
        current_datetime%second   = restart_second

        CALL date_to_time(current_datetime) ! fill date time structure
      END IF

    ELSE
      ! In an initial run, current date/time is, naturally, the initial date/time
      current_datetime = ini_datetime

    END IF !lrestart

    CALL message(' ',' ')
    CALL message(' ',' ')
    CALL message(routine,'Current date and time')
    CALL message(routine,'---------------------')
    CALL print_datetime_all(current_datetime)  ! print all date and time components

    ! 5.3 End date and time, and length of integration
    !     Here we define "nsteps" as the number of time steps THIS integration
    !     will last, regardless of the restart status.

    !
    IF (run_day/=0 .OR. run_hour/=0 .OR. run_minute/=0 .OR. run_second/=0.0_wp) THEN
      IF (run_day    < 0    ) CALL finish(routine,'"run_day" must not be negative')
      IF (run_hour   < 0    ) CALL finish(routine,'"run_hour" must not be negative')
      IF (run_minute < 0    ) CALL finish(routine,'"run_minute" must not be negative')
      IF (run_second < 0._wp) CALL finish(routine,'"run_second" must not be negative')
      !
      end_datetime = current_datetime
      CALL add_time(run_second,run_minute,run_hour,run_day,end_datetime)
      !
      cur_datetime_calsec = (REAL(current_datetime%calday,wp)+current_datetime%caltime) &
        &                   *REAL(current_datetime%daylen,wp)
      end_datetime_calsec = (REAL(end_datetime%calday,wp)+end_datetime%caltime) &
        &                   *REAL(end_datetime%daylen,wp)
      nsteps=INT((end_datetime_calsec-cur_datetime_calsec)/dtime)
      !
    ELSE
      ! compute nsteps from current_datetime, end_datetime and dtime
      end_datetime%calendar = calendar
      end_datetime%year     = end_year
      end_datetime%month    = end_month
      end_datetime%day      = end_day
      end_datetime%hour     = end_hour
      end_datetime%minute   = end_minute
      end_datetime%second   = end_second
      CALL date_to_time      (end_datetime) ! fill date time structure
      !
      cur_datetime_calsec = (REAL(current_datetime%calday,wp)+current_datetime%caltime) &
        &                   *REAL(current_datetime%daylen,wp)
      end_datetime_calsec = (REAL(end_datetime%calday,wp)+end_datetime%caltime) &
        &                   *REAL(end_datetime%daylen,wp)
      IF (end_datetime_calsec < cur_datetime_calsec) &
        & CALL finish(routine,'The end date and time must not be '// &
        &            'before the current date and time')
      !
    END IF
    !..................................................................
    !
    CALL message(' ',' ')
    CALL message(routine,'End date and time')
    CALL message(routine,'-----------------')
    CALL print_datetime_all(end_datetime)  ! print all date and time components

    CALL message(' ',' ')
    CALL message(routine,'Length of restart cycle')
    CALL message(routine,'-----------------------')
    WRITE(message_text,'(a,f10.2,a,f16.10,a)') &
         &'dt_restart :',dt_restart,' seconds =', dt_restart/86400._wp, ' days'
    CALL message(routine,message_text)

    CALL message(' ',' ')
    CALL message(routine,'Length of this run')
    CALL message(routine,'------------------')
    WRITE (message_text,'(a,f7.2)') 'dtime [s] :',dtime
    CALL message(routine,message_text)
    WRITE (message_text,'(a,i7)')   'nsteps    :',nsteps
    CALL message(routine,message_text)
    CALL message(' ',' ')


    !-----------------------------------------------------
    ! Store the namelist for restart
    !-----------------------------------------------------
    funit = open_tmpfile()
    WRITE(funit,NML=time_ctl)
    CALL store_and_close_namelist(funit, 'time_ctl')

    ! write the contents of the namelist to an ASCII file
    IF(p_pe == p_io) WRITE(nnml_output,nml=time_ctl)

 END SUBROUTINE time_nml_setup

END MODULE mo_time_nml
