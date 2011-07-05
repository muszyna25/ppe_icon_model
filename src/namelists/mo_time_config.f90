!>
!! <Short description of module for listings and indices>
!!
!! <Describe the concepts of the procedures and algorithms used in the module.>
!! <Details of procedures are documented below with their definitions.>
!! <Include any applicable external references inline as module::procedure,>
!! <external_procedure(), or by using @see.>
!! <Don't forget references to literature.>
!!
!! @author <name, affiliation>
!! @author <name, affiliation>
!!
!!
!! @par Revision History
!! <Description of activity> by <name, affiliation> (<YYYY-MM-DD>)
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
MODULE mo_time_config

  USE mo_kind,                  ONLY: wp
  USE mo_time_nml,              ONLY: nml_calendar, nml_ini_datetime,&
                                   &  nml_end_datetime, nml_dt_restart
  USE mo_exception,             ONLY: message, message_text, finish
  USE mo_datetime,              ONLY: t_datetime, &
                                   & date_to_time, add_time, print_datetime_all
  USE mo_master_nml,            ONLY: lrestart
  USE mo_io_restart_attributes, ONLY: get_restart_attribute
  USE mo_io_restart_namelist,ONLY: open_and_restore_namelist, close_tmpfile,&
                                  & open_tmpfile, store_and_close_namelist

  USE mo_mpi,                    ONLY: p_pe, p_io

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  !>
  !! Derived type containing variables for time control. 
  !!
  TYPE :: t_time_config

    INTEGER  :: calendar      !< calendar type 

    TYPE(t_datetime) :: ini_datetime  !< Starting time of model integration
    TYPE(t_datetime) :: end_datetime  !< Ending   time of model integration
    TYPE(t_datetime) :: current_datetime  !< Current  time model time 

    REAL(wp) :: dt_restart    !< Length of restart cycle in seconds

  END TYPE t_time_config 

  !! HW Comment: the character-type variables containing ini_ and end_time
  !! in the format "YYYYMMDDTHHMMSSZ" should be namelist variables
  !! and used for computing ini/end_datetime in this type.

  !>
  !! The state variable that contains the information
  !!
  TYPE(t_time_config) :: time_config


CONTAINS

SUBROUTINE time_setup


  ! time information
  ! ----------------
  !
  ! calendar type
  INTEGER          :: calendar
  !
  ! initial date and time
  ! - namelist variables
  INTEGER          :: ini_year, ini_month, ini_day
  INTEGER          :: ini_hour, ini_minute
  REAL(wp)         :: ini_second
  ! - data and time structure
  TYPE(t_datetime) :: ini_datetime
  !
  ! current model time, not a namelist variable
  TYPE(t_datetime) :: current_datetime
  !
  ! end date and time
  ! - namelist variables
  INTEGER          :: end_year, end_month, end_day
  INTEGER          :: end_hour, end_minute
  REAL(wp)         :: end_second
  ! - data and time structure
  TYPE(t_datetime) :: end_datetime
  !

  INTEGER:: calendar_old
   INTEGER  :: ini_year_old, ini_month_old, ini_day_old, ini_hour_old, ini_minute_old
   INTEGER  :: restart_year, restart_month, restart_day, restart_hour, restart_minute
   REAL(wp) :: ini_second_old
   REAL(wp) :: restart_second
   REAL(wp) :: cur_datetime_calsec, end_datetime_calsec, length_sec
  ! restart interval
  ! ----------------
  REAL(wp)         :: dt_restart          ! [s] length of a restart cycle 


   ! initial date and time

    !------------------------------------------------------------------------                  
    ! 2. If this is a resumed integration...
    !------------------------------------------------------------------------                  
    IF (lrestart) THEN                                                                 

   ! 2.1 Overwrite the defaults above by values in the restart file

      funit = open_and_restore_namelist('time_ctl')
      READ(funit,NML=run_ctl)
      CALL close_tmpfile(funit)


      ! 2.2 Save the calendar and initial date/time of the old run

      calendar_old    = nml_calendar
      ini_datetime    = nml_ini_datetime
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



    CALL date_to_time(end_datetime) ! fill date time structure
    CALL message(' ',' ')
    CALL message(routine,'End date and time')
    CALL message(routine,'-----------------')
    CALL print_datetime_all(end_datetime)  ! print all date and time components



    t_time_config%dt_restart =  nml_dt_restart
    CALL message(' ',' ')
    CALL message(routine,'Length of restart cycle')
    CALL message(routine,'-----------------------')
    WRITE(message_text,'(a,f10.2,a,f16.10,a)') &
         &'dt_restart :',dt_restart,' seconds =', dt_restart/86400._wp, ' days'
    CALL message(routine,message_text)


    !-----------------------------------------------------
    ! Store the namelist for restart
    !-----------------------------------------------------
    funit = open_tmpfile()
    WRITE(funit,NML=time_ctl)
    CALL store_and_close_namelist(funit, 'time_ctl')


END SUBROUTINE time_setup

END MODULE mo_time_config

