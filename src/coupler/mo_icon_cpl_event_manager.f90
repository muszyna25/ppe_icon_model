!>
!! Event Manager for coupling events and other
!!
!! <Description>
!!
!! Routines are used to define events and check for events.
!!
!! call event_init () to initialise the whole event handling
!! call event_add  (event_id) to add a new event to the list
!! event_check ( event_id) TRUE/FALSE to check whether it is
!!  time for some action.
!! Note that the module itself does not have any programmed upper
!! limit wrt the number of events which can be handled.
!! 
!! @author Rene Redler, Max-Planck Institute for Meteorology, Germany
!!
!!
!! @par Revision History
!! first implementation by Rene Redler (2011-02-28)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_icon_cpl_event_manager

  USE mo_time_config, ONLY : time_config
  USE mo_io_config, ONLY   : dt_checkpoint
  USE mo_time_base, ONLY   : t_julian_date, add_time
  USE mo_icon_cpl,  ONLY   : initial_date
  USE mo_exception, ONLY   : finish

  IMPLICIT NONE

  PRIVATE

  TYPE t_event
     LOGICAL             :: l_is_initialised
     INTEGER             :: delta_time
     INTEGER             :: time_step
     INTEGER             :: elapsed_time
     INTEGER             :: event_time
     INTEGER             :: restart_time
     INTEGER             :: check_time
     INTEGER             :: lag
     TYPE(t_julian_date) :: current_date
  END TYPE t_event

  INTEGER, PARAMETER     :: event_inc = 8
  INTEGER                :: number_of_events
  TYPE(t_event), POINTER :: events(:) => NULL()

  INTEGER                :: ierr, id

  PUBLIC :: event_init, event_add, event_check, event_update, t_event, events

CONTAINS

  SUBROUTINE event_init

    number_of_events = event_inc

    IF ( .NOT. ASSOCIATED(events) ) &
    ALLOCATE ( events(number_of_events), STAT = ierr )
    IF ( ierr > 0 ) THEN
       CALL finish("event_init","Error allocating events")
    ENDIF

    DO id = 1, number_of_events
       events(id)%l_is_initialised = .FALSE.
       events(id)%delta_time       = 0
       events(id)%time_step        = 0
       events(id)%elapsed_time     = 0
       events(id)%event_time       = 0
       events(id)%restart_time     = 0
       events(id)%check_time       = 0
       events(id)%lag              = 0
       events(id)%current_date     = initial_date
    ENDDO

  END SUBROUTINE event_init

  ! ---------------------------------------------------------------------

  SUBROUTINE event_add ( event_id,   & !< out 
                         delta_time, & !< in
                         time_step,  & !< in
                         lag        )  !< in

    INTEGER, INTENT(in)          :: time_step  !< time step between two consecutive calls
    INTEGER, INTENT(in)          :: delta_time !< interval at which events shall ocure
    INTEGER, INTENT(in)          :: lag        !< lag times time_step
    INTEGER, INTENT(out)         :: event_id   !< returned hanle

    INTEGER :: seconds
    INTEGER :: days

    INTEGER                      :: new_dim

    TYPE(t_event), POINTER       :: new_events(:)

    ! check for free event handles is the existing list

    DO id = 1, number_of_events
       IF ( .NOT. events(id)%l_is_initialised ) EXIT
    ENDDO

    event_id = id

    IF ( id > number_of_events ) THEN

       ! ----------------------------------------------------------------
       ! if old list is full extend the event list
       ! ----------------------------------------------------------------
       
       new_dim = number_of_events + event_inc

       ALLOCATE ( new_events(new_dim), STAT = ierr )

       IF ( ierr > 0 ) THEN
         CALL finish("event_add","Error allocating events")
       ENDIF

       ! ----------------------------------------------------------------
       ! copy old list
       ! ----------------------------------------------------------------
       
       new_events(1:number_of_events) = events(1:number_of_events)

       ! ----------------------------------------------------------------
       ! initialise new part
       ! ----------------------------------------------------------------

       new_events(number_of_events+1:new_dim)%l_is_initialised = .FALSE.
       new_events(number_of_events+1:new_dim)%delta_time       = 0
       new_events(number_of_events+1:new_dim)%time_step        = 0
       new_events(number_of_events+1:new_dim)%elapsed_time     = 0
       new_events(number_of_events+1:new_dim)%event_time       = 0
       new_events(number_of_events+1:new_dim)%restart_time     = 0
       new_events(number_of_events+1:new_dim)%check_time       = 0
       new_events(number_of_events+1:new_dim)%lag              = 0
       new_events(number_of_events+1:new_dim)%current_date     = initial_date

       ! ----------------------------------------------------------------
       ! reset pointer
       ! ----------------------------------------------------------------

       DEALLOCATE ( events, STAT = ierr )

       IF ( ierr > 0 ) THEN
         CALL finish("event_add"," Error deallocating events")
       ENDIF

       events => new_events

       ! ----------------------------------------------------------------
       ! increase counter
       ! ----------------------------------------------------------------

       number_of_events = new_dim

    ENDIF

    ! -------------------------------------------------------------------
    ! set the event id
    ! -------------------------------------------------------------------

    events(id)%l_is_initialised = .TRUE.
    events(id)%delta_time       = delta_time
    events(id)%time_step        = time_step
    events(id)%restart_time     = NINT(time_config%dt_restart)
    events(id)%check_time       = NINT(dt_checkpoint)

    IF ( lag > 0 ) THEN

       ! fast-forward internal event by lag coupling time steps

       events(id)%event_time    = lag * time_step
       events(id)%elapsed_time  = lag * time_step
       events(id)%lag           = lag

       ! Update date and time for this event

       seconds = lag * events(id)%time_step
       days    = 0

       CALL add_time ( days, seconds, events(id)%current_date )

    ENDIF

  END SUBROUTINE event_add

  ! ---------------------------------------------------------------------

  LOGICAL FUNCTION event_check ( event_id ) RESULT(l_action)

    INTEGER, INTENT(in) :: event_id

    INTEGER :: seconds
    INTEGER :: days

    l_action = .FALSE.

    IF ( events(event_id)%event_time == events(event_id)%delta_time ) THEN
       events(event_id)%event_time = 0
       l_action = .TRUE.
    ENDIF

  END FUNCTION event_check

  ! ---------------------------------------------------------------------

  SUBROUTINE event_update ( event_id )

    INTEGER, INTENT(in) :: event_id

    INTEGER :: seconds
    INTEGER :: days

    events(event_id)%event_time = &
    events(event_id)%event_time + events(event_id)%time_step

    events(event_id)%elapsed_time = &
    events(event_id)%elapsed_time + events(event_id)%time_step

    seconds = events(event_id)%time_step
    days    = 0

    ! Update date and time for the next event

    CALL add_time ( days, seconds, events(event_id)%current_date )

  END SUBROUTINE event_update

END MODULE mo_icon_cpl_event_manager