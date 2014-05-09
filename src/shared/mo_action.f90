!>
!! Routines for defining, initializing and performing action events
!!
!! Routines for defining, initializing and performing action events
!!
!! *****************************************************************
!!            RECIPE FOR CREATING A NEW ACTION EVENT
!! *****************************************************************
!! 1. Define new actionID and action-object (see 'List of available actions' 
!!    in 'mo_action')
!! 2. Assign new action to variables of your choice.
!!    E.g. add the following code snippet to an add_var/add_ref of your choice:
!!    action_list=actions(new_action(ACTION_XXX,'PTXXH'), new_action(...), ...)
!!    ACTION_XXX is the actionID defined in step 1, and PTXXH is the 
!!    interval at which the action should be triggered  
!! 3. Add another 'action_collect_vars' CALL to 'action_init', in order to assign
!!    fields to your particular action. I.e. this is the reverse operation 
!!    of assigning actions to fields as done in step 2.
!! 4. Write your own action-Routine 
!!    (see e.g. routine reset_action for actionID=ACTION_RESET)
!! 5. Add CALL of your newly defined action-routine (see step 4) to the code 
!!   (location depends on your needs)   
!!
!! @author Daniel Reinert, DWD
!!
!!
!! @par Revision History
!! Initial revision by Daniel Reinert, DWD (2014-01-09)
!!
!! @par Copyright
!! 2002-2010 by DWD and MPI-M
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
MODULE mo_action

  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: message, message_text
  USE mtime,                 ONLY: event, newEvent, datetime, newDatetime,    &
    &                              isCurrentEventActive, deallocateDatetime,  &
    &                              MAX_DATETIME_STR_LEN, PROLEPTIC_GREGORIAN, &
    &                              MAX_EVENTNAME_STR_LEN, newTimedelta,       &
    &                              timedelta, newTimedelta, deallocateTimedelta
  USE mo_run_config,         ONLY: msg_level
  USE mo_mtime_extensions,   ONLY: get_datetime_string, getTimeDeltaFromDateTime, &
    &                              getPTStringFromMS
  USE mo_action_types,       ONLY: t_var_action
  USE mo_time_config,        ONLY: time_config
  USE mo_var_list,           ONLY: nvar_lists, var_lists
  USE mo_linked_list,        ONLY: t_list_element
  USE mo_var_list_element,   ONLY: t_var_list_element


  IMPLICIT NONE

  PRIVATE

  ! TYPES
  PUBLIC  :: t_var_element_ptr
  PUBLIC  :: t_action_obj

  ! PROCEDURES/FUNCTIONS
  PUBLIC  :: action_init
  PUBLIC  :: reset_action

  ! VARIABLES
  PUBLIC :: reset_act


  ! List of available actions
  !
  INTEGER, PARAMETER, PUBLIC :: ACTION_RESET = 1   ! re-set field to 0


  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  INTEGER, PARAMETER :: NMAX_VARS = 50  ! maximum number of fields that can be 
                                        ! assigned to a single action

  ! base type for generating an array of pointers of type t_var_list_element
  !
  TYPE t_var_element_ptr
    TYPE(t_var_list_element), POINTER :: p
    TYPE(event)             , POINTER :: event     ! event from mtime library
  END TYPE t_var_element_ptr


  ! general definition of an action
  !
  TYPE t_action_obj
    INTEGER                    :: actionID                    ! action ID
    TYPE(t_var_element_ptr)    :: var_element_ptr(NMAX_VARS)  ! assigned variables
    INTEGER                    :: var_action_index(NMAX_VARS) ! index in var_element_ptr(10)%action

    INTEGER                    :: nvars                 ! number of variables for which 
                                                        ! this action is to be performed
  END TYPE t_action_obj



  ! actions
  !
  TYPE(t_action_obj) :: reset_act  ! action which resets field resetval%rval

CONTAINS

  !>
  !! Initialize actions by assigning variables to existing actions
  !!
  !! When generating a new field via add_var, it is possible to assign various 
  !! actions to this field. Here, we go the other way around. For each 
  !! available action, we loop over all fields and check, whether this 
  !! action must be performed for the particular field. If this is the case, 
  !! this field is assigned to the action.
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2014-01-13)
  !!
  SUBROUTINE action_init()

    ! 1. ACTION_RESET
    !
    CALL action_collect_vars(actionID=ACTION_RESET, act_obj=reset_act)

    ! 2. Add here the next collect_action CALL for ACTION_XXX
    !

  END SUBROUTINE action_init


  !> 
  !! Assign variables to specific actions/initialize action-object.
  !! 
  !! Loop over all variables and collect the variables names
  !! corresponding to the action @p action%actionID
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2014-01-13)
  !!
  SUBROUTINE action_collect_vars(actionID, act_obj)
    INTEGER            , INTENT(IN)   :: actionID
    TYPE(t_action_obj) , INTENT(OUT)  :: act_obj

    ! local variables
    CHARACTER(*), PARAMETER :: routine = TRIM("mo_var_list:action_collect_vars")
    INTEGER :: i, iact
    INTEGER :: nvars
    TYPE(t_list_element), POINTER :: element
    TYPE(t_var_action)  , POINTER :: action_list
    CHARACTER(LEN=2)                    :: str_actionID
    CHARACTER(LEN=128)                  :: intvl             ! action interval [PTnH]
    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: sim_start, sim_end
    CHARACTER(LEN=MAX_EVENTNAME_STR_LEN):: event_name
  !-------------------------------------------------------------------------

    ! init nvars
    nvars = 0

    ! compute sim_start, sim_end in a formate appropriate for mtime
    CALL get_datetime_string(sim_start, time_config%ini_datetime)
    CALL get_datetime_string(sim_end,   time_config%end_datetime)

    ! store actionID
    act_obj%actionID = actionID

    ! loop over all variable lists and variables
    !
    DO i = 1,nvar_lists
      element => NULL()

      LOOPVAR : DO
        IF(.NOT.ASSOCIATED(element)) THEN
          element => var_lists(i)%p%first_list_element
        ELSE
          element => element%next_list_element
        ENDIF
        IF(.NOT.ASSOCIATED(element)) EXIT LOOPVAR

        ! point to variable specific action list
        action_list => element%field%info%action_list

        ! Loop over all variable-specific actions
        !
        LOOPACTION: DO iact = 1,action_list%n_actions
 
          ! If the variable specific action fits, assign variable 
          ! to the corresponding action.
          IF (action_list%action(iact)%actionID == actionID) THEN

            ! Add field to action object
            nvars = nvars + 1
            act_obj%var_element_ptr(nvars)%p => element%field
            act_obj%var_action_index(nvars) = iact


            ! Create event for this specific field
            intvl = action_list%action(iact)%intvl
            write(str_actionID,'(i2)') actionID 
            event_name = 'act_ID'//TRIM(str_actionID)//'_'//TRIM(intvl)
            act_obj%var_element_ptr(nvars)%event =>newEvent(TRIM(event_name),TRIM(sim_start), &
              &                                       TRIM(sim_start), TRIM(sim_end), &
              &                                       TRIM(intvl))

!DR            action_list%action(iact)%event =>newEvent(TRIM(event_name),TRIM(sim_start), &
!DR              &                                       TRIM(sim_start), TRIM(sim_end), &
!DR              &                                       TRIM(intvl))

          END IF
        ENDDO  LOOPACTION ! loop over variable-specific actions

        IF(ASSOCIATED(action_list)) action_list => NULL()

      ENDDO LOOPVAR ! loop over vlist "i"
    ENDDO ! i = 1,nvar_lists

    ! set nvars
    act_obj%nvars = nvars

    IF (msg_level >= 11) THEN
      WRITE(message_text,'(a,i2,a)') 'Variables assigned to action ',act_obj%actionID,' :'
      CALL message(TRIM(routine),message_text)
      DO i=1,act_obj%nvars
      WRITE(message_text,'(a)') TRIM(act_obj%var_element_ptr(i)%p%info%name)
      CALL message('',message_text)
      ENDDO
    ENDIF

  END SUBROUTINE action_collect_vars


  !>
  !! ACTION_RESET: Reset fields to field%info%resetval%rval
  !!
  !! If time has come, event is triggered. This particular event/action 
  !! is only about resetting fields to resetval%rval.
  !! Works on reset_act object.
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2014-01-10)
  !!
  SUBROUTINE reset_action(slack)

    REAL(wp), INTENT(IN)      :: slack     !< allowed slack for event triggering  [s]
    ! local variables
    INTEGER :: ivar                        !< loop index for fields
    INTEGER :: var_action_idx              !< Index of this particular action in 
                                           !< field-specific action list

    TYPE(t_var_list_element), POINTER :: & !< Pointer to particular field
      &  field
    TYPE(event)             , POINTER :: & !< Pointer to variable specific event info
      &  this_event 
    TYPE(datetime),           POINTER :: & !< Current date in mtime format
      &  mtime_date 
    TYPE(datetime),           POINTER :: & !< Model start date in mtime format
      &  mtime_inidate 
    TYPE(timedelta),          POINTER :: & !< for forecast time computation
      &  time_range

    LOGICAL :: isactive

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: mtime_cur_datetime
    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: mtime_ini_datetime

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: str_slack       ! slack as string
    TYPE(timedelta), POINTER :: p_slack                    ! slack in 'timedelta'-Format

    CHARACTER(*), PARAMETER :: routine = TRIM("mo_action:reset_action")

  !-------------------------------------------------------------------------

    ! compute current datetime in a format appropriate for mtime
    CALL get_datetime_string(mtime_cur_datetime, time_config%cur_datetime)
    mtime_date  => newDatetime(TRIM(mtime_cur_datetime)) 

    ! compute allowed slack in PT-Format
    ! Use factor 999 instead of 1000, since no open intervall is available
    ! needed [trigger_date, trigger_date + slack[
    ! used   [trigger_date, trigger_date + slack]
    CALL getPTStringFromMS(INT(999._wp*slack),str_slack)
    ! get slack in 'timedelta'-format appropriate for isCurrentEventActive
    p_slack => newTimedelta(str_slack)


    ! Loop over all fields attached to this action
    DO ivar = 1, reset_act%nvars
      var_action_idx = reset_act%var_action_index(ivar)

      field      => reset_act%var_element_ptr(ivar)%p
      this_event => reset_act%var_element_ptr(ivar)%event
!DR      this_event => field%info%action_list%action(var_action_idx)%event


      ! Check whether event-pointer is associated.
      IF (.NOT. ASSOCIATED(this_event)) THEN
        WRITE (message_text,'(a,i2,a,a,a)')                           &
             'WARNING: action event ', var_action_idx, ' of field ',  &
              TRIM(field%info%name),': Event-Ptr is disassociated!'
        CALL message(routine,message_text)
      ENDIF


      ! Note that a second call to isCurrentEventActive will lead to 
      ! a different result! Is this a bug or a feature?
      ! triggers in interval [trigger_date + slack]
      isactive = isCurrentEventActive(this_event,mtime_date, plus_slack=p_slack)


      ! Check wheter the action 'reset_act' should be triggered for variable
      ! under consideration
      IF (isactive) THEN

        ! store latest triggering date
        field%info%action_list%action(var_action_idx)%lastActive = TRIM(mtime_cur_datetime)


        ! compute ini datetime in a format appropriate for mtime
        CALL get_datetime_string(mtime_ini_datetime, time_config%ini_datetime)
        mtime_inidate  => newDatetime(TRIM(mtime_ini_datetime))
        time_range     => newTimedelta("P01D")
        CALL getTimeDeltaFromDateTime(mtime_date, mtime_inidate, time_range)

! preparations
!!$        ! store forecast time (necessary for proper GRIB2 encoding)
!!$        field%info%volatile%fcast_time = 86400*time_range%day + 3600*time_range%hour &
!!$          &                            + 60*time_range%minute + time_range%second


        IF (msg_level >= 12) THEN
          WRITE(message_text,'(a,i2,a,a,a,a)') 'action ',reset_act%actionID, &
            &  ' triggered for ', TRIM(field%info%name), ' at ', TRIM(mtime_cur_datetime)
          CALL message(TRIM(routine),message_text)
        ENDIF

        ! re-set field to 'resetval'
        ! so far, works only for fields of type REAL
        field%r_ptr = field%info%resetval%rval

        ! cleanup
        !
        CALL deallocateTimedelta(time_range)
        CALL deallocateDatetime(mtime_inidate)
      ENDIF

    ENDDO


    ! cleanup
    !
    CALL deallocateDatetime(mtime_date)

  END SUBROUTINE reset_action



END MODULE mo_action

