!>
!! Routines for defining, initializing and executing action events
!!
!! Routines for defining, initializing and executing action events
!!
!! *****************************************************************
!!            RECIPE FOR CREATING A NEW ACTION EVENT
!! *****************************************************************
!! 1. Define a new actionID in mo_action
!! 2. Assign new actionID to variables of your choice.
!!    E.g. add the following code snippet to an add_var/add_ref of your choice:
!!    action_list=actions(new_action(ACTION_XXX,'PTXXH'), new_action(...), ...)
!!    ACTION_XXX is the actionID defined in step 1, and PTXXH is the 
!!    interval at which the action should be triggered.  
!! 3. Create an extension of the abstract type t_action_obj and overwrite 
!!    the deferred procedure 'kernel' with your action-specific kernel-routine 
!!    (to be defined in step 5).
!! 4. Create a variable (object) of the type defined in step 3.
!! 5. Write your own action-Routine (action-kernel). This is the routine which actually does 
!!    the work. (see e.g. routine 'reset_kernel' for actionID=ACTION_RESET)
!! 6. Initialize the new action object by invoking the type-bound procedure 'initialize'.
!!    (CALL act_obj%initialize(actionID)). The actionID defines the specific action to be 
!!    initialized. By this, you assign all matching fields to your particular action. 
!!    I.e. this is the reverse operation of assigning actions to fields as done in step 2.
!! 7. Execute your newly defined action object at a suitable place by invoking the 
!!    type-bound procedure 'execute' (CALL act_obj%execute(slack)). 'Slack' is the user-defined 
!!    maximum allowed time mismatch for executing the action.  
!!    
!!
!! @author Daniel Reinert, DWD
!!
!!
!! @par Revision History
!! Initial revision by Daniel Reinert, DWD (2014-01-09)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_action

  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: message, message_text
  USE mo_impl_constants,     ONLY: vname_len
  USE mtime,                 ONLY: event, newEvent, datetime, newDatetime,    &
    &                              isCurrentEventActive, deallocateDatetime,  &
    &                              MAX_DATETIME_STR_LEN, PROLEPTIC_GREGORIAN, &
    &                              MAX_EVENTNAME_STR_LEN, timedelta,          &
    &                              newTimedelta, deallocateTimedelta
  USE mo_mtime_extensions,   ONLY: get_datetime_string, getPTStringFromMS,    &
    &                              getTriggeredPreviousEventAtDateTime
  USE mo_util_string,        ONLY: remove_duplicates
  USE mo_run_config,         ONLY: msg_level
  USE mo_action_types,       ONLY: t_var_action
  USE mo_time_config,        ONLY: time_config
  USE mo_var_list,           ONLY: nvar_lists, var_lists
  USE mo_linked_list,        ONLY: t_list_element
  USE mo_var_list_element,   ONLY: t_var_list_element


  IMPLICIT NONE

  PRIVATE

  ! VARIABLES/OBJECTS
  PUBLIC :: reset_act

  !!!! temporary workaround for gfortran 4.5 and potentially others !!!!!!
  PUBLIC :: action_init    ! wrapper for CALL reset_act%initialize
  PUBLIC :: reset_action   ! wrapper for CALL reset_act%execute


  INTEGER, PARAMETER :: NMAX_VARS = 50  ! maximum number of fields that can be 
                                        ! assigned to a single action


  ! List of available actions
  !
  INTEGER, PARAMETER, PUBLIC :: ACTION_RESET = 1   ! re-set field to 0


  ! type for generating an array of pointers of type t_var_list_element
  !
  TYPE t_var_element_ptr
    TYPE(t_var_list_element), POINTER :: p
    TYPE(event)             , POINTER :: event     ! event from mtime library
  END TYPE t_var_element_ptr


  ! base type for action objects
  !
  TYPE, abstract:: t_action_obj
    INTEGER                    :: actionID                    ! action ID
    TYPE(t_var_element_ptr)    :: var_element_ptr(NMAX_VARS)  ! assigned variables
    INTEGER                    :: var_action_index(NMAX_VARS) ! index in var_element_ptr(10)%action

    INTEGER                    :: nvars                 ! number of variables for which 
                                                        ! this action is to be performed
  CONTAINS

    PROCEDURE :: initialize => action_collect_vars  ! initialize action object
    PROCEDURE :: execute    => action_execute       ! execute action object
    !
    ! deferred routine for action specific kernel (to be defined in extended type)
    PROCEDURE(kernel), deferred :: kernel
  END TYPE t_action_obj
  !
  ! kernel interface
  abstract INTERFACE
    SUBROUTINE kernel(act_obj, ivar)
      IMPORT                    :: t_action_obj
      CLASS(t_action_obj)       :: act_obj
      INTEGER, INTENT(IN)       :: ivar
    END SUBROUTINE kernel
  END INTERFACE



  ! extension of the action base type for the purpose of creating objects of that type.
  !
  ! create specific type for reset-action
  !
  TYPE, extends(t_action_obj) :: t_reset_obj
  CONTAINS
    PROCEDURE :: kernel => reset_kernel     ! type-specific action kernel (to be defined by user)
  END TYPE t_reset_obj


  ! create action object
  !
  TYPE(t_reset_obj) :: reset_act  ! action which resets field to resetval%rval

CONTAINS

  !> 
  !! Initialize action object
  !!
  !! Assign variables to specific actions/initialize action-object.
  !! When generating a new field via add_var, it is possible to assign various 
  !! actions to this field. Here, we go the other way around. For a specific 
  !! action, we loop over all fields and check, whether this action must 
  !! be performed for the particular field. If this is the case, this field 
  !! is assigned to the action.
  !! The action to be initialized is identified via the actionID.
  !!
  !! Loop over all variables and collect the variables names
  !! corresponding to the action @p action%actionID
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2014-01-13)
  !!
  SUBROUTINE action_collect_vars(act_obj, actionID) 
    CLASS(t_action_obj)         :: act_obj
    INTEGER      , INTENT(IN)   :: actionID

    ! local variables
    INTEGER :: i, iact
    INTEGER :: nvars                        ! number of variables assigned to action

    TYPE(t_list_element), POINTER :: element
    TYPE(t_var_action)  , POINTER :: action_list
    CHARACTER(LEN=2)                    :: str_actionID
    CHARACTER(LEN=128)                  :: intvl             ! action interval [PTnH]
    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: sim_start, sim_end
    CHARACTER(LEN=MAX_EVENTNAME_STR_LEN):: event_name
    CHARACTER(LEN=vname_len)            :: varlist(NMAX_VARS)
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

          END IF
        ENDDO  LOOPACTION ! loop over variable-specific actions

        IF(ASSOCIATED(action_list)) action_list => NULL()

      ENDDO LOOPVAR ! loop over vlist "i"
    ENDDO ! i = 1,nvar_lists

    ! set nvars
    act_obj%nvars = nvars


    IF (msg_level >= 11) THEN
      ! remove duplicate variable names
      DO i=1,act_obj%nvars
        varlist(i) = TRIM(act_obj%var_element_ptr(i)%p%info%name)
      ENDDO
      CALL remove_duplicates(varlist,nvars)

      WRITE(message_text,'(a,i0,a)') 'Variables assigned to action ',act_obj%actionID,':'
      DO i=1, nvars
        IF (i==1) THEN
          WRITE(message_text,'(a,a,a)') TRIM(message_text), " ",  TRIM(varlist(i))
        ELSE
          WRITE(message_text,'(a,a,a)') TRIM(message_text), ", ", TRIM(varlist(i))
        END IF
      ENDDO
      CALL message('',message_text)
    ENDIF

  END SUBROUTINE action_collect_vars




  !>
  !! Execute action
  !!
  !! For each field attached to this action it is checked, whether the action should 
  !! be executed at the datetime given. This routine does not make any assumption 
  !! about the details of the action to be executed. The action itself is encapsulated 
  !! in the kernel-routine.
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2014-09-11)
  !!
  SUBROUTINE action_execute(act_obj, slack)
    !
    CLASS(t_action_obj)       :: act_obj
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

    LOGICAL :: isactive

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: mtime_cur_datetime
    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: str_slack       ! slack as string

    TYPE(timedelta), POINTER :: p_slack                    ! slack in 'timedelta'-Format

    CHARACTER(*), PARAMETER :: routine = TRIM("mo_action:execute")

    TYPE(datetime) :: lastTrigger_datetime  ! latest intended triggering date

  !-------------------------------------------------------------------------

    ! compute current datetime in a format appropriate for mtime
    CALL get_datetime_string(mtime_cur_datetime, time_config%cur_datetime)
    mtime_date  => newDatetime(TRIM(mtime_cur_datetime)) 

    ! compute allowed slack in PT-Format
    ! Use factor 999 instead of 1000, since no open interval is available
    ! needed [trigger_date, trigger_date + slack[
    ! used   [trigger_date, trigger_date + slack]
    CALL getPTStringFromMS(INT(999._wp*slack),str_slack)
    ! get slack in 'timedelta'-format appropriate for isCurrentEventActive
    p_slack => newTimedelta(str_slack)


    ! Loop over all fields attached to this action
    DO ivar = 1, act_obj%nvars
      var_action_idx = act_obj%var_action_index(ivar)

      field      => act_obj%var_element_ptr(ivar)%p
      this_event => act_obj%var_element_ptr(ivar)%event


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
      isactive = LOGICAL(isCurrentEventActive(this_event,mtime_date, plus_slack=p_slack))



      ! Check wheter the action should be triggered for variable
      ! under consideration
      IF (isactive) THEN

        ! store latest true triggering date
        field%info%action_list%action(var_action_idx)%lastActive = TRIM(mtime_cur_datetime)
        ! store latest intended triggering date
        CALL getTriggeredPreviousEventAtDateTime(this_event, lastTrigger_datetime)
        field%info%action_list%action(var_action_idx)%EventLastTriggerDate = lastTrigger_datetime 


        IF (msg_level >= 12) THEN
          WRITE(message_text,'(a,i2,a,a,a,a)') 'action ',act_obj%actionID, &
            &  ' triggered for ', TRIM(field%info%name), ' at ', TRIM(mtime_cur_datetime)
          CALL message(TRIM(routine),message_text)
        ENDIF

        ! perform action on element number 'ivar'
        CALL act_obj%kernel(ivar)

      ENDIF

    ENDDO


    ! cleanup
    !
    CALL deallocateDatetime(mtime_date)
    CALL deallocateTimedelta(p_slack)

  END SUBROUTINE action_execute


  !>
  !! Reset-action kernel
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2014-09-12)
  !!
  SUBROUTINE reset_kernel(act_obj, ivar)
    ! usually, we have "t_reset_obj" as PASS type for this deferred
    ! subroutine, however, the PGI 12.3 expects the base class type...
#if defined (__PGI) 
    CLASS (t_action_obj) :: act_obj
#else
    CLASS (t_reset_obj)  :: act_obj
#endif
    INTEGER, INTENT(IN) :: ivar    ! element number 

    ! re-set field to its pre-defined reset-value
    act_obj%var_element_ptr(ivar)%p%r_ptr = act_obj%var_element_ptr(ivar)%p%info%resetval%rval

  END SUBROUTINE reset_kernel


  !=================================================================================!
  !       WORKAROUND for GFORTRAN 4.5 and potentially other ancient compilers       !
  !=================================================================================!
  !
  ! wrapper for reset_act%initialize
  !
  SUBROUTINE action_init(actionID)
    INTEGER, INTENT(IN) :: actionID

    ! Initialize reset-Action, i.e. assign variables to action object
    CALL reset_act%initialize(actionID)
  END SUBROUTINE action_init

  !
  ! wrapper for reset_act%execute
  !
  SUBROUTINE reset_action(slack)
    REAL(wp), INTENT(IN) :: slack

    ! execute reset action
    CALL reset_act%execute(slack)
  END SUBROUTINE reset_action

END MODULE mo_action

