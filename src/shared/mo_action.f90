!>
!! Routines for defining, initializing and executing action events
!!
!! Routines for defining, initializing and executing action events
!!
!! *****************************************************************
!!            RECIPE FOR CREATING A NEW ACTION EVENT
!! *****************************************************************
!! 1. Define a new actionTyp in mo_action
!! 2. Assign new actionTyp to variables of your choice.
!!    E.g. add the following code snippet to an add_var/add_ref of your choice:
!!    action_list=actions(new_action(ACTION_XXX,'PTXXH'), new_action(...), ...)
!!    ACTION_XXX is the actionTyp defined in step 1, and PTXXH is the
!!    interval at which the action should be triggered.
!! 3. Create an extension of the abstract type t_action_obj and overwrite
!!    the deferred procedure 'kernel' with your action-specific kernel-routine
!!    (to be defined in step 5).
!! 4. Create a variable (object) of the type defined in step 3.
!! 5. Write your own action-Routine (action-kernel). This is the routine which actually does
!!    the work. (see e.g. routine 'reset_kernel' for actionTyp=ACTION_RESET)
!! 6. Initialize the new action object by invoking the type-bound procedure 'initialize'.
!!    (CALL act_obj%initialize(actionTyp)). The actiontyp defines the specific action to be
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

  USE mo_kind,               ONLY: wp, i8
  USE mo_mpi,                ONLY: my_process_is_stdio, p_pe, p_io, p_comm_work, p_bcast
  USE mo_exception,          ONLY: message, message_text, finish
  USE mo_impl_constants,     ONLY: vname_len
  USE mtime,                 ONLY: event, newEvent, datetime,                        &
    &                              isCurrentEventActive, deallocateDatetime,         &
    &                              MAX_DATETIME_STR_LEN,                             &
    &                              MAX_EVENTNAME_STR_LEN, timedelta,                 &
    &                              newTimedelta, deallocateTimedelta,                &
    &                              getTriggeredPreviousEventAtDateTime,              &
    &                              getPTStringFromMS, datetimetostring
  USE mo_util_mtime,         ONLY: is_event_active
  USE mo_util_string,        ONLY: remove_duplicates
  USE mo_util_table,         ONLY: initialize_table, finalize_table, add_table_column, &
    &                              set_table_entry, print_table, t_table
  USE mo_action_types,       ONLY: t_var_action, action_names, action_reset
  USE mo_grid_config,        ONLY: n_dom
  USE mo_run_config,         ONLY: msg_level
  USE mo_parallel_config,    ONLY: proc0_offloading
  USE mo_var_list,           ONLY: var_lists_apply
  USE mo_linked_list,        ONLY: t_var_list_intrinsic
  USE mo_var_list_element,   ONLY: t_var_list_element
  USE mo_var_metadata_types, ONLY: t_var_metadata
  USE mo_fortran_tools,      ONLY: init

  IMPLICIT NONE

  PRIVATE

  ! VARIABLES/OBJECTS
  PUBLIC :: reset_act

  ! PARAMETER
  PUBLIC :: action_names, action_reset

  INTEGER, PARAMETER :: NMAX_VARS = 100 ! maximum number of fields that can be
                                        ! assigned to a single action


  ! type for generating an array of pointers of type t_var_list_element
  !
  TYPE t_var_element_ptr
    TYPE(t_var_list_element), POINTER :: p
    TYPE(event)             , POINTER :: event     ! event from mtime library
    INTEGER                           :: patch_id  ! patch on which field lives
  END TYPE t_var_element_ptr


  ! base type for action objects
  !
  TYPE, abstract:: t_action_obj
    INTEGER                    :: actionTyp                   ! Type of action
    TYPE(t_var_element_ptr)    :: var_element_ptr(NMAX_VARS)  ! assigned variables
    INTEGER                    :: var_action_index(NMAX_VARS) ! index in var_element_ptr(10)%action

    INTEGER                    :: nvars                 ! number of variables for which
                                                        ! this action is to be performed
  CONTAINS

    PROCEDURE :: initialize => action_collect_vars  ! initialize action object
#if defined (__SX__) || defined (__NEC_VH__)
    PROCEDURE :: execute    => action_execute_SX    ! execute action object
#else
    PROCEDURE :: execute    => action_execute       ! execute action object
#endif
    PROCEDURE :: print_setup=> action_print_setup   ! Screen print out of action object setup
    !
    ! deferred routine for action specific kernel (to be defined in extended type)
    PROCEDURE(kernel), deferred :: kernel
  END TYPE t_action_obj
  !
  ! kernel interface
  ABSTRACT INTERFACE
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
  !! The action to be initialized is identified via the actionTyp.
  !!
  !! Loop over all variables and collect the variables names
  !! corresponding to the action @p action%actionTyp
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2014-01-13)
  !!
  SUBROUTINE action_collect_vars(act_obj, actionTyp)
    CLASS(t_action_obj)         :: act_obj
    INTEGER      , INTENT(IN)   :: actionTyp

    ! local variables
    INTEGER :: i, l, tlen, sep_len
    INTEGER :: nvars                        ! number of variables assigned to action

    TYPE(t_var_action)  , POINTER       :: action_list
    CHARACTER(LEN=vname_len)            :: varlist(NMAX_VARS)
    CHARACTER(len=2), PARAMETER :: sep = ', '
  !-------------------------------------------------------------------------

    ! init nvars
    act_obj%nvars = 0

    ! store actionTyp
    act_obj%actionTyp = actionTyp

    ! loop over all variable lists and variables
    !
    CALL var_lists_apply(action_collect_var, act_obj)


    IF (msg_level >= 11) THEN

      ! remove duplicate variable names
      nvars = act_obj%nvars
      DO i=1, nvars
        varlist(i) = act_obj%var_element_ptr(i)%p%info%name
      ENDDO
      CALL remove_duplicates(varlist, nvars)

      WRITE(message_text,'(3a)') 'Variables assigned to action ',TRIM(ACTION_NAMES(actionTyp)),':'
      tlen = LEN_TRIM(message_text)
      DO i = 1, nvars
        sep_len = MERGE(1,2,i==1)
        l = LEN_TRIM(varlist(i))
        WRITE(message_text(tlen+1:),'(2a)') sep(3-sep_len:2), varlist(i)(:l)
        tlen = tlen + sep_len + l
      ENDDO
      CALL message('',message_text)

      IF(my_process_is_stdio()) THEN
        CALL act_obj%print_setup()
      ENDIF
    ENDIF

  END SUBROUTINE action_collect_vars

  !> wrapper, see action_collect_var_ for more information
  SUBROUTINE action_collect_var(field, state, var_list)
    TYPE(t_var_list_element), TARGET :: field
    CLASS(*), TARGET :: state
    TYPE(t_var_list_intrinsic), INTENT(in) :: var_list

    SELECT TYPE(state)
    CLASS IS (t_action_obj)
      CALL action_collect_var_(field, state, var_list%patch_id)
    END SELECT
  END SUBROUTINE action_collect_var

  !> inspect single field for applicable actions, if so, transfer
  !! action_list items into act_obj fields in case of matching
  !! actionTyp member
  SUBROUTINE action_collect_var_(field, act_obj, patch_id)
    TYPE(t_var_list_element), TARGET, INTENT(in) :: field
    CLASS(t_action_obj), INTENT(inout) :: act_obj
    INTEGER, INTENT(in) :: patch_id
    INTEGER :: nvars, i, action_type
    CHARACTER(LEN=MAX_EVENTNAME_STR_LEN):: event_name

    nvars = act_obj%nvars
    action_type = act_obj%actionTyp
    ! point to variable specific action list
    ASSOCIATE(action_list => field%info%action_list)
      !
      ! Loop over all variable-specific actions
      !
      LOOPACTION: DO i = 1,action_list%n_actions

        ! If the variable specific action fits, assign variable
        ! to the corresponding action.
        IF (action_list%action(i)%actionTyp == action_type) THEN

          ! Add field to action object
          nvars = nvars + 1
          act_obj%var_element_ptr(nvars)%p => field
          act_obj%var_action_index(nvars) = i
          act_obj%var_element_ptr(nvars)%patch_id = patch_id


          ! Create event for this specific field
          WRITE(event_name,'(a,i2,2a)') 'act_TYP', action_type, &
               '_', TRIM(action_list%action(i)%intvl)
          act_obj%var_element_ptr(nvars)%event =>newEvent(                   &
            &                                 event_name,                    &
            &                                 action_list%action(i)%ref,  &
            &                                 action_list%action(i)%start,&
            &                                 action_list%action(i)%end  ,&
            &                                 action_list%action(i)%intvl)

        END IF
      ENDDO  LOOPACTION ! loop over variable-specific actions
    END ASSOCIATE
    act_obj%nvars = nvars
  END SUBROUTINE action_collect_var_

  !>
  !! Screen print out of action event setup
  !!
  !! Screen print out of action event setup.
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2015-01-06)
  !!
  SUBROUTINE action_print_setup (act_obj)

    CLASS(t_action_obj)  :: act_obj  !< action for which setup will be printed

    ! local variables
    TYPE(t_table)   :: table
    INTEGER         :: ivar            ! loop counter
    INTEGER         :: irow            ! row to fill
    INTEGER         :: var_action_idx  ! index of current action in variable-specific action list
    INTEGER         :: jg              ! patch loop counter
    CHARACTER(LEN=2):: str_patch_id
    !--------------------------------------------------------------------------

    ! table-based output
    CALL initialize_table(table)
    ! the latter is no longer mandatory
    CALL add_table_column(table, "VarName")
    CALL add_table_column(table, "PID")
    CALL add_table_column(table, "Ref date")
    CALL add_table_column(table, "Start date")
    CALL add_table_column(table, "End date")
    CALL add_table_column(table, "Interval")

    irow = 0
    ! print event info sorted by patch ID in ascending order
    DO jg = 1, n_dom
      DO ivar=1,act_obj%nvars

        IF (act_obj%var_element_ptr(ivar)%patch_id /= jg) CYCLE

        var_action_idx = act_obj%var_action_index(ivar)

        irow = irow + 1
        CALL set_table_entry(table,irow,"VarName", TRIM(act_obj%var_element_ptr(ivar)%p%info%name))
        write(str_patch_id,'(i2)')  act_obj%var_element_ptr(ivar)%patch_id
        CALL set_table_entry(table,irow,"PID", TRIM(str_patch_id))
        CALL set_table_entry(table,irow,"Ref date", &
          &  TRIM(act_obj%var_element_ptr(ivar)%p%info%action_list%action(var_action_idx)%ref))
        CALL set_table_entry(table,irow,"Start date", &
          &  TRIM(act_obj%var_element_ptr(ivar)%p%info%action_list%action(var_action_idx)%start))
        CALL set_table_entry(table,irow,"End date", &
          &  TRIM(act_obj%var_element_ptr(ivar)%p%info%action_list%action(var_action_idx)%end))
        CALL set_table_entry(table,irow,"Interval", &
          &  TRIM(act_obj%var_element_ptr(ivar)%p%info%action_list%action(var_action_idx)%intvl))
      ENDDO
    ENDDO  ! jg

    CALL print_table(table, opt_delimiter=' | ')
    CALL finalize_table(table)

    WRITE (0,*) " " ! newline
  END SUBROUTINE action_print_setup



  !>
  !! Execute action
  !!
  !! For each variable attached to this action it is checked, whether the action should
  !! be executed at the datetime given. This routine does not make any assumption
  !! about the details of the action to be executed. The action itself is encapsulated
  !! in the kernel-routine.
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2014-09-11)
  !!
  SUBROUTINE action_execute(act_obj, slack, mtime_date)
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

    CHARACTER(*), PARAMETER :: routine = TRIM("mo_action:")

  !-------------------------------------------------------------------------

    ! compute allowed slack in PT-Format
    ! Use factor 999 instead of 1000, since no open interval is available
    ! needed [trigger_date, trigger_date + slack[
    ! used   [trigger_date, trigger_date + slack]
    CALL getPTStringFromMS(INT(999.0_wp*slack,i8),str_slack)
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
      isactive = isCurrentEventActive(this_event,mtime_date, plus_slack=p_slack)

      ! Check whether the action should be triggered for variable
      ! under consideration
      IF (isactive) THEN

        ! store latest true triggering date
        CALL datetimeToString(mtime_date, mtime_cur_datetime)
        field%info%action_list%action(var_action_idx)%lastActive = TRIM(mtime_cur_datetime)
        ! store latest intended triggering date
        CALL getTriggeredPreviousEventAtDateTime(this_event, &
          field%info%action_list%action(var_action_idx)%EventLastTriggerDate)

        IF (msg_level >= 12) THEN
          WRITE(message_text,'(5a,i2,a,a)') 'action ',TRIM(ACTION_NAMES(act_obj%actionTyp)),&
            &  ' triggered for ', TRIM(field%info%name),' (PID ',               &
            &  act_obj%var_element_ptr(ivar)%patch_id,') at ',                  &
            &  TRIM(mtime_cur_datetime)
          CALL message(TRIM(routine),message_text)
        ENDIF

        ! perform action on element number 'ivar'
        CALL act_obj%kernel(ivar)

      ENDIF

    ENDDO


    ! cleanup
    !
    CALL deallocateTimedelta(p_slack)

  END SUBROUTINE action_execute


  !>
  !! Execute action
  !!
  !! !!! 'optimized' variant for SX-architecture                                   !!!
  !! !!! isCurrentEventActive is only executed by PE0 and broadcasted to the other !!!
  !! !!! workers on the vector engine.                                             !!!
  !!
  !! For each variable attached to this action it is checked, whether the action should
  !! be executed at the datetime given. This routine does not make any assumption
  !! about the details of the action to be executed. The action itself is encapsulated
  !! in the kernel-routine.
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2014-09-11)
  !!
  SUBROUTINE action_execute_SX(act_obj, slack, mtime_date)
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

    LOGICAL :: isactive(NMAX_VARS)

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: mtime_cur_datetime
    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: str_slack       ! slack as string

    TYPE(timedelta), POINTER :: p_slack                    ! slack in 'timedelta'-Format

    CHARACTER(*), PARAMETER :: routine = TRIM("mo_action:")

  !-------------------------------------------------------------------------

    ! init
    isactive(:) =.FALSE.

    ! compute allowed slack in PT-Format
    ! Use factor 999 instead of 1000, since no open interval is available
    ! needed [trigger_date, trigger_date + slack[
    ! used   [trigger_date, trigger_date + slack]
    IF (.NOT. proc0_offloading .OR. p_pe == p_io) THEN
      CALL getPTStringFromMS(INT(999.0_wp*slack,i8),str_slack)
      ! get slack in 'timedelta'-format appropriate for isCurrentEventActive
      p_slack => newTimedelta(str_slack)

      ! Check for all action events whether they are due at the datetime given.
      DO ivar = 1, act_obj%nvars

        var_action_idx = act_obj%var_action_index(ivar)

        this_event => act_obj%var_element_ptr(ivar)%event

        ! Check whether event-pointer is associated.
        IF (.NOT. ASSOCIATED(this_event)) THEN
          field => act_obj%var_element_ptr(ivar)%p
          WRITE (message_text,'(a,i2,a,a,a)')                           &
               'WARNING: action event ', var_action_idx, ' of field ',  &
               TRIM(field%info%name),': Event-Ptr is disassociated!'
          CALL message(routine,message_text)
        ENDIF

        ! Note that a second call to isCurrentEventActive will lead to
        ! a different result! Is this a bug or a feature?
        ! triggers in interval [trigger_date + slack]
        isactive(ivar) = is_event_active(this_event, mtime_date, proc0_offloading, plus_slack=p_slack, opt_lasync=.TRUE.)
      ENDDO

      ! cleanup
      !
      CALL deallocateTimedelta(p_slack)
    ENDIF

    ! broadcast
    IF (proc0_offloading) CALL p_bcast(isactive, p_io, p_comm_work)


    ! Loop over all fields attached to this action
    DO ivar = 1, act_obj%nvars

      ! Check whether the action should be triggered for variable
      ! under consideration
      IF (isactive(ivar)) THEN

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


        ! store latest true triggering date
        CALL datetimeToString(mtime_date, mtime_cur_datetime)
        field%info%action_list%action(var_action_idx)%lastActive = TRIM(mtime_cur_datetime)
        ! store latest intended triggering date
        CALL getTriggeredPreviousEventAtDateTime(this_event, &
          field%info%action_list%action(var_action_idx)%EventLastTriggerDate)

        IF (msg_level >= 12) THEN
          WRITE(message_text,'(5a,i2,a,a)') 'action ',TRIM(ACTION_NAMES(act_obj%actionTyp)),&
            &  ' triggered for ', TRIM(field%info%name),' (PID ',               &
            &  act_obj%var_element_ptr(ivar)%patch_id,') at ',                  &
            &  TRIM(mtime_cur_datetime)
          CALL message(TRIM(routine),message_text)
        ENDIF

        ! perform action on element number 'ivar'
        CALL act_obj%kernel(ivar)

      ENDIF

    ENDDO


  END SUBROUTINE action_execute_SX


  !>
  !! Reset-action kernel
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2014-09-12)
  !! Modification by Daniel Reinert, DWD (2014-12-17)
  !! - extend reset kernel to integer fields
  !!
  SUBROUTINE reset_kernel(act_obj, ivar)
    CLASS (t_reset_obj)  :: act_obj
    INTEGER, INTENT(IN) :: ivar    ! element number

    CHARACTER(len=*), PARAMETER :: routine = 'mo_action:reset_kernel'
    !-------------------------------------------------------------------

    ! re-set field to its pre-defined reset-value
    IF (ASSOCIATED(act_obj%var_element_ptr(ivar)%p%r_ptr)) THEN
!$OMP PARALLEL
      CALL init(act_obj%var_element_ptr(ivar)%p%r_ptr, &
           act_obj%var_element_ptr(ivar)%p%info%resetval%rval)
!$OMP END PARALLEL
    ELSE IF (ASSOCIATED(act_obj%var_element_ptr(ivar)%p%i_ptr)) THEN
!$OMP PARALLEL
      CALL init(act_obj%var_element_ptr(ivar)%p%i_ptr, &
           act_obj%var_element_ptr(ivar)%p%info%resetval%ival)
!$OMP END PARALLEL
    ELSE
      CALL finish (routine, 'Field not allocated for '//TRIM(act_obj%var_element_ptr(ivar)%p%info%name))
    ENDIF

  END SUBROUTINE reset_kernel

END MODULE mo_action

