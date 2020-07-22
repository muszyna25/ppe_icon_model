!>
!! Type definition for action events
!!
!! Type definition for action events
!!
!! @author Daniel Reinert, DWD
!!
!!
!! @par Revision History
!! Initial revision by daniel Reinert, DWD (2014-01-13)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_action_types

  USE mtime,                 ONLY: MAX_DATETIME_STR_LEN, datetime

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: actions ! TODO: disentagle varlists and usage of action stuff, so "actions" can move to mo_action...
  PUBLIC :: t_var_action_element, t_var_action

  INTEGER, PARAMETER :: NMAX_ACTION=5

  ! defines single variable specific action
  !
  TYPE t_var_action_element
    INTEGER                             :: actionTyp             ! Type of action
    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: start                 ! action start time (TimeStamp)
                                                                 ! [YYYY-MM-DDThh:mm:ss]
    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: end                   ! action end time   (TimeStamp)
                                                                 ! [YYYY-MM-DDThh:mm:ss]
    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: ref                   ! action reference time (TimeStamp)
                                                                 ! [YYYY-MM-DDThh:mm:ss]
    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: intvl                 ! action interval (duration)
                                                                 ! [PnYnMnDTnHnMn]
    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: lastActive            ! date of last triggering (TimeStamp)
                                                                 ! [YYYY-MM-DDThh:mm:ss]
    TYPE(datetime)                      :: EventLastTriggerDate  ! last intended trigger date.
                                                               ! Differs from lastActive in the sense 
                                                               ! that it is the intended trigger date, whereas 
                                                               ! lastActive is the TRUE trigger date. These two 
                                                               ! can differ by the allowed 'slack'.  
  END TYPE t_var_action_element

  ! defines list of variable specific actions
  !
  TYPE t_var_action
    TYPE(t_var_action_element) :: action(NMAX_ACTION)
    INTEGER                    :: n_actions           ! number of variable specific actions  
  END TYPE t_var_action

CONTAINS

  !>
  !! Generate list (array) of variable specific actions
  !!
  !! Generate list (array) of variable specific actions.
  !! Creates array 'action_list' of type t_var_action
  !
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2014-01-13)
  !!
  FUNCTION actions(a01, a02, a03, a04, a05)  RESULT(action_list)
    TYPE(t_var_action_element), INTENT(IN), OPTIONAL :: a01, a02, a03, a04, a05
    TYPE(t_var_action)             :: action_list
    INTEGER :: n_act             ! action counter
    ! create action list
    n_act = 0
    CALL add_action_item(a01)
    CALL add_action_item(a02)
    CALL add_action_item(a03)
    CALL add_action_item(a04)
    CALL add_action_item(a05)
    action_list%n_actions = n_act
  CONTAINS

    SUBROUTINE add_action_item(aX)
      TYPE(t_var_action_element), INTENT(IN), OPTIONAL :: aX

      IF (PRESENT(aX)) THEN
        n_act = n_act + 1
        action_list%action(n_act) = aX
      END IF
    END SUBROUTINE add_action_item
  END FUNCTION actions

END MODULE mo_action_types

