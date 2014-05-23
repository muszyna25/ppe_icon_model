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

  USE mtime,                 ONLY: MAX_DATETIME_STR_LEN

  IMPLICIT NONE

  PRIVATE

  ! TYPES
  PUBLIC  :: t_var_action_element
  PUBLIC  :: t_var_action


  CHARACTER(len=*), PARAMETER :: &
    &  version = '$Id$'

  INTEGER, PARAMETER :: NMAX_ACTION=5

  ! defines single variable specific action
  !
  TYPE t_var_action_element
    INTEGER                         :: actionID              ! action ID  
    CHARACTER(LEN=128)              :: intvl                 ! action interval [PTnH]
    CHARACTER(MAX_DATETIME_STR_LEN) :: lastActive            ! date of last triggering
  END TYPE t_var_action_element


  ! defines list of variable specific actions
  !
  TYPE t_var_action
    TYPE(t_var_action_element) :: action(NMAX_ACTION)
    INTEGER                    :: n_actions           ! number of variable specific actions  
  END TYPE t_var_action


END MODULE mo_action_types

