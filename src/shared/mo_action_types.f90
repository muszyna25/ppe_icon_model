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

