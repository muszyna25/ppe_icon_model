!>
!! @author G. Zaengl
!!
!! @par Revision History
!! Moved configure state from namelists/mo_prepicon_nml:
!! F. Prill, DWD (2012-01-31)
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
MODULE mo_prepicon_config


  IMPLICIT NONE

  PUBLIC :: i_oper_mode


  CHARACTER(len=*),PARAMETER,PRIVATE :: &
    &  version = '$Id$'

  ! ----------------------------------------------------------------------------
  ! 1.0 Namelist variables for the prep_icon preprocessing program
  ! ----------------------------------------------------------------------------
  !
  INTEGER  :: i_oper_mode   ! operation mode


END MODULE mo_prepicon_config
