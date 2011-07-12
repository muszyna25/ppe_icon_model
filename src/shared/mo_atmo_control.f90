!>
!! Contains tha patch structures and the interpolation states for the atmosphere model
!!
!! @par Revision History
!! Leonidas Linardakis, MPI-M, March 2010
!!   Moved the data structures from the control_model.f90
!!
!! @par Copyright
!! 2002-2008 by DWD and MPI-M
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
MODULE mo_atmo_control

  USE mo_model_domain,        ONLY: t_patch

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*),    PARAMETER           :: version = '$Id$'

  TYPE(t_patch),         TARGET, ALLOCATABLE :: p_patch_global(:), p_patch_subdiv(:)
  TYPE(t_patch),         POINTER             :: p_patch(:)

  PUBLIC :: p_patch_global, p_patch_subdiv, p_patch

  !-------------------------------------------------------------------------

END MODULE mo_atmo_control
