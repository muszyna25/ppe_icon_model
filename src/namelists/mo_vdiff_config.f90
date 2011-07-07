!>
!! Configuration of the Brinkop and Roeckner (1995) turbulent mixing scheme.
!!
!! @author Hui Wan, MPI-M
!!
!! @par Revision History
!! First version by Hui Wan, MPI (2010-07)
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
MODULE mo_vdiff_config

  IMPLICIT NONE
  PUBLIC
  CHARACTER(len=*), PARAMETER, PRIVATE :: version = '$Id$'

  !>
  !! Derived type
  !!
  TYPE t_vdiff_config

    LOGICAL :: lsfc_mom_flux   !< switch on/off surface momentum flux
    LOGICAL :: lsfc_heat_flux  !< switch on/off surface heat flux                                  
                               !< (sensible AND latent)

  END TYPE t_vdiff_config

  !>
  !! The configuration state (variable).
  !!
  TYPE(t_vdiff_config) :: vdiff_config

END MODULE mo_vdiff_config
