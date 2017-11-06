!>
!! Configuration of the Brinkop and Roeckner (1995) turbulent mixing scheme.
!!
!! @author Hui Wan, MPI-M
!!
!! @par Revision History
!! First version by Hui Wan, MPI (2010-07)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_vdiff_config

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_vdiff_config, vdiff_config

  !>
  !! Derived type
  !!
  TYPE t_vdiff_config

    LOGICAL :: lsfc_mom_flux   !< switch on/off surface momentum flux
    LOGICAL :: lsfc_heat_flux  !< switch on/off surface heat flux
                               !< (sensible AND latent)

    LOGICAL :: lsfc_co2_flux   !< switch for coupled co2

  END TYPE t_vdiff_config

  !>
  !! The configuration state (variable).
  !!
  TYPE(t_vdiff_config) :: vdiff_config

END MODULE mo_vdiff_config
