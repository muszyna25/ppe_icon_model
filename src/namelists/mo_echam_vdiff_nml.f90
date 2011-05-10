!>
!! Namelist for the configuration of the vertical diffusion
!!
!!
!! @par Revision History
!! Revision history in mo_echam_vdiff_params.f90 (r4300)
!! Modification by Constantin Junk, MPI-M (2011-05-05)
!! - moved echam_vdiff namelist variables and subroutine setup_vdiff
!!   from mo_echam_vdiff_params to namelists/mo_echam_vdiff_nml
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
MODULE mo_echam_vdiff_nml

  USE mo_io_units,            ONLY: nnml
  USE mo_exception,           ONLY: message, print_value
  USE mo_namelist,            ONLY: position_nml, POSITIONED

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: echam_vdiff_nml_setup
  PUBLIC :: echam_vdiff_ctl                             !< namelist

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  !--------------------------------------!
  ! echam_vdiff_nml namelist variables   !
  !--------------------------------------!

  LOGICAL,PUBLIC :: lsfc_mom_flux   !< switch on/off surface momentum flux
  LOGICAL,PUBLIC :: lsfc_heat_flux  !< switch on/off surface heat flux
                                    !< (sensible AND latent)

  NAMELIST/echam_vdiff_ctl/ lsfc_mom_flux, lsfc_heat_flux


CONTAINS
!-------------------------------------------------------------------------
!
!-------------------------------------------------------------------------
!
!
!>
!!   Set up vertical diffusion
!!
!!
!! @par Revision History
!!   Revision History in mo_echam_vdiff_params (r4300)
!!   Modification by Constantin Junk, MPI-M (2011-05-05)
!!   - renamed setup_vdiff to echam_vdiff_nml_setup
!!

  SUBROUTINE echam_vdiff_nml_setup

    INTEGER :: ist


    lsfc_mom_flux  = .TRUE.
    lsfc_heat_flux = .TRUE.

    ! Read namelist (every CPU does this)

    CALL position_nml('echam_vdiff_ctl',status=ist)
    SELECT CASE (ist)
    CASE (POSITIONED)
      READ (nnml, echam_vdiff_ctl)
    END SELECT

    ! Check validity; send values to stdout

    CALL message('','')
    CALL message('','------- namelist echam_vdiff_ctl --------')

    CALL print_value(' lsfc_mom_flux  ',lsfc_mom_flux)
    CALL print_value(' lsfc_heat_flux ',lsfc_heat_flux)

    CALL message('','---------------------------')
    CALL message('','')

  END SUBROUTINE echam_vdiff_nml_setup

END MODULE mo_echam_vdiff_nml
