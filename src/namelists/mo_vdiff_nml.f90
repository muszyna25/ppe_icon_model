!>
!! Namelist for configuring turbulent mixing parameterization
!!
!! @par Revision History
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
MODULE mo_vdiff_nml

  USE mo_vdiff_config,        ONLY: vdiff_config
  USE mo_io_units,            ONLY: nnml
  USE mo_namelist,            ONLY: position_nml, POSITIONED, open_nml, close_nml
  USE mo_master_nml,          ONLY: lrestart
  USE mo_io_restart_namelist, ONLY: open_tmpfile, store_and_close_namelist,  &
                                  & open_and_restore_namelist, close_tmpfile

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: read_vdiff_namelist

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  !--------------------
  ! namelist variables   
  !--------------------

  LOGICAL :: lsfc_mom_flux   !< switch on/off surface momentum flux
  LOGICAL :: lsfc_heat_flux  !< switch on/off surface heat flux
                             !< (sensible AND latent)

  NAMELIST /vdiff_nml/ lsfc_mom_flux, lsfc_heat_flux

CONTAINS
  !>
  !!
  SUBROUTINE read_vdiff_namelist( filename )

    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER :: ist, funit

    !----------------------------------------------------------------
    ! Default values
    !----------------------------------------------------------------
    lsfc_mom_flux  = .TRUE.
    lsfc_heat_flux = .TRUE.

    !----------------------------------------------------------------
    ! If this is a resumed integration, overwrite the defaults above 
    ! by values in the previous integration.
    !----------------------------------------------------------------
    IF (lrestart) THEN
      funit = open_and_restore_namelist('vdiff_nml')
      READ(funit,NML=vdiff_nml)
      CALL close_tmpfile(funit)
    END IF

    !---------------------------------------------------------------------
    ! Read user's (new) specifications (Done so far by all MPI processes)
    !---------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml('vdiff_nml',STATUS=ist)
    SELECT CASE (ist)
    CASE (POSITIONED)
      READ (nnml, vdiff_nml)
    END SELECT
    CALL close_nml

    !-----------------------------------------------------
    ! Store the namelist for restart
    !-----------------------------------------------------
    funit = open_tmpfile()
    WRITE(funit,NML=vdiff_nml)
    CALL store_and_close_namelist(funit, 'vdiff_nml')

    !-----------------------------------------------------
    ! Fill the configuration state
    !-----------------------------------------------------
    vdiff_config%lsfc_mom_flux  = lsfc_mom_flux 
    vdiff_config%lsfc_heat_flux = lsfc_heat_flux 

  END SUBROUTINE read_vdiff_namelist

END MODULE mo_vdiff_nml
