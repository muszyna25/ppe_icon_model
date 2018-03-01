!>
!! Namelist for configuring turbulent mixing parameterization
!!
!! @par Revision History
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_vdiff_nml

  USE mo_vdiff_config,        ONLY: vdiff_config
  USE mo_io_units,            ONLY: nnml
  USE mo_namelist,            ONLY: position_nml, POSITIONED, open_nml, close_nml
  USE mo_master_control,      ONLY: use_restart_namelists
  USE mo_restart_namelist,    ONLY: open_tmpfile, store_and_close_namelist,  &
                                  & open_and_restore_namelist, close_tmpfile
  USE mo_mpi,                 ONLY: my_process_is_stdio
  USE mo_nml_annotate,        ONLY: temp_defaults, temp_settings

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: read_vdiff_namelist

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
    INTEGER :: iunit

    !----------------------------------------------------------------
    ! Default values
    !----------------------------------------------------------------
    lsfc_mom_flux  = .TRUE.
    lsfc_heat_flux = .TRUE.

    !----------------------------------------------------------------
    ! If this is a resumed integration, overwrite the defaults above 
    ! by values in the previous integration.
    !----------------------------------------------------------------
    IF (use_restart_namelists()) THEN
      funit = open_and_restore_namelist('vdiff_nml')
      READ(funit,NML=vdiff_nml)
      CALL close_tmpfile(funit)
    END IF

    !---------------------------------------------------------------------
    ! Read user's (new) specifications (Done so far by all MPI processes)
    !---------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml('vdiff_nml',STATUS=ist)
    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, vdiff_nml)  ! write defaults to temporary text file
    END IF
    SELECT CASE (ist)
    CASE (POSITIONED)
      READ (nnml, vdiff_nml)                                      ! overwrite default settings
      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, vdiff_nml)  ! write settings to temporary text file
      END IF
    END SELECT
    CALL close_nml

    !-----------------------------------------------------
    ! Store the namelist for restart
    !-----------------------------------------------------
    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=vdiff_nml)
      CALL store_and_close_namelist(funit, 'vdiff_nml')
    ENDIF

    !-----------------------------------------------------
    ! Fill the configuration state
    !-----------------------------------------------------
    vdiff_config%lsfc_mom_flux  = lsfc_mom_flux 
    vdiff_config%lsfc_heat_flux = lsfc_heat_flux 

  END SUBROUTINE read_vdiff_namelist

END MODULE mo_vdiff_nml
