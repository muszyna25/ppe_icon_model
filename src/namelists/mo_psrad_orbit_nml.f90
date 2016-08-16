!>
!! Namelist for orbit parameters of orbit used with psrad radiation
!!
!! @par Revision History
!!
!! Created by S. Rast (2016-01-28)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_psrad_orbit_nml

  USE mo_kind,                ONLY: wp
  USE mo_psrad_orbit_config,  ONLY: psrad_orbit_config
  USE mo_io_units,            ONLY: nnml
  USE mo_namelist,            ONLY: position_nml, POSITIONED, open_nml, close_nml
  USE mo_master_config,       ONLY: isRestart
  USE mo_io_restart_namelist, ONLY: open_tmpfile, store_and_close_namelist,  &
                                  & open_and_restore_namelist, close_tmpfile
  USE mo_mpi,                 ONLY: my_process_is_stdio
  USE mo_nml_annotate,        ONLY: temp_defaults, temp_settings

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: read_psrad_orbit_namelist

CONTAINS
  !>
  !!
  SUBROUTINE read_psrad_orbit_namelist( filename )

    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER :: ist, funit
    INTEGER :: iunit

    !--------------------
    ! namelist variables   
    !--------------------

    REAL(wp) :: cecc           !< Eccentricity of Earth's Orbit
    REAL(wp) :: cobld          !< Obliquity of Earth [Deg]
    LOGICAL  :: l_orbvsop87    !< .TRUE. for VSOP87 orbit, 
                               !< .FALSE. for Kepler orbit
    LOGICAL  :: l_sph_symm_irr !< .TRUE. for globally averaged irradiation (RCE)
                               !< .FALSE. for lat (lon) dependent irradiation

    NAMELIST /psrad_orbit_nml/ cecc, cobld, l_orbvsop87, l_sph_symm_irr

    !----------------------------------------------------------------
    ! Default values
    !----------------------------------------------------------------
      cecc        =  0.016715_wp !< Eccentricity of Earth's Orbit
      cobld       =  23.44100_wp !< Obliquity of Earth [Deg]
      l_orbvsop87 = .TRUE.       !< (real) observed orbit, not idealized
      l_sph_symm_irr = .FALSE.   !< lat (lon) dependent irradiation
    !----------------------------------------------------------------
    ! If this is a resumed integration, overwrite the defaults above 
    ! by values in the previous integration.
    !----------------------------------------------------------------
    IF (isRestart()) THEN
      funit = open_and_restore_namelist('psrad_orbit_nml')
      READ(funit,NML=psrad_orbit_nml)
      CALL close_tmpfile(funit)
    END IF

    !---------------------------------------------------------------------
    ! Read user's (new) specifications (Done so far by all MPI processes)
    !---------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml('psrad_orbit_nml',STATUS=ist)
    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, psrad_orbit_nml)  ! write defaults to temporary text file
    END IF
    SELECT CASE (ist)
    CASE (POSITIONED)
      READ (nnml, psrad_orbit_nml)                                      ! overwrite default settings
      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, psrad_orbit_nml)  ! write settings to temporary text file
      END IF
    END SELECT
    CALL close_nml

    !-----------------------------------------------------
    ! Store the namelist for restart
    !-----------------------------------------------------
    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=psrad_orbit_nml)
      CALL store_and_close_namelist(funit, 'psrad_orbit_nml')
    ENDIF

    !-----------------------------------------------------
    ! Fill the configuration state
    !-----------------------------------------------------
    psrad_orbit_config%cecc        = cecc
    psrad_orbit_config%cobld       = cobld
    psrad_orbit_config%l_orbvsop87 = l_orbvsop87
    psrad_orbit_config%l_sph_symm_irr = l_sph_symm_irr
  END SUBROUTINE read_psrad_orbit_namelist

END MODULE mo_psrad_orbit_nml
