!>
!! This module provides parameters controlling online emission module.
!!
!! Authors:
!! Michael Jähn, EMPA, michael.jaehn@empa.ch
!! Michael Steiner, EMPA, michael.steiner@empa.ch
!!
!!
!! Current Code Owner: C2SM, Anne Roches
!!  phone:  
!!  email:  anne.roches@env.ethz.ch
!!
!!
!!
!!
MODULE mo_oem_nml

    USE mo_oem_config, ONLY: config_vertical_profile_nc   => vertical_profile_nc,   &
                         &   config_hour_of_day_nc        => hour_of_day_nc,        &
                         &   config_day_of_week_nc        => day_of_week_nc,        &
                         &   config_month_of_year_nc      => month_of_year_nc,      &
                         &   config_hour_of_year_nc       => hour_of_year_nc,       &
                         &   config_gridded_emissions_nc  => gridded_emissions_nc

  USE mo_kind,               ONLY: wp
  USE mo_impl_constants,     ONLY: MAX_CHAR_LENGTH
  USE mo_mpi,                ONLY: my_process_is_stdio
  USE mo_namelist,           ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_io_units,           ONLY: nnml, nnml_output
  USE mo_master_control,     ONLY: use_restart_namelists
  USE mo_nml_annotate,       ONLY: temp_defaults, temp_settings
  USE mo_io_units,           ONLY: filename_max

  IMPLICIT NONE
  PRIVATE
  PUBLIC:: read_oemctrl_namelist

  !-----------------------------------
  ! namelist variables and parameters
  !-----------------------------------
  !
  CHARACTER(LEN=filename_max) :: &
    &  vertical_profile_nc,      & !< name of the oae vertical profile
    &  hour_of_day_nc,           & !< name of the oae hour of day file
    &  day_of_week_nc,           & !< name of the oae day of week file
    &  month_of_year_nc,         & !< name of the oae month of year file
    &  hour_of_year_nc,          & !< name of the oae hour of year file
    &  gridded_emissions_nc        !< name of the oae gridded emission file

  REAL(wp) ::                    restart_init_time

  !
  NAMELIST /oemctrl_nml/ vertical_profile_nc,   &
    &                    hour_of_day_nc,        &
    &                    day_of_week_nc,        &
    &                    month_of_year_nc,      &
    &                    hour_of_year_nc,       &
    &                    gridded_emissions_nc

!==============================================================================
! Module procedure in "mo_oem_nml"
!==============================================================================

CONTAINS


  SUBROUTINE read_oemctrl_namelist( filename )

    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER :: istat, funit, iz_err
    INTEGER :: iunit

    iz_err = 0

    !0!CHARACTER(len=*), PARAMETER ::  &
    !0!  &  routine = 'mo_oem_nml:read_oemctrl_namelist'

    !-----------------------
    ! 1. default settings   
    !-----------------------

    vertical_profile_nc   = ''
    hour_of_day_nc        = ''
    day_of_week_nc        = ''
    month_of_year_nc      = ''
    hour_of_year_nc       = ''
    gridded_emissions_nc  = ''


    !--------------------------------------------------------------------
    ! 2. Read user's (new) specifications (Done so far by all MPI processes)
    !--------------------------------------------------------------------

    CALL open_nml(TRIM(filename))
    CALL position_nml ('oemctrl_nml', status=istat)
    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, oemctrl_nml)  ! write defaults to temporary text file
    END IF
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, oemctrl_nml)          ! overwrite default settings
      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, oemctrl_nml)  ! write settings to temporary text file
      END IF
    END SELECT
    CALL close_nml

    !----------------------------------------------------
    ! 3. Fill the configuration state
    !----------------------------------------------------

    config_vertical_profile_nc   = vertical_profile_nc
    config_hour_of_day_nc        = hour_of_day_nc
    config_day_of_week_nc        = day_of_week_nc
    config_month_of_year_nc      = month_of_year_nc
    config_hour_of_year_nc       = hour_of_year_nc
    config_gridded_emissions_nc  = gridded_emissions_nc

    !-----------------------------------------------------
    ! 4. write the contents of the namelist to an ASCII file
    !-----------------------------------------------------

    IF(my_process_is_stdio()) WRITE(nnml_output,nml=oemctrl_nml)


  END SUBROUTINE read_oemctrl_namelist


!------------------------------------------------------------------------------
! End of module mo_oem_nml
!------------------------------------------------------------------------------

END MODULE mo_oem_nml

