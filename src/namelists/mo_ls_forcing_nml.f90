!>
!! Namelist file for large-scale forcing terms
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
!!
MODULE mo_ls_forcing_nml

  USE mo_kind,                ONLY: wp
  USE mo_mpi,                 ONLY: my_process_is_stdio 
  USE mo_exception,           ONLY: message, finish, message_text
  USE mo_io_units,            ONLY: nnml, nnml_output
  USE mo_namelist,            ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_master_control,      ONLY: use_restart_namelists
  USE mo_io_restart_namelist, ONLY: open_tmpfile, store_and_close_namelist,  &
                                  & open_and_restore_namelist, close_tmpfile
  USE mo_impl_constants,      ONLY: MAX_CHAR_LENGTH, max_dom
  USE mo_run_config,          ONLY: ltestcase
  USE mo_atm_phy_nwp_config,  ONLY: atm_phy_nwp_config
  USE mo_echam_phy_config,    ONLY: echam_phy_config
  USE mo_grid_config,         ONLY: is_plane_torus
  USE mo_nml_annotate,        ONLY: temp_defaults, temp_settings

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: read_ls_forcing_namelist, is_ls_forcing, is_subsidence_moment, is_subsidence_heat, &
            is_advection, is_geowind, is_rad_forcing, is_theta

  LOGICAL  :: is_ls_forcing  !true if any forcing is on
  LOGICAL  :: is_subsidence_moment  !true if subsidence is on for u and v
  LOGICAL  :: is_subsidence_heat    !true if subsidence is on for thermodyn. variables
  LOGICAL  :: is_advection   !true if horizontal advective forcing is on for any variable
  LOGICAL  :: is_geowind     !true if geostophic wind is set 
  LOGICAL  :: is_rad_forcing !true if radiative forcing is on
  LOGICAL  :: is_theta       !true is forcings are in terms of theta
 
  NAMELIST/ls_forcing_nml/ is_subsidence_moment, is_subsidence_heat, is_advection, is_geowind, is_rad_forcing,  &
                           is_theta

CONTAINS
  !-------------------------------------------------------------------------
  !>
  !! Read Namelist for LS forcing
  !!
  !! This subroutine 
  !! - reads the Namelist 
  !! - sets default values
  !! - potentially overwrites the defaults by values used in a 
  !!   previous integration (if this is a resumed run)
  !! - reads the user's (new) specifications
  !! - stores the Namelist for restart
  !! - fills the configuration state 
  !!
  !! @par Revision History
  !!  by Anurag Dipankar, MPIM (2013-May-31)
  !!
  SUBROUTINE read_ls_forcing_namelist( filename )

    CHARACTER(LEN=*), INTENT(IN) :: filename 
    INTEGER :: istat, funit, iunit

    CHARACTER(len=*), PARAMETER ::  &
      &  routine = 'mo_ls_forcing_nml: read_ls_forcing_namelist'

    !-----------------------
    ! 1. default settings
    !-----------------------
    is_ls_forcing = .FALSE.
    is_subsidence_moment = .FALSE.
    is_subsidence_heat   = .FALSE.
    is_advection  = .FALSE.
    is_geowind    = .FALSE.
    is_rad_forcing = .FALSE.
    is_theta       = .FALSE.

    !------------------------------------------------------------------
    ! 2. If this is a resumed integration, overwrite the defaults above 
    !    by values used in the previous integration.
    ! ltestcase is added here because it was causing trouble for AMIP runs
    ! restarting with a file generated long ago
    !------------------------------------------------------------------
    IF (use_restart_namelists() .AND. .NOT.echam_phy_config%lamip) THEN 
      funit = open_and_restore_namelist('ls_forcing_nml')
      READ(funit,NML=ls_forcing_nml)
      CALL close_tmpfile(funit)
    END IF

    !------------------------------------------------------------------------
    ! 3. Read user's (new) specifications (Done so far by all MPI processes)
    !------------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('ls_forcing_nml', status=istat)
    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, ls_forcing_nml)   ! write defaults to temporary text file
    END IF
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, ls_forcing_nml)                                       ! overwrite default settings
      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, ls_forcing_nml)   ! write settings to temporary text file
      END IF
    END SELECT
    CALL close_nml

    !4. checks
    !If any of the forcing is ON turn on is_ls_forcing
    IF(is_subsidence_moment .OR. is_subsidence_heat .OR. is_advection .OR. is_geowind .OR. is_rad_forcing) &
        is_ls_forcing = .TRUE.

    IF(is_ls_forcing .AND. .NOT.ltestcase) &
        CALL message(TRIM(routine),'ltestcase is turned ON because is_ls_forcing is ON!')

    !Check for testcases with large-scale forcing
    IF(is_rad_forcing .AND. atm_phy_nwp_config(1)%inwp_radiation>0) &
        CALL finish(TRIM(routine),'both inwp_rad and rad_forcing are turned on!')

    IF(is_geowind .AND. .NOT.is_plane_torus) & 
         CALL finish(TRIM(routine),'is_geowind is only applicable if is_plane_torus is turned on!')  

    !-----------------------------------------------------
    ! 5. Store the namelist for restart
    !-----------------------------------------------------
    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=ls_forcing_nml)                    
      CALL store_and_close_namelist(funit,'ls_forcing_nml') 
    ENDIF

    ! 6. write the contents of the namelist to an ASCII file
    IF(my_process_is_stdio()) WRITE(nnml_output,nml=ls_forcing_nml)

  END SUBROUTINE read_ls_forcing_namelist

END MODULE mo_ls_forcing_nml
