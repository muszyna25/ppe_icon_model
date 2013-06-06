!>
!! Namelist file for large-scale forcing terms
!!        
!! @par Revision History
!!
!! @par Copyright
!! 2002-2006 by DWD and MPI-M
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
!!
MODULE mo_ls_forcing_nml

  USE mo_kind,                ONLY: wp
  USE mo_mpi,                 ONLY: my_process_is_stdio 
  USE mo_exception,           ONLY: message, finish, message_text
  USE mo_io_units,            ONLY: nnml, nnml_output
  USE mo_namelist,            ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_master_control,      ONLY: is_restart_run
  USE mo_io_restart_namelist, ONLY: open_tmpfile, store_and_close_namelist,  &
                                  & open_and_restore_namelist, close_tmpfile
  USE mo_impl_constants,      ONLY: MAX_CHAR_LENGTH, max_dom
  USE mo_run_config,          ONLY: ltestcase
  USE mo_atm_phy_nwp_config,  ONLY: atm_phy_nwp_config

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: read_ls_forcing_namelist, is_ls_forcing, is_subsidence, is_advection, &
            is_geowind, is_rad_forcing, is_theta

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  LOGICAL  :: is_ls_forcing  !true if any forcing is on
  LOGICAL  :: is_subsidence  !true if subsidence is on
  LOGICAL  :: is_advection   !true if horizontal advective forcing is on for any variable
  LOGICAL  :: is_geowind     !true if geostophic wind is set 
  LOGICAL  :: is_rad_forcing !true if radiative forcing is on
  LOGICAL  :: is_theta       !true is forcings are in terms of theta
 
  NAMELIST/ls_forcing_nml/ is_subsidence, is_advection, is_geowind, is_rad_forcing,  &
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
    INTEGER :: istat, funit, jg

    CHARACTER(len=*), PARAMETER ::  &
      &  routine = 'mo_ls_forcing_nml: read_ls_forcing_namelist'

    !-----------------------
    ! 1. default settings
    !-----------------------
    is_ls_forcing = .FALSE.
    is_subsidence = .FALSE.
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
    IF (is_restart_run() .AND. ltestcase) THEN 
      funit = open_and_restore_namelist('ls_forcing_nml')
      READ(funit,NML=ls_forcing_nml)
      CALL close_tmpfile(funit)
    END IF

    !------------------------------------------------------------------------
    ! 3. Read user's (new) specifications (Done so far by all MPI processes)
    !------------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('ls_forcing_nml', status=istat)
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, ls_forcing_nml)
    END SELECT
    CALL close_nml

    !4. checks
    !If any of the forcing is ON turn on is_ls_forcing
    IF(is_subsidence .OR. is_advection .OR. is_geowind .OR. is_rad_forcing) &
        is_ls_forcing = .TRUE.

    IF(is_ls_forcing .AND. .NOT.ltestcase) &
        CALL message(TRIM(routine),'ltestcase is turned ON because is_ls_forcing is ON!')

    !Check for testcases with large-scale forcing
    IF(is_rad_forcing .AND. atm_phy_nwp_config(1)%inwp_radiation>0) &
        CALL finish(TRIM(routine),'both inwp_rad and rad_forcing are turned on!')

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
