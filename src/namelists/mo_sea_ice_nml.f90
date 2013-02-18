!>
!! @brief namelist setup for the sea-ice model
!!
!! Namelist setup for the sea-ice model
!! <Details of procedures are documented below with their definitions.>
!! <Include any applicable external references inline as module::procedure,>
!! <external_procedure(), or by using @see.>
!! <Don't forget references to literature.>
!!
!! @author Einar Olason, MPI-M (2013-01-29)
!!
!!
!! @par Revision History
!! New file based on mo_lnd_jsbach_nml.f90 by Einar Olason, MPI-M (2013-01-29)
!!
!! @par Copyright
!! 2002-2012 by DWD and MPI-M
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
MODULE mo_sea_ice_nml

  USE mo_kind,                ONLY: wp
  USE mo_namelist,            ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_io_units,            ONLY: nnml, nnml_output
  USE mo_master_control,      ONLY: is_restart_run
  USE mo_io_restart_namelist, ONLY: open_tmpfile, store_and_close_namelist, &
                                  & open_and_restore_namelist, close_tmpfile
  USE mo_exception,           ONLY: finish, message
  USE mo_mpi,                 ONLY: my_process_is_stdio

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: read_sea_ice_namelist

  CHARACTER(len=*), PARAMETER, PRIVATE :: version = '$Id$'

  INTEGER,PUBLIC :: kice                !< Number of ice classes
  INTEGER,PUBLIC :: i_ice_therm         !< Thermodynamic model switch:
                                        ! 1: Zero layers
                                        ! 2: Winton's two layer model
                                        ! 3: Zero layers with analytical fluxes
                                        ! 4: Zero layers - only calculates surface temperature
  INTEGER,PUBLIC :: i_ice_albedo        !< Albedo model (no albedo model implemented yet)

  REAL(wp),PUBLIC :: hnull              !< Hibler's h_0 for new ice formation
  REAL(wp),PUBLIC :: hmin               !< Minimum ice thickness allowed in the model

  NAMELIST /sea_ice_nml/ kice, i_ice_therm, i_ice_albedo, hnull

CONTAINS
  !>
  !!
  SUBROUTINE read_sea_ice_namelist(filename)

    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER :: istat, funit

    CHARACTER(len=*), PARAMETER ::  &
      &  routine = 'mo_sea_ice_nml:read_sea_ice_namelist'

    !------------------------------------------------------------------
    ! Set default values
    !------------------------------------------------------------------
    kice        = 1
    i_ice_therm = 2
    i_ice_albedo= 999

    hnull       = 0.5_wp
    hmin        = hnull

    !------------------------------------------------------------------
    ! If this is a resumed integration, overwrite the defaults above
    ! by values used in the previous integration.
    !------------------------------------------------------------------
    IF (is_restart_run()) THEN
      funit = open_and_restore_namelist('sea_ice_nml')
      READ(funit,NML=sea_ice_nml)
      CALL close_tmpfile(funit)
    END IF

    !------------------------------------------------------------------
    ! Read user's (new) specifications (done by all MPI processes)
    !------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('sea_ice_nml', STATUS=istat)
    SELECT CASE (istat)
    CASE (positioned)
      READ (nnml, sea_ice_nml)
    END SELECT
    CALL close_nml

    !------------------------------------------------------------------
    ! Sanity Check
    !------------------------------------------------------------------
    IF (kice /= 1 ) THEN
      CALL finish(TRIM(routine), 'Currently, kice must be 1.')
    END IF

    IF (i_ice_therm < 1 .AND. i_ice_therm > 4) THEN
      CALL finish(TRIM(routine), 'ice therm must be between 1 and 4.')
    END IF

    IF (i_ice_albedo /= 999) THEN
      CALL message(TRIM(routine), 'only one albedo scheme implemented')
    END IF

    IF (hmin > hnull) THEN
      CALL message(TRIM(routine), 'hmin cannot be larger than hnull')
      CALL message(TRIM(routine), 'setting hmin to hnull')
      hmin = hnull
    ENDIF


    !------------------------------------------------------------------
    ! Store the namelist for restart
    !------------------------------------------------------------------
    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=sea_ice_nml)
      CALL store_and_close_namelist(funit, 'sea_ice_nml')
    ENDIF

    !------------------------------------------------------------------
    ! Write the namelist to an ASCII file
    !------------------------------------------------------------------
    IF ( my_process_is_stdio() ) WRITE(nnml_output,nml=sea_ice_nml)

    !------------------------------------------------------------------
    ! Fill the configuration state
    !------------------------------------------------------------------

  END SUBROUTINE read_sea_ice_namelist

END MODULE mo_sea_ice_nml
