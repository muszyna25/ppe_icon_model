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
  USE mo_nml_annotate,        ONLY: temp_defaults, temp_settings

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
  INTEGER,PUBLIC :: i_ice_dyn           !< Dynamical model switch:
                                        ! 0: No dynamics
                                        ! 1: FEM dynamics (from AWI)
  INTEGER,PUBLIC :: i_Qio_type          !< Methods to calculate ice-ocean heatflux
                                        ! 1: Proportional to ocean cell thickness (like MPI-OM)
                                        ! 2: Proportional to speed difference btwn. ice and ocean

  REAL(wp),PUBLIC :: hnull              !< Hibler's h_0 for new ice formation
  REAL(wp),PUBLIC :: hmin               !< Minimum ice thickness allowed in the model
  REAL(wp),PUBLIC :: ramp_wind          !< Time (in days) that the wind stress is increased over.
                                        !  This is only necessary for runs which start off from a
                                        !  still ocean, i.e. not re-start
  REAL(wp),PUBLIC :: hci_layer          !< Thickness of stabilizing constant heat capacity layer
  REAL(wp),PUBLIC :: leadclose_1        !< Hibler's leadclose parameter for lateral melting

  INTEGER         :: iunit

  NAMELIST /sea_ice_nml/ kice, i_ice_therm, i_ice_albedo, i_ice_dyn, hnull, hmin, ramp_wind, &
    &           i_Qio_type, hci_layer, leadclose_1

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
    i_ice_albedo= 1
    i_ice_dyn   = 0
    i_Qio_type  = 2

    hnull       = 0.5_wp
    hmin        = 0.05_wp
    hci_layer   = 0.10_wp
    leadclose_1 = 0.5_wp

    ramp_wind   = 10._wp

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
    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, sea_ice_nml)    ! write defaults to temporary text file
    END IF
    SELECT CASE (istat)
    CASE (positioned)
      READ (nnml, sea_ice_nml)                                        ! overwrite default settings
      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, sea_ice_nml)    ! write settings to temporary text file
      END IF
    END SELECT
    CALL close_nml

    !------------------------------------------------------------------
    ! Sanity Check
    !------------------------------------------------------------------
    IF (kice /= 1 ) THEN
      CALL finish(TRIM(routine), 'Currently, kice must be 1.')
    END IF

    IF (i_ice_therm < 1 .OR. i_ice_therm > 4) THEN
      CALL finish(TRIM(routine), 'i_ice_therm must be between 1 and 4.')
    END IF

    IF (i_ice_albedo < 1 .OR. i_ice_albedo > 2 ) THEN
      CALL finish(TRIM(routine), 'i_ice_albedo must be either 1 or 2.')
    END IF

    IF (i_ice_dyn < 0 .OR. i_ice_dyn > 1) THEN
      CALL finish(TRIM(routine), 'i_ice_dyn must be either 0 or 1.')
    END IF

    ! TODO: This can be changed when we start advecting T1 and T2
    IF (i_ice_dyn == 1 ) THEN
      CALL message(TRIM(routine), 'i_ice_therm set to 1 because i_ice_dyn is 1')
      i_ice_therm = 1
    ENDIF


    IF (i_Qio_type < 1 .OR. i_Qio_type > 2) THEN
      CALL finish(TRIM(routine), 'i_Qio_type must be either 1 or 2.')
    END IF

    IF (i_ice_dyn == 0) THEN
      CALL message(TRIM(routine), 'i_Qio_type set to 1 because i_ice_dyn is 0')
      i_Qio_type = 1
    ENDIF

    IF (hmin > hnull) THEN
      CALL message(TRIM(routine), 'hmin cannot be larger than hnull')
      CALL message(TRIM(routine), 'setting hmin to hnull')
      hmin = hnull
    ENDIF

    IF (ramp_wind <= 0) THEN
      CALL message(TRIM(routine), 'ramp_wind cannot be smaller than or equal to zero')
      CALL message(TRIM(routine), 'setting ramp_wind to TINY(1._wp)')
      ramp_wind = TINY(1._wp)
    ENDIF

    IF (hci_layer < 0) THEN
      CALL message(TRIM(routine), 'hci_layer < 0, setting it equal to zero')
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
