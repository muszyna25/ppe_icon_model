!>
!! Namelist for configuring cumulus convection parameterization in
!! the ECHAM physics package
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
MODULE mo_echam_conv_nml

  USE mo_echam_conv_config,   ONLY: echam_conv_config
  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message_text,finish
  USE mo_io_units,            ONLY: nnml
  USE mo_namelist,            ONLY: position_nml, POSITIONED, open_nml, close_nml
  USE mo_master_control,      ONLY: use_restart_namelists
  USE mo_restart_namelist,    ONLY: open_tmpfile, store_and_close_namelist, &
                                  & open_and_restore_namelist, close_tmpfile
  USE mo_mpi,                 ONLY: my_process_is_stdio
  USE mo_nml_annotate,        ONLY: temp_defaults, temp_settings

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: read_echam_conv_namelist

  !--------------------------------------------------------------
  ! Namelist variables 
  !
  LOGICAL  :: lmfmid    !< true when midlevel    convection is switched on
  LOGICAL  :: lmfpen    !< true when penetrative convection is switched on
  LOGICAL  :: lmfdd     !< true when cumulus downdraft      is switched on
  LOGICAL  :: lmfdudv   !< true when cumulus friction       is switched on
  !
  REAL(wp) :: entrscv   !< average entrainment rate for shallow convection
  REAL(wp) :: entrmid   !< average entrainment rate for midlevel convection
  REAL(wp) :: entrpen   !< average entrainment rate for penetrative convection
  REAL(wp) :: entrdd    !< average entrainment rate for cumulus downdrafts
  !
  REAL(wp) :: cprcon    !< coefficient for determining conversion from cloud water to rain
  REAL(wp) :: cmfctop   !< fractional convective mass flux across the top of cloud
  REAL(wp) :: cmfdeps   !< fractional convective mass flux for downdrafts at lfs
  !
  REAL(wp) :: cminbuoy  !< minimum excess buoyancy
  REAL(wp) :: cmaxbuoy  !< maximum excess buoyancy
  REAL(wp) :: cbfac     !< factor for std dev of virtual pot temp
  REAL(wp) :: centrmax  !< maximum entrainment/detrainment rate
  !
  REAL(wp) :: dlev_land !< minimum cloud pressure depth for precipitation over land
  REAL(wp) :: dlev_ocean!< minimum cloud pressure depth for precipitation over ocean
  !
  REAL(wp) :: cmftau    !< characteristic adjustment time scale (s)
  !
  NAMELIST /echam_conv_nml/                     &
       & lmfmid  , lmfpen  , lmfdd  , lmfdudv , &
       & entrscv , entrmid , entrpen, entrdd  , &
       & cprcon  , cmfctop , cmfdeps,           &
       & cminbuoy, cmaxbuoy, cbfac  , centrmax, &
       & dlev_land, dlev_ocean, cmftau

CONTAINS
  !>
  !! Read the convection namelist
  !!
  SUBROUTINE read_echam_conv_namelist( filename )

    CHARACTER(LEN=*),INTENT(IN) :: filename
    CHARACTER(LEN=*),PARAMETER :: &
    routine = 'mo_echam_conv_nml:read_echam_conv_namelist'

    INTEGER  :: ist, funit, iunit

    !------------------------------------------------------------
    ! Set default values
    !
    lmfmid   = .TRUE.
    lmfpen   = .TRUE.
    lmfdd    = .TRUE.
    lmfdudv  = .TRUE.
    !
    entrscv  = 3.0e-3_wp
    entrmid  = 1.0e-4_wp
    entrpen  = 1.0e-4_wp
    entrdd   = 2.0e-4_wp
    !
    cprcon   = 1.5e-4_wp
    cmfctop  = 0.23_wp
    cmfdeps  = 0.3_wp
    !
    cminbuoy = 0.2_wp
    cmaxbuoy = 1.0_wp
    cbfac    = 1.0_wp
    centrmax = 3.0e-4_wp
    !
    dlev_land  = 3.0e4_wp
    dlev_ocean = 1.5e4_wp
    !
    cmftau   = 3600.0_wp

    !-------------------------------------------------------------------
    ! If this is a resumed integration, overwrite the defaults above
    ! by values used in the previous integration.
    !
    IF (use_restart_namelists()) THEN
      funit = open_and_restore_namelist('echam_conv_nml')
      READ(funit,NML=echam_conv_nml)
      CALL close_tmpfile(funit)
    END IF

    !-------------------------------------------------------------------------
    ! Read user's (new) specifications. (Done so far by all MPI processes)
    !
    CALL open_nml(TRIM(filename))
    CALL position_nml('echam_conv_nml',STATUS=ist)
    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, echam_conv_nml)    ! write defaults to temporary text file
    END IF
    SELECT CASE (ist)
    CASE (POSITIONED)
      READ (nnml, echam_conv_nml)                                        ! overwrite default settings
      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, echam_conv_nml)    ! write settings to temporary text file
      END IF
    END SELECT
    CALL close_nml

    !------------------------------------------------------------
    ! Sanity check
    !

    !-----------------------------------------------------
    ! Store the namelist for restart
    !
    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=echam_conv_nml)
      CALL store_and_close_namelist(funit, 'echam_conv_nml')
    ENDIF

    !-----------------------------------------------------
    ! Fill configuration state
    !
    echam_conv_config% lmfmid   = lmfmid
    echam_conv_config% lmfpen   = lmfpen
    echam_conv_config% lmfdd    = lmfdd
    echam_conv_config% lmfdudv  = lmfdudv
    !
    echam_conv_config% entrscv  = entrscv
    echam_conv_config% entrmid  = entrmid
    echam_conv_config% entrpen  = entrpen
    echam_conv_config% entrdd   = entrdd
    !
    echam_conv_config% cprcon   = cprcon
    echam_conv_config% cmfctop  = cmfctop
    echam_conv_config% cmfdeps  = cmfdeps
    !
    echam_conv_config% cminbuoy = cminbuoy
    echam_conv_config% cmaxbuoy = cmaxbuoy
    echam_conv_config% cbfac    = cbfac
    echam_conv_config% centrmax = centrmax
    !
    echam_conv_config% dlev_land  = dlev_land
    echam_conv_config% dlev_ocean = dlev_ocean
    !
    echam_conv_config% cmftau   = cmftau

  END SUBROUTINE read_echam_conv_namelist
  !-------------

END MODULE mo_echam_conv_nml
