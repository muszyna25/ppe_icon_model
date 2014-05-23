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
  USE mo_master_control,      ONLY: is_restart_run
  USE mo_io_restart_namelist, ONLY: open_tmpfile, store_and_close_namelist, &
                                  & open_and_restore_namelist, close_tmpfile
  USE mo_mpi,                 ONLY: my_process_is_stdio
  USE mo_nml_annotate,        ONLY: temp_defaults, temp_settings

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: read_echam_conv_namelist

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  !--------------------------------------------------------------
  ! Namelist variables 
  !--------------------------------------------------------------

  INTEGER :: iconv     !< 1,2,3 for different convection schemes
  INTEGER :: ncvmicro  !< 0 or 1. Scheme for convective microphysics

  LOGICAL :: lmfpen    !< true when penetrative convection is switched on
!  LOGICAL :: lmfmid    !< true when midlevel    convection is switched on
!  LOGICAL :: lmfdd     !< true when cumulus downdraft      is switched on
!  LOGICAL :: lmfdudv   !< true when cumulus friction       is switched on



  NAMELIST/echam_conv_nml/ ncvmicro, iconv,   &
                           lmfpen 


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
    !------------------------------------------------------------
    ncvmicro = 0
    iconv    = 1

    lmfpen   = .TRUE.

    !-------------------------------------------------------------------
    ! If this is a resumed integration, overwrite the defaults above
    ! by values used in the previous integration.
    !-------------------------------------------------------------------
    IF (is_restart_run()) THEN
      funit = open_and_restore_namelist('echam_conv_nml')
      READ(funit,NML=echam_conv_nml)
      CALL close_tmpfile(funit)
    END IF

    !-------------------------------------------------------------------------
    ! Read user's (new) specifications. (Done so far by all MPI processes)
    !-------------------------------------------------------------------------
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
    !------------------------------------------------------------
    SELECT CASE (iconv)
    CASE(1,2,3)  !OK
    CASE default
      WRITE(message_text,'(a,i0,a)') 'iconv = ',iconv,' is not supported'
      CALL finish(TRIM(routine),message_text)
    END SELECT

    SELECT CASE(ncvmicro)
    CASE (0)  !OK
    CASE DEFAULT
      CALL finish(TRIM(routine),'ncvmicro > 0 not yet supported in ICON')
    END SELECT

    !-----------------------------------------------------
    ! Store the namelist for restart
    !-----------------------------------------------------
    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=echam_conv_nml)
      CALL store_and_close_namelist(funit, 'echam_conv_nml')
    ENDIF
    !-----------------------------------------------------
    ! Fill configuration state
    !-----------------------------------------------------
    echam_conv_config% iconv    = iconv
    echam_conv_config% ncvmicro = ncvmicro
    echam_conv_config% lmfpen   = lmfpen

  END SUBROUTINE read_echam_conv_namelist
  !-------------

END MODULE mo_echam_conv_nml
