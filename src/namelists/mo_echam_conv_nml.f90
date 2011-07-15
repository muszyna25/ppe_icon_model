!>
!! Namelist for configuring cumulus convection parameterization in
!! the ECHAM physics package
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
MODULE mo_echam_conv_nml

  USE mo_echam_conv_config,   ONLY: echam_conv_config
  USE mo_kind,                ONLY: wp
  USE mo_impl_constants,      ONLY: SUCCESS
  USE mo_exception,           ONLY: print_value,message,message_text,finish
  USE mo_io_units,            ONLY: nnml
  USE mo_namelist,            ONLY: position_nml, POSITIONED, open_nml, close_nml
  USE mo_master_nml,          ONLY: lrestart
  USE mo_io_restart_namelist, ONLY: open_tmpfile, store_and_close_namelist, &
                                  & open_and_restore_namelist, close_tmpfile

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
  LOGICAL :: lmfmid    !< true when midlevel    convection is switched on
  LOGICAL :: lmfdd     !< true when cumulus downdraft      is switched on
  LOGICAL :: lmfdudv   !< true when cumulus friction       is switched on

  REAL(wp) :: dlev     !< "zdlev" in subroutine "cuasc". 
                       !< Critical thickness (unit: Pa) necessary for the 
                       !< onset of convective precipitation
  REAL(wp) :: cmftau   !< characteristic adjustment time scale
                       !< (replaces "ztau" in "cumastr")
  REAL(wp) :: cmfctop  !< fractional convective mass flux across the top of cloud 
  REAL(wp) :: cprcon   !< coefficient for determining conversion
                       !< from cloud water to rain
  REAL(wp) :: cminbuoy !< minimum excess buoyancy
  REAL(wp) :: entrpen  !< entrainment rate for penetrative convection

  NAMELIST/echam_conv_nml/ ncvmicro, iconv,   &
                           lmfpen,   lmfmid,  &
                           lmfdd,    lmfdudv, &
                           dlev,     cmftau,  &
                           cmfctop,  cprcon,  &
                           cminbuoy, entrpen

CONTAINS
  !>
  !! Read the convection namelist
  !!
  SUBROUTINE read_echam_conv_namelist( filename )

    CHARACTER(LEN=*),INTENT(IN) :: filename
    CHARACTER(LEN=*),PARAMETER :: &
    routine = 'mo_echam_conv_nml:read_echam_conv_namelist'

    INTEGER  :: ist, funit

    !------------------------------------------------------------
    ! Set default values
    !------------------------------------------------------------
    ncvmicro = 0
    iconv    = 1

    lmfpen   = .TRUE.
    lmfmid   = .TRUE.
    lmfdd    = .TRUE.
    lmfdudv  = .TRUE.

    dlev     = 3.0E4_wp   ! 300 hPa
    cmftau   = 10800._wp  ! 3 hours
    cmfctop  = 0.3_wp

    cprcon   = 1.E-4_wp
    cminbuoy = 0.025_wp
    entrpen  = 1.0E-4_wp

    !-------------------------------------------------------------------
    ! If this is a resumed integration, overwrite the defaults above
    ! by values used in the previous integration.
    !-------------------------------------------------------------------
    IF (lrestart) THEN
      funit = open_and_restore_namelist('echam_conv_nml')
      READ(funit,NML=echam_conv_nml)
      CALL close_tmpfile(funit)
    END IF

    !-------------------------------------------------------------------------
    ! Read user's (new) specifications. (Done so far by all MPI processors)
    !-------------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml('echam_conv_nml',STATUS=ist)
    SELECT CASE (ist)
    CASE (POSITIONED)
      READ (nnml, echam_conv_nml)
    END SELECT
    CALL close_nml

    !------------------------------------------------------------
    ! Sanity check
    !------------------------------------------------------------
    CALL message('','')
    CALL message('','------- namelist echam_conv_nml --------')

    SELECT CASE (iconv)
    CASE(1); CALL message('','--- iconv = 1 -> Convection: Nordeng (default)')
    CASE(2); CALL message('','--- iconv = 2 -> Convection: Tiedtke')
    CASE(3); CALL message('','--- iconv = 3 -> Convection: Hybrid')
    CASE default
      WRITE(message_text,'(a,i0,a)') 'iconv = ',iconv,' is not supported'
      CALL finish(TRIM(routine),message_text)
    END SELECT

    SELECT CASE(ncvmicro)
    CASE (0); CALL message('','--- ncvmicro = 0')
    CASE DEFAULT
      CALL finish(TRIM(routine),'ncvmicro > 0 not yet supported in ICON')
    END SELECT

    CALL print_value(' lmfpen  ',lmfpen)
    CALL print_value(' lmfmid  ',lmfmid)
    CALL print_value(' lmfdd   ',lmfdd)
    CALL print_value(' lmfdudv ',lmfdudv)

    CALL print_value(' cmftau   ',cmftau)
    CALL print_value(' cmfctop  ',cmfctop)
    CALL print_value(' cprcon   ',cprcon)
    CALL print_value(' cminbuoy ',cminbuoy)
    CALL print_value(' entrpen  ',entrpen)
    CALL print_value(' dlev     ',dlev)

    CALL message('','---------------------------')
    CALL message('','')

    !-----------------------------------------------------
    ! Store the namelist for restart
    !-----------------------------------------------------
    funit = open_tmpfile()
    WRITE(funit,NML=echam_conv_nml)
    CALL store_and_close_namelist(funit, 'echam_conv_nml')

    !-----------------------------------------------------
    ! Fill configuration state
    !-----------------------------------------------------
    echam_conv_config% iconv    = iconv
    echam_conv_config% ncvmicro = ncvmicro
    echam_conv_config% lmfpen   = lmfpen
    echam_conv_config% lmfmid   = lmfmid
    echam_conv_config% lmfdd    = lmfdd
    echam_conv_config% lmfdudv  = lmfdudv

    echam_conv_config% dlev     = dlev
    echam_conv_config% cmftau   = cmftau
    echam_conv_config% cmfctop  = cmfctop
    echam_conv_config% cprcon   = cprcon
    echam_conv_config% cminbuoy = cminbuoy
    echam_conv_config% entrpen  = entrpen

  END SUBROUTINE read_echam_conv_namelist
  !-------------

END MODULE mo_echam_conv_nml
