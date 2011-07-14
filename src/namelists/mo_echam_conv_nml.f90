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
  PUBLIC  :: read_echam_conv_namelist
  PUBLIC  :: echam_conv_nml_setup 
  PUBLIC  :: cleanup_cuparam

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  !--------------------------------------------------------------
  ! Namelist variables 
  !--------------------------------------------------------------

  INTEGER :: nml_iconv     !< 1,2,3 for different convection schemes
  INTEGER :: nml_ncvmicro  !< 0 or 1. Scheme for convective microphysics

  LOGICAL :: nml_lmfpen    !< true when penetrative convection is switched on
  LOGICAL :: nml_lmfmid    !< true when midlevel    convection is switched on
  LOGICAL :: nml_lmfdd     !< true when cumulus downdraft      is switched on
  LOGICAL :: nml_lmfdudv   !< true when cumulus friction       is switched on

  REAL(wp) :: nml_dlev     !< "zdlev" in subroutine "cuasc". Critical thickness (unit: Pa)
                           !< necessary for the onset of convective precipitation
  REAL(wp) :: nml_cmftau   !< characteristic adjustment time scale
                           !< (replaces "ztau" in "cumastr")
  REAL(wp) :: nml_cmfctop  !< fractional convective mass flux across the top of cloud 
  REAL(wp) :: nml_cprcon   !< coefficient for determining conversion
                           !< from cloud water to rain
  REAL(wp) :: nml_cminbuoy !< minimum excess buoyancy
  REAL(wp) :: nml_entrpen  !< entrainment rate for penetrative convection

  NAMELIST/echam_conv_nml/ nml_ncvmicro, nml_iconv,   &
                           nml_lmfpen,   nml_lmfmid,  &
                           nml_lmfdd,    nml_lmfdudv, &
                           nml_dlev,     nml_cmftau,  &
                           nml_cmfctop,  nml_cprcon,  &
                           nml_cminbuoy, nml_entrpen

CONTAINS
  !>
  !! Read the convection namelist
  !!
  SUBROUTINE read_echam_conv_namelist( filename )

    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER  :: ist, funit
    CHARACTER(LEN=*),PARAMETER :: &
    routine = 'mo_echam_conv_nml:read_echam_conv_namelist'

    !------------------------------------------------------------
    ! Set default values
    !------------------------------------------------------------
    nml_ncvmicro = 0
    nml_iconv    = 1

    nml_lmfpen   = .TRUE.
    nml_lmfmid   = .TRUE.
    nml_lmfdd    = .TRUE.
    nml_lmfdudv  = .TRUE.

    nml_dlev     = 3.0E4_wp   ! 300 hPa
    nml_cmftau   = 10800._wp  ! 3 hours
    nml_cmfctop  = 0.3_wp

    nml_cprcon   = 1.E-4_wp
    nml_cminbuoy = 0.025_wp
    nml_entrpen  = 1.0E-4_wp

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

    SELECT CASE (nml_iconv)
    CASE(1); CALL message('','--- nml_iconv = 1 -> Convection: Nordeng (default)')
    CASE(2); CALL message('','--- nml_iconv = 2 -> Convection: Tiedtke')
    CASE(3); CALL message('','--- nml_iconv = 3 -> Convection: Hybrid')
    CASE default
      WRITE(message_text,'(a,i0,a)') 'nml_iconv = ',nml_iconv,' is not supported'
      CALL finish(TRIM(routine),message_text)
    END SELECT

    SELECT CASE(nml_ncvmicro)
    CASE (0); CALL message('','--- nml_ncvmicro = 0')
    CASE DEFAULT
      CALL finish(TRIM(routine),'nml_ncvmicro > 0 not yet supported in ICON')
    END SELECT

    CALL print_value(' nml_lmfpen  ',nml_lmfpen)
    CALL print_value(' nml_lmfmid  ',nml_lmfmid)
    CALL print_value(' nml_lmfdd   ',nml_lmfdd)
    CALL print_value(' nml_lmfdudv ',nml_lmfdudv)

    CALL print_value(' nml_cmftau   ',nml_cmftau)
    CALL print_value(' nml_cmfctop  ',nml_cmfctop)
    CALL print_value(' nml_cprcon   ',nml_cprcon)
    CALL print_value(' nml_cminbuoy ',nml_cminbuoy)
    CALL print_value(' nml_entrpen  ',nml_entrpen)
    CALL print_value(' nml_dlev     ',nml_dlev)

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
    echam_conv_config% iconv    = nml_iconv
    echam_conv_config% ncvmicro = nml_ncvmicro
    echam_conv_config% lmfpen   = nml_lmfpen
    echam_conv_config% lmfmid   = nml_lmfmid
    echam_conv_config% lmfdd    = nml_lmfdd
    echam_conv_config% lmfdudv  = nml_lmfdudv

    echam_conv_config% dlev     = nml_dlev
    echam_conv_config% cmftau   = nml_cmftau
    echam_conv_config% cmfctop  = nml_cmfctop
    echam_conv_config% cprcon   = nml_cprcon
    echam_conv_config% cminbuoy = nml_cminbuoy
    echam_conv_config% entrpen  = nml_entrpen

  END SUBROUTINE read_echam_conv_namelist
  !>
  !!
  SUBROUTINE echam_conv_nml_setup

  ! REAL(wp) :: za, zb
  ! REAL(wp) :: zp(nlev), zph(nlevp1)
  ! INTEGER  :: jk, ist


  !  !------------------------------------------------------------
  !  ! CJ: calculations from the former cuparam subroutine
  !  !------------------------------------------------------------

  !  ! Determine highest level *nmctop* for cloud base of midlevel convection
  !  ! assuming nmctop=9 (300 hPa) for the standard 19 level model

  !  ! half level pressure values, assuming 101320. Pa surface pressure

  !  DO jk=1,nlevp1
  !    za = vct(jk)
  !    zb = vct(jk+nvclev)
  !    zph(jk) = za + zb*101320.0_wp
  !  END DO

  !  ! full level pressure

  !  DO jk = 1, nlev
  !    zp(jk) = (zph(jk)+zph(jk+1))*0.5_wp
  !  END DO

  !  ! search for 300 hPa level

  !  DO jk = 1, nlev
  !    nmctop=jk
  !    IF(zp(jk).GE.30000.0_wp) EXIT
  !  END DO

  !  CALL print_value(&
  !    &'lowest model level for cloud base of mid level convection: nmctop = ',nmctop)

  !  ! evaporation coefficient for kuo0

  !  ALLOCATE( cevapcu(nlev),STAT=ist )
  !  IF (ist/=SUCCESS) CALL finish('cuparam','allocation of cevapcu failed')

  !  DO jk = 1,nlev
  !     cevapcu(jk) = 1.E3_wp/(38.3_wp*0.293_wp)*SQRT(ceta(jk))
  !     cevapcu(jk) = 1.93E-6_wp*261._wp*SQRT(cevapcu(jk))*0.5_wp/grav
  !  END DO

  END SUBROUTINE echam_conv_nml_setup
  !------------
  !>
  !!
  SUBROUTINE cleanup_cuparam

  ! INTEGER :: ist

  ! DEALLOCATE( cevapcu,STAT=ist )
  ! IF (ist/=SUCCESS) CALL finish('cuparam_cleanup','deallocation of cevapcu failed')

  END SUBROUTINE cleanup_cuparam
  !-------------
END MODULE mo_echam_conv_nml
