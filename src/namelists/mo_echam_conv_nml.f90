!>
!! Namelist for the configuration of the convection parameters
!!
!!
!! @par Revision History
!! Revision history in mo_echam_conv_params.f90 (r4370)
!! Modification by Constantin Junk, MPI-M (2011-05-11)
!! - moved namelist parameters, setup_convection and cuparam
!!   to new module mo_echam_conv_nml
!! - included subroutine cuparam in setup_convection
!! - renamed setup_convection echam_conv_nml_setup
!! - moved echam_vdiff namelist variables and subroutine setup_vdiff
!!   from mo_echam_vdiff_params to namelists/mo_echam_vdiff_nml
!!
!! @par Copyright
!! 2002-2010 by DWD and MPI-M
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

  USE mo_kind,                ONLY: wp
  USE mo_impl_constants,      ONLY: SUCCESS
  USE mo_namelist,            ONLY: position_nml, POSITIONED
  USE mo_echam_conv_config,   ONLY: echam_conv_config
  USE mo_io_units,            ONLY: nnml
  USE mo_master_nml,          ONLY: lrestart
  USE mo_io_restart_namelist, ONLY: open_tmpfile, store_and_close_namelist, &
                                  & open_and_restore_namelist, close_tmpfile
 !USE mo_exception,           ONLY: print_value,message,message_text,finish
 !USE mo_run_nml,             ONLY: nlev,nlevp1,nvclev
 !USE mo_physical_constants,  ONLY: grav
 !USE mo_vertical_coord_table,ONLY: vct,ceta

  IMPLICIT NONE
  PRIVATE

  PUBLIC  :: iconv,lconvmassfix                     !< parameters
  PUBLIC  :: lmfpen,lmfmid,lmfscv,lmfdd,lmfdudv     !< parameters
  PUBLIC  :: cmftau,cmfdeps,cmfcmin,cmfcmax,cmfctop !< parameters
  PUBLIC  :: centrmax,cbfac,cminbuoy,cmaxbuoy       !< parameters
  PUBLIC  :: entrpen,entrmid,entrscv,entrdd         !< parameters
  PUBLIC  :: cprcon,cevapcu                         !< parameters
  PUBLIC  :: nmctop, dlev                           !< parameters
 !PUBLIC  :: echam_conv_ctl                         !< namelist
  PUBLIC  :: read_echam_conv_namelist
  PUBLIC  :: echam_conv_nml_setup                   !< subroutine
  PUBLIC  :: cleanup_cuparam                        !< subroutine 

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  !--------------------------------------------------------------!
  ! echam_conv_nml namelist variables and auxiliary parameters
  !--------------------------------------------------------------!

  INTEGER :: nml_ncvmicro     !< 0 or 1. Scheme for convective microphysics
  INTEGER :: iconv        !< 1,2,3 for different convection schemes
  LOGICAL :: lconvmassfix !< aerosol mass fixer in convection

  LOGICAL :: lmfpen    !< true when penetrative convection is switched on
  LOGICAL :: lmfmid    !< true when midlevel    convection is switched on
  LOGICAL :: lmfscv    !< true when shallow     convection is switched on
  LOGICAL :: lmfdd     !< true when cumulus downdraft      is switched on
  LOGICAL :: lmfdudv   !< true when cumulus friction       is switched on

  REAL(wp) :: dlev     !< "zdlev" in subroutine "cuasc". Critical thickness (unit: Pa)
                       !< necessary for the onset of convective precipitation

  REAL(wp) :: cmftau   !< characteristic adjustment time scale
                       !< (replaces "ztau" in "cumastr"
  REAL(wp) :: cmfdeps  !< fractional convective mass flux for downdrafts at lfs
  REAL(wp) :: cmfctop  !< fractional convective mass flux across the top of cloud 
  REAL(wp) :: cmfcmin  !< minimum massflux value (for safety)
  REAL(wp) :: cmfcmax  !< maximum massflux value allowed for
  REAL(wp) :: centrmax !<

  REAL(wp) :: cmaxbuoy !< maximum excess buoyancy
  REAL(wp) :: cminbuoy !< minimum excess buoyancy

  REAL(wp) :: cbfac    !< factor for std dev of virtual pot temp

  REAL(wp) :: entrpen  !< entrainment rate for penetrative convection
  REAL(wp) :: entrmid  !< entrainment rate for midlevel convection
  REAL(wp) :: entrscv  !< entrainment rate for shallow convection
  REAL(wp) :: entrdd   !< entrainment rate for cumulus downdrafts

  REAL(wp) :: cprcon   !< coefficient for determining conversion
                       !< from cloud water to rain
  INTEGER :: nmctop    !< max. level for cloud base of mid level conv.

  REAL(wp),ALLOCATABLE :: cevapcu(:)  !< evaporation coefficient for kuo0
                                      !< In ECHAM6 it is defined in mo_physc2,
                                      !< allocated in subroutine alloc_mods,
                                      !< and initialized in subroutine iniphy.

 !INTEGER :: nml_nauto        !< 1 or 2. autoconversion scheme

  NAMELIST/echam_conv_ctl/ nml_ncvmicro,iconv,lconvmassfix,       &
    &                      lmfpen,lmfmid,lmfscv,lmfdd,lmfdudv,      &
    &                      cmftau,cmfctop,cprcon,cmfdeps,cminbuoy,  &
    &                      entrpen, entrmid,entrscv,entrdd,dlev
   !&                      nml_nauto

CONTAINS
  !>
  !! Read the convection namelist
  !!
  SUBROUTINE read_echam_conv_namelist()

    INTEGER  :: ist, funit

    !------------------------------------------------------------
    ! set up the default values for echam_conv_ctl
    !------------------------------------------------------------

    nml_ncvmicro     = 0
    iconv        = 1
    lmfpen       = .TRUE.
    lmfmid       = .TRUE.
    lmfscv       = .TRUE.
    lmfdd        = .TRUE.
    lmfdudv      = .TRUE.
    lconvmassfix = .FALSE.
    cmftau       = 10800._wp  ! 3 hours
    cmfctop      = 0.3_wp
    cmfdeps      = 0.3_wp     ! Fractional massflux for downdrafts at lfs
    cprcon       = 1.E-4_wp
    cminbuoy     = 0.025_wp
    entrpen      = 1.0E-4_wp  ! average entrainment rate for penetrative convection
    entrmid      = 1.0E-4_wp  ! average entrainment rate for midlevel convection
    entrscv      = 3.0E-4_wp  ! average entrainment rate for shallow convection
    entrdd       = 2.0E-4_wp  ! average entrainment rate for downdrafts
    dlev         = 3.0E4_wp   ! 300 hPa

   !nml_nauto        = 1

    ! Set default values of auxiliary parameters

    cmfcmin    = 1.E-10_wp  ! Minimum massflux value (for safety)
    cmfcmax    = 1.0_wp     ! Maximum massflux value allowed for updrafts etc
    cmaxbuoy   = 1.0_wp
    cbfac      = 1.0_wp
    centrmax   = 3.E-4_wp

    !----------------------------------------------------------------
    ! If this is a resumed integration, overwrite the defaults above
    ! by values in the previous integration.
    !----------------------------------------------------------------
    IF (lrestart) THEN
      funit = open_and_restore_namelist('echam_conv_ctl')
      READ(funit,NML=echam_conv_ctl)
      CALL close_tmpfile(funit)
    END IF

    !-------------------------------------------------------------------------
    ! 3. Read user's (new) specifications. (Done so far by all MPI processes)
    !-------------------------------------------------------------------------
    CALL position_nml('echam_conv_ctl',STATUS=ist)
    SELECT CASE (ist)
    CASE (POSITIONED)
      READ (nnml, echam_conv_ctl)
    END SELECT

    !-----------------------------------------------------
    ! Store the namelist for restart
    !-----------------------------------------------------
    funit = open_tmpfile()
    WRITE(funit,NML=echam_conv_ctl)
    CALL store_and_close_namelist(funit, 'echam_conv_ctl')

    !-----------------------------------------------------
    ! Fill configuration state
    !-----------------------------------------------------
    echam_conv_config% ncvmicro = nml_ncvmicro
   !echam_conv_config% nauto    = nml_nauto

  END SUBROUTINE read_echam_conv_namelist

  !>
  !!
  SUBROUTINE echam_conv_nml_setup

  ! REAL(wp) :: za, zb
  ! REAL(wp) :: zp(nlev), zph(nlevp1)
  ! INTEGER  :: jk, ist

  !  !------------------------------------------------------------
  !  ! check the consistency of the parameters
  !  !------------------------------------------------------------
  !  CALL message('','')
  !  CALL message('','------- namelist echam_conv_ctl --------')

  !  SELECT CASE(ncvmicro)
  !  CASE (0)
  !    CALL message('','--- ncvmicro = 0')
  !  CASE DEFAULT
  !    CALL finish('echam_conv_nml_setup','ncvmicro > 0 not yet supported in ICON')
  !  END SELECT

  !  SELECT CASE(nauto)
  !  CASE (0)
  !  CASE (1)
  !    CALL message('','--- nauto = 1 --> Beheng (1994) - ECHAM5 Standard')
  !  CASE (2)
  !    CALL message('','--- nauto = 2 --> Khairoutdinov and Kogan (2000)')
  !  CASE DEFAULT
  !    WRITE(message_text,'(a,i0,a)') 'nauto = ',nauto,' is not supported'
  !    CALL message('NAMELIST echam_conv_ctl',message_text)
  !    CALL finish('echam_conv_nml_setup','Run terminated')
  !  END SELECT

  !  SELECT CASE (iconv)
  !  CASE(1)
  !    CALL message('','--- iconv = 1 --> Convection: Nordeng (default)')
  !  CASE(2)
  !    CALL message('','--- iconv = 2 --> Convection: Tiedtke')
  !  CASE(3)
  !    CALL message('','--- iconv = 3 --> Convection: Hybrid')
  !  CASE default
  !    WRITE(message_text,'(a,i0,a)') 'iconv = ',iconv,' is not supported'
  !    CALL message('NAMELIST echam_conv_ctl',message_text)
  !    CALL finish('echam_conv_nml_setup','Run terminated')
  !  END SELECT

  !  CALL print_value(' lmfpen  ',lmfpen)
  !  CALL print_value(' lmfmid  ',lmfmid)
  !  CALL print_value(' lmfscv  ',lmfscv)
  !  CALL print_value(' lmfdd   ',lmfdd)
  !  CALL print_value(' lmfdudv ',lmfdudv)
  !  CALL print_value(' lconvmassfix = ',lconvmassfix)

  !  !cmftau = MIN(10800._wp,cmftau)
  !  CALL print_value(' cmftau   ',cmftau)
  !  CALL print_value(' cmfctop  ',cmfctop)
  !  CALL print_value(' cprcon   ',cprcon)
  !  CALL print_value(' cminbuoy ',cminbuoy)
  !  CALL print_value(' entrpen  ',entrpen)
  !  CALL print_value(' dlev     ',dlev)

  !  CALL message('','---------------------------')
  !  CALL message('','')


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

    INTEGER :: ist

    DEALLOCATE( cevapcu,STAT=ist )
    IF (ist/=SUCCESS) CALL finish('cuparam_cleanup','deallocation of cevapcu failed')

  END SUBROUTINE cleanup_cuparam
  !-------------
END MODULE mo_echam_conv_nml
