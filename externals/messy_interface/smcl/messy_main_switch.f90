! *************************************************************************
MODULE messy_main_switch
! *************************************************************************

  ! DEFINITION AND NAMELIST INPUT OF GLOBAL SUBMODEL SWITCHES
  ! NOTE: TO ADD A NEW SWITCH LOOK FOR
  !       '### ADD HERE'

  IMPLICIT NONE
  PUBLIC

  CHARACTER(LEN=*), PARAMETER :: modstr = 'switch'
  CHARACTER(LEN=*), PARAMETER :: modver = '1.0'

  ! GLOBAL SWITCHES
  !
  LOGICAL :: L_TIME_INFO  = .true.
  !
  LOGICAL :: USE_AEROPT   = .false.
  LOGICAL :: USE_AIRSEA   = .false.
  LOGICAL :: USE_BUFLY    = .false.
  LOGICAL :: USE_CLOUD    = .false.
  LOGICAL :: USE_CONVECT  = .false.
  LOGICAL :: USE_CVTRANS  = .false.
  LOGICAL :: USE_D14CO    = .false.
  LOGICAL :: USE_DDEP     = .false.
  LOGICAL :: USE_DRADON   = .false.
  LOGICAL :: USE_E4CHEM   = .false.
  LOGICAL :: USE_GMXE     = .false.
  LOGICAL :: USE_GWAVE    = .false.
  LOGICAL :: USE_H2O      = .false.
  LOGICAL :: USE_HETCHEM  = .false.
  LOGICAL :: USE_JVAL     = .false.
  LOGICAL :: USE_LNOX     = .false.
  LOGICAL :: USE_M7       = .false.
  LOGICAL :: USE_MADE     = .false.
  LOGICAL :: USE_MECCA1   = .false.
  LOGICAL :: USE_MECCA    = .false.
  LOGICAL :: USE_MEGAN    = .false.
  LOGICAL :: USE_MLOCEAN  = .false.
  LOGICAL :: USE_MMFORCE  = .false.
  LOGICAL :: USE_MSBM     = .false.
  LOGICAL :: USE_O3ORIG   = .false.
  LOGICAL :: USE_OFFEMIS  = .false.
  LOGICAL :: USE_ONEMIS   = .false.
  LOGICAL :: USE_ORBIT    = .false.
  LOGICAL :: USE_PHOTO    = .false.
  LOGICAL :: USE_PLUMEGAS = .false.
  LOGICAL :: USE_PSC      = .false.
  LOGICAL :: USE_PTRAC    = .false.
  LOGICAL :: USE_QBO      = .false.
  LOGICAL :: USE_RAD4ALL  = .false.
  LOGICAL :: USE_SATSIMS  = .false.
  LOGICAL :: USE_SCALC    = .false.
  LOGICAL :: USE_SCAV     = .false.
  LOGICAL :: USE_SCOUT    = .false.
  LOGICAL :: USE_SEDI     = .false.
  LOGICAL :: USE_S4D      = .false.
  LOGICAL :: USE_SORBIT   = .false.
  LOGICAL :: USE_SPACENOX = .false.
  LOGICAL :: USE_SPE      = .false.
  LOGICAL :: USE_TIMEPOS  = .false.
  LOGICAL :: USE_TNUDGE   = .false.
  LOGICAL :: USE_TREXP    = .false.
  LOGICAL :: USE_TROPOP   = .false.
  LOGICAL :: USE_VAHR     = .false.
  LOGICAL :: USE_VISO     = .false.
  LOGICAL :: USE_SUBMOD1  = .false.
  LOGICAL :: USE_SUBMOD2  = .false.
  LOGICAL :: USE_SUBMOD3  = .false.

  ! ### ADD HERE

  NAMELIST /CTRL/ L_TIME_INFO, &
     L_TIME_INFO,  USE_AEROPT,   USE_AIRSEA,   USE_BUFLY,      &
     USE_CLOUD,    USE_CONVECT,  USE_CVTRANS,  USE_D14CO,      &
     USE_DDEP,     USE_DRADON,   USE_E4CHEM,   USE_GMXE,       &
     USE_GWAVE,    USE_H2O,      USE_HETCHEM,  USE_JVAL,       &
     USE_LNOX,     USE_M7,       USE_MADE,     USE_MECCA,      &
     USE_MECCA1,   USE_MEGAN,    USE_MLOCEAN,  USE_MMFORCE,    &
     USE_MSBM,     USE_O3ORIG,   USE_OFFEMIS,  USE_ONEMIS,     &
     USE_ORBIT,    USE_PHOTO,    USE_PLUMEGAS, USE_PSC,        &
     USE_PTRAC,    USE_QBO,      USE_RAD4ALL,  USE_S4D,        &
     USE_SATSIMS,  USE_SCALC,    USE_SCAV,     USE_SCOUT,      &
     USE_SEDI,     USE_SORBIT,   USE_SPACENOX, USE_SPE,        &
     USE_SUBMOD1,  USE_SUBMOD2,  USE_SUBMOD3,  USE_TIMEPOS,    &
     USE_TNUDGE,   USE_TREXP,    USE_TROPOP,   USE_VAHR,       &
     USE_VISO
  ! ### ADD HERE

  !PUBLIC :: messy_main_read_nml_ctrl
  !PUBLIC :: switch_init
  ! op_pj_20110621+
  ! temprarily moved to messy_main_switch_bi.f90 to reduce dependencies
  ! (via USE_MLOCEAN) of BML to SMCL
!!$  PRIVATE :: put_submodel_att
  ! op_pj_20110621-

CONTAINS

  ! -------------------------------------------------------------
  SUBROUTINE messy_main_read_nml_ctrl(status, iou)

    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status ! error status
    INTEGER, INTENT(IN)  :: iou    ! logical I/O unit
    
    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'messy_main_read_nml_ctrl'
    LOGICAL                     :: lex          ! file exists ?
    INTEGER                     :: fstat        ! file status

    ! INITIALIZE
    status = 1  ! DEFAULT: ERROR

    CALL read_nml_open(lex, substr, iou, 'CTRL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CTRL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! NO ERROR

  END SUBROUTINE messy_main_read_nml_ctrl
  ! -------------------------------------------------------------

! *************************************************************************
END MODULE messy_main_switch
! *************************************************************************
