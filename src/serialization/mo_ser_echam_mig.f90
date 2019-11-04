!--------------------------------------------------------------------
!
! Serialization routine for ECHAM mig
!
!--------------------------------------------------------------------

MODULE mo_ser_echam_mig

  USE mo_kind,               ONLY: vp, wp
  USE mo_exception,          ONLY: warning
  USE mo_ser_common,         ONLY: init
  USE mtime,                 ONLY: datetimeToString, MAX_DATETIME_STR_LEN
  USE mo_time_config,        ONLY: time_config
  USE mo_echam_phy_memory,   ONLY: t_echam_phy_field, t_echam_phy_tend
  USE mo_run_config,         ONLY: iqv, iqc, iqi, iqr, iqs, iqg
  IMPLICIT NONE

  INTEGER, PARAMETER :: selected_block = 1

  PUBLIC :: serialize_mig_input
  PUBLIC :: serialize_mig_output

  CONTAINS

  SUBROUTINE serialize_mig_input(jg, jb, jcs, jce, nproma, nlev, field)
    INTEGER, INTENT(IN)    :: jg, jb, jcs, jce, nproma, nlev
    TYPE(t_echam_phy_field), POINTER, INTENT(INOUT) :: field
#if defined(SERIALIZE_ECHAM_MIG) || defined(SERIALIZE_ECHAM_ALL) || defined(SERIALIZE_ALL)

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: date
    LOGICAL, PARAMETER :: lonlyonce = .TRUE.
#if defined( SERIALIZE_CREATE_REFERENCE )
    LOGICAL, SAVE :: lenabled = .TRUE.
#else
    LOGICAL, SAVE :: lenabled = .FALSE.
#endif
    LOGICAL, SAVE :: lactive = .FALSE.

    !$ser verbatim IF (selected_block < 0 .OR. jb == selected_block) THEN
    !$ser verbatim   lactive = .TRUE.
    !$ser verbatim ENDIF
    !$ser verbatim IF (lenabled .and. lactive) THEN
    !$ser verbatim   CALL warning('SER:mo_ser_echam_mig:input','Serialization is active!')
#if defined( _OPENACC )
    !$ser verbatim   CALL warning('GPU:mo_ser_echam_mig:input','GPU HOST synchronization forced by serialization!')
    !$ser verbatim   !$ACC UPDATE HOST( field%ta ) IF( ASSOCIATED(field%ta) )
    !$ser verbatim   !$ACC UPDATE HOST( field%presm_old ) IF( ASSOCIATED(field%presm_old) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rho ) IF( ASSOCIATED(field%rho) )
    !$ser verbatim   !$ACC UPDATE HOST( field%qtrc ) IF( ASSOCIATED(field%qtrc) )
#endif
    !$ser verbatim   CALL datetimeToString(time_config%tc_current_date, date)
    !$ser verbatim   CALL init('icon')
    !$ser savepoint echam_mig-input jg=jg jb=jb jcs=jcs jce=jce nproma=nproma nlev=nlev date=TRIM(date)
    #if defined( SERIALIZE_CREATE_REFERENCE )
    !$ser mode write
#elif defined( SERIALIZE_PERTURB_REFERENCE )
    !$ser mode read-perturb
#elif defined( SERIALIZE_READ_REFERENCE )
    !$ser mode read
#else
#error SERIALIZATION MODE IS NOT SET
#endif
    !$ser data echam_mig_ta=field%ta(:,:,jb) IF (ASSOCIATED(field%ta))
    !$ser data echam_mig_presm_old=field%presm_old(:,:,jb) IF (ASSOCIATED(field%presm_old))
    !$ser data echam_mig_rho=field%rho(:,:,jb) IF (ASSOCIATED(field%rho))
    !$ser data echam_mig_qv=field%qtrc(:,:,jb,iqv) IF (ASSOCIATED(field%qtrc))
    !$ser data echam_mig_qc=field%qtrc(:,:,jb,iqc) IF (ASSOCIATED(field%qtrc))
    !$ser data echam_mig_qi=field%qtrc(:,:,jb,iqi) IF (ASSOCIATED(field%qtrc))
    !$ser data echam_mig_qr=field%qtrc(:,:,jb,iqr) IF (ASSOCIATED(field%qtrc))
    !$ser data echam_mig_qs=field%qtrc(:,:,jb,iqs) IF (ASSOCIATED(field%qtrc))
    !$ser data echam_mig_qg=field%qtrc(:,:,jb,iqg) IF (ASSOCIATED(field%qtrc))
    !$ser verbatim lactive = .FALSE.
    !$ser verbatim IF (lonlyonce) THEN
    !$ser verbatim   lenabled = .FALSE.
    !$ser verbatim END IF
#if defined( SERIALIZE_READ_REFERENCE ) || defined( SERIALIZE_PERTURB_REFERENCE )
#if defined( _OPENACC )
    !$ser verbatim CALL warning('GPU:mo_ser_echam_mig:input','GPU DEVICE synchronization forced by serialization!')
    !$ser verbatim   !$ACC UPDATE HOST( field%ta ) IF( ASSOCIATED(field%ta) )
    !$ser verbatim   !$ACC UPDATE HOST( field%presm_old ) IF( ASSOCIATED(field%presm_old) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rho ) IF( ASSOCIATED(field%rho) )
    !$ser verbatim   !$ACC UPDATE HOST( field%qtrc ) IF( ASSOCIATED(field%qtrc) )
#endif
#endif
    !$ser verbatim ENDIF

#endif
  END SUBROUTINE serialize_mig_input

  SUBROUTINE serialize_mig_after_satad1(jg, jb, jcs, jce, nproma, nlev, field, &
      &           xlta, xlqv, xlqc, xlqi, xlqr, xlqs, xlqg, zqnc, &
      &           zprec_r, zprec_s, zprec_g, zqrsflux)
    INTEGER, INTENT(IN)    :: jg, jb, jcs, jce, nproma, nlev
    TYPE(t_echam_phy_field), POINTER, INTENT(INOUT) :: field
    REAL(wp), DIMENSION(nproma, nlev), INTENT(INOUT) :: xlta, xlqv, xlqc, xlqi, &
      &                                                 xlqr, xlqs, xlqg, zqnc
    REAL(wp), DIMENSION(nproma, nlev), INTENT(INOUT) :: zprec_r, zprec_s, zprec_g, zqrsflux
#if defined(SERIALIZE_ECHAM_MIG) || defined(SERIALIZE_ECHAM_ALL) || defined(SERIALIZE_ALL)

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: date
    LOGICAL, PARAMETER :: lonlyonce = .TRUE.
#if defined( SERIALIZE_CREATE_REFERENCE )
    LOGICAL, SAVE :: lenabled = .TRUE.
#else
    LOGICAL, SAVE :: lenabled = .FALSE.
#endif
    LOGICAL, SAVE :: lactive = .FALSE.

    !$ser verbatim IF (selected_block < 0 .OR. jb == selected_block) THEN
    !$ser verbatim   lactive = .TRUE.
    !$ser verbatim ENDIF
    !$ser verbatim IF (lenabled .and. lactive) THEN
    !$ser verbatim   CALL warning('SER:mo_ser_echam_mig:after_satad1','Serialization is active!')
#if defined( _OPENACC )
    !$ser verbatim   CALL warning('GPU:mo_ser_echam_mig:after_satad1','GPU HOST synchronization forced by serialization!')
    !$ser verbatim   !$ACC UPDATE HOST( field%presm_old ) IF( ASSOCIATED(field%presm_old) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rho ) IF( ASSOCIATED(field%rho) )
    !$ser verbatim   !$ACC UPDATE HOST( xlta, xlqv, xlqc, xlqi, xlqr, xlqs, xlqg, zqnc )
    !$ser verbatim   !$ACC UPDATE HOST( zprec_r, zprec_s, zprec_g, zqrsflux )
#endif
    !$ser verbatim   CALL datetimeToString(time_config%tc_current_date, date)
    !$ser verbatim   CALL init('icon')
    !$ser savepoint echam_mig-after-satad1 jg=jg jb=jb jcs=jcs jce=jce nproma=nproma nlev=nlev date=TRIM(date)
    #if defined( SERIALIZE_CREATE_REFERENCE )
    !$ser mode write
#elif defined( SERIALIZE_PERTURB_REFERENCE )
    !$ser mode read-perturb
#elif defined( SERIALIZE_READ_REFERENCE )
    !$ser mode read
#else
#error SERIALIZATION MODE IS NOT SET
#endif
    !$ser data echam_mig_ta=xlta
    !$ser data echam_mig_presm_old=field%presm_old(:,:,jb) IF (ASSOCIATED(field%presm_old))
    !$ser data echam_mig_rho=field%rho(:,:,jb) IF (ASSOCIATED(field%rho))
    !$ser data echam_mig_qv=xlqv
    !$ser data echam_mig_qc=xlqc
    !$ser data echam_mig_qi=xlqi
    !$ser data echam_mig_qr=xlqr
    !$ser data echam_mig_qs=xlqs
    !$ser data echam_mig_qg=xlqg
    !$ser data echam_mig_qnc=zqnc
    !$ser data echam_mig_zprec_r=zprec_r
    !$ser data echam_mig_zprec_s=zprec_s
    !$ser data echam_mig_zprec_g=zprec_g
    !$ser data echam_mig_zqrsflux=zqrsflux
    !$ser verbatim lactive = .FALSE.
    !$ser verbatim IF (lonlyonce) THEN
    !$ser verbatim   lenabled = .FALSE.
    !$ser verbatim END IF
#if defined( SERIALIZE_READ_REFERENCE ) || defined( SERIALIZE_PERTURB_REFERENCE )
#if defined( _OPENACC )
    !$ser verbatim CALL warning('GPU:mo_ser_echam_mig:after_satad1','GPU DEVICE synchronization forced by serialization!')
    !$ser verbatim   !$ACC UPDATE HOST( field%presm_old ) IF( ASSOCIATED(field%presm_old) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rho ) IF( ASSOCIATED(field%rho) )
    !$ser verbatim   !$ACC UPDATE HOST( xlta, xlqv, xlqc, xlqi, xlqr, xlqs, xlqg, zqnc )
    !$ser verbatim   !$ACC UPDATE HOST( zprec_r, zprec_s, zprec_g, zqrsflux )
#endif
#endif
    !$ser verbatim ENDIF

#endif
  END SUBROUTINE serialize_mig_after_satad1

  SUBROUTINE serialize_mig_before_satad2(jg, jb, jcs, jce, nproma, nlev, field, &
      &           xlta, xlqv, xlqc, xlqi, xlqr, xlqs, xlqg, zqnc, &
      &           zprec_r, zprec_s, zprec_g, zqrsflux)
    INTEGER, INTENT(IN)    :: jg, jb, jcs, jce, nproma, nlev
    TYPE(t_echam_phy_field), POINTER, INTENT(INOUT) :: field
    REAL(wp), DIMENSION(nproma, nlev), INTENT(INOUT) :: xlta, xlqv, xlqc, xlqi, &
      &                                                 xlqr, xlqs, xlqg, zqnc
    REAL(wp), DIMENSION(nproma, nlev), INTENT(INOUT) :: zprec_r, zprec_s, zprec_g, zqrsflux
#if defined(SERIALIZE_ECHAM_MIG) || defined(SERIALIZE_ECHAM_ALL) || defined(SERIALIZE_ALL)

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: date
    LOGICAL, PARAMETER :: lonlyonce = .TRUE.
#if defined( SERIALIZE_CREATE_REFERENCE )
    LOGICAL, SAVE :: lenabled = .TRUE.
#else
    LOGICAL, SAVE :: lenabled = .FALSE.
#endif
    LOGICAL, SAVE :: lactive = .FALSE.

    !$ser verbatim IF (selected_block < 0 .OR. jb == selected_block) THEN
    !$ser verbatim   lactive = .TRUE.
    !$ser verbatim ENDIF
    !$ser verbatim IF (lenabled .and. lactive) THEN
    !$ser verbatim   CALL warning('SER:mo_ser_echam_mig:before_satad2','Serialization is active!')
#if defined( _OPENACC )
    !$ser verbatim   CALL warning('GPU:mo_ser_echam_mig:before_satad2','GPU HOST synchronization forced by serialization!')
    !$ser verbatim   !$ACC UPDATE HOST( field%presm_old ) IF( ASSOCIATED(field%presm_old) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rho ) IF( ASSOCIATED(field%rho) )
    !$ser verbatim   !$ACC UPDATE HOST( xlta, xlqv, xlqc, xlqi, xlqr, xlqs, xlqg, zqnc )
    !$ser verbatim   !$ACC UPDATE HOST( zprec_r, zprec_s, zprec_g, zqrsflux )
#endif
    !$ser verbatim   CALL datetimeToString(time_config%tc_current_date, date)
    !$ser verbatim   CALL init('icon')
    !$ser savepoint echam_mig-before-satad2 jg=jg jb=jb jcs=jcs jce=jce nproma=nproma nlev=nlev date=TRIM(date)
    #if defined( SERIALIZE_CREATE_REFERENCE )
    !$ser mode write
#elif defined( SERIALIZE_PERTURB_REFERENCE )
    !$ser mode read-perturb
#elif defined( SERIALIZE_READ_REFERENCE )
    !$ser mode read
#else
#error SERIALIZATION MODE IS NOT SET
#endif
    !$ser data echam_mig_ta=xlta
    !$ser data echam_mig_presm_old=field%presm_old(:,:,jb) IF (ASSOCIATED(field%presm_old))
    !$ser data echam_mig_rho=field%rho(:,:,jb) IF (ASSOCIATED(field%rho))
    !$ser data echam_mig_qv=xlqv
    !$ser data echam_mig_qc=xlqc
    !$ser data echam_mig_qi=xlqi
    !$ser data echam_mig_qr=xlqr
    !$ser data echam_mig_qs=xlqs
    !$ser data echam_mig_qg=xlqg
    !$ser data echam_mig_qnc=zqnc
    !$ser data echam_mig_zprec_r=zprec_r
    !$ser data echam_mig_zprec_s=zprec_s
    !$ser data echam_mig_zprec_g=zprec_g
    !$ser data echam_mig_zqrsflux=zqrsflux
    !$ser verbatim lactive = .FALSE.
    !$ser verbatim IF (lonlyonce) THEN
    !$ser verbatim   lenabled = .FALSE.
    !$ser verbatim END IF
#if defined( SERIALIZE_READ_REFERENCE ) || defined( SERIALIZE_PERTURB_REFERENCE )
#if defined( _OPENACC )
    !$ser verbatim CALL warning('GPU:mo_ser_echam_mig:before_satad2','GPU DEVICE synchronization forced by serialization!')
    !$ser verbatim   !$ACC UPDATE HOST( field%presm_old ) IF( ASSOCIATED(field%presm_old) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rho ) IF( ASSOCIATED(field%rho) )
    !$ser verbatim   !$ACC UPDATE HOST( xlta, xlqv, xlqc, xlqi, xlqr, xlqs, xlqg, zqnc )
    !$ser verbatim   !$ACC UPDATE HOST( zprec_r, zprec_s, zprec_g, zqrsflux )
#endif
#endif
    !$ser verbatim ENDIF

#endif
  END SUBROUTINE serialize_mig_before_satad2

  SUBROUTINE serialize_mig_output(jg, jb, jcs, jce, nproma, nlev, field)
    INTEGER, INTENT(IN)    :: jg, jb, jcs, jce, nproma, nlev
    TYPE(t_echam_phy_field), POINTER, INTENT(INOUT) :: field
#if defined(SERIALIZE_ECHAM_MIG) || defined(SERIALIZE_ECHAM_ALL) || defined(SERIALIZE_ALL)

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: date
    LOGICAL, PARAMETER :: lonlyonce = .TRUE.
    LOGICAL, SAVE :: lenabled = .TRUE.
    LOGICAL, SAVE :: lactive = .FALSE.

    !$ser verbatim IF (selected_block < 0 .OR. jb == selected_block) THEN
    !$ser verbatim   lactive = .TRUE.
    !$ser verbatim ENDIF
    !$ser verbatim IF (lenabled .and. lactive) THEN
    !$ser verbatim   CALL warning('SER:mo_ser_echam_mig:output','Serialization is active!')
#if defined( _OPENACC )
    !$ser verbatim   CALL warning('GPU:mo_ser_echam_mig:output','GPU HOST synchronization forced by serialization!')
    !$ser verbatim   !$ACC UPDATE HOST( field%ta ) IF( ASSOCIATED(field%ta) )
    !$ser verbatim   !$ACC UPDATE HOST( field%presm_old ) IF( ASSOCIATED(field%presm_old) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rho ) IF( ASSOCIATED(field%rho) )
    !$ser verbatim   !$ACC UPDATE HOST( field%qtrc ) IF( ASSOCIATED(field%qtrc) )
#endif
    !$ser verbatim   CALL datetimeToString(time_config%tc_current_date, date)
    !$ser verbatim   CALL init('icon')
    !$ser savepoint echam_mig-output jg=jg jb=jb jcs=jcs jce=jce nproma=nproma nlev=nlev date=TRIM(date)
    !$ser mode write
    !$ser data echam_mig_ta=field%ta(:,:,jb) IF (ASSOCIATED(field%ta))
    !$ser data echam_mig_presm_old=field%presm_old(:,:,jb) IF (ASSOCIATED(field%presm_old))
    !$ser data echam_mig_rho=field%rho(:,:,jb) IF (ASSOCIATED(field%rho))
    !$ser data echam_mig_qv=field%qtrc(:,:,jb,iqv) IF (ASSOCIATED(field%qtrc))
    !$ser data echam_mig_qc=field%qtrc(:,:,jb,iqc) IF (ASSOCIATED(field%qtrc))
    !$ser data echam_mig_qi=field%qtrc(:,:,jb,iqi) IF (ASSOCIATED(field%qtrc))
    !$ser data echam_mig_qr=field%qtrc(:,:,jb,iqr) IF (ASSOCIATED(field%qtrc))
    !$ser data echam_mig_qs=field%qtrc(:,:,jb,iqs) IF (ASSOCIATED(field%qtrc))
    !$ser data echam_mig_qg=field%qtrc(:,:,jb,iqg) IF (ASSOCIATED(field%qtrc))
    !$ser verbatim lactive = .FALSE.
    !$ser verbatim IF (lonlyonce) THEN
    !$ser verbatim   lenabled = .FALSE.
    !$ser verbatim END IF
    !$ser verbatim ENDIF

#endif
  END SUBROUTINE serialize_mig_output

END MODULE mo_ser_echam_mig
