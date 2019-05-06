!--------------------------------------------------------------------
!
! Serialization routine for ECHAM GWD
!
!--------------------------------------------------------------------

MODULE mo_ser_echam_gwd

  USE mo_kind,               ONLY: vp, wp
  USE mo_exception,          ONLY: warning
  USE mo_ser_common,         ONLY: init
  USE mtime,                 ONLY: datetimeToString, MAX_DATETIME_STR_LEN
  USE mo_time_config,        ONLY: time_config
  USE mo_echam_phy_memory,   ONLY: t_echam_phy_field, t_echam_phy_tend
  IMPLICIT NONE

  INTEGER, PARAMETER :: selected_block = 1

  PUBLIC :: serialize_gwd_input
  PUBLIC :: serialize_gwd_output

  CONTAINS

  SUBROUTINE serialize_gwd_input(jg, jb, jcs, jce, nproma, nlev, field, tend)
    INTEGER, INTENT(IN)    :: jg, jb, jcs, jce, nproma, nlev
    TYPE(t_echam_phy_field), POINTER, INTENT(INOUT) :: field
    TYPE(t_echam_phy_tend), POINTER, INTENT(INOUT)  :: tend

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: date
    LOGICAL, PARAMETER :: lonlyonce = .TRUE.
    LOGICAL, SAVE :: lenabled = .TRUE.
    LOGICAL, SAVE :: lactive = .FALSE.

    !$ser verbatim IF (selected_block < 0 .OR. jb == selected_block) THEN
    !$ser verbatim   lactive = .TRUE.
    !$ser verbatim ENDIF
    !$ser verbatim IF (lenabled .and. lactive) THEN
    !$ser verbatim   CALL warning('SER:mo_ser_echam_gwd:input','Serialization is active!')
#if defined( _OPENACC )
    !$ser verbatim   CALL warning('GPU:mo_ser_echam_gwd:input','GPU HOST synchronization forced by serialization!')
    !$ser verbatim   !$ACC UPDATE HOST( field%presi_old ) IF( ASSOCIATED(field%presi_old) )
    !$ser verbatim   !$ACC UPDATE HOST( field%presm_old ) IF( ASSOCIATED(field%presm_old) )
    !$ser verbatim   !$ACC UPDATE HOST( field%zh ) IF( ASSOCIATED(field%zh) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rho ) IF( ASSOCIATED(field%rho) )
    !$ser verbatim   !$ACC UPDATE HOST( field%mair ) IF( ASSOCIATED(field%mair) )
    !$ser verbatim   !$ACC UPDATE HOST( field%ta ) IF( ASSOCIATED(field%ta) )
    !$ser verbatim   !$ACC UPDATE HOST( field%ua ) IF( ASSOCIATED(field%ua) )
    !$ser verbatim   !$ACC UPDATE HOST( field%va ) IF( ASSOCIATED(field%va) )
    !$ser verbatim   !$ACC UPDATE HOST( field%q_gwd ) IF( ASSOCIATED(field%q_gwd) )
    !$ser verbatim   !$ACC UPDATE HOST( field%q_gwd_vi ) IF( ASSOCIATED(field%q_gwd_vi) )
    !$ser verbatim   !$ACC UPDATE HOST( field%qconv ) IF( ASSOCIATED(field%qconv) )
    !$ser verbatim   !$ACC UPDATE HOST( field%q_phy ) IF( ASSOCIATED(field%q_phy) )
    !$ser verbatim   !$ACC UPDATE HOST( field%q_phy_vi ) IF( ASSOCIATED(field%q_phy_vi) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%ua_gwd ) IF( ASSOCIATED(tend%ua_gwd) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%va_gwd ) IF( ASSOCIATED(tend%va_gwd) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%ta_gwd ) IF( ASSOCIATED(tend%ta_gwd) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%ta_phy ) IF( ASSOCIATED(tend%ta_phy) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%ua_phy ) IF( ASSOCIATED(tend%ua_phy) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%va_phy ) IF( ASSOCIATED(tend%va_phy) )
#endif
    !$ser verbatim   CALL datetimeToString(time_config%tc_current_date, date)
    !$ser verbatim   CALL init('icon')
    !$ser savepoint echam_gwd-input jg=jg jb=jb jcs=jcs jce=jce nproma=nproma nlev=nlev date=TRIM(date)
#if defined( SERIALIZE_CREATE_REFERENCE )
    !$ser mode write
#elif defined( SERIALIZE_PERTURB_REFERENCE )
    !$ser mode read-perturb
#elif defined( SERIALIZE_READ_REFERENCE )
    !$ser mode read
#else
#error SERIALIZATION MODE IS NOT SET
#endif
    !$ser data echam_gwd_presi_old=field%presi_old(:,:,jb) IF (ASSOCIATED(field%presi_old))
    !$ser data echam_gwd_presm_old=field%presm_old(:,:,jb) IF (ASSOCIATED(field%presm_old))
    !$ser data echam_gwd_zh=field%zh(:,:,jb) IF (ASSOCIATED(field%zh))
    !$ser data echam_gwd_rho=field%rho(:,:,jb) IF (ASSOCIATED(field%rho))
    !$ser data echam_gwd_mair=field%mair(:,:,jb) IF (ASSOCIATED(field%mair))
    !$ser data echam_gwd_ta=field%ta(:,:,jb) IF (ASSOCIATED(field%ta))
    !$ser data echam_gwd_ua=field%ua(:,:,jb) IF (ASSOCIATED(field%ua))
    !$ser data echam_gwd_va=field%va(:,:,jb) IF (ASSOCIATED(field%va))
    !$ser data echam_gwd_q_gwd=field%q_gwd(:,:,jb) IF (ASSOCIATED(field%q_gwd))
    !$ser data echam_gwd_q_gwd_vi=field%q_gwd_vi(:,jb) IF (ASSOCIATED(field%q_gwd_vi))
    !$ser data echam_gwd_qconv=field%qconv(:,:,jb) IF (ASSOCIATED(field%qconv))
    !$ser data echam_gwd_q_phy=field%q_phy(:,:,jb) IF (ASSOCIATED(field%q_phy))
    !$ser data echam_gwd_q_phy_vi=field%q_phy_vi(:,jb) IF (ASSOCIATED(field%q_phy_vi))
    !$ser data echam_gwd_ua_gwd=tend%ua_gwd(:,:,jb) IF (ASSOCIATED(tend%ua_gwd))
    !$ser data echam_gwd_va_gwd=tend%va_gwd(:,:,jb) IF (ASSOCIATED(tend%va_gwd))
    !$ser data echam_gwd_ta_gwd=tend%ta_gwd(:,:,jb) IF (ASSOCIATED(tend%ta_gwd))
    !$ser data echam_gwd_ta_phy=tend%ta_phy(:,:,jb) IF (ASSOCIATED(tend%ta_phy))
    !$ser data echam_gwd_ua_phy=tend%ua_phy(:,:,jb) IF (ASSOCIATED(tend%ua_phy))
    !$ser data echam_gwd_va_phy=tend%va_phy(:,:,jb) IF (ASSOCIATED(tend%va_phy))
    !$ser verbatim lactive = .FALSE.
    !$ser verbatim IF (lonlyonce) THEN
    !$ser verbatim   lenabled = .FALSE.
    !$ser verbatim END IF
#if defined( SERIALIZE_READ_REFERENCE ) || defined( SERIALIZE_PERTURB_REFERENCE )
#if defined( _OPENACC )
    !$ser verbatim CALL warning('GPU:mo_ser_echam_gwd:input','GPU DEVICE synchronization forced by serialization!')
    !$ser verbatim   !$ACC UPDATE DEVICE( field%presi_old ) IF( ASSOCIATED(field%presi_old) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%presm_old ) IF( ASSOCIATED(field%presm_old) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%zh ) IF( ASSOCIATED(field%zh) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%rho ) IF( ASSOCIATED(field%rho) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%mair ) IF( ASSOCIATED(field%mair) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%ta ) IF( ASSOCIATED(field%ta) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%ua ) IF( ASSOCIATED(field%ua) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%va ) IF( ASSOCIATED(field%va) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%q_gwd ) IF( ASSOCIATED(field%q_gwd) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%q_gwd_vi ) IF( ASSOCIATED(field%q_gwd_vi) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%qconv ) IF( ASSOCIATED(field%qconv) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%q_phy ) IF( ASSOCIATED(field%q_phy) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%q_phy_vi ) IF( ASSOCIATED(field%q_phy_vi) )
    !$ser verbatim   !$ACC UPDATE DEVICE( tend%ua_gwd ) IF( ASSOCIATED(tend%ua_gwd) )
    !$ser verbatim   !$ACC UPDATE DEVICE( tend%va_gwd ) IF( ASSOCIATED(tend%va_gwd) )
    !$ser verbatim   !$ACC UPDATE DEVICE( tend%ta_gwd ) IF( ASSOCIATED(tend%ta_gwd) )
    !$ser verbatim   !$ACC UPDATE DEVICE( tend%ta_phy ) IF( ASSOCIATED(tend%ta_phy) )
    !$ser verbatim   !$ACC UPDATE DEVICE( tend%ua_phy ) IF( ASSOCIATED(tend%ua_phy) )
    !$ser verbatim   !$ACC UPDATE DEVICE( tend%va_phy ) IF( ASSOCIATED(tend%va_phy) )
#endif
#endif
    !$ser verbatim ENDIF

  END SUBROUTINE serialize_gwd_input

  SUBROUTINE serialize_gwd_output(jg, jb, jcs, jce, nproma, nlev, field, tend)
    INTEGER, INTENT(IN)    :: jg, jb, jcs, jce, nproma, nlev
    TYPE(t_echam_phy_field), POINTER, INTENT(INOUT) :: field
    TYPE(t_echam_phy_tend), POINTER, INTENT(INOUT)  :: tend

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: date
    LOGICAL, PARAMETER :: lonlyonce = .TRUE.
    LOGICAL, SAVE :: lenabled = .TRUE.
    LOGICAL, SAVE :: lactive = .FALSE.

    !$ser verbatim IF (selected_block < 0 .OR. jb == selected_block) THEN
    !$ser verbatim   lactive = .TRUE.
    !$ser verbatim ENDIF
    !$ser verbatim IF (lenabled .and. lactive) THEN
    !$ser verbatim   CALL warning('SER:mo_ser_echam_gwd:output','Serialization is active!')
#if defined( _OPENACC )
    !$ser verbatim   CALL warning('GPU:mo_ser_echam_gwd:output','GPU HOST synchronization forced by serialization!')
    !$ser verbatim   !$ACC UPDATE HOST( field%presi_old ) IF( ASSOCIATED(field%presi_old) )
    !$ser verbatim   !$ACC UPDATE HOST( field%presm_old ) IF( ASSOCIATED(field%presm_old) )
    !$ser verbatim   !$ACC UPDATE HOST( field%zh ) IF( ASSOCIATED(field%zh) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rho ) IF( ASSOCIATED(field%rho) )
    !$ser verbatim   !$ACC UPDATE HOST( field%mair ) IF( ASSOCIATED(field%mair) )
    !$ser verbatim   !$ACC UPDATE HOST( field%ta ) IF( ASSOCIATED(field%ta) )
    !$ser verbatim   !$ACC UPDATE HOST( field%ua ) IF( ASSOCIATED(field%ua) )
    !$ser verbatim   !$ACC UPDATE HOST( field%va ) IF( ASSOCIATED(field%va) )
    !$ser verbatim   !$ACC UPDATE HOST( field%q_gwd ) IF( ASSOCIATED(field%q_gwd) )
    !$ser verbatim   !$ACC UPDATE HOST( field%q_gwd_vi ) IF( ASSOCIATED(field%q_gwd_vi) )
    !$ser verbatim   !$ACC UPDATE HOST( field%qconv ) IF( ASSOCIATED(field%qconv) )
    !$ser verbatim   !$ACC UPDATE HOST( field%q_phy ) IF( ASSOCIATED(field%q_phy) )
    !$ser verbatim   !$ACC UPDATE HOST( field%q_phy_vi ) IF( ASSOCIATED(field%q_phy_vi) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%ua_gwd ) IF( ASSOCIATED(tend%ua_gwd) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%va_gwd ) IF( ASSOCIATED(tend%va_gwd) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%ta_gwd ) IF( ASSOCIATED(tend%ta_gwd) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%ta_phy ) IF( ASSOCIATED(tend%ta_phy) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%ua_phy ) IF( ASSOCIATED(tend%ua_phy) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%va_phy ) IF( ASSOCIATED(tend%va_phy) )
#endif
    !$ser verbatim   CALL datetimeToString(time_config%tc_current_date, date)
    !$ser verbatim   CALL init('icon')
    !$ser savepoint echam_gwd-output jg=jg jb=jb jcs=jcs jce=jce nproma=nproma nlev=nlev date=TRIM(date) 
    !$ser mode write
    !$ser data echam_gwd_presi_old=field%presi_old(:,:,jb) IF (ASSOCIATED(field%presi_old))
    !$ser data echam_gwd_presm_old=field%presm_old(:,:,jb) IF (ASSOCIATED(field%presm_old))
    !$ser data echam_gwd_zh=field%zh(:,:,jb) IF (ASSOCIATED(field%zh))
    !$ser data echam_gwd_rho=field%rho(:,:,jb) IF (ASSOCIATED(field%rho))
    !$ser data echam_gwd_mair=field%mair(:,:,jb) IF (ASSOCIATED(field%mair))
    !$ser data echam_gwd_ta=field%ta(:,:,jb) IF (ASSOCIATED(field%ta))
    !$ser data echam_gwd_ua=field%ua(:,:,jb) IF (ASSOCIATED(field%ua))
    !$ser data echam_gwd_va=field%va(:,:,jb) IF (ASSOCIATED(field%va))
    !$ser data echam_gwd_q_gwd=field%q_gwd(:,:,jb) IF (ASSOCIATED(field%q_gwd))
    !$ser data echam_gwd_q_gwd_vi=field%q_gwd_vi(:,jb) IF (ASSOCIATED(field%q_gwd_vi))
    !$ser data echam_gwd_qconv=field%qconv(:,:,jb) IF (ASSOCIATED(field%qconv))
    !$ser data echam_gwd_q_phy=field%q_phy(:,:,jb) IF (ASSOCIATED(field%q_phy))
    !$ser data echam_gwd_q_phy_vi=field%q_phy_vi(:,jb) IF (ASSOCIATED(field%q_phy_vi))
    !$ser data echam_gwd_ua_gwd=tend%ua_gwd(:,:,jb) IF (ASSOCIATED(tend%ua_gwd))
    !$ser data echam_gwd_va_gwd=tend%va_gwd(:,:,jb) IF (ASSOCIATED(tend%va_gwd))
    !$ser data echam_gwd_ta_gwd=tend%ta_gwd(:,:,jb) IF (ASSOCIATED(tend%ta_gwd))
    !$ser data echam_gwd_ta_phy=tend%ta_phy(:,:,jb) IF (ASSOCIATED(tend%ta_phy))
    !$ser data echam_gwd_ua_phy=tend%ua_phy(:,:,jb) IF (ASSOCIATED(tend%ua_phy))
    !$ser data echam_gwd_va_phy=tend%va_phy(:,:,jb) IF (ASSOCIATED(tend%va_phy))
    !$ser verbatim lactive = .FALSE.
    !$ser verbatim IF (lonlyonce) THEN
    !$ser verbatim   lenabled = .FALSE.
    !$ser verbatim END IF
    !$ser verbatim ENDIF

  END SUBROUTINE serialize_gwd_output

END MODULE mo_ser_echam_gwd
