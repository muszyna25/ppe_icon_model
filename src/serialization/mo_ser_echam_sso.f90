!--------------------------------------------------------------------
!
! Serialization routine for ECHAM SSO
!
!--------------------------------------------------------------------

MODULE mo_ser_echam_sso

  USE mo_kind,               ONLY: vp, wp
  USE mo_exception,          ONLY: warning
  USE mo_ser_common,         ONLY: init
  USE mtime,                 ONLY: datetimeToString, MAX_DATETIME_STR_LEN
  USE mo_time_config,        ONLY: time_config
  USE mo_echam_phy_memory,   ONLY: t_echam_phy_field, t_echam_phy_tend
  IMPLICIT NONE

  INTEGER, PARAMETER :: selected_block = 1

  PUBLIC :: serialize_sso_input
  PUBLIC :: serialize_sso_output

  CONTAINS

  SUBROUTINE serialize_sso_input(jg, jb, jcs, jce, nproma, nlev, field, tend)
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
    !$ser verbatim   CALL warning('SER:mo_ser_echam_sso:input','Serialization is active!')
#if defined( _OPENACC )
    !$ser verbatim   CALL warning('GPU:mo_ser_echam_sso:input','GPU HOST synchronization forced by serialization!')
    !$ser verbatim   !$ACC UPDATE HOST( field%sftlf ) IF( ASSOCIATED(field%sftlf) )
    !$ser verbatim   !$ACC UPDATE HOST( field%coriol ) IF( ASSOCIATED(field%coriol) )
    !$ser verbatim   !$ACC UPDATE HOST( field%zf ) IF( ASSOCIATED(field%zf) )
    !$ser verbatim   !$ACC UPDATE HOST( field%zh ) IF( ASSOCIATED(field%zh) )
    !$ser verbatim   !$ACC UPDATE HOST( field%presi_old ) IF( ASSOCIATED(field%presi_old) )
    !$ser verbatim   !$ACC UPDATE HOST( field%presm_old ) IF( ASSOCIATED(field%presm_old) )
    !$ser verbatim   !$ACC UPDATE HOST( field%mair ) IF( ASSOCIATED(field%mair) )
    !$ser verbatim   !$ACC UPDATE HOST( field%ta ) IF( ASSOCIATED(field%ta) )
    !$ser verbatim   !$ACC UPDATE HOST( field%ua ) IF( ASSOCIATED(field%ua) )
    !$ser verbatim   !$ACC UPDATE HOST( field%va ) IF( ASSOCIATED(field%va) )
    !$ser verbatim   !$ACC UPDATE HOST( field%oromea ) IF( ASSOCIATED(field%oromea) )
    !$ser verbatim   !$ACC UPDATE HOST( field%orostd ) IF( ASSOCIATED(field%orostd) )
    !$ser verbatim   !$ACC UPDATE HOST( field%orosig ) IF( ASSOCIATED(field%orosig) )
    !$ser verbatim   !$ACC UPDATE HOST( field%orogam ) IF( ASSOCIATED(field%orogam) )
    !$ser verbatim   !$ACC UPDATE HOST( field%orothe ) IF( ASSOCIATED(field%orothe) )
    !$ser verbatim   !$ACC UPDATE HOST( field%oropic ) IF( ASSOCIATED(field%oropic) )
    !$ser verbatim   !$ACC UPDATE HOST( field%oroval ) IF( ASSOCIATED(field%oroval) )
    !$ser verbatim   !$ACC UPDATE HOST( field%u_stress_sso ) IF( ASSOCIATED(field%u_stress_sso) )
    !$ser verbatim   !$ACC UPDATE HOST( field%v_stress_sso ) IF( ASSOCIATED(field%v_stress_sso) )
    !$ser verbatim   !$ACC UPDATE HOST( field%dissipation_sso ) IF( ASSOCIATED(field%dissipation_sso) )
    !$ser verbatim   !$ACC UPDATE HOST( field%q_sso ) IF( ASSOCIATED(field%q_sso) )
    !$ser verbatim   !$ACC UPDATE HOST( field%q_sso_vi ) IF( ASSOCIATED(field%q_sso_vi) )
    !$ser verbatim   !$ACC UPDATE HOST( field%qconv ) IF( ASSOCIATED(field%qconv) )
    !$ser verbatim   !$ACC UPDATE HOST( field%q_phy ) IF( ASSOCIATED(field%q_phy) )
    !$ser verbatim   !$ACC UPDATE HOST( field%q_phy_vi ) IF( ASSOCIATED(field%q_phy_vi) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%ua_sso ) IF( ASSOCIATED(tend%ua_sso) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%va_sso ) IF( ASSOCIATED(tend%va_sso) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%ta_sso ) IF( ASSOCIATED(tend%ta_sso) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%ua_phy ) IF( ASSOCIATED(tend%ua_phy) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%va_phy ) IF( ASSOCIATED(tend%va_phy) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%ta_phy ) IF( ASSOCIATED(tend%ta_phy) )
#endif
    !$ser verbatim   CALL datetimeToString(time_config%tc_current_date, date)
    !$ser verbatim   CALL init('icon')
    !$ser savepoint echam_sso-input jg=jg jb=jb jcs=jcs jce=jce nproma=nproma nlev=nlev date=TRIM(date)
#if defined( SERIALIZE_CREATE_REFERENCE )
    !$ser mode write
#elif defined( SERIALIZE_PERTURB_REFERENCE )
    !$ser mode read-perturb
#elif defined( SERIALIZE_READ_REFERENCE )
    !$ser mode read
#else
#error SERIALIZATION MODE IS NOT SET
#endif
    !$ser data echam_sso_sftlf=field%sftlf(:,jb) IF (ASSOCIATED(field%sftlf))
    !$ser data echam_sso_coriol=field%coriol(:,jb) IF (ASSOCIATED(field%coriol))
    !$ser data echam_sso_zf=field%zf(:,:,jb) IF (ASSOCIATED(field%zf))
    !$ser data echam_sso_zh=field%zh(:,:,jb) IF (ASSOCIATED(field%zh))
    !$ser data echam_sso_presi_old=field%presi_old(:,:,jb) IF (ASSOCIATED(field%presi_old))
    !$ser data echam_sso_presm_old=field%presm_old(:,:,jb) IF (ASSOCIATED(field%presm_old))
    !$ser data echam_sso_mair=field%mair(:,:,jb) IF (ASSOCIATED(field%mair))
    !$ser data echam_sso_ta=field%ta(:,:,jb) IF (ASSOCIATED(field%ta))
    !$ser data echam_sso_ua=field%ua(:,:,jb) IF (ASSOCIATED(field%ua))
    !$ser data echam_sso_va=field%va(:,:,jb) IF (ASSOCIATED(field%va))
    !$ser data echam_sso_oromea=field%oromea(:,jb) IF (ASSOCIATED(field%oromea))
    !$ser data echam_sso_orostd=field%orostd(:,jb) IF (ASSOCIATED(field%orostd))
    !$ser data echam_sso_orosig=field%orosig(:,jb) IF (ASSOCIATED(field%orosig))
    !$ser data echam_sso_orogam=field%orogam(:,jb) IF (ASSOCIATED(field%orogam))
    !$ser data echam_sso_orothe=field%orothe(:,jb) IF (ASSOCIATED(field%orothe))
    !$ser data echam_sso_oropic=field%oropic(:,jb) IF (ASSOCIATED(field%oropic))
    !$ser data echam_sso_oroval=field%oroval(:,jb) IF (ASSOCIATED(field%oroval))
    !$ser data echam_sso_u_stress_sso=field%u_stress_sso(:,jb) IF (ASSOCIATED(field%u_stress_sso))
    !$ser data echam_sso_v_stress_sso=field%v_stress_sso(:,jb) IF (ASSOCIATED(field%v_stress_sso))
    !$ser data echam_sso_dissipation_sso=field%dissipation_sso(:,jb) IF (ASSOCIATED(field%dissipation_sso))
    !$ser data echam_sso_q_sso=field%q_sso(:,:,jb) IF (ASSOCIATED(field%q_sso))
    !$ser data echam_sso_q_sso_vi=field%q_sso_vi(:,jb) IF (ASSOCIATED(field%q_sso_vi))
    !$ser data echam_sso_qconv=field%qconv(:,:,jb) IF (ASSOCIATED(field%qconv))
    !$ser data echam_sso_q_phy=field%q_phy(:,:,jb) IF (ASSOCIATED(field%q_phy))
    !$ser data echam_sso_q_phy_vi=field%q_phy_vi(:,jb) IF (ASSOCIATED(field%q_phy_vi))
    !$ser data echam_sso_ua_sso=tend%ua_sso(:,:,jb) IF (ASSOCIATED(tend%ua_sso))
    !$ser data echam_sso_va_sso=tend%va_sso(:,:,jb) IF (ASSOCIATED(tend%va_sso))
    !$ser data echam_sso_ta_sso=tend%ta_sso(:,:,jb) IF (ASSOCIATED(tend%ta_sso))
    !$ser data echam_sso_ua_phy=tend%ua_phy(:,:,jb) IF (ASSOCIATED(tend%ua_phy))
    !$ser data echam_sso_va_phy=tend%va_phy(:,:,jb) IF (ASSOCIATED(tend%va_phy))
    !$ser data echam_sso_ta_phy=tend%ta_phy(:,:,jb) IF (ASSOCIATED(tend%ta_phy))
    !$ser verbatim lactive = .FALSE.
    !$ser verbatim IF (lonlyonce) THEN
    !$ser verbatim   lenabled = .FALSE.
    !$ser verbatim END IF
#if defined( SERIALIZE_READ_REFERENCE ) || defined( SERIALIZE_PERTURB_REFERENCE )
#if defined( _OPENACC )
    !$ser verbatim CALL warning('GPU:mo_ser_echam_sso:input','GPU DEVICE synchronization forced by serialization!')
    !$ser verbatim   !$ACC UPDATE DEVICE( field%sftlf ) IF( ASSOCIATED(field%sftlf) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%coriol ) IF( ASSOCIATED(field%coriol) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%zf ) IF( ASSOCIATED(field%zf) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%zh ) IF( ASSOCIATED(field%zh) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%presi_old ) IF( ASSOCIATED(field%presi_old) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%presm_old ) IF( ASSOCIATED(field%presm_old) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%mair ) IF( ASSOCIATED(field%mair) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%ta ) IF( ASSOCIATED(field%ta) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%ua ) IF( ASSOCIATED(field%ua) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%va ) IF( ASSOCIATED(field%va) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%oromea ) IF( ASSOCIATED(field%oromea) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%orostd ) IF( ASSOCIATED(field%orostd) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%orosig ) IF( ASSOCIATED(field%orosig) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%orogam ) IF( ASSOCIATED(field%orogam) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%orothe ) IF( ASSOCIATED(field%orothe) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%oropic ) IF( ASSOCIATED(field%oropic) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%oroval ) IF( ASSOCIATED(field%oroval) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%u_stress_sso ) IF( ASSOCIATED(field%u_stress_sso) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%v_stress_sso ) IF( ASSOCIATED(field%v_stress_sso) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%dissipation_sso ) IF( ASSOCIATED(field%dissipation_sso) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%q_sso ) IF( ASSOCIATED(field%q_sso) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%q_sso_vi ) IF( ASSOCIATED(field%q_sso_vi) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%qconv ) IF( ASSOCIATED(field%qconv) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%q_phy ) IF( ASSOCIATED(field%q_phy) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%q_phy_vi ) IF( ASSOCIATED(field%q_phy_vi) )
    !$ser verbatim   !$ACC UPDATE DEVICE( tend%ua_sso ) IF( ASSOCIATED(tend%ua_sso) )
    !$ser verbatim   !$ACC UPDATE DEVICE( tend%va_sso ) IF( ASSOCIATED(tend%va_sso) )
    !$ser verbatim   !$ACC UPDATE DEVICE( tend%ta_sso ) IF( ASSOCIATED(tend%ta_sso) )
    !$ser verbatim   !$ACC UPDATE DEVICE( tend%ua_phy ) IF( ASSOCIATED(tend%ua_phy) )
    !$ser verbatim   !$ACC UPDATE DEVICE( tend%va_phy ) IF( ASSOCIATED(tend%va_phy) )
    !$ser verbatim   !$ACC UPDATE DEVICE( tend%ta_phy ) IF( ASSOCIATED(tend%ta_phy) )
#endif
#endif
    !$ser verbatim ENDIF

  END SUBROUTINE serialize_sso_input

  SUBROUTINE serialize_sso_output(jg, jb, jcs, jce, nproma, nlev, field, tend)
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
    !$ser verbatim   CALL warning('SER:mo_ser_echam_sso:output','Serialization is active!')
#if defined( _OPENACC )
    !$ser verbatim   CALL warning('GPU:mo_ser_echam_sso:output','GPU HOST synchronization forced by serialization!')
    !$ser verbatim   !$ACC UPDATE HOST( field%sftlf ) IF( ASSOCIATED(field%sftlf) )
    !$ser verbatim   !$ACC UPDATE HOST( field%coriol ) IF( ASSOCIATED(field%coriol) )
    !$ser verbatim   !$ACC UPDATE HOST( field%zf ) IF( ASSOCIATED(field%zf) )
    !$ser verbatim   !$ACC UPDATE HOST( field%zh ) IF( ASSOCIATED(field%zh) )
    !$ser verbatim   !$ACC UPDATE HOST( field%presi_old ) IF( ASSOCIATED(field%presi_old) )
    !$ser verbatim   !$ACC UPDATE HOST( field%presm_old ) IF( ASSOCIATED(field%presm_old) )
    !$ser verbatim   !$ACC UPDATE HOST( field%mair ) IF( ASSOCIATED(field%mair) )
    !$ser verbatim   !$ACC UPDATE HOST( field%ta ) IF( ASSOCIATED(field%ta) )
    !$ser verbatim   !$ACC UPDATE HOST( field%ua ) IF( ASSOCIATED(field%ua) )
    !$ser verbatim   !$ACC UPDATE HOST( field%va ) IF( ASSOCIATED(field%va) )
    !$ser verbatim   !$ACC UPDATE HOST( field%oromea ) IF( ASSOCIATED(field%oromea) )
    !$ser verbatim   !$ACC UPDATE HOST( field%orostd ) IF( ASSOCIATED(field%orostd) )
    !$ser verbatim   !$ACC UPDATE HOST( field%orosig ) IF( ASSOCIATED(field%orosig) )
    !$ser verbatim   !$ACC UPDATE HOST( field%orogam ) IF( ASSOCIATED(field%orogam) )
    !$ser verbatim   !$ACC UPDATE HOST( field%orothe ) IF( ASSOCIATED(field%orothe) )
    !$ser verbatim   !$ACC UPDATE HOST( field%oropic ) IF( ASSOCIATED(field%oropic) )
    !$ser verbatim   !$ACC UPDATE HOST( field%oroval ) IF( ASSOCIATED(field%oroval) )
    !$ser verbatim   !$ACC UPDATE HOST( field%u_stress_sso ) IF( ASSOCIATED(field%u_stress_sso) )
    !$ser verbatim   !$ACC UPDATE HOST( field%v_stress_sso ) IF( ASSOCIATED(field%v_stress_sso) )
    !$ser verbatim   !$ACC UPDATE HOST( field%dissipation_sso ) IF( ASSOCIATED(field%dissipation_sso) )
    !$ser verbatim   !$ACC UPDATE HOST( field%q_sso ) IF( ASSOCIATED(field%q_sso) )
    !$ser verbatim   !$ACC UPDATE HOST( field%q_sso_vi ) IF( ASSOCIATED(field%q_sso_vi) )
    !$ser verbatim   !$ACC UPDATE HOST( field%qconv ) IF( ASSOCIATED(field%qconv) )
    !$ser verbatim   !$ACC UPDATE HOST( field%q_phy ) IF( ASSOCIATED(field%q_phy) )
    !$ser verbatim   !$ACC UPDATE HOST( field%q_phy_vi ) IF( ASSOCIATED(field%q_phy_vi) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%ua_sso ) IF( ASSOCIATED(tend%ua_sso) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%va_sso ) IF( ASSOCIATED(tend%va_sso) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%ta_sso ) IF( ASSOCIATED(tend%ta_sso) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%ua_phy ) IF( ASSOCIATED(tend%ua_phy) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%va_phy ) IF( ASSOCIATED(tend%va_phy) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%ta_phy ) IF( ASSOCIATED(tend%ta_phy) )
#endif
    !$ser verbatim   CALL datetimeToString(time_config%tc_current_date, date)
    !$ser verbatim   CALL init('icon')
    !$ser savepoint echam_sso-output jg=jg jb=jb jcs=jcs jce=jce nproma=nproma nlev=nlev date=TRIM(date) 
    !$ser mode write
    !$ser data echam_sso_sftlf=field%sftlf(:,jb) IF (ASSOCIATED(field%sftlf))
    !$ser data echam_sso_coriol=field%coriol(:,jb) IF (ASSOCIATED(field%coriol))
    !$ser data echam_sso_zf=field%zf(:,:,jb) IF (ASSOCIATED(field%zf))
    !$ser data echam_sso_zh=field%zh(:,:,jb) IF (ASSOCIATED(field%zh))
    !$ser data echam_sso_presi_old=field%presi_old(:,:,jb) IF (ASSOCIATED(field%presi_old))
    !$ser data echam_sso_presm_old=field%presm_old(:,:,jb) IF (ASSOCIATED(field%presm_old))
    !$ser data echam_sso_mair=field%mair(:,:,jb) IF (ASSOCIATED(field%mair))
    !$ser data echam_sso_ta=field%ta(:,:,jb) IF (ASSOCIATED(field%ta))
    !$ser data echam_sso_ua=field%ua(:,:,jb) IF (ASSOCIATED(field%ua))
    !$ser data echam_sso_va=field%va(:,:,jb) IF (ASSOCIATED(field%va))
    !$ser data echam_sso_oromea=field%oromea(:,jb) IF (ASSOCIATED(field%oromea))
    !$ser data echam_sso_orostd=field%orostd(:,jb) IF (ASSOCIATED(field%orostd))
    !$ser data echam_sso_orosig=field%orosig(:,jb) IF (ASSOCIATED(field%orosig))
    !$ser data echam_sso_orogam=field%orogam(:,jb) IF (ASSOCIATED(field%orogam))
    !$ser data echam_sso_orothe=field%orothe(:,jb) IF (ASSOCIATED(field%orothe))
    !$ser data echam_sso_oropic=field%oropic(:,jb) IF (ASSOCIATED(field%oropic))
    !$ser data echam_sso_oroval=field%oroval(:,jb) IF (ASSOCIATED(field%oroval))
    !$ser data echam_sso_u_stress_sso=field%u_stress_sso(:,jb) IF (ASSOCIATED(field%u_stress_sso))
    !$ser data echam_sso_v_stress_sso=field%v_stress_sso(:,jb) IF (ASSOCIATED(field%v_stress_sso))
    !$ser data echam_sso_dissipation_sso=field%dissipation_sso(:,jb) IF (ASSOCIATED(field%dissipation_sso))
    !$ser data echam_sso_q_sso=field%q_sso(:,:,jb) IF (ASSOCIATED(field%q_sso))
    !$ser data echam_sso_q_sso_vi=field%q_sso_vi(:,jb) IF (ASSOCIATED(field%q_sso_vi))
    !$ser data echam_sso_qconv=field%qconv(:,:,jb) IF (ASSOCIATED(field%qconv))
    !$ser data echam_sso_q_phy=field%q_phy(:,:,jb) IF (ASSOCIATED(field%q_phy))
    !$ser data echam_sso_q_phy_vi=field%q_phy_vi(:,jb) IF (ASSOCIATED(field%q_phy_vi))
    !$ser data echam_sso_ua_sso=tend%ua_sso(:,:,jb) IF (ASSOCIATED(tend%ua_sso))
    !$ser data echam_sso_va_sso=tend%va_sso(:,:,jb) IF (ASSOCIATED(tend%va_sso))
    !$ser data echam_sso_ta_sso=tend%ta_sso(:,:,jb) IF (ASSOCIATED(tend%ta_sso))
    !$ser data echam_sso_ua_phy=tend%ua_phy(:,:,jb) IF (ASSOCIATED(tend%ua_phy))
    !$ser data echam_sso_va_phy=tend%va_phy(:,:,jb) IF (ASSOCIATED(tend%va_phy))
    !$ser data echam_sso_ta_phy=tend%ta_phy(:,:,jb) IF (ASSOCIATED(tend%ta_phy))
    !$ser verbatim lactive = .FALSE.
    !$ser verbatim IF (lonlyonce) THEN
    !$ser verbatim   lenabled = .FALSE.
    !$ser verbatim END IF
    !$ser verbatim ENDIF

  END SUBROUTINE serialize_sso_output

END MODULE mo_ser_echam_sso
