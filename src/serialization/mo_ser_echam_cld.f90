!--------------------------------------------------------------------
!
! Serialization routine for ECHAM CLD
!
!--------------------------------------------------------------------

MODULE mo_ser_echam_cld

  USE mo_kind,               ONLY: vp, wp
  USE mo_exception,          ONLY: warning
  USE mo_ser_common,         ONLY: init
  USE mtime,                 ONLY: datetimeToString, MAX_DATETIME_STR_LEN
  USE mo_time_config,        ONLY: time_config
  USE mo_echam_phy_memory,   ONLY: t_echam_phy_field, t_echam_phy_tend
  IMPLICIT NONE

  INTEGER, PARAMETER :: selected_block = 1

  PUBLIC :: serialize_cld_input
  PUBLIC :: serialize_cld_output

  CONTAINS

  SUBROUTINE serialize_cld_input(jg, jb, jcs, jce, nproma, nlev, field, tend)
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
    !$ser verbatim   CALL warning('SER:mo_ser_echam_cld:input','Serialization is active!')
#if defined( _OPENACC )
    !$ser verbatim   CALL warning('GPU:mo_ser_echam_cld:input','GPU HOST synchronization forced by serialization!')
    !$ser verbatim   !$ACC UPDATE HOST( field%rtype ) IF( ASSOCIATED(field%rtype) )
    !$ser verbatim   !$ACC UPDATE HOST( field%ictop ) IF( ASSOCIATED(field%ictop) )
    !$ser verbatim   !$ACC UPDATE HOST( field%presm_old ) IF( ASSOCIATED(field%presm_old) )
    !$ser verbatim   !$ACC UPDATE HOST( field%dz ) IF( ASSOCIATED(field%dz) )
    !$ser verbatim   !$ACC UPDATE HOST( field%mref ) IF( ASSOCIATED(field%mref) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rho ) IF( ASSOCIATED(field%rho) )
    !$ser verbatim   !$ACC UPDATE HOST( field%cpair ) IF( ASSOCIATED(field%cpair) )
    !$ser verbatim   !$ACC UPDATE HOST( field%acdnc ) IF( ASSOCIATED(field%acdnc) )
    !$ser verbatim   !$ACC UPDATE HOST( field%ta ) IF( ASSOCIATED(field%ta) )
    !$ser verbatim   !$ACC UPDATE HOST( field%qtrc ) IF( ASSOCIATED(field%qtrc) )
    !$ser verbatim   !$ACC UPDATE HOST( field%aclc ) IF( ASSOCIATED(field%aclc) )
    !$ser verbatim   !$ACC UPDATE HOST( field%aclcov ) IF( ASSOCIATED(field%aclcov) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rsfl ) IF( ASSOCIATED(field%rsfl) )
    !$ser verbatim   !$ACC UPDATE HOST( field%ssfl ) IF( ASSOCIATED(field%ssfl) )
    !$ser verbatim   !$ACC UPDATE HOST( field%hur ) IF( ASSOCIATED(field%hur) )
    !$ser verbatim   !$ACC UPDATE HOST( field%q_cld ) IF( ASSOCIATED(field%q_cld) )
    !$ser verbatim   !$ACC UPDATE HOST( field%q_cld_vi ) IF( ASSOCIATED(field%q_cld_vi) )
    !$ser verbatim   !$ACC UPDATE HOST( field%qconv ) IF( ASSOCIATED(field%qconv) )
    !$ser verbatim   !$ACC UPDATE HOST( field%q_phy ) IF( ASSOCIATED(field%q_phy) )
    !$ser verbatim   !$ACC UPDATE HOST( field%q_phy_vi ) IF( ASSOCIATED(field%q_phy_vi) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%qtrc_cld ) IF( ASSOCIATED(tend%qtrc_cld) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%ta_cld ) IF( ASSOCIATED(tend%ta_cld) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%ta_phy ) IF( ASSOCIATED(tend%ta_phy) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%qtrc_phy ) IF( ASSOCIATED(tend%qtrc_phy) )
#endif
    !$ser verbatim   CALL datetimeToString(time_config%tc_current_date, date)
    !$ser verbatim   CALL init('icon')
    !$ser savepoint echam_cld-input jg=jg jb=jb jcs=jcs jce=jce nproma=nproma nlev=nlev date=TRIM(date)
#if defined( SERIALIZE_CREATE_REFERENCE )
    !$ser mode write
#elif defined( SERIALIZE_PERTURB_REFERENCE )
    !$ser mode read-perturb
#elif defined( SERIALIZE_READ_REFERENCE )
    !$ser mode read
#else
#error SERIALIZATION MODE IS NOT SET
#endif
    !$ser data echam_cld_rtype=field%rtype(:,jb) IF (ASSOCIATED(field%rtype))
    !$ser data echam_cld_ictop=field%ictop(:,jb) IF (ASSOCIATED(field%ictop))
    !$ser data echam_cld_presm_old=field%presm_old(:,:,jb) IF (ASSOCIATED(field%presm_old))
    !$ser data echam_cld_dz=field%dz(:,:,jb) IF (ASSOCIATED(field%dz))
    !$ser data echam_cld_mref=field%mref(:,:,jb) IF (ASSOCIATED(field%mref))
    !$ser data echam_cld_rho=field%rho(:,:,jb) IF (ASSOCIATED(field%rho))
    !$ser data echam_cld_cpair=field%cpair(:,:,jb) IF (ASSOCIATED(field%cpair))
    !$ser data echam_cld_acdnc=field%acdnc(:,:,jb) IF (ASSOCIATED(field%acdnc))
    !$ser data echam_cld_ta=field%ta(:,:,jb) IF (ASSOCIATED(field%ta))
    !$ser data echam_cld_qtrc=field%qtrc(:,:,jb,:) IF (ASSOCIATED(field%qtrc))
    !$ser data echam_cld_aclc=field%aclc(:,:,jb) IF (ASSOCIATED(field%aclc))
    !$ser data echam_cld_aclcov=field%aclcov(:,jb) IF (ASSOCIATED(field%aclcov))
    !$ser data echam_cld_rsfl=field%rsfl(:,jb) IF (ASSOCIATED(field%rsfl))
    !$ser data echam_cld_ssfl=field%ssfl(:,jb) IF (ASSOCIATED(field%ssfl))
    !$ser data echam_cld_hur=field%hur IF (ASSOCIATED(field%hur))
    !$ser data echam_cld_q_cld=field%q_cld IF (ASSOCIATED(field%q_cld))
    !$ser data echam_cld_q_cld_vi=field%q_cld_vi IF (ASSOCIATED(field%q_cld_vi))  
    !$ser data echam_cld_qconv=field%qconv IF (ASSOCIATED(field%qconv))
    !$ser data echam_cld_q_phy=field%q_phy IF (ASSOCIATED(field%q_phy))
    !$ser data echam_cld_q_phy_vi=field%q_phy_vi IF (ASSOCIATED(field%q_phy_vi))
    !$ser data echam_cld_qtrc_cld=tend%qtrc_cld(:,:,jb,:) IF (ASSOCIATED(tend%qtrc_cld))
    !$ser data echam_cld_ta_cld=tend%ta_cld(:,:,jb) IF (ASSOCIATED(tend%ta_cld))
    !$ser data echam_cld_ta_phy=tend%ta_phy(:,:,jb) IF (ASSOCIATED(tend%ta_phy))
    !$ser data echam_cld_qtrc_phy=tend%qtrc_phy(:,:,jb,:) IF (ASSOCIATED(tend%qtrc_phy))
    !$ser verbatim lactive = .FALSE.
    !$ser verbatim IF (lonlyonce) THEN
    !$ser verbatim   lenabled = .FALSE.
    !$ser verbatim END IF
#if defined( SERIALIZE_READ_REFERENCE ) || defined( SERIALIZE_PERTURB_REFERENCE )
#if defined( _OPENACC )
    !$ser verbatim CALL warning('GPU:mo_ser_echam_cld:input','GPU DEVICE synchronization forced by serialization!')
    !$ser verbatim   !$ACC UPDATE DEVICE( field%rtype ) IF( ASSOCIATED(field%rtype) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%ictop ) IF( ASSOCIATED(field%ictop) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%presm_old ) IF( ASSOCIATED(field%presm_old) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%dz ) IF( ASSOCIATED(field%dz) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%mref ) IF( ASSOCIATED(field%mref) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%rho ) IF( ASSOCIATED(field%rho) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%cpair ) IF( ASSOCIATED(field%cpair) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%acdnc ) IF( ASSOCIATED(field%acdnc) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%ta ) IF( ASSOCIATED(field%ta) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%qtrc ) IF( ASSOCIATED(field%qtrc) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%aclc ) IF( ASSOCIATED(field%aclc) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%aclcov ) IF( ASSOCIATED(field%aclcov) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%rsfl ) IF( ASSOCIATED(field%rsfl) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%ssfl ) IF( ASSOCIATED(field%ssfl) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%hur ) IF( ASSOCIATED(field%hur) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%q_cld ) IF( ASSOCIATED(field%q_cld) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%q_cld_vi ) IF( ASSOCIATED(field%q_cld_vi) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%qconv ) IF( ASSOCIATED(field%qconv) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%q_phy ) IF( ASSOCIATED(field%q_phy) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%q_phy_vi ) IF( ASSOCIATED(field%q_phy_vi) )
    !$ser verbatim   !$ACC UPDATE DEVICE( tend%qtrc_cld ) IF( ASSOCIATED(tend%qtrc_cld) )
    !$ser verbatim   !$ACC UPDATE DEVICE( tend%ta_cld ) IF( ASSOCIATED(tend%ta_cld) )
    !$ser verbatim   !$ACC UPDATE DEVICE( tend%ta_phy ) IF( ASSOCIATED(tend%ta_phy) )
    !$ser verbatim   !$ACC UPDATE DEVICE( tend%qtrc_phy ) IF( ASSOCIATED(tend%qtrc_phy) )
#endif
#endif
    !$ser verbatim ENDIF

  END SUBROUTINE serialize_cld_input

  SUBROUTINE serialize_cld_output(jg, jb, jcs, jce, nproma, nlev, field, tend)
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
    !$ser verbatim   CALL warning('SER:mo_ser_echam_cld:output','Serialization is active!')
#if defined( _OPENACC )
    !$ser verbatim   CALL warning('GPU:mo_ser_echam_cld:output','GPU HOST synchronization forced by serialization!')
    !$ser verbatim   !$ACC UPDATE HOST( field%rtype ) IF( ASSOCIATED(field%rtype) )
    !$ser verbatim   !$ACC UPDATE HOST( field%ictop ) IF( ASSOCIATED(field%ictop) )
    !$ser verbatim   !$ACC UPDATE HOST( field%presm_old ) IF( ASSOCIATED(field%presm_old) )
    !$ser verbatim   !$ACC UPDATE HOST( field%dz ) IF( ASSOCIATED(field%dz) )
    !$ser verbatim   !$ACC UPDATE HOST( field%mref ) IF( ASSOCIATED(field%mref) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rho ) IF( ASSOCIATED(field%rho) )
    !$ser verbatim   !$ACC UPDATE HOST( field%cpair ) IF( ASSOCIATED(field%cpair) )
    !$ser verbatim   !$ACC UPDATE HOST( field%acdnc ) IF( ASSOCIATED(field%acdnc) )
    !$ser verbatim   !$ACC UPDATE HOST( field%ta ) IF( ASSOCIATED(field%ta) )
    !$ser verbatim   !$ACC UPDATE HOST( field%qtrc ) IF( ASSOCIATED(field%qtrc) )
    !$ser verbatim   !$ACC UPDATE HOST( field%aclc ) IF( ASSOCIATED(field%aclc) )
    !$ser verbatim   !$ACC UPDATE HOST( field%aclcov ) IF( ASSOCIATED(field%aclcov) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rsfl ) IF( ASSOCIATED(field%rsfl) )
    !$ser verbatim   !$ACC UPDATE HOST( field%ssfl ) IF( ASSOCIATED(field%ssfl) )
    !$ser verbatim   !$ACC UPDATE HOST( field%hur ) IF( ASSOCIATED(field%hur) )
    !$ser verbatim   !$ACC UPDATE HOST( field%q_cld ) IF( ASSOCIATED(field%q_cld) )
    !$ser verbatim   !$ACC UPDATE HOST( field%q_cld_vi ) IF( ASSOCIATED(field%q_cld_vi) )
    !$ser verbatim   !$ACC UPDATE HOST( field%qconv ) IF( ASSOCIATED(field%qconv) )
    !$ser verbatim   !$ACC UPDATE HOST( field%q_phy ) IF( ASSOCIATED(field%q_phy) )
    !$ser verbatim   !$ACC UPDATE HOST( field%q_phy_vi ) IF( ASSOCIATED(field%q_phy_vi) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%qtrc_cld ) IF( ASSOCIATED(tend%qtrc_cld) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%ta_cld ) IF( ASSOCIATED(tend%ta_cld) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%ta_phy ) IF( ASSOCIATED(tend%ta_phy) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%qtrc_phy ) IF( ASSOCIATED(tend%qtrc_phy) )
#endif
    !$ser verbatim   CALL datetimeToString(time_config%tc_current_date, date)
    !$ser verbatim   CALL init('icon')
    !$ser savepoint echam_cld-output jg=jg jb=jb jcs=jcs jce=jce nproma=nproma nlev=nlev date=TRIM(date) 
    !$ser mode write
    !$ser data echam_cld_rtype=field%rtype(:,jb) IF (ASSOCIATED(field%rtype))
    !$ser data echam_cld_ictop=field%ictop(:,jb) IF (ASSOCIATED(field%ictop))
    !$ser data echam_cld_presm_old=field%presm_old(:,:,jb) IF (ASSOCIATED(field%presm_old))
    !$ser data echam_cld_dz=field%dz(:,:,jb) IF (ASSOCIATED(field%dz))
    !$ser data echam_cld_mref=field%mref(:,:,jb) IF (ASSOCIATED(field%mref))
    !$ser data echam_cld_rho=field%rho(:,:,jb) IF (ASSOCIATED(field%rho))
    !$ser data echam_cld_cpair=field%cpair(:,:,jb) IF (ASSOCIATED(field%cpair))
    !$ser data echam_cld_acdnc=field%acdnc(:,:,jb) IF (ASSOCIATED(field%acdnc))
    !$ser data echam_cld_ta=field%ta(:,:,jb) IF (ASSOCIATED(field%ta))
    !$ser data echam_cld_qtrc=field%qtrc(:,:,jb,:) IF (ASSOCIATED(field%qtrc))
    !$ser data echam_cld_aclc=field%aclc(:,:,jb) IF (ASSOCIATED(field%aclc))
    !$ser data echam_cld_aclcov=field%aclcov(:,jb) IF (ASSOCIATED(field%aclcov))
    !$ser data echam_cld_rsfl=field%rsfl(:,jb) IF (ASSOCIATED(field%rsfl))
    !$ser data echam_cld_ssfl=field%ssfl(:,jb) IF (ASSOCIATED(field%ssfl))
    !$ser data echam_cld_hur=field%hur IF (ASSOCIATED(field%hur))
    !$ser data echam_cld_q_cld=field%q_cld IF (ASSOCIATED(field%q_cld))
    !$ser data echam_cld_q_cld_vi=field%q_cld_vi IF (ASSOCIATED(field%q_cld_vi))  
    !$ser data echam_cld_qconv=field%qconv IF (ASSOCIATED(field%qconv))
    !$ser data echam_cld_q_phy=field%q_phy IF (ASSOCIATED(field%q_phy))
    !$ser data echam_cld_q_phy_vi=field%q_phy_vi IF (ASSOCIATED(field%q_phy_vi))
    !$ser data echam_cld_qtrc_cld=tend%qtrc_cld(:,:,jb,:) IF (ASSOCIATED(tend%qtrc_cld))
    !$ser data echam_cld_ta_cld=tend%ta_cld(:,:,jb) IF (ASSOCIATED(tend%ta_cld))
    !$ser data echam_cld_ta_phy=tend%ta_phy(:,:,jb) IF (ASSOCIATED(tend%ta_phy))
    !$ser data echam_cld_qtrc_phy=tend%qtrc_phy(:,:,jb,:) IF (ASSOCIATED(tend%qtrc_phy))
    !$ser verbatim lactive = .FALSE.
    !$ser verbatim IF (lonlyonce) THEN
    !$ser verbatim   lenabled = .FALSE.
    !$ser verbatim END IF
    !$ser verbatim ENDIF

  END SUBROUTINE serialize_cld_output

END MODULE mo_ser_echam_cld
