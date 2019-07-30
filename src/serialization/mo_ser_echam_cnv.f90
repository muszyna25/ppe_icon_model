!--------------------------------------------------------------------
!
! Serialization routine for ECHAM CNV
!
!--------------------------------------------------------------------

MODULE mo_ser_echam_cnv

  USE mo_kind,               ONLY: vp, wp
  USE mo_exception,          ONLY: warning
  USE mo_ser_common,         ONLY: init
  USE mtime,                 ONLY: datetimeToString, MAX_DATETIME_STR_LEN
  USE mo_time_config,        ONLY: time_config
  USE mo_echam_phy_memory,   ONLY: t_echam_phy_field, t_echam_phy_tend
  IMPLICIT NONE

  INTEGER, PARAMETER :: selected_block = 1

  PUBLIC :: serialize_cnv_input
  PUBLIC :: serialize_cnv_output

  CONTAINS

  SUBROUTINE serialize_cnv_input(jg, jb, jcs, jce, nproma, nlev, field, tend)
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
    !$ser verbatim   CALL warning('SER:mo_ser_echam_cnv:input','Serialization is active!')
#if defined( _OPENACC )
    !$ser verbatim   CALL warning('GPU:mo_ser_echam_cnv:input','GPU HOST synchronization forced by serialization!')
    !$ser verbatim   !$ACC UPDATE HOST( field%ta ) IF( ASSOCIATED(field%ta) )
    !$ser verbatim   !$ACC UPDATE HOST( field%qtrc ) IF( ASSOCIATED(field%qtrc) )
    !$ser verbatim   !$ACC UPDATE HOST( field%ua ) IF( ASSOCIATED(field%ua) )
    !$ser verbatim   !$ACC UPDATE HOST( field%va ) IF( ASSOCIATED(field%va) )
    !$ser verbatim   !$ACC UPDATE HOST( field%sftlf ) IF( ASSOCIATED(field%sftlf) )
    !$ser verbatim   !$ACC UPDATE HOST( field%zf ) IF( ASSOCIATED(field%zf) )
    !$ser verbatim   !$ACC UPDATE HOST( field%zh ) IF( ASSOCIATED(field%zh) )
    !$ser verbatim   !$ACC UPDATE HOST( field%mref ) IF( ASSOCIATED(field%mref) )
    !$ser verbatim   !$ACC UPDATE HOST( field%omega ) IF( ASSOCIATED(field%omega) )
    !$ser verbatim   !$ACC UPDATE HOST( field%evap ) IF( ASSOCIATED(field%evap) )
    !$ser verbatim   !$ACC UPDATE HOST( field%presm_new ) IF( ASSOCIATED(field%presm_new) )
    !$ser verbatim   !$ACC UPDATE HOST( field%presi_new ) IF( ASSOCIATED(field%presi_new) )
    !$ser verbatim   !$ACC UPDATE HOST( field%geom ) IF( ASSOCIATED(field%geom) )
    !$ser verbatim   !$ACC UPDATE HOST( field%geoi ) IF( ASSOCIATED(field%geoi) )
    !$ser verbatim   !$ACC UPDATE HOST( field%thvsig ) IF( ASSOCIATED(field%thvsig) )
    !$ser verbatim   !$ACC UPDATE HOST( field%ictop ) IF( ASSOCIATED(field%ictop) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rsfc ) IF( ASSOCIATED(field%rsfc) )
    !$ser verbatim   !$ACC UPDATE HOST( field%ssfc ) IF( ASSOCIATED(field%ssfc) )
    !$ser verbatim   !$ACC UPDATE HOST( field%con_dtrl ) IF( ASSOCIATED(field%con_dtrl) )
    !$ser verbatim   !$ACC UPDATE HOST( field%con_dtri ) IF( ASSOCIATED(field%con_dtri) )
    !$ser verbatim   !$ACC UPDATE HOST( field%con_iteqv ) IF( ASSOCIATED(field%con_iteqv) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rtype ) IF( ASSOCIATED(field%rtype) )
    !$ser verbatim   !$ACC UPDATE HOST( field%topmax ) IF( ASSOCIATED(field%topmax) )
    !$ser verbatim   !$ACC UPDATE HOST( field%q_cnv ) IF( ASSOCIATED(field%q_cnv) )
    !$ser verbatim   !$ACC UPDATE HOST( field%q_cnv_vi ) IF( ASSOCIATED(field%q_cnv_vi) )
    !$ser verbatim   !$ACC UPDATE HOST( field%qconv ) IF( ASSOCIATED(field%qconv) )
    !$ser verbatim   !$ACC UPDATE HOST( field%q_phy ) IF( ASSOCIATED(field%q_phy) )
    !$ser verbatim   !$ACC UPDATE HOST( field%q_phy_vi ) IF( ASSOCIATED(field%q_phy_vi) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%qtrc_dyn ) IF( ASSOCIATED(tend%qtrc_dyn) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%qtrc_phy ) IF( ASSOCIATED(tend%qtrc_phy) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%ua_cnv ) IF( ASSOCIATED(tend%ua_cnv) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%va_cnv ) IF( ASSOCIATED(tend%va_cnv) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%qtrc_cnv ) IF( ASSOCIATED(tend%qtrc_cnv) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%ta_cnv ) IF( ASSOCIATED(tend%ta_cnv) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%ta_phy ) IF( ASSOCIATED(tend%ta_phy) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%ua_phy ) IF( ASSOCIATED(tend%ua_phy) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%va_phy ) IF( ASSOCIATED(tend%va_phy) )
#endif
    !$ser verbatim   CALL datetimeToString(time_config%tc_current_date, date)
    !$ser verbatim   CALL init('icon')
    !$ser savepoint echam_cnv-input jg=jg jb=jb jcs=jcs jce=jce nproma=nproma nlev=nlev date=TRIM(date)
#if defined( SERIALIZE_CREATE_REFERENCE )
    !$ser mode write
#elif defined( SERIALIZE_PERTURB_REFERENCE )
    !$ser mode read-perturb
#elif defined( SERIALIZE_READ_REFERENCE )
    !$ser mode read
#else
#error SERIALIZATION MODE IS NOT SET
#endif
    !$ser data echam_cnv_ta=field%ta(:,:,jb) IF (ASSOCIATED(field%ta))
    !$ser data echam_cnv_qtrc=field%qtrc(:,:,jb,:) IF (ASSOCIATED(field%qtrc))
    !$ser data echam_cnv_ua=field%ua(:,:,jb) IF (ASSOCIATED(field%ua))
    !$ser data echam_cnv_va=field%va(:,:,jb) IF (ASSOCIATED(field%va))
    !$ser data echam_cnv_sftlf=field%sftlf(:,jb) IF (ASSOCIATED(field%sftlf))
    !$ser data echam_cnv_zf=field%zf(:,:,jb) IF (ASSOCIATED(field%zf))
    !$ser data echam_cnv_zh=field%zh(:,:,jb) IF (ASSOCIATED(field%zh))
    !$ser data echam_cnv_mref=field%mref(:,:,jb) IF (ASSOCIATED(field%mref))
    !$ser data echam_cnv_omega=field%omega(:,:,jb) IF (ASSOCIATED(field%omega))
    !$ser data echam_cnv_evap=field%evap(:,jb) IF (ASSOCIATED(field%evap))
    !$ser data echam_cnv_presm_new=field%presm_new(:,:,jb) IF (ASSOCIATED(field%presm_new))
    !$ser data echam_cnv_presi_new=field%presi_new(:,:,jb) IF (ASSOCIATED(field%presi_new))
    !$ser data echam_cnv_geom=field%geom(:,:,jb) IF (ASSOCIATED(field%geom))
    !$ser data echam_cnv_geoi=field%geoi(:,:,jb) IF (ASSOCIATED(field%geoi))
    !$ser data echam_cnv_thvsig=field%thvsig(:,jb) IF (ASSOCIATED(field%thvsig))
    !$ser data echam_cnv_ictop=field%ictop(:,jb) IF (ASSOCIATED(field%ictop))
    !$ser data echam_cnv_rsfc=field%rsfc(:,jb) IF (ASSOCIATED(field%rsfc))
    !$ser data echam_cnv_ssfc=field%ssfc(:,jb) IF (ASSOCIATED(field%ssfc))
    !$ser data echam_cnv_con_dtrl=field%con_dtrl(:,jb) IF (ASSOCIATED(field%con_dtrl))
    !$ser data echam_cnv_con_dtri=field%con_dtri(:,jb) IF (ASSOCIATED(field%con_dtri))
    !$ser data echam_cnv_con_iteqv=field%con_iteqv(:,jb) IF (ASSOCIATED(field%con_iteqv))
    !$ser data echam_cnv_rtype=field%rtype(:,jb) IF (ASSOCIATED(field%rtype))
    !$ser data echam_cnv_topmax=field%topmax(:,jb) IF (ASSOCIATED(field%topmax))
    !$ser data echam_cnv_q_cnv=field%q_cnv(:,:,jb) IF (ASSOCIATED(field%q_cnv))
    !$ser data echam_cnv_qcnv_vi=field%q_cnv_vi(:,jb) IF (ASSOCIATED(field%q_cnv_vi))
    !$ser data echam_cnv_qconv=field%qconv(:,:,jb) IF (ASSOCIATED(field%qconv))
    !$ser data echam_cnv_q_phy=field%q_phy(:,:,jb) IF (ASSOCIATED(field%q_phy))
    !$ser data echam_cnv_q_phy_vi=field%q_phy_vi(:,jb) IF (ASSOCIATED(field%q_phy_vi))
    !$ser data echam_cnv_qtrc_dyn=tend%qtrc_dyn(:,:,jb,:) IF (ASSOCIATED(tend%qtrc_dyn))
    !$ser data echam_cnv_qtrc_phy=tend%qtrc_phy(:,:,jb,:) IF (ASSOCIATED(tend%qtrc_phy))
    !$ser data echam_cnv_ua_cnv=tend%ua_cnv(:,:,jb) IF (ASSOCIATED(tend%ua_cnv))
    !$ser data echam_cnv_va_cnv=tend%va_cnv(:,:,jb) IF (ASSOCIATED(tend%va_cnv))
    !$ser data echam_cnv_qtrc_cnv=tend%qtrc_cnv(:,:,jb,:) IF (ASSOCIATED(tend%qtrc_cnv))
    !$ser data echam_cnv_ta_cnv=tend%ta_cnv(:,:,jb) IF (ASSOCIATED(tend%ta_cnv))
    !$ser data echam_cnv_ta_phy=tend%ta_phy(:,:,jb) IF (ASSOCIATED(tend%ta_phy))
    !$ser data echam_cnv_ua_phy=tend%ua_phy(:,:,jb) IF (ASSOCIATED(tend%ua_phy))
    !$ser data echam_cnv_va_phy=tend%va_phy(:,:,jb) IF (ASSOCIATED(tend%va_phy))
    !$ser verbatim lactive = .FALSE.
    !$ser verbatim IF (lonlyonce) THEN
    !$ser verbatim   lenabled = .FALSE.
    !$ser verbatim END IF
#if defined( SERIALIZE_READ_REFERENCE ) || defined( SERIALIZE_PERTURB_REFERENCE )
#if defined( _OPENACC )
    !$ser verbatim CALL warning('GPU:mo_ser_echam_cnv:input','GPU DEVICE synchronization forced by serialization!')
    !$ser verbatim   !$ACC UPDATE DEVICE( field%ta ) IF( ASSOCIATED(field%ta) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%qtrc ) IF( ASSOCIATED(field%qtrc) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%ua ) IF( ASSOCIATED(field%ua) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%va ) IF( ASSOCIATED(field%va) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%sftlf ) IF( ASSOCIATED(field%sftlf) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%zf ) IF( ASSOCIATED(field%zf) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%zh ) IF( ASSOCIATED(field%zh) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%mref ) IF( ASSOCIATED(field%mref) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%omega ) IF( ASSOCIATED(field%omega) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%evap ) IF( ASSOCIATED(field%evap) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%presm_new ) IF( ASSOCIATED(field%presm_new) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%presi_new ) IF( ASSOCIATED(field%presi_new) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%geom ) IF( ASSOCIATED(field%geom) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%geoi ) IF( ASSOCIATED(field%geoi) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%thvsig ) IF( ASSOCIATED(field%thvsig) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%ictop ) IF( ASSOCIATED(field%ictop) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%rsfc ) IF( ASSOCIATED(field%rsfc) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%ssfc ) IF( ASSOCIATED(field%ssfc) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%con_dtrl ) IF( ASSOCIATED(field%con_dtrl) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%con_dtri ) IF( ASSOCIATED(field%con_dtri) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%con_iteqv ) IF( ASSOCIATED(field%con_iteqv) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%rtype ) IF( ASSOCIATED(field%rtype) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%topmax ) IF( ASSOCIATED(field%topmax) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%q_cnv ) IF( ASSOCIATED(field%q_cnv) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%q_cnv_vi ) IF( ASSOCIATED(field%q_cnv_vi) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%qconv ) IF( ASSOCIATED(field%qconv) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%q_phy ) IF( ASSOCIATED(field%q_phy) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%q_phy_vi ) IF( ASSOCIATED(field%q_phy_vi) )
    !$ser verbatim   !$ACC UPDATE DEVICE( tend%qtrc_dyn ) IF( ASSOCIATED(tend%qtrc_dyn) )
    !$ser verbatim   !$ACC UPDATE DEVICE( tend%qtrc_phy ) IF( ASSOCIATED(tend%qtrc_phy) )
    !$ser verbatim   !$ACC UPDATE DEVICE( tend%ua_cnv ) IF( ASSOCIATED(tend%ua_cnv) )
    !$ser verbatim   !$ACC UPDATE DEVICE( tend%va_cnv ) IF( ASSOCIATED(tend%va_cnv) )
    !$ser verbatim   !$ACC UPDATE DEVICE( tend%qtrc_cnv ) IF( ASSOCIATED(tend%qtrc_cnv) )
    !$ser verbatim   !$ACC UPDATE DEVICE( tend%ta_cnv ) IF( ASSOCIATED(tend%ta_cnv) )
    !$ser verbatim   !$ACC UPDATE DEVICE( tend%ta_phy ) IF( ASSOCIATED(tend%ta_phy) )
    !$ser verbatim   !$ACC UPDATE DEVICE( tend%ua_phy ) IF( ASSOCIATED(tend%ua_phy) )
    !$ser verbatim   !$ACC UPDATE DEVICE( tend%va_phy ) IF( ASSOCIATED(tend%va_phy) )
#endif
#endif
    !$ser verbatim ENDIF

  END SUBROUTINE serialize_cnv_input

  SUBROUTINE serialize_cnv_output(jg, jb, jcs, jce, nproma, nlev, field, tend)
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
    !$ser verbatim   CALL warning('SER:mo_ser_echam_cnv:output','Serialization is active!')
#if defined( _OPENACC )
    !$ser verbatim   CALL warning('GPU:mo_ser_echam_cnv:output','GPU HOST synchronization forced by serialization!')
    !$ser verbatim   !$ACC UPDATE HOST(field%ta ) IF( ASSOCIATED(field%ta) )
    !$ser verbatim   !$ACC UPDATE HOST(field%qtrc ) IF( ASSOCIATED(field%qtrc) )
    !$ser verbatim   !$ACC UPDATE HOST(field%ua ) IF( ASSOCIATED(field%ua) )
    !$ser verbatim   !$ACC UPDATE HOST(field%va ) IF( ASSOCIATED(field%va) )
    !$ser verbatim   !$ACC UPDATE HOST(field%sftlf ) IF( ASSOCIATED(field%sftlf) )
    !$ser verbatim   !$ACC UPDATE HOST(field%zf ) IF( ASSOCIATED(field%zf) )
    !$ser verbatim   !$ACC UPDATE HOST(field%zh ) IF( ASSOCIATED(field%zh) )
    !$ser verbatim   !$ACC UPDATE HOST(field%mref ) IF( ASSOCIATED(field%mref) )
    !$ser verbatim   !$ACC UPDATE HOST(field%omega ) IF( ASSOCIATED(field%omega) )
    !$ser verbatim   !$ACC UPDATE HOST(field%evap ) IF( ASSOCIATED(field%evap) )
    !$ser verbatim   !$ACC UPDATE HOST(field%presm_new ) IF( ASSOCIATED(field%presm_new) )
    !$ser verbatim   !$ACC UPDATE HOST(field%presi_new ) IF( ASSOCIATED(field%presi_new) )
    !$ser verbatim   !$ACC UPDATE HOST(field%geom ) IF( ASSOCIATED(field%geom) )
    !$ser verbatim   !$ACC UPDATE HOST(field%geoi ) IF( ASSOCIATED(field%geoi) )
    !$ser verbatim   !$ACC UPDATE HOST(field%thvsig ) IF( ASSOCIATED(field%thvsig) )
    !$ser verbatim   !$ACC UPDATE HOST(field%ictop ) IF( ASSOCIATED(field%ictop) )
    !$ser verbatim   !$ACC UPDATE HOST(field%rsfc ) IF( ASSOCIATED(field%rsfc) )
    !$ser verbatim   !$ACC UPDATE HOST(field%ssfc ) IF( ASSOCIATED(field%ssfc) )
    !$ser verbatim   !$ACC UPDATE HOST(field%con_dtrl ) IF( ASSOCIATED(field%con_dtrl) )
    !$ser verbatim   !$ACC UPDATE HOST(field%con_dtri ) IF( ASSOCIATED(field%con_dtri) )
    !$ser verbatim   !$ACC UPDATE HOST(field%con_iteqv ) IF( ASSOCIATED(field%con_iteqv) )
    !$ser verbatim   !$ACC UPDATE HOST(field%rtype ) IF( ASSOCIATED(field%rtype) )
    !$ser verbatim   !$ACC UPDATE HOST(field%topmax ) IF( ASSOCIATED(field%topmax) )
    !$ser verbatim   !$ACC UPDATE HOST(field%q_cnv ) IF( ASSOCIATED(field%q_cnv) )
    !$ser verbatim   !$ACC UPDATE HOST(field%q_cnv_vi ) IF( ASSOCIATED(field%q_cnv_vi) )
    !$ser verbatim   !$ACC UPDATE HOST(field%qconv ) IF( ASSOCIATED(field%qconv) )
    !$ser verbatim   !$ACC UPDATE HOST(field%q_phy ) IF( ASSOCIATED(field%q_phy) )
    !$ser verbatim   !$ACC UPDATE HOST(field%q_phy_vi ) IF( ASSOCIATED(field%q_phy_vi) )
    !$ser verbatim   !$ACC UPDATE HOST(tend%qtrc_dyn ) IF( ASSOCIATED(tend%qtrc_dyn) )
    !$ser verbatim   !$ACC UPDATE HOST(tend%qtrc_phy ) IF( ASSOCIATED(tend%qtrc_phy) )
    !$ser verbatim   !$ACC UPDATE HOST(tend%ua_cnv ) IF( ASSOCIATED(tend%ua_cnv) )
    !$ser verbatim   !$ACC UPDATE HOST(tend%va_cnv ) IF( ASSOCIATED(tend%va_cnv) )
    !$ser verbatim   !$ACC UPDATE HOST(tend%qtrc_cnv ) IF( ASSOCIATED(tend%qtrc_cnv) )
    !$ser verbatim   !$ACC UPDATE HOST(tend%ta_cnv ) IF( ASSOCIATED(tend%ta_cnv) )
    !$ser verbatim   !$ACC UPDATE HOST(tend%ta_phy ) IF( ASSOCIATED(tend%ta_phy) )
    !$ser verbatim   !$ACC UPDATE HOST(tend%ua_phy ) IF( ASSOCIATED(tend%ua_phy) )
    !$ser verbatim   !$ACC UPDATE HOST(tend%va_phy ) IF( ASSOCIATED(tend%va_phy) )
#endif
    !$ser verbatim   CALL datetimeToString(time_config%tc_current_date, date)
    !$ser verbatim   CALL init('icon')
    !$ser savepoint echam_cnv-output jg=jg jb=jb jcs=jcs jce=jce nproma=nproma nlev=nlev date=TRIM(date) 
    !$ser mode write
    !$ser data echam_cnv_ta=field%ta(:,:,jb) IF (ASSOCIATED(field%ta))
    !$ser data echam_cnv_qtrc=field%qtrc(:,:,jb,:) IF (ASSOCIATED(field%qtrc))
    !$ser data echam_cnv_ua=field%ua(:,:,jb) IF (ASSOCIATED(field%ua))
    !$ser data echam_cnv_va=field%va(:,:,jb) IF (ASSOCIATED(field%va))
    !$ser data echam_cnv_sftlf=field%sftlf(:,jb) IF (ASSOCIATED(field%sftlf))
    !$ser data echam_cnv_zf=field%zf(:,:,jb) IF (ASSOCIATED(field%zf))
    !$ser data echam_cnv_zh=field%zh(:,:,jb) IF (ASSOCIATED(field%zh))
    !$ser data echam_cnv_mref=field%mref(:,:,jb) IF (ASSOCIATED(field%mref))
    !$ser data echam_cnv_omega=field%omega(:,:,jb) IF (ASSOCIATED(field%omega))
    !$ser data echam_cnv_evap=field%evap(:,jb) IF (ASSOCIATED(field%evap))
    !$ser data echam_cnv_presm_new=field%presm_new(:,:,jb) IF (ASSOCIATED(field%presm_new))
    !$ser data echam_cnv_presi_new=field%presi_new(:,:,jb) IF (ASSOCIATED(field%presi_new))
    !$ser data echam_cnv_geom=field%geom(:,:,jb) IF (ASSOCIATED(field%geom))
    !$ser data echam_cnv_geoi=field%geoi(:,:,jb) IF (ASSOCIATED(field%geoi))
    !$ser data echam_cnv_thvsig=field%thvsig(:,jb) IF (ASSOCIATED(field%thvsig))
    !$ser data echam_cnv_ictop=field%ictop(:,jb) IF (ASSOCIATED(field%ictop))
    !$ser data echam_cnv_rsfc=field%rsfc(:,jb) IF (ASSOCIATED(field%rsfc))
    !$ser data echam_cnv_ssfc=field%ssfc(:,jb) IF (ASSOCIATED(field%ssfc))
    !$ser data echam_cnv_con_dtrl=field%con_dtrl(:,jb) IF (ASSOCIATED(field%con_dtrl))
    !$ser data echam_cnv_con_dtri=field%con_dtri(:,jb) IF (ASSOCIATED(field%con_dtri))
    !$ser data echam_cnv_con_iteqv=field%con_iteqv(:,jb) IF (ASSOCIATED(field%con_iteqv))
    !$ser data echam_cnv_rtype=field%rtype(:,jb) IF (ASSOCIATED(field%rtype))
    !$ser data echam_cnv_topmax=field%topmax(:,jb) IF (ASSOCIATED(field%topmax))
    !$ser data echam_cnv_q_cnv=field%q_cnv(:,:,jb) IF (ASSOCIATED(field%q_cnv))
    !$ser data echam_cnv_qcnv_vi=field%q_cnv_vi(:,jb) IF (ASSOCIATED(field%q_cnv_vi))
    !$ser data echam_cnv_qconv=field%qconv(:,:,jb) IF (ASSOCIATED(field%qconv))
    !$ser data echam_cnv_q_phy=field%q_phy(:,:,jb) IF (ASSOCIATED(field%q_phy))
    !$ser data echam_cnv_q_phy_vi=field%q_phy_vi(:,jb) IF (ASSOCIATED(field%q_phy_vi))
    !$ser data echam_cnv_qtrc_dyn=tend%qtrc_dyn(:,:,jb,:) IF (ASSOCIATED(tend%qtrc_dyn))
    !$ser data echam_cnv_qtrc_phy=tend%qtrc_phy(:,:,jb,:) IF (ASSOCIATED(tend%qtrc_phy))
    !$ser data echam_cnv_ua_cnv=tend%ua_cnv(:,:,jb) IF (ASSOCIATED(tend%ua_cnv))
    !$ser data echam_cnv_va_cnv=tend%va_cnv(:,:,jb) IF (ASSOCIATED(tend%va_cnv))
    !$ser data echam_cnv_qtrc_cnv=tend%qtrc_cnv(:,:,jb,:) IF (ASSOCIATED(tend%qtrc_cnv))
    !$ser data echam_cnv_ta_cnv=tend%ta_cnv(:,:,jb) IF (ASSOCIATED(tend%ta_cnv))
    !$ser data echam_cnv_ta_phy=tend%ta_phy(:,:,jb) IF (ASSOCIATED(tend%ta_phy))
    !$ser data echam_cnv_ua_phy=tend%ua_phy(:,:,jb) IF (ASSOCIATED(tend%ua_phy))
    !$ser data echam_cnv_va_phy=tend%va_phy(:,:,jb) IF (ASSOCIATED(tend%va_phy))
    !$ser verbatim lactive = .FALSE.
    !$ser verbatim IF (lonlyonce) THEN
    !$ser verbatim   lenabled = .FALSE.
    !$ser verbatim END IF
    !$ser verbatim ENDIF

  END SUBROUTINE serialize_cnv_output

END MODULE mo_ser_echam_cnv
