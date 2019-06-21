!--------------------------------------------------------------------
!
! Serialization routine for ECHAM RHT
!
!--------------------------------------------------------------------

MODULE mo_ser_echam_rht

  USE mo_kind,               ONLY: vp, wp
  USE mo_exception,          ONLY: warning
  USE mo_ser_common,         ONLY: init
  USE mtime,                 ONLY: datetimeToString, MAX_DATETIME_STR_LEN
  USE mo_time_config,        ONLY: time_config
  USE mo_echam_phy_memory,   ONLY: t_echam_phy_field, t_echam_phy_tend
  IMPLICIT NONE

  INTEGER, PARAMETER :: selected_block = 1

  PUBLIC :: serialize_rht_input
  PUBLIC :: serialize_rht_output

  CONTAINS

  SUBROUTINE serialize_rht_input(jg, jb, jcs, jce, nproma, nlev, field, tend, emis_rad)
    INTEGER, INTENT(IN)    :: jg, jb, jcs, jce, nproma, nlev
    TYPE(t_echam_phy_field), POINTER, INTENT(INOUT) :: field
    TYPE(t_echam_phy_tend), POINTER, INTENT(INOUT)  :: tend
    REAL(kind=wp), INTENT(INOUT)                    :: emis_rad(:,:)

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: date
    LOGICAL, PARAMETER :: lonlyonce = .TRUE.
    LOGICAL, SAVE :: lenabled = .TRUE.
    LOGICAL, SAVE :: lactive = .FALSE.

    !$ser verbatim IF (selected_block < 0 .OR. jb == selected_block) THEN
    !$ser verbatim   lactive = .TRUE.
    !$ser verbatim ENDIF
    !$ser verbatim IF (lenabled .and. lactive) THEN
    !$ser verbatim   CALL warning('SER:mo_ser_echam_rht:input','Serialization is active!')
#if defined( _OPENACC )
    !$ser verbatim   CALL warning('GPU:mo_ser_echam_rht:input','GPU HOST synchronization forced by serialization!')
    !$ser verbatim   !$ACC UPDATE HOST( field%cosmu0 ) IF( ASSOCIATED(field%cosmu0) )
    !$ser verbatim   !$ACC UPDATE HOST( field%daylght_frc ) IF( ASSOCIATED(field%daylght_frc) )
    !$ser verbatim   !$ACC UPDATE HOST( field%ts_rad ) IF( ASSOCIATED(field%ts_rad) )
    !$ser verbatim   !$ACC UPDATE HOST( field%ts_rad_rt ) IF( ASSOCIATED(field%ts_rad_rt) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rsd_rt ) IF( ASSOCIATED(field%rsd_rt) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rsu_rt ) IF( ASSOCIATED(field%rsu_rt) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rsdcs_rt ) IF( ASSOCIATED(field%rsdcs_rt) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rsucs_rt ) IF( ASSOCIATED(field%rsucs_rt) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rld_rt ) IF( ASSOCIATED(field%rld_rt) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rlu_rt ) IF( ASSOCIATED(field%rlu_rt) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rldcs_rt ) IF( ASSOCIATED(field%rldcs_rt) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rlucs_rt ) IF( ASSOCIATED(field%rlucs_rt) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rvds_dir_rt ) IF( ASSOCIATED(field%rvds_dir_rt) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rpds_dir_rt ) IF( ASSOCIATED(field%rpds_dir_rt) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rnds_dir_rt ) IF( ASSOCIATED(field%rnds_dir_rt) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rvds_dif_rt ) IF( ASSOCIATED(field%rvds_dif_rt) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rpds_dif_rt ) IF( ASSOCIATED(field%rpds_dif_rt) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rnds_dif_rt ) IF( ASSOCIATED(field%rnds_dif_rt) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rvus_rt ) IF( ASSOCIATED(field%rvus_rt) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rpus_rt ) IF( ASSOCIATED(field%rpus_rt) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rnus_rt ) IF( ASSOCIATED(field%rnus_rt) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rsdt ) IF( ASSOCIATED(field%rsdt) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rsut ) IF( ASSOCIATED(field%rsut) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rsds ) IF( ASSOCIATED(field%rsds) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rsus ) IF( ASSOCIATED(field%rsus) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rsutcs ) IF( ASSOCIATED(field%rsutcs) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rsdscs ) IF( ASSOCIATED(field%rsdscs) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rsuscs ) IF( ASSOCIATED(field%rsuscs) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rvds_dir ) IF( ASSOCIATED(field%rvds_dir) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rpds_dir ) IF( ASSOCIATED(field%rpds_dir) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rnds_dir ) IF( ASSOCIATED(field%rnds_dir) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rvds_dif ) IF( ASSOCIATED(field%rvds_dif) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rpds_dif ) IF( ASSOCIATED(field%rpds_dif) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rnds_dif ) IF( ASSOCIATED(field%rnds_dif) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rvus ) IF( ASSOCIATED(field%rvus) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rpus ) IF( ASSOCIATED(field%rpus) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rnus ) IF( ASSOCIATED(field%rnus) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rlut ) IF( ASSOCIATED(field%rlut) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rlds ) IF( ASSOCIATED(field%rlds) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rlus ) IF( ASSOCIATED(field%rlus) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rlutcs ) IF( ASSOCIATED(field%rlutcs) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rldscs ) IF( ASSOCIATED(field%rldscs) )
    !$ser verbatim   !$ACC UPDATE HOST( field%q_rlw_nlev ) IF( ASSOCIATED(field%q_rlw_nlev) )
    !$ser verbatim   !$ACC UPDATE HOST( field%qconv ) IF( ASSOCIATED(field%qconv) )
    !$ser verbatim   !$ACC UPDATE HOST( field%q_rad ) IF( ASSOCIATED(field%q_rad) )
    !$ser verbatim   !$ACC UPDATE HOST( field%q_rad_vi ) IF( ASSOCIATED(field%q_rad_vi) )
    !$ser verbatim   !$ACC UPDATE HOST( field%q_rsw ) IF( ASSOCIATED(field%q_rsw) )
    !$ser verbatim   !$ACC UPDATE HOST( field%q_rsw_vi ) IF( ASSOCIATED(field%q_rsw_vi) )
    !$ser verbatim   !$ACC UPDATE HOST( field%q_rlw ) IF( ASSOCIATED(field%q_rlw) )
    !$ser verbatim   !$ACC UPDATE HOST( field%q_rlw_vi ) IF( ASSOCIATED(field%q_rlw_vi) )
    !$ser verbatim   !$ACC UPDATE HOST( field%q_phy ) IF( ASSOCIATED(field%q_phy) )
    !$ser verbatim   !$ACC UPDATE HOST( field%q_phy_vi ) IF( ASSOCIATED(field%q_phy_vi) )
    !$ser verbatim   !$ACC UPDATE HOST( field%q_rlw_nlev ) IF( ASSOCIATED(field%q_rlw_nlev) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%ta_rad ) IF( ASSOCIATED(tend%ta_rad) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%ta_rsw ) IF( ASSOCIATED(tend%ta_rsw) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%ta_rlw ) IF( ASSOCIATED(tend%ta_rlw) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%ta_phy ) IF( ASSOCIATED(tend%ta_phy) )
    !$ser verbatim   !NOACC UPDATE HOST( emis_rad ) ! NOT ON GPU
#endif
    !$ser verbatim   CALL datetimeToString(time_config%tc_current_date, date)
    !$ser verbatim   CALL init('icon')
    !$ser savepoint echam_rht-input jg=jg jb=jb jcs=jcs jce=jce nproma=nproma nlev=nlev date=TRIM(date)
#if defined( SERIALIZE_CREATE_REFERENCE )
    !$ser mode write
#elif defined( SERIALIZE_PERTURB_REFERENCE )
    !$ser mode read-perturb
#elif defined( SERIALIZE_READ_REFERENCE )
    !$ser mode read
#else
#error SERIALIZATION MODE IS NOT SET
#endif
    !$ser data echam_rht_cosmu0=field%cosmu0(:,jb) IF (ASSOCIATED(field%cosmu0))
    !$ser data echam_rht_daylght_frc=field%daylght_frc(:,jb) IF (ASSOCIATED(field%daylght_frc))
    !$ser data echam_rht_ts_rad=field%ts_rad(:,jb) IF (ASSOCIATED(field%ts_rad))
    !$ser data echam_rht_ts_rad_rt=field%ts_rad_rt(:,jb) IF (ASSOCIATED(field%ts_rad_rt))
    !$ser data echam_rht_rsd_rt=field%rsd_rt(:,:,jb) IF (ASSOCIATED(field%rsd_rt))
    !$ser data echam_rht_rsu_rt=field%rsu_rt(:,:,jb) IF (ASSOCIATED(field%rsu_rt))
    !$ser data echam_rht_rsdcs_rt=field%rsdcs_rt(:,:,jb) IF (ASSOCIATED(field%rsdcs_rt))
    !$ser data echam_rht_rsucs_rt=field%rsucs_rt(:,:,jb) IF (ASSOCIATED(field%rsucs_rt))
    !$ser data echam_rht_rld_rt=field%rld_rt(:,:,jb) IF (ASSOCIATED(field%rld_rt))
    !$ser data echam_rht_rlu_rt=field%rlu_rt(:,:,jb) IF (ASSOCIATED(field%rlu_rt))
    !$ser data echam_rht_rldcs_rt=field%rldcs_rt(:,:,jb) IF (ASSOCIATED(field%rldcs_rt))
    !$ser data echam_rht_rlucs_rt=field%rlucs_rt(:,:,jb) IF (ASSOCIATED(field%rlucs_rt))
    !$ser data echam_rht_rvds_dir_rt=field%rvds_dir_rt(:,jb) IF (ASSOCIATED(field%rvds_dir_rt))
    !$ser data echam_rht_rpds_dir_rt=field%rpds_dir_rt(:,jb) IF (ASSOCIATED(field%rpds_dir_rt))
    !$ser data echam_rht_rnds_dir_rt=field%rnds_dir_rt(:,jb) IF (ASSOCIATED(field%rnds_dir_rt))
    !$ser data echam_rht_rvds_dif_rt=field%rvds_dif_rt(:,jb) IF (ASSOCIATED(field%rvds_dif_rt))
    !$ser data echam_rht_rpds_dif_rt=field%rpds_dif_rt(:,jb) IF (ASSOCIATED(field%rpds_dif_rt))
    !$ser data echam_rht_rnds_dif_rt=field%rnds_dif_rt(:,jb) IF (ASSOCIATED(field%rnds_dif_rt))
    !$ser data echam_rht_rvus_rt=field%rvus_rt(:,jb) IF (ASSOCIATED(field%rvus_rt))
    !$ser data echam_rht_rpus_rt=field%rpus_rt(:,jb) IF (ASSOCIATED(field%rpus_rt))
    !$ser data echam_rht_rnus_rt=field%rnus_rt(:,jb) IF (ASSOCIATED(field%rnus_rt))
    !$ser data echam_rht_rsdt=field%rsdt(:,jb) IF (ASSOCIATED(field%rsdt))
    !$ser data echam_rht_rsut=field%rsut(:,jb) IF (ASSOCIATED(field%rsut))
    !$ser data echam_rht_rsds=field%rsds(:,jb) IF (ASSOCIATED(field%rsds))
    !$ser data echam_rht_rsus=field%rsus(:,jb) IF (ASSOCIATED(field%rsus))
    !$ser data echam_rht_rsutcs=field%rsutcs(:,jb) IF (ASSOCIATED(field%rsutcs))
    !$ser data echam_rht_rsdscs=field%rsdscs(:,jb) IF (ASSOCIATED(field%rsdscs))
    !$ser data echam_rht_rsuscs=field%rsuscs(:,jb) IF (ASSOCIATED(field%rsuscs))
    !$ser data echam_rht_rvds_dir=field%rvds_dir(:,jb) IF (ASSOCIATED(field%rvds_dir))
    !$ser data echam_rht_rpds_dir=field%rpds_dir(:,jb) IF (ASSOCIATED(field%rpds_dir))
    !$ser data echam_rht_rnds_dir=field%rnds_dir(:,jb) IF (ASSOCIATED(field%rnds_dir))
    !$ser data echam_rht_rvds_dif=field%rvds_dif(:,jb) IF (ASSOCIATED(field%rvds_dif))
    !$ser data echam_rht_rpds_dif=field%rpds_dif(:,jb) IF (ASSOCIATED(field%rpds_dif))
    !$ser data echam_rht_rnds_dif=field%rnds_dif(:,jb) IF (ASSOCIATED(field%rnds_dif))
    !$ser data echam_rht_rvus=field%rvus(:,jb) IF (ASSOCIATED(field%rvus))
    !$ser data echam_rht_rpus=field%rpus(:,jb) IF (ASSOCIATED(field%rpus))
    !$ser data echam_rht_rnus=field%rnus(:,jb) IF (ASSOCIATED(field%rnus))
    !$ser data echam_rht_rlut=field%rlut(:,jb) IF (ASSOCIATED(field%rlut))
    !$ser data echam_rht_rlds=field%rlds(:,jb) IF (ASSOCIATED(field%rlds))
    !$ser data echam_rht_rlus=field%rlus(:,jb) IF (ASSOCIATED(field%rlus))
    !$ser data echam_rht_rlutcs=field%rlutcs(:,jb) IF (ASSOCIATED(field%rlutcs))
    !$ser data echam_rht_rldscs=field%rldscs(:,jb) IF (ASSOCIATED(field%rldscs))
    !$ser data echam_rht_q_rlw_nlev=field%q_rlw_nlev(:,jb) IF (ASSOCIATED(field%q_rlw_nlev))
    !$ser data echam_rht_qconv=field%qconv(:,:,jb) IF (ASSOCIATED(field%qconv))
    !$ser data echam_rht_q_rad=field%q_rad(:,:,jb) IF (ASSOCIATED(field%q_rad))
    !$ser data echam_rht_q_rad_vi=field%q_rad_vi(:,jb) IF (ASSOCIATED(field%q_rad_vi))
    !$ser data echam_rht_q_rsw=field%q_rsw(:,:,jb) IF (ASSOCIATED(field%q_rsw))
    !$ser data echam_rht_q_rsw_vi=field%q_rsw_vi(:,jb) IF (ASSOCIATED(field%q_rsw_vi))
    !$ser data echam_rht_q_rlw=field%q_rlw(:,:,jb) IF (ASSOCIATED(field%q_rlw))
    !$ser data echam_rht_q_rlw_vi=field%q_rlw_vi(:,jb) IF (ASSOCIATED(field%q_rlw_vi))
    !$ser data echam_rht_q_phy=field%q_phy(:,:,jb) IF (ASSOCIATED(field%q_phy))
    !$ser data echam_rht_q_phy_vi=field%q_phy_vi(:,jb) IF (ASSOCIATED(field%q_phy_vi))
    !$ser data echam_rht_ta_rad=tend%ta_rad(:,:,jb) IF (ASSOCIATED(tend%ta_rad))
    !$ser data echam_rht_ta_rsw=tend%ta_rsw(:,:,jb) IF (ASSOCIATED(tend%ta_rsw))
    !$ser data echam_rht_ta_rlw=tend%ta_rlw(:,:,jb) IF (ASSOCIATED(tend%ta_rlw))
    !$ser data echam_rht_ta_phy=tend%ta_phy(:,:,jb) IF (ASSOCIATED(tend%ta_phy))
    !$ser data echam_rht_ta=tend%ta(:,:,jb) IF (ASSOCIATED(tend%ta))
    !$ser data echam_rht_emis_rad=emis_rad(:,jb)
    !$ser verbatim lactive = .FALSE.
    !$ser verbatim IF (lonlyonce) THEN
    !$ser verbatim   lenabled = .FALSE.
    !$ser verbatim END IF
#if defined( SERIALIZE_READ_REFERENCE ) || defined( SERIALIZE_PERTURB_REFERENCE )
#if defined( _OPENACC )
    !$ser verbatim CALL warning('GPU:mo_ser_echam_rht:input','GPU DEVICE synchronization forced by serialization!')
    !$ser verbatim   !$ACC UPDATE DEVICE( field%cosmu0 ) IF( ASSOCIATED(field%cosmu0) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%daylght_frc ) IF( ASSOCIATED(field%daylght_frc) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%ts_rad ) IF( ASSOCIATED(field%ts_rad) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%ts_rad_rt ) IF( ASSOCIATED(field%ts_rad_rt) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%rsd_rt ) IF( ASSOCIATED(field%rsd_rt) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%rsu_rt ) IF( ASSOCIATED(field%rsu_rt) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%rsdcs_rt ) IF( ASSOCIATED(field%rsdcs_rt) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%rsucs_rt ) IF( ASSOCIATED(field%rsucs_rt) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%rld_rt ) IF( ASSOCIATED(field%rld_rt) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%rlu_rt ) IF( ASSOCIATED(field%rlu_rt) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%rldcs_rt ) IF( ASSOCIATED(field%rldcs_rt) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%rlucs_rt ) IF( ASSOCIATED(field%rlucs_rt) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%rvds_dir_rt ) IF( ASSOCIATED(field%rvds_dir_rt) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%rpds_dir_rt ) IF( ASSOCIATED(field%rpds_dir_rt) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%rnds_dir_rt ) IF( ASSOCIATED(field%rnds_dir_rt) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%rvds_dif_rt ) IF( ASSOCIATED(field%rvds_dif_rt) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%rpds_dif_rt ) IF( ASSOCIATED(field%rpds_dif_rt) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%rnds_dif_rt ) IF( ASSOCIATED(field%rnds_dif_rt) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%rvus_rt ) IF( ASSOCIATED(field%rvus_rt) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%rpus_rt ) IF( ASSOCIATED(field%rpus_rt) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%rnus_rt ) IF( ASSOCIATED(field%rnus_rt) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%rsdt ) IF( ASSOCIATED(field%rsdt) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%rsut ) IF( ASSOCIATED(field%rsut) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%rsds ) IF( ASSOCIATED(field%rsds) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%rsus ) IF( ASSOCIATED(field%rsus) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%rsutcs ) IF( ASSOCIATED(field%rsutcs) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%rsdscs ) IF( ASSOCIATED(field%rsdscs) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%rsuscs ) IF( ASSOCIATED(field%rsuscs) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%rvds_dir ) IF( ASSOCIATED(field%rvds_dir) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%rpds_dir ) IF( ASSOCIATED(field%rpds_dir) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%rnds_dir ) IF( ASSOCIATED(field%rnds_dir) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%rvds_dif ) IF( ASSOCIATED(field%rvds_dif) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%rpds_dif ) IF( ASSOCIATED(field%rpds_dif) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%rnds_dif ) IF( ASSOCIATED(field%rnds_dif) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%rvus ) IF( ASSOCIATED(field%rvus) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%rpus ) IF( ASSOCIATED(field%rpus) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%rnus ) IF( ASSOCIATED(field%rnus) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%rlut ) IF( ASSOCIATED(field%rlut) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%rlds ) IF( ASSOCIATED(field%rlds) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%rlus ) IF( ASSOCIATED(field%rlus) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%rlutcs ) IF( ASSOCIATED(field%rlutcs) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%rldscs ) IF( ASSOCIATED(field%rldscs) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%q_rlw_nlev ) IF( ASSOCIATED(field%q_rlw_nlev) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%qconv ) IF( ASSOCIATED(field%qconv) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%q_rad ) IF( ASSOCIATED(field%q_rad) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%q_rad_vi ) IF( ASSOCIATED(field%q_rad_vi) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%q_rsw ) IF( ASSOCIATED(field%q_rsw) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%q_rsw_vi ) IF( ASSOCIATED(field%q_rsw_vi) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%q_rlw ) IF( ASSOCIATED(field%q_rlw) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%q_rlw_vi ) IF( ASSOCIATED(field%q_rlw_vi) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%q_phy ) IF( ASSOCIATED(field%q_phy) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%q_phy_vi ) IF( ASSOCIATED(field%q_phy_vi) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%q_rlw_nlev ) IF( ASSOCIATED(field%q_rlw_nlev) )
    !$ser verbatim   !$ACC UPDATE DEVICE( tend%ta_rad ) IF( ASSOCIATED(tend%ta_rad) )
    !$ser verbatim   !$ACC UPDATE DEVICE( tend%ta_rsw ) IF( ASSOCIATED(tend%ta_rsw) )
    !$ser verbatim   !$ACC UPDATE DEVICE( tend%ta_rlw ) IF( ASSOCIATED(tend%ta_rlw) )
    !$ser verbatim   !$ACC UPDATE DEVICE( tend%ta_phy ) IF( ASSOCIATED(tend%ta_phy) )
    !$ser verbatim   !NOACC UPDATE DEVICE( emis_rad ) ! NOT ON GPU
#endif
#endif
    !$ser verbatim ENDIF

  END SUBROUTINE serialize_rht_input

  SUBROUTINE serialize_rht_output(jg, jb, jcs, jce, nproma, nlev, field, tend, emis_rad)
    INTEGER, INTENT(IN)    :: jg, jb, jcs, jce, nproma, nlev
    TYPE(t_echam_phy_field), POINTER, INTENT(INOUT) :: field
    TYPE(t_echam_phy_tend), POINTER, INTENT(INOUT)  :: tend
    REAL(kind=wp), INTENT(INOUT)                    :: emis_rad(:,:)

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: date
    LOGICAL, PARAMETER :: lonlyonce = .TRUE.
    LOGICAL, SAVE :: lenabled = .TRUE.
    LOGICAL, SAVE :: lactive = .FALSE.

    !$ser verbatim IF (selected_block < 0 .OR. jb == selected_block) THEN
    !$ser verbatim   lactive = .TRUE.
    !$ser verbatim ENDIF
    !$ser verbatim IF (lenabled .and. lactive) THEN
    !$ser verbatim   CALL warning('SER:mo_ser_echam_rht:output','Serialization is active!')
#if defined( _OPENACC )
    !$ser verbatim   CALL warning('GPU:mo_ser_echam_rht:output','GPU HOST synchronization forced by serialization!')
    !$ser verbatim   !$ACC UPDATE HOST( field%cosmu0 ) IF( ASSOCIATED(field%cosmu0) )
    !$ser verbatim   !$ACC UPDATE HOST( field%daylght_frc ) IF( ASSOCIATED(field%daylght_frc) )
    !$ser verbatim   !$ACC UPDATE HOST( field%ts_rad ) IF( ASSOCIATED(field%ts_rad) )
    !$ser verbatim   !$ACC UPDATE HOST( field%ts_rad_rt ) IF( ASSOCIATED(field%ts_rad_rt) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rsd_rt ) IF( ASSOCIATED(field%rsd_rt) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rsu_rt ) IF( ASSOCIATED(field%rsu_rt) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rsdcs_rt ) IF( ASSOCIATED(field%rsdcs_rt) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rsucs_rt ) IF( ASSOCIATED(field%rsucs_rt) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rld_rt ) IF( ASSOCIATED(field%rld_rt) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rlu_rt ) IF( ASSOCIATED(field%rlu_rt) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rldcs_rt ) IF( ASSOCIATED(field%rldcs_rt) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rlucs_rt ) IF( ASSOCIATED(field%rlucs_rt) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rvds_dir_rt ) IF( ASSOCIATED(field%rvds_dir_rt) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rpds_dir_rt ) IF( ASSOCIATED(field%rpds_dir_rt) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rnds_dir_rt ) IF( ASSOCIATED(field%rnds_dir_rt) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rvds_dif_rt ) IF( ASSOCIATED(field%rvds_dif_rt) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rpds_dif_rt ) IF( ASSOCIATED(field%rpds_dif_rt) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rnds_dif_rt ) IF( ASSOCIATED(field%rnds_dif_rt) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rvus_rt ) IF( ASSOCIATED(field%rvus_rt) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rpus_rt ) IF( ASSOCIATED(field%rpus_rt) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rnus_rt ) IF( ASSOCIATED(field%rnus_rt) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rsdt ) IF( ASSOCIATED(field%rsdt) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rsut ) IF( ASSOCIATED(field%rsut) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rsds ) IF( ASSOCIATED(field%rsds) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rsus ) IF( ASSOCIATED(field%rsus) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rsutcs ) IF( ASSOCIATED(field%rsutcs) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rsdscs ) IF( ASSOCIATED(field%rsdscs) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rsuscs ) IF( ASSOCIATED(field%rsuscs) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rvds_dir ) IF( ASSOCIATED(field%rvds_dir) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rpds_dir ) IF( ASSOCIATED(field%rpds_dir) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rnds_dir ) IF( ASSOCIATED(field%rnds_dir) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rvds_dif ) IF( ASSOCIATED(field%rvds_dif) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rpds_dif ) IF( ASSOCIATED(field%rpds_dif) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rnds_dif ) IF( ASSOCIATED(field%rnds_dif) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rvus ) IF( ASSOCIATED(field%rvus) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rpus ) IF( ASSOCIATED(field%rpus) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rnus ) IF( ASSOCIATED(field%rnus) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rlut ) IF( ASSOCIATED(field%rlut) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rlds ) IF( ASSOCIATED(field%rlds) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rlus ) IF( ASSOCIATED(field%rlus) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rlutcs ) IF( ASSOCIATED(field%rlutcs) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rldscs ) IF( ASSOCIATED(field%rldscs) )
    !$ser verbatim   !$ACC UPDATE HOST( field%q_rlw_nlev ) IF( ASSOCIATED(field%q_rlw_nlev) )
    !$ser verbatim   !$ACC UPDATE HOST( field%qconv ) IF( ASSOCIATED(field%qconv) )
    !$ser verbatim   !$ACC UPDATE HOST( field%q_rad ) IF( ASSOCIATED(field%q_rad) )
    !$ser verbatim   !$ACC UPDATE HOST( field%q_rad_vi ) IF( ASSOCIATED(field%q_rad_vi) )
    !$ser verbatim   !$ACC UPDATE HOST( field%q_rsw ) IF( ASSOCIATED(field%q_rsw) )
    !$ser verbatim   !$ACC UPDATE HOST( field%q_rsw_vi ) IF( ASSOCIATED(field%q_rsw_vi) )
    !$ser verbatim   !$ACC UPDATE HOST( field%q_rlw ) IF( ASSOCIATED(field%q_rlw) )
    !$ser verbatim   !$ACC UPDATE HOST( field%q_rlw_vi ) IF( ASSOCIATED(field%q_rlw_vi) )
    !$ser verbatim   !$ACC UPDATE HOST( field%q_phy ) IF( ASSOCIATED(field%q_phy) )
    !$ser verbatim   !$ACC UPDATE HOST( field%q_phy_vi ) IF( ASSOCIATED(field%q_phy_vi) )
    !$ser verbatim   !$ACC UPDATE HOST( field%q_rlw_nlev ) IF( ASSOCIATED(field%q_rlw_nlev) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%ta_rad ) IF( ASSOCIATED(tend%ta_rad) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%ta_rsw ) IF( ASSOCIATED(tend%ta_rsw) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%ta_rlw ) IF( ASSOCIATED(tend%ta_rlw) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%ta_phy ) IF( ASSOCIATED(tend%ta_phy) )
    !$ser verbatim   !NOACC UPDATE HOST( emis_rad ) ! NOT ON GPU
#endif
    !$ser verbatim   CALL datetimeToString(time_config%tc_current_date, date)
    !$ser verbatim   CALL init('icon')
    !$ser savepoint echam_rht-output jg=jg jb=jb jcs=jcs jce=jce nproma=nproma nlev=nlev date=TRIM(date) 
    !$ser mode write
    !$ser data echam_rht_cosmu0=field%cosmu0(:,jb) IF (ASSOCIATED(field%cosmu0))
    !$ser data echam_rht_daylght_frc=field%daylght_frc(:,jb) IF (ASSOCIATED(field%daylght_frc))
    !$ser data echam_rht_ts_rad=field%ts_rad(:,jb) IF (ASSOCIATED(field%ts_rad))
    !$ser data echam_rht_ts_rad_rt=field%ts_rad_rt(:,jb) IF (ASSOCIATED(field%ts_rad_rt))
    !$ser data echam_rht_rsd_rt=field%rsd_rt(:,:,jb) IF (ASSOCIATED(field%rsd_rt))
    !$ser data echam_rht_rsu_rt=field%rsu_rt(:,:,jb) IF (ASSOCIATED(field%rsu_rt))
    !$ser data echam_rht_rsdcs_rt=field%rsdcs_rt(:,:,jb) IF (ASSOCIATED(field%rsdcs_rt))
    !$ser data echam_rht_rsucs_rt=field%rsucs_rt(:,:,jb) IF (ASSOCIATED(field%rsucs_rt))
    !$ser data echam_rht_rld_rt=field%rld_rt(:,:,jb) IF (ASSOCIATED(field%rld_rt))
    !$ser data echam_rht_rlu_rt=field%rlu_rt(:,:,jb) IF (ASSOCIATED(field%rlu_rt))
    !$ser data echam_rht_rldcs_rt=field%rldcs_rt(:,:,jb) IF (ASSOCIATED(field%rldcs_rt))
    !$ser data echam_rht_rlucs_rt=field%rlucs_rt(:,:,jb) IF (ASSOCIATED(field%rlucs_rt))
    !$ser data echam_rht_rvds_dir_rt=field%rvds_dir_rt(:,jb) IF (ASSOCIATED(field%rvds_dir_rt))
    !$ser data echam_rht_rpds_dir_rt=field%rpds_dir_rt(:,jb) IF (ASSOCIATED(field%rpds_dir_rt))
    !$ser data echam_rht_rnds_dir_rt=field%rnds_dir_rt(:,jb) IF (ASSOCIATED(field%rnds_dir_rt))
    !$ser data echam_rht_rvds_dif_rt=field%rvds_dif_rt(:,jb) IF (ASSOCIATED(field%rvds_dif_rt))
    !$ser data echam_rht_rpds_dif_rt=field%rpds_dif_rt(:,jb) IF (ASSOCIATED(field%rpds_dif_rt))
    !$ser data echam_rht_rnds_dif_rt=field%rnds_dif_rt(:,jb) IF (ASSOCIATED(field%rnds_dif_rt))
    !$ser data echam_rht_rvus_rt=field%rvus_rt(:,jb) IF (ASSOCIATED(field%rvus_rt))
    !$ser data echam_rht_rpus_rt=field%rpus_rt(:,jb) IF (ASSOCIATED(field%rpus_rt))
    !$ser data echam_rht_rnus_rt=field%rnus_rt(:,jb) IF (ASSOCIATED(field%rnus_rt))
    !$ser data echam_rht_rsdt=field%rsdt(:,jb) IF (ASSOCIATED(field%rsdt))
    !$ser data echam_rht_rsut=field%rsut(:,jb) IF (ASSOCIATED(field%rsut))
    !$ser data echam_rht_rsds=field%rsds(:,jb) IF (ASSOCIATED(field%rsds))
    !$ser data echam_rht_rsus=field%rsus(:,jb) IF (ASSOCIATED(field%rsus))
    !$ser data echam_rht_rsutcs=field%rsutcs(:,jb) IF (ASSOCIATED(field%rsutcs))
    !$ser data echam_rht_rsdscs=field%rsdscs(:,jb) IF (ASSOCIATED(field%rsdscs))
    !$ser data echam_rht_rsuscs=field%rsuscs(:,jb) IF (ASSOCIATED(field%rsuscs))
    !$ser data echam_rht_rvds_dir=field%rvds_dir(:,jb) IF (ASSOCIATED(field%rvds_dir))
    !$ser data echam_rht_rpds_dir=field%rpds_dir(:,jb) IF (ASSOCIATED(field%rpds_dir))
    !$ser data echam_rht_rnds_dir=field%rnds_dir(:,jb) IF (ASSOCIATED(field%rnds_dir))
    !$ser data echam_rht_rvds_dif=field%rvds_dif(:,jb) IF (ASSOCIATED(field%rvds_dif))
    !$ser data echam_rht_rpds_dif=field%rpds_dif(:,jb) IF (ASSOCIATED(field%rpds_dif))
    !$ser data echam_rht_rnds_dif=field%rnds_dif(:,jb) IF (ASSOCIATED(field%rnds_dif))
    !$ser data echam_rht_rvus=field%rvus(:,jb) IF (ASSOCIATED(field%rvus))
    !$ser data echam_rht_rpus=field%rpus(:,jb) IF (ASSOCIATED(field%rpus))
    !$ser data echam_rht_rnus=field%rnus(:,jb) IF (ASSOCIATED(field%rnus))
    !$ser data echam_rht_rlut=field%rlut(:,jb) IF (ASSOCIATED(field%rlut))
    !$ser data echam_rht_rlds=field%rlds(:,jb) IF (ASSOCIATED(field%rlds))
    !$ser data echam_rht_rlus=field%rlus(:,jb) IF (ASSOCIATED(field%rlus))
    !$ser data echam_rht_rlutcs=field%rlutcs(:,jb) IF (ASSOCIATED(field%rlutcs))
    !$ser data echam_rht_rldscs=field%rldscs(:,jb) IF (ASSOCIATED(field%rldscs))
    !$ser data echam_rht_q_rlw_nlev=field%q_rlw_nlev(:,jb) IF (ASSOCIATED(field%q_rlw_nlev))
    !$ser data echam_rht_qconv=field%qconv(:,:,jb) IF (ASSOCIATED(field%qconv))
    !$ser data echam_rht_q_rad=field%q_rad(:,:,jb) IF (ASSOCIATED(field%q_rad))
    !$ser data echam_rht_q_rad_vi=field%q_rad_vi(:,jb) IF (ASSOCIATED(field%q_rad_vi))
    !$ser data echam_rht_q_rsw=field%q_rsw(:,:,jb) IF (ASSOCIATED(field%q_rsw))
    !$ser data echam_rht_q_rsw_vi=field%q_rsw_vi(:,jb) IF (ASSOCIATED(field%q_rsw_vi))
    !$ser data echam_rht_q_rlw=field%q_rlw(:,:,jb) IF (ASSOCIATED(field%q_rlw))
    !$ser data echam_rht_q_rlw_vi=field%q_rlw_vi(:,jb) IF (ASSOCIATED(field%q_rlw_vi))
    !$ser data echam_rht_q_phy=field%q_phy(:,:,jb) IF (ASSOCIATED(field%q_phy))
    !$ser data echam_rht_q_phy_vi=field%q_phy_vi(:,jb) IF (ASSOCIATED(field%q_phy_vi))
    !$ser data echam_rht_ta_rad=tend%ta_rad(:,:,jb) IF (ASSOCIATED(tend%ta_rad))
    !$ser data echam_rht_ta_rsw=tend%ta_rsw(:,:,jb) IF (ASSOCIATED(tend%ta_rsw))
    !$ser data echam_rht_ta_rlw=tend%ta_rlw(:,:,jb) IF (ASSOCIATED(tend%ta_rlw))
    !$ser data echam_rht_ta_phy=tend%ta_phy(:,:,jb) IF (ASSOCIATED(tend%ta_phy))
    !$ser data echam_rht_ta=tend%ta(:,:,jb) IF (ASSOCIATED(tend%ta))
    !$ser data echam_rht_emis_rad=emis_rad(:,jb)
    !$ser verbatim lactive = .FALSE.
    !$ser verbatim IF (lonlyonce) THEN
    !$ser verbatim   lenabled = .FALSE.
    !$ser verbatim END IF
    !$ser verbatim ENDIF

  END SUBROUTINE serialize_rht_output

END MODULE mo_ser_echam_rht
