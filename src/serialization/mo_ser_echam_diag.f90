!--------------------------------------------------------------------
!
! Serialization routine for ECHAM diag
!
!--------------------------------------------------------------------

MODULE mo_ser_echam_diag

  USE mo_kind,               ONLY: vp, wp
  USE mo_exception,          ONLY: warning
  USE mo_ser_common,         ONLY: init
  USE mtime,                 ONLY: datetimeToString, MAX_DATETIME_STR_LEN
  USE mo_time_config,        ONLY: time_config
  USE mo_echam_phy_memory,   ONLY: t_echam_phy_field, t_echam_phy_tend
  IMPLICIT NONE

  INTEGER, PARAMETER :: selected_block = 1

  PUBLIC :: serialize_initialize_input
  PUBLIC :: serialize_initialize_output
  PUBLIC :: serialize_fractions_input
  PUBLIC :: serialize_fractions_output
  PUBLIC :: serialize_droplet_number_input
  PUBLIC :: serialize_droplet_number_output
  PUBLIC :: serialize_cpair_cvair_qconv_input
  PUBLIC :: serialize_cpair_cvair_qconv_output
  PUBLIC :: serialize_finalize_input
  PUBLIC :: serialize_finalize_output

  CONTAINS

  SUBROUTINE serialize_initialize_input(jg, jb, jcs, jce, nproma, nlev, field)
    INTEGER, INTENT(IN)    :: jg, jb, jcs, jce, nproma, nlev
    TYPE(t_echam_phy_field), POINTER, INTENT(INOUT) :: field

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: date
    LOGICAL, PARAMETER :: lonlyonce = .TRUE.
    LOGICAL, SAVE :: lenabled = .TRUE.
    LOGICAL, SAVE :: lactive = .FALSE.

    !$ser verbatim IF (selected_block < 0 .OR. jb == selected_block) THEN
    !$ser verbatim   lactive = .TRUE.
    !$ser verbatim ENDIF
    !$ser verbatim IF (lenabled .and. lactive) THEN
    !$ser verbatim   CALL warning('SER:mo_ser_echam_diag_ini:input','Serialization is active!')
#if defined( _OPENACC )
    !$ser verbatim   CALL warning('GPU:mo_ser_echam_diag_ini:input','GPU HOST synchronization forced by serialization!')
    !$ser verbatim   !$ACC UPDATE HOST( field%q_phy ) IF( ASSOCIATED(field%q_phy) )
    !$ser verbatim   !$ACC UPDATE HOST( field%q_phy_vi ) IF( ASSOCIATED(field%q_phy_vi) )
#endif
    !$ser verbatim   CALL datetimeToString(time_config%tc_current_date, date)
    !$ser verbatim   CALL init('icon')
    !$ser savepoint echam_diag_ini-input jg=jg jb=jb jcs=jcs jce=jce nproma=nproma nlev=nlev date=TRIM(date)
#if defined( SERIALIZE_CREATE_REFERENCE )
    !$ser mode write
#elif defined( SERIALIZE_PERTURB_REFERENCE )
    !$ser mode read-perturb
#elif defined( SERIALIZE_READ_REFERENCE )
    !$ser mode read
#else
#error SERIALIZATION MODE IS NOT SET
#endif
    !$ser data echam_diag_ini_q_phy=field%q_phy(:,:,jb) IF (ASSOCIATED(field%q_phy))
    !$ser data echam_diag_ini_q_phy_vi=field%q_phy_vi(:,jb) IF (ASSOCIATED(field%q_phy_vi))
    !$ser verbatim lactive = .FALSE.
    !$ser verbatim IF (lonlyonce) THEN
    !$ser verbatim   lenabled = .FALSE.
    !$ser verbatim END IF
#if defined( SERIALIZE_READ_REFERENCE ) || defined( SERIALIZE_PERTURB_REFERENCE )
#if defined( _OPENACC )
    !$ser verbatim CALL warning('GPU:mo_ser_echam_diag_ini:input','GPU DEVICE synchronization forced by serialization!')
    !$ser verbatim !$ACC UPDATE DEVICE( field%q_phy ) IF( ASSOCIATED(field%q_phy) )
    !$ser verbatim !$ACC UPDATE DEVICE( field%q_phy_vi ) IF( ASSOCIATED(field%q_phy_vi) )
#endif
#endif
    !$ser verbatim ENDIF

  END SUBROUTINE serialize_initialize_input

  SUBROUTINE serialize_initialize_output(jg, jb, jcs, jce, nproma, nlev, field)
    INTEGER, INTENT(IN)    :: jg, jb, jcs, jce, nproma, nlev
    TYPE(t_echam_phy_field), POINTER, INTENT(INOUT) :: field

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: date
    LOGICAL, PARAMETER :: lonlyonce = .TRUE.
    LOGICAL, SAVE :: lenabled = .TRUE.
    LOGICAL, SAVE :: lactive = .FALSE.

    !$ser verbatim IF (selected_block < 0 .OR. jb == selected_block) THEN
    !$ser verbatim   lactive = .TRUE.
    !$ser verbatim ENDIF
    !$ser verbatim IF (lenabled .and. lactive) THEN
    !$ser verbatim   CALL warning('SER:mo_ser_echam_diag_ini:output','Serialization is active!')
#if defined( _OPENACC )
    !$ser verbatim   CALL warning('GPU:mo_ser_echam_diag_ini:output','GPU HOST synchronization forced by serialization!')
    !$ser verbatim   !$ACC UPDATE HOST( field%q_phy ) IF( ASSOCIATED(field%q_phy) )
    !$ser verbatim   !$ACC UPDATE HOST( field%q_phy_vi ) IF( ASSOCIATED(field%q_phy_vi) )
#endif
    !$ser verbatim   CALL datetimeToString(time_config%tc_current_date, date)
    !$ser verbatim   CALL init('icon')
    !$ser savepoint echam_diag_ini-output jg=jg jb=jb jcs=jcs jce=jce nproma=nproma nlev=nlev date=TRIM(date) 
    !$ser mode write
    !$ser data echam_diag_ini_q_phy=field%q_phy(:,:,jb) IF (ASSOCIATED(field%q_phy))
    !$ser data echam_diag_ini_q_phy_vi=field%q_phy_vi(:,jb) IF (ASSOCIATED(field%q_phy_vi))
    !$ser verbatim lactive = .FALSE.
    !$ser verbatim IF (lonlyonce) THEN
    !$ser verbatim   lenabled = .FALSE.
    !$ser verbatim END IF
    !$ser verbatim ENDIF

  END SUBROUTINE serialize_initialize_output

  SUBROUTINE serialize_fractions_input(jg, jb, jcs, jce, nproma, nlev, field)
    INTEGER, INTENT(IN)    :: jg, jb, jcs, jce, nproma, nlev
    TYPE(t_echam_phy_field), POINTER, INTENT(INOUT) :: field

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: date
    LOGICAL, PARAMETER :: lonlyonce = .TRUE.
    LOGICAL, SAVE :: lenabled = .TRUE.
    LOGICAL, SAVE :: lactive = .FALSE.

    !$ser verbatim IF (selected_block < 0 .OR. jb == selected_block) THEN
    !$ser verbatim   lactive = .TRUE.
    !$ser verbatim ENDIF
    !$ser verbatim IF (lenabled .and. lactive) THEN
    !$ser verbatim   CALL warning('SER:mo_ser_echam_diag_fra:input','Serialization is active!')
#if defined( _OPENACC )
    !$ser verbatim   CALL warning('GPU:mo_ser_echam_diag_fra:input','GPU HOST synchronization forced by serialization!')
    !$ser verbatim   !$ACC UPDATE HOST( field%lsmask ) IF( ASSOCIATED(field%lsmask) )
    !$ser verbatim   !$ACC UPDATE HOST( field%alake ) IF( ASSOCIATED(field%alake) )
    !$ser verbatim   !$ACC UPDATE HOST( field%seaice ) IF( ASSOCIATED(field%seaice) )
    !$ser verbatim   !$ACC UPDATE HOST( field%lake_ice_frc ) IF( ASSOCIATED(field%lake_ice_frc) )
    !$ser verbatim   !$ACC UPDATE HOST( field%ts_tile ) IF( ASSOCIATED(field%ts_tile) )
    !$ser verbatim   !$ACC UPDATE HOST( field%frac_tile ) IF( ASSOCIATED(field%frac_tile) )
#endif
    !$ser verbatim   CALL datetimeToString(time_config%tc_current_date, date)
    !$ser verbatim   CALL init('icon')
    !$ser savepoint echam_diag_fra-input jg=jg jb=jb jcs=jcs jce=jce nproma=nproma nlev=nlev date=TRIM(date)
#if defined( SERIALIZE_CREATE_REFERENCE )
    !$ser mode write
#elif defined( SERIALIZE_PERTURB_REFERENCE )
    !$ser mode read-perturb
#elif defined( SERIALIZE_READ_REFERENCE )
    !$ser mode read
#else
#error SERIALIZATION MODE IS NOT SET
#endif
    !$ser data echam_diag_fra_lsmask=field%lsmask(:,jb) IF (ASSOCIATED(field%lsmask))
    !$ser data echam_diag_fra_alake=field%alake(:,jb) IF (ASSOCIATED(field%alake))
    !$ser data echam_diag_fra_seaice=field%seaice(:,jb) IF (ASSOCIATED(field%seaice))
    !$ser data echam_diag_fra_lake_ice_frc=field%lake_ice_frc(:,jb) IF (ASSOCIATED(field%lake_ice_frc))
    !$ser data echam_diag_fra_ts_tile=field%ts_tile(:,jb,:) IF (ASSOCIATED(field%ts_tile))
    !$ser data echam_diag_fra_frac_tile=field%frac_tile(:,jb,:) IF (ASSOCIATED(field%frac_tile))
    !$ser verbatim lactive = .FALSE.
    !$ser verbatim IF (lonlyonce) THEN
    !$ser verbatim   lenabled = .FALSE.
    !$ser verbatim END IF
#if defined( SERIALIZE_READ_REFERENCE ) || defined( SERIALIZE_PERTURB_REFERENCE )
#if defined( _OPENACC )
    !$ser verbatim CALL warning('GPU:mo_ser_echam_diag_fra:input','GPU DEVICE synchronization forced by serialization!')
    !$ser verbatim !$ACC UPDATE DEVICE( field%lsmask ) IF( ASSOCIATED(field%lsmask) )
    !$ser verbatim !$ACC UPDATE DEVICE( field%alake ) IF( ASSOCIATED(field%alake) )
    !$ser verbatim !$ACC UPDATE DEVICE( field%seaice ) IF( ASSOCIATED(field%seaice) )
    !$ser verbatim !$ACC UPDATE DEVICE( field%lake_ice_frc ) IF( ASSOCIATED(field%lake_ice_frc) )
    !$ser verbatim !$ACC UPDATE DEVICE( field%ts_tile ) IF( ASSOCIATED(field%ts_tile) )
    !$ser verbatim !$ACC UPDATE DEVICE( field%frac_tile ) IF( ASSOCIATED(field%frac_tile) )
#endif
#endif
    !$ser verbatim ENDIF

  END SUBROUTINE serialize_fractions_input

  SUBROUTINE serialize_fractions_output(jg, jb, jcs, jce, nproma, nlev, field)
    INTEGER, INTENT(IN)    :: jg, jb, jcs, jce, nproma, nlev
    TYPE(t_echam_phy_field), POINTER, INTENT(INOUT) :: field

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: date
    LOGICAL, PARAMETER :: lonlyonce = .TRUE.
    LOGICAL, SAVE :: lenabled = .TRUE.
    LOGICAL, SAVE :: lactive = .FALSE.

    !$ser verbatim IF (selected_block < 0 .OR. jb == selected_block) THEN
    !$ser verbatim   lactive = .TRUE.
    !$ser verbatim ENDIF
    !$ser verbatim IF (lenabled .and. lactive) THEN
    !$ser verbatim   CALL warning('SER:mo_ser_echam_diag_fra:output','Serialization is active!')
#if defined( _OPENACC )
    !$ser verbatim   CALL warning('GPU:mo_ser_echam_diag_fra:output','GPU HOST synchronization forced by serialization!')
    !$ser verbatim   !$ACC UPDATE HOST( field%lsmask ) IF( ASSOCIATED(field%lsmask) )
    !$ser verbatim   !$ACC UPDATE HOST( field%alake ) IF( ASSOCIATED(field%alake) )
    !$ser verbatim   !$ACC UPDATE HOST( field%seaice ) IF( ASSOCIATED(field%seaice) )
    !$ser verbatim   !$ACC UPDATE HOST( field%lake_ice_frc ) IF( ASSOCIATED(field%lake_ice_frc) )
    !$ser verbatim   !$ACC UPDATE HOST( field%ts_tile ) IF( ASSOCIATED(field%ts_tile) )
    !$ser verbatim   !$ACC UPDATE HOST( field%frac_tile ) IF( ASSOCIATED(field%frac_tile) )
#endif
    !$ser verbatim   CALL datetimeToString(time_config%tc_current_date, date)
    !$ser verbatim   CALL init('icon')
    !$ser savepoint echam_diag_fra-output jg=jg jb=jb jcs=jcs jce=jce nproma=nproma nlev=nlev date=TRIM(date) 
    !$ser mode write
    !$ser data echam_diag_fra_lsmask=field%lsmask(:,jb) IF (ASSOCIATED(field%lsmask))
    !$ser data echam_diag_fra_alake=field%alake(:,jb) IF (ASSOCIATED(field%alake))
    !$ser data echam_diag_fra_seaice=field%seaice(:,jb) IF (ASSOCIATED(field%seaice))
    !$ser data echam_diag_fra_lake_ice_frc=field%lake_ice_frc(:,jb) IF (ASSOCIATED(field%lake_ice_frc))
    !$ser data echam_diag_fra_ts_tile=field%ts_tile(:,jb,:) IF (ASSOCIATED(field%ts_tile))
    !$ser data echam_diag_fra_frac_tile=field%frac_tile(:,jb,:) IF (ASSOCIATED(field%frac_tile))
    !$ser verbatim lactive = .FALSE.
    !$ser verbatim IF (lonlyonce) THEN
    !$ser verbatim   lenabled = .FALSE.
    !$ser verbatim END IF
    !$ser verbatim ENDIF

  END SUBROUTINE serialize_fractions_output

  SUBROUTINE serialize_droplet_number_input(jg, jb, jcs, jce, nproma, nlev, field)
    INTEGER, INTENT(IN)    :: jg, jb, jcs, jce, nproma, nlev
    TYPE(t_echam_phy_field), POINTER, INTENT(INOUT) :: field

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: date
    LOGICAL, PARAMETER :: lonlyonce = .TRUE.
    LOGICAL, SAVE :: lenabled = .TRUE.
    LOGICAL, SAVE :: lactive = .FALSE.

    !$ser verbatim IF (selected_block < 0 .OR. jb == selected_block) THEN
    !$ser verbatim   lactive = .TRUE.
    !$ser verbatim ENDIF
    !$ser verbatim IF (lenabled .and. lactive) THEN
    !$ser verbatim   CALL warning('SER:mo_ser_echam_diag_dpn:input','Serialization is active!')
#if defined( _OPENACC )
    !$ser verbatim   CALL warning('GPU:mo_ser_echam_diag_dpn:input','GPU HOST synchronization forced by serialization!')
    !$ser verbatim   !$ACC UPDATE HOST( field%sftlf ) IF( ASSOCIATED(field%sftlf) )
    !$ser verbatim   !$ACC UPDATE HOST( field%sftgif ) IF( ASSOCIATED(field%sftgif) )
    !$ser verbatim   !$ACC UPDATE HOST( field%presm_old ) IF( ASSOCIATED(field%presm_old) )
    !$ser verbatim   !$ACC UPDATE HOST( field%acdnc ) IF( ASSOCIATED(field%acdnc) )
#endif
    !$ser verbatim   CALL datetimeToString(time_config%tc_current_date, date)
    !$ser verbatim   CALL init('icon')
    !$ser savepoint echam_diag_dpn-input jg=jg jb=jb jcs=jcs jce=jce nproma=nproma nlev=nlev date=TRIM(date)
#if defined( SERIALIZE_CREATE_REFERENCE )
    !$ser mode write
#elif defined( SERIALIZE_PERTURB_REFERENCE )
    !$ser mode read-perturb
#elif defined( SERIALIZE_READ_REFERENCE )
    !$ser mode read
#else
#error SERIALIZATION MODE IS NOT SET
#endif
    !$ser data echam_diag_dpn_sftlf=field%sftlf(:,jb) IF (ASSOCIATED(field%sftlf))
    !$ser data echam_diag_dpn_sftgif=field%sftgif(:,jb) IF (ASSOCIATED(field%sftgif))
    !$ser data echam_diag_dpn_presm_old=field%presm_old(:,:,jb) IF (ASSOCIATED(field%presm_old))
    !$ser data echam_diag_dpn_acdnc=field%acdnc(:,:,jb) IF (ASSOCIATED(field%acdnc))
    !$ser verbatim lactive = .FALSE.
    !$ser verbatim IF (lonlyonce) THEN
    !$ser verbatim   lenabled = .FALSE.
    !$ser verbatim END IF
#if defined( SERIALIZE_READ_REFERENCE ) || defined( SERIALIZE_PERTURB_REFERENCE )
#if defined( _OPENACC )
    !$ser verbatim CALL warning('GPU:mo_ser_echam_diag_dpn:input','GPU DEVICE synchronization forced by serialization!')
    !$ser verbatim !$ACC UPDATE DEVICE( field%sftlf ) IF( ASSOCIATED(field%sftlf) )
    !$ser verbatim !$ACC UPDATE DEVICE( field%sftgif ) IF( ASSOCIATED(field%sftgif) )
    !$ser verbatim !$ACC UPDATE DEVICE( field%presm_old ) IF( ASSOCIATED(field%presm_old) )
    !$ser verbatim !$ACC UPDATE DEVICE( field%acdnc ) IF( ASSOCIATED(field%acdnc) )
#endif
#endif
    !$ser verbatim ENDIF

  END SUBROUTINE serialize_droplet_number_input

  SUBROUTINE serialize_droplet_number_output(jg, jb, jcs, jce, nproma, nlev, field)
    INTEGER, INTENT(IN)    :: jg, jb, jcs, jce, nproma, nlev
    TYPE(t_echam_phy_field), POINTER, INTENT(INOUT) :: field

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: date
    LOGICAL, PARAMETER :: lonlyonce = .TRUE.
    LOGICAL, SAVE :: lenabled = .TRUE.
    LOGICAL, SAVE :: lactive = .FALSE.

    !$ser verbatim IF (selected_block < 0 .OR. jb == selected_block) THEN
    !$ser verbatim   lactive = .TRUE.
    !$ser verbatim ENDIF
    !$ser verbatim IF (lenabled .and. lactive) THEN
    !$ser verbatim   CALL warning('SER:mo_ser_echam_diag_dpn:output','Serialization is active!')
#if defined( _OPENACC )
    !$ser verbatim   CALL warning('GPU:mo_ser_echam_diag_dpn:output','GPU HOST synchronization forced by serialization!')
    !$ser verbatim   !$ACC UPDATE HOST( field%sftlf ) IF( ASSOCIATED(field%sftlf) )
    !$ser verbatim   !$ACC UPDATE HOST( field%sftgif ) IF( ASSOCIATED(field%sftgif) )
    !$ser verbatim   !$ACC UPDATE HOST( field%presm_old ) IF( ASSOCIATED(field%presm_old) )
    !$ser verbatim   !$ACC UPDATE HOST( field%acdnc ) IF( ASSOCIATED(field%acdnc) )
#endif
    !$ser verbatim   CALL datetimeToString(time_config%tc_current_date, date)
    !$ser verbatim   CALL init('icon')
    !$ser savepoint echam_diag_dpn-output jg=jg jb=jb jcs=jcs jce=jce nproma=nproma nlev=nlev date=TRIM(date) 
    !$ser mode write
    !$ser data echam_diag_dpn_sftlf=field%sftlf(:,jb) IF (ASSOCIATED(field%sftlf))
    !$ser data echam_diag_dpn_sftgif=field%sftgif(:,jb) IF (ASSOCIATED(field%sftgif))
    !$ser data echam_diag_dpn_presm_old=field%presm_old(:,:,jb) IF (ASSOCIATED(field%presm_old))
    !$ser data echam_diag_dpn_acdnc=field%acdnc(:,:,jb) IF (ASSOCIATED(field%acdnc))
    !$ser verbatim lactive = .FALSE.
    !$ser verbatim IF (lonlyonce) THEN
    !$ser verbatim   lenabled = .FALSE.
    !$ser verbatim END IF
    !$ser verbatim ENDIF

  END SUBROUTINE serialize_droplet_number_output

  SUBROUTINE serialize_cpair_cvair_qconv_input(jg, jb, jcs, jce, nproma, nlev, field)
    INTEGER, INTENT(IN)    :: jg, jb, jcs, jce, nproma, nlev
    TYPE(t_echam_phy_field), POINTER, INTENT(INOUT) :: field

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: date
    LOGICAL, PARAMETER :: lonlyonce = .TRUE.
    LOGICAL, SAVE :: lenabled = .TRUE.
    LOGICAL, SAVE :: lactive = .FALSE.

    !$ser verbatim IF (selected_block < 0 .OR. jb == selected_block) THEN
    !$ser verbatim   lactive = .TRUE.
    !$ser verbatim ENDIF
    !$ser verbatim IF (lenabled .and. lactive) THEN
    !$ser verbatim   CALL warning('SER:mo_ser_echam_diag_ccq:input','Serialization is active!')
#if defined( _OPENACC )
    !$ser verbatim   CALL warning('GPU:mo_ser_echam_diag_ccq:input','GPU HOST synchronization forced by serialization!')
    !$ser verbatim   !$ACC UPDATE HOST( field%cpair ) IF( ASSOCIATED(field%cpair) )
    !$ser verbatim   !$ACC UPDATE HOST( field%qtrc ) IF( ASSOCIATED(field%qtrc) )
    !$ser verbatim   !$ACC UPDATE HOST( field%cvair ) IF( ASSOCIATED(field%cvair) )
    !$ser verbatim   !$ACC UPDATE HOST( field%qconv ) IF( ASSOCIATED(field%qconv) )
    !$ser verbatim   !$ACC UPDATE HOST( field%mair ) IF( ASSOCIATED(field%mair) )
#endif
    !$ser verbatim   CALL datetimeToString(time_config%tc_current_date, date)
    !$ser verbatim   CALL init('icon')
    !$ser savepoint echam_diag_ccq-input jg=jg jb=jb jcs=jcs jce=jce nproma=nproma nlev=nlev date=TRIM(date)
#if defined( SERIALIZE_CREATE_REFERENCE )
    !$ser mode write
#elif defined( SERIALIZE_PERTURB_REFERENCE )
    !$ser mode read-perturb
#elif defined( SERIALIZE_READ_REFERENCE )
    !$ser mode read
#else
#error SERIALIZATION MODE IS NOT SET
#endif
    !$ser data echam_diag_ccq_cpair=field%cpair(:,:,jb) IF (ASSOCIATED(field%cpair))
    !$ser data echam_diag_ccq_qtrc=field%qtrc(:,:,jb,:) IF (ASSOCIATED(field%qtrc))
    !$ser data echam_diag_ccq_cvair=field%cvair(:,:,jb) IF (ASSOCIATED(field%cvair))
    !$ser data echam_diag_ccq_qconv=field%qconv(:,:,jb) IF (ASSOCIATED(field%qconv))
    !$ser data echam_diag_ccq_mair=field%mair(:,:,jb) IF (ASSOCIATED(field%mair))
    !$ser verbatim lactive = .FALSE.
    !$ser verbatim IF (lonlyonce) THEN
    !$ser verbatim   lenabled = .FALSE.
    !$ser verbatim END IF
#if defined( SERIALIZE_READ_REFERENCE ) || defined( SERIALIZE_PERTURB_REFERENCE )
#if defined( _OPENACC )
    !$ser verbatim CALL warning('GPU:mo_ser_echam_diag_ccq:input','GPU DEVICE synchronization forced by serialization!')
    !$ser verbatim !$ACC UPDATE DEVICE( field%cpair ) IF( ASSOCIATED(field%cpair) )
    !$ser verbatim !$ACC UPDATE DEVICE( field%qtrc ) IF( ASSOCIATED(field%qtrc) )
    !$ser verbatim !$ACC UPDATE DEVICE( field%cvair ) IF( ASSOCIATED(field%cvair) )
    !$ser verbatim !$ACC UPDATE DEVICE( field%qconv ) IF( ASSOCIATED(field%qconv) )
    !$ser verbatim !$ACC UPDATE DEVICE( field%mair ) IF( ASSOCIATED(field%mair) )
#endif
#endif
    !$ser verbatim ENDIF

  END SUBROUTINE serialize_cpair_cvair_qconv_input

  SUBROUTINE serialize_cpair_cvair_qconv_output(jg, jb, jcs, jce, nproma, nlev, field)
    INTEGER, INTENT(IN)    :: jg, jb, jcs, jce, nproma, nlev
    TYPE(t_echam_phy_field), POINTER, INTENT(INOUT) :: field

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: date
    LOGICAL, PARAMETER :: lonlyonce = .TRUE.
    LOGICAL, SAVE :: lenabled = .TRUE.
    LOGICAL, SAVE :: lactive = .FALSE.

    !$ser verbatim IF (selected_block < 0 .OR. jb == selected_block) THEN
    !$ser verbatim   lactive = .TRUE.
    !$ser verbatim ENDIF
    !$ser verbatim IF (lenabled .and. lactive) THEN
    !$ser verbatim   CALL warning('SER:mo_ser_echam_diag_ccq:output','Serialization is active!')
#if defined( _OPENACC )
    !$ser verbatim   CALL warning('GPU:mo_ser_echam_diag_ccq:output','GPU HOST synchronization forced by serialization!')
    !$ser verbatim   !$ACC UPDATE HOST( field%cpair ) IF( ASSOCIATED(field%cpair) )
    !$ser verbatim   !$ACC UPDATE HOST( field%qtrc ) IF( ASSOCIATED(field%qtrc) )
    !$ser verbatim   !$ACC UPDATE HOST( field%cvair ) IF( ASSOCIATED(field%cvair) )
    !$ser verbatim   !$ACC UPDATE HOST( field%qconv ) IF( ASSOCIATED(field%qconv) )
    !$ser verbatim   !$ACC UPDATE HOST( field%mair ) IF( ASSOCIATED(field%mair) )
#endif
    !$ser verbatim   CALL datetimeToString(time_config%tc_current_date, date)
    !$ser verbatim   CALL init('icon')
    !$ser savepoint echam_diag_ccq-output jg=jg jb=jb jcs=jcs jce=jce nproma=nproma nlev=nlev date=TRIM(date) 
    !$ser mode write
    !$ser data echam_diag_ccq_cpair=field%cpair(:,:,jb) IF (ASSOCIATED(field%cpair))
    !$ser data echam_diag_ccq_qtrc=field%qtrc(:,:,jb,:) IF (ASSOCIATED(field%qtrc))
    !$ser data echam_diag_ccq_cvair=field%cvair(:,:,jb) IF (ASSOCIATED(field%cvair))
    !$ser data echam_diag_ccq_qconv=field%qconv(:,:,jb) IF (ASSOCIATED(field%qconv))
    !$ser data echam_diag_ccq_mair=field%mair(:,:,jb) IF (ASSOCIATED(field%mair))
    !$ser verbatim lactive = .FALSE.
    !$ser verbatim IF (lonlyonce) THEN
    !$ser verbatim   lenabled = .FALSE.
    !$ser verbatim END IF
    !$ser verbatim ENDIF

  END SUBROUTINE serialize_cpair_cvair_qconv_output

  SUBROUTINE serialize_finalize_input(jg, jb, jcs, jce, nproma, nlev, field, tend)
    INTEGER, INTENT(IN)    :: jg, jb, jcs, jce, nproma, nlev
    TYPE(t_echam_phy_field), POINTER, INTENT(INOUT) :: field
    TYPE(t_echam_phy_tend), POINTER, INTENT(INOUT) :: tend

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: date
    LOGICAL, PARAMETER :: lonlyonce = .TRUE.
    LOGICAL, SAVE :: lenabled = .TRUE.
    LOGICAL, SAVE :: lactive = .FALSE.

    !$ser verbatim IF (selected_block < 0 .OR. jb == selected_block) THEN
    !$ser verbatim   lactive = .TRUE.
    !$ser verbatim ENDIF
    !$ser verbatim IF (lenabled .and. lactive) THEN
    !$ser verbatim   CALL warning('SER:mo_ser_echam_diag_fin:input','Serialization is active!')
#if defined( _OPENACC )
    !$ser verbatim   CALL warning('GPU:mo_ser_echam_diag_fin:input','GPU HOST synchronization forced by serialization!')
    !$ser verbatim   !$ACC UPDATE HOST( field%pr ) IF( ASSOCIATED(field%pr) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rsfl ) IF( ASSOCIATED(field%rsfl) )
    !$ser verbatim   !$ACC UPDATE HOST( field%ssfl ) IF( ASSOCIATED(field%ssfl) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rsfc ) IF( ASSOCIATED(field%rsfc) )
    !$ser verbatim   !$ACC UPDATE HOST( field%ssfc ) IF( ASSOCIATED(field%ssfc) )
    !$ser verbatim   !$ACC UPDATE HOST( field%cpair ) IF( ASSOCIATED(field%cpair) )
    !$ser verbatim   !$ACC UPDATE HOST( field%cvair ) IF( ASSOCIATED(field%cvair) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%ta_phy ) IF( ASSOCIATED(tend%ta_phy) )
#endif
    !$ser verbatim   CALL datetimeToString(time_config%tc_current_date, date)
    !$ser verbatim   CALL init('icon')
    !$ser savepoint echam_diag_fin-input jg=jg jb=jb jcs=jcs jce=jce nproma=nproma nlev=nlev date=TRIM(date)
#if defined( SERIALIZE_CREATE_REFERENCE )
    !$ser mode write
#elif defined( SERIALIZE_PERTURB_REFERENCE )
    !$ser mode read-perturb
#elif defined( SERIALIZE_READ_REFERENCE )
    !$ser mode read
#else
#error SERIALIZATION MODE IS NOT SET
#endif
    !$ser data echam_diag_fin_pr=field%pr(:,jb) IF (ASSOCIATED(field%pr))
    !$ser data echam_diag_fin_rsfl=field%rsfl(:,jb) IF (ASSOCIATED(field%rsfl))
    !$ser data echam_diag_fin_ssfl=field%ssfl(:,jb) IF (ASSOCIATED(field%ssfl))
    !$ser data echam_diag_fin_rsfc=field%rsfc(:,jb) IF (ASSOCIATED(field%rsfc))
    !$ser data echam_diag_fin_ssfc=field%ssfc(:,jb) IF (ASSOCIATED(field%ssfc))
    !$ser data echam_diag_fin_cpair=field%cpair(:,:,jb) IF (ASSOCIATED(field%cpair))
    !$ser data echam_diag_fin_cvair=field%cvair(:,:,jb) IF (ASSOCIATED(field%cvair))
    !$ser data echam_diag_fin_ta_phy=tend%ta_phy(:,:,jb) IF (ASSOCIATED(tend%ta_phy))
    !$ser verbatim lactive = .FALSE.
    !$ser verbatim IF (lonlyonce) THEN
    !$ser verbatim   lenabled = .FALSE.
    !$ser verbatim END IF
#if defined( SERIALIZE_READ_REFERENCE ) || defined( SERIALIZE_PERTURB_REFERENCE )
#if defined( _OPENACC )
    !$ser verbatim CALL warning('GPU:mo_ser_echam_diag_fin:input','GPU DEVICE synchronization forced by serialization!')
    !$ser verbatim !$ACC UPDATE DEVICE( field%pr ) IF( ASSOCIATED(field%pr) )
    !$ser verbatim !$ACC UPDATE DEVICE( field%rsfl ) IF( ASSOCIATED(field%rsfl) )
    !$ser verbatim !$ACC UPDATE DEVICE( field%ssfl ) IF( ASSOCIATED(field%ssfl) )
    !$ser verbatim !$ACC UPDATE DEVICE( field%rsfc ) IF( ASSOCIATED(field%rsfc) )
    !$ser verbatim !$ACC UPDATE DEVICE( field%ssfc ) IF( ASSOCIATED(field%ssfc) )
    !$ser verbatim !$ACC UPDATE DEVICE( field%cpair ) IF( ASSOCIATED(field%cpair) )
    !$ser verbatim !$ACC UPDATE DEVICE( field%cvair ) IF( ASSOCIATED(field%cvair) )
    !$ser verbatim !$ACC UPDATE DEVICE( tend%ta_phy ) IF( ASSOCIATED(tend%ta_phy) )
#endif
#endif
    !$ser verbatim ENDIF

  END SUBROUTINE serialize_finalize_input

  SUBROUTINE serialize_finalize_output(jg, jb, jcs, jce, nproma, nlev, field, tend)
    INTEGER, INTENT(IN)    :: jg, jb, jcs, jce, nproma, nlev
    TYPE(t_echam_phy_field), POINTER, INTENT(INOUT) :: field
    TYPE(t_echam_phy_tend), POINTER, INTENT(INOUT) :: tend

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: date
    LOGICAL, PARAMETER :: lonlyonce = .TRUE.
    LOGICAL, SAVE :: lenabled = .TRUE.
    LOGICAL, SAVE :: lactive = .FALSE.

    !$ser verbatim IF (selected_block < 0 .OR. jb == selected_block) THEN
    !$ser verbatim   lactive = .TRUE.
    !$ser verbatim ENDIF
    !$ser verbatim IF (lenabled .and. lactive) THEN
    !$ser verbatim   CALL warning('SER:mo_ser_echam_diag_fin:output','Serialization is active!')
#if defined( _OPENACC )
    !$ser verbatim   CALL warning('GPU:mo_ser_echam_diag_fin:output','GPU HOST synchronization forced by serialization!')
    !$ser verbatim   !$ACC UPDATE HOST( field%pr ) IF( ASSOCIATED(field%pr) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rsfl ) IF( ASSOCIATED(field%rsfl) )
    !$ser verbatim   !$ACC UPDATE HOST( field%ssfl ) IF( ASSOCIATED(field%ssfl) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rsfc ) IF( ASSOCIATED(field%rsfc) )
    !$ser verbatim   !$ACC UPDATE HOST( field%ssfc ) IF( ASSOCIATED(field%ssfc) )
    !$ser verbatim   !$ACC UPDATE HOST( field%cpair ) IF( ASSOCIATED(field%cpair) )
    !$ser verbatim   !$ACC UPDATE HOST( field%cvair ) IF( ASSOCIATED(field%cvair) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%ta_phy ) IF( ASSOCIATED(tend%ta_phy) )
#endif
    !$ser verbatim   CALL datetimeToString(time_config%tc_current_date, date)
    !$ser verbatim   CALL init('icon')
    !$ser savepoint echam_diag_fin-output jg=jg jb=jb jcs=jcs jce=jce nproma=nproma nlev=nlev date=TRIM(date) 
    !$ser mode write
    !$ser data echam_diag_fin_pr=field%pr(:,jb) IF (ASSOCIATED(field%pr))
    !$ser data echam_diag_fin_rsfl=field%rsfl(:,jb) IF (ASSOCIATED(field%rsfl))
    !$ser data echam_diag_fin_ssfl=field%ssfl(:,jb) IF (ASSOCIATED(field%ssfl))
    !$ser data echam_diag_fin_rsfc=field%rsfc(:,jb) IF (ASSOCIATED(field%rsfc))
    !$ser data echam_diag_fin_ssfc=field%ssfc(:,jb) IF (ASSOCIATED(field%ssfc))
    !$ser data echam_diag_fin_cpair=field%cpair(:,:,jb) IF (ASSOCIATED(field%cpair))
    !$ser data echam_diag_fin_cvair=field%cvair(:,:,jb) IF (ASSOCIATED(field%cvair))
    !$ser data echam_diag_fin_ta_phy=tend%ta_phy(:,:,jb) IF (ASSOCIATED(tend%ta_phy))
    !$ser verbatim lactive = .FALSE.
    !$ser verbatim IF (lonlyonce) THEN
    !$ser verbatim   lenabled = .FALSE.
    !$ser verbatim END IF
    !$ser verbatim ENDIF

  END SUBROUTINE serialize_finalize_output

END MODULE mo_ser_echam_diag
