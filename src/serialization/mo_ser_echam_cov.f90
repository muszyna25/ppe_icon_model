!--------------------------------------------------------------------
!
! Serialization routine for ECHAM COV
!
!--------------------------------------------------------------------

MODULE mo_ser_echam_cov

  USE mo_kind,               ONLY: vp, wp
  USE mo_exception,          ONLY: warning
  USE mo_ser_common,         ONLY: init
  USE mtime,                 ONLY: datetimeToString, MAX_DATETIME_STR_LEN
  USE mo_time_config,        ONLY: time_config
  USE mo_echam_phy_memory,   ONLY: t_echam_phy_field, t_echam_phy_tend
  IMPLICIT NONE

  INTEGER, PARAMETER :: selected_block = 1

  PUBLIC :: serialize_cov_input
  PUBLIC :: serialize_cov_output

  CONTAINS

  SUBROUTINE serialize_cov_input(jg, jb, jcs, jce, nproma, nlev, field)
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
    !$ser verbatim   CALL warning('SER:mo_ser_echam_cov:input','Serialization is active!')
#if defined( _OPENACC )
    !$ser verbatim   CALL warning('GPU:mo_ser_echam_cov:input','GPU HOST synchronization forced by serialization!')
    !$ser verbatim   !$ACC UPDATE HOST( field%rtype ) IF( ASSOCIATED(field%rtype) )
    !$ser verbatim   !$ACC UPDATE HOST( field%frac_tile ) IF( ASSOCIATED(field%frac_tile) )
    !$ser verbatim   !$ACC UPDATE HOST( field%zf ) IF( ASSOCIATED(field%zf) )
    !$ser verbatim   !$ACC UPDATE HOST( field%presi_old ) IF( ASSOCIATED(field%presi_old) )
    !$ser verbatim   !$ACC UPDATE HOST( field%presm_old ) IF( ASSOCIATED(field%presm_old) )
    !$ser verbatim   !$ACC UPDATE HOST( field%ta ) IF( ASSOCIATED(field%ta) )
    !$ser verbatim   !$ACC UPDATE HOST( field%qtrc ) IF( ASSOCIATED(field%qtrc) )
    !$ser verbatim   !$ACC UPDATE HOST( field%aclc ) IF( ASSOCIATED(field%aclc) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rintop ) IF( ASSOCIATED(field%rintop) )
#endif
    !$ser verbatim   CALL datetimeToString(time_config%tc_current_date, date)
    !$ser verbatim   CALL init('icon')
    !$ser savepoint echam_cov-input jg=jg jb=jb jcs=jcs jce=jce nproma=nproma nlev=nlev date=TRIM(date)
#if defined( SERIALIZE_CREATE_REFERENCE )
    !$ser mode write
#elif defined( SERIALIZE_PERTURB_REFERENCE )
    !$ser mode read-perturb
#elif defined( SERIALIZE_READ_REFERENCE )
    !$ser mode read
#else
#error SERIALIZATION MODE IS NOT SET
#endif
    !$ser data echam_cov_rtype=field%rtype(:,jb) IF (ASSOCIATED(field%rtype))
    !$ser data echam_cov_frac_tile=field%frac_tile(:,jb,:) IF (ASSOCIATED(field%frac_tile))
    !$ser data echam_cov_zf=field%zf(:,:,jb) IF (ASSOCIATED(field%zf))
    !$ser data echam_cov_presi_old=field%presi_old(:,:,jb) IF (ASSOCIATED(field%presi_old))
    !$ser data echam_cov_presm_old=field%presm_old(:,:,jb) IF (ASSOCIATED(field%presm_old))
    !$ser data echam_cov_ta=field%ta(:,:,jb) IF (ASSOCIATED(field%ta))
    !$ser data echam_cov_qtrc=field%qtrc(:,:,jb,:) IF (ASSOCIATED(field%qtrc))
    !$ser data echam_cov_aclc=field%aclc(:,:,jb) IF (ASSOCIATED(field%aclc))
    !$ser data echam_cov_rintop=field%rintop(:,jb) IF (ASSOCIATED(field%rintop))
    !$ser verbatim lactive = .FALSE.
    !$ser verbatim IF (lonlyonce) THEN
    !$ser verbatim   lenabled = .FALSE.
    !$ser verbatim END IF
#if defined( SERIALIZE_READ_REFERENCE ) || defined( SERIALIZE_PERTURB_REFERENCE )
#if defined( _OPENACC )
    !$ser verbatim CALL warning('GPU:mo_ser_echam_cov:input','GPU DEVICE synchronization forced by serialization!')
    !$ser verbatim   !$ACC UPDATE DEVICE( field%rtype ) IF( ASSOCIATED(field%rtype) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%frac_tile ) IF( ASSOCIATED(field%frac_tile) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%zf ) IF( ASSOCIATED(field%zf) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%presi_old ) IF( ASSOCIATED(field%presi_old) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%presm_old ) IF( ASSOCIATED(field%presm_old) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%ta ) IF( ASSOCIATED(field%ta) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%qtrc ) IF( ASSOCIATED(field%qtrc) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%aclc ) IF( ASSOCIATED(field%aclc) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%rintop ) IF( ASSOCIATED(field%rintop) )
#endif
#endif
    !$ser verbatim ENDIF

  END SUBROUTINE serialize_cov_input

  SUBROUTINE serialize_cov_output(jg, jb, jcs, jce, nproma, nlev, field)
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
    !$ser verbatim   CALL warning('SER:mo_ser_echam_cov:output','Serialization is active!')
#if defined( _OPENACC )
    !$ser verbatim   CALL warning('GPU:mo_ser_echam_cov:output','GPU HOST synchronization forced by serialization!')
    !$ser verbatim   !$ACC UPDATE HOST( field%rtype ) IF( ASSOCIATED(field%rtype) )
    !$ser verbatim   !$ACC UPDATE HOST( field%frac_tile ) IF( ASSOCIATED(field%frac_tile) )
    !$ser verbatim   !$ACC UPDATE HOST( field%zf ) IF( ASSOCIATED(field%zf) )
    !$ser verbatim   !$ACC UPDATE HOST( field%presi_old ) IF( ASSOCIATED(field%presi_old) )
    !$ser verbatim   !$ACC UPDATE HOST( field%presm_old ) IF( ASSOCIATED(field%presm_old) )
    !$ser verbatim   !$ACC UPDATE HOST( field%ta ) IF( ASSOCIATED(field%ta) )
    !$ser verbatim   !$ACC UPDATE HOST( field%qtrc ) IF( ASSOCIATED(field%qtrc) )
    !$ser verbatim   !$ACC UPDATE HOST( field%aclc ) IF( ASSOCIATED(field%aclc) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rintop ) IF( ASSOCIATED(field%rintop) )
#endif
    !$ser verbatim   CALL datetimeToString(time_config%tc_current_date, date)
    !$ser verbatim   CALL init('icon')
    !$ser savepoint echam_cov-output jg=jg jb=jb jcs=jcs jce=jce nproma=nproma nlev=nlev date=TRIM(date) 
    !$ser mode write
    !$ser data echam_cov_rtype=field%rtype(:,jb) IF (ASSOCIATED(field%rtype))
    !$ser data echam_cov_frac_tile=field%frac_tile(:,jb,:) IF (ASSOCIATED(field%frac_tile))
    !$ser data echam_cov_zf=field%zf(:,:,jb) IF (ASSOCIATED(field%zf))
    !$ser data echam_cov_presi_old=field%presi_old(:,:,jb) IF (ASSOCIATED(field%presi_old))
    !$ser data echam_cov_presm_old=field%presm_old(:,:,jb) IF (ASSOCIATED(field%presm_old))
    !$ser data echam_cov_ta=field%ta(:,:,jb) IF (ASSOCIATED(field%ta))
    !$ser data echam_cov_qtrc=field%qtrc(:,:,jb,:) IF (ASSOCIATED(field%qtrc))
    !$ser data echam_cov_aclc=field%aclc(:,:,jb) IF (ASSOCIATED(field%aclc))
    !$ser data echam_cov_rintop=field%rintop(:,jb) IF (ASSOCIATED(field%rintop))
    !$ser verbatim lactive = .FALSE.
    !$ser verbatim IF (lonlyonce) THEN
    !$ser verbatim   lenabled = .FALSE.
    !$ser verbatim END IF
    !$ser verbatim ENDIF

  END SUBROUTINE serialize_cov_output

END MODULE mo_ser_echam_cov
