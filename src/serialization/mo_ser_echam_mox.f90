!--------------------------------------------------------------------
!
! Serialization routine for ECHAM MOX
!
!--------------------------------------------------------------------

MODULE mo_ser_echam_mox

  USE mo_kind,               ONLY: vp, wp
  USE mo_exception,          ONLY: warning
  USE mo_ser_common,         ONLY: init
  USE mtime,                 ONLY: datetimeToString, MAX_DATETIME_STR_LEN
  USE mo_time_config,        ONLY: time_config
  USE mo_echam_phy_memory,   ONLY: t_echam_phy_field, t_echam_phy_tend
  IMPLICIT NONE

  INTEGER, PARAMETER :: selected_block = 1

  PUBLIC :: serialize_mox_input
  PUBLIC :: serialize_mox_output

  CONTAINS

  SUBROUTINE serialize_mox_input(jg, jb, jcs, jce, nproma, nlev, field, tend)
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
    !$ser verbatim   CALL warning('SER:mo_ser_echam_mox:input','Serialization is active!')
#if defined( _OPENACC )
    !$ser verbatim   CALL warning('GPU:mo_ser_echam_mox:input','GPU HOST synchronization forced by serialization!')
    !$ser verbatim   !$ACC UPDATE HOST( field%presm_old ) IF(ASSOCIATED( field%presm_old) )
    !$ser verbatim   !$ACC UPDATE HOST( field%qtrc ) IF( ASSOCIATED(field%qtrc) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%qtrc_mox ) IF( ASSOCIATED(tend%qtrc_mox) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%qtrc_phy ) IF( ASSOCIATED(tend%qtrc_phy) )
#endif
    !$ser verbatim   CALL datetimeToString(time_config%tc_current_date, date)
    !$ser verbatim   CALL init('icon')
    !$ser savepoint echam_mox-input jg=jg jb=jb jcs=jcs jce=jce nproma=nproma nlev=nlev date=TRIM(date)
#if defined( SERIALIZE_CREATE_REFERENCE )
    !$ser mode write
#elif defined( SERIALIZE_PERTURB_REFERENCE )
    !$ser mode read-perturb
#elif defined( SERIALIZE_READ_REFERENCE )
    !$ser mode read
#else
#error SERIALIZATION MODE IS NOT SET
#endif
    !$ser data echam_mox_presm_old=field%presm_old(:,:,jb) IF (ASSOCIATED(field%presm_old))
    !$ser data echam_mox_qtrc=field%qtrc(:,:,jb,:) IF (ASSOCIATED(field%qtrc))
    !$ser data echam_mox_qtrc_mox=tend%qtrc_mox(:,:,jb,:) IF (ASSOCIATED(tend%qtrc_mox))
    !$ser data echam_mox_qtrc_phy=tend%qtrc_phy(:,:,jb,:) IF (ASSOCIATED(tend%qtrc_phy))
    !$ser verbatim lactive = .FALSE.
    !$ser verbatim IF (lonlyonce) THEN
    !$ser verbatim   lenabled = .FALSE.
    !$ser verbatim END IF
#if defined( SERIALIZE_READ_REFERENCE ) || defined( SERIALIZE_PERTURB_REFERENCE )
#if defined( _OPENACC )
    !$ser verbatim CALL warning('GPU:mo_ser_echam_mox:input','GPU DEVICE synchronization forced by serialization!')
    !$ser verbatim   !$ACC UPDATE DEVICE( field%presm_old ) IF(ASSOCIATED( field%presm_old) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%qtrc ) IF( ASSOCIATED(field%qtrc) )
    !$ser verbatim   !$ACC UPDATE DEVICE( tend%qtrc_mox ) IF( ASSOCIATED(tend%qtrc_mox) )
    !$ser verbatim   !$ACC UPDATE DEVICE( tend%qtrc_phy ) IF( ASSOCIATED(tend%qtrc_phy) )
#endif
#endif
    !$ser verbatim ENDIF

  END SUBROUTINE serialize_mox_input

  SUBROUTINE serialize_mox_output(jg, jb, jcs, jce, nproma, nlev, field, tend)
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
    !$ser verbatim   CALL warning('SER:mo_ser_echam_mox:output','Serialization is active!')
#if defined( _OPENACC )
    !$ser verbatim   CALL warning('GPU:mo_ser_echam_mox:output','GPU HOST synchronization forced by serialization!')
    !$ser verbatim   !$ACC UPDATE HOST( field%presm_old ) IF(ASSOCIATED( field%presm_old) )
    !$ser verbatim   !$ACC UPDATE HOST( field%qtrc ) IF( ASSOCIATED(field%qtrc) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%qtrc_mox ) IF( ASSOCIATED(tend%qtrc_mox) )
    !$ser verbatim   !$ACC UPDATE HOST( tend%qtrc_phy ) IF( ASSOCIATED(tend%qtrc_phy) )
#endif
    !$ser verbatim   CALL datetimeToString(time_config%tc_current_date, date)
    !$ser verbatim   CALL init('icon')
    !$ser savepoint echam_mox-output jg=jg jb=jb jcs=jcs jce=jce nproma=nproma nlev=nlev date=TRIM(date) 
    !$ser mode write
    !$ser data echam_mox_presm_old=field%presm_old(:,:,jb) IF (ASSOCIATED(field%presm_old))
    !$ser data echam_mox_qtrc=field%qtrc(:,:,jb,:) IF (ASSOCIATED(field%qtrc))
    !$ser data echam_mox_qtrc_mox=tend%qtrc_mox(:,:,jb,:) IF (ASSOCIATED(tend%qtrc_mox))
    !$ser data echam_mox_qtrc_phy=tend%qtrc_phy(:,:,jb,:) IF (ASSOCIATED(tend%qtrc_phy))
    !$ser verbatim lactive = .FALSE.
    !$ser verbatim IF (lonlyonce) THEN
    !$ser verbatim   lenabled = .FALSE.
    !$ser verbatim END IF
    !$ser verbatim ENDIF

  END SUBROUTINE serialize_mox_output

END MODULE mo_ser_echam_mox
