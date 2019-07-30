!--------------------------------------------------------------------
!
! Serialization routine for ECHAM WMO
!
!--------------------------------------------------------------------

MODULE mo_ser_echam_wmo

  USE mo_kind,               ONLY: vp, wp
  USE mo_exception,          ONLY: warning
  USE mo_ser_common,         ONLY: init
  USE mtime,                 ONLY: datetimeToString, MAX_DATETIME_STR_LEN
  USE mo_time_config,        ONLY: time_config
  USE mo_echam_phy_memory,   ONLY: t_echam_phy_field, t_echam_phy_tend
  IMPLICIT NONE

  INTEGER, PARAMETER :: selected_block = 1

  PUBLIC :: serialize_wmo_input
  PUBLIC :: serialize_wmo_output

  CONTAINS

  SUBROUTINE serialize_wmo_input(jg, jb, jcs, jce, nproma, nlev, field)
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
    !$ser verbatim   CALL warning('SER:mo_ser_echam_wmo:input','Serialization is active!')
#if defined( _OPENACC )
    !$ser verbatim   CALL warning('GPU:mo_ser_echam_wmo:input','GPU HOST synchronization forced by serialization!')
    !$ser verbatim   !$ACC UPDATE HOST( field%ta ) IF( ASSOCIATED(field%ta) )
    !$ser verbatim   !$ACC UPDATE HOST( field%presm_old ) IF( ASSOCIATED(field%presm_old) )
    !$ser verbatim   !$ACC UPDATE HOST( field%ptp ) IF( ASSOCIATED(field%ptp) )
#endif
    !$ser verbatim   CALL datetimeToString(time_config%tc_current_date, date)
    !$ser verbatim   CALL init('icon')
    !$ser savepoint echam_wmo-input jg=jg jb=jb jcs=jcs jce=jce nproma=nproma nlev=nlev date=TRIM(date)
#if defined( SERIALIZE_CREATE_REFERENCE )
    !$ser mode write
#elif defined( SERIALIZE_PERTURB_REFERENCE )
    !$ser mode read-perturb
#elif defined( SERIALIZE_READ_REFERENCE )
    !$ser mode read
#else
#error SERIALIZATION MODE IS NOT SET
#endif
    !$ser data echam_wmo_ta=field%ta(:,:,jb) IF (ASSOCIATED(field%ta))
    !$ser data echam_wmo_presm_old=field%presm_old(:,:,jb) IF (ASSOCIATED(field%presm_old))
    !$ser data echam_wmo_ptp=field%ptp(:,jb) IF (ASSOCIATED(field%ptp))
    !$ser verbatim lactive = .FALSE.
    !$ser verbatim IF (lonlyonce) THEN
    !$ser verbatim   lenabled = .FALSE.
    !$ser verbatim END IF
#if defined( SERIALIZE_READ_REFERENCE ) || defined( SERIALIZE_PERTURB_REFERENCE )
#if defined( _OPENACC )
    !$ser verbatim CALL warning('GPU:mo_ser_echam_wmo:input','GPU DEVICE synchronization forced by serialization!')
    !$ser verbatim   !$ACC UPDATE DEVICE( field%ta ) IF( ASSOCIATED(field%ta) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%presm_old ) IF( ASSOCIATED(field%presm_old) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%ptp ) IF( ASSOCIATED(field%ptp) )
#endif
#endif
    !$ser verbatim ENDIF

  END SUBROUTINE serialize_wmo_input

  SUBROUTINE serialize_wmo_output(jg, jb, jcs, jce, nproma, nlev, field)
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
    !$ser verbatim   CALL warning('SER:mo_ser_echam_wmo:output','Serialization is active!')
#if defined( _OPENACC )
    !$ser verbatim   CALL warning('GPU:mo_ser_echam_wmo:output','GPU HOST synchronization forced by serialization!')
    !$ser verbatim   !$ACC UPDATE HOST( field%ta ) IF( ASSOCIATED(field%ta) )
    !$ser verbatim   !$ACC UPDATE HOST( field%presm_old ) IF( ASSOCIATED(field%presm_old) )
    !$ser verbatim   !$ACC UPDATE HOST( field%ptp ) IF( ASSOCIATED(field%ptp) )
#endif
    !$ser verbatim   CALL datetimeToString(time_config%tc_current_date, date)
    !$ser verbatim   CALL init('icon')
    !$ser savepoint echam_wmo-output jg=jg jb=jb jcs=jcs jce=jce nproma=nproma nlev=nlev date=TRIM(date) 
    !$ser mode write
    !$ser data echam_wmo_ta=field%ta(:,:,jb) IF (ASSOCIATED(field%ta))
    !$ser data echam_wmo_presm_old=field%presm_old(:,:,jb) IF (ASSOCIATED(field%presm_old))
    !$ser data echam_wmo_ptp=field%ptp(:,jb) IF (ASSOCIATED(field%ptp))
    !$ser verbatim lactive = .FALSE.
    !$ser verbatim IF (lonlyonce) THEN
    !$ser verbatim   lenabled = .FALSE.
    !$ser verbatim END IF
    !$ser verbatim ENDIF

  END SUBROUTINE serialize_wmo_output

END MODULE mo_ser_echam_wmo
