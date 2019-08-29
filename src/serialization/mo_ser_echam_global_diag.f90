!--------------------------------------------------------------------
!
! Serialization routine for ECHAM global_diagnostic
!
!--------------------------------------------------------------------

MODULE mo_ser_echam_global_diag

  USE mo_kind,               ONLY: vp, wp
  USE mo_exception,          ONLY: warning
  USE mo_ser_common,         ONLY: init
  USE mtime,                 ONLY: datetimeToString, MAX_DATETIME_STR_LEN
  USE mo_time_config,        ONLY: time_config
  USE mo_echam_phy_memory,   ONLY: t_echam_phy_field, t_echam_phy_tend
  IMPLICIT NONE

  PUBLIC :: serialize_input
  PUBLIC :: serialize_output

  CONTAINS

  SUBROUTINE serialize_input(jg, field)
    INTEGER, INTENT(IN)                     :: jg
    TYPE(t_echam_phy_field), INTENT(INOUT)  :: field
#if defined(SERIALIZE_ECHAM_GLOBAL_DIAG) || defined(SERIALIZE_ECHAM_ALL) || defined(SERIALIZE_ALL)

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: date
    LOGICAL, PARAMETER :: lonlyonce = .TRUE.
#if defined( SERIALIZE_CREATE_REFERENCE )
    LOGICAL, SAVE :: lenabled = .TRUE.
#else
    LOGICAL, SAVE :: lenabled = .FALSE.
#endif
    LOGICAL, SAVE :: lactive = .TRUE.

    !$ser verbatim IF (lenabled .and. lactive) THEN
    !$ser verbatim   CALL warning('SER:mo_ser_echam_global_diag:input','Serialization is active!')
#if defined( _OPENACC )
    !$ser verbatim   CALL warning('GPU:mo_ser_echam_global_diag:input','GPU HOST synchronization forced by serialization!')
    !$ser verbatim   !$ACC UPDATE HOST( field%tas ) IF( ASSOCIATED(field%tas) )
    !$ser verbatim   !$ACC UPDATE HOST( field%tas_gmean ) IF( ASSOCIATED(field%tas_gmean) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rsdt ) IF( ASSOCIATED(field%rsdt) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rsdt_gmean ) IF( ASSOCIATED(field%rsdt_gmean) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rsut ) IF( ASSOCIATED(field%rsut) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rsut_gmean ) IF( ASSOCIATED(field%rsut_gmean) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rlut ) IF( ASSOCIATED(field%rlut) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rlut_gmean ) IF( ASSOCIATED(field%rlut_gmean) )
    !$ser verbatim   !$ACC UPDATE HOST( field%pr ) IF( ASSOCIATED(field%pr) )
    !$ser verbatim   !$ACC UPDATE HOST( field%prec_gmean ) IF( ASSOCIATED(field%prec_gmean) )
    !$ser verbatim   !$ACC UPDATE HOST( field%evap ) IF( ASSOCIATED(field%evap) )
    !$ser verbatim   !$ACC UPDATE HOST( field%evap_gmean ) IF( ASSOCIATED(field%evap_gmean) )
    !$ser verbatim   !$ACC UPDATE HOST( field%radtop_gmean ) IF( ASSOCIATED(field%radtop_gmean) )
    !$ser verbatim   !$ACC UPDATE HOST( field%sftof ) IF( ASSOCIATED(field%sftof) )
    !$ser verbatim   !$ACC UPDATE HOST( field%fwfoce_gmean ) IF( ASSOCIATED(field%fwfoce_gmean) )
#endif
    !$ser verbatim   CALL datetimeToString(time_config%tc_current_date, date)
    !$ser verbatim   CALL init('icon')
    !$ser savepoint echam_global_diag-input jg=jg date=TRIM(date)
#if defined( SERIALIZE_CREATE_REFERENCE )
    !$ser mode write
#elif defined( SERIALIZE_PERTURB_REFERENCE )
    !$ser mode read-perturb
#elif defined( SERIALIZE_READ_REFERENCE )
    !$ser mode read
#else
#error SERIALIZATION MODE IS NOT SET
#endif
    !$ser data echam_global_diag_tas=field%tas IF (ASSOCIATED(field%tas))
    !$ser data echam_global_diag_tas_gmean=field%tas_gmean IF (ASSOCIATED(field%tas_gmean))
    !$ser data echam_global_diag_rsdt=field%rsdt IF (ASSOCIATED(field%rsdt))
    !$ser data echam_global_diag_rsdt_gmean=field%rsdt_gmean IF (ASSOCIATED(field%rsdt_gmean))
    !$ser data echam_global_diag_rsut=field%rsut IF (ASSOCIATED(field%rsut))
    !$ser data echam_global_diag_rsut_gmean=field%rsut_gmean IF (ASSOCIATED(field%rsut_gmean))
    !$ser data echam_global_diag_rlut=field%rlut IF (ASSOCIATED(field%rlut))
    !$ser data echam_global_diag_rlut_gmean=field%rlut_gmean IF (ASSOCIATED(field%rlut_gmean))
    !$ser data echam_global_diag_pr=field%pr IF (ASSOCIATED(field%pr))
    !$ser data echam_global_diag_prec_gmean=field%prec_gmean IF (ASSOCIATED(field%prec_gmean))
    !$ser data echam_global_diag_evap=field%evap IF (ASSOCIATED(field%evap))
    !$ser data echam_global_diag_evap_gmean=field%evap_gmean IF (ASSOCIATED(field%evap_gmean))
    !$ser data echam_global_diag_radtop_gmean=field%radtop_gmean IF (ASSOCIATED(field%radtop_gmean))
    !$ser data echam_global_diag_sftof=field%sftof IF (ASSOCIATED(field%sftof))
    !$ser data echam_global_diag_fwfoce_gmean=field%fwfoce_gmean IF (ASSOCIATED(field%fwfoce_gmean))
    !$ser verbatim IF (lonlyonce) THEN
    !$ser verbatim   lactive = .FALSE.
    !$ser verbatim   lenabled = .FALSE.
    !$ser verbatim END IF
#if defined( SERIALIZE_READ_REFERENCE ) || defined( SERIALIZE_PERTURB_REFERENCE )
#if defined( _OPENACC )
    !$ser verbatim CALL warning('GPU:mo_ser_echam_global_diag:input','GPU DEVICE synchronization forced by serialization!')
    !$ser verbatim   !$ACC UPDATE DEVICE( field%tas ) IF( ASSOCIATED(field%tas) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%tas_gmean ) IF( ASSOCIATED(field%tas_gmean) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%rsdt ) IF( ASSOCIATED(field%rsdt) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%rsdt_gmean ) IF( ASSOCIATED(field%rsdt_gmean) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%rsut ) IF( ASSOCIATED(field%rsut) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%rsut_gmean ) IF( ASSOCIATED(field%rsut_gmean) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%rlut ) IF( ASSOCIATED(field%rlut) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%rlut_gmean ) IF( ASSOCIATED(field%rlut_gmean) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%pr ) IF( ASSOCIATED(field%pr) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%prec_gmean ) IF( ASSOCIATED(field%prec_gmean) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%evap ) IF( ASSOCIATED(field%evap) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%evap_gmean ) IF( ASSOCIATED(field%evap_gmean) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%radtop_gmean ) IF( ASSOCIATED(field%radtop_gmean) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%sftof ) IF( ASSOCIATED(field%sftof) )
    !$ser verbatim   !$ACC UPDATE DEVICE( field%fwfoce_gmean ) IF( ASSOCIATED(field%fwfoce_gmean) )
#endif
#endif
    !$ser verbatim ENDIF

#endif
  END SUBROUTINE serialize_input

  SUBROUTINE serialize_output(jg, field)
    INTEGER, INTENT(IN)                     :: jg
    TYPE(t_echam_phy_field), INTENT(INOUT)  :: field
#if defined(SERIALIZE_ECHAM_GLOBAL_DIAG) || defined(SERIALIZE_ECHAM_ALL) || defined(SERIALIZE_ALL)

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: date
    LOGICAL, PARAMETER :: lonlyonce = .TRUE.
    LOGICAL, SAVE :: lenabled = .TRUE.
    LOGICAL, SAVE :: lactive = .TRUE.

    !$ser verbatim IF (lenabled .and. lactive) THEN
    !$ser verbatim   CALL warning('SER:mo_ser_echam_global_diag:output','Serialization is active!')
#if defined( _OPENACC )
    !$ser verbatim   CALL warning('GPU:mo_ser_echam_global_diag:output','GPU HOST synchronization forced by serialization!')
    !$ser verbatim   !$ACC UPDATE HOST( field%tas ) IF( ASSOCIATED(field%tas) )
    !$ser verbatim   !$ACC UPDATE HOST( field%tas_gmean ) IF( ASSOCIATED(field%tas_gmean) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rsdt ) IF( ASSOCIATED(field%rsdt) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rsdt_gmean ) IF( ASSOCIATED(field%rsdt_gmean) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rsut ) IF( ASSOCIATED(field%rsut) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rsut_gmean ) IF( ASSOCIATED(field%rsut_gmean) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rlut ) IF( ASSOCIATED(field%rlut) )
    !$ser verbatim   !$ACC UPDATE HOST( field%rlut_gmean ) IF( ASSOCIATED(field%rlut_gmean) )
    !$ser verbatim   !$ACC UPDATE HOST( field%pr ) IF( ASSOCIATED(field%pr) )
    !$ser verbatim   !$ACC UPDATE HOST( field%prec_gmean ) IF( ASSOCIATED(field%prec_gmean) )
    !$ser verbatim   !$ACC UPDATE HOST( field%evap ) IF( ASSOCIATED(field%evap) )
    !$ser verbatim   !$ACC UPDATE HOST( field%evap_gmean ) IF( ASSOCIATED(field%evap_gmean) )
    !$ser verbatim   !$ACC UPDATE HOST( field%radtop_gmean ) IF( ASSOCIATED(field%radtop_gmean) )
    !$ser verbatim   !$ACC UPDATE HOST( field%sftof ) IF( ASSOCIATED(field%sftof) )
    !$ser verbatim   !$ACC UPDATE HOST( field%fwfoce_gmean ) IF( ASSOCIATED(field%fwfoce_gmean) )
#endif
    !$ser verbatim   CALL datetimeToString(time_config%tc_current_date, date)
    !$ser verbatim   CALL init('icon')
    !$ser savepoint echam_global_diag-output jg=jg date=TRIM(date)
    !$ser mode write
    !$ser data echam_global_diag_tas=field%tas IF (ASSOCIATED(field%tas))
    !$ser data echam_global_diag_tas_gmean=field%tas_gmean IF (ASSOCIATED(field%tas_gmean))
    !$ser data echam_global_diag_rsdt=field%rsdt IF (ASSOCIATED(field%rsdt))
    !$ser data echam_global_diag_rsdt_gmean=field%rsdt_gmean IF (ASSOCIATED(field%rsdt_gmean))
    !$ser data echam_global_diag_rsut=field%rsut IF (ASSOCIATED(field%rsut))
    !$ser data echam_global_diag_rsut_gmean=field%rsut_gmean IF (ASSOCIATED(field%rsut_gmean))
    !$ser data echam_global_diag_rlut=field%rlut IF (ASSOCIATED(field%rlut))
    !$ser data echam_global_diag_rlut_gmean=field%rlut_gmean IF (ASSOCIATED(field%rlut_gmean))
    !$ser data echam_global_diag_pr=field%pr IF (ASSOCIATED(field%pr))
    !$ser data echam_global_diag_prec_gmean=field%prec_gmean IF (ASSOCIATED(field%prec_gmean))
    !$ser data echam_global_diag_evap=field%evap IF (ASSOCIATED(field%evap))
    !$ser data echam_global_diag_evap_gmean=field%evap_gmean IF (ASSOCIATED(field%evap_gmean))
    !$ser data echam_global_diag_radtop_gmean=field%radtop_gmean IF (ASSOCIATED(field%radtop_gmean))
    !$ser data echam_global_diag_sftof=field%sftof IF (ASSOCIATED(field%sftof))
    !$ser data echam_global_diag_fwfoce_gmean=field%fwfoce_gmean IF (ASSOCIATED(field%fwfoce_gmean))
    !$ser verbatim IF (lonlyonce) THEN
    !$ser verbatim   lactive = .FALSE.
    !$ser verbatim   lenabled = .FALSE.
    !$ser verbatim END IF
    !$ser verbatim ENDIF

#endif
  END SUBROUTINE serialize_output

END MODULE mo_ser_echam_global_diag
