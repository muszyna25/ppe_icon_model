!--------------------------------------------------------------------
!
! Serialization routine for routine diagnose_pres_temp
!
!--------------------------------------------------------------------

MODULE mo_ser_diagnose_pres_temp

  USE mo_kind,               ONLY: vp, wp
  USE mo_exception,          ONLY: warning
  USE mo_ser_common,         ONLY: init
  USE mtime,                 ONLY: datetimeToString, MAX_DATETIME_STR_LEN
  USE mo_time_config,        ONLY: time_config
  USE mo_nonhydro_types,     ONLY: t_nh_prog, t_nh_diag
  USE mo_model_domain,       ONLY: t_patch
  IMPLICIT NONE

  PUBLIC :: serialize_input
  PUBLIC :: serialize_output

  CONTAINS

  SUBROUTINE serialize_input(jg, pt_prog, pt_prog_rcf, pt_diag, pt_patch)
    INTEGER, INTENT(IN) :: jg
    TYPE(t_nh_prog),    INTENT(INOUT) :: pt_prog      !!the prognostic variables
    TYPE(t_nh_prog),    INTENT(INOUT) :: pt_prog_rcf  !!the prognostic variables which are
                                                      !! treated with reduced calling frequency
    TYPE(t_nh_diag),    INTENT(INOUT) :: pt_diag      !!the diagnostic variables
    TYPE(t_patch),      INTENT(IN)    :: pt_patch     ! Patch
#if defined(SERIALIZE_ECHAM_ICONAM) || defined(SERIALIZE_ECHAM_ALL) || defined(SERIALIZE_ALL)

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: date
    LOGICAL, PARAMETER :: lonlyonce = .TRUE.
#if defined( SERIALIZE_CREATE_REFERENCE )
    LOGICAL, SAVE :: lenabled = .TRUE.
#else
    LOGICAL, SAVE :: lenabled = .FALSE.
#endif
    LOGICAL, SAVE :: lactive = .TRUE.

    !$ser verbatim IF (lenabled .and. lactive) THEN
    !$ser verbatim   CALL warning('SER:mo_ser_diagnose_pres_temp:input','Serialization is active!')
#if defined( _OPENACC )
    !$ser verbatim   CALL warning('GPU:mo_ser_diagnose_pres_temp:input','GPU HOST synchronization forced by serialization!')
    !$ser verbatim   !$ACC UPDATE HOST( pt_diag%tempv(:,:,:) ) IF( ASSOCIATED(pt_diag%tempv) )
    !$ser verbatim   !$ACC UPDATE HOST( pt_diag%temp(:,:,:) ) IF( ASSOCIATED(pt_diag%temp) )
    !$ser verbatim   !$ACC UPDATE HOST( pt_prog%exner(:,:,:) ) IF( ASSOCIATED(pt_prog%exner) )
    !$ser verbatim   !$ACC UPDATE HOST( pt_prog%theta_v(:,:,:) ) IF( ASSOCIATED(pt_prog%theta_v) )
    !$ser verbatim   !$ACC UPDATE HOST( pt_prog_rcf%tracer(:,:,:,:) ) IF( ASSOCIATED(pt_prog_rcf%tracer) )
#endif
    !$ser verbatim CALL datetimeToString(time_config%tc_current_date, date)
    !$ser verbatim CALL init('icon')
    !$ser savepoint echam_diagnose_pres_temp-input jg=jg date=TRIM(date)
#if defined SERIALIZE_CREATE_REFERENCE
    !$ser mode write
#elif defined SERIALIZE_PERTURB_REFERENCE
    !$ser mode read-perturb
#elif defined SERIALIZE_READ_REFERENCE
    !$ser mode read
#else
#error SERIALIZATION MODE IS NOT SET
#endif
    !$ser data nhs_ddyn_tempv=pt_diag%tempv(:,:,:)
    !$ser data nhs_ddyn_temp=pt_diag%temp(:,:,:)
    !$ser data nhs_ddyn_exner=pt_prog%exner(:,:,:)
    !$ser data nhs_ddyn_theta_v=pt_prog%theta_v(:,:,:)
    !$ser data nhs_ddyn_tracer_iqv=pt_prog_rcf%tracer(:,:,:,:)
    !$ser verbatim IF (lonlyonce) THEN
    !$ser verbatim   lactive = .FALSE.
    !$ser verbatim   lenabled = .FALSE.
    !$ser verbatim END IF
#if defined( SERIALIZE_READ_REFERENCE ) || defined( SERIALIZE_PERTURB_REFERENCE )
#if defined( _OPENACC )
    !$ser verbatim CALL warning('GPU:mo_ser_diagnose_pres_temp:input','GPU DEVICE synchronization forced by serialization!')
    !$ser verbatim   !$ACC UPDATE DEVICE( pt_diag%tempv(:,:,:) ) IF( ASSOCIATED(pt_diag%tempv) )
    !$ser verbatim   !$ACC UPDATE DEVICE( pt_diag%temp(:,:,:) ) IF( ASSOCIATED(pt_diag%temp) )
    !$ser verbatim   !$ACC UPDATE DEVICE( pt_prog%exner(:,:,:) ) IF( ASSOCIATED(pt_prog%exner) )
    !$ser verbatim   !$ACC UPDATE DEVICE( pt_prog%theta_v(:,:,:) ) IF( ASSOCIATED(pt_prog%theta_v) )
    !$ser verbatim   !$ACC UPDATE DEVICE( pt_prog_rcf%tracer(:,:,:,:) ) IF( ASSOCIATED(pt_prog_rcf%tracer) )
#endif
#endif
    !$ser verbatim ENDIF

#endif
  END SUBROUTINE serialize_input

  SUBROUTINE serialize_output(jg, pt_diag)
    INTEGER, INTENT(IN) :: jg
    TYPE(t_nh_diag),    INTENT(INOUT) :: pt_diag      !!the diagnostic variables
#if defined(SERIALIZE_ECHAM_ICONAM) || defined(SERIALIZE_ECHAM_ALL) || defined(SERIALIZE_ALL)

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: date
    LOGICAL, PARAMETER :: lonlyonce = .TRUE.
    LOGICAL, SAVE :: lenabled = .TRUE.
    LOGICAL, SAVE :: lactive = .TRUE.

    !$ser verbatim IF (lenabled .and. lactive) THEN
    !$ser verbatim   CALL warning('SER:mo_ser_diagnose_pres_temp:output','Serialization is active!')
#if defined( _OPENACC )
    !$ser verbatim   CALL warning('GPU:mo_ser_diagnose_pres_temp:output','GPU HOST synchronization forced by serialization!')
    !$ser verbatim   !$ACC UPDATE HOST( pt_diag%tempv(:,:,:) ) IF( ASSOCIATED(pt_diag%tempv) )
    !$ser verbatim   !$ACC UPDATE HOST( pt_diag%temp(:,:,:) ) IF( ASSOCIATED(pt_diag%temp) )
    !$ser verbatim   !$ACC UPDATE HOST( pt_diag%temp_ifc(:,:,:) ) IF( ASSOCIATED(pt_diag%temp_ifc) )
    !$ser verbatim   !$ACC UPDATE HOST( pt_diag%pres_ifc(:,:,:) ) IF( ASSOCIATED(pt_diag%pres_ifc) )
    !$ser verbatim   !$ACC UPDATE HOST( pt_diag%pres_sfc(:,:) ) IF( ASSOCIATED(pt_diag%pres_sfc) )
    !$ser verbatim   !$ACC UPDATE HOST( pt_diag%pres(:,:,:) ) IF( ASSOCIATED(pt_diag%pres) )
    !$ser verbatim   !$ACC UPDATE HOST( pt_diag%dpres_mc(:,:,:) ) IF( ASSOCIATED(pt_diag%dpres_mc) )
#endif
    !$ser verbatim CALL datetimeToString(time_config%tc_current_date, date)
    !$ser verbatim CALL init('icon')
    !$ser savepoint echam_diagnose_pres_temp-output jg=jg date=date
    !$ser mode write
    !$ser data nhs_ddyn_tempv=pt_diag%tempv(:,:,:)
    !$ser data nhs_ddyn_temp=pt_diag%temp(:,:,:)
    !$ser data nhs_ddyn_temp_ifc=pt_diag%temp_ifc(:,:,:)
    !$ser data nhs_ddyn_pres_ifc=pt_diag%pres_ifc(:,:,:)
    !$ser data nhs_ddyn_pres_sfc=pt_diag%pres_sfc(:,:)
    !$ser data nhs_ddyn_pres=pt_diag%pres(:,:,:)
    !$ser data nhs_ddyn_dpres_mc=pt_diag%dpres_mc(:,:,:)
    !$ser verbatim IF (lonlyonce) THEN
    !$ser verbatim   lactive = .FALSE.
    !$ser verbatim   lenabled = .FALSE.
    !$ser verbatim END IF
    !$ser verbatim END IF

#endif
  END SUBROUTINE serialize_output

END MODULE mo_ser_diagnose_pres_temp
