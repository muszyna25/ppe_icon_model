!--------------------------------------------------------------------
!
! Serialization routine for ECHAM nsurf_diag
!
!--------------------------------------------------------------------

MODULE mo_ser_echam_nsurf_diag

  USE mo_kind,               ONLY: vp, wp
  USE mo_ser_common,         ONLY: init
  USE mtime,                 ONLY: datetimeToString, MAX_DATETIME_STR_LEN
  USE mo_time_config,        ONLY: time_config
  USE mo_echam_phy_memory,   ONLY: t_echam_phy_field
  USE mo_run_config,         ONLY: nlev, nlevp1, iqv
  IMPLICIT NONE

  LOGICAL :: writeIn = .FALSE.
  LOGICAL :: writeOut = .FALSE.
  LOGICAL :: serializeStepIn = .TRUE.
  LOGICAL :: serializeStepOut = .TRUE.
  LOGICAL, PARAMETER :: singleStepIn = .TRUE.
  LOGICAL, PARAMETER :: singleStepOut = .TRUE.
  INTEGER, PARAMETER :: singleBlock = 7

  PUBLIC :: serialize_input
  PUBLIC :: serialize_output

  CONTAINS

  SUBROUTINE serialize_input(jb, kproma, kbdim, ksfc_type, idx_lnd, field,     &
                             pxm1_col, pcptgz_col, pcpt_tile, pbn_tile, pbhn_tile,     &
                             pbh_tile, pbm_tile, pri_tile)
    INTEGER, INTENT(IN)    :: jb, kproma
    INTEGER, INTENT(INOUT) :: kbdim, ksfc_type, idx_lnd
    TYPE(t_echam_phy_field),POINTER, INTENT(INOUT) :: field
    REAL(wp),INTENT(INOUT) :: &
      pxm1_col(:),            &
      pcptgz_col(:),          &
      pcpt_tile(:,:),         &
      pbn_tile(:,:),          &
      pbhn_tile(:,:),         &
      pbh_tile(:,:),          &
      pbm_tile(:,:),          &
      pri_tile(:,:)

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: date


    !$ser verbatim IF (singleBlock < 0 .OR. jb == singleBlock) THEN
    !$ser verbatim   writeIn = .TRUE.
    !$ser verbatim ENDIF
    !$ser verbatim IF (serializeStepIn .and. writeIn) THEN
    !$ser verbatim   CALL datetimeToString(time_config%tc_current_date, date)
    !$ser verbatim   CALL init('echam_nsurf_diag')
    !$ser savepoint echam_nsurf_diag-input jb=jb kproma=kproma nlev=nlev nlevp1=nlevp1 iqv=iqv &
    !$ser&          date=TRIM(date)
#if defined SERIALIZE_CREATE_REFERENCE 
    !$ser mode write
#elif defined SERIALIZE_PERTURB_REFERENCE
    !$ser mode read-perturb
#elif defined SERIALIZE_READ_REFERENCE
    !$ser mode read
#else
#error SERIALIZATION MODE IS NOT SET
#endif 
    !$ser data kbdim=kbdim                                  &
    !$ser&     ksfc_type=ksfc_type                          &
    !$ser&     idx_lnd=idx_lnd                              &
    !$ser&     pfrc=field%frac_tile(:,jb,:)                 &
    !$ser&     pqm1_col=field%qtrc(:,nlev,jb,iqv)           &
    !$ser&     ptm1_col=field%ta(:,nlev,jb)                 &
    !$ser&     papm1_col=field%presm_old(:,nlev,jb)         &
    !$ser&     paphm1_col=field%presi_old(:,nlevp1,jb)      &
    !$ser&     pxm1_col=pxm1_col                            &
    !$ser&     pum1_col=field%ua(:,nlev,jb)                 &
    !$ser&     pvm1_col=field%va(:,nlev,jb)                 &
    !$ser&     pocu=field%ocu(:,jb)                         &
    !$ser&     pocv=field%ocv(:,jb)                         &
    !$ser&     pzf_col=field%zf(:,nlev,jb)                  &
    !$ser&     pzs_col=field%zh(:,nlev+1,jb)                &
    !$ser&     pcptgz_col=pcptgz_col                        &
    !$ser&     pcpt_tile=pcpt_tile                          &
    !$ser&     pbn_tile=pbn_tile                            &
    !$ser&     pbhn_tile=pbhn_tile                          &
    !$ser&     pbh_tile=pbh_tile                            &
    !$ser&     pbm_tile=pbm_tile                            &
    !$ser&     pri_tile=pri_tile                            &
    !$ser&     ptasmax=field%tasmax(:,jb)                   &
    !$ser&     ptasmin=field%tasmin(:,jb)
    !$ser verbatim writeIn = .FALSE.
    !$ser verbatim IF (singleStepIn) THEN
    !$ser verbatim   serializeStepIn = .FALSE.
    !$ser verbatim END IF
    !$ser verbatim ENDIF

  END SUBROUTINE serialize_input

  SUBROUTINE serialize_output(jb, kproma, kbdim, ksfc_type, idx_lnd, field)
    INTEGER, INTENT(IN)    :: jb, kproma
    INTEGER, INTENT(INOUT) :: kbdim, ksfc_type, idx_lnd
    TYPE(t_echam_phy_field), POINTER, INTENT(INOUT) :: field

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: date

    !$ser verbatim IF (singleBlock < 0 .OR. jb == singleBlock) THEN
    !$ser verbatim   writeOut = .TRUE.
    !$ser verbatim ENDIF
    !$ser verbatim IF (serializeStepOut .and. writeOut) THEN
    !$ser verbatim   CALL datetimeToString(time_config%tc_current_date, date)
    !$ser verbatim   CALL init('echam_nsurf_diag')
    !$ser savepoint echam_nsurf_diag-output jb=jb kproma=kproma date=TRIM(date)
    !$ser mode write
    !$ser data kbdim=kbdim                                  &
    !$ser&     ksfc_type=ksfc_type                          &
    !$ser&     idx_lnd=idx_lnd                              &
    !$ser&     psfcWind_gbm=field%sfcWind(:,jb)             &
    !$ser&     ptas_gbm=field%tas(:,jb)                     &
    !$ser&     pdew2_gbm=field%dew2(:,jb)                   &
    !$ser&     puas_gbm=field%uas(:,jb)                     &
    !$ser&     pvas_gbm=field%vas(:,jb)                     &
    !$ser&     ptasmax=field%tasmax(:,jb)                   &
    !$ser&     ptasmin=field%tasmin(:,jb)                   &
    !$ser&     psfcWind_tile=field%sfcWind_tile(:,jb,:)     &
    !$ser&     ptas_tile=field%tas_tile(:,jb,:)             &
    !$ser&     pdew2_tile=field%dew2_tile(:,jb,:)           &
    !$ser&     puas_tile=field%uas_tile(:,jb,:)             &
    !$ser&     pvas_tile=field%vas_tile(:,jb,:)
    !$ser verbatim writeOut = .FALSE.
    !$ser verbatim IF (singleStepOut) THEN
    !$ser verbatim   serializeStepOut = .FALSE.
    !$ser verbatim END IF
    !$ser verbatim ENDIF

  END SUBROUTINE serialize_output

END MODULE mo_ser_echam_nsurf_diag