!--------------------------------------------------------------------
!
! Serialization routine for ECHAM vdiff down
!
!--------------------------------------------------------------------

MODULE mo_ser_echam_vdiff_down

  USE mo_kind,               ONLY: vp, wp
  USE mo_ser_common,         ONLY: init
  USE mtime,                 ONLY: datetimeToString, MAX_DATETIME_STR_LEN
  USE mo_time_config,        ONLY: time_config
  USE mo_echam_phy_memory,   ONLY: t_echam_phy_field
  USE mo_run_config,         ONLY: iqv, iqc, iqi, iqt
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

  SUBROUTINE serialize_input(jg, jb, jcs, kproma, kbdim, klev, klevm1, klevp1, ktrac, &
                             ksfc_type, idx_wtr, idx_ice, idx_lnd, pdtime,            &
                             field, pxm1, pxt_emis, pthvvar, pxvar, pwstar,           &
                             pqsat_tile, phdtcbl, pri, pri_tile, pmixlen, pcfm,       &
                             pcfm_tile, pcfh, pcfh_tile, pcfv, pcftotte, pcfthv,      &
                             aa, aa_btm, bb, bb_btm, pfactor_sfc, pcpt_tile,          &
                             pcptgz, pzthvvar, pztottevn)
    INTEGER, INTENT(IN)    :: jg, jb, jcs, kproma
    INTEGER, INTENT(IN)    :: kbdim, klev, klevm1, klevp1, ktrac
    INTEGER, INTENT(IN)    :: ksfc_type, idx_wtr, idx_ice, idx_lnd
    REAL(wp),INTENT(IN)    :: pdtime
    TYPE(t_echam_phy_field), POINTER, INTENT(INOUT) :: field
    REAL(wp),INTENT(INOUT) :: &
      & pxm1(:,:),         &
      & pxt_emis(:,:),     &
      & pthvvar(:,:),      &
      & pxvar(:,:),        &
      & pwstar(:),         &
      & pqsat_tile(:,:),   &
      & phdtcbl(:),        &
      & pri(:,:),          &
      & pri_tile(:,:),     &
      & pmixlen(:,:),      &
      & pcfm(:,:),         &
      & pcfm_tile(:,:),    &
      & pcfh(:,:),         &
      & pcfh_tile(:,:),    &
      & pcfv(:,:),         &
      & pcftotte (:,:),    &
      & pcfthv (:,:),      &
      & aa(:,:,:,:),       &
      & aa_btm(:,:,:,:),   &
      & bb(:,:,:),         &
      & bb_btm(:,:,:),     &
      & pfactor_sfc(:),    &
      & pcpt_tile(:,:),    &
      & pcptgz(:,:),       &
      & pzthvvar(:,:),     &
      & pztottevn(:,:)

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: date

    !$ser verbatim IF (singleBlock < 0 .OR. jb == singleBlock) THEN
    !$ser verbatim   writeIn = .TRUE.
    !$ser verbatim ENDIF
    !$ser verbatim IF (serializeStepIn .and. writeIn) THEN
    !$ser verbatim   CALL datetimeToString(time_config%tc_current_date, date)
    !$ser verbatim   CALL init('echam_vdiff_down')
    !$ser savepoint echam_vdiff_down-input jg=jg jb=jb jcs=jcs kproma=kproma &
    !$ser&          kbdim=kbdim klev=klev klevm1=klevm1 klevp1=klevp1 &
    !$ser&          ktrac=ktrac ksfc_type=ksfc_type idx_wtr=idx_wtr &
    !$ser&          idx_ice=idx_ice idx_lnd=idx_lnd iqv=iqv iqc=iqc iqi=iqi &
    !$ser&         iqt=iqt pdtime=pdtime date=TRIM(date)
#if defined SERIALIZE_CREATE_REFERENCE
    !$ser mode write
#elif defined SERIALIZE_PERTURB_REFERENCE
    !$ser mode read-perturb
#elif defined SERIALIZE_READ_REFERENCE
    !$ser mode read
#else
#error SERIALIZATION MODE IS NOT SET
#endif
    !$ser data pcoriol=field%coriol(:,jb)           &
    !$ser&     pzf=field%zf(:,:,jb)                 &
    !$ser&     pzh=field%zh(:,:,jb)                 &
    !$ser&     pfrc=field%frac_tile(:,jb,:)         &
    !$ser&     ptsfc_tile=field%ts_tile(:,jb,:)     &
    !$ser&     pocu=field%ocu(:,jb)                 &
    !$ser&     pocv=field%ocv(:,jb)                 &
    !$ser&     ppsfc=field%presi_old(:,klevp1,jb)   &
    !$ser&     pum1=field%ua(:,:,jb)                &
    !$ser&     pvm1=field%va(:,:,jb)                &
    !$ser&     ptm1=field%ta(:,:,jb)                &
    !$ser&     pqm1=field%qtrc(:,:,jb,iqv)          &
    !$ser&     pxlm1=field%qtrc(:,:,jb,iqc)         &
    !$ser&     pxim1=field%qtrc(:,:,jb,iqi)         &
    !$ser&     pxm1=pxm1                            &
    !$ser&     pxtm1=field%qtrc(:,:,jb,iqt:)        &
    !$ser&     pmair=field%mair(:,:,jb)             &
    !$ser&     pmref=field%mref(:,:,jb)             &
    !$ser&     paphm1=field%presi_old(:,:,jb)       &
    !$ser&     papm1=field%presm_old(:,:,jb)        &
    !$ser&     ptvm1=field%tv(:,:,jb)               &
    !$ser&     paclc=field%aclc(:,:,jb)             &
    !$ser&     pxt_emis=pxt_emis                    &
    !$ser&     pthvvar=pthvvar                      &
    !$ser&     pxvar=pxvar                          &
    !$ser&     pz0m_tile=field%z0m_tile(:,jb,:)     &
    !$ser&     ptottem1=field%tottem1(:,:,jb)       &
    !$ser&     pustar=field%ustar(:,jb)             &
    !$ser&     pwstar=pwstar                        &
    !$ser&     pwstar_tile=field%wstar_tile(:,jb,:) &
    !$ser&     pqsat_tile=pqsat_tile                &
    !$ser&     phdtcbl=phdtcbl                      &
    !$ser&     pri=pri                              &
    !$ser&     pri_tile=pri_tile                    &
    !$ser&     pmixlen=pmixlen                      &
    !$ser&     pcfm=pcfm                            &
    !$ser&     pcfm_tile=pcfm_tile                  &
    !$ser&     pcfh=pcfh                            &
    !$ser&     pcfh_tile=pcfh_tile                  &
    !$ser&     pcfv=pcfv                            &
    !$ser&     pcftotte=pcftotte                    &
    !$ser&     pcfthv=pcfthv                        &
    !$ser&     aa=aa                                &
    !$ser&     aa_btm=aa_btm                        &
    !$ser&     bb=bb                                &
    !$ser&     bb_btm=bb_btm                        &
    !$ser&     pfactor_sfc=pfactor_sfc              &
    !$ser&     pcpt_tile=pcpt_tile                  &
    !$ser&     pcptgz=pcptgz                        &
    !$ser&     pzthvvar=pzthvvar                    &
    !$ser&     pthvsig=field%thvsig(:,jb)           &
    !$ser&     pztottevn=pztottevn                  &
    !$ser&     pcsat=field%csat(:,jb)               &
    !$ser&     pcair=field%cair(:,jb)               &
    !$ser&     paz0lh=field%z0h_lnd(:,jb)
    !$ser verbatim writeIn = .FALSE.
    !$ser verbatim IF (singleStepIn) THEN
    !$ser verbatim   serializeStepIn = .FALSE.
    !$ser verbatim END IF
    !$ser verbatim ENDIF

  END SUBROUTINE serialize_input

  SUBROUTINE serialize_output(jg, jb, jcs, kproma, kbdim, klev, klevm1, klevp1, ktrac, &
                             ksfc_type, idx_wtr, idx_ice, idx_lnd, pdtime,             &
                             field, pwstar, pqsat_tile, phdtcbl, pri, pri_tile,        &
                             pmixlen, pcfm, pcfm_tile, pcfh, pcfh_tile, pcfv,          &
                             pcftotte, pcfthv, aa, aa_btm, bb, bb_btm,                 &
                             pfactor_sfc, pcpt_tile, pcptgz, pzthvvar,                 &
                             pztottevn, pch_tile, pbn_tile, pbhn_tile,                 &
                             pbm_tile, pbh_tile)
    INTEGER, INTENT(IN)     :: jg, jb, jcs, kproma
    INTEGER, INTENT(IN)     :: kbdim, klev, klevm1, klevp1, ktrac
    INTEGER, INTENT(IN)     :: ksfc_type, idx_wtr, idx_ice, idx_lnd
    REAL(wp), INTENT(IN)    :: pdtime
    TYPE(t_echam_phy_field) ,POINTER, INTENT(INOUT) :: field
    REAL(wp), INTENT(INOUT) :: &
      pwstar(:),            &
      pqsat_tile(:,:),      &
      phdtcbl(:),           &
      pri(:,:),             &
      pri_tile(:,:),        &
      pmixlen(:,:),         &
      pcfm(:,:),            &
      pcfm_tile(:,:),       &
      pcfh(:,:),            &
      pcfh_tile(:,:),       &
      pcfv(:,:),            &
      pcftotte (:,:),       &
      pcfthv (:,:),         &
      aa(:,:,:,:),          &
      aa_btm(:,:,:,:),      &
      bb(:,:,:),            &
      bb_btm(:,:,:),        &
      pfactor_sfc(:),       &
      pcpt_tile(:,:),       &
      pcptgz(:,:),          &
      pzthvvar(:,:),        &
      pztottevn(:,:),       &
      pch_tile(:,:),        &
      pbn_tile(:,:),        &
      pbhn_tile(:,:),       &
      pbm_tile(:,:),        &
      pbh_tile(:,:)


    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: date

    !$ser verbatim IF (singleBlock < 0 .OR. jb == singleBlock) THEN
    !$ser verbatim   writeOut = .TRUE.
    !$ser verbatim ENDIF
    !$ser verbatim IF (serializeStepOut .and. writeOut) THEN
    !$ser verbatim   CALL datetimeToString(time_config%tc_current_date, date)
    !$ser verbatim   CALL init('echam_vdiff_down')
    !$ser savepoint echam_vdiff_down-output jg=jg jb=jb jcs=jcs kproma=kproma &
    !$ser&          kbdim=kbdim klev=klev klevm1=klevm1 klevp1=klevp1 &
    !$ser&          ktrac=ktrac ksfc_type=ksfc_type idx_wtr=idx_wtr &
    !$ser&          idx_ice=idx_ice idx_lnd=idx_lnd iqv=iqv iqc=iqc iqi=iqi &
    !$ser&          iqt=iqt pdtime=pdtime date=TRIM(date)
    !$ser mode write
    !$ser data pustar=field%ustar(:,jb)             &
    !$ser&     pwstar=pwstar                        &
    !$ser&     pwstar_tile=field%wstar_tile(:,jb,:) &
    !$ser&     pqsat_tile=pqsat_tile                &
    !$ser&     phdtcbl=phdtcbl                      &
    !$ser&     pri=pri                              &
    !$ser&     pri_tile=pri_tile                    &
    !$ser&     pmixlen=pmixlen                      &
    !$ser&     pcfm=pcfm                            &
    !$ser&     pcfm_tile=pcfm_tile                  &
    !$ser&     pcfh=pcfh                            &
    !$ser&     pcfh_tile=pcfh_tile                  &
    !$ser&     pcfv=pcfv                            &
    !$ser&     pcftotte=pcftotte                    &
    !$ser&     pcfthv=pcfthv                        &
    !$ser&     aa=aa                                &
    !$ser&     aa_btm=aa_btm                        &
    !$ser&     bb=bb                                &
    !$ser&     bb_btm=bb_btm                        &
    !$ser&     pfactor_sfc=pfactor_sfc              &
    !$ser&     pcpt_tile=pcpt_tile                  &
    !$ser&     pcptgz=pcptgz                        &
    !$ser&     pzthvvar=pzthvvar                    &
    !$ser&     pthvsig=field%thvsig(:,jb)           &
    !$ser&     pztottevn=pztottevn                  &
    !$ser&     pch_tile=pch_tile                    &
    !$ser&     pbn_tile=pbn_tile                    &
    !$ser&     pbhn_tile=pbhn_tile                  &
    !$ser&     pbm_tile=pbm_tile                    &
    !$ser&     pbh_tile=pbh_tile
    !$ser verbatim writeOut = .FALSE.
    !$ser verbatim IF (singleStepOut) THEN
    !$ser verbatim   serializeStepOut = .FALSE.
    !$ser verbatim END IF
    !$ser verbatim ENDIF

  END SUBROUTINE serialize_output

END MODULE mo_ser_echam_vdiff_down
