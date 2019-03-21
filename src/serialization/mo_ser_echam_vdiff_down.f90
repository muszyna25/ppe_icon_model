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
    !$ser verbatim   CALL init('icon')
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
    !$ser data echam_vd_pcoriol=field%coriol(:,jb)           &
    !$ser&     echam_vd_pzf=field%zf(:,:,jb)                 &
    !$ser&     echam_vd_pzh=field%zh(:,:,jb)                 &
    !$ser&     echam_vd_pfrc=field%frac_tile(:,jb,:)         &
    !$ser&     echam_vd_ptsfc_tile=field%ts_tile(:,jb,:)     &
    !$ser&     echam_vd_pocu=field%ocu(:,jb)                 &
    !$ser&     echam_vd_pocv=field%ocv(:,jb)                 &
    !$ser&     echam_vd_ppsfc=field%presi_old(:,klevp1,jb)   &
    !$ser&     echam_vd_pum1=field%ua(:,:,jb)                &
    !$ser&     echam_vd_pvm1=field%va(:,:,jb)                &
    !$ser&     echam_vd_ptm1=field%ta(:,:,jb)                &
    !$ser&     echam_vd_pqm1=field%qtrc(:,:,jb,iqv)          &
    !$ser&     echam_vd_pxlm1=field%qtrc(:,:,jb,iqc)         &
    !$ser&     echam_vd_pxim1=field%qtrc(:,:,jb,iqi)         &
    !$ser&     echam_vd_pxm1=pxm1                            &
    !$ser&     echam_vd_pxtm1=field%qtrc(:,:,jb,iqt:)        &
    !$ser&     echam_vd_pmair=field%mair(:,:,jb)             &
    !$ser&     echam_vd_pmref=field%mref(:,:,jb)             &
    !$ser&     echam_vd_paphm1=field%presi_old(:,:,jb)       &
    !$ser&     echam_vd_papm1=field%presm_old(:,:,jb)        &
    !$ser&     echam_vd_ptvm1=field%tv(:,:,jb)               &
    !$ser&     echam_vd_paclc=field%aclc(:,:,jb)             &
    !$ser&     echam_vd_pxt_emis=pxt_emis                    &
    !$ser&     echam_vd_pthvvar=pthvvar                      &
    !$ser&     echam_vd_pxvar=pxvar                          &
    !$ser&     echam_vd_pz0m_tile=field%z0m_tile(:,jb,:)     &
    !$ser&     echam_vd_ptottem1=field%tottem1(:,:,jb)       &
    !$ser&     echam_vd_pustar=field%ustar(:,jb)             &
    !$ser&     echam_vd_pwstar=pwstar                        &
    !$ser&     echam_vd_pwstar_tile=field%wstar_tile(:,jb,:) &
    !$ser&     echam_vd_pqsat_tile=pqsat_tile                &
    !$ser&     echam_vd_phdtcbl=phdtcbl                      &
    !$ser&     echam_vd_pri=pri                              &
    !$ser&     echam_vd_pri_tile=pri_tile                    &
    !$ser&     echam_vd_pmixlen=pmixlen                      &
    !$ser&     echam_vd_pcfm=pcfm                            &
    !$ser&     echam_vd_pcfm_tile=pcfm_tile                  &
    !$ser&     echam_vd_pcfh=pcfh                            &
    !$ser&     echam_vd_pcfh_tile=pcfh_tile                  &
    !$ser&     echam_vd_pcfv=pcfv                            &
    !$ser&     echam_vd_pcftotte=pcftotte                    &
    !$ser&     echam_vd_pcfthv=pcfthv                        &
    !$ser&     echam_vd_aa=aa                                &
    !$ser&     echam_vd_aa_btm=aa_btm                        &
    !$ser&     echam_vd_bb=bb                                &
    !$ser&     echam_vd_bb_btm=bb_btm                        &
    !$ser&     echam_vd_pfactor_sfc=pfactor_sfc              &
    !$ser&     echam_vd_pcpt_tile=pcpt_tile                  &
    !$ser&     echam_vd_pcptgz=pcptgz                        &
    !$ser&     echam_vd_pzthvvar=pzthvvar                    &
    !$ser&     echam_vd_pthvsig=field%thvsig(:,jb)           &
    !$ser&     echam_vd_pztottevn=pztottevn                  &
    !$ser&     echam_vd_pcsat=field%csat(:,jb)               &
    !$ser&     echam_vd_pcair=field%cair(:,jb)               &
    !$ser&     echam_vd_paz0lh=field%z0h_lnd(:,jb)
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
    !$ser verbatim   CALL init('icon')
    !$ser savepoint echam_vdiff_down-output jg=jg jb=jb jcs=jcs kproma=kproma &
    !$ser&          kbdim=kbdim klev=klev klevm1=klevm1 klevp1=klevp1 &
    !$ser&          ktrac=ktrac ksfc_type=ksfc_type idx_wtr=idx_wtr &
    !$ser&          idx_ice=idx_ice idx_lnd=idx_lnd iqv=iqv iqc=iqc iqi=iqi &
    !$ser&          iqt=iqt pdtime=pdtime date=TRIM(date)
    !$ser mode write
    !$ser data echam_vd_pustar=field%ustar(:,jb)             &
    !$ser&     echam_vd_pwstar=pwstar                        &
    !$ser&     echam_vd_pwstar_tile=field%wstar_tile(:,jb,:) &
    !$ser&     echam_vd_pqsat_tile=pqsat_tile                &
    !$ser&     echam_vd_phdtcbl=phdtcbl                      &
    !$ser&     echam_vd_pri=pri                              &
    !$ser&     echam_vd_pri_tile=pri_tile                    &
    !$ser&     echam_vd_pmixlen=pmixlen                      &
    !$ser&     echam_vd_pcfm=pcfm                            &
    !$ser&     echam_vd_pcfm_tile=pcfm_tile                  &
    !$ser&     echam_vd_pcfh=pcfh                            &
    !$ser&     echam_vd_pcfh_tile=pcfh_tile                  &
    !$ser&     echam_vd_pcfv=pcfv                            &
    !$ser&     echam_vd_pcftotte=pcftotte                    &
    !$ser&     echam_vd_pcfthv=pcfthv                        &
    !$ser&     echam_vd_aa=aa                                &
    !$ser&     echam_vd_aa_btm=aa_btm                        &
    !$ser&     echam_vd_bb=bb                                &
    !$ser&     echam_vd_bb_btm=bb_btm                        &
    !$ser&     echam_vd_pfactor_sfc=pfactor_sfc              &
    !$ser&     echam_vd_pcpt_tile=pcpt_tile                  &
    !$ser&     echam_vd_pcptgz=pcptgz                        &
    !$ser&     echam_vd_pzthvvar=pzthvvar                    &
    !$ser&     echam_vd_pthvsig=field%thvsig(:,jb)           &
    !$ser&     echam_vd_pztottevn=pztottevn                  &
    !$ser&     echam_vd_pch_tile=pch_tile                    &
    !$ser&     echam_vd_pbn_tile=pbn_tile                    &
    !$ser&     echam_vd_pbhn_tile=pbhn_tile                  &
    !$ser&     echam_vd_pbm_tile=pbm_tile                    &
    !$ser&     echam_vd_pbh_tile=pbh_tile
    !$ser verbatim writeOut = .FALSE.
    !$ser verbatim IF (singleStepOut) THEN
    !$ser verbatim   serializeStepOut = .FALSE.
    !$ser verbatim END IF
    !$ser verbatim ENDIF

  END SUBROUTINE serialize_output

END MODULE mo_ser_echam_vdiff_down
