!--------------------------------------------------------------------
!
! Serialization routine for ECHAM vdiff up
!
!--------------------------------------------------------------------

MODULE mo_ser_echam_vdiff_up

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

  SUBROUTINE serialize_input(jb, jcs, kproma, kbdim, klev, klevm1, ktrac, ksfc_type, &
                             idx_wtr, pdtime, field, pcfm_tile, aa, pcptgz, pztottevn, bb, &
                             pzthvvar, pxvar, pkedisp)
    INTEGER, INTENT(IN)    :: jb, jcs, kproma, kbdim, klev, klevm1, ktrac, ksfc_type, idx_wtr
    REAL(wp),INTENT(IN)    :: pdtime
    TYPE(t_echam_phy_field) ,POINTER, INTENT(INOUT) :: field
    REAL(wp),INTENT(INOUT) :: &
      pcfm_tile(:,:),       &
      aa(:,:,:,:),          &
      pcptgz(:,:),          &
      pztottevn(:,:),       &
      bb(:,:,:),            &
      pzthvvar(:,:),        &
      pxvar(:,:),           &
      pkedisp(:)

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: date

    !$ser verbatim IF (singleBlock < 0 .OR. jb == singleBlock) THEN
    !$ser verbatim   writeIn = .TRUE.
    !$ser verbatim ENDIF
    !$ser verbatim IF (serializeStepIn .and. writeIn) THEN
    !$ser verbatim   CALL datetimeToString(time_config%tc_current_date, date)
    !$ser verbatim   CALL init('icon')
    !$ser savepoint echam_vdiff_up-input jb=jb jcs=jcs kproma=kproma kbdim=kbdim klev=klev &
    !$ser&          klevm1=klevm1 ktrac=ktrac ksfc_type=ksfc_type idx_wtr=idx_wtr &
    !$ser&          pdtime=pdtime iqv=iqv iqc=iqc iqi=iqi iqt=iqt date=TRIM(date)
#if defined SERIALIZE_CREATE_REFERENCE
    !$ser mode write
#elif defined SERIALIZE_PERTURB_REFERENCE
    !$ser mode read-perturb
#elif defined SERIALIZE_READ_REFERENCE
    !$ser mode read
#else
#error SERIALIZATION MODE IS NOT SET
#endif
    !$ser data echam_vu_pfrc=field%frac_tile(:,jb,:)         &
    !$ser&     echam_vu_pcfm_tile=pcfm_tile                  &
    !$ser&     echam_vu_aa=aa                                &
    !$ser&     echam_vu_pcptgz=pcptgz                        &
    !$ser&     echam_vu_pum1=field%ua(:,:,jb)                &
    !$ser&     echam_vu_pvm1=field%va(:,:,jb)                &
    !$ser&     echam_vu_ptm1=field%ta(:,:,jb)                &
    !$ser&     echam_vu_pmair=field%mair(:,:,jb)             &
    !$ser&     echam_vu_pmref=field%mref(:,:,jb)             &
    !$ser&     echam_vu_pqm1=field%qtrc(:,:,jb,iqv)          &
    !$ser&     echam_vu_pxlm1=field%qtrc(:,:,jb,iqc)         &
    !$ser&     echam_vu_pxim1=field%qtrc(:,:,jb,iqi)         &
    !$ser&     echam_vu_pxtm1=field%qtrc(:,:,jb,iqt:)        &
    !$ser&     echam_vu_pgeom1=field%geom(:,:,jb)            &
    !$ser&     echam_vu_pztottevn=pztottevn                  &
    !$ser&     echam_vu_bb=bb                                &
    !$ser&     echam_vu_pzthvvar=pzthvvar                    &
    !$ser&     echam_vu_pxvar=pxvar                          &
    !$ser&     echam_vu_pz0m_tile=field%z0m_tile(:,jb,:)
    !$ser verbatim writeIn = .FALSE.
    !$ser verbatim IF (singleStepIn) THEN
    !$ser verbatim   serializeStepIn = .FALSE.
    !$ser verbatim END IF
    !$ser verbatim ENDIF

  END SUBROUTINE serialize_input

  SUBROUTINE serialize_output(jb, jcs, kproma, kbdim, klev, klevm1, ktrac, &
                              ksfc_type, idx_wtr, pdtime, field, bb, pxvar, &
                              pkedisp, pute_vdf, pvte_vdf, pq_vdf, tend_qtrc_vdf, &
                              pthvvar)
    INTEGER, INTENT(IN)    :: jb, jcs, kproma, kbdim, klev, klevm1, ktrac, ksfc_type, idx_wtr
    REAL(wp),INTENT(IN)    :: pdtime
    TYPE(t_echam_phy_field) ,POINTER, INTENT(INOUT) :: field
    REAL(wp),INTENT(INOUT) :: &
      bb(:,:,:),            &
      pxvar(:,:),           &
      pkedisp(:),           &
      pute_vdf(:,:),        &
      pvte_vdf(:,:),        &
      pq_vdf(:,:),          &
      tend_qtrc_vdf(:,:,:), &
      pthvvar(:,:)

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: date

    !$ser verbatim IF (singleBlock < 0 .OR. jb == singleBlock) THEN
    !$ser verbatim   writeOut = .TRUE.
    !$ser verbatim ENDIF
    !$ser verbatim IF (serializeStepOut .and. writeOut) THEN
    !$ser verbatim CALL datetimeToString(time_config%tc_current_date, date)
    !$ser verbatim CALL init('icon')
    !$ser savepoint echam_vdiff_up-output jb=jb jcs=jcs kproma=kproma kbdim=kbdim klev=klev &
    !$ser&          klevm1=klevm1 ktrac=ktrac ksfc_type=ksfc_type idx_wtr=idx_wtr &
    !$ser&          pdtime=pdtime iqv=iqv iqc=iqc iqi=iqi iqt=iqt date=TRIM(date)
    !$ser mode write
    !$ser data echam_vu_bb=bb                                  &
    !$ser&     echam_vu_pxvar=pxvar                            &
    !$ser&     echam_vu_pz0m_tile=field%z0m_tile(:,jb,:)       &
    !$ser&     echam_vu_pkedisp=pkedisp                        &
    !$ser&     echam_vu_pute_vdf=pute_vdf                      &
    !$ser&     echam_vu_pvte_vdf=pvte_vdf                      &
    !$ser&     echam_vu_pq_vdf=pq_vdf                          &
    !$ser&     echam_vu_pqte_vdf=tend_qtrc_vdf(:,:,iqv)        &
    !$ser&     echam_vu_pxlte_vdf=tend_qtrc_vdf(:,:,iqc)       &
    !$ser&     echam_vu_pxite_vdf=tend_qtrc_vdf(:,:,iqi)       &
    !$ser&     echam_vu_pxtte_vdf=tend_qtrc_vdf(:,:,iqt:)      &
    !$ser&     echam_vu_pz0m=field%z0m(:,jb)                   &
    !$ser&     echam_vu_pthvvar=pthvvar                        &
    !$ser&     echam_vu_ptotte=field%totte(:,:,jb)             &
    !$ser&     echam_vu_psh_vdiff=field%sh_vdiff(:,jb)         &
    !$ser&     echam_vu_pqv_vdiff=field%qv_vdiff(:,jb)
    !$ser verbatim writeOut = .FALSE.
    !$ser verbatim IF (singleStepOut) THEN
    !$ser verbatim   serializeStepOut = .FALSE.
    !$ser verbatim END IF
    !$ser verbatim ENDIF

  END SUBROUTINE serialize_output

END MODULE mo_ser_echam_vdiff_up
