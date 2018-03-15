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
  USE mo_echam_phy_memory,   ONLY: t_echam_phy_field, t_echam_phy_tend
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

  SUBROUTINE serialize_input(jb, kproma, kbdim, klev, klevm1, ktrac, ksfc_type, &
                             idx_wtr, pdtime, field, aa, pcptgz, pztottevn, bb, &
                             pzthvvar, pxvar)
    INTEGER, INTENT(IN)    :: jb, kproma
    INTEGER, INTENT(INOUT) :: kbdim, klev, klevm1, ktrac
    INTEGER, INTENT(INOUT) :: ksfc_type, idx_wtr
    REAL(wp),INTENT(IN)    :: pdtime
    TYPE(t_echam_phy_field) ,POINTER, INTENT(INOUT) :: field
    REAL(wp),INTENT(INOUT) :: &
      aa(:,:,:,:),         &
      pcptgz(:,:),         &
      pztottevn(:,:),      &
      bb(:,:,:),           &
      pzthvvar(:,:),       & 
      pxvar(:,:)

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: date

    !$ser verbatim IF (singleBlock < 0 .OR. jb == singleBlock) THEN
    !$ser verbatim   writeIn = .TRUE.
    !$ser verbatim ENDIF
    !$ser verbatim IF (serializeStepIn .and. writeIn) THEN
    !$ser verbatim   CALL datetimeToString(time_config%tc_current_date, date)
    !$ser verbatim   CALL init('echam_vdiff_up')
    !$ser savepoint echam_vdiff_up-input jb=jb kproma=kproma iqv=iqv iqc=iqc &
    !$ser&          iqi=iqi iqt=iqt pdtime=pdtime date=TRIM(date)
#if defined SERIALIZE_CREATE_REFERENCE 
    !$ser mode write
#elif defined SERIALIZE_PERTURB_REFERENCE
    !$ser mode read-perturb
#elif defined SERIALIZE_READ_REFERENCE
    !$ser mode read
#else
#error SERIALIZATION MODE IS NOT SET
#endif 
    !$ser data kbdim=kbdim                          &
    !$ser&     klev=klev                            &
    !$ser&     klevm1=klevm1                        &
    !$ser&     ktrac=ktrac                          &
    !$ser&     ksfc_type=ksfc_type                  &
    !$ser&     idx_wtr=idx_wtr                      &
    !$ser&     pfrc=field%frac_tile(:,jb,:)         &
    !$ser&     pcfm_tile=field%cfm_tile(:,jb,:)     &
    !$ser&     aa=aa                                &
    !$ser&     pcptgz=pcptgz                        &
    !$ser&     pum1=field%ua(:,:,jb)                &
    !$ser&     pvm1=field%va(:,:,jb)                &
    !$ser&     ptm1=field%ta(:,:,jb)                &
    !$ser&     pmair=field%mair(:,:,jb)             &
    !$ser&     pmref=field%mref(:,:,jb)             &
    !$ser&     pqm1=field%qtrc(:,:,jb,iqv)          &
    !$ser&     pxlm1=field%qtrc(:,:,jb,iqc)         &
    !$ser&     pxim1,=field%qtrc(:,:,jb,iqi)        &
    !$ser&     pxtm1=field%qtrc(:,:,jb,iqt:)        &
    !$ser&     pgeom1=field%geom(:,:,jb)            &
    !$ser&     pztottevn=pztottevn                  &
    !$ser&     bb=bb                                &
    !$ser&     pzthvvar=pzthvvar                    &
    !$ser&     pxvar=pxvar                          &
    !$ser&     pz0m_tile=field%z0m_tile(:,jb,:) 
    !$ser verbatim writeIn = .FALSE.
    !$ser verbatim IF (singleStepIn) THEN
    !$ser verbatim   serializeStepIn = .FALSE.
    !$ser verbatim END IF
    !$ser verbatim ENDIF

  END SUBROUTINE serialize_input

  SUBROUTINE serialize_output(jb, kproma, kbdim, klev, klevm1, ktrac, &
                              ksfc_type, idx_wtr, pdtime, field, bb, pxvar, &
                              tend, pthvvar)
    INTEGER, INTENT(IN)    :: jb, kproma
    INTEGER, INTENT(INOUT) :: kbdim, klev, klevm1, ktrac
    INTEGER, INTENT(INOUT) :: ksfc_type, idx_wtr
    REAL(wp),INTENT(IN)    :: pdtime
    TYPE(t_echam_phy_field) ,POINTER, INTENT(INOUT) :: field
    TYPE(t_echam_phy_tend)  ,POINTER, INTENT(INOUT) :: tend
    REAL(wp),INTENT(INOUT) :: &
      bb(:,:,:),           &
      pxvar(:,:),          &
      pthvvar(:,:)

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: date

    !$ser verbatim IF (singleBlock < 0 .OR. jb == singleBlock) THEN
    !$ser verbatim   writeOut = .TRUE.
    !$ser verbatim ENDIF
    !$ser verbatim IF (serializeStepOut .and. writeOut) THEN
    !$ser verbatim CALL datetimeToString(time_config%tc_current_date, date)
    !$ser verbatim CALL init('echam_vdiff_up')
    !$ser savepoint echam_vdiff_up-output jb=jb kproma=kproma iqv=iqv iqc=iqc &
    !$ser&          iqi=iqi iqt=iqt pdtime=pdtime date=TRIM(date)
    !$ser mode write
    !$ser data kbdim=kbdim                            &
    !$ser&     klev=klev                              &
    !$ser&     klevm1=klevm1                          &
    !$ser&     ktrac=ktrac                            &
    !$ser&     ksfc_type=ksfc_type                    &
    !$ser&     idx_wtr=idx_wtr                        &
    !$ser&     bb=bb                                  &
    !$ser&     pxvar=pxvar                            &
    !$ser&     pz0m_tile=field%z0m_tile(:,jb,:)       &
    !$ser&     pkedisp=field%kedisp(:,jb)             &
    !$ser&     pute_vdf=tend%ua_vdf(:,:,jb)           &
    !$ser&     pvte_vdf=tend%va_vdf(:,:,jb)           &
    !$ser&     pq_vdf=field%q_vdf(:,:,jb)             &
    !$ser&     pqte_vdf=tend%qtrc_vdf(:,:,jb,iqv)     &
    !$ser&     pxlte_vdf=tend%qtrc_vdf(:,:,jb,iqc)    &
    !$ser&     pxite_vdf=tend%qtrc_vdf(:,:,jb,iqi)    &
    !$ser&     pxtte_vdf=tend%qtrc_vdf(:,:,jb,iqt:)   &
    !$ser&     pz0m=field%z0m(:,jb)                   &
    !$ser&     pthvvar=pthvvar                        &
    !$ser&     ptotte=field%totte(:,:,jb)             &
    !$ser&     psh_vdiff=field%sh_vdiff(:,jb)         &
    !$ser&     pqv_vdiff=field%qv_vdiff(:,jb)
    !$ser verbatim writeOut = .FALSE.
    !$ser verbatim IF (singleStepOut) THEN
    !$ser verbatim   serializeStepOut = .FALSE.
    !$ser verbatim END IF
    !$ser verbatim ENDIF

  END SUBROUTINE serialize_output

END MODULE mo_ser_echam_vdiff_up
