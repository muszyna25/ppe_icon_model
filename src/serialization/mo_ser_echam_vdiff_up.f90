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
  USE mo_echam_phy_memory,    ONLY: t_echam_phy_field
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
                             idx_wtr, pdtime, field)
    INTEGER, INTENT(IN) :: jb
    INTEGER, INTENT(IN) :: kproma, kbdim, klev, klevm1, ktrac
    INTEGER, INTENT(IN) :: ksfc_type, idx_wtr
    REAL(wp),INTENT(IN) :: pdtime
    TYPE(t_echam_phy_field) ,POINTER :: field

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: date

    !$ser verbatim IF (singleBlock < 0 .OR. jb == singleBlock) THEN
    !$ser verbatim   writeIn = .TRUE.
    !$ser verbatim ENDIF
    !$ser verbatim IF (serializeStepIn .and. writeIn) THEN
    !$ser verbatim   CALL datetimeToString(time_config%tc_current_date, date)
    !$ser verbatim   CALL init('echam_vdiff_up')
    !$ser savepoint echam_vdiff_up-input jb=jb date=TRIM(date)
#if defined SERIALIZE_CREATE_REFERENCE 
    !$ser mode write
#elif defined SERIALIZE_PERTURB_REFERENCE
    !$ser mode read-perturb
#elif defined SERIALIZE_READ_REFERENCE
    !$ser mode read
#else
#error SERIALIZATION MODE IS NOT SET
#endif 
    !$ser data kproma=kproma                        &
    !$ser&     kbdim=kbdim                          &
    !$ser&     klev=klev                            &
    !$ser&     klevm1=klevm1                        &
    !$ser&     ktrac=ktrac                          &
    !$ser&     ksfc_type=ksfc_type                  &
    !$ser&     idx_wtr=idx_wtr                      &
    !$ser&     pdtime=pdtime
    !$ser verbatim writeIn = .FALSE.
    !$ser verbatim IF (singleStepIn) THEN
    !$ser verbatim   serializeStepIn = .FALSE.
    !$ser verbatim END IF
    !$ser verbatim ENDIF

  END SUBROUTINE serialize_input

  SUBROUTINE serialize_output(jb, kproma, kbdim, klev, klevm1, ktrac, &
                              ksfc_type, idx_wtr, pdtime, field)
    INTEGER, INTENT(IN) :: jb
    INTEGER, INTENT(IN) :: kproma, kbdim, klev, klevm1, ktrac
    INTEGER, INTENT(IN) :: ksfc_type, idx_wtr
    REAL(wp),INTENT(IN) :: pdtime
    TYPE(t_echam_phy_field) ,POINTER :: field

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: date

    !$ser verbatim IF (singleBlock < 0 .OR. jb == singleBlock) THEN
    !$ser verbatim   writeOut = .TRUE.
    !$ser verbatim ENDIF
    !$ser verbatim IF (serializeStepOut .and. writeOut) THEN
    !$ser verbatim CALL datetimeToString(time_config%tc_current_date, date)
    !$ser verbatim CALL init('echam_vdiff_up')
    !$ser savepoint echam_vdiff_up-output jb=jb date=TRIM(date)
    !$ser mode write
    !$ser data kproma=kproma                        &
    !$ser&     kbdim=kbdim                          &
    !$ser&     klev=klev                            &
    !$ser&     klevm1=klevm1                        &
    !$ser&     ktrac=ktrac                          &
    !$ser&     ksfc_type=ksfc_type                  &
    !$ser&     idx_wtr=idx_wtr                      &
    !$ser&     pdtime=pdtime
    !$ser verbatim writeOut = .FALSE.
    !$ser verbatim IF (singleStepOut) THEN
    !$ser verbatim   serializeStepOut = .FALSE.
    !$ser verbatim END IF
    !$ser verbatim ENDIF

  END SUBROUTINE serialize_output

END MODULE mo_ser_echam_vdiff_up
