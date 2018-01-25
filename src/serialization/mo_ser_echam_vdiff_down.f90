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
  USE mo_echam_phy_memory,    ONLY: t_echam_phy_field
  IMPLICIT NONE

  LOGICAL :: writeIn = .TRUE.
  LOGICAL :: writeOut = .TRUE.

  PUBLIC :: serialize_input
  PUBLIC :: serialize_output

  CONTAINS

  SUBROUTINE serialize_input(jblock, kproma, kbdim, klev, klevm1, klevp1, ktrac, &
                             ksfc_type, idx_wtr, idx_ice, idx_lnd, pdtime, &
                             field)
    INTEGER, INTENT(IN) :: jblock
    INTEGER, INTENT(IN) :: kproma, kbdim, klev, klevm1, klevp1, ktrac
    INTEGER, INTENT(IN) :: ksfc_type, idx_wtr, idx_ice, idx_lnd
    REAL(wp),INTENT(IN) :: pdtime
    TYPE(t_echam_phy_field) ,POINTER :: field

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: date

    !$ser verbatim IF (writeIn) THEN
    !$ser verbatim call datetimeToString(time_config%tc_current_date, date)
    !$ser verbatim call init('echam_vdiff_down')
    !$ser savepoint echam_vdiff_down-input jblock=jblock date=TRIM(date)
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
    !$ser&     klevp1=klevp1                        &
    !$ser&     ktrac=ktrac                          &
    !$ser&     ksfc_type=ksfc_type                  &
    !$ser&     idx_wtr=idx_wtr                      &
    !$ser&     idx_ice=idx_ice                      &
    !$ser&     idx_lnd=idx_lnd                      &
    !$ser&     pdtime=pdtime
    !NOser verbatim writeIn = .FALSE.
    !$ser verbatim ENDIF

  END SUBROUTINE serialize_input

  SUBROUTINE serialize_output(jblock, kproma, kbdim, klev, klevm1, klevp1, ktrac, &
                             ksfc_type, idx_wtr, idx_ice, idx_lnd, pdtime, &
                             field)
    INTEGER, INTENT(IN) :: jblock
    INTEGER, INTENT(IN) :: kproma, kbdim, klev, klevm1, klevp1, ktrac
    INTEGER, INTENT(IN) :: ksfc_type, idx_wtr, idx_ice, idx_lnd
    REAL(wp),INTENT(IN) :: pdtime
    TYPE(t_echam_phy_field) ,POINTER :: field

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: date

    !$ser verbatim if (writeOut) then
    !$ser verbatim call datetimeToString(time_config%tc_current_date, date)
    !$ser verbatim call init('echam_vdiff_down')
    !$ser savepoint echam_vdiff_down-output jblock=jblock date=TRIM(date)
    !$ser mode write
    !$ser data kproma=kproma                        &
    !$ser&     kbdim=kbdim                          &
    !$ser&     klev=klev                            &
    !$ser&     klevm1=klevm1                        &
    !$ser&     klevp1=klevp1                        &
    !$ser&     ktrac=ktrac                          &
    !$ser&     ksfc_type=ksfc_type                  &
    !$ser&     idx_wtr=idx_wtr                      &
    !$ser&     idx_ice=idx_ice                      &
    !$ser&     idx_lnd=idx_lnd                      &
    !$ser&     pdtime=pdtime
    !NOser verbatim writeOut = .FALSE.
    !$ser verbatim endif

  END SUBROUTINE serialize_output

END MODULE mo_ser_echam_vdiff_down
