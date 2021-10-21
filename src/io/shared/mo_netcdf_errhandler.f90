MODULE mo_netcdf_errhandler

  USE mo_exception,          ONLY: finish, message_text, message, em_warn

  IMPLICIT NONE
  PRIVATE

  INCLUDE 'netcdf.inc'

  PUBLIC :: nf

CONTAINS

  SUBROUTINE nf(errstat, routine, warnonly, silent)

    INTEGER, INTENT(IN)           :: errstat
    CHARACTER(*), INTENT(IN) :: routine
    LOGICAL, INTENT(IN), OPTIONAL :: warnonly, silent
    LOGICAL :: lwarnonly

    IF(PRESENT(silent)) THEN
      IF (silent) RETURN
    END IF
    lwarnonly = .FALSE.
    IF(PRESENT(warnonly)) lwarnonly = .TRUE.
    IF (errstat .NE. nf_noerr) THEN
      IF (lwarnonly) THEN
        CALL message(TRIM(routine)//' netCDF error', nf_strerror(errstat), &
          & level=em_warn)
      ELSE
        CALL finish(TRIM(routine)//' netCDF error', nf_strerror(errstat))
      ENDIF
    ENDIF
  END SUBROUTINE nf

END MODULE mo_netcdf_errhandler
