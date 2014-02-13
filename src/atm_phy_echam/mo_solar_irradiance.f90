!>
!! Preliminary read and time interpolation of solar irradiance data
!!
!! This is  a clone of the respective ECHAM routine
!!
!! Read spectrally resolved solar irradiance yearly, apply primitive time
!! interpolation to monthly mean values and apply
!!
!! J.S. Rast, MPI, March 2010, original version for echam6
!! L. Kornblueh, MPI, April 2013, adapted as temporary reader in ICON
!!
MODULE mo_solar_irradiance

  USE mo_kind,                 ONLY: dp
  USE mo_exception,            ONLY: finish, message, message_text
  USE mo_datetime,             ONLY: t_datetime, add_time, date_to_time, idaylen, rdaylen
  USE mo_run_config,           ONLY: dtime
  USE mo_read_netcdf_parallel, ONLY: p_nf_open, p_nf_inq_dimid, p_nf_inq_dimlen, &
       &                             p_nf_inq_varid, p_nf_get_vara_double, p_nf_close, &
       &                             nf_read, nf_noerr, nf_strerror, p_nf_get_var_int
  USE mo_time_interpolation_weights,ONLY: t_wi_limm

  IMPLICIT NONE

  PRIVATE

  REAL(dp), POINTER :: tsi_radt_m(:) => NULL(), tsi_m(:) => NULL()
  REAL(dp), POINTER :: ssi_radt_m(:,:) => NULL()

  INTEGER, ALLOCATABLE :: ssi_years(:)
  INTEGER, ALLOCATABLE :: ssi_months(:)

  CHARACTER(len=*), PARAMETER :: ssi_fn = 'bc_ssi.nc'

  PUBLIC :: read_ssi_bc, ssi_time_interpolation

  REAL(dp), SAVE :: wgt1, wgt2
  INTEGER, SAVE  :: nmw1, nmw2

  LOGICAL, SAVE :: lread_solar = .TRUE., lread_solar_radt = .TRUE.
  INTEGER, SAVE :: last_year = -999999, last_year_radt = -999999

CONTAINS

  SUBROUTINE read_ssi_bc(year, lradt)
    INTEGER, INTENT(in) :: year
    LOGICAL, INTENT(in) :: lradt ! lradt=.true.: read data for radiation time step
                                 ! the radiation time step may be in the future
                                 ! lradt=.false.: read data for heating rates only
                                 ! data needed at every integration time step

    INTEGER :: ncid, ndimid, nvarid
    INTEGER :: ssi_time_entries
    INTEGER :: ssi_numwl
    INTEGER :: start(2), cnt(2)

    INTEGER :: idx, ym1, yp1, first_year, first_month

    ym1 = year - 1
    yp1 = year + 1

    IF (lradt) THEN
       IF (last_year_radt /= year ) lread_solar_radt=.TRUE.
    ELSE
       IF (last_year /= year) lread_solar=.TRUE.
    ENDIF

    IF (.NOT.lread_solar_radt .AND. .NOT.lread_solar) THEN
      CALL message('','Solar irradiance data already read ...')
      RETURN
    ENDIF

    CALL nf_check(p_nf_open('bc_ssi.nc', nf_read, ncid))

    CALL nf_check(p_nf_inq_dimid(ncid, 'time', ndimid))
    CALL nf_check(p_nf_inq_dimlen(ncid, ndimid, ssi_time_entries))
    CALL nf_check(p_nf_inq_dimid(ncid, 'numwl', ndimid))
    CALL nf_check(p_nf_inq_dimlen(ncid, ndimid, ssi_numwl))

    ALLOCATE (ssi_years(ssi_time_entries))
    ALLOCATE (ssi_months(ssi_time_entries))

    IF (lradt) THEN
       IF (.NOT.(ASSOCIATED(tsi_radt_m))) ALLOCATE(tsi_radt_m(0:13))
       IF (.NOT.(ASSOCIATED(ssi_radt_m))) ALLOCATE(ssi_radt_m(ssi_numwl,0:13))
    ELSE
       IF (.NOT.(ASSOCIATED(tsi_m)))      ALLOCATE(tsi_m(0:13))
    END IF

    CALL nf_check(p_nf_inq_varid(ncid, 'year', nvarid))
    CALL nf_check(p_nf_get_var_int (ncid, nvarid, ssi_years))
    CALL nf_check(p_nf_inq_varid(ncid, 'month', nvarid))
    CALL nf_check(p_nf_get_var_int (ncid, nvarid, ssi_months))

    first_year = ssi_years(1)
    first_month = ssi_months(1)

    ! not adding 1 in calculating the offset leads to an index to December of ym1
    idx = 12*(year - first_year)
    IF (idx < 1) THEN
      CALL finish('','No solar irradiance data available for the requested year')
    END IF

    CALL nf_check(p_nf_inq_varid (ncid, 'TSI', nvarid))
    start(1) = idx
    cnt(1) = 14
    IF (lradt) THEN
       CALL nf_check(p_nf_get_vara_double(ncid, nvarid, start, cnt, tsi_radt_m))
       CALL nf_check(p_nf_inq_varid (ncid, 'SSI', nvarid))
       start(1) = 1;   cnt(1) = ssi_numwl;
       start(2) = idx; cnt(2) = 14;
       CALL nf_check(p_nf_get_vara_double (ncid, nvarid, start, cnt, ssi_radt_m))
       lread_solar_radt=.FALSE.
    ELSE
       CALL nf_check(p_nf_get_vara_double(ncid, nvarid, start, cnt, tsi_m))
       lread_solar=.FALSE.
    END IF

    CALL nf_check(p_nf_close(ncid))

    DEALLOCATE(ssi_years)
    DEALLOCATE(ssi_months)

  END SUBROUTINE read_ssi_bc

  SUBROUTINE ssi_time_interpolation(wi, lradt, tsi, ssi)
    TYPE(t_wi_limm), INTENT(in)     :: wi
    LOGICAL, INTENT(in)             :: lradt
    REAL(dp), INTENT(out)           :: tsi
    REAL(dp), INTENT(out), OPTIONAL :: ssi(:)
    CHARACTER(len=14)               :: ctsi

    IF (lradt) THEN
      IF (.NOT.PRESENT(ssi)) THEN
        CALL finish ('ssi_time_interplation of mo_solar_irradiance', &
                     'Interpolation to radiation time step needs ssi',exit_no=1)
      END IF
      tsi    = wi%wgt1 * tsi_radt_m(wi%inm1) + wi%wgt2 * tsi_radt_m(wi%inm2)
      ssi(:) = wi%wgt1*ssi_radt_m(:,wi%inm1) + wi%wgt2*ssi_radt_m(:,wi%inm2)
      WRITE(ctsi,'(F14.8)') tsi
      CALL message('','Interpolated total solar irradiance and spectral ' &
                      //'bands for radiation transfer, tsi= '//ctsi)

    ELSE
      IF (PRESENT(ssi)) THEN
        CALL message ('ssi_time_interplation of mo_solar_irradiance', &
                     'Interpolation of ssi not necessary')
      END IF
      tsi    = wi%wgt1 * tsi_m(wi%inm1) + wi%wgt2 * tsi_m(wi%inm2)
    END IF
       
  END SUBROUTINE ssi_time_interpolation


  SUBROUTINE nf_check(iret)
    INTEGER, INTENT(in) :: iret
    IF (iret /= nf_noerr) THEN
      CALL finish('mo_solar_irradiance', nf_strerror(iret))
    ENDIF
  END SUBROUTINE nf_check

END MODULE mo_solar_irradiance
