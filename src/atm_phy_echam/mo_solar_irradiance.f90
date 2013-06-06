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

  IMPLICIT NONE

  PRIVATE

  REAL(dp), POINTER :: local_tsi(:) => NULL()
  REAL(dp), POINTER :: local_ssi(:,:) => NULL()

  INTEGER, ALLOCATABLE :: ssi_years(:)
  INTEGER, ALLOCATABLE :: ssi_months(:)

  CHARACTER(len=*), PARAMETER :: ssi_fn = 'bc_ssi.nc'

  PUBLIC :: read_ssi_bc
  PUBLIC :: ssi_time_weights
  PUBLIC :: ssi_time_interpolation

  REAL(dp), SAVE :: wgt1, wgt2
  INTEGER, SAVE  :: nmw1, nmw2

  LOGICAL, SAVE :: ssi_file_read = .FALSE.
  INTEGER, SAVE :: last_year = 0

CONTAINS

  SUBROUTINE read_ssi_bc(year)
    INTEGER, INTENT(in) :: year

    INTEGER :: ncid, ndimid, nvarid
    INTEGER :: ssi_time_entries
    INTEGER :: ssi_numwl
    INTEGER :: start(2), cnt(2)

    INTEGER :: idx, ym1, yp1, first_year, first_month

    ym1 = year - 1
    yp1 = year + 1

    IF (last_year /= year) THEN
      ssi_file_read = .FALSE.
    ENDIF

    IF (ssi_file_read) THEN
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

    ALLOCATE (local_tsi(0:13))
    ALLOCATE (local_ssi(ssi_numwl,0:13))

    CALL nf_check(p_nf_inq_varid(ncid, 'year', nvarid))
    CALL nf_check(p_nf_get_var_int (ncid, nvarid, ssi_years))
    CALL nf_check(p_nf_inq_varid(ncid, 'month', nvarid))
    CALL nf_check(p_nf_get_var_int (ncid, nvarid, ssi_months))

    first_year = ssi_years(1)
    first_month = ssi_months(1)

    ! not adding 1 in calculating the offset leads to an index to December of ym1
    idx = 12*(year - first_year)
    if (idx < 1) then
      CALL finish('','No solar irradiance data available for the requested year')
    endif

    CALL nf_check(p_nf_inq_varid (ncid, 'TSI', nvarid))
    start(1) = idx
    cnt(1) = 14
    CALL nf_check(p_nf_get_vara_double(ncid, nvarid, start, cnt, local_tsi))

    CALL nf_check(p_nf_inq_varid (ncid, 'SSI', nvarid))
    start(1) = 1;   cnt(1) = ssi_numwl;
    start(2) = idx; cnt(2) = 14;
    CALL nf_check(p_nf_get_vara_double (ncid, nvarid, start, cnt, local_ssi))

    ssi_file_read = .TRUE.

    CALL nf_check(p_nf_close(ncid))

    DEALLOCATE(ssi_years)
    DEALLOCATE(ssi_months)

  END SUBROUTINE read_ssi_bc

  SUBROUTINE ssi_time_weights(current_date)

    TYPE(t_datetime), INTENT(in) :: current_date

    ! calculates weighting factores for AMIP sst and sea ice

    TYPE(t_datetime) :: next_date

    TYPE(t_datetime) :: date_monm1, date_monp1
    INTEGER   :: yr, mo, dy, hr, mn, se
    INTEGER   :: isec
    INTEGER   :: imp1, imm1, imlenm1, imlen, imlenp1
    REAL (dp) :: zsec, zdayl
    REAL (dp) :: zmohlf, zmohlfp1, zmohlfm1

    ! time of next timestep and split
    !
    next_date = current_date
    CALL add_time(dtime,0,0,0,next_date)
    CALL date_to_time(next_date)

    yr = next_date%year
    mo = next_date%month
    dy = next_date%day
    hr = next_date%hour
    mn = next_date%minute
    se = INT(next_date%second)

    ! month index for AMIP data  (0..13)
    imp1 = mo+1
    imm1 = mo-1

    ! determine length of months and position within current month

    date_monm1%calendar = next_date%calendar
    IF (imm1 ==  0) THEN
      date_monm1%year = yr-1;  date_monm1%month = 12;   date_monm1%day = 1;
    ELSE
      date_monm1%year = yr;    date_monm1%month = imm1; date_monm1%day = 1;
    ENDIF
    date_monm1%hour = 0;   date_monm1%minute = 0; date_monm1%second   = 0;
    CALL date_to_time(date_monm1)

    date_monp1%calendar = next_date%calendar
    IF (imp1 == 13) THEN
      date_monp1%year = yr+1;  date_monp1%month = 1;    date_monp1%day = 1;
    ELSE
      date_monp1%year = yr;    date_monp1%month = imp1; date_monp1%day = 1;
    ENDIF
    date_monp1%hour     = 0;   date_monp1%minute   = 0;    date_monp1%second   = 0;
    CALL date_to_time(date_monp1)

    imlenm1 = date_monm1%monlen
    imlen   = next_date%monlen
    imlenp1 = date_monp1%monlen

    zdayl    = rdaylen
    zmohlfm1 = imlenm1*zdayl*0.5_dp
    zmohlf   = imlen  *zdayl*0.5_dp
    zmohlfp1 = imlenp1*zdayl*0.5_dp

    ! weighting factors for first/second half of month

    nmw1   = mo

    ! seconds in the present month
    isec = (dy-1) * idaylen + INT(next_date%daysec)
    zsec = REAL(isec,dp)

    IF(zsec <= zmohlf) THEN                     ! first part of month
      wgt1   = (zmohlfm1+zsec)/(zmohlfm1+zmohlf)
      wgt2   = 1.0_dp-wgt1
      nmw2   = imm1
    ELSE                                        ! second part of month
      wgt2   = (zsec-zmohlf)/(zmohlf+zmohlfp1)
      wgt1   = 1.0_dp-wgt2
      nmw2   = imp1
    ENDIF

  END SUBROUTINE ssi_time_weights

  SUBROUTINE ssi_time_interpolation(tsi, ssi)
    REAL(dp), INTENT(out) :: tsi
    REAL(dp), INTENT(out) :: ssi(:)

    tsi    = wgt1 * local_tsi(nmw1) + wgt2 * local_tsi(nmw2)
    ssi(:) = wgt1 * local_ssi(:,nmw1) + wgt2 * local_ssi(:,nmw2)

    CALL message('','Interpolated total solar irradiance and spectral bands of TSI.')

  END SUBROUTINE ssi_time_interpolation

  SUBROUTINE nf_check(iret)
    INTEGER, INTENT(in) :: iret
    IF (iret /= nf_noerr) THEN
      CALL finish('mo_solar_irradiance', nf_strerror(iret))
    ENDIF
  END SUBROUTINE nf_check

END MODULE mo_solar_irradiance
