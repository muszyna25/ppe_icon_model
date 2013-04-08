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
  USE mo_read_netcdf_parallel, ONLY: p_nf_open, p_nf_inq_dimid, p_nf_inq_dimlen, &
       &                             p_nf_inq_varid, p_nf_get_var_double, p_nf_close, &
       &                             nf_read, nf_noerr, nf_strerror, p_nf_get_var_int

  IMPLICIT NONE

  PRIVATE

  REAL(dp), POINTER :: tsi(:) => NULL()
  REAL(dp), POINTER :: ssi(:,:) => NULL()

  INTEGER, ALLOCATABLE :: ssi_years(:)

  CHARACTER(len=*), PARAMETER :: ssi_fn = 'bc_ssi.nc'

  PUBLIC :: read_ssi_bc

  LOGICAL, SAVE :: ssi_file_read = .FALSE.

CONTAINS

  SUBROUTINE read_ssi_bc

    INTEGER :: ncid, ndimid, nvarid
    INTEGER :: ssi_no_years
    INTEGER :: ssi_numwl

    IF (ssi_file_read) THEN
      CALL message('','Solar irradiance data already read ...')
      RETURN
    ENDIF

    CALL nf_check(p_nf_open('bc_ssi.nc', nf_read, ncid))
    CALL nf_check(p_nf_inq_dimid (ncid, 'time', ndimid))
    CALL nf_check(p_nf_inq_dimlen (ncid, ndimid, ssi_no_years))
    CALL nf_check(p_nf_inq_dimid (ncid, 'numwl', ndimid))
    CALL nf_check(p_nf_inq_dimlen (ncid, ndimid, ssi_numwl))

    ALLOCATE (ssi_years(ssi_no_years))
    ALLOCATE (tsi(0:13))
    ALLOCATE (ssi(ssi_numwl,0:13))

    CALL nf_check(p_nf_inq_varid(ncid, 'time', nvarid))
    CALL nf_check(p_nf_get_var_int (ncid, nvarid, ssi_years))

    CALL nf_check(p_nf_inq_varid (ncid, 'TSI', nvarid))
    CALL nf_check(p_nf_get_var_double (ncid, nvarid, tsi))

    CALL nf_check(p_nf_inq_varid (ncid, 'SSI', nvarid))
    CALL nf_check(p_nf_get_var_double (ncid, nvarid, ssi))

    ssi_file_read = .TRUE.

    CALL nf_check(p_nf_close(ncid))

  END SUBROUTINE read_ssi_bc

  SUBROUTINE nf_check(iret)
    INTEGER, INTENT(in) :: iret
    IF (iret /= nf_noerr) THEN
      CALL finish('mo_solar_irradiance', nf_strerror(iret))
    ENDIF
  END SUBROUTINE nf_check

END MODULE mo_solar_irradiance
