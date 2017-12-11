MODULE mo_psrad_io

  USE mo_psrad_general, ONLY : dp, finish

  USE mo_netcdf_parallel, ONLY: p_nf_open, p_nf_close, &
    &                           p_nf_inq_varid,        &
    &                           p_nf_get_vara_double,  &
    &                           nf_read, nf_noerr

  PUBLIC :: psrad_io_open, psrad_io_close, psrad_io_copy_double

CONTAINS

  SUBROUTINE psrad_io_open(filename, fileid)
    CHARACTER(len=*), INTENT(IN) :: filename
    INTEGER, INTENT(OUT) :: fileid
    INTEGER :: nf_status

    nf_status = p_nf_open(filename, nf_read, fileid)
    IF (nf_status /= nf_noerr) THEN
      CALL finish('psrad_io_open', 'p_nf_close /= nf_noerr')
    ENDIF
  END SUBROUTINE psrad_io_open

  SUBROUTINE psrad_io_close(fileid)
    INTEGER, INTENT(IN) :: fileid
    INTEGER :: nf_status ! BUG: generates warning, but won't compile...

    nf_status = p_nf_close(fileid) !...unless this return value is stored
    IF (nf_status /= nf_noerr) THEN
      CALL finish('psrad_io_close', 'p_nf_close /= nf_noerr')
    ENDIF
  END SUBROUTINE psrad_io_close

  SUBROUTINE psrad_io_copy_double(fileid, varname, src, tgt, buffer)
    INTEGER, INTENT(IN) :: fileid
    CHARACTER(len=*), INTENT(IN) :: varname
    INTEGER, INTENT(IN) :: src(*), tgt(*)
    REAL(dp), INTENT(OUT) :: buffer(*)
    INTEGER :: nf_status, varid

    nf_status = p_nf_inq_varid(fileid, varname, varid)
    IF (nf_status /= nf_noerr) THEN
      CALL finish('psrad_io_copy_double', 'p_nf_inq_varid /= nf_noerr')
    ENDIF
    nf_status = p_nf_get_vara_double(fileid, varid, src, tgt, buffer)
    IF (nf_status /= nf_noerr) THEN
      CALL finish('psrad_io_copy_double', 'p_nf_get_vara_double /= nf_noerr')
    ENDIF

  END SUBROUTINE psrad_io_copy_double
END MODULE mo_psrad_io
