MODULE mo_psrad_flat_data

  USE mo_psrad_general, ONLY: wp, max_minor_species, nbndlw, &
    nbndsw, npressure, ngas
  IMPLICIT NONE

  INCLUDE 'netcdf.inc'

PUBLIC

  TYPE dope
    INTEGER :: o, s, n
  END TYPE
  TYPE(dope), DIMENSION(2,nbndlw) :: lw_kmajor, lw_h2oref, lw_planck, &
    lw_kgas(max_minor_species,2,nbndlw)
  TYPE(dope), DIMENSION(2,nbndsw) :: sw_kmajor, sw_h2oref, sw_kgas, sw_rayl, &
    sw_sfluxref(nbndsw)
  TYPE(dope) :: lw_info, sw_info, lw_band_info(nbndlw), sw_band_info(nbndsw)
  INTEGER :: lw_count, sw_count

  INTEGER :: flat_data_size
  REAL(wp), ALLOCATABLE :: flat_data(:)

  PUBLIC :: setup_flat_data
CONTAINS

  SUBROUTINE setup_flat_data(pressure_scale, droplet_scale)

    REAL(wp), INTENT(IN) :: pressure_scale, droplet_scale
    CALL read_flat_data(pressure_scale, droplet_scale)
!$ACC ENTER DATA COPYIN(flat_data)

  END SUBROUTINE setup_flat_data

  SUBROUTINE count_dope(n_dope)
    INTEGER, INTENT(OUT) :: n_dope
    n_dope = 2 + &
      nbndlw * (2 * (3 + max_minor_species) + 1 ) + &
      nbndsw * (2 * 4 + 2)
  END SUBROUTINE count_dope

  SUBROUTINE read_dope(d, offset, o, s, n)
    TYPE(dope), INTENT(OUT) :: d
    INTEGER, INTENT(INOUT) :: offset
    INTEGER, INTENT(IN) :: o(:), s(:), n(:)
    d%o = o(offset)
    d%s = s(offset)
    d%n = n(offset)
    offset = offset + 1
  END SUBROUTINE read_dope

  SUBROUTINE unpack_dope(o, s, n)
    INTEGER, INTENT(IN) :: o(:), s(:), n(:)
    INTEGER :: i, j, k, offset

    offset = 1
    CALL read_dope(lw_info, offset, o, s, n)
    CALL read_dope(sw_info, offset, o, s, n)
    DO j = 1, nbndlw
      DO i = 1, 2
        CALL read_dope(lw_kmajor(i,j), offset, o, s, n)
        CALL read_dope(lw_h2oref(i,j), offset, o, s, n)
        CALL read_dope(lw_planck(i,j), offset, o, s, n)
        DO k = 1, max_minor_species
          CALL read_dope(lw_kgas(k,i,j), offset, o, s, n)
        ENDDO
      ENDDO
      CALL read_dope(lw_band_info(j), offset, o, s, n)
    ENDDO
    DO j = 1, nbndsw
      DO i = 1, 2
        CALL read_dope(sw_kmajor(i,j), offset, o, s, n)
        CALL read_dope(sw_h2oref(i,j), offset, o, s, n)
        CALL read_dope(sw_kgas(i,j), offset, o, s, n)
        CALL read_dope(sw_rayl(i,j), offset, o, s, n)
      ENDDO
      CALL read_dope(sw_sfluxref(j), offset, o, s, n)
      CALL read_dope(sw_band_info(j), offset, o, s, n)
    ENDDO
  END SUBROUTINE unpack_dope

  SUBROUTINE read_flat_data(pressure_scale, droplet_scale)
    USE mo_psrad_lrtm_kgs, ONLY: lrtm_count_flat, lrtm_unpack_data
    USE mo_psrad_general, ONLY: finish
    REAL(wp), INTENT(IN) :: pressure_scale, droplet_scale

    INTEGER :: n_dope 
    INTEGER, ALLOCATABLE :: flat_dope_o(:), flat_dope_s(:), flat_dope_n(:)
    INTEGER :: n_lrtm
    REAL(wp), ALLOCATABLE :: flat_lrtm(:)
    INTEGER, PARAMETER :: n_integer = 2, n_real = 2
    INTEGER :: flat_integer(n_integer), ret, nf_id, idummy, fid
    REAL(wp) :: flat_real(n_real), dims(5)
    INTEGER :: idx_real, idx_integer, idx_lrtm, idx_flat_data, &
      idx_dope_o, idx_dope_s, idx_dope_n, nf_length, i

    CHARACTER(len=16) :: dim_names(5) = (/ &
      'n_real          ', &
      'n_integer       ', &
      'n_lrtm          ', &
      'n_dope          ', &
      'flat_data_size  '/), &
      nf_name

    CALL count_dope(n_dope)
    ALLOCATE(flat_dope_o(n_dope))
    ALLOCATE(flat_dope_s(n_dope))
    ALLOCATE(flat_dope_n(n_dope))

    CALL lrtm_count_flat(n_lrtm)
    ALLOCATE(flat_lrtm(n_lrtm))

    dims(1:4) = (/n_real, n_integer, n_lrtm, n_dope/)

    ret = nf_open('lsdata.nc', nf_nowrite, fid)
    DO i = 1, 5
      ret = nf_inq_dim(fid, i, nf_name, nf_length)
      IF (TRIM(nf_name) /= TRIM(dim_names(i))) THEN
        WRITE(*,'(a,i1,a)') 'While reading lsdata.nc: dimension ', i, &
          ' is "' // TRIM(nf_name) // '", expected "' // &
          TRIM(dim_names(i)) // '"'
        STOP
      ENDIF
      IF (i == 5) EXIT
      IF (nf_length /= dims(i)) THEN
        WRITE(*,'(a,i1,a,i5,a,i5)') 'While reading lsdata.nc: dimension ', i, &
          ' has size ', nf_length, ', expected ', dims(i)
        STOP
      ENDIF
    ENDDO
    flat_data_size = nf_length
    ALLOCATE(flat_data(flat_data_size))

    ret = nf_inq_varid(fid, 'flat_real', idx_real)
    ret = nf_inq_varid(fid, 'flat_integer', idx_integer)
    ret = nf_inq_varid(fid, 'flat_lrtm', idx_lrtm)
    ret = nf_inq_varid(fid, 'flat_dope_o', idx_dope_o)
    ret = nf_inq_varid(fid, 'flat_dope_s', idx_dope_s)
    ret = nf_inq_varid(fid, 'flat_dope_n', idx_dope_n)
    ret = nf_inq_varid(fid, 'flat_data', idx_flat_data)

    ret = nf_get_vara_double(fid, idx_real, [1], [n_real], flat_real)
    ret = nf_get_vara_int(fid, idx_integer, [1], [n_integer], flat_integer)
    ret = nf_get_vara_double(fid, idx_lrtm, [1], [n_lrtm], flat_lrtm)
    ret = nf_get_vara_int(fid, idx_dope_o, [1], [n_dope], flat_dope_o)
    ret = nf_get_vara_int(fid, idx_dope_s, [1], [n_dope], flat_dope_s)
    ret = nf_get_vara_int(fid, idx_dope_n, [1], [n_dope], flat_dope_n)
    ret = nf_get_vara_double(fid, idx_flat_data, [1], [flat_data_size], &
      flat_data)

    ret = nf_close(fid)

    IF ((pressure_scale - flat_real(1)) / flat_real(1) > 1e-14) THEN
      CALL finish('psrad_flat_data','Initialization value for pressure ' // &
        'scale incompatible with the one used for generating data.')
    ENDIF
    IF ((droplet_scale - flat_real(2)) / flat_real(2) > 1e-14) THEN
      CALL finish('psrad_flat_data','Initialization value for droplet ' // &
        'scale incompatible with the one used for generating data.')
    ENDIF
    CALL unpack_dope(flat_dope_o, flat_dope_s, flat_dope_n)
    CALL lrtm_unpack_data(flat_lrtm)

  END SUBROUTINE read_flat_data

END MODULE mo_psrad_flat_data

