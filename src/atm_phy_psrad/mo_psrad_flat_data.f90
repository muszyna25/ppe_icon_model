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

#ifdef PSRAD_WITH_LEGACY
  INTERFACE flatten
    MODULE PROCEDURE flatten1
    MODULE PROCEDURE flatten2
    MODULE PROCEDURE flatten3
  END INTERFACE flatten
#endif

  PUBLIC :: setup_flat_data
CONTAINS

  SUBROUTINE setup_flat_data(pressure_scale, droplet_scale)

    REAL(wp), INTENT(IN) :: pressure_scale, droplet_scale
#ifdef PSRAD_WITH_LEGACY
    INTEGER :: offset

    offset = 0
    CALL zero_all_dope
    CALL measure_lw_data(offset)
    CALL measure_sw_data(offset)

    flat_data_size = lw_info%s + sw_info%s
    ALLOCATE(flat_data(flat_data_size))

    CALL flatten_lw_data
    CALL flatten_sw_data
    
    CALL write_flat_data(pressure_scale, droplet_scale)
#else
    CALL read_flat_data(pressure_scale, droplet_scale)
#endif
!$acc enter data copyin(flat_data)

  END SUBROUTINE setup_flat_data

#ifdef PSRAD_WITH_LEGACY

  SUBROUTINE zero_all_dope

    TYPE(dope), PARAMETER :: zero = dope(0,0,0)

    lw_kmajor = zero
    lw_h2oref = zero
    lw_planck = zero
    lw_kgas = zero
    sw_kmajor = zero
    sw_h2oref = zero
    sw_kgas = zero
    sw_rayl = zero
    sw_sfluxref = zero
    lw_info = zero
    sw_info = zero
    lw_band_info = zero
    sw_band_info = zero
    lw_count = 0
    sw_count = 0
  END SUBROUTINE zero_all_dope

  SUBROUTINE measure_lw_data(offset)

    USE mo_psrad_lrtm_kgs, ONLY: minor_species, ngpt, nsp, nh2oref, &
      nsp_fraction

    INTEGER, INTENT(INOUT) :: offset
    INTEGER :: band, atm, which, M, N, i, j, gas, gas2

    lw_info%o = offset
    DO band = 1,nbndlw
      lw_band_info(band)%o = offset
      DO atm = 1,2
        M = nsp_fraction(atm,band)
        IF (M /= 0) THEN
          lw_planck(atm,band) = dope(offset, M, M*ngpt(band))
          offset = offset + lw_planck(atm,band)%n
          lw_count = lw_count + 1
        ENDIF
        M = nsp(atm,band)
        IF (M /= 0) THEN
          lw_kmajor(atm,band) = dope(offset, 5 * M * npressure(atm), &
            5 * M * npressure(atm) * ngpt(band))
          offset = offset + lw_kmajor(atm,band)%n
          lw_count = lw_count + 1
        ENDIF
      ENDDO

      DO which = 1,2
        M = nh2oref(which,band)
        lw_h2oref(which,band) = dope(offset, M, M*ngpt(band))
        offset = offset + lw_h2oref(which,band)%n
        lw_count = lw_count + 1
      END DO

      DO atm = 1,2
        DO i = 1, max_minor_species
          gas = minor_species(i,atm,band)%gas
          IF (gas == 0) EXIT
          IF (gas <= ngas) THEN
            M = minor_species(i,atm,band)%M
            N = minor_species(i,atm,band)%N
            IF (N == 0) THEN
              lw_kgas(i,atm,band) = dope(offset, M, M*ngpt(band))
              offset = offset + lw_kgas(i,atm,band)%n
              lw_count = lw_count + 1
            ELSE
              lw_kgas(i,atm,band) = dope(offset, M*N, M*N*ngpt(band))
              offset = offset + lw_kgas(i,atm,band)%n
              lw_count = lw_count + 1
            ENDIF
          ELSE 
              IF (atm == 2) THEN
                DO j = 1,max_minor_species
                  gas2 = minor_species(j,1,band)%gas
                  IF (gas2 == 0) EXIT
                  IF (gas == gas2) THEN
                    lw_kgas(i,2,band) = lw_kgas(j,1,band)
                    EXIT
                  ENDIF
                ENDDO
                IF (lw_kgas(i,2,band)%n /= 0) CYCLE
              ENDIF
              lw_kgas(i,atm,band) = dope(offset, 1, ngpt(band))
              offset = offset + lw_kgas(i,atm,band)%n
              lw_count = lw_count + 1
          ENDIF
        ENDDO
      ENDDO
    ENDDO
    lw_info%s = offset
  END SUBROUTINE measure_lw_data

  SUBROUTINE measure_sw_data(offset)

    USE mo_psrad_srtm_kgs, ONLY: minor_Species, ngpt, rayl_type, &
      nsfluxref, nsp, nh2oref

    INTEGER, INTENT(INOUT) :: offset
    INTEGER :: band, atm, which, gas, M

    sw_info%o = offset
    DO band = 1,nbndsw
      sw_band_info(band)%o = offset
      DO atm = 1,2
        M = nsp(atm,band)
        IF (M /= 0) THEN
          sw_kmajor(atm,band) = dope(offset, 5 * M * npressure(atm), &
           5 * M * npressure(atm) * ngpt(band))
          offset = offset + sw_kmajor(atm,band)%n
          sw_count = sw_count + 1
        ENDIF
      ENDDO
      DO which = 1,2
        M = nh2oref(which,band) ! for which == 2, may be 3 or 4
        IF (M /= 0) THEN
          !NOTE: Must asways allocate 4 or get oorange.
          IF (M == 3) M = 4
          sw_h2oref(which,band) = dope(offset, M, M * ngpt(band) )
          offset = offset + sw_h2oref(which,band)%n
          sw_count = sw_count + 1
        ENDIF
      ENDDO
      DO atm = 1,2
        gas = minor_species(atm,band)
        IF (gas /= 0) THEN
          sw_kgas(atm,band) = dope(offset, 1, ngpt(band))
          offset = offset + sw_kgas(atm,band)%n
          sw_count = sw_count + 1
        END IF
      ENDDO
      DO atm = 1,2
        M = rayl_type(atm,band)
        IF (M == 0) THEN
          IF (atm == 1) THEN
            sw_rayl(atm,band) = dope(offset,1,1)
            offset = offset + 1
            sw_count = sw_count + 1
          ELSE
            sw_rayl(2,band) = sw_rayl(1,band)
          ENDIF
        ELSEIF (M == 1) THEN
          IF (atm == 1 .or. rayl_type(1,band) > 1) THEN
            sw_rayl(atm,band) = dope(offset, 1, ngpt(band))
            offset = offset + sw_rayl(atm,band)%n
            sw_count = sw_count + 1
          ELSE 
            sw_rayl(2,band) = sw_rayl(1,band)
          ENDIF
        ELSE
          sw_rayl(atm,band) = dope(offset, M, M*ngpt(band))
          offset = offset + sw_rayl(atm,band)%n
          sw_count = sw_count + 1
        ENDIF
      ENDDO
      M = nsfluxref(band)
      sw_sfluxref(band) = dope(offset, M, M*ngpt(band))
      offset = offset + sw_sfluxref(band)%n
      sw_count = sw_count + 1
    ENDDO
    sw_info%s = offset
  END SUBROUTINE measure_sw_data

  SUBROUTINE flatten_lw_data

    USE mo_psrad_legacy_data, ONLY: lw_kmajor_red, lw_h2oref_red, &
      lw_kgas_red, lw_planck_red
    USE mo_psrad_lrtm_kgs, ONLY: minor_species, nsp, nsp_fraction
    INTEGER :: band, atm, which, M, gas
    DO band = 1,nbndlw
      DO atm = 1,2 ! lower/upper
        M = nsp(atm,band)
        IF (M /= 0) THEN
          CALL flatten(lw_kmajor_red(atm,band)%v, flat_data, &
            lw_kmajor(atm,band)%o, lw_kmajor(atm,band)%n)
        ENDIF
      ENDDO
      DO which = 1,2 ! self/foreign
        CALL flatten(lw_h2oref_red(which,band)%v, flat_data, &
          lw_h2oref(which,band)%o, lw_h2oref(which,band)%n)
      ENDDO
      DO atm = 1,2 ! lower/upper
        M = nsp_fraction(atm,band)
        IF (M /= 0) THEN
          CALL flatten(lw_planck_red(atm,band)%v, flat_data, &
            lw_planck(atm,band)%o, lw_planck(atm,band)%n)
        ENDIF
        DO which= 1, max_minor_species
          gas = minor_species(which,atm,band)%gas
          IF (gas == 0) EXIT
          IF (gas <= ngas) THEN
            CALL flatten(lw_kgas_red(which,atm,band)%v, flat_data, &
              lw_kgas(which,atm,band)%o, lw_kgas(which,atm,band)%n)
          ELSEIF (atm == 1) THEN
            CALL flatten(lw_kgas_red(which,1,band)%v, flat_data, &
              lw_kgas(which,1,band)%o, lw_kgas(which,1,band)%n)
          ENDIF
        ENDDO
      ENDDO
    ENDDO
  END SUBROUTINE flatten_lw_data

  SUBROUTINE flatten_sw_data

    USE mo_psrad_legacy_data, ONLY: sw_kmajor_red, sw_h2oref_red, &
      sw_kgas_red, sw_rayl0, sw_rayl1_red, sw_rayl2_red, &
      sw_sfluxref_red
    USE mo_psrad_srtm_kgs, ONLY: nsp, nh2oref, minor_species, &
      rayl_type, nsfluxref
    INTEGER :: band, M, gas, atm, which

    DO band = 1,nbndsw
      DO atm = 1,2 ! lower/upper
        M = nsp(atm,band)
        IF (M /= 0) THEN
          CALL flatten(sw_kmajor_red(atm,band)%v, flat_data, &
            sw_kmajor(atm,band)%o, sw_kmajor(atm,band)%n)
        ENDIF
      END DO
      DO which = 1,2 ! self/foreign
        M = nh2oref(which,band);
        IF (M /= 0) THEN
          CALL flatten(sw_h2oref_red(which,band)%v, flat_data, &
            sw_h2oref(which,band)%o, sw_h2oref(which,band)%n)
        ENDIF
      ENDDO
      DO atm = 1,2 ! lower/upper
        gas = minor_species(atm,band)
        IF (gas /= 0) THEN
          CALL flatten(sw_kgas_red(atm,band)%v, flat_data, &
            sw_kgas(atm,band)%o, sw_kgas(atm,band)%n)
        ENDIF
      END DO
      DO atm = 1,2 ! lower/upper
        M = rayl_type(atm,band)
        IF (M == 0) THEN
          flat_data(sw_rayl(atm,band)%o+1) = sw_rayl0(band)
          IF (atm == 1) EXIT
        ELSEIF (M == 1) THEN
          CALL flatten(sw_rayl1_red(band)%v, flat_data, &
            sw_rayl(atm,band)%o, sw_rayl(atm,band)%n)
          IF (atm == 1) EXIT
        ELSE
          CALL flatten(sw_rayl2_red(atm,band)%v, flat_data, &
            sw_rayl(atm,band)%o, sw_rayl(atm,band)%n)
        ENDIF
      ENDDO
      IF (nsfluxref(band) /= 0) THEN
        CALL flatten(sw_sfluxref_red(band)%v, flat_data, &
          sw_sfluxref(band)%o, sw_sfluxref(band)%n)
      ENDIF
    ENDDO
  END SUBROUTINE flatten_sw_data

  SUBROUTINE flatten1(src, tgt, o, s)
    REAL(wp), INTENT(IN) :: src(:)
    REAL(wp), INTENT(INOUT) :: tgt(:)
    INTEGER, INTENT(IN) :: o, s
    INTEGER :: n1
    n1 = size(src,1)
    IF (n1 /= s) THEN
      WRITE(*,'(a,i3,a,i5)') 'Size mismatch in flatten2: ', &
        n1, ' /= ', s
    ENDIF
    tgt(o+1:o+s) = src(1:s)
  END SUBROUTINE flatten1

  SUBROUTINE flatten2(src, tgt, o, s)
    REAL(wp), INTENT(IN) :: src(:,:)
    REAL(wp), INTENT(INOUT) :: tgt(:)
    INTEGER, INTENT(IN) :: o, s
    INTEGER :: n1, n2
    n1 = size(src,1)
    n2 = size(src,2)
    IF (n1 * n2 /= s) THEN
      WRITE(*,'(a,i4,a,i4,a,i5)') 'Size mismatch in flatten2: ', &
        n1, ' * ', n2, ' /= ', s
    ENDIF
    tgt(o+1:o+s) = RESHAPE(src, SHAPE = (/n1 * n2/))
  END SUBROUTINE flatten2

  SUBROUTINE flatten3(src, tgt, o, s)
    REAL(wp), INTENT(IN) :: src(:,:,:)
    REAL(wp), INTENT(INOUT) :: tgt(:)
    INTEGER, INTENT(IN) :: o, s
    INTEGER :: n1, n2, n3
    n1 = size(src,1)
    n2 = size(src,2)
    n3 = size(src,3)
    IF (n1 * n2 * n3 /= s) THEN
      WRITE(*,'(a,i4,a,i4,a,i4,a,i5)') 'Size mismatch in flatten3: ', &
        n1, ' * ', n2, ' * ', n3, ' /= ', s
    ENDIF
    tgt(o+1:o+s) = RESHAPE(src, SHAPE = (/n1 * n2 * n3/))
  END SUBROUTINE flatten3

  SUBROUTINE flatten4(src, tgt, o, s)
    REAL(wp), INTENT(IN) :: src(:,:,:,:)
    REAL(wp), INTENT(INOUT) :: tgt(:)
    INTEGER, INTENT(IN) :: o, s
    INTEGER :: n1, n2, n3, n4
    n1 = size(src,1)
    n2 = size(src,2)
    n3 = size(src,3)
    n4 = size(src,4)
    IF (n1 * n2 * n3 * n4 /= s) THEN
      WRITE(*,'(a,i4,a,i4,a,i4,a,i4,a,i5)') 'Size mismatch in flatten3: ', &
        n1, ' * ', n2, ' * ', n3, ' * ', n4, ' /= ', s
    ENDIF
    tgt(o+1:o+s) = RESHAPE(src, SHAPE = (/n1 * n2 * n3 * n4/))
  END SUBROUTINE flatten4
#endif

  SUBROUTINE count_dope(n_dope)
    INTEGER, INTENT(OUT) :: n_dope
    n_dope = 2 + &
      nbndlw * (2 * (3 + max_minor_species) + 1 ) + &
      nbndsw * (2 * 4 + 2)
  END SUBROUTINE count_dope

#ifdef PSRAD_WITH_LEGACY

  SUBROUTINE write_dope(d, offset, o, s, n)
    TYPE(dope), INTENT(IN) :: d
    INTEGER, INTENT(INOUT) :: offset
    INTEGER, INTENT(INOUT) :: o(:), s(:), n(:)
    o(offset) = d%o
    s(offset) = d%s
    n(offset) = d%n
    offset = offset + 1
  END SUBROUTINE write_dope

  SUBROUTINE pack_dope(o, s, n)
    INTEGER, INTENT(OUT) :: o(:), s(:), n(:)
    INTEGER :: i, j, k, offset

    offset = 1
    CALL write_dope(lw_info, offset, o, s, n)
    CALL write_dope(sw_info, offset, o, s, n)
    DO j = 1, nbndlw
      DO i = 1, 2
        CALL write_dope(lw_kmajor(i,j), offset, o, s, n)
        CALL write_dope(lw_h2oref(i,j), offset, o, s, n)
        CALL write_dope(lw_planck(i,j), offset, o, s, n)
        DO k = 1, max_minor_species
          CALL write_dope(lw_kgas(k,i,j), offset, o, s, n)
        ENDDO
      ENDDO
      CALL write_dope(lw_band_info(j), offset, o, s, n)
    ENDDO
    DO j = 1, nbndsw
      DO i = 1, 2
        CALL write_dope(sw_kmajor(i,j), offset, o, s, n)
        CALL write_dope(sw_h2oref(i,j), offset, o, s, n)
        CALL write_dope(sw_kgas(i,j), offset, o, s, n)
        CALL write_dope(sw_rayl(i,j), offset, o, s, n)
      ENDDO
      CALL write_dope(sw_sfluxref(j), offset, o, s, n)
      CALL write_dope(sw_band_info(j), offset, o, s, n)
    ENDDO
  END SUBROUTINE pack_dope

  SUBROUTINE write_flat_data(pressure_scale, droplet_scale)
    USE mo_psrad_lrtm_kgs, ONLY: lrtm_count_flat, lrtm_pack_data
    REAL(wp), INTENT(IN) :: pressure_scale, droplet_scale

    INTEGER :: n_dope 
    INTEGER, ALLOCATABLE :: flat_dope_o(:), flat_dope_s(:), flat_dope_n(:)
    INTEGER :: n_lrtm
    REAL(wp), ALLOCATABLE :: flat_lrtm(:)
    INTEGER, PARAMETER :: n_integer = 2, n_real = 2
    INTEGER :: flat_integer(n_integer), ret, idummy, fid
    REAL(wp) :: flat_real(n_real)
    INTEGER :: idx_real, idx_integer, idx_lrtm, idx_flat_data, &
      idx_dope_o, idx_dope_s, idx_dope_n

    flat_integer = (/lw_count, sw_count/)
    flat_real = (/pressure_scale, droplet_scale/)
    
    CALL count_dope(n_dope)
    ALLOCATE(flat_dope_o(n_dope))
    ALLOCATE(flat_dope_s(n_dope))
    ALLOCATE(flat_dope_n(n_dope))
    CALL pack_dope(flat_dope_o, flat_dope_s, flat_dope_n)

    CALL lrtm_count_flat(n_lrtm)
    ALLOCATE(flat_lrtm(n_lrtm))
    CALL lrtm_pack_data(flat_lrtm)

    ret = nf_create('lsdata.nc', nf_clobber, fid)

    ret = nf_def_dim(fid, 'n_real', n_real, idummy)
    ret = nf_def_dim(fid, 'n_integer', n_integer, idummy)
    ret = nf_def_dim(fid, 'n_lrtm', n_lrtm, idummy)
    ret = nf_def_dim(fid, 'n_dope', n_dope, idummy)
    ret = nf_def_dim(fid, 'flat_data_size', flat_data_size, idummy)

    ret = nf_def_var(fid, 'flat_real', nf_double, 1, 1, idx_real)
    ret = nf_def_var(fid, 'flat_integer', nf_int, 1, 2, idx_integer)
    ret = nf_def_var(fid, 'flat_lrtm', nf_double, 1, 3, idx_lrtm)
    ret = nf_def_var(fid, 'flat_dope_o', nf_int, 1, 4, idx_dope_o)
    ret = nf_def_var(fid, 'flat_dope_s', nf_int, 1, 4, idx_dope_s)
    ret = nf_def_var(fid, 'flat_dope_n', nf_int, 1, 4, idx_dope_n)
    ret = nf_def_var(fid, 'flat_data', nf_double, 1, 5, idx_flat_data)

    ret = nf_enddef(fid)

    ret = nf_put_vara_double(fid, idx_real, [1], [n_real], flat_real)
    ret = nf_put_vara_int(fid, idx_integer, [1], [n_integer], flat_integer)
    ret = nf_put_vara_double(fid, idx_lrtm, [1], [n_lrtm], flat_lrtm)
    ret = nf_put_vara_int(fid, idx_dope_o, [1], [n_dope], flat_dope_o)
    ret = nf_put_vara_int(fid, idx_dope_s, [1], [n_dope], flat_dope_s)
    ret = nf_put_vara_int(fid, idx_dope_n, [1], [n_dope], flat_dope_n)
    ret = nf_put_vara_double(fid, idx_flat_data, [1], [flat_data_size], &
      flat_data)

    ret = nf_close(fid)
    
  END SUBROUTINE write_flat_data

#else

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
!FIXME: workaround for buggy library on Aurora testbed
#ifndef __NEC__
      CALL finish('psrad_flat_data','Initialization value for pressure scale incompatible with the one used for generating data.')
#else
      WRITE(*,*) 'psrad_flat_data: Initialization value for pressure scale incompatible with the one used for generating data.'
      STOP
#endif
    ENDIF
    IF ((droplet_scale - flat_real(2)) / flat_real(2) > 1e-14) THEN
!FIXME: workaround for buggy library on Aurora testbed
#ifndef __NEC__
      CALL finish('psrad_flat_data','Initialization value for droplet scale incompatible with the one used for generating data.')
#else
      WRITE (*,*) 'psrad_flat_data: Initialization value for droplet scale incompatible with the one used for generating data.'
      STOP
#endif
    ENDIF
    CALL unpack_dope(flat_dope_o, flat_dope_s, flat_dope_n)
    CALL lrtm_unpack_data(flat_lrtm)

  END SUBROUTINE read_flat_data

#endif

END MODULE mo_psrad_flat_data

