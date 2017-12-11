!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_psrad_srtm_netcdf

  USE mo_psrad_general, ONLY : wp, finish, nbndsw, ngas, &
    ngpt_orig, npressure, ptr1, ptr2, ptr3, ptr4
  USE mo_psrad_io, ONLY: psrad_io_open, psrad_io_close, &
    read=>psrad_io_copy_double
  USE mo_psrad_srtm_kgs, ONLY : nh2oref, nsp, nsfluxref, rayl_type, &
    rayl0, minor_species, minor_species_missing_data


  IMPLICIT NONE

  PRIVATE
  PUBLIC :: srtm_read

  INTEGER, PARAMETER :: Tdiff = 5
  !matching psrad_general (IDENTICAL to lrtm_netcdf...)
  ! ih2o  = 1, ico2  = 2, ich4  = 3, io2   = 4, 
  ! io3   = 5, in2o = 6, ico = 7, in2 = 8, ngas = 8
  INTEGER, PARAMETER :: netcdf_idx(ngas) = (/ &
    6, 7, 11, 12, 8, 9, 10, 1/)

  INTEGER, PARAMETER :: gPointSetNumber = 1
  INTEGER :: fileid     !< id number of netcdf file
  REAL(wp), ALLOCATABLE :: transposed_data(:,:)

CONTAINS 

  SUBROUTINE srtm_read(kmajor3_in, kmajor4_in, h2oref_in, kgas_in, &
      rayl1_in, rayl2_in, sfluxref_in)
    TYPE(ptr3), DIMENSION(2,nbndsw), INTENT(INOUT) :: kmajor3_in
    TYPE(ptr4), DIMENSION(2,nbndsw), INTENT(INOUT) :: kmajor4_in
    TYPE(ptr2), DIMENSION(2,nbndsw), INTENT(INOUT) :: h2oref_in
    TYPE(ptr1), DIMENSION(2,nbndsw), INTENT(INOUT) :: kgas_in
    TYPE(ptr1), INTENT(INOUT) :: rayl1_in(nbndsw)
    TYPE(ptr2), INTENT(INOUT) :: rayl2_in(2,nbndsw)
    TYPE(ptr2), DIMENSION(nbndsw), INTENT(INOUT) :: sfluxref_in

    CALL psrad_io_open('rrtmg_sw.nc', fileid)
    IF (fileid == 0) THEN
      CALL finish('mo_psrad_srtm_netcdf/srtm_read', &
        'File rrtmg_sw.nc cannot be opened')
    END IF

    ALLOCATE(transposed_data(ngpt_orig,9))
    CALL read_key_species(kmajor3_in, kmajor4_in)
    CALL read_water(h2oref_in)
    CALL read_gases(kgas_in)
    CALL read_solar_source(sfluxref_in)
    CALL read_rayleigh(rayl1_in, rayl2_in)
    DEALLOCATE(transposed_data)

    CALL psrad_io_close(fileid)

  END SUBROUTINE srtm_read

  SUBROUTINE read_key_species(kmajor3_in, kmajor4_in)
    TYPE(ptr3), DIMENSION(2,nbndsw), INTENT(INOUT) :: kmajor3_in
    TYPE(ptr4), DIMENSION(2,nbndsw), INTENT(INOUT) :: kmajor4_in

    CHARACTER(len=*), PARAMETER :: var_name(2) = (/ &
      'KeySpeciesAbsorptionCoefficientsLowerAtmos', &
      'KeySpeciesAbsorptionCoefficientsUpperAtmos'/)
    INTEGER :: band, atm, n
    DO band = 1, nbndsw
    DO atm = 1,2
      n = nsp(atm,band)
      IF (n /= 0) THEN
        IF (n == 1) THEN
          CALL read(fileid, var_name(atm), &
            (/1,1,1,1,band,gPointSetNumber/), &
            (/1,Tdiff,npressure(atm),ngpt_orig,1,1/), &
            kmajor3_in(atm,band)%v)
        ELSE
          CALL read(fileid, var_name(atm), &
            (/1,1,1,1,band,gPointSetNumber/), &
            (/n,Tdiff,npressure(atm),ngpt_orig,1,1/), &
            kmajor4_in(atm,band)%v)
        ENDIF
      ENDIF 
    ENDDO
    ENDDO
  END SUBROUTINE read_key_species

  SUBROUTINE read_water(h2oref_in)
    TYPE(ptr2), DIMENSION(2,nbndsw), INTENT(INOUT) :: h2oref_in
    CHARACTER(len=*), PARAMETER :: var_name(3) = (/ &
      'H2OSelfAbsorptionCoefficients             ', &
      'H2OForeignAbsorptionCoefficientsLowerAtmos', &
      'H2OForeignAbsorptionCoefficientsUpperAtmos'/)
    INTEGER :: band, which, n, hack_for_split_data

    hack_for_split_data = 0
    DO band = 1,nbndsw
      DO which = 1,2
        n = nh2oref(which,band)
        IF (n /= 0) THEN
          ! IF (which == 1) THEN ! n == 10
          IF (n == 4) THEN
            hack_for_split_data = 1
            n = 3
          END IF
          CALL read(fileid,trim(var_name(which)), &
            (/1,1,band,gPointSetNumber/), &
            (/n,ngpt_orig,1,1/), &
            h2oref_in(which,band)%v(1:n,:))
        ENDIF
        ! TODO: fix input data! Notice:
        ! the 2 in the first position of the 2nd argument
        IF (hack_for_split_data == 1) THEN
          CALL read(fileid,trim(var_name(3)), &
            (/2,1,band,gPointSetNumber/), &
            (/1,ngpt_orig,1,1/), &
            h2oref_in(2,band)%v(4,:))
          hack_for_split_data =0
        ENDIF
      ENDDO
    ENDDO
  END SUBROUTINE read_water

  SUBROUTINE read_gases(kgas_in)
    TYPE(ptr1), DIMENSION(2,nbndsw), INTENT(INOUT) :: kgas_in
    CHARACTER(len=*), PARAMETER :: var_name(2) = (/ &
      'AbsorptionCoefficientsLowerAtmos', &
      'AbsorptionCoefficientsUpperAtmos'/)
    INTEGER :: band, atm, gas
    DO band = 1,nbndsw
      DO atm = 1,2
        gas = minor_species(atm,band)
        IF (gas /= 0) THEN
          IF (minor_species_missing_data(atm,band) == 0 .and. &
            band /= 7) THEN
            CALL read(fileid,var_name(atm), &
              (/1,1,1,netcdf_idx(gas),band,gPointSetNumber/), &
              (/1,1,ngpt_orig,1,1,1/), &
              kgas_in(atm,band)%v)
          ENDIF
        ENDIF
      ENDDO
    ENDDO
  END SUBROUTINE read_gases


  SUBROUTINE read_solar_source(sfluxref_in)
    TYPE(ptr2), DIMENSION(nbndsw), INTENT(INOUT) :: sfluxref_in
    CHARACTER(len=*), PARAMETER :: var_name(2) = (/ &
      'SolarSourceFunctionLowerAtmos', &
      'SolarSourceFunctionUpperAtmos'/)
    INTEGER, PARAMETER :: which_source_variable(nbndsw) = (/ &
      1,2,1,1,1, 1,1,1,1,1, 1,1,2,1/)
    INTEGER :: band, m, which

    DO band = 1,nbndsw
      m = nsfluxref(band)
      which = which_source_variable(band)
      CALL read(fileid, var_name(which), &
        (/1,1,band,gPointSetNumber/), &
        (/ngpt_orig,m,1,1/), transposed_data(:,1:m))
        sfluxref_in(band)%v = TRANSPOSE(transposed_data(:,1:m))
    ENDDO
  END SUBROUTINE read_solar_source

  SUBROUTINE read_rayleigh(rayl1_in, rayl2_in)
    TYPE(ptr1), INTENT(INOUT) :: rayl1_in(nbndsw)
    TYPE(ptr2), INTENT(INOUT) :: rayl2_in(2,nbndsw)
    CHARACTER(len=*), PARAMETER :: var_name(2) = (/ &
      'RayleighExtinctionCoefficientsLowerAtmos', &
      'RayleighExtinctionCoefficientsUpperAtmos'/)
    INTEGER :: band, atm, m
    DO band = 1,nbndsw
      DO atm = 1,2
        m = rayl_type(atm,band)
        IF (m == 0) THEN
          CALL read(fileid, var_name(atm), &
            (/1,1,band,gPointSetNumber/), &
            (/1,1,1,1/), rayl0(band))
        ELSEIF (m == 1) THEN
          CALL read(fileid, var_name(atm), &
            (/1,1,band,gPointSetNumber/), &
            (/ngpt_orig,1,1,1/), rayl1_in(band)%v)
        ELSE
          CALL read(fileid, var_name(atm), &
            (/1,1,band,gPointSetNumber/), &
            (/ngpt_orig,m,1,1/), transposed_data(:,1:m))

            rayl2_in(atm,band)%v = TRANSPOSE(transposed_data(:,1:m))
            CYCLE
        ENDIF
        !NOTE: all bands share rayl between upper/lower, except 9th
        EXIT 
      END DO
    END DO
  END SUBROUTINE read_rayleigh

END MODULE mo_psrad_srtm_netcdf

