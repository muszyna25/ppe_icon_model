!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_psrad_lrtm_netcdf

  USE mo_psrad_general, ONLY: wp, finish, ngas, ncfc, nbndlw, &
    ih2o, ico2, io3, in2o, io2, ich4, ico, ngpt_orig, npressure, &
    cfc_offset, max_minor_species, ten20inv, ptr2, ptr3, ptr4
  USE mo_psrad_io, ONLY: psrad_io_open, psrad_io_close, &
    read=>psrad_io_copy_double
  USE mo_psrad_lrtm_kgs, ONLY : nsp, nsp_fraction, minor_species

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: lrtm_read

  INTEGER, PARAMETER :: keylower = 9, keyupper = 5, Tdiff = 5, ps = 59, &
    plower = 13, pupper = 47, Tself = 10, Tforeign = 4, pforeign = 4, &
    T = 19, Tplanck = 181, GPoint = 16, GPointSet = 2

  !matching psrad_general
  ! ih2o  = 1, ico2  = 2, ich4  = 3, io2   = 4, 
  ! io3   = 5, in2o = 6, ico = 7, in2 = 8, ngas = 8
  INTEGER, PARAMETER :: netcdf_idx(ngas) = (/ &
    6, 7, 11, 12, 8, 9, 10, 1/)
  !  iccl4 = cfc_offset+1, icfc11 = cfc_offset+2, &
  !  icfc12 = cfc_offset+3, icfc22 = cfc_offset+4, &
  INTEGER, PARAMETER :: netcdf_cfc_idx(ncfc) = (/ &
    2, 3, 4, 5/)

  INTEGER, PARAMETER :: gPointSetNumber = 1
  INTEGER :: fileid

CONTAINS 

  SUBROUTINE lrtm_read(kmajor3_in, kmajor4_in, h2oref_in, kgas_in, &
    planck_fraction2_in)

    USE mo_psrad_lrtm_kgs, ONLY: chi_mls, totplanck, totplanck16

    TYPE(ptr3), DIMENSION(2,nbndlw), INTENT(INOUT) :: kmajor3_in
    TYPE(ptr4), DIMENSION(2,nbndlw), INTENT(INOUT) :: kmajor4_in
    TYPE(ptr2), DIMENSION(2,nbndlw), INTENT(INOUT) :: h2oref_in
    TYPE(ptr3), DIMENSION(max_minor_species,2,nbndlw), &
      INTENT(INOUT) :: kgas_in
    TYPE(ptr2), DIMENSION(2,nbndlw), INTENT(INOUT) :: planck_fraction2_in


    CALL psrad_io_open('rrtmg_lw.nc', fileid)
    IF (fileid == 0) THEN
      CALL finish('mo_psrad_lrtm_netcdf/lrtm_read', 'File rrtmg_lw.nc cannot be opened')
    END IF

    chi_mls = 0
    CALL read(fileid, 'AbsorberAmountMLS', &
      (/6, 1/), &
      (/7, 59 /), & ! TODO
      chi_mls(1:7,:))
    !NOTE: scaling data instead of lots of runtime scaling/descalings
    chi_mls((/ih2o, ico2, io3, in2o, io2, ich4, ico/),:) = &
      ten20inv * chi_mls(1:7,:)

    CALL read(fileid, 'IntegratedPlanckFunction', &
      (/1, 1/), &
      (/SIZE(totplanck,1), SIZE(totplanck,2)/), &
      totplanck)
    CALL read(fileid, 'IntegratedPlanckFunctionBand16', &
      (/ 1 /), &
      (/SIZE(totplanck16)/), &
      totplanck16)
    CALL read_gases(kgas_in)
    CALL read_water(h2oref_in)
    CALL read_planck_fraction(planck_fraction2_in)
    CALL read_major_species(kmajor3_in, kmajor4_in)

    CALL psrad_io_close(fileid)

  END SUBROUTINE lrtm_read 

  SUBROUTINE read_major_species(kmajor3_in, kmajor4_in)
    TYPE(ptr3), DIMENSION(2,nbndlw), INTENT(INOUT) :: kmajor3_in
    TYPE(ptr4), DIMENSION(2,nbndlw), INTENT(INOUT) :: kmajor4_in

    CHARACTER(len=*), PARAMETER :: var_name(2) = (/ &
      'KeySpeciesAbsorptionCoefficientsLowerAtmos', &
      'KeySpeciesAbsorptionCoefficientsUpperAtmos'/)
    INTEGER :: band, atm, n
    DO band = 1, nbndlw
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
  END SUBROUTINE read_major_species

  SUBROUTINE read_water(h2oref_in)
    TYPE(ptr2), DIMENSION(2,nbndlw), INTENT(INOUT) :: h2oref_in
    CHARACTER(len=*), PARAMETER :: var_name(2) = (/ &
      'H20SelfAbsorptionCoefficients   ', &
      'H20ForeignAbsorptionCoefficients'/)
    INTEGER, PARAMETER :: key(2) = (/10, 4/)
    INTEGER :: band, which
    DO band = 1,nbndlw
    DO which = 1,2
      CALL read(fileid,trim(var_name(which)), &
        (/1,1,band,gPointSetNumber/), &
        (/key(which),ngpt_orig,1,1/), &
        h2oref_in(which,band)%v)
    ENDDO
    ENDDO
  END SUBROUTINE read_water

  SUBROUTINE read_gases(kgas_in)
    TYPE(ptr3), DIMENSION(max_minor_species,2,nbndlw), &
      INTENT(INOUT) :: kgas_in
    CHARACTER(len=*), PARAMETER :: var_name(2) = (/ &
      'AbsorptionCoefficientsLowerAtmos', &
      'AbsorptionCoefficientsUpperAtmos'/)
    INTEGER, PARAMETER :: key(2) = (/9, 5/)
    INTEGER :: band, atm, i, gas, N
    DO band = 1,nbndlw
      DO atm = 1,2
        DO i = 1, max_minor_species
          gas = minor_species(i,atm,band)%gas
          IF (gas == 0) EXIT
          IF (gas <= ngas) THEN
            N = minor_species(i,atm,band)%N
            IF (N == 0) THEN
              CALL read(fileid,var_name(atm), &
                (/1,1,1,netcdf_idx(gas),band,gPointSetNumber/), &
                (/1,T,ngpt_orig,1,1,1/), &
                kgas_in(i,atm,band)%v(:,:,1))
            ELSE
              CALL read(fileid,var_name(atm), &
                (/1,1,1,netcdf_idx(gas),band,gPointSetNumber/), &
                (/key(atm),T,ngpt_orig,1,1,1/), &
                kgas_in(i,atm,band)%v)
            ENDIF
          ELSEIF (atm == 1) THEN
            gas = gas - cfc_offset
            CALL read(fileid,var_name(1), &
              (/1,1,1,netcdf_cfc_idx(gas),band,gPointSetNumber/), &
              (/1,1,ngpt_orig,1,1,1/), &
              kgas_in(i,1,band)%v(:,1,1))
          ENDIF
        ENDDO
      ENDDO
    ENDDO
  END SUBROUTINE read_gases

  SUBROUTINE read_planck_fraction(planck_fraction2_in)
    TYPE(ptr2), DIMENSION(2,nbndlw), INTENT(INOUT) :: planck_fraction2_in
    CHARACTER(len=*), PARAMETER :: var_name(2) = (/ &
      'PlanckFractionLowerAtmos', &
      'PlanckFractionUpperAtmos'/)
    INTEGER :: band, atm, n
    REAL(wp) :: linear(ngpt_orig*MAXVAL(nsp_fraction))
    DO band = 1, nbndlw
    DO atm = 1,2
      n = nsp_fraction(atm,band)
      IF (n /= 0) THEN
        IF (n == 1) THEN
          CALL read(fileid,var_name(atm), &
            (/1,1,band,gPointSetNumber/), &
            (/ngpt_orig,1,1,1/), &
            planck_fraction2_in(atm,band)%v(1,:))
        ELSE
          CALL read(fileid,var_name(atm), &
            (/1,1,band,gPointSetNumber/), &
            (/ngpt_orig,n,1,1/), &
            linear)
            planck_fraction2_in(atm,band)%v(1:n,:) = TRANSPOSE(RESHAPE(&
              linear, SHAPE =(/ngpt_orig,n/)))
        ENDIF
      ENDIF
    ENDDO
    ENDDO
  END SUBROUTINE read_planck_fraction

END MODULE mo_psrad_lrtm_netcdf
