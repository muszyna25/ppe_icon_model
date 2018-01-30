!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_psrad_lrtm_netcdf

  USE mo_psrad_lrtm_kgs, ONLY : no, kmajor3_in, kmajor4_in, &
    planck_fraction1_in, planck_fraction2_in, &
    kgas2_in, kgas3_in, kgas2_list, kgas3_list, h2oref_in, cfc_in, &
    npressure, nsp_species, nsp_fraction, cfc_list
  USE mo_psrad_general, ONLY: wp, finish, ngas, ncfc, nbndlw, &
    ih2o, ico2, io3, in2o, io2, ich4, ico!, in2
  USE mo_psrad_io, ONLY: psrad_io_open, &
    read=>psrad_io_copy_double

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

  SUBROUTINE read_gases
    CHARACTER(len=*), PARAMETER :: var_name(2) = (/ &
      'AbsorptionCoefficientsLowerAtmos', &
      'AbsorptionCoefficientsUpperAtmos'/)
    INTEGER, PARAMETER :: key(2) = (/9, 5/)
    INTEGER :: band, atm, igas, M
    DO band = 1,nbndlw
      DO igas = 1, ncfc
        IF (cfc_list(igas,band) /= 0) THEN
          CALL read(fileid,var_name(1), &
            (/1,1,1,netcdf_cfc_idx(igas),band,gPointSetNumber/), &
            (/1,1,no,1,1,1/), &
            cfc_in(igas,band)%v)
        END IF
      END DO
      DO atm = 1,2
      DO igas = 1, ngas
        M = kgas2_list(igas,atm,band)
        IF (M /= 0) THEN
          CALL read(fileid,var_name(atm), &
            (/1,1,1,netcdf_idx(igas),band,gPointSetNumber/), &
            (/1,T,no,1,1,1/), &
            kgas2_in(igas,atm,band)%v)
        ENDIF
        M = kgas3_list(1,igas,atm,band)
        IF (M /= 0) THEN
          CALL read(fileid,var_name(atm), &
            (/1,1,1,netcdf_idx(igas),band,gPointSetNumber/), &
            (/key(atm),T,no,1,1,1/), &
            kgas3_in(igas,atm,band)%v)
        ENDIF
      ENDDO
      ENDDO
    ENDDO
  END SUBROUTINE read_gases

  SUBROUTINE read_water
    CHARACTER(len=*), PARAMETER :: var_name(2) = (/ &
      'H20SelfAbsorptionCoefficients   ', &
      'H20ForeignAbsorptionCoefficients'/)
    INTEGER, PARAMETER :: key(2) = (/10, 4/)
    INTEGER :: band, which
    DO band = 1,nbndlw
    DO which = 1,2
      CALL read(fileid,trim(var_name(which)), &
        (/1,1,band,gPointSetNumber/), &
        (/key(which),no,1,1/), &
        h2oref_in(1:key(which),:,which,band))
    ENDDO
    ENDDO
  END SUBROUTINE read_water

  SUBROUTINE read_planck_fraction
    CHARACTER(len=*), PARAMETER :: var_name(2) = (/ &
      'PlanckFractionLowerAtmos', &
      'PlanckFractionUpperAtmos'/)
    INTEGER :: band, atm, n
    REAL(wp) :: linear(no*MAXVAL(nsp_fraction))
    DO band = 1, nbndlw
    DO atm = 1,2
      n = nsp_fraction(atm,band)
      IF (n /= 0) THEN
        IF (n == 1) THEN
          CALL read(fileid,var_name(atm), &
            (/1,1,band,gPointSetNumber/), &
            (/no,1,1,1/), &
            planck_fraction1_in(:,atm,band))
        ELSE
          CALL read(fileid,var_name(atm), &
            (/1,1,band,gPointSetNumber/), &
            (/no,n,1,1/), &
            linear)
            planck_fraction2_in(atm,band)%v(1:n,:) = TRANSPOSE(RESHAPE(&
              linear, SHAPE =(/no,n/)))
        ENDIF
      ENDIF
    ENDDO
    ENDDO
  END SUBROUTINE read_planck_fraction

  SUBROUTINE read_key_species
    CHARACTER(len=*), PARAMETER :: var_name(2) = (/ &
      'KeySpeciesAbsorptionCoefficientsLowerAtmos', &
      'KeySpeciesAbsorptionCoefficientsUpperAtmos'/)
    INTEGER :: band, atm, n
    DO band = 1, nbndlw
    DO atm = 1,2
      n = nsp_species(atm,band)
      IF (n /= 0) THEN
        IF (n == 1) THEN
          CALL read(fileid, var_name(atm), &
            (/1,1,1,1,band,gPointSetNumber/), &
            (/1,Tdiff,npressure(atm),no,1,1/), &
            kmajor3_in(atm,band)%v)
        ELSE
          CALL read(fileid, var_name(atm), &
            (/1,1,1,1,band,gPointSetNumber/), &
            (/n,Tdiff,npressure(atm),no,1,1/), &
            kmajor4_in(atm,band)%v)
        ENDIF
      ENDIF 
    ENDDO
    ENDDO
  END SUBROUTINE read_key_species


  SUBROUTINE lrtm_read

    USE mo_psrad_lrtm_kgs, ONLY: chi_mls, totplanck, totplanck16

    CALL psrad_io_open('rrtmg_lw.nc', fileid)
    IF (fileid == 0) THEN
      CALL finish('mo_psrad_lrtm_netcdf/lrtm_read', 'File rrtmg_lw.nc cannot be opened')
    END IF

    chi_mls = 0
    CALL read(fileid, 'AbsorberAmountMLS', &
      (/6, 1/), &
      (/7, 59 /), & ! TODO
      chi_mls(1:7,:))
    chi_mls((/ih2o, ico2, io3, in2o, io2, ich4, ico/),:) = chi_mls(1:7,:)

    CALL read(fileid, 'IntegratedPlanckFunction', &
      (/1, 1/), &
      (/SIZE(totplanck,1), SIZE(totplanck,2)/), &
      totplanck)
    CALL read(fileid, 'IntegratedPlanckFunctionBand16', &
      (/ 1 /), &
      (/SIZE(totplanck16)/), &
      totplanck16)
    CALL read_gases
    CALL read_water
    CALL read_planck_fraction
    CALL read_key_species

  END SUBROUTINE lrtm_read 

END MODULE mo_psrad_lrtm_netcdf
