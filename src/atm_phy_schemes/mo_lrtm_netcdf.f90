!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.

MODULE mo_lrtm_netcdf

  USE mo_exception,            ONLY: finish

  USE mo_netcdf_parallel, ONLY: p_nf_open, p_nf_close, &
    &                           p_nf_inq_varid,        &
    &                           p_nf_get_vara_double,  &
    &                           nf_read, nf_noerr

  IMPLICIT NONE

  PRIVATE


  PUBLIC :: lrtm_read,                              &
    &       lw_kgb01, lw_kgb02, lw_kgb03, lw_kgb04, &
    &       lw_kgb05, lw_kgb06, lw_kgb07, lw_kgb08, &
    &       lw_kgb09, lw_kgb10, lw_kgb11, lw_kgb12, &
    &       lw_kgb13, lw_kgb14, lw_kgb15, lw_kgb16

  INTEGER, PARAMETER :: maxAbsorberNameLength =  5, &
    &                   Absorber              = 12

  CHARACTER(len = maxAbsorberNameLength), PARAMETER :: &
    &  AbsorberNames(Absorber) = (/'N2   ',            &
    &                              'CCL4 ',            &
    &                              'CFC11',            &
    &                              'CFC12',            &
    &                              'CFC22',            &
    &                              'H2O  ',            &
    &                              'CO2  ',            &
    &                              'O3   ',            &
    &                              'N2O  ',            &
    &                              'CO   ',            &
    &                              'CH4  ',            &
    &                              'O2   '/)

  INTEGER, PARAMETER :: &
    &  keylower  = 9,   &
    &  keyupper  = 5,   &
    &  Tdiff     = 5,   &
    &  ps        = 59,  &
    &  plower    = 13,  &
    &  pupper    = 47,  &
    &  Tself     = 10,  &
    &  Tforeign  = 4,   &
    &  pforeign  = 4,   &
    &  T         = 19,  &
    &  Tplanck   = 181, &
    &  band      = 16,  &
    &  GPoint    = 16,  &
    &  GPointSet = 2

  INTEGER, PARAMETER :: gPointSetNumber = 1

  INTEGER :: fileid     !< id number of netcdf file
  INTEGER :: varid      !< id number of variable in netcdf file
  INTEGER :: nf_status  !< return status of netcdf function

CONTAINS

  !=============================================================================

  INTEGER FUNCTION AbsorberIndex(AbsorberName)
    CHARACTER(len=*), INTENT(in)  :: AbsorberName

    INTEGER :: m

    AbsorberIndex = -1
    DO m = 1, Absorber
      IF (TRIM(AbsorberNames(m)) == TRIM(AbsorberName)) THEN
        AbsorberIndex = m
      END IF
    END DO

    IF (AbsorberIndex == -1) THEN
      CALL finish('Absorber name index lookup failed.')
    END IF

  END FUNCTION AbsorberIndex

  !=============================================================================

  SUBROUTINE lrtm_read(data_filename)

    !> NetCDF file containing longwave absorption coefficients and other data
    !> for RRTMG_LW k-distribution model
    CHARACTER (LEN=*), INTENT(IN) :: data_filename

    nf_status = p_nf_open(TRIM(data_filename), nf_read, fileid)

    IF (nf_status /= nf_noerr) THEN
      CALL finish('mo_lrtm_netcdf/lrtm_read', 'File '//TRIM(data_filename)//' cannot be opened')
    END IF

    CALL lw_kgb01  ! molecular absorption coefficients
    CALL lw_kgb02
    CALL lw_kgb03
    CALL lw_kgb04
    CALL lw_kgb05
    CALL lw_kgb06
    CALL lw_kgb07
    CALL lw_kgb08
    CALL lw_kgb09
    CALL lw_kgb10
    CALL lw_kgb11
    CALL lw_kgb12
    CALL lw_kgb13
    CALL lw_kgb14
    CALL lw_kgb15
    CALL lw_kgb16

    nf_status = p_nf_close(fileid)

  END SUBROUTINE lrtm_read

  SUBROUTINE lw_kgb01

    USE mo_rrlw_kg01, ONLY : fracrefao, fracrefbo, kao, kbo, kao_mn2, kbo_mn2, selfrefo, forrefo, no1

    INTEGER, PARAMETER :: bandNumber = 1
    INTEGER, PARAMETER :: numGPoints = no1


    nf_status = p_nf_inq_varid(fileid, 'PlanckFractionLowerAtmos', varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                      &
      &                              (/1,1,bandNumber,gPointSetNumber/), &
      &                              (/numGPoints,1,1,1/),               &
      &                              fracrefao)

    nf_status = p_nf_inq_varid(fileid,'PlanckFractionUpperAtmos',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                      &
      &                              (/1,1,bandNumber,gPointSetNumber/), &
      &                              (/numGPoints,1,1,1/),               &
      &                              fracrefbo)

    nf_status = p_nf_inq_varid(fileid,'KeySpeciesAbsorptionCoefficientsLowerAtmos', varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                          &
      &                              (/1,1,1,1,bandNumber,gPointSetNumber/), &
      &                              (/1,Tdiff,plower,numGPoints,1,1/),      &
      &                              kao)

    nf_status = p_nf_inq_varid(fileid,'KeySpeciesAbsorptionCoefficientsUpperAtmos', varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                          &
      &                              (/1,1,1,1,bandNumber,gPointSetNumber/), &
      &                              (/1,Tdiff,pupper,numGPoints,1,1/),      &
      &                              kbo)

    nf_status = p_nf_inq_varid(fileid,'H20SelfAbsorptionCoefficients',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                      &
      &                              (/1,1,bandNumber,gPointSetNumber/), &
      &                              (/Tself,numGPoints,1,1/),           &
      &                              selfrefo)

    nf_status = p_nf_inq_varid(fileid,'H20ForeignAbsorptionCoefficients',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                      &
      &                              (/1,1,bandNumber,gPointSetNumber/), &
      &                              (/Tforeign,numGPoints,1,1/),        &
      &                              forrefo)

    nf_status = p_nf_inq_varid(fileid,'AbsorptionCoefficientsLowerAtmos',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                                            &
      &                              (/1,1,1,AbsorberIndex('N2'),bandNumber,gPointSetNumber/), &
      &                              (/1,T,numGPoints,1,1,1/),                                 &
      &                              kao_mn2)

    nf_status = p_nf_inq_varid(fileid,'AbsorptionCoefficientsUpperAtmos',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                                            &
      &                              (/1,1,1,AbsorberIndex('N2'),bandNumber,gPointSetNumber/), &
      &                              (/1,T,numGPoints,1,1,1/),                                 &
      &                              kbo_mn2)

  END SUBROUTINE lw_kgb01

  SUBROUTINE lw_kgb02

    USE mo_rrlw_kg02, ONLY : fracrefao, fracrefbo, kao, kbo, selfrefo, forrefo, no2

    INTEGER, PARAMETER :: bandNumber = 2
    INTEGER, PARAMETER :: numGPoints = no2


    nf_status = p_nf_inq_varid(fileid,'PlanckFractionLowerAtmos',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                      &
      &                              (/1,1,bandNumber,gPointSetNumber/), &
      &                              (/numGPoints,1,1,1/),               &
      &                              fracrefao)

    nf_status = p_nf_inq_varid(fileid,'PlanckFractionUpperAtmos',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                      &
      &                              (/1,1,bandNumber,gPointSetNumber/), &
      &                              (/numGPoints,1,1,1/),               &
      &                              fracrefbo)

    nf_status = p_nf_inq_varid(fileid,'KeySpeciesAbsorptionCoefficientsLowerAtmos',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                          &
      &                              (/1,1,1,1,bandNumber,gPointSetNumber/), &
      &                              (/1,Tdiff,plower,numGPoints,1,1/),      &
      &                              kao)

    nf_status = p_nf_inq_varid(fileid,'KeySpeciesAbsorptionCoefficientsUpperAtmos',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                          &
      &                              (/1,1,1,1,bandNumber,gPointSetNumber/), &
      &                              (/1,Tdiff,pupper,numGPoints,1,1/),      &
      &                              kbo)

    nf_status = p_nf_inq_varid(fileid,'H20SelfAbsorptionCoefficients',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                      &
      &                              (/1,1,bandNumber,gPointSetNumber/), &
      &                              (/Tself,numGPoints,1,1/),           &
      &                              selfrefo)

    nf_status = p_nf_inq_varid(fileid,'H20ForeignAbsorptionCoefficients',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                      &
      &                              (/1,1,bandNumber,gPointSetNumber/), &
      &                              (/Tforeign,numGPoints,1,1/),        &
      &                              forrefo)

  END SUBROUTINE lw_kgb02

  SUBROUTINE lw_kgb03

    USE mo_rrlw_kg03, ONLY : fracrefao, fracrefbo, kao, kbo, kao_mn2o, kbo_mn2o, selfrefo, forrefo,no3

    INTEGER, PARAMETER :: bandNumber = 3
    INTEGER, PARAMETER :: numGPoints = no3


    nf_status = p_nf_inq_varid(fileid,'PlanckFractionLowerAtmos',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                      &
      &                              (/1,1,bandNumber,gPointSetNumber/), &
      &                              (/numGPoints,keylower,1,1/),        &
      &                              fracrefao)

    nf_status = p_nf_inq_varid(fileid,'PlanckFractionUpperAtmos',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                      &
      &                              (/1,1,bandNumber,gPointSetNumber/), &
      &                              (/numGPoints,keyupper,1,1/),        &
      &                              fracrefbo)

    nf_status = p_nf_inq_varid(fileid,'KeySpeciesAbsorptionCoefficientsLowerAtmos',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                            &
      &                              (/1,1,1,1,bandNumber,gPointSetNumber/),   &
      &                              (/keylower,Tdiff,plower,numGPoints,1,1/), &
      &                              kao)

    nf_status = p_nf_inq_varid(fileid,'KeySpeciesAbsorptionCoefficientsUpperAtmos',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                            &
      &                              (/1,1,1,1,bandNumber,gPointSetNumber/),   &
      &                              (/keyupper,Tdiff,pupper,numGPoints,1,1/), &
      &                              kbo)

    nf_status = p_nf_inq_varid(fileid,'H20SelfAbsorptionCoefficients',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                      &
      &                              (/1,1,bandNumber,gPointSetNumber/), &
      &                              (/Tself,numGPoints,1,1/),           &
      &                              selfrefo)

    nf_status = p_nf_inq_varid(fileid,'H20ForeignAbsorptionCoefficients',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                      &
      &                              (/1,1,bandNumber,gPointSetNumber/), &
      &                              (/Tforeign,numGPoints,1,1/),        &
      &                              forrefo)

    nf_status = p_nf_inq_varid(fileid,'AbsorptionCoefficientsLowerAtmos',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                                             &
      &                              (/1,1,1,AbsorberIndex('N2O'),bandNumber,gPointSetNumber/), &
      &                              (/keylower,T,numGPoints,1,1,1/),                           &
      &                              kao_mn2o)

    nf_status = p_nf_inq_varid(fileid,'AbsorptionCoefficientsUpperAtmos',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                                             &
      &                              (/1,1,1,AbsorberIndex('N2O'),bandNumber,gPointSetNumber/), &
      &                              (/keyupper,T,numGPoints,1,1,1/),                           &
      &                              kbo_mn2o)

  END SUBROUTINE lw_kgb03

  SUBROUTINE lw_kgb04

    USE mo_rrlw_kg04, ONLY : fracrefao, fracrefbo, kao, kbo, selfrefo, forrefo, no4

    INTEGER, PARAMETER :: bandNumber = 4
    INTEGER, PARAMETER :: numGPoints = no4


    nf_status = p_nf_inq_varid(fileid,'PlanckFractionLowerAtmos',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                      &
      &                              (/1,1,bandNumber,gPointSetNumber/), &
      &                              (/numGPoints,keylower,1,1/),        &
      &                              fracrefao)

    nf_status = p_nf_inq_varid(fileid,'PlanckFractionUpperAtmos',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                      &
      &                              (/1,1,bandNumber,gPointSetNumber/), &
      &                              (/numGPoints,keyupper,1,1/),        &
      &                              fracrefbo(:,1:5))

    nf_status = p_nf_inq_varid(fileid,'KeySpeciesAbsorptionCoefficientsLowerAtmos',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                            &
      &                              (/1,1,1,1,bandNumber,gPointSetNumber/),   &
      &                              (/keylower,Tdiff,plower,numGPoints,1,1/), &
      &                              kao)

    nf_status = p_nf_inq_varid(fileid,'KeySpeciesAbsorptionCoefficientsUpperAtmos',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                            &
      &                              (/1,1,1,1,bandNumber,gPointSetNumber/),   &
      &                              (/keyupper,Tdiff,pupper,numGPoints,1,1/), &
      &                              kbo)

    nf_status = p_nf_inq_varid(fileid,'H20SelfAbsorptionCoefficients',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                      &
      &                              (/1,1,bandNumber,gPointSetNumber/), &
      &                              (/Tself,numGPoints,1,1/),           &
      &                              selfrefo)

    nf_status = p_nf_inq_varid(fileid,'H20ForeignAbsorptionCoefficients',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                      &
      &                              (/1,1,bandNumber,gPointSetNumber/), &
      &                              (/Tforeign,numGPoints,1,1/),        &
      &                              forrefo)

  END SUBROUTINE lw_kgb04

  SUBROUTINE lw_kgb05

    USE mo_rrlw_kg05, ONLY : fracrefao, fracrefbo, kao, kbo, kao_mo3, selfrefo, forrefo, ccl4o, no5

    INTEGER, PARAMETER :: bandNumber = 5
    INTEGER, PARAMETER :: numGPoints = no5


    nf_status = p_nf_inq_varid(fileid,'PlanckFractionLowerAtmos',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                      &
      &                              (/1,1,bandNumber,gPointSetNumber/), &
      &                              (/numGPoints,keylower,1,1/),        &
      &                              fracrefao)

    nf_status = p_nf_inq_varid(fileid,'PlanckFractionUpperAtmos',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                      &
      &                              (/1,1,bandNumber,gPointSetNumber/), &
      &                              (/numGPoints,keyupper,1,1/),        &
      &                              fracrefbo)

    nf_status = p_nf_inq_varid(fileid,'KeySpeciesAbsorptionCoefficientsLowerAtmos',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                            &
      &                              (/1,1,1,1,bandNumber,gPointSetNumber/),   &
      &                              (/keylower,Tdiff,plower,numGPoints,1,1/), &
      &                              kao)

    nf_status = p_nf_inq_varid(fileid,'KeySpeciesAbsorptionCoefficientsUpperAtmos',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                            &
      &                              (/1,1,1,1,bandNumber,gPointSetNumber/),   &
      &                              (/keyupper,Tdiff,pupper,numGPoints,1,1/), &
      &                              kbo)

    nf_status = p_nf_inq_varid(fileid,'H20SelfAbsorptionCoefficients',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                      &
      &                              (/1,1,bandNumber,gPointSetNumber/), &
      &                              (/Tself,numGPoints,1,1/),           &
      &                              selfrefo)

    nf_status = p_nf_inq_varid(fileid,'H20ForeignAbsorptionCoefficients',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                      &
      &                              (/1,1,bandNumber,gPointSetNumber/), &
      &                              (/Tforeign,numGPoints,1,1/),        &
      &                              forrefo)

    nf_status = p_nf_inq_varid(fileid,'AbsorptionCoefficientsLowerAtmos',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                                            &
      &                              (/1,1,1,AbsorberIndex('O3'),bandNumber,gPointSetNumber/), &
      &                              (/keylower,T,numGPoints,1,1,1/),                          &
      &                              kao_mo3)

    nf_status = p_nf_inq_varid(fileid,'AbsorptionCoefficientsLowerAtmos',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                                              &
      &                              (/1,1,1,AbsorberIndex('CCL4'),bandNumber,gPointSetNumber/), &
      &                              (/1,1,numGPoints,1,1,1/),                                   &
      &                              ccl4o)

  END SUBROUTINE lw_kgb05

  SUBROUTINE lw_kgb06

    USE mo_rrlw_kg06, ONLY : fracrefao, kao, kao_mco2, selfrefo, forrefo, cfc11adjo, cfc12o, no6

    INTEGER, PARAMETER :: bandNumber = 6
    INTEGER, PARAMETER :: numGPoints = no6


    nf_status = p_nf_inq_varid(fileid,'PlanckFractionLowerAtmos',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                      &
      &                              (/1,1,bandNumber,gPointSetNumber/), &
      &                              (/numGPoints,1,1,1/),               &
      &                              fracrefao)

    nf_status = p_nf_inq_varid(fileid,'KeySpeciesAbsorptionCoefficientsLowerAtmos',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                          &
      &                              (/1,1,1,1,bandNumber,gPointSetNumber/), &
      &                              (/1,Tdiff,plower,numGPoints,1,1/),      &
      &                              kao)

    nf_status = p_nf_inq_varid(fileid,'H20SelfAbsorptionCoefficients',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                      &
      &                              (/1,1,bandNumber,gPointSetNumber/), &
      &                              (/Tself,numGPoints,1,1/),           &
      &                              selfrefo)

    nf_status = p_nf_inq_varid(fileid,'H20ForeignAbsorptionCoefficients',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                      &
      &                              (/1,1,bandNumber,gPointSetNumber/), &
      &                              (/Tforeign,numGPoints,1,1/),        &
      &                              forrefo)

    nf_status = p_nf_inq_varid(fileid,'AbsorptionCoefficientsLowerAtmos',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                                             &
      &                              (/1,1,1,AbsorberIndex('CO2'),bandNumber,gPointSetNumber/), &
      &                              (/1,T,numGPoints,1,1,1/),                                  &
      &                              kao_mco2)

    nf_status = p_nf_inq_varid(fileid,'AbsorptionCoefficientsLowerAtmos',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                                               &
      &                              (/1,1,1,AbsorberIndex('CFC11'),bandNumber,gPointSetNumber/), &
      &                              (/1,1,numGPoints,1,1,1/),                                    &
      &                              cfc11adjo)

    nf_status = p_nf_inq_varid(fileid,'AbsorptionCoefficientsLowerAtmos',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                                               &
      &                              (/1,1,1,AbsorberIndex('CFC12'),bandNumber,gPointSetNumber/), &
      &                              (/1,1,numGPoints,1,1,1/),                                    &
      &                              cfc12o)

  END SUBROUTINE lw_kgb06

  SUBROUTINE lw_kgb07

    USE mo_rrlw_kg07, ONLY : fracrefao, fracrefbo, kao, kbo, kao_mco2, kbo_mco2, selfrefo, forrefo,no7

    INTEGER, PARAMETER :: bandNumber = 7
    INTEGER, PARAMETER :: numGPoints = no7


    nf_status = p_nf_inq_varid(fileid,'PlanckFractionLowerAtmos',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                      &
      &                              (/1,1,bandNumber,gPointSetNumber/), &
      &                              (/numGPoints,keylower,1,1/),        &
      &                              fracrefao)

    nf_status = p_nf_inq_varid(fileid,'PlanckFractionUpperAtmos',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                      &
      &                              (/1,1,bandNumber,gPointSetNumber/), &
      &                              (/numGPoints,1,1,1/),               &
      &                              fracrefbo)

    nf_status = p_nf_inq_varid(fileid,'KeySpeciesAbsorptionCoefficientsLowerAtmos',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                            &
      &                              (/1,1,1,1,bandNumber,gPointSetNumber/),   &
      &                              (/keylower,Tdiff,plower,numGPoints,1,1/), &
      &                              kao)

    nf_status = p_nf_inq_varid(fileid,'KeySpeciesAbsorptionCoefficientsUpperAtmos',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                          &
      &                              (/1,1,1,1,bandNumber,gPointSetNumber/), &
      &                              (/1,Tdiff,pupper,numGPoints,1,1/),      &
      &                              kbo)

    nf_status = p_nf_inq_varid(fileid,'H20SelfAbsorptionCoefficients',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                      &
      &                              (/1,1,bandNumber,gPointSetNumber/), &
      &                              (/Tself,numGPoints,1,1/),           &
      &                              selfrefo)

    nf_status = p_nf_inq_varid(fileid,'H20ForeignAbsorptionCoefficients',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                      &
      &                              (/1,1,bandNumber,gPointSetNumber/), &
      &                              (/Tforeign,numGPoints,1,1/),        &
      &                              forrefo)

    nf_status = p_nf_inq_varid(fileid,'AbsorptionCoefficientsLowerAtmos',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                                             &
      &                              (/1,1,1,AbsorberIndex('CO2'),bandNumber,gPointSetNumber/), &
      &                              (/keylower,T,numGPoints,1,1,1/),                           &
      &                              kao_mco2)

    nf_status = p_nf_inq_varid(fileid,'AbsorptionCoefficientsUpperAtmos',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                                             &
      &                              (/1,1,1,AbsorberIndex('CO2'),bandNumber,gPointSetNumber/), &
      &                              (/1,T,numGPoints,1,1,1/),                                  &
      &                              kbo_mco2)

  END SUBROUTINE lw_kgb07

  SUBROUTINE lw_kgb08

    USE mo_rrlw_kg08, ONLY : fracrefao, fracrefbo, kao, kao_mco2, kao_mn2o, kao_mo3, &
      &                      kbo, kbo_mco2, kbo_mn2o, selfrefo, forrefo, cfc12o, cfc22adjo, no8

    INTEGER, PARAMETER :: bandNumber = 8
    INTEGER, PARAMETER :: numGPoints = no8


    nf_status = p_nf_inq_varid(fileid,'PlanckFractionLowerAtmos',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                      &
      &                              (/1,1,bandNumber,gPointSetNumber/), &
      &                              (/numGPoints,1,1,1/),               &
      &                              fracrefao)

    nf_status = p_nf_inq_varid(fileid,'PlanckFractionUpperAtmos',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                      &
      &                              (/1,1,bandNumber,gPointSetNumber/), &
      &                              (/numGPoints,1,1,1/),               &
      &                              fracrefbo)

    nf_status = p_nf_inq_varid(fileid,'KeySpeciesAbsorptionCoefficientsLowerAtmos',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                          &
      &                              (/1,1,1,1,bandNumber,gPointSetNumber/), &
      &                              (/1,Tdiff,plower,numGPoints,1,1/),      &
      &                              kao)

    nf_status = p_nf_inq_varid(fileid,'KeySpeciesAbsorptionCoefficientsUpperAtmos',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                          &
      &                              (/1,1,1,1,bandNumber,gPointSetNumber/), &
      &                              (/1,Tdiff,pupper,numGPoints,1,1/),      &
      &                              kbo)

    nf_status = p_nf_inq_varid(fileid,'H20SelfAbsorptionCoefficients',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                      &
      &                              (/1,1,bandNumber,gPointSetNumber/), &
      &                              (/Tself,numGPoints,1,1/),           &
      &                              selfrefo)

    nf_status = p_nf_inq_varid(fileid,'H20ForeignAbsorptionCoefficients',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                      &
      &                              (/1,1,bandNumber,gPointSetNumber/), &
      &                              (/Tforeign,numGPoints,1,1/),        &
      &                              forrefo)

    nf_status = p_nf_inq_varid(fileid,'AbsorptionCoefficientsLowerAtmos',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                                            &
      &                              (/1,1,1,AbsorberIndex('O3'),bandNumber,gPointSetNumber/), &
      &                              (/1,T,numGPoints,1,1,1/),                                 &
      &                              kao_mo3)

    nf_status = p_nf_inq_varid(fileid,'AbsorptionCoefficientsLowerAtmos',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                                             &
      &                              (/1,1,1,AbsorberIndex('CO2'),bandNumber,gPointSetNumber/), &
      &                              (/1,T,numGPoints,1,1,1/),                                  &
      &                              kao_mco2)

    nf_status = p_nf_inq_varid(fileid,'AbsorptionCoefficientsUpperAtmos',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                                             &
      &                              (/1,1,1,AbsorberIndex('CO2'),bandNumber,gPointSetNumber/), &
      &                              (/1,T,numGPoints,1,1,1/),                                  &
      &                              kbo_mco2)

    nf_status = p_nf_inq_varid(fileid,'AbsorptionCoefficientsLowerAtmos',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                                             &
      &                              (/1,1,1,AbsorberIndex('N2O'),bandNumber,gPointSetNumber/), &
      &                              (/1,T,numGPoints,1,1,1/),                                  &
      &                              kao_mn2o)

    nf_status = p_nf_inq_varid(fileid,'AbsorptionCoefficientsUpperAtmos',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                                             &
      &                              (/1,1,1,AbsorberIndex('N2O'),bandNumber,gPointSetNumber/), &
      &                              (/1,T,numGPoints,1,1,1/),                                  &
      &                              kbo_mn2o)

    nf_status = p_nf_inq_varid(fileid,'AbsorptionCoefficientsLowerAtmos',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                                               &
      &                              (/1,1,1,AbsorberIndex('CFC12'),bandNumber,gPointSetNumber/), &
      &                              (/1,1,numGPoints,1,1,1/),                                    &
      &                              cfc12o)

    nf_status = p_nf_inq_varid(fileid,'AbsorptionCoefficientsLowerAtmos',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                                               &
      &                              (/1,1,1,AbsorberIndex('CFC22'),bandNumber,gPointSetNumber/), &
      &                              (/1,1,numGPoints,1,1,1/),                                    &
      &                              cfc22adjo)

  END SUBROUTINE lw_kgb08

  SUBROUTINE lw_kgb09

    USE mo_rrlw_kg09, ONLY : fracrefao, fracrefbo, kao, kbo, kao_mn2o, kbo_mn2o, selfrefo, forrefo,no9

    INTEGER, PARAMETER :: bandNumber = 9
    INTEGER, PARAMETER :: numGPoints = no9


    nf_status = p_nf_inq_varid(fileid,'PlanckFractionLowerAtmos',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                      &
      &                              (/1,1,bandNumber,gPointSetNumber/), &
      &                              (/numGPoints,keylower,1,1/),        &
      &                              fracrefao)

    nf_status = p_nf_inq_varid(fileid,'PlanckFractionUpperAtmos',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                      &
      &                              (/1,1,bandNumber,gPointSetNumber/), &
      &                              (/numGPoints,1,1,1/),               &
      &                              fracrefbo)

    nf_status = p_nf_inq_varid(fileid,'KeySpeciesAbsorptionCoefficientsLowerAtmos',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                            &
      &                              (/1,1,1,1,bandNumber,gPointSetNumber/),   &
      &                              (/keylower,Tdiff,plower,numGPoints,1,1/), &
      &                              kao)

    nf_status = p_nf_inq_varid(fileid,'KeySpeciesAbsorptionCoefficientsUpperAtmos',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                          &
      &                              (/1,1,1,1,bandNumber,gPointSetNumber/), &
      &                              (/1,Tdiff,pupper,numGPoints,1,1/),      &
      &                              kbo)

    nf_status = p_nf_inq_varid(fileid,'H20SelfAbsorptionCoefficients',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                      &
      &                              (/1,1,bandNumber,gPointSetNumber/), &
      &                              (/Tself,numGPoints,1,1/),           &
      &                              selfrefo)

    nf_status = p_nf_inq_varid(fileid,'H20ForeignAbsorptionCoefficients',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                      &
      &                              (/1,1,bandNumber,gPointSetNumber/), &
      &                              (/Tforeign,numGPoints,1,1/),        &
      &                              forrefo)

    nf_status = p_nf_inq_varid(fileid,'AbsorptionCoefficientsLowerAtmos',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                                             &
      &                              (/1,1,1,AbsorberIndex('N2O'),bandNumber,gPointSetNumber/), &
      &                              (/keylower,T,numGPoints,1,1,1/),                           &
      &                              kao_mn2o)

    nf_status = p_nf_inq_varid(fileid,'AbsorptionCoefficientsUpperAtmos',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                                             &
      &                              (/1,1,1,AbsorberIndex('N2O'),bandNumber,gPointSetNumber/), &
      &                              (/1,T,numGPoints,1,1,1/),                                  &
      &                              kbo_mn2o)

  END SUBROUTINE lw_kgb09

  SUBROUTINE lw_kgb10

    USE mo_rrlw_kg10, ONLY : fracrefao, fracrefbo, kao, kbo, selfrefo, forrefo, no10

    INTEGER, PARAMETER :: bandNumber = 10
    INTEGER, PARAMETER :: numGPoints = no10


    nf_status = p_nf_inq_varid(fileid,'PlanckFractionLowerAtmos',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                      &
      &                              (/1,1,bandNumber,gPointSetNumber/), &
      &                              (/numGPoints,1,1,1/),               &
      &                              fracrefao)

    nf_status = p_nf_inq_varid(fileid,'PlanckFractionUpperAtmos',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                      &
      &                              (/1,1,bandNumber,gPointSetNumber/), &
      &                              (/numGPoints,1,1,1/),               &
      &                              fracrefbo)

    nf_status = p_nf_inq_varid(fileid,'KeySpeciesAbsorptionCoefficientsLowerAtmos',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                          &
      &                              (/1,1,1,1,bandNumber,gPointSetNumber/), &
      &                              (/1,Tdiff,plower,numGPoints,1,1/),      &
      &                              kao)

    nf_status = p_nf_inq_varid(fileid,'KeySpeciesAbsorptionCoefficientsUpperAtmos',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                          &
      &                              (/1,1,1,1,bandNumber,gPointSetNumber/), &
      &                              (/1,Tdiff,pupper,numGPoints,1,1/),      &
      &                              kbo)

    nf_status = p_nf_inq_varid(fileid,'H20SelfAbsorptionCoefficients',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                      &
      &                              (/1,1,bandNumber,gPointSetNumber/), &
      &                              (/Tself,numGPoints,1,1/),           &
      &                              selfrefo)

    nf_status = p_nf_inq_varid(fileid,'H20ForeignAbsorptionCoefficients',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                      &
      &                              (/1,1,bandNumber,gPointSetNumber/), &
      &                              (/Tforeign,numGPoints,1,1/),        &
      &                              forrefo)

  END SUBROUTINE lw_kgb10

  SUBROUTINE lw_kgb11
    USE mo_rrlw_kg11, ONLY : fracrefao, fracrefbo, kao, kbo, kao_mo2, kbo_mo2, selfrefo, forrefo, no11

    INTEGER, PARAMETER :: bandNumber = 11
    INTEGER, PARAMETER :: numGPoints = no11


    nf_status = p_nf_inq_varid(fileid,'PlanckFractionLowerAtmos',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                      &
      &                              (/1,1,bandNumber,gPointSetNumber/), &
      &                              (/numGPoints,1,1,1/),               &
      &                              fracrefao)

    nf_status = p_nf_inq_varid(fileid,'PlanckFractionUpperAtmos',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                      &
      &                              (/1,1,bandNumber,gPointSetNumber/), &
      &                              (/numGPoints,1,1,1/),               &
      &                              fracrefbo)

    nf_status = p_nf_inq_varid(fileid,'KeySpeciesAbsorptionCoefficientsLowerAtmos',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                          &
      &                              (/1,1,1,1,bandNumber,gPointSetNumber/), &
      &                              (/1,Tdiff,plower,numGPoints,1,1/),      &
      &                              kao)

    nf_status = p_nf_inq_varid(fileid,'KeySpeciesAbsorptionCoefficientsUpperAtmos',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                          &
      &                              (/1,1,1,1,bandNumber,gPointSetNumber/), &
      &                              (/1,Tdiff,pupper,numGPoints,1,1/),      &
      &                              kbo)

    nf_status = p_nf_inq_varid(fileid,'H20SelfAbsorptionCoefficients',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                      &
      &                              (/1,1,bandNumber,gPointSetNumber/), &
      &                              (/Tself,numGPoints,1,1/),           &
      &                              selfrefo)

    nf_status = p_nf_inq_varid(fileid,'H20ForeignAbsorptionCoefficients',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                      &
      &                              (/1,1,bandNumber,gPointSetNumber/), &
      &                              (/Tforeign,numGPoints,1,1/),        &
      &                              forrefo)


    nf_status = p_nf_inq_varid(fileid,'AbsorptionCoefficientsLowerAtmos',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                                            &
      &                              (/1,1,1,AbsorberIndex('O2'),bandNumber,gPointSetNumber/), &
      &                              (/1,T,numGPoints,1,1,1/),                                 &
      &                              kao_mo2)

    nf_status = p_nf_inq_varid(fileid,'AbsorptionCoefficientsUpperAtmos',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                                            &
      &                              (/1,1,1,AbsorberIndex('O2'),bandNumber,gPointSetNumber/), &
      &                              (/1,T,numGPoints,1,1,1/),                                 &
      &                              kbo_mo2)

  END SUBROUTINE lw_kgb11

  SUBROUTINE lw_kgb12

    USE mo_rrlw_kg12, ONLY : fracrefao, kao, selfrefo, forrefo, no12

    INTEGER, PARAMETER :: bandNumber = 12
    INTEGER, PARAMETER :: numGPoints = no12


    nf_status = p_nf_inq_varid(fileid,'PlanckFractionLowerAtmos',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                      &
      &                              (/1,1,bandNumber,gPointSetNumber/), &
      &                              (/numGPoints,keylower,1,1/),        &
      &                              fracrefao)

    nf_status = p_nf_inq_varid(fileid,'KeySpeciesAbsorptionCoefficientsLowerAtmos',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                            &
      &                              (/1,1,1,1,bandNumber,gPointSetNumber/),   &
      &                              (/keylower,Tdiff,plower,numGPoints,1,1/), &
      &                              kao)

    nf_status = p_nf_inq_varid(fileid,'H20SelfAbsorptionCoefficients',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                      &
      &                              (/1,1,bandNumber,gPointSetNumber/), &
      &                              (/Tself,numGPoints,1,1/),           &
      &                              selfrefo)

    nf_status = p_nf_inq_varid(fileid,'H20ForeignAbsorptionCoefficients',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                      &
      &                              (/1,1,bandNumber,gPointSetNumber/), &
      &                              (/Tforeign,numGPoints,1,1/),        &
      &                              forrefo)

  END SUBROUTINE lw_kgb12

  SUBROUTINE lw_kgb13

    USE mo_rrlw_kg13, ONLY : fracrefao, fracrefbo, kao, kao_mco2, kao_mco, kbo_mo3, &
      &                      selfrefo, forrefo, no13

    INTEGER, PARAMETER :: bandNumber = 13
    INTEGER, PARAMETER :: numGPoints = no13


    nf_status = p_nf_inq_varid(fileid,'PlanckFractionLowerAtmos',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                      &
      &                              (/1,1,bandNumber,gPointSetNumber/), &
      &                              (/numGPoints,keylower,1,1/),        &
      &                              fracrefao)

    nf_status = p_nf_inq_varid(fileid,'PlanckFractionUpperAtmos',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                       &
      &                              (/1,1,bandNumber,gPointSetNumber/),  &
      &                              (/numGPoints,1,1,1/),                &
      &                              fracrefbo)

    nf_status = p_nf_inq_varid(fileid,'KeySpeciesAbsorptionCoefficientsLowerAtmos',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                            &
      &                              (/1,1,1,1,bandNumber,gPointSetNumber/),   &
      &                              (/keylower,Tdiff,plower,numGPoints,1,1/), &
      &                              kao)

    nf_status = p_nf_inq_varid(fileid,'H20SelfAbsorptionCoefficients',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                      &
      &                              (/1,1,bandNumber,gPointSetNumber/), &
      &                              (/Tself,numGPoints,1,1/),           &
      &                              selfrefo)

    nf_status = p_nf_inq_varid(fileid,'H20ForeignAbsorptionCoefficients',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                      &
      &                              (/1,1,bandNumber,gPointSetNumber/), &
      &                              (/Tforeign,numGPoints,1,1/),        &
      &                              forrefo)

    nf_status = p_nf_inq_varid(fileid,'AbsorptionCoefficientsUpperAtmos',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                                            &
      &                              (/1,1,1,AbsorberIndex('O3'),bandNumber,gPointSetNumber/), &
      &                              (/1,T,numGPoints,1,1,1/),                                 &
      &                              kbo_mo3)

    nf_status = p_nf_inq_varid(fileid,'AbsorptionCoefficientsLowerAtmos',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                                             &
      &                              (/1,1,1,AbsorberIndex('CO2'),bandNumber,gPointSetNumber/), &
      &                              (/keylower,T,numGPoints,1,1,1/),                           &
      &                              kao_mco2)

    nf_status = p_nf_inq_varid(fileid,'AbsorptionCoefficientsLowerAtmos',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                                            &
      &                              (/1,1,1,AbsorberIndex('CO'),bandNumber,gPointSetNumber/), &
      &                              (/keylower,T,numGPoints,1,1,1/),                          &
      &                              kao_mco)

  END SUBROUTINE lw_kgb13

  SUBROUTINE lw_kgb14

    USE mo_rrlw_kg14, ONLY : fracrefao, fracrefbo, kao, kbo, selfrefo, forrefo, no14

    INTEGER, PARAMETER :: bandNumber = 14
    INTEGER, PARAMETER :: numGPoints = no14


    nf_status = p_nf_inq_varid(fileid,'PlanckFractionLowerAtmos',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                      &
      &                              (/1,1,bandNumber,gPointSetNumber/), &
      &                              (/numGPoints,1,1,1/),               &
      &                              fracrefao)

    nf_status = p_nf_inq_varid(fileid,'PlanckFractionUpperAtmos',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                      &
      &                              (/1,1,bandNumber,gPointSetNumber/), &
      &                              (/numGPoints,1,1,1/),               &
      &                              fracrefbo)

    nf_status = p_nf_inq_varid(fileid,'KeySpeciesAbsorptionCoefficientsLowerAtmos',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                          &
      &                              (/1,1,1,1,bandNumber,gPointSetNumber/), &
      &                              (/1,Tdiff,plower,numGPoints,1,1/),      &
      &                              kao)

    nf_status = p_nf_inq_varid(fileid,'KeySpeciesAbsorptionCoefficientsUpperAtmos',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                          &
      &                              (/1,1,1,1,bandNumber,gPointSetNumber/), &
      &                              (/1,Tdiff,pupper,numGPoints,1,1/),      &
      &                              kbo)

    nf_status = p_nf_inq_varid(fileid,'H20SelfAbsorptionCoefficients',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                      &
      &                              (/1,1,bandNumber,gPointSetNumber/), &
      &                              (/Tself,numGPoints,1,1/),           &
      &                              selfrefo)

    nf_status = p_nf_inq_varid(fileid,'H20ForeignAbsorptionCoefficients',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                      &
      &                              (/1,1,bandNumber,gPointSetNumber/), &
      &                              (/Tforeign,numGPoints,1,1/),        &
      &                              forrefo)

  END SUBROUTINE lw_kgb14

  SUBROUTINE lw_kgb15

    USE mo_rrlw_kg15, ONLY : fracrefao, kao, kao_mn2, selfrefo, forrefo, no15

    INTEGER, PARAMETER :: bandNumber = 15
    INTEGER, PARAMETER :: numGPoints = no15


    nf_status = p_nf_inq_varid(fileid,'PlanckFractionLowerAtmos',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                      &
      &                              (/1,1,bandNumber,gPointSetNumber/), &
      &                              (/numGPoints,keylower,1,1/),        &
      &                              fracrefao)

    nf_status = p_nf_inq_varid(fileid,'KeySpeciesAbsorptionCoefficientsLowerAtmos',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                            &
      &                              (/1,1,1,1,bandNumber,gPointSetNumber/),   &
      &                              (/keylower,Tdiff,plower,numGPoints,1,1/), &
      &                              kao)

    nf_status = p_nf_inq_varid(fileid,'H20SelfAbsorptionCoefficients',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                      &
      &                              (/1,1,bandNumber,gPointSetNumber/), &
      &                              (/Tself,numGPoints,1,1/),           &
      &                              selfrefo)

    nf_status = p_nf_inq_varid(fileid,'H20ForeignAbsorptionCoefficients',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                      &
      &                              (/1,1,bandNumber,gPointSetNumber/), &
      &                              (/Tforeign,numGPoints,1,1/),        &
      &                              forrefo)

    nf_status = p_nf_inq_varid(fileid,'AbsorptionCoefficientsLowerAtmos',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                                            &
      &                              (/1,1,1,AbsorberIndex('N2'),bandNumber,gPointSetNumber/), &
      &                              (/keylower,T,numGPoints,1,1,1/),                          &
      &                              kao_mn2)

  END SUBROUTINE lw_kgb15

  SUBROUTINE lw_kgb16

    USE mo_rrlw_kg16, ONLY : fracrefao, fracrefbo, kao, kbo, selfrefo, forrefo, no16

    INTEGER, PARAMETER :: bandNumber = 16
    INTEGER, PARAMETER :: numGPoints = no16


    nf_status = p_nf_inq_varid(fileid,'PlanckFractionLowerAtmos',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                      &
      &                              (/1,1,bandNumber,gPointSetNumber/), &
      &                              (/numGPoints,keylower,1,1/),        &
      &                              fracrefao)

    nf_status = p_nf_inq_varid(fileid,'PlanckFractionUpperAtmos',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                      &
      &                              (/1,1,bandNumber,gPointSetNumber/), &
      &                              (/numGPoints,1,1,1/),               &
      &                              fracrefbo)

    nf_status = p_nf_inq_varid(fileid,'KeySpeciesAbsorptionCoefficientsLowerAtmos',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                            &
      &                              (/1,1,1,1,bandNumber,gPointSetNumber/),   &
      &                              (/keylower,Tdiff,plower,numGPoints,1,1/), &
      &                              kao)

    nf_status = p_nf_inq_varid(fileid,'KeySpeciesAbsorptionCoefficientsUpperAtmos',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                          &
      &                              (/1,1,1,1,bandNumber,gPointSetNumber/), &
      &                              (/1,Tdiff,pupper,numGPoints,1,1/),      &
      &                              kbo)

    nf_status = p_nf_inq_varid(fileid,'H20SelfAbsorptionCoefficients',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                      &
      &                              (/1,1,bandNumber,gPointSetNumber/), &
      &                              (/Tself,numGPoints,1,1/),           &
      &                              selfrefo)

    nf_status = p_nf_inq_varid(fileid,'H20ForeignAbsorptionCoefficients',varid)
    nf_status = p_nf_get_vara_double(fileid, varid,                      &
      &                              (/1,1,bandNumber,gPointSetNumber/), &
      &                              (/Tforeign,numGPoints,1,1/),        &
      &                              forrefo)

  END SUBROUTINE lw_kgb16

END MODULE mo_lrtm_netcdf
