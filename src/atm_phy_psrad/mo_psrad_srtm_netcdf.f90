!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_psrad_srtm_netcdf

  USE mo_psrad_general, ONLY : wp, finish
  USE mo_psrad_io, ONLY: psrad_io_open, psrad_io_close, &
    read=>psrad_io_copy_double
  USE mo_psrad_srtm_kgs, ONLY : no

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: srtm_read

  INTEGER, PARAMETER :: keylower = 9, keyupper = 5, Tdiff = 5, ps = 59, &
    plower = 13, pupper = 47, Tself = 10, T = 19, band = 14, pforeign = 4, &
    Tforeignlower = 3, Tforeignupper = 2, GPoint = 16, GPointSet = 2 
  INTEGER, PARAMETER :: N2_idx = 1, CCL4_idx = 2, CFC11_idx = 3, &
    CFC12_idx = 4, CFC22_idx = 5, H2O_idx = 6, CO2_idx = 7, O3_idx = 8, &
    N2O_idx = 9, CO_idx = 10, CH4_idx = 11, O2_idx = 12

  INTEGER, PARAMETER :: gPointSetNumber = 1
  INTEGER :: fileid     !< id number of netcdf file

CONTAINS 

  SUBROUTINE srtm_read

    CALL psrad_io_open('rrtmg_sw.nc', fileid)
    IF (fileid == 0) THEN
      CALL finish('mo_psrad_srtm_netcdf/srtm_read', &
        'File rrtmg_sw.nc cannot be opened')
    END IF

    CALL sw_kgb16  ! molecular absorption coefficients
    CALL sw_kgb17
    CALL sw_kgb18
    CALL sw_kgb19
    CALL sw_kgb20
    CALL sw_kgb21
    CALL sw_kgb22
    CALL sw_kgb23
    CALL sw_kgb24
    CALL sw_kgb25
    CALL sw_kgb26
    CALL sw_kgb27
    CALL sw_kgb28
    CALL sw_kgb29

    CALL psrad_io_close(fileid)

  END SUBROUTINE srtm_read

  SUBROUTINE sw_kgb16

    USE mo_psrad_srtm_kgs,  ONLY : sfluxrefo=>sfluxrefo16, kao=>kao16, &
      kbo=>kbo16, selfrefo=>selfrefo16, forrefo=>forrefo16, rayl=>rayl16

    INTEGER, PARAMETER :: bandNumber = 1, numGPoints = no(1)

    REAL(kind=wp) :: ncrayl(1)

    CALL read(fileid,'SolarSourceFunctionLowerAtmos', &
      (/1,1,bandNumber,gPointSetNumber/), &
      (/numGPoints,1,1,1/), sfluxrefo)

    CALL read(fileid,'RayleighExtinctionCoefficientsLowerAtmos', &
      (/1,1,bandNumber,gPointSetNumber/), &
      (/1,1,1,1/), ncrayl)

    CALL read(fileid,'KeySpeciesAbsorptionCoefficientsLowerAtmos', &
      (/1,1,1,1,bandNumber,gPointSetNumber/), &
      (/keylower,Tdiff,plower,numGPoints,1,1/), kao)

    CALL read(fileid,'KeySpeciesAbsorptionCoefficientsUpperAtmos', &
      (/1,1,1,1,bandNumber,gPointSetNumber/), &
      (/1,Tdiff,pupper,numGPoints,1,1/), kbo)

    CALL read(fileid,'H2OSelfAbsorptionCoefficients', &
      (/1,1,bandNumber,gPointSetNumber/), &
      (/Tself,numGPoints,1,1/), selfrefo)

    CALL read(fileid,'H2OForeignAbsorptionCoefficientsLowerAtmos', &
      (/1,1,bandNumber,gPointSetNumber/), &
      (/Tforeignlower,numGPoints,1,1/), forrefo)

    rayl = ncrayl(1)
    
  END SUBROUTINE sw_kgb16

  SUBROUTINE sw_kgb17

    USE mo_psrad_srtm_kgs,  ONLY : sfluxrefo=>sfluxrefo17, kao=>kao17, &
      kbo=>kbo17, selfrefo=>selfrefo17, forrefo=>forrefo17, rayl=>rayl17

    INTEGER, PARAMETER :: bandNumber = 2
    INTEGER, PARAMETER :: numGPoints = no(2)

    REAL(kind=wp) :: ncrayl(1)

    CALL read(fileid,'SolarSourceFunctionUpperAtmos', &
      (/1,1,bandNumber,gPointSetNumber/), &
      (/numGPoints,keyupper,1,1/), sfluxrefo)

    CALL read(fileid,'RayleighExtinctionCoefficientsLowerAtmos', &
      (/1,1,bandNumber,gPointSetNumber/), &
      (/1,1,1,1/), ncrayl)

    CALL read(fileid,'KeySpeciesAbsorptionCoefficientsLowerAtmos', &
      (/1,1,1,1,bandNumber,gPointSetNumber/), &
      (/keylower,Tdiff,plower,numGPoints,1,1/), kao)

    CALL read(fileid,'KeySpeciesAbsorptionCoefficientsUpperAtmos', &
      (/1,1,1,1,bandNumber,gPointSetNumber/), &
      (/keyupper,Tdiff,pupper,numGPoints,1,1/), kbo)

    CALL read(fileid,'H2OSelfAbsorptionCoefficients', &
      (/1,1,bandNumber,gPointSetNumber/), &
      (/Tself,numGPoints,1,1/), selfrefo)

    CALL read(fileid,'H2OForeignAbsorptionCoefficientsLowerAtmos', &
      (/1,1,bandNumber,gPointSetNumber/), &
      (/Tforeignlower,numGPoints,1,1/), forrefo(1:3,:))

    CALL read(fileid,'H2OForeignAbsorptionCoefficientsUpperAtmos', &
      (/2,1,bandNumber,gPointSetNumber/), &
      (/1,numGPoints,1,1/), forrefo(4,:))

    rayl = ncrayl(1)
    
  END SUBROUTINE sw_kgb17

  SUBROUTINE sw_kgb18

    USE mo_psrad_srtm_kgs,  ONLY : sfluxrefo=>sfluxrefo18, kao=>kao18, &
      kbo=>kbo18, selfrefo=>selfrefo18, forrefo=>forrefo18, rayl=>rayl18

    INTEGER, PARAMETER :: bandNumber = 3, numGPoints = no(3)

    REAL(kind=wp) :: ncrayl(1)

    CALL read(fileid,'SolarSourceFunctionLowerAtmos', &
      (/1,1,bandNumber,gPointSetNumber/), &
      (/numGPoints,keylower,1,1/), sfluxrefo)

    CALL read(fileid,'RayleighExtinctionCoefficientsLowerAtmos', &
      (/1,1,bandNumber,gPointSetNumber/), &
      (/1,1,1,1/), ncrayl)

    CALL read(fileid,'KeySpeciesAbsorptionCoefficientsLowerAtmos', &
      (/1,1,1,1,bandNumber,gPointSetNumber/), &
      (/keylower,Tdiff,plower,numGPoints,1,1/), kao)

    CALL read(fileid,'KeySpeciesAbsorptionCoefficientsUpperAtmos', &
      (/1,1,1,1,bandNumber,gPointSetNumber/), &
      (/1,Tdiff,pupper,numGPoints,1,1/), kbo)

    CALL read(fileid,'H2OSelfAbsorptionCoefficients', &
      (/1,1,bandNumber,gPointSetNumber/), &
      (/Tself,numGPoints,1,1/), selfrefo)

    CALL read(fileid,'H2OForeignAbsorptionCoefficientsLowerAtmos', &
      (/1,1,bandNumber,gPointSetNumber/), &
      (/Tforeignlower,numGPoints,1,1/), forrefo)

    rayl = ncrayl(1)
    
  END SUBROUTINE sw_kgb18

  SUBROUTINE sw_kgb19

    USE mo_psrad_srtm_kgs,  ONLY : sfluxrefo=>sfluxrefo19, kao=>kao19, &
      kbo=>kbo19, selfrefo=>selfrefo19, forrefo=>forrefo19, rayl=>rayl19

    INTEGER, PARAMETER :: bandNumber = 4, numGPoints = no(4)

    REAL(kind=wp) :: ncrayl(1)

    CALL read(fileid,'SolarSourceFunctionLowerAtmos', &
      (/1,1,bandNumber,gPointSetNumber/), &
      (/numGPoints,keylower,1,1/), sfluxrefo)

    CALL read(fileid,'RayleighExtinctionCoefficientsLowerAtmos', &
      (/1,1,bandNumber,gPointSetNumber/), &
      (/1,1,1,1/), ncrayl)

    CALL read(fileid,'KeySpeciesAbsorptionCoefficientsLowerAtmos', &
      (/1,1,1,1,bandNumber,gPointSetNumber/), &
      (/keylower,Tdiff,plower,numGPoints,1,1/), kao)

    CALL read(fileid,'KeySpeciesAbsorptionCoefficientsUpperAtmos', &
      (/1,1,1,1,bandNumber,gPointSetNumber/), &
      (/1,Tdiff,pupper,numGPoints,1,1/), kbo)

    CALL read(fileid,'H2OSelfAbsorptionCoefficients', &
      (/1,1,bandNumber,gPointSetNumber/), &
      (/Tself,numGPoints,1,1/), selfrefo)

    CALL read(fileid,'H2OForeignAbsorptionCoefficientsLowerAtmos', &
      (/1,1,bandNumber,gPointSetNumber/), &
      (/Tforeignlower,numGPoints,1,1/), forrefo)

    rayl = ncrayl(1)

  END SUBROUTINE sw_kgb19

  SUBROUTINE sw_kgb20

    USE mo_psrad_srtm_kgs,  ONLY : sfluxrefo=>sfluxrefo20, kao=>kao20, &
      kbo=>kbo20, selfrefo=>selfrefo20, forrefo=>forrefo20, rayl=>rayl20, &
      absch4o=>absch4o20

    INTEGER, PARAMETER :: bandNumber = 5, numGPoints = no(5)

    REAL(kind=wp) :: ncrayl(1)

    CALL read(fileid,'SolarSourceFunctionLowerAtmos', &
      (/1,1,bandNumber,gPointSetNumber/), &
      (/numGPoints,1,1,1/), sfluxrefo)

    CALL read(fileid,'RayleighExtinctionCoefficientsLowerAtmos', &
      (/1,1,bandNumber,gPointSetNumber/), &
      (/1,1,1,1/), ncrayl)

    CALL read(fileid,'KeySpeciesAbsorptionCoefficientsLowerAtmos', &
      (/1,1,1,1,bandNumber,gPointSetNumber/), &
      (/1,Tdiff,plower,numGPoints,1,1/), kao)

    CALL read(fileid,'KeySpeciesAbsorptionCoefficientsUpperAtmos', &
      (/1,1,1,1,bandNumber,gPointSetNumber/), &
      (/1,Tdiff,pupper,numGPoints,1,1/), kbo)

    CALL read(fileid,'H2OSelfAbsorptionCoefficients', &
      (/1,1,bandNumber,gPointSetNumber/), &
      (/Tself,numGPoints,1,1/), selfrefo)

    CALL read(fileid,'H2OForeignAbsorptionCoefficientsLowerAtmos', &
      (/1,1,bandNumber,gPointSetNumber/), &
      (/Tforeignlower,numGPoints,1,1/), forrefo(1:3,:))

    CALL read(fileid,'H2OForeignAbsorptionCoefficientsUpperAtmos', &
      (/2,1,bandNumber,gPointSetNumber/), &
      (/1,numGPoints,1,1/), forrefo(4,:))

    CALL read(fileid,'AbsorptionCoefficientsLowerAtmos', &
      (/1,1,1,CH4_idx,bandNumber,gPointSetNumber/), &
      (/1,1,numGPoints,1,1,1/), absch4o)

    rayl = ncrayl(1)
    
  END SUBROUTINE sw_kgb20

  SUBROUTINE sw_kgb21

    USE mo_psrad_srtm_kgs,  ONLY : sfluxrefo=>sfluxrefo21, kao=>kao21, &
      kbo=>kbo21, selfrefo=>selfrefo21, forrefo=>forrefo21, rayl=>rayl21

    INTEGER, PARAMETER :: bandNumber = 6, numGPoints = no(6)

    REAL(kind=wp) :: ncrayl(1)

    CALL read(fileid,'SolarSourceFunctionLowerAtmos', &
      (/1,1,bandNumber,gPointSetNumber/), &
      (/numGPoints,keylower,1,1/), sfluxrefo)

    CALL read(fileid,'RayleighExtinctionCoefficientsLowerAtmos', &
      (/1,1,bandNumber,gPointSetNumber/), &
      (/1,1,1,1/), ncrayl)

    CALL read(fileid,'KeySpeciesAbsorptionCoefficientsLowerAtmos', &
      (/1,1,1,1,bandNumber,gPointSetNumber/), &
      (/keylower,Tdiff,plower,numGPoints,1,1/), kao)

    CALL read(fileid,'KeySpeciesAbsorptionCoefficientsUpperAtmos', &
      (/1,1,1,1,bandNumber,gPointSetNumber/), &
      (/keyupper,Tdiff,pupper,numGPoints,1,1/), kbo)

    CALL read(fileid,'H2OSelfAbsorptionCoefficients', &
      (/1,1,bandNumber,gPointSetNumber/), &
      (/Tself,numGPoints,1,1/), selfrefo )

    CALL read(fileid,'H2OForeignAbsorptionCoefficientsLowerAtmos', &
      (/1,1,bandNumber,gPointSetNumber/), &
      (/Tforeignlower,numGPoints,1,1/),  forrefo(1:3,:))

    CALL read(fileid,'H2OForeignAbsorptionCoefficientsUpperAtmos', &
      (/2,1,bandNumber,gPointSetNumber/), &
      (/1,numGPoints,1,1/), forrefo(4,:))

    rayl = ncrayl(1)
    
  END SUBROUTINE sw_kgb21

  SUBROUTINE sw_kgb22

    USE mo_psrad_srtm_kgs,  ONLY : sfluxrefo=>sfluxrefo22, kao=>kao22, &
      kbo=>kbo22, selfrefo=>selfrefo22, forrefo=>forrefo22, rayl=>rayl22

    INTEGER, PARAMETER :: bandNumber = 7, numGPoints = no(7)

    REAL(kind=wp) :: ncrayl(1)

    CALL read(fileid,'SolarSourceFunctionLowerAtmos', &
      (/1,1,bandNumber,gPointSetNumber/), &
      (/numGPoints,keylower,1,1/), sfluxrefo)
    CALL read(fileid,'RayleighExtinctionCoefficientsLowerAtmos', &
      (/1,1,bandNumber,gPointSetNumber/), &
      (/1,1,1,1/), ncrayl)
    CALL read(fileid,'KeySpeciesAbsorptionCoefficientsLowerAtmos', &
      (/1,1,1,1,bandNumber,gPointSetNumber/), &
      (/keylower,Tdiff,plower,numGPoints,1,1/), kao)
    CALL read(fileid,'KeySpeciesAbsorptionCoefficientsUpperAtmos', &
      (/1,1,1,1,bandNumber,gPointSetNumber/), &
      (/1,Tdiff,pupper,numGPoints,1,1/), kbo)
    CALL read(fileid,'H2OSelfAbsorptionCoefficients', &
      (/1,1,bandNumber,gPointSetNumber/), &
      (/Tself,numGPoints,1,1/), selfrefo)
    CALL read(fileid,'H2OForeignAbsorptionCoefficientsLowerAtmos', &
      (/1,1,bandNumber,gPointSetNumber/), &
      (/Tforeignlower,numGPoints,1,1/), forrefo)

    rayl = ncrayl(1)

    
  END SUBROUTINE sw_kgb22

  SUBROUTINE sw_kgb23

    USE mo_psrad_srtm_kgs,  ONLY : sfluxrefo=>sfluxrefo23, kao=>kao23, &
      selfrefo=>selfrefo23, forrefo=>forrefo23, raylo=>raylo23

    INTEGER, PARAMETER :: bandNumber = 8, numGPoints = no(8)

    CALL read(fileid,'SolarSourceFunctionLowerAtmos', &
      (/1,1,bandNumber,gPointSetNumber/), &
      (/numGPoints,1,1,1/), sfluxrefo)

    CALL read(fileid,'RayleighExtinctionCoefficientsLowerAtmos', &
      (/1,1,bandNumber,gPointSetNumber/), &
      (/numGPoints,1,1,1/), raylo)

    CALL read(fileid,'KeySpeciesAbsorptionCoefficientsLowerAtmos', &
      (/1,1,1,1,bandNumber,gPointSetNumber/), &
      (/1,Tdiff,plower,numGPoints,1,1/), kao)

    CALL read(fileid,'H2OSelfAbsorptionCoefficients', &
      (/1,1,bandNumber,gPointSetNumber/), &
      (/Tself,numGPoints,1,1/), selfrefo)

    CALL read(fileid,'H2OForeignAbsorptionCoefficientsLowerAtmos', &
      (/1,1,bandNumber,gPointSetNumber/), &
      (/Tforeignlower,numGPoints,1,1/), forrefo)

  END SUBROUTINE sw_kgb23

  SUBROUTINE sw_kgb24

    USE mo_psrad_srtm_kgs,  ONLY : sfluxrefo=>sfluxrefo24, kao=>kao24, &
      kbo=>kbo24, selfrefo=>selfrefo24, forrefo=>forrefo24, raylao=>raylao24, &
      raylbo=>raylbo24, abso3ao=>abso3ao24, abso3bo=>abso3bo24

    INTEGER, PARAMETER :: bandNumber = 9, numGPoints = no(9)

    CALL read(fileid,'SolarSourceFunctionLowerAtmos', &
      (/1,1,bandNumber,gPointSetNumber/), &
      (/numGPoints,keylower,1,1/), sfluxrefo)

    CALL read(fileid,'RayleighExtinctionCoefficientsLowerAtmos', &
      (/1,1,bandNumber,gPointSetNumber/), &
      (/numGPoints,keylower,1,1/), raylao)

    CALL read(fileid,'RayleighExtinctionCoefficientsUpperAtmos', &
      (/1,1,bandNumber,gPointSetNumber/), &
      (/numGPoints,1,1,1/), raylbo)

    CALL read(fileid,'KeySpeciesAbsorptionCoefficientsLowerAtmos', &
      (/1,1,1,1,bandNumber,gPointSetNumber/), &
      (/keylower,Tdiff,plower,numGPoints,1,1/), kao)

    CALL read(fileid,'KeySpeciesAbsorptionCoefficientsUpperAtmos', &
      (/1,1,1,1,bandNumber,gPointSetNumber/), &
      (/1,Tdiff,pupper,numGPoints,1,1/), kbo)

    CALL read(fileid,'H2OSelfAbsorptionCoefficients', &
      (/1,1,bandNumber,gPointSetNumber/), &
      (/Tself,numGPoints,1,1/), selfrefo)

    CALL read(fileid,'H2OForeignAbsorptionCoefficientsLowerAtmos', &
      (/1,1,bandNumber,gPointSetNumber/), &
      (/Tforeignlower,numGPoints,1,1/), forrefo)

    CALL read(fileid,'AbsorptionCoefficientsLowerAtmos', &
      (/1,1,1,O3_idx,bandNumber,gPointSetNumber/), &
      (/1,1,numGPoints,1,1,1/), abso3ao)
    CALL read(fileid,'AbsorptionCoefficientsUpperAtmos', &
      (/1,1,1,O3_idx,bandNumber,gPointSetNumber/), &
      (/1,1,numGPoints,1,1,1/), abso3bo)
    
  END SUBROUTINE sw_kgb24

  SUBROUTINE sw_kgb25

    USE mo_psrad_srtm_kgs,  ONLY : sfluxrefo=>sfluxrefo25, kao=>kao25, &
      raylo=>raylo25, abso3ao=>abso3ao25, abso3bo=>abso3bo25

    INTEGER, PARAMETER :: bandNumber = 10, numGPoints = no(10)

    CALL read(fileid,'SolarSourceFunctionLowerAtmos', &
      (/1,1,bandNumber,gPointSetNumber/), &
      (/numGPoints,1,1,1/), sfluxrefo)

    CALL read(fileid,'RayleighExtinctionCoefficientsLowerAtmos', &
      (/1,1,bandNumber,gPointSetNumber/), &
      (/numGPoints,1,1,1/), raylo)

    CALL read(fileid,'KeySpeciesAbsorptionCoefficientsLowerAtmos', &
      (/1,1,1,1,bandNumber,gPointSetNumber/), &
      (/1,Tdiff,plower,numGPoints,1,1/), kao)

    CALL read(fileid,'AbsorptionCoefficientsLowerAtmos', &
      (/1,1,1,O3_idx,bandNumber,gPointSetNumber/), &
      (/1,1,numGPoints,1,1,1/), abso3ao)

    CALL read(fileid,'AbsorptionCoefficientsUpperAtmos', &
      (/1,1,1,O3_idx,bandNumber,gPointSetNumber/), &
      (/1,1,numGPoints,1,1,1/), abso3bo)

  END SUBROUTINE sw_kgb25

  SUBROUTINE sw_kgb26

    USE mo_psrad_srtm_kgs,  ONLY : sfluxrefo=>sfluxrefo26, raylo=>raylo26

    INTEGER, PARAMETER :: bandNumber = 11, numGPoints = no(11)

    CALL read(fileid,'SolarSourceFunctionLowerAtmos', &
      (/1,1,bandNumber,gPointSetNumber/), &
      (/numGPoints,1,1,1/), sfluxrefo)

    CALL read(fileid,'RayleighExtinctionCoefficientsLowerAtmos', &
      (/1,1,bandNumber,gPointSetNumber/), &
      (/numGPoints,1,1,1/), raylo)

  END SUBROUTINE sw_kgb26

  SUBROUTINE sw_kgb27

    USE mo_psrad_srtm_kgs,  ONLY : sfluxrefo=>sfluxrefo27, kao=>kao27, &
      kbo=>kbo27, raylo=>raylo27

    INTEGER, PARAMETER :: bandNumber = 12, numGPoints = no(12)

    CALL read(fileid,'SolarSourceFunctionLowerAtmos', &
      (/1,1,bandNumber,gPointSetNumber/), &
      (/numGPoints,1,1,1/), sfluxrefo)

    CALL read(fileid,'RayleighExtinctionCoefficientsLowerAtmos', &
      (/1,1,bandNumber,gPointSetNumber/), &
      (/numGPoints,1,1,1/), raylo)

    CALL read(fileid,'KeySpeciesAbsorptionCoefficientsLowerAtmos', &
      (/1,1,1,1,bandNumber,gPointSetNumber/), &
      (/1,Tdiff,plower,numGPoints,1,1/), kao)

    CALL read(fileid,'KeySpeciesAbsorptionCoefficientsUpperAtmos', &
      (/1,1,1,1,bandNumber,gPointSetNumber/), &
      (/1,Tdiff,pupper,numGPoints,1,1/), kbo)

  END SUBROUTINE sw_kgb27

  SUBROUTINE sw_kgb28

    USE mo_psrad_srtm_kgs,  ONLY : sfluxrefo=>sfluxrefo28, kao=>kao28, &
      kbo=>kbo28, rayl=>rayl28

    INTEGER, PARAMETER :: bandNumber = 13, numGPoints = no(13)  

    REAL(kind=wp) :: ncrayl(1)

    CALL read(fileid,'SolarSourceFunctionUpperAtmos', &
      (/1,1,bandNumber,gPointSetNumber/), &
      (/numGPoints,keyupper,1,1/), sfluxrefo)

    CALL read(fileid,'RayleighExtinctionCoefficientsLowerAtmos', &
      (/1,1,bandNumber,gPointSetNumber/), &
      (/1,1,1,1/), ncrayl)

    CALL read(fileid,'KeySpeciesAbsorptionCoefficientsLowerAtmos', &
      (/1,1,1,1,bandNumber,gPointSetNumber/), &
      (/keylower,Tdiff,plower,numGPoints,1,1/), kao)

    CALL read(fileid,'KeySpeciesAbsorptionCoefficientsUpperAtmos', &
      (/1,1,1,1,bandNumber,gPointSetNumber/), &
      (/keyupper,Tdiff,pupper,numGPoints,1,1/), kbo)

    rayl = ncrayl(1)

  END SUBROUTINE sw_kgb28

  SUBROUTINE sw_kgb29

    USE mo_psrad_srtm_kgs, ONLY : sfluxrefo=>sfluxrefo29, kao=>kao29, &
      kbo=>kbo29, selfrefo=>selfrefo29, forrefo=>forrefo29, &
      absh2oo=>absh2oo29, absco2o=>absco2o29, rayl=>rayl29

    INTEGER, PARAMETER :: bandNumber = 14, numGPoints = no(14)

    REAL(kind=wp) :: ncrayl(1)

    CALL read(fileid,'SolarSourceFunctionLowerAtmos', &
      (/1,1,bandNumber,gPointSetNumber/), &
      (/numGPoints,1,1,1/), sfluxrefo)

    CALL read(fileid,'RayleighExtinctionCoefficientsLowerAtmos', &
      (/1,1,bandNumber,gPointSetNumber/), &
      (/1,1,1,1/), ncrayl)

    CALL read(fileid,'KeySpeciesAbsorptionCoefficientsLowerAtmos', &
      (/1,1,1,1,bandNumber,gPointSetNumber/), &
      (/1,Tdiff,plower,numGPoints,1,1/), kao)

    CALL read(fileid,'KeySpeciesAbsorptionCoefficientsUpperAtmos', &
      (/1,1,1,1,bandNumber,gPointSetNumber/), &
      (/1,Tdiff,pupper,numGPoints,1,1/), kbo)

    CALL read(fileid,'H2OSelfAbsorptionCoefficients', &
      (/1,1,bandNumber,gPointSetNumber/), &
      (/Tself,numGPoints,1,1/), selfrefo)

    CALL read(fileid,'H2OForeignAbsorptionCoefficientsLowerAtmos', &
      (/1,1,bandNumber,gPointSetNumber/), &
      (/Tforeignlower,numGPoints,1,1/), forrefo(1:3,:))

    CALL read(fileid,'H2OForeignAbsorptionCoefficientsUpperAtmos', &
      (/2,1,bandNumber,gPointSetNumber/), &
      (/1,numGPoints,1,1/), forrefo(4,:))

    CALL read(fileid,'AbsorptionCoefficientsLowerAtmos', &
      (/1,1,1,H2O_idx,bandNumber,gPointSetNumber/), &
      (/1,1,numGPoints,1,1,1/), absh2oo)

    CALL read(fileid,'AbsorptionCoefficientsLowerAtmos', &
      (/1,1,1,CO2_idx,bandNumber,gPointSetNumber/), &
      (/1,1,numGPoints,1,1,1/), absco2o)

    rayl = ncrayl(1)

  END SUBROUTINE sw_kgb29

END MODULE mo_psrad_srtm_netcdf

