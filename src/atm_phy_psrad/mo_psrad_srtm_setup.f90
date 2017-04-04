!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_psrad_srtm_setup

  USE mo_psrad_general, ONLY : wp, nbndsw, mg, ngptsw
  USE mo_psrad_srtm_netcdf,    ONLY : srtm_read

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: ngb, ngc, ngs, nsp, ioff, mult, lay_ref, lay_sref_mask, &
    ssi_default, ssi_preind, ssi_amip, ssi_RCEdiurnOn, ssi_RCEdiurnOff, &
    ssi_cmip6_picontrol, &
    wavenum1, wavenum2, delwave, setup_srtm
  !     Arrays for the g-point reduction from 224 to 112 for the 16 LW bands:
  !     This mapping from 224 to 112 points has been carefully selected to 
  !     minimize the effect on the resulting fluxes and cooling rates, and
  !     caution should be used if the mapping is modified.  The full 224
  !     g-point set can be restored with ngpt=224, ngc=16*16, ngn=224*1., etc.
  ! ngpt - The total number of new g-points
  !       for each band.  
  !       new g-point in each band.
  ! wt      

  INTEGER, PARAMETER :: &
    ng(nbndsw) = (/16,16,16,16,16,16,16,16,16,16,16,16, 16, 16/), &
    !NOTE: nsp and ioff were set at runtime!!!
    nsp(2,nbndsw) = RESHAPE((/&
      9, 9, 9, 9, 1, 9, 9, 1, 9, 1, 0, 1, 9, 1, &
      1, 5, 1, 1, 1, 5, 1, 0, 1, 0, 0, 1, 5, 1/), &
      SHAPE = (/2,nbndsw/), ORDER = (/2,1/)), &
    ioff(2,nbndsw) = RESHAPE((/&
      9, 9, 9, 9, 1, 9, 9, 1, 9, 1, 1, 1, 9, 1, &
      1, 5, 1, 1, 1, 5, 1, 1, 1, 1, 1, 1, 5, 1/), &
      SHAPE = (/2,nbndsw/), ORDER = (/2,1/)), &

  ! ngc - The number of new g-points in each band
    ngc(nbndsw) = (/6,12, 8, 8,10,10, 2,10, 8, 6, 6, 8,  6, 12/), &
  ! ngs - The cumulative sum of new g-points for each band
    ngs(nbndsw) = (/6,18,26,34,44,54,56,66,74,80,86,94,100,112/), &
  ! ngm - The index of each new g-point relative to the original 16 g-points 
    ngm(nbndsw*mg) = (/ &
      1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, & ! band 16
      1, 2, 3, 4, 5, 6, 6, 7, 8, 8, 9,10,10,11,12,12, & ! band 17
      1, 2, 3, 4, 5, 5, 6, 6, 7, 7, 7, 7, 8, 8, 8, 8, & ! band 18
      1, 2, 3, 4, 5, 5, 6, 6, 7, 7, 7, 7, 8, 8, 8, 8, & ! band 19
      1, 2, 3, 4, 5, 6, 7, 8, 9, 9,10,10,10,10,10,10, & ! band 20
      1, 2, 3, 4, 5, 6, 7, 8, 9, 9,10,10,10,10,10,10, & ! band 21
      1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, & ! band 22
      1, 1, 2, 2, 3, 4, 5, 6, 7, 8, 9, 9,10,10,10,10, & ! band 23
      1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, & ! band 24
      1, 2, 3, 3, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, & ! band 25
      1, 2, 3, 3, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, & ! band 26
      1, 2, 3, 4, 5, 6, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, & ! band 27
      1, 2, 3, 3, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, & ! band 28
      1, 2, 3, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9,10,11,12 /), & ! band 29
  ! ngn - The number of original g-points that are combined to make each 
    ngn(ngptsw) = (/ &
       2,2,2,2,4,4,             & ! band 16
       1,1,1,1,1,2,1,2,1,2,1,2, & ! band 17
       1,1,1,1,2,2,4,4,         & ! band 18
       1,1,1,1,2,2,4,4,         & ! band 19
       1,1,1,1,1,1,1,1,2,6,     & ! band 20
       1,1,1,1,1,1,1,1,2,6,     & ! band 21
       8,8,                     & ! band 22
       2,2,1,1,1,1,1,1,2,4,     & ! band 23
       2,2,2,2,2,2,2,2,         & ! band 24
       1,1,2,2,4,6,             & ! band 25
       1,1,2,2,4,6,             & ! band 26
       1,1,1,1,1,1,4,6,         & ! band 27
       1,1,2,2,4,6,             & ! band 28
       1,1,1,1,2,2,2,2,1,1,1,1 /), & ! band 29
  ! ngb - The band index for each new g-point.
    ngb(ngptsw) = (/ &
       1,1,1,1,1,1,             & ! band 1
       2,2,2,2,2,2,2,2,2,2,2,2, & ! band 2
       3,3,3,3,3,3,3,3,         & ! band 3
       4,4,4,4,4,4,4,4,         & ! band 4
       5,5,5,5,5,5,5,5,5,5,     & ! band 5
       6,6,6,6,6,6,6,6,6,6,     & ! band 6
       7,7,                     & ! band 7
       8,8,8,8,8,8,8,8,8,8,     & ! band 8
       9,9,9,9,9,9,9,9,         & ! band 9
       10,10,10,10,10,10,       & ! band 10
       11,11,11,11,11,11,       & ! band 11
       12,12,12,12,12,12,12,12, & ! band 12
       13,13,13,13,13,13,       & ! band 13
       14,14,14,14,14,14,14,14,14,14,14,14 /) ! band 14

  ! Shortwave spectral band limits (wavenumbers)
  REAL(wp), PARAMETER :: wavenum1(nbndsw) = (/ &
       2600._wp, 3250._wp, 4000._wp, 4650._wp, 5150._wp, 6150._wp, 7700._wp, &
       8050._wp,12850._wp,16000._wp,22650._wp,29000._wp,38000._wp,  820._wp/)
  REAL(wp), PARAMETER :: wavenum2(nbndsw) = (/ &
       3250., 4000., 4650., 5150., 6150., 7700., 8050., &
       12850.,16000.,22650.,29000.,38000.,50000., 2600./), &
    delwave(nbndsw)  = (/ &
       650.,  750.,  650.,  500., 1000., 1550.,  350., &
       4800., 3150., 6650., 6350., 9000.,12000., 1780./), &
    mult(2,2,nbndsw) = RESHAPE((/&
      252.131_wp, 0.0_wp, 8.0_wp, 8.0_wp, &
      0.364641_wp, 0.364641_wp, 8.0_wp, 4.0_wp, &
      38.9589_wp, 0.0_wp, 8.0_wp, 0.0_wp, &
      5.49281_wp, 0.0_wp, 8.0_wp, 0.0_wp, &
      0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
      0.0045321_wp, 0.0045321_wp, 8.0_wp, 4.0_wp, &
      !NOTE: some awesome banana peel there
      0.022708_wp * 1.6_wp, 0.0_wp, 8.0_wp, 0.0_wp, &
      0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
      0.124692_wp, 0.0_wp, 8.0_wp, 0.0_wp, &
      0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
      0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
      0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
      6.67029e-07_wp, 6.67029e-07_wp, 8.0_wp, 4.0_wp, &
      0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp/), &
      SHAPE = (/2,2,nbndsw/), ORDER = (/2,1,3/)), &
    !NOTE: The zero was missing in the original!
    lay_ref(nbndsw) = (/-18, 30, 6, 3, -3, 8, 2, -6, 1, 0, -1, -32, 58, -49/)

  INTEGER, PARAMETER :: lay_sref_mask(nbndsw) = &
    (/1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1/)

  REAL(wp) :: wt(mg) = (/ & ! RRTM weights for 16 g-points.
       0.1527534276_wp, 0.1491729617_wp, 0.1420961469_wp, &
       0.1316886544_wp, 0.1181945205_wp, 0.1019300893_wp, &
       0.0832767040_wp, 0.0626720116_wp, 0.0424925000_wp, &
       0.0046269894_wp, 0.0038279891_wp, 0.0030260086_wp, &
       0.0022199750_wp, 0.0014140010_wp, 0.0005330000_wp, &
       0.0000750000_wp /)

  REAL(wp) :: rwgt(nbndsw*mg)

  REAL(wp), PARAMETER :: &
    ssi_default(14) =  (/ & !< SRTM default solar flux (W/m2) in 14 SW bands
      12.1095682699999987_wp, 20.3650825429849398_wp, &
      23.7297328613475429_wp, 22.4276934179066849_wp, &
      55.6266126199999960_wp, 1.0293153082385436E+02_wp, &
      24.2936128100154392_wp, 3.4574251380000004E+02_wp, &
      2.1818712729860400E+02_wp, 3.4719231470000005E+02_wp, &
      1.2949501812000000E+02_wp, 50.1522503011060650_wp, &
      3.0799387010047838_wp, 12.8893773299999985_wp/), &
      ! sum of 14 bands is: 1.3682223735968237E+03
    ssi_amip(14) =  (/ & ! solar flux (W/m2) in 14 SW bands for AMIP-type 
    !...CMIP5 simulation (average from 1979-1988)
      11.95053_wp, 20.14766_wp, 23.40394_wp, 22.09458_wp, 55.41401_wp, &
      102.5134_wp, 24.69814_wp, 347.5362_wp, 217.2925_wp, 343.4221_wp, &
      129.403_wp, 47.14264_wp, 3.172126_wp, 13.18075_wp /), &
    ! sum of 14 bands is: 1361.371
    ssi_preind(14) =  (/ & !< solar flux (W/m2) in 14 SW bands for 
    ! ... preindutrial CMIP5 simulation (average from 1944-1856)
      11.95005_wp, 20.14612_wp, 23.40302_wp, 22.09443_wp, 55.41679_wp, &
      102.512_wp , 24.69536_wp, 347.4719_wp, 217.2217_wp, 343.2816_wp, &
      129.3001_wp, 47.07624_wp, 3.130199_wp, 13.17521_wp /), &
    ! sum of 14 bands is: 1360.875

    ! ssi_RCEdiurnOn added (diurnal cycle on)
    ssi_RCEdiurnOn(14) =  (/ & !< solar flux (W/m2) in 14 SW bands for
    !... RCE simulations with diurnal cycle. global mean 
    ! insolation = 340.3 W/m2
      9.386766_wp,  15.82535_wp,  18.38306_wp,   17.3546_wp,  43.52597_wp,  &
      80.52106_wp,  19.39961_wp,  272.9788_wp,  170.6764_wp,  269.7473_wp,  &
      101.642_wp,   37.02906_wp,  2.491606_wp,  10.35307_wp/), &
    ! sum of 14 bands is: 1069.315 for diurnal cycle on.

    ! ssi_RCEdiurnOFF added (diurnal cycle off)
    ssi_RCEdiurnOFF(14) = (/ & !< solar flux (W/m2) in 14 SW bands for
    ! ... RCE simulations with diurnal cycle switched off. global mean 
    ! insolation = 340.3 W/m2 rescaled from ssi_amip above, with constant 
    ! factor of app. 0.3183092
      3.803964_wp,  6.413186_wp,  7.44969_wp,   7.032908_wp,  17.63879_wp, &
      32.63096_wp,  7.861645_wp,  110.624_wp,   69.16621_wp,  109.3144_wp, &
      41.19017_wp,  15.00594_wp, 1.009717_wp,   4.195554_wp  /), &
    ! sum of 14 bands is: 433.3371
    ssi_cmip6_picontrol(14) =  (/ & !< solar flux (W/m2) in 14 SW bands for
    ! ...preindustrial CMIP6 simulation (average from 1850-1873)
    & 12.02503_wp, 20.24537_wp, 23.69633_wp, 22.42093_wp, 55.91312_wp,  &
    & 103.5685_wp, 24.46918_wp, 346.3545_wp, 217.1642_wp, 344.9984_wp,  &
    & 127.7391_wp, 45.95287_wp, 2.957935_wp, 13.2384_wp /)
    ! sum of 14 bands is: 1360.744

CONTAINS

  SUBROUTINE setup_srtm
    ! Lookup tables are computed for use in the SW radiative transfer, 
    ! and input absorption coefficient data for each spectral band are 
    ! reduced from 224 g-point intervals to 112.

    INTEGER :: ibnd, igc, ig, ind, ipr
    INTEGER :: igcsm, iprsm
    REAL(wp) :: wtsum, wtsm(mg)

    CALL srtm_read

    ! Perform g-point reduction from 16 per band (224 total points) to
    ! a band dependent number (112 total points) for all absorption
    ! coefficient input data and Planck fraction input data.
    ! Compute relative weighting for new g-point combinations.

    igcsm = 0
    DO ibnd = 1,nbndsw
      iprsm = 0
      IF (ngc(ibnd).LT.mg) THEN
        DO igc = 1,ngc(ibnd)
          igcsm = igcsm + 1
          wtsum = 0.
          DO ipr = 1, ngn(igcsm)
            iprsm = iprsm + 1
            wtsum = wtsum + wt(iprsm)
          ENDDO
          wtsm(igc) = wtsum
        ENDDO
        DO ig = 1, ng(ibnd)
          ind = (ibnd-1)*mg + ig
          rwgt(ind) = wt(ig)/wtsm(ngm(ind))
        ENDDO
      ELSE
        DO ig = 1, ng(ibnd)
          igcsm = igcsm + 1
          ind = (ibnd-1)*mg + ig
          rwgt(ind) = 1.0_wp
        ENDDO
      ENDIF
    ENDDO

    ! Reduce g-points for absorption coefficient data in each LW spectral band.
    CALL cmbgb16s
    CALL cmbgb17
    CALL cmbgb18
    CALL cmbgb19
    CALL cmbgb20
    CALL cmbgb21
    CALL cmbgb22
    CALL cmbgb23
    CALL cmbgb24
    CALL cmbgb25
    CALL cmbgb26
    CALL cmbgb27
    CALL cmbgb28
    CALL cmbgb29

  END SUBROUTINE setup_srtm

  !  The subroutines CMBGB16->CMBGB29 input the absorption coefficient
  !  data for each band, which are defined for 16 g-points and 14 spectral
  !  bands. The data are combined with appropriate weighting following the
  !  g-point mapping arrays specified in RRTMG_SW_INIT.  Solar source 
  !  function data in array SFLUXREF are combined without weighting.  All
  !  g-point reduced data are put into new arrays for use in RRTMG_SW.
  !

  !  band 16:  2600-3250 cm-1 (low key- h2o,ch4; high key - ch4)
  SUBROUTINE cmbgb16s

    USE mo_psrad_srtm_kgs, ONLY : kao=>kao16, kbo=>kbo16, selfrefo=>selfrefo16, forrefo=>forrefo16, sfluxrefo=>sfluxrefo16, &
         ka=>ka16, kb=>kb16, selfref=>selfref16, forref=>forref16, sfluxref=>sfluxref16

    INTEGER :: jn, jt, jp, igc, ipr, iprsm
    REAL(wp) :: sumk, sumf


    DO jn = 1,9
      DO jt = 1,5
        DO jp = 1,13
          iprsm = 0
          DO igc = 1,ngc(1)
            sumk = 0.
            DO ipr = 1, ngn(igc)
              iprsm = iprsm + 1
              sumk = sumk + kao(jn,jt,jp,iprsm)*rwgt(iprsm)
            ENDDO
            ka(jn,jt,jp,igc) = sumk
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    DO jt = 1,5
      DO jp = 1,47
        iprsm = 0
        DO igc = 1,ngc(1)
          sumk = 0.
          DO ipr = 1, ngn(igc)
            iprsm = iprsm + 1
            sumk = sumk + kbo(jt,jp,iprsm)*rwgt(iprsm)
          ENDDO
          kb(jt,jp,igc) = sumk
        ENDDO
      ENDDO
    ENDDO

    DO jt = 1,10
      iprsm = 0
      DO igc = 1,ngc(1)
        sumk = 0.
        DO ipr = 1, ngn(igc)
          iprsm = iprsm + 1
          sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm)
        ENDDO
        selfref(jt,igc) = sumk
      ENDDO
    ENDDO

    DO jt = 1,3
      iprsm = 0
      DO igc = 1,ngc(1)
        sumk = 0.
        DO ipr = 1, ngn(igc)
          iprsm = iprsm + 1
          sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm)
        ENDDO
        forref(jt,igc) = sumk
      ENDDO
    ENDDO

    iprsm = 0
    DO igc = 1,ngc(1)
      sumf = 0.
      DO ipr = 1, ngn(igc)
        iprsm = iprsm + 1
        sumf = sumf + sfluxrefo(iprsm)
      ENDDO
      sfluxref(igc) = sumf
    ENDDO

  END SUBROUTINE cmbgb16s

  !     band 17:  3250-4000 cm-1 (low - h2o,co2; high - h2o,co2)
  SUBROUTINE cmbgb17

    USE mo_psrad_srtm_kgs, ONLY : kao=>kao17, kbo=>kbo17, selfrefo=>selfrefo17, forrefo=>forrefo17, sfluxrefo=>sfluxrefo17, &
         ka=>ka17, kb=>kb17, selfref=>selfref17, forref=>forref17, sfluxref=>sfluxref17

    INTEGER :: jn, jt, jp, igc, ipr, iprsm
    REAL(wp) :: sumk, sumf

    DO jn = 1,9
      DO jt = 1,5
        DO jp = 1,13
          iprsm = 0
          DO igc = 1,ngc(2)
            sumk = 0.
            DO ipr = 1, ngn(ngs(1)+igc)
              iprsm = iprsm + 1
              sumk = sumk + kao(jn,jt,jp,iprsm)*rwgt(iprsm+16)
            ENDDO
            ka(jn,jt,jp,igc) = sumk
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    DO jn = 1,5
      DO jt = 1,5
        DO jp = 1,47
          iprsm = 0
          DO igc = 1,ngc(2)
            sumk = 0.
            DO ipr = 1, ngn(ngs(1)+igc)
              iprsm = iprsm + 1
              sumk = sumk + kbo(jn,jt,jp,iprsm)*rwgt(iprsm+16)
            ENDDO
            kb(jn,jt,jp,igc) = sumk
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    DO jt = 1,10
      iprsm = 0
      DO igc = 1,ngc(2)
        sumk = 0.
        DO ipr = 1, ngn(ngs(1)+igc)
          iprsm = iprsm + 1
          sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+16)
        ENDDO
        selfref(jt,igc) = sumk
      ENDDO
    ENDDO

    DO jt = 1,4
      iprsm = 0
      DO igc = 1,ngc(2)
        sumk = 0.
        DO ipr = 1, ngn(ngs(1)+igc)
          iprsm = iprsm + 1
          sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+16)
        ENDDO
        forref(jt,igc) = sumk
      ENDDO
    ENDDO

    DO jp = 1,5
      iprsm = 0
      DO igc = 1,ngc(2)
        sumf = 0.
        DO ipr = 1, ngn(ngs(1)+igc)
          iprsm = iprsm + 1
          sumf = sumf + sfluxrefo(iprsm,jp)
        ENDDO
        sfluxref(igc,jp) = sumf
      ENDDO
    ENDDO

  END SUBROUTINE cmbgb17

  !     band 18:  4000-4650 cm-1 (low - h2o,ch4; high - ch4)
  SUBROUTINE cmbgb18

    USE mo_psrad_srtm_kgs, ONLY : kao=>kao18, kbo=>kbo18, selfrefo=>selfrefo18, forrefo=>forrefo18, sfluxrefo=>sfluxrefo18, &
         ka=>ka18, kb=>kb18, selfref=>selfref18, forref=>forref18, sfluxref=>sfluxref18

    INTEGER :: jn, jt, jp, igc, ipr, iprsm
    REAL(wp) :: sumk, sumf

    DO jn = 1,9
      DO jt = 1,5
        DO jp = 1,13
          iprsm = 0
          DO igc = 1,ngc(3)
            sumk = 0.
            DO ipr = 1, ngn(ngs(2)+igc)
              iprsm = iprsm + 1
              sumk = sumk + kao(jn,jt,jp,iprsm)*rwgt(iprsm+32)
            ENDDO
            ka(jn,jt,jp,igc) = sumk
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    DO jt = 1,5
      DO jp = 1,47
        iprsm = 0
        DO igc = 1,ngc(3)
          sumk = 0.
          DO ipr = 1, ngn(ngs(2)+igc)
            iprsm = iprsm + 1
            sumk = sumk + kbo(jt,jp,iprsm)*rwgt(iprsm+32)
          ENDDO
          kb(jt,jp,igc) = sumk
        ENDDO
      ENDDO
    ENDDO

    DO jt = 1,10
      iprsm = 0
      DO igc = 1,ngc(3)
        sumk = 0.
        DO ipr = 1, ngn(ngs(2)+igc)
          iprsm = iprsm + 1
          sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+32)
        ENDDO
        selfref(jt,igc) = sumk
      ENDDO
    ENDDO

    DO jt = 1,3
      iprsm = 0
      DO igc = 1,ngc(3)
        sumk = 0.
        DO ipr = 1, ngn(ngs(2)+igc)
          iprsm = iprsm + 1
          sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+32)
        ENDDO
        forref(jt,igc) = sumk
      ENDDO
    ENDDO

    DO jp = 1,9
      iprsm = 0
      DO igc = 1,ngc(3)
        sumf = 0.
        DO ipr = 1, ngn(ngs(2)+igc)
          iprsm = iprsm + 1
          sumf = sumf + sfluxrefo(iprsm,jp)
        ENDDO
        sfluxref(igc,jp) = sumf
      ENDDO
    ENDDO

  END SUBROUTINE cmbgb18

  !     band 19:  4650-5150 cm-1 (low - h2o,co2; high - co2)
  SUBROUTINE cmbgb19

    USE mo_psrad_srtm_kgs, ONLY : kao=>kao19, kbo=>kbo19, selfrefo=>selfrefo19, forrefo=>forrefo19, sfluxrefo=>sfluxrefo19, &
         ka=>ka19, kb=>kb19, selfref=>selfref19, forref=>forref19, sfluxref=>sfluxref19

    INTEGER :: jn, jt, jp, igc, ipr, iprsm
    REAL(wp) :: sumk, sumf

    DO jn = 1,9
      DO jt = 1,5
        DO jp = 1,13
          iprsm = 0
          DO igc = 1,ngc(4)
            sumk = 0.
            DO ipr = 1, ngn(ngs(3)+igc)
              iprsm = iprsm + 1
              sumk = sumk + kao(jn,jt,jp,iprsm)*rwgt(iprsm+48)
            ENDDO
            ka(jn,jt,jp,igc) = sumk
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    DO jt = 1,5
      DO jp = 1,47
        iprsm = 0
        DO igc = 1,ngc(4)
          sumk = 0.
          DO ipr = 1, ngn(ngs(3)+igc)
            iprsm = iprsm + 1
            sumk = sumk + kbo(jt,jp,iprsm)*rwgt(iprsm+48)
          ENDDO
          kb(jt,jp,igc) = sumk
        ENDDO
      ENDDO
    ENDDO

    DO jt = 1,10
      iprsm = 0
      DO igc = 1,ngc(4)
        sumk = 0.
        DO ipr = 1, ngn(ngs(3)+igc)
          iprsm = iprsm + 1
          sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+48)
        ENDDO
        selfref(jt,igc) = sumk
      ENDDO
    ENDDO

    DO jt = 1,3
      iprsm = 0
      DO igc = 1,ngc(4)
        sumk = 0.
        DO ipr = 1, ngn(ngs(3)+igc)
          iprsm = iprsm + 1
          sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+48)
        ENDDO
        forref(jt,igc) = sumk
      ENDDO
    ENDDO

    DO jp = 1,9
      iprsm = 0
      DO igc = 1,ngc(4)
        sumf = 0.
        DO ipr = 1, ngn(ngs(3)+igc)
          iprsm = iprsm + 1
          sumf = sumf + sfluxrefo(iprsm,jp)
        ENDDO
        sfluxref(igc,jp) = sumf
      ENDDO
    ENDDO

  END SUBROUTINE cmbgb19

  !     band 20:  5150-6150 cm-1 (low - h2o; high - h2o)
  SUBROUTINE cmbgb20

    USE mo_psrad_srtm_kgs, ONLY : kao=>kao20, kbo=>kbo20, selfrefo=>selfrefo20, forrefo=>forrefo20, sfluxrefo=>sfluxrefo20, &
      absch4o=>absch4o20, ka=>ka20, kb=>kb20, selfref=>selfref20, forref=>forref20, sfluxref=>sfluxref20, absch4=>absch420

    INTEGER :: jt, jp, igc, ipr, iprsm
    REAL(wp) :: sumk, sumf1, sumf2

    DO jt = 1,5
      DO jp = 1,13
        iprsm = 0
        DO igc = 1,ngc(5)
          sumk = 0.
          DO ipr = 1, ngn(ngs(4)+igc)
            iprsm = iprsm + 1
            sumk = sumk + kao(jt,jp,iprsm)*rwgt(iprsm+64)
          ENDDO
          ka(jt,jp,igc) = sumk
        ENDDO
      ENDDO
      DO jp = 1,47
        iprsm = 0
        DO igc = 1,ngc(5)
          sumk = 0.
          DO ipr = 1, ngn(ngs(4)+igc)
            iprsm = iprsm + 1
            sumk = sumk + kbo(jt,jp,iprsm)*rwgt(iprsm+64)
          ENDDO
          kb(jt,jp,igc) = sumk
        ENDDO
      ENDDO
    ENDDO

    DO jt = 1,10
      iprsm = 0
      DO igc = 1,ngc(5)
        sumk = 0.
        DO ipr = 1, ngn(ngs(4)+igc)
          iprsm = iprsm + 1
          sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+64)
        ENDDO
        selfref(jt,igc) = sumk
      ENDDO
    ENDDO

    DO jt = 1,4
      iprsm = 0
      DO igc = 1,ngc(5)
        sumk = 0.
        DO ipr = 1, ngn(ngs(4)+igc)
          iprsm = iprsm + 1
          sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+64)
        ENDDO
        forref(jt,igc) = sumk
      ENDDO
    ENDDO

    iprsm = 0
    DO igc = 1,ngc(5)
      sumf1 = 0.
      sumf2 = 0.
      DO ipr = 1, ngn(ngs(4)+igc)
        iprsm = iprsm + 1
        sumf1 = sumf1 + sfluxrefo(iprsm)
        sumf2 = sumf2 + absch4o(iprsm)*rwgt(iprsm+64)
      ENDDO
      sfluxref(igc) = sumf1
      absch4(igc) = sumf2
    ENDDO

  END SUBROUTINE cmbgb20

  !     band 21:  6150-7700 cm-1 (low - h2o,co2; high - h2o,co2)
  SUBROUTINE cmbgb21

    USE mo_psrad_srtm_kgs, ONLY : kao=>kao21, kbo=>kbo21, selfrefo=>selfrefo21, forrefo=>forrefo21, sfluxrefo=>sfluxrefo21, &
         ka=>ka21, kb=>kb21, selfref=>selfref21, forref=>forref21, sfluxref=>sfluxref21

    INTEGER :: jn, jt, jp, igc, ipr, iprsm
    REAL(wp) :: sumk, sumf


    DO jn = 1,9
      DO jt = 1,5
        DO jp = 1,13
          iprsm = 0
          DO igc = 1,ngc(6)
            sumk = 0.
            DO ipr = 1, ngn(ngs(5)+igc)
              iprsm = iprsm + 1
              sumk = sumk + kao(jn,jt,jp,iprsm)*rwgt(iprsm+80)
            ENDDO
            ka(jn,jt,jp,igc) = sumk
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    DO jn = 1,5
      DO jt = 1,5
        DO jp = 1,47
          iprsm = 0
          DO igc = 1,ngc(6)
            sumk = 0.
            DO ipr = 1, ngn(ngs(5)+igc)
              iprsm = iprsm + 1
              sumk = sumk + kbo(jn,jt,jp,iprsm)*rwgt(iprsm+80)
            ENDDO
            kb(jn,jt,jp,igc) = sumk
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    DO jt = 1,10
      iprsm = 0
      DO igc = 1,ngc(6)
        sumk = 0.
        DO ipr = 1, ngn(ngs(5)+igc)
          iprsm = iprsm + 1
          sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+80)
        ENDDO
        selfref(jt,igc) = sumk
      ENDDO
    ENDDO

    DO jt = 1,4
      iprsm = 0
      DO igc = 1,ngc(6)
        sumk = 0.
        DO ipr = 1, ngn(ngs(5)+igc)
          iprsm = iprsm + 1
          sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+80)
        ENDDO
        forref(jt,igc) = sumk
      ENDDO
    ENDDO

    DO jp = 1,9
      iprsm = 0
      DO igc = 1,ngc(6)
        sumf = 0.
        DO ipr = 1, ngn(ngs(5)+igc)
          iprsm = iprsm + 1
          sumf = sumf + sfluxrefo(iprsm,jp)
        ENDDO
        sfluxref(igc,jp) = sumf
      ENDDO
    ENDDO

  END SUBROUTINE cmbgb21

  !     band 22:  7700-8050 cm-1 (low - h2o,o2; high - o2)
  SUBROUTINE cmbgb22

    USE mo_psrad_srtm_kgs, ONLY : kao=>kao22, kbo=>kbo22, selfrefo=>selfrefo22, forrefo=>forrefo22, sfluxrefo=>sfluxrefo22, &
         ka=>ka22, kb=>kb22, selfref=>selfref22, forref=>forref22, sfluxref=>sfluxref22

    INTEGER :: jn, jt, jp, igc, ipr, iprsm
    REAL(wp) :: sumk, sumf

    DO jn = 1,9
      DO jt = 1,5
        DO jp = 1,13
          iprsm = 0
          DO igc = 1,ngc(7)
            sumk = 0.
            DO ipr = 1, ngn(ngs(6)+igc)
              iprsm = iprsm + 1
              sumk = sumk + kao(jn,jt,jp,iprsm)*rwgt(iprsm+96)
            ENDDO
            ka(jn,jt,jp,igc) = sumk
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    DO jt = 1,5
      DO jp = 1,47
        iprsm = 0
        DO igc = 1,ngc(7)
          sumk = 0.
          DO ipr = 1, ngn(ngs(6)+igc)
            iprsm = iprsm + 1
            sumk = sumk + kbo(jt,jp,iprsm)*rwgt(iprsm+96)
          ENDDO
          kb(jt,jp,igc) = sumk
        ENDDO
      ENDDO
    ENDDO

    DO jt = 1,10
      iprsm = 0
      DO igc = 1,ngc(7)
        sumk = 0.
        DO ipr = 1, ngn(ngs(6)+igc)
          iprsm = iprsm + 1
          sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+96)
        ENDDO
        selfref(jt,igc) = sumk
      ENDDO
    ENDDO

    DO jt = 1,3
      iprsm = 0
      DO igc = 1,ngc(7)
        sumk = 0.
        DO ipr = 1, ngn(ngs(6)+igc)
          iprsm = iprsm + 1
          sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+96)
        ENDDO
        forref(jt,igc) = sumk
      ENDDO
    ENDDO

    DO jp = 1,9
      iprsm = 0
      DO igc = 1,ngc(7)
        sumf = 0.
        DO ipr = 1, ngn(ngs(6)+igc)
          iprsm = iprsm + 1
          sumf = sumf + sfluxrefo(iprsm,jp)
        ENDDO
        sfluxref(igc,jp) = sumf
      ENDDO
    ENDDO

  END SUBROUTINE cmbgb22

  !     band 23:  8050-12850 cm-1 (low - h2o; high - nothing)
  SUBROUTINE cmbgb23

    USE mo_psrad_srtm_kgs, ONLY : kao=>kao23, selfrefo=>selfrefo23, forrefo=>forrefo23, sfluxrefo=>sfluxrefo23, raylo=>raylo23, &
         ka=>ka23, selfref=>selfref23, forref=>forref23, sfluxref=>sfluxref23, rayl=>rayl23

    INTEGER :: jt, jp, igc, ipr, iprsm
    REAL(wp) :: sumk, sumf1, sumf2

    DO jt = 1,5
      DO jp = 1,13
        iprsm = 0
        DO igc = 1,ngc(8)
          sumk = 0.
          DO ipr = 1, ngn(ngs(7)+igc)
            iprsm = iprsm + 1
            sumk = sumk + kao(jt,jp,iprsm)*rwgt(iprsm+112)
          ENDDO
          ka(jt,jp,igc) = sumk
        ENDDO
      ENDDO
    ENDDO

    DO jt = 1,10
      iprsm = 0
      DO igc = 1,ngc(8)
        sumk = 0.
        DO ipr = 1, ngn(ngs(7)+igc)
          iprsm = iprsm + 1
          sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+112)
        ENDDO
        selfref(jt,igc) = sumk
      ENDDO
    ENDDO

    DO jt = 1,3
      iprsm = 0
      DO igc = 1,ngc(8)
        sumk = 0.
        DO ipr = 1, ngn(ngs(7)+igc)
          iprsm = iprsm + 1
          sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+112)
        ENDDO
        forref(jt,igc) = sumk
      ENDDO
    ENDDO

    iprsm = 0
    DO igc = 1,ngc(8)
      sumf1 = 0.
      sumf2 = 0.
      DO ipr = 1, ngn(ngs(7)+igc)
        iprsm = iprsm + 1
        sumf1 = sumf1 + sfluxrefo(iprsm)
        sumf2 = sumf2 + raylo(iprsm)*rwgt(iprsm+112)
      ENDDO
      sfluxref(igc) = sumf1
      rayl(igc) = sumf2
    ENDDO

  END SUBROUTINE cmbgb23

  !     band 24:  12850-16000 cm-1 (low - h2o,o2; high - o2)
  SUBROUTINE cmbgb24

    USE mo_psrad_srtm_kgs, ONLY : kao=>kao24, kbo=>kbo24, selfrefo=>selfrefo24, forrefo=>forrefo24, sfluxrefo=>sfluxrefo24, &
         abso3ao=>abso3ao24, abso3bo=>abso3bo24, raylao=>raylao24, raylbo=>raylbo24, &
         ka=>ka24, kb=>kb24, selfref=>selfref24, forref=>forref24, sfluxref=>sfluxref24, &
         abso3a=>abso3a24, abso3b=>abso3b24, rayla=>rayla24, raylb=>raylb24

    INTEGER :: jn, jt, jp, igc, ipr, iprsm
    REAL(wp) :: sumk, sumf1, sumf2, sumf3

    DO jn = 1,9
      DO jt = 1,5
        DO jp = 1,13
          iprsm = 0
          DO igc = 1,ngc(9)
            sumk = 0.
            DO ipr = 1, ngn(ngs(8)+igc)
              iprsm = iprsm + 1
              sumk = sumk + kao(jn,jt,jp,iprsm)*rwgt(iprsm+128)
            ENDDO
            ka(jn,jt,jp,igc) = sumk
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    DO jt = 1,5
      DO jp = 1,47
        iprsm = 0
        DO igc = 1,ngc(9)
          sumk = 0.
          DO ipr = 1, ngn(ngs(8)+igc)
            iprsm = iprsm + 1
            sumk = sumk + kbo(jt,jp,iprsm)*rwgt(iprsm+128)
          ENDDO
          kb(jt,jp,igc) = sumk
        ENDDO
      ENDDO
    ENDDO

    DO jt = 1,10
      iprsm = 0
      DO igc = 1,ngc(9)
        sumk = 0.
        DO ipr = 1, ngn(ngs(8)+igc)
          iprsm = iprsm + 1
          sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+128)
        ENDDO
        selfref(jt,igc) = sumk
      ENDDO
    ENDDO

    DO jt = 1,3
      iprsm = 0
      DO igc = 1,ngc(9)
        sumk = 0.
        DO ipr = 1, ngn(ngs(8)+igc)
          iprsm = iprsm + 1
          sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+128)
        ENDDO
        forref(jt,igc) = sumk
      ENDDO
    ENDDO

    iprsm = 0
    DO igc = 1,ngc(9)
      sumf1 = 0.
      sumf2 = 0.
      sumf3 = 0.
      DO ipr = 1, ngn(ngs(8)+igc)
        iprsm = iprsm + 1
        sumf1 = sumf1 + raylbo(iprsm)*rwgt(iprsm+128)
        sumf2 = sumf2 + abso3ao(iprsm)*rwgt(iprsm+128)
        sumf3 = sumf3 + abso3bo(iprsm)*rwgt(iprsm+128)
      ENDDO
      raylb(igc) = sumf1
      abso3a(igc) = sumf2
      abso3b(igc) = sumf3
    ENDDO

    DO jp = 1,9
      iprsm = 0
      DO igc = 1,ngc(9)
        sumf1 = 0.
        sumf2 = 0.
        DO ipr = 1, ngn(ngs(8)+igc)
          iprsm = iprsm + 1
          sumf1 = sumf1 + sfluxrefo(iprsm,jp)
          sumf2 = sumf2 + raylao(iprsm,jp)*rwgt(iprsm+128)
        ENDDO
        sfluxref(igc,jp) = sumf1
        rayla(igc,jp) = sumf2
      ENDDO
    ENDDO

  END SUBROUTINE cmbgb24

  !     band 25:  16000-22650 cm-1 (low - h2o; high - nothing)
  SUBROUTINE cmbgb25

    USE mo_psrad_srtm_kgs, ONLY : kao=>kao25, sfluxrefo=>sfluxrefo25, &
         abso3ao=>abso3ao25, abso3bo=>abso3bo25, raylo=>raylo25, &
         ka=>ka25, sfluxref=>sfluxref25, &
         abso3a=>abso3a25, abso3b=>abso3b25, rayl=>rayl25

    INTEGER :: jt, jp, igc, ipr, iprsm
    REAL(wp) :: sumk, sumf1, sumf2, sumf3, sumf4

    DO jt = 1,5
      DO jp = 1,13
        iprsm = 0
        DO igc = 1,ngc(10)
          sumk = 0.
          DO ipr = 1, ngn(ngs(9)+igc)
            iprsm = iprsm + 1
            sumk = sumk + kao(jt,jp,iprsm)*rwgt(iprsm+144)
          ENDDO
          ka(jt,jp,igc) = sumk
        ENDDO
      ENDDO
    ENDDO

    iprsm = 0
    DO igc = 1,ngc(10)
      sumf1 = 0.
      sumf2 = 0.
      sumf3 = 0.
      sumf4 = 0.
      DO ipr = 1, ngn(ngs(9)+igc)
        iprsm = iprsm + 1
        sumf1 = sumf1 + sfluxrefo(iprsm)
        sumf2 = sumf2 + abso3ao(iprsm)*rwgt(iprsm+144)
        sumf3 = sumf3 + abso3bo(iprsm)*rwgt(iprsm+144)
        sumf4 = sumf4 + raylo(iprsm)*rwgt(iprsm+144)
      ENDDO
      sfluxref(igc) = sumf1
      abso3a(igc) = sumf2
      abso3b(igc) = sumf3
      rayl(igc) = sumf4
    ENDDO

  END SUBROUTINE cmbgb25

  !     band 26:  22650-29000 cm-1 (low - nothing; high - nothing)
  SUBROUTINE cmbgb26

    USE mo_psrad_srtm_kgs, ONLY : sfluxrefo=>sfluxrefo26, raylo=>raylo26, &
         sfluxref=>sfluxref26, rayl=>rayl26

    INTEGER :: igc, ipr, iprsm
    REAL(wp) :: sumf1, sumf2

    iprsm = 0
    DO igc = 1,ngc(11)
      sumf1 = 0.
      sumf2 = 0.
      DO ipr = 1, ngn(ngs(10)+igc)
        iprsm = iprsm + 1
        sumf1 = sumf1 + raylo(iprsm)*rwgt(iprsm+160)
        sumf2 = sumf2 + sfluxrefo(iprsm)
      ENDDO
      rayl(igc) = sumf1
      sfluxref(igc) = sumf2
    ENDDO

  END SUBROUTINE cmbgb26

  !     band 27:  29000-38000 cm-1 (low - o3; high - o3)
  SUBROUTINE cmbgb27

    USE mo_psrad_srtm_kgs, ONLY : kao=>kao27, kbo=>kbo27, sfluxrefo=>sfluxrefo27, raylo=>raylo27, &
         ka=>ka27, kb=>kb27, sfluxref=>sfluxref27, rayl=>rayl27

    INTEGER :: jt, jp, igc, ipr, iprsm
    REAL(wp) :: sumk, sumf1, sumf2

    DO jt = 1,5
      DO jp = 1,13
        iprsm = 0
        DO igc = 1,ngc(12)
          sumk = 0.
          DO ipr = 1, ngn(ngs(11)+igc)
            iprsm = iprsm + 1
            sumk = sumk + kao(jt,jp,iprsm)*rwgt(iprsm+176)
          ENDDO
          ka(jt,jp,igc) = sumk
        ENDDO
      ENDDO
      DO jp = 1,47
        iprsm = 0
        DO igc = 1,ngc(12)
          sumk = 0.
          DO ipr = 1, ngn(ngs(11)+igc)
            iprsm = iprsm + 1
            sumk = sumk + kbo(jt,jp,iprsm)*rwgt(iprsm+176)
          ENDDO
          kb(jt,jp,igc) = sumk
        ENDDO
      ENDDO
    ENDDO

    iprsm = 0
    DO igc = 1,ngc(12)
      sumf1 = 0.
      sumf2 = 0.
      DO ipr = 1, ngn(ngs(11)+igc)
        iprsm = iprsm + 1
        sumf1 = sumf1 + sfluxrefo(iprsm)
        sumf2 = sumf2 + raylo(iprsm)*rwgt(iprsm+176)
      ENDDO
      sfluxref(igc) = sumf1
      rayl(igc) = sumf2
    ENDDO

  END SUBROUTINE cmbgb27

  !     band 28:  38000-50000 cm-1 (low - o3,o2; high - o3,o2)
  SUBROUTINE cmbgb28

    USE mo_psrad_srtm_kgs, ONLY : kao=>kao28, kbo=>kbo28, sfluxrefo=>sfluxrefo28, &
         ka=>ka28, kb=>kb28, sfluxref=>sfluxref28

    INTEGER :: jn, jt, jp, igc, ipr, iprsm
    REAL(wp) :: sumk, sumf


    DO jn = 1,9
      DO jt = 1,5
        DO jp = 1,13
          iprsm = 0
          DO igc = 1,ngc(13)
            sumk = 0.
            DO ipr = 1, ngn(ngs(12)+igc)
              iprsm = iprsm + 1
              sumk = sumk + kao(jn,jt,jp,iprsm)*rwgt(iprsm+192)
            ENDDO
            ka(jn,jt,jp,igc) = sumk
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    DO jn = 1,5
      DO jt = 1,5
        DO jp = 1,47
          iprsm = 0
          DO igc = 1,ngc(13)
            sumk = 0.
            DO ipr = 1, ngn(ngs(12)+igc)
              iprsm = iprsm + 1
              sumk = sumk + kbo(jn,jt,jp,iprsm)*rwgt(iprsm+192)
            ENDDO
            kb(jn,jt,jp,igc) = sumk
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    DO jp = 1,5
      iprsm = 0
      DO igc = 1,ngc(13)
        sumf = 0.
        DO ipr = 1, ngn(ngs(12)+igc)
          iprsm = iprsm + 1
          sumf = sumf + sfluxrefo(iprsm,jp)
        ENDDO
        sfluxref(igc,jp) = sumf
      ENDDO
    ENDDO

  END SUBROUTINE cmbgb28

  !     band 29:  820-2600 cm-1 (low - h2o; high - co2)
  SUBROUTINE cmbgb29

    USE mo_psrad_srtm_kgs, ONLY : kao=>kao29, kbo=>kbo29, selfrefo=>selfrefo29, forrefo=>forrefo29, sfluxrefo=>sfluxrefo29, &
         absh2oo=>absh2oo29, absco2o=>absco2o29, ka=>ka29, kb=>kb29, &
         selfref=>selfref29, forref=>forref29, sfluxref=>sfluxref29, absh2o=>absh2o29, absco2=>absco229

    INTEGER :: jt, jp, igc, ipr, iprsm
    REAL(wp) :: sumk, sumf1, sumf2, sumf3

    DO jt = 1,5
      DO jp = 1,13
        iprsm = 0
        DO igc = 1,ngc(14)
          sumk = 0.
          DO ipr = 1, ngn(ngs(13)+igc)
            iprsm = iprsm + 1
            sumk = sumk + kao(jt,jp,iprsm)*rwgt(iprsm+208)
          ENDDO
          ka(jt,jp,igc) = sumk
        ENDDO
      ENDDO
      DO jp = 1,47
        iprsm = 0
        DO igc = 1,ngc(14)
          sumk = 0.
          DO ipr = 1, ngn(ngs(13)+igc)
            iprsm = iprsm + 1
            sumk = sumk + kbo(jt,jp,iprsm)*rwgt(iprsm+208)
          ENDDO
          kb(jt,jp,igc) = sumk
        ENDDO
      ENDDO
    ENDDO

    DO jt = 1,10
      iprsm = 0
      DO igc = 1,ngc(14)
        sumk = 0.
        DO ipr = 1, ngn(ngs(13)+igc)
          iprsm = iprsm + 1
          sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+208)
        ENDDO
        selfref(jt,igc) = sumk
      ENDDO
    ENDDO

    DO jt = 1,4
      iprsm = 0
      DO igc = 1,ngc(14)
        sumk = 0.
        DO ipr = 1, ngn(ngs(13)+igc)
          iprsm = iprsm + 1
          sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+208)
        ENDDO
        forref(jt,igc) = sumk
      ENDDO
    ENDDO

    iprsm = 0
    DO igc = 1,ngc(14)
      sumf1 = 0.
      sumf2 = 0.
      sumf3 = 0.
      DO ipr = 1, ngn(ngs(13)+igc)
        iprsm = iprsm + 1
        sumf1 = sumf1 + sfluxrefo(iprsm)
        sumf2 = sumf2 + absco2o(iprsm)*rwgt(iprsm+208)
        sumf3 = sumf3 + absh2oo(iprsm)*rwgt(iprsm+208)
      ENDDO
      sfluxref(igc) = sumf1
      absco2(igc) = sumf2
      absh2o(igc) = sumf3
    ENDDO

  END SUBROUTINE cmbgb29

END MODULE mo_psrad_srtm_setup

