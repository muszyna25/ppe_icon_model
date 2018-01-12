
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_psrad_lrtm_setup

  USE mo_psrad_general, ONLY : wp, mg, nbndlw, ngptlw, maxperband, &
    ih2o, ico2, io3, in2o, io2, ich4, ico, in2, &
    cfc_offset, iccl4, icfc11, icfc12, icfc22, &
    ih2oco2, ih2oo3, ih2on2o, ih2och4, &
    in2oco2, io3co2
  USE mo_psrad_lrtm_netcdf, ONLY : lrtm_read
  USE mo_psrad_lrtm_kgs

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: ngb, ngs, ngc, delwave, setup_lrtm, wavenum1, wavenum2
  !
  ! spectra information that is entered at run time
  !
  REAL(wp) :: rwgt(nbndlw*mg) !< Weights for combining original gpts into reduced gpts

  INTEGER, PARAMETER :: &
    ngb(ngptlw) = (/& ! The band index for each new g-interval
      1,1,1,1,1,1,1,1,1,1, &                       ! band 1
      2,2,2,2,2,2,2,2,2,2,2,2, &                   ! band 2
      3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3, &           ! band 3
      4,4,4,4,4,4,4,4,4,4,4,4,4,4, &               ! band 4
      5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5, &           ! band 5
      6,6,6,6,6,6,6,6, &                           ! band 6
      7,7,7,7,7,7,7,7,7,7,7,7, &                   ! band 7
      8,8,8,8,8,8,8,8, &                           ! band 8
      9,9,9,9,9,9,9,9,9,9,9,9, &                   ! band 9
      10,10,10,10,10,10, &                         ! band 10
      11,11,11,11,11,11,11,11, &                   ! band 11
      12,12,12,12,12,12,12,12, &                   ! band 12
      13,13,13,13, &                               ! band 13
      14,14, &                                     ! band 14
      15,15, &                                     ! band 15
      16,16/), &                                   ! band 16
    ngs(nbndlw) = (/& ! The cumulative sum of new g-intervals for each band
      10,22,38,52,68,76,88,96,108,114,122,130,134,136,138,140/), &
    ngc(nbndlw) = (/& ! The number of new g-intervals in each band
      10,12,16,14,16,8,12,8,12,6,8,8,4,2,2,2/)
  REAL(wp), PARAMETER :: &
    delwave(nbndlw) = (/& ! Spectral band width in wavenumbers
      340._wp, 150._wp, 130._wp,  70._wp, 120._wp, 160._wp, &
      100._wp, 100._wp, 210._wp,  90._wp, 320._wp, 280._wp, &
      170._wp, 130._wp, 220._wp, 650._wp/)

  !PRIVATE:
  INTEGER, PARAMETER :: &
    ngm(nbndlw*mg) = (/& ! The index of each new gpt relative to the orignal
      1,2,3,3,4,4,5,5,6,6,7,7,8,8,9,10, &          ! band 1
      1,2,3,4,5,6,7,8,9,9,10,10,11,11,12,12, &     ! band 2
      1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, &    ! band 3
      1,2,3,4,5,6,7,8,9,10,11,12,13,14,14,14, &    ! band 4
      1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, &    ! band 5
      1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8, &           ! band 6
      1,1,2,2,3,4,5,6,7,8,9,10,11,11,12,12, &      ! band 7
      1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8, &           ! band 8
      1,2,3,4,5,6,7,8,9,9,10,10,11,11,12,12, &     ! band 9
      1,1,2,2,3,3,4,4,5,5,5,5,6,6,6,6, &           ! band 10
      1,2,3,3,4,4,5,5,6,6,7,7,7,8,8,8, &           ! band 11
      1,2,3,4,5,5,6,6,7,7,7,7,8,8,8,8, &           ! band 12
      1,1,1,2,2,2,3,3,3,3,4,4,4,4,4,4, &           ! band 13
      1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2, &           ! band 14
      1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2, &           ! band 15
      1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2/), &         ! band 16
    ngn(ngptlw) = (/& ! The number of original gs combined to make new pts
      1,1,2,2,2,2,2,2,1,1, &                       ! band 1
      1,1,1,1,1,1,1,1,2,2,2,2, &                   ! band 2
      1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &           ! band 3
      1,1,1,1,1,1,1,1,1,1,1,1,1,3, &               ! band 4
      1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &           ! band 5
      2,2,2,2,2,2,2,2, &                           ! band 6
      2,2,1,1,1,1,1,1,1,1,2,2, &                   ! band 7
      2,2,2,2,2,2,2,2, &                           ! band 8
      1,1,1,1,1,1,1,1,2,2,2,2, &                   ! band 9
      2,2,2,2,4,4, &                               ! band 10
      1,1,2,2,2,2,3,3, &                           ! band 11
      1,1,1,1,2,2,4,4, &                           ! band 12
      3,3,4,6, &                                   ! band 13
      8,8, &                                       ! band 14
      8,8, &                                       ! band 15
      4,12/)                                    ! band 16

  ! RRTM weights for the original 16 g-intervals
  REAL(wp), PARAMETER :: wt(mg) = (/& 
      0.1527534276_wp, 0.1491729617_wp, 0.1420961469_wp, &
      0.1316886544_wp, 0.1181945205_wp, 0.1019300893_wp, &
      0.0832767040_wp, 0.0626720116_wp, 0.0424925000_wp, &
      0.0046269894_wp, 0.0038279891_wp, 0.0030260086_wp, &
      0.0022199750_wp, 0.0014140010_wp, 0.0005330000_wp, &
      0.0000750000_wp/)
  REAL(wp), PARAMETER :: wavenum1(nbndlw) = (/ & !< Spectral band lower boundary in wavenumbers
    10._wp, 350._wp, 500._wp, 630._wp, 700._wp, 820._wp, &
    980._wp,1080._wp,1180._wp,1390._wp,1480._wp,1800._wp, &
    2080._wp,2250._wp,2380._wp,2600._wp/)
  REAL(wp), PARAMETER :: wavenum2(nbndlw) = (/ & !< Spectral band upper boundary in wavenumbers
    350._wp, 500._wp, 630._wp, 700._wp, 820._wp, 980._wp, &
    1080._wp,1180._wp,1390._wp,1480._wp,1800._wp,2080._wp, &
    2250._wp,2380._wp,2600._wp,3250._wp/)

CONTAINS

  ! Compute relative weighting for new g-point combinations.
  SUBROUTINE init_rwgt
    INTEGER :: ibnd, igc, ig, ind, ipr, igcsm, iprsm
    REAL(wp) :: wtsum, wtsm(mg)
    igcsm = 0
    DO ibnd = 1,nbndlw
      iprsm = 0
      IF (ngc(ibnd).LT.mg) THEN
        DO igc = 1,ngc(ibnd) 
          igcsm = igcsm + 1
          wtsum = 0._wp
          DO ipr = 1, ngn(igcsm)
            iprsm = iprsm + 1
            wtsum = wtsum + wt(iprsm)
          ENDDO
          wtsm(igc) = wtsum
        ENDDO
        DO ig = 1, no
          ind = (ibnd-1)*mg + ig
          rwgt(ind) = wt(ig)/wtsm(ngm(ind))
        ENDDO
      ELSE
        DO ig = 1, no
          igcsm = igcsm + 1
          ind = (ibnd-1)*mg + ig
          rwgt(ind) = 1.0_wp
        ENDDO
      ENDIF
    ENDDO
  END SUBROUTINE init_rwgt

  SUBROUTINE allocate_all
    INTEGER :: band, atm, igas, M, N
    DO band = 1,nbndlw
    DO atm = 1,2
      IF (planck_fraction_interpolation_layer(atm,band) == 0) THEN
        nullify(planck_fraction2(atm,band)%v, &
          planck_fraction2_delta(atm,band)%v, &
          planck_fraction2_in(atm,band)%v)
      END IF
      nullify(kmajor(atm,band)%v, kmajor3_in(atm,band)%v, &
        kmajor4_in(atm,band)%v)
      DO igas = 1, ngas
        nullify(kgas2(igas,atm,band)%v)
        nullify(kgas2_delta(igas,atm,band)%v)
        nullify(kgas2_in(igas,atm,band)%v)
        nullify(kgas3(igas,atm,band)%v)
        nullify(kgas3_delta(igas,atm,band)%v)
        nullify(kgas3_in(igas,atm,band)%v)
      ENDDO
    ENDDO
    ENDDO
    DO band = 1,nbndlw
      DO igas = 1,ncfc
        IF (cfc_list(igas,band) /= 0) THEN
          allocate(cfc(igas,band)%v(maxperband), &
            cfc_in(igas,band)%v(maxperband))
        ELSE
          nullify(cfc(igas,band)%v, cfc_in(igas,band)%v)
        END IF
      END DO
      DO atm = 1,2
        IF (planck_fraction_interpolation_layer(atm,band) /= 0) THEN
          allocate(planck_fraction2(atm,band)%v(9,maxperband), &
            planck_fraction2_delta(atm,band)%v(9,maxperband), &
            planck_fraction2_in(atm,band)%v(9,maxperband))
        END IF
        n = nsp_species(atm,band)
        IF (n /= 0) THEN
          IF (n == 1) THEN
            allocate(kmajor3_in(atm,band)%v(5,npressure(atm),maxperband))
          ELSE
            allocate(kmajor4_in(atm,band)%v(n,5,npressure(atm),maxperband))
          ENDIF
          allocate(kmajor(atm,band)%v(5 * n * npressure(atm),maxperband))
        ENDIF
        DO igas = 1, ngas
          M = kgas2_list(igas,atm,band)
          IF (M /= 0) THEN
            allocate(kgas2(igas,atm,band)%v(M,maxperband))
            allocate(kgas2_delta(igas,atm,band)%v(M,maxperband))
            allocate(kgas2_in(igas,atm,band)%v(M,maxperband))
          ENDIF
          M = kgas3_list(1,igas,atm,band)
          IF (M /= 0) THEN
            N = kgas3_list(2,igas,atm,band)
            allocate(kgas3(igas,atm,band)%v(M,N,maxperband))
            allocate(kgas3_delta(igas,atm,band)%v(M,N,maxperband))
            allocate(kgas3_in(igas,atm,band)%v(M,N,maxperband))
          ENDIF
        ENDDO
      ENDDO
    ENDDO
  END SUBROUTINE allocate_all

  SUBROUTINE reduce_all
    INTEGER :: band, igas, m, n, which
    DO band = 1,nbndlw
      DO igas = 1,ncfc
        if (cfc_list(igas,band) /= 0) THEN
          if (cfc_list(igas,band) == 1) THEN
            CALL reduce0(band, cfc_in(igas,band)%v, cfc(igas,band)%v)
          ELSE
            CALL reduce0_no_wt(band, cfc_in(igas,band)%v, cfc(igas,band)%v)
          END IF
        END IF
      END DO
      DO which = 1,2 ! self/foreign
        CALL reduce1(nh2oref(which), band, &
          h2oref_in(:,:,which,band), h2oref(:,:,which,band))
      ENDDO
      DO which = 1,2 ! lower/upper
        n = nsp_species(which,band)
        IF (n /= 0) THEN
          IF (n == 1) THEN
            CALL reduce2_pack(5,npressure(which),band, &
              kmajor3_in(which,band)%v, kmajor(which,band)%v)
          ELSE
            CALL reduce3_pack(n,5,npressure(which),band, &
              kmajor4_in(which,band)%v, kmajor(which,band)%v)
          ENDIF
        ENDIF
        n = nsp_fraction(which,band)
        IF (n /= 0) THEN
          IF (n == 1) THEN
            CALL reduce0_no_wt(band, planck_fraction1_in(:,which,band), &
              planck_fraction1(:,which,band))
          ELSE
            CALL reduce1_no_wt(n,band, &
                planck_fraction2_in(which,band)%v, &
                planck_fraction2(which,band)%v)
            planck_fraction2_delta(which,band)%v = 0
            planck_fraction2_delta(which,band)%v(1:n-1,:) = &
              planck_fraction2(which,band)%v(2:n,:) - &
              planck_fraction2(which,band)%v(1:n-1,:)
          ENDIF
        ENDIF
        DO igas = 1, ngas
          m = kgas2_list(igas,which,band)
          IF (m /= 0) THEN
            CALL reduce1(m, band, kgas2_in(igas,which,band)%v, &
              kgas2(igas,which,band)%v)
            kgas2_delta(igas,which,band)%v = 0
            kgas2_delta(igas,which,band)%v(1:m-1,:) = &
              kgas2(igas,which,band)%v(2:m,:) - &
              kgas2(igas,which,band)%v(1:m-1,:)
          ENDIF
          m = kgas3_list(1,igas,which,band)
          IF (m /= 0) THEN
            n = kgas3_list(2,igas,which,band)
            CALL reduce2(m, n, band, kgas3_in(igas,which,band)%v, &
              kgas3(igas,which,band)%v)
            kgas3_delta(igas,which,band)%v = 0
            kgas3_delta(igas,which,band)%v(1:m-1,:,:) = &
              kgas3(igas,which,band)%v(2:m,:,:) - &
              kgas3(igas,which,band)%v(1:m-1,:,:)
          ENDIF
        ENDDO
      ENDDO
    ENDDO
    DO which = 1,2 ! self/foreign
      h2oref_delta(1:nh2oref(which)-1,:,which,:) = &
        h2oref(2:nh2oref(which),:,which,:) - &
        h2oref(1:nh2oref(which)-1,:,which,:)
    ENDDO
  END SUBROUTINE reduce_all

  ! Calculate reference ratio to be used in calculation of Planck
  ! fraction in lower/upper atmosphere.
  SUBROUTINE fill_tables
    INTEGER :: i, j, k
    DO k = 1,59
    DO j = 1,ngas
    DO i = 1,ngas
    !from taumol03:
    !refrat_planck_a = chi_mls(ih2o,9)/chi_mls(ico2,9) !  P = 212.725 mb
    !refrat_planck_b = chi_mls(ih2o,13)/chi_mls(ico2,13) !  P = 95.58 mb
    !refrat_m_a = chi_mls(ih2o,3)/chi_mls(ico2,3) !  P = 706.270mb
    !refrat_m_b = chi_mls(ih2o,13)/chi_mls(ico2,13) !  P = 95.58 mb 
    !from taumol04:
    !refrat_planck_a = chi_mls(ih2o,11)/chi_mls(ico2,11) ! P =   142.5940 mb
    !refrat_planck_b = chi_mls(io3,13)/chi_mls(ico2,13) ! P = 95.58350 mb
    !from taumol05:
    !refrat_planck_a = chi_mls(ih2o,5)/chi_mls(ico2,5) ! P = 473.420 mb
    !refrat_planck_b = chi_mls(io3,43)/chi_mls(ico2,43) ! P = 0.2369 mb
    !refrat_m_a = chi_mls(ih2o,7)/chi_mls(ico2,7) ! P = 317.3480
    !from taumol07:
    !refrat_planck_a = chi_mls(ih2o,3)/chi_mls(io3,3) ! P = 706.2620 mb
    !refrat_m_a = chi_mls(ih2o,3)/chi_mls(io3,3) ! P = 706.2720 mb
    !from taumol09:
    !refrat_planck_a = chi_mls(ih2o,9)/chi_mls(ich4,9) ! P = 212 mb
    !refrat_m_a = chi_mls(ih2o,3)/chi_mls(ich4,3) ! P = 706.272 mb 
    !from taumol12:
    !refrat_planck_a = chi_mls(ih2o,10)/chi_mls(ico2,10) ! P =   174.164 mb 
    !from taumol13:
    !refrat_planck_a = chi_mls(ih2o,5)/chi_mls(in2o,5) ! P = 473.420 mb (Level 5)
    !refrat_m_a = chi_mls(ih2o,1)/chi_mls(in2o,1) ! P = 1053. (Level 1)
    !refrat_m_a3 = chi_mls(ih2o,3)/chi_mls(in2o,3) ! P = 706. (Level 3)
    !from taumol15:
    !refrat_planck_a = chi_mls(in2o,1)/chi_mls(ico2,1) ! P = 1053. mb (Level 1)
    !refrat_m_a = chi_mls(in2o,1)/chi_mls(ico2,1) ! P = 1053.
    !from taumol16:
    !refrat_planck_a = chi_mls(ih2o,6)/chi_mls(ich4,6) ! P = 387. mb (Level 6)

      IF (chi_mls(j,k) /= 0) THEN
        planck_ratio(i,j,k) = chi_mls(i,k) / chi_mls(j,k) 
      ELSE
        planck_ratio(i,j,k) = 0
      ENDIF
    END DO
    END DO
    END DO

    reaction_table = 0
    reaction_table(ih2o,ico2) = ih2oco2
    reaction_table(ih2o,io3) = ih2oo3
    reaction_table(ih2o,in2o) = ih2on2o
    reaction_table(ih2o,ich4) = ih2och4
    reaction_table(in2o,ico2) = in2oco2
    reaction_table(io3,ico2) = io3co2

  END SUBROUTINE fill_tables

  SUBROUTINE fill_in_missing
    INTEGER :: band, present, absent
    DO band = 1,nbndlw
      absent = missing_planck_fraction(band)
      IF (absent /= 0) THEN
        present = 3 - absent
        IF (planck_fraction_interpolation_layer(present,band) == 0) THEN
          planck_fraction1(:,absent,band) = planck_fraction1(:,present,band)
        ELSE
          planck_fraction2(absent,band)%v = planck_fraction2(present,band)%v
        END IF
      END IF
    END DO
  END SUBROUTINE fill_in_missing

  SUBROUTINE setup_tau_info
    INTEGER :: band, atm
    n_key_species = 0
    n_minor_species = 0
    stratosphere_fudge_flag = 0
    DO band = 1, nbndlw
      DO atm = 1, 2
        n_key_species(atm,band) = COUNT(key_species(:,atm,band) /= 0)
        n_minor_species(atm,band) = COUNT(minor_species(:,atm,band) /= 0)
      END DO
      IF (ANY(stratosphere_fudge(:,band) /= 0)) THEN
        stratosphere_fudge_flag(band) = 1
      ENDIF
    END DO
  END SUBROUTINE setup_tau_info

  SUBROUTINE setup_lrtm

    kgas2_list = 0
    kgas2_list(in2,1,1) = 19
    kgas2_list(in2,2,1) = 19
    kgas2_list(ico2,1,6) = 19
    kgas2_list(ico2,2,7) = 19
    kgas2_list(ico2,1,8) = 19
    kgas2_list(in2o,1,8) = 19
    kgas2_list(io3,1,8) = 19
    kgas2_list(ico2,2,8) = 19
    kgas2_list(in2o,2,8) = 19
    kgas2_list(in2o,2,8) = 19
    kgas2_list(in2o,2,9) = 19
    kgas2_list(io2,1,11) = 19
    kgas2_list(io2,2,11) = 19
    kgas2_list(io3,2,13) = 19

    kgas3_list = 0
    kgas3_list(:,in2o,1,3) = (/9,19/)
    kgas3_list(:,in2o,2,3) = (/5,19/)
    kgas3_list(:,io3,1,5) = (/9,19/)
    kgas3_list(:,ico2,1,7) = (/9,19/)
    kgas3_list(:,in2o,1,9) = (/9,19/)
    kgas3_list(:,ico2,1,13) = (/9,19/)
    kgas3_list(:,ico,1,13) = (/9,19/)
    kgas3_list(:,in2,1,15) = (/9,19/)

    cfc_list = 0
    cfc_list(iccl4-cfc_offset,5) = 1
    cfc_list(icfc11-cfc_offset,6) = 1
    cfc_list(icfc12-cfc_offset,6) = 1
    cfc_list(icfc12-cfc_offset,8) = 1
    cfc_list(icfc22-cfc_offset,8) = 1


    skip_atmosphere = 0
    skip_atmosphere(2,12) = 1
    skip_atmosphere(2,15) = 1

    missing_planck_fraction = 0
    missing_planck_fraction(6) = 2

    key_species = 0
    minor_species = 0
    minor_species_scale = 0
    minor_species_fudge = 0 !(4,max_minor_species,2,nbndlw)
    minor_species_interpolation_layer = 0
    planck_fraction_interpolation_layer = 0
    stratosphere_fudge = 0
    ! band 1:  10-350 cm-1 (low key - h2o; low minor - n2)
    !                      (high key - h2o; high minor - n2)
    !   old:   10-250 cm-1 (low - h2o; high - h2o)
    ! Minor gas mapping levels:
    !     lower - n2, p = 142.5490 mbar, t = 215.70 k
    !     upper - n2, p = 142.5490 mbar, t = 215.70 k
    key_species(1,1,1) = ih2o
    minor_species(1,1,1) = in2
    minor_species_scale(1,1,1) = 1
    key_species(1,2,1) = ih2o
    minor_species(1,2,1) = in2
    minor_species_scale(1,2,1) = 1
    ! band 2:  350-500 cm-1 (low key - h2o; high key - h2o)
    !   old:   250-500 cm-1 (low - h2o; high - h2o)
    key_species(1,1,2) = ih2o
    key_species(1,2,2) = ih2o
    ! band 3:  500-630 cm-1 (low key - h2o,co2; low minor - n2o)
    !                       (high key - h2o,co2; high minor - n2o)
    !    old:  500-630 cm-1 (low - h2o,co2; high - h2o,co2)
    ! Minor gas mapping levels:
    !     lower - n2o, p = 706.272 mbar, t = 278.94 k
    !     upper - n2o, p = 95.58 mbar, t = 215.7 k
    key_species(:,1,3) = (/ih2o, ico2/)
    minor_species(1,1,3) = in2o
    minor_species_fudge(:,1,1,3) = (/1.5_wp, 0.5_wp, 0.65_wp, 0.0_wp/)
    minor_species_interpolation_layer(1,1,3) = 3
    planck_fraction_interpolation_layer(1,3) = 9
    key_species(:,2,3) = (/ih2o, ico2/)
    minor_species(1,2,3) = in2o
    minor_species_fudge(:,1,2,3) = (/1.5_wp, 0.5_wp, 0.65_wp, 0.0_wp/)
    minor_species_interpolation_layer(1,2,3) = 13
    planck_fraction_interpolation_layer(2,3) = 13
    ! band 4:  630-700 cm-1 (low key - h2o,co2; high key - o3,co2)
    !    old:  630-700 cm-1 (low - h2o,co2; high - o3,co2)
    key_species(:,1,4) = (/ih2o, ico2/)
    planck_fraction_interpolation_layer(1,4) = 11
    key_species(:,2,4) = (/io3, ico2/)
    planck_fraction_interpolation_layer(2,4) = 13
    stratosphere_fudge(1:ngc(4),4) = &
      (/ 1., 1., 1., 1., 1., 1., 1., .92, .88, 1.07, 1.1, .99, .88, .943 /) 
    ! band 5:  700-820 cm-1 (low key - h2o,co2; low minor - o3, ccl4)
    !                       (high key - o3,co2)
    !    old:  700-820 cm-1 (low - h2o,co2; high - o3,co2)
    ! Minor gas mapping level :
    !     lower - o3, p = 317.34 mbar, t = 240.77 k
    !     lower - ccl4
    key_species(:,1,5) = (/ih2o,ico2/)
    minor_species(1:2,1,5) = (/io3,iccl4/)
    minor_species_interpolation_layer(1,1,5) = 7
    planck_fraction_interpolation_layer(1,5) = 5
    key_species(:,2,5) = (/io3,ico2/)
    planck_fraction_interpolation_layer(2,5) = 43
    !band 6:  820-980 cm-1 (low key - h2o; low minor - co2)
    !                      (high key - nothing; high minor - cfc11, cfc12)
    !   old:  820-980 cm-1 (low - h2o; high - nothing)
    ! NOTE: code has cfc11/12 in lower atmosphere as well.
    ! Minor gas mapping level:
    !     lower - co2, p = 706.2720 mb, t = 294.2 k
    !     upper - cfc11, cfc12
    key_species(1,1,6) = ih2o
    minor_species(1:3,1,6) = (/ico2,icfc11,icfc12/)
    minor_species_fudge(:,1,1,6) = (/3.0_wp, 2.0_wp, 0.77_wp, 0.0_wp/)
    minor_species(1:2,2,6) = (/icfc11, icfc12/)
    ! band 7:  98051080 cm-1 (low key - h2o,o3; low minor - co2)
    !                        (high key - o3; high minor - co2)
    !    old:  980-1080 cm-1 (low - h2o,o3; high - o3)
    ! Minor gas mapping level :
    !     lower - co2, p = 706.2620 mbar, t= 278.94 k
    !     upper - co2, p = 12.9350 mbar, t = 234.01 k
    key_species(:,1,7) = (/ih2o, io3/)
    minor_species(1,1,7) = ico2
    minor_species_fudge(:,1,1,7) = (/3.0_wp, 3.0_wp, 0.79_wp, 0.0_wp/)
    minor_species_interpolation_layer(1,1,7) = 3
    planck_fraction_interpolation_layer(1,7) = 3
    key_species(1,2,7) = io3
    minor_species(1,2,7) = ico2
    minor_species_fudge(:,1,2,7) = (/3.0_wp, 2.0_wp, 0.79_wp, 0.0_wp/)
    stratosphere_fudge(1:ngc(7),7) = &
      (/ 1., 1., 1., 1., 1., .92, .88, 1.07, 1.1, .99, .855, 1. /) 
    ! band 8:  1080-1180 cm-1 (low key - h2o; low minor - co2,o3,n2o)
    !                         (high key - o3; high minor - co2, n2o)
    ! NOTE: Code contained cfc12 and cfc22 as well!!!
    !    old:  1080-1180 cm-1 (low (i.e.>~300mb) - h2o; high - o3)
    ! Minor gas mapping level:
    !     lower - co2, p = 1053.63 mb, t = 294.2 k
    !     lower - o3,  p = 317.348 mb, t = 240.77 k
    !     lower - n2o, p = 706.2720 mb, t= 278.94 k
    !     lower - cfc12,cfc11
    !     upper - co2, p = 35.1632 mb, t = 223.28 k
    !     upper - n2o, p = 8.716e-2 mb, t = 226.03 k
    key_species(1,1,8) = ih2o
    minor_species(1:5,1,8) = (/ico2, io3, in2o, icfc12, icfc22/)
    minor_species_fudge(:,1,1,8) = (/3.0_wp, 2.0_wp, 0.65_wp, 0.0_wp/)
    key_species(1,2,8) = io3
    minor_species(1:4,2,8) = (/ico2, in2o, icfc12, icfc22/)
    minor_species_fudge(:,1,2,8) = (/3.0_wp, 2.0_wp, 0.65_wp, 0.0_wp/)
    ! band 9:  1180-1390 cm-1 (low key - h2o,ch4; low minor - n2o)
    !                         (high key - ch4; high minor - n2o)
    !    old:  1180-1390 cm-1 (low - h2o,ch4; high - ch4)
    ! Minor gas mapping level :
    !     lower - n2o, p = 706.272 mbar, t = 278.94 k
    !     upper - n2o, p = 95.58 mbar, t = 215.7 k
    key_species(:,1,9) = (/ih2o, ich4/)
    minor_species(1,1,9) = in2o
    minor_species_fudge(:,1,1,9) = (/1.5_wp, 0.5_wp, 0.65_wp, 0.0_wp/)
    minor_species_interpolation_layer(1,1,9) = 3
    planck_fraction_interpolation_layer(1,9) = 9
    key_species(1,2,9) = ich4
    minor_species(1,2,9) = in2o
    minor_species_fudge(:,1,2,9) = (/1.5_wp, 0.5_wp, 0.65_wp, 0.0_wp/)
    ! band 10:  1390-1480 cm-1 (low key - h2o; high key - h2o)
    !     old:  1390-1480 cm-1 (low - h2o; high - h2o)
    key_species(1,1,10) = ih2o
    key_species(1,2,10) = ih2o
    ! band 11:  1480-1800 cm-1 (low - h2o; low minor - o2)
    !                          (high key - h2o; high minor - o2)
    !     old:  1480-1800 cm-1 (low - h2o; low minor - o2)
    !                          (high key - h2o; high minor - o2)
    ! Minor gas mapping level :
    !     lower - o2, p = 706.2720 mbar, t = 278.94 k
    !     upper - o2, p = 4.758820 mbarm t = 250.85 k
    key_species(1,1,11) = ih2o
    minor_species(1,1,11) = io2
    minor_species_scale(1,1,11) = 2
    key_species(1,2,11) = ih2o
    minor_species(1,2,11) = io2
    minor_species_scale(1,2,11) = 2
    ! band 12:  1800-2080 cm-1 (low - h2o,co2; high - nothing)
    !     old:  1800-2080 cm-1 (low - h2o,co2; high - nothing)
    key_species(:,1,12) = (/ih2o, ico2/)
    planck_fraction_interpolation_layer(1,12) = 10
    ! band 13:  2080-2250 cm-1 (low key - h2o,n2o; high minor - o3 minor)
    !     old:  2080-2250 cm-1 (low - h2o,n2o; high - nothing)
    ! NOTE: comment did not reflect code, see minor species for lower atm
    ! Minor gas mapping levels :
    !     lower - co2, p = 1053.63 mb, t = 294.2 k
    !     lower - co, p = 706 mb, t = 278.94 k
    !     upper - o3, p = 95.5835 mb, t = 215.7 k
    key_species(:,1,13) = (/ih2o,in2o/)
    minor_species(1:2,1,13) = (/ico2,ico/)
    minor_species_fudge(:,1,1,13) = (/3.0_wp, 2.0_wp, 0.68_wp, 3.55e-4_wp/)
    minor_species_interpolation_layer(1:2,1,13) = (/1,3/)
    planck_fraction_interpolation_layer(1,13) = 5
    minor_species(1,2,13) = io3
    ! band 14:  2250-2380 cm-1 (low - co2; high - co2)
    !     old:  2250-2380 cm-1 (low - co2; high - co2)
    key_species(1,1,14) = ico2
    key_species(1,2,14) = ico2
    ! band 15:  2380-2600 cm-1 (low - n2o,co2; low minor - n2)
    !                          (high - nothing)
    !     old:  2380-2600 cm-1 (low - n2o,co2; high - nothing)
    ! Minor gas mapping level : 
    !     Lower - Nitrogen Continuum, P = 1053., T = 294.
    key_species(:,1,15) = (/in2o, ico2/)
    minor_species(1,1,15) = in2
    minor_species_scale(1,1,15) = 3
    minor_species_interpolation_layer(1,1,15) = 1
    planck_fraction_interpolation_layer(1,15) = 1
    ! band 16:  2600-3250 cm-1 (low key- h2o,ch4; high key - ch4)
    !     old:  2600-3000 cm-1 (low - h2o,ch4; high - nothing)
    key_species(:,1,16) = (/ih2o, ich4/)
    planck_fraction_interpolation_layer(1,16) = 6
    key_species(1,2,16) = ich4

    h2o_absorption_flag = RESHAPE((/ &
      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, & ! self,lower
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, & ! self,upper
      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, & ! foreign, lower
      1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0/), & !foreign, upper
      SHAPE = (/2,2,nbndlw/), ORDER = (/3,2,1/))

    pressure_dependent_tau_correction = 0
    pressure_dependent_tau_correction(1,1) = 1
    pressure_dependent_tau_correction(2,1) = 2
    pressure_dependent_tau_correction(1,2) = 3


    CALL setup_tau_info
    CALL allocate_all
    CALL init_rwgt
    CALL lrtm_read
    CALL reduce_all
    CALL fill_in_missing
    CALL fill_tables

  END SUBROUTINE setup_lrtm

  !  The subroutines CMBGB1->CMBGB16 input the absorption coefficient
  !  data for each band, which are defined for 16 g-points and 16 spectral
  !  bands. The data are combined with appropriate weighting following the
  !  g-point mapping arrays specified in RRTMINIT.  Plank fraction data
  !  in arrays FRACREFA and FRACREFB are combined without weighting.  All
  !  g-point reduced data are put into new arrays for use in RRTM.
  SUBROUTINE reduce0(band, ref, k)
    INTEGER, INTENT(IN) :: band
    REAL(wp), INTENT(IN) :: ref(:)
    REAL(wp), INTENT(OUT) :: k(:)

    INTEGER :: ipr, iprsm, iprsoff, igc, ngsoff
    REAL(wp) :: sumk

    k = 0
    iprsoff = (band-1) * 16
    if (band == 1) then
      ngsoff = 0
    else
      ngsoff = ngs(band-1)
    end if
    iprsm = 0
    DO igc = 1,ngc(band)
      sumk = 0.0_wp
      DO ipr = 1, ngn(ngsoff + igc)
        iprsm = iprsm + 1
        sumk = sumk + ref(iprsm)*rwgt(iprsm+iprsoff)
      ENDDO
      k(igc) = sumk
    ENDDO
  END SUBROUTINE reduce0

  SUBROUTINE reduce0_no_wt(band, ref, k)
    INTEGER, INTENT(IN) :: band
    REAL(wp), INTENT(IN) :: ref(:)
    REAL(wp), INTENT(OUT) :: k(:)

    INTEGER :: ipr, iprsm, igc, ngsoff
    REAL(wp) :: sumk

    k = 0
    if (band == 1) then
      ngsoff = 0
    else
      ngsoff = ngs(band-1)
    end if
    iprsm = 0
    DO igc = 1,ngc(band)
      sumk = 0.0_wp
      DO ipr = 1, ngn(ngsoff + igc)
        iprsm = iprsm + 1
        sumk = sumk + ref(iprsm)
      ENDDO
      k(igc) = sumk
    ENDDO
  END SUBROUTINE reduce0_no_wt


  SUBROUTINE reduce1_no_wt(nt, band, ref, k)
    INTEGER, INTENT(IN) :: nt, band
    REAL(wp), INTENT(IN) :: ref(:,:) !NOTE: Was ref(no,nt), because of input
    REAL(wp), INTENT(INOUT) :: k(:,:) !fracref[ab] => max(5,9)!

    INTEGER :: jt
    DO jt = 1,nt
      CALL reduce0_no_wt(band, ref(jt,:), k(jt,:))
    ENDDO
  END SUBROUTINE reduce1_no_wt


  SUBROUTINE reduce1(nt, band, ref, k)
    INTEGER, INTENT(IN) :: nt, band
    REAL(wp), INTENT(IN) :: ref(:,:)
    REAL(wp), INTENT(OUT) :: k(:,:)

    INTEGER :: jt

    DO jt = 1,nt
      CALL reduce0(band, ref(jt,:), k(jt,:))
    ENDDO
  END SUBROUTINE reduce1

  SUBROUTINE reduce2(nt, np, band, ref, k)
    INTEGER, INTENT(IN) :: nt, np, band
    REAL(wp), INTENT(IN) :: ref(:,:,:)
    REAL(wp), INTENT(OUT) :: k(:,:,:)

    INTEGER :: jt, jp

    DO jt = 1,nt
      DO jp = 1,np
        CALL reduce0(band, ref(jt,jp,:), k(jt,jp,:))
      ENDDO
    ENDDO
  END SUBROUTINE reduce2

  SUBROUTINE reduce2_pack(nt, np, band, ref, k)
    INTEGER, INTENT(IN) :: nt, np, band
    REAL(wp), INTENT(IN) :: ref(:,:,:)
    REAL(wp), INTENT(OUT) :: k(:,:)
    REAL(wp) :: tmp(nt,np,size(ref,3))

    CALL reduce2(nt, np, band, ref, tmp)
    k = RESHAPE(tmp, SHAPE = (/nt*np,size(ref,3)/))
  END SUBROUTINE reduce2_pack

  SUBROUTINE reduce3(nn, nt, np, band, ref, k)
    INTEGER, INTENT(IN) :: nn, nt, np, band
    REAL(wp), INTENT(IN) :: ref(:,:,:,:)
    REAL(wp), INTENT(OUT) :: k(:,:,:,:)

    INTEGER :: jn, jt, jp

    DO jn = 1,nn
      DO jt = 1,nt
        DO jp = 1,np
          CALL reduce0(band, ref(jn,jt,jp,:), k(jn,jt,jp,:))
        ENDDO
      ENDDO
    ENDDO
  END SUBROUTINE reduce3

  SUBROUTINE reduce3_pack(nn, nt, np, band, ref, k)
    INTEGER, INTENT(IN) :: nn, nt, np, band
    REAL(wp), INTENT(IN) :: ref(:,:,:,:)
    REAL(wp), INTENT(OUT) :: k(:,:)
    REAL(wp) :: tmp(nn,nt,np,size(ref,4))

    CALL reduce3(nn, nt, np, band, ref, tmp)
    k = RESHAPE(tmp, SHAPE = (/nn*nt*np,size(ref,4)/))
  END SUBROUTINE reduce3_pack


END MODULE mo_psrad_lrtm_setup
