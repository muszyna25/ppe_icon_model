! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
!     path:      $Source: /storm/rc1/cvsroot/rc/rrtmg_lw/src/rrtmg_lw_taumol.f90,v $
!     author:    $Author: mike $
!     revision:  $Revision: 1.7 $
!     created:   $Date: 2009/10/20 15:08:37 $
!

! *     TAUMOL
! *
! *     This file contains the subroutines TAUGBn (where n goes from
! *     1 to 16).  TAUGBn calculates the optical depths and Planck fractions
! *     per g-value and layer for band n.
! *
! *  Output:  optical depths (unitless)
! *           fractions needed to compute Planck functions at every layer
! *               and g-value
! *
! *     COMMON /TAUGCOM/  TAUG(MXLAY,MG)
! *     COMMON /PLANKG/   FRACS(MXLAY,MG)
! *
! *  Input
! *
! *     COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)
! *     COMMON /PRECISE/  ONEMINUS
! *     COMMON /PROFILE/  NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
! *     &                 PZ(0:MXLAY),TZ(0:MXLAY)
! *     COMMON /PROFDATA/ LAYTROP,
! *    &                  COLH2O(MXLAY),COLCO2(MXLAY),COLO3(MXLAY),
! *    &                  COLN2O(MXLAY),COLCO(MXLAY),COLCH4(MXLAY),
! *    &                  COLO2(MXLAY)
! *     COMMON /INTFAC/   FAC00(MXLAY),FAC01(MXLAY),
! *    &                  FAC10(MXLAY),FAC11(MXLAY)
! *     COMMON /INTIND/   JP(MXLAY),JT(MXLAY),JT1(MXLAY)
! *     COMMON /SELF/     SELFFAC(MXLAY), SELFFRAC(MXLAY), INDSELF(MXLAY)
! *
! *     Description:
! *     NG(IBAND) - number of g-values in band IBAND
! *     NSPA(IBAND) - for the lower atmosphere, the number of reference
! *                   atmospheres that are stored for band IBAND per
! *                   pressure level and temperature.  Each of these
! *                   atmospheres has different relative amounts of the
! *                   key species for the band (i.e. different binary
! *                   species parameters).
! *     NSPB(IBAND) - same for upper atmosphere
! *     ONEMINUS - since problems are caused in some cases by interpolation
! *                parameters equal to or greater than 1, for these cases
! *                these parameters are set to this value, slightly < 1.
! *     PAVEL - layer pressures (mb)
! *     TAVEL - layer temperatures (degrees K)
! *     PZ - level pressures (mb)
! *     TZ - level temperatures (degrees K)
! *     LAYTROP - layer at which switch is made from one combination of
! *               key species to another
! *     COLH2O, COLCO2, COLO3, COLN2O, COLCH4 - column amounts of water
! *               vapor,carbon dioxide, ozone, nitrous ozide, methane,
! *               respectively (molecules/cm**2)
! *     FACij(LAY) - for layer LAY, these are factors that are needed to
! *                  compute the interpolation factors that multiply the
! *                  appropriate reference k-values.  A value of 0 (1) for
! *                  i,j indicates that the corresponding factor multiplies
! *                  reference k-value for the lower (higher) of the two
! *                  appropriate temperatures, and altitudes, respectively.
! *     JP - the index of the lower (in altitude) of the two appropriate
! *          reference pressure levels needed for interpolation
! *     JT, JT1 - the indices of the lower of the two appropriate reference
! *               temperatures needed for interpolation (for pressure
! *               levels JP and JP+1, respectively)
! *     SELFFAC - scale factor needed for water vapor self-continuum, equals
! *               (water vapor density)/(atmospheric density at 296K and
! *               1013 mb)
! *     SELFFRAC - factor needed for temperature interpolation of reference
! *                water vapor self-continuum data
! *     INDSELF - index of the lower of the two appropriate reference
! *               temperatures needed for the self-continuum interpolation
! *     FORFAC  - scale factor needed for water vapor foreign-continuum.
! *     FORFRAC - factor needed for temperature interpolation of reference
! *                water vapor foreign-continuum data
! *     INDFOR  - index of the lower of the two appropriate reference
! *               temperatures needed for the foreign-continuum interpolation
! *
! *  Data input
! *     COMMON /Kn/ KA(NSPA(n),5,13,MG), KB(NSPB(n),5,13:59,MG), SELFREF(10,MG),
! *                 FORREF(4,MG), KA_M'MGAS', KB_M'MGAS'
! *        (note:  n is the band number,'MGAS' is the species name of the minor
! *         gas)
! *
! *     Description:
! *     KA - k-values for low reference atmospheres (key-species only)
! *          (units: cm**2/molecule)
! *     KB - k-values for high reference atmospheres (key-species only)
! *          (units: cm**2/molecule)
! *     KA_M'MGAS' - k-values for low reference atmosphere minor species
! *          (units: cm**2/molecule)
! *     KB_M'MGAS' - k-values for high reference atmosphere minor species
! *          (units: cm**2/molecule)
! *     SELFREF - k-values for water vapor self-continuum for reference
! *               atmospheres (used below LAYTROP)
! *               (units: cm**2/molecule)
! *     FORREF  - k-values for water vapor foreign-continuum for reference
! *               atmospheres (used below/above LAYTROP)
! *               (units: cm**2/molecule)
! *
! *     DIMENSION ABSA(65*NSPA(n),MG), ABSB(235*NSPB(n),MG)
! *     EQUIVALENCE (KA,ABSA),(KB,ABSB)
#include "mod1.inc"
MODULE mo_psrad_lrtm_gas_optics

  USE mo_psrad_general, ONLY : wp, nbndlw, ngptlw, &
    ngas, ncfc, cfc_offset, nreact, maxperband
  USE mo_psrad_lrtm_setup, ONLY : ngs, ngc !ngb
  USE mo_psrad_lrtm_kgs, ONLY : chi_mls, &
    h2oref, h2oref_delta, kgas2_list, kgas2, kgas2_delta, &
    kgas3_list, kgas3, kgas3_delta, &
    kmajor, nsp_species, nsp_fraction, cfc, &
    nsp_species_broken_16, &
    planck_fraction1, planck_fraction2, planck_fraction2_delta, &
    planck_ratio, planck_fraction_interpolation_layer, &
    skip_atmosphere, &
    key_species, reaction_table, h2o_absorption_flag, &
    minor_species, minor_species_scale, n_minor_species, &
    minor_species_fudge, &
    minor_species_interpolation_layer, &
    pressure_dependent_tau_correction, &
    stratosphere_fudge, stratosphere_fudge_flag

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: gas_optics_lw

  REAL(wp), PARAMETER :: oneminus = 1.0_wp - 1.0e-06_wp

  INTEGER, PARAMETER :: fracs_mult(2,nbndlw) = MAX(0,nsp_fraction-1)

CONTAINS

  SUBROUTINE gas_optics_lw(jcs, kproma, kbdim, klev, play, &
    wx, coldry, laytrop, jp1, iabs, gases, colbrd, &
    fac, reaction_ratio, &
    h2o_factor, h2o_fraction, h2o_index, &
    minorfrac, scaleminor, scaleminorn2, indminor, fracs_ret, tau_ret)

    INTEGER, INTENT(in) :: jcs, kproma, kbdim, klev 
    INTEGER, INTENT(in) :: laytrop(KBDIM), & ! tropopause layer index
      iabs(KBDIM,2,2,klev)
    REAL(wp), DIMENSION(KBDIM,klev,2), INTENT(IN) :: h2o_factor,h2o_fraction
    INTEGER, DIMENSION(KBDIM,klev,2), INTENT(IN) :: h2o_index
    INTEGER, DIMENSION(KBDIM,klev), INTENT(in) :: &
      jp1, indminor
    REAL(wp), INTENT(IN) :: &
      wx(KBDIM,ncfc,klev) ! cross-section amounts (mol/cm2)
    REAL(wp), DIMENSION(KBDIM,klev), INTENT(IN) :: &
      play, & ! (klev) layer pressures [mb]
      coldry, & ! (klev) column amount (dry air)
      colbrd, &
      minorfrac, &
      scaleminor, scaleminorn2
    REAL(wp), INTENT(IN) :: fac(KBDIM,2,2,klev), gases(KBDIM,klev,ngas), &
      reaction_ratio(KBDIM,2,klev,nreact)
    REAL(wp) :: tau_upper(KBDIM,klev,ngptlw), fracs_upper(KBDIM,klev,ngptlw)

    ! Output arrays have size (klev)
    REAL(wp), DIMENSION(KBDIM,klev,ngptlw), INTENT(OUT) :: &
      fracs_ret, & ! planck fractions
      tau_ret ! gaseous optical depth 

    INTEGER :: gpt, lay
    INTEGER :: atm_range(2,2), merge_range(2), range2(2)
    LOGICAL :: must_merge

    merge_range = (/minval(laytrop(jcs:kproma)), maxval(laytrop(jcs:kproma))/)
    must_merge = merge_range(1) < merge_range(2)
    IF (must_merge) THEN
      atm_range = RESHAPE((/&
        1, merge_range(2), &
        merge_range(2)+1, klev/), SHAPE=(/2,2/))
      merge_range(1) = merge_range(1)+1
    ELSE
      atm_range = RESHAPE((/&
        1, merge_range(1), &
        merge_range(1)+1, klev/), SHAPE=(/2,2/))
    END IF

    !CALL do_atmosphere(1, atm_range(:,1), tau_ret, fracs_ret)
    !CALL do_atmosphere(2, atm_range(:,2), tau_ret, fracs_ret)

    range2(1) = atm_range(1,1)
    range2(2) = atm_range(2,1)
    CALL do_atmosphere(jcs, kproma, kbdim, klev, 1, range2, fac, gases, reaction_ratio, &
         &             iabs, h2o_factor, h2o_fraction, h2o_index, jp1, indminor,                &
         &             wx, play, coldry, colbrd, minorfrac, scaleminor, scaleminorn2,           &
         &             tau_ret, fracs_ret)
    !CALL do_atmosphere(jcs, kproma, kbdim, klev, 2, range2, fac, gases, reaction_ratio, &
         !&             iabs, h2o_factor, h2o_fraction, h2o_index, jp1, indminor,                &
         !&             wx, play, coldry, colbrd, minorfrac, scaleminor, scaleminorn2,           &
         !&             tau_ret, fracs_ret)

    range2(1) = atm_range(1,2)
    range2(2) = atm_range(2,2)
    !CALL do_atmosphere(jcs, kproma, kbdim, klev, 1, range2, fac, gases, reaction_ratio, &
         !&             iabs, h2o_factor, h2o_fraction, h2o_index, jp1, indminor,                &
         !&             wx, play, coldry, colbrd, minorfrac, scaleminor, scaleminorn2,           &
         !&             tau_ret, fracs_ret)
    CALL do_atmosphere(jcs, kproma, kbdim, klev, 2, range2, fac, gases, reaction_ratio, &
         &             iabs, h2o_factor, h2o_fraction, h2o_index, jp1, indminor,                &
         &             wx, play, coldry, colbrd, minorfrac, scaleminor, scaleminorn2,           &
         &             tau_ret, fracs_ret)


    IF (must_merge) THEN
      !CALL do_atmosphere(2, merge_range, tau_upper, fracs_upper)
      CALL do_atmosphere(jcs, kproma, kbdim, klev, 2, merge_range, fac, gases, reaction_ratio, &
         &               iabs, h2o_factor, h2o_fraction, h2o_index, jp1, indminor,             &
         &               wx, play, coldry, colbrd, minorfrac, scaleminor, scaleminorn2,        &
         &               tau_upper, fracs_upper)
      DO gpt = 1,ngptlw
      DO lay = merge_range(1),merge_range(2)
        WHERE (laytrop(jcs:kproma) < lay)
          tau_ret(jcs:kproma,lay,gpt) = tau_upper(jcs:kproma,lay,gpt)
          fracs_ret(jcs:kproma,lay,gpt) = fracs_upper(jcs:kproma,lay,gpt)
        END WHERE
      END DO
      END DO
    END IF

  END SUBROUTINE gas_optics_lw
    !CONTAINS

    !SUBROUTINE do_atmosphere(atm, range_in, tau, fracs)
    SUBROUTINE do_atmosphere(jcs, kproma, kbdim, klev, atm, zrange_in, &
                             fac, gases, reaction_ratio, iabs, h2o_factor, h2o_fraction, h2o_index, &
                             jp1, indminor, &
                             wx, play, coldry, colbrd, minorfrac, scaleminor, scaleminorn2, &
                             tau, fracs)

      INTEGER, INTENT(IN) :: jcs, kproma, kbdim, klev
      INTEGER, INTENT(IN) :: atm, zrange_in(2)

      REAL(wp), INTENT(IN) :: fac(KBDIM,2,2,klev), gases(KBDIM,klev,ngas), &
                              reaction_ratio(KBDIM,2,klev,nreact)
      REAL(wp), DIMENSION(KBDIM,klev,2), INTENT(IN) :: h2o_factor, h2o_fraction
      INTEGER, INTENT(in) :: iabs(KBDIM,2,2,klev)
      INTEGER, DIMENSION(KBDIM,klev,2), INTENT(IN) :: h2o_index
      INTEGER, DIMENSION(KBDIM,klev), INTENT(in) :: jp1, indminor
      REAL(wp), INTENT(IN) :: wx(KBDIM,ncfc,klev) ! cross-section amounts (mol/cm2)
      REAL(wp), DIMENSION(KBDIM,klev), INTENT(IN) :: &
                play, & ! (klev) layer pressures [mb]
                coldry, & ! (klev) column amount (dry air)
                colbrd, &
                minorfrac, &
                scaleminor, scaleminorn2

      REAL(wp), DIMENSION(KBDIM,klev,ngptlw), INTENT(INOUT) :: fracs, tau

      REAL(wp) :: minorscale(KBDIM,klev), ratio
      INTEGER :: gpt_low, gpt_high, gpt_in, &
                 band, igas1, igas2, ireaction, gpt, which, igas_minor, zrange(2)

      zrange(1) = min(max(zrange_in(1),1),klev)
      zrange(2) = min(max(zrange_in(2),1),klev)

      gpt_high = 0
      DO band = 1,nbndlw
        gpt_low = gpt_high+1
        gpt_high = ngs(band)
        !n_in_band = ngc(band)
        IF (skip_atmosphere(atm,band) /= 0) THEN
          tau(:,zrange(1):zrange(2), gpt_low:gpt_high) = 0.
          fracs(:,zrange(1):zrange(2), gpt_low:gpt_high) = 0.
        ELSE
          IF (planck_fraction_interpolation_layer(atm,band) == 0) THEN
            gpt_in = 0
            DO gpt = gpt_low,gpt_high
              gpt_in = gpt_in+1
              fracs(:,zrange(1):zrange(2), gpt) = &
                planck_fraction1(gpt_in,atm,band)
            END DO
          ELSE
            igas1 = key_species(1,atm,band)
            igas2 = key_species(2,atm,band)
            ratio = planck_ratio(igas1, igas2, &
              planck_fraction_interpolation_layer(atm,band))
            CALL get_planck_fractions_interp(jcs, kproma, kbdim, klev, &
              gases(:,:,igas1), ratio, gases(:,:,igas2), &
              fracs(:,:,gpt_low:gpt_high), atm, band, &
              zrange)
          END IF
          IF (associated(kmajor(atm,band)%v)) THEN
            IF (key_species(2,atm,band) == 0) THEN
              igas1 = key_species(1,atm,band)
              CALL get_tau_simple(jcs, kproma, kbdim, klev, atm, band, &
                gpt_low, gpt_high, &
                gases(:,:,igas1), fac, zrange, &
                iabs, kmajor(atm,band)%v, tau)
            ELSE
              igas1 = key_species(1,atm,band)
              igas2 = key_species(2,atm,band)
              ireaction = reaction_table(igas1, igas2)
              IF (atm == 1) THEN
                CALL get_tau_major_lower(jcs, kproma, kbdim, klev, atm, band, &
                  zrange, gpt_low, gpt_high, &
                  gases(:,:,igas1), reaction_ratio(:,:,:,ireaction), &
                  gases(:,:,igas2), &
                  fac, iabs, tau)
              ELSE
                CALL get_tau_major_upper(jcs, kproma, kbdim, klev, atm, band, &
                  zrange, gpt_low, gpt_high, &
                  gases(:,:,igas1), reaction_ratio(:,:,:,ireaction), &
                  gases(:,:,igas2), &
                  fac, iabs, kmajor(atm,band)%v, tau)
              ENDIF
            END IF
          ELSE
            tau(:,zrange(1):zrange(2), gpt_low:gpt_high) = 0
          END IF
        END IF
        DO which = 1,2 ! self/foreign
          IF (h2o_absorption_flag(which,atm,band) == 1) THEN
            CALL get_tau_minor(jcs, kproma, kbdim, klev, zrange, &
              gpt_low, gpt_high, h2o_factor(:,:,which), &
              h2o_fraction(:,:,which), h2o_index(:,:,which), &
              10, &
              h2oref(:,:,which,band), h2oref_delta(:,:,which,band), &
              tau)
          END IF
        END DO
        DO which = 1,n_minor_species(atm,band)
          igas_minor = minor_species(which,atm,band)
          IF (igas_minor <= ngas) THEN
            IF (key_species(2,atm,band) /= 0) THEN
              igas1 = key_species(1,atm,band)
              igas2 = key_species(2,atm,band)
              ratio = planck_ratio(igas1, igas2, &
                minor_species_interpolation_layer(which,atm,band))
            ENDIF
            CALL fill_scale(jcs, kproma, kbdim, klev, which, zrange, atm, &
              band, gases, jp1, coldry, colbrd, scaleminor, &
              scaleminorn2, minorscale)
            IF (key_species(2,atm,band) == 0) THEN
              CALL get_tau_minor(jcs, kproma, kbdim, klev, zrange, &
                gpt_low, gpt_high, minorscale, minorfrac, &
                indminor, kgas2_list(igas_minor,atm,band), &
                kgas2(igas_minor,atm,band)%v, &
                kgas2_delta(igas_minor,atm,band)%v, &
                tau)
            ELSE
              CALL get_tau_minor_spec(jcs, kproma, kbdim, klev, &
                zrange, gpt_low, gpt_high, &
                fracs_mult(atm,band), &
                kgas3_list(:,igas_minor,atm,band), &
                gases(:,:,igas1), ratio, &
                gases(:,:,igas2), minorscale, minorfrac, indminor, &
                kgas3(igas_minor,atm,band)%v, &
                kgas3_delta(igas_minor,atm,band)%v, &
                tau)
            ENDIF
          ELSE !(igas_minor >= first_cfc)
            igas_minor = igas_minor - cfc_offset
            CALL get_tau_cfc(jcs, kproma, kbdim, klev, zrange, &
              gpt_low, gpt_high, igas_minor, wx, &
              cfc(igas_minor,band)%v, tau)
          ENDIF
        ENDDO
        IF (pressure_dependent_tau_correction(atm,band) /= 0) THEN
          CALL pressure_correct_tau(jcs, kproma, kbdim, klev, zrange, gpt_low, &
            gpt_high, pressure_dependent_tau_correction(atm,band), &
            play, tau)
        END IF
        IF (atm == 2 .and. stratosphere_fudge_flag(band) /= 0) THEN
          CALL stratosphere_correction(jcs, kproma, kbdim, klev, band, zrange, &
            gpt_low, gpt_high, tau)
        END IF
      END DO
    END SUBROUTINE do_atmosphere

  !END SUBROUTINE gas_optics_lw

  SUBROUTINE spec_index_1d(jcs, kproma, kbdim, gas1, ratio, gas2, m, &
    parm, comb, js, fs)
    INTEGER, INTENT(IN) :: jcs, kproma, kbdim, m
    REAL(wp), INTENT(IN) :: gas1(KBDIM), ratio, gas2(KBDIM)
    INTEGER, INTENT(OUT) :: js(KBDIM)
    REAL(wp), INTENT(OUT) :: parm(KBDIM), comb(KBDIM), fs(KBDIM)
    REAL(wp) :: mult(KBDIM)

    comb(jcs:kproma) = gas1(jcs:kproma) + ratio * gas2(jcs:kproma)
    parm(jcs:kproma) = MIN(oneminus, gas1(jcs:kproma)/comb(jcs:kproma))
    mult(jcs:kproma) = m * parm(jcs:kproma)
    js(jcs:kproma) = 1 + INT(mult(jcs:kproma))
    fs(jcs:kproma) = MOD1(mult(jcs:kproma))
  END SUBROUTINE spec_index_1d

  SUBROUTINE spec_index_2d(jcs, kproma, kbdim, gas1, ratio, gas2, m, iabs, nsp, &
    parm, comb, js, fs)
    INTEGER, INTENT(IN) :: jcs, kproma, kbdim, iabs(KBDIM), nsp, m
    REAL(wp), INTENT(IN) :: gas1(KBDIM), ratio(KBDIM), gas2(KBDIM)
    INTEGER, INTENT(OUT) :: js(KBDIM)
    REAL(wp), INTENT(OUT) :: parm(KBDIM), comb(KBDIM), fs(KBDIM)
    REAL(wp) :: mult(KBDIM)

    comb(jcs:kproma) = gas1(jcs:kproma) + ratio(jcs:kproma) * gas2(jcs:kproma)
    parm(jcs:kproma) = MIN(oneminus, gas1(jcs:kproma)/comb(jcs:kproma))
    mult(jcs:kproma) = m * parm(jcs:kproma)
    js(jcs:kproma) = 1 + INT(mult(jcs:kproma)) + iabs(jcs:kproma) * nsp 
    fs(jcs:kproma) = MOD1(mult(jcs:kproma))
  END SUBROUTINE spec_index_2d

  SUBROUTINE pressure_correct_tau(jcs, kproma, kbdim, klev, zrange, gpt_low, gpt_high, &
    which, play, tau)
    INTEGER, INTENT(IN) :: jcs, kproma, kbdim, klev, zrange(2), gpt_low, gpt_high, which
    REAL(wp), INTENT(IN) :: play(KBDIM,klev)
    REAL(wp), INTENT(INOUT) :: tau(KBDIM,klev,ngptlw)
    INTEGER :: lay, gpt

    SELECT CASE(which)
      CASE(1)
      DO gpt = gpt_low,gpt_high
      DO lay = zrange(1),zrange(2)
        WHERE(play(jcs:kproma,lay) .LT. 250._wp)
          tau(jcs:kproma,lay,gpt) = tau(jcs:kproma,lay,gpt) * &
            (1._wp - 0.15_wp * (250._wp-play(jcs:kproma,lay)) / 154.4_wp)
        END WHERE
      ENDDO
      ENDDO
      CASE(2)
      DO gpt = gpt_low,gpt_high
      DO lay = zrange(1),zrange(2)
        tau(jcs:kproma,lay,gpt) = tau(jcs:kproma,lay,gpt) * &
          (1._wp - 0.15_wp * (play(jcs:kproma,lay) / 95.6_wp))
      ENDDO
      ENDDO
      CASE(3)
      DO gpt = gpt_low,gpt_high
      DO lay = zrange(1),zrange(2)
        tau(jcs:kproma,lay,gpt) = tau(jcs:kproma,lay,gpt) * &
          (1._wp - .05_wp * (play(jcs:kproma,lay) - 100._wp) / 900._wp)
      ENDDO
      ENDDO
    END SELECT
  END SUBROUTINE pressure_correct_tau


  SUBROUTINE get_tau_cfc(jcs, kproma, kbdim, klev, zrange, &
    gpt_low, gpt_high, icfc, wx, ref, tau)
    INTEGER, INTENT(IN) :: jcs, kproma, kbdim, klev, zrange(2), &
      gpt_low, gpt_high, icfc
    REAL(wp), INTENT(IN) :: &
      wx(KBDIM,ncfc,klev), & ! cross-section amounts (mol/cm2)
      ref(maxperband)
    REAL(wp), INTENT(INOUT) :: tau(KBDIM,klev,ngptlw)
    INTEGER :: gpt, gpt_in, lay
    gpt_in = 0
    DO gpt = gpt_low,gpt_high
      gpt_in = gpt_in+1
      DO lay = zrange(1),zrange(2)
        tau(jcs:kproma,lay,gpt) = tau(jcs:kproma,lay,gpt) + &
          wx(jcs:kproma,icfc,lay) * ref(gpt_in)
      END DO
    END DO
  END SUBROUTINE

  SUBROUTINE stratosphere_correction(jcs, kproma, kbdim, klev, band, zrange, &
    gpt_low, gpt_high, tau)
    INTEGER, INTENT(IN) :: jcs, kproma, kbdim, klev, band, zrange(2), gpt_low, gpt_high
    REAL(wp), INTENT(INOUT) :: tau(KBDIM,klev,ngptlw)
    INTEGER :: gpt, gpt_in
    gpt_in = 0
    DO gpt = gpt_low,gpt_high
      gpt_in = gpt_in+1
      tau(jcs:kproma,zrange(1):zrange(2),gpt) = tau(jcs:kproma,zrange(1):zrange(2),gpt) * &
        stratosphere_fudge(gpt_in,band)
    ENDDO
  END SUBROUTINE stratosphere_correction

  SUBROUTINE get_tau_minor_spec(jcs, kproma, kbdim, klev, zrange, &
    gpt_low, gpt_high, m, dim, gas1, ratio, gas2, &
    scale, fraction, index, ref, delta, tau)
    INTEGER, INTENT(IN) :: jcs, kproma, kbdim, klev, zrange(2), &
      gpt_low, gpt_high, m, dim(2)
    REAL(wp), INTENT(IN) :: ratio
    REAL(wp), DIMENSION(KBDIM,klev), INTENT(IN) :: gas1, gas2, &
      scale, fraction
    REAL(wp), DIMENSION(dim(1),dim(2),maxperband), INTENT(IN) :: &
      ref, delta
    INTEGER, DIMENSION(KBDIM,klev), INTENT(IN) :: index
    REAL(wp), INTENT(INOUT) :: tau(KBDIM,klev,ngptlw)
    INTEGER :: lay, gpt_in, gpt, j(KBDIM), k
    REAL(wp), DIMENSION(KBDIM) :: a, b, f, comb, parm

    gpt_in = 0
    DO gpt = gpt_low,gpt_high
      gpt_in = gpt_in+1
      DO lay = zrange(1),zrange(2)
        CALL spec_index_1d(jcs, kproma, kbdim, gas1(:,lay), ratio, gas2(:,lay), m, &
          parm, comb, j, f)
        DO k = jcs, kproma
          a(k) = ref(j(k),index(k,lay),gpt_in) + &
            f(k) * delta(j(k),index(k,lay),gpt_in)
          b(k) = ref(j(k),index(k,lay)+1,gpt_in) + &
            f(k) * delta(j(k),index(k,lay)+1,gpt_in)
        END DO
        tau(jcs:kproma,lay,gpt) = tau(jcs:kproma,lay,gpt) + &
          scale(jcs:kproma,lay) * (a(jcs:kproma) + fraction(jcs:kproma,lay) * &
          (b(jcs:kproma) - a(jcs:kproma)));
      END DO
    END DO
  END SUBROUTINE get_tau_minor_spec

  SUBROUTINE fill_scale(jcs, kproma, kbdim, klev, which, zrange, atm, band, gases, &
    jp1, coldry, colbrd, scaleminor, scaleminorn2, minorscale)
    INTEGER, INTENT(IN) :: jcs, kproma, kbdim, klev, which, atm, band, zrange(2), &
      jp1(KBDIM,klev)
    REAL(wp), DIMENSION(KBDIM,klev,ngas), INTENT(IN) :: gases
    REAL(wp), DIMENSION(KBDIM,klev), INTENT(IN) :: coldry, colbrd, &
      scaleminor, scaleminorn2
    REAL(wp), INTENT(OUT) :: minorscale(KBDIM,klev)
    INTEGER :: igas, lay
    REAL(wp) :: threshold, a, b, ugly_hack, x(KBDIM), chi(KBDIM)
    REAL(wp), PARAMETER :: fudge = 1e20_wp, fudgem1 = 1e-20_wp

    SELECT CASE(minor_species_scale(which,atm,band))
      CASE(1)
        minorscale(jcs:kproma,zrange(1):zrange(2)) = &
          colbrd(jcs:kproma,zrange(1):zrange(2)) * scaleminorn2(jcs:kproma,zrange(1):zrange(2))
      CASE(2)
        igas = minor_species(which,atm,band)
        minorscale(jcs:kproma,zrange(1):zrange(2)) = &
          gases(jcs:kproma,zrange(1):zrange(2),igas) * scaleminor(jcs:kproma,zrange(1):zrange(2))
      CASE(3)
        minorscale(jcs:kproma,zrange(1):zrange(2)) = &
          colbrd(jcs:kproma,zrange(1):zrange(2)) * scaleminor(jcs:kproma,zrange(1):zrange(2))
      CASE DEFAULT
        igas = minor_species(which,atm,band)
        IF (minor_species_fudge(1,which,atm,band) == 0) THEN
           minorscale(jcs:kproma,zrange(1):zrange(2)) = gases(jcs:kproma,zrange(1):zrange(2),igas)
        ELSE
  !  In atmospheres where the amount of N2O/CO2 is too great to be considered
  !  a minor species, adjust the column amount by an empirical 
  !  factor to obtain the proper contribution.
          threshold = minor_species_fudge(1,which,atm,band)
          a = minor_species_fudge(2,which,atm,band)
          b = minor_species_fudge(3,which,atm,band)
          ugly_hack = minor_species_fudge(4,which,atm,band)
          IF (ugly_hack /= 0) THEN
            chi = ugly_hack
          END IF
          DO lay = zrange(1),zrange(2)
            IF (ugly_hack == 0) THEN
              chi(jcs:kproma) = chi_mls(igas,jp1(jcs:kproma,lay))
            ENDIF
            x(jcs:kproma) = fudge * (gases(jcs:kproma,lay,igas) / coldry(jcs:kproma,lay)) / &
                            chi(jcs:kproma)
            WHERE (x(jcs:kproma) > threshold)
              minorscale(jcs:kproma,lay) = (a + (x(jcs:kproma)-a)**b) * &
                chi(jcs:kproma) * fudgem1 * coldry(jcs:kproma,lay)
            ELSEWHERE
              minorscale(jcs:kproma,lay) = gases(jcs:kproma,lay,igas)
            END WHERE
          END DO
        END IF
    END SELECT
  END SUBROUTINE fill_scale

  SUBROUTINE get_tau_minor(jcs, kproma, kbdim, klev, zrange, gpt_low, gpt_high, &
    factor, fraction, index, nref, ref, delta, tau)
    INTEGER, INTENT(IN) :: jcs, kproma, kbdim, klev, zrange(2), &
      gpt_low, gpt_high, nref
    REAL(wp), DIMENSION(KBDIM,klev), INTENT(IN) :: factor, fraction
    REAL(wp), DIMENSION(nref,maxperband), INTENT(IN) :: ref, delta
    INTEGER, DIMENSION(KBDIM,klev), INTENT(IN) :: index
    REAL(wp), INTENT(INOUT) :: tau(KBDIM,klev,ngptlw)
    INTEGER :: lay, ind(KBDIM), gpt_in, gpt

    gpt_in = 0
    DO gpt = gpt_low, gpt_high
      gpt_in = gpt_in+1
      DO lay = zrange(1),zrange(2)
          ind(jcs:kproma) = index(jcs:kproma,lay)
          tau(jcs:kproma,lay,gpt) = tau(jcs:kproma,lay,gpt) + &
            factor(jcs:kproma,lay) * (&
              ref(ind(jcs:kproma),gpt_in) + &
              fraction(jcs:kproma,lay) * &
              delta(ind(jcs:kproma),gpt_in))
      END DO
    END DO
  END SUBROUTINE get_tau_minor

  SUBROUTINE get_tau_major_lower(jcs, kproma, kbdim, klev, atm, band, &
    zrange, gpt_low, gpt_high, &
    gas1, reaction_ratio, gas2, fac, iabs, tau)
    INTEGER, INTENT(IN) :: jcs, kproma, kbdim, klev, atm, band, zrange(2), &
      gpt_low, gpt_high, iabs(KBDIM,2,2,klev)
    REAL(wp), INTENT(IN) :: gas1(KBDIM,klev), reaction_ratio(KBDIM,2,klev), &
      gas2(KBDIM,klev), &
      fac(KBDIM,2,2,klev)
    REAL(wp), INTENT(INOUT) :: tau(KBDIM,klev,ngptlw)

    REAL(wp), DIMENSION(KBDIM) :: speccomb, specparm, fs, p, p4
    REAL(wp) :: fk(KBDIM), acc(KBDIM)
    INTEGER :: js(KBDIM), i, j(KBDIM), gpt, gpt_in, lay
    INTEGER, PARAMETER :: stencil_6(6,2) = RESHAPE((/ &
      0, 9, 1, 10, 2, 11, &
      1, 10, 0, 9, -1, 8/), SHAPE =(/6,2/)), &

      stencil_4(4) = (/1,10,0,9/)
    INTEGER :: mask(KBDIM)

    gpt_in = 0
    DO gpt = gpt_low, gpt_high
      gpt_in = gpt_in+1
      DO lay = zrange(1),zrange(2)
        DO i = 1,2
          CALL spec_index_2d(jcs, kproma, kbdim, gas1(:,lay), &
            reaction_ratio(:,i,lay), gas2(:,lay), &
            8, iabs(:,i,atm,lay), nsp_species(atm,band), &
            specparm, speccomb, js, fs)
          mask = 0
          WHERE (specparm(jcs:kproma) < 0.125_wp) mask(jcs:kproma) = -1
          WHERE (specparm(jcs:kproma) > 0.875_wp) mask(jcs:kproma) = 1
          WHERE (mask(jcs:kproma) == 0) 
            acc(jcs:kproma) = fs(jcs:kproma) * fac(jcs:kproma,1,i,lay) * &
              kmajor(atm,band)%v(js(jcs:kproma) + stencil_4(1),gpt_in)
            acc(jcs:kproma) = acc(jcs:kproma) + fs(jcs:kproma) * fac(jcs:kproma,2,i,lay) * &
              kmajor(atm,band)%v(js(jcs:kproma) + stencil_4(2),gpt_in)
            fk(jcs:kproma) = 1._wp - fs(jcs:kproma)
            acc(jcs:kproma) = acc(jcs:kproma) + fk(jcs:kproma) * fac(jcs:kproma,1,i,lay) * &
              kmajor(atm,band)%v(js(jcs:kproma) + stencil_4(3),gpt_in)
            acc(jcs:kproma) = acc(jcs:kproma) + fk(jcs:kproma) * fac(jcs:kproma,2,i,lay) * &
              kmajor(atm,band)%v(js(jcs:kproma) + stencil_4(4),gpt_in)
          ENDWHERE
          IF (ANY(mask(jcs:kproma) /= 0)) THEN
            WHERE (mask(jcs:kproma) < 0)
              j(jcs:kproma) = 1
              p(jcs:kproma) = fs(jcs:kproma) - 1
            ELSEWHERE
              j(jcs:kproma) = 2
              p(jcs:kproma) = -fs(jcs:kproma) 
            ENDWHERE
            WHERE (mask(jcs:kproma) /= 0)
              p4(jcs:kproma) = p(jcs:kproma)**4
              acc(jcs:kproma) = p4(jcs:kproma) * fac(jcs:kproma,1,i,lay) * &
                kmajor(atm,band)%v(js(jcs:kproma) + stencil_6(1,j(jcs:kproma)),gpt_in)
              acc(jcs:kproma) = acc(jcs:kproma) + p4(jcs:kproma) * fac(jcs:kproma,2,i,lay) * &
                kmajor(atm,band)%v(js(jcs:kproma) + stencil_6(2,j(jcs:kproma)),gpt_in)
              fk(jcs:kproma) = 1 - p(jcs:kproma) - 2.0_wp*p4(jcs:kproma)
              acc(jcs:kproma) = acc(jcs:kproma) + fk(jcs:kproma) * fac(jcs:kproma,1,i,lay) * &
                kmajor(atm,band)%v(js(jcs:kproma) + stencil_6(3,j(jcs:kproma)),gpt_in)
              acc(jcs:kproma) = acc(jcs:kproma) + fk(jcs:kproma) * fac(jcs:kproma,2,i,lay) * &
                kmajor(atm,band)%v(js(jcs:kproma) + stencil_6(4,j(jcs:kproma)),gpt_in)
              fk(jcs:kproma) = p(jcs:kproma) + p4(jcs:kproma)
              acc(jcs:kproma) = acc(jcs:kproma) + fk(jcs:kproma) * fac(jcs:kproma,1,i,lay) * &
                kmajor(atm,band)%v(js(jcs:kproma) + stencil_6(5,j(jcs:kproma)),gpt_in)
              acc(jcs:kproma) = acc(jcs:kproma) + fk(jcs:kproma) * fac(jcs:kproma,2,i,lay) * &
                kmajor(atm,band)%v(js(jcs:kproma) + stencil_6(6,j(jcs:kproma)),gpt_in)
            ENDWHERE
          END IF
          IF (i == 1) THEN
            tau(jcs:kproma,lay,gpt) = speccomb(jcs:kproma) * acc(jcs:kproma)
          ELSE
            tau(jcs:kproma,lay,gpt) = tau(jcs:kproma,lay,gpt) + speccomb(jcs:kproma) * acc(jcs:kproma)
          ENDIF
        END DO
      END DO
    END DO
  END SUBROUTINE get_tau_major_lower

  SUBROUTINE get_tau_major_upper(jcs, kproma, kbdim, klev, atm, band, &
    zrange, gpt_low, gpt_high, &
    gas1, reaction_ratio, gas2, fac, iabs, ref, tau)
    INTEGER, INTENT(IN) :: jcs, kproma, kbdim, klev, atm, band, zrange(2), &
      gpt_low, gpt_high, iabs(KBDIM,2,2,klev)
    REAL(wp), INTENT(IN) :: gas1(KBDIM,klev), reaction_ratio(KBDIM,2,klev), &
      gas2(KBDIM,klev), fac(KBDIM,2,2,klev), ref(:,:)
    REAL(wp), INTENT(INOUT) :: tau(KBDIM,klev,ngptlw)

    REAL(wp), DIMENSION(KBDIM) :: speccomb, specparm, fs
    REAL(wp) :: acc(KBDIM)
    INTEGER :: js(KBDIM), i, gpt, gpt_in, lay
    INTEGER, PARAMETER :: stencil_4(4) = (/1,6,0,5/)

    gpt_in = 0
    DO gpt = gpt_low, gpt_high
      gpt_in = gpt_in+1
      DO lay = zrange(1),zrange(2)
        DO i = 1,2
          CALL spec_index_2d(jcs, kproma, kbdim, gas1(:,lay), &
            reaction_ratio(:,i,lay), gas2(:,lay), &
            4, iabs(:,i,atm,lay), nsp_species(atm,band), &
            specparm, speccomb, js, fs)
          acc(jcs:kproma) = fs(jcs:kproma) * ( &
              fac(jcs:kproma,1,i,lay) * ref(js(jcs:kproma) + stencil_4(1),gpt_in) + &
              fac(jcs:kproma,2,i,lay) * ref(js(jcs:kproma) + stencil_4(2),gpt_in)) + &
            (1._wp - fs(jcs:kproma)) * ( &
              fac(jcs:kproma,1,i,lay) * ref(js(jcs:kproma) + stencil_4(3),gpt_in) + &
              fac(jcs:kproma,2,i,lay) * ref(js(jcs:kproma) + stencil_4(4),gpt_in))
          IF (i == 1) THEN
            tau(jcs:kproma,lay,gpt) = speccomb(jcs:kproma) * acc(jcs:kproma)
          ELSE
            tau(jcs:kproma,lay,gpt) = tau(jcs:kproma,lay,gpt) + speccomb(jcs:kproma) * acc(jcs:kproma)
          ENDIF
        END DO
      END DO
    END DO
  END SUBROUTINE get_tau_major_upper

  SUBROUTINE get_tau_simple(jcs, kproma, kbdim, klev, atm, band, &
    gpt_low, gpt_high, &
    gas, fac, zrange, iabs, ref, tau)
    INTEGER, INTENT(IN) :: jcs, kproma, kbdim, klev, atm, band, zrange(2), &
      iabs(KBDIM,2,2,klev), gpt_low, gpt_high
    REAL(wp), INTENT(IN) :: gas(KBDIM,klev), &
      fac(KBDIM,2,2,klev), ref(:,:)
    REAL(wp), INTENT(INOUT) :: tau(KBDIM,klev,ngptlw)
    INTEGER :: ind(KBDIM), gpt, gpt_in, lay 
    REAL(wp) :: acc(KBDIM)

    gpt_in = 0
    DO gpt = gpt_low, gpt_high
      gpt_in = gpt_in+1
      DO lay = zrange(1),zrange(2)
        ind(jcs:kproma) = iabs(jcs:kproma,1,atm,lay) * nsp_species_broken_16(atm,band) + 1
        acc(jcs:kproma) = fac(jcs:kproma,1,1,lay) * ref(ind(jcs:kproma),gpt_in)
        acc(jcs:kproma) = acc(jcs:kproma) + fac(jcs:kproma,2,1,lay) * ref(ind(jcs:kproma)+1,gpt_in)
        ind(jcs:kproma) = iabs(jcs:kproma,2,atm,lay) * nsp_species_broken_16(atm,band) + 1
        acc(jcs:kproma) = acc(jcs:kproma) + fac(jcs:kproma,1,2,lay) * ref(ind(jcs:kproma),gpt_in)
        acc(jcs:kproma) = acc(jcs:kproma) + fac(jcs:kproma,2,2,lay) * ref(ind(jcs:kproma)+1,gpt_in)
        tau(jcs:kproma,lay,gpt) = gas(jcs:kproma,lay) * acc(jcs:kproma)
      END DO
    END DO
  END SUBROUTINE get_tau_simple

  SUBROUTINE get_planck_fractions_interp(jcs, kproma, kbdim, klev, gas1, &
    ratio, gas2, fracs, atm, band, zrange)
    INTEGER, INTENT(IN) :: jcs, kproma, kbdim, klev, atm, band, zrange(2)
    REAL(wp), INTENT(IN) :: gas1(KBDIM,klev), gas2(KBDIM,klev), &
      ratio
    REAL(wp), INTENT(INOUT) :: fracs(KBDIM,klev,ngc(band))
    INTEGER :: lay, gpt
    REAL(wp), DIMENSION(KBDIM) :: speccomb_planck, specparm_planck, &
      specmult_planck, fpl
    INTEGER :: jpl(KBDIM)

    DO gpt = 1, ngc(band)
    DO lay = zrange(1),zrange(2)
      speccomb_planck(jcs:kproma) = gas1(jcs:kproma,lay) + ratio * gas2(jcs:kproma,lay)
      specparm_planck(jcs:kproma) = MIN(oneminus, gas1(jcs:kproma,lay)/speccomb_planck(jcs:kproma))
      specmult_planck(jcs:kproma) = fracs_mult(atm,band) * specparm_planck(jcs:kproma)
      jpl(jcs:kproma)= 1 + INT(specmult_planck(jcs:kproma))
      fpl(jcs:kproma) = MOD1(specmult_planck(jcs:kproma))
      fracs(jcs:kproma,lay,gpt) = &
        planck_fraction2(atm,band)%v(jpl(jcs:kproma),gpt) + fpl(jcs:kproma) * &
          planck_fraction2_delta(atm,band)%v(jpl(jcs:kproma),gpt)
    END DO
    END DO
  END SUBROUTINE get_planck_fractions_interp
END MODULE mo_psrad_lrtm_gas_optics

