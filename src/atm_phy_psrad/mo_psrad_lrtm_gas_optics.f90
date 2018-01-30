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

MODULE mo_psrad_lrtm_gas_optics

  USE mo_psrad_general, ONLY: wp, nbndlw, ngptlw, ncfc, ngas, cfc_offset, &
    oneminus, finish, max_minor_species
  USE mo_psrad_lrtm_kgs, ONLY: major_species, minor_species, &
    h2o_absorption_flag, ngpt, fracs_mult, nsp, skip_atmosphere, &
    planck_fraction_interpolation_layer, planck_ratio, &
    pressure_dependent_tau_correction, stratosphere_fudge_idx, &
    minor_species_fudge, chi_mls, stratosphere_fudge
#ifdef PSRAD_DEVEL
  USE mo_psrad_dump, ONLY: toggle_lw_sw, save_tau_major, save_tau_inc, &
    save_correction
#endif
  USE mo_psrad_flat_data, ONLY: flat_data, lw_planck, lw_kmajor, lw_h2oref, &
    lw_kgas
  USE mo_psrad_gas_optics, ONLY: get_tau_major_combined, &
    get_tau_major_simple, get_tau_minor, get_tau_gas
#ifdef PSRAD_NO_INLINE
  USE mo_psrad_gas_optics, ONLY: spec_index_1d
#endif
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: gas_optics_lw

CONTAINS

  SUBROUTINE gas_optics_lw(kproma,kbdim, klev, jp, fac, iabs, &
    laytrop, gases, h2o_factor, h2o_fraction, h2o_index, &
    play, wx, coldry, colbrd, mixing_ratio, &
    minorfrac, scaleminor, scaleminorn2, indminor, fracs_ret, tau_ret)

    USE mo_psrad_general, ONLY: upwards
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: kproma, kbdim, klev 
    INTEGER, INTENT(IN) :: laytrop(KBDIM), & ! tropopause layer index
      iabs(KBDIM,2,2,klev)
    REAL(wp), DIMENSION(KBDIM,klev,2), INTENT(IN) :: h2o_factor,h2o_fraction
    INTEGER, DIMENSION(KBDIM,klev,2), INTENT(IN) :: h2o_index
    INTEGER, DIMENSION(KBDIM,klev), INTENT(IN) :: jp
    REAL(wp), INTENT(IN) :: fac(KBDIM,2,2,klev)
    REAL(wp), TARGET, INTENT(IN) :: gases(:,:,:), mixing_ratio(:,:,:,:)
    INTEGER, DIMENSION(KBDIM,klev), INTENT(IN) :: indminor
    REAL(wp), INTENT(IN) :: &
      wx(KBDIM,klev,ncfc) ! cross-section amounts (mol/cm2)
    REAL(wp), DIMENSION(KBDIM,klev), INTENT(IN) :: &
      play, & ! (klev) layer pressures [mb] 
      coldry, & ! (klev) column amount (dry air)
      colbrd, &
      minorfrac, &
      scaleminor 
    REAL(wp), TARGET, INTENT(IN) :: scaleminorn2(:,:)
    ! Output arrays have size (klev)
    REAL(wp), DIMENSION(KBDIM,klev,ngptlw), INTENT(INOUT) :: &
      fracs_ret, & ! planck fractions
      tau_ret ! gaseous optical depth 

    INTEGER, PARAMETER :: max_merge_size = 5
    REAL(wp), DIMENSION(KBDIM,max_merge_size,ngptlw) :: &
      tau_merge, fracs_merge

    INTEGER :: gpt, i, j, k, merge_size 
    INTEGER :: atm_range(2,2), merge_range(2)
    REAL(wp), TARGET :: actual_corrected_gas(KBDIM,klev)
    CHARACTER(len=128) :: msg

#ifdef PSRAD_DEVEL
    INTEGER :: dump_range(2)
    CALL toggle_lw_sw(.true.)
#endif

    merge_range = (/MINVAL(laytrop(1:kproma)), MAXVAL(laytrop(1:kproma))/)
    merge_size = merge_range(2) - merge_range(1)
    IF (merge_size > 0) THEN
      IF (merge_size > max_merge_size) THEN
        WRITE(msg,'(a,i3,a,i3)') 'Merge_size ', merge_size, &
          'above hardcoded limit of ', max_merge_size
        CALL finish('lrtm_gas_optics', msg, 1)
      ENDIF
      tau_merge(:,1:merge_size,:) = 0
      fracs_merge(:,1:merge_size,:) = 0
      IF (upwards) THEN
        atm_range = RESHAPE((/&
          1, merge_range(2), &
          merge_range(2)+1, klev/), SHAPE=(/2,2/))
      ELSE
        atm_range = RESHAPE((/&
          merge_range(1), klev, &
          1, merge_range(1)-1/), SHAPE=(/2,2/))
      ENDIF
    ELSE
      IF (upwards) THEN
        atm_range = RESHAPE((/&
          1, merge_range(1), &
          merge_range(1)+1, klev/), SHAPE=(/2,2/))
      ELSE
        atm_range = RESHAPE((/&
          merge_range(1), klev, &
          1, merge_range(1)-1/), SHAPE=(/2,2/))
      ENDIF
    END IF
    IF (upwards) THEN
      merge_range(1) = merge_range(1) + 1
    ELSE
      merge_range(2) = merge_range(2) - 1
    ENDIF

#ifdef PSRAD_DEVEL
    IF (upwards) THEN
      dump_range = (/1, merge_range(2)/)
    ELSE
      dump_range = (/merge_range(1), klev/)
    ENDIF
#endif
    CALL do_atmosphere(1, atm_range(:,1), tau_ret, fracs_ret, 0)
#ifdef PSRAD_DEVEL
    IF (upwards) THEN
      dump_range = (/merge_range(2)+1, klev/)
    ELSE
      dump_range = (/1, merge_range(1)-1/)
    ENDIF
#endif
    CALL do_atmosphere(2, atm_range(:,2), tau_ret, fracs_ret, 0)
    IF (merge_size > 0) THEN
#ifdef PSRAD_DEVEL
      dump_range = merge_range
#endif
      CALL do_atmosphere(2, merge_range, &
        tau_merge, fracs_merge, 1-merge_range(1))
      DO gpt = 1,ngptlw
        j = 0
        IF (upwards) THEN
          DO k = merge_range(1),merge_range(2)
            j = j+1
            DO i = 1, kproma
              IF (laytrop(i) < k) THEN
                tau_ret(i,k,gpt) = tau_merge(i,j,gpt)
                fracs_ret(i,k,gpt) = fracs_merge(i,j,gpt)
              ENDIF
            ENDDO
          ENDDO
        ELSE
          DO k = merge_range(1),merge_range(2)
            j = j+1
            DO i = 1, kproma
              IF (laytrop(i) > k) THEN
                tau_ret(i,k,gpt) = tau_merge(i,j,gpt)
                fracs_ret(i,k,gpt) = fracs_merge(i,j,gpt)
              ENDIF
            ENDDO
          ENDDO
        ENDIF
      ENDDO
    ENDIF

  CONTAINS

    SUBROUTINE do_atmosphere(atm, range, tau, fracs, o)
      INTEGER, INTENT(IN) :: atm, range(2), o
      REAL(wp), DIMENSION(:,:,:), INTENT(INOUT) :: &
        fracs, tau
      REAL(wp) :: ratio
      REAL(wp), POINTER :: corrected_gas(:,:)
      INTEGER :: gpt_range(2), gpt_in, & !n_in_band, &
        band, igas1, igas2, imixing, gpt, which, igas_minor

      gpt_range = 0
      DO band = 1,nbndlw
        gpt_range = (/1,ngpt(band)/) + gpt_range(2)
        IF (skip_atmosphere(atm,band) /= 0) THEN
          tau(1:kproma,range(1)+o:range(2)+o, gpt_range(1):gpt_range(2)) = 0
          fracs(1:kproma,range(1)+o:range(2)+o, gpt_range(1):gpt_range(2)) = 0
        ELSE
          IF (planck_fraction_interpolation_layer(atm,band) == 0) THEN
            gpt_in = 0
            DO gpt = gpt_range(1),gpt_range(2)
              gpt_in = gpt_in+1
              fracs(1:kproma,range(1)+o:range(2)+o, gpt) = &
                flat_data(lw_planck(atm,band)%o + gpt_in)
            END DO
          ELSE
            igas1 = major_species(1,atm,band)
            igas2 = major_species(2,atm,band)
            ratio = planck_ratio(igas1, igas2, &
              planck_fraction_interpolation_layer(atm,band))
            CALL get_planck_fractions_interp(kproma, kbdim, klev, atm, &
              band, range, gases(:,:,igas1), ratio, gases(:,:,igas2), &
              lw_planck(atm,band)%o, lw_planck(atm,band)%s, &
              fracs(:,:,gpt_range(1):gpt_range(2)), o)
          END IF
          IF (major_species(2,atm,band) /= 0) THEN
            igas1 = major_species(1,atm,band)
            igas2 = major_species(2,atm,band)
            imixing = major_species(3,atm,band)
            CALL get_tau_major_combined(kproma, kbdim, klev, atm, &
              range, gpt_range, fracs_mult(atm,band), nsp(atm,band), &
              gases(:,:,igas1), gases(:,:,igas2), &
              0.0_wp, imixing, mixing_ratio, &
              fac, iabs, lw_kmajor(atm,band)%o, lw_kmajor(atm,band)%s, &
              atm==1, tau, o)
          ELSEIF (major_species(1,atm,band) /= 0) THEN
            igas1 = major_species(1,atm,band)
            CALL get_tau_major_simple(kproma, kbdim, klev, atm, band, &
              gpt_range, gases(:,:,igas1), fac, range, iabs, &
              lw_kmajor(atm,band)%o, lw_kmajor(atm,band)%s, tau, o)
          ELSE
            tau(1:kproma,range(1)+o:range(2)+o, gpt_range(1):gpt_range(2)) = 0
          ENDIF
#ifdef PSRAD_DEVEL
          CALL save_tau_major(kproma, kbdim, laytrop, dump_range, atm, gpt_range, tau, o)
#endif
          DO which = 1,2 ! self/foreign
            IF (h2o_absorption_flag(which,atm,band) == 1) THEN
              CALL get_tau_minor(kproma, kbdim, klev, laytrop, range, atm, &
                gpt_range, h2o_factor(:,:,which), &
                h2o_fraction(:,:,which), h2o_index(:,:,which), &
                lw_h2oref(which,band)%o, lw_h2oref(which,band)%s, &
                tau, o, which)
            END IF
          END DO
          DO which = 1, max_minor_species
            igas_minor = minor_species(which,atm,band)%gas
            IF (igas_minor == 0) THEN
              EXIT
            ENDIF
            IF (igas_minor <= ngas) THEN
              IF (major_species(2,atm,band) /= 0) THEN
                igas1 = major_species(1,atm,band)
                igas2 = major_species(2,atm,band)
                ratio = planck_ratio(igas1, igas2, &
                  minor_species(which,atm,band)%interpolation_layer)
              ENDIF
              IF (minor_species(which,atm,band)%scaling_type == 0) THEN
                corrected_gas => gases(:,:,minor_species(which,atm,band)%gas)
              ELSE
                corrected_gas => actual_corrected_gas
                CALL correct_gas(kproma, kbdim, klev, which, range, atm, &
                  band, gases, jp, coldry, colbrd, scaleminor, &
                  scaleminorn2, corrected_gas)
              ENDIF
              IF (major_species(2,atm,band) == 0) THEN
                CALL get_tau_minor(kproma, kbdim, klev, laytrop, range, atm, &
                  gpt_range, corrected_gas, minorfrac, &
                  indminor, lw_kgas(which,atm,band)%o, &
                  lw_kgas(which,atm,band)%s, tau, o, which+2)
              ELSE
                CALL get_tau_minor_spec(kproma, kbdim, klev, &
                  laytrop, range, atm, gpt_range, &
                  fracs_mult(atm,band), &
                  gases(:,:,igas1), ratio, &
                  gases(:,:,igas2), corrected_gas, minorfrac, indminor, &
                  lw_kgas(which,atm,band)%o, lw_kgas(which,atm,band)%s, &
                  tau, o, which+2)
              ENDIF
            ELSE !(igas_minor >= first_cfc)
              CALL get_tau_gas(kproma, kbdim, laytrop, range, atm, &
                gpt_range, wx(:,:,igas_minor-cfc_offset), &
                lw_kgas(which,atm,band)%o, tau, o, which+2)
            ENDIF
          ENDDO
          IF (pressure_dependent_tau_correction(atm,band) /= 0) THEN
            CALL pressure_correct_tau(kproma, kbdim, klev, laytrop, range, atm, &
              gpt_range, pressure_dependent_tau_correction(atm,band), &
              play, tau, o)
          END IF
          IF (atm == 2 .and. stratosphere_fudge_idx(band) /= 0) THEN
            CALL stratosphere_correction(kproma, kbdim, laytrop, atm, band, &
              range, gpt_range, tau, o)
          END IF
        ENDIF
      END DO
    END SUBROUTINE do_atmosphere

  END SUBROUTINE gas_optics_lw


  SUBROUTINE correct_gas(kproma, kbdim, klev, which, range, atm, band, &
    gases, jp, coldry, colbrd, scaleminor, scaleminorn2, corrected_gas)
    INTEGER, INTENT(IN) :: kproma, kbdim, klev, which, atm, band, range(2), &
      jp(KBDIM,klev)
    REAL(wp), TARGET, INTENT(IN) :: gases(:,:,:)
    REAL(wp), DIMENSION(KBDIM,klev), INTENT(IN) :: coldry, colbrd, &
      scaleminor
    REAL(wp), TARGET, INTENT(IN) :: scaleminorn2(:,:)
    REAL(wp), POINTER, INTENT(INOUT) :: corrected_gas(:,:)
    INTEGER :: igas, i, lay, ifudge
    REAL(wp) :: threshold, a, b, ugly_hack, x, chi

    SELECT CASE(minor_species(which,atm,band)%scaling_type)
      CASE(1)
        corrected_gas => scaleminorn2
        RETURN
      CASE(2)
        igas = minor_species(which,atm,band)%gas
        corrected_gas(1:kproma,range(1):range(2)) = &
          gases(1:kproma,range(1):range(2),igas) * &
          scaleminor(1:kproma,range(1):range(2)) 
        RETURN
      CASE(3)
        corrected_gas(1:kproma,range(1):range(2)) = &
          colbrd(1:kproma,range(1):range(2)) * &
          scaleminor(1:kproma,range(1):range(2))
        RETURN
      CASE(4)
        igas = minor_species(which,atm,band)%gas
        ifudge = minor_species(which,atm,band)%fudge_index
!  In atmospheres where the amount of N2O/CO2 is too great to be considered
!  a minor species, adjust the column amount by an empirical 
!  factor to obtain the proper contribution.
        threshold = minor_species_fudge(1,ifudge)
        a = minor_species_fudge(2,ifudge)
        b = minor_species_fudge(3,ifudge)
        ugly_hack = minor_species_fudge(4,ifudge)
        IF (ugly_hack /= 0) THEN
          chi = ugly_hack
        END IF
        DO lay = range(1),range(2)
        DO i = 1, kproma
          IF (ugly_hack == 0) THEN
            chi = chi_mls(igas,jp(i,lay)+1)
          ENDIF
          x = (gases(i,lay,igas) / coldry(i,lay)) / chi
          IF (x > threshold) THEN
            corrected_gas(i,lay) = (a + (x-a)**b) * &
              chi * coldry(i,lay)
          ELSE
            corrected_gas(i,lay) = gases(i,lay,igas)
          ENDIF
        END DO
        END DO
        RETURN
    END SELECT
  END SUBROUTINE correct_gas


  SUBROUTINE get_tau_minor_spec(kproma, kbdim, klev, laytrop, range, atm, &
    gpt_range, m, gas1, ratio, gas2, &
    scale, fraction, index, ref_o, ref_s, tau, o, dump_index)
    INTEGER, INTENT(IN) :: kproma, kbdim, klev, laytrop(KBDIM), range(2), atm, &
      gpt_range(2), m, ref_o, ref_s, o, dump_index
    REAL(wp), INTENT(IN) :: ratio
    REAL(wp), DIMENSION(KBDIM,klev), INTENT(IN) :: gas1, gas2, &
      scale, fraction
    INTEGER, DIMENSION(KBDIM,klev), INTENT(IN) :: index
    REAL(wp), INTENT(INOUT) :: tau(:,:,:)
    INTEGER :: i, j, lay, gpt, ptr_base, js(KBDIM)
    REAL(wp) :: a, b, fs(KBDIM)
#ifdef PSRAD_DEVEL
    INTEGER :: gpt_in
    REAL(wp) :: tau_inc(KBDIM)
#endif

    ptr_base = ref_o
#ifdef PSRAD_DEVEL
    gpt_in = 0
#endif
    DO gpt = gpt_range(1),gpt_range(2)
#ifdef PSRAD_DEVEL
      gpt_in = gpt_in+1
#endif
      DO lay = range(1),range(2)
        CALL spec_index_1d(kproma, KBDIM, gas1(:,lay), ratio, gas2(:,lay), &
          m, js, fs)
        js(1:kproma) = js(1:kproma) + ptr_base
        DO i = 1, kproma
          j = js(i) + (index(i,lay)-1) * (m+1)
          a = flat_data(j) + fs(i) * (flat_data(j+1)- flat_data(j))
          j = j + (m+1)
          b = flat_data(j) + fs(i) * (flat_data(j+1)- flat_data(j))
#ifdef PSRAD_DEVEL
          tau_inc(i) = scale(i,lay) * (a + fraction(i,lay) * (b - a));
          tau(i,lay+o,gpt) = tau(i,lay+o,gpt) + tau_inc(i)
#else
          tau(i,lay+o,gpt) = tau(i,lay+o,gpt) + &
            scale(i,lay) * (a + fraction(i,lay) * (b - a));
#endif
        END DO
#ifdef PSRAD_DEVEL
        CALL save_tau_inc(kproma, kbdim, laytrop, lay, atm, gpt, dump_index, tau_inc)
#endif
      END DO
      ptr_base = ptr_base + ref_s
    END DO
  END SUBROUTINE get_tau_minor_spec

  SUBROUTINE pressure_correct_tau(kproma, kbdim, klev, laytrop, range, atm, &
    gpt_range, which, play, tau, o)

    USE mo_psrad_lrtm_kgs, ONLY: pa1, pa2, pa3, pb1, pb2, pc1, pc2, pc3

    INTEGER, INTENT(IN) :: kproma, kbdim, klev, laytrop(KBDIM), &
      range(2), atm, gpt_range(2), which, o
    REAL(wp), INTENT(IN) :: play(KBDIM,klev)
    REAL(wp), INTENT(INOUT) :: tau(:,:,:)
    INTEGER :: i, lay, gpt
#ifdef PSRAD_DEVEL
    REAL(wp) :: correction(KBDIM)
#endif

    SELECT CASE(which)
      CASE(1)
      DO gpt = gpt_range(1),gpt_range(2)
      DO lay = range(1),range(2)
        DO i = 1, kproma
#ifdef PSRAD_DEVEL
          IF (play(i,lay) < pa2) THEN
            correction(i) = 1._wp - &
              pa1 * (pa2 - play(i,lay)) / pa3
            tau(i,lay+o,gpt) = tau(i,lay+o,gpt) * correction(i)
          ELSE
            correction(i) = 1
          ENDIF
#else
          IF (play(i,lay) < pa2 ) THEN
            tau(i,lay+o,gpt) = tau(i,lay+o,gpt) * &
              (1._wp - pa1 * (pa2 - play(i,lay)) / pa3)
          ENDIF
#endif
        ENDDO
#ifdef PSRAD_DEVEL
        CALL save_correction(kproma, kbdim, laytrop, (/lay,lay/), atm, gpt, correction)
#endif
      ENDDO
      ENDDO

      CASE(2)
      DO gpt = gpt_range(1),gpt_range(2)
      DO lay = range(1),range(2)
        DO i = 1, kproma
#ifdef PSRAD_DEVEL
          correction(i) = (1._wp - pb1 * (play(i,lay) / pb2))
          tau(i,lay+o,gpt) = tau(i,lay+o,gpt) * correction(i)
#else
          tau(i,lay+o,gpt) = tau(i,lay+o,gpt) * &
            (1._wp - pb1 * (play(i,lay) / pb2))
#endif
        ENDDO
#ifdef PSRAD_DEVEL
        CALL save_correction(kproma, kbdim, laytrop, (/lay,lay/), atm, gpt, correction)
#endif
      ENDDO
      ENDDO

      CASE(3)
      DO gpt = gpt_range(1),gpt_range(2)
      DO lay = range(1),range(2)
        DO i = 1, kproma
#ifdef PSRAD_DEVEL
          correction(i) = &
            (1._wp - pc1 * (play(i,lay) - pc2) / pc3)
          tau(i,lay+o,gpt) = tau(i,lay+o,gpt) * correction(i)
#else
          tau(i,lay+o,gpt) = tau(i,lay+o,gpt) * &
            (1._wp - pc1 * (play(i,lay) - pc2) / pc3)
#endif
        ENDDO
#ifdef PSRAD_DEVEL
        CALL save_correction(kproma, kbdim, laytrop, (/lay,lay/), atm, gpt, correction)
#endif
      ENDDO
      ENDDO
    END SELECT
  END SUBROUTINE pressure_correct_tau

  SUBROUTINE stratosphere_correction(kproma, kbdim, laytrop, atm, band, range, &
    gpt_range, tau, o)
    INTEGER, INTENT(IN) :: kproma, kbdim, laytrop(KBDIM), atm, band, &
      range(2), gpt_range(2), o
    REAL(wp), INTENT(INOUT) :: tau(:,:,:)
    INTEGER :: gpt, gpt_in, idx
#ifdef PSRAD_DEVEL
    REAL(wp), DIMENSION(KBDIM) :: correction
#endif
    gpt_in = 0
    idx = stratosphere_fudge_idx(band)
    DO gpt = gpt_range(1),gpt_range(2)
      gpt_in = gpt_in+1
#ifdef PSRAD_DEVEL
      correction(1:kproma) = stratosphere_fudge(gpt_in,idx)
      CALL save_correction(kproma, kbdim, laytrop, range, atm, gpt, correction)
#endif
      tau(1:kproma,range(1)+o:range(2)+o,gpt) = &
        tau(1:kproma,range(1)+o:range(2)+o,gpt) * &
        stratosphere_fudge(gpt_in,idx)
    ENDDO
  END SUBROUTINE stratosphere_correction


  SUBROUTINE get_planck_fractions_interp(kproma, kbdim, klev, atm, band, &
    range, gas1, ratio, gas2, ref_o, ref_s, fracs, o)
    INTEGER, INTENT(IN) :: kproma, kbdim, klev, atm, band, range(2), &
      ref_o, ref_s, o
    REAL(wp), INTENT(IN) :: gas1(KBDIM,klev), gas2(KBDIM,klev), &
      ratio
    REAL(wp), INTENT(INOUT) :: fracs(:,:,:)
    INTEGER :: lay, gpt, j(KBDIM), i, ptr_base
    REAL(wp) :: f(KBDIM)

    ptr_base = ref_o
    DO gpt = 1, ngpt(band)
      DO lay = range(1),range(2)
        CALL spec_index_1d(kproma, KBDIM, gas1(:,lay), ratio, gas2(:,lay), &
          fracs_mult(atm,band), j, f)
        j(1:kproma) = j(1:kproma) + ptr_base
        DO i = 1, kproma
          fracs(i,lay+o,gpt) = flat_data(j(i)) + &
            f(i) * (flat_data(j(i)+1) - flat_data(j(i)))
        END DO
      END DO
      ptr_base = ptr_base + ref_s
    END DO
  END SUBROUTINE get_planck_fractions_interp

#ifndef PSRAD_NO_INLINE
!DEC$ ATTRIBUTES FORCEINLINE :: spec_index_1d
  SUBROUTINE spec_index_1d(kproma, kbdim, gas1, ratio, gas2, m, js, fs)
    INTEGER, INTENT(IN) :: kproma, kbdim, m
    REAL(wp), INTENT(IN) :: gas1(KBDIM), ratio, gas2(KBDIM)
    INTEGER, INTENT(INOUT) :: js(KBDIM)
    REAL(wp), INTENT(INOUT) :: fs(KBDIM)
    REAL(wp) :: mult(KBDIM), parm(KBDIM), comb(KBDIM)

    comb(1:kproma) = gas1(1:kproma) + ratio * gas2(1:kproma)
    parm(1:kproma) = MIN(oneminus, gas1(1:kproma)/comb(1:kproma))
    mult(1:kproma) = m * parm(1:kproma)
    js(1:kproma) = INT(mult(1:kproma))
    fs(1:kproma) = mult(1:kproma) - js(1:kproma)
    js(1:kproma) = js(1:kproma) + 1
  END SUBROUTINE spec_index_1d
#endif

END MODULE mo_psrad_lrtm_gas_optics

