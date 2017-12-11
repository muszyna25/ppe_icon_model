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

  USE mo_psrad_general
  USE mo_psrad_lrtm_kgs

  USE mo_psrad_gas_optics, ONLY: get_tau_major_combined, &
    get_tau_major_simple, get_tau_minor, get_tau_gas
#ifdef PSRAD_NO_INLINE
  USE mo_psrad_gas_optics, ONLY : spec_index_1d
#endif
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: gas_optics_lw

CONTAINS

  SUBROUTINE gas_optics_lw(kproma,kbdim, klev, jp, fac, iabs, &
    laytrop, gases, h2o_factor, h2o_fraction, h2o_index, &
    play, wx, coldry, colbrd, mixing_ratio, &
    minorfrac, scaleminor, scaleminorn2, indminor, fracs_ret, tau_ret)

    USE mo_psrad_general, ONLY : upwards
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: kproma, kbdim, klev 
    INTEGER, INTENT(IN) :: laytrop(KBDIM), & ! tropopause layer index
      iabs(KBDIM,2,2,klev)
    REAL(wp), DIMENSION(KBDIM,klev,2), INTENT(IN) :: h2o_factor,h2o_fraction
    INTEGER, DIMENSION(KBDIM,klev,2), INTENT(IN) :: h2o_index
    INTEGER, DIMENSION(KBDIM,klev), INTENT(IN) :: jp
    REAL(wp), INTENT(IN) :: fac(KBDIM,2,2,klev)
    REAL(wp), POINTER, INTENT(IN) :: gases(:,:,:), mixing_ratio(:,:,:,:)
    INTEGER, DIMENSION(KBDIM,klev), INTENT(IN) :: indminor
    REAL(wp), INTENT(IN) :: &
      wx(KBDIM,klev,ncfc) ! cross-section amounts (mol/cm2)
    REAL(wp), DIMENSION(KBDIM,klev), INTENT(IN) :: &
      play, & ! (klev) layer pressures [mb] 
      coldry, & ! (klev) column amount (dry air)
      colbrd, &
      minorfrac, &
      scaleminor 
    REAL(wp), POINTER, INTENT(IN) :: scaleminorn2(:,:)
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

    CALL do_atmosphere(1, atm_range(:,1), tau_ret, fracs_ret, 0)
    CALL do_atmosphere(2, atm_range(:,2), tau_ret, fracs_ret, 0)
    IF (merge_size > 0) THEN
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
                planck_fraction2(atm,band)%v(1,gpt_in)
            END DO
          ELSE
            igas1 = major_species(1,atm,band)
            igas2 = major_species(2,atm,band)
            ratio = planck_ratio(igas1, igas2, &
              planck_fraction_interpolation_layer(atm,band))
            CALL get_planck_fractions_interp(kproma, kbdim, atm, &
              band, range, gases(:,:,igas1), ratio, gases(:,:,igas2), &
              planck_fraction2(atm,band)%v, &
              fracs(:,:,gpt_range(1):gpt_range(2)), o)
          END IF
          IF (major_species(2,atm,band) /= 0) THEN
            igas1 = major_species(1,atm,band)
            igas2 = major_species(2,atm,band)
            imixing = major_species(3,atm,band)
            CALL get_tau_major_combined(kproma, kbdim, atm, &
              range, gpt_range, fracs_mult(atm,band), nsp(atm,band), &
              gases(:,:,igas1), gases(:,:,igas2), &
              0.0_wp, imixing, mixing_ratio, &
              fac, iabs, kmajor(atm,band)%v, atm==1, tau, o)
          ELSEIF (major_species(1,atm,band) /= 0) THEN
            igas1 = major_species(1,atm,band)
            CALL get_tau_major_simple(kproma, atm, band, &
              gpt_range, &
              gases(:,:,igas1), fac, range, &
              iabs, kmajor(atm,band)%v, tau, o)
          ELSE
            tau(1:kproma,range(1)+o:range(2)+o, gpt_range(1):gpt_range(2)) = 0
          ENDIF
          DO which = 1,2 ! self/foreign
            IF (h2o_absorption_flag(which,atm,band) == 1) THEN
              CALL get_tau_minor(kproma, range, &
                gpt_range, h2o_factor(:,:,which), &
                h2o_fraction(:,:,which), h2o_index(:,:,which), &
                h2oref(which,band)%v, tau, o)
            END IF
          END DO
          DO which = 1,n_minor_species(atm,band)
            igas_minor = minor_species(which,atm,band)%gas
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
                CALL correct_gas(kproma, which, range, atm, &
                  band, gases, jp, coldry, colbrd, scaleminor, &
                  scaleminorn2, corrected_gas)
              ENDIF
              IF (major_species(2,atm,band) == 0) THEN
                CALL get_tau_minor(kproma, range, &
                  gpt_range, corrected_gas, minorfrac, &
                  indminor, kgas(which,atm,band)%v, tau, o)
              ELSE
                CALL get_tau_minor_spec(kproma, kbdim, &
                  range, gpt_range, &
                  fracs_mult(atm,band), &
                  gases(:,:,igas1), ratio, &
                  gases(:,:,igas2), corrected_gas, minorfrac, indminor, &
                  kgas(which,atm,band)%v, tau, o)
              ENDIF
            ELSE !(igas_minor >= first_cfc)
              CALL get_tau_gas(kproma, range, &
                gpt_range, wx(:,:,igas_minor-cfc_offset), &
                kgas(which,atm,band)%v(1,:), tau, o)
            ENDIF
          ENDDO
          IF (pressure_dependent_tau_correction(atm,band) /= 0) THEN
            CALL pressure_correct_tau(kproma, range, &
              gpt_range, pressure_dependent_tau_correction(atm,band), &
              play, tau, o)
          END IF
          IF (atm == 2 .and. stratosphere_fudge_flag(band) /= 0) THEN
            CALL stratosphere_correction(kproma, band, &
              range, gpt_range, tau, o)
          END IF
        ENDIF
      END DO
    END SUBROUTINE do_atmosphere

  END SUBROUTINE gas_optics_lw


  SUBROUTINE correct_gas(kproma, which, range, atm, band, &
    gases, jp, coldry, colbrd, scaleminor, scaleminorn2, corrected_gas)
    INTEGER, INTENT(IN) :: kproma, which, atm, band, range(2), &
      jp(:,:)
    REAL(wp), POINTER, INTENT(IN) :: gases(:,:,:)
    REAL(wp), DIMENSION(:,:), INTENT(IN) :: coldry, colbrd, &
      scaleminor
    REAL(wp), POINTER, INTENT(IN) :: scaleminorn2(:,:)
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


  SUBROUTINE get_tau_minor_spec(kproma, kbdim, range, &
    gpt_range, m, gas1, ratio, gas2, &
    scale, fraction, index, ref, tau, o)
    INTEGER, INTENT(IN) :: kproma, kbdim, range(2), &
      gpt_range(2), m, o
    REAL(wp), INTENT(IN) :: ratio
    REAL(wp), DIMENSION(:,:), INTENT(IN) :: gas1, gas2, &
      scale, fraction
    REAL(wp), DIMENSION(:,:), INTENT(IN) :: ref
    INTEGER, DIMENSION(:,:), INTENT(IN) :: index
    REAL(wp), INTENT(INOUT) :: tau(:,:,:)
    INTEGER :: i, j, lay, gpt_in, gpt, js(kbdim)
    REAL(wp) :: a, b, fs(kbdim)

    gpt_in = 0
    DO gpt = gpt_range(1),gpt_range(2)
      gpt_in = gpt_in+1
      DO lay = range(1),range(2)
        CALL spec_index_1d(kproma, kbdim, gas1(:,lay), ratio, gas2(:,lay), &
          m, js, fs)
        DO i = 1, kproma
          j = js(i) + (index(i,lay)-1) * (m+1)
          a = ref(j,gpt_in) + fs(i) * (ref(j+1,gpt_in)- ref(j,gpt_in))
          j = j + (m+1)
          b = ref(j,gpt_in) + fs(i) * (ref(j+1,gpt_in)- ref(j,gpt_in))
          tau(i,lay+o,gpt) = tau(i,lay+o,gpt) + &
            scale(i,lay) * (a + fraction(i,lay) * (b - a));
        END DO
      END DO
    END DO
  END SUBROUTINE get_tau_minor_spec

  SUBROUTINE pressure_correct_tau(kproma, range, &
    gpt_range, which, play, tau, o)

    USE mo_psrad_lrtm_setup, ONLY : pa1, pa2, pa3, pb1, pb2, pc1, pc2, pc3

    INTEGER, INTENT(IN) :: kproma, &
      range(2), gpt_range(2), which, o
    REAL(wp), INTENT(IN) :: play(:,:)
    REAL(wp), INTENT(INOUT) :: tau(:,:,:)
    INTEGER :: i, lay, gpt

    SELECT CASE(which)
      CASE(1)
      DO gpt = gpt_range(1),gpt_range(2)
      DO lay = range(1),range(2)
        DO i = 1, kproma
          IF (play(i,lay) < pa2 ) THEN
            tau(i,lay+o,gpt) = tau(i,lay+o,gpt) * &
              (1._wp - pa1 * (pa2 - play(i,lay)) / pa3)
          ENDIF
        ENDDO
      ENDDO
      ENDDO

      CASE(2)
      DO gpt = gpt_range(1),gpt_range(2)
      DO lay = range(1),range(2)
        DO i = 1, kproma
          tau(i,lay+o,gpt) = tau(i,lay+o,gpt) * &
            (1._wp - pb1 * (play(i,lay) / pb2))
        ENDDO
      ENDDO
      ENDDO

      CASE(3)
      DO gpt = gpt_range(1),gpt_range(2)
      DO lay = range(1),range(2)
        DO i = 1, kproma
          tau(i,lay+o,gpt) = tau(i,lay+o,gpt) * &
            (1._wp - pc1 * (play(i,lay) - pc2) / pc3)
        ENDDO
      ENDDO
      ENDDO
    END SELECT
  END SUBROUTINE pressure_correct_tau

  SUBROUTINE stratosphere_correction(kproma, band, range, &
    gpt_range, tau, o)
    INTEGER, INTENT(IN) :: kproma, band, &
      range(2), gpt_range(2), o
    REAL(wp), INTENT(INOUT) :: tau(:,:,:)
    INTEGER :: gpt, gpt_in
    gpt_in = 0
    DO gpt = gpt_range(1),gpt_range(2)
      gpt_in = gpt_in+1
      tau(1:kproma,range(1)+o:range(2)+o,gpt) = &
        tau(1:kproma,range(1)+o:range(2)+o,gpt) * &
        stratosphere_fudge(gpt_in,band)
    ENDDO
  END SUBROUTINE stratosphere_correction


  SUBROUTINE get_planck_fractions_interp(kproma, kbdim, atm, band, &
    range, gas1, ratio, gas2, ref, fracs, o)
    INTEGER, INTENT(IN) :: kproma, kbdim, atm, band, range(2), o
    REAL(wp), INTENT(IN) :: gas1(:,:), gas2(:,:), &
      ratio, ref(:,:)
    REAL(wp), INTENT(INOUT) :: fracs(:,:,:)
    INTEGER :: lay, gpt, j(KBDIM), i
    REAL(wp) :: f(KBDIM)

    DO gpt = 1, ngpt(band)
    DO lay = range(1),range(2)
      CALL spec_index_1d(kproma, KBDIM, gas1(:,lay), ratio, gas2(:,lay), &
        fracs_mult(atm,band), j, f)
      DO i = 1, kproma
        fracs(i,lay+o,gpt) = ref(j(i),gpt) + &
          f(i) * (ref(j(i)+1,gpt) - ref(j(i),gpt))
      END DO
    END DO
    END DO
  END SUBROUTINE get_planck_fractions_interp

  SUBROUTINE spec_index_1d(kproma, kbdim, gas1, ratio, gas2, m, js, fs)
    INTEGER, INTENT(IN) :: kproma, kbdim, m
    REAL(wp), INTENT(IN) :: gas1(KBDIM), ratio, gas2(KBDIM)
    INTEGER, INTENT(INOUT) :: js(KBDIM)
    REAL(wp), INTENT(INOUT) :: fs(KBDIM)
    REAL(wp) :: mult(KBDIM), parm(KBDIM), comb(KBDIM)

    INTEGER :: i

    DO i = 1, kproma
      comb(i) = gas1(i) + ratio * gas2(i)
      parm(i) = MIN(oneminus, gas1(i)/comb(i))
      mult(i) = m * parm(i)
      js(i) = INT(mult(i))
      fs(i) = mult(i) - js(i)
      js(i) = js(i) + 1
    ENDDO
  END SUBROUTINE spec_index_1d


END MODULE mo_psrad_lrtm_gas_optics

