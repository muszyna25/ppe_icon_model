!> !! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_psrad_srtm_gas_optics

  USE mo_psrad_general, ONLY: wp, nbndsw, ngptsw, oneminus, io2, ih2o, &
    oneminus, dummy4, finish
  USE mo_psrad_srtm_kgs, ONLY: rayl_type, nsfluxref, &
    major_species, minor_species, h2o_absorption_flag, ngpt, fracs_mult, nsp
  USE mo_psrad_flat_data, ONLY: flat_data, sw_kmajor, sw_h2oref, sw_kgas, &
    sw_sfluxref, sw_rayl
  USE mo_psrad_gas_optics, ONLY: get_tau_major_combined, &
    get_tau_major_simple, get_tau_minor, get_tau_gas
#ifdef PSRAD_NO_INLINE
  USE mo_psrad_gas_optics, ONLY: spec_index_1d
#endif
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: gas_optics_sw

  REAL(wp), PARAMETER :: mixing_ratio(2,nbndsw) = RESHAPE((/&
      252.131_wp, 0.0_wp, &
      0.364641_wp, 0.364641_wp, &
      38.9589_wp, 0.0_wp, &
      5.49281_wp, 0.0_wp, &
      0.0_wp, 0.0_wp, &

      0.0045321_wp, 0.0045321_wp, &
      0.022708_wp * 1.6_wp, 0.0_wp, & !NOTE: Notice 1.6 factor (band 7)
      0.0_wp, 0.0_wp, &
      0.124692_wp, 0.0_wp, &
      0.0_wp, 0.0_wp, &

      0.0_wp, 0.0_wp, &
      0.0_wp, 0.0_wp, &
      6.67029e-07_wp, 6.67029e-07_wp, &
      0.0_wp, 0.0_wp/), SHAPE=(/2,nbndsw/))

CONTAINS

  SUBROUTINE gas_optics_sw(kproma, kbdim, klev, &
    jp, fac, iabs, laytrop, &
    gases, h2o_factor, h2o_fraction, h2o_index, colmol, &
    sflx_zen, tau_gas_ret, tau_aer_ret)

    USE mo_psrad_general, ONLY: upwards

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: kproma, kbdim, klev 
    INTEGER, INTENT(IN) :: laytrop(KBDIM), & ! tropopause layer index
      iabs(KBDIM,2,2,klev)
    REAL(wp), DIMENSION(KBDIM,klev,2), INTENT(IN) :: h2o_factor,h2o_fraction
    INTEGER, DIMENSION(KBDIM,klev,2), INTENT(IN) :: h2o_index
    INTEGER, DIMENSION(KBDIM,klev), INTENT(in) :: jp
    REAL(wp), INTENT(IN) :: fac(KBDIM,2,2,klev), colmol(:,:)
    REAL(wp), TARGET, INTENT(IN) :: gases(:,:,:)
    REAL(wp), INTENT(INOUT) :: sflx_zen(:,:), tau_gas_ret(:,:,:), &
      tau_aer_ret(:,:,:)

    INTEGER :: gpt, i, j, k, merge_size
    INTEGER, PARAMETER :: max_merge_size = 5
    REAL(wp), DIMENSION(KBDIM,klev,ngptsw) :: &
      tau_gas_merge, tau_aer_merge
    INTEGER :: atm_range(2,2), merge_range(2)
    CHARACTER(len=128) :: msg

    merge_range = (/MINVAL(laytrop(1:kproma)), MAXVAL(laytrop(1:kproma))/)
    merge_size = merge_range(2) - merge_range(1)
    IF (merge_size > 0) THEN
      IF (merge_size > max_merge_size) THEN
        WRITE(msg,'(a,i3,a,i3)') 'Merge_size ', merge_size, &
          'above hardcoded limit of ', max_merge_size
        CALL finish('lrtm_gas_optics', msg, 1)
      ENDIF
      tau_gas_merge(:,1:merge_size,:) = 0
      tau_aer_merge(:,1:merge_size,:) = 0
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
    CALL do_atmosphere(1, atm_range(:,1), tau_gas_ret, tau_aer_ret, 0)
    CALL do_atmosphere(2, atm_range(:,2), tau_gas_ret, tau_aer_ret, 0)
    IF (merge_size > 0) THEN
      CALL do_atmosphere(2, merge_range, &
        tau_gas_merge, tau_aer_merge, 1-merge_range(1))
      DO gpt = 1,ngptsw
      j = 0
      DO k = merge_range(1),merge_range(2)
        j = j+1
        IF (upwards) THEN
          DO i = 1, kproma
            IF (laytrop(i) < k) THEN
              tau_gas_ret(i,k,gpt) = tau_gas_merge(i,j,gpt)
              tau_aer_ret(i,k,gpt) = tau_aer_merge(i,j,gpt)
            ENDIF
          END DO
        ELSE
          DO i = 1, kproma
            IF (laytrop(i) > k) THEN
              tau_gas_ret(i,k,gpt) = tau_gas_merge(i,j,gpt)
              tau_aer_ret(i,k,gpt) = tau_aer_merge(i,j,gpt)
            ENDIF
          END DO
        ENDIF
      END DO
      END DO
    END IF

    CALL get_solar_flux(kproma, KBDIM, klev, &
      jp, laytrop, gases, sflx_zen)

  CONTAINS

    SUBROUTINE do_atmosphere(atm, range, tau_gas, tau_aer, o)
      INTEGER, INTENT(IN) :: atm, range(2), o
      REAL(wp), DIMENSION(:,:,:), INTENT(INOUT) :: &
        tau_gas, tau_aer
      INTEGER :: gpt_range(2), band, igas1, igas2, which, igas_minor

      gpt_range = 0
      DO band = 1,nbndsw
        gpt_range = (/1,ngpt(band)/) + gpt_range(2)
        IF (major_species(2,atm,band) /= 0) THEN
          igas1 = major_species(1,atm,band)
          igas2 = major_species(2,atm,band)
          CALL get_tau_major_combined(kproma, kbdim, klev, atm, &
            range, gpt_range, fracs_mult(atm,band), nsp(atm,band), &
            gases(:,:,igas1), gases(:,:,igas2), &
            mixing_ratio(atm,band), 0, dummy4, &
            fac, iabs, sw_kmajor(atm,band)%o, sw_kmajor(atm,band)%s, &
            .false., tau_gas, o)
        ELSEIF (major_species(1,atm,band) /= 0) THEN
          igas1 = major_species(1,atm,band)
          CALL get_tau_major_simple(kproma, kbdim, klev, atm, band, &
            gpt_range, gases(:,:,igas1), fac, range, iabs, &
            sw_kmajor(atm,band)%o, sw_kmajor(atm,band)%s, tau_gas, o)
        ELSE
          tau_gas(1:kproma,range(1)+o:range(2)+o, gpt_range(1):gpt_range(2)) = 0
        ENDIF
        IF (h2o_absorption_flag(atm,band) > 0) THEN
          DO which = 1,2 ! self/foreign
            CALL get_tau_minor(kproma, kbdim, klev, laytrop, range, atm, &
              gpt_range, h2o_factor(:,:,which), &
              h2o_fraction(:,:,which), h2o_index(:,:,which), &
              sw_h2oref(which,band)%o, sw_h2oref(which,band)%s, &
              tau_gas, o, which)
          END DO
        END IF
        igas_minor = minor_species(atm,band)
        IF (igas_minor /= 0) THEN
          CALL get_tau_gas(kproma, kbdim, laytrop, range, atm, &
            gpt_range, gases(:,:,igas_minor), &
            sw_kgas(atm,band)%o, tau_gas, o, 3)
        ENDIF
        CALL get_tau_aer(kproma, kbdim, atm, band, laytrop, range, &
          gpt_range, colmol, gases, tau_aer, o)
      END DO
    END SUBROUTINE do_atmosphere

  END SUBROUTINE gas_optics_sw

  SUBROUTINE get_tau_aer(kproma, kbdim, atm, band, laytrop, range, &
    gpt_range, colmol, gases, tau_aer, o)
    INTEGER, INTENT(IN) :: kproma, kbdim, atm, band, laytrop(KBDIM), &
      range(2), gpt_range(2), o
    REAL(wp), INTENT(IN) :: colmol(:,:), gases(:,:,:)
    REAL(wp), INTENT(INOUT) :: tau_aer(:,:,:)
    INTEGER :: i, lay, gpt, gpt_in, js(KBDIM)
    REAL(wp) :: fs(KBDIM)

    IF (rayl_type(atm,band) == 0) THEN
      gpt_in = 0
      DO gpt = gpt_range(1),gpt_range(2)
        gpt_in = gpt_in+1
        tau_aer(1:kproma,range(1)+o:range(2)+o,gpt) = &
          colmol(1:kproma,range(1):range(2)) * &
          flat_data(sw_rayl(atm,band)%o+1)
      ENDDO
    ELSEIF (rayl_type(atm,band) == 1) THEN
      gpt_in = 0
      DO gpt = gpt_range(1),gpt_range(2)
        gpt_in = gpt_in+1
        tau_aer(1:kproma,range(1)+o:range(2)+o,gpt) = &
          colmol(1:kproma,range(1):range(2)) * &
          flat_data(sw_rayl(atm,band)%o+gpt_in)
      END DO
    ELSE !this is a hack for band 9, atm 1
      DO lay = range(1),range(2)
        CALL spec_index_1d(kproma, kbdim, gases(:,lay,ih2o), &
          mixing_ratio(1,9), gases(:,lay,io2), fracs_mult(1,9), js, fs)
        gpt_in = 0
        js = js + sw_rayl(atm,band)%o
        DO gpt = gpt_range(1),gpt_range(2)
          gpt_in = gpt_in+1
          !THIS IMPLIES the branch where fs was set was taken!
          DO i = 1, kproma
            tau_aer(i,lay,gpt) = colmol(i,lay) * &
              (flat_data(js(i)) + fs(i) * &
              (flat_data(js(i)+1) - flat_data(js(i))))
          ENDDO
          js = js + sw_rayl(atm,band)%s
        ENDDO
      ENDDO
    ENDIF
  END SUBROUTINE get_tau_aer

! FIXME: desperately in need of clearing up, kept as is for backward
! compatibility
  SUBROUTINE get_solar_flux(kproma, kbdim, klev, jp, laytrop, &
    gases, sflx_zen)
    USE mo_psrad_general, ONLY: upwards
    INTEGER, INTENT(IN) :: kproma, kbdim, klev, jp(KBDIM,klev), laytrop(KBDIM)
    REAL(wp), INTENT(IN) :: gases(:,:,:)
    REAL(wp), INTENT(INOUT) :: sflx_zen(:,:)
    INTEGER  :: i, gpt, gpt_in, jk, js_save, band, atm, &
      jk_below, jk_above, laysolfr, ptr
    REAL(wp) :: specparm, specmult, speccomb, fs_save
    INTEGER :: jFROM, jTO, jINC

    INTEGER, PARAMETER :: &
      !NOTE: negative values make no sense. This will be compared
      !to 1 <= jp(:,:) <= 58
      !(/-18,30,6,3,-3, 8,2,-6,1,-2, -1,-32,58,-49/), &
      pref_effective(nbndsw) = (/0,30,6,3,0, 8,2,0,1,0, 0,0,58,0/), &
      laysolfr_ref(nbndsw) = & !NOTE: originally set at runtime with 1=>klev
      (/1,1,0,0,0, 0,0,0,0,0, 0,1,1,1/)

    gpt = 0
    DO band = 1,nbndsw
      ptr = sw_sfluxref(band)%o
      DO gpt_in = 1,ngpt(band)
        gpt = gpt + 1
        DO i = 1,kproma
          js_save = 0
          IF (upwards) THEN
            laysolfr = laytrop(i) + laysolfr_ref(band)
          ELSE
            laysolfr = laytrop(i) - laysolfr_ref(band)
          ENDIF
          DO atm = 1,2
            IF (major_species(2,atm,band) > 0) THEN
              IF (atm == 1) THEN
                IF (upwards) THEN
                  jFROM = 1
                  jTO = laytrop(i)
                  jINC = 1
                ELSE
                  jFROM = klev
                  jTO = laytrop(i)
                  jINC = -1
                ENDIF
              ELSE
                IF (upwards) THEN
                  jFROM = laytrop(i)+1
                  jTO = klev
                  jINC = 1
                ELSE
                  jFROM = laytrop(i)-1
                  jTO = 1
                  jINC = -1
                ENDIF
              ENDIF
              DO jk = jFROM, jTO, jINC
              !NOTE: This code was not guarded against atm/band
              ! with only one major species!
                IF (atm == 1) THEN
                  jk_below = jk
                  jk_above = jk+jINC
                ELSE
                  jk_below = jk-jINC
                  jk_above = jk
                ENDIF
                IF (jp(i,jk_below) < pref_effective(band) .and. &
                  jp(i,jk_above) >= pref_effective(band)) THEN
                  IF (upwards) THEN
                    laysolfr = MIN(jk+jINC,laytrop(i)) 
                    IF (laysolfr <= laytrop(i) .and. &
                      laysolfr_ref(band) /= 0) THEN
                      laysolfr = jk
                    ENDIF
                  ELSE
                    laysolfr = MAX(jk+jINC,laytrop(i)) 
                    IF (laysolfr >= laytrop(i) .and. &
                      laysolfr_ref(band) /= 0) THEN
                      laysolfr = jk
                    ENDIF
                  ENDIF
                ENDIF
                IF (jk == laysolfr .and. pref_effective(band) >= 0) THEN
                  speccomb = gases(i,jk,major_species(1, atm, band)) 
                  specparm = speccomb + mixing_ratio(atm,band) * &
                    gases(i,jk,major_species(2,atm,band))
                  specmult = &
                    fracs_mult(atm,band) * MIN(oneminus,speccomb/specparm)
                  fs_save = specmult - INT(specmult)
                  js_save = 1  + INT(specmult)
                ENDIF
              ENDDO !jk
            ENDIF
          ENDDO !atm
          ! if an effective layer is specified for the solar flux then replace 
          ! the default solar flux at zenith with one interpolated using a value
          ! at an effective layer -- an RRTM special feature ;-)
          IF (js_save > 0 .and. nsfluxref(band) > 1) THEN
            sflx_zen(i,gpt) = flat_data(ptr+js_save) + fs_save * ( &
              flat_data(ptr+js_save+1) -flat_data(ptr+js_save)) 
          ELSE
            sflx_zen(i,gpt) = flat_data(ptr+1)
          ENDIF
        END DO !i
        ptr = ptr + sw_sfluxref(band)%s
      END DO !gpt_in
    END DO !band

  END SUBROUTINE get_solar_flux

#ifndef PSRAD_NO_INLINE
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
#endif

END MODULE mo_psrad_srtm_gas_optics
