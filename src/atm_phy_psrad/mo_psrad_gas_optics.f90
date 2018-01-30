MODULE mo_psrad_gas_optics

  USE mo_psrad_general, ONLY: wp, oneminus, &
    ngas, ih2o, ico2, ich4, io2, io3, in2o, ico, &
    nmixture, ih2oco2, ih2oo3, ih2on2o, ih2och4, in2oco2, io3co2, &
    ten20
  USE mo_psrad_flat_data, ONLY: flat_data

  IMPLICIT NONE

  PUBLIC :: spec_index_2d, get_tau_major_simple, get_tau_major_combined, &
    get_tau_minor, get_tau_gas, precomputation
#ifdef PSRAD_NO_INLINE
  PUBLIC :: spec_index_1d
#endif

  INTEGER, PARAMETER :: mixture(3,nmixture) = RESHAPE((/ &
    ih2oco2, ih2o, ico2, &
    ih2oo3, ih2o, io3, & ! Needed only in lower atmos (plog > 4.56_wp) 
    ih2on2o, ih2o, in2o, &
    ih2och4, ih2o, ich4, &
    in2oco2, in2o, ico2, &
    io3co2, io3, ico2/), & ! Needed only in upper atmos (plog <= 4.56_wp) 
    SHAPE=(/3,nmixture/))


CONTAINS

  SUBROUTINE precomputation(kproma, kbdim, klev, srtm_only, play, tlay, &
    coldry, gases, laytrop, jp, iabs, colbrd, colmol, fac, ratio, &
    h2o_factor, h2o_fraction, h2o_index, &
    minorfrac, scaleminor, scaleminorn2, indminor)

    USE mo_psrad_general, ONLY: ngas, upwards
    USE mo_psrad_setup, ONLY: jp_coeff, jp_mult, jp_fudge, log_pref, tref, &
      n_ref, tropopause_ref, stpfac, pressure_scale
    USE mo_psrad_lrtm_kgs, ONLY: chi_mls

    INTEGER, INTENT(IN) :: kproma, kbdim, klev
    LOGICAL, INTENT(IN) :: srtm_only
    REAL(wp), INTENT(IN) :: &
      play(KBDIM,klev), & ! layer pressures (mb) 
      tlay(KBDIM,klev), & ! layer temperatures (K)
      coldry(KBDIM,klev) ! dry air column density (mol/cm2)
    REAL(wp), INTENT(INOUT) :: &
      gases(KBDIM,klev,ngas) !< molecular amounts (mol/cm-2) (klev,ngas)

    INTEGER, INTENT(INOUT) :: laytrop(KBDIM), & !< tropopause layer index
      iabs(KBDIM,2,2,klev)
    INTEGER, DIMENSION(KBDIM,klev), INTENT(INOUT) :: jp, &
      indminor
    REAL(wp), DIMENSION(KBDIM,klev), INTENT(INOUT) :: colbrd, colmol, &
      minorfrac, scaleminor, scaleminorn2
    REAL(wp), INTENT(INOUT) :: fac(KBDIM,2,2,klev), &
      ratio(KBDIM,2,klev,nmixture)
    REAL(wp), DIMENSION(KBDIM,klev,2), INTENT(INOUT) :: h2o_factor,h2o_fraction
    INTEGER, DIMENSION(KBDIM,klev,2), INTENT(INOUT) :: h2o_index
    
    INTEGER :: i, j, k
    INTEGER, DIMENSION(KBDIM,klev) :: jt, jt1
    REAL(wp), DIMENSION(KBDIM,klev) :: plog, ft, ft1, water, scalefac 
    REAL(wp) :: fp

    !  Find the two reference pressures on either side of the
    !  layer pressure.  Store them in JP and JP1.  Store in FP the
    !  fraction of the difference (in ln(pressure)) between these
    !  two values that the layer pressure lies.
    DO j = 1, klev
    DO i = 1, kproma
      plog(i,j) = LOG(play(i,j))
    ENDDO
    ENDDO
    IF (upwards) THEN
      !NOTE: this branch was not properly tested
      laytrop = 0
      DO j = 1, klev
        WHERE(plog(1:kproma,j) > log_pref(tropopause_ref))
          laytrop(1:kproma) = laytrop(1:kproma) + 1
        ENDWHERE
      ENDDO
    ELSE
      laytrop = 1
      DO j = 1, klev
        WHERE(plog(1:kproma,j) <= log_pref(tropopause_ref))
          laytrop(1:kproma) = laytrop(1:kproma) + 1
        ENDWHERE
      ENDDO
    ENDIF
    DO j = 1, klev
    DO i = 1, kproma
      jp(i,j) = MIN(n_ref-1, MAX(1, &
        INT(jp_coeff - jp_mult * (plog(i,j) + jp_fudge)) ))
    ENDDO
    ENDDO
    !  Determine, for each reference pressure (JP and JP1), which
    !  reference temperature (these are different for each  
    !  reference pressure) is nearest the layer temperature but does
    !  not exceed it.  Store these indices in JT and JT1, resp.
    !  Store in FT (resp. FT1) the fraction of the way between JT
    !  (JT1) and the next highest reference temperature that the 
    !  layer temperature falls.
    DO j = 1, klev
    DO i = 1, kproma
      jt(i,j) = MIN(4,MAX(1, &
        INT(3._wp + (tlay(i,j) - tref(jp(i,j)))/15._wp) ))
      jt1(i,j) = MIN(4,MAX(1, &
        INT(3._wp + (tlay(i,j) - tref(jp(i,j)+1))/15._wp) ))
    ENDDO 
    ENDDO 
    DO j = 1, klev
    DO i = 1, kproma
      ft(i,j) = ((tlay(i,j)-tref(jp(i,j)))/15._wp) - float(jt(i,j)-3)
      ft1(i,j) = ((tlay(i,j)-tref(jp(i,j)+1))/15._wp) - float(jt1(i,j)-3)
    ENDDO 
    ENDDO 

    !  We have now isolated the layer ln pressure and temperature,
    !  between two reference pressures and two reference temperatures 
    !  (for each reference pressure).  We multiply the pressure 
    !  fraction FP with the appropriate temperature fractions to get 
    !  the factors that will be needed for the interpolation that yields
    !  the optical depths (performed in routines TAUGBn for band n).`
    DO j = 1, klev
    DO i = 1, kproma
      iabs(i,1,1,j) = MIN(11,jp(i,j)-1)*5+(jt(i,j)-1)
      iabs(i,2,1,j) = MIN(12,jp(i,j))*5+(jt1(i,j)-1)
      iabs(i,1,2,j) = MAX(0,jp(i,j)-13)*5+(jt(i,j)-1)
      iabs(i,2,2,j) = MAX(1,jp(i,j)-12)*5+(jt1(i,j)-1)
    ENDDO
    ENDDO
    DO j = 1, klev
    DO i = 1, kproma
      fp = jp_mult * (log_pref(jp(i,j)) - plog(i,j))
      fac(i,2,2,j) = fp * ft1(i,j)
      fac(i,1,2,j) = fp * (1._wp - ft1(i,j))
      fp = 1. - fp
      fac(i,2,1,j) = fp * ft(i,j)
      fac(i,1,1,j) = fp * (1._wp - ft(i,j))
    ENDDO
    ENDDO
    
    ! Tropopause defined in terms of pressure (~100 hPa)
    ! We're looking for the first layer (counted from the bottom) at which 
    ! the pressure reaches or falls below this value

    ! Calculate needed column amounts. Only a few ratios are used in the 
    ! upper atmosphere but masking may be less efficient
    ! BUG: if size mismatch cannot be determined at compile time,
    ! no exception is thrown at runtime either!
    DO k = ih2o,ico
    DO j = 1, klev
    DO i = 1, kproma
      IF (gases(i,j,k) <= 0._wp) THEN
        gases(i,j,k) = 1.e-12_wp * coldry(i,j)
      ENDIF
    ENDDO
    ENDDO
    ENDDO

    DO j = 1, klev
    DO i = 1, kproma
      colmol(i,j) = coldry(i,j) + gases(i,j,ih2o)
      water(i,j) = gases(i,j,ih2o) / coldry(i,j)
      scalefac(i,j) = play(i,j) * stpfac / tlay(i,j)
    ENDDO
    ENDDO
    !NOTE: ITS A TRAP!!! MOVED INTO SRTM_COEFFS
    ! PSRAD comment:
    !!! Water vapor continuum broadening factors are used differently 
    !!! in LW and SW? 
    !h2o_factor(:,:,1) = h2o_factor(:,:,1) * gases(:,:,ih2o)
    !h2o_factor(:,:,2) = h2o_factor(:,:,2) * gases(:,:,ih2o)
    ! Interpolation coefficients 
    DO j = 1, klev
    DO i = 1, kproma
      h2o_factor(i,j,2) = gases(i,j,ih2o) * scalefac(i,j) / (1._wp+water(i,j))
    ENDDO
    ENDDO
    ! Set up factors needed to separately include the water vapor
    ! self-continuum in the calculation of absorption coefficient.
    DO j = 1, klev
    DO i = 1, kproma
      h2o_factor(i,j,1)  = water(i,j) * h2o_factor(i,j,2)
    ENDDO
    ENDDO

    !NOTE: BEAUTIFUL CODE!!!
    !NOTE: This will break the code by allowing indexes up to 4. Previous code
    ! used uninitialized arrays. Fixed by adding ghost point in data...
    DO j = 1, klev
    DO i = 1, kproma
      fp = (332.0_wp-tlay(i,j)) / 36.0_wp
    ! PSRAD comment:
    ! If the pressure is less than ~100mb, perform a different set of species
    ! interpolations.
      IF (plog(i,j) <= log_pref(tropopause_ref)) THEN
        h2o_index(i,j,2) = 3
        !h2o_fraction(i,j,2) = (tlay(i,j)-188.0_wp)/36.0_wp - 1.0_wp
        h2o_fraction(i,j,2) = (tlay(i,j)-224.0_wp)/36.0_wp
      ELSE
        h2o_index(i,j,2) = MIN(2, MAX(1, INT(fp)))
        h2o_fraction(i,j,2) = fp - FLOAT(h2o_index(i,j,2))
      ENDIF
    ENDDO
    ENDDO

    ! PSRAD comment:
    ! In RRTMG code, this calculation is done only in the lower atmosphere 
    ! (plog > 4.56) 
    DO j = 1, klev
    DO i = 1, kproma
      fp = (tlay(i,j)-188.0_wp)/7.2_wp
      h2o_index(i,j,1) = MIN(9, MAX(1, INT(fp)-7))
      h2o_fraction(i,j,1) = fp - FLOAT(h2o_index(i,j,1) + 7)
    ENDDO
    ENDDO

    IF (srtm_only) RETURN

    ! PSRAD comment:
    ! Broadening gases -- the number of molecules per cm^2 of all gases 
    ! not specified explicitly (water is excluded) 
    DO j = 1, klev
    DO i = 1, kproma
      colbrd(i,j) = coldry(i,j) - &
        SUM(gases(i,j,(/ico2,io3,in2o,ich4,io2/)))
    ENDDO
    ENDDO
    ! PSRAD comment:
    !  Set up factors needed to separately include the minor gases
    !  in the calculation of absorption coefficient
    DO j = 1, klev
    DO i = 1, kproma
      scaleminor(i,j) = play(i,j) / (pressure_scale * tlay(i,j))
      scaleminorn2(i,j) = colbrd(i,j)**2 * scaleminor(i,j) / colmol(i,j)
    ENDDO
    ENDDO
    DO j = 1, klev
    DO i = 1, kproma
      fp = (tlay(i,j)-180.8_wp) / 7.2_wp
      indminor(i,j) = MIN(18, MAX(1, INT(fp)))
      minorfrac(i,j) = fp - FLOAT(indminor(i,j))
    ENDDO
    ENDDO
    
    ! PSRAD comment:
    ! Setup reference ratio to be used in calculation of binary species 
    ! parameter.
    ratio = ten20
    DO k = 1,nmixture
      DO j = 1,klev
      DO i = 1, kproma
        IF (chi_mls(mixture(3,k), jp(i,j)) /= 0) THEN
          ratio(i,1,j,mixture(1,k)) = chi_mls(mixture(2,k), jp(i,j)) /&
            chi_mls(mixture(3,k), jp(i,j))
        ENDIF
        IF (chi_mls(mixture(3,k), jp(i,j)+1) /= 0) THEN
          ratio(i,2,j,mixture(1,k)) = chi_mls(mixture(2,k), jp(i,j)+1) /&
            chi_mls(mixture(3,k), jp(i,j)+1)
        ENDIF
      ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE precomputation

#ifdef PSRAD_NO_INLINE
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

  SUBROUTINE spec_index_2d(kproma, kbdim, gas1, ratio, gas2, m, iabs, nsp, &
    parm, comb, js, fs)
    INTEGER, INTENT(IN) :: kproma, kbdim, iabs(KBDIM), nsp, m
    REAL(wp), INTENT(IN) :: gas1(KBDIM), ratio(KBDIM), gas2(KBDIM)
    INTEGER, INTENT(INOUT) :: js(KBDIM)
    REAL(wp), INTENT(INOUT) :: parm(KBDIM), comb(KBDIM), fs(KBDIM)
    REAL(wp) :: mult(KBDIM)

    INTEGER :: i

    DO i = 1, kproma
      comb(i) = gas1(i) + ratio(i) * gas2(i)
      parm(i) = MIN(oneminus, gas1(i)/comb(i))
      mult(i) = m * parm(i)
      js(i) = INT(mult(i)) 
      fs(i) = mult(i) - js(i)
      js(i) = js(i) + 1 + iabs(i) * nsp 
    ENDDO
  END SUBROUTINE spec_index_2d

  SUBROUTINE get_tau_major_combined(kproma, kbdim, klev, atm, &
    range, gpt_range, m, nsp, gas1, gas2, &
    constant_ratio, ratio_index, variable_ratio, &
    fac, iabs, ref_o, ref_s, fancy, tau, o)
    INTEGER, INTENT(IN) :: kproma, kbdim, klev, atm, range(2), &
      gpt_range(2), iabs(KBDIM,2,2,klev), ref_o, ref_s, o
    REAL(wp), INTENT(IN) :: gas1(KBDIM,klev), constant_ratio, &
      gas2(KBDIM,klev), fac(KBDIM,2,2,klev)
    INTEGER, INTENT(IN) :: m, nsp, ratio_index
    LOGICAL, INTENT(IN) :: fancy
    REAL(wp), TARGET, INTENT(IN) :: variable_ratio(:,:,:,:)
    REAL(wp), INTENT(INOUT) :: tau(:,:,:)
    REAL(wp), POINTER :: current_ratio(:)

    REAL(wp), DIMENSION(KBDIM) :: speccomb, specparm, fs
    REAL(wp) :: fk, acc, p, p4
    REAL(wp), TARGET :: array_ratio(KBDIM)
    INTEGER :: js(KBDIM), j, gpt, gpt_in, i, pass, lay

    !NOTE: These are derived from nsp(atm,band) being always 9 (+- 1...)
    INTEGER, PARAMETER :: stencil_6(6,2) = RESHAPE((/ &
      0, 9, 1, 10, 2, 11, &
      1, 10, 0, 9, -1, 8/), SHAPE =(/6,2/))

    INTEGER :: mask(KBDIM)
    LOGICAL :: any_non_fancy
    REAL(wp) :: tau_inc(KBDIM)

    IF (ratio_index == 0) THEN
      current_ratio => array_ratio
      array_ratio = constant_ratio
    ENDIF
    tau(1:kproma,range(1)+o:range(2)+o,gpt_range(1):gpt_range(2)) = 0.0_wp
    mask = 0
    any_non_fancy = .true.
    DO lay = range(1),range(2)
      DO pass = 1,2
        IF (ratio_index /= 0) THEN
          current_ratio => variable_ratio(:,pass,lay,ratio_index)
        ENDIF
        CALL spec_index_2d(kproma, kbdim, gas1(:,lay), &
          current_ratio, gas2(:,lay), &
          m, iabs(:,pass,atm,lay), nsp, &
          specparm, speccomb, js, fs)
        gpt_in = 0
        js(1:kproma) = js(1:kproma) + ref_o
        DO gpt = gpt_range(1),gpt_range(2)
          gpt_in = gpt_in+1
          IF (fancy) THEN
            mask = 0
            WHERE (specparm(1:kproma) < 0.125_wp) mask(1:kproma) = -1
            WHERE (specparm(1:kproma) > 0.875_wp) mask(1:kproma) = 1
            any_non_fancy = ANY(mask(1:kproma) == 0)
          ENDIF
          IF (any_non_fancy) THEN
            DO i = 1,kproma
              tau_inc(i) = speccomb(i) * &
                (fs(i) * ( &
                  fac(i,1,pass,lay) * flat_data(js(i) + 1) + &
                  fac(i,2,pass,lay) * flat_data(js(i) + nsp+1)) + &
                (1._wp - fs(i)) * ( &
                  fac(i,1,pass,lay) * flat_data(js(i) + 0) + &
                  fac(i,2,pass,lay) * flat_data(js(i) + nsp)))
            ENDDO
          ENDIF
          IF (fancy .and. ANY(mask(1:kproma) /= 0)) THEN
            DO i = 1,kproma
              IF (mask(i) == 0) CYCLE
              IF (mask(i) < 0) THEN
                j = 1
                p = fs(i) - 1
              ELSE
                j = 2
                p = -fs(i)
              ENDIF
              p4 = p**4
              acc = p4 * ( &
                fac(i,1,pass,lay) * flat_data(js(i) + stencil_6(1,j)) + &
                fac(i,2,pass,lay) * flat_data(js(i) + stencil_6(2,j)))
              fk = 1 - p - 2.0_wp*p4
              acc = acc + fk * ( &
                fac(i,1,pass,lay) * flat_data(js(i) + stencil_6(3,j)) + &
                fac(i,2,pass,lay) * flat_data(js(i) + stencil_6(4,j)))
              fk = p + p4
              acc = acc + fk * ( &
                fac(i,1,pass,lay) * flat_data(js(i) + stencil_6(5,j)) + &
                fac(i,2,pass,lay) * flat_data(js(i) + stencil_6(6,j)))
              tau_inc(i) = speccomb(i) * acc
            ENDDO
          ENDIF
          tau(1:kproma,lay+o,gpt) = &
            tau(1:kproma,lay+o,gpt) + tau_inc(1:kproma)
          js(1:kproma) = js(1:kproma) + ref_s
        ENDDO
      ENDDO
    ENDDO
  END SUBROUTINE get_tau_major_combined

  SUBROUTINE get_tau_major_simple(kproma, kbdim, klev, atm, band, &
    gpt_range, gas, fac, range, iabs, ref_o, ref_s, tau, o)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: kproma, kbdim, klev, atm, band, range(2), &
      iabs(KBDIM,2,2,klev), gpt_range(2), ref_o, ref_s, o
    REAL(wp), INTENT(IN) :: gas(KBDIM,klev), &
      fac(KBDIM,2,2,klev)
    REAL(wp), INTENT(INOUT) :: tau(:,:,:)
    INTEGER :: gpt, gpt_in, lay, i, idx1, idx2, ptr

    ! The two lines marked below should read
    !ind = iabs(:,1,atm,lay) * nsp(atm,band) + 1
    !ind = iabs(:,2,atm,lay) * nsp(atm,band) + 1
    ! but nsp(atm,band) is always 1 when this function is called.
    !... except for the case of broken band 16, which is handled below
    ! for backward compatibility
    ! duplicated.
    INTEGER :: m
    m = 1
    !FIXME: data for band 16 atm 2 is not homogeneous!
    IF (atm == 2 .and. band == 16) m = 0 
    gpt_in = 0
    ptr = ref_o
    DO gpt = gpt_range(1),gpt_range(2)
      gpt_in = gpt_in+1
      DO lay = range(1),range(2)
        DO i = 1,kproma
          idx1 = iabs(i,1,atm,lay) * m + 1
          idx2 = iabs(i,2,atm,lay) * m + 1
          tau(i,lay+o,gpt) = gas(i,lay) * ( &
            fac(i,1,1,lay) * flat_data(ptr+idx1) + &
            fac(i,2,1,lay) * flat_data(ptr+idx1+1) + &
            fac(i,1,2,lay) * flat_data(ptr+idx2) + &
            fac(i,2,2,lay) * flat_data(ptr+idx2+1))
        ENDDO
      ENDDO
      ptr = ptr + ref_s
    ENDDO
  END SUBROUTINE get_tau_major_simple

  SUBROUTINE get_tau_minor(kproma, kbdim, klev, laytrop, range, atm, &
    gpt_range, factor, fraction, index, ref_o, ref_s, tau, o, save_target)
    INTEGER, INTENT(IN) :: kproma, kbdim, klev, laytrop(KBDIM), &
      range(2), atm, gpt_range(2), ref_o, ref_s, o, save_target
    REAL(wp), DIMENSION(KBDIM,klev), INTENT(IN) :: factor, fraction
    INTEGER, DIMENSION(KBDIM,klev), INTENT(IN) :: index
    REAL(wp), INTENT(INOUT) :: tau(:,:,:)
    INTEGER :: lay, gpt_in, gpt, i, idx, ptr
#ifdef PSRAD_DEVEL
    REAL(wp) :: tau_inc(KBDIM)
#endif

    ptr = ref_o
    gpt_in = 0
    DO gpt = gpt_range(1),gpt_range(2)
      gpt_in = gpt_in+1
      DO lay = range(1),range(2)
        ! This loop is explicit to work around bug in NAG 6.1 Build 6106 
        DO i = 1,kproma
          idx = index(i,lay)
#ifdef PSRAD_DEVEL
          tau_inc(i) = &
#else
          tau(i,lay+o,gpt) = tau(i,lay+o,gpt) + &
#endif
            factor(i,lay) * (&
              flat_data(ptr+idx) + &
              fraction(i,lay) * &
              (flat_data(ptr+idx+1) - &
              flat_data(ptr+idx)))
#ifdef PSRAD_DEVEL
          tau(i,lay+o,gpt) = tau(i,lay+o,gpt) + tau_inc(i)
#endif
        ENDDO
#ifdef PSRAD_DEVEL
        CALL save_tau_inc(kproma, KBDIM, laytrop, lay, atm, gpt, save_target, &
          tau_inc)
#endif
      ENDDO
      ptr = ptr + ref_s
    ENDDO
  END SUBROUTINE get_tau_minor

  SUBROUTINE get_tau_gas(kproma, kbdim, laytrop, range, atm, &
    gpt_range, gas, ref_o, tau, o, save_target)

    INTEGER, INTENT(IN) :: kproma, kbdim, laytrop(KBDIM), range(2), atm, &
      gpt_range(2), ref_o, o, save_target
    REAL(wp), INTENT(IN) :: gas(:,:)
    REAL(wp), INTENT(INOUT) :: tau(:,:,:)
    INTEGER :: gpt, gpt_in, i, lay, ptr
#ifdef PSRAD_DEVEL
    REAL(wp) :: tau_inc(KBDIM)
#endif

    gpt_in = 0
    DO gpt = gpt_range(1),gpt_range(2)
      gpt_in = gpt_in+1
      ptr = ref_o + gpt_in
      DO lay = range(1),range(2)
        DO i = 1,kproma
#ifdef PSRAD_DEVEL
          tau_inc(i) = gas(i,lay) * flat_data(ptr)
          tau(i,lay+o,gpt) = tau(i,lay+o,gpt) + tau_inc(i)
#else
          tau(i,lay+o,gpt) = tau(i,lay+o,gpt) + &
            gas(i,lay) * flat_data(ptr)
#endif
        ENDDO
#ifdef PSRAD_DEVEL
        CALL save_tau_inc(kproma, KBDIM, laytrop, lay, atm, gpt, save_target, &
          tau_inc)
#endif
      ENDDO
    ENDDO
  END SUBROUTINE get_tau_gas

END MODULE mo_psrad_gas_optics
