MODULE mo_psrad_gas_optics

  USE mo_psrad_general, ONLY : wp, oneminus, &
    ngas, ih2o, ico2, ich4, io2, io3, in2o, ico, &
    nmixture, ih2oco2, ih2oo3, ih2on2o, ih2och4, in2oco2, io3co2, &
    ten20, ten20inv

  PUBLIC :: spec_index_2d, get_tau_major_simple, get_tau_major_combined, &
    get_tau_minor, get_tau_gas, precomputation

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
    coldry, wkl, laytrop, jp, iabs, gases, colbrd, colmol, fac, ratio, &
    h2o_factor, h2o_fraction, h2o_index, &
    minorfrac, scaleminor, scaleminorn2, indminor)

    USE mo_psrad_general, ONLY: ngas, upwards
    USE mo_psrad_setup, ONLY : jp_coeff, jp_mult, jp_fudge, log_pref, tref, &
      n_ref, tropopause_ref, stpfac, pressure_scale
    USE mo_psrad_lrtm_kgs, ONLY : chi_mls

    INTEGER, INTENT(IN) :: kproma, kbdim, klev
    LOGICAL, INTENT(IN) :: srtm_only
    REAL(wp), INTENT(IN) :: &
      play(kbdim,klev), & ! layer pressures (mb) 
      tlay(kbdim,klev), & ! layer temperatures (K)
      coldry(kbdim,klev), & ! dry air column density (mol/cm2)
      wkl(kbdim,klev,ngas) !< molecular amounts (mol/cm-2) (klev,ngas)

    INTEGER, INTENT(INOUT) :: laytrop(kbdim), & !< tropopause layer index
      iabs(kbdim,2,2,klev)
    INTEGER, DIMENSION(kbdim,klev), INTENT(INOUT) :: jp, &
      indminor
    REAL(wp), DIMENSION(kbdim,klev), INTENT(INOUT) :: colbrd, colmol, &
      minorfrac, scaleminor, scaleminorn2
    REAL(wp), INTENT(INOUT) :: fac(kbdim,2,2,klev), gases(kbdim,klev,ngas), &
      ratio(kbdim,2,klev,nmixture)
    REAL(wp), DIMENSION(kbdim,klev,2), INTENT(INOUT) :: h2o_factor,h2o_fraction
    INTEGER, DIMENSION(kbdim,klev,2), INTENT(INOUT) :: h2o_index
    
    INTEGER :: i, j, k
    INTEGER, DIMENSION(kbdim,klev) :: jt, jt1
    REAL(wp), DIMENSION(kbdim,klev) :: plog, ft, ft1, water, scalefac, wbroad
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
      DO i = 1, kproma
        laytrop(i) = COUNT(plog(i,:) > log_pref(tropopause_ref))
      ENDDO
    ELSE
      DO i = 1, kproma
        laytrop(i) = COUNT(plog(i,:) <= log_pref(tropopause_ref)) + 1
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
    gases(1:kproma,:,ih2o:ico) = ten20inv * wkl(1:kproma,:,ih2o:ico)
    DO k = ih2o,ico
    DO j = 1, klev
    DO i = 1, kproma
      IF (wkl(i,j,k) <= 0._wp) THEN
        gases(i,j,k) = 1.e-32_wp * coldry(i,j)
      ENDIF
    ENDDO
    ENDDO
    ENDDO
    colmol(1:kproma,:) = ten20inv * coldry(1:kproma,:) + gases(1:kproma,:,ih2o)

    water(1:kproma,:) = wkl(1:kproma,:,ih2o)/coldry(1:kproma,:)
    scalefac(1:kproma,:) = play(1:kproma,:) * stpfac / tlay(1:kproma,:)
    ! FIXME?
    ! Original comment:
    !!! Water vapor continuum broadening factors are used differently 
    !!! in LW and SW? 
    !h2o_factor(:,:,1) = h2o_factor(:,:,1) * gases(:,:,ih2o)
    !h2o_factor(:,:,2) = h2o_factor(:,:,2) * gases(:,:,ih2o)
    ! Interpolation coefficients 
    h2o_factor(1:kproma,:,2) = &
      gases(1:kproma,:,ih2o) * scalefac(1:kproma,:) / &
      (1._wp+water(1:kproma,:))
    ! Set up factors needed to separately include the water vapor
    ! self-continuum in the calculation of absorption coefficient.
    h2o_factor(1:kproma,:,1)  = water(1:kproma,:) * h2o_factor(1:kproma,:,2)

    ! FIXME? Is the data actually of size 3?
    ! Original comment:
    !NOTE: This will break the code by allowing indexes up to 4. Previous code
    ! used uninitialized arrays. Fixed by adding ghost point in data...
    DO j = 1, klev
    DO i = 1, kproma
      fp = (332.0_wp-tlay(i,j)) / 36.0_wp
    ! If the pressure is less than ~100mb, perform a different set of species
    ! interpolations.
      IF (plog(i,j) <= log_pref(tropopause_ref)) THEN
        h2o_index(i,j,2) = 3
        h2o_fraction(i,j,2) = (tlay(i,j)-188.0_wp)/36.0_wp - 1.0_wp
      ELSE
        h2o_index(i,j,2) = MIN(2, MAX(1, INT(fp)))
        h2o_fraction(i,j,2) = fp - FLOAT(h2o_index(i,j,2))
      ENDIF
    ENDDO
    ENDDO

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

    ! Broadening gases -- the number of molecules per cm^2 of all gases 
    ! not specified explicitly (water is excluded) 
    DO j = 1, klev
    DO i = 1, kproma
      wbroad(i,j) = coldry(i,j) - &
        SUM(wkl(i,j,(/ico2,io3,in2o,ich4,io2/)))
      colbrd(i,j) = ten20inv * wbroad(i,j)
    ENDDO
    ENDDO
    !  Set up factors needed to separately include the minor gases
    !  in the calculation of absorption coefficient
    DO j = 1, klev
    DO i = 1, kproma
      scaleminor(i,j) = play(i,j) / (pressure_scale * tlay(i,j))
      scaleminorn2(i,j) = colbrd(i,j) * scaleminor(i,j) * &
        (wbroad(i,j) / (coldry(i,j)+wkl(i,j,ih2o)))
    ENDDO
    ENDDO
    DO j = 1, klev
    DO i = 1, kproma
      fp = (tlay(i,j)-180.8_wp) / 7.2_wp
      indminor(i,j) = MIN(18, MAX(1, INT(fp)))
      minorfrac(i,j) = fp - FLOAT(indminor(i,j))
    ENDDO
    ENDDO
    
    !  Setup reference ratio to be used in calculation of binary species 
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

  SUBROUTINE spec_index_2d(kproma, gas1, ratio, gas2, m, iabs, nsp, &
    parm, comb, js, fs)
    INTEGER, INTENT(IN) :: kproma, iabs(:), nsp, m
    REAL(wp), INTENT(IN) :: gas1(:), ratio(:), gas2(:)
    INTEGER, INTENT(INOUT) :: js(:)
    REAL(wp), INTENT(INOUT) :: parm(:), comb(:), fs(:)

    REAL(wp) :: mult
    INTEGER :: i

    DO i = 1, kproma
      comb(i) = gas1(i) + ratio(i) * gas2(i)
      parm(i) = MIN(oneminus, gas1(i)/comb(i))
      mult = m * parm(i)
      js(i) = INT(mult) 
      fs(i) = mult - js(i)
      js(i) = js(i) + 1 + iabs(i) * nsp 
    ENDDO
  END SUBROUTINE spec_index_2d

  SUBROUTINE get_tau_major_combined(kproma, kbdim, atm, &
    range, gpt_range, m, nsp, gas1, gas2, &
    constant_ratio, ratio_index, variable_ratio, &
    fac, iabs, ref, fancy, tau, o)
    INTEGER, INTENT(IN) :: kproma, kbdim, atm, range(2), &
      gpt_range(2), iabs(:,:,:,:), o
    REAL(wp), INTENT(IN) :: gas1(:,:), constant_ratio, &
      gas2(:,:), fac(:,:,:,:), ref(:,:)
    INTEGER, INTENT(IN) :: m, nsp, ratio_index
    LOGICAL, INTENT(IN) :: fancy
    REAL(wp), POINTER, INTENT(IN) :: variable_ratio(:,:,:,:)
    REAL(wp), INTENT(INOUT) :: tau(:,:,:)
    REAL(wp), POINTER :: current_ratio(:)

    REAL(wp), DIMENSION(kbdim) :: speccomb, specparm, fs
    REAL(wp) :: fk, acc, p, p4
    REAL(wp), TARGET :: array_ratio(kbdim)
    INTEGER :: js(kbdim), j, gpt, gpt_in, i, pass, lay

    ! These are derived from nsp(atm,band) being always 9 (+- 1...)
    INTEGER, PARAMETER :: stencil_6(6,2) = RESHAPE((/ &
      0, 9, 1, 10, 2, 11, &
      1, 10, 0, 9, -1, 8/), SHAPE =(/6,2/))

    INTEGER :: mask(kbdim)
    LOGICAL :: any_non_fancy
    REAL(wp) :: tau_inc(kbdim)

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
        CALL spec_index_2d(kproma, gas1(:,lay), &
          current_ratio, gas2(:,lay), &
          m, iabs(:,pass,atm,lay), nsp, &
          specparm, speccomb, js, fs)
        gpt_in = 0
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
                  fac(i,1,pass,lay) * ref(js(i) + 1,gpt_in) + &
                  fac(i,2,pass,lay) * ref(js(i) + nsp+1,gpt_in)) + &
                (1._wp - fs(i)) * ( &
                  fac(i,1,pass,lay) * ref(js(i) + 0,gpt_in) + &
                  fac(i,2,pass,lay) * ref(js(i) + nsp,gpt_in)))
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
                fac(i,1,pass,lay) * ref(js(i) + stencil_6(1,j),gpt_in) + &
                fac(i,2,pass,lay) * ref(js(i) + stencil_6(2,j),gpt_in))
              fk = 1 - p - 2.0_wp*p4
              acc = acc + fk * ( &
                fac(i,1,pass,lay) * ref(js(i) + stencil_6(3,j),gpt_in) + &
                fac(i,2,pass,lay) * ref(js(i) + stencil_6(4,j),gpt_in))
              fk = p + p4
              acc = acc + fk * ( &
                fac(i,1,pass,lay) * ref(js(i) + stencil_6(5,j),gpt_in) + &
                fac(i,2,pass,lay) * ref(js(i) + stencil_6(6,j),gpt_in))
              tau_inc(i) = speccomb(i) * acc
            ENDDO
          ENDIF
          tau(1:kproma,lay+o,gpt) = &
            tau(1:kproma,lay+o,gpt) + tau_inc(1:kproma)
        ENDDO
      ENDDO
    ENDDO
  END SUBROUTINE get_tau_major_combined

  SUBROUTINE get_tau_major_simple(kproma, atm, band, &
    gpt_range, gas, fac, range, iabs, ref, tau, o)
    INTEGER, INTENT(IN) :: kproma, atm, band, range(2), &
      iabs(:,:,:,:), gpt_range(2), o
    REAL(wp), INTENT(IN) :: gas(:,:), &
      fac(:,:,:,:), ref(:,:)
    REAL(wp), INTENT(INOUT) :: tau(:,:,:)
    INTEGER :: gpt, gpt_in, lay, i, idx1, idx2 

    ! The two lines marked below should read
    !ind = iabs(:,1,atm,lay) * nsp(atm,band) + 1
    !ind = iabs(:,2,atm,lay) * nsp(atm,band) + 1
    ! but nsp(atm,band) is always 1 when this function is called.
    !... except for the case of broken band 16, which is handled below
    ! for backward compatibility
    ! duplicated.
    INTEGER :: m

    m = 1
    IF (atm == 2 .and. band == 16) m = 0 
    gpt_in = 0
    DO gpt = gpt_range(1),gpt_range(2)
      gpt_in = gpt_in+1
      DO lay = range(1),range(2)
        DO i = 1,kproma
          idx1 = iabs(i,1,atm,lay) * m + 1
          idx2 = iabs(i,2,atm,lay) * m + 1
          tau(i,lay+o,gpt) = gas(i,lay) * ( &
            fac(i,1,1,lay) * ref(idx1,gpt_in) + &
            fac(i,2,1,lay) * ref(idx1+1,gpt_in) + &
            fac(i,1,2,lay) * ref(idx2,gpt_in) + &
            fac(i,2,2,lay) * ref(idx2+1,gpt_in))
        ENDDO
      ENDDO
    ENDDO
  END SUBROUTINE get_tau_major_simple

  SUBROUTINE get_tau_minor(kproma, range, gpt_range, &
    factor, fraction, index, ref, tau, o)
    INTEGER, INTENT(IN) :: kproma, range(2), gpt_range(2), o
    REAL(wp), DIMENSION(:,:), INTENT(IN) :: factor, fraction
    REAL(wp), DIMENSION(:,:), INTENT(IN) :: ref
    INTEGER, DIMENSION(:,:), INTENT(IN) :: index
    REAL(wp), INTENT(INOUT) :: tau(:,:,:)
    INTEGER :: lay, gpt_in, gpt, i, idx

    gpt_in = 0
    DO gpt = gpt_range(1),gpt_range(2)
      gpt_in = gpt_in+1
      DO lay = range(1),range(2)
        ! This loop is explicit to work around bug in NAG 6.1 Build 6106 
        DO i = 1,kproma
          idx = index(i,lay)
          tau(i,lay+o,gpt) = tau(i,lay+o,gpt) + &
            factor(i,lay) * (&
              ref(idx,gpt_in) + &
              fraction(i,lay) * &
              (ref(idx+1,gpt_in) - &
              ref(idx,gpt_in)))
        ENDDO
      ENDDO
    ENDDO
  END SUBROUTINE get_tau_minor

  SUBROUTINE get_tau_gas(kproma, range, &
    gpt_range, gas, ref, tau, o)

    INTEGER, INTENT(IN) :: kproma, range(2), &
      gpt_range(2), o
    REAL(wp), INTENT(IN) :: gas(:,:), ref(:)
    REAL(wp), INTENT(INOUT) :: tau(:,:,:)
    INTEGER :: gpt, gpt_in, i, lay

    gpt_in = 0
    DO gpt = gpt_range(1),gpt_range(2)
      gpt_in = gpt_in+1
      DO lay = range(1),range(2)
        DO i = 1,kproma
          tau(i,lay+o,gpt) = tau(i,lay+o,gpt) + &
            gas(i,lay) * ref(gpt_in)
        ENDDO
      ENDDO
    ENDDO
  END SUBROUTINE get_tau_gas

END MODULE mo_psrad_gas_optics
