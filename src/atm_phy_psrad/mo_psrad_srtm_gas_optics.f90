!> !! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_psrad_srtm_gas_optics

  USE mo_psrad_general, ONLY: wp, nbndsw, ngptsw
  USE mo_psrad_srtm_setup,  ONLY: ngb, ngs, nsp
  USE mo_psrad_srtm_kgs

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: gpt_taumol, ih2o, ich4, ico2, io2, io3

  LOGICAL :: initialized = .FALSE.
  REAL(wp), PARAMETER :: oneminus = 1.0_wp - 1.0e-6_wp
  REAL(wp) :: &
    rayl_24_js(9,ngptsw), &
    abs_ab(1175,2,ngptsw) , &
    xk_mnr(2,ngptsw) , &
    ref_self(10,ngptsw) , &
    ref_frgn(4,ngptsw) , &
    ref_sflx(9,ngptsw) , &
    s_wght1(2,nbndsw)   , &
    s_wght2(2,nbndsw)   , &
    rayl(ngptsw)
  INTEGER  :: igas(4,2,nbndsw), ioff(2,nbndsw), &
    layreffr(nbndsw), laysolfr_ref(nbndsw)
! DEAD CODE/BUG??? What were these for???
!  REAL(wp) :: s_wght0(2,nbndsw)
!  INTEGER  :: laysolfr

  INTEGER, PARAMETER :: nogas = 0, ih2o  = 1, ico2  = 2, ich4  = 3, &
    io2   = 4, io3   = 5

CONTAINS

  SUBROUTINE gpt_taumol(jcs, kproma,kbdim, klev, jp, jt, jt1, laytrop, indself, indfor, &
    gases, colmol, fac, selffac, selffrac, forfac, &
    forfrac, sflx_zen, taug, taur)

    INTEGER, INTENT(IN) :: jcs, kproma, kbdim, klev, jp(:,:), jt(:,:), jt1(:,:), laytrop(:), &
      indself(:,:), indfor(:,:)

    REAL(wp), INTENT(IN) :: colmol(:,:), gases(:,:,:), &
      selffac(:,:), selffrac(:,:), forfac(:,:), forfrac(:,:)
    REAL(wp), INTENT(IN) :: fac(kbdim,2,2,klev)

    REAL(wp), INTENT(INOUT) :: sflx_zen(:,:), taug(:,:,:), taur(:,:,:)

    INTEGER, PARAMETER :: noffset(2) = (/1,61/)
    INTEGER  :: jl, ig, jk, js, ii, ib, itrop, jk1, jk2, laysolfr, &
      i000, i001, i010, i011, i100, i101, i110, i111
    REAL(wp) :: fs, gas_mjr, specmult, gas_key, fs_wght, &
      x000, x001, x010, x011, x100, x101, x110, x111

    IF (.NOT.initialized) CALL setup_taumol(klev)

    DO jl = jcs,kproma
    DO ig = 1,ngptsw
      js = 0
      ib = ngb(ig)
      ! IF (laysolfr_ref(ib) == 0) laysolfr = laytrop
      laysolfr = laytrop(jl)
      DO jk = 1,klev
        itrop = MAX(0,MIN(1,jk-laytrop(jl))) + 1    ! itrop = 1 below 100 hPa, 2 above

        gas_key = 0.0_wp
        IF (igas(1,itrop,ib) > 0) gas_key = gases(jl,jk,igas(1, itrop, ib)) 
        gas_mjr = gas_key

        fs = 0.0_wp
        ii = (jp(jl,jk)*5 - noffset(itrop) )*nsp(itrop,ib) + 1 
        IF (igas(2,itrop,ib) > 0) THEN
          gas_mjr = gas_key + s_wght1(itrop,ib)*gases(jl,jk,igas(2,itrop,ib))
          specmult = s_wght2(itrop,ib) * MIN(oneminus,gas_key/gas_mjr)
          fs   = specmult - AINT(specmult)
          ii   = ii + INT(specmult)
        END IF

        i000 = ii + (jt(jl,jk) - 5)*nsp(itrop,ib)
        i001 = ii + (jt1(jl,jk))*nsp(itrop,ib)
        i100 = i000+1
        i101 = i001+1
        i110 = i100+ioff(itrop,ib)
        i010 = i000+ioff(itrop,ib)
        i111 = i101+ioff(itrop,ib)
        i011 = i001+ioff(itrop,ib)

        x100 = fac(jl,1,1,jk) * fs
        x000 = fac(jl,1,1,jk) - x100
        x101 = fac(jl,1,2,jk) * fs
        x001 = fac(jl,1,2,jk) - x101
        x110 = fac(jl,2,1,jk) * fs
        x010 = fac(jl,2,1,jk) - x110
        x111 = fac(jl,2,2,jk) * fs
        x011 = fac(jl,2,2,jk) - x111

        taug(jl,jk,ig) = 0.0_wp

        IF (igas(1,itrop,ib) > 0) THEN
          taug(jl,jk,ig) = taug(jl,jk,ig) + gas_mjr * ( &
            x000 * abs_ab(i000,itrop,ig) + &
            x001 * abs_ab(i001,itrop,ig) + &
            x010 * abs_ab(i010,itrop,ig) + &
            x011 * abs_ab(i011,itrop,ig) )
        END IF

        IF (igas(2,itrop,ib) > 0) THEN
          taug(jl,jk,ig) = taug(jl,jk,ig) + gas_mjr * ( &
             x100 * abs_ab(i100,itrop,ig) + &
             x101 * abs_ab(i101,itrop,ig) + &
             x110 * abs_ab(i110,itrop,ig) + &
             x111 * abs_ab(i111,itrop,ig) )
        END IF

        IF (igas(3,itrop,ib) > 0) THEN
          taug(jl,jk,ig) = taug(jl,jk,ig) + gases(jl,jk,igas(3, itrop, ib)) * ( &
            forfac(jl,jk) * (ref_frgn(indfor(jl,jk),ig) + forfrac(jl,jk) * &
            (ref_frgn(indfor(jl,jk)+1,ig) - ref_frgn(indfor(jl,jk),ig))))
        END IF

        IF (igas(3,itrop,ib) > 0 .AND. selffac(jl,jk) > EPSILON(1.0_wp)) THEN
          taug(jl,jk,ig) = taug(jl,jk,ig) + gases(jl,jk,igas(3, itrop, ib)) * ( &
            selffac(jl,jk) * (ref_self(indself(jl,jk),ig) + selffrac(jl,jk) * &
            (ref_self(indself(jl,jk)+1,ig) - ref_self(indself(jl,jk),ig))))
        END IF

        IF (igas(4,itrop,ib) > 0) THEN
          taug(jl,jk,ig) = taug(jl,jk,ig) + gases(jl,jk,igas(4, itrop, ib)) * xk_mnr(itrop,ig) 
        END IF

        taur(jl,jk,ig) = colmol(jl,jk)*rayl(ig)

        IF (ib == 9 .AND. itrop == 1) THEN
          js = 1 + INT(specmult)
          taur(jl,jk,ig) = colmol(jl,jk) * (rayl_24_js(js,ig) + &
               fs * (rayl_24_js(js+1,ig) - rayl_24_js(js,ig)))
          js = 0
        END IF

        jk1 = jk+1-itrop
        jk2 = jk1+1
        IF (jp(jl,jk1) < layreffr(ib) .AND. jp(jl,jk2) >= layreffr(ib))               &
             laysolfr = MIN(jk+1,laytrop(jl))
        IF (laysolfr <=  laytrop(jl) .and. laysolfr_ref(ib) > laytrop(jl))              &
             laysolfr = jk

        IF (jk == laysolfr .AND. layreffr(ib) >= 0) THEN
          js      = 1  + INT(specmult)
          fs_wght = fs
        END IF
      ENDDO
      !
      ! if an effective layer is specified for the solar flux then replace 
      ! the default solar flux at zenith with one interpolated using a value
      ! at an effective layer -- an RRTM special feature ;-)
      !
      sflx_zen(jl,ig) = ref_sflx(1,ig)
      IF (js > 0) sflx_zen(jl,ig) =                                             &
           ref_sflx(js,ig) + fs_wght*(ref_sflx(js+1,ig) -ref_sflx(js,ig)) 
    END DO
    END DO

  END SUBROUTINE gpt_taumol

  SUBROUTINE setup_taumol(klev)

    INTEGER, INTENT (IN) :: klev

    INTEGER :: ig, jg, ibeg, iband, ib
    !  Prepare band 16:  2600-3250 cm-1 (low - h2o,ch4; high - ch4)
    ibeg  = 1
    iband = ngb(ibeg)
    ib = iband 
    ioff(:,ib) = (/9,1/)
    igas(:,:,ib)  = RESHAPE ( SOURCE = (/ih2o, ich4, ih2o, nogas, &
      ich4, nogas, nogas, nogas/), SHAPE = (/4,2/) )
    s_wght1(:,ib) = (/252.131_wp, 0.0_wp/)
    s_wght2(:,ib) = (/8.0_wp, 8.0_wp/)
    layreffr(ib) = -18
    laysolfr_ref(ib) = klev

    DO jg=1,ng(1)
      ig = jg + ibeg - 1
      rayl(ig) = rayl16
      abs_ab(1:585,1,ig) = absa16(:,jg)
      abs_ab(1:235,2,ig) = absb16(:,jg)
      ref_self(1:10,ig) = selfref16(:,jg)
      ref_frgn(1:3 ,ig) = forref16(:,jg)
      ref_sflx(1,ig) = sfluxref16(jg)
    END DO

    !  Prepare band 17:  3250-4000 cm-1 (low - h2o,co2; high - h2o,co2)
    ibeg  = ngs(ib) + 1
    iband = ngb(ibeg)
    ib = iband 
    ioff(:,ib) = (/9,5/)
    igas(:,:,ib)  = RESHAPE ( SOURCE = (/ih2o,ico2,ih2o,nogas,  &
      ih2o,ico2,ih2o,nogas/), SHAPE = (/4,2/) )
    s_wght1(:,ib) = (/0.364641_wp, 0.364641_wp/)
    s_wght2(:,ib) = (/8.0_wp, 4.0_wp/)
    layreffr(ib)  = 30
    laysolfr_ref(ib)  = klev

    DO jg=1,ng(2)
      ig = jg + ibeg - 1
      rayl(ig) = rayl17
      abs_ab(1:585,1,ig) = absa17(:,jg)
      abs_ab(1:1175,2,ig) = absb17(:,jg)
      ref_self(1:10,ig) = selfref17(:,jg)
      ref_frgn(1:4 ,ig) = forref17(:,jg)
      ref_sflx(1:5,ig) = sfluxref17(jg,:)
    END DO

    !  Prepare band 18:  4000-4650 cm-1 (low - h2o,ch4; high - ch4)
    ibeg  = ngs(ib) + 1
    iband = ngb(ibeg)
    ib = iband 
    ioff(:,ib)   = (/9,1/)
    igas(:,:,ib) = RESHAPE ( SOURCE = (/ih2o,ich4,ih2o,nogas,&
      ich4,nogas,nogas,nogas/), SHAPE = (/4,2/) )
    s_wght1(:,ib) = (/38.9589_wp, 0.0_wp/)
    s_wght2(:,ib) = (/8.0_wp, 0.0_wp/)
    layreffr(ib) = 6
    laysolfr_ref(ib) = 0

    DO jg=1,ng(3)
      ig = jg + ibeg - 1
      rayl(ig) = rayl18
      abs_ab(1:585,1,ig) = absa18(:,jg)
      abs_ab(1:235,2,ig) = absb18(:,jg)
      ref_self(1:10,ig) = selfref18(:,jg)
      ref_frgn(1:3 ,ig) = forref18(:,jg)
      ref_sflx(1:9,ig) = sfluxref18(jg,:)
    END DO

    !  Prepare band 19: 4650-5150 cm-1 (low - h2o,co2; high - co2)
    ibeg  = ngs(ib) + 1
    iband = ngb(ibeg)
    ib = iband 
    ioff(:,ib) = (/9,1/)
    igas(:,:,ib) = RESHAPE(SOURCE = (/ih2o,ico2,ih2o,nogas,    &
      ico2,nogas,nogas,nogas/), SHAPE = (/4,2/) )
    s_wght1(:,ib) = (/5.49281_wp, 0.0_wp/)
    s_wght2(:,ib) = (/8.0_wp, 0.0_wp/)
    layreffr(ib)  = 3
    laysolfr_ref(ib)  = 0

    DO jg=1,ng(4)
      ig = jg + ibeg - 1
      rayl(ig) = rayl19
      abs_ab(1:585,1,ig) = absa19(:,jg)
      abs_ab(1:235,2,ig) = absb19(:,jg)
      ref_self(1:10,ig) = selfref19(:,jg)
      ref_frgn(1:3 ,ig) = forref19(:,jg)
      ref_sflx(1:9,ig) = sfluxref19(jg,:)
    END DO

    !  Prepare band 20: 5150-6150 cm-1 (low - h2o; high - h2o)
    ibeg  = ngs(ib) + 1
    iband = ngb(ibeg)
    ib = iband 
    ioff(:,ib) = (/1,1/)
    igas(:,:,ib) = RESHAPE(SOURCE = (/ih2o,nogas,ih2o,ich4,  &
      ih2o,nogas,ih2o,ich4/), SHAPE = (/4,2/) )
    s_wght1(:,ib) = (/0.0_wp, 0.0_wp/)
    s_wght2(:,ib) = (/0.0_wp, 0.0_wp/)
    layreffr(ib) = -3
    laysolfr_ref(ib) = 0

    DO jg=1,ng(5)
      ig = jg + ibeg - 1
      rayl(ig) = rayl20
      xk_mnr(1,ig) = absch420(jg)
      xk_mnr(2,ig) = absch420(jg)
      abs_ab(1:65,1,ig) = absa20(:,jg)
      abs_ab(1:235,2,ig) = absb20(:,jg)
      ref_self(1:10,ig) = selfref20(:,jg)
      ref_frgn(1:4 ,ig) = forref20(:,jg)
      ref_sflx(1,ig) = sfluxref20(jg)
    END DO

    !  Prepare band 21: 6150-7700 cm-1 (low - h2o,co2; high - h2o,co2)
    ibeg  = ngs(ib) + 1
    iband = ngb(ibeg)
    ib = iband 
    ioff(:,ib) = (/9,5/)
    igas(:,:,ib)  = RESHAPE ( SOURCE = (/ih2o,ico2,ih2o,nogas, &
      ih2o,ico2,ih2o,nogas/), SHAPE = (/4,2/) )
    s_wght1(:,ib) = (/0.0045321_wp, 0.0045321_wp/)
    s_wght2(:,ib) = (/8.0_wp, 4.0_wp/)
    layreffr(ib)  = 8
    laysolfr_ref(ib) = 0

    DO jg=1,ng(6)
      ig = jg + ibeg - 1
      rayl(ig) = rayl21
      abs_ab(1:585,1,ig) = absa21(:,jg)
      abs_ab(1:1175,2,ig) = absb21(:,jg)
      ref_self(1:10,ig) = selfref21(:,jg)
      ref_frgn(1:4 ,ig) = forref21(:,jg)
      ref_sflx(1:9,ig) = sfluxref21(jg,:)
    END DO

    ! Prepare band 22: 6150-7700 cm-1 (low - h2o,co2; high - h2o,co2)
    ! In this band the ratio of total O2 band intensity (lines 
    ! and Mate continuum) to O2 band intensity (line only) is 1.6 and is used
    ! to adjust the optical depths since the k's include only lines.  This is
    ! done by multiplying swght1 by 1.6 and abs_ab(:,2,:) by 1.6 to account for 
    ! this difference in below and above 100 hPa respectively.  Also note that
    ! the minor gas does not have a g-point dependent absorption because the
    ! o2 abosrption is a continuum effect
    ibeg  = ngs(ib) + 1
    iband = ngb(ibeg)
    ib = iband 
    ioff(:,ib) = (/9,1/)
    igas(:,:,ib) = RESHAPE(SOURCE = (/ih2o,io2,ih2o,io2,    &
      io2,nogas,nogas,io2/), SHAPE = (/4,2/) )
    s_wght1(:,ib) = (/0.022708_wp, 0.0_wp/)*1.6_wp
    s_wght2(:,ib) = (/8.0_wp, 0.0_wp/)
    layreffr(ib) = 2
    laysolfr_ref(ib) = 0

    DO jg=1,ng(7)
      ig = jg + ibeg - 1
      rayl(ig) = rayl22
      xk_mnr(1,ig) = 4.35e-4_wp/(350.0_wp*2.0_wp) 
      xk_mnr(2,ig) = 4.35e-4_wp/(350.0_wp*2.0_wp) 
      abs_ab(1:585,1,ig) = absa22(:,jg)
      abs_ab(1:235,2,ig) = absb22(:,jg)*1.6_wp
      ref_self(1:10,ig) = selfref22(:,jg)
      ref_frgn(1:3 ,ig) = forref22(:,jg)
      ref_sflx(1:9,ig) = sfluxref22(jg,:)
    END DO

    !  Prepare band 23:  8050-12850 cm-1 (low - h2o; high - nothing)
    !  Average Giver et al. correction factor for this band is 1.029
    !  and is implemented by multiplying the abosrption coefficients
    ibeg  = ngs(ib) + 1
    iband = ngb(ibeg)
    ib = iband 
    ioff(:,ib) = (/1,1/)
    igas(:,:,ib) = RESHAPE ( SOURCE = (/ih2o,nogas,ih2o,nogas,     &
      nogas,nogas,nogas,nogas/), SHAPE = (/4,2/) )
    s_wght1(:,ib) = (/0.0_wp, 0.0_wp/)
    s_wght2(:,ib) = (/0.0_wp, 0.0_wp/)
    layreffr(ib)  = -6
    laysolfr_ref(ib)  = 0

    DO jg=1,ng(8)
      ig = jg + ibeg - 1
      rayl(ig) = rayl23(jg)
      abs_ab(1:65,1,ig) = absa23(:,jg) * 1.029_wp
      ref_self(1:10,ig) = selfref23(:,jg)
      ref_frgn(1:3 ,ig) = forref23(:,jg)
      ref_sflx(1,ig) = sfluxref23(jg)
    END DO
    
    !  Prepare band 24:  12850-16000 cm-1 (low - h2o,o2; high - o2)
    ibeg  = ngs(ib) + 1
    iband = ngb(ibeg)
    ib = iband 
    ioff(:,ib) = (/9,1/)
    igas(:,:,ib) = RESHAPE(SOURCE = (/ih2o,io2,ih2o,io3,     &
      io2,nogas,nogas,io3/), SHAPE = (/4,2/) )
    s_wght1(:,ib) = (/0.124692_wp, 0.0_wp/)
    s_wght2(:,ib) = (/8.0_wp, 0.0_wp/)
    layreffr(ib) = 1
    laysolfr_ref(ib) = 0

    DO jg=1,ng(9)
      ig = jg + ibeg - 1
      rayl_24_js(:,ig)  = rayla24(jg,:)
      rayl(ig) = raylb24(jg)
      xk_mnr(1,ig) = abso3a24(jg)
      xk_mnr(2,ig) = abso3b24(jg)
      abs_ab(1:585,1,ig) = absa24(:,jg)
      abs_ab(1:235,2,ig) = absb24(:,jg)
      ref_self(1:10,ig) = selfref24(:,jg)
      ref_frgn(1:3 ,ig) = forref24(:,jg)
      ref_sflx(1:9,ig) = sfluxref24(jg,:)
    END DO

    !  Prepare band 25: 16000-22650 cm-1 (low - h2o; high - nothing)
    ibeg  = ngs(ib) + 1
    iband = ngb(ibeg)
    ib = iband 
    ioff(:,ib) = (/1,1/)
    igas(:,:,ib) = RESHAPE(SOURCE = (/ih2o,nogas,nogas,io3,  &
      nogas,nogas,nogas,io3/), SHAPE = (/4,2/) )
    s_wght1(:,ib) = (/0.0_wp, 0.0_wp/)
    s_wght2(:,ib) = (/0.0_wp, 0.0_wp/)
    layreffr(ib)  = -2
    laysolfr_ref(ib)  = 0

    DO jg=1,ng(10)
      ig = jg + ibeg - 1
      rayl(ig) = rayl25(jg)
      xk_mnr(1,ig) = abso3a25(jg)
      xk_mnr(2,ig) = abso3b25(jg)
      abs_ab(1:65,1,ig) = absa25(:,jg)
      ref_sflx(1,ig) = sfluxref25(jg)
    END DO

    !  Prepare band 26: 22650-29000 cm-1 (low - nothing; high - nothing)
    ibeg  = ngs(ib) + 1
    iband = ngb(ibeg)
    ib = iband 
    ioff(:,ib) = (/1,1/)
    igas(:,:,ib) = RESHAPE(SOURCE = (/nogas,nogas,nogas,nogas,  &
      nogas,nogas,nogas,nogas/), SHAPE = (/4,2/) )
    s_wght1(:,ib) = (/0.0_wp, 0.0_wp/)
    s_wght2(:,ib) = (/0.0_wp, 0.0_wp/)
    layreffr(ib)  = -1
    laysolfr_ref(ib)  = 0

    DO jg=1,ng(11)
      ig = jg + ibeg - 1
      rayl(ig) = rayl26(jg)
      ref_sflx(1,ig)  = sfluxref26(jg)
    END DO

    !  Prepare band 27:  29000-38000 cm-1 (low - o3; high - o3)
    ! Kurucz solar source function
    ! The values in sfluxref were obtained using the "low resolution"
    ! version of the Kurucz solar source function.  For unknown reasons,
    ! the total irradiance in this band differs from the corresponding
    ! total in the "high-resolution" version of the Kurucz function.
    ! Therefore, these values are scaled below by the factor SCALEKUR.
    ibeg  = ngs(ib) + 1
    iband = ngb(ibeg)
    ib = iband 
    ioff(:,ib) = (/1,1/)
    igas(:,:,ib) = RESHAPE(SOURCE = (/io3,nogas,nogas,nogas,  &
      io3,nogas,nogas,nogas/), SHAPE = (/4,2/) )
    s_wght1(:,ib) = (/0.0_wp, 0.0_wp/)
    s_wght2(:,ib) = (/0.0_wp, 0.0_wp/)
    layreffr(ib)  = -32
    laysolfr_ref(ib)  = klev

    DO jg=1,ng(12)
      ig = jg + ibeg - 1
      rayl(ig) = rayl27(jg)
      abs_ab(1:65,1,ig) = absa27(:,jg)
      abs_ab(1:235,2,ig) = absb27(:,jg)
      ref_sflx(1,ig) = sfluxref27(jg) * 50.15_wp/48.37_wp
    END DO

    !  Prepare band 28:   38000-50000 cm-1 (low - o3,o2; high - o3,o2)
    ibeg  = ngs(ib) + 1
    iband = ngb(ibeg)
    ib = iband 
    ioff(:,ib) = (/9,5/)
    igas(:,:,ib) = RESHAPE(SOURCE = (/io3,io2,nogas,nogas,  &
      io3,io2,nogas,nogas/), SHAPE = (/4,2/) )
    s_wght1(:,ib) = (/6.67029e-07_wp, 6.67029e-07_wp/)
    s_wght2(:,ib) = (/8.0_wp, 4.0_wp/)
    layreffr(ib)  = 58
    laysolfr_ref(ib)  = klev

    DO jg=1,ng(13)
      ig = jg + ibeg - 1
      rayl(ig) = rayl28
      abs_ab(1:585,1,ig) = absa28(:,jg)
      abs_ab(1:1175,2,ig) = absb28(:,jg)
      ref_sflx(1:5,ig) = sfluxref28(jg,:)
    END DO

    !  Prepare band 29:   820-2600 cm-1 (low - h2o; high - co2)
    ibeg  = ngs(ib) + 1
    iband = ngb(ibeg)
    ib = iband
    ioff(:,ib) = (/1,1/)
    igas(:,:,ib) = RESHAPE(SOURCE = (/ih2o,nogas,ih2o,ico2,  &
      ico2,nogas,nogas,ih2o/), SHAPE = (/4,2/) )
    s_wght1(:,ib) = (/0.0_wp, 0.0_wp/)
    s_wght2(:,ib) = (/0.0_wp, 0.0_wp/)
    layreffr(ib)  = -49
    laysolfr_ref(ib)  = klev

    DO jg=1,ng(14)
      ig = jg + ibeg - 1
      rayl(ig) = rayl29
      abs_ab(1:65,1,ig) = absa29(:,jg)
      abs_ab(1:235,2,ig) = absb29(:,jg)
      xk_mnr(1,ig) = absco229(jg)
      xk_mnr(2,ig) = absh2o29(jg)
      ref_self(1:10,ig) = selfref29(:,jg)
      ref_frgn(1:4 ,ig) = forref29(:,jg)
      ref_sflx(1,ig) = sfluxref29(jg)
    END DO

    initialized = .TRUE.
  END SUBROUTINE setup_taumol

END MODULE mo_psrad_srtm_gas_optics
