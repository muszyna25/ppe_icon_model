!  USE mo_psrad_srtm_driver,   ONLY : srtm

MODULE mo_psrad_c_binding

  USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_int, c_double, c_ptr, c_f_pointer
  IMPLICIT NONE

CONTAINS

  SUBROUTINE psrad_set_amw(val) BIND(C)
    USE mo_psrad_general, ONLY : amw
    REAL(c_double), VALUE, INTENT(IN) :: val
    amw = val
  END SUBROUTINE psrad_set_amw
    
  REAL(c_double) FUNCTION psrad_get_amw() BIND(c)
    USE mo_psrad_general, ONLY : amw
    psrad_get_amw = amw
  END FUNCTION

  SUBROUTINE psrad_set_amd(val) BIND(C)
    USE mo_psrad_general, ONLY : amd
    REAL(c_double), VALUE, INTENT(IN) :: val
    amd = val
  END SUBROUTINE psrad_set_amd

  REAL(c_double) FUNCTION psrad_get_amd() BIND(c)
    USE mo_psrad_general, ONLY : amd
    psrad_get_amd = amd
  END FUNCTION

  SUBROUTINE psrad_set_grav(val) BIND(C)
    USE mo_psrad_general, ONLY : grav
    REAL(c_double), VALUE, INTENT(IN) :: val
    grav = val
  END SUBROUTINE psrad_set_grav

  REAL(c_double) FUNCTION psrad_get_grav() BIND(c)
    USE mo_psrad_general, ONLY : grav
    psrad_get_grav = grav
  END FUNCTION

  SUBROUTINE psrad_set_pi(val) BIND(C)
    USE mo_psrad_general, ONLY : pi
    REAL(c_double), VALUE, INTENT(IN) :: val
    pi = val
  END SUBROUTINE psrad_set_pi

  REAL(c_double) FUNCTION psrad_get_pi() BIND(c)
    USE mo_psrad_general, ONLY : pi
    psrad_get_pi = pi
  END FUNCTION

  INTEGER(c_int) FUNCTION psrad_get_mg() BIND(c)
    USE mo_psrad_general, ONLY : mg
    psrad_get_mg = mg
  END FUNCTION

  INTEGER(c_int) FUNCTION psrad_get_maxxsec() BIND(c)
    USE mo_psrad_general, ONLY : ncfc
    psrad_get_maxxsec = ncfc
  END FUNCTION

  INTEGER(c_int) FUNCTION psrad_get_ngas() BIND(c)
    USE mo_psrad_general, ONLY : ngas
    psrad_get_ngas = ngas
  END FUNCTION

  INTEGER(c_int) FUNCTION psrad_get_jpband() BIND(c)
    USE mo_psrad_general, ONLY : jpband
    psrad_get_jpband = jpband
  END FUNCTION

  INTEGER(c_int) FUNCTION psrad_get_nbndsw() BIND(c)
    USE mo_psrad_general, ONLY : nbndsw
    psrad_get_nbndsw = nbndsw
  END FUNCTION

  INTEGER(c_int) FUNCTION psrad_get_ngptsw() BIND(c)
    USE mo_psrad_general, ONLY : ngptsw
    psrad_get_ngptsw = ngptsw
  END FUNCTION

  INTEGER(c_int) FUNCTION psrad_get_nbndlw() BIND(c)
    USE mo_psrad_general, ONLY : nbndlw
    psrad_get_nbndlw = nbndlw
  END FUNCTION

  INTEGER(c_int) FUNCTION psrad_get_ngptlw() BIND(c)
    USE mo_psrad_general, ONLY : ngptlw
    psrad_get_ngptlw = ngptlw
  END FUNCTION

  SUBROUTINE psrad_cloud_optics(laglac, laland, kproma, kbdim, klev, ktype, &
    icldlyr, zlwp, ziwp, zlwc, ziwc, zcdnc, tau_lw, tau_sw, &
    omg, asy, re_droplets2d, re_crystals2d) BIND(c)
    USE mo_psrad_general, ONLY : nbndsw, nbndlw
    USE mo_psrad_cloud_optics, ONLY : cloud_optics
    INTEGER(c_int), VALUE, INTENT (IN) :: kproma, kbdim, klev
    INTEGER(c_int), INTENT (IN) :: ktype(KBDIM), icldlyr(KBDIM,klev), &
       laglac(KBDIM), laland(KBDIM) 
    REAL (c_double), INTENT (IN) :: zlwp(KBDIM,klev), ziwp(KBDIM,klev), &
       zcdnc(KBDIM,klev), zlwc(KBDIM,klev), ziwc(KBDIM,klev)

    REAL (c_double), INTENT (OUT) :: tau_lw(KBDIM,klev,nbndlw), &
         tau_sw(KBDIM,klev,nbndsw), omg(KBDIM,klev,nbndsw), &
         asy(KBDIM,klev,nbndsw), re_droplets2d(KBDIM,klev), &
         re_crystals2d(KBDIM,klev)

    LOGICAL, ALLOCATABLE, SAVE :: laglac_f(:), laland_f(:)
    IF (.not. allocated(laglac_f)) THEN
      ALLOCATE(laglac_f(KBDIM), laland_f(KBDIM))
    END IF
    laglac_f(1:KBDIM) = laglac(1:KBDIM) .ne. 0
    laland_f(1:KBDIM) = laland(1:KBDIM) .ne. 0

    CALL cloud_optics(laglac_f, laland_f, kproma, KBDIM, klev, ktype, &
      icldlyr, zlwp, ziwp, zlwc, ziwc, zcdnc, &
      tau_lw, tau_sw, omg, asy, re_droplets2d, re_crystals2d)
    
  END SUBROUTINE

  SUBROUTINE psrad_basic_setup() BIND(c)
    USE mo_psrad_cloud_optics,  ONLY : setup_cloud_optics
    USE mo_psrad_lrtm_setup,  ONLY : setup_lrtm
    USE mo_psrad_srtm_setup,  ONLY : setup_srtm
    CALL setup_cloud_optics
    CALL setup_lrtm
    CALL setup_srtm
  END SUBROUTINE


  REAL(c_double) FUNCTION psrad_get_ssi_default(i) BIND(c)
    USE mo_psrad_srtm_setup,  ONLY : ssi_default
    INTEGER(c_int), VALUE, INTENT(IN) :: i

    psrad_get_ssi_default = ssi_default(i)
  END FUNCTION

  INTEGER(c_int) FUNCTION psrad_get_num_gpts_lw() BIND(c)
    USE mo_psrad_general, ONLY : ngptlw
    psrad_get_num_gpts_lw = ngptlw
  END FUNCTION

  INTEGER(c_int) FUNCTION psrad_get_num_gpts_sw() BIND(c)
    USE mo_psrad_general, ONLY : ngptsw
    psrad_get_num_gpts_sw = ngptsw
  END FUNCTION

  SUBROUTINE psrad_lrtm(kproma, kbdim, klev, play, psfc, tlay, tlev, tsfc, &
    wkl, wx, coldry, emis, cldfr, taucld, tauaer, rnseeds, seed_size, &
    uflx, dflx, uflxc, dflxc) BIND(c)

    USE mo_psrad_lrtm_driver,   ONLY : lrtm
    USE mo_psrad_general, ONLY : ncfc, nbndlw, ngas

    INTEGER(c_int), VALUE, INTENT(IN) :: kbdim, kproma, klev, seed_size

    INTEGER(c_int), INTENT(INOUT) :: rnseeds(KBDIM, seed_size)
    REAL(c_double), INTENT(IN) :: play(KBDIM,klev), psfc(KBDIM), &
      tlay(KBDIM,klev), tlev(KBDIM,klev+1), tsfc(KBDIM), &
      wkl(KBDIM,klev,ngas), wx(KBDIM,ncfc,klev), &
      coldry(KBDIM,klev), emis(KBDIM,nbndlw), cldfr(KBDIM,klev), &
      taucld(KBDIM,klev,nbndlw), tauaer(KBDIM,klev,nbndlw)
    REAL(c_double), INTENT(OUT) :: uflx(KBDIM,klev+1), dflx(KBDIM,klev+1), &
      uflxc(KBDIM,klev+1), dflxc(KBDIM,klev+1)

    CALL lrtm(kproma, KBDIM, klev, play, psfc, tlay, tlev, tsfc, &
      wkl, wx, coldry, emis, cldfr, taucld, tauaer, rnseeds, &
      uflx, dflx, uflxc, dflxc)

  END SUBROUTINE

  SUBROUTINE psrad_srtm(kproma, kbdim, klev, play, tlay, wkl, coldry, &
    asdir, asdif, aldir, aldif, prmu0, daylight_frc, &
    ssi_factor, psctm, &
    cld_frc, cld_tau_sw, cld_cg_sw, cld_piz_sw, aer_tau_sw, aer_cg_sw, &
    aer_piz_sw, rnseeds, seed_size, flxd_sw, flxu_sw, flxd_sw_clr, &
    flxu_sw_clr, vis_dn_dir_sfc, par_dn_dir_sfc, nir_dn_dir_sfc, &
    vis_dn_dff_sfc, par_dn_dff_sfc, nir_dn_dff_sfc, vis_up_sfc, par_up_sfc, &
    nir_up_sfc) BIND(C)

    USE mo_psrad_srtm_driver, ONLY : srtm
    USE mo_psrad_general, ONLY : nbndsw, ngas

    INTEGER(c_int), VALUE, INTENT(in) :: kproma, kbdim, klev, seed_size
    REAL(c_double), VALUE, INTENT(in) :: psctm

    INTEGER(c_int), INTENT(INOUT) :: rnseeds(KBDIM, seed_size)
    REAL(c_double), DIMENSION(KBDIM,klev), INTENT(IN) :: play, &
      tlay, coldry, cld_frc
    REAL(c_double), DIMENSION(KBDIM), INTENT(IN) :: asdir, aldir, asdif, &
      aldif, prmu0
    REAL(c_double), INTENT(IN) ::  wkl(KBDIM,klev,ngas), &
      daylight_frc(KBDIM), ssi_factor(nbndsw)

    REAL(c_double), DIMENSION(KBDIM,klev,nbndsw), INTENT(IN) :: &
      cld_tau_sw, cld_cg_sw, cld_piz_sw, aer_tau_sw, aer_cg_sw, aer_piz_sw

    REAL(c_double), DIMENSION(KBDIM,klev+1), INTENT(OUT) :: flxd_sw, &
      flxu_sw, flxd_sw_clr, flxu_sw_clr
    
    REAL(c_double), DIMENSION(KBDIM), INTENT(OUT) :: vis_dn_dir_sfc, &
      par_dn_dir_sfc, nir_dn_dir_sfc, vis_dn_dff_sfc, par_dn_dff_sfc, &
      nir_dn_dff_sfc, vis_up_sfc, par_up_sfc, nir_up_sfc

    CALL srtm(kproma, KBDIM, klev, play, tlay, wkl, coldry, &
      asdir, asdif, aldir, aldif, prmu0, daylight_frc, &
      ssi_factor, psctm, &
      cld_frc, cld_tau_sw, cld_cg_sw, cld_piz_sw, aer_tau_sw, aer_cg_sw, &
      aer_piz_sw, rnseeds, flxd_sw, flxu_sw, flxd_sw_clr, flxu_sw_clr, &
      vis_dn_dir_sfc, par_dn_dir_sfc, nir_dn_dir_sfc, vis_dn_dff_sfc, &
      par_dn_dff_sfc, nir_dn_dff_sfc, vis_up_sfc, par_up_sfc, &
      nir_up_sfc)

  END SUBROUTINE
END MODULE mo_psrad_c_binding
