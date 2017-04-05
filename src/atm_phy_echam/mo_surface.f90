!>
!!
!! @author Hui Wan, MPI-M
!!
!! @par Revision History
!! First version by Hui Wan (2011-08)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_surface

  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: finish
  USE mo_surface_diag,      ONLY: wind_stress, surface_fluxes
  USE mo_echam_vdiff_params,ONLY: tpfac2
  USE mo_echam_phy_config,  ONLY: phy_config => echam_phy_config
  USE mo_vdiff_solver,      ONLY: ih, iqv, iu, iv, imh, imqv, imuv, &
                                & nmatrix, nvar_vdiff,              &
                                & matrix_to_richtmyer_coeff
  USE mo_coupling_config,   ONLY: is_coupled_run

#ifndef __NO_JSBACH__
  USE mo_jsb_interface,     ONLY: jsbach_interface
  USE mo_radiation_config,  ONLY: mmr_co2      ! This should be here only temporarily
#endif
  USE mo_echam_sfc_indices, ONLY: nsfc_type
  USE mo_echam_phy_memory,  ONLY: cdimissval
#ifndef __NO_ICON_OCEAN__
  USE mo_sea_ice,           ONLY: ice_fast
  USE mo_ml_ocean,            ONLY: ml_ocean
#endif
  USE mo_physical_constants,ONLY: rhos, rhoi, Tf, alf, albedoW, zemiss_def, stbo, tmelt
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: update_surface

CONTAINS
  !>
  !!
  !!
  SUBROUTINE update_surface( lsfc_heat_flux, lsfc_mom_flux,     &! in
                           & pdtime, psteplen,                  &! in
                           & jg,                                &! in
                           & kproma, kbdim,                     &! in
                           & kice,                              &! in
                           & klev, ksfc_type,                   &! in
                           & idx_wtr, idx_ice, idx_lnd,         &! in
                           & pfrc,                              &! in
                           & pcfh_tile, pcfm_tile,              &! in
                           & pfac_sfc, pocu, pocv,              &! in
                           & aa, aa_btm, bb, bb_btm,            &! inout
                           & pcpt_tile, pqsat_tile,             &! inout
                           & ptsfc_tile,                        &! inout
                           & pu_stress_gbm, pv_stress_gbm,      &! out
                           & plhflx_gbm, pshflx_gbm,            &! out
                           & pevap_gbm,                         &! out
                           & pu_stress_tile,   pv_stress_tile,  &! out
                           & plhflx_tile, pshflx_tile,          &! out
                           & dshflx_dT_tile,                    &! out
                           & pevap_tile,                        &! out
                           !! optional
                           & nblock,                            &! in
                           & lsm,                               &! in
                           & alake,                             &! in
                           & pu,                                &! in
                           & pv,                                &! in
                           & ptemp,                             &! in
                           & pq,                                &! in
                           & prsfl,                             &! in
                           & prsfc,                             &! in
                           & pssfl,                             &! in
                           & pssfc,                             &! in
                           & plw,                               &! inout
                           & plw_down,                          &! in
                           & psw,                               &! inout
                           & pswvis, pswnir, pswpar_down,       &! in
                           & pvisdff, pnirdff, ppardff,         &! in
                           & presi_old,                         &! in
                           & pcosmu0,                           &! in
                           & pch_tile,                          &! in
                           !! for JSBACH
                           & pcsat,                             &! inout
                           & pcair,                             &! inout
                           & q_snocpymlt,                       &! out
                           !
                           & z0h_lnd, z0m_tile,                 &! out
                           & albvisdir, albnirdir, albvisdif, albnirdif, &! inout
                           & albvisdir_tile,                    &! inout
                           & albnirdir_tile,                    &! inout
                           & albvisdif_tile,                    &! inout
                           & albnirdif_tile,                    &! inout
                           & albedo, albedo_tile,               &! inout
                           & ptsfc,                             &! out
                           & ptsfc_rad,                         &! out
                           & pswflx_tile, plwflx_tile,          &! out
                           & lake_ice_frc,                      &! out
                           !! Sea ice
                           & Tsurf,                             &! inout
                           & T1,                                &! inout
                           & T2,                                &! inout
                           & hi,                                &! in
                           & hs,                                &! inout
                           & Qtop,                              &! out
                           & Qbot,                              &! out
                           & conc,                              &! in
                           & albvisdir_ice, albvisdif_ice,      &! inout
                           & albnirdir_ice, albnirdif_ice)       ! inout

    LOGICAL, INTENT(IN) :: lsfc_heat_flux, lsfc_mom_flux
    REAL(wp),INTENT(IN) :: pdtime, psteplen
    INTEGER, INTENT(IN) :: jg
    INTEGER, INTENT(IN) :: kproma, kbdim
    INTEGER, INTENT(IN) :: klev, ksfc_type
    INTEGER, INTENT(IN) :: idx_wtr, idx_ice, idx_lnd
    REAL(wp),INTENT(IN) :: pfrc      (kbdim,ksfc_type)
    REAL(wp),INTENT(IN) :: pcfh_tile (kbdim,ksfc_type)
    REAL(wp),INTENT(IN) :: pcfm_tile (kbdim,ksfc_type)
    REAL(wp),INTENT(IN) :: pfac_sfc  (kbdim)
    REAL(wp),INTENT(IN) :: pocu      (kbdim)
    REAL(wp),INTENT(IN) :: pocv      (kbdim)
    REAL(wp),INTENT(INOUT) :: aa     (kbdim,klev,3,nmatrix)
    REAL(wp),INTENT(INOUT) :: aa_btm (kbdim,3,ksfc_type,imh:imqv)
    REAL(wp),INTENT(INOUT) :: bb     (kbdim,klev,nvar_vdiff)
    REAL(wp),INTENT(INOUT) :: bb_btm (kbdim,ksfc_type,ih:iqv)
    REAL(wp),INTENT(INOUT) :: pcpt_tile (kbdim,ksfc_type)
    REAL(wp),INTENT(INOUT) :: pqsat_tile(kbdim,ksfc_type)
    REAL(wp),INTENT(INOUT) :: ptsfc_tile (kbdim,ksfc_type)

    REAL(wp),INTENT(INOUT)   :: pu_stress_gbm (kbdim)              ! OUT
    REAL(wp),INTENT(INOUT)   :: pv_stress_gbm (kbdim)              ! OUT
    REAL(wp),INTENT(INOUT)   ::    plhflx_gbm (kbdim)              ! OUT
    REAL(wp),INTENT(INOUT)   ::    pshflx_gbm (kbdim)              ! OUT
    REAL(wp),INTENT(INOUT)   ::     pevap_gbm (kbdim)              ! OUT

    REAL(wp),INTENT(INOUT)   :: pu_stress_tile (kbdim,ksfc_type)   ! OUT
    REAL(wp),INTENT(INOUT)   :: pv_stress_tile (kbdim,ksfc_type)   ! OUT
    REAL(wp),INTENT(INOUT)   ::    plhflx_tile (kbdim,ksfc_type)   ! OUT
    REAL(wp),INTENT(INOUT)   ::    pshflx_tile (kbdim,ksfc_type)   ! OUT
    REAL(wp),INTENT(INOUT)   :: dshflx_dT_tile (kbdim,ksfc_type)   ! OUT
    REAL(wp),INTENT(INOUT)   ::     pevap_tile (kbdim,ksfc_type)   ! OUT

    !! JSBACH input
    INTEGER, OPTIONAL,INTENT(IN) :: nblock
    REAL(wp),OPTIONAL,INTENT(IN) :: lsm(kbdim)
    REAL(wp),OPTIONAL,INTENT(IN) :: alake(kbdim)
    REAL(wp),OPTIONAL,INTENT(IN) :: pu        (kbdim)              ! zonal wind lowest level
    REAL(wp),OPTIONAL,INTENT(IN) :: pv        (kbdim)              ! meridional wind lowest level
    REAL(wp),OPTIONAL,INTENT(IN) :: ptemp     (kbdim)              ! temperature of lowest atmospheric level
    REAL(wp),OPTIONAL,INTENT(IN) :: pq        (kbdim)              ! humidity of lowest atmospheric level
    REAL(wp),OPTIONAL,INTENT(IN) :: prsfl     (kbdim)              ! rain large scale
    REAL(wp),OPTIONAL,INTENT(IN) :: prsfc     (kbdim)              ! rain convective
    REAL(wp),OPTIONAL,INTENT(IN) :: pssfl     (kbdim)              ! snow large scale
    REAL(wp),OPTIONAL,INTENT(IN) :: pssfc     (kbdim)              ! snow convective
    REAL(wp),OPTIONAL,INTENT(IN) :: plw_down  (kbdim)              ! downward surface longwave flux [W/m2]
    REAL(wp),OPTIONAL,INTENT(IN) :: pswvis    (kbdim)              ! net surface shortwave flux in VIS [W/m2]
    REAL(wp),OPTIONAL,INTENT(IN) :: pswnir    (kbdim)              ! net surface shortwave flux in NIR [W/m2]
    REAL(wp),OPTIONAL,INTENT(IN) :: pswpar_down(kbdim)             ! downward surface shortwave flux in PAR [W/m2]
    REAL(wp),OPTIONAL,INTENT(IN) :: pvisdff   (kbdim)              ! diffuse fraction of VIS shortwave flux
    REAL(wp),OPTIONAL,INTENT(IN) :: pnirdff   (kbdim)              ! diffuse fraction of NIR shortwave flux
    REAL(wp),OPTIONAL,INTENT(IN) :: ppardff   (kbdim)              ! diffuse fraction of PAR shortwave flux
    REAL(wp),OPTIONAL,INTENT(IN) :: presi_old (kbdim,klev+1)       ! half level pressure
    REAL(wp),OPTIONAL,INTENT(IN) :: pcosmu0   (kbdim)              ! cos of zenith angle
    REAL(wp),OPTIONAL,INTENT(IN) :: pch_tile  (kbdim,ksfc_type)
    !! JSBACH output
    REAL(wp),OPTIONAL,INTENT(INOUT) :: pcsat(kbdim)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: pcair(kbdim)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: q_snocpymlt(kbdim)  ! OUT
    REAL(wp),OPTIONAL,INTENT(INOUT) :: z0h_lnd(kbdim), z0m_tile(kbdim,ksfc_type)  ! OUT
    !
    REAL(wp),OPTIONAL,INTENT(INOUT) :: albvisdir_tile(kbdim,ksfc_type)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: albnirdir_tile(kbdim,ksfc_type)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: albvisdif_tile(kbdim,ksfc_type)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: albnirdif_tile(kbdim,ksfc_type)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: albedo(kbdim)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: albvisdir(kbdim), albvisdif(kbdim)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: albnirdir(kbdim), albnirdif(kbdim)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: albedo_tile(kbdim,ksfc_type)
    REAL(wp),OPTIONAL,INTENT(OUT)   :: ptsfc    (kbdim) ! OUT
    REAL(wp),OPTIONAL,INTENT(OUT)   :: ptsfc_rad(kbdim) ! OUT
    REAL(wp),OPTIONAL,INTENT(INOUT) :: plw       (kbdim)            ! INOUT net surface longwave flux [W/m2]
    REAL(wp),OPTIONAL,INTENT(INOUT) :: psw       (kbdim)            ! IOUT net surface shortwave flux [W/m2]
    REAL(wp),OPTIONAL,INTENT(INOUT) :: pswflx_tile(kbdim,ksfc_type) ! OUT
    REAL(wp),OPTIONAL,INTENT(INOUT) :: plwflx_tile(kbdim,ksfc_type) ! OUT
    REAL(wp),OPTIONAL,INTENT(OUT)   :: lake_ice_frc(kbdim)          ! OUT
    !! Sea ice
    INTEGER,          INTENT(IN)    :: kice ! Number of ice thickness classes
    REAL(wp),OPTIONAL,INTENT(INOUT) :: Tsurf(kbdim,kice)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: T1   (kbdim,kice)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: T2   (kbdim,kice)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: hi   (kbdim,kice)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: hs   (kbdim,kice)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: Qtop (kbdim,kice) ! OUT
    REAL(wp),OPTIONAL,INTENT(INOUT) :: Qbot (kbdim,kice) ! OUT
    REAL(wp),OPTIONAL,INTENT(IN)    :: conc (kbdim,kice)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: albvisdir_ice(kbdim,kice)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: albvisdif_ice(kbdim,kice)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: albnirdir_ice(kbdim,kice)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: albnirdif_ice(kbdim,kice)

! locals

    INTEGER  :: jsfc, jk, jkm1, im, k
    REAL(wp) :: se_sum(kbdim), qv_sum(kbdim), wgt_sum(kbdim), wgt(kbdim)
    REAL(wp) :: zca(kbdim,ksfc_type), zcs(kbdim,ksfc_type)
    REAL(wp) :: zfrc_oce(kbdim)

    REAL(wp) :: zen_h (kbdim,ksfc_type)
    REAL(wp) :: zfn_h (kbdim,ksfc_type)
    REAL(wp) :: zen_qv(kbdim,ksfc_type)
    REAL(wp) :: zfn_qv(kbdim,ksfc_type)

    REAL(wp) :: &
      & sat_surface_specific_humidity(kbdim), dry_static_energy(kbdim),                             &
      & ztsfc_lnd(kbdim), ztsfc_lnd_eff(kbdim), ztsfc_wtr(kbdim), ztsfc_ice(kbdim),                 &
      & zswvis_down(kbdim), zswnir_down(kbdim), zsw_down(kbdim),                                    &
      & zswvisdif_down(kbdim), zswvisdir_down(kbdim), zswnirdif_down(kbdim), zswnirdir_down(kbdim), &
      & zalbvis(kbdim), zalbnir(kbdim)

    REAL(wp) :: zgrnd_hflx(kbdim,ksfc_type), zgrnd_hcap(kbdim,ksfc_type), ztsfc(kbdim), &
      & zevap_wtr(kbdim), zevap_ice(kbdim), zlhflx_wtr(kbdim), zlhflx_ice(kbdim),       &
      & zshflx_wtr(kbdim), zshflx_ice(kbdim), zalbvisdir_ice(kbdim)

    !REAL(wp) :: zt2s_conv(kbdim,ksfc_type)

    ! Sea ice
    REAL(wp) :: Tfw(kbdim)
    REAL(wp) :: swflx_ice(kbdim,kice), nonsolar_ice(kbdim,kice), dnonsolardT(kbdim,kice), conc_sum(kbdim)

    LOGICAL :: mask(kbdim)

#if defined(__NO_JSBACH__) || defined (__NO_ICON_OCEAN__)
   CHARACTER(len=*), PARAMETER :: method_name='mo_surface:update_surface'
#endif

    ! save old grid box mean surface temperature
    ztsfc(:) = 0._wp
    DO jsfc=1,ksfc_type
      ztsfc(1:kproma) = ztsfc(1:kproma) + pfrc(1:kproma,jsfc) * ptsfc_tile(1:kproma,jsfc)
    ENDDO

    ! Compute factor for conversion temperature to dry static energy
    !DO jsfc=1,ksfc_type
    !  zt2s_conv(1:kproma,jsfc) = pcpt_tile(1:kproma,jsfc) / ptsfc_tile(1:kproma,jsfc)
    !END DO

    !===================================================================
    ! BEFORE CALLING land/ocean/ice model
    !===================================================================
    ! Compute wind stress at the old time step.
    ! At this point bb(:,klev,iu) = u_klev(t)/tpfac1 (= udif in echam)
    !               bb(:,klev,iv) = v_klev(t)/tpfac1 (= vdif in echam)

    CALL wind_stress( lsfc_mom_flux, psteplen,             &! in
                    & kproma, kbdim, ksfc_type,            &! in
                    & pfrc, pcfm_tile, pfac_sfc,           &! in
                    & bb(:,klev,iu), bb(:,klev,iv),        &! in
                    & pu_stress_gbm,  pv_stress_gbm,       &! out
                    & pu_stress_tile, pv_stress_tile       )! out

    ! Compute downward shortwave surface fluxes
    zswvisdif_down(1:kproma)  =          pvisdff(1:kproma)  * pswvis(1:kproma) / (1._wp - albvisdif(1:kproma))
    zswvisdir_down(1:kproma)  = (1._wp - pvisdff(1:kproma)) * pswvis(1:kproma) / (1._wp - albvisdir(1:kproma))
    zswvis_down(1:kproma)     = zswvisdif_down(1:kproma) + zswvisdir_down(1:kproma)
    zswnirdif_down(1:kproma)  =          pnirdff(1:kproma)  * pswnir(1:kproma) / (1._wp - albnirdif(1:kproma))
    zswnirdir_down(1:kproma)  = (1._wp - pnirdff(1:kproma)) * pswnir(1:kproma) / (1._wp - albnirdir(1:kproma))
    zswnir_down(1:kproma)     = zswnirdif_down(1:kproma) + zswnirdir_down(1:kproma)
    zsw_down(1:kproma)        = zswvis_down(1:kproma) + zswnir_down(1:kproma)

    ! Turbulent transport of moisture:
    ! - finish matrix set up;
    ! - perform bottom level elimination;
    ! - convert matrix entries to Richtmyer-Morton coefficients
    IF (idx_lnd <= ksfc_type) THEN
      CALL matrix_to_richtmyer_coeff( kproma, kbdim, klev, ksfc_type, idx_lnd, &! in
        & aa(:,:,:,imh:imqv), bb(:,:,ih:iqv),      &! in
        & aa_btm, bb_btm,                          &! inout
        & zen_h, zfn_h, zen_qv, zfn_qv,            &! out
        & pcair = pcair(:),                        &! in
        & pcsat = pcsat(:))                         ! in
    ELSE
      CALL matrix_to_richtmyer_coeff( kproma, kbdim, klev, ksfc_type, idx_lnd, &! in
        & aa(:,:,:,imh:imqv), bb(:,:,ih:iqv),      &! in
        & aa_btm, bb_btm,                          &! inout
        & zen_h, zfn_h, zen_qv, zfn_qv             )! out
    END IF

    ! Set defaults
    zca(1:kproma,:) = 1._wp
    zcs(1:kproma,:) = 1._wp

    !===========================================================================
    ! Land surface
    !===========================================================================

    z0h_lnd(:)       = 0._wp
    q_snocpymlt(:)   = 0._wp

    IF (idx_lnd <= ksfc_type) THEN

      ! If land is present, JSBACH is currently the only surface scheme supported by ECHAM physcis package
#ifndef __NO_JSBACH__

      sat_surface_specific_humidity(:) = 0._wp
      dry_static_energy(:) = 0._wp

      ztsfc_lnd(:)        = 280._wp
      ztsfc_lnd_eff(:)    = 280._wp
      z0m_tile(:,idx_lnd) = 0._wp

      IF (phy_config%llake) THEN
        CALL jsbach_interface ( jg, nblock, 1, kproma, pdtime, psteplen,                   & ! in
          & t_air             = ptemp(1:kproma),                                           & ! in
          & q_air             = pq(1:kproma),                                              & ! in
          & rain              = prsfl(1:kproma) + prsfc(1:kproma),                         & ! in
          & snow              = pssfl(1:kproma) + pssfc(1:kproma),                         & ! in
          & wind_air          = SQRT(pu(1:kproma)**2 + pv(1:kproma)**2),                   & ! in
          & wind_10m          = SQRT(pu(1:kproma)**2 + pv(1:kproma)**2),                   & ! in, temporary
          & lw_srf_down       = plw_down(1:kproma),                                        & ! in
          & swvis_srf_down    = zswvis_down(1:kproma),                                     & ! in
          & swnir_srf_down    = zswnir_down(1:kproma),                                     & ! in
          & swpar_srf_down    = pswpar_down(1:kproma),                                     & ! in
          & frac_par_diffuse  = ppardff(1:kproma),                                         & ! in
          & press_srf         = presi_old(1:kproma,klev+1),                                & ! in
          & drag_srf          = pfac_sfc(1:kproma) * pcfh_tile(1:kproma,idx_lnd),          & ! in
          & t_acoef           = zen_h(1:kproma, idx_lnd),                                  & ! in
          & t_bcoef           = zfn_h(1:kproma, idx_lnd),                                  & ! in
          & q_acoef           = zen_qv(1:kproma, idx_lnd),                                 & ! in
          & q_bcoef           = zfn_qv(1:kproma, idx_lnd),                                 & ! in
          & pch               = MERGE(pch_tile(1:kproma,idx_lnd),1._wp,lsm(1:kproma)>0._wp),  & ! in
          & cos_zenith_angle  = pcosmu0(1:kproma),                                         & ! in
          & CO2_air           = SPREAD(mmr_co2, DIM=1, NCOPIES=kproma),                    & ! in
          & t_srf             = ztsfc_lnd(1:kproma),                                       & ! out (T_s^(n+1)) surface temp (filtered, if Asselin)
          & t_eff_srf         = ztsfc_lnd_eff(1:kproma),                                   & ! out (T_s^eff) surface temp (effective, for longwave rad)
          & qsat_srf          = sat_surface_specific_humidity(1:kproma),                   & ! out
          & s_srf             = dry_static_energy(1:kproma),                               & ! out (s_s^star, for vertical diffusion scheme)
          & fact_q_air        = pcair(1:kproma),                                           & ! out
          & fact_qsat_srf     = pcsat(1:kproma),                                           & ! out
          & evapotrans        = pevap_tile(1:kproma, idx_lnd),                             & ! out
          & latent_hflx       = plhflx_tile(1:kproma, idx_lnd),                            & ! out
          & sensible_hflx     = pshflx_tile(1:kproma, idx_lnd),                            & ! out
          & grnd_hflx         = zgrnd_hflx(1:kproma, idx_lnd),                             & ! out
          & grnd_hcap         = zgrnd_hcap(1:kproma, idx_lnd),                             & ! out
          & rough_h_srf       = z0h_lnd(1:kproma),                                         & ! out
          & rough_m_srf       = z0m_tile(1:kproma, idx_lnd),                               & ! out
          & q_snocpymlt       = q_snocpymlt(1:kproma),                                     & ! out
          & alb_vis_dir       = albvisdir_tile(1:kproma, idx_lnd),                         & ! out
          & alb_nir_dir       = albnirdir_tile(1:kproma, idx_lnd),                         & ! out
          & alb_vis_dif       = albvisdif_tile(1:kproma, idx_lnd),                         & ! out
          & alb_nir_dif       = albnirdif_tile(1:kproma, idx_lnd),                         & ! out
          !
          & drag_wtr          = pfac_sfc(1:kproma) * pcfh_tile(1:kproma,idx_wtr),          & ! in
          & drag_ice          = pfac_sfc(1:kproma) * pcfh_tile(1:kproma,idx_ice),          & ! in
          & t_acoef_wtr       = zen_h(1:kproma, idx_wtr),                                  & ! in
          & t_bcoef_wtr       = zfn_h(1:kproma, idx_wtr),                                  & ! in
          & q_acoef_wtr       = zen_qv(1:kproma, idx_wtr),                                 & ! in
          & q_bcoef_wtr       = zfn_qv(1:kproma, idx_wtr),                                 & ! in
          & t_acoef_ice       = zen_h(1:kproma, idx_ice),                                  & ! in
          & t_bcoef_ice       = zfn_h(1:kproma, idx_ice),                                  & ! in
          & q_acoef_ice       = zen_qv(1:kproma, idx_ice),                                 & ! in
          & q_bcoef_ice       = zfn_qv(1:kproma, idx_ice),                                 & ! in
          & albedo_lwtr       = albedo_tile(1:kproma, idx_wtr),                            & ! in
          & t_lwtr            = ztsfc_wtr(1:kproma),                                       & ! out
          & evapo_wtr         = zevap_wtr(1:kproma),                                       & ! out
          & latent_hflx_wtr   = zlhflx_wtr(1:kproma),                                      & ! out
          & sensible_hflx_wtr = zshflx_wtr(1:kproma),                                      & ! out
          & t_lice            = ztsfc_ice(1:kproma),                                       & ! out
          & evapo_ice         = zevap_ice(1:kproma),                                       & ! out
          & latent_hflx_ice   = zlhflx_ice(1:kproma),                                      & ! out
          & sensible_hflx_ice = zshflx_ice(1:kproma),                                      & ! out
          & albedo_lice       = zalbvisdir_ice(1:kproma),                                  & ! out
          & ice_fract_lake    = lake_ice_frc(1:kproma)                                     & ! out
          )
      ELSE
        CALL jsbach_interface ( jg, nblock, 1, kproma, pdtime, psteplen,                   & ! in
          & t_air             = ptemp(1:kproma),                                           & ! in
          & q_air             = pq(1:kproma),                                              & ! in
          & rain              = prsfl(1:kproma) + prsfc(1:kproma),                         & ! in
          & snow              = pssfl(1:kproma) + pssfc(1:kproma),                         & ! in
          & wind_air          = SQRT(pu(1:kproma)**2 + pv(1:kproma)**2),                   & ! in
          & wind_10m          = SQRT(pu(1:kproma)**2 + pv(1:kproma)**2),                   & ! in, temporary
          & lw_srf_down       = plw_down(1:kproma),                                        & ! in
          & swvis_srf_down    = zswvis_down(1:kproma),                                     & ! in
          & swnir_srf_down    = zswnir_down(1:kproma),                                     & ! in
          & swpar_srf_down    = pswpar_down(1:kproma),                                     & ! in
          & frac_par_diffuse  = ppardff(1:kproma),                                         & ! in
          & press_srf         = presi_old(1:kproma,klev+1),                                & ! in
          & drag_srf          = pfac_sfc(1:kproma) * pcfh_tile(1:kproma,idx_lnd),          & ! in
          & t_acoef           = zen_h(1:kproma, idx_lnd),                                  & ! in
          & t_bcoef           = zfn_h(1:kproma, idx_lnd),                                  & ! in
          & q_acoef           = zen_qv(1:kproma, idx_lnd),                                 & ! in
          & q_bcoef           = zfn_qv(1:kproma, idx_lnd),                                 & ! in
          & pch               = MERGE(pch_tile(1:kproma,idx_lnd),1._wp,lsm(1:kproma)>0._wp),  & ! in
          & cos_zenith_angle  = pcosmu0(1:kproma),                                         & ! in
          & CO2_air           = SPREAD(mmr_co2, DIM=1, NCOPIES=kproma),                    & ! in
          & t_srf             = ztsfc_lnd(1:kproma),                                       & ! out (T_s^(n+1)) surface temp (filtered, if Asselin)
          & t_eff_srf         = ztsfc_lnd_eff(1:kproma),                                   & ! out (T_s^eff) surface temp (effective, for longwave rad)
          & qsat_srf          = sat_surface_specific_humidity(1:kproma),                   & ! out
          & s_srf             = dry_static_energy(1:kproma),                               & ! out (s_s^star, for vertical diffusion scheme)
          & fact_q_air        = pcair(1:kproma),                                           & ! out
          & fact_qsat_srf     = pcsat(1:kproma),                                           & ! out
          & evapotrans        = pevap_tile(1:kproma, idx_lnd),                             & ! out
          & latent_hflx       = plhflx_tile(1:kproma, idx_lnd),                            & ! out
          & sensible_hflx     = pshflx_tile(1:kproma, idx_lnd),                            & ! out
          & grnd_hflx         = zgrnd_hflx(1:kproma, idx_lnd),                             & ! out
          & grnd_hcap         = zgrnd_hcap(1:kproma, idx_lnd),                             & ! out
          & rough_h_srf       = z0h_lnd(1:kproma),                                         & ! out
          & rough_m_srf       = z0m_tile(1:kproma, idx_lnd),                               & ! out
          & q_snocpymlt       = q_snocpymlt(1:kproma),                                     & ! out
          & alb_vis_dir       = albvisdir_tile(1:kproma, idx_lnd),                         & ! out
          & alb_nir_dir       = albnirdir_tile(1:kproma, idx_lnd),                         & ! out
          & alb_vis_dif       = albvisdif_tile(1:kproma, idx_lnd),                         & ! out
          & alb_nir_dif       = albnirdif_tile(1:kproma, idx_lnd)                          & ! out
          )
      END IF

      ptsfc_tile(1:kproma,idx_lnd) = ztsfc_lnd(1:kproma)
      IF (phy_config%llake) THEN
        WHERE (alake(1:kproma) > 0._wp)
          ptsfc_tile    (1:kproma, idx_wtr) = ztsfc_wtr     (1:kproma)
          ptsfc_tile    (1:kproma, idx_ice) = ztsfc_ice     (1:kproma)
          pevap_tile    (1:kproma, idx_wtr) = zevap_wtr     (1:kproma)
          plhflx_tile   (1:kproma, idx_wtr) = zlhflx_wtr    (1:kproma)
          pshflx_tile   (1:kproma, idx_wtr) = zshflx_wtr    (1:kproma)
          pevap_tile    (1:kproma, idx_ice) = zevap_ice     (1:kproma)
          plhflx_tile   (1:kproma, idx_ice) = zlhflx_ice    (1:kproma)
          pshflx_tile   (1:kproma, idx_ice) = zshflx_ice    (1:kproma)
          albvisdir_tile(1:kproma, idx_ice) = zalbvisdir_ice(1:kproma)
          albvisdif_tile(1:kproma, idx_ice) = zalbvisdir_ice(1:kproma)
          albnirdir_tile(1:kproma, idx_ice) = zalbvisdir_ice(1:kproma)
          albnirdif_tile(1:kproma, idx_ice) = zalbvisdir_ice(1:kproma)
        ELSEWHERE
          lake_ice_frc(1:kproma) = 0._wp
        ENDWHERE
      END IF
      pcpt_tile (1:kproma,idx_lnd) = dry_static_energy(1:kproma)
      pqsat_tile(1:kproma,idx_lnd) = sat_surface_specific_humidity(1:kproma)

      ! Set the evapotranspiration coefficients, to be used later in
      ! blending and in diagnosing surface fluxes.
      !
      zca(1:kproma,idx_lnd) = pcair(1:kproma)
      zcs(1:kproma,idx_lnd) = pcsat(1:kproma)

#else
      CALL finish(method_name, "The JSBACH component is not activated")
#endif
    END IF

    !===========================================================================
    ! Ocean model
    !===========================================================================
    IF (idx_wtr <= ksfc_type) THEN

#ifndef __NO_ICON_OCEAN__
      IF (phy_config%lmlo) THEN
        ztsfc_wtr(:) = ptsfc_tile(:,idx_wtr)
        CALL ml_ocean ( kbdim, 1, kproma, pdtime, & !< in
          & pahflw=plhflx_tile(:,idx_wtr),        & !< in
          & pahfsw=pshflx_tile(:,idx_wtr),        & !< in
          & ptrflw=plw(:),                        & !< in
          & psoflw=pswvis(:) + pswnir(:),         & !< in
          & ptsw=ztsfc_wtr(:) )                     ! inout
        WHERE (alake(1:kproma) < EPSILON(1._wp))
          ptsfc_tile(1:kproma, idx_wtr) = ztsfc_wtr(1:kproma)
        END WHERE
      END IF
#endif

      ! Albedo model for the ocean (including open water of lakes)
      ! TBD: This should be replaced by routine mo_surface_ocean:update_albedo_ocean from ECHAM6.2
      albvisdir_tile(1:kproma,idx_wtr) = albedoW
      albvisdif_tile(1:kproma,idx_wtr) = albedoW
      albnirdir_tile(1:kproma,idx_wtr) = albedoW
      albnirdif_tile(1:kproma,idx_wtr) = albedoW

    END IF

    !===========================================================================
    ! Sea-ice model (thermodynamic)
    !===========================================================================

    IF (idx_ice <= ksfc_type .AND. phy_config%lice) THEN

#ifndef __NO_ICON_OCEAN__
      ! LL This is a temporary solution,
      ! we should restrcure ice thermodynamics in a more stand-alone way

      ! For explicit coupling to ice:

      ! Freezing point of sea-water
      Tfw = Tf

      ! ECHAM has no tiles for SW & LW and this is how it's solved there
      ! Net shortwave on all bands.
      ! Net longwave - we don't have tiles yet
      ! First all ice classes
      DO k=1,kice
        swflx_ice(1:kproma,k) = &
          & zswvisdif_down(1:kproma) * (1._wp - albvisdif_ice(1:kproma,k)) + &
          & zswvisdir_down(1:kproma) * (1._wp - albvisdir_ice(1:kproma,k)) + &
          & zswnirdif_down(1:kproma) * (1._wp - albnirdif_ice(1:kproma,k)) + &
          & zswnirdir_down(1:kproma) * (1._wp - albnirdir_ice(1:kproma,k))

        nonsolar_ice(1:kproma,k) = &
          zemiss_def * (plw_down(1:kproma) - stbo * (Tsurf(1:kproma,k)+tmelt)**4) &  ! longwave net
          & + plhflx_tile(1:kproma,idx_ice) + pshflx_tile(1:kproma,idx_ice)

        dnonsolardT(1:kproma,k) = -4._wp * zemiss_def * stbo * (Tsurf(1:kproma,k)+tmelt)**3

      ENDDO

      CALL ice_fast(1, kproma, kbdim, kice, pdtime, &
        &   Tsurf,              &
        &   T1,                 &
        &   T2,                 &
        &   hi,                 &
        &   hs,                 &
        &   Qtop,               &
        &   Qbot,               &
        &   swflx_ice,          &
        &   nonsolar_ice,       &
        &   dnonsolardT,        &
        &   Tfw,                &
        &   albvisdir_ice,      &
        &   albvisdif_ice,      &
        &   albnirdir_ice,      &
        &   albnirdif_ice )

      ! Let it snow and melt and grow
      !      DO k=1,kice
      !        WHERE ( hi(:,1) > 0._wp )
      !          hs(:,k) = hs(:,k) + (pssfl + pssfc)*pdtime/rhos &
      !            &   - MIN( Qtop(:,k)*pdtime/( alf*rhos ), hs(:,k) )
      !        ENDWHERE
      !        hi(:,k) = hi(:,k) - MIN( Qbot(:,k)*pdtime/( alf*rhoi ), hi(:,k) )
      !        WHERE ( hs(:,1) <= 0._wp )
      !          hi(:,k) = hi(:,k) - MIN( Qtop(:,k)*pdtime/( alf*rhoi ), hi(:,k) )
      !        ENDWHERE
      !      ENDDO
      !      hi(:,:) = max( hi(:,:), 0._wp )
      ! Let it snow in AMIP
      IF ( phy_config%lamip ) THEN
        DO k=1,kice
          ! Snowfall on ice - no ice => no snow
          WHERE ( hi(1:kproma,k) > 0._wp )
            ! Snow only falls when it's below freezing
            WHERE ( Tsurf(1:kproma,k) < 0._wp )
              hs(1:kproma,k) = hs(1:kproma,k) + (pssfl(1:kproma) + pssfc(1:kproma))*pdtime/rhos
            ENDWHERE
            ! Snow melt
            hs(1:kproma,k) = hs(1:kproma,k) - MIN( Qtop(1:kproma,k)*pdtime/( alf*rhos ), hs(1:kproma,k) )
          ELSEWHERE
            hs(1:kproma,k) = 0._wp
          ENDWHERE
        ENDDO
      ENDIF

      conc_sum(1:kproma) = SUM(conc(1:kproma,:),2)
      WHERE (alake(1:kproma) < EPSILON(1._wp))
        ! Average the albedo.
!        albvisdir_tile(1:kproma,idx_ice) = 0._wp
!        albvisdif_tile(1:kproma,idx_ice) = 0._wp
!        albnirdir_tile(1:kproma,idx_ice) = 0._wp
!        albnirdif_tile(1:kproma,idx_ice) = 0._wp
        WHERE (conc_sum(1:kproma) > 1.e-6_wp)
          albvisdir_tile(1:kproma,idx_ice) = SUM( conc(1:kproma,:) * albvisdir_ice(1:kproma,:), 2 ) / conc_sum(1:kproma)
          albvisdif_tile(1:kproma,idx_ice) = SUM( conc(1:kproma,:) * albvisdif_ice(1:kproma,:), 2 ) / conc_sum(1:kproma)
          albnirdir_tile(1:kproma,idx_ice) = SUM( conc(1:kproma,:) * albnirdir_ice(1:kproma,:), 2 ) / conc_sum(1:kproma)
          albnirdif_tile(1:kproma,idx_ice) = SUM( conc(1:kproma,:) * albnirdif_ice(1:kproma,:), 2 ) / conc_sum(1:kproma)

          ! Set the tile temperature
          ptsfc_tile(1:kproma,idx_ice) = Tsurf(1:kproma,1) + tmelt
        END WHERE
      END WHERE

      ! Compute new dry static energy
      ! (Switched off for now, should be used for implicit coupling)
      !pcpt_tile(1:kproma,idx_ice) = ptsfc_tile(1:kproma,idx_ice) * zt2s_conv(1:kproma,idx_ice)

#else
    ! __NO_ICON_OCEAN__
      CALL finish(method_name, "The ice process requires the ICON_OCEAN component")
#endif
    ENDIF ! lice

    !===================================================================
    ! AFTER CALLING land/ocean/ice model
    !===================================================================
    ! Turbulent transport of moisture and dry static energy:
    ! Get solution of the two variables on the lowest model level.
    !-------------------------------------------------------------------
    ! - Over individual tiles
    !   For echam developers: relationship to "update_surface" of echam6:
    !   bb_btm(:,jsfc,ih) : tpfac2*land%ztklevl, tpfac2*ice%ztklevi, tpfac2*ocean%ztklevw
    !   bb_btm(:,jsfc,iqv): tpfac2*land%zqklevl, tpfac2*ice%zqklevi, tpfac2*ocean%zqklevw

    DO jsfc = 1,ksfc_type
       bb_btm(1:kproma,jsfc,ih)  = tpfac2*(    zen_h (1:kproma,jsfc) &
                                 &         *pcpt_tile(1:kproma,jsfc) &
                                 &         +   zfn_h (1:kproma,jsfc) )

       bb_btm(1:kproma,jsfc,iqv) = tpfac2*(    zen_qv(1:kproma,jsfc) &
                                 &        *pqsat_tile(1:kproma,jsfc) &
                                 &        +    zfn_qv(1:kproma,jsfc) )
    END DO

    ! - Grid box mean
    !   For echam developers: relationship to "update_surface" of echam6:
    !   bb(:,klev,ih) : ztdif_new
    !   bb(:,klev,iqv): zqdif_new

     se_sum(1:kproma) = 0._wp    ! sum of weighted solution
     qv_sum(1:kproma) = 0._wp    ! sum of weighted solution
    wgt_sum(1:kproma) = 0._wp    ! sum of weights

    DO jsfc = 1,ksfc_type
           wgt(1:kproma) = pfrc(1:kproma,jsfc)
       wgt_sum(1:kproma) = wgt_sum(1:kproma) + wgt(1:kproma)
        se_sum(1:kproma) = se_sum(1:kproma) + bb_btm(1:kproma,jsfc,ih ) * wgt(1:kproma)
        qv_sum(1:kproma) = qv_sum(1:kproma) + bb_btm(1:kproma,jsfc,iqv) * wgt(1:kproma)
    ENDDO

    IF (lsfc_heat_flux) THEN
      bb(1:kproma,klev,ih ) = se_sum(1:kproma)/wgt_sum(1:kproma)
      bb(1:kproma,klev,iqv) = qv_sum(1:kproma)/wgt_sum(1:kproma)
    ELSE
      jsfc = 1
      bb(1:kproma,klev,ih ) = bb_btm(1:kproma,jsfc,ih )
      bb(1:kproma,klev,iqv) = bb_btm(1:kproma,jsfc,iqv)
    END IF

    !-------------------------------------------------------------------
    ! Turbulent transport of u and v: adjust the right-hand side vector,
    ! then perform the bottom level elimination to get the solution
    !-------------------------------------------------------------------
    ! Add additional terms to the r.h.s. of the velocity equations
    ! to take into account ocean currents.
    ! Note that in subroutine rhs_setup the constant tpfac2 has been
    ! multiplied to the r.h.s. array bb. Thus the additional terms here
    ! need to be scaled by the same factor.

    IF (idx_wtr.LE.ksfc_type) THEN   ! Open water is considered
      IF (idx_ice.LE.ksfc_type) THEN ! Sea ice is also considered
        zfrc_oce(1:kproma) = pfrc(1:kproma,idx_wtr)+pfrc(1:kproma,idx_ice)
      ELSE ! only open water
        zfrc_oce(1:kproma) = pfrc(1:kproma,idx_wtr)
      ENDIF
      bb(1:kproma,klev,iu) =   bb(1:kproma,klev,iu)                   &
                           & - pocu(1:kproma)*zfrc_oce(1:kproma)*tpfac2
      bb(1:kproma,klev,iv) =   bb(1:kproma,klev,iv)                   &
                           & - pocv(1:kproma)*zfrc_oce(1:kproma)*tpfac2
    ENDIF

    ! Bottom level elimination

    im   = imuv
    jk   = klev    ! Bottom level index
    jkm1 = jk - 1

    aa(1:kproma,jk,2,im) =  aa(1:kproma,jk,2,im)                      &
                         & -aa(1:kproma,jk,1,im)*aa(1:kproma,jkm1,3,im)
    aa(1:kproma,jk,3,im) =  aa(1:kproma,jk,3,im)/aa(1:kproma,jk,2,im)

    bb(1:kproma,jk,iu) = (bb(1:kproma,jk,iu)                         &
                       & -aa(1:kproma,jk,1,im)*bb(1:kproma,jkm1,iu)) &
                       & /aa(1:kproma,jk,2,im)

    bb(1:kproma,jk,iv) = (bb(1:kproma,jk,iv)                         &
                       & -aa(1:kproma,jk,1,im)*bb(1:kproma,jkm1,iv)) &
                       & /aa(1:kproma,jk,2,im)

   !-------------------------------------------------------------------
   ! Various diagnostics
   !-------------------------------------------------------------------

   CALL surface_fluxes( lsfc_heat_flux, psteplen,             &! in
                      & kproma, kbdim, ksfc_type,             &! in
                      & idx_wtr, idx_ice, idx_lnd, ih, iqv,   &! in
                      & pfrc, alake,                          &! in
                      & pcfh_tile, pfac_sfc,                  &! in
                      & pcpt_tile, pqsat_tile,                &! in
                      & zca, zcs, bb_btm(:,:,ih:iqv),         &! in
                      & plhflx_gbm, pshflx_gbm,               &! out
                      & pevap_gbm,                            &! out
                      & plhflx_tile, pshflx_tile,             &! inout
                      & dshflx_dT_tile,                       &! out
                      & pevap_tile                            &! inout
                      & )


    DO jsfc=1,ksfc_type
      zalbvis(:) = 0._wp
      WHERE(zswvis_down(1:kproma) > 0._wp)
        zalbvis(1:kproma) = &
          & (albvisdir_tile(1:kproma,jsfc) * zswvisdir_down(1:kproma) + albvisdif_tile(1:kproma,jsfc) * zswvisdif_down(1:kproma)) &
          & / zswvis_down(1:kproma)
      END WHERE
      zalbnir(:) = 0._wp
      WHERE(zswnir_down(1:kproma) > 0._wp)
        zalbnir(1:kproma) = &
          & (albnirdir_tile(1:kproma,jsfc) * zswnirdir_down(1:kproma) + albnirdif_tile(1:kproma,jsfc) * zswnirdif_down(1:kproma)) &
          & / zswnir_down(1:kproma)
      END WHERE
      WHERE(zsw_down(1:kproma) > 0._wp)
        albedo_tile(1:kproma,jsfc) = &
          & (zalbvis(1:kproma) * zswvis_down(1:kproma) + zalbnir(1:kproma) * zswnir_down(1:kproma)) &
          & / zsw_down(1:kproma)
      END WHERE
    END DO

    ! calculate grid box mean surface temperature
    ptsfc(:) = 0._wp
    DO jsfc=1,ksfc_type
      ptsfc(1:kproma) = ptsfc(1:kproma) + pfrc(1:kproma,jsfc) * ptsfc_tile(1:kproma,jsfc)
    ENDDO

    ! calculate grid box mean radiative temperature for use in radiation
    ptsfc_rad(:) = 0._wp
    DO jsfc=1,ksfc_type
      ptsfc_rad(1:kproma) = ptsfc_rad(1:kproma) + pfrc(1:kproma,jsfc) * ptsfc_tile(1:kproma,jsfc)**4
    ENDDO
    ptsfc_rad(1:kproma) = ptsfc_rad(1:kproma)**0.25_wp

    ! Compute lw and sw surface radiation fluxes on tiles
    DO jsfc=1,ksfc_type
      IF (jsfc == idx_lnd) THEN
        plwflx_tile(1:kproma,jsfc) = zemiss_def * (plw_down(1:kproma) - stbo * ztsfc_lnd_eff(1:kproma)**4)
      ELSE
        plwflx_tile(1:kproma,jsfc) = zemiss_def * (plw_down(1:kproma) - stbo * ptsfc_tile(1:kproma,jsfc)**4)
      END IF
    END DO

    DO jsfc=1,ksfc_type
      pswflx_tile(1:kproma,jsfc) = &
        & zswvisdif_down(1:kproma) * (1._wp - albvisdif_tile(1:kproma,jsfc)) + &
        & zswvisdir_down(1:kproma) * (1._wp - albvisdir_tile(1:kproma,jsfc)) + &
        & zswnirdif_down(1:kproma) * (1._wp - albnirdif_tile(1:kproma,jsfc)) + &
        & zswnirdir_down(1:kproma) * (1._wp - albnirdir_tile(1:kproma,jsfc))
    END DO

    ! Merge sw and lw surface fluxes
    ! This includes the update of the lw flux on land due to the new surface temperature where only part
    ! of the net radiation was used (due to the Taylor truncation in the surface energy balance)
    plw(:) = 0._wp
    psw(:) = 0._wp
    DO jsfc=1,ksfc_type
      plw(1:kproma) = plw(1:kproma) + pfrc(1:kproma,jsfc) * plwflx_tile(1:kproma,jsfc)
      psw(1:kproma) = psw(1:kproma) + pfrc(1:kproma,jsfc) * pswflx_tile(1:kproma,jsfc)
    END DO

    ! Merge surface albedos
    albvisdir(:) = 0._wp
    albvisdif(:) = 0._wp
    albnirdir(:) = 0._wp
    albnirdif(:) = 0._wp
    albedo   (:) = 0._wp
    DO jsfc=1,nsfc_type
      albvisdir(1:kproma) = albvisdir(1:kproma) + pfrc(1:kproma,jsfc) * albvisdir_tile(1:kproma,jsfc)
      albvisdif(1:kproma) = albvisdif(1:kproma) + pfrc(1:kproma,jsfc) * albvisdif_tile(1:kproma,jsfc)
      albnirdir(1:kproma) = albnirdir(1:kproma) + pfrc(1:kproma,jsfc) * albnirdir_tile(1:kproma,jsfc)
      albnirdif(1:kproma) = albnirdif(1:kproma) + pfrc(1:kproma,jsfc) * albnirdif_tile(1:kproma,jsfc)
      albedo   (1:kproma) = albedo   (1:kproma) + pfrc(1:kproma,jsfc) * albedo_tile   (1:kproma,jsfc)
    END DO

    ! Mask out tiled variables
    ! wtr and ice tiles only for uncoupled case ... coupled atmo/ocean runs yield different results if wtr/ice is
    ! masked out over land.
    DO jsfc=1,ksfc_type
      mask(:) = .FALSE.
      IF (jsfc == idx_lnd) mask(1:kproma) = pfrc(1:kproma,jsfc) == 0._wp
      IF (.NOT. is_coupled_run()) THEN
        IF (jsfc == idx_wtr .OR. jsfc == idx_ice) &
          & mask(1:kproma) = pfrc(1:kproma,idx_wtr) == 0._wp .AND. pfrc(1:kproma,idx_ice) == 0._wp &
            &          .AND. alake(1:kproma) < EPSILON(1._wp)
      END IF
      WHERE (mask(1:kproma))
!        ptsfc_tile     (1:kproma,jsfc) = cdimissval
!        pqsat_tile     (1:kproma,jsfc) = cdimissval
!        pswflx_tile    (1:kproma,jsfc) = cdimissval
!        plwflx_tile    (1:kproma,jsfc) = cdimissval
!        pevap_tile     (1:kproma,jsfc) = cdimissval
!        pshflx_tile    (1:kproma,jsfc) = cdimissval
!        plhflx_tile    (1:kproma,jsfc) = cdimissval
!        albedo_tile    (1:kproma,jsfc) = cdimissval
!        albvisdir_tile (1:kproma,jsfc) = cdimissval
!        albvisdif_tile (1:kproma,jsfc) = cdimissval
!        albnirdir_tile (1:kproma,jsfc) = cdimissval
!        albnirdif_tile (1:kproma,jsfc) = cdimissval
        pu_stress_tile (1:kproma,jsfc) = cdimissval
        pv_stress_tile (1:kproma,jsfc) = cdimissval
        dshflx_dT_tile (1:kproma,jsfc) = cdimissval
        z0m_tile       (1:kproma,jsfc) = cdimissval
      END WHERE
      IF (jsfc == idx_lnd) THEN
        WHERE (mask(1:kproma))
          z0h_lnd(1:kproma) = cdimissval
        END WHERE
      END IF
    END DO


    END SUBROUTINE update_surface
  !-------------

END MODULE mo_surface
