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
  USE mo_physical_constants,ONLY: grav, Tf, alf, albedoW, zemiss_def, stbo, tmelt, rhos!!$, rhoi
  USE mo_mpi_phy_config,    ONLY: mpi_phy_config
  USE mo_echam_phy_memory,  ONLY: cdimissval
  USE mo_echam_vdiff_params,ONLY: tpfac2
  USE mo_vdiff_config,      ONLY: vdiff_config
  USE mo_vdiff_solver,      ONLY: ih, iqv, iu, iv, imh, imqv, imuv, &
                                & nmatrix, nvar_vdiff,              &
                                & matrix_to_richtmyer_coeff
  USE mo_surface_diag,      ONLY: wind_stress, surface_fluxes
#ifndef __NO_JSBACH__
  USE mo_jsb_interface,     ONLY: jsbach_interface
  USE mo_radiation_config,  ONLY: mmr_co2      ! This should be here only temporarily
#endif
  USE mo_echam_sfc_indices, ONLY: nsfc_type
#ifndef __NO_ICON_OCEAN__
  USE mo_sea_ice,           ONLY: ice_fast
  USE mo_ml_ocean,          ONLY: ml_ocean
#endif
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: update_surface

CONTAINS
  !>
  !!
  !!
  SUBROUTINE update_surface( jg,                                &! in
                           & kproma, kbdim,                     &! in
                           & kice,                              &! in
                           & klev, ksfc_type,                   &! in
                           & idx_wtr, idx_ice, idx_lnd,         &! in
                           & pdtime,                            &! in
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
                           & pevap_tile,                        &! out
                           & pco2nat,                           &! out
                           !! optional
                           & nblock,                            &! in
                           & lsm,                               &! in
                           & alake,                             &! in
                           & pu,                                &! in
                           & pv,                                &! in
                           & ptemp,                             &! in
                           & pq,                                &! in
                           & pco2,                              &! in
                           & prsfl,                             &! in
                           & prsfc,                             &! in
                           & pssfl,                             &! in
                           & pssfc,                             &! in
                           & rlds,                              &! in
                           & rlus,                              &! inout
                           & rsds,                              &! in
                           & rsus,                              &! in
                           !
                           & rvds_dir,                          &! in
                           & rpds_dir,                          &! in
                           & rnds_dir,                          &! in
                           & rvds_dif,                          &! in
                           & rpds_dif,                          &! in
                           & rnds_dif,                          &! in
                           !
                           & ps,                                &! in
                           & pcosmu0,                           &! in
                           & pch_tile,                          &! in
                           !! for JSBACH
                           & pcsat,                             &! inout
                           & pcair,                             &! inout
                           & q_snocpymlt,                       &! out
                           !
                           & z0m_tile, z0h_lnd,                 &! out
                           & albvisdir, albnirdir, albvisdif, albnirdif, &! inout
                           & albvisdir_tile,                    &! inout
                           & albnirdir_tile,                    &! inout
                           & albvisdif_tile,                    &! inout
                           & albnirdif_tile,                    &! inout
                           & albedo, albedo_tile,               &! inout
                           & pco2_flux_tile,                    &! inout
                           & ptsfc,                             &! out
                           & ptsfc_rad,                         &! out
                           & rsns_tile, rlns_tile,              &! out
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

    REAL(wp),INTENT(IN) :: pdtime
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

    REAL(wp),INTENT(OUT)   :: pu_stress_gbm (kbdim)
    REAL(wp),INTENT(OUT)   :: pv_stress_gbm (kbdim)
    REAL(wp),INTENT(OUT)   ::    plhflx_gbm (kbdim)
    REAL(wp),INTENT(OUT)   ::    pshflx_gbm (kbdim)
    REAL(wp),INTENT(OUT)   ::     pevap_gbm (kbdim)

    REAL(wp),INTENT(OUT)   :: pu_stress_tile (kbdim,ksfc_type)
    REAL(wp),INTENT(OUT)   :: pv_stress_tile (kbdim,ksfc_type)
    REAL(wp),INTENT(INOUT) :: plhflx_tile (kbdim,ksfc_type)   ! OUT
    REAL(wp),INTENT(INOUT) :: pshflx_tile (kbdim,ksfc_type)   ! OUT
    REAL(wp),INTENT(OUT)   :: pevap_tile (kbdim,ksfc_type)

    !! JSBACH input
    INTEGER, OPTIONAL,INTENT(IN) :: nblock
    REAL(wp),OPTIONAL,INTENT(IN) :: lsm(kbdim)
    REAL(wp),OPTIONAL,INTENT(IN) :: alake(kbdim)
    REAL(wp),OPTIONAL,INTENT(IN) :: pu        (kbdim)              ! zonal wind lowest level
    REAL(wp),OPTIONAL,INTENT(IN) :: pv        (kbdim)              ! meridional wind lowest level
    REAL(wp),OPTIONAL,INTENT(IN) :: ptemp     (kbdim)              ! temperature of lowest atmospheric level
    REAL(wp),OPTIONAL,INTENT(IN) :: pq        (kbdim)              ! humidity of lowest atmospheric level
    REAL(wp),OPTIONAL,INTENT(IN) :: pco2      (kbdim)              ! co2 of lowest atmospheric level
    REAL(wp),OPTIONAL,INTENT(IN) :: prsfl     (kbdim)              ! rain large scale
    REAL(wp),OPTIONAL,INTENT(IN) :: prsfc     (kbdim)              ! rain convective
    REAL(wp),OPTIONAL,INTENT(IN) :: pssfl     (kbdim)              ! snow large scale
    REAL(wp),OPTIONAL,INTENT(IN) :: pssfc     (kbdim)              ! snow convective
    REAL(wp),OPTIONAL,INTENT(IN) :: rlds      (kbdim)              ! downward surface  longwave flux [W/m2]
    REAL(wp),OPTIONAL,INTENT(IN) :: rsds      (kbdim)              ! downward surface shortwave flux [W/m2]
    
    REAL(wp),INTENT(IN) :: rvds_dir(kbdim)        ! all-sky   vis. dir. downward flux at current   time [W/m2]
    REAL(wp),INTENT(IN) :: rpds_dir(kbdim)        ! all-sky   par  dir. downward flux at current   time [W/m2]
    REAL(wp),INTENT(IN) :: rnds_dir(kbdim)        ! all-sky   nir  dir. downward flux at current   time [W/m2]
    REAL(wp),INTENT(IN) :: rvds_dif(kbdim)        ! all-sky   vis. dif. downward flux at current   time [W/m2]
    REAL(wp),INTENT(IN) :: rpds_dif(kbdim)        ! all-sky   par  dif. downward flux at current   time [W/m2]
    REAL(wp),INTENT(IN) :: rnds_dif(kbdim)        ! all-sky   nir  dif. downward flux at current   time [W/m2]

    REAL(wp),OPTIONAL,INTENT(IN) :: ps        (kbdim)              ! surface pressure
    REAL(wp),OPTIONAL,INTENT(IN) :: pcosmu0   (kbdim)              ! cos of zenith angle
    REAL(wp),OPTIONAL,INTENT(IN) :: pch_tile  (kbdim,ksfc_type)
    !! JSBACH output
    REAL(wp),OPTIONAL,INTENT(INOUT) :: pcsat(kbdim)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: pcair(kbdim)
    REAL(wp),OPTIONAL,INTENT(OUT)   :: q_snocpymlt(kbdim)
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
    REAL(wp),OPTIONAL,INTENT(OUT)   :: pco2nat  (kbdim)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: pco2_flux_tile(kbdim,ksfc_type)
    REAL(wp),OPTIONAL,INTENT(OUT)   :: ptsfc    (kbdim)
    REAL(wp),OPTIONAL,INTENT(OUT)   :: ptsfc_rad(kbdim)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: rlus     (kbdim)           ! INOUT upward surface  longwave flux [W/m2]
    REAL(wp),OPTIONAL,INTENT(IN)    :: rsus     (kbdim)           ! IN upward surface shortwave flux [W/m2]
    REAL(wp),OPTIONAL,INTENT(INOUT) :: rsns_tile(kbdim,ksfc_type) ! shortwave net flux at surface on tiles
    REAL(wp),OPTIONAL,INTENT(INOUT) :: rlns_tile(kbdim,ksfc_type) ! longwave net flux at surface on tiles
    REAL(wp),OPTIONAL,INTENT(OUT)   :: lake_ice_frc(kbdim)        ! fraction of ice on lakes
    !! Sea ice
    INTEGER,          INTENT(IN)    :: kice ! Number of ice thickness classes
    REAL(wp),OPTIONAL,INTENT(INOUT) :: Tsurf(kbdim,kice)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: T1   (kbdim,kice) ! for coupled ocean only
    REAL(wp),OPTIONAL,INTENT(INOUT) :: T2   (kbdim,kice) ! for coupled ocean only
    REAL(wp),OPTIONAL,INTENT(INOUT) :: hi   (kbdim,kice) ! for coupled ocean only
    REAL(wp),OPTIONAL,INTENT(INOUT) :: hs   (kbdim,kice) ! for coupled ocean only
    REAL(wp),OPTIONAL,INTENT(INOUT) :: Qtop (kbdim,kice) ! OUT
    REAL(wp),OPTIONAL,INTENT(INOUT) :: Qbot (kbdim,kice) ! OUT
    REAL(wp),OPTIONAL,INTENT(IN)    :: conc (kbdim,kice) ! for coupled ocean only
    REAL(wp),OPTIONAL,INTENT(INOUT) :: albvisdir_ice(kbdim,kice)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: albvisdif_ice(kbdim,kice)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: albnirdir_ice(kbdim,kice)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: albnirdif_ice(kbdim,kice)

! locals

    INTEGER  :: loidx  (kbdim,ksfc_type) !< counter for masks
    INTEGER  :: is     (ksfc_type)       !< counter for masks

    INTEGER  :: jsfc, jk, jkm1, im, k, jl, jls, js
    REAL(wp) :: se_sum(kbdim), qv_sum(kbdim), wgt_sum(kbdim), wgt(kbdim)
    REAL(wp) :: zca(kbdim,ksfc_type), zcs(kbdim,ksfc_type)
    REAL(wp) :: zfrc_oce(kbdim)

    REAL(wp) :: zen_h (kbdim,ksfc_type)
    REAL(wp) :: zfn_h (kbdim,ksfc_type)
    REAL(wp) :: zen_qv(kbdim,ksfc_type)
    REAL(wp) :: zfn_qv(kbdim,ksfc_type)

    REAL(wp) ::                                                    &
      & zlhflx_lnd(kbdim), zlhflx_lwtr(kbdim), zlhflx_lice(kbdim), &
      & zshflx_lnd(kbdim), zshflx_lwtr(kbdim), zshflx_lice(kbdim), &
      & zevap_lnd(kbdim), zevap_lwtr(kbdim), zevap_lice(kbdim),    &
      & sat_surface_specific_humidity(kbdim),                      &
      & dry_static_energy(kbdim),                                  &
      & ztsfc_lnd(kbdim), ztsfc_lnd_eff(kbdim),                    &
      & ztsfc_wtr(kbdim), ztsfc_lwtr(kbdim), ztsfc_lice(kbdim),    &
      & rvds(kbdim), rnds(kbdim), rpds(kbdim),                     &
      & rsns(kbdim), rlns(kbdim), frac_par_diffuse(kbdim),         &
      & zalbvis(kbdim), zalbnir(kbdim),                            &
      & zalbedo_lwtr(kbdim), zalbedo_lice(kbdim)

    REAL(wp) :: zgrnd_hflx(kbdim,ksfc_type), zgrnd_hcap(kbdim,ksfc_type)

    !REAL(wp) :: zt2s_conv(kbdim,ksfc_type)

    ! Sea ice
    REAL(wp) :: Tfw(kbdim)
    REAL(wp) :: swflx_ice(kbdim,kice), nonsolar_ice(kbdim,kice), dnonsolardT(kbdim,kice), conc_sum(kbdim)

    LOGICAL :: mask(kbdim)

   CHARACTER(len=*), PARAMETER :: method_name='mo_surface:update_surface'

    ! check for masks
    !
    DO jsfc = 1,ksfc_type
      is(jsfc) = 0
      DO jl = 1,kproma
        IF(pfrc(jl,jsfc).GT.0.0_wp) THEN
          is(jsfc) = is(jsfc) + 1
          loidx(is(jsfc),jsfc) = jl
        ENDIF
      ENDDO
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

    IF (vdiff_config%lsfc_mom_flux) THEN
       CALL wind_stress( kproma, kbdim, ksfc_type,            &! in
            &            pdtime,                              &! in
            &            pfrc, pcfm_tile, pfac_sfc,           &! in
            &            bb(:,klev,iu), bb(:,klev,iv),        &! in
            &            pu_stress_gbm,  pv_stress_gbm,       &! out
            &            pu_stress_tile, pv_stress_tile       )! out
    ELSE
       pu_stress_tile(:,:) = 0._wp
       pv_stress_tile(:,:) = 0._wp
       pu_stress_gbm (:)   = 0._wp
       pv_stress_gbm (:)   = 0._wp
    END IF

    ! Compute downward shortwave surface fluxes
    rvds(1:kproma)      = rvds_dif(1:kproma) + rvds_dir(1:kproma)
    rnds(1:kproma)      = rnds_dif(1:kproma) + rnds_dir(1:kproma)
    rpds(1:kproma)      = rpds_dif(1:kproma) + rpds_dir(1:kproma)

    ! Turbulent transport of moisture:
    ! - finish matrix set up;
    ! - perform bottom level elimination;
    ! - convert matrix entries to Richtmyer-Morton coefficients
    IF (idx_lnd <= ksfc_type) THEN
      CALL matrix_to_richtmyer_coeff( jg, kproma, kbdim, klev, ksfc_type, idx_lnd, &! in
        & aa(:,:,:,imh:imqv), bb(:,:,ih:iqv),      &! in
        & aa_btm, bb_btm,                          &! inout
        & zen_h, zfn_h, zen_qv, zfn_qv,            &! out
        & pcair = pcair(:),                        &! in
        & pcsat = pcsat(:))                         ! in
    ELSE
      CALL matrix_to_richtmyer_coeff( jg, kproma, kbdim, klev, ksfc_type, idx_lnd, &! in
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
    
    zlhflx_lnd(:)    = 0._wp
    zlhflx_lwtr(:)   = 0._wp
    zlhflx_lice(:)   = 0._wp
    zshflx_lnd(:)    = 0._wp
    zshflx_lwtr(:)   = 0._wp
    zshflx_lice(:)   = 0._wp
    zevap_lnd(:)     = 0._wp
    zevap_lwtr(:)    = 0._wp
    zevap_lice(:)    = 0._wp
    z0h_lnd(:)       = 0._wp
    q_snocpymlt(:)   = 0._wp
    lake_ice_frc(:)  = 0._wp

    IF (idx_lnd <= ksfc_type) THEN

      ! If land is present, JSBACH is currently the only surface scheme supported by ECHAM physcis package
#ifndef __NO_JSBACH__

      sat_surface_specific_humidity(:) = 0._wp
      dry_static_energy(:) = 0._wp
      ztsfc_lnd(:)         = 0._wp
      ztsfc_lnd_eff(:)     = 0._wp
      ztsfc_lwtr(:)        = 0._wp
      ztsfc_lice(:)        = 0._wp
      z0m_tile(:,idx_lnd)  = 0._wp

      WHERE (rpds(1:kproma) > 0._wp)
        frac_par_diffuse(1:kproma) = rpds_dif(1:kproma) / rpds(1:kproma)
      ELSE WHERE
        frac_par_diffuse(1:kproma) = 0._wp
      END WHERE

      IF (mpi_phy_config(jg)%llake) THEN
        CALL jsbach_interface ( jg, nblock, 1, kproma, pdtime, pdtime,                     & ! in
          & t_air             = ptemp(1:kproma),                                           & ! in
          & q_air             = pq(1:kproma),                                              & ! in
          & rain              = prsfl(1:kproma) + prsfc(1:kproma),                         & ! in
          & snow              = pssfl(1:kproma) + pssfc(1:kproma),                         & ! in
          & wind_air          = SQRT(pu(1:kproma)**2 + pv(1:kproma)**2),                   & ! in
          ! @todo: use real 10m wind
          & wind_10m          = SQRT(pu(1:kproma)**2 + pv(1:kproma)**2),                   & ! in, temporary
          & lw_srf_down       = rlds(1:kproma),                                            & ! in
          & swvis_srf_down    = rvds(1:kproma),                                            & ! in
          & swnir_srf_down    = rnds(1:kproma),                                            & ! in
          & swpar_srf_down    = rpds(1:kproma),                                            & ! in
          & frac_par_diffuse  = frac_par_diffuse(1:kproma),                                & ! in
          & press_srf         = ps(1:kproma),                                              & ! in
          & drag_srf          = grav*pfac_sfc(1:kproma) * pcfh_tile(1:kproma,idx_lnd),     & ! in
          & t_acoef           = zen_h(1:kproma, idx_lnd),                                  & ! in
          & t_bcoef           = zfn_h(1:kproma, idx_lnd),                                  & ! in
          & q_acoef           = zen_qv(1:kproma, idx_lnd),                                 & ! in
          & q_bcoef           = zfn_qv(1:kproma, idx_lnd),                                 & ! in
          & pch               = MERGE(pch_tile(1:kproma,idx_lnd),1._wp,lsm(1:kproma)>0._wp),  & ! in
          & cos_zenith_angle  = pcosmu0(1:kproma),                                         & ! in
          & CO2_air           = SPREAD(mmr_co2, DIM=1, NCOPIES=kproma),                    & ! in
          ! & CO2_air           = pco2(1:kproma),                                            & ! in
          & t_srf             = ztsfc_lnd(1:kproma),                                       & ! out (T_s^(n+1)) surface temp (filtered, if Asselin)
                                                                                             ! (filtered, if Asselin)
          & t_eff_srf         = ztsfc_lnd_eff(1:kproma),                                   & ! out (T_s^eff) surface temp (effective, for longwave rad)
                                                                                             ! (effective, for longwave rad)
          & qsat_srf          = sat_surface_specific_humidity(1:kproma),                   & ! out
          & s_srf             = dry_static_energy(1:kproma),                               & ! out (s_s^star, for vertical diffusion scheme)
          & fact_q_air        = pcair(1:kproma),                                           & ! out
          & fact_qsat_srf     = pcsat(1:kproma),                                           & ! out
          & evapotrans        = zevap_lnd(1:kproma),                                       & ! out
          & latent_hflx       = zlhflx_lnd(1:kproma),                                      & ! out
          & sensible_hflx     = zshflx_lnd(1:kproma),                                      & ! out
          & grnd_hflx         = zgrnd_hflx(1:kproma, idx_lnd),                             & ! out
          & grnd_hcap         = zgrnd_hcap(1:kproma, idx_lnd),                             & ! out
          & rough_h_srf       = z0h_lnd(1:kproma),                                         & ! out
          & rough_m_srf       = z0m_tile(1:kproma, idx_lnd),                               & ! out
          & q_snocpymlt       = q_snocpymlt(1:kproma),                                     & ! out
          & alb_vis_dir       = albvisdir_tile(1:kproma, idx_lnd),                         & ! out
          & alb_nir_dir       = albnirdir_tile(1:kproma, idx_lnd),                         & ! out
          & alb_vis_dif       = albvisdif_tile(1:kproma, idx_lnd),                         & ! out
          & alb_nir_dif       = albnirdif_tile(1:kproma, idx_lnd),                         & ! out
          & co2_flux          = pco2_flux_tile(1:kproma, idx_lnd),                         & ! out
          !
          & drag_wtr          = grav*pfac_sfc(1:kproma) * pcfh_tile(1:kproma,idx_wtr),     & ! in
          & drag_ice          = grav*pfac_sfc(1:kproma) * pcfh_tile(1:kproma,idx_ice),     & ! in
          & t_acoef_wtr       = zen_h(1:kproma, idx_wtr),                                  & ! in
          & t_bcoef_wtr       = zfn_h(1:kproma, idx_wtr),                                  & ! in
          & q_acoef_wtr       = zen_qv(1:kproma, idx_wtr),                                 & ! in
          & q_bcoef_wtr       = zfn_qv(1:kproma, idx_wtr),                                 & ! in
          & t_acoef_ice       = zen_h(1:kproma, idx_ice),                                  & ! in
          & t_bcoef_ice       = zfn_h(1:kproma, idx_ice),                                  & ! in
          & q_acoef_ice       = zen_qv(1:kproma, idx_ice),                                 & ! in
          & q_bcoef_ice       = zfn_qv(1:kproma, idx_ice),                                 & ! in
          & t_lwtr            = ztsfc_lwtr(1:kproma),                                      & ! out
          & evapo_wtr         = zevap_lwtr(1:kproma),                                      & ! out
          & latent_hflx_wtr   = zlhflx_lwtr(1:kproma),                                     & ! out
          & sensible_hflx_wtr = zshflx_lwtr(1:kproma),                                     & ! out
          & albedo_lwtr       = zalbedo_lwtr(1:kproma),                                    & ! out
          & t_lice            = ztsfc_lice(1:kproma),                                      & ! out
          & evapo_ice         = zevap_lice(1:kproma),                                      & ! out
          & latent_hflx_ice   = zlhflx_lice(1:kproma),                                     & ! out
          & sensible_hflx_ice = zshflx_lice(1:kproma),                                     & ! out
          & albedo_lice       = zalbedo_lice(1:kproma),                                    & ! out
          & ice_fract_lake    = lake_ice_frc(1:kproma)                                     & ! out
          )
      ELSE
        CALL jsbach_interface ( jg, nblock, 1, kproma, pdtime, pdtime,                    & ! in
          & t_air            = ptemp(1:kproma),                                           & ! in
          & q_air            = pq(1:kproma),                                              & ! in
          & rain             = prsfl(1:kproma) + prsfc(1:kproma),                         & ! in
          & snow             = pssfl(1:kproma) + pssfc(1:kproma),                         & ! in
          & wind_air         = SQRT(pu(1:kproma)**2 + pv(1:kproma)**2),                   & ! in
          ! @todo: use real 10m wind
          & wind_10m         = SQRT(pu(1:kproma)**2 + pv(1:kproma)**2),                   & ! in, temporary
          & lw_srf_down      = rlds(1:kproma),                                            & ! in
          & swvis_srf_down   = rvds(1:kproma),                                            & ! in
          & swnir_srf_down   = rnds(1:kproma),                                            & ! in
          & swpar_srf_down   = rpds(1:kproma),                                            & ! in
          & frac_par_diffuse = frac_par_diffuse(1:kproma),                                & ! in
          & press_srf        = ps(1:kproma),                                              & ! in
          & drag_srf         = grav*pfac_sfc(1:kproma) * pcfh_tile(1:kproma,idx_lnd),     & ! in
          & t_acoef          = zen_h(1:kproma, idx_lnd),                                  & ! in
          & t_bcoef          = zfn_h(1:kproma, idx_lnd),                                  & ! in
          & q_acoef          = zen_qv(1:kproma, idx_lnd),                                 & ! in
          & q_bcoef          = zfn_qv(1:kproma, idx_lnd),                                 & ! in
          & pch              = MERGE(pch_tile(1:kproma,idx_lnd),1._wp,lsm(1:kproma)>0._wp),  & ! in
          & cos_zenith_angle = pcosmu0(1:kproma),                                         & ! in
          & CO2_air          = SPREAD(mmr_co2, DIM=1, NCOPIES=kproma),                    & ! in
          ! & CO2_air          = pco2(1:kproma),                                            & ! in
          & t_srf            = ztsfc_lnd(1:kproma),                                       & ! out (T_s^(n+1)) surface temp 
                                                                                            ! (filtered, if Asselin)
          & t_eff_srf        = ztsfc_lnd_eff(1:kproma),                                   & ! out (T_s^eff) surface temp 
                                                                                            ! (effective, for longwave rad)
          & qsat_srf         = sat_surface_specific_humidity(1:kproma),                   & ! out
          & s_srf            = dry_static_energy(1:kproma),                               & ! out (s_s^star, for vert. diff. scheme)
          & fact_q_air       = pcair(1:kproma),                                           & ! out
          & fact_qsat_srf    = pcsat(1:kproma),                                           & ! out
          & evapotrans       = zevap_lnd(1:kproma),                                       & ! out
          & latent_hflx      = zlhflx_lnd(1:kproma),                                      & ! out
          & sensible_hflx    = zshflx_lnd(1:kproma),                                      & ! out
          & grnd_hflx        = zgrnd_hflx(1:kproma, idx_lnd),                             & ! out
          & grnd_hcap        = zgrnd_hcap(1:kproma, idx_lnd),                             & ! out
          & rough_h_srf      = z0h_lnd(1:kproma),                                         & ! out
          & rough_m_srf      = z0m_tile(1:kproma, idx_lnd),                               & ! out
          & q_snocpymlt      = q_snocpymlt(1:kproma),                                     & ! out
          & alb_vis_dir      = albvisdir_tile(1:kproma, idx_lnd),                         & ! out
          & alb_nir_dir      = albnirdir_tile(1:kproma, idx_lnd),                         & ! out
          & alb_vis_dif      = albvisdif_tile(1:kproma, idx_lnd),                         & ! out
          & alb_nir_dif      = albnirdif_tile(1:kproma, idx_lnd),                         & ! out
          & co2_flux         = pco2_flux_tile(1:kproma, idx_lnd)                          & ! out
        )
      END IF

      ! preliminary, dummy values
      pco2_flux_tile(1:kproma, idx_ice) =  0._wp

      ptsfc_tile(1:kproma,idx_lnd) = ztsfc_lnd(1:kproma)
      pcpt_tile (1:kproma,idx_lnd) = dry_static_energy(1:kproma)
      pqsat_tile(1:kproma,idx_lnd) = sat_surface_specific_humidity(1:kproma)
      IF (mpi_phy_config(jg)%llake) THEN
        IF (idx_wtr <= ksfc_type) THEN
          WHERE (alake(1:kproma) > 0._wp)
            ptsfc_tile    (1:kproma, idx_wtr) = ztsfc_lwtr   (1:kproma)
            albvisdir_tile(1:kproma, idx_wtr) = zalbedo_lwtr (1:kproma)
            albvisdif_tile(1:kproma, idx_wtr) = zalbedo_lwtr (1:kproma)
            albnirdir_tile(1:kproma, idx_wtr) = zalbedo_lwtr (1:kproma)
            albnirdif_tile(1:kproma, idx_wtr) = zalbedo_lwtr (1:kproma)
          END WHERE
        END IF
        IF (idx_ice <= ksfc_type) THEN
          WHERE (alake(1:kproma) > 0._wp)
            ptsfc_tile    (1:kproma, idx_ice) = ztsfc_lice   (1:kproma)
            albvisdir_tile(1:kproma, idx_ice) = zalbedo_lice (1:kproma)
            albvisdif_tile(1:kproma, idx_ice) = zalbedo_lice (1:kproma)
            albnirdir_tile(1:kproma, idx_ice) = zalbedo_lice (1:kproma)
            albnirdif_tile(1:kproma, idx_ice) = zalbedo_lice (1:kproma)
          ELSEWHERE
            lake_ice_frc(1:kproma) = 0._wp
          ENDWHERE
        END IF
      END IF

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

      rsns(1:kproma)      = rsds(1:kproma) - rsus(1:kproma)
      rlns(1:kproma)      = rlds(1:kproma) - rlus(1:kproma)

      IF (mpi_phy_config(jg)%lmlo) THEN
        CALL ml_ocean ( kbdim, 1, kproma, pdtime, &
          & pahflw=plhflx_tile(:,idx_wtr),        & ! dependency on kproma has to be checked
          & pahfsw=pshflx_tile(:,idx_wtr),        & ! dependency on kproma has to be checked
          & ptrflw=rlns(:),                       &
          & psoflw=rsns(:),                       &
          & ptsw=ztsfc_wtr(:) )                     ! out
        WHERE (alake(1:kproma) < EPSILON(1._wp))
          ptsfc_tile(1:kproma, idx_wtr) = ztsfc_wtr(1:kproma)
        END WHERE
      END IF
#endif

      ! Albedo model for the ocean
      ! TBD: This should be replaced by routine mo_surface_ocean:update_albedo_ocean from ECHAM6.2
      WHERE (alake(1:kproma) < EPSILON(1._wp))
        albvisdir_tile(1:kproma,idx_wtr) = albedoW
        albvisdif_tile(1:kproma,idx_wtr) = albedoW
        albnirdir_tile(1:kproma,idx_wtr) = albedoW
        albnirdif_tile(1:kproma,idx_wtr) = albedoW
      END WHERE

    END IF

    !===========================================================================
    ! Sea-ice model (thermodynamic)
    !===========================================================================

    IF (idx_ice <= ksfc_type .AND. mpi_phy_config(jg)%lice) THEN

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
          & rvds_dif(1:kproma) * (1._wp - albvisdif_ice(1:kproma,k)) + &
          & rvds_dir(1:kproma) * (1._wp - albvisdir_ice(1:kproma,k)) + &
          & rnds_dif(1:kproma) * (1._wp - albnirdif_ice(1:kproma,k)) + &
          & rnds_dir(1:kproma) * (1._wp - albnirdir_ice(1:kproma,k))

        nonsolar_ice(1:kproma,k) = &
          zemiss_def * (rlds(1:kproma) - stbo * (Tsurf(1:kproma,k)+tmelt)**4) &  ! longwave net
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
      IF ( mpi_phy_config(jg)%lamip ) THEN
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

      ! Average the albedo.
      conc_sum(1:kproma) = SUM(conc(1:kproma,:),2)
      WHERE (alake(1:kproma) < EPSILON(1._wp))
        albvisdir_tile(1:kproma,idx_ice) = 0._wp
        albvisdif_tile(1:kproma,idx_ice) = 0._wp
        albnirdir_tile(1:kproma,idx_ice) = 0._wp
        albnirdif_tile(1:kproma,idx_ice) = 0._wp
        WHERE (conc_sum(1:kproma) > 1.e-6_wp)
          albvisdir_tile(1:kproma,idx_ice) = SUM( conc(1:kproma,:) * albvisdir_ice(1:kproma,:), 2 ) / conc_sum(1:kproma)
          albvisdif_tile(1:kproma,idx_ice) = SUM( conc(1:kproma,:) * albvisdif_ice(1:kproma,:), 2 ) / conc_sum(1:kproma)
          albnirdir_tile(1:kproma,idx_ice) = SUM( conc(1:kproma,:) * albnirdir_ice(1:kproma,:), 2 ) / conc_sum(1:kproma)
          albnirdif_tile(1:kproma,idx_ice) = SUM( conc(1:kproma,:) * albnirdif_ice(1:kproma,:), 2 ) / conc_sum(1:kproma)

          ! Set the tile temperature, convert back to K
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

    ! calculate grid box mean surface of co2
    pco2nat(:) = 0._wp
    DO jsfc=1,ksfc_type
      pco2nat(1:kproma) = pco2nat(1:kproma) + pfrc(1:kproma,jsfc) * pco2_flux_tile(1:kproma,jsfc)
    ENDDO
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

    IF (vdiff_config%lsfc_heat_flux) THEN
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

    IF (vdiff_config%lsfc_heat_flux) THEN
       CALL surface_fluxes( kproma, kbdim, ksfc_type,             &! in
            &               idx_wtr, idx_ice, idx_lnd, ih, iqv,   &! in
            &               pdtime,                               &! in
            &               pfrc, alake, pcfh_tile, pfac_sfc,     &! in
            &               pcpt_tile, pqsat_tile,                &! in
            &               zca, zcs, bb_btm(:,:,ih:iqv),         &! in
            &               zlhflx_lnd, zlhflx_lwtr, zlhflx_lice, &! in
            &               zshflx_lnd, zshflx_lwtr, zshflx_lice, &! in
            &               zevap_lnd, zevap_lwtr, zevap_lice,    &! in
            &               plhflx_gbm, pshflx_gbm,               &! out
            &               pevap_gbm,                            &! out
            &               plhflx_tile, pshflx_tile,             &! out
            &               pevap_tile )                           ! out
    ELSE
       plhflx_tile(:,:) = 0._wp
       pshflx_tile(:,:) = 0._wp
       pevap_tile (:,:) = 0._wp
       plhflx_gbm (:)   = 0._wp
       pshflx_gbm (:)   = 0._wp
       pevap_gbm  (:)   = 0._wp
    END IF

    DO jsfc=1,ksfc_type
      zalbvis(:) = 0._wp
      WHERE(rvds(1:kproma) > 0._wp)
        zalbvis(1:kproma) = &
          & (albvisdir_tile(1:kproma,jsfc) * rvds_dir(1:kproma) + albvisdif_tile(1:kproma,jsfc) * rvds_dif(1:kproma)) &
          & / rvds(1:kproma)
      END WHERE
      zalbnir(:) = 0._wp
      WHERE(rnds(1:kproma) > 0._wp)
        zalbnir(1:kproma) = &
          & (albnirdir_tile(1:kproma,jsfc) * rnds_dir(1:kproma) + albnirdif_tile(1:kproma,jsfc) * rnds_dif(1:kproma)) &
          & / rnds(1:kproma)
      END WHERE
      albedo_tile(:,jsfc) = 0._wp
      WHERE(rsds(1:kproma) > 0._wp)
        albedo_tile(1:kproma,jsfc) = &
          & (zalbvis(1:kproma) * rvds(1:kproma) + zalbnir(1:kproma) * rnds(1:kproma)) &
          & / rsds(1:kproma)
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
      DO jls = 1,is(jsfc)
        ! set index
        js=loidx(jls,jsfc)
        IF (jsfc == idx_lnd) THEN
          rlns_tile(js,jsfc) = zemiss_def * (rlds(js) - stbo * ztsfc_lnd_eff(js)**4)
        ELSE
          rlns_tile(js,jsfc) = zemiss_def * (rlds(js) - stbo * ptsfc_tile(js,jsfc)**4)
        END IF

        rsns_tile(js,jsfc) = rvds_dif(js) * (1._wp - albvisdif_tile(js,jsfc)) + &
                           & rvds_dir(js) * (1._wp - albvisdir_tile(js,jsfc)) + &
                           & rnds_dif(js) * (1._wp - albnirdif_tile(js,jsfc)) + &
                           & rnds_dir(js) * (1._wp - albnirdir_tile(js,jsfc))
      END DO
    END DO

    ! Merge sw and lw surface fluxes
    ! This includes the update of the lw flux on land due to the new surface temperature where only part
    ! of the net radiation was used (due to the Taylor truncation in the surface energy balance)
    rlns(:) = 0._wp
    DO jsfc=1,ksfc_type
      DO jls = 1,is(jsfc)
        ! set index
        js=loidx(jls,jsfc)
        rlns(js) = rlns(js) + pfrc(js,jsfc) * rlns_tile(js,jsfc)
      END DO
    END DO
    rlus(1:kproma) = rlds(1:kproma) -rlns(1:kproma)

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
    DO jsfc=1,ksfc_type
      mask(:) = .FALSE.
      !
      mask(1:kproma) = pfrc(1:kproma,jsfc) <= 0._wp
      !
      WHERE (mask(1:kproma))
        pqsat_tile     (1:kproma,jsfc) = cdimissval
        albedo_tile    (1:kproma,jsfc) = cdimissval
        albvisdir_tile (1:kproma,jsfc) = cdimissval
        albvisdif_tile (1:kproma,jsfc) = cdimissval
        albnirdir_tile (1:kproma,jsfc) = cdimissval
        albnirdif_tile (1:kproma,jsfc) = cdimissval
      !  rsns_tile      (1:kproma,jsfc) = cdimissval
      !  rlns_tile      (1:kproma,jsfc) = cdimissval
        pevap_tile     (1:kproma,jsfc) = cdimissval
      !  pshflx_tile    (1:kproma,jsfc) = cdimissval
      !  plhflx_tile    (1:kproma,jsfc) = cdimissval
        ptsfc_tile     (1:kproma,jsfc) = cdimissval
      END WHERE
      ! land only
      IF (jsfc == idx_lnd) THEN
        WHERE (mask(1:kproma))
          z0h_lnd        (1:kproma)      = cdimissval
          z0m_tile       (1:kproma,jsfc) = cdimissval
          rsns_tile      (1:kproma,jsfc) = cdimissval
          rlns_tile      (1:kproma,jsfc) = cdimissval
          pshflx_tile    (1:kproma,jsfc) = cdimissval
          plhflx_tile    (1:kproma,jsfc) = cdimissval
        END WHERE
      END IF
    END DO

    !----------------------------------------------------------------------------
    ! For consistency z0m_tile for ice is masked out here
    !----------------------------------------------------------------------------
    IF (idx_ice<=ksfc_type) THEN  ! ice surface exists in the simulation
      mask(1:kproma) = pfrc(1:kproma,idx_ice) == 0._wp
      WHERE (mask(1:kproma))
        z0m_tile(1:kproma,idx_ice) = cdimissval
      ELSEWHERE
        ! z0m for ice is not calculated yet, so in case the ice surface changes
        ! it is set to the initial value again
        z0m_tile(1:kproma,idx_ice) = 1.e-3_wp
      ENDWHERE
    ENDIF

  !---------------------------------------------------------------------------
  !
  !---------------------------------------------------------------------------

    END SUBROUTINE update_surface
  !-------------

END MODULE mo_surface
