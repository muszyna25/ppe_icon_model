!>
!!
!! @author Hui Wan, MPI-M
!!
!! @par Revision History
!! First version by Hui Wan (2011-08)
!!
!! @par Copyright
!! 2002-2011 by DWD and MPI-M
!! This software is provided for non-commercial use only.
!! See the LICENSE and the WARRANTY conditions.
!!
!! @par License
!! The use of ICON is hereby granted free of charge for an unlimited time,
!! provided the following rules are accepted and applied:
!! <ol>
!! <li> You may use or modify this code for your own non commercial and non
!!    violent purposes.
!! <li> The code may not be re-distributed without the consent of the authors.
!! <li> The copyright notice and statement of authorship must appear in all
!!    copies.
!! <li> You accept the warranty conditions (see WARRANTY).
!! <li> In case you intend to use the code commercially, we oblige you to sign
!!    an according license agreement with DWD and MPI-M.
!! </ol>
!!
!! @par Warranty
!! This code has been tested up to a certain level. Defects and weaknesses,
!! which may be included in the code, do not establish any warranties by the
!! authors.
!! The authors do not make any warranty, express or implied, or assume any
!! liability or responsibility for the use, acquisition or application of this
!! software.
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
#ifndef __NO_JSBACH__
  USE mo_jsb_interface,     ONLY: jsbach_interface
#endif
  USE mo_icoham_sfc_indices,ONLY: nsfc_type
#ifndef __NO_ICON_OCEAN__
  USE mo_sea_ice,           ONLY: ice_fast
#endif
  USE mo_physical_constants,ONLY: rhos, rhoi, Tf, alf, albedoW, zemiss_def, stbo, tmelt
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: update_surface

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

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
                           & pu,                                &! in
                           & pv,                                &! in
                           & ptemp,                             &! in
                           & pq,                                &! in
                           & prsfl,                             &! in
                           & prsfc,                             &! in
                           & pssfl,                             &! in
                           & pssfc,                             &! in
                           & pemterall,                         &! in
                           & ptrsolall,                         &! in
                           & presi_old,                         &! in
                           & pcosmu0,                           &! in
                           & pch_tile,                          &! in
                           !! for JSBACH
                           & pcsat,                             &! inout
                           & pcair,                             &! inout
                           & zhsoil,                            &! out
                           & tte_corr,                          &! out
                           & z0h_lnd, z0m_lnd,                  &! out
                           & albvisdir,                         &! inout
                           & albnirdir,                         &! inout
                           & albvisdif,                         &! inout
                           & albnirdif,                         &! inout
                           & surface_temperature_rad,           &! out
                           & surface_temperature_eff,           &! out
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
                           & albnirdir_ice, albnirdif_ice,      &! inout
                           & albvisdir_wtr, albvisdif_wtr,      &! inout
                           & albnirdir_wtr, albnirdif_wtr,      &! inout
                           & pswflx_tile, plwflx_tile)             ! out (for coupling)

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
    REAL(wp),OPTIONAL,INTENT(IN) :: pu        (kbdim)              ! zonal wind lowest level
    REAL(wp),OPTIONAL,INTENT(IN) :: pv        (kbdim)              ! meridional wind lowest level
    REAL(wp),OPTIONAL,INTENT(IN) :: ptemp     (kbdim)              ! temperature of lowest atmospheric level
    REAL(wp),OPTIONAL,INTENT(IN) :: pq        (kbdim)              ! humidity of lowest atmospheric level
    REAL(wp),OPTIONAL,INTENT(IN) :: prsfl     (kbdim)              ! rain large scale
    REAL(wp),OPTIONAL,INTENT(IN) :: prsfc     (kbdim)              ! rain convective
    REAL(wp),OPTIONAL,INTENT(IN) :: pssfl     (kbdim)              ! snow large scale
    REAL(wp),OPTIONAL,INTENT(IN) :: pssfc     (kbdim)              ! snow convective
    REAL(wp),OPTIONAL,INTENT(IN) :: pemterall (kbdim)              ! net surface longwave flux [W/m2]
    REAL(wp),OPTIONAL,INTENT(IN) :: ptrsolall (kbdim)              ! net surface shortwave flux [W/m2]
    REAL(wp),OPTIONAL,INTENT(IN) :: presi_old (kbdim,klev+1)       ! half level pressure
    REAL(wp),OPTIONAL,INTENT(IN) :: pcosmu0   (kbdim)              ! cos of zenith angle
    REAL(wp),OPTIONAL,INTENT(IN) :: pch_tile  (kbdim,ksfc_type)
    !! JSBACH output
    REAL(wp),OPTIONAL,INTENT(INOUT) :: pcsat(kbdim)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: pcair(kbdim)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: zhsoil(kbdim)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: tte_corr(kbdim)  ! OUT
    REAL(wp),OPTIONAL,INTENT(INOUT) :: z0h_lnd(kbdim), z0m_lnd(kbdim)  ! OUT
    REAL(wp),OPTIONAL,INTENT(INOUT) :: albvisdir(kbdim)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: albnirdir(kbdim)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: albvisdif(kbdim)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: albnirdif(kbdim)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: surface_temperature_rad(kbdim) ! OUT
    REAL(wp),OPTIONAL,INTENT(INOUT) :: surface_temperature_eff(kbdim) ! OUT
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
    REAL(wp),OPTIONAL,INTENT(INOUT) :: albvisdir_wtr(kbdim)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: albvisdif_wtr(kbdim)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: albnirdir_wtr(kbdim)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: albnirdif_wtr(kbdim)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: pswflx_tile(kbdim,ksfc_type) ! OUT
    REAL(wp),OPTIONAL,INTENT(INOUT) :: plwflx_tile(kbdim,ksfc_type) ! OUT

! locals

    INTEGER  :: ilsm(kbdim)
    INTEGER  :: jsfc, jk, jkm1, im, k
    REAL(wp) :: se_sum(kbdim), qv_sum(kbdim), wgt_sum(kbdim), wgt(kbdim)
    REAL(wp) :: zca(kbdim,ksfc_type), zcs(kbdim,ksfc_type)
    REAL(wp) :: zfrc_oce(kbdim)

    REAL(wp) :: zen_h (kbdim,ksfc_type)
    REAL(wp) :: zfn_h (kbdim,ksfc_type)
    REAL(wp) :: zen_qv(kbdim,ksfc_type)
    REAL(wp) :: zfn_qv(kbdim,ksfc_type)

    REAL(wp) :: &
      & lwup(kbdim), evapotranspiration(kbdim),                        &
      & surface_temperature(kbdim), surface_temperature_last(kbdim),   &
      & sat_surface_specific_humidity(kbdim), dry_static_energy(kbdim)

    ! Sea ice
    REAL(wp) :: Tfw(kbdim), LWin(kbdim)
    REAL(wp) :: swflx_ice(kbdim,kice), nonsolar_ice(kbdim,kice), dnonsolardT(kbdim,kice)

#ifdef __NO_JSBACH__
   CHARACTER(len=*), PARAMETER :: method_name='mo_surface:update_surface'
#endif
#ifdef __ICON_OCEAN__
   CHARACTER(len=*), PARAMETER :: method_name='mo_surface:update_surface'
#endif

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

    ! Turbulent transport of moisture:
    ! - finish matrix set up;
    ! - perform bottom level elimination;
    ! - convert matrix entries to Richtmyer-Morton coefficients

    IF (phy_config%ljsbach) THEN
#ifndef __NO_JSBACH__

      CALL matrix_to_richtmyer_coeff( kproma, kbdim, klev, ksfc_type, idx_lnd, &! in
                                    & aa(:,:,:,imh:imqv), bb(:,:,ih:iqv),      &! in
                                    & aa_btm, bb_btm,                          &! inout
                                    & zen_h, zfn_h, zen_qv, zfn_qv,            &! out
                                    & pcair = pcair(:),                        &! in
                                    & pcsat = pcsat(:))                         ! in

      lwup(1:kproma) = zemiss_def * stbo * ptsfc_tile(1:kproma,idx_lnd)**4

      surface_temperature(:) = 0._wp
      surface_temperature_last(:) = 0._wp
      sat_surface_specific_humidity(:) = 0._wp
      dry_static_energy(:) = 0._wp
      evapotranspiration(:) = 0._wp
      surface_temperature_last(1:kproma) = ptsfc_tile(1:kproma,idx_lnd)

      z0m_lnd(:) = 0._wp
      z0h_lnd(:) = 0._wp
      tte_corr(:) = 0._wp

      WHERE (lsm(:) > 0.5_wp)
        ilsm(:) = 1
      ELSEWHERE
        ilsm(:) = 0
      ENDWHERE

      CALL jsbach_interface ( jg, nblock, 1, kproma, pdtime, psteplen,                  & ! in
        & t_air            = ptemp(1:kproma),                                           & ! in
        & q_air            = pq(1:kproma),                                              & ! in
        & rain             = prsfl(1:kproma) + prsfc(1:kproma),                         & ! in
        & snow             = pssfl(1:kproma) + pssfc(1:kproma),                         & ! in
        & wind_air         = SQRT(pu(1:kproma)**2 + pv(1:kproma)**2),                   & ! in
        & wind_10m         = SQRT(pu(1:kproma)**2 + pv(1:kproma)**2),                   & ! in, temporary
        & lwrad_srf_down   = pemterall(1:kproma) + lwup(1:kproma),                      & ! in
        & swrad_srf_net    = ptrsolall(1:kproma),                                       & ! in
        & press_srf        = presi_old(1:kproma,klev+1),                                & ! in
        & drag_srf         = pfac_sfc(1:kproma) * pcfh_tile(1:kproma,idx_lnd),          & ! in
        & t_acoef          = zen_h(1:kproma, idx_lnd),                                  & ! in
        & t_bcoef          = zfn_h(1:kproma, idx_lnd),                                  & ! in
        & q_acoef          = zen_qv(1:kproma, idx_lnd),                                 & ! in
        & q_bcoef          = zfn_qv(1:kproma, idx_lnd),                                 & ! in
        & pch              = MERGE(pch_tile(1:kproma,idx_lnd),1._wp,ilsm(1:kproma)>0),  & ! in
        & cos_zenith_angle = pcosmu0(1:kproma),                                         & ! in
        & t_srf            = surface_temperature(1:kproma),                             & ! out
        & t_rad_srf        = surface_temperature_rad(1:kproma),                         & ! out
        & qsat_srf         = sat_surface_specific_humidity(1:kproma),                   & ! out
        & s_srf            = dry_static_energy(1:kproma),                               & ! out
        & fact_q_air       = pcair(1:kproma),                                           & ! out
        & fact_qsat_srf    = pcsat(1:kproma),                                           & ! out
        & rel_hum_srf      = zhsoil(1:kproma),                                          & ! out
        & evapotrans       = evapotranspiration(1:kproma),                              & ! out
        & latent_hflx      = plhflx_tile(1:kproma, idx_lnd),                            & ! out
        & sensible_hflx    = pshflx_tile(1:kproma, idx_lnd),                            & ! out
        & rough_h_srf      = z0h_lnd(1:kproma),                                         & ! out
        & rough_m_srf      = z0m_lnd(1:kproma),                                         & ! out
        & tte_corr         = tte_corr(1:kproma),                                        & ! out
        & alb_vis_dir      = albvisdir(1:kproma),                                       & ! out
        & alb_nir_dir      = albnirdir(1:kproma),                                       & ! out
        & alb_vis_dif      = albvisdif(1:kproma),                                       & ! out
        & alb_nir_dif      = albnirdif(1:kproma)                                        & ! out
        )

      WHERE (lsm(1:kproma) > 0.5_wp)
        ptsfc_tile(1:kproma,idx_lnd) = surface_temperature(1:kproma)
        pcpt_tile (1:kproma,idx_lnd) = dry_static_energy(1:kproma)
        pqsat_tile(1:kproma,idx_lnd) = sat_surface_specific_humidity(1:kproma)
      END WHERE

      tte_corr(1:kproma) = tte_corr(1:kproma) / (presi_old(1:kproma,klev+1) - presi_old(1:kproma,klev))

      ! calculate effective temperature for use in radheat
      WHERE (lsm(1:kproma) > 0.5_wp)
        surface_temperature_eff(1:kproma) = (surface_temperature_last(1:kproma) ** 3 *  &
                                            (4._wp*surface_temperature_rad(1:kproma) -  &
                                            3._wp * surface_temperature_last(1:kproma)))**0.25
      ELSEWHERE
        surface_temperature_eff(1:kproma) = ptsfc_tile(1:kproma,idx_wtr)
        surface_temperature_rad(1:kproma) = ptsfc_tile(1:kproma,idx_wtr)
      ENDWHERE
#else
      CALL finish(method_name, "The JSBACH component is not activated")
#endif
    ELSE
      CALL matrix_to_richtmyer_coeff( kproma, kbdim, klev, ksfc_type, idx_lnd, &! in
                                    & aa(:,:,:,imh:imqv), bb(:,:,ih:iqv),      &! in
                                    & aa_btm, bb_btm,                          &! inout
                                    & zen_h, zfn_h, zen_qv, zfn_qv             )! out

    END IF

    ! Set the evapotranspiration coefficients, to be used later in
    ! blending and in diagnoising surface fluxes.
    !
    zca(1:kproma,:) = 1._wp
    zcs(1:kproma,:) = 1._wp

    IF (idx_lnd <= ksfc_type .AND. phy_config%ljsbach) THEN
      zca(1:kproma,idx_lnd) = pcair(1:kproma)
      zcs(1:kproma,idx_lnd) = pcsat(1:kproma)
    END IF

    ! CALL sea_ice_thermodynamics ?

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
           wgt(1:kproma) = pfrc(1:kproma,jsfc)*pcfh_tile(1:kproma,jsfc)*pfac_sfc(1:kproma)
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
                      & kproma, kbdim, klev, ksfc_type,       &! in
                      & idx_wtr, idx_ice, idx_lnd, ih, iqv,   &! in
                      & pfrc, pcfh_tile, pfac_sfc,            &! in
                      & pcpt_tile, ptsfc_tile, pqsat_tile,    &! in
                      & zca, zcs, bb(:,:,ih:iqv),             &! in
                      & plhflx_gbm, pshflx_gbm,               &! out
                      & pevap_gbm,                            &! out
                      & plhflx_tile, pshflx_tile,             &! inout
                      & dshflx_dT_tile,                       &! out
                      & pevap_tile,                           &! out
                      & evapotranspiration)                    ! in (optional)


    IF(phy_config%lice) THEN
#ifndef __NO_ICON_OCEAN__
! LL This is a temporary solution,
! we should restrcure ice thermodynamics in a more stand-alone way

! For explicit coupling to ice:
    IF ( idx_ice <= nsfc_type ) THEN
! Freezing point of sea-water
      Tfw = Tf
! ECHAM has no tiles for SW & LW and this is how it's solved there
! Net shortwave on all bands. 0.25*ptrsolall to be replaced by field%vissfc, etc.
! We need to divide with the box albedo because there are no tiles yet
! Net longwave - we don't have tiles yet
      LWin(1:kproma) = pemterall(1:kproma)
      DO jsfc=1,nsfc_type
        LWin(1:kproma) = LWin(1:kproma) + zemiss_def * stbo * &
          & pfrc(1:kproma,jsfc) * ptsfc_tile(1:kproma,jsfc)**4
      ENDDO
! First all ice classes
      DO k=1,kice
!        swflx_ice(1:kproma,k) = &
!          &   0.28_wp * ptrsolall(1:kproma) / (1._wp-albvisdir(1:kproma)) * (1._wp-albvisdir_ice(1:kproma,k)) + &
!          &   0.24_wp * ptrsolall(1:kproma) / (1._wp-albvisdif(1:kproma)) * (1._wp-albvisdif_ice(1:kproma,k)) + &
!          &   0.31_wp * ptrsolall(1:kproma) / (1._wp-albnirdir(1:kproma)) * (1._wp-albnirdir_ice(1:kproma,k)) + &
!          &   0.17_wp * ptrsolall(1:kproma) / (1._wp-albnirdif(1:kproma)) * (1._wp-albnirdif_ice(1:kproma,k))
        swflx_ice(1:kproma,k) = ptrsolall(1:kproma) * (                                   &
          &   0.28_wp / (1._wp-albvisdir(1:kproma)) * (1._wp-albvisdir_ice(1:kproma,k)) + &
          &   0.24_wp / (1._wp-albvisdif(1:kproma)) * (1._wp-albvisdif_ice(1:kproma,k)) + &
          &   0.31_wp / (1._wp-albnirdir(1:kproma)) * (1._wp-albnirdir_ice(1:kproma,k)) + &
          &   0.17_wp / (1._wp-albnirdif(1:kproma)) * (1._wp-albnirdif_ice(1:kproma,k)) )

        nonsolar_ice(1:kproma,k) = LWin(1:kproma) - zemiss_def * stbo * (Tsurf(1:kproma,k)+tmelt)**4 &
          &     + plhflx_tile(1:kproma,idx_ice) + pshflx_tile(1:kproma,idx_ice)

        dnonsolardT(1:kproma,k) = -4._wp * zemiss_def * stbo * (Tsurf(1:kproma,k)+tmelt)**3

      ENDDO
! Then open water
      IF ( PRESENT(plwflx_tile) ) THEN
        plwflx_tile(1:kproma,idx_wtr) = LWin(1:kproma) - zemiss_def * stbo * ptsfc_tile(1:kproma,idx_wtr)**4
        plwflx_tile(1:kproma,idx_ice) = LWin(1:kproma) - zemiss_def * stbo * (Tsurf(1:kproma,1)+tmelt)**4
      ENDIF

      IF ( PRESENT(pswflx_tile) ) THEN
        pswflx_tile(1:kproma,idx_wtr) = ptrsolall(1:kproma) * (                                 &
          &        0.28_wp / (1._wp-albvisdir(1:kproma)) * (1._wp-albvisdir_wtr(1:kproma)) +    &
          &        0.24_wp / (1._wp-albvisdif(1:kproma)) * (1._wp-albvisdif_wtr(1:kproma)) +    &
          &        0.31_wp / (1._wp-albnirdir(1:kproma)) * (1._wp-albnirdir_wtr(1:kproma)) +    &
          &        0.17_wp / (1._wp-albnirdif(1:kproma)) * (1._wp-albnirdif_wtr(1:kproma)) )
        pswflx_tile(1:kproma,idx_ice) = swflx_ice(1:kproma,1)
      ENDIF

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

! Albedo model for the ocean
      albvisdir_wtr(1:kproma) = albedoW
      albvisdif_wtr(1:kproma) = albedoW
      albnirdir_wtr(1:kproma) = albedoW
      albnirdif_wtr(1:kproma) = albedoW

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
          WHERE ( hi(:,k) > 0._wp )
            ! Snow only falls when it's below freezing
            WHERE ( Tsurf(:,k) < 0._wp )
              hs(:,k) = hs(:,k) + (pssfl + pssfc)*pdtime/rhos 
            ENDWHERE
            ! Snow melt
            hs(:,k) = hs(:,k) - MIN( Qtop(:,k)*pdtime/( alf*rhos ), hs(:,k) )
          ELSEWHERE
            hs(:,k) = 0._wp
          ENDWHERE
        ENDDO
      ENDIF
! Average the albedo.
      IF ( idx_lnd <= nsfc_type ) THEN

        WHERE ( pfrc(1:kproma,idx_lnd) < 1 )
          albvisdir(1:kproma) = pfrc(1:kproma,idx_wtr) * albvisdir_wtr(1:kproma) + &
            & SUM( conc(1:kproma,:) * albvisdir_ice(1:kproma,:), 2 )
          albvisdif(1:kproma) = pfrc(1:kproma,idx_wtr) * albvisdif_wtr(1:kproma) + &
            & SUM( conc(1:kproma,:) * albvisdif_ice(1:kproma,:), 2 )
          albnirdir(1:kproma) = pfrc(1:kproma,idx_wtr) * albnirdir_wtr(1:kproma) + &
            & SUM( conc(1:kproma,:) * albnirdir_ice(1:kproma,:), 2 )
          albnirdif(1:kproma) = pfrc(1:kproma,idx_wtr) * albnirdif_wtr(1:kproma) + &
            & SUM( conc(1:kproma,:) * albnirdif_ice(1:kproma,:), 2 )
        ENDWHERE

      ELSE

        albvisdir(1:kproma) = pfrc(1:kproma,idx_wtr) * albvisdir_wtr(1:kproma) + &
          & SUM( conc(1:kproma,:) * albvisdir_ice(1:kproma,:), 2 )
        albvisdif(1:kproma) = pfrc(1:kproma,idx_wtr) * albvisdif_wtr(1:kproma) + &
          & SUM( conc(1:kproma,:) * albvisdif_ice(1:kproma,:), 2 )
        albnirdir(1:kproma) = pfrc(1:kproma,idx_wtr) * albnirdir_wtr(1:kproma) + &
          & SUM( conc(1:kproma,:) * albnirdir_ice(1:kproma,:), 2 )
        albnirdif(1:kproma) = pfrc(1:kproma,idx_wtr) * albnirdif_wtr(1:kproma) + &
          & SUM( conc(1:kproma,:) * albnirdif_ice(1:kproma,:), 2 )
!        albvisdir(1:kproma) = pfrc(:,idx_wtr)*albvisdir_wtr + SUM(conc(:,:)*albvisdir_ice(:,:), 2)
!        albvisdif(1:kproma) = pfrc(:,idx_wtr)*albvisdif_wtr + SUM(conc(:,:)*albvisdif_ice(:,:), 2)
!        albnirdir(1:kproma) = pfrc(:,idx_wtr)*albnirdir_wtr + SUM(conc(:,:)*albnirdir_ice(:,:), 2)
!        albnirdif(1:kproma) = pfrc(:,idx_wtr)*albnirdif_wtr + SUM(conc(:,:)*albnirdif_ice(:,:), 2)

      ENDIF
! Set the tile temperature
      ptsfc_tile(1:kproma,idx_ice) = Tsurf(1:kproma,1) + tmelt
    ENDIF
#else
    ! not __ICON_OCEAN__
      CALL finish(method_name, "The ice process requires the ICON_OCEAN component")
#endif
    ENDIF ! lice

  END SUBROUTINE update_surface
  !-------------

END MODULE mo_surface
