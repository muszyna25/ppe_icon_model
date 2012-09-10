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
  USE mo_echam_phy_config,    ONLY: phy_config => echam_phy_config
  USE mo_vdiff_solver,      ONLY: ih, iqv, iu, iv, imh, imqv, imuv, &
                                & nmatrix, nvar_vdiff,              &
                                & matrix_to_richtmyer_coeff
  USE mo_jsbach_interface_icon,ONLY: jsbach_inter_1d
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: update_surface

CONTAINS
  !>
  !!
  !!
  SUBROUTINE update_surface( lsfc_heat_flux, lsfc_mom_flux,     &! in
                           & pdtime, psteplen,                  &! in
                           & kproma, kbdim,                     &! in
                           & klev, ksfc_type,                   &! in
                           & idx_wtr, idx_ice, idx_lnd,         &! in
                           & pfrc,                              &! in
                           & pcfh_tile, pcfm_tile,              &! in
                           & pfac_sfc, pocu, pocv,              &! in
                           & aa, aa_btm, bb, bb_btm,            &! inout
                           & pcpt_tile, pqsat_tile,             &! inout
                           & ptsfc_tile,                        &! inout
                           & pu_stress_gbm_ac, pv_stress_gbm_ac,&! inout
                           & plhflx_gbm_ac, pshflx_gbm_ac,      &! inout
                           & pevap_gbm_ac,  dshflx_dT_ac_tile,  &! inout
                           & pu_stress_tile,   pv_stress_tile,  &! out
                           & plhflx_tile, pshflx_tile,          &! out
                           & dshflx_dT_tile,                    &! out
                           & pevap_tile, pevap_gbm,             &! out
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
                           !! added for testing JSBACH (hydrology)
                           & pcsat,                             &! inout
                           & pcair,                             &! inout
                           & csat_transpiration,                &! inout
                           & moisture1,                         &! inout
                           & moisture2,                         &! inout
                           & moisture3,                         &! inout
                           & moisture4,                         &! inout
                           & moisture5,                         &! inout
                           & moisture_all,                      &! inout
                           & sat_surface_specific_humidity,     &! inout
                           & skin_reservoir,                    &! inout
                           & snow_fract,                        &! inout
                           & snow,                              &! inout
                           & snow_canopy,                       &! inout
                           & snow_melt,                         &! inout
                           & snow_acc,                          &! inout
                           & snow_melt_acc,                     &! inout
                           & glacier_runoff_acc,                &! inout
                           & runoff_acc,                        &! inout
                           & drainage_acc,                      &! inout
                           !! added for testing JSBACH (energy balance)
                           & surface_temperature,               &! inout
                           & surface_temperature_old,           &! inout
                           & c_soil_temperature1,               &! inout
                           & c_soil_temperature2,               &! inout
                           & c_soil_temperature3,               &! inout
                           & c_soil_temperature4,               &! inout
                           & c_soil_temperature5,               &! inout
                           & d_soil_temperature1,               &! inout
                           & d_soil_temperature2,               &! inout
                           & d_soil_temperature3,               &! inout
                           & d_soil_temperature4,               &! inout
                           & d_soil_temperature5,               &! inout
                           & soil_temperature1,                 &! inout
                           & soil_temperature2,                 &! inout
                           & soil_temperature3,                 &! inout
                           & soil_temperature4,                 &! inout
                           & soil_temperature5,                 &! inout
                           & heat_capacity,                     &! inout
                           & ground_heat_flux,                  &! inout
                           & swnet,                             &! inout
                           & time_steps_soil,                   &! inout
                           & evapotranspiration,                &! out
                           & surface_temperature_rad,           &! out
                           & surface_temperature_eff            &! out
                           )

    LOGICAL, INTENT(IN) :: lsfc_heat_flux, lsfc_mom_flux
    REAL(wp),INTENT(IN) :: pdtime, psteplen
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
    REAL(wp),INTENT(INOUT) :: pu_stress_gbm_ac (kbdim)
    REAL(wp),INTENT(INOUT) :: pv_stress_gbm_ac (kbdim)
    REAL(wp),INTENT(INOUT) :: plhflx_gbm_ac (kbdim)
    REAL(wp),INTENT(INOUT) :: pshflx_gbm_ac (kbdim)
    REAL(wp),INTENT(INOUT) :: pevap_gbm_ac (kbdim)
    REAL(wp),INTENT(INOUT) :: dshflx_dT_ac_tile(kbdim,ksfc_type)

    REAL(wp),INTENT(INOUT) :: pu_stress_tile (kbdim,ksfc_type) ! practically out
    REAL(wp),INTENT(INOUT) :: pv_stress_tile (kbdim,ksfc_type) ! practically out
    REAL(wp),INTENT(INOUT) ::    plhflx_tile (kbdim,ksfc_type) ! practically out
    REAL(wp),INTENT(INOUT) ::    pshflx_tile (kbdim,ksfc_type) ! practically out
    REAL(wp),INTENT(OUT)   :: dshflx_dT_tile (kbdim,ksfc_type)
    REAL(wp),INTENT(OUT)   :: pevap_tile (kbdim,ksfc_type)
    REAL(wp),INTENT(OUT)   :: pevap_gbm  (kbdim)

    !! added as JSBACH input
    INTEGER, OPTIONAL,INTENT(IN) :: nblock
    INTEGER, OPTIONAL,INTENT(IN) :: lsm(kbdim)
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
    REAL(wp),OPTIONAL,INTENT(IN) :: presi_old (kbdim)              ! surface pressure
    REAL(wp),OPTIONAL,INTENT(IN) :: pcosmu0   (kbdim)              ! cos of zenith angle
    REAL(wp),OPTIONAL,INTENT(IN) :: pch_tile  (kbdim,ksfc_type)
    !! added for testing JSBACH (hydrology)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: pcsat(kbdim)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: pcair(kbdim)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: csat_transpiration(kbdim)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: moisture1(kbdim)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: moisture2(kbdim)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: moisture3(kbdim)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: moisture4(kbdim)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: moisture5(kbdim)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: moisture_all(kbdim)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: sat_surface_specific_humidity(kbdim)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: skin_reservoir(kbdim)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: snow_fract(kbdim)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: snow(kbdim)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: snow_canopy(kbdim)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: snow_melt(kbdim)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: snow_acc(kbdim)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: snow_melt_acc(kbdim)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: glacier_runoff_acc(kbdim)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: runoff_acc(kbdim)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: drainage_acc(kbdim)
    !! added for testing JSBACH (energy balance)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: surface_temperature(kbdim)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: surface_temperature_old(kbdim)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: c_soil_temperature1(kbdim)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: c_soil_temperature2(kbdim)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: c_soil_temperature3(kbdim)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: c_soil_temperature4(kbdim)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: c_soil_temperature5(kbdim)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: d_soil_temperature1(kbdim)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: d_soil_temperature2(kbdim)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: d_soil_temperature3(kbdim)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: d_soil_temperature4(kbdim)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: d_soil_temperature5(kbdim)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: soil_temperature1(kbdim)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: soil_temperature2(kbdim)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: soil_temperature3(kbdim)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: soil_temperature4(kbdim)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: soil_temperature5(kbdim)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: heat_capacity(kbdim)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: ground_heat_flux(kbdim)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: swnet(kbdim)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: time_steps_soil(kbdim)
    REAL(wp),OPTIONAL,INTENT(OUT)   :: evapotranspiration(kbdim)
    REAL(wp),OPTIONAL,INTENT(OUT)   :: surface_temperature_rad(kbdim)
    REAL(wp),OPTIONAL,INTENT(OUT)   :: surface_temperature_eff(kbdim)

! locals

    LOGICAL  :: lfland(kbdim)
    REAL(wp)  :: surface_temperature_last(kbdim)
    INTEGER  :: jsfc, jk, jkm1, im
    REAL(wp) :: se_sum(kbdim), qv_sum(kbdim), wgt_sum(kbdim), wgt(kbdim)
    REAL(wp) :: zca(kbdim,ksfc_type), zcs(kbdim,ksfc_type)
    REAL(wp) :: zfrc_oce(kbdim)

    REAL(wp) :: zen_h (kbdim,ksfc_type)
    REAL(wp) :: zfn_h (kbdim,ksfc_type)
    REAL(wp) :: zen_qv(kbdim,ksfc_type)
    REAL(wp) :: zfn_qv(kbdim,ksfc_type)

    REAL(wp) :: lwup(kbdim)

    !===================================================================
    ! BEFORE CALLING land/ocean/ice model
    !===================================================================
    ! Compute wind stress at the old time step.
    ! At this point bb(:,klev,iu) = u_klev(t)/tpfac1 (= udif in echam)
    !               bb(:,klev,iv) = v_klev(t)/tpfac1 (= vdif in echam)

    CALL wind_stress( lsfc_mom_flux, pdtime, psteplen,     &! in
                    & kproma, kbdim, ksfc_type,            &! in
                    & pfrc, pcfm_tile, pfac_sfc,           &! in
                    & bb(:,klev,iu), bb(:,klev,iv),        &! in
                    & pu_stress_gbm_ac, pv_stress_gbm_ac,  &! inout
                    & pu_stress_tile,   pv_stress_tile     )! out

    ! Turbulent transport of moisture:
    ! - finish matrix set up;
    ! - perform bottom level elimination;
    ! - convert matrix entries to Richtmyer-Morton coefficients



    IF (phy_config%ljsbach) THEN
    CALL matrix_to_richtmyer_coeff( kproma, kbdim, klev, ksfc_type, idx_lnd, &! in
                                  & aa(:,:,:,imh:imqv), bb(:,:,ih:iqv),      &! in
                                  & aa_btm, bb_btm,                          &! inout
                                  & zen_h, zfn_h, zen_qv, zfn_qv,            &! out
                                  & pcair = pcair(:),                        &! in
                                  & pcsat = pcsat(:))                         ! in

    lfland(1:kproma) = .true. !! TR for testing
!!$ TR long wave upward radiation for water (aqua planet setup)
    lwup(1:kproma) = 0.996_wp * 5.67e-8_wp * ptsfc_tile(1:kproma,idx_lnd)**4
    surface_temperature_last(1:kproma) = surface_temperature(1:kproma)

    CALL jsbach_inter_1d (kdim = kproma,                                   &
                          kblock = nblock,                                 &
                          kland = COUNT(lfland(1:kproma)),                 &
                          mask_land = lfland(1:kproma),                    &
                          delta_time = pdtime,                             &
                          time_step_len = psteplen,                        &
                          time_steps_soil = time_steps_soil(1:kproma),     &
                          wind = SQRT(pu(1:kproma)**2 + pv(1:kproma)**2),  &
!!$ TR                          wind10 =     ! comment out to simplify for testing, 10 m wind speed will be implemented later
!!$ TR                                       ! set wind10 in update soil to wind for testing
                          temp_air = ptemp(1:kproma),                      &
                          qair     = pq(1:kproma),                         &
                          precip_rain = prsfl(1:kproma) + prsfc(1:kproma), &
                          precip_snow = pssfl(1:kproma) + pssfc(1:kproma), &
                          lwdown   = pemterall(1:kproma) + lwup(1:kproma), &
                          swdown   = ptrsolall(1:kproma),                  & !! ATTENTION !! replaces sw_vis_net + sw_nir_net
                          pressure = presi_old(1:kproma),                  &
!!$ TR                          czenith  = pcosmu0(1:kproma),                    & ! comment out to simplify for testing
!!$ TR                          CO2_concentration = 0.00055_wp,                  & ! comment out to simplify for testing, mass mixing ratio
                          cdrag = pfac_sfc(1:kproma) * pcfh_tile(1:kproma,idx_lnd), &
                          etAcoef = zen_h(1:kproma,idx_lnd),               &
                          etBcoef = zfn_h(1:kproma,idx_lnd),               &
                          eqAcoef = zen_qv(1:kproma,idx_lnd),              &
                          eqBcoef = zfn_qv(1:kproma,idx_lnd),              &
                          p_echam_zchl = pch_tile(1:kproma,idx_lnd),       & ! intent in
                           !! added for testing JSBACH (hydrology)
                          cair = pcair(1:kproma),                          & ! intent out
                          csat = pcsat(1:kproma),                          & ! intent out
!!$ TR                          zhsoil = ,                                      & ! intent out
                          csat_transpiration= csat_transpiration(1:kproma),& ! intent out
                          moisture1 = moisture1(1:kproma),                 & ! intent out
                          moisture2 = moisture2(1:kproma),                 & ! intent out
                          moisture3 = moisture3(1:kproma),                 & ! intent out
                          moisture4 = moisture4(1:kproma),                 & ! intent out
                          moisture5 = moisture5(1:kproma),                 & ! intent out
                          moisture_all = moisture_all(1:kproma),           & ! intent out
                          sat_surface_specific_humidity = sat_surface_specific_humidity(1:kproma),& ! intent out
                          skin_reservoir = skin_reservoir(1:kproma),       & ! intent out
                          snow_fract = snow_fract(1:kproma),               & ! intent out
                          snow = snow(1:kproma),                           & ! intent out
                          snow_canopy = snow_canopy(1:kproma),             & ! intent out
                          snow_melt = snow_melt(1:kproma),                 & ! intent out
                          snow_acc = snow_acc(1:kproma),                   & ! intent out
                          snow_melt_acc = snow_melt_acc(1:kproma),         & ! intent out
                          glacier_runoff_acc = glacier_runoff_acc(1:kproma), & ! intent out
                          runoff_acc = runoff_acc(1:kproma),               & ! intent out
                          drainage_acc = drainage_acc(1:kproma),           & ! intent out
                           !! added for testing JSBACH (energy balance)
                          surface_temperature = surface_temperature(1:kproma), &
                          surface_temperature_old = surface_temperature_old(1:kproma), &
                          surface_temperature_rad = surface_temperature_rad(1:kproma), &
                          evapotranspiration = evapotranspiration(1:kproma), &
                          c_soil_temperature1 = c_soil_temperature1(1:kproma), &
                          c_soil_temperature2 = c_soil_temperature2(1:kproma), &
                          c_soil_temperature3 = c_soil_temperature3(1:kproma), &
                          c_soil_temperature4 = c_soil_temperature4(1:kproma), &
                          c_soil_temperature5 = c_soil_temperature5(1:kproma), &
                          d_soil_temperature1 = d_soil_temperature1(1:kproma), &
                          d_soil_temperature2 = d_soil_temperature2(1:kproma), &
                          d_soil_temperature3 = d_soil_temperature3(1:kproma), &
                          d_soil_temperature4 = d_soil_temperature4(1:kproma), &
                          d_soil_temperature5 = d_soil_temperature5(1:kproma), &
                          soil_temperature1 = soil_temperature1(1:kproma), &
                          soil_temperature2 = soil_temperature2(1:kproma), &
                          soil_temperature3 = soil_temperature3(1:kproma), &
                          soil_temperature4 = soil_temperature4(1:kproma), &
                          soil_temperature5 = soil_temperature5(1:kproma), &
                          heat_capacity = heat_capacity(1:kproma), &
                          ground_heat_flux = ground_heat_flux(1:kproma), &
                          swnet = swnet(1:kproma) &
                          )


    ! calculate effective temperature for use in radheat
    WHERE (lsm(1:kproma) > 0) 
      surface_temperature_eff(1:kproma) = (surface_temperature_last(1:kproma) ** 3 *  &
                                          (4._wp*surface_temperature_rad(1:kproma) -  &
                                          3._wp * surface_temperature_last(1:kproma)))**0.25
    ELSEWHERE
      surface_temperature_eff(1:kproma) = ptsfc_tile(1:kproma,idx_wtr)
      surface_temperature_rad(1:kproma) = ptsfc_tile(1:kproma,idx_wtr)
    ENDWHERE

    ELSE ! not ljsbach

    CALL matrix_to_richtmyer_coeff( kproma, kbdim, klev, ksfc_type, idx_lnd, &! in
                                  & aa(:,:,:,imh:imqv), bb(:,:,ih:iqv),      &! in
                                  & aa_btm, bb_btm,                          &! inout
                                  & zen_h, zfn_h, zen_qv, zfn_qv             )! out

    END IF ! ljsbach

    ! Set the evapotranspiration coefficients, to be used later in
    ! blending and in diagnoising surface fluxes.

    zca(1:kproma,:) = 1._wp
    zcs(1:kproma,:) = 1._wp

    IF (idx_lnd<=ksfc_type) THEN
      IF (.NOT. phy_config%ljsbach) CALL finish('mo_surface','land surface JSBACH not chosen')
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
           wgt(1:kproma) =  pfrc(1:kproma,jsfc)*pcfh_tile(1:kproma,jsfc)*pfac_sfc(1:kproma)
       wgt_sum(1:kproma) = wgt_sum(1:kproma) + wgt(1:kproma)
        se_sum(1:kproma) =  se_sum(1:kproma) + bb_btm(1:kproma,jsfc,ih)*wgt(1:kproma)
        qv_sum(1:kproma) =  qv_sum(1:kproma) + bb_btm(1:kproma,jsfc,iqv) &
                         &                       *wgt(1:kproma)*zca(1:kproma,jsfc)
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

!!$ TR couple surface temperature of water
!!$ TR    WRITE (*,*) 'lsm',lsm(1:kproma)
    IF (phy_config%ljsbach) THEN
      ptsfc_tile(1:kproma,idx_lnd) = surface_temperature(1:kproma)
    END IF ! ljsbach
   !-------------------------------------------------------------------
   ! Various diagnostics
   !-------------------------------------------------------------------

   CALL surface_fluxes( lsfc_heat_flux, pdtime, psteplen,     &! in
                      & kproma, kbdim, klev, ksfc_type,       &! in
                      & idx_wtr, idx_ice, idx_lnd, ih, iqv,   &! in
                      & pfrc, pcfh_tile, pfac_sfc,            &! in
                      & pcpt_tile, ptsfc_tile, pqsat_tile,    &! in
                      & zca, zcs, bb(:,:,ih:iqv),             &! in
                      & plhflx_gbm_ac, pshflx_gbm_ac,         &! inout
                      & pevap_gbm_ac,  dshflx_dT_ac_tile,     &! inout
                      & plhflx_tile, pshflx_tile,             &! inout (practically out)
                      & dshflx_dT_tile,                       &! out
                      & pevap_tile, pevap_gbm,                &! out
                      & evapotranspiration)                    ! in (optional)


  END SUBROUTINE update_surface
  !-------------

END MODULE mo_surface
