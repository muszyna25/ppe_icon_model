!>
!! @brief Subroutine interface_echam_vdf calls the vertical diffusion and the surface schemes.
!!
!! @author Hui Wan, MPI-M
!! @author Marco Giorgetta, MPI-M
!!
!! @par Revision History
!!  Original version from ECHAM6 (revision 2028)
!!  Modified for ICOHAM by Hui Wan and Marco Giorgetta (2010)
!!  Modified for ICONAM by Marco Giorgetta (2014)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!

MODULE mo_interface_echam_vdf

  !$ser verbatim USE mo_ser_echam_vdiff_down, ONLY: serialize_vdiff_down_input => serialize_input, &
  !$ser&                                            serialize_vdiff_down_output => serialize_output
  !$ser verbatim USE mo_ser_echam_vdiff_up, ONLY: serialize_vdiff_up_input => serialize_input, &
  !$ser&                                          serialize_vdiff_up_output => serialize_output
  !$ser verbatim USE mo_ser_echam_update_surface, ONLY: serialize_update_surface_input => serialize_input, &
  !$ser&                                                serialize_update_surface_output => serialize_output
  !$ser verbatim USE mo_ser_echam_nsurf_diag, ONLY: serialize_nsurf_diag_input => serialize_input, &
  !$ser&                                            serialize_nsurf_diag_output => serialize_output

  USE mo_kind                ,ONLY: wp
  USE mtime                  ,ONLY: datetime

  USE mo_exception           ,ONLY: finish, warning

  USE mo_echam_phy_config    ,ONLY: echam_phy_config
  USE mo_echam_phy_memory    ,ONLY: t_echam_phy_field, prm_field, &
    &                               t_echam_phy_tend,  prm_tend

  USE mo_timer               ,ONLY: ltimer, timer_start, timer_stop, timer_vdf, &
    &                               timer_vdf_dn, timer_vdf_sf, timer_vdf_up

  USE mo_echam_rad_config    ,ONLY: echam_rad_config
  USE mo_ccycle_config       ,ONLY: ccycle_config
  USE mo_physical_constants  ,ONLY: amco2, amd
  USE mo_bc_greenhouse_gases ,ONLY: ghg_co2mmr

  USE mo_run_config          ,ONLY: ntracer, iqv, iqc, iqi, iqt, ico2
  USE mo_vdiff_downward_sweep,ONLY: vdiff_down
  USE mo_vdiff_upward_sweep  ,ONLY: vdiff_up
  USE mo_vdiff_solver        ,ONLY: nvar_vdiff, nmatrix, imh, imqv, ih_vdiff=>ih, iqv_vdiff=>iqv
  
  USE mo_echam_sfc_indices   ,ONLY: nsfc_type, iwtr, iice, ilnd
  USE mo_surface             ,ONLY: update_surface
  USE mo_surface_diag        ,ONLY: nsurf_diag

  IMPLICIT NONE
  PRIVATE
  PUBLIC  :: interface_echam_vdf

CONTAINS

  SUBROUTINE interface_echam_vdf(jg, jb, jcs, jce     ,&
       &                         nproma,nlev          ,& 
       &                         is_in_sd_ed_interval ,&
       &                         is_active            ,&
       &                         datetime_old         ,&
       &                         pdtime               )

    ! Arguments
    !
    INTEGER                 ,INTENT(in) :: jg,jb,jcs,jce
    INTEGER                 ,INTENT(in) :: nproma,nlev
    LOGICAL                 ,INTENT(in) :: is_in_sd_ed_interval
    LOGICAL                 ,INTENT(in) :: is_active
    TYPE(datetime)          ,POINTER    :: datetime_old
    REAL(wp)                ,INTENT(in) :: pdtime

    ! Pointers
    !
    LOGICAL                 ,POINTER    :: lparamcpl
    INTEGER                 ,POINTER    :: fc_vdf
    REAL(wp)                ,POINTER    :: vmr_co2
    TYPE(t_echam_phy_field) ,POINTER    :: field
    TYPE(t_echam_phy_tend)  ,POINTER    :: tend

    ! Local variables
    !
    INTEGER  :: nlevm1, nlevp1
    INTEGER  :: ntrac
    !
    REAL(wp) :: zxt_emis(nproma,ntracer-iqt+1)   !< tracer tendency due to surface emission

    REAL(wp) :: zcpt_sfc_tile(nproma,nsfc_type)  !< dry static energy at surface
   
    ! Coefficient matrices and right-hand-side vectors for the turbulence solver
    ! _btm refers to the lowest model level (i.e., full level "klev", not the surface)
    REAL(wp) :: zri_tile(nproma,nsfc_type)          !< Richardson number
    REAL(wp) :: zaa    (nproma,nlev,3,nmatrix)      !< coeff. matrices, all variables
    REAL(wp) :: zaa_btm(nproma,3,nsfc_type,imh:imqv)!< last row of coeff. matrix of heat and moisture
    REAL(wp) :: zbb    (nproma,nlev,nvar_vdiff)             !< r.h.s., all variables
    REAL(wp) :: zbb_btm(nproma,nsfc_type,ih_vdiff:iqv_vdiff)!< last row of r.h.s. of heat and moisture

    ! Temporary arrays used by VDIFF, JSBACH
    REAL(wp) :: zfactor_sfc(nproma)
    REAL(wp) :: zco2       (nproma)          !< co2 value passed on to jsbach

    REAL(wp) :: zqx      (nproma,nlev)       !< total cloud condensate
    REAL(wp) :: zcptgz   (nproma,nlev)       !< dry static energy
    REAL(wp) :: zthvvar  (nproma,nlev)       !< intermediate value of thvvar
    REAL(wp) :: dummy    (nproma,nlev)       !< to replace thvvar
    REAL(wp) :: dummyx   (nproma,nlev)       !< to replace xvar
    REAL(wp) :: ztottevn (nproma,nlev)       !< intermediate value of TTE
    REAL(wp) :: zch_tile (nproma,nsfc_type)  !< for "nsurf_diag"
!!$    REAL(wp) :: zchn_tile(nproma,nsfc_type)  !< for "nsurf_diag"
!!$    REAL(wp) :: zcdn_tile(nproma,nsfc_type)  !< for "nsurf_diag"
!!$    REAL(wp) :: zcfnc_tile(nproma,nsfc_type) !< for "nsurf_diag"
    REAL(wp) :: zbn_tile (nproma,nsfc_type)  !< for "nsurf_diag"
    REAL(wp) :: zbhn_tile(nproma,nsfc_type)  !< for "nsurf_diag"
    REAL(wp) :: zbm_tile (nproma,nsfc_type)  !< for "nsurf_diag"
    REAL(wp) :: zbh_tile (nproma,nsfc_type)  !< for "nsurf_diag"

    REAL(wp) :: zq_snocpymlt(nproma)         !< heating by melting of snow on the canopy [W/m2]
    !                                        !  which warms the lowermost atmospheric layer (JSBACH)
    
    REAL(wp) :: zq_rlw_impl (nproma)         !< additional heating by LW rad. due to implicit coupling
    !                                        !  in surface energy balance [W/m2]

    REAL(wp) :: mmr_co2

    IF (ltimer) CALL timer_start(timer_vdf)

    ! associate pointers
    lparamcpl => echam_phy_config(jg)%lparamcpl
    fc_vdf    => echam_phy_config(jg)%fc_vdf
    vmr_co2   => echam_rad_config(jg)%vmr_co2
    field     => prm_field(jg)
    tend      => prm_tend (jg)

    nlevm1 = nlev-1
    nlevp1 = nlev+1
    ntrac  = ntracer-iqt+1  ! number of tracers excluding water vapour and hydrometeors

!!$    ! Emission of aerosols or other tracers (not implemented yet)
!!$    IF (ntrac>0) THEN
!!$       CALL tracer_emission()
!!$    ENDIF
    !
    ! default setting for all tracers
    zxt_emis(jcs:jce,:) = 0._wp
    !

    ! set emissions of co2, if any (hardcoded, co2-tracer 5 eq. zxt_emis(:,2))

    IF (ntrac>=2) THEN
      IF (ccycle_config(jg)%iccy_co2conc .EQ. 1 .AND. ccycle_config(jg)%iccy_co2flux .EQ. 2) THEN
        zxt_emis(jcs:jce,2) = field%fco2nat(jcs:jce,jb)
      ELSE
        CALL warning('echam_vdf','co2 emissions not possible in this setup. Please check your settings!')
      END IF
    END IF
!!$    ! Dry deposition of aerosols or other tracers (not implemented yet)
!!$    CALL dry_deposition()

    !
    ! setting of co2 for jsbach
    !
    mmr_co2 = vmr_co2   * amco2/amd

    IF(ccycle_config(jg)%iccy_co2conc .EQ. 0) THEN   !  default
      zco2(:) = SPREAD(mmr_co2, DIM=1, NCOPIES=nproma)
    ELSE IF (ccycle_config(jg)%iccy_co2conc .EQ. 2 .OR. ccycle_config(jg)%iccy_co2conc .EQ. 4) THEN
      zco2(:) = SPREAD(ghg_co2mmr, DIM=1, NCOPIES=nproma)
    ELSE IF (ccycle_config(jg)%iccy_co2conc .EQ. 1) THEN
      zco2(:) = field% qtrc(:,nlev,jb,ico2)
    ELSE
      CALL finish('echam_vdf: setting of co2 not recommended. Please check your setting!')
    END IF

    IF ( is_in_sd_ed_interval ) THEN
       !
       IF ( is_active ) THEN
          ! Set dummy values to zero to prevent invalid floating point operations:
          dummy (:,:)=0._wp
          dummyx(:,:)=0._wp 
          !
          zqx(jcs:jce,:) =  prm_field(jg)%qtrc(jcs:jce,:,jb,iqc) &
               &           +prm_field(jg)%qtrc(jcs:jce,:,jb,iqi)
          !
          IF (ltimer) CALL timer_start(timer_vdf_dn)
          !
          ! Turbulent mixing, part I:
          ! - computation of exchange coefficients in the atmosphere and at the surface;
          ! - build up the tridiagonal linear algebraic system;
          ! - downward sweep (Gaussian elimination from top till level nlev-1)
          !
          !----------------------------------------------------------------------------------------
          ! Serialbox2 input fields serialization
          !$ser verbatim call serialize_vdiff_down_input(jb, jg, jce, nproma,&
          !$ser verbatim   nlev, nlevm1, nlevp1, ntrac, nsfc_type, iwtr, iice,&
          !$ser verbatim   ilnd, pdtime, field, zqx, zxt_emis, dummy, dummyx,&
          !$ser verbatim   zri_tile,zaa,zaa_btm,zbb,zbb_btm,zfactor_sfc,&
          !$ser verbatim   zcpt_sfc_tile, zcptgz, zthvvar, ztottevn)
          !
          CALL vdiff_down(jg,                              &! in
               &          jce, nproma, nlev, nlevm1,nlevp1,&! in
               &          ntrac, nsfc_type,                &! in
               &          iwtr, iice, ilnd,                &! in, indices of different surface types
               &          pdtime,                          &! in, time step
               &          field%coriol(:,jb),              &! in, Coriolis parameter
               &          field%   zf(:,:,jb),             &! in, geopot. height above sea level, full level
               &          field%   zh(:,:,jb),             &! in, geopot. height above sea level, half level
               &          field%frac_tile(:,jb,:),         &! in, area fraction of each sfc type
               &          field% ts_tile(:,jb,:),          &! in, surface temperature
               &          field% ocu (:,jb),               &! in, ocean sfc velocity, u-component
               &          field% ocv (:,jb),               &! in, ocean sfc velocity, v-component
               &          field% presi_old(:,nlevp1,jb),   &! in, sfc pressure
               &          field%   ua(:,:,jb),             &! in, um1
               &          field%   va(:,:,jb),             &! in, vm1
               &          field%   ta(:,:,jb),             &! in, tm1
               &          field% qtrc(:,:,jb,iqv),         &! in, qm1
               &          field% qtrc(:,:,jb,iqc),         &! in, xlm1
               &          field% qtrc(:,:,jb,iqi),         &! in, xim1
               &                  zqx(:,:),                &! in, xlm1 + xim1
               &          field% qtrc(:,:,jb,iqt:),        &! in, xtm1
               &          field% mair(:,:,jb),             &! in, moist     air mass
               &          field% mref(:,:,jb),             &! in, reference air mass
               &          field% presi_old(:,:,jb),        &! in, aphm1
               &          field% presm_old(:,:,jb),        &! in, apm1
               &          field%   tv(:,:,jb),             &! in, virtual temperaturea
               &          field% aclc(:,:,jb),             &! in, cloud fraction
               &          zxt_emis,                        &! in, zxtems
               &          dummy(:,:),                      &! in, variance of theta_v at step t-dt
               &          dummyx(:,:),                     &! in
               &          field% z0m_tile(:,jb,:),         &! in
               &          field%  tottem1(:,:,jb),         &! in, TTE at step t-dt
               &          field%  ustar(:,  jb),           &! inout
               &          field%  wstar(:,  jb),           &! out, convective velocity scale
               &          field%  wstar_tile(:,jb,:),      &! inout, convective velocity scale (each sfc type)
               &          field% qs_sfc_tile(:,jb,:),      &! out, sfc specific humidity at saturation
               &          field%  hdtcbl(:,jb),            &! out, for output
               &          field%      ri (:,:,jb),         &! out, for output
               &          zri_tile (:,:),                  &! out, for nsurf_diag
               &          field%  mixlen (:,:,jb),         &! out, for output
               &          field% cfm     (:,:,jb),         &! out, for output
               &          field% cfm_tile(:,jb,:),         &! out, for output and "vdiff_up"
               &          field% cfh     (:,:,jb),         &! out, for output
               &          field% cfh_tile(:,jb,:),         &! out, for output and "vdiff_up"
               &          field% cfv     (:,:,jb),         &! out, for output
               &          field% cftotte (:,:,jb),         &! out, for output
               &          field% cfthv   (:,:,jb),         &! out, for output
               &          zaa, zaa_btm, zbb, zbb_btm,      &! out, for "vdiff_up"
               &          zfactor_sfc(:),                  &! out, for "vdiff_up"
               &          zcpt_sfc_tile(:,:),              &! out, for "vdiff_up"
               &          zcptgz(:,:),                     &! out, for "vdiff_up"
               &          zthvvar(:,:),                    &! out, for "vdiff_up"
               &          field%   thvsig(:,  jb),         &! out, for "cucall"
               &          ztottevn (:,:),                  &! out, for "vdiff_up"
               &          zch_tile(:,:),                   &! out, for "nsurf_diag"
!!$               &          zchn_tile(:,:),                  &! out, for "nsurf_diag"
!!$               &          zcdn_tile(:,:),                  &! out, for "nsurf_diag"
!!$               &          zcfnc_tile(:,:),                 &! out, for "nsurf_diag"
               &          zbn_tile(:,:),                   &! out, for "nsurf_diag"
               &          zbhn_tile(:,:),                  &! out, for "nsurf_diag"
               &          zbm_tile(:,:),                   &! out, for "nsurf_diag"
               &          zbh_tile(:,:),                   &! out, for "nsurf_diag"
               &          pcsat = field% csat(:,jb),       &! in, optional, area fraction with wet land surface
               &          pcair = field% cair(:,jb),       &! in, optional, area fraction with wet land surface (air)
               &          paz0lh = field% z0h_lnd(:,jb))    ! in, optional, roughness length for heat over land
          !
          !----------------------------------------------------------------------------------------
          ! Serialbox2 output fields serialization
          !$ser verbatim call serialize_vdiff_down_output(jb, jg, jce, nproma,&
          !$ser verbatim   nlev, nlevm1, nlevp1, ntrac, nsfc_type, iwtr, iice,&
          !$ser verbatim   ilnd, pdtime, field, zri_tile, zaa, zaa_btm, zbb,&
          !$ser verbatim   zbb_btm, zfactor_sfc, zcpt_sfc_tile, zcptgz,&
          !$ser verbatim   zthvvar, ztottevn, zch_tile, zbn_tile, zbhn_tile,&
          !$ser verbatim   zbm_tile, zbh_tile)
          !
          IF (ltimer) CALL timer_stop(timer_vdf_dn)
          !
          !
          ! Surface processes that provide time-dependent lower boundary
          ! condition for wind, temperature, tracer concentraion, etc.
          !
          field% lhflx_tile(jcs:jce,jb,:) = 0._wp
          field% shflx_tile(jcs:jce,jb,:) = 0._wp
          field% evap_tile (jcs:jce,jb,:) = 0._wp
          !
          IF (ltimer) CALL timer_start(timer_vdf_sf)
          !
          !----------------------------------------------------------------------------------------
          ! Serialbox2 input fields serialization
          !$ser verbatim call serialize_update_surface_input(jb, jg, jce, nproma,&
          !$ser verbatim   nlev, nsfc_type, iwtr, iice, ilnd,&
          !$ser verbatim   pdtime, field, zfactor_sfc, zaa, zaa_btm, zbb, zbb_btm,&
          !$ser verbatim   zcpt_sfc_tile, zco2, zch_tile)
          !
          CALL update_surface(jg, jce, nproma, field%kice,                    &! in
               &              nlev, nsfc_type,                                &! in
               &              iwtr, iice, ilnd,                               &! in, indices of surface types
               &              pdtime,                                         &! in, time step
               &              field%frac_tile(:,jb,:),                        &! in, area fraction
               &              field% cfh_tile(:,jb,:),                        &! in, from "vdiff_down"
               &              field% cfm_tile(:,jb,:),                        &! in, from "vdiff_down"
               &              zfactor_sfc(:),                                 &! in, from "vdiff_down"
               &              field% ocu (:,jb),                              &! in, ocean sfc velocity, u-component
               &              field% ocv (:,jb),                              &! in, ocean sfc velocity, v-component
               &              zaa, zaa_btm, zbb, zbb_btm,                     &! inout
               &              zcpt_sfc_tile(:,:),                             &! inout, from "vdiff_down", for "vdiff_up"
               &              field%qs_sfc_tile(:,jb,:),                      &! inout, from "vdiff_down", for "vdiff_up"
               &              field% ts_tile(:,jb,:),                         &! inout
               &              field%u_stress    (:,  jb),                     &! out
               &              field%v_stress    (:,  jb),                     &! out
               &              field% lhflx      (:,  jb),                     &! out
               &              field% shflx      (:,  jb),                     &! out
               &              field%  evap      (:,  jb),                     &! out, for "cucall"
               &              field%u_stress_tile  (:,jb,:),                  &! out
               &              field%v_stress_tile  (:,jb,:),                  &! out
               &              field% lhflx_tile    (:,jb,:),                  &! out
               &              field% shflx_tile    (:,jb,:),                  &! out
               &              field%  evap_tile    (:,jb,:),                  &! out
               &              field%  fco2nat      (:,  jb),                  &! inout
               &              nblock = jb,                                    &! in
               &              lsm = field%lsmask(:,jb),                       &!< in, land-sea mask
               &              alake = field%alake(:,jb),                      &! in, lake fraction
               &              pu    = field% ua(:,nlev,jb),                   &! in, um1
               &              pv    = field% va(:,nlev,jb),                   &! in, vm1
               &              ptemp = field% ta(:,nlev,jb),                   &! in, tm1
               &              pq = field% qtrc(:,nlev,jb,iqv),                &! in, qm1
               &              pco2 = zco2,                                    &! in, co2, lowest level or fixed value
               &              prsfl = field% rsfl(:,jb),                      &! in, rain surface large scale (from cloud)
               &              prsfc = field% rsfc(:,jb),                      &! in, rain surface concective (from cucall)
               &              pssfl = field% ssfl(:,jb),                      &! in, snow surface large scale (from cloud)
               &              pssfc = field% ssfc(:,jb),                      &! in, snow surface concective (from cucall)
               &              rlds        = field% rlds (:,jb),               &! in,  downward surface  longwave flux [W/m2]
               &              rlus        = field% rlus (:,jb),               &! inout, upward surface  longwave flux [W/m2]
               &              rsds        = field% rsds (:,jb),               &! in,  downward surface shortwave flux [W/m2]
               &              rsus        = field% rsus (:,jb),               &! in,  upward surface shortwave flux [W/m2]
               !
               &              rvds_dir   = field%rvds_dir   (:,jb),           &! in, all-sky downward direct visible radiation at surface
               &              rpds_dir   = field%rpds_dir   (:,jb),           &! in, all-sky downward direct PAR     radiation at surface
               &              rnds_dir   = field%rnds_dir   (:,jb),           &! in, all-sky downward direct near-IR radiation at surface
               &              rvds_dif   = field%rvds_dif   (:,jb),           &! in, all-sky downward diffuse visible radiation at surface
               &              rpds_dif   = field%rpds_dif   (:,jb),           &! in, all-sky downward diffuse PAR     radiation at surface
               &              rnds_dif   = field%rnds_dif   (:,jb),           &! in, all-sky downward diffuse near-IR radiation at surface
               !
               &              ps = field% presi_old(:,nlevp1,jb),             &! in, paphm1, half level pressure
               &              pcosmu0 = field% cosmu0(:,jb),                  &! in, amu0_x, cos of zenith angle
               &              pch_tile = zch_tile(:,:),                       &! in, from "vdiff_down" for JSBACH
               &              pcsat = field%csat(:,jb),                       &! inout, area fraction with wet land surface
               &              pcair = field%cair(:,jb),                       &! inout, area fraction with wet land surface (air)
               &              q_snocpymlt = zq_snocpymlt(:),                  &! out, heating  by melting snow on the canopy [W/m2]
               &              z0m_tile = field% z0m_tile(:,jb,:),             &! inout, roughness length for momentum over tiles
               &              z0h_lnd  = field% z0h_lnd (:,jb),               &! out, roughness length for heat over land
               &              albvisdir      = field% albvisdir     (:,jb)  , &! inout
               &              albnirdir      = field% albnirdir     (:,jb)  , &! inout
               &              albvisdif      = field% albvisdif     (:,jb)  , &! inout
               &              albnirdif      = field% albnirdif     (:,jb)  , &! inout
               &              albvisdir_tile = field% albvisdir_tile(:,jb,:), &! inout
               &              albnirdir_tile = field% albnirdir_tile(:,jb,:), &! inout
               &              albvisdif_tile = field% albvisdif_tile(:,jb,:), &! inout
               &              albnirdif_tile = field% albnirdif_tile(:,jb,:), &! inout
               &              albedo         = field% albedo        (:,jb)  , &! inout
               &              albedo_tile    = field% albedo_tile(:,jb,:),    &! inout
               &              pco2_flux_tile = field% co2_flux_tile(:,jb,:),  &! inout
               &              ptsfc     = field%ts    (:,jb),                 &! out
               &              ptsfc_rad = field%ts_rad(:,jb),                 &! out
               &              rlns_tile = field%lwflxsfc_tile(:,jb,:),        &! out (for coupling)
               &              rsns_tile = field%swflxsfc_tile(:,jb,:),        &! out (for coupling)
               &              lake_ice_frc = field%lake_ice_frc(:,jb),        &! out
               &              Tsurf = field% Tsurf(:,:,jb),                   &! inout, for sea ice
               &              T1    = field% T1   (:,:,jb),                   &! inout, for sea ice
               &              T2    = field% T2   (:,:,jb),                   &! inout, for sea ice
               &              hi    = field% hi   (:,:,jb),                   &! in, for sea ice
               &              hs    = field% hs   (:,:,jb),                   &! in, for sea ice
               &              conc  = field% conc (:,:,jb),                   &! in, for sea ice
               &              Qtop  = field% Qtop (:,:,jb),                   &! out, for sea ice
               &              Qbot  = field% Qbot (:,:,jb),                   &! out, for sea ice
               &              albvisdir_ice = field% albvisdir_ice(:,:,jb),   &! inout ice albedos
               &              albnirdir_ice = field% albnirdir_ice(:,:,jb),   &! inout
               &              albvisdif_ice = field% albvisdif_ice(:,:,jb),   &! inout
               &              albnirdif_ice = field% albnirdif_ice(:,:,jb))    ! inout
          !
          !----------------------------------------------------------------------------------------
          ! Serialbox2 output fields serialization
          !$ser verbatim call serialize_update_surface_output(jb, jg, jce, nproma,&
          !$ser verbatim   nlev, nsfc_type, iwtr, iice, ilnd,&
          !$ser verbatim   pdtime, field, zaa, zaa_btm, zbb, zbb_btm,&
          !$ser verbatim   zcpt_sfc_tile, zq_snocpymlt)
          !
          IF (ltimer) CALL timer_stop(timer_vdf_sf)
          !
          !
          !
          ! Turbulent mixing, part II:
          ! - Elimination for the lowest model level using boundary conditions
          !   provided by the surface model(s);
          ! - Back substitution to get solution of the tridiagonal system;
          ! - Compute tendencies and additional diagnostics.
          !
          IF (ltimer) CALL timer_start(timer_vdf_up)
          !
          !----------------------------------------------------------------------------------------
          ! Serialbox2 input fields serialization
          !$ser verbatim call serialize_vdiff_up_input(jb, jce, nproma, nlev,&
          !$ser verbatim   nlevm1, ntrac, nsfc_type, iwtr, pdtime, field, zaa,&
          !$ser verbatim   zcptgz, ztottevn, zbb, zthvvar, dummyx)
          !
          CALL vdiff_up(jce, nproma, nlev, nlevm1,       &! in
               &        ntrac, nsfc_type,                &! in
               &        iwtr,                            &! in, indices of different sfc types
               &        pdtime,                          &! in, time steps
               &        field%frac_tile(:,jb,:),         &! in, area fraction of each sfc type
               &        field% cfm_tile(:,jb,:),         &! in
               &        zaa,                             &! in, from "vdiff_down"
               &         zcptgz(:,:),                    &! in, from "vdiff_down"
               &        field%   ua(:,:,jb),             &! in, um1
               &        field%   va(:,:,jb),             &! in, vm1
               &        field%   ta(:,:,jb),             &! in, tm1
               &        field% mair(:,:,jb),             &! in, moist     air mass [kg/m2]
               &        field% mref(:,:,jb),             &! in, reference air mass [kg/m2]
               &        field% qtrc(:,:,jb,iqv),         &! in, qm1
               &        field% qtrc(:,:,jb,iqc),         &! in, xlm1
               &        field% qtrc(:,:,jb,iqi),         &! in, xim1
               &        field% qtrc(:,:,jb,iqt:),        &! in, xtm1
               &        field% geom(:,:,jb),             &! in, pgeom1 = geopotential above ground
               &             ztottevn(:,:),              &! in, TTE at intermediate time step
               &        zbb,                             &! in
               &        zthvvar(:,:),                    &! inout
               &        dummyx(:,:),                     &! inout
               &        field% z0m_tile(:,jb,:),         &! inout
               &        field% kedisp(:,  jb),           &! out, vert. integr. diss. kin. energy [W/m2]
               &         tend%   ua_vdf(:,:,jb),         &! out
               &         tend%   va_vdf(:,:,jb),         &! out
               &        field% q_vdf   (:,:,jb),         &! out   heating W/m2
               &         tend% qtrc_vdf(:,:,jb,iqv),     &! out
               &         tend% qtrc_vdf(:,:,jb,iqc),     &! out
               &         tend% qtrc_vdf(:,:,jb,iqi),     &! out
               &         tend% qtrc_vdf(:,:,jb,iqt:),    &! out
               &        field%   z0m   (:,  jb),         &! out, for the next step
               &        dummy(:,:),                      &! 
               &        field%      totte(:,:,jb),       &! out
               &        field%   sh_vdiff(:,  jb),       &! out, for energy diagnostic
               &        field%   qv_vdiff(:,  jb)        )! out, for energy diagnostic
          !
          !----------------------------------------------------------------------------------------
          ! Serialbox2 output fields serialization
          !$ser verbatim call serialize_vdiff_up_output(jb, jce, nproma, nlev,&
          !$ser verbatim   nlevm1, ntrac, nsfc_type, iwtr, pdtime, field, zbb,&
          !$ser verbatim   dummyx, tend, dummy)
          !
          IF (ltimer) CALL timer_stop(timer_vdf_up)
          !
       END IF
       !
       !
       ! Surface effect on the lowermost layer
       !
       IF (echam_phy_config(jg)%ljsb) THEN
          !
          ! convert    heating
          ! zq_snocpymlt = heating for melting of snow on canopy
          !              = cooling of atmosphere --> negative sign
          tend% ta_sfc(jcs:jce,jb)      = -zq_snocpymlt(jcs:jce) * field% qconv(jcs:jce,nlev,jb)
          !
          ! accumulate heating
          field% q_phy(jcs:jce,nlev,jb) = field% q_phy(jcs:jce,nlev,jb) - zq_snocpymlt(jcs:jce)
          !
          ! accumulate tendencies for later updating the model state
          SELECT CASE(fc_vdf)
          CASE(0)
             ! diagnostic, do not use tendency
          CASE(1)
             ! use tendency to update the model state
             tend% ta_phy(jcs:jce,nlev,jb) = tend% ta_phy(jcs:jce,nlev,jb) + tend% ta_sfc(jcs:jce,jb)
!!$          CASE(2)
!!$             ! use tendency as forcing in the dynamics
!!$             ...
          END SELECT
          !
          ! update physics state for input to the next physics process
          IF (lparamcpl) THEN
             field% ta(jcs:jce,nlev,jb) = field% ta(jcs:jce,nlev,jb) + tend% ta_sfc(jcs:jce,jb)*pdtime
          END IF
          !
       END IF
       !
       !
       ! Vertical diffusion efect on the atmospheric column
       !
       ! convert    heating
       tend% ta_vdf(jcs:jce,:,jb) = field% q_vdf(jcs:jce,:,jb) * field% qconv(jcs:jce,:,jb)
       !
       ! accumulate heating
       field% q_phy(jcs:jce,:,jb) = field% q_phy(jcs:jce,:,jb) + field% q_vdf(jcs:jce,:,jb)
       !
       ! accumulate tendencies for later updating the model state
       SELECT CASE(fc_vdf)
       CASE(0)
          ! diagnostic, do not use tendency
       CASE(1)
          ! use tendency to update the model state
          tend%   ua_phy(jcs:jce,:,jb)      = tend%   ua_phy(jcs:jce,:,jb)      + tend%   ua_vdf(jcs:jce,:,jb)
          tend%   va_phy(jcs:jce,:,jb)      = tend%   va_phy(jcs:jce,:,jb)      + tend%   va_vdf(jcs:jce,:,jb)
          tend%   ta_phy(jcs:jce,:,jb)      = tend%   ta_phy(jcs:jce,:,jb)      + tend%   ta_vdf(jcs:jce,:,jb)
          tend% qtrc_phy(jcs:jce,:,jb,iqv)  = tend% qtrc_phy(jcs:jce,:,jb,iqv)  + tend% qtrc_vdf(jcs:jce,:,jb,iqv)
          tend% qtrc_phy(jcs:jce,:,jb,iqc)  = tend% qtrc_phy(jcs:jce,:,jb,iqc)  + tend% qtrc_vdf(jcs:jce,:,jb,iqc)
          tend% qtrc_phy(jcs:jce,:,jb,iqi)  = tend% qtrc_phy(jcs:jce,:,jb,iqi)  + tend% qtrc_vdf(jcs:jce,:,jb,iqi)
          tend% qtrc_phy(jcs:jce,:,jb,iqt:) = tend% qtrc_phy(jcs:jce,:,jb,iqt:) + tend% qtrc_vdf(jcs:jce,:,jb,iqt:)
!!$       CASE(2)
!!$          ! use tendency as forcing in the dynamics
!!$          ...
       END SELECT
       !
       ! update physics state for input to the next physics process
       IF (lparamcpl) THEN
          field%   ua(jcs:jce,:,jb)      = field%   ua(jcs:jce,:,jb)      + tend%   ua_vdf(jcs:jce,:,jb)     *pdtime
          field%   va(jcs:jce,:,jb)      = field%   va(jcs:jce,:,jb)      + tend%   va_vdf(jcs:jce,:,jb)     *pdtime
          field%   ta(jcs:jce,:,jb)      = field%   ta(jcs:jce,:,jb)      + tend%   ta_vdf(jcs:jce,:,jb)     *pdtime
          field% qtrc(jcs:jce,:,jb,iqv)  = field% qtrc(jcs:jce,:,jb,iqv)  + tend% qtrc_vdf(jcs:jce,:,jb,iqv) *pdtime
          field% qtrc(jcs:jce,:,jb,iqc)  = field% qtrc(jcs:jce,:,jb,iqc)  + tend% qtrc_vdf(jcs:jce,:,jb,iqc) *pdtime
          field% qtrc(jcs:jce,:,jb,iqi)  = field% qtrc(jcs:jce,:,jb,iqi)  + tend% qtrc_vdf(jcs:jce,:,jb,iqi) *pdtime
          field% qtrc(jcs:jce,:,jb,iqt:) = field% qtrc(jcs:jce,:,jb,iqt:) + tend% qtrc_vdf(jcs:jce,:,jb,iqt:)*pdtime
       END IF
       !
       !
       IF (echam_phy_config(jg)%lnew) THEN
       ! Correction related to implicitness, due to the fact that surface model only used
       ! part of longwave radiation to compute new surface temperature
       ! 
       zq_rlw_impl(jcs:jce) =                                                   &
            &  ( (field%rld_rt(jcs:jce,nlev,jb)-field%rlu_rt(jcs:jce,nlev,jb))  & ! ( rln  from "radiation", at top of layer nlev
            &   -(field%rlds  (jcs:jce,jb)     -field%rlus  (jcs:jce,jb)     )) & !  -rlns from "radheating" and "update_surface")
            & -field%q_rlw(jcs:jce,nlev,jb)                                       ! -old heating in layer nlev from "radheating"
       !
       ! convert    heating
       tend%ta_rlw_impl(jcs:jce,jb) = zq_rlw_impl(jcs:jce) * field% qconv(jcs:jce,nlev,jb)
       !
       ! accumulate heating
       field% q_phy(jcs:jce,nlev,jb) = field% q_phy(jcs:jce,nlev,jb) + zq_rlw_impl(jcs:jce)
       !
       ! accumulate tendencies for later updating the model state
       SELECT CASE(fc_vdf)
       CASE(0)
          ! diagnostic, do not use tendency
       CASE(1)
          ! use tendency to update the model state
          tend%ta_phy(jcs:jce,nlev,jb) = tend%ta_phy(jcs:jce,nlev,jb) + tend%ta_rlw_impl(jcs:jce,jb)
!!$       CASE(2)
!!$          ! use tendency as forcing in the dynamics
!!$          ...
       END SELECT
       !
       ! update physics state for input to the next physics process
       IF (lparamcpl) THEN
          field% ta(jcs:jce,nlev,jb) = field% ta(jcs:jce,nlev,jb) + tend% ta_rlw_impl(jcs:jce,jb)*pdtime
       END IF
       !
       END IF

       ! 2-tl-scheme
       field% tottem1(jcs:jce,:,jb) = field% totte (jcs:jce,:,jb)


       ! Turbulent mixing, part III:
       ! - Further diagnostics.
       !----------------------------------------------------------------------------------------
       ! Serialbox2 input fields serialization
       !$ser verbatim call serialize_nsurf_diag_input(jb, jce, nproma, nsfc_type,&
       !$ser verbatim   ilnd, field, zqx(:,nlev), zcptgz(:,nlev), zcpt_sfc_tile, zbn_tile,&
       !$ser verbatim   zbhn_tile, zbh_tile, zbm_tile, zri_tile)
       !
       CALL nsurf_diag(jce, nproma, nsfc_type,          &! in
            &          ilnd,                            &! in
            &          field%frac_tile(:,jb,:),         &! in
            &          field%  qtrc(:,nlev,jb,iqv),     &! in humidity qm1
            &          field%    ta(:,nlev,jb),         &! in tm1
            &          field% presm_old(:,nlev,jb),     &! in, apm1
            &          field% presi_old(:,nlevp1,jb),   &! in, aphm1
            &                  zqx(:,nlev),             &! in, xlm1 + xim1
            &          field%   ua(:,nlev,jb),          &! in, um1
            &          field%   va(:,nlev,jb),          &! in, vm1
            &          field% ocu (:,jb),               &! in, ocean sfc velocity, u-component
            &          field% ocv (:,jb),               &! in, ocean sfc velocity, v-component
            &          field% zf  (:,nlev  ,jb),        &! in, height of lowermost full level (m)
            &          field% zh  (:,nlev+1,jb),        &! in, surface height    (m)
            &          zcptgz(:,nlev),                  &! in dry static energy
            &          zcpt_sfc_tile(:,:),              &! in dry static energy
            &          zbn_tile(:,:),                   &! in for diagnostic
            &          zbhn_tile(:,:),                  &! in for diagnostic
            &          zbh_tile(:,:),                   &! in for diagnostic
            &          zbm_tile(:,:),                   &! in for diagnostic
            &          zri_tile(:,:),                   &! in 
            &          field%sfcWind(:,  jb),           &! out 10m windspeed
            &          field%    tas(:,  jb),           &! out temperature in 2m
            &          field%   dew2(:,  jb),           &! out dew point temperature in 2m
            &          field%    uas(:,  jb),           &! out zonal wind in 10m
            &          field%    vas(:,  jb),           &! out meridional wind in 10m
            &          field%tasmax (:,  jb),           &! out max 2m temperature
            &          field%tasmin (:,  jb),           &! out min 2m temperature
            &          field%sfcWind_tile(:,jb,:),      &! out 10m windspeed on tiles
            &          field%    tas_tile(:,jb,:),      &! out temperature in 2m on tiles
            &          field%   dew2_tile(:,jb,:),      &! out dew point temperature in 2m on tiles
            &          field%    uas_tile(:,jb,:),      &! out zonal wind in 10m on tiles
            &          field%    vas_tile(:,jb,:)       )! out meridional wind in 10m on tiles
       !
       !----------------------------------------------------------------------------------------
       ! Serialbox2 output fields serialization
       !$ser verbatim call serialize_nsurf_diag_output(jb, jce, nproma,&
       !$ser verbatim   nsfc_type, ilnd, field)


       !
       !
    ELSE
       !
       !
       ! vdiff_down
       field% ustar          (jcs:jce,  jb  ) = 0.0_wp
       field% wstar          (jcs:jce,  jb  ) = 0.0_wp
       field% wstar_tile     (jcs:jce,  jb,:) = 0.0_wp
       field% qs_sfc_tile    (jcs:jce,  jb,:) = 0.0_wp
       field% hdtcbl         (jcs:jce,  jb  ) = 0.0_wp
       field% ri             (jcs:jce,:,jb  ) = 0.0_wp
       field% mixlen         (jcs:jce,:,jb  ) = 0.0_wp
       field% cfm            (jcs:jce,:,jb  ) = 0.0_wp
       field% cfm_tile       (jcs:jce,  jb,:) = 0.0_wp
       field% cfh            (jcs:jce,:,jb  ) = 0.0_wp
       field% cfh_tile       (jcs:jce,  jb,:) = 0.0_wp
       field% cfv            (jcs:jce,:,jb  ) = 0.0_wp
       field% cftotte        (jcs:jce,:,jb  ) = 0.0_wp
       field% cfthv          (jcs:jce,:,jb  ) = 0.0_wp
       field% thvsig         (jcs:jce,  jb  ) = 0.0_wp
       !
       ! update_surface
       field% ts_tile        (jcs:jce,  jb,:) = 0.0_wp
       field% u_stress       (jcs:jce,  jb  ) = 0.0_wp
       field% v_stress       (jcs:jce,  jb  ) = 0.0_wp
       field% lhflx          (jcs:jce,  jb  ) = 0.0_wp
       field% shflx          (jcs:jce,  jb  ) = 0.0_wp
       field%  evap          (jcs:jce,  jb  ) = 0.0_wp
       field% u_stress_tile  (jcs:jce,  jb,:) = 0.0_wp
       field% v_stress_tile  (jcs:jce,  jb,:) = 0.0_wp
       field% lhflx_tile     (jcs:jce,  jb,:) = 0.0_wp
       field% shflx_tile     (jcs:jce,  jb,:) = 0.0_wp
       field%  evap_tile     (jcs:jce,  jb,:) = 0.0_wp
       field% rlus           (jcs:jce,  jb  ) = 0.0_wp
       field% rsus           (jcs:jce,  jb  ) = 0.0_wp
       field% csat           (jcs:jce,  jb  ) = 0.0_wp
       field% cair           (jcs:jce,  jb  ) = 0.0_wp
       field% z0h_lnd        (jcs:jce,  jb  ) = 0.0_wp
       field% albvisdir      (jcs:jce,  jb  ) = 0.0_wp
       field% albnirdir      (jcs:jce,  jb  ) = 0.0_wp
       field% albvisdif      (jcs:jce,  jb  ) = 0.0_wp
       field% albnirdif      (jcs:jce,  jb  ) = 0.0_wp
       field% albvisdir_tile (jcs:jce,  jb,:) = 0.0_wp
       field% albnirdir_tile (jcs:jce,  jb,:) = 0.0_wp
       field% albvisdif_tile (jcs:jce,  jb,:) = 0.0_wp
       field% albnirdif_tile (jcs:jce,  jb,:) = 0.0_wp
       field% albedo         (jcs:jce,  jb  ) = 0.0_wp
       field% albedo_tile    (jcs:jce,  jb,:) = 0.0_wp
       field% co2_flux_tile  (jcs:jce,  jb,:) = 0.0_wp
       field% ts             (jcs:jce,  jb  ) = 0.0_wp
       field% ts_rad         (jcs:jce,  jb  ) = 0.0_wp
       field% lwflxsfc_tile  (jcs:jce,  jb,:) = 0.0_wp
       field% swflxsfc_tile  (jcs:jce,  jb,:) = 0.0_wp
       !
       field% Tsurf          (jcs:jce,:,jb  ) = 0.0_wp
       field% T1             (jcs:jce,:,jb  ) = 0.0_wp
       field% T2             (jcs:jce,:,jb  ) = 0.0_wp
       field% Qtop           (jcs:jce,:,jb  ) = 0.0_wp
       field% Qbot           (jcs:jce,:,jb  ) = 0.0_wp
       field% albvisdir_ice  (jcs:jce,:,jb  ) = 0.0_wp
       field% albnirdir_ice  (jcs:jce,:,jb  ) = 0.0_wp
       field% albvisdif_ice  (jcs:jce,:,jb  ) = 0.0_wp
       field% albnirdif_ice  (jcs:jce,:,jb  ) = 0.0_wp
       !
       tend% ta_sfc          (jcs:jce,  jb  ) = 0.0_wp
       !
       ! vdiff_up
       field% totte          (jcs:jce,:,jb  ) = 0.0_wp
       field% z0m            (jcs:jce,  jb  ) = 0.0_wp
       field% z0m_tile       (jcs:jce,  jb,:) = 0.0_wp
       field% kedisp         (jcs:jce,  jb  ) = 0.0_wp
       field% sh_vdiff       (jcs:jce,  jb  ) = 0.0_wp
       field% qv_vdiff       (jcs:jce,  jb  ) = 0.0_wp
       !
       tend%   ua_vdf        (jcs:jce,:,jb  ) = 0.0_wp
       tend%   va_vdf        (jcs:jce,:,jb  ) = 0.0_wp
       tend%   ta_vdf        (jcs:jce,:,jb  ) = 0.0_wp
       tend% qtrc_vdf        (jcs:jce,:,jb,:) = 0.0_wp
       !
    END IF

    ! disassociate pointers
    NULLIFY(lparamcpl)
    NULLIFY(fc_vdf)
    NULLIFY(vmr_co2)
    NULLIFY(field)
    NULLIFY(tend)

    IF (ltimer) CALL timer_stop(timer_vdf)

  END SUBROUTINE interface_echam_vdf

END MODULE mo_interface_echam_vdf
