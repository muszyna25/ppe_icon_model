!>
!! @brief Subroutine echam_phy_main calls all the parameterization schemes
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
#if defined __xlC__ && !defined NOXLFPROCESS
@PROCESS HOT
@PROCESS SPILLSIZE(5000)
#endif
!OCL NOALIAS

MODULE mo_echam_phy_main

  USE mo_kind,                ONLY: wp
  USE mtime,                  ONLY: datetime
  !
  USE mo_math_constants,      ONLY: pi
  USE mo_physical_constants,  ONLY: cpd, cpv, cvd, cvv, amd, amo3
  !
  USE mo_run_config,          ONLY: ntracer, nlev, nlevm1, nlevp1,    &
    &                               iqv, iqc, iqi, iqt, io3
  !
  USE mo_echam_phy_config,    ONLY: echam_phy_config
  USE mo_echam_phy_memory,    ONLY: t_echam_phy_field, prm_field,     &
    &                               t_echam_phy_tend,  prm_tend
  !
  USE mo_cover,               ONLY: cover
  !
  USE mo_ext_data_state,      ONLY: ext_data
  USE mo_psrad_radiation_parameters, ONLY: psctm
  USE mo_psrad_radiation,     ONLY: psrad_radiation
  USE mo_radheating,          ONLY: radheating
  !
  USE mo_vdiff_config,        ONLY: vdiff_config
  USE mo_vdiff_downward_sweep,ONLY: vdiff_down
  USE mo_vdiff_upward_sweep,  ONLY: vdiff_up
  USE mo_vdiff_solver,        ONLY: nvar_vdiff, nmatrix, imh, imqv,   &
    &                               ih_vdiff=>ih, iqv_vdiff=>iqv
  !
  USE mo_echam_sfc_indices,   ONLY: nsfc_type, iwtr, iice, ilnd
  USE mo_surface,             ONLY: update_surface
  USE mo_surface_diag,        ONLY: nsurf_diag
  !
  USE mo_bcs_time_interpolation, ONLY: t_time_interpolation_weights, &
    &                                  calculate_time_interpolation_weights
  USE mo_lcariolle_types,     ONLY: avi, t_time_interpolation
  !
  USE mo_gw_hines,            ONLY: gw_hines
  USE mo_ssortns,             ONLY: ssodrag
  !
  USE mo_echam_conv_config,   ONLY: echam_conv_config
  USE mo_cumastr,             ONLY: cumastr
  !
  USE mo_echam_cloud_config,  ONLY: echam_cloud_config
  USE mo_cloud,               ONLY: cloud
  !
  USE mo_timer,               ONLY: ltimer, timer_start, timer_stop,                &
    &                               timer_cover, timer_radiation, timer_radheat,    &
    &                               timer_vdiff_down, timer_surface,timer_vdiff_up, &
    &                               timer_gw_hines, timer_ssodrag,                  &
    &                               timer_convection, timer_cloud

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: echam_phy_main

CONTAINS
  !>
  !!
  SUBROUTINE echam_phy_main( jg,jb,jcs,jce,nbdim,      &
    &                        this_datetime,pdtime,     &
    &                        ltrig_rad                 )

    INTEGER         ,INTENT(IN) :: jg             !< grid level/domain index
    INTEGER         ,INTENT(IN) :: jb             !< block index
    INTEGER         ,INTENT(IN) :: jcs, jce       !< start/end column index within this block
    INTEGER         ,INTENT(IN) :: nbdim          !< size of this block

    TYPE(datetime), POINTER     :: this_datetime  !< date and time
    REAL(wp)        ,INTENT(IN) :: pdtime         !< time step

    LOGICAL         ,INTENT(IN) :: ltrig_rad      !< perform radiative transfer computation

    ! Local variables

    TYPE(t_echam_phy_field),   POINTER :: field
    TYPE(t_echam_phy_tend) ,   POINTER :: tend

    REAL(wp) :: zlat_deg(nbdim)           !< latitude in deg N

    REAL(wp) :: zvmixtau   (nbdim,nlev)   !< timescale of mixing for vertical turbulence
    REAL(wp) :: zqtvar_prod(nbdim,nlev)   !< production rate of total water variance
                                          !< due to turbulence. Computed in "vdiff",
                                          !< used by "cloud"
    INTEGER  :: itype(nbdim)              !< type of convection
    INTEGER  :: ictop (nbdim)             !< from massflux

    REAL(wp) :: zfrl (nbdim)              !< fraction of land in the grid box
    REAL(wp) :: zfrw (nbdim)              !< fraction of water (without ice) in the grid point
    REAL(wp) :: zfri (nbdim)              !< fraction of ice in the grid box
    REAL(wp) :: zfrc (nbdim,nsfc_type)    !< zfrl, zfrw, zfrc combined
    REAL(wp) :: zri_tile(nbdim,nsfc_type) !< Richardson number

!    REAL(wp) :: zcvcbot(nbdim)
!    REAL(wp) :: zwcape (nbdim)

    REAL(wp) :: zta    (nbdim,nlev)         !< provisional temperature         [K]
    REAL(wp) :: zqtrc  (nbdim,nlev,ntracer) !< provisional mass mixing ratios  [kg/kg]
    REAL(wp) :: zua    (nbdim,nlev)         !< provisional zonal      wind     [m/s]
    REAL(wp) :: zva    (nbdim,nlev)         !< provisional meridional wind     [m/s]

    REAL(wp) :: zcpair (nbdim,nlev)       !< specific heat of moist air at const. pressure [J/K/kg]
    REAL(wp) :: zcvair (nbdim,nlev)       !< specific heat of moist air at const. volume   [J/K/kg]
    REAL(wp) :: zconv  (nbdim,nlev)       !< conversion factor q-->dT/dt       [(K/s)/(W/m2)]
    REAL(wp) :: zdelp  (nbdim,nlev)       !< layer thickness in pressure coordinate  [Pa]

    REAL(wp) :: zq_phy (nbdim,nlev)       !< heating by whole ECHAM physics    [W/m2]
    REAL(wp) :: zq_rsw (nbdim,nlev)       !< heating by short wave radiation   [W/m2]
    REAL(wp) :: zq_rlw (nbdim,nlev)       !< heating by long  wave radiation   [W/m2]
    REAL(wp) :: zq_rlw_impl (nbdim)       !< additional heating by LW rad. due to impl. coupling in surface energy balance [W/m2]
!!$    REAL(wp) :: zq_vdf (nbdim,nlev)       !< heating by vertical diffusion     [W/m2]
    REAL(wp) :: zq_sso (nbdim,nlev)       !< heating by subgrid scale orogr.   [W/m2]
    REAL(wp) :: zq_gwh (nbdim,nlev)       !< heating by atm. gravity waves     [W/m2]
    REAL(wp) :: zq_cnv (nbdim,nlev)       !< heating by convection             [W/m2]
    REAL(wp) :: zq_cld (nbdim,nlev)       !< heating by stratiform clouds      [W/m2]

    INTEGER  :: ihpbl  (nbdim)            !< location of PBL top given as vertical level index
    REAL(wp) :: zxt_emis(nbdim,ntracer-iqt+1)  !< tracer tendency due to surface emission
                                               !< and dry deposition. "zxtems" in ECHAM5
    INTEGER  :: jk
    INTEGER  :: jks   !< start index for vertical loops
    INTEGER  :: nc    !< number of cells/columns from (jce-jcs+1)
    INTEGER  :: jc
    INTEGER  :: ntrac !< # of tracers excluding water vapour and hydrometeors
                      !< (handled by sub-models, e.g., chemical species)

    ! Coefficient matrices and right-hand-side vectors for the turbulence solver
    ! _btm refers to the lowest model level (i.e., full level "klev", not the surface)

    REAL(wp) :: zaa    (nbdim,nlev,3,nmatrix)       !< coeff. matrices, all variables
    REAL(wp) :: zaa_btm(nbdim,3,nsfc_type,imh:imqv) !< last row of coeff. matrix of heat and moisture
    REAL(wp) :: zbb    (nbdim,nlev,nvar_vdiff)  !< r.h.s., all variables
    REAL(wp) :: zbb_btm(nbdim,nsfc_type,ih_vdiff:iqv_vdiff) !< last row of r.h.s. of heat and moisture

    ! Temporary arrays used by VDIFF, JSBACH

    REAL(wp) :: zfactor_sfc(nbdim)
    REAL(wp) :: zcpt_sfc_tile(nbdim,nsfc_type)  !< dry static energy at surface

    REAL(wp) :: zcptgz   (nbdim,nlev) !< dry static energy
    REAL(wp) :: zrhoh    (nbdim,nlev) !< air density at half levels
    REAL(wp) :: zqshear  (nbdim,nlev) !<
    REAL(wp) :: zthvvar  (nbdim,nlev) !< intermediate value of thvvar
    REAL(wp) :: ztkevn   (nbdim,nlev) !< intermediate value of tke
    REAL(wp) :: zch_tile (nbdim,nsfc_type)   !<  for "nsurf_diag"
!    REAL(wp) :: zchn_tile(nbdim,nsfc_type)   !<  for "nsurf_diag"
!    REAL(wp) :: zcdn_tile(nbdim,nsfc_type)   !<  for "nsurf_diag"
!    REAL(wp) :: zcfnc_tile(nbdim,nsfc_type)  !<  for "nsurf_diag"
    REAL(wp) :: zbn_tile (nbdim,nsfc_type)   !<  for "nsurf_diag"
    REAL(wp) :: zbhn_tile(nbdim,nsfc_type)   !<  for "nsurf_diag"
    REAL(wp) :: zbm_tile (nbdim,nsfc_type)   !<  for "nsurf_diag"
    REAL(wp) :: zbh_tile (nbdim,nsfc_type)   !<  for "nsurf_diag"

    REAL(wp) :: ztte_corr(nbdim)      !< tte correction for snow melt over land (JSBACH)

    ! Temporary variables for Cariolle scheme (ozone)
    REAL(wp)    :: do3dt(nbdim,nlev)
    TYPE(t_time_interpolation) :: time_interpolation
    EXTERNAL       lcariolle_lat_intp_li, lcariolle_pres_intp_li
    TYPE(t_time_interpolation_weights) :: current_time_interpolation_weights
 
    ! Temporary array used by GW_HINES

    REAL(wp) :: zdis_gwh(nbdim,nlev)  !<  out, energy dissipation rate [J/s/kg]


    ! Temporary array used by SSODRAG

    REAL(wp) :: zdis_sso(nbdim,nlev)  !<  out, energy dissipation rate [J/s/kg]


    ! Temporary array used by convection

    REAL(wp) :: zqtrc_cnd(nbdim,nlev) !<  cloud condensate mixing ratio [kg/kg]
    REAL(wp) :: ztend_qv(nbdim,nlev)  !<  moisture tendency from dynamics and physics before convection
    REAL(wp) :: ztop(nbdim)           !<  convective cloud top pressure [Pa]

    ! Temporary variables used for cloud droplet number concentration

    REAL(wp) :: zprat, zn1, zn2, zcdnc
    LOGICAL  :: lland(nbdim), lglac(nbdim)

    ! number of cells/columns from index jcs to jce
    nc = jce-jcs+1

    ! start index for vertical loops
    jks=1

    ! 1. Associate pointers

    field  => prm_field(jg)
    tend   => prm_tend (jg)

    ! provisionally copy the incoming tedencies

    tend%   ta_phy (jcs:jce,:,jb)   = tend%   ta (jcs:jce,:,jb)
    tend%   ua_phy (jcs:jce,:,jb)   = tend%   ua (jcs:jce,:,jb)
    tend%   va_phy (jcs:jce,:,jb)   = tend%   va (jcs:jce,:,jb)
    tend% qtrc_phy (jcs:jce,:,jb,:) = tend% qtrc (jcs:jce,:,jb,:)

    ! initialize physics heating
    zq_phy(:,:) = 0._wp

    ! 2. local switches and parameters

    ntrac = ntracer-iqt+1  !# of tracers excluding water vapour and hydrometeors

    !------------------------------------------------------------
    ! 3. COMPUTE SOME FIELDS NEEDED BY THE PHYSICAL ROUTINES.
    !------------------------------------------------------------

    ! 3.2b Specific heat of moist air
    !
    zcpair  (:,:) = cpd+(cpv-cpd)*field%qtrc(:,:,jb,iqv)
    zcvair  (:,:) = cvd+(cvv-cvd)*field%qtrc(:,:,jb,iqv)
    !
    zconv   (:,:) = 1._wp/(field%mair(:,:,jb)*zcpair(:,:))

    DO jk = 1,nlev
      !
      ! 3.2 Thickness of model layer in pressure coordinate
      !
      zdelp   (:,jk) = field% presi_old (:,jk+1,jb) - field% presi_old (:,jk,jb)
      !
    END DO

    ! 3.3 Weighting factors for fractional surface coverage
    !     Accumulate ice portion for diagnostics

    DO jc=jcs,jce

      ! fraction of land in the grid box. lsmask: land-sea mask, 1.= land

      ! TBD: use fractional mask here
      zfrl(jc) = field% lsmask(jc,jb)

      ! fraction of sea/lake in the grid box
      ! * (1. - fraction of sea ice in the sea/lake part of the grid box)
      ! => fraction of open water in the grid box

      zfrw(jc) = (1._wp-zfrl(jc))*(1._wp-field%seaice(jc,jb))

      ! fraction of sea ice in the grid box
      zfri(jc) = 1._wp-zfrl(jc)-zfrw(jc)
    END DO

    ! 3.4 Merge three pieces of information into one array for vdiff
    IF (ilnd.LE.nsfc_type) zfrc(jcs:jce,ilnd) = zfrl(jcs:jce)
    IF (iwtr.LE.nsfc_type) zfrc(jcs:jce,iwtr) = zfrw(jcs:jce)
    IF (iice.LE.nsfc_type) zfrc(jcs:jce,iice) = zfri(jcs:jce)

    !---------------------------------------------------------------------
    ! 3.9 DETERMINE TROPOPAUSE HEIGHT AND MASS BUDGETS
    !     (Needed only for sub-models. Note: sequence of arguments
    !      different from the original ECHAM6)
    !---------------------------------------------------------------------
    !
    !CALL WMO_tropopause( jce, nbdim, nlev,         &! in
    !                   & ncctop, nccbot, lresum,   &! in
    !                   & field% ta(:,:,jb),        &! in
    !                   & field% presm_old(:,:,jb), &! in
    !                   & field% tropo(:,jb),       &! out for diagnostics
    !                   & itrpwmo, itrpwmop1        )! out for submodel

    !---------------------------------------------------------------------
    ! 3.12 INITIALISATION OF CLOUD DROPLET NUMBER CONCENTRATION 
    !      (1/M**3) USED IN RADLSW AND CLOUD
    !---------------------------------------------------------------------

    DO jc=jcs,jce
      lland(jc) = field%lfland(jc,jb)
      lglac(jc) = lland(jc).AND.field%glac(jc,jb).GT.0._wp
    END DO

    DO jk = 1,nlev
      DO jc = jcs,jce
        !
        zprat=(MIN(8._wp,80000._wp/field%presm_old(jc,jk,jb)))**2

        IF (lland(jc).AND.(.NOT.lglac(jc))) THEN
          zn1= echam_cloud_config% cn1lnd
          zn2= echam_cloud_config% cn2lnd
        ELSE
          zn1= echam_cloud_config% cn1sea
          zn2= echam_cloud_config% cn2sea
        ENDIF
        IF (field%presm_old(jc,jk,jb).LT.80000._wp) THEN
          zcdnc=1.e6_wp*(zn1+(zn2-zn1)*(EXP(1._wp-zprat)))
        ELSE
          zcdnc=zn2*1.e6_wp
        ENDIF
        field% acdnc(jc,jk,jb) = zcdnc
        !
      END DO
    END DO

    !-------------------------------------------------------------------
    ! 3.13 DIAGNOSE CURRENT CLOUD COVER
    !-------------------------------------------------------------------

    itype(jcs:jce) = NINT(field%rtype(jcs:jce,jb))

    IF (echam_phy_config%lcond) THEN
      IF (ltimer) CALL timer_start(timer_cover)

      CALL cover( jce, nbdim, jks,          &! in
        &         nlev, nlevp1,             &! in
        &         itype,  zfrw, zfri,       &! in
        &         field% zf(:,:,jb),        &! in
        &         field% presi_old(:,:,jb), &! in
        &         field% presm_old(:,:,jb), &! in
        &         field%  ta(:,:,jb),       &! in    tm1
        &         field%  qtrc(:,:,jb,iqv), &! in    qm1
        &         field%  qtrc(:,:,jb,iqi), &! in    xim1
        &         field%  aclc(:,:,jb),     &! out   (for "radiation" and "vdiff_down")
        &         field% rintop(:,  jb)    ) ! out   (for output)

      IF (ltimer) CALL timer_stop(timer_cover)
    ENDIF ! lcond

    !-------------------------------------------------------------------
    ! 4. RADIATION PARAMETERISATION
    !-------------------------------------------------------------------
    IF (echam_phy_config%lrad) THEN

       ! 4.1 RADIATIVE TRANSFER
       !-----------------------
       IF (ltrig_rad) THEN

          ! store ts_rad of this radiatiative transfer timestep in ts_rad_rt,
          ! so that it can be reused in radheat in the other timesteps
          field%ts_rad_rt(jcs:jce,jb) = field%ts_rad(jcs:jce,jb)

        IF (ltimer) CALL timer_start(timer_radiation)

        CALL psrad_radiation(                       &
        & jg                                       ,&!< in  domain index
        & jb                                       ,&!< in  block index
        & kproma         = jce                     ,&!< in  end index for loop over block
        & kbdim          = nbdim                   ,&!< in  dimension of block over cells
        & klev           = nlev                    ,&!< in  number of full levels = number of layers
        & klevp1         = nlevp1                  ,&!< in  number of half levels = number of layer interfaces
        & ktype          = itype(:)                ,&!< in  type of convection
        & loland         = lland                   ,&!< in  land-sea mask. (logical)
        & loglac         = lglac                   ,&!< in  glacier mask (logical)
        & this_datetime  = this_datetime           ,&!< in  actual time step
        & pcos_mu0       = field%cosmu0_rt(:,jb)   ,&!< in  solar zenith angle
        & alb_vis_dir    = field%albvisdir(:,jb)   ,&!< in  surface albedo for visible range, direct
        & alb_nir_dir    = field%albnirdir(:,jb)   ,&!< in  surface albedo for near IR range, direct
        & alb_vis_dif    = field%albvisdif(:,jb)   ,&!< in  surface albedo for visible range, diffuse
        & alb_nir_dif    = field%albnirdif(:,jb)   ,&!< in  surface albedo for near IR range, diffuse
        & tk_sfc         = field%ts_rad_rt(:,jb)   ,&!< in  grid box mean surface temperature
        & zf             = field%zf(:,:,jb)        ,&!< in  geometric height at full level      [m]
        & zh             = field%zh(:,:,jb)        ,&!< in  geometric height at half level      [m]
        & dz             = field%dz(:,:,jb)        ,&!< in  geometric height thickness of layer [m]
        & pp_hl          = field%presi_old(:,:,jb) ,&!< in  pressure at half levels at t-dt [Pa]
        & pp_fl          = field%presm_old(:,:,jb) ,&!< in  pressure at full levels at t-dt [Pa]
        & tk_fl          = field%ta(:,:,jb)        ,&!< in  tk_fl  = temperature at full level at t-dt
        & xm_dry         = field%mdry(:,:,jb)      ,&!< in  dry air mass in layer [kg/m2]
        & xm_trc         = field%mtrc(:,:,jb,:)    ,&!< in  tracer  mass in layer [kg/m2]
        & xm_ozn         = field%o3(:,:,jb)        ,&!< inout  ozone  mass mixing ratio [kg/kg]
        !
        & cdnc           = field% acdnc(:,:,jb)    ,&!< in   cloud droplet number conc
        & cld_frc        = field% aclc(:,:,jb)     ,&!< in   cloud fraction [m2/m2]
        & cld_cvr        = field%aclcov(:,jb)      ,&!< out  total cloud cover
        !
        & lw_dnw_clr     = field%rldcs_rt(:,:,jb)  ,&!< out  Clear-sky net longwave  at all levels
        & lw_upw_clr     = field%rlucs_rt(:,:,jb)  ,&!< out  Clear-sky net longwave  at all levels
        & sw_dnw_clr     = field%rsdcs_rt(:,:,jb)  ,&!< out  Clear-sky net shortwave at all levels
        & sw_upw_clr     = field%rsucs_rt(:,:,jb)  ,&!< out  Clear-sky net shortwave at all levels
        & lw_dnw         = field%rld_rt  (:,:,jb)  ,&!< out  All-sky net longwave  at all levels
        & lw_upw         = field%rlu_rt  (:,:,jb)  ,&!< out  All-sky net longwave  at all levels
        & sw_dnw         = field%rsd_rt  (:,:,jb)  ,&!< out  All-sky net longwave  at all levels
        & sw_upw         = field%rsu_rt  (:,:,jb)  ,&!< out  All-sky net longwave  at all levels
        !
        & vis_dn_dir_sfc = field%rvds_dir_rt(:,jb) ,&!< out  all-sky downward direct visible radiation at surface
        & par_dn_dir_sfc = field%rpds_dir_rt(:,jb) ,&!< all-sky downward direct PAR     radiation at surface
        & nir_dn_dir_sfc = field%rnds_dir_rt(:,jb) ,&!< all-sky downward direct near-IR radiation at surface
        & vis_dn_dff_sfc = field%rvds_dif_rt(:,jb) ,&!< all-sky downward diffuse visible radiation at surface
        & par_dn_dff_sfc = field%rpds_dif_rt(:,jb) ,&!< all-sky downward diffuse PAR     radiation at surface
        & nir_dn_dff_sfc = field%rnds_dif_rt(:,jb) ,&!< all-sky downward diffuse near-IR radiation at surface
        & vis_up_sfc     = field%rvus_rt    (:,jb) ,&!< all-sky upward visible radiation at surface
        & par_up_sfc     = field%rpus_rt    (:,jb) ,&!< all-sky upward PAR     radiation at surfac
        & nir_up_sfc     = field%rnus_rt    (:,jb) ) !< all-sky upward near-IR radiation at surface

        
        IF (ltimer) CALL timer_stop(timer_radiation)

      END IF ! ltrig_rad

      ! 4.2 RADIATIVE HEATING
      !----------------------

      ! radheat first computes the shortwave and longwave radiation for the current time step from transmissivity and
      ! the longwave flux at the radiation time step and, from there, the radiative heating due to sw and lw radiation.
      ! If radiation is called every time step, the longwave flux is not changed.

      IF (ltimer) CALL timer_start(timer_radheat)

      CALL radheating (                                &
        !
        ! input
        ! -----
        !
        & jcs        = jcs                            ,&! loop start index
        & jce        = jce                            ,&! loop end index
        & kbdim      = nbdim                          ,&! dimension size
        & klev       = nlev                           ,&! vertical dimension size
        & klevp1     = nlevp1                         ,&! vertical dimension size
        !
        & rsdt0      = psctm                          ,&! toa incident shortwave radiation for sun in zenith
        & cosmu0     = field%cosmu0    (:,jb)         ,&! solar zenith angle at current time
        !
        & emiss      = ext_data(jg)%atm%emis_rad(:,jb),&! lw sfc emissivity
        & tsr        = field%ts_rad (:,jb)            ,&! radiative surface temperature at current   time [K]
        & tsr_rt     = field%ts_rad_rt(:,jb)          ,&! radiative surface temperature at radiation time [K]
        !
        & rsd_rt     = field%rsd_rt           (:,:,jb),&! all-sky   shortwave downward flux at radiation time [W/m2]
        & rsu_rt     = field%rsu_rt           (:,:,jb),&! all-sky   shortwave upward   flux at radiation time [W/m2]
        !
        & rsdcs_rt   = field%rsdcs_rt         (:,:,jb),&! clear-sky shortwave downward flux at radiation time [W/m2]
        & rsucs_rt   = field%rsucs_rt         (:,:,jb),&! clear-sky shortwave upward   flux at radiation time [W/m2]
        !
        & rld_rt     = field%rld_rt           (:,:,jb),&! all-sky   longwave  downward flux at radiation time [W/m2]
        & rlu_rt     = field%rlu_rt           (:,:,jb),&! all-sky   longwave  upward   flux at radiation time [W/m2]
        !
        & rldcs_rt   = field%rldcs_rt         (:,:,jb),&! clear-sky longwave  downward flux at radiation time [W/m2]
        & rlucs_rt   = field%rlucs_rt         (:,:,jb),&! clear-sky longwave  upward   flux at radiation time [W/m2]
        !
        & rvds_dir_rt= field%rvds_dir_rt        (:,jb),&!< out  all-sky downward direct visible radiation at surface
        & rpds_dir_rt= field%rpds_dir_rt        (:,jb),&!< out  all-sky downward direct PAR     radiation at surface
        & rnds_dir_rt= field%rnds_dir_rt        (:,jb),&!< out  all-sky downward direct near-IR radiation at surface
        & rvds_dif_rt= field%rvds_dif_rt        (:,jb),&!< out  all-sky downward diffuse visible radiation at surface
        & rpds_dif_rt= field%rpds_dif_rt        (:,jb),&!< out  all-sky downward diffuse PAR     radiation at surface
        & rnds_dif_rt= field%rnds_dif_rt        (:,jb),&!< out  all-sky downward diffuse near-IR radiation at surface
        & rvus_rt    = field%rvus_rt            (:,jb),&!< out  all-sky upward visible radiation at surface
        & rpus_rt    = field%rpus_rt            (:,jb),&!< out  all-sky upward PAR     radiation at surfac
        & rnus_rt    = field%rnus_rt            (:,jb),&!< out  all-sky upward near-IR radiation at surface
        !
        ! output
        ! ------
        !
        & rsdt       = field%rsdt               (:,jb),&! all-sky   shortwave downward flux at current   time [W/m2]
        & rsut       = field%rsut               (:,jb),&! all-sky   shortwave upward   flux at current   time [W/m2]
        & rsds       = field%rsds               (:,jb),&! all-sky   shortwave downward flux at current   time [W/m2]
        & rsus       = field%rsus               (:,jb),&! all-sky   shortwave upward   flux at current   time [W/m2]
        !
        & rsutcs     = field%rsutcs             (:,jb),&! clear-sky shortwave upward   flux at current   time [W/m2]
        & rsdscs     = field%rsdscs             (:,jb),&! clear-sky shortwave downward flux at current   time [W/m2]
        & rsuscs     = field%rsuscs             (:,jb),&! clear-sky shortwave upward   flux at current   time [W/m2]
        !
        & rvds_dir   = field%rvds_dir           (:,jb),&!< out  all-sky downward direct visible radiation at surface
        & rpds_dir   = field%rpds_dir           (:,jb),&!< out  all-sky downward direct PAR     radiation at surface
        & rnds_dir   = field%rnds_dir           (:,jb),&!< out  all-sky downward direct near-IR radiation at surface
        & rvds_dif   = field%rvds_dif           (:,jb),&!< out  all-sky downward diffuse visible radiation at surface
        & rpds_dif   = field%rpds_dif           (:,jb),&!< out  all-sky downward diffuse PAR     radiation at surface
        & rnds_dif   = field%rnds_dif           (:,jb),&!< out  all-sky downward diffuse near-IR radiation at surface
        & rvus       = field%rvus               (:,jb),&!< out  all-sky upward visible radiation at surface
        & rpus       = field%rpus               (:,jb),&!< out  all-sky upward PAR     radiation at surfac
        & rnus       = field%rnus               (:,jb),&!< out  all-sky upward near-IR radiation at surface
        !
        & rlut       = field%rlut               (:,jb),&! all-sky   longwave  upward   flux at current   time [W/m2]
        & rlds       = field%rlds               (:,jb),&! all-sky   longwave  downward flux at current   time [W/m2]
        & rlus       = field%rlus               (:,jb),&! all-sky   longwave  upward   flux at current   time [W/m2]
        !
        & rlutcs     = field%rlutcs             (:,jb),&! clear-sky longwave  upward   flux at current   time [W/m2]
        & rldscs     = field%rldscs             (:,jb),&! clear-sky longwave  downward flux at current   time [W/m2]
        !
        & q_rsw      = zq_rsw                   (:,:) ,&! rad. heating by SW           [W/m2]
        & q_rlw      = zq_rlw                   (:,:) ) ! rad. heating by LW           [W/m2]

      IF (ltimer) CALL timer_stop(timer_radheat)

      ! heating accumulated
      zq_phy(jcs:jce,:) = zq_phy(jcs:jce,:) + zq_rsw(jcs:jce,:) + zq_rlw(jcs:jce,:)
      
      ! tendencies
      tend%ta_rsw(jcs:jce,:,jb) = zq_rsw(jcs:jce,:) * zconv(jcs:jce,:)
      tend%ta_rlw(jcs:jce,:,jb) = zq_rlw(jcs:jce,:) * zconv(jcs:jce,:)

      ! tendencies accumulated
      tend% ta(jcs:jce,:,jb) = tend% ta     (jcs:jce,:,jb) &
        &                    + tend% ta_rsw (jcs:jce,:,jb) &
        &                    + tend% ta_rlw (jcs:jce,:,jb)

    END IF ! lrad

    !-------------------------------------------------------------------
    ! 5. BOUNDARY LAYER AND SURFACE PROCESSES
    !-------------------------------------------------------------------
    ! Note: In ECHAM this part is located between "CALL radiation" and
    !       "CALL radheat".
    !
    ! 5.1 Emission of aerosols or other tracers. Not implemented yet.

      IF (ntrac>0) THEN
        !CALL tracer_emission()
        zxt_emis(jcs:jce,:) = 0._wp
      ENDIF
    !
    ! 5.2 Dry deposition of aerosols or other tracers. Not implemented yet.
    ! CALL dry_deposition()
    !

    ! 5.3 Turbulent mixing, part I:
    !     computation of exchange coefficients in the atmosphere and at the surface;
    !     build up the tridiagonal linear algebraic system;
    !     downward sweep (Gaussian elimination from top till level nlev-1)

    IF (echam_phy_config%lvdiff) THEN
      IF (ltimer) CALL timer_start(timer_vdiff_down)


      CALL vdiff_down( vdiff_config%lsfc_mom_flux,      &! in
                     & vdiff_config%lsfc_heat_flux,     &! in
                     & jce, nbdim, nlev, nlevm1, nlevp1,&! in
                     & ntrac, nsfc_type,                &! in
                     & iwtr, iice, ilnd,                &! in, indices of different surface types
                     & pdtime,                          &! in, time step
                     & field%coriol(:,jb),              &! in, Coriolis parameter
                     & zfrc(:,:),                       &! in, area fraction of each sfc type
                     & field% ts_tile(:,jb,:),          &! in, surface temperature
                     & field% ocu (:,jb),               &! in, ocean sfc velocity, u-component
                     & field% ocv (:,jb),               &! in, ocean sfc velocity, v-component
                     & field% presi_old(:,nlevp1,jb),   &! in, sfc pressure
                     & field%   ua(:,:,jb),             &! in, um1
                     & field%   va(:,:,jb),             &! in, vm1
                     & field%   ta(:,:,jb),             &! in, tm1
                     & field% qtrc(:,:,jb,iqv),         &! in, qm1
                     & field% qtrc(:,:,jb,iqc),         &! in, xlm1
                     & field% qtrc(:,:,jb,iqi),         &! in, xim1
                     & field%   qx(:,:,jb),             &! in, xlm1 + xim1
                     & field% qtrc(:,:,jb,iqt:),        &! in, xtm1
                     & field% presi_old(:,:,jb),        &! in, aphm1
                     & field% presm_old(:,:,jb),        &! in, apm1
                     & zdelp(:,:),                      &! in, layer thickness [Pa]
                     & field% geom(:,:,jb),             &! in, pgeom1 = geopotential above ground
                     & field% geoi(:,:,jb),             &! in, pgeohm1 = half-level geopotential
                     & field%   tv(:,:,jb),             &! in, virtual temperaturea
                     & field% aclc(:,:,jb),             &! in, cloud fraction
                     & zxt_emis,                        &! in, zxtems
                     & field% thvvar(:,:,jb),           &! in, variance of theta_v at step t-dt
                     & field%   xvar(:,:,jb),           &! in
                     & field% z0m_tile(:,jb,:),         &! in
                     & field%  tkem1(:,:,jb),           &! in, TKE at step t-dt
                     & field%  ustar(:,  jb),           &! inout
                     & field%  wstar(:,  jb),           &! out, convective velocity scale
                     & field%  wstar_tile(:,jb,:),      &! inout, convective velocity scale (each sfc type)
                     & field% qs_sfc_tile(:,jb,:),      &! out, sfc specific humidity at saturation
                     & ihpbl(:),                        &! out, for "vdiff_up"
                     & field%    ghpbl(:,jb),           &! out, for output
                     & field%      ri (:,:,jb),         &! out, for output
                     & zri_tile (:,:),                  &! out, for nsurf_diag
                     & field%  mixlen (:,:,jb),         &! out, for output
                     & field% cfm     (:,:,jb),         &! out, for output
                     & field% cfm_tile(:,jb,:),         &! out, for output and "vdiff_up"
                     & field% cfh     (:,:,jb),         &! out, for output
                     & field% cfh_tile(:,jb,:),         &! out, for output and "vdiff_up"
                     & field% cfv     (:,:,jb),         &! out, for output
                     & field% cftke   (:,:,jb),         &! out, for output
                     & field% cfthv   (:,:,jb),         &! out, for output
                     & zaa, zaa_btm, zbb, zbb_btm,      &! out, for "vdiff_up"
                     & zfactor_sfc(:),                  &! out, for "vdiff_up"
                     & zcpt_sfc_tile(:,:),              &! out, for "vdiff_up"
                     & zcptgz(:,:), zrhoh(:,:),         &! out, for "vdiff_up"
                     & zqshear(:,:),                    &! out, for "vdiff_up"
                     & zthvvar(:,:),                    &! out, for "vdiff_up"
                     & field%   thvsig(:,  jb),         &! out, for "cucall"
                     & ztkevn (:,:),                    &! out, for "vdiff_up"
                     & zch_tile(:,:),                   &! out, for "nsurf_diag"
!                     & zchn_tile(:,:),                  &! out, for "nsurf_diag"
!                     & zcdn_tile(:,:),                  &! out, for "nsurf_diag"
!                     & zcfnc_tile(:,:),                 &! out, for "nsurf_diag"
                     & zbn_tile(:,:),                   &! out, for "nsurf_diag"
                     & zbhn_tile(:,:),                  &! out, for "nsurf_diag"
                     & zbm_tile(:,:),                   &! out, for "nsurf_diag"
                     & zbh_tile(:,:),                   &! out, for "nsurf_diag"
                     & pcsat = field% csat(:,jb),       &! in, optional, area fraction with wet land surface
                     & pcair = field% cair(:,jb),       &! in, optional, area fraction with wet land surface (air)
                     & paz0lh = field% z0h_lnd(:,jb))     ! in, optional, roughness length for heat over land

      IF (ltimer) CALL timer_stop(timer_vdiff_down)

    ! 5.4 Surface processes that provide time-dependent lower boundary
    !     condition for wind, temperature, tracer concentraion, etc.
    !KF To avoid
        field% lhflx_tile(:,jb,:) = 0._wp
        field% shflx_tile(:,jb,:) = 0._wp
        field% evap_tile (:,jb,:) = 0._wp

        IF (ltimer) CALL timer_start(timer_surface)

        CALL update_surface( vdiff_config%lsfc_heat_flux,  &! in
          & vdiff_config%lsfc_mom_flux,   &! in
          & pdtime,                       &! in, time step
          & jg, jce, nbdim, field%kice,   &! in
          & nlev, nsfc_type,              &! in
          & iwtr, iice, ilnd,             &! in, indices of surface types
          & zfrc(:,:),                    &! in, area fraction
          & field% cfh_tile(:,jb,:),      &! in, from "vdiff_down"
          & field% cfm_tile(:,jb,:),      &! in, from "vdiff_down"
          & zfactor_sfc(:),               &! in, from "vdiff_down"
          & field% ocu (:,jb),            &! in, ocean sfc velocity, u-component
          & field% ocv (:,jb),            &! in, ocean sfc velocity, v-component
          & zaa, zaa_btm, zbb, zbb_btm,   &! inout
          & zcpt_sfc_tile(:,:),           &! inout, from "vdiff_down", for "vdiff_up"
          & field%qs_sfc_tile(:,jb,:),    &! inout, from "vdiff_down", for "vdiff_up"
          & field% ts_tile(:,jb,:),       &! inout
          & field%u_stress    (:,  jb),   &! out
          & field%v_stress    (:,  jb),   &! out
          & field% lhflx      (:,  jb),   &! out
          & field% shflx      (:,  jb),   &! out
          & field%  evap      (:,  jb),   &! out, for "cucall"
          & field%u_stress_tile  (:,jb,:),   &! out
          & field%v_stress_tile  (:,jb,:),   &! out
          & field% lhflx_tile    (:,jb,:),   &! out
          & field% shflx_tile    (:,jb,:),   &! out
          & field% dshflx_dT_tile(:,jb,:),   &! out for Sea ice
          & field%  evap_tile    (:,jb,:),   &! out
                                !! optional
          & nblock = jb,                  &! in
          & lsm = field%lsmask(:,jb), &!< in, land-sea mask
          & pu    = field% ua(:,nlev,jb), &! in, um1
          & pv    = field% va(:,nlev,jb), &! in, vm1
          & ptemp = field% ta(:,nlev,jb), &! in, tm1
          & pq = field% qtrc(:,nlev,jb,iqv),  &! in, qm1
          & prsfl = field% rsfl(:,jb),    &! in, rain surface large scale (from cloud)
          & prsfc = field% rsfc(:,jb),    &! in, rain surface concective (from cucall)
          & pssfl = field% ssfl(:,jb),    &! in, snow surface large scale (from cloud)
          & pssfc = field% ssfc(:,jb),    &! in, snow surface concective (from cucall)
          & rlds        = field% rlds (:,jb), &! in,  downward surface  longwave flux [W/m2]
          & rlus        = field% rlus (:,jb), &! inout, upward surface  longwave flux [W/m2]
          & rsds        = field% rsds (:,jb), &! in,  downward surface shortwave flux [W/m2]
          & rsus        = field% rsus (:,jb), &! inout, upward surface shortwave flux [W/m2]
          !
          & rvds_dir   = field%rvds_dir   (:,jb), &! in, all-sky downward direct visible radiation at surface
          & rpds_dir   = field%rpds_dir   (:,jb), &! in, all-sky downward direct PAR     radiation at surface
          & rnds_dir   = field%rnds_dir   (:,jb), &! in, all-sky downward direct near-IR radiation at surface
          & rvds_dif   = field%rvds_dif   (:,jb), &! in, all-sky downward diffuse visible radiation at surface
          & rpds_dif   = field%rpds_dif   (:,jb), &! in, all-sky downward diffuse PAR     radiation at surface
          & rnds_dif   = field%rnds_dif   (:,jb), &! in, all-sky downward diffuse near-IR radiation at surface
          !
          & presi_old = field% presi_old(:,:,jb),&! in, paphm1, half level pressure
          & pcosmu0 = field% cosmu0(:,jb),&! in, amu0_x, cos of zenith angle
          & pch_tile = zch_tile(:,:),     &! in, from "vdiff_down" for JSBACH
          & pcsat = field%csat(:,jb),      &! inout, area fraction with wet land surface
          & pcair = field%cair(:,jb),      &! inout, area fraction with wet land surface (air)
          & tte_corr = ztte_corr(:),       &! out, tte correction for snow melt over land
          & z0m_tile = field% z0m_tile(:,jb,:), &! inout, roughness length for momentum over tiles
          & z0h_lnd  = field% z0h_lnd (:,jb),   &! out, roughness length for heat over land
          & albvisdir      = field% albvisdir     (:,jb)  ,                    &! inout
          & albnirdir      = field% albnirdir     (:,jb)  ,                    &! inout
          & albvisdif      = field% albvisdif     (:,jb)  ,                    &! inout
          & albnirdif      = field% albnirdif     (:,jb)  ,                    &! inout
          & albvisdir_tile = field% albvisdir_tile(:,jb,:),                    &! inout
          & albnirdir_tile = field% albnirdir_tile(:,jb,:),                    &! inout
          & albvisdif_tile = field% albvisdif_tile(:,jb,:),                    &! inout
          & albnirdif_tile = field% albnirdif_tile(:,jb,:),                    &! inout
          & albedo         = field% albedo        (:,jb)  ,                    &! inout
          & albedo_tile    = field% albedo_tile(:,jb,:),                       &! inout
          & ptsfc     = field%ts    (:,jb),                        &! out
          & ptsfc_rad = field%ts_rad(:,jb),                        &! out
          & rlns_tile = field%lwflxsfc_tile(:,jb,:),               &! out (for coupling)
          & rsns_tile = field%swflxsfc_tile(:,jb,:),               &! out (for coupling)
          & Tsurf = field% Tsurf(:,:,jb),  &! inout, for sea ice
          & T1    = field% T1   (:,:,jb),  &! inout, for sea ice
          & T2    = field% T2   (:,:,jb),  &! inout, for sea ice
          & hi    = field% hi   (:,:,jb),  &! in, for sea ice
          & hs    = field% hs   (:,:,jb),  &! in, for sea ice
          & conc  = field% conc (:,:,jb),  &! in, for sea ice
          & Qtop  = field% Qtop (:,:,jb),  &! out, for sea ice
          & Qbot  = field% Qbot (:,:,jb),  &! out, for sea ice
          & albvisdir_ice = field% albvisdir_ice(:,:,jb), &! inout ice albedos
          & albnirdir_ice = field% albnirdir_ice(:,:,jb), &! inout
          & albvisdif_ice = field% albvisdif_ice(:,:,jb), &! inout
          & albnirdif_ice = field% albnirdif_ice(:,:,jb))  ! inout

        IF (ltimer) CALL timer_stop(timer_surface)

    ! 5.5 Turbulent mixing, part II:
    !     - Elimination for the lowest model level using boundary conditions
    !       provided by the surface model(s);
    !     - Back substitution to get solution of the tridiagonal system;
    !     - Compute tendencies and additional diagnostics.

      IF (ltimer) CALL timer_start(timer_vdiff_up)

      CALL vdiff_up( jce, nbdim, nlev, nlevm1, nlevp1,&! in
                   & ntrac, nsfc_type,                &! in
                   & iwtr,                            &! in, indices of different sfc types
                   & pdtime,                          &! in, time steps
                   & zfrc(:,:),                       &! in, area fraction of each sfc type
                   & field% cfm_tile(:,jb,:),         &! in
                   & zaa,                             &! in, from "vdiff_down"
                   &   ihpbl(:),                      &! in, from "vdiff_down"
                   &  zcptgz(:,:),                    &! in, from "vdiff_down"
                   &   zrhoh(:,:),                    &! in, from "vdiff_down"
                   & zqshear(:,:),                    &! in, from "vdiff_down"
                   & field%   ua(:,:,jb),             &! in, um1
                   & field%   va(:,:,jb),             &! in, vm1
                   & field%   ta(:,:,jb),             &! in, tm1
                   & field% qtrc(:,:,jb,iqv),         &! in, qm1
                   & field% qtrc(:,:,jb,iqc),         &! in, xlm1
                   & field% qtrc(:,:,jb,iqi),         &! in, xim1
                   & field% qtrc(:,:,jb,iqt:),        &! in, xtm1
                   & zdelp(:,:),                      &! in, layer thickness [Pa]
                   & field% geom(:,:,jb),             &! in, pgeom1 = geopotential above ground
                   &      ztkevn(:,:),                &! in, tke at intermediate time step
                   & field%tkem1(:,:,jb),             &! in, TKE at step t-dt
                   & ztte_corr(:),                    &! in
                   & zbb,                             &! inout
                   & zthvvar(:,:),                    &! inout
                   & field%   xvar(:,:,jb),           &! inout
                   & field% z0m_tile(:,jb,:),         &! inout
                   & field% kedisp(:,  jb),           &! inout, "vdis" in ECHAM
                   &  tend%   ua_vdf(:,:,jb),         &! out
                   &  tend%   va_vdf(:,:,jb),         &! out
                   &  tend%   ta_vdf(:,:,jb),         &! out
!!$                   &          zq_vdf(:,:),            &! out   heating W/m2
                   &  tend% qtrc_vdf(:,:,jb,iqv),     &! out
                   &  tend% qtrc_vdf(:,:,jb,iqc),     &! out
                   &  tend% qtrc_vdf(:,:,jb,iqi),     &! out
                   &  tend% qtrc_vdf(:,:,jb,iqt:),    &! out
                   &  zqtvar_prod,                    &! out, for "cloud" ("zvdiffp" in echam)
                   &  zvmixtau,                       &! out, for "cloud"
                   & field%   z0m   (:,  jb),         &! out, for the next step
                   & field%   thvvar(:,:,jb),         &! out, for the next step
                   & field%   thvsig(:,  jb),         &! out, for "cucall"
                   & field%      tke(:,:,jb),         &! out
                   & field%   sh_vdiff(:,  jb),       &! out, for energy diagnostic
                   & field%   qv_vdiff(:,  jb)        )! out, for energy diagnostic

      IF (ltimer) CALL timer_stop(timer_vdiff_up)

!!$      ! heating accumulated
!!$      zq_phy(jcs:jce,:) = zq_phy(jcs:jce,:) + zq_vdf(jcs:jce,:)
!!$
!!$      ! tendency
!!$      tend% temp_vdf(jcs:jce,:,jb) = zq_vdf(jcs:jce,:)*zconv(jcs:jce,:)
!!$
      ! tendencies accumulated
      tend%   ua(jcs:jce,:,jb)      = tend%   ua(jcs:jce,:,jb)      + tend%   ua_vdf(jcs:jce,:,jb)
      tend%   va(jcs:jce,:,jb)      = tend%   va(jcs:jce,:,jb)      + tend%   va_vdf(jcs:jce,:,jb)
      tend%   ta(jcs:jce,:,jb)      = tend%   ta(jcs:jce,:,jb)      + tend%   ta_vdf(jcs:jce,:,jb)
      tend% qtrc(jcs:jce,:,jb,iqv)  = tend% qtrc(jcs:jce,:,jb,iqv)  + tend% qtrc_vdf(jcs:jce,:,jb,iqv)
      tend% qtrc(jcs:jce,:,jb,iqc)  = tend% qtrc(jcs:jce,:,jb,iqc)  + tend% qtrc_vdf(jcs:jce,:,jb,iqc)
      tend% qtrc(jcs:jce,:,jb,iqi)  = tend% qtrc(jcs:jce,:,jb,iqi)  + tend% qtrc_vdf(jcs:jce,:,jb,iqi)
      tend% qtrc(jcs:jce,:,jb,iqt:) = tend% qtrc(jcs:jce,:,jb,iqt:) + tend% qtrc_vdf(jcs:jce,:,jb,iqt:)

!    ! TIME FILTER FOR TURBULENT KINETIC ENERGY
!
!    IF(.NOT.lstart) THEN
!      zeps=eps
!    ELSE
!      zeps=0._wp
!    END IF
!    DO 397 jk=ktdia,klev
!      DO 396 jl=1,kproma
!        ptkem1(jl,jk)=ptkem(jl,jk)                                    &
!                  +zeps*(ptkem1(jl,jk)-2._wp*ptkem(jl,jk)+ptke(jl,jk))
!        ptkem(jl,jk)=ptke(jl,jk)
!396   END DO
!397 END DO

      ! 2-tl-scheme
      field% tkem1(jcs:jce,:,jb) = field% tke  (jcs:jce,:,jb)

    ! 5.6 Turbulent mixing, part III:
    !     - Further diagnostics.

    CALL nsurf_diag( jce, nbdim, nsfc_type,           &! in
                   & ilnd,                            &! in
                   & zfrc(:,:),                       &! in
                   & field%  qtrc(:,nlev,jb,iqv),     &! in humidity qm1
                   & field%    ta(:,nlev,jb),         &! in tm1
                   & field% presm_old(:,nlev,jb),     &! in, apm1
                   & field% presi_old(:,nlevp1,jb),   &! in, aphm1
                   & field%   qx(:,nlev,jb),          &! in, xlm1 + xim1
                   & field%   ua(:,nlev,jb),          &! in, um1
                   & field%   va(:,nlev,jb),          &! in, vm1
                   & field% ocu (:,jb),               &! in, ocean sfc velocity, u-component
                   & field% ocv (:,jb),               &! in, ocean sfc velocity, v-component
                   & field%  geom(:,nlev,jb),         &! in geopotential above surface
                   & zcptgz(:,nlev),                  &! in dry static energy
                   & zcpt_sfc_tile(:,:),              &! in dry static energy
                   & zbn_tile(:,:),                   &! in for diagnostic
                   & zbhn_tile(:,:),                  &! in for diagnostic
                   & zbh_tile(:,:),                   &! in for diagnostic
                   & zbm_tile(:,:),                   &! in for diagnostic
                   & zri_tile(:,:),                   &! in 
                   & field%sfcWind(:,  jb),           &! out 10m windspeed
                   & field%    tas(:,  jb),           &! out temperature in 2m
                   & field%   dew2(:,  jb),           &! out dew point temperature in 2m
                   & field%    uas(:,  jb),           &! out zonal wind in 10m
                   & field%    vas(:,  jb),           &! out meridional wind in 10m
                   & field%tasmax (:,  jb),           &! out max 2m temperature
                   & field%tasmin (:,  jb),           &! out min 2m temperature
                   & field%sfcWind_tile(:,jb,:),      &! out 10m windspeed on tiles
                   & field%    tas_tile(:,jb,:),      &! out temperature in 2m on tiles
                   & field%   dew2_tile(:,jb,:),      &! out dew point temperature in 2m on tiles
                   & field%    uas_tile(:,jb,:),      &! out zonal wind in 10m on tiles
                   & field%    vas_tile(:,jb,:)       )! out meridional wind in 10m on tiles

    ELSE
      zvmixtau   (jcs:jce,:) = 0._wp
      zqtvar_prod(jcs:jce,:) = 0._wp

    ENDIF !lvdiff

    IF (echam_phy_config%lrad) THEN

      ! Heating due to the fact that surface model only used part of longwave radiation to compute new surface temperature
      zq_rlw_impl(jcs:jce) =                                                  &
        &   ( (field%rld_rt(jcs:jce,nlev,jb)-field%rlu_rt(jcs:jce,nlev,jb))   & ! (( rlns from "radiation"
        &    -(field%rlds  (jcs:jce,jb)     -field%rlus  (jcs:jce,jb)     ))  & !   -rlns from "radheating" and "update_surface")
        &  -zq_rlw(jcs:jce,nlev)                                                ! old heating from radheat

      ! Heating accumulated
      zq_phy(jcs:jce,nlev) = zq_phy(jcs:jce,nlev) + zq_rlw_impl(jcs:jce)

      ! Tendency
      tend%ta_rlw_impl(jcs:jce,jb) = zq_rlw_impl(jcs:jce) * zconv(jcs:jce,nlev)

      ! Tendencies accumulated
      tend%ta(jcs:jce,nlev,jb) = tend%ta(jcs:jce,nlev,jb) + tend%ta_rlw_impl(jcs:jce,jb)

    END IF

    IF (echam_phy_config%lcariolle) THEN
      avi%tmprt(jcs:jce,:)=field%ta(jcs:jce,:,jb)
      avi%vmr2molm2(jcs:jce,:)=field%mdry(jcs:jce,:,jb)/amd*1.e3_wp
      avi%pres(jcs:jce,:)=field%presm_old(jcs:jce,:,jb)
      avi%cell_center_lat(jcs:jce)=field%clat(jcs:jce,jb)
      avi%lday(jcs:jce)=field%cosmu0(jcs:jce,jb)>1.e-3_wp
      avi%ldown=.TRUE.
      current_time_interpolation_weights = calculate_time_interpolation_weights(this_datetime)
      time_interpolation%imonth1=current_time_interpolation_weights%month1_index
      time_interpolation%imonth2=current_time_interpolation_weights%month2_index
      time_interpolation%weight1=current_time_interpolation_weights%weight1
      time_interpolation%weight2=current_time_interpolation_weights%weight2
      avi%o3_vmr(jcs:jce,:)=field%qtrc(jcs:jce,:,jb,io3)*amd/amo3
      CALL lcariolle_do3dt(                                                    &
         & jcs,                    jce,                nbdim,                  &
         & nlev,                   time_interpolation, lcariolle_lat_intp_li,  &
         & lcariolle_pres_intp_li, avi,                do3dt                   )
      tend% qtrc(jcs:jce,:,jb,io3) = tend% qtrc(jcs:jce,:,jb,io3) + do3dt(jcs:jce,:)*amo3/amd
    END IF

    !-------------------------------------------------------------------
    ! 6. ATMOSPHERIC GRAVITY WAVES
    !-------------------------------------------------------------------

    ! 6.1   CALL SUBROUTINE GW_HINES

    IF (echam_phy_config%lgw_hines) THEN

      zlat_deg(jcs:jce) = field% clat(jcs:jce,jb) * 180._wp/pi

      IF (ltimer) call timer_start(timer_gw_hines)

      CALL gw_hines ( jg                       ,&
        &             nbdim                    ,&
        &             jcs                      ,&
        &             jce                      ,&
        &             nc                       ,&
        &             nlev                     ,&
        &             field% presi_old(:,:,jb) ,&
        &             field% presm_old(:,:,jb) ,&
        &             field%   ta(:,:,jb)      ,&
        &             field%   ua(:,:,jb)      ,&
        &             field%   va(:,:,jb)      ,&
        &             zlat_deg(:)              ,&
!!$        &             aprflux(:,krow)          ,&
        &             zdis_gwh(:,:)            ,&
        &             tend%   ua_gwh(:,:,jb)   ,&
        &             tend%   va_gwh(:,:,jb) )

      IF (ltimer) call timer_stop(timer_gw_hines)

      ! heating
      zq_gwh(jcs:jce,:) = zdis_gwh(jcs:jce,:) * field%mair(jcs:jce,:,jb)

      ! heating accumulated
      zq_phy(jcs:jce,:) = zq_phy(jcs:jce,:) + zq_gwh(jcs:jce,:)

      ! tendency
      tend% ta_gwh(jcs:jce,:,jb) = zq_gwh(jcs:jce,:)*zconv(jcs:jce,:)

      ! tendencies accumulated
      tend%   ta(jcs:jce,:,jb) = tend%   ta(jcs:jce,:,jb) + tend%   ta_gwh(jcs:jce,:,jb)
      tend%   ua(jcs:jce,:,jb) = tend%   ua(jcs:jce,:,jb) + tend%   ua_gwh(jcs:jce,:,jb)
      tend%   va(jcs:jce,:,jb) = tend%   va(jcs:jce,:,jb) + tend%   va_gwh(jcs:jce,:,jb)

    END IF !lgw_hines


    ! 6.2   CALL SUBROUTINE SSODRAG

    IF (echam_phy_config%lssodrag) THEN

      IF (ltimer) call timer_start(timer_ssodrag)

       CALL ssodrag( nc                                        ,& ! in,  number of cells/columns in loop (jce-jcs+1)
                     nbdim                                     ,& ! in,  dimension of block of cells/columns
                     nlev                                      ,& ! in,  number of levels
                     !
                     field% clat(:,jb)                         ,& ! in,  Latitude in radians
                     pdtime                                    ,& ! in,  time step length
                     !
                     field% presi_old(:,:,jb)                  ,& ! in,  p at half levels
                     field% presm_old(:,:,jb)                  ,& ! in,  p at full levels
                     field% geom(:,:,jb)                       ,& ! in,  geopotential above surface (t-dt)
                     field%   ta(:,:,jb)                       ,& ! in,  T
                     field%   ua(:,:,jb)                       ,& ! in,  u
                     field%   va(:,:,jb)                       ,& ! in,  v
                     !
                     field% oromea(:,jb)                       ,& ! in,  Mean Orography (m)
                     field% orostd(:,jb)                       ,& ! in,  SSO standard deviation (m)
                     field% orosig(:,jb)                       ,& ! in,  SSO slope
                     field% orogam(:,jb)                       ,& ! in,  SSO Anisotropy
                     field% orothe(:,jb)                       ,& ! in,  SSO Angle
                     field% oropic(:,jb)                       ,& ! in,  SSO Peaks elevation (m)
                     field% oroval(:,jb)                       ,& ! in,  SSO Valleys elevation (m)
                     !
                     field% u_stress_sso(:,jb)                 ,& ! out, u-gravity wave stress
                     field% v_stress_sso(:,jb)                 ,& ! out, v-gravity wave stress
                     field% dissipation_sso(:,jb)              ,& ! out, dissipation by gravity wave drag
                     !
                     zdis_sso(:,:)                             ,& ! out, energy dissipation rate
                     tend%   ua_sso(:,:,jb)                    ,& ! out, tendency of zonal wind
                     tend%   va_sso(:,:,jb)                     ) ! out, tendency of meridional wind

      IF (ltimer) call timer_stop(timer_ssodrag)

      ! heating
      zq_sso(jcs:jce,:) = zdis_sso(jcs:jce,:) * field%mair(jcs:jce,:,jb)

      ! heating accumulated
      zq_phy(jcs:jce,:) = zq_phy(jcs:jce,:) + zq_sso(jcs:jce,:)

      ! tendency
      tend% ta_sso(jcs:jce,:,jb) = zq_sso(jcs:jce,:)*zconv(jcs:jce,:)

      ! tendencies accumulated
      tend%   ta(jcs:jce,:,jb) = tend%   ta(jcs:jce,:,jb) + tend%   ta_sso(jcs:jce,:,jb)
      tend%   ua(jcs:jce,:,jb) = tend%   ua(jcs:jce,:,jb) + tend%   ua_sso(jcs:jce,:,jb)
      tend%   va(jcs:jce,:,jb) = tend%   va(jcs:jce,:,jb) + tend%   va_sso(jcs:jce,:,jb)

    END IF ! SSODRAG

    !-------------------------------------------------------------------

    ! Update physics state for input to next parameterization
    zta  (:,:)   =     field% ta  (:,:,jb)   + pdtime*tend% ta  (:,:,jb)
    zqtrc(:,:,:) = MAX(field% qtrc(:,:,jb,:) + pdtime*tend% qtrc(:,:,jb,:), 0.0_wp)
    zua  (:,:)   =     field% ua  (:,:,jb)   + pdtime*tend% ua  (:,:,jb)
    zva  (:,:)   =     field% va  (:,:,jb)   + pdtime*tend% va  (:,:,jb)

    !-------------------------------------------------------------------
    ! 7. CONVECTION PARAMETERISATION
    !-------------------------------------------------------------------

    zqtrc_cnd(:,:) = zqtrc(:,:,iqc) + zqtrc(:,:,iqi)
    ztend_qv(:,:)  = tend%qtrc_dyn(:,:,jb,iqv) + tend%qtrc_phy(:,:,jb,iqv)

    IF (echam_phy_config%lconv) THEN

      IF (ltimer) CALL timer_start(timer_convection)

      CALL cumastr(jce, nbdim,                   &! in
        &          nlev, nlevp1, nlevm1,         &! in
        &          pdtime,                       &! in
        &          field% zf       (:,:,jb),     &! in
        &          field% zh       (:,:,jb),     &! in
        &          field% mdry     (:,:,jb),     &! in
        &                zta       (:,:),        &! in
        &                zqtrc     (:,:,   iqv), &! in
        &                zqtrc_cnd (:,:),        &! in
        &                zua       (:,:),        &! in
        &                zva       (:,:),        &! in
        &          ntrac,                        &! in
        &          field% lfland   (:,  jb),     &! in
        &                zqtrc     (:,:,   iqt:),&! in
        &          field% omega    (:,:,jb),     &! in
        &          field% evap     (:,  jb),     &! in
        &          field% presm_new(:,:,jb),     &! in
        &          field% presi_new(:,:,jb),     &! in
        &          field% geom     (:,:,jb),     &! in
        &          field% geoi     (:,:,jb),     &! in
        &                ztend_qv  (:,:),        &! in
        &          field% thvsig   (:,  jb),     &! in
        &          echam_conv_config%cevapcu,    &! in
        &          itype,                        &! out
        &          ictop,                        &! out
        &          field% rsfc     (:,  jb),     &! out
        &          field% ssfc     (:,  jb),     &! out
        &          field% con_dtrl (:,jb),       &! out
        &          field% con_dtri (:,jb),       &! out
        &          field% con_iteqv(:,jb),       &! out
        &                   zq_cnv (:,:),        &! out
        &           tend%   ua_cnv (:,:,jb),     &! out
        &           tend%   va_cnv (:,:,jb),     &! out
        &           tend% qtrc_cnv (:,:,jb,iqv), &! out
        &           tend% qtrc_cnv (:,:,jb,iqt:),&! out
        &           tend% qtrc_cnv (:,:,jb,iqc), &! out
        &           tend% qtrc_cnv (:,:,jb,iqi), &! out
        &                ztop      (:)           )! out

      IF (ltimer) CALL timer_stop(timer_convection)

      ! store convection type as real value
      field% rtype(:,jb) = REAL(itype(:),wp)

      ! keep minimum conv. cloud top pressure (= max. conv. cloud top height) of this output interval
      field% topmax(:,jb) = MIN(field% topmax(:,jb),ztop(:))

      ! heating accumulated
      zq_phy(:,:) = zq_phy(:,:) + zq_cnv(:,:)

      ! tendency
      tend% ta_cnv(:,:,jb) = zq_cnv(:,:)*zconv(:,:)

      ! tendencies accumulated
      tend%   ua(:,:,jb)      = tend%   ua(:,:,jb)      + tend%   ua_cnv(:,:,jb)
      tend%   va(:,:,jb)      = tend%   va(:,:,jb)      + tend%   va_cnv(:,:,jb)
      tend%   ta(:,:,jb)      = tend%   ta(:,:,jb)      + tend%   ta_cnv(:,:,jb)
      tend% qtrc(:,:,jb,iqv)  = tend% qtrc(:,:,jb,iqv)  + tend% qtrc_cnv(:,:,jb,iqv)
      tend% qtrc(:,:,jb,iqc)  = tend% qtrc(:,:,jb,iqc)  + tend% qtrc_cnv(:,:,jb,iqc)
      tend% qtrc(:,:,jb,iqi)  = tend% qtrc(:,:,jb,iqi)  + tend% qtrc_cnv(:,:,jb,iqi)
      tend% qtrc(:,:,jb,iqt:) = tend% qtrc(:,:,jb,iqt:) + tend% qtrc_cnv(:,:,jb,iqt:)

    ELSE ! NECESSARY COMPUTATIONS IF MASSFLUX IS BY-PASSED

      ictop(:)   = nlev-1
      itype(:)   = 0

    ENDIF !lconv

    !-------------------------------------------------------------
    ! Update provisional physics state
    !
!!$    field% ta  (:,:,jb)     = field% ta  (:,:,jb)     + tend% ta  (:,:,jb)    *pdtime
!!$    field% qtrc(:,:,jb,iqv) = field% qtrc(:,:,jb,iqv) + tend% qtrc(:,:,jb,iqv)*pdtime
    field% qtrc(:,:,jb,iqc) = field% qtrc(:,:,jb,iqc) + tend% qtrc(:,:,jb,iqc)*pdtime
    field% qtrc(:,:,jb,iqi) = field% qtrc(:,:,jb,iqi) + tend% qtrc(:,:,jb,iqi)*pdtime
    !
    !-------------------------------------------------------------


    !-------------------------------------------------------------
    ! 7. LARGE SCALE CONDENSATION.
    !-------------------------------------------------------------
    IF(echam_phy_config%lcond) THEN

      !IF (lcotra) CALL get_col_pol( tend%ta(:,:,jb),tend%qtrc(:,:,jb,iqv),jb )

      IF (ltimer) CALL timer_start(timer_cloud)

      CALL cloud(jce, nbdim, jks, nlev,        &! in
        &        pdtime,                       &! in
        &        ictop,                        &! in (from "cucall")
        &        field% presm_old(:,:,jb),     &! in
        &        field% dz       (:,:,jb),     &! in
        &        field% mdry     (:,:,jb),     &! in
        &        field% rho      (:,:,jb),     &! in
        &               zcpair   (:,:),        &! in
        &        field% acdnc    (:,:,jb),     &! in  acdnc
        &        field% ta       (:,:,jb),     &! in  tm1
        &        field% qtrc     (:,:,jb,iqv), &! in  qm1
        &        field% qtrc     (:,:,jb,iqc), &! in  xlm1
        &        field% qtrc     (:,:,jb,iqi), &! in  xim1
        &         tend% ta       (:,:,jb),     &! in  tte
        &         tend% qtrc     (:,:,jb,iqv), &! in  qte
        !
        &        itype,                        &! inout
        &        field% aclc     (:,:,jb),     &! inout
        !
        &        field% aclcov   (:,  jb),     &! out
        &        field% rsfl     (:,  jb),     &! out
        &        field% ssfl     (:,  jb),     &! out
        &        field% relhum   (:,:,jb),     &! out
        &               zq_cld   (:,:),        &! out
        &         tend% qtrc_cld (:,:,jb,iqv), &! out
        &         tend% qtrc_cld (:,:,jb,iqc), &! out
        &         tend% qtrc_cld (:,:,jb,iqi)  )! out

      IF (ltimer) CALL timer_stop(timer_cloud)

      field% rtype(:,jb) = REAL(itype(:),wp)

      ! heating accumulated
      zq_phy(:,:) = zq_phy(:,:) + zq_cld(:,:)

      ! tendency
      tend% ta_cld(:,:,jb) = zq_cld(:,:)*zconv(:,:)

      ! tendencies accumulated
      tend%   ta(:,:,jb)      = tend%   ta(:,:,jb)      + tend%   ta_cld(:,:,jb)
      tend% qtrc(:,:,jb,iqv)  = tend% qtrc(:,:,jb,iqv)  + tend% qtrc_cld(:,:,jb,iqv)
      tend% qtrc(:,:,jb,iqc)  = tend% qtrc(:,:,jb,iqc)  + tend% qtrc_cld(:,:,jb,iqc)
      tend% qtrc(:,:,jb,iqi)  = tend% qtrc(:,:,jb,iqi)  + tend% qtrc_cld(:,:,jb,iqi)

    ENDIF !lcond

    ! KF accumulate fields for diagnostics

    !  total precipitation flux
       field% totprec (:,jb)     =  field% rsfl (:,jb) & ! rain large scale
            &                      +field% ssfl (:,jb) & ! snow large scale
            &                      +field% rsfc (:,jb) & ! rain convection
            &                      +field% ssfc (:,jb)   ! snow convection

    ! Now compute tendencies from physics alone

    tend%   ta_phy (:,:,jb)   = tend%   ta (:,:,jb)   - tend%   ta_phy (:,:,jb)
    tend%   ua_phy (:,:,jb)   = tend%   ua (:,:,jb)   - tend%   ua_phy (:,:,jb)
    tend%   va_phy (:,:,jb)   = tend%   va (:,:,jb)   - tend%   va_phy (:,:,jb)
    tend% qtrc_phy (:,:,jb,:) = tend% qtrc (:,:,jb,:) - tend% qtrc_phy (:,:,jb,:)

    ! And convert constant pressure temperture tendency to constant volume temperature tendency
    tend% ta_phy (:,:,jb) = tend% ta_phy(:,:,jb)*zcpair(:,:)/zcvair(:,:)

    ! Done. Disassociate pointers.
    NULLIFY(field,tend)

  END SUBROUTINE echam_phy_main
  !-------------

END MODULE mo_echam_phy_main
