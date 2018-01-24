#if (defined (__GNUC__) || defined(__SUNPRO_F95) || defined(__SX__))
#define HAVE_F95
#endif

!>
!! Description:  Contains the data structures
!!  to store the physical model state and other auxiliary variables
!!  in order to run the ECHAM physics.
!!  This module should be an analogon to 'mo_nonhydro_types.f90'

!!  TODO/To think about:
!     - should physics be called before or after dynamics?
!     - allocate fluxes at edges instead at the centers?
!     - horizontal/vertical tracer flux (reconstruct q'v_n' into q'u' and q'v') ?
!     - provide the "virt_inc" with meaning
!     - where to provide the lat/lon info for radiation?
!     - how to implement the echam-modules - rewriting them or "capsulate"?
!     - revision of fields if there are needed or tp be replaced
!     - fill the physics tendency construction/destruction subroutine
!     - later implement already calculated icon gradients for echam physics
!     - think about variables for flexible time steps
!!
!! @author Kristina Froehlich, DWD
!! @author Marco Giorgetta, MPI-M
!!
!!
!! @par Revision History
!! Initial  by Kristina Froehlich (2009-06-10)
!! Memory allocation method changed from explicit allocation to Luis' 
!! infrastructure by Kristina Froehlich (MPI-M, 2011-04-27)
!! Added clch, clcm, clcl, hbas_con, htop_con by Helmut Frank (DWD, 2013-01-17)
!! Added hzerocl and gusts                    by Helmut Frank (DWD, 2013-03-13)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_nwp_phy_types

  USE mo_kind,                ONLY: wp, vp2
  USE mo_fortran_tools,       ONLY: t_ptr_2d3d,t_ptr_tracer

  IMPLICIT NONE
  PRIVATE



  !public interface
  !
  !types
  PUBLIC :: t_nwp_phy_diag
  PUBLIC :: t_nwp_phy_tend



  !> derived data type for synthetic satellite images
  TYPE t_rttov_image
    REAL(wp), POINTER :: p(:,:)      ! pointer to 2D image
  END TYPE t_rttov_image


  !
  !!data structure defining model states
  !
  !!diagnostic variables
  !

  TYPE t_nwp_phy_diag

    TYPE(t_ptr_2d3d),ALLOCATABLE :: tot_ptr(:)  !< pointer array: one pointer for each tot var (grid+subgrid)
    TYPE(t_ptr_2d3d),ALLOCATABLE :: tci_ptr(:)  !< pointer array: total column-integrated values
    TYPE(t_ptr_2d3d),ALLOCATABLE :: tav_ptr(:)  !< pointer array: average of tci

    TYPE(t_ptr_2d3d),ALLOCATABLE :: cfm_ptr(:)  !< pointer array: average of cfm
    TYPE(t_ptr_2d3d),ALLOCATABLE :: cfh_ptr(:)  !< pointer array: average of cfh
    TYPE(t_ptr_2d3d),ALLOCATABLE :: z0m_ptr(:)  !< pointer array: average of z0m
    TYPE(t_ptr_2d3d),ALLOCATABLE :: albdif_t_ptr(:)   !< pointer array: tile-specific albedo (shortwave)
    TYPE(t_ptr_2d3d),ALLOCATABLE :: albvisdif_t_ptr(:)!< pointer array: tile-specific albedo (UV/visible)
    TYPE(t_ptr_2d3d),ALLOCATABLE :: albnirdif_t_ptr(:)!< pointer array: tile-specific albedo (NIR)
    TYPE(t_ptr_2d3d),ALLOCATABLE :: swflxsfc_t_ptr(:) !< pointer array: shortwave net flux at surface
    TYPE(t_ptr_2d3d),ALLOCATABLE :: lwflxsfc_t_ptr(:) !< pointer array: longwave net flux at surface
    TYPE(t_ptr_2d3d),ALLOCATABLE :: tcm_t_ptr(:) !< pointer array: turbulent transfer coefficients for momentum
    TYPE(t_ptr_2d3d),ALLOCATABLE :: tch_t_ptr(:) !< pointer array: turbulent transfer coefficients for heat
    TYPE(t_ptr_2d3d),ALLOCATABLE :: tfv_t_ptr(:) !< pointer array: laminar reduction factor for evaporation
    TYPE(t_ptr_2d3d),ALLOCATABLE :: tvm_t_ptr(:) !< pointer array: turbulent transfer velocity for momentum
    TYPE(t_ptr_2d3d),ALLOCATABLE :: tvh_t_ptr(:) !< pointer array: turbulent transfer velocity for heat
    TYPE(t_ptr_2d3d),ALLOCATABLE :: tkr_t_ptr(:) !< pointer array: turbulent reference surface diffusion coefficient
    TYPE(t_ptr_2d3d),ALLOCATABLE :: gz0_t_ptr(:) !< pointer array: roughness length * gravity

    TYPE(t_ptr_2d3d),ALLOCATABLE :: tvs_s_t_ptr(:)  !< pointer array: turbulent velocity scale at surface
    TYPE(t_ptr_2d3d),ALLOCATABLE :: tkvm_s_t_ptr(:) !< pointer array: exchange coefficient for momentum at surface
    TYPE(t_ptr_2d3d),ALLOCATABLE :: tkvh_s_t_ptr(:) !< pointer array: exchange coefficient for heat at surface
    TYPE(t_ptr_2d3d),ALLOCATABLE :: u_10m_t_ptr(:)  !< pointer array: zonal wind at 10m
    TYPE(t_ptr_2d3d),ALLOCATABLE :: v_10m_t_ptr(:)  !< pointer array: meridional wind at 10m
    TYPE(t_ptr_2d3d),ALLOCATABLE :: shfl_s_t_ptr(:) !< pointer array: surface sensible heat flux 
    TYPE(t_ptr_2d3d),ALLOCATABLE :: lhfl_s_t_ptr(:) !< pointer array: surface latent heat flux
    TYPE(t_ptr_2d3d),ALLOCATABLE :: umfl_s_t_ptr(:) !< pointer array: u-momentum flux at the surface
    TYPE(t_ptr_2d3d),ALLOCATABLE :: vmfl_s_t_ptr(:) !< pointer array: v-momentum flux at the surface
    TYPE(t_ptr_2d3d),ALLOCATABLE :: qhfl_s_t_ptr(:) !< pointer array: surface moisture flux
    TYPE(t_ptr_2d3d),ALLOCATABLE :: lhfl_bs_t_ptr(:)!< pointer array: lhf from bare soil
    TYPE(t_ptr_2d3d),ALLOCATABLE :: lhfl_pl_t_ptr(:)!< pointer array: lhf from plants
    TYPE(t_ptr_2d3d),ALLOCATABLE :: aerosol_ptr(:)  !< pointer array: prognostic vertically integrated aerosol optical depth

    REAL(wp), POINTER          &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
      , CONTIGUOUS             &
#endif
      &  ::                    &
      &   rain_gsp_rate(:,:),  & !! grid-scale surface rain rate                         [kg/m2/s]
      &   snow_gsp_rate(:,:),  & !! grid_scale surface snow rate                         [kg/m2/s]
      &   ice_gsp_rate(:,:),   & !! grid_scale surface ice rate                          [kg/m2/s]
      &   graupel_gsp_rate(:,:),&!! grid_scale surface graupel rate                      [kg/m2/s]
      &   hail_gsp_rate(:,:),  & !! grid_scale surface hail rate                         [kg/m2/s]
      &   rain_con_rate(:,:),  & !! convective surface rain rate                         [kg/m2/s]
      &   snow_con_rate(:,:),  & !! convective surface snow_rate                         [kg/m2/s]
      &   rain_con_rate_3d(:,:,:),  & !! 3d convective rain rate (convection scheme)     [kg/m2/s]
      &   snow_con_rate_3d(:,:,:),  & !! 3d convective snow_rate (convection scheme)     [kg/m2/s]
      &   rain_edmf_rate_3d(:,:,:), & !! 3d convective rain rate (EDMF scheme)           [kg/m2/s]
      &   snow_edmf_rate_3d(:,:,:), & !! 3d convective snow_rate (EDMF scheme)           [kg/m2/s]
      &   rain_gsp(:,:),       & !! accumulated grid-scale surface rain                  [kg/m2]
      &   snow_gsp(:,:),       & !! accumulated grid_scale surface snow                  [kg/m2]
      &   ice_gsp(:,:),        & !! accumulated grid_scale surface ice                   [kg/m2]
      &   hail_gsp(:,:),       & !! accumulated grid_scale surface hail                  [kg/m2]
      &   graupel_gsp(:,:),    & !! accumulated grid_scale surface graupel               [kg/m2]
      &   rain_con(:,:),       & !! accumulated convective surface rain                  [kg/m2]
      &   snow_con(:,:),       & !! accumulated convective surface snow                  [kg/m2]
      &   tot_prec(:,:),       & !! accumulated grid-scale plus convective surface       [kg/m2]
                                 !! total precipitation
      &   tot_prec_rate_avg(:,:),   & !! average since model start of                    [kg/m2/s]
                                 !! grid-scale plus convective surface 
                                 !! total precipitation rate
      &   con_prec_rate_avg(:,:),   & !! average since model start of                    [kg/m2/s]
                                 !! convective surface precipitation rate
      &   gsp_prec_rate_avg(:,:),   & !! average since model start of                    [kg/m2/s]
                                 !! grid-scale surface precipitation rate
!     the following precipitation variables *0 are accumulated only to the previous call of ww_diagnostics
      &   rain_gsp0(:,:),      & !! accumulated grid-scale surface rain                  [kg/m2]
      &   snow_gsp0(:,:),      & !! accumulated grid_scale surface snow                  [kg/m2]
      &   rain_con0(:,:),      & !! accumulated convective surface rain                  [kg/m2]
      &   snow_con0(:,:),      & !! accumulated convective surface snow                  [kg/m2]
      &   acdnc(:,:,:),        & !! cloud droplet number concentration                   [1/m**3]
      &   cloud_num(:,:),      & !! 2D cloud droplet number concentration for simple aerosol-cloud coupling [1/m**3]
      &   cape    (:,:),       & !! convective available energy
      &   cape_ml (:,:),       & !! convective available energy of mean surface layer parcel
      &   cin_ml  (:,:),       & !! convective inhibition of mean surface layer parcel
      &   con_gust(:,:),       & !! convective gusts near surface
      &   con_udd(:,:,:,:),    & !!(nproma,nlev,nblks,8) convective up/downdraft fields
                                 !! 1= convective updraft mass flux (pmfu)
                                 !! 2= convective downdraft mass flux (pmfd)
                                 !! 3= updraft   detrainment rate  (pmfude_rate)
                                 !! 4= downdraft   detrainment rate (pmfdde_rate)
                                 !! 5= temperature in updraft region (ptu)
                                 !! 6= humidity in updraft region (pqu)
                                 !! 7= condensate in updraft region (plu)
      &  rain_upd(:,:),        & !! total precipitation produced in updrafts [kg/m2/s]
      &  hzerocl(:,:),         & !! height of 0 deg C level [m]
      &  shfl_s(:,:),          & !! sensible heat flux (surface) ( W/m2)
      &  shfl_s_t(:,:,:),      & !! sensible heat flux (surface) ( W/m2)
      &  lhfl_s(:,:),          & !! latent   heat flux (surface) ( W/m2)
      &  lhfl_s_t(:,:,:),      & !! latent   heat flux (surface) ( W/m2)
      &  lhfl_bs(:,:),         & !! latent heat flux from bare soil evap. (surface) ( W/m2)
      &  lhfl_bs_t(:,:,:),     & !! latent heat flux from bare soil evap. (surface) ( W/m2)
      &  lhfl_pl(:,:,:),       & !! latent heat flux from plants                    ( W/m2)
      &  lhfl_pl_t(:,:,:,:),   & !! latent heat flux from plants                    ( W/m2)
      &  qhfl_s(:,:),          & !!      moisture flux (surface) ( Kg/m2/s)
                                 !!      = evaporation rate at surface
      &  qhfl_s_t(:,:,:),      & !! moisture flux (surface)                         ( Kg/m2/s)
                                 !!      = evaporation rate at surface
      &  ashfl_s(:,:),         & !! average or accumulated since model start of shfl_s [W/m2]
      &  alhfl_s(:,:),         & !! average or accumulated since model start of lhfl_s [W/m2]
      &  aqhfl_s(:,:),         & !! average since model start of qhfl_s ( Kg/m2/s) 
                                 !! = average of evaporation rate at surface
      &  alhfl_bs(:,:),        & !! average or accumulated since model start of lhfl_bs [W/m2]
      &  alhfl_pl(:,:,:),      & !! average or accumulated since model start of lhfl_pl [W/m2]
      &  clc(:,:,:),           & !! cloud cover  
      &  clct(:,:),            & !! total cloud cover  
      &  clch(:,:),            & !! cloud cover of high-level clouds
      &  clcm(:,:),            & !! cloud cover of mid-level clouds
      &  clcl(:,:),            & !! cloud cover of low-level clouds
      &  cldepth(:,:),         & !! modified cloud depth for media
      &  clct_mod(:,:),        & !! modified total cloud cover for media
      &  hbas_con(:,:),        & !! height of base of convection [m]
      &  htop_con(:,:),        & !! height of top of convection [m]
      &  htop_dc(:,:),         & !! height above msl of the top of dry convection [m]
      &  tot_cld(:,:,:,:),     & !! total cloud variables (qv,qc,qi)
      &  tot_cld_vi(:,:,:),    & !! vertically integrated tot_cld (qv,qc,qi), including vertical 
                                 !! integrals of qr and qs 
      &  tot_cld_vi_avg(:,:,:),& !! average since model start of the 
                                 !! vertically integrated tot_cld (qv,qc,qi)
      &  clct_avg(:,:),        & !! average since model start of the total cloud cover  
      &  cosmu0(:,:),          & !! cosine of solar zenith angle
      &  albdif(:,:),          & !! Shortwave albedo for diffuse radiation  (0.3-5.0um)
      &  albvisdif(:,:),       & !! UV visible albedo for diffuse radiation (0.3-0.7um)
      &  albvisdir(:,:),       & !! UV visible albedo for direct radiation  (0.3-0.7um)
      &  albnirdif(:,:),       & !! near IR albedo for diffuse radiation    (0.7-5.0um)
      &  albnirdir(:,:),       & !! near IR albedo for direct radiation     (0.7-5.0um)
      &  albdif_t(:,:,:),      & !! tile-based shortwave albedo for diffuse radiation  (0.3-5.0um)
      &  albvisdif_t(:,:,:),   & !! tile-based UV visible albedo for diffuse radiation (0.3-0.7um)
      &  albnirdif_t(:,:,:),   & !! tile-based near IR albedo for diffuse radiation (0.3-0.7um)
      &  vio3(:,:),            & !! vertically integrated ozone amount (Pa O3)
      &  hmo3(:,:),            & !! height of O3 maximum (Pa)
      &  flxdwswtoa(:,:),      & !! downward shortwave flux at TOA [W/m2]
      &  tsfctrad(:,:),        & !! surface temperature at trad [K]
      &  lwflxall(:,:,:),      & !! longwave net flux           [W/m2]
      &  lwflxsfc(:,:),        & !! longwave net flux at surface [W/m2]
      &  lwflx_up_sfc(:,:),    & !! longwave upward flux at surface [W/m2]
      &  lwflx_up_sfc_rs(:,:), & !! longwave upward flux at surface calculated at radiation time steps [W/m2]
      &  lwflxsfc_t(:,:,:),    & !! tile-based longwave net flux at surface [W/m2]
      &  trsolall(:,:,:),      & !! shortwave net tranmissivity (i.e. net flux normalized by irradiance) []
      &  trsolclr_sfc(:,:),    & !! clear-sky shortwave net tranmissivity at the surface
      &  trsol_up_toa(:,:),    & !! normalized shortwave upward flux at the top of the atmosphere
      &  trsol_up_sfc(:,:),    & !! normalized shortwave upward flux at the surface
      &  trsol_par_sfc(:,:),   & !! normalized downward photosynthetically active flux at the surface
      &  trsol_dn_sfc_diff(:,:),& !! normalized shortwave diffuse downward radiative flux at the surface
      &  swflx_up_toa(:,:),    & !! shortwave upward flux at the top of the atmosphere [W/m2]
      &  swflx_up_sfc(:,:),    & !! shortwave upward flux at the surface [W/m2]
      &  swflx_par_sfc(:,:),   & !! shortwave downward photosynthetically active flux at the surface [W/m2]
      &  aswflx_par_sfc(:,:),  & !! shortwave downward photosynthetically active flux at the surface [W/m2]
                                 !! accumulated or mean since model start
      &  swflx_dn_sfc_diff(:,:),& !! shortwave diffuse downward radiative flux at the surface [W/m2]
      &  swflxsfc(:,:),        & !! shortwave net flux at surface [W/m2]
      &  swflxsfc_t(:,:,:),    & !! tile-based shortwave net flux at surface [W/m2]
      &  swflxtoa(:,:),        & !! shortwave net flux at toa [W/m2]
      &  lwflxsfc_a(:,:),      & !! Surface net thermal radiation [W/m2], accumulated or mean since model start
      &  swflxsfc_a(:,:),      & !! Surface net solar radiation [W/m2], accumulated or mean since model start
      &  lwflxtoa_a(:,:),      & !! TOA net thermal radiation [W/m2], accumulated or mean since model start
      &  swflxtoa_a(:,:),      & !! shortwave net flux at toa [W/m2], accumulated or mean since model start
      &  asod_t    (:,:),      & !! Top down solar radiation  [W/m2], accumulated or mean since model start
      &  asou_t    (:,:),      & !! Top up solar radiation  [W/m2], accumulated or mean since model start
      &  athd_s    (:,:),      & !! Surface down thermal radiation [W/m2], accumulated or mean since model start
      &  athu_s    (:,:),      & !! Surface up thermal radiation [W/m2], accumulated or mean since model start
      &  asodird_s (:,:),      & !! Surface down solar direct rad. [W/m2], accumulated or mean since model start 
      &  asodifd_s (:,:),      & !! Surface down solar diff. rad. [W/m2], accumulated or mean since model start 
      &  asodifu_s (:,:),      & !! Surface up solar diff. rad. [W/m2], accumulated or mean since model start 
                                 !! _a means average values if lflux_avg=.TRUE.
                                 !! and accumulated values if lflux_avg=.FALSE., default is .FALSE.
      &  snowlmt     (:,:),    & !! height of snowfall limit above MSL
      &  drag_u_grid (:,:),    & !! zonal resolved surface stress [N/m2]
      &  drag_v_grid (:,:),    & !! meridional resolved surface stress [N/m2]
      &  adrag_u_grid(:,:),    & !! zonal resolved surface stress, accumulated or mean since model start
      &  adrag_v_grid(:,:),    & !! meridional resolved surface stress, accumulated or mean since model start
      &  str_u_sso   (:,:),    & !! zonal sso surface stress [N/m2]
      &  str_v_sso   (:,:),    & !! meridional sso surface stress [N/m2]
      &  astr_u_sso  (:,:),    & !! zonal sso surface stress, accumulated or mean since model start
      &  astr_v_sso  (:,:)       !! meridional sso surface stress, accumulated or mean since model start



    !> Parameter fields for turbulence
    REAL(wp), POINTER      &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
      , CONTIGUOUS         &
#endif
      ::                   &
      rcld(:,:,:)     ,    & !> standard deviation of the saturation deficit    --
      tcm(:,:)        ,    & !! turbulent transfer coefficients for momentum    --
      tch(:,:)        ,    & !! turbulent transfer coefficients for heat        --
      tfm(:,:)        ,    & !! factor of laminar transfer of momentum          --
      tfh(:,:)        ,    & !! factor of laminar transfer of scalars           --
      tfv(:,:)        ,    & !! laminar reduction factor for evaporation        --
      tvm(:,:)        ,    & !! turbulent transfer velocity for momentum      (m/s)
      tvh(:,:)        ,    & !! factor of laminar transfer of scalars           --
      tkr(:,:)        ,    & !! turbulent reference surface diffusion coeff.  (m2/s) (Ustar*kap*z0)
      tkred_sfc(:,:)  ,    & !! reduction factor for minimum diffusion coefficients near the surface
      pat_len(:,:)    ,    & !! length scale of sub-grid scale roughness elements (m)
      gz0(:,:),            & !! roughness length * g of the vertically not
                             !! resolved canopy                               (m2/s2)
      tkvm(:,:,:),         & !! turbulent diffusion coefficients for momentum (m/s2 )
      tkvh(:,:,:),         & !! turbulent diffusion coefficients for heat     (m/s2 )
      t_2m(:,:)       ,    & !! temperature in 2m                             (  K  )
      t_2m_land(:,:)  ,    & !! temperature in 2m (land tiles only)           (  K  )
      tmax_2m(:,:)    ,    & !! maximum temperature in 2m (for specified timerange) ( K )
      tmin_2m(:,:)    ,    & !! minimum temperature in 2m (for specified timerange) ( K )
      qv_2m (:,:)     ,    & !! specific water vapor content in 2m            (kg/kg)
      td_2m (:,:)     ,    & !! dew-point in 2m                               (  K  )
      rh_2m (:,:)     ,    & !! relative humidity in 2m                       (  %  )
      td_2m_land (:,:),    & !! dew-point in 2m (land tiles only)             (  K  )
      rh_2m_land (:,:),    & !! relative humidity in 2m  (land tiles only)    (  %  )
      u_10m (:,:)     ,    & !! zonal wind in 10m                             ( m/s )
      v_10m (:,:)     ,    & !! meridional wind in 10m                        ( m/s )
      sp_10m(:,:)     ,    & !! wind speed in 10m                             ( m/s )
      dyn_gust(:,:)   ,    & !! dynamic gust at 10m                           ( m/s )
      gust10(:,:)     ,    & !! max. gust at 10m                              ( m/s )
      edr   (:,:,:)    ,   & !! eddy dissipation rate
      tcm_t(:,:,:)     ,   & !! turbulent transfer coefficients for momentum    --
      tch_t(:,:,:)     ,   & !! turbulent transfer coefficients for heat        --
      tfv_t(:,:,:)     ,   & !! laminar reduction factor for evaporation        --
      tvm_t(:,:,:)     ,   & !! turbulent transfer velocity for momentum      ( m/s )
      tvh_t(:,:,:)     ,   & !! turbulent transfer velocity for heat          ( m/s )
      tkr_t(:,:,:)     ,   & !! turbulent reference surface diffusion coeff.  ( m2/s) (Ustar*kap*z0)
      gz0_t(:,:,:)     ,   & !! roughness length * g                          (m2/s2)
      tvs_s_t(:,:,:)   ,   & !! surface turbulence velocity scale (SQRT(2*TKE)) (m/s)
                             !! (tile based)
      tkvm_s_t(:,:,:)  ,   & !! surface turbulent diffusion coefficients for momentum (m/s2)
                             !! (tile based)
      tkvh_s_t(:,:,:)  ,   & !! surface turbulent diffusion coefficients for heat (m/s2)
                             !! (tile based)
      u_10m_t(:,:,:)   ,   & !! zonal wind at 10m                             ( m/s )
      v_10m_t(:,:,:)   ,   & !! meridional wind at 10m                        ( m/s )
      umfl_s_t(:,:,:)  ,   & !! u-momentum flux at the surface (tile based)    (N/m2)
      vmfl_s_t(:,:,:)  ,   & !! v-momentum flux at the surface (tile based)    (N/m2)
      umfl_s(:,:)      ,   & !! u-momentum flux at the surface                 (N/m2)
      vmfl_s(:,:)      ,   & !! v-momentum flux at the surface                 (N/m2)
      aumfl_s(:,:)     ,   & !! u-momentum flux at the surface (N/m2), accumulated or mean since model start
      avmfl_s(:,:)           !! v-momentum flux at the surface (N/m2), accumulated or mean since model start
                             !! a means average values if lflux_avg=.TRUE.
                             !! and accumulated values if lflux_avg=.FALSE., default is .FALSE.

    ! need only for EDMF
    REAL(wp), POINTER       &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
      , CONTIGUOUS          &
#endif
      & ::                  &
      & z0m     (:,:)       !< aerodynamic roughness length

    !> Diagnostics for LES turbulence
    REAL(wp), POINTER      &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
      , CONTIGUOUS         &
#endif
      ::                   &
      z_pbl(:,:)     ,     & !> Boundary layer height  (m)
      bruvais(:,:,:) ,     & !> Brunt Vaisala Frequency
      mech_prod(:,:,:),    & !> Mechanical production/loss term in TKE equation
      t_cbase(:,:),        & !>cloud base temperature
      p_cbase(:,:),        & !>cloud base pressure
      t_ctop(:,:),         & !>cloud top temperature
      p_ctop(:,:) !         & !>cloud top pressure
      !cld_opt_thck(

    ! for old aerosol climatology from COSMO (to be used with inwp_radiation==2)
    REAL(wp), POINTER       &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
      , CONTIGUOUS          &
#endif
      & ::                  &
      & aersea  (:,:),      &
      & aerlan  (:,:),      &
      & aerurb  (:,:),      &
      & aerdes  (:,:)

    ! time-interpolated values for Tegen aerosol climatology (needed as state fields for coupling with microphysics and convection)
    REAL(wp), POINTER       &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
      , CONTIGUOUS          &
#endif
      & ::                  &
      & aercl_ss  (:,:),    &
      & aercl_or  (:,:),    &
      & aercl_bc  (:,:),    &
      & aercl_su  (:,:),    &
      & aercl_du  (:,:),    &
      & aerosol   (:,:,:)

    INTEGER, POINTER        &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
      , CONTIGUOUS          &
#endif
      & ::                  &
      &  mbas_con(:,:),     & !< cloud base level index
      &  mtop_con(:,:),     & !< cloud top  level index
      &  ktype   (:,:),     & !< Type of convection
      &  k850    (:,:),     & !< level index that corrsponds to the height 
                              !< of the standard atmosphere 850hPa level above ground
      &  k950    (:,:),     & !< level index that corresponds to the height 
                              !< of the standard atmosphere 950hPa level above ground
      &  k800    (:,:),     & !< level index that corresponds to the height 
                              !< of the standard atmosphere 800hPa level above ground
      &  k400    (:,:),     & !< level index that corresponds to the height 
                              !< of the standard atmosphere 400hPa level above ground
      &  ktop_envel(:,:),   & !< level index of upper boundary of SSO envelope layer
      &  iww     (:,:)        !< significant weather

    REAL(wp), POINTER :: tropics_mask(:,:)      !< mask field that is 1 in the tropics and 0 in the extratropics
    REAL(wp), POINTER :: innertropics_mask(:,:) !< mask field that is 1 in the inner tropics and 0 elsewhere

    LOGICAL, POINTER        &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
      , CONTIGUOUS          &
#endif
      & ::                  &
      & locum     (:,:),    & !< convective  activity indicator
      & ldshcv    (:,:)       !< shallow convection indicator

    !> (Optional:) Additional diagnostic fields:
    REAL(wp), POINTER ::  &
      rh(:,:,:),          &   !> relative humidity
      pv(:,:,:)               !> potential vorticity


    ! Buffer field needed when vertical nesting is combined with a reduced radiation
    ! grid and latm_above_top = .TRUE.
    REAL(wp), POINTER :: buffer_rrg(:,:,:)

    ! Buffer field needed for RTTOV calculations on a vertical nested grid
    REAL(wp), POINTER :: buffer_rttov(:,:,:)

    ! pointer to satellite images (all images in one array):
    REAL(wp), POINTER    :: synsat_arr(:,:,:)

    ! pointers to satellite images (list of 2D slices)
    TYPE (t_rttov_image), ALLOCATABLE :: synsat_image(:)

    !> Special 1D and 0D diagnostics for LES runs
    REAL(wp), ALLOCATABLE :: &
      turb_diag_1dvar(:,:), turb_diag_0dvar(:)  

  END TYPE t_nwp_phy_diag
  !
  ! !---tendencies of type global!
  !
  TYPE t_nwp_phy_tend

    REAL(wp), POINTER           &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
      , CONTIGUOUS              &
#endif
      ::                        &
      ddt_temp_radsw  (:,:,:)  ,& !! Temp-tendency from shortwave radiation
      ddt_temp_radlw  (:,:,:)  ,& !! Temp-tendency from longwave radiation
      ddt_temp_turb   (:,:,:)  ,& !! Temp-tendency from turbulence
      ddt_temp_gscp   (:,:,:)  ,& !! Temp-tendency from microphysics
      ddt_u_turb      (:,:,:)  ,& !! ZonalW-tendency from turbulence
      ddt_u_pconv     (:,:,:)  ,& !! ZonalW-tendency from convective prec
      ddt_v_turb      (:,:,:)  ,& !! MeridW-tendency from turbulence
      ddt_w_turb      (:,:,:)  ,& !! VertW-tendency from turbulence
      ddt_v_pconv     (:,:,:)  ,& !! MeridW-tendency from convective prec
      ddt_tracer_turb (:,:,:,:),& !! Hydromet-tendency from turbulence
      ddt_tracer_pconv(:,:,:,:),& !! Hydromet-tendency from convective prec
      ddt_tke_pconv   (:,:,:)  ,& !! TKE tendency from convective prec
      ddt_tke_hsh     (:,:,:)  ,& !! TKE tendency from horizontal shear
      ddt_tke         (:,:,:)     !! tendency for turbulent velocity scale [m/s^2]

    REAL(vp2), POINTER           &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
      , CONTIGUOUS              &
#endif
      ::                        &
      ddt_temp_drag   (:,:,:)  ,& !! Temp-tendency from sso + gravity-wave drag + Rayleigh friction
      ddt_temp_pconv  (:,:,:)  ,& !! Temp-tendency from convective prec
      ddt_tracer_gscp (:,:,:,:),& !! Hydromet-tendency from microphysics
      ddt_u_gwd       (:,:,:)  ,& !! ZonalW-tendency from gravity wave drag
      ddt_u_sso       (:,:,:)  ,& !! ZonalW-tendency from sso drag
      ddt_v_gwd       (:,:,:)  ,& !! MeridW-tendency from gravity wave drag
      ddt_v_sso       (:,:,:)     !! MeridW-tendency from sso drag


    !Anurag Dipankar, MPIM (2013-May-31)
    !Large-scale tendencies for idealized testcases (nlev)
    REAL(wp), ALLOCATABLE ::   &
      ddt_u_ls       (:), &     !! LS tendency for u 
      ddt_v_ls       (:), &     !! LS tendency for v 
      ddt_temp_ls    (:), &     !! LS tendency for temp 
      ddt_tracer_ls  (:,:)      !! LS tendency for tracer

    TYPE(t_ptr_2d3d),ALLOCATABLE ::  &
      &  tracer_turb_ptr(:)    ,& !< pointer array: one pointer for each component
      &  tracer_conv_ptr(:)    ,& !< pointer array: one pointer for each component
      &  tracer_gscp_ptr(:)       !< pointer array: one pointer for each component

    TYPE(t_ptr_tracer), ALLOCATABLE :: conv_tracer_tend(:,:) !< pointer for chemical tracer conv. tend.

    TYPE(t_ptr_tracer), ALLOCATABLE :: turb_tracer_tend(:,:) !< pointer for chemical tracer turb. tend.

  END TYPE t_nwp_phy_tend

END MODULE mo_nwp_phy_types
