#if (defined (__GNUC__) || defined(__SUNPRO_F95) || defined(__SX__))
#define HAVE_F95
#endif
!>

!! Type definition for the dynamical core of ICONAM.
!!
!! @author Almut Gassmann (MPI-M)
!! @author Daniel Reinert (DWD)
!! @author Guenther Zaengl (DWD)
!!
!! @par Revision History
!! Initial release by Daniel Reinert, DWD (2012-02-07)
!! - Moved here from mo_nonhydro_state to avoid circular dependencies
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!!
MODULE mo_nonhydro_types

  USE mo_kind,                 ONLY: wp, vp
  USE mo_fortran_tools,        ONLY: t_ptr_2d3d, t_ptr_2d3d_vp, t_ptr_tracer
  USE mo_var_list,             ONLY: t_var_list_ptr

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_nh_prog             ! state vector of prognostic variables (type)
  PUBLIC :: t_nh_diag             ! state vector of diagnostic variables (type)
                                  ! on p- and/or z-levels
  PUBLIC :: t_nh_ref              ! state vector of reference state (type)
  PUBLIC :: t_nh_metrics          ! state vector of metrics variables (type)
  PUBLIC :: t_nh_state            ! state vector of nonhydrostatic variables (type)
  PUBLIC :: t_nh_state_lists      ! lists for state vector of nonhydrostatic variables (type)



  ! prognostic variables state vector
  TYPE t_nh_prog

    REAL(wp), POINTER    &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
    , CONTIGUOUS         &
#endif
      ::                 &
      w(:,:,:),          & !> orthogonal vertical wind (nproma,nlevp1,nblks_c)     [m/s]
      vn(:,:,:),         & !! orthogonal normal wind (nproma,nlev,nblks_e)         [m/s]
      rho(:,:,:),        & !! density (nproma,nlev,nblks_c)                     [kg/m^3]
      exner(:,:,:),      & !! Exner pressure (nproma,nlev,nblks_c)                   [-]
      theta_v(:,:,:),    & !! virtual potential temperature (nproma,nlev,nblks_c)    [K]
      tracer(:,:,:,:),   & !! tracer concentration (nproma,nlev,nblks_c,ntracer) [kg/kg]
      tke   (:,:,:)      & !! turbulent kinetic energy                         [m^2/s^2]
        => NULL()          !! (defined on half levels) with 2 time levels
    TYPE(t_ptr_2d3d),ALLOCATABLE :: tracer_ptr(:)  !< pointer array: one pointer for each tracer
    TYPE(t_ptr_tracer),ALLOCATABLE :: conv_tracer(:,:)  
    TYPE(t_ptr_tracer),ALLOCATABLE :: turb_tracer(:,:)  
  END TYPE t_nh_prog


  ! diagnostic variables state vector
  TYPE t_nh_diag

    REAL(wp), POINTER       &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
    , CONTIGUOUS            &
#endif
    &  ::                   &
    ! a) variables needed for intermediate storage and physics-dynamics coupling
    &  u(:,:,:),            & ! zonal wind (nproma,nlev,nblks_c)               [m/s]
    &  v(:,:,:),            & ! meridional wind (nproma,nlev,nblks_c)          [m/s]
    &  omega_z(:,:,:),      & ! relative vertical vorticity at dual grid
                              ! (nproma,nlev,nblks_v)                          [1/s]
    &  vor(:,:,:),          & ! relative vertical vorticity interpolated to cells
                              ! (nproma,nlev,nblks_c)                          [1/s]
    &  ddt_tracer_adv(:,:,:,:), &! advective tendency of tracers          [kg/kg/s]
    &  tracer_vi(:,:,:),    & ! vertically integrated tracers( mass related ones only) [kg/m**2]
    &  tracer_vi_avg(:,:,:),& ! average (since model start) of tracer_vi [kg/m**2]
    &  exner_pr(:,:,:),     & ! exner pressure perturbation, saved from previous step (nproma,nlev,nblks_c)
    &  temp(:,:,:),         & ! temperature (nproma,nlev,nblks_c)                 [K]
    &  tempv(:,:,:),        & ! virtual temperature (nproma,nlev,nblks_c)         [K]
    &  temp_ifc(:,:,:),     & ! temperature at half levels (nproma,nlevp1,nblks_c)[K]
    &  pres(:,:,:),         & ! pressure (nproma,nlev,nblks_c)                  [Pa]
    &  pres_ifc(:,:,:),     & ! pressure at interfaces (nproma,nlevp1,nblks_c)  [Pa]
    &  pres_sfc(:,:),       & ! diagnosed surface pressure (nproma,nblks_c)     [Pa]
    &  pres_sfc_old(:,:),   & ! diagnosed surface pressure at previous timestep (nproma,nblks_c) [Pa]
    &  ddt_pres_sfc(:,:),   & ! current time tendency of diagnosed surface pressure (nproma,nblks_c) [Pa/s]
    &  dpres_mc(:,:,:),     & ! pressure thickness at masspoints(nproma,nlevp,nblks_c)  [Pa]
    &  hfl_tracer(:,:,:,:), & ! horizontal tracer flux at edges             [kg/m/s]
                              ! (nproma,nlev,nblks_e,ntracer)
    &  vfl_tracer(:,:,:,:), & ! vertical tracer flux at cells               [kg/m/s]
                              ! (nproma,nlevp1,nblks_c,ntracer)
    &  div(:,:,:),          & ! divergence(nproma,nlev,nblks_c)     [1/s]
    &  mass_fl_e(:,:,:),    & ! horizontal mass flux at edges (nproma,nlev,nblks_e) [kg/m/s]
    &  rho_ic(:,:,:),       & ! density at half levels (nproma,nlevp1,nblks_c)     [kg/m^3]
    &  theta_v_ic(:,:,:),   & ! theta_v at half levels (nproma,nlevp1,nblks_c)         [K]
    &  airmass_now(:,:,:),  & ! mass of air in layer at physics time step now [kg/m^2]
    &  airmass_new(:,:,:),  & ! mass of air in layer at physics time step new [kg/m^2]

    !
    ! b) variables needed for grid nesting
    &  grf_tend_vn(:,:,:),  & ! vn tendency field for use in grid refinement
                              ! (nproma,nlev,nblks_e)                        [m/s^2]
    &  grf_tend_w(:,:,:),   & ! w tendency field for use in grid refinement
                              ! (nproma,nlevp1,nblks_c)                      [m/s^2]
    &  grf_tend_rho(:,:,:), & ! rho tendency field for use in grid refinement
                              ! (nproma,nlev,nblks_c)                     [kg/m^3/s]
    &  grf_tend_mflx(:,:,:),& ! rho*vn tendency field for use in grid refinement
                              ! (nproma,nlev,nblks_e)                     [kg/m^2/s^2]
    &  grf_bdy_mflx(:,:,:), & ! rho*vn boundary field for use in grid refinement
                              ! (nlev,npoints,2)                            [kg/m^2/s^2]
    &  grf_tend_thv(:,:,:), & ! theta_v tendency field for use in grid refinement
                              ! (nproma,nlev,nblks_c)                          [K/s]
    &  grf_tend_tracer(:,:,:,:), & ! tracer tendency field for use in grid refinement
                                   ! (nproma,nlev,nblks_c,ntracer)          [kg/kg/s]
    &  dvn_ie_int(:,:),         & ! Storage field for vertical nesting: vn at parent interface level
    &  dvn_ie_ubc(:,:),         & ! Storage field for vertical nesting: vn at child upper boundary
    &  w_int(:,:,:),            & ! Storage field for vertical nesting: w at parent interface level
    &  w_ubc(:,:,:),            & ! Storage field for vertical nesting: 
                                  ! average w plus time tendency at child upper boundary
    &  theta_v_ic_int(:,:,:),   & ! Storage field for vertical nesting: theta at parent interface level
    &  theta_v_ic_ubc(:,:,:),   & ! Storage field for vertical nesting: 
                                  ! average theta plus time tendency at child upper boundary
    &  rho_ic_int(:,:,:),       & ! Storage field for vertical nesting: rho at parent interface level
    &  rho_ic_ubc(:,:,:),       & ! Storage field for vertical nesting: 
                                  ! average rho plus time tendency at child upper boundary
    &  mflx_ic_int(:,:,:),      & ! Storage field for vertical nesting: mass flux at parent interface level
    &  mflx_ic_ubc(:,:,:),      & ! Storage field for vertical nesting: 
                                  ! average mass flux plus time tendency at child upper boundary

    !
    ! c) variables derived from analysis increments
    &  t2m_bias (:,:),       & !! filtered T2M bias from surface analysis [K]
    &  rh_avginc(:,:),       & !! time-averaged/filtered RH increments from DA at lowest model level
    &  t_avginc(:,:),        & !! time-averaged/filtered T increments from DA at lowest model level
    &  t_wgt_avginc(:,:),    & !! time-averaged/filtered T increments from DA at lowest model level, weighted with COS(local time)
    &  p_avginc(:,:),        & !! time-averaged/filtered P increments from DA at lowest model level

    !
    ! d) storage variables for time-averaged first-guess output
    &  u_avg    (:,:,:),    & ! normal velocity average          [m/s]
    &  v_avg    (:,:,:),    & ! normal velocity average          [m/s]
    &  pres_avg (:,:,:),    & ! exner average                    [-]
    &  temp_avg (:,:,:),    & ! moist density average            [kg/m^3]
    &  qv_avg   (:,:,:),    & ! specific humidity average        [kg/kg]

    !
    ! e) optional diagnostics
    &  pres_msl(:,:),       & ! diagnosed mean sea level pressure (nproma,nblks_c)  [Pa]
    &  omega(:,:,:),        & ! vertical velocity ( omega=dp/dt )           [Pa/s]
    &  vor_u(:,:,:),        & ! zonal component of relative vorticity
                              ! (nproma,nlev,nblks_c)            [s-1]
    &  vor_v(:,:,:),        & ! meridional component of relative vorticity
                              ! (nproma,nlev,nblks_c)            [s-1]
    &  bvf2(:,:,:),         & ! square of Brunt-Vaisala frequency
                              ! (nproma,nlev,nblks_c)            [s-2]
    &  parcelfreq2(:,:,:)   & ! square of general parcel oscillation frequency
                              ! (nproma,nlev,nblks_c)            [s-2]
    &  => NULL()

    ! d) variables that are in single precision when "__MIXED_PRECISION" is defined
    REAL(vp), POINTER       &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
    , CONTIGUOUS            &
#endif
    &  ::                   &
    ! analysis increments
    &  vn_incr   (:,:,:),   & ! normal velocity increment        [m/s]
    &  exner_incr(:,:,:),   & ! exner inrement                   [-]
    &  rho_incr  (:,:,:),   & ! moist density increment          [kg/m^3]
    &  rhov_incr (:,:,:),   & ! water vapour partial density increment [kg/m^3]
    &  rhoc_incr (:,:,:),   & ! cloud water partial density increment [kg/m^3]
    &  rhoi_incr (:,:,:),   & ! cloud ice partial density increment [kg/m^3]
    &  rhor_incr (:,:,:),   & ! rain partial density increment [kg/m^3]
    &  rhos_incr (:,:,:),   & ! snow partial density increment [kg/m^3]
    &  rhog_incr (:,:,:),   & ! graupel partial density increment [kg/m^3]
    &  rhoh_incr (:,:,:),   & ! hail partial density increment [kg/m^3]
    &  rhonc_incr (:,:,:),   & ! cloud water number density increment [1/m^3]
    &  rhoni_incr (:,:,:),   & ! cloud ice number density increment [1/m^3]
    &  rhonr_incr (:,:,:),   & ! rain number density increment [1/m^3]
    &  rhons_incr (:,:,:),   & ! snow number density increment [1/m^3]
    &  rhong_incr (:,:,:),   & ! graupel number density increment [1/m^3]
    &  rhonh_incr (:,:,:),   & ! hail number density increment [1/m^3]
    !
    ! tendencies, physics increments and derived velocity fields
    &  vt(:,:,:),           & ! tangential wind (nproma,nlev,nblks_e)          [m/s]
    &  ddt_exner_phy(:,:,:),& ! exner pressure tendency from physical forcing 
                              ! (nproma,nlev,nblks_c)                     [1/s]
    &  ddt_vn_phy(:,:,:),   & ! normal wind tendency from forcing
                              ! (nproma,nlev,nblks_e)                          [m/s^2]
    &  exner_dyn_incr(:,:,:), & ! exner pres dynamics increment (nproma,nlev,nblks_c)
    &  vn_ie(:,:,:),        & ! normal wind at half levels (nproma,nlevp1,nblks_e)   [m/s]
    &  w_concorr_c(:,:,:),  & ! contravariant vert correction (nproma,nlevp1,nblks_c)[m/s]
    &  mass_fl_e_sv(:,:,:), & ! storage field for horizontal mass flux at edges (nproma,nlev,nblks_e) [kg/m/s]
    &  ddt_vn_apc_pc(:,:,:,:), & ! normal wind tendency from advection plus coriolis, for predictor (p) and corrector (c) steps
                              ! (nproma,nlev,nblks_e,1:3)                    [m/s^2]
    &  ddt_vn_cor_pc(:,:,:,:), & ! normal wind tendency from coriolis,  for predictor (p) and corrector (c) steps
                              ! (nproma,nlev,nblks_e,1:3)                    [m/s^2]
    &  ddt_w_adv_pc(:,:,:,:),  & ! vert. wind tendency from advection,  for predictor (p) and corrector (c) steps
    ! fields for 3D elements in turbdiff
    &  div_ic(:,:,:),       & ! divergence at half levels(nproma,nlevp1,nblks_c)     [1/s]
    &  hdef_ic(:,:,:),      & ! horizontal wind field deformation (nproma,nlevp1,nblks_c)     [1/s^2]
    &  dwdx(:,:,:),         & ! zonal gradient of vertical wind speed (nproma,nlevp1,nblks_c)     [1/s]
    &  dwdy(:,:,:)          & ! meridional gradient of vertical wind speed (nproma,nlevp1,nblks_c)     [1/s]
    &  => NULL()              ! (nproma,nlevp1,nblks_c,1:3)                  [m/s^2]

    REAL(wp), POINTER       &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
    &  , CONTIGUOUS    &
#endif
    &  ::              &
    ! wind tendencies in dynamics [m/s^2]
    &  ddt_vn_dyn  (:,:,:) => NULL() ,& ! vn   total = sum of the following contributions
    &  ddt_ua_dyn  (:,:,:) => NULL() ,& ! ua
    &  ddt_va_dyn  (:,:,:) => NULL() ,& ! va
    &  ddt_vn_dmp  (:,:,:) => NULL() ,& ! vn   divergence damping
    &  ddt_ua_dmp  (:,:,:) => NULL() ,& ! ua
    &  ddt_va_dmp  (:,:,:) => NULL() ,& ! va
    &  ddt_vn_hdf  (:,:,:) => NULL() ,& ! vn   horizontal diffusion
    &  ddt_ua_hdf  (:,:,:) => NULL() ,& ! ua
    &  ddt_va_hdf  (:,:,:) => NULL() ,& ! va
    &  ddt_vn_adv  (:,:,:) => NULL() ,& ! vn   advection
    &  ddt_ua_adv  (:,:,:) => NULL() ,& ! ua
    &  ddt_va_adv  (:,:,:) => NULL() ,& ! va
    &  ddt_vn_cor  (:,:,:) => NULL() ,& ! vn   Coriolis effect
    &  ddt_ua_cor  (:,:,:) => NULL() ,& ! ua
    &  ddt_va_cor  (:,:,:) => NULL() ,& ! va
    &  ddt_vn_pgr  (:,:,:) => NULL() ,& ! vn   pressure gradient
    &  ddt_ua_pgr  (:,:,:) => NULL() ,& ! ua
    &  ddt_va_pgr  (:,:,:) => NULL() ,& ! va
    &  ddt_vn_phd  (:,:,:) => NULL() ,& ! vn   physics applied in dynamics
    &  ddt_ua_phd  (:,:,:) => NULL() ,& ! ua
    &  ddt_va_phd  (:,:,:) => NULL() ,& ! va
    &  ddt_vn_cen  (:,:,:) => NULL() ,& ! vn   centrifugal acceleration in deep atmosphere
    &  ddt_ua_cen  (:,:,:) => NULL() ,& ! ua
    &  ddt_va_cen  (:,:,:) => NULL() ,& ! va
    &  ddt_vn_iau  (:,:,:) => NULL() ,& ! vn   incremental analysis update
    &  ddt_ua_iau  (:,:,:) => NULL() ,& ! ua
    &  ddt_va_iau  (:,:,:) => NULL() ,& ! va
    &  ddt_vn_ray  (:,:,:) => NULL() ,& ! vn   Rayleigh damping
    &  ddt_ua_ray  (:,:,:) => NULL() ,& ! ua
    &  ddt_va_ray  (:,:,:) => NULL() ,& ! va
    &  ddt_vn_grf  (:,:,:) => NULL() ,& ! vn   grid refinement
    &  ddt_ua_grf  (:,:,:) => NULL() ,& ! ua
    &  ddt_va_grf  (:,:,:) => NULL()    ! va

    LOGICAL :: &
    ! flags indicating whether (=.TRUE.) or not (=.FALSE.) the wind tendency pointers are associated with memory
    &  ddt_vn_dyn_is_associated = .FALSE. ,& ! vn   total = sum of the following contributions
    &  ddt_ua_dyn_is_associated = .FALSE. ,& ! ua
    &  ddt_va_dyn_is_associated = .FALSE. ,& ! va
    &  ddt_vn_dmp_is_associated = .FALSE. ,& ! vn   divergence damping
    &  ddt_ua_dmp_is_associated = .FALSE. ,& ! ua
    &  ddt_va_dmp_is_associated = .FALSE. ,& ! va
    &  ddt_vn_hdf_is_associated = .FALSE. ,& ! vn   horizontal diffusion
    &  ddt_ua_hdf_is_associated = .FALSE. ,& ! ua
    &  ddt_va_hdf_is_associated = .FALSE. ,& ! va
    &  ddt_vn_adv_is_associated = .FALSE. ,& ! vn   advection
    &  ddt_ua_adv_is_associated = .FALSE. ,& ! ua
    &  ddt_va_adv_is_associated = .FALSE. ,& ! va
    &  ddt_vn_cor_is_associated = .FALSE. ,& ! vn   Coriolis effect
    &  ddt_ua_cor_is_associated = .FALSE. ,& ! ua
    &  ddt_va_cor_is_associated = .FALSE. ,& ! va
    &  ddt_vn_pgr_is_associated = .FALSE. ,& ! vn   pressure gradient
    &  ddt_ua_pgr_is_associated = .FALSE. ,& ! ua
    &  ddt_va_pgr_is_associated = .FALSE. ,& ! va
    &  ddt_vn_phd_is_associated = .FALSE. ,& ! vn   physics applied in dynamics
    &  ddt_ua_phd_is_associated = .FALSE. ,& ! ua
    &  ddt_va_phd_is_associated = .FALSE. ,& ! va
    &  ddt_vn_cen_is_associated = .FALSE. ,& ! vn   centrifugal acceleration in deep atmosphere
    &  ddt_ua_cen_is_associated = .FALSE. ,& ! ua
    &  ddt_va_cen_is_associated = .FALSE. ,& ! va
    &  ddt_vn_iau_is_associated = .FALSE. ,& ! vn   incremental analysis update
    &  ddt_ua_iau_is_associated = .FALSE. ,& ! ua
    &  ddt_va_iau_is_associated = .FALSE. ,& ! va
    &  ddt_vn_ray_is_associated = .FALSE. ,& ! vn   Rayleigh damping
    &  ddt_ua_ray_is_associated = .FALSE. ,& ! ua
    &  ddt_va_ray_is_associated = .FALSE. ,& ! va
    &  ddt_vn_grf_is_associated = .FALSE. ,& ! vn   grid refinement
    &  ddt_ua_grf_is_associated = .FALSE. ,& ! ua
    &  ddt_va_grf_is_associated = .FALSE.    ! va

    REAL(vp), POINTER       & ! single precision if "__MIXED_PRECISION" is defined
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
    , CONTIGUOUS            &
#endif
    &  ::                   &
    &  ddt_temp_dyn(:,:,:)  & ! rediagnosed temperature tendency from dynamics [K/s]
    &  => NULL()

    INTEGER, POINTER ::     &
    &  nsteps_avg(:)        &  ! number of time steps summed up for averaging
    &  => NULL()

    REAL(wp), POINTER ::    & !
     &  extra_2d(:,:,:)  ,  & !> extra debug output in 2d and
     &  extra_3d(:,:,:,:)   & !!                       3d
     &  => NULL()

    REAL(vp) :: max_vcfl_dyn=0._vp  ! maximum vertical CFL number in dynamical core

    TYPE(t_ptr_2d3d),ALLOCATABLE ::   &
      &  ddt_grf_trc_ptr(:),   &  !< pointer array: one pointer for each tracer
      &  hfl_trc_ptr    (:),   &  !< pointer array: one pointer for each tracer
      &  vfl_trc_ptr    (:),   &  !< pointer array: one pointer for each tracer
      &  ddt_trc_adv_ptr(:),   &  !< pointer array: one pointer for each tracer
      &  tracer_vi_ptr  (:),   &  !< pointer array: one pointer for each tracer
      &  tracer_vi_avg_ptr(:), &  !< pointer array: one pointer for each tracer
      &  extra_2d_ptr   (:),   &
      &  extra_3d_ptr   (:)

    TYPE(t_ptr_2d3d_vp),ALLOCATABLE ::   &
      &  ddt_vn_apc_pc_ptr (:),&  !< pointer array: one pointer for each step: predictor (p) and corrector (c)
      &  ddt_vn_cor_pc_ptr (:),&  !< pointer array: one pointer for each step: predictor (p) and corrector (c)
      &  ddt_w_adv_pc_ptr  (:)    !< pointer array: one pointer for each step: predictor (p) and corrector (c)

  END TYPE t_nh_diag


  TYPE t_nh_ref
    REAL(wp), POINTER ::    &
      vn_ref  (:,:,:),      & !! orthogonal normal wind (nproma,nlev,nblks_e)      [m/s]
      w_ref   (:,:,:)       & !> orthogonal vertical wind (nproma,nlevp1,nblks_c)  [m/s]
      => NULL()
  END TYPE t_nh_ref


  TYPE t_nh_metrics

    ! Variables that are always in double precision
    REAL(wp), POINTER      &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
    , CONTIGUOUS           &
#endif
     ::                    &
     ! a) General geometric quantities
     !
     z_ifc(:,:,:)        , & ! geometric height at the vertical interface of cells (nproma,nlevp1,nblks_c)
     z_mc(:,:,:)         , & ! geometric height at full levels (nproma,nlev,nblks_c)
     ddqz_z_full(:,:,:)  , & ! functional determinant of the metrics [sqrt(gamma)] (nproma,nlev,nblks_c)
     geopot(:,:,:)       , & ! geopotential at cell center (nproma,nlev,nblks_c)
     geopot_agl(:,:,:)   , & ! geopotential above ground level at cell center (nproma,nlev,nblks_c)
     geopot_agl_ifc(:,:,:),& ! geopotential above ground level at interfaces and cell center (nproma,nlevp1,nblks_c)
     dgeopot_mc(:,:,:)   , & ! geopotential at cell center (nproma,nlev,nblks_c)
     !
     ! b) Specific fields for the dynamical core
     !
     rayleigh_w(:)       , & ! Rayleigh damping on the vertical velocity
     rayleigh_vn(:)      , & ! Rayleigh damping on the normal velocity
     enhfac_diffu(:)     , & ! Enhancement factor for nabla4 background diffusion
     scalfac_dd3d(:)     , & ! Scaling factor for 3D divergence damping terms
     hmask_dd3d(:,:)     , & ! Horizontal mask field for 3D divergence damping terms (nproma,nblks_e)
     vwind_expl_wgt(:,:)  , & ! explicit weight in vertical wind solver (nproma,nblks_c)
     vwind_impl_wgt(:,:)  , & ! implicit weight in vertical wind solver (nproma,nblks_c)
     !
     ! c) Fields for truly horizontal temperature diffusion
     !
     zd_intcoef(:,:) , & 
     zd_geofac(:,:)  , & 
     zd_e2cell(:,:)  , & 
     zd_diffcoef(:)  , & 
     !
     ! d) Fields for LES Model : Anurag Dipankar, MPIM (2013-04)
     !
     ! Vertical grid related
     inv_ddqz_z_full_e(:,:,:)  , & 
     inv_ddqz_z_full_v(:,:,:)  , & 
     inv_ddqz_z_half(:,:,:)    , & 
     inv_ddqz_z_half_e(:,:,:)  , & 
     inv_ddqz_z_half_v(:,:,:)  , & 
     wgtfac_v(:,:,:)           , & 
     ! Mixing length for Smagorinsky model
     mixing_length_sq(:,:,:)   , & 
     !
     ! e) Other stuff
     !
     ! Mask field for mountain or upper slope points
     mask_mtnpoints(:,:) , & ! 
     mask_mtnpoints_g(:,:) , & ! 
     ! slope angle and azimuth (used for slope-dependent radiation)
     slope_angle(:,:)  , & ! [rad]
     slope_azimuth(:,:), & ! [rad]; zero means south-facing slope
     ! Area of subdomain for which feedback is performed; dim: (nlev)
     fbk_dom_volume(:)       &
     => NULL()

    ! Variables that are in single precision when "__MIXED_PRECISION" is defined
    REAL(vp), POINTER      &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
    , CONTIGUOUS           &
#endif
     ::                    &
     ! a) Layer thicknesses
     !
     ddxn_z_full(:,:,:)    , & ! slope of the terrain in normal direction (nproma,nlev,nblks_e)
     ddxn_z_full_c(:,:,:)  , & ! slope of the terrain in normal direction (nproma,nlev,nblks_c)
     ddxn_z_full_v(:,:,:)  , & ! slope of the terrain in normal direction (nproma,nlev,nblks_v)
     ddxn_z_half_e(:,:,:)  , & ! slope of the terrain in normal direction (nproma,nlev,nblks_e)
     ddxn_z_half_c(:,:,:)  , & ! slope of the terrain in normal direction (nproma,nlev,nblks_c)
     ddxt_z_full(:,:,:)    , & ! slope of the terrain in tangential direction (nproma,nlev,nblks_e)
     ddxt_z_full_c(:,:,:)  , & ! slope of the terrain in tangential direction (nproma,nlev,nblks_c)
     ddxt_z_full_v(:,:,:)  , & ! slope of the terrain in tangential direction (nproma,nlev,nblks_v)
     ddxt_z_half_e(:,:,:)  , & ! slope of the terrain in tangential direction (nproma,nlev,nblks_e)
     ddxt_z_half_c(:,:,:)  , & ! slope of the terrain in tangential direction (nproma,nlev,nblks_c)
     ddxt_z_half_v(:,:,:)  , & ! slope of the terrain in tangential direction (nproma,nlev,nblks_v)
     ddqz_z_full_e(:,:,:)  , & ! functional determinant of the metrics [sqrt(gamma)] (nproma,nlev,nblks_e)
     ddqz_z_half(:,:,:)    , & ! functional determinant of the metrics [sqrt(gamma)] (nproma,nlevp1,nblks_c)
     inv_ddqz_z_full(:,:,:), & ! Inverse layer thickness of full levels (nproma,nlev,nblks_c)
     !
     ! b) Interpolation coefficients
     !
     wgtfac_c(:,:,:)      , & ! weighting factor for interpolation from full to half levels (nproma,nlevp1,nblks_c)
     wgtfac_e(:,:,:)      , & ! weighting factor for interpolation from full to half levels (nproma,nlevp1,nblks_e)
     wgtfacq_c(:,:,:)     , & ! weighting factor for quadratic interpolation to surface (nproma,3,nblks_c)
     wgtfacq_e(:,:,:)     , & ! weighting factor for quadratic interpolation to surface (nproma,3,nblks_e)
     wgtfacq1_c(:,:,:)    , & ! weighting factor for quadratic interpolation to model top (nproma,3,nblks_c)
     wgtfacq1_e(:,:,:)    , & ! weighting factor for quadratic interpolation to model top (nproma,3,nblks_e)
     coeff_gradekin(:,:,:), & ! Coefficients for improved discretization of horizontal kinetic energy gradient (nproma,2,nblks_e)
     coeff1_dwdz(:,:,:)   , & 
     coeff2_dwdz(:,:,:)   , & ! Coefficients for second-order accurate dw/dz term (nproma,nlev,nblks_c)
     zdiff_gradp(:,:,:,:) , & ! Height differences between local edge point and neighbor cell points used for
                              ! pressure gradient computation (2,nproma,nlev,nblks_e)
     coeff_gradp(:,:,:,:) , & ! Interpolation coefficients for cubic interpolation of Exner pressure (8,nproma,nlev,nblks_e)
     exner_exfac(:,:,:)   , & ! extrapolation factor for Exner pressure (slope-dependent for stability optimization) 
     !
     ! c) Fields for reference atmosphere
     !
     theta_ref_mc(:,:,:) , & 
     theta_ref_me(:,:,:) , & 
     theta_ref_ic(:,:,:) , & 
     tsfc_ref(:,:)       , & 
     exner_ref_mc(:,:,:) , & 
     rho_ref_mc  (:,:,:) , &  
     rho_ref_me  (:,:,:) , & 
     d_exner_dz_ref_ic(:,:,:), & 
     d2dexdz2_fac1_mc(:,:,:) , & 
     d2dexdz2_fac2_mc(:,:,:) , &
     !
     ! Correction term needed to use perturbation density for lateral boundary nudging
     ! (note: this field is defined on the local parent grid in case of MPI parallelization)
     rho_ref_corr(:,:,:) , &
     !
     ! d) other stuff
     !
     pg_exdist (:)         &  ! extrapolation distance needed for igradp_method = 3
     => NULL()

    INTEGER, POINTER          &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
    , CONTIGUOUS              &
#endif
     ::                       &
     vertidx_gradp(:,:,:,:) , &  ! Vertical index of neighbor points needed for Taylor-expansion-based 
                                 ! pressure gradient (2,nproma,nlev,nblks_e)
     !
     ! Fields for truly horizontal temperature diffusion
     zd_indlist(:,:) , & 
     zd_blklist(:,:) , & 
     zd_edgeidx(:,:) , & 
     zd_edgeblk(:,:) , & 
     zd_vertidx(:,:) , & 
     !
     ! Fields for igradp_method = 3
     pg_edgeidx(:) , & 
     pg_edgeblk(:) , & 
     pg_vertidx(:) , & 
     !
     ! Index lists for grid points on which lateral boundary nudging is applied
     nudge_c_idx(:) , & 
     nudge_e_idx(:) , & 
     nudge_c_blk(:) , & 
     nudge_e_blk(:) , & 
     !
     ! Index lists and mask fields needed to minimize the number of halo communications
     ! a) index lists for halo points belonging to the nest boundary region
     bdy_halo_c_idx(:) , & 
     bdy_halo_c_blk(:) , & 
     ! b) index lists for halo points belonging to a nest overlap region
     !    the additional dimension is n_childdom
     ovlp_halo_c_dim(:)   , &
     ovlp_halo_c_idx(:,:) , & 
     ovlp_halo_c_blk(:,:) , & 
     ! c) index lists for mass fluxes at lateral nest boundary (including the required halo points)
     bdy_mflx_e_idx(:) , & 
     bdy_mflx_e_blk(:)   &
     => NULL()

     !
     ! Vertically varying nudging coefficient: nudgecoeff_vert(nlev)
    REAL(wp), POINTER &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
    , CONTIGUOUS           &
#endif
     :: nudgecoeff_vert(:)

    ! Upper atmosphere/deep atmosphere
    !
    REAL(wp), POINTER     &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
      , CONTIGUOUS        &
#endif
     ::                   &
     zgpot_ifc(:,:,:)   , & ! geopotential height of the grid layer interfaces (nproma,nlevp1,nblks_c)
     zgpot_mc(:,:,:)    , & ! geopotential height of cell centers (nproma,nlev,nblks_c)
     dzgpot_mc(:,:,:)     & ! geopotential layer thickness (nproma,nlev,nblks_c)   
     => NULL()
    !
    REAL(wp), POINTER      &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
      , CONTIGUOUS         &
#endif
      ::                   &
      deepatmo_t1mc(:,:),  & ! metrical modification factors for full levels (nlev,nitem)
      deepatmo_t1ifc(:,:), & ! metrical modification factors for half levels (nlevp1,nitem)
      deepatmo_t2mc(:,:)   & ! metrical modification factors for full levels (nitem,nlev)
      => NULL()    


   ! Corresponding scalar list dimensions
   INTEGER  :: zd_listdim  ! for truly horizontal temperature diffusion
   INTEGER  :: pg_listdim  ! for igradp_method = 3
   INTEGER  :: nudge_c_dim, nudge_e_dim ! for grid points on which lateral boundary nudging is applied
   INTEGER  :: bdy_halo_c_dim ! for halo points belonging to the nest boundary region
   INTEGER  :: bdy_mflx_e_dim ! for mass fluxes at lateral nest boundary


   ! Finally, a mask field that excludes boundary halo points
   LOGICAL,  POINTER :: mask_prog_halo_c(:,:) => NULL()

  END TYPE t_nh_metrics


!-------------------------------------------------------------------------
!                      STATE VECTORS AND LISTS
!-------------------------------------------------------------------------
  TYPE t_nh_state

    !array of prognostic states at different timelevels
    TYPE(t_nh_prog),  ALLOCATABLE :: prog(:)       !< shape: (timelevels)
    TYPE(t_nh_diag)    :: diag
    TYPE(t_nh_ref)     :: ref
    TYPE(t_nh_metrics) :: metrics

  END TYPE t_nh_state

  TYPE t_nh_state_lists

    ! array of prognostic state lists at different timelevels
    ! splitting this out of t_nh_state allows for a deep copy 
    ! of the p_nh_state variable to accelerator devices with OpenACC
    TYPE(t_var_list_ptr), ALLOCATABLE :: prog_list(:)  !< shape: (timelevels)
    TYPE(t_var_list_ptr)   :: diag_list
    TYPE(t_var_list_ptr)   :: ref_list
    TYPE(t_var_list_ptr)   :: metrics_list
    TYPE(t_var_list_ptr), ALLOCATABLE :: tracer_list(:) !< shape: (timelevels)

  END TYPE t_nh_state_lists


END MODULE mo_nonhydro_types




