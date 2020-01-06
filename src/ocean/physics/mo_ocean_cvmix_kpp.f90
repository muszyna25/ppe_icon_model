!>
!! Provide an implementation of the ocean physics.
!!
!! Provide an implementation of the physical parameters and characteristics
!! for the hydrostatic ocean model.
!!
!! This module defines the interface to the KPP scheme of the CVMix library (cvmix_kpp.f90).
!! The interface calculates all necessary additional fields that are required as input to the
!! KPP scheme. The interface follows closely that of MOM6.
!!
!! A detailed description can be found in the CVMix documentation:
!!   Griffies, S.M., Levy, M., Adcroft, A. J., Danabasoglu, R., Hallberg, R. W.,
!!     Jacobsen, D., Large, W., Ringler, T. D.: Theory and numerics of the
!!     Community Ocean Vertical Mixing (CVMix) Project. Tech. rep., available
!!     at: https://github.com/CVMix/CVMix-description, 2013.
!!
!! reference: LMD94: Large, W. G., McWilliams, J. C., Doney, S. C.: Oceanic vertical mixing: A review and a model
!!                     with a nonlocal boundary layer parameterization. Rev. Geophys., 21, 363-403,
!!                     https://doi.org/10.1029./94RG01872,1994.
!!
!! Main subroutines:
!! setup_kpp(): this prepares all parameters for the KPP scheme.
!! calc_kpp(): this prepares all input fields and calls calc_coeff_kpp() from the CVMix library.
!!
!!
!! @author Oliver Gutjahr, MPI-M
!!
!! @par Revision History
!!  Original version by Oliver Gutjahr, MPI-M (04-2019)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!----------------------------
#include "omp_definitions.inc"
!----------------------------
MODULE mo_ocean_cvmix_kpp

  USE mo_kind,                           ONLY:  &
      wp                                          ! precision

  USE cvmix_kpp,                         ONLY:  & ! CVMix-library:
      cvmix_init_kpp                           ,& !
      cvmix_coeffs_kpp                         ,& ! 
      cvmix_kpp_compute_OBL_depth              ,& !
      cvmix_kpp_compute_turbulent_scales       ,& ! calculates vertical turbulent velocity scales (momentum, tracer)
      cvmix_kpp_compute_bulk_Richardson        ,& !
      cvmix_kpp_compute_unresolved_shear       ,& !
      cvmix_kpp_params_type                    ,& !
      cvmix_kpp_compute_kOBL_depth                !

  USE cvmix_math,                        ONLY:  &
      cvmix_math_evaluate_cubic                   ! function to calculate G(sigma)

  USE mo_ocean_nml,                      ONLY:  & ! --> parameters from namelist
      n_zlev                                   ,& ! number of vertical levels
      velocity_VerticalDiffusion_background    ,& ! background mixing value for viscosity
      Temperature_VerticalDiffusion_background ,& ! background mixing value for temperature
      Salinity_VerticalDiffusion_background    ,& ! background mixing value for salinity
      no_tracer                                ,& ! ???
      tracer_convection_MixingCoefficient      ,& ! enhanced diffusivity constant during convective conditions
      RichardsonDiffusion_threshold            ,& ! Threshold for Ri-number below convection parameterisation is applied
      convection_InstabilityThreshold          ,& ! Threshold for density gradient below convection parameterisation is applied
      tracer_RichardsonCoeff                   ,& ! ???
      velocity_RichardsonCoeff                 ,& ! ???
      use_reduced_mixing_under_ice             ,& ! flag if wind mixing under sea ice is active 
      EOS_TYPE                                 ,& ! equation of state type
      LinearThermoExpansionCoefficient         ,& ! alpha
      LinearHalineContractionCoefficient       ,& ! beta
      ReferencePressureIndbars                 ,& ! reference pressure in bar
      OceanReferenceDensity                    ,& ! reference ocean density (rho_0)
      Ri_crit                                  ,& ! critical bulk-Rickardson number (Rib) used to diagnose OBL depth
      vonKarman                                ,& ! von-Karman constant (dimensionless)
      cs                                       ,& ! parameter for computing velocity scale function (dimensionless)
      cs2                                      ,& ! parameter for multiplying by non-local term
      minOBLdepth                              ,& ! if non-zero, sets the minimum depth for the OBL (m)
      minVtsqr                                 ,& ! minimum for the squared unresolved velocity used in Rib CVMix calculation (m2/s2)
      langmuir_Efactor                         ,& ! Langmuir enhancement factor for turb. vertical velocity scale (w_s)
      surf_layer_ext                           ,& ! fraction of OBL depth considered as the surface layer (nondim) (epsilon=0.1)
      enhance_diffusion                        ,& ! True => add enhanced diffusivity at base of boundary layer (not recommended)
      computeEkman                             ,& ! True => compute Ekman depth limit for OBLdepth (not tested) 
      computeMoninObukhov                      ,& ! True => compute Monin-Obukhov limit for OBLdepth (not tested)
      llangmuirEF                              ,& ! True => apply Langmuir enhancement factor to w_s (not tested)
      lenhanced_entr                           ,& ! True => enhance entrainment by adding Stokes shear to unresolved shear (not tested)
      CS_is_one                                ,& ! True => forces cs2 (cs in NLT computation) to equal 1 for parabolic non-local option.
      interpType                               ,& ! Type of interpolation in determining OBL depth: linear,quadratic,cubic
      MatchTechnique                           ,& ! Method used in CVMix for setting diffusivity and NLT profile functions
      internal_mix_scheme                      ,& ! Ri-number dependent mixing scheme below the OBL: 'PP' or 'KPP'
      lnonlocal_trans                          ,& ! True => non-local transport (NLT) terms are applied (disable only for testing)
      fixedOBLdepth                            ,& ! True => fixed OBL depth at fixedOBLdepth_value (only for testing)
      lconvection                              ,& ! True => convection param. 'enhanced diff.' is called below the OBL
      lnl_trans_under_sea_ice                  ,& ! True => non-local transport (NLT) tendencies are calculated below sea ice
      diag_WS                                  ,& ! True => diagnose turb. velocity scale (w_s) with correct OBL depth (not required?)
      diag_vtsqr                               ,& ! True => vert. turb. shear acting on OBL depth will be diagnosed (not working yet)
      diag_G                                   ,& ! True => diagnose non-dimensional shape function G(sigma)=sigma(1-sigma)**2
      min_thickness                            ,& ! min. thickness (m) to avoid division by small numbers in vicinity of vanished layers
      deepOBLoffset                            ,& ! If non-zero, a distance from the bottom that the OBL cannot penetrate through (m)
      fixedOBLdepth_value                      ,& ! value for the fixed OBL depth when fixedOBLdepth=.True.
      SW_METHOD                                ,& ! Sets method for using shortwave radiation in surface buoyancy flux
      NLT_shape                                ,& ! Use a different shape function (G) for non-local transport (NLT)
      KPP_nu_zero                              ,& ! leading coefficient of shear mixing formula (m^2/s; default= 5e-3)
      KPP_Ri_zero                              ,& ! critical Richardson number value (0.7 in LMD94)
      KPP_loc_exp                              ,& ! Exponent of unitless factor of diffusities (3.0 in LMD94)
      PP_nu_zero                               ,& ! leading coefficient of shear mixing in PP formula (m^2/s; default= 1e-2)
      PP_alpha                                 ,& ! coefficient in PP scheme
      PP_loc_exp                              !,& ! coefficient in PP scheme 


  USE mo_ocean_physics_types,            ONLY:  & ! --> contains variables for physics
      t_ho_params                                 ! ???
  USE mo_parallel_config,                ONLY:  &
      nproma
  USE mo_model_domain,                   ONLY:  &
      t_patch                                  ,& ! ???
      t_patch_3d                                  ! ???
  USE mo_impl_constants,                 ONLY:  &
      max_char_length                             ! maximum number of characters
  USE mo_util_dbg_prnt,                  ONLY:  & ! 
      dbg_print                                ,& ! ???
      debug_print_MaxMinMean                      ! ???
  USE mo_ocean_types,                    ONLY:  &
      t_hydro_ocean_state                      ,& ! ???
      t_onEdges_Pointer_3d_wp                  ,& ! ???
      t_onCells_HalfLevels_Pointer_wp          ,& ! ???
      t_operator_coeff                            ! ??? 
  USE mo_ocean_state,                    ONLY:  &
      oce_config                                  ! ???
  USE mo_physical_constants,             ONLY:  & !--> physical constants
      grav                                     ,& ! gravitational constant
      sal_ref                                  ,& ! reference salinity
      clw                                         ! specific heat capacity of liquid water [J/K/kg]
  USE mo_dynamics_config,                ONLY:  &  
      nold                                        ! old timestep
  USE mo_run_config,                     ONLY:  & ! --> not used 
      dtime                                       ! numerical timestep
  USE mo_linked_list,                    ONLY:  &
      t_var_list                                  ! ???
  !USE mo_zaxis_type,                     ONLY:  &
  !    za_depth_below_sea                       ,& ! ???
  !    za_depth_below_sea_half                  ,& ! ???
  !    za_surface                                  ! ???
  USE mo_grid_subset,                    ONLY:  &
      t_subset_range                           ,& ! ???
      get_index_range                             ! ???
  USE mo_sync,                           ONLY:  &
      sync_c                                   ,& ! ???
      sync_e                                   ,& ! ???
      sync_v                                   ,& ! ???
      sync_patch_array                         ,& ! ???
      global_max                               ,& ! ???
      sync_patch_array_mult                       ! ???
  USE mo_ocean_thermodyn,                ONLY:  &
      calculate_density_onColumn                  ! function to calculate in-situ density over 1d-column
  USE mo_math_constants,                 ONLY:  &
      dbl_eps                                     ! almost zero
  USE mo_statistics,                     ONLY:  &
      global_minmaxmean
  USE mo_sea_ice_types,                  ONLY:  &
      t_sea_ice                                ,& ! sea ice temperature ?
      t_atmos_fluxes                              ! atmospheric fluxes over sea ice ?
  USE mo_ocean_surface_types,            ONLY:  &
      t_ocean_surface                             ! contains p_oce_sfc 
  USE mo_ocean_thermodyn,                ONLY:  &
      calc_neutralslope_coeff_func_onColumn    ,& ! calculates alpha and beta (1968)
      calc_neutralslope_coeff_func_onColumn_UNESCO ! calculates alpha and beta (UNESCO)


IMPLICIT NONE
  !PRIVATE

  PUBLIC :: calc_kpp
  PUBLIC :: setup_kpp

  !> CVmix parameters
  !type(CVmix_kpp_params_type), public, pointer :: KPP_params => NULL()
  !type KPP_CS
  !-------------------------------------------------
  ! Comments on and options for some namelist parameters
  !-------------------------------------------------
  !
  !  surface_layer_ext=0.1          Defines the depth of the surface layer (SL) by surface_layer_ext*OBLdepth. 
  !                                 Usually this value is never changed in literature.
  !
  !  Ri_crit=0.3                    Critical bulk Richardson number. Defines OBL depth where RiB>Ri_crit. Commonly
  !                                 used values are in the range 0.25-0.7.    
  ! 
  !  vonKarman=0.40                 von-Karman constant (dimensionless). Universally, 0.40 is used, but 0.41 is quite 
  !                                 common in atmospheric applications.
  !
  !  NLT_shape="CVMIX"              Overwrites NLT shape function if not NLT_shape="CVMIX"
  !                                 Allowed values are:
  !                                   CVMIX     - Uses the profiles from CVmix specified by MATCH_TECHNIQUE
  !                                   LINEAR    - A linear profile, G(sigma) = 1-sigma
  !                                   PARABOLIC - A parablic profile, G(sigma) = (1-sigma)^2
  !                                   CUBIC     - A cubic profile, G(sigma) =  (1-sigma)^2(1+2*sigma)
  !                                   CUBIC_LMD - The original KPP profile
  !                                   default='CVMIX'.
  !
  !  cs2=6.32739901508              Parameter for multiplying by non-local term. This is active for NLT_SHAPE_CUBIC_LMD only.
  !                                 Note: MOM6 recommends NLT_shape="PARABOLIC", which
  !                                 results in deeper OBL and less spurios extremes
  !
  !  MatchTechnique="SimpleShapes"  Method used in CVMix for setting diffusivity and NLT profile functions:
  !                                   SimpleShapes      - G(sigma) = sigma*(1-sigma)^2 for both diffusivity and NLT
  !                                   MatchGradient     - G(sigma) = sigma*(1-sigma)^2 for NLT; diffusivity profile from matching
  !                                   MatchBoth         - match gradient for both diffusivity and 
  !                                   ParabolicNonLocal - G(sigma) = sigma*(1-sigma)^2 for diffusivity;
  !                                                       G(sigma) = (1-sigma)^2 for NLT
  !
  !
  !  SW_METHOD="SW_METHOD_ALL_SW"   Sets method for using shortwave radiation in surface buoyancy flux (Bf)
  !                                 Alternatives:
  !                                    SW_METHOD_ALL_SW -
  !                                    SW_METHOD_MXL_SW -
  !                                    SW_METHOD_LV1_SW -
  !
  !  interpType="cubic"             Type of interpolation in determining OBL depth: 
  !                                    linear
  !                                    quadratic
  !                                    cubic  
  !  
  !  internal_mix_scheme="PP"       Ri-number dependent mixing scheme below the OBL: 
  !                                    PP  - Pacanowski & Philander (PP), 1981 scheme
  !                                    KPP - scheme from Large et al. (1994) = LMD94 (not fully tested yet!)
  !
  !  CS_is_one=.FALSE.              If True, forces cs2 (cs in NLT computation) to equal 1 for parabolic non-local option.
  !
  !  lconvection=.TRUE.             If True, convection below the OBL is parameterized as enhanced diffusivity (but not viscosity)
  !
  !  lnl_trans_under_sea_ice=.TRUE. If True, non-local transport tendencies are calculated below sea ice. Might be removed in favour of
  !                                 'use_reduced_mixing_under_ice' flag.
  !  
  !  use_reduced_mixing_under_ice=.TRUE.   This affects whether u* is reduced with increasing sea ice cover: u*=u* * (1-A)^2.
  !
  !  enhance_diffusion=.FALSE.      True => add enhanced diffusivity at base of boundary layer (not recommended).
  !
  !  computeEkman=.FALSE.           True => compute Ekman depth limit for OBLdepth.
  !  computeMoninObukhov=.FALSE.    True => compute Monin-Obukhov limit for OBLdepth.
  !  llangmuirEF=.FALSE.            True => apply Langmuir enhancement factor (langmuir_Efactor) to w_s/w_m; 
  !                                         w_s=w_s*langmuir_Efactor and w_m=w_m*langmuir_Efactor.  
  !  lenhanced_entr=.FALSE.         True => enhance entrainment by adding Stokes shear.
  !
  !  diag_WS=.FALSE.                If True, the turbulent velocity scale for tracers (w_s) will be recalculated with correct OBL depth.
  !
  !  diag_vtsqr=.TRUE.              If True, vertical turbulent shear acting on OBL depth will be diagnosed (not working a.t.m).
  !  
  !  diag_G=.TRUE.                  The shape function G(sigma) is not available as output from CVMix, if set to True, 
  !                                 G(sigma)=sigma(1-sigma)**2 is diagnosed. Note if NLT_SHAPE!="CVMIX" then this diagnosed G will differ 
  !                                 from the G that is used to calculate NLT.
  !
  !  minVtsqr=1e-10                 Min for the squared unresolved velocity used in Rib CVMix calculation (m2/s2), see p.XX in the CVMiX documentation. 
  !  -------------------------------------------------------------------------------------
  !  > Namelists that are only used for testing and debugging (and can be deleted later):
  !  
  !  lnonlocal_trans=.TRUE.         If True, non-local transport terms are calculated, otherwise they are zero
  !  minOBLdepth=0.0                If non-zero, set the minimum depth for the OBL (m). Note that the minimum OBLdepth is always the depth of the first
  !                                 interface level below the surface.  
  !  deepOBLoffset=0.0              If non-zero, is a distance from the bottom that the OBL can not penetrate through (m).  
  !  fixedOBLdepth=.FALSE.          If True, will fix the OBL depth at 'fixedOBLdepth_value'.
  !----------------------------------------------------------------------------------------------------------------------------

  !> CVmix variables
     ! pointer for convenience 
     REAL(wp), POINTER :: dz(:,:,:)            ! cell vertical distances
     REAL(wp), POINTER :: dzi(:,:,:)           ! 1/dz
     REAL(wp), POINTER :: dzw(:,:,:)           ! cell thicknesses
     REAL(wp), POINTER :: temp(:,:,:)          ! temperature
     REAL(wp), POINTER :: salt(:,:,:)          ! salinity  
     REAL(wp), POINTER :: Av_old(:,:,:)        ! Viscosity of previous time step 
     REAL(wp), POINTER :: kv_old(:,:,:)        ! Diffusivity of previous time step

     INTEGER, POINTER :: kbot(:,:)             ! model level of bottom layer 

    ! pointer for KPP variables - 2d
    REAL(wp), POINTER :: OBLdepth(:,:)               ! Depth (positive) of OBL (m)     
    REAL(wp), POINTER :: Usurf(:,:)                  ! Surface layer averaged velocity (m/s)
    REAL(wp), POINTER :: ustar(:,:)                  ! friction velocity (ms-1)  
    REAL(wp), POINTER :: SLdepth_2D(:,:)             ! surface layer depth (m)
    REAL(wp), POINTER :: Coriolis(:,:)               ! Coriolis parameter at cell centre
    REAL(wp), POINTER :: kOBL_2D(:,:)                ! model layer of maximum mixed layer depth
    REAL(wp), POINTER :: net_heat(:,:)               ! net heat flux (Wms-1)
    REAL(wp), POINTER :: net_salt(:,:)               ! net freshwater flux (?)

    ! pointer for KPP variables - 3d
    REAL(wp), POINTER :: avo_kpp(:,:,:)              ! viscosity from KPP scheme
    REAL(wp), POINTER :: dvo_heat_kpp(:,:,:)         ! diffusivity for heat KPP scheme
    REAL(wp), POINTER :: dvo_salt_kpp(:,:,:)         ! diffusivity for salt KPP scheme
    REAL(wp), POINTER :: richardson_no(:,:,:)        ! Bulk Richardson number
    REAL(wp), POINTER :: BulkRi(:,:,:)               ! non-local bulk Richardson number for each layer (dimensionless)
    REAL(wp), POINTER :: dRho(:,:,:)                 ! non-local density gradient
    REAL(wp), POINTER :: Uz2(:,:,:)                  ! non-local square of bulk difference in resolved velocity (m2/s2)
    REAL(wp), POINTER :: Nbuoy(:,:,:)                ! buoyancy frequency (unit?)
    REAL(wp), POINTER :: N2(:,:,:)                   ! squared buoyancy frequency
    REAL(wp), POINTER :: WS_cntr(:,:,:)              ! turbulent velocity scale for scalars (m/s) at cell centres
    REAL(wp), POINTER :: WS(:,:,:)                   ! turbulent velocity scale for scalars (m/s)
    REAL(wp), POINTER :: WM(:,:,:)                   ! turbulent velocity scale for momentum (m/s)
    REAL(wp), POINTER :: stab_KPP(:,:,:)             ! vertical stability (drho/dz) 
    REAL(wp), POINTER :: cellHeight_KPP(:,:,:)       ! depths of cell center
    REAL(wp), POINTER :: iFaceHeight_KPP(:,:,:)      ! depths of interfaces
    REAL(wp), POINTER :: BuoyFlux_3D(:,:,:)          ! buoyancy flux (m2/s3)
    REAL(wp), POINTER :: Tsurf(:,:,:)                ! temperature averaged over surface layer
    REAL(wp), POINTER :: Ssurf(:,:,:)                ! salinity averaged over surface layer
    REAL(wp), POINTER :: G_3D(:,:,:)                 ! shape function within OBL: G(sigma)
    REAL(wp), POINTER :: Vt2(:,:,:)                  ! vertical turbulent velocity shear
    REAL(wp), POINTER :: nl_trans_tend_heat(:,:,:)   ! non-local heat transport tendency: cs2 * G(sigma) * Qh
    REAL(wp), POINTER :: nl_trans_tend_salt(:,:,:)   ! non-local salt transport tendency: cd2 * G(sigma) * Qs
    REAL(wp), POINTER :: nonLocalTransHeat(:,:,:)    ! non-local heat transport: cs2 * G(sigma)
    REAL(wp), POINTER :: nonLocalTransScalar(:,:,:)  ! non-local scalar transport: cs2 * G(sigma)

    REAL(wp), POINTER :: vert_density_grad(:,:,:)    ! vertical density gradient 
    !type(KPP_CS), pointer       :: KPP_params       !< Control structure

    TYPE(t_subset_range), POINTER         :: cells_in_domain

    ! debugging variables
    CHARACTER(LEN=12)       :: str_module    = 'ocePhys_KPP'  ! Output of module for 1 line debug
    INTEGER :: idt_src       = 0                              ! Level of detail for 1 line debug
    

CONTAINS

  !-------------------------------------------------------------------------
  !>
  !! Initialise the KPP scheme
  !!
  !! This routine sets all required parameters for internal calculations with
  !! routine cvmix_coeffs_kpp
  !!  !!
  !! @par Revision History
  !! Initial release by Oliver Gutjahr, MPI-M (2019-04)
  !-------------------------------------------------------------------------
  SUBROUTINE setup_kpp()
    !---------------------------
    ! initialises the KPP scheme
    !---------------------------
    IF (MatchTechnique .EQ. 'ParabolicNonLocal') THEN
       ! This forces Cs2 (Cs in non-local computation) to equal 1 for parabolic
       ! non-local option.
       !  May be used during CVmix initialization.
       Cs_is_one=.TRUE.
    ENDIF

    !========================================================================
    ! CALL initialisation routine of KPP from CVMix library
    !========================================================================

    ! Call the cvmix subroutine to initialise all required namelists (checks
    ! whether they are in the allowed range).
    ! Note: the namelist-parameters for KPP could be stored in the pointer
    ! 'KPP_params' (not working a.t.m.).
    CALL cvmix_init_kpp( Ri_crit        = Ri_crit,             &
                         minOBLdepth    = minOBLdepth,         &
                         minVtsqr       = minVtsqr,            &
                         vonKarman      = vonKarman,           &
                         surf_layer_ext = surf_layer_ext,      &
                         interp_type    = interpType,          &
                         lEkman         = computeEkman,        &
                         lMonOb         = computeMoninObukhov, &
                         MatchTechnique = MatchTechnique,      &
                         lenhanced_diff = enhance_diffusion   ,&
                         llangmuirEF    = llangmuirEF         ,&
                         lenhanced_entr = lenhanced_entr      ,&
                         lnonzero_surf_nonlocal = Cs_is_one   )!,& 
                         !cvmix_kpp_params_user = KPP_params)

  END SUBROUTINE setup_kpp


!========================================================================
!========================================================================
  !-------------------------------------------------------------------------
  !>
  !! Run the KPP scheme
  !!
  !! This routine is a wrapper that calculates all additional fields required
  !! for the KPP scheme and finally runs the routine cvmix_coeffs_kpp to
  !! calculate the vertical mixing coefficients in the mixed layer.
  !! 
  !! Below the mixed layer there are two options: either using a 
  !! Richardson-number dependent scheme of Large et al. (1994), or using the 
  !! PP-scheme. Enhanced diffusivity for convection is active by default.
  !!
  !! Note 1: the nonlocal transport terms are not applied in this routine, has
  !!         to be done somewhere else.
  !! Note 2: there is no matching of the mixing coefficients at the base of the
  !!         mixed layer by default.
  !! 
  !! @par Revision History
  !! Initial release by Oliver Gutjahr, MPI-M (2019-04)
  !-------------------------------------------------------------------------

SUBROUTINE calc_kpp(patch_3d, ocean_state, params_oce, atmos_fluxes, p_oce_sfc, concsum)
  TYPE(t_patch_3d ),TARGET, INTENT(in)  :: patch_3d
  !TYPE(t_subset_range), POINTER         :: edges_in_domain 
  TYPE(t_hydro_ocean_state), TARGET     :: ocean_state
  TYPE(t_atmos_fluxes)                  :: atmos_fluxes
  TYPE(t_ho_params), INTENT(inout)      :: params_oce
  TYPE (t_ocean_surface), INTENT(IN)    :: p_oce_sfc
  REAL(wp), TARGET                      :: concsum(:,:) ! t_sea_ice%concs
  TYPE(t_patch),POINTER                 :: patch_2D


    ! Begin calculation of KPP
    ! 
    ! -- Part 1: within mixed layer --
    !
    ! 1) Calculate the surface layer depth, averaged surface layer quantities, and vertical shear
    ! 2) Calculate buoyancy flux (m2/s3)
    ! 3) Compute in-situ density and buoyancy frequency (N2)
    ! 4) Calculate friction velocity (ustar) at surface (m/s)
    ! 5) Calculate the turbulent velocity scales w_s and w_m at interface depths
    ! 6) Calculate non-local bulk Richardson number
    ! 7) Compute OBL depth (m)
    ! 8) Call CVMix/KPP to obtain OBL diffusivities (k), OBL depth (h), and non-dim. shape function (G)
    ! 9) Update non-local transport terms and non-local tendency terms (if NLT_shape != "CVMIX")
    !10) Diagnostics for KPP: WS,WM, G, and Vt2
    !
    ! -- Part 2: below mixed layer --
    !
    !11) Calculate diffusivities, viscosities below the OBL
    !12) Ensure a background values
    !13) Convection below the OBL

    !-- Part 3: updating the mixing coefficients
    !
    !14) Update viscosity and diffusivities

    !========================================================================
    ! Renaming stuff  
    !======================================================================== 
    ! grid infos
    dz  => patch_3d%p_patch_1d(1)%prism_center_dist_c      ! depth of centres
    dzi => patch_3d%p_patch_1d(1)%inv_prism_center_dist_c  ! 1/dz
    dzw => patch_3d%p_patch_1d(1)%prism_thick_c            ! cell thicknesses
    kbot => patch_3d%p_patch_1d(1)%dolic_c(:,:)            ! last water cell
    patch_2D => patch_3d%p_patch_2d(1)
    cells_in_domain => patch_2D%cells%in_domain

    ! variables (maybe move into subroutines where it is used?)
    temp => ocean_state%p_prog(nold(1))%tracer(:,:,:,1)
    salt => ocean_state%p_prog(nold(1))%tracer(:,:,:,2)
    Av_old => params_oce%a_veloc_v(:,:,:)
    kv_old => params_oce%a_tracer_v(:,:,:,1)

    ! necessary for mlotst diagnostic
    vert_density_grad => ocean_state%p_diag%zgrad_rho

    ! KPP cvmix variables
    avo_kpp => params_oce%cvmix_params%avo_kpp(:,:,:)
    dvo_heat_kpp => params_oce%cvmix_params%dvo_heat_kpp(:,:,:)
    dvo_salt_kpp => params_oce%cvmix_params%dvo_salt_kpp(:,:,:)

    cellHeight_KPP => params_oce%cvmix_params%cellHeight_KPP(:,:,:)
    iFaceHeight_KPP => params_oce%cvmix_params%iFaceHeight_KPP(:,:,:)
    Tsurf => params_oce%cvmix_params%Tsurf(:,:,:)
    Ssurf => params_oce%cvmix_params%Ssurf(:,:,:)
    Usurf => params_oce%cvmix_params%Usurf(:,:)
    Uz2 => params_oce%cvmix_params%Uz2(:,:,:)
    Coriolis => params_oce%cvmix_params%Coriolis(:,:)
    Coriolis = patch_2D%cells%f_c
    SLdepth_2D => params_oce%cvmix_params%SLdepth_2D(:,:)
    BuoyFlux_3D => params_oce%cvmix_params%BuoyFlux_3D(:,:,:)
    stab_KPP => params_oce%cvmix_params%stab_KPP(:,:,:)
    dRho => params_oce%cvmix_params%dRho(:,:,:)
    N2 => params_oce%cvmix_params%N2(:,:,:)
    Nbuoy => params_oce%cvmix_params%Nbuoy(:,:,:)
    richardson_no => params_oce%cvmix_params%richardson_no(:,:,:)
    ustar => params_oce%cvmix_params%ustar(:,:)
    BulkRi => params_oce%cvmix_params%BulkRi(:,:,:)
    WS_cntr =>  params_oce%cvmix_params%WS_cntr(:,:,:)
    nonLocalTransHeat => params_oce%cvmix_params%nonLocalTransHeat(:,:,:)
    nonLocalTransScalar => params_oce%cvmix_params%nonLocalTransScalar(:,:,:)
    OBLdepth => params_oce%cvmix_params%OBLdepth(:,:)
    kOBL_2D  => params_oce%cvmix_params%kOBL_2D(:,:)
    WS => params_oce%cvmix_params%WS(:,:,:)
    WM => params_oce%cvmix_params%WM(:,:,:)
    nl_trans_tend_heat => params_oce%cvmix_params%nl_trans_tend_heat(:,:,:)
    nl_trans_tend_salt => params_oce%cvmix_params%nl_trans_tend_salt(:,:,:)
    G_3D => params_oce%cvmix_params%G_3D(:,:,:)
    Vt2 => params_oce%cvmix_params%Vt2(:,:,:)
    net_heat => params_oce%cvmix_params%net_heat(:,:)
    net_salt => params_oce%cvmix_params%net_salt(:,:)

    !========================================================================
    ! Prepare the atmospheric surface fluxes
    !========================================================================
    ! Qh: net heat flux (Wm-2 -> K m s-1): 1/(rho0*cp) * Qh
    net_heat = p_oce_sfc%HeatFlux_Total   ! Wm-2

    ! use this if coupled run?
    !net_heat => atmos_fluxes%HeatFlux_Total
    
    ! change unit of net heat flux
    net_heat = net_heat / (OceanReferenceDensity * clw) !K m s-1
   
    ! Qs: freshwater flux (kg m-2 s-1 -> psu m s-1)
    net_salt = p_oce_sfc%FrshFlux_TotalOcean         !psu m s-1

!------------------------------------------------------------------------
! PART 1: within mixed layer
!------------------------------------------------------------------------

    !========================================================================
    ! 1) Calculate the surface layer depth (SLdepth_2D), averaged surface layer quantities
    !     (Tsurf,Ssurf,U-sl), and square of bulk difference in resolved velocity (Uz2)
    !========================================================================
    ! calculates: SLdepth_2D
    !========================================================================
    ! FIXME: make Tsurf,Ssurf 2D variables?
    CALL calc_surface_layer_averages(patch_3d, ocean_state, params_oce,   &
               cellHeight_KPP, iFaceHeight_KPP, Coriolis, temp, salt,     &
               Tsurf, Ssurf, Usurf, Uz2, SLdepth_2D)

    !========================================================================
    ! 2) Calculate buoyancy flux (m2/s3) for WS in 5)
    !========================================================================
    ! calculates: BuoyFlux_3D
    !========================================================================
    CALL calc_buoyancy_flux(patch_3d, ocean_state, params_oce, temp, salt, &
                            net_heat, net_salt, BuoyFlux_3D)

    !========================================================================
    ! 3) compute in-situ density and buoyancy frequency
    !========================================================================
    ! calculates: dRho, richardson_no, Nbuoy, N2, stab_KPP
    !========================================================================
    CALL calc_insitu_density_and_N2(patch_3d, ocean_state, params_oce,    &
                temp, salt, Tsurf, Ssurf, dRho, richardson_no, Nbuoy, N2, &
                stab_KPP )

      ! update ocean state density gradient (necessary for mlotst diagnostic?
      vert_density_grad = stab_KPP 
   
    !========================================================================
    ! 4) Calculate friction velocity (ustar) at surface (m/s) --> WS in 5)
    !========================================================================
    ! calculates: ustar
    !========================================================================
    CALL calc_ustar(patch_3d, params_oce, atmos_fluxes, concsum, ustar)

    !========================================================================
    ! 5) Calculate the turbulent velocity scales w_s and w_m (WS) 
    !     at the cell centers
    !========================================================================
    ! calculates: WS_cntr
    !========================================================================
    CALL calc_turb_velocity_scales(patch_3d, params_oce,                  &
                cellHeight_KPP, BuoyFlux_3D, ustar, WS_cntr)

    !========================================================================
    ! 6) Calculate Bulk Richardson number (BulkRi)
    !========================================================================
    ! calculates: BulkRi
    !========================================================================
    CALL calc_bulk_Richardson_number(patch_3d, params_oce,                &
                cellHeight_KPP, dRho, Uz2, WS_cntr, Nbuoy, BulkRi)

    !========================================================================
    ! 7) Compute OBL depth (h) (units: m)
    !========================================================================
    ! calculates: OBLdepth, kOBL_2D
    !========================================================================
    CALL calc_OBLdepth(patch_3d, params_oce,                              &
               cellHeight_KPP, iFaceHeight_KPP, BulkRi, BuoyFlux_3D,      &
               ustar, Coriolis, OBLdepth, kOBL_2D)

    !========================================================================
    ! 8) Call CVMix/KPP to obtain OBL diffusivities (k), OBL depth (h), 
    ! and non-dim. shape function (G), and non-local transport: nlt=cs2*G
    ! cs2=6.32739901508
    ! note: main calculation, calculates also the turbulent vertical velocity
    !       scales and the non-local transport terms, assuming 
    !       NLT_shape="CVMIX". The non-local transports are not applied in this
    !       routine!
    !========================================================================
    ! calculates: avo_kpp, dvo_heat_kpp, dvo_salt_kpp, WM, WS, 
    !             nonLocalTransHeat, nonLocalTransScalar
    ! Note: it already updates the diffusivities in params_oce%a_veloc_v and
    !       params_oce%a_tracer_v   

    !FIXME: debugging fields
    !CALL check_minmaxmean_2D('u*',ustar(:,:),cells_in_domain,str_module,idt_src)
    !CALL check_minmaxmean_2D('SL_Bf',BuoyFlux_3D(:,1,:),cells_in_domain,str_module,idt_src)
    !CALL check_minmaxmean_2D('OBL',OBLdepth(:,:),cells_in_domain,str_module,idt_src)
    !CALL check_minmaxmean_2D('SLdepth',SLdepth_2d(:,:),cells_in_domain,str_module,idt_src)

    !CALL check_minmaxmean_3D('c_depth',cellHeight_KPP(:,:,:),cells_in_domain,str_module,idt_src)
    !CALL check_minmaxmean_3D('i_depth',iFaceHeight_KPP(:,:,:),cells_in_domain,str_module,idt_src)
    !CALL check_minmaxmean_3D('Bf',BuoyFlux_3D(:,:,:),cells_in_domain,str_module,idt_src)

    CALL calc_diffusivities_within_OBL(patch_3d, params_oce,              &
                BuoyFlux_3D, cellHeight_KPP, iFaceHeight_KPP, OBLdepth,   &
                kOBL_2D, ustar, avo_kpp, dvo_heat_kpp, dvo_salt_kpp,      &
                WM, WS, nonLocalTransHeat, nonLocalTransScalar)

    !========================================================================
    ! 9) (Re-) Calculate nonlocal transport terms and non-local tendency terms
    !========================================================================
    ! Adjusts the non-local transport terms if  NLT_shape!="CVMIX" and updates
    ! them, otherwise nothing will happen in calc_nonlocal_trans_heat_salt()
    !
    ! ! calculate the non-local transports (T,S): nlt=cs2*G
    ! note1: nlt is just the factor of the non-local transport terms 
    !        (nlt*Qh and nlt*Qs)
    ! note2: here G might be exchanged by user option
    !
    ! KPP_NonLocalTransport_temp() and KPP_NonLocalTransport_saln() calculate  
    ! then the non-local transport tendencies for temperature and salinity that
    ! are added to the tracer diffusion equation (mo_ocean_tracer.f90).
    !========================================================================
    ! calculates: nonLocalTransHeat, nonLocalTransScalar
    !========================================================================    
    CALL calc_nonlocal_trans_heat_salt(patch_3d, params_oce,               &
                iFaceHeight_KPP, BuoyFlux_3D, OBLdepth, nonLocalTransHeat, &
                nonLocalTransScalar)

    ! calculate the non-local tendencies (T,S): 
    !   dT/dt = ( nlt(k-1) - nlt(k) ) / dzw(k) * Qh
    !   dS/dt = ( nlt(k-1) - nlt(k) ) / dzw(k) * Qs
    ! NOTE: the tendencies have to be applied outside the KPP scheme!
    CALL KPP_NonLocalTransport_temp(patch_3d, params_oce, concsum,         &
                    nonLocalTransHeat, net_heat, nl_trans_tend_heat)

    CALL KPP_NonLocalTransport_saln(patch_3d, params_oce, concsum,         &
                    nonLocalTransScalar, net_salt, nl_trans_tend_salt)


    !========================================================================
    ! 10) Diagnostics for KPP
    !========================================================================

    ! vertical turbulent velocity scale (scalars)
    ! calculates: WM, WS
    IF (diag_WS) THEN 
      CALL calc_diag_turb_velocity_scales(patch_3d, BuoyFlux_3D,           &
                 iFaceHeight_KPP, OBLdepth, ustar, WM, WS)
    ENDIF

    ! non-dimensional shape function G=sigma(1-sigma)**2
    ! calculates: G_3D
    IF (diag_G) THEN
      CALL diagnose_G(patch_3d, params_oce,                                &
                 iFaceHeight_KPP, OBLdepth, kOBL_2D, G_3D)
    ENDIF

    ! unresolved turbulent vertical shear
    ! calculates: Vt2
    IF (diag_vtsqr) THEN
      CALL diagnose_vtsqr(patch_3d, params_oce,                            &
                 cellHeight_KPP, WS_cntr, Nbuoy, Vt2)
    ENDIF

!------------------------------------------------------------------------ 
! PART 2: below mixed layer
!------------------------------------------------------------------------

    !========================================================================
    ! 11) Calculate diffusivities, viscosities below the OBL
    !========================================================================
    ! calculates: avo_kpp, dvo_heat_kpp and dvo_salt_kpp
    ! Note: there is no distinguishing of tracers below the OBL, so that
    !       dvo_heat_kpp = dvo_salt_kpp
    ! Note2: it already updates the diffusivities in params_oce%a_veloc_v !
    !        and params_oce%a_tracer_v
    !========================================================================
    CALL calc_diffusivities_below_OBL(patch_3d, params_oce, richardson_no, &
             kOBL_2D, avo_kpp, dvo_heat_kpp, dvo_salt_kpp)
 
    !========================================================================
    ! 12) ensure a background value for avo/dvo
    !========================================================================
   
    CALL ensure_background_values(patch_3d, params_oce, & 
            avo_kpp, dvo_heat_kpp, dvo_salt_kpp)

    !========================================================================
    ! 13) Apply convection parameterization below the OBL
    !========================================================================

    IF (lconvection) THEN

     ! apply the enhanced diffusion parameterisation for convection.
     ! Note: only for diffusivities, viscosity is unaffected
     !       it updates dvo_heat_kpp and dvo_salt_kpp
     CALL calc_enhanced_vert_diff_below_OBL(patch_3d, params_oce,          &
                      stab_KPP, kOBL_2D, dvo_heat_kpp, dvo_salt_kpp)

    ENDIF

!------------------------------------------------------------------------
! PART 3: update the mixing coefficients
!------------------------------------------------------------------------

    !========================================================================
    ! 14) Update viscosity and diffusivities
    !========================================================================
    ! updates params_oce%a_veloc_v and params_oce%a_tracer_v
    ! FIXME: the latter only for tracers 1 (heat) and 2 (salinity) a.t.m.
    !========================================================================
    CALL update_mixing_coefficients(patch_3d, params_oce,                 &
                avo_kpp, dvo_heat_kpp, dvo_salt_kpp)



  END SUBROUTINE calc_kpp

!========================================================================
!========================================================================

  !-------------------------------------------------------------------------
  !>
  !! This routine calculates the averaged surface layer quantities T_SL, S_SL, U_SL, 
  !! and dU^2=(U-U_SL)^2, which is required to calculate the modified bulk
  !! Richardson number (RiB). The depth of the surface layer (SL) is defined as 
  !! surface_layer_ext*h. As default surface_layer_ext=0.1, so that the SL is 
  !! 10% of the ocean surface boundary layer (OBL), h.
  !! This routine follows the Column Sampling Method of the CVMix documentation,
  !! see eq. (8.173) and section 8.5.7.2 on page 73. It assumes that the SL
  !! depth equals the present cell depth.
  !!
  !!
  !! @par Revision History
  !! Initial release by Oliver Gutjahr, MPI-M (2019-04)
  !-------------------------------------------------------------------------
  SUBROUTINE calc_surface_layer_averages(patch_3d, ocean_state, params_oce, &
               cellHeight_KPP, iFaceHeight_KPP, Coriolis, temp, salt, Tsurf, Ssurf, &
               Usurf, Uz2, SLdepth_2D)
    TYPE(t_patch_3d ),TARGET, INTENT(in)  :: patch_3d
    TYPE(t_patch),POINTER                 :: patch_2D                        
    TYPE(t_subset_range), POINTER         :: all_cells !,edges_in_domain
    TYPE(t_hydro_ocean_state), TARGET     :: ocean_state
    TYPE(t_ho_params), INTENT(inout)      :: params_oce

    ! required variables
    REAL(wp),DIMENSION(:,:,:),INTENT(IN)    :: &
      temp                                    ,& ! temperature (K)
      salt                                       ! salinity (psu)

    REAL(wp),DIMENSION(:,:),INTENT(IN)      :: &
      Coriolis                                   ! Coriolis (s-1)

    REAL(wp),DIMENSION(:,:,:),INTENT(INOUT) :: &     
      cellHeight_KPP                          ,& ! absolute depth of cell centres (m)
      iFaceHeight_KPP                         ,& ! absolute depth of cell interfaces (m) 
      Tsurf                                   ,& ! surface-layer averaged temperature (K)
      Ssurf                                   ,& ! surface-layer averaged salinity (psu)
      Uz2                                        ! Square of bulk difference in resolved velocity: 
                                                 !   (U(k) - Usurf)**2 (m^2/s-2)

    REAL(wp),DIMENSION(:,:),INTENT(INOUT)   :: &
      Usurf                                   ,& ! surface-layer averaged velocity (m/s)
      SLdepth_2D                                 ! surface-layer depth (m) 

    ! local variables
    REAL(wp)                               ::  &
      SLdepth_0d                              ,& ! surface-layer depth (m)
      surfHtemp                               ,& ! dummy: level temperature (degC)
      surfHsalt                               ,& ! dummy: level salinity (psu)
      surfHu                                  ,& ! dummy: level velocity (m/s)
      hTot                                    ,& ! SL thickness (m)
      delH                                       ! 

    ! loop variables
    INTEGER :: jc, blockNo, jk, ksfc, ktmp
    INTEGER :: start_index, end_index
    INTEGER :: levels

    REAL(wp),POINTER                       :: &
      depth_CellMiddle(:,:,:)                ,&  ! depth of cell centres (m)
      depth_CellInterface(:,:,:)                 ! depth of cell interfaces (m)

    ! grid information
    patch_2D => patch_3d%p_patch_2d(1) 
    all_cells => patch_3d%p_patch_2d(1)%cells%ALL
    depth_CellMiddle => patch_3d%p_patch_1d(1)%depth_CellMiddle
    depth_CellInterface => patch_3d%p_patch_1d(1)%depth_CellInterface
    dzw => patch_3d%p_patch_1d(1)%prism_thick_c
    kbot => patch_3d%p_patch_1d(1)%dolic_c(:,:) 
    
    !========================================================================
    ! 1) Calculate the surface layer depth, averaged surface layer quantities,
    !     and non-local shear of vertical velocity
    !
    !    cellHeight_KPP  --> cell center depths (n_zlev)
    !    iFaceHeight_KPP --> interface depths   (n_zlev+1) 
    !    Uz2(i,j,k) --> square of bulk difference in resolved velocity (m2/s2) required for Bulk_Ri_3D in 7)
    !    Tsurf(i,j) --> required for dRho in 4)
    !    Ssurf(i,j) --> required for dRho in 4)
    !    Usurf(i,j) -->
    !    SLdepth_2D(i,j) -> diagnostic: surface layer depth 
    !
    ! note: this routine follows the Column Sampling Method of the CVMix docu
    !       eq. (8.173) and section 8.5.7.2 on page 73. It assumes that the 
    !       surface layer depth equals the present cell depth.
    !========================================================================
 
    ! initialise surface layer depth variable
    SLdepth_2D(:,:) = 0.0_wp

    DO blockNo = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blockNo, start_index, end_index)
      DO jc = start_index, end_index
        !levels       = patch_3d%p_patch_1d(1)%dolic_e(jc, blockNo)   !at edge 
        !levels       = patch_3d%p_patch_1d(1)%dolic_c(jc, blockNo)   !at centre
        !iFaceHeight_KPP(jc,1,blockNo) = 0.0_wp

        ! only for ocean cells
        IF (kbot(jc,blockNo)>0) THEN

        DO jk=1,n_zlev
        !DO jk=1,levels  !levels seem to cause an error later?
        !DO jk=1,kbot(jc,blockNo) 
          !========================================================================
          ! 1.1) depth below surface (cell centre and interface points) in
          ! meters. Negative in ocean.
          !========================================================================
          ! patch3d wet_c --> similar to 'weto'
          ! prism_thick_c
          !
          ! Option A): displacement of sea surface height is considered
          ! note: depth_CellMiddle and depth_CellInterface vary slightly over
          ! the spatial domain. This is because of the the sea surface height
          ! change.
          ! FIXME: not working correctly: causes wrong OBL depths...
!          cellHeight_KPP(jc,jk,blockNo)  = -depth_CellMiddle(jc,jk,blockNo)        ! depth at tracer point (cell middle)
!          iFaceHeight_KPP(jc,jk,blockNo) = -depth_CellInterface(jc,jk,blockNo)   ! depth at interface point (cell interface)

          !Option B): displacement of sea suface height not considered 
          cellHeight_KPP(jc,jk,blockNo)  = -1.0_wp * patch_3d%p_patch_1d(1)%zlev_m(jk)
          iFaceHeight_KPP(jc,jk,blockNo) = -1.0_wp * patch_3d%p_patch_1d(1)%zlev_i(jk)

          !========================================================================
          ! 1.2) surface layer depth (m) and associated model level
          !========================================================================
          ! get surface layer depth (z-axis from top to bottom!) 
          ! Assumption is 10% of total water depth, but at least as deep as
          ! first interface level below surface.
        
          ! we simply assume 10% of total water depth
          SLdepth_0d = surf_layer_ext * MAX( MAX( -cellHeight_KPP(jc,jk,blockNo), -iFaceHeight_KPP(jc,2,blockNo) ), minOBLdepth )
          ! debugging: fixed SL depth
          !SLdepth_0d = -1.0_wp * surf_layer_ext * iFaceHeight_KPP(jc,5,blockNo)

          !FIXME: if 1.1A) is used, limit water column to depth of kbot in case of 1.1A)
!          IF (jk .GE. kbot(jc,blockNo)) THEN
!            SLdepth_0d = surf_layer_ext * MAX( MAX(-cellHeight_KPP(jc,kbot(jc,blockNo),blockNo), -iFaceHeight_KPP(jc,2,blockNo) ), minOBLdepth )
!          END IF

          ! find model layer (ksfc) for cell where "surface layer" sits,
          ! corresponds to Eq.(8.173) of the CVMix documentation
          ksfc = jk
          DO ktmp = 1,jk
            IF (-1.0_wp*iFaceHeight_KPP(jc,ktmp+1,blockNo) >= SLdepth_0d) THEN
              ksfc = ktmp
              EXIT
            ENDIF
          ENDDO !ktmp

          !========================================================================
          ! 1.3) Average temperature, salinity, and velocity over surface layer
          !========================================================================
          ! initialise values
          surfHtemp = 0.0_wp
          surfHsalt = 0.0_wp
          surfHu    = 0.0_wp
          !hTot      = 0.0_wp
          hTot      = dbl_eps

          DO ktmp = 1,ksfc
            ! SLdepth can be between cell interfaces
            delH = MIN( MAX(0.0_wp, SLdepth_0D - hTot), dzw(jc,ktmp,blockNo) )
            ! surface layer thickness
            hTot = hTot + delH

            ! surface layer averaged fields,
            ! individual layers are weighted by the layer thickness
            surfHtemp = surfHtemp + temp(jc,ktmp,blockNo) * delH 
            surfHsalt = surfHsalt + salt(jc,ktmp,blockNo) * delH 

            ! ocean_state%p_diag%p_vn%x has dimension x(3), therefore use SUM()
            ! p_vn in centre of box and %x has 3 kartesian components 
            ! (vn would be on edges)
            surfHu = SUM(ocean_state%p_diag%p_vn(jc,ktmp,blockNo)%x) * delH
 
          ENDDO ! ktmp

          ! Surface layer averages (T, S, and velocity)
          ! FIXME: Tsurf and Ssurf are 3d-fields here, may change to 2d field
          ! later
          Tsurf(jc,jk,blockNo) = surfHTemp / hTot
          Ssurf(jc,jk,blockNo) = surfHSalt / hTot
          Usurf(jc,blockNo) = surfHu / hTot

          !------------------------------------------------------------------------
          ! Square of bulk difference in resolved velocity
          ! [Uz2] = m2s-2
          !------------------------------------------------------------------------
          Uz2(jc,jk,blockNo) =  dbl_eps +  SUM( (ocean_state%p_diag%p_vn(jc,jk,blockNo)%x &
                                                - Usurf(jc,blockNo) )**2)

        ENDDO  !jk 

        !========================================================================
        ! 1.4) store surface layer depth (diagnostic)
        !========================================================================
 
        SLdepth_2D(jc,blockNo) = SLdepth_0d

       END IF !ocean cell

      ENDDO
    ENDDO

  END SUBROUTINE calc_surface_layer_averages

!========================================================================
!========================================================================
  ! This routine calculates the buoyancy flux as Bf(k) = g / rho0 * ( alpha * Qheat(k)/cp + beta(k)*Qsalt) (m2/s3)
  ! note: cp is already contained in 'net_heat'
  ! Bf is required to calculate the modified bulk Richardson number (RiB). The RiB is crucial
  ! in defining the depth of the ocean surface boundary layer (OBL), where RiB>Ri_crit.

  SUBROUTINE calc_buoyancy_flux(patch_3d, ocean_state, params_oce, temp, salt, &
                                net_heat, net_salt, BuoyFlux_3D)
    TYPE(t_patch_3d ),TARGET, INTENT(in)  :: patch_3d
    TYPE(t_subset_range), POINTER         :: all_cells
    TYPE(t_hydro_ocean_state), TARGET     :: ocean_state
    TYPE(t_ho_params), INTENT(inout)      :: params_oce

    ! loop variables
    INTEGER :: jc, blockNo, jk
    INTEGER :: start_index, end_index
    INTEGER :: levels

    ! required variables
    REAL(wp),DIMENSION(:,:),INTENT(IN)    ::   &
      net_heat                                ,& ! net heat flux at surface (Km s-1)
      net_salt                                   ! net freshwater flux at surface (ppm)

    REAL(wp),DIMENSION(:,:,:),INTENT(IN)    :: &
      temp                                    ,& ! temperature (K)
      salt                                       ! salinity (psu)

    REAL(wp),DIMENSION(:,:,:),INTENT(INOUT) :: &
      BuoyFlux_3D                                ! buoyancy flux (m3 s-2)


    REAL(wp) :: z_grav_rho
    REAL(wp) :: neutral_coeff(1:n_zlev, 2)!, salinityColumn(1:n_zlev)


    !========================================================================
    ! 2) Calculate buoyancy flux (m2/s3) for WS in 6) 
    !========================================================================
    ! alpha = neutral_coeff(:,1) ! K-1
    ! beta  = neutral_coeff(:,2) ! psu-1
    !-------------------------------------------------------------------------
    z_grav_rho = grav/OceanReferenceDensity
    all_cells => patch_3d%p_patch_2d(1)%cells%ALL
    kbot => patch_3d%p_patch_1d(1)%dolic_c(:,:)
 
    ! FIXME: in case of absent salinity tracer (necessary?)
    !salinityColumn(1:n_zlev) = sal_ref

    DO blockNo = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blockNo, start_index, end_index)
      DO jc = start_index, end_index
        levels       = patch_3d%p_patch_1d(1)%dolic_e(jc, blockNo)   ! 

        ! only for ocean cells
        IF (kbot(jc,blockNo)>0) THEN

        ! calculate alpha and beta
        IF(EOS_TYPE/=1)THEN
          !Nonlinear EOS, slope coefficients are calculated via the McDougall-method
          ! as in FESOM
          ! FIXME: use salinity instead of salinityColumn?
          neutral_coeff = calc_neutralslope_coeff_func_onColumn(          &
              !& temp(jc,1:levels,blockNo), salinityColumn(1:levels),  &
              & temp(jc,1:levels,blockNo), salt(jc,1:levels,blockNo), &
              !& depth_cellinterface(jc,2:levels+1,blockNo),levels)
              & patch_3D%p_patch_1d(1)%depth_cellinterface(jc,2:levels+1,blockNo),levels) 

       ELSEIF(EOS_TYPE==1)THEN
          !Linear EOS: slope coefficients are equal to EOS-coefficients
          neutral_coeff(:,1) = LinearThermoExpansionCoefficient
          neutral_coeff(:,2) = LinearHalineContractionCoefficient
       ENDIF


       ! loop over depth levels
       DO jk=1,levels

          ! bouyancy_flux = g/rho0 * (alpha*Qheat/cp + beta*Qsalt) (m2/s3)
          ! -> FIXME: is something like ddpo needed?
          ! heat_flux & water_flux: positive up (?)

          !FIXME: check whether sign is correct, or use -z_grav_rho * ... ? 
          ! heat/fw flux negative up (?)
          ! net_heat is in unit K m s-1
          ! net_salt is in unit psu m s-1
          ! alpha is in unit K-1
          ! beta is in unit psu-1
          ! BuoyFlux_3D is in unit m2 s-3
          !---> net_heat was already multiplied by g/(rho0*Cp)
          !BuoyFlux_3D(jc,jk,blockNo) = z_grav_rho * ( neutral_coeff(jk,1)* net_heat(jc,blockNo)/clw   &           
           BuoyFlux_3D(jc,jk,blockNo) = grav * ( neutral_coeff(jk,1) * net_heat(jc,blockNo) &
                                                   +  neutral_coeff(jk,2)* net_salt(jc,blockNo) )

          ! FIXME: add shortwave radiation absorption to buoyancy flux?
          ! 
          ! ...

       ENDDO !jk

       END IF ! ocean only
      ENDDO !jc
    ENDDO !blockNo


  END SUBROUTINE calc_buoyancy_flux

!========================================================================
!========================================================================
  ! This routine calculates the insitu density gradient (stab_KPP), N2,
  ! and the local Richardson number (Ri), which is required to define convective conditions
  ! below the ocean surface boundary layer (OBL). 

  SUBROUTINE calc_insitu_density_and_N2(patch_3d, ocean_state, params_oce, &
                temp, salt, Tsurf, Ssurf, dRho, richardson_no, Nbuoy, N2, stab_KPP )
    TYPE(t_patch_3d ),TARGET, INTENT(in)  :: patch_3d
    TYPE(t_subset_range), POINTER         :: all_cells
    TYPE(t_hydro_ocean_state), TARGET     :: ocean_state
    TYPE(t_ho_params), INTENT(inout)      :: params_oce

    ! loop variables
    INTEGER :: jc, blockNo, jk
    INTEGER :: start_index, end_index
    INTEGER :: levels

    ! required variables
    REAL(wp),DIMENSION(:,:,:),INTENT(IN)    :: &
      temp                                    ,& ! temperature (K)
      salt                                    ,& ! salinity (psu)
      Tsurf                                   ,& ! surface-layer averaged temperature (K)
      Ssurf                                      ! surface-layer averaged salinity (psu)

    REAL(wp),DIMENSION(:,:,:),INTENT(INOUT) :: &
      stab_KPP                                ,& ! vertical stability ()
      dRho                                    ,& ! vertical density gradient ()
      Nbuoy                                   ,& ! buoyancy frequency () 
      N2                                      ,& ! squared buoyancy frequency ()
      richardson_no


    REAL(wp) :: z_grav_rho                       ! g/rho0 ()
    REAL(wp) :: vert_velocity_shear

    ! local variables - 1d
    REAL(wp) :: rho_up(n_zlev)
    REAL(wp) :: rho_down(n_zlev)
    REAL(wp) :: rho_sl(n_zlev)       
    REAL(wp) :: pressure(n_zlev)

    !========================================================================
    ! 3.1) vertical stability (non-local dRho for BulkRi_3D in 6))
    !========================================================================

    z_grav_rho = grav/OceanReferenceDensity 
    all_cells => patch_3d%p_patch_2d(1)%cells%ALL
    levels = n_zlev
    kbot => patch_3d%p_patch_1d(1)%dolic_c(:,:)

    ! initialise vectors/fields
    rho_up(:) = 0.0_wp
    rho_down(:) = 0.0_wp
    stab_KPP(:,:,:) = 0.0_wp
    dRho(:,:,:) = 0.0_wp
    N2(:,:,:) = 0.0_wp
    Nbuoy(:,:,:) = 0.0_wp

    DO blockNo = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blockNo, start_index, end_index)

      ! initialise Richardson number
      richardson_no(:,:,blockNo) = 0.0_wp

      DO jc = start_index, end_index
        !levels       = patch_3d%p_patch_1d(1)%dolic_e(jc, blockNo)   !

        ! only for ocean cells
        IF (kbot(jc,blockNo)>0) THEN

        pressure(1:levels) = patch_3d%p_patch_1d(1)%depth_CellInterface(jc, 1:levels, blockNo) * ReferencePressureIndbars 

        ! averaged surface layer density
        rho_sl(2:levels) = calculate_density_onColumn(&
            & Tsurf(jc,2:levels,blockNo), &
            & Ssurf(jc,2:levels,blockNo), &
            & pressure(2:levels), levels-1)   ! FIXME: correct pressure(2:levels)?

        ! at k-1
        rho_up(1:levels-1)  = calculate_density_onColumn(&
            & temp(jc,1:levels-1,blockNo), &
            & salt(jc,1:levels-1,blockNo), &
            & pressure(2:levels), levels-1)
        ! at k 
        rho_down(2:levels)  = calculate_density_onColumn(&
            & temp(jc,2:levels,blockNo), &
            & salt(jc,2:levels,blockNo), &
            & pressure(2:levels), levels-1)
  
        DO jk = 2, levels
          ! non-local rho difference: difference to surface layer averaged density
          dRho(jc,jk,blockNo) = rho_down(jk) - rho_sl(jk)

          !========================================================================
          ! 3.2) buoyancy frequency based on local stability (1:ke+1) for BulkRi_3D
          ! in 6)
          !========================================================================
          ! vertical stability
          stab_KPP(jc,jk,blockNo) = dRho(jc,jk,blockNo) * dzi(jc,jk,blockNo)

          ! vertical stability and buoyancy frequency
          N2(jc,jk,blockNo) = z_grav_rho*stab_KPP(jc,jk,blockNo)

          ! used for unresolved shear calculation
          Nbuoy(jc,jk,blockNo)  = SQRT( MAX( N2(jc,jk,blockNo), 0.0_wp) )

          ! division by dz**2 is omitted in this calculation of velocity shear: shear = (d_vn)**2
          vert_velocity_shear = dbl_eps + &
              & SUM((ocean_state%p_diag%p_vn(jc,jk-1,blockNo)%x - ocean_state%p_diag%p_vn(jc,jk,blockNo)%x)**2)

          ! Richardson-number
          ! Ri = g/OceanReferenceDensity * dz * d_rho/(d_vn)**2
          ! note: Ri is not restricted to positive values! Important for mixing
          ! below the mixed layer (internal_mix_scheme ='KPP')
          richardson_no(jc,jk,blockNo) = patch_3d%p_patch_1d(1)%prism_center_dist_c(jc,jk,blockNo) * z_grav_rho * &
              &                           (rho_down(jk) - rho_up(jk-1)) / vert_velocity_shear

        ENDDO ! levels

        END IF ! ocean only

      ENDDO !jc
    ENDDO !blockNo

  END SUBROUTINE calc_insitu_density_and_N2

!========================================================================
!========================================================================
  ! This routine calculates the surface friction velocity u* = sqrt(tau / rho0).
  ! u* is the main quantity that defines the turbulent vertical velocity scales.
 
  SUBROUTINE calc_ustar(patch_3d, params_oce, atmos_fluxes, concsum, ustar)
    TYPE(t_patch_3d ),TARGET, INTENT(IN)  :: patch_3d
    TYPE(t_subset_range), POINTER         :: all_cells
    TYPE(t_atmos_fluxes)                  :: atmos_fluxes
    TYPE(t_ho_params), INTENT(INOUT)      :: params_oce
    REAL(wp), TARGET                      :: concsum(:,:) ! t_sea_ice%concsum

    ! required variables
    REAL(wp), DIMENSION(:,:), INTENT(INOUT) :: ustar  ! surface friction velocity (m s-1) 

    ! loop variables
    INTEGER :: jc, blockNo
    INTEGER :: start_index, end_index

    REAL(wp), DIMENSION(nproma,patch_3d%p_patch_2d(1)%alloc_cell_blocks) ::    &
      tau_abs                                    ! wind stress at ocean surface (m s-1)

    ! calculation
    all_cells => patch_3d%p_patch_2d(1)%cells%ALL 
    kbot => patch_3d%p_patch_1d(1)%dolic_c(:,:)

    DO blockNo = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blockNo, start_index, end_index)
      DO jc = start_index, end_index

        ! only for ocean cells
        !IF (kbot(jc,blockNo)>0) THEN

        ! wind_stress at ocean surface (m/s): u* = sqrt( tau/rho0 )
        ! note: stress_xw and stress_yw are only defined for ocean grid cells?
        tau_abs(jc,blockNo) = SQRT( (atmos_fluxes%stress_xw(jc,blockNo))**2     &
                                  + (atmos_fluxes%stress_yw(jc,blockNo))**2  )

        IF (use_reduced_mixing_under_ice) THEN
         ! reduce ustar below sea ice
         ustar(jc,blockNo) = (1.0_wp - concsum(jc,blockNo))**2 * &
                             SQRT( tau_abs(jc,blockNo) / OceanReferenceDensity )
        ELSE
         ! wind stress ignores presence of sea ice
         ustar(jc,blockNo) = SQRT( tau_abs(jc,blockNo) / OceanReferenceDensity )
        ENDIF

        !END IF !ocean only

      ENDDO !jc
    ENDDO !blockNo
 
  END SUBROUTINE calc_ustar

!========================================================================
!========================================================================
  ! This routine calculates the turbulent vertical velocity scales for momentum
  ! , w_m, and for scalars/tracers, w_s. These velocity scales are required to
  ! calculate the viscosity (K=h*w_m*G) and the diffusivities (K=h*w_s*G).
  ! Optionally, if ' =.TRUE.' the scales are modified by including a 
  ! Langmuir factor 'langmuir_Efactor'that represents Langmuir turbulence. 
  ! The modified scales are then: w_m=w_m*langmuir_Efactor and w_s=w_s*langmuir_Efactor. 
  ! However, it was not tested yet so that 'langmuir_Efactor=1.0' as default.
  ! The turbulent velocity scales mainly depend on the surface friction velocity
  ! u* and on the surface buoyancy forcing (Bf). The depth of the surface layer
  ! (surface_layer_ext*h) further defines the depth at which the profiles of w_m and w_s become
  ! constant with depth. As default 'surface_layer_ext=0.1'.

  SUBROUTINE calc_turb_velocity_scales(patch_3d, params_oce, &
                cellHeight_KPP, BuoyFlux_3D, ustar, WS_cntr)
    TYPE(t_patch_3d ),TARGET, INTENT(in)  :: patch_3d
    TYPE(t_subset_range), POINTER         :: all_cells
    TYPE(t_ho_params), INTENT(inout)      :: params_oce

    ! required variables
    REAL(wp), DIMENSION(:,:,:), INTENT(IN) ::      &
      cellHeight_KPP                              ,& ! depth of cell centres (m)
      BuoyFlux_3D                                    ! buoyancy flux (m2 s-3)
 
    REAL(wp), DIMENSION(:,:), INTENT(IN) ::        &
      ustar                                          ! surface friction velocity (m/s)

    REAL(wp), DIMENSION(:,:,:), INTENT(INOUT) ::   &
      WS_cntr                                        ! vertical turb. velocity scale for momentum

    ! loop variables
    INTEGER :: jc, blockNo
    INTEGER :: start_index, end_index

    ! parameters from namelist:
    !  surf_layer_ext
    !  langmuir_Efactor

    !========================================================================
    ! 5) Calculate the turbulent velocity scales w_s and w_m (WS) at the cell
    ! centers
    !========================================================================
    ! Note that if sigma > surf_layer_ext, then CVmix_kpp_compute_turbulent_scales
    ! computes w_s and w_m velocity scale at sigma=surf_layer_ext. So we
    ! only pass sigma=surf_layer_ext for this calculation.
    ! NOTE: if llangmuirEF=.TRUE.:  w_m = w_m * langmuir_Efactor
    !                               w_s = w_s * langmuir_Efactor
    ! surface_layer_ext = Fraction of OBL depth considered in the surface layer
    ! (nondim)
    ! SL depth is then calculated by surf_layer_ext * -cellHeight

    all_cells => patch_3d%p_patch_2d(1)%cells%ALL
    kbot => patch_3d%p_patch_1d(1)%dolic_c(:,:)

    DO blockNo = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blockNo, start_index, end_index)
      DO jc = start_index, end_index

        ! only for ocean cells
        IF (kbot(jc,blockNo)>0) THEN

        ! output from subroutine is: WS_cntr
        CALL cvmix_kpp_compute_turbulent_scales(                  &
                 sigma_coord     = surf_layer_ext,                & ! (in)  sigma: Normalized surface layer depth
                 OBL_depth       = -cellHeight_KPP(jc,:,blockNo), & ! (in)  Assume OBL depth (m) = -cellHeight(k)
                 surf_buoy_force = BuoyFlux_3D(jc,:,blockNo),     & ! (in)  Buoyancy flux at surface (m2/s3) (note: 1d-vector)
                 surf_fric_vel   = ustar(jc,blockNo),             & ! (in)  Turbulent friction velocity at surface (m/s)
                 langmuir_Efactor= langmuir_Efactor,              & ! (in)  Enhancement factor due to Langmuir circulation
                 w_s             = WS_cntr(jc,:,blockNo) ) ! ,         & ! (out) Turbulent velocity scale profile (m/s)
                 !CVmix_kpp_params_user = KPP_params               & ! (in) KPP parameter list
             !)

        END IF ! ocean only
      ENDDO
    ENDDO

  END SUBROUTINE calc_turb_velocity_scales

!========================================================================
!========================================================================
  ! This routine calculates the modified bulk Richardson-number (RiB(k)) for 
  ! each model layer k. In the definition of this modified RiB, all quantities are
  ! expressed as local differences to the surface layer averaged quantities, e.g. 
  ! dT(k) = T(k) - T_SL. In the denominator, the shear term dU(k)^2=(U(k)-U_SL is)^2
  ! extended by the unresolved turbulent vertical shear (Vt_sqr). This term is
  ! optional input, but currently disabled. In that case, CVMix calculates
  ! Vt_sqr internally. Vt_sqr can be diagnosed with the routine 'diagnose_vtsqr,
  !' if the namelist 'diag_vtsqr'=".TRUE.". In addition, RiB is also optionally 
  ! modified by the surface Stokes drift velocity, which can be provided as input 
  ! to the routine, but is disabled a.t.m.

  SUBROUTINE calc_bulk_Richardson_number(patch_3d, params_oce, &
                cellHeight_KPP, dRho, Uz2, WS_cntr, Nbuoy, BulkRi)
    TYPE(t_patch_3d ),TARGET, INTENT(in)  :: patch_3d
    TYPE(t_subset_range), POINTER         :: all_cells
    TYPE(t_ho_params), INTENT(inout)      :: params_oce

    ! loop variables
    INTEGER :: jc, blockNo 
    INTEGER :: start_index, end_index

    ! required variables
    REAL(wp), DIMENSION(:,:,:), INTENT(IN) ::         &
      cellHeight_KPP                                 ,& ! depth of cell centres (m)
      dRho                                           ,& ! vertical density gradient ()
      Uz2                                            ,& ! vertical velocity shear (U(k) - Usurf)
      WS_cntr                                        ,& ! vertical turb. velocity scale for scalars
      Nbuoy                                             ! buoyancy frequency ()

    REAL(wp), DIMENSION(:,:,:), INTENT(INOUT) ::      &
      BulkRi                                            ! Bulk-Richardson number 

    REAL(wp) :: z_grav_rho

    ! optional fields:
    !   stokes_drift(:,:)
    !   Vt_sqr_cntr(:,:,:) 
 
    !========================================================================
    ! 6) Calculate Bulk Richardson number (BulkRi_3D), CVMix docu p.73
    !========================================================================
    ! computed for each cell in a column,
    ! assuming OBLdepth = depth of present grid cell. After Rib(k) is
    ! known for the column, then CVMix interpolates to find
    ! the actual OBLdepth (RiB>Ri_crit). This approach avoids need to iterate
    ! on the OBLdepth calculation. It follows that used in MOM6
    !------------------------------------------------------------------------
    ! NOTE: stokes_drift (optional) is not used at the moment
    !       Vt_sqr_cntr is not specified, cvmix calculates it internally

    z_grav_rho = grav/OceanReferenceDensity
    all_cells => patch_3d%p_patch_2d(1)%cells%ALL
    kbot => patch_3d%p_patch_1d(1)%dolic_c(:,:)

    DO blockNo = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blockNo, start_index, end_index)
      DO jc = start_index, end_index
        !levels = patch_3d%p_patch_1d(1)%dolic_e(jc, blockNo) 
        ! only for ocean cells
        IF (kbot(jc,blockNo)>0) THEN

        ! Bulk-Richardson number at cell center
        BulkRi(jc,:,blockNo) = cvmix_kpp_compute_bulk_Richardson(          &
                    zt_cntr         = cellHeight_KPP(jc,:,blockNo),        & ! (in) Depth of cell center (m)
                    delta_buoy_cntr = z_grav_rho*dRho(jc,:,blockNo),       & ! (in) Bulk buoyancy difference, Br-B(z) (1/s)
                    delta_Vsqr_cntr = Uz2(jc,:,blockNo),                   & ! (in) Square of resolved velocity difference (m2/s2)
                    ws_cntr         = WS_cntr(jc,:,blockNo),               & ! (in) Turbulent velocity scale profile (m/s)
                    !stokes_drift   = stokes_drift(jc,blockNo),            & ! (in, optional) surface Stokes drift velocity (m/s)
                    !Vt_sqr_cntr    = Vt_sqr_cntr(jc,:,blockNo),           & ! (in, optional) squared unresolved shear term (m2/s2)
                    N_iface         = Nbuoy(jc,:,blockNo) )!,              & ! (in) Buoyancy frequency on interfaces (n_zlev+1) (1/s)
                    !CVmix_kpp_params_user = KPP_params)                      ! (in) KPP parameter list                        
        END IF ! ocean only
      ENDDO
    ENDDO


  END SUBROUTINE calc_bulk_Richardson_number

!========================================================================
!========================================================================
  ! Computes the ocean surface boundary layer (OBL) depth, h, which is the main quantity
  ! to control the behaviour of the KPP scheme. The OBL depth is required to
  ! calculate K = h * w_s * G. The OBL depth is defined as the (linear
  ! interpolated) depth where the local modified bulk Ri exceeds the critical Ri
  ! number 'Ri_crit'.

  SUBROUTINE calc_OBLdepth(patch_3d, params_oce, &
               cellHeight_KPP, iFaceHeight_KPP, BulkRi, BuoyFlux_3D, ustar, Coriolis, OBLdepth, kOBL_2D)
    TYPE(t_patch_3d ),TARGET, INTENT(in)  :: patch_3d
    TYPE(t_patch),POINTER                 :: patch_2D
    TYPE(t_subset_range), POINTER         :: all_cells
    TYPE(t_ho_params), INTENT(inout)      :: params_oce

    ! loop variables
    INTEGER :: jc, blockNo
    INTEGER :: start_index, end_index

    ! required variables
    REAL(wp), DIMENSION(:,:,:), INTENT(IN) ::         &
      cellHeight_KPP                                 ,& ! depth of cell centres (m)
      iFaceHeight_KPP                                ,& ! depth of cell interfaces (m)
      BulkRi                                         ,& ! Bulk-Richardson number
      BuoyFlux_3D                                       ! Buoyancy flux (m2 s-3)

    REAL(wp), DIMENSION(:,:), INTENT(IN) ::           &
      Coriolis                                       ,& ! Coriolis parameter (s-1)
      ustar                                             ! surface friction velocity (m/s)

    REAL(wp), DIMENSION(:,:), INTENT(INOUT) ::        &
      OBLdepth                                       ,& ! Ocean bounday layer depth (m)
      kOBL_2D                                           ! model layer at OBL depth (float)

    REAL(wp)                                    ::    &
      zBottomMinusOffset                                !

    !========================================================================
    ! 8) Compute OBL depth (m)
    !========================================================================
    ! This is only used in kpp_compute_OBL_depth to limit
    ! h to Monin-Obukov (default is false, ie. not used)
    ! note: if BulkRi_1d equals size of iFaceHeight, then the calculations use
    !       iFaceHeight, otherwise if it equals cellHeight, then the
    !       cellHeights are used.

    patch_2D => patch_3d%p_patch_2d(1) 
    all_cells => patch_3d%p_patch_2d(1)%cells%ALL 
    kbot => patch_3d%p_patch_1d(1)%dolic_c(:,:)

    DO blockNo = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blockNo, start_index, end_index)
      DO jc = start_index, end_index

       ! only for ocean cells
        IF (kbot(jc,blockNo)>0) THEN

        CALL cvmix_kpp_compute_OBL_depth(            &
               Ri_bulk    = BulkRi(jc,:,blockNo),           & ! (in) Bulk Richardson number dim=(n_zlev+1)
               zw_iface   = iFaceHeight_KPP(jc,:,blockNo),  & ! (in) Height of interfaces (m) dim=(n_zlev+1)
               OBL_depth  = OBLdepth(jc,blockNo),           & ! (out) OBL depth (m) dim=1
               kOBL_depth = kOBL_2D(jc,blockNo),            & ! (out) level (+fraction) of OBL extent dim=1
               zt_cntr    = cellHeight_KPP(jc,:,blockNo),   & ! (in) Height of cell centers (m) dim=(n_zlev)
               surf_fric  = ustar(jc,blockNo),              & ! (in) Turbulent friction velocity at surface (m/s) dim=1
               surf_buoy  = BuoyFlux_3D(jc,1,blockNo),      & ! (in) Buoyancy flux at surface (m2/s3) dim=1
               Coriolis   = Coriolis(jc,blockNo) )!,           & ! (in) Coriolis parameter (1/s) dim=1
               !CVmix_kpp_params_user = KPP_params)            ! (in) KPP parameter list                         


        ! A hack to avoid KPP reaching the bottom (only if deepOBLoffset is set to > 0m). 
        ! It was needed during development  in MOM6 because KPP was unable to handle 
        ! vanishingly small layers near the bottom
        IF (deepOBLoffset > 0.0_wp) THEN
          zBottomMinusOffset = iFaceHeight_KPP(jc,n_zlev+1,blockNo) &
               + MIN(deepOBLoffset,-0.1_wp*iFaceHeight_KPP(jc,n_zlev+1,blockNo))
          OBLdepth(jc,blockNo) = MIN( OBLdepth(jc,blockNo), -zBottomMinusOffset )
        END IF

        ! FIXME: use this only for debugging
        IF(fixedOBLdepth) THEN
          OBLdepth(jc,blockNo) = fixedOBLdepth_value
        END IF

        ! apply some constraints on OBLdepth

        ! no shallower than top layer
        OBLdepth(jc,blockNo) = MAX( OBLdepth(jc,blockNo), -iFaceHeight_KPP(jc,2,blockNo) )

        ! no deeper than bottom
        !OBLdepth(jc,blockNo) = MIN( OBLdepth(jc,blockNo), -iFaceHeight_KPP(jc,blockNo,n_zlev+1) )

        ! no deeper than kbot
        OBLdepth(jc,blockNo) = MIN( OBLdepth(jc,blockNo), -iFaceHeight_KPP(jc,kbot(jc,blockNo),blockNo) )

        ! prevent negative/too shallow OBL depths
        OBLdepth(jc,blockNo) = MAX( OBLdepth(jc,blockNo), minOBLdepth )

        ! model level of OBL depth (note: float number)
        kOBL_2D(jc,blockNo)  = cvmix_kpp_compute_kOBL_depth( iFaceHeight_KPP(jc,:,blockNo) &
             ,cellHeight_KPP(jc,:,blockNo), OBLdepth(jc,blockNo) )

        ! safety for kOBL: level of OBL must not exceed bottom level
        IF (kOBL_2D(jc,blockNo) .GT. kbot(jc,blockNo)) THEN
          kOBL_2D(jc,blockNo) = kbot(jc,blockNo)  ! or kbot(jc,blockNo) - 1 ???
        END IF
      
        END IF ! kbot
 
      ENDDO
    ENDDO

  END SUBROUTINE calc_OBLdepth

!========================================================================
!========================================================================
  ! This is the main routine to calculate the viscosity (K=h*w_m*G) and diffusivities
  ! (K=h*w_s*G), this K is also used for the non-local transport (NLT) terms,
  ! which are also calculated in this routine. In this routine the dimensionless
  ! uniform shape fuction G is always G(sigma)=sigma*(1-sigma)**2. 
  ! Note, however, that the NLT terms might be overwritten in the routine 
  ! 'calc_nonlocal_trans_heat_salt' if NLT_shape!="CVMIX", which would modify
  ! the shape of G.
  SUBROUTINE calc_diffusivities_within_OBL(patch_3d, params_oce,        &
                BuoyFlux_3D, cellHeight_KPP, iFaceHeight_KPP, OBLdepth, &
                kOBL_2D, ustar, avo_kpp, dvo_heat_kpp, dvo_salt_kpp,    &
                WM, WS, nonLocalTransHeat, nonLocalTransScalar)
    TYPE(t_patch_3d ),TARGET, INTENT(in)  :: patch_3d
    TYPE(t_subset_range), POINTER         :: all_cells
    TYPE(t_ho_params), INTENT(inout)      :: params_oce

    ! required variables
    REAL(wp), DIMENSION(:,:,:), INTENT(IN) ::         &
      BuoyFlux_3D                                    ,& ! buoyancy flux (m2 s-3)
      cellHeight_KPP                                 ,& ! depth of cell centres (m)
      iFaceHeight_KPP                                   ! depth of cell interfaces (m)

    REAL(wp), DIMENSION(:,:), INTENT(IN) ::           &
      OBLdepth                                       ,& ! ocean boundary layer depth (m)
      kOBL_2D                                        ,& ! model layer at OBL depth (float)
      ustar                                             ! surface friction velocity (m/s)

    REAL(wp), DIMENSION(:,:,:), INTENT(INOUT) ::      &
      avo_kpp                                        ,& ! vertical viscosity (m2 s-1)
      dvo_heat_kpp                                   ,& ! vertical diffusivity of heat (m2 s-1)
      dvo_salt_kpp                                   ,& ! vertical diffusivity of salt (m2 s-1)
      WM                                             ,& ! vertical turb. velocity scale for momentum
      WS                                             ,& ! vertical turb. velocity scale for scalar (tracer)
      nonLocalTransHeat                              ,& ! non-local transport term for heat
      nonLocalTransScalar                               ! non-local transport term for scalar (tracer)
      
    ! loop variables
    INTEGER :: jc, blockNo, jk
    INTEGER :: start_index, end_index
    !INTEGER :: levels

    ! local variables

    INTEGER     ::  &
      nlev         ,&
      max_nlev

    REAL(wp), DIMENSION(nproma,patch_3d%p_patch_2d(1)%alloc_cell_blocks)                 ::  &
      surfBuoyFlux                                        ! surface buyancy flux()

    REAL(wp), DIMENSION(nproma,n_zlev+1,patch_3d%p_patch_2d(1)%alloc_cell_blocks)        ::  &
      Kviscosity                                      !,& ! vertical viscosity at interfaces (m2/s) for CVMIX
      !G                                                  ! non-dimensional shape function  

    REAL(wp),DIMENSION(nproma,n_zlev+1,patch_3d%p_patch_2d(1)%alloc_cell_blocks,2)       ::  &
      Kdiffusivity                                    , & ! vertical diffusivity for heat/salt at interfaces (m2/s)
      nonLocalTrans                                       ! non-local transport for heat/salt at interfaces (non-dimensional)


    all_cells => patch_3d%p_patch_2d(1)%cells%ALL
    nlev = n_zlev
    max_nlev = nlev

    ! initialise surfBuoyFlux
    surfBuoyFlux = BuoyFlux_3D(:,1,:)  ! first level

    kbot => patch_3d%p_patch_1d(1)%dolic_c(:,:)

    DO blockNo = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blockNo, start_index, end_index)
      DO jc = start_index, end_index
        !levels       = patch_3d%p_patch_1d(1)%dolic_e(jc, blockNo)

        ! only for ocean cells
        IF (kbot(jc,blockNo)>0) THEN

       !========================================================================
        ! 9.1) Add option for use of surface buoyancy flux with total shortwave
        ! radiation flux.
        !========================================================================
        IF (SW_METHOD .EQ. "SW_METHOD_ALL_SW") THEN
          ! Use all shortwave radiation
          surfBuoyFlux(jc,blockNo) = BuoyFlux_3D(jc,1,blockNo)
        ELSEIF (SW_METHOD .EQ. "SW_METHOD_MXL_SW") THEN
          ! Use shortwave radiation absorbed in mixing layer
          ! FIXME: use INT(kOBL_2D) or FLOOR(kOBL_2D) ?  
          surfBuoyFlux(jc,blockNo)  = BuoyFlux_3D(jc,1,blockNo) - BuoyFlux_3D(jc,INT(kOBL_2D(jc,blockNo)),blockNo)
        ELSEIF (SW_METHOD .EQ. "SW_METHOD_LV1_SW") THEN
          ! Use shortwave radiation absorbed in layer 1
          surfBuoyFlux(jc,blockNo)  = BuoyFlux_3D(jc,1,blockNo) - BuoyFlux_3D(jc,2,blockNo)
        ENDIF

        !========================================================================
        ! 9.2) If option "MatchBoth" is selected in CVMix, ICON should be capable
        ! of matching.
        ! Unlike LMD94, we do not match to interior diffusivities. If using the
        ! original LMD94 shape function, not matching is equivalent to matching to
        ! a zero diffusivity.
        ! Note: matching is not recommended by CVMix --> might create extremes in 
        ! tracer field if G > 1 at mld base, causing a larger redistribution of
        ! the surface fluxes than what is available at the surface.
        !========================================================================

        IF (.NOT. (MatchTechnique.EQ.'MatchBoth')) THEN
          Kdiffusivity(jc,:,blockNo,:) = 0.0_wp ! Diffusivities for heat and salt (m2/s)
          Kviscosity(jc,:,blockNo)     = 0.0_wp ! Viscosity (m2/s)
        ELSE
          Kdiffusivity(jc,:,blockNo,1) = dvo_heat_kpp(jc,:,blockNo)
          Kdiffusivity(jc,:,blockNo,2) = dvo_salt_kpp(jc,:,blockNo)
          Kviscosity(jc,:,blockNo)     = avo_kpp(jc,:,blockNo)
        ENDIF


        !========================================================================
        ! 9.3) Calculate the turbulent diffusion coefficients
        !========================================================================
        ! NOTE: the turbulent velocity scale is calculated anew with the correct
        ! OBL depth within the CALL
        ! NOTE: ! they dont make use of 'LangmuirEnhancementFactor' -> always
        ! default is used ?

        CALL cvmix_coeffs_kpp(                                   &
                     Mdiff_out  = Kviscosity(jc,:,blockNo),      & ! (inout) Total viscosity (m2/s) (max_nlev+1)
                     Tdiff_out  = Kdiffusivity(jc,:,blockNo,1),  & ! (inout) Total heat diffusivity (m2/s) (max_nlev+1)
                     Sdiff_out  = Kdiffusivity(jc,:,blockNo,2),  & ! (inout) Total salt diffusivity (m2/s) (max_nlev+1)
                     zw         = iFaceHeight_KPP(jc,:,blockNo), & ! (in)    Height of interfaces (m) (max_nlev+1)
                     zt         = cellHeight_KPP(jc,:,blockNo),  & ! (in)    Height of level centers (m) (max_nlev)
                     old_Mdiff  = Kviscosity(jc,:,blockNo),      & ! (in)    Original viscosity (m2/s) (max_nlev+1)
                     old_Tdiff  = Kdiffusivity(jc,:,blockNo,1),  & ! (in)    Original heat diffusivity (m2/s) (max_nlev+1)
                     old_Sdiff  = Kdiffusivity(jc,:,blockNo,2),  & ! (in)    Original salt diffusivity (m2/s) (max_nlev+1)
                     OBL_depth  = OBLdepth(jc,blockNo),          & ! (in)    OBL depth (m) (scalar)
                     kOBL_depth = kOBL_2D(jc,blockNo),           & ! (in)    level (+fraction) of OBL extent (scalar)
                     Tnonlocal  = nonLocalTrans(jc,:,blockNo,1), & ! (out)   Non-local heat transport (non-dimensional) (max_nlev+1)
                     Snonlocal  = nonLocalTrans(jc,:,blockNo,2), & ! (out)   Non-local salt transport (non-dimensional) (max_nlev+1)
                     surf_fric  = ustar(jc,blockNo),             & ! (in)    Turbulent friction velocity at surface (m/s) (scalar)
                     surf_buoy  = surfBuoyFlux(jc,blockNo),      & ! (in)    Buoyancy flux at surface (m2/s3) (scalar)
                     nlev       = nlev,                          & ! (in)    Number of levels to compute coeffs for (scalar)
                     max_nlev   = max_nlev,                      & ! (in)    maximum vertical levels (scalar)
                     langmuir_Efactor = langmuir_Efactor,        & ! (in)    Enhancement factor due to Langmuir circulation (scalar)
                     !CVmix_kpp_params_user = KPP_params,         & ! (in)    KPP parameter list 
                     w_m        = WM(jc,:,blockNo),              & ! (out)   turbulent velocity scale for momentum (nlev+1)
                     w_s        = WS(jc,:,blockNo)  )              ! (out)   turbulent velocity scale for tracer (nlev+1)
                     !G=G(jc,:,blockNo)             )              ! (out)   shape function G(sigma)

        ! viscosity and diffusivity at cell centres
        avo_kpp(jc,:,blockNo) = Kviscosity(jc,:,blockNo)

        dvo_heat_kpp(jc,:,blockNo) = Kdiffusivity(jc,:,blockNo,1)
        dvo_salt_kpp(jc,:,blockNo) = Kdiffusivity(jc,:,blockNo,2)

        ! KPP_NonLocalTransport_temp and KPP_NonLocalTransport_saln
        nonLocalTransHeat(jc,:,blockNo)   = nonLocalTrans(jc,:,blockNo,1) ! temp
        nonLocalTransScalar(jc,:,blockNo) = nonLocalTrans(jc,:,blockNo,2) ! salt

        END IF !ocean only

      ENDDO
    ENDDO

  END SUBROUTINE calc_diffusivities_within_OBL

!========================================================================
!========================================================================
  ! This routine calculates the viscosity and diffusivities below the OBL. 
  ! Currently, two options are implemented 'internal_mix_scheme=PP' calls
  ! the Pacanowski & Philander (1981) scheme and 'internal_mix_scheme=KPP' 
  ! calls the Large et al. (1994) Ri-dependent scheme.
  ! If 'convection='.TRUE.' the diffusivities are enhanced to parameterize
  ! convection. The viscosity is not enhanced.

  SUBROUTINE calc_diffusivities_below_OBL(patch_3d, params_oce, & 
                richardson_no, kOBL_2D, avo_kpp, dvo_heat_kpp, dvo_salt_kpp)
    TYPE(t_patch_3d ),TARGET, INTENT(in)  :: patch_3d
    !TYPE(t_patch),POINTER                 :: patch_2D
    TYPE(t_subset_range), POINTER         :: all_cells
    TYPE(t_ho_params), INTENT(inout)      :: params_oce

    ! required variables
    REAL(wp), DIMENSION(:,:,:), INTENT(IN) ::         &
      richardson_no

    REAL(wp), DIMENSION(:,:), INTENT(IN) ::           &
      kOBL_2D

    REAL(wp), DIMENSION(:,:,:), INTENT(INOUT) ::      & 
      avo_kpp                                        ,&
      dvo_heat_kpp                                   ,&
      dvo_salt_kpp 
 
    ! loop variables
    INTEGER :: jc, blockNo, tracer_index, kml
    INTEGER :: start_index, end_index

    INTEGER :: kobl  ! first model level below mixed layer

    INTEGER :: cell_1_idx, cell_1_block, cell_2_idx, cell_2_block
    REAL(wp) :: denom 

    ! parameters from namelist
    !   internal_mix_scheme
    !   PP_alpha, PP_nu_zero, PP_loc_exp
    !   KPP_nu_zero, KPP_Ri_zero, KPP_loc_exp 

    !--------------------------------------------------------------------------
    ! Richardson-dependent mixing below the OBL (either 'PP' or 'KPP')
    !--------------------------------------------------------------------------
    ! PP:  PP scheme 
    ! KPP: Ri-dependent scheme as described in LMD94

    !patch_2D => patch_3d%p_patch_2d(1)
    all_cells => patch_3d%p_patch_2d(1)%cells%ALL
    kbot => patch_3d%p_patch_1d(1)%dolic_c(:,:)

    DO blockNo = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blockNo, start_index, end_index)
      DO jc = start_index, end_index
       
        ! only for ocean cells
        IF (kbot(jc,blockNo)>0) THEN

        ! index of first model level below mixed layer
        kobl = INT(kOBL_2D(jc,blockNo)) + 1

         SELECT CASE (internal_mix_scheme)

          CASE('PP')
            !========================================================================
            ! A) Pacanowski-Philander (PP81)
            !========================================================================
              DO kml=kobl,kbot(jc,blockNo) 
              !DO kml=kobl,kbot(jc,blockNo) + 1 (as in cvmix?)
 
!               ! viscosity
!               avo_kpp(jc,kml,blockNo) = PP_nu_zero / ((1.0_wp + PP_alpha *          &
!                                         richardson_no(jc,kml,blockNo))**PP_loc_exp) &
!                                         + params_oce%a_veloc_v_back

               IF (richardson_no(jc,kml,blockNo) .GT. 0.0_wp) THEN
                 denom = 1.0_wp + PP_alpha * richardson_no(jc,kml,blockNo)
               ELSE
                 ! Treat non-negative Richardson number as Ri = 0
                 denom = 1.0_wp
               END IF


               ! viscosity
               avo_kpp(jc,kml,blockNo) = PP_nu_zero / (denom)**PP_loc_exp + params_oce%a_veloc_v_back
 
               ! FIXME: note that below OBL, there is no distinguishing between
               ! tracers, they all use the same diffusivity

               ! diffusivity
               DO tracer_index=1,no_tracer

                 dvo_heat_kpp(jc,kml,blockNo) = PP_nu_zero / (denom**3) + params_oce%a_tracer_v_back(tracer_index)
                 dvo_salt_kpp(jc,kml,blockNo) = dvo_heat_kpp(jc,kml,blockNo)

               ENDDO

             
             ENDDO


         CASE('KPP')
            !========================================================================
            ! B) KPP shear parameterization below the mixed layer (Large et
            ! al.,1994)
            ! FIXME: maybe change to call subroutine   subroutine
            ! cvmix_coeffs_shear_low(Mdiff_out, Tdiff_out, RICH, nlev,         &
            !                         max_nlev, CVmix_shear_params_user)
            ! where CVmix_shear_params_user%mixing_scheme='KPP'
            !========================================================================
              DO kml=kobl,kbot(jc,blockNo)
              !DO kml=kobl,kbot(jc,blockNo) + 1 (as in cvmix?)

                ! viscosity
                DO tracer_index=1,no_tracer

                  IF(richardson_no(jc,kml,blockNo) .LT. 0.0_wp) THEN
                    !params_oce%a_tracer_v(jc,kml,blockNo,tracer_index) = KPP_nu_zero 
                    dvo_heat_kpp(jc,kml,blockNo) = KPP_nu_zero

                  ELSEIF(richardson_no(jc,kml,blockNo) .LT. KPP_Ri_zero) THEN
                    !params_oce%a_tracer_v(jc,kml,blockNo,tracer_index) = KPP_nu_zero * (1.0_wp - &
                    !                                                     (richardson_no(jc,kml,blockNo) &
                    !                                                      / KPP_Ri_zero)**2)**KPP_loc_exp
                    dvo_heat_kpp(jc,kml,blockNo) = KPP_nu_zero * (1.0_wp - &
                                                   (richardson_no(jc,kml,blockNo) &
                                                   / KPP_Ri_zero)**2)**KPP_loc_exp
                    dvo_salt_kpp(jc,kml,blockNo) = dvo_heat_kpp(jc,kml,blockNo) 
                  ELSE ! Ri_g >= Ri_zero
                    !params_oce%a_tracer_v(jc,kml,blockNo,tracer_index) = 0.0_wp
                    dvo_heat_kpp(jc,kml,blockNo) = 0.0_wp
                    dvo_salt_kpp(jc,kml,blockNo) = dvo_heat_kpp(jc,kml,blockNo)

                  END IF

                ENDDO

                ! set viscosity equal to heat diffusivity
                ! FIXME: how to deal with more than one tracer here?
                ! Is this how it is done in Large et al. (1994)?
                !avo_kpp(jc,kml,blockNo) = params_oce%a_tracer_v(jc,kml,blockNo,1)
                 avo_kpp(jc,kml,blockNo) = dvo_heat_kpp(jc,kml,blockNo)

              ENDDO


         END SELECT

          END IF ! ocean only
       ENDDO
     ENDDO

  END SUBROUTINE calc_diffusivities_below_OBL

!========================================================================
!========================================================================
  ! This routine overwrites the non-local transport (NLT) terms, if
  ! NLT_shape!="CVMIX". 
  ! The CVMix code has yet to update for these options, so we compute it here.
  ! Note that nonLocalTrans = Cs * G(sigma) (LMD94 notation), with Cs =
  ! 6.32739901508.
  ! IF NLT_shape="CVMIX", then nothing will happen.
  SUBROUTINE calc_nonlocal_trans_heat_salt(patch_3d, params_oce,           &
                iFaceHeight_KPP, BuoyFlux_3D, OBLdepth, nonLocalTransHeat, &
                nonLocalTransScalar)
    TYPE(t_patch_3d ),TARGET, INTENT(in)  :: patch_3d
    TYPE(t_subset_range), POINTER         :: all_cells
    TYPE(t_ho_params), INTENT(inout)      :: params_oce

    ! required variables
    REAL(wp), DIMENSION(:,:,:), INTENT(IN) ::         &
      BuoyFlux_3D                                    ,&
      iFaceHeight_KPP

    REAL(wp), DIMENSION(:,:), INTENT(IN) ::           &
      OBLdepth

    REAL(wp), DIMENSION(:,:,:), INTENT(INOUT) ::      &
      nonLocalTransHeat                              ,&
      nonLocalTransScalar                

    ! local variables 
    REAL(wp)                                     ::  &
      sigma                                            ! non-dimensional depth in OBL

    REAL(wp) :: surfBuoyFlux(nproma,patch_3d%p_patch_2d(1)%alloc_cell_blocks) ! surface buoyancy flux (unit?)

    ! non-local transport terms (heat, scalar) (unit?)
    REAL(wp) :: nonLocalTrans(nproma,n_zlev+1,patch_3d%p_patch_2d(1)%alloc_cell_blocks,2)

    ! loop variables
    INTEGER :: jc, blockNo, jk
    INTEGER :: start_index, end_index
    INTEGER :: levels

    ! store first level of buoyancy flux
    surfBuoyFlux= BuoyFlux_3D(:,1,:)

    ! might be overwritten in the following if not NLT_shape="CVMIX" is used.
    nonLocalTrans(:,:,:,1) = nonLocalTransHeat(:,:,:)
    nonLocalTrans(:,:,:,2) = nonLocalTransScalar(:,:,:)
    
    !========================================================================
    ! 12) Calculate nonlocal transport terms
    !========================================================================
    ! Over-write CVMix NLT shape function with one of the following choices.
    !
    ! Start do-loop at k=2, since k=1 is ocean surface (sigma=0)
    ! and we do not wish to double-count the surface forcing.
    ! Only compute nonlocal transport for 0 <= sigma <= 1.
    ! MOM6 recommended shape is the parabolic; it gives deeper boundary layer
    ! and no spurious extrema.
    !
    ! --> the non-local transport is then applied in additional subroutines
    !     see subroutines KPP_NonLocalTransport_temp and
    !     KPP_NonLocalTransport_saln below.
    kbot => patch_3d%p_patch_1d(1)%dolic_c(:,:)
    all_cells => patch_3d%p_patch_2d(1)%cells%ALL

    DO blockNo = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blockNo, start_index, end_index)
      DO jc = start_index, end_index
 
      levels       = patch_3d%p_patch_1d(1)%dolic_e(jc, blockNo)

      ! only for ocean cells
      IF (kbot(jc,blockNo)>0) THEN


      ! only if surface buoyancy flux is negative
      IF (surfBuoyFlux(jc,blockNo) < 0.0_wp) THEN

        ! nonLocalTrans is only changed if NLT_shape != "CVMIX"
        IF (NLT_shape == "CUBIC") THEN
          DO jk = 2, levels
            sigma = MIN(1.0_wp,-iFaceHeight_KPP(jc,jk,blockNo)/OBLdepth(jc,blockNo))
            nonLocalTrans(jc,jk,blockNo,1) = (1.0_wp - sigma)**2 * (1.0_wp + 2.0_wp*sigma)!*cs2
            nonLocalTrans(jc,jk,blockNo,2) = nonLocalTrans(jc,jk,blockNo,1)
          ENDDO
        ELSEIF (NLT_shape == "PARABOLIC") THEN
          DO jk = 2, levels
            sigma = MIN(1.0_wp,-iFaceHeight_KPP(jc,jk,blockNo)/OBLdepth(jc,blockNo))
            nonLocalTrans(jc,jk,blockNo,1) = (1.0_wp - sigma)**2 !*cs2
            nonLocalTrans(jc,jk,blockNo,2) = nonLocalTrans(jc,jk,blockNo,1)
          ENDDO
        ELSEIF (NLT_shape == "LINEAR") THEN
          DO jk = 2, levels
            sigma = MIN(1.0_wp,-iFaceHeight_KPP(jc,jk,blockNo)/OBLdepth(jc,blockNo))
            nonLocalTrans(jc,jk,blockNo,1) = (1.0_wp - sigma)!*cs2
            nonLocalTrans(jc,jk,blockNo,2) = nonLocalTrans(jc,jk,blockNo,1)
          ENDDO
        ELSEIF (NLT_shape == "CUBIC_LMD") THEN
          ! Sanity check (should agree with CVMix result using simple matching)
          DO jk = 2, levels
            sigma = MIN(1.0_wp,-iFaceHeight_KPP(jc,jk,blockNo)/OBLdepth(jc,blockNo))
            nonLocalTrans(jc,jk,blockNo,1) = cs2 * sigma*(1.0_wp -sigma)**2
            nonLocalTrans(jc,jk,blockNo,2)= nonLocalTrans(jc,jk,blockNo,1)
          ENDDO
        ENDIF

      ENDIF

      ! Save non-local transport terms in module variables
      ! KPP_NonLocalTransport_temp and KPP_NonLocalTransport_saln
      nonLocalTransHeat(jc,:,blockNo)   = nonLocalTrans(jc,:,blockNo,1) ! temp
      nonLocalTransScalar(jc,:,blockNo) = nonLocalTrans(jc,:,blockNo,2) ! salt

      END IF ! ocean only

    ENDDO
  ENDDO

  !------------------------------------------------------------------------
  ! Apply non-local transport of heat and salt
  !------------------------------------------------------------------------

  ! non-local contribution is applied in tracer_rhs directly
  ! (see mo_ocean_tracer.f90) 


  END SUBROUTINE calc_nonlocal_trans_heat_salt


!========================================================================
!========================================================================
  ! This routine calculates the non-local transport tendency for temperature.
  ! Note that the tendency is applied in the r.h.s. of the tracer diffusion
  ! equation elsewhere in the code.

  SUBROUTINE KPP_NonLocalTransport_temp(patch_3d, params_oce, concsum, &
                    nonLocalTransHeat, net_heat, nl_trans_tend_heat)
    TYPE(t_patch_3d ),TARGET, INTENT(in)  :: patch_3d
    TYPE(t_subset_range), POINTER         :: all_cells
    TYPE(t_ho_params), INTENT(inout)      :: params_oce
    REAL(wp), TARGET                      :: concsum(:,:) ! t_sea_ice%concs

    ! required variables
    REAL(wp), DIMENSION(:,:,:), INTENT(IN) ::         &
      nonLocalTransHeat                              

    REAL(wp), DIMENSION(:,:), INTENT(IN) ::           &
      net_heat

    REAL(wp), DIMENSION(:,:,:), INTENT(INOUT) ::      &
      nl_trans_tend_heat


    ! loop variables
    INTEGER :: jc, blockNo, jk
    INTEGER :: start_index, end_index
    INTEGER :: levels

    ! This routine calculates the non-local transport tendency for temperature

    kbot => patch_3d%p_patch_1d(1)%dolic_c(:,:)
    dzw => patch_3d%p_patch_1d(1)%prism_thick_c
    all_cells => patch_3d%p_patch_2d(1)%cells%ALL

    DO blockNo = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blockNo, start_index, end_index)
      DO jc = start_index, end_index

        levels = patch_3d%p_patch_1d(1)%dolic_e(jc, blockNo)

        ! only for ocean cells
        IF (kbot(jc,blockNo)>0) THEN


        IF (.NOT. lnl_trans_under_sea_ice ) THEN
          ! no non-local transports under sea ice
          IF (concsum(jc,blockNo) .GT. 0.0_wp ) CYCLE
        ENDIF


        DO jk = 1, levels

          ! net_heat already contains clw
          nl_trans_tend_heat(jc,jk,blockNo) = -(nonLocalTransHeat(jc,jk+1,blockNo) - nonLocalTransHeat(jc,jk,blockNo) ) / &
                                               (dzw(jc,jk,blockNo)+dbl_eps) * net_heat(jc,blockNo)

       ENDDO !jk

       END IF ! ocean only

      ENDDO !jc
    ENDDO !blockNo



  END SUBROUTINE KPP_NonLocalTransport_temp 
  
!========================================================================
!========================================================================
  ! This routine calculates the non-local transport tendency for salinity.
  ! Note that the tendency is applied in the r.h.s. of the tracer diffusion
  ! equation elsewhere in the code.

  SUBROUTINE KPP_NonLocalTransport_saln(patch_3d, params_oce, concsum, &
               nonLocalTransScalar, net_salt, nl_trans_tend_salt)
    TYPE(t_patch_3d ),TARGET, INTENT(in)  :: patch_3d
    TYPE(t_subset_range), POINTER         :: all_cells
    TYPE(t_ho_params), INTENT(inout)      :: params_oce
    REAL(wp), TARGET                      :: concsum(:,:) ! t_sea_ice%concs

    ! required variables
    REAL(wp), DIMENSION(:,:,:), INTENT(IN) ::         &
      nonLocalTransScalar

    REAL(wp), DIMENSION(:,:), INTENT(IN) ::           &
      net_salt

    REAL(wp), DIMENSION(:,:,:), INTENT(INOUT) ::      &
      nl_trans_tend_salt

    ! loop variables
    INTEGER :: jc, blockNo, jk
    INTEGER :: start_index, end_index
    INTEGER :: levels

    ! This routine calculates the non-local transport tendency for salinity
    kbot => patch_3d%p_patch_1d(1)%dolic_c(:,:)
    dzw => patch_3d%p_patch_1d(1)%prism_thick_c
    all_cells => patch_3d%p_patch_2d(1)%cells%ALL

    DO blockNo = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blockNo, start_index, end_index)
      DO jc = start_index, end_index

        levels = patch_3d%p_patch_1d(1)%dolic_e(jc, blockNo)

        ! only for ocean cells
        IF (kbot(jc,blockNo)>0) THEN

        IF (.NOT. lnl_trans_under_sea_ice ) THEN
          ! no non-local transports under sea ice
          IF (concsum(jc,blockNo) .GT. 0.0_wp ) CYCLE
        ENDIF

  
        DO jk = 1, levels

          nl_trans_tend_salt(jc,jk,blockNo) = -( nonLocalTransScalar(jc,jk+1,blockNo) - nonLocalTransScalar(jc,jk,blockNo) ) / &
                                               (dzw(jc,jk,blockNo)+dbl_eps) * net_salt(jc,blockNo)

        ENDDO !jk

        END IF ! ocean only

      ENDDO !jc
    ENDDO !blockNo

  END SUBROUTINE KPP_NonLocalTransport_saln

!========================================================================
!========================================================================
  ! This routine calculates the vertical turbulent velocity scales, for momentum
  ! (w_m) and for tracers (w_s). These scales are required to calculate
  ! K=h*w_s*G.
  SUBROUTINE calc_diag_turb_velocity_scales(patch_3d, BuoyFlux_3D, iFaceHeight_KPP, &
                OBLdepth, ustar, WM, WS)
    TYPE(t_patch_3d ),TARGET, INTENT(in)  :: patch_3d
    TYPE(t_subset_range), POINTER         :: all_cells

    ! required variables
    REAL(wp), DIMENSION(:,:,:), INTENT(IN) ::         &
      BuoyFlux_3D                                    ,&
     !cellHeight_KPP
      iFaceHeight_KPP                                

    REAL(wp), DIMENSION(:,:), INTENT(IN) ::           & 
      OBLdepth                                       ,&
      ustar                                          

    REAL(wp), DIMENSION(:,:,:), INTENT(INOUT) ::      &
      WM                                             ,&
      WS


    ! loop variables
    INTEGER :: jc, blockNo
    INTEGER :: start_index, end_index
    !INTEGER :: levels

    REAL(wp), DIMENSION(nproma,patch_3d%p_patch_2d(1)%alloc_cell_blocks)::  &
      surfBuoyFlux                                        ! surface buyancy flux()

    ! parameter from namelist
    !  langmuir_Efactor

    !========================================================================
    ! Diagnostic calculation of turbulent vertical velocity scale for tracer
    ! (WS)
    !========================================================================
    ! recompute ws for diagnostics, now that we in fact know boundary
    ! layer depth.
    ! Note: the correct ws is calculated within cvmix_coeffs_kpp() 
    all_cells => patch_3d%p_patch_2d(1)%cells%ALL
    kbot => patch_3d%p_patch_1d(1)%dolic_c(:,:)
    
    surfBuoyFlux= BuoyFlux_3D(:,1,:)
 

    DO blockNo = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blockNo, start_index, end_index)
      DO jc = start_index, end_index

        ! only for ocean cells
        IF (kbot(jc,blockNo)>0) THEN

        ! output from subroutine is: WS
        CALL cvmix_kpp_compute_turbulent_scales(                      &
               sigma_coord     = -iFaceHeight_KPP(jc,:,blockNo)/OBLdepth(jc,blockNo), & ! (in) Normalized depth (z_nlev+1)
               !sigma_coord     = -cellHeight_KPP(jc,:,blockNo)/OBLdepth(jc,blockNo),  & ! (in) Normalized depth (z_nlev)
               OBL_depth       = OBLdepth(jc,blockNo),                & ! (in) OBL depth (m)
               surf_buoy_force = surfBuoyFlux(jc,blockNo),            & ! (in) surfacve buoyancy flux (m2/s3)
               surf_fric_vel   = ustar(jc,blockNo),                   & ! (in) Turbulent friction velocity at surface (m/s)
               langmuir_Efactor= langmuir_Efactor,                    & ! (in) Langmuir factor
               w_m             = WM(jc,:,blockNo),                    & ! (out) Turb. velocity scale for momentum (m/s) (z_nlev+1) 
               w_s             = WS(jc,:,blockNo) ) !,                    & ! (out) Turb. velocity scale for tracer (m/s) (z_nlev+1)
             !  CVmix_kpp_params_user = KPP_params                     &
             !)
        END IF ! ocean only
      ENDDO
    ENDDO

  END SUBROUTINE calc_diag_turb_velocity_scales

!========================================================================
!========================================================================
  ! This routine diagnoses the non-dimensional shape function G(sigma) =
  ! sigma*(1-sigma)**2. G(sigma) is used to calculate the local K = h*ws*G
  ! and for the non-local transport (NLT) terms. However, for the NLT terms
  ! G can have a different form if NLT_shape!="CVMIX", then the diagnosed form here
  ! is only correct for the local K.

  SUBROUTINE diagnose_G(patch_3d, params_oce, &
               iFaceHeight_KPP, OBLdepth, kOBL_2D, G_3D)
    TYPE(t_patch_3d ),TARGET, INTENT(in)  :: patch_3d
    TYPE(t_subset_range), POINTER         :: all_cells
    TYPE(t_ho_params), INTENT(inout)      :: params_oce

    ! required variables
    REAL(wp), DIMENSION(:,:,:), INTENT(IN) ::         &
      iFaceHeight_KPP

    REAL(wp), DIMENSION(:,:), INTENT(IN) ::           &
      OBLdepth                                       ,&
      kOBL_2D

    REAL(wp), DIMENSION(:,:,:), INTENT(INOUT) ::      &
      G_3D

    ! loop variables
    INTEGER :: jc, blockNo, jk
    INTEGER :: start_index, end_index

    ! local variables
    REAL(wp), DIMENSION(4)                     ::  &
      Mshape                                      !,&
      !Tshape,                                      &
      !Sshape
    
 
    ! diagnose the non-dimensional shape function G(sigma) = sigma*(1-sigma)**2
    Mshape(1) = 0.0_wp
    Mshape(2) = 1.0_wp
    Mshape(3) = -2.0_wp
    Mshape(4) = 1.0_wp
    !Tshape = Mshape
    !Sshape = Mshape
    all_cells => patch_3d%p_patch_2d(1)%cells%ALL
    kbot => patch_3d%p_patch_1d(1)%dolic_c(:,:)
 
    DO blockNo = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blockNo, start_index, end_index)
      DO jc = start_index, end_index

        ! only for ocean cells
        IF (kbot(jc,blockNo)>0) THEN

        DO jk=2,FLOOR(kOBL_2D(jc,blockNo))
         G_3D(jc,jk,blockNo)=cvmix_math_evaluate_cubic(Mshape,MIN(1.0_wp,-iFaceHeight_KPP(jc,jk,blockNo)/OBLdepth(jc,blockNo)))
        ENDDO


         END IF ! ocean only
      ENDDO
    ENDDO

  END SUBROUTINE diagnose_G

!========================================================================
!========================================================================
  
  SUBROUTINE diagnose_vtsqr(patch_3d, params_oce,     &
              cellHeight_KPP, WS_cntr, Nbuoy, Vt2)
    TYPE(t_patch_3d ),TARGET, INTENT(in)  :: patch_3d
    TYPE(t_subset_range), POINTER         :: all_cells
    TYPE(t_ho_params), INTENT(inout)      :: params_oce

    ! required variables
    REAL(wp), DIMENSION(:,:,:), INTENT(IN) ::         &
      cellHeight_KPP                                 ,&
      WS_cntr                                        ,&
      Nbuoy

    REAL(wp), DIMENSION(:,:,:), INTENT(INOUT) ::      &
     Vt2 

    ! loop variables
    INTEGER :: jc, blockNo
    INTEGER :: start_index, end_index
    !INTEGER :: levels

    all_cells => patch_3d%p_patch_2d(1)%cells%ALL
    kbot => patch_3d%p_patch_1d(1)%dolic_c(:,:)

    DO blockNo = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blockNo, start_index, end_index)
      DO jc = start_index, end_index

        ! only for ocean cells
        IF (kbot(jc,blockNo)>0) THEN

          ! compute unresolved squared velocity for diagnostics
          ! Note if Vt2 = MIN(Vt2, minVtsqr)
          Vt2(jc,:,blockNo) = cvmix_kpp_compute_unresolved_shear( &
                     zt_cntr = -cellHeight_KPP(jc,:,blockNo),     & ! Depth of cell center (m) (z_nlev)
                     ws_cntr = WS_cntr(jc,:,blockNo),             & ! Turbulent velocity scale profile, at centers (m/s)
                     N_iface = Nbuoy(jc,:,blockNo) )!,            & ! Buoyancy frequency at interface (1/s)
                     !cvmix_kpp_params_user = KPP_params ) ! KPP parameters

         END IF ! ocean only
      ENDDO
    ENDDO

  END SUBROUTINE diagnose_vtsqr

!========================================================================
!========================================================================

  SUBROUTINE calc_enhanced_vert_diff_below_OBL(patch_3d, params_oce, & 
                      stab_KPP, kOBL_2D, dvo_heat_kpp, dvo_salt_kpp)
    ! This is a parameterisation of convection by greatly enhanced vertical diffusion. 
    ! In case of the convection parametrizations ioconv = 1 and ioconv = 4, the vertical eddy
    ! diffusivities for tracers are increased to the namelist parameter
    ! convection_InstabilityThreshold.
    ! Note: the eddy viscosity is not enhanced.
    TYPE(t_patch_3d ),TARGET, INTENT(in)  :: patch_3d
    TYPE(t_subset_range), POINTER         :: all_cells
    TYPE(t_ho_params), INTENT(inout)      :: params_oce

    ! required variables
    REAL(wp), DIMENSION(:,:,:), INTENT(IN) ::         &
      stab_KPP                                       

    REAL(wp), DIMENSION(:,:), INTENT(IN) ::           &
      kOBL_2D    

    REAL(wp), DIMENSION(:,:,:), INTENT(INOUT) ::      &
      dvo_heat_kpp                                   ,&
      dvo_salt_kpp       

    ! local variables
    REAL(wp) :: diffusion_weight

    ! loop variables
    INTEGER :: jc, blockNo, tracer_index, kml
    INTEGER :: start_index, end_index
    !INTEGER :: levels
    INTEGER :: kobl  ! first model level below mixed layer
  
    all_cells => patch_3d%p_patch_2d(1)%cells%ALL
    kbot => patch_3d%p_patch_1d(1)%dolic_c(:,:)
 
    DO blockNo = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blockNo, start_index, end_index)
      DO jc = start_index, end_index

        ! only for ocean cells
        IF (kbot(jc,blockNo)>0) THEN

        ! index of first model level below mixed layer
        kobl = INT(kOBL_2D(jc,blockNo)) + 1

        DO kml=kobl,kbot(jc,blockNo)

          !convection
          IF ( stab_KPP(jc,kml,blockNo) <= convection_InstabilityThreshold) THEN
            ! no convection for velocity
            !avo_kpp(jc,kml,blockNo) = MAX(cdvocon * (convection_InstabilityThreshold - stab_KPP(jc,kml,blockNo)) &
            !     / (convection_InstabilityThreshold + ABS(stabio(i, j, k))), avo_conv(i,j,k))
            dvo_heat_kpp(jc,kml,blockNo) = tracer_convection_MixingCoefficient
            dvo_salt_kpp(jc,kml,blockNo) = tracer_convection_MixingCoefficient
    
          ELSE

            IF ( stab_KPP(jc,kml,blockNo) < RichardsonDiffusion_threshold) THEN
            ! loop over tracer

              !params_oce%a_tracer_v(jc,kml,blockNo,tracer_index) = MAX(tracer_convection_MixingCoefficient &
              !     * (convection_InstabilityThreshold - stab_KPP(jc,kml,blockNo)) &
              !     / (convection_InstabilityThreshold + ABS(stab_KPP(jc,kml,blockNo))), &
              !     params_oce%a_tracer_v(jc,kml,blockNo,tracer_index))

              diffusion_weight =  &
                  & (stab_KPP(jc,kml,blockNo) - convection_InstabilityThreshold) / &
                  & (RichardsonDiffusion_threshold - convection_InstabilityThreshold)

              dvo_heat_kpp(jc,kml,blockNo) = &
                  &  tracer_convection_MixingCoefficient * (1.0_wp - diffusion_weight) +&
                  & diffusion_weight * dvo_heat_kpp(jc,kml,blockNo)

              dvo_salt_kpp(jc,kml,blockNo) = &
                  &  tracer_convection_MixingCoefficient * (1.0_wp - diffusion_weight) +&
                  & diffusion_weight * dvo_salt_kpp(jc,kml,blockNo)

              

              !dvo_heat_kpp(jc,kml,blockNo)  = MAX(tracer_convection_MixingCoefficient &
              !     * (convection_InstabilityThreshold - stab_KPP(jc,kml,blockNo)) &
              !     / (convection_InstabilityThreshold + ABS(stab_KPP(jc,kml,blockNo))), &
              !     dvo_heat_kpp(jc,kml,blockNo))
   
              !dvo_salt_kpp(jc,kml,blockNo) = dvo_heat_kpp(jc,kml,blockNo)
 
           END IF 
            
          END IF !convection

        ENDDO ! kml



        END IF ! ocean only

      ENDDO !jc
    ENDDO !blockNo

  
  END SUBROUTINE calc_enhanced_vert_diff_below_OBL

!========================================================================
!========================================================================

  SUBROUTINE update_mixing_coefficients(patch_3d, params_oce, &
                avo_kpp, dvo_heat_kpp, dvo_salt_kpp)
    ! updates the viscosity and diffusivity fields
    ! Note: it is not neccessary to distinguish within and below OBL
    !       because this was already accounted for in the calculation
    !       of the mixing coefficients
    TYPE(t_patch_3d ),TARGET, INTENT(in)  :: patch_3d
    TYPE(t_patch),POINTER                 :: patch_2D
    TYPE(t_subset_range), POINTER         :: all_cells,edges_in_domain
    TYPE(t_ho_params), INTENT(inout)      :: params_oce

    ! required variables
    REAL(wp), DIMENSION(:,:,:), INTENT(IN) ::         &
      avo_kpp                                        ,&
      dvo_heat_kpp                                   ,&
      dvo_salt_kpp                                  

    !REAL(wp), DIMENSION(:,:), INTENT(IN) ::           &
    !  kOBL_2D

    ! loop variables
    INTEGER :: jc, jk, blockNo, tracer_index, kml, je
    INTEGER :: start_index, end_index
    INTEGER :: levels 

    INTEGER :: cell_1_idx, cell_1_block, cell_2_idx, cell_2_block

    patch_2D => patch_3d%p_patch_2d(1)
    all_cells => patch_3d%p_patch_2d(1)%cells%ALL
    edges_in_domain => patch_2D%edges%in_domain
 

    ! interpolate vert. visosity from cell center to edges
    params_oce%a_veloc_v = 0.0
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_index, end_index)
      DO je = start_index, end_index
        levels       = patch_3d%p_patch_1d(1)%dolic_e(je, blockNo)
        cell_1_idx   = patch_2D%edges%cell_idx(je,blockNo,1)
        cell_1_block = patch_2D%edges%cell_blk(je,blockNo,1)
        cell_2_idx   = patch_2D%edges%cell_idx(je,blockNo,2)
        cell_2_block = patch_2D%edges%cell_blk(je,blockNo,2)

        ! update viscosity 
        !params_oce%a_veloc_v = 0.0
        DO jk = 2, levels
          params_oce%a_veloc_v(je,jk,blockNo) = &
          & 0.5_wp * (    avo_kpp(cell_1_idx,jk,cell_1_block) &
          &             + avo_kpp(cell_2_idx,jk,cell_2_block) )
        END DO

      END DO
    END DO

    ! write kpp vert. diffusivity to vert tracer diffusivities 
    DO blockNo = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blockNo, start_index, end_index)
      DO jc = start_index, end_index
        ! update diffusivity

        ! loop over tracer
        !DO tracer_index=1,no_tracer
        params_oce%a_tracer_v(jc,:,blockNo,1) = dvo_heat_kpp(jc,:,blockNo)
        params_oce%a_tracer_v(jc,:,blockNo,2) = dvo_salt_kpp(jc,:,blockNo)
        !END DO 
      END DO
    END DO 

  END SUBROUTINE update_mixing_coefficients

!========================================================================
!========================================================================

  SUBROUTINE ensure_background_values(patch_3d, params_oce, &
                   avo_kpp, dvo_heat_kpp, dvo_salt_kpp)

    TYPE(t_patch_3d ),TARGET, INTENT(in)  :: patch_3d
    TYPE(t_patch),POINTER                 :: patch_2D
    TYPE(t_subset_range), POINTER         :: all_cells
    TYPE(t_ho_params), INTENT(inout)      :: params_oce

    ! required variables
    REAL(wp), DIMENSION(:,:,:), INTENT(INOUT) ::      &
      avo_kpp                                        ,&
      dvo_heat_kpp                                   ,&
      dvo_salt_kpp

    !REAL(wp), DIMENSION(:,:), INTENT(IN) ::           &
    !  kOBL_2D

    ! loop variables
    INTEGER :: jc, jk, blockNo
    INTEGER :: start_index, end_index
    INTEGER :: levels

    INTEGER :: cell_1_idx, cell_1_block, cell_2_idx, cell_2_block

    patch_2D => patch_3d%p_patch_2d(1)
    all_cells => patch_3d%p_patch_2d(1)%cells%ALL
    kbot => patch_3d%p_patch_1d(1)%dolic_c(:,:)

    DO blockNo = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blockNo, start_index, end_index)
      DO jc = start_index, end_index
        !levels = patch_3d%p_patch_1d(1)%dolic_c(jc, blockNo)

        ! set background values if necessary     
        DO jk = 2,kbot(jc,blockNo)
          avo_kpp(jc,jk,blockNo) = MAX( avo_kpp(jc,jk,blockNo), velocity_VerticalDiffusion_background )
          dvo_heat_kpp(jc,jk,blockNo) = MAX( dvo_heat_kpp(jc,jk,blockNo), Temperature_VerticalDiffusion_background  )
          dvo_salt_kpp(jc,jk,blockNo) = MAX( dvo_salt_kpp(jc,jk,blockNo), Salinity_VerticalDiffusion_background )
        END DO

        ! FIXME: set first layer to zero
        avo_kpp(jc,1,blockNo) = 0.0_wp
        dvo_heat_kpp(jc,1,blockNo) = 0.0_wp
        dvo_salt_kpp(jc,1,blockNo) = 0.0_wp
    
      END DO
    END DO  

  END SUBROUTINE ensure_background_values

!========================================================================
!========================================================================

  SUBROUTINE check_minmaxmean_2D(info_text,dummy2D,in_subset,str_module,idt_src)
    CHARACTER(*)        :: info_text     ! text to print
    REAL(wp)            :: dummy2D(:,:)  ! input field
    CHARACTER(LEN=12)   :: str_module    ! Output of module for 1 line debug
    INTEGER             :: idt_src       ! Level of detail for 1 line debug;
                                         ! output print level (1-5, fix)
    REAL(wp)            :: minmaxmean(3) ! array to store min/max/mean value

    TYPE(t_subset_range), POINTER :: in_subset

    ! get min/max/mean value of field
    minmaxmean(:) = global_minmaxmean(values = dummy2D(:,:),in_subset=in_subset)

    ! print min/max/mean value
    CALL debug_print_MaxMinMean(info_text, minmaxmean, str_module, idt_src)

  END SUBROUTINE check_minmaxmean_2D

  SUBROUTINE check_minmaxmean_3D(info_text,dummy3D,in_subset,str_module,idt_src)
    CHARACTER(*)        :: info_text     ! text to print
    REAL(wp)            :: dummy3D(:,:,:)  ! input field
    CHARACTER(LEN=12)   :: str_module    ! Output of module for 1 line debug
    INTEGER             :: idt_src       ! Level of detail for 1 line debug;
                                         ! output print level (1-5, fix)
    REAL(wp)            :: minmaxmean(3) ! array to store min/max/mean value

    TYPE(t_subset_range), POINTER :: in_subset

    ! get min/max/mean value of field
    minmaxmean(:) = global_minmaxmean(values = dummy3D(:,:,:),in_subset=in_subset)

    ! print min/max/mean value
    CALL debug_print_MaxMinMean(info_text, minmaxmean, str_module, idt_src)

  END SUBROUTINE check_minmaxmean_3D



END MODULE mo_ocean_cvmix_kpp

