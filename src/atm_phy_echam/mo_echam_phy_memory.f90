#if (defined (__GNUC__) || defined(__SUNPRO_F95) || defined(__SX__))
#define HAVE_F95
#endif
!>
!! Data types and variables used by the ECHAM6 physics package.
!!
!! This module contains
!! <ol>
!! <li> definition of data types for organising the physical quantities in the
!!      ECHAM physics package,
!! <li> the actual variables that are declared of these types, and
!! <li> subroutines for (de-)allocating memory for the variables.
!! </ol>
!! This module has a functionality similar to mo_memory_g3b in ECHAM,
!! but uses derived data types in order to allow for local refinement.
!!
!! @author Hui≤ Wan (MPI-M)
!! @author Marco Giorgetta (MPI-M)
!! @author Kristina Froehlich (DWD, MPI-M)
!! @author Luis Kornblueh (MPI-M)
!!
!! @par Revision History
!! First version by Hui Wan, Marco Giorgetta and Kristina Froehlich, 2010-10-28.
!! Memory allocation method changed from explicit allocation to Luis' 
!! infrastructure by Hui Wan (MPI-M, 2011-04-24)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_echam_phy_memory

  USE mo_kind,                ONLY: dp, wp
  USE mo_impl_constants,      ONLY: SUCCESS, MAX_CHAR_LENGTH,  & 
    &                               VINTP_METHOD_PRES,         &
    &                               VINTP_METHOD_LIN,          &
    &                               VINTP_METHOD_LIN_NLEVP1
  USE mo_cdi_constants,       ONLY: GRID_UNSTRUCTURED_CELL,    &
    &                               GRID_CELL
  USE mo_exception,           ONLY: message, finish
  USE mo_fortran_tools,       ONLY: t_ptr_2d, t_ptr_3d
  USE mo_parallel_config,     ONLY: nproma
  USE mo_io_config,           ONLY: lnetcdf_flt64_output
  USE mo_name_list_output_config,   ONLY: first_output_name_list, &
    &                                     is_variable_in_output
  USE mtime,                  ONLY: OPERATOR(>)
  USE mo_time_config,         ONLY: time_config
  USE mo_echam_phy_config,    ONLY: echam_phy_tc, dt_zero
  USE mo_echam_sfc_indices,   ONLY: nsfc_type, csfc
  USE mo_model_domain,        ONLY: t_patch

  USE mo_linked_list,         ONLY: t_var_list
  USE mo_var_list,            ONLY: default_var_list_settings, &
    &                               add_var, add_ref,          &
    &                               new_var_list,              &
    &                               delete_var_list
  USE mo_var_metadata,        ONLY: create_vert_interp_metadata, vintp_types, &
    &                               new_action, actions
  USE mo_action,              ONLY: ACTION_RESET
  USE mo_cf_convention,       ONLY: t_cf_var
  USE mo_grib2,               ONLY: t_grib2_var, grib2_var
  USE mo_cdi,                 ONLY: DATATYPE_PACK16, DATATYPE_PACK24,  &
    &                               DATATYPE_FLT32,  DATATYPE_FLT64,   &
    &                               GRID_UNSTRUCTURED, GRID_LONLAT,    &
    &                               TSTEP_INSTANT, TSTEP_CONSTANT,     &
    &                               TSTEP_MIN, TSTEP_MAX,              &
    &                               cdiInqMissval, DATATYPE_INT
  USE mo_zaxis_type,          ONLY: ZA_REFERENCE, ZA_REFERENCE_HALF,         &
    &                               ZA_SURFACE, ZA_GENERIC_ICE
  USE mo_sea_ice_nml,         ONLY: kice
    !ART
  USE mo_art_config,          ONLY: ctracer_art
  USE mo_run_config,          ONLY: lart

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: prm_field, prm_tend                         !< variables
  PUBLIC :: prm_field_list, prm_tend_list               !< variable lists
  PUBLIC :: construct_echam_phy_state                   !< subroutine
  PUBLIC :: destruct_echam_phy_state                    !< subroutines
  PUBLIC :: t_echam_phy_field, t_echam_phy_tend         !< derived types

  PUBLIC :: cdimissval

  CHARACTER(len=*), PARAMETER :: thismodule = 'mo_echam_phy_memory'

  !!--------------------------------------------------------------------------
  !!                               DATA TYPES
  !!--------------------------------------------------------------------------
  !>
  !! Derived data type: t_echam_phy_field
  !!
  !! This structure contains two kinds of components:
  !! <ol>
  !! <li> quantities involved in the parameterisation scheme but not in the
  !!      dynamical core, e.g., cloud cover, 10-m wind speed;
  !! <li> atmospheric state variables shared by dynamics and physics, e.g.
  !!      wind, temperature, pressure and tracer concentrations.
  !!      At each time step, the dynamical core provides initial values of
  !!      these quantites on the dynamics grid, which are then interpolated
  !!      to the physics grid and passed on to the physics package.
  !!      In the physics package, these variables may be updated once or
  !!      even more times, depending on the actual numerical schemes chosen
  !!      for the physics-dynamics coupling and the coupling between
  !!      individual parameterisation schemes;
  !! </ol>
  !!
  !! All components are arrays of one of the following shapes:
  !! <ol>
  !! <li> (nproma,           nblks_phy)
  !! <li> (nproma, nlev_phy, nblks_phy)
  !! <li> (nproma, nlev_phy, nblks_phy, ntracers)
  !! </ol>
  !! Currently the physics grid has the same spatial resolution as the
  !! dynamics grid, but is unstaggered. This means
  !!
  !!    nlev_phy = nlev
  !!   nblks_phy = patch%nblks_c
  !!
  !! In the long run, the physics grid and dynamics grid may differ in
  !! horizontal and/or vertical resolution, or even in shape.

  TYPE t_echam_phy_field

    ! Metrics
    REAL(wp),POINTER ::     &
      ! horizontal grid
      & clon      (:,:)=>NULL(),    &!< [rad]   longitude at cell center
      & clat      (:,:)=>NULL(),    &!< [rad]   longitude at cell center
      & areacella (:,:)=>NULL(),    &!< [m2]    atmosphere grid-cell area
      ! vertical grid
      & zh        (:,:,:)=>NULL(),  &!< [m]     geometric height at half levels
      & zf        (:,:,:)=>NULL(),  &!< [m]     geometric height at full levels
      & dz        (:,:,:)=>NULL()    !< [m]     geometric height thickness of layer
      
    ! Meteorology and tracers
    REAL(wp),POINTER ::     &
      & ua        (:,:,:)=>NULL(),  &!< [m/s]   zonal wind
      & va        (:,:,:)=>NULL(),  &!< [m/s]   meridional wind
      & vor       (:,:,:)=>NULL(),  &!< [1/s]   relative vorticity
      & ta        (:,:,:)=>NULL(),  &!< [K]     temperature          (tm1  of memory_g1a in ECHAM)
      & tv        (:,:,:)=>NULL(),  &!< [K]     virtual temperature  (tvm1 of memory_g1a in ECHAM)
      & qtrc      (:,:,:,:)=>NULL(),&!< [kg/kg] tracer concentration (qm1, xlm1, xim1 of memory_g1a in ECHAM)
      & mtrc      (:,:,:,:)=>NULL(),&!< [kg/m2] tracer content
      & mtrcvi    (:,:,:)=>NULL(),  &!< [kg/m2] tracer content, vertically integrated through the atmospheric column
      & prw       (:,:)=>NULL(),    &!< [kg/m2] water vapor content, vertically integrated through the atmospheric column
      & cllvi     (:,:)=>NULL(),    &!< [kg/m2] cloud water content, vertically integrated through the atmospheric column
      & clivi     (:,:)=>NULL(),    &!< [kg/m2] cloud ice   content, vertically integrated through the atmospheric column
      & rho       (:,:,:)=>NULL(),  &!< [kg/m3] air density
      & mh2o      (:,:,:)=>NULL(),  &!< [kg/m2] h2o content (vap+liq+ice)
      & mair      (:,:,:)=>NULL(),  &!< [kg/m2] air content
      & mdry      (:,:,:)=>NULL(),  &!< [kg/m2] dry air content
      & mref      (:,:,:)=>NULL(),  &!< [kg/m2] reference air content
      & xref      (:,:,:)=>NULL(),  &!< []      ratio mair/mdry
      & mh2ovi    (:,:)=>NULL(),    &!< [kg/m2] h2o content, vertically integrated through the atmospheric column
      & mairvi    (:,:)=>NULL(),    &!< [kg/m2] air content, vertically integrated through the atmospheric column
      & mdryvi    (:,:)=>NULL(),    &!< [kg/m2] dry air content, vertically integrated through the atmospheric column
      & mrefvi    (:,:)=>NULL(),    &!< [kg/m2] reference air content, vertically integrated through the atmospheric column
      & omega     (:,:,:)=>NULL(),  &!< [Pa/s]  vertical velocity in pressure coord. ("vervel" in ECHAM)
      & geoi      (:,:,:)=>NULL(),  &!< [m2/s2] geopotential above ground at half levels (vertical interfaces)
      & geom      (:,:,:)=>NULL(),  &!< [m2/s2] geopotential above ground at full levels (layer ave. or mid-point value)
      & presi_old (:,:,:)=>NULL(),  &!< [Pa]    pressure at half levels at time step "old"
      & presm_old (:,:,:)=>NULL(),  &!< [Pa]    pressure at full levels at time step "old"
      & presi_new (:,:,:)=>NULL(),  &!< [Pa]    pressure at half levels at time step "new"
      & presm_new (:,:,:)=>NULL()    !< [Pa]    pressure at full levels at time step "new"

    TYPE(t_ptr_3d),ALLOCATABLE :: qtrc_ptr(:)
    TYPE(t_ptr_3d),ALLOCATABLE :: mtrc_ptr(:)
    TYPE(t_ptr_2d),ALLOCATABLE :: mtrcvi_ptr(:)

    ! Radiation
    REAL(wp),POINTER ::       &
      !
      & cosmu0      (:,  :)=>NULL(),  &!< [ ]    cos of zenith angle mu0 for radiative heating  calculation
      & cosmu0_rt   (:,  :)=>NULL(),  &!< [ ]    cos of zenith angle mu0 at radiation time step
      & daylght_frc (:,  :)=>NULL(),  &!< [ ]    daylight fraction at each grid point
      & daylght_frc_rt(:,:)=>NULL(),  &!< [ ]    daylight fraction at each grid point for radiation time step
      !
      ! shortwave fluxes
      ! - through the atmosphere at the radiation time step
      & rsd_rt      (:,:,:)=>NULL(),  &!< [W/m2] downwelling shortwave radiation
      & rsu_rt      (:,:,:)=>NULL(),  &!< [W/m2] upwelling   shortwave radiation
      & rsdcs_rt    (:,:,:)=>NULL(),  &!< [W/m2] downwelling clear-sky shortwave radiation
      & rsucs_rt    (:,:,:)=>NULL(),  &!< [W/m2] upwelling   clear-sky shortwave radiation
      ! - at the top of the atmosphere at all times
      & rsdt        (:,  :)=>NULL(),  &!< [W/m2] toa incident shortwave radiation
      & rsut        (:,  :)=>NULL(),  &!< [W/m2] toa outgoing shortwave radiation
      & rsutcs      (:,  :)=>NULL(),  &!< [W/m2] toa outgoing clear-sky shortwave radiation
      ! - at the surface at all times
      & rsds        (:,  :)=>NULL(),  &!< [W/m2] surface downwelling shortwave radiation
      & rsus        (:,  :)=>NULL(),  &!< [W/m2] surface upwelling   shortwave radiation
      & rsdscs      (:,  :)=>NULL(),  &!< [W/m2] surface downwelling clear-sky shortwave radiation
      & rsuscs      (:,  :)=>NULL(),  &!< [W/m2] surface upwelling   clear-sky shortwave radiation
      !
      ! shortwave flux components at the surface
      ! - at radiation times
      & rvds_dir_rt (:,  :)=>NULL(),  &!< [W/m2] surface downwelling direct  visible            radiation
      & rpds_dir_rt (:,  :)=>NULL(),  &!< [W/m2] surface downwelling direct  photosynth. active radiation
      & rnds_dir_rt (:,  :)=>NULL(),  &!< [W/m2] surface downwelling direct  near-infrared      radiation
      & rvds_dif_rt (:,  :)=>NULL(),  &!< [W/m2] surface downwelling diffuse visible            radiation
      & rpds_dif_rt (:,  :)=>NULL(),  &!< [W/m2] surface downwelling diffuse photosynth. active radiation
      & rnds_dif_rt (:,  :)=>NULL(),  &!< [W/m2] surface downwelling diffuse near-infrared      radiation
      & rvus_rt     (:,  :)=>NULL(),  &!< [W/m2] surface   upwelling         visible            radiation
      & rpus_rt     (:,  :)=>NULL(),  &!< [W/m2] surface   upwelling         photosynth. active radiation
      & rnus_rt     (:,  :)=>NULL(),  &!< [W/m2] surface   upwelling         near-infrared      radiation
      ! - at all times
      & rvds_dir    (:,  :)=>NULL(),  &!< [W/m2] surface downwelling direct  visible            radiation
      & rpds_dir    (:,  :)=>NULL(),  &!< [W/m2] surface downwelling direct  photosynth. active radiation
      & rnds_dir    (:,  :)=>NULL(),  &!< [W/m2] surface downwelling direct  near-infrared      radiation
      & rvds_dif    (:,  :)=>NULL(),  &!< [W/m2] surface downwelling diffuse visible            radiation
      & rpds_dif    (:,  :)=>NULL(),  &!< [W/m2] surface downwelling diffuse photosynth. active radiation
      & rnds_dif    (:,  :)=>NULL(),  &!< [W/m2] surface downwelling diffuse near-infrared      radiation
      & rvus        (:,  :)=>NULL(),  &!< [W/m2] surface   upwelling         visible            radiation
      & rpus        (:,  :)=>NULL(),  &!< [W/m2] surface   upwelling         photosynth. active radiation
      & rnus        (:,  :)=>NULL(),  &!< [W/m2] surface   upwelling         near-infrared      radiation
      !
      ! longwave fluxes
      ! - through the atmosphere at the radiation time step
      & rld_rt      (:,:,:)=>NULL(),  &!< [W/m2] downwelling longwave radiation
      & rlu_rt      (:,:,:)=>NULL(),  &!< [W/m2] upwelling   longwave radiation
      & rldcs_rt    (:,:,:)=>NULL(),  &!< [W/m2] downwelling clear-sky longwave radiation
      & rlucs_rt    (:,:,:)=>NULL(),  &!< [W/m2] upwelling   clear-sky longwave radiation
      ! - at the top of the atmosphere at all times
      & rlut        (:,  :)=>NULL(),  &!< [W/m2] toa outgoing longwave radiation
      & rlutcs      (:,  :)=>NULL(),  &!< [W/m2] toa outgoing clear-sky longwave radiation
      ! - at the surface at all times
      & rlds        (:,  :)=>NULL(),  &!< [W/m2] surface downwelling longwave radiation
      & rlus        (:,  :)=>NULL(),  &!< [W/m2] surface upwelling   longwave radiation
      & rldscs      (:,  :)=>NULL(),  &!< [W/m2] surface downwelling clear-sky longwave radiation
      & rluscs      (:,  :)=>NULL(),  &!< [W/m2] surface downwelling clear-sky longwave radiation
      !
      & o3          (:,:,:)    !< temporary set ozone mass mixing ratio  
    ! aerosol optical properties
    REAL(wp),POINTER ::      &
      & aer_aod_533 (:,:,:)=>NULL(),  &!< aerosol optical depth at 533 nm
      & aer_ssa_533 (:,:,:)=>NULL(),  &!< aerosol single scattering albedo at 533 nm
      & aer_asy_533 (:,:,:)=>NULL(),  &!< aerosol asymmetry factor at 533 nm
      & aer_aod_2325(:,:,:)=>NULL(),  &!< aerosol optical depth at 2325 nm
      & aer_ssa_2325(:,:,:)=>NULL(),  &!< aerosol single scattering albedo at 2325 nm
      & aer_asy_2325(:,:,:)=>NULL(),  &!< aerosol asymmetry factor at 2325 nm
      & aer_aod_9731(:,:,:)=>NULL()    !< effective aerosol optical depth at 9731 nm
            !< the last quantity is in the thermal wavelength ranch, 
            !< the first lie in the solar spectrum
    ! Cloud and precipitation
    REAL(wp),POINTER ::     &
      & aclc      (:,:,:)=>NULL(),  &!< [m2/m2] cloud area fractional
      & aclcov    (:,  :)=>NULL(),  &!< [m2/m2] total cloud cover
      & acdnc     (:,:,:)=>NULL(),  &!< cloud droplet number concentration [1/m**3]
      & hur       (:,:,:)=>NULL(),  &!< relative humidity
      & rsfl      (:,  :)=>NULL(),  &!< sfc rain flux, large scale [kg m-2 s-1]
      & rsfc      (:,  :)=>NULL(),  &!< sfc rain flux, convective  [kg m-2 s-1]
      & ssfl      (:,  :)=>NULL(),  &!< sfc snow flux, large scale [kg m-2 s-1]
      & ssfc      (:,  :)=>NULL(),  &!< sfc snow flux, convective  [kg m-2 s-1]
      & pr        (:,  :)=>NULL()    !< precipitation flux         [kg m-2 s-1]

    ! Tropopause
    REAL(wp),POINTER ::     &
      & ptp       (:,  :)=>NULL()    !< tropopause air pressure [Pa]

    REAL(wp),POINTER ::     &
      & rintop (:,  :)=>NULL(),     &!< low lever inversion, computed by "cover" (memory_g3b)
      & rtype  (:,  :)=>NULL(),     &!< type of convection 0...3. (in memory_g3b in ECHAM)
      & topmax (:,  :)=>NULL()       !< maximum height of convective cloud tops [Pa] (memory_g3b)
    INTEGER ,POINTER ::     &
      & ictop  (:,  :)=>NULL()       !< level index of cnovective cloud top

    ! Vertical diffusion
    REAL(wp),POINTER ::     &
      & thvsig (:,  :)=>NULL()       !< Std. dev. of virtual potential temperature at the upper
                                     !< interface of the lowest model layer.
                                     !< Computed in "vdiff" and used by "cucall".

    REAL(wp),POINTER ::     &
      & siced  (:,  :)=>NULL(),     &!< ice depth
      & alb    (:,  :)=>NULL(),     &!< surface background albedo
      & seaice (:,  :)=>NULL()       !< sea ice as read in from amip input

    ! Energy and moisture budget related diagnostic variables
    REAL(wp),POINTER ::     &
      & cpair    (:,:,:)=>NULL(),   &!< specific heat of air at constant pressure [J/kg/K]
      & cvair    (:,:,:)=>NULL(),   &!< specific heat of air at constant volume   [J/kg/K]
      & qconv    (:,:,:)=>NULL(),   &!< convert heating to temp tend. [(K/s)/(W/m^2)]
      !
      & q_phy    (:,:,:)=>NULL(),   &!< layer heating by physics [W/m^2]
      & q_phy_vi (:,  :)=>NULL(),   &!< vertically integrated heating by physics [W/m^2]
      !
      & q_rad    (:,:,:)=>NULL(),   &!< Layer heating by LW+SW radiation
      & q_rad_vi (:,  :)=>NULL(),   &!< Vertically integrated heating by LW+SW radiation
      & q_rlw    (:,:,:)=>NULL(),   &!< Layer heating by LW radiation
      & q_rlw_vi (:,  :)=>NULL(),   &!< Vertically integrated heating by LW radiation
      & q_rsw    (:,:,:)=>NULL(),   &!< Layer heating by SW radiation
      & q_rsw_vi (:,  :)=>NULL(),   &!< Vertically integrated heating by SW radiation
      & q_vdf    (:,:,:)=>NULL(),   &!< Layer heating by vertical diffusion
      & q_vdf_vi (:,  :)=>NULL(),   &!< Vertically integrated heating by vertical diffusion
      & q_cnv    (:,:,:)=>NULL(),   &!< Layer heating by convection
      & q_cnv_vi (:,  :)=>NULL(),   &!< Vertically integrated heating by convection
      & q_cld    (:,:,:)=>NULL(),   &!< Layer heating by cloud processes
      & q_cld_vi (:,  :)=>NULL(),   &!< Vertically integrated heating by cloud processes
      & q_gwd    (:,:,:)=>NULL(),   &!< Layer heating by atmospheric gravity wave dissipation
      & q_gwd_vi (:,  :)=>NULL(),   &!< Vertically integrated heating by atmospheric gravity wave dissipation
      & q_sso    (:,:,:)=>NULL(),   &!< Layer heating by orographic gravity wave dissipation
      & q_sso_vi (:,  :)=>NULL(),   &!< Vertically integrated heating by orographic gravity wave dissipation
      !
      & sh_vdiff (:,  :)=>NULL(),   &!< sensible heat flux of vdiff
      & qv_vdiff (:,  :)=>NULL(),   &!< qv flux of vdiff
      & con_dtrl (:,  :)=>NULL(),   &!< detrainment of liquid from convection
      & con_dtri (:,  :)=>NULL(),   &!< detrainment of ice from convection 
      & con_iteqv(:,  :)=>NULL()     !< v. int. tendency of water vapor within convection

    ! orography
    REAL(wp),POINTER :: &
      & oromea (:,  :)=>NULL(),     &!< Orographic mean elevation
      & orostd (:,  :)=>NULL(),     &!< Orographic standard deviation
      & orosig (:,  :)=>NULL(),     &!< Orographic slope
      & orogam (:,  :)=>NULL(),     &!< Orographic anisotropy
      & orothe (:,  :)=>NULL(),     &!< Orographic angle
      & oropic (:,  :)=>NULL(),     &!< Orographic peacks elevation
      & oroval (:,  :)=>NULL()       !< Orographic valleys elevation

    ! JSBACH
    REAL(wp),POINTER :: &
      & ts_rad     (:,  :)=>NULL(),  &!< [K] radiative sfc. temperature for use in radiation
      & ts_rad_rt  (:,  :)=>NULL(),  &!< [K] radiative sfc. temperature at radiation time
      & csat       (:,  :)=>NULL(),  &!<
      & cair       (:,  :)=>NULL(),  &!<
      & q_snocpymlt(:,  :)=>NULL(),  &!< [W/m2] heating used to melt snow on the canopy
      & q_rlw_impl (:,  :)=>NULL(),  &!< [W/m2] heating correction due to implicit land surface coupling
      & q_rlw_nlev (:,  :)=>NULL()    !< [W/m2] heating in the lowest layer

    ! CO2
    REAL(wp),POINTER :: &
      & co2_flux_tile   (:,:,:)=>NULL(),  &!< CO2 flux on tiles (land, ocean)
      & fco2nat         (:,  :)=>NULL()    !< Surface Carbon Mass Flux into the Atmosphere Due to Natural Sources

    TYPE(t_ptr_2d),ALLOCATABLE :: co2_flux_tile_ptr(:)

    ! Sea ice.
    ! See also sea_ice/thermodyn/mo_sea_ice_types.f90
    INTEGER              :: kice  ! Number of ice-thickness classes
    REAL(wp),POINTER     ::     &
      & Tsurf   (:,:,:)=>NULL(),        &! Ice surface temperature [degC]
      & T1      (:,:,:)=>NULL(),        &! Temperature of upper ice layer [degC]
      & T2      (:,:,:)=>NULL(),        &! Temperature of lower ice layer [degC]
      & hi      (:,:,:)=>NULL(),        &! Ice thickness [m]
      & hs      (:,:,:)=>NULL(),        &! Snow thickness on ice [m]
      & Qtop    (:,:,:)=>NULL(),        &! Energy flux available for surface melting [W/m^2]
      & Qbot    (:,:,:)=>NULL(),        &! Energy flux at ice-ocean interface [W/m^2]
      & conc    (:,:,:)=>NULL(),        &! Ice concentration [0,1]
      & albvisdir_ice(:,:,:)=>NULL(),   &! Ice surface albedo for visible range, direct
      & albvisdif_ice(:,:,:)=>NULL(),   &! Ice surface albedo for visible range, diffuse
      & albnirdir_ice(:,:,:)=>NULL(),   &! Ice surface albedo for near IR range, direct
      & albnirdif_ice(:,:,:)=>NULL()     ! Ice surface albedo for near IR range, diffuse

    ! Sub grid scale orographic effects (sso)
    REAL(wp),POINTER ::         &
      & u_stress_sso   (:,:)=>NULL(),   &!< Zonal gravity wave stress
      & v_stress_sso   (:,:)=>NULL(),   &!< Meridional gravity wave stress
      & dissipation_sso(:,:)=>NULL()     !< Dissipation of orographic waves

    ! Turbulence
    REAL(wp),POINTER ::     &
      & totte       (:,:,:)=>NULL(),  &!< total turbulent energy at step n+1
      & tottem0     (:,:,:)=>NULL(),  &!< total turbulent energy at step n
      & tottem1     (:,:,:)=>NULL()    !< total turbulent energy at step n-1

    ! need only for vdiff ++++
    REAL(wp),POINTER ::     &
      & ri_atm    (:,:,:)=>NULL(),  &!< moist Richardson number at layer interfaces
      & mixlen    (:,:,:)=>NULL()    !< mixing length at layer interfaces


    REAL(wp),POINTER ::      &
      & cfm     (:,:,:)=>NULL(),     &!< turbulent exchange coefficient
      & cfm_tile(:,:,:)=>NULL(),     &!< turbulent exchange coefficient
      & cfh     (:,:,:)=>NULL(),     &!< turbulent exchange coefficient
      & cfh_tile(:,:,:)=>NULL(),     &!< turbulent exchange coefficient
      & cfv     (:,:,:)=>NULL(),     &!< turbulent exchange coefficient
      & cftotte (:,:,:)=>NULL(),     &!< turbulent exchange coefficient
      & cfthv   (:,:,:)=>NULL()       !< turbulent exchange coefficient

    TYPE(t_ptr_2d),ALLOCATABLE :: cfm_tile_ptr(:)
    TYPE(t_ptr_2d),ALLOCATABLE :: cfh_tile_ptr(:)

    REAL(wp),POINTER ::     &
      & coriol(:,:)=>NULL(),        &!< Coriolis parameter
      & hdtcbl (:,:)=>NULL(),       &!< height of the top of the atmospheric dry convective boundary layer
      & z0m_tile(:,:,:)=>NULL(),    &!< aerodynamic roughness length (over each surface type)
      & z0m   (:,:)=>NULL(),        &!< aerodynamic roughness length (grid box mean)
      & z0h_lnd(:,:)=>NULL(),       &!< roughness length for heat (over land)
      & ustar (:,:)=>NULL(),        &!<
      & wstar (:,:)=>NULL(),        &!< convective velocity scale
      & wstar_tile(:,:,:)=>NULL(),  &!< convective velocity scale (over each surface type)
      & kedisp(:,:)=>NULL(),        &!< time-mean (or integrated?) vertically integrated dissipation of kinetic energy
      & ocu   (:,:)=>NULL(),        &!< eastward  velocity of ocean surface current
      & ocv   (:,:)=>NULL()          !< northward velocity of ocean surface current

      !
    REAL(wp),POINTER ::     &
      ! net fluxes at TOA and surface
      & swflxsfc_tile(:,:,:)=>NULL(),  &!< [ W/m2] shortwave net flux at surface
      & lwflxsfc_tile(:,:,:)=>NULL(),  &!< [ W/m2] longwave net flux at surface
      & dlwflxsfc_dT (:,:)  =>NULL()       !< [ W/m2/K] longwave net flux temp tend at surface

    TYPE(t_ptr_2d),ALLOCATABLE :: swflxsfc_tile_ptr(:)
    TYPE(t_ptr_2d),ALLOCATABLE :: lwflxsfc_tile_ptr(:)

    TYPE(t_ptr_2d),ALLOCATABLE :: z0m_tile_ptr(:)
    TYPE(t_ptr_2d),ALLOCATABLE :: wstar_tile_ptr(:)

    ! need only for vdiff ----

    ! Surface variables

    REAL(wp),POINTER :: &
      & orog(:,:)=>NULL(),          &!< surface altitude [m]
      & sftlf (:,:)=>NULL(),        &!< cell area fraction occupied by land including lakes (1. = land and/or lakes only, 0. = ocean only)
      & sftgif(:,:)=>NULL(),        &!< cell area fraction occupied by land ice             (1. = land ice only, 0. = no land ice)
      & sftof (:,:)=>NULL(),        &!< cell area fraction occupied by ocean                (1. = ocean only, 0. = land and/or lakes only)
      & lsmask(:,:)=>NULL(),        &!< cell area fraction occupied by land excluding lakes (1. = land, 0. = ocean or lake only) 
      & alake (:,:)=>NULL(),        &!< cell area fraction occupied by lakes
      & glac  (:,:)=>NULL(),        &!< land area fraction that is glaciated
      & lake_ice_frc(:,:)=>NULL(),  &!< lake area fraction that is ice covered
      & icefrc(:,:)=>NULL(),        &!< ice cover given as the fraction of grid box (friac  in memory_g3b)
      & ts_tile(:,:,:)=>NULL(),     &!< surface temperature over land/water/ice
      & ts     (:,  :)=>NULL(),     &!< surface temperature, grid box mean
      & qs_sfc_tile(:,:,:)=>NULL()   !< saturation specific humidity at surface 

    TYPE(t_ptr_2d),ALLOCATABLE :: ts_tile_ptr(:)
    TYPE(t_ptr_2d),ALLOCATABLE :: qs_sfc_tile_ptr(:)

    ! Surface albedo
    REAL(wp),POINTER :: &
      & albvisdir_tile (:,:,:)=>NULL(),  &!< [ ] surface albedo over tiles for visible range, direct
      & albvisdif_tile (:,:,:)=>NULL(),  &!< [ ] surface albedo over tiles for visible range, diffuse
      & albnirdir_tile (:,:,:)=>NULL(),  &!< [ ] surface albedo over tiles for near-IR range, direct
      & albnirdif_tile (:,:,:)=>NULL(),  &!< [ ] surface albedo over tiles for near-IR range, diffuse
      & albedo_tile    (:,:,:)=>NULL(),  &!< [ ] surface albedo over tiles
      & albvisdir      (:,:  )=>NULL(),  &!< [ ] surface albedo for visible range, direct, grid-box mean
      & albvisdif      (:,:  )=>NULL(),  &!< [ ] surface albedo for visible range, diffuse, grid-box mean
      & albnirdir      (:,:  )=>NULL(),  &!< [ ] surface albedo for near-IR range, direct, grid-box mean
      & albnirdif      (:,:  )=>NULL(),  &!< [ ] surface albedo for near-IR range, diffuse, grid-box mean
      & albedo         (:,:  )=>NULL()    !< [ ] surface albedo, grid-box mean

    TYPE(t_ptr_2d),ALLOCATABLE :: albvisdir_tile_ptr(:), albvisdif_tile_ptr(:), &
      & albnirdir_tile_ptr(:), albnirdif_tile_ptr(:), albedo_tile_ptr(:)

    REAL(wp),POINTER :: &
      & lhflx     (:,  :)=>NULL(),    &!< grid box mean latent   heat flux at surface
      & shflx     (:,  :)=>NULL(),    &!< grid box mean sensible heat flux at surface
      & evap      (:,  :)=>NULL(),    &!< grid box mean evaporation at surface
      & lhflx_tile(:,:,:)=>NULL(),    &!< latent   heat flux at surface on tiles
      & shflx_tile(:,:,:)=>NULL(),    &!< sensible heat flux at surface on tiles
      & evap_tile (:,:,:)=>NULL(),    &!< evaporation at surface on tiles
      & frac_tile (:,:,:)=>NULL()      !< surface fraction of tiles:
                                       !  - fraction of land without lakes
                                       !  - fraction of ice covered water in the grid box, for sea and lakes
                                       !  - fraction of open water in the grid box, for sea and lakes

    TYPE(t_ptr_2d),ALLOCATABLE :: lhflx_tile_ptr(:)
    TYPE(t_ptr_2d),ALLOCATABLE :: shflx_tile_ptr(:)
    TYPE(t_ptr_2d),ALLOCATABLE :: evap_tile_ptr(:)
    TYPE(t_ptr_2d),ALLOCATABLE :: frac_tile_ptr(:)

    REAL(wp),POINTER :: &
      & u_stress     (:,  :)=>NULL(), &!< grid box mean wind stress
      & v_stress     (:,  :)=>NULL(), &!< grid box mean wind stress
      & u_stress_tile(:,:,:)=>NULL(), &!< wind stress on tiles
      & v_stress_tile(:,:,:)=>NULL()   !< wind stress on tiles

    TYPE(t_ptr_2d),ALLOCATABLE :: u_stress_tile_ptr(:)
    TYPE(t_ptr_2d),ALLOCATABLE :: v_stress_tile_ptr(:)

    ! Near surface diagnostics (2m temp; 2m dew point temp; 10m wind)
    !
    REAL(wp),POINTER ::        &
      & sfcwind     (:,  :)=>NULL(),   &!< grid box mean 10 m wind
      & uas         (:,  :)=>NULL(),   &!< grid box mean 10m u-velocity
      & vas         (:,  :)=>NULL(),   &!< grid box mean 10m v-velocity
      & tas         (:,  :)=>NULL(),   &!< grid box mean 2m temperature
      & dew2        (:,  :)=>NULL(),   &!< grid box mean 2m dew point temperature
      & tasmax      (:,  :)=>NULL(),   &!< grid box mean maximum 2m temperature
      & tasmin      (:,  :)=>NULL(),   &!< grid box mean minimum 2m temperature
      & sfcwind_tile(:,:,:)=>NULL(),   &!< 10 m wind on tiles
      & uas_tile    (:,:,:)=>NULL(),   &!< 10m u-velocity on tiles
      & vas_tile    (:,:,:)=>NULL(),   &!< 10m v-velocity on tiles
      & tas_tile    (:,:,:)=>NULL(),   &!< 2m temperature on tiles
      & dew2_tile   (:,:,:)=>NULL()     !< 2m dew point temperature on tiles

    TYPE(t_ptr_2d),ALLOCATABLE :: sfcwind_tile_ptr(:)
    TYPE(t_ptr_2d),ALLOCATABLE :: uas_tile_ptr(:)
    TYPE(t_ptr_2d),ALLOCATABLE :: vas_tile_ptr(:)
    TYPE(t_ptr_2d),ALLOCATABLE :: tas_tile_ptr(:)
    TYPE(t_ptr_2d),ALLOCATABLE :: dew2_tile_ptr(:)

    ! global diagnostics
    REAL(wp),POINTER ::       &
      !
      & tas_gmean    (:)=>NULL(),      &!< [K] global mean 2m-temperature
      & rsdt_gmean   (:)=>NULL(),      &!< [W/m2] global mean toa incident shortwave radiation
      & rsut_gmean   (:)=>NULL(),      &!< [W/m2] global mean toa outgoing shortwave radiation
      & rlut_gmean   (:)=>NULL(),      &!< [W/m2] global mean toa outgoing longwave radiation
      & prec_gmean   (:)=>NULL(),      &!< [kg/m2/s] global mean precipitation flux
      & evap_gmean   (:)=>NULL(),      &!< [kg/m2/s] global mean evaporation flux
      & radtop_gmean (:)=>NULL(),      &!< [W/m2] global mean toa total radiation, derived variable
      & fwfoce_gmean (:)=>NULL(),      &!< [kg/m2/s] global mean freshwater flux over ocean area, derived variable
      & icefrc_gmean (:)=>NULL()!,      &!< global mean ice cover given as the fraction of grid box
   
    ! coupling to HAMOCC lcpl_co2_atmoce
    REAL(wp),POINTER :: &
      & co2mmr(:,:) =>NULL(),  &  !< co2 mixing ratio
      & co2flux(:,:)=>NULL()      !< co2 flux

  END TYPE t_echam_phy_field

  !>
  !! Data type containing the tendencies returned by the individual parameterizations
  !!
  !! The components are arrays of one of the following shapes:
  !! <ol>
  !! <li> (nproma, nlev_phy, nblks_phy)
  !! <li> (nproma, nlev_phy, nblks_phy, ntracers)
  !! </ol>
  !!
  TYPE t_echam_phy_tend

    REAL(wp), POINTER ::   &
      !
      ! tendency due to all processes
      !
      &   ua     (:,:,:)=>NULL()  , & !< [m/s2]    u-wind
      &   va     (:,:,:)=>NULL()  , & !< [m/s2]    v-wind
      &   ta     (:,:,:)=>NULL()  , & !< [K/s]     temperature
      & qtrc     (:,:,:,:)=>NULL(), & !< [kg/kg/s] tracer mass mixing ratio
      !
      ! tendency due to resolved dynamics
      !
      &   ua_dyn (:,:,:)=>NULL()  , & !< [m/s2]    u-wind
      &   va_dyn (:,:,:)=>NULL()  , & !< [m/s2]    v-wind
      &   ta_dyn (:,:,:)=>NULL()  , & !< [K/s]     temperature (for const. volume)
      & qtrc_dyn (:,:,:,:)=>NULL(), & !< [kg/kg/s] tracer mass mixing ratio
      !
      ! tendency due to parameterized processes
      !
      &   ua_phy (:,:,:)=>NULL()  , & !< [m/s2]    u-wind
      &   va_phy (:,:,:)=>NULL()  , & !< [m/s2]    v-wind
      &   ta_phy (:,:,:)=>NULL()  , & !< [K/s]     temperature (for const. volume)
      & qtrc_phy (:,:,:,:)=>NULL(), & !< [kg/kg/s] tracer mass mixing ratio
      & mtrc_phy (:,:,:,:)=>NULL(), & !< [kg/m2/s] tracer mass
      & mtrcvi_phy(:,:,  :)=>NULL(),& !< [kg/m2/s] tracer content, vertically integrated through the atmospheric column
      !
      ! cloud microphysics
      !
      &   ta_cld (:,:,:)=>NULL()  , & !< temperature tendency due to large scale cloud processes (for const. pressure)
      & qtrc_cld (:,:,:,:)=>NULL(), & !< tracer tendency  due to large scale cloud processes
      !
      ! cumulus convection
      !
      &   ta_cnv (:,:,:)=>NULL(),   & !< temperature tendency due to convective cloud processes (for const. pressure)
      &   ua_cnv (:,:,:)=>NULL(),   & !< u-wind tendency due to convective cloud processes
      &   va_cnv (:,:,:)=>NULL(),   & !< v-wind tendency due to convective cloud processes
      & qtrc_cnv (:,:,:,:)=>NULL(), & !< tracer tendency due to convective cloud processes
      !
      ! vertical diffusion ("vdiff")
      !
      &   ta_vdf (:,:,:)=>NULL()  , & !< temperature tendency due to vertical diffusion (for const. pressure)
      &   ua_vdf (:,:,:)=>NULL()  , & !< u-wind tendency due to vertical diffusion
      &   va_vdf (:,:,:)=>NULL()  , & !< v-wind tendency due to vertical diffusion
      & qtrc_vdf (:,:,:,:)=>NULL(), & !< tracer tendency due to vertical diffusion
      !
      ! surface scheme
      !
      &   ta_sfc (:,:)=>NULL()  , & !< temperature tendency in lowermost layer due to surface processes (for const. pressure)
      !
      ! Hines param. for atmospheric gravity waves
      !
      &   ua_gwd (:,:,:)=>NULL()  , & !< u-wind tendency due to non-orographic gravity waves
      &   va_gwd (:,:,:)=>NULL()  , & !< v-wind tendency due to non-orographic gravity waves
      &   ta_gwd (:,:,:)=>NULL()  , & !< temperature tendency due to non-orographic gravity waves (for const. pressure)
      !
      ! subgrid scale orographic (sso) blocking and gravity wave drag
      !
      &   ua_sso (:,:,:)=>NULL()  , & !< u-wind tendency due to sub grid scale orography
      &   va_sso (:,:,:)=>NULL()  , & !< v-wind tendency due to sub grid scale orography
      &   ta_sso (:,:,:)=>NULL()  , & !< temperature tendency due to sub grid scale orography (for const. pressure)
      !
      ! radiation
      !
      &   ta_rsw (:,:,:)=>NULL()  , & !< temperature due to shortwave radiation (for const. pressure)
      &   ta_rlw (:,:,:)=>NULL()  , & !< temperature due to longwave  radiation  (for const. pressure)
      &   ta_rad (:,:,:)=>NULL()  , & !< temperature due to SW + LW   radiation (for const. pressure)
      &   ta_rlw_impl(:,:)=>NULL(), & !< temperature tendency due to LW rad. due to implicit land surface temperature change
      !                          (for const. pressure)
      !
      ! methane oxidation
      ! 
      & qtrc_mox (:,:,:,:)=>NULL(), & !< tracer mass mixing ratio (in fact that of water vapour)
      !                                  due to methane oxidation and H2O photolysis
      !
      ! Cariolle linearized ozone
      ! 
      & qtrc_car (:,:,:,:)=>NULL()    !< tracer mass mixing ratio (in fact that of ozone)
      !                                  due to Cariolle's linearized ozone chemistry

    TYPE(t_ptr_3d),ALLOCATABLE :: qtrc_ptr(:)
    TYPE(t_ptr_3d),ALLOCATABLE :: qtrc_dyn_ptr(:)
    TYPE(t_ptr_3d),ALLOCATABLE :: qtrc_phy_ptr(:)
    TYPE(t_ptr_3d),ALLOCATABLE :: qtrc_cld_ptr(:)
    TYPE(t_ptr_3d),ALLOCATABLE :: qtrc_cnv_ptr(:)
    TYPE(t_ptr_3d),ALLOCATABLE :: qtrc_vdf_ptr(:)
    TYPE(t_ptr_3d),ALLOCATABLE :: qtrc_mox_ptr(:)
    TYPE(t_ptr_3d),ALLOCATABLE :: qtrc_car_ptr(:)
              
    TYPE(t_ptr_3d),ALLOCATABLE :: mtrc_phy_ptr(:)
    TYPE(t_ptr_2d),ALLOCATABLE :: mtrcvi_phy_ptr(:)

  END TYPE t_echam_phy_tend

  !!--------------------------------------------------------------------------
  !!                          STATE VARIABLES 
  !!--------------------------------------------------------------------------
  !! The variable names have the prefix "prm_" in order to emphasize that they
  !! are defined for and used in parameterisations.

  TYPE(t_echam_phy_field),ALLOCATABLE,TARGET :: prm_field(:)  !< shape: (n_dom)
  TYPE(t_echam_phy_tend ),ALLOCATABLE,TARGET :: prm_tend (:)  !< shape: (n_dom)

  !!--------------------------------------------------------------------------
  !!                          VARIABLE LISTS
  !!--------------------------------------------------------------------------
  TYPE(t_var_list),ALLOCATABLE :: prm_field_list(:)  !< shape: (n_dom)
  TYPE(t_var_list),ALLOCATABLE :: prm_tend_list (:)  !< shape: (n_dom)

  REAL(dp), SAVE :: cdimissval
 
CONTAINS


  !!--------------------------------------------------------------------------
  !!                SUBROUTINES FOR BUILDING AND DELETING VARIABLE LISTS 
  !!--------------------------------------------------------------------------
  !>
  !! Top-level procedure for building the physics state
  !!
  SUBROUTINE construct_echam_phy_state( patch_array, ntracer )

    TYPE(t_patch),INTENT(IN) :: patch_array(:)
    INTEGER,INTENT(IN) :: ntracer
    CHARACTER(len=13) :: listname_f
    CHARACTER(len=12) :: listname_t
    CHARACTER(len=MAX_CHAR_LENGTH) :: ctracer(ntracer) !< tracer acronyms
    INTEGER :: ndomain, jg, ist, nblks, nlev, jtrc

    !---

    CALL message(thismodule,'Construction of ECHAM physics state started.')


    IF (lart) THEN
        ctracer = ctracer_art
    ELSE

    ! Define tracer names to be used in the construction of tracer related variables
    !
      DO jtrc = 1,ntracer
        !
        ! specific names
        SELECT CASE (jtrc)
        CASE(1)
          ctracer(jtrc) = 'hus'
        CASE(2)
          ctracer(jtrc) = 'clw'
        CASE(3)
          ctracer(jtrc) = 'cli'
        CASE(4)
          ctracer(jtrc) = 'o3'
        CASE(5)
          ctracer(jtrc) = 'co2'
        CASE(6)
          ctracer(jtrc) = 'ch4'
        CASE(7)
          ctracer(jtrc) = 'n2o'
        CASE DEFAULT
          WRITE(ctracer(jtrc),'(a1,i0)') 't',jtrc
        END SELECT
      END DO

    ENDIF

    cdimissval = cdiInqMissval()

    ! Allocate pointer arrays prm_field and prm_tend, 
    ! as well as the corresponding list arrays.

    ndomain = SIZE(patch_array)

    ALLOCATE( prm_field(ndomain), prm_tend(ndomain), STAT=ist)
    IF (ist/=SUCCESS) CALL finish(thismodule, &
      &'allocation of prm_field/tend array failed')

    ALLOCATE( prm_field_list(ndomain), prm_tend_list(ndomain), STAT=ist)
    IF (ist/=SUCCESS) CALL finish(thismodule, &
      &'allocation of prm_field/tend list array failed')

    ! Build a field list and a tendency list for each grid level.
    ! This includes memory allocation. 

    DO jg = 1,ndomain

      nblks = patch_array(jg)%nblks_c
      nlev  = patch_array(jg)%nlev

      WRITE(listname_f,'(a,i2.2)') 'prm_field_D',jg
      CALL new_echam_phy_field_list( jg, nproma, nlev, nblks, ntracer, ctracer, &
                                   & nsfc_type, listname_f, '',             &
                                   & prm_field_list(jg), prm_field(jg)          )

      WRITE(listname_t,'(a,i2.2)') 'prm_tend_D',jg
      CALL new_echam_phy_tend_list( jg, nproma, nlev, nblks, ntracer, ctracer, &
                                  & listname_t, 'tend_',                   &
                                  & prm_tend_list(jg), prm_tend(jg)            )
    ENDDO
    CALL message(thismodule,'Construction of ECHAM physics state finished.')

  END SUBROUTINE construct_echam_phy_state
  !--------------------------------------------------------------------


  !--------------------------------------------------------------------
  !>
  !! Release memory used by the state variable arrays and list arrays
  !!
  SUBROUTINE destruct_echam_phy_state

    INTEGER :: ndomain  !< total # of grid levels/domains
    INTEGER :: jg       !< grid level/domain index
    INTEGER :: ist      !< system status code

    !---
    CALL message(thismodule,'Destruction of ECHAM physics state started.')

    ndomain = SIZE(prm_field)

    DO jg = 1,ndomain
      CALL delete_var_list( prm_field_list(jg) )
      CALL delete_var_list( prm_tend_list (jg) )
    ENDDO

    DEALLOCATE( prm_field_list, prm_tend_list, STAT=ist )
    IF (ist/=SUCCESS) CALL finish(thismodule, &
      & 'deallocation of prm_field/tend list array failed')

    DEALLOCATE( prm_field, prm_tend, STAT=ist )
    IF (ist/=SUCCESS) CALL finish(thismodule, &
      & 'deallocation of prm_field/tend array failed')

    CALL message(thismodule,'Destruction of ECHAM physics state finished.')

  END SUBROUTINE destruct_echam_phy_state
  !--------------------------------------------------------------------


  !--------------------------------------------------------------------
  !>
  SUBROUTINE new_echam_phy_field_list( jg, kproma, klev, kblks, ktracer, &
                                     & ctracer, ksfc_type, listname,     &
                                     & prefix, field_list, field         )

    INTEGER,INTENT(IN) :: jg !> patch ID
    INTEGER,INTENT(IN) :: kproma, klev, kblks, ktracer, ksfc_type  !< dimension sizes

    CHARACTER(len=*)              ,INTENT(IN) :: listname, prefix
    CHARACTER(len=MAX_CHAR_LENGTH),INTENT(IN) :: ctracer(ktracer) !< tracer acronyms


    TYPE(t_var_list),       INTENT(INOUT) :: field_list
    TYPE(t_echam_phy_field),INTENT(INOUT) :: field

    ! Local variables

    CHARACTER(len=MAX_CHAR_LENGTH) :: varname
    LOGICAL :: contvar_is_in_output

    TYPE(t_cf_var)    ::    cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER :: shape2d(2), shape3d(3), shapesfc(3), shapeice(3), shape3d_layer_interfaces(3)
!0!    INTEGER :: shape4d(4)
    INTEGER :: ibits, iextbits
    INTEGER :: datatype_flt
    INTEGER :: jsfc, jtrc, tlen

    ibits = DATATYPE_PACK16
    iextbits = DATATYPE_PACK24

    IF ( lnetcdf_flt64_output ) THEN
      datatype_flt = DATATYPE_FLT64
    ELSE
      datatype_flt = DATATYPE_FLT32
    ENDIF

    shape2d  = (/kproma,       kblks/)
    shape3d  = (/kproma, klev, kblks/)
    shapesfc = (/kproma, kblks, ksfc_type/)
    shape3d_layer_interfaces = (/kproma,klev+1,kblks/)


    ! Register a field list and apply default settings

    CALL new_var_list( field_list, listname, patch_id=jg )
    CALL default_var_list_settings( field_list,                &
                                  & lrestart=.TRUE.  )

    !------------------------------
    ! Metrics
    !------------------------------

    cf_desc    = t_cf_var('cell_longitude', 'rad',                              &
                &         'cell center longitude',                              &
                &         datatype_flt)
    grib2_desc = grib2_var(0,191,2, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'clon', field%clon,                       &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,      &
                & ldims=shape2d,                                                &
                & lrestart = .FALSE.,                                           &
                & isteptype=TSTEP_CONSTANT  )

    cf_desc    = t_cf_var('cell_latitude', 'rad',                               &
                &         'cell center latitude',                               &
                &         datatype_flt)
    grib2_desc = grib2_var(0,191,1, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'clat', field%clat,                       &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,      &
                & ldims=shape2d,                                                &
                & lrestart = .FALSE.,                                           &
                & isteptype=TSTEP_CONSTANT  )

    cf_desc    = t_cf_var('cell_area', 'm2',                                    &
                &         'Atmosphere Grid-Cell Area',                          &
                &         datatype_flt)
    grib2_desc = grib2_var(255,255,255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'areacella', field%areacella,             &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,      &
                & ldims=shape2d,                                                &
                & lrestart = .FALSE.,                                           &
                & isteptype=TSTEP_CONSTANT  )

    cf_desc    = t_cf_var('geometric_height_at_half_level', 'm',                &
                &         'Geometric height at half level in physics',          &
                &         datatype_flt)
    grib2_desc = grib2_var(0, 3, 6, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'zh_phy', field%zh,                       &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,  &
                & ldims=shape3d_layer_interfaces,                               &
                & vert_interp=create_vert_interp_metadata(                      &
                &   vert_intp_type=vintp_types("P","Z","I"),                    &
                &   vert_intp_method=VINTP_METHOD_LIN_NLEVP1 ),                 &
                & lrestart = .FALSE.,                                           &
                & isteptype=TSTEP_CONSTANT  )

    cf_desc    = t_cf_var('geometric_height_at_full_level', 'm',                &
                &         'Geometric height at full level in physics',          &
                &         datatype_flt)
    grib2_desc = grib2_var(0, 3, 6, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'zf_phy', field%zf,                       &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,       &
                & ldims=shape3d,                                                &
                & vert_interp=create_vert_interp_metadata(                      &
                &   vert_intp_type=vintp_types("P","Z","I"),                    &
                &   vert_intp_method=VINTP_METHOD_LIN ),                        &
                & lrestart = .FALSE.,                                           &
                & isteptype=TSTEP_CONSTANT  )

    cf_desc    = t_cf_var('geometric_height_thickness', 'm',                    &
                &         'Geometric height thickness in physics',              &
                &         datatype_flt)
    grib2_desc = grib2_var(0, 3, 6, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'dz_phy', field%dz,                       &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,       &
                & ldims=shape3d,                                                &
                & vert_interp=create_vert_interp_metadata(                      &
                &   vert_intp_type=vintp_types("P","Z","I"),                    &
                &   vert_intp_method=VINTP_METHOD_LIN ),                        &
                & lrestart = .FALSE.,                                           &
                & isteptype=TSTEP_CONSTANT  )


    !------------------------------
    ! Meteorological quantities
    !------------------------------

    ! &       field% ua        (nproma,nlev  ,nblks),          &
    cf_desc    = t_cf_var('eastward_wind', 'm s-1', 'u-component of wind in physics', datatype_flt)
    grib2_desc = grib2_var(0, 2, 2, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'ua_phy', field%ua,                                &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc, ldims=shape3d, &
                & lrestart = .FALSE.,                                                    &
                & vert_interp = &
                &   create_vert_interp_metadata( vert_intp_type=vintp_types("P","Z","I") ) )

    ! &       field% va        (nproma,nlev  ,nblks),          &
    cf_desc    = t_cf_var('northward_wind', 'm s-1', 'v-component of wind in physics', datatype_flt)
    grib2_desc = grib2_var(0, 2, 3, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'va_phy', field%va,                                &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc, ldims=shape3d, &
                & lrestart = .FALSE.,                                                    &
                & vert_interp = &
                &   create_vert_interp_metadata( vert_intp_type=vintp_types("P","Z","I") ) )

    ! &       field% vor       (nproma,nlev  ,nblks),          &
    cf_desc    = t_cf_var('vorticity', 's-1', 'relative vorticity in physics', datatype_flt)
    grib2_desc = grib2_var(0, 2, 12, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'vor_phy', field%vor,                     &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc, ldims=shape3d, &
                & lrestart = .FALSE.,                                           &
                & vert_interp=create_vert_interp_metadata(                      &
                &   vert_intp_type=vintp_types("P","Z","I"),                    &
                &   vert_intp_method=VINTP_METHOD_LIN ) )

    ! &       field% ta        (nproma,nlev  ,nblks),          &
    cf_desc    = t_cf_var('temperature', 'K', 'temperature in physics', datatype_flt)
    grib2_desc = grib2_var(0, 0, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'ta_phy', field%ta,                       &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc, ldims=shape3d, &
                & lrestart = .FALSE.,                                           &
                & vert_interp=create_vert_interp_metadata(                      &
                &   vert_intp_type=vintp_types("P","Z","I"),                    &
                &   vert_intp_method=VINTP_METHOD_LIN ) )

    ! &       field% tv        (nproma,nlev  ,nblks),          &
    cf_desc    = t_cf_var('virtual_temperature', 'K', 'virtual temperature in physics', datatype_flt)
    grib2_desc = grib2_var(0,0,1, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'tv_phy', field%tv,                       &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc, ldims=shape3d, &
                & lrestart = .FALSE.,                                           &
                & vert_interp=create_vert_interp_metadata(                      &
                &   vert_intp_type=vintp_types("P","Z","I"),                    &
                &   vert_intp_method=VINTP_METHOD_LIN ) )

    ! OZONE 
    ! &       field% o3        (nproma,nlev  ,nblks),          &
    cf_desc    = t_cf_var('ozone', 'kg/kg', 'ozone mixing ratio', datatype_flt)
    grib2_desc = grib2_var(0,14,1, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'tro3', field%o3,                         &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc, ldims=shape3d, &
                & lrestart = .FALSE.,                                           &
                & vert_interp=create_vert_interp_metadata(                      &
                &   vert_intp_type=vintp_types("P","Z","I"),                    &
                &   vert_intp_method=VINTP_METHOD_LIN ) )

    ! aerosol optical properties
    ! at 533 nm
    cf_desc    = t_cf_var('aer_aod_533','-','aerosol optical depth at 533 nm', &
                & datatype_flt)
    grib2_desc = grib2_var(0,20,102, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'aer_aod_533', field%aer_aod_533,        &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,      &
                & ldims=shape3d,                                               &
                & lrestart = .FALSE.,                                          &
                & vert_interp=create_vert_interp_metadata(                     &
                &   vert_intp_type=vintp_types("P","Z","I"),                   &
                &   vert_intp_method=VINTP_METHOD_LIN,                         &
                &   l_extrapol=.FALSE. ) )
    cf_desc    = t_cf_var('aer_ssa_533','-',                                   &
                & 'aerosol single scattering albedo at 533 nm', datatype_flt)
    grib2_desc = grib2_var(0,20,103, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'aer_ssa_533', field%aer_ssa_533,        &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,      &
                & ldims=shape3d,                                               &
                & lrestart = .FALSE.,                                          &
                & vert_interp=create_vert_interp_metadata(                     &
                &   vert_intp_type=vintp_types("P","Z","I"),                   &
                &   vert_intp_method=VINTP_METHOD_LIN,                         &
                &   l_extrapol=.FALSE. ) )
    cf_desc    = t_cf_var('aer_asy_533','-',                                   &
                & 'aerosol asymmetry factor at 533 nm', datatype_flt)
    grib2_desc = grib2_var(0,20,104, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'aer_asy_533', field%aer_asy_533,        &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,      &
                & ldims=shape3d,                                               &
                & lrestart = .FALSE.,                                          &
                & vert_interp=create_vert_interp_metadata(                     &
                &   vert_intp_type=vintp_types("P","Z","I"),                   &
                &   vert_intp_method=VINTP_METHOD_LIN ) )
    ! at 2325 nm
    cf_desc    = t_cf_var('aer_aod_2325','-',                                  &
                & 'aerosol optical depth at 2325 nm', datatype_flt)
    grib2_desc = grib2_var(0,20,102, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'aer_aod_2325', field%aer_aod_2325,      &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,      &
                & ldims=shape3d,                                               &
                & lrestart = .FALSE.,                                          &
                & vert_interp=create_vert_interp_metadata(                     &
                &   vert_intp_type=vintp_types("P","Z","I"),                   &
                &   vert_intp_method=VINTP_METHOD_LIN,                         &
                &   l_extrapol=.FALSE. ) )
    cf_desc    = t_cf_var('aer_ssa_2325','-',                                  &
                & 'aerosol single scattering albedo at 2325 nm', datatype_flt)
    grib2_desc = grib2_var(0,20,103, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'aer_ssa_2325', field%aer_ssa_2325,      &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,      &
                & ldims=shape3d,                                               &
                & lrestart = .FALSE.,                                          &
                & vert_interp=create_vert_interp_metadata(                     &
                &   vert_intp_type=vintp_types("P","Z","I"),                   &
                &   vert_intp_method=VINTP_METHOD_LIN,                         &
                &   l_extrapol=.FALSE. ) )
    cf_desc    = t_cf_var('aer_asy_2325','-',                                  &
                & 'aerosol asymmetry factor at 2325 nm', datatype_flt)
    grib2_desc = grib2_var(0,20,104, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'aer_asy_2325', field%aer_asy_2325,      &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,      &
                & ldims=shape3d,                                               &
                & lrestart = .FALSE.,                                          &
                & vert_interp=create_vert_interp_metadata(                     &
                &   vert_intp_type=vintp_types("P","Z","I"),                   &
                &   vert_intp_method=VINTP_METHOD_LIN,                         &
                &   l_extrapol=.FALSE. ) )
    ! at 9731 nm
    cf_desc    = t_cf_var('aer_aod_9731','-',                                  &
                & 'effective aerosol optical depth at 9731 nm', datatype_flt)
    grib2_desc = grib2_var(0,20,102, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'aer_aod_9731', field%aer_aod_9731,      &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,      &
                & ldims=shape3d,                                               &
                & lrestart = .FALSE.,                                          &
                & vert_interp=create_vert_interp_metadata(                     &
                &   vert_intp_type=vintp_types("P","Z","I"),                   &
                &   vert_intp_method=VINTP_METHOD_LIN,                         &
                &   l_extrapol=.FALSE. ) )

    ! &       field% qtrc      (nproma,nlev  ,nblks,ntracer),  &
    CALL add_var( field_list, prefix//'qtrc_phy', field%qtrc,                  &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                           &
                & t_cf_var('mass_fraction_of_tracer_in_air', 'kg kg-1',        &
                &          'mass fraction of tracer in air (physics)',         &
                &          datatype_flt),                                      &
                & grib2_var(0,20,2, ibits, GRID_UNSTRUCTURED, GRID_CELL),      &
                & ldims = (/kproma,klev,kblks,ktracer/),                       &
                & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.         )

    ! &       field% mtrc      (nproma,nlev  ,nblks,ntracer),  &
    CALL add_var( field_list, prefix//'mtrc_phy', field%mtrc,                  &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                           &
                & t_cf_var('mass_of_tracer_in_air', 'kg m-2',                  &
                &          'mass of tracer in air (physics)',                  &
                &          datatype_flt),                                      &
                & grib2_var(0,20,2, ibits, GRID_UNSTRUCTURED, GRID_CELL),      &
                & ldims = (/kproma,klev,kblks,ktracer/),                       &
                & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.         )

    ! &       field% mtrcvi      (nproma,nblks,ntracer),  &
    CALL add_var( field_list, prefix//'mtrcvi_phy', field%mtrcvi,              &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                          &
                & t_cf_var('atmosphere_tracer_content', 'kg m-2',              &
                &          'tracer path (physics)',                            &
                &          datatype_flt),                                      &
                & grib2_var(0,20,2, ibits, GRID_UNSTRUCTURED, GRID_CELL),      &
                & ldims = (/kproma,kblks,ktracer/),                            &
                & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.         )

    ALLOCATE(field%qtrc_ptr(ktracer))
    ALLOCATE(field%mtrc_ptr(ktracer))
    ALLOCATE(field%mtrcvi_ptr(ktracer))
    
    ! Generic references to single tracers
    DO jtrc = 1,ktracer
      tlen = LEN_TRIM(ctracer(jtrc))
      CALL add_ref( field_list, prefix//'qtrc_phy',                            &
                  & prefix//'q'//ctracer(jtrc)(1:tlen)//'_phy', field%qtrc_ptr(jtrc)%p, &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                         &
                  & t_cf_var('mass_fraction_of_'//ctracer(jtrc)(1:tlen)//'_in_air', &
                  &          'kg kg-1',                                        &
                  &          'mass fraction of '//ctracer(jtrc)(1:tlen)//' in air (physics)', &
                  &          datatype_flt),                                    &
                  & grib2_var(0,20,2, ibits, GRID_UNSTRUCTURED, GRID_CELL),    &
                  & ref_idx=jtrc, ldims=(/kproma,klev,kblks/),                 &
                  & lrestart = .FALSE.,                                        &
                  & vert_interp=create_vert_interp_metadata(                   &
                  &             vert_intp_type=vintp_types("P","Z","I"),       & 
                  &             vert_intp_method=VINTP_METHOD_LIN,             &
                  &             l_loglin=.FALSE.,                              &
                  &             l_extrapol=.TRUE., l_pd_limit=.FALSE.,         &
                  &             lower_limit=0._wp )                            )
      CALL add_ref( field_list, prefix//'mtrc_phy',                            &
                  & prefix//'m'//ctracer(jtrc)(1:tlen)//'_phy', field%mtrc_ptr(jtrc)%p, &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                         &
                  & t_cf_var('mass_of_'//ctracer(jtrc)(1:tlen)//'_in_air',       &
                  &          'kg m-2',                                         &
                  &          'mass of '//ctracer(jtrc)(1:tlen)//' in air (physics)', &
                  &          datatype_flt),                                    &
                  & grib2_var(0,20,2, ibits, GRID_UNSTRUCTURED, GRID_CELL),    &
                  & ref_idx=jtrc, ldims=(/kproma,klev,kblks/),                 &
                  & lrestart = .FALSE.,                                        &
                  & vert_interp=create_vert_interp_metadata(                   &
                  &             vert_intp_type=vintp_types("P","Z","I"),       & 
                  &             vert_intp_method=VINTP_METHOD_LIN,             &
                  &             l_loglin=.FALSE.,                              &
                  &             l_extrapol=.TRUE., l_pd_limit=.FALSE.,         &
                  &             lower_limit=0._wp )                            )
      CALL add_ref( field_list, prefix//'mtrcvi_phy',                          &
                  & prefix//'m'//ctracer(jtrc)(1:tlen)//'vi_phy', field%mtrcvi_ptr(jtrc)%p, &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                        &
                  & t_cf_var('atmosphere_'//ctracer(jtrc)(1:tlen)//'_content',   &
                  &          'kg m-2', ctracer(jtrc)(1:tlen)//' path (physics)', &
                  &          datatype_flt),                                    &
                  & grib2_var(0,20,2, ibits, GRID_UNSTRUCTURED, GRID_CELL),    &
                  & ref_idx=jtrc, ldims=(/kproma,kblks/),                      &
                  & lrestart = .FALSE.                                         )
    END DO                                                                                

    ! Specific references for tracers with specific names
    DO jtrc = 1,ktracer
      tlen = LEN_TRIM(ctracer(jtrc))
       IF ( ctracer(jtrc)(1:tlen) == 'hus' ) THEN
          CALL add_ref( field_list, prefix//'mtrcvi_phy',                      &
                  & prefix//'prw', field%prw,                                  &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                        &
                  & t_cf_var('total_vapour', 'kg m-2', 'vertically integrated water vapour', &
                  &          datatype_flt),                                    &
                  & grib2_var(0,1,64, ibits, GRID_UNSTRUCTURED, GRID_CELL),    &
                  & ref_idx=jtrc, ldims=(/kproma,kblks/),                      &
                  & lrestart = .FALSE.                                         )
       END IF
       IF ( ctracer(jtrc)(1:tlen) == 'clw' ) THEN
          CALL add_ref( field_list, prefix//'mtrcvi_phy',                      &
                  & prefix//'cllvi', field%cllvi,                              &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                        &
                  & t_cf_var('total_cloud_water', 'kg m-2', 'vertically integrated cloud water', &
                  &          datatype_flt),                                    &
                  & grib2_var(0,1,69, ibits, GRID_UNSTRUCTURED, GRID_CELL),    &
                  & ref_idx=jtrc, ldims=(/kproma,kblks/),                      &
                  & lrestart = .FALSE.                                         )
       END IF
       IF ( ctracer(jtrc)(1:tlen) == 'cli' ) THEN
          CALL add_ref( field_list, prefix//'mtrcvi_phy',                      &
                  & prefix//'clivi', field%clivi,                              &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                        &
                  & t_cf_var('total_cloud_ice', 'kg m-2', 'vertically integrated cloud ice', &
                  &          datatype_flt),                                    &
                  & grib2_var(0,1,70, ibits, GRID_UNSTRUCTURED, GRID_CELL),    &
                  & ref_idx=jtrc, ldims=(/kproma,kblks/),                      &
                  & lrestart = .FALSE.                                         )
       END IF
    END DO

    ! &       field% rho        (nproma,nlev  ,nblks),          &
    cf_desc    = t_cf_var('air_density', 'kg m-3', 'density of air',           &
         &                datatype_flt)
    grib2_desc = grib2_var(0,3,10, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'rho_phy', field%rho,                    &
         &        GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,      &
         &        ldims=shape3d, lrestart = .FALSE.,                           &
         &        vert_interp=create_vert_interp_metadata(                     &
         &                    vert_intp_type=vintp_types("P","Z","I"),         & 
         &                    vert_intp_method=VINTP_METHOD_LIN,               &
         &                    l_loglin=.FALSE.,                                &
         &                    l_extrapol=.TRUE., l_pd_limit=.FALSE.,           &
         &                    lower_limit=0._wp  ) )

    ! &       field% mh2o        (nproma,nlev  ,nblks),          &
    cf_desc    = t_cf_var('h2o_mass', 'kg m-2', 'h2o (vap+liq+ice) mass in layer', &
         &                datatype_flt)
    grib2_desc = grib2_var(0,1,21, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'mh2o_phy', field%mh2o,                  &
         &        GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,      &
         &        ldims=shape3d, lrestart = .FALSE.,                           &
         &        vert_interp=create_vert_interp_metadata(                     &
         &                    vert_intp_type=vintp_types("P","Z","I"),         & 
         &                    vert_intp_method=VINTP_METHOD_LIN,               &
         &                    l_loglin=.FALSE.,                                &
         &                    l_extrapol=.TRUE., l_pd_limit=.FALSE.,           &
         &                    lower_limit=0._wp  ) )

    ! &       field% mh2ovi     (nproma,nblks),          &
    cf_desc    = t_cf_var('atmosphere_h2o_content', 'kg m-2', 'h2o (vap+liq+ice) path (physics)', &
         &                datatype_flt)
    grib2_desc = grib2_var(0,1,64, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'mh2ovi_phy', field%mh2ovi,              &
         &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                          &
         &        cf_desc, grib2_desc,                                         &
         &        ldims=shape2d,                                               &
         &        lrestart = .FALSE.,                                          &
         &        isteptype=TSTEP_INSTANT )

    ! &       field% mair        (nproma,nlev  ,nblks),          &
    cf_desc    = t_cf_var('air_mass', 'kg m-2', 'air mass in layer', &
         &                datatype_flt)
    grib2_desc = grib2_var(0,1,21, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'mair_phy', field%mair,                  &
         &        GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,      &
         &        ldims=shape3d, lrestart = .FALSE., initval=1.0_wp,           &
         &        vert_interp=create_vert_interp_metadata(                     &
         &                    vert_intp_type=vintp_types("P","Z","I"),         & 
         &                    vert_intp_method=VINTP_METHOD_LIN,               &
         &                    l_loglin=.FALSE.,                                &
         &                    l_extrapol=.TRUE., l_pd_limit=.FALSE.,           &
         &                    lower_limit=0._wp  ) )

    ! &       field% mairvi     (nproma,nblks),          &
    cf_desc    = t_cf_var('atmosphere_air_content', 'kg m-2', 'air path (physics)', &
         &                datatype_flt)
    grib2_desc = grib2_var(0,1,64, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'mairvi_phy', field%mairvi,              &
         &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                          &
         &        cf_desc, grib2_desc,                                         &
         &        ldims=shape2d,                                               &
         &        lrestart = .FALSE.,                                          &
         &        isteptype=TSTEP_INSTANT )

    ! &       field% mdry        (nproma,nlev  ,nblks),          &
    cf_desc    = t_cf_var('dry_air_mass', 'kg m-2', 'dry air mass in layer', &
         &                datatype_flt)
    grib2_desc = grib2_var(0,1,21, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'mdry_phy', field%mdry,                  &
         &        GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,      &
         &        ldims=shape3d, lrestart = .FALSE., initval=1.0_wp,           &
         &        vert_interp=create_vert_interp_metadata(                     &
         &                    vert_intp_type=vintp_types("P","Z","I"),         & 
         &                    vert_intp_method=VINTP_METHOD_LIN,               &
         &                    l_loglin=.FALSE.,                                &
         &                    l_extrapol=.TRUE., l_pd_limit=.FALSE.,           &
         &                    lower_limit=0._wp  ) )

    ! &       field% mdryvi     (nproma,nblks),          &
    cf_desc    = t_cf_var('atmosphere_dry_air_content', 'kg m-2', 'dry air path (physics)', &
         &                datatype_flt)
    grib2_desc = grib2_var(0,1,64, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'mdryvi_phy', field%mdryvi,              &
         &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                          &
         &        cf_desc, grib2_desc,                                         &
         &        ldims=shape2d,                                               &
         &        lrestart = .FALSE.,                                          &
         &        isteptype=TSTEP_INSTANT )

    ! &       field% mref        (nproma,nlev  ,nblks),          &
    cf_desc    = t_cf_var('ref_air_mass', 'kg m-2', 'ref air mass in layer', &
         &                datatype_flt)
    grib2_desc = grib2_var(0,1,21, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'mref_phy', field%mref,                  &
         &        GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,      &
         &        ldims=shape3d, lrestart = .FALSE., initval=1.0_wp,           &
         &        vert_interp=create_vert_interp_metadata(                     &
         &                    vert_intp_type=vintp_types("P","Z","I"),         & 
         &                    vert_intp_method=VINTP_METHOD_LIN,               &
         &                    l_loglin=.FALSE.,                                &
         &                    l_extrapol=.TRUE., l_pd_limit=.FALSE.,           &
         &                    lower_limit=0._wp  ) )

    ! &       field% mrefvi     (nproma,nblks),          &
    cf_desc    = t_cf_var('atmosphere_ref_air_content', 'kg m-2', 'ref air path (physics)', &
         &                datatype_flt)
    grib2_desc = grib2_var(0,1,64, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'mrefvi_phy', field%mrefvi,              &
         &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                          &
         &        cf_desc, grib2_desc,                                         &
         &        ldims=shape2d,                                               &
         &        lrestart = .FALSE.,                                          &
         &        isteptype=TSTEP_INSTANT )

    ! &       field% xref        (nproma,nlev  ,nblks),          &
    cf_desc    = t_cf_var('ratio_mair_mdry', '', 'ratio mair/mdry', &
         &                datatype_flt)
    grib2_desc = grib2_var(0,1,255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'xref_phy', field%xref,                  &
         &        GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,      &
         &        ldims=shape3d, lrestart = .FALSE., initval=1.0_wp,           &
         &        vert_interp=create_vert_interp_metadata(                     &
         &                    vert_intp_type=vintp_types("P","Z","I"),         & 
         &                    vert_intp_method=VINTP_METHOD_LIN,               &
         &                    l_loglin=.FALSE.,                                &
         &                    l_extrapol=.TRUE., l_pd_limit=.FALSE.,           &
         &                    lower_limit=0._wp  ) )

    ! &       field% omega     (nproma,nlev  ,nblks),          &
    cf_desc    = t_cf_var('vertical_velocity', 'Pa s-1', 'vertical velocity in physics', datatype_flt)
    grib2_desc = grib2_var(0,2,8, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'omega_phy', field%omega,                &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,      &
                & ldims=shape3d, lrestart = .FALSE.,                           &
                & vert_interp=create_vert_interp_metadata(                     &
                &             vert_intp_type=vintp_types("P","Z","I"),         &
                &             vert_intp_method=VINTP_METHOD_LIN,               &
                &             l_loglin=.FALSE., l_extrapol=.FALSE.) )

    ! &       field% geom      (nproma,nlev  ,nblks),          &
    cf_desc    = t_cf_var('geopotential', 'm2 s-2', 'geopotential above surface', datatype_flt)
    grib2_desc = grib2_var(0, 3, 4, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'gpsm', field%geom,                      &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,      &
                & ldims=shape3d, lrestart = .FALSE.,                           &
                & isteptype=TSTEP_CONSTANT,                                    &
                & vert_interp=create_vert_interp_metadata(                     &
                &             vert_intp_type=vintp_types("P","Z","I"),         &
                &             vert_intp_method=VINTP_METHOD_LIN,               &
                &             l_extrapol=.TRUE., l_pd_limit=.FALSE.) )

    ! &       field% presm_old (nproma,nlev  ,nblks),          &
    cf_desc    = t_cf_var('pressure', 'Pa', 'pressure at old time step', datatype_flt)
    grib2_desc = grib2_var(0, 3, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'pam_old', field%presm_old,              &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,      &
                & ldims=shape3d, lrestart = .FALSE.,                           &
                & vert_interp=create_vert_interp_metadata(                     &
                &             vert_intp_type=vintp_types("Z","I"),             &
                &             vert_intp_method=VINTP_METHOD_PRES ) )

    ! &       field% presm_new (nproma,nlev  ,nblks),          &
    cf_desc    = t_cf_var('pressure', 'Pa', 'pressure at new time step', datatype_flt)
    grib2_desc = grib2_var(0, 3, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'pam_new', field%presm_new,              &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,      &
                & ldims=shape3d, lrestart = .FALSE.,                           &
                & vert_interp=create_vert_interp_metadata(                     &
                &             vert_intp_type=vintp_types("Z","I"),             &
                &             vert_intp_method=VINTP_METHOD_PRES ) )


    !-- Variables defined at layer interfaces --
    ! &       field% geoi      (nproma,nlevp1,nblks),          &
    cf_desc    = t_cf_var('geopotential', 'm2 s-2', 'geopotential above surface', datatype_flt)
    grib2_desc = grib2_var(0, 3, 4, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'gpsi', field%geoi,                        &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,   &
                & ldims=shape3d_layer_interfaces, lrestart = .FALSE.,            &
                & isteptype=TSTEP_CONSTANT,                                      &
                & vert_interp=create_vert_interp_metadata(                       &
                &   vert_intp_type=vintp_types("P","Z","I"),                     &
                &   vert_intp_method=VINTP_METHOD_LIN_NLEVP1 ) )

    ! &       field% presi_old (nproma,nlevp1,nblks),          &
    cf_desc    = t_cf_var('pressure', 'Pa', 'pressure at old time step', datatype_flt)
    grib2_desc = grib2_var(0, 3, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'pai_old', field%presi_old,               &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,  &
                & ldims=shape3d_layer_interfaces, lrestart = .FALSE.,           &
                & vert_interp=create_vert_interp_metadata(                      &
                &   vert_intp_type=vintp_types("Z","I"),                        &
                &   vert_intp_method=VINTP_METHOD_LIN_NLEVP1 ) )

    ! &       field% presi_new (nproma,nlevp1,nblks),          &
    cf_desc    = t_cf_var('pressure', 'Pa', 'pressure at new time step', datatype_flt)
    grib2_desc = grib2_var(0, 3, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'pai_new', field%presi_new,               &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,  &
                & ldims=shape3d_layer_interfaces, lrestart = .FALSE.,           &
                & vert_interp=create_vert_interp_metadata(                      &
                &   vert_intp_type=vintp_types("Z","I"),                        &
                &   vert_intp_method=VINTP_METHOD_LIN_NLEVP1 ) )

    !------------------
    ! Radiation
    !------------------
    !

    cf_desc    = t_cf_var( 'cosmu0'                                      , &
         &                 ''                                            , &
         &                 'cosine of the zenith angle for rad. heating' , &
         &                 datatype_flt                                  )
    grib2_desc = grib2_var(192,214,1, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(field_list, prefix//'cosmu0' , field%cosmu0, &
         &       GRID_UNSTRUCTURED_CELL , ZA_SURFACE        , &
         &       cf_desc , grib2_desc                       , &
         &       lrestart = .FALSE.                         , &
         &       ldims=shape2d                              )

    cf_desc    = t_cf_var( 'cosmu0_rt'                                    , &
         &                 ''                                             , &
         &                 'cosine of the zenith angle at radiation time' , &
         &                 datatype_flt                                   )
    grib2_desc = grib2_var(192,214,1, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(field_list, prefix//'cosmu0_rt', field%cosmu0_rt, &
         &       GRID_UNSTRUCTURED_CELL , ZA_SURFACE             , &
         &       cf_desc , grib2_desc                            , &
         &       lrestart = .TRUE.                               , &
         &       ldims=shape2d                                   )

    cf_desc    = t_cf_var( 'daylght_frc', &
         &                 ''           , &
         &                 ''           , &
         &                 datatype_flt )
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(field_list, prefix//'daylght_frc', field%daylght_frc, &
         &       GRID_UNSTRUCTURED_CELL , ZA_SURFACE                 , &
         &       cf_desc, grib2_desc                                 , &
         &       lrestart = .FALSE.                                  , &
         &       ldims=shape2d                                       )

    cf_desc    = t_cf_var( 'daylght_frc_rt' , &
         &                 ''               , &
         &                 ''               , &
         &                 datatype_flt     )
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(field_list, prefix//'daylght_frc_rt', field%daylght_frc_rt, &
         &       GRID_UNSTRUCTURED_CELL , ZA_SURFACE                       , &
         &       cf_desc, grib2_desc                                       , &
         &       lrestart = .TRUE.                                         , &
         &       ldims=shape2d                                             )

    IF ( echam_phy_tc(jg)%dt_rad > dt_zero ) THEN
       !
       ! shortwave fluxes
       ! - through the atmosphere
       !
       cf_desc    = t_cf_var('downwelling_shortwave_flux_in_air', &
            &                'W m-2'                            , &
            &                'downwelling shortwave radiation'  , &
            &                datatype_flt                       )
       grib2_desc = grib2_var(0,4,7, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var(field_list, prefix//'rsd', field%rsd_rt     , &
            &       GRID_UNSTRUCTURED_CELL   , ZA_REFERENCE_HALF   , &
            &       cf_desc, grib2_desc                         , &
            &       lrestart = .TRUE.                           , &
            &       ldims=shape3d_layer_interfaces              , &
            &       vert_interp=create_vert_interp_metadata       &
            &         (vert_intp_type=vintp_types("P","Z","I") ,  &
            &          vert_intp_method=VINTP_METHOD_LIN_NLEVP1))

       cf_desc    = t_cf_var('upwelling_shortwave_flux_in_air', &
            &                'W m-2'                          , &
            &                'upwelling shortwave radiation'  , &
            &                datatype_flt                     )
       grib2_desc = grib2_var(0,4,8, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var(field_list, prefix//'rsu' , field%rsu_rt    , &
            &       GRID_UNSTRUCTURED_CELL    , ZA_REFERENCE_HALF  , &
            &       cf_desc, grib2_desc                         , &
            &       lrestart = .TRUE.                           , &
            &       ldims=shape3d_layer_interfaces              , &
            &       vert_interp=create_vert_interp_metadata       &
            &         (vert_intp_type=vintp_types("P","Z","I") ,  &
            &          vert_intp_method=VINTP_METHOD_LIN_NLEVP1))

       cf_desc    = t_cf_var('downwelling_shortwave_flux_in_air_assuming_clear_sky', &
            &                'W m-2'                                               , &
            &                'downwelling clear-sky shortwave radiation'           , &
            &                datatype_flt                                          )
       grib2_desc = grib2_var(0,4,7, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var(field_list, prefix//'rsdcs' , field%rsdcs_rt , &
            &       GRID_UNSTRUCTURED_CELL    , ZA_REFERENCE_HALF   , &
            &       cf_desc, grib2_desc                          , &
            &       lrestart = .TRUE.                            , &
            &       ldims=shape3d_layer_interfaces               , &
            &       vert_interp=create_vert_interp_metadata        &
            &         (vert_intp_type=vintp_types("P","Z","I") ,   &
            &          vert_intp_method=VINTP_METHOD_LIN_NLEVP1) )

       cf_desc    = t_cf_var('upwelling_shortwave_flux_in_air_assuming_clear_sky', &
            &                'W m-2'                                             , &
            &                'upwelling clear-sky shortwave radiation'           , &
            &                datatype_flt                                        )
       grib2_desc = grib2_var(0,4,8, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var(field_list, prefix//'rsucs' , field%rsucs_rt , &
            &       GRID_UNSTRUCTURED_CELL    , ZA_REFERENCE_HALF   , &
            &       cf_desc, grib2_desc                          , &
            &       lrestart = .TRUE.                            , &
            &       ldims=shape3d_layer_interfaces               , &
            &       vert_interp=create_vert_interp_metadata        &
            &         (vert_intp_type=vintp_types("P","Z","I") ,   &
            &          vert_intp_method=VINTP_METHOD_LIN_NLEVP1) )

       ! - at the top of the atmosphere
       !
       cf_desc    = t_cf_var('toa_incoming_shortwave_flux'     , &
            &                'W m-2'                           , &
            &                'toa incident shortwave radiation', &
            &                datatype_flt                      )
       grib2_desc = grib2_var(0,4,7, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var(field_list, prefix//'rsdt', field%rsdt, &
            &       GRID_UNSTRUCTURED_CELL    , ZA_SURFACE, &
            &       cf_desc, grib2_desc                   , &
            &       lrestart = .FALSE.                    , &
            &       ldims=shape2d                         )

       cf_desc    = t_cf_var('toa_outgoing_shortwave_flux'     , &
            &                'W m-2'                           , &
            &                'toa outgoing shortwave radiation', &
            &                datatype_flt                      )
       grib2_desc = grib2_var(0,4,8, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var(field_list, prefix//'rsut', field%rsut, &
            &       GRID_UNSTRUCTURED_CELL    , ZA_SURFACE, &
            &       cf_desc, grib2_desc                   , &
            &       lrestart = .FALSE.                    , &
            &       ldims=shape2d                         )

       cf_desc    = t_cf_var('toa_outgoing_shortwave_flux_assuming_clear_sky', &
            &                'W m-2'                                         , &
            &                'toa outgoing clear-sky shortwave radiation'    , &
            &                datatype_flt                                    )
       grib2_desc = grib2_var(0,4,8, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var(field_list, prefix//'rsutcs', field%rsutcs, &
            &       GRID_UNSTRUCTURED_CELL      , ZA_SURFACE  , &
            &       cf_desc, grib2_desc                       , &
            &       lrestart = .FALSE.                        , &
            &       ldims=shape2d                             )

    END IF

    ! - at the surface (also used in update surface)
    !
    cf_desc    = t_cf_var('surface_downwelling_shortwave_flux_in_air', &
         &                'W m-2'                                    , &
         &                'surface downwelling shortwave radiation'  , &
         &                datatype_flt                               )
    grib2_desc = grib2_var(0,4,7, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(field_list, prefix//'rsds', field%rsds, &
         &       GRID_UNSTRUCTURED_CELL    , ZA_SURFACE, &
         &       cf_desc, grib2_desc                   , &
         &       lrestart = .FALSE.                    , &
         &       ldims=shape2d                         )

    cf_desc    = t_cf_var('surface_upwelling_shortwave_flux_in_air', &
         &                'W m-2'                                  , &
         &                'surface upwelling shortwave radiation'  , &
         &                datatype_flt                             )
    grib2_desc = grib2_var(0,4,8, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(field_list, prefix//'rsus', field%rsus, &
         &       GRID_UNSTRUCTURED_CELL    , ZA_SURFACE, &
         &       cf_desc, grib2_desc                   , &
         &       lrestart = .FALSE.                    , &
         &       ldims=shape2d                         )


    IF ( echam_phy_tc(jg)%dt_rad > dt_zero ) THEN
       !
       cf_desc    = t_cf_var('surface_downwelling_shortwave_flux_in_air_assuming_clear_sky', &
            &                'W m-2'                                                       , &
            &                'surface downwelling clear-sky shortwave radiation'           , &
            &                datatype_flt                                                  )
       grib2_desc = grib2_var(0,4,7, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var(field_list, prefix//'rsdscs', field%rsdscs, &
            &       GRID_UNSTRUCTURED_CELL      , ZA_SURFACE  , &
            &       cf_desc, grib2_desc                       , &
            &       lrestart = .FALSE.                        , &
            &       ldims=shape2d                             )

       cf_desc    = t_cf_var('surface_upwelling_shortwave_flux_in_air_assuming_clear_sky', &
            &                'W m-2'                                                     , &
            &                'surface upwelling clear-sky shortwave radiation'           , &
            &                datatype_flt                                                )
       grib2_desc = grib2_var(0,4,8, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var(field_list, prefix//'rsuscs', field%rsuscs, &
            &       GRID_UNSTRUCTURED_CELL      , ZA_SURFACE  , &
            &       cf_desc, grib2_desc                       , &
            &       lrestart = .FALSE.                        , &
            &       ldims=shape2d                             )

       !-----------------------------------------------------------------------------------
       ! shortwave flux components at the surface
       ! - at radiation times
       cf_desc    = t_cf_var('surface_downwelling_direct_visible_flux_in_air_at_rad_time'    , &
            &                'W m-2'                                                         , &
            &                'surface downwelling direct visible radiation at radiation time', &
            &                datatype_flt                                                    )
       grib2_desc = grib2_var(0,4,7, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var(field_list, prefix//'rvds_dir_rt', field%rvds_dir_rt, &
            &       GRID_UNSTRUCTURED_CELL           , ZA_SURFACE       , &
            &       cf_desc, grib2_desc                                 , &
            &       lrestart = .TRUE.                                   , &
            &       ldims=shape2d                                       )

       cf_desc    = t_cf_var('surface_downwelling_direct_par_flux_in_air_at_rad_time'                          , &
            &                'W m-2'                                                                           , &
            &                'surface downwelling direct photosynthetically active radiation at radiation time', &
            &                datatype_flt                                                                      )
       grib2_desc = grib2_var(0,4,7, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var(field_list, prefix//'rpds_dir_rt', field%rpds_dir_rt, &
            &       GRID_UNSTRUCTURED_CELL           , ZA_SURFACE       , &
            &       cf_desc, grib2_desc                                 , &
            &       lrestart = .TRUE.                                   , &
            &       ldims=shape2d                                       )

       cf_desc    = t_cf_var('surface_downwelling_direct_nearir_flux_in_air_at_rad_time'           , &
            &                'W m-2'                                                               , &
            &                'surface downwelling direct near infrared radiation at radiation time', &
            &                datatype_flt                                                          )
       grib2_desc = grib2_var(0,4,7, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var(field_list, prefix//'rnds_dir_rt', field%rnds_dir_rt, &
            &       GRID_UNSTRUCTURED_CELL           , ZA_SURFACE       , &
            &       cf_desc, grib2_desc                                 , &
            &       lrestart = .TRUE.                                   , &
            &       ldims=shape2d                                       )


       cf_desc    = t_cf_var('surface_downwelling_diffuse_visible_flux_in_air_at_rad_time'    , &
            &                'W m-2'                                                          , &
            &                'surface downwelling diffuse visible radiation at radiation time', &
            &                datatype_flt                                                     )
       grib2_desc = grib2_var(0,4,7, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var(field_list, prefix//'rvds_dif_rt', field%rvds_dif_rt, &
            &       GRID_UNSTRUCTURED_CELL        , ZA_SURFACE          , &
            &       cf_desc, grib2_desc                                 , &
            &       lrestart = .TRUE.                                   , &
            &       ldims=shape2d                                       )

       cf_desc    = t_cf_var('surface_downwelling_diffuse_par_flux_in_air_at_rad_time'                          , &
            &                'W m-2'                                                                            , &
            &                'surface downwelling diffuse photosynthetically active radiation at radiation time', &
            &                datatype_flt                                                                       )
       grib2_desc = grib2_var(0,4,7, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var(field_list, prefix//'rpds_dif_rt', field%rpds_dif_rt, &
            &       GRID_UNSTRUCTURED_CELL           , ZA_SURFACE       , &
            &       cf_desc, grib2_desc                                 , &
            &       lrestart = .TRUE.                                   , &
            &       ldims=shape2d                                       )

       cf_desc    = t_cf_var('surface_downwelling_diffuse_nearir_flux_in_air_at_rad_time'           , &
            &                'W m-2'                                                                , &
            &                'surface downwelling diffuse near infrared radiation at radiation time', &
            &                datatype_flt                                                           )
       grib2_desc = grib2_var(0,4,7, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var(field_list, prefix//'rnds_dif_rt', field%rnds_dif_rt, &
            &       GRID_UNSTRUCTURED_CELL           , ZA_SURFACE       , &
            &       cf_desc, grib2_desc                                 , &
            &       lrestart = .TRUE.                                   , &
            &       ldims=shape2d                                       )


       cf_desc    = t_cf_var('surface_upwelling_visible_flux_in_air_at_rad_time'    , &
            &                'W m-2'                                                , &
            &                'surface upwelling visible radiation at radiation time', &
            &                datatype_flt                                           )
       grib2_desc = grib2_var(0,4,7, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var(field_list, prefix//'rvus_rt', field%rvus_rt, &
            &       GRID_UNSTRUCTURED_CELL       , ZA_SURFACE   , &
            &       cf_desc, grib2_desc                         , &
            &       lrestart = .TRUE.                           , &
            &       ldims=shape2d                               )

       cf_desc    = t_cf_var('surface_upwelling_par_flux_in_air_at_rad_time'                          , &
            &                'W m-2'                                                                  , &
            &                'surface upwelling photosynthetically active radiation at radiation time', &
            &                datatype_flt                                                             )
       grib2_desc = grib2_var(0,4,7, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var(field_list, prefix//'rpus_rt', field%rpus_rt, &
            &       GRID_UNSTRUCTURED_CELL       , ZA_SURFACE   , &
            &       cf_desc, grib2_desc                         , &
            &       lrestart = .TRUE.                           , &
            &       ldims=shape2d                               )

       cf_desc    = t_cf_var('surface_upwelling_nearir_flux_in_air_at_rad_time'           , &
            &                'W m-2'                                                      , &
            &                'surface upwelling near infrared radiation at radiation time', &
            &                datatype_flt                                                 )
       grib2_desc = grib2_var(0,4,7, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var(field_list, prefix//'rnus_rt', field%rnus_rt, &
            &       GRID_UNSTRUCTURED_CELL       , ZA_SURFACE   , &
            &       cf_desc, grib2_desc                         , &
            &       lrestart = .TRUE.                           , &
            &       ldims=shape2d                               )

    END IF

    !-----------------------------------------------------------------------------------------
    ! shortwave flux components at the surface (also used in update_surface)
    ! - at all times
    cf_desc    = t_cf_var('surface_downwelling_direct_visible_flux_in_air', &
         &                'W m-2'                                         , &
         &                'surface downwelling direct visible radiation'  , &
         &                datatype_flt                                    )
    grib2_desc = grib2_var(0,4,7, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(field_list, prefix//'rvds_dir', field%rvds_dir, &
         &       GRID_UNSTRUCTURED_CELL        , ZA_SURFACE    , &
         &       cf_desc, grib2_desc                           , &
         &       lrestart = .FALSE.                            , &
         &       ldims=shape2d                                 )

    cf_desc    = t_cf_var('surface_downwelling_direct_par_flux_in_air'                    , &
         &                'W m-2'                                                         , &
         &                'surface downwelling direct photosynthetically active radiation', &
         &                datatype_flt                                                    )
    grib2_desc = grib2_var(0,4,7, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(field_list, prefix//'rpds_dir', field%rpds_dir, &
         &       GRID_UNSTRUCTURED_CELL        , ZA_SURFACE    , &
         &       cf_desc, grib2_desc                           , &
         &       lrestart = .FALSE.                            , &
         &       ldims=shape2d                                 )

    cf_desc    = t_cf_var('surface_downwelling_direct_nearir_flux_in_air'     , &
         &                'W m-2'                                             , &
         &                'surface downwelling direct near infrared radiation', &
         &                datatype_flt                                        )
    grib2_desc = grib2_var(0,4,7, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(field_list, prefix//'rnds_dir', field%rnds_dir, &
         &       GRID_UNSTRUCTURED_CELL        , ZA_SURFACE    , &
         &       cf_desc, grib2_desc                           , &
         &       lrestart = .FALSE.                            , &
         &       ldims=shape2d                                 )


    cf_desc    = t_cf_var('surface_downwelling_diffuse_visible_flux_in_air', &
         &                'W m-2'                                          , &
         &                'surface downwelling diffuse visible radiation'  , &
         &                datatype_flt                                     )
    grib2_desc = grib2_var(0,4,7, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(field_list, prefix//'rvds_dif', field%rvds_dif, &
         &       GRID_UNSTRUCTURED_CELL        , ZA_SURFACE    , &
         &       cf_desc, grib2_desc                           , &
         &       lrestart = .FALSE.                            , &
         &       ldims=shape2d                                 )

    cf_desc    = t_cf_var('surface_downwelling_diffuse_par_flux_in_air'                    , &
         &                'W m-2'                                                          , &
         &                'surface downwelling diffuse photosynthetically active radiation', &
         &                datatype_flt                                                     )
    grib2_desc = grib2_var(0,4,7, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(field_list, prefix//'rpds_dif', field%rpds_dif, &
         &       GRID_UNSTRUCTURED_CELL        , ZA_SURFACE    , &
         &       cf_desc, grib2_desc                           , &
         &       lrestart = .FALSE.                            , &
         &       ldims=shape2d                                 )

    cf_desc    = t_cf_var('surface_downwelling_diffuse_nearir_flux_in_air'     , &
         &                'W m-2'                                              , &
         &                'surface downwelling diffuse near infrared radiation', &
         &                datatype_flt                                         )
    grib2_desc = grib2_var(0,4,7, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(field_list, prefix//'rnds_dif', field%rnds_dif, &
         &       GRID_UNSTRUCTURED_CELL        , ZA_SURFACE    , &
         &       cf_desc, grib2_desc                           , &
         &       lrestart = .FALSE.                            , &
         &       ldims=shape2d                                 )


    IF ( echam_phy_tc(jg)%dt_rad > dt_zero ) THEN
       !
       cf_desc    = t_cf_var('surface_upwelling_visible_flux_in_air', &
            &                'W m-2'                                , &
            &                'surface upwelling visible radiation'  , &
            &                datatype_flt                           )
       grib2_desc = grib2_var(0,4,7, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var(field_list, prefix//'rvus', field%rvus, &
            &       GRID_UNSTRUCTURED_CELL    , ZA_SURFACE, &
            &       cf_desc, grib2_desc                   , &
            &       lrestart = .FALSE.                    , &
            &       ldims=shape2d                         )

       cf_desc    = t_cf_var('surface_upwelling_par_flux_in_air'                    , &
            &                'W m-2'                                                , &
            &                'surface upwelling photosynthetically active radiation', &
            &                datatype_flt                                           )
       grib2_desc = grib2_var(0,4,7, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var(field_list, prefix//'rpus', field%rpus, &
            &       GRID_UNSTRUCTURED_CELL    , ZA_SURFACE, &
            &       cf_desc, grib2_desc                   , &
            &       lrestart = .FALSE.                    , &
            &       ldims=shape2d                         )

       cf_desc    = t_cf_var('surface_upwelling_nearir_flux_in_air'     , &
            &                'W m-2'                                    , &
            &                'surface upwelling near infrared radiation', &
            &                datatype_flt                               )
       grib2_desc = grib2_var(0,4,7, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var(field_list, prefix//'rnus', field%rnus, &
            &       GRID_UNSTRUCTURED_CELL    , ZA_SURFACE, &
            &       cf_desc, grib2_desc                   , &
            &       lrestart = .FALSE.                    , &
            &       ldims=shape2d                         )
       !---------------------------------------------------------

    END IF
    !
    !------------------
    !

    ! longwave  fluxes
    ! - through the atmosphere
    !
    ! (rld_rt and rlu_rt are also needed in update_surface)

    cf_desc    = t_cf_var('downwelling_longwave_flux_in_air', &
         &                'W m-2'                           , &
         &                'downwelling longwave radiation'  , &
         &                datatype_flt                       )
    grib2_desc = grib2_var(0,5,3, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(field_list, prefix//'rld' , field%rld_rt     , &
         &       GRID_UNSTRUCTURED_CELL    , ZA_REFERENCE_HALF   , &
         &       cf_desc, grib2_desc                          , &
         &       lrestart = .TRUE.                            , &
         &       ldims=shape3d_layer_interfaces               , &
         &       vert_interp=create_vert_interp_metadata        &
         &         (vert_intp_type=vintp_types("P","Z","I") ,   &
         &          vert_intp_method=VINTP_METHOD_LIN_NLEVP1) )

    cf_desc    = t_cf_var('upwelling_longwave_flux_in_air', &
         &                'W m-2'                         , &
         &                'upwelling longwave radiation'  , &
         &                datatype_flt                     )
    grib2_desc = grib2_var(0,5,4, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(field_list, prefix//'rlu' , field%rlu_rt     , &
         &       GRID_UNSTRUCTURED_CELL    , ZA_REFERENCE_HALF   , &
         &       cf_desc, grib2_desc                          , &
         &       lrestart = .TRUE.                            , &
         &       ldims=shape3d_layer_interfaces               , &
         &       vert_interp=create_vert_interp_metadata        &
         &         (vert_intp_type=vintp_types("P","Z","I") ,   &
         &          vert_intp_method=VINTP_METHOD_LIN_NLEVP1) )

    IF ( echam_phy_tc(jg)%dt_rad > dt_zero ) THEN
       !
       cf_desc    = t_cf_var('downwelling_longwave_flux_in_air_assuming_clear_sky', &
            &                'W m-2'                                              , &
            &                'downwelling clear-sky longwave radiation'           , &
            &                datatype_flt                                         )
       grib2_desc = grib2_var(0,5,3, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var(field_list, prefix//'rldcs' , field%rldcs_rt , &
            &       GRID_UNSTRUCTURED_CELL    , ZA_REFERENCE_HALF   , &
            &       cf_desc, grib2_desc                          , &
            &       lrestart = .TRUE.                            , &
            &       ldims=shape3d_layer_interfaces               , &
            &       vert_interp=create_vert_interp_metadata        &
            &         (vert_intp_type=vintp_types("P","Z","I") ,   &
            &          vert_intp_method=VINTP_METHOD_LIN_NLEVP1) )

       cf_desc    = t_cf_var('upwelling_longwave_flux_in_air_assuming_clear_sky', &
            &                'W m-2'                                            , &
            &                'upwelling clear-sky longwave radiation'           , &
            &                datatype_flt                                       )
       grib2_desc = grib2_var(0,5,4, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var(field_list, prefix//'rlucs' , field%rlucs_rt , &
            &       GRID_UNSTRUCTURED_CELL    , ZA_REFERENCE_HALF   , &
            &       cf_desc, grib2_desc                          , &
            &       lrestart = .TRUE.                            , &
            &       ldims=shape3d_layer_interfaces               , &
            &       vert_interp=create_vert_interp_metadata        &
            &         (vert_intp_type=vintp_types("P","Z","I") ,   &
            &          vert_intp_method=VINTP_METHOD_LIN_NLEVP1) )

       ! - at the top of the atmosphere
       !
       cf_desc    = t_cf_var('toa_outgoing_longwave_flux'     , &
            &                'W m-2'                          , &
            &                'toa outgoing longwave radiation', &
            &                datatype_flt                     )
       grib2_desc = grib2_var(0,5,4, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var(field_list, prefix//'rlut', field%rlut, &
            &       GRID_UNSTRUCTURED_CELL    , ZA_SURFACE, &
            &       cf_desc, grib2_desc                   , &
            &       lrestart = .FALSE.                    , &
            &       ldims=shape2d                         )

       cf_desc    = t_cf_var('toa_outgoing_longwave_flux_assuming_clear_sky', &
            &                'W m-2'                                        , &
            &                'toa outgoing clear-sky longwave radiation'    , &
            &                datatype_flt                                   )
       grib2_desc = grib2_var(0,5,4, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var(field_list, prefix//'rlutcs', field%rlutcs, &
            &       GRID_UNSTRUCTURED_CELL      , ZA_SURFACE  , &
            &       cf_desc, grib2_desc                       , &
            &       lrestart = .FALSE.                        , &
            &       ldims=shape2d                             )

    END IF

    ! - at the surface (also used in update_surface)
    !
    cf_desc    = t_cf_var('surface_downwelling_longwave_flux_in_air', &
         &                'W m-2'                                   , &
         &                'surface downwelling longwave radiation'  , &
         &                datatype_flt                               )
    grib2_desc = grib2_var(0,5,3, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(field_list, prefix//'rlds', field%rlds, &
         &       GRID_UNSTRUCTURED_CELL    , ZA_SURFACE, &
         &       cf_desc, grib2_desc                   , &
         &       lrestart = .FALSE.                    , &
         &       ldims=shape2d                         )

    cf_desc    = t_cf_var('surface_upwelling_longwave_flux_in_air', &
         &                'W m-2'                                 , &
         &                'surface upwelling longwave radiation'  , &
         &                datatype_flt                             )
    grib2_desc = grib2_var(0,5,4, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(field_list, prefix//'rlus', field%rlus, &
         &       GRID_UNSTRUCTURED_CELL    , ZA_SURFACE, &
         &       cf_desc, grib2_desc                   , &
         &       lrestart = .FALSE.                    , &
         &       ldims=shape2d                         )

    IF ( echam_phy_tc(jg)%dt_rad > dt_zero ) THEN
       !
       cf_desc    = t_cf_var('surface_downwelling_longwave_flux_in_air_assuming_clear_sky', &
            &                'W m-2'                                                      , &
            &                'surface downwelling clear-sky longwave radiation'           , &
            &                datatype_flt                                                 )
       grib2_desc = grib2_var(0,5,3, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var(field_list, prefix//'rldscs', field%rldscs, &
            &       GRID_UNSTRUCTURED_CELL      , ZA_SURFACE  , &
            &       cf_desc, grib2_desc                       , &
            &       lrestart = .FALSE.                        , &
            &       ldims=shape2d                             )

       cf_desc    = t_cf_var('surface_upwelling_longwave_flux_in_air_assuming_clear_sky', &
            &                'W m-2'                                                    , &
            &                'surface upwelling clear-sky longwave radiation'           , &
            &                datatype_flt                                               )
       grib2_desc = grib2_var(0,5,4, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var(field_list, prefix//'rluscs', field%rluscs, &
            &       GRID_UNSTRUCTURED_CELL      , ZA_SURFACE  , &
            &       cf_desc, grib2_desc                       , &
            &       lrestart = .FALSE.                        , &
            &       ldims=shape2d                             )

       !
    END IF

    !
    !------------------
    !


    cf_desc    = t_cf_var('drlns_dT', 'W m-2 K-1', 'longwave net flux T-derivative at surface', &
         &                datatype_flt)
    grib2_desc = grib2_var(0, 5, 5, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'drlns_dT', field%dlwflxsfc_dT, &
         &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                 &
         &        cf_desc, grib2_desc,                                &
         &        ldims=shape2d,                                      &
         &        lrestart = .FALSE.,                                 &
         &        isteptype=TSTEP_INSTANT )

    cf_desc    = t_cf_var('siced', 'm', 'sea ice thickness', datatype_flt)
    grib2_desc = grib2_var(10,2,1, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'sit', field%siced,                  &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
                & lrestart = .TRUE., ldims=shape2d )

    cf_desc    = t_cf_var('alb', '', 'surface albedo from external file', datatype_flt)
    grib2_desc = grib2_var(0,19,1, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'alb', field%alb,                    &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
                & lrestart = .FALSE., ldims=shape2d )

    cf_desc    = t_cf_var('ts_rad', 'K', 'radiative surface temperature', datatype_flt)
    grib2_desc = grib2_var(0,0,17, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'ts_rad', field%ts_rad,              &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
                & lrestart = .TRUE., ldims=shape2d )

    cf_desc    = t_cf_var('ts_rad_rt', 'K', 'radiative surface temperature at rad. time', datatype_flt)
    grib2_desc = grib2_var(0,0,17, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'ts_rad_rt', field%ts_rad_rt,        &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
                & lrestart = .TRUE., ldims=shape2d )

    IF (echam_phy_tc(jg)%dt_vdf > time_config%tc_dt_dyn(jg) .OR.                            &
      & is_variable_in_output(first_output_name_list, var_name=prefix//'q_snocpymlt')) THEN
       cf_desc    = t_cf_var('q_snocpymlt', 'W/m2', 'heating for snow melt on canopy', datatype_flt)
       grib2_desc = grib2_var(255,255,255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var( field_list, prefix//'q_snocpymlt', field%q_snocpymlt,    &
                   & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
                   & lrestart = .FALSE., ldims=shape2d )
    END IF

    IF (is_variable_in_output(first_output_name_list, var_name=prefix//'q_rlw_impl')) THEN
       cf_desc    = t_cf_var('q_rlw_impl', 'W/m2', 'heating correction due to implicit land surface coupling', datatype_flt)
       grib2_desc = grib2_var(255,255,255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var( field_list, prefix//'q_rlw_impl', field%q_rlw_impl,    &
                   & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
                   & lrestart = .FALSE., ldims=shape2d )
    END IF

    cf_desc    = t_cf_var('q_rlw_nlev', 'W/m2', 'LW heating in the lowest layer', datatype_flt)
    grib2_desc = grib2_var(255,255,255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'q_rlw_nlev', field%q_rlw_nlev,    &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
                & lrestart = .FALSE., ldims=shape2d )
    !
    !------------------
    !
    ! CO2

    cf_desc = t_cf_var('fco2nat', 'kg m-2 s-1',                                &
                & 'Surface Carbon Mass Flux into the Atmosphere Due to Natural Sources', datatype_flt)
    grib2_desc = grib2_var(255,255,255, ibits, GRID_UNSTRUCTURED, GRID_CELL)

    CALL add_var( field_list, prefix//'fco2nat', field%fco2nat,                &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,     &
                & lrestart = .TRUE., initval =  0.0_wp, ldims=shape2d )

    ! &       field% co2_flux_tile(nproma,nblks,nsfc_type), &
    CALL add_var( field_list, prefix//'co2_flux_tile', field%co2_flux_tile,         &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                          &
                & t_cf_var('co2_flux_tile',  'kg m-2 s-1',                     &
                & 'surface_upward_mass_flux_of_carbon_dioxide', datatype_flt), &
                & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED,GRID_CELL),&
                & ldims=shapesfc, initval=0.0_wp,                              &
                & lcontainer=.TRUE., lrestart=.FALSE.      )

    ALLOCATE(field%co2_flux_tile_ptr(ksfc_type))

    DO jsfc = 1,ksfc_type

      CALL add_ref( field_list, prefix//'co2_flux_tile',                         &
                  & prefix//'co2_flux_'//csfc(jsfc), field%co2_flux_tile_ptr(jsfc)%p,      &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                          &
                  & t_cf_var('co2_flux_'//csfc(jsfc), 'kg m-2 s-1',              &
                  & 'surface_upward_mass_flux_of_carbon_dioxide', datatype_flt), &
                  & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED,GRID_CELL),&
                  & lrestart=.TRUE., ldims=shape2d,  initval=0.0_wp,             &
                  & lmiss=.TRUE., missval=cdimissval )

    END DO

    !
    !------------------
    !

    
    ! Parameterized topography
    !
    IF ( echam_phy_tc(jg)%dt_sso > dt_zero ) THEN
       !
       cf_desc    = t_cf_var('surface_height_above_sea_level', 'm',   &
                   &         'Mean height of orography above sea level', datatype_flt)
       grib2_desc = grib2_var(0,3,6, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var( field_list, prefix//'oromea', field%oromea,              &
                   & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
                   & lrestart = .FALSE., ldims=shape2d,                       &
                   & isteptype=TSTEP_CONSTANT )
       !
       cf_desc    = t_cf_var('standard_deviation_of_height', 'm',     &
                   &         'Standard deviation of height above sea level of sub-grid scale orography', &
                   &         datatype_flt)
       grib2_desc = grib2_var(0,3,20, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var( field_list, prefix//'orostd', field%orostd,              &
                   & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
                   & lrestart = .FALSE., ldims=shape2d,                       &
                   & isteptype=TSTEP_CONSTANT )
       !
       cf_desc    = t_cf_var('slope_of_terrain', '-',                 &
                   &         'Slope of sub-gridscale orography', datatype_flt)
       grib2_desc = grib2_var(0,3,22, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var( field_list, prefix//'orosig', field%orosig,      &
                   & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
                   & lrestart = .FALSE., ldims=shape2d,                       &
                   & isteptype=TSTEP_CONSTANT )
       !
       cf_desc    = t_cf_var('anisotropy_factor', '-',                &
                   &         'Anisotropy of sub-gridscale orography', datatype_flt)
       grib2_desc = grib2_var(0,3,24, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var( field_list, prefix//'orogam', field%orogam,              &
                   & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
                   & lrestart = .FALSE., ldims=shape2d,                       &
                   & isteptype=TSTEP_CONSTANT )
       !
       cf_desc    = t_cf_var('angle_of_principal_axis', 'radians',    &
                   &         'Angle of sub-gridscale orography', datatype_flt)
       grib2_desc = grib2_var(0,3,21, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var( field_list, prefix//'orothe', field%orothe,              &
                   & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
                   & lrestart = .FALSE., ldims=shape2d,                       &
                   & isteptype=TSTEP_CONSTANT )
       !
       cf_desc    = t_cf_var('height_of_peaks', 'm', 'Height above sea level of peaks', datatype_flt)
       grib2_desc = grib2_var(0,3,6, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var( field_list, prefix//'oropic', field%oropic,              &
                   & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
                   & lrestart = .FALSE., ldims=shape2d,                       &
                   & isteptype=TSTEP_CONSTANT )
       !
       cf_desc    = t_cf_var('height_of_valleys', 'm', 'Height above sea level of valleys', datatype_flt)
       grib2_desc = grib2_var(0,3,6, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var( field_list, prefix//'oroval', field%oroval,              &
                   & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
                   & lrestart = .FALSE., ldims=shape2d,                       &
                   & isteptype=TSTEP_CONSTANT )
       !
    END IF
    !
    !------------------
    !

    cf_desc    = t_cf_var('csat', '', '', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'csat', field%csat,                  &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
                & lrestart = .TRUE., initval =  1.0_wp, ldims=shape2d )

    cf_desc    = t_cf_var('cair', '', '', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'cair', field%cair,                  &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
                & lrestart = .TRUE., initval =  1.0_wp, ldims=shape2d )

    !-------------------------
    ! Sea ice
    !-------------------------

    field%kice = kice ! Number of thickness classes - always 1, as of yet
    shapeice = (/kproma, field%kice, kblks/)

    CALL add_var( field_list, prefix//'ts_icecl', field%Tsurf ,               &
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE,                      &
      &          t_cf_var('ts_icecl', 'C','surface temperature',datatype_flt),&
      &          grib2_var(10,2,8, ibits, GRID_UNSTRUCTURED, GRID_CELL),      &
      &          ldims=shapeice, lrestart=.TRUE.)
    CALL add_var( field_list, prefix//'t1_icecl', field%T1 ,                  &
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE,                      &
      &          t_cf_var('t1_icecl','C','Temperature upper layer',datatype_flt), &
      &          grib2_var(10,2,8, ibits, GRID_UNSTRUCTURED, GRID_CELL),      &
      &          ldims=shapeice, lrestart=.TRUE.)
    CALL add_var( field_list, prefix//'t2_icecl', field%T2 ,                  &
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE,                      &
      &          t_cf_var('t2_icecl','C','Temperature lower layer', datatype_flt),&
      &          grib2_var(10,2,8, ibits, GRID_UNSTRUCTURED, GRID_CELL),      &
      &          ldims=shapeice, lrestart=.TRUE.)
    CALL add_var( field_list, prefix//'sit_icecl', field%hi ,                 &
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE,                      &
      &          t_cf_var('sit_icecl', 'm', 'ice thickness', datatype_flt),   &
      &          grib2_var(10,2,1, ibits, GRID_UNSTRUCTURED, GRID_CELL),      &
      &          ldims=shapeice, lrestart=.TRUE.)
    CALL add_var( field_list, prefix//'hs_icecl', field%hs ,                  &
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE,                      &
      &          t_cf_var('hs_icecl', 'm', 'snow thickness', datatype_flt),   &
      &          grib2_var(10,2,255, ibits, GRID_UNSTRUCTURED, GRID_CELL),    &
      &          ldims=shapeice, lrestart=.TRUE.)
    CALL add_var( field_list, prefix//'qtop_icecl', field%Qtop ,              &
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE,                      &
      &          t_cf_var('qtop_icecl', 'W/m^2', 'Energy flux available for surface melting', &
      &                   datatype_flt),                                      &
      &          grib2_var(10,2,255, ibits, GRID_UNSTRUCTURED, GRID_CELL),    &
      &          ldims=shapeice, lrestart=.FALSE.)
    CALL add_var( field_list, prefix//'qbot_icecl', field%Qbot ,              &
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE,                      &
      &          t_cf_var('qbot_icecl', 'W/m^2', 'Energy flux at ice-ocean interface', datatype_flt),&
      &          grib2_var(10,2,255, ibits, GRID_UNSTRUCTURED, GRID_CELL),    &
      &          ldims=shapeice, lrestart=.FALSE.)


    CALL add_var( field_list, prefix//'sic_icecl', field%conc ,               &
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE,                      &
      &          t_cf_var('sic_icecl', '', 'ice concentration in each ice class', datatype_flt),&
      &          grib2_var(10,2,0, ibits, GRID_UNSTRUCTURED, GRID_CELL),      &
      &          ldims=shapeice, lrestart=.TRUE.)

    ! &       field% albvisdir_ice (nproma,field%kice,nblks),          &
    cf_desc    = t_cf_var('albvisdir_icecl', '', 'ice albedo VIS direct', datatype_flt)
    grib2_desc = grib2_var(192,128,15, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'albvisdir_icecl', field%albvisdir_ice,  &
                & GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, cf_desc, grib2_desc, &
                & ldims=shapeice, lrestart=.TRUE. )

    ! &       field% albvisdif_ice (nproma,field%kice,nblks),          &
    cf_desc    = t_cf_var('albvisdif_icecl', '', 'ice albedo VIS diffuse', datatype_flt)
    grib2_desc = grib2_var(0,19,222, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'albvisdif_icecl', field%albvisdif_ice,  &
                & GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, cf_desc, grib2_desc, &
         & ldims=shapeice, lrestart=.TRUE. )

    ! &       field% albnirdir_ice (nproma,field%kice,nblks),          &
    cf_desc    = t_cf_var('albnirdir_icecl', '', 'ice albedo NIR direct', datatype_flt)
    grib2_desc = grib2_var(192,128,17, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'albnirdir_icecl', field%albnirdir_ice,  &
                & GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, cf_desc, grib2_desc, &
                & ldims=shapeice, lrestart=.TRUE. )

    ! &       field% albnirdif_ice (nproma,field%kice,nblks),          &
    cf_desc    = t_cf_var('albnirdif_icecl', '', 'ice albedo NIR diffuse', datatype_flt)
    grib2_desc = grib2_var(0, 19, 223, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'albnirdif_icecl', field%albnirdif_ice,  &
                & GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, cf_desc, grib2_desc, &
                & ldims=shapeice, lrestart=.TRUE. )


    !-------------------------
    ! Cloud and precipitation
    !-------------------------
    cf_desc    = t_cf_var('cl', 'm2 m-2', 'cloud area fraction', datatype_flt)
    grib2_desc = grib2_var(0,6,22, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'cl', field%aclc,                                  &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc, ldims=shape3d, &
                & lrestart = .FALSE.,                                                    &
                & vert_interp=create_vert_interp_metadata(                               &
                &             vert_intp_type=vintp_types("P","Z","I"),                   &
                &             vert_intp_method=VINTP_METHOD_LIN,                         &
                &             l_loglin=.FALSE.,                                          &
                &             l_extrapol=.FALSE., l_pd_limit=.FALSE.,                    &
                &             lower_limit=0._wp ) )

    cf_desc    = t_cf_var('clt', 'm2 m-2', &
               & 'total cloud cover', datatype_flt)
    grib2_desc = grib2_var(0,6,1, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'clt', field%aclcov,       &
         &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE,            &
         &        cf_desc, grib2_desc,                           &
         &        ldims=shape2d,                                 &
         &        lrestart = .FALSE.,                            &
         &        isteptype=TSTEP_INSTANT )

    ! &       field% acdnc  (nproma,nlev  ,nblks), &
    cf_desc    = t_cf_var('acdnc', 'm-3', 'cloud droplet number concentration', datatype_flt)
    grib2_desc = grib2_var(0,6,28, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'acdnc', field%acdnc,                              &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc, ldims=shape3d, &
                & lrestart = .FALSE.,                                                    &
                & vert_interp=create_vert_interp_metadata(                               &
                &             vert_intp_type=vintp_types("P","Z","I"),                   &
                &             vert_intp_method=VINTP_METHOD_LIN,                         &
                &             l_loglin=.FALSE.,                                          &
                &             l_extrapol=.TRUE., l_pd_limit=.FALSE.,                     &
                &             lower_limit=0._wp ) )

    IF (is_variable_in_output(first_output_name_list, var_name=prefix//'hur')) THEN
       cf_desc    = t_cf_var('relative_humidity', '', 'relative humidity', datatype_flt)
       grib2_desc = grib2_var(0, 1, 1, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var( field_list, prefix//'hur', field%hur   ,                               &
                   & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc, ldims=shape3d, &
                   & lrestart = .FALSE.,                                                    &
                   & vert_interp=create_vert_interp_metadata(                               &
                   &             vert_intp_type=vintp_types("P","Z","I"),                   &
                   &             vert_intp_method=VINTP_METHOD_LIN,                         &
                   &             l_loglin=.FALSE.,                                          &
                   &             l_extrapol=.FALSE., l_pd_limit=.TRUE.,                     &
                   &             lower_limit=0._wp ) )
    END IF

    cf_desc    = t_cf_var('prlr', 'kg m-2 s-1',    &
               & 'large-scale precipitation flux (water)', datatype_flt)
    grib2_desc = grib2_var(0,1,77, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'prlr', field%rsfl,        &
         &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE,            &
         &        cf_desc, grib2_desc,                           &
         &        ldims=shape2d,                                 &
         &        lrestart = .TRUE.,                             &
         &        isteptype=TSTEP_INSTANT )

    cf_desc    = t_cf_var('prcr', 'kg m-2 s-1',    &
               & 'convective precipitation flux (water)', datatype_flt)
    grib2_desc = grib2_var(0,1,76, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'prcr', field%rsfc,        &
         &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE,            &
         &        cf_desc, grib2_desc,                           &
         &        ldims=shape2d,                                 &
         &        lrestart = .TRUE.,                             &
         &        isteptype=TSTEP_INSTANT )

    cf_desc    = t_cf_var('prls', 'kg m-2 s-1',    &
               & 'large-scale precipitation flux (snow)', datatype_flt)
    grib2_desc = grib2_var(0,1,59, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'prls', field%ssfl,        &
         &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE,            &
         &        cf_desc, grib2_desc,                           &
         &        ldims=shape2d,                                 &
         &        lrestart = .TRUE.,                             &
         &        isteptype=TSTEP_INSTANT )

    cf_desc    = t_cf_var('prcs', 'kg m-2 s-1',    &
               & 'convective precipitation flux (snow)', datatype_flt)
    grib2_desc = grib2_var(0,1,58, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'prcs', field%ssfc,        &
         &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE,            &
         &        cf_desc, grib2_desc,                           &
         &        ldims=shape2d,                                 &
         &        lrestart = .TRUE.,                             &
         &        isteptype=TSTEP_INSTANT )

    cf_desc    = t_cf_var('pr', 'kg m-2 s-1',                    &
         &                'precipitation flux',                  &
         &                datatype_flt)
    grib2_desc = grib2_var(0, 1, 52, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'pr', field%pr,            &
         &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE,            &
         &        cf_desc, grib2_desc,                           &
         &        ldims=shape2d,                                 &
         &        lrestart = .FALSE.,                            &
         &        isteptype=TSTEP_INSTANT )

    ! &       field% rintop (nproma,       nblks), &
    cf_desc    = t_cf_var('rintop', '', '', datatype_flt)
    grib2_desc = grib2_var(0,6,255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'rintop', field%rintop,              &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
                & lrestart = .FALSE., ldims=shape2d )

    ! &       field% rtype  (nproma,       nblks), &
    cf_desc    = t_cf_var('convection_type', '', 'convection_type (0...3)', datatype_flt)
    grib2_desc = grib2_var(0,6,255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'rtype', field%rtype,                &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
                & lrestart = .TRUE., ldims=shape2d )

    IF (is_variable_in_output(first_output_name_list, var_name=prefix//'topmax')) THEN
       cf_desc    = t_cf_var('topmax', 'Pa', 'maximum height of convective cloud tops', &
            &                datatype_flt)
       grib2_desc = grib2_var(0,6,255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var( field_list, prefix//'topmax', field%topmax,              &
                   & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
                   & lrestart = .FALSE., ldims=shape2d,                       &
                   & initval =  99999.0_wp, resetval =  99999.0_wp,           &
                   & isteptype=TSTEP_MIN,                                     &
                   & action_list=actions(new_action(ACTION_RESET,"P1D"))      )
    END IF

    ! &       field% ictop  (nproma,       nblks), &
    cf_desc    = t_cf_var('ictop', '-', 'level index of convective cloud tops', &
         &                datatype_int)
    grib2_desc = grib2_var(0,6,255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'ictop', field%ictop,                &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
                & lrestart = .FALSE., ldims=shape2d,                       &
                & isteptype=TSTEP_INSTANT )

    ! &       field% totte  (nproma,nlev  ,nblks), &
    cf_desc    = t_cf_var('total_turbulent_energy', 'J kg-1', 'total turbulent energy', &
         &                datatype_flt)
    grib2_desc = grib2_var(0,19,11, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'totte', field%totte,                &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,  &
                & lrestart = .FALSE., initval = 1.e-4_wp, ldims=shape3d,   &
                & vert_interp=create_vert_interp_metadata(                 &
                &             vert_intp_type=vintp_types("P","Z","I"),     &
                &             vert_intp_method=VINTP_METHOD_LIN,           &
                &             l_loglin=.FALSE.,                            &
                &             l_extrapol=.TRUE., l_pd_limit=.FALSE.,       &
                &             lower_limit=0._wp ) )

    ! &       field% thvsig (nproma,       nblks), &
    cf_desc    = t_cf_var('thvsig', 'K', '', datatype_flt)
    grib2_desc = grib2_var(0,19,255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'thvsig', field%thvsig,              &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
                & lrestart = .FALSE., initval = 1.e-2_wp, ldims=shape2d )

    !---------------------------
    ! WMO tropopause
    !---------------------------

    ! &       field% ptp (nproma,       nblks), &
    cf_desc    = t_cf_var('ptp', 'Pa', 'tropopause air pressure', datatype_flt)
    grib2_desc = grib2_var(0,6,1, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'ptp', field%ptp,          &
         &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE,            &
         &        cf_desc, grib2_desc,                           &
         &        ldims=shape2d,                                 &
         &        lrestart = .TRUE., initval = 20000.0_wp,       &
         &        isteptype=TSTEP_INSTANT )

    !---------------------------
    ! Variables for energy diagnostic of echam6 physics
    !---------------------------

    CALL add_var( field_list, prefix//'cpair', field%cpair,                       &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                           &
                & t_cf_var('cpair', 'J/kg/K',                                     &
                &          'specific heat of air at constant pressure',           &
                &          datatype_flt),                                         &
                & grib2_var(0,0,255, ibits, GRID_UNSTRUCTURED, GRID_CELL),        &
                & ldims=shape3d,                                                  &
                & lrestart = .FALSE.,                                             &
                & vert_interp=create_vert_interp_metadata(                        &
                &   vert_intp_type=vintp_types("P","Z","I"),                      &
                &   vert_intp_method=VINTP_METHOD_LIN ) )

    CALL add_var( field_list, prefix//'cvair', field%cvair,                       &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                           &
                & t_cf_var('cvair', 'J/kg/K',                                     &
                &          'specific heat of air at constant colume',             &
                &          datatype_flt),                                         &
                & grib2_var(0,0,255, ibits, GRID_UNSTRUCTURED, GRID_CELL),        &
                & ldims=shape3d,                                                  &
                & lrestart = .FALSE.,                                             &
                & vert_interp=create_vert_interp_metadata(                        &
                &   vert_intp_type=vintp_types("P","Z","I"),                      &
                &   vert_intp_method=VINTP_METHOD_LIN ) )

    CALL add_var( field_list, prefix//'qconv', field%qconv,                       &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                           &
                & t_cf_var('qconv', '(K/s)/(W/m2)',                               &
                &          'conv. factor layer heating to temp. tendency',        &
                &          datatype_flt),                                         &
                & grib2_var(0,0,255, ibits, GRID_UNSTRUCTURED, GRID_CELL),        &
                & ldims=shape3d,                                                  &
                & lrestart = .FALSE.,                                             &
                & vert_interp=create_vert_interp_metadata(                        &
                &   vert_intp_type=vintp_types("P","Z","I"),                      &
                &   vert_intp_method=VINTP_METHOD_LIN ) )

    IF (is_variable_in_output(first_output_name_list, var_name=prefix//'q_phy')) THEN
       CALL add_var( field_list, prefix//'q_phy', field%q_phy,                       &
                   & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                           &
                   & t_cf_var('q_phy', 'W m-2',                                      &
                   &          'layer heating by physics',                            &
                   &          datatype_flt),                                         &
                   & grib2_var(0,0,255, ibits, GRID_UNSTRUCTURED, GRID_CELL),        &
                   & ldims=shape3d,                                                  &
                   & lrestart = .FALSE.,                                             &
                   & vert_interp=create_vert_interp_metadata(                        &
                   &   vert_intp_type=vintp_types("P","Z","I"),                      &
                   &   vert_intp_method=VINTP_METHOD_LIN ) )
    END IF

    IF (is_variable_in_output(first_output_name_list, var_name=prefix//'q_phy_vi')) THEN
       CALL add_var( field_list, prefix//'q_phy_vi', field%q_phy_vi,                 &
                   & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                             &
                   & t_cf_var('q_phy_vi', 'W m-2',                                   &
                   &          'vert. integr. heating by physics',                    &
                   &          datatype_flt),                                         &
                   & grib2_var(0,0,255, ibits, GRID_UNSTRUCTURED, GRID_CELL),        &
                   & ldims=shape2d,                                                  &
                   & lrestart = .FALSE. )
    END IF

    IF ( echam_phy_tc(jg)%dt_rad > dt_zero ) THEN
       !
       IF (is_variable_in_output(first_output_name_list, var_name=prefix//'q_rad')) THEN
          CALL add_var( field_list, prefix//'q_rlw', field%q_rad,                       &
                      & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                           &
                      & t_cf_var('q_rad', 'W m-2',                                      &
                      &          'layer heating by LW+SW radiation',                    &
                      &          datatype_flt),                                         &
                      & grib2_var(0,0,255, ibits, GRID_UNSTRUCTURED, GRID_CELL),        &
                      & ldims=shape3d,                                                  &
                      & lrestart = .FALSE.,                                             &
                      & vert_interp=create_vert_interp_metadata(                        &
                      &   vert_intp_type=vintp_types("P","Z","I"),                      &
                      &   vert_intp_method=VINTP_METHOD_LIN ) )
       END IF
       !
       IF (is_variable_in_output(first_output_name_list, var_name=prefix//'q_rad_vi')) THEN
          CALL add_var( field_list, prefix//'q_rad_vi', field%q_rad_vi,                 &
                      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                             &
                      & t_cf_var('q_rlw_vi', 'W m-2',                                   &
                      &          'vert. integr. heating by LW+SW radiation',            &
                      &          datatype_flt),                                         &
                      & grib2_var(0,0,255, ibits, GRID_UNSTRUCTURED, GRID_CELL),        &
                      & ldims=shape2d,                                                  &
                      & lrestart = .FALSE. )
       END IF
       !
       IF (is_variable_in_output(first_output_name_list, var_name=prefix//'q_rlw')) THEN
          CALL add_var( field_list, prefix//'q_rlw', field%q_rlw,                       &
                      & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                           &
                      & t_cf_var('q_rlw', 'W m-2',                                      &
                      &          'layer heating by LW radiation',                       &
                      &          datatype_flt),                                         &
                      & grib2_var(0,0,255, ibits, GRID_UNSTRUCTURED, GRID_CELL),        &
                      & ldims=shape3d,                                                  &
                      & lrestart = .FALSE.,                                             &
                      & vert_interp=create_vert_interp_metadata(                        &
                      &   vert_intp_type=vintp_types("P","Z","I"),                      &
                      &   vert_intp_method=VINTP_METHOD_LIN ) )
       END IF
       !
       IF (is_variable_in_output(first_output_name_list, var_name=prefix//'q_rlw_vi')) THEN
          CALL add_var( field_list, prefix//'q_rlw_vi', field%q_rlw_vi,                 &
                      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                             &
                      & t_cf_var('q_rlw_vi', 'W m-2',                                   &
                      &          'vert. integr. heating by LW radiation',               &
                      &          datatype_flt),                                         &
                      & grib2_var(0,0,255, ibits, GRID_UNSTRUCTURED, GRID_CELL),        &
                      & ldims=shape2d,                                                  &
                      & lrestart = .FALSE. )
       END IF
       !
       IF (is_variable_in_output(first_output_name_list, var_name=prefix//'q_rsw')) THEN
          CALL add_var( field_list, prefix//'q_rsw', field%q_rsw,                       &
                      & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                           &
                      & t_cf_var('q_rsw', 'W m-2',                                      &
                      &          'layer heating by SW radiation',                       &
                      &          datatype_flt),                                         &
                      & grib2_var(0,0,255, ibits, GRID_UNSTRUCTURED, GRID_CELL),        &
                      & ldims=shape3d,                                                  &
                      & lrestart = .FALSE.,                                             &
                      & vert_interp=create_vert_interp_metadata(                        &
                      &   vert_intp_type=vintp_types("P","Z","I"),                      &
                      &   vert_intp_method=VINTP_METHOD_LIN ) )
       END IF
       !
       IF (is_variable_in_output(first_output_name_list, var_name=prefix//'q_rsw_vi')) THEN
          CALL add_var( field_list, prefix//'q_rsw_vi', field%q_rsw_vi,                 &
                      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                             &
                      & t_cf_var('q_rsw_vi', 'W m-2',                                   &
                      &          'vert. integr. heating by SW radiation',               &
                      &          datatype_flt),                                         &
                      & grib2_var(0,0,255, ibits, GRID_UNSTRUCTURED, GRID_CELL),        &
                      & ldims=shape2d,                                                  &
                      & lrestart = .FALSE. )
       END IF
       !
    END IF

    IF ( echam_phy_tc(jg)%dt_vdf > dt_zero ) THEN
       !
       IF (echam_phy_tc(jg)%dt_vdf > time_config%tc_dt_dyn(jg) .OR.                     &
         & is_variable_in_output(first_output_name_list, var_name=prefix//'q_vdf')) THEN
          CALL add_var( field_list, prefix//'q_vdf', field%q_vdf,                       &
                      & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                           &
                      & t_cf_var('q_vdf', 'W m-2',                                      &
                      &          'layer heating by vertical diffusion',                 &
                      &          datatype_flt),                                         &
                      & grib2_var(0,0,255, ibits, GRID_UNSTRUCTURED, GRID_CELL),        &
                      & ldims=shape3d,                                                  &
                      & lrestart = .FALSE.,                                             &
                      & vert_interp=create_vert_interp_metadata(                        &
                      &   vert_intp_type=vintp_types("P","Z","I"),                      &
                      &   vert_intp_method=VINTP_METHOD_LIN ) )
       END IF
       !
       IF (is_variable_in_output(first_output_name_list, var_name=prefix//'q_vdf_vi')) THEN
          CALL add_var( field_list, prefix//'q_vdf_vi', field%q_vdf_vi,                 &
                      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                             &
                      & t_cf_var('q_vdf_vi', 'W m-2',                                   &
                      &          'vert. integr. heating by vertical diffusion',         &
                      &          datatype_flt),                                         &
                      & grib2_var(0,0,255, ibits, GRID_UNSTRUCTURED, GRID_CELL),        &
                      & ldims=shape2d,                                                  &
                      & lrestart = .FALSE. )
       END IF
       !
    END IF

    IF ( echam_phy_tc(jg)%dt_cnv > dt_zero ) THEN
       !
       IF (echam_phy_tc(jg)%dt_cnv > time_config%tc_dt_dyn(jg) .OR.                     &
         & is_variable_in_output(first_output_name_list, var_name=prefix//'q_cnv')) THEN
          CALL add_var( field_list, prefix//'q_cnv', field%q_cnv,                       &
                      & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                           &
                      & t_cf_var('q_cnv', 'W m-2',                                      &
                      &          'layer heating by vertical diffusion',                 &
                      &          datatype_flt),                                         &
                      & grib2_var(0,0,255, ibits, GRID_UNSTRUCTURED, GRID_CELL),        &
                      & ldims=shape3d,                                                  &
                      & lrestart = .FALSE.,                                             &
                      & vert_interp=create_vert_interp_metadata(                        &
                      &   vert_intp_type=vintp_types("P","Z","I"),                      &
                      &   vert_intp_method=VINTP_METHOD_LIN ) )
       END IF
       !
       IF (is_variable_in_output(first_output_name_list, var_name=prefix//'q_cnv_vi')) THEN
          CALL add_var( field_list, prefix//'q_cnv_vi', field%q_cnv_vi,                 &
                      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                             &
                      & t_cf_var('q_cnv_vi', 'W m-2',                                   &
                      &          'vert. integr. heating by vertical diffusion',         &
                      &          datatype_flt),                                         &
                      & grib2_var(0,0,255, ibits, GRID_UNSTRUCTURED, GRID_CELL),        &
                      & ldims=shape2d,                                                  &
                      & lrestart = .FALSE. )
       END IF
       !
    END IF

    IF ( echam_phy_tc(jg)%dt_cld > dt_zero ) THEN
       !
       IF (echam_phy_tc(jg)%dt_cld > time_config%tc_dt_dyn(jg) .OR.                     &
         & is_variable_in_output(first_output_name_list, var_name=prefix//'q_cld')) THEN
          CALL add_var( field_list, prefix//'q_cld', field%q_cld,                       &
                      & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                           &
                      & t_cf_var('q_cld', 'W m-2',                                      &
                      &          'layer heating by vertical diffusion',                 &
                      &          datatype_flt),                                         &
                      & grib2_var(0,0,255, ibits, GRID_UNSTRUCTURED, GRID_CELL),        &
                      & ldims=shape3d,                                                  &
                      & lrestart = .FALSE.,                                             &
                      & vert_interp=create_vert_interp_metadata(                        &
                      &   vert_intp_type=vintp_types("P","Z","I"),                      &
                      &   vert_intp_method=VINTP_METHOD_LIN ) )
       END IF
       !
       IF (is_variable_in_output(first_output_name_list, var_name=prefix//'q_cld_vi')) THEN
          CALL add_var( field_list, prefix//'q_cld_vi', field%q_cld_vi,                 &
                      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                             &
                      & t_cf_var('q_cld_vi', 'W m-2',                                   &
                      &          'vert. integr. heating by vertical diffusion',         &
                      &          datatype_flt),                                         &
                      & grib2_var(0,0,255, ibits, GRID_UNSTRUCTURED, GRID_CELL),        &
                      & ldims=shape2d,                                                  &
                      & lrestart = .FALSE. )
       END IF
       !
    END IF

    IF ( echam_phy_tc(jg)%dt_gwd > dt_zero ) THEN
       !
       IF (echam_phy_tc(jg)%dt_gwd > time_config%tc_dt_dyn(jg) .OR.                     &
         & is_variable_in_output(first_output_name_list, var_name=prefix//'q_gwd')) THEN
          CALL add_var( field_list, prefix//'q_gwd', field%q_gwd,                       &
                      & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                           &
                      & t_cf_var('q_gwd', 'W m-2',                                      &
                      &          'layer heating by atm. gravity wave drag',             &
                      &          datatype_flt),                                         &
                      & grib2_var(0,0,255, ibits, GRID_UNSTRUCTURED, GRID_CELL),        &
                      & ldims=shape3d,                                                  &
                      & lrestart = .FALSE.,                                             &
                      & vert_interp=create_vert_interp_metadata(                        &
                      &   vert_intp_type=vintp_types("P","Z","I"),                      &
                      &   vert_intp_method=VINTP_METHOD_LIN ) )
       END IF
       !
       IF (is_variable_in_output(first_output_name_list, var_name=prefix//'q_gwd_vi')) THEN
          CALL add_var( field_list, prefix//'q_gwd_vi', field%q_gwd_vi,                 &
                      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                             &
                      & t_cf_var('q_gwd_vi', 'W m-2',                                   &
                      &          'vert. integr. heating by atm. gravity wave drag',     &
                      &          datatype_flt),                                         &
                      & grib2_var(0,0,255, ibits, GRID_UNSTRUCTURED, GRID_CELL),        &
                      & ldims=shape2d,                                                  &
                      & lrestart = .FALSE. )
       END IF
       !
    END IF

    IF ( echam_phy_tc(jg)%dt_sso > dt_zero ) THEN
       !
       IF (echam_phy_tc(jg)%dt_sso > time_config%tc_dt_dyn(jg) .OR.                     &
         & is_variable_in_output(first_output_name_list, var_name=prefix//'q_sso')) THEN
          CALL add_var( field_list, prefix//'q_sso', field%q_sso,                       &
                      & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                           &
                      & t_cf_var('q_sso', 'W m-2',                                      &
                      &          'layer heating by atm. gravity wave drag',             &
                      &          datatype_flt),                                         &
                      & grib2_var(0,0,255, ibits, GRID_UNSTRUCTURED, GRID_CELL),        &
                      & ldims=shape3d,                                                  &
                      & lrestart = .FALSE.,                                             &
                      & vert_interp=create_vert_interp_metadata(                        &
                      &   vert_intp_type=vintp_types("P","Z","I"),                      &
                      &   vert_intp_method=VINTP_METHOD_LIN ) )
       END IF
       !
       IF (is_variable_in_output(first_output_name_list, var_name=prefix//'q_sso_vi')) THEN
          CALL add_var( field_list, prefix//'q_sso_vi', field%q_sso_vi,                 &
                      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                             &
                      & t_cf_var('q_sso_vi', 'W m-2',                                   &
                      &          'vert. integr. heating by atm. gravity wave drag',     &
                      &          datatype_flt),                                         &
                      & grib2_var(0,0,255, ibits, GRID_UNSTRUCTURED, GRID_CELL),        &
                      & ldims=shape2d,                                                  &
                      & lrestart = .FALSE. )
       END IF
       !
    END IF

       cf_desc    = t_cf_var('sh_vdiff','J m-2 s-1', '', datatype_flt)
       grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var( field_list, prefix//'sh_vdiff', field%sh_vdiff,          &
                   & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
                   & lrestart = .FALSE., ldims=shape2d )
       cf_desc    = t_cf_var('con_dtrl','?', '', datatype_flt)
       grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var( field_list, prefix//'con_dtrl', field%con_dtrl,          &
                   & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
                   & lrestart = .FALSE., ldims=shape2d )
       cf_desc    = t_cf_var('con_dtri','?', '', datatype_flt)
       grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var( field_list, prefix//'con_dtri', field%con_dtri,          &
                   & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
                   & lrestart = .FALSE., ldims=shape2d )
       cf_desc    = t_cf_var('con_iteqv','?', '', datatype_flt)
       grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var( field_list, prefix//'con_iteqv', field%con_iteqv,        &
                   & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
                   & lrestart = .FALSE., ldims=shape2d )
       cf_desc    = t_cf_var('qv_vdiff','kg/m^2/s', '', datatype_flt)
       grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var( field_list, prefix//'qv_vdiff', field%qv_vdiff,          &
                   & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
                   & lrestart = .FALSE., ldims=shape2d )

    IF ( echam_phy_tc(jg)%dt_sso > dt_zero ) THEN
       !
       !---------------------------
       ! Sub grid scale orographic effects (sso)
       !---------------------------
       IF (is_variable_in_output(first_output_name_list, var_name=prefix//'tauu_sso')) THEN
          CALL add_var( field_list, prefix//'tauu_sso', field%u_stress_sso,             &
                      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                             &
                      & t_cf_var('u_stress_sso', 'N m-2',                               &
                      &          'zonal stress from subgrid scale orographic drag',     &
                      &          datatype_flt),                                         &
                      & grib2_var(0,2,17, ibits, GRID_UNSTRUCTURED, GRID_CELL),         &
                      & ldims=shape2d,                                                  &
                      & lrestart = .FALSE.,                                             &
                      & isteptype=TSTEP_INSTANT                                         )
       END IF
       !
       IF (is_variable_in_output(first_output_name_list, var_name=prefix//'tauv_sso')) THEN
          CALL add_var( field_list, prefix//'tauv_sso', field%v_stress_sso,             &
                      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                             &
                      & t_cf_var('v_stress_sso', 'N m-2',                               &
                      &          'meridional stress from subgrid scale orographic drag',&
                      &          datatype_flt),                                         &
                      & grib2_var(0,2,18, ibits, GRID_UNSTRUCTURED, GRID_CELL),         &
                      & ldims=shape2d,                                                  &
                      & lrestart = .FALSE.,                                             &
                      & isteptype=TSTEP_INSTANT                                         )
       END IF
       !
       IF (is_variable_in_output(first_output_name_list, var_name=prefix//'diss_sso')) THEN
          CALL add_var( field_list, prefix//'diss_sso', field%dissipation_sso,          &
                      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                             &
                      & t_cf_var('dissipation_sso', '',                                 &
                      &          'dissipation of orographic waves',                     &
                      &          datatype_flt),                                         &
                      & grib2_var(0,0,255, ibits, GRID_UNSTRUCTURED, GRID_CELL),        &
                      & ldims=shape2d,                                                  &
                      & lrestart = .FALSE.,                                             &
                      & isteptype=TSTEP_INSTANT                                         )
       END IF
    !
    END IF

    !--------------------
    ! Turbulence
    !--------------------

     ! shape2d  = (/kproma,            kblks/)
     ! shape3d  = (/kproma, klev,      kblks/)
     !shapesfc = (/kproma, ksfc_type, kblks/)
     ! shapesfc = (/kproma, kblks, ksfc_type/)

      IF (is_variable_in_output(first_output_name_list, var_name=prefix//'ri_atm')) THEN
         cf_desc    = t_cf_var('richardson_number', ' ', 'moist Richardson number', datatype_flt)
         grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
         CALL add_var( field_list, prefix//'ri_atm', field%ri_atm,              &
                   & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc, &
                   & lrestart = .FALSE., ldims=shape3d )
      END IF

      IF (is_variable_in_output(first_output_name_list, var_name=prefix//'mixlen')) THEN
         cf_desc    = t_cf_var('mixing_length', 'm', 'mixing_length', datatype_flt)
         grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
         CALL add_var( field_list, prefix//'mixlen', field%mixlen,           &
                   & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc, &
                   & lrestart = .FALSE., initval = -999._wp, ldims=shape3d )
      END IF

      ! &       field% tottem0 (nproma,nlev,nblks), &
      cf_desc    = t_cf_var('totte', 'm2 s-2', 'TTE at step t', datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( field_list, prefix//'tottem0', field%tottem0,              &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc, &
                  & lrestart = .FALSE., initval = 1.e-4_wp, ldims=shape3d )

      ! &       field% tottem1  (nproma,nlev,nblks), &
      cf_desc    = t_cf_var('totte', 'm2 s-2', 'TTE at step t-dt', datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( field_list, prefix//'tottem1', field%tottem1,              &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc, &
                  & lrestart = .TRUE., initval = 1.e-4_wp, ldims=shape3d )

      
      IF (is_variable_in_output(first_output_name_list, var_name=prefix//'cfm')) THEN
         cf_desc    = t_cf_var('turb_exchng_coeff_momentum', '', '', datatype_flt)
         grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
         CALL add_var( field_list, prefix//'cfm', field%cfm,                      &
                     & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc, &
                     & lrestart = .FALSE., ldims=shape3d )
      END IF
      !
      contvar_is_in_output = .FALSE.
      DO jsfc = 1,ksfc_type
         varname=prefix//'cfm_'//csfc(jsfc)
         IF (is_variable_in_output(first_output_name_list, var_name=TRIM(varname))) THEN
            contvar_is_in_output = .TRUE.
         END IF
      END DO
      !
      IF (contvar_is_in_output) THEN
         CALL add_var( field_list, prefix//'cfm_tile', field%cfm_tile,              &
                     & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                          &
                     & t_cf_var('turb_exchng_coeff_momentum', '', '', datatype_flt),&
                     & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED,GRID_CELL),&
                     & ldims=shapesfc, lmiss=.TRUE., missval=cdimissval,            &
                     & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.         )
         ALLOCATE(field%cfm_tile_ptr(ksfc_type))
      END IF
      !
      DO jsfc = 1,ksfc_type
         varname=prefix//'cfm_'//csfc(jsfc)
         IF (is_variable_in_output(first_output_name_list, var_name=TRIM(varname))) THEN
            CALL add_ref( field_list, prefix//'cfm_tile',                              &
                        & TRIM(varname), field%cfm_tile_ptr(jsfc)%p,                   &
                        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                          &
                        & t_cf_var('turb_exchng_coeff_momentum_'//csfc(jsfc), '', '', datatype_flt),&
                        & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED,GRID_CELL),&
                        & lrestart=.FALSE., ldims=shape2d,                             &
                        & lmiss=.TRUE., missval=cdimissval )
         END IF
      END DO


      IF (is_variable_in_output(first_output_name_list, var_name=prefix//'cfh')) THEN
         cf_desc    = t_cf_var('turb_exchng_coeff_heat', '', '', datatype_flt)
         grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
         CALL add_var( field_list, prefix//'cfh', field%cfh,                      &
                     & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc, &
                     & lrestart = .FALSE., ldims=shape3d )
      END IF
      !
      contvar_is_in_output = .FALSE.
      DO jsfc = 1,ksfc_type
         varname=prefix//'cfh_'//csfc(jsfc)
         IF (is_variable_in_output(first_output_name_list, var_name=TRIM(varname))) THEN
            contvar_is_in_output = .TRUE.
         END IF
      END DO
      !
      IF (contvar_is_in_output) THEN
         CALL add_var( field_list, prefix//'cfh_tile', field%cfh_tile,              &
                     & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                          &
                     & t_cf_var('turb_exchng_coeff_heat', '', '', datatype_flt),    &
                     & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED,GRID_CELL),&
                     & ldims=shapesfc, lmiss=.TRUE., missval=cdimissval,            &
                     & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.         )
         ALLOCATE(field%cfh_tile_ptr(ksfc_type))
      END IF
      !
      DO jsfc = 1,ksfc_type
         varname=prefix//'cfh_'//csfc(jsfc)
         IF (is_variable_in_output(first_output_name_list, var_name=TRIM(varname))) THEN
            CALL add_ref( field_list, prefix//'cfh_tile',                              &
                        & TRIM(varname), field%cfh_tile_ptr(jsfc)%p,                   &
                        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                          &
                        & t_cf_var('turb_exchng_coeff_heat_'//csfc(jsfc), '', '',      &
                        &          datatype_flt),                                      &
                        & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED,GRID_CELL),&
                        & lrestart = .FALSE., ldims=shape2d,                           &
                        & lmiss=.TRUE., missval=cdimissval )
         END IF
      END DO


      IF (is_variable_in_output(first_output_name_list, var_name=prefix//'cfv')) THEN
         cf_desc    = t_cf_var('turb_exchng_coeff_water_var', '', '', datatype_flt)
         grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
         CALL add_var( field_list, prefix//'cfv', field%cfv,                      &
                     & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc, &
                     & lrestart = .FALSE., ldims=shape3d )
      END IF

      IF (is_variable_in_output(first_output_name_list, var_name=prefix//'cfv')) THEN
         cf_desc    = t_cf_var('turb_exchng_coeff_totte', '', '', datatype_flt)
         grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
         CALL add_var( field_list, prefix//'cftotte', field%cftotte,              &
                     & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc, &
                     & lrestart = .FALSE., ldims=shape3d )
      END IF

      IF (is_variable_in_output(first_output_name_list, var_name=prefix//'cfthv')) THEN
         cf_desc    = t_cf_var('turb_exchng_coeff_thv', '', '', datatype_flt)
         grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
         CALL add_var( field_list, prefix//'cfthv', field%cfthv,                  &
                     & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc, &
                     & lrestart = .FALSE., ldims=shape3d )
      END IF

      cf_desc    = t_cf_var('Coriolis_param', 's-1', 'Coriolis parameter', datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( field_list, prefix//'coriol', field%coriol,              &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
                  & lrestart = .FALSE., ldims=shape2d )

      IF (is_variable_in_output(first_output_name_list, var_name=prefix//'hdtcbl')) THEN
         cf_desc    = t_cf_var('height_pbl_top', 'm', 'height of PBL top', datatype_flt)
         grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
         CALL add_var( field_list, prefix//'hdtcbl', field%hdtcbl,              &
                     & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
                     & lrestart = .FALSE., ldims=shape2d )
      END IF

      !-----------------------------------
      ! &       field% z0m(nproma,nblks), &
      cf_desc    = t_cf_var('z0m', '', 'aerodynamic roughness length, grid box mean', &
           &                datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( field_list, prefix//'z0m', field%z0m,                    &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
                  & lrestart = .FALSE., ldims=shape2d )

      ! &       field% z0m_tile(nproma,nblks,nsfc_type), &
      CALL add_var( field_list, prefix//'z0m_tile', field%z0m_tile,              &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                          &
                  & t_cf_var('z0m_tile', '', 'aerodynamic roughness length', datatype_flt),&
                  & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED,GRID_CELL),&
                  & ldims=shapesfc, lmiss=.TRUE., missval=cdimissval,            &
                  & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.         )

      ALLOCATE(field%z0m_tile_ptr(ksfc_type))
      DO jsfc = 1,ksfc_type
        CALL add_ref( field_list, prefix//'z0m_tile',                              &
                    & prefix//'z0m_'//csfc(jsfc), field%z0m_tile_ptr(jsfc)%p,      &
                    & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                          &
                    & t_cf_var('z0m_'//csfc(jsfc), '','', datatype_flt),           &
                    & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED,GRID_CELL),&
                    & lrestart = .TRUE., ldims=shape2d,                            &
                    & lmiss=.TRUE., missval=cdimissval )
      END DO

   


      ! &        field% z0h_lnd(nproma, nblks), &
      cf_desc    = t_cf_var('z0h_lnd', '', 'roughness length heat, land', &
        &                datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( field_list, prefix//'z0h_lnd', field%z0h_lnd,            &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
                  & lrestart = .TRUE., ldims=shape2d,                        &
                  & lmiss=.TRUE., missval=cdimissval )

      !-----------------------------------

      ! &       field% ustar  (nproma,nblks),                &
      cf_desc    = t_cf_var('friction_velocity', 'm s-1', 'friction velocity', datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( field_list, prefix//'ustar', field%ustar,                &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
                  & lrestart = .TRUE., initval = 1._wp, ldims=shape2d )

      IF (is_variable_in_output(first_output_name_list, var_name=prefix//'wstar')) THEN
         cf_desc    = t_cf_var('conv_velocity_scale', 'm s-1', 'convective velocity scale', datatype_flt)
         grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
         CALL add_var( field_list, prefix//'wstar', field%wstar,                &
                     & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
                     & lrestart = .FALSE., ldims=shape2d )
      END IF

      ! &       field% wstar_tile(nproma,nblks,nsfc_type), &
      CALL add_var( field_list, prefix//'wstar_tile', field%wstar_tile,          &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                          &
                  & t_cf_var('wstar_tile', '', 'convective velocity scale', datatype_flt),&
                  & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED,GRID_CELL),&
                  & ldims=shapesfc, lcontainer=.TRUE.,                           &
                  & lrestart=.FALSE., loutput=.FALSE.                            )

      ALLOCATE(field%wstar_tile_ptr(ksfc_type))
      DO jsfc = 1,ksfc_type
        CALL add_ref( field_list, prefix//'wstar_tile',                            &
                    & prefix//'wstar_'//csfc(jsfc), field%wstar_tile_ptr(jsfc)%p,  &
                    & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                          &
                    & t_cf_var('wstar_'//csfc(jsfc), '','', datatype_flt),         &
                    & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED,GRID_CELL),&
                    & lrestart=.TRUE., ldims=shape2d                               )
      END DO


      IF (is_variable_in_output(first_output_name_list, var_name=prefix//'kedisp')) THEN
         cf_desc    = t_cf_var('vert_int_dissip_kin_energy', 'W/m2',            &
                     &         'vert. integr. dissip. kin. energy', datatype_flt)
         grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
         CALL add_var( field_list, prefix//'kedisp', field%kedisp,              &
                     & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
                     & lrestart=.FALSE., ldims=shape2d )
      END IF

      ! &       field% ocu    (nproma,nblks),                &
      cf_desc    = t_cf_var('ocean_sfc_u', 'm/s', 'u-component of ocean current/ice', datatype_flt)
      grib2_desc = grib2_var(10,1,2, iextbits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( field_list, prefix//'ocu', field%ocu, &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE,   &
        &           cf_desc, grib2_desc, ldims=shape2d,   &
        &           lrestart=.TRUE. )

      ! &       field% ocv    (nproma,nblks),                &
      cf_desc    = t_cf_var('ocean_sfc_v', 'm/s', 'v-component of ocean current/ice', datatype_flt)
      grib2_desc = grib2_var(10,1,3, iextbits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( field_list, prefix//'ocv', field%ocv, &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE,   &
        &           cf_desc, grib2_desc, ldims=shape2d,   &
        &           lrestart=.TRUE. )

    !-----------------------
    ! Surface
    !-----------------------

    cf_desc    = t_cf_var('surface_altitude', 'm',   &
                &         'surface altitude', datatype_flt)
    grib2_desc = grib2_var(2,0,7, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'orog', field%orog,                  &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
                & lrestart = .FALSE., ldims=shape2d,                       &
                & isteptype=TSTEP_CONSTANT )

    cf_desc    = t_cf_var('land_area_fraction', 'm2/m2',   &
                &         'cell area fraction occupied by land including lakes', datatype_flt)
    grib2_desc = grib2_var(2,0,0, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'sftlf', field%sftlf,                &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
                & lrestart = .FALSE., ldims=shape2d,                       &
                & isteptype=TSTEP_CONSTANT )

    cf_desc    = t_cf_var('land_ice_area_fraction', 'm2/m2',   &
                &         'cell area fraction occupied by land ice', datatype_flt)
    grib2_desc = grib2_var(255,255,255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'sftgif', field%sftgif,              &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
                & lrestart = .FALSE., ldims=shape2d,                       &
                & isteptype=TSTEP_CONSTANT )

    cf_desc    = t_cf_var('ocean_area_fraction', 'm2/m2',   &
                &         'cell area fraction occupied by ocean', datatype_flt)
    grib2_desc = grib2_var(2,0,255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'sftof', field%sftof,                &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
                & lrestart = .FALSE., ldims=shape2d,                       &
                & isteptype=TSTEP_CONSTANT )


    ! &       field% lsmask (nproma, nblks),                 &
    cf_desc    = t_cf_var('land_cover', '', 'land cover', datatype_flt)
    grib2_desc = grib2_var(1, 2, 8, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'land', field%lsmask,              &
              & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
              & lrestart=.FALSE., ldims=shape2d )

    ! &       field% glac   (nproma, nblks),                 &
    cf_desc    = t_cf_var('glacier_cover', '', 'fraction of land covered by glaciers', &
         &                datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'glac', field%glac,                &
              & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
              & lrestart=.FALSE., ldims=shape2d )

    ! &       field% seaice (nproma, nblks),                 &
    cf_desc    = t_cf_var('sea_ice_cover', '', 'fraction of ocean covered by sea ice', &
         &                datatype_flt)
    grib2_desc = grib2_var(10,2,0, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'sic', field%seaice,    &
      &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE,         &
      &           cf_desc, grib2_desc, ldims=shape2d,         &
      &           lrestart=.TRUE. )

    ! &       field% alake (nproma, nblks),                 &
    cf_desc    = t_cf_var('alake', '', 'fraction of lakes', &
         &                datatype_flt)
    grib2_desc = grib2_var(1,2,2, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'alake', field%alake,              &
              & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
              & lrestart=.FALSE., ldims=shape2d )

    ! &       field% lake_ice_frc (nproma, nblks),                 &
    cf_desc    = t_cf_var('lake_ice_frc', '', 'fraction of ice on lakes', & 
         &                datatype_flt)
    grib2_desc = grib2_var(1,2,7, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'lake_ice_frc', field%lake_ice_frc,  &
              & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,   &
              & initval=0._wp, lrestart=.TRUE., ldims=shape2d )

    !-----------------------------------
    ! &       field% ts(nproma,nblks), &
    cf_desc    = t_cf_var('surface_temperature', 'K', 'surface temperature', datatype_flt)
    grib2_desc = grib2_var(0,0,0, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'ts', field%ts,                    &
              & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
              & lrestart=.TRUE., ldims=shape2d )

    ! &       field% ts_tile(nproma,nblks,nsfc_type), &
    CALL add_var( field_list, prefix//'ts_tile', field%ts_tile,                &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                          &
                & t_cf_var('_tile', 'K', 'surface temperature on tiles', datatype_flt), &
                & grib2_var(0,0,0, ibits, GRID_UNSTRUCTURED, GRID_CELL),       &
                & ldims=shapesfc, lmiss=.TRUE., missval=cdimissval,            &
                & lcontainer=.TRUE., lrestart=.FALSE.         )

    ALLOCATE(field%ts_tile_ptr(ksfc_type))
    DO jsfc = 1,ksfc_type
      CALL add_ref( field_list, prefix//'ts_tile',                                   &
                  & prefix//'ts_'//csfc(jsfc), field%ts_tile_ptr(jsfc)%p,            &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                              &
                  & t_cf_var('ts_'//csfc(jsfc), 'K',                                 &
                  &          'surface temperature on '//csfc(jsfc), datatype_flt),   &
                  & grib2_var(0,0,0, ibits, GRID_UNSTRUCTURED, GRID_CELL),           &
                  & lrestart=.TRUE., ldims=shape2d,                                  &
                  & lmiss=.TRUE., missval=cdimissval )
    END DO
    !-----------------------------------

    contvar_is_in_output = .FALSE.
    DO jsfc = 1,ksfc_type
       varname=prefix//'qs_sfc_'//csfc(jsfc)
       IF (is_variable_in_output(first_output_name_list, var_name=TRIM(varname))) THEN
          contvar_is_in_output = .TRUE.
       END IF
    END DO
    !
    IF (contvar_is_in_output) THEN
       CALL add_var( field_list, prefix//'qs_sfc_tile', field%qs_sfc_tile,        &
                   & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                          &
                   & t_cf_var('qs_sfc_tile', '', '', datatype_flt),               &
                   & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED,GRID_CELL),&
                   & ldims=shapesfc, lmiss=.TRUE., missval=cdimissval,            &
                   & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.         )
       ALLOCATE(field%qs_sfc_tile_ptr(ksfc_type))
    END IF
    !
    DO jsfc = 1,ksfc_type
       varname=prefix//'qs_sfc_'//csfc(jsfc)
       IF (is_variable_in_output(first_output_name_list, var_name=TRIM(varname))) THEN
          CALL add_ref( field_list, prefix//'qs_sfc_tile',                           &
                      & TRIM(varname), field%qs_sfc_tile_ptr(jsfc)%p,                &
                      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                          &
                      & t_cf_var('qs_sfc_'//csfc(jsfc), '', '', datatype_flt),       &
                      & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED,GRID_CELL),&
                      & lrestart=.FALSE., ldims=shape2d,                             &
                      & lmiss=.TRUE., missval=cdimissval )
       END IF
    END DO

    !-----------------------------------
    ! &       field% albedo (nproma,nblks),          &
    cf_desc    = t_cf_var('albedo', '', 'surface albedo', datatype_flt)
    grib2_desc = grib2_var(0,19,1, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'albedo', field%albedo,              &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
                & lrestart=.FALSE., ldims=shape2d )

    ! &       field% albvisdir (nproma,nblks),          &
    cf_desc    = t_cf_var('albvisdir', '', 'albedo VIS direct', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'albvisdir', field%albvisdir,        &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
                & lrestart=.TRUE., ldims=shape2d  )

    ! &       field% albvisdif (nproma,nblks),          &
    cf_desc    = t_cf_var('albvisdif', '', 'albedo VIS diffuse', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'albvisdif', field%albvisdif,        &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
                & lrestart=.TRUE., ldims=shape2d  )

    ! &       field% albnirdir (nproma,nblks),          &
    cf_desc    = t_cf_var('albnirdir', '', 'albedo NIR direct', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'albnirdir', field%albnirdir,        &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
                & lrestart=.TRUE., ldims=shape2d  )

    ! &       field% albnirdif (nproma,nblks),          &
    cf_desc    = t_cf_var('albnirdif', '', 'albedo NIR diffuse', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'albnirdif', field%albnirdif,        &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
                & lrestart=.TRUE., ldims=shape2d  )

    ! &       field% albvisdir_tile (nproma,nblks,nsfc_type),          &
    cf_desc    = t_cf_var('albvisdir_tile', '', 'albedo VIS direct', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'albvisdir_tile', field%albvisdir_tile,             &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,                &
                & ldims=shapesfc, lmiss=.TRUE., missval=cdimissval,                       &
                & lcontainer=.TRUE., lrestart=.FALSE.         )

    ! &       field% albvisdif_tile (nproma,nblks,nsfc_type),          &
    cf_desc    = t_cf_var('albvisdif_tile', '', 'albedo VIS diffuse', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'albvisdif_tile', field%albvisdif_tile,             &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,                &
                & ldims=shapesfc, lmiss=.TRUE., missval=cdimissval,                       &
                & lcontainer=.TRUE., lrestart=.FALSE.         )

    ! &       field% albnirdir_tile (nproma,nblks,nsfc_type),          &
    cf_desc    = t_cf_var('albnirdir_tile', '', 'albedo NIR direct', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'albnirdir_tile', field%albnirdir_tile,             &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,                &
                & ldims=shapesfc, lmiss=.TRUE., missval=cdimissval,                       &
                & lcontainer=.TRUE., lrestart=.FALSE.         )

    ! &       field% albnirdif_tile (nproma,nblks,nsfc_type),          &
    cf_desc    = t_cf_var('albnirdif_tile', '', 'albedo NIR diffuse', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'albnirdif_tile', field%albnirdif_tile,             &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,                &
                & ldims=shapesfc, lmiss=.TRUE., missval=cdimissval,                       &
                & lcontainer=.TRUE., lrestart=.FALSE.         )

    ! &       field% albedo_tile (nproma,nblks,nsfc_type),          &
    cf_desc    = t_cf_var('albedo_tile', '', 'albedo', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'albedo_tile', field%albedo_tile,                   &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,                &
                & ldims=shapesfc, lmiss=.TRUE., missval=cdimissval,                       &
                & lcontainer=.TRUE., lrestart=.FALSE.         )

    ALLOCATE(field%albvisdir_tile_ptr(ksfc_type), field%albvisdif_tile_ptr(ksfc_type), &
             field%albnirdir_tile_ptr(ksfc_type), field%albnirdif_tile_ptr(ksfc_type), &
             field%albedo_tile_ptr(ksfc_type)                                          )

    DO jsfc = 1,ksfc_type
      CALL add_ref( field_list, prefix//'albvisdir_tile',                              &
                  & prefix//'albvisdir_'//csfc(jsfc), field%albvisdir_tile_ptr(jsfc)%p,&
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                &
                  & t_cf_var('albvisdir_'//csfc(jsfc), '', '', datatype_flt),          &
                  & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),     &
                  & lrestart=.FALSE., ldims=shape2d,                                   &
                  & lmiss=.TRUE., missval=cdimissval )
      CALL add_ref( field_list, prefix//'albvisdif_tile',                              &
                  & prefix//'albvisdif_'//csfc(jsfc), field%albvisdif_tile_ptr(jsfc)%p,&
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                &
                  & t_cf_var('albvisdif_'//csfc(jsfc), '', '', datatype_flt),          &
                  & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),     &
                  & lrestart=.FALSE., ldims=shape2d,                                   &
                  & lmiss=.TRUE., missval=cdimissval )
      CALL add_ref( field_list, prefix//'albnirdir_tile',                              &
                  & prefix//'albnirdir_'//csfc(jsfc), field%albnirdir_tile_ptr(jsfc)%p,&
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                &
                  & t_cf_var('albnirdir_'//csfc(jsfc), '', '', datatype_flt),          &
                  & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),     &
                  & lrestart=.FALSE., ldims=shape2d,                                   &
                  & lmiss=.TRUE., missval=cdimissval )
      CALL add_ref( field_list, prefix//'albnirdif_tile',                              &
                  & prefix//'albnirdif_'//csfc(jsfc), field%albnirdif_tile_ptr(jsfc)%p,&
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                &
                  & t_cf_var('albnirdif_'//csfc(jsfc), '', '', datatype_flt),          &
                  & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),     &
                  & lrestart=.FALSE., ldims=shape2d,                                   &
                  & lmiss=.FALSE., missval=cdimissval )
      CALL add_ref( field_list, prefix//'albedo_tile',                                 &
                  & prefix//'albedo_'//csfc(jsfc), field%albedo_tile_ptr(jsfc)%p,      &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                &
                  & t_cf_var('albedo_'//csfc(jsfc), '', '', datatype_flt),             &
                  & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),     &
                  & lrestart=.FALSE., ldims=shape2d,                                   &
                  & lmiss=.TRUE., missval=cdimissval )
    END DO

    !---------------------------
    ! Surface fluxes
    !---------------------------
    ! gridbox mean
    CALL add_var( field_list, prefix//'co2flux', field%co2flux,           &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                     &
                & t_cf_var('co2flux', 'kg m-2 s-1', 'co2 flux',           &
                & datatype_flt),                                          &
                & grib2_var(255,255,255,iextbits, GRID_UNSTRUCTURED, GRID_CELL),&
                & ldims=shape2d,                                          &
                & lrestart = .FALSE.,                                     &
                & isteptype=TSTEP_INSTANT                                 )


    CALL add_var( field_list, prefix//'evspsbl', field%evap,              &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                     &
                & t_cf_var('evap', 'kg m-2 s-1', 'evaporation',           &
                & datatype_flt),                                          &
                & grib2_var(0,1,6,iextbits, GRID_UNSTRUCTURED, GRID_CELL),&
                & ldims=shape2d,                                          &
                & lrestart = .FALSE.,                                     &
                & isteptype=TSTEP_INSTANT                                 )

    CALL add_var( field_list, prefix//'hfls', field%lhflx,                &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                     &
                & t_cf_var('lhflx', 'W m-2 ', 'latent heat flux',         &
                & datatype_flt),                                          &
                & grib2_var(0,0,10, ibits, GRID_UNSTRUCTURED, GRID_CELL), &
                & ldims=shape2d,                                          &
                & lrestart = .FALSE.,                                     &
                & isteptype=TSTEP_INSTANT                                 )

    CALL add_var( field_list, prefix//'hfss', field%shflx,                &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                     &
                & t_cf_var('shflx', 'W m-2 ', 'sensible heat flux',       &
                & datatype_flt),                                          &
                & grib2_var(0,0,11, ibits, GRID_UNSTRUCTURED, GRID_CELL), &
                & ldims=shape2d,                                          &
                & lrestart = .FALSE.,                                     &
                & isteptype=TSTEP_INSTANT                                 )

    !---------------------------------
    ! values on tiles

    CALL add_var( field_list, prefix//'rsns_tile',field%swflxsfc_tile,    &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                     &
                & t_cf_var('rsns_tile', 'W m-2',                          &
                &          'shortwave net flux at surface on tiles',      &
                &          datatype_flt),                                 &
                & grib2_var(0,4,9, ibits, GRID_UNSTRUCTURED, GRID_CELL),  &
                & ldims=shapesfc, lmiss=.TRUE., missval=cdimissval,       &
                & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.    )

    CALL add_var( field_list, prefix//'rlns_tile',field%lwflxsfc_tile,    &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                     &
                & t_cf_var('rlns_tile', 'W m-2',                          &
                &          'longwave net flux at surface on tiles',       &
                &          datatype_flt),                                 &
                & grib2_var(0,5,5, ibits, GRID_UNSTRUCTURED, GRID_CELL),  &
                & ldims=shapesfc, lmiss=.TRUE., missval=cdimissval,       &
                & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.    )

    CALL add_var( field_list, prefix//'evspsbl_tile', field%evap_tile,    &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                     &
                & t_cf_var('evspsbl_tile', 'kg m-2 s-1',                  &
                &          'evaporation on tiles', datatype_flt),         &
                & grib2_var(0,1,6, ibits, GRID_UNSTRUCTURED, GRID_CELL),  &
                & ldims=shapesfc,                                         &
                & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.    )

    CALL add_var( field_list, prefix//'hfls_tile', field%lhflx_tile,      &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                     &
                & t_cf_var('hfls_tile', 'W m-2',                          &
                &          'latent heat flux on tiles', datatype_flt),    &
                & grib2_var(0,0,10, ibits, GRID_UNSTRUCTURED, GRID_CELL), &
                & ldims=shapesfc, lmiss=.TRUE., missval=cdimissval,       &
                & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.    )

    CALL add_var( field_list, prefix//'hfss_tile', field%shflx_tile,      &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                     &
                & t_cf_var('hfss_tile', 'W m-2',                          &
                &          'sensible heat flux on tiles', datatype_flt),  &
                & grib2_var(0,0,11, ibits, GRID_UNSTRUCTURED, GRID_CELL), &
                & ldims=shapesfc, lmiss=.TRUE., missval=cdimissval,       &
                & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.    )

    CALL add_var( field_list, prefix//'frac_tile', field%frac_tile,       &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                     &
                & t_cf_var('frac_tile', '%',                              &
                &          'surface fraction of tiles', datatype_flt),    &
                & grib2_var(255,255,255, ibits, GRID_UNSTRUCTURED, GRID_CELL), &
                & ldims=shapesfc, lmiss=.TRUE., missval=cdimissval,       &
                & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.    )

    ALLOCATE(field%swflxsfc_tile_ptr(ksfc_type))
    ALLOCATE(field%lwflxsfc_tile_ptr(ksfc_type))
    ALLOCATE(field%evap_tile_ptr(ksfc_type))
    ALLOCATE(field%lhflx_tile_ptr(ksfc_type))
    ALLOCATE(field%shflx_tile_ptr(ksfc_type))
    ALLOCATE(field%frac_tile_ptr(ksfc_type))

    DO jsfc = 1,ksfc_type

      CALL add_ref( field_list, prefix//'rsns_tile',                                   &
                  & prefix//'rsns_'//csfc(jsfc), field%swflxsfc_tile_ptr(jsfc)%p,      &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                &
                  & t_cf_var('rsns_'//csfc(jsfc), 'W m-2',                             &
                  &          'shortwave net flux at surface on tile '//csfc(jsfc),     &
                  &          datatype_flt),                                            &
                  & grib2_var(0,4,9, ibits, GRID_UNSTRUCTURED, GRID_CELL),             &
                  & lrestart=.FALSE., ldims=shape2d,                                   &
                  & lmiss=.TRUE., missval=cdimissval )

      CALL add_ref( field_list, prefix//'rlns_tile',                                   &
                  & prefix//'rlns_'//csfc(jsfc), field%lwflxsfc_tile_ptr(jsfc)%p,      &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                &
                  & t_cf_var('rlns_'//csfc(jsfc), 'W m-2',                             &
                  &          'longwave net flux at surface on tile '//csfc(jsfc),      &
                  &          datatype_flt),                                            &
                  & grib2_var(0,5,5, ibits, GRID_UNSTRUCTURED, GRID_CELL),             &
                  & lrestart=.FALSE., ldims=shape2d,                                   &
                  & lmiss=.TRUE., missval=cdimissval )

      CALL add_ref( field_list, prefix//'evspsbl_tile',                                &
                  & prefix//'evspsbl_'//csfc(jsfc), field%evap_tile_ptr(jsfc)%p,       &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                &
                  & t_cf_var('evspsbl_'//csfc(jsfc), 'kg m-2 s-1',                     &
                  &          'evaporation on tile '//csfc(jsfc), datatype_flt),        &
                  & grib2_var(0,1,6, ibits, GRID_UNSTRUCTURED, GRID_CELL),             &
                  & lrestart=.FALSE., ldims=shape2d,                                   &
                  & lmiss=.TRUE., missval=cdimissval )

      CALL add_ref( field_list, prefix//'hfls_tile',                                   &
                  & prefix//'hfls_'//csfc(jsfc), field%lhflx_tile_ptr(jsfc)%p,         &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                &
                  & t_cf_var('hfls_'//csfc(jsfc), 'W m-2',                             &
                  &          'latent heat flux on tile '//csfc(jsfc), datatype_flt),   &
                  & grib2_var(0,0,10, ibits, GRID_UNSTRUCTURED, GRID_CELL),            &
                  & lrestart=.FALSE., ldims=shape2d,                                   &
                  & lmiss=.TRUE., missval=cdimissval )

      CALL add_ref( field_list, prefix//'hfss_tile',                                   &
                  & prefix//'hfss_'//csfc(jsfc), field%shflx_tile_ptr(jsfc)%p,         &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                &
                  & t_cf_var('hfss_'//csfc(jsfc), 'W m-2',                             &
                  &          'sensible heat flux on tile '//csfc(jsfc),datatype_flt),  &
                  & grib2_var(0,0,11, ibits, GRID_UNSTRUCTURED, GRID_CELL),            &
                  & lrestart=.FALSE., ldims=shape2d,                                   &
                  & lmiss=.TRUE., missval=cdimissval )

      CALL add_ref( field_list, prefix//'frac_tile',                                   &
                  & prefix//'frac_'//csfc(jsfc), field%frac_tile_ptr(jsfc)%p,          &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                &
                  & t_cf_var('frac_'//csfc(jsfc), '%',                                 &
                  &          'surface fraction of tile '//csfc(jsfc),                  &
                  &          datatype_flt),                                            &
                  & grib2_var(255,255,255, ibits, GRID_UNSTRUCTURED, GRID_CELL),       &
                  &  ldims=shape2d, lmiss=.TRUE., missval=cdimissval )

    END DO

    !-----------------------------------------
    ! wind stress, grid box mean
    !-----------------------------------------

    CALL add_var( field_list, prefix//'tauu', field%u_stress        ,           &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                           &
                & t_cf_var('u_stress','N m-2','u-momentum flux at the surface', &
                &          datatype_flt),                                       &
                & grib2_var(0,2,17, ibits, GRID_UNSTRUCTURED, GRID_CELL),       &
                & ldims=shape2d,                                                &
                & lrestart = .FALSE.,                                           &
                & isteptype=TSTEP_INSTANT                                       )

    CALL add_var( field_list, prefix//'tauv', field%v_stress,                   &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                           &
                & t_cf_var('v_stress','N m-2','v-momentum flux at the surface', &
                &          datatype_flt),                                       &
                & grib2_var(0,2,18, ibits, GRID_UNSTRUCTURED, GRID_CELL),       &
                & ldims=shape2d,                                                &
                & lrestart = .FALSE.,                                           &
                & isteptype=TSTEP_INSTANT                                       )

    ! wind stress, instantaneous tile values 

    CALL add_var( field_list, prefix//'tauu_tile', field%u_stress_tile,         &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                           &
                & t_cf_var('u_stress_tile', 'N m-2',                            &
                &          'u-momentum flux at the surface on tiles',           &
                &          datatype_flt),                                       &
                & grib2_var(0,2,17, ibits, GRID_UNSTRUCTURED, GRID_CELL),       &
                & ldims=shapesfc, lmiss=.TRUE., missval=cdimissval,             &
                & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.          )

    CALL add_var( field_list, prefix//'tauv_tile', field%v_stress_tile,         &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                           &
                & t_cf_var('v_stress_tile', 'N m-2',                            &
                &          'v-momentum flux at the surface on tiles',           &
                &          datatype_flt),                                       &
                & grib2_var(0,2,18, ibits, GRID_UNSTRUCTURED, GRID_CELL),       &
                & ldims=shapesfc, lmiss=.TRUE., missval=cdimissval,             &
                & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.          )

    ALLOCATE(field%u_stress_tile_ptr(ksfc_type))
    ALLOCATE(field%v_stress_tile_ptr(ksfc_type))

    DO jsfc = 1,ksfc_type

      CALL add_ref( field_list, prefix//'tauu_tile',                                &
                  & prefix//'tauu_'//csfc(jsfc), field%u_stress_tile_ptr(jsfc)%p,   &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                             &
                  & t_cf_var('u_stress_'//csfc(jsfc), 'N m-2',                      &
                  &          'u-momentum flux at the surface on tile '//csfc(jsfc), &
                  &          datatype_flt),                                         &
                  & grib2_var(0,2,17, ibits, GRID_UNSTRUCTURED, GRID_CELL),         &
                  & lrestart=.FALSE., ldims=shape2d,                                &
                  & lmiss=.TRUE., missval=cdimissval )

      CALL add_ref( field_list, prefix//'tauv_tile',                                &
                  & prefix//'tauv_'//csfc(jsfc), field%v_stress_tile_ptr(jsfc)%p,   &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                             &
                  & t_cf_var('v_stress_'//csfc(jsfc), 'N m-2',                      &
                  &          'v-momentum flux at the surface on tile '//csfc(jsfc), &
                  &          datatype_flt),                                         &
                  & grib2_var(0,2,18, ibits, GRID_UNSTRUCTURED, GRID_CELL),         &
                  & lrestart=.FALSE., ldims=shape2d,                                &
                  & lmiss=.TRUE., missval=cdimissval )
    END DO

    !-----------------------------------------
    ! near surface diagnostics, grid box mean
    !-----------------------------------------

    CALL add_var( field_list, prefix//'co2mmr', field%co2mmr,                   &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                           &
                & t_cf_var('CO2 MR','kg kg-1','co2 mixing ratio',               &
                &          datatype_flt),                                       &
                & grib2_var(255,255,255, ibits, GRID_UNSTRUCTURED, GRID_CELL),  &
                & ldims=shape2d,                                                &
                & lrestart = .FALSE.,                                           &
                & isteptype=TSTEP_INSTANT                                       )

    CALL add_var( field_list, prefix//'sfcwind', field%sfcwind,                 &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                           &
                & t_cf_var('sfcwind','m s-1','10m windspeed',                   &
                &          datatype_flt),                                       &
                & grib2_var(0,2,1, ibits, GRID_UNSTRUCTURED, GRID_CELL),        &
                & ldims=shape2d,                                                &
                & lrestart = .FALSE.,                                           &
                & isteptype=TSTEP_INSTANT                                       )

    CALL add_var( field_list, prefix//'uas', field%uas,                         &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                           &
                & t_cf_var('uas','m s-1','zonal wind in 10m',                   &
                &          datatype_flt),                                       &
                & grib2_var(0,2,2, ibits, GRID_UNSTRUCTURED, GRID_CELL),        &
                & ldims=shape2d,                                                &
                & lrestart = .FALSE.,                                           &
                & isteptype=TSTEP_INSTANT                                       )

    CALL add_var( field_list, prefix//'vas', field%vas,                         &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                           &
                & t_cf_var('vas','m s-1','meridional wind in 10m',              &
                &          datatype_flt),                                       &
                & grib2_var(0,2,3, ibits, GRID_UNSTRUCTURED, GRID_CELL),        &
                & ldims=shape2d,                                                &
                & lrestart = .FALSE.,                                           &
                & isteptype=TSTEP_INSTANT                                       )

    CALL add_var( field_list, prefix//'tas', field%tas,                         &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                           &
                & t_cf_var('tas','K','temperature in 2m',                       &
                &          datatype_flt),                                       &
                & grib2_var(0,0,0, ibits, GRID_UNSTRUCTURED, GRID_CELL),        &
                & ldims=shape2d,                                                &
                & lrestart = .FALSE.,                                           &
                & isteptype=TSTEP_INSTANT                                       )

    CALL add_var( field_list, prefix//'dew2', field%dew2,                       &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                           &
                & t_cf_var('dew2','K','dew point temperature in 2m',            &
                &          datatype_flt),                                       &
                & grib2_var(0,0,6, ibits, GRID_UNSTRUCTURED, GRID_CELL),        &
                & ldims=shape2d,                                                &
                & lrestart = .FALSE.,                                           &
                & isteptype=TSTEP_INSTANT                                       )

    CALL add_var( field_list, prefix//'tasmax', field%tasmax,                   &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                           &
                & t_cf_var('tasmax','K','maximum 2m temperature',               &
                &          datatype_flt),                                       &
                & grib2_var(0,0,4, ibits, GRID_UNSTRUCTURED, GRID_CELL),        &
                & ldims=shape2d,                                                &
                & lrestart = .FALSE.,                                           &
                & initval = -99._wp, resetval = -99._wp,                        &
                & isteptype=TSTEP_MAX,                                          &
                & action_list=actions(new_action(ACTION_RESET,"P1D"))           )

    CALL add_var( field_list, prefix//'tasmin', field%tasmin,                   &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                           &
                & t_cf_var('tasmin','K','minimum 2m temperature',               &
                &          datatype_flt),                                       &
                & grib2_var(0,0,5, ibits, GRID_UNSTRUCTURED, GRID_CELL),        &
                & ldims=shape2d,                                                &
                & lrestart = .FALSE.,                                           &
                & initval = 999._wp, resetval = 999._wp,                        &
                & isteptype=TSTEP_MIN,                                          &
                & action_list=actions(new_action(ACTION_RESET,"P1D"))           )

    !--------------------------------------
    ! near surface diagnostics, tile values
    !--------------------------------------

    CALL add_var( field_list, prefix//'sfcwind_tile', field%sfcwind_tile,       &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                           &
                & t_cf_var('sfcwind_tile','m s-1','10m windspeed on tiles',     &
                &          datatype_flt),                                       &
                & grib2_var(0,2,1, ibits, GRID_UNSTRUCTURED, GRID_CELL),        &
                & ldims=shapesfc,                                               &
                & lcontainer=.TRUE., lrestart=.FALSE.,                          &
                & isteptype=TSTEP_INSTANT                                       )

    CALL add_var( field_list, prefix//'uas_tile', field%uas_tile,               &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                           &
                & t_cf_var('uas_tile','m s-1','zonal wind in 10m on tiles',     &
                &          datatype_flt),                                       &
                & grib2_var(0,2,2, ibits, GRID_UNSTRUCTURED, GRID_CELL),        &
                & ldims=shapesfc,                                               &
                & lcontainer=.TRUE., lrestart=.FALSE.,                          &
                & isteptype=TSTEP_INSTANT                                       )

    CALL add_var( field_list, prefix//'vas_tile', field%vas_tile,               &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                           &
                & t_cf_var('vas_tile','m s-1','meridional wind in 10m on tiles',&
                &          datatype_flt),                                       &
                & grib2_var(0,2,3, ibits, GRID_UNSTRUCTURED, GRID_CELL),        &
                & ldims=shapesfc,                                               &
                & lcontainer=.TRUE., lrestart=.FALSE.,                          &
                & isteptype=TSTEP_INSTANT                                       )

    CALL add_var( field_list, prefix//'tas_tile', field%tas_tile,               &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                           &
                & t_cf_var('tas_tile','K','temperature in 2m on tiles',         &
                &          datatype_flt),                                       &
                & grib2_var(0,0,0, ibits, GRID_UNSTRUCTURED, GRID_CELL),        &
                & ldims=shapesfc,                                               &
                & lcontainer=.TRUE., lrestart=.FALSE.,                          &
                & isteptype=TSTEP_INSTANT                                       )

    CALL add_var( field_list, prefix//'dew2_tile', field%dew2_tile,             &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                           &
                & t_cf_var('dew2_tile','K','dew point temperature in 2m on tiles',&
                &          datatype_flt),                                       &
                & grib2_var(0,0,6, ibits, GRID_UNSTRUCTURED, GRID_CELL),        &
                & ldims=shapesfc,                                               &
                & lcontainer=.TRUE., lrestart=.FALSE.,                          &
                & isteptype=TSTEP_INSTANT                                       )


    ALLOCATE(field%sfcwind_tile_ptr(ksfc_type))
    ALLOCATE(field%uas_tile_ptr(ksfc_type))
    ALLOCATE(field%vas_tile_ptr(ksfc_type))
    ALLOCATE(field%tas_tile_ptr(ksfc_type))
    ALLOCATE(field%dew2_tile_ptr(ksfc_type))

    DO jsfc = 1,ksfc_type

      CALL add_ref( field_list, prefix//'sfcwind_tile',                             &
                  & prefix//'sfcwind_'//csfc(jsfc), field%sfcwind_tile_ptr(jsfc)%p, &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                             &
                  & t_cf_var('sfcwind_'//csfc(jsfc), 'm s-1',                       &
                  &          '10m windspeed on tile '//csfc(jsfc),                  &
                  &          datatype_flt),                                         &
                  & grib2_var(0,2,1, ibits, GRID_UNSTRUCTURED, GRID_CELL),          &
                  & lrestart=.FALSE., ldims=shape2d,                                &
                  & lmiss=.TRUE., missval=cdimissval )

      CALL add_ref( field_list, prefix//'uas_tile',                                 &
                  & prefix//'uas_'//csfc(jsfc), field%uas_tile_ptr(jsfc)%p,         &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                             &
                  & t_cf_var('uas_'//csfc(jsfc), 'm s-1',                           &
                  &          'zonal wind in 10m on tile '//csfc(jsfc),              &
                  &          datatype_flt),                                         &
                  & grib2_var(0,2,2, ibits, GRID_UNSTRUCTURED, GRID_CELL),          &
                  & lrestart=.FALSE., ldims=shape2d,                                &
                  & lmiss=.TRUE., missval=cdimissval )

      CALL add_ref( field_list, prefix//'vas_tile',                                 &
                  & prefix//'vas_'//csfc(jsfc), field%vas_tile_ptr(jsfc)%p,         &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                             &
                  & t_cf_var('vas_'//csfc(jsfc), 'm s-1',                           &
                  &          'meridional wind in 10m on tile '//csfc(jsfc),         &
                  &          datatype_flt),                                         &
                  & grib2_var(0,2,3, ibits, GRID_UNSTRUCTURED, GRID_CELL),          &
                  & lrestart=.FALSE., ldims=shape2d,                                &
                  & lmiss=.TRUE., missval=cdimissval )

      CALL add_ref( field_list, prefix//'tas_tile',                                 &
                  & prefix//'tas_'//csfc(jsfc), field%tas_tile_ptr(jsfc)%p,         &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                             &
                  & t_cf_var('tas_'//csfc(jsfc), 'K',                               &
                  &          'temperature in 2m on tile '//csfc(jsfc),              &
                  &          datatype_flt),                                         &
                  & grib2_var(0,0,0, ibits, GRID_UNSTRUCTURED, GRID_CELL),          &
                  & lrestart=.FALSE., ldims=shape2d,                                &
                  & lmiss=.TRUE., missval=cdimissval )

      CALL add_ref( field_list, prefix//'dew2_tile',                                &
                  & prefix//'dew2_'//csfc(jsfc), field%dew2_tile_ptr(jsfc)%p,       &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                             &
                  & t_cf_var('dew2_'//csfc(jsfc), 'K',                              &
                  &          'dew point temperature in 2m on tile '//csfc(jsfc),    &
                  &          datatype_flt),                                         &
                  & grib2_var(0,0,6, ibits, GRID_UNSTRUCTURED, GRID_CELL),          &
                  & lrestart=.FALSE., ldims=shape2d,                                &
                  & lmiss=.TRUE., missval=cdimissval )

    END DO

    ! global diagnostics
    cf_desc    = t_cf_var('tas_gmean', 'K', 'global mean temperature at 2m', datatype_flt,'tas_gmean')
    grib2_desc = grib2_var(255,255,255, ibits, GRID_UNSTRUCTURED, GRID_LONLAT)
    CALL add_var( field_list, prefix//'tas_gmean', field%tas_gmean,              &
                & GRID_LONLAT, ZA_SURFACE, cf_desc, grib2_desc, &
                & lrestart = .FALSE., ldims=(/1/) )

    cf_desc    = t_cf_var('rsdt_gmean', 'W m-2', 'global mean toa incident shortwave radiation', datatype_flt,'rsdt_gmean')
    grib2_desc = grib2_var(255,255,255, ibits, GRID_UNSTRUCTURED, GRID_LONLAT)
    CALL add_var( field_list, prefix//'rsdt_gmean', field%rsdt_gmean,            &
                & GRID_LONLAT, ZA_SURFACE, cf_desc, grib2_desc, &
                & lrestart = .FALSE., ldims=(/1/) )

    cf_desc    = t_cf_var('rsut_gmean', 'W m-2', 'global mean toa outgoing shortwave radiation', datatype_flt,'rsut_gmean')
    grib2_desc = grib2_var(255,255,255, ibits, GRID_UNSTRUCTURED, GRID_LONLAT)
    CALL add_var( field_list, prefix//'rsut_gmean', field%rsut_gmean,            &
                & GRID_LONLAT, ZA_SURFACE, cf_desc, grib2_desc, &
                & lrestart = .FALSE., ldims=(/1/) )

    cf_desc    = t_cf_var('rlut_gmean', 'W m-2', 'global mean toa outgoing longwave radiation', datatype_flt,'rlut_gmean')
    grib2_desc = grib2_var(255,255,255, ibits, GRID_UNSTRUCTURED, GRID_LONLAT)
    CALL add_var( field_list, prefix//'rlut_gmean', field%rlut_gmean,            &
                & GRID_LONLAT, ZA_SURFACE, cf_desc, grib2_desc, &
                & lrestart = .FALSE., ldims=(/1/) )

    cf_desc    = t_cf_var('prec_gmean', 'kg m-2 s-1', 'global mean precipitation flux', datatype_flt,'prec_gmean')
    grib2_desc = grib2_var(255,255,255, ibits, GRID_UNSTRUCTURED, GRID_LONLAT)
    CALL add_var( field_list, prefix//'prec_gmean', field%prec_gmean,            &
                & GRID_LONLAT, ZA_SURFACE, cf_desc, grib2_desc, &
                & lrestart = .FALSE., ldims=(/1/) )

    cf_desc    = t_cf_var('evap_gmean', 'kg m-2 s-1', 'global mean evaporation flux', datatype_flt,'evap_gmean')
    grib2_desc = grib2_var(255,255,255, ibits, GRID_UNSTRUCTURED, GRID_LONLAT)
    CALL add_var( field_list, prefix//'evap_gmean', field%evap_gmean,            &
                & GRID_LONLAT, ZA_SURFACE, cf_desc, grib2_desc, &
                & lrestart = .FALSE., ldims=(/1/) )

!   derived variable
    cf_desc    = t_cf_var('radtop_gmean', 'W m-2', 'global mean toa total radiation', datatype_flt,'radtop_gmean')
    grib2_desc = grib2_var(255,255,255, ibits, GRID_UNSTRUCTURED, GRID_LONLAT)
    CALL add_var( field_list, prefix//'radtop_gmean', field%radtop_gmean,            &
                & GRID_LONLAT, ZA_SURFACE, cf_desc, grib2_desc, &
                & lrestart = .FALSE., ldims=(/1/) )
!   derived variable
    cf_desc    = t_cf_var('fwfoce_gmean', 'kg m-2 s-1', 'mean surface freshwater flux over ocean surface', &
                & datatype_flt,'fwfoce_gmean')
    grib2_desc = grib2_var(255,255,255, ibits, GRID_UNSTRUCTURED, GRID_LONLAT)
    CALL add_var( field_list, prefix//'fwfoce_gmean', field%fwfoce_gmean,            &
                & GRID_LONLAT, ZA_SURFACE, cf_desc, grib2_desc, &
                & lrestart = .FALSE., ldims=(/1/) )

! icefrc not allocated in atmosphere
!   cf_desc    = t_cf_var('icefrc_gmean', 'frac', 'global mean ice cover of grid box', datatype_flt,'icefrc_gmean')
!   grib2_desc = grib2_var(255,255,255, ibits, GRID_UNSTRUCTURED, GRID_LONLAT)
!   CALL add_var( field_list, prefix//'icefrc_gmean', field%icefrc_gmean,            &
!               & GRID_LONLAT, ZA_SURFACE, cf_desc, grib2_desc, &
!               & lrestart = .FALSE., ldims=(/1/) )

  END SUBROUTINE new_echam_phy_field_list
  !-------------
  !>
  !!
  !!
  SUBROUTINE new_echam_phy_tend_list( jg, kproma, klev, kblks, ktracer, &
                                    & ctracer, listname, prefix,        &
                                    & tend_list, tend )

    INTEGER,INTENT(IN) :: jg !> patch ID
    INTEGER,INTENT(IN) :: kproma, klev, kblks, ktracer  !< dimension sizes

    CHARACTER(len=*)              ,INTENT(IN) :: listname, prefix
    CHARACTER(len=MAX_CHAR_LENGTH),INTENT(IN) :: ctracer(ktracer) !< tracer acronyms

    TYPE(t_var_list)      ,INTENT(INOUT) :: tend_list
    TYPE(t_echam_phy_tend),INTENT(INOUT) :: tend

    ! Local variables

    CHARACTER(len=MAX_CHAR_LENGTH) :: varname
    LOGICAL :: contvar_is_in_output

    TYPE(t_cf_var)    ::    cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER :: shape2d(2), shape3d(3), shape_trc(4)
    INTEGER :: ibits, jtrc, tlen
    INTEGER :: datatype_flt
    !------------------------------

    ibits = DATATYPE_PACK16 ! "entropy" of horizontal slice

    IF ( lnetcdf_flt64_output ) THEN
      datatype_flt = DATATYPE_FLT64
    ELSE
      datatype_flt = DATATYPE_FLT32
    ENDIF

    shape2d   = (/kproma, kblks/)
    shape3d   = (/kproma, klev, kblks/)
    shape_trc = (/kproma, klev, kblks, ktracer/)

    CALL new_var_list( tend_list, listname, patch_id=jg )
    CALL default_var_list_settings( tend_list, lrestart=.FALSE. )

    !------------------------------
    ! Temperature tendencies
    !------------------------------
    ! &       tend% ta      (nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('temperature_tendency', 'K s-1',                               &
                &         'temperature tendency (cv)',                                   &
                &         datatype_flt)
    grib2_desc = grib2_var(0,0,255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( tend_list, prefix//'ta', tend%ta,                                      &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc,grib2_desc,ldims=shape3d,&
                & vert_interp=create_vert_interp_metadata(                               &
                &   vert_intp_type=vintp_types("P","Z","I"),                             &
                &   vert_intp_method=VINTP_METHOD_LIN,                                   &
                &   l_extrapol=.FALSE. ) )

    ! &       tend% ta_dyn  (nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('temperature_tendency_dyn', 'K s-1',                           &
                &         'temperature tendency due to  due to resolved dynamics (cv)',  &
                &         datatype_flt)
    grib2_desc = grib2_var(0,0,255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( tend_list, prefix//'ta_dyn', tend%  ta_dyn,                            &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc,grib2_desc,ldims=shape3d,&
                & vert_interp=create_vert_interp_metadata(                               &
                &   vert_intp_type=vintp_types("P","Z","I"),                             &
                &   vert_intp_method=VINTP_METHOD_LIN,                                   &
                &   l_extrapol=.FALSE. ) )

    ! &       tend% ta_phy  (nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('temperature_tendency_phy', 'K s-1',                           &
                &         'temperature tendency due to parameterized processes (cv)',    &
                &         datatype_flt)
    grib2_desc = grib2_var(0,0,255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( tend_list, prefix//'ta_phy', tend%  ta_phy,                            &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc,grib2_desc,ldims=shape3d,&
                & vert_interp=create_vert_interp_metadata(                               &
                &   vert_intp_type=vintp_types("P","Z","I"),                             &
                &   vert_intp_method=VINTP_METHOD_LIN,                                   &
                &   l_extrapol=.FALSE. ) )

    IF ( echam_phy_tc(jg)%dt_rad > dt_zero ) THEN
       !
       IF (is_variable_in_output(first_output_name_list, var_name=prefix//'ta_rsw')) THEN
          cf_desc    = t_cf_var('temperature_tendency_rsw', 'K s-1',                           &
                      &         'temperature tendency due to shortwave radiation (cp)',        &
                      &         datatype_flt)
          grib2_desc = grib2_var(0,0,255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
          CALL add_var( tend_list, prefix//'ta_rsw', tend%  ta_rsw,                            &
                      & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc,grib2_desc,ldims=shape3d,&
                      & vert_interp=create_vert_interp_metadata(                               &
                      &   vert_intp_type=vintp_types("P","Z","I"),                             &
                      &   vert_intp_method=VINTP_METHOD_LIN,                                   &
                      &   l_extrapol=.FALSE. ) )
       END IF
       !
       IF (is_variable_in_output(first_output_name_list, var_name=prefix//'ta_rlw')) THEN
          cf_desc    = t_cf_var('temperature_tendency_rlw', 'K s-1',                           &
                      &         'temperature tendency due to longwave radiation (cp)',         &
                      &         datatype_flt)
          grib2_desc = grib2_var(0,0,255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
          CALL add_var( tend_list, prefix//'ta_rlw', tend%  ta_rlw,                            &
                      & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc,grib2_desc,ldims=shape3d,&
                      & vert_interp=create_vert_interp_metadata(                               &
                      &   vert_intp_type=vintp_types("P","Z","I"),                             &
                      &   vert_intp_method=VINTP_METHOD_LIN,                                   &
                      &   l_extrapol=.FALSE. ) )
       END IF
       !
       IF (is_variable_in_output(first_output_name_list, var_name=prefix//'ta_rad')) THEN
          cf_desc    = t_cf_var('temperature_tendency_rad', 'K s-1',                           &
                      &         'temperature tendency due to radiation (cp)',                  &
                      &         datatype_flt)
          grib2_desc = grib2_var(0,0,255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
          CALL add_var( tend_list, prefix//'ta_rad', tend%  ta_rlw,                            &
                      & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc,grib2_desc,ldims=shape3d,&
                      & vert_interp=create_vert_interp_metadata(                               &
                      &   vert_intp_type=vintp_types("P","Z","I"),                             &
                      &   vert_intp_method=VINTP_METHOD_LIN,                                   &
                      &   l_extrapol=.FALSE. ) )
       END IF
       !
    END IF

    IF (is_variable_in_output(first_output_name_list, var_name=prefix//'ta_rlw_impl')) THEN
       cf_desc    = t_cf_var('temperature_tendency_rlw_impl', 'K s-1',                      &
                   &         'temperature tendency due to LW rad. due to implicit land surface temperature change (cp)', &
                   &         datatype_flt)
       grib2_desc = grib2_var(0,0,255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var( tend_list, prefix//'ta_rlw_impl', tend%  ta_rlw_impl,                  &
                   & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,               &
                   & ldims=(/kproma,kblks/))
    END IF

    IF ( echam_phy_tc(jg)%dt_cld > dt_zero ) THEN
       !
       IF (is_variable_in_output(first_output_name_list, var_name=prefix//'ta_cld')) THEN
          cf_desc    = t_cf_var('temperature_tendency_cloud', 'K s-1',                         &
                      &         'temperature tendency due to large scale cloud processes (cp)',&
                      &         datatype_flt)
          grib2_desc = grib2_var(0,0,255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
          CALL add_var( tend_list, prefix//'ta_cld', tend%  ta_cld,                            &
                      & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc,grib2_desc,ldims=shape3d,&
                      & vert_interp=create_vert_interp_metadata(                               &
                      &   vert_intp_type=vintp_types("P","Z","I"),                             &
                      &   vert_intp_method=VINTP_METHOD_LIN,                                   &
                      &   l_extrapol=.FALSE. ) )
       END IF
       !
    END IF

    IF ( echam_phy_tc(jg)%dt_cnv > dt_zero ) THEN
       !
       IF (is_variable_in_output(first_output_name_list, var_name=prefix//'ta_cnv')) THEN
          cf_desc    = t_cf_var('temperature_tendency_convective', 'K s-1',                    &
                      &         'temperature tendency due to convective cloud processes (cp)', &
                      &         datatype_flt)
          grib2_desc = grib2_var(0,0,255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
          CALL add_var( tend_list, prefix//'ta_cnv', tend%  ta_cnv,                            &
                      & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,cf_desc,grib2_desc,ldims=shape3d, &
                      & vert_interp=create_vert_interp_metadata(                               &
                      &   vert_intp_type=vintp_types("P","Z","I"),                             &
                      &   vert_intp_method=VINTP_METHOD_LIN,                                   &
                      &   l_extrapol=.FALSE. ) )
       END IF
       !
    END IF

    IF ( echam_phy_tc(jg)%dt_vdf > dt_zero ) THEN
       !
       IF (is_variable_in_output(first_output_name_list, var_name=prefix//'ta_vdf')) THEN
          cf_desc    = t_cf_var('temperature_tendency_turbulent', 'K s-1',                     &
                      &         'temperature tendency due to vertical diffusion (cp)',         &
                      &         datatype_flt)
          grib2_desc = grib2_var(0,0,255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
          CALL add_var( tend_list, prefix//'ta_vdf', tend%  ta_vdf,                            &
                      & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc,grib2_desc,ldims=shape3d,&
                      & vert_interp=create_vert_interp_metadata(                               &
                      &   vert_intp_type=vintp_types("P","Z","I"),                             &
                      &   vert_intp_method=VINTP_METHOD_LIN,                                   &
                      &   l_extrapol=.FALSE. ) )
       END IF
       !
       IF (is_variable_in_output(first_output_name_list, var_name=prefix//'ta_sfc')) THEN
          cf_desc    = t_cf_var('temperature_tendency_surface',   'K s-1',                     &
                      &         'temperature tendency due to surface porcesses (cp)',          &
                      &         datatype_flt)
          grib2_desc = grib2_var(0,0,255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
          CALL add_var( tend_list, prefix//'ta_sfc', tend%  ta_sfc,                            &
                      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                    &
                      & cf_desc, grib2_desc, ldims=shape2d )
       END IF
       !
    END IF

    IF ( echam_phy_tc(jg)%dt_gwd > dt_zero ) THEN
       !
       IF (is_variable_in_output(first_output_name_list, var_name=prefix//'ta_gwd')) THEN
          cf_desc    = t_cf_var('temperature_tendency_Hines_gw', 'K s-1',                      &
                      &         'temperature tendency due to non-orogr. gravity waves (cp)',   &
                      &         datatype_flt)
          grib2_desc = grib2_var(0,0,255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
          CALL add_var( tend_list, prefix//'ta_gwd', tend%  ta_gwd,                            &
                      & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc,grib2_desc,ldims=shape3d,&
                      & vert_interp=create_vert_interp_metadata(                               &
                      &   vert_intp_type=vintp_types("P","Z","I"),                             &
                      &   vert_intp_method=VINTP_METHOD_LIN,                                   &
                      &   l_extrapol=.FALSE. ) )
       END IF
       !
    END IF

    IF ( echam_phy_tc(jg)%dt_sso > dt_zero ) THEN
       !
       IF (echam_phy_tc(jg)%dt_sso > time_config%tc_dt_dyn(jg) .OR.                            &
         & is_variable_in_output(first_output_name_list, var_name=prefix//'ta_sso')) THEN
          cf_desc    = t_cf_var('temperature_tendency_sso', 'K s-1',                           &
                      &         'temperature tendency due to sub grid scale orography (cp)',   &
                      &         datatype_flt)
          grib2_desc = grib2_var(0,0,255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
          CALL add_var( tend_list, prefix//'ta_sso', tend%  ta_sso,                            &
                      & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc,grib2_desc,ldims=shape3d,&
                      & vert_interp=create_vert_interp_metadata(                               &
                      &   vert_intp_type=vintp_types("P","Z","I"),                             &
                      &   vert_intp_method=VINTP_METHOD_LIN,                                   &
                      &   l_extrapol=.FALSE. ) )
       END IF
       !
    END IF

    !------------------------------
    ! U-wind tendencies
    !------------------------------
    ! &       tend%    ua     (nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('u_wind_tendency', 'm s-2',                                    &
                &         'u-wind tendency',                                             &
                &         datatype_flt)
    grib2_desc = grib2_var(0,2,255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( tend_list, prefix//'ua', tend%ua,                                      &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc,grib2_desc,ldims=shape3d,&
                & vert_interp=create_vert_interp_metadata(                               &
                &   vert_intp_type=vintp_types("P","Z","I"),                             &
                &   vert_intp_method=VINTP_METHOD_LIN,                                   &
                &   l_extrapol=.FALSE. ) )

    ! &       tend%    ua_dyn (nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('u_wind_tendency_dyn', 'm s-2',                                &
                &         'u-wind tendency due to resolved dynamics',                    &
                &         datatype_flt)
    grib2_desc = grib2_var(0,2,255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( tend_list, prefix//'ua_dyn', tend%ua_dyn,                              &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc,grib2_desc,ldims=shape3d,&
                & vert_interp=create_vert_interp_metadata(                               &
                &   vert_intp_type=vintp_types("P","Z","I"),                             &
                &   vert_intp_method=VINTP_METHOD_LIN,                                   &
                &   l_extrapol=.FALSE. ) )

    ! &       tend%    ua_phy (nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('u_wind_tendency_phy', 'm s-2',                                &
                &         'u-wind tendency due to parameterized processes',              &
                &         datatype_flt)
    grib2_desc = grib2_var(0,2,255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( tend_list, prefix//'ua_phy', tend%ua_phy,                              &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc,grib2_desc,ldims=shape3d,&
                & vert_interp=create_vert_interp_metadata(                               &
                &   vert_intp_type=vintp_types("P","Z","I"),                             &
                &   vert_intp_method=VINTP_METHOD_LIN,                                   &
                &   l_extrapol=.FALSE. ) )

    IF ( echam_phy_tc(jg)%dt_cnv > dt_zero ) THEN
       !
       IF (echam_phy_tc(jg)%dt_cnv > time_config%tc_dt_dyn(jg) .OR.                            &
         & is_variable_in_output(first_output_name_list, var_name=prefix//'ua_cnv')) THEN
          cf_desc    = t_cf_var('u_wind_tendency_convective', 'm s-2',                         &
                      &         'u-wind tendency due to convective cloud processes',           &
                      &         datatype_flt)
          grib2_desc = grib2_var(0,2,255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
          CALL add_var( tend_list, prefix//'ua_cnv', tend%ua_cnv,                              &
                      & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc,grib2_desc,ldims=shape3d,&
                      & vert_interp=create_vert_interp_metadata(                               &
                      &   vert_intp_type=vintp_types("P","Z","I"),                             &
                      &   vert_intp_method=VINTP_METHOD_LIN,                                   &
                      &   l_extrapol=.FALSE. ) )
       END IF
       !
    END IF

    IF ( echam_phy_tc(jg)%dt_vdf > dt_zero ) THEN
       !
       IF (echam_phy_tc(jg)%dt_vdf > time_config%tc_dt_dyn(jg) .OR.                            &
         & is_variable_in_output(first_output_name_list, var_name=prefix//'ua_vdf')) THEN
          cf_desc    = t_cf_var('u_wind_tendency_turbulent', 'm s-2',                          &
                      &         'u-wind tendency due to vertical diffusion',                   &
                      &         datatype_flt)
          grib2_desc = grib2_var(0,2,255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
          CALL add_var( tend_list, prefix//'ua_vdf', tend%ua_vdf,                              &
                      & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc,grib2_desc,ldims=shape3d,&
                      & vert_interp=create_vert_interp_metadata(                               &
                      &   vert_intp_type=vintp_types("P","Z","I"),                             &
                      &   vert_intp_method=VINTP_METHOD_LIN,                                   &
                      &   l_extrapol=.FALSE. ) )
       END IF
       !
    END IF

    IF ( echam_phy_tc(jg)%dt_gwd > dt_zero ) THEN
       !
       IF (echam_phy_tc(jg)%dt_gwd > time_config%tc_dt_dyn(jg) .OR.                            &
         & is_variable_in_output(first_output_name_list, var_name=prefix//'ua_gwd')) THEN
          cf_desc    = t_cf_var('u_wind_tendency_nonoro_gw', 'm s-2',                          &
                      &         'u-wind tendency due to non-orographic gravity waves',         &
                      &         datatype_flt)
          grib2_desc = grib2_var(0,2,255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
          CALL add_var( tend_list, prefix//'ua_gwd', tend%ua_gwd,                              &
                      & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc,grib2_desc,ldims=shape3d,&
                      & vert_interp=create_vert_interp_metadata(                               &
                      &   vert_intp_type=vintp_types("P","Z","I"),                             &
                      &   vert_intp_method=VINTP_METHOD_LIN,                                   &
                      &   l_extrapol=.FALSE. ) )
       END IF
       !
    END IF

    IF ( echam_phy_tc(jg)%dt_sso > dt_zero ) THEN
       !
       IF (echam_phy_tc(jg)%dt_sso > time_config%tc_dt_dyn(jg) .OR.                            &
         & is_variable_in_output(first_output_name_list, var_name=prefix//'ua_sso')) THEN
          cf_desc    = t_cf_var('u_wind_tendency_sso', 'm s-2',                                &
                      &         'u-wind tendency due to sub grid scale orography',             &
                      &         datatype_flt)
          grib2_desc = grib2_var(0,2,255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
          CALL add_var( tend_list, prefix//'ua_sso', tend%ua_sso,                              &
                      & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc,grib2_desc,ldims=shape3d,&
                      & vert_interp=create_vert_interp_metadata(                               &
                      &   vert_intp_type=vintp_types("P","Z","I"),                             &
                      &   vert_intp_method=VINTP_METHOD_LIN,                                   &
                      &   l_extrapol=.FALSE. ) )
       END IF
       !
    END IF

    !------------------------------
    ! V-wind tendencies
    !------------------------------
    ! &       tend%    va     (nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('v_wind_tendency', 'm s-2',                                    &
                &         'v-wind tendency',                                             &
                &         datatype_flt)
    grib2_desc = grib2_var(0,2,255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( tend_list, prefix//'va', tend%va,                                      &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc,grib2_desc,ldims=shape3d,&
                & vert_interp=create_vert_interp_metadata(                               &
                &   vert_intp_type=vintp_types("P","Z","I"),                             &
                &   vert_intp_method=VINTP_METHOD_LIN,                                   &
                &   l_extrapol=.FALSE. ) )

    ! &       tend%    va_dyn (nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('v_wind_tendency_dyn', 'm s-2',                                &
                &         'v-wind tendency due to resolved dynamics',                    &
                &         datatype_flt)
    grib2_desc = grib2_var(0,2,255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( tend_list, prefix//'va_dyn', tend%va_dyn,                              &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc,grib2_desc,ldims=shape3d,&
                & vert_interp=create_vert_interp_metadata(                               &
                &   vert_intp_type=vintp_types("P","Z","I"),                             &
                &   vert_intp_method=VINTP_METHOD_LIN,                                   &
                &   l_extrapol=.FALSE. ) )

    ! &       tend%    va_phy (nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('v_wind_tendency_phy', 'm s-2',                                &
                &         'v-wind tendency due to parameterized processes',              &
                &         datatype_flt)
    grib2_desc = grib2_var(0,2,255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( tend_list, prefix//'va_phy', tend%va_phy,                              &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc,grib2_desc,ldims=shape3d,&
                & vert_interp=create_vert_interp_metadata(                               &
                &   vert_intp_type=vintp_types("P","Z","I"),                             &
                &   vert_intp_method=VINTP_METHOD_LIN,                                   &
                &   l_extrapol=.FALSE. ) )

    IF ( echam_phy_tc(jg)%dt_cnv > dt_zero ) THEN
       !
       IF (echam_phy_tc(jg)%dt_cnv > time_config%tc_dt_dyn(jg) .OR.                            &
         & is_variable_in_output(first_output_name_list, var_name=prefix//'va_cnv')) THEN
          cf_desc    = t_cf_var('v_wind_tendency', 'm s-2',                                    &
                      &         'v-wind tendency due to convective cloud processes',           &
                      &         datatype_flt)
          grib2_desc = grib2_var(0,2,255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
          CALL add_var( tend_list, prefix//'va_cnv', tend%va_cnv,                              &
                      & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc,grib2_desc,ldims=shape3d,&
                      & vert_interp=create_vert_interp_metadata(                               &
                      &   vert_intp_type=vintp_types("P","Z","I"),                             &
                      &   vert_intp_method=VINTP_METHOD_LIN,                                   &
                      &   l_extrapol=.FALSE. ) )
       END IF
       !
    END IF

    IF ( echam_phy_tc(jg)%dt_vdf > dt_zero ) THEN
       !
       IF (echam_phy_tc(jg)%dt_vdf > time_config%tc_dt_dyn(jg) .OR.                            &
         & is_variable_in_output(first_output_name_list, var_name=prefix//'va_vdf')) THEN
          cf_desc    = t_cf_var('v_wind_tendency_turbulent', 'm s-2',                          &
                      &         'v-wind tendency due to vertical diffusion',                   &
                      &         datatype_flt)
          grib2_desc = grib2_var(0,2,255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
          CALL add_var( tend_list, prefix//'va_vdf', tend%va_vdf,                              &
                      & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc,grib2_desc,ldims=shape3d,&
                      & vert_interp=create_vert_interp_metadata(                               &
                      &   vert_intp_type=vintp_types("P","Z","I"),                             &
                      &   vert_intp_method=VINTP_METHOD_LIN,                                   &
                      &   l_extrapol=.FALSE. ) )
       END IF
       !
    END IF

    IF ( echam_phy_tc(jg)%dt_gwd > dt_zero ) THEN
       !
       IF (echam_phy_tc(jg)%dt_gwd > time_config%tc_dt_dyn(jg) .OR.                            &
         & is_variable_in_output(first_output_name_list, var_name=prefix//'va_gwd')) THEN
          cf_desc    = t_cf_var('v_wind_tendency_Hines_gw', 'm s-2',                           &
                      &         'v-wind tendency due to non-orographic gravity waves',         &
                      &         datatype_flt)
          grib2_desc = grib2_var(0,2,255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
          CALL add_var( tend_list, prefix//'va_gwd', tend%va_gwd,                              &
                      & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc,grib2_desc,ldims=shape3d,&
                      & vert_interp=create_vert_interp_metadata(                               &
                      &   vert_intp_type=vintp_types("P","Z","I"),                             &
                      &   vert_intp_method=VINTP_METHOD_LIN,                                   &
                      &   l_extrapol=.FALSE. ) )
       END IF
       !
    END IF

    IF ( echam_phy_tc(jg)%dt_sso > dt_zero ) THEN
       !
       IF (echam_phy_tc(jg)%dt_sso > time_config%tc_dt_dyn(jg) .OR.                            &
         & is_variable_in_output(first_output_name_list, var_name=prefix//'va_sso')) THEN
          cf_desc    = t_cf_var('v_wind_tendency_sso', 'm s-2',                                &
                      &         'v-wind tendency due to sub grid scale orography',             &
                      &         datatype_flt)
          grib2_desc = grib2_var(0,2,255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
          CALL add_var( tend_list, prefix//'va_sso', tend%va_sso,                              &
                      & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc,grib2_desc,ldims=shape3d,&
                      & vert_interp=create_vert_interp_metadata(                               &
                      &   vert_intp_type=vintp_types("P","Z","I"),                             &
                      &   vert_intp_method=VINTP_METHOD_LIN,                                   &
                      &   l_extrapol=.FALSE. ) )
       END IF
       !
    END IF

    !-------------------
    ! Tracer tendencies
    !-------------------
    ! Tracer arrays for (model) internal use                                               

    CALL add_var( tend_list, prefix//'qtrc', tend%qtrc,                        &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                        &
                & t_cf_var('tend_qtrc', 'kg kg-1 s-1',                         &
                &          'tendency of mass mixing ratio of tracers',         &
                &          datatype_flt),                                      &
                & grib2_var(0,20,2, ibits, GRID_UNSTRUCTURED, GRID_CELL),      &
                & ldims = shape_trc,                                           &
                & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.         )
    ALLOCATE(tend% qtrc_ptr(ktracer))

    CALL add_var( tend_list, prefix//'qtrc_dyn', tend%qtrc_dyn,                &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                        &
                & t_cf_var('tend_qtrc_dyn', 'kg kg-1 s-1',                     &
                &          'tendency of mass mixing ratio of tracers '//       &
                &          'due to resolved dynamics',                         &
                &          datatype_flt),                                      &
                & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED,GRID_CELL),&
                & ldims = shape_trc,                                           &
                & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.         )
    ALLOCATE(tend% qtrc_dyn_ptr(ktracer))

    CALL add_var( tend_list, prefix//'qtrc_phy', tend%qtrc_phy,                &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                        &
                & t_cf_var('tend_qtrc_phy', 'kg kg-1 s-1',                     &
                &          'tendency of mass mixing ratio of tracers '//       &
                &          'due to parameterized processes',                   &
                &          datatype_flt),                                      &
                & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED,GRID_CELL),&
                & ldims = shape_trc,                                           &
                & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.         )
    ALLOCATE(tend% qtrc_phy_ptr(ktracer))

    IF ( echam_phy_tc(jg)%dt_cld > dt_zero ) THEN
       !
       contvar_is_in_output = .FALSE.
       DO jtrc = 1,ktracer
          tlen = LEN_TRIM(ctracer(jtrc))
          varname=prefix//'q'//ctracer(jtrc)(1:tlen)//'_cld'
          IF (is_variable_in_output(first_output_name_list, var_name=TRIM(varname))) THEN
             contvar_is_in_output = .TRUE.
          END IF
       END DO
       !
       IF (echam_phy_tc(jg)%dt_cld > time_config%tc_dt_dyn(jg) .OR.                  &
         & contvar_is_in_output) THEN
          CALL add_var( tend_list, prefix//'qtrc_cld', tend%qtrc_cld,                &
                      & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                        &
                      & t_cf_var('tend_qtrc_cld', 'kg kg-1 s-1',                     &
                      &          'tendency of mass mixing ratio of tracers '//       &
                      &          'due to large scale cloud processes',               &
                      &          datatype_flt),                                      &           
                      & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED,GRID_CELL),&
                      & ldims = shape_trc,                                           &
                      & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.         )
          ALLOCATE(tend% qtrc_cld_ptr(ktracer))
       END IF
       !
       DO jtrc = 1,ktracer
          !
          tlen = LEN_TRIM(ctracer(jtrc))
          varname=prefix//'q'//ctracer(jtrc)(1:tlen)//'_cld'
          !
          IF (is_variable_in_output(first_output_name_list, var_name=TRIM(varname))) THEN
             CALL add_ref( tend_list, prefix//'qtrc_cld',                                        &
                         & TRIM(varname), tend%qtrc_cld_ptr(jtrc)%p,                             &
                         & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                                 &
                         & t_cf_var('tend_q'//ctracer(jtrc)(1:tlen)//'_cld', 'kg kg-1 s-1',      &
                         &          'tendency of mass mixing ratio of tracer '//                 &
                         &          ctracer(jtrc)(1:tlen)//                                      &
                         &          ' due to large scale cloud processes',                       &
                         &          datatype_flt),                                               &
                         & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),        &
                         & ref_idx=jtrc, ldims=(/kproma,klev,kblks/),                            &
                         & vert_interp=create_vert_interp_metadata(                              &
                         &             vert_intp_type=vintp_types("P","Z","I"),                  &
                         &             vert_intp_method=VINTP_METHOD_LIN )                       )
          END IF
          !
       END DO
       !
    END IF

    IF ( echam_phy_tc(jg)%dt_cnv > dt_zero ) THEN
       !
       contvar_is_in_output = .FALSE.
       DO jtrc = 1,ktracer
          tlen = LEN_TRIM(ctracer(jtrc))
          varname=prefix//'q'//ctracer(jtrc)(1:tlen)//'_cnv'
          IF (is_variable_in_output(first_output_name_list, var_name=TRIM(varname))) THEN
             contvar_is_in_output = .TRUE.
          END IF
       END DO
       !
       IF (echam_phy_tc(jg)%dt_cnv > time_config%tc_dt_dyn(jg) .OR.                  &
         & contvar_is_in_output) THEN
          CALL add_var( tend_list, prefix//'qtrc_cnv', tend%qtrc_cnv,                &
                      & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                        &
                      & t_cf_var('tend_qtrc_cnv', 'kg kg-1 s-1',                     &
                      &          'tendency of mass mixing ratio of tracers '//       &
                      &          'due to convective cloud processes',                &
                      &          datatype_flt),                                      &           
                      & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED,GRID_CELL),&
                      & ldims = shape_trc,                                           &
                      & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.         )
          ALLOCATE(tend% qtrc_cnv_ptr(ktracer))
       END IF
       !
       DO jtrc = 1,ktracer
          !
          tlen = LEN_TRIM(ctracer(jtrc))
          varname=prefix//'q'//ctracer(jtrc)(1:tlen)//'_cnv'
          !
          IF (is_variable_in_output(first_output_name_list, var_name=TRIM(varname))) THEN
             CALL add_ref( tend_list, prefix//'qtrc_cnv',                                        &
                         & TRIM(varname), tend%qtrc_cnv_ptr(jtrc)%p,                             &
                         & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                                 &
                         & t_cf_var('tend_q'//ctracer(jtrc)(1:tlen)//'_cnv', 'kg kg-1 s-1',      &
                         &          'tendency of mass mixing ratio of tracer '//                 &
                         &          ctracer(jtrc)(1:tlen)//                                      &
                         &          ' due to convective cloud processes',                        &
                         &          datatype_flt),                                               &
                         & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),        &
                         & ref_idx=jtrc, ldims=(/kproma,klev,kblks/),                            &
                         & vert_interp=create_vert_interp_metadata(                              &
                         &             vert_intp_type=vintp_types("P","Z","I"),                  &
                         &             vert_intp_method=VINTP_METHOD_LIN )                       )
          END IF
          !
       END DO
       !
    END IF


    IF ( echam_phy_tc(jg)%dt_vdf > dt_zero ) THEN
       !
       contvar_is_in_output = .FALSE.
       DO jtrc = 1,ktracer
          tlen = LEN_TRIM(ctracer(jtrc))
          varname=prefix//'q'//ctracer(jtrc)(1:tlen)//'_vdf'
          IF (is_variable_in_output(first_output_name_list, var_name=TRIM(varname))) THEN
             contvar_is_in_output = .TRUE.
          END IF
       END DO
       !
       IF (echam_phy_tc(jg)%dt_vdf > time_config%tc_dt_dyn(jg) .OR.                  &
         & contvar_is_in_output) THEN
          CALL add_var( tend_list, prefix//'qtrc_vdf', tend%qtrc_vdf,                &
                      & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                        &
                      & t_cf_var('tend_qtrc_vdf', 'kg kg-1 s-1',                     &
                      &          'tendency of mass mixing ratio of tracers '//       &
                      &          'due to vertical diffusion',                        &
                      &          datatype_flt),                                      &           
                      & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED,GRID_CELL),&
                      & ldims = shape_trc,                                           &
                      & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.         )
          ALLOCATE(tend% qtrc_vdf_ptr(ktracer))
       END IF
       !
       DO jtrc = 1,ktracer
          !
          tlen = LEN_TRIM(ctracer(jtrc))
          varname=prefix//'q'//ctracer(jtrc)(1:tlen)//'_vdf'
          !
          IF (is_variable_in_output(first_output_name_list, var_name=TRIM(varname))) THEN
             CALL add_ref( tend_list, prefix//'qtrc_vdf',                                        &
                         & TRIM(varname), tend%qtrc_vdf_ptr(jtrc)%p,                             &
                         & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                                 &
                         & t_cf_var('tend_q'//ctracer(jtrc)(1:tlen)//'_vdf', 'kg kg-1 s-1',      &
                         &          'tendency of mass mixing ratio of tracer '//                 &
                         &          ctracer(jtrc)(1:tlen)//                                      &
                         &          ' due to vertical diffusion',                                &
                         &          datatype_flt),                                               &
                         & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),        &
                         & ref_idx=jtrc, ldims=(/kproma,klev,kblks/),                            &
                         & vert_interp=create_vert_interp_metadata(                              &
                         &             vert_intp_type=vintp_types("P","Z","I"),                  &
                         &             vert_intp_method=VINTP_METHOD_LIN )                       )
          END IF
          !
       END DO
       !
    END IF

    IF ( (echam_phy_tc(jg)%dt_mox > dt_zero) ) THEN
       !
       contvar_is_in_output = .FALSE.
       DO jtrc = 1,ktracer
          tlen = LEN_TRIM(ctracer(jtrc))
          varname=prefix//'q'//ctracer(jtrc)(1:tlen)//'_mox'
          IF (is_variable_in_output(first_output_name_list, var_name=TRIM(varname))) THEN
             contvar_is_in_output = .TRUE.
          END IF
       END DO
       !
       IF (echam_phy_tc(jg)%dt_mox > time_config%tc_dt_dyn(jg) .OR.                  &
         & contvar_is_in_output) THEN
          CALL add_var( tend_list, prefix//'qtrc_mox', tend%qtrc_mox,                &
                      & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                        &
                      & t_cf_var('tend_qtrc_mox', 'kg kg-1 s-1',                     &
                      &          'tendency of mass mixing ratio of tracers '//       &
                      &          'due to methane ox. and H2O photolysis',            &
                      &          datatype_flt),                                      &           
                      & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED,GRID_CELL),&
                      & ldims = shape_trc,                                           &
                      & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.         )
          ALLOCATE(tend% qtrc_mox_ptr(ktracer))
       END IF
       !
       DO jtrc = 1,ktracer
          !
          tlen = LEN_TRIM(ctracer(jtrc))
          varname=prefix//'q'//ctracer(jtrc)(1:tlen)//'_mox'
          !
          IF (is_variable_in_output(first_output_name_list, var_name=TRIM(varname))) THEN
             CALL add_ref( tend_list, prefix//'qtrc_mox',                                        &
                         & TRIM(varname), tend%qtrc_mox_ptr(jtrc)%p,                             &
                         & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                                 &
                         & t_cf_var('tend_q'//TRIM(ctracer(jtrc))//'_mox', 'kg kg-1 s-1',        &
                         &          'tendency of mass mixing ratio of tracer '//                 &
                         &          TRIM(ctracer(jtrc))//                                        &
                         &          ' due to methane oxidation and H2O photolysis',              &
                         &          datatype_flt),                                               &
                         & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),        &
                         & ref_idx=jtrc, ldims=(/kproma,klev,kblks/),                            &
                         & vert_interp=create_vert_interp_metadata(                              &
                         &             vert_intp_type=vintp_types("P","Z","I"),                  &
                         &             vert_intp_method=VINTP_METHOD_LIN )                       )
          END IF
          !
       END DO
       !
    END IF

    IF ( (echam_phy_tc(jg)%dt_car > dt_zero) ) THEN
       !
       contvar_is_in_output = .FALSE.
       DO jtrc = 1,ktracer
          tlen = LEN_TRIM(ctracer(jtrc))
          varname=prefix//'q'//ctracer(jtrc)(1:tlen)//'_car'
          IF (is_variable_in_output(first_output_name_list, var_name=TRIM(varname))) THEN
             contvar_is_in_output = .TRUE.
          END IF
       END DO
       !
       IF (echam_phy_tc(jg)%dt_car > time_config%tc_dt_dyn(jg) .OR.                  &
         & contvar_is_in_output) THEN
          CALL add_var( tend_list, prefix//'qtrc_car', tend%qtrc_car,                &
                      & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                        &
                      & t_cf_var('tend_qtrc_car', 'kg kg-1 s-1',                     &
                      &          'tendency of mass mixing ratio of tracers '//       &
                      &          'due to linearized ozone chemistry (Cariolle)',     &
                      &          datatype_flt),                                      &           
                      & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED,GRID_CELL),&
                      & ldims = shape_trc,                                           &
                      & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.         )
          ALLOCATE(tend% qtrc_car_ptr(ktracer))
       END IF
       !
       DO jtrc = 1,ktracer
          !
          tlen = LEN_TRIM(ctracer(jtrc))
          varname=prefix//'q'//ctracer(jtrc)(1:tlen)//'_car'
          !
          IF (is_variable_in_output(first_output_name_list, var_name=TRIM(varname))) THEN
             CALL add_ref( tend_list, prefix//'qtrc_car',                                        &
                         & TRIM(varname), tend%qtrc_car_ptr(jtrc)%p,                             &
                         & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                                 &
                         & t_cf_var('tend_q'//TRIM(ctracer(jtrc))//'_car', 'kg kg-1 s-1',        &
                         &          'tendency of mass mixing ratio of tracer '//                 &
                         &          TRIM(ctracer(jtrc))//                                        &
                         &          ' due to linearized ozone chemistry (Cariolle)',             &
                         &          datatype_flt),                                               &
                         & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),        &
                         & ref_idx=jtrc, ldims=(/kproma,klev,kblks/),                            &
                         & vert_interp=create_vert_interp_metadata(                              &
                         &             vert_intp_type=vintp_types("P","Z","I"),                  &
                         &             vert_intp_method=VINTP_METHOD_LIN )                       )
          END IF
          !
       END DO
       !
    END IF

    CALL add_var( tend_list, prefix//'mtrc_phy', tend%mtrc_phy,                &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                        &
                & t_cf_var('tend_mtrc_phy', 'kg m-2 s-1',                      &
                &          'tendency of tracer mass '//                        &
                &          'due to parameterized processes',                   &
                &          datatype_flt),                                      &
                & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED,GRID_CELL),&
                & ldims = shape_trc,                                           &
                & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.         )
    ALLOCATE(tend% mtrc_phy_ptr(ktracer))

    CALL add_var( tend_list, prefix//'mtrcvi_phy', tend%mtrcvi_phy,            &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                          &
                & t_cf_var('tend_mtrcvi_phy', 'kg m-2 s-1',                    &
                &          'tendency of path of tracers '//                    &
                &          'due to parameterized processes',                   &
                &          datatype_flt),                                      &
                & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED,GRID_CELL),&
                & ldims = (/kproma,kblks,ktracer/),                            &
                & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.         )
    ALLOCATE(tend% mtrcvi_phy_ptr(ktracer))

    ! Referrence to individual tracer, for I/O

    DO jtrc = 1,ktracer
      tlen = LEN_TRIM(ctracer(jtrc))
      CALL add_ref( tend_list, prefix//'qtrc',                                            &
                  & prefix//'q'//ctracer(jtrc)(1:tlen), tend%qtrc_ptr(jtrc)%p,            &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                                 &
                  & t_cf_var('tend_q'//ctracer(jtrc)(1:tlen), 'kg kg-1 s-1',              &
                  &          'tendency of mass mixing ratio of tracer '//                 &
                  &          ctracer(jtrc)(1:tlen),                                       &
                  &          datatype_flt),                                               &
                  & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),        &
                  & ref_idx=jtrc, ldims=(/kproma,klev,kblks/),                            &
                  & vert_interp=create_vert_interp_metadata(                              &
                  &             vert_intp_type=vintp_types("P","Z","I"),                  &
                  &             vert_intp_method=VINTP_METHOD_LIN )                       )

      CALL add_ref( tend_list, prefix//'qtrc_dyn',                                        &
                  & prefix//'q'//ctracer(jtrc)(1:tlen)//'_dyn', tend%qtrc_dyn_ptr(jtrc)%p,&
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                                 &
                  & t_cf_var('tend_q'//ctracer(jtrc)(1:tlen)//'_dyn', 'kg kg-1 s-1',      &
                  &          'tendency of mass mixing ratio of tracer '//                 &
                  &          ctracer(jtrc)(1:tlen)//                                      &
                  &          ' due to resolved dynamics',                                 &
                  &          datatype_flt),                                               &
                  & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),        &
                  & ref_idx=jtrc, ldims=(/kproma,klev,kblks/),                            &
                  & vert_interp=create_vert_interp_metadata(                              &
                  &             vert_intp_type=vintp_types("P","Z","I"),                  &
                  &             vert_intp_method=VINTP_METHOD_LIN )                       )

      CALL add_ref( tend_list, prefix//'qtrc_phy',                                        &
                  & prefix//'q'//ctracer(jtrc)(1:tlen)//'_phy', tend%qtrc_phy_ptr(jtrc)%p,&
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                                 &
                  & t_cf_var('tend_q'//ctracer(jtrc)(1:tlen)//'_phy', 'kg kg-1 s-1',      &
                  &          'tendency of mass mixing ratio of tracer '//                 &
                  &          ctracer(jtrc)(1:tlen)//                                      &
                  &          ' due to parameterized processes',                           &
                  &          datatype_flt),                                               &
                  & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),        &
                  & ref_idx=jtrc, ldims=(/kproma,klev,kblks/),                            &
                  & vert_interp=create_vert_interp_metadata(                              &
                  &             vert_intp_type=vintp_types("P","Z","I"),                  &
                  &             vert_intp_method=VINTP_METHOD_LIN )                       )

      CALL add_ref( tend_list, prefix//'mtrc_phy',                                        &
                  & prefix//'m'//ctracer(jtrc)(1:tlen)//'_phy', tend%mtrc_phy_ptr(jtrc)%p,&
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                                 &
                  & t_cf_var('tend_m'//ctracer(jtrc)(1:tlen)//'_phy', 'kg m-2 s-1',       &
                  &          'tendency of '//ctracer(jtrc)(1:tlen)//                      &
                  &          ' mass due to parameterized processes',                      &
                  &          datatype_flt),                                               &
                  & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),        &
                  & ref_idx=jtrc, ldims=(/kproma,klev,kblks/),                            &
                  & vert_interp=create_vert_interp_metadata(                              &
                  &             vert_intp_type=vintp_types("P","Z","I"),                  &
                  &             vert_intp_method=VINTP_METHOD_LIN )                       )

      CALL add_ref( tend_list, prefix//'mtrcvi_phy',                                      &
                  & prefix//'m'//ctracer(jtrc)(1:tlen)//'vi_phy', tend%mtrcvi_phy_ptr(jtrc)%p, &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                   &
                  & t_cf_var('tend_m'//ctracer(jtrc)(1:tlen)//'vi_phy', 'kg m-2 s-1',     &
                  &          'tendency of '//ctracer(jtrc)(1:tlen)//                      &
                  &          ' path due to parameterized processes',                      &
                  &          datatype_flt),                                               &
                  & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),        &
                  & ref_idx=jtrc, ldims=(/kproma,kblks/)                                  )

    END DO

  END SUBROUTINE new_echam_phy_tend_list
  !-------------

END MODULE mo_echam_phy_memory
