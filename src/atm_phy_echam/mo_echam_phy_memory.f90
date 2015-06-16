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
!! @author Hui Wan (MPI-M)
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

  USE mo_kind,                ONLY: wp
  USE mo_impl_constants,      ONLY: SUCCESS, MAX_CHAR_LENGTH,  & 
    &                               VINTP_METHOD_PRES,         &
    &                               VINTP_METHOD_LIN,          &
    &                               VINTP_METHOD_LIN_NLEVP1
  USE mo_exception,           ONLY: message, finish
  USE mo_parallel_config,     ONLY: nproma
  USE mo_echam_sfc_indices,   ONLY: nsfc_type, csfc
  USE mo_model_domain,        ONLY: t_patch

  USE mo_linked_list,         ONLY: t_var_list
  USE mo_var_list,            ONLY: default_var_list_settings, &
    &                               add_var, add_ref,          &
    &                               new_var_list,              &
    &                               delete_var_list
  USE mo_var_metadata,        ONLY: create_vert_interp_metadata, vintp_types
  USE mo_cf_convention,       ONLY: t_cf_var
  USE mo_grib2,               ONLY: t_grib2_var
  USE mo_cdi_constants,       ONLY: GRID_REFERENCE,                    &
    &                               GRID_UNSTRUCTURED_CELL, GRID_CELL, &
    &                               ZA_HYBRID, ZA_HYBRID_HALF,         &
    &                               ZA_SURFACE, ZA_GENERIC_ICE,        &
    &                               DATATYPE_PACK16, DATATYPE_PACK24,  &
    &                               DATATYPE_FLT32,                    &
    &                               TSTEP_INSTANT, TSTEP_AVG
  USE mo_sea_ice_nml,         ONLY: kice


  IMPLICIT NONE
  PRIVATE

  PUBLIC :: prm_field, prm_tend                         !< variables
  PUBLIC :: prm_field_list, prm_tend_list               !< variable lists
  PUBLIC :: construct_echam_phy_state                   !< subroutine
  PUBLIC :: destruct_echam_phy_state                    !< subroutines
  PUBLIC :: t_echam_phy_field, t_echam_phy_tend         !< derived types

#ifdef HAVE_F95
  PUBLIC :: t_ptr2d, t_ptr3d
#endif
  CHARACTER(len=*), PARAMETER :: thismodule = 'mo_echam_phy_memory'

  !!--------------------------------------------------------------------------
  !!                               DATA TYPES
  !!--------------------------------------------------------------------------
  !>
  !! Derived data types for building pointer arrays
  !!
  TYPE t_ptr2d
    REAL(wp),POINTER :: p(:,:)    ! pointer to 2D (spatial) array
  END TYPE t_ptr2d

  TYPE t_ptr3d
    REAL(wp),POINTER :: p(:,:,:)  ! pointer to 3D (spatial) array
  END TYPE t_ptr3d

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

    ! Meteorology and tracers
    REAL(wp),POINTER ::     &
      & u         (:,:,:),  &!< [m/s]   zonal wind
      & v         (:,:,:),  &!< [m/s]   meridional wind
      & vor       (:,:,:),  &!< [1/s]   relative vorticity
      & temp      (:,:,:),  &!< [K]     temperature          (tm1  of memory_g1a in ECHAM)
      & tv        (:,:,:),  &!< [K]     virtual temperature  (tvm1 of memory_g1a in ECHAM)
      & q         (:,:,:,:),&!< [kg/kg] tracer concentration (qm1, xlm1, xim1 of memory_g1a in ECHAM)
      & qx        (:,:,:),  &!< [kg/kg] total concentration of hydrometeors
      & omega     (:,:,:),  &!< [Pa/s]  vertical velocity in pressure coord. ("vervel" in ECHAM)
      & geoi      (:,:,:),  &!< [m2/s2] geopotential at half levels (vertical interfaces)
      & geom      (:,:,:),  &!< [m2/s2] geopotential at full levels (layer ave. or mid-point value)
      & presi_old (:,:,:),  &!< [Pa]    pressure at half levels at time step "old"
      & presm_old (:,:,:),  &!< [Pa]    pressure at full levels at time step "old"
      & presi_new (:,:,:),  &!< [Pa]    pressure at half levels at time step "new"
      & presm_new (:,:,:)    !< [Pa]    pressure at full levels at time step "new"

    TYPE(t_ptr3d),ALLOCATABLE :: q_ptr(:)


    ! Radiation
    REAL(wp),POINTER ::       &
      !
      ! insolation at TOA
      & cosmu0      (:,  :),  &!< [ ]    cos of zenith angle mu0 for radiative heating  calculation
      & flxdwswtoa  (:,  :),  &!< [W/m2] downward shortwave flux at TOA
      !
      ! shortwave surface fluxes (updated every time step)
      & vissfc      (:,  :),  &!< [ ]    net shortwave radiation in VIS (positive downward)
      & nirsfc      (:,  :),  &!< [ ]    net shortwave radiation in NIR (positive downward)
      & parsfcdn    (:,  :),  &!< [ ]    downward shortwave radiation in PAR
      !
      ! fractions of net surface shortwave flux (updated every radiation time step)
      & visfrcsfc   (:,  :),  &!< [ ]    visible fraction of net surface shortwave flux 
      !
      ! fractions of diffuse shortwave surface radiation (updated every radiation time step)
      & visdffsfc   (:,  :),  &!< [ ]    diffuse fraction in VIS downward flux
      & nirdffsfc   (:,  :),  &!< [ ]    diffuse fraction in NIR downward flux
      & pardffsfc   (:,  :),  &!< [ ]    diffuse fraction in PAR downward flux

      ! shortwave net transmissivity (updated at radiation time steps)
      & swtrmclr    (:,:,:),  &!< [ ]    net shortwave transmissivity  , clear-sky, positive downward
      & swtrmall    (:,:,:),  &!< [ ]    net shortwave transmissivity  , all-sky,   positive downward
      & partrmdnsfc (:,  :),  &!< [ ]    downward shortwave transmissivity in PAR, all-sky
      !
      ! longwave net fluxes (updated at radiation time steps)
      & lwflxclr    (:,:,:),  &!< [W/m2] net longwave flux,            clear-sky, positive downward
      & lwflxall    (:,:,:),  &!< [W/m2] net longwave flux,            all-sky,   positive downward
      & lwflxupsfc  (:,  :),  &!< [W/m2] upward longwave flux at surface, all-sky
      !
      & o3          (:,:,:)    !< temporary set ozone mass mixing ratio  
    ! aerosol optical properties
    REAL(wp), POINTER ::      &
      & aer_aod_533 (:,:,:),  &!< aerosol optical depth at 533 nm
      & aer_ssa_533 (:,:,:),  &!< aerosol single scattering albedo at 533 nm
      & aer_asy_533 (:,:,:),  &!< aerosol asymmetry factor at 533 nm
      & aer_aod_2325(:,:,:),  &!< aerosol optical depth at 2325 nm
      & aer_ssa_2325(:,:,:),  &!< aerosol single scattering albedo at 2325 nm
      & aer_asy_2325(:,:,:),  &!< aerosol asymmetry factor at 2325 nm
      & aer_aod_9731(:,:,:)    !< effective aerosol optical depth at 9731 nm
            !< the last quantity is in the thermal wavelength ranch, 
            !< the first lie in the solar spectrum
    ! Cloud and precipitation
    REAL(wp),POINTER ::     &
      & aclc      (:,:,:),  &!< [m2/m2] cloud area fractional
      & aclcov    (:,  :),  &!< [m2/m2] total cloud cover
      & acdnc     (:,:,:),  &!< cloud droplet number concentration [1/m**3]
      & xvar      (:,:,:),  &!< variance of total water amount qv+qi+ql [kg/kg] (memory_g3b)
      & xskew     (:,:,:),  &!< skewness of total water amount qv+qi+ql [kg/kg]
      & relhum    (:,:,:),  &!< relative humidity (relhum of memory_g3b in ECHAM)
      & rsfl      (:,  :),  &!< sfc rain flux, large scale [kg m-2 s-1]
      & rsfc      (:,  :),  &!< sfc rain flux, convective  [kg m-2 s-1]
      & ssfl      (:,  :),  &!< sfc snow flux, large scale [kg m-2 s-1]
      & ssfc      (:,  :),  &!< sfc snow flux, convective  [kg m-2 s-1]
      & totprec   (:,  :),  &!< total precipitation flux,[kg m-2 s-1]
      & totprec_avg(:, :),  &!< (time ave)  total precipitation flux,[kg m-2 s-1]
      & qvi       (:,  :),  &!< vertically integrated water vapor [kg/m**2]
      & xlvi      (:,  :),  &!< vertically integrated cloud water [kg/m**2]
      & xivi      (:,  :)    !< vertically integrated cloud ice   [kg/m**2]

    REAL(wp),POINTER :: &
      & rintop (:,  :),     &!< low lever inversion, computed by "cover" (memory_g3b)
      & rtype  (:,  :),     &!< type of convection 0...3. (in memory_g3b in ECHAM)
      & topmax (:,  :),     &!< maximum height of convective cloud tops [Pa] (memory_g3b)
      & thvsig (:,  :)       !< Std. dev. of virtual potential temperature at the upper
                             !< interface of the lowest model layer.
                             !< Computed in "vdiff" by getting the square root of
                             !< thvvar(:,nlev-1,:). Used by "cucall".

    REAL(wp),POINTER :: &
      & siced  (:,  :),     &!< ice depth
      & alake  (:,  :),     &!< lake mask
      & alb    (:,  :),     &!< surface background albedo
      & seaice (:,  :)       !< sea ice as read in from amip input

    ! Energy and moisture budget related diagnostic variables
    REAL(wp),POINTER :: &
      & sh_vdiff (:,  :),   &!< sensible heat flux of vdiff
      & qv_vdiff (:,  :),   &!< qv flux of vdiff
      & ch_concloud (:,  :),&!< condensational heating, convection, large scale clouds
      & con_dtrl (:,  :),   &!< detrainment of liquid from convection
      & con_dtri (:,  :),   &!< detrainment of ice from convection 
      & con_iteqv (:,  :),  &!< v. int. tendency of water vapor within convection
      & cld_dtrl (:,  :),   &!< entrainment of liquid from convection to cloud
      & cld_dtri (:,  :),   &!< entrainment of ice from convection to cloud
      & cld_iteq (:,  :)    !< v. int. tendency of qv,qc, and qi within cloud

    ! orography
    REAL(wp),POINTER :: &
      & oromea (:,  :),     &!< Orographic mean elevation
      & orostd (:,  :),     &!< Orographic standard deviation
      & orosig (:,  :),     &!< Orographic slope
      & orogam (:,  :),     &!< Orographic anisotropy
      & orothe (:,  :),     &!< Orographic angle
      & oropic (:,  :),     &!< Orographic peacks elevation
      & oroval (:,  :)       !< Orographic valleys elevation

    ! JSBACH
    REAL(wp),POINTER :: &
      & tsfc_rad  (:,  :),  &!< [K] radiative sfc. temperature for use in radiation
      & tsfc_radt (:,  :),  &!< [K] radiative sfc. temperature at radiation time
      & csat      (:,  :),  &!<
      & cair      (:,  :)    !<

    ! Sea ice.
    ! See also atm_oce_lnd_interface/mo_sea_ice_types.f90
    INTEGER              :: kice  ! Number of ice-thickness classes
    REAL(wp),POINTER     ::     &
      & Tsurf   (:,:,:),        & ! Ice surface temperature [degC]
      & T1      (:,:,:),        & ! Temperature of upper ice layer [degC]
      & T2      (:,:,:),        & ! Temperature of lower ice layer [degC]
      & hi      (:,:,:),        & ! Ice thickness [m]
      & hs      (:,:,:),        & ! Snow thickness on ice [m]
      & Qtop    (:,:,:),        & ! Energy flux available for surface melting [W/m^2]
      & Qbot    (:,:,:),        & ! Energy flux at ice-ocean interface [W/m^2]
      & conc    (:,:,:),        & ! Ice concentration [0,1]
      & albvisdir_ice(:,:,:),   & ! Ice surface albedo for visible range, direct
      & albvisdif_ice(:,:,:),   & ! Ice surface albedo for visible range, diffuse
      & albnirdir_ice(:,:,:),   & ! Ice surface albedo for near IR range, direct
      & albnirdif_ice(:,:,:)      ! Ice surface albedo for near IR range, diffuse

    ! Orographic wave drag (ssodrag)

    REAL(wp),POINTER ::     &
      & u_stress_sso   (:,:),  &! < Zonal gravity wave stress
      & v_stress_sso   (:,:),  &! < Meridional gravity wave stress
      & dissipation_sso(:,:)    ! < Dissipation of orographic waves

    ! Turbulence

    REAL(wp),POINTER ::     &
      & tke       (:,:,:),  &!< turbulent kinetik energy at step n+1
      & tkem0     (:,:,:),  &!< turbulent kinetik energy at step n
      & tkem1     (:,:,:)    !< turbulent kinetik energy at step n-1

    ! need only for vdiff ++++
    REAL(wp),POINTER ::     &
      & ri        (:,:,:),  &!< moist Richardson number at layer interfaces
      & mixlen    (:,:,:),  &!< mixing length at layer interfaces
      & thvvar    (:,:,:)    !< variance of virtual potential temperature at layer interfaces.
                             !< Computed in "vdiff" by solving a prognostic equation of
                             !< the variance. Used for getting "thvsig".

    REAL(wp),POINTER ::      &
      & cfm     (:,:,:),     &!< turbulent exchange coefficient
      & cfm_tile(:,:,:),     &!< turbulent exchange coefficient
      & cfh     (:,:,:),     &!< turbulent exchange coefficient
      & cfh_tile(:,:,:),     &!< turbulent exchange coefficient
      & cfv     (:,:,:),     &!< turbulent exchange coefficient
      & cftke   (:,:,:),     &!< turbulent exchange coefficient
      & cfthv   (:,:,:)       !< turbulent exchange coefficient

    TYPE(t_ptr2d),ALLOCATABLE :: cfm_tile_ptr(:)
    TYPE(t_ptr2d),ALLOCATABLE :: cfh_tile_ptr(:)

    REAL(wp),POINTER ::     &
      & coriol(:,:),        &!< Coriolis parameter, needed for diagnosing PBL height.
      & ghpbl (:,:),        &!< geopotential of the top of the atmospheric boundary layer
      & z0m_tile(:,:,:),    &!< aerodynamic roughness length (over each surface type)
      & z0m   (:,:),        &!< aerodynamic roughness length (grid box mean)
      & z0h_lnd(:,:),       &!< roughness length for heat (over land)
      & ustar (:,:),        &!<
      & kedisp(:,:),        &!< time-mean (or integrated?) vertically integrated dissipation of kinetic energy
      & ocu   (:,:),        &!< eastward  velocity of ocean surface current
      & ocv   (:,:)          !< northward velocity of ocean surface current

      !
    REAL(wp),POINTER ::     &
      ! net fluxes at TOA and surface
      & swflxsfc    (:,:),     &!< [ W/m2] shortwave net flux at surface
      & swflxsfc_tile(:,:,:),  &!< [ W/m2] shortwave net flux at surface
      & lwflxsfc    (:,:),     &!< [ W/m2] longwave net flux at surface
      & lwupflxsfc    (:,:),   &!< [ W/m2] longwave upward flux at surface
      & lwflxsfc_tile(:,:,:),  &!< [ W/m2] longwave net flux at surface
      & dlwflxsfc_dT(:,:),     &!< [ W/m2/K] longwave net flux temp tend at surface
      & swflxtoa    (:,:),     &!< [ W/m2] shortwave net flux at TOA 
      & lwflxtoa    (:,:)       !< [ W/m2] shortwave net flux at TOA

    TYPE(t_ptr2d),ALLOCATABLE :: swflxsfc_tile_ptr(:)
    TYPE(t_ptr2d),ALLOCATABLE :: lwflxsfc_tile_ptr(:)

    TYPE(t_ptr2d),ALLOCATABLE :: z0m_tile_ptr(:)

    ! need only for vdiff ----

    ! Surface variables

    LOGICAL, POINTER :: &
      & lfland(:,:),        &!< .TRUE. when fraction of land > 0.
      & lfglac(:,:)          !< .TRUE. when fraction of glaciated land > 0.

    REAL(wp),POINTER :: &
      & lsmask(:,:),        &!< land-sea mask. (1. = land, 0. = sea/lakes) (slm in memory_g3b)
      & glac  (:,:),        &!< fraction of land covered by glaciers (glac in memory_g3b)
      & icefrc(:,:),        &!< ice cover given as the fraction of grid box (friac  in memory_g3b)
      & tsfc_tile (:,:,:),  &!< surface temperature over land/water/ice (tsw/l/i in memory_g3b)
      & tsfc      (:,  :),  &!< surface temperature, grid box mean
      & qs_sfc_tile(:,:,:)   !< saturation specitifc humidity at surface 

    TYPE(t_ptr2d),ALLOCATABLE ::   tsfc_tile_ptr(:)
    TYPE(t_ptr2d),ALLOCATABLE :: qs_sfc_tile_ptr(:)

    ! Surface albedo
    REAL(wp),POINTER :: &
      & albvisdir_tile (:,:,:),  &!< [ ] surface albedo over tiles for visible range, direct
      & albvisdif_tile (:,:,:),  &!< [ ] surface albedo over tiles for visible range, diffuse
      & albnirdir_tile (:,:,:),  &!< [ ] surface albedo over tiles for near-IR range, direct
      & albnirdif_tile (:,:,:),  &!< [ ] surface albedo over tiles for near-IR range, diffuse
      & albedo_tile    (:,:,:),  &!< [ ] surface albedo over tiles
      & albvisdir      (:,:  ),  &!< [ ] surface albedo for visible range, direct, grid-box mean
      & albvisdif      (:,:  ),  &!< [ ] surface albedo for visible range, diffuse, grid-box mean
      & albnirdir      (:,:  ),  &!< [ ] surface albedo for near-IR range, direct, grid-box mean
      & albnirdif      (:,:  ),  &!< [ ] surface albedo for near-IR range, diffuse, grid-box mean
      & albedo         (:,:  )    !< [ ] surface albedo, grid-box mean

    TYPE(t_ptr2d),ALLOCATABLE :: albvisdir_tile_ptr(:), albvisdif_tile_ptr(:), &
      & albnirdir_tile_ptr(:), albnirdif_tile_ptr(:), albedo_tile_ptr(:)

    REAL(wp),POINTER :: &
      & lhflx     (:,  :),    &!< grid box mean latent   heat flux at surface
      & shflx     (:,  :),    &!< grid box mean sensible heat flux at surface
      & evap      (:,  :),    &!< grid box mean evaporation at surface
      & lhflx_tile(:,:,:),    &!< latent   heat flux at surface on tiles
      & shflx_tile(:,:,:),    &!< sensible heat flux at surface on tiles
      & evap_tile(:,:,:),     &!< evaporation at surface on tiles
      & dshflx_dT_tile(:,:,:)  !< temp tendency of SHF at surface on tiles

    TYPE(t_ptr2d),ALLOCATABLE :: lhflx_tile_ptr(:)
    TYPE(t_ptr2d),ALLOCATABLE :: shflx_tile_ptr(:)
    TYPE(t_ptr2d),ALLOCATABLE :: evap_tile_ptr(:)
    TYPE(t_ptr2d),ALLOCATABLE :: dshflx_dT_tile_ptr(:)

    REAL(wp),POINTER :: &
      & u_stress     (:,  :), &!< grid box mean wind stress
      & v_stress     (:,  :), &!< grid box mean wind stress
      & u_stress_tile(:,:,:), &!< wind stress on tiles
      & v_stress_tile(:,:,:)   !< wind stress on tiles

    TYPE(t_ptr2d),ALLOCATABLE :: u_stress_tile_ptr(:)
    TYPE(t_ptr2d),ALLOCATABLE :: v_stress_tile_ptr(:)

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
      ! all processes
      !
      &    u     (:,:,:)  , & !< u-wind tendency
      &    v     (:,:,:)  , & !< v-wind tendency
      & temp     (:,:,:)  , & !< temperature tendency
      &    q     (:,:,:,:), & !< tracer tendency
      !
      ! resolved dynamcis
      !
      &    u_dyn (:,:,:)  , & !< u-wind tendency due to resolved dynamics
      &    v_dyn (:,:,:)  , & !< v-wind tendency due to resolved dynamics
      & temp_dyn (:,:,:)  , & !< temperature tendency due to resolved dynamics
      &    q_dyn (:,:,:,:), & !< tracer tendency due to resolved dynamics
      !
      ! all parameterized processes
      !
      &    u_phy (:,:,:)  , & !< u-wind tendency due to parameterized processes
      &    v_phy (:,:,:)  , & !< v-wind tendency due to parameterized processes
      & temp_phy (:,:,:)  , & !< temperature tendency due to parameterized processes
      &    q_phy (:,:,:,:), & !< tracer tendency due to parameterized processes
      !
      ! cloud microphysics
      !
      & temp_cld (:,:,:)  , & !< temperature tendency due to large scale cloud processes
      &    q_cld (:,:,:,:), & !< tracer tendency  due to large scale cloud processes
      !
      ! cumulus convection
      !
      & temp_cnv (:,:,:),   & !< temperature tendency due to convective cloud processes
      &    u_cnv (:,:,:),   & !< u-wind tendency due to convective cloud processes
      &    v_cnv (:,:,:),   & !< v-wind tendency due to convective cloud processes
      &    q_cnv (:,:,:,:), & !< tracer tendency due to convective cloud processes
      &    xl_dtr(:,:,:),   & !< cloud water tendency due to detrainment from convective clouds
      &    xi_dtr(:,:,:),   & !< cloud ice tendency due to detrainment from convective clouds
      !
      ! vertical diffusion ("vdiff")
      !
      & temp_vdf (:,:,:)  , & !< temperature tendency due to vertical diffusion
      &    u_vdf (:,:,:)  , & !< u-wind tendency due to vertical diffusion
      &    v_vdf (:,:,:)  , & !< v-wind tendency due to vertical diffusion
      &    q_vdf (:,:,:,:), & !< tracer tendency due to vertical diffusion
      !
      ! Hines param. for atmospheric gravity waves
      !
      & u_gwh    (:,:,:)  , & !< u-wind tendency due to non-orographic gravity waves
      & v_gwh    (:,:,:)  , & !< v-wind tendency due to non-orographic gravity waves
      & temp_gwh (:,:,:)  , & !< temperature tendency due to non-orographic gravity waves
      !
      ! subgrid scale orographic (sso) blocking and gravity wave drag
      !
      & u_sso    (:,:,:)  , & !< u-wind tendency due to sub grid scale orography
      & v_sso    (:,:,:)  , & !< v-wind tendency due to sub grid scale orography
      & temp_sso (:,:,:)  , & !< temperature tendency due to sub grid scale orography
      !
      ! radiation
      !
      & temp_rsw (:,:,:)  , & !< temperature due to shortwave radiation
      & temp_rlw (:,:,:)  , & !< temperature due to longwave radiation
      & temp_rlw_impl(:,:)    !< temperature tendency due to LW rad. due to implicit land surface temperature change

    TYPE(t_ptr3d),ALLOCATABLE ::     q_ptr(:)
    TYPE(t_ptr3d),ALLOCATABLE :: q_dyn_ptr(:)
    TYPE(t_ptr3d),ALLOCATABLE :: q_phy_ptr(:)
    TYPE(t_ptr3d),ALLOCATABLE :: q_cld_ptr(:)
    TYPE(t_ptr3d),ALLOCATABLE :: q_cnv_ptr(:)
    TYPE(t_ptr3d),ALLOCATABLE :: q_vdf_ptr(:)

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
 
CONTAINS


  !!--------------------------------------------------------------------------
  !!                SUBROUTINES FOR BUILDING AND DELETING VARIABLE LISTS 
  !!--------------------------------------------------------------------------
  !>
  !! Top-level procedure for building the physics state
  !!
  SUBROUTINE construct_echam_phy_state( ntracer, patch_array )

    INTEGER,INTENT(IN) :: ntracer
    TYPE(t_patch),INTENT(IN) :: patch_array(:)
    CHARACTER(len=MAX_CHAR_LENGTH) :: listname
    CHARACTER(len=MAX_CHAR_LENGTH) :: ctracer(ntracer) !< tracer acronyms
    INTEGER :: ndomain, jg, ist, nblks, nlev

    !---

    CALL message(TRIM(thismodule),'Construction of ECHAM physics state started.')

    ! Stop if ntracer/=3 (to be generalized later)
    IF (ntracer /= 3) CALL finish( TRIM(thismodule)//'construct_echam_phy_state', &
      &                            'Currently does not work for ntracer /= 3'     )

    ctracer(1:3) = (/'hus','clw','cli'/)

    ! Allocate pointer arrays prm_field and prm_tend, 
    ! as well as the corresponding list arrays.

    ndomain = SIZE(patch_array)

    ALLOCATE( prm_field(ndomain), prm_tend(ndomain), STAT=ist)
    IF (ist/=SUCCESS) CALL finish(TRIM(thismodule), &
      &'allocation of prm_field/tend array failed')

    ALLOCATE( prm_field_list(ndomain), prm_tend_list(ndomain), STAT=ist)
    IF (ist/=SUCCESS) CALL finish(TRIM(thismodule), &
      &'allocation of prm_field/tend list array failed')

    ! Build a field list and a tendency list for each grid level.
    ! This includes memory allocation. 

    DO jg = 1,ndomain

      nblks = patch_array(jg)%nblks_c
      nlev  = patch_array(jg)%nlev

      WRITE(listname,'(a,i2.2)') 'prm_field_D',jg
      CALL new_echam_phy_field_list( jg, nproma, nlev, nblks, ntracer, ctracer, &
                                   & nsfc_type, TRIM(listname), 'prm_',         &
                                   & prm_field_list(jg), prm_field(jg)          )

      WRITE(listname,'(a,i2.2)') 'prm_tend_D',jg
      CALL new_echam_phy_tend_list( jg, nproma, nlev, nblks, ntracer, ctracer, &
                                  & TRIM(listname), 'prm_tend_',               &
                                  & prm_tend_list(jg), prm_tend(jg)            )
    ENDDO
    CALL message(TRIM(thismodule),'Construction of ECHAM physics state finished.')

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
    CALL message(TRIM(thismodule),'Destruction of ECHAM physics state started.')

    ndomain = SIZE(prm_field)

    DO jg = 1,ndomain
      CALL delete_var_list( prm_field_list(jg) )
      CALL delete_var_list( prm_tend_list (jg) )
    ENDDO

    DEALLOCATE( prm_field_list, prm_tend_list, STAT=ist )
    IF (ist/=SUCCESS) CALL finish(TRIM(thismodule), &
      & 'deallocation of prm_field/tend list array failed')

    DEALLOCATE( prm_field, prm_tend, STAT=ist )
    IF (ist/=SUCCESS) CALL finish(TRIM(thismodule), &
      & 'deallocation of prm_field/tend array failed')

    CALL message(TRIM(thismodule),'Destruction of ECHAM physics state finished.')

  END SUBROUTINE destruct_echam_phy_state
  !--------------------------------------------------------------------


  !--------------------------------------------------------------------
  !>
  SUBROUTINE new_echam_phy_field_list( k_jg, kproma, klev, kblks, ktracer, &
                                     & ctracer, ksfc_type, listname,       &
                                     & prefix, field_list, field           )

    INTEGER,INTENT(IN) :: k_jg !> patch ID
    INTEGER,INTENT(IN) :: kproma, klev, kblks, ktracer, ksfc_type  !< dimension sizes

    CHARACTER(len=*)              ,INTENT(IN) :: listname, prefix
    CHARACTER(len=MAX_CHAR_LENGTH),INTENT(IN) :: ctracer(ktracer) !< tracer acronyms


    TYPE(t_var_list),       INTENT(INOUT) :: field_list
    TYPE(t_echam_phy_field),INTENT(INOUT) :: field

    ! Local variables

    TYPE(t_cf_var)    ::    cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER :: shape2d(2), shape3d(3), shapesfc(3), shapeice(3), shape3d_layer_interfaces(3)
!0!    INTEGER :: shape4d(4)
    INTEGER :: ibits, iextbits, jsfc, jtrc

    ibits = DATATYPE_PACK16
    iextbits = DATATYPE_PACK24

    shape2d  = (/kproma,       kblks/)
    shape3d  = (/kproma, klev, kblks/)
    shapesfc = (/kproma, kblks, ksfc_type/)
    shape3d_layer_interfaces = (/kproma,klev+1,kblks/)


    ! Register a field list and apply default settings

    CALL new_var_list( field_list, TRIM(listname), patch_id=k_jg )
    CALL default_var_list_settings( field_list,                &
                                  & lrestart=.TRUE.  )

    !------------------------------
    ! Meteorological quantities
    !------------------------------

    ! &       field% u         (nproma,nlev  ,nblks),          &
    cf_desc    = t_cf_var('eastward_wind', 'm s-1', 'u-component of wind', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 2, 2, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'ua', field%u,                                     &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d, &
                & vert_interp = &
                &   create_vert_interp_metadata( vert_intp_type=vintp_types("P","Z","I") ) )

    ! &       field% v         (nproma,nlev  ,nblks),          &
    cf_desc    = t_cf_var('northward_wind', 'm s-1', 'v-component of wind', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 2, 3, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'va', field%v,                                     &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d, &
                & vert_interp = &
                &   create_vert_interp_metadata( vert_intp_type=vintp_types("P","Z","I") ) )

    ! &       field% vor       (nproma,nlev  ,nblks),          &
    cf_desc    = t_cf_var('vorticity', 's-1', 'relative vorticity', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 2, 12, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'vor', field%vor,                         &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d, &
                & vert_interp=create_vert_interp_metadata(                      &
                &   vert_intp_type=vintp_types("P","Z","I"),                    &
                &   vert_intp_method=VINTP_METHOD_LIN ) )

    ! &       field% temp      (nproma,nlev  ,nblks),          &
    cf_desc    = t_cf_var('temperature', 'K', 'temperature', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 0, 0, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'ta', field%temp,                         &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d, &
                & vert_interp=create_vert_interp_metadata(                      &
                &   vert_intp_type=vintp_types("P","Z","I"),                    &
                &   vert_intp_method=VINTP_METHOD_LIN ) )

    ! &       field% tv        (nproma,nlev  ,nblks),          &
    cf_desc    = t_cf_var('virtual_temperature', 'K', 'virtual temperature', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0,0,1, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'tv', field%tv,                           &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d, &
                & vert_interp=create_vert_interp_metadata(                      &
                &   vert_intp_type=vintp_types("P","Z","I"),                    &
                &   vert_intp_method=VINTP_METHOD_LIN ) )

    ! OZONE 
    ! &       field% o3        (nproma,nlev  ,nblks),          &
    cf_desc    = t_cf_var('ozone', 'kg/kg', 'ozone mixing ratio', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0,14,1, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'tro3', field%o3,                         &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d, &
                & vert_interp=create_vert_interp_metadata(                      &
                &   vert_intp_type=vintp_types("P","Z","I"),                    &
                &   vert_intp_method=VINTP_METHOD_LIN ) )

    ! aerosol optical properties
    ! at 533 nm
    cf_desc    = t_cf_var('aer_aod_533','-','aerosol optical depth at 533 nm', &
                & DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0,20,102, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'aer_aod_533', field%aer_aod_533,        &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc,      &
                & ldims=shape3d,                                               &
                & vert_interp=create_vert_interp_metadata(                     &
                &   vert_intp_type=vintp_types("P","Z","I"),                   &
                &   vert_intp_method=VINTP_METHOD_LIN,                         &
                &   l_extrapol=.FALSE. ) )
    cf_desc    = t_cf_var('aer_ssa_533','-',                                   &
                & 'aerosol single scattering albedo at 533 nm', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0,20,103, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'aer_ssa_533', field%aer_ssa_533,        &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc,      &
                & ldims=shape3d,                                               &
                & vert_interp=create_vert_interp_metadata(                     &
                &   vert_intp_type=vintp_types("P","Z","I"),                   &
                &   vert_intp_method=VINTP_METHOD_LIN,                         &
                &   l_extrapol=.FALSE. ) )
    cf_desc    = t_cf_var('aer_asy_533','-',                                   &
                & 'aerosol asymmetry factor at 533 nm', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0,20,104, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'aer_asy_533', field%aer_asy_533,        &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc,      &
                & ldims=shape3d,                                               &
                & vert_interp=create_vert_interp_metadata(                     &
                &   vert_intp_type=vintp_types("P","Z","I"),                   &
                &   vert_intp_method=VINTP_METHOD_LIN ) )
    ! at 2325 nm
    cf_desc    = t_cf_var('aer_aod_2325','-',                                  &
                & 'aerosol optical depth at 2325 nm', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0,20,102, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'aer_aod_2325', field%aer_aod_2325,      &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc,      &
                & ldims=shape3d,                                               &
                & vert_interp=create_vert_interp_metadata(                     &
                &   vert_intp_type=vintp_types("P","Z","I"),                   &
                &   vert_intp_method=VINTP_METHOD_LIN,                         &
                &   l_extrapol=.FALSE. ) )
    cf_desc    = t_cf_var('aer_ssa_2325','-',                                  &
                & 'aerosol single scattering albedo at 2325 nm', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0,20,103, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'aer_ssa_2325', field%aer_ssa_2325,      &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc,      &
                & ldims=shape3d,                                               &
                & vert_interp=create_vert_interp_metadata(                     &
                &   vert_intp_type=vintp_types("P","Z","I"),                   &
                &   vert_intp_method=VINTP_METHOD_LIN,                         &
                &   l_extrapol=.FALSE. ) )
    cf_desc    = t_cf_var('aer_asy_2325','-',                                  &
                & 'aerosol asymmetry factor at 2325 nm', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0,20,104, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'aer_asy_2325', field%aer_asy_2325,      &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc,      &
                & ldims=shape3d,                                               &
                & vert_interp=create_vert_interp_metadata(                     &
                &   vert_intp_type=vintp_types("P","Z","I"),                   &
                &   vert_intp_method=VINTP_METHOD_LIN,                         &
                &   l_extrapol=.FALSE. ) )
    ! at 9731 nm
    cf_desc    = t_cf_var('aer_aod_9731','-',                                  &
                & 'effective aerosol optical depth at 9731 nm', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0,20,102, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'aer_aod_9731', field%aer_aod_9731,      &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc,      &
                & ldims=shape3d,                                               &
                & vert_interp=create_vert_interp_metadata(                     &
                &   vert_intp_type=vintp_types("P","Z","I"),                   &
                &   vert_intp_method=VINTP_METHOD_LIN,                         &
                &   l_extrapol=.FALSE. ) )

    ! &       field% q         (nproma,nlev  ,nblks,ntracer),  &
    CALL add_var( field_list, prefix//'tracer', field%q,                       &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID,                           &
                & t_cf_var('tracer', 'kg kg-1', 'mass mixing ratio of tracers',&
                &          DATATYPE_FLT32),                                    &
                & t_grib2_var(0,20,2, ibits, GRID_REFERENCE, GRID_CELL),       &
                & ldims = (/kproma,klev,kblks,ktracer/),                       &
                & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.         )

    ALLOCATE(field%q_ptr(ktracer))
    DO jtrc = 1,ktracer
      CALL add_ref( field_list, prefix//'tracer',                              &
                  & prefix//TRIM(ctracer(jtrc)), field%q_ptr(jtrc)%p,          &
                  & GRID_UNSTRUCTURED_CELL, ZA_HYBRID,                         &
                  & t_cf_var(TRIM(ctracer(jtrc)), 'kg kg-1',                   &
                  &          'mass mixing ratio of tracer '//TRIM(ctracer(jtrc)), &
                  &          DATATYPE_FLT32),                                  &
                  & t_grib2_var(0,20,2, ibits, GRID_REFERENCE, GRID_CELL),     &
                  & ldims=(/kproma,klev,kblks/),                               &
                  & vert_interp=create_vert_interp_metadata(                   &
                  &             vert_intp_type=vintp_types("P","Z","I"),       & 
                  &             vert_intp_method=VINTP_METHOD_LIN,             &
                  &             l_loglin=.FALSE.,                              &
                  &             l_extrapol=.TRUE., l_pd_limit=.FALSE.,         &
                  &             lower_limit=0._wp )                            )
    END DO                                                                                

    ! &       field% qx        (nproma,nlev  ,nblks),          &
    cf_desc    = t_cf_var('condensated_water', 'kg kg-1', 'cloud water + cloud ice', &
         &                DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0,1,21, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'qx', field%qx,                               &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d, &
                & vert_interp=create_vert_interp_metadata(                          &
                &             vert_intp_type=vintp_types("P","Z","I"),              & 
                &             vert_intp_method=VINTP_METHOD_LIN,                    &
                &             l_loglin=.FALSE.,                                     &
                &             l_extrapol=.TRUE., l_pd_limit=.FALSE.,                &
                &             lower_limit=0._wp  ) )

    ! &       field% omega     (nproma,nlev  ,nblks),          &
    cf_desc    = t_cf_var('vertical_velocity', 'Pa s-1', 'vertical velocity', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0,2,8, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'omega', field%omega,                     &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d, &
                & vert_interp=create_vert_interp_metadata(                       &
                &             vert_intp_type=vintp_types("P","Z","I"),           &
                &             vert_intp_method=VINTP_METHOD_LIN,                 &
                &             l_loglin=.FALSE., l_extrapol=.FALSE.) )

    ! &       field% geom      (nproma,nlev  ,nblks),          &
    cf_desc    = t_cf_var('geopotential', 'm2 s-2', 'geopotential above surface', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 3, 4, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'gpsm', field%geom,                       &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d, &
                & vert_interp=create_vert_interp_metadata(                      &
                &             vert_intp_type=vintp_types("P","Z","I"),          &
                &             vert_intp_method=VINTP_METHOD_LIN,                &
                &             l_extrapol=.TRUE., l_pd_limit=.FALSE.) )

    ! &       field% presm_old (nproma,nlev  ,nblks),          &
    cf_desc    = t_cf_var('pressure', 'Pa', 'pressure at old time step', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 3, 0, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'pam_old', field%presm_old,               &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d, &
                & vert_interp=create_vert_interp_metadata(                      &
                &             vert_intp_type=vintp_types("Z","I"),              &
                &             vert_intp_method=VINTP_METHOD_PRES ) )

    ! &       field% presm_new (nproma,nlev  ,nblks),          &
    cf_desc    = t_cf_var('pressure', 'Pa', 'pressure at new time step', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 3, 0, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'pam_new', field%presm_new,               &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d, &
                & vert_interp=create_vert_interp_metadata(                      &
                &             vert_intp_type=vintp_types("Z","I"),              &
                &             vert_intp_method=VINTP_METHOD_PRES ) )


    !-- Variables defined at layer interfaces --
    ! &       field% geoi      (nproma,nlevp1,nblks),          &
    cf_desc    = t_cf_var('geopotential', 'm2 s-2', 'geopotential above surface', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 3, 4, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'gpsi', field%geoi,                        &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID_HALF, cf_desc, grib2_desc,   &
                & ldims=shape3d_layer_interfaces,                                &
                & vert_interp=create_vert_interp_metadata(                       &
                &   vert_intp_type=vintp_types("P","Z","I"),                     &
                &   vert_intp_method=VINTP_METHOD_LIN_NLEVP1 ) )

    ! &       field% presi_old (nproma,nlevp1,nblks),          &
    cf_desc    = t_cf_var('pressure', 'Pa', 'pressure at old time step', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 3, 0, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'pai_old', field%presi_old,             &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID_HALF, cf_desc, grib2_desc,  &
                & ldims=shape3d_layer_interfaces,                               &
                & vert_interp=create_vert_interp_metadata(                      &
                &   vert_intp_type=vintp_types("Z","I"),                        &
                &   vert_intp_method=VINTP_METHOD_LIN_NLEVP1 ) )

    ! &       field% presi_new (nproma,nlevp1,nblks),          &
    cf_desc    = t_cf_var('pressure', 'Pa', 'pressure at new time step', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 3, 0, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'pai_new', field%presi_new,             &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID_HALF, cf_desc, grib2_desc,  &
                & ldims=shape3d_layer_interfaces,                               &
                & vert_interp=create_vert_interp_metadata(                      &
                &   vert_intp_type=vintp_types("Z","I"),                        &
                &   vert_intp_method=VINTP_METHOD_LIN_NLEVP1 ) )

    !------------------
    ! Radiation
    !------------------
    ! 2D variables

   !ALLOCATE( field% cosmu0    (nproma,       nblks),          &
    cf_desc    = t_cf_var('cosmu0', '', 'cosine of the zenith angle', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(192,214,1, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'cosmu0', field%cosmu0,                   &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    ! &       field% flxdwswtoa(nproma,       nblks),          &
    cf_desc    = t_cf_var('flxdwswtoa', 'W m-2',                                  &
                &         'downward shortwave flux at the top of the atmosphere', &
                &         DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0,4,7, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'rsdt', field%flxdwswtoa,                   &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    ! &       field% visfrcsfc    (nproma,       nblks),          &
    cf_desc    = t_cf_var('visfrcsfc', '', 'visible fraction of sfc net sw', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0,4,255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'visfrcsfc', field%visfrcsfc,             &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    ! &       field% vissfc    (nproma,       nblks),          &
    cf_desc    = t_cf_var('vissfc', 'W m-2', 'net visible flux at surface', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0,4,255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'vissfc', field%vissfc,                   &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    ! &       field% visdffsfc (nproma,       nblks),          &
    cf_desc    = t_cf_var('visdffsfc', 'W m-2', 'net visible diffuse flux at surface', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0,4,255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'visdffsfc', field%visdffsfc,             &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    ! &       field% nirsfc    (nproma,       nblks),          &
    cf_desc    = t_cf_var('nirsfc', 'W m-2', 'net near-IR flux at surface', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0,4,255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'nirsfc', field%nirsfc,                   &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    ! &       field% nirdffsfc (nproma,       nblks),          &
    cf_desc    = t_cf_var('nirdffsfc', 'W m-2', 'net near-IR diffuse flux at surface', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0,4,255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'nirdffsfc', field%nirdffsfc,             &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    ! &       field% parsfcdn    (nproma,       nblks),          &
    cf_desc    = t_cf_var('parsfcdn', 'W m-2',                                            &
                &         'downward photoysnthetically active radiation flux at surface', &
                &         DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0,4,255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'parsfcdn', field%parsfcdn,               &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    ! &       field% pardffsfc (nproma,       nblks),          &
    cf_desc    = t_cf_var('pardffsfc', 'W m-2',                                              &
                &         'net photoysnthetically active radiation diffuse flux at surface', &
                &         DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0,4,255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'pardffsfc', field%pardffsfc,             &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    ! &       field% partrmdnsfc (nproma,       nblks),        &
    cf_desc    = t_cf_var('partrmdnsfc', 'W m-2',                                                   &
                &         'downward photoysnthetically active radiation transmissivity at surface', &
                &          DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0,4,255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'partrmdnsfc', field%partrmdnsfc,         &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    cf_desc    = t_cf_var('rsns', 'W m-2', ' shortwave net flux at surface', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 4, 9, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'rsns', field%swflxsfc,    &
         &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE,            &
         &        cf_desc, grib2_desc,                           &
         &        ldims=shape2d,                                 &
         &        lrestart = .FALSE.,                            &
         &        isteptype=TSTEP_INSTANT )
        
    cf_desc    = t_cf_var('rsnt', 'W m-2', ' shortwave net flux at TOA', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 4, 9, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'rsnt', field%swflxtoa,    &
         &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE,            &
         &        cf_desc, grib2_desc,                           &
         &        ldims=shape2d,                                 &
         &        lrestart = .FALSE.,                            &
         &        isteptype=TSTEP_INSTANT )
        
    cf_desc    = t_cf_var('rlns', 'W m-2', 'longwave net flux at surface', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 5, 5, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'rlns', field%lwflxsfc,    &
         &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE,            &
         &        cf_desc, grib2_desc,                           &
         &        ldims=shape2d,                                 &
         &        lrestart = .FALSE.,                            &
         &        isteptype=TSTEP_INSTANT )

    cf_desc    = t_cf_var('rlus', 'W m-2', 'longwave upward flux at surface', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 5, 5, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'rlus', field%lwupflxsfc,&
         &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE,            &
         &        cf_desc, grib2_desc,                           &
         &        ldims=shape2d,                                 &
         &        lrestart = .FALSE.,                            &
         &        isteptype=TSTEP_INSTANT )

    cf_desc    = t_cf_var('rlus_radt', 'W m-2', 'longwave upward flux at surface at rad time', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 5, 5, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'rlus_radt', field%lwflxupsfc,&
         &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE,            &
         &        cf_desc, grib2_desc,                           &
         &        ldims=shape2d,                                 &
         &        lrestart = .TRUE.,                             &
         &        isteptype=TSTEP_INSTANT )

    cf_desc    = t_cf_var('drlns_dT', 'W m-2 K-1', 'longwave net flux T-derivative at surface', &
         &                DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 5, 5, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'drlns_dT', field%dlwflxsfc_dT, &
         &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                     &
         &        cf_desc, grib2_desc,                                    &
         &        ldims=shape2d,                                          &
         &        lrestart = .FALSE.,                                     &
         &        isteptype=TSTEP_INSTANT )

    cf_desc    = t_cf_var('rlnt', 'W m-2', 'longwave net flux at TOA', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 5, 5, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'rlnt', field%lwflxtoa,    &
         &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE,            &
         &        cf_desc, grib2_desc,                           &
         &        ldims=shape2d,                                 &
         &        lrestart = .FALSE.,                            &
         &        isteptype=TSTEP_INSTANT )

    cf_desc    = t_cf_var('siced', 'm', 'sea ice thickness', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(10,2,1, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'sit', field%siced,      &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    cf_desc    = t_cf_var('alb', '', 'surface albedo from external file', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0,19,1, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'alb', field%alb,      &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    cf_desc    = t_cf_var('surface_height_above_sea_level', 'm',   &
                &         'Mean height above sea level of orography', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0,3,6, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'oromea', field%oromea,      &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    cf_desc    = t_cf_var('standard_deviation_of_height', 'm',     &
                &         'Standard deviation of height above sea level of sub-grid scale orography', &
                &         DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0,3,20, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'orostd', field%orostd,      &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    cf_desc    = t_cf_var('slope_of_terrain', '-',                 &
                &         'Slope of sub-gridscale orography', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0,3,22, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'orosig', field%orosig,      &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    cf_desc    = t_cf_var('anisotropy_factor', '-',                &
                &         'Anisotropy of sub-gridscale orography', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0,3,24, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'orogam', field%orogam,      &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    cf_desc    = t_cf_var('angle_of_principal_axis', 'radians',    &
                &         'Angle of sub-gridscale orography', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0,3,21, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'orothe', field%orothe,      &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    cf_desc    = t_cf_var('height_of_peaks', 'm', 'Height above sea level of peaks', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0,3,6, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'oropic', field%oropic,      &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    cf_desc    = t_cf_var('height_of_valleys', 'm', 'Height above sea level of valleys', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0,3,6, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'oroval', field%oroval,      &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    cf_desc    = t_cf_var('ts_rad', 'K', 'radiative surface temperature', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0,0,17, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'ts_rad', field%tsfc_rad,    &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    cf_desc    = t_cf_var('ts_radt', 'K', 'radiative surface temperature at rad. time', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0,0,17, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'ts_radt', field%tsfc_radt,  &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    cf_desc    = t_cf_var('csat', '', '', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'csat', field%csat,      &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    cf_desc    = t_cf_var('cair', '', '', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'cair', field%cair,      &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    !-------------------------
    ! Sea ice
    !-------------------------

    field%kice = kice ! Number of thickness classes - always 1, as of yet
    shapeice = (/kproma, field%kice, kblks/)

    CALL add_var( field_list, prefix//'ts_icecl', field%Tsurf ,               &
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE,                      &
      &          t_cf_var('ts_icecl', 'C', 'surface temperature',DATATYPE_FLT32),&
      &          t_grib2_var(10,2,8, ibits, GRID_REFERENCE, GRID_CELL),       &
      &          ldims=shapeice)
    CALL add_var( field_list, prefix//'t1_icecl', field%T1 ,                  &
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE,                      &
      &          t_cf_var('t1_icecl','C','Temperature upper layer',DATATYPE_FLT32), &
      &          t_grib2_var(10,2,8, ibits, GRID_REFERENCE, GRID_CELL),       &
      &          ldims=shapeice, lrestart=.TRUE.)
    CALL add_var( field_list, prefix//'t2_icecl', field%T2 ,                  &
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE,                      &
      &          t_cf_var('t2_icecl','C','Temperature lower layer', DATATYPE_FLT32),&
      &          t_grib2_var(10,2,8, ibits, GRID_REFERENCE, GRID_CELL),       &
      &          ldims=shapeice, lrestart=.TRUE.)
    CALL add_var( field_list, prefix//'sit_icecl', field%hi ,                 &
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE,                      &
      &          t_cf_var('sit_icecl', 'm', 'ice thickness', DATATYPE_FLT32), &
      &          t_grib2_var(10,2,1, ibits, GRID_REFERENCE, GRID_CELL),       &
      &          ldims=shapeice, lrestart=.TRUE.)
    CALL add_var( field_list, prefix//'hs_icecl', field%hs ,                  &
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE,                      &
      &          t_cf_var('hs_icecl', 'm', 'snow thickness', DATATYPE_FLT32), &
      &          t_grib2_var(10,2,255, ibits, GRID_REFERENCE, GRID_CELL),&
      &          ldims=shapeice, lrestart=.TRUE.)
    CALL add_var( field_list, prefix//'qtop_icecl', field%Qtop ,              &
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE,                      &
      &          t_cf_var('qtop_icecl', 'W/m^2', 'Energy flux available for surface melting', &
      &                   DATATYPE_FLT32),                                    &
      &          t_grib2_var(10,2,255, ibits, GRID_REFERENCE, GRID_CELL),&
      &          ldims=shapeice, lrestart=.FALSE.)
    CALL add_var( field_list, prefix//'qbot_icecl', field%Qbot ,              &
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE,                      &
      &          t_cf_var('qbot_icecl', 'W/m^2', 'Energy flux at ice-ocean interface', DATATYPE_FLT32),&
      &          t_grib2_var(10,2,255, ibits, GRID_REFERENCE, GRID_CELL),&
      &          ldims=shapeice, lrestart=.FALSE.)


    CALL add_var( field_list, prefix//'sic_icecl', field%conc ,               &
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE,                      &
      &          t_cf_var('sic_icecl', '', 'ice concentration in each ice class', DATATYPE_FLT32),&
      &          t_grib2_var(10,2,0, ibits, GRID_REFERENCE, GRID_CELL),&
      &          ldims=shapeice, lrestart=.TRUE.)

    ! &       field% albvisdir_ice (nproma,field%kice,nblks),          &
    cf_desc    = t_cf_var('albvisdir_icecl', '', 'ice albedo VIS direct', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(192,128,15, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'albvisdir_icecl', field%albvisdir_ice,             &
                & GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, cf_desc, grib2_desc, ldims=shapeice )

    ! &       field% albvisdif_ice (nproma,field%kice,nblks),          &
    cf_desc    = t_cf_var('albvisdif_icecl', '', 'ice albedo VIS diffuse', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0,19,222, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'albvisdif_icecl', field%albvisdif_ice,             &
                & GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, cf_desc, grib2_desc, ldims=shapeice )

    ! &       field% albnirdir_ice (nproma,field%kice,nblks),          &
    cf_desc    = t_cf_var('albnirdir_icecl', '', 'ice albedo NIR direct', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(192,128,17, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'albnirdir_icecl', field%albnirdir_ice,             &
                & GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, cf_desc, grib2_desc, ldims=shapeice )

    ! &       field% albnirdif_ice (nproma,field%kice,nblks),          &
    cf_desc    = t_cf_var('albnirdif_icecl', '', 'ice albedo NIR diffuse', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 19, 223, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'albnirdif_icecl', field%albnirdif_ice,             &
                & GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, cf_desc, grib2_desc, ldims=shapeice )


    !---- 3D variables defined at layer interfaces ----

    ! &       field% lwflxclr  (nproma,nlevp1,nblks),          &
    cf_desc    = t_cf_var('lwflxclr', 'W/m2', 'net longwave flux, clear-sky, positive downward', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0,5,6, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'lwflxclr', field%lwflxclr,               &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID_HALF, cf_desc, grib2_desc,  &
                & ldims=shape3d_layer_interfaces,                               &
                & vert_interp=create_vert_interp_metadata(                      &
                &   vert_intp_type=vintp_types("P","Z","I"),                    &
                &   vert_intp_method=VINTP_METHOD_LIN_NLEVP1 ) )

    ! &       field% lwflxall  (nproma,nlevp1,nblks),          &
    cf_desc    = t_cf_var('lwflxall', 'W/m2', 'net longwave flux, all-sky, positive downward', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0,5,5, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'lwflxall', field%lwflxall,               &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID_HALF, cf_desc, grib2_desc,  &
                & ldims=shape3d_layer_interfaces,                               &
                & vert_interp=create_vert_interp_metadata(                      &
                &   vert_intp_type=vintp_types("P","Z","I"),                    &
                &   vert_intp_method=VINTP_METHOD_LIN_NLEVP1 ) )

    ! &       field% swtrmclr  (nproma,nlevp1,nblks),          &
    cf_desc    = t_cf_var('swtrmclr', '', 'shortwave transmissivity, clear-sky', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0,4,255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'swtrmclr', field%swtrmclr,               &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID_HALF, cf_desc, grib2_desc,  &
                & ldims=shape3d_layer_interfaces,                               &
                & vert_interp=create_vert_interp_metadata(                      &
                &   vert_intp_type=vintp_types("P","Z","I"),                    &
                &   vert_intp_method=VINTP_METHOD_LIN_NLEVP1 ) )

    ! &       field% swtrmall  (nproma,nlevp1,nblks),          &
    cf_desc    = t_cf_var('swtrmall', '', 'shortwave transmissivity, all-sky', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0,4,255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'swtrmall', field%swtrmall,               &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID_HALF, cf_desc, grib2_desc,  &
                & ldims=shape3d_layer_interfaces,                               &
                & vert_interp=create_vert_interp_metadata(                      &
                &   vert_intp_type=vintp_types("P","Z","I"),                    &
                &   vert_intp_method=VINTP_METHOD_LIN_NLEVP1 ) )


    !-------------------------
    ! Cloud and precipitation
    !-------------------------
    cf_desc    = t_cf_var('cl', 'm2 m-2', 'cloud area fraction', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0,6,22, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'cl', field%aclc,                                  &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d, &
                & vert_interp=create_vert_interp_metadata(                               &
                &             vert_intp_type=vintp_types("P","Z","I"),                   &
                &             vert_intp_method=VINTP_METHOD_LIN,                         &
                &             l_loglin=.FALSE.,                                          &
                &             l_extrapol=.FALSE., l_pd_limit=.FALSE.,                    &
                &             lower_limit=0._wp ) )

    cf_desc    = t_cf_var('clt', 'm2 m-2', &
               & 'total cloud cover', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0,6,1, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'clt', field%aclcov,       &
         &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE,            &
         &        cf_desc, grib2_desc,                           &
         &        ldims=shape2d,                                 &
         &        lrestart = .FALSE.,                            &
         &        isteptype=TSTEP_INSTANT )

    ! &       field% acdnc  (nproma,nlev  ,nblks), &
    cf_desc    = t_cf_var('acdnc', 'm-3', 'cloud droplet number concentration', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0,6,28, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'acdnc', field%acdnc,                     &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d, &
                & vert_interp=create_vert_interp_metadata(                               &
                &             vert_intp_type=vintp_types("P","Z","I"),                   &
                &             vert_intp_method=VINTP_METHOD_LIN,                         &
                &             l_loglin=.FALSE.,                                          &
                &             l_extrapol=.TRUE., l_pd_limit=.FALSE.,                     &
                &             lower_limit=0._wp ) )

    ! &       field% xvar   (nproma,nlev  ,nblks), &
    cf_desc    = t_cf_var('variance_of_total_water', '', 'subgrid variance of total water', &
         &                DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0,6,255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'xvar', field%xvar,                       &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d, &
                & vert_interp=create_vert_interp_metadata(                               &
                &             vert_intp_type=vintp_types("P","Z","I"),                   &
                &             vert_intp_method=VINTP_METHOD_LIN,                         &
                &             l_loglin=.FALSE.,                                          &
                &             l_extrapol=.TRUE., l_pd_limit=.FALSE.,                     &
                &             lower_limit=0._wp ) )

    ! &       field% xskew  (nproma,nlev  ,nblks), &
    cf_desc    = t_cf_var('skewness_of_total_water', '', 'skewness of total water', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0,6,255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'xskew', field%xskew,                     &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d, &
                & vert_interp=create_vert_interp_metadata(                               &
                &             vert_intp_type=vintp_types("P","Z","I"),                   &
                &             vert_intp_method=VINTP_METHOD_LIN,                         &
                &             l_loglin=.FALSE.,                                          &
                &             l_extrapol=.TRUE., l_pd_limit=.FALSE. ) )

    ! &       field% relhum (nproma,nlev  ,nblks), &
    cf_desc    = t_cf_var('relative_humidity', '', 'relative humidity', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 1, 1, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'hur', field%relhum,                               &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d, &
                & vert_interp=create_vert_interp_metadata(                               &
                &             vert_intp_type=vintp_types("P","Z","I"),                   &
                &             vert_intp_method=VINTP_METHOD_LIN,                         &
                &             l_loglin=.FALSE.,                                          &
                &             l_extrapol=.FALSE., l_pd_limit=.TRUE.,                     &
                &             lower_limit=0._wp ) )

    cf_desc    = t_cf_var('prlr', 'kg m-2 s-1',    &
               & 'large-scale precipitation flux (water)', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0,1,77, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'prlr', field%rsfl,        &
         &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE,            &
         &        cf_desc, grib2_desc,                           &
         &        ldims=shape2d,                                 &
         &        lrestart = .TRUE.,                             &
         &        isteptype=TSTEP_INSTANT )

    cf_desc    = t_cf_var('prcr', 'kg m-2 s-1',    &
               & 'convective precipitation flux (water)', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0,1,76, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'prcr', field%rsfc,        &
         &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE,            &
         &        cf_desc, grib2_desc,                           &
         &        ldims=shape2d,                                 &
         &        lrestart = .TRUE.,                             &
         &        isteptype=TSTEP_INSTANT )

    cf_desc    = t_cf_var('prls', 'kg m-2 s-1',    &
               & 'large-scale precipitation flux (snow)', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0,1,59, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'prls', field%ssfl,        &
         &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE,            &
         &        cf_desc, grib2_desc,                           &
         &        ldims=shape2d,                                 &
         &        lrestart = .TRUE.,                             &
         &        isteptype=TSTEP_INSTANT )

    cf_desc    = t_cf_var('prcs', 'kg m-2 s-1',    &
               & 'convective precipitation flux (snow)', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0,1,58, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'prcs', field%ssfc,        &
         &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE,            &
         &        cf_desc, grib2_desc,                           &
         &        ldims=shape2d,                                 &
         &        lrestart = .TRUE.,                             &
         &        isteptype=TSTEP_INSTANT )

    cf_desc    = t_cf_var('TOTPREC', 'kg m-2 s-1',               &
         &                'total precipitation flux',            &
         &                DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 1, 52, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'pr', field%totprec,       &
         &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE,            &
         &        cf_desc, grib2_desc,                           &
         &        ldims=shape2d,                                 &
         &        lrestart = .FALSE.,                            &
         &        isteptype=TSTEP_INSTANT )

    cf_desc    = t_cf_var('TOTPREC_AVG', 'kg m-2 s-1',              &
         &                'time averaged total precipitation flux', &
         &       DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 1, 52, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'totprec_avg', field%totprec_avg, &
         &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                   &
         &        cf_desc, grib2_desc,                                  &
         &        ldims=shape2d,                                        &
         &        lrestart = .TRUE.,                                    &
         &        isteptype=TSTEP_AVG )

    cf_desc    = t_cf_var('total_vapour', 'kg m-2', 'vertically integrated water vapour', &
         &                DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0,1,64, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'prw', field%qvi,          &
         &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE,            &
         &        cf_desc, grib2_desc,                           &
         &        ldims=shape2d,                                 &
         &        lrestart = .FALSE.,                            &
         &        isteptype=TSTEP_INSTANT )

    cf_desc    = t_cf_var('total_cloud_water', 'kg m-2',&
               & 'vertically integrated cloud water', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0,1,69, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'cllvi', field%xlvi,       &
         &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE,            &
         &        cf_desc, grib2_desc,                           &
         &        ldims=shape2d,                                 &
         &        lrestart = .FALSE.,                            &
         &        isteptype=TSTEP_INSTANT )

    cf_desc    = t_cf_var('total_cloud_ice', 'kg m-2',&
               & 'vertically integrated cloud ice', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0,1,70, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'clivi', field%xivi,       &
         &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE,            &
         &        cf_desc, grib2_desc,                           &
         &        ldims=shape2d,                                 &
         &        lrestart = .FALSE.,                            &
         &        isteptype=TSTEP_INSTANT )

    ! &       field% rintop (nproma,       nblks), &
    cf_desc    = t_cf_var('rintop', '', '', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0,6,255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'rintop', field%rintop,                   &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    ! &       field% rtype  (nproma,       nblks), &
    cf_desc    = t_cf_var('convection_type', '', 'convection_type (0...3)', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0,6,255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'rtype', field%rtype,                     &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    ! &       field% topmax (nproma,       nblks), &
    cf_desc    = t_cf_var('topmax', 'Pa', 'maximum height of convective cloud tops', &
         &                DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0,6,255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'topmax', field%topmax,                   &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    ! &       field% tke    (nproma,nlev  ,nblks), &
    cf_desc    = t_cf_var('turbulent_kinetic_energy', 'J kg-1', 'turbulent kinetic energy', &
         &                DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0,19,11, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'tke', field%tke,                         &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d, &
                & vert_interp=create_vert_interp_metadata(                               &
                &             vert_intp_type=vintp_types("P","Z","I"),                   &
                &             vert_intp_method=VINTP_METHOD_LIN,                         &
                &             l_loglin=.FALSE.,                                          &
                &             l_extrapol=.TRUE., l_pd_limit=.FALSE.,                     &
                &             lower_limit=0._wp ) )

    ! &       field% thvsig (nproma,       nblks), &
    cf_desc    = t_cf_var('thvsig', 'K', '', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0,19,255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'thvsig', field%thvsig,                   &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    !---------------------------
    ! Variables for energy diagnostic of echam6 physics
    !---------------------------

       cf_desc    = t_cf_var('sh_vdiff','J m-2 s-1', '', DATATYPE_FLT32)
       grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
       CALL add_var( field_list, prefix//'sh_vdiff', field%sh_vdiff,             &
                   & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )
       cf_desc    = t_cf_var('ch_concloud','J m-2 s-1', '', DATATYPE_FLT32)
       grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
       CALL add_var( field_list, prefix//'ch_concloud', field%ch_concloud,       &
                   & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )
       cf_desc    = t_cf_var('con_dtrl','?', '', DATATYPE_FLT32)
       grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
       CALL add_var( field_list, prefix//'con_dtrl', field%con_dtrl,       &
                   & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )
       cf_desc    = t_cf_var('con_dtri','?', '', DATATYPE_FLT32)
       grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
       CALL add_var( field_list, prefix//'con_dtri', field%con_dtri,       &
                   & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )
       cf_desc    = t_cf_var('con_iteqv','?', '', DATATYPE_FLT32)
       grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
       CALL add_var( field_list, prefix//'con_iteqv', field%con_iteqv,       &
                   & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )
       cf_desc    = t_cf_var('qv_vdiff','kg/m^2/s', '', DATATYPE_FLT32)
       grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
       CALL add_var( field_list, prefix//'qv_vdiff', field%qv_vdiff,             &
                   & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )
       cf_desc    = t_cf_var('cld_dtrl','?', '', DATATYPE_FLT32)
       grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
       CALL add_var( field_list, 'cld_dtrl', field%cld_dtrl,       &
                   & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )
       cf_desc    = t_cf_var('cld_dtri','?', '', DATATYPE_FLT32)
       grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
       CALL add_var( field_list, 'cld_dtri', field%cld_dtri,       &
                   & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )
       cf_desc    = t_cf_var('cld_iteq','?', '', DATATYPE_FLT32)
       grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
       CALL add_var( field_list, 'cld_iteq', field%cld_iteq,       &
                   & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    !---------------------------
    ! Orographic wave drag diagnostics
    !---------------------------
    CALL add_var( field_list, prefix//'tauu_sso', field%u_stress_sso,             &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                             &
                & t_cf_var('u_stress_sso', 'N m-2',                               &
                &          'zonal stress from subgrid scale orographic drag',     &
                &          DATATYPE_FLT32),                                       &
                & t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),   &
                & ldims=shape2d,                                                  &
                & lrestart = .FALSE.,                                             &
                & isteptype=TSTEP_INSTANT                                         )

    CALL add_var( field_list, prefix//'tauv_sso', field%v_stress_sso,             &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                             &
                & t_cf_var('v_stress_sso', 'N m-2',                               &
                &          'meridional stress from subgrid scale orographic drag',&
                &          DATATYPE_FLT32),                                       &
                & t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),   &
                & ldims=shape2d,                                                  &
                & lrestart = .FALSE.,                                             &
                & isteptype=TSTEP_INSTANT                                         )

    CALL add_var( field_list, prefix//'diss_sso', field%dissipation_sso,          &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                             &
                & t_cf_var('dissipation_sso', '',                                 &
                &          'dissipation of orographic waves',                     &
                &          DATATYPE_FLT32),                                       &
                & t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),   &
                & ldims=shape2d,                                                  &
                & lrestart = .FALSE.,                                             &
                & isteptype=TSTEP_INSTANT                                         )

    !--------------------
    ! Turbulence
    !--------------------

     ! shape2d  = (/kproma,            kblks/)
     ! shape3d  = (/kproma, klev,      kblks/)
     !shapesfc = (/kproma, ksfc_type, kblks/)
     ! shapesfc = (/kproma, kblks, ksfc_type/)

     !ALLOCATE( field% ri     (nproma,nlev,nblks), &
      cf_desc    = t_cf_var('richardson_number', ' ', 'moist Richardson number', DATATYPE_FLT32)
      grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( field_list, prefix//'ri', field%ri,                         &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

      ! &       field% mixlen (nproma,nlev,nblks), &
      cf_desc    = t_cf_var('mixing_length', 'm', 'mixing_length', DATATYPE_FLT32)
      grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( field_list, prefix//'mixlen', field%mixlen,                 &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

      ! &       field% thvvar (nproma,nlev,nblks), &
      cf_desc    = t_cf_var('thvvar', 'K2',                           &
                 & 'subgrid variance of virtual potential temperature', DATATYPE_FLT32)
      grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( field_list, prefix//'thvvar', field%thvvar,                 &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

      ! &       field% tkem0  (nproma,nlev,nblks), &
      cf_desc    = t_cf_var('tke', 'm2 s-2', 'TKE at step t', DATATYPE_FLT32)
      grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( field_list, prefix//'tkem0', field%tkem0,                   &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

      ! &       field% tkem1  (nproma,nlev,nblks), &
      cf_desc    = t_cf_var('tke', 'm2 s-2', 'TKE at step t-dt', DATATYPE_FLT32)
      grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( field_list, prefix//'tkem1', field%tkem1,                   &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

     !---------
     !ALLOCATE( field% cfm    (nproma,nlev,     nblks), &
      cf_desc    = t_cf_var('turb_exchng_coeff_momentum', '', '', DATATYPE_FLT32)
      grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( field_list, prefix//'cfm', field%cfm,                       &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

      ! &       field% cfm_tile(nproma,nblks,nsfc_type), &
      CALL add_var( field_list, prefix//'cfm_tile', field%cfm_tile,              &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                          &
                  & t_cf_var('turb_exchng_coeff_momentum', '', '', DATATYPE_FLT32), &
                  & t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),&
                  & ldims=shapesfc,                                              &
                  & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.         )

      ALLOCATE(field%cfm_tile_ptr(ksfc_type))
      DO jsfc = 1,ksfc_type
        CALL add_ref( field_list, prefix//'cfm_tile',                              &
                    & prefix//'cfm_'//csfc(jsfc), field%cfm_tile_ptr(jsfc)%p,      &
                    & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                          &
                    & t_cf_var('turb_exchng_coeff_momentum_'//csfc(jsfc), '', '', DATATYPE_FLT32),&
                    & t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),&
                    & ldims=shape2d )
      END DO

      !---------
      ! &       field% cfh    (nproma,nlev,     nblks), &
      cf_desc    = t_cf_var('turb_exchng_coeff_heat', '', '', DATATYPE_FLT32)
      grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( field_list, prefix//'cfh', field%cfh,                       &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

      ! &       field% cfh_tile(nproma,nblks,nsfc_type), &
      CALL add_var( field_list, prefix//'cfh_tile', field%cfh_tile,              &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                          &
                  & t_cf_var('turb_exchng_coeff_heat', '', '', DATATYPE_FLT32),  &
                  & t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),&
                  & ldims=shapesfc,                                              &
                  & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.         )

      ALLOCATE(field%cfh_tile_ptr(ksfc_type))
      DO jsfc = 1,ksfc_type
        CALL add_ref( field_list, prefix//'cfh_tile',                              &
                    & prefix//'cfh_'//csfc(jsfc), field%cfh_tile_ptr(jsfc)%p,      &
                    & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                          &
                    & t_cf_var('turb_exchng_coeff_heat_'//csfc(jsfc), '', '',      &
                    &          DATATYPE_FLT32),                                    &
                    & t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),&
                    & ldims=shape2d )
      END DO

      !---------
      ! &       field% cfv    (nproma,nlev,     nblks), &
      cf_desc    = t_cf_var('turb_exchng_coeff_water_var', '', '', DATATYPE_FLT32)
      grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( field_list, prefix//'cfv', field%cfv,                       &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

      ! &       field% cftke  (nproma,nlev,     nblks), &
      cf_desc    = t_cf_var('turb_exchng_coeff_tke', '', '', DATATYPE_FLT32)
      grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( field_list, prefix//'cftke', field%cftke,                   &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

      ! &       field% cfthv  (nproma,nlev,     nblks)  )
      cf_desc    = t_cf_var('turb_exchng_coeff_thv', '', '', DATATYPE_FLT32)
      grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( field_list, prefix//'cfthv', field%cfthv,                   &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

     !ALLOCATE( field% coriol (nproma,nblks),                &
      cf_desc    = t_cf_var('Coriolis_param', 's-1', 'Coriolis parameter', DATATYPE_FLT32)
      grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( field_list, prefix//'coriol', field%coriol,                 &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

      ! &       field% ghpbl  (nproma,nblks),                &
      cf_desc    = t_cf_var('geopot_pbl_top', '', 'geopotential of PBL top', DATATYPE_FLT32)
      grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( field_list, prefix//'ghpbl', field%ghpbl,                   &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

      !-----------------------------------
      ! &       field% z0m(nproma,nblks), &
      cf_desc    = t_cf_var('z0m', '', 'aerodynamic roughness length, grid box mean', &
           &                DATATYPE_FLT32)
      grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( field_list, prefix//'z0m', field%z0m,                       &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

      ! &       field% z0m_tile(nproma,nblks,nsfc_type), &
      CALL add_var( field_list, prefix//'z0m_tile', field%z0m_tile,              &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                          &
                  & t_cf_var('z0m_tile', '', 'aerodynamic roughness length', DATATYPE_FLT32),&
                  & t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),&
                  & ldims=shapesfc,                                              &
                  & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.         )

      ALLOCATE(field%z0m_tile_ptr(ksfc_type))
      DO jsfc = 1,ksfc_type
        CALL add_ref( field_list, prefix//'z0m_tile',                              &
                    & prefix//'z0m_'//csfc(jsfc), field%z0m_tile_ptr(jsfc)%p,      &
                    & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                          &
                    & t_cf_var('z0m_'//csfc(jsfc), '','', DATATYPE_FLT32),         &
                    & t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),&
                    & ldims=shape2d )
      END DO

      ! &        field% z0h_lnd(nproma, nblks), &
      cf_desc    = t_cf_var('z0h_lnd', '', 'roughness length heat, land', &
        &                DATATYPE_FLT32)
      grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( field_list, prefix//'z0h_lnd', field%z0h_lnd,                  &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

      !-----------------------------------

      ! &       field% ustar  (nproma,nblks),                &
      cf_desc    = t_cf_var('friction_velocity', 'm s-1', 'friction velocity', DATATYPE_FLT32)
      grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( field_list, prefix//'ustar', field%ustar,                   &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

      ! &       field% kedisp (nproma,nblks),                &
      cf_desc    = t_cf_var('KE dissipation rate', '', '', DATATYPE_FLT32)
      grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( field_list, prefix//'kedisp', field%kedisp,                 &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

      ! &       field% ocu    (nproma,nblks),                &
      cf_desc    = t_cf_var('ocean_sfc_u', 'm/s', 'u-component of ocean current/ice', DATATYPE_FLT32)
      grib2_desc = t_grib2_var(10,1,2, iextbits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( field_list, prefix//'ocu', field%ocu, &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE,   &
        &           cf_desc, grib2_desc, ldims=shape2d,   &
        &           lrestart=.TRUE. )

      ! &       field% ocv    (nproma,nblks),                &
      cf_desc    = t_cf_var('ocean_sfc_v', 'm/s', 'v-component of ocean current/ice', DATATYPE_FLT32)
      grib2_desc = t_grib2_var(10,1,3, iextbits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( field_list, prefix//'ocv', field%ocv, &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE,   &
        &           cf_desc, grib2_desc, ldims=shape2d,   &
        &           lrestart=.TRUE. )

    !-----------------------
    ! Surface
    !-----------------------
   !ALLOCATE( field% lfland (nproma, nblks),                 &
    cf_desc    = t_cf_var('lfland', '', '', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'lfland', field%lfland,        &
              & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                  &
              & cf_desc, grib2_desc, ldims=shape2d, lrestart=.FALSE. )

    ! &       field% lfglac (nproma, nblks),                 &
    cf_desc    = t_cf_var('lfglac', '', '', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'lfglac', field%lfglac,         &
              & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                   &
              & cf_desc, grib2_desc, ldims=shape2d, lrestart=.FALSE. )

  ! ALLOCATE( field% lfland (kproma, kblks), &
  !         & field% lfglac (kproma, kblks)  ) 

    ! &       field% lsmask (nproma, nblks),                 &
    cf_desc    = t_cf_var('land_cover', '', 'land cover', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(2, 0, 0, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'land', field%lsmask,                   &
              & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    ! &       field% glac   (nproma, nblks),                 &
    cf_desc    = t_cf_var('glacier_cover', '', 'fraction of land covered by glaciers', &
         &                DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'glac', field%glac,                     &
              & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    ! &       field% seaice (nproma, nblks),                 &
    cf_desc    = t_cf_var('sea_ice_cover', '', 'fraction of ocean covered by sea ice', &
         &                DATATYPE_FLT32)
    grib2_desc = t_grib2_var(10,2,0, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'sic', field%seaice,    &
      &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE,         &
      &           cf_desc, grib2_desc, ldims=shape2d,         &
      &           lrestart=.TRUE. )

    ! &       field% alake (nproma, nblks),                 &
    cf_desc    = t_cf_var('alake', '', 'fraction of lakes', &
         &                DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'alake', field%alake,                 &
              & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    ! &       field% icefrc (nproma, nblks),                 &
    cf_desc    = t_cf_var('ice_cover', '', 'ice cover given as fraction of grid box', & 
         &                DATATYPE_FLT32)
    grib2_desc = t_grib2_var(10,2,0, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'icefrc', field%icefrc,                 &
              & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    !-----------------------------------
    ! &       field% tsfc(nproma,nblks), &
    cf_desc    = t_cf_var('surface_temperature', 'K', 'surface temperature', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0,0,0, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'ts', field%tsfc,                       &
              & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    ! &       field% tsfc_tile(nproma,nblks,nsfc_type), &
    CALL add_var( field_list, prefix//'ts_tile', field%tsfc_tile,              &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                          &
                & t_cf_var('_tile', 'K', 'surface temperature on tiles', DATATYPE_FLT32), &
                & t_grib2_var(0,0,0, ibits, GRID_REFERENCE, GRID_CELL),        &
                & ldims=shapesfc,                                              &
                & lcontainer=.TRUE., lrestart=.FALSE.         )

    ALLOCATE(field%tsfc_tile_ptr(ksfc_type))
    DO jsfc = 1,ksfc_type
      CALL add_ref( field_list, prefix//'ts_tile',                                   &
                  & prefix//'ts_'//csfc(jsfc), field%tsfc_tile_ptr(jsfc)%p,          &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                              &
                  & t_cf_var('ts_'//csfc(jsfc), 'K',                                 &
                  &          'surface temperature on '//csfc(jsfc), DATATYPE_FLT32), &
                  & t_grib2_var(0,0,0, ibits, GRID_REFERENCE, GRID_CELL),            &
                  & ldims=shape2d, lrestart=.TRUE. )
    END DO
    !-----------------------------------

    ! &       field% qs_sfc_tile (nproma,nblks,nsfc_type), &
    CALL add_var( field_list, prefix//'qs_sfc_tile', field%qs_sfc_tile,        &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                          &
                & t_cf_var('qs_sfc_tile', '', '', DATATYPE_FLT32),             &
                & t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),&
                & ldims=shapesfc,                                              &
                & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.         )

    ALLOCATE(field%qs_sfc_tile_ptr(ksfc_type))
    DO jsfc = 1,ksfc_type
      CALL add_ref( field_list, prefix//'qs_sfc_tile',                           &
                  & prefix//'qs_sfc_'//csfc(jsfc), field%qs_sfc_tile_ptr(jsfc)%p,&
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                          &
                  & t_cf_var('qs_sfc_'//csfc(jsfc), '', '', DATATYPE_FLT32),     &
                  & t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),&
                  & ldims=shape2d )
    END DO

    !-----------------------------------
    ! &       field% albedo (nproma,nblks),          &
    cf_desc    = t_cf_var('albedo', '', 'surface albedo', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'albedo', field%albedo,                            &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    ! &       field% albvisdir (nproma,nblks),          &
    cf_desc    = t_cf_var('albvisdir', '', 'albedo VIS direct', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'albvisdir', field%albvisdir,                       &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d  )

    ! &       field% albvisdif (nproma,nblks),          &
    cf_desc    = t_cf_var('albvisdif', '', 'albedo VIS diffuse', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'albvisdif', field%albvisdif,                       &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d  )

    ! &       field% albnirdir (nproma,nblks),          &
    cf_desc    = t_cf_var('albnirdir', '', 'albedo NIR direct', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'albnirdir', field%albnirdir,                       &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d  )

    ! &       field% albnirdif (nproma,nblks),          &
    cf_desc    = t_cf_var('albnirdif', '', 'albedo NIR diffuse', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'albnirdif', field%albnirdif,                       &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d  )

    ! &       field% albvisdir_tile (nproma,nblks,nsfc_type),          &
    cf_desc    = t_cf_var('albvisdir_tile', '', 'albedo VIS direct', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'albvisdir_tile', field%albvisdir_tile,             &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shapesfc,&
                & lcontainer=.TRUE., lrestart=.FALSE.         )

    ! &       field% albvisdif_tile (nproma,nblks,nsfc_type),          &
    cf_desc    = t_cf_var('albvisdif_tile', '', 'albedo VIS diffuse', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'albvisdif_tile', field%albvisdif_tile,             &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shapesfc,&
                & lcontainer=.TRUE., lrestart=.FALSE.         )

    ! &       field% albnirdir_tile (nproma,nblks,nsfc_type),          &
    cf_desc    = t_cf_var('albnirdir_tile', '', 'albedo NIR direct', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'albnirdir_tile', field%albnirdir_tile,             &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shapesfc,&
                & lcontainer=.TRUE., lrestart=.FALSE.         )

    ! &       field% albnirdif_tile (nproma,nblks,nsfc_type),          &
    cf_desc    = t_cf_var('albnirdif_tile', '', 'albedo NIR diffuse', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'albnirdif_tile', field%albnirdif_tile,             &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shapesfc,&
                & lcontainer=.TRUE., lrestart=.FALSE.         )

    ! &       field% albedo_tile (nproma,nblks,nsfc_type),          &
    cf_desc    = t_cf_var('albedo_tile', '', 'albedo', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'albedo_tile', field%albedo_tile,                   &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shapesfc,&
                & lcontainer=.TRUE., lrestart=.FALSE.         )

    ALLOCATE(field%albvisdir_tile_ptr(ksfc_type), field%albvisdif_tile_ptr(ksfc_type), &
             field%albnirdir_tile_ptr(ksfc_type), field%albnirdif_tile_ptr(ksfc_type), &
             field%albedo_tile_ptr(ksfc_type)                                          )

    DO jsfc = 1,ksfc_type
      CALL add_ref( field_list, prefix//'albvisdir_tile',                              &
                  & prefix//'albvisdir_'//csfc(jsfc), field%albvisdir_tile_ptr(jsfc)%p,&
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                &
                  & t_cf_var('albvisdir_'//csfc(jsfc), '', '', DATATYPE_FLT32),        &
                  & t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),      &
                  & ldims=shape2d )
      CALL add_ref( field_list, prefix//'albvisdif_tile',                              &
                  & prefix//'albvisdif_'//csfc(jsfc), field%albvisdif_tile_ptr(jsfc)%p,&
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                &
                  & t_cf_var('albvisdif_'//csfc(jsfc), '', '', DATATYPE_FLT32),        &
                  & t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),      &
                  & ldims=shape2d )
      CALL add_ref( field_list, prefix//'albnirdir_tile',                              &
                  & prefix//'albnirdir_'//csfc(jsfc), field%albnirdir_tile_ptr(jsfc)%p,&
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                &
                  & t_cf_var('albnirdir_'//csfc(jsfc), '', '', DATATYPE_FLT32),        &
                  & t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),      &
                  & ldims=shape2d )
      CALL add_ref( field_list, prefix//'albnirdif_tile',                              &
                  & prefix//'albnirdif_'//csfc(jsfc), field%albnirdif_tile_ptr(jsfc)%p,&
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                &
                  & t_cf_var('albnirdif_'//csfc(jsfc), '', '', DATATYPE_FLT32),        &
                  & t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),      &
                  & ldims=shape2d )
      CALL add_ref( field_list, prefix//'albedo_tile',                                 &
                  & prefix//'albedo_'//csfc(jsfc), field%albedo_tile_ptr(jsfc)%p,      &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                &
                  & t_cf_var('albedo_'//csfc(jsfc), '', '', DATATYPE_FLT32),           &
                  & t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),      &
                  & ldims=shape2d )
    END DO

    !---------------------------
    ! Surface fluxes
    !---------------------------
    ! gridbox mean

    CALL add_var( field_list, prefix//'evspsbl', field%evap,              &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                     &
                & t_cf_var('evap', 'kg m-2 s-1', 'evaporation',           &
                & DATATYPE_FLT32),                                        &
                & t_grib2_var(0,1,6,iextbits, GRID_REFERENCE, GRID_CELL), &
                & ldims=shape2d,                                          &
                & lrestart = .FALSE.,                                     &
                & isteptype=TSTEP_INSTANT                                 )

    CALL add_var( field_list, prefix//'hfls', field%lhflx,                &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                     &
                & t_cf_var('lhflx', 'W m-2 ', 'latent heat flux',         &
                & DATATYPE_FLT32),                                        &
                & t_grib2_var(0,0,10, ibits, GRID_REFERENCE, GRID_CELL),  &
                & ldims=shape2d,                                          &
                & lrestart = .FALSE.,                                     &
                & isteptype=TSTEP_INSTANT                                 )

    CALL add_var( field_list, prefix//'hfss', field%shflx,                &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                     &
                & t_cf_var('shflx', 'W m-2 ', 'sensible heat flux',       &
                & DATATYPE_FLT32),                                        &
                & t_grib2_var(0,0,11, ibits, GRID_REFERENCE, GRID_CELL),  &
                & ldims=shape2d,                                          &
                & lrestart = .FALSE.,                                     &
                & isteptype=TSTEP_INSTANT                                 )

    !---------------------------------
    ! values on tiles

    CALL add_var( field_list, prefix//'rsns_tile',field%swflxsfc_tile,    &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                     &
                & t_cf_var('rsns_tile', 'W m-2',                          &
                &          'shortwave net flux at surface on tiles',      &
                &          DATATYPE_FLT32),                               &
                & t_grib2_var(0,4,9, ibits, GRID_REFERENCE, GRID_CELL),   &
                & ldims=shapesfc,                                         &
                & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.    )

    CALL add_var( field_list, prefix//'rlns_tile',field%lwflxsfc_tile,    &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                     &
                & t_cf_var('rlns_tile', 'W m-2',                          &
                &          'longwave net flux at surface on tiles',       &
                &          DATATYPE_FLT32),                               &
                & t_grib2_var(0,5,5, ibits, GRID_REFERENCE, GRID_CELL),   &
                & ldims=shapesfc,                                         &
                & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.    )

    CALL add_var( field_list, prefix//'evspsbl_tile', field%evap_tile,    &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                     &
                & t_cf_var('evspsbl_tile', 'kg m-2 s-1',                  &
                &          'evaporation on tiles', DATATYPE_FLT32),       &
                & t_grib2_var(0,1,6, ibits, GRID_REFERENCE, GRID_CELL),   &
                & ldims=shapesfc,                                         &
                & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.    )

    CALL add_var( field_list, prefix//'hfls_tile', field%lhflx_tile,      &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                     &
                & t_cf_var('hfls_tile', 'W m-2',                          &
                &          'latent heat flux on tiles', DATATYPE_FLT32),  &
                & t_grib2_var(0,0,10, ibits, GRID_REFERENCE, GRID_CELL),  &
                & ldims=shapesfc,                                         &
                & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.    )

    CALL add_var( field_list, prefix//'hfss_tile', field%shflx_tile,      &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                     &
                & t_cf_var('hfss_tile', 'W m-2',                          &
                &          'sensible heat flux on tiles', DATATYPE_FLT32),&
                & t_grib2_var(0,0,11, ibits, GRID_REFERENCE, GRID_CELL),  &
                & ldims=shapesfc,                                         &
                & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.    )

    CALL add_var( field_list,prefix//'dhfss_dT_tile',field%dshflx_dT_tile,&
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                     &
                & t_cf_var('dhfss_dT_tile', 'W m-2 K-1',                  &
                &          'temp tend of SHF on tiles', DATATYPE_FLT32),  &
                & t_grib2_var(0,0,11, ibits, GRID_REFERENCE, GRID_CELL),  &
                & ldims=shapesfc,                                         &
                & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.    )


    ALLOCATE(field%swflxsfc_tile_ptr(ksfc_type))
    ALLOCATE(field%lwflxsfc_tile_ptr(ksfc_type))
    ALLOCATE(field%evap_tile_ptr(ksfc_type))
    ALLOCATE(field%lhflx_tile_ptr(ksfc_type))
    ALLOCATE(field%shflx_tile_ptr(ksfc_type))
    ALLOCATE(field%dshflx_dT_tile_ptr(ksfc_type))

    DO jsfc = 1,ksfc_type

      CALL add_ref( field_list, prefix//'rsns_tile',                                   &
                  & prefix//'rsns_'//csfc(jsfc), field%swflxsfc_tile_ptr(jsfc)%p,      &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                &
                  & t_cf_var('rsns_'//csfc(jsfc), 'W m-2',                             &
                  &          'shortwave net flux at surface on tile '//csfc(jsfc),     &
                  &          DATATYPE_FLT32),                                          &
                  & t_grib2_var(0,4,9, ibits, GRID_REFERENCE, GRID_CELL),              &
                  & ldims=shape2d                                                      )

      CALL add_ref( field_list, prefix//'rlns_tile',                                   &
                  & prefix//'rlns_'//csfc(jsfc), field%lwflxsfc_tile_ptr(jsfc)%p,      &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                &
                  & t_cf_var('rlns_'//csfc(jsfc), 'W m-2',                             &
                  &          'longwave net flux at surface on tile '//csfc(jsfc),      &
                  &          DATATYPE_FLT32),                                          &
                  & t_grib2_var(0,5,5, ibits, GRID_REFERENCE, GRID_CELL),              &
                  & ldims=shape2d                                                      )

      CALL add_ref( field_list, prefix//'evspsbl_tile',                                &
                  & prefix//'evspsbl_'//csfc(jsfc), field%evap_tile_ptr(jsfc)%p,       &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                &
                  & t_cf_var('evspsbl_'//csfc(jsfc), 'kg m-2 s-1',                     &
                  &          'evaporation on tile '//csfc(jsfc), DATATYPE_FLT32),      &
                  & t_grib2_var(0,1,6, ibits, GRID_REFERENCE, GRID_CELL),              &
                  & ldims=shape2d                                                      )

      CALL add_ref( field_list, prefix//'hfls_tile',                                   &
                  & prefix//'hfls_'//csfc(jsfc), field%lhflx_tile_ptr(jsfc)%p,         &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                &
                  & t_cf_var('hfls_'//csfc(jsfc), 'W m-2',                             &
                  &          'latent heat flux on tile '//csfc(jsfc), DATATYPE_FLT32), &
                  & t_grib2_var(0,0,10, ibits, GRID_REFERENCE, GRID_CELL),             &
                  & ldims=shape2d                                                      )

      CALL add_ref( field_list, prefix//'hfss_tile',                                   &
                  & prefix//'hfss_'//csfc(jsfc), field%shflx_tile_ptr(jsfc)%p,         &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                &
                  & t_cf_var('hfss_'//csfc(jsfc), 'W m-2',                             &
                  &          'sensible heat flux on tile '//csfc(jsfc),DATATYPE_FLT32),&
                  & t_grib2_var(0,0,11, ibits, GRID_REFERENCE, GRID_CELL),             &
                  & ldims=shape2d                                                      )

      CALL add_ref( field_list, prefix//'dhfss_dT_tile',                               &
                  & prefix//'dhfss_dT_'//csfc(jsfc), field%dshflx_dT_tile_ptr(jsfc)%p, &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                &
                  & t_cf_var('dhfss_dT_'//csfc(jsfc), 'W m-2 K-1',                     &
                  &          'temp tend of SHF on tile '//csfc(jsfc), DATATYPE_FLT32), &
                  & t_grib2_var(0,0,11, ibits, GRID_REFERENCE, GRID_CELL),             &
                  & ldims=shape2d                                                      )
    END DO

    !-----------------------------------------
    ! wind stress, grid box mean
    !-----------------------------------------

    CALL add_var( field_list, prefix//'tauu', field%u_stress        ,           &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                           &
                & t_cf_var('u_stress','N m-2','u-momentum flux at the surface', &
                &          DATATYPE_FLT32),                                     &
                & t_grib2_var(0,2,17, ibits, GRID_REFERENCE, GRID_CELL),        &
                & ldims=shape2d,                                                &
                & lrestart = .FALSE.,                                           &
                & isteptype=TSTEP_INSTANT                                       )

    CALL add_var( field_list, prefix//'tauv', field%v_stress,                   &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                           &
                & t_cf_var('v_stress','N m-2','v-momentum flux at the surface', &
                &          DATATYPE_FLT32),                                     &
                & t_grib2_var(0,2,18, ibits, GRID_REFERENCE, GRID_CELL),        &
                & ldims=shape2d,                                                &
                & lrestart = .FALSE.,                                           &
                & isteptype=TSTEP_INSTANT                                       )

    ! wind stress, instantaneous tile values 

    CALL add_var( field_list, prefix//'tauu_tile', field%u_stress_tile,         &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                           &
                & t_cf_var('u_stress_tile', 'N m-2',                            &
                &          'u-momentum flux at the surface on tiles',           &
                &          DATATYPE_FLT32),                                     &
                & t_grib2_var(0,2,17, ibits, GRID_REFERENCE, GRID_CELL),        &
                & ldims=shapesfc,                                               &
                & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.          )

    CALL add_var( field_list, prefix//'tauv_tile', field%v_stress_tile,         &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                           &
                & t_cf_var('v_stress_tile', 'N m-2',                            &
                &          'v-momentum flux at the surface on tiles',           &
                &          DATATYPE_FLT32),                                     &
                & t_grib2_var(0,2,18, ibits, GRID_REFERENCE, GRID_CELL),        &
                & ldims=shapesfc,                                               &
                & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.          )

    ALLOCATE(field%u_stress_tile_ptr(ksfc_type))
    ALLOCATE(field%v_stress_tile_ptr(ksfc_type))

    DO jsfc = 1,ksfc_type

      CALL add_ref( field_list, prefix//'tauu_tile',                                &
                  & prefix//'tauu_'//csfc(jsfc), field%u_stress_tile_ptr(jsfc)%p,   &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                             &
                  & t_cf_var('u_stress_'//csfc(jsfc), 'N m-2',                      &
                  &          'u-momentum flux at the surface on tile '//csfc(jsfc), &
                  &          DATATYPE_FLT32),                                       &
                  & t_grib2_var(0,2,17, ibits, GRID_REFERENCE, GRID_CELL),          &
                  & ldims=shape2d                                                   )

      CALL add_ref( field_list, prefix//'tauv_tile',                                &
                  & prefix//'tauv_'//csfc(jsfc), field%v_stress_tile_ptr(jsfc)%p,   &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                             &
                  & t_cf_var('v_stress_'//csfc(jsfc), 'N m-2',                      &
                  &          'v-momentum flux at the surface on tile '//csfc(jsfc), &
                  &          DATATYPE_FLT32),                                       &
                  & t_grib2_var(0,2,18, ibits, GRID_REFERENCE, GRID_CELL),          &
                  & ldims=shape2d                                                   )
    END DO

  END SUBROUTINE new_echam_phy_field_list
  !-------------
  !>
  !!
  !!
  SUBROUTINE new_echam_phy_tend_list( k_jg, kproma, klev, kblks, ktracer, &
                                    & ctracer, listname, prefix,          &
                                    & tend_list, tend )

    INTEGER,INTENT(IN) :: k_jg !> patch ID
    INTEGER,INTENT(IN) :: kproma, klev, kblks, ktracer  !< dimension sizes

    CHARACTER(len=*)              ,INTENT(IN) :: listname, prefix
    CHARACTER(len=MAX_CHAR_LENGTH),INTENT(IN) :: ctracer(ktracer) !< tracer acronyms

    TYPE(t_var_list)      ,INTENT(INOUT) :: tend_list
    TYPE(t_echam_phy_tend),INTENT(INOUT) :: tend

    ! Local variables

    TYPE(t_cf_var)    ::    cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER :: shape3d(3), shape_trc(4)
    INTEGER :: ibits, jtrc
    !------------------------------

    ibits = DATATYPE_PACK16 ! "entropy" of horizontal slice

    shape3d   = (/kproma, klev, kblks/)
    shape_trc = (/kproma, klev, kblks, ktracer/)

    CALL new_var_list( tend_list, TRIM(listname), patch_id=k_jg )
    CALL default_var_list_settings( tend_list, lrestart=.FALSE. )

    !------------------------------
    ! Temperature tendencies
    !------------------------------
    ! &       tend% temp      (nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('temperature_tendency', 'K s-1',                               &
                &         'temperature tendency',                                        &
                &         DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0,0,255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( tend_list, prefix//'ta', tend%temp,                                    &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d, &
                & vert_interp=create_vert_interp_metadata(                               &
                &   vert_intp_type=vintp_types("P","Z","I"),                             &
                &   vert_intp_method=VINTP_METHOD_LIN,                                   &
                &   l_extrapol=.FALSE. ) )

    ! &       tend% temp_dyn  (nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('temperature_tendency_dyn', 'K s-1',                           &
                &         'temperature tendency due to  due to resolved dynamics',       &
                &         DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0,0,255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( tend_list, prefix//'ta_dyn', tend%temp_dyn,                            &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d, &
                & vert_interp=create_vert_interp_metadata(                               &
                &   vert_intp_type=vintp_types("P","Z","I"),                             &
                &   vert_intp_method=VINTP_METHOD_LIN,                                   &
                &   l_extrapol=.FALSE. ) )

    ! &       tend% temp_phy  (nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('temperature_tendency_phy', 'K s-1',                           &
                &         'temperature tendency due to parameterized processes',         &
                &         DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0,0,255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( tend_list, prefix//'ta_phy', tend%temp_phy,                            &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d, &
                & vert_interp=create_vert_interp_metadata(                               &
                &   vert_intp_type=vintp_types("P","Z","I"),                             &
                &   vert_intp_method=VINTP_METHOD_LIN,                                   &
                &   l_extrapol=.FALSE. ) )

    ! &       tend% temp_rsw(nproma,nlev,nblks),            &
    cf_desc    = t_cf_var('temperature_tendency_rsw', 'K s-1',                           &
                &         'temperature tendency due to shortwave radiation',             &
                &         DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0,0,255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( tend_list, prefix//'ta_rsw', tend%temp_rsw,                            &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d, &
                & vert_interp=create_vert_interp_metadata(                               &
                &   vert_intp_type=vintp_types("P","Z","I"),                             &
                &   vert_intp_method=VINTP_METHOD_LIN,                                   &
                &   l_extrapol=.FALSE. ) )

    ! &       tend% temp_rlw(nproma,nlev,nblks),            &
    cf_desc    = t_cf_var('temperature_tendency_rlw', 'K s-1',                           &
                &         'temperature tendency due to longwave radiation',              &
                &         DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0,0,255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( tend_list, prefix//'ta_rlw', tend%temp_rlw,                            &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d, &
                & vert_interp=create_vert_interp_metadata(                               &
                &   vert_intp_type=vintp_types("P","Z","I"),                             &
                &   vert_intp_method=VINTP_METHOD_LIN,                                   &
                &   l_extrapol=.FALSE. ) )

    ! &       tend% temp_rlw_impl(nproma,nblks),            &
    cf_desc    = t_cf_var('temperature_tendency_rlw_impl', 'K s-1',                      &
                &         'temperature tendency due to LW rad. due to implicit land surface temperature change', &
                &         DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0,0,255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( tend_list, prefix//'ta_rlw_impl', tend%temp_rlw_impl,                  &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,               &
                & ldims=(/kproma,kblks/))

    ! &       tend% temp_cld  (nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('temperature_tendency_cloud', 'K s-1',                         &
                &         'temperature tendency due to large scale cloud processes',     &
                &         DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0,0,255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( tend_list, prefix//'ta_cld', tend%temp_cld,                            &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d, &
                & vert_interp=create_vert_interp_metadata(                               &
                &   vert_intp_type=vintp_types("P","Z","I"),                             &
                &   vert_intp_method=VINTP_METHOD_LIN,                                   &
                &   l_extrapol=.FALSE. ) )

    ! &       tend% temp_cnv  (nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('temperature_tendency_convective', 'K s-1',                    &
                &         'temperature tendency due to convective cloud processes',      &
                &         DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0,0,255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( tend_list, prefix//'ta_cnv', tend%temp_cnv,                            &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d, &
                & vert_interp=create_vert_interp_metadata(                               &
                &   vert_intp_type=vintp_types("P","Z","I"),                             &
                &   vert_intp_method=VINTP_METHOD_LIN,                                   &
                &   l_extrapol=.FALSE. ) )

    ! &       tend% temp_vdf  (nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('temperature_tendency_turbulent', 'K s-1',                     &
                &         'temperature tendency due to vertical diffusion',              &
                &         DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0,0,255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( tend_list, prefix//'ta_vdf', tend%temp_vdf,                            &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d, &
                & vert_interp=create_vert_interp_metadata(                               &
                &   vert_intp_type=vintp_types("P","Z","I"),                             &
                &   vert_intp_method=VINTP_METHOD_LIN,                                   &
                &   l_extrapol=.FALSE. ) )

    ! &       tend% temp_gwh  (nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('temperature_tendency_Hines_gw', 'K s-1',                      &
                &         'temperature tendency due to non-orographic gravity waves',    &
                &         DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0,0,255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( tend_list, prefix//'ta_gwh', tend%temp_gwh,                            &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d, &
                & vert_interp=create_vert_interp_metadata(                               &
                &   vert_intp_type=vintp_types("P","Z","I"),                             &
                &   vert_intp_method=VINTP_METHOD_LIN,                                   &
                &   l_extrapol=.FALSE. ) )

    ! &       tend% temp_sso  (nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('temperature_tendency_sso', 'K s-1',                           &
                &         'temperature tendency due to sub grid scale orography',        &
                &         DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0,0,255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( tend_list, prefix//'ta_sso', tend%temp_sso,                            &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d, &
                & vert_interp=create_vert_interp_metadata(                               &
                &   vert_intp_type=vintp_types("P","Z","I"),                             &
                &   vert_intp_method=VINTP_METHOD_LIN,                                   &
                &   l_extrapol=.FALSE. ) )

    !------------------------------
    ! U-wind tendencies
    !------------------------------
    ! &       tend%    u      (nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('u_wind_tendency', 'm s-2',                                    &
                &         'u-wind tendency',                                             &
                &         DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0,2,255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( tend_list, prefix//'ua', tend%u,                                       &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d, &
                & vert_interp=create_vert_interp_metadata(                               &
                &   vert_intp_type=vintp_types("P","Z","I"),                             &
                &   vert_intp_method=VINTP_METHOD_LIN,                                   &
                &   l_extrapol=.FALSE. ) )

    ! &       tend%    u_dyn  (nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('u_wind_tendency_dyn', 'm s-2',                                &
                &         'u-wind tendency due to resolved dynamics',                    &
                &         DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0,2,255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( tend_list, prefix//'ua_dyn', tend%u_dyn,                               &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d, &
                & vert_interp=create_vert_interp_metadata(                               &
                &   vert_intp_type=vintp_types("P","Z","I"),                             &
                &   vert_intp_method=VINTP_METHOD_LIN,                                   &
                &   l_extrapol=.FALSE. ) )

    ! &       tend%    u_phy  (nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('u_wind_tendency_phy', 'm s-2',                                &
                &         'u-wind tendency due to parameterized processes',              &
                &         DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0,2,255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( tend_list, prefix//'ua_phy', tend%u_phy,                               &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d, &
                & vert_interp=create_vert_interp_metadata(                               &
                &   vert_intp_type=vintp_types("P","Z","I"),                             &
                &   vert_intp_method=VINTP_METHOD_LIN,                                   &
                &   l_extrapol=.FALSE. ) )

    ! &       tend%    u_cnv  (nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('u_wind_tendency_convective', 'm s-2',                         &
                &         'u-wind tendency due to convective cloud processes',           &
                &         DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0,2,255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( tend_list, prefix//'ua_cnv', tend%u_cnv,                               &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d, &
                & vert_interp=create_vert_interp_metadata(                               &
                &   vert_intp_type=vintp_types("P","Z","I"),                             &
                &   vert_intp_method=VINTP_METHOD_LIN,                                   &
                &   l_extrapol=.FALSE. ) )

    ! &       tend%    u_vdf  (nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('u_wind_tendency_turbulent', 'm s-2',                          &
                &         'u-wind tendency due to vertical diffusion',                   &
                &         DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0,2,255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( tend_list, prefix//'ua_vdf', tend%u_vdf,                               &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d, &
                & vert_interp=create_vert_interp_metadata(                               &
                &   vert_intp_type=vintp_types("P","Z","I"),                             &
                &   vert_intp_method=VINTP_METHOD_LIN,                                   &
                &   l_extrapol=.FALSE. ) )

    ! &       tend%    u_gwh  (nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('u_wind_tendency_nonoro_gw', 'm s-2',                          &
                &         'u-wind tendency due to non-orographic gravity waves',         &
                &         DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0,2,255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( tend_list, prefix//'ua_gwh', tend%u_gwh,                               &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d, &
                & vert_interp=create_vert_interp_metadata(                               &
                &   vert_intp_type=vintp_types("P","Z","I"),                             &
                &   vert_intp_method=VINTP_METHOD_LIN,                                   &
                &   l_extrapol=.FALSE. ) )

    ! &       tend%    u_sso  (nproma,nlev,nblks),          &   
    cf_desc    = t_cf_var('u_wind_tendency_sso', 'm s-2',                                &
                &         'u-wind tendency due to sub grid scale orography',             &
                &         DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0,2,255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( tend_list, prefix//'ua_sso', tend%u_sso,                               &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d, &
                & vert_interp=create_vert_interp_metadata(                               &
                &   vert_intp_type=vintp_types("P","Z","I"),                             &
                &   vert_intp_method=VINTP_METHOD_LIN,                                   &
                &   l_extrapol=.FALSE. ) )

    !------------------------------
    ! V-wind tendencies
    !------------------------------
    ! &       tend%    v      (nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('v_wind_tendency', 'm s-2',                                    &
                &         'v-wind tendency',                                             &
                &         DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0,2,255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( tend_list, prefix//'va', tend%v,                                       &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d, &
                & vert_interp=create_vert_interp_metadata(                               &
                &   vert_intp_type=vintp_types("P","Z","I"),                             &
                &   vert_intp_method=VINTP_METHOD_LIN,                                   &
                &   l_extrapol=.FALSE. ) )

    ! &       tend%    v_dyn  (nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('v_wind_tendency_dyn', 'm s-2',                                &
                &         'v-wind tendency due to resolved dynamics',                    &
                &         DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0,2,255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( tend_list, prefix//'va_dyn', tend%v_dyn,                               &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d, &
                & vert_interp=create_vert_interp_metadata(                               &
                &   vert_intp_type=vintp_types("P","Z","I"),                             &
                &   vert_intp_method=VINTP_METHOD_LIN,                                   &
                &   l_extrapol=.FALSE. ) )

    ! &       tend%    v_phy  (nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('v_wind_tendency_phy', 'm s-2',                                &
                &         'v-wind tendency due to parameterized processes',              &
                &         DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0,2,255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( tend_list, prefix//'va_phy', tend%v_phy,                               &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d, &
                & vert_interp=create_vert_interp_metadata(                               &
                &   vert_intp_type=vintp_types("P","Z","I"),                             &
                &   vert_intp_method=VINTP_METHOD_LIN,                                   &
                &   l_extrapol=.FALSE. ) )

    ! &       tend%    v_cnv  (nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('v_wind_tendency', 'm s-2',                                    &
                &         'v-wind tendency due to convective cloud processes',           &
                &         DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0,2,255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( tend_list, prefix//'va_cnv', tend%v_cnv,                               &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d, &
                & vert_interp=create_vert_interp_metadata(                               &
                &   vert_intp_type=vintp_types("P","Z","I"),                             &
                &   vert_intp_method=VINTP_METHOD_LIN,                                   &
                &   l_extrapol=.FALSE. ) )

    ! &       tend%    v_vdf  (nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('v_wind_tendency_turbulent', 'm s-2',                          &
                &         'v-wind tendency due to vertical diffusion',                   &
                &         DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0,2,255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( tend_list, prefix//'va_vdf', tend%v_vdf,                               &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d, &
                & vert_interp=create_vert_interp_metadata(                               &
                &   vert_intp_type=vintp_types("P","Z","I"),                             &
                &   vert_intp_method=VINTP_METHOD_LIN,                                   &
                &   l_extrapol=.FALSE. ) )

    ! &       tend%    v_gwh  (nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('v_wind_tendency_Hines_gw', 'm s-2',                           &
                &         'v-wind tendency due to non-orographic gravity waves',         &
                &         DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0,2,255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( tend_list, prefix//'va_gwh', tend%v_gwh,                               &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d, &
                & vert_interp=create_vert_interp_metadata(                               &
                &   vert_intp_type=vintp_types("P","Z","I"),                             &
                &   vert_intp_method=VINTP_METHOD_LIN,                                   &
                &   l_extrapol=.FALSE. ) )

    ! &       tend%    v_sso  (nproma,nlev,nblks),          &              
    cf_desc    = t_cf_var('v_wind_tendency_sso', 'm s-2',                                &
                &         'v-wind tendency due to sub grid scale orography',             &
                &         DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0,2,255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( tend_list, prefix//'va_sso', tend%v_sso,                               &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d, &
                & vert_interp=create_vert_interp_metadata(                               &
                &   vert_intp_type=vintp_types("P","Z","I"),                             &
                &   vert_intp_method=VINTP_METHOD_LIN,                                   &
                &   l_extrapol=.FALSE. ) )

    !------------------------------
    ! Detrainment
    !------------------------------

    ! &       tend%    xl_dtr  (nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('liquid_detrain_rate', 'kg kg-1 s-1',                          &
                &         'cloud water tendency due to detrainment from convective clouds', & 
                &         DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0,6,255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( tend_list, prefix//'clw_dtr', tend%xl_dtr,                              &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d, &
                & vert_interp=create_vert_interp_metadata(                               &
                &   vert_intp_type=vintp_types("P","Z","I"),                             &
                &   vert_intp_method=VINTP_METHOD_LIN,                                   &
                &   l_extrapol=.FALSE. ) )

   ! &       tend%    xi_dtr  (nproma,nlev,nblks),          &
   cf_desc    = t_cf_var('ice_detrain_rate', 'kg kg-1 s-1',                              &
                &        'cloud ice tendency due to detrainment from convective clouds', & 
                &         DATATYPE_FLT32)
   grib2_desc = t_grib2_var(0,6,255, ibits, GRID_REFERENCE, GRID_CELL)
   CALL add_var( tend_list, prefix//'cli_dtr', tend%xi_dtr,                               &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d, &
                & vert_interp=create_vert_interp_metadata(                               &
                &   vert_intp_type=vintp_types("P","Z","I"),                             &
                &   vert_intp_method=VINTP_METHOD_LIN,                                   &
                &   l_extrapol=.FALSE. ) )

    !-------------------
    ! Tracer tendencies
    !-------------------
    ! Tracer arrays for (model) internal use                                               

    CALL add_var( tend_list, prefix//'tracer', tend%q,                         &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID,                           &
                & t_cf_var('tend_tracer', 'kg kg-1 s-1',                       &
                &          'tendency of mass mixing ratio of tracers',         &
                &          DATATYPE_FLT32),                                    &
                & t_grib2_var(0,20,2, ibits, GRID_REFERENCE, GRID_CELL),       &
                & ldims = shape_trc,                                           &
                & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.         )

    CALL add_var( tend_list, prefix//'tracer_dyn', tend%q_dyn,                 &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID,                           &
                & t_cf_var('tend_tracer_dyn', 'kg kg-1 s-1',                   &
                &          'tendency of mass mixing ratio of tracers '//       &
                &          'due to resolved dynamics',                         &
                &          DATATYPE_FLT32),                                    &
                & t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),&
                & ldims = shape_trc,                                           &
                & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.         )

    CALL add_var( tend_list, prefix//'tracer_phy', tend%q_phy,                 &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID,                           &
                & t_cf_var('tend_tracer_phy', 'kg kg-1 s-1',                   &
                &          'tendency of mass mixing ratio of tracers '//       &
                &          'due to parameterized processes',                   &
                &          DATATYPE_FLT32),                                    &
                & t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),&
                & ldims = shape_trc,                                           &
                & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.         )

    CALL add_var( tend_list, prefix//'tracer_cld', tend%q_cld,                 &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID,                           &
                & t_cf_var('tend_tracer_cld', 'kg kg-1 s-1',                   &
                &          'tendency of mass mixing ratio of tracers '//       &
                &          'due to large scale cloud processes',               &
                &          DATATYPE_FLT32),                                    &           
                & t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),&
                & ldims = shape_trc,                                           &
                & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.         )

    CALL add_var( tend_list, prefix//'tracer_cnv', tend%q_cnv,                 &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID,                           &
                & t_cf_var('tend_tracer_cnv', 'kg kg-1 s-1',                   &
                &          'tendency of mass mixing ratio of tracers '//       &
                &          'due to convective cloud processes',                &
                &          DATATYPE_FLT32),                                    &           
                & t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),&
                & ldims = shape_trc,                                           &
                & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.         )

    CALL add_var( tend_list, prefix//'tracer_vdf', tend%q_vdf,                 &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID,                           &
                & t_cf_var('tend_tracer_vdf', 'kg kg-1 s-1',                   &
                &          'tendency of mass mixing ratio of tracers '//       &
                &          'due to vertical diffusion',                        &
                &          DATATYPE_FLT32),                                    &           
                & t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),&
                & ldims = shape_trc,                                           &
                & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.         )

    ! Referrence to individual tracer, for I/O

    ALLOCATE(tend%     q_ptr(ktracer))
    ALLOCATE(tend% q_dyn_ptr(ktracer))
    ALLOCATE(tend% q_phy_ptr(ktracer))
    ALLOCATE(tend% q_cld_ptr(ktracer))
    ALLOCATE(tend% q_cnv_ptr(ktracer))
    ALLOCATE(tend% q_vdf_ptr(ktracer))

    DO jtrc = 1,ktracer

      CALL add_ref( tend_list, prefix//'tracer',                                          &
                  & prefix//TRIM(ctracer(jtrc)), tend%q_ptr(jtrc)%p,                      &
                  & GRID_UNSTRUCTURED_CELL, ZA_HYBRID,                                    &
                  & t_cf_var('tend_'//TRIM(ctracer(jtrc)), 'kg kg-1 s-1',                 &
                  &          'tendency of mass mixing ratio of tracer '//                 &
                  &          TRIM(ctracer(jtrc)),                                         &
                  &          DATATYPE_FLT32),                                             &
                  & t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),         &
                  & ldims=(/kproma,klev,kblks/),                                          &
                  & vert_interp=create_vert_interp_metadata(                              &
                  &             vert_intp_type=vintp_types("P","Z","I"),                  &
                  &             vert_intp_method=VINTP_METHOD_LIN )                       )

      CALL add_ref( tend_list, prefix//'tracer_dyn',                                      &
                  & prefix//TRIM(ctracer(jtrc))//'_dyn', tend%q_dyn_ptr(jtrc)%p,          &
                  & GRID_UNSTRUCTURED_CELL, ZA_HYBRID,                                    &
                  & t_cf_var('tend_'//TRIM(ctracer(jtrc))//'_dyn', 'kg kg-1 s-1',         &
                  &          'tendency of mass mixing ratio of tracer '//                 &
                  &          TRIM(ctracer(jtrc))//                                        &
                  &          ' due to resolved dynamics',                                 &
                  &          DATATYPE_FLT32),                                             &
                  & t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),         &
                  & ldims=(/kproma,klev,kblks/),                                          &
                  & vert_interp=create_vert_interp_metadata(                              &
                  &             vert_intp_type=vintp_types("P","Z","I"),                  &
                  &             vert_intp_method=VINTP_METHOD_LIN )                       )

      CALL add_ref( tend_list, prefix//'tracer_phy',                                      &
                  & prefix//TRIM(ctracer(jtrc))//'_phy', tend%q_phy_ptr(jtrc)%p,          &
                  & GRID_UNSTRUCTURED_CELL, ZA_HYBRID,                                    &
                  & t_cf_var('tend_'//TRIM(ctracer(jtrc))//'_phy', 'kg kg-1 s-1',         &
                  &          'tendency of mass mixing ratio of tracer '//                 &
                  &          TRIM(ctracer(jtrc))//                                        &
                  &          ' due to parameterized processes',                           &
                  &          DATATYPE_FLT32),                                             &
                  & t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),         &
                  & ldims=(/kproma,klev,kblks/),                                          &
                  & vert_interp=create_vert_interp_metadata(                              &
                  &             vert_intp_type=vintp_types("P","Z","I"),                  &
                  &             vert_intp_method=VINTP_METHOD_LIN )                       )

      CALL add_ref( tend_list, prefix//'tracer_cld',                                      &
                  & prefix//TRIM(ctracer(jtrc))//'_cld', tend%q_cld_ptr(jtrc)%p,          &
                  & GRID_UNSTRUCTURED_CELL, ZA_HYBRID,                                    &
                  & t_cf_var('tend_'//TRIM(ctracer(jtrc))//'_cld', 'kg kg-1 s-1',         &
                  &          'tendency of mass mixing ratio of tracer '//                 &
                  &          TRIM(ctracer(jtrc))//                                        &
                  &          ' due to large scale cloud processes',                       &
                  &          DATATYPE_FLT32),                                             &
                  & t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),         &
                  & ldims=(/kproma,klev,kblks/),                                          &
                  & vert_interp=create_vert_interp_metadata(                              &
                  &             vert_intp_type=vintp_types("P","Z","I"),                  &
                  &             vert_intp_method=VINTP_METHOD_LIN )                       )

      CALL add_ref( tend_list, prefix//'tracer_cnv',                                      &
                  & prefix//TRIM(ctracer(jtrc))//'_cnv', tend%q_cnv_ptr(jtrc)%p,          &
                  & GRID_UNSTRUCTURED_CELL, ZA_HYBRID,                                    &
                  & t_cf_var('tend_'//TRIM(ctracer(jtrc))//'_cnv', 'kg kg-1 s-1',         &
                  &          'tendency of mass mixing ratio of tracer '//                 &
                  &          TRIM(ctracer(jtrc))//                                        &
                  &          ' due to convective cloud processes',                        &
                  &          DATATYPE_FLT32),                                             &
                  & t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),         &
                  & ldims=(/kproma,klev,kblks/),                                          &
                  & vert_interp=create_vert_interp_metadata(                              &
                  &             vert_intp_type=vintp_types("P","Z","I"),                  &
                  &             vert_intp_method=VINTP_METHOD_LIN )                       )

      CALL add_ref( tend_list, prefix//'tracer_vdf',                                      &
                  & prefix//TRIM(ctracer(jtrc))//'_vdf', tend%q_vdf_ptr(jtrc)%p,          &
                  & GRID_UNSTRUCTURED_CELL, ZA_HYBRID,                                    &
                  & t_cf_var('tend_'//TRIM(ctracer(jtrc))//'_vdf', 'kg kg-1 s-1',         &
                  &          'tendency of mass mixing ratio of tracer '//                 &
                  &          TRIM(ctracer(jtrc))//                                        &
                  &          ' due to vertical diffusion',                                &
                  &          DATATYPE_FLT32),                                             &
                  & t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),         &
                  & ldims=(/kproma,klev,kblks/),                                          &
                  & vert_interp=create_vert_interp_metadata(                              &
                  &             vert_intp_type=vintp_types("P","Z","I"),                  &
                  &             vert_intp_method=VINTP_METHOD_LIN )                       )
    END DO

  END SUBROUTINE new_echam_phy_tend_list
  !-------------

END MODULE mo_echam_phy_memory
