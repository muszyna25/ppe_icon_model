#if (defined (__GNUC__) || defined(__SUNPRO_F95) || defined(__SX__))
#define HAVE_F95
#endif
!>
!! Data types and variables used by the ECHAM6 physics package implemented
!! in the ICOHAM model.
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
!!      violent purposes.
!! <li> The code may not be re-distributed without the consent of the authors.
!! <li> The copyright notice and statement of authorship must appear in all
!!      copies.
!! <li> You accept the warranty conditions (see WARRANTY).
!! <li> In case you intend to use the code commercially, we oblige you to sign
!!      an according license agreement with DWD and MPI-M.
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
MODULE mo_echam_phy_memory

  USE mo_kind,                ONLY: wp
  USE mo_impl_constants,      ONLY: SUCCESS, MAX_CHAR_LENGTH
  USE mo_exception,           ONLY: message, finish
  USE mo_parallel_config,     ONLY: nproma
  USE mo_advection_config,    ONLY: advection_config
  USE mo_icoham_sfc_indices,  ONLY: nsfc_type
  USE mo_echam_phy_config,    ONLY: get_lvdiff, get_ljsbach, get_lamip
  USE mo_model_domain,        ONLY: t_patch

  USE mo_linked_list,         ONLY: t_var_list
  USE mo_var_list,            ONLY: default_var_list_settings, &
    &                               add_var, add_ref,          &
    &                               new_var_list,              &
    &                               delete_var_list
  USE mo_cf_convention,       ONLY: t_cf_var
  USE mo_grib2,               ONLY: t_grib2_var
  USE mo_cdi_constants,       ONLY: GRID_REFERENCE,                    &
    &                               GRID_UNSTRUCTURED_CELL, GRID_CELL, &
    &                               ZA_HYBRID, ZA_HYBRID_HALF,         &
    &                               ZA_SURFACE, ZA_GENERIC_ICE,        &
    &                               DATATYPE_PACK16, DATATYPE_PACK24,  &
    &                               DATATYPE_FLT32, FILETYPE_NC2,      &
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

  CHARACTER(len=*), PARAMETER :: version = '$Id$'
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
      ! shortwave surface albedo
      & albvisdir   (:,  :),  &!< [ ]    surface albedo for visible range, direct
      & albvisdif   (:,  :),  &!< [ ]    surface albedo for visible range, diffuse
      & albnirdir   (:,  :),  &!< [ ]    surface albedo for near IR range, direct
      & albnirdif   (:,  :),  &!< [ ]    surface albedo for near IR range, diffuse
      !
      ! shortwave net surface fluxes
      & vissfc      (:,  :),  &!< [ ]    solar transmissivity in VIS, net downward
      & visdffsfc   (:,  :),  &!< [ ]    diffuse fraction in VIS net downw. flux (?)
      & nirsfc      (:,  :),  &!< [ ]    solar transmissivity in NIR, net downward
      & nirdffsfc   (:,  :),  &!< [ ]    diffuse fraction in NIR net downw. flux (?)
      & parsfc      (:,  :),  &!< [ ]    solar transmissivity in PAR, downward
      & pardffsfc   (:,  :),  &!< [ ]    diffuse fraction in PAR net downw. flux (?)
      !
      ! shortwave net transmissivity
      & trsolclr    (:,:,:),  &!< [ ]    solar transmissivity  , clear sky, net downward
      & trsolall    (:,:,:),  &!< [ ]    solar transmissivity  , all   sky, net downward
      !
      ! longwave net fluxes
      & emterclr    (:,:,:),  &!< [W/m2] terrestrial emissivity, clear sky, net downward
      & emterall    (:,:,:),  & !< [W/m2] terrestrial emissivity, all   sky, net downward
      & o3          (:,:,:)     !< temporary set ozone mass mixing ratio  
    ! aerosol optical properties
    REAL(wp), POINTER ::      &
      & aer_aod_533 (:,:,:),  & !< aerosol optical depth at 533 nm
      & aer_ssa_533 (:,:,:),  & !< aerosol single scattering albedo at 533 nm
      & aer_asy_533 (:,:,:),  & !< aerosol asymmetry factor at 533 nm
      & aer_aod_2325(:,:,:),  & !< aerosol optical depth at 2325 nm
      & aer_ssa_2325(:,:,:),  & !< aerosol single scattering albedo at 2325 nm
      & aer_asy_2325(:,:,:),  & !< aerosol asymmetry factor at 2325 nm
      & aer_aod_9731(:,:,:)       !< effective aerosol optical depth at 9731 nm
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

    ! AMIP-sst and ice, preliminary 3 compartiments for water/ice/land for surface temp.
    REAL(wp),POINTER :: &
      & tsurfw (:,  :),     &!< sst as read in from amip input (==tsw)
      & tsurfi (:,  :),     &!< ice surface temperature
      & tsurfl (:,  :),     &!< land surface temperature
      & siced  (:,  :),     &!< ice depth
      & alake  (:,  :),     &!< lake mask
      & alb    (:,  :),     &!< surface background albedo
      & seaice (:,  :)       !< sea ice as read in from amip input

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
      & surface_temperature_rad  (:,  :),  &!< [K] radiative sfc. temperature
      & surface_temperature_eff  (:,  :),  &!< [K] effective sfc. temperature
      & zhsoil                   (:,  :),  &!< rel. humidity of land surface
      & csat                     (:,  :),  &!<
      & cair                     (:,  :)    !<

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
      & albnirdif_ice(:,:,:),   & ! Ice surface albedo for near IR range, diffuse
      & albvisdir_wtr(:,  :),   & ! Ocean surface albedo for visible range, direct
      & albvisdif_wtr(:,  :),   & ! Ocean surface albedo for visible range, diffuse
      & albnirdir_wtr(:,  :),   & ! Ocean surface albedo for near IR range, direct
      & albnirdif_wtr(:,  :)      ! Ocean surface albedo for near IR range, diffuse

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
      & swflxsfc    (:,:),  &!< [ W/m2] shortwave net flux at surface
      & swflxsfc_tile(:,:,:),  &!< [ W/m2] shortwave net flux at surface
      & lwflxsfc    (:,:),  &!< [ W/m2] longwave net flux at surface
      & lwflxsfc_tile(:,:,:),  &!< [ W/m2] longwave net flux at surface
      & dlwflxsfc_dT(:,:),  &!< [ W/m2/K] longwave net flux temp tend at surface
      & swflxtoa    (:,:),  &!< [ W/m2] shortwave net flux at TOA 
      & lwflxtoa    (:,:)    !< [ W/m2] shortwave net flux at TOA

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
!      & seaice(:,:),        &!< ice cover given as the fraction of (1- slm) (seaice in memory_g3b)
      & icefrc(:,:),        &!< ice cover given as the fraction of grid box (friac  in memory_g3b)
      & tsfc_tile (:,:,:),  &!< surface temperature over land/water/ice (tsw/l/i in memory_g3b)
      & tsfc      (:,  :),  &!< surface temperature, grid box mean
      & qs_sfc_tile(:,:,:)   !< saturation specitifc humidity at surface 

    TYPE(t_ptr2d),ALLOCATABLE ::   tsfc_tile_ptr(:)
    TYPE(t_ptr2d),ALLOCATABLE :: qs_sfc_tile_ptr(:)

    REAL(wp),POINTER :: &
      & lhflx     (:,  :),    &!< grid box mean latent   heat flux at surface 
      & shflx     (:,  :),    &!< grid box mean sensible heat flux at surface 
      & evap      (:,  :),    &!< grid box mean evaporation at surface 
      & lhflx_tile(:,:,:),    &!< (instantaneous) latent   heat flux at surface 
      & shflx_tile(:,:,:),    &!< (instantaneous) sensible heat flux at surface 
      & evap_tile(:,:,:),     &!< (instantaneous) evaporation at surface 
      & dshflx_dT_tile(:,:,:)  !< (instantaneous) temp tendency of SHF at surface

    TYPE(t_ptr2d),ALLOCATABLE :: lhflx_tile_ptr(:)
    TYPE(t_ptr2d),ALLOCATABLE :: shflx_tile_ptr(:)
    TYPE(t_ptr2d),ALLOCATABLE :: evap_tile_ptr(:)
    TYPE(t_ptr2d),ALLOCATABLE :: dshflx_dT_tile_ptr(:)

    REAL(wp),POINTER :: &
      & u_stress     (:,  :), &!< grid box mean wind stress 
      & v_stress     (:,  :), &!< grid box mean wind stress 
      & u_stress_tile(:,:,:), &!< (instantaneous) wind stress 
      & v_stress_tile(:,:,:)   !< (instantaneous) wind stress 

    TYPE(t_ptr2d),ALLOCATABLE :: u_stress_tile_ptr(:)
    TYPE(t_ptr2d),ALLOCATABLE :: v_stress_tile_ptr(:)

!!$    ! Variables for debugging
!!$    REAL(wp),POINTER :: &
!!$      ! 2d arrays
!!$      & debug_2d_1(:,:),    &!< 2d variable for debug purposes
!!$      & debug_2d_2(:,:),    &!< 2d variable for debug purposes
!!$      & debug_2d_3(:,:),    &!< 2d variable for debug purposes
!!$      & debug_2d_4(:,:),    &!< 2d variable for debug purposes
!!$      & debug_2d_5(:,:),    &!< 2d variable for debug purposes
!!$      & debug_2d_6(:,:),    &!< 2d variable for debug purposes
!!$      & debug_2d_7(:,:),    &!< 2d variable for debug purposes
!!$      & debug_2d_8(:,:),    &!< 2d variable for debug purposes
!!$      ! 3d arrays
!!$      & debug_3d_1(:,:,:),  &!< 3d variable for debug purposes
!!$      & debug_3d_2(:,:,:),  &!< 3d variable for debug purposes
!!$      & debug_3d_3(:,:,:),  &!< 3d variable for debug purposes
!!$      & debug_3d_4(:,:,:),  &!< 3d variable for debug purposes
!!$      & debug_3d_5(:,:,:),  &!< 3d variable for debug purposes
!!$      & debug_3d_6(:,:,:),  &!< 3d variable for debug purposes
!!$      & debug_3d_7(:,:,:),  &!< 3d variable for debug purposes
!!$      & debug_3d_8(:,:,:)    !< 3d variable for debug purposes

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
      &    u       (:,:,:)  , & !< accumulated tendency
      &    v       (:,:,:)  , & !< accumulated tendency
      & temp       (:,:,:)  , & !< accumulated tendency
      &    q       (:,:,:,:), & !< accumulated tendency
      !
      ! all physics processes
      !
      &    u_phy   (:,:,:)  , & !< accumulated tendency
      &    v_phy   (:,:,:)  , & !< accumulated tendency
      & temp_phy   (:,:,:)  , & !< accumulated tendency
      &    q_phy   (:,:,:,:), & !< accumulated tendency
      !
      ! cloud microphysics
      !
      & temp_cld   (:,:,:)  , & !< temperature tendency from cloud microphysical processes
      &    q_cld   (:,:,:,:), & !< tracer tendency from cloud microphysical process
      !
      ! cumulus convection
      !
      & temp_cnv    (:,:,:),   & !< temperature tendency from cumulus convection
      &    u_cnv    (:,:,:),   & !< u-wind tendency from cumulus convection
      &    v_cnv    (:,:,:),   & !< v-wind tendency from cumulus convection
      &    q_cnv    (:,:,:,:), & !< tracer tendency from cumulus convection
      &    xl_dtr   (:,:,:),   & !< cloud liquid tendency due to detrainment (memory_g3b:xtecl)
      &    xi_dtr   (:,:,:),   & !< cloud ice tendency due to detrainment (memory_g3b:xteci)
      !
      ! vertical turbulent mixing ("vdiff")
      !
      & temp_vdf   (:,:,:)  , & !< temperature tendency due to turbulent mixing
      &    u_vdf   (:,:,:)  , & !< u-wind tendency due to turbulent mixing
      &    v_vdf   (:,:,:)  , & !< v-wind tendency due to turbulent mixing
      &    q_vdf   (:,:,:,:), & !< tracer tendency due to turbulent mixing
      !
      ! Hines param. for atmospheric gravity waves
      !
      & u_gwh      (:,:,:)  , & !< u-wind tendency from Hines gravity wave param.
      & v_gwh      (:,:,:)  , & !< v-wind tendency from Hines gravity wave param.
      & temp_gwh   (:,:,:)  , & !< temperature tendency from Hines gravity wave param.
      !
      ! subgrid scale orographic (sso) blocking and gravity wave drag
      !
      & u_sso      (:,:,:)  , & !< u-wind tendency from sso drag
      & v_sso      (:,:,:)  , & !< v-wind tendency from sso drag
      & temp_sso   (:,:,:)  , & !< temperature tendency from sso drag
      !
      ! radiation
      !
      & temp_radsw (:,:,:)  , & !< temperature tendency from radiation
      & temp_radlw (:,:,:)      !< temperature tendency from radiation

    TYPE(t_ptr3d),ALLOCATABLE ::     q_ptr(:)
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
    CHARACTER(len=MAX_CHAR_LENGTH) :: & !< list of tracers to initialize
    &  ctracer_list
    INTEGER :: ndomain, jg, ist, nblks, nlev

    !---

    CALL message(TRIM(thismodule),'Construction of ECHAM physics state started.')

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

      ! get ctracer_list
      ctracer_list = advection_config(jg)%ctracer_list

      WRITE(listname,'(a,i2.2)') 'prm_field_D',jg
      CALL new_echam_phy_field_list( jg, nproma, nlev, nblks, ntracer, ctracer_list, &
                                   & nsfc_type, TRIM(listname), 'prm_',              &
                                   & prm_field_list(jg), prm_field(jg)          )

      WRITE(listname,'(a,i2.2)') 'prm_tend_D',jg
      CALL new_echam_phy_tend_list( jg, nproma, nlev, nblks, ntracer, ctracer_list, &
                                  & TRIM(listname), 'prm_tend_',                    &
                                  & prm_tend_list(jg), prm_tend(jg)             )
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
  SUBROUTINE new_echam_phy_field_list( k_jg, kproma, klev, kblks, ktracer,      &
                                     & ctracer_list, ksfc_type, listname, &
                                     & prefix, field_list, field          )

    INTEGER,INTENT(IN) :: k_jg !> patch ID
    INTEGER,INTENT(IN) :: kproma, klev, kblks, ktracer, ksfc_type  !< dimension sizes

    CHARACTER(len=*),INTENT(IN)    :: listname, prefix
    CHARACTER(len=MAX_CHAR_LENGTH) :: & !< list of tracers to initialize
    &  ctracer_list

    TYPE(t_var_list),       INTENT(INOUT) :: field_list
    TYPE(t_echam_phy_field),INTENT(INOUT) :: field

    ! Local variables

    TYPE(t_cf_var)    ::    cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER :: shape2d(2), shape3d(3), shapesfc(3), shapeice(3), shape3d_layer_interfaces(3)
!0!    INTEGER :: shape4d(4)
    INTEGER :: ibits, iextbits, jsfc, jtrc

    CHARACTER(LEN=1) :: csfc

    ibits = DATATYPE_PACK16
    iextbits = DATATYPE_PACK24

    shape2d  = (/kproma,       kblks/)
    shape3d  = (/kproma, klev, kblks/)
    shapesfc = (/kproma, kblks, ksfc_type/)
    shape3d_layer_interfaces = (/kproma,klev+1,kblks/)


    ! Register a field list and apply default settings

    CALL new_var_list( field_list, TRIM(listname), patch_id=k_jg )
    CALL default_var_list_settings( field_list,                &
                                  & lrestart=.TRUE.,           &
                                  & restart_type=FILETYPE_NC2  )

    !------------------------------
    ! Meteorological quantities
    !------------------------------

    ! &       field% u         (nproma,nlev  ,nblks),          &
    cf_desc    = t_cf_var('eastward_wind', 'm s-1', 'u-component of wind', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 2, 2, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'u', field%u,                             &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

    ! &       field% v         (nproma,nlev  ,nblks),          &
    cf_desc    = t_cf_var('northward_wind', 'm s-1', 'v-component of wind', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 2, 3, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'v', field%v,                             &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

    ! &       field% vor       (nproma,nlev  ,nblks),          &
    cf_desc    = t_cf_var('vorticity', 's-1', 'relative vorticity', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 2, 12, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'vo', field%vor,                          &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

    ! &       field% temp      (nproma,nlev  ,nblks),          &
    cf_desc    = t_cf_var('temperature', 'K', 'temperature', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 0, 0, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'t', field%temp,                          &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

    ! &       field% tv        (nproma,nlev  ,nblks),          &
    cf_desc    = t_cf_var('virtual_temperature', 'K', 'virtual temperature', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'vtemp', field%tv,                        &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

    ! OZONE 
    ! &       field% o3        (nproma,nlev  ,nblks),          &
    cf_desc    = t_cf_var('ozone', 'g/g', 'ozone mixing ratio', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'o3', field%o3,                        &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

    ! aerosol optical properties
    ! at 533 nm
    cf_desc    = t_cf_var('aer_aod_533','-','aerosol optical depth at 533 nm', &
                & DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'aer_aod_533', field%aer_aod_533,        &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc,      &
                & ldims=shape3d )
    cf_desc    = t_cf_var('aer_ssa_533','-',                                   &
                & 'aerosol single scattering albedo at 533 nm', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'aer_ssa_533', field%aer_ssa_533,        &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc,      &
                & ldims=shape3d )
    cf_desc    = t_cf_var('aer_asy_533','-',                                   &
                & 'aerosol asymmetry factor at 533 nm', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'aer_asy_533', field%aer_asy_533,        &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc,      &
                & ldims=shape3d )
    ! at 2325 nm
    cf_desc    = t_cf_var('aer_aod_2325','-',                                  &
                & 'aerosol optical depth at 2325 nm', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'aer_aod_2325', field%aer_aod_2325,      &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc,      &
                & ldims=shape3d )
    cf_desc    = t_cf_var('aer_ssa_2325','-',                                  &
                & 'aerosol single scattering albedo at 2325 nm', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'aer_ssa_2325', field%aer_ssa_2325,      &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc,      &
                & ldims=shape3d )
    cf_desc    = t_cf_var('aer_asy_2325','-',                                  &
                & 'aerosol asymmetry factor at 2325 nm', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'aer_asy_2325', field%aer_asy_2325,      &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc,      &
                & ldims=shape3d )
    ! at 9731 nm
    cf_desc    = t_cf_var('aer_aod_9731','-',                                  &
                & 'effective aerosol optical depth at 9731 nm', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'aer_aod_9731', field%aer_aod_9731,      &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc,      &
                & ldims=shape3d )
    ! &       field% q         (nproma,nlev  ,nblks,ntracer),  &
    CALL add_var( field_list, prefix//'q', field%q,                            &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID,                           &
                & t_cf_var('q', 'kg kg-1', '', DATATYPE_FLT32),                &
                & t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),&
                & ldims = (/kproma,klev,kblks,ktracer/),                       &
                & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.         )

    ALLOCATE(field%q_ptr(ktracer))                                                   
    DO jtrc = 1,ktracer                                                                   
      CALL add_ref( field_list, prefix//'q',                                     &
                  & prefix//'q_'//ctracer_list(jtrc:jtrc), field%q_ptr(jtrc)%p,  &
                  & GRID_UNSTRUCTURED_CELL, ZA_HYBRID,                           &
                  & t_cf_var('q_'//ctracer_list(jtrc:jtrc), 'kg kg-1', '', DATATYPE_FLT32), &
                  & t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),&
                  & ldims=(/kproma,klev,kblks/))                                        
    END DO                                                                                

    ! &       field% qx        (nproma,nlev  ,nblks),          &
    cf_desc    = t_cf_var('condensated_water', 'kg kg-1', 'cloud water + cloud ice', &
         &                DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'qx', field%qx,                           &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

    ! &       field% omega     (nproma,nlev  ,nblks),          &
    cf_desc    = t_cf_var('vertical_velocity', 'Pa s-1', 'vertical velocity', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'omega', field%omega,                     &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

    ! &       field% geom      (nproma,nlev  ,nblks),          &
    cf_desc    = t_cf_var('geopotential', 'm2 s-2', 'geopotential above surface', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 3, 5, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'ghm', field%geom,                        &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

    ! &       field% presm_old (nproma,nlev  ,nblks),          &
    cf_desc    = t_cf_var('pressure', 'Pa', 'pressure at old time step', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 3, 0, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'presm_old', field%presm_old,             &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

    ! &       field% presm_new (nproma,nlev  ,nblks),          &
    cf_desc    = t_cf_var('pressure', 'Pa', 'pressure at new time step', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 3, 0, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'presm_new', field%presm_new,             &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d )


    !-- Variables defined at layer interfaces --
    ! &       field% geoi      (nproma,nlevp1,nblks),          &
    cf_desc    = t_cf_var('geopotential', 'm2 s-2', 'geopotential above surface', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 3, 5, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'ghi', field%geoi,                         &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID_HALF, cf_desc, grib2_desc, ldims=shape3d_layer_interfaces )

    ! &       field% presi_old (nproma,nlevp1,nblks),          &
    cf_desc    = t_cf_var('pressure', 'Pa', 'pressure at old time step', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 3, 0, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'presi_old', field%presi_old,             &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID_HALF, cf_desc, grib2_desc, ldims=shape3d_layer_interfaces )

    ! &       field% presi_new (nproma,nlevp1,nblks),          &
    cf_desc    = t_cf_var('pressure', 'Pa', 'pressure at new time step', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 3, 0, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'presi_new', field%presi_new,             &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID_HALF, cf_desc, grib2_desc, ldims=shape3d_layer_interfaces  )

    !------------------
    ! Radiation
    !------------------
    ! 2D variables

   !ALLOCATE( field% cosmu0    (nproma,       nblks),          &
    cf_desc    = t_cf_var('cosmu0', '', '', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'cosmu0', field%cosmu0,                   &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    ! &       field% flxdwswtoa(nproma,       nblks),          &
    cf_desc    = t_cf_var('flxdwswtoa', '', '', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'flxdwswtoa', field%flxdwswtoa,           &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    ! &       field% albvisdir (nproma,       nblks),          &
    cf_desc    = t_cf_var('albvisdir', '', '', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'albvisdir', field%albvisdir,             &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    ! &       field% albvisdif (nproma,       nblks),          &
    cf_desc    = t_cf_var('albvisdif', '', '', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'albvisdif', field%albvisdif,             &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    ! &       field% albnirdir (nproma,       nblks),          &
    cf_desc    = t_cf_var('albnirdir', '', '', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'albnirdir', field%albnirdir,             &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    ! &       field% albnirdif (nproma,       nblks),          &
    cf_desc    = t_cf_var('albnirdif', '', '', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'albnirdif', field%albnirdif,             &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    ! &       field% vissfc    (nproma,       nblks),          &
    cf_desc    = t_cf_var('vissfc', '', '', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'vissfc', field%vissfc,                   &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    ! &       field% visdffsfc (nproma,       nblks),          &
    cf_desc    = t_cf_var('visdffsfc', '', '', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'visdffsfc', field%visdffsfc,             &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    ! &       field% nirsfc    (nproma,       nblks),          &
    cf_desc    = t_cf_var('nirsfc', '', '', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'nirsfc', field%nirsfc,                   &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    ! &       field% nirdffsfc (nproma,       nblks),          &
    cf_desc    = t_cf_var('nirdffsfc', '', '', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'nirdffsfc', field%nirdffsfc,             &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    ! &       field% parsfc    (nproma,       nblks),          &
    cf_desc    = t_cf_var('parsfc', '', '', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'parsfc', field%parsfc,                   &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    ! &       field% pardffsfc (nproma,       nblks),          &
    cf_desc    = t_cf_var('pardffsfc', '', '', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'pardffsfc', field%pardffsfc,             &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    cf_desc    = t_cf_var('swflxsfc', 'W m-2', ' shortwave net flux at surface', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 4, 9, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'swflxsfc', field%swflxsfc,&
         &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE,            &
         &        cf_desc, grib2_desc,                           &
         &        ldims=shape2d,                                 &
         &        lrestart = .FALSE.,                            &
         &        isteptype=TSTEP_INSTANT )
        
    cf_desc    = t_cf_var('swflxtoa', 'W m-2', ' shortwave net flux at TOA', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 4, 9, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'swflxtoa', field%swflxtoa,&
         &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE,            &
         &        cf_desc, grib2_desc,                           &
         &        ldims=shape2d,                                 &
         &        lrestart = .FALSE.,                            &
         &        isteptype=TSTEP_INSTANT )
        
    cf_desc    = t_cf_var('lwflxsfc', 'W m-2', 'longwave net flux at surface', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 5, 5, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'lwflxsfc', field%lwflxsfc,&
         &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE,            &
         &        cf_desc, grib2_desc,                           &
         &        ldims=shape2d,                                 &
         &        lrestart = .FALSE.,                            &
         &        isteptype=TSTEP_INSTANT )

    cf_desc    = t_cf_var('dlwflxsfc_dT', 'W m-2 K-1', 'longwave net flux T-derivative at surface', &
         &                DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 5, 5, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'dlwflxsfc_dT', field%dlwflxsfc_dT, &
         &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                     &
         &        cf_desc, grib2_desc,                                    &
         &        ldims=shape2d,                                          &
         &        lrestart = .FALSE.,                                     &
         &        isteptype=TSTEP_INSTANT )

    cf_desc    = t_cf_var('lwflxtoa', 'W m-2', 'longwave net flux at TOA', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 5, 5, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'lwflxtoa', field%lwflxtoa,&
         &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE,            &
         &        cf_desc, grib2_desc,                           &
         &        ldims=shape2d,                                 &
         &        lrestart = .FALSE.,                            &
         &        isteptype=TSTEP_INSTANT )

    IF (get_lamip()) THEN
    cf_desc    = t_cf_var('tsfc_wtr', 'K', 'surface temperature over water', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'tsfc_wtr', field%tsurfw,      &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    cf_desc    = t_cf_var('tsfc_ice', 'K', 'surface temperature over ice', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'tsfc_ice', field%tsurfi,      &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    cf_desc    = t_cf_var('tsfc_lnd', 'K', 'surface temperature over land', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'tsfc_lnd', field%tsurfl,      &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    cf_desc    = t_cf_var('siced', '', '', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'siced', field%siced,      &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    cf_desc    = t_cf_var('alb', '', '', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'alb', field%alb,      &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    cf_desc    = t_cf_var('oromea', '', '', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'oromea', field%oromea,      &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    cf_desc    = t_cf_var('orostd', '', '', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'orostd', field%orostd,      &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    cf_desc    = t_cf_var('orosig', '', '', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'orosig', field%orosig,      &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    cf_desc    = t_cf_var('orogam', '', '', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'orogam', field%orogam,      &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    cf_desc    = t_cf_var('orothe', '', '', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'orothe', field%orothe,      &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    cf_desc    = t_cf_var('oropic', '', '', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'oropic', field%oropic,      &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    cf_desc    = t_cf_var('oroval', '', '', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'oroval', field%oroval,      &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )
    END IF

    IF (get_ljsbach()) THEN

    cf_desc    = t_cf_var('tsfc_rad', '', '', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'tsfc_rad', field%surface_temperature_rad, &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    cf_desc    = t_cf_var('tsfc_eff', '', '', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'tsfc_eff', field%surface_temperature_eff, &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    cf_desc    = t_cf_var('zhsoil', '', '', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'zhsoil', field%zhsoil,      &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    cf_desc    = t_cf_var('csat', '', '', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'csat', field%csat,      &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    cf_desc    = t_cf_var('cair', '', '', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'cair', field%cair,      &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    END IF ! ljsbach

    !-------------------------
    ! Sea ice
    !-------------------------

    field%kice = kice ! Number of thickness classes - always 1, as of yet
    shapeice = (/kproma, field%kice, kblks/)

    CALL add_var( field_list, prefix//'Tsurf', field%Tsurf ,                  &
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE,                      &
      &          t_cf_var('Tsurf', 'C', 'surface temperature', DATATYPE_FLT32),&
      &          t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),&
      &          ldims=shapeice)
    CALL add_var( field_list, prefix//'T1', field%T1 ,                        &
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE,                      &
      &          t_cf_var('T1', 'C', 'Temperature upper layer', DATATYPE_FLT32),&
      &          t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),&
      &          ldims=shapeice)
    CALL add_var( field_list, prefix//'T2', field%T2 ,                        &
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE,                      &
      &          t_cf_var('T2', 'C', 'Temperature lower layer', DATATYPE_FLT32),&
      &          t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),&
      &          ldims=shapeice)
    CALL add_var( field_list, prefix//'hi', field%hi ,                        &
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE,                      &
      &          t_cf_var('hi', 'm', 'ice thickness', DATATYPE_FLT32),        &
      &          t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),&
      &          ldims=shapeice)
    CALL add_var( field_list, prefix//'hs', field%hs ,                        &
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE,                      &
      &          t_cf_var('hs', 'm', 'snow thickness', DATATYPE_FLT32),       &
      &          t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),&
      &          ldims=shapeice)
    CALL add_var( field_list, prefix//'Qtop', field%Qtop ,                    &
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE,                      &
      &          t_cf_var('Qtop', 'W/m^2', 'Energy flux available for surface melting', &
      &                   DATATYPE_FLT32),                                    &
      &          t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),&
      &          ldims=shapeice, lrestart=.FALSE.)
    CALL add_var( field_list, prefix//'Qbot', field%Qbot ,                    &
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE,                      &
      &          t_cf_var('Qbot', 'W/m^2', 'Energy flux at ice-ocean interface', DATATYPE_FLT32),&
      &          t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),&
      &          ldims=shapeice, lrestart=.FALSE.)


    CALL add_var( field_list, prefix//'conc', field%conc ,                    &
      &          GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE,                      &
      &          t_cf_var('conc', '', 'ice concentration in each ice class', DATATYPE_FLT32),&
      &          t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),&
      &          ldims=shapeice)

    ! &       field% albvisdir_ice (nproma,field%kice,nblks),          &
    cf_desc    = t_cf_var('albvisdir_ice', '', 'ice albedo VIS direct', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'albvisdir_ice', field%albvisdir_ice,             &
                & GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, cf_desc, grib2_desc, ldims=shapeice )

    ! &       field% albvisdif_ice (nproma,field%kice,nblks),          &
    cf_desc    = t_cf_var('albvisdif_ice', '', 'ice albedo VIS diffuse', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'albvisdif_ice', field%albvisdif_ice,             &
                & GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, cf_desc, grib2_desc, ldims=shapeice )

    ! &       field% albnirdir_ice (nproma,field%kice,nblks),          &
    cf_desc    = t_cf_var('albnirdir_ice', '', 'ice albedo NIR direct', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'albnirdir_ice', field%albnirdir_ice,             &
                & GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, cf_desc, grib2_desc, ldims=shapeice )

    ! &       field% albnirdif_ice (nproma,field%kice,nblks),          &
    cf_desc    = t_cf_var('albnirdif_ice', '', 'ice albedo NIR direct', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'albnirdif_ice', field%albnirdif_ice,             &
                & GRID_UNSTRUCTURED_CELL, ZA_GENERIC_ICE, cf_desc, grib2_desc, ldims=shapeice )


    !---- 3D variables defined at layer interfaces ----

    ! &       field% emterclr  (nproma,nlevp1,nblks),          &
    cf_desc    = t_cf_var('emterclr', '', '', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'emterclr', field%emterclr,               &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID_HALF, cf_desc, grib2_desc, ldims=shape3d_layer_interfaces )

    ! &       field% emterall  (nproma,nlevp1,nblks),          &
    cf_desc    = t_cf_var('emterall', '', '', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'emterall', field%emterall,               &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID_HALF, cf_desc, grib2_desc, ldims=shape3d_layer_interfaces )

    ! &       field% trsolclr  (nproma,nlevp1,nblks),          &
    cf_desc    = t_cf_var('trsolclr', '', '', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'trsolclr', field%trsolclr,               &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID_HALF, cf_desc, grib2_desc, ldims=shape3d_layer_interfaces )

    ! &       field% trsolall  (nproma,nlevp1,nblks),          &
    cf_desc    = t_cf_var('trsolall', '', '', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'trsolall', field%trsolall,               &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID_HALF, cf_desc, grib2_desc, ldims=shape3d_layer_interfaces )


    !-------------------------
    ! Cloud and precipitation
    !-------------------------
    cf_desc    = t_cf_var('ACLC', 'm2 m-2', 'cloud area fraction, instantaneous', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0,6,1, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'aclc', field%aclc,                       &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

    cf_desc    = t_cf_var('ACLCOV', 'm2 m-2', &
               & 'total cloud cover', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0,6, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'aclcov', field%aclcov,    &
         &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE,            &
         &        cf_desc, grib2_desc,                           &
         &        ldims=shape2d,                                 &
         &        lrestart = .FALSE.,                            &
         &        isteptype=TSTEP_INSTANT )

    ! &       field% acdnc  (nproma,nlev  ,nblks), &
    cf_desc    = t_cf_var('acdnc', '', '', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'acdnc', field%acdnc,                     &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

    ! &       field% xvar   (nproma,nlev  ,nblks), &
    cf_desc    = t_cf_var('variance_of_total_water', '', 'subgrid variance of total water', &
         &                DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'xvar', field%xvar,                       &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

    ! &       field% xskew  (nproma,nlev  ,nblks), &
    cf_desc    = t_cf_var('skewness_of_total_water', '', 'skewness of total water', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'xskew', field%xskew,                     &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

    ! &       field% relhum (nproma,nlev  ,nblks), &
    cf_desc    = t_cf_var('relative_humidity', '', 'relative humidity', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 1, 1, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'r', field%relhum,                        &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

    cf_desc    = t_cf_var('RSFL', 'kg m-2 s-1',    &
               & 'instantaneous large-scale precipitation flux (water)', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'rsfl', field%rsfl,        &
         &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE,            &
         &        cf_desc, grib2_desc,                           &
         &        ldims=shape2d,                                 &
         &        lrestart = .TRUE.,                             &
         &        isteptype=TSTEP_INSTANT )

    cf_desc    = t_cf_var('RSFC', 'kg m-2 s-1',    &
               & 'instantaneous convective precipitation flux (water)', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'rsfc', field%rsfc,        &
         &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE,            &
         &        cf_desc, grib2_desc,                           &
         &        ldims=shape2d,                                 &
         &        lrestart = .TRUE.,                             &
         &        isteptype=TSTEP_INSTANT )

    cf_desc    = t_cf_var('SSFL', 'kg m-2 s-1',    &
               & 'instantaneous large-scale precipitation flux (snow)', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'ssfl', field%ssfl,        &
         &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE,            &
         &        cf_desc, grib2_desc,                           &
         &        ldims=shape2d,                                 &
         &        lrestart = .TRUE.,                             &
         &        isteptype=TSTEP_INSTANT )

    cf_desc    = t_cf_var('SSFC', 'kg m-2 s-1',    &
               & 'instantaneous convective precipitation flux (snow)', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'ssfc', field%ssfc,        &
         &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE,            &
         &        cf_desc, grib2_desc,                           &
         &        ldims=shape2d,                                 &
         &        lrestart = .TRUE.,                             &
         &        isteptype=TSTEP_INSTANT )

    cf_desc    = t_cf_var('TOTPREC', 'kg m-2 s-1',                  &
         &                'instantaneous total precipitation flux', &
         &                DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 1, 52, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'totprec', field%totprec,  &
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
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'qvi', field%qvi,          &
         &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE,            &
         &        cf_desc, grib2_desc,                           &
         &        ldims=shape2d,                                 &
         &        lrestart = .FALSE.,                            &
         &        isteptype=TSTEP_INSTANT )

    cf_desc    = t_cf_var('total_cloud_water', 'kg m-2',&
               & 'vertically integrated cloud water', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'xlvi', field%xlvi,        &
         &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE,            &
         &        cf_desc, grib2_desc,                           &
         &        ldims=shape2d,                                 &
         &        lrestart = .FALSE.,                            &
         &        isteptype=TSTEP_INSTANT )

    cf_desc    = t_cf_var('total_cloud_ice', 'kg m-2',&
               & 'vertically integrated cloud ice', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'xivi', field%xivi,                       &
         &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE,            &
         &        cf_desc, grib2_desc,                           &
         &        ldims=shape2d,                                 &
         &        lrestart = .FALSE.,                            &
         &        isteptype=TSTEP_INSTANT )

    ! &       field% rintop (nproma,       nblks), &
    cf_desc    = t_cf_var('rintop', '', '', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'rintop', field%rintop,                   &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    ! &       field% rtype  (nproma,       nblks), &
    cf_desc    = t_cf_var('convection_type', '', 'convection_type (0...3)', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'rtype', field%rtype,                     &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    ! &       field% topmax (nproma,       nblks), &
    cf_desc    = t_cf_var('topmax', 'Pa', 'maximum height of convective cloud tops', &
         &                DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'topmax', field%topmax,                   &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    ! &       field% tke    (nproma,nlev  ,nblks), &
    cf_desc    = t_cf_var('turbulent_kinetic_energy', 'm2 s-2', 'turbulent kinetic energy', &
         &                DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'tke', field%tke,                         &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

    ! &       field% thvsig (nproma,       nblks), &
    cf_desc    = t_cf_var('thvsig', 'K', '', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'thvsig', field%thvsig,                   &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    !---------------------------
    ! Orographic wave drag diagnostics
    !---------------------------
    CALL add_var( field_list, prefix//'u_stress_sso', field%u_stress_sso,         &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                             &
                & t_cf_var('u_stress_sso', 'N/m2',                                &
                &          'zonal stress from subgrid scale orographic drag',     &
                &          DATATYPE_FLT32),                                       &
                & t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),   &
                & ldims=shape2d,                                                  &
                & lrestart = .FALSE.,                                             &
                & isteptype=TSTEP_INSTANT                                         )

    CALL add_var( field_list, prefix//'v_stress_sso', field%v_stress_sso,         &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                             &
                & t_cf_var('v_stress_sso', 'N/m2',                                &
                &          'meridional stress from subgrid scale orographic drag',&
                &          DATATYPE_FLT32),                                       &
                & t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),   &
                & ldims=shape2d,                                                  &
                & lrestart = .FALSE.,                                             &
                & isteptype=TSTEP_INSTANT                                         )

    CALL add_var( field_list, prefix//'dissipation_sso', field%dissipation_sso,   &
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
    IF (get_lvdiff()) THEN

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

      ! &       field% cfm_tile(nproma,nsfc_type,nblks), &
      CALL add_var( field_list, prefix//'cfm_tile', field%cfm_tile,              &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                          &
                  & t_cf_var('turb_exchng_coeff_momentum', '', '', DATATYPE_FLT32), &
                  & t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),&
                  & ldims=shapesfc,                                              &
                  & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.         )

      ALLOCATE(field%cfm_tile_ptr(ksfc_type))
      DO jsfc = 1,ksfc_type
        WRITE(csfc,'(i1)') jsfc 
        CALL add_ref( field_list, prefix//'cfm_tile',                              &
                    & prefix//'cfm_tile_'//csfc, field%cfm_tile_ptr(jsfc)%p,       &
                    & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                          &
                    & t_cf_var('turb_exchng_coeff_momentum_'//csfc, '', '', DATATYPE_FLT32),&
                    & t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),&
                    & ldims=shape2d )
      END DO

      !---------
      ! &       field% cfh    (nproma,nlev,     nblks), &
      cf_desc    = t_cf_var('turb_exchng_coeff_heat', '', '', DATATYPE_FLT32)
      grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( field_list, prefix//'cfh', field%cfh,                       &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

      ! &       field% cfh_tile(nproma,nsfc_type,nblks), &
      CALL add_var( field_list, prefix//'cfh_tile', field%cfh_tile,              &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                          &
                  & t_cf_var('turb_exchng_coeff_heat', '', '', DATATYPE_FLT32),  &
                  & t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),&
                  & ldims=shapesfc,                                              &
                  & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.         )

      ALLOCATE(field%cfh_tile_ptr(ksfc_type))
      DO jsfc = 1,ksfc_type
        WRITE(csfc,'(i1)') jsfc 
        CALL add_ref( field_list, prefix//'cfh_tile',                              &
                    & prefix//'cfh_tile_'//csfc, field%cfh_tile_ptr(jsfc)%p,       &
                    & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                          &
                    & t_cf_var('turb_exchng_coeff_heat_'//csfc, '', '', DATATYPE_FLT32), &
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

      ! &       field% z0m_tile(nproma,nsfc_type,nblks), &
      CALL add_var( field_list, prefix//'z0m_tile', field%z0m_tile,              &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                          &
                  & t_cf_var('z0m_tile', '', 'aerodynamic roughness length', DATATYPE_FLT32),&
                  & t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),&
                  & ldims=shapesfc,                                              &
                  & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.         )

      ALLOCATE(field%z0m_tile_ptr(ksfc_type))
      DO jsfc = 1,ksfc_type
        WRITE(csfc,'(i1)') jsfc 
        CALL add_ref( field_list, prefix//'z0m_tile',                              &
                    & prefix//'z0m_tile_'//csfc, field%z0m_tile_ptr(jsfc)%p,       &
                    & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                          &
                    & t_cf_var('z0m_tile_'//csfc, '','', DATATYPE_FLT32),          &
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
      cf_desc    = t_cf_var('fricktion_velocity', 'm s-1', 'friction velocity', DATATYPE_FLT32)
      grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( field_list, prefix//'ustar', field%ustar,                   &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

      ! &       field% kedisp (nproma,nblks),                &
      cf_desc    = t_cf_var('KE dissipation rate', '', '', DATATYPE_FLT32)
      grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( field_list, prefix//'kedisp', field%kedisp,                 &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

      ! &       field% ocu    (nproma,nblks),                &
      cf_desc    = t_cf_var('ocean_sfc_u', '', '', DATATYPE_FLT32)
      grib2_desc = t_grib2_var(255, 255, 255, iextbits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( field_list, prefix//'ocu', field%ocu,                       &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

      ! &       field% ocv    (nproma,nblks),                &
      cf_desc    = t_cf_var('ocean_sfc_v', '', '', DATATYPE_FLT32)
      grib2_desc = t_grib2_var(255, 255, 255, iextbits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( field_list, prefix//'ocv', field%ocv,                       &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

      ! &       field% albvisdir_wtr (nproma,       nblks),          &
      cf_desc    = t_cf_var('albvisdir_wtr', '', 'ocean albedo VIS direct', DATATYPE_FLT32)
      grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( field_list, prefix//'albvisdir_wtr', field%albvisdir_wtr,             &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

      ! &       field% albvisdif_wtr (nproma,       nblks),          &
      cf_desc    = t_cf_var('albvisdif_wtr', '', 'ocean albedo VIS diffuse', DATATYPE_FLT32)
      grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( field_list, prefix//'albvisdif_wtr', field%albvisdif_wtr,             &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

      ! &       field% albnirdir_wtr (nproma,       nblks),          &
      cf_desc    = t_cf_var('albnirdir_wtr', '', 'ocean albedo NIR direct', DATATYPE_FLT32)
      grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( field_list, prefix//'albnirdir_wtr', field%albnirdir_wtr,             &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

      ! &       field% albnirdif_wtr (nproma,       nblks),          &
      cf_desc    = t_cf_var('albnirdif_wtr', '', 'ocean albedo NIR direct', DATATYPE_FLT32)
      grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( field_list, prefix//'albnirdif_wtr', field%albnirdif_wtr,             &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    ENDIF ! get_lvdiff

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
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'seaice', field%seaice,                 &
              & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    ! &       field% alake (nproma, nblks),                 &
    cf_desc    = t_cf_var('alake', '', 'fraction of lakes', &
         &                DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'alake', field%alake,                 &
              & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    ! &       field% icefrc (nproma, nblks),                 &
    cf_desc    = t_cf_var('ice_cover', '', 'ice cover given as fraction of grid box', & 
         &                DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'icefrc', field%icefrc,                 &
              & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    !-----------------------------------
    ! &       field% tsfc(nproma,nblks), &
    cf_desc    = t_cf_var('surface_temperature', '', 'surface temperature', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, prefix//'tsfc', field%tsfc,                     &
              & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    ! &       field% tsfc_tile(nproma,nsfc_type,nblks), &
    CALL add_var( field_list, prefix//'tsfc_tile', field%tsfc_tile,            &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                          &
                & t_cf_var('tsfc_tile', '', 'skin temperature', DATATYPE_FLT32), &
                & t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),&
                & ldims=shapesfc,                                              &
!                & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.         )
                & lcontainer=.TRUE., lrestart=.FALSE.         )

    ALLOCATE(field%tsfc_tile_ptr(ksfc_type))
    DO jsfc = 1,ksfc_type
      WRITE(csfc,'(i1)') jsfc 
      CALL add_ref( field_list, prefix//'tsfc_tile',                             &
                  & prefix//'tsfc_tile_'//csfc, field%tsfc_tile_ptr(jsfc)%p,     &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                          &
                  & t_cf_var('tsfc_tile_'//csfc, '', '', DATATYPE_FLT32),        &
                  & t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),&
                  & ldims=shape2d )
    END DO
    !-----------------------------------

    ! &       field% qs_sfc_tile (nproma,nsfc_type, nblks), &
    CALL add_var( field_list, prefix//'qs_sfc_tile', field%qs_sfc_tile,        &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                          &
                & t_cf_var('qs_sfc_tile', '', '', DATATYPE_FLT32),             &
                & t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),&
                & ldims=shapesfc,                                              &
                & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.         )

    ALLOCATE(field%qs_sfc_tile_ptr(ksfc_type))
    DO jsfc = 1,ksfc_type
      WRITE(csfc,'(i1)') jsfc 
      CALL add_ref( field_list, prefix//'qs_sfc_tile',                           &
                  & prefix//'qs_sfc_tile_'//csfc, field%qs_sfc_tile_ptr(jsfc)%p, &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                          &
                  & t_cf_var('qs_sfc_tile_'//csfc, '', '', DATATYPE_FLT32),      &
                  & t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),&
                  & ldims=shape2d )
    END DO

    !---------------------------
    ! Surface fluxes
    !---------------------------
    ! Averaged gridbox mean

    CALL add_var( field_list, prefix//'evap', field%evap,                 &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                     &
                & t_cf_var('evap', 'kg m-2', 'evaporation',               &
                & DATATYPE_FLT32),                                        &
                & t_grib2_var(2,0,6,iextbits, GRID_REFERENCE, GRID_CELL), &
                & ldims=shape2d,                                          &
                & lrestart = .FALSE.,                                     &
                & isteptype=TSTEP_INSTANT                                 )

    CALL add_var( field_list, prefix//'lhflx', field%lhflx,               &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                     &
                & t_cf_var('lhflx', 'W m-2 ', 'latent heat flux',         &
                & DATATYPE_FLT32),                                        &
                & t_grib2_var(2, 0, 6, ibits, GRID_REFERENCE, GRID_CELL), &
                & ldims=shape2d,                                          &
                & lrestart = .FALSE.,                                     &
                & isteptype=TSTEP_INSTANT                                 )

    CALL add_var( field_list, prefix//'shflx', field%shflx,               &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                     &
                & t_cf_var('shflx', 'W m-2 ', 'sensible heat flux',       &
                & DATATYPE_FLT32),                                        &
                & t_grib2_var(2, 0, 6, ibits, GRID_REFERENCE, GRID_CELL), &
                & ldims=shape2d,                                          &
                & lrestart = .FALSE.,                                     &
                & isteptype=TSTEP_INSTANT                                 )

    !---------------------------------
    ! Instantaneous values over tiles

    CALL add_var( field_list, prefix//'swflxsfc_tile', field%swflxsfc_tile,&
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                     &
                & t_cf_var('swflxsfc_tile', '', '', DATATYPE_FLT32),      &
                & t_grib2_var(0, 4, 9, ibits, GRID_REFERENCE, GRID_CELL), &
                & ldims=shapesfc,                                         &
                & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.    )

    CALL add_var( field_list, prefix//'lwflxsfc_tile', field%lwflxsfc_tile,&
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                     &
                & t_cf_var('lwflxsfc_tile', '', '', DATATYPE_FLT32),      &
                & t_grib2_var(0, 5, 5, ibits, GRID_REFERENCE, GRID_CELL), &
                & ldims=shapesfc,                                         &
                & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.    )

    CALL add_var( field_list, prefix//'evap_tile', field%evap_tile,       &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                     &
                & t_cf_var('evap_tile', '', '', DATATYPE_FLT32),          &
                & t_grib2_var(2, 0, 6, ibits, GRID_REFERENCE, GRID_CELL), &
                & ldims=shapesfc,                                         &
                & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.    )

    CALL add_var( field_list, prefix//'lhflx_tile', field%lhflx_tile,     &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                     &
                & t_cf_var('lhflx_tile', 'W m-2', 'latent heat flux', DATATYPE_FLT32), &
                & t_grib2_var(2, 0, 6, ibits, GRID_REFERENCE, GRID_CELL), &
                & ldims=shapesfc,                                         &
                & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.    )

    CALL add_var( field_list, prefix//'shflx_tile', field%shflx_tile,     &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                     &
                & t_cf_var('shflx_tile', 'W m-2', 'sensible heat flux', DATATYPE_FLT32),  &
                & t_grib2_var(2, 0, 6, ibits, GRID_REFERENCE, GRID_CELL), &
                & ldims=shapesfc,                                         &
                & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.    )

    CALL add_var( field_list, prefix//'dshflx_dT_tile', field%dshflx_dT_tile,&
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                     &
                & t_cf_var('dshflx_dT_tile', 'W m-2', 'temp tend of SHF', DATATYPE_FLT32),  &
                & t_grib2_var(2, 0, 6, ibits, GRID_REFERENCE, GRID_CELL), &
                & ldims=shapesfc,                                         &
                & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.    )


    ALLOCATE(field%swflxsfc_tile_ptr(ksfc_type))
    ALLOCATE(field%lwflxsfc_tile_ptr(ksfc_type))
    ALLOCATE(field%evap_tile_ptr(ksfc_type))
    ALLOCATE(field%lhflx_tile_ptr(ksfc_type))
    ALLOCATE(field%shflx_tile_ptr(ksfc_type))
    ALLOCATE(field%dshflx_dT_tile_ptr(ksfc_type))

    DO jsfc = 1,ksfc_type
      WRITE(csfc,'(i1)') jsfc 

      CALL add_ref( field_list, prefix//'swflxsfc_tile',                            &
                  & prefix//'swflxsfc_tile_'//csfc, field%swflxsfc_tile_ptr(jsfc)%p,&
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                             &
                  & t_cf_var('swflxsfc_tile_'//csfc, '', '', DATATYPE_FLT32),       &
                  & t_grib2_var(2,0,6, ibits, GRID_REFERENCE, GRID_CELL),           &
                  & ldims=shape2d )

      CALL add_ref( field_list, prefix//'lwflxsfc_tile',                            &
                  & prefix//'lwflxsfc_tile_'//csfc, field%lwflxsfc_tile_ptr(jsfc)%p,&
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                             &
                  & t_cf_var('lwflxsfc_tile_'//csfc, '', '', DATATYPE_FLT32),       &
                  & t_grib2_var(2,0,6, ibits, GRID_REFERENCE, GRID_CELL),           &
                  & ldims=shape2d )

      CALL add_ref( field_list, prefix//'evap_tile',                             &
                  & prefix//'evap_tile_'//csfc, field%evap_tile_ptr(jsfc)%p,     &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                          &
                  & t_cf_var('evap_tile_'//csfc, '', '', DATATYPE_FLT32),        &
                  & t_grib2_var(2,0,6, ibits, GRID_REFERENCE, GRID_CELL),        &
                  & ldims=shape2d )

      CALL add_ref( field_list, prefix//'lhflx_tile',                            &
                  & prefix//'lhflx_tile_'//csfc, field%lhflx_tile_ptr(jsfc)%p,   &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                          &
                  & t_cf_var('lhflx_tile_'//csfc, '', '', DATATYPE_FLT32),       &
                  & t_grib2_var(2,0,6, ibits, GRID_REFERENCE, GRID_CELL),        &
                  & ldims=shape2d )

      CALL add_ref( field_list, prefix//'shflx_tile',                            &
                  & prefix//'shflx_tile_'//csfc, field%shflx_tile_ptr(jsfc)%p,   &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                          &
                  & t_cf_var('shflx_tile_'//csfc, '', '', DATATYPE_FLT32),       &
                  & t_grib2_var(2,0,6, ibits, GRID_REFERENCE, GRID_CELL),        &
                  & ldims=shape2d )

      CALL add_ref( field_list, prefix//'dshflx_dT_tile',                             &
                  & prefix//'dshflx_dT_tile_'//csfc, field%dshflx_dT_tile_ptr(jsfc)%p,&
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                               &
                  & t_cf_var('dshflx_dT_tile_'//csfc, '', '', DATATYPE_FLT32),        &
                  & t_grib2_var(2,0,6, ibits, GRID_REFERENCE, GRID_CELL),             &
                  & ldims=shape2d )
    END DO

    !-----------------------------------------
    ! wind stress, accumulated grid box mean
    !-----------------------------------------

    CALL add_var( field_list, prefix//'u_stress', field%u_stress        ,       &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                           &
                & t_cf_var('u_stress', 'N m-2', 'surface wind stress',          &
                &          DATATYPE_FLT32),                                     &
                & t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL), &
                & ldims=shape2d,                                                &
                & lrestart = .FALSE.,                                           &
                & isteptype=TSTEP_INSTANT )

    CALL add_var( field_list, prefix//'v_stress', field%v_stress,               &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                           &
                & t_cf_var('v_stress', 'N m-2', 'surface wind stress',          &
                &          DATATYPE_FLT32),                                     &
                & t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL), &
                & ldims=shape2d,                                                &
                & lrestart = .FALSE.,                                           &
                & isteptype=TSTEP_INSTANT )

    ! wind stress, instantaneous tile values 

    CALL add_var( field_list, prefix//'u_stress_tile', field%u_stress_tile,    &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                          &
                & t_cf_var('u_stress_tile', '', '', DATATYPE_FLT32),           &
                & t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),&
                & ldims=shapesfc,                                              &
                & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.         )

    CALL add_var( field_list, prefix//'v_stress_tile', field%v_stress_tile,    &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                          &
                & t_cf_var('v_stress_tile', '', '', DATATYPE_FLT32),           &
                & t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),&
                & ldims=shapesfc,                                              &
                & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.         )

    ALLOCATE(field%u_stress_tile_ptr(ksfc_type))
    ALLOCATE(field%v_stress_tile_ptr(ksfc_type))

    DO jsfc = 1,ksfc_type
      WRITE(csfc,'(i1)') jsfc 

      CALL add_ref( field_list, prefix//'u_stress_tile',                            &
                  & prefix//'u_stress_tile_'//csfc, field%u_stress_tile_ptr(jsfc)%p,&
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                             &
                  & t_cf_var('u_stress_tile_'//csfc, '', '', DATATYPE_FLT32),       &
                  & t_grib2_var(255,255,255, ibits, GRID_REFERENCE, GRID_CELL),     &
                  & ldims=shape2d )

      CALL add_ref( field_list, prefix//'v_stress_tile',                            &
                  & prefix//'v_stress_tile_'//csfc, field%v_stress_tile_ptr(jsfc)%p,&
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                             &
                  & t_cf_var('v_stress_tile_'//csfc, '', '', DATATYPE_FLT32),       &
                  & t_grib2_var(255,255,255, ibits, GRID_REFERENCE, GRID_CELL),     &
                  & ldims=shape2d )
    END DO

  END SUBROUTINE new_echam_phy_field_list
  !-------------
  !>
  !!
  !!
  SUBROUTINE new_echam_phy_tend_list( k_jg, kproma, klev, kblks, ktracer,   &
                                    & ctracer_list, listname, prefix,       &
                                    & tend_list, tend )

    INTEGER,INTENT(IN) :: k_jg !> patch ID
    INTEGER,INTENT(IN) :: kproma, klev, kblks, ktracer  !< dimension sizes

    CHARACTER(len=*),INTENT(IN)    :: listname, prefix
    CHARACTER(len=MAX_CHAR_LENGTH) :: & !< list of tracers to initialize
    &  ctracer_list

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
    cf_desc    = t_cf_var('temperature_tendency', 'K s-1', '', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( tend_list, prefix//'temp', tend%temp,                    &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

    ! &       tend% temp_phy  (nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('temperature_tendency_phy', 'K s-1', '', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( tend_list, prefix//'temp_phy', tend%temp_phy,            &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

    ! &       tend% temp_radsw(nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('temperature_tendency_radsw', 'K s-1', '', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( tend_list, prefix//'temp_radsw', tend%temp_radsw,        &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

    ! &       tend% temp_radlw(nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('temperature_tendency_radlw', 'K s-1', '', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( tend_list, prefix//'temp_radlw', tend%temp_radlw,        &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

    ! &       tend% temp_cld  (nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('temperature_tendency_cloud', 'K s-1', '', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( tend_list, prefix//'temp_cld', tend%temp_cld,            &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

    ! &       tend% temp_cnv  (nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('temperature_tendency_convective', 'K s-1', '', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( tend_list, prefix//'temp_cnv', tend%temp_cnv,            &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

    ! &       tend% temp_vdf  (nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('temperature_tendency_turbulent', 'K s-1', '', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( tend_list, prefix//'temp_vdf', tend%temp_vdf,            &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

    ! &       tend% temp_gwh  (nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('temperature_tendency_Hines_gw', 'K s-1', '', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( tend_list, prefix//'temp_gwh', tend%temp_gwh,            &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

    ! &       tend% temp_sso  (nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('temperature_tendency_sso', 'K s-1', '', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( tend_list, prefix//'temp_sso', tend%temp_sso,            &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

    !------------------------------
    ! U-wind tendencies
    !------------------------------
    ! &       tend%    u      (nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('u_wind_tendency', 'm s-2', '', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( tend_list, prefix//'u', tend%u,                          &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

    ! &       tend%    u_phy  (nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('u_wind_tendency_phy', 'm s-2', '', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( tend_list, prefix//'u_phy', tend%u_phy,                  &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

    ! &       tend%    u_cnv  (nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('u_wind_tendency_convective', 'm s-2', '', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( tend_list, prefix//'u_cnv', tend%u_cnv,                  &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

    ! &       tend%    u_vdf  (nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('u_wind_tendency_turbulent', 'm s-2', '', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( tend_list, prefix//'u_vdf', tend%u_vdf,                  &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

    ! &       tend%    u_gwh  (nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('u_wind_tendency_Hines_gw', 'm s-2', '', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( tend_list, prefix//'u_gwh', tend%u_gwh,                  &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

    ! &       tend%    u_sso  (nproma,nlev,nblks),          &   
    cf_desc    = t_cf_var('u_wind_tendency_sso', 'm s-2', '', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( tend_list, prefix//'u_sso', tend%u_sso,                  &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

    !------------------------------
    ! V-wind tendencies
    !------------------------------
    ! &       tend%    v      (nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('v_wind_tendency', 'm s-2', '', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( tend_list, prefix//'v', tend%v,                          &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

    ! &       tend%    v_phy  (nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('v_wind_tendency_phy', 'm s-2', '', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( tend_list, prefix//'v_phy', tend%v_phy,                          &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

    ! &       tend%    v_cnv  (nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('v_wind_tendency', 'm s-2', '', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( tend_list, prefix//'v_cnv', tend%v_cnv,                  &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

    ! &       tend%    v_vdf  (nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('v_wind_tendency_turbulent', 'm s-2', '', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( tend_list, prefix//'v_vdf', tend%v_vdf,                  &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

    ! &       tend%    v_gwh  (nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('v_wind_tendency_Hines_gw', 'm s-2', '', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( tend_list, prefix//'v_gwh', tend%v_gwh,                  &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

    ! &       tend%    v_sso  (nproma,nlev,nblks),          &              
    cf_desc    = t_cf_var('v_wind_tendency_sso', 'm s-2', '', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( tend_list, prefix//'v_sso', tend%v_sso,                  &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

    !------------------------------
    ! Detrainment
    !------------------------------

    ! &       tend%    xl_dtr  (nproma,nlev,nblks),          &
    cf_desc    = t_cf_var('liquid_detrain_rate', 's-1', '', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( tend_list, prefix//'xl_dtr', tend%xl_dtr,                  &
    & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

   ! &       tend%    xi_dtr  (nproma,nlev,nblks),          &
   cf_desc    = t_cf_var('ice_detrain_rate', 's-1', '', DATATYPE_FLT32)
   grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
   CALL add_var( tend_list, prefix//'xi_dtr', tend%xi_dtr,                  &
   & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

    !-------------------
    ! Tracer tendencies
    !-------------------
    ! Tracer arrays for (model) internal use                                               
                                                                                            
    CALL add_var( tend_list, prefix//'q', tend%q,                              &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID,                           &
                & t_cf_var('tend_q', 's-1', 'tracer tendency', DATATYPE_FLT32),&
                & t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),&
                & ldims = shape_trc,                                           &
                & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.         )

    CALL add_var( tend_list, prefix//'q_phy', tend%q_phy,                      &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID,                           &
                & t_cf_var('tend_q_phy', 's-1', 'tracer tendency phyiscs',     &
                & DATATYPE_FLT32),                                             &
                & t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),&
                & ldims = shape_trc,                                           &
                & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.         )

    CALL add_var( tend_list, prefix//'q_cld', tend%q_cld,                         &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID,                              &
                & t_cf_var('tend_q_cld', 's-1', 'tracer tendency condensational', &
                & DATATYPE_FLT32),                                                &           
                & t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),   &
                & ldims = shape_trc,                                              &
                & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.            )

    CALL add_var( tend_list, prefix//'q_cnv', tend%q_cnv,                      &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID,                           &
                & t_cf_var('tend_q_cnv', 's-1', 'tracer tendency convective',  &
                &          DATATYPE_FLT32),                                    &           
                & t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),&
                & ldims = shape_trc,                                           &
                & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.         )

    CALL add_var( tend_list, prefix//'q_vdf', tend%q_vdf,                      &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID,                           &
                & t_cf_var('tend_q_vdf', 's-1', 'tracer tendency turbulent',   &
                &          DATATYPE_FLT32),                                    &           
                & t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),&
                & ldims = shape_trc,                                           &
                & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.         )

    ! Referrence to individual tracer, for I/O                                            
                                                                                            
    ALLOCATE(tend%     q_ptr(ktracer))
    ALLOCATE(tend% q_phy_ptr(ktracer))
    ALLOCATE(tend% q_cld_ptr(ktracer))
    ALLOCATE(tend% q_cnv_ptr(ktracer))
    ALLOCATE(tend% q_vdf_ptr(ktracer))

    DO jtrc = 1,ktracer                                                                   
                                                                                          
      CALL add_ref( tend_list, prefix//'q',                                        &       
                  & prefix//'q'//ctracer_list(jtrc:jtrc), tend%q_ptr(jtrc)%p,      &       
                  & GRID_UNSTRUCTURED_CELL, ZA_HYBRID,                             &       
                  & t_cf_var('tend_q'//ctracer_list(jtrc:jtrc), 's-1', '',         &
                  &          DATATYPE_FLT32),                                      &       
                  & t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),  &       
                  & ldims=(/kproma,klev,kblks/))                                        
                                                                                          
      CALL add_ref( tend_list, prefix//'q_phy',                                           &       
                  & prefix//'q'//ctracer_list(jtrc:jtrc)//'_phy', tend%q_phy_ptr(jtrc)%p, &       
                  & GRID_UNSTRUCTURED_CELL, ZA_HYBRID,                                    &       
                  & t_cf_var('tend_q'//ctracer_list(jtrc:jtrc)//'_phy', 's-1', '',        &
                  &          DATATYPE_FLT32),                                             &       
                  & t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),         &       
                  & ldims=(/kproma,klev,kblks/))                                        
                                                                                          
      CALL add_ref( tend_list, prefix//'q_cld',                                           &       
                  & prefix//'q'//ctracer_list(jtrc:jtrc)//'_cld', tend%q_cld_ptr(jtrc)%p, &       
                  & GRID_UNSTRUCTURED_CELL, ZA_HYBRID,                                    &       
                  & t_cf_var('tend_q'//ctracer_list(jtrc:jtrc)//'_cld', 's-1', '',        &
                  &          DATATYPE_FLT32),                                             &       
                  & t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),         &       
                  & ldims=(/kproma,klev,kblks/))                                        

      CALL add_ref( tend_list, prefix//'q_cnv',                                           &       
                  & prefix//'q'//ctracer_list(jtrc:jtrc)//'_cnv', tend%q_cnv_ptr(jtrc)%p, &
                  & GRID_UNSTRUCTURED_CELL, ZA_HYBRID,                                    &
                  & t_cf_var('tend_q'//ctracer_list(jtrc:jtrc)//'_cnv', 's-1', '',        &
                  &          DATATYPE_FLT32),                                             &
                  & t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),         &
                  & ldims=(/kproma,klev,kblks/))                                        

      CALL add_ref( tend_list, prefix//'q_vdf',                                           &
                  & prefix//'q'//ctracer_list(jtrc:jtrc)//'_vdf', tend%q_vdf_ptr(jtrc)%p, &
                  & GRID_UNSTRUCTURED_CELL, ZA_HYBRID,                                    &
                  & t_cf_var('tend_q'//ctracer_list(jtrc:jtrc)//'_vdf', 's-1', '',        & 
                  &          DATATYPE_FLT32),                                             &
                  & t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),         &
                  & ldims=(/kproma,klev,kblks/))                                        
    END DO                                                                                

  END SUBROUTINE new_echam_phy_tend_list
  !-------------

END MODULE mo_echam_phy_memory
