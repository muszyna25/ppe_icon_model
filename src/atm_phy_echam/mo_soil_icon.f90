!+ Definitions for soil parameters
! 
Module mo_soil_icon

  ! 
  ! Description: 
  !   Variables that hold soil parameters for JSBACH
  ! 
  ! Current Code Owner: jsbach_admin
  ! 
  ! History: 
  !  
  ! Version   Date        Comment 
  ! -------   ----        ------- 
  ! 0.1       2001/07/01  Original code. Reiner Schnur
  ! 
  ! Code Description: 
  !   Language:           Fortran 90. 
  !   Software Standards: "European Standards for Writing and  
  !     Documenting Exchangeable Fortran 90 Code". 
  ! 
  ! Modules used: 
  ! 
!!$ TR  USE mo_netCDF, ONLY: FILE_INFO, NF_MAX_NAME, BELOWSUR, SOILLEV
!!$ TR  USE mo_jsbach, ONLY: debug, test_stream
!!$ TR  USE mo_jsbach_grid, ONLY: kstart, kend, nidx
!!$ TR  USE mo_linked_list, ONLY : t_stream
!!$ TR  USE mo_mpi, ONLY: p_parallel, p_parallel_io, p_bcast, p_io, p_pe, p_nprocs
!!$ TR  USE mo_io_units, ONLY : nout, nerr
  USE mo_kind,   ONLY: wp 
!!$ TR  USE mo_filename, ONLY : path_limit
  USE mo_exception, ONLY: finish, message !! for testing , message_text
!!$ TR  USE mo_test, ONLY: test
  USE mo_jsbach_hydrology, ONLY: update_surf_down, update_surf_up

  IMPLICIT NONE

  PRIVATE
!!$ TR  PUBLIC :: soil_type, soil_param_type, init_soil, update_soil, soil_diagnostics, get_soil_diag
  PUBLIC :: update_soil_icon

  REAL(wp), PARAMETER :: cdel(5) = (/0.065_wp,0.254_wp,0.913_wp,2.902_wp,5.700_wp/) ! thicknesses of soil layers [m] for 5 layer scheme (should be namelist var?)

CONTAINS

  SUBROUTINE update_soil_icon(ntiles, nidx, &
       delta_time, time_step_len, &
       time_steps_soil, &
!!$ TR       nidx, surface, soil, soil_param, useDynveg, &
       canopy_conductance_max, &
!!$ TR       root_depth, &
       lai, cdrag, t_Acoef, t_Bcoef, q_Acoef, q_Bcoef, air_temperature, air_moisture, &
       surface_pressure, windspeed, wind10, rad_longwave_down, rad_shortwave_net, &
       precip_rain, precip_snow, p_echam_zchl, &
!!$ TR       cair, csat, p_echam_zchl, zhsoil_avg, tte_corr_avg, canopy_conductance_limited, &
!!$ TR       glac_runoff_evap, surf_runoff_hd, drainage_hd &
       !! output for testing (hydrology)
       cair, csat, csat_transpiration, &
       moisture1, moisture2, moisture3, moisture4, moisture5, moisture_all, &
       sat_surface_specific_humidity, skin_reservoir, &
       snow_fract, snow, snow_canopy, snow_melt, snow_acc, snow_melt_acc, &
       glacier_runoff_acc, runoff_acc, drainage_acc, canopy_snow_fract, &
       !! output for testing (energy balance)
       surface_temperature, surface_temperature_old, &
       surface_temperature_rad, &
       evapotranspiration, &
       c_soil_temperature1, &
       c_soil_temperature2, &
       c_soil_temperature3, &
       c_soil_temperature4, &
       c_soil_temperature5, &
       d_soil_temperature1, &
       d_soil_temperature2, &
       d_soil_temperature3, &
       d_soil_temperature4, &
       d_soil_temperature5, &
       soil_temperature1, &
       soil_temperature2, &
       soil_temperature3, &
       soil_temperature4, &
       soil_temperature5, &
       heat_capacity,     &
       ground_heat_flux   &
       )

!!$ TR   USE mo_land_surface, ONLY: land_surface_type
!!$ TR   USE mo_jsbach_constants, ONLY : Gravity, RhoH2O, SpecificHeatDryAirConstPressure, &
!!$ TR        SpecificHeatVaporConstPressure, Emissivity, StefanBoltzmann, &
!!$ TR        LatentHeatVaporization, LatentHeatSublimation, GasConstantDryAir, &
!!$ TR        GasConstantWaterVapor, zsigfac, zepsec, zqsncr
!!$ TR   USE mo_constants,    ONLY: tmelt
!!$ TR    USE mo_atmosphere,   ONLY: sat_specific_humidity
!!$ TR   USE mo_utils,        ONLY: average_tiles

!!$ TR#ifndef STANDALONE
!!$ TR   USE mo_physc2, ONLY: cvdifts       ! Factor for time step weighting in ECHAM5
!!$ TR#endif
!!$ TR   USE mo_time_control, ONLY : lstart, delta_time, time_step_len
!!$ TR   USE mo_semi_impl, ONLY: eps
!!$ TR   USE mo_param_switches, ONLY: lsurf

    INTEGER, INTENT(in) :: ntiles
    INTEGER, INTENT(in) :: nidx
    REAL(wp), INTENT(in) :: delta_time
    REAL(wp), INTENT(in) :: time_step_len
!!$ TR    INTEGER, INTENT(in) :: nidx
!!$ TR    TYPE(land_surface_type), INTENT(in) :: surface
!!$ TR    TYPE(soil_type),    INTENT(inout) :: soil
!!$ TR    TYPE(soil_param_type), INTENT(in) :: soil_param
!!$ TR    LOGICAL,               INTENT(in) :: useDynveg
    REAL(wp), DIMENSION(:,:), INTENT(in) :: &  ! Dimension (nidx,ntiles)
         canopy_conductance_max                   ! Unstressed canopy resistance
    REAL(wp), DIMENSION(:,:), INTENT(in) :: &  ! Dimension (nidx,ntiles)
         lai
!!$ TR    REAL(wp), DIMENSION(:,:), INTENT(inout) :: &  ! Dimension (nidx,ntiles)
!!$ TR         canopy_snow, &                             ! Snow depth in canopy
!!$ TR    REAL(dp), DIMENSION(:,:), INTENT(out) :: &  ! Dimension (nidx,ntiles)
!!$ TR         canopy_conductance_limited                  ! resistance limited by water avaiability (water stress)
    REAL(wp), DIMENSION(:), INTENT(in) :: &    ! Dimension nidx
         cdrag, &
         t_Acoef, t_Bcoef, &
         q_Acoef, q_Bcoef, &
         air_temperature, &    ! Temperature at lowest atmospheric level
         air_moisture, &       ! Specific humidity at lowest atmospheric level
         surface_pressure, &   ! Surface pressure
         windspeed, &          ! Wind speed
         wind10, &             ! 10m wind speed
         rad_longwave_down, &  ! Longwave radiation down
         rad_shortwave_net, &  ! Net shortwave radiation
         precip_rain, &        ! Rainfall
         precip_snow, &        ! Snowfall
         p_echam_zchl

!!$ TR    REAL(dp), INTENT(in) :: root_depth(:,:,:)
    !! output for testing (hydrology)
    REAL(wp), DIMENSION(:), INTENT(inout) :: cair                !! area fraction with wet surface
    REAL(wp), DIMENSION(:), INTENT(inout) :: csat                !! area fraction with wet surface (air)
    REAL(wp), DIMENSION(:), INTENT(inout) :: csat_transpiration
    REAL(wp), DIMENSION(:), INTENT(inout) :: moisture1
    REAL(wp), DIMENSION(:), INTENT(inout) :: moisture2
    REAL(wp), DIMENSION(:), INTENT(inout) :: moisture3
    REAL(wp), DIMENSION(:), INTENT(inout) :: moisture4
    REAL(wp), DIMENSION(:), INTENT(inout) :: moisture5
    REAL(wp), DIMENSION(:), INTENT(inout) :: moisture_all
    REAL(wp), DIMENSION(:), INTENT(inout) :: sat_surface_specific_humidity
    REAL(wp), DIMENSION(:), INTENT(inout) :: skin_reservoir
    REAL(wp), DIMENSION(:), INTENT(inout) :: snow_fract
    REAL(wp), DIMENSION(:), INTENT(inout) :: snow
    REAL(wp), DIMENSION(:), INTENT(inout) :: snow_canopy
    REAL(wp), DIMENSION(:), INTENT(inout) :: snow_melt
    REAL(wp), DIMENSION(:), INTENT(inout) :: snow_acc
    REAL(wp), DIMENSION(:), INTENT(inout) :: snow_melt_acc
    REAL(wp), DIMENSION(:), INTENT(inout) :: glacier_runoff_acc
    REAL(wp), DIMENSION(:), INTENT(inout) :: runoff_acc
    REAL(wp), DIMENSION(:), INTENT(inout) :: drainage_acc
    REAL(wp), DIMENSION(:,:), INTENT(out) :: canopy_snow_fract
    !! output for testing (energy balance)
    REAL(wp), DIMENSION(:), INTENT(inout) ::  surface_temperature  ! Dimension (nidx)
    REAL(wp), DIMENSION(:), INTENT(inout) ::  surface_temperature_old  ! Dimension (nidx)
    REAL(wp), DIMENSION(:), INTENT(out)   ::  surface_temperature_rad  ! Dimension (nidx)
    REAL(wp), DIMENSION(:), INTENT(out)   ::  evapotranspiration
    REAL(wp), DIMENSION(:), INTENT(inout) ::  c_soil_temperature1
    REAL(wp), DIMENSION(:), INTENT(inout) ::  c_soil_temperature2
    REAL(wp), DIMENSION(:), INTENT(inout) ::  c_soil_temperature3
    REAL(wp), DIMENSION(:), INTENT(inout) ::  c_soil_temperature4
    REAL(wp), DIMENSION(:), INTENT(inout) ::  c_soil_temperature5
    REAL(wp), DIMENSION(:), INTENT(inout) ::  d_soil_temperature1
    REAL(wp), DIMENSION(:), INTENT(inout) ::  d_soil_temperature2
    REAL(wp), DIMENSION(:), INTENT(inout) ::  d_soil_temperature3
    REAL(wp), DIMENSION(:), INTENT(inout) ::  d_soil_temperature4
    REAL(wp), DIMENSION(:), INTENT(inout) ::  d_soil_temperature5
    REAL(wp), DIMENSION(:), INTENT(inout) ::  soil_temperature1
    REAL(wp), DIMENSION(:), INTENT(inout) ::  soil_temperature2
    REAL(wp), DIMENSION(:), INTENT(inout) ::  soil_temperature3
    REAL(wp), DIMENSION(:), INTENT(inout) ::  soil_temperature4
    REAL(wp), DIMENSION(:), INTENT(inout) ::  soil_temperature5
    REAL(wp), DIMENSION(:), INTENT(inout) ::  heat_capacity
    REAL(wp), DIMENSION(:), INTENT(inout) ::  ground_heat_flux
    REAL(wp), DIMENSION(:), INTENT(inout) ::  time_steps_soil
!!$ TR    REAL(dp), DIMENSION(:), INTENT(out) :: &  ! Dimension nidx
!!$ TR         cair, csat, zhsoil_avg ,tte_corr_avg, &
!!$ TR         glac_runoff_evap, surf_runoff_hd, drainage_hd ! INPUT for HD-Model in coupled ocean model
    ! Local variables
    REAL(wp), DIMENSION(nidx) :: &
         net_radiation, &                        ! Net radiation at surface
         dry_static_energy, &                    ! Surface dry static energy (= C_p * T )
         dry_static_energy_new, &                ! New dry static energy
         surface_qsat_new, &                     ! New surface saturated specific humidity
         transpiration, &                        ! transpiration by plants through stomata
!!$ TR         air_qsat, &                             ! Saturated specific humidity at lowest atmospheric level
         evapotranspiration_no_snow_skin, &      ! Evapotranspiration without that from snow and the skin reservoir
         zdqsl, &                                ! Sensitivity of saturated surface specific humidity wrt temperature
         zcpq, &                                 ! Conversion factor for humidity from dry static energy
!!$ TR         snow_avg, &                             ! Snow [m water equivalent] averaged over all tiles
         melt_water_excess                       ! water from snow melting which exceeds skin reservoir and infiltrates soil
!!$ TR         csat_transpiration                      ! fraction of grid box that contributes to transpiration
                                                 ! (considered to be completely wet)
!!$ TR    REAL(wp), DIMENSION(nidx) :: &
!!$ TR         sat_surf_specific_hum_avg, &
!!$ TR         ground_heat_flux_avg, &
!!$ TR         heat_capacity_avg, &
!!$ TR         surface_temperature_avg, soil_temperature_avg, &
!!$ TR         soil_moisture_avg

    REAL(wp), DIMENSION(nidx,ntiles) :: &
         glacier_depth, &
         skin_reservoir_max, &          
         wet_skin_fract, &
         glacier_precip_minus_evap, &       ! P-E for glaciers [m]
         surface_runoff, drainage, &        ! Surface runoff and drainage for HD model [m]
         qsat_fact, qair_fact, &            ! Factors for implicit coupling
         air_dry_static_energy_new, &
         air_moisture_new, &
         evaporation_pot, &                 ! Potential evaporation
         canopy_resistance, &               ! Water-limited canopy resistance
         canopy_conductance_limited, &      ! water limited canopy conductance
         soil_moisture_root, &              ! Soil moisture in root zone
         soil_moisture_root_max, &          ! Field capacity in root zone
         water_stress_factor, &             ! Water stress factor (=1: no stress, =0: infinite stress)
!!$ TR         relative_humidity, &               ! Relative humidity (Eq. 3.3.2.9 ECHAM3 Manual)
!!$ TR         zhsoil, &
         tte_corr, &
!!$ TR         sat_surface_specific_hum_old, &
         qsat_veg, qair_veg
!!$ TR         csat_tiles, cair_tiles, &
!!$ TR         qsat_transpiration, csat_transpiration_tiles
!!$ TR#ifdef __PGI
!!$ TR    REAL(wp), DIMENSION(nidx,soil%ntiles) :: moisture_max
!!$ TR#endif
    REAL(wp), DIMENSION(nidx,1,ntiles) :: soil_moisture ! substitute for soil%moisture
    REAL(wp), DIMENSION(nidx,5)        :: soil_layer_moisture ! substitute for soil%layer_moisture
    REAL(wp), DIMENSION(nidx,1,ntiles) :: root_depth
    REAL(wp), DIMENSION(nidx,1) :: soil_depth
    REAL(wp), DIMENSION(nidx,1) :: soil_MaxMoisture
!!$ TR    LOGICAL, DIMENSION(nidx,soil%ntiles) :: &
!!$ TR         soil_mask                          ! True if not glacier and not lake
!
    REAL(wp) :: zcons30, ztpfac2, ztpfac3, vtmpc2
!!$ TR    REAL(wp) :: hlp1
    INTEGER :: itile, i, j

! 5 layer scheme variables and variables for extra water balance diagnostics
    REAL(wp) :: FieldCapacity(nidx, 5)                ! Field capacity of soil layer i [m]
    INTEGER  :: jllog                                 ! Switch/Index of water balance test grid box
    REAL(wp) :: skin_res_evap_tile(nidx,ntiles)       ! Skin reservoir evaporation from update_surf_up per tile [m]
    REAL(wp) :: skin_reservoire_evap(nidx)            ! Total-Skin reservoir evaporation [m]
    REAL(wp) :: snow_evap(nidx)                       ! Total snow evaporation [kg/m**2/s]
    REAL(wp) :: snow_evap_tile(nidx,ntiles)           ! Evaporation over snow per tile [kg/m**2/s]
    REAL(wp) :: layer_moisture_tile(nidx, 5, ntiles)  ! Temporary field to distribute layer_moisture to different tiles [m]
    REAL(wp) :: RedEvapotransp(nidx,ntiles)           ! Evapotransp. without snow and the skin reservoir evap. per tile [m]
    REAL(wp) :: reduced_evap(nidx,ntiles)             ! Diagnostic evaporation reduction obtained from 5 layer scheme
                                                      ! --> is currently distributed between the soil layers [kg/m**2/s]
    REAL(wp) :: soil_water_pre(nidx)                  ! Soil water storages of previous time step
    REAL(wp) :: water_balance(nidx)                   ! Soil water balance within time step [mm]
    REAL(wp) :: dum(nidx)                             ! Array for various temporary data
!   The atmosphere sees only the average skin reservoir regarding evapotranspiration. 
!   Thus, tiling makes no sense and leads to errors in Back-calculations of bare soil evap.
    REAL(wp) :: SkinReservoire(nidx)                  ! Mean skin reservoir [m]
    REAL(wp) :: SkinReservoire_tile(nidx,ntiles) ! Mean skin reservoir distributed on tiles [m]
    REAL(wp) :: SnowDepthCan(nidx)                    ! Mean snow on canopy [m]
    REAL(wp) :: WetSkin_frac(nidx)                    ! Mean wet skin fraction
    REAL(wp) :: MaxSkinResCap(nidx)                   ! Mean maximum skin reservoir capacity[m]
    REAL(wp) :: RootDepth(nidx)
    REAL(wp) :: SoilDepth(nidx)
    REAL(wp) :: hyd_cond_sat(nidx)
    REAL(wp) :: PotMoisture(nidx)
    REAL(wp) :: bclapp(nidx)
    REAL(wp) :: SoilPorosity(nidx)
    REAL(wp) :: FieldCapacity_param(nidx)
    REAL(wp) :: WiltingPoint(nidx)
    REAL(wp) :: PoreSizeIndex(nidx)

!!$ TR#ifdef STANDALONE
!!$ TR    REAL(wp), PARAMETER :: cvdifts = 1.0_dp
!!$ TR#endif

    !! local variables for testing
    LOGICAL, DIMENSION(nidx,ntiles) :: is_present,is_glacier
    REAL(wp), DIMENSION(nidx) :: sat_surface_specific_hum_old
    REAL(wp), DIMENSION(nidx) :: surface_temperature_unfiltered
    INTEGER  ::  ntsoil,nsoil
    REAL(wp), DIMENSION(nidx) :: ThermalDiffusivity
    REAL(wp), DIMENSION(nidx) :: VolHeatCapacity
    REAL(wp)  ::  cvdifts
    REAL(wp)  ::  eps
    REAL(wp)  ::  Emissivity
    REAL(wp)  ::  StefanBoltzmann
    REAL(wp)  ::  Gravity
    REAL(wp)  ::  SpecificHeatDryAirConstPressure
    REAL(wp)  ::  SpecificHeatVaporConstPressure
    REAL(wp)  ::  LatentHeatVaporization
    REAL(wp)  ::  RhoH2O, tmelt
    REAL(wp)  ::  moist_crit_fract, moist_wilt_fract
    REAL(wp)  ::  skin_res_max, zqsncr, zepsec, zsigfac
    REAL(wp), DIMENSION(nidx)  ::  oro_std_dev
    REAL(wp)  ::  crit_snow_depth
    REAL(wp), DIMENSION(nidx) :: veg_ratio_max


    REAL(wp), DIMENSION(nidx,5)  ::  c_soil_temperature
    REAL(wp), DIMENSION(nidx,5)  ::  d_soil_temperature
    REAL(wp), DIMENSION(nidx,5)  ::  soil_temperature

    ntsoil = 5   ! number of thermal soil layers
    nsoil = 5    ! number of hydrological soil layers
    ThermalDiffusivity(1:nidx) = 7.e-7_wp
    VolHeatCapacity(1:nidx) = 2.e+6_wp
    cvdifts = 1.5_wp
    eps = 0.1_wp
    Emissivity = 0.996_wp
    StefanBoltzmann = 5.67e-8_wp
    Gravity = 9.80665_wp
    SpecificHeatDryAirConstPressure = 1005.46_wp
    SpecificHeatVaporConstPressure = 1869.46_wp
    LatentHeatVaporization = 2.5008e6_wp
    RhoH2O = 1000._wp            !! Density of liquid water (kg/m^3)
    tmelt = 273.15_wp
!!$ TR    ntiles = soil%ntiles
!!$ TR    kidx0 = kstart
!!$ TR    kidx1 = kend
!!$ TR    nroot_zones = SIZE(root_depth,2)

    zcons30 = 1._wp / (cvdifts*Gravity*time_step_len)
    ztpfac2 = 1._wp / cvdifts
    ztpfac3 = 1._wp - ztpfac2
    vtmpc2 = SpecificHeatVaporConstPressure / SpecificHeatDryAirConstPressure - 1._wp

    moist_crit_fract = 0.75_wp
    moist_wilt_fract = 0.35_wp
    skin_res_max     = 2.e-4_wp
    zqsncr  = 0.95_wp   ! inverse of equivalent water height when snow is considered to completely cover the ground
    zepsec  = 1.e-12_wp
    zsigfac = 0.15_wp
    crit_snow_depth  = 5.85036e-3_wp
    oro_std_dev(1:nidx) = 140.606_wp

!!$ TR    soil_mask = .NOT. surface%is_glacier(kidx0:kidx1,:) .AND. .NOT. surface%is_lake(kidx0:kidx1,:) &
!!$ TR         .AND. surface%is_present(kidx0:kidx1,:)
    !----------------------------------------------------------------------------------------
!!$ TR    IF(lstart) THEN
!!$ TR       csat(1:nidx) = .5_dp
!!$ TR       cair(1:nidx) = .5_dp
!!$ TR       csat_transpiration(1:nidx) = .5_dp
!!$ TR    ELSE
!!$ TR       csat(1:nidx) = soil%csat(kidx0:kidx1)
!!$ TR       cair(1:nidx) = soil%cair(kidx0:kidx1)
!!$ TR       csat_transpiration(1:nidx) = soil%csat_transpiration(kidx0:kidx1)
!!$ TR    END IF

    !----------------------------------------------------------------------------------------
    ! Soil moisture in root zone
!!$ TR   DO itile=1,ntiles
!!$ TR      soil_moisture_root(:,itile) = &
!!$ TR           calc_moist_root_zone(soil%moisture(kidx0:kidx1,:,itile), &
!!$ TR                                soil_param%Depth(kidx0:kidx1,:), root_depth(:,:,itile))
!!$ TR      soil_moisture_root_max(:,itile) = &
!!$ TR           calc_moist_root_zone(soil_param%MaxMoisture(kidx0:kidx1,:), &
!!$ TR                               soil_param%Depth(kidx0:kidx1,:), root_depth(:,:,itile))
!!$ TR   END DO
!!$ TR   soil_moisture_root = MERGE(soil_moisture_root, 0._dp, surface%is_vegetation(kidx0:kidx1,:))
!!$ TR   soil_moisture_root_max = MERGE(soil_moisture_root_max, 0._dp, surface%is_vegetation(kidx0:kidx1,:))

    root_depth(:,:,:) = 1._wp
    soil_depth(:,:) = 1._wp
    soil_MaxMoisture(:,:) = 0.5_wp
    veg_ratio_max(:) = 1._wp ! for testing

    soil_layer_moisture(1:nidx,1) = moisture1(1:nidx)
    soil_layer_moisture(1:nidx,2) = moisture2(1:nidx)
    soil_layer_moisture(1:nidx,3) = moisture3(1:nidx)
    soil_layer_moisture(1:nidx,4) = moisture4(1:nidx)
    soil_layer_moisture(1:nidx,5) = moisture5(1:nidx)

    DO itile=1,ntiles
       soil_moisture(1:nidx,1,itile) = moisture_all(1:nidx)
       soil_moisture_root(1:nidx,itile) =                               &
            calc_moist_root_zone(soil_moisture(1:nidx,:,itile), &
                                 soil_depth(1:nidx,:), root_depth(:,:,itile))
       soil_moisture_root_max(1:nidx,itile) =                           &
            calc_moist_root_zone(soil_MaxMoisture(1:nidx,:), &
                                 soil_depth(1:nidx,:), root_depth(:,:,itile))
    END DO

    !------------------------------------------------------------------------------------------

    ! Surface dry static energy 
    !    ... (loop 331 in vdiff)
!!$ TR    DO itile=1,ntiles
!!$ TR       soil%sat_surface_specific_humidity(kidx0:kidx1,itile) = &
!!$ TR            sat_specific_humidity(soil%surface_temperature(kidx0:kidx1,itile), surface_pressure(1:nidx))
!!$ TR       IF (ANY(soil%sat_surface_specific_humidity(kidx0:kidx1,itile) > 0.99_dp*HUGE(1._dp))) THEN
!!$ TR          CALL message('sat_specific_humidity', 'lookup table overflow')
!!$ TR       ENDIF
!!$ TR       sat_surface_specific_hum_old(1:nidx,itile) = soil%sat_surface_specific_humidity(kidx0:kidx1,itile)
!!$ TR    END DO

    sat_surface_specific_humidity(1:nidx) = &
         sat_specific_humidity(surface_temperature(1:nidx), surface_pressure(1:nidx))
    IF (ANY(sat_surface_specific_humidity(1:nidx) > 0.99_wp*HUGE(1._wp))) THEN
       CALL message('sat_specific_humidity', 'lookup table overflow')
    ENDIF
    sat_surface_specific_hum_old(1:nidx) = sat_surface_specific_humidity(1:nidx)

!!$ TR    CALL average_tiles(soil%surface_temperature(kidx0:kidx1,:),surface%is_present(kidx0:kidx1,:) &
!!$ TR                       .AND. .NOT. surface%is_lake(kidx0:kidx1,:), &
!!$ TR                       surface%cover_fract(kidx0:kidx1,:), surface_temperature_avg(:))
!!$ TR    CALL average_tiles(soil%sat_surface_specific_humidity(kidx0:kidx1,:),surface%is_present(kidx0:kidx1,:) &
!!$ TR         .AND. .NOT. surface%is_lake(kidx0:kidx1,:), &
!!$ TR         surface%cover_fract(kidx0:kidx1,:),sat_surf_specific_hum_avg(:))

    ! Modify minimum canopy resistance (no water stress) according to water limitation in the soil root zone 
    canopy_resistance(:,:) = 1.e20_wp
    water_stress_factor(:,:) = 0._wp
    canopy_conductance_limited(:,:) = 1._wp/1.e20_wp

    DO itile=1,ntiles
       DO i=1,nidx
!!$ TR    j=kidx0+i-1
!!$ TR    IF (surface%is_vegetation(j,itile)) THEN
       water_stress_factor(i,itile) = calc_water_stress_factor(soil_moisture_root(i,itile), &
                                      soil_moisture_root_max(i,itile), moist_crit_fract, &
                                      moist_wilt_fract)
       IF (water_stress_factor(i,itile) > EPSILON(1._wp) .AND. &
           canopy_conductance_max(i,itile) > EPSILON(1._wp) .AND. &
           air_moisture(i) .LE. sat_surface_specific_humidity(i)) THEN
          canopy_resistance(i,itile) = 1._wp / (canopy_conductance_max(i,itile) * &
             water_stress_factor(i,itile) + 1.e-20_wp)
       ELSE
          canopy_resistance(i,itile) = 1.e20_wp
       END IF
       canopy_conductance_limited(i,itile) = 1._wp/MAX(canopy_resistance(i,itile),1.e-20_wp)
!!$ TR    END IF
       END DO
    END DO

    ! Sensitivity of saturated surface specific humidity to temperature
!!$ TR    zdqsl = (sat_specific_humidity(surface_temperature_avg(:) + 0.001_dp, surface_pressure(1:nidx)) - &
!!$ TR             sat_surf_specific_hum_avg(:)) * 1000._dp
    zdqsl = (sat_specific_humidity(surface_temperature(1:nidx) + 0.001_wp, &
             surface_pressure(1:nidx)) - sat_surface_specific_humidity(:)) * 1000._wp

    ! Compute ECHAM-type skin reservoir and snow fraction (bare soil and soil under canopy)
!!$ TR    skin_reservoir_max(:,:) = MERGE(soil_options%SkinReservoirMax * (1._dp + lai(:,:) *               &
!!$ TR                              SPREAD(surface%veg_ratio_max(kidx0:kidx1),DIM=2,ncopies=ntiles)), 0._dp, &
!!$ TR                              .NOT. surface%is_glacier(kidx0:kidx1,:))
    skin_reservoir_max(:,:) = skin_res_max * (1._wp + lai(:,:) * &
                              SPREAD(veg_ratio_max(1:nidx),DIM=2,ncopies=ntiles))

    evaporation_pot(:,:) = 0._wp
    DO  itile=1,ntiles
      DO  i=1,nidx
!!$ TR       j=kidx0+i-1
        IF (skin_reservoir_max(i,itile) > EPSILON(1._wp)) THEN
           canopy_snow_fract(i,itile) = MIN(1._wp, snow_canopy(i) / skin_reservoir_max(i,itile))
        ELSE
           canopy_snow_fract(i,itile) = 0._wp
        END IF
!!$ TR        IF (soil_mask(i,itile)) THEN
!!$ TR           wet_skin_fract(i,itile) = MERGE(MIN(1._dp, soil%skin_reservoir(j,itile) / skin_reservoir_max(i,itile)), &
!!$ TR                0.0_dp, soil_options%SkinReservoirMax > EPSILON(1._dp))
           wet_skin_fract(i,itile) = MIN(1._wp, skin_reservoir(i) / skin_reservoir_max(i,itile))
          
           ! Roesch et al, 2002, Climate Dynamics
!!$ TR           soil%snow_fract(j,itile) = zqsncr * TANH(soil%snow(j,itile) * 100._dp) *            &
!!$ TR                SQRT(soil%snow(j,itile) * 1000._dp / (soil%snow(j,itile) * 1000._dp + zepsec + &
!!$ TR                zsigfac * surface%oro_std_dev(j)))
           snow_fract(i) = zqsncr * TANH(snow(i) * 100._wp) *            &
                SQRT(snow(i) * 1000._wp / (snow(i) * 1000._wp + zepsec + &
                zsigfac * oro_std_dev(i)))           

!!$ TR           IF (soil%snow_fract(j,itile) < EPSILON(1._dp) .AND. canopy_snow_fract(i,itile) > EPSILON(1._dp)) &
!!$ TR                soil%snow_fract(j,itile) = canopy_snow_fract(i,itile)
           IF (snow_fract(i) < EPSILON(1._wp) .AND. canopy_snow_fract(i,itile) > EPSILON(1._wp)) &
                snow_fract(i) = canopy_snow_fract(i,itile)
           ! Modify snow cover if snow loss during the time step due to potential evaporation is larger than 
           ! equivalent snow water content from soil and canopy; same for skin reservoir
           ! Potential evaporation using old values of air and surface humidity
!!$ TR           evaporation_pot(i,itile) = zcons30 * cdrag(i) * &
!!$ TR                (soil%sat_surface_specific_humidity(j,itile) - air_moisture(i))
           evaporation_pot(i,itile) = zcons30 * cdrag(i) * &
                (sat_surface_specific_humidity(i) - air_moisture(i))
!!$ TR           IF (soil%snow_fract(j,itile) > 0._dp) THEN
!!$ TR              soil%snow_fract(j,itile) = soil%snow_fract(j,itile) /                              &
!!$ TR                   MAX(1._dp, soil%snow_fract(j,itile) * evaporation_pot(i,itile) * delta_time / &
!!$ TR                   (RhoH2O*(soil%snow(j,itile) + canopy_snow(i,itile))))
!!$ TR           END IF
           IF (snow_fract(i) > 0._wp) THEN
              snow_fract(i) = snow_fract(i) / &
                   MAX(1._wp, snow_fract(i) * evaporation_pot(i,itile) * delta_time / &
                   (RhoH2O * (snow(i) + snow_canopy(i))))
           END IF
!!$ TR           IF (wet_skin_fract(i,itile) > 0._dp) THEN
!!$ TR              wet_skin_fract(i,itile) = wet_skin_fract(i,itile) /                                          &
!!$ TR                   MAX(1._dp, (1._dp - soil%snow_fract(j,itile)) * evaporation_pot(i,itile) * delta_time / &
!!$ TR                   (RhoH2O * MAX(EPSILON(1._dp),soil%skin_reservoir(j,itile))))
!!$ TR           END IF
           IF (wet_skin_fract(i,itile) > 0._wp) THEN
              wet_skin_fract(i,itile) = wet_skin_fract(i,itile) / &
                   MAX(1._wp, (1._wp - snow_fract(i)) * evaporation_pot(i,itile) * delta_time / &
                   (RhoH2O * MAX(EPSILON(1._wp),skin_reservoir(i))))
           END IF
           
!!$ TR        ELSE IF (surface%is_glacier(j,itile)) THEN
!!$ TR           wet_skin_fract(i,itile)  = 0.0_dp
!!$ TR           soil%snow_fract(j,itile) = 1.0_dp
!!$ TR        ELSE
!!$ TR           wet_skin_fract(i,itile)  = 0.0_dp
!!$ TR           soil%snow_fract(j,itile) = 0.0_dp
!!$ TR        END IF
      END DO
    END DO

    !----------------------------------------------------------------------------------------------------------------------

    ! Note: at the moment, surface_temperature and sat_surface_specific_humidity are the same for all tiles in a grid box
!!$ TR    dry_static_energy(:) = surface_temperature_avg(:) * SpecificHeatDryAirConstPressure * &
!!$ TR         (1._dp+ vtmpc2 * ( csat(:) * sat_surf_specific_hum_avg(:) + (1._dp - cair(:)) * air_moisture(:)))
    dry_static_energy(1:nidx) = surface_temperature(1:nidx) * SpecificHeatDryAirConstPressure * &
        (1._wp + vtmpc2 * ( csat(1:nidx) * sat_surface_specific_humidity(1:nidx) + &
        (1._wp - cair(1:nidx)) * air_moisture(1:nidx)))
    
    ! Conversion factor for humidity from dry static energy
    zcpq(1:nidx) = dry_static_energy(1:nidx) / surface_temperature(1:nidx)

    net_radiation(1:nidx) = rad_shortwave_net(1:nidx) + &
         rad_longwave_down(1:nidx)  - Emissivity * StefanBoltzmann * &
         surface_temperature(1:nidx)**4

    ! Compute new surface temperature and moisture
    dry_static_energy_new = 0._wp
    surface_qsat_new      = 0._wp
!! ------------------------------------------------------
!!$ TR    CALL average_tiles(soil%ground_heat_flux(kidx0:kidx1,:),surface%is_present(kidx0:kidx1,:) &
!!$ TR         .AND. .NOT. surface%is_lake(kidx0:kidx1,:), &
!!$ TR         surface%cover_fract(kidx0:kidx1,:),ground_heat_flux_avg(:))
!!$ TR    CALL average_tiles(soil%heat_capacity(kidx0:kidx1,:),surface%is_present(kidx0:kidx1,:) &
!!$ TR         .AND. .NOT. surface%is_lake(kidx0:kidx1,:) , &
!!$ TR         surface%cover_fract(kidx0:kidx1,:),heat_capacity_avg(:))

!!$ TR    IF (.NOT. lstart) THEN
!!$ TR       CALL update_surfacetemp(nidx, zcpq(:),                                         &
!!$ TR            t_Acoef(:), t_Bcoef(:), q_Acoef(:), q_Bcoef(:),                                  &
!!$ TR            dry_static_energy(:), sat_surf_specific_hum_avg(:), zdqsl(:),            &
!!$ TR            net_radiation(:), ground_heat_flux_avg(:),                                       &
!!$ TR            zcons30*cdrag(:), cair(:), csat(:),                                          &
!!$ TR            SUM(surface%cover_fract(kidx0:kidx1,:) * soil%snow_fract(kidx0:kidx1,:), DIM=2), &
!!$ TR            heat_capacity_avg(:),                                                            &
!!$ TR            dry_static_energy_new(:), surface_qsat_new(:))
!!$ TR    ELSE
!!$ TR       dry_static_energy_new(:) = dry_static_energy(:)
!!$ TR       surface_qsat_new(:) = sat_surf_specific_hum_avg(:)
!!$ TR    ENDIF

!!$ TR    IF (.NOT. lstart) THEN !! lstart ???
       IF (time_steps_soil(1) > 0.5_wp) THEN
       CALL update_surfacetemp_icon(nidx, time_step_len,                              &
            cvdifts, zcpq(:),                                                         &
            t_Acoef(:), t_Bcoef(:), q_Acoef(:), q_Bcoef(:),                           &
            dry_static_energy(:), sat_surface_specific_humidity(:), zdqsl(:),         &
            net_radiation(:), ground_heat_flux(:),                                    &
            zcons30*cdrag(:), cair(:), csat(:),                                       &
!!$ TR            SUM(surface%cover_fract(kidx0:kidx1,:) * &
!!$ TR            soil%snow_fract(kidx0:kidx1,:), DIM=2), &
            snow_fract(:),                                                            &
            heat_capacity(:),                                                         &
            dry_static_energy_new(:), surface_qsat_new(:))           
       ELSE
          dry_static_energy_new(:) = dry_static_energy(:)
!!$ TR          surface_qsat_new(:) = sat_surf_specific_hum_avg(:)
          surface_qsat_new(:) = sat_surface_specific_humidity(:)
       ENDIF

    ! New land surface temperature and surface saturated humidity
!!$ TR    DO itile=1,ntiles
!!$ TR       IF (.NOT. lstart) THEN
!!$ TR          soil%sat_surface_specific_humidity(kidx0:kidx1,itile) = surface_qsat_new(:)
!!$ TR       END IF
!!$ TR       soil%dry_static_energy_new(kidx0:kidx1,itile) = dry_static_energy_new(:)
!!$ TR    END DO
    IF (time_steps_soil(1) > 0.5_wp) THEN
       sat_surface_specific_humidity(:) = surface_qsat_new(:)
    END IF 

    ! New unfiltered land surface temperature at t+dt (\tilde(X)^(t+1))
!!$ TR    soil%surface_temperature_unfiltered(kidx0:kidx1,:) = &
!!$ TR         SPREAD((ztpfac2 * dry_static_energy_new(:) + ztpfac3 * dry_static_energy(:)) / zcpq, NCOPIES=ntiles, DIM=2)
    surface_temperature_unfiltered(1:nidx) = &
         (ztpfac2 * dry_static_energy_new(1:nidx) + &
         ztpfac3 * dry_static_energy(1:nidx)) / zcpq(1:nidx)

    ! Correction for snowmelt
!!$ TR    CALL average_tiles(soil%snow(kidx0:kidx1,:), surface%is_present(kidx0:kidx1,:) .AND. .NOT. surface%is_lake(kidx0:kidx1,:) , &
!!$ TR         surface%cover_fract(kidx0:kidx1,:), snow_avg(:))
!!$ TR    WHERE (snow_avg(:) > soil_options%CriticalSnowDepth .OR. ANY(surface%is_glacier(kidx0:kidx1,:),DIM=2)) &
!!$ TR       dry_static_energy_new(:) = (MIN(soil%surface_temperature_unfiltered(kidx0:kidx1,1),tmelt)                &
!!$ TR       * zcpq(:) - ztpfac3 * dry_static_energy(:)) / ztpfac2
    ! Correction for snowmelt
    WHERE (snow(:) > crit_snow_depth) &
       dry_static_energy_new(:) = (MIN(surface_temperature_unfiltered(1:nidx),tmelt)  &
       * zcpq(1:nidx) - ztpfac3 * dry_static_energy(1:nidx)) / ztpfac2

    DO itile=1,ntiles
       ! Compute temperature and moisture at lowest atmospheric level by back-substitution
       air_dry_static_energy_new(:,itile) = t_Acoef(:) * dry_static_energy_new(:) + t_Bcoef(:)
!!$ TR       air_moisture_new(:,itile) = q_Acoef(:) * soil%sat_surface_specific_humidity(kidx0:kidx1,itile) + q_Bcoef(:)
       air_moisture_new(:,itile) = q_Acoef(:) * sat_surface_specific_humidity(:) + q_Bcoef(:)
    END DO

!!$ TR    DO itile=1,ntiles
!!$ TR       soil%radiative_temperature(kidx0:kidx1,itile) = dry_static_energy_new(:) / zcpq(:)
!!$ TR    END DO
    surface_temperature_rad(1:nidx) = dry_static_energy_new(1:nidx) / zcpq(1:nidx)

    ! Compute sensible heat flux
!!$ TR    DO itile=1,ntiles
!!$ TR       soil%sensible_heat_flux(kidx0:kidx1,itile) = zcons30 * cdrag(:) *         &
!!$ TR         (air_dry_static_energy_new(:,itile) - dry_static_energy_new(:) -        &
!!$ TR         SpecificHeatDryAirConstPressure * vtmpc2 * surface_temperature_avg(:) * &
!!$ TR         (cair(:) * air_moisture_new(:,itile) - csat(:) *                        &
!!$ TR         soil%sat_surface_specific_humidity(kidx0:kidx1,itile)))
!!$ TR    END DO

    ! Compute evaporation and latent heat flux
!!$ TR    DO itile=1,ntiles
!!$ TR       ! Evapotranspiration
!!$ TR       soil%evapotranspiration(kidx0:kidx1,itile) = zcons30 * cdrag(:)  &
!!$ TR       * (cair(:) * air_moisture_new(:,itile) - csat(:)  * soil%sat_surface_specific_humidity(kidx0:kidx1,itile))
!!$ TR       ! Transpiration
!!$ TR       soil%transpiration(kidx0:kidx1,itile) = zcons30 * cdrag(:)  &
!!$ TR       * csat_transpiration(:) * (air_moisture_new(:,itile) - soil%sat_surface_specific_humidity(kidx0:kidx1,itile))
!!$ TR       ! Potential evaporation
!!$ TR       soil%evaporation_pot(kidx0:kidx1,itile) = zcons30 * cdrag(:)  &
!!$ TR       * (air_moisture_new(:,itile) - soil%sat_surface_specific_humidity(kidx0:kidx1,itile))
!!$ TR       ! Evaporation over snow
!!$ TR       IF ((nsoil==5) .OR. (ldiag_soil)) &
!!$ TR         snow_evap_tile(1:nidx,itile) = soil%snow_fract(kidx0:kidx1,itile) * soil%evaporation_pot(kidx0:kidx1,itile)
!!$ TR    END DO
    ! Latent heat flux
!!$ TR    soil%latent_heat_flux(kidx0:kidx1,:) = LatentHeatVaporization  * soil%evapotranspiration(kidx0:kidx1,:) &
!!$ TR    + (LatentHeatSublimation - LatentHeatVaporization) * soil%snow_fract(kidx0:kidx1,:) * soil%evaporation_pot(kidx0:kidx1,:)

    DO itile=1,ntiles
       ! Evapotranspiration
       evapotranspiration(1:nidx) = zcons30 * cdrag(1:nidx) * &
       (cair(1:nidx) * air_moisture_new(1:nidx,itile) - &
       csat(1:nidx)  * sat_surface_specific_humidity(1:nidx))
       ! Transpiration
       transpiration(1:nidx) = zcons30 * cdrag(1:nidx) * &
       csat_transpiration(1:nidx) * (air_moisture_new(1:nidx,itile) - &
       sat_surface_specific_humidity(1:nidx))
       ! Potential evaporation
       evaporation_pot(1:nidx,itile) = zcons30 * cdrag(1:nidx) * &
       (air_moisture_new(1:nidx,itile) - sat_surface_specific_humidity(1:nidx))
    END DO

!   Initial values for 5 layer scheme
!!$ TR    IF (nsoil==5) THEN
!   *** Note that the atmosphere does not see tiles, thus the evaporative fluxes are
!          valid for the gridbox average. Therefore, separating the water fluxes per tile is 
!          leading to errors in the water balance, especially for the 5 layer hydrology scheme.
!!$ TR      CALL average_tiles(wet_skin_fract(1:nidx, :),                                       &
!!$ TR           surface%is_present(kidx0:kidx1,:) .AND. .NOT. surface%is_lake(kidx0:kidx1,:) , &
!!$ TR           surface%cover_fract(kidx0:kidx1,:), WetSkin_frac(1:nidx) )
!!$ TR      CALL average_tiles(skin_reservoir_max(1:nidx, :),                                   &
!!$ TR           surface%is_present(kidx0:kidx1,:) .AND. .NOT. surface%is_lake(kidx0:kidx1,:) , &
!!$ TR           surface%cover_fract(kidx0:kidx1,:), MaxSkinResCap(1:nidx) )
!!$ TR    ENDIF

!!$ TR    IF ((nsoil==5) .OR. (ldiag_soil)) THEN
    IF (nsoil == 5) THEN
!!$ TR      CALL average_tiles(soil%skin_reservoir(kidx0:kidx1,:),                              &
!!$ TR           surface%is_present(kidx0:kidx1,:) .AND. .NOT. surface%is_lake(kidx0:kidx1,:) , &
!!$ TR           surface%cover_fract(kidx0:kidx1,:), SkinReservoire(1:nidx) )
       SkinReservoire(1:nidx) = skin_reservoir(1:nidx)
!
!   *** Save sum of water storages from previous time step
!!$ TR      CALL average_tiles(canopy_snow(1:nidx, :),                                          &
!!$ TR           surface%is_present(kidx0:kidx1,:) .AND. .NOT. surface%is_lake(kidx0:kidx1,:) , &
!!$ TR           surface%cover_fract(kidx0:kidx1,:), SnowDepthCan(1:nidx) )
!!$ TR      IF (nsoil == 5) THEN
!!$ TR        DO i=kidx0, kidx1
!!$ TR          soil_water_pre(i - kidx0 + 1) = SUM( soil%layer_moisture(i,1:soil%nsoil) )
!!$ TR        ENDDO
        DO itile=1,ntiles
           SkinReservoire_tile(1:nidx, itile) = SkinReservoire(1:nidx)
        ENDDO
!!$ TR      ELSE
!!$ TR        CALL average_tiles(soil%moisture(kidx0:kidx1,1,:), surface%is_present(kidx0:kidx1,:) &
!!$ TR             .AND. .NOT. surface%is_lake(kidx0:kidx1,:), surface%cover_fract(kidx0:kidx1,:), &
!!$ TR             soil_water_pre(1:nidx)  )
!!$ TR      ENDIF
!!$ TR      soil_water_pre(1:nidx) = soil_water_pre(1:nidx) + SkinReservoire(1:nidx) + SnowDepthCan(1:nidx) + snow_avg(1:nidx)
    ENDIF
    !
    !--------------------------------------------------------------------------------------------------------
    ! 
    tte_corr(:,:)                  = 0._wp
    glacier_precip_minus_evap(:,:) = 0._wp
    melt_water_excess(:)           = 0._wp
    surface_runoff(:,:)            = 0._wp
    drainage(:,:)                  = 0._wp
    glacier_depth(:,:)             = 0._wp
    is_present(:,:) = .TRUE.
    is_glacier(:,:) = .FALSE.

    DO itile=1,ntiles

!!$ TR       IF (nsoil==5) SkinReservoire(1:nidx) = soil%skin_reservoir(kidx0:kidx1,itile)
       IF (nsoil==5) SkinReservoire(1:nidx) = skin_reservoir(1:nidx)

!!$ TR       CALL update_surf_down( nidx,                                                  &
!!$ TR            air_moisture(1:nidx),                                                    &
!!$ TR            soil%surface_temperature_unfiltered(kidx0:kidx1,itile), wind10(1:nidx),  &
!!$ TR            air_temperature(1:nidx), soil%skin_reservoir(kidx0:kidx1,itile),         &
!!$ TR            wet_skin_fract(1:nidx,itile), skin_reservoir_max(1:nidx,itile),          &
!!$ TR            soil%snow(kidx0:kidx1,itile), soil%snow_fract(kidx0:kidx1,itile),        &
!!$ TR            canopy_snow(1:nidx,itile), soil%glacier_depth(kidx0:kidx1,itile),        &
!!$ TR            soil%heat_capacity(kidx0:kidx1,itile),                                   &
!!$ TR            soil%evapotranspiration(kidx0:kidx1,itile),                              &
!!$ TR            soil%evaporation_pot(kidx0:kidx1,itile),                                 &
!!$ TR            evapotranspiration_no_snow_skin(1:nidx),                                 &
!!$ TR            rain, snow,                                                              &
!!$ TR            surface%is_present(kidx0:kidx1,itile),                                   &
!!$ TR            surface%is_glacier(kidx0:kidx1,itile),                                   &
!!$ TR            glacier_precip_minus_evap(:,itile),                                      &
!!$ TR            soil%snow_acc(kidx0:kidx1,itile), soil%snow_melt_acc(kidx0:kidx1,itile), &
!!$ TR            soil%glacier_runoff_acc(kidx0:kidx1,itile), melt_water_excess(1:nidx),   &
!!$ TR            tte_corr(1:nidx,itile),soil%snow_melt(kidx0:kidx1,itile))

       CALL update_surf_down( nidx, delta_time, time_steps_soil(1:nidx),        &
            air_moisture(1:nidx),                                               &
            surface_temperature_unfiltered(1:nidx), wind10(1:nidx),             &
            air_temperature(1:nidx), skin_reservoir(1:nidx),                    &
            wet_skin_fract(1:nidx,itile), skin_reservoir_max(1:nidx,itile),     &
            snow(1:nidx), snow_fract(1:nidx),                                   &
            snow_canopy(1:nidx), glacier_depth(1:nidx,itile),                   &
            heat_capacity(1:nidx),                                              &
            evapotranspiration(1:nidx),                                         &
            evaporation_pot(1:nidx,itile),                                      &
            evapotranspiration_no_snow_skin(1:nidx),                            &
            precip_rain(1:nidx), precip_snow(1:nidx),                           &
            is_present(1:nidx,itile),                                           &
            is_glacier(1:nidx,itile),                                           &
            glacier_precip_minus_evap(1:nidx,itile),                            &
            snow_acc(1:nidx), snow_melt_acc(1:nidx),                            &
            glacier_runoff_acc(1:nidx), melt_water_excess(1:nidx),              &
            tte_corr(1:nidx,itile), snow_melt(1:nidx))

!!$ TR       CALL update_soiltemp(nidx,soil%ntsoil                                                                   &
!!$ TR                         , soil%surface_temperature_unfiltered(kidx0:kidx1,itile),soil%snow(kidx0:kidx1,itile) &
!!$ TR                         , soil_param%ThermalDiffusivity(kidx0:kidx1)                                          &
!!$ TR                         , soil_param%VolHeatCapacity(kidx0:kidx1)                                             &
!!$ TR                         , soil%c_soil_temperature(kidx0:kidx1,:,itile)                                        &
!!$ TR                         , soil%d_soil_temperature(kidx0:kidx1,:,itile)                                        &
!!$ TR                         , soil%soil_temperature(kidx0:kidx1,:,itile)                                          &
!!$ TR                         , soil%heat_capacity(kidx0:kidx1,itile),soil%ground_heat_flux(kidx0:kidx1,itile)      &
!!$ TR                         , surface%is_present(kidx0:kidx1,itile)                                               &
!!$ TR                         , surface%is_glacier(kidx0:kidx1,itile) )

       !! for JSBACH testing
       c_soil_temperature(1:nidx,1) = c_soil_temperature1(1:nidx)
       c_soil_temperature(1:nidx,2) = c_soil_temperature2(1:nidx)
       c_soil_temperature(1:nidx,3) = c_soil_temperature3(1:nidx)
       c_soil_temperature(1:nidx,4) = c_soil_temperature4(1:nidx)
       c_soil_temperature(1:nidx,5) = c_soil_temperature5(1:nidx)
       d_soil_temperature(1:nidx,1) = d_soil_temperature1(1:nidx)
       d_soil_temperature(1:nidx,2) = d_soil_temperature2(1:nidx)
       d_soil_temperature(1:nidx,3) = d_soil_temperature3(1:nidx)
       d_soil_temperature(1:nidx,4) = d_soil_temperature4(1:nidx)
       d_soil_temperature(1:nidx,5) = d_soil_temperature5(1:nidx)
       soil_temperature(1:nidx,1) = soil_temperature1(1:nidx)
       soil_temperature(1:nidx,2) = soil_temperature2(1:nidx)
       soil_temperature(1:nidx,3) = soil_temperature3(1:nidx)
       soil_temperature(1:nidx,4) = soil_temperature4(1:nidx)
       soil_temperature(1:nidx,5) = soil_temperature5(1:nidx)

       CALL update_soiltemp_icon(nidx,ntsoil                        &
                      , delta_time                                  &
                      , time_steps_soil(1:nidx)                     &
                      , surface_temperature_unfiltered(1:nidx)      &
                      , ThermalDiffusivity(1:nidx)                  &
                      , VolHeatCapacity(1:nidx)                     &
                      , c_soil_temperature(1:nidx,1:ntsoil)         &
                      , d_soil_temperature(1:nidx,1:ntsoil)         &
                      , soil_temperature(1:nidx,1:ntsoil)           &
                      , heat_capacity(1:nidx)                       &
                      , ground_heat_flux(1:nidx)                    &
                      )

!
!      update_surf_up is called for each tile, 5 layer hydrology routine is without tiles
!               as this is incompatible with REMO
!               --> Dummy Separation into tiles before update_surf_up and tile averaging
!                   afterwards using new temporary field zlayer_moisture_tile
       IF (nsoil == 5) THEN
!!$ TR         layer_moisture_tile(1:nidx, :, itile) = soil%layer_moisture(kidx0:kidx1, :)
!!$ TR         SkinReservoire_tile(1:nidx,    itile) = SkinReservoire_tile(1:nidx, itile) + &
!!$ TR                                                 soil%skin_reservoir(kidx0:kidx1,itile) - SkinReservoire(1:nidx)
         layer_moisture_tile(1:nidx, :, itile) = soil_layer_moisture(1:nidx, :)
         SkinReservoire_tile(1:nidx,    itile) = SkinReservoire_tile(1:nidx, itile) + &
                                                 skin_reservoir(1:nidx) - SkinReservoire(1:nidx)

          ! Note: surface_runoff and drainage have to be passed to the runoff model once this is part of jsbach

!!$ TR         CALL update_surf_up(nidx, itile, nsoil,                                          &
!!$ TR            SkinReservoire_tile(1:nidx, itile), WetSkin_frac(1:nidx), MaxSkinResCap(1:nidx),    &
!!$ TR            soil%moisture(kidx0:kidx1,1,itile), soil_param%MaxMoisture(kidx0:kidx1,1),                    &
!!$ TR            soil%soil_temperature(kidx0:kidx1,1,itile), soil%snow_fract(kidx0:kidx1,itile),               &
!!$ TR            surface%oro_std_dev(kidx0:kidx1),                                                             &
!!$ TR            soil%evapotranspiration(kidx0:kidx1,itile), soil%evaporation_pot(kidx0:kidx1,itile),          &
!!$ TR            evapotranspiration_no_snow_skin(1:nidx),                                                      &
!!$ TR            rain(1:nidx),                                                                                 &
!!$ TR            surface%is_present(kidx0:kidx1,itile),                                                        &
!!$ TR            surface%is_glacier(kidx0:kidx1,itile),                                                        &
!!$ TR            soil%runoff_acc(kidx0:kidx1,itile), soil%drainage_acc(kidx0:kidx1,itile),                     &
!!$ TR            surface_runoff(1:nidx,itile), drainage(1:nidx,itile), melt_water_excess(1:nidx),              &
!!$ TR            skin_res_evap_tile(1:nidx,itile),                                                             &
!!$ TR            ! 5 layer hydrology fields (optional)      
!!$ TR            cdel,                                  &
!!$ TR            layer_moisture_tile(1:nidx, :, itile), &            
!!$ TR            FieldCapacity(1:nidx, :),              &            
!!$ TR            soil_param%RootDepth(kidx0:kidx1),     &            
!!$ TR            soil_param%SoilDepth(kidx0:kidx1),     &            
!!$ TR            soil_param%hyd_cond_sat(kidx0:kidx1),  &            
!!$ TR            soil_param%PotMoisture(kidx0:kidx1),   &            
!!$ TR            soil_param%bclapp(kidx0:kidx1),        &            
!!$ TR            soil_param%SoilPorosity(kidx0:kidx1),  &            
!!$ TR            soil%transpiration(kidx0:kidx1,itile), &
!!$ TR            reduced_evap(1:nidx,itile),            &
!!$ TR            soil_param%FieldCapacity(kidx0:kidx1), &            
!!$ TR            soil_param%WiltingPoint(kidx0:kidx1),  &            
!!$ TR            soil_param%PoreSizeIndex(kidx0:kidx1), &            
!!$ TR            jllog                                  &
!!$ TR           )
!!$ TR global values for soil properties to test JSBACH
         RootDepth(1:nidx) = 2.65355_wp
         SoilDepth(1:nidx) = 2.72603_wp
         hyd_cond_sat(1:nidx) = 4.60588e-06_wp
         PotMoisture(1:nidx) = 0.275715_wp
         bclapp(1:nidx) = 8.58824_wp
         SoilPorosity(1:nidx) = 0.441515_wp
         FieldCapacity_param(1:nidx) = 0.318517_wp
         WiltingPoint = 0.185245_wp
         PoreSizeIndex = 0.20747_wp

         CALL update_surf_up(nidx, itile, delta_time, nsoil,                                 &
            SkinReservoire_tile(1:nidx, itile),                                              &
            wet_skin_fract(1:nidx,itile), skin_reservoir_max(1:nidx,itile),                  &
            soil_moisture(1:nidx,1,itile), soil_MaxMoisture(1:nidx,1),                       &
            soil_temperature(1:nidx,1), snow_fract(1:nidx),                                  &
            oro_std_dev(1:nidx),                                                             &
            evapotranspiration(1:nidx), evaporation_pot(1:nidx,itile),                       &
            evapotranspiration_no_snow_skin(1:nidx),                                         &
            precip_rain(1:nidx),                                                             &
            is_present(1:nidx,itile),                                                        &
            is_glacier(1:nidx,itile),                                                        &
            runoff_acc(1:nidx), drainage_acc(1:nidx),                                        &
            surface_runoff(1:nidx,itile), drainage(1:nidx,itile), melt_water_excess(1:nidx), &
            skin_res_evap_tile(1:nidx,itile),                                                &
            ! 5 layer hydrology fields (optional)      
            cdel(1:nsoil),                         &
            layer_moisture_tile(1:nidx, :, itile), &            
            FieldCapacity(1:nidx,1:nsoil),         &            
            RootDepth(1:nidx),                     &            
            SoilDepth(1:nidx),                     &            
            hyd_cond_sat(1:nidx),                  &            
            PotMoisture(1:nidx),                   &            
            bclapp(1:nidx),                        &            
            SoilPorosity(1:nidx),                  &            
            transpiration(1:nidx),                 &
            reduced_evap(1:nidx,itile),            &
            FieldCapacity_param(1:nidx),           &            
            WiltingPoint(1:nidx),                  &            
            PoreSizeIndex(1:nidx),                 &            
            jllog                                  &
           )

!!$ TR       ELSE ! 1 layer soil model below
!!$ TR
!!$ TR         CALL update_surf_up(nidx, itile, nsoil,                                                          &
!!$ TR            soil%skin_reservoir(kidx0:kidx1,itile), wet_skin_fract(:,itile), skin_reservoir_max(:,itile), &
!!$ TR            soil%moisture(kidx0:kidx1,1,itile), soil_param%MaxMoisture(kidx0:kidx1,1),                    &
!!$ TR            soil%soil_temperature(kidx0:kidx1,1,itile), soil%snow_fract(kidx0:kidx1,itile),               &
!!$ TR            surface%oro_std_dev(kidx0:kidx1),                                                             &
!!$ TR            soil%evapotranspiration(kidx0:kidx1,itile), soil%evaporation_pot(kidx0:kidx1,itile),          &
!!$ TR            evapotranspiration_no_snow_skin(1:nidx),                                                      &
!!$ TR            rain(1:nidx),                                                                                 &
!!$ TR            surface%is_present(kidx0:kidx1,itile),                                                        &
!!$ TR            surface%is_glacier(kidx0:kidx1,itile),                                                        &
!!$ TR            soil%runoff_acc(kidx0:kidx1,itile), soil%drainage_acc(kidx0:kidx1,itile),                     &
!!$ TR            surface_runoff(1:nidx,itile), drainage(1:nidx,itile), melt_water_excess(1:nidx),              &
!!$ TR            skin_res_evap_tile(1:nidx,itile)                                                              &
!!$ TR           )
       ENDIF

!!$ TR       IF ((nsoil == 5) .OR. (ldiag_soil)) RedEvapotransp(1:nidx,itile) = evapotranspiration_no_snow_skin(1:nidx)
       IF (nsoil == 5) RedEvapotransp(1:nidx,itile) = &
          evapotranspiration_no_snow_skin(1:nidx)

    END DO                  ! End of loop over tiles

!!$ TR    CALL average_tiles(tte_corr, surface%is_present(kidx0:kidx1,:) .AND. .NOT. surface%is_lake(kidx0:kidx1,:) , &
!!$ TR         surface%cover_fract(kidx0:kidx1,:), tte_corr_avg)
!!$ TR    CALL average_tiles(soil%moisture(kidx0:kidx1,1,:), surface%is_present(kidx0:kidx1,:)              &
!!$ TR                      .AND. .NOT. surface%is_lake(kidx0:kidx1,:), surface%cover_fract(kidx0:kidx1,:), &
!!$ TR                      soil_moisture_avg)
!!$ TR    CALL average_tiles(soil%snow(kidx0:kidx1,:), surface%is_present(kidx0:kidx1,:)                    &
!!$ TR                      .AND. .NOT. surface%is_lake(kidx0:kidx1,:), surface%cover_fract(kidx0:kidx1,:), &
!!$ TR                      snow_avg)
!!$ TR    CALL average_tiles(soil%surface_temperature_unfiltered(kidx0:kidx1,:), surface%is_present(kidx0:kidx1,:) &
!!$ TR                      .AND. .NOT. surface%is_lake(kidx0:kidx1,:), surface%cover_fract(kidx0:kidx1,:),        &
!!$ TR                      surface_temperature_avg)
!!$ TR    IF ((nsoil==5) .OR. (ldiag_soil)) THEN
!!$ TR      CALL average_tiles(skin_res_evap_tile(1:nidx,:), surface%is_present(kidx0:kidx1,:)                &
!!$ TR                        .AND. .NOT. surface%is_lake(kidx0:kidx1,:), surface%cover_fract(kidx0:kidx1,:), &
!!$ TR                        skin_reservoire_evap)
!!$ TR      CALL average_tiles(RedEvapotransp(1:nidx,:), surface%is_present(kidx0:kidx1,:)                    &
!!$ TR                        .AND. .NOT. surface%is_lake(kidx0:kidx1,:), surface%cover_fract(kidx0:kidx1,:), &
!!$ TR                        evapotranspiration_no_snow_skin)
!!$ TR      CALL average_tiles(snow_evap_tile(1:nidx,:), surface%is_present(kidx0:kidx1,:)                    &
!!$ TR                        .AND. .NOT. surface%is_lake(kidx0:kidx1,:), surface%cover_fract(kidx0:kidx1,:), &
!!$ TR                        snow_evap)
!!$ TR    ENDIF
    IF (nsoil==5) THEN
!!$ TR      DO i=1, soil%nsoil
!!$ TR        CALL average_tiles(layer_moisture_tile(1:nidx,i,:), surface%is_present(kidx0:kidx1,:)           &
!!$ TR                        .AND. .NOT. surface%is_lake(kidx0:kidx1,:), surface%cover_fract(kidx0:kidx1,:), &
!!$ TR                        soil%layer_moisture(kidx0:kidx1,i) )
!!$ TR      END DO
      moisture1(1:nidx) = layer_moisture_tile(1:nidx,1,1)
      moisture2(1:nidx) = layer_moisture_tile(1:nidx,2,1)
      moisture3(1:nidx) = layer_moisture_tile(1:nidx,3,1)
      moisture4(1:nidx) = layer_moisture_tile(1:nidx,4,1)
      moisture5(1:nidx) = layer_moisture_tile(1:nidx,5,1)
      moisture_all(1:nidx) = soil_moisture(1:nidx,1,1)

      ! update the skin_reservoir per tile in a consistent way
!!$ TR      CALL average_tiles(SkinReservoire_tile(1:nidx, :),                                  &
!!$ TR           surface%is_present(kidx0:kidx1,:) .AND. .NOT. surface%is_lake(kidx0:kidx1,:) , &
!!$ TR           surface%cover_fract(kidx0:kidx1,:), SkinReservoire(1:nidx) )
      DO itile=1,ntiles
        DO i=1, nidx
!!$ TR          j=kidx0+i-1
!!$ TR          IF (MaxSkinResCap(i) > EPSILON(1._dp)) THEN
!!$ TR            soil%skin_reservoir(j, itile) = MIN(SkinReservoire(i)/MaxSkinResCap(i), 1._dp) * skin_reservoir_max(i, itile)
!!$ TR          ELSE    
!!$ TR            soil%skin_reservoir(j, itile) = 0._dp
          IF (skin_reservoir_max(i,itile) > EPSILON(1._wp)) THEN
            skin_reservoir(i) = MIN(SkinReservoire_tile(i,itile),skin_reservoir_max(i,itile))
          ELSE
            skin_reservoir(i) = 0._wp
          ENDIF
        ENDDO
      ENDDO
    ENDIF

    ! Test that soil moisture has still a positive value
!!$ TR    IF (MINVAL(soil_moisture_avg(:)) < 0._dp .AND. flag_soil_moisture) THEN
!!$ TR      CALL message('update_soil','Warning: Negative soil moisture occurred!')
!!$ TR      flag_soil_moisture = .false.
!!$ TR    END IF
!!$ TR    soil_moisture_avg(:) = MAX(EPSILON(1._dp),soil_moisture_avg(:))
!!$ TR    DO itile=1,ntiles
!!$ TR      soil%moisture(kidx0:kidx1,1,itile) = soil_moisture_avg(:)
!!$ TR      soil%snow(kidx0:kidx1,itile) = snow_avg(:)
!!$ TR      soil%surface_temperature_unfiltered(kidx0:kidx1,itile) = surface_temperature_avg(:)
!!$ TR      soil%soil_temperature(kidx0:kidx1,1,itile) = surface_temperature_avg(:)
    soil_temperature(1:nidx,1) = surface_temperature_unfiltered(1:nidx)
!!$ TR    END DO
!!$ TR    DO isoil=2,soil%ntsoil
!!$ TR       CALL average_tiles(soil%soil_temperature(kidx0:kidx1,isoil,:), surface%is_present(kidx0:kidx1,:) &
!!$ TR                    .AND. .NOT. surface%is_lake(kidx0:kidx1,:), surface%cover_fract(kidx0:kidx1,:), &
!!$ TR                    soil_temperature_avg)
!!$ TR       DO itile=1,ntiles
!!$ TR          soil%soil_temperature(kidx0:kidx1,isoil,itile) = soil_temperature_avg
!!$ TR       END DO
!!$ TR    END DO

    !! for JSBACH testing
    c_soil_temperature1(1:nidx) = c_soil_temperature(1:nidx,1)
    c_soil_temperature2(1:nidx) = c_soil_temperature(1:nidx,2)
    c_soil_temperature3(1:nidx) = c_soil_temperature(1:nidx,3)
    c_soil_temperature4(1:nidx) = c_soil_temperature(1:nidx,4)
    c_soil_temperature5(1:nidx) = c_soil_temperature(1:nidx,5)
    d_soil_temperature1(1:nidx) = d_soil_temperature(1:nidx,1)
    d_soil_temperature2(1:nidx) = d_soil_temperature(1:nidx,2)
    d_soil_temperature3(1:nidx) = d_soil_temperature(1:nidx,3)
    d_soil_temperature4(1:nidx) = d_soil_temperature(1:nidx,4)
    d_soil_temperature5(1:nidx) = d_soil_temperature(1:nidx,5)
    soil_temperature1(1:nidx) = soil_temperature(1:nidx,1)
    soil_temperature2(1:nidx) = soil_temperature(1:nidx,2)
    soil_temperature3(1:nidx) = soil_temperature(1:nidx,3)
    soil_temperature4(1:nidx) = soil_temperature(1:nidx,4)
    soil_temperature5(1:nidx) = soil_temperature(1:nidx,5)

        ! Time filter for surface temperature
!!$ TR    IF (lsurf) THEN
!!$ TR       IF (.NOT. lstart) THEN  !! lstart ???
!!$ TR          WHERE (surface%is_present(kidx0:kidx1,:))
         IF (time_steps_soil(1) > 0.5_wp) THEN
             surface_temperature(1:nidx) = surface_temperature_old(1:nidx) +     &
                  eps * (surface_temperature(1:nidx)                             &
                  - 2._wp * surface_temperature_old(1:nidx)                      &
                  + surface_temperature_unfiltered(1:nidx))             
             surface_temperature_old(1:nidx) = surface_temperature_unfiltered(1:nidx)
!!$ TR          END WHERE
         ELSE
            surface_temperature(1:nidx) = surface_temperature_old(1:nidx)
         END IF
!!$ TR    ELSE
!!$ TR       soil%surface_temperature_old(kidx0:kidx1,:) = soil%surface_temperature(kidx0:kidx1,:)
!!$ TR       soil%surface_temperature_unfiltered(kidx0:kidx1,:) = soil%surface_temperature(kidx0:kidx1,:)
!!$ TR    END IF

    !--------------------------------------------------------------------------------------------------------
    ! Re-calculate qsat,qair and zhsoil with updated moisture parameters due to time mismatch with old echam
    !--------------------------------------------------------------------------------------------------------

    DO  itile=1,ntiles
      DO  i=1,nidx
!!$ TR       j=kidx0+i-1
!!$ TR        IF (soil_mask(i,itile)) THEN
!!$ TR           wet_skin_fract(i,itile) = MERGE(MIN(1._dp, soil%skin_reservoir(j,itile) / skin_reservoir_max(i,itile)), &
!!$ TR                0.0_dp, soil_options%SkinReservoirMax > EPSILON(1._dp))
           wet_skin_fract(i,itile) = MIN(1._wp, skin_reservoir(i) / skin_reservoir_max(i,itile))
          
           ! Roesch et al, 2002, Climate Dynamics
!!$ TR           soil%snow_fract(j,itile) = zqsncr * TANH(soil%snow(j,itile) * 100._dp) *            &
!!$ TR                SQRT(soil%snow(j,itile) * 1000._dp / (soil%snow(j,itile) * 1000._dp + zepsec + &
!!$ TR                zsigfac * surface%oro_std_dev(j)))
           snow_fract(i) = zqsncr * TANH(snow(i) * 100._wp) *            &
                SQRT(snow(i) * 1000._wp / (snow(i) * 1000._wp + zepsec + &
                zsigfac * oro_std_dev(i)))           

!!$ TR           IF (soil%snow_fract(j,itile) < EPSILON(1._dp) .AND. canopy_snow_fract(i,itile) > EPSILON(1._dp)) &
!!$ TR                soil%snow_fract(j,itile) = canopy_snow_fract(i,itile)
           IF (snow_fract(i) < EPSILON(1._wp) .AND. canopy_snow_fract(i,itile) > EPSILON(1._wp)) &
                snow_fract(i) = canopy_snow_fract(i,itile)
           ! Modify snow cover if snow loss during the time step due to potential evaporation is larger than 
           ! equivalent snow water content from soil and canopy; same for skin reservoir
           ! Potential evaporation using old values of air and surface humidity
!!$ TR           evaporation_pot(i,itile) = zcons30 * cdrag(i) * &
!!$ TR                (sat_surface_specific_hum_old(i,itile) - air_moisture(i))
           evaporation_pot(i,itile) = zcons30 * cdrag(i) * &
                (sat_surface_specific_hum_old(i) - air_moisture(i))
!!$ TR      IF (soil%snow_fract(j,itile) > 0._dp .AND. soil%snow(j,itile) + canopy_snow(i,itile) > 0._dp) THEN
!!$ TR              soil%snow_fract(j,itile) = soil%snow_fract(j,itile) /                              &
!!$ TR                   MAX(1._dp, soil%snow_fract(j,itile) * evaporation_pot(i,itile) * delta_time / &
!!$ TR                   (RhoH2O*(soil%snow(j,itile) + canopy_snow(i,itile))))
!!$ TR           END IF

           IF (snow_fract(i) > 0._wp .AND. snow(i) + snow_canopy(i) > 0._wp) THEN
              snow_fract(i) = snow_fract(i) / &
                   MAX(1._wp, snow_fract(i) * evaporation_pot(i,itile) * delta_time / &
                   (RhoH2O * (snow(i) + snow_canopy(i))))
           END IF
!!$ TR           IF (wet_skin_fract(i,itile) > 0._dp) THEN
!!$ TR              wet_skin_fract(i,itile) = wet_skin_fract(i,itile) /                                          &
!!$ TR                   MAX(1._dp, (1._dp - soil%snow_fract(j,itile)) * evaporation_pot(i,itile) * delta_time / &
!!$ TR                   (RhoH2O * MAX(EPSILON(1._dp),soil%skin_reservoir(j,itile))))
!!$ TR           END IF
           IF (wet_skin_fract(i,itile) > 0._wp) THEN
              wet_skin_fract(i,itile) = wet_skin_fract(i,itile) / &
                   MAX(1._wp, (1._wp - snow_fract(i)) * evaporation_pot(i,itile) * delta_time / &
                   (RhoH2O * MAX(EPSILON(1._wp),skin_reservoir(i))))
           END IF
           
!!$ TR        ELSE IF (surface%is_glacier(j,itile)) THEN
!!$ TR           wet_skin_fract(i,itile)  = 0.0_dp
!!$ TR           soil%snow_fract(j,itile) = 1.0_dp
!!$ TR        ELSE
!!$ TR           wet_skin_fract(i,itile)  = 0.0_dp
!!$ TR           soil%snow_fract(j,itile) = 0.0_dp
!!$ TR        END IF
      END DO
    END DO

!!$    relative_humidity = 0.0_dp
!!$
!!$    ! Calculate relative humidity from water content in first soil layer
!!$#ifdef __PGI
!!$    moisture_max = SPREAD(soil_param%MaxMoisture(kidx0:kidx1,1), NCOPIES=ntiles, DIM=2)
!!$    relative_humidity = &
!!$         calc_relative_humidity(soil%moisture(kidx0:kidx1,1,:), &
!!$         moisture_max, &
!!$         soil_options%MoistureFractWilting)
!!$#else
!!$    relative_humidity = &
!!$         calc_relative_humidity(soil%moisture(kidx0:kidx1,1,:), &
!!$         SPREAD(soil_param%MaxMoisture(kidx0:kidx1,1), NCOPIES=ntiles, DIM=2), &
!!$         soil_options%MoistureFractWilting)
!!$#endif
!!$
!!$    !------------------------------------------------------------------------------------------------------
!!$
    qsat_fact(:,:) = 0._wp
    qair_fact(:,:) = 0._wp
    qsat_veg(:,:) = 0._wp
    qair_veg(:,:) = 0._wp
!!$    qsat_transpiration(:,:) = 0._dp
!!$    zhsoil = 0._dp
!!$
!!$    DO itile=1,ntiles
!!$       DO i=1,nidx
!!$          IF (soil%moisture(j,1,itile) > & 
!!$              (soil_options%MoistureFractWilting * soil_param%MaxMoisture(j,1))) THEN
!!$             qsat_veg(i,itile) = soil%snow_fract(j,itile) + (1._dp - soil%snow_fract(j,itile)) *   &
!!$                                 (wet_skin_fract(i,itile) + (1._dp - wet_skin_fract(i,itile)) /    &
!!$                                 (1._dp + p_echam_zchl(i) * canopy_resistance(i,itile) *           &
!!$                                 MAX(1.0_dp,windspeed(i))))
!!$          qsat_transpiration(i,itile) = (1._dp - soil%snow_fract(j,itile)) *                    &
!!$                                        (1._dp - wet_skin_fract(i,itile)) /                     &
!!$                                        (1._dp + p_echam_zchl(i) * canopy_resistance(i,itile) * &
!!$                                        MAX(1.0_dp,windspeed(i)))
!!$          ELSE
!!$             qsat_veg(i,itile) = soil%snow_fract(j,itile) + (1._dp - soil%snow_fract(j,itile)) * &
!!$                                 wet_skin_fract(i,itile)
!!$             qsat_transpiration(i,itile) = 0._dp
!!$          END IF
!!$          qair_veg(i,itile) = qsat_veg(i,itile)
!!$       END DO
!!$    END DO
!!$        
!!$    WHERE (relative_humidity(1:nidx,:) > SPREAD(air_moisture(1:nidx), NCOPIES=ntiles, DIM=2) / &
!!$         SPREAD(sat_surf_specific_hum_avg(1:nidx), NCOPIES=ntiles, DIM=2) .AND. &
!!$         relative_humidity(1:nidx,:) > 1.e-10_dp )
!!$
!!$       qsat_fact(:,:) = soil%snow_fract(kidx0:kidx1,:) +                                 & 
!!$                        (1._dp - soil%snow_fract(kidx0:kidx1,:)) *                       &
!!$                        (wet_skin_fract(1:nidx,:) + (1._dp - wet_skin_fract(1:nidx,:)) * &
!!$                        relative_humidity(1:nidx,:))
!!$       qair_fact(:,:) = 1._dp
!!$
!!$    ELSEWHERE (surface%is_present(kidx0:kidx1,:))
!!$
!!$       qsat_fact(:,:) = soil%snow_fract(kidx0:kidx1,:) + &
!!$                        (1._dp - soil%snow_fract(kidx0:kidx1,:)) * wet_skin_fract(1:nidx,:)
!!$       qair_fact(:,:) = qsat_fact(:,:)
!!$
!!$    END WHERE
!!$
!!$    WHERE(SPREAD(air_moisture(1:nidx), NCOPIES=ntiles, DIM=2) > &
!!$         SPREAD(sat_surf_specific_hum_avg(1:nidx), NCOPIES=ntiles, DIM=2))
!!$       qsat_fact(:,:) = 1._dp
!!$       qair_fact(:,:) = 1._dp
!!$    END WHERE
!!$  
!!$    csat_tiles(1:nidx,:)=Surface%veg_ratio(kidx0:kidx1,:) * qsat_veg(1:nidx,:) +            &
!!$                         (1._dp - Surface%veg_ratio(kidx0:kidx1,:)) * qsat_fact(1:nidx,:)
!!$    cair_tiles(1:nidx,:)=Surface%veg_ratio(kidx0:kidx1,:) * qair_veg(1:nidx,:) +            &
!!$                         (1._dp - Surface%veg_ratio(kidx0:kidx1,:)) * qair_fact(1:nidx,:)
!!$    csat_transpiration_tiles(1:nidx,:) = Surface%veg_ratio(kidx0:kidx1,:) * qsat_transpiration(1:nidx,:)
!!$
!!$    ! ECHAM5 compatibility: one surface temperature for whole grid box
!!$
!!$    CALL average_tiles(csat_tiles(1:nidx,:), surface%is_present(kidx0:kidx1,:) .AND. .NOT.                &
!!$                       surface%is_lake(kidx0:kidx1,:), surface%cover_fract(kidx0:kidx1,:), csat(1:nidx))
!!$    CALL average_tiles(cair_tiles(1:nidx,:), surface%is_present(kidx0:kidx1,:) .AND. .NOT.                &
!!$                       surface%is_lake(kidx0:kidx1,:), surface%cover_fract(kidx0:kidx1,:), cair(1:nidx))
!!$    CALL average_tiles(csat_transpiration_tiles(1:nidx,:), surface%is_present(kidx0:kidx1,:) .AND. .NOT.  &
!!$                       surface%is_lake(kidx0:kidx1,:), surface%cover_fract(kidx0:kidx1,:), csat_transpiration(1:nidx))
!!$
!!$    WHERE (SPREAD(air_moisture(1:nidx), NCOPIES=ntiles, DIM=2) < sat_surface_specific_hum_old(1:nidx,:))
!!$       zhsoil = soil%snow_fract(kidx0:kidx1,:) +                                 & 
!!$                (1._dp - soil%snow_fract(kidx0:kidx1,:)) *                       &
!!$                (wet_skin_fract(1:nidx,:) + (1._dp - wet_skin_fract(1:nidx,:)) * &
!!$                relative_humidity(1:nidx,:))
!!$    ELSEWHERE
!!$       zhsoil = 1._dp
!!$    END WHERE
!!$    CALL average_tiles(zhsoil, surface%is_present(kidx0:kidx1,:) .AND. .NOT. surface%is_lake(kidx0:kidx1,:) , &
!!$         surface%cover_fract(kidx0:kidx1,:), zhsoil_avg)
!!$
!!$    soil%csat(kidx0:kidx1) = csat(1:nidx)
!!$    soil%cair(kidx0:kidx1) = cair(1:nidx)
!!$    soil%csat_transpiration(kidx0:kidx1) = csat_transpiration(1:nidx)
    csat(1:nidx) = 1._wp
    cair(1:nidx) = 1._wp
    csat_transpiration(1:nidx) = 1._wp
    !------------------------------------------------------------------------------------------------------------
    ! END RECALC QSAT QAir and zhsoil
    !------------------------------------------------------------------------------------------------------------
    ! Note: glacier_precip_minus_evap is in m water equivalent; to convert to kg/m^2s multiply by RhoH2O/delta_time,
    !       and for accumulation this is again multiplied by delta_time, i.e. just multiply by RhoH2O
    ! Accumulate P-E for glaciers
!!$ TR    WHERE (surface%is_present(kidx0:kidx1,:))
!!$ TR       soil%glacier_precip_minus_evap_acc(kidx0:kidx1,:) = soil%glacier_precip_minus_evap_acc(kidx0:kidx1,:) + &
!!$ TR                                                           glacier_precip_minus_evap * RhoH2O
!!$ TR    END WHERE

    ! INPUT for HD-Model in coupled case:----------
!!$ TR    glac_runoff_evap = 0._dp
!!$ TR    surf_runoff_hd = 0._dp
!!$ TR    drainage_hd = 0._dp

!!$ TR    CALL average_tiles(glacier_precip_minus_evap(1:nidx,1:ntiles), surface%is_present(kidx0:kidx1,1:ntiles), & 
!!$ TR                       surface%cover_fract(kidx0:kidx1,1:ntiles), glac_runoff_evap(1:nidx))
!!$ TR    CALL average_tiles(surface_runoff(1:nidx,1:ntiles), surface%is_present(kidx0:kidx1,1:ntiles), &
!!$ TR                       surface%cover_fract(kidx0:kidx1,1:ntiles), surf_runoff_hd(1:nidx))
!!$ TR    CALL average_tiles(drainage(1:nidx,1:ntiles), surface%is_present(kidx0:kidx1,1:ntiles), &
!!$ TR                       surface%cover_fract(kidx0:kidx1,1:ntiles), drainage_hd(1:nidx))

!!$ TR    IF (useDynveg) THEN
!!$ TR       air_qsat(1:nidx) = sat_specific_humidity(surface_temperature_avg(1:nidx), surface_pressure(1:nidx))
!!$ TR       IF (ANY(air_qsat(1:nidx) > 0.99_dp*HUGE(1._dp))) &
!!$ TR          CALL message('sat_specific_humidity', 'lookup table overflow')
!!$ TR       soil%relative_humidity_air(kidx0:kidx1) = 100._dp * air_moisture(1:nidx) * &
!!$ TR         ((1._dp - GasConstantDryAir / GasConstantWaterVapor) * &
!!$ TR         air_qsat(1:nidx) + (GasConstantDryAir / GasConstantWaterVapor)) / &
!!$ TR         (air_qsat(1:nidx) * ((1._dp - GasConstantDryAir / GasConstantWaterVapor) * air_moisture(1:nidx) + &
!!$ TR         (GasConstantDryAir / GasConstantWaterVapor)))
!!$ TR    END IF

!!$ TR    DO i=1,nidx
!!$ TR       j=kidx0+i-1
!!$ TR       IF (MAXVAL(soil%snow_fract(j,:)) > EPSILON(1._dp)) THEN
!!$ TR          hlp1 = MIN(1._dp,EXP(5000._dp * ((1._dp / tmelt) - (1._dp / surface_temperature_avg(i)))))
!!$ TR          soil%snow_age(j) = MAX(0.0_dp,(soil%snow_age(j) + &
!!$ TR                             ((hlp1**10._dp + hlp1 + 0.3_dp) * delta_time * 1.e-06_dp)) * &
!!$ TR                             (1._dp - (snow(i) * delta_time * 0.5_dp))) 
!!$ TR       ELSE
!!$ TR          soil%snow_age(j) = 0._dp
!!$ TR       END IF
!!$ TR    END DO

!!$ TR set time step counter
    time_steps_soil(1:nidx) = time_steps_soil(1:nidx) + 1.0_wp

  END SUBROUTINE update_soil_icon

SUBROUTINE update_surfacetemp_icon(klon, time_step_len,              &
  &            cvdifts, pcp,                                         &
  &            pescoe, pfscoe, peqcoe, pfqcoe,                       &
  &            psold, pqsold, pdqsold,                               &
  &            pnetrad, pgrdfl,                                      &
  &            pcfh, pcair, pcsat, pfracsu, pgrdcap,                 &
  &            psnew, pqsnew)


!!$ TR  USE mo_time_control,     ONLY: time_step_len
!!$ TR#ifndef STANDALONE
!!$ TR    USE mo_physc2, ONLY: cvdifts
!!$ TR#endif
!!$ TR  USE mo_radiation_parameters,        ONLY: cemiss
!!$ TR  USE mo_constants,        ONLY: stbo, cpd, vtmpc2, alv, als

  INTEGER,  INTENT(in)    :: klon
  REAL(wp), INTENT(in)    :: time_step_len
  REAL(wp), INTENT(in)    :: cvdifts
  REAL(wp),     INTENT(in)    :: pcp(klon), pfscoe(klon), pescoe(klon), pfqcoe(klon), peqcoe(klon)
  REAL(wp),     INTENT(in)    :: psold(klon), pqsold(klon), pdqsold(klon)
  REAL(wp),     INTENT(in)    :: pnetrad(klon), pgrdfl(klon)
  REAL(wp),     INTENT(in)    :: pcfh(klon), pcair(klon), pcsat(klon), pfracsu(klon)
  REAL(wp),     INTENT(in)    :: pgrdcap(klon)
  REAL(wp),     INTENT(out)   :: psnew(klon), pqsnew(klon)
!
  INTEGER :: jl
  REAL(wp) :: zcolin(klon), zcohfl(klon), zcoind(klon), zicp(klon), zca(klon), zcs(klon)
  REAL(wp) :: pdt, pemi, pboltz
  REAL(wp) :: pc16, platev, platsu
  REAL(wp) :: ztpfac1, ztmst
  REAL(wp) :: cpd, cpv

!!$ TR#ifdef STANDALONE
!!$ TR    REAL(wp), PARAMETER :: cvdifts = 1.0_wp
!!$ TR#endif


!-------------------------------------------------------------------------------------
! Constants
  cpd = 1005.46_wp                         ! specific heat of dry air at constant pressure in J/K/kg
  cpv = 1869.46_wp                         ! specific heat of water vapour at constant pressure in J/K/kg
  ztpfac1 = cvdifts
  ztmst   = time_step_len
  pdt     = ztpfac1*ztmst                  ! zcons29 in 'old' vdiff
  pemi    = 0.996_wp                       ! emissivity: compare wirh cemiss in echam/mo_radiation_parameters.f90
  pboltz  = 5.67E-8_wp                     ! Stefan Boltzman constant: compare with stbo in echam/mo_constants.f90
  pc16    = cpd * (cpv / cpd - 1.0_wp)     ! cpd * vtmpc2 : compare with echam/mo_constants.f90
  platev  = 2.5008e6_wp                    ! latent heat for vaporisation in J/kg: compare with alv in echam/mo_constants.f90
  platsu  = 2.8345e6_wp                    ! latent heat for sublimation in J/kg: compare with als in echam/mo_constants.f90

!************************************************************************************
!
     zicp(:) = 1._wp / pcp(:)
!
     zca(:)    = platsu * pfracsu(:) +  platev * (pcair(:) - pfracsu(:))
     zcs(:)    = platsu * pfracsu(:) +  platev * (pcsat(:) - pfracsu(:))
!
     zcolin(:) = pgrdcap(:)*zicp(:) +                                             &
          &      pdt * (zicp(:) * 4._wp * pemi * pboltz * ((zicp(:) * psold(:))**3) -     &
          &      pcfh(:) * (zca(:) * peqcoe(:) - zcs(:) -                      &
          &      zicp(:) * pc16 * psold(:) *                        &
          &      (pcair(:) * peqcoe(:) - pcsat(:)))*        &
          &      zicp(:) * pdqsold(:))
!
     zcohfl(:) = -pdt * pcfh(:) * (pescoe(:) - 1._wp)
!
     zcoind(:) = pdt * (pnetrad(:) + pcfh(:) * pfscoe(:) +  pcfh(:)*          &
          &      ((zca(:) * peqcoe(:) - zcs(:)) * pqsold(:) + zca(:) * pfqcoe(:) -   &
          &      zicp(:) * pc16 * psold(:) *                                &
          &      ((pcair(:) * peqcoe(:) - pcsat(:)) * pqsold(:) +    &
          &      pcair(:) * pfqcoe(:))) + pgrdfl(:))
!
    psnew(:)  = (zcolin(:) * psold(:) + zcoind(:)) / (zcolin(:) + zcohfl(:))
    pqsnew(:) = pqsold(:) + zicp(:) * pdqsold(:) * (psnew(:) - psold(:))

END SUBROUTINE update_surfacetemp_icon

SUBROUTINE update_soiltemp_icon(nidx,nsoil                   &
                         , delta_time                        &
                         , time_steps_soil                   &
                         , pts                               &
!!$ TR                         , psn                               &
                         , psodif, prgcgn                    &
                         , pgrndc, pgrndd                    &
                         , ptsoil                            &
                         , pgrndcapc, pgrndhflx              &
!!$ TR                         , lmask, ldglac)
                         )
!
!   AUTHOR:  FREDERIC HOURDIN     30/01/92
!
!            ADAPTED TO THE LMD-GCM BY JAN POLCHER  26/02/92
!            ADAPTED TO THE ECHAM-GCM BY JAN-PETER SCHULZ, MPI  03/02/96
!
!            J.-P. SCHULZ   MPI - OCTOBER 1997 :
!               ROUTINE USED FOR IMPLEMENTATION OF AN IMPLICIT
!               COUPLING BETWEEN LAND SURFACE AND ATMOSPHERE IN THE
!               ECHAM4 GCM.
!            U.SCHLESE DKRZ - NOVEMBER 1999  MODIFIED FOR ECHAM5
!            U.Schlese DKRZ - February 2000  new soil temperatures
!            L Kornblueh, MPI, January 2003, removed MERGE
!            
!            Adapted to JSBACH by Thomas Raddatz, Mai 2004
!            Adapted to ICON  by Thomas Raddatz, Sep 2011
!   OBJECTIVE:  COMPUTATION OF:
!               THE GROUND TEMPERATURE EVOLUTION
!               THE GROUND SPECIFIC HEAT "CAPCAL"
!               THE SURFACE DIFFUSIVE FLUX FROM GROUND "F0"
!
!
!   METHOD:  IMPLICIT TIME INTEGRATION
!
!   CONSECUTIVES GROUND TEMPERATURES ARE RELATED BY:
!           T(K+1) = C(K) + D(K)*T(K)  (1)
!   THE COEFFICIENTS C (=pgrndc) AND D (=pgrndd) ARE COMPUTED AT THE
!   T-DT TIME-STEP.
!   ROUTINE STRUCTURE:
!   1)NEW TEMPERATURES ARE COMPUTED  USING (1)
!   2)C AND D COEFFICIENTS ARE COMPUTED FROM THE NEW TEMPERATURE
!     PROFILE FOR THE T+DT TIME-STEP
!   3)THE COEFFICIENTS A AND B ARE COMPUTED WHERE THE DIFFUSIVE
!     FLUXES AT THE T+DT TIME-STEP IS GIVEN BY
!            FDIFF = A + B TS(T+DT)
!     OR     FDIFF = F0 + CAPCAL (TS(T+DT)-TS(T))/DT
!            WITH F0 = A + B (TS(T))
!                 CAPCAL = B*DT
!
!     ------------------------------------------------------------------
!
!   DECLARATIONS:
!
!!$ TR USE mo_jsbach_constants   , ONLY: RhoH2O
!!$ TR USE mo_time_control       , ONLY: lstart
!
!-----------------------------------------------------------------------
!  ARGUMENTS
!
  INTEGER, Intent(in)  ::  nidx                 !! length of the vector
  INTEGER, Intent(in)  ::  nsoil                !! number of soil layers (fixed to 5)
  REAL(wp), Intent(in)     ::  delta_time       !! time step
  REAL(wp), Intent(in)     ::  time_steps_soil(nidx)   !! number of time steps since initialisation of the soil
  REAL(wp), Intent(in)     ::  pts(nidx)            !! surface temperature at top of soil [K]
!!$ TR  REAL(wp), Intent(in)     ::  psn(nidx)            !! equivalent snow depth [m water]
  REAL(wp), Intent(in)     ::  psodif(nidx)         !! soil temperature diffusivity [m^2/s]
  REAL(wp), Intent(in)     ::  prgcgn(nidx)         !! soil heat capacity [J/m^3K]
  REAL(wp), Intent(inout)  ::  pgrndc(nidx,nsoil)   !!
  REAL(wp), Intent(inout)  ::  pgrndd(nidx,nsoil)   !!
  REAL(wp), Intent(inout)  ::  ptsoil(nidx,nsoil)   !! soil temperature [K]
  REAL(wp), Intent(out)    ::  pgrndcapc(nidx)      !!
  REAL(wp), Intent(out)    ::  pgrndhflx(nidx)      !! ground heat flux
!!$ TR  LOGICAL, Intent(in)  ::  lmask(nidx)
!!$ TR  LOGICAL, Intent(in)  ::  ldglac(nidx)         !! glacier mask
!
!     ------------------------------------------------------------------
!
!  local Variables
!
  INTEGER :: jk
  REAL(wp) :: zso_cond(nidx), zso_capa(nidx)
  REAL(wp) :: z1(nidx)
  REAL(wp) :: zd1(nsoil)
  REAL(wp) :: zdz1(nidx,nsoil),   zdz2(nidx,nsoil)
  REAL(wp) :: zkappa(nidx,nsoil), zcapa(nidx,nsoil)
  REAL(wp) :: zsnow_h(nidx), zx1(nidx), zx2(nidx)
  REAL(wp) :: cdel(nsoil),cmid(nsoil)
  REAL(wp) :: zrici, zdifiz, zsn_cond, zsn_dens, zsn_capa
!
!     ------------------------------------------------------------------
!
!*    1.  SPECIFYING THE DEPTHS OF THE TEMPERATURE LEVELS.
!
!*    1.1 SOME CONSTANTS.
!
  zrici = 2.09e+06_wp                                 !! volumetric heat capacity of ice [j/m**3/k]
  zdifiz = 12.e-07_wp                                 !! temperature diffusivity of ice  [m**2/s]
  zsn_cond = 0.31_wp                                  !! snow thermal conductivity [j/s/m/k]
  zsn_dens = 330.0_wp                                 !! snow density              [kg/m**3]
  zsn_capa = 634500.0_wp                              !! snow  heat capacity   [j/m**3/k]
  cdel = (/0.065_wp,0.254_wp,0.913_wp,2.902_wp,5.700_wp/)         !! thicknesses of soil layers [m]
  cmid = (/0.0325_wp,0.192_wp,0.7755_wp,2.683_wp,6.984_wp/)       !! depth of mids of soil layers [m]
!
!*    1.2 COMPUTING SOME USEFUL CONSTANTS.
!
  DO jk = 1,nsoil-1
     zd1(jk) = 1._wp / (cmid(jk+1) - cmid(jk))
  END DO
!
!*    1.3 COMPUTE OF THE SOIL THERMAL CONDUCTIVITY [J/S/M/K] FROM
!*        THE SOIL TEMPERATURE DIFFUSIVITY [M**2/S].
!
!!$ TR  WHERE (lmask(:) .AND. ldglac(:))
!!$ TR    zso_capa(:) = zrici
!!$ TR    zso_cond(:) = zso_capa(:) * zdifiz
!!$ TR  ELSEWHERE (lmask(:))
  zso_capa(:) = prgcgn(:)
  zso_cond(:) = zso_capa(:) * psodif(:)
!
!*    1.4 PRE-SET THERMAL CONDUCTIVITY AT ALL LEVELS.
!
!!$ TR  DO jk = 1,nsoil
!!$ TR     WHERE (lmask)
!!$ TR        zkappa(:,jk) = zso_cond(:)
!!$ TR        zcapa(:,jk)  = zso_capa(:)
!!$ TR     END WHERE
!!$ TR  END DO

  DO jk = 1,nsoil
     zkappa(:,jk) = zso_cond(:)
     zcapa(:,jk)  = zso_capa(:)
  END DO

!
!   --------------------------------------------------------------
!   COMPUTATION OF THE GROUND TEMPERATURES USING THE CGRD AND DGRD
!   COEFFICIENTS COMPUTED AT THE PREVIOUS TIME-STEP
!   --------------------------------------------------------------
!
!!$ TR  IF(.NOT.lstart) THEN
!   Upper layer
!
!!$ TR    ptsoil(:,1) = pts(:)
    WHERE (time_steps_soil(:) > 0.5_wp) ptsoil(:,1) = pts(:)
!
!   Deeper layers
!
    DO jk = 1,nsoil-1
!!$ TR      WHERE(lmask) ptsoil(:,jk+1) = pgrndc(:,jk) + pgrndd(:,jk) * ptsoil(:,jk)
      WHERE (time_steps_soil(:) > 0.5_wp) ptsoil(:,jk+1) = pgrndc(:,jk) + &
        pgrndd(:,jk) * ptsoil(:,jk)
    END DO
!!$ TR  END IF
!
!   ---------------------------------------------------------------
!   COMPUTATION OF THE CGRD AND DGRD COEFFICIENTS FOR THE NEXT STEP
!   ---------------------------------------------------------------
!
!!$ TR  WHERE (lmask) 
!!$ TR    zsnow_h(:) = psn(:) * RhoH2O / zsn_dens
     zsnow_h(:) = 0._wp
     
!
!*       Special treatment for first layer
!
     WHERE ( zsnow_h(:) > cmid(2) )
        zcapa(:,1) = zsn_capa
        zkappa(:,1) = zsn_cond
     ELSEWHERE( zsnow_h(:) > 0.0_wp .AND. zsnow_h(:) <= cmid(2) )
        zx1 = zsnow_h(:) / cmid(2)
        zx2 = ( cmid(2) - zsnow_h(:)) / cmid(2)
        zcapa(:,1) = zx1 * zsn_capa + zx2 * zso_capa(:)
        zkappa(:,1) = 1._wp / ( zx1 / zsn_cond + zx2 / zso_cond(:) )
     ELSEWHERE
        zcapa(:,1) = zso_capa(:)
        zkappa(:,1) = zso_cond(:)
     ENDWHERE
!!$ TR  END WHERE
!
  DO jk = 2, nsoil - 2
!!$ TR    WHERE (lmask)
       WHERE ( zsnow_h(:) > cmid(jk+1) )
          zcapa(:,jk) = zsn_capa
          zkappa(:,jk) = zsn_cond
       ELSEWHERE ( zsnow_h(:) > cmid(jk) .AND. zsnow_h(:) <= cmid(jk+1) )
          zx1 = (zsnow_h(:) - cmid(jk)) * zd1(jk)
          zx2 = ( cmid(jk+1) - zsnow_h(:)) * zd1(jk)
          zcapa(:,jk) = zx1 * zsn_capa + zx2 * zso_capa(:)
          zkappa(:,jk) = 1._wp / ( zx1 / zsn_cond + zx2 / zso_cond(:) )
       ELSEWHERE
          zcapa(:,jk) = zso_capa(:)
          zkappa(:,jk) = zso_cond(:)
       ENDWHERE
!!$ TR    END WHERE
  END DO
!
  DO jk=1,nsoil
!!$ TR    WHERE (lmask) zdz2(:,jk) = zcapa(:,jk) * cdel(jk) / delta_time
     zdz2(:,jk) = zcapa(:,jk) * cdel(jk) / delta_time
  END DO
!
  DO jk=1,nsoil-1
!!$ TR    WHERE (lmask) zdz1(:,jk) = zd1(jk) * zkappa(:,jk)
     zdz1(:,jk) = zd1(jk) * zkappa(:,jk)
  END DO
!
!!$ TR  WHERE (lmask)
     z1(:) = zdz2(:,nsoil) + zdz1(:,nsoil-1)
     pgrndc(:,nsoil-1) = zdz2(:,nsoil) * ptsoil(:,nsoil) / z1(:)
     pgrndd(:,nsoil-1) = zdz1(:,nsoil-1) / z1(:)
!!$ TR  END WHERE
!
  DO jk=nsoil-1,2,-1
!!$ TR     WHERE (lmask)
        z1(:) = 1._wp / (zdz2(:,jk) + zdz1(:,jk-1) + zdz1(:,jk) * (1._wp - pgrndd(:,jk)))
        pgrndc(:,jk-1) = (ptsoil(:,jk) * zdz2(:,jk) + zdz1(:,jk) * pgrndc(:,jk)) * z1(:)
        pgrndd(:,jk-1) = zdz1(:,jk-1) * z1(:)
!!$ TR     END WHERE
  END DO
!
!   ---------------------------------------------------------
!   COMPUTATION OF THE SURFACE DIFFUSIVE FLUX FROM GROUND AND
!   CALORIFIC CAPACITY OF THE GROUND:
!   ---------------------------------------------------------
!
!!$ TR  WHERE (lmask)
     pgrndhflx(:) = zdz1(:,1) * (pgrndc(:,1) + (pgrndd(:,1) - 1._wp) * ptsoil(:,1))
     pgrndcapc(:) = (zdz2(:,1) * delta_time + delta_time * (1._wp - pgrndd(:,1)) * zdz1(:,1))
!!$ TR  ELSEWHERE
!!$ TR     pgrndhflx = 0._wp
!!$ TR     pgrndcapc = 0._wp
!!$ TR  END WHERE

END SUBROUTINE update_soiltemp_icon

  FUNCTION calc_moist_root_zone(soil_moisture, soil_depth, root_depth) RESULT(root_moisture)

    ! Calculate moisture in root zone
    ! This assumes for now that depths of root zones are equal to depths of soil moisture layers, 
    ! except the last root zone

    REAL(wp), INTENT(in) :: soil_moisture(:,:)    ! nland x nsoil
    REAL(wp), INTENT(in) :: soil_depth(:,:)       ! nland x nsoil
    REAL(wp), INTENT(in) :: root_depth(:,:)       ! nland x nroot_zones
    REAL(wp)             :: root_moisture(SIZE(soil_moisture,1))

    INTEGER :: n_soil, n_root

    n_soil = SIZE(soil_moisture,2)
    n_root = SIZE(root_depth,2)

    root_moisture(:) = SUM(soil_moisture(:,1:n_root-1), DIM=2)
    root_moisture(:) = root_moisture(:) + soil_moisture(:,n_root) * root_depth(:,n_root) / &
                       soil_depth(:,n_root)

  END FUNCTION calc_moist_root_zone

  ELEMENTAL FUNCTION calc_water_stress_factor(moisture, moisture_max, fract_critical, &
                                              fract_wilting) RESULT(stress)

    REAL(wp), INTENT(in)  :: moisture
    REAL(wp), INTENT(in)  :: moisture_max
    REAL(wp), INTENT(in)  :: fract_critical, fract_wilting
    REAL(wp)              :: stress

    REAL(wp) :: moisture_critical, moisture_wilting

    moisture_critical = fract_critical * moisture_max
    moisture_wilting  = fract_wilting  * moisture_max

    stress = MAX(0._wp, MIN(1._wp, (moisture - moisture_wilting) / &
             (moisture_critical - moisture_wilting)))

  END FUNCTION calc_water_stress_factor

!!$  FUNCTION calc_relative_humidity(moisture, moisture_max, fract_wilting) RESULT(rel_humidity)
!!$
!!$    USE mo_constants,          ONLY: api
!!$
!!$    REAL(dp), INTENT(in)  :: moisture(:,:)
!!$    REAL(dp), INTENT(in)  :: moisture_max(:,:)
!!$    REAL(dp), INTENT(in)  :: fract_wilting
!!$    REAL(dp)              :: rel_humidity(SIZE(moisture,1),SIZE(moisture,2))
!!$
!!$    REAL(dp) moisture_limit(SIZE(moisture,1),SIZE(moisture,2)), evap_stop(SIZE(moisture,1),SIZE(moisture,2))
!!$
!!$    evap_stop      = MIN(0.1_dp, moisture_max)
!!$    moisture_limit = moisture_max - evap_stop
!!$
!!$    WHERE (moisture > moisture_limit .AND. moisture > fract_wilting*moisture_max)
!!$       rel_humidity = 0.5_dp * (1._dp - COS((moisture-moisture_limit) * api / evap_stop))
!!$    ELSEWHERE
!!$       rel_humidity = 0._dp
!!$    END WHERE
!!$
!!$  END FUNCTION calc_relative_humidity

  FUNCTION calc_relative_humidity_upper(moisture, moisture_max, ntiles) RESULT(rel_humidity)

    REAL(wp), INTENT(in)  :: moisture(:)
    REAL(wp), INTENT(in)  :: moisture_max(:)
    INTEGER, INTENT(in)   :: ntiles              !! Number of soil tiles
    REAL(wp)              :: rel_humidity(SIZE(moisture), ntiles)

    INTEGER               :: i
    REAL(wp)              :: moisture_upper(SIZE(moisture))
    REAL(wp)              :: rhum(SIZE(moisture))
    REAL(wp),PARAMETER    :: api = 3.141592653589793_wp ! pi 

    moisture_upper = DMIN1(moisture, moisture_max)
    WHERE (moisture_upper > 0._wp)
       rhum = 0.5_wp * (1._wp - COS((moisture_upper) * api / moisture_max))
    ELSEWHERE
       rhum = 0._wp
    END WHERE
    DO i = 1, ntiles
      rel_humidity(:,i) = rhum(:)
    ENDDO
  
  END FUNCTION calc_relative_humidity_upper

  FUNCTION sat_specific_humidity(temp, pressure) RESULT(qsat)
    !
    ! Returns saturation specific humidity for given temperature and pressure (at some atmospheric level or at surface)
    ! Uses Eq. 2.27 of Fundamentals of Atmospheric Modeling for the saturated case, but the saturation vapor pressure
    ! of water over a liquid surface resp. ice (see pp 32-34 of Ref.) is computed as in ECHAM5 (mo_convect_tables)
    !
    USE mo_convect_tables, ONLY: tlucua, &   ! Table for water vapor pressure multiplied by R_d/R_v
                                 jptlucu1, jptlucu2, lookuperror

    REAL(wp), INTENT(in) :: temp(:)      ! Air temperature at level [K]
    REAL(wp), INTENT(in) :: pressure(:)  ! Pressure at level [Pa]
    REAL(wp)             :: qsat(SIZE(temp,1))

    REAL(wp) :: water_vapor_pressure  ! Water vapor pressure [Pa]
    REAL(wp) :: vtmpc1
    INTEGER :: it(SIZE(temp,1))
    REAL(wp) :: tluc(SIZE(temp,1))

    vtmpc1 = 0.6077686814143877_wp
    it = NINT(temp*1000._wp)
    tluc = tlucua(it)

    WHERE (it >= jptlucu1 .AND. it <= jptlucu2)
       qsat = tluc / (pressure - vtmpc1*tluc)
    ELSEWHERE
       qsat = HUGE(1._wp)
    END WHERE

  END FUNCTION sat_specific_humidity

END MODULE mo_soil_icon

