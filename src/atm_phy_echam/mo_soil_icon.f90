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

  IMPLICIT NONE

  PRIVATE
!!$ TR  PUBLIC :: soil_type, soil_param_type, init_soil, update_soil, soil_diagnostics, get_soil_diag
  PUBLIC :: update_soil_icon

CONTAINS

  SUBROUTINE update_soil_icon(ntiles, kidx, &
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
       moisture1, moisture2, moisture3, moisture4, moisture5, &
       sat_surface_specific_humidity, skin_reservoir, &
       snow_fract, snow, snow_canopy, snow_melt, snow_acc, snow_melt_acc, &
       glacier_runoff_acc, runoff_acc, drainage_acc, &
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
    INTEGER, INTENT(in) :: kidx
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
!!$ TR         canopy_snow_fract
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
    REAL(wp), DIMENSION(kidx) :: &
         net_radiation, &                        ! Net radiation at surface
         dry_static_energy, &                    ! Surface dry static energy (= C_p * T )
         dry_static_energy_new, &                ! New dry static energy
         surface_qsat_new, &                     ! New surface saturated specific humidity
         canopy_snow_fract, &                    !
!!$ TR         air_qsat, &                             ! Saturated specific humidity at lowest atmospheric level
!!$ TR         evapotranspiration_no_snow_skin, &      ! Evapotranspiration without that from snow and the skin reservoir
         zdqsl, &                                ! Sensitivity of saturated surface specific humidity wrt temperature
         zcpq                                   ! Conversion factor for humidity from dry static energy
!!$ TR         snow_avg, &                             ! Snow [m water equivalent] averaged over all tiles
!!$ TR         melt_water_excess, &                    ! water from snow melting which exceeds skin reservoir and infiltrates soil
!!$ TR         csat_transpiration                      ! fraction of grid box that contributes to transpiration
                                                 ! (considered to be completely wet)
!!$ TR    REAL(wp), DIMENSION(kidx) :: &
!!$ TR         sat_surf_specific_hum_avg, &
!!$ TR         ground_heat_flux_avg, &
!!$ TR         heat_capacity_avg, &
!!$ TR         surface_temperature_avg, soil_temperature_avg, &
!!$ TR         soil_moisture_avg

    REAL(wp), DIMENSION(kidx,ntiles) :: &
!!$ TR         skin_reservoir_max, &          
!!$ TR         wet_skin_fract, &
!!$ TR         glacier_precip_minus_evap, &       ! P-E for glaciers [m]
!!$ TR         surface_runoff, drainage, &        ! Surface runoff and drainage for HD model [m]
         qsat_fact, qair_fact, &            ! Factors for implicit coupling
         air_dry_static_energy_new, &
         air_moisture_new, &
!!$ TR         evaporation_pot, &                 ! Potential evaporation
         canopy_resistance, &               ! Water-limited canopy resistance
         canopy_conductance_limited, &      ! water limited canopy conductance
         soil_moisture_root, &              ! Soil moisture in root zone
         soil_moisture_root_max, &          ! Field capacity in root zone
         water_stress_factor, &             ! Water stress factor (=1: no stress, =0: infinite stress)
!!$ TR         relative_humidity, &               ! Relative humidity (Eq. 3.3.2.9 ECHAM3 Manual)
!!$ TR         zhsoil, tte_corr, &
!!$ TR         sat_surface_specific_hum_old, &
         qsat_veg, qair_veg
!!$ TR         csat_tiles, cair_tiles, &
!!$ TR         qsat_transpiration, csat_transpiration_tiles
!!$ TR#ifdef __PGI
!!$ TR    REAL(wp), DIMENSION(kidx,soil%ntiles) :: moisture_max
!!$ TR#endif
    REAL(wp), DIMENSION(kidx,5,ntiles) :: soil_moisture ! substitute for soil%moisture
    REAL(wp), DIMENSION(kidx,5,ntiles) :: root_depth
    REAL(wp), DIMENSION(kidx,5) :: soil_depth
    REAL(wp), DIMENSION(kidx,5) :: soil_MaxMoisture
!!$ TR    LOGICAL, DIMENSION(kidx,soil%ntiles) :: &
!!$ TR         soil_mask                          ! True if not glacier and not lake

    REAL(wp) :: zcons30, ztpfac2, ztpfac3, vtmpc2
!!$ TR    REAL(wp) :: hlp1
    INTEGER :: nidx, itile, i, j

!!$ TR#ifdef STANDALONE
!!$ TR    REAL(wp), PARAMETER :: cvdifts = 1.0_dp
!!$ TR#endif

    !! local variables for testing
    REAL(wp), DIMENSION(kidx) :: sat_surface_specific_hum_old
    REAL(wp), DIMENSION(kidx) :: surface_temperature_unfiltered
    INTEGER  ::  ntsoil
    REAL(wp), DIMENSION(kidx) :: ThermalDiffusivity
    REAL(wp), DIMENSION(kidx) :: VolHeatCapacity
    REAL(wp)  ::  cvdifts
    REAL(wp)  ::  eps
    REAL(wp)  ::  Emissivity
    REAL(wp)  ::  StefanBoltzmann
    REAL(wp)  ::  Gravity
    REAL(wp)  ::  SpecificHeatDryAirConstPressure
    REAL(wp)  ::  SpecificHeatVaporConstPressure
    REAL(wp)  ::  moist_crit_fract, moist_wilt_fract


    REAL(wp), DIMENSION(kidx,5)  ::  c_soil_temperature
    REAL(wp), DIMENSION(kidx,5)  ::  d_soil_temperature
    REAL(wp), DIMENSION(kidx,5)  ::  soil_temperature

    EXTERNAL update_surfacetemp_icon, update_soiltemp_icon

    nidx = kidx
    ntsoil = 5
    ThermalDiffusivity(1:nidx) = 7.e-7_wp
    VolHeatCapacity(1:nidx) = 2.e+6_wp
    cvdifts = 1.5_wp
    eps = 0.1_wp
    Emissivity = 0.996_wp
    StefanBoltzmann = 5.67e-8_wp
    Gravity = 9.80665_wp
    SpecificHeatDryAirConstPressure = 1005.46_wp
    SpecificHeatVaporConstPressure = 1869.46_wp

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

    root_depth(:,:,:) = 0.1_wp
    soil_depth(:,:) = 0.1_wp
    soil_MaxMoisture(:,:) = 0.1_wp

    DO itile=1,ntiles
       soil_moisture(1:nidx,1,itile) = moisture1(1:nidx)
       soil_moisture(1:nidx,2,itile) = moisture2(1:nidx)
       soil_moisture(1:nidx,3,itile) = moisture3(1:nidx)
       soil_moisture(1:nidx,4,itile) = moisture4(1:nidx)
       soil_moisture(1:nidx,5,itile) = moisture5(1:nidx)
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

!!$ TR    evaporation_pot = 0._dp
!!$ TR    do itile=1,ntiles
!!$ TR    do i=1,nidx
!!$ TR    j=kidx0+i-1
!!$ .   IF (skin_reservoir_max(i,itile) > EPSILON(1._dp)) THEN
!!$ .       canopy_snow_fract(i,itile) = MIN(1._dp, canopy_snow(i,itile) / skin_reservoir_max(i,itile))
!!$ .    ELSE
!!$       canopy_snow_fract(i,itile) = 0._dp
!!$    END IF
!!$    IF (soil_mask(i,itile)) THEN
!!$       wet_skin_fract(i,itile) = MERGE(MIN(1._dp, soil%skin_reservoir(j,itile) / skin_reservoir_max(i,itile)), &
!!$            0.0_dp, soil_options%SkinReservoirMax > EPSILON(1._dp))
!!$      
!!$       ! Roesch et al, 2002, Climate Dynamics
!!$       soil%snow_fract(j,itile) = zqsncr * TANH(soil%snow(j,itile) * 100._dp) * &
!!$            SQRT(soil%snow(j,itile) * 1000._dp / (soil%snow(j,itile) * 1000._dp + zepsec + &
!!$            zsigfac * surface%oro_std_dev(j)))
!!$       IF (soil%snow_fract(j,itile) .LT. EPSILON(1._dp) .AND. canopy_snow_fract(i,itile) .GE. EPSILON(1._dp)) &
!!$            soil%snow_fract(j,itile) = canopy_snow_fract(i,itile)
!!$       ! Modify snow cover if snow loss during the time step due to potential evaporation is larger than 
!!$       ! equivalent snow water content from soil and canopy; same for skin reservoir
!!$       ! Potential evaporation using old values of air and surface humidity
!!$       evaporation_pot(i,itile) = zcons30 * cdrag(i) * &
!!$            (soil%sat_surface_specific_humidity(j,itile) - air_moisture(i))
!!$       IF (soil%snow_fract(j,itile) > 0._dp) THEN
!!$          soil%snow_fract(j,itile) = soil%snow_fract(j,itile) / &
!!$               MAX(1._dp, soil%snow_fract(j,itile) * evaporation_pot(i,itile) * delta_time / &
!!$               (RhoH2O*(soil%snow(j,itile) + canopy_snow(i,itile))))
!!$       END IF
!!$       IF (wet_skin_fract(i,itile) > 0._dp) THEN
!!$          wet_skin_fract(i,itile) = wet_skin_fract(i,itile) / &
!!$               MAX(1._dp, (1._dp - soil%snow_fract(j,itile)) * evaporation_pot(i,itile) * delta_time / &
!!$               (RhoH2O * MAX(EPSILON(1._dp),soil%skin_reservoir(j,itile))))
!!$       END IF
!!$       
!!$    ELSE IF (surface%is_glacier(j,itile)) THEN
!!$       wet_skin_fract(i,itile) = 0.0_dp
!!$       soil%snow_fract(j,itile) = 1.0_dp
!!$ .    ELSE
!!$ .       wet_skin_fract(i,itile) = 0.0_dp
!!$ .       soil%snow_fract(j,itile) = 0.0_dp
!!$ TR    END IF
!!$ TR    END DO
!!$ TR    END do

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

    DO itile=1,ntiles
       ! Compute temperature and moisture at lowest atmospheric level by back-substitution
       air_dry_static_energy_new(:,itile) = t_Acoef(:) * dry_static_energy_new(:) + t_Bcoef(:)
!!$ TR       air_moisture_new(:,itile) = q_Acoef(:) * soil%sat_surface_specific_humidity(kidx0:kidx1,itile) + q_Bcoef(:)
       air_moisture_new(:,itile) = q_Acoef(:) * surface_qsat_new(:) + q_Bcoef(:)
    END DO

!!$ TR    DO itile=1,ntiles
!!$ TR       soil%radiative_temperature(kidx0:kidx1,itile) = dry_static_energy_new(:) / zcpq(:)
!!$ TR    END DO
    surface_temperature_rad(1:nidx) = dry_static_energy_new(1:nidx) / zcpq(1:nidx)

    ! Compute sensible heat flux
!!$ TR    DO itile=1,ntiles
!!$ TR       soil%sensible_heat_flux(kidx0:kidx1,itile) = zcons30 * cdrag(:) *              &
!!$ TR         (air_dry_static_energy_new(:,itile) - dry_static_energy_new(:) -             &
!!$ TR         SpecificHeatDryAirConstPressure * vtmpc2 * surface_temperature_avg(:) *      &
!!$ TR         (cair(:) * air_moisture_new(:,itile) - csat(:) *                             &
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
!!$ TR    END DO
!!$ TR    ! Latent heat flux
!!$ TR    soil%latent_heat_flux(kidx0:kidx1,:) = LatentHeatVaporization  * soil%evapotranspiration(kidx0:kidx1,:) &
!!$ TR    + (LatentHeatSublimation - LatentHeatVaporization) * soil%snow_fract(kidx0:kidx1,:) * soil%evaporation_pot(kidx0:kidx1,:)

    evapotranspiration(1:nidx) = zcons30 * cdrag(1:nidx) * &
      (cair(1:nidx) * (q_Acoef(1:nidx) * surface_qsat_new(1:nidx) + q_Bcoef(1:nidx)) - &
      csat(1:nidx)  * surface_qsat_new(1:nidx))
    !
    !--------------------------------------------------------------------------------------------------------
    ! 
!!$ TR    tte_corr = 0._dp
!!$ TR    glacier_precip_minus_evap = 0._dp
!!$ TR    melt_water_excess = 0._dp
!!$ TR    surface_runoff = 0._dp
!!$ TR    drainage = 0._dp

!!$ TR    DO itile=1,ntiles
!!$ TR       CALL update_surf_down( nidx, delta_time,                                                             &
!!$ TR            air_moisture(1:nidx),                                                                           &
!!$ .            soil%surface_temperature_unfiltered(kidx0:kidx1,itile), wind10(1:nidx), air_temperature(1:nidx), &
!!$ .            soil%skin_reservoir(kidx0:kidx1,itile), wet_skin_fract(1:nidx,itile), skin_reservoir_max(1:nidx,itile), &
!!$ .            soil%snow(kidx0:kidx1,itile), soil%snow_fract(kidx0:kidx1,itile),                                &
!!$            canopy_snow(1:nidx,itile), soil%glacier_depth(kidx0:kidx1,itile),                                &
!!$            soil%heat_capacity(kidx0:kidx1,itile),                                                           &
!!$            soil%evapotranspiration(kidx0:kidx1,itile), soil%evaporation_pot(kidx0:kidx1,itile),             &
!!$            evapotranspiration_no_snow_skin(1:nidx),                                                         &
!!$            precip_rain, precip_snow,                                                                        &
!!$            surface%is_present(kidx0:kidx1,itile),                                                           &
!!$            surface%is_glacier(kidx0:kidx1,itile),                                                           &
!!$            glacier_precip_minus_evap(:,itile),                                                              &
!!$            soil%snow_acc(kidx0:kidx1,itile), soil%snow_melt_acc(kidx0:kidx1,itile),                         &
!!$            soil%glacier_runoff_acc(kidx0:kidx1,itile), melt_water_excess(1:nidx),                           &
!!$            tte_corr(1:nidx,itile),soil%snow_melt(kidx0:kidx1,itile))
!!$
!!$       CALL update_soiltemp(nidx,soil%ntsoil                                                               &
!!$                         , soil%surface_temperature_unfiltered(kidx0:kidx1,itile),soil%snow(kidx0:kidx1,itile) &
!!$                         , soil_param%ThermalDiffusivity(kidx0:kidx1)                                      &
!!$                         , soil_param%VolHeatCapacity(kidx0:kidx1)                                         &
!!$                         , soil%c_soil_temperature(kidx0:kidx1,:,itile)                                    &
!!$                         , soil%d_soil_temperature(kidx0:kidx1,:,itile)                                    &
!!$                         , soil%soil_temperature(kidx0:kidx1,:,itile)                                      &
!!$                         , soil%heat_capacity(kidx0:kidx1,itile),soil%ground_heat_flux(kidx0:kidx1,itile)  &
!!$                         , surface%is_present(kidx0:kidx1,itile)                                           &
!!$                         , surface%is_glacier(kidx0:kidx1,itile))
!!$
!!$       CALL update_surf_up(nidx, itile, delta_time,                                                       &
!!$            soil%skin_reservoir(kidx0:kidx1,itile), wet_skin_fract(:,itile), skin_reservoir_max(:,itile), &
!!$            soil%moisture(kidx0:kidx1,1,itile), soil_param%MaxMoisture(kidx0:kidx1,1),                    &
!!$            soil%soil_temperature(kidx0:kidx1,1,itile), soil%snow_fract(kidx0:kidx1,itile),               &
!!$            surface%oro_std_dev(kidx0:kidx1),                                                             &
!!$            soil%evapotranspiration(kidx0:kidx1,itile), soil%evaporation_pot(kidx0:kidx1,itile),          &
!!$            evapotranspiration_no_snow_skin(1:nidx),                                                      &
!!$            precip_rain(1:nidx),                                                                          &
!!$            surface%is_present(kidx0:kidx1,itile),                                                        &
!!$            surface%is_glacier(kidx0:kidx1,itile),                                                        &
!!$            soil%runoff_acc(kidx0:kidx1,itile), soil%drainage_acc(kidx0:kidx1,itile),                     &
!!$            surface_runoff(1:nidx,itile), drainage(1:nidx,itile), melt_water_excess(1:nidx))
!!$       ! Note: surface_runoff and drainage have to be passed to the runoff model once this is part of jsbach
!!$    END DO
!!$    CALL average_tiles(tte_corr, surface%is_present(kidx0:kidx1,:) .AND. .NOT. surface%is_lake(kidx0:kidx1,:) , &
!!$         surface%cover_fract(kidx0:kidx1,:), tte_corr_avg)
!!$    CALL average_tiles(soil%moisture(kidx0:kidx1,1,:), surface%is_present(kidx0:kidx1,:) &
!!$                      .AND. .NOT. surface%is_lake(kidx0:kidx1,:), surface%cover_fract(kidx0:kidx1,:), &
!!$                      soil_moisture_avg)
!!$ .    CALL average_tiles(soil%snow(kidx0:kidx1,:), surface%is_present(kidx0:kidx1,:) &
!!$ .                      .AND. .NOT. surface%is_lake(kidx0:kidx1,:), surface%cover_fract(kidx0:kidx1,:), &
!!$ .                      snow_avg)
!!$ TR    CALL average_tiles(soil%surface_temperature_unfiltered(kidx0:kidx1,:), surface%is_present(kidx0:kidx1,:) &
!!$ TR                      .AND. .NOT. surface%is_lake(kidx0:kidx1,:), surface%cover_fract(kidx0:kidx1,:), &
!!$ TR                      surface_temperature_avg)
    !! for testing
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

!!$ TR CALL update_surf_down from current cosmos-landveg

    CALL update_soiltemp_icon(nidx,ntsoil                           &
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

!!$ TR CALL update_surf_up from current cosmos-landveg

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
!!$ TR    END DO
!!$ TR    DO isoil=2,soil%ntsoil
!!$ TR       CALL average_tiles(soil%soil_temperature(kidx0:kidx1,isoil,:), surface%is_present(kidx0:kidx1,:) &
!!$ TR                    .AND. .NOT. surface%is_lake(kidx0:kidx1,:), surface%cover_fract(kidx0:kidx1,:), &
!!$ TR                    soil_temperature_avg)
!!$ TR       DO itile=1,ntiles
!!$ TR          soil%soil_temperature(kidx0:kidx1,isoil,itile) = soil_temperature_avg
!!$ TR       END DO
!!$ TR    END DO


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

!!$    DO itile=1,ntiles
!!$    DO i=1,nidx
!!$    j=kidx0+i-1
!!$    IF (soil_mask(i,itile)) THEN
!!$       wet_skin_fract(i,itile) = MERGE(MIN(1._dp, soil%skin_reservoir(j,itile) / skin_reservoir_max(i,itile)), &
!!$            0.0_dp, soil_options%SkinReservoirMax > EPSILON(1._dp))
!!$      
!!$       ! Roesch et al, 2002, Climate Dynamics
!!$       soil%snow_fract(j,itile) = zqsncr * TANH(soil%snow(j,itile) * 100._dp) * &
!!$            SQRT(soil%snow(j,itile)*1000. / (soil%snow(j,itile)*1000._dp + zepsec + &
!!$            zsigfac * surface%oro_std_dev(j)))
!!$
!!$       IF (soil%snow_fract(j,itile) .LT. EPSILON(1._dp) .AND. canopy_snow_fract(i,itile) .GE. EPSILON(1._dp)) &
!!$            soil%snow_fract(j,itile) = canopy_snow_fract(i,itile)
!!$       ! Modify snow cover if snow loss during the time step due to potential evaporation is larger than
!!$       ! equivalent snow water content from soil and canopy; same for skin reservoir
!!$       ! Potential evaporation using old values of air and surface humidity
!!$       evaporation_pot(i,itile) = zcons30 * cdrag(i) * &
!!$            (sat_surface_specific_hum_old(i,itile) - air_moisture(i))
!!$       IF (soil%snow_fract(j,itile) > 0._dp .AND. soil%snow(j,itile) + canopy_snow(i,itile) > 0._dp) THEN
!!$          soil%snow_fract(j,itile) = soil%snow_fract(j,itile) / &
!!$               MAX(1._dp, soil%snow_fract(j,itile) * evaporation_pot(i,itile) * delta_time / &
!!$               (RhoH2O*(soil%snow(j,itile) + canopy_snow(i,itile))))
!!$       END IF
!!$
!!$       IF (wet_skin_fract(i,itile) > 0._dp) THEN
!!$          wet_skin_fract(i,itile) = wet_skin_fract(i,itile) / &
!!$               MAX(1._dp, (1._dp - soil%snow_fract(j,itile)) * evaporation_pot(i,itile) * delta_time / &
!!$               (RhoH2O * MAX(EPSILON(1._dp),soil%skin_reservoir(j,itile))))
!!$       END IF
!!$       
!!$    ELSE IF (surface%is_glacier(j,itile)) THEN
!!$       wet_skin_fract(i,itile) = 0.0_dp
!!$       soil%snow_fract(j,itile) = 1.0_dp
!!$    ELSE
!!$       wet_skin_fract(i,itile) = 0.0_dp
!!$       soil%snow_fract(j,itile) = 0.0_dp
!!$    END IF
!!$    END DO
!!$    END DO
!!$
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

