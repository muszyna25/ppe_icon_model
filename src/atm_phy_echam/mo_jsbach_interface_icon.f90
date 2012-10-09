!> 
!! JSBACH interface adapted for use with ICON
!!
!! @author Thomas Raddatz, MPI-M
!!
!! @par Revision History
!! First version by Thomas Raddatz (2011-09)
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
!!    violent purposes.
!! <li> The code may not be re-distributed without the consent of the authors.
!! <li> The copyright notice and statement of authorship must appear in all
!!    copies.
!! <li> You accept the warranty conditions (see WARRANTY).
!! <li> In case you intend to use the code commercially, we oblige you to sign
!!    an according license agreement with DWD and MPI-M.
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
MODULE mo_jsbach_interface_icon

  USE mo_exception,         ONLY: finish, message
  USE mo_kind,              ONLY: wp
  USE mo_soil_icon,         ONLY: update_soil_icon
  USE mo_canopy,            ONLY: unstressed_canopy_cond_par
  USE mo_land_surface,      ONLY: update_albedo_echam5

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: jsbach_inter_1d

CONTAINS
  !>
  !!
  !!
  SUBROUTINE jsbach_inter_1d ( &
       ! INPUT
       kdim, &                            !! Length of vectors
       kblock, &                          !! Index of block
       kland, &                           !! Number of land points in vectors ( = SUM(lmask) )
       mask_land, &                       !! land sea mask
       delta_time, &                      !! time step
       time_step_len, &                   !! 2 * delta_time for leapfrog
       time_steps_soil, &                 !! number of time steps since initialisation of the soil
       wind, &
!!$ TR wind10, &
       temp_air, &
       qair, &
       precip_rain, &
       precip_snow, &
       lwdown, &
       swdown, &                      !! introduced for testing to replace sw_vis_net + sw_nir_net
!!$ TR       sw_vis_net, &                      !! net solar radiation in the visible band
!!$ TR       sw_vis_frac_diffuse, &             !! fraction of diffuse solar radiation contained in sw_vis_net
!!$ TR       sw_nir_net, &                      !! net solar radiation in the near infrared band
!!$ TR       sw_nir_frac_diffuse, &             !! fraction of diffuse solar radiation contained in sw_nir_net
!!$ TR       sw_par_down, &                     !! downward PAR
!!$ TR       sw_par_frac_diffuse, &             !! fraction of diffuse solar radiation contained in sw_par_down
       pressure, &
       czenith, &
!!$ TR       CO2_concentration, &
       cdrag, &
       etAcoef, &
       etBcoef, &
       eqAcoef, &
       eqBcoef, &
       p_echam_zchl, &
       albedo_vis_soil, &
       albedo_nir_soil, &
       albedo_vis_canopy, &
       albedo_nir_canopy, &
       albedo_background, &
       forest_fract, &
       !! inout for testing (hydrology)
       cair, &
       csat, &
       csat_transpiration, &
       albvisdir, &
       albnirdir, &
       albvisdif, &
       albnirdif, &
       moisture1, &
       moisture2, &
       moisture3, &
       moisture4, &
       moisture5, &
       moisture_all, &
       sat_surface_specific_humidity, &
       skin_reservoir, &
       snow_fract, &
       snow, &
       snow_canopy, &
       snow_melt, &
       snow_acc, &
       snow_melt_acc, &
       glacier_runoff_acc, &
       runoff_acc, &
       drainage_acc, &
       !! inout for testing (energy balance)
       surface_temperature, &
       surface_temperature_old, &
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
       heat_capacity, &
       ground_heat_flux, &
       swnet             &
       )

       ! Output
       ! surface and soil temperature for testing
!!$ TR       evap_act, evap_pot, sensible, latent, &
!!$ TR       CO2_flux_npp, CO2_flux_soilresp, CO2_flux_herbivory, CO2_flux_dynveg, &
!!$ TR       CO2_emission_lcc, CO2_emission_harvest, &
!!$ TR       tsoil_rad, temp_soil_new, qsurf, &
!!$ TR       albedo, albedo_vis, albedo_nir, emis, &
!!$ TR       z0h, z0m, &
!!$ TR       cdrag, &                    !! shifted up as these are input variables  
!!$ TR       etAcoef, etBcoef, eqAcoef, eqBcoef, &
!!$ TR       zhsoil, echam_zchl, &
!!$ TR       mask_land, mode, surf_dry_static_energy, &
!!$ TR       kblock, soil_wetness, snow_depth, &
!!$ TR       runoff, drainage, skin_res, &
!!$ TR       tte_corr, &
!!$ TR       glac_runoff_evap, surf_runoff_hd , drainage_hd, &
!!$ TR       glacier_depth, snow_melt_acc, glacier_p_minus_e_acc, snow_acc, glacier_runoff_acc, &
!!$ TR       wsmax, &
!!$ TR       snow_melt &  ! for dry deposition modules

    ! Subroutine arguments
    INTEGER,            INTENT(in)    :: kdim                     !! Length of vectors
    INTEGER, OPTIONAL,  INTENT(in)    :: kblock                   !! Index of block  in domain that is to be
                                                                  !! processed. If missing: one call for whole domain
    INTEGER,            INTENT(in)    :: kland                    !! Number of land points in vectors
    INTEGER, OPTIONAL,  INTENT(in)    :: mask_land(kdim)          !! Land-sea mask (land includes glaciers)
    REAL(wp), INTENT(in)              :: delta_time               !! time step
    REAL(wp), INTENT(in)              :: time_step_len            !! 2*delta_time for leapfrog
    REAL(wp), OPTIONAL, INTENT(inout) :: time_steps_soil(kdim)    !! number of time steps since initialisation of the soil
    REAL(wp), OPTIONAL, INTENT(in)    :: wind(kdim)               !! Lowest level wind speed [m/s]
!!$ TR    REAL(wp), OPTIONAL, INTENT(in)    :: wind10(kdim)             !! 10m wind speed [m/s] (for update_surface_down)
    REAL(wp), OPTIONAL, INTENT(in)    :: temp_air(kdim)           !! Lowest level air temperature [Kelvin]
    REAL(wp), OPTIONAL, INTENT(in)    :: qair(kdim)               !! Lowest level specific humidity
    REAL(wp), OPTIONAL, INTENT(in)    :: precip_rain(kdim)        !! Precipitation as rain [kg/(m^2 s)]
    REAL(wp), OPTIONAL, INTENT(in)    :: precip_snow(kdim)        !! Precipitation as snow [kg/(m^2 s)]
    REAL(wp), OPTIONAL, INTENT(in)    :: lwdown(kdim)             !! Downward longwave flux
    REAL(wp), OPTIONAL, INTENT(in)    :: swdown(kdim)             !! Downward shortwave flux
!!$ TR   REAL(wp), OPTIONAL, INTENT(in)    :: sw_vis_net(kdim)          !! net surface visible radiation [W/m^2]
!!$ TR   REAL(wp), OPTIONAL, INTENT(in)    :: sw_vis_frac_diffuse(kdim) !! fraction of diffuse radiation contained in sw_vis_net
!!$ TR   REAL(wp), OPTIONAL, INTENT(in)    :: sw_nir_net(kdim)          !! net surface NIR radiation [W/m^2]
!!$ TR   REAL(wp), OPTIONAL, INTENT(in)    :: sw_nir_frac_diffuse(kdim) !! fraction of diffuse radiation contained in sw_nir_net
!!$ TR   REAL(wp), OPTIONAL, INTENT(in)    :: sw_par_down(kdim)         !! downward surface PAR [W/m^2]
!!$ TR   REAL(wp), OPTIONAL, INTENT(in)    :: sw_par_frac_diffuse(kdim) !! fraction of diffuse radiation contained in sw_par_net
    REAL(wp), OPTIONAL, INTENT(in)    :: pressure(kdim)            !! Surface pressure
   REAL(wp), OPTIONAL, INTENT(in)    :: czenith(kdim)              !! Cosine of solar zenith angle
!!$ TR   REAL(wp), OPTIONAL, INTENT(in)    :: CO2_concentration(kdim)  !! Atmospheric CO2 concentration [kg(CO2)/kg(air)]
    REAL(wp), OPTIONAL, INTENT(in)    :: cdrag(kdim)              !! Surface drag
    REAL(wp), OPTIONAL, INTENT(in)    :: etAcoef(kdim)            !! Richtmeyer Morton coeff. temperature
    REAL(wp), OPTIONAL, INTENT(in)    :: etBcoef(kdim)            !! 
    REAL(wp), OPTIONAL, INTENT(in)    :: eqAcoef(kdim)            !! Richtmeyer Morton coeff. humidity
    REAL(wp), OPTIONAL, INTENT(in)    :: eqBcoef(kdim)            !!
    REAL(wp), OPTIONAL, INTENT(in)    :: p_echam_zchl(kdim)       !!
    REAL(wp), OPTIONAL, INTENT(in)    :: albedo_vis_soil(kdim)    !!
    REAL(wp), OPTIONAL, INTENT(in)    :: albedo_nir_soil(kdim)    !!
    REAL(wp), OPTIONAL, INTENT(in)    :: albedo_vis_canopy(kdim)  !!
    REAL(wp), OPTIONAL, INTENT(in)    :: albedo_nir_canopy(kdim)  !!
    REAL(wp), OPTIONAL, INTENT(in)    :: albedo_background(kdim)  !!
    REAL(wp), OPTIONAL, INTENT(in)    :: forest_fract(kdim)       !!
    !! inout for testing (hydrology
    REAL(wp), OPTIONAL, INTENT(inout) :: cair(kdim)                !! area fraction with wet surface
    REAL(wp), OPTIONAL, INTENT(inout) :: csat(kdim)                !! area fraction with wet surface (air)
    REAL(wp), OPTIONAL, INTENT(inout) :: csat_transpiration(kdim)
    REAL(wp), OPTIONAL, INTENT(inout) :: albvisdir(kdim)
    REAL(wp), OPTIONAL, INTENT(inout) :: albnirdir(kdim)
    REAL(wp), OPTIONAL, INTENT(inout) :: albvisdif(kdim)
    REAL(wp), OPTIONAL, INTENT(inout) :: albnirdif(kdim)
    REAL(wp), OPTIONAL, INTENT(inout) :: moisture1(kdim)
    REAL(wp), OPTIONAL, INTENT(inout) :: moisture2(kdim)
    REAL(wp), OPTIONAL, INTENT(inout) :: moisture3(kdim)
    REAL(wp), OPTIONAL, INTENT(inout) :: moisture4(kdim)
    REAL(wp), OPTIONAL, INTENT(inout) :: moisture5(kdim)
    REAL(wp), OPTIONAL, INTENT(inout) :: moisture_all(kdim)
    REAL(wp), OPTIONAL, INTENT(inout) :: sat_surface_specific_humidity(kdim)
    REAL(wp), OPTIONAL, INTENT(inout) :: skin_reservoir(kdim)
    REAL(wp), OPTIONAL, INTENT(inout) :: snow_fract(kdim)
    REAL(wp), OPTIONAL, INTENT(inout) :: snow(kdim)
    REAL(wp), OPTIONAL, INTENT(inout) :: snow_canopy(kdim)
    REAL(wp), OPTIONAL, INTENT(inout) :: snow_melt(kdim)
    REAL(wp), OPTIONAL, INTENT(inout) :: snow_acc(kdim)
    REAL(wp), OPTIONAL, INTENT(inout) :: snow_melt_acc(kdim)
    REAL(wp), OPTIONAL, INTENT(inout) :: glacier_runoff_acc(kdim)
    REAL(wp), OPTIONAL, INTENT(inout) :: runoff_acc(kdim)
    REAL(wp), OPTIONAL, INTENT(inout) :: drainage_acc(kdim)
    !! inout for testing (energy balance)
    REAL(wp), OPTIONAL, INTENT(inout) :: surface_temperature(kdim) !! surface_temperature
    REAL(wp), OPTIONAL, INTENT(inout) :: surface_temperature_old(kdim) !! old surface_temperature
    REAL(wp), OPTIONAL, INTENT(out) :: surface_temperature_rad(kdim) !! old surface_temperature
    REAL(wp), OPTIONAL, INTENT(out) :: evapotranspiration(kdim) !! evapotranspiration
    REAL(wp), OPTIONAL, INTENT(inout) :: c_soil_temperature1(kdim) !! soil temperature parameter
    REAL(wp), OPTIONAL, INTENT(inout) :: c_soil_temperature2(kdim) !! soil temperature parameter
    REAL(wp), OPTIONAL, INTENT(inout) :: c_soil_temperature3(kdim) !! soil temperature parameter
    REAL(wp), OPTIONAL, INTENT(inout) :: c_soil_temperature4(kdim) !! soil temperature parameter
    REAL(wp), OPTIONAL, INTENT(inout) :: c_soil_temperature5(kdim) !! soil temperature parameter
    REAL(wp), OPTIONAL, INTENT(inout) :: d_soil_temperature1(kdim) !! soil temperature parameter
    REAL(wp), OPTIONAL, INTENT(inout) :: d_soil_temperature2(kdim) !! soil temperature parameter
    REAL(wp), OPTIONAL, INTENT(inout) :: d_soil_temperature3(kdim) !! soil temperature parameter
    REAL(wp), OPTIONAL, INTENT(inout) :: d_soil_temperature4(kdim) !! soil temperature parameter
    REAL(wp), OPTIONAL, INTENT(inout) :: d_soil_temperature5(kdim) !! soil temperature parameter
    REAL(wp), OPTIONAL, INTENT(inout) :: soil_temperature1(kdim) !! soil temperature
    REAL(wp), OPTIONAL, INTENT(inout) :: soil_temperature2(kdim) !! soil temperature
    REAL(wp), OPTIONAL, INTENT(inout) :: soil_temperature3(kdim) !! soil temperature
    REAL(wp), OPTIONAL, INTENT(inout) :: soil_temperature4(kdim) !! soil temperature
    REAL(wp), OPTIONAL, INTENT(inout) :: soil_temperature5(kdim) !! soil temperature
    REAL(wp), OPTIONAL, INTENT(inout) :: heat_capacity(kdim)     !! heat capacity
    REAL(wp), OPTIONAL, INTENT(inout) :: ground_heat_flux(kdim)  !! ground heat flux
    REAL(wp), OPTIONAL, INTENT(inout) :: swnet(kdim)  !! net surface solar flux


    !! Local declarations for packed (gathered) input fields
    !!$ TR JSBACH testing: replace kdim with kland
    REAL(wp), DIMENSION(kdim) ::   zwind
!! TR    REAL(wp), DIMENSION(kdim) ::   zwind10
    REAL(wp), DIMENSION(kdim) ::   ztemp_air
    REAL(wp), DIMENSION(kdim) ::   zqair
    REAL(wp), DIMENSION(kdim) ::   zprecip_rain         !! Packed rain
    REAL(wp), DIMENSION(kdim) ::   zprecip_snow         !! Packed snow
    REAL(wp), DIMENSION(kdim) ::   zlwdown                            !! Packed downward longwave radiation
    REAL(wp), DIMENSION(kdim) ::   zswdown                            !! Packed downward shortwave radiation
!!$ TR    REAL(wp), DIMENSION(kdim) ::   zsw_vis_net                        !! Packed radiation from visible band [W/m^2]
!!$ TR    REAL(wp), DIMENSION(kdim) ::   zsw_vis_frac_diffuse               !! Packed fraction of diffuse visible radiation
!!$ TR    REAL(wp), DIMENSION(kdim) ::   zsw_nir_net                        !! Packed radiation from NIR band [W/m^2]
!!$ TR    REAL(wp), DIMENSION(kdim) ::   zsw_nir_frac_diffuse               !! Packed fraction of diffuse NIR radiation
!!$ TR    REAL(wp), DIMENSION(kdim) ::   zsw_par_down                       !! Packed surface downward PAR [W/m^2]
!!$ TR    REAL(wp), DIMENSION(kdim) ::   zsw_par_frac_diffuse               !! Packed fraction of diffuse PAR
    REAL(wp), DIMENSION(kdim) ::   zpressure                          !! Packed surface pressure
    REAL(wp), DIMENSION(kdim) ::   zcdrag
    REAL(wp), DIMENSION(kdim) ::   zetAcoef
    REAL(wp), DIMENSION(kdim) ::   zetBcoef
    REAL(wp), DIMENSION(kdim) ::   zeqAcoef
    REAL(wp), DIMENSION(kdim) ::   zeqBcoef
    !!$ TR inout for JSBACH testing 
    REAL(wp), DIMENSION(kdim) ::   zsurface_temperature
    REAL(wp), DIMENSION(kdim) ::   zsurface_temperature_old
    REAL(wp), DIMENSION(kdim) ::   zsurface_temperature_rad
    REAL(wp), DIMENSION(kdim) ::   zevapotranspiration
    REAL(wp), DIMENSION(kdim) ::   zc_soil_temperature1
    REAL(wp), DIMENSION(kdim) ::   zc_soil_temperature2
    REAL(wp), DIMENSION(kdim) ::   zc_soil_temperature3
    REAL(wp), DIMENSION(kdim) ::   zc_soil_temperature4
    REAL(wp), DIMENSION(kdim) ::   zc_soil_temperature5
    REAL(wp), DIMENSION(kdim) ::   zd_soil_temperature1
    REAL(wp), DIMENSION(kdim) ::   zd_soil_temperature2
    REAL(wp), DIMENSION(kdim) ::   zd_soil_temperature3
    REAL(wp), DIMENSION(kdim) ::   zd_soil_temperature4
    REAL(wp), DIMENSION(kdim) ::   zd_soil_temperature5
    REAL(wp), DIMENSION(kdim) ::   zsoil_temperature1
    REAL(wp), DIMENSION(kdim) ::   zsoil_temperature2
    REAL(wp), DIMENSION(kdim) ::   zsoil_temperature3
    REAL(wp), DIMENSION(kdim) ::   zsoil_temperature4
    REAL(wp), DIMENSION(kdim) ::   zsoil_temperature5
    REAL(wp), DIMENSION(kdim) ::   zheat_capacity
    REAL(wp), DIMENSION(kdim) ::   zground_heat_flux
    REAL(wp), DIMENSION(kdim) ::   zswnet
    REAL(wp), DIMENSION(kdim) ::   ztime_steps_soil

    REAL(wp), DIMENSION(kdim) ::   zsky_view_factor

    !! Other local declarations
    LOGICAL,  DIMENSION(kdim)  ::   mask                               !! Land mask or all true's
    LOGICAL, ALLOCATABLE, DIMENSION(:,:) :: is_glacier
    REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: albedo_echam5
    REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: canopy_snow_fract

    !! local declarations for testing
    REAL(wp), ALLOCATABLE, DIMENSION(:,:)  :: canopy_conductance
!!$ TR    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)  :: root_depth
    REAL(wp), ALLOCATABLE, DIMENSION(:,:)  :: lai
    REAL(wp), ALLOCATABLE, DIMENSION(:)  :: zwind10
!!$ TR   REAL(wp), ALLOCATABLE, DIMENSION(:)  :: swnet

    INTEGER :: nidx, ntiles !!$ TR, nsoil, isoil
    INTEGER :: itile

!!$ TR    INTEGER                 ::   i, kidx0, kidx1, itile

!!$ TR    IF (PRESENT(mask_land)) THEN
!!$ TR       IF (COUNT(mask_land(1:kdim)) /= kland) THEN
!!$ TR          CALL message('jsbach_interface ','kland, COUNT(mask_land) = '// &
!!$ TR                            int2string(kland)//', '//int2string(COUNT(mask_land(1:kdim))))
!!$ TR          CALL finish('jsbach_interface ','wrong mask_land')
!!$ TR       END IF
!!$ TR    END IF

!!$ TR    CALL update_current_call(kdim, kland, kblock=kblock, mask=mask_land)

    nidx = kland
    ntiles = 1
!!$ TR    nsoil = 1

!!$ TR    kidx0 = kstart  ! kstart by use from mo_jsbach_grid
!!$ TR    kidx1 = kend    ! kend by use from mo_jsbach_grid

!!$ TR    IF (PRESENT(mask_land)) THEN
!!$ TR       mask = mask_land
!!$ TR   ELSE
       mask = .TRUE.
!!$ TR    ENDIF
!!$ TR    csat(:) = 1._wp
!!$ TR    cair(:) = 1._wp

    !Generate local packed forcing arrays for each domain (processor)
    IF (PRESENT(wind))         zwind       = PACK(wind, MASK=mask)
!!$ TR    IF (PRESENT(wind10))       zwind10     = PACK(wind10, MASK=mask)
    IF (PRESENT(temp_air)) ztemp_air    = PACK(temp_air, MASK=mask)
    IF (PRESENT(qair)) zqair        = PACK(qair, MASK=mask)
    IF (PRESENT(precip_rain))  zprecip_rain = PACK(precip_rain, MASK=mask)
    IF (PRESENT(precip_snow))  zprecip_snow = PACK(precip_snow, MASK=mask)
    IF (PRESENT(lwdown))       zlwdown      = PACK(lwdown, MASK=mask)
    IF (PRESENT(swdown))       zswdown      = PACK(swdown, MASK=mask)
!!$ TR    IF (PRESENT(sw_vis_net))   zsw_vis_net = PACK(sw_vis_net, MASK=mask)
!!$ TR    IF (PRESENT(sw_vis_frac_diffuse)) zsw_vis_frac_diffuse = PACK(sw_vis_frac_diffuse, MASK=mask)
!!$ TR    IF (PRESENT(sw_nir_net))   zsw_nir_net = PACK(sw_nir_net, MASK=mask)
!!$ TR    IF (PRESENT(sw_nir_frac_diffuse)) zsw_nir_frac_diffuse = PACK(sw_nir_frac_diffuse, MASK=mask)
!!$ TR    IF (PRESENT(sw_par_down))  zsw_par_down = PACK(sw_par_down, MASK=mask)
!!$ TR    IF (PRESENT(sw_par_frac_diffuse)) zsw_par_frac_diffuse = PACK(sw_par_frac_diffuse, MASK=mask)
    IF (PRESENT(pressure)) zpressure    = PACK(pressure, MASK=mask)
!!$ TR    IF (PRESENT(czenith))      cos_zenith(kidx0:kidx1) = PACK(czenith, MASK=mask)
!!$ TR    IF (PRESENT(CO2_concentration)) zCO2         = PACK(CO2_concentration, MASK=mask)
    IF (PRESENT(cdrag))        zcdrag       = PACK(cdrag, MASK=mask)
    IF (PRESENT(etAcoef))      zetAcoef     = PACK(etAcoef, MASK=mask)
    IF (PRESENT(etBcoef))      zetBcoef     = PACK(etBcoef, MASK=mask)
    IF (PRESENT(eqAcoef))      zeqAcoef     = PACK(eqAcoef, MASK=mask)
    IF (PRESENT(eqBcoef))      zeqBcoef     = PACK(eqBcoef, MASK=mask)
    !! inout for testing
    IF (PRESENT(surface_temperature)) zsurface_temperature = &
       PACK(surface_temperature, MASK=mask)
    IF (PRESENT(surface_temperature_old)) zsurface_temperature_old = &
       PACK(surface_temperature_old, MASK=mask)
    IF (PRESENT(surface_temperature_rad)) zsurface_temperature_rad = &
       PACK(surface_temperature_rad, MASK=mask)
    IF (PRESENT(evapotranspiration)) zevapotranspiration = &
       PACK(evapotranspiration, MASK=mask)
    IF (PRESENT(c_soil_temperature1)) zc_soil_temperature1 = &
       PACK(c_soil_temperature1, MASK=mask)
    IF (PRESENT(c_soil_temperature2)) zc_soil_temperature2 = &
       PACK(c_soil_temperature2, MASK=mask)
    IF (PRESENT(c_soil_temperature3)) zc_soil_temperature3 = &
       PACK(c_soil_temperature3, MASK=mask)
    IF (PRESENT(c_soil_temperature4)) zc_soil_temperature4 = &
       PACK(c_soil_temperature4, MASK=mask)
    IF (PRESENT(c_soil_temperature5)) zc_soil_temperature5 = &
       PACK(c_soil_temperature5, MASK=mask)
    IF (PRESENT(d_soil_temperature1)) zd_soil_temperature1 = &
       PACK(d_soil_temperature1, MASK=mask)
    IF (PRESENT(d_soil_temperature2)) zd_soil_temperature2 = &
       PACK(d_soil_temperature2, MASK=mask)
    IF (PRESENT(d_soil_temperature3)) zd_soil_temperature3 = &
       PACK(d_soil_temperature3, MASK=mask)
    IF (PRESENT(d_soil_temperature4)) zd_soil_temperature4 = &
       PACK(d_soil_temperature4, MASK=mask)
    IF (PRESENT(d_soil_temperature5)) zd_soil_temperature5 = &
       PACK(d_soil_temperature5, MASK=mask)
    IF (PRESENT(soil_temperature1)) zsoil_temperature1 = PACK(soil_temperature1, MASK=mask)
    IF (PRESENT(soil_temperature2)) zsoil_temperature2 = PACK(soil_temperature2, MASK=mask)
    IF (PRESENT(soil_temperature3)) zsoil_temperature3 = PACK(soil_temperature3, MASK=mask)
    IF (PRESENT(soil_temperature4)) zsoil_temperature4 = PACK(soil_temperature4, MASK=mask)
    IF (PRESENT(soil_temperature5)) zsoil_temperature5 = PACK(soil_temperature5, MASK=mask)
    IF (PRESENT(heat_capacity)) zheat_capacity = PACK(heat_capacity, MASK=mask)
    IF (PRESENT(ground_heat_flux)) zground_heat_flux = PACK(ground_heat_flux, MASK=mask)
    IF (PRESENT(swnet)) zswnet = PACK(swnet, MASK=mask)
    IF (PRESENT(time_steps_soil)) ztime_steps_soil = PACK(time_steps_soil, MASK=mask)

    ALLOCATE(is_glacier(kdim,ntiles))
    ALLOCATE(albedo_echam5(kdim,ntiles))
    ALLOCATE(canopy_snow_fract(kdim,ntiles))
    ALLOCATE(canopy_conductance(kdim,ntiles))
!!$ TR    ALLOCATE(root_depth(nidx,ntiles,nsoil))
    !!$ JSBACH testing: replace kdim by nidx
    ALLOCATE(lai(kdim,ntiles))
    ALLOCATE(zwind10(kdim))
!!$    ALLOCATE(swnet(kdim))

    !! for testing
!!$ TR    canopy_conductance(:,:) = 1._wp               ! set to 1 for numerical reasons
!!$ TR    root_depth(:,:,:) = 1._wp                     ! root depth equals the depth of the soil bucket (as in the current JSBACH)
    lai(:,:) = 5._wp                              ! no leaves, no transpiration
    zwind10(:) = zwind(:) * 0.8_wp

!!$ TR preliminary calculation of albedo
!!$    zsky_view_factor(:) = EXP(-lai(:,1)/2._wp)
!!$    albedo_vis(:) = (zsky_view_factor(:) * albedo_vis_soil(:)) + &
!!$                    ((1._wp - zsky_view_factor(:)) * albedo_vis_canopy(:))
!!$    albedo_nir(:) = (zsky_view_factor(:) * albedo_nir_soil(:)) + &
!!$                    ((1._wp - zsky_view_factor(:)) * albedo_nir_canopy(:))
!!$ TR conversion of swdown in swnet is not accurate and preliminary for testing
    zswnet(:) = zswdown(:) * (1._wp - (albvisdir(:) + albnirdir(:)) / 2._wp)

    DO itile=1,ntiles
       CALL unstressed_canopy_cond_par(lai(1:nidx,itile), zswdown(:) / 2._wp, &
                  canopy_conductance(1:nidx,itile))
    END DO

!!$ TR    ALLOCATE(root_depth(nidx,ntiles,nsoil))    
    CALL update_soil_icon(ntiles, kdim, &
                     delta_time, time_step_len, &
                     ztime_steps_soil, &
                     canopy_conductance(:,:), &
!!$ TR                     root_depth(:,:,:), &
                     lai(:,:), zcdrag, &
                     zetAcoef, zetBcoef, zeqAcoef, zeqBcoef, &
                     ztemp_air, zqair, &
                     zpressure, zwind, zwind10, &
                     zlwdown, zswnet, zprecip_rain, zprecip_snow, &
                     p_echam_zchl, &
                     !! output for testing (hydrology)
                     cair, &
                     csat, &
                     csat_transpiration, &
                     moisture1, &
                     moisture2, &
                     moisture3, &
                     moisture4, &
                     moisture5, &
                     moisture_all, &
                     sat_surface_specific_humidity, &
                     skin_reservoir, &
                     snow_fract, &
                     snow, &
                     snow_canopy, &
                     snow_melt, &
                     snow_acc, &
                     snow_melt_acc, &
                     glacier_runoff_acc, &
                     runoff_acc, &
                     drainage_acc, &
                     canopy_snow_fract, &
                     !! output for testing (energy balance)
                     zsurface_temperature, &
                     zsurface_temperature_old, &
                     zsurface_temperature_rad, &
                     zevapotranspiration, &
                     zc_soil_temperature1, &
                     zc_soil_temperature2, &
                     zc_soil_temperature3, &
                     zc_soil_temperature4, &
                     zc_soil_temperature5, &
                     zd_soil_temperature1, &
                     zd_soil_temperature2, &
                     zd_soil_temperature3, &
                     zd_soil_temperature4, &
                     zd_soil_temperature5, &
                     zsoil_temperature1, &
                     zsoil_temperature2, &
                     zsoil_temperature3, &
                     zsoil_temperature4, &
                     zsoil_temperature5, &
                     zheat_capacity, &
                     zground_heat_flux &
                     )
    zswnet(:) = zlwdown(:)
                                          
!!$ TR                     p_echam_zchl, zzhsoil, ztte_corr, &
!!$ TR                     theLand%Vegetation%canopy_conductance_limited(kidx0:kidx1,:), &
!!$ TR                     zglac_runoff_evap, zsurf_runoff_hd, zdrainage_hd)

    is_glacier(:,:) = .false.

    CALL update_albedo_echam5(nidx, ntiles, is_glacier(:,:), forest_fract(:), &
           zsurface_temperature(:),                            &
           SPREAD(snow_fract(:),DIM=2, NCOPIES=ntiles),        &
           SPREAD(albedo_background(:),DIM=2, NCOPIES=ntiles), &
           canopy_snow_fract  (:,:),                           &
           lai                (:,:),                           &
           albedo_echam5(:,:))

    WHERE (mask_land(:) > 0)
      albvisdir(:) = albedo_echam5(:,1)
      albnirdir(:) = albedo_echam5(:,1)
      albvisdif(:) = albedo_echam5(:,1)
      albnirdif(:) = albedo_echam5(:,1)
    ENDWHERE

    !! for testing
!!$ TR    zhsoil() = 0._wp                              ! dry surface
    
    IF (PRESENT(surface_temperature)) surface_temperature(:) = &
        UNPACK(zsurface_temperature, mask, 0.0_wp)
    IF (PRESENT(surface_temperature_old)) surface_temperature_old(:) = &
        UNPACK(zsurface_temperature_old, mask, 0.0_wp)
    IF (PRESENT(surface_temperature_rad)) surface_temperature_rad(:) = &
        UNPACK(zsurface_temperature_rad, mask, 0.0_wp)
    IF (PRESENT(evapotranspiration)) evapotranspiration(:) = &
        UNPACK(zevapotranspiration, mask, 0.0_wp)
    IF (PRESENT(c_soil_temperature1)) c_soil_temperature1(:) = &
        UNPACK(zc_soil_temperature1, mask, 0.0_wp)
    IF (PRESENT(c_soil_temperature2)) c_soil_temperature2(:) = &
        UNPACK(zc_soil_temperature2, mask, 0.0_wp)
    IF (PRESENT(c_soil_temperature3)) c_soil_temperature3(:) = &
        UNPACK(zc_soil_temperature3, mask, 0.0_wp)
    IF (PRESENT(c_soil_temperature4)) c_soil_temperature4(:) = &
        UNPACK(zc_soil_temperature4, mask, 0.0_wp)
    IF (PRESENT(c_soil_temperature5)) c_soil_temperature5(:) = &
        UNPACK(zc_soil_temperature5, mask, 0.0_wp)
    IF (PRESENT(d_soil_temperature1)) d_soil_temperature1(:) = &
        UNPACK(zd_soil_temperature1, mask, 0.0_wp)
    IF (PRESENT(d_soil_temperature2)) d_soil_temperature2(:) = &
        UNPACK(zd_soil_temperature2, mask, 0.0_wp)
    IF (PRESENT(d_soil_temperature3)) d_soil_temperature3(:) = &
        UNPACK(zd_soil_temperature3, mask, 0.0_wp)
    IF (PRESENT(d_soil_temperature4)) d_soil_temperature4(:) = &
        UNPACK(zd_soil_temperature4, mask, 0.0_wp)
    IF (PRESENT(d_soil_temperature5)) d_soil_temperature5(:) = &
        UNPACK(zd_soil_temperature5, mask, 0.0_wp)
    IF (PRESENT(soil_temperature1)) soil_temperature1(:) = UNPACK(zsoil_temperature1, mask, 0.0_wp)
    IF (PRESENT(soil_temperature2)) soil_temperature2(:) = UNPACK(zsoil_temperature2, mask, 0.0_wp)
    IF (PRESENT(soil_temperature3)) soil_temperature3(:) = UNPACK(zsoil_temperature3, mask, 0.0_wp)
    IF (PRESENT(soil_temperature4)) soil_temperature4(:) = UNPACK(zsoil_temperature4, mask, 0.0_wp)
    IF (PRESENT(soil_temperature5)) soil_temperature5(:) = UNPACK(zsoil_temperature5, mask, 0.0_wp)
    IF (PRESENT(heat_capacity)) heat_capacity(:) = UNPACK(zheat_capacity, mask, 0.0_wp)
    IF (PRESENT(ground_heat_flux)) ground_heat_flux(:) = UNPACK(zground_heat_flux, mask, 0.0_wp)
    IF (PRESENT(swnet)) swnet(:) = UNPACK(zswnet, mask, 0.0_wp)
    IF (PRESENT(time_steps_soil)) time_steps_soil(:) = UNPACK(ztime_steps_soil, mask, 0.0_wp)

    !! deallocate variables for testing
    DEALLOCATE(is_glacier)
    DEALLOCATE(albedo_echam5)
    DEALLOCATE(canopy_snow_fract)
    DEALLOCATE(canopy_conductance)
!!$ TR    DEALLOCATE(root_depth)
    DEALLOCATE(lai)
    DEALLOCATE(zwind10)
!!$ TR    DEALLOCATE(swnet)

   END SUBROUTINE jsbach_inter_1d

END MODULE mo_jsbach_interface_icon
