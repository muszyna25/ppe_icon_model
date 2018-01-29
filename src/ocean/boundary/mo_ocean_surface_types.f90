!>
!! This module provides definition of ocean surface module types that describe:
!!  1) ocean surface fluxes : t_ocean_surface
!!  2) OMIP fluxes          : t_atmos_fluxes
!! ----------------------------------------------------------------------------------------
!!
!!
!! @author Stephan Lorenz, MPI
!!
!! @par Revision History
!!  Initial release by Stephan Lorenz, MPI-M (2015-04)
!!  Modified        by Vladimir Lapin, MPI-M (2017-04)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_ocean_surface_types
  USE mo_kind,                ONLY: wp
  USE mo_math_types,          ONLY: t_cartesian_coordinates

  IMPLICIT NONE
  PRIVATE

  PUBLIC  :: t_ocean_surface
  PUBLIC  :: t_atmos_for_ocean

! ---------------------------------------------------------------------------------------

  ! --------------------
  ! Ocean surface types
  ! --------------------
  !
  TYPE t_ocean_surface

    ! The forcing is specified as fluxes at the air-sea interface defined on cell-centers
    ! dimension: (nproma, nblks_c)
    REAL(wp), POINTER ::   &
      &  Wind_Speed_10m            (:,:), & ! wind speed in 10m height                                  [m/s]
      &  TopBC_WindStress_u        (:,:), & ! forcing of zonal component of velocity equation           [Pa]
      &  TopBC_WindStress_v        (:,:), & ! forcing of meridional component of velocity equation      [Pa]
      &  SST                       (:,:), & ! sea surface temperature                                   [C]
      &  SSS                       (:,:), & ! sea surface salinity                                      [psu]
      &  cellThicknessUnderIce     (:,:), & ! thickness of freeboard, open water below ice              [m]

      ! heat fluxes
      &  HeatFlux_Total            (:,:), & ! sum of forcing surface heat flux                          [W/m2]
      &  HeatFlux_Shortwave        (:,:), & ! shortwave heat flux for penetration into deeper layers    [W/m2]
      ! auxillary heat fluxes
      &  HeatFlux_LongWave         (:,:), & ! surface long wave heat flux                               [W/m2]
      &  HeatFlux_Sensible         (:,:), & ! surface sensible heat flux                                [W/m2]
      &  HeatFlux_Latent           (:,:), & ! surface latent heat flux                                  [W/m2]

      ! freshwater fluxes
      &  FrshFlux_TotalIce         (:,:), & ! forcing surface freshwater flux due to sea ice change     [m/s]
      &  FrshFlux_VolumeTotal      (:,:), & ! sum of forcing volume flux including relaxation           [m/s]
      &  FrshFlux_IceSalt          (:,:), & ! salt volume flux due to sea ice change                    [psu*m/s]
      ! auxillary freshwater fluxes
      &  FrshFlux_Precipitation    (:,:), & ! total precipitation flux                                  [m/s]
      &  FrshFlux_SnowFall         (:,:), & ! total snow flux                                           [m/s]
      &  FrshFlux_Evaporation      (:,:), & ! evaporation flux                                          [m/s]
      &  FrshFlux_Runoff           (:,:), & ! river runoff flux                                         [m/s]
      &  FrshFlux_TotalSalt        (:,:), & ! sum of forcing surface freshwater flux from BC            [m/s]
      &  FrshFlux_TotalOcean       (:,:), & ! forcing surface freshwater flux at open ocean             [m/s]
      &  FrshFlux_VolumeIce        (:,:), & ! forcing volume flux for height equation under sea ice     [m/s]

      ! relaxaton
      &  data_surfRelax_Temp       (:,:), & ! contains data to which temperature is relaxed             [C]
      &  data_surfRelax_Salt       (:,:), & ! contains data to which salinity is relaxed                [psu]
      &  TopBC_Temp_vdiff          (:,:), & ! forcing of temperature in vertical diffusion equation     [C*m/s]
      &  TopBC_Salt_vdiff          (:,:), & ! forcing of salinity in vertical diffusion equation        [psu*m/s]
      &  TempFlux_Relax            (:,:), & ! temperature tracer flux due to relaxation                 [C/s]
      &  SaltFlux_Relax            (:,:), & ! salinity tracer flux due to relaxation                    [psu/s]
      &  HeatFlux_Relax            (:,:), & ! surface heat flux due to relaxation                       [W/m2]
      &  FrshFlux_Relax            (:,:), & ! surface freshwater flux due to relaxation                 [m/s]

      !  accumulation variables - comments see above
      &  Wind_Speed_10m_acc             (:,:),  &
      &  TopBC_WindStress_u_acc         (:,:),  &
      &  TopBC_WindStress_v_acc         (:,:),  &
!      &  cellThicknessUnderIce_acc      (:,:),  &

      &  HeatFlux_Total_acc             (:,:),  &
      &  HeatFlux_Shortwave_acc         (:,:),  &
      &  HeatFlux_LongWave_acc          (:,:),  &
      &  HeatFlux_Sensible_acc          (:,:),  &
      &  HeatFlux_Latent_acc            (:,:),  &

      &  FrshFlux_TotalIce_acc          (:,:),  &
      &  FrshFlux_VolumeTotal_acc       (:,:),  &
      &  FrshFlux_Precipitation_acc     (:,:),  &
      &  FrshFlux_SnowFall_acc          (:,:),  &
      &  FrshFlux_Evaporation_acc       (:,:),  &
      &  FrshFlux_Runoff_acc            (:,:),  &
      &  FrshFlux_TotalSalt_acc         (:,:),  &
      &  FrshFlux_TotalOcean_acc        (:,:),  &
      &  FrshFlux_VolumeIce_acc         (:,:),  &

      ! Are these relaxation-related accumulated fields really necessary?
      &  data_surfRelax_Temp_acc        (:,:),  &
      &  data_surfRelax_Salt_acc        (:,:),  &
      &  TopBC_Temp_vdiff_acc           (:,:),  &
      &  TopBC_Salt_vdiff_acc           (:,:),  &
!      &  TempFlux_Relax_acc             (:,:),  &
!      &  SaltFlux_Relax_acc             (:,:),  &
      &  HeatFlux_Relax_acc             (:,:),  &
      &  FrshFlux_Relax_acc             (:,:)

    TYPE(t_cartesian_coordinates), & ! wind forcing with cartesian vector, located at cell centers
      & ALLOCATABLE :: TopBC_WindStress_cc(:,:)

  END TYPE t_ocean_surface

  ! global type variables
  TYPE(t_ocean_surface), PUBLIC, TARGET :: v_oce_sfc

! ---------------------------------------------------------------------------------------

  !---------------------------------------
  !------------  OMIP forcing ------------
  !---------------------------------------
  ! OMIP representation of atmosphere state for driving the ocean model.
  ! These fields are transformed via bulk fomulas into atmospheric fluxes.
  ! The fluxes are then used to set the oceans surface boundary conditions.

  TYPE t_atmos_for_ocean

    ! on cell-centers, dimension: (nproma, nblks_c)
    REAL(wp), POINTER :: &
      & tafo                     (:,:), & ! 2 m air temperature                              [C]
      & ftdew                    (:,:), & ! 2 m dew-point temperature                        [K]
      & fclou                    (:,:), & ! Fractional cloud cover
      & fu10                     (:,:), & ! 10 m wind speed                                  [m/s]
      & fswr                     (:,:), & ! Incoming surface solar radiation                 [W/m]
      & pao                      (:,:), & ! Surface atmospheric pressure                     [hPa]
      & u                        (:,:), & ! wind in reference height                         [m/s]
      & v                        (:,:), &
!      & precip                   (:,:), & ! precipitation rate                               [m/s]
!      & evap                     (:,:), & ! evaporation   rate                               [m/s]
!      & runoff                   (:,:), & ! river runoff  rate                               [m/s]
      & topBoundCond_windStress_u(:,:), &
      & topBoundCond_windStress_v(:,:), &
      & FrshFlux_Precipitation   (:,:), &
      & FrshFlux_Runoff          (:,:), &
      & data_surfRelax_Salt      (:,:), &
      & data_surfRelax_Temp      (:,:)

  END TYPE t_atmos_for_ocean

! ---------------------------------------------------------------------------------------

END MODULE mo_ocean_surface_types

