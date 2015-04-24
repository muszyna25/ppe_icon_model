!>
!! Provide an implementation of the types for the surface module
!!
!! Provide an implementation of the types for the surface module
!! used between the atmosphere and the hydrostatic ocean model.
!!
!! @author Stephan Lorenz, MPI
!!
!! @par Revision History
!!  Initial release by Stephan Lorenz, MPI-M (2015-04)
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
  USE mo_math_utilities,      ONLY: t_cartesian_coordinates

  IMPLICIT NONE
  PRIVATE


  ! Definition of forcing types for ocean surface module
  ! public types
  PUBLIC  :: t_ocean_surface_fluxes
  PUBLIC  :: t_ptr2d



  TYPE t_ptr2d
    REAL(wp),POINTER :: p(:,:)  ! pointer to 2D (spatial) array
  END TYPE t_ptr2d

  !------  Definition of surface flux type---------------------

  ! These fluxes will successively replace the fluxes defined in type t_sfc_flx,
  ! they are at the end exclusively used by module mo_ocean_surface

  TYPE t_ocean_surface_fluxes

    ! The forcing is specified as fluxes at the air-sea interface defined on cell-centers
    ! dimension: (nproma, nblks_c)
    REAL(wp), POINTER ::   &
      &  topBoundCond_windStress_u (:,:), & ! forcing of zonal component of velocity equation           [Pa]
      &  topBoundCond_windStress_v (:,:), & ! forcing of meridional component of velocity equation      [Pa]
      &  HeatFlux_ShortWave        (:,:), & ! surface short wave heat flux                              [W/m2]
      &  HeatFlux_LongWave         (:,:), & ! surface long wave heat flux                               [W/m2]
      &  HeatFlux_Sensible         (:,:), & ! surface sensible heat flux                                [W/m2]
      &  HeatFlux_Latent           (:,:), & ! surface latent heat flux                                  [W/m2]
      &  HeatFlux_Total            (:,:), & ! sum of forcing surface heat flux                          [W/m2]
      &  FrshFlux_Precipitation    (:,:), & ! total precipitation flux                                  [m/s]
      &  FrshFlux_SnowFall         (:,:), & ! total snow flux                                           [m/s]
      &  FrshFlux_Evaporation      (:,:), & ! evaporation flux                                          [m/s]
      &  FrshFlux_Runoff           (:,:), & ! river runoff flux                                         [m/s]
      &  FrshFlux_TotalSalt        (:,:), & ! sum of forcing surface freshwater flux from BC            [m/s]
      &  FrshFlux_TotalOcean       (:,:), & ! forcing surface freshwater flux at open ocean             [m/s]
      &  FrshFlux_TotalIce         (:,:), & ! forcing surface freshwater flux under sea ice             [m/s]
      &  FrshFlux_VolumeIce        (:,:), & ! forcing volume flux for height equation under sea ice     [m/s]
      &  FrshFlux_VolumeTotal      (:,:), & ! sum of forcing volume flux including relaxation           [m/s]
      &  topBoundCond_Temp_vdiff   (:,:), & ! forcing of temperature in vertical diffusion equation     [K*m/s]
      &  topBoundCond_Salt_vdiff   (:,:), & ! forcing of salinity in vertical diffusion equation        [psu*m/s]
      &  data_surfRelax_Temp(:,:),        & ! contains data to which temperature is relaxed             [K]
      &  data_surfRelax_Salt(:,:),        & ! contains data to which salinity is relaxed                [psu]
      &  HeatFlux_Relax            (:,:), & ! surface heat flux due to relaxation                       [W/m2]
      &  FrshFlux_Relax            (:,:), & ! surface freshwater flux due to relaxation                 [m/s]
      &  TempFlux_Relax            (:,:), & ! temperature tracer flux due to relaxation                 [K/s]
      &  SaltFlux_Relax            (:,:), & ! salinity tracer flux due to relaxation                    [psu/s]
      !
      !  accumulations variables - comments see above
      &  topBoundCond_windStress_u_acc  (:,:),  &
      &  topBoundCond_windStress_v_acc  (:,:),  &
      &  HeatFlux_ShortWave_acc         (:,:),  &
      &  HeatFlux_LongWave_acc          (:,:),  &
      &  HeatFlux_Sensible_acc          (:,:),  &
      &  HeatFlux_Latent_acc            (:,:),  &
      &  HeatFlux_Total_acc             (:,:),  &
      &  FrshFlux_Precipitation_acc     (:,:),  &
      &  FrshFlux_SnowFall_acc          (:,:),  &
      &  FrshFlux_Evaporation_acc       (:,:),  &
      &  FrshFlux_Runoff_acc            (:,:),  &
      &  FrshFlux_TotalSalt_acc         (:,:),  &
      &  FrshFlux_TotalOcean_acc        (:,:),  &
      &  FrshFlux_TotalIce_acc          (:,:),  &
      &  FrshFlux_VolumeIce_acc         (:,:),  &
      &  FrshFlux_VolumeTotal_acc       (:,:),  &
      &  topBoundCond_Temp_vdiff_acc    (:,:),  &
      &  topBoundCond_Salt_vdiff_acc    (:,:),  &
      &  HeatFlux_Relax_acc             (:,:),  &
      &  FrshFlux_Relax_acc             (:,:),  &
      &  TempFlux_Relax_acc             (:,:),  &
      &  SaltFlux_Relax_acc             (:,:),  &
      &  data_surfRelax_Temp_acc        (:,:),  &
      &  data_surfRelax_Salt_acc        (:,:),  &
      !
      !
      &  cellThicknessUnderIce          (:,:)

    TYPE(t_cartesian_coordinates), & ! wind forcing with cartesian vector, located at cell centers
      & ALLOCATABLE :: topBoundCond_windStress_cc(:,:)

    TYPE(t_ptr2d),ALLOCATABLE :: tracer_ptr(:)  !< pointer array: one pointer for each tracer
  END TYPE t_ocean_surface_fluxes

  ! global type variables
  TYPE(t_ocean_surface_fluxes), PUBLIC, TARGET :: v_ocean_sfc_flx

END MODULE mo_ocean_surface_types

