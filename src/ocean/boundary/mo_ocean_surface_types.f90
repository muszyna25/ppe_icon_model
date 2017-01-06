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
  PUBLIC  :: t_ocean_surface


  !------  Definition of surface flux type---------------------

  ! These fluxes will successively replace the fluxes defined in type t_sfc_flx,
  ! they are at the end exclusively used by module mo_ocean_surface

  TYPE t_ocean_surface

    ! The forcing is specified as fluxes at the air-sea interface defined on cell-centers
    ! dimension: (nproma, nblks_c)
    REAL(wp), POINTER ::   &
      &  TopBC_WindStress_u        (:,:), & ! forcing of zonal component of velocity equation           [Pa]
      &  TopBC_WindStress_v        (:,:), & ! forcing of meridional component of velocity equation      [Pa]
      &  Wind_Speed_10m            (:,:), & ! wind speed in 10m height                                  [m/s]
      &  HeatFlux_Total            (:,:), & ! sum of forcing surface heat flux                          [W/m2]
      &  HeatFlux_Shortwave        (:,:), & ! shortwave heat flux for penetration into deeper layers    [W/m2]
      &  FrshFlux_TotalIce         (:,:), & ! forcing surface freshwater flux due to sea ice change     [m/s]
      &  FrshFlux_VolumeTotal      (:,:), & ! sum of forcing volume flux including relaxation           [m/s]
      &  SST                       (:,:), & ! sea surface temperature                                   [C]
      &  SSS                       (:,:), & ! sea surface salinity                                      [psu]
      &  data_surfRelax_Temp(:,:),        & ! contains data to which temperature is relaxed             [C]
      &  data_surfRelax_Salt(:,:),        & ! contains data to which salinity is relaxed                [psu]
      &  HeatFlux_Relax            (:,:), & ! surface heat flux due to relaxation                       [W/m2]
      &  FrshFlux_Relax            (:,:), & ! surface freshwater flux due to relaxation                 [m/s]
      &  cellThicknessUnderIce     (:,:), & ! thickness of freeboard, open water below ice              [m]
      !
      !  accumulation variables - comments see above
      &  TopBC_WindStress_u_acc         (:,:),  &
      &  TopBC_WindStress_v_acc         (:,:),  &
      &  Wind_Speed_10m_acc             (:,:),  &
      &  HeatFlux_Total_acc             (:,:),  &
      &  HeatFlux_Shortwave_acc         (:,:),  &
      &  FrshFlux_TotalIce_acc          (:,:),  &
      &  FrshFlux_VolumeTotal_acc       (:,:),  &
      &  SST_acc                        (:,:),  &
      &  SSS_acc                        (:,:),  &
      &  HeatFlux_Relax_acc             (:,:),  &
      &  FrshFlux_Relax_acc             (:,:),  &
      &  cellThicknessUnderIce_acc      (:,:)

    TYPE(t_cartesian_coordinates), & ! wind forcing with cartesian vector, located at cell centers
      & ALLOCATABLE :: TopBC_WindStress_cc(:,:)

  ! TYPE(t_ptr2d),ALLOCATABLE :: tracer_ptr(:)  !< pointer array: one pointer for each tracer
  END TYPE t_ocean_surface

  ! global type variables
  TYPE(t_ocean_surface), PUBLIC, TARGET :: v_oce_sfc

END MODULE mo_ocean_surface_types

