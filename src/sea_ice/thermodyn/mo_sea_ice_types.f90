!>
!! Provide an implementation of the sea-ice model.
!!
!! Provide an implementation of the parameters of the surface module (sea ice)
!! used between the atmopshere and the hydrostatic ocean model.
!!
!! @author
!!
!! @par Revision History
!!
!! @par Copyright
!! 2002-2007 by DWD and MPI-M
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
MODULE mo_sea_ice_types
  USE mo_kind,                ONLY: wp
  USE mo_math_utilities,      ONLY: t_cartesian_coordinates

  IMPLICIT NONE
  PRIVATE


  ! Definition of forcing types
  ! public types
  PUBLIC  :: t_sea_ice
  PUBLIC  :: t_sea_ice_acc
  PUBLIC  :: t_sfc_flx
  PUBLIC  :: t_atmos_fluxes
  PUBLIC  :: t_atmos_for_ocean
  PUBLIC  :: t_ptr2d



  TYPE t_ptr2d
    REAL(wp),POINTER :: p(:,:)  ! pointer to 2D (spatial) array
  END TYPE t_ptr2d
  !------  Definition of surface flux type---------------------
  TYPE t_sfc_flx

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
      &  forc_hfrelax              (:,:), & ! diagnosed surface heat flux due to relaxation (EMPTY)     [m/s]
      &  forc_fwrelax              (:,:), & ! diagnosed surface freshwater flux due to relaxation       [m/s]
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
      &  forc_hfrelax_acc     (:,:),     &
      &  forc_fwrelax_acc               (:,:),  &
      &  data_surfRelax_Temp_acc          (:,:),   &
      &  data_surfRelax_Salt_acc(:,:)

    TYPE(t_cartesian_coordinates), & ! wind forcing with cartesian vector, located at cell centers
      & ALLOCATABLE :: topBoundCond_windStress_cc(:,:)

    TYPE(t_ptr2d),ALLOCATABLE :: tracer_ptr(:)  !< pointer array: one pointer for each tracer
  END TYPE t_sfc_flx

  ! global type variables
  TYPE(t_sfc_flx), PUBLIC, TARGET :: v_sfc_flx

  !------  Definition of representation of atm state in ocean model---
  !
  !representation of atmosphere in ocean model. Data are coming either from
  !atmosphere model or from file. These fields are transformed via bulk fomulas
  !into atmospheric fluxes, the fluxes are then used to set the oceans surface
  !boundary conditions
  ! On cells
  TYPE t_atmos_for_ocean

    REAL(wp), ALLOCATABLE :: &
      & tafo(:,:),             &  ! 2 m air temperature                              [C]
      & ftdew(:,:),            &  ! 2 m dew-point temperature                        [K]
      & fclou(:,:),            &  ! Fractional cloud cover
      & fu10(:,:) ,            &  ! 10 m wind speed                                  [m/s]
      & fswr(:,:),             &  ! Incoming surface solar radiation                 [W/m]
      & pao(:,:),              &  ! Surface atmospheric pressure                     [hPa]
      & u(:,:),                &  ! wind in reference height                         [m/s]
      & v(:,:),                &
      & precip(:,:),           &  ! precipitation rate                               [m/s]
      & evap  (:,:),           &  ! evaporation   rate                               [m/s]
      & runoff(:,:),           &  ! river runoff  rate                               [m/s]
      & topBoundCond_windStress_u(:,:), &
      & topBoundCond_windStress_v(:,:), &
      & FrshFlux_Precipitation(:,:),    &
      & FrshFlux_Evaporation(:,:),    &
      & FrshFlux_TotalOcean(:,:),    &
      & FrshFlux_Runoff(:,:),           &
      & data_surfRelax_Temp(:,:)

  END TYPE t_atmos_for_ocean



  !------  Definition of forcing---------------------
  TYPE t_atmos_fluxes

    REAL(wp), ALLOCATABLE ::   &
      & sens    (:,:,:),           & ! Sensible heat flux at ice surface           [W/m2]
      & lat     (:,:,:),           & ! Latent heat flux at ice surface             [W/m2]
      & LWout   (:,:,:),           & ! outgoing LW radiation flux at ice surface   [W/m2]
      & LWnet   (:,:,:),           & ! net LW radiation flux at ice surface        [W/m2]
      & SWnet   (:,:,:),           & ! net SW radiation flux over ice              [W/m2]
      & bot     (:,:,:),           & ! Ocean heat flux at ice bottom               [W/m2]
      & dsensdT (:,:,:),           & ! d sensible Flux / d T_surf                  [W/m2/K]
      & dlatdT  (:,:,:),           & ! d latent Flux / d T_surf                    [W/m2/K]
      & dLWdT   (:,:,:),           & ! d radiation Flux / d T_surf                 [W/m2/K]
      & stress_x(:,:),             & ! Wind stress at the ice surface              [Pa]
      & stress_y(:,:)                ! Wind stress at the ice surface              [Pa]

    REAL(wp), ALLOCATABLE ::   &
      & rprecw (:,:),             & ! liquid precipitation rate                   [m/s]
      & rpreci (:,:),             & ! solid  precipitation rate                   [m/s]
      & sensw  (:,:),             & ! Sensible heat flux over water               [W/m2]
      & latw   (:,:),             & ! Latent heat flux over water                 [W/m2]
      & LWoutw (:,:),             & ! outgoing LW radiation flux over water       [W/m2]
      & LWnetw (:,:),             & ! net LW radiation flux over water            [W/m2]
      & SWnetw (:,:),             & ! net SW radiation flux over water            [W/m2]
      & LWin   (:,:),             & ! incoming LW radiation flux                  [W/m2]
      & stress_xw(:,:),           & ! Wind stress at the ocean surface            [Pa]
      & stress_yw(:,:)              ! Wind stress at the ocean surface            [Pa]

! Albedos
    REAL(wp), POINTER ::     &
      & albvisdir (:,:,:),      & ! VIS direct/paralell (ice)
      & albvisdif (:,:,:),      & ! VIS diffuse (ice)
      & albnirdir (:,:,:),      & ! NIR direct/paralell (ice)
      & albnirdif (:,:,:),      & ! NIR diffuse (ice)
      & albvisdirw(:,:),        & ! VIS direct/paralell (ocean)
      & albvisdifw(:,:),        & ! VIS diffuse (ocean)
      & albnirdirw(:,:),        & ! NIR direct/paralell (ocean)
      & albnirdifw(:,:)           ! NIR diffuse (ocean)

    INTEGER ::     counter

    REAL(wp), ALLOCATABLE ::   &
      &  topBoundCond_windStress_u(:,:),     & ! forcing of zonal component of velocity equation,
      &  topBoundCond_windStress_v(:,:),     & ! forcing of meridional component of velocity equation,
      &  HeatFlux_ShortWave       (:,:),     & ! surface short wave heat flux                              [W/m2]
      &  HeatFlux_LongWave        (:,:),     & ! surface long wave heat flux                               [W/m2]
      &  HeatFlux_Sensible        (:,:),     & ! surface sensible heat flux                                [W/m2]
      &  HeatFlux_Latent          (:,:),     & ! surface latent heat flux                                  [W/m2]
      &  FrshFlux_Precipitation      (:,:),     & ! total precipitation flux                                  [m/s]
      &  FrshFlux_SnowFall        (:,:),     & ! total snow flux                                           [m/s]
      &  FrshFlux_Evaporation        (:,:),     & ! evaporation flux                                          [m/s]
      &  FrshFlux_Runoff      (:,:)!     & ! river runoff flux                                         [m/s]

    TYPE(t_cartesian_coordinates), & ! wind forcing with cartesian vector, located at cell centers
      & ALLOCATABLE :: topBoundCond_windStress_cc(:,:)

  END TYPE t_atmos_fluxes

  TYPE t_sea_ice_acc
  ! accumulated fields of the sea-ice state
    REAL(wp), POINTER :: &
      & hi         (:,:,:)       ,   & ! Ice thickness                                 [m]
      & hs         (:,:,:)       ,   & ! Snow thickness                                [m]
      & conc       (:,:,:)             ! ice concentration in each ice class
    REAL(wp), POINTER :: &
      & u(:,:)          ,      & ! Zonal velocity on cell centre (diagnostic)    [m/s]
      & v(:,:)                   ! Meridional velocity on cell centre (diagn.)   [m/s]
  END TYPE t_sea_ice_acc
  TYPE t_sea_ice

  ! The description of the sea-ice state, defined on cell-centers
  ! dimension: (nproma, nblks_c)

    REAL(wp), POINTER :: &
      & alb        (:,:,:)       ,   & ! Albedo of snow-ice system
      & Tsurf      (:,:,:)       ,   & ! Surface temperature                           [C]
      & T1         (:,:,:)       ,   & ! Temperature upper layer                       [C]
      & T2         (:,:,:)       ,   & ! Temperature lower layer                       [C]
      & E1         (:,:,:)       ,   & ! Energy content upper layer                    [Jm/kg]
      & E2         (:,:,:)       ,   & ! Energy content lower layer                    [Jm/kg]
      & vol        (:,:,:)       ,   & ! Ice volume                                    [m^3]
      & vols       (:,:,:)       ,   & ! Snow volume                                   [m^3]
      & hi         (:,:,:)       ,   & ! Ice thickness                                 [m]
      & hs         (:,:,:)       ,   & ! Snow thickness                                [m]
      & hiold      (:,:,:)       ,   & ! Ice thickness at previous time step           [m]
      & hsold      (:,:,:)       ,   & ! Snow thickness at previous time step          [m]
      & Qtop       (:,:,:)       ,   & ! Energy flux available for surface melting     [W/m2]
      & Qbot       (:,:,:)       ,   & ! Energy flux at ice-ocean interface            [W/m2]
      & heatocei   (:,:,:)       ,   & ! Energy to ocean when all ice is melted        [J]
      & snow_to_ice(:,:,:)       ,   & ! amount of snow that is transformed to ice     [m]
      & surfmelt   (:,:,:)       ,   & ! surface melt water running into ocean         [m]
      & surfmeltT  (:,:,:)       ,   & ! Mean temperature of surface melt water        [C]
      & evapwi     (:,:,:)       ,   & ! amount of evaporated water if no ice left     [kg/m2]
      & conc       (:,:,:)             ! ice concentration in each ice class

    REAL(wp), POINTER :: &
      & u_prog(:,:)     ,      & ! Zonal velocity (prognostic, rotated grid)     [m/s]
      & v_prog(:,:)     ,      & ! Meridional velocity (prognostic, rotated)     [m/s]
      & u(:,:)          ,      & ! Zonal velocity on cell centre (diagnostic)    [m/s]
      & v(:,:)          ,      & ! Meridional velocity on cell centre (diagn.)   [m/s]
      & vn_e(:,:)       ,      & ! Edge normal velocity (diagnostic)             [m/s]
      & concSum(:,:)    ,      & ! Total ice concentration within a grid cell
      & newice(:,:)     ,      & ! New ice growth in open water                  [m]
      & zUnderIce(:,:)           ! water in upper ocean grid cell below ice      [m]

    INTEGER ::  kice           ! Number of ice-thickness classes

    REAL(wp), ALLOCATABLE ::  hi_lim(:)   ! Thickness limits

    TYPE(t_sea_ice_acc)   :: acc

  END TYPE t_sea_ice



  ! global type variables
  TYPE(t_sea_ice),PUBLIC, SAVE, TARGET :: v_sea_ice



END MODULE mo_sea_ice_types
