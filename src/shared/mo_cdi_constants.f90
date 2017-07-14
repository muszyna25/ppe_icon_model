!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
MODULE mo_cdi_constants

  USE mo_cdi, ONLY: ZAXIS_ALTITUDE, ZAXIS_ATMOSPHERE, ZAXIS_CLOUD_BASE, ZAXIS_CLOUD_TOP,        &
    &               ZAXIS_DEPTH_BELOW_LAND, ZAXIS_DEPTH_BELOW_SEA, ZAXIS_GENERIC, ZAXIS_HEIGHT, &
    &               ZAXIS_HYBRID, ZAXIS_HYBRID_HALF, ZAXIS_ISENTROPIC, ZAXIS_ISOTHERM_ZERO,     &
    &               ZAXIS_LAKE_BOTTOM, ZAXIS_MEANSEA, ZAXIS_MIX_LAYER, ZAXIS_PRESSURE,          &
    &               ZAXIS_REFERENCE, ZAXIS_SEDIMENT_BOTTOM_TW, ZAXIS_SURFACE, ZAXIS_TOA


  IMPLICIT NONE

  PUBLIC

  INTEGER, PARAMETER :: GRID_CELL   = 1
  INTEGER, PARAMETER :: GRID_VERTEX = 2
  INTEGER, PARAMETER :: GRID_EDGE   = 3

  INTEGER, PARAMETER :: GRID_UNSTRUCTURED_CELL = 1
  INTEGER, PARAMETER :: GRID_UNSTRUCTURED_VERT = 2
  INTEGER, PARAMETER :: GRID_UNSTRUCTURED_EDGE = 3

  INTEGER, PARAMETER :: GRID_REGULAR_LONLAT    = 45

! Intermediate fix, until GRID_REFERENCE has completely been removed from ICON
  INTEGER, PARAMETER :: GRID_REFERENCE = 9


  ! Parameters for naming all Z-axes created with zaxisCreate.
  ! Moreover, add_var/add_ref makes use of these parameters instead 
  ! of the ZAXIS types defined in cdi.inc. This makes it much easier 
  ! to distinguish vertical axes of the same type (e.g. ZAXIS_HEIGHT), 
  ! which may have different level heights or numbers. 
  !  
  ENUM, BIND(C)
    ENUMERATOR ::                  &
      ZA_SURFACE = 1             , &
      ! Atmosphere
      ZA_CLOUD_BASE              , &
      ZA_CLOUD_TOP               , &
      ZA_ISOTHERM_ZERO           , &
      ZA_REFERENCE               , &
      ZA_REFERENCE_HALF          , &
      ZA_REFERENCE_HALF_HHL      , &
      ZA_HYBRID                  , &
      ZA_HYBRID_HALF             , &
      ZA_HYBRID_HALF_HHL         , &
      ZA_DEPTH_BELOW_LAND        , &
      ZA_DEPTH_BELOW_LAND_P1     , &
      ZA_SNOW                    , &
      ZA_SNOW_HALF               , &
      ZA_PRESSURE                , &
      ZA_HEIGHT                  , &
      ZA_HEIGHT_2M               , &
      ZA_HEIGHT_10M              , &
      ZA_ALTITUDE                , &
      ZA_MEANSEA                 , &
      ZA_ISENTROPIC              , &
      ZA_TOA                     , &
      ZA_PRESSURE_800            , &
      ZA_PRESSURE_400            , &
      ZA_PRESSURE_0              , &
      ZA_DEPTH_RUNOFF_S          , &
      ZA_DEPTH_RUNOFF_G          , &
      ! Lake 
      ZA_LAKE_BOTTOM             , &
      ZA_LAKE_BOTTOM_HALF        , &
      ZA_MIX_LAYER               , &
      ZA_SEDIMENT_BOTTOM_TW_HALF , &
      ! Ocean
      ZA_DEPTH_BELOW_SEA         , &
      ZA_DEPTH_BELOW_SEA_HALF    , &
      ZA_GENERIC_ICE             , &
      ! HAMOCC sediment
      ZA_OCEAN_SEDIMENT          , &
      ! Volcanic Ash parameters  , needed for ICON-ART
      ZA_PRES_FL_SFC_200         , &
      ZA_PRES_FL_200_350         , &
      ZA_PRES_FL_350_550         , &
      ZA_PRES_FL_SFC_100         , &
      ZA_PRES_FL_100_245         , &
      ZA_PRES_FL_245_390         , &
      ZA_PRES_FL_390_530         , &
      ZA_ATMOSPHERE              , &
      !
      ZA_HEIGHT_2M_LAYER         , &
      !
      !
      !
      ! Max for all the constants above, used to dimension arrays with
      ! one entry for each
      ZA_COUNT
  END ENUM

  ! Array mapping our internal ZA_... constants to the respective CDI
  ! zaxis types.
  !
  ! Note that the ordering of this array is crucial, since it
  ! represents the relationship between the ZA_xxxx constants above
  ! and the CDI-internal constants.
  !
  INTEGER, PARAMETER, PUBLIC :: cdi_zaxis_types(ZA_COUNT-1) = &
    [ ZAXIS_SURFACE            , & ! ZA_SURFACE
    & ZAXIS_CLOUD_BASE         , & ! ZA_CLOUD_BASE
    & ZAXIS_CLOUD_TOP          , & ! ZA_CLOUD_TOP
    & ZAXIS_ISOTHERM_ZERO      , & ! ZA_ISOTHERM_ZERO
    & ZAXIS_REFERENCE          , & ! ZA_REFERENCE
    & ZAXIS_REFERENCE          , & ! ZA_REFERENCE_HALF
    & ZAXIS_REFERENCE          , & ! ZA_REFERENCE_HALF_HHL
    & ZAXIS_HYBRID             , & ! ZA_HYBRID
    & ZAXIS_HYBRID_HALF        , & ! ZA_HYBRID_HALF
    & ZAXIS_HYBRID_HALF        , & ! ZA_HYBRID_HALF_HHL
    & ZAXIS_DEPTH_BELOW_LAND   , & ! ZA_DEPTH_BELOW_LAND
    & ZAXIS_DEPTH_BELOW_LAND   , & ! ZA_DEPTH_BELOW_LAND_P1
    & ZAXIS_GENERIC            , & ! ZA_SNOW
    & ZAXIS_GENERIC            , & ! ZA_SNOW_HALF
    & ZAXIS_PRESSURE           , & ! ZA_PRESSURE
    & ZAXIS_HEIGHT             , & ! ZA_HEIGHT
    & ZAXIS_HEIGHT             , & ! ZA_HEIGHT_2M
    & ZAXIS_HEIGHT             , & ! ZA_HEIGHT_10M
    & ZAXIS_ALTITUDE           , & ! ZA_ALTITUDE
    & ZAXIS_MEANSEA            , & ! ZA_MEANSEA
    & ZAXIS_ISENTROPIC         , & ! ZA_ISENTROPIC
    & ZAXIS_TOA                , & ! ZA_TOA
    & ZAXIS_PRESSURE           , & ! ZA_PRESSURE_800
    & ZAXIS_PRESSURE           , & ! ZA_PRESSURE_400
    & ZAXIS_PRESSURE           , & ! ZA_PRESSURE_0
    & ZAXIS_DEPTH_BELOW_LAND   , & ! ZA_DEPTH_RUNOFF_S
    & ZAXIS_DEPTH_BELOW_LAND   , & ! ZA_DEPTH_RUNOFF_G
    & ZAXIS_LAKE_BOTTOM        , & ! ZA_LAKE_BOTTOM
    & ZAXIS_LAKE_BOTTOM        , & ! ZA_LAKE_BOTTOM_HALF
    & ZAXIS_MIX_LAYER          , & ! ZA_MIX_LAYER
    & ZAXIS_SEDIMENT_BOTTOM_TW , & ! ZA_SEDIMENT_BOTTOM_TW_HALF
    & ZAXIS_DEPTH_BELOW_SEA    , & ! ZA_DEPTH_BELOW_SEA
    & ZAXIS_DEPTH_BELOW_SEA    , & ! ZA_DEPTH_BELOW_SEA_HALF
    & ZAXIS_GENERIC            , & ! ZA_GENERIC_ICE
    & ZAXIS_GENERIC            , & ! ZA_OCEAN_SEDIMENT
    & ZAXIS_PRESSURE           , & ! ZA_PRES_FL_SFC_200
    & ZAXIS_PRESSURE           , & ! ZA_PRES_FL_200_350
    & ZAXIS_PRESSURE           , & ! ZA_PRES_FL_350_550
    & ZAXIS_PRESSURE           , & ! ZA_PRES_FL_SFC_100
    & ZAXIS_PRESSURE           , & ! ZA_PRES_FL_100_245
    & ZAXIS_PRESSURE           , & ! ZA_PRES_FL_245_390
    & ZAXIS_PRESSURE           , & ! ZA_PRES_FL_390_530
    & ZAXIS_ATMOSPHERE         , & ! ZA_ATMOSPHERE
    & ZAXIS_HEIGHT               & ! ZA_HEIGHT_2M_LAYER
    & ]

CONTAINS

  PURE FUNCTION is_2d_field(izaxis)
    LOGICAL :: is_2d_field
    INTEGER, INTENT(IN) :: izaxis

    is_2d_field = (izaxis == ZA_surface)                 .OR.  &
      &           (izaxis == ZA_cloud_base)              .OR.  &
      &           (izaxis == ZA_cloud_top)               .OR.  &
      &           (izaxis == ZA_isotherm_zero)           .OR.  &
      &           (izaxis == ZA_height_2m)               .OR.  &
      &           (izaxis == ZA_height_10m)              .OR.  &
      &           (izaxis == ZA_meansea)                 .OR.  &
      &           (izaxis == ZA_TOA)                     .OR.  &
      &           (izaxis == ZA_LAKE_BOTTOM)             .OR.  &
      &           (izaxis == ZA_LAKE_BOTTOM_HALF)        .OR.  &
      &           (izaxis == ZA_MIX_LAYER)               .OR.  &
      &           (izaxis == ZA_SEDIMENT_BOTTOM_TW_HALF) .OR.  &
      &           (izaxis == ZA_PRESSURE_800)            .OR.  &
      &           (izaxis == ZA_PRESSURE_400)            .OR.  &
      &           (izaxis == ZA_PRESSURE_0)              .OR.  &
      &           (izaxis == ZA_DEPTH_RUNOFF_S)          .OR.  &
      &           (izaxis == ZA_DEPTH_RUNOFF_G)          .OR.  &
      &           (izaxis == ZA_HEIGHT_2M_LAYER)

    IF (.NOT.is_2d_field) THEN
       is_2d_field = (izaxis == ZA_PRES_FL_SFC_200)      .OR.  &
                     (izaxis == ZA_PRES_FL_200_350)      .OR.  &
                     (izaxis == ZA_PRES_FL_350_550)      .OR.  &
                     (izaxis == ZA_PRES_FL_SFC_100)      .OR.  &
                     (izaxis == ZA_PRES_FL_100_245)      .OR.  &
                     (izaxis == ZA_PRES_FL_245_390)      .OR.  &
                     (izaxis == ZA_PRES_FL_390_530)      .OR.  &
                     (izaxis == ZA_ATMOSPHERE)
     END IF
  END FUNCTION is_2d_field

END MODULE mo_cdi_constants
