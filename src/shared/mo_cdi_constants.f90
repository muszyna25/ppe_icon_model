MODULE mo_cdi_constants

  IMPLICIT NONE

  PUBLIC

  INCLUDE 'cdi.inc'

  INTEGER, PARAMETER :: GRID_CELL   = 1
  INTEGER, PARAMETER :: GRID_VERTEX = 2
  INTEGER, PARAMETER :: GRID_EDGE   = 3

  INTEGER, PARAMETER :: GRID_UNSTRUCTURED_CELL = 1
  INTEGER, PARAMETER :: GRID_UNSTRUCTURED_VERT = 2
  INTEGER, PARAMETER :: GRID_UNSTRUCTURED_EDGE = 3

  INTEGER, PARAMETER :: GRID_REGULAR_LONLAT    = 4



  ! Parameters for naming all Z-axes created with zaxisCreate.
  ! Moreover, add_var/add_ref makes use of these parameters instead 
  ! of the ZAXIS types defined in cdi.inc. This makes it much easier 
  ! to distinguish vertical axes of the same type (e.g. ZAXIS_HEIGHT), 
  ! which may have different level heights or numbers. 
  !  
  INTEGER, PARAMETER, PUBLIC      :: ZA_SURFACE             =  1
  ! Atmosphere
  INTEGER, PARAMETER, PUBLIC      :: ZA_CLOUD_BASE          =  2
  INTEGER, PARAMETER, PUBLIC      :: ZA_CLOUD_TOP           =  3
  INTEGER, PARAMETER, PUBLIC      :: ZA_ISOTHERM_ZERO       =  4
  INTEGER, PARAMETER, PUBLIC      :: ZA_REFERENCE           =  5
  INTEGER, PARAMETER, PUBLIC      :: ZA_REFERENCE_HALF      =  6
  INTEGER, PARAMETER, PUBLIC      :: ZA_HYBRID              =  7
  INTEGER, PARAMETER, PUBLIC      :: ZA_HYBRID_HALF         =  8
  INTEGER, PARAMETER, PUBLIC      :: ZA_HYBRID_HALF_HHL     =  9
  INTEGER, PARAMETER, PUBLIC      :: ZA_DEPTH_BELOW_LAND    = 10
  INTEGER, PARAMETER, PUBLIC      :: ZA_DEPTH_BELOW_LAND_P1 = 11
  INTEGER, PARAMETER, PUBLIC      :: ZA_SNOW                = 12
  INTEGER, PARAMETER, PUBLIC      :: ZA_SNOW_HALF           = 13
  INTEGER, PARAMETER, PUBLIC      :: ZA_PRESSURE            = 14
  INTEGER, PARAMETER, PUBLIC      :: ZA_HEIGHT              = 15
  INTEGER, PARAMETER, PUBLIC      :: ZA_HEIGHT_2M           = 16
  INTEGER, PARAMETER, PUBLIC      :: ZA_HEIGHT_10M          = 17
  INTEGER, PARAMETER, PUBLIC      :: ZA_ALTITUDE            = 18
  INTEGER, PARAMETER, PUBLIC      :: ZA_MEANSEA             = 19
  INTEGER, PARAMETER, PUBLIC      :: ZA_ISENTROPIC          = 20
  INTEGER, PARAMETER, PUBLIC      :: ZA_TOA                 = 21
  INTEGER, PARAMETER, PUBLIC      :: ZA_PRESSURE_800        = 22
  INTEGER, PARAMETER, PUBLIC      :: ZA_PRESSURE_400        = 23
  INTEGER, PARAMETER, PUBLIC      :: ZA_PRESSURE_0          = 24
  INTEGER, PARAMETER, PUBLIC      :: ZA_DEPTH_RUNOFF_S      = 25
  INTEGER, PARAMETER, PUBLIC      :: ZA_DEPTH_RUNOFF_G      = 26
  ! Ocean
  INTEGER, PARAMETER, PUBLIC      :: ZA_DEPTH_BELOW_SEA     = 27
  INTEGER, PARAMETER, PUBLIC      :: ZA_DEPTH_BELOW_SEA_HALF= 28
  INTEGER, PARAMETER, PUBLIC      :: ZA_GENERIC_ICE         = 29

CONTAINS

  PURE FUNCTION is_2d_field(izaxis)
    LOGICAL :: is_2d_field
    INTEGER, INTENT(IN) :: izaxis

    is_2d_field = (izaxis == ZA_surface)       .OR.  &
      &           (izaxis == ZA_cloud_base)    .OR.  &
      &           (izaxis == ZA_cloud_top)     .OR.  &
      &           (izaxis == ZA_isotherm_zero) .OR.  &
      &           (izaxis == ZA_height_2m)     .OR.  &
      &           (izaxis == ZA_height_10m)    .OR.  &
      &           (izaxis == ZA_meansea)       .OR.  &
      &           (izaxis == ZA_TOA)           .OR.  &
      &           (izaxis == ZA_PRESSURE_800)  .OR.  &
      &           (izaxis == ZA_PRESSURE_400)  .OR.  &
      &           (izaxis == ZA_PRESSURE_0)
  END FUNCTION is_2d_field

END MODULE mo_cdi_constants
