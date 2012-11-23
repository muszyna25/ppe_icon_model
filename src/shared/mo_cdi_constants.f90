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
  INTEGER, PARAMETER, PUBLIC      :: ZA_HYBRID              =  2
  INTEGER, PARAMETER, PUBLIC      :: ZA_HYBRID_HALF         =  3
  INTEGER, PARAMETER, PUBLIC      :: ZA_DEPTH_BELOW_LAND    =  4
  INTEGER, PARAMETER, PUBLIC      :: ZA_DEPTH_BELOW_LAND_P1 =  5
  INTEGER, PARAMETER, PUBLIC      :: ZA_GENERIC_SNOW        =  6
  INTEGER, PARAMETER, PUBLIC      :: ZA_GENERIC_SNOW_P1     =  7
  INTEGER, PARAMETER, PUBLIC      :: ZA_PRESSURE            =  8
  INTEGER, PARAMETER, PUBLIC      :: ZA_HEIGHT              =  9
  INTEGER, PARAMETER, PUBLIC      :: ZA_HEIGHT_2M           = 10
  INTEGER, PARAMETER, PUBLIC      :: ZA_HEIGHT_10M          = 11
  INTEGER, PARAMETER, PUBLIC      :: ZA_ALTITUDE            = 12
  INTEGER, PARAMETER, PUBLIC      :: ZA_MEANSEA             = 13
  INTEGER, PARAMETER, PUBLIC      :: ZA_ISENTROPIC          = 14
  INTEGER, PARAMETER, PUBLIC      :: ZA_TOA                 = 15
  ! Ocean
  INTEGER, PARAMETER, PUBLIC      :: ZA_DEPTH_BELOW_SEA     = 16
  INTEGER, PARAMETER, PUBLIC      :: ZA_DEPTH_BELOW_SEA_HALF= 17
  INTEGER, PARAMETER, PUBLIC      :: ZA_GENERIC_ICE         = 18

CONTAINS

  PURE FUNCTION is_2d_field(izaxis)
    LOGICAL :: is_2d_field
    INTEGER, INTENT(IN) :: izaxis

    is_2d_field = (izaxis == ZA_surface)    .OR.  &
      &           (izaxis == ZA_height_2m)  .OR.  &
      &           (izaxis == ZA_height_10m) .OR.  &
      &           (izaxis == ZA_meansea)
  END FUNCTION is_2d_field

END MODULE mo_cdi_constants
