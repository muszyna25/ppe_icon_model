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
  INTEGER, PARAMETER, PUBLIC      :: ZA_surface             =  1
  ! Atmosphere
  INTEGER, PARAMETER, PUBLIC      :: ZA_hybrid              =  2
  INTEGER, PARAMETER, PUBLIC      :: ZA_hybrid_half         =  3
  INTEGER, PARAMETER, PUBLIC      :: ZA_depth_below_land    =  4
  INTEGER, PARAMETER, PUBLIC      :: ZA_depth_below_land_p1 =  5
  INTEGER, PARAMETER, PUBLIC      :: ZA_generic_snow        =  6
  INTEGER, PARAMETER, PUBLIC      :: ZA_generic_snow_p1     =  7
  INTEGER, PARAMETER, PUBLIC      :: ZA_pressure            =  8
  INTEGER, PARAMETER, PUBLIC      :: ZA_height              =  9
  INTEGER, PARAMETER, PUBLIC      :: ZA_height_2m           = 10
  INTEGER, PARAMETER, PUBLIC      :: ZA_height_10m          = 11
  INTEGER, PARAMETER, PUBLIC      :: ZA_altitude            = 12
  INTEGER, PARAMETER, PUBLIC      :: ZA_meansea             = 13
  INTEGER, PARAMETER, PUBLIC      :: ZA_isentropic          = 14
  INTEGER, PARAMETER, PUBLIC      :: ZA_TOA                 = 15
  ! Ocean
  INTEGER, PARAMETER, PUBLIC      :: ZA_depth_below_sea     = 16
  INTEGER, PARAMETER, PUBLIC      :: ZA_depth_below_sea_half= 17
  INTEGER, PARAMETER, PUBLIC      :: ZA_generic_ice         = 18



END MODULE mo_cdi_constants
