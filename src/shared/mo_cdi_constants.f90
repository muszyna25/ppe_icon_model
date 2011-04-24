MODULE mo_cdi_constants

  IMPLICIT NONE

  PUBLIC

#ifdef HAVE_CDI
  INCLUDE 'cdi.inc'
#else
  ! just a small extract for testing without CDI - it'll be mandatory to 
  ! have CDI a little bit later ...

  INTEGER, PARAMETER :: GRID_UNSTRUCTURED      =  9
  INTEGER, PARAMETER :: GRID_REFERENCE         = 15
  
  INTEGER, PARAMETER :: ZAXIS_HYBRID           =  2
  INTEGER, PARAMETER :: ZAXIS_HYBRID_HALF      =  3
  
  INTEGER, PARAMETER :: CALENDAR_PROLEPTIC     =  1
  
  INTEGER, PARAMETER :: COMPRESS_SZIP          =  1
  
  INTEGER, PARAMETER :: FILETYPE_GRB2          =  2
  INTEGER, PARAMETER :: FILETYPE_NC2           =  4
  INTEGER, PARAMETER :: FILETYPE_NC4           =  5

#endif

  INTEGER, PARAMETER :: GRID_CELL   = 1
  INTEGER, PARAMETER :: GRID_VERTEX = 2
  INTEGER, PARAMETER :: GRID_EDGE   = 3

END MODULE mo_cdi_constants
