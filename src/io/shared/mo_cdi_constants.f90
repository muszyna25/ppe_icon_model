!> Some CDI-specific constants
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_cdi_constants

  USE mo_zaxis_type, ONLY: ZA_SURFACE, ZA_DEPTH_BELOW_LAND
  PUBLIC

  !------------------------------------------------!
  !  CDI constants for horizontal grid
  !------------------------------------------------!

  INTEGER, PARAMETER :: GRID_CELL   = 1
  INTEGER, PARAMETER :: GRID_VERTEX = 2
  INTEGER, PARAMETER :: GRID_EDGE   = 3

  INTEGER, PARAMETER :: GRID_UNSTRUCTURED_CELL = 1
  INTEGER, PARAMETER :: GRID_UNSTRUCTURED_VERT = 2
  INTEGER, PARAMETER :: GRID_UNSTRUCTURED_EDGE = 3

  INTEGER, PARAMETER :: GRID_REGULAR_LONLAT    = 45

  ! Some constants and variables are re-published
  ! (backward-compatibility needed for JSBACH adapter)
  !
  ! Intermediate fix, until GRID_REFERENCE has completely been removed
  ! from ICON
  INTEGER, PARAMETER :: GRID_REFERENCE = 9

  PUBLIC :: ZA_SURFACE, ZA_DEPTH_BELOW_LAND

END MODULE mo_cdi_constants
