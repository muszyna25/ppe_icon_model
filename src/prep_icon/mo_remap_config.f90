!! Definition of constants for remapping algorithm.
!!
MODULE mo_remap_config

  IMPLICIT NONE
  PUBLIC

  ! Maximum size of array-structured stencil. Every entry that does
  ! not fit into this array will be appended to a sequential list
  ! (with some performance losses)
  INTEGER, PARAMETER :: MAX_NSTENCIL = 40

  ! vertex-neighbor cells: stencil size
  INTEGER, PARAMETER :: N_VNB_STENCIL = 13

  ! level of output verbosity
  INTEGER, PARAMETER :: dbg_level =  10

  ! max. name string length   
  INTEGER, PARAMETER :: MAX_NAME_LENGTH = 32

END MODULE mo_remap_config
