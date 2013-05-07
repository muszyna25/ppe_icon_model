!! Definition of constants for remapping algorithm.
!!
MODULE mo_remap_config

  IMPLICIT NONE
  PUBLIC

  ! Conservative remapping: Maximum size of array-structured
  ! stencil. Every entry that does not fit into this array will be
  ! appended to a sequential list (with some performance losses)
#ifdef __SX__
  INTEGER, PARAMETER :: MAX_NSTENCIL_CONS = 50
#else
  INTEGER, PARAMETER :: MAX_NSTENCIL_CONS = 32
#endif


  ! RBF interpolation: Maximum size of coefficient stencil
  INTEGER, PARAMETER :: MAX_NSTENCIL_RBF = 4

  ! vertex-neighbor cells: stencil size
  INTEGER, PARAMETER :: N_VNB_STENCIL_ICON = 13

  ! level of output verbosity
#ifdef __ICON__
  INTEGER            :: dbg_level =  1
#else
  INTEGER            :: dbg_level =  3
#endif

  ! max. name string length
  INTEGER, PARAMETER :: MAX_NAME_LENGTH = 128

  ! meta data ("config state")
  INTEGER, PARAMETER :: MAX_INPUT_FIELDS = 50 !< maximum number of INPUT fields
  INTEGER, PARAMETER :: MAX_NZAXIS       = 15 !< maximum number of z-axis objects

  ! parallelization
  INTEGER, PARAMETER :: MIN_NFOREIGN = 50000  !< Initial list size for inter-process comm.

  ! RBF constants
  INTEGER, PARAMETER :: rbf_vec_dim   = 4     !< size of the RBF stencil

END MODULE mo_remap_config
