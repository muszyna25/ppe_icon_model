!! Meta-data type definitions for ICON variables.
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
MODULE mo_var_metadata_types

  USE mo_kind,                  ONLY: dp, wp, sp
  USE mo_impl_constants,        ONLY: VARNAME_LEN
  USE mo_grib2,                 ONLY: t_grib2_var
  USE mo_action_types,          ONLY: t_var_action
  USE mo_cf_convention,         ONLY: t_cf_var
  USE mo_tracer_metadata_types, ONLY: t_tracer_meta
  USE mo_model_domain,          ONLY: t_subset_range
  USE mo_var_groups,            ONLY: MAX_GROUPS

  IMPLICIT NONE

  PRIVATE

  !> module name string
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_var_metadata_types'


  ! ---------------------------------------------------------------
  ! CONSTANTS
  ! ---------------------------------------------------------------


  ! list of vertical interpolation types
  ! 
  ! A variable can have any combination of this which means that it
  ! can be interpolated vertically in these different ways.
  CHARACTER(len=1), PARAMETER :: VINTP_TYPE_LIST(3) = &
    (/ "Z",  "P", "I" /)


  ! list of available post-op's (small arithmetic operations on
  ! fields). The implementation is placed in "mo_post_op.f90".o
  INTEGER, PARAMETER, PUBLIC   :: POST_OP_NONE      = -1  !< trivial post-op ("do nothing")
  INTEGER, PARAMETER, PUBLIC   :: POST_OP_SCALE     =  1  !< multiply by scalar factor "arg1"
  INTEGER, PARAMETER, PUBLIC   :: POST_OP_RHO       =  2  !< multiply by rho to get densities instead
  INTEGER, PARAMETER, PUBLIC   :: POST_OP_LUC       =  3  !< convert landuse classes from internal values 
                                                          !< to GRIB2 values (table 4.243) and vice versa. 

  ! list of available variable classes
  INTEGER, PARAMETER, PUBLIC :: CLASS_DEFAULT       = 0
  INTEGER, PARAMETER, PUBLIC :: CLASS_TILE          = 1   !< variable contains tile-specific information
  INTEGER, PARAMETER, PUBLIC :: CLASS_TILE_LAND     = 2   !< variable contains tile-specific information
                                                          !< but is restricted to land-tiles only
  INTEGER, PARAMETER, PUBLIC :: CLASS_SYNSAT        = 3
  INTEGER, PARAMETER, PUBLIC :: CLASS_CHEM          = 4   !< atmospheric chemical constituent (PDT 40)
  INTEGER, PARAMETER, PUBLIC :: CLASS_CHEM_STAT     = 5   !< atmospheric chemical constituent (PDT 42)
                                                          !< statistical process
  INTEGER, PARAMETER, PUBLIC :: CLASS_CHEM_OPTP     = 6   !< atmospheric chemical constituent (PDT 48)
                                                          !< optical properties
  INTEGER, PARAMETER, PUBLIC :: CLASS_DISTR         = 7   !< variable based on a distribuition function (PDT 57)
  INTEGER, PARAMETER, PUBLIC :: CLASS_DISTR_STAT    = 8   !< variable based on a distribuition function (PDT 40467)
                                                          !< statistical process


  ! ---------------------------------------------------------------
  ! META-DATA TYPE DEFINITIONS
  ! ---------------------------------------------------------------

  TYPE t_union_vals
    REAL(dp) :: rval
    REAL(sp) :: sval
    INTEGER  :: ival
    LOGICAL  :: lval
  END type t_union_vals


  !> data specific for pz-level interpolation.
  TYPE t_vert_interp_meta
    ! meta data containing the groups to which a variable belongs
    LOGICAL  :: vert_intp_type(SIZE(VINTP_TYPE_LIST))
    INTEGER  :: vert_intp_method
    LOGICAL  :: l_hires_intp, l_restore_fricred, l_loglin, &
         &      l_extrapol, l_satlimit, l_restore_pbldev,  &
         &      l_pd_limit, l_restore_sfcinv, l_hires_corr
    REAL(wp) :: lower_limit, extrapol_dist
  END TYPE t_vert_interp_meta


  !> data specific for horizontal interpolation.
  TYPE t_hor_interp_meta
    INTEGER :: hor_intp_type ! NONE/RBF/Nearest-Neighbor/...
    INTEGER :: fallback_type ! replaces "hor_intp_type" if this is not feasible
    INTEGER :: lonlat_id     ! lon-lat grid (ID in global list)
  END TYPE t_hor_interp_meta


  !> This type defines small arithmetic operations ("post-ops") as
  !  post-processing tasks.
  !
  !  These post-processing tasks are restricted to point-wise
  !  operations (no halo synchronization) of a single field, like
  !  value scaling.
  !
  !  @note The "post-ops" are performed at output time and DO NOT
  !        MODIFY THE FIELD ITSELF.
  !
  TYPE t_post_op_meta
    INTEGER                    :: ipost_op_type         !< type of post-processing operation
    !
    LOGICAL                    :: lnew_cf
    TYPE(t_cf_var)             :: new_cf                !< CF information of modified field
    LOGICAL                    :: lnew_grib2
    TYPE(t_grib2_var)          :: new_grib2             !< GRIB2 information of modified field
    !
    TYPE(t_union_vals)         :: arg1                  !< post-op argument (e.g. scaling factor)
  END TYPE t_post_op_meta


  TYPE t_var_metadata
    !
    INTEGER                    :: key                   ! hash value of name
    CHARACTER(len=VARNAME_LEN) :: name                  ! variable name
    INTEGER                    :: var_class             ! variable type
    !                                                   ! 0: CLASS_DEFAULT, 1: CLASS_TILE, ... 
    INTEGER                    :: data_type             ! variable data type: REAL_T, SINGLE_T, INT_T, BOOL_T
    !
    TYPE(t_cf_var)             :: cf                    ! CF convention information 
    TYPE(t_grib2_var)          :: grib2                 ! GRIB2 related information
    !
    LOGICAL                    :: allocated             ! allocation status
    INTEGER                    :: ndims                 ! number of dimensions used
    INTEGER                    :: used_dimensions(5)    ! final dimensions of variable
    ! 
    LOGICAL                    :: lrestart              ! write field to restart
    LOGICAL                    :: loutput               ! write field to output
    INTEGER                    :: isteptype             ! Type of statistical processing
    !                                         
    TYPE(t_union_vals)         :: resetval              ! reset value for accumulated fields
    LOGICAL                    :: lrestart_cont         ! continue if not in restart file     
    LOGICAL                    :: lrestart_read         ! field has been set from restart file
    TYPE(t_union_vals)         :: initval               ! value if not in restart file
    !     
    LOGICAL                    :: lcontainer            ! true, if this is a container
    LOGICAL                    :: lcontained            ! true, if this is in a container
    INTEGER                    :: ncontained            ! index in container
    INTEGER                    :: maxcontained          ! container size   
    INTEGER                    :: var_ref_pos           ! for containers: dimension index for references
    !
    INTEGER                    :: hgrid                 ! CDI horizontal grid type
    INTEGER                    :: vgrid                 ! CDI vertical grid type
    TYPE(t_subset_range)       :: subset             ! subset for latter field access
    !
    INTEGER                    :: tlev_source           ! Information where to find the actual
    !                                                     timelevel for timelevel dependent variables:        
    !                                                      = 0 : nnow
    !                                                      = 1 : nnow_rcf
    !                                                      ... more may follow
    !
    INTEGER                    :: cdiVarID
    INTEGER                    :: cdiGridID
    !
    ! Metadata for "post-ops" (small arithmetic operations)
    !
    TYPE(t_post_op_meta)       :: post_op               !<  "post-op" (small arithmetic operations) for this variable
    !
    ! Metadata for "actions" (regularly triggered events)
    !
    TYPE(t_var_action)         :: action_list
    !
    ! Metadata for vertical/horizontal interpolation
    !
    ! Note that setting these parameters to non-default values does
    ! not mean that interpolation is actually performed for this
    ! variables (this is controlled by namelist settings) but only
    ! that this is possible!
    !
    TYPE(t_vert_interp_meta)   :: vert_interp 
    TYPE(t_hor_interp_meta)    :: hor_interp 
    !
    ! meta data containing the groups to which a variable belongs
    LOGICAL :: in_group(MAX_GROUPS)

    ! Flag: defines, if this field is updated by the internal
    ! post-processing scheduler
    INTEGER :: l_pp_scheduler_task

    ! Metadata for missing value masking

    LOGICAL                    :: lmiss          ! flag: true, if variable should be initialized with missval
    TYPE(t_union_vals)         :: missval        ! missing value
    LOGICAL                    :: lmask_boundary ! flag: true, if interpolation zone should be masked *in output*

  END TYPE t_var_metadata

  ! The type t_var_metadata_dynamic is (in contrast to t_var_metadata) not transfered to the output PE.
  ! This allows for dynamical objects inside t_var_metadata_dynamic like pointers or allocatables.
  TYPE t_var_metadata_dynamic
    CLASS(t_tracer_meta), POINTER       :: tracer      ! Tracer-specific metadata
  END TYPE t_var_metadata_dynamic


  PUBLIC :: VINTP_TYPE_LIST
  PUBLIC :: VARNAME_LEN
  PUBLIC :: MAX_GROUPS

  PUBLIC :: t_union_vals
  PUBLIC :: t_var_metadata
  PUBLIC :: t_var_metadata_dynamic
  PUBLIC :: t_vert_interp_meta
  PUBLIC :: t_hor_interp_meta
  PUBLIC :: t_post_op_meta

END MODULE mo_var_metadata_types
