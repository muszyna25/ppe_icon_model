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
  USE mo_impl_constants,        ONLY: TLEV_NNOW, VARNAME_LEN, &
    & VINTP_METHOD_LIN, HINTP_TYPE_LONLAT_RBF
  USE mo_grib2,                 ONLY: t_grib2_var
  USE mo_action_types,          ONLY: t_var_action
  USE mo_cf_convention,         ONLY: t_cf_var
  USE mo_tracer_metadata_types, ONLY: t_tracer_meta
  USE mo_model_domain,          ONLY: t_subset_range
  USE mo_var_groups,            ONLY: MAX_GROUPS
  USE ISO_C_BINDING,            ONLY: C_F_POINTER, C_LOC
#ifdef __PGI
  USE ISO_C_BINDING,            ONLY: C_SIZEOF
#endif
  USE mo_cdi,                   ONLY: TSTEP_INSTANT, CDI_UNDEFID

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
  INTEGER, PARAMETER, PUBLIC   :: POST_OP_LIN2DBZ   =  4  !< convert linear values to dbz: dbzval = 10*log10(val)

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
  INTEGER, PARAMETER, PUBLIC :: CLASS_DISTR_STAT    = 8   !< variable based on a distribuition function (PDT 67)
                                                          !< statistical process


  ! ---------------------------------------------------------------
  ! META-DATA TYPE DEFINITIONS
  ! ---------------------------------------------------------------

  TYPE t_union_vals
    REAL(dp) :: rval = 0._dp
    REAL(sp) :: sval = 0._sp
    INTEGER  :: ival = 0
    LOGICAL  :: lval = .FALSE.
  END type t_union_vals

  !> data specific for pz-level interpolation.
  TYPE t_vert_interp_meta
    ! meta data containing the groups to which a variable belongs
    LOGICAL  :: vert_intp_type(SIZE(VINTP_TYPE_LIST)) = .FALSE.
    INTEGER  :: vert_intp_method                      = VINTP_METHOD_LIN
    LOGICAL  :: l_hires_intp                          = .FALSE., &
         &      l_restore_fricred                     = .FALSE., &
         &      l_loglin                              = .FALSE., &
         &      l_extrapol                            = .TRUE.,  &
         &      l_satlimit                            = .FALSE., &
         &      l_restore_pbldev                      = .FALSE., &
         &      l_pd_limit                            = .FALSE., &
         &      l_restore_sfcinv, l_hires_corr
    REAL(wp) :: lower_limit = 0._wp, extrapol_dist
  END TYPE t_vert_interp_meta


  !> data specific for horizontal interpolation.
  TYPE t_hor_interp_meta
    INTEGER :: hor_intp_type = HINTP_TYPE_LONLAT_RBF ! NONE/RBF/Nearest-Neighbor/...
    INTEGER :: fallback_type = HINTP_TYPE_LONLAT_RBF ! replaces "hor_intp_type" if this is not feasible
    INTEGER :: lonlat_id     = 0 ! lon-lat grid (ID in global list)
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
    INTEGER                    :: ipost_op_type = POST_OP_NONE !< type of post-processing operation
    !
    LOGICAL                    :: lnew_cf       = .FALSE.
    TYPE(t_cf_var)             :: new_cf        = t_cf_var('', '', '', -1) !< CF information of modified field
    LOGICAL                    :: lnew_grib2    = .FALSE.
    TYPE(t_grib2_var)          :: new_grib2 !< GRIB2 information of modified field
    !
    TYPE(t_union_vals)         :: arg1 !< post-op argument (e.g. scaling factor)
  END TYPE t_post_op_meta


  TYPE t_var_metadata
    !
    INTEGER                    :: key         = 0              ! hash value of name
    CHARACTER(len=VARNAME_LEN) :: name        = ''             ! variable name
    INTEGER                    :: var_class   = CLASS_DEFAULT  ! variable type
    !                                                   ! 0: CLASS_DEFAULT, 1: CLASS_TILE, ... 
    INTEGER                    :: data_type             ! variable data type: REAL_T, SINGLE_T, INT_T, BOOL_T
    !
    TYPE(t_cf_var)             :: cf          = t_cf_var('', '', '', -1)  ! CF convention information 
    TYPE(t_grib2_var)          :: grib2  ! GRIB2 related information
    !
    LOGICAL                    :: allocated   = .FALSE.        ! allocation status
    INTEGER                    :: ndims       = 0              ! number of dimensions used
    INTEGER                    :: used_dimensions(5) = 0     ! final dimensions of variable
    ! 
    LOGICAL                    :: lrestart              ! write field to restart
    LOGICAL                    :: loutput     = .TRUE.         ! write field to output
    INTEGER                    :: isteptype   = TSTEP_INSTANT  ! Type of statistical processing
    !                                          
    TYPE(t_union_vals)         :: resetval                     ! reset value for accumulated fields
    LOGICAL                    :: lrestart_cont = .FALSE.      ! continue if not in restart file     
    LOGICAL                    :: lrestart_read = .FALSE.      ! field has been set from restart file
    TYPE(t_union_vals)         :: initval                      ! value if not in restart file
    !     
    LOGICAL                    :: lcontainer   = .FALSE.       ! true, if this is a container
    LOGICAL                    :: lcontained   = .FALSE.       ! true, if this is in a container
    INTEGER                    :: ncontained   = 0             ! index in container
    INTEGER                    :: maxcontained = 0             ! container size   
    INTEGER                    :: var_ref_pos  = -1            ! for containers: dimension index for references
    !
    INTEGER                    :: hgrid        = -1            ! CDI horizontal grid type
    INTEGER                    :: vgrid        = -1            ! CDI vertical grid type
    TYPE(t_subset_range)       :: subset                ! subset for latter field access
    INTEGER, POINTER           :: dom                   ! pointer to the variable list
    !
    INTEGER                    :: tlev_source  = TLEV_NNOW     ! Information where to find the actual
    !                                                     timelevel for timelevel dependent variables:        
    !                                                      = 0 : nnow
    !                                                      = 1 : nnow_rcf
    !                                                      ... more may follow
    !
    INTEGER                    :: cdiVarID     = CDI_UNDEFID
    INTEGER                    :: cdiGridID    = CDI_UNDEFID
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
    INTEGER :: l_pp_scheduler_task             = 0

    ! Metadata for missing value masking

    LOGICAL                    :: lmiss          ! flag: true, if variable should be initialized with missval
    TYPE(t_union_vals)         :: missval        ! missing value
    LOGICAL                    :: lmask_boundary ! flag: true, if interpolation zone should be masked *in output*

    ! Index of tracer in tracer and in diagnostics container
    INTEGER                    :: idx_tracer          !< index of tracer in tracer container
    INTEGER                    :: idx_diag            !< index of tracer in diagnostics container
    LOGICAL                    :: lopenacc     = .FALSE.   ! Variable exists on GPU

  CONTAINS
    PROCEDURE :: finalize => t_var_metadata_finalize   !< destructor
  END TYPE t_var_metadata

  TYPE t_var_metadata_ptr
    TYPE(t_var_metadata), POINTER :: p => NULL()
  END TYPE t_var_metadata_ptr
  ! The type t_var_metadata_dynamic is (in contrast to t_var_metadata) not transfered to the output PE.
  ! This allows for dynamical objects inside t_var_metadata_dynamic like pointers or allocatables.
  TYPE t_var_metadata_dynamic
    CLASS(t_tracer_meta), POINTER       :: tracer      ! Tracer-specific metadata

  CONTAINS
    PROCEDURE :: finalize => t_var_metadata_dynamic_finalize   !< destructor
  END TYPE t_var_metadata_dynamic


  PUBLIC :: VINTP_TYPE_LIST
  PUBLIC :: VARNAME_LEN
  PUBLIC :: MAX_GROUPS

  PUBLIC :: t_union_vals
  PUBLIC :: t_var_metadata, t_var_metadata_ptr, var_metadata_get_size
  PUBLIC :: var_metadata_toBinary, var_metadata_fromBinary
  PUBLIC :: t_var_metadata_dynamic
  PUBLIC :: t_vert_interp_meta
  PUBLIC :: t_hor_interp_meta
  PUBLIC :: t_post_op_meta

CONTAINS

  SUBROUTINE t_var_metadata_finalize(this)
    CLASS(t_var_metadata), INTENT(INOUT) :: this
    ! nothing to be done (yet)
  END SUBROUTINE t_var_metadata_finalize

  !-------------------------------------------------------------------------------------------------
  !> @return size of a single variable's info object
  !
  !  @author F. Prill, DWD
  !
  INTEGER FUNCTION var_metadata_get_size()
    TYPE(t_var_metadata) :: info
#ifndef __PGI
    INTEGER, PARAMETER :: info_size = STORAGE_SIZE(info) / STORAGE_SIZE(1)
#else
    INTEGER, SAVE :: info_size = -1

    IF (info_size .EQ. -1) info_size = C_SIZEOF(info) / C_SIZEOF(1)
#endif
    var_metadata_get_size = info_size
  END FUNCTION var_metadata_get_size

  FUNCTION var_metadata_toBinary(info) RESULT(binptr)
    TYPE(t_var_metadata), TARGET, INTENT(IN) :: info
    INTEGER, POINTER :: binptr(:)

    CALL C_F_POINTER(C_LOC(info), binptr, (/ var_metadata_get_size() /))
  END FUNCTION var_metadata_toBinary

  FUNCTION var_metadata_fromBinary(bin) RESULT(infoptr)
    INTEGER, TARGET, INTENT(IN) :: bin(:)
    TYPE(t_var_metadata), POINTER :: infoptr

    NULLIFY(infoptr)
    IF (SIZE(bin) .EQ. var_metadata_get_size()) &
      & CALL C_F_POINTER(C_LOC(bin), infoptr)
  END FUNCTION var_metadata_fromBinary

  SUBROUTINE t_var_metadata_dynamic_finalize(this)
    CLASS(t_var_metadata_dynamic), INTENT(INOUT) :: this
    IF (ASSOCIATED(this%tracer)) THEN
      DEALLOCATE(this%tracer)
    END IF
  END SUBROUTINE t_var_metadata_dynamic_finalize

END MODULE mo_var_metadata_types
