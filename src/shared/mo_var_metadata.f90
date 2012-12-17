MODULE mo_var_metadata

  USE mo_kind,           ONLY: dp, wp
  USE mo_grib1,          ONLY: t_grib1_var
  USE mo_grib2,          ONLY: t_grib2_var
  USE mo_cf_convention,  ONLY: t_cf_var

  IMPLICIT NONE

  PRIVATE

  ! maximum string length for variable names
  INTEGER, PARAMETER :: varname_len = 32

  ! list of variable groups
  CHARACTER(len=varname_len), PARAMETER :: var_groups(16) = &
    (/ "ALL                   ",  &
    &  "ATMO_ML_VARS          ",  &
    &  "ATMO_PL_VARS          ",  &
    &  "ATMO_ZL_VARS          ",  &
    &  "NH_PROG_VARS          ",  &
    &  "ATMO_DERIVED_VARS     ",  &
    &  "RAD_VARS              ",  &
    &  "PRECIP_VARS           ",  &
    &  "CLOUD_DIAG            ",  &
    &  "PBL_VARS              ",  &
    &  "PHYS_TENDENCIES       ",  &
    &  "LAND_VARS             ",  &
    &  "LAND_TILE_VARS        ",  &
    &  "MULTISNOW_VARS        ",  &
    &  "ADDITIONAL_PRECIP_VARS",  &
    &  "SNOW_VARS             "/)

  TYPE t_union_vals
    REAL(dp) :: rval
    INTEGER  :: ival
    LOGICAL  :: lval
  END type t_union_vals


  TYPE t_tracer_meta
    !
    LOGICAL :: lis_tracer         ! this is a tracer field (TRUE/FALSE)
    !  
    INTEGER :: ihadv_tracer       ! method for horizontal transport
    INTEGER :: ivadv_tracer       ! method for vertical transport
    !
    LOGICAL :: lturb_tracer       ! turbulent transport (TRUE/FALSE)
    LOGICAL :: lsed_tracer        ! sedimentation (TRUE/FALSE)
    LOGICAL :: ldep_tracer        ! dry deposition (TRUE/FALSE)  
    LOGICAL :: lconv_tracer       ! convection  (TRUE/FALSE)
    LOGICAL :: lwash_tracer       ! washout (TRUE/FALSE)
    !
    REAL(wp) :: rdiameter_tracer   ! particle diameter in m
    REAL(wp) :: rrho_tracer        ! particle density in kg m^-3
  END TYPE t_tracer_meta


  !> data specific for pz-level interpolation.
  TYPE t_vert_interp_meta
    INTEGER  :: vert_intp_type, vert_intp_method
    LOGICAL  :: l_hires_intp, l_restore_fricred, l_loglin, &
         &      l_extrapol, l_satlimit, l_restore_pbldev,  &
         &      l_pd_limit, l_restore_sfcinv, l_hires_corr
    REAL(wp) :: lower_limit, extrapol_dist
  END TYPE t_vert_interp_meta

  !> data specific for horizontal interpolation.
  TYPE t_hor_interp_meta
    INTEGER :: hor_intp_type ! NONE/LONLAT
    INTEGER :: lonlat_id     ! lon-lat grid (ID in global list)
  END TYPE t_hor_interp_meta


  TYPE t_var_metadata
    !
    INTEGER                    :: key                   ! hash value of name
    CHARACTER(len=varname_len) :: name                  ! variable name  
    !
    TYPE(t_cf_var)             :: cf                    ! CF convention information 
    TYPE(t_grib1_var)          :: grib1                 ! GRIB1 related information
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
    LOGICAL                    :: lmiss                 ! missing value flag
    TYPE(t_union_vals)         :: missval               ! missing value
    LOGICAL                    :: lrestart_cont         ! continue if not in restart file     
    LOGICAL                    :: lrestart_read         ! field has been set from restart file
    TYPE(t_union_vals)         :: initval               ! value if not in restart file
    !     
    LOGICAL                    :: lcontainer            ! true, if this is a container
    LOGICAL                    :: lcontained            ! true, if this is in a container
    INTEGER                    :: ncontained            ! index in container   
    !
    INTEGER                    :: hgrid                 ! CDI horizontal grid type
    INTEGER                    :: vgrid                 ! CDI vertical grid type
    !
    INTEGER                    :: tlev_source           ! Information where to find the actual
    !                                                     timelevel for timelevel dependent variables:        
    !                                                      = 0 : nnow
    !                                                      = 1 : nnow_rcf
    !                                                      ... more may follow
    !
    INTEGER                    :: cdiVarID
    INTEGER                    :: cdiVarID_2            ! for 2nd vector component in LatLon interpolation
    INTEGER                    :: cdiGridID
    INTEGER                    :: cdiZaxisID
    INTEGER                    :: cdiDataType
    INTEGER                    :: cdiTimeID             ! CDI time mode (TIME_VARIABLE/TIME_CONSTANT)
    !
    TYPE(t_tracer_meta)        :: tracer                ! metadata for tracer fields
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
    LOGICAL :: in_group(SIZE(var_groups))

    ! Flag: defines, if this field is updated by the internal
    ! post-processing scheduler
    INTEGER :: l_pp_scheduler_task

  END TYPE t_var_metadata

  PUBLIC :: t_union_vals
  PUBLIC :: t_var_metadata
  PUBLIC :: t_tracer_meta
  PUBLIC :: t_vert_interp_meta
  PUBLIC :: t_hor_interp_meta
  PUBLIC :: varname_len
  PUBLIC :: var_groups

END MODULE mo_var_metadata
