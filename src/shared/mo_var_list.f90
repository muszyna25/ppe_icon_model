! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
MODULE mo_var_list

#include <icon_contiguous_defines.h>

#if defined (__INTEL_COMPILER) || defined (__PGI) || defined (NAGFOR)
#ifdef VARLIST_INITIZIALIZE_WITH_NAN
  USE, INTRINSIC :: ieee_features
  USE, INTRINSIC :: ieee_arithmetic
  USE, INTRINSIC :: ieee_exceptions
#endif
#endif

  USE mo_kind,             ONLY: sp, dp, i8
  USE mo_cf_convention,    ONLY: t_cf_var
  USE mo_grib2,            ONLY: t_grib2_var, grib2_var
  USE mo_var_groups,       ONLY: var_groups_dyn, groups
  USE mo_var_metadata_types,ONLY: t_var_metadata, t_union_vals,     &
    & t_var_metadata_dynamic, t_var_metadata_ptr, &
    &                            t_vert_interp_meta,                &
    &                            t_hor_interp_meta, t_post_op_meta, &
    &                            MAX_GROUPS, VINTP_TYPE_LIST,       &
    &                            CLASS_TILE, CLASS_TILE_LAND
  USE mo_var_metadata,     ONLY: create_vert_interp_metadata,       &
    &                            create_hor_interp_metadata
  USE mo_tracer_metadata,  ONLY: create_tracer_metadata
  USE mo_tracer_metadata_types, ONLY: t_tracer_meta
  USE mo_var_list_element, ONLY: t_var_list_element
  USE mo_linked_list, ONLY: t_var_list, t_list_element, append_list_element
  USE mo_exception,        ONLY: message, finish, message_text
  USE mo_util_texthash,    ONLY: text_hash_c
  USE mo_util_string,      ONLY: toupper, tolower
  USE mo_impl_constants,   ONLY: VARNAME_LEN, TIMELEVEL_SUFFIX,     &
    &                            STR_HINTP_TYPE, MAX_TIME_LEVELS,   &
    &                            REAL_T, SINGLE_T, BOOL_T, INT_T
  USE mo_fortran_tools,    ONLY: init_contiguous_dp, init_contiguous_sp, &
    &                            init_contiguous_i4, init_contiguous_l
  USE mo_action_types,     ONLY: t_var_action
  USE mo_io_config,        ONLY: restart_file_type

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: print_var_list
  PUBLIC :: print_memory_use

  PUBLIC :: default_var_list_settings ! set default settings for a whole list

  PUBLIC :: add_var                   ! create/allocate a new var_list list entry
  PUBLIC :: add_ref                   ! create/reference a new var_list list entry
  PUBLIC :: get_var                   ! obtain reference to existing list entry

  PUBLIC :: get_var_name              ! return plain variable name (without timelevel)
  PUBLIC :: get_var_timelevel         ! return variable timelevel (or "-1")
  PUBLIC :: get_var_tileidx           ! return variable tile index
  PUBLIC :: get_var_list_element_info ! return a copy of the metadata for a var_list element
  PUBLIC :: get_tracer_info_dyn_by_idx! return a copy of the dynamic metadata of a certain tracer
  PUBLIC :: get_timelevel_string      ! return the default string with timelevel encoded
  PUBLIC :: get_varname_with_timelevel! join varname with timelevel string

  PUBLIC :: fget_var_list_element_r1d
  PUBLIC :: fget_var_list_element_r2d
  PUBLIC :: fget_var_list_element_r3d
  PUBLIC :: fget_var_list_element_s1d
  PUBLIC :: fget_var_list_element_s2d
  PUBLIC :: fget_var_list_element_s3d

  PUBLIC :: find_list_element   ! find an element in the list

 INTERFACE add_var  ! create a new list entry
    MODULE PROCEDURE add_var_list_element_5d
    MODULE PROCEDURE add_var_list_element_r4d
    MODULE PROCEDURE add_var_list_element_r3d
    MODULE PROCEDURE add_var_list_element_r2d
    MODULE PROCEDURE add_var_list_element_r1d
    MODULE PROCEDURE add_var_list_element_s4d
    MODULE PROCEDURE add_var_list_element_s3d
    MODULE PROCEDURE add_var_list_element_s2d
    MODULE PROCEDURE add_var_list_element_s1d
    MODULE PROCEDURE add_var_list_element_i4d
    MODULE PROCEDURE add_var_list_element_i3d
    MODULE PROCEDURE add_var_list_element_i2d
    MODULE PROCEDURE add_var_list_element_i1d
    MODULE PROCEDURE add_var_list_element_l4d
    MODULE PROCEDURE add_var_list_element_l3d
    MODULE PROCEDURE add_var_list_element_l2d
    MODULE PROCEDURE add_var_list_element_l1d
  END INTERFACE add_var

  INTERFACE add_ref
    MODULE PROCEDURE add_var_list_reference_r3d
    MODULE PROCEDURE add_var_list_reference_r2d
    MODULE PROCEDURE add_var_list_reference_s3d
    MODULE PROCEDURE add_var_list_reference_s2d
    MODULE PROCEDURE add_var_list_reference_i2d
  END INTERFACE add_ref

  INTERFACE get_var  ! obtain reference to a list entry
    MODULE PROCEDURE get_var_list_element_r5d
    MODULE PROCEDURE get_var_list_element_r4d
    MODULE PROCEDURE get_var_list_element_r3d
    MODULE PROCEDURE get_var_list_element_r2d
    MODULE PROCEDURE get_var_list_element_r1d
    MODULE PROCEDURE get_var_list_element_s5d
    MODULE PROCEDURE get_var_list_element_s4d
    MODULE PROCEDURE get_var_list_element_s3d
    MODULE PROCEDURE get_var_list_element_s2d
    MODULE PROCEDURE get_var_list_element_s1d
    MODULE PROCEDURE get_var_list_element_i5d
    MODULE PROCEDURE get_var_list_element_i4d
    MODULE PROCEDURE get_var_list_element_i3d
    MODULE PROCEDURE get_var_list_element_i2d
    MODULE PROCEDURE get_var_list_element_i1d
    MODULE PROCEDURE get_var_list_element_l5d
    MODULE PROCEDURE get_var_list_element_l4d
    MODULE PROCEDURE get_var_list_element_l3d
    MODULE PROCEDURE get_var_list_element_l2d
    MODULE PROCEDURE get_var_list_element_l1d
  END INTERFACE get_var

  CHARACTER(*), PARAMETER :: modname = "mo_var_list"

CONTAINS
  !------------------------------------------------------------------------------------------------
  !> @return Plain variable name (i.e. without TIMELEVEL_SUFFIX)
  !
  FUNCTION get_var_name(var)
    CHARACTER(LEN=VARNAME_LEN) :: get_var_name
    TYPE(t_var_list_element)   :: var
    INTEGER :: idx

    idx = INDEX(var%info%name,TIMELEVEL_SUFFIX)
    IF (idx .EQ. 0) THEN
      get_var_name = var%info%name
    ELSE
      get_var_name = var%info%name(1:idx-1)
    END IF
  END FUNCTION get_var_name

  !------------------------------------------------------------------------------------------------
  ! construct varname  with timelevel
  !
  CHARACTER(LEN=VARNAME_LEN) FUNCTION get_varname_with_timelevel(varname,timelevel)
    CHARACTER(LEN=VARNAME_LEN), INTENT(IN) :: varname
    INTEGER, INTENT(IN)        :: timelevel

    get_varname_with_timelevel = TRIM(varname)//get_timelevel_string(timelevel)
  END FUNCTION get_varname_with_timelevel

  !------------------------------------------------------------------------------------------------
  ! construct string for timelevel encoding into variable names
  !
  FUNCTION get_timelevel_string(timelevel) RESULT(suffix)
    INTEGER, INTENT(IN) :: timelevel
    CHARACTER(len=4) :: suffix

    WRITE(suffix,'("'//TIMELEVEL_SUFFIX//'",i1)') timelevel
  END FUNCTION get_timelevel_string

  !------------------------------------------------------------------------------------------------
  !> @return time level (extracted from time level suffix) or "-1"
  !
  INTEGER FUNCTION get_var_timelevel(info) RESULT(tl)
    TYPE(t_var_metadata), INTENT(IN) :: info
    CHARACTER(LEN=*), PARAMETER :: routine = modname//':get_var_timelevel'

    tl = INDEX(info%name,TIMELEVEL_SUFFIX)
    IF (tl .EQ. 0) THEN
      tl = -1
    ELSE
      tl = ICHAR(info%name(tl+3:tl+3)) - ICHAR('0')
      IF (tl .LE. 0 .OR. tl .GT. MAX_TIME_LEVELS) &
        & CALL finish(routine, 'Illegal time level in '//TRIM(info%name))
    END IF
  END FUNCTION get_var_timelevel

  ! return logical if a variable name has a timelevel encoded
  LOGICAL FUNCTION has_time_level(varname)
    CHARACTER(*), INTENT(IN) :: varname

    has_time_level = (0 .EQ. INDEX(varname,TIMELEVEL_SUFFIX))
  END FUNCTION

  !------------------------------------------------------------------------------------------------
  !> @return tile index (extracted from tile index suffix "t_") or "-1"
  !
  INTEGER FUNCTION get_var_tileidx(varname) RESULT(tidx)
    CHARACTER(LEN=*), INTENT(IN) :: varname
    CHARACTER(LEN=*), PARAMETER :: routine = modname//':get_var_tileidx'

    tidx = INDEX(varname,'_t_')
    IF (tidx .NE. 0) THEN
      tidx = ICHAR(varname(+3:tidx+3)) - ICHAR('0')
      IF (tidx .LE. 0) &
        CALL finish(routine, 'Illegal time level in '//TRIM(varname))
    END IF
  END FUNCTION get_var_tileidx

  !------------------------------------------------------------------------------------------------
  !
  ! Change parameters of an already existent output var_list
  !
  SUBROUTINE set_var_list (this_list, output_type, restart_type,   &
      & post_suf, rest_suf, init_suf, loutput, lrestart, linitial, &
      & patch_id, vlevel_type, filename, compression_type, model_type)
    TYPE(t_var_list), INTENT(inout)        :: this_list      ! output var_list to change
    INTEGER,          INTENT(in), OPTIONAL :: output_type, restart_type   ! 'GRIB' or 'NetCDF'
    CHARACTER(len=*), INTENT(in), OPTIONAL :: post_suf, rest_suf, init_suf ! suffix of output/restart/initial file
    LOGICAL,          INTENT(in), OPTIONAL :: loutput, lrestart, linitial  ! in standard output/restart/initial file
    INTEGER,          INTENT(in), OPTIONAL :: patch_id       ! patch ID
    INTEGER,          INTENT(in), OPTIONAL :: vlevel_type    ! 1/2/3 for model/pres./height levels
    CHARACTER(len=*), INTENT(in), OPTIONAL :: filename       ! name of output file
    INTEGER,          INTENT(in), OPTIONAL :: compression_type
    CHARACTER(len=*), INTENT(in), OPTIONAL :: model_type     ! output file associated

    this_list%p%restart_type = restart_file_type
    IF (PRESENT(output_type))      this_list%p%output_type      = output_type
    IF (PRESENT(restart_type))     this_list%p%restart_type     = restart_type
    IF (PRESENT(post_suf))         this_list%p%post_suf         = post_suf
    IF (PRESENT(rest_suf))         this_list%p%rest_suf         = rest_suf
    IF (PRESENT(init_suf))         this_list%p%init_suf         = init_suf
    IF (PRESENT(loutput))          this_list%p%loutput          = loutput
    IF (PRESENT(lrestart))         this_list%p%lrestart         = lrestart
    IF (PRESENT(linitial))         this_list%p%linitial         = linitial
    IF (PRESENT(patch_id))         this_list%p%patch_id         = patch_id
    IF (PRESENT(vlevel_type))      this_list%p%vlevel_type      = vlevel_type
    IF (PRESENT(filename))         this_list%p%filename         = filename
    IF (PRESENT(compression_type)) this_list%p%compression_type =  compression_type
    IF (PRESENT(model_type))       this_list%p%model_type       = model_type
  END SUBROUTINE set_var_list
  !------------------------------------------------------------------------------------------------
  !
  ! Set default meta data of output var_list
  !
  SUBROUTINE default_var_list_settings (this_list, filename, loutput, &
    & lrestart, linitial, post_suf, rest_suf, init_suf, output_type,  &
    & restart_type, compression_type, model_type)
    TYPE(t_var_list), INTENT(INOUT)        :: this_list        ! output var_list
    LOGICAL,          INTENT(IN), OPTIONAL :: loutput, lrestart, linitial  ! in standard output/restart/initial file
    CHARACTER(len=*), INTENT(IN), OPTIONAL :: filename         ! name of output file
    INTEGER,          INTENT(IN), OPTIONAL :: output_type, restart_type   ! 'GRIB' or 'NetCDF'
    CHARACTER(len=*), INTENT(IN), OPTIONAL :: post_suf, rest_suf, init_suf ! suffix of output/restart/initial file
    INTEGER,          INTENT(IN), OPTIONAL :: compression_type ! compression type
    CHARACTER(len=*), INTENT(IN), OPTIONAL :: model_type       ! output file associated

    CALL set_var_list(this_list, output_type=output_type, restart_type=restart_type,    &
      & post_suf=post_suf, rest_suf=rest_suf, init_suf=init_suf, loutput=loutput,       &
      & lrestart=lrestart, linitial=linitial, filename=filename, model_type=model_type, &
      & compression_type=compression_type)
  END SUBROUTINE default_var_list_settings
  !------------------------------------------------------------------------------------------------
  !
  ! Get a copy of the metadata concerning a var_list element
  !
  SUBROUTINE get_var_list_element_info (this_list, name, info)
    !
    TYPE(t_var_list),     INTENT(in)  :: this_list    ! list
    CHARACTER(len=*),     INTENT(in)  :: name         ! name of variable
    TYPE(t_var_metadata), INTENT(out) :: info         ! variable meta data
    !
    TYPE(t_list_element), POINTER :: element
    !
    element => find_list_element (this_list, name)
    IF (ASSOCIATED (element)) THEN
      info = element%field%info
    ENDIF
    !
  END SUBROUTINE get_var_list_element_info
  !------------------------------------------------------------------------------------------------
  !
  ! Get a copy of the dynamic metadata concerning a var_list element by index of the element
  !
  SUBROUTINE get_tracer_info_dyn_by_idx (this_list, ncontained, info_dyn)
    !    
    TYPE(t_var_list),             INTENT(in)  :: this_list    ! list
    INTEGER,                      INTENT(in)  :: ncontained   ! index of variable in container
    TYPE(t_var_metadata_dynamic), INTENT(out) :: info_dyn     ! dynamic variable meta data
    !    
    TYPE(t_list_element), POINTER :: element
    !    
    element => find_tracer_by_index (this_list, ncontained)
    IF (ASSOCIATED (element)) THEN 
      info_dyn = element%field%info_dyn
    ENDIF
    !    
  END SUBROUTINE get_tracer_info_dyn_by_idx

  SUBROUTINE default_var_list_metadata(this_info, this_list)
    TYPE(t_var_metadata), INTENT(out) :: this_info
    TYPE(t_var_list), INTENT(in)      :: this_list
    !
    this_info%grib2               = grib2_var(-1, -1, -1, -1, -1, -1)
    this_info%lrestart            = this_list%p%lrestart
    this_info%lmiss               = this_list%p%lmiss
    this_info%lmask_boundary      = this_list%p%lmask_boundary
    this_info%vert_interp         = create_vert_interp_metadata()
    this_info%hor_interp          = create_hor_interp_metadata()
    this_info%in_group(:)         = groups()
  END SUBROUTINE default_var_list_metadata
  !------------------------------------------------------------------------------------------------
  !
  ! Set parameters of list element already created
  ! (private routine within this module)
  !
  ! Set each parameter in data type var_metadata if the respective
  ! optional parameter is present.
  !
  SUBROUTINE set_var_metadata (info,                                           &
         &                     name, hgrid, vgrid, cf, grib2, ldims,           &
         &                     loutput, lcontainer, lrestart, lrestart_cont,   &
         &                     initval, isteptype, resetval, lmiss, missval,   &
         &                     tlev_source, vert_interp,                       &
         &                     hor_interp, in_group, verbose,                  &
         &                     l_pp_scheduler_task, post_op, action_list,      &
         &                     var_class, data_type, idx_tracer, idx_diag,     &   
         &                     lopenacc)
    !
    TYPE(t_var_metadata),    INTENT(inout)        :: info          ! memory info struct.
    CHARACTER(len=*),        INTENT(in), OPTIONAL :: name          ! variable name
    INTEGER,                 INTENT(in), OPTIONAL :: hgrid         ! horizontal grid type used
    INTEGER,                 INTENT(in), OPTIONAL :: vgrid         ! vertical grid type used
    TYPE(t_cf_var),          INTENT(in), OPTIONAL :: cf            ! CF convention
    TYPE(t_grib2_var),       INTENT(in), OPTIONAL :: grib2         ! GRIB2
    INTEGER,                 INTENT(in)           :: ldims(:)      ! used dimensions
    LOGICAL,                 INTENT(in), OPTIONAL :: loutput       ! into output var_list
    LOGICAL,                 INTENT(in), OPTIONAL :: lcontainer    ! true if container
    LOGICAL,                 INTENT(in), OPTIONAL :: lrestart      ! restart file flag
    LOGICAL,                 INTENT(in), OPTIONAL :: lrestart_cont ! continue on restart
    TYPE(t_union_vals),      INTENT(in), OPTIONAL :: initval       ! value if var not available
    INTEGER,                 INTENT(in), OPTIONAL :: isteptype     ! type of statistical processing
    TYPE(t_union_vals),      INTENT(in), OPTIONAL :: resetval      ! reset value
    LOGICAL,                 INTENT(in), OPTIONAL :: lmiss         ! missing value flag
    TYPE(t_union_vals),      INTENT(in), OPTIONAL :: missval       ! missing value
    INTEGER,                 INTENT(in), OPTIONAL :: tlev_source   ! actual TL for TL dependent vars
    TYPE(t_vert_interp_meta),INTENT(in), OPTIONAL :: vert_interp   ! vertical interpolation metadata
    TYPE(t_hor_interp_meta), INTENT(in), OPTIONAL :: hor_interp    ! horizontal interpolation metadata
    LOGICAL, INTENT(in), OPTIONAL :: in_group(:)          ! groups to which a variable belongs
    LOGICAL,                 INTENT(in), OPTIONAL :: verbose
    INTEGER,                 INTENT(in), OPTIONAL :: l_pp_scheduler_task ! .TRUE., if field is updated by pp scheduler
    TYPE(t_post_op_meta),    INTENT(in), OPTIONAL :: post_op       !< "post-op" (small arithmetic operations) for this variable
    TYPE(t_var_action),      INTENT(in), OPTIONAL :: action_list   !< regularly triggered events
    INTEGER,                 INTENT(in), OPTIONAL :: var_class     ! variable class/species
    INTEGER,                 INTENT(IN), OPTIONAL :: data_type     ! variable data type
    INTEGER,                 INTENT(IN), OPTIONAL :: idx_tracer    ! index of tracer in tracer container 
    INTEGER,                 INTENT(IN), OPTIONAL :: idx_diag      ! index of tracer in diagnostics container 
    LOGICAL,                 INTENT(IN), OPTIONAL :: lopenacc      ! variable data type
    !
    LOGICAL :: lverbose
    !
    ! set flags from optional parameters
    lverbose = .FALSE.
    IF (PRESENT(name)) THEN
      info%name      = name
      info%key = text_hash_c(TRIM(name))
    END IF
    IF (PRESENT(verbose))       lverbose             = verbose
    IF (PRESENT(data_type))     info%data_type       = data_type
    ! set components describing the 'Content of the field'
    IF (PRESENT(var_class))     info%var_class       = var_class
    IF (PRESENT(cf))            info%cf              = cf
    IF (PRESENT(grib2))         info%grib2           = grib2
    IF (PRESENT(hgrid))         info%hgrid           = hgrid
    IF (PRESENT(vgrid))         info%vgrid           = vgrid
    info%used_dimensions = 1
    info%used_dimensions(1:SIZE(ldims)) = ldims
    IF (PRESENT(loutput))       info%loutput         = loutput
    IF (PRESENT(lcontainer))    info%lcontainer      = lcontainer
    IF (info%lcontainer) THEN
      info%ncontained   =  0
      info%var_ref_pos  = -1 ! UNDEFINED
    END IF
    IF (PRESENT(resetval))      info%resetval      = resetval
    IF (PRESENT(isteptype))     info%isteptype     = isteptype
    IF (PRESENT(lmiss))         info%lmiss         = lmiss
    IF (PRESENT(missval))       info%missval       = missval
    IF (PRESENT(lrestart))      info%lrestart      = lrestart
    IF (PRESENT(lrestart_cont)) info%lrestart_cont = lrestart_cont
    IF (PRESENT(initval))       info%initval       = initval
    IF (PRESENT(tlev_source))   info%tlev_source   = tlev_source
    ! set flags concerning vertical interpolation
    IF (PRESENT(vert_interp))   info%vert_interp   = vert_interp
    IF (PRESENT(hor_interp))    info%hor_interp    = hor_interp
    IF (PRESENT(in_group)) &
      & info%in_group(:SIZE(in_group)) = in_group(:)
    IF (PRESENT(l_pp_scheduler_task)) &
      & info%l_pp_scheduler_task = l_pp_scheduler_task
    IF (PRESENT(post_op))       info%post_op       = post_op
    IF (PRESENT(action_list))   info%action_list   = action_list
    ! indices of tracer in tracer container and in diagnostic container
    IF (PRESENT(idx_tracer))    info%idx_tracer    = idx_tracer
    IF (PRESENT(idx_diag))      info%idx_diag      = idx_diag
    IF (PRESENT(lopenacc))      info%lopenacc      = lopenacc
    ! perform consistency checks on variable's meta-data:
    CALL check_metadata_consistency(info)
    !LK    IF (lverbose) CALL print_var_metadata (info)
  END SUBROUTINE set_var_metadata


  !------------------------------------------------------------------------------------------------
  !
  ! Set dynamic metadata, i.e. polymorphic tracer metadata
  ! (private routine within this module)
  !
  SUBROUTINE set_var_metadata_dyn(this_info_dyn,tracer_info)
    TYPE(t_var_metadata_dynamic),INTENT(INOUT) :: this_info_dyn
    CLASS(t_tracer_meta), INTENT(IN), OPTIONAL :: tracer_info
    CLASS(t_tracer_meta), POINTER :: tmp
    
    IF (PRESENT(tracer_info)) THEN
      ALLOCATE(this_info_dyn%tracer, source=tracer_info)
    ELSE
      ALLOCATE(t_tracer_meta :: this_info_dyn%tracer)
      tmp => this_info_dyn%tracer
      SELECT TYPE(tmp)
      TYPE IS(t_tracer_meta)
        tmp = create_tracer_metadata(lis_tracer=.FALSE.)
      END SELECT
    ENDIF
  END SUBROUTINE set_var_metadata_dyn

  !------------------------------------------------------------------------------------------------
  !
  ! Create a list new entry
  !
  ! Specific routines for pointers of different rank
  !
  ! create (allocate) a new table entry
  ! optionally obtain pointer to 5d-field
  ! optionally overwrite default meta data
  !
  SUBROUTINE add_var_list_element_5d(ndims, data_type, this_list, name,         &
    &   hgrid, vgrid, cf, grib2, ldims, new_list_element, loutput, lcontainer,  &
    &   lrestart, lrestart_cont, isteptype, lmiss, tlev_source,                 &
    &   info, vert_interp, hor_interp, in_group, verbose,                       &
    &   l_pp_scheduler_task, post_op, action_list, tracer_info,                 &
    &   p5_r, p5_s, p5_i, p5_l, initval_r, initval_s, initval_i, initval_l,     &
    &   resetval_r, resetval_s, resetval_i, resetval_l,                         &
    &   missval_r, missval_s, missval_i, missval_l, var_class, lopenacc )

    INTEGER,                 INTENT(IN)           :: ndims                        ! used dimensions (1...5)
    INTEGER,                 INTENT(IN)           :: data_type
    TYPE(t_var_list),        INTENT(inout)        :: this_list                    ! list
    CHARACTER(len=*),        INTENT(in)           :: name                         ! name of variable
    INTEGER,                 INTENT(in)           :: hgrid                        ! horizontal grid type used
    INTEGER,                 INTENT(in)           :: vgrid                        ! vertical grid type used
    TYPE(t_cf_var),          INTENT(in)           :: cf                           ! CF related metadata
    TYPE(t_grib2_var),       INTENT(in)           :: grib2                        ! GRIB2 related metadata
    INTEGER,                 INTENT(in)           :: ldims(5)                     ! local dimensions
    TYPE(t_list_element),    POINTER              :: new_list_element             ! pointer to new var list element
    LOGICAL,                 INTENT(in), OPTIONAL :: loutput                      ! output flag
    LOGICAL,                 INTENT(in), OPTIONAL :: lcontainer                   ! container flag
    LOGICAL,                 INTENT(in), OPTIONAL :: lrestart                     ! restart flag
    LOGICAL,                 INTENT(in), OPTIONAL :: lrestart_cont                ! continue restart if var not available
    INTEGER,                 INTENT(in), OPTIONAL :: isteptype                    ! type of statistical processing
    LOGICAL,                 INTENT(in), OPTIONAL :: lmiss                        ! missing value flag
    INTEGER,                 INTENT(in), OPTIONAL :: tlev_source                  ! actual TL for TL dependent vars
    TYPE(t_var_metadata),    POINTER,    OPTIONAL :: info                         ! returns reference to metadata
    REAL(dp),         CONTIGUOUS_TARGET, OPTIONAL :: p5_r(:,:,:,:,:)              ! provided pointer
    REAL(sp),         CONTIGUOUS_TARGET, OPTIONAL :: p5_s(:,:,:,:,:)              ! provided pointer
    INTEGER,          CONTIGUOUS_TARGET, OPTIONAL :: p5_i(:,:,:,:,:)              ! provided pointer
    LOGICAL,          CONTIGUOUS_TARGET, OPTIONAL :: p5_l(:,:,:,:,:)              ! provided pointer
    TYPE(t_vert_interp_meta),INTENT(in), OPTIONAL :: vert_interp                  ! vertical interpolation metadata
    TYPE(t_hor_interp_meta), INTENT(in), OPTIONAL :: hor_interp                   ! horizontal interpolation metadata
    LOGICAL,                 INTENT(in), OPTIONAL :: in_group(:)                  ! groups to which a variable belongs
    LOGICAL,                 INTENT(in), OPTIONAL :: verbose                      ! print information
    INTEGER,                 INTENT(in), OPTIONAL :: l_pp_scheduler_task          ! .TRUE., if field is updated by pp scheduler
    TYPE(t_post_op_meta),    INTENT(IN), OPTIONAL :: post_op                      ! "post-op" (small arithmetic operations) for this variable
    TYPE(t_var_action),      INTENT(IN), OPTIONAL :: action_list                  ! regularly triggered events
    CLASS(t_tracer_meta),    INTENT(in), OPTIONAL :: tracer_info                  ! tracer meta data
    REAL(dp),                INTENT(in), OPTIONAL :: initval_r                    ! value if var not available
    REAL(sp),                INTENT(in), OPTIONAL :: initval_s                    ! value if var not available
    INTEGER,                 INTENT(in), OPTIONAL :: initval_i                    ! value if var not available
    LOGICAL,                 INTENT(in), OPTIONAL :: initval_l                    ! value if var not available
    REAL(dp),                INTENT(in), OPTIONAL :: resetval_r                   ! reset value (after accumulation)
    REAL(sp),                INTENT(in), OPTIONAL :: resetval_s                   ! reset value (after accumulation)
    INTEGER,                 INTENT(in), OPTIONAL :: resetval_i                   ! reset value (after accumulation)
    LOGICAL,                 INTENT(in), OPTIONAL :: resetval_l                   ! reset value (after accumulation)
    REAL(dp),                INTENT(in), OPTIONAL :: missval_r                    ! missing value
    REAL(sp),                INTENT(in), OPTIONAL :: missval_s                    ! missing value
    INTEGER,                 INTENT(in), OPTIONAL :: missval_i                    ! missing value
    LOGICAL,                 INTENT(in), OPTIONAL :: missval_l                    ! missing value
    INTEGER,                 INTENT(in), OPTIONAL :: var_class                    !< variable type/species
    LOGICAL,                 INTENT(in), OPTIONAL :: lopenacc                     ! create variable on GPU

    ! local variables
    TYPE(t_union_vals) :: missval, initval, resetval, ivals
    INTEGER :: idims(5), istat
    LOGICAL :: referenced, is_restart_var
    CHARACTER(LEN = *), PARAMETER :: routine = modname//":add_var_list_element_5d"

    ! consistency check for restart and output
    TYPE(t_list_element), POINTER :: duplicate

    ! Check for a variable of the same name. This consistency check
    ! only makes sense for single-domain setups which, in addition,
    ! must not use internal post-processing (lon-lat or vertically
    ! interpolated output).
!TODO: find way without cycle-dependency
!    IF (msg_level > 20) THEN
!      duplicate => find_element(name)
!      IF (ASSOCIATED(duplicate)) THEN
!        CALL message('ADD_VAR:','Found double entry for varname:'//TRIM(name))
!        NULLIFY(duplicate)
!      ENDIF
!    END IF

    is_restart_var = this_list%p%lrestart
    IF (PRESENT(lrestart)) THEN
      is_restart_var = lrestart
      IF (.NOT. this_list%p%lrestart .AND. lrestart) THEN
        CALL finish(routine, 'for list '//TRIM(this_list%p%name)//' restarting not enabled, '// &
                           & 'but restart of '//TRIM(name)//' requested.')
      ENDIF
    ENDIF
    IF (is_restart_var .AND. (.NOT. ANY(data_type == (/REAL_T, SINGLE_T, INT_T/)))) THEN
      CALL finish(routine, 'unsupported data_type for "'//TRIM(NAME)//'": '// &
        & 'data_type of restart variables must be floating-point or integer type.')
    END IF

    ! add list entry

    CALL append_list_element (this_list, new_list_element)
    CALL default_var_list_metadata(new_list_element%field%info, this_list)

    ! init local fields

    missval = new_list_element%field%info%missval
    initval = new_list_element%field%info%initval
    resetval= new_list_element%field%info%resetval

    ! and set meta data
    referenced = ANY([PRESENT(p5_r), PRESENT(p5_s), PRESENT(p5_i), PRESENT(p5_l)])
    IF (PRESENT(missval_r))  missval%rval  = missval_r
    IF (PRESENT(missval_s))  missval%sval  = missval_s
    IF (PRESENT(missval_i))  missval%ival  = missval_i
    IF (PRESENT(missval_l))  missval%lval  = missval_l
    IF (PRESENT(initval_r))  initval%rval  = initval_r
    IF (PRESENT(initval_s))  initval%sval  = initval_s
    IF (PRESENT(initval_i))  initval%ival  = initval_i
    IF (PRESENT(initval_l))  initval%lval  = initval_l
    IF (PRESENT(resetval_r)) resetval%rval = resetval_r
    IF (PRESENT(resetval_s)) resetval%sval = resetval_s
    IF (PRESENT(resetval_i)) resetval%ival = resetval_i
    IF (PRESENT(resetval_l)) resetval%lval = resetval_l
    CALL set_var_metadata( new_list_element%field%info,                      &
         name=name, hgrid=hgrid, vgrid=vgrid, cf=cf, grib2=grib2,            &
         ldims=ldims(1:ndims), loutput=loutput, lcontainer=lcontainer,       &
         lrestart=lrestart, lrestart_cont=lrestart_cont, initval=initval,    &
         isteptype=isteptype, resetval=resetval, lmiss=lmiss,                &
         missval=missval, tlev_source=tlev_source,                           &
         vert_interp=vert_interp, hor_interp=hor_interp, in_group=in_group,  &
         verbose=verbose, l_pp_scheduler_task=l_pp_scheduler_task,           &
         post_op=post_op, action_list=action_list, var_class=var_class,      &
         data_type=data_type, lopenacc=lopenacc )
    ! set dynamic metadata, i.e. polymorphic tracer metadata
    CALL set_var_metadata_dyn (new_list_element%field%info_dyn,              &
                               tracer_info=tracer_info)

    new_list_element%field%info%ndims                    = ndims
    new_list_element%field%info%used_dimensions(1:ndims) = ldims(1:ndims)
    new_list_element%field%info%dom                      => this_list%p%patch_id
    NULLIFY(new_list_element%field%r_ptr, new_list_element%field%s_ptr, &
      &     new_list_element%field%i_ptr, new_list_element%field%l_ptr)
    IF (.NOT. referenced) THEN
      idims(1:ndims)    = new_list_element%field%info%used_dimensions(1:ndims)
      idims((ndims+1):) = 1
      SELECT CASE(data_type)
      CASE (REAL_T)
        new_list_element%field%var_base_size    = 8
        ALLOCATE(new_list_element%field%r_ptr(idims(1), idims(2), idims(3), idims(4), idims(5)), STAT=istat)
        !$ACC ENTER DATA CREATE( new_list_element%field%r_ptr ) IF( new_list_element%field%info%lopenacc )
      CASE (SINGLE_T)
        new_list_element%field%var_base_size    = 4
        ALLOCATE(new_list_element%field%s_ptr(idims(1), idims(2), idims(3), idims(4), idims(5)), STAT=istat)
        !$ACC ENTER DATA CREATE( new_list_element%field%s_ptr ) IF( new_list_element%field%info%lopenacc )
      CASE (INT_T)
        new_list_element%field%var_base_size    = 4
        ALLOCATE(new_list_element%field%i_ptr(idims(1), idims(2), idims(3), idims(4), idims(5)), STAT=istat)
        !$ACC ENTER DATA CREATE( new_list_element%field%i_ptr ) IF( new_list_element%field%info%lopenacc )
      CASE (BOOL_T)
        new_list_element%field%var_base_size    = 4
        ALLOCATE(new_list_element%field%l_ptr(idims(1), idims(2), idims(3), idims(4), idims(5)), STAT=istat)
        !$ACC ENTER DATA CREATE( new_list_element%field%l_ptr ) IF( new_list_element%field%info%lopenacc )
      END SELECT
      IF (istat /= 0) &
        CALL finish(routine, 'allocation of array '//TRIM(name)//' failed')
      this_list%p%memory_used = this_list%p%memory_used &
           + INT(new_list_element%field%var_base_size, i8) &
           & * INT(PRODUCT(idims(1:5)),i8)
    ELSE
      SELECT CASE(data_type)
      CASE (REAL_T)
        new_list_element%field%r_ptr => p5_r
      CASE (SINGLE_T)
        new_list_element%field%s_ptr => p5_s
      CASE (INT_T)
        new_list_element%field%i_ptr => p5_i
      CASE (BOOL_T)
        new_list_element%field%l_ptr => p5_l
      END SELECT
    ENDIF
    new_list_element%field%info%allocated = .TRUE.
    IF(PRESENT(info)) info => new_list_element%field%info

    ! initialize the new array
#if    defined (VARLIST_INITIZIALIZE_WITH_NAN) \
    && (defined (__INTEL_COMPILER) || defined (__PGI) || defined (NAGFOR))
    ivals%rval = ieee_value(ptr, ieee_signaling_nan)
    ivals%sval = ieee_value(ptr, ieee_signaling_nan)
#endif
    IF (ANY([PRESENT(initval_r), PRESENT(initval_s), PRESENT(initval_i), PRESENT(initval_l)])) THEN
      ivals = initval
    ELSE IF (ANY([PRESENT(missval_r), PRESENT(missval_s), PRESENT(missval_i), PRESENT(missval_l)])) THEN
      ivals = missval
    END IF
    SELECT CASE(data_type)
    CASE (REAL_T)
      !$OMP PARALLEL
      CALL init_contiguous_dp(new_list_element%field%r_ptr, SIZE(new_list_element%field%r_ptr), ivals%rval)
      !$OMP END PARALLEL
      !$ACC UPDATE DEVICE( new_list_element%field%r_ptr ) IF( new_list_element%field%info%lopenacc )
    CASE (SINGLE_T)
      !$OMP PARALLEL
      CALL init_contiguous_sp(new_list_element%field%s_ptr, SIZE(new_list_element%field%s_ptr), ivals%sval)
      !$OMP END PARALLEL
      !$ACC UPDATE DEVICE( new_list_element%field%s_ptr ) IF( new_list_element%field%info%lopenacc )
    CASE (INT_T)
      !$OMP PARALLEL
      CALL init_contiguous_i4(new_list_element%field%i_ptr, SIZE(new_list_element%field%i_ptr), ivals%ival)
      !$OMP END PARALLEL
      !$ACC UPDATE DEVICE( new_list_element%field%i_ptr ) IF( new_list_element%field%info%lopenacc )
    CASE (BOOL_T)
      !$OMP PARALLEL
      CALL init_contiguous_l(new_list_element%field%l_ptr, SIZE(new_list_element%field%l_ptr), ivals%lval)
      !$OMP END PARALLEL
      !$ACC UPDATE DEVICE( new_list_element%field%l_ptr ) IF( new_list_element%field%info%lopenacc )
    END SELECT
  END SUBROUTINE add_var_list_element_5d

  !------------------------------------------------------------------------------------------------
  ! create (allocate) a new table entry
  ! optionally obtain pointer to 4d-field
  ! optionally overwrite default meta data
  !
  SUBROUTINE add_var_list_element_r4d(this_list, name, ptr,       &
    &   hgrid, vgrid, cf, grib2, ldims, loutput, lcontainer,      &
    &   lrestart, lrestart_cont, initval, isteptype,              &
    &   resetval, lmiss, missval, tlev_source, info, p5,          &
    &   vert_interp, hor_interp, in_group, verbose, new_element,  &
    &   l_pp_scheduler_task, post_op, action_list, var_class,     &
    &   lopenacc)
    !
    TYPE(t_var_list),        INTENT(inout)        :: this_list                    ! list
    CHARACTER(len=*),        INTENT(in)           :: name                         ! name of variable
    REAL(dp),       POINTER, INTENT(OUT)          :: ptr(:,:,:,:)                 ! reference to field
    INTEGER,                 INTENT(in)           :: hgrid                        ! horizontal grid type used
    INTEGER,                 INTENT(in)           :: vgrid                        ! vertical grid type used
    TYPE(t_cf_var),          INTENT(in)           :: cf                           ! CF related metadata
    TYPE(t_grib2_var),       INTENT(in)           :: grib2                        ! GRIB2 related metadata
    INTEGER,                 INTENT(in)           :: ldims(:)                     ! local dimensions
    LOGICAL,                 INTENT(in), OPTIONAL :: loutput                      ! output flag
    LOGICAL,                 INTENT(in), OPTIONAL :: lcontainer                   ! container flag
    LOGICAL,                 INTENT(in), OPTIONAL :: lrestart                     ! restart flag
    LOGICAL,                 INTENT(in), OPTIONAL :: lrestart_cont                ! continue restart if var not available
    REAL(dp),                INTENT(in), OPTIONAL :: initval                      ! value if var not available
    INTEGER,                 INTENT(in), OPTIONAL :: isteptype                    ! type of statistical processing
    REAL(dp),                INTENT(in), OPTIONAL :: resetval                     ! reset value (after accumulation)
    LOGICAL,                 INTENT(in), OPTIONAL :: lmiss                        ! missing value flag
    REAL(dp),                INTENT(in), OPTIONAL :: missval                      ! missing value
    INTEGER,                 INTENT(in), OPTIONAL :: tlev_source                  ! actual TL for TL dependent vars
    TYPE(t_var_metadata),    POINTER,    OPTIONAL :: info                         ! returns reference to metadata
    REAL(dp),         CONTIGUOUS_TARGET, OPTIONAL :: p5(:,:,:,:,:)                ! provided pointer
    TYPE(t_vert_interp_meta),INTENT(in), OPTIONAL :: vert_interp                  ! vertical interpolation metadata
    TYPE(t_hor_interp_meta), INTENT(in), OPTIONAL :: hor_interp                   ! horizontal interpolation metadata
    LOGICAL,                 INTENT(in), OPTIONAL :: in_group(:)                  ! groups to which a variable belongs
    LOGICAL,                 INTENT(in), OPTIONAL :: verbose                      ! print information
    TYPE(t_list_element),    POINTER,    OPTIONAL :: new_element                  ! pointer to new var list element
    INTEGER,                 INTENT(in), OPTIONAL :: l_pp_scheduler_task          ! .TRUE., if field is updated by pp scheduler
    TYPE(t_post_op_meta),    INTENT(IN), OPTIONAL :: post_op                      ! "post-op" (small arithmetic operations) for this variable
    TYPE(t_var_action),      INTENT(IN), OPTIONAL :: action_list                  ! regularly triggered events
    INTEGER,                 INTENT(in), OPTIONAL :: var_class                    !< variable type/species
    LOGICAL,                 INTENT(in), OPTIONAL :: lopenacc                     ! create variable on GPU
    ! local variables
    TYPE(t_list_element), POINTER :: element
    INTEGER                       :: idims(5), ndims

    ndims = 4
    idims(1:ndims) = ldims(1:ndims)
    CALL add_var_list_element_5d(ndims, REAL_T, this_list, name,           &
      &   hgrid, vgrid, cf, grib2, idims, element, loutput, lcontainer,    &
      &   lrestart, lrestart_cont, isteptype, lmiss, tlev_source, info,    &
      &   vert_interp, hor_interp, in_group, verbose,                      &
      &   l_pp_scheduler_task, post_op, action_list, p5_r=p5,              &
      &   initval_r=initval, resetval_r=resetval, missval_r=missval,       &
      &   var_class=var_class, lopenacc=lopenacc)
    ptr => element%field%r_ptr(:,:,:,:,1)
    IF (PRESENT(new_element))  new_element => element
  END SUBROUTINE add_var_list_element_r4d


  !------------------------------------------------------------------------------------------------
  ! create (allocate) a new table entry
  ! optionally obtain pointer to 3d-field
  ! optionally overwrite default meta data
  !
  SUBROUTINE add_var_list_element_r3d(this_list, name, ptr,            &
    &   hgrid, vgrid, cf, grib2, ldims, loutput, lcontainer,           &
    &   lrestart, lrestart_cont, initval, isteptype,                   &
    &   resetval, lmiss, missval, tlev_source, tracer_info, info, p5,  &
    &   vert_interp, hor_interp, in_group, verbose, new_element,       &
    &   l_pp_scheduler_task, post_op, action_list, var_class, lopenacc)
    !
    TYPE(t_var_list),        INTENT(inout)        :: this_list                    ! list
    CHARACTER(len=*),        INTENT(in)           :: name                         ! name of variable
    REAL(dp),       POINTER, INTENT(OUT)          :: ptr(:,:,:)                   ! reference to field
    INTEGER,                 INTENT(in)           :: hgrid                        ! horizontal grid type used
    INTEGER,                 INTENT(in)           :: vgrid                        ! vertical grid type used
    TYPE(t_cf_var),          INTENT(in)           :: cf                           ! CF related metadata
    TYPE(t_grib2_var),       INTENT(in)           :: grib2                        ! GRIB2 related metadata
    INTEGER,                 INTENT(in)           :: ldims(:)                     ! local dimensions
    LOGICAL,                 INTENT(in), OPTIONAL :: loutput                      ! output flag
    LOGICAL,                 INTENT(in), OPTIONAL :: lcontainer                   ! container flag
    LOGICAL,                 INTENT(in), OPTIONAL :: lrestart                     ! restart flag
    LOGICAL,                 INTENT(in), OPTIONAL :: lrestart_cont                ! continue restart if var not available
    REAL(dp),                INTENT(in), OPTIONAL :: initval                      ! value if var not available
    INTEGER,                 INTENT(in), OPTIONAL :: isteptype                    ! type of statistical processing
    REAL(dp),                INTENT(in), OPTIONAL :: resetval                     ! reset value (after accumulation)
    LOGICAL,                 INTENT(in), OPTIONAL :: lmiss                        ! missing value flag
    REAL(dp),                INTENT(in), OPTIONAL :: missval                      ! missing value
    INTEGER,                 INTENT(in), OPTIONAL :: tlev_source                  ! actual TL for TL dependent vars
    CLASS(t_tracer_meta),    INTENT(in), OPTIONAL :: tracer_info                  ! tracer meta data
    TYPE(t_var_metadata),    POINTER,    OPTIONAL :: info                         ! returns reference to metadata
    REAL(dp),         CONTIGUOUS_TARGET, OPTIONAL :: p5(:,:,:,:,:)                ! provided pointer
    TYPE(t_vert_interp_meta),INTENT(in), OPTIONAL :: vert_interp                  ! vertical interpolation metadata
    TYPE(t_hor_interp_meta), INTENT(in), OPTIONAL :: hor_interp                   ! horizontal interpolation metadata
    LOGICAL,                 INTENT(in), OPTIONAL :: in_group(:)                  ! groups to which a variable belongs
    LOGICAL,                 INTENT(in), OPTIONAL :: verbose                      ! print information
    TYPE(t_list_element),    POINTER,    OPTIONAL :: new_element                  ! pointer to new var list element
    INTEGER,                 INTENT(in), OPTIONAL :: l_pp_scheduler_task          ! .TRUE., if field is updated by pp scheduler
    TYPE(t_post_op_meta),    INTENT(IN), OPTIONAL :: post_op                      ! "post-op" (small arithmetic operations) for this variable
    TYPE(t_var_action),      INTENT(IN), OPTIONAL :: action_list                  ! regularly triggered events
    INTEGER,                 INTENT(in), OPTIONAL :: var_class                    !< variable type/species
    LOGICAL,                 INTENT(in), OPTIONAL :: lopenacc                     ! create variable on GPU
    ! local variables
    TYPE(t_list_element), POINTER :: element
    INTEGER                       :: idims(5), ndims

    ndims = 3
    idims(1:ndims) = ldims(1:ndims)
    CALL add_var_list_element_5d(ndims, REAL_T, this_list, name,               &
      &   hgrid, vgrid, cf, grib2, idims, element, loutput, lcontainer,        &
      &   lrestart, lrestart_cont, isteptype, lmiss, tlev_source,              &
      &   info, vert_interp, hor_interp, in_group, verbose,                    &
      &   l_pp_scheduler_task, post_op, action_list, tracer_info=tracer_info,  &
      &   p5_r=p5, initval_r=initval, resetval_r=resetval, missval_r=missval,  &
      &   var_class=var_class, lopenacc=lopenacc)
    ptr => element%field%r_ptr(:,:,:,1,1)
    IF (PRESENT(new_element))  new_element => element
  END SUBROUTINE add_var_list_element_r3d


  !------------------------------------------------------------------------------------------------
  ! create (allocate) a new table entry
  ! optionally obtain pointer to 2d-field
  ! optionally overwrite default meta data
  !
  SUBROUTINE add_var_list_element_r2d(this_list, name, ptr,            &
    &   hgrid, vgrid, cf, grib2, ldims, loutput, lcontainer,           &
    &   lrestart, lrestart_cont, initval, isteptype,                   &
    &   resetval, lmiss, missval, tlev_source, tracer_info, info, p5,  &
    &   vert_interp, hor_interp, in_group, verbose, new_element,       &
    &   l_pp_scheduler_task, post_op, action_list, var_class, lopenacc)
    !
    TYPE(t_var_list),        INTENT(inout)        :: this_list                    ! list
    CHARACTER(len=*),        INTENT(in)           :: name                         ! name of variable
    REAL(dp),       POINTER, INTENT(OUT)          :: ptr(:,:)                     ! reference to field
    INTEGER,                 INTENT(in)           :: hgrid                        ! horizontal grid type used
    INTEGER,                 INTENT(in)           :: vgrid                        ! vertical grid type used
    TYPE(t_cf_var),          INTENT(in)           :: cf                           ! CF related metadata
    TYPE(t_grib2_var),       INTENT(in)           :: grib2                        ! GRIB2 related metadata
    INTEGER,                 INTENT(in)           :: ldims(:)                     ! local dimensions
    LOGICAL,                 INTENT(in), OPTIONAL :: loutput                      ! output flag
    LOGICAL,                 INTENT(in), OPTIONAL :: lcontainer                   ! container flag
    LOGICAL,                 INTENT(in), OPTIONAL :: lrestart                     ! restart flag
    LOGICAL,                 INTENT(in), OPTIONAL :: lrestart_cont                ! continue restart if var not available
    REAL(dp),                INTENT(in), OPTIONAL :: initval                      ! value if var not available
    INTEGER,                 INTENT(in), OPTIONAL :: isteptype                    ! type of statistical processing
    REAL(dp),                INTENT(in), OPTIONAL :: resetval                     ! reset value (after accumulation)
    LOGICAL,                 INTENT(in), OPTIONAL :: lmiss                        ! missing value flag
    REAL(dp),                INTENT(in), OPTIONAL :: missval                      ! missing value
    INTEGER,                 INTENT(in), OPTIONAL :: tlev_source                  ! actual TL for TL dependent vars
    CLASS(t_tracer_meta),    INTENT(in), OPTIONAL :: tracer_info                  ! tracer meta data
    TYPE(t_var_metadata),    POINTER,    OPTIONAL :: info                         ! returns reference to metadata
    REAL(dp),         CONTIGUOUS_TARGET, OPTIONAL :: p5(:,:,:,:,:)                ! provided pointer
    TYPE(t_vert_interp_meta),INTENT(in), OPTIONAL :: vert_interp                  ! vertical interpolation metadata
    TYPE(t_hor_interp_meta), INTENT(in), OPTIONAL :: hor_interp                   ! horizontal interpolation metadata
    LOGICAL,                 INTENT(in), OPTIONAL :: in_group(:)                  ! groups to which a variable belongs
    LOGICAL,                 INTENT(in), OPTIONAL :: verbose                      ! print information
    TYPE(t_list_element),    POINTER,    OPTIONAL :: new_element                  ! pointer to new var list element
    INTEGER,                 INTENT(in), OPTIONAL :: l_pp_scheduler_task          ! .TRUE., if field is updated by pp scheduler
    TYPE(t_post_op_meta),    INTENT(IN), OPTIONAL :: post_op                      ! "post-op" (small arithmetic operations) for this variable
    TYPE(t_var_action),      INTENT(IN), OPTIONAL :: action_list                  ! regularly triggered events
    INTEGER,                 INTENT(in), OPTIONAL :: var_class                    !< variable type/species
    LOGICAL,                 INTENT(in), OPTIONAL :: lopenacc                     ! create variable on GPU
    ! local variables
    TYPE(t_list_element), POINTER :: element
    INTEGER                       :: idims(5), ndims

    ndims = 2
    idims(1:ndims) = ldims(1:ndims)
    CALL add_var_list_element_5d(ndims, REAL_T, this_list, name,               &
      &   hgrid, vgrid, cf, grib2, idims, element, loutput, lcontainer,        &
      &   lrestart, lrestart_cont, isteptype, lmiss, tlev_source,              &
      &   info, vert_interp, hor_interp, in_group, verbose,                    &
      &   l_pp_scheduler_task, post_op, action_list, tracer_info=tracer_info,  &
      &   p5_r=p5, initval_r=initval, resetval_r=resetval, missval_r=missval,  &
      &   var_class=var_class, lopenacc=lopenacc)
    ptr => element%field%r_ptr(:,:,1,1,1)
    IF (PRESENT(new_element))  new_element => element
  END SUBROUTINE add_var_list_element_r2d


  !------------------------------------------------------------------------------------------------
  ! create (allocate) a new table entry
  ! optionally obtain pointer to 1d-field
  ! optionally overwrite default meta data
  !
  SUBROUTINE add_var_list_element_r1d(this_list, name, ptr,       &
    &   hgrid, vgrid, cf, grib2, ldims, loutput, lcontainer,      &
    &   lrestart, lrestart_cont, initval, isteptype,              &
    &   resetval, lmiss, missval, tlev_source, info, p5,          &
    &   vert_interp, hor_interp, in_group, verbose, new_element,  &
    &   l_pp_scheduler_task, post_op, action_list, var_class,     &
    &   lopenacc)
    !
    TYPE(t_var_list),        INTENT(inout)        :: this_list                    ! list
    CHARACTER(len=*),        INTENT(in)           :: name                         ! name of variable
    REAL(dp),       POINTER, INTENT(OUT)          :: ptr(:)                       ! reference to field
    INTEGER,                 INTENT(in)           :: hgrid                        ! horizontal grid type used
    INTEGER,                 INTENT(in)           :: vgrid                        ! vertical grid type used
    TYPE(t_cf_var),          INTENT(in)           :: cf                           ! CF related metadata
    TYPE(t_grib2_var),       INTENT(in)           :: grib2                        ! GRIB2 related metadata
    INTEGER,                 INTENT(in)           :: ldims(:)                     ! local dimensions
    LOGICAL,                 INTENT(in), OPTIONAL :: loutput                      ! output flag
    LOGICAL,                 INTENT(in), OPTIONAL :: lcontainer                   ! container flag
    LOGICAL,                 INTENT(in), OPTIONAL :: lrestart                     ! restart flag
    LOGICAL,                 INTENT(in), OPTIONAL :: lrestart_cont                ! continue restart if var not available
    REAL(dp),                INTENT(in), OPTIONAL :: initval                      ! value if var not available
    INTEGER,                 INTENT(in), OPTIONAL :: isteptype                    ! type of statistical processing
    REAL(dp),                INTENT(in), OPTIONAL :: resetval                     ! reset value (after accumulation)
    LOGICAL,                 INTENT(in), OPTIONAL :: lmiss                        ! missing value flag
    REAL(dp),                INTENT(in), OPTIONAL :: missval                      ! missing value
    INTEGER,                 INTENT(in), OPTIONAL :: tlev_source                  ! actual TL for TL dependent vars
    TYPE(t_var_metadata),    POINTER,    OPTIONAL :: info                         ! returns reference to metadata
    REAL(dp),         CONTIGUOUS_TARGET, OPTIONAL :: p5(:,:,:,:,:)                ! provided pointer
    TYPE(t_vert_interp_meta),INTENT(in), OPTIONAL :: vert_interp                  ! vertical interpolation metadata
    TYPE(t_hor_interp_meta), INTENT(in), OPTIONAL :: hor_interp                   ! horizontal interpolation metadata
    LOGICAL,                 INTENT(in), OPTIONAL :: in_group(:)                  ! groups to which a variable belongs
    LOGICAL,                 INTENT(in), OPTIONAL :: verbose                      ! print information
    TYPE(t_list_element),    POINTER,    OPTIONAL :: new_element                  ! pointer to new var list element
    INTEGER,                 INTENT(in), OPTIONAL :: l_pp_scheduler_task          ! .TRUE., if field is updated by pp scheduler
    TYPE(t_post_op_meta),    INTENT(IN), OPTIONAL :: post_op                      ! "post-op" (small arithmetic operations) for this variable
    TYPE(t_var_action),      INTENT(IN), OPTIONAL :: action_list                  ! regularly triggered events
    INTEGER,                 INTENT(in), OPTIONAL :: var_class                    !< variable type/species
    LOGICAL,                 INTENT(in), OPTIONAL :: lopenacc                     ! create variable on GPU
    ! local variables
    TYPE(t_list_element), POINTER :: element
    INTEGER                       :: idims(5), ndims

    ndims = 1
    idims(1:ndims) = ldims(1:ndims)
    CALL add_var_list_element_5d(ndims, REAL_T, this_list, name,        &
      &   hgrid, vgrid, cf, grib2, idims, element, loutput, lcontainer, &
      &   lrestart, lrestart_cont, isteptype, lmiss, tlev_source, info, &
      &   vert_interp, hor_interp, in_group, verbose,                   &
      &   l_pp_scheduler_task, post_op, action_list, p5_r=p5,           &
      &   initval_r=initval, resetval_r=resetval, missval_r=missval,    &
      &   var_class=var_class,lopenacc=lopenacc)
    ptr => element%field%r_ptr(:,1,1,1,1)
    IF (PRESENT(new_element))  new_element => element

  END SUBROUTINE add_var_list_element_r1d

  !------------------------------------------------------------------------------------------------
  ! create (allocate) a new table entry
  ! optionally obtain pointer to 4d-field
  ! optionally overwrite default meta data
  !
  SUBROUTINE add_var_list_element_s4d(this_list, name, ptr,       &
    &   hgrid, vgrid, cf, grib2, ldims, loutput, lcontainer,      &
    &   lrestart, lrestart_cont, initval, isteptype,              &
    &   resetval, lmiss, missval, tlev_source, info, p5,          &
    &   vert_interp, hor_interp, in_group, verbose, new_element,  &
    &   l_pp_scheduler_task, post_op, action_list, var_class,     &
    &   lopenacc)
    !
    TYPE(t_var_list),        INTENT(inout)        :: this_list                    ! list
    CHARACTER(len=*),        INTENT(in)           :: name                         ! name of variable
    REAL(sp),       POINTER, INTENT(OUT)          :: ptr(:,:,:,:)                 ! reference to field
    INTEGER,                 INTENT(in)           :: hgrid                        ! horizontal grid type used
    INTEGER,                 INTENT(in)           :: vgrid                        ! vertical grid type used
    TYPE(t_cf_var),          INTENT(in)           :: cf                           ! CF related metadata
    TYPE(t_grib2_var),       INTENT(in)           :: grib2                        ! GRIB2 related metadata
    INTEGER,                 INTENT(in)           :: ldims(:)                     ! local dimensions
    LOGICAL,                 INTENT(in), OPTIONAL :: loutput                      ! output flag
    LOGICAL,                 INTENT(in), OPTIONAL :: lcontainer                   ! container flag
    LOGICAL,                 INTENT(in), OPTIONAL :: lrestart                     ! restart flag
    LOGICAL,                 INTENT(in), OPTIONAL :: lrestart_cont                ! continue restart if var not available
    REAL(sp),                INTENT(in), OPTIONAL :: initval                      ! value if var not available
    INTEGER,                 INTENT(in), OPTIONAL :: isteptype                    ! type of statistical processing
    REAL(sp),                INTENT(in), OPTIONAL :: resetval                     ! reset value (after accumulation)
    LOGICAL,                 INTENT(in), OPTIONAL :: lmiss                        ! missing value flag
    REAL(sp),                INTENT(in), OPTIONAL :: missval                      ! missing value
    INTEGER,                 INTENT(in), OPTIONAL :: tlev_source                  ! actual TL for TL dependent vars
    TYPE(t_var_metadata),    POINTER,    OPTIONAL :: info                         ! returns reference to metadata
    REAL(sp),         CONTIGUOUS_TARGET, OPTIONAL :: p5(:,:,:,:,:)                ! provided pointer
    TYPE(t_vert_interp_meta),INTENT(in), OPTIONAL :: vert_interp                  ! vertical interpolation metadata
    TYPE(t_hor_interp_meta), INTENT(in), OPTIONAL :: hor_interp                   ! horizontal interpolation metadata
    LOGICAL,                 INTENT(in), OPTIONAL :: in_group(:)                  ! groups to which a variable belongs
    LOGICAL,                 INTENT(in), OPTIONAL :: verbose                      ! print information
    TYPE(t_list_element),    POINTER,    OPTIONAL :: new_element                  ! pointer to new var list element
    INTEGER,                 INTENT(in), OPTIONAL :: l_pp_scheduler_task          ! .TRUE., if field is updated by pp scheduler
    TYPE(t_post_op_meta),    INTENT(IN), OPTIONAL :: post_op                      ! "post-op" (small arithmetic operations) for this variable
    TYPE(t_var_action),      INTENT(IN), OPTIONAL :: action_list                  ! regularly triggered events
    INTEGER,                 INTENT(in), OPTIONAL :: var_class                    !< variable type/species
    LOGICAL,                 INTENT(in), OPTIONAL :: lopenacc                     ! create variable on GPU
    ! local variables
    TYPE(t_list_element), POINTER :: element
    INTEGER                       :: idims(5), ndims

    ndims = 4
    idims(1:ndims) = ldims(1:ndims)
    CALL add_var_list_element_5d(ndims, SINGLE_T, this_list, name,         &
      &   hgrid, vgrid, cf, grib2, idims, element, loutput, lcontainer,    &
      &   lrestart, lrestart_cont, isteptype, lmiss, tlev_source, info,    &
      &   vert_interp, hor_interp, in_group, verbose,                      &
      &   l_pp_scheduler_task, post_op, action_list, p5_s=p5,              &
      &   initval_s=initval, resetval_s=resetval, missval_s=missval,       &
      &   var_class=var_class, lopenacc=lopenacc)
    ptr => element%field%s_ptr(:,:,:,:,1)
    IF (PRESENT(new_element))  new_element => element
  END SUBROUTINE add_var_list_element_s4d


  !------------------------------------------------------------------------------------------------
  ! create (allocate) a new table entry
  ! optionally obtain pointer to 3d-field
  ! optionally overwrite default meta data
  !
  SUBROUTINE add_var_list_element_s3d(this_list, name, ptr,            &
    &   hgrid, vgrid, cf, grib2, ldims, loutput, lcontainer,           &
    &   lrestart, lrestart_cont, initval, isteptype,                   &
    &   resetval, lmiss, missval, tlev_source, tracer_info, info, p5,  &
    &   vert_interp, hor_interp, in_group, verbose, new_element,       &
    &   l_pp_scheduler_task, post_op, action_list, var_class, lopenacc)
    !
    TYPE(t_var_list),        INTENT(inout)        :: this_list                    ! list
    CHARACTER(len=*),        INTENT(in)           :: name                         ! name of variable
    REAL(sp),       POINTER, INTENT(OUT)          :: ptr(:,:,:)                   ! reference to field
    INTEGER,                 INTENT(in)           :: hgrid                        ! horizontal grid type used
    INTEGER,                 INTENT(in)           :: vgrid                        ! vertical grid type used
    TYPE(t_cf_var),          INTENT(in)           :: cf                           ! CF related metadata
    TYPE(t_grib2_var),       INTENT(in)           :: grib2                        ! GRIB2 related metadata
    INTEGER,                 INTENT(in)           :: ldims(:)                     ! local dimensions
    LOGICAL,                 INTENT(in), OPTIONAL :: loutput                      ! output flag
    LOGICAL,                 INTENT(in), OPTIONAL :: lcontainer                   ! container flag
    LOGICAL,                 INTENT(in), OPTIONAL :: lrestart                     ! restart flag
    LOGICAL,                 INTENT(in), OPTIONAL :: lrestart_cont                ! continue restart if var not available
    REAL(sp),                INTENT(in), OPTIONAL :: initval                      ! value if var not available
    INTEGER,                 INTENT(in), OPTIONAL :: isteptype                    ! type of statistical processing
    REAL(sp),                INTENT(in), OPTIONAL :: resetval                     ! reset value (after accumulation)
    LOGICAL,                 INTENT(in), OPTIONAL :: lmiss                        ! missing value flag
    REAL(sp),                INTENT(in), OPTIONAL :: missval                      ! missing value
    INTEGER,                 INTENT(in), OPTIONAL :: tlev_source                  ! actual TL for TL dependent vars
    CLASS(t_tracer_meta),    INTENT(in), OPTIONAL :: tracer_info                  ! tracer meta data
    TYPE(t_var_metadata),    POINTER,    OPTIONAL :: info                         ! returns reference to metadata
    REAL(sp),         CONTIGUOUS_TARGET, OPTIONAL :: p5(:,:,:,:,:)                ! provided pointer
    TYPE(t_vert_interp_meta),INTENT(in), OPTIONAL :: vert_interp                  ! vertical interpolation metadata
    TYPE(t_hor_interp_meta), INTENT(in), OPTIONAL :: hor_interp                   ! horizontal interpolation metadata
    LOGICAL,                 INTENT(in), OPTIONAL :: in_group(:)                  ! groups to which a variable belongs
    LOGICAL,                 INTENT(in), OPTIONAL :: verbose                      ! print information
    TYPE(t_list_element),    POINTER,    OPTIONAL :: new_element                  ! pointer to new var list element
    INTEGER,                 INTENT(in), OPTIONAL :: l_pp_scheduler_task          ! .TRUE., if field is updated by pp scheduler
    TYPE(t_post_op_meta),    INTENT(IN), OPTIONAL :: post_op                      ! "post-op" (small arithmetic operations) for this variable
    TYPE(t_var_action),      INTENT(IN), OPTIONAL :: action_list                  ! regularly triggered events
    INTEGER,                 INTENT(in), OPTIONAL :: var_class                    !< variable type/species
    LOGICAL,                 INTENT(in), OPTIONAL :: lopenacc                     ! create variable on GPU
    ! local variables
    TYPE(t_list_element), POINTER :: element
    INTEGER                       :: idims(5), ndims

    ndims = 3
    idims(1:ndims) = ldims(1:ndims)
    CALL add_var_list_element_5d(ndims, SINGLE_T, this_list, name,             &
      &   hgrid, vgrid, cf, grib2, idims, element, loutput, lcontainer,        &
      &   lrestart, lrestart_cont, isteptype, lmiss, tlev_source,              &
      &   info, vert_interp, hor_interp, in_group, verbose,                    &
      &   l_pp_scheduler_task, post_op, action_list, tracer_info=tracer_info,  &
      &   p5_s=p5, initval_s=initval, resetval_s=resetval, missval_s=missval,  &
      &   var_class=var_class, lopenacc=lopenacc)
    ptr => element%field%s_ptr(:,:,:,1,1)
    IF (PRESENT(new_element))  new_element => element
  END SUBROUTINE add_var_list_element_s3d


  !------------------------------------------------------------------------------------------------
  ! create (allocate) a new table entry
  ! optionally obtain pointer to 2d-field
  ! optionally overwrite default meta data
  !
  SUBROUTINE add_var_list_element_s2d(this_list, name, ptr,            &
    &   hgrid, vgrid, cf, grib2, ldims, loutput, lcontainer,           &
    &   lrestart, lrestart_cont, initval, isteptype,                   &
    &   resetval, lmiss, missval, tlev_source, tracer_info, info, p5,  &
    &   vert_interp, hor_interp, in_group, verbose, new_element,       &
    &   l_pp_scheduler_task, post_op, action_list, var_class, lopenacc)
    !
    TYPE(t_var_list),        INTENT(inout)        :: this_list                    ! list
    CHARACTER(len=*),        INTENT(in)           :: name                         ! name of variable
    REAL(sp),       POINTER, INTENT(OUT)          :: ptr(:,:)                     ! reference to field
    INTEGER,                 INTENT(in)           :: hgrid                        ! horizontal grid type used
    INTEGER,                 INTENT(in)           :: vgrid                        ! vertical grid type used
    TYPE(t_cf_var),          INTENT(in)           :: cf                           ! CF related metadata
    TYPE(t_grib2_var),       INTENT(in)           :: grib2                        ! GRIB2 related metadata
    INTEGER,                 INTENT(in)           :: ldims(:)                     ! local dimensions
    LOGICAL,                 INTENT(in), OPTIONAL :: loutput                      ! output flag
    LOGICAL,                 INTENT(in), OPTIONAL :: lcontainer                   ! container flag
    LOGICAL,                 INTENT(in), OPTIONAL :: lrestart                     ! restart flag
    LOGICAL,                 INTENT(in), OPTIONAL :: lrestart_cont                ! continue restart if var not available
    REAL(sp),                INTENT(in), OPTIONAL :: initval                      ! value if var not available
    INTEGER,                 INTENT(in), OPTIONAL :: isteptype                    ! type of statistical processing
    REAL(sp),                INTENT(in), OPTIONAL :: resetval                     ! reset value (after accumulation)
    LOGICAL,                 INTENT(in), OPTIONAL :: lmiss                        ! missing value flag
    REAL(sp),                INTENT(in), OPTIONAL :: missval                      ! missing value
    INTEGER,                 INTENT(in), OPTIONAL :: tlev_source                  ! actual TL for TL dependent vars
    CLASS(t_tracer_meta),    INTENT(in), OPTIONAL :: tracer_info                  ! tracer meta data
    TYPE(t_var_metadata),    POINTER,    OPTIONAL :: info                         ! returns reference to metadata
    REAL(sp),         CONTIGUOUS_TARGET, OPTIONAL :: p5(:,:,:,:,:)                ! provided pointer
    TYPE(t_vert_interp_meta),INTENT(in), OPTIONAL :: vert_interp                  ! vertical interpolation metadata
    TYPE(t_hor_interp_meta), INTENT(in), OPTIONAL :: hor_interp                   ! horizontal interpolation metadata
    LOGICAL,                 INTENT(in), OPTIONAL :: in_group(:)                  ! groups to which a variable belongs
    LOGICAL,                 INTENT(in), OPTIONAL :: verbose                      ! print information
    TYPE(t_list_element),    POINTER,    OPTIONAL :: new_element                  ! pointer to new var list element
    INTEGER,                 INTENT(in), OPTIONAL :: l_pp_scheduler_task          ! .TRUE., if field is updated by pp scheduler
    TYPE(t_post_op_meta),    INTENT(IN), OPTIONAL :: post_op                      ! "post-op" (small arithmetic operations) for this variable
    TYPE(t_var_action),      INTENT(IN), OPTIONAL :: action_list                  ! regularly triggered events
    INTEGER,                 INTENT(in), OPTIONAL :: var_class                    !< variable type/species
    LOGICAL,                 INTENT(in), OPTIONAL :: lopenacc                     ! create variable on GPU
    ! local variables
    TYPE(t_list_element), POINTER :: element
    INTEGER                       :: idims(5), ndims

    ndims = 2
    idims(1:ndims) = ldims(1:ndims)
    CALL add_var_list_element_5d(ndims, SINGLE_T, this_list, name,             &
      &   hgrid, vgrid, cf, grib2, idims, element, loutput, lcontainer,        &
      &   lrestart, lrestart_cont, isteptype, lmiss, tlev_source,              &
      &   info, vert_interp, hor_interp, in_group, verbose,                    &
      &   l_pp_scheduler_task, post_op, action_list, tracer_info=tracer_info,  &
      &   p5_s=p5, initval_s=initval, resetval_s=resetval, missval_s=missval,  &
      &   var_class=var_class, lopenacc=lopenacc)
    ptr => element%field%s_ptr(:,:,1,1,1)
    IF (PRESENT(new_element))  new_element => element
  END SUBROUTINE add_var_list_element_s2d


  !------------------------------------------------------------------------------------------------
  ! create (allocate) a new table entry
  ! optionally obtain pointer to 1d-field
  ! optionally overwrite default meta data
  !
  SUBROUTINE add_var_list_element_s1d(this_list, name, ptr,       &
    &   hgrid, vgrid, cf, grib2, ldims, loutput, lcontainer,      &
    &   lrestart, lrestart_cont, initval, isteptype,              &
    &   resetval, lmiss, missval, tlev_source, info, p5,          &
    &   vert_interp, hor_interp, in_group, verbose, new_element,  &
    &   l_pp_scheduler_task, post_op, action_list, var_class,     &
    &   lopenacc)
    !
    TYPE(t_var_list),        INTENT(inout)        :: this_list                    ! list
    CHARACTER(len=*),        INTENT(in)           :: name                         ! name of variable
    REAL(sp),       POINTER, INTENT(OUT)          :: ptr(:)                       ! reference to field
    INTEGER,                 INTENT(in)           :: hgrid                        ! horizontal grid type used
    INTEGER,                 INTENT(in)           :: vgrid                        ! vertical grid type used
    TYPE(t_cf_var),          INTENT(in)           :: cf                           ! CF related metadata
    TYPE(t_grib2_var),       INTENT(in)           :: grib2                        ! GRIB2 related metadata
    INTEGER,                 INTENT(in)           :: ldims(:)                     ! local dimensions
    LOGICAL,                 INTENT(in), OPTIONAL :: loutput                      ! output flag
    LOGICAL,                 INTENT(in), OPTIONAL :: lcontainer                   ! container flag
    LOGICAL,                 INTENT(in), OPTIONAL :: lrestart                     ! restart flag
    LOGICAL,                 INTENT(in), OPTIONAL :: lrestart_cont                ! continue restart if var not available
    REAL(sp),                INTENT(in), OPTIONAL :: initval                      ! value if var not available
    INTEGER,                 INTENT(in), OPTIONAL :: isteptype                    ! type of statistical processing
    REAL(sp),                INTENT(in), OPTIONAL :: resetval                     ! reset value (after accumulation)
    LOGICAL,                 INTENT(in), OPTIONAL :: lmiss                        ! missing value flag
    REAL(sp),                INTENT(in), OPTIONAL :: missval                      ! missing value
    INTEGER,                 INTENT(in), OPTIONAL :: tlev_source                  ! actual TL for TL dependent vars
    TYPE(t_var_metadata),    POINTER,    OPTIONAL :: info                         ! returns reference to metadata
    REAL(sp),         CONTIGUOUS_TARGET, OPTIONAL :: p5(:,:,:,:,:)                ! provided pointer
    TYPE(t_vert_interp_meta),INTENT(in), OPTIONAL :: vert_interp                  ! vertical interpolation metadata
    TYPE(t_hor_interp_meta), INTENT(in), OPTIONAL :: hor_interp                   ! horizontal interpolation metadata
    LOGICAL,                 INTENT(in), OPTIONAL :: in_group(:)                  ! groups to which a variable belongs
    LOGICAL,                 INTENT(in), OPTIONAL :: verbose                      ! print information
    TYPE(t_list_element),    POINTER,    OPTIONAL :: new_element                  ! pointer to new var list element
    INTEGER,                 INTENT(in), OPTIONAL :: l_pp_scheduler_task          ! .TRUE., if field is updated by pp scheduler
    TYPE(t_post_op_meta),    INTENT(IN), OPTIONAL :: post_op                      ! "post-op" (small arithmetic operations) for this variable
    TYPE(t_var_action),      INTENT(IN), OPTIONAL :: action_list                  ! regularly triggered events
    INTEGER,                 INTENT(in), OPTIONAL :: var_class                    !< variable type/species
    LOGICAL,                 INTENT(in), OPTIONAL :: lopenacc                     ! create variable on GPU
    ! local variables
    TYPE(t_list_element), POINTER :: element
    INTEGER                       :: idims(5), ndims

    ndims = 1
    idims(1:ndims) = ldims(1:ndims)
    CALL add_var_list_element_5d(ndims, SINGLE_T, this_list, name,      &
      &   hgrid, vgrid, cf, grib2, idims, element, loutput, lcontainer, &
      &   lrestart, lrestart_cont, isteptype, lmiss, tlev_source, info, &
      &   vert_interp, hor_interp, in_group, verbose,                   &
      &   l_pp_scheduler_task, post_op, action_list, p5_s=p5,           &
      &   initval_s=initval, resetval_s=resetval, missval_s=missval,    &
      &   var_class=var_class, lopenacc=lopenacc)
    ptr => element%field%s_ptr(:,1,1,1,1)
    IF (PRESENT(new_element))  new_element => element

  END SUBROUTINE add_var_list_element_s1d

  !------------------------------------------------------------------------------------------------
  ! create (allocate) a new table entry
  ! optionally obtain pointer to 4d-field
  ! optionally overwrite default meta data
  !
  SUBROUTINE add_var_list_element_i4d(this_list, name, ptr,       &
    &   hgrid, vgrid, cf, grib2, ldims, loutput, lcontainer,      &
    &   lrestart, lrestart_cont, initval, isteptype,              &
    &   resetval, lmiss, missval, tlev_source, info, p5,          &
    &   vert_interp, hor_interp, in_group, verbose, new_element,  &
    &   l_pp_scheduler_task, post_op, action_list, var_class,     &
    &   lopenacc)
    !
    TYPE(t_var_list),        INTENT(inout)        :: this_list                    ! list
    CHARACTER(len=*),        INTENT(in)           :: name                         ! name of variable
    INTEGER,        POINTER, INTENT(OUT)          :: ptr(:,:,:,:)                 ! reference to field
    INTEGER,                 INTENT(in)           :: hgrid                        ! horizontal grid type used
    INTEGER,                 INTENT(in)           :: vgrid                        ! vertical grid type used
    TYPE(t_cf_var),          INTENT(in)           :: cf                           ! CF related metadata
    TYPE(t_grib2_var),       INTENT(in)           :: grib2                        ! GRIB2 related metadata
    INTEGER,                 INTENT(in)           :: ldims(:)                     ! local dimensions
    LOGICAL,                 INTENT(in), OPTIONAL :: loutput                      ! output flag
    LOGICAL,                 INTENT(in), OPTIONAL :: lcontainer                   ! container flag
    LOGICAL,                 INTENT(in), OPTIONAL :: lrestart                     ! restart flag
    LOGICAL,                 INTENT(in), OPTIONAL :: lrestart_cont                ! continue restart if var not available
    INTEGER,                 INTENT(in), OPTIONAL :: initval                      ! value if var not available
    INTEGER,                 INTENT(in), OPTIONAL :: isteptype                    ! type of statistical processing
    INTEGER,                 INTENT(in), OPTIONAL :: resetval                     ! reset value (after accumulation)
    LOGICAL,                 INTENT(in), OPTIONAL :: lmiss                        ! missing value flag
    INTEGER,                 INTENT(in), OPTIONAL :: missval                      ! missing value
    INTEGER,                 INTENT(in), OPTIONAL :: tlev_source                  ! actual TL for TL dependent vars
    TYPE(t_var_metadata),    POINTER,    OPTIONAL :: info                         ! returns reference to metadata
    INTEGER,          CONTIGUOUS_TARGET, OPTIONAL :: p5(:,:,:,:,:)                ! provided pointer
    TYPE(t_vert_interp_meta),INTENT(in), OPTIONAL :: vert_interp                  ! vertical interpolation metadata
    TYPE(t_hor_interp_meta), INTENT(in), OPTIONAL :: hor_interp                   ! horizontal interpolation metadata
    LOGICAL,                 INTENT(in), OPTIONAL :: in_group(:)                  ! groups to which a variable belongs
    LOGICAL,                 INTENT(in), OPTIONAL :: verbose                      ! print information
    TYPE(t_list_element),    POINTER,    OPTIONAL :: new_element                  ! pointer to new var list element
    INTEGER,                 INTENT(in), OPTIONAL :: l_pp_scheduler_task          ! .TRUE., if field is updated by pp scheduler
    TYPE(t_post_op_meta),    INTENT(IN), OPTIONAL :: post_op                      ! "post-op" (small arithmetic operations) for this variable
    TYPE(t_var_action),      INTENT(IN), OPTIONAL :: action_list                  ! regularly triggered events
    INTEGER,                 INTENT(in), OPTIONAL :: var_class                    !< variable type/species
    LOGICAL,                 INTENT(in), OPTIONAL :: lopenacc                     ! create variable on GPU
    ! local variables
    TYPE(t_list_element), POINTER :: element
    INTEGER                       :: idims(5), ndims

    ndims = 4
    idims(1:ndims) = ldims(1:ndims)
    CALL add_var_list_element_5d(ndims, INT_T, this_list, name,         &
      &   hgrid, vgrid, cf, grib2, idims, element, loutput, lcontainer, &
      &   lrestart, lrestart_cont, isteptype, lmiss, tlev_source, info, &
      &   vert_interp, hor_interp, in_group, verbose,                   &
      &   l_pp_scheduler_task, post_op, action_list, p5_i=p5,           &
      &   initval_i=initval, resetval_i=resetval, missval_i=missval,    &
      &   var_class=var_class, lopenacc=lopenacc)
    ptr => element%field%i_ptr(:,:,:,:,1)
    IF (PRESENT(new_element))  new_element => element
  END SUBROUTINE add_var_list_element_i4d


  !------------------------------------------------------------------------------------------------
  ! create (allocate) a new table entry
  ! optionally obtain pointer to 3d-field
  ! optionally overwrite default meta data
  !
  SUBROUTINE add_var_list_element_i3d(this_list, name, ptr,       &
    &   hgrid, vgrid, cf, grib2, ldims, loutput, lcontainer,      &
    &   lrestart, lrestart_cont, initval, isteptype,              &
    &   resetval, lmiss, missval, tlev_source, info, p5,          &
    &   vert_interp, hor_interp, in_group, verbose, new_element,  &
    &   l_pp_scheduler_task, post_op, action_list, var_class,     &
    &   lopenacc)
    !
    TYPE(t_var_list),        INTENT(inout)        :: this_list                    ! list
    CHARACTER(len=*),        INTENT(in)           :: name                         ! name of variable
    INTEGER,        POINTER, INTENT(OUT)          :: ptr(:,:,:)                   ! reference to field
    INTEGER,                 INTENT(in)           :: hgrid                        ! horizontal grid type used
    INTEGER,                 INTENT(in)           :: vgrid                        ! vertical grid type used
    TYPE(t_cf_var),          INTENT(in)           :: cf                           ! CF related metadata
    TYPE(t_grib2_var),       INTENT(in)           :: grib2                        ! GRIB2 related metadata
    INTEGER,                 INTENT(in)           :: ldims(:)                     ! local dimensions
    LOGICAL,                 INTENT(in), OPTIONAL :: loutput                      ! output flag
    LOGICAL,                 INTENT(in), OPTIONAL :: lcontainer                   ! container flag
    LOGICAL,                 INTENT(in), OPTIONAL :: lrestart                     ! restart flag
    LOGICAL,                 INTENT(in), OPTIONAL :: lrestart_cont                ! continue restart if var not available
    INTEGER,                 INTENT(in), OPTIONAL :: initval                      ! value if var not available
    INTEGER,                 INTENT(in), OPTIONAL :: isteptype                    ! type of statistical processing
    INTEGER,                 INTENT(in), OPTIONAL :: resetval                     ! reset value (after accumulation)
    LOGICAL,                 INTENT(in), OPTIONAL :: lmiss                        ! missing value flag
    INTEGER,                 INTENT(in), OPTIONAL :: missval                      ! missing value
    INTEGER,                 INTENT(in), OPTIONAL :: tlev_source                  ! actual TL for TL dependent vars
    TYPE(t_var_metadata),    POINTER,    OPTIONAL :: info                         ! returns reference to metadata
    INTEGER,          CONTIGUOUS_TARGET, OPTIONAL :: p5(:,:,:,:,:)                ! provided pointer
    TYPE(t_vert_interp_meta),INTENT(in), OPTIONAL :: vert_interp                  ! vertical interpolation metadata
    TYPE(t_hor_interp_meta), INTENT(in), OPTIONAL :: hor_interp                   ! horizontal interpolation metadata
    LOGICAL,                 INTENT(in), OPTIONAL :: in_group(:)                  ! groups to which a variable belongs
    LOGICAL,                 INTENT(in), OPTIONAL :: verbose                      ! print information
    TYPE(t_list_element),    POINTER,    OPTIONAL :: new_element                  ! pointer to new var list element
    INTEGER,                 INTENT(in), OPTIONAL :: l_pp_scheduler_task          ! .TRUE., if field is updated by pp scheduler
    TYPE(t_post_op_meta),    INTENT(IN), OPTIONAL :: post_op                      ! "post-op" (small arithmetic operations) for this variable
    TYPE(t_var_action),      INTENT(IN), OPTIONAL :: action_list                  ! regularly triggered events
    INTEGER,                 INTENT(in), OPTIONAL :: var_class                    !< variable type/species
    LOGICAL,                 INTENT(in), OPTIONAL :: lopenacc                     ! create variable on GPU
    ! local variables
    TYPE(t_list_element), POINTER :: element
    INTEGER                       :: idims(5), ndims

    ndims = 3
    idims(1:ndims) = ldims(1:ndims)
    CALL add_var_list_element_5d(ndims, INT_T, this_list, name,         &
      &   hgrid, vgrid, cf, grib2, idims, element, loutput, lcontainer, &
      &   lrestart, lrestart_cont, isteptype, lmiss, tlev_source, info, &
      &   vert_interp, hor_interp, in_group, verbose,                   &
      &   l_pp_scheduler_task, post_op, action_list, p5_i=p5,           &
      &   initval_i=initval, resetval_i=resetval, missval_i=missval,    &
      &   var_class=var_class, lopenacc=lopenacc)
    ptr => element%field%i_ptr(:,:,:,1,1)
    IF (PRESENT(new_element))  new_element => element
  END SUBROUTINE add_var_list_element_i3d


  !------------------------------------------------------------------------------------------------
  ! create (allocate) a new table entry
  ! optionally obtain pointer to 2d-field
  ! optionally overwrite default meta data
  !
  SUBROUTINE add_var_list_element_i2d(this_list, name, ptr,       &
    &   hgrid, vgrid, cf, grib2, ldims, loutput, lcontainer,      &
    &   lrestart, lrestart_cont, initval, isteptype,              &
    &   resetval, lmiss, missval, tlev_source, info, p5,          &
    &   vert_interp, hor_interp, in_group, verbose, new_element,  &
    &   l_pp_scheduler_task, post_op, action_list, var_class,     &
    &   lopenacc)
    !
    TYPE(t_var_list),        INTENT(inout)        :: this_list                    ! list
    CHARACTER(len=*),        INTENT(in)           :: name                         ! name of variable
    INTEGER,        POINTER, INTENT(OUT)          :: ptr(:,:)                     ! reference to field
    INTEGER,                 INTENT(in)           :: hgrid                        ! horizontal grid type used
    INTEGER,                 INTENT(in)           :: vgrid                        ! vertical grid type used
    TYPE(t_cf_var),          INTENT(in)           :: cf                           ! CF related metadata
    TYPE(t_grib2_var),       INTENT(in)           :: grib2                        ! GRIB2 related metadata
    INTEGER,                 INTENT(in)           :: ldims(:)                     ! local dimensions
    LOGICAL,                 INTENT(in), OPTIONAL :: loutput                      ! output flag
    LOGICAL,                 INTENT(in), OPTIONAL :: lcontainer                   ! container flag
    LOGICAL,                 INTENT(in), OPTIONAL :: lrestart                     ! restart flag
    LOGICAL,                 INTENT(in), OPTIONAL :: lrestart_cont                ! continue restart if var not available
    INTEGER,                 INTENT(in), OPTIONAL :: initval                      ! value if var not available
    INTEGER,                 INTENT(in), OPTIONAL :: isteptype                    ! type of statistical processing
    INTEGER,                 INTENT(in), OPTIONAL :: resetval                     ! reset value (after accumulation)
    LOGICAL,                 INTENT(in), OPTIONAL :: lmiss                        ! missing value flag
    INTEGER,                 INTENT(in), OPTIONAL :: missval                      ! missing value
    INTEGER,                 INTENT(in), OPTIONAL :: tlev_source                  ! actual TL for TL dependent vars
    TYPE(t_var_metadata),    POINTER,    OPTIONAL :: info                         ! returns reference to metadata
    INTEGER,          CONTIGUOUS_TARGET, OPTIONAL :: p5(:,:,:,:,:)                ! provided pointer
    TYPE(t_vert_interp_meta),INTENT(in), OPTIONAL :: vert_interp                  ! vertical interpolation metadata
    TYPE(t_hor_interp_meta), INTENT(in), OPTIONAL :: hor_interp                   ! horizontal interpolation metadata
    LOGICAL,                 INTENT(in), OPTIONAL :: in_group(:)                  ! groups to which a variable belongs
    LOGICAL,                 INTENT(in), OPTIONAL :: verbose                      ! print information
    TYPE(t_list_element),    POINTER,    OPTIONAL :: new_element                  ! pointer to new var list element
    INTEGER,                 INTENT(in), OPTIONAL :: l_pp_scheduler_task          ! .TRUE., if field is updated by pp scheduler
    TYPE(t_post_op_meta),    INTENT(IN), OPTIONAL :: post_op                      ! "post-op" (small arithmetic operations) for this variable
    TYPE(t_var_action),      INTENT(IN), OPTIONAL :: action_list                  ! regularly triggered events
    INTEGER,                 INTENT(in), OPTIONAL :: var_class                    !< variable type/species
    LOGICAL,                 INTENT(in), OPTIONAL :: lopenacc                     ! create variable on GPU
    ! local variables
    TYPE(t_list_element), POINTER :: element
    INTEGER                       :: idims(5), ndims

    ndims = 2
    idims(1:ndims) = ldims(1:ndims)
    CALL add_var_list_element_5d(ndims, INT_T, this_list, name,         &
      &   hgrid, vgrid, cf, grib2, idims, element, loutput, lcontainer, &
      &   lrestart, lrestart_cont, isteptype, lmiss, tlev_source, info, &
      &   vert_interp, hor_interp, in_group, verbose,                   &
      &   l_pp_scheduler_task, post_op, action_list, p5_i=p5,           &
      &   initval_i=initval, resetval_i=resetval, missval_i=missval,    &
      &   var_class=var_class, lopenacc=lopenacc)
    ptr => element%field%i_ptr(:,:,1,1,1)
    IF (PRESENT(new_element))  new_element => element
  END SUBROUTINE add_var_list_element_i2d


  !------------------------------------------------------------------------------------------------
  ! create (allocate) a new table entry
  ! optionally obtain pointer to 1d-field
  ! optionally overwrite default meta data
  !
  SUBROUTINE add_var_list_element_i1d(this_list, name, ptr,       &
    &   hgrid, vgrid, cf, grib2, ldims, loutput, lcontainer,      &
    &   lrestart, lrestart_cont, initval, isteptype,              &
    &   resetval, lmiss, missval, tlev_source, info, p5,          &
    &   vert_interp, hor_interp, in_group, verbose, new_element,  &
    &   l_pp_scheduler_task, post_op, action_list, var_class,     &
    &   lopenacc)
    !
    TYPE(t_var_list),        INTENT(inout)        :: this_list                    ! list
    CHARACTER(len=*),        INTENT(in)           :: name                         ! name of variable
    INTEGER,        POINTER, INTENT(OUT)          :: ptr(:)                       ! reference to field
    INTEGER,                 INTENT(in)           :: hgrid                        ! horizontal grid type used
    INTEGER,                 INTENT(in)           :: vgrid                        ! vertical grid type used
    TYPE(t_cf_var),          INTENT(in)           :: cf                           ! CF related metadata
    TYPE(t_grib2_var),       INTENT(in)           :: grib2                        ! GRIB2 related metadata
    INTEGER,                 INTENT(in)           :: ldims(:)                     ! local dimensions
    LOGICAL,                 INTENT(in), OPTIONAL :: loutput                      ! output flag
    LOGICAL,                 INTENT(in), OPTIONAL :: lcontainer                   ! container flag
    LOGICAL,                 INTENT(in), OPTIONAL :: lrestart                     ! restart flag
    LOGICAL,                 INTENT(in), OPTIONAL :: lrestart_cont                ! continue restart if var not available
    INTEGER,                 INTENT(in), OPTIONAL :: initval                      ! value if var not available
    INTEGER,                 INTENT(in), OPTIONAL :: isteptype                    ! type of statistical processing
    INTEGER,                 INTENT(in), OPTIONAL :: resetval                     ! reset value (after accumulation)
    LOGICAL,                 INTENT(in), OPTIONAL :: lmiss                        ! missing value flag
    INTEGER,                 INTENT(in), OPTIONAL :: missval                      ! missing value
    INTEGER,                 INTENT(in), OPTIONAL :: tlev_source                  ! actual TL for TL dependent vars
    TYPE(t_var_metadata),    POINTER,    OPTIONAL :: info                         ! returns reference to metadata
    INTEGER,          CONTIGUOUS_TARGET, OPTIONAL :: p5(:,:,:,:,:)                ! provided pointer
    TYPE(t_vert_interp_meta),INTENT(in), OPTIONAL :: vert_interp                  ! vertical interpolation metadata
    TYPE(t_hor_interp_meta), INTENT(in), OPTIONAL :: hor_interp                   ! horizontal interpolation metadata
    LOGICAL,                 INTENT(in), OPTIONAL :: in_group(:)                  ! groups to which a variable belongs
    LOGICAL,                 INTENT(in), OPTIONAL :: verbose                      ! print information
    TYPE(t_list_element),    POINTER,    OPTIONAL :: new_element                  ! pointer to new var list element
    INTEGER,                 INTENT(in), OPTIONAL :: l_pp_scheduler_task          ! .TRUE., if field is updated by pp scheduler
    TYPE(t_post_op_meta),    INTENT(IN), OPTIONAL :: post_op                      ! "post-op" (small arithmetic operations) for this variable
    TYPE(t_var_action),      INTENT(IN), OPTIONAL :: action_list                  ! regularly triggered events
    INTEGER,                 INTENT(in), OPTIONAL :: var_class                    !< variable type/species
    LOGICAL,                 INTENT(in), OPTIONAL :: lopenacc                     ! create variable on GPU
    ! local variables
    TYPE(t_list_element), POINTER :: element
    INTEGER                       :: idims(5), ndims

    ndims = 1
    idims(1:ndims) = ldims(1:ndims)
    CALL add_var_list_element_5d(ndims, INT_T, this_list, name,         &
      &   hgrid, vgrid, cf, grib2, idims, element, loutput, lcontainer, &
      &   lrestart, lrestart_cont, isteptype, lmiss, tlev_source, info, &
      &   vert_interp, hor_interp, in_group, verbose,                   &
      &   l_pp_scheduler_task, post_op, action_list, p5_i=p5,           &
      &   initval_i=initval, resetval_i=resetval, missval_i=missval,    &
      &   var_class=var_class, lopenacc=lopenacc)
    ptr => element%field%i_ptr(:,1,1,1,1)
    IF (PRESENT(new_element))  new_element => element
  END SUBROUTINE add_var_list_element_i1d

  !------------------------------------------------------------------------------------------------
  ! create (allocate) a new table entry
  ! optionally obtain pointer to 4d-field
  ! optionally overwrite default meta data
  !
  SUBROUTINE add_var_list_element_l4d(this_list, name, ptr,       &
    &   hgrid, vgrid, cf, grib2, ldims, loutput, lcontainer,      &
    &   lrestart, lrestart_cont, initval, isteptype,              &
    &   resetval, lmiss, missval, tlev_source, info, p5,          &
    &   vert_interp, hor_interp, in_group, verbose, new_element,  &
    &   l_pp_scheduler_task, post_op, action_list, var_class,     &
    &   lopenacc)
    !
    TYPE(t_var_list),        INTENT(inout)        :: this_list                    ! list
    CHARACTER(len=*),        INTENT(in)           :: name                         ! name of variable
    LOGICAL,        POINTER, INTENT(OUT)          :: ptr(:,:,:,:)                 ! reference to field
    INTEGER,                 INTENT(in)           :: hgrid                        ! horizontal grid type used
    INTEGER,                 INTENT(in)           :: vgrid                        ! vertical grid type used
    TYPE(t_cf_var),          INTENT(in)           :: cf                           ! CF related metadata
    TYPE(t_grib2_var),       INTENT(in)           :: grib2                        ! GRIB2 related metadata
    INTEGER,                 INTENT(in)           :: ldims(:)                     ! local dimensions
    LOGICAL,                 INTENT(in), OPTIONAL :: loutput                      ! output flag
    LOGICAL,                 INTENT(in), OPTIONAL :: lcontainer                   ! container flag
    LOGICAL,                 INTENT(in), OPTIONAL :: lrestart                     ! restart flag
    LOGICAL,                 INTENT(in), OPTIONAL :: lrestart_cont                ! continue restart if var not available
    LOGICAL,                 INTENT(in), OPTIONAL :: initval                      ! value if var not available
    INTEGER,                 INTENT(in), OPTIONAL :: isteptype                    ! type of statistical processing
    LOGICAL,                 INTENT(in), OPTIONAL :: resetval                     ! reset value (after accumulation)
    LOGICAL,                 INTENT(in), OPTIONAL :: lmiss                        ! missing value flag
    LOGICAL,                 INTENT(in), OPTIONAL :: missval                      ! missing value
    INTEGER,                 INTENT(in), OPTIONAL :: tlev_source                  ! actual TL for TL dependent vars
    TYPE(t_var_metadata),    POINTER,    OPTIONAL :: info                         ! returns reference to metadata
    LOGICAL,          CONTIGUOUS_TARGET, OPTIONAL :: p5(:,:,:,:,:)                ! provided pointer
    TYPE(t_vert_interp_meta),INTENT(in), OPTIONAL :: vert_interp                  ! vertical interpolation metadata
    TYPE(t_hor_interp_meta), INTENT(in), OPTIONAL :: hor_interp                   ! horizontal interpolation metadata
    LOGICAL,                 INTENT(in), OPTIONAL :: in_group(:)                  ! groups to which a variable belongs
    LOGICAL,                 INTENT(in), OPTIONAL :: verbose                      ! print information
    TYPE(t_list_element),    POINTER,    OPTIONAL :: new_element                  ! pointer to new var list element
    INTEGER,                 INTENT(in), OPTIONAL :: l_pp_scheduler_task          ! .TRUE., if field is updated by pp scheduler
    TYPE(t_post_op_meta),    INTENT(IN), OPTIONAL :: post_op                      ! "post-op" (small arithmetic operations) for this variable
    TYPE(t_var_action),      INTENT(IN), OPTIONAL :: action_list                  ! regularly triggered events
    INTEGER,                 INTENT(in), OPTIONAL :: var_class                    !< variable type/species
    LOGICAL,                 INTENT(in), OPTIONAL :: lopenacc                     ! create variable on GPU
    ! local variables
    TYPE(t_list_element), POINTER :: element
    INTEGER                       :: idims(5), ndims

    ndims = 4
    idims(1:ndims) = ldims(1:ndims)
    CALL add_var_list_element_5d(ndims, BOOL_T, this_list, name,        &
      &   hgrid, vgrid, cf, grib2, idims, element, loutput, lcontainer, &
      &   lrestart, lrestart_cont, isteptype, lmiss, tlev_source, info, &
      &   vert_interp, hor_interp, in_group, verbose,                   &
      &   l_pp_scheduler_task, post_op, action_list, p5_l=p5,           &
      &   initval_l=initval, resetval_l=resetval, missval_l=missval,    &
      &   var_class=var_class, lopenacc=lopenacc)
    ptr => element%field%l_ptr(:,:,:,:,1)
    IF (PRESENT(new_element))  new_element => element
  END SUBROUTINE add_var_list_element_l4d


  !------------------------------------------------------------------------------------------------
  ! create (allocate) a new table entry
  ! optionally obtain pointer to 3d-field
  ! optionally overwrite default meta data
  !
  SUBROUTINE add_var_list_element_l3d(this_list, name, ptr,       &
    &   hgrid, vgrid, cf, grib2, ldims, loutput, lcontainer,      &
    &   lrestart, lrestart_cont, initval, isteptype,              &
    &   resetval, lmiss, missval, tlev_source, info, p5,          &
    &   vert_interp, hor_interp, in_group, verbose, new_element,  &
    &   l_pp_scheduler_task, post_op, action_list, var_class,     &
    &   lopenacc)
    !
    TYPE(t_var_list),        INTENT(inout)        :: this_list                    ! list
    CHARACTER(len=*),        INTENT(in)           :: name                         ! name of variable
    LOGICAL,        POINTER, INTENT(OUT)          :: ptr(:,:,:)                   ! reference to field
    INTEGER,                 INTENT(in)           :: hgrid                        ! horizontal grid type used
    INTEGER,                 INTENT(in)           :: vgrid                        ! vertical grid type used
    TYPE(t_cf_var),          INTENT(in)           :: cf                           ! CF related metadata
    TYPE(t_grib2_var),       INTENT(in)           :: grib2                        ! GRIB2 related metadata
    INTEGER,                 INTENT(in)           :: ldims(:)                     ! local dimensions
    LOGICAL,                 INTENT(in), OPTIONAL :: loutput                      ! output flag
    LOGICAL,                 INTENT(in), OPTIONAL :: lcontainer                   ! container flag
    LOGICAL,                 INTENT(in), OPTIONAL :: lrestart                     ! restart flag
    LOGICAL,                 INTENT(in), OPTIONAL :: lrestart_cont                ! continue restart if var not available
    LOGICAL,                 INTENT(in), OPTIONAL :: initval                      ! value if var not available
    INTEGER,                 INTENT(in), OPTIONAL :: isteptype                    ! type of statistical processing
    LOGICAL,                 INTENT(in), OPTIONAL :: resetval                     ! reset value (after accumulation)
    LOGICAL,                 INTENT(in), OPTIONAL :: lmiss                        ! missing value flag
    LOGICAL,                 INTENT(in), OPTIONAL :: missval                      ! missing value
    INTEGER,                 INTENT(in), OPTIONAL :: tlev_source                  ! actual TL for TL dependent vars
    TYPE(t_var_metadata),    POINTER,    OPTIONAL :: info                         ! returns reference to metadata
    LOGICAL,          CONTIGUOUS_TARGET, OPTIONAL :: p5(:,:,:,:,:)                ! provided pointer
    TYPE(t_vert_interp_meta),INTENT(in), OPTIONAL :: vert_interp                  ! vertical interpolation metadata
    TYPE(t_hor_interp_meta), INTENT(in), OPTIONAL :: hor_interp                   ! horizontal interpolation metadata
    LOGICAL,                 INTENT(in), OPTIONAL :: in_group(:)                  ! groups to which a variable belongs
    LOGICAL,                 INTENT(in), OPTIONAL :: verbose                      ! print information
    TYPE(t_list_element),    POINTER,    OPTIONAL :: new_element                  ! pointer to new var list element
    INTEGER,                 INTENT(in), OPTIONAL :: l_pp_scheduler_task          ! .TRUE., if field is updated by pp scheduler
    TYPE(t_post_op_meta),    INTENT(IN), OPTIONAL :: post_op                      ! "post-op" (small arithmetic operations) for this variable
    TYPE(t_var_action),      INTENT(IN), OPTIONAL :: action_list                  ! regularly triggered events
    INTEGER,                 INTENT(in), OPTIONAL :: var_class                    !< variable type/species
    LOGICAL,                 INTENT(in), OPTIONAL :: lopenacc                     ! create variable on GPU
    ! local variables
    TYPE(t_list_element), POINTER :: element
    INTEGER                       :: idims(5), ndims

    ndims = 3
    idims(1:ndims) = ldims(1:ndims)
    CALL add_var_list_element_5d(ndims, BOOL_T, this_list, name,        &
      &   hgrid, vgrid, cf, grib2, idims, element, loutput, lcontainer, &
      &   lrestart, lrestart_cont, isteptype, lmiss, tlev_source, info, &
      &   vert_interp, hor_interp, in_group, verbose,                   &
      &   l_pp_scheduler_task, post_op, action_list, p5_l=p5,           &
      &   initval_l=initval, resetval_l=resetval, missval_l=missval,    &
      &   var_class=var_class, lopenacc=lopenacc)
    ptr => element%field%l_ptr(:,:,:,1,1)
    IF (PRESENT(new_element))  new_element => element
  END SUBROUTINE add_var_list_element_l3d


  !------------------------------------------------------------------------------------------------
  ! create (allocate) a new table entry
  ! optionally obtain pointer to 2d-field
  ! optionally overwrite default meta data
  !
  SUBROUTINE add_var_list_element_l2d(this_list, name, ptr,       &
    &   hgrid, vgrid, cf, grib2, ldims, loutput, lcontainer,      &
    &   lrestart, lrestart_cont, initval, isteptype,              &
    &   resetval, lmiss, missval, tlev_source, info, p5,          &
    &   vert_interp, hor_interp, in_group, verbose, new_element,  &
    &   l_pp_scheduler_task, post_op, action_list, var_class,     &
    &   lopenacc)
    !
    TYPE(t_var_list),        INTENT(inout)        :: this_list                    ! list
    CHARACTER(len=*),        INTENT(in)           :: name                         ! name of variable
    LOGICAL,        POINTER, INTENT(OUT)          :: ptr(:,:)                     ! reference to field
    INTEGER,                 INTENT(in)           :: hgrid                        ! horizontal grid type used
    INTEGER,                 INTENT(in)           :: vgrid                        ! vertical grid type used
    TYPE(t_cf_var),          INTENT(in)           :: cf                           ! CF related metadata
    TYPE(t_grib2_var),       INTENT(in)           :: grib2                        ! GRIB2 related metadata
    INTEGER,                 INTENT(in)           :: ldims(:)                     ! local dimensions
    LOGICAL,                 INTENT(in), OPTIONAL :: loutput                      ! output flag
    LOGICAL,                 INTENT(in), OPTIONAL :: lcontainer                   ! container flag
    LOGICAL,                 INTENT(in), OPTIONAL :: lrestart                     ! restart flag
    LOGICAL,                 INTENT(in), OPTIONAL :: lrestart_cont                ! continue restart if var not available
    LOGICAL,                 INTENT(in), OPTIONAL :: initval                      ! value if var not available
    INTEGER,                 INTENT(in), OPTIONAL :: isteptype                    ! type of statistical processing
    LOGICAL,                 INTENT(in), OPTIONAL :: resetval                     ! reset value (after accumulation)
    LOGICAL,                 INTENT(in), OPTIONAL :: lmiss                        ! missing value flag
    LOGICAL,                 INTENT(in), OPTIONAL :: missval                      ! missing value
    INTEGER,                 INTENT(in), OPTIONAL :: tlev_source                  ! actual TL for TL dependent vars
    TYPE(t_var_metadata),    POINTER,    OPTIONAL :: info                         ! returns reference to metadata
    LOGICAL,          CONTIGUOUS_TARGET, OPTIONAL :: p5(:,:,:,:,:)                ! provided pointer
    TYPE(t_vert_interp_meta),INTENT(in), OPTIONAL :: vert_interp                  ! vertical interpolation metadata
    TYPE(t_hor_interp_meta), INTENT(in), OPTIONAL :: hor_interp                   ! horizontal interpolation metadata
    LOGICAL,                 INTENT(in), OPTIONAL :: in_group(:)                  ! groups to which a variable belongs
    LOGICAL,                 INTENT(in), OPTIONAL :: verbose                      ! print information
    TYPE(t_list_element),    POINTER,    OPTIONAL :: new_element                  ! pointer to new var list element
    INTEGER,                 INTENT(in), OPTIONAL :: l_pp_scheduler_task          ! .TRUE., if field is updated by pp scheduler
    TYPE(t_post_op_meta),    INTENT(IN), OPTIONAL :: post_op                      ! "post-op" (small arithmetic operations) for this variable
    TYPE(t_var_action),      INTENT(IN), OPTIONAL :: action_list                  ! regularly triggered events
    INTEGER,                 INTENT(in), OPTIONAL :: var_class                    !< variable type/species
    LOGICAL,                 INTENT(in), OPTIONAL :: lopenacc                     ! create variable on GPU
    ! local variables
    TYPE(t_list_element), POINTER :: element
    INTEGER                       :: idims(5), ndims

    ndims = 2
    idims(1:ndims) = ldims(1:ndims)
    CALL add_var_list_element_5d(ndims, BOOL_T, this_list, name,        &
      &   hgrid, vgrid, cf, grib2, idims, element, loutput, lcontainer, &
      &   lrestart, lrestart_cont, isteptype, lmiss, tlev_source, info, &
      &   vert_interp, hor_interp, in_group, verbose,                   &
      &   l_pp_scheduler_task, post_op, action_list, p5_l=p5,           &
      &   initval_l=initval, resetval_l=resetval, missval_l=missval,    &
      &   var_class=var_class, lopenacc=lopenacc)
    ptr => element%field%l_ptr(:,:,1,1,1)
    IF (PRESENT(new_element))  new_element => element
  END SUBROUTINE add_var_list_element_l2d


  !------------------------------------------------------------------------------------------------
  ! create (allocate) a new table entry
  ! optionally obtain pointer to 1d-field
  ! optionally overwrite default meta data
  !
  SUBROUTINE add_var_list_element_l1d(this_list, name, ptr,       &
    &   hgrid, vgrid, cf, grib2, ldims, loutput, lcontainer,      &
    &   lrestart, lrestart_cont, initval, isteptype,              &
    &   resetval, lmiss, missval, tlev_source, info, p5,          &
    &   vert_interp, hor_interp, in_group, verbose, new_element,  &
    &   l_pp_scheduler_task, post_op, action_list, var_class,     &
    &   lopenacc)
    !
    TYPE(t_var_list),        INTENT(inout)        :: this_list                    ! list
    CHARACTER(len=*),        INTENT(in)           :: name                         ! name of variable
    LOGICAL,        POINTER, INTENT(OUT)          :: ptr(:)                       ! reference to field
    INTEGER,                 INTENT(in)           :: hgrid                        ! horizontal grid type used
    INTEGER,                 INTENT(in)           :: vgrid                        ! vertical grid type used
    TYPE(t_cf_var),          INTENT(in)           :: cf                           ! CF related metadata
    TYPE(t_grib2_var),       INTENT(in)           :: grib2                        ! GRIB2 related metadata
    INTEGER,                 INTENT(in)           :: ldims(:)                     ! local dimensions
    LOGICAL,                 INTENT(in), OPTIONAL :: loutput                      ! output flag
    LOGICAL,                 INTENT(in), OPTIONAL :: lcontainer                   ! container flag
    LOGICAL,                 INTENT(in), OPTIONAL :: lrestart                     ! restart flag
    LOGICAL,                 INTENT(in), OPTIONAL :: lrestart_cont                ! continue restart if var not available
    LOGICAL,                 INTENT(in), OPTIONAL :: initval                      ! value if var not available
    INTEGER,                 INTENT(in), OPTIONAL :: isteptype                    ! type of statistical processing
    LOGICAL,                 INTENT(in), OPTIONAL :: resetval                     ! reset value (after accumulation)
    LOGICAL,                 INTENT(in), OPTIONAL :: lmiss                        ! missing value flag
    LOGICAL,                 INTENT(in), OPTIONAL :: missval                      ! missing value
    INTEGER,                 INTENT(in), OPTIONAL :: tlev_source                  ! actual TL for TL dependent vars
    TYPE(t_var_metadata),    POINTER,    OPTIONAL :: info                         ! returns reference to metadata
    LOGICAL,          CONTIGUOUS_TARGET, OPTIONAL :: p5(:,:,:,:,:)                ! provided pointer
    TYPE(t_vert_interp_meta),INTENT(in), OPTIONAL :: vert_interp                  ! vertical interpolation metadata
    TYPE(t_hor_interp_meta), INTENT(in), OPTIONAL :: hor_interp                   ! horizontal interpolation metadata
    LOGICAL,                 INTENT(in), OPTIONAL :: in_group(:)                  ! groups to which a variable belongs
    LOGICAL,                 INTENT(in), OPTIONAL :: verbose                      ! print information
    TYPE(t_list_element),    POINTER,    OPTIONAL :: new_element                  ! pointer to new var list element
    INTEGER,                 INTENT(in), OPTIONAL :: l_pp_scheduler_task          ! .TRUE., if field is updated by pp scheduler
    TYPE(t_post_op_meta),    INTENT(IN), OPTIONAL :: post_op                      ! "post-op" (small arithmetic operations) for this variable
    TYPE(t_var_action),      INTENT(IN), OPTIONAL :: action_list                  ! regularly triggered events
    INTEGER,                 INTENT(in), OPTIONAL :: var_class                    !< variable type/species
    LOGICAL,                 INTENT(in), OPTIONAL :: lopenacc                     ! create variable on GPU
    ! local variables
    TYPE(t_list_element), POINTER :: element
    INTEGER                       :: idims(5), ndims

    ndims = 1
    idims(1:ndims) = ldims(1:ndims)
    CALL add_var_list_element_5d(ndims, BOOL_T, this_list, name,        &
      &   hgrid, vgrid, cf, grib2, idims, element, loutput, lcontainer, &
      &   lrestart, lrestart_cont, isteptype, lmiss, tlev_source, info, &
      &   vert_interp, hor_interp, in_group, verbose,                   &
      &   l_pp_scheduler_task, post_op, action_list, p5_l=p5,           &
      &   initval_l=initval, resetval_l=resetval, missval_l=missval,    &
      &   var_class=var_class, lopenacc=lopenacc)
    ptr => element%field%l_ptr(:,1,1,1,1)
    IF (PRESENT(new_element))  new_element => element
  END SUBROUTINE add_var_list_element_l1d


  !================================================================================================
  !------------------------------------------------------------------------------------------------
  !
  ! Get element of a list
  !
  ! Specific routines for pointers of different rank
  !
  !================================================================================================
  ! REAL SECTION ----------------------------------------------------------------------------------
  !
  ! obtain pointer to 4d-field
  !
  SUBROUTINE get_var_list_element_r5d (this_list, name, ptr)
    TYPE(t_var_list), INTENT(in) :: this_list      ! list
    CHARACTER(len=*), INTENT(in) :: name           ! name of variable
    REAL(dp),         POINTER    :: ptr(:,:,:,:,:) ! reference to allocated field
    !
    TYPE(t_list_element), POINTER :: element
    !
    element => find_list_element (this_list, name)
    NULLIFY (ptr)
    IF (ASSOCIATED (element)) ptr => element%field%r_ptr(:,:,:,:,:)
    !
  END SUBROUTINE get_var_list_element_r5d
  !------------------------------------------------------------------------------------------------
  !
  ! obtain pointer to 4d-field
  !
  SUBROUTINE get_var_list_element_r4d (this_list, name, ptr)
    TYPE(t_var_list), INTENT(in) :: this_list    ! list
    CHARACTER(len=*), INTENT(in) :: name         ! name of variable
    REAL(dp),         POINTER    :: ptr(:,:,:,:) ! reference to allocated field
    !
    TYPE(t_list_element), POINTER :: element
    !
    element => find_list_element (this_list, name)
    NULLIFY (ptr)
    IF (ASSOCIATED (element)) ptr => element%field%r_ptr(:,:,:,:,1)
    !
  END SUBROUTINE get_var_list_element_r4d
  !------------------------------------------------------------------------------------------------
  !
  ! obtain pointer to 3d-field
  !
  SUBROUTINE get_var_list_element_r3d (this_list, name, ptr)
    TYPE(t_var_list), INTENT(in) :: this_list    ! list
    CHARACTER(len=*), INTENT(in) :: name         ! name of variable
    REAL(dp),         POINTER    :: ptr(:,:,:)   ! reference to allocated field
    !
    TYPE(t_list_element), POINTER :: element
    element => find_list_element (this_list, name)
    NULLIFY (ptr)
    IF (ASSOCIATED (element)) ptr => element%field%r_ptr(:,:,:,1,1)
    !
  END SUBROUTINE get_var_list_element_r3d
  !------------------------------------------------------------------------------------------------
  !
  ! obtain pointer to 2d-field
  !
  SUBROUTINE get_var_list_element_r2d (this_list, name, ptr)
    TYPE(t_var_list), INTENT(in) :: this_list  ! list
    CHARACTER(len=*), INTENT(in) :: name       ! name of variable
    REAL(dp),         POINTER    :: ptr(:,:)   ! reference to allocated field
    !
    TYPE(t_list_element), POINTER :: element
    !
    element => find_list_element (this_list, name)
    NULLIFY (ptr)
    IF (ASSOCIATED (element)) ptr => element%field%r_ptr(:,:,1,1,1)
    !
  END SUBROUTINE get_var_list_element_r2d
  !------------------------------------------------------------------------------------------------
  !
  ! obtain pointer to 1d-field
  !
  SUBROUTINE get_var_list_element_r1d (this_list, name, ptr)
    TYPE(t_var_list), INTENT(in) :: this_list ! list
    CHARACTER(len=*), INTENT(in) :: name      ! name of variable
    REAL(dp),         POINTER    :: ptr(:)    ! reference to allocated field
    !
    TYPE(t_list_element), POINTER :: element
    element => find_list_element (this_list, name)
    NULLIFY (ptr)
    IF (ASSOCIATED (element)) ptr => element%field%r_ptr(:,1,1,1,1)
    !
  END SUBROUTINE get_var_list_element_r1d


  ! Obtain pointer to 2D REAL field
  !
  FUNCTION fget_var_list_element_r1d (this_list, name) RESULT(ptr)
    TYPE(t_var_list), INTENT(in) :: this_list   ! list
    CHARACTER(len=*), INTENT(in) :: name        ! name of variable
    REAL(dp),         POINTER    :: ptr(:)      ! reference to allocated field
    !
    TYPE(t_list_element), POINTER :: element
    element => find_list_element (this_list, name)
    NULLIFY (ptr)
    IF (element%field%info%lcontained) THEN
      IF (ASSOCIATED (element)) ptr => element%field%r_ptr(:,element%field%info%ncontained,1,1,1)
    ELSE
      IF (ASSOCIATED (element)) ptr => element%field%r_ptr(:,1,1,1,1)
    ENDIF
    !
  END FUNCTION fget_var_list_element_r1d


  ! Obtain pointer to 2D REAL field
  !
  FUNCTION fget_var_list_element_r2d (this_list, name) RESULT(ptr)
    TYPE(t_var_list), INTENT(in) :: this_list   ! list
    CHARACTER(len=*), INTENT(in) :: name        ! name of variable
    REAL(dp),         POINTER    :: ptr(:,:)    ! reference to allocated field
    !
    TYPE(t_list_element), POINTER :: element
    element => find_list_element (this_list, name)
    NULLIFY (ptr)
    IF (element%field%info%lcontained) THEN
      IF (ASSOCIATED (element)) ptr => element%field%r_ptr(:,:,element%field%info%ncontained,1,1)
    ELSE
      IF (ASSOCIATED (element)) ptr => element%field%r_ptr(:,:,1,1,1)
    ENDIF
    !
  END FUNCTION fget_var_list_element_r2d


  ! Obtain pointer to 3D REAL field
  !
  FUNCTION fget_var_list_element_r3d (this_list, name) RESULT(ptr)
    TYPE(t_var_list), INTENT(in) :: this_list    ! list
    CHARACTER(len=*), INTENT(in) :: name         ! name of variable
    REAL(dp),         POINTER    :: ptr(:,:,:)   ! reference to allocated field
    !
    TYPE(t_list_element), POINTER :: element
    element => find_list_element (this_list, name)
    NULLIFY (ptr)
    IF (element%field%info%lcontained) THEN
      IF (ASSOCIATED (element)) ptr => element%field%r_ptr(:,:,:,element%field%info%ncontained,1)
    ELSE
      IF (ASSOCIATED (element)) ptr => element%field%r_ptr(:,:,:,1,1)
    ENDIF
    !
  END FUNCTION fget_var_list_element_r3d

  !================================================================================================
  ! REAL(sp) SECTION ----------------------------------------------------------------------------------
  !
  ! obtain pointer to 4d-field
  !
  SUBROUTINE get_var_list_element_s5d (this_list, name, ptr)
    TYPE(t_var_list), INTENT(in) :: this_list      ! list
    CHARACTER(len=*), INTENT(in) :: name           ! name of variable
    REAL(sp),         POINTER    :: ptr(:,:,:,:,:) ! reference to allocated field
    !
    TYPE(t_list_element), POINTER :: element
    !
    element => find_list_element (this_list, name)
    NULLIFY (ptr)
    IF (ASSOCIATED (element)) ptr => element%field%s_ptr(:,:,:,:,:)
    !
  END SUBROUTINE get_var_list_element_s5d
  !------------------------------------------------------------------------------------------------
  !
  ! obtain pointer to 4d-field
  !
  SUBROUTINE get_var_list_element_s4d (this_list, name, ptr)
    TYPE(t_var_list), INTENT(in) :: this_list    ! list
    CHARACTER(len=*), INTENT(in) :: name         ! name of variable
    REAL(sp),         POINTER    :: ptr(:,:,:,:) ! reference to allocated field
    !
    TYPE(t_list_element), POINTER :: element
    !
    element => find_list_element (this_list, name)
    NULLIFY (ptr)
    IF (ASSOCIATED (element)) ptr => element%field%s_ptr(:,:,:,:,1)
    !
  END SUBROUTINE get_var_list_element_s4d
  !------------------------------------------------------------------------------------------------
  !
  ! obtain pointer to 3d-field
  !
  SUBROUTINE get_var_list_element_s3d (this_list, name, ptr)
    TYPE(t_var_list), INTENT(in) :: this_list    ! list
    CHARACTER(len=*), INTENT(in) :: name         ! name of variable
    REAL(sp),         POINTER    :: ptr(:,:,:)   ! reference to allocated field
    !
    TYPE(t_list_element), POINTER :: element
    element => find_list_element (this_list, name)
    NULLIFY (ptr)
    IF (ASSOCIATED (element)) ptr => element%field%s_ptr(:,:,:,1,1)
    !
  END SUBROUTINE get_var_list_element_s3d
  !------------------------------------------------------------------------------------------------
  !
  ! obtain pointer to 2d-field
  !
  SUBROUTINE get_var_list_element_s2d (this_list, name, ptr)
    TYPE(t_var_list), INTENT(in) :: this_list  ! list
    CHARACTER(len=*), INTENT(in) :: name       ! name of variable
    REAL(sp),         POINTER    :: ptr(:,:)   ! reference to allocated field
    !
    TYPE(t_list_element), POINTER :: element
    !
    element => find_list_element (this_list, name)
    NULLIFY (ptr)
    IF (ASSOCIATED (element)) ptr => element%field%s_ptr(:,:,1,1,1)
    !
  END SUBROUTINE get_var_list_element_s2d
  !------------------------------------------------------------------------------------------------
  !
  ! obtain pointer to 1d-field
  !
  SUBROUTINE get_var_list_element_s1d (this_list, name, ptr)
    TYPE(t_var_list), INTENT(in) :: this_list ! list
    CHARACTER(len=*), INTENT(in) :: name      ! name of variable
    REAL(sp),         POINTER    :: ptr(:)    ! reference to allocated field
    !
    TYPE(t_list_element), POINTER :: element
    element => find_list_element (this_list, name)
    NULLIFY (ptr)
    IF (ASSOCIATED (element)) ptr => element%field%s_ptr(:,1,1,1,1)
    !
  END SUBROUTINE get_var_list_element_s1d


  ! Obtain pointer to 2D REAL field
  !
  FUNCTION fget_var_list_element_s1d (this_list, name) RESULT(ptr)
    TYPE(t_var_list), INTENT(in) :: this_list   ! list
    CHARACTER(len=*), INTENT(in) :: name        ! name of variable
    REAL(sp),         POINTER    :: ptr(:)      ! reference to allocated field
    !
    TYPE(t_list_element), POINTER :: element
    element => find_list_element (this_list, name)
    NULLIFY (ptr)
    IF (element%field%info%lcontained) THEN
      IF (ASSOCIATED (element)) ptr => element%field%s_ptr(:,element%field%info%ncontained,1,1,1)
    ELSE
      IF (ASSOCIATED (element)) ptr => element%field%s_ptr(:,1,1,1,1)
    ENDIF
    !
  END FUNCTION fget_var_list_element_s1d


  ! Obtain pointer to 2D REAL field
  !
  FUNCTION fget_var_list_element_s2d (this_list, name) RESULT(ptr)
    TYPE(t_var_list), INTENT(in) :: this_list   ! list
    CHARACTER(len=*), INTENT(in) :: name        ! name of variable
    REAL(sp),         POINTER    :: ptr(:,:)    ! reference to allocated field
    !
    TYPE(t_list_element), POINTER :: element
    element => find_list_element (this_list, name)
    NULLIFY (ptr)
    IF (element%field%info%lcontained) THEN
      IF (ASSOCIATED (element)) ptr => element%field%s_ptr(:,:,element%field%info%ncontained,1,1)
    ELSE
      IF (ASSOCIATED (element)) ptr => element%field%s_ptr(:,:,1,1,1)
    ENDIF
    !
  END FUNCTION fget_var_list_element_s2d


  ! Obtain pointer to 3D REAL field
  !
  FUNCTION fget_var_list_element_s3d (this_list, name) RESULT(ptr)
    TYPE(t_var_list), INTENT(in) :: this_list    ! list
    CHARACTER(len=*), INTENT(in) :: name         ! name of variable
    REAL(sp),         POINTER    :: ptr(:,:,:)   ! reference to allocated field
    !
    TYPE(t_list_element), POINTER :: element
    element => find_list_element (this_list, name)
    NULLIFY (ptr)
    IF (element%field%info%lcontained) THEN
      IF (ASSOCIATED (element)) ptr => element%field%s_ptr(:,:,:,element%field%info%ncontained,1)
    ELSE
      IF (ASSOCIATED (element)) ptr => element%field%s_ptr(:,:,:,1,1)
    ENDIF
    !
  END FUNCTION fget_var_list_element_s3d


  !================================================================================================
  ! INTEGER SECTION -------------------------------------------------------------------------------
  !
  ! obtain pointer to 5d-field
  !
  SUBROUTINE get_var_list_element_i5d (this_list, name, ptr)
    TYPE(t_var_list), INTENT(in) :: this_list      ! list
    CHARACTER(len=*), INTENT(in) :: name           ! name of variable
    INTEGER,          POINTER    :: ptr(:,:,:,:,:) ! reference to allocated field
    !
    TYPE(t_list_element), POINTER :: element
    !
    element => find_list_element (this_list, name)
    NULLIFY (ptr)
    IF (ASSOCIATED (element)) ptr => element%field%i_ptr(:,:,:,:,:)
    !
  END SUBROUTINE get_var_list_element_i5d
  !------------------------------------------------------------------------------------------------
  !
  ! obtain pointer to 4d-field
  !
  SUBROUTINE get_var_list_element_i4d (this_list, name, ptr)
    TYPE(t_var_list), INTENT(in) :: this_list    ! list
    CHARACTER(len=*), INTENT(in) :: name         ! name of variable
    INTEGER,          POINTER    :: ptr(:,:,:,:) ! reference to allocated field
    !
    TYPE(t_list_element), POINTER :: element
    !
    element => find_list_element (this_list, name)
    NULLIFY (ptr)
    IF (ASSOCIATED (element)) ptr => element%field%i_ptr(:,:,:,:,1)
    !
  END SUBROUTINE get_var_list_element_i4d
  !------------------------------------------------------------------------------------------------
  !
  ! obtain pointer to 3d-field
  !
  SUBROUTINE get_var_list_element_i3d (this_list, name, ptr)
    TYPE(t_var_list), INTENT(in) :: this_list    ! list
    CHARACTER(len=*), INTENT(in) :: name         ! name of variable
    INTEGER,          POINTER    :: ptr(:,:,:)   ! reference to allocated field
    !
    TYPE(t_list_element), POINTER :: element
    element => find_list_element (this_list, name)
    NULLIFY (ptr)
    IF (ASSOCIATED (element)) ptr => element%field%i_ptr(:,:,:,1,1)
    !
  END SUBROUTINE get_var_list_element_i3d
  !------------------------------------------------------------------------------------------------
  !
  ! obtain pointer to 2d-field
  !
  SUBROUTINE get_var_list_element_i2d (this_list, name, ptr)
    TYPE(t_var_list), INTENT(in) :: this_list  ! list
    CHARACTER(len=*), INTENT(in) :: name       ! name of variable
    INTEGER,          POINTER    :: ptr(:,:)   ! reference to allocated field
    !
    TYPE(t_list_element), POINTER :: element
    !
    element => find_list_element (this_list, name)
    NULLIFY (ptr)
    IF (ASSOCIATED (element)) ptr => element%field%i_ptr(:,:,1,1,1)
    !
  END SUBROUTINE get_var_list_element_i2d
  !------------------------------------------------------------------------------------------------
  !
  ! obtain pointer to 1d-field
  !
  SUBROUTINE get_var_list_element_i1d (this_list, name, ptr)
    TYPE(t_var_list), INTENT(in) :: this_list ! list
    CHARACTER(len=*), INTENT(in) :: name      ! name of variable
    INTEGER,          POINTER    :: ptr(:)    ! reference to allocated field
    !
    TYPE(t_list_element), POINTER :: element
    element => find_list_element (this_list, name)
    NULLIFY (ptr)
    IF (ASSOCIATED (element)) ptr => element%field%i_ptr(:,1,1,1,1)
    !
  END SUBROUTINE get_var_list_element_i1d
  !================================================================================================
  ! LOGICAL SECTION -------------------------------------------------------------------------------
  !
  ! obtain pointer to 5d-field
  !
  SUBROUTINE get_var_list_element_l5d (this_list, name, ptr)
    TYPE(t_var_list), INTENT(in) :: this_list      ! list
    CHARACTER(len=*), INTENT(in) :: name           ! name of variable
    LOGICAL,          POINTER    :: ptr(:,:,:,:,:) ! reference to allocated field
    !
    TYPE(t_list_element), POINTER :: element
    !
    element => find_list_element (this_list, name)
    NULLIFY (ptr)
    IF (ASSOCIATED (element)) ptr => element%field%l_ptr(:,:,:,:,:)
    !
  END SUBROUTINE get_var_list_element_l5d
  !------------------------------------------------------------------------------------------------
  !
  ! obtain pointer to 4d-field
  !
  SUBROUTINE get_var_list_element_l4d (this_list, name, ptr)
    TYPE(t_var_list), INTENT(in) :: this_list    ! list
    CHARACTER(len=*), INTENT(in) :: name         ! name of variable
    LOGICAL,          POINTER    :: ptr(:,:,:,:) ! reference to allocated field
    !
    TYPE(t_list_element), POINTER :: element
    !
    element => find_list_element (this_list, name)
    NULLIFY (ptr)
    IF (ASSOCIATED (element)) ptr => element%field%l_ptr(:,:,:,:,1)
    !
  END SUBROUTINE get_var_list_element_l4d
  !------------------------------------------------------------------------------------------------
  !
  ! obtain pointer to 3d-field
  !
  SUBROUTINE get_var_list_element_l3d (this_list, name, ptr)
    TYPE(t_var_list), INTENT(in) :: this_list    ! list
    CHARACTER(len=*), INTENT(in) :: name         ! name of variable
    LOGICAL,          POINTER    :: ptr(:,:,:)   ! reference to allocated field
    !
    TYPE(t_list_element), POINTER :: element
    element => find_list_element (this_list, name)
    NULLIFY (ptr)
    IF (ASSOCIATED (element)) ptr => element%field%l_ptr(:,:,:,1,1)
    !
  END SUBROUTINE get_var_list_element_l3d
  !------------------------------------------------------------------------------------------------
  !
  ! obtain pointer to 2d-field
  !
  SUBROUTINE get_var_list_element_l2d (this_list, name, ptr)
    TYPE(t_var_list), INTENT(in) :: this_list  ! list
    CHARACTER(len=*), INTENT(in) :: name       ! name of variable
    LOGICAL,          POINTER    :: ptr(:,:)   ! reference to allocated field
    !
    TYPE(t_list_element), POINTER :: element
    !
    element => find_list_element (this_list, name)
    NULLIFY (ptr)
    IF (ASSOCIATED (element)) ptr => element%field%l_ptr(:,:,1,1,1)
    !
  END SUBROUTINE get_var_list_element_l2d
  !------------------------------------------------------------------------------------------------
  !
  ! obtain pointer to 1d-field
  !
  SUBROUTINE get_var_list_element_l1d (this_list, name, ptr)
    TYPE(t_var_list), INTENT(in) :: this_list ! list
    CHARACTER(len=*), INTENT(in) :: name      ! name of variable
    LOGICAL,          POINTER    :: ptr(:)    ! reference to allocated field
    !
    TYPE(t_list_element), POINTER :: element
    element => find_list_element (this_list, name)
    NULLIFY (ptr)
    IF (ASSOCIATED (element)) ptr => element%field%l_ptr(:,1,1,1,1)
    !
  END SUBROUTINE get_var_list_element_l1d
  !================================================================================================
  !------------------------------------------------------------------------------------------------
  !
  ! Create a refernce to a list entry
  !
  ! Specific routines for pointers of different rank
  !
  !================================================================================================
  ! REAL SECTION ----------------------------------------------------------------------------------
  !
  ! create (allocate) a new table entry
  ! reference to an existing pointer to 3d-field
  ! optionally overwrite some default meta data
  !
  SUBROUTINE add_var_list_reference_r3d (this_list, target_name, name, ptr,                      &
       &                                 hgrid, vgrid, cf, grib2, ref_idx, ldims, loutput,       &
       &                                 lrestart, lrestart_cont, initval, isteptype,            &
       &                                 resetval, lmiss, missval, tlev_source, tracer_info,     &
       &                                 info, vert_interp, hor_interp, in_group, verbose,       &
       &                                 new_element, l_pp_scheduler_task, post_op, action_list, &
       &                                 opt_var_ref_pos, var_class)
    !
    TYPE(t_var_list),        INTENT(inout)           :: this_list
    CHARACTER(len=*),        INTENT(in)              :: target_name
    CHARACTER(len=*),        INTENT(in)              :: name
    REAL(dp), POINTER                                :: ptr(:,:,:)
    INTEGER,                 INTENT(in)              :: hgrid                      ! horizontal grid type used
    INTEGER,                 INTENT(in)              :: vgrid                      ! vertical grid type used
    TYPE(t_cf_var),          INTENT(in)              :: cf                         ! CF related metadata
    TYPE(t_grib2_var),       INTENT(in)              :: grib2                      ! GRIB2 related metadata
    INTEGER,                 INTENT(in)              :: ref_idx                    ! idx of slice to be referenced
    INTEGER,                 INTENT(in)              :: ldims(3)                   ! local dimensions, for checking
    LOGICAL,                 INTENT(in),    OPTIONAL :: loutput                    ! output flag
    LOGICAL,                 INTENT(in),    OPTIONAL :: lrestart                   ! restart flag
    LOGICAL,                 INTENT(in),    OPTIONAL :: lrestart_cont              ! continue restart if var not available
    REAL(dp),                INTENT(in),    OPTIONAL :: initval                    ! value if var not available
    INTEGER,                 INTENT(in),    OPTIONAL :: isteptype                  ! type of statistical processing
    REAL(dp),                INTENT(in),    OPTIONAL :: resetval                   ! reset value (after accumulation)
    LOGICAL,                 INTENT(in),    OPTIONAL :: lmiss                      ! missing value flag
    REAL(dp),                INTENT(in),    OPTIONAL :: missval                    ! missing value
    INTEGER,                 INTENT(in),    OPTIONAL :: tlev_source                ! actual TL for TL dependent vars
    CLASS(t_tracer_meta),    INTENT(in),    OPTIONAL :: tracer_info                ! tracer meta data
    TYPE(t_var_metadata), POINTER,          OPTIONAL :: info                       ! returns reference to metadata
    TYPE(t_vert_interp_meta),INTENT(in),    OPTIONAL :: vert_interp                ! vertical interpolation metadata
    TYPE(t_hor_interp_meta), INTENT(in),    OPTIONAL :: hor_interp                 ! horizontal interpolation metadata
    LOGICAL,                 INTENT(in),    OPTIONAL :: in_group(:)                ! groups to which a variable belongs
    LOGICAL,                 INTENT(in),    OPTIONAL :: verbose
    TYPE(t_list_element), POINTER,          OPTIONAL :: new_element                ! pointer to new var list element
    INTEGER,                 INTENT(in),    OPTIONAL :: l_pp_scheduler_task        ! .TRUE., if field is updated by pp scheduler
    TYPE(t_post_op_meta),    INTENT(IN),    OPTIONAL :: post_op                    !< "post-op" (small arithmetic operations) for this variable
    TYPE(t_var_action),      INTENT(IN),    OPTIONAL :: action_list                !< regularly triggered events
    INTEGER,                 INTENT(IN),    OPTIONAL :: opt_var_ref_pos            !< (optional:) position of container index
    INTEGER,                 INTENT(in),    OPTIONAL :: var_class                  !< variable type/species
    ! local variables
    CHARACTER(*), PARAMETER :: routine = modname//"::add_var_list_reference_r3d"
    !
    TYPE(t_list_element), POINTER :: target_element
    TYPE(t_var_metadata), POINTER :: target_info, ref_info
    TYPE(t_list_element), POINTER :: new_list_element
    TYPE(t_union_vals)            :: missvalt, initvalt, resetvalt
    INTEGER                       :: ndims, var_ref_pos, dim_indices(5), index
    LOGICAL :: in_group_new(MAX_GROUPS)             ! groups to which a variable belongs
                                                    ! (for taking into account tile groups)
    !
    ndims = 3
    target_element => find_list_element (this_list, target_name)
    target_info    => target_element%field%info
    IF (.NOT. ASSOCIATED(target_element%field%r_ptr))  CALL finish(routine, TRIM(name)//' not created.')

    !
    ! The parameter "var_ref_pos" contains the dimension index which
    ! points to the reference slice. Usually, this is "ndims+1", such
    ! that 3D slices, e.g., are stored in a 4D array as (:,:,:,1),
    ! (:,:,:,2), (:,:,:,3), etc.
    IF (PRESENT(opt_var_ref_pos)) THEN
      var_ref_pos    = opt_var_ref_pos
      IF (.NOT. target_info%lcontainer) &
        &  CALL finish(routine, "Container index does not make sense: Target is not a container variable!")
      IF ((target_info%var_ref_pos /= var_ref_pos) .AND. (target_info%var_ref_pos /= -1)) THEN
        CALL finish(routine, "Container index does not match the previously set value!")
      END IF
      target_info%var_ref_pos = var_ref_pos
    ELSE
      var_ref_pos    = ndims + 1
    END IF
    SELECT CASE(var_ref_pos)
    CASE (1)
      dim_indices    = (/ 2, 3, 4, 0, 0 /)
    CASE (2)
      dim_indices    = (/ 1, 3, 4, 0, 0 /)
    CASE (3)
      dim_indices    = (/ 1, 2, 4, 0, 0 /)
    CASE (4)
      dim_indices    = (/ 1, 2, 3, 0, 0 /)
    CASE DEFAULT
      CALL finish(routine, "Internal error!")
    END SELECT

    IF (target_info%lcontainer) THEN
      ! Counting the number of existing references is deactivated, if the slice index
      ! to be referenced is given explicitly.
!!$      IF ( PRESENT(ref_idx) ) THEN
        target_info%ncontained = target_info%ncontained+1
        ! only check validity of given slice index
        IF ( (ref_idx > SIZE(target_element%field%r_ptr, var_ref_pos)) .OR. (ref_idx < 1)) THEN
          WRITE (message_text, *) &
            &  'Slice idx ', ref_idx, ' for ', TRIM(name), &
            &  ' out of allowable range [1,',SIZE(target_element%field%r_ptr, var_ref_pos),']'
          CALL finish(routine, message_text)
        ENDIF
!!$      ELSE
!!$        target_info%ncontained = target_info%ncontained+1
!!$        IF (SIZE(target_element%field%r_ptr, var_ref_pos) < target_info%ncontained) THEN
!!$          WRITE (message_text, *) &
!!$            &  TRIM(name), ' exceeds the number of predefined entries in container:', &
!!$            &  SIZE(target_element%field%r_ptr, var_ref_pos)
!!$          CALL finish(routine, message_text)
!!$        ENDIF
!!$      ENDIF
      IF ( ANY(ldims(1:ndims) /=  target_info%used_dimensions(dim_indices(1:ndims))) ) THEN
        CALL finish(routine, TRIM(name)//' dimensions requested and available differ.')
      ENDIF
    ENDIF
    !
    ! add list entry
    !
    CALL append_list_element (this_list, new_list_element)
    IF (PRESENT(new_element)) new_element=>new_list_element
    ref_info => new_list_element%field%info
    CALL default_var_list_metadata(ref_info, this_list)

    !
    ! init local fields
    !
    missvalt  = ref_info%missval
    initvalt  = ref_info%initval
    resetvalt = ref_info%resetval
    IF (PRESENT(missval))  missvalt%rval  = missval
    IF (PRESENT(initval))  initvalt%rval  = initval
    IF (PRESENT(resetval)) resetvalt%rval = resetval
    !
    CALL set_var_metadata (new_list_element%field%info,                      &
         name=name, hgrid=hgrid, vgrid=vgrid,                                &
         cf=cf, grib2=grib2, ldims=ldims, loutput=loutput,                   &
         lrestart=lrestart, lrestart_cont=lrestart_cont, initval=initvalt,   &
         isteptype=isteptype, resetval=resetvalt, lmiss=lmiss,               &
         missval=missvalt, tlev_source=tlev_source,                          &
         vert_interp=vert_interp, hor_interp=hor_interp,                     &
         in_group=in_group, verbose=verbose,                                 &
         l_pp_scheduler_task=l_pp_scheduler_task,                            &
         post_op=post_op, action_list=action_list, var_class=var_class,      &
         data_type=REAL_T )
    ! set dynamic metadata, i.e. polymorphic tracer metadata
    CALL set_var_metadata_dyn (new_list_element%field%info_dyn,              &
                               tracer_info=tracer_info)

    ref_info%ndims = ndims
    ref_info%used_dimensions(:)       = 0
    ref_info%used_dimensions(1:ndims) = target_element%field%info%used_dimensions(dim_indices(1:ndims))

    index = 1
    !
    IF (PRESENT(var_class)) THEN
      IF ( ANY((/CLASS_TILE, CLASS_TILE_LAND/) == var_class)) THEN
        ! automatically add tile to its variable specific tile-group
        CALL var_groups_dyn%add(group_name=target_name, in_group_new=in_group_new, opt_in_group=in_group)
        !
        ! update in_group metainfo
        new_list_element%field%info%in_group(:) = in_group_new(:)
      ENDIF
    END IF
    !
    IF (target_info%lcontainer) THEN
      ref_info%lcontained                   = .TRUE.
      ref_info%used_dimensions(ndims+1)     = 1
      ref_info%var_ref_pos                  = var_ref_pos
      !
      ref_info%maxcontained = SIZE(target_element%field%r_ptr,var_ref_pos)
      !
!!$      IF ( PRESENT(ref_idx) ) THEN
        ref_info%ncontained = ref_idx
!!$      ELSE
!!$        ref_info%ncontained = target_info%ncontained
!!$      ENDIF
      index = ref_info%ncontained
    ENDIF
    SELECT CASE(var_ref_pos)
    CASE(1)
      ptr => target_element%field%r_ptr(index,:,:,:,1)
    CASE(2)
      ptr => target_element%field%r_ptr(:,index,:,:,1)
    CASE(3)
      ptr => target_element%field%r_ptr(:,:,index,:,1)
    CASE(4)
      ptr => target_element%field%r_ptr(:,:,:,index,1)
    CASE default
      CALL finish(routine, "internal error!")
    END SELECT
    new_list_element%field%r_ptr => target_element%field%r_ptr
    !
    IF (.NOT. ASSOCIATED(new_list_element%field%r_ptr)) THEN
      WRITE (0,*) 'problem with association of ptr for '//TRIM(name)
    ENDIF
    IF(PRESENT(info)) info => new_list_element%field%info
    IF (PRESENT(lmiss)) THEN
      ptr = new_list_element%field%info%missval%rval
    ELSE
      ptr = 0.0_dp
    END IF
  END SUBROUTINE add_var_list_reference_r3d

  !------------------------------------------------------------------------------------------------
  !
  ! create (allocate) a new table entry
  ! reference to an existing pointer to 2d-field
  ! optionally overwrite some default meta data
  !
  SUBROUTINE add_var_list_reference_r2d (this_list, target_name, name, ptr,                      &
       &                                 hgrid, vgrid, cf, grib2, ref_idx, ldims, loutput,       &
       &                                 lrestart, lrestart_cont, initval, isteptype,            &
       &                                 resetval, lmiss, missval, tlev_source, tracer_info,     &
       &                                 info, vert_interp, hor_interp, in_group,                &
       &                                 verbose, new_element, l_pp_scheduler_task,              &
       &                                 post_op, action_list, opt_var_ref_pos, var_class,       &
       &                                 idx_tracer, idx_diag) 

    TYPE(t_var_list),        INTENT(inout)        :: this_list
    CHARACTER(len=*),        INTENT(in)           :: target_name
    CHARACTER(len=*),        INTENT(in)           :: name
    REAL(dp), POINTER                             :: ptr(:,:)
    INTEGER,                 INTENT(in)           :: hgrid                       ! horizontal grid type used
    INTEGER,                 INTENT(in)           :: vgrid                       ! vertical grid type used
    TYPE(t_cf_var),          INTENT(in)           :: cf                          ! CF related metadata
    TYPE(t_grib2_var),       INTENT(in)           :: grib2                       ! GRIB2 related metadata
    INTEGER,                 INTENT(in)           :: ref_idx                     ! idx of slice to be referenced
    INTEGER,                 INTENT(in)           :: ldims(2)                    ! local dimensions, for checking
    LOGICAL,                 INTENT(in), OPTIONAL :: loutput                     ! output flag
    LOGICAL,                 INTENT(in), OPTIONAL :: lrestart                    ! restart flag
    LOGICAL,                 INTENT(in), OPTIONAL :: lrestart_cont               ! continue restart if var not available
    REAL(dp),                INTENT(in), OPTIONAL :: initval                     ! value if var not available
    INTEGER,                 INTENT(in), OPTIONAL :: isteptype                   ! type of statistical processing
    REAL(dp),                INTENT(in), OPTIONAL :: resetval                    ! reset value (after accumulation)
    LOGICAL,                 INTENT(in), OPTIONAL :: lmiss                       ! missing value flag
    REAL(dp),                INTENT(in), OPTIONAL :: missval                     ! missing value
    INTEGER,                 INTENT(in), OPTIONAL :: tlev_source                 ! actual TL for TL dependent vars
    CLASS(t_tracer_meta),    INTENT(in), OPTIONAL :: tracer_info                 ! tracer meta data
    TYPE(t_var_metadata), POINTER,       OPTIONAL :: info                        ! returns reference to metadata
    TYPE(t_vert_interp_meta),INTENT(in), OPTIONAL :: vert_interp                 ! vertical interpolation metadata
    TYPE(t_hor_interp_meta), INTENT(in), OPTIONAL :: hor_interp                  ! horizontal interpolation metadata
    LOGICAL,                 INTENT(in), OPTIONAL :: in_group(:)                 ! groups to which a variable belongs
    LOGICAL,                 INTENT(in), OPTIONAL :: verbose
    TYPE(t_list_element), POINTER,       OPTIONAL :: new_element                 ! pointer to new var list element
    INTEGER,                 INTENT(in), OPTIONAL :: l_pp_scheduler_task         ! .TRUE., if field is updated by pp scheduler
    TYPE(t_post_op_meta),    INTENT(IN), OPTIONAL :: post_op                     !< "post-op" (small arithmetic operations) for this variable
    TYPE(t_var_action),      INTENT(IN), OPTIONAL :: action_list                 !< regularly triggered events
    INTEGER,                 INTENT(IN), OPTIONAL :: opt_var_ref_pos             !< (optional:) position of container index
    INTEGER,                 INTENT(in), OPTIONAL :: var_class                   !< variable type/species
    INTEGER,                 INTENT(IN), OPTIONAL :: idx_tracer                  !< index of tracer in tracer container 
    INTEGER,                 INTENT(IN), OPTIONAL :: idx_diag                    !< index of tracer in diagnostics container 

    ! local variables
    CHARACTER(*), PARAMETER :: routine = modname//"::add_var_list_reference_r2d"
    !
    TYPE(t_list_element), POINTER :: target_element
    TYPE(t_var_metadata), POINTER :: target_info, ref_info
    TYPE(t_list_element), POINTER :: new_list_element
    TYPE(t_union_vals)            :: missvalt, initvalt, resetvalt
    INTEGER                       :: ndims, var_ref_pos, dim_indices(5), index
    LOGICAL :: in_group_new(MAX_GROUPS)             ! groups to which a variable belongs
                                                    ! (for taking into account tile groups)
    !
    ndims = 2

    target_element => find_list_element (this_list, target_name)
    target_info => target_element%field%info
    IF (.NOT. ASSOCIATED(target_element%field%r_ptr))  CALL finish(routine, TRIM(name)//' not created.')
    !
    ! The parameter "var_ref_pos" contains the dimension index which
    ! points to the reference slice. Usually, this is "ndims+1", such
    ! that 3D slices, e.g., are stored in a 4D array as (:,:,:,1),
    ! (:,:,:,2), (:,:,:,3), etc.
    IF (PRESENT(opt_var_ref_pos)) THEN
      var_ref_pos    = opt_var_ref_pos
      IF (.NOT. target_info%lcontainer) &
        &  CALL finish(routine, "Container index does not make sense: Target is not a container variable!")
      IF ((target_info%var_ref_pos /= var_ref_pos) .AND. (target_info%var_ref_pos /= -1)) THEN
        CALL finish(routine, "Container index does not match the previously set value!")
      END IF
      target_info%var_ref_pos = var_ref_pos
    ELSE
      var_ref_pos    = ndims + 1
    END IF
    SELECT CASE(var_ref_pos)
    CASE (1)
      dim_indices    = (/ 2, 3, 0, 0, 0 /)
    CASE (2)
      dim_indices    = (/ 1, 3, 0, 0, 0 /)
    CASE (3)
      dim_indices    = (/ 1, 2, 0, 0, 0 /)
    CASE DEFAULT
      CALL finish(routine, "Internal error!")
    END SELECT

    IF (target_info%lcontainer) THEN
      ! Counting the number of existing references is deactivated, if the slice index
      ! to be referenced is given explicitly.
!!$      IF ( PRESENT(ref_idx) ) THEN
        target_info%ncontained = target_info%ncontained+1
        ! only check validity of given slice index
        IF ( (ref_idx > SIZE(target_element%field%r_ptr, var_ref_pos)) .OR. (ref_idx < 1)) THEN
          WRITE (message_text, *) &
            &  'Slice idx ', ref_idx, ' for ', TRIM(name), &
            &  ' out of allowable range [1,',SIZE(target_element%field%r_ptr, var_ref_pos),']'
          CALL finish(routine, message_text)
        ENDIF
!!$      ELSE
!!$        target_info%ncontained = target_info%ncontained+1
!!$        IF (SIZE(target_element%field%r_ptr, var_ref_pos) < target_info%ncontained) THEN
!!$          WRITE (message_text, *) &
!!$            &  TRIM(name), ' exceeds the number of predefined entries in container:', &
!!$            &  SIZE(target_element%field%r_ptr, var_ref_pos)
!!$          CALL finish(routine, message_text)
!!$        ENDIF
!!$      ENDIF
      IF (ANY(ldims(1:ndims) /=  target_info%used_dimensions(dim_indices(1:ndims)))) THEN
        CALL finish(routine, TRIM(name)//' dimensions requested and available differ.')
      ENDIF
    ENDIF
    !
    ! add list entry
    !
    CALL append_list_element (this_list, new_list_element)
    IF (PRESENT(new_element)) new_element=>new_list_element
    ref_info => new_list_element%field%info
    CALL default_var_list_metadata(ref_info, this_list)
    !
    ! init local fields
    !
    missvalt  = ref_info%missval
    initvalt  = ref_info%initval
    resetvalt = ref_info%resetval
    IF (PRESENT(missval))  missvalt%rval  = missval
    IF (PRESENT(initval))  initvalt%rval  = initval
    IF (PRESENT(resetval)) resetvalt%rval = resetval
    !
    CALL set_var_metadata (new_list_element%field%info,                      &
         name=name, hgrid=hgrid, vgrid=vgrid,                                &
         cf=cf, grib2=grib2, ldims=ldims, loutput=loutput,                   &
         lrestart=lrestart, lrestart_cont=lrestart_cont, initval=initvalt,   &
         isteptype=isteptype, resetval=resetvalt, lmiss=lmiss,               &
         missval=missvalt, tlev_source=tlev_source,                          &
         vert_interp=vert_interp, hor_interp=hor_interp,                     &
         in_group=in_group, verbose=verbose,                                 &
         l_pp_scheduler_task=l_pp_scheduler_task,                            &
         post_op=post_op, action_list=action_list, var_class=var_class,      &
         data_type=REAL_T, idx_tracer=idx_tracer, idx_diag=idx_diag)     
    ! set dynamic metadata, i.e. polymorphic tracer metadata
    CALL set_var_metadata_dyn (new_list_element%field%info_dyn,              &
                               tracer_info=tracer_info)

    ref_info%ndims = ndims
    ref_info%used_dimensions(:)       = 0
    ref_info%used_dimensions(1:ndims) = target_element%field%info%used_dimensions(dim_indices(1:ndims))

    index = 1
    !
    IF (PRESENT(var_class)) THEN
      IF ( ANY((/CLASS_TILE, CLASS_TILE_LAND/) == var_class)) THEN
        ! automatically add tile to its variable specific tile-group
        CALL var_groups_dyn%add(group_name=target_name, in_group_new=in_group_new, opt_in_group=in_group)
        !
        ! update in_group metainfo
        new_list_element%field%info%in_group(:) = in_group_new(:)
      ENDIF
    END IF

    IF (target_info%lcontainer) THEN
      ref_info%lcontained                   = .TRUE.
      ref_info%used_dimensions(ndims+1)     = 1
      ref_info%var_ref_pos                  = var_ref_pos
      !
      ref_info%maxcontained = SIZE(target_element%field%r_ptr,var_ref_pos)
      !
!!$      IF ( PRESENT(ref_idx) ) THEN
        ref_info%ncontained = ref_idx
!!$      ELSE
!!$        ref_info%ncontained = target_info%ncontained
!!$      ENDIF
      index = ref_info%ncontained
    ENDIF
    SELECT CASE(var_ref_pos)
    CASE(1)
      ptr => target_element%field%r_ptr(index,:,:,1,1)
    CASE(2)
      ptr => target_element%field%r_ptr(:,index,:,1,1)
    CASE(3)
      ptr => target_element%field%r_ptr(:,:,index,1,1)
    CASE default
      CALL finish(routine, "internal error!")
    END SELECT
    new_list_element%field%r_ptr => target_element%field%r_ptr
    !
    IF (.NOT. ASSOCIATED(new_list_element%field%r_ptr)) THEN
      WRITE (0,*) 'problem with association of ptr for '//TRIM(name)
    ENDIF
    !
    IF(PRESENT(info)) info => new_list_element%field%info
    IF (PRESENT(lmiss)) THEN
      ptr = new_list_element%field%info%missval%rval
    ELSE
      ptr = 0.0_dp
    END IF
  END SUBROUTINE add_var_list_reference_r2d

  SUBROUTINE add_var_list_reference_s3d (this_list, target_name, name, ptr,                      &
       &                                 hgrid, vgrid, cf, grib2, ref_idx, ldims, loutput,       &
       &                                 lrestart, lrestart_cont, initval, isteptype,            &
       &                                 resetval, lmiss, missval, tlev_source, tracer_info,     &
       &                                 info, vert_interp, hor_interp, in_group, verbose,       &
       &                                 new_element, l_pp_scheduler_task, post_op, action_list, &
       &                                 opt_var_ref_pos, var_class)
    !
    TYPE(t_var_list),        INTENT(inout)           :: this_list
    CHARACTER(len=*),        INTENT(in)              :: target_name
    CHARACTER(len=*),        INTENT(in)              :: name
    REAL(sp), POINTER                                :: ptr(:,:,:)
    INTEGER,                 INTENT(in)              :: hgrid                      ! horizontal grid type used
    INTEGER,                 INTENT(in)              :: vgrid                      ! vertical grid type used
    TYPE(t_cf_var),          INTENT(in)              :: cf                         ! CF related metadata
    TYPE(t_grib2_var),       INTENT(in)              :: grib2                      ! GRIB2 related metadata
    INTEGER,                 INTENT(in)              :: ref_idx                    ! idx of slice to be referenced
    INTEGER,                 INTENT(in)              :: ldims(3)                   ! local dimensions, for checking
    LOGICAL,                 INTENT(in),    OPTIONAL :: loutput                    ! output flag
    LOGICAL,                 INTENT(in),    OPTIONAL :: lrestart                   ! restart flag
    LOGICAL,                 INTENT(in),    OPTIONAL :: lrestart_cont              ! continue restart if var not available
    REAL(sp),                INTENT(in),    OPTIONAL :: initval                    ! value if var not available
    INTEGER,                 INTENT(in),    OPTIONAL :: isteptype                  ! type of statistical processing
    REAL(sp),                INTENT(in),    OPTIONAL :: resetval                   ! reset value (after accumulation)
    LOGICAL,                 INTENT(in),    OPTIONAL :: lmiss                      ! missing value flag
    REAL(sp),                INTENT(in),    OPTIONAL :: missval                    ! missing value
    INTEGER,                 INTENT(in),    OPTIONAL :: tlev_source                ! actual TL for TL dependent vars
    CLASS(t_tracer_meta),    INTENT(in),    OPTIONAL :: tracer_info                ! tracer meta data
    TYPE(t_var_metadata), POINTER,          OPTIONAL :: info                       ! returns reference to metadata
    TYPE(t_vert_interp_meta),INTENT(in),    OPTIONAL :: vert_interp                ! vertical interpolation metadata
    TYPE(t_hor_interp_meta), INTENT(in),    OPTIONAL :: hor_interp                 ! horizontal interpolation metadata
    LOGICAL,                 INTENT(in),    OPTIONAL :: in_group(:)                ! groups to which a variable belongs
    LOGICAL,                 INTENT(in),    OPTIONAL :: verbose
    TYPE(t_list_element), POINTER,          OPTIONAL :: new_element                ! pointer to new var list element
    INTEGER,                 INTENT(in),    OPTIONAL :: l_pp_scheduler_task        ! .TRUE., if field is updated by pp scheduler
    TYPE(t_post_op_meta),    INTENT(IN),    OPTIONAL :: post_op                    !< "post-op" (small arithmetic operations) for this variable
    TYPE(t_var_action),      INTENT(IN),    OPTIONAL :: action_list                !< regularly triggered events
    INTEGER,                 INTENT(IN),    OPTIONAL :: opt_var_ref_pos            !< (optional:) position of container index
    INTEGER,                 INTENT(in),    OPTIONAL :: var_class                  !< variable type/species
    ! local variables
    CHARACTER(*), PARAMETER :: routine = modname//"::add_var_list_reference_s3d"
    !
    TYPE(t_list_element), POINTER :: target_element
    TYPE(t_var_metadata), POINTER :: target_info, ref_info
    TYPE(t_list_element), POINTER :: new_list_element
    TYPE(t_union_vals)            :: missvalt, initvalt, resetvalt
    INTEGER                       :: ndims, var_ref_pos, dim_indices(5), index
    LOGICAL :: in_group_new(MAX_GROUPS)             ! groups to which a variable belongs
                                                    ! (for taking into account tile groups)
    !
    ndims = 3
    target_element => find_list_element (this_list, target_name)
    target_info    => target_element%field%info
    IF (.NOT. ASSOCIATED(target_element%field%s_ptr))  CALL finish(routine, TRIM(name)//' not created.')

    !
    ! The parameter "var_ref_pos" contains the dimension index which
    ! points to the reference slice. Usually, this is "ndims+1", such
    ! that 3D slices, e.g., are stored in a 4D array as (:,:,:,1),
    ! (:,:,:,2), (:,:,:,3), etc.
    IF (PRESENT(opt_var_ref_pos)) THEN
      var_ref_pos    = opt_var_ref_pos
      IF (.NOT. target_info%lcontainer) &
        &  CALL finish(routine, "Container index does not make sense: Target is not a container variable!")
      IF ((target_info%var_ref_pos /= var_ref_pos) .AND. (target_info%var_ref_pos /= -1)) THEN
        CALL finish(routine, "Container index does not match the previously set value!")
      END IF
      target_info%var_ref_pos = var_ref_pos
    ELSE
      var_ref_pos    = ndims + 1
    END IF
    SELECT CASE(var_ref_pos)
    CASE (1)
      dim_indices    = (/ 2, 3, 4, 0, 0 /)
    CASE (2)
      dim_indices    = (/ 1, 3, 4, 0, 0 /)
    CASE (3)
      dim_indices    = (/ 1, 2, 4, 0, 0 /)
    CASE (4)
      dim_indices    = (/ 1, 2, 3, 0, 0 /)
    CASE DEFAULT
      CALL finish(routine, "Internal error!")
    END SELECT

    IF (target_info%lcontainer) THEN
      ! Counting the number of existing references is deactivated, if the slice index
      ! to be referenced is given explicitly.
!!$      IF ( PRESENT(ref_idx) ) THEN
        target_info%ncontained = target_info%ncontained+1
        ! only check validity of given slice index
        IF ( (ref_idx > SIZE(target_element%field%s_ptr, var_ref_pos)) .OR. (ref_idx < 1)) THEN
          WRITE (message_text, *) &
            &  'Slice idx ', ref_idx, ' for ', TRIM(name), &
            &  ' out of allowable range [1,',SIZE(target_element%field%s_ptr, var_ref_pos),']'
          CALL finish(routine, message_text)
        ENDIF
!!$      ELSE
!!$        target_info%ncontained = target_info%ncontained+1
!!$        IF (SIZE(target_element%field%s_ptr, var_ref_pos) < target_info%ncontained) THEN
!!$          WRITE (message_text, *) &
!!$            &  TRIM(name), ' exceeds the number of predefined entries in container:', &
!!$            &  SIZE(target_element%field%s_ptr, var_ref_pos)
!!$          CALL finish(routine, message_text)
!!$        ENDIF
!!$      ENDIF
      IF ( ANY(ldims(1:ndims) /=  target_info%used_dimensions(dim_indices(1:ndims))) ) THEN
        CALL finish(routine, TRIM(name)//' dimensions requested and available differ.')
      ENDIF
    ENDIF
    !
    ! add list entry
    !
    CALL append_list_element (this_list, new_list_element)
    IF (PRESENT(new_element)) new_element=>new_list_element
    ref_info => new_list_element%field%info
    CALL default_var_list_metadata(ref_info, this_list)

    !
    ! init local fields
    !
    missvalt  = ref_info%missval
    initvalt  = ref_info%initval
    resetvalt = ref_info%resetval
    IF (PRESENT(missval))  missvalt%sval  = missval
    IF (PRESENT(initval))  initvalt%sval  = initval
    IF (PRESENT(resetval)) resetvalt%sval = resetval
    CALL set_var_metadata (new_list_element%field%info,                      &
         name=name, hgrid=hgrid, vgrid=vgrid,                                &
         cf=cf, grib2=grib2, ldims=ldims, loutput=loutput,                   &
         lrestart=lrestart, lrestart_cont=lrestart_cont, initval=initvalt,   &
         isteptype=isteptype, resetval=resetvalt, lmiss=lmiss,               &
         missval=missvalt, tlev_source=tlev_source,                          &
         vert_interp=vert_interp, hor_interp=hor_interp,                     &
         in_group=in_group, verbose=verbose,                                 &
         l_pp_scheduler_task=l_pp_scheduler_task,                            &
         post_op=post_op, action_list=action_list, var_class=var_class,      &
         data_type=SINGLE_T )
    ! set dynamic metadata, i.e. polymorphic tracer metadata
    CALL set_var_metadata_dyn (new_list_element%field%info_dyn,              &
                               tracer_info=tracer_info)

    ref_info%ndims = ndims
    ref_info%used_dimensions(:)       = 0
    ref_info%used_dimensions(1:ndims) = target_element%field%info%used_dimensions(dim_indices(1:ndims))

    index = 1
    !
    IF (PRESENT(var_class)) THEN
      IF ( ANY((/CLASS_TILE, CLASS_TILE_LAND/) == var_class)) THEN
        ! automatically add tile to its variable specific tile-group
        CALL var_groups_dyn%add(group_name=target_name, in_group_new=in_group_new, opt_in_group=in_group)
        !
        ! update in_group metainfo
        new_list_element%field%info%in_group(:) = in_group_new(:)
      ENDIF
    END IF
    !
    IF (target_info%lcontainer) THEN
      ref_info%lcontained                   = .TRUE.
      ref_info%used_dimensions(ndims+1)     = 1
      ref_info%var_ref_pos                  = var_ref_pos
      !
      ref_info%maxcontained = SIZE(target_element%field%s_ptr,var_ref_pos)
      !
!!$      IF ( PRESENT(ref_idx) ) THEN
        ref_info%ncontained = ref_idx
!!$      ELSE
!!$        ref_info%ncontained = target_info%ncontained
!!$      ENDIF
      index = ref_info%ncontained
    ENDIF
    SELECT CASE(var_ref_pos)
    CASE(1)
      ptr => target_element%field%s_ptr(index,:,:,:,1)
    CASE(2)
      ptr => target_element%field%s_ptr(:,index,:,:,1)
    CASE(3)
      ptr => target_element%field%s_ptr(:,:,index,:,1)
    CASE(4)
      ptr => target_element%field%s_ptr(:,:,:,index,1)
    CASE default
      CALL finish(routine, "internal error!")
    END SELECT
    new_list_element%field%s_ptr => target_element%field%s_ptr
    !
    IF (.NOT. ASSOCIATED(new_list_element%field%s_ptr)) THEN
      WRITE (0,*) 'problem with association of ptr for '//TRIM(name)
    ENDIF
    !
    IF(PRESENT(info)) info => new_list_element%field%info
    IF (PRESENT(lmiss)) THEN
      ptr = new_list_element%field%info%missval%sval
    ELSE
      ptr = 0.0_sp
    END IF
  END SUBROUTINE add_var_list_reference_s3d

  !------------------------------------------------------------------------------------------------
  !
  ! create (allocate) a new table entry
  ! reference to an existing pointer to 2d-field
  ! optionally overwrite some default meta data
  !
  SUBROUTINE add_var_list_reference_s2d (this_list, target_name, name, ptr,                      &
       &                                 hgrid, vgrid, cf, grib2, ref_idx, ldims, loutput,       &
       &                                 lrestart, lrestart_cont, initval, isteptype,            &
       &                                 resetval, lmiss, missval, tlev_source, tracer_info,     &
       &                                 info, vert_interp, hor_interp, in_group,                &
       &                                 verbose, new_element, l_pp_scheduler_task,              &
       &                                 post_op, action_list, opt_var_ref_pos, var_class)

    TYPE(t_var_list),        INTENT(inout)        :: this_list
    CHARACTER(len=*),        INTENT(in)           :: target_name
    CHARACTER(len=*),        INTENT(in)           :: name
    REAL(sp), POINTER                             :: ptr(:,:)
    INTEGER,                 INTENT(in)           :: hgrid                       ! horizontal grid type used
    INTEGER,                 INTENT(in)           :: vgrid                       ! vertical grid type used
    TYPE(t_cf_var),          INTENT(in)           :: cf                          ! CF related metadata
    TYPE(t_grib2_var),       INTENT(in)           :: grib2                       ! GRIB2 related metadata
    INTEGER,                 INTENT(in)           :: ref_idx                     ! idx of slice to be referenced
    INTEGER,                 INTENT(in)           :: ldims(2)                    ! local dimensions, for checking
    LOGICAL,                 INTENT(in), OPTIONAL :: loutput                     ! output flag
    LOGICAL,                 INTENT(in), OPTIONAL :: lrestart                    ! restart flag
    LOGICAL,                 INTENT(in), OPTIONAL :: lrestart_cont               ! continue restart if var not available
    REAL(sp),                INTENT(in), OPTIONAL :: initval                     ! value if var not available
    INTEGER,                 INTENT(in), OPTIONAL :: isteptype                   ! type of statistical processing
    REAL(sp),                INTENT(in), OPTIONAL :: resetval                    ! reset value (after accumulation)
    LOGICAL,                 INTENT(in), OPTIONAL :: lmiss                       ! missing value flag
    REAL(sp),                INTENT(in), OPTIONAL :: missval                     ! missing value
    INTEGER,                 INTENT(in), OPTIONAL :: tlev_source                 ! actual TL for TL dependent vars
    CLASS(t_tracer_meta),    INTENT(in), OPTIONAL :: tracer_info                 ! tracer meta data
    TYPE(t_var_metadata), POINTER,       OPTIONAL :: info                        ! returns reference to metadata
    TYPE(t_vert_interp_meta),INTENT(in), OPTIONAL :: vert_interp                 ! vertical interpolation metadata
    TYPE(t_hor_interp_meta), INTENT(in), OPTIONAL :: hor_interp                  ! horizontal interpolation metadata
    LOGICAL,                 INTENT(in), OPTIONAL :: in_group(:)                 ! groups to which a variable belongs
    LOGICAL,                 INTENT(in), OPTIONAL :: verbose
    TYPE(t_list_element), POINTER,       OPTIONAL :: new_element                 ! pointer to new var list element
    INTEGER,                 INTENT(in), OPTIONAL :: l_pp_scheduler_task         ! .TRUE., if field is updated by pp scheduler
    TYPE(t_post_op_meta),    INTENT(IN), OPTIONAL :: post_op                     !< "post-op" (small arithmetic operations) for this variable
    TYPE(t_var_action),      INTENT(IN), OPTIONAL :: action_list                 !< regularly triggered events
    INTEGER,                 INTENT(IN), OPTIONAL :: opt_var_ref_pos             !< (optional:) position of container index
    INTEGER,                 INTENT(in), OPTIONAL :: var_class                   !< variable type/species
    ! local variables
    CHARACTER(*), PARAMETER :: routine = modname//"::add_var_list_reference_s2d"
    !
    TYPE(t_list_element), POINTER :: target_element
    TYPE(t_var_metadata), POINTER :: target_info, ref_info
    TYPE(t_list_element), POINTER :: new_list_element
    TYPE(t_union_vals)            :: missvalt, initvalt, resetvalt
    INTEGER                       :: ndims, var_ref_pos, dim_indices(5), index
    LOGICAL :: in_group_new(MAX_GROUPS)             ! groups to which a variable belongs
                                                    ! (for taking into account tile groups)
    !
    ndims = 2

    target_element => find_list_element (this_list, target_name)
    target_info => target_element%field%info
    IF (.NOT. ASSOCIATED(target_element%field%s_ptr))  CALL finish(routine, TRIM(name)//' not created.')
    !
    ! The parameter "var_ref_pos" contains the dimension index which
    ! points to the reference slice. Usually, this is "ndims+1", such
    ! that 3D slices, e.g., are stored in a 4D array as (:,:,:,1),
    ! (:,:,:,2), (:,:,:,3), etc.
    IF (PRESENT(opt_var_ref_pos)) THEN
      var_ref_pos    = opt_var_ref_pos
      IF (.NOT. target_info%lcontainer) &
        &  CALL finish(routine, "Container index does not make sense: Target is not a container variable!")
      IF ((target_info%var_ref_pos /= var_ref_pos) .AND. (target_info%var_ref_pos /= -1)) THEN
        CALL finish(routine, "Container index does not match the previously set value!")
      END IF
      target_info%var_ref_pos = var_ref_pos
    ELSE
      var_ref_pos    = ndims + 1
    END IF
    SELECT CASE(var_ref_pos)
    CASE (1)
      dim_indices    = (/ 2, 3, 0, 0, 0 /)
    CASE (2)
      dim_indices    = (/ 1, 3, 0, 0, 0 /)
    CASE (3)
      dim_indices    = (/ 1, 2, 0, 0, 0 /)
    CASE DEFAULT
      CALL finish(routine, "Internal error!")
    END SELECT

    IF (target_info%lcontainer) THEN
      ! Counting the number of existing references is deactivated, if the slice index
      ! to be referenced is given explicitly.
!!$      IF ( PRESENT(ref_idx) ) THEN
        target_info%ncontained = target_info%ncontained+1
        ! only check validity of given slice index
        IF ( (ref_idx > SIZE(target_element%field%s_ptr, var_ref_pos)) .OR. (ref_idx < 1)) THEN
          WRITE (message_text, *) &
            &  'Slice idx ', ref_idx, ' for ', TRIM(name), &
            &  ' out of allowable range [1,',SIZE(target_element%field%s_ptr, var_ref_pos),']'
          CALL finish(routine, message_text)
        ENDIF
!!$      ELSE
!!$        target_info%ncontained = target_info%ncontained+1
!!$        IF (SIZE(target_element%field%s_ptr, var_ref_pos) < target_info%ncontained) THEN
!!$          WRITE (message_text, *) &
!!$            &  TRIM(name), ' exceeds the number of predefined entries in container:', &
!!$            &  SIZE(target_element%field%s_ptr, var_ref_pos)
!!$          CALL finish(routine, message_text)
!!$        ENDIF
!!$      ENDIF
      IF (ANY(ldims(1:ndims) /=  target_info%used_dimensions(dim_indices(1:ndims)))) THEN
        CALL finish(routine, TRIM(name)//' dimensions requested and available differ.')
      ENDIF
    ENDIF
    !
    ! add list entry
    !
    CALL append_list_element (this_list, new_list_element)
    IF (PRESENT(new_element)) new_element=>new_list_element
    ref_info => new_list_element%field%info
    CALL default_var_list_metadata(ref_info, this_list)
    !
    ! init local fields
    !
    missvalt  = ref_info%missval
    initvalt  = ref_info%initval
    resetvalt = ref_info%resetval
    IF (PRESENT(missval))  missvalt%sval  = missval
    IF (PRESENT(initval))  initvalt%sval  = initval
    IF (PRESENT(resetval)) resetvalt%sval = resetval
    CALL set_var_metadata (new_list_element%field%info,                      &
         name=name, hgrid=hgrid, vgrid=vgrid,                                &
         cf=cf, grib2=grib2, ldims=ldims, loutput=loutput,                   &
         lrestart=lrestart, lrestart_cont=lrestart_cont, initval=initvalt,   &
         isteptype=isteptype, resetval=resetvalt, lmiss=lmiss,               &
         missval=missvalt, tlev_source=tlev_source,                          &
         vert_interp=vert_interp, hor_interp=hor_interp,                     &
         in_group=in_group, verbose=verbose,                                 &
         l_pp_scheduler_task=l_pp_scheduler_task,                            &
         post_op=post_op, action_list=action_list, var_class=var_class,      &
         data_type=SINGLE_T )
    ! set dynamic metadata, i.e. polymorphic tracer metadata
    CALL set_var_metadata_dyn (new_list_element%field%info_dyn,              &
                               tracer_info=tracer_info)

    ref_info%ndims = ndims
    ref_info%used_dimensions(:)       = 0
    ref_info%used_dimensions(1:ndims) = target_element%field%info%used_dimensions(dim_indices(1:ndims))

    index = 1
    !
    IF (PRESENT(var_class)) THEN
      IF ( ANY((/CLASS_TILE, CLASS_TILE_LAND/) == var_class)) THEN
        ! automatically add tile to its variable specific tile-group
        CALL var_groups_dyn%add(group_name=target_name, in_group_new=in_group_new, opt_in_group=in_group)
        !
        ! update in_group metainfo
        new_list_element%field%info%in_group(:) = in_group_new(:)
      ENDIF
    END IF

    IF (target_info%lcontainer) THEN
      ref_info%lcontained                   = .TRUE.
      ref_info%used_dimensions(ndims+1)     = 1
      ref_info%var_ref_pos                  = var_ref_pos
      !
      ref_info%maxcontained = SIZE(target_element%field%s_ptr,var_ref_pos)
      !
!!$      IF ( PRESENT(ref_idx) ) THEN
        ref_info%ncontained = ref_idx
!!$      ELSE
!!$        ref_info%ncontained = target_info%ncontained
!!$      ENDIF
      index = ref_info%ncontained
    ENDIF
    SELECT CASE(var_ref_pos)
    CASE(1)
      ptr => target_element%field%s_ptr(index,:,:,1,1)
    CASE(2)
      ptr => target_element%field%s_ptr(:,index,:,1,1)
    CASE(3)
      ptr => target_element%field%s_ptr(:,:,index,1,1)
    CASE default
      CALL finish(routine, "internal error!")
    END SELECT
    new_list_element%field%s_ptr => target_element%field%s_ptr
    !
    IF (.NOT. ASSOCIATED(new_list_element%field%s_ptr)) THEN
      WRITE (0,*) 'problem with association of ptr for '//TRIM(name)
    ENDIF
    !
    IF(PRESENT(info)) info => new_list_element%field%info
    IF (PRESENT(lmiss)) THEN
      ptr = new_list_element%field%info%missval%sval
    ELSE
      ptr = 0.0_sp
    END IF
  END SUBROUTINE add_var_list_reference_s2d

  ! INTEGER SECTION ----------------------------------------------------------------------------------
  !
  ! create (allocate) a new table entry
  ! reference to an existing pointer to 3d-field
  ! optionally overwrite some default meta data
  !
  SUBROUTINE add_var_list_reference_i2d (this_list, target_name, name, ptr,                      &
       &                                 hgrid, vgrid, cf, grib2, ref_idx, ldims, loutput,       &
       &                                 lrestart, lrestart_cont, initval, isteptype,            &
       &                                 resetval, lmiss, missval, tlev_source, tracer_info,     &
       &                                 info, vert_interp, hor_interp, in_group, verbose,       &
       &                                 new_element, l_pp_scheduler_task, post_op, action_list, &
       &                                 opt_var_ref_pos, var_class)
    !
    TYPE(t_var_list),        INTENT(inout)        :: this_list
    CHARACTER(len=*),        INTENT(in)           :: target_name
    CHARACTER(len=*),        INTENT(in)           :: name
    INTEGER, POINTER                              :: ptr(:,:)
    INTEGER,                 INTENT(in)           :: hgrid                        ! horizontal grid type used
    INTEGER,                 INTENT(in)           :: vgrid                        ! vertical grid type used
    TYPE(t_cf_var),          INTENT(in)           :: cf                           ! CF related metadata
    TYPE(t_grib2_var),       INTENT(in)           :: grib2                        ! GRIB2 related metadata
    INTEGER,                 INTENT(in)           :: ref_idx                      ! idx of slice to be referenced
    INTEGER,                 INTENT(in)           :: ldims(2)                     ! local dimensions, for checking
    LOGICAL,                 INTENT(in), OPTIONAL :: loutput                      ! output flag
    LOGICAL,                 INTENT(in), OPTIONAL :: lrestart                     ! restart flag
    LOGICAL,                 INTENT(in), OPTIONAL :: lrestart_cont                ! continue restart if var not available
    INTEGER,                 INTENT(in), OPTIONAL :: initval                      ! value if var not available
    INTEGER,                 INTENT(in), OPTIONAL :: isteptype                    ! type of statistical processing
    INTEGER,                 INTENT(in), OPTIONAL :: resetval                     ! reset value (after accumulation)
    LOGICAL,                 INTENT(in), OPTIONAL :: lmiss                        ! missing value flag
    INTEGER,                 INTENT(in), OPTIONAL :: missval                      ! missing value
    INTEGER,                 INTENT(in), OPTIONAL :: tlev_source                  ! actual TL for TL dependent vars
    CLASS(t_tracer_meta),    INTENT(in), OPTIONAL :: tracer_info                  ! tracer meta data
    TYPE(t_var_metadata), POINTER,       OPTIONAL :: info                         ! returns reference to metadata
    TYPE(t_vert_interp_meta),INTENT(in), OPTIONAL :: vert_interp                  ! vertical interpolation metadata
    TYPE(t_hor_interp_meta), INTENT(in), OPTIONAL :: hor_interp                   ! horizontal interpolation metadata
    LOGICAL,                 INTENT(in), OPTIONAL :: in_group(:)                  ! groups to which a variable belongs
    LOGICAL,                 INTENT(in), OPTIONAL :: verbose
    TYPE(t_list_element), POINTER,       OPTIONAL :: new_element                  ! pointer to new var list element
    INTEGER,                 INTENT(in), OPTIONAL :: l_pp_scheduler_task          ! .TRUE., if field is updated by pp scheduler
    TYPE(t_post_op_meta),    INTENT(IN), OPTIONAL :: post_op                      !< "post-op" (small arithmetic operations) for this variable
    TYPE(t_var_action),      INTENT(IN), OPTIONAL :: action_list                  !< regularly triggered events
    INTEGER,                 INTENT(IN), OPTIONAL :: opt_var_ref_pos              !< (optional:) position of container index
    INTEGER,                 INTENT(in), OPTIONAL :: var_class                    !< variable type/species
    ! local variables
    CHARACTER(*), PARAMETER :: routine = modname//"::add_var_list_reference_i2d"
    !
    TYPE(t_list_element), POINTER :: target_element
    TYPE(t_var_metadata), POINTER :: target_info, ref_info
    TYPE(t_list_element), POINTER :: new_list_element
    TYPE(t_union_vals)            :: missvalt, initvalt, resetvalt
    INTEGER                       :: ndims, var_ref_pos, dim_indices(5), index
    LOGICAL :: in_group_new(MAX_GROUPS)             ! groups to which a variable belongs
                                                    ! (for taking into account tile groups)
    !
    ndims = 2

    target_element => find_list_element (this_list, target_name)
    target_info => target_element%field%info
    IF (.NOT. ASSOCIATED(target_element%field%i_ptr))  CALL finish(routine, TRIM(name)//' not created.')
    !
    ! The parameter "var_ref_pos" contains the dimension index which
    ! points to the reference slice. Usually, this is "ndims+1", such
    ! that 3D slices, e.g., are stored in a 4D array as (:,:,:,1),
    ! (:,:,:,2), (:,:,:,3), etc.
    IF (PRESENT(opt_var_ref_pos)) THEN
      var_ref_pos    = opt_var_ref_pos
      IF (.NOT. target_info%lcontainer) &
        &  CALL finish(routine, "Container index does not make sense: Target is not a container variable!")
      IF ((target_info%var_ref_pos /= var_ref_pos) .AND. (target_info%var_ref_pos /= -1)) THEN
        CALL finish(routine, "Container index does not match the previously set value!")
      END IF
      target_info%var_ref_pos = var_ref_pos
    ELSE
      var_ref_pos    = ndims + 1
    END IF
    SELECT CASE(var_ref_pos)
    CASE (1)
      dim_indices    = (/ 2, 3, 0, 0, 0 /)
    CASE (2)
      dim_indices    = (/ 1, 3, 0, 0, 0 /)
    CASE (3)
      dim_indices    = (/ 1, 2, 0, 0, 0 /)
    CASE DEFAULT
      CALL finish(routine, "Internal error!")
    END SELECT

    IF (target_info%lcontainer) THEN
      ! Counting the number of existing references is deactivated, if the slice index
      ! to be referenced is given explicitly.
!!$      IF ( PRESENT(ref_idx) ) THEN
        target_info%ncontained = target_info%ncontained+1
        ! only check validity of given slice index
        IF ( (ref_idx > SIZE(target_element%field%i_ptr, var_ref_pos)) .OR. (ref_idx < 1)) THEN
          WRITE (message_text, *) &
            &  'Slice idx ', ref_idx, ' for ', TRIM(name), &
            &  ' out of allowable range [1,',SIZE(target_element%field%i_ptr, var_ref_pos),']'
          CALL finish(routine, message_text)
        ENDIF
!!$      ELSE
!!$        target_info%ncontained = target_info%ncontained+1
!!$        IF (SIZE(target_element%field%i_ptr, var_ref_pos) < target_info%ncontained) THEN
!!$          WRITE (message_text, *) &
!!$            &  TRIM(name), ' exceeds the number of predefined entries in container:', &
!!$            &  SIZE(target_element%field%i_ptr, var_ref_pos)
!!$          CALL finish(routine, message_text)
!!$        ENDIF
!!$      ENDIF
      IF (any(ldims(1:ndims) /=  target_info%used_dimensions(dim_indices(1:ndims)))) THEN
        CALL finish(routine, TRIM(name)//' dimensions requested and available differ.')
      ENDIF
    ENDIF
    !
    ! add list entry
    !
    CALL append_list_element (this_list, new_list_element)
    IF (PRESENT(new_element)) new_element=>new_list_element
    ref_info => new_list_element%field%info
    CALL default_var_list_metadata(ref_info, this_list)
    !
    ! init local fields
    !
    missvalt  = ref_info%missval
    initvalt  = ref_info%initval
    resetvalt = ref_info%resetval
    IF (PRESENT(missval))  missvalt%ival  = missval
    IF (PRESENT(initval))  initvalt%ival  = initval
    IF (PRESENT(resetval)) resetvalt%ival = resetval
    CALL set_var_metadata (new_list_element%field%info,                      &
         name=name, hgrid=hgrid, vgrid=vgrid,                                &
         cf=cf, grib2=grib2, ldims=ldims, loutput=loutput,                   &
         lrestart=lrestart, lrestart_cont=lrestart_cont, initval=initvalt,   &
         isteptype=isteptype, resetval=resetvalt, lmiss=lmiss,               &
         missval=missvalt, tlev_source=tlev_source,                          &
         vert_interp=vert_interp, hor_interp=hor_interp,                     &
         in_group=in_group, verbose=verbose,                                 &
         l_pp_scheduler_task=l_pp_scheduler_task,                            &
         post_op=post_op, action_list=action_list, var_class=var_class,      &
         data_type=INT_T )
    ! set dynamic metadata, i.e. polymorphic tracer metadata
    CALL set_var_metadata_dyn (new_list_element%field%info_dyn,              &
                               tracer_info=tracer_info)
    !
    ref_info%ndims = ndims
    ref_info%used_dimensions(:)       = 0
    ref_info%used_dimensions(1:ndims) = target_element%field%info%used_dimensions(dim_indices(1:ndims))

    IF (PRESENT(var_class)) THEN
      IF ( ANY((/CLASS_TILE, CLASS_TILE_LAND/) == var_class)) THEN
        ! automatically add tile to its variable specific tile-group
        CALL var_groups_dyn%add(group_name=target_name, in_group_new=in_group_new, opt_in_group=in_group)
        !
        ! update in_group metainfo
        new_list_element%field%info%in_group(:) = in_group_new(:)
      ENDIF
    END IF
    !
    index = 1
    IF (target_info%lcontainer) THEN
      ref_info%lcontained                   = .TRUE.
      ref_info%used_dimensions(ndims+1)     = 1
      ref_info%var_ref_pos                  = var_ref_pos
      !
      ref_info%maxcontained = SIZE(target_element%field%i_ptr,var_ref_pos)
      !
!!$      IF ( PRESENT(ref_idx) ) THEN
        ref_info%ncontained = ref_idx
!!$      ELSE
!!$        ref_info%ncontained = target_info%ncontained
!!$      ENDIF
      index = ref_info%ncontained
    ENDIF
    SELECT CASE(var_ref_pos)
    CASE(1)
      ptr => target_element%field%i_ptr(index,:,:,1,1)
    CASE(2)
      ptr => target_element%field%i_ptr(:,index,:,1,1)
    CASE(3)
      ptr => target_element%field%i_ptr(:,:,index,1,1)
    CASE default
      CALL finish(routine, "internal error!")
    END SELECT
    new_list_element%field%i_ptr => target_element%field%i_ptr
    !
    IF (.NOT. ASSOCIATED(new_list_element%field%i_ptr)) THEN
      WRITE (0,*) 'problem with association of ptr for '//TRIM(name)
    ENDIF
    !
    IF(PRESENT(info)) info => new_list_element%field%info
    IF (PRESENT(lmiss)) THEN
      ptr = new_list_element%field%info%missval%ival
    ELSE
      ptr = 0
    END IF
  END SUBROUTINE add_var_list_reference_i2d

  !================================================================================================
  !------------------------------------------------------------------------------------------------
  !
  ! perform consistency checks on variable's meta-data.
  !
  SUBROUTINE check_metadata_consistency(info)
    TYPE(t_var_metadata), INTENT(IN) :: info  ! variable meta data
    CHARACTER(LEN=*), PARAMETER :: routine = modname//':check_metadata_consistency'

    IF (info%lrestart .AND. info%lcontainer) THEN
      CALL finish(routine//' - '//TRIM(info%name), &
        &         'Container variables are not restartable! Use var references instead.')
    END IF
    ! ... put other consistency checks here ...
  END SUBROUTINE check_metadata_consistency
  !------------------------------------------------------------------------------------------------
  !
  ! Print routines for control output and debuggung
  !
  SUBROUTINE print_memory_use (this_list, ldetailed)
    TYPE(t_var_list) ,INTENT(in) :: this_list ! list
    LOGICAL, INTENT(in), OPTIONAL :: ldetailed
    !
    IF (PRESENT(ldetailed)) THEN
      WRITE (message_text,'(a32,a,a,i10,a,i4,a)')                    &
           TRIM(this_list%p%name), '-buffer: ',                      &
           'Memory in use: ', this_list%p%memory_used, ' bytes in ', &
           this_list%p%list_elements, ' fields.'
    ELSE
      WRITE (message_text,'(a32,a,a,i10,a,i6,a)')                         &
           TRIM(this_list%p%name), '-buffer: ',                           &
           'Memory in use: ', this_list%p%memory_used/1024_i8, ' kb in ', &
           this_list%p%list_elements, ' fields.'
    ENDIF
    CALL message('',message_text)
    !
  END SUBROUTINE print_memory_use
  !------------------------------------------------------------------------------------------------
  !
  ! print current memory table
  !
  SUBROUTINE print_var_list (this_list, lshort)
    TYPE(t_var_list),  INTENT(in) :: this_list ! list
    LOGICAL, OPTIONAL :: lshort
    !
    TYPE(t_list_element), POINTER :: this_list_element
    CHARACTER(len=32) :: dimension_text, dtext,keytext
    INTEGER :: i, igrp, ivintp_type
    CHARACTER(len=4) :: localMode = '----'
    LOGICAL :: short

    short = .FALSE.
    IF (PRESENT(lshort)) short = lshort
    CALL message('','')
    CALL message('','')
    CALL message('','Status of variable list '//TRIM(this_list%p%name)//':')
    CALL message('','')
    !
    this_list_element => this_list%p%first_list_element
    !
    DO WHILE (ASSOCIATED(this_list_element))
      !
      IF (short) THEN

        IF (this_list_element%field%info%name /= '' .AND. &
             .NOT. this_list_element%field%info%lcontainer) THEN
          IF (this_list_element%field%info%lrestart) localMode(1:1) = 'r'
          IF (this_list_element%field%info%lcontained) localMode(2:2) = 't'
          SELECT CASE (this_list_element%field%info%isteptype)
          CASE (1)
            localMode(3:3) = 'i'
          CASE (2)
            localMode(3:3) = 'm'
          CASE (3)
            localMode(3:3) = 'a'
          END SELECT
          SELECT CASE (this_list_element%field%info%hgrid)
          CASE (1)
            localMode(4:4) = 'c'
          CASE (2)
            localMode(4:4) = 'v'
          CASE (3)
            localMode(4:4) = 'e'
          END SELECT

          WRITE(message_text, '(a4,3i4,a24,a40)') localMode,                                 &
               &                              this_list_element%field%info%grib2%discipline, &
               &                              this_list_element%field%info%grib2%category,   &
               &                              this_list_element%field%info%grib2%number,     &
               &                              TRIM(this_list_element%field%info%name),       &
               &                              TRIM(this_list_element%field%info%cf%standard_name)
          CALL message('', message_text)

          localMode = '----'
        ENDIF

      ELSE

      IF (this_list_element%field%info%name /= '' .AND. &
           .NOT. this_list_element%field%info%lcontainer) THEN
        !
        WRITE (message_text,'(a,a)')       &
             'Table entry name                            : ', &
             TRIM(this_list_element%field%info%name)
        CALL message('', message_text)
        WRITE (keytext,'(i32.1)') this_list_element%field%info%key

        WRITE (message_text,'(a,a)')       &
             'Key entry                                   : ', &
             TRIM(keytext)
        CALL message('', message_text)
        !
        IF (ASSOCIATED(this_list_element%field%r_ptr) .OR. &
          & ASSOCIATED(this_list_element%field%s_ptr) .OR. &
          & ASSOCIATED(this_list_element%field%i_ptr) .OR. &
          & ASSOCIATED(this_list_element%field%l_ptr)) THEN
          CALL message ('','Pointer status                              : in use.')
          dimension_text = '('
          DO i = 1, this_list_element%field%info%ndims
            WRITE(dtext,'(i0)') this_list_element%field%info%used_dimensions(i)
            IF (this_list_element%field%info%ndims == i) THEN
              dimension_text = TRIM(dimension_text)//TRIM(dtext)//')'
            ELSE
              dimension_text = TRIM(dimension_text)//TRIM(dtext)//','
            ENDIF
          ENDDO
          WRITE (message_text,'(a,a)') &
               'Local field dimensions                      : ', TRIM(dimension_text)
          CALL message('', message_text)
        ELSE
          CALL message('', 'Pointer status                              : not in use.')
        ENDIF
        !
        WRITE (message_text,'(a,3i4)') &
             'Assigned GRIB discipline/category/parameter : ', &
             this_list_element%field%info%grib2%discipline,    &
             this_list_element%field%info%grib2%category,      &
             this_list_element%field%info%grib2%number
        CALL message('', message_text)
        !
        WRITE (message_text,'(a,a,a,a)')                          &
             'CF convention standard name/unit            : ',    &
             TRIM(this_list_element%field%info%cf%standard_name), &
             '     ',                                             &
             TRIM(this_list_element%field%info%cf%units)
        CALL message('', message_text)
        !
        WRITE (message_text,'(2a)') &
             'CF convention long name                     : ', &
             TRIM(this_list_element%field%info%cf%long_name)
        !
        IF (this_list_element%field%info%lcontained) THEN
          CALL message('', 'Field is in a container                     : yes.')
          WRITE (message_text,'(a,i2)')                        &
             ' Index in container                          : ',&
             this_list_element%field%info%ncontained
          CALL message('', message_text)
        ELSE
          CALL message('', 'Field is in a container                     : no.')
          WRITE (message_text,'(a)')                           &
             ' Index in container                          : --'
          CALL message('', message_text)
        ENDIF
        !
        WRITE (message_text,'(a,i2)')                          &
             ' horizontal grid type used (C=1,V=2,E=3)     : ',&
             this_list_element%field%info%hgrid
        CALL message('', message_text)
        !
        WRITE (message_text,'(a,i2)')                          &
             ' vertical grid type used (see cdilib.c)      : ',&
             this_list_element%field%info%vgrid
        CALL message('', message_text)
        !
        WRITE (message_text,'(a,i2)')                          &
             ' type of stat. processing (I=1,AVG=2,ACC=3...: ',&
             this_list_element%field%info%isteptype
        CALL message('', message_text)
        !
        IF (this_list_element%field%info%lmiss) THEN
          IF (ASSOCIATED(this_list_element%field%r_ptr)) THEN
            WRITE (message_text,'(a,e20.12)')      &
                 'Missing value                               : ', &
                 this_list_element%field%info%missval%rval
          ELSE IF (ASSOCIATED(this_list_element%field%s_ptr)) THEN
            WRITE (message_text,'(a,e20.12)')      &
                 'Missing value                               : ', &
                 this_list_element%field%info%missval%sval
          ELSE IF (ASSOCIATED(this_list_element%field%i_ptr)) THEN
            WRITE (message_text,'(a,i8)')      &
                 'Missing value                               : ', &
                 this_list_element%field%info%missval%ival
          ELSE IF (ASSOCIATED(this_list_element%field%l_ptr)) THEN
            WRITE (message_text,'(a,l8)')      &
                 'Missing value                               : ', &
                 this_list_element%field%info%missval%lval
          ENDIF
          CALL message('', message_text)
        ELSE
          CALL message('', 'Missing values                              : off.')
        ENDIF
        !
        IF (this_list_element%field%info%lrestart) THEN
          CALL message('', 'Added to restart                            : yes.')
        ELSE
          CALL message('', 'Added to Restart                            : no.')
        ENDIF
        !
        IF (this_list_element%field%info_dyn%tracer%lis_tracer) THEN
          CALL message('', 'Tracer field                                : yes.')

          IF (this_list_element%field%info_dyn%tracer%lfeedback) THEN
            CALL message('', 'Child-to-parent feedback                  : yes.')
          ELSE
            CALL message('', 'Child-to-parent feedback                  : no.')
          ENDIF

          WRITE (message_text,'(a,3i3)') &
             'Horizontal transport method                 : ', &
             this_list_element%field%info_dyn%tracer%ihadv_tracer
          CALL message('', message_text)

          WRITE (message_text,'(a,3i3)') &
             'Vertical transport method                   : ', &
             this_list_element%field%info_dyn%tracer%ivadv_tracer
          CALL message('', message_text)

          IF (this_list_element%field%info_dyn%tracer%lturb_tracer) THEN
            CALL message('', 'Turbulent transport                         : yes.')
          ELSE
            CALL message('', 'Turbulent transport                         : no.')
          ENDIF

        ELSE
          CALL message('', 'Tracer field                                : no.')
        ENDIF !lis_tracer

        ! print variable class/species
        WRITE (message_text,'(a,i2)')       &
             'Variable class/species                      : ', &
             this_list_element%field%info%var_class
        CALL message('', message_text)

        !
        ! print groups, to which this variable belongs:
        IF (ANY(this_list_element%field%info%in_group(:))) THEN
          WRITE (message_text,'(a)')  'Variable group(s)                           :'
          DO igrp=1,SIZE(this_list_element%field%info%in_group)
            IF (this_list_element%field%info%in_group(igrp)) THEN
              IF (igrp == 1) THEN
                message_text = TRIM(message_text)//" "//TRIM(var_groups_dyn%name(igrp))
              ELSE
                message_text = TRIM(message_text)//", "//TRIM(var_groups_dyn%name(igrp))
              END IF
            ENDIF
          END DO
          CALL message('', message_text)
        END IF

        !
        ! print horizontal and vertical interpolation method(s):
        WRITE (message_text,'(a)')  &
          &  'Horizontal interpolation                    : '//  &
          &  TRIM(STR_HINTP_TYPE(this_list_element%field%info%hor_interp%hor_intp_type))
        CALL message('', message_text)

        LOOP_VINTP_TYPES : DO ivintp_type=1,SIZE(VINTP_TYPE_LIST)
          IF (this_list_element%field%info%vert_interp%vert_intp_type(ivintp_type)) THEN
            WRITE (message_text,'(a)')  &
              &  'Vertical interpolation                      : '//  &
              &  toupper(TRIM(VINTP_TYPE_LIST(ivintp_type)))
            CALL message('', message_text)
          END IF
        END DO LOOP_VINTP_TYPES
        CALL message('', '')
      ENDIF

      ENDIF
      !
      ! select next element in linked list
      !
      this_list_element => this_list_element%next_list_element
    ENDDO

    !
  END SUBROUTINE print_var_list

  !------------------------------------------------------------------------------------------------
  LOGICAL FUNCTION elementFoundByName(key2look4,name2look4,name_has_time_level,element,case_insensitive)
    INTEGER, INTENT(in) :: key2look4
    CHARACTER(len=*),   INTENT(in) :: name2look4
    TYPE(t_list_element), INTENT(in) :: element
    LOGICAL, INTENT(in) :: name_has_time_level, case_insensitive

    ! go forward only if both variables have NO or THE SAME timelevel
    IF (name_has_time_level .NEQV. has_time_level(element%field%info%name)) THEN
      elementFoundByName = .FALSE.
      RETURN
    ENDIF

    IF (case_insensitive) THEN
      elementFoundByName &
        = tolower(name2look4) == tolower(get_var_name(element%field))
    ELSE
      ! fixme: unless perfect hashing can be employed, this
      ! might create false positives
      elementFoundByName = key2look4 == element%field%info%key
    END IF
  END FUNCTION elementFoundByName
  !-----------------------------------------------------------------------------
  
  ! Should be overloaded to be able to search for the different information 
  ! In the proposed structure for the linked list, in the example only
  ! A character string is used so it is straight forward only one find
  !
  FUNCTION find_list_element (this_list, name, opt_hgrid, opt_caseInsensitive) RESULT(element)
    TYPE(t_var_list),   INTENT(in) :: this_list
    CHARACTER(len=*),   INTENT(in) :: name
    INTEGER, OPTIONAL              :: opt_hgrid
    LOGICAL, OPTIONAL              :: opt_caseInsensitive
    TYPE(t_list_element), POINTER :: element
    INTEGER :: key,hgrid
    LOGICAL :: name_has_time_level, case_insensitive

    case_insensitive = .FALSE.
    IF (PRESENT(opt_caseInsensitive)) case_insensitive = opt_caseInsensitive
    hgrid = -1
    IF (PRESENT(opt_hgrid)) hgrid = opt_hgrid
    key = text_hash_c(TRIM(name))
    name_has_time_level = has_time_level(name)
    element => this_list%p%first_list_element
    DO WHILE (ASSOCIATED(element))
      IF (-1 == hgrid .OR. hgrid == element%field%info%hgrid) THEN
        IF (elementFoundByName(key,name,name_has_time_level,&
          &                    element,case_insensitive)) RETURN
      ENDIF
      element => element%next_list_element
    ENDDO
  END FUNCTION find_list_element
  
  !------------------------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  !
  ! Overloaded to search for a tracer by its index (ncontained)
  !
  FUNCTION find_tracer_by_index (this_list, ncontained, opt_hgrid) RESULT(this_list_element)
    TYPE(t_var_list),   INTENT(in) :: this_list
    INTEGER,            INTENT(in) :: ncontained
    INTEGER, OPTIONAL              :: opt_hgrid
    TYPE(t_list_element), POINTER  :: this_list_element
    INTEGER :: hgrid

    hgrid = -1
    IF (PRESENT(opt_hgrid)) hgrid = opt_hgrid
    this_list_element => this_list%p%first_list_element
    DO WHILE (ASSOCIATED(this_list_element))
      IF (this_list_element%field%info_dyn%tracer%lis_tracer) THEN
        IF(ncontained == this_list_element%field%info%ncontained) THEN
          IF (-1 == hgrid) THEN
            RETURN
          ELSE
            IF (hgrid == this_list_element%field%info%hgrid) RETURN
          ENDIF
        ENDIF
      ENDIF
      this_list_element => this_list_element%next_list_element
    ENDDO
    NULLIFY (this_list_element)
  END FUNCTION find_tracer_by_index
END MODULE mo_var_list
