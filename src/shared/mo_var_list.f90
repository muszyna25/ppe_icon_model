! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
MODULE mo_var_list

#include <icon_contiguous_defines.h>
#include <omp_definitions.inc>

#if defined (__INTEL_COMPILER) || defined (__PGI) || defined (NAGFOR)
#ifdef VARLIST_INITIZIALIZE_WITH_NAN
  USE, INTRINSIC :: ieee_arithmetic
#endif
#endif

  USE mo_kind,             ONLY: sp, dp, i8
  USE mo_cf_convention,    ONLY: t_cf_var
  USE mo_grib2,            ONLY: t_grib2_var, grib2_var
  USE mo_var_groups,       ONLY: var_groups_dyn, groups, MAX_GROUPS
  USE mo_var_metadata_types,ONLY: t_var_metadata, t_union_vals,     &
    &                            t_vert_interp_meta, CLASS_TILE,    &
    &                            t_hor_interp_meta, t_post_op_meta, &
    &                            VINTP_TYPE_LIST, CLASS_TILE_LAND
  USE mo_var_metadata,     ONLY: create_vert_interp_metadata,       &
    &                            create_hor_interp_metadata, &
    &                            set_var_metadata, set_var_metadata_dyn, &
    &                            get_var_timelevel
  USE mo_tracer_metadata_types, ONLY: t_tracer_meta
  USE mo_var_list_element, ONLY: t_var_list_element, level_type_ml
  USE mo_exception,        ONLY: message, finish, message_text
  USE mo_util_texthash,    ONLY: text_hash_c
  USE mo_util_string,      ONLY: toupper, tolower
  USE mo_impl_constants,   ONLY: STR_HINTP_TYPE, REAL_T, SINGLE_T, BOOL_T, INT_T
  USE mo_fortran_tools,    ONLY: init_contiguous_dp, init_contiguous_sp, &
    &                            init_contiguous_i4, init_contiguous_l
  USE mo_action_types,     ONLY: t_var_action

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: print_var_list, add_var, add_ref
  PUBLIC :: find_list_element, delete_list, append_list_element
  PUBLIC :: t_var_list_ptr, t_list_element, sel_var_list

  ! t_list_element provides the entry to the actual information 
  ! and a reference to the next element in the list
  TYPE t_list_element
    TYPE(t_var_list_element)      :: field
    TYPE(t_list_element), POINTER :: next_list_element => NULL()
  END TYPE t_list_element

  TYPE t_var_list
    INTEGER                       :: id = -1
    CHARACTER(len=128)            :: name = ''          ! stream name
    TYPE(t_list_element), POINTER :: first_list_element => NULL() ! reference to first
    INTEGER(i8)                   :: memory_used = 0_i8 ! memory allocated
    INTEGER                       :: list_elements = 0  ! allocated elements
    LOGICAL :: loutput = .TRUE., lrestart = .FALSE., linitial = .FALSE.
    CHARACTER(len=256)            :: filename = ''      ! name of file
    CHARACTER(len=8) :: post_suf = '', rest_suf = '', init_suf = '' ! suffixes
    LOGICAL                       :: first = .FALSE.    ! first var_list in file
    INTEGER :: output_type = -1, restart_type = -1, compression_type = -1 ! CDI file/compression type
    LOGICAL :: restart_opened = .FALSE., output_opened = .FALSE. ! true, if restart/output file opened
    CHARACTER(len=8)              :: model_type = 'atm' ! store model type (default is 'atm' for reasons)
    INTEGER                       :: patch_id = -1      ! ID of patch to which list variables belong
    INTEGER                       :: vlevel_type = level_type_ml ! 1: modellevels, 2: pressure levels, 3: height levels
    INTEGER                       :: nvars = 0
    LOGICAL                    :: lmiss = .FALSE.         ! flag: true, if variables should be initialized with missval
    LOGICAL                    :: lmask_boundary = .TRUE. ! flag: true, if interpolation zone should be masked *in output*
  END TYPE t_var_list

  TYPE t_var_list_ptr
    TYPE(t_var_list), POINTER :: p => NULL()
  END type t_var_list_ptr

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

  CHARACTER(*), PARAMETER :: modname = "mo_var_list"

CONTAINS

  FUNCTION sel_var_list(obj) RESULT(vl)
    CLASS(*), POINTER, INTENT(IN) :: obj
    TYPE(t_var_list), POINTER :: vl

    NULLIFY(vl)
    IF (ASSOCIATED(obj)) THEN
      SELECT TYPE(obj)
      TYPE IS(t_var_list)
        vl => obj
      END SELECT
    END IF
  END FUNCTION sel_var_list

  !-----------------------------------------------------------------------------
  ! remove all elements of a linked list
  ! check if all elements are removed
  SUBROUTINE delete_list(this_list)
    TYPE(t_var_list_ptr), INTENT(INOUT) :: this_list
    TYPE(t_list_element), POINTER   :: this, next

    next => this_list%p%first_list_element
    DO WHILE (ASSOCIATED(next))
      this => next
      next => this%next_list_element
      IF (this%field%info%allocated) THEN
        IF (ASSOCIATED(this%field%r_ptr)) THEN
          !$ACC EXIT DATA DELETE(this%field%r_ptr) IF(this%field%info%lopenacc)
          DEALLOCATE (this%field%r_ptr)
        ELSE IF (ASSOCIATED(this%field%s_ptr)) THEN
          !$ACC EXIT DATA DELETE(this%field%s_ptr) IF(this%field%info%lopenacc)
          DEALLOCATE (this%field%s_ptr)
        ELSE IF (ASSOCIATED(this%field%i_ptr)) THEN
          !$ACC EXIT DATA DELETE(this%field%i_ptr) IF(this%field%info%lopenacc)
          DEALLOCATE (this%field%i_ptr)
        ELSE IF (ASSOCIATED(this%field%l_ptr)) THEN
          !$ACC EXIT DATA DELETE(this%field%l_ptr) IF(this%field%info%lopenacc)
          DEALLOCATE (this%field%l_ptr)
        ENDIF
      ENDIF
      DEALLOCATE (this)
    END DO
  END SUBROUTINE delete_list

  !-----------------------------------------------------------------------------
  ! add a list element to the linked list
  SUBROUTINE append_list_element(this_list, new_element)
    TYPE(t_var_list_ptr), INTENT(INOUT) :: this_list
    TYPE(t_list_element), POINTER, INTENT(OUT) :: new_element
    TYPE(t_list_element), POINTER :: cur_element

    IF (.NOT.ASSOCIATED(this_list%p%first_list_element)) THEN
      ! insert as first element if list is empty
      ALLOCATE(this_list%p%first_list_element)
      new_element => this_list%p%first_list_element
    ELSE
      ! loop over list elements to find position
      cur_element => this_list%p%first_list_element
      DO WHILE (ASSOCIATED(cur_element%next_list_element))
        cur_element => cur_element%next_list_element
      ENDDO
      ! insert element
      ALLOCATE(new_element)
      new_element%next_list_element => cur_element%next_list_element
      cur_element%next_list_element => new_element
    END IF
    this_list%p%list_elements = this_list%p%list_elements+1
    this_list%p%nvars = this_list%p%nvars + 1
  END SUBROUTINE append_list_element

  SUBROUTINE inherit_var_list_metadata(this_info, this_list)
    TYPE(t_var_metadata), INTENT(out) :: this_info
    TYPE(t_var_list_ptr), INTENT(in) :: this_list

    this_info%grib2               = grib2_var(-1, -1, -1, -1, -1, -1)
    this_info%lrestart            = this_list%p%lrestart
    this_info%lmiss               = this_list%p%lmiss
    this_info%lmask_boundary      = this_list%p%lmask_boundary
    this_info%vert_interp         = create_vert_interp_metadata()
    this_info%hor_interp          = create_hor_interp_metadata()
    this_info%in_group(:)         = groups()
  END SUBROUTINE inherit_var_list_metadata

  !------------------------------------------------------------------------------------------------
  ! Create a list new entry
  SUBROUTINE add_var_list_element_5d(data_type, this_list, varname,         &
    &   hgrid, vgrid, cf, grib2, ldims, new_elem, loutput, lcontainer,  &
    &   lrestart, lrestart_cont, isteptype, lmiss, tlev_source,                 &
    &   info, vert_interp, hor_interp, in_group, verbose,                       &
    &   l_pp_scheduler_task, post_op, action_list, tracer_info,                 &
    &   p5_r, p5_s, p5_i, p5_l, initval_r, initval_s, initval_i, initval_l,     &
    &   resetval_r, resetval_s, resetval_i, resetval_l, new_element,           &
    &   missval_r, missval_s, missval_i, missval_l, var_class, lopenacc )
    INTEGER, INTENT(IN) :: data_type, hgrid, vgrid, ldims(:)
    TYPE(t_var_list_ptr), INTENT(INOUT) :: this_list
    CHARACTER(*), INTENT(IN) :: varname
    TYPE(t_cf_var), INTENT(IN) :: cf
    TYPE(t_grib2_var), INTENT(IN) :: grib2
    TYPE(t_list_element), POINTER, INTENT(OUT) :: new_elem
    LOGICAL, INTENT(IN), OPTIONAL :: loutput, lcontainer, lrestart, lrestart_cont, &
      & lmiss, in_group(:), verbose, initval_l, resetval_l, missval_l, lopenacc
    INTEGER, INTENT(IN), OPTIONAL :: isteptype, tlev_source, l_pp_scheduler_task, &
      & initval_i, resetval_i, missval_i, var_class
    TYPE(t_var_metadata), POINTER, OPTIONAL :: info
    REAL(dp), CONTIGUOUS_TARGET, OPTIONAL :: p5_r(:,:,:,:,:)
    REAL(sp), CONTIGUOUS_TARGET, OPTIONAL :: p5_s(:,:,:,:,:)
    INTEGER, CONTIGUOUS_TARGET, OPTIONAL :: p5_i(:,:,:,:,:)
    LOGICAL, CONTIGUOUS_TARGET, OPTIONAL :: p5_l(:,:,:,:,:)
    TYPE(t_vert_interp_meta), INTENT(IN), OPTIONAL :: vert_interp
    TYPE(t_hor_interp_meta), INTENT(IN), OPTIONAL :: hor_interp
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op
    TYPE(t_var_action), INTENT(IN), OPTIONAL :: action_list
    CLASS(t_tracer_meta), INTENT(IN), OPTIONAL :: tracer_info
    REAL(dp), INTENT(IN), OPTIONAL :: initval_r, resetval_r, missval_r
    REAL(sp), INTENT(IN), OPTIONAL :: initval_s, resetval_s, missval_s
    TYPE(t_list_element), POINTER, INTENT(OUT), OPTIONAL :: new_element
    TYPE(t_union_vals) :: missval, initval, resetval, ivals
    INTEGER :: d(5), istat, ndims
    LOGICAL :: referenced, is_restart_var
    CHARACTER(LEN = *), PARAMETER :: routine = modname//":add_var_list_element_5d"

    ndims = SIZE(ldims)
    ! Check for a variable of the same name in this list
    ! This consistency check only makes sense inside individual lists.
    ! For single-domain setups and/or when using internal post-processing 
    ! (e.g. lon-lat or vertically interpolated output)  
    ! duplicate names may exist in different lists
    IF (ASSOCIATED(find_list_element(this_list, varname))) &
      & CALL finish(routine, "duplicate var-name ("//TRIM(varname)//") exists in var_list ("//TRIM(this_list%p%name)//")")
    is_restart_var = this_list%p%lrestart
    IF (PRESENT(lrestart)) THEN
      is_restart_var = lrestart
      IF (.NOT. this_list%p%lrestart .AND. lrestart) THEN
        CALL finish(routine, 'for list '//TRIM(this_list%p%name)//' restarting not enabled, '// &
                           & 'but restart of '//TRIM(varname)//' requested.')
      ENDIF
    ENDIF
    IF (is_restart_var .AND. (.NOT. ANY(data_type == (/REAL_T, SINGLE_T, INT_T/)))) THEN
      CALL finish(routine, 'unsupported data_type for "'//TRIM(varname)//'": '// &
        & 'data_type of restart variables must be floating-point or integer type.')
    END IF
    ! add list entry
    CALL append_list_element(this_list, new_elem)
    CALL inherit_var_list_metadata(new_elem%field%info, this_list)
    ! init local fields
    missval = new_elem%field%info%missval
    initval = new_elem%field%info%initval
    resetval= new_elem%field%info%resetval
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
    CALL set_var_metadata( new_elem%field%info,                      &
         name=varname, hgrid=hgrid, vgrid=vgrid, cf=cf, grib2=grib2,            &
         ldims=ldims(1:ndims), loutput=loutput, lcontainer=lcontainer,       &
         lrestart=lrestart, lrestart_cont=lrestart_cont, initval=initval,    &
         isteptype=isteptype, resetval=resetval, lmiss=lmiss,                &
         missval=missval, tlev_source=tlev_source,                           &
         vert_interp=vert_interp, hor_interp=hor_interp, in_group=in_group,  &
         verbose=verbose, l_pp_scheduler_task=l_pp_scheduler_task,           &
         post_op=post_op, action_list=action_list, var_class=var_class,      &
         data_type=data_type, lopenacc=lopenacc )
    ! set dynamic metadata, i.e. polymorphic tracer metadata
    CALL set_var_metadata_dyn (new_elem%field%info_dyn, tracer_info=tracer_info)
    new_elem%field%info%ndims = ndims
    new_elem%field%info%used_dimensions(1:ndims) = ldims(1:ndims)
    new_elem%field%info%dom => this_list%p%patch_id
    IF(PRESENT(info)) info => new_elem%field%info
    NULLIFY(new_elem%field%r_ptr, new_elem%field%s_ptr, &
      &     new_elem%field%i_ptr, new_elem%field%l_ptr)
    d(1:ndims)    = new_elem%field%info%used_dimensions(1:ndims)
    d((ndims+1):) = 1
#if    defined (VARLIST_INITIZIALIZE_WITH_NAN) \
    && (defined (__INTEL_COMPILER) || defined (__PGI) || defined (NAGFOR))
    ivals%rval = ieee_value(ptr, ieee_signaling_nan)
    ivals%sval = ieee_value(ptr, ieee_signaling_nan)
#endif
    IF (ANY([PRESENT(initval_r), PRESENT(initval_s), PRESENT(initval_i), PRESENT(initval_l)])) THEN
      ivals = initval
    ELSE IF (PRESENT(lmiss)) THEN
      ivals = missval
    END IF
    SELECT CASE(data_type)
    CASE(REAL_T)
      IF (referenced) THEN
        new_elem%field%r_ptr => p5_r
      ELSE
        new_elem%field%var_base_size = 8
        ALLOCATE(new_elem%field%r_ptr(d(1), d(2), d(3), d(4), d(5)), STAT=istat)
        IF (istat /= 0) CALL finish(routine, 'allocation of array '//TRIM(varname)//' failed')
        !$ACC ENTER DATA CREATE(new_elem%field%r_ptr) IF(new_elem%field%info%lopenacc)
      END IF
      !ICON_OMP PARALLEL
      CALL init_contiguous_dp(new_elem%field%r_ptr, PRODUCT(d(1:5)), ivals%rval)
      !ICON_OMP END PARALLEL
      !$ACC UPDATE DEVICE(new_elem%field%r_ptr) IF(new_elem%field%info%lopenacc)
    CASE(SINGLE_T)
      IF (referenced) THEN
        new_elem%field%s_ptr => p5_s
      ELSE
        new_elem%field%var_base_size = 4
        ALLOCATE(new_elem%field%s_ptr(d(1), d(2), d(3), d(4), d(5)), STAT=istat)
        IF (istat /= 0) CALL finish(routine, 'allocation of array '//TRIM(varname)//' failed')
        !$ACC ENTER DATA CREATE(new_elem%field%s_ptr) IF(new_elem%field%info%lopenacc)
      END IF
      !ICON_OMP PARALLEL
      CALL init_contiguous_sp(new_elem%field%s_ptr, PRODUCT(d(1:5)), ivals%sval)
      !ICON_OMP END PARALLEL
      !$ACC UPDATE DEVICE(new_elem%field%s_ptr) IF(new_elem%field%info%lopenacc)
    CASE(INT_T)
      IF (referenced) THEN
        new_elem%field%i_ptr => p5_i
      ELSE
        new_elem%field%var_base_size = 4
        ALLOCATE(new_elem%field%i_ptr(d(1), d(2), d(3), d(4), d(5)), STAT=istat)
        IF (istat /= 0) CALL finish(routine, 'allocation of arrayb'//TRIM(varname)//' failed')
        !$ACC ENTER DATA CREATE(new_elem%field%i_ptr)bIF(new_elem%field%info%lopenacc)
      END IF
      !ICON_OMP PARALLEL
      CALL init_contiguous_i4(new_elem%field%i_ptr, PRODUCT(d(1:5)), ivals%ival)
      !ICON_OMP END PARALLEL
      !$ACC UPDATE DEVICE(new_elem%field%i_ptr) IF(new_elem%field%info%lopenacc)
    CASE(BOOL_T)
      IF (referenced) THEN
        new_elem%field%l_ptr => p5_l
      ELSE
        new_elem%field%var_base_size = 4
        ALLOCATE(new_elem%field%l_ptr(d(1), d(2), d(3), d(4), d(5)), STAT=istat)
        IF (istat /= 0) CALL finish(routine, 'allocation of array '//TRIM(varname)//' failed')
        !$ACC ENTER DATA CREATE(new_elem%field%l_ptr) IF(new_elem%field%info%lopenacc)
      END IF
      !ICON_OMP PARALLEL
      CALL init_contiguous_l(new_elem%field%l_ptr, PRODUCT(d(1:5)), ivals%lval)
      !ICON_OMP END PARALLEL
      !$ACC UPDATE DEVICE(new_elem%field%l_ptr) IF(new_elem%field%info%lopenacc)
    END SELECT
    IF (.NOT.referenced) this_list%p%memory_used = this_list%p%memory_used + &
      & INT(new_elem%field%var_base_size, i8) * INT(PRODUCT(d(1:5)),i8)
    new_elem%field%info%allocated = .TRUE.
    IF (PRESENT(new_element)) new_element => new_elem
  END SUBROUTINE add_var_list_element_5d

  SUBROUTINE add_var_list_element_r4d(this_list, varname, ptr, hgrid, vgrid, &
    & cf, grib2, ldims, loutput, lcontainer, lrestart, lrestart_cont,     &
    & initval, isteptype, resetval, lmiss, missval, tlev_source, info,    &
    & p5, vert_interp, hor_interp, in_group, verbose, new_element,        &
    & l_pp_scheduler_task, post_op, action_list, var_class, lopenacc)
    TYPE(t_var_list_ptr), INTENT(INOUT) :: this_list
    CHARACTER(len=*), INTENT(IN) :: varname
    REAL(dp), POINTER, INTENT(OUT) :: ptr(:,:,:,:)
    INTEGER, INTENT(IN) :: hgrid, vgrid, ldims(4)
    TYPE(t_cf_var), INTENT(IN) :: cf
    TYPE(t_grib2_var), INTENT(IN) :: grib2
    LOGICAL, INTENT(IN), OPTIONAL :: loutput, lcontainer, lrestart, &
      & lrestart_cont, lmiss, in_group(:), verbose, lopenacc
    INTEGER, INTENT(IN), OPTIONAL :: isteptype, l_pp_scheduler_task, &
      & tlev_source, var_class
    REAL(dp), INTENT(IN), OPTIONAL :: initval, resetval, missval
    TYPE(t_var_metadata), POINTER, OPTIONAL :: info
    REAL(dp), CONTIGUOUS_TARGET, OPTIONAL :: p5(:,:,:,:,:)
    TYPE(t_vert_interp_meta), INTENT(IN), OPTIONAL :: vert_interp
    TYPE(t_hor_interp_meta), INTENT(IN), OPTIONAL :: hor_interp
    TYPE(t_list_element), POINTER, OPTIONAL :: new_element
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op
    TYPE(t_var_action), INTENT(IN), OPTIONAL :: action_list
    TYPE(t_list_element), POINTER :: element

    CALL add_var_list_element_5d(REAL_T, this_list, varname,           &
      &   hgrid, vgrid, cf, grib2, ldims, element, loutput, lcontainer,    &
      &   lrestart, lrestart_cont, isteptype, lmiss, tlev_source, info,    &
      &   vert_interp, hor_interp, in_group, verbose,                      &
      &   l_pp_scheduler_task, post_op, action_list, p5_r=p5,              &
      &   initval_r=initval, resetval_r=resetval, missval_r=missval,       &
      &   var_class=var_class, lopenacc=lopenacc, new_element=new_element)
    ptr => element%field%r_ptr(:,:,:,:,1)
  END SUBROUTINE add_var_list_element_r4d

  SUBROUTINE add_var_list_element_r3d(this_list, varname, ptr, hgrid, vgrid, &
    & cf, grib2, ldims, loutput, lcontainer, lrestart, lrestart_cont,     &
    & initval, isteptype, resetval, lmiss, missval, tlev_source, info,    &
    & p5, vert_interp, hor_interp, in_group, verbose, new_element, tracer_info, &
    & l_pp_scheduler_task, post_op, action_list, var_class, lopenacc)
    TYPE(t_var_list_ptr), INTENT(INOUT) :: this_list
    CHARACTER(len=*), INTENT(IN) :: varname
    REAL(dp), POINTER, INTENT(OUT) :: ptr(:,:,:)
    INTEGER, INTENT(IN) :: hgrid, vgrid, ldims(3)
    TYPE(t_cf_var), INTENT(IN) :: cf
    TYPE(t_grib2_var), INTENT(IN) :: grib2
    LOGICAL, INTENT(IN), OPTIONAL :: loutput, lcontainer, lrestart, &
      & lrestart_cont, lmiss, in_group(:), verbose, lopenacc
    INTEGER, INTENT(IN), OPTIONAL :: isteptype, l_pp_scheduler_task, &
      & tlev_source, var_class
    REAL(dp), INTENT(IN), OPTIONAL :: initval, resetval, missval
    TYPE(t_var_metadata), POINTER, OPTIONAL :: info
    CLASS(t_tracer_meta), INTENT(in), OPTIONAL :: tracer_info
    REAL(dp), CONTIGUOUS_TARGET, OPTIONAL :: p5(:,:,:,:,:)
    TYPE(t_vert_interp_meta), INTENT(IN), OPTIONAL :: vert_interp
    TYPE(t_hor_interp_meta), INTENT(IN), OPTIONAL :: hor_interp
    TYPE(t_list_element), POINTER, OPTIONAL :: new_element
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op
    TYPE(t_var_action), INTENT(IN), OPTIONAL :: action_list
    TYPE(t_list_element), POINTER :: element

    CALL add_var_list_element_5d(REAL_T, this_list, varname,               &
      &   hgrid, vgrid, cf, grib2, ldims, element, loutput, lcontainer,        &
      &   lrestart, lrestart_cont, isteptype, lmiss, tlev_source,              &
      &   info, vert_interp, hor_interp, in_group, verbose,                    &
      &   l_pp_scheduler_task, post_op, action_list, tracer_info=tracer_info,  &
      &   p5_r=p5, initval_r=initval, resetval_r=resetval, missval_r=missval,  &
      &   var_class=var_class, lopenacc=lopenacc, new_element=new_element)
    ptr => element%field%r_ptr(:,:,:,1,1)
  END SUBROUTINE add_var_list_element_r3d

  SUBROUTINE add_var_list_element_r2d(this_list, varname, ptr, hgrid, vgrid, &
    & cf, grib2, ldims, loutput, lcontainer, lrestart, lrestart_cont,     &
    & initval, isteptype, resetval, lmiss, missval, tlev_source, info,    &
    & p5, vert_interp, hor_interp, in_group, verbose, new_element, tracer_info, &
    & l_pp_scheduler_task, post_op, action_list, var_class, lopenacc)
    TYPE(t_var_list_ptr), INTENT(INOUT) :: this_list
    CHARACTER(len=*), INTENT(IN) :: varname
    REAL(dp), POINTER, INTENT(OUT) :: ptr(:,:)
    INTEGER, INTENT(IN) :: hgrid, vgrid, ldims(2)
    TYPE(t_cf_var), INTENT(IN) :: cf
    TYPE(t_grib2_var), INTENT(IN) :: grib2
    LOGICAL, INTENT(IN), OPTIONAL :: loutput, lcontainer, lrestart, &
      & lrestart_cont, lmiss, in_group(:), verbose, lopenacc
    INTEGER, INTENT(IN), OPTIONAL :: isteptype, l_pp_scheduler_task, &
      & tlev_source, var_class
    REAL(dp), INTENT(IN), OPTIONAL :: initval, resetval, missval
    TYPE(t_var_metadata), POINTER, OPTIONAL :: info
    CLASS(t_tracer_meta), INTENT(IN), OPTIONAL :: tracer_info
    REAL(dp), CONTIGUOUS_TARGET, OPTIONAL :: p5(:,:,:,:,:)
    TYPE(t_vert_interp_meta), INTENT(IN), OPTIONAL :: vert_interp
    TYPE(t_hor_interp_meta), INTENT(IN), OPTIONAL :: hor_interp
    TYPE(t_list_element), POINTER, OPTIONAL :: new_element
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op
    TYPE(t_var_action), INTENT(IN), OPTIONAL :: action_list
    TYPE(t_list_element), POINTER :: element

    CALL add_var_list_element_5d(REAL_T, this_list, varname,               &
      &   hgrid, vgrid, cf, grib2, ldims, element, loutput, lcontainer,        &
      &   lrestart, lrestart_cont, isteptype, lmiss, tlev_source,              &
      &   info, vert_interp, hor_interp, in_group, verbose,                    &
      &   l_pp_scheduler_task, post_op, action_list, tracer_info=tracer_info,  &
      &   p5_r=p5, initval_r=initval, resetval_r=resetval, missval_r=missval,  &
      &   var_class=var_class, lopenacc=lopenacc, new_element=new_element)
    ptr => element%field%r_ptr(:,:,1,1,1)
  END SUBROUTINE add_var_list_element_r2d

  SUBROUTINE add_var_list_element_r1d(this_list, varname, ptr, hgrid, vgrid, &
    & cf, grib2, ldims, loutput, lcontainer, lrestart, lrestart_cont,     &
    & initval, isteptype, resetval, lmiss, missval, tlev_source, info,    &
    & p5, vert_interp, hor_interp, in_group, verbose, new_element,        &
    & l_pp_scheduler_task, post_op, action_list, var_class, lopenacc)
    TYPE(t_var_list_ptr), INTENT(INOUT) :: this_list
    CHARACTER(len=*), INTENT(IN) :: varname
    REAL(dp), POINTER, INTENT(OUT) :: ptr(:)
    INTEGER, INTENT(IN) :: hgrid, vgrid, ldims(1)
    TYPE(t_cf_var), INTENT(IN) :: cf
    TYPE(t_grib2_var), INTENT(IN) :: grib2
    LOGICAL, INTENT(IN), OPTIONAL :: loutput, lcontainer, lrestart, &
      & lrestart_cont, lmiss, in_group(:), verbose, lopenacc
    INTEGER, INTENT(IN), OPTIONAL :: isteptype, l_pp_scheduler_task, &
      & tlev_source, var_class
    REAL(dp), INTENT(IN), OPTIONAL :: initval, resetval, missval
    TYPE(t_var_metadata), POINTER, OPTIONAL :: info
    REAL(dp), CONTIGUOUS_TARGET, OPTIONAL :: p5(:,:,:,:,:)
    TYPE(t_vert_interp_meta), INTENT(IN), OPTIONAL :: vert_interp
    TYPE(t_hor_interp_meta), INTENT(IN), OPTIONAL :: hor_interp
    TYPE(t_list_element), POINTER, OPTIONAL :: new_element
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op
    TYPE(t_var_action), INTENT(IN), OPTIONAL :: action_list
    TYPE(t_list_element), POINTER :: element

    CALL add_var_list_element_5d(REAL_T, this_list, varname,        &
      &   hgrid, vgrid, cf, grib2, ldims, element, loutput, lcontainer, &
      &   lrestart, lrestart_cont, isteptype, lmiss, tlev_source, info, &
      &   vert_interp, hor_interp, in_group, verbose,                   &
      &   l_pp_scheduler_task, post_op, action_list, p5_r=p5,           &
      &   initval_r=initval, resetval_r=resetval, missval_r=missval,    &
      &   var_class=var_class,lopenacc=lopenacc, new_element=new_element)
    ptr => element%field%r_ptr(:,1,1,1,1)
  END SUBROUTINE add_var_list_element_r1d

  SUBROUTINE add_var_list_element_s4d(this_list, varname, ptr, hgrid, vgrid, &
    & cf, grib2, ldims, loutput, lcontainer, lrestart, lrestart_cont,     &
    & initval, isteptype, resetval, lmiss, missval, tlev_source, info,    &
    & p5, vert_interp, hor_interp, in_group, verbose, new_element,        &
    & l_pp_scheduler_task, post_op, action_list, var_class, lopenacc)
    TYPE(t_var_list_ptr), INTENT(INOUT) :: this_list
    CHARACTER(len=*), INTENT(IN) :: varname
    REAL(sp), POINTER, INTENT(OUT) :: ptr(:,:,:,:)
    INTEGER, INTENT(IN) :: hgrid, vgrid, ldims(4)
    TYPE(t_cf_var), INTENT(IN) :: cf
    TYPE(t_grib2_var), INTENT(IN) :: grib2
    LOGICAL, INTENT(IN), OPTIONAL :: loutput, lcontainer, lrestart, &
      & lrestart_cont, lmiss, in_group(:), verbose, lopenacc
    INTEGER, INTENT(IN), OPTIONAL :: isteptype, l_pp_scheduler_task, &
      & tlev_source, var_class
    REAL(sp), INTENT(IN), OPTIONAL :: initval, resetval, missval
    TYPE(t_var_metadata), POINTER, OPTIONAL :: info
    REAL(sp), CONTIGUOUS_TARGET, OPTIONAL :: p5(:,:,:,:,:)
    TYPE(t_vert_interp_meta), INTENT(IN), OPTIONAL :: vert_interp
    TYPE(t_hor_interp_meta), INTENT(IN), OPTIONAL :: hor_interp
    TYPE(t_list_element), POINTER, OPTIONAL :: new_element
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op
    TYPE(t_var_action), INTENT(IN), OPTIONAL :: action_list
    TYPE(t_list_element), POINTER :: element

    CALL add_var_list_element_5d(SINGLE_T, this_list, varname,         &
      &   hgrid, vgrid, cf, grib2, ldims, element, loutput, lcontainer,    &
      &   lrestart, lrestart_cont, isteptype, lmiss, tlev_source, info,    &
      &   vert_interp, hor_interp, in_group, verbose,                      &
      &   l_pp_scheduler_task, post_op, action_list, p5_s=p5,              &
      &   initval_s=initval, resetval_s=resetval, missval_s=missval,       &
      &   var_class=var_class, lopenacc=lopenacc, new_element=new_element)
    ptr => element%field%s_ptr(:,:,:,:,1)
  END SUBROUTINE add_var_list_element_s4d

  SUBROUTINE add_var_list_element_s3d(this_list, varname, ptr, hgrid, vgrid, &
    & cf, grib2, ldims, loutput, lcontainer, lrestart, lrestart_cont,     &
    & initval, isteptype, resetval, lmiss, missval, tlev_source, info,    &
    & p5, vert_interp, hor_interp, in_group, verbose, new_element, tracer_info, &
    & l_pp_scheduler_task, post_op, action_list, var_class, lopenacc)
    TYPE(t_var_list_ptr), INTENT(INOUT) :: this_list
    CHARACTER(len=*), INTENT(IN) :: varname
    REAL(sp), POINTER, INTENT(OUT) :: ptr(:,:,:)
    INTEGER, INTENT(IN) :: hgrid, vgrid, ldims(:)
    TYPE(t_cf_var), INTENT(IN) :: cf
    TYPE(t_grib2_var), INTENT(IN) :: grib2
    LOGICAL, INTENT(IN), OPTIONAL :: loutput, lcontainer, lrestart, &
      & lrestart_cont, lmiss, in_group(:), verbose, lopenacc
    INTEGER, INTENT(IN), OPTIONAL :: isteptype, l_pp_scheduler_task, &
      & tlev_source, var_class
    REAL(sp), INTENT(IN), OPTIONAL :: initval, resetval, missval
    TYPE(t_var_metadata), POINTER, OPTIONAL :: info
    CLASS(t_tracer_meta),    INTENT(in), OPTIONAL :: tracer_info
    REAL(sp), CONTIGUOUS_TARGET, OPTIONAL :: p5(:,:,:,:,:)
    TYPE(t_vert_interp_meta), INTENT(IN), OPTIONAL :: vert_interp
    TYPE(t_hor_interp_meta), INTENT(IN), OPTIONAL :: hor_interp
    TYPE(t_list_element), POINTER, OPTIONAL :: new_element
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op
    TYPE(t_var_action), INTENT(IN), OPTIONAL :: action_list
    TYPE(t_list_element), POINTER :: element

    CALL add_var_list_element_5d(SINGLE_T, this_list, varname,             &
      &   hgrid, vgrid, cf, grib2, ldims, element, loutput, lcontainer,        &
      &   lrestart, lrestart_cont, isteptype, lmiss, tlev_source,              &
      &   info, vert_interp, hor_interp, in_group, verbose,                    &
      &   l_pp_scheduler_task, post_op, action_list, tracer_info=tracer_info,  &
      &   p5_s=p5, initval_s=initval, resetval_s=resetval, missval_s=missval,  &
      &   var_class=var_class, lopenacc=lopenacc, new_element=new_element)
    ptr => element%field%s_ptr(:,:,:,1,1)
  END SUBROUTINE add_var_list_element_s3d

  SUBROUTINE add_var_list_element_s2d(this_list, varname, ptr, hgrid, vgrid, &
    & cf, grib2, ldims, loutput, lcontainer, lrestart, lrestart_cont,     &
    & initval, isteptype, resetval, lmiss, missval, tlev_source, info,    &
    & p5, vert_interp, hor_interp, in_group, verbose, new_element, tracer_info, &
    & l_pp_scheduler_task, post_op, action_list, var_class, lopenacc)
    TYPE(t_var_list_ptr), INTENT(INOUT) :: this_list
    CHARACTER(len=*), INTENT(IN) :: varname
    REAL(sp), POINTER, INTENT(OUT) :: ptr(:,:)
    INTEGER, INTENT(IN) :: hgrid, vgrid, ldims(2)
    TYPE(t_cf_var), INTENT(IN) :: cf
    TYPE(t_grib2_var), INTENT(IN) :: grib2
    LOGICAL, INTENT(IN), OPTIONAL :: loutput, lcontainer, lrestart, &
      & lrestart_cont, lmiss, in_group(:), verbose, lopenacc
    INTEGER, INTENT(IN), OPTIONAL :: isteptype, l_pp_scheduler_task, &
      & tlev_source, var_class
    REAL(sp), INTENT(IN), OPTIONAL :: initval, resetval, missval
    TYPE(t_var_metadata), POINTER, OPTIONAL :: info
    CLASS(t_tracer_meta),    INTENT(in), OPTIONAL :: tracer_info
    REAL(sp), CONTIGUOUS_TARGET, OPTIONAL :: p5(:,:,:,:,:)
    TYPE(t_vert_interp_meta), INTENT(IN), OPTIONAL :: vert_interp
    TYPE(t_hor_interp_meta), INTENT(IN), OPTIONAL :: hor_interp
    TYPE(t_list_element), POINTER, OPTIONAL :: new_element
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op
    TYPE(t_var_action), INTENT(IN), OPTIONAL :: action_list
    TYPE(t_list_element), POINTER :: element

    CALL add_var_list_element_5d(SINGLE_T, this_list, varname,             &
      &   hgrid, vgrid, cf, grib2, ldims, element, loutput, lcontainer,        &
      &   lrestart, lrestart_cont, isteptype, lmiss, tlev_source,              &
      &   info, vert_interp, hor_interp, in_group, verbose,                    &
      &   l_pp_scheduler_task, post_op, action_list, tracer_info=tracer_info,  &
      &   p5_s=p5, initval_s=initval, resetval_s=resetval, missval_s=missval,  &
      &   var_class=var_class, lopenacc=lopenacc, new_element=new_element)
    ptr => element%field%s_ptr(:,:,1,1,1)
  END SUBROUTINE add_var_list_element_s2d

  SUBROUTINE add_var_list_element_s1d(this_list, varname, ptr, hgrid, vgrid, &
    & cf, grib2, ldims, loutput, lcontainer, lrestart, lrestart_cont,     &
    & initval, isteptype, resetval, lmiss, missval, tlev_source, info,    &
    & p5, vert_interp, hor_interp, in_group, verbose, new_element,        &
    & l_pp_scheduler_task, post_op, action_list, var_class, lopenacc)
    TYPE(t_var_list_ptr), INTENT(INOUT) :: this_list
    CHARACTER(len=*), INTENT(IN) :: varname
    REAL(sp), POINTER, INTENT(OUT) :: ptr(:)
    INTEGER, INTENT(IN) :: hgrid, vgrid, ldims(1)
    TYPE(t_cf_var), INTENT(IN) :: cf
    TYPE(t_grib2_var), INTENT(IN) :: grib2
    LOGICAL, INTENT(IN), OPTIONAL :: loutput, lcontainer, lrestart, &
      & lrestart_cont, lmiss, in_group(:), verbose, lopenacc
    INTEGER, INTENT(IN), OPTIONAL :: isteptype, l_pp_scheduler_task, &
      & tlev_source, var_class
    REAL(sp), INTENT(IN), OPTIONAL :: initval, resetval, missval
    TYPE(t_var_metadata), POINTER, OPTIONAL :: info
    REAL(sp), CONTIGUOUS_TARGET, OPTIONAL :: p5(:,:,:,:,:)
    TYPE(t_vert_interp_meta), INTENT(IN), OPTIONAL :: vert_interp
    TYPE(t_hor_interp_meta), INTENT(IN), OPTIONAL :: hor_interp
    TYPE(t_list_element), POINTER, OPTIONAL :: new_element
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op
    TYPE(t_var_action), INTENT(IN), OPTIONAL :: action_list
    TYPE(t_list_element), POINTER :: element

    CALL add_var_list_element_5d(SINGLE_T, this_list, varname,      &
      &   hgrid, vgrid, cf, grib2, ldims, element, loutput, lcontainer, &
      &   lrestart, lrestart_cont, isteptype, lmiss, tlev_source, info, &
      &   vert_interp, hor_interp, in_group, verbose,                   &
      &   l_pp_scheduler_task, post_op, action_list, p5_s=p5,           &
      &   initval_s=initval, resetval_s=resetval, missval_s=missval,    &
      &   var_class=var_class, lopenacc=lopenacc, new_element=new_element)
    ptr => element%field%s_ptr(:,1,1,1,1)
  END SUBROUTINE add_var_list_element_s1d

  SUBROUTINE add_var_list_element_i4d(this_list, varname, ptr, hgrid, vgrid, &
    & cf, grib2, ldims, loutput, lcontainer, lrestart, lrestart_cont,     &
    & initval, isteptype, resetval, lmiss, missval, tlev_source, info,    &
    & p5, vert_interp, hor_interp, in_group, verbose, new_element,        &
    & l_pp_scheduler_task, post_op, action_list, var_class, lopenacc)
    TYPE(t_var_list_ptr), INTENT(INOUT) :: this_list
    CHARACTER(len=*), INTENT(IN) :: varname
    INTEGER, POINTER, INTENT(OUT) :: ptr(:,:,:,:)
    INTEGER, INTENT(IN) :: hgrid, vgrid, ldims(4)
    TYPE(t_cf_var), INTENT(IN) :: cf
    TYPE(t_grib2_var), INTENT(IN) :: grib2
    LOGICAL, INTENT(IN), OPTIONAL :: loutput, lcontainer, lrestart, &
      & lrestart_cont, lmiss, in_group(:), verbose, lopenacc
    INTEGER, INTENT(IN), OPTIONAL :: isteptype, l_pp_scheduler_task, &
      & tlev_source, var_class, initval, resetval, missval
    TYPE(t_var_metadata), POINTER, OPTIONAL :: info
    INTEGER, CONTIGUOUS_TARGET, OPTIONAL :: p5(:,:,:,:,:)
    TYPE(t_vert_interp_meta), INTENT(IN), OPTIONAL :: vert_interp
    TYPE(t_hor_interp_meta), INTENT(IN), OPTIONAL :: hor_interp
    TYPE(t_list_element), POINTER, OPTIONAL :: new_element
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op
    TYPE(t_var_action), INTENT(IN), OPTIONAL :: action_list
    TYPE(t_list_element), POINTER :: element

    CALL add_var_list_element_5d(INT_T, this_list, varname,         &
      &   hgrid, vgrid, cf, grib2, ldims, element, loutput, lcontainer, &
      &   lrestart, lrestart_cont, isteptype, lmiss, tlev_source, info, &
      &   vert_interp, hor_interp, in_group, verbose,                   &
      &   l_pp_scheduler_task, post_op, action_list, p5_i=p5,           &
      &   initval_i=initval, resetval_i=resetval, missval_i=missval,    &
      &   var_class=var_class, lopenacc=lopenacc, new_element=new_element)
    ptr => element%field%i_ptr(:,:,:,:,1)
  END SUBROUTINE add_var_list_element_i4d

  SUBROUTINE add_var_list_element_i3d(this_list, varname, ptr, hgrid, vgrid, &
    & cf, grib2, ldims, loutput, lcontainer, lrestart, lrestart_cont,     &
    & initval, isteptype, resetval, lmiss, missval, tlev_source, info,    &
    & p5, vert_interp, hor_interp, in_group, verbose, new_element,        &
    & l_pp_scheduler_task, post_op, action_list, var_class, lopenacc)
    TYPE(t_var_list_ptr), INTENT(INOUT) :: this_list
    CHARACTER(len=*), INTENT(IN) :: varname
    INTEGER, POINTER, INTENT(OUT) :: ptr(:,:,:)
    INTEGER, INTENT(IN) :: hgrid, vgrid, ldims(3)
    TYPE(t_cf_var), INTENT(IN) :: cf
    TYPE(t_grib2_var), INTENT(IN) :: grib2
    LOGICAL, INTENT(IN), OPTIONAL :: loutput, lcontainer, lrestart, &
      & lrestart_cont, lmiss, in_group(:), verbose, lopenacc
    INTEGER, INTENT(IN), OPTIONAL :: isteptype, l_pp_scheduler_task, &
      & tlev_source, var_class, initval, resetval, missval
    TYPE(t_var_metadata), POINTER, OPTIONAL :: info
    INTEGER, CONTIGUOUS_TARGET, OPTIONAL :: p5(:,:,:,:,:)
    TYPE(t_vert_interp_meta), INTENT(IN), OPTIONAL :: vert_interp
    TYPE(t_hor_interp_meta), INTENT(IN), OPTIONAL :: hor_interp
    TYPE(t_list_element), POINTER, OPTIONAL :: new_element
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op
    TYPE(t_var_action), INTENT(IN), OPTIONAL :: action_list
    TYPE(t_list_element), POINTER :: element

    CALL add_var_list_element_5d(INT_T, this_list, varname,         &
      &   hgrid, vgrid, cf, grib2, ldims, element, loutput, lcontainer, &
      &   lrestart, lrestart_cont, isteptype, lmiss, tlev_source, info, &
      &   vert_interp, hor_interp, in_group, verbose,                   &
      &   l_pp_scheduler_task, post_op, action_list, p5_i=p5,           &
      &   initval_i=initval, resetval_i=resetval, missval_i=missval,    &
      &   var_class=var_class, lopenacc=lopenacc, new_element=new_element)
    ptr => element%field%i_ptr(:,:,:,1,1)
  END SUBROUTINE add_var_list_element_i3d

  SUBROUTINE add_var_list_element_i2d(this_list, varname, ptr, hgrid, vgrid, &
    & cf, grib2, ldims, loutput, lcontainer, lrestart, lrestart_cont,     &
    & initval, isteptype, resetval, lmiss, missval, tlev_source, info,    &
    & p5, vert_interp, hor_interp, in_group, verbose, new_element,        &
    & l_pp_scheduler_task, post_op, action_list, var_class, lopenacc)
    TYPE(t_var_list_ptr), INTENT(INOUT) :: this_list
    CHARACTER(len=*), INTENT(IN) :: varname
    INTEGER, POINTER, INTENT(OUT) :: ptr(:,:)
    INTEGER, INTENT(IN) :: hgrid, vgrid, ldims(2)
    TYPE(t_cf_var), INTENT(IN) :: cf
    TYPE(t_grib2_var), INTENT(IN) :: grib2
    LOGICAL, INTENT(IN), OPTIONAL :: loutput, lcontainer, lrestart, &
      & lrestart_cont, lmiss, in_group(:), verbose, lopenacc
    INTEGER, INTENT(IN), OPTIONAL :: isteptype, l_pp_scheduler_task, &
      & tlev_source, var_class, initval, resetval, missval
    TYPE(t_var_metadata), POINTER, OPTIONAL :: info
    INTEGER, CONTIGUOUS_TARGET, OPTIONAL :: p5(:,:,:,:,:)
    TYPE(t_vert_interp_meta), INTENT(IN), OPTIONAL :: vert_interp
    TYPE(t_hor_interp_meta), INTENT(IN), OPTIONAL :: hor_interp
    TYPE(t_list_element), POINTER, OPTIONAL :: new_element
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op
    TYPE(t_var_action), INTENT(IN), OPTIONAL :: action_list
    TYPE(t_list_element), POINTER :: element

    CALL add_var_list_element_5d(INT_T, this_list, varname,         &
      &   hgrid, vgrid, cf, grib2, ldims, element, loutput, lcontainer, &
      &   lrestart, lrestart_cont, isteptype, lmiss, tlev_source, info, &
      &   vert_interp, hor_interp, in_group, verbose,                   &
      &   l_pp_scheduler_task, post_op, action_list, p5_i=p5,           &
      &   initval_i=initval, resetval_i=resetval, missval_i=missval,    &
      &   var_class=var_class, lopenacc=lopenacc, new_element=new_element)
    ptr => element%field%i_ptr(:,:,1,1,1)
  END SUBROUTINE add_var_list_element_i2d

  SUBROUTINE add_var_list_element_i1d(this_list, varname, ptr, hgrid, vgrid, &
    & cf, grib2, ldims, loutput, lcontainer, lrestart, lrestart_cont,     &
    & initval, isteptype, resetval, lmiss, missval, tlev_source, info,    &
    & p5, vert_interp, hor_interp, in_group, verbose, new_element,        &
    & l_pp_scheduler_task, post_op, action_list, var_class, lopenacc)
    TYPE(t_var_list_ptr), INTENT(INOUT) :: this_list
    CHARACTER(len=*), INTENT(IN) :: varname
    INTEGER, POINTER, INTENT(OUT) :: ptr(:)
    INTEGER, INTENT(IN) :: hgrid, vgrid, ldims(1)
    TYPE(t_cf_var), INTENT(IN) :: cf
    TYPE(t_grib2_var), INTENT(IN) :: grib2
    LOGICAL, INTENT(IN), OPTIONAL :: loutput, lcontainer, lrestart, &
      & lrestart_cont, lmiss, in_group(:), verbose, lopenacc
    INTEGER, INTENT(IN), OPTIONAL :: isteptype, l_pp_scheduler_task, &
      & tlev_source, var_class, initval, resetval, missval
    TYPE(t_var_metadata), POINTER, OPTIONAL :: info
    INTEGER, CONTIGUOUS_TARGET, OPTIONAL :: p5(:,:,:,:,:)
    TYPE(t_vert_interp_meta), INTENT(IN), OPTIONAL :: vert_interp
    TYPE(t_hor_interp_meta), INTENT(IN), OPTIONAL :: hor_interp
    TYPE(t_list_element), POINTER, OPTIONAL :: new_element
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op
    TYPE(t_var_action), INTENT(IN), OPTIONAL :: action_list
    TYPE(t_list_element), POINTER :: element

    CALL add_var_list_element_5d(INT_T, this_list, varname,         &
      &   hgrid, vgrid, cf, grib2, ldims, element, loutput, lcontainer, &
      &   lrestart, lrestart_cont, isteptype, lmiss, tlev_source, info, &
      &   vert_interp, hor_interp, in_group, verbose,                   &
      &   l_pp_scheduler_task, post_op, action_list, p5_i=p5,           &
      &   initval_i=initval, resetval_i=resetval, missval_i=missval,    &
      &   var_class=var_class, lopenacc=lopenacc, new_element=new_element)
    ptr => element%field%i_ptr(:,1,1,1,1)
  END SUBROUTINE add_var_list_element_i1d

  SUBROUTINE add_var_list_element_l4d(this_list, varname, ptr, hgrid, vgrid, &
    & cf, grib2, ldims, loutput, lcontainer, lrestart, lrestart_cont,     &
    & initval, isteptype, resetval, lmiss, missval, tlev_source, info,    &
    & p5, vert_interp, hor_interp, in_group, verbose, new_element,        &
    & l_pp_scheduler_task, post_op, action_list, var_class, lopenacc)
    TYPE(t_var_list_ptr), INTENT(INOUT) :: this_list
    CHARACTER(len=*), INTENT(IN) :: varname
    LOGICAL, POINTER, INTENT(OUT) :: ptr(:,:,:,:)
    INTEGER, INTENT(IN) :: hgrid, vgrid, ldims(4)
    TYPE(t_cf_var), INTENT(IN) :: cf
    TYPE(t_grib2_var), INTENT(IN) :: grib2
    LOGICAL, INTENT(IN), OPTIONAL :: loutput, lcontainer, lrestart, &
      & lrestart_cont, initval, resetval, lmiss, missval, in_group(:), &
      & verbose, lopenacc
    INTEGER, INTENT(IN), OPTIONAL :: isteptype, l_pp_scheduler_task, &
      & tlev_source, var_class
    TYPE(t_var_metadata), POINTER, OPTIONAL :: info
    LOGICAL, CONTIGUOUS_TARGET, OPTIONAL :: p5(:,:,:,:,:)
    TYPE(t_vert_interp_meta), INTENT(IN), OPTIONAL :: vert_interp
    TYPE(t_hor_interp_meta), INTENT(IN), OPTIONAL :: hor_interp
    TYPE(t_list_element), POINTER, OPTIONAL :: new_element
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op
    TYPE(t_var_action), INTENT(IN), OPTIONAL :: action_list
    TYPE(t_list_element), POINTER :: element

    CALL add_var_list_element_5d(BOOL_T, this_list, varname,        &
      &   hgrid, vgrid, cf, grib2, ldims, element, loutput, lcontainer, &
      &   lrestart, lrestart_cont, isteptype, lmiss, tlev_source, info, &
      &   vert_interp, hor_interp, in_group, verbose,                   &
      &   l_pp_scheduler_task, post_op, action_list, p5_l=p5,           &
      &   initval_l=initval, resetval_l=resetval, missval_l=missval,    &
      &   var_class=var_class, lopenacc=lopenacc, new_element=new_element)
    ptr => element%field%l_ptr(:,:,:,:,1)
  END SUBROUTINE add_var_list_element_l4d

  SUBROUTINE add_var_list_element_l3d(this_list, varname, ptr, hgrid, vgrid, &
    & cf, grib2, ldims, loutput, lcontainer, lrestart, lrestart_cont,     &
    & initval, isteptype, resetval, lmiss, missval, tlev_source, info,    &
    & p5, vert_interp, hor_interp, in_group, verbose, new_element,        &
    & l_pp_scheduler_task, post_op, action_list, var_class, lopenacc)
    TYPE(t_var_list_ptr), INTENT(INOUT) :: this_list
    CHARACTER(len=*), INTENT(IN) :: varname
    LOGICAL, POINTER, INTENT(OUT) :: ptr(:,:,:)
    INTEGER, INTENT(IN) :: hgrid, vgrid, ldims(3)
    TYPE(t_cf_var), INTENT(IN) :: cf
    TYPE(t_grib2_var), INTENT(IN) :: grib2
    LOGICAL, INTENT(IN), OPTIONAL :: loutput, lcontainer, lrestart, &
      & lrestart_cont, initval, resetval, lmiss, missval, in_group(:), &
      & verbose, lopenacc
    INTEGER, INTENT(IN), OPTIONAL :: isteptype, l_pp_scheduler_task, &
      & tlev_source, var_class
    TYPE(t_var_metadata), POINTER, OPTIONAL :: info
    LOGICAL, CONTIGUOUS_TARGET, OPTIONAL :: p5(:,:,:,:,:)
    TYPE(t_vert_interp_meta), INTENT(IN), OPTIONAL :: vert_interp
    TYPE(t_hor_interp_meta), INTENT(IN), OPTIONAL :: hor_interp
    TYPE(t_list_element), POINTER, OPTIONAL :: new_element
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op
    TYPE(t_var_action), INTENT(IN), OPTIONAL :: action_list
    TYPE(t_list_element), POINTER :: element

    CALL add_var_list_element_5d(BOOL_T, this_list, varname,        &
      &   hgrid, vgrid, cf, grib2, ldims, element, loutput, lcontainer, &
      &   lrestart, lrestart_cont, isteptype, lmiss, tlev_source, info, &
      &   vert_interp, hor_interp, in_group, verbose,                   &
      &   l_pp_scheduler_task, post_op, action_list, p5_l=p5,           &
      &   initval_l=initval, resetval_l=resetval, missval_l=missval,    &
      &   var_class=var_class, lopenacc=lopenacc, new_element=new_element)
    ptr => element%field%l_ptr(:,:,:,1,1)
  END SUBROUTINE add_var_list_element_l3d

  SUBROUTINE add_var_list_element_l2d(this_list, varname, ptr, hgrid, vgrid, &
    & cf, grib2, ldims, loutput, lcontainer, lrestart, lrestart_cont,     &
    & initval, isteptype, resetval, lmiss, missval, tlev_source, info,    &
    & p5, vert_interp, hor_interp, in_group, verbose, new_element,        &
    & l_pp_scheduler_task, post_op, action_list, var_class, lopenacc)
    TYPE(t_var_list_ptr), INTENT(INOUT) :: this_list
    CHARACTER(len=*), INTENT(IN) :: varname
    LOGICAL, POINTER, INTENT(OUT) :: ptr(:,:)
    INTEGER, INTENT(IN) :: hgrid, vgrid, ldims(2)
    TYPE(t_cf_var), INTENT(IN) :: cf
    TYPE(t_grib2_var), INTENT(IN) :: grib2
    LOGICAL, INTENT(IN), OPTIONAL :: loutput, lcontainer, lrestart, &
      & lrestart_cont, initval, resetval, lmiss, missval, in_group(:), &
      & verbose, lopenacc
    INTEGER, INTENT(IN), OPTIONAL :: isteptype, l_pp_scheduler_task, &
      & tlev_source, var_class
    TYPE(t_var_metadata), POINTER, OPTIONAL :: info
    LOGICAL, CONTIGUOUS_TARGET, OPTIONAL :: p5(:,:,:,:,:)
    TYPE(t_vert_interp_meta), INTENT(IN), OPTIONAL :: vert_interp
    TYPE(t_hor_interp_meta), INTENT(IN), OPTIONAL :: hor_interp
    TYPE(t_list_element), POINTER, OPTIONAL :: new_element
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op
    TYPE(t_var_action), INTENT(IN), OPTIONAL :: action_list
    TYPE(t_list_element), POINTER :: element

    CALL add_var_list_element_5d(BOOL_T, this_list, varname,        &
      &   hgrid, vgrid, cf, grib2, ldims, element, loutput, lcontainer, &
      &   lrestart, lrestart_cont, isteptype, lmiss, tlev_source, info, &
      &   vert_interp, hor_interp, in_group, verbose,                   &
      &   l_pp_scheduler_task, post_op, action_list, p5_l=p5,           &
      &   initval_l=initval, resetval_l=resetval, missval_l=missval,    &
      &   var_class=var_class, lopenacc=lopenacc, new_element=new_element)
    ptr => element%field%l_ptr(:,:,1,1,1)
  END SUBROUTINE add_var_list_element_l2d

  SUBROUTINE add_var_list_element_l1d(this_list, varname, ptr, hgrid, vgrid, &
    & cf, grib2, ldims, loutput, lcontainer, lrestart, lrestart_cont,     &
    & initval, isteptype, resetval, lmiss, missval, tlev_source, info,    &
    & p5, vert_interp, hor_interp, in_group, verbose, new_element,        &
    & l_pp_scheduler_task, post_op, action_list, var_class, lopenacc)
    TYPE(t_var_list_ptr), INTENT(INOUT) :: this_list
    CHARACTER(len=*), INTENT(IN) :: varname
    LOGICAL, POINTER, INTENT(OUT) :: ptr(:)
    INTEGER, INTENT(IN) :: hgrid, vgrid, ldims(1)
    TYPE(t_cf_var), INTENT(IN) :: cf
    TYPE(t_grib2_var), INTENT(IN) :: grib2
    LOGICAL, INTENT(IN), OPTIONAL :: loutput, lcontainer, lrestart, &
      & lrestart_cont, initval, resetval, lmiss, missval, in_group(:), &
      & verbose, lopenacc
    INTEGER, INTENT(IN), OPTIONAL :: isteptype, l_pp_scheduler_task, &
      & tlev_source, var_class
    TYPE(t_var_metadata), POINTER, OPTIONAL :: info
    LOGICAL, CONTIGUOUS_TARGET, OPTIONAL :: p5(:,:,:,:,:)
    TYPE(t_vert_interp_meta), INTENT(IN), OPTIONAL :: vert_interp
    TYPE(t_hor_interp_meta), INTENT(IN), OPTIONAL :: hor_interp
    TYPE(t_list_element), POINTER, OPTIONAL :: new_element
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op
    TYPE(t_var_action), INTENT(IN), OPTIONAL :: action_list
    TYPE(t_list_element), POINTER :: element

    CALL add_var_list_element_5d(BOOL_T, this_list, varname,        &
      &   hgrid, vgrid, cf, grib2, ldims, element, loutput, lcontainer, &
      &   lrestart, lrestart_cont, isteptype, lmiss, tlev_source, info, &
      &   vert_interp, hor_interp, in_group, verbose,                   &
      &   l_pp_scheduler_task, post_op, action_list, p5_l=p5,           &
      &   initval_l=initval, resetval_l=resetval, missval_l=missval,    &
      &   var_class=var_class, lopenacc=lopenacc, new_element=new_element)
    ptr => element%field%l_ptr(:,1,1,1,1)
  END SUBROUTINE add_var_list_element_l1d

  SUBROUTINE add_var_list_reference_util(target_element, new_list_element,      &
    & this_list, target_name, refname, hgrid, vgrid, cf, grib2, ref_idx, ldims, &
    & dtype, icontainer, loutput, lrestart, lrestart_cont, isteptype, lmiss,    &
    & tlev_source, tracer_info, info, vert_interp, hor_interp, in_group,        &
    & verbose, new_element, l_pp_scheduler_task, post_op, action_list, idx_diag,&
    & var_class, opt_var_ref_pos, initval_r, initval_s, initval_i, missval_r,   &
    & missval_s, missval_i, resetval_r, resetval_s, resetval_i, idx_tracer)
    TYPE(t_var_list_ptr), INTENT(INOUT) :: this_list
    TYPE(t_list_element), INTENT(OUT), POINTER :: target_element, new_list_element
    CHARACTER(len=*), INTENT(IN) :: target_name, refname
    INTEGER, INTENT(IN) :: hgrid, vgrid, ref_idx, ldims(:), dtype
    INTEGER, INTENT(OUT) :: icontainer
    TYPE(t_cf_var), INTENT(IN) :: cf
    TYPE(t_grib2_var), INTENT(IN) :: grib2
    LOGICAL, INTENT(IN), OPTIONAL :: loutput, lrestart, lrestart_cont, &
      & lmiss, in_group(:), verbose
    INTEGER, INTENT(IN), OPTIONAL :: isteptype, tlev_source, var_class, &
      & l_pp_scheduler_task, opt_var_ref_pos, idx_tracer, idx_diag
    REAL(dp), INTENT(IN), OPTIONAL :: initval_r, resetval_r, missval_r
    REAL(sp), INTENT(IN), OPTIONAL :: initval_s, resetval_s, missval_s
    INTEGER, INTENT(IN), OPTIONAL :: initval_i, resetval_i, missval_i
    CLASS(t_tracer_meta), INTENT(IN), OPTIONAL :: tracer_info
    TYPE(t_var_metadata), POINTER, OPTIONAL :: info
    TYPE(t_vert_interp_meta),INTENT(IN), OPTIONAL :: vert_interp
    TYPE(t_hor_interp_meta), INTENT(IN), OPTIONAL :: hor_interp
    TYPE(t_list_element), POINTER, OPTIONAL :: new_element
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op
    TYPE(t_var_action), INTENT(IN), OPTIONAL :: action_list
    TYPE(t_var_metadata), POINTER :: target_info, ref_info
    TYPE(t_union_vals) :: missvalt, initvalt, resetvalt
    CHARACTER(*), PARAMETER :: routine = modname//":add_var_list_reference_util"
    INTEGER :: var_ref_pos, ndims, di(5), di3, max_ref
    LOGICAL :: in_group_new(MAX_GROUPS)

    ndims = SIZE(ldims)
    target_element => find_list_element(this_list, target_name)
    target_info => target_element%field%info
    IF (PRESENT(opt_var_ref_pos)) THEN
      var_ref_pos = opt_var_ref_pos
      IF (.NOT. target_info%lcontainer) &
        &  CALL finish(routine, "invalid container index: Target is not a container variable!")
      IF ((target_info%var_ref_pos /= var_ref_pos) .AND. &
        & (target_info%var_ref_pos /= -1)) THEN
        CALL finish(routine, "Container index does not match the previously set value!")
      END IF
      target_info%var_ref_pos = var_ref_pos
    ELSE
      var_ref_pos = ndims + 1
    END IF
    di3 = MERGE(4, 0, ndims.EQ.3)
    IF (.NOT. ANY(NDIMS .EQ. (/ 2, 3 /))) CALL finish(routine, "Internal error!")
    SELECT CASE(var_ref_pos)
    CASE(1)
      di = (/ 2, 3, di3, 0, 0 /)
    CASE(2)
      di = (/ 1, 3, di3, 0, 0 /)
    CASE(3)
      di = (/ 1, 2, di3, 0, 0 /)
    CASE(4)
      IF (NDIMS.EQ.2) CALL finish(routine, "Internal error!")
      di = (/ 1, 2, 3, 0, 0 /)
    CASE DEFAULT
      CALL finish(routine, "Internal error!")
    END SELECT
    IF (target_info%lcontainer) THEN
      max_ref = 0
      IF (ASSOCIATED(target_element%field%r_ptr)) THEN
        max_ref = SIZE(target_element%field%r_ptr, var_ref_pos)
      ELSE IF (ASSOCIATED(target_element%field%s_ptr)) THEN
        max_ref = SIZE(target_element%field%s_ptr, var_ref_pos)
      ELSE IF (ASSOCIATED(target_element%field%i_ptr)) THEN
        max_ref = SIZE(target_element%field%i_ptr, var_ref_pos)
      ELSE IF (ASSOCIATED(target_element%field%l_ptr)) THEN
        max_ref = SIZE(target_element%field%l_ptr, var_ref_pos)
      END IF
      ! Counting the number of existing references is deactivated, 
      ! if the slice index to be referenced is given explicitly.
        target_info%ncontained = target_info%ncontained+1
        ! only check validity of given slice index
        IF ( (ref_idx > max_ref) .OR. (ref_idx < 1)) THEN
          WRITE (message_text, "(2(a,i3),a)") 'Slice idx ', ref_idx, ' for ' // &
            & TRIM(refname) // ' out of allowable range [1,',max_ref,']'
          CALL finish(routine, message_text)
        ENDIF
      IF (ANY(ldims(1:ndims) /= target_info%used_dimensions(di(1:ndims)))) &
        & CALL finish(routine, TRIM(refname)//' dimensions requested and available differ.')
    ENDIF
    ! add list entry
    CALL append_list_element (this_list, new_list_element)
    IF (PRESENT(new_element)) new_element=>new_list_element
    new_list_element%field%ref_to => target_element%field
    ref_info => new_list_element%field%info
    CALL inherit_var_list_metadata(ref_info, this_list)
    ! init local fields
    missvalt  = ref_info%missval
    initvalt  = ref_info%initval
    resetvalt = ref_info%resetval
    IF (PRESENT(missval_r))  missvalt%rval  = missval_r
    IF (PRESENT(missval_s))  missvalt%sval  = missval_s
    IF (PRESENT(missval_i))  missvalt%ival  = missval_i
    IF (PRESENT(initval_r))  initvalt%rval  = initval_r
    IF (PRESENT(initval_s))  initvalt%sval  = initval_s
    IF (PRESENT(initval_i))  initvalt%ival  = initval_i
    IF (PRESENT(resetval_r)) resetvalt%rval = resetval_r
    IF (PRESENT(resetval_s)) resetvalt%sval = resetval_s
    IF (PRESENT(resetval_i)) resetvalt%ival = resetval_i
    CALL set_var_metadata (ref_info, name=refname, hgrid=hgrid, vgrid=vgrid, &
      & cf=cf, grib2=grib2, ldims=ldims, loutput=loutput, lrestart=lrestart, &
      & lrestart_cont=lrestart_cont, initval=initvalt, isteptype=isteptype,  &
      & resetval=resetvalt, lmiss=lmiss, missval=missvalt, verbose=verbose,  &
      & tlev_source=tlev_source, vert_interp=vert_interp, in_group=in_group, &
      & hor_interp=hor_interp, l_pp_scheduler_task=l_pp_scheduler_task,      &
      & post_op=post_op, action_list=action_list, var_class=var_class,       &
      & idx_tracer=idx_tracer, idx_diag=idx_diag, data_type=dtype)
    ! set dynamic metadata, i.e. polymorphic tracer metadata
    CALL set_var_metadata_dyn(new_list_element%field%info_dyn, tracer_info=tracer_info)
    ref_info%ndims = ndims
    ref_info%used_dimensions(:) = 0
    ref_info%used_dimensions(1:ndims) = target_info%used_dimensions(di(1:ndims))
    IF (PRESENT(var_class)) THEN
      IF ( ANY((/CLASS_TILE, CLASS_TILE_LAND/) == var_class)) THEN
        ! automatically add tile to its variable specific tile-group
        CALL var_groups_dyn%add(group_name=target_info%name, &
          & in_group_new=in_group_new, opt_in_group=in_group)
        ! update in_group metainfo
        ref_info%in_group(:) = in_group_new(:)
      ENDIF
    END IF
    IF (target_info%lcontainer) THEN
      ref_info%lcontained                   = .TRUE.
      ref_info%used_dimensions(ndims+1)     = 1
      ref_info%var_ref_pos = var_ref_pos
      ref_info%maxcontained = max_ref
      ref_info%ncontained = ref_idx
    ENDIF
    icontainer = MERGE(ref_info%ncontained, 1, target_info%lcontainer)
    IF(PRESENT(info)) info => ref_info
  END SUBROUTINE add_var_list_reference_util

  SUBROUTINE add_var_list_reference_r3d(this_list, target_name, refname, ptr,       &
    & hgrid, vgrid, cf, grib2, ref_idx, ldims, loutput, lrestart, lrestart_cont, &
    & initval, isteptype, resetval, lmiss, missval, tlev_source, tracer_info,    &
    & info, vert_interp, hor_interp, in_group, verbose, new_element,             &
    & l_pp_scheduler_task, post_op, action_list, opt_var_ref_pos, var_class)
    TYPE(t_var_list_ptr), INTENT(INOUT) :: this_list
    CHARACTER(len=*), INTENT(IN) :: target_name, refname
    REAL(dp), POINTER :: ptr(:,:,:)
    INTEGER, INTENT(IN) :: hgrid, vgrid, ref_idx, ldims(3)
    TYPE(t_cf_var), INTENT(IN) :: cf
    TYPE(t_grib2_var), INTENT(IN) :: grib2
    LOGICAL, INTENT(IN), OPTIONAL :: loutput, lrestart, lrestart_cont, &
      & lmiss, in_group(:), verbose
    INTEGER, INTENT(IN), OPTIONAL :: isteptype, tlev_source, var_class, &
      & l_pp_scheduler_task, opt_var_ref_pos
    REAL(dp), INTENT(IN), OPTIONAL :: initval, resetval, missval
    CLASS(t_tracer_meta), INTENT(IN), OPTIONAL :: tracer_info
    TYPE(t_var_metadata), POINTER, OPTIONAL :: info
    TYPE(t_vert_interp_meta),INTENT(IN), OPTIONAL :: vert_interp
    TYPE(t_hor_interp_meta), INTENT(IN), OPTIONAL :: hor_interp
    TYPE(t_list_element), POINTER, OPTIONAL :: new_element
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op
    TYPE(t_var_action), INTENT(IN), OPTIONAL :: action_list
    CHARACTER(*), PARAMETER :: routine = modname//"::add_var_list_reference_r3d"
    TYPE(t_list_element), POINTER :: target_element, new_list_element
    INTEGER :: icontainer

    CALL add_var_list_reference_util(target_element, new_list_element,     &
      & this_list, target_name, refname, hgrid, vgrid, cf, grib2, ref_idx, &
      & ldims, REAL_T, icontainer, loutput=loutput, lrestart=lrestart,     &
      & lrestart_cont=lrestart_cont, isteptype=isteptype, lmiss=lmiss,     &
      & tlev_source=tlev_source, tracer_info=tracer_info, info=info,       &
      & vert_interp=vert_interp, hor_interp=hor_interp, in_group=in_group, &
      & verbose=verbose, new_element=new_element,                          &
      & l_pp_scheduler_task=l_pp_scheduler_task, post_op=post_op,          &
      & action_list=action_list, var_class=var_class,                      &
      & opt_var_ref_pos=opt_var_ref_pos, initval_r=initval,                &
      & missval_r=missval, resetval_r=resetval)
    IF (.NOT. ASSOCIATED(target_element%field%r_ptr)) &
      & CALL finish(routine, TRIM(refname)//' not created.')
    SELECT CASE(new_list_element%field%info%var_ref_pos)
    CASE(1)
      ptr => target_element%field%r_ptr(icontainer,:,:,:,1)
    CASE(2)
      ptr => target_element%field%r_ptr(:,icontainer,:,:,1)
    CASE(3)
      ptr => target_element%field%r_ptr(:,:,icontainer,:,1)
    CASE(4)
      ptr => target_element%field%r_ptr(:,:,:,icontainer,1)
    CASE default
      CALL finish(routine, "internal error!")
    END SELECT
    new_list_element%field%r_ptr => target_element%field%r_ptr
    IF (.NOT. ASSOCIATED(new_list_element%field%r_ptr)) &
      & WRITE (0,*) 'problem with association of ptr for '//TRIM(refname)
    IF (PRESENT(lmiss)) THEN
      ptr = new_list_element%field%info%missval%rval
    ELSE
      ptr = 0.0_dp
    END IF
  END SUBROUTINE add_var_list_reference_r3d

  SUBROUTINE add_var_list_reference_r2d(this_list, target_name, refname, ptr,    &
    & hgrid, vgrid, cf, grib2, ref_idx, ldims, loutput, lrestart, lrestart_cont, &
    & initval, isteptype, resetval, lmiss, missval, tlev_source, tracer_info,    &
    & info, vert_interp, hor_interp, in_group, verbose, new_element,             &
    & l_pp_scheduler_task, post_op, action_list, opt_var_ref_pos, var_class,     &
    & idx_tracer, idx_diag)
    TYPE(t_var_list_ptr), INTENT(INOUT) :: this_list
    CHARACTER(len=*), INTENT(IN) :: target_name, refname
    REAL(dp), POINTER :: ptr(:,:)
    INTEGER, INTENT(IN) :: hgrid, vgrid, ref_idx, ldims(2)
    TYPE(t_cf_var), INTENT(IN) :: cf                  
    TYPE(t_grib2_var), INTENT(IN) :: grib2               
    LOGICAL, INTENT(IN), OPTIONAL :: loutput, lrestart, lrestart_cont, &
      & lmiss, in_group(:), verbose
    INTEGER, INTENT(IN), OPTIONAL :: isteptype, tlev_source, var_class, &
      & l_pp_scheduler_task, opt_var_ref_pos, idx_tracer, idx_diag
    REAL(dp), INTENT(IN), OPTIONAL :: initval, resetval, missval
    CLASS(t_tracer_meta), INTENT(IN), OPTIONAL :: tracer_info         
    TYPE(t_var_metadata), POINTER, OPTIONAL :: info                
    TYPE(t_vert_interp_meta),INTENT(IN), OPTIONAL :: vert_interp         
    TYPE(t_hor_interp_meta), INTENT(IN), OPTIONAL :: hor_interp          
    TYPE(t_list_element), POINTER, OPTIONAL :: new_element         
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op            
    TYPE(t_var_action), INTENT(IN), OPTIONAL :: action_list         
    CHARACTER(*), PARAMETER :: routine = modname//"::add_var_list_reference_r2d"
    TYPE(t_list_element), POINTER :: target_element, new_list_element
    INTEGER :: icontainer

    CALL add_var_list_reference_util(target_element, new_list_element,     &
      & this_list, target_name, refname, hgrid, vgrid, cf, grib2, ref_idx, &
      & ldims, REAL_T, icontainer, loutput=loutput, lrestart=lrestart,     &
      & lrestart_cont=lrestart_cont, isteptype=isteptype, lmiss=lmiss,     &
      & tlev_source=tlev_source, tracer_info=tracer_info, info=info,       &
      & vert_interp=vert_interp, hor_interp=hor_interp, in_group=in_group, &
      & verbose=verbose, new_element=new_element,                          &
      & l_pp_scheduler_task=l_pp_scheduler_task, post_op=post_op,          &
      & action_list=action_list, var_class=var_class, idx_diag=idx_diag,   &
      & opt_var_ref_pos=opt_var_ref_pos, initval_r=initval,                &
      & missval_r=missval, resetval_r=resetval, idx_tracer=idx_tracer)
    IF (.NOT. ASSOCIATED(target_element%field%r_ptr)) &
      & CALL finish(routine, TRIM(refname)//' not created.')
    SELECT CASE(new_list_element%field%info%var_ref_pos)
    CASE(1)
      ptr => target_element%field%r_ptr(icontainer,:,:,1,1)
    CASE(2)
      ptr => target_element%field%r_ptr(:,icontainer,:,1,1)
    CASE(3)
      ptr => target_element%field%r_ptr(:,:,icontainer,1,1)
    CASE default
      CALL finish(routine, "internal error!")
    END SELECT
    new_list_element%field%r_ptr => target_element%field%r_ptr
    IF (.NOT. ASSOCIATED(new_list_element%field%r_ptr)) &
      & WRITE (0,*) 'problem with association of ptr for '//TRIM(refname)
    IF (PRESENT(lmiss)) THEN
      ptr = new_list_element%field%info%missval%rval
    ELSE
      ptr = 0.0_dp
    END IF
  END SUBROUTINE add_var_list_reference_r2d

  SUBROUTINE add_var_list_reference_s3d(this_list, target_name, refname, ptr,    &
    & hgrid, vgrid, cf, grib2, ref_idx, ldims, loutput, lrestart, lrestart_cont, &
    & initval, isteptype, resetval, lmiss, missval, tlev_source, tracer_info,    &
    & info, vert_interp, hor_interp, in_group, verbose, new_element,             &
    & l_pp_scheduler_task, post_op, action_list, opt_var_ref_pos, var_class)
    TYPE(t_var_list_ptr), INTENT(INOUT) :: this_list
    CHARACTER(len=*), INTENT(IN) :: target_name, refname
    REAL(sp), POINTER :: ptr(:,:,:)
    INTEGER, INTENT(IN) :: hgrid, vgrid, ref_idx, ldims(3)
    TYPE(t_cf_var), INTENT(IN) :: cf
    TYPE(t_grib2_var), INTENT(IN) :: grib2
    LOGICAL, INTENT(IN), OPTIONAL :: loutput, lrestart, lrestart_cont, &
      & lmiss, in_group(:), verbose
    INTEGER, INTENT(IN), OPTIONAL :: isteptype, tlev_source, var_class, &
      & l_pp_scheduler_task, opt_var_ref_pos
    REAL(sp), INTENT(IN), OPTIONAL :: initval, resetval, missval
    CLASS(t_tracer_meta), INTENT(IN), OPTIONAL :: tracer_info
    TYPE(t_var_metadata), POINTER, OPTIONAL :: info
    TYPE(t_vert_interp_meta),INTENT(IN), OPTIONAL :: vert_interp
    TYPE(t_hor_interp_meta), INTENT(IN), OPTIONAL :: hor_interp
    TYPE(t_list_element), POINTER, OPTIONAL :: new_element
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op
    TYPE(t_var_action), INTENT(IN), OPTIONAL :: action_list
    CHARACTER(*), PARAMETER :: routine = modname//"::add_var_list_reference_s3d"
    TYPE(t_list_element), POINTER :: target_element, new_list_element
    INTEGER :: icontainer

    CALL add_var_list_reference_util(target_element, new_list_element,     &
      & this_list, target_name, refname, hgrid, vgrid, cf, grib2, ref_idx, &
      & ldims, SINGLE_T, icontainer, loutput=loutput, lrestart=lrestart,   &
      & lrestart_cont=lrestart_cont, isteptype=isteptype, lmiss=lmiss,     &
      & tlev_source=tlev_source, tracer_info=tracer_info, info=info,       &
      & vert_interp=vert_interp, hor_interp=hor_interp, in_group=in_group, &
      & verbose=verbose, new_element=new_element,                          &
      & l_pp_scheduler_task=l_pp_scheduler_task, post_op=post_op,          &
      & action_list=action_list, var_class=var_class,                      &
      & opt_var_ref_pos=opt_var_ref_pos, initval_s=initval,                &
      & missval_s=missval, resetval_s=resetval)
    IF (.NOT. ASSOCIATED(target_element%field%s_ptr)) &
      & CALL finish(routine, TRIM(refname)//' not created.')
    SELECT CASE(new_list_element%field%info%var_ref_pos)
    CASE(1)
      ptr => target_element%field%s_ptr(icontainer,:,:,:,1)
    CASE(2)
      ptr => target_element%field%s_ptr(:,icontainer,:,:,1)
    CASE(3)
      ptr => target_element%field%s_ptr(:,:,icontainer,:,1)
    CASE(4)
      ptr => target_element%field%s_ptr(:,:,:,icontainer,1)
    CASE default
      CALL finish(routine, "internal error!")
    END SELECT
    new_list_element%field%s_ptr => target_element%field%s_ptr
    IF (.NOT. ASSOCIATED(new_list_element%field%s_ptr)) &
      & WRITE (0,*) 'problem with association of ptr for '//TRIM(refname)
    IF (PRESENT(lmiss)) THEN
      ptr = new_list_element%field%info%missval%sval
    ELSE
      ptr = 0.0_sp
    END IF
  END SUBROUTINE add_var_list_reference_s3d

  SUBROUTINE add_var_list_reference_s2d(this_list, target_name, refname, ptr,    &
    & hgrid, vgrid, cf, grib2, ref_idx, ldims, loutput, lrestart, lrestart_cont, &
    & initval, isteptype, resetval, lmiss, missval, tlev_source, tracer_info,    &
    & info, vert_interp, hor_interp, in_group, verbose, new_element,             &
    & l_pp_scheduler_task, post_op, action_list, opt_var_ref_pos, var_class)
    TYPE(t_var_list_ptr), INTENT(INOUT) :: this_list
    CHARACTER(len=*), INTENT(IN) :: target_name, refname
    REAL(sp), POINTER :: ptr(:,:)
    INTEGER, INTENT(IN) :: hgrid, vgrid, ref_idx, ldims(2)
    TYPE(t_cf_var), INTENT(IN) :: cf
    TYPE(t_grib2_var), INTENT(IN) :: grib2
    LOGICAL, INTENT(IN), OPTIONAL :: loutput, lrestart, lrestart_cont, &
      & lmiss, in_group(:), verbose
    INTEGER, INTENT(IN), OPTIONAL :: isteptype, tlev_source, var_class, &
      & l_pp_scheduler_task, opt_var_ref_pos
    REAL(sp), INTENT(IN), OPTIONAL :: initval, resetval, missval
    CLASS(t_tracer_meta), INTENT(IN), OPTIONAL :: tracer_info
    TYPE(t_var_metadata), POINTER, OPTIONAL :: info
    TYPE(t_vert_interp_meta),INTENT(IN), OPTIONAL :: vert_interp
    TYPE(t_hor_interp_meta), INTENT(IN), OPTIONAL :: hor_interp
    TYPE(t_list_element), POINTER, OPTIONAL :: new_element
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op
    TYPE(t_var_action), INTENT(IN), OPTIONAL :: action_list
    CHARACTER(*), PARAMETER :: routine = modname//"::add_var_list_reference_s2d"
    TYPE(t_list_element), POINTER :: target_element, new_list_element
    INTEGER :: icontainer

    CALL add_var_list_reference_util(target_element, new_list_element,     &
      & this_list, target_name, refname, hgrid, vgrid, cf, grib2, ref_idx, &
      & ldims, SINGLE_T, icontainer, loutput=loutput, lrestart=lrestart,   &
      & lrestart_cont=lrestart_cont, isteptype=isteptype, lmiss=lmiss,     &
      & tlev_source=tlev_source, tracer_info=tracer_info, info=info,       &
      & vert_interp=vert_interp, hor_interp=hor_interp, in_group=in_group, &
      & verbose=verbose, new_element=new_element,                          &
      & l_pp_scheduler_task=l_pp_scheduler_task, post_op=post_op,          &
      & action_list=action_list, var_class=var_class,                      &
      & opt_var_ref_pos=opt_var_ref_pos, initval_s=initval,                &
      & missval_s=missval, resetval_s=resetval)
    IF (.NOT. ASSOCIATED(target_element%field%s_ptr)) &
      & CALL finish(routine, TRIM(refname)//' not created.')
    SELECT CASE(new_list_element%field%info%var_ref_pos)
    CASE(1)
      ptr => target_element%field%s_ptr(icontainer,:,:,1,1)
    CASE(2)
      ptr => target_element%field%s_ptr(:,icontainer,:,1,1)
    CASE(3)
      ptr => target_element%field%s_ptr(:,:,icontainer,1,1)
    CASE default
      CALL finish(routine, "internal error!")
    END SELECT
    new_list_element%field%s_ptr => target_element%field%s_ptr
    IF (.NOT. ASSOCIATED(new_list_element%field%s_ptr)) &
      & WRITE (0,*) 'problem with association of ptr for '//TRIM(refname)
    IF (PRESENT(lmiss)) THEN
      ptr = new_list_element%field%info%missval%sval
    ELSE
      ptr = 0.0_sp
    END IF
  END SUBROUTINE add_var_list_reference_s2d

  SUBROUTINE add_var_list_reference_i2d(this_list, target_name, refname, ptr,    &
    & hgrid, vgrid, cf, grib2, ref_idx, ldims, loutput, lrestart, lrestart_cont, &
    & initval, isteptype, resetval, lmiss, missval, tlev_source, tracer_info,    &
    & info, vert_interp, hor_interp, in_group, verbose, new_element,             &
    & l_pp_scheduler_task, post_op, action_list, opt_var_ref_pos, var_class)
    TYPE(t_var_list_ptr), INTENT(INOUT) :: this_list
    CHARACTER(len=*), INTENT(IN) :: target_name, refname
    INTEGER(sp), POINTER :: ptr(:,:)
    INTEGER, INTENT(IN) :: hgrid, vgrid, ref_idx, ldims(2)
    TYPE(t_cf_var), INTENT(IN) :: cf
    TYPE(t_grib2_var), INTENT(IN) :: grib2
    LOGICAL, INTENT(IN), OPTIONAL :: loutput, lrestart, lrestart_cont, &
      & lmiss, in_group(:), verbose
    INTEGER, INTENT(IN), OPTIONAL :: isteptype, tlev_source, var_class, &
      & l_pp_scheduler_task, opt_var_ref_pos, initval, resetval, missval
    CLASS(t_tracer_meta), INTENT(IN), OPTIONAL :: tracer_info
    TYPE(t_var_metadata), POINTER, OPTIONAL :: info
    TYPE(t_vert_interp_meta),INTENT(IN), OPTIONAL :: vert_interp
    TYPE(t_hor_interp_meta), INTENT(IN), OPTIONAL :: hor_interp
    TYPE(t_list_element), POINTER, OPTIONAL :: new_element
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op
    TYPE(t_var_action), INTENT(IN), OPTIONAL :: action_list
    CHARACTER(*), PARAMETER :: routine = modname//"::add_var_list_reference_i2d"
    TYPE(t_list_element), POINTER :: target_element, new_list_element
    INTEGER :: icontainer

    CALL add_var_list_reference_util(target_element, new_list_element,     &
      & this_list, target_name, refname, hgrid, vgrid, cf, grib2, ref_idx, &
      & ldims, INT_T, icontainer, loutput=loutput, lrestart=lrestart,      &
      & lrestart_cont=lrestart_cont, isteptype=isteptype, lmiss=lmiss,     &
      & tlev_source=tlev_source, tracer_info=tracer_info, info=info,       &
      & vert_interp=vert_interp, hor_interp=hor_interp, in_group=in_group, &
      & verbose=verbose, new_element=new_element,                          &
      & l_pp_scheduler_task=l_pp_scheduler_task, post_op=post_op,          &
      & action_list=action_list, var_class=var_class,                      &
      & opt_var_ref_pos=opt_var_ref_pos, initval_i=initval,                &
      & missval_i=missval, resetval_i=resetval)
    IF (.NOT. ASSOCIATED(target_element%field%i_ptr)) &
      & CALL finish(routine, TRIM(refname)//' not created.')
    SELECT CASE(new_list_element%field%info%var_ref_pos)
    CASE(1)
      ptr => target_element%field%i_ptr(icontainer,:,:,1,1)
    CASE(2)
      ptr => target_element%field%i_ptr(:,icontainer,:,1,1)
    CASE(3)
      ptr => target_element%field%i_ptr(:,:,icontainer,1,1)
    CASE default
      CALL finish(routine, "internal error!")
    END SELECT
    new_list_element%field%i_ptr => target_element%field%i_ptr
    IF (.NOT. ASSOCIATED(new_list_element%field%i_ptr)) &
      & WRITE (0,*) 'problem with association of ptr for '//TRIM(refname)
    IF (PRESENT(lmiss)) THEN
      ptr = new_list_element%field%info%missval%ival
    ELSE
      ptr = 0
    END IF
  END SUBROUTINE add_var_list_reference_i2d

  SUBROUTINE print_var_list (this_list, lshort)
    TYPE(t_var_list_ptr),  INTENT(in) :: this_list ! list
    LOGICAL, OPTIONAL :: lshort
    TYPE(t_var_metadata), POINTER :: info
    TYPE(t_list_element), POINTER :: this_list_element
    CHARACTER(len=32) :: dimension_text, dtext,keytext
    INTEGER :: i, igrp, ivintp_type, ndims
    CHARACTER(len=4) :: localMode
    LOGICAL :: short
    LOGICAL, POINTER :: in_group(:)

    short = .FALSE.
    IF (PRESENT(lshort)) short = lshort
    CALL message('','')
    CALL message('','')
    CALL message('','Status of variable list '//TRIM(this_list%p%name)//':')
    CALL message('','')
    this_list_element => this_list%p%first_list_element
    DO WHILE (ASSOCIATED(this_list_element))
      info => this_list_element%field%info
      IF (short) THEN
        localMode = '----'
        IF (info%name /= '' .AND. .NOT. info%lcontainer) THEN
          IF (info%lrestart) localMode(1:1) = 'r'
          IF (info%lcontained) localMode(2:2) = 't'
          SELECT CASE (info%isteptype)
          CASE (1)
            localMode(3:3) = 'i'
          CASE (2)
            localMode(3:3) = 'm'
          CASE (3)
            localMode(3:3) = 'a'
          END SELECT
          SELECT CASE (info%hgrid)
          CASE (1)
            localMode(4:4) = 'c'
          CASE (2)
            localMode(4:4) = 'v'
          CASE (3)
            localMode(4:4) = 'e'
          CASE (45)
            localMode(4:4) = 'L'
          END SELECT
          WRITE(message_text, '(a4,3i4,a24,a40)') localMode, &
            & info%grib2%discipline, info%grib2%category,   &
            & info%grib2%number, TRIM(info%name), TRIM(info%cf%standard_name)
          CALL message('', message_text)
        ENDIF
      ELSE
        IF (info%name /= '' .AND. .NOT. info%lcontainer) THEN
          WRITE (message_text,'(a,a)')       &
               'Table entry name                            : ', TRIM(info%name)
          CALL message('', message_text)
          WRITE (keytext,'(i32.1)') info%key
          WRITE (message_text,'(a,a)')       &
               'Key entry                                   : ', TRIM(keytext)
          CALL message('', message_text)
          IF (ASSOCIATED(this_list_element%field%r_ptr) .OR. &
            & ASSOCIATED(this_list_element%field%s_ptr) .OR. &
            & ASSOCIATED(this_list_element%field%i_ptr) .OR. &
            & ASSOCIATED(this_list_element%field%l_ptr)) THEN
            CALL message ('','Pointer status                              : in use.')
            dimension_text = '('
            ndims = info%ndims
            DO i = 1, ndims
              WRITE(dtext,'(i0)') info%used_dimensions(i)
              IF (info%ndims == i) THEN
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
          WRITE (message_text,'(a,3i4)') &
            & 'Assigned GRIB discipline/category/parameter : ', &
            & info%grib2%discipline, info%grib2%category, info%grib2%number
          CALL message('', message_text)
          WRITE (message_text,'(a,a,a,a)')                          &
               'CF convention standard name/unit            : ',    &
               TRIM(info%cf%standard_name), '     ', TRIM(info%cf%units)
          CALL message('', message_text)
          WRITE (message_text,'(2a)') &
            & 'CF convention long name                     : ', TRIM(info%cf%long_name)
          IF (info%lcontained) THEN
            CALL message('', 'Field is in a container                     : yes.')
            WRITE (message_text,'(a,i2)')                        &
               ' Index in container                          : ', info%ncontained
            CALL message('', message_text)
          ELSE
            CALL message('', 'Field is in a container                     : no.')
            WRITE (message_text,'(a)')                           &
               ' Index in container                          : --'
            CALL message('', message_text)
          ENDIF
          WRITE (message_text,'(a,i2)')                          &
            & ' horizontal grid type used (C=1,V=2,E=3)     : ', info%hgrid
          CALL message('', message_text)
          WRITE (message_text,'(a,i2)')                          &
            & ' vertical grid type used (see cdilib.c)      : ', info%vgrid
          CALL message('', message_text)
          WRITE (message_text,'(a,i2)')                          &
            & ' type of stat. processing (I=1,AVG=2,ACC=3...: ', info%isteptype
          CALL message('', message_text)
          IF (info%lmiss) THEN
            IF (ASSOCIATED(this_list_element%field%r_ptr)) THEN
              WRITE (message_text,'(a,e20.12)')      &
                & 'Missing value                               : ', info%missval%rval
            ELSE IF (ASSOCIATED(this_list_element%field%s_ptr)) THEN
              WRITE (message_text,'(a,e20.12)')      &
                & 'Missing value                               : ', info%missval%sval
            ELSE IF (ASSOCIATED(this_list_element%field%i_ptr)) THEN
              WRITE (message_text,'(a,i8)')      &
                & 'Missing value                               : ', info%missval%ival
            ELSE IF (ASSOCIATED(this_list_element%field%l_ptr)) THEN
              WRITE (message_text,'(a,l8)')      &
                & 'Missing value                               : ', info%missval%lval
            ENDIF
            CALL message('', message_text)
          ELSE
            CALL message('', 'Missing values                              : off.')
          ENDIF
          CALL message('', 'Added to restart                            : ' // &
            & MERGE("yes.", " no.", info%lrestart))
          CALL message('', 'Tracer field                                : ' // &
            & MERGE("yes.", " no.", this_list_element%field%info_dyn%tracer%lis_tracer))
          IF (this_list_element%field%info_dyn%tracer%lis_tracer) THEN
            CALL message('', 'Child-to-parent feedback                  : ' // &
              & MERGE("yes.", " no.", this_list_element%field%info_dyn%tracer%lfeedback))
            WRITE (message_text,'(a,3i3)') 'Horizontal transport method                 : ', &
               this_list_element%field%info_dyn%tracer%ihadv_tracer
            CALL message('', message_text)
            WRITE (message_text,'(a,3i3)') 'Vertical transport method                   : ', &
               this_list_element%field%info_dyn%tracer%ivadv_tracer
            CALL message('', message_text)
            CALL message('', 'Turbulent transport                         : ' // &
              & MERGE("yes.", " no.", this_list_element%field%info_dyn%tracer%lturb_tracer))
          ENDIF !lis_tracer
          ! print variable class/species
          WRITE (message_text,'(a,i2)')       &
            & 'Variable class/species                      : ', info%var_class
          CALL message('', message_text)
          ! print groups, to which this variable belongs:
          in_group => info%in_group
          IF (ANY(in_group(:))) THEN
            WRITE (message_text,'(a)')  'Variable group(s)                           :'
            DO igrp=1,SIZE(in_group)
              IF (in_group(igrp)) THEN
                IF (igrp == 1) THEN
                  message_text = TRIM(message_text)//" "//TRIM(var_groups_dyn%name(igrp))
                ELSE
                  message_text = TRIM(message_text)//", "//TRIM(var_groups_dyn%name(igrp))
                END IF
              ENDIF
            END DO
            CALL message('', message_text)
          END IF
          ! print horizontal and vertical interpolation method(s):
          WRITE (message_text,'(a)')  &
            &  'Horizontal interpolation                    : '//  &
            &  TRIM(STR_HINTP_TYPE(info%hor_interp%hor_intp_type))
          CALL message('', message_text)
          LOOP_VINTP_TYPES : DO ivintp_type=1,SIZE(VINTP_TYPE_LIST)
            IF (info%vert_interp%vert_intp_type(ivintp_type)) THEN
              WRITE (message_text,'(a)')  &
                &  'Vertical interpolation                      : '//  &
                &  toupper(VINTP_TYPE_LIST(ivintp_type))
              CALL message('', message_text)
            END IF
          END DO LOOP_VINTP_TYPES
          CALL message('', '')
        END IF
      END IF
      ! select next element in linked list
      this_list_element => this_list_element%next_list_element
    ENDDO
  END SUBROUTINE print_var_list

  !-----------------------------------------------------------------------------
  ! Should be overloaded to be able to search for the different information 
  ! In the proposed structure for the linked list, in the example only
  ! A character string is used so it is straight forward only one find
  FUNCTION find_list_element (this_list, name, opt_hgrid, opt_with_tl, opt_output) RESULT(element)
    TYPE(t_var_list_ptr),   INTENT(in) :: this_list
    CHARACTER(len=*),   INTENT(in) :: name
    INTEGER, OPTIONAL              :: opt_hgrid
    LOGICAL, OPTIONAL              :: opt_with_tl, opt_output
    TYPE(t_list_element), POINTER :: element
    INTEGER :: key, hgrid, time_lev
    LOGICAL :: with_tl, omit_output, with_output

    with_tl = .TRUE.
    IF (PRESENT(opt_with_tl)) with_tl = opt_with_tl
    hgrid = -1
    IF (PRESENT(opt_hgrid)) hgrid = opt_hgrid
    IF (with_tl) THEN
      key = text_hash_c(TRIM(name))
    ELSE
      key = text_hash_c(tolower(name))
    END IF
    with_output = .TRUE.
    omit_output = .NOT.PRESENT(opt_output)
    IF (.NOT.omit_output) with_output = opt_output
    time_lev = get_var_timelevel(name)
    element => this_list%p%first_list_element
    DO WHILE (ASSOCIATED(element))
      IF (-1 == hgrid .OR. hgrid == element%field%info%hgrid) THEN
        IF (MERGE(.TRUE., with_output .EQV. element%field%info%loutput, omit_output)) THEN
          IF (time_lev .EQ. get_var_timelevel(element%field%info%name)) THEN
            IF (key .EQ. MERGE(element%field%info%key, &
              &                element%field%info%key_notl, with_tl)) EXIT
          END IF
        END IF
      END IF
      element => element%next_list_element
    ENDDO
  END FUNCTION find_list_element
  
END MODULE mo_var_list
