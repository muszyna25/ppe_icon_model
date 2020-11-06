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
    &                            get_var_timelevel, get_var_name
  USE mo_tracer_metadata_types, ONLY: t_tracer_meta
  USE mo_var,              ONLY: t_var, t_var_ptr, level_type_ml
  USE mo_exception,        ONLY: message, finish, message_text
  USE mo_util_texthash,    ONLY: text_hash_c
  USE mo_util_string,      ONLY: toupper, tolower
  USE mo_impl_constants,   ONLY: STR_HINTP_TYPE, REAL_T, SINGLE_T, BOOL_T, INT_T
  USE mo_fortran_tools,    ONLY: init_contiguous_dp, init_contiguous_sp, &
    &                            init_contiguous_i4, init_contiguous_l
  USE mo_action_types,     ONLY: t_var_action

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: add_var, add_ref, find_list_element
  PUBLIC :: t_var_list_ptr

  TYPE :: t_var_list
    CHARACTER(len=256) :: filename = ''
    CHARACTER(len=128) :: name = ''
    INTEGER(i8) :: memory_used = 0_i8
    CHARACTER(len=8) :: post_suf = '', rest_suf = '', init_suf = '', &
      & model_type = 'atm' ! model type (default is 'atm' for reasons)
    INTEGER :: patch_id = -1, nvars = 0, vlevel_type = level_type_ml, &
      & output_type = -1, restart_type = -1, compression_type = -1
    LOGICAL :: loutput = .TRUE., lrestart = .FALSE., linitial = .FALSE., &
      & restart_opened = .FALSE., output_opened = .FALSE., lmiss = .FALSE., &
      & lmask_boundary = .TRUE. , first = .FALSE.
    INTEGER, ALLOCATABLE :: tl(:), hgrid(:), key(:), key_notl(:)
    LOGICAL, ALLOCATABLE :: lout(:)
    TYPE(t_var_ptr), ALLOCATABLE :: vl(:)
  END TYPE t_var_list

  TYPE :: t_var_list_ptr
    TYPE(t_var_list), POINTER :: p => NULL()
  CONTAINS
    PROCEDURE :: register => register_list_element
    PROCEDURE :: delete => delete_list
    PROCEDURE :: print => print_var_list
  END TYPE t_var_list_ptr

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

  ! remove all elements of a linked list
  SUBROUTINE delete_list(this)
    CLASS(t_var_list_ptr), INTENT(INOUT) :: this
    INTEGER :: i, n

    IF (ASSOCIATED(this%p)) THEN
      IF (ALLOCATED(this%p%vl)) THEN
        n = SIZE(this%p%vl)
        DO i = 1, n
          IF (ASSOCIATED(this%p%vl(i)%p)) THEN
            IF (this%p%vl(i)%p%info%allocated) THEN
              SELECT CASE(this%p%vl(i)%p%info%data_type)
              CASE(REAL_T)
!$ACC EXIT DATA DELETE(this%p%vl(i)%p%r_ptr) IF(this%p%vl(i)%p%info%lopenacc)
                DEALLOCATE(this%p%vl(i)%p%r_ptr)
              CASE(SINGLE_T)
!$ACC EXIT DATA DELETE(this%p%vl(i)%p%s_ptr) IF(this%p%vl(i)%p%info%lopenacc)
                DEALLOCATE(this%p%vl(i)%p%s_ptr)
              CASE(INT_T)
!$ACC EXIT DATA DELETE(this%p%vl(i)%p%i_ptr) IF(this%p%vl(i)%p%info%lopenacc)
                DEALLOCATE(this%p%vl(i)%p%i_ptr)
              CASE(BOOL_T)
!$ACC EXIT DATA DELETE(this%p%vl(i)%p%l_ptr) IF(this%p%vl(i)%p%info%lopenacc)
                DEALLOCATE(this%p%vl(i)%p%l_ptr)
              END SELECT
            END IF
            DEALLOCATE(this%p%vl(i)%p)
          END IF
        END DO
      END IF
    END IF
  END SUBROUTINE delete_list

  !-----------------------------------------------------------------------------
  ! add a list element to the linked list
  SUBROUTINE register_list_element(this, varp)
    CLASS(t_var_list_ptr), INTENT(INOUT) :: this
    TYPE(t_var), INTENT(IN), POINTER :: varp
    INTEGER :: iv, na , nv
    TYPE(t_var_ptr), ALLOCATABLE :: vtmp(:)
    CHARACTER(*), PARAMETER :: routine = modname//":register_list_element"
    INTEGER, ALLOCATABLE :: itmp1(:), itmp2(:), itmp3(:), itmp4(:)
    LOGICAL, ALLOCATABLE :: ltmp(:)

    IF (.NOT.ASSOCIATED(this%p)) CALL finish(routine, "not a valid var_list")
    IF (.NOT.ASSOCIATED(varp)) CALL finish(routine, "not a valid var")
    na = 0
    nv = this%p%nvars
    IF (nv .EQ. 0) THEN
      na = 16
    ELSE IF (SIZE(this%p%vl) .EQ. nv) THEN
      na = nv + MAX(8, nv / 8)
    END IF
    IF (na .GT. 0) THEN
      ALLOCATE(itmp1(na), itmp2(na), itmp3(na), itmp4(na), ltmp(na), vtmp(na))
      IF (nv .GT. 0) THEN
        itmp1(1:nv) = this%p%tl(1:nv)
        itmp2(1:nv) = this%p%hgrid(1:nv)
        itmp3(1:nv) = this%p%key(1:nv)
        itmp4(1:nv) = this%p%key_notl(1:nv)
        ltmp(1:nv) = this%p%lout(1:nv)
        DO iv = 1, nv
          vtmp(iv)%p => this%p%vl(iv)%p
        END DO
      END IF
      CALL MOVE_ALLOC(itmp1, this%p%tl)
      CALL MOVE_ALLOC(itmp2, this%p%hgrid)
      CALL MOVE_ALLOC(itmp3, this%p%key)
      CALL MOVE_ALLOC(itmp4, this%p%key_notl)
      CALL MOVE_ALLOC(ltmp, this%p%lout)
      CALL MOVE_ALLOC(vtmp, this%p%vl)
    END IF
    nv = nv + 1
    this%p%vl(nv)%p => varp
    this%p%tl(nv) = get_var_timelevel(varp%info%name)
    this%p%hgrid(nv) = varp%info%hgrid
    this%p%key(nv) = text_hash_c(TRIM(varp%info%name))
    this%p%key_notl(nv) = text_hash_c(tolower(get_var_name(varp%info)))
    this%p%lout(nv) = varp%info%loutput
    this%p%nvars = nv
  END SUBROUTINE register_list_element

  SUBROUTINE inherit_var_list_metadata(this, info)
    CLASS(t_var_list_ptr), INTENT(IN) :: this
    TYPE(t_var_metadata), INTENT(OUT) :: info

    info%grib2               = grib2_var(-1, -1, -1, -1, -1, -1)
    info%lrestart            = this%p%lrestart
    info%lmiss               = this%p%lmiss
    info%lmask_boundary      = this%p%lmask_boundary
    info%vert_interp         = create_vert_interp_metadata()
    info%hor_interp          = create_hor_interp_metadata()
    info%in_group(:)         = groups()
  END SUBROUTINE inherit_var_list_metadata

  !------------------------------------------------------------------------------------------------
  ! Create a list new entry
  SUBROUTINE add_var_list_element_5d(data_type, list, varname, hgrid, vgrid, &
    & cf, grib2, ldims, new_elem, loutput, lcontainer, lrestart,             &
    & lrestart_cont, isteptype, lmiss, tlev_source, info, vert_interp,       &
    & hor_interp, in_group, l_pp_scheduler_task, post_op, action_list,       &
    & tracer_info, p5_r, p5_s, p5_i, p5_l, initval_r, initval_s, initval_i,  &
    & initval_l, resetval_r, resetval_s, resetval_i, resetval_l, new_element,&
    & missval_r, missval_s, missval_i, missval_l, var_class, lopenacc)
    INTEGER, INTENT(IN) :: data_type, hgrid, vgrid, ldims(:)
    TYPE(t_var_list_ptr), INTENT(INOUT) :: list
    CHARACTER(*), INTENT(IN) :: varname
    TYPE(t_cf_var), INTENT(IN) :: cf
    TYPE(t_grib2_var), INTENT(IN) :: grib2
    TYPE(t_var), POINTER, INTENT(OUT) :: new_elem
    LOGICAL, INTENT(IN), OPTIONAL :: loutput, lcontainer, lrestart, lrestart_cont, &
      & lmiss, in_group(:), initval_l, resetval_l, missval_l, lopenacc
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
    TYPE(t_var), POINTER, INTENT(OUT), OPTIONAL :: new_element
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
    IF (ASSOCIATED(find_list_element(list, varname))) &
      & CALL finish(routine, "duplicate entry ("//TRIM(varname)//") in var_list ("//TRIM(list%p%name)//")")
    is_restart_var = list%p%lrestart
    IF (PRESENT(lrestart)) THEN
      is_restart_var = lrestart
      IF (.NOT.list%p%lrestart .AND. lrestart) &
        & CALL finish(routine, 'for list '//TRIM(list%p%name)//' restarting not enabled, '// &
                           & 'but restart of '//TRIM(varname)//' requested.')
    ENDIF
    IF (is_restart_var .AND. (.NOT. ANY(data_type == (/REAL_T, SINGLE_T, INT_T/)))) &
      & CALL finish(routine, 'unsupported data_type for "'//TRIM(varname)//'": '// &
        & 'data_type of restart variables must be floating-point or integer type.')
    ALLOCATE(new_elem)
    CALL inherit_var_list_metadata(list, new_elem%info)
    ! init local fields
    missval = new_elem%info%missval
    initval = new_elem%info%initval
    resetval= new_elem%info%resetval
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
    CALL set_var_metadata(new_elem%info, ldims, name=varname,       &
      & hgrid=hgrid, vgrid=vgrid, cf=cf, grib2=grib2, loutput=loutput,       &
      & lcontainer=lcontainer, lrestart=lrestart, missval=missval,           &
      & lrestart_cont=lrestart_cont, initval=initval, isteptype=isteptype,   &
      & resetval=resetval, tlev_source=tlev_source, vert_interp=vert_interp, &
      & hor_interp=hor_interp, l_pp_scheduler_task=l_pp_scheduler_task,      &
      & post_op=post_op, action_list=action_list, var_class=var_class,       &
      & data_type=data_type, lopenacc=lopenacc, lmiss=lmiss, in_group=in_group)
    ! set dynamic metadata, i.e. polymorphic tracer metadata
    CALL set_var_metadata_dyn (new_elem%info_dyn, tracer_info=tracer_info)
    new_elem%info%ndims = ndims
    new_elem%info%used_dimensions(1:ndims) = ldims(1:ndims)
    new_elem%info%dom => list%p%patch_id
    IF(PRESENT(info)) info => new_elem%info
    NULLIFY(new_elem%r_ptr, new_elem%s_ptr, new_elem%i_ptr, new_elem%l_ptr)
    d(1:ndims)    = new_elem%info%used_dimensions(1:ndims)
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
        new_elem%r_ptr => p5_r
      ELSE
        new_elem%var_base_size = 8
        ALLOCATE(new_elem%r_ptr(d(1), d(2), d(3), d(4), d(5)), STAT=istat)
        IF (istat /= 0) CALL finish(routine, 'allocation of array '//TRIM(varname)//' failed')
        !$ACC ENTER DATA CREATE(new_elem%r_ptr) IF(new_elem%info%lopenacc)
      END IF
      !ICON_OMP PARALLEL
      CALL init_contiguous_dp(new_elem%r_ptr, PRODUCT(d(1:5)), ivals%rval)
      !ICON_OMP END PARALLEL
      !$ACC UPDATE DEVICE(new_elem%r_ptr) IF(new_elem%info%lopenacc)
    CASE(SINGLE_T)
      IF (referenced) THEN
        new_elem%s_ptr => p5_s
      ELSE
        new_elem%var_base_size = 4
        ALLOCATE(new_elem%s_ptr(d(1), d(2), d(3), d(4), d(5)), STAT=istat)
        IF (istat /= 0) CALL finish(routine, 'allocation of array '//TRIM(varname)//' failed')
        !$ACC ENTER DATA CREATE(new_elem%s_ptr) IF(new_elem%info%lopenacc)
      END IF
      !ICON_OMP PARALLEL
      CALL init_contiguous_sp(new_elem%s_ptr, PRODUCT(d(1:5)), ivals%sval)
      !ICON_OMP END PARALLEL
      !$ACC UPDATE DEVICE(new_elem%s_ptr) IF(new_elem%info%lopenacc)
    CASE(INT_T)
      IF (referenced) THEN
        new_elem%i_ptr => p5_i
      ELSE
        new_elem%var_base_size = 4
        ALLOCATE(new_elem%i_ptr(d(1), d(2), d(3), d(4), d(5)), STAT=istat)
        IF (istat /= 0) CALL finish(routine, 'allocation of arrayb'//TRIM(varname)//' failed')
        !$ACC ENTER DATA CREATE(new_elem%i_ptr) IF(new_elem%info%lopenacc)
      END IF
      !ICON_OMP PARALLEL
      CALL init_contiguous_i4(new_elem%i_ptr, PRODUCT(d(1:5)), ivals%ival)
      !ICON_OMP END PARALLEL
      !$ACC UPDATE DEVICE(new_elem%i_ptr) IF(new_elem%info%lopenacc)
    CASE(BOOL_T)
      IF (referenced) THEN
        new_elem%l_ptr => p5_l
      ELSE
        new_elem%var_base_size = 4
        ALLOCATE(new_elem%l_ptr(d(1), d(2), d(3), d(4), d(5)), STAT=istat)
        IF (istat /= 0) CALL finish(routine, 'allocation of array '//TRIM(varname)//' failed')
        !$ACC ENTER DATA CREATE(new_elem%l_ptr) IF(new_elem%info%lopenacc)
      END IF
      !ICON_OMP PARALLEL
      CALL init_contiguous_l(new_elem%l_ptr, PRODUCT(d(1:5)), ivals%lval)
      !ICON_OMP END PARALLEL
      !$ACC UPDATE DEVICE(new_elem%l_ptr) IF(new_elem%info%lopenacc)
    END SELECT
    CALL register_list_element(list, new_elem)
    IF (.NOT.referenced) list%p%memory_used = list%p%memory_used + &
      & INT(new_elem%var_base_size, i8) * INT(PRODUCT(d(1:5)),i8)
    new_elem%info%allocated = .TRUE.
    IF (PRESENT(new_element)) new_element => new_elem
  END SUBROUTINE add_var_list_element_5d

  SUBROUTINE add_var_list_element_r4d(this_list, varname, ptr, hgrid, vgrid, &
    & cf, grib2, ldims, loutput, lcontainer, lrestart, lrestart_cont,     &
    & initval, isteptype, resetval, lmiss, missval, tlev_source, info,    &
    & p5, vert_interp, hor_interp, in_group, new_element,        &
    & l_pp_scheduler_task, post_op, action_list, var_class, lopenacc)
    TYPE(t_var_list_ptr), INTENT(INOUT) :: this_list
    CHARACTER(*), INTENT(IN) :: varname
    REAL(dp), POINTER, INTENT(OUT) :: ptr(:,:,:,:)
    INTEGER, INTENT(IN) :: hgrid, vgrid, ldims(4)
    TYPE(t_cf_var), INTENT(IN) :: cf
    TYPE(t_grib2_var), INTENT(IN) :: grib2
    LOGICAL, INTENT(IN), OPTIONAL :: loutput, lcontainer, lrestart, &
      & lrestart_cont, lmiss, in_group(:), lopenacc
    INTEGER, INTENT(IN), OPTIONAL :: isteptype, l_pp_scheduler_task, &
      & tlev_source, var_class
    REAL(dp), INTENT(IN), OPTIONAL :: initval, resetval, missval
    TYPE(t_var_metadata), POINTER, OPTIONAL :: info
    REAL(dp), CONTIGUOUS_TARGET, OPTIONAL :: p5(:,:,:,:,:)
    TYPE(t_vert_interp_meta), INTENT(IN), OPTIONAL :: vert_interp
    TYPE(t_hor_interp_meta), INTENT(IN), OPTIONAL :: hor_interp
    TYPE(t_var), POINTER, OPTIONAL :: new_element
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op
    TYPE(t_var_action), INTENT(IN), OPTIONAL :: action_list
    TYPE(t_var), POINTER :: element

    CALL add_var_list_element_5d(REAL_T, this_list, varname, hgrid, vgrid, &
      & cf, grib2, ldims, element, loutput, lcontainer, lrestart,          &
      & lrestart_cont, isteptype, lmiss, tlev_source, info, vert_interp,   &
      & hor_interp, in_group, l_pp_scheduler_task, post_op, action_list,   &
      & p5_r=p5, initval_r=initval, resetval_r=resetval, missval_r=missval,&
      & var_class=var_class, lopenacc=lopenacc, new_element=new_element)
    ptr => element%r_ptr(:,:,:,:,1)
  END SUBROUTINE add_var_list_element_r4d

  SUBROUTINE add_var_list_element_r3d(this_list, varname, ptr, hgrid, vgrid, &
    & cf, grib2, ldims, loutput, lcontainer, lrestart, lrestart_cont,     &
    & initval, isteptype, resetval, lmiss, missval, tlev_source, info,    &
    & p5, vert_interp, hor_interp, in_group, new_element, tracer_info, &
    & l_pp_scheduler_task, post_op, action_list, var_class, lopenacc)
    TYPE(t_var_list_ptr), INTENT(INOUT) :: this_list
    CHARACTER(*), INTENT(IN) :: varname
    REAL(dp), POINTER, INTENT(OUT) :: ptr(:,:,:)
    INTEGER, INTENT(IN) :: hgrid, vgrid, ldims(3)
    TYPE(t_cf_var), INTENT(IN) :: cf
    TYPE(t_grib2_var), INTENT(IN) :: grib2
    LOGICAL, INTENT(IN), OPTIONAL :: loutput, lcontainer, lrestart, &
      & lrestart_cont, lmiss, in_group(:), lopenacc
    INTEGER, INTENT(IN), OPTIONAL :: isteptype, l_pp_scheduler_task, &
      & tlev_source, var_class
    REAL(dp), INTENT(IN), OPTIONAL :: initval, resetval, missval
    TYPE(t_var_metadata), POINTER, OPTIONAL :: info
    CLASS(t_tracer_meta), INTENT(in), OPTIONAL :: tracer_info
    REAL(dp), CONTIGUOUS_TARGET, OPTIONAL :: p5(:,:,:,:,:)
    TYPE(t_vert_interp_meta), INTENT(IN), OPTIONAL :: vert_interp
    TYPE(t_hor_interp_meta), INTENT(IN), OPTIONAL :: hor_interp
    TYPE(t_var), POINTER, OPTIONAL :: new_element
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op
    TYPE(t_var_action), INTENT(IN), OPTIONAL :: action_list
    TYPE(t_var), POINTER :: element

    CALL add_var_list_element_5d(REAL_T, this_list, varname, hgrid, vgrid, &
      & cf, grib2, ldims, element, loutput, lcontainer, lrestart,          &
      & lrestart_cont, isteptype, lmiss, tlev_source, info, vert_interp,   &
      & hor_interp, in_group, l_pp_scheduler_task, post_op, action_list,   &
      & p5_r=p5, initval_r=initval, resetval_r=resetval, missval_r=missval,&
      & var_class=var_class, lopenacc=lopenacc, new_element=new_element,   &
      & tracer_info=tracer_info)
    ptr => element%r_ptr(:,:,:,1,1)
  END SUBROUTINE add_var_list_element_r3d

  SUBROUTINE add_var_list_element_r2d(this_list, varname, ptr, hgrid, vgrid, &
    & cf, grib2, ldims, loutput, lcontainer, lrestart, lrestart_cont,     &
    & initval, isteptype, resetval, lmiss, missval, tlev_source, info,    &
    & p5, vert_interp, hor_interp, in_group, new_element, tracer_info, &
    & l_pp_scheduler_task, post_op, action_list, var_class, lopenacc)
    TYPE(t_var_list_ptr), INTENT(INOUT) :: this_list
    CHARACTER(*), INTENT(IN) :: varname
    REAL(dp), POINTER, INTENT(OUT) :: ptr(:,:)
    INTEGER, INTENT(IN) :: hgrid, vgrid, ldims(2)
    TYPE(t_cf_var), INTENT(IN) :: cf
    TYPE(t_grib2_var), INTENT(IN) :: grib2
    LOGICAL, INTENT(IN), OPTIONAL :: loutput, lcontainer, lrestart, &
      & lrestart_cont, lmiss, in_group(:), lopenacc
    INTEGER, INTENT(IN), OPTIONAL :: isteptype, l_pp_scheduler_task, &
      & tlev_source, var_class
    REAL(dp), INTENT(IN), OPTIONAL :: initval, resetval, missval
    TYPE(t_var_metadata), POINTER, OPTIONAL :: info
    CLASS(t_tracer_meta), INTENT(IN), OPTIONAL :: tracer_info
    REAL(dp), CONTIGUOUS_TARGET, OPTIONAL :: p5(:,:,:,:,:)
    TYPE(t_vert_interp_meta), INTENT(IN), OPTIONAL :: vert_interp
    TYPE(t_hor_interp_meta), INTENT(IN), OPTIONAL :: hor_interp
    TYPE(t_var), POINTER, OPTIONAL :: new_element
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op
    TYPE(t_var_action), INTENT(IN), OPTIONAL :: action_list
    TYPE(t_var), POINTER :: element

    CALL add_var_list_element_5d(REAL_T, this_list, varname, hgrid, vgrid, &
      & cf, grib2, ldims, element, loutput, lcontainer, lrestart,          &
      & lrestart_cont, isteptype, lmiss, tlev_source, info, vert_interp,   &
      & hor_interp, in_group, l_pp_scheduler_task, post_op, action_list,   &
      & p5_r=p5, initval_r=initval, resetval_r=resetval, missval_r=missval,&
      & var_class=var_class, lopenacc=lopenacc, new_element=new_element,   &
      & tracer_info=tracer_info)
    ptr => element%r_ptr(:,:,1,1,1)
  END SUBROUTINE add_var_list_element_r2d

  SUBROUTINE add_var_list_element_r1d(this_list, varname, ptr, hgrid, vgrid, &
    & cf, grib2, ldims, loutput, lcontainer, lrestart, lrestart_cont,     &
    & initval, isteptype, resetval, lmiss, missval, tlev_source, info,    &
    & p5, vert_interp, hor_interp, in_group, new_element,        &
    & l_pp_scheduler_task, post_op, action_list, var_class, lopenacc)
    TYPE(t_var_list_ptr), INTENT(INOUT) :: this_list
    CHARACTER(*), INTENT(IN) :: varname
    REAL(dp), POINTER, INTENT(OUT) :: ptr(:)
    INTEGER, INTENT(IN) :: hgrid, vgrid, ldims(1)
    TYPE(t_cf_var), INTENT(IN) :: cf
    TYPE(t_grib2_var), INTENT(IN) :: grib2
    LOGICAL, INTENT(IN), OPTIONAL :: loutput, lcontainer, lrestart, &
      & lrestart_cont, lmiss, in_group(:), lopenacc
    INTEGER, INTENT(IN), OPTIONAL :: isteptype, l_pp_scheduler_task, &
      & tlev_source, var_class
    REAL(dp), INTENT(IN), OPTIONAL :: initval, resetval, missval
    TYPE(t_var_metadata), POINTER, OPTIONAL :: info
    REAL(dp), CONTIGUOUS_TARGET, OPTIONAL :: p5(:,:,:,:,:)
    TYPE(t_vert_interp_meta), INTENT(IN), OPTIONAL :: vert_interp
    TYPE(t_hor_interp_meta), INTENT(IN), OPTIONAL :: hor_interp
    TYPE(t_var), POINTER, OPTIONAL :: new_element
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op
    TYPE(t_var_action), INTENT(IN), OPTIONAL :: action_list
    TYPE(t_var), POINTER :: element

    CALL add_var_list_element_5d(REAL_T, this_list, varname, hgrid, vgrid, &
      & cf, grib2, ldims, element, loutput, lcontainer, lrestart,          &
      & lrestart_cont, isteptype, lmiss, tlev_source, info, vert_interp,   &
      & hor_interp, in_group, l_pp_scheduler_task, post_op, action_list,   &
      & p5_r=p5, initval_r=initval, resetval_r=resetval, missval_r=missval,&
      & var_class=var_class, lopenacc=lopenacc, new_element=new_element)
    ptr => element%r_ptr(:,1,1,1,1)
  END SUBROUTINE add_var_list_element_r1d

  SUBROUTINE add_var_list_element_s4d(this_list, varname, ptr, hgrid, vgrid, &
    & cf, grib2, ldims, loutput, lcontainer, lrestart, lrestart_cont,     &
    & initval, isteptype, resetval, lmiss, missval, tlev_source, info,    &
    & p5, vert_interp, hor_interp, in_group, new_element,        &
    & l_pp_scheduler_task, post_op, action_list, var_class, lopenacc)
    TYPE(t_var_list_ptr), INTENT(INOUT) :: this_list
    CHARACTER(*), INTENT(IN) :: varname
    REAL(sp), POINTER, INTENT(OUT) :: ptr(:,:,:,:)
    INTEGER, INTENT(IN) :: hgrid, vgrid, ldims(4)
    TYPE(t_cf_var), INTENT(IN) :: cf
    TYPE(t_grib2_var), INTENT(IN) :: grib2
    LOGICAL, INTENT(IN), OPTIONAL :: loutput, lcontainer, lrestart, &
      & lrestart_cont, lmiss, in_group(:), lopenacc
    INTEGER, INTENT(IN), OPTIONAL :: isteptype, l_pp_scheduler_task, &
      & tlev_source, var_class
    REAL(sp), INTENT(IN), OPTIONAL :: initval, resetval, missval
    TYPE(t_var_metadata), POINTER, OPTIONAL :: info
    REAL(sp), CONTIGUOUS_TARGET, OPTIONAL :: p5(:,:,:,:,:)
    TYPE(t_vert_interp_meta), INTENT(IN), OPTIONAL :: vert_interp
    TYPE(t_hor_interp_meta), INTENT(IN), OPTIONAL :: hor_interp
    TYPE(t_var), POINTER, OPTIONAL :: new_element
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op
    TYPE(t_var_action), INTENT(IN), OPTIONAL :: action_list
    TYPE(t_var), POINTER :: element

    CALL add_var_list_element_5d(SINGLE_T, this_list, varname, hgrid, vgrid,&
      & cf, grib2, ldims, element, loutput, lcontainer, lrestart,          &
      & lrestart_cont, isteptype, lmiss, tlev_source, info, vert_interp,   &
      & hor_interp, in_group, l_pp_scheduler_task, post_op, action_list,   &
      & p5_s=p5, initval_s=initval, resetval_s=resetval, missval_s=missval,&
      & var_class=var_class, lopenacc=lopenacc, new_element=new_element)
    ptr => element%s_ptr(:,:,:,:,1)
  END SUBROUTINE add_var_list_element_s4d

  SUBROUTINE add_var_list_element_s3d(this_list, varname, ptr, hgrid, vgrid, &
    & cf, grib2, ldims, loutput, lcontainer, lrestart, lrestart_cont,     &
    & initval, isteptype, resetval, lmiss, missval, tlev_source, info,    &
    & p5, vert_interp, hor_interp, in_group, new_element, tracer_info, &
    & l_pp_scheduler_task, post_op, action_list, var_class, lopenacc)
    TYPE(t_var_list_ptr), INTENT(INOUT) :: this_list
    CHARACTER(*), INTENT(IN) :: varname
    REAL(sp), POINTER, INTENT(OUT) :: ptr(:,:,:)
    INTEGER, INTENT(IN) :: hgrid, vgrid, ldims(:)
    TYPE(t_cf_var), INTENT(IN) :: cf
    TYPE(t_grib2_var), INTENT(IN) :: grib2
    LOGICAL, INTENT(IN), OPTIONAL :: loutput, lcontainer, lrestart, &
      & lrestart_cont, lmiss, in_group(:), lopenacc
    INTEGER, INTENT(IN), OPTIONAL :: isteptype, l_pp_scheduler_task, &
      & tlev_source, var_class
    REAL(sp), INTENT(IN), OPTIONAL :: initval, resetval, missval
    TYPE(t_var_metadata), POINTER, OPTIONAL :: info
    CLASS(t_tracer_meta),    INTENT(in), OPTIONAL :: tracer_info
    REAL(sp), CONTIGUOUS_TARGET, OPTIONAL :: p5(:,:,:,:,:)
    TYPE(t_vert_interp_meta), INTENT(IN), OPTIONAL :: vert_interp
    TYPE(t_hor_interp_meta), INTENT(IN), OPTIONAL :: hor_interp
    TYPE(t_var), POINTER, OPTIONAL :: new_element
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op
    TYPE(t_var_action), INTENT(IN), OPTIONAL :: action_list
    TYPE(t_var), POINTER :: element

    CALL add_var_list_element_5d(SINGLE_T, this_list, varname, hgrid, vgrid,&
      & cf, grib2, ldims, element, loutput, lcontainer, lrestart,          &
      & lrestart_cont, isteptype, lmiss, tlev_source, info, vert_interp,   &
      & hor_interp, in_group, l_pp_scheduler_task, post_op, action_list,   &
      & p5_s=p5, initval_s=initval, resetval_s=resetval, missval_s=missval,&
      & var_class=var_class, lopenacc=lopenacc, new_element=new_element,   &
      & tracer_info=tracer_info)
    ptr => element%s_ptr(:,:,:,1,1)
  END SUBROUTINE add_var_list_element_s3d

  SUBROUTINE add_var_list_element_s2d(this_list, varname, ptr, hgrid, vgrid, &
    & cf, grib2, ldims, loutput, lcontainer, lrestart, lrestart_cont,     &
    & initval, isteptype, resetval, lmiss, missval, tlev_source, info,    &
    & p5, vert_interp, hor_interp, in_group, new_element, tracer_info, &
    & l_pp_scheduler_task, post_op, action_list, var_class, lopenacc)
    TYPE(t_var_list_ptr), INTENT(INOUT) :: this_list
    CHARACTER(*), INTENT(IN) :: varname
    REAL(sp), POINTER, INTENT(OUT) :: ptr(:,:)
    INTEGER, INTENT(IN) :: hgrid, vgrid, ldims(2)
    TYPE(t_cf_var), INTENT(IN) :: cf
    TYPE(t_grib2_var), INTENT(IN) :: grib2
    LOGICAL, INTENT(IN), OPTIONAL :: loutput, lcontainer, lrestart, &
      & lrestart_cont, lmiss, in_group(:), lopenacc
    INTEGER, INTENT(IN), OPTIONAL :: isteptype, l_pp_scheduler_task, &
      & tlev_source, var_class
    REAL(sp), INTENT(IN), OPTIONAL :: initval, resetval, missval
    TYPE(t_var_metadata), POINTER, OPTIONAL :: info
    CLASS(t_tracer_meta),    INTENT(in), OPTIONAL :: tracer_info
    REAL(sp), CONTIGUOUS_TARGET, OPTIONAL :: p5(:,:,:,:,:)
    TYPE(t_vert_interp_meta), INTENT(IN), OPTIONAL :: vert_interp
    TYPE(t_hor_interp_meta), INTENT(IN), OPTIONAL :: hor_interp
    TYPE(t_var), POINTER, OPTIONAL :: new_element
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op
    TYPE(t_var_action), INTENT(IN), OPTIONAL :: action_list
    TYPE(t_var), POINTER :: element

    CALL add_var_list_element_5d(SINGLE_T, this_list, varname, hgrid, vgrid,&
      & cf, grib2, ldims, element, loutput, lcontainer, lrestart,          &
      & lrestart_cont, isteptype, lmiss, tlev_source, info, vert_interp,   &
      & hor_interp, in_group, l_pp_scheduler_task, post_op, action_list,   &
      & p5_s=p5, initval_s=initval, resetval_s=resetval, missval_s=missval,&
      & var_class=var_class, lopenacc=lopenacc, new_element=new_element,   &
      & tracer_info=tracer_info)
    ptr => element%s_ptr(:,:,1,1,1)
  END SUBROUTINE add_var_list_element_s2d

  SUBROUTINE add_var_list_element_s1d(this_list, varname, ptr, hgrid, vgrid, &
    & cf, grib2, ldims, loutput, lcontainer, lrestart, lrestart_cont,     &
    & initval, isteptype, resetval, lmiss, missval, tlev_source, info,    &
    & p5, vert_interp, hor_interp, in_group, new_element,        &
    & l_pp_scheduler_task, post_op, action_list, var_class, lopenacc)
    TYPE(t_var_list_ptr), INTENT(INOUT) :: this_list
    CHARACTER(*), INTENT(IN) :: varname
    REAL(sp), POINTER, INTENT(OUT) :: ptr(:)
    INTEGER, INTENT(IN) :: hgrid, vgrid, ldims(1)
    TYPE(t_cf_var), INTENT(IN) :: cf
    TYPE(t_grib2_var), INTENT(IN) :: grib2
    LOGICAL, INTENT(IN), OPTIONAL :: loutput, lcontainer, lrestart, &
      & lrestart_cont, lmiss, in_group(:), lopenacc
    INTEGER, INTENT(IN), OPTIONAL :: isteptype, l_pp_scheduler_task, &
      & tlev_source, var_class
    REAL(sp), INTENT(IN), OPTIONAL :: initval, resetval, missval
    TYPE(t_var_metadata), POINTER, OPTIONAL :: info
    REAL(sp), CONTIGUOUS_TARGET, OPTIONAL :: p5(:,:,:,:,:)
    TYPE(t_vert_interp_meta), INTENT(IN), OPTIONAL :: vert_interp
    TYPE(t_hor_interp_meta), INTENT(IN), OPTIONAL :: hor_interp
    TYPE(t_var), POINTER, OPTIONAL :: new_element
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op
    TYPE(t_var_action), INTENT(IN), OPTIONAL :: action_list
    TYPE(t_var), POINTER :: element

    CALL add_var_list_element_5d(SINGLE_T, this_list, varname, hgrid, vgrid,& 
      & cf, grib2, ldims, element, loutput, lcontainer, lrestart,          &
      & lrestart_cont, isteptype, lmiss, tlev_source, info, vert_interp,   &
      & hor_interp, in_group, l_pp_scheduler_task, post_op, action_list,   &
      & p5_s=p5, initval_s=initval, resetval_s=resetval, missval_s=missval,&
      & var_class=var_class, lopenacc=lopenacc, new_element=new_element)
    ptr => element%s_ptr(:,1,1,1,1)
  END SUBROUTINE add_var_list_element_s1d

  SUBROUTINE add_var_list_element_i4d(this_list, varname, ptr, hgrid, vgrid, &
    & cf, grib2, ldims, loutput, lcontainer, lrestart, lrestart_cont,     &
    & initval, isteptype, resetval, lmiss, missval, tlev_source, info,    &
    & p5, vert_interp, hor_interp, in_group, new_element,        &
    & l_pp_scheduler_task, post_op, action_list, var_class, lopenacc)
    TYPE(t_var_list_ptr), INTENT(INOUT) :: this_list
    CHARACTER(*), INTENT(IN) :: varname
    INTEGER, POINTER, INTENT(OUT) :: ptr(:,:,:,:)
    INTEGER, INTENT(IN) :: hgrid, vgrid, ldims(4)
    TYPE(t_cf_var), INTENT(IN) :: cf
    TYPE(t_grib2_var), INTENT(IN) :: grib2
    LOGICAL, INTENT(IN), OPTIONAL :: loutput, lcontainer, lrestart, &
      & lrestart_cont, lmiss, in_group(:), lopenacc
    INTEGER, INTENT(IN), OPTIONAL :: isteptype, l_pp_scheduler_task, &
      & tlev_source, var_class, initval, resetval, missval
    TYPE(t_var_metadata), POINTER, OPTIONAL :: info
    INTEGER, CONTIGUOUS_TARGET, OPTIONAL :: p5(:,:,:,:,:)
    TYPE(t_vert_interp_meta), INTENT(IN), OPTIONAL :: vert_interp
    TYPE(t_hor_interp_meta), INTENT(IN), OPTIONAL :: hor_interp
    TYPE(t_var), POINTER, OPTIONAL :: new_element
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op
    TYPE(t_var_action), INTENT(IN), OPTIONAL :: action_list
    TYPE(t_var), POINTER :: element

    CALL add_var_list_element_5d(INT_T, this_list, varname, hgrid, vgrid,  & 
      & cf, grib2, ldims, element, loutput, lcontainer, lrestart,          &
      & lrestart_cont, isteptype, lmiss, tlev_source, info, vert_interp,   &
      & hor_interp, in_group, l_pp_scheduler_task, post_op, action_list,   &
      & p5_i=p5, initval_i=initval, resetval_i=resetval, missval_i=missval,&
      & var_class=var_class, lopenacc=lopenacc, new_element=new_element)
    ptr => element%i_ptr(:,:,:,:,1)
  END SUBROUTINE add_var_list_element_i4d

  SUBROUTINE add_var_list_element_i3d(this_list, varname, ptr, hgrid, vgrid, &
    & cf, grib2, ldims, loutput, lcontainer, lrestart, lrestart_cont,     &
    & initval, isteptype, resetval, lmiss, missval, tlev_source, info,    &
    & p5, vert_interp, hor_interp, in_group, new_element,        &
    & l_pp_scheduler_task, post_op, action_list, var_class, lopenacc)
    TYPE(t_var_list_ptr), INTENT(INOUT) :: this_list
    CHARACTER(*), INTENT(IN) :: varname
    INTEGER, POINTER, INTENT(OUT) :: ptr(:,:,:)
    INTEGER, INTENT(IN) :: hgrid, vgrid, ldims(3)
    TYPE(t_cf_var), INTENT(IN) :: cf
    TYPE(t_grib2_var), INTENT(IN) :: grib2
    LOGICAL, INTENT(IN), OPTIONAL :: loutput, lcontainer, lrestart, &
      & lrestart_cont, lmiss, in_group(:), lopenacc
    INTEGER, INTENT(IN), OPTIONAL :: isteptype, l_pp_scheduler_task, &
      & tlev_source, var_class, initval, resetval, missval
    TYPE(t_var_metadata), POINTER, OPTIONAL :: info
    INTEGER, CONTIGUOUS_TARGET, OPTIONAL :: p5(:,:,:,:,:)
    TYPE(t_vert_interp_meta), INTENT(IN), OPTIONAL :: vert_interp
    TYPE(t_hor_interp_meta), INTENT(IN), OPTIONAL :: hor_interp
    TYPE(t_var), POINTER, OPTIONAL :: new_element
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op
    TYPE(t_var_action), INTENT(IN), OPTIONAL :: action_list
    TYPE(t_var), POINTER :: element

    CALL add_var_list_element_5d(INT_T, this_list, varname, hgrid, vgrid,  & 
      & cf, grib2, ldims, element, loutput, lcontainer, lrestart,          &
      & lrestart_cont, isteptype, lmiss, tlev_source, info, vert_interp,   &
      & hor_interp, in_group, l_pp_scheduler_task, post_op, action_list,   &
      & p5_i=p5, initval_i=initval, resetval_i=resetval, missval_i=missval,&
      & var_class=var_class, lopenacc=lopenacc, new_element=new_element)
    ptr => element%i_ptr(:,:,:,1,1)
  END SUBROUTINE add_var_list_element_i3d

  SUBROUTINE add_var_list_element_i2d(this_list, varname, ptr, hgrid, vgrid, &
    & cf, grib2, ldims, loutput, lcontainer, lrestart, lrestart_cont,     &
    & initval, isteptype, resetval, lmiss, missval, tlev_source, info,    &
    & p5, vert_interp, hor_interp, in_group, new_element,        &
    & l_pp_scheduler_task, post_op, action_list, var_class, lopenacc)
    TYPE(t_var_list_ptr), INTENT(INOUT) :: this_list
    CHARACTER(*), INTENT(IN) :: varname
    INTEGER, POINTER, INTENT(OUT) :: ptr(:,:)
    INTEGER, INTENT(IN) :: hgrid, vgrid, ldims(2)
    TYPE(t_cf_var), INTENT(IN) :: cf
    TYPE(t_grib2_var), INTENT(IN) :: grib2
    LOGICAL, INTENT(IN), OPTIONAL :: loutput, lcontainer, lrestart, &
      & lrestart_cont, lmiss, in_group(:), lopenacc
    INTEGER, INTENT(IN), OPTIONAL :: isteptype, l_pp_scheduler_task, &
      & tlev_source, var_class, initval, resetval, missval
    TYPE(t_var_metadata), POINTER, OPTIONAL :: info
    INTEGER, CONTIGUOUS_TARGET, OPTIONAL :: p5(:,:,:,:,:)
    TYPE(t_vert_interp_meta), INTENT(IN), OPTIONAL :: vert_interp
    TYPE(t_hor_interp_meta), INTENT(IN), OPTIONAL :: hor_interp
    TYPE(t_var), POINTER, OPTIONAL :: new_element
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op
    TYPE(t_var_action), INTENT(IN), OPTIONAL :: action_list
    TYPE(t_var), POINTER :: element

    CALL add_var_list_element_5d(INT_T, this_list, varname, hgrid, vgrid,  & 
      & cf, grib2, ldims, element, loutput, lcontainer, lrestart,          &
      & lrestart_cont, isteptype, lmiss, tlev_source, info, vert_interp,   &
      & hor_interp, in_group, l_pp_scheduler_task, post_op, action_list,   &
      & p5_i=p5, initval_i=initval, resetval_i=resetval, missval_i=missval,&
      & var_class=var_class, lopenacc=lopenacc, new_element=new_element)
    ptr => element%i_ptr(:,:,1,1,1)
  END SUBROUTINE add_var_list_element_i2d

  SUBROUTINE add_var_list_element_i1d(this_list, varname, ptr, hgrid, vgrid, &
    & cf, grib2, ldims, loutput, lcontainer, lrestart, lrestart_cont,     &
    & initval, isteptype, resetval, lmiss, missval, tlev_source, info,    &
    & p5, vert_interp, hor_interp, in_group, new_element,        &
    & l_pp_scheduler_task, post_op, action_list, var_class, lopenacc)
    TYPE(t_var_list_ptr), INTENT(INOUT) :: this_list
    CHARACTER(*), INTENT(IN) :: varname
    INTEGER, POINTER, INTENT(OUT) :: ptr(:)
    INTEGER, INTENT(IN) :: hgrid, vgrid, ldims(1)
    TYPE(t_cf_var), INTENT(IN) :: cf
    TYPE(t_grib2_var), INTENT(IN) :: grib2
    LOGICAL, INTENT(IN), OPTIONAL :: loutput, lcontainer, lrestart, &
      & lrestart_cont, lmiss, in_group(:), lopenacc
    INTEGER, INTENT(IN), OPTIONAL :: isteptype, l_pp_scheduler_task, &
      & tlev_source, var_class, initval, resetval, missval
    TYPE(t_var_metadata), POINTER, OPTIONAL :: info
    INTEGER, CONTIGUOUS_TARGET, OPTIONAL :: p5(:,:,:,:,:)
    TYPE(t_vert_interp_meta), INTENT(IN), OPTIONAL :: vert_interp
    TYPE(t_hor_interp_meta), INTENT(IN), OPTIONAL :: hor_interp
    TYPE(t_var), POINTER, OPTIONAL :: new_element
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op
    TYPE(t_var_action), INTENT(IN), OPTIONAL :: action_list
    TYPE(t_var), POINTER :: element

    CALL add_var_list_element_5d(INT_T, this_list, varname, hgrid, vgrid,  &
      & cf, grib2, ldims, element, loutput, lcontainer, lrestart,          &
      & lrestart_cont, isteptype, lmiss, tlev_source, info, vert_interp,   &
      & hor_interp, in_group, l_pp_scheduler_task, post_op, action_list,   &
      & p5_i=p5, initval_i=initval, resetval_i=resetval, missval_i=missval,&
      & var_class=var_class, lopenacc=lopenacc, new_element=new_element)
    ptr => element%i_ptr(:,1,1,1,1)
  END SUBROUTINE add_var_list_element_i1d

  SUBROUTINE add_var_list_element_l4d(this_list, varname, ptr, hgrid, vgrid, &
    & cf, grib2, ldims, loutput, lcontainer, lrestart, lrestart_cont,     &
    & initval, isteptype, resetval, lmiss, missval, tlev_source, info,    &
    & p5, vert_interp, hor_interp, in_group, new_element,        &
    & l_pp_scheduler_task, post_op, action_list, var_class, lopenacc)
    TYPE(t_var_list_ptr), INTENT(INOUT) :: this_list
    CHARACTER(*), INTENT(IN) :: varname
    LOGICAL, POINTER, INTENT(OUT) :: ptr(:,:,:,:)
    INTEGER, INTENT(IN) :: hgrid, vgrid, ldims(4)
    TYPE(t_cf_var), INTENT(IN) :: cf
    TYPE(t_grib2_var), INTENT(IN) :: grib2
    LOGICAL, INTENT(IN), OPTIONAL :: loutput, lcontainer, lrestart, &
      & lrestart_cont, initval, resetval, lmiss, missval, in_group(:), &
      & lopenacc
    INTEGER, INTENT(IN), OPTIONAL :: isteptype, l_pp_scheduler_task, &
      & tlev_source, var_class
    TYPE(t_var_metadata), POINTER, OPTIONAL :: info
    LOGICAL, CONTIGUOUS_TARGET, OPTIONAL :: p5(:,:,:,:,:)
    TYPE(t_vert_interp_meta), INTENT(IN), OPTIONAL :: vert_interp
    TYPE(t_hor_interp_meta), INTENT(IN), OPTIONAL :: hor_interp
    TYPE(t_var), POINTER, OPTIONAL :: new_element
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op
    TYPE(t_var_action), INTENT(IN), OPTIONAL :: action_list
    TYPE(t_var), POINTER :: element

    CALL add_var_list_element_5d(BOOL_T, this_list, varname, hgrid, vgrid, & 
      & cf, grib2, ldims, element, loutput, lcontainer, lrestart,          &
      & lrestart_cont, isteptype, lmiss, tlev_source, info, vert_interp,   &
      & hor_interp, in_group, l_pp_scheduler_task, post_op, action_list,   &
      & p5_l=p5, initval_l=initval, resetval_l=resetval, missval_l=missval,&
      & var_class=var_class, lopenacc=lopenacc, new_element=new_element)
    ptr => element%l_ptr(:,:,:,:,1)
  END SUBROUTINE add_var_list_element_l4d

  SUBROUTINE add_var_list_element_l3d(this_list, varname, ptr, hgrid, vgrid, &
    & cf, grib2, ldims, loutput, lcontainer, lrestart, lrestart_cont,     &
    & initval, isteptype, resetval, lmiss, missval, tlev_source, info,    &
    & p5, vert_interp, hor_interp, in_group, new_element,        &
    & l_pp_scheduler_task, post_op, action_list, var_class, lopenacc)
    TYPE(t_var_list_ptr), INTENT(INOUT) :: this_list
    CHARACTER(*), INTENT(IN) :: varname
    LOGICAL, POINTER, INTENT(OUT) :: ptr(:,:,:)
    INTEGER, INTENT(IN) :: hgrid, vgrid, ldims(3)
    TYPE(t_cf_var), INTENT(IN) :: cf
    TYPE(t_grib2_var), INTENT(IN) :: grib2
    LOGICAL, INTENT(IN), OPTIONAL :: loutput, lcontainer, lrestart, &
      & lrestart_cont, initval, resetval, lmiss, missval, in_group(:), &
      & lopenacc
    INTEGER, INTENT(IN), OPTIONAL :: isteptype, l_pp_scheduler_task, &
      & tlev_source, var_class
    TYPE(t_var_metadata), POINTER, OPTIONAL :: info
    LOGICAL, CONTIGUOUS_TARGET, OPTIONAL :: p5(:,:,:,:,:)
    TYPE(t_vert_interp_meta), INTENT(IN), OPTIONAL :: vert_interp
    TYPE(t_hor_interp_meta), INTENT(IN), OPTIONAL :: hor_interp
    TYPE(t_var), POINTER, OPTIONAL :: new_element
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op
    TYPE(t_var_action), INTENT(IN), OPTIONAL :: action_list
    TYPE(t_var), POINTER :: element

    CALL add_var_list_element_5d(BOOL_T, this_list, varname, hgrid, vgrid, & 
      & cf, grib2, ldims, element, loutput, lcontainer, lrestart,          &
      & lrestart_cont, isteptype, lmiss, tlev_source, info, vert_interp,   &
      & hor_interp, in_group, l_pp_scheduler_task, post_op, action_list,   &
      & p5_l=p5, initval_l=initval, resetval_l=resetval, missval_l=missval,&
      & var_class=var_class, lopenacc=lopenacc, new_element=new_element)
    ptr => element%l_ptr(:,:,:,1,1)
  END SUBROUTINE add_var_list_element_l3d

  SUBROUTINE add_var_list_element_l2d(this_list, varname, ptr, hgrid, vgrid, &
    & cf, grib2, ldims, loutput, lcontainer, lrestart, lrestart_cont,     &
    & initval, isteptype, resetval, lmiss, missval, tlev_source, info,    &
    & p5, vert_interp, hor_interp, in_group, new_element,        &
    & l_pp_scheduler_task, post_op, action_list, var_class, lopenacc)
    TYPE(t_var_list_ptr), INTENT(INOUT) :: this_list
    CHARACTER(*), INTENT(IN) :: varname
    LOGICAL, POINTER, INTENT(OUT) :: ptr(:,:)
    INTEGER, INTENT(IN) :: hgrid, vgrid, ldims(2)
    TYPE(t_cf_var), INTENT(IN) :: cf
    TYPE(t_grib2_var), INTENT(IN) :: grib2
    LOGICAL, INTENT(IN), OPTIONAL :: loutput, lcontainer, lrestart, &
      & lrestart_cont, initval, resetval, lmiss, missval, in_group(:), &
      & lopenacc
    INTEGER, INTENT(IN), OPTIONAL :: isteptype, l_pp_scheduler_task, &
      & tlev_source, var_class
    TYPE(t_var_metadata), POINTER, OPTIONAL :: info
    LOGICAL, CONTIGUOUS_TARGET, OPTIONAL :: p5(:,:,:,:,:)
    TYPE(t_vert_interp_meta), INTENT(IN), OPTIONAL :: vert_interp
    TYPE(t_hor_interp_meta), INTENT(IN), OPTIONAL :: hor_interp
    TYPE(t_var), POINTER, OPTIONAL :: new_element
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op
    TYPE(t_var_action), INTENT(IN), OPTIONAL :: action_list
    TYPE(t_var), POINTER :: element

    CALL add_var_list_element_5d(BOOL_T, this_list, varname, hgrid, vgrid, & 
      & cf, grib2, ldims, element, loutput, lcontainer, lrestart,          &
      & lrestart_cont, isteptype, lmiss, tlev_source, info, vert_interp,   &
      & hor_interp, in_group, l_pp_scheduler_task, post_op, action_list,   &
      & p5_l=p5, initval_l=initval, resetval_l=resetval, missval_l=missval,&
      & var_class=var_class, lopenacc=lopenacc, new_element=new_element)
    ptr => element%l_ptr(:,:,1,1,1)
  END SUBROUTINE add_var_list_element_l2d

  SUBROUTINE add_var_list_element_l1d(this_list, varname, ptr, hgrid, vgrid, &
    & cf, grib2, ldims, loutput, lcontainer, lrestart, lrestart_cont,     &
    & initval, isteptype, resetval, lmiss, missval, tlev_source, info,    &
    & p5, vert_interp, hor_interp, in_group, new_element,        &
    & l_pp_scheduler_task, post_op, action_list, var_class, lopenacc)
    TYPE(t_var_list_ptr), INTENT(INOUT) :: this_list
    CHARACTER(*), INTENT(IN) :: varname
    LOGICAL, POINTER, INTENT(OUT) :: ptr(:)
    INTEGER, INTENT(IN) :: hgrid, vgrid, ldims(1)
    TYPE(t_cf_var), INTENT(IN) :: cf
    TYPE(t_grib2_var), INTENT(IN) :: grib2
    LOGICAL, INTENT(IN), OPTIONAL :: loutput, lcontainer, lrestart, &
      & lrestart_cont, initval, resetval, lmiss, missval, in_group(:), &
      & lopenacc
    INTEGER, INTENT(IN), OPTIONAL :: isteptype, l_pp_scheduler_task, &
      & tlev_source, var_class
    TYPE(t_var_metadata), POINTER, OPTIONAL :: info
    LOGICAL, CONTIGUOUS_TARGET, OPTIONAL :: p5(:,:,:,:,:)
    TYPE(t_vert_interp_meta), INTENT(IN), OPTIONAL :: vert_interp
    TYPE(t_hor_interp_meta), INTENT(IN), OPTIONAL :: hor_interp
    TYPE(t_var), POINTER, OPTIONAL :: new_element
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op
    TYPE(t_var_action), INTENT(IN), OPTIONAL :: action_list
    TYPE(t_var), POINTER :: element

    CALL add_var_list_element_5d(BOOL_T, this_list, varname, hgrid, vgrid, & 
      & cf, grib2, ldims, element, loutput, lcontainer, lrestart,          &
      & lrestart_cont, isteptype, lmiss, tlev_source, info, vert_interp,   &
      & hor_interp, in_group, l_pp_scheduler_task, post_op, action_list,   &
      & p5_l=p5, initval_l=initval, resetval_l=resetval, missval_l=missval,&
      & var_class=var_class, lopenacc=lopenacc, new_element=new_element)
    ptr => element%l_ptr(:,1,1,1,1)
  END SUBROUTINE add_var_list_element_l1d

  SUBROUTINE add_var_list_reference_util(target_element, new_list_element,      &
    & this_list, target_name, refname, hgrid, vgrid, cf, grib2, ref_idx, ldims, &
    & dtype, icontainer, loutput, lrestart, lrestart_cont, isteptype, lmiss,    &
    & tlev_source, tracer_info, info, vert_interp, hor_interp, in_group,        &
    & new_element, l_pp_scheduler_task, post_op, action_list, idx_diag,&
    & var_class, opt_var_ref_pos, initval_r, initval_s, initval_i, missval_r,   &
    & missval_s, missval_i, resetval_r, resetval_s, resetval_i, idx_tracer)
    TYPE(t_var_list_ptr), INTENT(INOUT) :: this_list
    TYPE(t_var), INTENT(OUT), POINTER :: target_element, new_list_element
    CHARACTER(*), INTENT(IN) :: target_name, refname
    INTEGER, INTENT(IN) :: hgrid, vgrid, ref_idx, ldims(:), dtype
    INTEGER, INTENT(OUT) :: icontainer
    TYPE(t_cf_var), INTENT(IN) :: cf
    TYPE(t_grib2_var), INTENT(IN) :: grib2
    LOGICAL, INTENT(IN), OPTIONAL :: loutput, lrestart, lrestart_cont, &
      & lmiss, in_group(:)
    INTEGER, INTENT(IN), OPTIONAL :: isteptype, tlev_source, var_class, &
      & l_pp_scheduler_task, opt_var_ref_pos, idx_tracer, idx_diag
    REAL(dp), INTENT(IN), OPTIONAL :: initval_r, resetval_r, missval_r
    REAL(sp), INTENT(IN), OPTIONAL :: initval_s, resetval_s, missval_s
    INTEGER, INTENT(IN), OPTIONAL :: initval_i, resetval_i, missval_i
    CLASS(t_tracer_meta), INTENT(IN), OPTIONAL :: tracer_info
    TYPE(t_var_metadata), POINTER, OPTIONAL :: info
    TYPE(t_vert_interp_meta),INTENT(IN), OPTIONAL :: vert_interp
    TYPE(t_hor_interp_meta), INTENT(IN), OPTIONAL :: hor_interp
    TYPE(t_var), POINTER, OPTIONAL :: new_element
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op
    TYPE(t_var_action), INTENT(IN), OPTIONAL :: action_list
    TYPE(t_var_metadata), POINTER :: target_info, ref_info
    TYPE(t_union_vals) :: missvalt, initvalt, resetvalt
    CHARACTER(*), PARAMETER :: routine = modname//":add_var_list_reference_util"
    INTEGER :: var_ref_pos, ndims, di(5), di3, max_ref
    LOGICAL :: in_group_new(MAX_GROUPS)

    ndims = SIZE(ldims)
    target_element => find_list_element(this_list, target_name)
    target_info => target_element%info
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
      IF (ASSOCIATED(target_element%r_ptr)) THEN
        max_ref = SIZE(target_element%r_ptr, var_ref_pos)
      ELSE IF (ASSOCIATED(target_element%s_ptr)) THEN
        max_ref = SIZE(target_element%s_ptr, var_ref_pos)
      ELSE IF (ASSOCIATED(target_element%i_ptr)) THEN
        max_ref = SIZE(target_element%i_ptr, var_ref_pos)
      ELSE IF (ASSOCIATED(target_element%l_ptr)) THEN
        max_ref = SIZE(target_element%l_ptr, var_ref_pos)
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
    ALLOCATE(new_list_element)
    IF (PRESENT(new_element)) new_element=>new_list_element
    new_list_element%ref_to => target_element
    ref_info => new_list_element%info
    CALL inherit_var_list_metadata(this_list, ref_info)
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
    CALL set_var_metadata(ref_info, ldims, name=refname, hgrid=hgrid, &
      & vgrid=vgrid, cf=cf, grib2=grib2, loutput=loutput, lrestart=lrestart,   &
      & missval=missvalt, lrestart_cont=lrestart_cont, initval=initvalt,       &
      & isteptype=isteptype, resetval=resetvalt, tlev_source=tlev_source,      &
      & vert_interp=vert_interp, hor_interp=hor_interp, post_op=post_op,       &
      & l_pp_scheduler_task=l_pp_scheduler_task, action_list=action_list,      &
      & var_class=var_class, data_type=dtype, lmiss=lmiss, in_group=in_group,  &
      & idx_tracer=idx_tracer, idx_diag=idx_diag)
    ! set dynamic metadata, i.e. polymorphic tracer metadata
    CALL set_var_metadata_dyn(new_list_element%info_dyn, tracer_info=tracer_info)
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
    CALL register_list_element(this_list, new_list_element)
  END SUBROUTINE add_var_list_reference_util

  SUBROUTINE add_var_list_reference_r3d(this_list, target_name, refname, ptr,       &
    & hgrid, vgrid, cf, grib2, ref_idx, ldims, loutput, lrestart, lrestart_cont, &
    & initval, isteptype, resetval, lmiss, missval, tlev_source, tracer_info,    &
    & info, vert_interp, hor_interp, in_group, new_element,             &
    & l_pp_scheduler_task, post_op, action_list, opt_var_ref_pos, var_class)
    TYPE(t_var_list_ptr), INTENT(INOUT) :: this_list
    CHARACTER(*), INTENT(IN) :: target_name, refname
    REAL(dp), POINTER :: ptr(:,:,:)
    INTEGER, INTENT(IN) :: hgrid, vgrid, ref_idx, ldims(3)
    TYPE(t_cf_var), INTENT(IN) :: cf
    TYPE(t_grib2_var), INTENT(IN) :: grib2
    LOGICAL, INTENT(IN), OPTIONAL :: loutput, lrestart, lrestart_cont, &
      & lmiss, in_group(:)
    INTEGER, INTENT(IN), OPTIONAL :: isteptype, tlev_source, var_class, &
      & l_pp_scheduler_task, opt_var_ref_pos
    REAL(dp), INTENT(IN), OPTIONAL :: initval, resetval, missval
    CLASS(t_tracer_meta), INTENT(IN), OPTIONAL :: tracer_info
    TYPE(t_var_metadata), POINTER, OPTIONAL :: info
    TYPE(t_vert_interp_meta),INTENT(IN), OPTIONAL :: vert_interp
    TYPE(t_hor_interp_meta), INTENT(IN), OPTIONAL :: hor_interp
    TYPE(t_var), POINTER, OPTIONAL :: new_element
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op
    TYPE(t_var_action), INTENT(IN), OPTIONAL :: action_list
    CHARACTER(*), PARAMETER :: routine = modname//"::add_var_list_reference_r3d"
    TYPE(t_var), POINTER :: target_element, new_list_element
    INTEGER :: icontainer

    CALL add_var_list_reference_util(target_element, new_list_element,     &
      & this_list, target_name, refname, hgrid, vgrid, cf, grib2, ref_idx, &
      & ldims, REAL_T, icontainer, loutput=loutput, lrestart=lrestart,      &
      & lrestart_cont=lrestart_cont, isteptype=isteptype, lmiss=lmiss,     &
      & tlev_source=tlev_source, tracer_info=tracer_info, info=info,       &
      & vert_interp=vert_interp, hor_interp=hor_interp, in_group=in_group, &
      & new_element=new_element, l_pp_scheduler_task=l_pp_scheduler_task,  &
      & post_op=post_op, action_list=action_list, var_class=var_class,     &
      & opt_var_ref_pos=opt_var_ref_pos, initval_r=initval,                &
      & missval_r=missval, resetval_r=resetval)
    IF (.NOT. ASSOCIATED(target_element%r_ptr)) &
      & CALL finish(routine, TRIM(refname)//' not created.')
    SELECT CASE(new_list_element%info%var_ref_pos)
    CASE(1)
      ptr => target_element%r_ptr(icontainer,:,:,:,1)
    CASE(2)
      ptr => target_element%r_ptr(:,icontainer,:,:,1)
    CASE(3)
      ptr => target_element%r_ptr(:,:,icontainer,:,1)
    CASE(4)
      ptr => target_element%r_ptr(:,:,:,icontainer,1)
    CASE default
      CALL finish(routine, "internal error!")
    END SELECT
    new_list_element%r_ptr => target_element%r_ptr
    IF (.NOT. ASSOCIATED(new_list_element%r_ptr)) &
      & WRITE (0,*) 'problem with association of ptr for '//TRIM(refname)
    IF (PRESENT(lmiss)) THEN
      ptr = new_list_element%info%missval%rval
    ELSE
      ptr = 0.0_dp
    END IF
  END SUBROUTINE add_var_list_reference_r3d

  SUBROUTINE add_var_list_reference_r2d(this_list, target_name, refname, ptr,    &
    & hgrid, vgrid, cf, grib2, ref_idx, ldims, loutput, lrestart, lrestart_cont, &
    & initval, isteptype, resetval, lmiss, missval, tlev_source, tracer_info,    &
    & info, vert_interp, hor_interp, in_group, new_element,             &
    & l_pp_scheduler_task, post_op, action_list, opt_var_ref_pos, var_class,     &
    & idx_tracer, idx_diag)
    TYPE(t_var_list_ptr), INTENT(INOUT) :: this_list
    CHARACTER(*), INTENT(IN) :: target_name, refname
    REAL(dp), POINTER :: ptr(:,:)
    INTEGER, INTENT(IN) :: hgrid, vgrid, ref_idx, ldims(2)
    TYPE(t_cf_var), INTENT(IN) :: cf                  
    TYPE(t_grib2_var), INTENT(IN) :: grib2               
    LOGICAL, INTENT(IN), OPTIONAL :: loutput, lrestart, lrestart_cont, &
      & lmiss, in_group(:)
    INTEGER, INTENT(IN), OPTIONAL :: isteptype, tlev_source, var_class, &
      & l_pp_scheduler_task, opt_var_ref_pos, idx_tracer, idx_diag
    REAL(dp), INTENT(IN), OPTIONAL :: initval, resetval, missval
    CLASS(t_tracer_meta), INTENT(IN), OPTIONAL :: tracer_info         
    TYPE(t_var_metadata), POINTER, OPTIONAL :: info                
    TYPE(t_vert_interp_meta),INTENT(IN), OPTIONAL :: vert_interp         
    TYPE(t_hor_interp_meta), INTENT(IN), OPTIONAL :: hor_interp          
    TYPE(t_var), POINTER, OPTIONAL :: new_element         
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op            
    TYPE(t_var_action), INTENT(IN), OPTIONAL :: action_list         
    CHARACTER(*), PARAMETER :: routine = modname//"::add_var_list_reference_r2d"
    TYPE(t_var), POINTER :: target_element, new_list_element
    INTEGER :: icontainer

    CALL add_var_list_reference_util(target_element, new_list_element,     &
      & this_list, target_name, refname, hgrid, vgrid, cf, grib2, ref_idx, &
      & ldims, REAL_T, icontainer, loutput=loutput, lrestart=lrestart,      &
      & lrestart_cont=lrestart_cont, isteptype=isteptype, lmiss=lmiss,     &
      & tlev_source=tlev_source, tracer_info=tracer_info, info=info,       &
      & vert_interp=vert_interp, hor_interp=hor_interp, in_group=in_group, &
      & new_element=new_element, l_pp_scheduler_task=l_pp_scheduler_task,  &
      & post_op=post_op, action_list=action_list, var_class=var_class,     &
      & opt_var_ref_pos=opt_var_ref_pos, initval_r=initval,                &
      & missval_r=missval, resetval_r=resetval, idx_tracer=idx_tracer,     &
      & idx_diag=idx_diag)
    IF (.NOT. ASSOCIATED(target_element%r_ptr)) &
      & CALL finish(routine, TRIM(refname)//' not created.')
    SELECT CASE(new_list_element%info%var_ref_pos)
    CASE(1)
      ptr => target_element%r_ptr(icontainer,:,:,1,1)
    CASE(2)
      ptr => target_element%r_ptr(:,icontainer,:,1,1)
    CASE(3)
      ptr => target_element%r_ptr(:,:,icontainer,1,1)
    CASE default
      CALL finish(routine, "internal error!")
    END SELECT
    new_list_element%r_ptr => target_element%r_ptr
    IF (.NOT. ASSOCIATED(new_list_element%r_ptr)) &
      & WRITE (0,*) 'problem with association of ptr for '//TRIM(refname)
    IF (PRESENT(lmiss)) THEN
      ptr = new_list_element%info%missval%rval
    ELSE
      ptr = 0.0_dp
    END IF
  END SUBROUTINE add_var_list_reference_r2d

  SUBROUTINE add_var_list_reference_s3d(this_list, target_name, refname, ptr,    &
    & hgrid, vgrid, cf, grib2, ref_idx, ldims, loutput, lrestart, lrestart_cont, &
    & initval, isteptype, resetval, lmiss, missval, tlev_source, tracer_info,    &
    & info, vert_interp, hor_interp, in_group, new_element,             &
    & l_pp_scheduler_task, post_op, action_list, opt_var_ref_pos, var_class)
    TYPE(t_var_list_ptr), INTENT(INOUT) :: this_list
    CHARACTER(*), INTENT(IN) :: target_name, refname
    REAL(sp), POINTER :: ptr(:,:,:)
    INTEGER, INTENT(IN) :: hgrid, vgrid, ref_idx, ldims(3)
    TYPE(t_cf_var), INTENT(IN) :: cf
    TYPE(t_grib2_var), INTENT(IN) :: grib2
    LOGICAL, INTENT(IN), OPTIONAL :: loutput, lrestart, lrestart_cont, &
      & lmiss, in_group(:)
    INTEGER, INTENT(IN), OPTIONAL :: isteptype, tlev_source, var_class, &
      & l_pp_scheduler_task, opt_var_ref_pos
    REAL(sp), INTENT(IN), OPTIONAL :: initval, resetval, missval
    CLASS(t_tracer_meta), INTENT(IN), OPTIONAL :: tracer_info
    TYPE(t_var_metadata), POINTER, OPTIONAL :: info
    TYPE(t_vert_interp_meta),INTENT(IN), OPTIONAL :: vert_interp
    TYPE(t_hor_interp_meta), INTENT(IN), OPTIONAL :: hor_interp
    TYPE(t_var), POINTER, OPTIONAL :: new_element
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op
    TYPE(t_var_action), INTENT(IN), OPTIONAL :: action_list
    CHARACTER(*), PARAMETER :: routine = modname//"::add_var_list_reference_s3d"
    TYPE(t_var), POINTER :: target_element, new_list_element
    INTEGER :: icontainer

    CALL add_var_list_reference_util(target_element, new_list_element,     &
      & this_list, target_name, refname, hgrid, vgrid, cf, grib2, ref_idx, &
      & ldims, SINGLE_T, icontainer, loutput=loutput, lrestart=lrestart,   &
      & lrestart_cont=lrestart_cont, isteptype=isteptype, lmiss=lmiss,     &
      & tlev_source=tlev_source, tracer_info=tracer_info, info=info,       &
      & vert_interp=vert_interp, hor_interp=hor_interp, in_group=in_group, &
      & new_element=new_element, l_pp_scheduler_task=l_pp_scheduler_task,  &
      & post_op=post_op, action_list=action_list, var_class=var_class,     &
      & opt_var_ref_pos=opt_var_ref_pos, initval_s=initval,                &
      & missval_s=missval, resetval_s=resetval)
    IF (.NOT. ASSOCIATED(target_element%s_ptr)) &
      & CALL finish(routine, TRIM(refname)//' not created.')
    SELECT CASE(new_list_element%info%var_ref_pos)
    CASE(1)
      ptr => target_element%s_ptr(icontainer,:,:,:,1)
    CASE(2)
      ptr => target_element%s_ptr(:,icontainer,:,:,1)
    CASE(3)
      ptr => target_element%s_ptr(:,:,icontainer,:,1)
    CASE(4)
      ptr => target_element%s_ptr(:,:,:,icontainer,1)
    CASE default
      CALL finish(routine, "internal error!")
    END SELECT
    new_list_element%s_ptr => target_element%s_ptr
    IF (.NOT. ASSOCIATED(new_list_element%s_ptr)) &
      & WRITE (0,*) 'problem with association of ptr for '//TRIM(refname)
    IF (PRESENT(lmiss)) THEN
      ptr = new_list_element%info%missval%sval
    ELSE
      ptr = 0.0_sp
    END IF
  END SUBROUTINE add_var_list_reference_s3d

  SUBROUTINE add_var_list_reference_s2d(this_list, target_name, refname, ptr,    &
    & hgrid, vgrid, cf, grib2, ref_idx, ldims, loutput, lrestart, lrestart_cont, &
    & initval, isteptype, resetval, lmiss, missval, tlev_source, tracer_info,    &
    & info, vert_interp, hor_interp, in_group, new_element,             &
    & l_pp_scheduler_task, post_op, action_list, opt_var_ref_pos, var_class)
    TYPE(t_var_list_ptr), INTENT(INOUT) :: this_list
    CHARACTER(*), INTENT(IN) :: target_name, refname
    REAL(sp), POINTER :: ptr(:,:)
    INTEGER, INTENT(IN) :: hgrid, vgrid, ref_idx, ldims(2)
    TYPE(t_cf_var), INTENT(IN) :: cf
    TYPE(t_grib2_var), INTENT(IN) :: grib2
    LOGICAL, INTENT(IN), OPTIONAL :: loutput, lrestart, lrestart_cont, &
      & lmiss, in_group(:)
    INTEGER, INTENT(IN), OPTIONAL :: isteptype, tlev_source, var_class, &
      & l_pp_scheduler_task, opt_var_ref_pos
    REAL(sp), INTENT(IN), OPTIONAL :: initval, resetval, missval
    CLASS(t_tracer_meta), INTENT(IN), OPTIONAL :: tracer_info
    TYPE(t_var_metadata), POINTER, OPTIONAL :: info
    TYPE(t_vert_interp_meta),INTENT(IN), OPTIONAL :: vert_interp
    TYPE(t_hor_interp_meta), INTENT(IN), OPTIONAL :: hor_interp
    TYPE(t_var), POINTER, OPTIONAL :: new_element
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op
    TYPE(t_var_action), INTENT(IN), OPTIONAL :: action_list
    CHARACTER(*), PARAMETER :: routine = modname//"::add_var_list_reference_s2d"
    TYPE(t_var), POINTER :: target_element, new_list_element
    INTEGER :: icontainer

    CALL add_var_list_reference_util(target_element, new_list_element,     &
      & this_list, target_name, refname, hgrid, vgrid, cf, grib2, ref_idx, &
      & ldims, SINGLE_T, icontainer, loutput=loutput, lrestart=lrestart,   &
      & lrestart_cont=lrestart_cont, isteptype=isteptype, lmiss=lmiss,     &
      & tlev_source=tlev_source, tracer_info=tracer_info, info=info,       &
      & vert_interp=vert_interp, hor_interp=hor_interp, in_group=in_group, &
      & new_element=new_element, l_pp_scheduler_task=l_pp_scheduler_task,  &
      & post_op=post_op, action_list=action_list, var_class=var_class,     &
      & opt_var_ref_pos=opt_var_ref_pos, initval_s=initval,                &
      & missval_s=missval, resetval_s=resetval)
    IF (.NOT. ASSOCIATED(target_element%s_ptr)) &
      & CALL finish(routine, TRIM(refname)//' not created.')
    SELECT CASE(new_list_element%info%var_ref_pos)
    CASE(1)
      ptr => target_element%s_ptr(icontainer,:,:,1,1)
    CASE(2)
      ptr => target_element%s_ptr(:,icontainer,:,1,1)
    CASE(3)
      ptr => target_element%s_ptr(:,:,icontainer,1,1)
    CASE default
      CALL finish(routine, "internal error!")
    END SELECT
    new_list_element%s_ptr => target_element%s_ptr
    IF (.NOT. ASSOCIATED(new_list_element%s_ptr)) &
      & WRITE (0,*) 'problem with association of ptr for '//TRIM(refname)
    IF (PRESENT(lmiss)) THEN
      ptr = new_list_element%info%missval%sval
    ELSE
      ptr = 0.0_sp
    END IF
  END SUBROUTINE add_var_list_reference_s2d

  SUBROUTINE add_var_list_reference_i2d(this_list, target_name, refname, ptr,    &
    & hgrid, vgrid, cf, grib2, ref_idx, ldims, loutput, lrestart, lrestart_cont, &
    & initval, isteptype, resetval, lmiss, missval, tlev_source, tracer_info,    &
    & info, vert_interp, hor_interp, in_group, new_element,             &
    & l_pp_scheduler_task, post_op, action_list, opt_var_ref_pos, var_class)
    TYPE(t_var_list_ptr), INTENT(INOUT) :: this_list
    CHARACTER(*), INTENT(IN) :: target_name, refname
    INTEGER, POINTER :: ptr(:,:)
    INTEGER, INTENT(IN) :: hgrid, vgrid, ref_idx, ldims(2)
    TYPE(t_cf_var), INTENT(IN) :: cf
    TYPE(t_grib2_var), INTENT(IN) :: grib2
    LOGICAL, INTENT(IN), OPTIONAL :: loutput, lrestart, lrestart_cont, &
      & lmiss, in_group(:)
    INTEGER, INTENT(IN), OPTIONAL :: isteptype, tlev_source, var_class, &
      & l_pp_scheduler_task, opt_var_ref_pos, initval, resetval, missval
    CLASS(t_tracer_meta), INTENT(IN), OPTIONAL :: tracer_info
    TYPE(t_var_metadata), POINTER, OPTIONAL :: info
    TYPE(t_vert_interp_meta),INTENT(IN), OPTIONAL :: vert_interp
    TYPE(t_hor_interp_meta), INTENT(IN), OPTIONAL :: hor_interp
    TYPE(t_var), POINTER, OPTIONAL :: new_element
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op
    TYPE(t_var_action), INTENT(IN), OPTIONAL :: action_list
    CHARACTER(*), PARAMETER :: routine = modname//"::add_var_list_reference_i2d"
    TYPE(t_var), POINTER :: target_element, new_list_element
    INTEGER :: icontainer

    CALL add_var_list_reference_util(target_element, new_list_element,     &
      & this_list, target_name, refname, hgrid, vgrid, cf, grib2, ref_idx, &
      & ldims, INT_T, icontainer, loutput=loutput, lrestart=lrestart,      &
      & lrestart_cont=lrestart_cont, isteptype=isteptype, lmiss=lmiss,     &
      & tlev_source=tlev_source, tracer_info=tracer_info, info=info,       &
      & vert_interp=vert_interp, hor_interp=hor_interp, in_group=in_group, &
      & new_element=new_element, l_pp_scheduler_task=l_pp_scheduler_task,  &
      & post_op=post_op, action_list=action_list, var_class=var_class,     &
      & opt_var_ref_pos=opt_var_ref_pos, initval_i=initval,                &
      & missval_i=missval, resetval_i=resetval)
    IF (.NOT. ASSOCIATED(target_element%i_ptr)) &
      & CALL finish(routine, TRIM(refname)//' not created.')
    SELECT CASE(new_list_element%info%var_ref_pos)
    CASE(1)
      ptr => target_element%i_ptr(icontainer,:,:,1,1)
    CASE(2)
      ptr => target_element%i_ptr(:,icontainer,:,1,1)
    CASE(3)
      ptr => target_element%i_ptr(:,:,icontainer,1,1)
    CASE default
      CALL finish(routine, "internal error!")
    END SELECT
    new_list_element%i_ptr => target_element%i_ptr
    IF (.NOT. ASSOCIATED(new_list_element%i_ptr)) &
      & WRITE (0,*) 'problem with association of ptr for '//TRIM(refname)
    IF (PRESENT(lmiss)) THEN
      ptr = new_list_element%info%missval%ival
    ELSE
      ptr = 0
    END IF
  END SUBROUTINE add_var_list_reference_i2d

  SUBROUTINE print_var_list(this, lshort)
    CLASS(t_var_list_ptr), INTENT(in) :: this
    LOGICAL, OPTIONAL :: lshort
    TYPE(t_var_metadata), POINTER :: info
    TYPE(t_var), POINTER :: le
    CHARACTER(len=32) :: dimension_text, dtext,keytext
    INTEGER :: i, igrp, ivintp_type, ndims, j
    CHARACTER(len=4) :: localMode
    LOGICAL :: short
    LOGICAL, POINTER :: in_group(:)

    short = .FALSE.
    IF (PRESENT(lshort)) short = lshort
    CALL message('','')
    CALL message('','')
    CALL message('','Status of variable list '//TRIM(this%p%name)//':')
    CALL message('','')
    DO j = 1, this%p%nvars
      le => this%p%vl(j)%p
      info => le%info
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
          message_text = 'Table entry name                            : '//TRIM(info%name)
          CALL message('', message_text)
          WRITE (keytext,'(i32.1)') this%p%key(j)
          WRITE (message_text,'(a,a)')       &
               'Key entry                                   : ', TRIM(keytext)
          CALL message('', message_text)
          IF (ASSOCIATED(le%r_ptr) .OR. ASSOCIATED(le%s_ptr) .OR. &
            & ASSOCIATED(le%i_ptr) .OR. ASSOCIATED(le%l_ptr)) THEN
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
            IF (ASSOCIATED(le%r_ptr)) THEN
              WRITE (message_text,'(a,e20.12)')      &
                & 'Missing value                               : ', info%missval%rval
            ELSE IF (ASSOCIATED(le%s_ptr)) THEN
              WRITE (message_text,'(a,e20.12)')      &
                & 'Missing value                               : ', info%missval%sval
            ELSE IF (ASSOCIATED(le%i_ptr)) THEN
              WRITE (message_text,'(a,i8)')      &
                & 'Missing value                               : ', info%missval%ival
            ELSE IF (ASSOCIATED(le%l_ptr)) THEN
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
            & MERGE("yes.", " no.", le%info_dyn%tracer%lis_tracer))
          IF (le%info_dyn%tracer%lis_tracer) THEN
            CALL message('', 'Child-to-parent feedback                  : ' // &
              & MERGE("yes.", " no.", le%info_dyn%tracer%lfeedback))
            WRITE (message_text,'(a,3i3)') 'Horizontal transport method                 : ', &
              & le%info_dyn%tracer%ihadv_tracer
            CALL message('', message_text)
            WRITE (message_text,'(a,3i3)') 'Vertical transport method                   : ', &
              & le%info_dyn%tracer%ivadv_tracer
            CALL message('', message_text)
            CALL message('', 'Turbulent transport                         : ' // &
              & MERGE("yes.", " no.", le%info_dyn%tracer%lturb_tracer))
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
    ENDDO
  END SUBROUTINE print_var_list

  !-----------------------------------------------------------------------------
  ! Should be overloaded to be able to search for the different information 
  ! In the proposed structure for the linked list, in the example only
  ! A character string is used so it is straight forward only one find
  FUNCTION find_list_element(this, vname, opt_hgrid, opt_with_tl, opt_output) RESULT(element)
    TYPE(t_var_list_ptr), INTENT(IN) :: this
    CHARACTER(*), INTENT(IN) :: vname
    INTEGER, OPTIONAL, INTENT(IN) :: opt_hgrid
    LOGICAL, OPTIONAL, INTENT(IN) :: opt_with_tl, opt_output
    TYPE(t_var), POINTER :: element
    INTEGER :: key, hgrid, time_lev, iv
    LOGICAL :: with_tl, omit_output, with_output

    NULLIFY(element)
    with_tl = .TRUE.
    IF (PRESENT(opt_with_tl)) with_tl = opt_with_tl
    hgrid = -1
    IF (PRESENT(opt_hgrid)) hgrid = opt_hgrid
    IF (with_tl) THEN
      key = text_hash_c(TRIM(vname))
    ELSE
      key = text_hash_c(tolower(vname))
    END IF
    with_output = .TRUE.
    omit_output = .NOT.PRESENT(opt_output)
    IF (.NOT.omit_output) with_output = opt_output
    time_lev = get_var_timelevel(vname)
    DO iv = 1, this%p%nvars
      IF (-1 .NE. hgrid .AND. this%p%hgrid(iv) .NE. hgrid) CYCLE
      IF (.NOT.MERGE(.TRUE., with_output .EQV. this%p%lout(iv), omit_output)) &
        CYCLE
      IF (time_lev .NE. this%p%tl(iv)) CYCLE
      IF (key .NE. MERGE(this%p%key(iv), this%p%key_notl(iv), with_tl)) CYCLE
      element => this%p%vl(iv)%p
    ENDDO
  END FUNCTION find_list_element
  
END MODULE mo_var_list
