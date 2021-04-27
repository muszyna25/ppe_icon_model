! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
#include "icon_contiguous_defines.h"
MODULE mo_var_list

#if defined (__INTEL_COMPILER) || defined (__PGI) || defined (NAGFOR)
#ifdef VARLIST_INITIZIALIZE_WITH_NAN
  USE, INTRINSIC :: ieee_features
  USE, INTRINSIC :: ieee_arithmetic
  USE, INTRINSIC :: ieee_exceptions
#endif
#endif

  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_null_char, C_SIZE_T

  USE mo_kind,             ONLY: sp, wp, i8
  USE mo_cdi,              ONLY: TSTEP_INSTANT,                     &
       &                         CDI_UNDEFID
  USE mo_mpi,              ONLY: my_process_is_stdio, p_bcast
  USE mo_cf_convention,    ONLY: t_cf_var
  USE mo_grib2,            ONLY: t_grib2_var, grib2_var
  USE mo_run_config,       ONLY: msg_level
  USE mo_var_groups,       ONLY: var_groups_dyn, groups
  USE mo_var_metadata_types,ONLY: t_var_metadata, t_union_vals,     &
    &                            t_var_metadata_dynamic,            &
    &                            t_vert_interp_meta,                &
    &                            t_hor_interp_meta,                 &
    &                            MAX_GROUPS, VINTP_TYPE_LIST,       &
    &                            t_post_op_meta,                    &
    &                            CLASS_DEFAULT, CLASS_TILE,         &
    &                            CLASS_TILE_LAND
  USE mo_var_metadata,     ONLY: create_vert_interp_metadata,       &
    &                            create_hor_interp_metadata,        &
    &                            post_op, actions
  USE mo_var_groups,       ONLY: groups
  USE mo_tracer_metadata,  ONLY: create_tracer_metadata
  USE mo_tracer_metadata_types,ONLY: t_tracer_meta
  USE mo_var_list_element, ONLY: t_var_list_element, t_p_var_list_element, &
    &                            level_type_ml
  USE mo_linked_list,      ONLY: t_var_list, t_list_element,        &
       &                         new_list, delete_list,             &
       &                         append_list_element,               &
       &                         t_var_list_intrinsic
  USE mo_exception,        ONLY: message, message_text, finish
  USE mo_util_hash,        ONLY: util_hashword
  USE mo_util_string,      ONLY: remove_duplicates, toupper,        &
    &                            pretty_print_string_list, tolower, &
    &                            difference, find_trailing_number,  &
    &                            lowcase
  USE mo_impl_constants,   ONLY: max_var_lists, vname_len,          &
    &                            max_var_list_name_len,             &
    &                            STR_HINTP_TYPE, MAX_TIME_LEVELS,   &
    &                            TLEV_NNOW, REAL_T, SINGLE_T,       &
    &                            BOOL_T, INT_T, SUCCESS,            &
    &                            VARNAME_LEN,                       &
    &                            TIMELEVEL_SUFFIX
  USE mo_cdi_constants,    ONLY: GRID_UNSTRUCTURED_CELL,            &
    &                            GRID_REGULAR_LONLAT
  USE mo_fortran_tools,    ONLY: assign_if_present, &
    &                            init_contiguous_dp, init_contiguous_sp, &
    &                            init_contiguous_i4, init_contiguous_l
  USE mo_action_types,     ONLY: t_var_action
  USE mo_io_config,        ONLY: restart_file_type
  USE mo_packed_message,   ONLY: t_PackedMessage, kPackOp, kUnpackOp
  USE mo_util_sort,        ONLY: quicksort
#ifdef DEBUG_MVSTREAM
  USE mo_mpi,              ONLY: my_process_is_stdio
  USE mo_util_string,      ONLY: int2string
  USE self_assert,         ONLY: print_summary
#endif

  IMPLICIT NONE

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_var_list'

  PRIVATE

  PUBLIC :: new_var_list              ! get a pointer to a new output var_list
  PUBLIC :: delete_var_list           ! delete an output var_list
  PUBLIC :: delete_var_lists          ! delete all output var_lists
  PUBLIC :: get_var_list              ! get a pointer to an existing output var_list
  PUBLIC :: set_var_list              ! set default parameters of an output var_list
  PUBLIC :: print_var
  PUBLIC :: print_var_list
  PUBLIC :: print_all_var_lists
  PUBLIC :: print_memory_use

  PUBLIC :: default_var_list_settings ! set default settings for a whole list
  PUBLIC :: collect_group

  !PUBLIC :: var_lists                 ! vector of output var_lists
  PUBLIC :: nvar_lists                ! number of output var_lists defined so far
  PUBLIC :: max_var_lists

  PUBLIC :: add_var                   ! create/allocate a new var_list list entry
  PUBLIC :: add_var_list_reference
  PUBLIC :: add_ref                   ! create/reference a new var_list list entry
  PUBLIC :: get_var                   ! obtain reference to existing list entry
  PUBLIC :: get_all_var_names         ! obtain a list of variables names

  PUBLIC :: get_var_name              ! return plain variable name (without timelevel)
  PUBLIC :: get_var_timelevel         ! return variable timelevel (or "-1")
  PUBLIC :: get_var_tileidx           ! return variable tile index
  PUBLIC :: get_var_container
  PUBLIC :: get_var_list_element_info ! return a copy of the metadata for a var_list element
  PUBLIC :: get_tracer_info_dyn_by_idx! return a copy of the dynamic metadata of a certain tracer
  PUBLIC :: get_timelevel_string      ! return the default string with timelevel encoded
  PUBLIC :: get_varname_with_timelevel! join varname with timelevel string

  PUBLIC :: total_number_of_variables ! returns total number of defined variables
  PUBLIC :: get_model_types
  PUBLIC :: get_restart_vars
  PUBLIC :: var_lists_apply
  PUBLIC :: fget_var_list_element_r1d
  PUBLIC :: fget_var_list_element_r2d
  PUBLIC :: fget_var_list_element_r3d
  PUBLIC :: fget_var_list_element_s1d
  PUBLIC :: fget_var_list_element_s2d
  PUBLIC :: fget_var_list_element_s3d

  PUBLIC :: find_list_element   ! find an element in the list
  PUBLIC :: find_element   ! find an element in the list

  PUBLIC :: varlistPacker

  PUBLIC :: print_group_details

  INTERFACE find_element
    MODULE PROCEDURE find_list_element
    MODULE PROCEDURE find_element_from_all
    MODULE PROCEDURE find_field_from_selected
    MODULE PROCEDURE find_field_from_selected2
  END INTERFACE find_element

  INTERFACE print_var_list
    MODULE PROCEDURE print_var_list
    MODULE PROCEDURE print_var_list_intrinsic
  END INTERFACE print_var_list

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

  INTERFACE add_var_list_reference
    MODULE PROCEDURE add_var_list_reference
    MODULE PROCEDURE add_var_list_reference_named
  END INTERFACE add_var_list_reference

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

  INTERFACE struct_assign_if_present  ! purely internal
    MODULE PROCEDURE assign_if_present_cf
    MODULE PROCEDURE assign_if_present_grib2
    MODULE PROCEDURE assign_if_present_union
    MODULE PROCEDURE assign_if_present_tracer_meta
    MODULE PROCEDURE assign_if_present_vert_interp
    MODULE PROCEDURE assign_if_present_hor_interp
    MODULE PROCEDURE assign_if_present_post_op
    MODULE PROCEDURE assign_if_present_action_list
  END INTERFACE struct_assign_if_present

  ABSTRACT INTERFACE
    FUNCTION list_select(var_list, state) RESULT(is_selected)
      IMPORT :: t_var_list_intrinsic
      LOGICAL :: is_selected
      TYPE(t_var_list_intrinsic), INTENT(in) :: var_list
      CLASS(*), TARGET :: state
    END FUNCTION list_select
    FUNCTION var_filter(field, state, var_list) RESULT(is_selected)
      IMPORT :: t_var_list_element, t_var_list_intrinsic
      LOGICAL :: is_selected
      TYPE(t_var_list_element), INTENT(in) :: field
      CLASS(*), TARGET :: state
      TYPE(t_var_list_intrinsic), INTENT(in) :: var_list
    END FUNCTION var_filter
    SUBROUTINE var_apply(field, state, var_list)
      IMPORT :: t_var_list_element, t_var_list_intrinsic
      TYPE(t_var_list_element), TARGET :: field
      CLASS(*), TARGET :: state
      TYPE(t_var_list_intrinsic), INTENT(in) :: var_list
    END SUBROUTINE var_apply
  END INTERFACE

  INTEGER,                  SAVE :: nvar_lists     =   0      ! var_lists allocated so far
  !
  TYPE(t_var_list), TARGET, SAVE :: var_lists(max_var_lists)  ! memory buffer array
  !
#ifndef NOMPI
  PUBLIC :: replicate_var_lists
#endif
CONTAINS
  !------------------------------------------------------------------------------------------------
  !
  ! Create a new memory buffer / output var_list
  ! Get a pointer to the new var_list
  !
  SUBROUTINE new_var_list (this_list, name, output_type, restart_type,      &
       &                   post_suf, rest_suf, init_suf, loutput, lrestart, &
       &                   linitial, patch_id, vlevel_type)
    !
    TYPE(t_var_list), INTENT(inout)        :: this_list    ! anchor
    CHARACTER(len=*), INTENT(in)           :: name         ! name of output var_list
    INTEGER,          INTENT(in), OPTIONAL :: output_type  ! 'GRIB2' or 'NetCDF'
    INTEGER,          INTENT(in), OPTIONAL :: restart_type ! 'GRIB2' or 'NetCDF'
    CHARACTER(len=*), INTENT(in), OPTIONAL :: post_suf     ! suffix of output file
    CHARACTER(len=*), INTENT(in), OPTIONAL :: rest_suf     ! suffix of restart file
    CHARACTER(len=*), INTENT(in), OPTIONAL :: init_suf     ! suffix of initial file
    LOGICAL,          INTENT(in), OPTIONAL :: loutput      ! write to  output file
    LOGICAL,          INTENT(in), OPTIONAL :: lrestart     ! write to restart file
    LOGICAL,          INTENT(in), OPTIONAL :: linitial     ! read from initial file
    INTEGER,          INTENT(in), OPTIONAL :: patch_id     ! patch ID
    INTEGER,          INTENT(in), OPTIONAL :: vlevel_type  ! 1/2/3 for model/pres./height levels
    !
    INTEGER :: i
    !
    ! look, if name exists already in list
    !
    DO i = 1, nvar_lists
      IF (var_lists(i)%p%name == name) THEN
        CALL finish('new_list', 'output var_list '//TRIM(name)//' already used.')
      ENDIF
    ENDDO
    !
    this_list%p => NULL()
    !
    ! - check, if there is an entry without name in the existing vector
    !
    DO i = 1, nvar_lists
      IF (var_lists(i)%p%name == '') THEN
        this_list%p => var_lists(i)%p
        EXIT
      ENDIF
    END DO
    !
    ! - if not successful, append to vector of lists
    !
    IF(.NOT. ASSOCIATED(this_list%p)) THEN
      nvar_lists = nvar_lists + 1
      IF (nvar_lists > max_var_lists) THEN
        CALL finish('new_list', &
             &      'var_lists container overflow, increase "max_var_lists" in mo_var_list.f90')
      ENDIF
    ENDIF
    !
    CALL new_list (var_lists(nvar_lists))
    !
    ! connect anchor and backbone by referencing
    !
    this_list%p => var_lists(nvar_lists)%p
    !
    ! set default list characteristics
    !
    this_list%p%name     = name
    this_list%p%post_suf = '_'//name
    this_list%p%rest_suf = this_list%p%post_suf
    this_list%p%init_suf = this_list%p%post_suf
    this_list%p%loutput  = .TRUE.
    !
    ! set non-default list characteristics
    !
    this_list%p%restart_type = restart_file_type

    CALL assign_if_present(this_list%p%output_type,  output_type)
    CALL assign_if_present(this_list%p%restart_type, restart_type)
    CALL assign_if_present(this_list%p%post_suf,     post_suf)
    CALL assign_if_present(this_list%p%rest_suf,     rest_suf)
    CALL assign_if_present(this_list%p%init_suf,     init_suf)
    CALL assign_if_present(this_list%p%loutput,      loutput)
    CALL assign_if_present(this_list%p%lrestart,     lrestart)
    CALL assign_if_present(this_list%p%linitial,     linitial)
    CALL assign_if_present(this_list%p%patch_id,     patch_id)
    CALL assign_if_present(this_list%p%vlevel_type,  vlevel_type)
    !
    CALL message('','')
    CALL message('','adding new var_list '//TRIM(name))
    !
  END SUBROUTINE new_var_list
  !------------------------------------------------------------------------------------------------
  !
  ! Get a reference to a memory buffer/output var_list
  !
  SUBROUTINE get_var_list (this_list, name)
    !
    TYPE(t_var_list), INTENT(out), POINTER    :: this_list ! pointer
    CHARACTER(len=*), INTENT(in) :: name      ! name of output var_list
    !
    INTEGER :: i
    !
    NULLIFY (this_list)
    !
    DO i = 1, nvar_lists
      IF (var_lists(i)%p%name == name) THEN
        this_list => var_lists(i)
        EXIT
      ENDIF
    END DO
    !
  END SUBROUTINE get_var_list

  !------------------------------------------------------------------------------------------------
  !
  ! @return total number of (non-container) variables
  !
  FUNCTION total_number_of_variables(opt_list_select, opt_state, opt_var_select)
    INTEGER :: total_number_of_variables
    PROCEDURE(list_select), OPTIONAL :: opt_list_select
    PROCEDURE(var_filter), OPTIONAL :: opt_var_select
    CLASS(*), OPTIONAL :: opt_state
    ! local variables
    INTEGER :: i
    TYPE(t_list_element), POINTER :: element

    total_number_of_variables = 0
    !- loop over variables

    ! Note that there may be several variables with different time
    ! levels, we just add unconditionally all
    DO i = 1,nvar_lists
      IF (PRESENT(opt_list_select)) THEN
        IF (.NOT. opt_list_select(var_lists(i)%p, opt_state)) CYCLE
      END IF
      element => var_lists(i)%p%first_list_element
      IF (.NOT. PRESENT(opt_var_select)) THEN
        LOOPVAR : DO WHILE (ASSOCIATED(element))
          ! Do not count element if it is a container
          total_number_of_variables = total_number_of_variables &
               + MERGE(1, 0, .NOT. element%field%info%lcontainer)
          element => element%next_list_element
        ENDDO LOOPVAR ! loop over vlist "i"
      ELSE
        LOOPVAR_CUSTOM_FILTER : DO WHILE (ASSOCIATED(element))
          ! Do not count element if it is a container
          total_number_of_variables = total_number_of_variables &
            &  + MERGE(1, 0, opt_var_select(element%field, opt_state, &
            &                               var_lists(i)%p))
          element => element%next_list_element
        ENDDO LOOPVAR_CUSTOM_FILTER ! loop over vlist "i"
      END IF
    ENDDO ! i = 1,nvar_lists
  END FUNCTION total_number_of_variables

  !> builds array of all model type names seen in any variable list
  SUBROUTINE get_model_types(model_types, patch_id)
    CHARACTER(len=8), ALLOCATABLE, INTENT(out) :: model_types(:)
    INTEGER, intent(in) :: patch_id
    INTEGER :: i, nmt, ierror
    CHARACTER(len=8) :: model_types_(nvar_lists)
    CHARACTER(LEN=*), PARAMETER :: routine = modname//":get_model_types"

    nmt = 0
    DO i = 1, nvar_lists
      IF (var_lists(i)%p%patch_id == patch_id &
          .AND. ALL(model_types_(1:nmt) /= var_lists(i)%p%model_type)) THEN
        nmt = nmt + 1
        model_types_(nmt) = var_lists(i)%p%model_type
      END IF
    END DO
    ALLOCATE(model_types(nmt), STAT=ierror)
    IF (ierror /= SUCCESS) CALL finish(routine, "memory allocation failed")
    model_types(1:nmt) = model_types_(1:nmt)
  END SUBROUTINE get_model_types

  !> used as a filter to only select variable lists that match a given
  !! patch id and model type name
  LOGICAL FUNCTION list_is_in_patch_restart(varlist, patch_id, model_type)
    TYPE(t_var_list), INTENT(IN) :: varlist
    INTEGER, INTENT(IN) :: patch_id
    CHARACTER(LEN = *), INTENT(IN) :: model_type

    list_is_in_patch_restart = varlist%p%lrestart &
      .AND. varlist%p%patch_id == patch_id &
      .AND. varlist%p%model_type == model_type
  END FUNCTION list_is_in_patch_restart

  ! compute the number of  restart variables for the given logical patch.
  FUNCTION number_of_restart_variables(patch_id, modelType) RESULT(fld_cnt)
    INTEGER :: fld_cnt
    INTEGER, INTENT(in) :: patch_id
    CHARACTER(LEN = *), INTENT(IN) :: modelType
    INTEGER :: i
    TYPE(t_list_element), POINTER :: element

    fld_cnt = 0
    DO i = 1, nvar_lists
      IF (list_is_in_patch_restart(var_lists(i), patch_id, modelType)) THEN
        ! check, if the list has valid restart fields
        element => var_lists(i)%p%first_list_element
        DO WHILE (ASSOCIATED(element))
          fld_cnt = fld_cnt + MERGE(1, 0, element%field%info%lrestart)
          element => element%next_list_element
        END DO
      END IF
    ENDDO
  END FUNCTION number_of_restart_variables

  !> collect references to all variables which match patch id and
  !! model type name into the var_data argument
  SUBROUTINE get_restart_vars(var_data, patch_id, modelType, out_restartType)
    TYPE(t_p_var_list_element), ALLOCATABLE, INTENT(out) :: var_data(:)
    INTEGER, INTENT(in) :: patch_id
    CHARACTER(LEN = *), INTENT(IN) :: modelType
    INTEGER, OPTIONAL, INTENT(OUT) :: out_restartType
    INTEGER :: varCount, varIndex, error, curList, restartType
    TYPE(t_list_element), POINTER :: element
    CHARACTER(LEN = *), PARAMETER :: routine = modname//":get_restart_vars"

    ! init. main variables
    restartType = -1
    ! counts number of restart variables for this file (logical patch ident)
    varCount = number_of_restart_variables(patch_id, modelType)
    ! allocate the array of restart variables
    ALLOCATE(var_data(varCount), STAT = error)
    IF(error /= SUCCESS) CALL finish(routine, "memory allocation failed")
    ! fill the array of restart variables
    varIndex = 1
    DO curList = 1, nvar_lists
      IF(list_is_in_patch_restart(var_lists(curList), patch_id, modelType)) THEN
        IF(restartType == -1) THEN
          restartType = var_lists(curList)%p%restart_type
        ELSE IF(restartType /= var_lists(curList)%p%restart_type) THEN
          CALL finish(routine, "assertion failed: var_lists contains inconsistent restart_type values")
        END IF
        ! check, if the list has valid restart fields
        element => var_lists(curList)%p%first_list_element
        DO WHILE (ASSOCIATED(element))
          IF (element%field%info%lrestart) THEN
            var_data(varIndex)%p => element%field
            varIndex = varIndex + 1
          END IF
          element => element%next_list_element
        END DO
      END IF
    END DO
    IF(varIndex /= varCount + 1) &
         CALL finish(routine, "assertion failed: wrong restart variable count")
    IF (PRESENT(out_restartType)) out_restartType = restartType
  END SUBROUTINE get_restart_vars

  !> loop over all variable lists and call fn for each variable list field
  !!
  !! @param fn callback to evaluate for each field
  !! @param opt_list_select optional predicate function that needs to
  !! evaluate to false for variable lists whose variables are to be
  !! excluded from evaluation.
  !! @param state is passed to both fn and opt_list_select (if present)
  SUBROUTINE var_lists_apply(fn, state, opt_list_select)
    PROCEDURE(var_apply) :: fn
    CLASS(*), TARGET :: state
    PROCEDURE(list_select), OPTIONAL :: opt_list_select
    TYPE(t_list_element), POINTER :: element
    INTEGER :: i

    DO i = 1,nvar_lists
      IF (PRESENT(opt_list_select)) THEN
        IF (.NOT. opt_list_select(var_lists(i)%p, state)) CYCLE
      END IF
      element => var_lists(i)%p%first_list_element
      DO WHILE (ASSOCIATED(element))
        CALL fn(element%field, state, var_lists(i)%p)
        element => element%next_list_element
      ENDDO
    ENDDO
  END SUBROUTINE var_lists_apply

  !-------------------------------------------------------------------------
  !
  ! Get a list of variable names matching a given criterion.
  !
  SUBROUTINE get_all_var_names (varlist, ivar, opt_vlevel_type,           &
    &                           opt_vert_intp_type,                       &
    &                           opt_hor_intp_type, opt_lcontainer,        &
    &                           opt_loutput, opt_patch_id)
    CHARACTER(LEN=vname_len), INTENT(INOUT) :: varlist(:)
    INTEGER,                  INTENT(OUT)   :: ivar
    LOGICAL, OPTIONAL, INTENT(IN)  :: opt_vert_intp_type(SIZE(VINTP_TYPE_LIST))
    INTEGER, OPTIONAL, INTENT(IN)  :: opt_vlevel_type, opt_patch_id,   &
      &                               opt_hor_intp_type
    LOGICAL, OPTIONAL, INTENT(IN)  :: opt_lcontainer, opt_loutput
    ! local variables
    INTEGER                       :: i, hor_intp_type_match
    LOGICAL                       :: lcontainer, loutput_matters, loutput, &
         vert_intp_type(SIZE(VINTP_TYPE_LIST)), hor_intp_matters
    TYPE(t_list_element), POINTER :: element
    TYPE(t_var_metadata), POINTER :: info

    ! clear result
    DO i=1,SIZE(varlist)
      varlist(i) = ' '
    END DO

    ! default values for some criteria:
    lcontainer = .FALSE.
    CALL assign_if_present(lcontainer, opt_lcontainer)

    loutput_matters = PRESENT(opt_loutput)
    IF (loutput_matters) THEN
      loutput = opt_loutput
    ELSE
      loutput = .FALSE.
    END IF
    hor_intp_matters= PRESENT(opt_hor_intp_type)
    IF (hor_intp_matters) THEN
      hor_intp_type_match = opt_hor_intp_type
    ELSE
      hor_intp_type_match = -1
    END IF
    IF (PRESENT(opt_vert_intp_type)) THEN
      vert_intp_type = opt_vert_intp_type
    ELSE
      vert_intp_type = .FALSE.
    END IF
    !- loop over variables
    ivar = 0
    ! Note that there may be several variables with different time
    ! levels, we just add unconditionally all
    LOOP_VARLISTS : DO i = 1,nvar_lists
      IF (PRESENT(opt_patch_id)) THEN
        IF (var_lists(i)%p%patch_id /= opt_patch_id) CYCLE LOOP_VARLISTS
      END IF
      IF (PRESENT(opt_vlevel_type)) THEN
        IF(var_lists(i)%p%vlevel_type /= opt_vlevel_type) CYCLE LOOP_VARLISTS
      END IF
      IF (PRESENT(opt_loutput)) THEN
        ! Skip var_lists for which loutput .NEQV. opt_loutput
        IF (opt_loutput .NEQV. var_lists(i)%p%loutput) CYCLE LOOP_VARLISTS
      END IF
      element => var_lists(i)%p%first_list_element
      LOOPVAR : DO WHILE(ASSOCIATED(element))
        info => element%field%info
        ! Do not inspect element if it is a container
        IF ((info%lcontainer .EQV. lcontainer) &
             ! Do not inspect element if "loutput=.false."
             .AND. ((.NOT. loutput_matters) .OR. (loutput .EQV. info%loutput)) &
             ! Do not inspect element if it does not contain info for
             ! horizontal interpolation
             .AND. (.NOT. hor_intp_matters &
             &      .OR. info%hor_interp%hor_intp_type == hor_intp_type_match) &
             ) THEN
          ! Do not inspect element if it does not contain matching info for
          ! vertical interpolation
          IF (ALL(.NOT. vert_intp_type(:) .OR. &
               & info%vert_interp%vert_intp_type(:))) THEN
            ivar = ivar+1
            ! assign without time level suffix:
            varlist(ivar) = get_var_name(element%field)
          END IF
        END IF
        element => element%next_list_element
      ENDDO LOOPVAR ! loop over vlist "i"
    ENDDO LOOP_VARLISTS ! i = 1,nvar_lists

    CALL remove_duplicates(varlist, ivar)

  END SUBROUTINE get_all_var_names


  !------------------------------------------------------------------------------------------------
  !> @return Plain variable name (i.e. without TIMELEVEL_SUFFIX)
  !
  FUNCTION get_var_name(var)
    CHARACTER(LEN=VARNAME_LEN) :: get_var_name
    TYPE(t_var_list_element)   :: var
    ! local variable
    INTEGER :: idx

    idx = INDEX(var%info%name,TIMELEVEL_SUFFIX)
    IF (idx==0) THEN
      get_var_name = var%info%name
    ELSE
      get_var_name = var%info%name(1:idx-1)
    END IF
  END FUNCTION get_var_name

  !------------------------------------------------------------------------------------------------
  ! construct varname  with timelevel
  !
  FUNCTION get_varname_with_timelevel(varname,timelevel)
    CHARACTER(LEN=VARNAME_LEN) :: varname
    INTEGER, INTENT(IN)        :: timelevel

    CHARACTER(LEN=VARNAME_LEN) :: get_varname_with_timelevel

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
  FUNCTION get_var_timelevel(info)
    INTEGER :: get_var_timelevel
    TYPE(t_var_metadata), INTENT(IN) :: info
    ! local variable
    CHARACTER(LEN=*), PARAMETER :: routine = modname//':get_var_timelevel'
    INTEGER :: idx

    idx = INDEX(info%name,TIMELEVEL_SUFFIX)
    IF (idx == 0) THEN
      get_var_timelevel = -1
      RETURN
    END IF

    ! Get time level
    get_var_timelevel = ICHAR(info%name(idx+3:idx+3)) - ICHAR('0')
    IF(get_var_timelevel<=0 .OR. get_var_timelevel>MAX_TIME_LEVELS) &
      CALL finish(routine, 'Illegal time level in '//TRIM(info%name))
  END FUNCTION get_var_timelevel

  ! return logical if a variable name has a timelevel encoded
  LOGICAL FUNCTION has_time_level(varname)
    CHARACTER(LEN=*) :: varname
    CHARACTER(LEN=*), PARAMETER :: routine = modname//':has_time_level'
    INTEGER :: idx

    idx = INDEX(varname,TIMELEVEL_SUFFIX)
    has_time_level = (0 .EQ. idx)
  END FUNCTION

  !------------------------------------------------------------------------------------------------
  !> @return tile index (extracted from tile index suffix "t_") or "-1"
  !
  FUNCTION get_var_tileidx(varname)
    INTEGER :: get_var_tileidx
    CHARACTER(LEN=*) :: varname
    ! local variable
    CHARACTER(LEN=*), PARAMETER :: routine = modname//':get_var_tileidx'
    INTEGER :: idx

    idx = INDEX(varname,'_t_')
    IF (idx == 0) THEN
      get_var_tileidx = 0
      RETURN
    END IF

    ! Get time level
    get_var_tileidx = ICHAR(varname(idx+3:idx+3)) - ICHAR('0')
    IF (get_var_tileidx<=0) &
      CALL finish(routine, 'Illegal time level in '//TRIM(varname))
  END FUNCTION get_var_tileidx

  !------------------------------------------------------------------------------------------------
  !> @return container variable
  !
  FUNCTION get_var_container(patch_id, contained_elt)  RESULT(res)
    TYPE(t_var_list_element), POINTER :: res
    INTEGER, INTENT(IN) :: patch_id
    TYPE(t_var_list_element), TARGET, INTENT(IN) :: contained_elt

    ! local variables
    TYPE(t_list_element), POINTER :: element
    INTEGER                       :: i

    ! Unfortunately, there does not (yet) exist a link (pointer)
    ! between the contained element and the variable
    ! container. Therefore we need to loop over all variables and
    ! compare pointers.

    res => contained_elt
    VARLIST_LOOP : DO i = 1, nvar_lists
      IF (var_lists(i)%p%patch_id /= patch_id) CYCLE

      element => var_lists(i)%p%first_list_element
      DO WHILE (ASSOCIATED(element))
        IF (element%field%info%ncontained > 0) THEN
          IF (ASSOCIATED(element%field%r_ptr, contained_elt%r_ptr) .OR. &
            & ASSOCIATED(element%field%s_ptr, contained_elt%s_ptr) .OR. &
            & ASSOCIATED(element%field%i_ptr, contained_elt%i_ptr) .OR. &
            & ASSOCIATED(element%field%l_ptr, contained_elt%l_ptr)) THEN
            res => element%field
            EXIT VARLIST_LOOP
          END IF
        END IF

        element => element%next_list_element
      END DO
    END DO VARLIST_LOOP
  END FUNCTION get_var_container


  !------------------------------------------------------------------------------------------------
  !
  ! Change parameters of an already existent output var_list
  !
  SUBROUTINE set_var_list (this_list, output_type, restart_type,  &
       &                   post_suf, rest_suf, init_suf, loutput, &
       &                   lrestart, linitial, patch_id,          &
       &                   vlevel_type)
    !
    TYPE(t_var_list), INTENT(inout)        :: this_list      ! output var_list to change
    INTEGER,          INTENT(in), OPTIONAL :: output_type    ! 'GRIB' or 'NetCDF'
    INTEGER,          INTENT(in), OPTIONAL :: restart_type   ! 'GRIB' or 'NetCDF'
    CHARACTER(len=*), INTENT(in), OPTIONAL :: post_suf       ! suffix of output  file
    CHARACTER(len=*), INTENT(in), OPTIONAL :: rest_suf       ! suffix of restart file
    CHARACTER(len=*), INTENT(in), OPTIONAL :: init_suf       ! suffix of initial file
    LOGICAL,          INTENT(in), OPTIONAL :: loutput        ! in standard output file
    LOGICAL,          INTENT(in), OPTIONAL :: lrestart       ! in standard restartfile
    LOGICAL,          INTENT(in), OPTIONAL :: linitial       ! in standard initialfile
    INTEGER,          INTENT(in), OPTIONAL :: patch_id       ! patch ID
    INTEGER,          INTENT(in), OPTIONAL :: vlevel_type    ! 1/2/3 for model/pres./height levels
    !
    this_list%p%restart_type = restart_file_type

    CALL assign_if_present(this_list%p%output_type,  output_type)
    CALL assign_if_present(this_list%p%restart_type, restart_type)
    CALL assign_if_present(this_list%p%post_suf,     post_suf)
    CALL assign_if_present(this_list%p%rest_suf,     rest_suf)
    CALL assign_if_present(this_list%p%init_suf,     init_suf)
    CALL assign_if_present(this_list%p%loutput,      loutput)
    CALL assign_if_present(this_list%p%lrestart,     lrestart)
    CALL assign_if_present(this_list%p%linitial,     linitial)
    CALL assign_if_present(this_list%p%patch_id,     patch_id)
    CALL assign_if_present(this_list%p%vlevel_type,  vlevel_type)
    !
  END SUBROUTINE set_var_list
  !------------------------------------------------------------------------------------------------
  !
  ! Delete an output var_list, nullify the associated pointer
  !
  SUBROUTINE delete_var_list(this_list)
    !
    TYPE(t_var_list) :: this_list
    !
    IF (ASSOCIATED(this_list%p)) THEN
      CALL delete_list(this_list)
      DEALLOCATE(this_list%p)
    ENDIF
    !
  END SUBROUTINE delete_var_list
  !------------------------------------------------------------------------------------------------
  !
  ! Delete all output var_lists
  !
  SUBROUTINE delete_var_lists
    !
    TYPE(t_var_list), POINTER :: this_list
    !
    INTEGER :: i
    !
    DO i = 1, nvar_lists
      this_list => var_lists(i)
      CALL delete_var_list (this_list)
    END DO
    !
  END SUBROUTINE delete_var_lists
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
  !------------------------------------------------------------------------------------------------
  !
  ! Set default meta data of output var_list
  !
  SUBROUTINE default_var_list_settings (this_list,                                   &
       &                                filename,                                    &
       &                                loutput, lrestart, linitial,                 &
       &                                post_suf, rest_suf, init_suf,                &
       &                                output_type, restart_type, compression_type, &
       &                                model_type)
    !
    TYPE(t_var_list),   INTENT(inout)        :: this_list        ! output var_list
    LOGICAL,            INTENT(in), OPTIONAL :: loutput          ! to output
    LOGICAL,            INTENT(in), OPTIONAL :: lrestart         ! from/to restart
    LOGICAL,            INTENT(in), OPTIONAL :: linitial         ! from initial
    CHARACTER(len=*),   INTENT(in), OPTIONAL :: filename         ! name of output file
    CHARACTER(len=*),   INTENT(in), OPTIONAL :: post_suf         ! suffix of output  file
    CHARACTER(len=*),   INTENT(in), OPTIONAL :: rest_suf         ! suffix of restart file
    CHARACTER(len=*),   INTENT(in), OPTIONAL :: init_suf         ! suffix of initial file
    INTEGER,            INTENT(in), OPTIONAL :: output_type      ! output file type
    INTEGER,            INTENT(in), OPTIONAL :: restart_type     ! restart file type
    INTEGER,            INTENT(in), OPTIONAL :: compression_type ! compression type
    CHARACTER(len=*),   INTENT(in), OPTIONAL :: model_type       ! output file associated
    !
    this_list%p%restart_type = restart_file_type

    CALL assign_if_present (this_list%p%loutput,          loutput)
    CALL assign_if_present (this_list%p%lrestart,         lrestart)
    CALL assign_if_present (this_list%p%linitial,         linitial)
    CALL assign_if_present (this_list%p%filename,         filename)
    CALL assign_if_present (this_list%p%post_suf,         post_suf)
    CALL assign_if_present (this_list%p%rest_suf,         rest_suf)
    CALL assign_if_present (this_list%p%init_suf,         init_suf)
    CALL assign_if_present (this_list%p%output_type,      output_type)
    CALL assign_if_present (this_list%p%restart_type,     restart_type)
    CALL assign_if_present (this_list%p%compression_type, compression_type)
    CALL assign_if_present (this_list%p%model_type,       model_type)
    !
  END SUBROUTINE default_var_list_settings
  !------------------------------------------------------------------------------------------------
  SUBROUTINE default_var_list_metadata(this_info, this_list)
    !> memory info structure
    TYPE(t_var_metadata), INTENT(out) :: this_info
    !
    !> output var_list
    TYPE(t_var_list), INTENT(in)      :: this_list
    !

    this_info%key                 = 0
    this_info%name                = ''
    this_info%var_class           = CLASS_DEFAULT
    !
    this_info%cf                  = t_cf_var('', '', '', -1)
    this_info%grib2               = grib2_var(-1, -1, -1, -1, -1, -1)
    !
    this_info%allocated           = .FALSE.
    this_info%ndims               = 0
    this_info%used_dimensions(:)  = 0
    !
    ! RJ: Set default loutput to .TRUE., regardless of this_list%p%loutput
    this_info%loutput             = .TRUE.
    this_info%isteptype           = TSTEP_INSTANT
    this_info%resetval            = t_union_vals( 0.0_wp, 0.0_sp, 0, .FALSE.)
    this_info%lrestart            = this_list%p%lrestart
    this_info%lrestart_cont       = .FALSE.
    this_info%lrestart_read       = .FALSE.

    this_info%lmiss               = this_list%p%lmiss
    this_info%missval             = t_union_vals( 0.0_wp, 0.0_sp, 0, .FALSE.)
    this_info%lmask_boundary      = this_list%p%lmask_boundary

    this_info%initval             = t_union_vals( 0.0_wp, 0.0_sp, 0, .FALSE.)
    !
    this_info%lcontainer          = .FALSE.
    this_info%lcontained          = .FALSE.
    this_info%ncontained          = 0
    this_info%var_ref_pos         = -1 ! UNDEFINED
    this_info%maxcontained        = 0
    !
    this_info%hgrid               = -1
    this_info%vgrid               = -1
    !
    this_info%tlev_source         = TLEV_NNOW
    !
    this_info%cdiVarID            = CDI_UNDEFID
    this_info%cdiGridID           = CDI_UNDEFID
    !
    this_info%vert_interp         = create_vert_interp_metadata()
    this_info%hor_interp          = create_hor_interp_metadata()
    !
    this_info%post_op             = post_op()
    !
    this_info%in_group(:)         = groups()
    !
    this_info%action_list         = actions()
    !
    this_info%l_pp_scheduler_task = 0
    !
    this_info%lopenacc            = .FALSE.

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
    !
    lverbose = .FALSE.
    CALL assign_if_present (lverbose, verbose)
    !
    ! set components describing the 'Content of the field'
    CALL assign_if_present(info%data_type, data_type)
    !
    ! set components describing the 'Content of the field'
    !
    CALL assign_if_present (info%name,  name)
    CALL assign_if_present (info%var_class,  var_class)
    CALL struct_assign_if_present (info%cf,    cf)
    CALL struct_assign_if_present (info%grib2, grib2)
    !
    ! hash variable name for fast search
    !
    IF (PRESENT(name)) info%key = util_hashword(name//c_null_char, int(LEN_TRIM(name), C_SIZE_T), 0)
    !
!!$    CALL assign_if_present (info%used_dimensions(1:SIZE(ldims)), ldims)
    CALL assign_if_present (info%used_dimensions, ldims)
    !
    ! set grid type
    !
    CALL assign_if_present (info%hgrid, hgrid)
    CALL assign_if_present (info%vgrid, vgrid)
    !
    ! set flags concerning I/O
    !
    CALL assign_if_present (info%loutput,       loutput)
    CALL assign_if_present (info%lcontainer,    lcontainer)
    IF (info%lcontainer) THEN
      info%ncontained   =  0
      info%var_ref_pos  = -1 ! UNDEFINED
    END IF
    CALL struct_assign_if_present(info%resetval,resetval)
    CALL assign_if_present(info%isteptype, isteptype)
    CALL assign_if_present(info%lmiss, lmiss)
    CALL struct_assign_if_present(info%missval, missval)
    CALL assign_if_present(info%lrestart, lrestart)
    CALL assign_if_present(info%lrestart_cont, lrestart_cont)
    CALL struct_assign_if_present(info%initval, initval)
    CALL assign_if_present(info%tlev_source, tlev_source)
    !
    ! set flags concerning vertical interpolation
    CALL struct_assign_if_present (info%vert_interp,   vert_interp )

    ! set flags concerning horizontal interpolation
    CALL struct_assign_if_present (info%hor_interp,    hor_interp )

    ! set meta data containing the groups to which a variable belongs
    IF (PRESENT(in_group)) THEN
      info%in_group(1:SIZE(in_group)) = in_group(:)
    END IF

    CALL assign_if_present (info%l_pp_scheduler_task, l_pp_scheduler_task)

    CALL struct_assign_if_present (info%post_op, post_op)

    CALL struct_assign_if_present (info%action_list, action_list)

    ! indices of tracer in tracer container and in diagnostic container
    CALL assign_if_present (info%idx_tracer, idx_tracer)
    CALL assign_if_present (info%idx_diag, idx_diag)

    ! Create data on GPU
    CALL assign_if_present(info%lopenacc, lopenacc)

    ! perform consistency checks on variable's meta-data:
    CALL check_metadata_consistency(info)

    !
    ! printout (optional)
    !

    !LK    IF (lverbose) CALL print_var_metadata (info)
    !
  END SUBROUTINE set_var_metadata


  !------------------------------------------------------------------------------------------------
  !
  ! Set dynamic metadata, i.e. polymorphic tracer metadata
  ! (private routine within this module)
  !
  SUBROUTINE set_var_metadata_dyn(this_info_dyn,tracer_info)
    TYPE(t_var_metadata_dynamic),INTENT(INOUT) :: this_info_dyn
    CLASS(t_tracer_meta),INTENT(IN),OPTIONAL :: tracer_info

    CALL assign_if_present_tracer_meta(this_info_dyn%tracer,tracer_info)

  END SUBROUTINE set_var_metadata_dyn


  ! Auxiliary routine: initialize array, REAL(wp) variant
  SUBROUTINE init_array_r5d(ptr, linit, initval, lmiss, missval)
    REAL(wp),           POINTER     :: ptr(:,:,:,:,:)      ! pointer to field
    LOGICAL,            INTENT(IN)  :: linit, lmiss
    TYPE(t_union_vals), INTENT(IN)  :: initval, missval    ! optional initialization value
    REAL(wp) :: init_val

    IF (linit) THEN
      init_val = initval%rval
    ELSE IF (lmiss) THEN
      init_val = missval%rval
    ELSE
#if    defined (VARLIST_INITIZIALIZE_WITH_NAN) \
    && (defined (__INTEL_COMPILER) || defined (__PGI) || defined (NAGFOR))
      init_val = ieee_value(ptr, ieee_signaling_nan)
#else
      init_val = 0.0_wp
#endif
    END IF
!$omp parallel
    CALL init_contiguous_dp(ptr, SIZE(ptr), init_val)
!$omp end parallel
  END SUBROUTINE init_array_r5d

  ! Auxiliary routine: initialize array, REAL(sp) variant
  SUBROUTINE init_array_s5d(ptr, linit, initval, lmiss, missval)
    REAL(sp),           POINTER     :: ptr(:,:,:,:,:)      ! pointer to field
    LOGICAL,            INTENT(IN)  :: linit, lmiss
    TYPE(t_union_vals), INTENT(IN)  :: initval, missval    ! optional initialization value
    REAL(sp) :: init_val

    IF (linit) THEN
      init_val = initval%sval
    ELSE IF (lmiss) THEN
      init_val = missval%sval
    ELSE
#if    defined (VARLIST_INITIZIALIZE_WITH_NAN) \
    && (defined (__INTEL_COMPILER) || defined (__PGI) || defined (NAGFOR))
      init_val = ieee_value(ptr, ieee_signaling_nan)
#else
      init_val = 0.0_sp
#endif
    END IF
!$omp parallel
    CALL init_contiguous_sp(ptr, SIZE(ptr), init_val)
!$omp end parallel
  END SUBROUTINE init_array_s5d


  ! Auxiliary routine: initialize array, INTEGER variant
  SUBROUTINE init_array_i5d(ptr, linit, initval, lmiss, missval)
    INTEGER, POINTER                :: ptr(:,:,:,:,:)      ! pointer to field
    LOGICAL,            INTENT(IN)  :: linit, lmiss
    TYPE(t_union_vals), INTENT(IN)  :: initval, missval    ! optional initialization value
    INTEGER :: init_val

    IF (linit) THEN
      init_val = initval%ival
    ELSE IF (lmiss) THEN
      init_val = missval%ival
    ELSE
      init_val = 0
    END IF
!$omp parallel
    CALL init_contiguous_i4(ptr, SIZE(ptr), init_val)
!$omp end parallel
  END SUBROUTINE init_array_i5d


  ! Auxiliary routine: initialize array, REAL(wp) variant
  SUBROUTINE init_array_l5d(ptr, linit, initval, lmiss, missval)
    LOGICAL, POINTER                :: ptr(:,:,:,:,:)      ! pointer to field
    LOGICAL,            INTENT(IN)  :: linit, lmiss
    TYPE(t_union_vals), INTENT(IN)  :: initval, missval    ! optional initialization value
    LOGICAL :: init_val

    IF (linit) THEN
      init_val = initval%lval
    ELSE IF (lmiss) THEN
      init_val = missval%lval
    ELSE
      init_val = .FALSE.
    END IF
!$omp parallel
    CALL init_contiguous_l(ptr, SIZE(ptr), init_val)
!$omp end parallel
  END SUBROUTINE init_array_l5d


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
    REAL(wp),                TARGET &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
      , CONTIGUOUS &
#endif
      ,    OPTIONAL :: p5_r(:,:,:,:,:)              ! provided pointer
    REAL(sp),                TARGET &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
      , CONTIGUOUS &
#endif
      ,    OPTIONAL :: p5_s(:,:,:,:,:)              ! provided pointer
    INTEGER,                 TARGET &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
      , CONTIGUOUS &
#endif
      ,    OPTIONAL :: p5_i(:,:,:,:,:)              ! provided pointer
    LOGICAL,                 TARGET &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
      , CONTIGUOUS &
#endif
      ,    OPTIONAL :: p5_l(:,:,:,:,:)              ! provided pointer
    TYPE(t_vert_interp_meta),INTENT(in), OPTIONAL :: vert_interp                  ! vertical interpolation metadata
    TYPE(t_hor_interp_meta), INTENT(in), OPTIONAL :: hor_interp                   ! horizontal interpolation metadata
    LOGICAL,                 INTENT(in), OPTIONAL :: in_group(:)                  ! groups to which a variable belongs
    LOGICAL,                 INTENT(in), OPTIONAL :: verbose                      ! print information
    INTEGER,                 INTENT(in), OPTIONAL :: l_pp_scheduler_task          ! .TRUE., if field is updated by pp scheduler
    TYPE(t_post_op_meta),    INTENT(IN), OPTIONAL :: post_op                      ! "post-op" (small arithmetic operations) for this variable
    TYPE(t_var_action),      INTENT(IN), OPTIONAL :: action_list                  ! regularly triggered events
    CLASS(t_tracer_meta),    INTENT(in), OPTIONAL :: tracer_info                  ! tracer meta data
    REAL(wp),                INTENT(in), OPTIONAL :: initval_r                    ! value if var not available
    REAL(sp),                INTENT(in), OPTIONAL :: initval_s                    ! value if var not available
    INTEGER,                 INTENT(in), OPTIONAL :: initval_i                    ! value if var not available
    LOGICAL,                 INTENT(in), OPTIONAL :: initval_l                    ! value if var not available
    REAL(wp),                INTENT(in), OPTIONAL :: resetval_r                   ! reset value (after accumulation)
    REAL(sp),                INTENT(in), OPTIONAL :: resetval_s                   ! reset value (after accumulation)
    INTEGER,                 INTENT(in), OPTIONAL :: resetval_i                   ! reset value (after accumulation)
    LOGICAL,                 INTENT(in), OPTIONAL :: resetval_l                   ! reset value (after accumulation)
    REAL(wp),                INTENT(in), OPTIONAL :: missval_r                    ! missing value
    REAL(sp),                INTENT(in), OPTIONAL :: missval_s                    ! missing value
    INTEGER,                 INTENT(in), OPTIONAL :: missval_i                    ! missing value
    LOGICAL,                 INTENT(in), OPTIONAL :: missval_l                    ! missing value
    INTEGER,                 INTENT(in), OPTIONAL :: var_class                    !< variable type/species
    LOGICAL,                 INTENT(in), OPTIONAL :: lopenacc                     ! create variable on GPU

    ! local variables
    TYPE(t_union_vals) :: missval, initval, resetval
    INTEGER :: idims(5), istat
    LOGICAL :: referenced, is_restart_var
    CHARACTER(LEN = *), PARAMETER :: routine = modname//":add_var_list_element_5d"

    ! consistency check for restart and output
    TYPE(t_list_element), POINTER :: duplicate

    ! Check for a variable of the same name. This consistency check
    ! only makes sense for single-domain setups which, in addition,
    ! must not use internal post-processing (lon-lat or vertically
    ! interpolated output).
    IF (msg_level > 20) THEN
      duplicate => find_element(name)
      IF (ASSOCIATED(duplicate)) THEN
        CALL message('ADD_VAR:','Found double entry for varname:'//TRIM(name))
        NULLIFY(duplicate)
      ENDIF
    END IF

    is_restart_var = this_list%p%lrestart
    CALL assign_if_present(is_restart_var, lrestart)

    IF (PRESENT(lrestart)) THEN
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

    IF (PRESENT(p5_r) .OR. PRESENT(p5_s) .OR. PRESENT(p5_i) .OR. PRESENT(p5_l)) THEN
      referenced = .TRUE.
      new_list_element%field%info%allocated = .TRUE.
    ELSE
      referenced = .FALSE.
    ENDIF

    CALL assign_if_present(missval%rval, missval_r)
    CALL assign_if_present(missval%sval, missval_s)
    CALL assign_if_present(missval%ival, missval_i)
    CALL assign_if_present(missval%lval, missval_l)

    CALL assign_if_present(initval%rval, initval_r)
    CALL assign_if_present(initval%sval, initval_s)
    CALL assign_if_present(initval%ival, initval_i)
    CALL assign_if_present(initval%lval, initval_l)

    CALL assign_if_present(resetval%rval, resetval_r)
    CALL assign_if_present(resetval%sval, resetval_s)
    CALL assign_if_present(resetval%ival, resetval_i)
    CALL assign_if_present(resetval%lval, resetval_l)

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
    new_list_element%field%info%dom                      = this_list%p%patch_id
    IF (.NOT. referenced) THEN
      idims(1:ndims)    = new_list_element%field%info%used_dimensions(1:ndims)
      idims((ndims+1):) = 1

      NULLIFY(new_list_element%field%r_ptr)
      NULLIFY(new_list_element%field%s_ptr)
      NULLIFY(new_list_element%field%i_ptr)
      NULLIFY(new_list_element%field%l_ptr)
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

      IF (istat /= 0) THEN
        CALL finish(routine, 'allocation of array '//TRIM(name)//' failed')
      ELSE
        new_list_element%field%info%allocated = .TRUE.
      ENDIF
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

    IF(PRESENT(info)) info => new_list_element%field%info

    ! initialize the new array
    SELECT CASE(data_type)
    CASE (REAL_T)
      CALL init_array_r5d(new_list_element%field%r_ptr, linit=PRESENT(initval_r), initval=initval, &
        &                 lmiss=PRESENT(lmiss), missval=missval)
      !$ACC UPDATE DEVICE( new_list_element%field%r_ptr ) IF( new_list_element%field%info%lopenacc )
    CASE (SINGLE_T)
      CALL init_array_s5d(new_list_element%field%s_ptr, linit=PRESENT(initval_s), initval=initval, &
        &                 lmiss=PRESENT(lmiss), missval=missval)
      !$ACC UPDATE DEVICE( new_list_element%field%s_ptr ) IF( new_list_element%field%info%lopenacc )
    CASE (INT_T)
      CALL init_array_i5d(new_list_element%field%i_ptr, linit=PRESENT(initval_i), initval=initval, &
        &                 lmiss=PRESENT(lmiss), missval=missval)
      !$ACC UPDATE DEVICE( new_list_element%field%i_ptr ) IF( new_list_element%field%info%lopenacc )
    CASE (BOOL_T)
      CALL init_array_l5d(new_list_element%field%l_ptr, linit=PRESENT(initval_l), initval=initval, &
        &                 lmiss=PRESENT(lmiss), missval=missval)
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
    REAL(wp),                POINTER              :: ptr(:,:,:,:)                 ! reference to field
    INTEGER,                 INTENT(in)           :: hgrid                        ! horizontal grid type used
    INTEGER,                 INTENT(in)           :: vgrid                        ! vertical grid type used
    TYPE(t_cf_var),          INTENT(in)           :: cf                           ! CF related metadata
    TYPE(t_grib2_var),       INTENT(in)           :: grib2                        ! GRIB2 related metadata
    INTEGER,                 INTENT(in)           :: ldims(:)                     ! local dimensions
    LOGICAL,                 INTENT(in), OPTIONAL :: loutput                      ! output flag
    LOGICAL,                 INTENT(in), OPTIONAL :: lcontainer                   ! container flag
    LOGICAL,                 INTENT(in), OPTIONAL :: lrestart                     ! restart flag
    LOGICAL,                 INTENT(in), OPTIONAL :: lrestart_cont                ! continue restart if var not available
    REAL(wp),                INTENT(in), OPTIONAL :: initval                      ! value if var not available
    INTEGER,                 INTENT(in), OPTIONAL :: isteptype                    ! type of statistical processing
    REAL(wp),                INTENT(in), OPTIONAL :: resetval                     ! reset value (after accumulation)
    LOGICAL,                 INTENT(in), OPTIONAL :: lmiss                        ! missing value flag
    REAL(wp),                INTENT(in), OPTIONAL :: missval                      ! missing value
    INTEGER,                 INTENT(in), OPTIONAL :: tlev_source                  ! actual TL for TL dependent vars
    TYPE(t_var_metadata),    POINTER,    OPTIONAL :: info                         ! returns reference to metadata
    REAL(wp),                TARGET &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
      , CONTIGUOUS &
#endif
      ,    OPTIONAL :: p5(:,:,:,:,:)                ! provided pointer
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
    REAL(wp),                POINTER              :: ptr(:,:,:)                   ! reference to field
    INTEGER,                 INTENT(in)           :: hgrid                        ! horizontal grid type used
    INTEGER,                 INTENT(in)           :: vgrid                        ! vertical grid type used
    TYPE(t_cf_var),          INTENT(in)           :: cf                           ! CF related metadata
    TYPE(t_grib2_var),       INTENT(in)           :: grib2                        ! GRIB2 related metadata
    INTEGER,                 INTENT(in)           :: ldims(:)                     ! local dimensions
    LOGICAL,                 INTENT(in), OPTIONAL :: loutput                      ! output flag
    LOGICAL,                 INTENT(in), OPTIONAL :: lcontainer                   ! container flag
    LOGICAL,                 INTENT(in), OPTIONAL :: lrestart                     ! restart flag
    LOGICAL,                 INTENT(in), OPTIONAL :: lrestart_cont                ! continue restart if var not available
    REAL(wp),                INTENT(in), OPTIONAL :: initval                      ! value if var not available
    INTEGER,                 INTENT(in), OPTIONAL :: isteptype                    ! type of statistical processing
    REAL(wp),                INTENT(in), OPTIONAL :: resetval                     ! reset value (after accumulation)
    LOGICAL,                 INTENT(in), OPTIONAL :: lmiss                        ! missing value flag
    REAL(wp),                INTENT(in), OPTIONAL :: missval                      ! missing value
    INTEGER,                 INTENT(in), OPTIONAL :: tlev_source                  ! actual TL for TL dependent vars
    CLASS(t_tracer_meta),    INTENT(in), OPTIONAL :: tracer_info                  ! tracer meta data
    TYPE(t_var_metadata),    POINTER,    OPTIONAL :: info                         ! returns reference to metadata
    REAL(wp),                TARGET &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
      , CONTIGUOUS &
#endif
      ,    OPTIONAL :: p5(:,:,:,:,:)                ! provided pointer
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
    REAL(wp),                POINTER              :: ptr(:,:)                     ! reference to field
    INTEGER,                 INTENT(in)           :: hgrid                        ! horizontal grid type used
    INTEGER,                 INTENT(in)           :: vgrid                        ! vertical grid type used
    TYPE(t_cf_var),          INTENT(in)           :: cf                           ! CF related metadata
    TYPE(t_grib2_var),       INTENT(in)           :: grib2                        ! GRIB2 related metadata
    INTEGER,                 INTENT(in)           :: ldims(:)                     ! local dimensions
    LOGICAL,                 INTENT(in), OPTIONAL :: loutput                      ! output flag
    LOGICAL,                 INTENT(in), OPTIONAL :: lcontainer                   ! container flag
    LOGICAL,                 INTENT(in), OPTIONAL :: lrestart                     ! restart flag
    LOGICAL,                 INTENT(in), OPTIONAL :: lrestart_cont                ! continue restart if var not available
    REAL(wp),                INTENT(in), OPTIONAL :: initval                      ! value if var not available
    INTEGER,                 INTENT(in), OPTIONAL :: isteptype                    ! type of statistical processing
    REAL(wp),                INTENT(in), OPTIONAL :: resetval                     ! reset value (after accumulation)
    LOGICAL,                 INTENT(in), OPTIONAL :: lmiss                        ! missing value flag
    REAL(wp),                INTENT(in), OPTIONAL :: missval                      ! missing value
    INTEGER,                 INTENT(in), OPTIONAL :: tlev_source                  ! actual TL for TL dependent vars
    CLASS(t_tracer_meta),    INTENT(in), OPTIONAL :: tracer_info                  ! tracer meta data
    TYPE(t_var_metadata),    POINTER,    OPTIONAL :: info                         ! returns reference to metadata
    REAL(wp),                TARGET &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
      , CONTIGUOUS &
#endif
      ,    OPTIONAL :: p5(:,:,:,:,:)                ! provided pointer
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
    REAL(wp),                POINTER              :: ptr(:)                       ! reference to field
    INTEGER,                 INTENT(in)           :: hgrid                        ! horizontal grid type used
    INTEGER,                 INTENT(in)           :: vgrid                        ! vertical grid type used
    TYPE(t_cf_var),          INTENT(in)           :: cf                           ! CF related metadata
    TYPE(t_grib2_var),       INTENT(in)           :: grib2                        ! GRIB2 related metadata
    INTEGER,                 INTENT(in)           :: ldims(:)                     ! local dimensions
    LOGICAL,                 INTENT(in), OPTIONAL :: loutput                      ! output flag
    LOGICAL,                 INTENT(in), OPTIONAL :: lcontainer                   ! container flag
    LOGICAL,                 INTENT(in), OPTIONAL :: lrestart                     ! restart flag
    LOGICAL,                 INTENT(in), OPTIONAL :: lrestart_cont                ! continue restart if var not available
    REAL(wp),                INTENT(in), OPTIONAL :: initval                      ! value if var not available
    INTEGER,                 INTENT(in), OPTIONAL :: isteptype                    ! type of statistical processing
    REAL(wp),                INTENT(in), OPTIONAL :: resetval                     ! reset value (after accumulation)
    LOGICAL,                 INTENT(in), OPTIONAL :: lmiss                        ! missing value flag
    REAL(wp),                INTENT(in), OPTIONAL :: missval                      ! missing value
    INTEGER,                 INTENT(in), OPTIONAL :: tlev_source                  ! actual TL for TL dependent vars
    TYPE(t_var_metadata),    POINTER,    OPTIONAL :: info                         ! returns reference to metadata
    REAL(wp),                TARGET &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
      , CONTIGUOUS &
#endif
      ,    OPTIONAL :: p5(:,:,:,:,:)                ! provided pointer
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
    REAL(sp),                POINTER              :: ptr(:,:,:,:)                 ! reference to field
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
    REAL(sp),                TARGET &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
      , CONTIGUOUS &
#endif
      ,    OPTIONAL :: p5(:,:,:,:,:)                ! provided pointer
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
    REAL(sp),                POINTER              :: ptr(:,:,:)                   ! reference to field
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
    REAL(sp),                TARGET &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
      , CONTIGUOUS &
#endif
      ,    OPTIONAL :: p5(:,:,:,:,:)                ! provided pointer
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
    REAL(sp),                POINTER              :: ptr(:,:)                     ! reference to field
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
    REAL(sp),                TARGET &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
      , CONTIGUOUS &
#endif
      ,    OPTIONAL :: p5(:,:,:,:,:)                ! provided pointer
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
    REAL(sp),                POINTER              :: ptr(:)                       ! reference to field
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
    REAL(sp),                TARGET &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
      , CONTIGUOUS &
#endif
      ,    OPTIONAL :: p5(:,:,:,:,:)                ! provided pointer
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
    INTEGER,                 POINTER              :: ptr(:,:,:,:)                 ! reference to field
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
    INTEGER,                 TARGET &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
      , CONTIGUOUS &
#endif
      ,    OPTIONAL :: p5(:,:,:,:,:)                ! provided pointer
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
    INTEGER,                 POINTER              :: ptr(:,:,:)                   ! reference to field
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
    INTEGER,                 TARGET &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
      , CONTIGUOUS &
#endif
      ,    OPTIONAL :: p5(:,:,:,:,:)                ! provided pointer
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
    INTEGER,                 POINTER              :: ptr(:,:)                     ! reference to field
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
    INTEGER,                 TARGET &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
      , CONTIGUOUS &
#endif
      ,    OPTIONAL :: p5(:,:,:,:,:)                ! provided pointer
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
    INTEGER,                 POINTER              :: ptr(:)                       ! reference to field
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
    INTEGER,                 TARGET &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
      , CONTIGUOUS &
#endif
      ,    OPTIONAL :: p5(:,:,:,:,:)                ! provided pointer
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
    LOGICAL,                 POINTER              :: ptr(:,:,:,:)                 ! reference to field
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
    LOGICAL,                 TARGET &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
      , CONTIGUOUS &
#endif
      ,    OPTIONAL :: p5(:,:,:,:,:)                ! provided pointer
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
    LOGICAL,                 POINTER              :: ptr(:,:,:)                   ! reference to field
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
    LOGICAL,                 TARGET &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
      , CONTIGUOUS &
#endif
      ,    OPTIONAL :: p5(:,:,:,:,:)                ! provided pointer
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
    LOGICAL,                 POINTER              :: ptr(:,:)                     ! reference to field
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
    LOGICAL,                 TARGET &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
      , CONTIGUOUS &
#endif
      ,    OPTIONAL :: p5(:,:,:,:,:)                ! provided pointer
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
    LOGICAL,                 POINTER              :: ptr(:)                       ! reference to field
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
    LOGICAL,                 TARGET &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
      , CONTIGUOUS &
#endif
      ,    OPTIONAL :: p5(:,:,:,:,:)                ! provided pointer
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
    REAL(wp),         POINTER    :: ptr(:,:,:,:,:) ! reference to allocated field
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
    REAL(wp),         POINTER    :: ptr(:,:,:,:) ! reference to allocated field
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
    REAL(wp),         POINTER    :: ptr(:,:,:)   ! reference to allocated field
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
    REAL(wp),         POINTER    :: ptr(:,:)   ! reference to allocated field
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
    REAL(wp),         POINTER    :: ptr(:)    ! reference to allocated field
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
    REAL(wp),         POINTER    :: ptr(:)      ! reference to allocated field
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
    REAL(wp),         POINTER    :: ptr(:,:)    ! reference to allocated field
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
    REAL(wp),         POINTER    :: ptr(:,:,:)   ! reference to allocated field
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
    REAL(wp), POINTER                                :: ptr(:,:,:)
    INTEGER,                 INTENT(in)              :: hgrid                      ! horizontal grid type used
    INTEGER,                 INTENT(in)              :: vgrid                      ! vertical grid type used
    TYPE(t_cf_var),          INTENT(in)              :: cf                         ! CF related metadata
    TYPE(t_grib2_var),       INTENT(in)              :: grib2                      ! GRIB2 related metadata
    INTEGER,                 INTENT(in)              :: ref_idx                    ! idx of slice to be referenced
    INTEGER,                 INTENT(in)              :: ldims(3)                   ! local dimensions, for checking
    LOGICAL,                 INTENT(in),    OPTIONAL :: loutput                    ! output flag
    LOGICAL,                 INTENT(in),    OPTIONAL :: lrestart                   ! restart flag
    LOGICAL,                 INTENT(in),    OPTIONAL :: lrestart_cont              ! continue restart if var not available
    REAL(wp),                INTENT(in),    OPTIONAL :: initval                    ! value if var not available
    INTEGER,                 INTENT(in),    OPTIONAL :: isteptype                  ! type of statistical processing
    REAL(wp),                INTENT(in),    OPTIONAL :: resetval                   ! reset value (after accumulation)
    LOGICAL,                 INTENT(in),    OPTIONAL :: lmiss                      ! missing value flag
    REAL(wp),                INTENT(in),    OPTIONAL :: missval                    ! missing value
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
    missvalt = ref_info%missval
    initvalt = ref_info%initval
    resetvalt= ref_info%resetval
    !
    CALL assign_if_present(missvalt%rval,  missval)
    CALL assign_if_present(initvalt%rval,  initval)
    CALL assign_if_present(resetvalt%rval, resetval)
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
    !
    IF(PRESENT(info)) info => new_list_element%field%info
    !
    IF (PRESENT(lmiss)) THEN
      ptr = new_list_element%field%info%missval%rval
    END IF
    !
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
    REAL(wp), POINTER                             :: ptr(:,:)
    INTEGER,                 INTENT(in)           :: hgrid                       ! horizontal grid type used
    INTEGER,                 INTENT(in)           :: vgrid                       ! vertical grid type used
    TYPE(t_cf_var),          INTENT(in)           :: cf                          ! CF related metadata
    TYPE(t_grib2_var),       INTENT(in)           :: grib2                       ! GRIB2 related metadata
    INTEGER,                 INTENT(in)           :: ref_idx                     ! idx of slice to be referenced
    INTEGER,                 INTENT(in)           :: ldims(2)                    ! local dimensions, for checking
    LOGICAL,                 INTENT(in), OPTIONAL :: loutput                     ! output flag
    LOGICAL,                 INTENT(in), OPTIONAL :: lrestart                    ! restart flag
    LOGICAL,                 INTENT(in), OPTIONAL :: lrestart_cont               ! continue restart if var not available
    REAL(wp),                INTENT(in), OPTIONAL :: initval                     ! value if var not available
    INTEGER,                 INTENT(in), OPTIONAL :: isteptype                   ! type of statistical processing
    REAL(wp),                INTENT(in), OPTIONAL :: resetval                    ! reset value (after accumulation)
    LOGICAL,                 INTENT(in), OPTIONAL :: lmiss                       ! missing value flag
    REAL(wp),                INTENT(in), OPTIONAL :: missval                     ! missing value
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
    missvalt = ref_info%missval
    initvalt = ref_info%initval
    resetvalt= ref_info%resetval
    !
    CALL assign_if_present(missvalt%rval,  missval)
    CALL assign_if_present(initvalt%rval,  initval)
    CALL assign_if_present(resetvalt%rval, resetval)
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
    !
    IF (PRESENT(lmiss)) THEN
      ptr = new_list_element%field%info%missval%rval
    END IF
    !
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
    missvalt = ref_info%missval
    initvalt = ref_info%initval
    resetvalt= ref_info%resetval
    !
    CALL assign_if_present(missvalt%sval,  missval)
    CALL assign_if_present(initvalt%sval,  initval)
    CALL assign_if_present(resetvalt%sval, resetval)
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
    !
    IF (PRESENT(lmiss)) THEN
      ptr = new_list_element%field%info%missval%sval
    END IF
    !
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
    missvalt = ref_info%missval
    initvalt = ref_info%initval
    resetvalt= ref_info%resetval
    !
    CALL assign_if_present(missvalt%sval,  missval)
    CALL assign_if_present(initvalt%sval,  initval)
    CALL assign_if_present(resetvalt%sval, resetval)
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
    !
    IF (PRESENT(lmiss)) THEN
      ptr = new_list_element%field%info%missval%sval
    END IF
    !
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
    missvalt = ref_info%missval
    initvalt = ref_info%initval
    resetvalt= ref_info%resetval
    !
    CALL assign_if_present(missvalt%ival,  missval)
    CALL assign_if_present(initvalt%ival,  initval)
    CALL assign_if_present(resetvalt%ival, resetval)
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
    !
    IF (PRESENT(lmiss)) THEN
      ptr = new_list_element%field%info%missval%ival
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



  !================================================================================================
  !------------------------------------------------------------------------------------------------
  !
  ! add supplementary fields to a different var list (eg. geopotential, surface pressure, ...)
  !
  SUBROUTINE add_var_list_reference_named(to_var_list, name, &
       from_var_list_name, loutput, bit_precision, in_group)
    TYPE(t_var_list), TARGET, INTENT(inout)  :: to_var_list
    CHARACTER(len=*), INTENT(in)             :: name
    CHARACTER(len=*), INTENT(in)             :: from_var_list_name
    LOGICAL,          INTENT(in),   OPTIONAL :: loutput
    INTEGER,          INTENT(in),   OPTIONAL :: bit_precision
    LOGICAL,          INTENT(in),   OPTIONAL :: in_group(MAX_GROUPS)  ! groups to which a variable belongs
    !
    TYPE(t_var_list), POINTER :: from_var_list
    !
    CALL get_var_list(from_var_list, from_var_list_name)
    IF (ASSOCIATED(from_var_list)) THEN
      CALL add_var_list_reference(to_var_list, name, from_var_list, &
        &                         loutput, bit_precision, in_group)
    END IF
    !
  END SUBROUTINE add_var_list_reference_named

  !------------------------------------------------------------------------------------------------
  !
  ! add supplementary fields to a different var list (eg. geopotential, surface pressure, ...)
  !
  SUBROUTINE add_var_list_reference(to_var_list, name, from_var_list, loutput, bit_precision, in_group)
    TYPE(t_var_list), INTENT(inout)          :: to_var_list
    CHARACTER(len=*), INTENT(in)             :: name
    TYPE(t_var_list), TARGET, INTENT(in)     :: from_var_list
    LOGICAL,          INTENT(in),   OPTIONAL :: loutput
    INTEGER,          INTENT(in),   OPTIONAL :: bit_precision
    LOGICAL,          INTENT(in),   OPTIONAL :: in_group(MAX_GROUPS)  ! groups to which a variable belongs
    !
    TYPE(t_var_list_element), POINTER :: source
    TYPE(t_list_element),     POINTER :: element
    !
    element => find_list_element(from_var_list, name)
    IF (ASSOCIATED(element)) THEN
      source => element%field
      CALL append_list_element(to_var_list, element)
      element%field                = source
      element%field%info%allocated = .FALSE.
      element%field%info%lrestart  = .FALSE.
      CALL assign_if_present(element%field%info%loutput, loutput)
      CALL assign_if_present(element%field%info%grib2%bits, bit_precision)
      if (present(in_group)) then
        element%field%info%in_group(:)=in_group(:)
      end if
    ENDIF
    !
  END SUBROUTINE add_var_list_reference

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
    CALL print_var_list_intrinsic(this_list%p, lshort)
    !
  END SUBROUTINE print_var_list
  !------------------------------------------------------------------------------------------------
  !
  ! print current memory table
  !
  SUBROUTINE print_var_list_intrinsic(this_list, lshort)
    TYPE(t_var_list_intrinsic),  INTENT(in) :: this_list ! list
    LOGICAL, OPTIONAL :: lshort
    !
    TYPE(t_list_element), POINTER :: this_list_element

    CALL message('','')
    CALL message('','')
    CALL message('','Status of variable list '//TRIM(this_list%name)//':')
    CALL message('','')
    !
    this_list_element => this_list%first_list_element
    DO WHILE (ASSOCIATED(this_list_element))
      CALL print_var(this_list_element%field, lshort)
      this_list_element => this_list_element%next_list_element
    ENDDO
  END SUBROUTINE print_var_list_intrinsic

  !---------------------------------------------------------------------
  !
  ! print single variable
  !
  SUBROUTINE print_var(field, lshort)
    TYPE(t_var_list_element),  INTENT(in) :: field ! list
    LOGICAL, OPTIONAL :: lshort
    !
    LOGICAL :: short, lfirst
    CHARACTER(len=32) :: dimension_text, dtext
    INTEGER :: i, igrp, ivintp_type, tlen, alen
    CHARACTER(len=4) :: localMode

    localMode = '----'
    short = .FALSE.
    CALL assign_if_present(short,lshort)
    IF (short) THEN
      IF (field%info%name /= '' .AND. .NOT. field%info%lcontainer) THEN
        IF (field%info%lrestart) localMode(1:1) = 'r'
        IF (field%info%lcontained) localMode(2:2) = 't'
        SELECT CASE (field%info%isteptype)
        CASE (1)
          localMode(3:3) = 'i'
        CASE (2)
          localMode(3:3) = 'm'
        CASE (3)
          localMode(3:3) = 'a'
        END SELECT
        SELECT CASE (field%info%hgrid)
        CASE (1)
          localMode(4:4) = 'c'
        CASE (2)
          localMode(4:4) = 'v'
        CASE (3)
          localMode(4:4) = 'e'
        END SELECT

        WRITE(message_text, '(a4,3i4,a24,a40)') localMode,                                 &
           &                              field%info%grib2%discipline, &
           &                              field%info%grib2%category,   &
           &                              field%info%grib2%number,     &
           &                              TRIM(field%info%name),       &
           &                              TRIM(field%info%cf%standard_name)
        CALL message('', message_text)

      ENDIF

    ELSE

      IF (field%info%name /= '' .AND. .NOT. field%info%lcontainer) THEN
        !
        message_text = 'Table entry name                            : ' &
          // field%info%name
        CALL message('', message_text)
        WRITE (message_text,'(a,i32.1)') &
          'Key entry                                   : ', field%info%key

        CALL message('', message_text)
        !
        IF (ASSOCIATED(field%r_ptr) .OR. &
          & ASSOCIATED(field%s_ptr) .OR. &
          & ASSOCIATED(field%i_ptr) .OR. &
          & ASSOCIATED(field%l_ptr)) THEN
          CALL message ('','Pointer status                              : in use.')
          dimension_text = '('
          tlen = 1
          DO i = 1, field%info%ndims
            WRITE(dtext,'(i0)') field%info%used_dimensions(i)
            alen = LEN_TRIM(dtext)
            dimension_text(tlen+1:tlen+alen) = dtext(1:alen)
            tlen = tlen + alen + 1
            dimension_text(tlen:tlen) = MERGE(',', ')', field%info%ndims /= i)
          ENDDO
          message_text = 'Local field dimensions                      : ' &
               // dimension_text(1:tlen)
          CALL message('', message_text)
        ELSE
          CALL message('', 'Pointer status                              : not in use.')
        ENDIF
        !
        WRITE (message_text,'(a,3i4)') &
             'Assigned GRIB discipline/category/parameter : ', &
             field%info%grib2%discipline,    &
             field%info%grib2%category,      &
             field%info%grib2%number
        CALL message('', message_text)
        !
        WRITE (message_text,'(a,a,a,a)')                          &
             'CF convention standard name/unit            : ',    &
             TRIM(field%info%cf%standard_name), &
             '     ',                                             &
             TRIM(field%info%cf%units)
        CALL message('', message_text)
        !
        message_text = 'CF convention long name                     : ' &
             // field%info%cf%long_name
        !
        IF (field%info%lcontained) THEN
          CALL message('', 'Field is in a container                     : yes.')
          WRITE (message_text,'(a,i2)')                        &
             ' Index in container                          : ',&
             field%info%ncontained
          CALL message('', message_text)
        ELSE
          CALL message('', 'Field is in a container                     : no.')
          CALL message('', ' Index in container                          : --')
        ENDIF
        !
        WRITE (message_text,'(a,i2)')                          &
             ' horizontal grid type used (C=1,V=2,E=3)     : ',&
             field%info%hgrid
        CALL message('', message_text)
        !
        WRITE (message_text,'(a,i2)')                          &
             ' vertical grid type used (see cdilib.c)      : ',&
             field%info%vgrid
        CALL message('', message_text)
        !
        WRITE (message_text,'(a,i2)')                          &
             ' type of stat. processing (I=1,AVG=2,ACC=3...: ',&
             field%info%isteptype
        CALL message('', message_text)
        !
        IF (field%info%lmiss) THEN
          IF (ASSOCIATED(field%r_ptr)) THEN
            WRITE (message_text,'(a,e20.12)')      &
                 'Missing value                               : ', &
                 field%info%missval%rval
          ELSE IF (ASSOCIATED(field%s_ptr)) THEN
            WRITE (message_text,'(a,e20.12)')      &
                 'Missing value                               : ', &
                 field%info%missval%sval
          ELSE IF (ASSOCIATED(field%i_ptr)) THEN
            WRITE (message_text,'(a,i8)')      &
                 'Missing value                               : ', &
                 field%info%missval%ival
          ELSE IF (ASSOCIATED(field%l_ptr)) THEN
            WRITE (message_text,'(a,l8)')      &
                 'Missing value                               : ', &
                 field%info%missval%lval
          ENDIF
          CALL message('', message_text)
        ELSE
          CALL message('', 'Missing values                              : off.')
        ENDIF
        !
        IF (field%info%lrestart) THEN
          CALL message('', 'Added to restart                            : yes.')
        ELSE
          CALL message('', 'Added to Restart                            : no.')
        ENDIF
        !
        IF (field%info_dyn%tracer%lis_tracer) THEN
          CALL message('', 'Tracer field                                : yes.')

          IF (field%info_dyn%tracer%lfeedback) THEN
            CALL message('', 'Child-to-parent feedback                  : yes.')
          ELSE
            CALL message('', 'Child-to-parent feedback                  : no.')
          ENDIF

          WRITE (message_text,'(a,3i3)') &
             'Horizontal transport method                 : ', &
             field%info_dyn%tracer%ihadv_tracer
          CALL message('', message_text)

          WRITE (message_text,'(a,3i3)') &
             'Vertical transport method                   : ', &
             field%info_dyn%tracer%ivadv_tracer
          CALL message('', message_text)

          IF (field%info_dyn%tracer%lturb_tracer) THEN
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
             field%info%var_class
        CALL message('', message_text)

        !
        ! print groups, to which this variable belongs:
        IF (ANY(field%info%in_group(:))) THEN
          message_text = 'Variable group(s)                           :'
          tlen = LEN_TRIM(message_text)
          lfirst = .TRUE.
          DO igrp=1,SIZE(field%info%in_group)
            IF (field%info%in_group(igrp)) THEN
              message_text(tlen+1:tlen+2) = MERGE("  ", ", ", lfirst)
              tlen = tlen + 1 + MERGE(1, 0, .NOT. lfirst)
              alen = LEN_TRIM(var_groups_dyn%name(igrp))
              message_text(tlen+1:tlen+alen) = var_groups_dyn%name(igrp)(1:alen)
              tlen = tlen + alen
            ENDIF
          END DO
          CALL message('', message_text)
        END IF

        !
        ! print horizontal and vertical interpolation method(s):
        message_text = 'Horizontal interpolation                    : ' &
             // STR_HINTP_TYPE(field%info%hor_interp%hor_intp_type)
        CALL message('', message_text)

        LOOP_VINTP_TYPES : DO ivintp_type=1,SIZE(VINTP_TYPE_LIST)
          IF (field%info%vert_interp%vert_intp_type(ivintp_type)) THEN
            WRITE (message_text,'(a)')  &
              &  'Vertical interpolation                      : '//  &
              &  toupper(VINTP_TYPE_LIST(ivintp_type))
            CALL message('', message_text)
          END IF
        END DO LOOP_VINTP_TYPES
        CALL message('', '')
      ENDIF

    ENDIF
  END SUBROUTINE print_var


  !------------------------------------------------------------------------------------------------
  !
  ! print all var lists
  !
  SUBROUTINE print_all_var_lists(lshort)
    LOGICAL, OPTIONAL, INTENT(in) :: lshort
    INTEGER :: i
    DO i=1,nvar_lists
      CALL print_var_list(var_lists(i), lshort)
    END DO
  END SUBROUTINE print_all_var_lists
  !------------------------------------------------------------------------------------------------
  !
  ! print current stat table
  !
  SUBROUTINE print_sinfo (this_list)
    TYPE(t_var_list),  INTENT(in) :: this_list
    !
    WRITE (message_text,'(a16,a)') TRIM(this_list%p%name), '-buffer: '
    CALL message('',message_text)
    CALL message('','')
    CALL message('','')
    CALL message('','Statistic of base memory:')
    CALL message('','')
    !
    !LK    CALL print_sinfo_list (this_list)
    !
  END SUBROUTINE print_sinfo



  !> Loops over all variables and collects the variables names
  !  corresponding to the group @p grp_name
  !
  SUBROUTINE collect_group(grp_name, var_name, nvars,       &
    &                      loutputvars_only, lremap_lonlat, &
    &                      opt_vlevel_type, opt_dom_id,     &
    &                      opt_lquiet)
    CHARACTER(LEN=*),           INTENT(IN)    :: grp_name
    CHARACTER(LEN=VARNAME_LEN), INTENT(OUT)   :: var_name(:)
    INTEGER,                    INTENT(OUT)   :: nvars
    ! loutputvars_only: If set to .TRUE. all variables in the group
    ! which have the the loutput flag equal to .FALSE. are skipped.
    LOGICAL,                    INTENT(IN)    :: loutputvars_only
    ! lremap_lonlat: If set to .TRUE. only variables in the group
    ! which can be interpolated onto lon-lat grids are considered.
    LOGICAL,                    INTENT(IN)    :: lremap_lonlat

    ! 1: model levels, 2: pressure levels, 3: height level
    INTEGER, OPTIONAL,          INTENT(IN)    :: opt_vlevel_type
    ! (optional:) domain id
    INTEGER, OPTIONAL,          INTENT(IN)    :: opt_dom_id
    ! (optional:) quiet mode (no log output)
    LOGICAL, OPTIONAL,          INTENT(IN)    :: opt_lquiet

    ! local variables
    CHARACTER(*), PARAMETER :: routine = modname//":collect_group"
    INTEGER                       :: i, grp_id, llmsg_len
    TYPE(t_list_element), POINTER :: element
    TYPE(t_var_metadata), POINTER :: info
    CHARACTER(LEN=VARNAME_LEN)    :: name
    CHARACTER(len=*), PARAMETER   :: llmsg = " lon-lat"
    LOGICAL                       :: lquiet, verbose, skip

    nvars  = 0
    grp_id = var_groups_dyn%group_id(grp_name)
    lquiet = .FALSE.
    IF (PRESENT(opt_lquiet))  lquiet = opt_lquiet
    verbose = .NOT. lquiet

    ! loop over all variable lists and variables
    DO i = 1,nvar_lists
      IF (PRESENT(opt_vlevel_type)) THEN
        IF (var_lists(i)%p%vlevel_type /= opt_vlevel_type) CYCLE
      ENDIF
      IF (PRESENT(opt_dom_id)) THEN
        ! do not inspect variable list if its domain does not match:
        IF (var_lists(i)%p%patch_id /= opt_dom_id)  CYCLE
      END IF

      element => var_lists(i)%p%first_list_element
      LOOPVAR : DO WHILE (ASSOCIATED(element))
        info => element%field%info
        ! Do not inspect element if it is a container
        IF (.NOT. info%lcontainer .AND. info%in_group(grp_id)) THEN
          name = get_var_name(element%field)
          llmsg_len = 0
          ! Skip element if we need only output variables:
          skip = loutputvars_only .AND. &
            & ((.NOT. info%loutput) .OR. (.NOT. var_lists(i)%p%loutput))

          IF (lremap_lonlat) THEN
            skip = skip .OR. info%hgrid /= GRID_UNSTRUCTURED_CELL
            llmsg_len = LEN(llmsg)
          ELSE
            ! If no lon-lat interpolation is requested for this output file,
            ! skip all variables of this kind:
            skip = skip .OR. &
                 (loutputvars_only .AND. (info%hgrid == GRID_REGULAR_LONLAT))
          END IF

          IF (.NOT. skip) THEN
            nvars = nvars + 1
            var_name(nvars) = name
          ELSE IF (verbose) THEN
            CALL message(routine, "Skipping variable "//TRIM(name)//" for " &
                 //llmsg(1:llmsg_len)//"output.")
          END IF
        END IF
        element => element%next_list_element
      ENDDO LOOPVAR ! loop over vlist "i"
    ENDDO ! i = 1,nvar_lists

    CALL remove_duplicates(var_name, nvars)

  END SUBROUTINE collect_group


  !------------------------------------------------------------------------------------------------
  SUBROUTINE assign_if_present_cf (y,x)
    TYPE(t_cf_var), INTENT(inout)        :: y
    TYPE(t_cf_var), INTENT(in) ,OPTIONAL :: x
    IF (PRESENT(x)) y = x
  END SUBROUTINE assign_if_present_cf
  !------------------------------------------------------------------------------------------------
  SUBROUTINE assign_if_present_grib2 (y,x)
    TYPE(t_grib2_var), INTENT(inout)        :: y
    TYPE(t_grib2_var) ,INTENT(in) ,OPTIONAL :: x
    IF (PRESENT(x)) y = x
  END SUBROUTINE assign_if_present_grib2
  !------------------------------------------------------------------------------------------------
  SUBROUTINE assign_if_present_union (y,x)
    TYPE(t_union_vals), INTENT(inout)        :: y
    TYPE(t_union_vals) ,INTENT(in) ,OPTIONAL :: x
    IF (PRESENT(x)) y = x
  END SUBROUTINE assign_if_present_union
  !------------------------------------------------------------------------------------------------
  SUBROUTINE assign_if_present_tracer_meta (y,x)
    CLASS(t_tracer_meta), POINTER, INTENT(out) :: y
    CLASS(t_tracer_meta) ,INTENT(in) ,OPTIONAL :: x
    IF (PRESENT(x)) THEN
      ALLOCATE(y, source=x)
    ELSE
      ALLOCATE(t_tracer_meta :: y)
      SELECT TYPE(y)
        TYPE IS(t_tracer_meta)
          y = create_tracer_metadata(lis_tracer=.FALSE.)
      END SELECT
    ENDIF
  END SUBROUTINE assign_if_present_tracer_meta
  !------------------------------------------------------------------------------------------------
  SUBROUTINE assign_if_present_vert_interp (y,x)
    TYPE(t_vert_interp_meta), INTENT(inout)        :: y
    TYPE(t_vert_interp_meta) ,INTENT(in) ,OPTIONAL :: x
    IF (PRESENT(x)) y = x
  END SUBROUTINE assign_if_present_vert_interp
  !------------------------------------------------------------------------------------------------
  SUBROUTINE assign_if_present_hor_interp (y,x)
    TYPE(t_hor_interp_meta), INTENT(inout)        :: y
    TYPE(t_hor_interp_meta) ,INTENT(in) ,OPTIONAL :: x
    IF (PRESENT(x)) y = x
  END SUBROUTINE assign_if_present_hor_interp
  !------------------------------------------------------------------------------------------------
  SUBROUTINE assign_if_present_post_op (y,x)
    TYPE(t_post_op_meta), INTENT(inout)        :: y
    TYPE(t_post_op_meta) ,INTENT(in) ,OPTIONAL :: x
    IF (PRESENT(x)) y = x
  END SUBROUTINE assign_if_present_post_op
  !------------------------------------------------------------------------------------------------
  SUBROUTINE assign_if_present_action_list (y,x)
    TYPE(t_var_action), INTENT(inout)        :: y
    TYPE(t_var_action) ,INTENT(in) ,OPTIONAL :: x
    IF (PRESENT(x)) y = x
  END SUBROUTINE assign_if_present_action_list
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
    !
    TYPE(t_var_list),   INTENT(in) :: this_list
    CHARACTER(len=*),   INTENT(in) :: name
    INTEGER, OPTIONAL              :: opt_hgrid
    LOGICAL, OPTIONAL              :: opt_caseInsensitive
    !
    TYPE(t_list_element), POINTER :: element
    INTEGER :: key,hgrid
    LOGICAL :: name_has_time_level
    LOGICAL :: case_insensitive
    case_insensitive = .FALSE.
    CALL assign_if_present(case_insensitive, opt_caseInsensitive)

    hgrid = -1
    CALL assign_if_present(hgrid,opt_hgrid)
    !
    key = util_hashword(name, INT(LEN_TRIM(name), C_SIZE_T), 0)
    name_has_time_level = has_time_level(name)
    !
    element => this_list%p%first_list_element
    DO WHILE (ASSOCIATED(element))
      IF (-1 == hgrid .OR. hgrid == element%field%info%hgrid) THEN
        IF (elementFoundByName(key,name,name_has_time_level,&
          &                    element,case_insensitive)) RETURN
      ENDIF
      element => element%next_list_element
    ENDDO
    !
  END FUNCTION find_list_element
  
  !------------------------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  !
  ! Overloaded to search for a tracer by its index (ncontained)
  !
  FUNCTION find_tracer_by_index (this_list, ncontained, opt_hgrid) RESULT(this_list_element)
    !
    TYPE(t_var_list),   INTENT(in) :: this_list
    INTEGER,            INTENT(in) :: ncontained
    INTEGER, OPTIONAL              :: opt_hgrid
    !
    TYPE(t_list_element), POINTER  :: this_list_element
    INTEGER :: key,hgrid

    hgrid = -1
    CALL assign_if_present(hgrid,opt_hgrid)
    !
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
    !
    NULLIFY (this_list_element)
    !
  END FUNCTION find_tracer_by_index

  !-----------------------------------------------------------------------------
  !
  ! Find named list element across all known variable lists
  !
  FUNCTION find_element_from_all (name, opt_patch_id, opt_hgrid, opt_caseInsensitive,opt_returnList) RESULT(element)
    CHARACTER(len=*),   INTENT(in) :: name
    INTEGER, OPTIONAL              :: opt_patch_id
    INTEGER, OPTIONAL              :: opt_hgrid
    LOGICAL, OPTIONAL              :: opt_caseInsensitive
    TYPE(t_var_list), POINTER, OPTIONAL     :: opt_returnList

    TYPE(t_list_element), POINTER :: element
    INTEGER :: i,patch_id

    patch_id = 1
    CALL assign_if_present(patch_id,opt_patch_id)

    DO i=1,nvar_lists
      IF ( patch_id /= var_lists(i)%p%patch_id ) CYCLE
      element => find_list_element(var_lists(i),name,opt_hgrid,opt_caseInsensitive)
      IF (ASSOCIATED (element)) THEN
#ifdef DEBUG_MVSTREAM
        if (my_process_is_stdio()) call &
            & print_summary('destination PATCHID:'//TRIM(int2string(patch_id)))
#endif
        IF (PRESENT(opt_returnList)) opt_returnList => var_lists(i)
        RETURN
      END IF
    END DO
  END FUNCTION! find_element_from_all_lists

  !-----------------------------------------------------------------------------
  !
  ! Find named list element across selected known variable lists
  !
  FUNCTION find_field_from_selected(name, list_is_selected, opt_state, &
       opt_caseInsensitive) RESULT(field)
    CHARACTER(len=*),   INTENT(in) :: name
    PROCEDURE(list_select) :: list_is_selected
    CLASS(*), TARGET, OPTIONAL :: opt_state
    LOGICAL, OPTIONAL :: opt_caseInsensitive

    TYPE(t_list_element), POINTER :: element
    TYPE(t_var_list_element), POINTER :: field
    INTEGER :: i

    NULLIFY(field)
    DO i=1,nvar_lists
      IF (list_is_selected(var_lists(i)%p, opt_state)) THEN
        element => find_list_element(var_lists(i), name, &
             opt_caseInsensitive=opt_caseInsensitive)
        IF (ASSOCIATED (element)) THEN
          field => element%field
          RETURN
        END IF
      END IF
    END DO
  END FUNCTION find_field_from_selected

  !-----------------------------------------------------------------------------
  !
  ! Find named list element across selected known variable lists
  !
  FUNCTION find_field_from_selected2(var_is_selected, list_is_selected, &
       state) RESULT(field)
    PROCEDURE(var_filter) :: var_is_selected
    PROCEDURE(list_select) :: list_is_selected
    CLASS(*), TARGET :: state

    TYPE(t_list_element), POINTER :: element
    TYPE(t_var_list_element), POINTER :: field
    INTEGER :: i

    NULLIFY(field)
    DO i=1,nvar_lists
      IF (list_is_selected(var_lists(i)%p, state)) THEN
        element => var_lists(i)%p%first_list_element
        DO WHILE (ASSOCIATED(element))
          IF (var_is_selected(element%field, state, var_lists(i)%p)) THEN
            field => element%field
            RETURN
          END IF
          element => element%next_list_element
        ENDDO
      END IF
    END DO
  END FUNCTION find_field_from_selected2

  !-----------------------------------------------------------------------------
  !
  ! (Un)pack the var_lists
  ! This IS needed for the restart modules that need to communicate the var_lists from the worker PEs to dedicated restart PEs.
  !
  SUBROUTINE varlistPacker(operation, packedMessage)
    INTEGER, VALUE :: operation
    TYPE(t_PackedMessage), INTENT(INOUT) :: packedMessage

    INTEGER :: info_size, iv, nv, nelems, patch_id, restart_type, vlevel_type, n, ierrstat
    INTEGER, ALLOCATABLE            :: info_storage(:)
    TYPE(t_list_element), POINTER   :: element, newElement
    TYPE(t_var_metadata)            :: info
    TYPE(t_var_list)                :: p_var_list
    CHARACTER(LEN=max_var_list_name_len) :: var_list_name
    CHARACTER(LEN=32)               :: model_type
    LOGICAL                         :: lrestart

    CHARACTER(LEN=*), PARAMETER :: routine = modname//':varlistPacker'

    ! delete old var lists
    IF(operation == kUnpackOp) CALL delete_var_lists

    ! get the size - in default INTEGER words - which is needed to
    ! hold the contents of TYPE(t_var_metadata)
    info_size = SIZE(TRANSFER(info, (/ 0 /)))
    ALLOCATE(info_storage(info_size), STAT=ierrstat)
    IF(ierrstat /= SUCCESS) CALL finish(routine, "memory allocation failure")

    ! get the number of var lists
    nv = nvar_lists
    CALL packedMessage%packer(operation, nv)

    ! for each var list, get its components
    DO iv = 1, nv
        IF(operation == kPackOp) THEN
            ! copy the values needed for the new_var_list() CALL to local variables
            lrestart = var_lists(iv)%p%lrestart
            var_list_name = var_lists(iv)%p%name
            model_type = var_lists(iv)%p%model_type
            patch_id = var_lists(iv)%p%patch_id
            restart_type = var_lists(iv)%p%restart_type
            vlevel_type = var_lists(iv)%p%vlevel_type

            ! count the number of variable restart entries
            element => var_lists(iv)%p%first_list_element
            nelems = 0
            DO WHILE (ASSOCIATED(element))
                IF(element%field%info%lrestart) nelems = nelems+1
                element => element%next_list_element
            END DO
        END IF
        CALL packedMessage%packer(operation, lrestart)
        CALL packedMessage%packer(operation, var_list_name)
        CALL packedMessage%packer(operation, model_type)
        CALL packedMessage%packer(operation, patch_id)
        CALL packedMessage%packer(operation, restart_type)
        CALL packedMessage%packer(operation, vlevel_type)
        CALL packedMessage%packer(operation, nelems)

        IF(.NOT. lrestart) CYCLE  ! transfer only a restart var_list
        IF(nelems == 0) CYCLE ! check if there are valid restart fields

        IF(operation == kPackOp) THEN
            element => var_lists(iv)%p%first_list_element
            DO WHILE (ASSOCIATED(element))
                IF(element%field%info%lrestart) THEN
                    info_storage = TRANSFER(element%field%info, (/ 0 /))
                    CALL packedMessage%packer(operation, info_storage)
                END IF
                element => element%next_list_element
            END DO
        END IF

        IF(operation == kUnpackOp) THEN
            ! create var list
            CALL new_var_list(p_var_list, var_list_name, patch_id=patch_id, restart_type=restart_type, vlevel_type=vlevel_type, &
                             &lrestart=.TRUE.)
            p_var_list%p%model_type = TRIM(model_type)

            ! insert elements into var list
            DO n = 1, nelems
                ! ALLOCATE a new element
                ALLOCATE(newElement, STAT=ierrstat)
                IF(ierrstat /= SUCCESS) CALL finish(routine, "memory allocation failure")
                IF(n == 1) THEN   ! the first element pointer needs to be stored IN a different variable than the later pointers (there are no double pointers IN FORTRAN...)
                    p_var_list%p%first_list_element => newElement
                ELSE
                    element%next_list_element => newElement
                END IF
                element => newElement
                element%next_list_element => NULL()

                ! these pointers don't make sense on the restart PEs, NULLIFY them
                NULLIFY(element%field%r_ptr, element%field%s_ptr, element%field%i_ptr, element%field%l_ptr)
                element%field%var_base_size = 0 ! Unknown here

                ! set info structure from binary representation in info_storage
                CALL packedMessage%packer(operation, info_storage)
                element%field%info = TRANSFER(info_storage, info)
            END DO
        END IF

    END DO
  END SUBROUTINE varlistPacker


  !>  Detailed print-out of variable groups.
  !
  SUBROUTINE print_group_details(idom, opt_latex_fmt, opt_reduce_trailing_num, opt_skip_trivial)
    INTEGER, INTENT(IN)           :: idom          !< domain ID
    LOGICAL, INTENT(IN), OPTIONAL :: opt_latex_fmt !< Flag: .TRUE., if output shall be formatted for LaTeX
    LOGICAL, INTENT(IN), OPTIONAL :: opt_reduce_trailing_num !< Flag: replace trailing numbers by "*"
    LOGICAL, INTENT(IN), OPTIONAL :: opt_skip_trivial        !< Flag: skip empty of single-entry groups
    ! local variables
    CHARACTER(*), PARAMETER :: routine = modname//"::print_group_details"
    CHARACTER(len=VARNAME_LEN), ALLOCATABLE :: group_names(:)
    CHARACTER(LEN=VARNAME_LEN), ALLOCATABLE :: grp_vars(:), grp_vars_output(:)
    INTEGER                                 :: ngrp_vars, ngrp_vars_output, ierrstat, i, j
    LOGICAL                                 :: latex_fmt, reduce_trailing_num, skip_trivial

    latex_fmt = .FALSE.
    IF (PRESENT(opt_latex_fmt))  latex_fmt = opt_latex_fmt 
    reduce_trailing_num = .FALSE.
    IF (PRESENT(opt_reduce_trailing_num))  reduce_trailing_num = opt_reduce_trailing_num 
    skip_trivial = .FALSE.
    IF (PRESENT(opt_skip_trivial))  skip_trivial = opt_skip_trivial 

    IF (latex_fmt) THEN
      WRITE (0,*) " "
      WRITE (0,*) "% ---------------------------------------"
      WRITE (0,'(a,i0,a)') " % Variable group info (for domain #", idom, "):"
      WRITE (0,*) "% ---------------------------------------"
      WRITE (0,*) "% "
      WRITE (0,*) "% LaTeX formatted output, requires suitable environment 'varlist' and"
      WRITE (0,*) "% macros 'varname' and 'grpname'."
      WRITE (0,*) " "
    ELSE
      WRITE (0,*) " "
      WRITE (0,*) "---------------------------------------"
      WRITE (0,'(a,i0,a)') " Variable group info (for domain #", idom, "):"
      WRITE (0,*) "---------------------------------------"
      WRITE (0,*) " "
    END IF

    group_names = var_groups_dyn%alphabetical_list()
    IF (latex_fmt) THEN
      WRITE (0,*) "% List of groups:"
      CALL pretty_print_string_list(group_names, opt_prefix=" %    ")
    ELSE
      WRITE (0,*) "List of groups:"
      CALL pretty_print_string_list(group_names, opt_prefix="    ")
    END IF
    
    ! temporary variables needed for variable group parsing
    i = total_number_of_variables()
    ALLOCATE(grp_vars_output(i), grp_vars(i), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
    
    DO i=1,SIZE(group_names)
      
      ! for each group we collect the contained variables two times
      ! - the first time we collect *all* model level variables on
      ! the triangular grid, and the second time we collect only
      ! those variables which are available for output.
      CALL collect_group(group_names(i), grp_vars, ngrp_vars,    &
        &               loutputvars_only = .FALSE.,              &
        &               lremap_lonlat    = .FALSE.,              &
        &               opt_vlevel_type  = level_type_ml,        &
        &               opt_dom_id       = idom,                 &
        &               opt_lquiet       = .TRUE.)
     
      IF (ngrp_vars > 0) THEN
        
        CALL lowcase(grp_vars(1:ngrp_vars))
        ! (optionally) replace trailing numbers by "*"
        IF (reduce_trailing_num) &
             CALL star_out_trailing_num(ngrp_vars, grp_vars)
        CALL quicksort(grp_vars(1:ngrp_vars))

        IF ((skip_trivial) .AND. (ngrp_vars <= 1))  CYCLE
               
        CALL collect_group(group_names(i), grp_vars_output, ngrp_vars_output,    &
          &               loutputvars_only = .TRUE.,               &
          &               lremap_lonlat    = .FALSE.,              &
          &               opt_vlevel_type  = level_type_ml,        &
          &               opt_dom_id       = idom,                 &
          &               opt_lquiet       = .TRUE.)

        CALL lowcase(grp_vars_output(1:ngrp_vars_output))
        ! (optionally) replace trailing numbers by "*"
        IF (reduce_trailing_num) &
          CALL star_out_trailing_num(ngrp_vars_output, grp_vars_output)
        CALL quicksort(grp_vars_output(1:ngrp_vars_output))

        IF (latex_fmt) THEN
          DO j=1,ngrp_vars
            grp_vars(j) = "\varname{"//TRIM(grp_vars(j))//"}"
          END DO
          DO j=1,ngrp_vars_output
            grp_vars_output(j) = "\varname{"//TRIM(grp_vars_output(j))//"}"
          END DO
        END IF
        
        WRITE (0,*) ' '
        IF (latex_fmt) THEN
          WRITE (0,*) "\begin{varlist}{\grpname{"//TRIM(group_names(i))//"}}"
        ELSE
          WRITE (0,*) 'GROUP "', TRIM(group_names(i)), '":'
        END IF
        CALL pretty_print_string_list(grp_vars_output(1:ngrp_vars_output), opt_prefix="    ")
        
        CALL difference(grp_vars, ngrp_vars, grp_vars_output, ngrp_vars_output)
        IF (ngrp_vars > 0) THEN
          WRITE (0,*) " "
          WRITE (0,*) "   Non-output variables:"
          WRITE (0,*) " "

          CALL pretty_print_string_list(grp_vars(1:ngrp_vars), opt_prefix="    ")
        END IF
        
      ELSE

        IF ((skip_trivial) .AND. (ngrp_vars <= 1))  CYCLE

        WRITE (0,*) ' '
        IF (latex_fmt) THEN
          WRITE (0,*) "\begin{varlist}{\grpname{"//TRIM(group_names(i))//"}}"
        ELSE
          WRITE (0,*) 'GROUP "', TRIM(group_names(i)), '":'
        END IF
        WRITE (0,*) "   -- empty --"
        
      END IF

      IF (latex_fmt) THEN
        WRITE (0,*) "\end{varlist}"
      END IF
    END DO
    WRITE (0,*) " "
    
    DEALLOCATE(group_names, grp_vars, grp_vars_output, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')

  END SUBROUTINE print_group_details

  SUBROUTINE star_out_trailing_num(ngrp_vars, grp_vars)
    INTEGER, INTENT(inout) :: ngrp_vars
    CHARACTER(LEN=VARNAME_LEN), INTENT(inout) :: grp_vars(ngrp_vars)
    INTEGER :: slen(ngrp_vars), j, k, t
    LOGICAL :: lfound
    ngrp_vars = SIZE(grp_vars)
    DO j=1,ngrp_vars
      slen(j) = find_trailing_number(grp_vars(j))
    END DO
    search_loop: DO j=1,ngrp_vars
      t = slen(j)
      IF (t /= -1) THEN
        ! replace trailing number if any other variable of "this
        ! kind" exists:
        lfound = .FALSE.
        DO k=j+1,ngrp_vars
          IF (t == slen(k)) THEN
            IF (grp_vars(j)(1:t-1) == grp_vars(k)(1:t-1)) THEN
              grp_vars(k)(t:) = "*"
              slen(k) = -1
              lfound = .TRUE.
            END IF
          END IF
        END DO
        IF (lfound) grp_vars(j)(t:) = "*"
      END IF
    END DO search_loop
    CALL remove_duplicates(grp_vars, ngrp_vars)
  END SUBROUTINE star_out_trailing_num

#ifndef NOMPI
  !> mirrors all variable meta-data in var_lists on another group of processes
  !! is_sender must be consistently set for the sending group
  !! (typically workers)
  SUBROUTINE replicate_var_lists(intercomm, bcast_root, is_sender)
    INTEGER, INTENT(in) :: intercomm, bcast_root
    LOGICAL, INTENT(in) :: is_sender

    ! var_list_name should have at least the length of var_list names
    ! (although this doesn't matter as long as it is big enough for every name)
    CHARACTER(LEN=max_var_list_name_len) :: var_list_name
    INTEGER :: info_size, iv, nv,nelems, n, ierror, list_info(4)
    INTEGER, ALLOCATABLE          :: info_storage(:,:)
    TYPE(t_list_element), POINTER :: element
    TYPE(t_var_list)              :: p_var_list
    TYPE(t_var_metadata)          :: info
    CHARACTER(len=*), PARAMETER :: routine = modname//'replicate_var_lists'
    !---------------------------------------------------------------------
    ! Replicate variable lists

    ! Get the size - in default INTEGER words - which is needed to
    ! hold the contents of TYPE(t_var_metadata)
    info_size = SIZE(TRANSFER(info, (/ 0 /)))

    ! Get the number of var_lists
    IF (is_sender) THEN
      nv = nvar_lists
      list_info(1) = nv
      list_info(2) = 0
      DO iv = 1, nv
        element => var_lists(iv)%p%first_list_element
        nelems = 0
        DO WHILE (ASSOCIATED(element))
          nelems = nelems+1
          element => element%next_list_element
        ENDDO
        IF (nelems > list_info(2)) list_info(2) = nelems
      END DO
    END IF
    CALL p_bcast(list_info(1:2), bcast_root, intercomm)
    nv = list_info(1)

    ALLOCATE(info_storage(info_size, list_info(2)), STAT=ierror)
    IF (ierror /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")

    ! For each var list, get its components
    DO iv = 1, nv

      ! Send name
      IF (is_sender) var_list_name = var_lists(iv)%p%name
      CALL p_bcast(var_list_name, bcast_root, intercomm)

      IF (is_sender) THEN

        ! Count the number of variable entries
        element => var_lists(iv)%p%first_list_element
        nelems = 0
        DO WHILE (ASSOCIATED(element))
          nelems = nelems+1
          element => element%next_list_element
        ENDDO

        ! Gather the components needed for name list I/O and send them.
        ! Please note that not the complete list is replicated, unneeded
        ! entries are left away!

        list_info(1) = nelems
        list_info(2) = var_lists(iv)%p%patch_id
        list_info(3) = var_lists(iv)%p%vlevel_type
        list_info(4) = MERGE(1,0,var_lists(iv)%p%loutput)

      ENDIF

      ! Send basic info:

      CALL p_bcast(list_info, bcast_root, intercomm)

      IF (.NOT. is_sender) THEN
        nelems = list_info(1)
        ! Create var list
        CALL new_var_list(p_var_list, var_list_name, patch_id=list_info(2), &
                          vlevel_type=list_info(3), loutput=(list_info(4)/=0) )
      ENDIF

      ! Get the binary representation of all info members of the variables
      ! of the list and send it to the receiver.
      ! Using the Fortran TRANSFER intrinsic may seem like a hack,
      ! but it has the advantage that it is completely independet of the
      ! actual declaration if TYPE(t_var_metadata).
      ! Thus members may added to or removed from TYPE(t_var_metadata)
      ! without affecting the code below and we don't have an additional
      ! cross dependency between TYPE(t_var_metadata) and this module.

      IF (is_sender) THEN
        element => var_lists(iv)%p%first_list_element
        nelems = 0
        DO WHILE (ASSOCIATED(element))
          nelems = nelems+1
          info_storage(:,nelems) = TRANSFER(element%field%info, (/ 0 /))
          element => element%next_list_element
        ENDDO
      ENDIF

      ! Send binary representation of all info members

      CALL p_bcast(info_storage, bcast_root, intercomm)

      IF (.NOT. is_sender) THEN
        IF (nelems > 0) THEN
          ! Insert elements into var list
          ALLOCATE(var_lists(iv)%p%first_list_element)
          element => var_lists(iv)%p%first_list_element

          DO n = 1, nelems-1
            ! Nullify all pointers in element%field,
            ! reconstruct these later if needed
            NULLIFY(element%field%r_ptr, element%field%s_ptr, &
              &     element%field%i_ptr, element%field%l_ptr)
            element%field%var_base_size = 0 ! Unknown here

            ! Set info structure from binary representation in info_storage
            element%field%info = TRANSFER(info_storage(:, n), info)
            ALLOCATE(element%next_list_element)
            element => element%next_list_element
          ENDDO
          element%field%info = TRANSFER(info_storage(:, nelems), info)
          NULLIFY(element%next_list_element, &
            &     element%field%r_ptr, element%field%s_ptr, &
            &     element%field%i_ptr, element%field%l_ptr)
          element%field%var_base_size = 0 ! Unknown here
        ELSE
          NULLIFY(var_lists(iv)%p%first_list_element)
        END IF
      ENDIF

    ENDDO

    DEALLOCATE(info_storage, STAT=ierror)
    IF (ierror /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")

  END SUBROUTINE replicate_var_lists
#endif

END MODULE mo_var_list
