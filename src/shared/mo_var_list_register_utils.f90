! @par Copyright and Li/cense
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
MODULE mo_var_list_register_utils
  USE mo_var_groups,       ONLY: var_groups_dyn, MAX_GROUPS
  USE mo_var_metadata_types, ONLY: t_var_metadata
  USE mo_var_metadata,     ONLY: get_var_name
  USE mo_var, ONLY: level_type_ml, t_var
  USE mo_var_list,         ONLY: find_list_element, t_var_list_ptr
  USE mo_exception,        ONLY: finish
  USE mo_util_string,      ONLY: remove_duplicates, pretty_print_string_list, &
    &                            lowcase, difference, find_trailing_number
  USE mo_impl_constants,   ONLY: vname_len, SUCCESS
  USE mo_cdi_constants,    ONLY: GRID_UNSTRUCTURED_CELL, GRID_REGULAR_LONLAT
  USE mo_util_sort,        ONLY: quicksort
  USE mo_var_list_register, ONLY: vlr_get, t_vl_register_iter

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: vlr_add_vref, vlr_find, vlr_print_vls
  PUBLIC :: vlr_print_groups, vlr_group

  CHARACTER(*), PARAMETER :: modname = "mo_var_list_register_utils"

CONTAINS
  ! add supplementary fields to a different var list (eg. geopotential, surface pressure, ...)
  SUBROUTINE vlr_add_vref(to_vl, vname, from_vl, loutput, bit_precision, in_group)
    TYPE(t_var_list_ptr), INTENT(INOUT) :: to_vl
    CHARACTER(*), INTENT(IN) :: vname, from_vl
    LOGICAL, INTENT(in), OPTIONAL :: loutput, in_group(MAX_GROUPS)
    INTEGER, INTENT(in), OPTIONAL :: bit_precision
    TYPE(t_var_list_ptr) :: vlp_from
    TYPE(t_var), POINTER :: n_e, o_e => NULL()

    CALL vlr_get(vlp_from, from_vl)
    IF (ASSOCIATED(vlp_from%p)) o_e => find_list_element(vlp_from, vname)
    IF (ASSOCIATED(o_e)) THEN
      ALLOCATE(n_e, SOURCE=o_e)
      n_e%info%allocated = .FALSE.
      n_e%info%lrestart  = .FALSE.
      IF (PRESENT(loutput))       n_e%info%loutput     = loutput
      IF (PRESENT(bit_precision)) n_e%info%grib2%bits  = bit_precision
      IF (PRESENT(in_group))      n_e%info%in_group(:) = in_group(:)
      CALL to_vl%register(n_e)
    ENDIF
  END SUBROUTINE vlr_add_vref

  ! print all var lists
  SUBROUTINE vlr_print_vls(lshort)
    LOGICAL, OPTIONAL, INTENT(IN) :: lshort
    TYPE(t_vl_register_iter) :: iter

    DO WHILE(iter%next())
      CALL iter%cur%print(lshort=lshort)
    END DO
  END SUBROUTINE vlr_print_vls

  !> Loops over all variables and collects the variables names
  !  corresponding to the group @p grp_name
  SUBROUTINE vlr_group(grp_name, var_name, nvars,       &
    &                      loutputvars_only, lremap_lonlat, &
    &                      opt_vlevel_type, opt_dom_id)
    CHARACTER(*), INTENT(IN) :: grp_name
    CHARACTER(LEN=vname_len), INTENT(OUT) :: var_name(:)
    INTEGER, INTENT(OUT) :: nvars
    LOGICAL, INTENT(IN) :: loutputvars_only, lremap_lonlat
    INTEGER, OPTIONAL, INTENT(IN) :: opt_vlevel_type, opt_dom_id
    TYPE(t_vl_register_iter) :: iter
    CHARACTER(*), PARAMETER :: routine = modname//":collect_group"
    INTEGER :: grp_id, i
    TYPE(t_var_metadata), POINTER :: info

    nvars  = 0
    grp_id = var_groups_dyn%group_id(grp_name)
    DO WHILE(iter%next())
      IF (PRESENT(opt_vlevel_type)) THEN
        IF (iter%cur%p%vlevel_type /= opt_vlevel_type) CYCLE
      ENDIF
      IF (PRESENT(opt_dom_id)) THEN
        IF (iter%cur%p%patch_id /= opt_dom_id)  CYCLE
      END IF
      DO i = 1, iter%cur%p%nvars
        info => iter%cur%p%vl(i)%p%info
        IF (info%lcontainer .OR. .NOT.info%in_group(grp_id)) CYCLE
        IF (loutputvars_only .AND. (.NOT.info%loutput .OR. &
          & .NOT.iter%cur%p%loutput)) CYCLE
        IF (lremap_lonlat) THEN
          IF (info%hgrid .NE. GRID_UNSTRUCTURED_CELL) CYCLE
        ELSE
          IF (loutputvars_only .AND. (info%hgrid == GRID_REGULAR_LONLAT)) &
            & CYCLE
        END IF
        nvars = nvars + 1
        var_name(nvars) = get_var_name(info)
      ENDDO ! loop over vlist "i"
    ENDDO ! i = 1, SIZE(var_lists)
    CALL remove_duplicates(var_name, nvars)
  END SUBROUTINE vlr_group

  ! Find named list element accross all knows variable lists
  FUNCTION vlr_find(vname, opt_patch_id, opt_hgrid, opt_list, opt_cs, opt_output) RESULT(element)
    CHARACTER(*), INTENT(in) :: vname
    INTEGER, OPTIONAL, INTENT(IN) :: opt_patch_id, opt_hgrid
    TYPE(t_var_list_ptr), OPTIONAL, INTENT(OUT) :: opt_list
    LOGICAL, OPTIONAL, INTENT(IN) :: opt_cs, opt_output
    TYPE(t_var), POINTER :: element
    INTEGER :: patch_id
    LOGICAL :: output, check_output
    TYPE(t_vl_register_iter) :: iter

    patch_id = 1
    check_output = PRESENT(opt_output)
    IF (check_output) output = opt_output
    IF (PRESENT(opt_patch_id)) patch_id = opt_patch_id
    DO WHILE(iter%next())
      IF ( patch_id /= iter%cur%p%patch_id ) CYCLE
      IF (check_output) THEN
        IF (output .NEQV. iter%cur%p%loutput) CYCLE
      END IF
      element => find_list_element(iter%cur, vname, opt_hgrid, opt_cs, opt_output)
      IF (ASSOCIATED (element)) THEN
        IF (PRESENT(opt_list)) opt_list%p => iter%cur%p
        EXIT
      END IF
    END DO
  END FUNCTION vlr_find

  !>  Detailed print-out of variable groups.
  SUBROUTINE vlr_print_groups(idom, opt_latex_fmt, opt_reduce_trailing_num, opt_skip_trivial)
    INTEGER, INTENT(IN)           :: idom          !< domain ID
    LOGICAL, INTENT(IN), OPTIONAL :: opt_latex_fmt !< Flag: .TRUE., if output shall be formatted for LaTeX
    LOGICAL, INTENT(IN), OPTIONAL :: opt_reduce_trailing_num !< Flag: replace trailing numbers by "*"
    LOGICAL, INTENT(IN), OPTIONAL :: opt_skip_trivial        !< Flag: skip empty of single-entry groups
    CHARACTER(*), PARAMETER :: routine = modname//"::print_group_details"
    CHARACTER(LEN=vname_len), ALLOCATABLE :: grp_names(:), grp_vars(:), grp_vars_out(:)
    INTEGER :: ngrp_vars, ngrp_vars_out, ierrstat, i, j, lx_l
    LOGICAL :: latex_fmt, reduce_trailing_num, skip_trivial
    CHARACTER(LEN=2), PARAMETER :: lx = ' %'
    TYPE(t_vl_register_iter) :: iter

    latex_fmt = .FALSE.
    IF (PRESENT(opt_latex_fmt)) latex_fmt = opt_latex_fmt
    lx_l = MERGE(0, 2, latex_fmt)
    reduce_trailing_num = .TRUE.
    IF (PRESENT(opt_reduce_trailing_num)) reduce_trailing_num = opt_reduce_trailing_num
    skip_trivial = .TRUE.
    IF (PRESENT(opt_skip_trivial)) skip_trivial = opt_skip_trivial
    WRITE (0,*) " "
    WRITE (0,*) lx(1:lx_l)//"---------------------------------------"
    WRITE (0,'(a,i0,a)') lx(1:lx_l)//" Variable group info (for domain #", idom, "):"
    WRITE (0,*) lx(1:lx_l)//"---------------------------------------"
    IF (latex_fmt) THEN
      WRITE (0,*) "% "
      WRITE (0,*) "% LaTeX formatted output, requires suitable environment 'varlist' and"
      WRITE (0,*) "% macros 'varname' and 'grpname'."
    END IF
    WRITE (0,*) " "
    grp_names = var_groups_dyn%alphabetical_list()
    WRITE (0,*) lx(1:lx_l)//"List of groups:"
    CALL pretty_print_string_list(grp_names, opt_prefix=lx(1:lx_l)//"    ")
    ! temporary variables needed for variable group parsing
    i = 0
    DO WHILE(iter%next())
      DO j = 1, iter%cur%p%nvars
        i = i + MERGE(0, 1, iter%cur%p%vl(j)%p%info%lcontainer)
      END DO
    END DO
    ALLOCATE(grp_vars(i), grp_vars_out(i), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
    DO i = 1, SIZE(grp_names)
      ! for each group we collect the contained variables two times
      ! - the first time we collect *all* model level variables on
      ! the triangular grid, and the second time we collect only
      ! those variables which are available for output.
      CALL compose_grp_varlist(grp_names(i), .FALSE., grp_vars, ngrp_vars)
      IF ((skip_trivial) .AND. (ngrp_vars <= 1))  CYCLE
      WRITE (0,*) ' '
      IF (latex_fmt) THEN
        WRITE (0,"(3a)") "\begin{varlist}{\grpname{", TRIM(grp_names(i)), "}}"
      ELSE
        WRITE (0,"(3a)") 'GROUP "', TRIM(grp_names(i)), '":'
      END IF
      IF (ngrp_vars > 0) THEN
        CALL compose_grp_varlist(grp_names(i), .TRUE., grp_vars_out, ngrp_vars_out)
        CALL pretty_print_string_list(grp_vars_out(1:ngrp_vars_out), opt_prefix="    ")
        CALL difference(grp_vars, ngrp_vars, grp_vars_out, ngrp_vars_out)
        IF (ngrp_vars > 0) THEN
          WRITE (0,*) " "
          WRITE (0,*) "   Non-output variables:"
          WRITE (0,*) " "
          CALL pretty_print_string_list(grp_vars(1:ngrp_vars), opt_prefix="    ")
        END IF
      ELSE
        WRITE (0,*) "   -- empty --"
      END IF
      IF (latex_fmt) WRITE (0,*) "\end{varlist}"
    END DO
    WRITE (0,*) " "
  CONTAINS

    SUBROUTINE compose_grp_varlist(grpname, lout, vargrp, nvargrp)
      CHARACTER(*), INTENT(IN) :: grpname
      LOGICAL, INTENT(IN) :: lout
      CHARACTER(LEN=vname_len), INTENT(INOUT) :: vargrp(:)
      INTEGER, INTENT(INOUT) :: nvargrp
      INTEGER, ALLOCATABLE :: slen(:)
      INTEGER :: jj, kk
      LOGICAL :: lfound

      CALL vlr_group(grpname, vargrp, nvargrp, opt_dom_id=idom, &
        & loutputvars_only=lout, lremap_lonlat=.FALSE., opt_vlevel_type=level_type_ml)
      IF (nvargrp .LE. 0) RETURN
      CALL lowcase(vargrp(1:nvargrp))
      IF (reduce_trailing_num) THEN
        ALLOCATE(slen(nvargrp))
        slen = [(find_trailing_number(vargrp(jj))-1, jj = 1, nvargrp)]
        DO jj = 1, nvargrp
          IF (slen(jj) .EQ. -2) CYCLE
          ! replace trailing number if any other variable of "this kind" exists:
          lfound = .FALSE.
          DO kk = 1, nvargrp
            IF (lfound) EXIT
            IF (jj .EQ. kk) CYCLE
            lfound = (vargrp(jj)(1:slen(jj)) == vargrp(kk)(1:slen(kk)))
          END DO
          IF (lfound) vargrp(jj) = vargrp(jj)(1:slen(jj))//"*"
        END DO
        CALL remove_duplicates(vargrp(1:nvargrp), nvargrp)
      END IF
      CALL quicksort(vargrp(1:nvargrp))
      IF (latex_fmt) THEN
        DO jj = 1, nvargrp
          vargrp(jj) = "\varname{" // TRIM(vargrp(jj)) // "}"
        END DO
      END IF
    END SUBROUTINE compose_grp_varlist
  END SUBROUTINE vlr_print_groups

END MODULE mo_var_list_register_utils
