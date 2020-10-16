! @par Copyright and Li/cense
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
MODULE mo_var_list_register
  USE mo_var_groups,       ONLY: var_groups_dyn, MAX_GROUPS
  USE mo_var_metadata_types,ONLY: t_var_metadata, var_metadata_fromBinary, &
    & var_metadata_toBinary
  USE mo_var_metadata,     ONLY: get_var_name
  USE mo_var, ONLY: level_type_ml, t_var
  USE mo_var_list,         ONLY: find_list_element, t_var_list_ptr, sel_var_list
  USE mo_exception,        ONLY: message, finish
  USE mo_util_string,      ONLY: remove_duplicates, pretty_print_string_list, &
    &                            tolower, difference, find_trailing_number
  USE mo_impl_constants,   ONLY: vname_len, SUCCESS
  USE mo_cdi_constants,    ONLY: GRID_UNSTRUCTURED_CELL, GRID_REGULAR_LONLAT
  USE mo_io_config,        ONLY: restart_file_type
  USE mo_packed_message,   ONLY: t_PackedMessage, kPackOp, kUnpackOp
  USE mo_util_sort,        ONLY: quicksort
  USE mo_util_texthash,    ONLY: text_hash, text_isEqual
#ifdef __PGI
  USE mo_util_texthash,    ONLY: t_char_workaround
#endif
  USE mo_hash_table,       ONLY: t_HashTable, hashTable_make, t_HashIterator

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: vl_register, t_vl_register_iter

  TYPE t_vl_register_iter
    PRIVATE
    TYPE(t_hashIterator) :: hash_iter
#ifdef HAVE_CDI_ORDERING_DEFECT
    INTEGER :: cur_idx = -1
#endif
    TYPE(t_var_list_ptr), PUBLIC :: cur
    LOGICAL :: anew = .TRUE.
  CONTAINS
#ifdef HAVE_CDI_ORDERING_DEFECT
    PROCEDURE, PUBLIC :: next_workaround => iter_next_workaround
#endif
    PROCEDURE, PUBLIC :: next => iter_next
  END TYPE t_vl_register_iter

#ifdef HAVE_CDI_ORDERING_DEFECT
  INTEGER, SAVE :: vl_new_counter = 0,  &
    &              vl_del_counter = 0,  &
    &              vl_nc_lastup   = -1, &
    &              vl_dc_lastup   = -1
  TYPE(t_var_list_ptr), ALLOCATABLE :: vl_vec_reordered(:)
#endif

  TYPE t_var_list_store
    PRIVATE
    LOGICAL, PUBLIC :: is_init = .FALSE.
    TYPE(t_hashTable), POINTER :: storage
  CONTAINS
    PROCEDURE, PUBLIC :: new => new_var_list
    PROCEDURE, PUBLIC :: delete => delete_var_list
    PROCEDURE, PUBLIC :: get => get_var_list
    PROCEDURE, PUBLIC, NOPASS :: find_var_all => find_var_all
    PROCEDURE, PUBLIC, NOPASS :: print_all => print_all_var_lists
    PROCEDURE, PUBLIC, NOPASS :: collect_group => collect_group
    PROCEDURE, PUBLIC, NOPASS :: print_group => print_group_details
    PROCEDURE, PUBLIC, NOPASS :: n_var => total_number_of_variables
    PROCEDURE, PUBLIC :: packer => varlistPacker
    PROCEDURE, PUBLIC :: new_var_ref => add_var_reference
  END TYPE t_var_list_store

  TYPE(t_var_list_store) :: vl_register

  CHARACTER(*), PARAMETER :: modname = 'mo_var_list_store'

CONTAINS
  !------------------------------------------------------------------------------------------------
  ! Create a new memory buffer / output var_list
  ! Get a pointer to the new var_list
  SUBROUTINE new_var_list(this, list, vlname, output_type, restart_type,    &
    & post_suf, rest_suf, init_suf, loutput, lrestart, linitial, patch_id,          &
    & vlevel_type, model_type, filename, compression_type)
    CLASS(t_var_list_store), INTENT(INOUT) :: this
    TYPE(t_var_list_ptr), INTENT(OUT) :: list    ! anchor
    CHARACTER(*), INTENT(IN) :: vlname         ! name of output var_list
    INTEGER, INTENT(IN), OPTIONAL :: output_type, restart_type, patch_id, vlevel_type, compression_type
    CHARACTER(len=*), INTENT(IN), OPTIONAL :: post_suf, rest_suf, init_suf, model_type, filename
    LOGICAL,          INTENT(IN), OPTIONAL :: loutput, lrestart, linitial  ! in standard output/restart/initial file
    INTEGER :: vln_len

    vln_len = LEN_TRIM(vlname)
    CALL message('','')
    CALL message('','adding new var_list '//vlname(1:vln_len))
    IF (.NOT.this%is_init) THEN
      this%storage => hashTable_make(text_hash, text_isEqual)
      this%is_init = .TRUE.
    END IF
    CALL get_var_list(this, list, vlname)
    IF (ASSOCIATED(list%p)) CALL finish('new_var_list', ' >'//vlname(1:vln_len)//'< already in use.')
    ALLOCATE(list%p)
#ifdef HAVE_CDI_ORDERING_DEFECT
    vl_new_counter = vl_new_counter + 1
#endif
    CALL put_var_list(tolower(vlname(1:vln_len)))
    ! set default list characteristics
#ifdef HAVE_CDI_ORDERING_DEFECT
    list%p%id       = vl_new_counter
#endif
    list%p%name     = vlname(1:vln_len)
    list%p%post_suf = '_'//vlname(1:vln_len)
    list%p%rest_suf = list%p%post_suf
    list%p%init_suf = list%p%post_suf
    ! set non-default list characteristics
    list%p%restart_type = restart_file_type
    IF (PRESENT(output_type))      list%p%output_type      = output_type
    IF (PRESENT(restart_type))     list%p%restart_type     = restart_type
    IF (PRESENT(post_suf))         list%p%post_suf         = post_suf
    IF (PRESENT(rest_suf))         list%p%rest_suf         = rest_suf
    IF (PRESENT(init_suf))         list%p%init_suf         = init_suf
    IF (PRESENT(loutput))          list%p%loutput          = loutput
    IF (PRESENT(lrestart))         list%p%lrestart         = lrestart
    IF (PRESENT(linitial))         list%p%linitial         = linitial
    IF (PRESENT(patch_id))         list%p%patch_id         = patch_id
    IF (PRESENT(vlevel_type))      list%p%vlevel_type      = vlevel_type
    IF (PRESENT(filename))         list%p%filename         = filename
    IF (PRESENT(compression_type)) list%p%compression_type = compression_type
    IF (PRESENT(model_type))       list%p%model_type       = model_type
  CONTAINS

    SUBROUTINE put_var_list(key) ! Cray needs this to happen in an extra subroutine
      CHARACTER(*), INTENT(IN) :: key
      CLASS(*), POINTER :: keyObj, valObj
#ifdef __PGI
      TYPE(t_char_workaround), POINTER :: key_p

      ALLOCATE(key_p)
      ALLOCATE(CHARACTER(LEN=LEN(key)) :: key_p%c)
      WRITE(key_p%c, "(a)") key
      keyObj => key_p
#else
      ALLOCATE(keyObj, SOURCE=key)
#endif
      valObj => list%p
      CALL this%storage%setEntry(keyObj, valObj)
    END SUBROUTINE put_var_list
  END SUBROUTINE new_var_list

  !------------------------------------------------------------------------------------------------
  ! Get a reference to a memory buffer/output var_list
  SUBROUTINE get_var_list(this, list, vlname)
    CLASS(t_var_list_store), INTENT(IN) :: this
    TYPE(t_var_list_ptr), INTENT(OUT) :: list ! pointer
    CHARACTER(*), INTENT(IN) :: vlname
    CLASS(*), POINTER :: valObj

    valObj => get_valObj(tolower(vlname)) ! Cray needs this to happen in an extra function
    IF (ASSOCIATED(valObj)) list%p => sel_var_list(valObj)
  CONTAINS

    FUNCTION get_valObj(key) RESULT(valObj)
      CHARACTER(*), INTENT(IN), TARGET :: key
      CLASS(*), POINTER :: keyObj, valObj

      keyObj => key
      valObj => this%storage%getEntry(keyObj)
    END FUNCTION get_valObj
  END SUBROUTINE get_var_list

#ifdef HAVE_CDI_ORDERING_DEFECT
  LOGICAL FUNCTION iter_next_workaround(this) RESULT(valid)
    CLASS(t_vl_register_iter), INTENT(INOUT) :: this
    CLASS(*), POINTER :: keyObj, valObj

    valid = .false.
    IF (vl_register%is_init) THEN
      IF (vl_new_counter .NE. vl_nc_lastup .OR. &
        & vl_del_counter .NE. vl_dc_lastup) THEN
        IF (.NOT.this%anew) &
          & CALL finish("vl_iter_next", "impossible, go away!")
        CALL construct_ordered_vl_vector()
        vl_nc_lastup = vl_new_counter
        vl_dc_lastup = vl_del_counter
      END IF
      IF (this%anew) THEN
        this%cur_idx = 0
        this%anew = .false.
      END IF
      this%cur_idx = this%cur_idx + 1
      IF (this%cur_idx .LE. SIZE(vl_vec_reordered)) THEN
        valid = .true.
        this%cur%p => vl_vec_reordered(this%cur_idx)%p
      ELSE
        NULLIFY(this%cur%p)
        this%anew = .true.
      END IF
    END IF
  END FUNCTION iter_next_workaround

  SUBROUTINE construct_ordered_vl_vector()
    TYPE(t_var_list_ptr), ALLOCATABLE :: tmp(:)
    INTEGER :: max_id, i_vl, i, n_vl
    TYPE(t_vl_register_iter) :: iter

    IF (ALLOCATED(vl_vec_reordered)) DEALLOCATE(vl_vec_reordered)
    n_vl = 0
    max_id = 0
    DO WHILE(iter_next(iter))
      n_vl = n_vl + 1
      max_id = MAX(max_id, iter%cur%p%id)
    ENDDO
    ALLOCATE(tmp(max_id), vl_vec_reordered(n_vl))
    DO WHILE(iter_next(iter))
      tmp(iter%cur%p%id)%p => iter%cur%p
    ENDDO
    i_vl = 0
    DO i = 1, max_id
      IF (ASSOCIATED(tmp(i)%p)) THEN
        i_vl = i_vl + 1
        vl_vec_reordered(i_vl)%p => tmp(i)%p
      END IF
    END DO
    IF (i_vl .NE. n_vl) CALL finish("construct_ordered_vl_vector", "inconsistency")
  END SUBROUTINE construct_ordered_vl_vector
#endif

  LOGICAL FUNCTION iter_next(this) RESULT(valid)
    CLASS(t_vl_register_iter), INTENT(INOUT) :: this
    CLASS(*), POINTER :: keyObj, valObj

    valid = .false.
    IF (vl_register%is_init) THEN
      IF (this%anew) CALL this%hash_iter%init(vl_register%storage)
      valid = this%hash_iter%nextEntry(keyObj, valObj)
      IF (valid) THEN
        this%cur%p => sel_var_list(valObj)
      ELSE ! prepare for the next loop ...
        NULLIFY(this%cur%p)   
        this%anew = .true.
        CALL this%hash_iter%reset()
      END IF
    END IF
  END FUNCTION iter_next

  !------------------------------------------------------------------------------------------------
  ! @return total number of (non-container) variables
  INTEGER FUNCTION total_number_of_variables() RESULT(nvar)
    TYPE(t_vl_register_iter) :: iter
    INTEGER :: i

    nvar = 0
    ! Note that there may be several variables with different time
    ! levels, we just add unconditionally all
    DO WHILE(iter_next(iter))
      DO i = 1, iter%cur%p%nvars
        IF (.NOT.iter%cur%p%vl(i)%p%info%lcontainer) &
          nvar = nvar + 1
      END DO
    ENDDO
  END FUNCTION total_number_of_variables

  !------------------------------------------------------------------------------------------------
  ! Delete an output var_list, nullify the associated pointer
  SUBROUTINE delete_var_list(this, list)
    CLASS(t_var_list_store), INTENT(INOUT) :: this
    TYPE(t_var_list_ptr), INTENT(INOUT) :: list

    IF (ASSOCIATED(list%p)) & ! Cray needs this to happen in an extra subroutine
      & CALL remove_var_list(tolower(list%p%name))
  CONTAINS

    SUBROUTINE remove_var_list(key)
      CHARACTER(*), INTENT(IN), TARGET :: key
      CLASS(*), POINTER :: keyObj

      keyObj => key
      CALL list%delete()
      CALL this%storage%removeEntry(keyObj)
#ifdef HAVE_CDI_ORDERING_DEFECT
      vl_del_counter = vl_del_counter + 1
#endif
    END SUBROUTINE remove_var_list
  END SUBROUTINE delete_var_list

  !------------------------------------------------------------------------------------------------
  ! Delete all output var_lists
  SUBROUTINE delete_var_lists(this)
    CLASS(t_var_list_store), INTENT(INOUT) :: this
    TYPE(t_vl_register_iter) :: iter

    DO WHILE(iter_next(iter))
      CALL delete_var_list(this, iter%cur)
    END DO
  END SUBROUTINE delete_var_lists

  !------------------------------------------------------------------------------------------------
  ! add supplementary fields to a different var list (eg. geopotential, surface pressure, ...)
  SUBROUTINE add_var_reference(this, to_vl, vname, from_vl, loutput, bit_precision, in_group)
    CLASS(t_var_list_store), INTENT(INOUT) :: this
    TYPE(t_var_list_ptr), INTENT(INOUT) :: to_vl
    CHARACTER(*), INTENT(IN) :: vname, from_vl
    LOGICAL, INTENT(in), OPTIONAL :: loutput, in_group(MAX_GROUPS)
    INTEGER, INTENT(in), OPTIONAL :: bit_precision
    TYPE(t_var_list_ptr) :: vlp_from
    TYPE(t_var), POINTER :: n_e, o_e => NULL()

    CALL get_var_list(this, vlp_from, from_vl)
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
  END SUBROUTINE add_var_reference

  !------------------------------------------------------------------------------------------------
  ! print all var lists
  SUBROUTINE print_all_var_lists(lshort)
    LOGICAL, OPTIONAL, INTENT(IN) :: lshort
    TYPE(t_vl_register_iter) :: iter

    DO WHILE(iter_next(iter))
      CALL iter%cur%print(lshort=lshort)
    END DO
  END SUBROUTINE print_all_var_lists

  !> Loops over all variables and collects the variables names
  !  corresponding to the group @p grp_name
  SUBROUTINE collect_group(grp_name, var_name, nvars,       &
    &                      loutputvars_only, lremap_lonlat, &
    &                      opt_vlevel_type, opt_dom_id,     &
    &                      opt_lquiet)
    CHARACTER(*),             INTENT(IN)    :: grp_name
    CHARACTER(LEN=vname_len), INTENT(INOUT) :: var_name(:)
    INTEGER,                  INTENT(OUT)   :: nvars
    ! loutputvars_only: If set to .TRUE. all variables in the group
    ! which have the the loutput flag equal to .FALSE. are skipped.
    LOGICAL,                  INTENT(IN)    :: loutputvars_only
    ! lremap_lonlat: If set to .TRUE. only variables in the group
    ! which can be interpolated onto lon-lat grids are considered.
    LOGICAL,                  INTENT(IN)    :: lremap_lonlat
    ! 1: model levels, 2: pressure levels, 3: height level
    INTEGER, OPTIONAL,        INTENT(IN)    :: opt_vlevel_type, opt_dom_id
    LOGICAL, OPTIONAL,        INTENT(IN)    :: opt_lquiet
    TYPE(t_vl_register_iter) :: iter
    CHARACTER(*), PARAMETER :: routine = modname//":collect_group", llmsg = " lon-lat"
    INTEGER :: grp_id, llmsg_len, i
    TYPE(t_var_metadata), POINTER :: info
    CHARACTER(LEN=vname_len) :: vname
    LOGICAL                       :: lquiet, verbose, skip

    nvars  = 0
    grp_id = var_groups_dyn%group_id(grp_name)
    lquiet = .FALSE.
    IF (PRESENT(opt_lquiet))  lquiet = opt_lquiet
    verbose = .NOT. lquiet
    ! loop over all variable lists and variables
#ifdef HAVE_CDI_ORDERING_DEFECT 
    DO WHILE(iter_next_workaround(iter))
#else
    DO WHILE(iter_next(iter))
#endif
      IF (PRESENT(opt_vlevel_type)) THEN
        IF (iter%cur%p%vlevel_type /= opt_vlevel_type) CYCLE
      ENDIF
      IF (PRESENT(opt_dom_id)) THEN
        ! do not inspect variable list if its domain does not match:
        IF (iter%cur%p%patch_id /= opt_dom_id)  CYCLE
      END IF
      LOOPVAR : DO i = 1, iter%cur%p%nvars
        info => iter%cur%p%vl(i)%p%info
        ! Do not inspect element if it is a container
        IF (.NOT. info%lcontainer .AND. info%in_group(grp_id)) THEN
          vname = get_var_name(info)
          llmsg_len = 0
          ! Skip element if we need only output variables:
          skip = loutputvars_only .AND. &
            & (.NOT.info%loutput .OR. .NOT.iter%cur%p%loutput)

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
            var_name(nvars) = vname
          ELSE IF (verbose) THEN
            CALL message(routine, "Skipping variable "//TRIM(vname)//" for " &
                 //llmsg(1:llmsg_len)//"output.")
          END IF
        END IF
      ENDDO LOOPVAR ! loop over vlist "i"
    ENDDO ! i = 1, SIZE(var_lists)
    
    CALL remove_duplicates(var_name, nvars)
  END SUBROUTINE collect_group

  !-----------------------------------------------------------------------------
  ! Find named list element accross all knows variable lists
  FUNCTION find_var_all(vname, opt_patch_id, opt_hgrid, opt_list, opt_cs, opt_output) RESULT(element)
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
    DO WHILE(iter_next(iter))
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
  END FUNCTION find_var_all

  !-----------------------------------------------------------------------------
  ! (Un)pack the var_lists
  ! This IS needed for the restart modules that need to communicate the
  ! var_lists from the worker PEs to dedicated restart PEs.
  SUBROUTINE varlistPacker(this, operation, pmsg, restart_only, nv_all)
    CLASS(t_var_list_store), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: operation
    TYPE(t_PackedMessage), INTENT(INOUT) :: pmsg
    LOGICAL, INTENT(IN) :: restart_only
    INTEGER, INTENT(OUT), OPTIONAL :: nv_all
    INTEGER :: ivl, nvl, iv, nv, nv_al, patch_id, restart_type, vlevel_type, ierrstat
#ifdef HAVE_CDI_ORDERING_DEFECT
    INTEGER :: id
#endif
    INTEGER, ALLOCATABLE :: info_buf(:)
    TYPE(t_var), POINTER :: elem
    CHARACTER(LEN=128) :: var_list_name
    CHARACTER(LEN=32) :: model_type
    LOGICAL :: lrestart, loutput
    TYPE(t_var_list_ptr) :: vlp
    TYPE(t_vl_register_iter) :: iter
    CHARACTER(LEN=*), PARAMETER :: routine = modname//':varlistPacker'

    IF (operation .EQ. kUnpackOp) CALL delete_var_lists(this)
    nvl = 0
    IF (this%is_init) nvl = this%storage%getEntryCount()
    CALL pmsg%packer(operation, nvl)
    nv_al = 0
    DO ivl = 1, nvl
      IF(operation .EQ. kPackOp) THEN
        IF (.NOT.iter_next(iter)) CALL finish(routine, "inconsistency")
        vlp%p => iter%cur%p
        ! copy the values needed for the new_var_list() CALL to local variables
        lrestart      = vlp%p%lrestart
        loutput       = vlp%p%loutput
        var_list_name = vlp%p%name
#ifdef HAVE_CDI_ORDERING_DEFECT
        id            = vlp%p%id
#endif
        model_type    = vlp%p%model_type
        patch_id      = vlp%p%patch_id
        restart_type  = vlp%p%restart_type
        vlevel_type   = vlp%p%vlevel_type
        nv = vlp%p%nvars
        IF (restart_only) THEN
          nv = 0
          DO iv = 1, vlp%p%nvars
            IF(vlp%p%vl(iv)%p%info%lrestart) nv = nv + 1
          END DO
        END IF
        nv_al = nv_al + nv
      END IF
      CALL pmsg%packer(operation, lrestart)
      CALL pmsg%packer(operation, loutput)
      CALL pmsg%packer(operation, var_list_name)
#ifdef HAVE_CDI_ORDERING_DEFECT
      CALL pmsg%packer(operation, id)
#endif
      CALL pmsg%packer(operation, model_type)
      CALL pmsg%packer(operation, patch_id)
      CALL pmsg%packer(operation, restart_type)
      CALL pmsg%packer(operation, vlevel_type)
      CALL pmsg%packer(operation, nv)
      IF (restart_only .AND. .NOT. lrestart) CYCLE  ! transfer only a restart var_list
      IF (nv .EQ. 0) CYCLE ! check if there are valid restart fields
      IF (operation .EQ. kPackOp) THEN
        DO iv = 1, vlp%p%nvars
          elem => vlp%p%vl(iv)%p
          IF(elem%info%lrestart .OR. .NOT.restart_only) THEN
            info_buf = var_metadata_toBinary(elem%info)
            CALL pmsg%pack(info_buf)
            CALL pmsg%pack(vlp%p%tl(iv))
            CALL pmsg%pack(vlp%p%hgrid(iv))
            CALL pmsg%pack(vlp%p%key(iv))
            CALL pmsg%pack(vlp%p%key_notl(iv))
            CALL pmsg%pack(vlp%p%lout(iv))
          END IF
        END DO
      ELSE
        ! create var list
        CALL new_var_list(this, vlp, var_list_name, patch_id=patch_id, &
          & model_type=model_type, restart_type=restart_type, &
          & vlevel_type=vlevel_type, lrestart=lrestart, loutput=loutput)
#ifdef HAVE_CDI_ORDERING_DEFECT
        vlp%p%id = id ! override !
#endif
        vlp%p%nvars = nv
        ! insert elements into var list
        ALLOCATE(vlp%p%vl(nv), vlp%p%tl(nv), vlp%p%hgrid(nv), &
          & vlp%p%key(nv), vlp%p%key_notl(nv), vlp%p%lout(nv))
        DO iv = 1, nv
          ALLOCATE(vlp%p%vl(iv)%p, STAT=ierrstat)
          IF (ierrstat /= SUCCESS) CALL finish(routine, "memory allocation failure")
          elem => vlp%p%vl(iv)%p
          NULLIFY(elem%r_ptr, elem%s_ptr, elem%i_ptr, elem%l_ptr)
          elem%var_base_size = 0 ! Unknown here
          CALL pmsg%unpack(info_buf)
          elem%info = var_metadata_fromBinary(info_buf)
          CALL pmsg%unpack(vlp%p%tl(iv))
          CALL pmsg%unpack(vlp%p%hgrid(iv))
          CALL pmsg%unpack(vlp%p%key(iv))
          CALL pmsg%unpack(vlp%p%key_notl(iv))
          CALL pmsg%unpack(vlp%p%lout(iv))
        END DO
      END IF
    END DO
    IF (operation .EQ. kPackOp) THEN
      IF (iter_next(iter)) CALL finish(routine, "inconsistency -- second kind")
    END IF 
    CALL pmsg%packer(operation, nv_al)
    IF (PRESENT(nv_all)) nv_all = nv_al
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
    CHARACTER(LEN=vname_len), ALLOCATABLE :: grp_names(:), grp_vars(:), grp_vars_out(:)
    INTEGER, ALLOCATABLE :: slen(:)
    INTEGER :: ngrp_vars, ngrp_vars_out, ierrstat, i, j, k, t
    LOGICAL :: latex_fmt, reduce_trailing_num, skip_trivial, lfound

    latex_fmt = .FALSE.
    IF (PRESENT(opt_latex_fmt)) latex_fmt = opt_latex_fmt 
    reduce_trailing_num = .FALSE.
    IF (PRESENT(opt_reduce_trailing_num)) reduce_trailing_num = opt_reduce_trailing_num 
    skip_trivial = .FALSE.
    IF (PRESENT(opt_skip_trivial)) skip_trivial = opt_skip_trivial 
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
    grp_names = var_groups_dyn%alphabetical_list()
    IF (latex_fmt) THEN
      WRITE (0,*) "% List of groups:"
      CALL pretty_print_string_list(grp_names, opt_prefix=" %    ")
    ELSE
      WRITE (0,*) "List of groups:"
      CALL pretty_print_string_list(grp_names, opt_prefix="    ")
    END IF
    ! temporary variables needed for variable group parsing
    i = total_number_of_variables()
    ALLOCATE(grp_vars_out(i), grp_vars(i), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
    DO i = 1, SIZE(grp_names)
      ! for each group we collect the contained variables two times
      ! - the first time we collect *all* model level variables on
      ! the triangular grid, and the second time we collect only
      ! those variables which are available for output.
      CALL collect_group(grp_names(i), grp_vars, ngrp_vars,    &
        &               loutputvars_only = .FALSE.,              &
        &               lremap_lonlat    = .FALSE.,              &
        &               opt_vlevel_type  = level_type_ml,        &
        &               opt_dom_id       = idom,                 &
        &               opt_lquiet       = .TRUE.)
      IF (ngrp_vars > 0) THEN
        DO j=1,ngrp_vars
          grp_vars(j) = tolower(grp_vars(j))
        END DO
        ! (optionally) replace trailing numbers by "*"
        IF (reduce_trailing_num) THEN
          ALLOCATE(slen(ngrp_vars))
          DO j=1,ngrp_vars
            slen(j) = find_trailing_number(grp_vars(j))
          END DO
          DO j=1,ngrp_vars
            t = slen(j)
            IF (t /= -1) THEN
              ! replace trailing number if any other variable of "this
              ! kind" exists:
              lfound = .FALSE.
              INNER_LOOP1 : DO k=1,ngrp_vars
                IF (j==k)  CYCLE INNER_LOOP1
                IF (grp_vars(j)(1:(t-1)) == grp_vars(k)(1:(slen(k)-1))) lfound = .TRUE.
              END DO INNER_LOOP1
              IF (lfound) THEN
                grp_vars(j) = grp_vars(j)(1:(t-1))//"*"
              END IF
            END IF
          END DO
          CALL remove_duplicates(grp_vars(1:ngrp_vars), ngrp_vars)
          DEALLOCATE(slen)
        END IF
        CALL quicksort(grp_vars(1:ngrp_vars))
        IF ((skip_trivial) .AND. (ngrp_vars <= 1))  CYCLE
        CALL collect_group(grp_names(i), grp_vars_out, ngrp_vars_out,    &
          &               loutputvars_only = .TRUE.,               &
          &               lremap_lonlat    = .FALSE.,              &
          &               opt_vlevel_type  = level_type_ml,        &
          &               opt_dom_id       = idom,                 &
          &               opt_lquiet       = .TRUE.)
        
        DO j=1,ngrp_vars_out
          grp_vars_out(j) = tolower(grp_vars_out(j))
        END DO
        ! (optionally) replace trailing numbers by "*"
        IF (reduce_trailing_num) THEN
          ALLOCATE(slen(ngrp_vars_out))
          DO j=1,ngrp_vars_out
            slen(j) = find_trailing_number(grp_vars_out(j))
          END DO
          DO j=1,ngrp_vars_out
            t = slen(j)
            IF (t /= -1) THEN
              ! replace trailing number if any other variable of "this
              ! kind" exists:
              lfound = .FALSE.
              INNER_LOOP2 : DO k=1,ngrp_vars_out
                IF (j==k)  CYCLE INNER_LOOP2
                IF (grp_vars_out(j)(1:(t-1)) == grp_vars_out(k)(1:(slen(k)-1))) lfound = .TRUE.
              END DO INNER_LOOP2
              IF (lfound) grp_vars_out(j) = grp_vars_out(j)(1:(t-1))//"*"
            END IF
          END DO
          CALL remove_duplicates(grp_vars_out(1:ngrp_vars_out), ngrp_vars_out)
          DEALLOCATE(slen)
        END IF
        CALL quicksort(grp_vars_out(1:ngrp_vars_out))
        IF (latex_fmt) THEN
          DO j=1,ngrp_vars
            grp_vars(j) = "\varname{"//TRIM(grp_vars(j))//"}"
          END DO
          DO j=1,ngrp_vars_out
            grp_vars_out(j) = "\varname{"//TRIM(grp_vars_out(j))//"}"
          END DO
        END IF
        WRITE (0,*) ' '
        IF (latex_fmt) THEN
          WRITE (0,*) "\begin{varlist}{\grpname{"//TRIM(grp_names(i))//"}}"
        ELSE
          WRITE (0,*) 'GROUP "', TRIM(grp_names(i)), '":'
        END IF
        CALL pretty_print_string_list(grp_vars_out(1:ngrp_vars_out), opt_prefix="    ")
        CALL difference(grp_vars, ngrp_vars, grp_vars_out, ngrp_vars_out)
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
          WRITE (0,*) "\begin{varlist}{\grpname{"//TRIM(grp_names(i))//"}}"
        ELSE
          WRITE (0,*) 'GROUP "', TRIM(grp_names(i)), '":'
        END IF
        WRITE (0,*) "   -- empty --"
      END IF
      IF (latex_fmt) WRITE (0,*) "\end{varlist}"
    END DO
    WRITE (0,*) " "
    DEALLOCATE(grp_names, grp_vars, grp_vars_out, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')
  END SUBROUTINE print_group_details

END MODULE mo_var_list_register
