! @par Copyright and License
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
  USE mo_var_list_element, ONLY: level_type_ml
  USE mo_var_list,         ONLY: find_list_element, print_var_list, delete_list, &
    & t_var_list_ptr, t_list_element, append_list_element, sel_var_list
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
  USE mo_hash_table,            ONLY: t_HashTable, hashTable_make, t_HashIterator
#ifdef DEBUG_MVSTREAM
  USE mo_util_string,      ONLY: int2string
  USE mo_mpi,              ONLY: my_process_is_stdio
  USE self_assert,         ONLY: print_summary
#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: vl_register, vl_iter

  TYPE t_var_list_iterator
    PRIVATE
    TYPE(t_hashIterator) :: iter
    TYPE(t_var_list_ptr), PUBLIC :: cur
    LOGICAL :: anew = .TRUE.
  CONTAINS
    PROCEDURE, PUBLIC :: next => iter_next
    PROCEDURE, PUBLIC :: reset => iter_reset
  END TYPE t_var_list_iterator

  TYPE t_var_list_store
    PRIVATE
    LOGICAL, PUBLIC :: is_init = .FALSE.
    TYPE(t_hashTable), POINTER :: storage
    TYPE(t_var_list_iterator) :: it
  CONTAINS
    PROCEDURE, PUBLIC :: new => new_var_list
    PROCEDURE, PUBLIC :: delete => delete_var_list
    PROCEDURE, PUBLIC :: get => get_var_list
    PROCEDURE, PUBLIC :: find_var_all => find_var_all
    PROCEDURE, PUBLIC :: print_all => print_all_var_lists
    PROCEDURE, PUBLIC :: collect_group => collect_group
    PROCEDURE, PUBLIC :: print_group => print_group_details
    PROCEDURE, PUBLIC :: n_var => total_number_of_variables
    PROCEDURE, PUBLIC :: packer => varlistPacker
    PROCEDURE, PUBLIC :: new_var_ref => add_var_reference
  END TYPE t_var_list_store

  TYPE(t_var_list_store) :: vl_register
  TYPE(t_var_list_iterator) :: vl_iter
  INTEGER, SAVE :: vl_new_counter = 0

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_var_list_store'

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
    vl_new_counter = vl_new_counter + 1
    CALL put_var_list(tolower(vlname(1:vln_len)))
    ! set default list characteristics
    list%p%id       = vl_new_counter
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

  LOGICAL FUNCTION iter_next(this) RESULT(valid)
    CLASS(t_var_list_iterator), INTENT(INOUT) :: this
    CLASS(*), POINTER :: keyObj, valObj

    valid = .false.
    IF (vl_register%is_init) THEN
      IF (this%anew) CALL this%iter%init(vl_register%storage)
      valid = this%iter%nextEntry(keyObj, valObj)
      IF (valid) this%cur%p => sel_var_list(valObj)
      IF (.NOT.valid) CALL iter_reset(this)
    END IF
  END FUNCTION iter_next

  SUBROUTINE iter_reset(this)
    CLASS(t_var_list_iterator), INTENT(INOUT) :: this

    NULLIFY(this%cur%p)
    this%anew = .TRUE.
    CALL this%iter%reset()
  END SUBROUTINE iter_reset

  !------------------------------------------------------------------------------------------------
  ! @return total number of (non-container) variables
  INTEGER FUNCTION total_number_of_variables(this) RESULT(nvar)
    CLASS(t_var_list_store), INTENT(INOUT) :: this
    TYPE(t_list_element), POINTER :: element

    nvar = 0
    ! Note that there may be several variables with different time
    ! levels, we just add unconditionally all
    DO WHILE(iter_next(this%it))
      element => this%it%cur%p%first_list_element
      DO WHILE (ASSOCIATED(element))
        ! Do not count element if it is a container
        IF (.NOT.element%field%info%lcontainer) nvar = nvar + 1
        element => element%next_list_element
      ENDDO
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
      CALL delete_list(list)
      CALL this%storage%removeEntry(keyObj)
    END SUBROUTINE remove_var_list
  END SUBROUTINE delete_var_list

  !------------------------------------------------------------------------------------------------
  ! Delete all output var_lists
  SUBROUTINE delete_var_lists(this)
    CLASS(t_var_list_store), INTENT(INOUT) :: this
    
    DO WHILE(iter_next(this%it))
      CALL delete_var_list(this, this%it%cur)
    END DO
  END SUBROUTINE delete_var_lists

  !------------------------------------------------------------------------------------------------
  ! add supplementary fields to a different var list (eg. geopotential, surface pressure, ...)
  SUBROUTINE add_var_reference(this, to_var_list, vname, from_var_list, loutput, bit_precision, in_group)
    CLASS(t_var_list_store), INTENT(INOUT)   :: this
    TYPE(t_var_list_ptr), INTENT(inout)      :: to_var_list
    CHARACTER(len=*), INTENT(in)             :: vname, from_var_list
    LOGICAL,          INTENT(in),   OPTIONAL :: loutput, in_group(MAX_GROUPS)
    INTEGER,          INTENT(in),   OPTIONAL :: bit_precision
    TYPE(t_var_list_ptr) :: vlp_from
    TYPE(t_list_element),     POINTER :: n_list_e, o_list_e => NULL()
    !
    CALL get_var_list(this, vlp_from, from_var_list)
    IF (ASSOCIATED(vlp_from%p)) o_list_e => find_list_element(vlp_from, vname)
    IF (ASSOCIATED(o_list_e)) THEN
      CALL append_list_element(to_var_list, n_list_e)
      n_list_e%field                = o_list_e%field
      n_list_e%field%info%allocated = .FALSE.
      n_list_e%field%info%lrestart  = .FALSE.
      IF (PRESENT(loutput))       n_list_e%field%info%loutput     = loutput
      IF (PRESENT(bit_precision)) n_list_e%field%info%grib2%bits  = bit_precision
      IF (PRESENT(in_group))      n_list_e%field%info%in_group(:) = in_group(:)
    ENDIF
  END SUBROUTINE add_var_reference

  !------------------------------------------------------------------------------------------------
  ! print all var lists
  SUBROUTINE print_all_var_lists(this, lshort)
    CLASS(t_var_list_store), INTENT(INOUT) :: this
    LOGICAL, OPTIONAL, INTENT(IN) :: lshort

    DO WHILE(iter_next(this%it))
      CALL print_var_list(this%it%cur, lshort=lshort)
    END DO
  END SUBROUTINE print_all_var_lists

  !> Loops over all variables and collects the variables names
  !  corresponding to the group @p grp_name
  SUBROUTINE collect_group(this, grp_name, var_name, nvars,       &
    &                      loutputvars_only, lremap_lonlat, &
    &                      opt_vlevel_type, opt_dom_id,     &
    &                      opt_lquiet)
    CLASS(t_var_list_store), INTENT(INOUT)  :: this
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
    ! local variables
    CHARACTER(*), PARAMETER :: routine = modname//":collect_group", llmsg = " lon-lat"
    INTEGER                       :: grp_id, llmsg_len
    TYPE(t_list_element), POINTER :: element
    TYPE(t_var_metadata), POINTER :: info
    CHARACTER(LEN=vname_len)      :: name
    LOGICAL                       :: lquiet, verbose, skip

    nvars  = 0
    grp_id = var_groups_dyn%group_id(grp_name)
    lquiet = .FALSE.
    IF (PRESENT(opt_lquiet))  lquiet = opt_lquiet
    verbose = .NOT. lquiet

    ! loop over all variable lists and variables
    DO WHILE(iter_next(this%it))
      IF (PRESENT(opt_vlevel_type)) THEN
        IF (this%it%cur%p%vlevel_type /= opt_vlevel_type) CYCLE
      ENDIF
      IF (PRESENT(opt_dom_id)) THEN
        ! do not inspect variable list if its domain does not match:
        IF (this%it%cur%p%patch_id /= opt_dom_id)  CYCLE
      END IF

      element => this%it%cur%p%first_list_element
      LOOPVAR : DO WHILE (ASSOCIATED(element))
        info => element%field%info
        ! Do not inspect element if it is a container
        IF (.NOT. info%lcontainer .AND. info%in_group(grp_id)) THEN
          name = get_var_name(info)
          llmsg_len = 0
          ! Skip element if we need only output variables:
          skip = loutputvars_only .AND. &
            & (.NOT.info%loutput .OR. .NOT.this%it%cur%p%loutput)

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
    ENDDO ! i = 1, SIZE(var_lists)
    
    CALL remove_duplicates(var_name, nvars)
  END SUBROUTINE collect_group
  !-----------------------------------------------------------------------------
  !
  ! Find named list element accross all knows variable lists
  !
  FUNCTION find_var_all(this, vname, opt_patch_id, opt_hgrid, opt_list, opt_cs, opt_output) RESULT(element)
    CLASS(t_var_list_store), INTENT(INOUT) :: this
    CHARACTER(*), INTENT(in) :: vname
    INTEGER, OPTIONAL, INTENT(IN) :: opt_patch_id, opt_hgrid
    TYPE(t_var_list_ptr), OPTIONAL, INTENT(OUT) :: opt_list
    LOGICAL, OPTIONAL, INTENT(IN) :: opt_cs, opt_output
    TYPE(t_list_element), POINTER :: element
    INTEGER :: patch_id
    LOGICAL :: output, check_output

    patch_id = 1
    check_output = PRESENT(opt_output)
    IF (check_output) output = opt_output
    IF (PRESENT(opt_patch_id)) patch_id = opt_patch_id
    DO WHILE(iter_next(this%it))
      IF ( patch_id /= this%it%cur%p%patch_id ) CYCLE
      IF (check_output) THEN
        IF (output .NEQV. this%it%cur%p%loutput) CYCLE
      END IF
      element => find_list_element(this%it%cur, vname, opt_hgrid, opt_cs, opt_output)
      IF (ASSOCIATED (element)) THEN
#ifdef DEBUG_MVSTREAM
        if (my_process_is_stdio()) call &
            & print_summary('destination PATCHID:'//TRIM(int2string(patch_id)))
#endif
        IF (PRESENT(opt_list)) opt_list%p => this%it%cur%p
        EXIT
      END IF
    END DO
    CALL iter_reset(this%it)
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
    INTEGER :: iv, nv, nelems, nelems_all, patch_id, restart_type, vlevel_type, n, ierrstat, id
    INTEGER, ALLOCATABLE :: info_buf(:)
    TYPE(t_list_element), POINTER   :: element, newElement
    CHARACTER(LEN=128)              :: var_list_name
    CHARACTER(LEN=32)               :: model_type
    LOGICAL                         :: lrestart, loutput
    TYPE(t_var_list_ptr) :: vlp
    CHARACTER(LEN=*), PARAMETER :: routine = modname//':varlistPacker'

    IF (operation .EQ. kUnpackOp) CALL delete_var_lists(this)
    nv = 0
    IF (this%is_init) nv = this%storage%getEntryCount()
    CALL pmsg%packer(operation, nv)
    nelems_all = 0
    DO iv = 1, nv
      IF(operation .EQ. kPackOp) THEN
        IF (.NOT.iter_next(this%it)) CALL finish(routine, "inconsistency")
        vlp%p => this%it%cur%p
        ! copy the values needed for the new_var_list() CALL to local variables
        lrestart      = vlp%p%lrestart
        loutput       = vlp%p%loutput
        var_list_name = vlp%p%name
        id            = vlp%p%id
        model_type    = vlp%p%model_type
        patch_id      = vlp%p%patch_id
        restart_type  = vlp%p%restart_type
        vlevel_type   = vlp%p%vlevel_type
        element     =>  vlp%p%first_list_element
        nelems = 0
        DO WHILE (ASSOCIATED(element))
          IF(element%field%info%lrestart .OR. .NOT.restart_only) nelems = nelems + 1
          element => element%next_list_element
        END DO
        nelems_all = nelems_all + nelems
      END IF
      CALL pmsg%packer(operation, lrestart)
      CALL pmsg%packer(operation, loutput)
      CALL pmsg%packer(operation, var_list_name)
      CALL pmsg%packer(operation, id)
      CALL pmsg%packer(operation, model_type)
      CALL pmsg%packer(operation, patch_id)
      CALL pmsg%packer(operation, restart_type)
      CALL pmsg%packer(operation, vlevel_type)
      CALL pmsg%packer(operation, nelems)
      IF (restart_only .AND. .NOT. lrestart) CYCLE  ! transfer only a restart var_list
      IF (nelems .EQ. 0) CYCLE ! check if there are valid restart fields
      IF (operation .EQ. kPackOp) THEN
        element => vlp%p%first_list_element
        DO WHILE (ASSOCIATED(element))
          IF(element%field%info%lrestart .OR. .NOT.restart_only) THEN
            info_buf = var_metadata_toBinary(element%field%info)
            CALL pmsg%pack(info_buf)
          END IF
          element => element%next_list_element
        END DO
      ELSE
        ! create var list
        CALL new_var_list(this, vlp, var_list_name, patch_id=patch_id, &
          & model_type=model_type, restart_type=restart_type, &
          & vlevel_type=vlevel_type, lrestart=lrestart, loutput=loutput)
        vlp%p%id = id ! override !
        ! insert elements into var list
        DO n = 1, nelems
          ALLOCATE(newElement, STAT=ierrstat)
          IF (ierrstat /= SUCCESS) CALL finish(routine, "memory allocation failure")
          IF (n .EQ. 1) THEN
            vlp%p%first_list_element => newElement
          ELSE
            element%next_list_element => newElement
          END IF
          element => newElement
          NULLIFY(element%next_list_element ,element%field%r_ptr, element%field%s_ptr, &
            &     element%field%i_ptr, element%field%l_ptr)
          element%field%var_base_size = 0 ! Unknown here
          CALL pmsg%unpack(info_buf)
          element%field%info = var_metadata_fromBinary(info_buf)
        END DO
      END IF
    END DO
    IF (operation .EQ. kPackOp) THEN
      IF (iter_next(this%it)) CALL finish(routine, "inconsistency -- second kind")
    END IF 
    CALL pmsg%packer(operation, nelems_all)
    IF (PRESENT(nv_all)) nv_all = nelems_all
  END SUBROUTINE varlistPacker

  !>  Detailed print-out of variable groups.
  !
  SUBROUTINE print_group_details(this, idom, opt_latex_fmt, opt_reduce_trailing_num, opt_skip_trivial)
    CLASS(t_var_list_store), INTENT(INOUT)   :: this
    INTEGER, INTENT(IN)           :: idom          !< domain ID
    LOGICAL, INTENT(IN), OPTIONAL :: opt_latex_fmt !< Flag: .TRUE., if output shall be formatted for LaTeX
    LOGICAL, INTENT(IN), OPTIONAL :: opt_reduce_trailing_num !< Flag: replace trailing numbers by "*"
    LOGICAL, INTENT(IN), OPTIONAL :: opt_skip_trivial        !< Flag: skip empty of single-entry groups
    ! local variables
    CHARACTER(*), PARAMETER :: routine = modname//"::print_group_details"
    CHARACTER(len=vname_len), ALLOCATABLE :: group_names(:)
    CHARACTER(LEN=vname_len), ALLOCATABLE :: grp_vars(:), grp_vars_output(:)
    INTEGER,                    ALLOCATABLE :: slen(:)
    INTEGER                                 :: ngrp_vars, ngrp_vars_output, ierrstat, i, j, k, t
    LOGICAL                                 :: latex_fmt, reduce_trailing_num, skip_trivial, lfound

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
    i = total_number_of_variables(this)
    ALLOCATE(grp_vars_output(i), grp_vars(i), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
    
    DO i=1,SIZE(group_names)
      
      ! for each group we collect the contained variables two times
      ! - the first time we collect *all* model level variables on
      ! the triangular grid, and the second time we collect only
      ! those variables which are available for output.
      CALL collect_group(this, group_names(i), grp_vars, ngrp_vars,    &
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
               
        CALL collect_group(this, group_names(i), grp_vars_output, ngrp_vars_output,    &
          &               loutputvars_only = .TRUE.,               &
          &               lremap_lonlat    = .FALSE.,              &
          &               opt_vlevel_type  = level_type_ml,        &
          &               opt_dom_id       = idom,                 &
          &               opt_lquiet       = .TRUE.)
        
        DO j=1,ngrp_vars_output
          grp_vars_output(j) = tolower(grp_vars_output(j))
        END DO
        ! (optionally) replace trailing numbers by "*"
        IF (reduce_trailing_num) THEN
          ALLOCATE(slen(ngrp_vars_output))
          DO j=1,ngrp_vars_output
            slen(j) = find_trailing_number(grp_vars_output(j))
          END DO
          DO j=1,ngrp_vars_output
            t = slen(j)
            IF (t /= -1) THEN
              ! replace trailing number if any other variable of "this
              ! kind" exists:
              lfound = .FALSE.
              INNER_LOOP2 : DO k=1,ngrp_vars_output
                IF (j==k)  CYCLE INNER_LOOP2
                IF (grp_vars_output(j)(1:(t-1)) == grp_vars_output(k)(1:(slen(k)-1))) lfound = .TRUE.
              END DO INNER_LOOP2
              IF (lfound) THEN
                grp_vars_output(j) = grp_vars_output(j)(1:(t-1))//"*"
              END IF
            END IF
          END DO
          CALL remove_duplicates(grp_vars_output(1:ngrp_vars_output), ngrp_vars_output)
          DEALLOCATE(slen)
        END IF
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

END MODULE mo_var_list_register
