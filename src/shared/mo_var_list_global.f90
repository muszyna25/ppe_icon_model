! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
MODULE mo_var_list_global

  USE mo_var_groups,       ONLY: var_groups_dyn, MAX_GROUPS
  USE mo_var_metadata_types,ONLY: t_var_metadata, VINTP_TYPE_LIST, &
    & var_metadata_fromBinary, var_metadata_toBinary
  USE mo_var_list_element, ONLY: t_var_list_element, level_type_ml
  USE mo_var_list,         ONLY: find_list_element, print_var_list, get_var_name
  USE mo_linked_list,      ONLY: t_var_list, t_list_element,        &
       &                         delete_list, append_list_element
  USE mo_exception,        ONLY: message, finish
  USE mo_util_string,      ONLY: remove_duplicates,        &
    &                            pretty_print_string_list, tolower, &
    &                            difference, find_trailing_number
  USE mo_impl_constants,   ONLY: vname_len, VARNAME_LEN, MAX_TIME_LEVELS,   &
    &                            SUCCESS, TIMELEVEL_SUFFIX
  USE mo_cdi_constants, ONLY: GRID_UNSTRUCTURED_CELL, GRID_REGULAR_LONLAT
  USE mo_io_config,        ONLY: restart_file_type
  USE mo_packed_message,   ONLY: t_PackedMessage, kPackOp, kUnpackOp
  USE mo_util_sort,        ONLY: quicksort
  USE mo_key_value_store,  ONLY: t_key_value_store
#ifdef DEBUG_MVSTREAM
  USE mo_util_string,      ONLY: int2string
  USE mo_mpi,              ONLY: my_process_is_stdio
  USE self_assert,         ONLY: print_summary
#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: new_var_list              ! get a pointer to a new output var_list
  PUBLIC :: delete_var_list           ! delete an output var_list
  PUBLIC :: delete_var_lists          ! delete all output var_lists
  PUBLIC :: get_var_list              ! get a pointer to an existing output var_list
  PUBLIC :: print_all_var_lists
  PUBLIC :: collect_group
  PUBLIC :: var_lists                 ! vector of output var_lists
  PUBLIC :: nvar_lists                ! number of output var_lists defined so far
  PUBLIC :: get_all_var_names         ! obtain a list of variables names
  PUBLIC :: total_number_of_variables ! returns total number of defined variables
  PUBLIC :: find_var_global   ! find an element in the list
  PUBLIC :: varlistPacker
  PUBLIC :: print_group_details
  PUBLIC :: add_var_list_reference

  INTEGER :: nvar_lists     =   0      ! var_lists allocated so far
  TYPE(t_var_list), ALLOCATABLE, TARGET :: var_lists(:)  ! memory buffer array
  TYPE(t_key_value_store) :: var_lists_map
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_var_list'

CONTAINS
  !------------------------------------------------------------------------------------------------
  !
  ! Create a new memory buffer / output var_list
  ! Get a pointer to the new var_list
  !
  SUBROUTINE new_var_list (this_list, vlname, output_type, restart_type,      &
       &                   post_suf, rest_suf, init_suf, loutput, lrestart, &
       &                   linitial, patch_id, vlevel_type)
    TYPE(t_var_list), INTENT(inout)        :: this_list    ! anchor
    CHARACTER(len=*), INTENT(in)           :: vlname         ! name of output var_list
    INTEGER,          INTENT(in), OPTIONAL :: output_type, restart_type   ! 'GRIB' or 'NetCDF'
    CHARACTER(len=*), INTENT(in), OPTIONAL :: post_suf, rest_suf, init_suf ! suffix of output/restart/initial file
    LOGICAL,          INTENT(in), OPTIONAL :: loutput, lrestart, linitial  ! in standard output/restart/initial file
    INTEGER,          INTENT(in), OPTIONAL :: patch_id     ! patch ID
    INTEGER,          INTENT(in), OPTIONAL :: vlevel_type  ! 1/2/3 for model/pres./height levels
    INTEGER :: i, ierr, nvl_used
    TYPE(t_var_list), ALLOCATABLE :: tmp(:)
    !
    CALL message('','')
    CALL message('','adding new var_list '//TRIM(vlname))
    IF (.NOT.var_lists_map%is_init) CALL var_lists_map%init(.FALSE.)
    CALL var_lists_map%get(vlname, i, ierr)
    IF (ierr .EQ. 0) CALL finish('new_var_list', ' >'//TRIM(vlname)//'< already in use.')
    this_list%p => NULL()
    nvl_used = var_lists_map%getEntryCount()
    IF (nvar_lists .NE. nvl_used) CALL finish('new_var_list', "inconsistent element counts")
    IF (.NOT.ALLOCATED(var_lists)) THEN
      ALLOCATE(var_lists(12))
    ELSE IF (nvl_used .GE. SIZE(var_lists)) THEN
      ALLOCATE(tmp(SIZE(var_lists) + 4))
      FORALL(i = 1:nvl_used) tmp(i)%p => var_lists(i)%p
      CALL MOVE_ALLOC(tmp, var_lists) 
    END IF
    nvl_used = nvl_used + 1
    nvar_lists = nvl_used
    ALLOCATE(var_lists(nvl_used)%p)
    this_list%p => var_lists(nvl_used)%p
    CALL var_lists_map%put(vlname, nvar_lists)
    ! set default list characteristics
    this_list%p%name     = vlname
    this_list%p%post_suf = '_'//TRIM(vlname)
    this_list%p%rest_suf = this_list%p%post_suf
    this_list%p%init_suf = this_list%p%post_suf
    this_list%p%loutput  = .TRUE.
    ! set non-default list characteristics
    CALL set_var_list(this_list, output_type=output_type,                &
      & restart_type=restart_type, post_suf=post_suf, rest_suf=rest_suf, &
      & init_suf=init_suf, loutput=loutput, lrestart=lrestart,           &
      & linitial=linitial, patch_id=patch_id, vlevel_type=vlevel_type)
  END SUBROUTINE new_var_list
  !------------------------------------------------------------------------------------------------
  !
  ! Get a reference to a memory buffer/output var_list
  !
  SUBROUTINE get_var_list (this_list, vlname)
    TYPE(t_var_list), POINTER, INTENT(OUT) :: this_list ! pointer
    CHARACTER(len=*), INTENT(IN) :: vlname      ! name of output var_list
    INTEGER :: i, ierr

    CALL var_lists_map%get(vlname, i, ierr)
    NULLIFY(this_list)
    IF (ierr .EQ. 0) this_list => var_lists(i)
  END SUBROUTINE get_var_list

  !------------------------------------------------------------------------------------------------
  !
  ! @return total number of (non-container) variables
  !
  FUNCTION total_number_of_variables()
    INTEGER :: total_number_of_variables
    INTEGER :: i
    TYPE(t_list_element), POINTER :: element

    total_number_of_variables = 0
    ! Note that there may be several variables with different time
    ! levels, we just add unconditionally all
    DO i = 1,nvar_lists
      element => var_lists(i)%p%first_list_element
      LOOPVAR : DO WHILE (ASSOCIATED(element))
        ! Do not count element if it is a container
        total_number_of_variables = total_number_of_variables &
             + MERGE(1, 0, .NOT. element%field%info%lcontainer)
        element => element%next_list_element
      ENDDO LOOPVAR ! loop over vlist "i"
    ENDDO ! i = 1,nvar_lists
  END FUNCTION total_number_of_variables

  !------------------------------------------------------------------------------------------------
  !
  ! Get a list of variable names matching a given criterion.
  !
  SUBROUTINE get_all_var_names (varlist, ivar, opt_vlevel_type,                          &
    &                           opt_vert_intp_type,                                      &
    &                           opt_hor_intp_type, opt_lcontainer,                       &
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
    FORALL(i = 1:SIZE(varlist)) varlist(i) = ' '
    ! default values for some criteria:
    lcontainer = .FALSE.
    loutput = .FALSE.
    hor_intp_type_match = -1
    vert_intp_type = .FALSE.
    IF (PRESENT(opt_lcontainer)) lcontainer = opt_lcontainer
    loutput_matters = PRESENT(opt_loutput)
    IF (loutput_matters) loutput = opt_loutput
    hor_intp_matters= PRESENT(opt_hor_intp_type)
    IF (hor_intp_matters) hor_intp_type_match = opt_hor_intp_type
    IF (PRESENT(opt_vert_intp_type)) vert_intp_type = opt_vert_intp_type
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
      IF (loutput_matters) THEN
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
  ! Delete an output var_list, nullify the associated pointer
  !
  SUBROUTINE delete_var_list(this_list)
    TYPE(t_var_list) :: this_list
    !
    IF (ASSOCIATED(this_list%p)) THEN
      CALL var_lists_map%remove(TRIM(this_list%p%name))
      CALL delete_list(this_list)
      DEALLOCATE(this_list%p)
    ENDIF
  END SUBROUTINE delete_var_list
  !------------------------------------------------------------------------------------------------
  !
  ! Delete all output var_lists
  !
  SUBROUTINE delete_var_lists
    TYPE(t_var_list), POINTER :: this_list
    INTEGER :: i
    !
    DO i = 1, nvar_lists
      this_list => var_lists(i)
      IF (ASSOCIATED(this_list%p)) THEN
        CALL delete_list(this_list)
        DEALLOCATE(this_list%p)
      END IF
    END DO
    CALL var_lists_map%destruct()
  END SUBROUTINE delete_var_lists

  !================================================================================================
  !------------------------------------------------------------------------------------------------
  !
  ! add supplementary fields to a different var list (eg. geopotential, surface pressure, ...)
  !
  SUBROUTINE add_var_list_reference (to_var_list, vname, from_var_list, loutput, bit_precision, in_group)
    TYPE(t_var_list), INTENT(inout)          :: to_var_list
    CHARACTER(len=*), INTENT(in)             :: vname, from_var_list
    LOGICAL,          INTENT(in),   OPTIONAL :: loutput, in_group(MAX_GROUPS)
    INTEGER,          INTENT(in),   OPTIONAL :: bit_precision
    TYPE(t_list_element),     POINTER :: n_list_e, o_list_e => NULL()
    INTEGER :: i, ierr
    !
    CALL var_lists_map%get(from_var_list, i, ierr)
    IF (ierr .EQ. 0) o_list_e => find_list_element(var_lists(i), vname)
    IF (ASSOCIATED(o_list_e)) THEN
      CALL append_list_element(to_var_list, n_list_e)
      n_list_e%field                = o_list_e%field
      n_list_e%field%info%allocated = .FALSE.
      n_list_e%field%info%lrestart  = .FALSE.
      IF (PRESENT(loutput))       n_list_e%field%info%loutput     = loutput
      IF (PRESENT(bit_precision)) n_list_e%field%info%grib2%bits  = bit_precision
      IF (PRESENT(in_group))      n_list_e%field%info%in_group(:) = in_group(:)
    ENDIF
  END SUBROUTINE add_var_list_reference
  !------------------------------------------------------------------------------------------------
  !
  ! print all var lists
  !
  SUBROUTINE print_all_var_lists(lshort)
    LOGICAL, OPTIONAL, INTENT(IN) :: lshort
    INTEGER :: i

    DO i = 1, nvar_lists
      CALL print_var_list(var_lists(i), lshort=lshort)
    END DO
  END SUBROUTINE print_all_var_lists

  !> Loops over all variables and collects the variables names
  !  corresponding to the group @p grp_name
  !
  SUBROUTINE collect_group(grp_name, var_name, nvars,       &
    &                      loutputvars_only, lremap_lonlat, &
    &                      opt_vlevel_type, opt_dom_id,     &
    &                      opt_lquiet)
    CHARACTER(LEN=*),           INTENT(IN)    :: grp_name
    CHARACTER(LEN=VARNAME_LEN), INTENT(INOUT) :: var_name(:)
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
  !-----------------------------------------------------------------------------
  !
  ! Find named list element accross all knows variable lists
  !
  FUNCTION find_var_global(vname, opt_patch_id, opt_hgrid, opt_caseInsensitive,opt_returnList) RESULT(element)
    CHARACTER(len=*),   INTENT(in) :: vname
    INTEGER, OPTIONAL              :: opt_patch_id, opt_hgrid
    LOGICAL, OPTIONAL              :: opt_caseInsensitive
    TYPE(t_var_list), POINTER, OPTIONAL :: opt_returnList
    TYPE(t_list_element), POINTER :: element
    INTEGER :: i, patch_id

    patch_id = 1
    IF (PRESENT(opt_patch_id)) patch_id = opt_patch_id
    DO i=1,nvar_lists
      IF ( patch_id /= var_lists(i)%p%patch_id ) CYCLE
      element => find_list_element(var_lists(i),vname,opt_hgrid,opt_caseInsensitive)
      IF (ASSOCIATED (element)) THEN
#ifdef DEBUG_MVSTREAM
        if (my_process_is_stdio()) call &
            & print_summary('destination PATCHID:'//TRIM(int2string(patch_id)))
#endif
        IF (PRESENT(opt_returnList)) opt_returnList => var_lists(i)
        RETURN
      END IF
    END DO
  END FUNCTION find_var_global

  !-----------------------------------------------------------------------------
  !
  ! (Un)pack the var_lists
  ! This IS needed for the restart modules that need to communicate the
  ! var_lists from the worker PEs to dedicated restart PEs.
  SUBROUTINE varlistPacker(operation, packedMessage, restart_only, nv_all)
    INTEGER, INTENT(IN) :: operation
    TYPE(t_PackedMessage), INTENT(INOUT) :: packedMessage
    LOGICAL, INTENT(IN) :: restart_only
    INTEGER, INTENT(OUT), OPTIONAL :: nv_all
    INTEGER :: iv, nv, nelems, nelems_all, patch_id, restart_type, vlevel_type, n, ierrstat
    INTEGER, ALLOCATABLE :: info_buf(:)
    TYPE(t_list_element), POINTER   :: element, newElement
    TYPE(t_var_list)                :: p_var_list
    CHARACTER(LEN=128)              :: var_list_name
    CHARACTER(LEN=32)               :: model_type
    LOGICAL                         :: lrestart
    CHARACTER(LEN=*), PARAMETER :: routine = modname//':varlistPacker'

    IF(operation == kUnpackOp) CALL delete_var_lists
    nv = nvar_lists
    nelems_all = 0
    CALL packedMessage%packer(operation, nv)
    DO iv = 1, nv
      IF(operation == kPackOp) THEN
        ! copy the values needed for the new_var_list() CALL to local variables
        lrestart = var_lists(iv)%p%lrestart
        var_list_name = var_lists(iv)%p%name
        model_type = var_lists(iv)%p%model_type
        patch_id = var_lists(iv)%p%patch_id
        restart_type = var_lists(iv)%p%restart_type
        vlevel_type = var_lists(iv)%p%vlevel_type
        element => var_lists(iv)%p%first_list_element
        nelems = 0
        DO WHILE (ASSOCIATED(element))
          IF(element%field%info%lrestart .OR. .NOT.restart_only) nelems = nelems+1
          element => element%next_list_element
        END DO
        nelems_all = nelems_all + nelems
      END IF
      CALL packedMessage%packer(operation, lrestart)
      CALL packedMessage%packer(operation, var_list_name)
      CALL packedMessage%packer(operation, model_type)
      CALL packedMessage%packer(operation, patch_id)
      CALL packedMessage%packer(operation, restart_type)
      CALL packedMessage%packer(operation, vlevel_type)
      CALL packedMessage%packer(operation, nelems)
      IF (restart_only .AND. .NOT. lrestart) CYCLE  ! transfer only a restart var_list
      IF (nelems == 0) CYCLE ! check if there are valid restart fields
      IF (operation == kPackOp) THEN
        element => var_lists(iv)%p%first_list_element
        DO WHILE (ASSOCIATED(element))
          IF(element%field%info%lrestart .OR. .NOT.restart_only) THEN
            info_buf = var_metadata_toBinary(element%field%info)
            CALL packedMessage%pack(info_buf)
          END IF
          element => element%next_list_element
        END DO
      END IF
      IF (operation == kUnpackOp) THEN
        ! create var list
        CALL new_var_list(p_var_list, var_list_name, patch_id=patch_id, &
          & restart_type=restart_type, vlevel_type=vlevel_type, lrestart=lrestart)
        p_var_list%p%model_type = TRIM(model_type)
        NULLIFY(p_var_list%p%first_list_element)
        ! insert elements into var list
        DO n = 1, nelems
          ALLOCATE(newElement, STAT=ierrstat)
          IF (ierrstat /= SUCCESS) CALL finish(routine, "memory allocation failure")
          IF (n == 1) THEN   ! the first element pointer needs to be stored IN a different variable than the later pointers (there are no double pointers IN FORTRAN...)
            p_var_list%p%first_list_element => newElement
          ELSE
            element%next_list_element => newElement
          END IF
          element => newElement
          NULLIFY(element%next_list_element ,element%field%r_ptr, element%field%s_ptr, &
            &     element%field%i_ptr, element%field%l_ptr)
          element%field%var_base_size = 0 ! Unknown here
          CALL packedMessage%unpack(info_buf)
          element%field%info = var_metadata_fromBinary(info_buf)
        END DO
      END IF
    END DO
    CALL packedMessage%packer(operation, nelems_all)
    IF (PRESENT(nv_all)) nv_all = nelems_all
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
               
        CALL collect_group(group_names(i), grp_vars_output, ngrp_vars_output,    &
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

END MODULE mo_var_list_global
