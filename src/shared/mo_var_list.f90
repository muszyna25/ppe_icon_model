MODULE mo_var_list

  USE mo_kind,             ONLY: wp, i8
  USE mo_cdi_constants,    ONLY: DATATYPE_FLT64,                    &
       &                         DATATYPE_INT32,                    &
       &                         DATATYPE_INT8,                     &
       &                         TIME_VARIABLE,                     &
       &                         TSTEP_INSTANT,                     &
       &                         GRID_UNSTRUCTURED_CELL
  USE mo_cf_convention,    ONLY: t_cf_var
  USE mo_grib2,            ONLY: t_grib2_var
  USE mo_var_metadata,     ONLY: t_var_metadata, t_union_vals,      &
    &                            t_tracer_meta,                     &
    &                            t_vert_interp_meta,                &
    &                            t_hor_interp_meta,                 &
    &                            VARNAME_LEN, VAR_GROUPS,           &
    &                            VINTP_TYPE_LIST,                   &
    &                            t_post_op_meta, POST_OP_NONE
  USE mo_var_list_element, ONLY: t_var_list_element
  USE mo_linked_list,      ONLY: t_var_list, t_list_element,        &
       &                         new_list, delete_list,             &
       &                         append_list_element,               &
       &                         find_list_element,                 &
       &                         delete_list_element 
  USE mo_exception,        ONLY: message, message_text, finish
  USE mo_util_hash,        ONLY: util_hashword
  USE mo_util_string,      ONLY: remove_duplicates, toupper
  USE mo_impl_constants,   ONLY: VINTP_METHOD_LIN,                  &
    &                            HINTP_TYPE_LONLAT_RBF,             &
    &                            max_var_lists, vname_len,          &
    &                            STR_HINTP_TYPE, MAX_TIME_LEVELS
  USE mo_advection_config, ONLY: t_advection_config
  USE mo_fortran_tools,    ONLY: assign_if_present, t_ptr_2d3d
  
  USE mo_art_config,       ONLY: t_art_config, art_config

  IMPLICIT NONE

  PRIVATE
  
  PUBLIC :: new_var_list              ! get a pointer to a new output var_list
  PUBLIC :: delete_var_list           ! delete an output var_list
  PUBLIC :: delete_var_lists          ! delete all output var_lists
  PUBLIC :: get_var_list              ! get a pointer to an existing output var_list
  PUBLIC :: set_var_list              ! set default parameters of an output var_list
  PUBLIC :: create_tracer_metadata    ! create metadata for tracer variables
  PUBLIC :: create_vert_interp_metadata ! create metadata for vertical interpolation
  PUBLIC :: create_hor_interp_metadata  ! create metadata for horizontal interpolation
  PUBLIC :: print_var_list
  PUBLIC :: print_memory_use

  PUBLIC :: default_var_list_settings ! set default settings for a whole list

  PUBLIC :: var_lists                 ! vector of output var_lists
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

  PUBLIC :: total_number_of_variables ! returns total number of defined variables
  PUBLIC :: groups                    ! group array constructor
  PUBLIC :: post_op
  PUBLIC :: vintp_type_id             ! mapping: string specifier of vert. interp. typer -> ID
  PUBLIC :: vintp_types               ! vertical interpolation type (array) constructor
  PUBLIC :: collect_group             ! obtain variables in a group
  
  PUBLIC :: add_tracer_ref            ! add new tracer component

 INTERFACE add_var  ! create a new list entry
    MODULE PROCEDURE add_var_list_element_r5d
    MODULE PROCEDURE add_var_list_element_r4d
    MODULE PROCEDURE add_var_list_element_r3d 
    MODULE PROCEDURE add_var_list_element_r2d 
    MODULE PROCEDURE add_var_list_element_r1d 
    MODULE PROCEDURE add_var_list_element_i5d
    MODULE PROCEDURE add_var_list_element_i4d
    MODULE PROCEDURE add_var_list_element_i3d 
    MODULE PROCEDURE add_var_list_element_i2d 
    MODULE PROCEDURE add_var_list_element_i1d 
    MODULE PROCEDURE add_var_list_element_l5d
    MODULE PROCEDURE add_var_list_element_l4d
    MODULE PROCEDURE add_var_list_element_l3d 
    MODULE PROCEDURE add_var_list_element_l2d 
    MODULE PROCEDURE add_var_list_element_l1d 
  END INTERFACE add_var

  INTERFACE add_ref
    MODULE PROCEDURE add_var_list_reference_r3d
    MODULE PROCEDURE add_var_list_reference_r2d
  END INTERFACE add_ref


  INTERFACE add_tracer_ref
    MODULE PROCEDURE add_var_list_reference_tracer
  END INTERFACE add_tracer_ref


  INTERFACE get_var  ! obtain reference to a list entry
    MODULE PROCEDURE get_var_list_element_r5d
    MODULE PROCEDURE get_var_list_element_r4d
    MODULE PROCEDURE get_var_list_element_r3d
    MODULE PROCEDURE get_var_list_element_r2d
    MODULE PROCEDURE get_var_list_element_r1d
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
  END INTERFACE struct_assign_if_present
  
  INTEGER,                  SAVE :: nvar_lists     =   0      ! var_lists allocated so far
  !
  TYPE(t_var_list), TARGET, SAVE :: var_lists(max_var_lists)  ! memory buffer array
  !
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
    this_list%p%post_suf = '_'//TRIM(name)
    this_list%p%rest_suf = this_list%p%post_suf
    this_list%p%init_suf = this_list%p%post_suf
    this_list%p%loutput  = .TRUE.
    !
    ! set non-default list characteristics
    !
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
    TYPE(t_var_list), POINTER    :: this_list ! pointer
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
  ! @return total number of variables
  !
  FUNCTION total_number_of_variables()
    INTEGER :: total_number_of_variables
    ! local variables
    INTEGER :: i
    TYPE(t_list_element), POINTER :: element
    
    total_number_of_variables = 0
    !- loop over variables

    ! Note that there may be several variables with different time
    ! levels, we just add unconditionally all
    DO i = 1,nvar_lists
      element => NULL()
      LOOPVAR : DO
        IF(.NOT.ASSOCIATED(element)) THEN
          element => var_lists(i)%p%first_list_element
        ELSE
          element => element%next_list_element
        ENDIF
        IF(.NOT.ASSOCIATED(element)) EXIT LOOPVAR
        ! Do not inspect element if it is a container
        IF (element%field%info%lcontainer) CYCLE LOOPVAR
        
        total_number_of_variables = total_number_of_variables + 1
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
    INTEGER                       :: i, ivintp_type
    LOGICAL                       :: lcontainer
    TYPE(t_list_element), POINTER :: element
    TYPE(t_var_metadata), POINTER :: info

    ! clear result
    DO i=1,SIZE(varlist)
      varlist(i) = ' '
    END DO

    ! default values for some criteria:
    lcontainer = .FALSE.
    CALL assign_if_present(lcontainer, opt_lcontainer)

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
      element => NULL()
      LOOPVAR : DO
        IF(.NOT.ASSOCIATED(element)) THEN
          element => var_lists(i)%p%first_list_element
        ELSE
          element => element%next_list_element
        ENDIF
        IF(.NOT.ASSOCIATED(element)) EXIT LOOPVAR

        info => element%field%info
        ! Do not inspect element if it is a container
        IF (info%lcontainer .NEQV. lcontainer) CYCLE LOOPVAR
        ! Do not inspect element if "loutput=.false."
        IF (PRESENT(opt_loutput)) THEN
          IF (opt_loutput .NEQV. info%loutput) CYCLE LOOPVAR
        END IF
        ! Do not inspect element if it does not contain info for
        ! vertical interpolation
        IF (PRESENT(opt_vert_intp_type)) THEN
          LOOP_VINTP_TYPES : DO ivintp_type=1,SIZE(VINTP_TYPE_LIST)
            IF (opt_vert_intp_type(ivintp_type) .AND. &
              & .NOT. info%vert_interp%vert_intp_type(ivintp_type)) CYCLE LOOPVAR
          END DO LOOP_VINTP_TYPES
        END IF
        ! Do not inspect element if it does not contain info for
        ! horizontal interpolation
        IF (PRESENT(opt_hor_intp_type)) THEN
          IF (info%hor_interp%hor_intp_type /= opt_hor_intp_type) CYCLE LOOPVAR
        END IF

        ! Check for time level suffix:
        ivar = ivar+1
        varlist(ivar) = TRIM(get_var_name(element%field))
      ENDDO LOOPVAR ! loop over vlist "i"
    ENDDO LOOP_VARLISTS ! i = 1,nvar_lists

    CALL remove_duplicates(varlist, ivar)

  END SUBROUTINE get_all_var_names


  !------------------------------------------------------------------------------------------------
  !> @return Plain variable name (i.e. without time level suffix ".TL")
  !
  FUNCTION get_var_name(var)
    CHARACTER(LEN=VARNAME_LEN) :: get_var_name
    TYPE(t_var_list_element)   :: var
    ! local variable
    INTEGER :: idx

    idx = INDEX(var%info%name,'.TL')
    IF (idx==0) THEN
      get_var_name = TRIM(var%info%name)
    ELSE
      get_var_name = TRIM(var%info%name(1:idx-1))
    END IF
  END FUNCTION get_var_name


  !------------------------------------------------------------------------------------------------
  !> @return time level (extracted from time level suffix ".TL") or "-1"
  !
  FUNCTION get_var_timelevel(var)
    INTEGER :: get_var_timelevel
    TYPE(t_var_list_element)   :: var
    ! local variable
    CHARACTER(LEN=*), PARAMETER :: routine = 'mo_var_list:get_var_timelevel'
    INTEGER :: idx

    idx = INDEX(var%info%name,'.TL')
    IF (idx == 0) THEN
      get_var_timelevel = -1
      RETURN
    END IF

    ! Get time level
    get_var_timelevel = ICHAR(var%info%name(idx+3:idx+3)) - ICHAR('0')
    IF(get_var_timelevel<=0 .OR. get_var_timelevel>max_time_levels) &
      CALL finish(routine, 'Illegal time level in '//TRIM(var%info%name))
  END FUNCTION get_var_timelevel


  !------------------------------------------------------------------------------------------------
  !> @return tile index (extracted from tile index suffix "t_") or "-1"
  !
  FUNCTION get_var_tileidx(varname)
    INTEGER :: get_var_tileidx
    CHARACTER(LEN=*) :: varname
    ! local variable
    CHARACTER(LEN=*), PARAMETER :: routine = 'mo_var_list:get_var_tileidx'
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
  FUNCTION default_var_list_metadata(this_list) RESULT(this_info)
    !
    TYPE(t_var_metadata)         :: this_info        ! memory info structure
    !
    TYPE(t_var_list), INTENT(in) :: this_list        ! output var_list    
    !
    this_info%key                 = 0    
    this_info%name                = ''
    !
    this_info%cf                  = t_cf_var('', '', '', -1)
    this_info%grib2               = t_grib2_var(-1, -1, -1, -1, -1, -1)
    !
    this_info%allocated           = .FALSE.
    this_info%ndims               = 0
    this_info%used_dimensions(:)  = 0
    !
    ! RJ: Set default loutput to .TRUE., regardless of this_list%p%loutput
    this_info%loutput             = .TRUE.
    this_info%isteptype           = TSTEP_INSTANT
    this_info%resetval            = t_union_vals( 0.0_wp, 0, .FALSE.)
    this_info%lrestart            = this_list%p%lrestart
    this_info%lrestart_cont       = .FALSE.
    this_info%lrestart_read       = .FALSE.
    this_info%lmiss               = this_list%p%lmiss
    this_info%missval             = t_union_vals( 0.0_wp, 0, .FALSE.)
    this_info%initval             = t_union_vals( 0.0_wp, 0, .FALSE.)
    !
    this_info%lcontainer          = .FALSE.
    this_info%lcontained          = .FALSE.
    this_info%ncontained          = 0
    !
    this_info%hgrid               = -1
    this_info%vgrid               = -1
    !
    this_info%tlev_source         = 0
    !
    this_info%cdiVarID            = -1
    this_info%cdiVarID_2          = -1
    this_info%cdiGridID           = -1
    this_info%cdiZaxisID          = -1
    this_info%cdiDataType         = -1
    !
    this_info%tracer              = create_tracer_metadata()
    this_info%vert_interp         = create_vert_interp_metadata()
    this_info%hor_interp          = create_hor_interp_metadata()
    !
    this_info%post_op             = post_op()
    !
    this_info%in_group(:)         = groups()
    !
    this_info%cdiTimeID           = TIME_VARIABLE
    this_info%l_pp_scheduler_task = 0

  END FUNCTION default_var_list_metadata
  !------------------------------------------------------------------------------------------------
  !
  ! Set tracer parameters of list element already created
  ! (public routine. Can be used for two things:
  ! 1.) default settings: If used without any argument, it gives back a variable
  !     of type t_tracer_meta, containing the default settings.
  ! 2.) Setting of metadata: If used with arguments, it gives back a variable 
  !     of type t_tracer_meta, containing the default settings except for those components 
  !     which are given in the argument list.
  !
  ! Comment by DR: Maybe for the future one could define different sets of default values
  ! for different groups of ART species.
  ! 
  FUNCTION create_tracer_metadata(lis_tracer,tracer_class,                              &
    &                            ihadv_tracer, ivadv_tracer, lturb_tracer,              &
    &                            lsed_tracer, ldep_tracer, lconv_tracer,                &
    &                            lwash_tracer,rdiameter_tracer, rrho_tracer,            &
    &                            halflife_tracer,imis_tracer) RESULT(tracer_meta) 

    LOGICAL, INTENT(IN), OPTIONAL :: lis_tracer      ! this is a tracer field (TRUE/FALSE)
    CHARACTER(len=*), INTENT(IN), OPTIONAL :: tracer_class  ! type of tracer (cloud, volcash, radioact,...) 
    INTEGER, INTENT(IN), OPTIONAL :: ihadv_tracer    ! method for horizontal transport
    INTEGER, INTENT(IN), OPTIONAL :: ivadv_tracer    ! method for vertical transport
    LOGICAL, INTENT(IN), OPTIONAL :: lturb_tracer    ! turbulent transport (TRUE/FALSE)
    LOGICAL, INTENT(IN), OPTIONAL :: lsed_tracer     ! sedimentation (TRUE/FALSE)
    LOGICAL, INTENT(IN), OPTIONAL :: ldep_tracer     ! dry deposition (TRUE/FALSE)  
    LOGICAL, INTENT(IN), OPTIONAL :: lconv_tracer    ! convection  (TRUE/FALSE)
    LOGICAL, INTENT(IN), OPTIONAL :: lwash_tracer    ! washout (TRUE/FALSE)
    REAL(wp), INTENT(IN), OPTIONAL :: rdiameter_tracer! particle diameter in m
    REAL(wp), INTENT(IN), OPTIONAL :: rrho_tracer     ! particle density in kg m^-3
    REAL(wp), INTENT(IN), OPTIONAL :: halflife_tracer       ! radioactive half-life in s^-1
    INTEGER, INTENT(IN), OPTIONAL :: imis_tracer         ! IMIS number

    TYPE(t_tracer_meta) :: tracer_meta               ! tracer metadata

    ! lis_tracer
    IF ( PRESENT(lis_tracer) ) THEN
      tracer_meta%lis_tracer = lis_tracer
    ELSE
      tracer_meta%lis_tracer = .FALSE.
    ENDIF

    ! tracer_class
    IF ( PRESENT(tracer_class) ) THEN
      tracer_meta%tracer_class = tracer_class
    ELSE
      tracer_meta%tracer_class = "cloud"  ! USE cloud as standard
    ENDIF

    ! ihadv_tracer
    IF ( PRESENT(ihadv_tracer) ) THEN
      tracer_meta%ihadv_tracer = ihadv_tracer
    ELSE
      tracer_meta%ihadv_tracer = 2
    ENDIF

    ! ivadv_tracer
    IF ( PRESENT(ivadv_tracer) ) THEN
      tracer_meta%ivadv_tracer = ivadv_tracer
    ELSE
      tracer_meta%ivadv_tracer = 3
    ENDIF

    ! lturb_tracer  
    IF ( PRESENT(lturb_tracer) ) THEN
      tracer_meta%lturb_tracer = lturb_tracer
    ELSE
      tracer_meta%lturb_tracer = .FALSE.
    ENDIF

    ! lsed_tracer
    IF ( PRESENT(lsed_tracer) ) THEN
      tracer_meta%lsed_tracer = lsed_tracer
    ELSE
      tracer_meta%lsed_tracer = .FALSE.
    ENDIF

    ! ldep_tracer
    IF ( PRESENT(ldep_tracer) ) THEN
      tracer_meta%ldep_tracer = ldep_tracer
    ELSE
      tracer_meta%ldep_tracer = .FALSE.
    ENDIF

    ! lconv_tracer
    IF ( PRESENT(lconv_tracer) ) THEN
      tracer_meta%lconv_tracer = lconv_tracer
    ELSE
      tracer_meta%lconv_tracer = .FALSE.
    ENDIF

    ! lwash_tracer
    IF ( PRESENT(lwash_tracer) ) THEN
      tracer_meta%lwash_tracer = lwash_tracer
    ELSE
      tracer_meta%lwash_tracer = .FALSE.
    ENDIF

    ! rdiameter_tracer
    IF ( PRESENT(rdiameter_tracer) ) THEN
      tracer_meta%rdiameter_tracer = rdiameter_tracer
    ELSE
      tracer_meta%rdiameter_tracer = 1.0e-6_wp  ! particle diameter in m
    ENDIF

    ! rrho_tracer
    IF ( PRESENT(rrho_tracer) ) THEN
      tracer_meta%rrho_tracer = rrho_tracer
    ELSE
      tracer_meta%rrho_tracer = 1000.0_wp       ! particle density in kg m^-3
    ENDIF

    ! halflife_tracer
    IF ( PRESENT(halflife_tracer) ) THEN
      tracer_meta%halflife_tracer = halflife_tracer
    ELSE
      tracer_meta%halflife_tracer = 0.0_wp       ! half-life in s^-1
    ENDIF

    ! imis_tracer
    IF ( PRESENT(imis_tracer) ) THEN
      tracer_meta%imis_tracer = imis_tracer
    ELSE
      tracer_meta%imis_tracer = -999       ! IMIS number
    ENDIF

  END FUNCTION create_tracer_metadata

  
  !------------------------------------------------------------------------------------------------
  !
  ! Quasi-constructor for vertical interpolation meta data
  ! 
  ! Fills data structure with default values (unless set otherwise).
  FUNCTION create_vert_interp_metadata(vert_intp_type, vert_intp_method,                     &
    &  l_hires_intp, l_restore_fricred, l_loglin, l_extrapol, l_satlimit, l_restore_pbldev,  &
    &  l_pd_limit, l_restore_sfcinv, l_hires_corr, lower_limit, extrapol_dist)               &
    RESULT(vert_interp_meta)

    TYPE(t_vert_interp_meta) :: vert_interp_meta    
    LOGICAL, INTENT(IN), OPTIONAL      :: &
      &  vert_intp_type(SIZE(VINTP_TYPE_LIST))
    INTEGER, INTENT(IN), OPTIONAL      :: &
      &  vert_intp_method
    LOGICAL, INTENT(IN), OPTIONAL      :: &
      &  l_hires_intp, l_restore_fricred, l_loglin, &
      &  l_extrapol, l_satlimit, l_restore_pbldev,  &
      &  l_pd_limit, l_restore_sfcinv, l_hires_corr
    REAL(wp), INTENT(IN), OPTIONAL     :: &
      &  lower_limit, extrapol_dist

    ! set default values
    vert_interp_meta%vert_intp_type(:) = .FALSE.
    vert_interp_meta%vert_intp_method  = VINTP_METHOD_LIN
    vert_interp_meta%l_hires_intp      = .FALSE.
    vert_interp_meta%l_restore_fricred = .FALSE.
    vert_interp_meta%l_loglin          = .FALSE.
    vert_interp_meta%l_extrapol        = .TRUE.
    vert_interp_meta%l_satlimit        = .FALSE.
    vert_interp_meta%l_restore_pbldev  = .FALSE.
    vert_interp_meta%l_pd_limit        = .FALSE.
    vert_interp_meta%l_restore_sfcinv  = .FALSE.
    vert_interp_meta%l_hires_corr      = .FALSE.
    vert_interp_meta%lower_limit       = 0._wp
    vert_interp_meta%extrapol_dist     = 500._wp
    ! supersede with user definitions
    CALL assign_if_present(vert_interp_meta%vert_intp_type     , vert_intp_type    )
    CALL assign_if_present(vert_interp_meta%vert_intp_method   , vert_intp_method  )
    CALL assign_if_present(vert_interp_meta%l_hires_intp       , l_hires_intp      )
    CALL assign_if_present(vert_interp_meta%l_restore_fricred  , l_restore_fricred )
    CALL assign_if_present(vert_interp_meta%l_loglin           , l_loglin          )
    CALL assign_if_present(vert_interp_meta%l_extrapol         , l_extrapol        )
    CALL assign_if_present(vert_interp_meta%l_satlimit         , l_satlimit        )
    CALL assign_if_present(vert_interp_meta%l_restore_pbldev   , l_restore_pbldev  )
    CALL assign_if_present(vert_interp_meta%l_pd_limit         , l_pd_limit        )
    CALL assign_if_present(vert_interp_meta%l_restore_sfcinv   , l_restore_sfcinv  )
    CALL assign_if_present(vert_interp_meta%l_hires_corr       , l_hires_corr      )
    CALL assign_if_present(vert_interp_meta%lower_limit        , lower_limit       )
    CALL assign_if_present(vert_interp_meta%extrapol_dist      , extrapol_dist     )

  END FUNCTION create_vert_interp_metadata


  !------------------------------------------------------------------------------------------------
  !
  ! Quasi-constructor for horizontal interpolation meta data
  ! 
  ! Fills data structure with default values (unless set otherwise).
  FUNCTION create_hor_interp_metadata(hor_intp_type, lonlat_id)    &
    RESULT(hor_interp_meta)

    TYPE(t_hor_interp_meta) :: hor_interp_meta    
    INTEGER, INTENT(IN), OPTIONAL      :: &
      &  hor_intp_type, lonlat_id

    ! set default values
    hor_interp_meta%hor_intp_type    = HINTP_TYPE_LONLAT_RBF
    hor_interp_meta%lonlat_id        = 0 ! invalid ID

    ! supersede with user definitions
    CALL assign_if_present(hor_interp_meta%hor_intp_type, hor_intp_type)
    CALL assign_if_present(hor_interp_meta%lonlat_id,     lonlat_id)

  END FUNCTION create_hor_interp_metadata


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
         &                     initval, isteptype, resetval, lmiss,            &
         &                     missval, tlev_source, tracer_info, vert_interp, &
         &                     hor_interp, in_group, verbose, cdiTimeID,       &
         &                     l_pp_scheduler_task, post_op)
    !
    TYPE(t_var_metadata),    INTENT(inout)        :: info          ! memory info struct.
    CHARACTER(len=*),        INTENT(in), OPTIONAL :: name          ! variable name
    INTEGER,                 INTENT(in), OPTIONAL :: hgrid         ! horizontal grid type used 
    INTEGER,                 INTENT(in), OPTIONAL :: vgrid         ! vertical grid type used
    TYPE(t_cf_var),          INTENT(in), OPTIONAL :: cf            ! CF convention
    TYPE(t_grib2_var),       INTENT(in), OPTIONAL :: grib2         ! GRIB2
    INTEGER,                 INTENT(in), OPTIONAL :: ldims(:)      ! used dimensions 
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
    TYPE(t_tracer_meta),     INTENT(in), OPTIONAL :: tracer_info   ! tracer meta data
    TYPE(t_vert_interp_meta),INTENT(in), OPTIONAL :: vert_interp   ! vertical interpolation metadata
    TYPE(t_hor_interp_meta), INTENT(in), OPTIONAL :: hor_interp    ! horizontal interpolation metadata
    LOGICAL, INTENT(in), OPTIONAL :: in_group(SIZE(VAR_GROUPS))    ! groups to which a variable belongs
    LOGICAL,                 INTENT(in), OPTIONAL :: verbose
    INTEGER,                 INTENT(in), OPTIONAL :: cdiTimeID     ! CDI time mode (TIME_VARIABLE/TIME_CONSTANT)
    INTEGER,                 INTENT(in), OPTIONAL :: l_pp_scheduler_task ! .TRUE., if field is updated by pp scheduler
    TYPE(t_post_op_meta),    INTENT(IN), OPTIONAL :: post_op       !< "post-op" (small arithmetic operations) for this variable
    !
    LOGICAL :: lverbose
    !
    ! set flags from optional parameters
    !
    lverbose = .FALSE. 
    CALL assign_if_present (lverbose, verbose)
    !
    ! set components describing the 'Content of the field'
    !
    CALL assign_if_present (info%name,  name)
    CALL struct_assign_if_present (info%cf,    cf)
    CALL struct_assign_if_present (info%grib2, grib2)
    !
    ! hash variable name for fast search
    !
    IF (PRESENT(name)) info%key = util_hashword(name, LEN_TRIM(name), 0)
    !
    CALL assign_if_present (info%used_dimensions(1:SIZE(ldims)), ldims)
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
    IF (info%lcontainer) info%ncontained = 0
    CALL struct_assign_if_present (info%resetval,resetval)
    CALL assign_if_present (info%isteptype,     isteptype)
    CALL assign_if_present (info%lmiss,         lmiss)
    CALL struct_assign_if_present (info%missval,       missval)
    CALL assign_if_present (info%lrestart,      lrestart)
    CALL assign_if_present (info%lrestart_cont, lrestart_cont)
    CALL struct_assign_if_present (info%initval,       initval)
    CALL assign_if_present (info%tlev_source,   tlev_source)
    !
    ! set flags concerning tracer fields
    !
    CALL struct_assign_if_present (info%tracer,   tracer_info)
    !
    ! set flags concerning vertical interpolation
    CALL struct_assign_if_present (info%vert_interp,   vert_interp )

    ! set flags concerning horizontal interpolation
    CALL struct_assign_if_present (info%hor_interp,    hor_interp )

    ! set meta data containing the groups to which a variable belongs
    IF (PRESENT(in_group)) THEN
      info%in_group(:) = in_group(:)
    END IF

    CALL assign_if_present (info%cdiTimeID,   cdiTimeID)

    CALL assign_if_present (info%l_pp_scheduler_task, l_pp_scheduler_task)

    CALL struct_assign_if_present (info%post_op, post_op)
    
    !
    ! printout (optional)
    !
    !LK    IF (lverbose) CALL print_var_metadata (info)
    !
  END SUBROUTINE set_var_metadata
  !------------------------------------------------------------------------------------------------
  !
  ! Create a list new entry
  !
  ! Specific routines for pointers of different rank
  !
  !================================================================================================ 
  ! REAL SECTION ----------------------------------------------------------------------------------
  !
  ! create (allocate) a new table entry
  ! optionally obtain pointer to 5d-field
  ! optionally overwrite default meta data 
  !
  SUBROUTINE add_var_list_element_r5d(this_list, name, ptr,    &
       hgrid, vgrid, cf, grib2, ldims, loutput, lcontainer,    &
       lrestart, lrestart_cont, initval_r, isteptype,          &
       resetval_r, lmiss, missval_r, tlev_source, info, p5,    &
       vert_interp, hor_interp, in_group, verbose, new_element, &
       cdiTimeID, l_pp_scheduler_task, post_op)
    !
    TYPE(t_var_list),     INTENT(inout)        :: this_list           ! list
    CHARACTER(len=*),     INTENT(in)           :: name                ! name of variable
    REAL(wp),             POINTER              :: ptr(:,:,:,:,:)      ! reference to field
    INTEGER,              INTENT(in)           :: hgrid               ! horizontal grid type used
    INTEGER,              INTENT(in)           :: vgrid               ! vertical grid type used
    TYPE(t_cf_var),       INTENT(in)           :: cf                  ! CF related metadata
    TYPE(t_grib2_var),    INTENT(in)           :: grib2               ! GRIB2 related metadata
    INTEGER,              INTENT(in), OPTIONAL :: ldims(5)            ! local dimensions
    LOGICAL,              INTENT(in), OPTIONAL :: loutput             ! output flag
    LOGICAL,              INTENT(in), OPTIONAL :: lcontainer          ! container flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart            ! restart flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart_cont       ! continue restart if var not available
    REAL(wp),             INTENT(in), OPTIONAL :: initval_r           ! value if var not available
    INTEGER,              INTENT(in), OPTIONAL :: isteptype           ! type of statistical processing
    REAL(wp),             INTENT(in), OPTIONAL :: resetval_r          ! reset value (after accumulation)
    LOGICAL,              INTENT(in), OPTIONAL :: lmiss               ! missing value flag
    REAL(wp),             INTENT(in), OPTIONAL :: missval_r           ! missing value
    INTEGER,              INTENT(in), OPTIONAL :: tlev_source         ! actual TL for TL dependent vars
    TYPE(t_var_metadata), POINTER,    OPTIONAL :: info                ! returns reference to metadata
    REAL(wp),             POINTER,    OPTIONAL :: p5(:,:,:,:,:)       ! provided pointer
    TYPE(t_vert_interp_meta),INTENT(in), OPTIONAL :: vert_interp      ! vertical interpolation metadata
    TYPE(t_hor_interp_meta), INTENT(in), OPTIONAL :: hor_interp       ! horizontal interpolation metadata
    LOGICAL, INTENT(in), OPTIONAL :: in_group(SIZE(VAR_GROUPS))       ! groups to which a variable belongs
    LOGICAL,              INTENT(in), OPTIONAL :: verbose             ! print information
    TYPE(t_list_element), POINTER, OPTIONAL  :: new_element ! pointer to new var list element
    INTEGER,              INTENT(in), OPTIONAL :: cdiTimeID           ! CDI time mode (TIME_VARIABLE/TIME_CONSTANT)
    INTEGER,              INTENT(in), OPTIONAL :: l_pp_scheduler_task ! .TRUE., if field is updated by pp scheduler
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op             !< "post-op" (small arithmetic operations) for this variable
    !
    TYPE(t_list_element), POINTER :: new_list_element
    TYPE(t_union_vals) :: missval, initval, resetval
    INTEGER :: idims(5)
    INTEGER :: istat
    LOGICAL :: referenced
    !
    ! consistency check for restart and output
    !
    IF (PRESENT(lrestart)) THEN
      IF (.NOT. this_list%p%lrestart .AND. lrestart) THEN
        CALL finish('mo_var_list:add_var_list_element_r5d',                            &
             &      'for list '//TRIM(this_list%p%name)//' restarting not enabled, '// &
             &      'but restart of '//TRIM(name)//' requested.')
      ENDIF
    ENDIF
    !
    ! add list entry
    !
    CALL append_list_element (this_list, new_list_element)
    IF (PRESENT(new_element)) new_element=>new_list_element
    new_list_element%field%info = default_var_list_metadata(this_list)
    !
    ! and set meta data
    !
    IF (PRESENT(p5)) THEN
      referenced = .TRUE.
      new_list_element%field%info%allocated = .TRUE.
    ELSE
      referenced = .FALSE.
    ENDIF
    !
    CALL assign_if_present(missval%rval, missval_r)
    CALL assign_if_present(initval%rval, initval_r)
    CALL assign_if_present(resetval%rval, resetval_r)
    CALL set_var_metadata (new_list_element%field%info,                           &
         name=name, hgrid=hgrid, vgrid=vgrid,                                     & 
         cf=cf, grib2=grib2, ldims=ldims, loutput=loutput, lcontainer=lcontainer, &
         lrestart=lrestart, lrestart_cont=lrestart_cont, initval=initval,         &
         isteptype=isteptype, resetval=resetval, lmiss=lmiss,                     &
         missval=missval, tlev_source=tlev_source, vert_interp=vert_interp,       &
         hor_interp=hor_interp, in_group=in_group, verbose=verbose,               &
         cdiTimeID=cdiTimeID, l_pp_scheduler_task=l_pp_scheduler_task,            &
         post_op=post_op)
    !
    IF (.NOT. referenced) THEN
      CALL assign_if_present(new_list_element%field%info%used_dimensions(1:5), ldims(1:5))
      idims(1) = new_list_element%field%info%used_dimensions(1)
      idims(2) = new_list_element%field%info%used_dimensions(2)
      idims(3) = new_list_element%field%info%used_dimensions(3)
      idims(4) = new_list_element%field%info%used_dimensions(4)
      idims(5) = new_list_element%field%info%used_dimensions(5)
      new_list_element%field%info%ndims = 5
      new_list_element%field%info%cdiDataType = DATATYPE_FLT64
      new_list_element%field%var_base_size = 8
      ALLOCATE(new_list_element%field%r_ptr(idims(1), idims(2), idims(3), idims(4), idims(5)), &
           &   STAT=istat)
      IF (istat /= 0) THEN
        CALL finish('mo_var_list:add_var_list_element_r5d', &
             &      'allocation of array '//TRIM(name)//' failed')
      ELSE
        new_list_element%field%info%allocated = .TRUE.
      ENDIF
      NULLIFY(new_list_element%field%i_ptr)
      NULLIFY(new_list_element%field%l_ptr)
      this_list%p%memory_used = this_list%p%memory_used &
           +INT(new_list_element%field%var_base_size*SIZE(new_list_element%field%r_ptr),i8)
    ELSE
      new_list_element%field%r_ptr => p5
    ENDIF
    ptr => new_list_element%field%r_ptr(:,:,:,:,:)
    IF(PRESENT(info)) info => new_list_element%field%info
    !
    IF (PRESENT(lmiss)) THEN
      new_list_element%field%r_ptr = new_list_element%field%info%missval%rval
    ELSE
      new_list_element%field%r_ptr = 0.0_wp
    END IF
    ! 
    IF (PRESENT(initval_r)) THEN
      new_list_element%field%r_ptr = new_list_element%field%info%initval%rval
    ENDIF

  END SUBROUTINE add_var_list_element_r5d
  !------------------------------------------------------------------------------------------------
  ! create (allocate) a new table entry
  ! optionally obtain pointer to 4d-field
  ! optionally overwrite default meta data 
  !
  SUBROUTINE add_var_list_element_r4d(this_list, name, ptr,    &
       hgrid, vgrid, cf, grib2, ldims, loutput, lcontainer,    &
       lrestart, lrestart_cont, initval_r, isteptype,          &
       resetval_r, lmiss, missval_r, tlev_source, info, p5,    &
       vert_interp, hor_interp, in_group, verbose, new_element, &
       cdiTimeID, l_pp_scheduler_task, post_op)
    !
    TYPE(t_var_list),     INTENT(inout)        :: this_list           ! list
    CHARACTER(len=*),     INTENT(in)           :: name                ! name of variable
    REAL(wp),             POINTER              :: ptr(:,:,:,:)        ! reference to field
    INTEGER,              INTENT(in)           :: hgrid               ! horizontal grid type used
    INTEGER,              INTENT(in)           :: vgrid               ! vertical grid type used
    TYPE(t_cf_var),       INTENT(in)           :: cf                  ! CF related metadata
    TYPE(t_grib2_var),    INTENT(in)           :: grib2               ! GRIB2 related metadata
    INTEGER,              INTENT(in), OPTIONAL :: ldims(4)            ! local dimensions
    LOGICAL,              INTENT(in), OPTIONAL :: loutput             ! output flag
    LOGICAL,              INTENT(in), OPTIONAL :: lcontainer          ! container flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart            ! restart flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart_cont       ! continue restart if var not available
    REAL(wp),             INTENT(in), OPTIONAL :: initval_r           ! value if var not available
    INTEGER,              INTENT(in), OPTIONAL :: isteptype           ! type of statistical processing
    REAL(wp),             INTENT(in), OPTIONAL :: resetval_r          ! reset value (after accumulation)
    LOGICAL,              INTENT(in), OPTIONAL :: lmiss               ! missing value flag
    REAL(wp),             INTENT(in), OPTIONAL :: missval_r           ! missing value
    INTEGER,              INTENT(in), OPTIONAL :: tlev_source         ! actual TL for TL dependent vars
    TYPE(t_var_metadata), POINTER,    OPTIONAL :: info                ! returns reference to metadata
    REAL(wp),             POINTER,    OPTIONAL :: p5(:,:,:,:,:)       ! provided pointer
    TYPE(t_vert_interp_meta),INTENT(in), OPTIONAL :: vert_interp      ! vertical interpolation metadata
    TYPE(t_hor_interp_meta), INTENT(in), OPTIONAL :: hor_interp       ! horizontal interpolation metadata
    LOGICAL, INTENT(in), OPTIONAL :: in_group(SIZE(VAR_GROUPS))       ! groups to which a variable belongs
    LOGICAL,              INTENT(in), OPTIONAL :: verbose             ! print information
    TYPE(t_list_element), POINTER, OPTIONAL  :: new_element ! pointer to new var list element
    INTEGER,              INTENT(in), OPTIONAL :: cdiTimeID           ! CDI time mode (TIME_VARIABLE/TIME_CONSTANT)
    INTEGER,              INTENT(in), OPTIONAL :: l_pp_scheduler_task ! .TRUE., if field is updated by pp scheduler
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op             !< "post-op" (small arithmetic operations) for this variable
    !
    TYPE(t_list_element), POINTER :: new_list_element
    TYPE(t_union_vals) :: missval, initval, resetval
    INTEGER :: idims(5)
    INTEGER :: istat
    LOGICAL :: referenced
    !
    ! consistency check for restart and output
    !
    IF (PRESENT(lrestart)) THEN
      IF (.NOT. this_list%p%lrestart .AND. lrestart) THEN
        CALL finish('mo_var_list:add_var_list_element_r4d',                            &
             &      'for list '//TRIM(this_list%p%name)//' restarting not enabled, '// &
             &      'but restart of '//TRIM(name)//' requested.')
      ENDIF
    ENDIF
    !
    ! add list entry
    !
    CALL append_list_element (this_list, new_list_element)
    IF (PRESENT(new_element)) new_element=>new_list_element
    new_list_element%field%info = default_var_list_metadata(this_list)
    !
    ! and set meta data
    !
    IF (PRESENT(p5)) THEN
      referenced = .TRUE.
      new_list_element%field%info%allocated = .TRUE.
    ELSE
      referenced = .FALSE.
    ENDIF
    !
    CALL assign_if_present(missval%rval, missval_r)
    CALL assign_if_present(initval%rval, initval_r)
    CALL assign_if_present(resetval%rval, resetval_r)
    CALL set_var_metadata (new_list_element%field%info,                           &
         name=name, hgrid=hgrid, vgrid=vgrid,                                     & 
         cf=cf, grib2=grib2, ldims=ldims, loutput=loutput, lcontainer=lcontainer, &
         lrestart=lrestart, lrestart_cont=lrestart_cont, initval=initval,         &
         isteptype=isteptype, resetval=resetval, lmiss=lmiss,                     &
         missval=missval, tlev_source=tlev_source, vert_interp=vert_interp,       &
         hor_interp=hor_interp, in_group=in_group, verbose=verbose,               &
         cdiTimeID=cdiTimeID, l_pp_scheduler_task=l_pp_scheduler_task,            &
         post_op=post_op)
    !
    IF (.NOT. referenced) THEN
      CALL assign_if_present(new_list_element%field%info%used_dimensions(1:4), ldims(1:4))
      idims(1) = new_list_element%field%info%used_dimensions(1)
      idims(2) = new_list_element%field%info%used_dimensions(2)
      idims(3) = new_list_element%field%info%used_dimensions(3)
      idims(4) = new_list_element%field%info%used_dimensions(4)
      idims(5) = 1
      new_list_element%field%info%ndims = 4
      new_list_element%field%info%cdiDataType = DATATYPE_FLT64
      new_list_element%field%var_base_size = 8
      ALLOCATE(new_list_element%field%r_ptr(idims(1), idims(2), idims(3), idims(4), idims(5)), &
           &   STAT=istat)
      IF (istat /= 0) THEN
        CALL finish('mo_var_list:add_var_list_element_r4d', &
             &      'allocation of array '//TRIM(name)//' failed')
      ELSE
        new_list_element%field%info%allocated = .TRUE.
      ENDIF
      NULLIFY(new_list_element%field%i_ptr)
      NULLIFY(new_list_element%field%l_ptr)
      this_list%p%memory_used = this_list%p%memory_used &
           +INT(new_list_element%field%var_base_size*SIZE(new_list_element%field%r_ptr),i8)
    ELSE
      new_list_element%field%r_ptr => p5
    ENDIF
    ptr => new_list_element%field%r_ptr(:,:,:,:,1)
    IF(PRESENT(info)) info => new_list_element%field%info
    !
    IF (PRESENT(lmiss)) THEN
      new_list_element%field%r_ptr = new_list_element%field%info%missval%rval
    ELSE
      new_list_element%field%r_ptr = 0.0_wp
    END IF
    !
    IF (PRESENT(initval_r)) THEN
      new_list_element%field%r_ptr = new_list_element%field%info%initval%rval
    ENDIF
    !
  END SUBROUTINE add_var_list_element_r4d
  !------------------------------------------------------------------------------------------------
  !
  ! create (allocate) a new table entry
  ! optionally obtain pointer to 3d-field
  ! optionally overwrite default meta data 
  !
  SUBROUTINE add_var_list_element_r3d(this_list, name, ptr,    &
       hgrid, vgrid, cf, grib2, ldims, loutput, lcontainer,    &
       lrestart, lrestart_cont, initval_r, isteptype,          &
       resetval_r, lmiss, missval_r, tlev_source, info, p5,    &
       vert_interp, hor_interp, in_group, verbose, new_element, &
       cdiTimeID, l_pp_scheduler_task, post_op)
    !
    TYPE(t_var_list),     INTENT(inout)        :: this_list           ! list
    CHARACTER(len=*),     INTENT(in)           :: name                ! name of variable
    REAL(wp),             POINTER              :: ptr(:,:,:)          ! reference to field
    INTEGER,              INTENT(in)           :: hgrid               ! horizontal grid type used
    INTEGER,              INTENT(in)           :: vgrid               ! vertical grid type used
    TYPE(t_cf_var),       INTENT(in)           :: cf                  ! CF related metadata
    TYPE(t_grib2_var),    INTENT(in)           :: grib2               ! GRIB2 related metadata
    INTEGER,              INTENT(in), OPTIONAL :: ldims(3)            ! local dimensions
    LOGICAL,              INTENT(in), OPTIONAL :: loutput             ! output flag
    LOGICAL,              INTENT(in), OPTIONAL :: lcontainer          ! container flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart            ! restart flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart_cont       ! continue restart if var not available
    REAL(wp),             INTENT(in), OPTIONAL :: initval_r           ! value if var not available
    INTEGER,              INTENT(in), OPTIONAL :: isteptype           ! type of statistical processing
    REAL(wp),             INTENT(in), OPTIONAL :: resetval_r          ! reset value (after accumulation)
    LOGICAL,              INTENT(in), OPTIONAL :: lmiss               ! missing value flag
    REAL(wp),             INTENT(in), OPTIONAL :: missval_r           ! missing value
    INTEGER,              INTENT(in), OPTIONAL :: tlev_source         ! actual TL for TL dependent vars
    TYPE(t_var_metadata), POINTER,    OPTIONAL :: info                ! returns reference to metadata
    REAL(wp),             POINTER,    OPTIONAL :: p5(:,:,:,:,:)       ! provided pointer
    TYPE(t_vert_interp_meta),INTENT(in), OPTIONAL :: vert_interp      ! vertical interpolation metadata
    TYPE(t_hor_interp_meta), INTENT(in), OPTIONAL :: hor_interp       ! horizontal interpolation metadata
    LOGICAL, INTENT(in), OPTIONAL :: in_group(SIZE(VAR_GROUPS))       ! groups to which a variable belongs
    LOGICAL,              INTENT(in), OPTIONAL :: verbose             ! print information
    TYPE(t_list_element), POINTER, OPTIONAL  :: new_element ! pointer to new var list element
    INTEGER,              INTENT(in), OPTIONAL :: cdiTimeID           ! CDI time mode (TIME_VARIABLE/TIME_CONSTANT)
    INTEGER,              INTENT(in), OPTIONAL :: l_pp_scheduler_task ! .TRUE., if field is updated by pp scheduler
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op             !< "post-op" (small arithmetic operations) for this variable
    !
    TYPE(t_list_element), POINTER :: new_list_element
    TYPE(t_union_vals) :: missval, initval, resetval
    INTEGER :: idims(5)
    INTEGER :: istat
    LOGICAL :: referenced
    !
    ! consistency check for restart and output
    !
    IF (PRESENT(lrestart)) THEN
      IF (.NOT. this_list%p%lrestart .AND. lrestart) THEN
        CALL finish('mo_var_list:add_var_list_element_r3d',                            &
             &      'for list '//TRIM(this_list%p%name)//' restarting not enabled, '// &
             &      'but restart of '//TRIM(name)//' requested.')
      ENDIF
    ENDIF
    !
    ! add list entry
    !
    CALL append_list_element (this_list, new_list_element)
    IF (PRESENT(new_element)) new_element=>new_list_element
    new_list_element%field%info = default_var_list_metadata(this_list)
    !
    ! and set meta data
    !
    IF (PRESENT(p5)) THEN
      referenced = .TRUE.
      new_list_element%field%info%allocated = .TRUE.
    ELSE
      referenced = .FALSE.
    ENDIF
    !
    CALL assign_if_present(missval%rval, missval_r)
    CALL assign_if_present(initval%rval, initval_r)
    CALL assign_if_present(resetval%rval, resetval_r)
    CALL set_var_metadata (new_list_element%field%info,                           &
         name=name, hgrid=hgrid, vgrid=vgrid,                                     & 
         cf=cf, grib2=grib2, ldims=ldims, loutput=loutput, lcontainer=lcontainer, &
         lrestart=lrestart, lrestart_cont=lrestart_cont, initval=initval,         &
         isteptype=isteptype, resetval=resetval, lmiss=lmiss,                     &
         missval=missval, tlev_source=tlev_source, vert_interp=vert_interp,       &
         hor_interp=hor_interp, in_group=in_group, verbose=verbose,               &
         cdiTimeID=cdiTimeID, l_pp_scheduler_task=l_pp_scheduler_task,            &
         post_op=post_op)
    !
    IF (.NOT. referenced) THEN
      CALL assign_if_present(new_list_element%field%info%used_dimensions(1:3), ldims(1:3))
      idims(1) = new_list_element%field%info%used_dimensions(1)
      idims(2) = new_list_element%field%info%used_dimensions(2)
      idims(3) = new_list_element%field%info%used_dimensions(3)
      idims(4) = 1
      idims(5) = 1
      new_list_element%field%info%ndims = 3
      new_list_element%field%info%cdiDataType = DATATYPE_FLT64
      new_list_element%field%var_base_size = 8
      ALLOCATE (new_list_element%field%r_ptr(idims(1), idims(2), idims(3), idims(4), idims(5)), &
           &    STAT=istat)
      IF (istat /= 0) THEN
        CALL finish('mo_var_list:add_var_list_element_r3d', &
             &      'allocation of array '//TRIM(name)//' failed')
      ELSE
        new_list_element%field%info%allocated = .TRUE.
      ENDIF
      NULLIFY(new_list_element%field%i_ptr)
      NULLIFY(new_list_element%field%l_ptr)
      this_list%p%memory_used = this_list%p%memory_used &
           +INT(new_list_element%field%var_base_size*SIZE(new_list_element%field%r_ptr),i8)
    ELSE
      new_list_element%field%r_ptr => p5
    ENDIF
    ptr => new_list_element%field%r_ptr(:,:,:,1,1)
    IF(PRESENT(info)) info => new_list_element%field%info
    !
    IF (PRESENT(lmiss)) THEN
      new_list_element%field%r_ptr = new_list_element%field%info%missval%rval
    ELSE
      new_list_element%field%r_ptr = 0.0_wp
    END IF
    !
    IF (PRESENT(initval_r)) THEN
      new_list_element%field%r_ptr = new_list_element%field%info%initval%rval
    ENDIF
    !
  END SUBROUTINE add_var_list_element_r3d
  !------------------------------------------------------------------------------------------------
  !
  ! create (allocate) a new table entry
  ! optionally obtain pointer to 2d-field
  ! optionally overwrite default meta data 
  !
  SUBROUTINE add_var_list_element_r2d(this_list, name, ptr,    &
       hgrid, vgrid, cf, grib2, ldims, loutput, lcontainer,    &
       lrestart, lrestart_cont, initval_r, isteptype,          &
       resetval_r, lmiss, missval_r, tlev_source, info, p5,    &
       vert_interp, hor_interp, in_group, verbose, new_element, &
       cdiTimeID, l_pp_scheduler_task, post_op)
    !
    TYPE(t_var_list),     INTENT(inout)        :: this_list           ! list
    CHARACTER(len=*),     INTENT(in)           :: name                ! name of variable
    REAL(wp),             POINTER              :: ptr(:,:)            ! reference to field
    INTEGER,              INTENT(in)           :: hgrid               ! horizontal grid type used
    INTEGER,              INTENT(in)           :: vgrid               ! vertical grid type used
    TYPE(t_cf_var),       INTENT(in)           :: cf                  ! CF related metadata
    TYPE(t_grib2_var),    INTENT(in)           :: grib2               ! GRIB2 related metadata
    INTEGER,              INTENT(in), OPTIONAL :: ldims(2)            ! local dimensions
    LOGICAL,              INTENT(in), OPTIONAL :: loutput             ! output flag
    LOGICAL,              INTENT(in), OPTIONAL :: lcontainer          ! container flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart            ! restart flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart_cont       ! continue restart if var not available
    REAL(wp),             INTENT(in), OPTIONAL :: initval_r           ! value if var not available
    INTEGER,              INTENT(in), OPTIONAL :: isteptype           ! type of statistical processing
    REAL(wp),             INTENT(in), OPTIONAL :: resetval_r          ! reset value (after accumulation)
    LOGICAL,              INTENT(in), OPTIONAL :: lmiss               ! missing value flag
    REAL(wp),             INTENT(in), OPTIONAL :: missval_r           ! missing value
    INTEGER,              INTENT(in), OPTIONAL :: tlev_source         ! actual TL for TL dependent vars
    TYPE(t_var_metadata), POINTER,    OPTIONAL :: info                ! returns reference to metadata
    REAL(wp),             POINTER,    OPTIONAL :: p5(:,:,:,:,:)       ! provided pointer
    TYPE(t_vert_interp_meta),INTENT(in), OPTIONAL :: vert_interp      ! vertical interpolation metadata
    TYPE(t_hor_interp_meta), INTENT(in), OPTIONAL :: hor_interp       ! horizontal interpolation metadata
    LOGICAL, INTENT(in), OPTIONAL :: in_group(SIZE(VAR_GROUPS))       ! groups to which a variable belongs
    LOGICAL,              INTENT(in), OPTIONAL :: verbose             ! print information
    TYPE(t_list_element), POINTER, OPTIONAL  :: new_element ! pointer to new var list element
    INTEGER,              INTENT(in), OPTIONAL :: cdiTimeID           ! CDI time mode (TIME_VARIABLE/TIME_CONSTANT)
    INTEGER,              INTENT(in), OPTIONAL :: l_pp_scheduler_task ! .TRUE., if field is updated by pp scheduler
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op             !< "post-op" (small arithmetic operations) for this variable
    !
    TYPE(t_list_element), POINTER :: new_list_element
    TYPE(t_union_vals) :: missval, initval, resetval
    INTEGER :: idims(5)
    INTEGER :: istat
    LOGICAL :: referenced
    !
    ! consistency check for restart and output
    !
    IF (PRESENT(lrestart)) THEN
      IF (.NOT. this_list%p%lrestart .AND. lrestart) THEN
        CALL finish('mo_var_list:add_var_list_element_r2d',                            &
             &      'for list '//TRIM(this_list%p%name)//' restarting not enabled, '// &
             &      'but restart of '//TRIM(name)//' requested.')
      ENDIF
    ENDIF
    !
    ! add list entry
    !
    CALL append_list_element (this_list, new_list_element)
    IF (PRESENT(new_element)) new_element=>new_list_element
    new_list_element%field%info = default_var_list_metadata(this_list)
    !
    ! and set meta data
    !
    IF (PRESENT(p5)) THEN
      referenced = .TRUE.
      new_list_element%field%info%allocated = .TRUE.
    ELSE
      referenced = .FALSE.
    ENDIF
    !
    CALL assign_if_present(missval%rval, missval_r)
    CALL assign_if_present(initval%rval, initval_r)
    CALL assign_if_present(resetval%rval, resetval_r)
    CALL set_var_metadata (new_list_element%field%info,                           &
         name=name, hgrid=hgrid, vgrid=vgrid,                                     & 
         cf=cf, grib2=grib2, ldims=ldims, loutput=loutput, lcontainer=lcontainer, &
         lrestart=lrestart, lrestart_cont=lrestart_cont, initval=initval,         &
         isteptype=isteptype, resetval=resetval, lmiss=lmiss,                     &
         missval=missval, tlev_source=tlev_source, vert_interp=vert_interp,       &
         hor_interp=hor_interp, in_group=in_group, verbose=verbose,               &
         cdiTimeID=cdiTimeID, l_pp_scheduler_task=l_pp_scheduler_task,            &
         post_op=post_op)
    !
    IF (.NOT. referenced) THEN
      CALL assign_if_present(new_list_element%field%info%used_dimensions(1:2), ldims(1:2))
      idims(1) = new_list_element%field%info%used_dimensions(1)
      idims(2) = new_list_element%field%info%used_dimensions(2)
      idims(3) = 1
      idims(4) = 1
      idims(5) = 1
      new_list_element%field%info%ndims = 2
      new_list_element%field%info%cdiDataType = DATATYPE_FLT64
      new_list_element%field%var_base_size = 8
      ALLOCATE (new_list_element%field%r_ptr(idims(1), idims(2), idims(3), idims(4), idims(5)), &
           &    STAT=istat)
      IF (istat /= 0) THEN
        CALL finish('mo_var_list:add_var_list_element_r2d', &
             &      'allocation of array '//TRIM(name)//' failed')
      ELSE
        new_list_element%field%info%allocated = .TRUE.
      ENDIF
      NULLIFY(new_list_element%field%i_ptr)
      NULLIFY(new_list_element%field%l_ptr)
      this_list%p%memory_used = this_list%p%memory_used &
           +INT(new_list_element%field%var_base_size*SIZE(new_list_element%field%r_ptr),i8)
    ELSE
      new_list_element%field%r_ptr => p5
    ENDIF
    ptr => new_list_element%field%r_ptr(:,:,1,1,1)
    IF(PRESENT(info)) info => new_list_element%field%info
    !
    IF (PRESENT(lmiss)) THEN
      new_list_element%field%r_ptr = new_list_element%field%info%missval%rval
    ELSE
      new_list_element%field%r_ptr = 0.0_wp
    END IF
    !
    IF (PRESENT(initval_r)) THEN
      new_list_element%field%r_ptr = new_list_element%field%info%initval%rval
    ENDIF
    !
  END SUBROUTINE add_var_list_element_r2d
  !------------------------------------------------------------------------------------------------
  !
  ! create (allocate) a new table entry
  ! optionally obtain pointer to 1d-field
  ! optionally overwrite default meta data 
  !
  SUBROUTINE add_var_list_element_r1d(this_list, name, ptr,    &
       hgrid, vgrid, cf, grib2, ldims, loutput, lcontainer,    &
       lrestart, lrestart_cont, initval_r, isteptype,          &
       resetval_r, lmiss, missval_r, tlev_source, info, p5,    &
       vert_interp, hor_interp, in_group, verbose, new_element, &
       cdiTimeID, l_pp_scheduler_task, post_op)
    !
    TYPE(t_var_list),     INTENT(inout)        :: this_list           ! list
    CHARACTER(len=*),     INTENT(in)           :: name                ! name of variable
    REAL(wp),             POINTER              :: ptr(:)              ! reference to field
    INTEGER,              INTENT(in)           :: hgrid               ! horizontal grid type used
    INTEGER,              INTENT(in)           :: vgrid               ! vertical grid type used
    TYPE(t_cf_var),       INTENT(in)           :: cf                  ! CF related metadata
    TYPE(t_grib2_var),    INTENT(in)           :: grib2               ! GRIB2 related metadata
    INTEGER,              INTENT(in), OPTIONAL :: ldims(1)            ! local dimensions
    LOGICAL,              INTENT(in), OPTIONAL :: loutput             ! output flag
    LOGICAL,              INTENT(in), OPTIONAL :: lcontainer          ! container flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart            ! restart flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart_cont       ! continue restart if var not available
    REAL(wp),             INTENT(in), OPTIONAL :: initval_r           ! value if var not available
    INTEGER,              INTENT(in), OPTIONAL :: isteptype           ! type of statistical processing
    REAL(wp),             INTENT(in), OPTIONAL :: resetval_r          ! reset value (after accumulation)
    LOGICAL,              INTENT(in), OPTIONAL :: lmiss               ! missing value flag
    REAL(wp),             INTENT(in), OPTIONAL :: missval_r           ! missing value
    INTEGER,              INTENT(in), OPTIONAL :: tlev_source         ! actual TL for TL dependent vars
    TYPE(t_var_metadata), POINTER,    OPTIONAL :: info                ! returns reference to metadata
    REAL(wp),             POINTER,    OPTIONAL :: p5(:,:,:,:,:)       ! provided pointer
    TYPE(t_vert_interp_meta),INTENT(in), OPTIONAL :: vert_interp      ! vertical interpolation metadata
    TYPE(t_hor_interp_meta), INTENT(in), OPTIONAL :: hor_interp       ! horizontal interpolation metadata
    LOGICAL, INTENT(in), OPTIONAL :: in_group(SIZE(VAR_GROUPS))       ! groups to which a variable belongs
    LOGICAL,              INTENT(in), OPTIONAL :: verbose             ! print information
    TYPE(t_list_element), POINTER, OPTIONAL  :: new_element ! pointer to new var list element
    INTEGER,              INTENT(in), OPTIONAL :: cdiTimeID           ! CDI time mode (TIME_VARIABLE/TIME_CONSTANT)
    INTEGER,              INTENT(in), OPTIONAL :: l_pp_scheduler_task ! .TRUE., if field is updated by pp scheduler
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op             !< "post-op" (small arithmetic operations) for this variable
    !
    TYPE(t_list_element), POINTER :: new_list_element
    TYPE(t_union_vals) :: missval, initval, resetval
    INTEGER :: idims(5)
    INTEGER :: istat
    LOGICAL :: referenced
    !
    ! consistency check for restart and output
    !
    IF (PRESENT(lrestart)) THEN
      IF (.NOT. this_list%p%lrestart .AND. lrestart) THEN
        CALL finish('mo_var_list:add_var_list_element_r1d',                            &
             &      'for list '//TRIM(this_list%p%name)//' restarting not enabled, '// &
             &      'but restart of '//TRIM(name)//' requested.')
      ENDIF
    ENDIF
    !
    ! add list entry
    !
    CALL append_list_element (this_list, new_list_element)
    IF (PRESENT(new_element)) new_element=>new_list_element
    new_list_element%field%info = default_var_list_metadata(this_list)
    !
    ! and set meta data
    !
    IF (PRESENT(p5)) THEN
      referenced = .TRUE.
      new_list_element%field%info%allocated = .TRUE.
    ELSE
      referenced = .FALSE.
    ENDIF
    !
    CALL assign_if_present(missval%rval, missval_r)
    CALL assign_if_present(initval%rval, initval_r)
    CALL assign_if_present(resetval%rval, resetval_r)
    CALL set_var_metadata (new_list_element%field%info,                           &
         name=name, hgrid=hgrid, vgrid=vgrid,                                     & 
         cf=cf, grib2=grib2, ldims=ldims, loutput=loutput, lcontainer=lcontainer, &
         lrestart=lrestart, lrestart_cont=lrestart_cont, initval=initval,         &
         isteptype=isteptype,resetval=resetval, lmiss=lmiss,                      &
         missval=missval, tlev_source=tlev_source, vert_interp=vert_interp,       &
         hor_interp=hor_interp, in_group=in_group, verbose=verbose,               &
         cdiTimeID=cdiTimeID, l_pp_scheduler_task=l_pp_scheduler_task,            &
         post_op=post_op)
    !
    IF (.NOT. referenced) THEN
      CALL assign_if_present(new_list_element%field%info%used_dimensions(1:1), ldims(1:1))
      idims(1) = new_list_element%field%info%used_dimensions(1)
      idims(2) = 1
      idims(3) = 1
      idims(4) = 1
      idims(5) = 1
      new_list_element%field%info%ndims = 1
      new_list_element%field%info%cdiDataType = DATATYPE_FLT64
      new_list_element%field%var_base_size = 8
      ALLOCATE (new_list_element%field%r_ptr(idims(1), idims(2), idims(3), idims(4), idims(5)), &
           &    STAT=istat)
      IF (istat /= 0) THEN
        CALL finish('mo_var_list:add_var_list_element_r1d', &
             &      'allocation of array '//TRIM(name)//' failed')
      ELSE
        new_list_element%field%info%allocated = .TRUE.
      ENDIF
      NULLIFY(new_list_element%field%i_ptr)
      NULLIFY(new_list_element%field%l_ptr)
      this_list%p%memory_used = this_list%p%memory_used &
           +INT(new_list_element%field%var_base_size*SIZE(new_list_element%field%r_ptr),i8)
    ELSE
      new_list_element%field%r_ptr => p5
    ENDIF
    ptr => new_list_element%field%r_ptr(:,1,1,1,1)
    IF(PRESENT(info)) info => new_list_element%field%info
    !
    IF (PRESENT(lmiss)) THEN
      new_list_element%field%r_ptr = new_list_element%field%info%missval%rval
    ELSE
      new_list_element%field%r_ptr = 0.0_wp
    END IF
    !
    IF (PRESENT(initval_r)) THEN
      new_list_element%field%r_ptr = new_list_element%field%info%initval%rval
    ENDIF
    !
  END SUBROUTINE add_var_list_element_r1d
  !
  !================================================================================================
  ! INTEGER SECTION -------------------------------------------------------------------------------
  !
  ! create (allocate) a new table entry
  ! optionally obtain pointer to 4d-field
  ! optionally overwrite default meta data 
  !
  SUBROUTINE add_var_list_element_i5d(this_list, name, ptr,    &
       hgrid, vgrid, cf, grib2, ldims, loutput,                &
       lrestart, lrestart_cont, initval_i, isteptype,          &
       resetval_i, lmiss, missval_i, info, p5, hor_interp,     &
       vert_interp, in_group, verbose, new_element, cdiTimeID, &
       l_pp_scheduler_task, post_op)
    !
    TYPE(t_var_list),     INTENT(inout)        :: this_list           ! list
    CHARACTER(len=*),     INTENT(in)           :: name                ! name of variable
    INTEGER,              POINTER              :: ptr(:,:,:,:,:)      ! reference to field
    INTEGER,              INTENT(in)           :: hgrid               ! horizontal grid type used
    INTEGER,              INTENT(in)           :: vgrid               ! vertical grid type used
    TYPE(t_cf_var),       INTENT(in)           :: cf                  ! CF related metadata
    TYPE(t_grib2_var),    INTENT(in)           :: grib2               ! GRIB2 related metadata
    INTEGER,              INTENT(in), OPTIONAL :: ldims(5)            ! local dimensions
    LOGICAL,              INTENT(in), OPTIONAL :: loutput             ! output flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart            ! restart flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart_cont       ! continue restart if var not available
    INTEGER,              INTENT(in), OPTIONAL :: initval_i           ! value if var not available
    INTEGER,              INTENT(in), OPTIONAL :: isteptype           ! type of statistical processing
    INTEGER,              INTENT(in), OPTIONAL :: resetval_i          ! reset value (after accumulation)
    LOGICAL,              INTENT(in), OPTIONAL :: lmiss               ! missing value flag
    INTEGER,              INTENT(in), OPTIONAL :: missval_i           ! missing value
    TYPE(t_var_metadata), POINTER,    OPTIONAL :: info                ! returns reference to metadata
    INTEGER,              POINTER,    OPTIONAL :: p5(:,:,:,:,:)       ! provided pointer
    TYPE(t_vert_interp_meta),INTENT(in), OPTIONAL :: vert_interp      ! vertical interpolation metadata
    TYPE(t_hor_interp_meta), INTENT(in), OPTIONAL :: hor_interp       ! horizontal interpolation metadata
    LOGICAL, INTENT(in), OPTIONAL :: in_group(SIZE(VAR_GROUPS))       ! groups to which a variable belongs
    LOGICAL,              INTENT(in), OPTIONAL :: verbose             ! print information
    TYPE(t_list_element), POINTER, OPTIONAL  :: new_element ! pointer to new var list element
    INTEGER,              INTENT(in), OPTIONAL :: cdiTimeID           ! CDI time mode (TIME_VARIABLE/TIME_CONSTANT)
    INTEGER,              INTENT(in), OPTIONAL :: l_pp_scheduler_task ! .TRUE., if field is updated by pp scheduler
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op             !< "post-op" (small arithmetic operations) for this variable
    !
    TYPE(t_list_element), POINTER :: new_list_element
    TYPE(t_union_vals) :: missval, initval, resetval
    INTEGER :: idims(5)
    INTEGER :: istat
    LOGICAL :: referenced
    !
    ! consistency check for restart and output
    !
    IF (PRESENT(lrestart)) THEN
      IF (.NOT. this_list%p%lrestart .AND. lrestart) THEN
        CALL finish('mo_var_list:add_var_list_element_i5d',                            &
             &      'for list '//TRIM(this_list%p%name)//' restarting not enabled, '// &
             &      'but restart of '//TRIM(name)//' requested.')
      ENDIF
    ENDIF
    !
    ! add list entry
    !
    CALL append_list_element (this_list, new_list_element)
    IF (PRESENT(new_element)) new_element=>new_list_element
    new_list_element%field%info = default_var_list_metadata(this_list)
    !
    ! and set meta data
    !
    IF (PRESENT(p5)) THEN
      referenced = .TRUE.
      new_list_element%field%info%allocated = .TRUE.
    ELSE
      referenced = .FALSE.
    ENDIF
    !
    CALL assign_if_present(missval%ival, missval_i)
    CALL assign_if_present(initval%ival, initval_i)
    CALL assign_if_present(resetval%ival, resetval_i)
    CALL set_var_metadata (new_list_element%field%info,                   &
         name=name, hgrid=hgrid, vgrid=vgrid,                             & 
         cf=cf, grib2=grib2, ldims=ldims, loutput=loutput,                &
         lrestart=lrestart, lrestart_cont=lrestart_cont, initval=initval, &
         isteptype=isteptype, resetval=resetval, lmiss=lmiss,             & 
         missval=missval, vert_interp=vert_interp, hor_interp=hor_interp, &
         in_group=in_group, verbose=verbose, cdiTimeID=cdiTimeID,         &
         l_pp_scheduler_task=l_pp_scheduler_task,                         &
         post_op=post_op)
    !
    IF (.NOT. referenced) THEN
      CALL assign_if_present(new_list_element%field%info%used_dimensions(1:5), ldims(1:5))
      idims(1) = new_list_element%field%info%used_dimensions(1)
      idims(2) = new_list_element%field%info%used_dimensions(2)
      idims(3) = new_list_element%field%info%used_dimensions(3)
      idims(4) = new_list_element%field%info%used_dimensions(4)
      idims(5) = new_list_element%field%info%used_dimensions(5)
      new_list_element%field%info%ndims = 5
      new_list_element%field%info%cdiDataType = DATATYPE_INT32
      new_list_element%field%var_base_size = 4
      ALLOCATE (new_list_element%field%i_ptr(idims(1), idims(2), idims(3), idims(4), idims(5)), &
           &    STAT=istat)
      IF (istat /= 0) THEN
        CALL finish('mo_var_list:add_var_list_element_i5d', &
             &      'allocation of array '//TRIM(name)//' failed')
      ELSE
        new_list_element%field%info%allocated = .TRUE.
      ENDIF
      NULLIFY(new_list_element%field%r_ptr)
      NULLIFY(new_list_element%field%l_ptr)
      this_list%p%memory_used = this_list%p%memory_used &
           +INT(new_list_element%field%var_base_size*SIZE(new_list_element%field%i_ptr),i8)
    ELSE
      new_list_element%field%i_ptr => p5
    ENDIF
    ptr => new_list_element%field%i_ptr(:,:,:,:,:)
    IF(PRESENT(info)) info => new_list_element%field%info
    !
    IF (PRESENT(lmiss)) THEN
      new_list_element%field%i_ptr = new_list_element%field%info%missval%ival
    ELSE
      new_list_element%field%i_ptr = 0
    END IF
    !
    IF (PRESENT(initval_i)) THEN
      new_list_element%field%i_ptr = new_list_element%field%info%initval%ival
    ENDIF
    !
  END SUBROUTINE add_var_list_element_i5d
  !------------------------------------------------------------------------------------------------
  !
  ! create (allocate) a new table entry
  ! optionally obtain pointer to 4d-field
  ! optionally overwrite default meta data 
  !
  SUBROUTINE add_var_list_element_i4d(this_list, name, ptr,    &
       hgrid, vgrid, cf, grib2, ldims, loutput,                &
       lrestart, lrestart_cont, initval_i, isteptype,          &
       resetval_i, lmiss, missval_i, info, p5, vert_interp,    &
       hor_interp, in_group, verbose, new_element, cdiTimeID,  &
       l_pp_scheduler_task, post_op)
    !
    TYPE(t_var_list),     INTENT(inout)        :: this_list           ! list
    CHARACTER(len=*),     INTENT(in)           :: name                ! name of variable
    INTEGER,              POINTER              :: ptr(:,:,:,:)        ! reference to field
    INTEGER,              INTENT(in)           :: hgrid               ! horizontal grid type used
    INTEGER,              INTENT(in)           :: vgrid               ! vertical grid type used
    TYPE(t_cf_var),       INTENT(in)           :: cf                  ! CF related metadata
    TYPE(t_grib2_var),    INTENT(in)           :: grib2               ! GRIB2 related metadata
    INTEGER,              INTENT(in), OPTIONAL :: ldims(4)            ! local dimensions
    LOGICAL,              INTENT(in), OPTIONAL :: loutput             ! output flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart            ! restart flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart_cont       ! continue restart if var not available
    INTEGER,              INTENT(in), OPTIONAL :: initval_i           ! value if var not available
    INTEGER,              INTENT(in), OPTIONAL :: isteptype           ! type of statistical processing
    INTEGER,              INTENT(in), OPTIONAL :: resetval_i          ! reset value (after accumulation)
    LOGICAL,              INTENT(in), OPTIONAL :: lmiss               ! missing value flag
    INTEGER,              INTENT(in), OPTIONAL :: missval_i           ! missing value
    TYPE(t_var_metadata), POINTER,    OPTIONAL :: info                ! returns reference to metadata
    INTEGER,              POINTER,    OPTIONAL :: p5(:,:,:,:,:)       ! provided pointer
    TYPE(t_vert_interp_meta),INTENT(in), OPTIONAL :: vert_interp      ! vertical interpolation metadata
    TYPE(t_hor_interp_meta), INTENT(in), OPTIONAL :: hor_interp       ! horizontal interpolation metadata
    LOGICAL, INTENT(in), OPTIONAL :: in_group(SIZE(VAR_GROUPS))       ! groups to which a variable belongs
    LOGICAL,              INTENT(in), OPTIONAL :: verbose             ! print information
    TYPE(t_list_element), POINTER, OPTIONAL  :: new_element ! pointer to new var list element
    INTEGER,              INTENT(in), OPTIONAL :: cdiTimeID           ! CDI time mode (TIME_VARIABLE/TIME_CONSTANT)
    INTEGER,              INTENT(in), OPTIONAL :: l_pp_scheduler_task ! .TRUE., if field is updated by pp scheduler
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op             !< "post-op" (small arithmetic operations) for this variable
    !
    TYPE(t_list_element), POINTER :: new_list_element
    TYPE(t_union_vals) :: missval, initval, resetval
    INTEGER :: idims(5)
    INTEGER :: istat
    LOGICAL :: referenced
    !
    ! consistency check for restart and output
    !
    IF (PRESENT(lrestart)) THEN
      IF (.NOT. this_list%p%lrestart .AND. lrestart) THEN
        CALL finish('mo_var_list:add_var_list_element_i4d',                            &
             &      'for list '//TRIM(this_list%p%name)//' restarting not enabled, '// &
             &      'but restart of '//TRIM(name)//' requested.')
      ENDIF
    ENDIF
    !
    ! add list entry
    !
    CALL append_list_element (this_list, new_list_element)
    IF (PRESENT(new_element)) new_element=>new_list_element
    new_list_element%field%info = default_var_list_metadata(this_list)
    !
    ! and set meta data
    !
    IF (PRESENT(p5)) THEN
      referenced = .TRUE.
      new_list_element%field%info%allocated = .TRUE.
    ELSE
      referenced = .FALSE.
    ENDIF
    !
    CALL assign_if_present(missval%ival, missval_i)
    CALL assign_if_present(initval%ival, initval_i)
    CALL assign_if_present(resetval%ival, resetval_i)
    CALL set_var_metadata (new_list_element%field%info,                   &
         name=name, hgrid=hgrid, vgrid=vgrid,                             & 
         cf=cf, grib2=grib2, ldims=ldims, loutput=loutput,                &
         lrestart=lrestart, lrestart_cont=lrestart_cont, initval=initval, &
         isteptype=isteptype, resetval=resetval, lmiss=lmiss,             &
         missval=missval, vert_interp=vert_interp, hor_interp=hor_interp, &
         in_group=in_group, verbose=verbose, cdiTimeID=cdiTimeID,         &
         l_pp_scheduler_task=l_pp_scheduler_task,                         &
         post_op=post_op)
    !
    IF (.NOT. referenced) THEN
      CALL assign_if_present(new_list_element%field%info%used_dimensions(1:4), ldims(1:4))
      idims(1) = new_list_element%field%info%used_dimensions(1)
      idims(2) = new_list_element%field%info%used_dimensions(2)
      idims(3) = new_list_element%field%info%used_dimensions(3)
      idims(4) = new_list_element%field%info%used_dimensions(4)
      idims(5) = 1
      new_list_element%field%info%ndims = 4
      new_list_element%field%info%cdiDataType = DATATYPE_INT32
      new_list_element%field%var_base_size = 4
      ALLOCATE (new_list_element%field%i_ptr(idims(1), idims(2), idims(3), idims(4), idims(5)), &
           &    STAT=istat)
      IF (istat /= 0) THEN
        CALL finish('mo_var_list:add_var_list_element_i4d', &
             &      'allocation of array '//TRIM(name)//' failed')
      ELSE
        new_list_element%field%info%allocated = .TRUE.
      ENDIF
      NULLIFY(new_list_element%field%r_ptr)
      NULLIFY(new_list_element%field%l_ptr)
      this_list%p%memory_used = this_list%p%memory_used &
           +INT(new_list_element%field%var_base_size*SIZE(new_list_element%field%i_ptr),i8)
    ELSE
      new_list_element%field%i_ptr => p5
    ENDIF
    ptr => new_list_element%field%i_ptr(:,:,:,:,1)
    IF(PRESENT(info)) info => new_list_element%field%info
    !
    IF (PRESENT(lmiss)) THEN
      new_list_element%field%i_ptr = new_list_element%field%info%missval%ival
    ELSE
      new_list_element%field%i_ptr = 0
    END IF
    !
    IF (PRESENT(initval_i)) THEN
      new_list_element%field%i_ptr = new_list_element%field%info%initval%ival
    ENDIF
    !
  END SUBROUTINE add_var_list_element_i4d
  !------------------------------------------------------------------------------------------------
  !
  ! create (allocate) a new table entry
  ! optionally obtain pointer to 3d-field
  ! optionally overwrite default meta data 
  !
  SUBROUTINE add_var_list_element_i3d(this_list, name, ptr,    &
       hgrid, vgrid, cf, grib2, ldims, loutput,                &
       lrestart, lrestart_cont, initval_i, isteptype,          &
       resetval_i, lmiss, missval_i, info, p5, vert_interp,    &
       hor_interp, in_group, verbose, new_element, cdiTimeID,  &
       l_pp_scheduler_task, post_op)
    !
    TYPE(t_var_list),     INTENT(inout)        :: this_list           ! list
    CHARACTER(len=*),     INTENT(in)           :: name                ! name of variable
    INTEGER,              POINTER              :: ptr(:,:,:)          ! reference to field
    INTEGER,              INTENT(in)           :: hgrid               ! horizontal grid type used
    INTEGER,              INTENT(in)           :: vgrid               ! vertical grid type used
    TYPE(t_cf_var),       INTENT(in)           :: cf                  ! CF related metadata
    TYPE(t_grib2_var),    INTENT(in)           :: grib2               ! GRIB2 related metadata
    INTEGER,              INTENT(in), OPTIONAL :: ldims(3)            ! local dimensions
    LOGICAL,              INTENT(in), OPTIONAL :: loutput             ! output flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart            ! restart flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart_cont       ! continue restart if var not available
    INTEGER,              INTENT(in), OPTIONAL :: initval_i           ! value if var not available
    INTEGER,              INTENT(in), OPTIONAL :: isteptype           ! type of statistical processing
    INTEGER,              INTENT(in), OPTIONAL :: resetval_i          ! reset value (after accumulation)
    LOGICAL,              INTENT(in), OPTIONAL :: lmiss               ! missing value flag
    INTEGER,              INTENT(in), OPTIONAL :: missval_i           ! missing value
    TYPE(t_var_metadata), POINTER,    OPTIONAL :: info                ! returns reference to metadata
    INTEGER,              POINTER,    OPTIONAL :: p5(:,:,:,:,:)       ! provided pointer
    TYPE(t_vert_interp_meta),INTENT(in), OPTIONAL :: vert_interp      ! vertical interpolation metadata
    TYPE(t_hor_interp_meta), INTENT(in), OPTIONAL :: hor_interp       ! horizontal interpolation metadata
    LOGICAL, INTENT(in), OPTIONAL :: in_group(SIZE(VAR_GROUPS))       ! groups to which a variable belongs
    LOGICAL,              INTENT(in), OPTIONAL :: verbose             ! print information
    TYPE(t_list_element), POINTER, OPTIONAL  :: new_element ! pointer to new var list element
    INTEGER,              INTENT(in), OPTIONAL :: cdiTimeID           ! CDI time mode (TIME_VARIABLE/TIME_CONSTANT)
    INTEGER,              INTENT(in), OPTIONAL :: l_pp_scheduler_task ! .TRUE., if field is updated by pp scheduler
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op             !< "post-op" (small arithmetic operations) for this variable
    !
    TYPE(t_list_element), POINTER :: new_list_element
    TYPE(t_union_vals) :: missval, initval, resetval
    INTEGER :: idims(5)
    INTEGER :: istat
    LOGICAL :: referenced
    !
    ! consistency check for restart and output
    !
    IF (PRESENT(lrestart)) THEN
      IF (.NOT. this_list%p%lrestart .AND. lrestart) THEN
        CALL finish('mo_var_list:add_var_list_element_i3d',                            &
             &      'for list '//TRIM(this_list%p%name)//' restarting not enabled, '// &
             &      'but restart of '//TRIM(name)//' requested.')
      ENDIF
    ENDIF
    !
    ! add list entry
    !
    CALL append_list_element (this_list, new_list_element)
    IF (PRESENT(new_element)) new_element=>new_list_element
    new_list_element%field%info = default_var_list_metadata(this_list)
    !
    ! and set meta data
    !
    IF (PRESENT(p5)) THEN
      referenced = .TRUE.
      new_list_element%field%info%allocated = .TRUE.
    ELSE
      referenced = .FALSE.
    ENDIF
    !
    CALL assign_if_present(missval%ival, missval_i)
    CALL assign_if_present(initval%ival, initval_i)
    CALL assign_if_present(resetval%ival, resetval_i)
    CALL set_var_metadata (new_list_element%field%info,                   &
         name=name, hgrid=hgrid, vgrid=vgrid,                             & 
         cf=cf, grib2=grib2, ldims=ldims, loutput=loutput,                &
         lrestart=lrestart, lrestart_cont=lrestart_cont, initval=initval, &
         isteptype=isteptype, resetval=resetval, lmiss=lmiss,             &
         missval=missval, vert_interp=vert_interp, hor_interp=hor_interp, &
         in_group=in_group, verbose=verbose, cdiTimeID=cdiTimeID,         &
         l_pp_scheduler_task=l_pp_scheduler_task,                         &
         post_op=post_op)
    !
    IF (.NOT. referenced) THEN
      CALL assign_if_present(new_list_element%field%info%used_dimensions(1:3), ldims(1:3))
      idims(1) = new_list_element%field%info%used_dimensions(1)
      idims(2) = new_list_element%field%info%used_dimensions(2)
      idims(3) = new_list_element%field%info%used_dimensions(3)
      idims(4) = 1
      idims(5) = 1
      new_list_element%field%info%ndims = 3
      new_list_element%field%info%cdiDataType = DATATYPE_INT32
      new_list_element%field%var_base_size = 4
      ALLOCATE (new_list_element%field%i_ptr(idims(1), idims(2), idims(3), idims(4), idims(5)), &
           &    STAT=istat)
      IF (istat /= 0) THEN
        CALL finish('mo_var_list:add_var_list_element_i3d', &
             &      'allocation of array '//TRIM(name)//' failed')
      ELSE
        new_list_element%field%info%allocated = .TRUE.
      ENDIF
      NULLIFY(new_list_element%field%r_ptr)
      NULLIFY(new_list_element%field%l_ptr)
      this_list%p%memory_used = this_list%p%memory_used &
           +INT(new_list_element%field%var_base_size*SIZE(new_list_element%field%i_ptr),i8)
    ELSE
      new_list_element%field%i_ptr => p5
    ENDIF
    ptr => new_list_element%field%i_ptr(:,:,:,1,1)
    IF(PRESENT(info)) info => new_list_element%field%info
    !
    IF (PRESENT(lmiss)) THEN
      new_list_element%field%i_ptr = new_list_element%field%info%missval%ival
    ELSE
      new_list_element%field%i_ptr = 0
    END IF
    !
    IF (PRESENT(initval_i)) THEN
      new_list_element%field%i_ptr = new_list_element%field%info%initval%ival
    ENDIF
    !
  END SUBROUTINE add_var_list_element_i3d
  !------------------------------------------------------------------------------------------------
  !
  ! create (allocate) a new table entry
  ! optionally obtain pointer to 2d-field
  ! optionally overwrite default meta data 
  !
  SUBROUTINE add_var_list_element_i2d(this_list, name, ptr,    &
       hgrid, vgrid, cf, grib2, ldims, loutput,                &
       lrestart, lrestart_cont, initval_i, isteptype,          &
       resetval_i, lmiss, missval_i, info, p5, vert_interp,    &
       hor_interp, in_group, verbose, new_element, cdiTimeID,  &
       l_pp_scheduler_task, post_op)
    !
    TYPE(t_var_list),     INTENT(inout)        :: this_list           ! list
    CHARACTER(len=*),     INTENT(in)           :: name                ! name of variable
    INTEGER,              POINTER              :: ptr(:,:)            ! reference to field
    INTEGER,              INTENT(in)           :: hgrid               ! horizontal grid type used
    INTEGER,              INTENT(in)           :: vgrid               ! vertical grid type used
    TYPE(t_cf_var),       INTENT(in)           :: cf                  ! CF related metadata
    TYPE(t_grib2_var),    INTENT(in)           :: grib2               ! GRIB2 related metadata
    INTEGER,              INTENT(in), OPTIONAL :: ldims(2)            ! local dimensions
    LOGICAL,              INTENT(in), OPTIONAL :: loutput             ! output flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart            ! restart flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart_cont       ! continue restart if var not available
    INTEGER,              INTENT(in), OPTIONAL :: initval_i           ! value if var not available
    INTEGER,              INTENT(in), OPTIONAL :: isteptype           ! type of statistical processing
    INTEGER,              INTENT(in), OPTIONAL :: resetval_i          ! reset value (after accumulation)
    LOGICAL,              INTENT(in), OPTIONAL :: lmiss               ! missing value flag
    INTEGER,              INTENT(in), OPTIONAL :: missval_i           ! missing value
    TYPE(t_var_metadata), POINTER,    OPTIONAL :: info                ! returns reference to metadata
    INTEGER,              POINTER,    OPTIONAL :: p5(:,:,:,:,:)       ! provided pointer
    TYPE(t_vert_interp_meta),INTENT(in), OPTIONAL :: vert_interp      ! vertical interpolation metadata
    TYPE(t_hor_interp_meta), INTENT(in), OPTIONAL :: hor_interp       ! horizontal interpolation metadata
    LOGICAL, INTENT(in), OPTIONAL :: in_group(SIZE(VAR_GROUPS))       ! groups to which a variable belongs
    LOGICAL,              INTENT(in), OPTIONAL :: verbose             ! print information
    TYPE(t_list_element), POINTER, OPTIONAL  :: new_element ! pointer to new var list element
    INTEGER,              INTENT(in), OPTIONAL :: cdiTimeID           ! CDI time mode (TIME_VARIABLE/TIME_CONSTANT)
    INTEGER,              INTENT(in), OPTIONAL :: l_pp_scheduler_task ! .TRUE., if field is updated by pp scheduler
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op             !< "post-op" (small arithmetic operations) for this variable
    !
    TYPE(t_list_element), POINTER :: new_list_element
    TYPE(t_union_vals) :: missval, initval, resetval
    INTEGER :: idims(5)
    INTEGER :: istat
    LOGICAL :: referenced
    !
    ! consistency check for restart and output
    !
    IF (PRESENT(lrestart)) THEN
      IF (.NOT. this_list%p%lrestart .AND. lrestart) THEN
        CALL finish('mo_var_list:add_var_list_element_i2d',                            &
             &      'for list '//TRIM(this_list%p%name)//' restarting not enabled, '// &
             &      'but restart of '//TRIM(name)//' requested.')
      ENDIF
    ENDIF
    !
    ! add list entry
    !
    CALL append_list_element (this_list, new_list_element)
    IF (PRESENT(new_element)) new_element=>new_list_element
    new_list_element%field%info = default_var_list_metadata(this_list)
    !
    ! and set meta data
    !
    IF (PRESENT(p5)) THEN
      referenced = .TRUE.
      new_list_element%field%info%allocated = .TRUE.
    ELSE
      referenced = .FALSE.
    ENDIF
    !
    CALL assign_if_present(missval%ival, missval_i)
    CALL assign_if_present(initval%ival, initval_i)
    CALL assign_if_present(resetval%ival, resetval_i)
    CALL set_var_metadata (new_list_element%field%info,                   &
         name=name, hgrid=hgrid, vgrid=vgrid,                             & 
         cf=cf, grib2=grib2, ldims=ldims, loutput=loutput,                &
         lrestart=lrestart, lrestart_cont=lrestart_cont, initval=initval, &
         isteptype=isteptype, resetval=resetval, lmiss=lmiss,             &
         missval=missval, vert_interp=vert_interp, hor_interp=hor_interp, &
         in_group=in_group, verbose=verbose,                              &
         l_pp_scheduler_task=l_pp_scheduler_task,                         &
         post_op=post_op)
    !
    IF (.NOT. referenced) THEN
      CALL assign_if_present(new_list_element%field%info%used_dimensions(1:2), ldims(1:2))
      idims(1) = new_list_element%field%info%used_dimensions(1)
      idims(2) = new_list_element%field%info%used_dimensions(2)
      idims(3) = 1
      idims(4) = 1
      idims(5) = 1
      new_list_element%field%info%ndims = 2
      new_list_element%field%info%cdiDataType = DATATYPE_INT32
      new_list_element%field%var_base_size = 4
      ALLOCATE (new_list_element%field%i_ptr(idims(1), idims(2), idims(3), idims(4), idims(5)), &
           &    STAT=istat)
      IF (istat /= 0) THEN
        CALL finish('mo_var_list:add_var_list_element_i2d', &
             &      'allocation of array '//TRIM(name)//' failed')
      ELSE
        new_list_element%field%info%allocated = .TRUE.
      ENDIF
      NULLIFY(new_list_element%field%r_ptr)
      NULLIFY(new_list_element%field%l_ptr)
      this_list%p%memory_used = this_list%p%memory_used &
           +INT(new_list_element%field%var_base_size*SIZE(new_list_element%field%i_ptr),i8)
    ELSE
      new_list_element%field%i_ptr => p5
    ENDIF
    ptr => new_list_element%field%i_ptr(:,:,1,1,1)
    IF(PRESENT(info)) info => new_list_element%field%info
    !
    IF (PRESENT(lmiss)) THEN
      new_list_element%field%i_ptr = new_list_element%field%info%missval%ival
    ELSE
      new_list_element%field%i_ptr = 0
    END IF
    !
    IF (PRESENT(initval_i)) THEN
      new_list_element%field%i_ptr = new_list_element%field%info%initval%ival
    ENDIF
    !
  END SUBROUTINE add_var_list_element_i2d
  !------------------------------------------------------------------------------------------------
  !
  ! create (allocate) a new table entry
  ! optionally obtain pointer to 1d-field
  ! optionally overwrite default meta data 
  !
  SUBROUTINE add_var_list_element_i1d(this_list, name, ptr,    &
       hgrid, vgrid, cf, grib2, ldims, loutput,                &
       lrestart, lrestart_cont, initval_i, isteptype,          &
       resetval_i, lmiss, missval_i, info, p5, vert_interp,    &
       hor_interp, in_group, verbose, new_element, cdiTimeID,  &
       l_pp_scheduler_task, post_op)
    !
    TYPE(t_var_list),     INTENT(inout)        :: this_list           ! list
    CHARACTER(len=*),     INTENT(in)           :: name                ! name of variable
    INTEGER,              POINTER              :: ptr(:)              ! reference to field
    INTEGER,              INTENT(in)           :: hgrid               ! horizontal grid type used
    INTEGER,              INTENT(in)           :: vgrid               ! vertical grid type used
    TYPE(t_cf_var),       INTENT(in)           :: cf                  ! CF related metadata
    TYPE(t_grib2_var),    INTENT(in)           :: grib2               ! GRIB2 related metadata
    INTEGER,              INTENT(in), OPTIONAL :: ldims(1)            ! local dimensions
    LOGICAL,              INTENT(in), OPTIONAL :: loutput             ! output flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart            ! restart flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart_cont       ! continue restart if var not available
    INTEGER,              INTENT(in), OPTIONAL :: initval_i           ! value if var not available
    INTEGER,              INTENT(in), OPTIONAL :: isteptype           ! type of statistical processing
    INTEGER,              INTENT(in), OPTIONAL :: resetval_i          ! reset value (after accumulation)
    LOGICAL,              INTENT(in), OPTIONAL :: lmiss               ! missing value flag
    INTEGER,              INTENT(in), OPTIONAL :: missval_i           ! missing value
    TYPE(t_var_metadata), POINTER,    OPTIONAL :: info                ! returns reference to metadata
    INTEGER,              POINTER,    OPTIONAL :: p5(:,:,:,:,:)       ! provided pointer
    TYPE(t_vert_interp_meta),INTENT(in), OPTIONAL :: vert_interp      ! vertical interpolation metadata
    TYPE(t_hor_interp_meta), INTENT(in), OPTIONAL :: hor_interp       ! horizontal interpolation metadata
    LOGICAL, INTENT(in), OPTIONAL :: in_group(SIZE(VAR_GROUPS))       ! groups to which a variable belongs
    LOGICAL,              INTENT(in), OPTIONAL :: verbose             ! print information
    TYPE(t_list_element), POINTER, OPTIONAL  :: new_element ! pointer to new var list element
    INTEGER,              INTENT(in), OPTIONAL :: cdiTimeID           ! CDI time mode (TIME_VARIABLE/TIME_CONSTANT)
    INTEGER,              INTENT(in), OPTIONAL :: l_pp_scheduler_task ! .TRUE., if field is updated by pp scheduler
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op             !< "post-op" (small arithmetic operations) for this variable
    !
    TYPE(t_list_element), POINTER :: new_list_element
    TYPE(t_union_vals) :: missval, initval, resetval
    INTEGER :: idims(5)
    INTEGER :: istat
    LOGICAL :: referenced
    !
    ! consistency check for restart and output
    !
    IF (PRESENT(lrestart)) THEN
      IF (.NOT. this_list%p%lrestart .AND. lrestart) THEN
        CALL finish('mo_var_list:add_var_list_element_i1d',                            &
             &      'for list '//TRIM(this_list%p%name)//' restarting not enabled, '// &
             &      'but restart of '//TRIM(name)//' requested.')
      ENDIF
    ENDIF
    !
    ! add list entry
    !
    CALL append_list_element (this_list, new_list_element)
    IF (PRESENT(new_element)) new_element=>new_list_element
    new_list_element%field%info = default_var_list_metadata(this_list)
    !
    ! and set meta data
    !
    IF (PRESENT(p5)) THEN
      referenced = .TRUE.
      new_list_element%field%info%allocated = .TRUE.
    ELSE
      referenced = .FALSE.
    ENDIF
    !
    CALL assign_if_present(missval%ival, missval_i)
    CALL assign_if_present(initval%ival, initval_i)
    CALL assign_if_present(resetval%ival, resetval_i)
    CALL set_var_metadata (new_list_element%field%info,                   &
         name=name, hgrid=hgrid, vgrid=vgrid,                             & 
         cf=cf, grib2=grib2, ldims=ldims, loutput=loutput,                &
         lrestart=lrestart, lrestart_cont=lrestart_cont, initval=initval, &
         isteptype=isteptype, resetval=resetval, lmiss=lmiss,             &
         missval=missval, vert_interp=vert_interp, hor_interp=hor_interp, &
         in_group=in_group, verbose=verbose,                              &
         l_pp_scheduler_task=l_pp_scheduler_task,                         &
         post_op=post_op)
    !
    IF (.NOT. referenced) THEN
      CALL assign_if_present(new_list_element%field%info%used_dimensions(1:1), ldims(1:1))
      idims(1) = new_list_element%field%info%used_dimensions(1)
      idims(2) = 1
      idims(3) = 1
      idims(4) = 1
      idims(5) = 1
      new_list_element%field%info%ndims = 1
      new_list_element%field%info%cdiDataType = DATATYPE_INT32
      new_list_element%field%var_base_size = 4
      ALLOCATE (new_list_element%field%i_ptr(idims(1), idims(2), idims(3), idims(4), idims(5)), &
           &    STAT=istat)
      IF (istat /= 0) THEN
        CALL finish('mo_var_list:add_var_list_element_i1d', &
             &      'allocation of array '//TRIM(name)//' failed')
      ELSE
        new_list_element%field%info%allocated = .TRUE.
      ENDIF
      NULLIFY(new_list_element%field%r_ptr)
      NULLIFY(new_list_element%field%l_ptr)
      this_list%p%memory_used = this_list%p%memory_used &
           +INT(new_list_element%field%var_base_size*SIZE(new_list_element%field%i_ptr),i8)
    ELSE
      new_list_element%field%i_ptr => p5
    ENDIF
    ptr => new_list_element%field%i_ptr(:,1,1,1,1)
    IF(PRESENT(info)) info => new_list_element%field%info
    !
    IF (PRESENT(lmiss)) THEN
      new_list_element%field%i_ptr = new_list_element%field%info%missval%ival
    ELSE
      new_list_element%field%i_ptr = 0
    END IF
    !
    IF (PRESENT(initval_i)) THEN
      new_list_element%field%i_ptr = new_list_element%field%info%initval%ival
    ENDIF
    !
  END SUBROUTINE add_var_list_element_i1d
  !
  !================================================================================================
  ! LOGICAL SECTION -------------------------------------------------------------------------------
  !
  ! create (allocate) a new table entry
  ! optionally obtain pointer to 4d-field
  ! optionally overwrite default meta data 
  !
  SUBROUTINE add_var_list_element_l5d(this_list, name, ptr,    &
       hgrid, vgrid, cf, grib2, ldims, loutput,                &
       lrestart, lrestart_cont, initval_l, isteptype,          &
       resetval_l, lmiss, missval_l, info, p5, vert_interp,    &
       hor_interp, in_group, verbose, new_element, cdiTimeID,  &
       l_pp_scheduler_task, post_op)
    !
    TYPE(t_var_list),     INTENT(inout)        :: this_list           ! list
    CHARACTER(len=*),     INTENT(in)           :: name                ! name of variable
    LOGICAL,              POINTER              :: ptr(:,:,:,:,:)      ! reference to field
    INTEGER,              INTENT(in)           :: hgrid               ! horizontal grid type used
    INTEGER,              INTENT(in)           :: vgrid               ! vertical grid type used
    TYPE(t_cf_var),       INTENT(in)           :: cf                  ! CF related metadata
    TYPE(t_grib2_var),    INTENT(in)           :: grib2               ! GRIB2 related metadata
    INTEGER,              INTENT(in), OPTIONAL :: ldims(5)            ! local dimensions
    LOGICAL,              INTENT(in), OPTIONAL :: loutput             ! output flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart            ! restart flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart_cont       ! continue restart if var not available
    LOGICAL,              INTENT(in), OPTIONAL :: initval_l           ! value if var not available
    INTEGER,              INTENT(in), OPTIONAL :: isteptype           ! type of statistical processing
    LOGICAL,              INTENT(in), OPTIONAL :: resetval_l          ! reset value (after accumulation)
    LOGICAL,              INTENT(in), OPTIONAL :: lmiss               ! missing value flag
    LOGICAL,              INTENT(in), OPTIONAL :: missval_l           ! missing value
    TYPE(t_var_metadata), POINTER,    OPTIONAL :: info                ! returns reference to metadata
    LOGICAL,              POINTER,    OPTIONAL :: p5(:,:,:,:,:)       ! provided pointer
    TYPE(t_vert_interp_meta),INTENT(in), OPTIONAL :: vert_interp      ! vertical interpolation metadata
    TYPE(t_hor_interp_meta), INTENT(in), OPTIONAL :: hor_interp       ! horizontal interpolation metadata
    LOGICAL, INTENT(in), OPTIONAL :: in_group(SIZE(VAR_GROUPS))       ! groups to which a variable belongs
    LOGICAL,              INTENT(in), OPTIONAL :: verbose             ! print information
    TYPE(t_list_element), POINTER, OPTIONAL  :: new_element ! pointer to new var list element
    INTEGER,              INTENT(in), OPTIONAL :: cdiTimeID           ! CDI time mode (TIME_VARIABLE/TIME_CONSTANT)
    INTEGER,              INTENT(in), OPTIONAL :: l_pp_scheduler_task ! .TRUE., if field is updated by pp scheduler
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op             !< "post-op" (small arithmetic operations) for this variable
    !
    TYPE(t_list_element), POINTER :: new_list_element
    TYPE(t_union_vals) :: missval, initval, resetval
    INTEGER :: idims(5)
    INTEGER :: istat
    LOGICAL :: referenced
    !
    ! consistency check for restart and output
    !
    IF (PRESENT(lrestart)) THEN
      IF (.NOT. this_list%p%lrestart .AND. lrestart) THEN
        CALL finish('mo_var_list:add_var_list_element_l5d',                            &
             &      'for list '//TRIM(this_list%p%name)//' restarting not enabled, '// &
             &      'but restart of '//TRIM(name)//' requested.')
      ENDIF
    ENDIF
    !
    ! add list entry
    !
    CALL append_list_element (this_list, new_list_element)
    IF (PRESENT(new_element)) new_element=>new_list_element
    new_list_element%field%info = default_var_list_metadata(this_list)
    !
    ! and set meta data
    !
    IF (PRESENT(p5)) THEN
      referenced = .TRUE.
      new_list_element%field%info%allocated = .TRUE.
    ELSE
      referenced = .FALSE.
    ENDIF
    !
    CALL assign_if_present(missval%lval, missval_l)
    CALL assign_if_present(initval%lval, initval_l)
    CALL assign_if_present(resetval%lval, resetval_l)
    CALL set_var_metadata (new_list_element%field%info,                   &
         name=name, hgrid=hgrid, vgrid=vgrid,                             & 
         cf=cf, grib2=grib2, ldims=ldims, loutput=loutput,                &
         lrestart=lrestart, lrestart_cont=lrestart_cont, initval=initval, &
         isteptype=isteptype, resetval=resetval, lmiss=lmiss,             &
         missval=missval, vert_interp=vert_interp, hor_interp=hor_interp, &
         in_group=in_group, verbose=verbose,                              &
         l_pp_scheduler_task=l_pp_scheduler_task,                         &
         post_op=post_op)
    !
    IF (.NOT. referenced) THEN
      CALL assign_if_present(new_list_element%field%info%used_dimensions(1:4), ldims(1:4))
      idims(1) = new_list_element%field%info%used_dimensions(1)
      idims(2) = new_list_element%field%info%used_dimensions(2)
      idims(3) = new_list_element%field%info%used_dimensions(3)
      idims(4) = new_list_element%field%info%used_dimensions(4)
      idims(5) = new_list_element%field%info%used_dimensions(5)
      new_list_element%field%info%ndims = 5
      new_list_element%field%info%cdiDataType = DATATYPE_INT8      
      new_list_element%field%var_base_size = 4
      ALLOCATE (new_list_element%field%l_ptr(idims(1), idims(2), idims(3), idims(4), idims(5)), &
           &    STAT=istat)
      IF (istat /= 0) THEN
        CALL finish('mo_var_list:add_var_list_element_l5d', &
             &      'allocation of array '//TRIM(name)//' failed')
      ELSE
        new_list_element%field%info%allocated = .TRUE.
      ENDIF
      NULLIFY(new_list_element%field%r_ptr)
      NULLIFY(new_list_element%field%i_ptr)      
      this_list%p%memory_used = this_list%p%memory_used &
           +INT(new_list_element%field%var_base_size*SIZE(new_list_element%field%l_ptr),i8)
    ELSE
      new_list_element%field%l_ptr => p5
    ENDIF
    ptr => new_list_element%field%l_ptr(:,:,:,:,:)
    IF(PRESENT(info)) info => new_list_element%field%info
    !
    IF (PRESENT(lmiss)) THEN
      new_list_element%field%l_ptr = new_list_element%field%info%missval%lval
    ELSE
      new_list_element%field%l_ptr = .FALSE.
    END IF
    !
    IF (PRESENT(initval_l)) THEN
      new_list_element%field%l_ptr = new_list_element%field%info%initval%lval
    ENDIF
    !
  END SUBROUTINE add_var_list_element_l5d
  !------------------------------------------------------------------------------------------------
  !
  ! create (allocate) a new table entry
  ! optionally obtain pointer to 4d-field
  ! optionally overwrite default meta data 
  !
  SUBROUTINE add_var_list_element_l4d(this_list, name, ptr,    &
       hgrid, vgrid, cf, grib2, ldims, loutput,                &
       lrestart, lrestart_cont, initval_l, isteptype,          &
       resetval_l, lmiss, missval_l, info, p5, vert_interp,    &
       hor_interp, in_group, verbose, new_element, cdiTimeID,  &
       l_pp_scheduler_task, post_op)
    !
    TYPE(t_var_list),     INTENT(inout)        :: this_list           ! list
    CHARACTER(len=*),     INTENT(in)           :: name                ! name of variable
    LOGICAL,              POINTER              :: ptr(:,:,:,:)        ! reference to field
    INTEGER,              INTENT(in)           :: hgrid               ! horizontal grid type used
    INTEGER,              INTENT(in)           :: vgrid               ! vertical grid type used
    TYPE(t_cf_var),       INTENT(in)           :: cf                  ! CF related metadata
    TYPE(t_grib2_var),    INTENT(in)           :: grib2               ! GRIB2 related metadata
    INTEGER,              INTENT(in), OPTIONAL :: ldims(4)            ! local dimensions
    LOGICAL,              INTENT(in), OPTIONAL :: loutput             ! output flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart            ! restart flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart_cont       ! continue restart if var not available
    LOGICAL,              INTENT(in), OPTIONAL :: initval_l           ! value if var not available
    INTEGER,              INTENT(in), OPTIONAL :: isteptype           ! type of statistical processing
    LOGICAL,              INTENT(in), OPTIONAL :: resetval_l          ! reset value (after accumulation)
    LOGICAL,              INTENT(in), OPTIONAL :: lmiss               ! missing value flag
    LOGICAL,              INTENT(in), OPTIONAL :: missval_l           ! missing value
    TYPE(t_var_metadata), POINTER,    OPTIONAL :: info                ! returns reference to metadata
    LOGICAL,              POINTER,    OPTIONAL :: p5(:,:,:,:,:)       ! provided pointer
    TYPE(t_vert_interp_meta),INTENT(in), OPTIONAL :: vert_interp      ! vertical interpolation metadata
    TYPE(t_hor_interp_meta), INTENT(in), OPTIONAL :: hor_interp       ! horizontal interpolation metadata
    LOGICAL, INTENT(in), OPTIONAL :: in_group(SIZE(VAR_GROUPS))       ! groups to which a variable belongs
    LOGICAL,              INTENT(in), OPTIONAL :: verbose             ! print information
    TYPE(t_list_element), POINTER, OPTIONAL  :: new_element ! pointer to new var list element
    INTEGER,              INTENT(in), OPTIONAL :: cdiTimeID           ! CDI time mode (TIME_VARIABLE/TIME_CONSTANT)
    INTEGER,              INTENT(in), OPTIONAL :: l_pp_scheduler_task ! .TRUE., if field is updated by pp scheduler
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op             !< "post-op" (small arithmetic operations) for this variable
    !
    TYPE(t_list_element), POINTER :: new_list_element
    TYPE(t_union_vals) :: missval, initval, resetval
    INTEGER :: idims(5)
    INTEGER :: istat
    LOGICAL :: referenced
    !
    ! consistency check for restart and output
    !
    IF (PRESENT(lrestart)) THEN
      IF (.NOT. this_list%p%lrestart .AND. lrestart) THEN
        CALL finish('mo_var_list:add_var_list_element_l4d',                            &
             &      'for list '//TRIM(this_list%p%name)//' restarting not enabled, '// &
             &      'but restart of '//TRIM(name)//' requested.')
      ENDIF
    ENDIF
    !
    ! add list entry
    !
    CALL append_list_element (this_list, new_list_element)
    IF (PRESENT(new_element)) new_element=>new_list_element
    new_list_element%field%info = default_var_list_metadata(this_list)
    !
    ! and set meta data
    !
    IF (PRESENT(p5)) THEN
      referenced = .TRUE.
      new_list_element%field%info%allocated = .TRUE.
    ELSE
      referenced = .FALSE.
    ENDIF
    !
    CALL assign_if_present(missval%lval, missval_l)
    CALL assign_if_present(initval%lval, initval_l)
    CALL assign_if_present(resetval%lval, resetval_l)
    CALL set_var_metadata (new_list_element%field%info,                   &
         name=name, hgrid=hgrid, vgrid=vgrid,                             & 
         cf=cf, grib2=grib2, ldims=ldims, loutput=loutput,                &
         lrestart=lrestart, lrestart_cont=lrestart_cont, initval=initval, &
         isteptype=isteptype, resetval=resetval, lmiss=lmiss,             &
         missval=missval, vert_interp=vert_interp, hor_interp=hor_interp, &
         in_group=in_group, verbose=verbose,                              &
         l_pp_scheduler_task=l_pp_scheduler_task,                         &
         post_op=post_op)
    !
    IF (.NOT. referenced) THEN
      CALL assign_if_present(new_list_element%field%info%used_dimensions(1:4), ldims(1:4))
      idims(1) = new_list_element%field%info%used_dimensions(1)
      idims(2) = new_list_element%field%info%used_dimensions(2)
      idims(3) = new_list_element%field%info%used_dimensions(3)
      idims(4) = new_list_element%field%info%used_dimensions(4)
      idims(5) = 1
      new_list_element%field%info%ndims = 4
      new_list_element%field%info%cdiDataType = DATATYPE_INT8      
      new_list_element%field%var_base_size = 4
      ALLOCATE (new_list_element%field%l_ptr(idims(1), idims(2), idims(3), idims(4), idims(5)), &
           &    STAT=istat)
      IF (istat /= 0) THEN
        CALL finish('mo_var_list:add_var_list_element_l4d', &
             &      'allocation of array '//TRIM(name)//' failed')
      ELSE
        new_list_element%field%info%allocated = .TRUE.
      ENDIF
      NULLIFY(new_list_element%field%r_ptr)
      NULLIFY(new_list_element%field%i_ptr)      
      this_list%p%memory_used = this_list%p%memory_used &
           +INT(new_list_element%field%var_base_size*SIZE(new_list_element%field%l_ptr),i8)
    ELSE
      new_list_element%field%l_ptr => p5
    ENDIF
    ptr => new_list_element%field%l_ptr(:,:,:,:,1)
    IF(PRESENT(info)) info => new_list_element%field%info
    !
    IF (PRESENT(lmiss)) THEN
      new_list_element%field%l_ptr = new_list_element%field%info%missval%lval
    ELSE
      new_list_element%field%l_ptr = .FALSE.
    END IF
    !
    IF (PRESENT(initval_l)) THEN
      new_list_element%field%l_ptr = new_list_element%field%info%initval%lval
    ENDIF
    !
  END SUBROUTINE add_var_list_element_l4d
  !------------------------------------------------------------------------------------------------
  !
  ! create (allocate) a new table entry
  ! optionally obtain pointer to 3d-field
  ! optionally overwrite default meta data 
  !
  SUBROUTINE add_var_list_element_l3d(this_list, name, ptr,    &
       hgrid, vgrid, cf, grib2, ldims, loutput,                &
       lrestart, lrestart_cont, initval_l, isteptype,          &
       resetval_l, lmiss, missval_l, info, p5, vert_interp,    &
       hor_interp, in_group, verbose, new_element, cdiTimeID,  &
       l_pp_scheduler_task, post_op)
    !
    TYPE(t_var_list),     INTENT(inout)        :: this_list           ! list
    CHARACTER(len=*),     INTENT(in)           :: name                ! name of variable
    LOGICAL,              POINTER              :: ptr(:,:,:)          ! reference to field
    INTEGER,              INTENT(in)           :: hgrid               ! horizontal grid type used
    INTEGER,              INTENT(in)           :: vgrid               ! vertical grid type used
    TYPE(t_cf_var),       INTENT(in)           :: cf                  ! CF related metadata
    TYPE(t_grib2_var),    INTENT(in)           :: grib2               ! GRIB2 related metadata
    INTEGER,              INTENT(in), OPTIONAL :: ldims(3)            ! local dimensions
    LOGICAL,              INTENT(in), OPTIONAL :: loutput             ! output flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart            ! restart flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart_cont       ! continue restart if var not available
    LOGICAL,              INTENT(in), OPTIONAL :: initval_l           ! value if var not available
    INTEGER,              INTENT(in), OPTIONAL :: isteptype           ! type of statistical processing
    LOGICAL,              INTENT(in), OPTIONAL :: resetval_l          ! reset value (after accumulation)
    LOGICAL,              INTENT(in), OPTIONAL :: lmiss               ! missing value flag
    LOGICAL,              INTENT(in), OPTIONAL :: missval_l           ! missing value
    TYPE(t_var_metadata), POINTER,    OPTIONAL :: info                ! returns reference to metadata
    LOGICAL,              POINTER,    OPTIONAL :: p5(:,:,:,:,:)       ! provided pointer
    TYPE(t_vert_interp_meta),INTENT(in), OPTIONAL :: vert_interp      ! vertical interpolation metadata
    TYPE(t_hor_interp_meta), INTENT(in), OPTIONAL :: hor_interp       ! horizontal interpolation metadata
    LOGICAL, INTENT(in), OPTIONAL :: in_group(SIZE(VAR_GROUPS))       ! groups to which a variable belongs
    LOGICAL,              INTENT(in), OPTIONAL :: verbose             ! print information
    TYPE(t_list_element), POINTER, OPTIONAL  :: new_element ! pointer to new var list element
    INTEGER,              INTENT(in), OPTIONAL :: cdiTimeID           ! CDI time mode (TIME_VARIABLE/TIME_CONSTANT)
    INTEGER,              INTENT(in), OPTIONAL :: l_pp_scheduler_task ! .TRUE., if field is updated by pp scheduler
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op             !< "post-op" (small arithmetic operations) for this variable
    !
    TYPE(t_list_element), POINTER :: new_list_element
    TYPE(t_union_vals) :: missval, initval, resetval
    INTEGER :: idims(5)
    INTEGER :: istat
    LOGICAL :: referenced
    !
    ! consistency check for restart and output
    !
    IF (PRESENT(lrestart)) THEN
      IF (.NOT. this_list%p%lrestart .AND. lrestart) THEN
        CALL finish('mo_var_list:add_var_list_element_l3d',                            &
             &      'for list '//TRIM(this_list%p%name)//' restarting not enabled, '// &
             &      'but restart of '//TRIM(name)//' requested.')
      ENDIF
    ENDIF
    !
    ! add list entry
    !
    CALL append_list_element (this_list, new_list_element)
    IF (PRESENT(new_element)) new_element=>new_list_element
    new_list_element%field%info = default_var_list_metadata(this_list)
    !
    ! and set meta data
    !
    IF (PRESENT(p5)) THEN
      referenced = .TRUE.
      new_list_element%field%info%allocated = .TRUE.
    ELSE
      referenced = .FALSE.
    ENDIF
    !
    CALL assign_if_present(missval%lval, missval_l)
    CALL assign_if_present(initval%lval, initval_l)
    CALL assign_if_present(resetval%lval, resetval_l)
    CALL set_var_metadata (new_list_element%field%info,                   &
         name=name, hgrid=hgrid, vgrid=vgrid,                             & 
         cf=cf, grib2=grib2, ldims=ldims, loutput=loutput,                &
         lrestart=lrestart, lrestart_cont=lrestart_cont, initval=initval, &
         isteptype=isteptype, resetval=resetval, lmiss=lmiss,             &
         missval=missval, vert_interp=vert_interp, hor_interp=hor_interp, &
         in_group=in_group, verbose=verbose,                              &
         l_pp_scheduler_task=l_pp_scheduler_task,                         &
         post_op=post_op)
    !
    IF (.NOT. referenced) THEN
      CALL assign_if_present(new_list_element%field%info%used_dimensions(1:3), ldims(1:3))
      idims(1) = new_list_element%field%info%used_dimensions(1)
      idims(2) = new_list_element%field%info%used_dimensions(2)
      idims(3) = new_list_element%field%info%used_dimensions(3)
      idims(4) = 1
      idims(5) = 1
      new_list_element%field%info%ndims = 3
      new_list_element%field%info%cdiDataType = DATATYPE_INT8      
      new_list_element%field%var_base_size = 4
      ALLOCATE (new_list_element%field%l_ptr(idims(1), idims(2), idims(3), idims(4), idims(5)), &
           &    STAT=istat)
      IF (istat /= 0) THEN
        CALL finish('mo_var_list:add_var_list_element_l3d', &
             &      'allocation of array '//TRIM(name)//' failed')
      ELSE
        new_list_element%field%info%allocated = .TRUE.
      ENDIF
      NULLIFY(new_list_element%field%r_ptr)
      NULLIFY(new_list_element%field%i_ptr)      
      this_list%p%memory_used = this_list%p%memory_used &
           +INT(new_list_element%field%var_base_size*SIZE(new_list_element%field%l_ptr),i8)
    ELSE
      new_list_element%field%l_ptr => p5
    ENDIF
    ptr => new_list_element%field%l_ptr(:,:,:,1,1)
    IF(PRESENT(info)) info => new_list_element%field%info
    !
    IF (PRESENT(lmiss)) THEN
      new_list_element%field%l_ptr = new_list_element%field%info%missval%lval
    ELSE
      new_list_element%field%l_ptr = .FALSE.
    END IF
    !
    IF (PRESENT(initval_l)) THEN
      new_list_element%field%l_ptr = new_list_element%field%info%initval%lval
    ENDIF
    !
  END SUBROUTINE add_var_list_element_l3d
  !------------------------------------------------------------------------------------------------
  !
  ! create (allocate) a new table entry
  ! optionally obtain pointer to 2d-field
  ! optionally overwrite default meta data 
  !
  SUBROUTINE add_var_list_element_l2d(this_list, name, ptr,    &
       hgrid, vgrid, cf, grib2, ldims, loutput,                &
       lrestart, lrestart_cont, initval_l, isteptype,          &
       resetval_l, lmiss, missval_l, info, p5, vert_interp,    &
       hor_interp, in_group, verbose, new_element, cdiTimeID,  &
       l_pp_scheduler_task, post_op)
    !
    TYPE(t_var_list),     INTENT(inout)        :: this_list           ! list
    CHARACTER(len=*),     INTENT(in)           :: name                ! name of variable
    LOGICAL,              POINTER              :: ptr(:,:)            ! reference to field
    INTEGER,              INTENT(in)           :: hgrid               ! horizontal grid type used
    INTEGER,              INTENT(in)           :: vgrid               ! vertical grid type used
    TYPE(t_cf_var),       INTENT(in)           :: cf                  ! CF related metadata
    TYPE(t_grib2_var),    INTENT(in)           :: grib2               ! GRIB2 related metadata
    INTEGER,              INTENT(in), OPTIONAL :: ldims(2)            ! local dimensions
    LOGICAL,              INTENT(in), OPTIONAL :: loutput             ! output flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart            ! restart flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart_cont       ! continue restart if var not available
    LOGICAL,              INTENT(in), OPTIONAL :: initval_l           ! value if var not available
    INTEGER,              INTENT(in), OPTIONAL :: isteptype           ! type of statistical processing
    LOGICAL,              INTENT(in), OPTIONAL :: resetval_l          ! reset value (after accumulation)
    LOGICAL,              INTENT(in), OPTIONAL :: lmiss               ! missing value flag
    LOGICAL,              INTENT(in), OPTIONAL :: missval_l           ! missing value
    TYPE(t_var_metadata), POINTER,    OPTIONAL :: info                ! returns reference to metadata
    LOGICAL,              POINTER,    OPTIONAL :: p5(:,:,:,:,:)       ! provided pointer
    TYPE(t_vert_interp_meta),INTENT(in), OPTIONAL :: vert_interp      ! vertical interpolation metadata
    TYPE(t_hor_interp_meta), INTENT(in), OPTIONAL :: hor_interp       ! horizontal interpolation metadata
    LOGICAL, INTENT(in), OPTIONAL :: in_group(SIZE(VAR_GROUPS))       ! groups to which a variable belongs
    LOGICAL,              INTENT(in), OPTIONAL :: verbose             ! print information
    TYPE(t_list_element), POINTER, OPTIONAL  :: new_element ! pointer to new var list element
    INTEGER,              INTENT(in), OPTIONAL :: cdiTimeID           ! CDI time mode (TIME_VARIABLE/TIME_CONSTANT)
    INTEGER,              INTENT(in), OPTIONAL :: l_pp_scheduler_task ! .TRUE., if field is updated by pp scheduler
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op             !< "post-op" (small arithmetic operations) for this variable
    !
    TYPE(t_list_element), POINTER :: new_list_element
    TYPE(t_union_vals) :: missval, initval, resetval
    INTEGER :: idims(5)
    INTEGER :: istat
    LOGICAL :: referenced
    !
    ! consistency check for restart and output
    !
    IF (PRESENT(lrestart)) THEN
      IF (.NOT. this_list%p%lrestart .AND. lrestart) THEN
        CALL finish('mo_var_list:add_var_list_element_l2d',                            &
             &      'for list '//TRIM(this_list%p%name)//' restarting not enabled, '// &
             &      'but restart of '//TRIM(name)//' requested.')
      ENDIF
    ENDIF
    !
    ! add list entry
    !
    CALL append_list_element (this_list, new_list_element)
    IF (PRESENT(new_element)) new_element=>new_list_element
    new_list_element%field%info = default_var_list_metadata(this_list)
    !
    ! and set meta data
    !
    IF (PRESENT(p5)) THEN
      referenced = .TRUE.
      new_list_element%field%info%allocated = .TRUE.
    ELSE
      referenced = .FALSE.
    ENDIF
    !
    CALL assign_if_present(missval%lval, missval_l)
    CALL assign_if_present(initval%lval, initval_l)
    CALL assign_if_present(resetval%lval, resetval_l)
    CALL set_var_metadata (new_list_element%field%info,                   &
         name=name, hgrid=hgrid, vgrid=vgrid,                             & 
         cf=cf, grib2=grib2, ldims=ldims, loutput=loutput,                &
         lrestart=lrestart, lrestart_cont=lrestart_cont, initval=initval, &
         isteptype=isteptype, resetval=resetval, lmiss=lmiss,             &
         missval=missval, vert_interp=vert_interp, hor_interp=hor_interp, &
         in_group=in_group, verbose=verbose,                              &
         l_pp_scheduler_task=l_pp_scheduler_task,                         &
         post_op=post_op)
    !
    IF (.NOT. referenced) THEN
      CALL assign_if_present(new_list_element%field%info%used_dimensions(1:2), ldims(1:2))
      idims(1) = new_list_element%field%info%used_dimensions(1)
      idims(2) = new_list_element%field%info%used_dimensions(2)
      idims(3) = 1
      idims(4) = 1
      idims(5) = 1
      new_list_element%field%info%ndims = 2
      new_list_element%field%info%cdiDataType = DATATYPE_INT8      
      new_list_element%field%var_base_size = 4
      ALLOCATE (new_list_element%field%l_ptr(idims(1), idims(2), idims(3), idims(4), idims(5)), &
           &    STAT=istat)
      IF (istat /= 0) THEN
        CALL finish('mo_var_list:add_var_list_element_l2d', &
             &      'allocation of array '//TRIM(name)//' failed')
      ELSE
        new_list_element%field%info%allocated = .TRUE.
      ENDIF
      NULLIFY(new_list_element%field%r_ptr)
      NULLIFY(new_list_element%field%i_ptr)      
      this_list%p%memory_used = this_list%p%memory_used &
           +INT(new_list_element%field%var_base_size*SIZE(new_list_element%field%l_ptr),i8)
    ELSE
      new_list_element%field%l_ptr => p5
    ENDIF
    ptr => new_list_element%field%l_ptr(:,:,1,1,1)
    IF(PRESENT(info)) info => new_list_element%field%info
    !
    IF (PRESENT(lmiss)) THEN
      new_list_element%field%l_ptr = new_list_element%field%info%missval%lval
    ELSE
      new_list_element%field%l_ptr = .FALSE.
    END IF
    !
    IF (PRESENT(initval_l)) THEN
      new_list_element%field%l_ptr = new_list_element%field%info%initval%lval
    ENDIF
    !
  END SUBROUTINE add_var_list_element_l2d
  !------------------------------------------------------------------------------------------------
  !
  ! create (allocate) a new table entry
  ! optionally obtain pointer to 1d-field
  ! optionally overwrite default meta data 
  !
  SUBROUTINE add_var_list_element_l1d(this_list, name, ptr,    &
       hgrid, vgrid, cf, grib2, ldims, loutput,                &
       lrestart, lrestart_cont, initval_l, isteptype,          &
       resetval_l, lmiss, missval_l, info, p5, vert_interp,    &
       hor_interp, in_group, verbose, new_element, cdiTimeID,  &
       l_pp_scheduler_task, post_op)
    !
    TYPE(t_var_list),     INTENT(inout)        :: this_list           ! list
    CHARACTER(len=*),     INTENT(in)           :: name                ! name of variable
    LOGICAL,              POINTER              :: ptr(:)              ! reference to field
    INTEGER,              INTENT(in)           :: hgrid               ! horizontal grid type used
    INTEGER,              INTENT(in)           :: vgrid               ! vertical grid type used
    TYPE(t_cf_var),       INTENT(in)           :: cf                  ! CF related metadata
    TYPE(t_grib2_var),    INTENT(in)           :: grib2               ! GRIB2 related metadata
    INTEGER,              INTENT(in), OPTIONAL :: ldims(1)            ! local dimensions
    LOGICAL,              INTENT(in), OPTIONAL :: loutput             ! output flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart            ! restart flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart_cont       ! continue restart if var not available
    LOGICAL,              INTENT(in), OPTIONAL :: initval_l           ! value if var not available
    INTEGER,              INTENT(in), OPTIONAL :: isteptype           ! type of statistical processing
    LOGICAL,              INTENT(in), OPTIONAL :: resetval_l          ! reset value (after accumulation)
    LOGICAL,              INTENT(in), OPTIONAL :: lmiss               ! missing value flag
    LOGICAL,              INTENT(in), OPTIONAL :: missval_l           ! missing value
    TYPE(t_var_metadata), POINTER,    OPTIONAL :: info                ! returns reference to metadata
    LOGICAL,              POINTER,    OPTIONAL :: p5(:,:,:,:,:)       ! provided pointer
    TYPE(t_vert_interp_meta),INTENT(in), OPTIONAL :: vert_interp      ! vertical interpolation metadata
    TYPE(t_hor_interp_meta), INTENT(in), OPTIONAL :: hor_interp       ! horizontal interpolation metadata
    LOGICAL, INTENT(in), OPTIONAL :: in_group(SIZE(VAR_GROUPS))       ! groups to which a variable belongs
    LOGICAL,              INTENT(in), OPTIONAL :: verbose             ! print information
    TYPE(t_list_element), POINTER, OPTIONAL  :: new_element ! pointer to new var list element
    INTEGER,              INTENT(in), OPTIONAL :: cdiTimeID           ! CDI time mode (TIME_VARIABLE/TIME_CONSTANT)
    INTEGER,              INTENT(in), OPTIONAL :: l_pp_scheduler_task ! .TRUE., if field is updated by pp scheduler
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op             !< "post-op" (small arithmetic operations) for this variable
    !
    TYPE(t_list_element), POINTER :: new_list_element
    TYPE(t_union_vals) :: missval, initval, resetval
    INTEGER :: idims(5)
    INTEGER :: istat
    LOGICAL :: referenced
    !
    ! consistency check for restart and output
    !
    IF (PRESENT(lrestart)) THEN
      IF (.NOT. this_list%p%lrestart .AND. lrestart) THEN
        CALL finish('mo_var_list:add_var_list_element_l1d',                            &
             &      'for list '//TRIM(this_list%p%name)//' restarting not enabled, '// &
             &      'but restart of '//TRIM(name)//' requested.')
      ENDIF
    ENDIF
    !
    ! add list entry
    !
    CALL append_list_element (this_list, new_list_element)
    IF (PRESENT(new_element)) new_element=>new_list_element
    new_list_element%field%info = default_var_list_metadata(this_list)
    !
    ! and set meta data
    !
    IF (PRESENT(p5)) THEN
      referenced = .TRUE.
      new_list_element%field%info%allocated = .TRUE.
    ELSE
      referenced = .FALSE.
    ENDIF
    !
    CALL assign_if_present(missval%lval, missval_l)
    CALL assign_if_present(initval%lval, initval_l)
    CALL assign_if_present(resetval%lval, resetval_l)
    CALL set_var_metadata (new_list_element%field%info,                   &
         name=name, hgrid=hgrid, vgrid=vgrid,                             & 
         cf=cf, grib2=grib2, ldims=ldims, loutput=loutput,                &
         lrestart=lrestart, lrestart_cont=lrestart_cont, initval=initval, &
         isteptype=isteptype, resetval=resetval, lmiss=lmiss,             &
         missval=missval, vert_interp=vert_interp, hor_interp=hor_interp, &
         in_group=in_group, verbose=verbose, cdiTimeID=cdiTimeID,         &
         l_pp_scheduler_task=l_pp_scheduler_task,                         &
         post_op=post_op)
    !
    IF (.NOT. referenced) THEN
      CALL assign_if_present(new_list_element%field%info%used_dimensions(1:1), ldims(1:1))
      idims(1) = new_list_element%field%info%used_dimensions(1)
      idims(2) = 1
      idims(3) = 1
      idims(4) = 1
      idims(5) = 1
      new_list_element%field%info%ndims = 1
      new_list_element%field%info%cdiDataType = DATATYPE_INT8      
      new_list_element%field%var_base_size = 4
      ALLOCATE (new_list_element%field%l_ptr(idims(1), idims(2), idims(3), idims(4), idims(5)), &
           &    STAT=istat)
      IF (istat /= 0) THEN
        CALL finish('mo_var_list:add_var_list_element_l1d', &
             &      'allocation of array '//TRIM(name)//' failed')
      ELSE
        new_list_element%field%info%allocated = .TRUE.
      ENDIF
      NULLIFY(new_list_element%field%r_ptr)
      NULLIFY(new_list_element%field%i_ptr)      
      this_list%p%memory_used = this_list%p%memory_used &
           +INT(new_list_element%field%var_base_size*SIZE(new_list_element%field%l_ptr),i8)
    ELSE
      new_list_element%field%l_ptr => p5
    ENDIF
    ptr => new_list_element%field%l_ptr(:,1,1,1,1)
    IF(PRESENT(info)) info => new_list_element%field%info
    !
    IF (PRESENT(lmiss)) THEN
      new_list_element%field%l_ptr = new_list_element%field%info%missval%lval
    ELSE
      new_list_element%field%l_ptr = .FALSE.
    END IF
    !
    IF (PRESENT(initval_l)) THEN
      new_list_element%field%l_ptr = new_list_element%field%info%initval%lval
    ENDIF
    !
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
       &                                 lrestart, lrestart_cont, initval_r, isteptype,          &
       &                                 resetval_r, lmiss, missval_r, tlev_source, tracer_info, &
       &                                 info, vert_interp, hor_interp, in_group, verbose,       &
       &                                 new_element, cdiTimeID, l_pp_scheduler_task, post_op)
    !
    TYPE(t_var_list), INTENT(inout)            :: this_list
    CHARACTER(len=*), INTENT(in)               :: target_name
    CHARACTER(len=*), INTENT(in)               :: name
    REAL(wp), POINTER                          :: ptr(:,:,:)
    INTEGER,              INTENT(in)           :: hgrid               ! horizontal grid type used
    INTEGER,              INTENT(in)           :: vgrid               ! vertical grid type used
    TYPE(t_cf_var),       INTENT(in)           :: cf                  ! CF related metadata
    TYPE(t_grib2_var),    INTENT(in)           :: grib2               ! GRIB2 related metadata
    INTEGER,              INTENT(in), OPTIONAL :: ref_idx             ! idx of slice to be referenced
    INTEGER,              INTENT(in), OPTIONAL :: ldims(3)            ! local dimensions, for checking
    LOGICAL,              INTENT(in), OPTIONAL :: loutput             ! output flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart            ! restart flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart_cont       ! continue restart if var not available
    REAL(wp),             INTENT(in), OPTIONAL :: initval_r           ! value if var not available
    INTEGER,              INTENT(in), OPTIONAL :: isteptype           ! type of statistical processing
    REAL(wp),             INTENT(in), OPTIONAL :: resetval_r          ! reset value (after accumulation)
    LOGICAL,              INTENT(in), OPTIONAL :: lmiss               ! missing value flag
    REAL(wp),             INTENT(in), OPTIONAL :: missval_r           ! missing value
    INTEGER,              INTENT(in), OPTIONAL :: tlev_source         ! actual TL for TL dependent vars
    TYPE(t_tracer_meta),  INTENT(in), OPTIONAL :: tracer_info         ! tracer meta data
    TYPE(t_var_metadata), POINTER,    OPTIONAL :: info                ! returns reference to metadata
    TYPE(t_vert_interp_meta),INTENT(in), OPTIONAL :: vert_interp      ! vertical interpolation metadata
    TYPE(t_hor_interp_meta), INTENT(in), OPTIONAL :: hor_interp       ! horizontal interpolation metadata
    LOGICAL, INTENT(in), OPTIONAL :: in_group(SIZE(VAR_GROUPS))       ! groups to which a variable belongs
    LOGICAL,              INTENT(in), OPTIONAL :: verbose
    TYPE(t_list_element), POINTER, OPTIONAL  :: new_element ! pointer to new var list element
    INTEGER,              INTENT(in), OPTIONAL :: cdiTimeID           ! CDI time mode (TIME_VARIABLE/TIME_CONSTANT)
    INTEGER,              INTENT(in), OPTIONAL :: l_pp_scheduler_task ! .TRUE., if field is updated by pp scheduler
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op             !< "post-op" (small arithmetic operations) for this variable
    !
    REAL(wp), POINTER :: target_ptr3d(:,:,:)
    REAL(wp), POINTER :: target_ptr4d(:,:,:,:)
    !
    TYPE(t_list_element), POINTER :: target_element    
    TYPE(t_var_metadata), POINTER :: target_info, ref_info
    TYPE(t_list_element), POINTER :: new_list_element
    TYPE(t_union_vals) :: missval, initval, resetval
    !
    NULLIFY(target_ptr3d)
    NULLIFY(target_ptr4d)
    !
    target_element => find_list_element (this_list, target_name)
    target_info => target_element%field%info
    IF (target_info%lcontainer) THEN
      target_ptr4d => target_element%field%r_ptr(:,:,:,:,1)
      !
      IF (.NOT. ASSOCIATED(target_ptr4d)) THEN
        CALL finish('add_var_list_reference_r3d', &
             TRIM(name)//' not created.')
      ENDIF
      !
      ! Counting the number of existing references is deactivated, if the slice index 
      ! to be referenced is given explicitly.
      IF ( PRESENT(ref_idx) ) THEN
        ! only check validity of given slice index
        IF ( (ref_idx > SIZE(target_ptr4d,4)) .OR. (ref_idx < 1)) THEN
          WRITE (message_text, *) &
            &  'Slice idx ', ref_idx, ' for ', TRIM(name), &
            &  ' out of allowable range [1,',SIZE(target_ptr4d,4),']'      
          CALL finish('add_var_list_reference_r3d', message_text)
        ENDIF
      ELSE
        target_info%ncontained = target_info%ncontained+1
        IF (SIZE(target_ptr4d,4) < target_info%ncontained) THEN
          WRITE (message_text, *) &
            &  TRIM(name), ' exceeds the number of predefined entries in container:', &
            &  SIZE(target_ptr4d,4)      
          CALL finish('add_var_list_reference_r3d', message_text)
        ENDIF
      ENDIF
      IF ( ldims(1) /=  target_info%used_dimensions(1) .OR. &
           ldims(2) /=  target_info%used_dimensions(2) .OR. & 
           ldims(3) /=  target_info%used_dimensions(3) ) THEN
        CALL finish('add_var_list_reference_r3d', & 
             TRIM(name)//' dimensions requested and available differ.')      
      ENDIF
    ELSE
      target_ptr3d => target_element%field%r_ptr(:,:,:,1,1)
      !
      IF (.NOT. ASSOCIATED(target_ptr3d)) THEN
        CALL finish('add_var_list_reference_r3d', &
             TRIM(name)//' not created.')
      ENDIF
    ENDIF
    !
    ! add list entry
    !
    CALL append_list_element (this_list, new_list_element)
    IF (PRESENT(new_element)) new_element=>new_list_element
    ref_info => new_list_element%field%info
    ref_info = default_var_list_metadata(this_list)

    !
    CALL assign_if_present(missval%rval, missval_r)
    CALL assign_if_present(initval%rval, initval_r)
    CALL assign_if_present(resetval%rval, resetval_r)
    !
    CALL set_var_metadata (new_list_element%field%info,                     &
         name=name, hgrid=hgrid, vgrid=vgrid,                               & 
         cf=cf, grib2=grib2, ldims=ldims, loutput=loutput,                  &
         lrestart=lrestart, lrestart_cont=lrestart_cont, initval=initval,   &
         isteptype=isteptype, resetval=resetval, lmiss=lmiss,               &
         missval=missval, tlev_source=tlev_source, tracer_info=tracer_info, &
         vert_interp=vert_interp, hor_interp=hor_interp,                    &
         in_group=in_group, verbose=verbose, cdiTimeID=cdiTimeID,           &
         l_pp_scheduler_task=l_pp_scheduler_task,                           &
         post_op=post_op)
    !
    ref_info%ndims = 3
    ref_info%used_dimensions =  target_element%field%info%used_dimensions
    !
    IF (target_info%lcontainer) THEN
      ref_info%lcontained = .TRUE.
      ref_info%used_dimensions(4) = 1
      !
      IF ( PRESENT(ref_idx) ) THEN
        ref_info%ncontained = ref_idx
        ptr => target_element%field%r_ptr(:,:,:,ref_idx,1)
      ELSE
        ref_info%ncontained = target_info%ncontained
        ptr => target_element%field%r_ptr(:,:,:,target_info%ncontained,1)
      ENDIF
      !
      new_list_element%field%r_ptr => target_element%field%r_ptr
    ELSE
      new_list_element%field%r_ptr => target_element%field%r_ptr
      ptr => target_element%field%r_ptr(:,:,:,1,1)
    ENDIF
    !
    IF (.NOT. ASSOCIATED(new_list_element%field%r_ptr)) THEN
      write (0,*) 'problem with association of ptr for '//TRIM(name)
    ENDIF
    !
    IF(PRESENT(info)) info => new_list_element%field%info
    !
    IF (PRESENT(lmiss)) THEN
      ptr = new_list_element%field%info%missval%rval
    ELSE
      ptr = 0.0_wp
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
       &                                 lrestart, lrestart_cont, initval_r, isteptype,          &
       &                                 resetval_r, lmiss, missval_r, tlev_source, tracer_info, &
       &                                 info, vert_interp, hor_interp, in_group,                &
       &                                 verbose, new_element, cdiTimeID, l_pp_scheduler_task,   &
       &                                 post_op)

    TYPE(t_var_list), INTENT(inout)            :: this_list
    CHARACTER(len=*), INTENT(in)               :: target_name
    CHARACTER(len=*), INTENT(in)               :: name
    REAL(wp), POINTER                          :: ptr(:,:)
    INTEGER,              INTENT(in)           :: hgrid               ! horizontal grid type used
    INTEGER,              INTENT(in)           :: vgrid               ! vertical grid type used
    TYPE(t_cf_var),       INTENT(in)           :: cf                  ! CF related metadata
    TYPE(t_grib2_var),    INTENT(in)           :: grib2               ! GRIB2 related metadata
    INTEGER,              INTENT(in), OPTIONAL :: ref_idx             ! idx of slice to be referenced
    INTEGER,              INTENT(in), OPTIONAL :: ldims(2)            ! local dimensions, for checking
    LOGICAL,              INTENT(in), OPTIONAL :: loutput             ! output flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart            ! restart flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart_cont       ! continue restart if var not available
    REAL(wp),             INTENT(in), OPTIONAL :: initval_r           ! value if var not available
    INTEGER,              INTENT(in), OPTIONAL :: isteptype           ! type of statistical processing
    REAL(wp),             INTENT(in), OPTIONAL :: resetval_r          ! reset value (after accumulation)
    LOGICAL,              INTENT(in), OPTIONAL :: lmiss               ! missing value flag
    REAL(wp),             INTENT(in), OPTIONAL :: missval_r           ! missing value
    INTEGER,              INTENT(in), OPTIONAL :: tlev_source         ! actual TL for TL dependent vars
    TYPE(t_tracer_meta),  INTENT(in), OPTIONAL :: tracer_info         ! tracer meta data
    TYPE(t_var_metadata), POINTER,    OPTIONAL :: info                ! returns reference to metadata
    TYPE(t_vert_interp_meta),INTENT(in), OPTIONAL :: vert_interp      ! vertical interpolation metadata
    TYPE(t_hor_interp_meta), INTENT(in), OPTIONAL :: hor_interp       ! horizontal interpolation metadata
    LOGICAL, INTENT(in), OPTIONAL :: in_group(SIZE(VAR_GROUPS))       ! groups to which a variable belongs
    LOGICAL,              INTENT(in), OPTIONAL :: verbose
    TYPE(t_list_element), POINTER, OPTIONAL  :: new_element ! pointer to new var list element
    INTEGER,              INTENT(in), OPTIONAL :: cdiTimeID           ! CDI time mode (TIME_VARIABLE/TIME_CONSTANT)
    INTEGER,              INTENT(in), OPTIONAL :: l_pp_scheduler_task ! .TRUE., if field is updated by pp scheduler
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op             !< "post-op" (small arithmetic operations) for this variable
    !
    REAL(wp), POINTER :: target_ptr2d(:,:)
    REAL(wp), POINTER :: target_ptr3d(:,:,:)
    !
    TYPE(t_list_element), POINTER :: target_element    
    TYPE(t_var_metadata), POINTER :: target_info, ref_info
    TYPE(t_list_element), POINTER :: new_list_element
    TYPE(t_union_vals) :: missval, initval, resetval
    !
    NULLIFY(target_ptr2d)
    NULLIFY(target_ptr3d)
    !
    target_element => find_list_element (this_list, target_name)
    target_info => target_element%field%info
    IF (target_info%lcontainer) THEN
      target_ptr3d => target_element%field%r_ptr(:,:,:,1,1)
      !
      IF (.NOT. ASSOCIATED(target_ptr3d)) THEN
        CALL finish('add_var_list_reference_r2d', &
             TRIM(name)//' not created.')
      ENDIF
      !
      ! Counting the number of existing references is deactivated, if the slice index 
      ! to be referenced is given explicitly.
      IF ( PRESENT(ref_idx) ) THEN
        ! only check validity of given slice index
        IF ( (ref_idx > SIZE(target_ptr3d,3)) .OR. (ref_idx < 1)) THEN
          WRITE (message_text, *) &
            &  'Slice idx ', ref_idx, ' for ', TRIM(name), &
            &  ' out of allowable range [1,',SIZE(target_ptr3d,3),']'      
          CALL finish('add_var_list_reference_r2d', message_text)
        ENDIF
      ELSE
        target_info%ncontained = target_info%ncontained+1
        IF (SIZE(target_ptr3d,3) < target_info%ncontained) THEN
          WRITE (message_text, *) &
            &  TRIM(name), ' exceeds the number of predefined entries in container:', &
            &  SIZE(target_ptr3d,3)
          CALL finish('add_var_list_reference_r2d', message_text)
        ENDIF
      ENDIF
      IF ( ldims(1) /=  target_info%used_dimensions(1) .OR. &
           ldims(2) /=  target_info%used_dimensions(2) ) THEN
        CALL finish('add_var_list_reference_r2d', & 
             TRIM(name)//' dimensions requested and available differ.')      
      ENDIF
    ELSE
      target_ptr2d => target_element%field%r_ptr(:,:,1,1,1)
      !
      IF (.NOT. ASSOCIATED(target_ptr2d)) THEN
        CALL finish('add_var_list_reference_r2d', &
             TRIM(name)//' not created.')
      ENDIF
    ENDIF
    !
    ! add list entry
    !
    CALL append_list_element (this_list, new_list_element)
    IF (PRESENT(new_element)) new_element=>new_list_element
    ref_info => new_list_element%field%info
    ref_info = default_var_list_metadata(this_list)

    !
    CALL assign_if_present(missval%rval, missval_r)
    CALL assign_if_present(initval%rval, initval_r)
    CALL assign_if_present(resetval%rval, resetval_r)
    !
    CALL set_var_metadata (new_list_element%field%info,                     &
         name=name, hgrid=hgrid, vgrid=vgrid,                               & 
         cf=cf, grib2=grib2, ldims=ldims, loutput=loutput,                  &
         lrestart=lrestart, lrestart_cont=lrestart_cont, initval=initval,   &
         isteptype=isteptype, resetval=resetval, lmiss=lmiss,               &
         missval=missval, tlev_source=tlev_source, tracer_info=tracer_info, &
         vert_interp=vert_interp, hor_interp=hor_interp,                    &
         in_group=in_group, verbose=verbose, cdiTimeID=cdiTimeID,           &
         l_pp_scheduler_task=l_pp_scheduler_task,                           &
         post_op=post_op)
    !
    ref_info%ndims = 2
    ref_info%used_dimensions = target_element%field%info%used_dimensions
    !
    IF (target_info%lcontainer) THEN
      ref_info%lcontained = .TRUE.
      ref_info%used_dimensions(3) = 1
      !
      IF ( PRESENT(ref_idx) ) THEN
        ref_info%ncontained = ref_idx
        ptr => target_element%field%r_ptr(:,:,ref_idx,1,1)
      ELSE
        ref_info%ncontained = target_info%ncontained
        ptr => target_element%field%r_ptr(:,:,target_info%ncontained,1,1)
      ENDIF
      !
      new_list_element%field%r_ptr => target_element%field%r_ptr
    ELSE
      new_list_element%field%r_ptr => target_element%field%r_ptr
      ptr => target_element%field%r_ptr(:,:,1,1,1)
    ENDIF
    !
    IF (.NOT. ASSOCIATED(new_list_element%field%r_ptr)) THEN
      write (0,*) 'problem with association of ptr for '//TRIM(name)
    ENDIF
    !
    IF(PRESENT(info)) info => new_list_element%field%info
    !
    IF (PRESENT(lmiss)) THEN
      ptr = new_list_element%field%info%missval%rval
    ELSE
      ptr = 0.0_wp
    END IF
    !
  END SUBROUTINE add_var_list_reference_r2d


  !================================================================================================
  !------------------------------------------------------------------------------------------------
  !
  ! add supplementary fields to a different var list (eg. geopotential, surface pressure, ...)
  !
  SUBROUTINE add_var_list_reference (to_var_list, name, from_var_list, loutput, bit_precision, in_group)
    TYPE(t_var_list), INTENT(inout)          :: to_var_list
    CHARACTER(len=*), INTENT(in)             :: name
    CHARACTER(len=*), INTENT(in)             :: from_var_list
    LOGICAL,          INTENT(in),   OPTIONAL :: loutput
    INTEGER,          INTENT(in),   OPTIONAL :: bit_precision
    LOGICAL,          INTENT(in),   OPTIONAL :: in_group(SIZE(VAR_GROUPS))  ! groups to which a variable belongs
    !
    TYPE(t_var_list_element), POINTER :: source
    TYPE(t_list_element),     POINTER :: new_list_element
    !
    CALL locate (source, name, from_var_list)
    IF (ASSOCIATED(source)) THEN
      CALL append_list_element (to_var_list, new_list_element)
      new_list_element%field                = source
      new_list_element%field%info%allocated = .FALSE.
      new_list_element%field%info%lrestart  = .FALSE.
      CALL assign_if_present(new_list_element%field%info%loutput, loutput)
      CALL assign_if_present(new_list_element%field%info%grib2%bits, bit_precision)
      if (present(in_group)) then
        new_list_element%field%info%in_group(:)=in_group(:)
      end if
    ENDIF
    !
  CONTAINS
    !----------------------------------------------------------------------------------------------
    !
    ! find an entry
    !
    SUBROUTINE locate (element, name, in_var_list)
      TYPE(t_var_list_element), POINTER        :: element
      CHARACTER(len=*), INTENT(in)           :: name
      CHARACTER(len=*), INTENT(in), OPTIONAL :: in_var_list
      !
      INTEGER                     :: i
      TYPE(t_list_element), POINTER :: link
      !
      NULLIFY (element)
      !
      DO i = 1, nvar_lists
        IF (PRESENT(in_var_list)) THEN
          IF (in_var_list /= var_lists(i)%p%name) CYCLE
        ENDIF
        link => find_list_element (var_lists(i), name)
        IF (ASSOCIATED(link)) THEN
          element => link%field
          EXIT
        ENDIF
      END DO
      !
    END SUBROUTINE locate
    !
  END SUBROUTINE add_var_list_reference
  !------------------------------------------------------------------------------------------------
  !
  ! remove one element from the list
  ! the element is identified by its name
  !
  SUBROUTINE delete_var_list_element (this_list, name)
    TYPE(t_var_list), INTENT(inout) :: this_list
    CHARACTER(len=*), INTENT(in)    :: name   
    !
    TYPE(t_list_element), POINTER :: ptr
    !
    IF (this_list%p%first_list_element%field%info%name == name) THEN
      CALL delete_list_element (this_list, this_list%p%first_list_element)
      RETURN
    ELSE
      ptr => this_list%p%first_list_element
      DO
        IF (.NOT.ASSOCIATED (ptr%next_list_element)) EXIT
        IF (ptr%next_list_element%field%info%name == name) THEN
          CALL delete_list_element (this_list, ptr%next_list_element)
          EXIT
        ENDIF
        ptr => ptr%next_list_element
      END DO
    ENDIF
    !
  END SUBROUTINE delete_var_list_element
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
  SUBROUTINE print_var_list (this_list)
    TYPE(t_var_list),  INTENT(in) :: this_list ! list
    !
    TYPE(t_list_element), POINTER :: this_list_element
    CHARACTER(len=32) :: dimension_text, dtext
    INTEGER :: i, igrp, ivintp_type
   
    CALL message('','')
    CALL message('','')
    CALL message('','Status of variable list '//TRIM(this_list%p%name)//':')    
    CALL message('','')
    !
    this_list_element => this_list%p%first_list_element
    !
    DO WHILE (ASSOCIATED(this_list_element))
      ! 
      IF (this_list_element%field%info%name /= '' .AND. &
           .NOT. this_list_element%field%info%lcontainer) THEN
        !
        WRITE (message_text,'(a,a)')       &
             'Table entry name                            : ', &
             TRIM(this_list_element%field%info%name)
        CALL message('', message_text)
        !
        IF (ASSOCIATED(this_list_element%field%r_ptr) .OR. &
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
!DR Due to problems with insufficient record length on SX9 
        WRITE (0, * ) &
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
        IF (this_list_element%field%info%tracer%lis_tracer) THEN
          CALL message('', 'Tracer field                                : yes.')

          WRITE (message_text,'(a,a)') &
             'Tracer class                                : ', &
             this_list_element%field%info%tracer%tracer_class
          CALL message('', message_text)


          WRITE (message_text,'(a,3i3)') &
             'Horizontal transport method                 : ', &
             this_list_element%field%info%tracer%ihadv_tracer
          CALL message('', message_text)

          WRITE (message_text,'(a,3i3)') &
             'Vertical transport method                   : ', &
             this_list_element%field%info%tracer%ivadv_tracer
          CALL message('', message_text)

          IF (this_list_element%field%info%tracer%lturb_tracer) THEN
            CALL message('', 'Turbulent transport                         : yes.')
          ELSE
            CALL message('', 'Turbulent transport                         : no.')
          ENDIF

          IF (this_list_element%field%info%tracer%lsed_tracer) THEN
            CALL message('', 'Sedimentation                               : yes.')
          ELSE
            CALL message('', 'Sedimentation                               : no.')
          ENDIF

          IF (this_list_element%field%info%tracer%ldep_tracer) THEN
            CALL message('', 'Dry deposition                              : yes.')
          ELSE
            CALL message('', 'Dry deposition                              : no.')
          ENDIF

          IF (this_list_element%field%info%tracer%lconv_tracer) THEN
            CALL message('', 'Convection                                  : yes.')
          ELSE
            CALL message('', 'Convection                                  : no.')
          ENDIF

          IF (this_list_element%field%info%tracer%lwash_tracer) THEN
            CALL message('', 'Washout                                     : yes.')
          ELSE
            CALL message('', 'Washout                                     : no.')
          ENDIF

          WRITE (message_text,'(a,e18.12)') &
             'Particle diameter in m                      : ', &
             this_list_element%field%info%tracer%rdiameter_tracer
          CALL message('', message_text)

          WRITE (message_text,'(a,e18.12)') &
             'particle density in kg m^-3                 : ', &
             this_list_element%field%info%tracer%rrho_tracer
          CALL message('', message_text)

        WRITE (message_text,'(a,e18.12)') &
             'Radioactive half-life in s^-1                      : ', &
             this_list_element%field%info%tracer%halflife_tracer
          CALL message('', message_text)

        WRITE (message_text,'(a,i3)') &
             'IMIS number                      : ', &
             this_list_element%field%info%tracer%imis_tracer
          CALL message('', message_text)


        ELSE
          CALL message('', 'Tracer field                                : no.')
        ENDIF

        ! 
        ! print groups, to which this variable belongs:
        IF (ANY(this_list_element%field%info%in_group(:))) THEN
          WRITE (message_text,'(a)')  'Variable group(s)                           :'
          DO igrp=1,SIZE(var_groups)
            IF (this_list_element%field%info%in_group(igrp)) THEN
              IF (igrp == 1) THEN
                message_text = TRIM(message_text)//" "//TRIM(var_groups(igrp))
              ELSE
                message_text = TRIM(message_text)//", "//TRIM(var_groups(igrp))
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

      !
      ! select next element in linked list 
      !
      this_list_element => this_list_element%next_list_element
    ENDDO

    !    
  END SUBROUTINE print_var_list
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


  !------------------------------------------------------------------------------------------------
  ! HANDLING OF VERTICAL INTERPOLATION MODES
  !------------------------------------------------------------------------------------------------

  !> Implements a (somewhat randomly chosen) one-to-one mapping
  !  between a string and an integer ID number between 1 and
  !  MAX_VINTP_TYPES.
  !
  FUNCTION vintp_type_id(in_str)
    INTEGER                      :: vintp_type_id, ivintp_type
    CHARACTER(LEN=*), INTENT(IN) :: in_str
    CHARACTER(*), PARAMETER :: routine = TRIM("mo_var_list:vintp_type_id")

    vintp_type_id = 0
    LOOP_VINTP_TYPES : DO ivintp_type=1,SIZE(VINTP_TYPE_LIST)
      IF (toupper(TRIM(in_str)) == toupper(TRIM(VINTP_TYPE_LIST(ivintp_type)))) THEN
        vintp_type_id = ivintp_type
        EXIT LOOP_VINTP_TYPES
      END IF
    END DO LOOP_VINTP_TYPES
    ! paranoia:
    IF ((vintp_type_id < 1) .OR. (vintp_type_id > SIZE(VINTP_TYPE_LIST))) &
      &  CALL finish(routine, "Invalid vertical interpolation type!")
  END FUNCTION vintp_type_id


  !> Utility function with *a lot* of optional string parameters v1,
  !  v2, v3, v4, ...; mapping those onto a
  !  LOGICAL(DIMENSION=MAX_VAR_GROUPS) according to the "group_id"
  !  function.
  !
  FUNCTION vintp_types(v01, v02, v03, v04, v05, v06, v07, v08, v09, v10)
    LOGICAL :: vintp_types(SIZE(VINTP_TYPE_LIST))
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: &
      &   v01, v02, v03, v04, v05, v06, v07, v08, v09, v10
    
    vintp_types(:) = .FALSE.
    IF (PRESENT(v01)) vintp_types(vintp_type_id(v01)) = .TRUE.
    IF (PRESENT(v02)) vintp_types(vintp_type_id(v02)) = .TRUE.
    IF (PRESENT(v03)) vintp_types(vintp_type_id(v03)) = .TRUE.
    IF (PRESENT(v04)) vintp_types(vintp_type_id(v04)) = .TRUE.
    IF (PRESENT(v05)) vintp_types(vintp_type_id(v05)) = .TRUE.
    IF (PRESENT(v06)) vintp_types(vintp_type_id(v06)) = .TRUE.
    IF (PRESENT(v07)) vintp_types(vintp_type_id(v07)) = .TRUE.
    IF (PRESENT(v08)) vintp_types(vintp_type_id(v08)) = .TRUE.
    IF (PRESENT(v09)) vintp_types(vintp_type_id(v09)) = .TRUE.
    IF (PRESENT(v10)) vintp_types(vintp_type_id(v10)) = .TRUE.
  END FUNCTION vintp_types


  !------------------------------------------------------------------------------------------------
  ! HANDLING OF VARIABLE GROUPS
  !------------------------------------------------------------------------------------------------
  
  !> Implements a (somewhat randomly chosen) one-to-one mapping
  !  between a string and an integer ID number between 1 and
  !  MAX_VAR_GROUPS.
  !
  FUNCTION group_id(in_str)
    INTEGER                      :: group_id, igrp
    CHARACTER(LEN=*), INTENT(IN) :: in_str
    CHARACTER(*), PARAMETER :: routine = TRIM("mo_var_list:group_id")

    group_id = 0
    LOOP_GROUPS : DO igrp=1,SIZE(VAR_GROUPS)
      IF (toupper(TRIM(in_str)) == toupper(TRIM(VAR_GROUPS(igrp)))) THEN
        group_id = igrp
        EXIT LOOP_GROUPS
      END IF
    END DO LOOP_GROUPS
    ! paranoia:
    IF ((group_id < 1) .OR. (group_id > SIZE(VAR_GROUPS))) &
      &  CALL finish(routine, "Invalid group ID!")

  END FUNCTION group_id

  !> Utility function with *a lot* of optional string parameters g1,
  !  g2, g3, g4, ...; mapping those onto a
  !  LOGICAL(DIMENSION=MAX_VAR_GROUPS) according to the "group_id"
  !  function.
  !
  FUNCTION groups(g01, g02, g03, g04, g05, g06, g07, g08, g09, g10)
    LOGICAL :: groups(SIZE(VAR_GROUPS))
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: &
      &   g01, g02, g03, g04, g05, g06, g07, g08, g09, g10
    
    groups(:) = .FALSE.
    groups(group_id("ALL")) = .TRUE.
    IF (PRESENT(g01)) groups(group_id(g01)) = .TRUE.
    IF (PRESENT(g02)) groups(group_id(g02)) = .TRUE.
    IF (PRESENT(g03)) groups(group_id(g03)) = .TRUE.
    IF (PRESENT(g04)) groups(group_id(g04)) = .TRUE.
    IF (PRESENT(g05)) groups(group_id(g05)) = .TRUE.
    IF (PRESENT(g06)) groups(group_id(g06)) = .TRUE.
    IF (PRESENT(g07)) groups(group_id(g07)) = .TRUE.
    IF (PRESENT(g08)) groups(group_id(g08)) = .TRUE.
    IF (PRESENT(g09)) groups(group_id(g09)) = .TRUE.
    IF (PRESENT(g10)) groups(group_id(g10)) = .TRUE.
  END FUNCTION groups

  !> Loops over all variables and collects the variables names
  !  corresponding to the group @p grp_name
  !
  SUBROUTINE collect_group(grp_name, var_name, nvars, &
    &                      loutputvars_only, lremap_lonlat)
    CHARACTER(LEN=*),           INTENT(IN)    :: grp_name
    CHARACTER(LEN=VARNAME_LEN), INTENT(INOUT) :: var_name(:)
    INTEGER,                    INTENT(OUT)   :: nvars
    ! loutputvars_only: If set to .TRUE. all variables in the group
    ! which have the the loutput flag equal to .FALSE. are skipped.
    LOGICAL,                    INTENT(IN)    :: loutputvars_only
    ! lremap_lonlat: If set to .TRUE. only variables in the group
    ! which can be interpolated onto lon-lat grids are considered.
    LOGICAL,                    INTENT(IN)    :: lremap_lonlat

    ! local variables
    CHARACTER(*), PARAMETER :: routine = TRIM("mo_var_list:collect_group")
    INTEGER :: i, ivar, grp_id
    TYPE(t_list_element), POINTER :: element
    TYPE(t_var_metadata), POINTER :: info
    CHARACTER(LEN=VARNAME_LEN)    :: name
    
    nvars  = 0
    grp_id = group_id(grp_name)

    ! loop over all variable lists and variables
    DO i = 1,nvar_lists
      element => NULL()
      LOOPVAR : DO
        IF(.NOT.ASSOCIATED(element)) THEN
          element => var_lists(i)%p%first_list_element
        ELSE
          element => element%next_list_element
        ENDIF
        IF(.NOT.ASSOCIATED(element)) EXIT LOOPVAR
        info => element%field%info
        ! Do not inspect element if it is a container
        IF (info%lcontainer) CYCLE LOOPVAR
        
        IF (info%in_group(grp_id)) THEN
          name = TRIM(get_var_name(element%field))

          ! Skip element if we need only output variables:
          IF (loutputvars_only .AND. &
            & (.NOT. info%loutput) .OR. (.NOT. var_lists(i)%p%loutput)) THEN
            CALL message(routine, "Skipping variable "//TRIM(name)//" for output.")
            CYCLE LOOPVAR
          END IF

          IF (lremap_lonlat .AND. &
            & (info%hgrid /= GRID_UNSTRUCTURED_CELL)) THEN
            CALL message(routine, "Skipping variable "//TRIM(name)//" for lon-lat output.")
            CYCLE LOOPVAR
          END IF

          nvars = nvars + 1
          var_name(nvars) = name
        END IF
      ENDDO LOOPVAR ! loop over vlist "i"
    ENDDO ! i = 1,nvar_lists

    CALL remove_duplicates(var_name, nvars)

  END SUBROUTINE collect_group


  !> Utility function with *a lot* of optional string parameters g1,
  !  g2, g3, g4, ...; mapping those onto a
  !  LOGICAL(DIMENSION=MAX_VAR_GROUPS) according to the "group_id"
  !  function.
  !
  FUNCTION post_op(ipost_op_type, new_cf, new_grib2, arg1)
    TYPE(t_post_op_meta) :: post_op
    INTEGER,           INTENT(IN), OPTIONAL :: ipost_op_type    !< type of post-processing operation
    TYPE(t_cf_var),    INTENT(IN), OPTIONAL :: new_cf           !< CF information of modified field
    TYPE(t_grib2_var), INTENT(IN), OPTIONAL :: new_grib2        !< GRIB2 information of modified field
    REAL(wp),          INTENT(IN), OPTIONAL :: arg1             !< post-op argument (e.g. scaling factor)

    post_op%ipost_op_type = POST_OP_NONE
    post_op%lnew_cf       = .FALSE.
    post_op%lnew_grib2    = .FALSE.
    post_op%arg1          = 0._wp

    IF (PRESENT(ipost_op_type)) post_op%ipost_op_type = ipost_op_type
    IF (PRESENT(arg1))          post_op%arg1          = arg1
    IF (PRESENT(new_cf)) THEN
      post_op%lnew_cf = .TRUE.
      post_op%new_cf  = new_cf
    END IF
    IF (PRESENT(new_grib2)) THEN
      post_op%lnew_grib2 = .TRUE.
      post_op%new_grib2  = new_grib2
    END IF
  END FUNCTION post_op


  !------------------------------------------------------------------------------------------------
  SUBROUTINE assign_if_present_cf (y,x)
    TYPE(t_cf_var), INTENT(inout)        :: y
    TYPE(t_cf_var), INTENT(in) ,OPTIONAL :: x
    IF (.NOT.PRESENT(x)) RETURN
    y = x
  END SUBROUTINE assign_if_present_cf
  !------------------------------------------------------------------------------------------------
  SUBROUTINE assign_if_present_grib2 (y,x)
    TYPE(t_grib2_var), INTENT(inout)        :: y
    TYPE(t_grib2_var) ,INTENT(in) ,OPTIONAL :: x
    IF (.NOT.PRESENT(x)) RETURN
    y = x
  END SUBROUTINE assign_if_present_grib2
  !------------------------------------------------------------------------------------------------
  SUBROUTINE assign_if_present_union (y,x)
    TYPE(t_union_vals), INTENT(inout)        :: y
    TYPE(t_union_vals) ,INTENT(in) ,OPTIONAL :: x
    IF (.NOT.PRESENT(x)) RETURN
    y = x
  END SUBROUTINE assign_if_present_union
  !------------------------------------------------------------------------------------------------
  SUBROUTINE assign_if_present_tracer_meta (y,x)
    TYPE(t_tracer_meta), INTENT(inout)        :: y
    TYPE(t_tracer_meta) ,INTENT(in) ,OPTIONAL :: x
    IF (.NOT.PRESENT(x)) RETURN
    y = x
  END SUBROUTINE assign_if_present_tracer_meta
  !------------------------------------------------------------------------------------------------
  SUBROUTINE assign_if_present_vert_interp (y,x)
    TYPE(t_vert_interp_meta), INTENT(inout)        :: y
    TYPE(t_vert_interp_meta) ,INTENT(in) ,OPTIONAL :: x
    IF (.NOT.PRESENT(x)) RETURN
    y = x
  END SUBROUTINE assign_if_present_vert_interp
  !------------------------------------------------------------------------------------------------
  SUBROUTINE assign_if_present_hor_interp (y,x)
    TYPE(t_hor_interp_meta), INTENT(inout)        :: y
    TYPE(t_hor_interp_meta) ,INTENT(in) ,OPTIONAL :: x
    IF (.NOT.PRESENT(x)) RETURN
    y = x
  END SUBROUTINE assign_if_present_hor_interp
  !------------------------------------------------------------------------------------------------
  SUBROUTINE assign_if_present_post_op (y,x)
    TYPE(t_post_op_meta), INTENT(inout)        :: y
    TYPE(t_post_op_meta) ,INTENT(in) ,OPTIONAL :: x
    IF (.NOT.PRESENT(x)) RETURN
    y = x
  END SUBROUTINE assign_if_present_post_op
  !------------------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------------------
  !
  ! create (allocate) a new table entry
  ! reference to an existing pointer to a 3D tracer field
  ! optionally overwrite some default meta data
  !
  SUBROUTINE add_var_list_reference_tracer(this_list, target_name, tracer_name,  &
    &        tracer_idx, ptr_arr, cf, grib2, advconf,jg, ldims, loutput, lrestart,  &
    &        isteptype, tlev_source, vert_interp, hor_interp, in_group,          &
    &        lis_tracer,tracer_class,                                            &
    &        ihadv_tracer, ivadv_tracer, lturb_tracer, lsed_tracer,              &
    &        ldep_tracer, lconv_tracer, lwash_tracer, rdiameter_tracer,          &
    &        rrho_tracer,halflife_tracer,imis_tracer)

    TYPE(t_var_list)    , INTENT(inout)        :: this_list
    CHARACTER(len=*)    , INTENT(in)           :: target_name
    CHARACTER(len=*)    , INTENT(in)           :: tracer_name
    INTEGER             , INTENT(inout)        :: tracer_idx     ! index in 4D tracer container
    TYPE(t_ptr_2d3d)    , INTENT(inout)        :: ptr_arr(:)
    TYPE(t_cf_var)      , INTENT(in)           :: cf             ! CF related metadata
    TYPE(t_grib2_var)   , INTENT(in)           :: grib2          ! GRIB2 related metadata
    TYPE(t_advection_config), INTENT(inout)    :: advconf        ! adv configure state
    INTEGER             , INTENT(in), OPTIONAL :: jg             ! patch id
    INTEGER             , INTENT(in), OPTIONAL :: ldims(3)       ! local dimensions, for checking
    LOGICAL             , INTENT(in), OPTIONAL :: loutput        ! output flag
    LOGICAL             , INTENT(in), OPTIONAL :: lrestart       ! restart flag
    INTEGER,              INTENT(in), OPTIONAL :: isteptype      ! type of statistical processing
    INTEGER             , INTENT(in), OPTIONAL :: tlev_source    ! actual TL for TL dependent vars
    TYPE(t_vert_interp_meta),INTENT(in), OPTIONAL :: vert_interp ! vertical interpolation metadata
    TYPE(t_hor_interp_meta), INTENT(in), OPTIONAL :: hor_interp  ! horizontal interpolation metadata
    LOGICAL, INTENT(in), OPTIONAL :: in_group(SIZE(VAR_GROUPS))  ! groups to which a variable belongs
    LOGICAL             , INTENT(in), OPTIONAL :: lis_tracer     ! this is a tracer field (TRUE/FALSE)
    CHARACTER(len=*)    , INTENT(in), OPTIONAL :: tracer_class  ! type of tracer (cloud, volcash, radioact,...) 
    INTEGER             , INTENT(in), OPTIONAL :: ihadv_tracer   ! method for hor. transport
    INTEGER             , INTENT(in), OPTIONAL :: ivadv_tracer   ! method for vert. transport
    LOGICAL             , INTENT(in), OPTIONAL :: lturb_tracer   ! turbulent transport (TRUE/FALSE)
    LOGICAL             , INTENT(in), OPTIONAL :: lsed_tracer    ! sedimentation (TRUE/FALSE)
    LOGICAL             , INTENT(in), OPTIONAL :: ldep_tracer    ! dry deposition (TRUE/FALSE)
    LOGICAL             , INTENT(in), OPTIONAL :: lconv_tracer   ! convection  (TRUE/FALSE)
    LOGICAL             , INTENT(in), OPTIONAL :: lwash_tracer   ! washout (TRUE/FALSE)
    REAL(wp)            , INTENT(in), OPTIONAL :: rdiameter_tracer ! particle diameter in m
    REAL(wp)            , INTENT(in), OPTIONAL :: rrho_tracer    ! particle density in kg m^-3
    REAL(wp)            , INTENT(in), OPTIONAL :: halflife_tracer       ! radioactive half-life in s^-1
    INTEGER             , INTENT(in), OPTIONAL :: imis_tracer         ! IMIS number


    ! Local variables:
    TYPE(t_list_element), POINTER :: target_element
    TYPE(t_var_metadata), POINTER :: target_info
    TYPE(t_tracer_meta)           :: tracer_info    ! tracer meta data

    INTEGER :: zihadv_tracer, zivadv_tracer

    CHARACTER(*), PARAMETER :: routine = "add_tracer_ref"

  !------------------------------------------------------------------

    ! get pointer to target element (in this case 4D tracer container)
    target_element => find_list_element (this_list, target_name)
    ! get tracer field metadata
    target_info => target_element%field%info

    ! get index of current field in 4D container and set
    ! tracer index accordingly.
    ! Note that this will be repeated for each patch. Otherwise it may happen that
    ! some MPI-processes miss this assignment.
    !
    tracer_idx = target_info%ncontained+1  ! index in 4D tracer container

    WRITE (message_text,'(a,i3,a,a)')                                            &
      & "tracer index ", tracer_idx," assigned to tracer field ", TRIM(tracer_name)
    CALL message(TRIM(routine),message_text)



    ! set default values
    ! If ihadv_tracer or ivadv_tracer are not present, take respective value from the 
    ! advection_config state.
    ! If ihadv_tracer or ivadv_tracer are present, take those values, and overwrite 
    ! the respective values of the advection_config state. 
    !
    IF ( PRESENT( ihadv_tracer )) THEN
      zihadv_tracer = ihadv_tracer
      ! BE AWARE, THAT ivadv_tracer IS NOT SANITY-CHECKED. THIS OVERWRITES THE 
      ! SANITY CHECKED NAMELIST SETTINGS.
      advconf%ihadv_tracer(tracer_idx) = ihadv_tracer
    ELSE
      zihadv_tracer = advconf%ihadv_tracer(tracer_idx)
    ENDIF

    IF ( PRESENT( ivadv_tracer )) THEN
      zivadv_tracer = ivadv_tracer
      ! BE AWARE, THAT ivadv_tracer IS NOT SANITY-CHECKED. THIS OVERWRITES THE 
      ! SANITY CHECKED NAMELIST SETTINGS.
      advconf%ivadv_tracer(tracer_idx) = ivadv_tracer
    ELSE
      zivadv_tracer = advconf%ivadv_tracer(tracer_idx)
    ENDIF

    ! generate tracer metadata information
    !
    tracer_info = create_tracer_metadata(                                     &
      &                                  lis_tracer       = lis_tracer,       &
      &                                  tracer_class     = tracer_class,     &
      &                                  ihadv_tracer     = zihadv_tracer,    &
      &                                  ivadv_tracer     = zivadv_tracer,    &
      &                                  lturb_tracer     = lturb_tracer,     &
      &                                  lsed_tracer      = lsed_tracer,      &
      &                                  ldep_tracer      = ldep_tracer,      &
      &                                  lconv_tracer     = lconv_tracer,     &
      &                                  lwash_tracer     = lwash_tracer,     &
      &                                  rdiameter_tracer = rdiameter_tracer, &
      &                                  rrho_tracer      = rrho_tracer,      &
      &                                  halflife_tracer  = halflife_tracer,        &
      &                                  imis_tracer      = imis_tracer       )



    ! create new table entry reference including additional tracer metadata
    CALL add_ref( this_list, target_name, tracer_name, ptr_arr(tracer_idx)%p_3d, &
       &          target_info%hgrid, target_info%vgrid, cf, grib2,               &
       &          ldims=ldims, loutput=loutput, lrestart=lrestart,               &
       &          isteptype=isteptype, tlev_source=tlev_source,                  &
       &          vert_interp=vert_interp, hor_interp=hor_interp,                &
       &          tracer_info=tracer_info, in_group=in_group                     )

    ! Get the number of convection tracers  

    IF (art_config(jg)%lart_conv) THEN
         IF( lconv_tracer) THEN
      art_config(jg)%nconv_tracer = art_config(jg)%nconv_tracer + 1
         ELSE
          WRITE (message_text,*) 'WARNING. lart_conv=',&
          & art_config(jg)%lart_conv,&
          & 'lconv_tracer=',lconv_tracer
          CALL message('add_var_list_reference_tracer',message_text)
         ENDIF
    ENDIF

  END SUBROUTINE add_var_list_reference_tracer


END MODULE mo_var_list
