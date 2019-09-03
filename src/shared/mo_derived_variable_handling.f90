!>
!! Routines for handling proxy variables e.g. accumulation buffers
!!
MODULE mo_derived_variable_handling

  USE self_object
  USE self_vector_ref
  USE self_vector
  USE self_map
  USE self_assert
#ifdef DEBUG_MVSTREAM
  USE mo_util_dbg_prnt, ONLY: dbg_print
#endif

  USE mo_kind, ONLY: wp
  USE mo_model_domain, ONLY: t_patch
  USE mo_io_config, ONLY: lnetcdf_flt64_output
  USE mo_dynamics_config, ONLY: nnow, nnew, nold
  USE mo_statistics, ONLY: add_fields, max_fields, min_fields, assign_fields, add_sqr_fields
  USE mo_var_metadata_types, ONLY: VARNAME_LEN
  USE mo_impl_constants, ONLY: vname_len, SUCCESS, max_char_length, TLEV_NNOW, TLEV_NNEW,    &
    &                          TLEV_NNOW_RCF, TLEV_NNEW_RCF, REAL_T
  USE mo_cdi_constants, ONLY: GRID_UNSTRUCTURED_CELL, GRID_UNSTRUCTURED_EDGE, GRID_UNSTRUCTURED_VERT
  USE mo_name_list_output_types, ONLY: t_output_name_list
  USE mo_zaxis_type, ONLY:  ZA_OCEAN_SEDIMENT
  USE mo_mpi, ONLY: my_process_is_stdio
  USE mo_var_list_element, ONLY: level_type_ml, level_type_pl, level_type_hl, level_type_il
  USE mo_name_list_output_metadata, ONLY: metainfo_get_timelevel
  USE mo_var_list, ONLY: new_var_list,&
       total_number_of_variables, &
       get_var_name, default_var_list_settings, add_var, find_element, find_list_element, &
       & get_varname_with_timelevel, delete_var_list, print_all_var_lists
  USE mo_linked_list, ONLY: t_var_list, t_list_element
  USE mo_exception, ONLY: finish, message, message_text
  USE mtime, ONLY: MAX_DATETIME_STR_LEN, newEvent, event, isCurrentEventActive,&
    & newDatetime, datetime, eventToString, divideDatetimeDifferenceInSeconds, &
    & divisionquotientTimespan
  USE mo_util_string, ONLY: int2string
  USE mtime_datetime, ONLY: datetimeToString
  USE mo_output_event_types,  ONLY: t_sim_step_info
  USE mo_time_config,         ONLY: time_config
  USE mo_cdi,                 ONLY: DATATYPE_FLT32, DATATYPE_FLT64, GRID_LONLAT, GRID_ZONAL, TSTEP_CONSTANT

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_derived_variable_handling'

  TYPE(map), SAVE    :: meanMap, meanEvents, meanEventsActivity, meanVarCounter, meanPrognosticPointers
  TYPE(vector), SAVE :: meanPrognostics
  TYPE(map), SAVE    :: maxMap, maxEvents, maxEventsActivity, maxPrognosticPointers, maxVarCounter
  TYPE(vector), SAVE :: maxPrognostics
  TYPE(map), SAVE    :: minMap, minEvents, minEventsActivity, minVarCounter, minPrognosticPointers
  TYPE(vector), SAVE :: minPrognostics
  TYPE(map), SAVE    :: squareMap, squareEvents, squareEventsActivity, squareVarCounter, squarePrognosticPointers
  TYPE(vector), SAVE :: squarePrognostics
  CHARACTER(1), PARAMETER :: separator = achar(124)

  PUBLIC :: init_statistics_streams
  PUBLIC :: finish_statistics_streams
  PUBLIC :: update_statistics
  PUBLIC :: reset_statistics
  PUBLIC :: process_statistics_stream

  TYPE :: t_mvstream_pair
    TYPE(t_list_element), POINTER :: source
    TYPE(t_list_element), POINTER :: destination
    INTEGER                       :: counter = 0
  END TYPE t_mvstream_pair

  CHARACTER(LEN = MAX_CHAR_LENGTH), PARAMETER :: MEAN   = "mean"
  CHARACTER(LEN = MAX_CHAR_LENGTH), PARAMETER :: MAX    = "max"
  CHARACTER(LEN = MAX_CHAR_LENGTH), PARAMETER :: MIN    = "min"
  CHARACTER(LEN = MAX_CHAR_LENGTH), PARAMETER :: SQUARE = "square"
  !>
  ! wrapper for bind_c event type - needed because bind_c types are not allowed
  ! in "select type" - this is kind of a hack
  type :: t_event_wrapper
    type(event), pointer :: this
  end type

CONTAINS

  !>
  !! Create all needed lists and maps
  !!
  SUBROUTINE init_statistics_streams()
    CHARACTER(LEN=max_char_length) :: listname

    ! main map for automatical mean value computation:
    ! key: eventString
    ! value: self-vector of t_list_elements objects that belong to that event
    call meanMap%init(verbose=.false.)
    ! map for holding the events itself
    ! since the events cannot be saved a keys for a map because they are only
    ! C-pointers the same string representation of events is used as for
    ! "meanMap"
    call meanEvents%init(verbose=.false.)
    call meanEventsActivity%init(verbose=.false.)
    call meanVarCounter%init(verbose=.false.)
    call meanPrognostics%init(verbose=.false.)
    call meanPrognosticPointers%init(verbose=.false.)

    call maxMap%init(verbose=.false.)
    call maxEvents%init(verbose=.false.)
    call maxEventsActivity%init(verbose=.false.)
    call maxVarCounter%init(verbose=.false.)
    call maxPrognostics%init(verbose=.false.)
    call maxPrognosticPointers%init(verbose=.false.)

    call minMap%init(verbose=.false.)
    call minEvents%init(verbose=.false.)
    call minEventsActivity%init(verbose=.false.)
    call minVarCounter%init(verbose=.false.)
    call minPrognostics%init(verbose=.false.)
    call minPrognosticPointers%init(verbose=.false.)

    call squareMap%init(verbose=.false.)
    call squareEvents%init(verbose=.false.)
    call squareEventsActivity%init(verbose=.false.)
    call squareVarCounter%init(verbose=.false.)
    call squarePrognostics%init(verbose=.false.)
    call squarePrognosticPointers%init(verbose=.false.)
  END SUBROUTINE init_statistics_streams

  !>
  !! Print contents of a vector, just giving the name for t_list_elements
  !!
  !! Optional label is printed first, on a line by its own
  !!
  SUBROUTINE var_print(this, label)
    class(vector_ref) , INTENT(in) :: this
    CHARACTER(*), INTENT(in), OPTIONAL :: label

    TYPE(vector_iterator) :: my_iter
    CLASS(*), POINTER :: my_buffer
    INTEGER :: i

    IF (PRESENT(label)) PRINT *, label

    my_iter = this%iter()
    DO WHILE(my_iter%next(my_buffer))
      SELECT TYPE(my_buffer)
      TYPE is (t_list_element)
        PRINT *,'t_list_element:varname:',         trim(my_buffer%field%info%name)
      TYPE is (t_mvstream_pair)
        PRINT *,'t_mvstream_pair:source     :',trim(my_buffer%source%field%info%name)
        PRINT *,'t_mvstream_pair:destination:',trim(my_buffer%destination%field%info%name)
      CLASS default
        PRINT *,' default class print  :'
        print *,object_string (my_buffer)
      END SELECT
    END DO
  END SUBROUTINE var_print

  !>
  !! return the event string from given output name list, based on output_start, -end and -interval
  !!
  FUNCTION get_event_key(output_name_list) RESULT(event_key)
    TYPE(t_output_name_list) :: output_name_list
    CHARACTER(LEN=1000) :: event_key
    CHARACTER(LEN=1)    :: separator

    separator = '_'
    event_key = &
      &trim(output_name_list%output_start(1))//separator//&
      &trim(output_name_list%output_end(1))//separator//&
      &trim(output_name_list%output_interval(1))
  END FUNCTION get_event_key

  !>
  !! function to check if a varable is prognostic and saved for stats
  !!
  logical function is_statsPrognosticVariable(variable,statList)
    type(t_list_element), intent(in) :: variable
    type(vector), intent(in) :: statList

    is_statsPrognosticVariable = statList%includes(variable%field%info%name)
  end function

  !>
  !! return pointer to prognostic accumulation copy for given timelevel
  !!
  function get_prognostics_source_pointer(destinationVariable, timelevelIndex, pointerMap) result(sourcePointer)
    type(t_list_element), pointer :: sourcePointer
    type(t_list_element)          :: destinationVariable
    type(map)                     :: pointerMap

    CHARACTER(LEN=*), PARAMETER :: routine =  modname//"::get_prognostics_source_pointer"
    INTEGER , INTENT(IN) :: timelevelIndex

    class(*), pointer :: dummy

    dummy => pointerMap%get( &
      & get_varname_with_timelevel( &
      &   destinationVariable%field%info%name, &
      &   timelevelIndex) &
      & )

    select type (dummy)
    type is (t_list_element)
      sourcePointer => dummy
    class default
      call finish(routine,"Cound not find pointer to prognostics variable for given timelevel!")
    end select
  end function get_prognostics_source_pointer

  !>
  !! The current mean values are supported
  !!   * on the global domain (dom = 1)
  !!   * without stream partitioning
  !! the model should abort under these circumstances
  SUBROUTINE statStreamCrossCheck(p_onl)
    TYPE(t_output_name_list) :: p_onl

    CHARACTER(LEN=*), PARAMETER :: routine =  modname//"::statStreamCrossCheck"
    LOGICAL :: abort

    abort = .FALSE.

    IF (  p_onl%stream_partitions_ml > 1 .OR. &
      &   p_onl%stream_partitions_pl > 1 .OR. &
      &   p_onl%stream_partitions_hl > 1 .OR. &
      &   p_onl%stream_partitions_il > 1 ) abort = .TRUE.

    IF (abort) THEN
      call finish(routine,"meanValues are only supported on global domain 1 and without stream partitioning!")
    END IF
  END SUBROUTINE statStreamCrossCheck

  !>
  !! Delete internal mean value fields
  !!
  SUBROUTINE finish_statistics_streams()

#ifdef DEBUG_MVSTREAM
    IF (my_process_is_stdio()) CALL print_summary('destruct mean stream variables')
#endif

  END SUBROUTINE finish_statistics_streams

  integer function output_varlist_length(in_varlist)
    CHARACTER(LEN=vname_len), intent(in) :: in_varlist(:)

    output_varlist_length = 0
    DO
      IF (in_varlist(output_varlist_length+1) == ' ') EXIT
      output_varlist_length = output_varlist_length + 1
    END DO
  end function output_varlist_length

  SUBROUTINE fill_list_of_output_varnames(varlist,in_varlist,ntotal_vars)
    CHARACTER(LEN=VARNAME_LEN), INTENT(INOUT) :: varlist(:)
    CHARACTER(LEN=vname_len), POINTER, INTENT(IN) :: in_varlist(:)
    INTEGER, INTENT(IN)               :: ntotal_vars

    INTEGER :: output_variables

    ! count variables {{{
    output_variables = 0
    DO
      IF (in_varlist(output_variables+1) == ' ') EXIT
      output_variables = output_variables + 1
    END DO

    IF (output_variables > 0)  varlist(1:output_variables) = in_varlist(1:output_variables)
    varlist((output_variables+1):ntotal_vars) = " "
    ! }}}
  END SUBROUTINE fill_list_of_output_varnames

  SUBROUTINE print_src_dest_info(src_element, dest_element, dest_element_name,in_varlist,i)
    TYPE(t_list_element), POINTER :: src_element, dest_element
    CHARACTER(LEN=VARNAME_LEN) :: dest_element_name
    CHARACTER(LEN=vname_len), POINTER :: in_varlist(:)
    INTEGER,INTENT(IN) :: i

      if (my_process_is_stdio()) CALL print_summary('src(name)     :|'//trim(src_element%field%info%name)//'|',stderr=.true.)
      if (my_process_is_stdio()) CALL print_summary('src(grid)     :|'//int2string(src_element%field%info%hgrid)//'|', &
          & stderr=.true.)
      if (my_process_is_stdio()) CALL print_summary('src(dims:1)   :|'//int2string(src_element%field%info%used_dimensions(1))//'|',&
          & stderr=.true.)
      if (my_process_is_stdio()) CALL print_summary('src(dims:2)   :|'//int2string(src_element%field%info%used_dimensions(2))//'|',&
          & stderr=.true.)
      if (my_process_is_stdio()) CALL print_summary('src(dims:3)   :|'//int2string(src_element%field%info%used_dimensions(3))//'|',&
          & stderr=.true.)
      if (my_process_is_stdio()) CALL print_summary('src(dims:4)   :|'//int2string(src_element%field%info%used_dimensions(4))//'|',&
          & stderr=.true.)
      if (my_process_is_stdio()) CALL print_summary('src(dims:5)   :|'//int2string(src_element%field%info%used_dimensions(5))//'|',&
          & stderr=.true.)
      if (my_process_is_stdio()) CALL print_summary('varlist(name) :|'//trim(in_varlist(i))//'|',stderr=.true.)
      if (my_process_is_stdio()) CALL print_summary('new name      :|'//trim(dest_element_name)//'|',stderr=.true.)
      if (my_process_is_stdio()) CALL print_summary('new grid      :|'//int2string(dest_element%field%info%hgrid)//'|',&
          & stderr=.true.)
  END SUBROUTINE

  SUBROUTINE find_src_element(src_element, varlist_element, dest_element_name, &
          &  dom, src_list, prognosticsList, prognosticsPointerList)
    TYPE(t_list_element), POINTER :: src_element
    CHARACTER(LEN=VARNAME_LEN), INTENT(IN) :: varlist_element
    CHARACTER(LEN=VARNAME_LEN) :: dest_element_name
    INTEGER, INTENT(IN) :: dom
    TYPE(t_var_list), POINTER :: src_list
    TYPE(vector) :: prognosticsList
    TYPE(map)    :: prognosticsPointerList

    LOGICAL :: foundPrognostic
    INTEGER :: timelevel, timelevels(3)
    CHARACTER(LEN=*), PARAMETER :: routine =  modname//"::find_src_element"

    ! find existing source variable on all possible ICON grids with the identical name {{{
    src_element => find_element ( TRIM(varlist_element), &
        &                         opt_patch_id=dom, &
        &                         opt_hgrid=GRID_UNSTRUCTURED_CELL, &
        &                         opt_caseInsensitive=.true., &
        &                         opt_returnList=src_list)
    IF (.not. ASSOCIATED(src_element) ) &
        & src_element => find_element ( TRIM(varlist_element),&
        &                               opt_patch_id=dom, &
        &                               opt_hgrid=GRID_UNSTRUCTURED_EDGE, &
        &                               opt_caseInsensitive=.true., &
        &                               opt_returnList=src_list)
    IF (.not. ASSOCIATED(src_element) ) &
        & src_element => find_element ( TRIM(varlist_element),&
        &                               opt_patch_id=dom, &
        &                               opt_hgrid=GRID_UNSTRUCTURED_VERT, &
        &                               opt_caseInsensitive=.true., &
        &                               opt_returnList=src_list)
    IF (.not. ASSOCIATED(src_element) ) &
        & src_element => find_element ( TRIM(varlist_element),&
        &                               opt_patch_id=dom, &
        &                               opt_hgrid=GRID_LONLAT, &
        &                               opt_caseInsensitive=.true., &
        &                               opt_returnList=src_list)
    IF (.not. ASSOCIATED(src_element) ) &
        & src_element => find_element ( TRIM(varlist_element),&
        &                               opt_patch_id=dom, &
        &                               opt_hgrid=GRID_ZONAL, &
        &                               opt_caseInsensitive=.true., &
        &                               opt_returnList=src_list)
    ! }}}

    ! if not found: maybe it is a prognostic variable, so it has the
    ! time-level in its name
    ! ATTENTION: this is only a placeholder, because it catches the first match
    ! the correct pointer is the one with nnew(), but its value changes each timestep
    IF (.not. ASSOCIATED (src_element)) THEN
      foundPrognostic = .false.
      timelevels = (/nold(1),nnow(1),nnew(1)/)
      DO timelevel=1,3

#ifdef DEBUG_MVSTREAM
          if (my_process_is_stdio()) call print_error(get_varname_with_timelevel(varlist_element,timelevels(timelevel)))
#endif
          src_element => find_element(get_varname_with_timelevel(varlist_element,timelevels(timelevel)),&
              &                       opt_patch_id=dom, opt_returnList=src_list)
          if ( ASSOCIATED(src_element) ) then
#ifdef DEBUG_MVSTREAM
            if (my_process_is_stdio()) write(0,*)'found prognostic:',&
                & TRIM(get_varname_with_timelevel(varlist_element,timelevels(timelevel)))
#endif
            if ( .not. foundPrognostic ) then
              ! save the name of the original output variable if a prognosting version was found
              call prognosticsList%add(dest_element_name)
              foundPrognostic = .true.
#ifdef DEBUG_MVSTREAM
              if (my_process_is_stdio()) call print_error('prognosticsList%add():'//varlist_element)
#endif
            end if
            ! save the the pointers for all time levels of a prognostic variable
            ! these must be used for correct accumulation during the time loop
            call prognosticsPointerList%add(get_varname_with_timelevel(dest_element_name,timelevels(timelevel)),src_element)
          end if
      END DO
    END IF
    ! in case nothing appropriate could be found, throw an error
    IF (.not. ASSOCIATED (src_element)) THEN
      call finish(routine,'Could not find source variable:'//TRIM(varlist_element))
    END IF
  END SUBROUTINE

  SUBROUTINE setup_statistics_events_and_lists(statsMap, eventTag, statisticsVariablesForEvent, eventMap, eventsActivityMap,p_onl)
    TYPE(map) :: statsMap, eventMap, eventsActivityMap
    CHARACTER(LEN=100), intent(in) :: eventTag
    TYPE(vector_ref) :: statisticsVariablesForEvent
    TYPE (t_output_name_list), intent(in),target  :: p_onl

    class(*), pointer :: myBuffer
    type(t_event_wrapper) :: event_wrapper

    IF ( statsMap%has_key(eventTag) ) THEN
      myBuffer => statsMap%get(eventTag)
      select type (myBuffer)
      type is (vector_ref)
        ! use exiting list
        statisticsVariablesForEvent = myBuffer
      end select
    ELSE
      ! create new variable list
      call statisticsVariablesForEvent%init()

      ! wrap c-pointer for events into a fortran type - look weird, but does the trick
      event_wrapper = t_event_wrapper(this=newEvent(eventTag, &
        &                 p_onl%output_start(1), &
        &                 p_onl%output_start(1), &
        &                 p_onl%output_end(1), &
        &                 p_onl%output_interval(1) &
        &                ))

      ! collect the new event
      call eventMap%add(eventTag,event_wrapper)
    END IF
    ! initialize all events as in-active
    call eventsActivityMap%add(eventTag,.false.)
  END SUBROUTINE

  SUBROUTINE print_prognosticLists(varlist,progMap,progPointers)
    TYPE(vector_ref) :: varlist
    TYPE(map) :: progPointers, progMap

    CHARACTER(LEN=*), PARAMETER :: routine =  modname//"::print_prognosticLists"

    call var_print(varlist)
    call print_routine(routine, 'progMap%to_string()')
    call print_routine(routine, progMap%to_string())
    call print_routine(routine, progPointers%to_string())
    call var_print(progPointers%values())
  END SUBROUTINE

  !>
  !! copy varlist element (source_element) to target varlist (list) with given name
  !! lrestart = .true.
  !!
  FUNCTION copy_var_to_list(list,name,source_element,patch_2d) RESULT(dest_element)
    TYPE(t_var_list) :: list
    CHARACTER(LEN=VARNAME_LEN) :: name
    TYPE(t_list_element),POINTER :: source_element
    TYPE(t_patch),TARGET, INTENT(in)    :: patch_2d

    TYPE(t_list_element), POINTER :: dest_element

    integer :: dataType = -1
    CHARACTER(LEN=*), PARAMETER :: routine =  modname//"::copy_var_to_list"

    dataType = MERGE(DATATYPE_FLT64, DATATYPE_FLT32, lnetcdf_flt64_output)
#ifdef DEBUG_MVSTREAM
    if (my_process_is_stdio()) then
      call print_verbose("COPY variable:"//TRIM(name))
      call print_verbose("     INTO    :"//TRIM(list%p%name))
    endif
#endif
    CALL add_var(source_element%field%info%ndims, &
      & REAL_T, &
      & list, name, &
      & source_element%field%info%hgrid, source_element%field%info%vgrid, &
      & source_element%field%info%cf, source_element%field%info%grib2,  &
      & source_element%field%info%used_dimensions, &
      & dest_element, &
      & tlev_source=source_element%field%info%tlev_source, &
      & isteptype=source_element%field%info%isteptype, &
      & post_op=source_element%field%info%post_op, &
      & initval_r=source_element%field%info%resetval%rval, &
      & resetval_r=source_element%field%info%resetval%rval, &
      & lmiss=source_element%field%info%lmiss, &
      & missval_r=source_element%field%info%missval%rval, &
      & action_list=source_element%field%info%action_list, &
      & vert_interp=source_element%field%info%vert_interp, &
      & hor_interp=source_element%field%info%hor_interp, &
      & in_group=source_element%field%info%in_group, &
      & l_pp_scheduler_task=source_element%field%info%l_pp_scheduler_task, &
      & loutput=.TRUE., lrestart=.FALSE., &
      & var_class=source_element%field%info%var_class )

    ! add the subset for later accumulation on all types of horizontal grids
    select case (source_element%field%info%hgrid)
    case (GRID_UNSTRUCTURED_CELL)
      dest_element%field%info%subset = patch_2d%cells%owned
    case (GRID_UNSTRUCTURED_EDGE)
      dest_element%field%info%subset = patch_2d%edges%owned
    case (GRID_UNSTRUCTURED_VERT)
      dest_element%field%info%subset = patch_2d%verts%owned
    end select
  END FUNCTION copy_var_to_list

  !>
  !! return internal name for accumulation variables
  !!
  FUNCTION get_statistics_varname(varname,output_setup)
    CHARACTER(LEN=VARNAME_LEN)  :: varname
    type(t_output_name_list) :: output_setup

    CHARACTER(LEN=VARNAME_LEN)  :: get_statistics_varname

    get_statistics_varname = &
      &TRIM(varname)//separator//&
      &TRIM(output_setup%operation)//separator//&
      &TRIM(output_setup%output_interval(1))//separator//&
      &TRIM(output_setup%output_start(1))//separator//&
      &'DOM'//TRIM(int2string(output_setup%dom))

  END FUNCTION get_statistics_varname

  !>
  !! return internal name from accumulation variable name
  !!
  FUNCTION get_real_varname(mean_varname)
    CHARACTER(LEN=VARNAME_LEN)  :: mean_varname

    CHARACTER(LEN=VARNAME_LEN)  :: get_real_varname

    get_real_varname = mean_varname(1:INDEX(mean_varname,separator)-1)
  END FUNCTION get_real_varname
 
  SUBROUTINE process_mvstream(p_onl,i_typ, sim_step_info, patch_2d, &
           & statisticMap, &
           & statisticEvents, &
           & statisticEventsActivity, &
           & statisticVarCounter, &
           & statisticPrognostics, &
           & statisticPrognosticPointers, &
           & operation)
    TYPE (t_output_name_list), target  :: p_onl
    INTEGER                            :: i_typ
    TYPE (t_sim_step_info), INTENT(IN) :: sim_step_info
    TYPE(t_patch), INTENT(IN)          :: patch_2d

    type(map) :: statisticMap
    type(map) :: statisticEvents
    type(map) :: statisticEventsActivity
    type(map) :: statisticVarCounter
    type(vector) :: statisticPrognostics
    type(map) :: statisticPrognosticPointers
    CHARACTER(len=*), intent(in) :: operation

    !-----------------------------------------------------------------------------
    CHARACTER(LEN=vname_len), POINTER :: in_varlist(:)
    INTEGER :: ntotal_vars, output_variables,i,ierrstat, dataType
    INTEGER :: timelevel, timelevels(3)
    type(vector_ref) :: statisticVariables, prognosticVariables
    CHARACTER(LEN=100) :: eventKey
    type(t_event_wrapper) :: event_wrapper
    class(*), pointer :: myBuffer
    TYPE(t_list_element), POINTER :: src_element, dest_element
    CHARACTER(LEN=VARNAME_LEN), ALLOCATABLE :: varlist(:)
    CHARACTER(LEN=VARNAME_LEN) :: dest_element_name
    LOGICAL :: foundPrognostic
    TYPE(t_var_list), POINTER :: src_list
    CHARACTER(LEN=*), PARAMETER :: routine =  modname//"::process_mvstream"

#ifdef DEBUG_MVSTREAM
    if (my_process_is_stdio()) call print_routine(routine,'start')
#endif

    IF (trim(operation) .EQ. TRIM(p_onl%operation)) THEN

#ifdef DEBUG_MVSTREAM
      if (my_process_is_stdio()) call print_routine(routine,'found "'//trim(operation)//'" operation')
#endif
      call statStreamCrossCheck(p_onl)

      ntotal_vars = total_number_of_variables()
      ! temporary variables needed for variable group parsing
      ALLOCATE(varlist(ntotal_vars), STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')

      IF (i_typ == level_type_ml) in_varlist => p_onl%ml_varlist
      IF (i_typ == level_type_pl) in_varlist => p_onl%pl_varlist
      IF (i_typ == level_type_hl) in_varlist => p_onl%hl_varlist
      IF (i_typ == level_type_il) in_varlist => p_onl%il_varlist

      ! count variables {{{
      output_variables = 0
      DO
        IF (in_varlist(output_variables+1) == ' ') EXIT
        output_variables = output_variables + 1
      END DO

      IF (output_variables > 0)  varlist(1:output_variables) = in_varlist(1:output_variables)
      varlist((output_variables+1):ntotal_vars) = " "
      ! }}}

      ! uniq identifier for an event based on output start/end/interval
      eventKey = get_event_key(p_onl)
      ! this has the advantage that we can compute a uniq id without creating
      ! the event itself

#ifdef DEBUG_MVSTREAM
      if (my_process_is_stdio()) call print_summary('eventKey:'//trim(eventKey))
#endif

      ! fill main dictionary of variables for different events
       CALL setup_statistics_events_and_lists(statisticMap, &
                                           & eventKey, &
                                           & statisticVariables, &
                                           & statisticEvents, &
                                           & statisticEventsActivity, &
                                           & p_onl)

      ! create adhoc copies of all variables for later accumulation
      !
      ! varlist MUST by a list of add_var-identifiers, i.e. the given name
      DO i=1, output_variables
        ! collect data variables only
        ! variables names like 'grid:clon' which should be excluded
        IF ( INDEX(varlist(i),':') > 0) CYCLE

        ! check for already created meanStream variable (maybe from another output_nml with the same output_interval)
        ! names consist of original spot-value names PLUS event information (start + interval of output)
        ! TODO: unify with eventKey definition if possible
        dest_element_name = get_statistics_varname(varlist(i),p_onl)
        dest_element => find_element(trim(dest_element_name),opt_patch_id=p_onl%dom)
#ifdef DEBUG_MVSTREAM
        if (my_process_is_stdio()) call print_summary('destination variable NAME:'//TRIM(dest_element_name),stderr=.true.)
#endif
        IF (.not. ASSOCIATED(dest_element) ) THEN !not found -->> create a new on
          ! find existing source variable on all possible ICON grids with the identical name
          CALL find_src_element(src_element, varlist(i), dest_element_name, p_onl%dom, src_list, &
            & statisticPrognostics, statisticPrognosticPointers)
          ! avoid mean processing for instantaneous fields
          IF (TSTEP_CONSTANT .eq. src_element%field%info%isteptype) CYCLE

          ! add new mean variable, copy the meta-data from the existing variable
          ! 1. copy the source variable to destination pointer
          dest_element => copy_var_to_list(src_list,dest_element_name,src_element, patch_2d)

#ifdef DEBUG_MVSTREAM
          CALL print_src_dest_info(src_element, dest_element, dest_element_name, in_varlist, i)
#endif

          ! set output to double precission if necessary
          dest_element%field%info%cf%datatype = MERGE(DATATYPE_FLT64, DATATYPE_FLT32, lnetcdf_flt64_output)

          ! 2. update the nc-shortname to internal name of the source variable unless it is already set by the user
          IF("" == dest_element%field%info%cf%short_name) dest_element%field%info%cf%short_name = get_var_name(src_element%field)

          ! Collect variable pointers for source and destination in the same list {{{
          CALL statisticVariables%add(src_element)
          CALL statisticVariables%add(dest_element)
          ! collect the counter for each destination in a separate list
          CALL statisticVarCounter%add(dest_element%field%info%name,0)

        ! replace existince varname in output_nml with the meanStream Variable
#ifdef DEBUG_MVSTREAM
          if ( my_process_is_stdio()) CALL print_summary('dst(name)     :|'//trim(dest_element%field%info%name)//'|',&
              & stderr=.true.)
          if ( my_process_is_stdio()) CALL print_summary('dst(shortname):|'//trim(dest_element%field%info%cf%short_name)//'|',&
              & stderr=.true.)
#endif
        END IF
        in_varlist(i) = trim(dest_element%field%info%name)
      END DO
      CALL statisticMap%add(eventKey,statisticVariables)
#ifdef DEBUG_MVSTREAM
      if (my_process_is_stdio()) call print_error(routine//": statisticPrognostics%to_string()")
      if (my_process_is_stdio()) call print_error(statisticPrognostics%to_string())
#endif
      DEALLOCATE(varlist)
    ELSE
#ifdef DEBUG_MVSTREAM
      if (my_process_is_stdio()) call print_routine(routine,'NO "statistic" operation found')
#endif
    END IF

#ifdef DEBUG_MVSTREAM
    if (my_process_is_stdio()) then
      call print_prognosticLists(statisticVariables, statisticMap, statisticPrognosticPointers)
    endif
    if (my_process_is_stdio()) call print_routine(routine,'end',stderr=.true.)
#endif

  END SUBROUTINE process_mvstream
  SUBROUTINE process_statistics_stream(p_onl, i_typ, sim_step_info, patch_2d)
    TYPE (t_output_name_list), target  :: p_onl
    INTEGER                            :: i_typ
    TYPE (t_sim_step_info), INTENT(IN) :: sim_step_info
    TYPE(t_patch), INTENT(IN)          :: patch_2d

    CALL process_mvstream(p_onl,i_typ,sim_step_info, patch_2d, &
                        & meanMap, &
                        & meanEvents, &
                        & meanEventsActivity, &
                        & meanVarCounter, &
                        & meanPrognostics, &
                        & meanPrognosticPointers, &
                        & MEAN)
    CALL process_mvstream(p_onl,i_typ,sim_step_info, patch_2d, &
                        & maxMap, &
                        & maxEvents, &
                        & maxEventsActivity, &
                        & maxVarCounter, &
                        & maxPrognostics, &
                        & maxPrognosticPointers, &
                        & MAX)
    CALL process_mvstream(p_onl,i_typ,sim_step_info, patch_2d, &
                        & minMap, &
                        & minEvents, &
                        & minEventsActivity, &
                        & minVarCounter, &
                        & minPrognostics, &
                        & minPrognosticPointers, &
                        & MIN)
    CALL process_mvstream(p_onl,i_typ,sim_step_info, patch_2d, &
                        & squareMap, &
                        & squareEvents, &
                        & squareEventsActivity, &
                        & squareVarCounter, &
                        & squarePrognostics, &
                        & squarePrognosticPointers, &
                        & SQUARE)
  END SUBROUTINE process_statistics_stream

  ! methods needed for update statistics each timestep
  !>
  !! implement addition for source fields to the internal accumulation fields
  !!
  SUBROUTINE statistics_add(source, destination, counter)
    type(t_list_element) , INTENT(IN)    :: source
    type(t_list_element) , INTENT(INOUT) :: destination
    integer, intent(inout)               :: counter

    integer :: index

    index = -1

    IF (.NOT. ASSOCIATED(source%field%r_ptr)) THEN
      CALL add_sp()
    ELSE
      CALL add_wp()
    ENDIF
    counter                 = counter + 1

    CONTAINS
      SUBROUTINE add_wp()
        IF (source%field%info%lcontained) THEN
          index = source%field%info%ncontained

          ! tackle 4d containers
          SELECT CASE(source%field%info%var_ref_pos)
          CASE(1)
            destination%field%r_ptr(:,:,:,:,:) = &
              & destination%field%r_ptr (:,:,:,:,:) + source%field%r_ptr(index:index,:,:,:,:)
          CASE(2)
            destination%field%r_ptr(:,:,:,:,:) = &
              & destination%field%r_ptr (:,:,:,:,:) + source%field%r_ptr(:,index:index,:,:,:)
          CASE(3)
            destination%field%r_ptr(:,:,:,:,:) = &
              & destination%field%r_ptr (:,:,:,:,:) + source%field%r_ptr(:,:,index:index,:,:)
          CASE(4)
            destination%field%r_ptr(:,:,:,:,:) = &
              & destination%field%r_ptr (:,:,:,:,:) + source%field%r_ptr(:,:,:,index:index,:)
          END SELECT
        ELSE
          ! tackle 1d, 2d and 3d
          SELECT CASE(destination%field%info%ndims)
          CASE(3)
           !hack for sea ice variables which uses vertical level for ice class
            IF (1 == destination%field%info%used_dimensions(2)) THEN
              call add_fields(destination%field%r_ptr(:,:,:,1,1), &
                &             source%field%r_ptr(:,:,:,1,1), &
                &             destination%field%info%subset, &
                &             levels=1, &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            ELSE IF (ZA_OCEAN_SEDIMENT == destination%field%info%vgrid) THEN
              call add_fields(destination%field%r_ptr(:,:,:,1,1), &
                &             source%field%r_ptr(:,:,:,1,1), &
                &             destination%field%info%subset, &
                &             levels=destination%field%info%used_dimensions(2), &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            ELSE
              call add_fields(destination%field%r_ptr(:,:,:,1,1), &
                &             source%field%r_ptr(:,:,:,1,1), &
                &             destination%field%info%subset, &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            ENDIF
          CASE(2)
            IF (GRID_ZONAL .EQ. destination%field%info%hgrid) THEN
              call add_fields(destination%field%r_ptr(:,:,1,1,1), &
                &             source%field%r_ptr(:,:,1,1,1), &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            ELSE
              call add_fields(destination%field%r_ptr(:,:,1,1,1), &
                &             source%field%r_ptr(:,:,1,1,1), &
                &             destination%field%info%subset, &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            ENDIF
          CASE DEFAULT
            destination%field%r_ptr(:,:,:,:,:) = destination%field%r_ptr (:,:,:,:,:) + source%field%r_ptr(:,:,:,:,:)
          END SELECT
        ENDIF
      END SUBROUTINE add_wp
      SUBROUTINE add_sp() 
        IF (source%field%info%lcontained) then
          index = source%field%info%ncontained
          ! tackle 4d containers
          if (my_process_is_stdio()) print *,'sheeep:',shape(source%field%s_ptr)
          if (my_process_is_stdio()) print *,'deeeep:',shape(destination%field%s_ptr)
          if (my_process_is_stdio()) print *,'index:',trim(int2string(index))
          SELECT CASE(source%field%info%var_ref_pos)
          CASE(1)
            destination%field%r_ptr(:,:,:,:,:) = &
              & destination%field%r_ptr (:,:,:,:,:) + REAL(source%field%r_ptr(index:index,:,:,:,:),wp)
          CASE(2)
            destination%field%r_ptr(:,:,:,:,:) = &
              & destination%field%r_ptr (:,:,:,:,:) + REAL(source%field%r_ptr(:,index:index,:,:,:),wp)
          CASE(3)
            destination%field%r_ptr(:,:,:,:,:) = &
              & destination%field%r_ptr (:,:,:,:,:) + REAL(source%field%r_ptr(:,:,index:index,:,:) ,wp)
          CASE(4)
            destination%field%r_ptr(:,:,:,:,:) = &
              & destination%field%r_ptr (:,:,:,:,:) + REAL(source%field%r_ptr(:,:,:,index:index,:),wp)
          END SELECT
        ELSE
          ! tackle 1d, 2d and 3d
          SELECT CASE(destination%field%info%ndims)
          CASE(3)
           !hack for sea ice variables which uses vertical level for ice class
            IF (1 == destination%field%info%used_dimensions(2)) THEN
              call add_fields(destination%field%r_ptr(:,:,:,1,1), &
                &             source%field%s_ptr(:,:,:,1,1), &
                &             destination%field%info%subset, &
                &             levels=1, &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            ELSE IF (ZA_OCEAN_SEDIMENT == destination%field%info%vgrid) THEN
              call add_fields(destination%field%r_ptr(:,:,:,1,1), &
                &             source%field%s_ptr(:,:,:,1,1), &
                &             destination%field%info%subset, &
                &             levels=destination%field%info%used_dimensions(2), &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            ELSE
              call add_fields(destination%field%r_ptr(:,:,:,1,1), &
                &             source%field%s_ptr(:,:,:,1,1), &
                &             destination%field%info%subset, &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            ENDIF
          CASE(2)
            IF (GRID_ZONAL .EQ. destination%field%info%hgrid) THEN
              call add_fields(destination%field%r_ptr(:,:,1,1,1), &
                &             source%field%s_ptr(:,:,1,1,1), &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            ELSE
              call add_fields(destination%field%r_ptr(:,:,1,1,1), &
                &             source%field%s_ptr(:,:,1,1,1), &
                &             destination%field%info%subset, &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            ENDIF
          CASE DEFAULT
            destination%field%r_ptr(:,:,:,:,:) = destination%field%r_ptr (:,:,:,:,:) + REAL(source%field%s_ptr(:,:,:,:,:),wp)
          END SELECT
        ENDIF
      END SUBROUTINE add_sp
  END SUBROUTINE statistics_add

  ! methods needed for update statistics each timestep
  !>
  !! implement addition for source fields to the internal accumulation fields
  !!
  SUBROUTINE statistics_add_sqr(source, destination, counter)
    type(t_list_element) , INTENT(IN)    :: source
    type(t_list_element) , INTENT(INOUT) :: destination
    integer, intent(inout)               :: counter

    integer :: index

    index = -1

    IF (.NOT. ASSOCIATED(source%field%r_ptr)) THEN
      CALL add_sqr_sp()
    ELSE
      CALL add_sqr_wp()
    ENDIF
    counter                 = counter + 1

    CONTAINS
      SUBROUTINE add_sqr_wp()
        IF (source%field%info%lcontained) THEN
          index = source%field%info%ncontained

          ! tackle 4d containers
          SELECT CASE(source%field%info%var_ref_pos)
          CASE(1)
            destination%field%r_ptr(:,:,:,:,:) = &
              & destination%field%r_ptr (:,:,:,:,:) + source%field%r_ptr(index:index,:,:,:,:)**2
          CASE(2)
            destination%field%r_ptr(:,:,:,:,:) = &
              & destination%field%r_ptr (:,:,:,:,:) + source%field%r_ptr(:,index:index,:,:,:)**2
          CASE(3)
            destination%field%r_ptr(:,:,:,:,:) = &
              & destination%field%r_ptr (:,:,:,:,:) + source%field%r_ptr(:,:,index:index,:,:)**2
          CASE(4)
            destination%field%r_ptr(:,:,:,:,:) = &
              & destination%field%r_ptr (:,:,:,:,:) + source%field%r_ptr(:,:,:,index:index,:)**2
          END SELECT
        ELSE
          ! tackle 1d, 2d and 3d
          SELECT CASE(destination%field%info%ndims)
          CASE(3)
           !hack for sea ice variables which uses vertical level for ice class
            IF (1 == destination%field%info%used_dimensions(2)) THEN
              call add_sqr_fields(destination%field%r_ptr(:,:,:,1,1), &
                &             source%field%r_ptr(:,:,:,1,1), &
                &             destination%field%info%subset, &
                &             levels=1, &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            ELSE IF (ZA_OCEAN_SEDIMENT == destination%field%info%vgrid) THEN
              call add_sqr_fields(destination%field%r_ptr(:,:,:,1,1), &
                &             source%field%r_ptr(:,:,:,1,1), &
                &             destination%field%info%subset, &
                &             levels=destination%field%info%used_dimensions(2), &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            ELSE
              call add_sqr_fields(destination%field%r_ptr(:,:,:,1,1), &
                &             source%field%r_ptr(:,:,:,1,1), &
                &             destination%field%info%subset, &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            ENDIF
          CASE(2)
            IF (GRID_ZONAL .EQ. destination%field%info%hgrid) THEN
              call add_sqr_fields(destination%field%r_ptr(:,:,1,1,1), &
                &             source%field%r_ptr(:,:,1,1,1), &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            ELSE
              call add_sqr_fields(destination%field%r_ptr(:,:,1,1,1), &
                &             source%field%r_ptr(:,:,1,1,1), &
                &             destination%field%info%subset, &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            ENDIF
          CASE DEFAULT
            destination%field%r_ptr(:,:,:,:,:) = destination%field%r_ptr (:,:,:,:,:) + source%field%r_ptr(:,:,:,:,:)**2
          END SELECT
        ENDIF
      END SUBROUTINE add_sqr_wp
      SUBROUTINE add_sqr_sp()
        IF (source%field%info%lcontained) THEN
          index = source%field%info%ncontained

          ! tackle 4d containers
          SELECT CASE(source%field%info%var_ref_pos)
          CASE(1)
            destination%field%r_ptr(:,:,:,:,:) = &
              & destination%field%r_ptr (:,:,:,:,:) + REAL(source%field%r_ptr(index:index,:,:,:,:),wp)**2
          CASE(2)
            destination%field%r_ptr(:,:,:,:,:) = &
              & destination%field%r_ptr (:,:,:,:,:) + REAL(source%field%r_ptr(:,index:index,:,:,:),wp)**2
          CASE(3)
            destination%field%r_ptr(:,:,:,:,:) = &
              & destination%field%r_ptr (:,:,:,:,:) + REAL(source%field%r_ptr(:,:,index:index,:,:),wp)**2
          CASE(4)
            destination%field%r_ptr(:,:,:,:,:) = &
              & destination%field%r_ptr (:,:,:,:,:) + REAL(source%field%r_ptr(:,:,:,index:index,:),wp)**2
          END SELECT
        ELSE
          ! tackle 1d, 2d and 3d
          SELECT CASE(destination%field%info%ndims)
          CASE(3)
           !hack for sea ice variables which uses vertical level for ice class
            IF (1 == destination%field%info%used_dimensions(2)) THEN
              call add_sqr_fields(destination%field%r_ptr(:,:,:,1,1), &
                &             source%field%s_ptr(:,:,:,1,1), &
                &             destination%field%info%subset, &
                &             levels=1, &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            ELSE IF (ZA_OCEAN_SEDIMENT == destination%field%info%vgrid) THEN
              call add_sqr_fields(destination%field%r_ptr(:,:,:,1,1), &
                &             source%field%s_ptr(:,:,:,1,1), &
                &             destination%field%info%subset, &
                &             levels=destination%field%info%used_dimensions(2), &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            ELSE
              call add_sqr_fields(destination%field%r_ptr(:,:,:,1,1), &
                &             source%field%s_ptr(:,:,:,1,1), &
                &             destination%field%info%subset, &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            ENDIF
          CASE(2)
            IF (GRID_ZONAL .EQ. destination%field%info%hgrid) THEN
              call add_sqr_fields(destination%field%r_ptr(:,:,1,1,1), &
                &             source%field%s_ptr(:,:,1,1,1), &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            ELSE
              call add_sqr_fields(destination%field%r_ptr(:,:,1,1,1), &
                &             source%field%s_ptr(:,:,1,1,1), &
                &             destination%field%info%subset, &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            ENDIF
          CASE DEFAULT
            destination%field%r_ptr(:,:,:,:,:) = destination%field%r_ptr (:,:,:,:,:) + source%field%r_ptr(:,:,:,:,:)**2
          END SELECT
        ENDIF
      END SUBROUTINE add_sqr_sp
  END SUBROUTINE statistics_add_sqr

  SUBROUTINE statistics_assign(source, destination, counter)
    type(t_list_element) , INTENT(IN)    :: source
    type(t_list_element) , INTENT(INOUT) :: destination
    integer, intent(inout)               :: counter

    integer :: index

    index = -1

    IF (.NOT. ASSOCIATED(source%field%r_ptr)) THEN
      CALL assign_sp()
    ELSE
      CALL assign_wp()
    ENDIF
    counter                 = counter + 1

    CONTAINS
      SUBROUTINE assign_wp()
        IF (source%field%info%lcontained) THEN
          index = source%field%info%ncontained

          ! tackle 4d containers
          SELECT CASE(source%field%info%var_ref_pos)
          CASE(1)
            destination%field%r_ptr(:,:,:,:,:) = source%field%r_ptr(index:index,:,:,:,:)
          CASE(2)
            destination%field%r_ptr(:,:,:,:,:) = source%field%r_ptr(:,index:index,:,:,:)
          CASE(3)
            destination%field%r_ptr(:,:,:,:,:) = source%field%r_ptr(:,:,index:index,:,:)
          CASE(4)
            destination%field%r_ptr(:,:,:,:,:) = source%field%r_ptr(:,:,:,index:index,:)
          END SELECT
        ELSE
          ! tackle 1d, 2d and 3d
          SELECT CASE(destination%field%info%ndims)
          CASE(3)
           !hack for sea ice variables which uses vertical level for ice class
            IF (1 == destination%field%info%used_dimensions(2)) THEN
              call assign_fields(destination%field%r_ptr(:,:,:,1,1), &
                &             source%field%r_ptr(:,:,:,1,1), &
                &             destination%field%info%subset, &
                &             levels=1, &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            ELSE IF (ZA_OCEAN_SEDIMENT == destination%field%info%vgrid) THEN
              call assign_fields(destination%field%r_ptr(:,:,:,1,1), &
                &             source%field%r_ptr(:,:,:,1,1), &
                &             destination%field%info%subset, &
                &             levels=destination%field%info%used_dimensions(2), &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            ELSE
              call assign_fields(destination%field%r_ptr(:,:,:,1,1), &
                &             source%field%r_ptr(:,:,:,1,1), &
                &             destination%field%info%subset, &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            ENDIF
          CASE(2)
            IF (GRID_ZONAL .EQ. destination%field%info%hgrid) THEN
              call assign_fields(destination%field%r_ptr(:,:,1,1,1), &
                &             source%field%r_ptr(:,:,1,1,1), &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            ELSE
              call assign_fields(destination%field%r_ptr(:,:,1,1,1), &
                &             source%field%r_ptr(:,:,1,1,1), &
                &             destination%field%info%subset, &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            ENDIF
          CASE DEFAULT
            destination%field%r_ptr(:,:,:,:,:) = source%field%r_ptr(:,:,:,:,:)
          END SELECT
        ENDIF
      END SUBROUTINE assign_wp
      SUBROUTINE assign_sp()
        IF (source%field%info%lcontained) THEN
          index = source%field%info%ncontained

          ! tackle 4d containers
          SELECT CASE(source%field%info%var_ref_pos)
          CASE(1)
            destination%field%r_ptr(:,:,:,:,:) = REAL(source%field%r_ptr(index:index,:,:,:,:),wp)
          CASE(2)
            destination%field%r_ptr(:,:,:,:,:) = REAL(source%field%r_ptr(:,index:index,:,:,:),wp)
          CASE(3)
            destination%field%r_ptr(:,:,:,:,:) = REAL(source%field%r_ptr(:,:,index:index,:,:),wp)
          CASE(4)
            destination%field%r_ptr(:,:,:,:,:) = REAL(source%field%r_ptr(:,:,:,index:index,:),wp)
          END SELECT
        ELSE
          ! tackle 1d, 2d and 3d
          SELECT CASE(destination%field%info%ndims)
          CASE(3)
           !hack for sea ice variables which uses vertical level for ice class
            IF (1 == destination%field%info%used_dimensions(2)) THEN
              call assign_fields(destination%field%r_ptr(:,:,:,1,1), &
                &             source%field%s_ptr(:,:,:,1,1), &
                &             destination%field%info%subset, &
                &             levels=1, &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            ELSE IF (ZA_OCEAN_SEDIMENT == destination%field%info%vgrid) THEN
              call assign_fields(destination%field%r_ptr(:,:,:,1,1), &
                &             source%field%s_ptr(:,:,:,1,1), &
                &             destination%field%info%subset, &
                &             levels=destination%field%info%used_dimensions(2), &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            ELSE
              call assign_fields(destination%field%r_ptr(:,:,:,1,1), &
                &             source%field%s_ptr(:,:,:,1,1), &
                &             destination%field%info%subset, &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            ENDIF
          CASE(2)
            IF (GRID_ZONAL .EQ. destination%field%info%hgrid) THEN
              call assign_fields(destination%field%r_ptr(:,:,1,1,1), &
                &             source%field%s_ptr(:,:,1,1,1), &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            ELSE
              call assign_fields(destination%field%r_ptr(:,:,1,1,1), &
                &             source%field%s_ptr(:,:,1,1,1), &
                &             destination%field%info%subset, &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            ENDIF
          CASE DEFAULT
            destination%field%r_ptr(:,:,:,:,:) = REAL(source%field%r_ptr(:,:,:,:,:),wp)
          END SELECT
        ENDIF
      END SUBROUTINE assign_sp
  END SUBROUTINE statistics_assign

  SUBROUTINE statistics_max(source, destination, counter)
    type(t_list_element) , INTENT(IN)    :: source
    type(t_list_element) , INTENT(INOUT) :: destination
    integer, intent(inout)               :: counter

    integer :: index

    index = -1

    IF (.NOT. ASSOCIATED(source%field%r_ptr)) THEN
      CALL max_sp()
    ELSE
      CALL max_wp()
    ENDIF
    counter                 = counter + 1

    CONTAINS
      SUBROUTINE max_wp()
        IF (source%field%info%lcontained) THEN
          index = source%field%info%ncontained

          ! tackle 4d containers
          SELECT CASE(source%field%info%var_ref_pos)
          CASE(1)
            destination%field%r_ptr(:,:,:,:,:) = &
            & merge(destination%field%r_ptr (:,:,:,:,:), &
            &       source%field%r_ptr(index:index,:,:,:,:), &
            &       destination%field%r_ptr(:,:,:,:,:).gt.source%field%r_ptr(index:index,:,:,:,:))
          CASE(2)
            destination%field%r_ptr(:,:,:,:,:) = &
            & merge(destination%field%r_ptr (:,:,:,:,:), &
            &       source%field%r_ptr(:,index:index,:,:,:), &
            &       destination%field%r_ptr(:,:,:,:,:).gt.source%field%r_ptr(:,index:index,:,:,:))
          CASE(3)
            destination%field%r_ptr(:,:,:,:,:) = &
            & merge(destination%field%r_ptr (:,:,:,:,:), &
            &       source%field%r_ptr(:,:,index:index,:,:), &
            &       destination%field%r_ptr(:,:,:,:,:).gt.source%field%r_ptr(:,:,index:index,:,:))
          CASE(4)
            destination%field%r_ptr(:,:,:,:,:) = &
            & merge(destination%field%r_ptr (:,:,:,:,:), &
            &       source%field%r_ptr(:,:,:,index:index,:), &
            &       destination%field%r_ptr(:,:,:,:,:).gt.source%field%r_ptr(:,:,:,index:index,:))
          END SELECT
        ELSE
          ! tackle 1d, 2d and 3d
          SELECT CASE(destination%field%info%ndims)
          CASE(3)
           !hack for sea ice variables which uses vertical level for ice class
            IF (1 == destination%field%info%used_dimensions(2)) THEN
              call max_fields(destination%field%r_ptr(:,:,:,1,1), &
                &             source%field%r_ptr(:,:,:,1,1), &
                &             destination%field%info%subset, &
                &             levels=1, &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            ELSE IF (ZA_OCEAN_SEDIMENT == destination%field%info%vgrid) THEN
              call max_fields(destination%field%r_ptr(:,:,:,1,1), &
                &             source%field%r_ptr(:,:,:,1,1), &
                &             destination%field%info%subset, &
                &             levels=destination%field%info%used_dimensions(2), &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            ELSE
              call max_fields(destination%field%r_ptr(:,:,:,1,1), &
                &             source%field%r_ptr(:,:,:,1,1), &
                &             destination%field%info%subset, &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            ENDIF
          CASE(2)
            IF (GRID_ZONAL .EQ. destination%field%info%hgrid) THEN
              call max_fields(destination%field%r_ptr(:,:,1,1,1), &
                &             source%field%r_ptr(:,:,1,1,1), &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            ELSE
              call max_fields(destination%field%r_ptr(:,:,1,1,1), &
                &             source%field%r_ptr(:,:,1,1,1), &
                &             destination%field%info%subset, &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            ENDIF
          CASE DEFAULT
            destination%field%r_ptr(:,:,:,:,:) = &
                & merge(destination%field%r_ptr, &
                &       source%field%r_ptr, &
                &       destination%field%r_ptr.gt.source%field%r_ptr)
          END SELECT
        ENDIf
      END SUBROUTINE max_wp
      SUBROUTINE max_sp()
        IF (source%field%info%lcontained) THEN
          index = source%field%info%ncontained

          ! tackle 4d containers
          SELECT CASE(source%field%info%var_ref_pos)
          CASE(1)
            destination%field%r_ptr(:,:,:,:,:) = &
            & merge(destination%field%r_ptr (:,:,:,:,:), &
            &       REAL(source%field%s_ptr(index:index,:,:,:,:),wp), &
            &       destination%field%r_ptr(:,:,:,:,:).gt.REAL(source%field%r_ptr(index:index,:,:,:,:),wp))
          CASE(2)
            destination%field%r_ptr(:,:,:,:,:) = &
            & merge(destination%field%r_ptr (:,:,:,:,:), &
            &       REAL(source%field%s_ptr(:,index:index,:,:,:),wp), &
            &       destination%field%r_ptr(:,:,:,:,:).gt.REAL(source%field%r_ptr(:,index:index,:,:,:),wp))
          CASE(3)
            destination%field%r_ptr(:,:,:,:,:) = &
            & merge(destination%field%r_ptr (:,:,:,:,:), &
            &       REAL(source%field%s_ptr(:,:,index:index,:,:),wp), &
            &       destination%field%r_ptr(:,:,:,:,:).gt.REAL(source%field%r_ptr(:,:,index:index,:,:),wp))
          CASE(4)
            destination%field%r_ptr(:,:,:,:,:) = &
            & merge(destination%field%r_ptr (:,:,:,:,:), &
            &       REAL(source%field%s_ptr(:,:,:,index:index,:),wp), &
            &       destination%field%r_ptr(:,:,:,:,:).gt.REAL(source%field%r_ptr(:,:,:,index:index,:),wp))
          END SELECT
        ELSE
          ! tackle 1d, 2d and 3d
          SELECT CASE(destination%field%info%ndims)
          CASE(3)
           !hack for sea ice variables which uses vertical level for ice class
            IF (1 == destination%field%info%used_dimensions(2)) THEN
              call max_fields(destination%field%r_ptr(:,:,:,1,1), &
                &             source%field%s_ptr(:,:,:,1,1), &
                &             destination%field%info%subset, &
                &             levels=1, &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            ELSE IF (ZA_OCEAN_SEDIMENT == destination%field%info%vgrid) THEN
              call max_fields(destination%field%r_ptr(:,:,:,1,1), &
                &             source%field%s_ptr(:,:,:,1,1), &
                &             destination%field%info%subset, &
                &             levels=destination%field%info%used_dimensions(2), &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            ELSE
              call max_fields(destination%field%r_ptr(:,:,:,1,1), &
                &             source%field%s_ptr(:,:,:,1,1), &
                &             destination%field%info%subset, &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            ENDIF
          CASE(2)
            IF (GRID_ZONAL .EQ. destination%field%info%hgrid) THEN
              call max_fields(destination%field%r_ptr(:,:,1,1,1), &
                &             source%field%s_ptr(:,:,1,1,1), &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            ELSE
              call max_fields(destination%field%r_ptr(:,:,1,1,1), &
                &             source%field%s_ptr(:,:,1,1,1), &
                &             destination%field%info%subset, &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            ENDIF
          CASE DEFAULT
            destination%field%r_ptr(:,:,:,:,:) = &
                & merge(destination%field%r_ptr, &
                &       REAL(source%field%s_ptr,wp), &
                &       destination%field%r_ptr.gt.REAL(source%field%s_ptr,wp))
          END SELECT
        ENDIf
      END SUBROUTINE max_sp
  END SUBROUTINE statistics_max

  SUBROUTINE statistics_min(source, destination, counter)
    type(t_list_element) , INTENT(IN)    :: source
    type(t_list_element) , INTENT(INOUT) :: destination
    integer, intent(inout)               :: counter

    integer :: index

    index = -1
    IF (.NOT. ASSOCIATED(source%field%r_ptr)) THEN
      CALL min_sp()
    ELSE
      CALL min_wp()
    ENDIF
    counter                 = counter + 1

    CONTAINS
      SUBROUTINE min_wp()
        IF (source%field%info%lcontained) THEN
          index = source%field%info%ncontained

          ! tackle 4d containers
          SELECT CASE(source%field%info%var_ref_pos)
          CASE(1)
            destination%field%r_ptr(:,:,:,:,:) = &
            & merge(destination%field%r_ptr (:,:,:,:,:), &
            &       source%field%r_ptr(index:index,:,:,:,:), &
            &       destination%field%r_ptr(:,:,:,:,:).lt.source%field%r_ptr(index:index,:,:,:,:))
          CASE(2)
            destination%field%r_ptr(:,:,:,:,:) = &
            & merge(destination%field%r_ptr (:,:,:,:,:), &
            &       source%field%r_ptr(:,index:index,:,:,:), &
            &       destination%field%r_ptr(:,:,:,:,:).lt.source%field%r_ptr(:,index:index,:,:,:))
          CASE(3)
            destination%field%r_ptr(:,:,:,:,:) = &
            & merge(destination%field%r_ptr (:,:,:,:,:), &
            &       source%field%r_ptr(:,:,index:index,:,:), &
            &       destination%field%r_ptr(:,:,:,:,:).lt.source%field%r_ptr(:,:,index:index,:,:))
          CASE(4)
            destination%field%r_ptr(:,:,:,:,:) = &
            & merge(destination%field%r_ptr (:,:,:,:,:), &
            &       source%field%r_ptr(:,:,:,index:index,:), &
            &       destination%field%r_ptr(:,:,:,:,:).lt.source%field%r_ptr(:,:,:,index:index,:))
          END SELECT
        ELSE
          ! tackle 1d, 2d and 3d
          SELECT CASE(destination%field%info%ndims)
          CASE(3)
           !hack for sea ice variables which uses vertical level for ice class
            IF (1 == destination%field%info%used_dimensions(2)) THEN
              call min_fields(destination%field%r_ptr(:,:,:,1,1), &
                &             source%field%r_ptr(:,:,:,1,1), &
                &             destination%field%info%subset, &
                &             levels=1, &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            ELSE IF (ZA_OCEAN_SEDIMENT == destination%field%info%vgrid) THEN
              call min_fields(destination%field%r_ptr(:,:,:,1,1), &
                &             source%field%r_ptr(:,:,:,1,1), &
                &             destination%field%info%subset, &
                &             levels=destination%field%info%used_dimensions(2), &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            ELSE
              call min_fields(destination%field%r_ptr(:,:,:,1,1), &
                &             source%field%r_ptr(:,:,:,1,1), &
                &             destination%field%info%subset, &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            ENDIF
          CASE(2)
            IF (GRID_ZONAL .EQ. destination%field%info%hgrid) THEN
              call min_fields(destination%field%r_ptr(:,:,1,1,1), &
                &             source%field%r_ptr(:,:,1,1,1), &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            ELSE
              call min_fields(destination%field%r_ptr(:,:,1,1,1), &
                &             source%field%r_ptr(:,:,1,1,1), &
                &             destination%field%info%subset, &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            ENDIF
          CASE DEFAULT
            destination%field%r_ptr(:,:,:,:,:) = &
                & merge(destination%field%r_ptr, &
                &       source%field%r_ptr, &
                &       destination%field%r_ptr.lt.source%field%r_ptr)
          END SELECT
        ENDIF
      END SUBROUTINE min_wp
      SUBROUTINE min_sp()
        IF (source%field%info%lcontained) THEN
          index = source%field%info%ncontained

          ! tackle 4d containers
          SELECT CASE(source%field%info%var_ref_pos)
          CASE(1)
            destination%field%r_ptr(:,:,:,:,:) = &
            & merge(destination%field%r_ptr (:,:,:,:,:), &
            &       REAL(source%field%s_ptr(index:index,:,:,:,:),wp), &
            &       destination%field%r_ptr(:,:,:,:,:).lt.REAL(source%field%s_ptr(index:index,:,:,:,:),wp))
          CASE(2)
            destination%field%r_ptr(:,:,:,:,:) = &
            & merge(destination%field%r_ptr (:,:,:,:,:), &
            &       REAL(source%field%s_ptr(:,index:index,:,:,:),wp), &
            &       destination%field%r_ptr(:,:,:,:,:).lt.REAL(source%field%r_ptr(:,index:index,:,:,:),wp))
          CASE(3)
            destination%field%r_ptr(:,:,:,:,:) = &
            & merge(destination%field%r_ptr (:,:,:,:,:), &
            &       REAL(source%field%s_ptr(:,:,index:index,:,:),wp), &
            &       destination%field%r_ptr(:,:,:,:,:).lt.REAL(source%field%r_ptr(:,:,index:index,:,:),wp))
          CASE(4)
            destination%field%r_ptr(:,:,:,:,:) = &
            & merge(destination%field%r_ptr (:,:,:,:,:), &
            &       REAL(source%field%s_ptr(:,:,:,index:index,:),wp), &
            &       destination%field%r_ptr(:,:,:,:,:).lt.REAL(source%field%r_ptr(:,:,:,index:index,:),wp))
          END SELECT
        ELSE
          ! tackle 1d, 2d and 3d
          SELECT CASE(destination%field%info%ndims)
          CASE(3)
           !hack for sea ice variables which uses vertical level for ice class
            IF (1 == destination%field%info%used_dimensions(2)) THEN
              call min_fields(destination%field%r_ptr(:,:,:,1,1), &
                &             source%field%s_ptr(:,:,:,1,1), &
                &             destination%field%info%subset, &
                &             levels=1, &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            ELSE IF (ZA_OCEAN_SEDIMENT == destination%field%info%vgrid) THEN
              call min_fields(destination%field%r_ptr(:,:,:,1,1), &
                &             source%field%s_ptr(:,:,:,1,1), &
                &             destination%field%info%subset, &
                &             levels=destination%field%info%used_dimensions(2), &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            ELSE
              call min_fields(destination%field%r_ptr(:,:,:,1,1), &
                &             source%field%s_ptr(:,:,:,1,1), &
                &             destination%field%info%subset, &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            ENDIF
          CASE(2)
            IF (GRID_ZONAL .EQ. destination%field%info%hgrid) THEN
              call min_fields(destination%field%r_ptr(:,:,1,1,1), &
                &             source%field%s_ptr(:,:,1,1,1), &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            ELSE
              call min_fields(destination%field%r_ptr(:,:,1,1,1), &
                &             source%field%s_ptr(:,:,1,1,1), &
                &             destination%field%info%subset, &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            ENDIF
          CASE DEFAULT
            destination%field%r_ptr(:,:,:,:,:) = &
                & merge(destination%field%r_ptr, &
                &       REAL(source%field%r_ptr,wp), &
                &       destination%field%r_ptr.lt.REAL(source%field%r_ptr,wp))
          END SELECT
        ENDIF
      END SUBROUTINE min_sp
  END SUBROUTINE statistics_min
  !>
  !! Execute the accumulation forall internal variables and compute mean values
  !! if the corresponding event is active
  !!
  subroutine update_statistics()
    call update_mvstream(meanMap, meanEvents, meanEventsActivity, meanPrognostics, meanPrognosticPointers,  meanVarCounter,MEAN)
    call update_mvstream(maxMap, maxEvents, maxEventsActivity, maxPrognostics, maxPrognosticPointers,  maxVarCounter,MAX)
    call update_mvstream(minMap, minEvents, minEventsActivity, minPrognostics, minPrognosticPointers,  minVarCounter,MIN)
    call update_mvstream(squareMap, &
      &                  squareEvents, &
      &                  squareEventsActivity, &
      &                  squarePrognostics, &
      &                  squarePrognosticPointers, &
      &                  squareVarCounter, &
      &                  SQUARE)
  end subroutine update_statistics

  SUBROUTINE update_mvstream(statisticMap, statisticEvents, statisticEventsActivity, &
          & statisticPrognostics, statisticPrognosticPointers, &
          & statisticVarCounter, operation)
    TYPE(map) :: statisticMap
    TYPE(map) :: statisticEvents
    TYPE(map) :: statisticEventsActivity
    TYPE(vector) :: statisticPrognostics
    TYPE(map) :: statisticPrognosticPointers
    TYPE(map) :: statisticVarCounter
    CHARACTER(LEN=*), INTENT(IN) :: operation

    !----------------------------------------------------------------------------------
    INTEGER :: element_counter
    class(*),pointer :: varListForMeanEvent, statisticEventKey
    class(*),pointer :: sourceVariable, destinationVariable, sourceVariable4Prognostics
    class(*),pointer :: counter, statisticEvent, myItem
    class(*),pointer :: eventActive
    type(t_list_element), pointer :: source, destination
    type(vector_iterator) :: statisticMapIterator, statisticEventIterator
    class(*), pointer :: mvstream_pair

    TYPE(datetime), POINTER :: mtime_date
    logical :: isactive
    integer :: varcounter
    integer :: timelevel
#ifdef DEBUG_MVSTREAM
    character(len=max_datetime_str_len) :: datetimestring
#endif

    CHARACTER(LEN=*), PARAMETER :: routine =  modname//"::update_statistics"

#ifdef DEBUG_MVSTREAM
    if (my_process_is_stdio()) call print_routine(routine,'start',stderr=.true.)
#endif

    statisticMapIterator   = statisticMap%iter()
    statisticEventIterator = statisticEvents%iter()

    ! Check events first {{{
    ! this is necessary because of mtime internals
    isactive   = .false.
    mtime_date => newDatetime(time_config%tc_current_date)

#ifdef DEBUG_MVSTREAM
    CALL datetimeToString(mtime_date, datetimestring)
    if (my_process_is_stdio()) call print_summary('Current mtime timestamp:'//trim(datetimestring),stderr=.true.)
#endif

    ! Save results for (not so much) later
    do while (statisticEventIterator%next(myItem))
      select type (myItem)
      type is (map_item)
        statisticEventKey => myItem%key
        statisticEvent    => myItem%value
        select type (statisticEvent)
        type is (t_event_wrapper)
          isactive = LOGICAL(isCurrentEventActive(statisticEvent%this,mtime_date))
        end select
        call statisticEventsActivity%add(statisticEventKey,isactive)
      end select
    end do

#ifdef DEBUG_MVSTREAM
    if (my_process_is_stdio()) call print_error(statisticEventsActivity%to_string(),stderr=.true.)
#endif
    ! }}}

    do while (statisticMapIterator%next(myItem))

      select type (myItem)
      type is (map_item)

        statisticEventKey        => myItem%key
        varListForMeanEvent => myItem%value

#ifdef DEBUG_MVSTREAM
        if (my_process_is_stdio()) call print_summary(object_pointer_string(statisticEventKey)//"PERFORM ACCU",stderr=.true.) !TODO
#endif

        select type(varListForMeanEvent)
        class is (vector_ref)
          do element_counter=1,varListForMeanEvent%length(),2

#ifdef DEBUG_MVSTREAM
          if (my_process_is_stdio()) call print_routine("update_statistics",object_string(element_counter),stderr=.true.)
#endif
            sourceVariable      => varListForMeanEvent%at(element_counter)
            destinationVariable => varListForMeanEvent%at(element_counter+1)
            if (associated(sourceVariable) .and. associated(destinationVariable)) then
              select type (sourceVariable)
              type is (t_list_element)
                select type (destinationVariable)
                type is (t_list_element)
                  ! check for prognostics, pointer must be shifted according to given timelevelIndex {{{
                  if (is_statsPrognosticVariable(destinationVariable,statisticPrognostics)) then
                    ! find the correct pointer:
                    ! check f reduced calling freq. variables are used for
                    ! output? if true, theses variables should be used for
                    ! accumulation, too

#ifdef DEBUG_MVSTREAM
                    if (my_process_is_stdio()) &
                      & call print_error("destination var IS     prognostic:"//TRIM(destinationVariable%field%info%name),&
                      & stderr=.true.)
#endif

                    timelevel =  metainfo_get_timelevel(destinationVariable%field%info, &
                        & destinationVariable%field%info%dom)

                    source    => get_prognostics_source_pointer(destinationVariable, timelevel, statisticPrognosticPointers)

                  else
                  ! }}}
                    source    => sourceVariable
#ifdef DEBUG_MVSTREAM
                    if (my_process_is_stdio()) &
                      & call print_error("destination var IS NOT prognostic:"//TRIM(destinationVariable%field%info%name),&
                      & stderr=.true.)
#endif
                  end if
                  destination => destinationVariable
                  counter     => statisticVarCounter%get(destination%field%info%name)
                  select type (counter)
                  type is (integer)

#ifdef DEBUG_MVSTREAM
                    IF ( my_process_is_stdio() ) call print_summary('sourceName : '//trim(source%field%info%name),&
                        &stderr=.true.)
                    IF ( my_process_is_stdio() ) call print_summary('destName   : '//trim(destination%field%info%name),&
                        & stderr=.true.)
                    IF ( my_process_is_stdio() ) call print_summary('destNameOut: '//trim(destination%field%info%cf%short_name),&
                        &stderr=.true.)
                    IF ( my_process_is_stdio() ) write (0,*)'old counter: ',counter
#endif

                    ! update of internal fields according to the operator {{{
                    varcounter = counter !TODO work around for integer pointer, ugly
                    if (MEAN.EQ.operation) then
                      CALL statistics_add(source, destination, varcounter)

                    elseif (MAX.EQ.operation) then
                      if (0 .EQ. varcounter) then
                        CALL statistics_assign(source, destination, varcounter)
                      else
                        CALL statistics_max(source, destination, varcounter)
                      endif

                    elseif (MIN.EQ.operation) then
                      if (0 .EQ. varcounter) then
                        CALL statistics_assign(source, destination, varcounter)
                      else
                        CALL statistics_min(source, destination, varcounter)
                      endif

                    elseif (SQUARE.EQ.operation) then
                      CALL statistics_add_sqr(source, destination, varcounter)

                    else
                      CALL finish('update_mvstream','Found unknown operation:'//trim(operation))

                    endif
                    counter = varcounter
                    ! }}}
#ifdef DEBUG_MVSTREAM
                    IF ( my_process_is_stdio() )  write (0,*)'new counter: ',counter
#endif

                  end select

                  ! MEAN VALUE COMPUTAION {{{
                  ! check if the field will be written to disk this timestep {{{
                  eventActive => statisticEventsActivity%get(statisticEventKey)
                  select type (eventActive)
                  type is (logical)
                    isactive = eventActive
                    if ( isactive ) then

#ifdef DEBUG_MVSTREAM
                    if (my_process_is_stdio()) &
                        & CALL print_summary(" --> PERFORM MEAN VALUE COMP for"//trim(destination%field%info%name),stderr=.true.)
#endif

                      counter => statisticVarCounter%get(destination%field%info%name)
                      select type(counter)
                      type is (integer)
#ifdef DEBUG_MVSTREAM
                      IF ( my_process_is_stdio() ) write (0,*)' ------> MEAN VALUE counter:',counter
#endif
                      if (MEAN.EQ.operation .OR. SQUARE.EQ.operation) then
                        destination%field%r_ptr = destination%field%r_ptr / REAL(counter,wp)
                      endif
                      counter = 0
                      end select
                    end if
                  end select
                  ! }}}
                end select
                ! ! }}}
              class default
                call finish('update_statistics','Found unknown source variable type')
              end select
            else
              call finish(routine,'source or destination variable cannot be found')
            end if
          end do
        end select
      end select
    end do

#ifdef DEBUG_MVSTREAM
    if (my_process_is_stdio()) call print_routine(routine,'finish')
#endif

  END SUBROUTINE update_mvstream

  !>
  !! reset internal fields to zero, if the corresponding event is active
  !!
  ! min/max do not need extra reset due to JEFs implementation idea
  SUBROUTINE reset_statistics
    call reset_mvstream(meanMap, meanEventsActivity)
    call reset_mvstream(squareMap, squareEventsActivity)
  END SUBROUTINE reset_statistics


  SUBROUTINE reset_mvstream(statisticMap, statisticEventsActivity)
    type(map) :: statisticMap, statisticEventsActivity

    INTEGER :: element_counter
    type(t_list_element), pointer :: destination
    class(*),pointer :: varListForStatisticEvent, statisticEventKey
    class(*),pointer :: destinationVariable
    class(*),pointer :: statisticEvent, myItem
    class(*), pointer :: eventActive
    type(vector_iterator) :: statisticMapIterator

    CHARACTER(LEN=*), PARAMETER :: routine =  modname//"::reset_statistics"

#ifdef DEBUG_MVSTREAM
    if (my_process_is_stdio()) call print_routine(routine,'start')
#endif

    statisticMapIterator = statisticMap%iter()

    do while (statisticMapIterator%next(myItem))

#ifdef DEBUG_MVSTREAM
      if (my_process_is_stdio()) call print_error("started event loop",stderr=.true.)
#endif
      select type (myItem)
      type is (map_item)

        statisticEventKey   => myItem%key
        varListForStatisticEvent => myItem%value

        select type(varListForStatisticEvent)
        class is (vector_ref)

          do element_counter=1,varListForStatisticEvent%length(),2 !start at 2 because the event is at index 1

            destinationVariable => varListForStatisticEvent%at(element_counter+1)

#ifdef DEBUG_MVSTREAM
            if (my_process_is_stdio()) call print_error("got destinationVariable",stderr=.true.)
#endif

            if (associated(destinationVariable)) then

#ifdef DEBUG_MVSTREAM
            if (my_process_is_stdio()) call print_error("    destinationVariable is associated",stderr=.true.)
#endif

              select type (destinationVariable)
              type is (t_list_element)
                  destination => destinationVariable

#ifdef DEBUG_MVSTREAM
                  if (my_process_is_stdio()) call print_error("    destinationVariable is t_list_element",stderr=.true.)
#endif

                  eventActive => statisticEventsActivity%get(statisticEventKey)

#ifdef DEBUG_MVSTREAM
                  if (.not.associated(eventActive)) then
                    if (my_process_is_stdio()) call print_error("       eventActive not associated",stderr=.true.)
                  end if
#endif

                  select type (eventActive)
                  type is (logical)

#ifdef DEBUG_MVSTREAM
                    if (my_process_is_stdio()) call print_error("       eventActive is logical",stderr=.true.)
#endif

                    if ( LOGICAL(eventActive) ) then

#ifdef DEBUG_MVSTREAM
                      if (my_process_is_stdio()) call print_error("       eventActive is true",stderr=.true.)
                      if (my_process_is_stdio()) call print_error(object_string(statisticEventKey)//' : --> PERFORM RESET',&
                          & stderr=.true.)
if (my_process_is_stdio()) call print_error(object_string(statisticEventKey)//' : --> '//trim(destination%field%info%name),&
    & stderr=.true.)
#endif

                      destination%field%r_ptr = 0.0_wp ! take the neutral element of addition
#ifdef DEBUG_MVSTREAM
                    else
                      if (my_process_is_stdio()) call print_error("       eventActive is false",stderr=.true.)
#endif
                    end if
                  class default
#ifdef DEBUG_MVSTREAM
                      if (my_process_is_stdio()) call print_error("       eventActive has wrong type",stderr=.true.)
#endif
                  end select
              class default
#ifdef DEBUG_MVSTREAM
                  if (my_process_is_stdio()) call print_error("     destinationVariable is not t_list_element",stderr=.true.)
#endif
              end select
            else
#ifdef DEBUG_MVSTREAM
              if (my_process_is_stdio()) call print_error(routine//TRIM(": cannot find destination variable!"),stderr=.true.)
#endif
            end if
          end do
        end select
      end select
    end do

#ifdef DEBUG_MVSTREAM
    if (my_process_is_stdio()) call print_routine(routine,'finish',stderr=.true.)
#endif
  END SUBROUTINE

  SUBROUTINE reset_mean
    INTEGER :: element_counter
    type(t_list_element), pointer :: destination
    class(*),pointer :: varListForMeanEvent, meanEventKey
    class(*),pointer :: destinationVariable
    class(*),pointer :: meanEvent, myItem
    class(*), pointer :: eventActive
    type(vector_iterator) :: meanMapIterator

    CHARACTER(LEN=*), PARAMETER :: routine =  modname//"::reset_statistics"

#ifdef DEBUG_MVSTREAM
    if (my_process_is_stdio()) call print_routine(routine,'start')
#endif

    meanMapIterator   = meanMap%iter()

    do while (meanMapIterator%next(myItem))

#ifdef DEBUG_MVSTREAM
      if (my_process_is_stdio()) call print_error("started event loop",stderr=.true.)
#endif
      select type (myItem)
      type is (map_item)

        meanEventKey        => myItem%key
        varListForMeanEvent => myItem%value

        select type(varListForMeanEvent)
        class is (vector_ref)

          do element_counter=1,varListForMeanEvent%length(),2 !start at 2 because the event is at index 1

            destinationVariable => varListForMeanEvent%at(element_counter+1)

#ifdef DEBUG_MVSTREAM
            if (my_process_is_stdio()) call print_error("got destinationVariable",stderr=.true.)
#endif

            if (associated(destinationVariable)) then

#ifdef DEBUG_MVSTREAM
            if (my_process_is_stdio()) call print_error("    destinationVariable is associated",stderr=.true.)
#endif

              select type (destinationVariable)
              type is (t_list_element)
                  destination => destinationVariable

#ifdef DEBUG_MVSTREAM
                  if (my_process_is_stdio()) call print_error("    destinationVariable is t_list_element",stderr=.true.)
#endif

                  eventActive => meanEventsActivity%get(meanEventKey)

#ifdef DEBUG_MVSTREAM
                  if (.not.associated(eventActive)) then
                    if (my_process_is_stdio()) call print_error("       eventActive not associated",stderr=.true.)
                  end if
#endif

                  select type (eventActive)
                  type is (logical)

#ifdef DEBUG_MVSTREAM
                    if (my_process_is_stdio()) call print_error("       eventActive is logical",stderr=.true.)
#endif

                    if ( LOGICAL(eventActive) ) then

#ifdef DEBUG_MVSTREAM
                      if (my_process_is_stdio()) call print_error("       eventActive is true",stderr=.true.)
                      if (my_process_is_stdio()) call print_error(object_string(meanEventKey)//' : --> PERFORM RESET',&
                          & stderr=.true.)
if (my_process_is_stdio()) call print_error(object_string(meanEventKey)//' : --> '//trim(destination%field%info%name),&
    & stderr=.true.)
#endif

                      destination%field%r_ptr = 0.0_wp ! take the neutral element of addition
#ifdef DEBUG_MVSTREAM
                    else
                      if (my_process_is_stdio()) call print_error("       eventActive is false",stderr=.true.)
#endif
                    end if
                  class default
#ifdef DEBUG_MVSTREAM
                      if (my_process_is_stdio()) call print_error("       eventActive has wrong type",stderr=.true.)
#endif
                  end select
              class default
#ifdef DEBUG_MVSTREAM
                  if (my_process_is_stdio()) call print_error("     destinationVariable is not t_list_element",stderr=.true.)
#endif
              end select
            else
#ifdef DEBUG_MVSTREAM
              if (my_process_is_stdio()) call print_error(routine//TRIM(": cannot find destination variable!"),stderr=.true.)
#endif
            end if
          end do
        end select
      end select
    end do

#ifdef DEBUG_MVSTREAM
    if (my_process_is_stdio()) call print_routine(routine,'finish',stderr=.true.)
#endif

  END SUBROUTINE

END MODULE mo_derived_variable_handling

! vim:tw=0
