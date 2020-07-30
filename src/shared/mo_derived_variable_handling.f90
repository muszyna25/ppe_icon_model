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

  USE mo_kind, ONLY: wp, sp
  USE mo_model_domain, ONLY: t_patch
  USE mo_io_config, ONLY: lnetcdf_flt64_output
  USE mo_dynamics_config, ONLY: nnow, nnew, nold
  USE mo_statistics, ONLY: add_fields, max_fields, min_fields, assign_fields, add_sqr_fields
  USE mo_impl_constants, ONLY: vname_len, SUCCESS, max_char_length, TLEV_NNOW, TLEV_NNEW,    &
    &                          TLEV_NNOW_RCF, TLEV_NNEW_RCF, REAL_T
  USE mo_cdi_constants, ONLY: GRID_UNSTRUCTURED_CELL, GRID_UNSTRUCTURED_EDGE, GRID_UNSTRUCTURED_VERT, &
                              & GRID_ZONAL
  USE mo_name_list_output_types, ONLY: t_output_name_list
  USE mo_zaxis_type, ONLY:  ZA_OCEAN_SEDIMENT
  USE mo_mpi, ONLY: my_process_is_stdio
  USE mo_var_list_element, ONLY: level_type_ml, level_type_pl, level_type_hl, level_type_il
  USE mo_name_list_output_metadata, ONLY: metainfo_get_timelevel
  USE mo_var_list, ONLY: get_var_name, add_var, get_varname_with_timelevel
  USE mo_var_list_global, ONLY: find_element => find_var_global, new_var_list, &
    & total_number_of_variables, delete_var_list, print_all_var_lists
  USE mo_linked_list, ONLY: t_var_list, t_list_element
  USE mo_exception, ONLY: finish, message, message_text
  USE mtime, ONLY: MAX_DATETIME_STR_LEN, newEvent, event, isCurrentEventActive,&
    & newDatetime, datetime, eventToString, divideDatetimeDifferenceInSeconds, &
    & divisionquotientTimespan
  USE mo_util_string, ONLY: int2string
  USE mtime_datetime, ONLY: datetimeToString
  USE mo_output_event_types,  ONLY: t_sim_step_info
  USE mo_time_config,         ONLY: time_config
  USE mo_cdi,                 ONLY: DATATYPE_FLT32, DATATYPE_FLT64, GRID_LONLAT, TSTEP_CONSTANT
#ifdef _OPENACC
  USE mo_mpi,                 ONLY: i_am_accel_node
#endif

#include "add_var_acc_macro.inc"

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
  ! in "SELECT TYPE" - this is kind of a hack
  TYPE :: t_event_wrapper
    TYPE(event), pointer :: this
  end TYPE

CONTAINS

  !>
  !! Create all needed lists and maps
  !!
  SUBROUTINE init_statistics_streams()
    CHARACTER(LEN=max_char_length) :: listname

    ! main map for automatical mean value computation:
    ! key: eventString
    ! value: self-vector of t_list_elements objects that belong to that event
    CALL meanMap%init(verbose=.false.)
    ! map for holding the events itself
    ! since the events cannot be saved a keys for a map because they are only
    ! C-pointers the same string representation of events is used as for
    ! "meanMap"
    CALL meanEvents%init(verbose=.false.)
    CALL meanEventsActivity%init(verbose=.false.)
    CALL meanVarCounter%init(verbose=.false.)
    CALL meanPrognostics%init(verbose=.false.)
    CALL meanPrognosticPointers%init(verbose=.false.)

    CALL maxMap%init(verbose=.false.)
    CALL maxEvents%init(verbose=.false.)
    CALL maxEventsActivity%init(verbose=.false.)
    CALL maxVarCounter%init(verbose=.false.)
    CALL maxPrognostics%init(verbose=.false.)
    CALL maxPrognosticPointers%init(verbose=.false.)

    CALL minMap%init(verbose=.false.)
    CALL minEvents%init(verbose=.false.)
    CALL minEventsActivity%init(verbose=.false.)
    CALL minVarCounter%init(verbose=.false.)
    CALL minPrognostics%init(verbose=.false.)
    CALL minPrognosticPointers%init(verbose=.false.)

    CALL squareMap%init(verbose=.false.)
    CALL squareEvents%init(verbose=.false.)
    CALL squareEventsActivity%init(verbose=.false.)
    CALL squareVarCounter%init(verbose=.false.)
    CALL squarePrognostics%init(verbose=.false.)
    CALL squarePrognosticPointers%init(verbose=.false.)
  END SUBROUTINE init_statistics_streams

  !>
  !! Print contents of a vector, just giving the name for t_list_elements
  !!
  !! Optional label is printed first, on a line by its own
  !!
  SUBROUTINE var_print(this, label)
    CLASS(vector_ref) , INTENT(in) :: this
    CHARACTER(*), INTENT(in), OPTIONAL :: label

    TYPE(vector_iterator) :: my_iter
    CLASS(*), POINTER :: my_buffer
    INTEGER :: i

    IF (PRESENT(label)) PRINT *, label

    my_iter = this%iter()
    DO WHILE(my_iter%next(my_buffer))
      SELECT TYPE(my_buffer)
      TYPE IS (t_list_element)
        PRINT *,'t_list_element:varname:',         trim(my_buffer%field%info%name)
      TYPE IS (t_mvstream_pair)
        PRINT *,'t_mvstream_pair:source     :',trim(my_buffer%source%field%info%name)
        PRINT *,'t_mvstream_pair:destination:',trim(my_buffer%destination%field%info%name)
      CLASS DEFAULT
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
  LOGICAL FUNCTION is_statsPrognosticVariable(variable,statList)
    TYPE(t_list_element), INTENT(in) :: variable
    TYPE(vector), INTENT(in) :: statList

    is_statsPrognosticVariable = statList%includes(variable%field%info%name)
  END FUNCTION

  !>
  !! return pointer to prognostic accumulation copy for given timelevel
  !!
  FUNCTION get_prognostics_source_pointer(destinationVariable, timelevelIndex, pointerMap) result(sourcePointer)
    TYPE(t_list_element), pointer :: sourcePointer
    TYPE(t_list_element)          :: destinationVariable
    TYPE(map)                     :: pointerMap

    CHARACTER(LEN=*), PARAMETER :: routine =  modname//"::get_prognostics_source_pointer"
    INTEGER , INTENT(IN) :: timelevelIndex

    CLASS(*), pointer :: dummy

    dummy => pointerMap%get( &
      & get_varname_with_timelevel( &
      &   destinationVariable%field%info%name, &
      &   timelevelIndex) &
      & )

    SELECT TYPE (dummy)
    TYPE IS (t_list_element)
      sourcePointer => dummy
    CLASS DEFAULT
      CALL finish(routine,"Cound not find pointer to prognostics variable for given timelevel!")
    END SELECT
  END FUNCTION get_prognostics_source_pointer

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
      CALL finish(routine,"meanValues are only supported on global domain 1 and without stream partitioning!")
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

  INTEGER FUNCTION output_varlist_length(in_varlist)
    CHARACTER(LEN=vname_len), INTENT(in) :: in_varlist(:)

    output_varlist_length = 0
    DO
      IF (in_varlist(output_varlist_length+1) == ' ') EXIT
      output_varlist_length = output_varlist_length + 1
    END DO
  END FUNCTION output_varlist_length

  SUBROUTINE fill_list_of_output_varnames(varlist,in_varlist,ntotal_vars)
    CHARACTER(LEN=vname_len), INTENT(INOUT) :: varlist(:)
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
    CHARACTER(LEN=vname_len) :: dest_element_name
    CHARACTER(LEN=vname_len), POINTER :: in_varlist(:)
    INTEGER,INTENT(IN) :: i

      IF (my_process_is_stdio()) CALL print_summary('src(name)     :|'//trim(src_element%field%info%name)//'|',stderr=.true.)
      IF (my_process_is_stdio()) CALL print_summary('src(grid)     :|'//int2string(src_element%field%info%hgrid)//'|', &
          & stderr=.true.)
      IF (my_process_is_stdio()) CALL print_summary('src(dims:1)   :|'//int2string(src_element%field%info%used_dimensions(1))//'|',&
          & stderr=.true.)
      IF (my_process_is_stdio()) CALL print_summary('src(dims:2)   :|'//int2string(src_element%field%info%used_dimensions(2))//'|',&
          & stderr=.true.)
      IF (my_process_is_stdio()) CALL print_summary('src(dims:3)   :|'//int2string(src_element%field%info%used_dimensions(3))//'|',&
          & stderr=.true.)
      IF (my_process_is_stdio()) CALL print_summary('src(dims:4)   :|'//int2string(src_element%field%info%used_dimensions(4))//'|',&
          & stderr=.true.)
      IF (my_process_is_stdio()) CALL print_summary('src(dims:5)   :|'//int2string(src_element%field%info%used_dimensions(5))//'|',&
          & stderr=.true.)
      IF (my_process_is_stdio()) CALL print_summary('varlist(name) :|'//trim(in_varlist(i))//'|',stderr=.true.)
      IF (my_process_is_stdio()) CALL print_summary('new name      :|'//trim(dest_element_name)//'|',stderr=.true.)
      IF (my_process_is_stdio()) CALL print_summary('new grid      :|'//int2string(dest_element%field%info%hgrid)//'|',&
          & stderr=.true.)
  END SUBROUTINE

  SUBROUTINE find_src_element(src_element, varlist_element, dest_element_name, &
          &  dom, src_list, prognosticsList, prognosticsPointerList)
    TYPE(t_list_element), POINTER :: src_element
    CHARACTER(LEN=vname_len), INTENT(IN) :: varlist_element
    CHARACTER(LEN=vname_len) :: dest_element_name
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
    IF (.NOT. ASSOCIATED(src_element) ) &
        & src_element => find_element ( TRIM(varlist_element),&
        &                               opt_patch_id=dom, &
        &                               opt_hgrid=GRID_UNSTRUCTURED_EDGE, &
        &                               opt_caseInsensitive=.true., &
        &                               opt_returnList=src_list)
    IF (.NOT. ASSOCIATED(src_element) ) &
        & src_element => find_element ( TRIM(varlist_element),&
        &                               opt_patch_id=dom, &
        &                               opt_hgrid=GRID_UNSTRUCTURED_VERT, &
        &                               opt_caseInsensitive=.true., &
        &                               opt_returnList=src_list)
    IF (.NOT. ASSOCIATED(src_element) ) &
        & src_element => find_element ( TRIM(varlist_element),&
        &                               opt_patch_id=dom, &
        &                               opt_hgrid=GRID_LONLAT, &
        &                               opt_caseInsensitive=.true., &
        &                               opt_returnList=src_list)
    IF (.NOT. ASSOCIATED(src_element) ) &
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
    IF (.NOT. ASSOCIATED (src_element)) THEN
      foundPrognostic = .false.
      timelevels = (/nold(1),nnow(1),nnew(1)/)
      DO timelevel=1,3

#ifdef DEBUG_MVSTREAM
          IF (my_process_is_stdio()) CALL print_error(get_varname_with_timelevel(varlist_element,timelevels(timelevel)))
#endif
          src_element => find_element(get_varname_with_timelevel(varlist_element,timelevels(timelevel)),&
              &                       opt_patch_id=dom, opt_returnList=src_list)
          IF ( ASSOCIATED(src_element) ) THEN
#ifdef DEBUG_MVSTREAM
            IF (my_process_is_stdio()) write(0,*)'found prognostic:',&
                & TRIM(get_varname_with_timelevel(varlist_element,timelevels(timelevel)))
#endif
            IF ( .NOT. foundPrognostic ) THEN
              ! save the name of the original output variable if a prognosting version was found
              CALL prognosticsList%add(dest_element_name)
              foundPrognostic = .true.
#ifdef DEBUG_MVSTREAM
              IF (my_process_is_stdio()) CALL print_error('prognosticsList%add():'//varlist_element)
#endif
            END IF
            ! save the the pointers for all time levels of a prognostic variable
            ! these must be used for correct accumulation during the time loop
            CALL prognosticsPointerList%add(get_varname_with_timelevel(dest_element_name,timelevels(timelevel)),src_element)
          END IF
      END DO
    END IF
    ! in case nothing appropriate could be found, throw an error
    IF (.NOT. ASSOCIATED (src_element)) THEN
      CALL finish(routine,'Could not find source variable:'//TRIM(varlist_element))
    END IF
  END SUBROUTINE

  SUBROUTINE setup_statistics_events_and_lists(statsMap, eventTag, statisticsVariablesForEvent, eventMap, eventsActivityMap,p_onl)
    TYPE(map) :: statsMap, eventMap, eventsActivityMap
    CHARACTER(LEN=100), INTENT(in) :: eventTag
    TYPE(vector_ref) :: statisticsVariablesForEvent
    TYPE (t_output_name_list), INTENT(in),target  :: p_onl

    CLASS(*), pointer :: myBuffer
    TYPE(t_event_wrapper) :: event_wrapper

    IF ( statsMap%has_key(eventTag) ) THEN
      myBuffer => statsMap%get(eventTag)
      SELECT TYPE (myBuffer)
      TYPE IS (vector_ref)
        ! use exiting list
        statisticsVariablesForEvent = myBuffer
      END SELECT
    ELSE
      ! create new variable list
      CALL statisticsVariablesForEvent%init()

      ! wrap c-pointer for events into a fortran type - look weird, but does the trick
      event_wrapper = t_event_wrapper(this=newEvent(eventTag, &
        &                 p_onl%output_start(1), &
        &                 p_onl%output_start(1), &
        &                 p_onl%output_end(1), &
        &                 p_onl%output_interval(1) &
        &                ))

      ! collect the new event
      CALL eventMap%add(eventTag,event_wrapper)
    END IF
    ! initialize all events as in-active
    CALL eventsActivityMap%add(eventTag,.false.)
  END SUBROUTINE

  SUBROUTINE print_prognosticLists(varlist,progMap,progPointers)
    TYPE(vector_ref) :: varlist
    TYPE(map) :: progPointers, progMap

    CHARACTER(LEN=*), PARAMETER :: routine =  modname//"::print_prognosticLists"

    CALL var_print(varlist)
    CALL print_routine(routine, 'progMap%to_string()')
    CALL print_routine(routine, progMap%to_string())
    CALL print_routine(routine, progPointers%to_string())
    CALL var_print(progPointers%values())
  END SUBROUTINE

  !>
  !! copy varlist element (source_element) to target varlist (list) with given name
  !! lrestart = .true.
  !!
  FUNCTION copy_var_to_list(list,name,source_element,patch_2d) RESULT(dest_element)
    TYPE(t_var_list) :: list
    CHARACTER(LEN=vname_len) :: name
    TYPE(t_list_element),POINTER :: source_element
    TYPE(t_patch),TARGET, INTENT(in)    :: patch_2d

    TYPE(t_list_element), POINTER :: dest_element

    INTEGER :: dataType = -1
    CHARACTER(LEN=*), PARAMETER :: routine =  modname//"::copy_var_to_list"

    dataType = MERGE(DATATYPE_FLT64, DATATYPE_FLT32, lnetcdf_flt64_output)
#ifdef DEBUG_MVSTREAM
    IF (my_process_is_stdio()) THEN
      CALL print_verbose("COPY variable:"//TRIM(name))
      CALL print_verbose("     INTO    :"//TRIM(list%p%name))
    END IF
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
      & var_class=source_element%field%info%var_class, &
      & lopenacc = .TRUE. )
!!!!      & lopenacc = source_element%field%info%lopenacc )   ! Not clear why this is not .TRUE.
    __acc_attach( dest_element%field%r_ptr )

    ! add the subset for later accumulation on all types of horizontal grids
    SELECT case (source_element%field%info%hgrid)
    case (GRID_UNSTRUCTURED_CELL)
      dest_element%field%info%subset = patch_2d%cells%owned
    case (GRID_UNSTRUCTURED_EDGE)
      dest_element%field%info%subset = patch_2d%edges%owned
    case (GRID_UNSTRUCTURED_VERT)
      dest_element%field%info%subset = patch_2d%verts%owned
    END SELECT
  END FUNCTION copy_var_to_list

  !>
  !! return internal name for accumulation variables
  !!
  FUNCTION get_statistics_varname(varname,output_setup)
    CHARACTER(LEN=vname_len)  :: varname, get_statistics_varname
    type(t_output_name_list) :: output_setup

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
    CHARACTER(LEN=vname_len)  :: mean_varname, get_real_varname

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

    TYPE(map) :: statisticMap
    TYPE(map) :: statisticEvents
    TYPE(map) :: statisticEventsActivity
    TYPE(map) :: statisticVarCounter
    TYPE(vector) :: statisticPrognostics
    TYPE(map) :: statisticPrognosticPointers
    CHARACTER(len=*), INTENT(in) :: operation

    !-----------------------------------------------------------------------------
    CHARACTER(LEN=vname_len), POINTER :: in_varlist(:)
    INTEGER :: ntotal_vars, output_variables,i,ierrstat, dataType
    INTEGER :: timelevel, timelevels(3)
    TYPE(vector_ref) :: statisticVariables, prognosticVariables
    CHARACTER(LEN=100) :: eventKey
    TYPE(t_event_wrapper) :: event_wrapper
    CLASS(*), pointer :: myBuffer
    TYPE(t_list_element), POINTER :: src_element, dest_element
    CHARACTER(LEN=vname_len), ALLOCATABLE :: varlist(:)
    CHARACTER(LEN=vname_len) :: dest_element_name
    LOGICAL :: foundPrognostic
    TYPE(t_var_list), POINTER :: src_list
    CHARACTER(LEN=*), PARAMETER :: routine =  modname//"::process_mvstream"

#ifdef DEBUG_MVSTREAM
    IF (my_process_is_stdio()) CALL print_routine(routine,'start')
#endif

    IF (trim(operation) .EQ. TRIM(p_onl%operation)) THEN

#ifdef DEBUG_MVSTREAM
      IF (my_process_is_stdio()) CALL print_routine(routine,'found "'//trim(operation)//'" operation')
#endif
      CALL statStreamCrossCheck(p_onl)

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
      IF (my_process_is_stdio()) CALL print_summary('eventKey:'//trim(eventKey))
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
        IF (my_process_is_stdio()) CALL print_summary('destination variable NAME:'//TRIM(dest_element_name),stderr=.true.)
#endif
        IF (.NOT. ASSOCIATED(dest_element) ) THEN !not found -->> create a new on
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
          IF ( my_process_is_stdio()) CALL print_summary('dst(name)     :|'//trim(dest_element%field%info%name)//'|',&
              & stderr=.true.)
          IF ( my_process_is_stdio()) CALL print_summary('dst(shortname):|'//trim(dest_element%field%info%cf%short_name)//'|',&
              & stderr=.true.)
#endif
        END IF
        in_varlist(i) = trim(dest_element%field%info%name)
      END DO
      CALL statisticMap%add(eventKey,statisticVariables)
#ifdef DEBUG_MVSTREAM
      IF (my_process_is_stdio()) CALL print_error(routine//": statisticPrognostics%to_string()")
      IF (my_process_is_stdio()) CALL print_error(statisticPrognostics%to_string())
#endif
      DEALLOCATE(varlist)
    ELSE
#ifdef DEBUG_MVSTREAM
      IF (my_process_is_stdio()) CALL print_routine(routine,'NO "statistic" operation found')
#endif
    END IF

#ifdef DEBUG_MVSTREAM
    IF (my_process_is_stdio()) THEN
      CALL print_prognosticLists(statisticVariables, statisticMap, statisticPrognosticPointers)
    END IF
    IF (my_process_is_stdio()) CALL print_routine(routine,'end',stderr=.true.)
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
    TYPE(t_list_element) , INTENT(IN)    :: source
    TYPE(t_list_element) , INTENT(INOUT) :: destination
    INTEGER, INTENT(inout)               :: counter

    INTEGER :: index

    index = -1

    IF (.NOT. ASSOCIATED(source%field%r_ptr)) THEN
      CALL add_sp()
    ELSE
      CALL add_wp()
    END IF
    counter                 = counter + 1

    CONTAINS
      SUBROUTINE add_wp()
        REAL(wp), POINTER :: dst_ptr(:,:,:,:,:), src_ptr(:,:,:,:,:)
        dst_ptr => destination%field%r_ptr
        src_ptr => source%field%r_ptr
        IF (source%field%info%lcontained) THEN
          ! tackle 4d containers
          CALL add_fields(dst_ptr,src_ptr, source%field%info%var_ref_pos,source%field%info%ncontained)
        ELSE
          ! tackle 1d, 2d and 3d
          SELECT CASE(destination%field%info%ndims)
          CASE(3)
           !hack for sea ice variables which uses vertical level for ice class
            IF (1 == destination%field%info%used_dimensions(2)) THEN
              CALL add_fields(dst_ptr(:,:,:,1,1), &
                &             src_ptr(:,:,:,1,1), &
                &             destination%field%info%subset, &
                &             levels=1, &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            ELSE IF (ZA_OCEAN_SEDIMENT == destination%field%info%vgrid) THEN
              CALL add_fields(dst_ptr(:,:,:,1,1), &
                &             src_ptr(:,:,:,1,1), &
                &             destination%field%info%subset, &
                &             levels=destination%field%info%used_dimensions(2), &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            ELSE
              CALL add_fields(dst_ptr(:,:,:,1,1), &
                &             src_ptr(:,:,:,1,1), &
                &             destination%field%info%subset, &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            END IF
          CASE(2)
            IF (GRID_ZONAL .EQ. destination%field%info%hgrid) THEN
              CALL add_fields(dst_ptr(:,:,1,1,1), &
                &             src_ptr(:,:,1,1,1), &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            ELSE
              CALL add_fields(dst_ptr(:,:,1,1,1), &
                &             src_ptr(:,:,1,1,1), &
                &             destination%field%info%subset, &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            END IF
          CASE DEFAULT
            CALL add_fields(dst_ptr, src_ptr)
          END SELECT
        END IF

      END SUBROUTINE add_wp
      SUBROUTINE add_sp() 
        REAL(wp), POINTER :: dst_ptr(:,:,:,:,:)
        REAL(sp), POINTER :: src_ptr(:,:,:,:,:)
        dst_ptr => destination%field%r_ptr
        src_ptr => source%field%s_ptr
        IF (source%field%info%lcontained) THEN
          ! tackle 4d containers
          CALL add_fields(dst_ptr,src_ptr, source%field%info%var_ref_pos,source%field%info%ncontained)
        ELSE
          ! tackle 1d, 2d and 3d
          SELECT CASE(destination%field%info%ndims)
          CASE(3)
           !hack for sea ice variables which uses vertical level for ice class
            IF (1 == destination%field%info%used_dimensions(2)) THEN
              CALL add_fields(dst_ptr(:,:,:,1,1), &
                &             src_ptr(:,:,:,1,1), &
                &             destination%field%info%subset, &
                &             levels=1, &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            ELSE IF (ZA_OCEAN_SEDIMENT == destination%field%info%vgrid) THEN
              CALL add_fields(dst_ptr(:,:,:,1,1), &
                &             src_ptr(:,:,:,1,1), &
                &             destination%field%info%subset, &
                &             levels=destination%field%info%used_dimensions(2), &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            ELSE
              CALL add_fields(dst_ptr(:,:,:,1,1), &
                &             src_ptr(:,:,:,1,1), &
                &             destination%field%info%subset, &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            END IF
          CASE(2)
            IF (GRID_ZONAL .EQ. destination%field%info%hgrid) THEN
              CALL add_fields(dst_ptr(:,:,1,1,1), &
                &             src_ptr(:,:,1,1,1), &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            ELSE
              CALL add_fields(dst_ptr(:,:,1,1,1), &
                &             source%field%s_ptr(:,:,1,1,1), &
                &             destination%field%info%subset, &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            END IF
          CASE DEFAULT
            CALL add_fields(dst_ptr,src_ptr)
          END SELECT
        END IF
      END SUBROUTINE add_sp
  END SUBROUTINE statistics_add

  ! methods needed for update statistics each timestep
  !>
  !! implement addition for source fields to the internal accumulation fields
  !!
  SUBROUTINE statistics_add_sqr(source, destination, counter)
    TYPE(t_list_element) , INTENT(IN)    :: source
    TYPE(t_list_element) , INTENT(INOUT) :: destination
    INTEGER, INTENT(inout)               :: counter

    INTEGER :: index

    index = -1

    IF (.NOT. ASSOCIATED(source%field%r_ptr)) THEN
      CALL add_sqr_sp()
    ELSE
      CALL add_sqr_wp()
    END IF
    counter                 = counter + 1

    CONTAINS
      SUBROUTINE add_sqr_wp()
        REAL(wp), POINTER :: dst_ptr(:,:,:,:,:), src_ptr(:,:,:,:,:)
        dst_ptr => destination%field%r_ptr
        src_ptr => source%field%r_ptr
        IF (source%field%info%lcontained) THEN
          ! tackle 4d containers
          CALL add_sqr_fields(dst_ptr,src_ptr, source%field%info%var_ref_pos,source%field%info%ncontained)
        ELSE
          ! tackle 1d, 2d and 3d
          SELECT CASE(destination%field%info%ndims)
          CASE(3)
           !hack for sea ice variables which uses vertical level for ice class
            IF (1 == destination%field%info%used_dimensions(2)) THEN
              CALL add_sqr_fields(dst_ptr(:,:,:,1,1), &
                &             src_ptr(:,:,:,1,1), &
                &             destination%field%info%subset, &
                &             levels=1, &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            ELSE IF (ZA_OCEAN_SEDIMENT == destination%field%info%vgrid) THEN
              CALL add_sqr_fields(dst_ptr(:,:,:,1,1), &
                &             src_ptr(:,:,:,1,1), &
                &             destination%field%info%subset, &
                &             levels=destination%field%info%used_dimensions(2), &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            ELSE
              CALL add_sqr_fields(dst_ptr(:,:,:,1,1), &
                &             src_ptr(:,:,:,1,1), &
                &             destination%field%info%subset, &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            END IF
          CASE(2)
            IF (GRID_ZONAL .EQ. destination%field%info%hgrid) THEN
              CALL add_sqr_fields(dst_ptr(:,:,1,1,1), &
                &             src_ptr(:,:,1,1,1), &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            ELSE
              CALL add_sqr_fields(dst_ptr(:,:,1,1,1), &
                &             src_ptr(:,:,1,1,1), &
                &             destination%field%info%subset, &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            END IF
          CASE DEFAULT
            CALL add_sqr_fields(dst_ptr,src_ptr)
          END SELECT
        END IF
      END SUBROUTINE add_sqr_wp
      SUBROUTINE add_sqr_sp()
        REAL(wp), POINTER :: dst_ptr(:,:,:,:,:)
        REAL(sp), POINTER :: src_ptr(:,:,:,:,:)
        dst_ptr => destination%field%r_ptr
        src_ptr => source%field%s_ptr
        IF (source%field%info%lcontained) THEN
          ! tackle 4d containers
          CALL add_sqr_fields(dst_ptr,src_ptr, source%field%info%var_ref_pos,source%field%info%ncontained)
        ELSE
          ! tackle 1d, 2d and 3d
          SELECT CASE(destination%field%info%ndims)
          CASE(3)
           !hack for sea ice variables which uses vertical level for ice class
            IF (1 == destination%field%info%used_dimensions(2)) THEN
              CALL add_sqr_fields(dst_ptr(:,:,:,1,1), &
                &             source%field%s_ptr(:,:,:,1,1), &
                &             destination%field%info%subset, &
                &             levels=1, &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            ELSE IF (ZA_OCEAN_SEDIMENT == destination%field%info%vgrid) THEN
              CALL add_sqr_fields(dst_ptr(:,:,:,1,1), &
                &             source%field%s_ptr(:,:,:,1,1), &
                &             destination%field%info%subset, &
                &             levels=destination%field%info%used_dimensions(2), &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            ELSE
              CALL add_sqr_fields(dst_ptr(:,:,:,1,1), &
                &             source%field%s_ptr(:,:,:,1,1), &
                &             destination%field%info%subset, &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            END IF
          CASE(2)
            IF (GRID_ZONAL .EQ. destination%field%info%hgrid) THEN
              CALL add_sqr_fields(dst_ptr(:,:,1,1,1), &
                &             source%field%s_ptr(:,:,1,1,1), &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            ELSE
              CALL add_sqr_fields(dst_ptr(:,:,1,1,1), &
                &             source%field%s_ptr(:,:,1,1,1), &
                &             destination%field%info%subset, &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            END IF
          CASE DEFAULT
            CALL add_sqr_fields(dst_ptr,src_ptr)
          END SELECT
        END IF
      END SUBROUTINE add_sqr_sp
  END SUBROUTINE statistics_add_sqr

  SUBROUTINE statistics_assign(source, destination, counter)
    TYPE(t_list_element) , INTENT(IN)    :: source
    TYPE(t_list_element) , INTENT(INOUT) :: destination
    INTEGER, INTENT(inout)               :: counter

    INTEGER :: index

    index = -1

    IF (.NOT. ASSOCIATED(source%field%r_ptr)) THEN
      CALL assign_sp()
    ELSE
      CALL assign_wp()
    END IF
    counter                 = counter + 1

    CONTAINS
      SUBROUTINE assign_wp()
        REAL(wp), POINTER :: dst_ptr(:,:,:,:,:), src_ptr(:,:,:,:,:)
        dst_ptr => destination%field%r_ptr
        src_ptr => source%field%r_ptr
        IF (source%field%info%lcontained) THEN
          ! tackle 4d containers
          CALL assign_fields(dst_ptr,src_ptr, source%field%info%var_ref_pos,source%field%info%ncontained)
        ELSE
          ! tackle 1d, 2d and 3d
          SELECT CASE(destination%field%info%ndims)
          CASE(3)
           !hack for sea ice variables which uses vertical level for ice class
            IF (1 == destination%field%info%used_dimensions(2)) THEN
              CALL assign_fields(dst_ptr(:,:,:,1,1), &
                &             src_ptr(:,:,:,1,1), &
                &             destination%field%info%subset, &
                &             levels=1, &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            ELSE IF (ZA_OCEAN_SEDIMENT == destination%field%info%vgrid) THEN
              CALL assign_fields(dst_ptr(:,:,:,1,1), &
                &             src_ptr(:,:,:,1,1), &
                &             destination%field%info%subset, &
                &             levels=destination%field%info%used_dimensions(2), &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            ELSE
              CALL assign_fields(dst_ptr(:,:,:,1,1), &
                &             src_ptr(:,:,:,1,1), &
                &             destination%field%info%subset, &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            END IF
          CASE(2)
            IF (GRID_ZONAL .EQ. destination%field%info%hgrid) THEN
              CALL assign_fields(dst_ptr(:,:,1,1,1), &
                &             src_ptr(:,:,1,1,1), &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            ELSE
              CALL assign_fields(dst_ptr(:,:,1,1,1), &
                &             src_ptr(:,:,1,1,1), &
                &             destination%field%info%subset, &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            END IF
          CASE DEFAULT
            CALL assign_fields(dst_ptr,src_ptr)
          END SELECT
        END IF
      END SUBROUTINE assign_wp
      SUBROUTINE assign_sp()
        REAL(wp), POINTER :: dst_ptr(:,:,:,:,:)
        REAL(sp), POINTER :: src_ptr(:,:,:,:,:)
        dst_ptr => destination%field%r_ptr
        src_ptr => source%field%s_ptr
        IF (source%field%info%lcontained) THEN
          ! tackle 4d containers
          CALL assign_fields(dst_ptr,src_ptr, source%field%info%var_ref_pos,source%field%info%ncontained)
        ELSE
          ! tackle 1d, 2d and 3d
          SELECT CASE(destination%field%info%ndims)
          CASE(3)
           !hack for sea ice variables which uses vertical level for ice class
            IF (1 == destination%field%info%used_dimensions(2)) THEN
              CALL assign_fields(dst_ptr(:,:,:,1,1), &
                &             source%field%s_ptr(:,:,:,1,1), &
                &             destination%field%info%subset, &
                &             levels=1, &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            ELSE IF (ZA_OCEAN_SEDIMENT == destination%field%info%vgrid) THEN
              CALL assign_fields(dst_ptr(:,:,:,1,1), &
                &             source%field%s_ptr(:,:,:,1,1), &
                &             destination%field%info%subset, &
                &             levels=destination%field%info%used_dimensions(2), &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            ELSE
              CALL assign_fields(dst_ptr(:,:,:,1,1), &
                &             source%field%s_ptr(:,:,:,1,1), &
                &             destination%field%info%subset, &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            END IF
          CASE(2)
            IF (GRID_ZONAL .EQ. destination%field%info%hgrid) THEN
              CALL assign_fields(dst_ptr(:,:,1,1,1), &
                &             source%field%s_ptr(:,:,1,1,1), &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            ELSE
              CALL assign_fields(dst_ptr(:,:,1,1,1), &
                &             source%field%s_ptr(:,:,1,1,1), &
                &             destination%field%info%subset, &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            END IF
          CASE DEFAULT
            CALL assign_fields(dst_ptr,src_ptr)
          END SELECT
        END IF
      END SUBROUTINE assign_sp
  END SUBROUTINE statistics_assign

  SUBROUTINE statistics_max(source, destination, counter)
    TYPE(t_list_element) , INTENT(IN)    :: source
    TYPE(t_list_element) , INTENT(INOUT) :: destination
    INTEGER, INTENT(inout)               :: counter

    INTEGER :: index

    index = -1

    IF (.NOT. ASSOCIATED(source%field%r_ptr)) THEN
      CALL max_sp()
    ELSE
      CALL max_wp()
    END IF
    counter                 = counter + 1

    CONTAINS
      SUBROUTINE max_wp()
        REAL(wp), POINTER :: dst_ptr(:,:,:,:,:), src_ptr(:,:,:,:,:)
        dst_ptr => destination%field%r_ptr
        src_ptr => source%field%r_ptr
        IF (source%field%info%lcontained) THEN
          ! tackle 4d containers
          CALL max_fields(dst_ptr,src_ptr, source%field%info%var_ref_pos,source%field%info%ncontained)
        ELSE
          ! tackle 1d, 2d and 3d
          SELECT CASE(destination%field%info%ndims)
          CASE(3)
           !hack for sea ice variables which uses vertical level for ice class
            IF (1 == destination%field%info%used_dimensions(2)) THEN
              CALL max_fields(dst_ptr(:,:,:,1,1), &
                &             src_ptr(:,:,:,1,1), &
                &             destination%field%info%subset, &
                &             levels=1, &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            ELSE IF (ZA_OCEAN_SEDIMENT == destination%field%info%vgrid) THEN
              CALL max_fields(dst_ptr(:,:,:,1,1), &
                &             src_ptr(:,:,:,1,1), &
                &             destination%field%info%subset, &
                &             levels=destination%field%info%used_dimensions(2), &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            ELSE
              CALL max_fields(dst_ptr(:,:,:,1,1), &
                &             src_ptr(:,:,:,1,1), &
                &             destination%field%info%subset, &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            END IF
          CASE(2)
            IF (GRID_ZONAL .EQ. destination%field%info%hgrid) THEN
              CALL max_fields(dst_ptr(:,:,1,1,1), &
                &             src_ptr(:,:,1,1,1), &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            ELSE
              CALL max_fields(dst_ptr(:,:,1,1,1), &
                &             src_ptr(:,:,1,1,1), &
                &             destination%field%info%subset, &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            END IF
          CASE DEFAULT
            CALL max_fields(dst_ptr,src_ptr)
          END SELECT
        END IF
      END SUBROUTINE max_wp
      SUBROUTINE max_sp()
        REAL(wp), POINTER :: dst_ptr(:,:,:,:,:)
        REAL(sp), POINTER :: src_ptr(:,:,:,:,:)
        dst_ptr => destination%field%r_ptr
        src_ptr => source%field%s_ptr
        IF (source%field%info%lcontained) THEN
          ! tackle 4d containers
          CALL max_fields(dst_ptr,src_ptr, source%field%info%var_ref_pos,source%field%info%ncontained)
        ELSE
          ! tackle 1d, 2d and 3d
          SELECT CASE(destination%field%info%ndims)
          CASE(3)
           !hack for sea ice variables which uses vertical level for ice class
            IF (1 == destination%field%info%used_dimensions(2)) THEN
              CALL max_fields(dst_ptr(:,:,:,1,1), &
                &             source%field%s_ptr(:,:,:,1,1), &
                &             destination%field%info%subset, &
                &             levels=1, &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            ELSE IF (ZA_OCEAN_SEDIMENT == destination%field%info%vgrid) THEN
              CALL max_fields(dst_ptr(:,:,:,1,1), &
                &             source%field%s_ptr(:,:,:,1,1), &
                &             destination%field%info%subset, &
                &             levels=destination%field%info%used_dimensions(2), &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            ELSE
              CALL max_fields(dst_ptr(:,:,:,1,1), &
                &             source%field%s_ptr(:,:,:,1,1), &
                &             destination%field%info%subset, &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            END IF
          CASE(2)
            IF (GRID_ZONAL .EQ. destination%field%info%hgrid) THEN
              CALL max_fields(dst_ptr(:,:,1,1,1), &
                &             source%field%s_ptr(:,:,1,1,1), &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            ELSE
              CALL max_fields(dst_ptr(:,:,1,1,1), &
                &             source%field%s_ptr(:,:,1,1,1), &
                &             destination%field%info%subset, &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            END IF
          CASE DEFAULT
            CALL max_fields(dst_ptr,src_ptr)
          END SELECT
        END IF
      END SUBROUTINE max_sp
  END SUBROUTINE statistics_max

  SUBROUTINE statistics_min(source, destination, counter)
    TYPE(t_list_element) , INTENT(IN)    :: source
    TYPE(t_list_element) , INTENT(INOUT) :: destination
    INTEGER, INTENT(inout)               :: counter

    INTEGER :: index

    index = -1
    IF (.NOT. ASSOCIATED(source%field%r_ptr)) THEN
      CALL min_sp()
    ELSE
      CALL min_wp()
    END IF
    counter                 = counter + 1

    CONTAINS
      SUBROUTINE min_wp()
        REAL(wp), POINTER :: dst_ptr(:,:,:,:,:), src_ptr(:,:,:,:,:)
        dst_ptr => destination%field%r_ptr
        src_ptr => source%field%r_ptr
        IF (source%field%info%lcontained) THEN
          ! tackle 4d containers
          CALL min_fields(dst_ptr,src_ptr, source%field%info%var_ref_pos,source%field%info%ncontained)
        ELSE
          ! tackle 1d, 2d and 3d
          SELECT CASE(destination%field%info%ndims)
          CASE(3)
           !hack for sea ice variables which uses vertical level for ice class
            IF (1 == destination%field%info%used_dimensions(2)) THEN
              CALL min_fields(dst_ptr(:,:,:,1,1), &
                &             src_ptr(:,:,:,1,1), &
                &             destination%field%info%subset, &
                &             levels=1, &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            ELSE IF (ZA_OCEAN_SEDIMENT == destination%field%info%vgrid) THEN
              CALL min_fields(dst_ptr(:,:,:,1,1), &
                &             src_ptr(:,:,:,1,1), &
                &             destination%field%info%subset, &
                &             levels=destination%field%info%used_dimensions(2), &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            ELSE
              CALL min_fields(dst_ptr(:,:,:,1,1), &
                &             src_ptr(:,:,:,1,1), &
                &             destination%field%info%subset, &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            END IF
          CASE(2)
            IF (GRID_ZONAL .EQ. destination%field%info%hgrid) THEN
              CALL min_fields(dst_ptr(:,:,1,1,1), &
                &             src_ptr(:,:,1,1,1), &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            ELSE
              CALL min_fields(dst_ptr(:,:,1,1,1), &
                &             src_ptr(:,:,1,1,1), &
                &             destination%field%info%subset, &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            END IF
          CASE DEFAULT
            CALL min_fields(dst_ptr,src_ptr)
          END SELECT
        END IF
      END SUBROUTINE min_wp
      SUBROUTINE min_sp()
        REAL(wp), POINTER :: dst_ptr(:,:,:,:,:)
        REAL(sp), POINTER :: src_ptr(:,:,:,:,:)
        dst_ptr => destination%field%r_ptr
        src_ptr => source%field%s_ptr
        IF (source%field%info%lcontained) THEN
          ! tackle 4d containers
          CALL min_fields(dst_ptr,src_ptr, source%field%info%var_ref_pos,source%field%info%ncontained)
        ELSE
          ! tackle 1d, 2d and 3d
          SELECT CASE(destination%field%info%ndims)
          CASE(3)
           !hack for sea ice variables which uses vertical level for ice class
            IF (1 == destination%field%info%used_dimensions(2)) THEN
              CALL min_fields(dst_ptr(:,:,:,1,1), &
                &             source%field%s_ptr(:,:,:,1,1), &
                &             destination%field%info%subset, &
                &             levels=1, &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            ELSE IF (ZA_OCEAN_SEDIMENT == destination%field%info%vgrid) THEN
              CALL min_fields(dst_ptr(:,:,:,1,1), &
                &             source%field%s_ptr(:,:,:,1,1), &
                &             destination%field%info%subset, &
                &             levels=destination%field%info%used_dimensions(2), &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            ELSE
              CALL min_fields(dst_ptr(:,:,:,1,1), &
                &             source%field%s_ptr(:,:,:,1,1), &
                &             destination%field%info%subset, &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            END IF
          CASE(2)
            IF (GRID_ZONAL .EQ. destination%field%info%hgrid) THEN
              CALL min_fields(dst_ptr(:,:,1,1,1), &
                &             source%field%s_ptr(:,:,1,1,1), &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            ELSE
              CALL min_fields(dst_ptr(:,:,1,1,1), &
                &             source%field%s_ptr(:,:,1,1,1), &
                &             destination%field%info%subset, &
                &             has_missvals=destination%field%info%lmiss, &
                &             missval=destination%field%info%missval%rval)
            END IF
          CASE DEFAULT
            CALL min_fields(dst_ptr,src_ptr)
          END SELECT
        END IF
      END SUBROUTINE min_sp
  END SUBROUTINE statistics_min
  !>
  !! Execute the accumulation forall internal variables and compute mean values
  !! if the corresponding event is active
  !!
  SUBROUTINE update_statistics()
    CALL update_mvstream(meanMap, meanEvents, meanEventsActivity, meanPrognostics, meanPrognosticPointers,  meanVarCounter,MEAN)
    CALL update_mvstream(maxMap, maxEvents, maxEventsActivity, maxPrognostics, maxPrognosticPointers,  maxVarCounter,MAX)
    CALL update_mvstream(minMap, minEvents, minEventsActivity, minPrognostics, minPrognosticPointers,  minVarCounter,MIN)
    CALL update_mvstream(squareMap, &
      &                  squareEvents, &
      &                  squareEventsActivity, &
      &                  squarePrognostics, &
      &                  squarePrognosticPointers, &
      &                  squareVarCounter, &
      &                  SQUARE)
  END SUBROUTINE update_statistics

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
    CLASS(*),pointer :: varListForMeanEvent, statisticEventKey
    CLASS(*),pointer :: sourceVariable, destinationVariable, sourceVariable4Prognostics
    CLASS(*),pointer :: counter, statisticEvent, myItem
    CLASS(*),pointer :: eventActive
    TYPE(t_list_element), pointer :: source, destination
    TYPE(vector_iterator) :: statisticMapIterator, statisticEventIterator
    CLASS(*), pointer :: mvstream_pair

    TYPE(datetime), POINTER :: mtime_date
    LOGICAL :: isactive
    INTEGER :: varcounter
    INTEGER :: timelevel
#ifdef DEBUG_MVSTREAM
    character(len=max_datetime_str_len) :: datetimestring
#endif

    CHARACTER(LEN=*), PARAMETER :: routine =  modname//"::update_statistics"
    REAL(wp), POINTER :: dst_ptr(:,:,:,:,:)

#ifdef DEBUG_MVSTREAM
    IF (my_process_is_stdio()) CALL print_routine(routine,'start',stderr=.true.)
#endif

    statisticMapIterator   = statisticMap%iter()
    statisticEventIterator = statisticEvents%iter()

    ! Check events first {{{
    ! this is necessary because of mtime internals
    isactive   = .false.
    mtime_date => newDatetime(time_config%tc_current_date)

#ifdef DEBUG_MVSTREAM
    CALL datetimeToString(mtime_date, datetimestring)
    IF (my_process_is_stdio()) CALL print_summary('Current mtime timestamp:'//trim(datetimestring),stderr=.true.)
#endif

    ! Save results for (not so much) later
    DO WHILE (statisticEventIterator%next(myItem))
      SELECT TYPE (myItem)
      TYPE IS (map_item)
        statisticEventKey => myItem%key
        statisticEvent    => myItem%value
        SELECT TYPE (statisticEvent)
        TYPE IS (t_event_wrapper)
          isactive = LOGICAL(isCurrentEventActive(statisticEvent%this,mtime_date))
        END SELECT
        CALL statisticEventsActivity%add(statisticEventKey,isactive)
      END SELECT
    END DO

#ifdef DEBUG_MVSTREAM
    IF (my_process_is_stdio()) CALL print_error(statisticEventsActivity%to_string(),stderr=.true.)
#endif
    ! }}}

    DO WHILE (statisticMapIterator%next(myItem))

      SELECT TYPE (myItem)
      TYPE IS (map_item)

        statisticEventKey        => myItem%key
        varListForMeanEvent => myItem%value

#ifdef DEBUG_MVSTREAM
        IF (my_process_is_stdio()) CALL print_summary(object_pointer_string(statisticEventKey)//"PERFORM ACCU",stderr=.true.) !TODO
#endif

        SELECT TYPE(varListForMeanEvent)
        CLASS IS (vector_ref)
          DO element_counter=1,varListForMeanEvent%length(),2

#ifdef DEBUG_MVSTREAM
          IF (my_process_is_stdio()) CALL print_routine("update_statistics",object_string(element_counter),stderr=.true.)
#endif
            sourceVariable      => varListForMeanEvent%at(element_counter)
            destinationVariable => varListForMeanEvent%at(element_counter+1)
            IF (ASSOCIATED(sourceVariable) .and. ASSOCIATED(destinationVariable)) THEN
              SELECT TYPE (sourceVariable)
              TYPE IS (t_list_element)
                SELECT TYPE (destinationVariable)
                TYPE IS (t_list_element)
                  ! check for prognostics, pointer must be shifted according to given timelevelIndex {{{
                  IF (is_statsPrognosticVariable(destinationVariable,statisticPrognostics)) THEN
                    ! find the correct pointer:
                    ! check f reduced calling freq. variables are used for
                    ! output? if true, theses variables should be used for
                    ! accumulation, too

#ifdef DEBUG_MVSTREAM
                    IF (my_process_is_stdio()) &
                      & CALL print_error("destination var IS     prognostic:"//TRIM(destinationVariable%field%info%name),&
                      & stderr=.true.)
#endif

                    timelevel =  metainfo_get_timelevel(destinationVariable%field%info, &
                        & destinationVariable%field%info%dom)

                    source    => get_prognostics_source_pointer(destinationVariable, timelevel, statisticPrognosticPointers)

                  ELSE
                  ! }}}
                    source    => sourceVariable
#ifdef DEBUG_MVSTREAM
                    IF (my_process_is_stdio()) &
                      & CALL print_error("destination var IS NOT prognostic:"//TRIM(destinationVariable%field%info%name),&
                      & stderr=.true.)
#endif
                  END IF
                  destination => destinationVariable
                  counter     => statisticVarCounter%get(destination%field%info%name)
                  SELECT TYPE (counter)
                  TYPE IS (INTEGER)

#ifdef DEBUG_MVSTREAM
                    IF ( my_process_is_stdio() ) CALL print_summary('sourceName : '//trim(source%field%info%name),&
                        &stderr=.true.)
                    IF ( my_process_is_stdio() ) CALL print_summary('destName   : '//trim(destination%field%info%name),&
                        & stderr=.true.)
                    IF ( my_process_is_stdio() ) CALL print_summary('destNameOut: '//trim(destination%field%info%cf%short_name),&
                        &stderr=.true.)
                    IF ( my_process_is_stdio() ) write (0,*)'old counter: ',counter
#endif

                    ! update of internal fields according to the operator {{{
                    varcounter = counter !TODO work around for integer pointer, ugly
                    IF (MEAN.EQ.operation) THEN
                      CALL statistics_add(source, destination, varcounter)

                    ELSE IF (MAX.EQ.operation) THEN
                      IF (0 .EQ. varcounter) THEN
                        CALL statistics_assign(source, destination, varcounter)
                      ELSE
                        CALL statistics_max(source, destination, varcounter)
                      END IF

                    ELSE IF (MIN.EQ.operation) THEN
                      IF (0 .EQ. varcounter) THEN
                        CALL statistics_assign(source, destination, varcounter)
                      ELSE
                        CALL statistics_min(source, destination, varcounter)
                      END IF

                    ELSE IF (SQUARE.EQ.operation) THEN
                      CALL statistics_add_sqr(source, destination, varcounter)

                    ELSE
                      CALL finish('update_mvstream','Found unknown operation:'//trim(operation))

                    END IF
                    counter = varcounter
                    ! }}}
#ifdef DEBUG_MVSTREAM
                    IF ( my_process_is_stdio() )  write (0,*)'new counter: ',counter
#endif

                  END SELECT

! Subsequently the accumulated arrays may be updated on the HOST
!$ACC WAIT

                  ! MEAN VALUE COMPUTAION {{{
                  ! check if the field will be written to disk this timestep {{{
                  eventActive => statisticEventsActivity%get(statisticEventKey)
                  SELECT TYPE (eventActive)
                  TYPE IS (LOGICAL)
                    isactive = eventActive
                    IF ( isactive ) THEN

#ifdef DEBUG_MVSTREAM
                    IF (my_process_is_stdio()) &
                        & CALL print_summary(" --> PERFORM MEAN VALUE COMP for"//trim(destination%field%info%name),stderr=.true.)
#endif

                      counter => statisticVarCounter%get(destination%field%info%name)
                      SELECT TYPE(counter)
                      TYPE IS (INTEGER)
#ifdef DEBUG_MVSTREAM
                      IF ( my_process_is_stdio() ) write (0,*)' ------> MEAN VALUE counter:',counter
#endif
                      IF (MEAN.EQ.operation .OR. SQUARE.EQ.operation) THEN
                        dst_ptr => destination%field%r_ptr
                        varcounter = counter   ! PGI 19.9 cannot properly pass this scalar to Kernel
                        !$ACC KERNELS PRESENT(dst_ptr) IF (i_am_accel_node)
                        dst_ptr = dst_ptr / REAL(varcounter,wp)
                        !$ACC END KERNELS
                      END IF
                      counter = 0
                      END SELECT
                    END IF
                  END SELECT
                  ! }}}
                END SELECT
                ! ! }}}
              CLASS DEFAULT
                CALL finish('update_statistics','Found unknown source variable type')
              END SELECT
            ELSE
              CALL finish(routine,'source or destination variable cannot be found')
            END IF
          END DO
        END SELECT
      END SELECT
    END DO

#ifdef DEBUG_MVSTREAM
    IF (my_process_is_stdio()) CALL print_routine(routine,'finish')
#endif

  END SUBROUTINE update_mvstream

  !>
  !! reset internal fields to zero, if the corresponding event is active
  !!
  ! min/max do not need extra reset due to JEFs implementation idea
  SUBROUTINE reset_statistics
    CALL reset_mvstream(meanMap, meanEventsActivity)
    CALL reset_mvstream(squareMap, squareEventsActivity)
  END SUBROUTINE reset_statistics


  SUBROUTINE reset_mvstream(statisticMap, statisticEventsActivity)
    TYPE(map) :: statisticMap, statisticEventsActivity

    INTEGER :: element_counter
    TYPE(t_list_element), pointer :: destination
    CLASS(*),pointer :: varListForStatisticEvent, statisticEventKey
    CLASS(*),pointer :: destinationVariable
    CLASS(*),pointer :: statisticEvent, myItem
    CLASS(*), pointer :: eventActive
    TYPE(vector_iterator) :: statisticMapIterator

    CHARACTER(LEN=*), PARAMETER :: routine =  modname//"::reset_statistics"
    REAL(wp), POINTER :: dst_ptr(:,:,:,:,:)

#ifdef DEBUG_MVSTREAM
    IF (my_process_is_stdio()) CALL print_routine(routine,'start')
#endif

    statisticMapIterator = statisticMap%iter()

    DO WHILE (statisticMapIterator%next(myItem))

#ifdef DEBUG_MVSTREAM
      IF (my_process_is_stdio()) CALL print_error("started event loop",stderr=.true.)
#endif
      SELECT TYPE (myItem)
      TYPE IS (map_item)

        statisticEventKey   => myItem%key
        varListForStatisticEvent => myItem%value

        SELECT TYPE(varListForStatisticEvent)
        CLASS IS (vector_ref)

          DO element_counter=1,varListForStatisticEvent%length(),2 !start at 2 because the event is at index 1

            destinationVariable => varListForStatisticEvent%at(element_counter+1)

#ifdef DEBUG_MVSTREAM
            IF (my_process_is_stdio()) CALL print_error("got destinationVariable",stderr=.true.)
#endif

            IF (ASSOCIATED(destinationVariable)) THEN

#ifdef DEBUG_MVSTREAM
            IF (my_process_is_stdio()) CALL print_error("    destinationVariable is ASSOCIATED",stderr=.true.)
#endif

              SELECT TYPE (destinationVariable)
              TYPE IS (t_list_element)
                  destination => destinationVariable

#ifdef DEBUG_MVSTREAM
                  IF (my_process_is_stdio()) CALL print_error("    destinationVariable is t_list_element",stderr=.true.)
#endif

                  eventActive => statisticEventsActivity%get(statisticEventKey)

#ifdef DEBUG_MVSTREAM
                  IF (.NOT.ASSOCIATED(eventActive)) THEN
                    IF (my_process_is_stdio()) CALL print_error("       eventActive not ASSOCIATED",stderr=.true.)
                  END IF
#endif

                  SELECT TYPE (eventActive)
                  TYPE IS (LOGICAL)

#ifdef DEBUG_MVSTREAM
                    IF (my_process_is_stdio()) CALL print_error("       eventActive is LOGICAL",stderr=.true.)
#endif

                    IF ( LOGICAL(eventActive) ) THEN

#ifdef DEBUG_MVSTREAM
                      IF (my_process_is_stdio()) CALL print_error("       eventActive is true",stderr=.true.)
                      IF (my_process_is_stdio()) CALL print_error(object_string(statisticEventKey)//' : --> PERFORM RESET',&
                          & stderr=.true.)
IF (my_process_is_stdio()) CALL print_error(object_string(statisticEventKey)//' : --> '//trim(destination%field%info%name),&
    & stderr=.true.)
#endif
                      dst_ptr => destination%field%r_ptr
                      !$ACC KERNELS PRESENT(dst_ptr) IF (i_am_accel_node)
                      dst_ptr = 0.0_wp ! take the neutral element of addition
                      !$ACC END KERNELS
#ifdef DEBUG_MVSTREAM
                    ELSE
                      IF (my_process_is_stdio()) CALL print_error("       eventActive is false",stderr=.true.)
#endif
                    END IF
                  CLASS DEFAULT
#ifdef DEBUG_MVSTREAM
                      IF (my_process_is_stdio()) CALL print_error("       eventActive has wrong type",stderr=.true.)
#endif
                  END SELECT
              CLASS DEFAULT
#ifdef DEBUG_MVSTREAM
                  IF (my_process_is_stdio()) CALL print_error("     destinationVariable is not t_list_element",stderr=.true.)
#endif
              END SELECT
            ELSE
#ifdef DEBUG_MVSTREAM
              IF (my_process_is_stdio()) CALL print_error(routine//TRIM(": cannot find destination variable!"),stderr=.true.)
#endif
            END IF
          END DO
        END SELECT
      END SELECT
    END DO

#ifdef DEBUG_MVSTREAM
    IF (my_process_is_stdio()) CALL print_routine(routine,'finish',stderr=.true.)
#endif
  END SUBROUTINE

  SUBROUTINE reset_mean
    INTEGER :: element_counter
    TYPE(t_list_element), pointer :: destination
    CLASS(*),pointer :: varListForMeanEvent, meanEventKey
    CLASS(*),pointer :: destinationVariable
    CLASS(*),pointer :: meanEvent, myItem
    CLASS(*), pointer :: eventActive
    TYPE(vector_iterator) :: meanMapIterator

    CHARACTER(LEN=*), PARAMETER :: routine =  modname//"::reset_statistics"
    REAL(wp), POINTER :: dst_ptr(:,:,:,:,:)

#ifdef DEBUG_MVSTREAM
    IF (my_process_is_stdio()) CALL print_routine(routine,'start')
#endif

    meanMapIterator   = meanMap%iter()

    DO WHILE (meanMapIterator%next(myItem))

#ifdef DEBUG_MVSTREAM
      IF (my_process_is_stdio()) CALL print_error("started event loop",stderr=.true.)
#endif
      SELECT TYPE (myItem)
      TYPE IS (map_item)

        meanEventKey        => myItem%key
        varListForMeanEvent => myItem%value

        SELECT TYPE(varListForMeanEvent)
        CLASS IS (vector_ref)

          DO element_counter=1,varListForMeanEvent%length(),2 !start at 2 because the event is at index 1

            destinationVariable => varListForMeanEvent%at(element_counter+1)

#ifdef DEBUG_MVSTREAM
            IF (my_process_is_stdio()) CALL print_error("got destinationVariable",stderr=.true.)
#endif

            IF (ASSOCIATED(destinationVariable)) THEN

#ifdef DEBUG_MVSTREAM
            IF (my_process_is_stdio()) CALL print_error("    destinationVariable is ASSOCIATED",stderr=.true.)
#endif

              SELECT TYPE (destinationVariable)
              TYPE IS (t_list_element)
                  destination => destinationVariable

#ifdef DEBUG_MVSTREAM
                  IF (my_process_is_stdio()) CALL print_error("    destinationVariable is t_list_element",stderr=.true.)
#endif

                  eventActive => meanEventsActivity%get(meanEventKey)

#ifdef DEBUG_MVSTREAM
                  IF (.NOT. ASSOCIATED(eventActive)) THEN
                    IF (my_process_is_stdio()) CALL print_error("       eventActive not ASSOCIATED",stderr=.true.)
                  END IF
#endif

                  SELECT TYPE (eventActive)
                  TYPE IS (LOGICAL)

#ifdef DEBUG_MVSTREAM
                    IF (my_process_is_stdio()) CALL print_error("       eventActive is LOGICAL",stderr=.true.)
#endif

                    IF ( LOGICAL(eventActive) ) THEN

#ifdef DEBUG_MVSTREAM
                      IF (my_process_is_stdio()) CALL print_error("       eventActive is true",stderr=.true.)
                      IF (my_process_is_stdio()) CALL print_error(object_string(meanEventKey)//' : --> PERFORM RESET',&
                          & stderr=.true.)
IF (my_process_is_stdio()) CALL print_error(object_string(meanEventKey)//' : --> '//trim(destination%field%info%name),&
    & stderr=.true.)
#endif

                      dst_ptr => destination%field%r_ptr
                      !$ACC KERNELS PRESENT(dst_ptr) IF (i_am_accel_node)
                      dst_ptr = 0.0_wp ! take the neutral element of addition
                      !$ACC END KERNELS
#ifdef DEBUG_MVSTREAM
                    ELSE
                      IF (my_process_is_stdio()) CALL print_error("       eventActive is false",stderr=.true.)
#endif
                    END IF
                  CLASS DEFAULT
#ifdef DEBUG_MVSTREAM
                      IF (my_process_is_stdio()) CALL print_error("       eventActive has wrong type",stderr=.true.)
#endif
                  END SELECT
              CLASS DEFAULT
#ifdef DEBUG_MVSTREAM
                  IF (my_process_is_stdio()) CALL print_error("     destinationVariable is not t_list_element",stderr=.true.)
#endif
              END SELECT
            ELSE
#ifdef DEBUG_MVSTREAM
              IF (my_process_is_stdio()) CALL print_error(routine//TRIM(": cannot find destination variable!"),stderr=.true.)
#endif
            END IF
          END DO
        END SELECT
      END SELECT
    END DO

#ifdef DEBUG_MVSTREAM
    IF (my_process_is_stdio()) CALL print_routine(routine,'finish',stderr=.true.)
#endif

  END SUBROUTINE

END MODULE mo_derived_variable_handling

! vim:tw=0
