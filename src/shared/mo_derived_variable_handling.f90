!>
!! Routines for handling proxy variables e.g. accumulation buffers
!!
MODULE mo_derived_variable_handling

  USE self_object
  USE self_vector_ref
  USE self_vector
  USE self_map
  USE self_assert

  USE mo_kind, ONLY: wp
  USE mo_model_domain, ONLY: t_patch
  USE mo_io_config, ONLY: lnetcdf_flt64_output
  USE mo_dynamics_config, ONLY: nnow, nnew, nold
  USE mo_statistics, ONLY: add_fields
  USE mo_var_metadata_types, ONLY: VARNAME_LEN
  USE mo_impl_constants, ONLY: vname_len, SUCCESS, max_char_length, TLEV_NNOW, TLEV_NNEW,    &
    &                          TLEV_NNOW_RCF, TLEV_NNEW_RCF, REAL_T
  USE mo_cdi_constants, ONLY: GRID_UNSTRUCTURED_CELL, GRID_UNSTRUCTURED_EDGE, GRID_UNSTRUCTURED_VERT
  USE mo_name_list_output_types, ONLY: t_output_name_list
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
  USE mo_cdi,                 ONLY: DATATYPE_FLT32, DATATYPE_FLT64, GRID_LONLAT, GRID_ZONAL

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_derived_variable_handling'

  TYPE(map), SAVE    :: meanMap, meanEvents, meanEventsActivity, meanVarCounter, meanPrognosticPointers
  TYPE(vector), SAVE :: meanPrognostics
  TYPE(t_var_list)   :: mean_stream_list

  PUBLIC :: init_mean_stream
  PUBLIC :: finish_mean_stream
  PUBLIC :: perform_accumulation
  PUBLIC :: reset_accumulation
  PUBLIC :: process_mean_stream

  TYPE :: t_accumulation_pair
    TYPE(t_list_element), POINTER :: source
    TYPE(t_list_element), POINTER :: destination
    INTEGER                       :: counter = 0
  END TYPE t_accumulation_pair

  !>
  ! wrapper for bind_c event type - needed because bind_c types are not allowed
  ! in "select type" - this is kind of a hack
  type :: t_event_wrapper
    type(event), pointer :: this
  end type

CONTAINS

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
      TYPE is (t_accumulation_pair)
        PRINT *,'t_accumulation_pair:source     :',trim(my_buffer%source%field%info%name)
        PRINT *,'t_accumulation_pair:destination:',trim(my_buffer%destination%field%info%name)
      CLASS default
        PRINT *,' default class print  :'
        print *,object_string (my_buffer)
      END SELECT
    END DO
  END SUBROUTINE var_print

  !>
  !! function to check if a varable is prognostic and saved for accumulation
  !!
  logical function is_meanPrognosticVariable(variable)
    type(t_list_element), intent(in) :: variable
    is_meanPrognosticVariable = meanPrognostics%includes(variable%field%info%cf%short_name)
  end function is_meanPrognosticVariable

  !>
  !! return pointer to prognostic accumulation copy for given timelevel
  !!
  function get_prognostics_source_pointer(destinationVariable, timelevelIndex) result(sourcePointer)
    type(t_list_element), pointer :: sourcePointer
    type(t_list_element)          :: destinationVariable

    CHARACTER(LEN=*), PARAMETER :: routine =  modname//"::get_prognostics_source_pointer"
    INTEGER , INTENT(IN) :: timelevelIndex

    class(*), pointer :: dummy

    dummy => meanPrognosticPointers%get( &
      & get_varname_with_timelevel( &
      &   destinationVariable%field%info%cf%short_name, &
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
  SUBROUTINE meanStreamCrossCheck(p_onl)
    TYPE(t_output_name_list) :: p_onl

    CHARACTER(LEN=*), PARAMETER :: routine =  modname//"::meanStreamCrossCheck"
    LOGICAL :: abort

    abort = .FALSE.

    IF (MAXVAL(p_onl%dom) > 1) abort = .TRUE.

    IF (  p_onl%stream_partitions_ml > 1 .OR. &
      &   p_onl%stream_partitions_pl > 1 .OR. &
      &   p_onl%stream_partitions_hl > 1 .OR. &
      &   p_onl%stream_partitions_il > 1 ) abort = .TRUE.

    IF (abort) THEN
      call finish(routine,"meanValues are only supported on global domain 1 and without stream partitioning!")
    END IF
  END SUBROUTINE meanStreamCrossCheck

  !>
  !! Create a variable list
  !!
  SUBROUTINE init_mean_stream(patch_2d)
    TYPE(t_patch), TARGET, INTENT(in) :: patch_2d
    
    CHARACTER(LEN=max_char_length) :: listname
    
    ! main map for automatical mean value computation:
    ! key: eventString
    ! value: self-vector of t_list_elements objects that belong to that event
    meanMap            = map(verbose=.false.)

    ! map for holding the events itself
    ! since the events cannot be saved a keys for a map because they are only
    ! C-pointers the same string representation of events is used as for
    ! "meanMap"
    meanEvents         = map(verbose=.false.)
    meanEventsActivity = map(verbose=.false.)
    meanVarCounter     = map(verbose=.false.)
    meanPrognostics    = vector(verbose=.false.)
    meanPrognosticPointers    = map(verbose=.false.)

    listname = 'mean_stream_list'
    CALL new_var_list(mean_stream_list, listname, patch_id=patch_2d%id)
    CALL default_var_list_settings( mean_stream_list,lrestart=.FALSE.,loutput=.TRUE.)
  END SUBROUTINE init_mean_stream

  !>
  !! Delete internal mean value fields
  !!
  SUBROUTINE finish_mean_stream()

#ifdef DEBUG_MVSTREAM
    IF (my_process_is_stdio()) CALL print_summary('destruct mean stream variables')
#endif

    CALL delete_var_list(mean_stream_list)
  END SUBROUTINE finish_mean_stream

  integer function output_varlist_length(in_varlist)
    CHARACTER(LEN=vname_len), intent(in) :: in_varlist(:)

    output_varlist_length = 0
    DO
      IF (in_varlist(output_varlist_length+1) == ' ') EXIT
      output_varlist_length = output_varlist_length + 1
    END DO
  end function output_varlist_length

  !>
  !! Go through the output namelists and create events and accumulation fields if needed
  !!
  SUBROUTINE process_mean_stream(p_onl,i_typ, sim_step_info, patch_2d)
    TYPE (t_output_name_list), target  :: p_onl
    INTEGER                            :: i_typ
    TYPE (t_sim_step_info), INTENT(IN) :: sim_step_info
    TYPE(t_patch), INTENT(IN)          :: patch_2d

    CHARACTER(LEN=vname_len), POINTER :: in_varlist(:)
    INTEGER :: ntotal_vars, output_variables,i,ierrstat, dataType
    INTEGER :: timelevel, timelevels(3)
    type(vector_ref) :: meanVariables, prognosticVariables
    CHARACTER(LEN=100) :: eventKey
    type(t_event_wrapper) :: event_wrapper
    class(*), pointer :: myBuffer
    TYPE(t_list_element), POINTER :: src_element, dest_element
    CHARACTER(LEN=VARNAME_LEN), ALLOCATABLE :: varlist(:)
    CHARACTER(LEN=VARNAME_LEN) :: dest_element_name
    LOGICAL :: foundPrognostic
    CHARACTER(LEN=*), PARAMETER :: routine =  modname//"::process_mean_stream"

#ifdef DEBUG_MVSTREAM
    if (my_process_is_stdio()) call print_routine(routine,'start')
#endif

    IF ("mean" .EQ. TRIM(p_onl%operation)) THEN

#ifdef DEBUG_MVSTREAM
      if (my_process_is_stdio()) call print_routine(routine,'found "mean" operation')
#endif
      call meanStreamCrossCheck(p_onl)

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
      IF ( meanMap%has_key(eventKey) ) THEN
        myBuffer => meanMap%get(eventKey)
        select type (myBuffer)
        type is (vector_ref)
          ! use exiting list
          meanVariables = myBuffer
        end select
      ELSE
        ! create new variable list
        meanVariables = vector_ref()

        ! wrap c-pointer for events into a fortran type - look weird, but does the trick
        event_wrapper = t_event_wrapper(this=newEvent(eventKey, &
          &                 p_onl%output_start(1), &
          &                 p_onl%output_start(1), &
          &                 p_onl%output_end(1), &
          &                 p_onl%output_interval(1) &
          &                ))

        ! collect the new event
        call meanEvents%add(eventKey,event_wrapper)
      END IF
      ! initialize all events as in-active
      call meanEventsActivity%add(eventKey,.false.)

      ! create adhoc copies of all variables for later accumulation
      DO i=1, output_variables
        ! collect data variables only
        ! variables names like 'grid:clon' which should be excluded
        IF ( INDEX(varlist(i),':') > 0) CYCLE
 
        ! check for already created meanStream variable (maybe from another output_nml with the same output_interval)
        ! names consist of original spot-value names PLUS event information (start + interval of output)
        ! TODO: unify with eventKey definition if possible
        dest_element_name = get_accumulation_varname(varlist(i),p_onl)
        dest_element => find_list_element(mean_stream_list, trim(dest_element_name))
#ifdef DEBUG_MVSTREAM
        if (my_process_is_stdio()) call print_summary('destination variable NAME:'//TRIM(dest_element_name),stderr=.true.)
#endif
        IF (.not. ASSOCIATED(dest_element) ) THEN !not found -->> create a new on
          ! find existing source variable on all possible ICON grids with the identical name
          src_element => find_element ( TRIM(varlist(i)),opt_hgrid=GRID_UNSTRUCTURED_CELL,opt_caseInsensitive=.true.)
          IF (.not. ASSOCIATED(src_element) ) &
              & src_element => find_element ( TRIM(varlist(i)),opt_hgrid=GRID_UNSTRUCTURED_EDGE,opt_caseInsensitive=.true.)
          IF (.not. ASSOCIATED(src_element) ) &
              & src_element => find_element ( TRIM(varlist(i)),opt_hgrid=GRID_UNSTRUCTURED_VERT,opt_caseInsensitive=.true.)
          IF (.not. ASSOCIATED(src_element) ) &
              & src_element => find_element ( TRIM(varlist(i)),opt_hgrid=GRID_LONLAT,opt_caseInsensitive=.true.)
          IF (.not. ASSOCIATED(src_element) ) &
              & src_element => find_element ( TRIM(varlist(i)),opt_hgrid=GRID_ZONAL,opt_caseInsensitive=.true.)

          ! if not found: maybe it is a prognostic variable, so it has the
          ! time-level in its name
          ! ATTENTION: this is only a placeholder, because it catches the first match
          ! the correct pointer is the one with nnew(), but its value changes each timestep
          IF (.not. ASSOCIATED (src_element)) THEN
            foundPrognostic = .false.
            timelevels = (/nold(1),nnow(1),nnew(1)/)
            do timelevel=1,3

#ifdef DEBUG_MVSTREAM
              if (my_process_is_stdio()) call print_error(get_varname_with_timelevel(varlist(i),timelevels(timelevel)))
#endif
              src_element => find_element(get_varname_with_timelevel(varlist(i),timelevels(timelevel)))
              if ( ASSOCIATED(src_element) ) then
if (my_process_is_stdio()) write(0,*)'found prognostic:',TRIM(get_varname_with_timelevel(varlist(i),timelevels(timelevel)))
                if ( .not. foundPrognostic ) then
                  ! save the name of the original output variable if a prognosting version was found
                  call meanPrognostics%add(varlist(i))
                  foundPrognostic = .true.
                  if (my_process_is_stdio()) call print_error('meanPrognostics%add():'//varlist(i))
                end if
                ! save the the pointers for all time levels of a prognostic variable
                call meanPrognosticPointers%add(get_varname_with_timelevel(varlist(i),timelevels(timelevel)),src_element)
if (my_process_is_stdio()) write(0,*)'IS pROGNOSTIC:',TRIM(varlist(i))
              end if
            end do
          END IF
          ! in case nothing appropriate could be found, throw an error
          IF (.not. ASSOCIATED (src_element)) THEN
            call finish(routine,'Could not find source variable:'//TRIM(varlist(i)))
          END IF

          ! add new mean variable, copy the meta-data from the existing variable
          ! 1. copy the source variable to destination pointer
          dest_element => copy_var_to_list(mean_stream_list,dest_element_name,src_element, patch_2d)

#ifdef DEBUG_MVSTREAM
          if (my_process_is_stdio()) CALL print_summary('src(name)     :|'//trim(src_element%field%info%name)//'|',stderr=.true.)
          if (my_process_is_stdio()) CALL print_summary('src(grid)     :|'//int2string(src_element%field%info%hgrid)//'|', &
              & stderr=.true.)
          if (my_process_is_stdio()) CALL print_summary('varlist(name) :|'//trim(in_varlist(i))//'|',stderr=.true.)
          if (my_process_is_stdio()) CALL print_summary('new name      :|'//trim(dest_element_name)//'|',stderr=.true.)
          if (my_process_is_stdio()) CALL print_summary('new grid      :|'//int2string(dest_element%field%info%hgrid)//'|',&
              & stderr=.true.)
#endif

          ! set output to double precission if necessary
          dest_element%field%info%cf%datatype = MERGE(DATATYPE_FLT64, DATATYPE_FLT32, lnetcdf_flt64_output)

          ! 2. update the nc-shortname to internal name of the source variable
          dest_element%field%info%cf%short_name = get_var_name(src_element%field)

          ! Collect variable pointers for source and destination in the same list {{{
          CALL meanVariables%add(src_element)
          CALL meanVariables%add(dest_element)
          ! collect the counter for each destination in a separate list
          CALL meanVarCounter%add(dest_element%field%info%name,0)

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

      call meanMap%add(eventKey,meanVariables)
      if (my_process_is_stdio()) call print_error(routine//": meanPrognostics%to_string()")
      if (my_process_is_stdio()) call print_error(meanPrognostics%to_string())
    ELSE
#ifdef DEBUG_MVSTREAM
      if (my_process_is_stdio()) call print_routine(routine,'NO "mean" operation found')
#endif
    END IF

#ifdef DEBUG_MVSTREAM
    if (my_process_is_stdio()) then
        call var_print(meanVariables)
      call print_routine(routine, 'meanMap%to_string()')
      call print_routine(routine, meanMap%to_string())
      call print_routine(routine, 'meanVariables%to_string()')

      call print_routine(routine, meanPrognosticPointers%to_string())
        call var_print(meanPrognosticPointers%values())
      call print_routine(routine, 'meanPrognosticPointers%to_string()')
      call print_routine(routine, meanPrognosticPointers%to_string())
    endif
    if (my_process_is_stdio()) call print_routine(routine,'end',stderr=.true.)
#endif

  END SUBROUTINE process_mean_stream

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
    call print_summary("COPY variable:"//TRIM(name),stderr=.true.)
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

    ! copy user defined vertical axes
    dest_element%field%info%cdiZaxisID = source_element%field%info%cdiZaxisID

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
  FUNCTION get_accumulation_varname(varname,output_setup)
    CHARACTER(LEN=VARNAME_LEN)  :: varname
    type(t_output_name_list) :: output_setup

    CHARACTER(LEN=VARNAME_LEN)  :: get_accumulation_varname
    CHARACTER(LEN=1)            :: separator

    separator = '_'
    get_accumulation_varname = &
      &TRIM(varname)//separator//&
      &TRIM(output_setup%operation)//separator//&
      &TRIM(output_setup%output_interval(1))//separator//&
      &TRIM(output_setup%output_start(1))

  END FUNCTION get_accumulation_varname

  !>
  !! implement addition for source fields to the internal accumulation fields
  !!
  SUBROUTINE accumulation_add(source, destination, counter)
    type(t_list_element) , INTENT(IN)    :: source
    type(t_list_element) , INTENT(INOUT) :: destination
    integer, intent(inout)               :: counter

    integer :: index

    index = -1
    if (source%field%info%lcontained) then
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
    else
      ! tackle 1d, 2d and 3d
      SELECT CASE(destination%field%info%ndims)
      CASE(3)
        call add_fields(destination%field%r_ptr(:,:,:,1,1), &
          &             source%field%r_ptr(:,:,:,1,1), &
          &             destination%field%info%subset, &
          &             has_missvals=destination%field%info%lmiss, &
          &             missval=destination%field%info%missval%rval)
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
    endif
    counter                 = counter + 1
  END SUBROUTINE accumulation_add

  !>
  !! Execute the accumulation forall internal variables and compute mean values
  !! if the corresponding event is active
  !!
  SUBROUTINE perform_accumulation(timelevelIndex, timelevelIndex_rcf)
    INTEGER :: timelevelIndex
    INTEGER :: timelevelIndex_rcf

    INTEGER :: element_counter
    class(*),pointer :: varListForMeanEvent, meanEventKey
    class(*),pointer :: sourceVariable, destinationVariable, sourceVariable4Prognostics
    class(*),pointer :: counter, meanEvent, myItem
    class(*),pointer :: eventActive
    type(t_list_element), pointer :: source, destination
    type(vector_iterator) :: meanMapIterator, meanEventIterator
    TYPE(datetime), POINTER :: mtime_date 
    logical :: isactive
    integer :: varcounter
    integer :: timelevel
#ifdef DEBUG_MVSTREAM
    character(len=max_datetime_str_len) :: datetimestring
#endif

    CHARACTER(LEN=*), PARAMETER :: routine =  modname//"::perform_accumulation"

#ifdef DEBUG_MVSTREAM
    if (my_process_is_stdio()) call print_routine(routine,'start',stderr=.true.)
#endif

    meanMapIterator   = meanMap%iter()
    meanEventIterator = meanEvents%iter()

    ! Check events first {{{
    ! this is necessary because of mtime internals
    isactive   = .false.
    mtime_date => newDatetime(time_config%tc_current_date)

#ifdef DEBUG_MVSTREAM
    CALL datetimeToString(mtime_date, datetimestring)
    if (my_process_is_stdio()) call print_summary('Current mtime timestamp:'//trim(datetimestring),stderr=.true.)
#endif

    ! Save results for (not so much) later
    do while (meanEventIterator%next(myItem))
      select type (myItem)
      type is (map_item)
        meanEventKey => myItem%key
        meanEvent    => myItem%value
        select type (meanEvent)
        type is (t_event_wrapper)
          isactive = LOGICAL(isCurrentEventActive(meanEvent%this,mtime_date))
        end select
        call meanEventsActivity%add(meanEventKey,isactive)
      end select
    end do

#ifdef DEBUG_MVSTREAM
    if (my_process_is_stdio()) call print_error(meanEventsActivity%to_string(),stderr=.true.)
#endif
    ! }}}

    do while (meanMapIterator%next(myItem))

      select type (myItem)
      type is (map_item)

        meanEventKey        => myItem%key
        varListForMeanEvent => myItem%value

#ifdef DEBUG_MVSTREAM
        if (my_process_is_stdio()) call print_summary(object_pointer_string(meanEventKey)//"PERFORM ACCU",stderr=.true.) !TODO
#endif

        select type(varListForMeanEvent)
        class is (vector_ref)
          do element_counter=1,varListForMeanEvent%length(),2

#ifdef DEBUG_MVSTREAM
          if (my_process_is_stdio()) call print_routine("perform_accumulation",object_string(element_counter),stderr=.true.)
#endif

            sourceVariable      => varListForMeanEvent%at(element_counter)
            destinationVariable => varListForMeanEvent%at(element_counter+1)

            if (associated(sourceVariable) .and. associated(destinationVariable)) then
              select type (sourceVariable)
              type is (t_list_element)
                select type (destinationVariable)
                type is (t_list_element)
                  ! check for prognostics, pointer must be shifted according to given timelevelIndex {{{
                  if (is_meanPrognosticVariable(destinationVariable)) then
                    ! find the correct pointer:
                    ! check f reduced calling freq. variables are used for
                    ! output? if true, theses variables should be used for
                    ! accumulation, too

#ifdef DEBUG_MVSTREAM
                    if (my_process_is_stdio()) &
                      & call print_error("show meanPrognostic:"//TRIM(destinationVariable%field%info%name),stderr=.true.)
#endif

                    timelevel =  metainfo_get_timelevel(destinationVariable%field%info, 1)
                    source    => get_prognostics_source_pointer (destinationVariable, timelevel)
                    
                  else
                  ! }}}
                    source    => sourceVariable
                  end if
                  destination => destinationVariable
                  counter     => meanVarCounter%get(destination%field%info%name)
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

                    ! FIELD ACCUMULATION {{
                    varcounter = counter !TODO work around for integer pointer, ugly
                    CALL accumulation_add(source, destination, varcounter)
                    counter = varcounter
                    ! }}}
#ifdef DEBUG_MVSTREAM
                    IF ( my_process_is_stdio() )  write (0,*)'new counter: ',counter
#endif

                  end select

                  ! MEAN VALUE COMPUTAION {{{
                  ! check if the field will be written to disk this timestep {{{
                  eventActive => meanEventsActivity%get(meanEventKey)
                  select type (eventActive)
                  type is (logical)
                    isactive = eventActive
                    if ( isactive ) then

#ifdef DEBUG_MVSTREAM
                    if (my_process_is_stdio()) CALL print_summary(" --> PERFORM MEAN VALUE COMP!!!!",stderr=.true.)
#endif

                      counter => meanVarCounter%get(destination%field%info%name)
                      select type(counter)
                      type is (integer)
                        destination%field%r_ptr = destination%field%r_ptr / REAL(counter,wp)
                        counter = 0
                      end select
                    end if
                  end select
                  ! }}}
                end select
                ! ! }}}
              class default
                call finish('perform_accumulation','Found unknown source variable type')
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

  END SUBROUTINE perform_accumulation

  !>
  !! reset internal fields to zero, if the corresponding event is active
  !!
  SUBROUTINE reset_accumulation
    INTEGER :: element_counter
    type(t_list_element), pointer :: destination
    class(*),pointer :: varListForMeanEvent, meanEventKey
    class(*),pointer :: destinationVariable
    class(*),pointer :: meanEvent, myItem
    class(*), pointer :: eventActive
    type(vector_iterator) :: meanMapIterator

    CHARACTER(LEN=*), PARAMETER :: routine =  modname//"::reset_accumulation"

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
#endif

                      destination%field%r_ptr = 0.0_wp
#ifdef DEBUG_MVSTREAM
                    else
                      if (my_process_is_stdio()) call print_error("       eventActive is false",stderr=.true.)
#endif
                    end if
                  class default
                      if (my_process_is_stdio()) call print_error("       eventActive has wrong type",stderr=.true.)
                  end select
              class default
                  if (my_process_is_stdio()) call print_error("     destinationVariable is not t_list_element",stderr=.true.)
              end select
            else
              if (my_process_is_stdio()) call print_error(routine//TRIM(": cannot find destination variable!"),stderr=.true.)
            end if
          end do
        end select 
      end select 
    end do

#ifdef DEBUG_MVSTREAM
    if (my_process_is_stdio()) call print_routine(routine,'finish',stderr=.true.)
#endif

  END SUBROUTINE reset_accumulation

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

END MODULE mo_derived_variable_handling

! vim:tw=0
