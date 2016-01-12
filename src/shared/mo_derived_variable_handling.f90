!>
!! Routines for handling proxy variables e.g. accumulation buffers
!!
MODULE mo_derived_variable_handling

  USE self_object
  USE self_vector_ref
  USE self_vector
  USE self_map
  USE self_assert

  USE mo_kind, ONLY: wp,i8
  USE mo_model_domain, ONLY: t_patch
  USE mo_parallel_config,     ONLY: nproma, use_dp_mpi2io
  USE mo_io_config, ONLY: lnetcdf_flt64_output
  USE mo_run_config, ONLY: tc_dt_model
  USE mo_dynamics_config, ONLY: nnow, nnew, nold
  USE mo_master_config, ONLY: tc_startdate
  USE mo_ocean_nml,           ONLY: n_zlev
  USE mo_statistics, ONLY: add_fields
  USE mo_var_metadata_types, ONLY: varname_len
  USE mo_impl_constants, ONLY: vname_len, success, max_char_length, TLEV_NNOW, TLEV_NNEW, TLEV_NNOW_RCF, TLEV_NNEW_RCF
  USE mo_name_list_output_types, ONLY: t_output_name_list
  USE mo_mpi, ONLY: my_process_is_stdio
  USE mo_var_list_element, ONLY: level_type_ml, level_type_pl, level_type_hl, level_type_il
  USE mo_name_list_output_config, ONLY: first_output_name_list
  USE mo_name_list_output_metadata, ONLY: metainfo_get_timelevel
  USE mo_var_list, ONLY: nvar_lists, max_var_lists, var_lists, new_var_list,&
       total_number_of_variables, collect_group, get_var_timelevel,&
       get_var_name, default_var_list_settings, add_var, REAL_T, find_element, find_list_element, &
       & get_varname_with_timelevel
  USE mo_linked_list, ONLY: t_var_list, t_list_element
  USE mo_util_string, ONLY: tolower
  USE mo_exception, ONLY: finish, message, message_text
  USE mtime, ONLY: MAX_DATETIME_STR_LEN, newEvent, event, isCurrentEventActive,&
    & newDatetime, datetime, eventToString, divideDatetimeDifferenceInSeconds, &
    & divisionquotientTimespan
  USE mo_mtime_extensions,                  ONLY: get_datetime_string
  USE mo_output_event_types, ONLY: t_sim_step_info
  USE mo_time_config,         ONLY: time_config
  USE mo_cdi,                  ONLY: DATATYPE_FLT32, DATATYPE_FLT64
  USE mo_cdi_constants, ONLY: GRID_UNSTRUCTURED_CELL, GRID_UNSTRUCTURED_EDGE, GRID_UNSTRUCTURED_VERT

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_derived_variable_handling'

  TYPE(map), SAVE    :: meanMap, meanEvents, meanEventsActivity, meanVarCounter, meanPrognosticPointers
  TYPE(vector), SAVE :: meanPrognostics
  TYPE(t_var_list)   :: mean_stream_list
  INTEGER, PARAMETER :: ntotal = 10

  PUBLIC :: init_mean_stream
  PUBLIC :: finish_mean_stream
  PUBLIC :: mean_stream_list
  PUBLIC :: copy_var_to_list
  PUBLIC :: perform_accumulation
  PUBLIC :: reset_accumulation
  PUBLIC :: process_mean_stream

  TYPE :: t_accumulation_pair
    TYPE(t_list_element), POINTER :: source
    TYPE(t_list_element), POINTER :: destination
    INTEGER                       :: counter = 0
  END TYPE t_accumulation_pair

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

  logical function is_meanPrognosticVariable(variable)
    type(t_list_element), intent(in) :: variable
    is_meanPrognosticVariable = meanPrognostics%includes(variable%field%info%cf%short_name)
  end function is_meanPrognosticVariable

  function get_prognstics_source_pointer(destinationVariable, timelevelIndex) result(sourcePointer)
    type(t_list_element), pointer :: sourcePointer
    type(t_list_element)          :: destinationVariable
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
    end select
  end function get_prognstics_source_pointer


  !>
  !! Create a variable list
  !!
  SUBROUTINE init_mean_stream(patch_2d)
    TYPE(t_patch), TARGET, INTENT(in) :: patch_2d
    
    CHARACTER(LEN=max_char_length) :: listname
    
    meanMap            = map(verbose=.false.)
    meanEvents         = map(verbose=.false.)
    meanEventsActivity = map(verbose=.false.)
    meanVarCounter     = map(verbose=.false.)
    meanPrognostics    = vector(verbose=.false.)
    meanPrognosticPointers    = map(verbose=.false.)

    WRITE(listname,'(a)')  'mean_stream_list'
    CALL new_var_list(mean_stream_list, listname, patch_id=patch_2d%id)
    CALL default_var_list_settings( mean_stream_list,lrestart=.FALSE.,loutput=.TRUE., model_type='oce' )
  END SUBROUTINE init_mean_stream

  !>
  !!
  !!
  SUBROUTINE finish_mean_stream()
    IF (my_process_is_stdio()) THEN
    CALL print_summary(&
         '==== FINISH MAP ===================================================')
    print *,meanMap%to_string()
    CALL print_summary(&
         '===================================================================')
    END IF
  END SUBROUTINE finish_mean_stream

  integer function output_varlist_length(in_varlist)
    CHARACTER(LEN=vname_len), intent(in) :: in_varlist(:)

    output_varlist_length = 0
    DO
      IF (in_varlist(output_varlist_length+1) == ' ') EXIT
      output_varlist_length = output_varlist_length + 1
    END DO
  end function output_varlist_length

  SUBROUTINE process_mean_stream(p_onl,i_typ, sim_step_info, patch_2d)
    TYPE (t_output_name_list), target :: p_onl
    INTEGER :: i_typ
    TYPE (t_sim_step_info), INTENT(IN) :: sim_step_info
    TYPE(t_patch), INTENT(IN) :: patch_2d

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

    IF ("mean" .EQ. TRIM(p_onl%operation)) THEN

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

      eventKey = get_event_key(p_onl)
      if (my_process_is_stdio()) write (0,*)'eventKey:',trim(eventKey)
      IF ( meanMap%has_key(eventKey) ) THEN
        myBuffer => meanMap%get(eventKey)
        select type (myBuffer)
        type is (vector_ref)
          meanVariables = myBuffer
        end select
      ELSE
        meanVariables = vector_ref()

        event_wrapper = t_event_wrapper(this=newEvent(eventKey, &
          &                 sim_step_info%sim_start, &
          &                 p_onl%output_start(1), &
          &                 p_onl%output_end(1), &
          &                 p_onl%output_interval(1) &
          &                ))


        call meanEvents%add(eventKey,event_wrapper)
      END IF
      call meanEventsActivity%add(eventKey,.false.)

      ! create adhoc copies of all variables for later accumulation
      DO i=1, output_variables
        ! collect data variables only, there variables names like
        ! 'grid:clon' which should be excluded
      print *,varlist(i)
        IF ( INDEX(varlist(i),':') > 0) CYCLE
      print *,varlist(i)
 
        ! check for already create meanStream variable (maybe from another output_nml with the same output_interval)
        dest_element_name = get_accumulation_varname(varlist(i),p_onl)
              call print_summary('CHECK NAME:'//TRIM(dest_element_name))
        dest_element => find_list_element(mean_stream_list, trim(dest_element_name))
        IF (.not. ASSOCIATED(dest_element) ) THEN !not found -->> create a new on
          ! find existing source variable
          src_element => find_element ( TRIM(varlist(i)))
          IF (.not. ASSOCIATED (src_element)) THEN
            ! try to find timelevel variables 
            foundPrognostic = .false.
            timelevels = (/nold(1),nnow(1),nnew(1)/)
            do timelevel=1,3
              call print_error(get_varname_with_timelevel(varlist(i),timelevels(timelevel)))
              src_element => find_element(get_varname_with_timelevel(varlist(i),timelevels(timelevel)))
              if ( ASSOCIATED(src_element) ) then
                if ( .not. foundPrognostic ) then
                  call meanPrognostics%add(varlist(i))
                  foundPrognostic = .true.
                end if
                call meanPrognosticPointers%add(get_varname_with_timelevel(varlist(i),timelevels(timelevel)),src_element)
              end if
            end do
          END IF
          IF (.not. ASSOCIATED (src_element)) THEN
            call finish(routine,'Could not find source variable:'//TRIM(varlist(i)))
          END IF
      if ( my_process_is_stdio())CALL print_summary('src(name)     :|'//trim(src_element%field%info%name)//'|')
      if ( my_process_is_stdio())CALL print_summary('varlist(name) :|'//trim(in_varlist(i))//'|')
      if ( my_process_is_stdio())CALL print_summary('new name      :|'//trim(dest_element_name)//'|')
          ! add new variable, copy the meta-data from the existing variable
          ! 1. copy the source variable to destination pointer
          dest_element => copy_var_to_list(mean_stream_list,dest_element_name,src_element, patch_2d)

          ! set output to double precission if necessary
          dest_element%field%info%cf%datatype = MERGE(DATATYPE_FLT64, DATATYPE_FLT32, lnetcdf_flt64_output)
      if ( my_process_is_stdio()) print *,'copy_var to list CALLED'
          ! 2. update the nc-shortname to internal name of the source variable
          dest_element%field%info%cf%short_name = get_var_name(src_element%field)

          ! Collect variable pointers {{{
          CALL meanVariables%add(src_element)
      print *,'src added'
          CALL meanVariables%add(dest_element)
          CALL meanVarCounter%add(dest_element%field%info%name,0)
      print *,'dst added'
        ! replace existince varname in output_nml with the meanStream Variable
      if ( my_process_is_stdio())CALL print_summary('dst(name)     :|'//trim(dest_element%field%info%name)//'|')
      if ( my_process_is_stdio())CALL print_summary('dst(shortname):|'//trim(dest_element%field%info%cf%short_name)//'|')
        END IF
        in_varlist(i) = trim(dest_element%field%info%name)
      END DO
      print *,'meanVariables num:',meanVariables%length()
      call meanMap%add(eventKey,meanVariables)
    ELSE
      RETURN
    END IF
          if (my_process_is_stdio())  CALL print_error(meanVarCounter%to_string())
  END SUBROUTINE process_mean_stream

  FUNCTION copy_var_to_list(list,name,source_element,patch_2d) RESULT(dest_element)
    TYPE(t_var_list) :: list
    CHARACTER(LEN=VARNAME_LEN) :: name
    TYPE(t_list_element),POINTER :: source_element
    TYPE(t_patch),TARGET, INTENT(in)    :: patch_2d

    TYPE(t_list_element), POINTER :: dest_element
    integer :: dataType = -1
    CHARACTER(LEN=*), PARAMETER :: routine =  modname//"::copy_var_to_list"

    dataType = MERGE(DATATYPE_FLT64, DATATYPE_FLT32, lnetcdf_flt64_output)

    call print_summary("COPY variable:"//TRIM(name))
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
      & initval_r=source_element%field%info%initval%rval, &
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

  FUNCTION get_accumulation_varname(varname,output_setup)
    CHARACTER(LEN=VARNAME_LEN)  :: varname
    type(t_output_name_list) :: output_setup

    CHARACTER(LEN=VARNAME_LEN)  :: get_accumulation_varname
    CHARACTER(LEN=1)            :: separator

    separator = '_'
    get_accumulation_varname = &
      &TRIM(varname)//separator//&
      &trim(output_setup%operation)//separator//&
      &TRIM(output_setup%output_interval(1))//separator//&
      &trim(output_setup%output_start(1))

  END FUNCTION get_accumulation_varname

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
      ! tackle 2d and 3d
      if ( 3 .eq. destination%field%info%ndims ) then
        call add_fields(destination%field%r_ptr(:,:,:,1,1), &
          &             source%field%r_ptr(:,:,:,1,1), &
          &             destination%field%info%subset, &
          &             has_missvals=destination%field%info%lmiss, &
          &             missval=destination%field%info%missval%rval)
      else
        call add_fields(destination%field%r_ptr(:,:,1,1,1), &
          &             source%field%r_ptr(:,:,1,1,1), &
          &             destination%field%info%subset, &
          &             has_missvals=destination%field%info%lmiss, &
          &             missval=destination%field%info%missval%rval)
      end if
!     destination%field%r_ptr(:,:,:,:,:) = destination%field%r_ptr (:,:,:,:,:) + source%field%r_ptr(:,:,:,:,:)
    endif
    counter                 = counter + 1
  END SUBROUTINE accumulation_add

  SUBROUTINE perform_accumulation(timelevelIndex, timelevelIndex_rcf)
    INTEGER :: timelevelIndex
    INTEGER :: timelevelIndex_rcf

    INTEGER :: element_counter
    class(*),pointer :: varListForMeanEvent, meanEventKey
    class(*),pointer :: sourceVariable, destinationVariable, sourceVariable4Prognostics
    class(*),pointer :: counter, meanEvent, myItem
    class(*), pointer :: eventActive
    type(t_list_element), pointer :: source, destination
    type(vector_iterator) :: meanMapIterator, meanEventIterator
    TYPE(datetime), POINTER :: mtime_date 
    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: mtime_cur_datetime
    logical :: isactive
    integer :: varcounter
    real(wp) :: counter_inverse
    integer :: timelevel

    TYPE(divisionquotienttimespan) :: quot
!   integer, pointer :: counter
    CHARACTER(LEN=*), PARAMETER :: routine =  modname//"::perform_accumulation"
    
    call print_routine(routine,'start')

    meanMapIterator   = meanMap%iter()
    meanEventIterator = meanEvents%iter()

    ! Check events first {{{
    isactive = .false.
    CALL get_datetime_string(mtime_cur_datetime, time_config%cur_datetime)
    mtime_date  => newDatetime(TRIM(mtime_cur_datetime)) 
          if (my_process_is_stdio()) call print_summary('Current mtime timestamp:'//trim(mtime_cur_datetime))
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
      if (my_process_is_stdio()) call print_error(meanEventsActivity%to_string()) !TODO
    ! }}}
!#ifdef T
    do while (meanMapIterator%next(myItem))

      select type (myItem)
      type is (map_item)

        meanEventKey        => myItem%key
        varListForMeanEvent => myItem%value

      if (my_process_is_stdio()) call print_summary(object_pointer_string(meanEventKey)//"PERFORM ACCU") !TODO
        select type(varListForMeanEvent)
        class is (vector_ref)
      !IF ( my_process_is_stdio() ) write(0,*)'type: vector' !TODO
          do element_counter=1,varListForMeanEvent%length(),2
!           call print_routine("perform_accumulation",object_string(element_counter))

            sourceVariable      => varListForMeanEvent%at(element_counter)
            destinationVariable => varListForMeanEvent%at(element_counter+1)
!           call print_routine("perform_accumulation","after calling 'AT'")

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

       call print_error("show meanPrognostic:"//TRIM(destinationVariable%field%info%name))
                    SELECT CASE (destinationVariable%field%info%tlev_source)
                    CASE(TLEV_NNOW);     timelevel = timelevelIndex
                    CASE(TLEV_NNEW);     timelevel = timelevelIndex
                    CASE(TLEV_NNOW_RCF); timelevel = timelevelIndex_rcf
                    CASE(TLEV_NNEW_RCF); timelevel = timelevelIndex_rcf
                    CASE DEFAULT;        timelevel = timelevelIndex
                    END SELECT
                    timelevel = metainfo_get_timelevel(destinationVariable%field%info, 1)
                    source         => get_prognstics_source_pointer (destinationVariable, timelevel)
                    
                  else
                  ! }}}
                  source      => sourceVariable
                  end if
                  destination => destinationVariable
                  counter     => meanVarCounter%get(destination%field%info%name)
                  select type (counter)
                  type is (integer)
      IF ( my_process_is_stdio() ) call print_summary('sourceName : '//trim(source%field%info%name))
      IF ( my_process_is_stdio() ) call print_summary('destName   : '//trim(destination%field%info%name))
      IF ( my_process_is_stdio() ) call print_summary('destNameOut: '//trim(destination%field%info%cf%short_name))
      IF ( my_process_is_stdio() )  write (0,*)'counter: ',counter
                    ! ACCUMULATION {{
                    varcounter = counter !TODO work around for integer pointer, ugly
                    CALL accumulation_add(source, destination, varcounter)
                    counter = varcounter
                    ! }}}
      IF ( my_process_is_stdio() )  write (0,*)'counter: ',counter
                  end select

                  ! MEAN VALUE COMPUTAION {{{
                  ! check if the field will be written to disk this timestep {{{
                  eventActive => meanEventsActivity%get(meanEventKey)
                  select type (eventActive)
                  type is (logical)
                    isactive = eventActive
                    if ( isactive ) then
      if (my_process_is_stdio()) CALL print_summary(" --------------->>>>  PERFORM MEAN VALUE COMP!!!!")

                      counter => meanVarCounter%get(destination%field%info%name)
                      select type(counter)
                      type is (integer)
                        counter_inverse = 1/REAL(counter,wp)
                        destination%field%r_ptr = destination%field%r_ptr * counter_inverse
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
!#endif
    call print_routine(routine,'finish')
  END SUBROUTINE perform_accumulation

  SUBROUTINE reset_accumulation
    INTEGER :: element_counter
    type(t_list_element), pointer :: destination
    class(*),pointer :: varListForMeanEvent, meanEventKey
    class(*),pointer :: destinationVariable
    class(*),pointer :: meanEvent, myItem
    class(*), pointer :: eventActive
    type(vector_iterator) :: meanMapIterator

    CHARACTER(LEN=*), PARAMETER :: routine =  modname//"::reset_accumulation"

         call print_routine(routine,'start')

    meanMapIterator   = meanMap%iter()

      if (my_process_is_stdio()) call print_error("got iterator")
    do while (meanMapIterator%next(myItem))

      if (my_process_is_stdio()) call print_error("started event loop")
      select type (myItem)
      type is (map_item)

        meanEventKey        => myItem%key
        varListForMeanEvent => myItem%value

        select type(varListForMeanEvent)
        class is (vector_ref)

          do element_counter=1,varListForMeanEvent%length(),2 !start at 2 because the event is at index 1

            destinationVariable => varListForMeanEvent%at(element_counter+1)
       if (my_process_is_stdio()) call print_error("got destinationVariable")
            if (associated(destinationVariable)) then
       if (my_process_is_stdio()) call print_error("    destinationVariable is associated")
              select type (destinationVariable)
              type is (t_list_element)
                  destination => destinationVariable
       if (my_process_is_stdio()) call print_error("    destinationVariable is t_list_element")
                  eventActive => meanEventsActivity%get(meanEventKey)
       if (my_process_is_stdio()) call print_error("       eventActive got       ")
                  if (.not.associated(eventActive)) then
                    if (my_process_is_stdio()) call print_error("       eventActive not associated")
                  end if
                  select type (eventActive)
                  type is (logical)
       if (my_process_is_stdio()) call print_error("       eventActive is logical")
                    if ( LOGICAL(eventActive) ) then
       if (my_process_is_stdio()) call print_error("       eventActive is true")
       if (my_process_is_stdio()) call print_error(object_string(meanEventKey)//' : ------------ >>>> PERFORM RESET')
                      destination%field%r_ptr = 0.0_wp
                    else
      if (my_process_is_stdio()) call print_error("       eventActive is false")
                    end if
                  class default
      if (my_process_is_stdio()) call print_error("       eventActive has wrong type")
                  end select
              class default
      if (my_process_is_stdio()) call print_error("     destinationVariable is not t_list_element")
              end select
            else
              call print_error(routine//TRIM(": cannot find destination variable!"))
            end if
          end do
        end select 
      end select 
    end do

         call print_routine(routine,'finish')
  END SUBROUTINE reset_accumulation
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
