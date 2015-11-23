!>
!! Routines for handling proxy variables e.g. accumulation buffers
!!
MODULE mo_derived_variable_handling

  USE self_object
  USE self_vector_ref
  USE self_vector
  USE self_map_ref
  USE self_map
  USE self_assert

  USE mo_kind, ONLY: wp,i8
  USE mo_model_domain, ONLY: t_patch
  USE mo_parallel_config,     ONLY: nproma, use_dp_mpi2io
  USE mo_run_config, ONLY: tc_dt_model
  USE mo_master_config, ONLY: tc_startdate
  USE mo_ocean_nml,           ONLY: n_zlev
  USE mo_var_metadata_types, ONLY: varname_len
  USE mo_impl_constants, ONLY: vname_len, success, max_char_length
  USE mo_name_list_output_types, ONLY: t_output_name_list
  USE mo_mpi, ONLY: my_process_is_stdio
  USE mo_var_list_element, ONLY: level_type_ml, level_type_pl, level_type_hl, level_type_il
  USE mo_name_list_output_config, ONLY: first_output_name_list
  USE mo_var_list, ONLY: nvar_lists, max_var_lists, var_lists, new_var_list,&
       total_number_of_variables, collect_group, get_var_timelevel,&
       get_var_name, default_var_list_settings, add_var, REAL_T, find_element, find_list_element
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

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_derived_variable_handling'

  TYPE(map), SAVE    :: meanMap, meanEvents, meanEventsActivity, meanVarCounter
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

  SUBROUTINE process_mean_stream(p_onl,i_typ, sim_step_info)
    TYPE (t_output_name_list), target :: p_onl
    INTEGER :: i_typ
    TYPE (t_sim_step_info), INTENT(IN) :: sim_step_info

    CHARACTER(LEN=vname_len), POINTER :: in_varlist(:)
    INTEGER :: ntotal_vars, output_variables,i,ierrstat, dataType
    type(vector_ref) :: meanVariables
    CHARACTER(LEN=100) :: eventKey
    type(t_event_wrapper) :: event_wrapper
    class(*), pointer :: myBuffer
    TYPE(t_list_element), POINTER :: src_element, dest_element
    TYPE(vector) :: keys, values 
    CHARACTER(LEN=VARNAME_LEN), ALLOCATABLE :: varlist(:)
    CHARACTER(LEN=VARNAME_LEN) :: dest_element_name
    TYPE(t_accumulation_pair) :: pair
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

      ! create adhoc copies of all variables for later accumulation
      DO i=1, output_variables
        ! collect data variables only, there variables names like
        ! 'grid:clon' which should be excluded
!TODO print *,varlist(i)
        IF ( INDEX(varlist(i),':') > 0) CYCLE
!TODO print *,varlist(i)
 
        ! check for already create meanStream variable (maybe from another output_nml with the same output_interval)
        dest_element_name = get_accumulation_varname(varlist(i),p_onl)
        call print_summary('CHECK NAME:'//TRIM(dest_element_name))
        dest_element => find_list_element(mean_stream_list, trim(dest_element_name))
        IF (.not. ASSOCIATED(dest_element) ) THEN !not found -->> create a new on
          ! find existing source variable
          src_element => find_element ( TRIM(varlist(i)))
          IF (.not. ASSOCIATED (src_element)) THEN
            call finish(routine,'Could not find source variable:'//TRIM(varlist(i)))
          end if
CALL print_summary('src(name)     :|'//trim(src_element%field%info%name)//'|')
CALL print_summary('varlist(name) :|'//trim(in_varlist(i))//'|')
CALL print_summary('new name      :|'//trim(dest_element_name)//'|')
          ! add new variable, copy the meta-data from the existing variable
          ! 1. copy the source variable to destination pointer
          dest_element => copy_var_to_list(mean_stream_list,dest_element_name,src_element)

          ! set output to double precission if necessary
          IF ( use_dp_mpi2io ) THEN                                                                                                                                                                                                                  
            dataType = DATATYPE_FLT64
          ELSE
            dataType = DATATYPE_FLT32
          ENDIF
          dest_element%field%info%cf%datatype = dataType
          if ( my_process_is_stdio()) then
            print *,'copy_var to list CALLED'
          endif
          ! 2. update the nc-shortname to internal name of the source variable
          dest_element%field%info%cf%short_name = src_element%field%info%name
          CALL meanVariables%add(src_element) ! source element comes first
!TODO print *,'src added'
          CALL meanVariables%add(dest_element)
          CALL meanVarCounter%add(dest_element%field%info%name,0)
!         pair = t_accumulation_pair(source=src_element, destination=dest_element, counter=0)
!         CALL meanVariables%add(t_accumulation_pair(source=src_element, destination=dest_element, counter=0))
!TODO print *,'dst added'
        ! replace existince varname in output_nml with the meanStream Variable
CALL print_summary('dst(name)     :|'//trim(dest_element%field%info%name)//'|')
CALL print_summary('dst(shortname):|'//trim(dest_element%field%info%cf%short_name)//'|')
        END IF
        in_varlist(i) = trim(dest_element%field%info%name)
      END DO
!TODO print *,'meanVariables num:',meanVariables%length()
      call meanMap%add(eventKey,meanVariables)
    ELSE
      RETURN
    END IF
    CALL print_error(meanVarCounter%to_string())
  END SUBROUTINE process_mean_stream
  FUNCTION copy_var_to_list(list,name,source_element) RESULT(dest_element)
    TYPE(t_var_list) :: list
    CHARACTER(LEN=VARNAME_LEN) :: name
    TYPE(t_list_element),POINTER :: source_element

    TYPE(t_list_element), POINTER :: dest_element
    integer :: dataType = -1
    CHARACTER(LEN=*), PARAMETER :: routine =  modname//"::copy_var_to_list"

!   CALL message(routine,'START')
    CALL add_var(source_element%field%info%ndims, REAL_T, &
      & list, name, &
      & source_element%field%info%hgrid, source_element%field%info%vgrid, &
      & source_element%field%info%cf, source_element%field%info%grib2,  &
      & source_element%field%info%used_dimensions, &
      & dest_element, &
      & post_op=source_element%field%info%post_op, &
      & action_list=source_element%field%info%action_list, &
      & vert_interp=source_element%field%info%vert_interp, &
      & hor_interp=source_element%field%info%hor_interp, &
      & in_group=source_element%field%info%in_group, &
      & l_pp_scheduler_task=source_element%field%info%l_pp_scheduler_task, &
      & loutput=.TRUE., lrestart=.FALSE., &
      & var_class=source_element%field%info%var_class )
!   CALL message(routine,'FINISH')
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

    destination%field%r_ptr(:,:,:,:,:) = destination%field%r_ptr (:,:,:,:,:)+ source%field%r_ptr(:,:,:,:,:)
    counter                 = counter + 1
  END SUBROUTINE accumulation_add

  SUBROUTINE perform_accumulation
    INTEGER :: k,i
    INTEGER :: element_counter
    class(*),pointer :: elements,check_src, check_pair, check_dest, meanEvent,counter
    type(t_list_element), pointer :: source, destination
    type(vector_iterator) :: value_iterator
    type(vector) :: values, keys
    TYPE(datetime), POINTER :: mtime_date 
    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: mtime_cur_datetime
    character(len=132) :: eventKey
    type(event),pointer :: event_pointer
    type(event) :: event_from_nml
    character(len=132) :: msg
    type(event),pointer :: e
    logical :: isactive
    integer :: varcounter

    class(*), pointer :: eventActive
    TYPE(divisionquotienttimespan) :: quot
!   integer, pointer :: counter
    TYPE(t_accumulation_pair), POINTER :: accumulation_pair
    CHARACTER(LEN=*), PARAMETER :: routine =  modname//"::perform_accumulation"
    
    values = meanMap%values()
    keys   = meanMap%keys()

    ! check events first {{{
    isactive = .false.
    CALL get_datetime_string(mtime_cur_datetime, time_config%cur_datetime)
    mtime_date  => newDatetime(TRIM(mtime_cur_datetime)) 
    if (my_process_is_stdio()) call print_summary('Current mtime timestamp:'//trim(mtime_cur_datetime))

    do i=1,keys%length()
      select type (eventString => keys%at(i))
      type is (character(*))
if (my_process_is_stdio()) call print_summary(eventString)
        meanEvent => meanEvents%get(eventString)
        select type (meanEvent)
        type is (t_event_wrapper)
          isactive = LOGICAL(isCurrentEventActive(meanEvent%this,mtime_date))
        end select
        call meanEventsActivity%add(eventString,isactive)
      end select
    end do
if (my_process_is_stdio()) call print_error(meanEventsActivity%to_string())
    ! }}}

IF ( my_process_is_stdio() ) write(0,*)'values%length = ',values%length() !TODO
    do i=1, values%length()
IF ( my_process_is_stdio() ) write(0,*)'perform_accumulation: i=',i !TODO
      elements => values%at(i)
      select type(elements)
      class is (vector_ref)
IF ( my_process_is_stdio() ) write(0,*)'type: vector' !TODO
        do element_counter=1,elements%length(),2 !start at 2 because the event is at index 1
          check_src => elements%at(element_counter)
          check_dest => elements%at(element_counter+1)
!         check_pair => elements%at(element_counter)
          if (associated(check_src)) then
            select type (check_src)
            type is (t_list_element)
              source      => check_src
            end select
          end if
!         if (associated(check_pair)) then
!           select type (check_pair)
!           type is (t_accumulation_pair)
!             source      => check_pair%source
!             destination => check_pair%destination
!             counter     => check_pair%counter
!           end select
!         end if
          if (associated(check_dest)) then
            select type (check_dest)
            type is (t_list_element)
              destination => check_dest
            end select
          end if
 
!         if (associated(check_pair)) then
!           select type (check_pair)
!           type is (t_accumulation_pair)
          if (associated(check_src)) then
            select type (check_src)
            type is (t_list_element)
              if (associated(check_dest)) then
              select type (check_dest)
              type is (t_list_element)
                counter => meanVarCounter%get(destination%field%info%name)
                select type(counter)
                type is (integer)
IF ( my_process_is_stdio() ) call print_summary('sourceName: '//trim(source%field%info%name))
IF ( my_process_is_stdio() ) call print_summary('destName: '//trim(destination%field%info%name))
IF ( my_process_is_stdio() )  write (0,*)'counter: ',counter
                varcounter = counter
                CALL accumulation_add(source, destination, varcounter)
                counter = varcounter
IF ( my_process_is_stdio() )  write (0,*)'counter: ',counter
                end select

                ! check if the field will be written to disk this timestep {{{
                select type (eventString => keys%at(i))
                type is (character(*))
                  if (my_process_is_stdio()) call print_summary(eventString)
                  ! check if the event is active wrt the current datatime {{{
!                 eventActive => meanEventsActivity%get(eventString)
                  select type (eventActive => meanEventsActivity%get(eventString))
                  type is (logical)
                    isactive = eventActive
                    if ( isactive ) then
if (my_process_is_stdio()) CALL print_summary(" --------------->>>>  PERFORM MEAN VALUE COMP!!!!")

                counter => meanVarCounter%get(destination%field%info%name)
                      select type(counter)
                      type is (integer)
                        destination%field%r_ptr = destination%field%r_ptr / (REAL(counter))
                        counter = 0
                      end select
                    end if
                  end select
                end select
                  ! }}}
              end select
                ! ! }}}
              end if
            class default
              call finish('perform_accumulation','Found unknown source variable type')
            end select
          end if
        end do
      end select 
    end do

  END SUBROUTINE perform_accumulation
  SUBROUTINE reset_accumulation
    INTEGER :: k,i
    INTEGER :: element_counter
    class(*),pointer :: elements,check_src, check_pair, check_dest, meanEvent
    type(t_list_element), pointer :: source, destination
    type(vector_iterator) :: value_iterator
    type(vector) :: values, keys
    TYPE(datetime), POINTER :: mtime_date 
    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: mtime_cur_datetime
    character(len=132) :: eventKey
    type(event),pointer :: event_pointer
    type(event) :: event_from_nml
    character(len=132) :: msg
    type(event),pointer :: e
    integer :: isactive

    class(*), pointer :: eventActive
    TYPE(divisionquotienttimespan) :: quot
!   integer, pointer :: counter
    TYPE(t_accumulation_pair), POINTER :: accumulation_pair
    CHARACTER(LEN=*), PARAMETER :: routine =  modname//"::reset_accumulation"
    
    values = meanMap%values()
    keys   = meanMap%keys()

    do i=1, values%length()
      elements => values%at(i)
      select type(elements)
      class is (vector_ref)
        do element_counter=1,elements%length(),2 !start at 2 because the event is at index 1
          check_dest => elements%at(element_counter+1)
          if (associated(check_dest)) then
            select type (check_dest)
            type is (t_list_element)
              destination => check_dest
              select type (eventString => keys%at(i))
              type is (character(*))
                if (my_process_is_stdio()) call print_summary(eventString)
                ! check if the event is active wrt the current datatime {{{
!                 eventActive => meanEventsActivity%get(eventString)
                select type (eventActive => meanEventsActivity%get(eventString))
                type is (logical)
                  if ( eventActive ) then
if (my_process_is_stdio()) call print_error(eventString//' : ------------ >>>> PERFORM RESET')
                  destination%field%r_ptr = 0.0_wp
                  end if
                end select
              end select
            end select
          end if
        end do
      end select 
    end do

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
