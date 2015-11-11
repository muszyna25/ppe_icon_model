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

  USE mo_kind, ONLY: wp
  USE mo_model_domain, ONLY: t_patch
  USE mo_parallel_config,     ONLY: nproma
  USE mo_ocean_nml,           ONLY: n_zlev
  USE mo_var_metadata_types, ONLY: varname_len
  USE mo_impl_constants, ONLY: vname_len, success, max_char_length
  USE mo_name_list_output_types, ONLY: t_output_name_list
  USE mo_mpi, ONLY: my_process_is_stdio
  USE mo_var_list_element, ONLY: level_type_ml, level_type_pl, level_type_hl, level_type_il
  USE mo_name_list_output_config, ONLY: first_output_name_list
  USE mo_var_list, ONLY: nvar_lists, max_var_lists, var_lists, new_var_list,&
       total_number_of_variables, collect_group, get_var_timelevel,&
       get_var_name, default_var_list_settings, add_var, REAL_T, find_element
  USE mo_linked_list, ONLY: t_var_list, t_list_element
  USE mo_util_string, ONLY: tolower
  USE mo_exception, ONLY: finish, message, message_text
  USE mtime, ONLY: MAX_DATETIME_STR_LEN, newEvent, event, isCurrentEventActive, newDatetime, datetime, eventToString
  USE mo_mtime_extensions,                  ONLY: get_datetime_string
  USE mo_output_event_types, ONLY: t_sim_step_info
  USE mo_time_config,         ONLY: time_config

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_derived_variable_handling'

  TYPE(map)       , SAVE :: meanMap, meanEvents
  TYPE(t_var_list)   :: mean_stream_list
  INTEGER, PARAMETER :: ntotal = 10

  PUBLIC :: init_mean_stream
  PUBLIC :: finish_mean_stream
  PUBLIC :: collect_meanstream_variables
  PUBLIC :: mean_stream_list
  PUBLIC :: copy_var_to_list
  PUBLIC :: perform_accumulation
  PUBLIC :: process_mean_stream

  TYPE :: t_accumulation_pair
    TYPE(t_list_element), POINTER :: source, destination
  END TYPE t_accumulation_pair

!!!  SUBROUTINE collect_target_variables()
!!!  END SUBROUTINE collect_target_variables
  
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
    
    integer :: i
    
    meanMap = map()

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

  SUBROUTINE process_mean_stream(p_onl,i_typ)
    TYPE (t_output_name_list), target :: p_onl
    INTEGER :: i_typ

    CHARACTER(LEN=vname_len), POINTER :: in_varlist(:)
    INTEGER :: ntotal_vars, output_variables,i,ierrstat
    type(vector_ref) :: meanVariables
    CHARACTER(LEN=100) :: eventKey
    class(*), pointer :: myBuffer
    TYPE(t_list_element), POINTER :: src_element, dest_element
    TYPE(vector) :: keys, values 
    CHARACTER(LEN=VARNAME_LEN), ALLOCATABLE :: varlist(:)
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

      ! collect mean variables {{{
      eventKey = get_event_key(p_onl)
      write (0,*)'eventKey:',trim(eventKey)
      IF ( meanMap%has_key(eventKey) ) THEN
        myBuffer => meanMap%get(eventKey)
        select type (myBuffer)
        type is (vector_ref)
          meanVariables = myBuffer
        end select
      ELSE
        meanVariables = vector_ref()
  !     e => newEvent( &
  !       &                 sim_step_info%sim_start, &
  !       &                 p_onl%output_start(1), &
  !       &                 p_onl%output_end(1), &
  !       &                 p_onl%output_interval(1) &
  !       &                )

  !     call meanEvents%add(eventKey,newEvent(eventKey, &
  !       &                 sim_step_info%sim_start, &
  !       &                 p_onl%output_start(1), &
  !       &                 p_onl%output_end(1), &
  !       &                 p_onl%output_interval(1) &
  !       &                ))
      END IF

      ! create adhoc copies of all variables for later accumulation
      DO i=1, output_variables
        ! collect data variables only, there variables names like
        ! 'grid:clon' which should be excluded
!TODO print *,varlist(i)
        IF ( INDEX(varlist(i),':') > 0) CYCLE
!TODO print *,varlist(i)
 
        ! find existing variable
        src_element => find_element ( TRIM(varlist(i)))
        IF (.not. ASSOCIATED (src_element)) THEN
          call finish(routine,'Could not find source variable:'//TRIM(varlist(i)))
        end if
CALL print_summary('src(name)     :|'//trim(src_element%field%info%name)//'|')
CALL print_summary('varlist(name) :|'//trim(in_varlist(i))//'|')
CALL print_summary('new name      :|'//trim(get_accumulation_varname(varlist(i),p_onl))//'|')
        ! add new variable, copy the meta-data from the existing variable
        ! 1. copy the source variable to destination pointer
        dest_element => copy_var_to_list(mean_stream_list,get_accumulation_varname(varlist(i),p_onl),src_element)
        if ( my_process_is_stdio()) then
          print *,'copy_var to list CALLED'
        endif
        ! 2. update the nc-shortname to internal name of the source variable
        dest_element%field%info%cf%short_name = src_element%field%info%name
        CALL meanVariables%add(src_element) ! source element comes first
!TODO print *,'src added'
        CALL meanVariables%add(dest_element)
!TODO print *,'dst added'
        ! replace existince varname in output_nml with the meanStream Variable
        in_varlist(i) = trim(dest_element%field%info%name)
CALL print_summary('dst(name)     :|'//trim(dest_element%field%info%name)//'|')
CALL print_summary('dst(shortname):|'//trim(dest_element%field%info%cf%short_name)//'|')
      END DO
!TODO print *,'meanVariables num:',meanVariables%length()
          call meanMap%add(eventKey,meanVariables)
    ELSE
      RETURN
    END IF
  END SUBROUTINE process_mean_stream
  FUNCTION copy_var_to_list(list,name,source_element) RESULT(dest_element)
    TYPE(t_var_list) :: list
    CHARACTER(LEN=VARNAME_LEN) :: name
    TYPE(t_list_element),POINTER :: source_element

    TYPE(t_list_element), POINTER :: dest_element
    CHARACTER(LEN=*), PARAMETER :: routine =  modname//"::copy_var_to_list"

!   CALL message(routine,'START')
    CALL add_var(source_element%field%info%ndims, REAL_T, &
      & list, name, &
      & source_element%field%info%hgrid, source_element%field%info%vgrid, &
      & source_element%field%info%cf, source_element%field%info%grib2,  &
      & source_element%field%info%used_dimensions, &
      & dest_element, &
      & post_op=source_element%field%info%post_op, &
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

  SUBROUTINE accumulation_add(source, destination)
    type(t_list_element) , INTENT(IN)    :: source
    type(t_list_element) , INTENT(INOUT) :: destination

    destination%field%r_ptr = destination%field%r_ptr + source%field%r_ptr
  END SUBROUTINE accumulation_add

  SUBROUTINE perform_accumulation
    INTEGER :: k,i
    INTEGER :: element_counter
    class(*),pointer :: elements,check_src, check_dest, meanEvent
    type(t_list_element), pointer :: source, destination
    type(vector_iterator) :: value_iterator
    type(vector) :: values, keys
    TYPE(datetime), POINTER :: mtime_date 
    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: mtime_cur_datetime
    character(len=132) :: eventKey
    class(*), pointer :: eventString
    type(event),pointer :: event_pointer
    type(event) :: event_from_nml
    character(len=132) :: msg
    type(event),pointer :: e
    logical :: isactive
    
    values = meanMap%values()
    keys   = meanMap%keys()

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
          if (associated(check_src)) then
            select type (check_src)
            type is (t_list_element)
              source      => check_src
            end select
          end if
          if (associated(check_dest)) then
            select type (check_dest)
            type is (t_list_element)
              destination => check_dest
            end select
          end if
 
          if (associated(check_src)) then
            select type (check_src)
            type is (t_list_element)
              !if (associated(check_dest)) then
              select type (check_dest)
              type is (t_list_element)
                IF ( my_process_is_stdio() ) write(0,*)'sourceName:',trim(source%field%info%name)
                IF ( my_process_is_stdio() ) write(0,*)'destName:',trim(destination%field%info%name)
                CALL accumulation_add(source, destination)
              class default
                call finish('perform_accumulation','Found unknown destination variable type')
              end select
              !end if
            class default
              call finish('perform_accumulation','Found unknown source variable type')
            end select
          end if
        end do
        eventString => keys%at(i)
       !do k=1,3
       !  e => meanEvents(k)
       !  call eventToString(e, msg)
       !  IF ( my_process_is_stdio() ) THEN
       !    write (0,*)' k  :', k  
       !    call print_summary(msg,stderr=.true.)
       !    if (msg == eventKey) then
       !      ! found the correnspondin even
       !      call print_summary("key found !!!",stderr=.true.)
       !      ! now check for activity
       !      CALL get_datetime_string(mtime_cur_datetime, time_config%cur_datetime)
       !      call print_summary(trim(mtime_cur_datetime))
       !      mtime_date  => newDatetime(TRIM(mtime_cur_datetime)) 
       !      isactive = LOGICAL(isCurrentEventActive(e,mtime_date))
       !      write (0,*)'---- isactive ----- ',isactive
       !    end if
       !  end if
       !end do
      end select 
    end do

  END SUBROUTINE perform_accumulation
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
