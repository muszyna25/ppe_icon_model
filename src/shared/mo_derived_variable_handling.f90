module mo_derived_variable_handling
  use self_vector
  use self_map
  use self_assert

  USE mo_var_metadata_types,                ONLY: VARNAME_LEN
  USE mo_impl_constants,                    ONLY: vname_len, SUCCESS
  USE mo_name_list_output_types,            ONLY: t_output_name_list
  USE mo_mpi,                               ONLY: my_process_is_stdio
  USE mo_var_list_element,                  ONLY: level_type_ml, level_type_pl, level_type_hl, level_type_il
  USE mo_name_list_output_config,           ONLY: first_output_name_list
  USE mo_var_list,                          ONLY: nvar_lists, max_var_lists, var_lists,           &
    &                                             new_var_list,                                   &
    &                                             total_number_of_variables, collect_group,       &
    &                                             get_var_timelevel, get_var_name
  USE mo_linked_list,                       ONLY: find_list_element
  USE mo_exception,                         ONLY: finish, message, message_text

  implicit none

  private

  type(map)   , save :: meanMap
  type(vector), save :: meanVariables
  integer, parameter :: ntotal = 1024
  CHARACTER(LEN=VARNAME_LEN) , SAVE, target :: varlist_buffer(ntotal)
  CHARACTER(LEN=VARNAME_LEN) , SAVE, target :: periods_buffer(ntotal)

  public :: init_mean_stream
  public :: collect_meanstream_variables

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_derived_variable_handling'

! subroutine collect_target_variables()
! end subroutine collect_target_variables
CONTAINS

  SUBROUTINE init_mean_stream()
    meanMap = map()
    meanVariables = vector()

    varlist_buffer = ''
    periods_buffer = ''
  end subroutine init_mean_stream
  SUBROUTINE collect_meanstream_variables()

    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::collect_meanStream_variables"
    !
    CHARACTER(LEN=VARNAME_LEN)              :: varname, mean_varname, message_text
    INTEGER                                 :: nvars, i_typ, ierrstat, i, ntotal_vars,j
    INTEGER                                 :: varlist_length, periods_counter, v_counter
    CHARACTER(LEN=vname_len),  POINTER      :: in_varlist(:)
    TYPE (t_output_name_list), POINTER      :: p_onl
    CHARACTER(LEN=VARNAME_LEN), ALLOCATABLE :: varlist(:)

    nvars           = 1
    periods_counter = 1
    v_counter       = 1
    ntotal_vars     = total_number_of_variables()
    ! temporary variables needed for variable group parsing
    ALLOCATE(varlist(ntotal_vars), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')

    ! -- loop over all output namelists
    p_onl => first_output_name_list

    DO
      IF (.NOT.ASSOCIATED(p_onl)) EXIT
      IF ("mean" .EQ. TRIM(p_onl%operation)) THEN
        WRITE(message_text,'(3a)') 'outputInterval: ',TRIM(p_onl%output_interval(1))
        CALL message('',message_text)

        DO i_typ = 1, 4
   
          IF (i_typ == level_type_ml) in_varlist => p_onl%ml_varlist
          IF (i_typ == level_type_pl) in_varlist => p_onl%pl_varlist
          IF (i_typ == level_type_hl) in_varlist => p_onl%hl_varlist
          IF (i_typ == level_type_il) in_varlist => p_onl%il_varlist
   
          varlist_length = SIZE(in_varlist)

          DO
            IF (in_varlist(nvars) == ' ') EXIT
            nvars = nvars + 1
          END DO
          nvars = nvars - 1

WRITE(message_text,FMT=*) 'nvars: ',nvars
CALL message('',message_text)
   
          IF (nvars > 0)  varlist(1:nvars) = in_varlist(1:nvars)
          varlist((nvars+1):ntotal_vars) = " "

          periods_buffer(periods_counter) = trim(p_onl%output_interval(1))

          IF ( meanMap%has_key(periods_buffer(periods_counter)) ) THEN
            call meanMap%get(trim(p_onl%output_interval(1)),meanVariables)
            !meanVariables = vector()
          ELSE
            meanVariables = vector()
          END IF
          IF (i_typ == level_type_ml) THEN
            do i=1,nvars
              ! collect data variables only
              if ( index(varlist(i),':') < 1 ) then
                j = (periods_counter-1)*nvars + i
                varlist_buffer(j) = trim(varlist(i))
                call print_green('var:'//trim(varlist_buffer(j))//'---')
                call meanVariables%addUniq(trim(varlist_buffer(j)))
              end if
            end do
            ! add stuff to map
IF ( my_process_is_stdio() ) THEN
  call print_aqua('collected variables:{{{{{{{{{{{{{{{{{{{{{{{{{{{')
  call meanVariables%print()
  call print_aqua('}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}')
END IF
            call meanMap%add(trim(p_onl%output_interval(1)),meanVariables)
          END IF
periods_counter = periods_counter + 1
        END DO
      END IF
      p_onl => p_onl%next
    END DO
    IF ( my_process_is_stdio() ) THEN
      call print_aqua('collected map {{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{')
      call meanMap%print()
      call print_aqua('}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}')
    END IF
    !
    !1. Collect variables uniq by the output interval
    !   this will allow collective events for all variables in this group
    !2. for each variable, store
    !      source pointer [got with find_element]
    !      target pointer: copy meta info the the source, but choose new name with '_m'
    !
    !
    !write(0,*)'varlist:',varlist
  END SUBROUTINE collect_meanStream_variables
end module mo_derived_variable_handling
