module derived_variable_handling
  use self_vector
  use self_map
  use self_assert

  USE mo_var_metadata_types,                ONLY: VARNAME_LEN
  USE mo_impl_constants,                    ONLY: vname_len
  USE mo_name_list_output_types,            ONLY: t_output_name_list

  implicit none
  private

  public :: collect_target_variables
  public :: collect_meanStream_variables

contains

  subroutine collect_target_variables
  end subroutine collect_target_variables

  SUBROUTINE collect_meanStream_variables()
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::collect_meanStream_variables"
    !
    CHARACTER(LEN=VARNAME_LEN), ALLOCATABLE :: varlist(:)
    CHARACTER(LEN=VARNAME_LEN)              :: varname, mean_varname, message_text
    INTEGER                                 :: nvars, i_typ, ierrstat, i, ntotal_vars
    INTEGER                                 :: varlist_length
    CHARACTER(LEN=vname_len),  POINTER      :: in_varlist(:)
    TYPE (t_output_name_list), POINTER      :: p_onl

    type(map) :: meanMap
    type(vector) :: meanVariables

    meanMap = map()
    meanVariables = vector()


    ntotal_vars = total_number_of_variables()
    ! temporary variables needed for variable group parsing
    ALLOCATE(varlist(ntotal_vars), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')

    nvars = 1
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
   
          IF (nvars > 0)  varlist(1:nvars) = in_varlist(1:nvars)
          varlist((nvars+1):ntotal_vars) = " "

          IF ( meanMap%has_key(trim(p_onl%output_interval(1))) ) THEN
            !call meanMap%get(trim(p_onl%output_interval(1)),meanVariables)
            meanVariables = vector()
          ELSE
            meanVariables = vector()
          END IF
          IF (i_typ == level_type_ml) THEN
            do i=1,nvars
            ! collect data variables only
            if ( index(varlist(i),':') < 1 ) then
              call print_green('var:'//trim(varlist(i))//'---')
              call meanVariables%addUniq(trim(varlist(i)))
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

        END DO
      END IF
      p_onl => p_onl%next
    END DO
    IF ( my_process_is_stdio() ) THEN
      call print_aqua('collected variables:{{{{{{{{{{{{{{{{{{{{{{{{{{{')
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
end module derived_variable_handling
