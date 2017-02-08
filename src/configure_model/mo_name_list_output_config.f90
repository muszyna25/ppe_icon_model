!>
!! @brief configuration setup for output name lists
!!
!! Configuration setup for output name lists
!!
!! @author R.Johanni, F. Prill
!!
!! @note This is only a preliminary implementation of a configure state
!!       for output name lists; based on R. Johanni's name list
!!       handling which was originally implemented in
!!       "shared/mo_name_list_output"
!!
!! @par Revision History
!! Moved configure state from shared/mo_name_list_output:
!! F. Prill, DWD (2012-01-26)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_name_list_output_config

  USE mo_kind,                  ONLY: wp
  USE mo_exception,             ONLY: finish
  USE mo_impl_constants,        ONLY: max_var_ml, max_var_pl, &
    &                                 max_var_hl, max_var_il, &
    &                                 MAX_TIME_LEVELS
  USE mo_cdi,                   ONLY: FILETYPE_GRB, FILETYPE_GRB2
  USE mo_var_metadata_types,    ONLY: t_var_metadata
  USE mo_util_string,           ONLY: toupper
  USE mo_name_list_output_types,ONLY: t_output_name_list, t_output_file, &
    &                                 t_var_desc

  IMPLICIT NONE

  PUBLIC :: is_grib_output,                                  &
    &       is_variable_in_output_nml, is_variable_in_output
  PUBLIC :: use_async_name_list_io
  PUBLIC :: first_output_name_list
  PUBLIC :: add_var_desc


  ! Constant defining how many variable entries are added when resizing array:
  INTEGER, PARAMETER :: NVARS_GROW = 10

  ! Pointer to a linked list of output name lists:
  TYPE(t_output_name_list), POINTER :: first_output_name_list => NULL()

  ! Flag whether async name_list I/O is used, it is set in the main program:
  LOGICAL :: use_async_name_list_io = .FALSE.

CONTAINS
  
  !-------------------------------------------------------------------------------------------------
  !>
  !! @return .TRUE. if one of the output namelists has been specified with GRIB output.
  FUNCTION is_grib_output() RESULT(retval)
    LOGICAL                           :: retval
    TYPE(t_output_name_list), POINTER :: p_onl

    retval = .FALSE.
    p_onl => first_output_name_list
    DO
      IF(.NOT.ASSOCIATED(p_onl)) EXIT
      retval = retval                           .OR.  &
        &      (p_onl%filetype == FILETYPE_GRB) .OR.  &
        &      (p_onl%filetype == FILETYPE_GRB2)
      p_onl => p_onl%next
    ENDDO
  END FUNCTION is_grib_output


  !-------------------------------------------------------------------------------------------------
  !>
  !! @return .TRUE. if output for a given variable is due for a given
  !!         namelist
  !!
  FUNCTION is_variable_in_output_nml(p_onl, var_name) RESULT(retval)
    LOGICAL :: retval

    TYPE(t_output_name_list), POINTER         :: p_onl      !< output name list
    CHARACTER(LEN=*), INTENT(IN) :: var_name   !< variable name
    ! local variables
    INTEGER :: ivar
   
    ! if a specific variable name has been provided, loop over the
    ! variables for this output file
    retval = .FALSE.
    DO ivar=1,max_var_ml
      IF (p_onl%ml_varlist(ivar) == ' ') CYCLE
      IF (toupper(TRIM(p_onl%ml_varlist(ivar))) == toupper(TRIM(var_name))) retval=.TRUE.
      IF (retval) EXIT
    END DO
    DO ivar=1,max_var_pl
      IF (p_onl%pl_varlist(ivar) == ' ') CYCLE
      IF (toupper(TRIM(p_onl%pl_varlist(ivar))) == toupper(TRIM(var_name))) retval=.TRUE.
      IF (retval) EXIT
    END DO
    DO ivar=1,max_var_hl
      IF (p_onl%hl_varlist(ivar) == ' ') CYCLE
      IF (toupper(TRIM(p_onl%hl_varlist(ivar))) == toupper(TRIM(var_name))) retval=.TRUE.
      IF (retval) EXIT
    END DO
    DO ivar=1,max_var_il
      IF (p_onl%il_varlist(ivar) == ' ') CYCLE
      IF (toupper(TRIM(p_onl%il_varlist(ivar))) == toupper(TRIM(var_name))) retval=.TRUE.
      IF (retval) EXIT
    END DO
  END FUNCTION is_variable_in_output_nml


  !-------------------------------------------------------------------------------------------------
  !>
  !! @return .TRUE. if output for a given variable is due in any
  !!         output namelist.
  !!
  FUNCTION is_variable_in_output(first_output_name_list, var_name) RESULT(retval)
    LOGICAL                           :: retval

    TYPE(t_output_name_list), POINTER          :: first_output_name_list   !< head output namelist list
    CHARACTER(LEN=*), INTENT(IN)  :: var_name   !< variable name
    ! local variables
    TYPE (t_output_name_list), POINTER :: p_onl

    p_onl => first_output_name_list
    retval = .FALSE.
    DO
      IF(.NOT.ASSOCIATED(p_onl)) EXIT
      IF (retval) EXIT
      retval = is_variable_in_output_nml(p_onl, var_name=var_name)
      p_onl => p_onl%next
    END DO
  END FUNCTION is_variable_in_output


  !------------------------------------------------------------------------------------------------
  !> Append variable descriptor to the end of a (dynamically growing) list
  !! 
  !! @author  F. Prill, DWD
  SUBROUTINE add_var_desc(p_of, var_desc)
    TYPE(t_output_file), INTENT(INOUT)        :: p_of       !< output file
    TYPE(t_var_desc),    INTENT(IN)           :: var_desc   !< variable descriptor
    ! local variables
    CHARACTER(*), PARAMETER :: routine = TRIM("mo_name_list_output_config:add_var_desc")
    INTEGER                       :: errstat, new_max_vars, i, ivar
    TYPE(t_var_desc), ALLOCATABLE :: tmp(:)

    ! increase number of variables currently in use:
    p_of%num_vars = p_of%num_vars + 1
    IF (p_of%num_vars > p_of%max_vars) THEN
      ! array full, enlarge and make a triangle copy:
      new_max_vars = p_of%max_vars + NVARS_GROW
      IF (p_of%max_vars > 0) THEN
        ALLOCATE(tmp(p_of%max_vars), STAT=errstat)
        IF (errstat /= 0)  CALL finish (routine, 'Error in ALLOCATE operation!')
        tmp(1:p_of%max_vars) = p_of%var_desc(1:p_of%max_vars)
        DEALLOCATE(p_of%var_desc, STAT=errstat)
        IF (errstat /= 0)  CALL finish (routine, 'Error in DEALLOCATE operation!')
      END IF

      ALLOCATE(p_of%var_desc(new_max_vars), STAT=errstat)
      IF (errstat /= 0)    CALL finish (routine, 'Error in ALLOCATE operation!')
      ! Nullify pointers in p_of%var_desc
      DO ivar=(p_of%max_vars+1),new_max_vars
        p_of%var_desc(ivar)%r_ptr => NULL()
        p_of%var_desc(ivar)%s_ptr => NULL()
        p_of%var_desc(ivar)%i_ptr => NULL()
        DO i = 1, max_time_levels
          p_of%var_desc(ivar)%tlev_rptr(i)%p => NULL()
          p_of%var_desc(ivar)%tlev_sptr(i)%p => NULL()
          p_of%var_desc(ivar)%tlev_iptr(i)%p => NULL()
        ENDDO
      END DO

      IF (p_of%max_vars > 0) THEN
        p_of%var_desc(1:p_of%max_vars) = tmp(1:p_of%max_vars)
        DEALLOCATE(tmp, STAT=errstat)
        IF (errstat /= 0)  CALL finish (routine, 'Error in DEALLOCATE operation!')
      END IF
      p_of%max_vars = new_max_vars
    END IF
    ! add new element to array
    p_of%var_desc(p_of%num_vars) = var_desc
  END SUBROUTINE add_var_desc


END MODULE mo_name_list_output_config
