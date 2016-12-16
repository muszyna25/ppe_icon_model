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
  USE mo_name_list_output_types,ONLY: t_output_name_list

  IMPLICIT NONE

  PUBLIC :: is_grib_output,                                  &
    &       is_variable_in_output
  PUBLIC :: use_async_name_list_io
  PUBLIC :: first_output_name_list

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
    CHARACTER(len=LEN(var_name)) :: var_name_u

    ! if a specific variable name has been provided, loop over the
    ! variables for this output file
    retval = .FALSE.
    var_name_u = toupper(var_name)
    DO ivar=1,max_var_ml
      IF (p_onl%ml_varlist(ivar) /= ' ') THEN
        retval = toupper(TRIM(p_onl%ml_varlist(ivar))) == var_name_u
        IF (retval) RETURN
      END IF
    END DO
    DO ivar=1,max_var_pl
      IF (p_onl%pl_varlist(ivar) /= ' ') THEN
        retval = toupper(TRIM(p_onl%pl_varlist(ivar))) == var_name_u
        IF (retval) RETURN
      END IF
    END DO
    DO ivar=1,max_var_hl
      IF (p_onl%hl_varlist(ivar) /= ' ') THEN
        retval = toupper(TRIM(p_onl%hl_varlist(ivar))) == var_name_u
        IF (retval) RETURN
      END IF
    END DO
    DO ivar=1,max_var_il
      IF (p_onl%il_varlist(ivar) /= ' ') THEN
        retval = toupper(TRIM(p_onl%il_varlist(ivar))) == var_name_u
        IF (retval) RETURN
      END IF
    END DO
  END FUNCTION is_variable_in_output_nml


  !-------------------------------------------------------------------------------------------------
  !>
  !! @return .TRUE. if output for a given variable is due in any
  !!         output namelist.
  !!
  FUNCTION is_variable_in_output(first_output_name_list, var_name) RESULT(retval)
    LOGICAL                           :: retval

    !> head output namelist list
    TYPE(t_output_name_list), POINTER :: first_output_name_list
    CHARACTER(LEN=*), INTENT(IN)  :: var_name   !< variable name
    ! local variables
    TYPE (t_output_name_list), POINTER :: p_onl
    INTEGER :: tlen

    tlen = LEN_TRIM(var_name)
    p_onl => first_output_name_list
    retval = .FALSE.
    DO WHILE (ASSOCIATED(p_onl) .AND. .NOT. retval)
      retval = is_variable_in_output_nml(p_onl, var_name=var_name(1:tlen))
      p_onl => p_onl%next
    END DO
  END FUNCTION is_variable_in_output

END MODULE mo_name_list_output_config
