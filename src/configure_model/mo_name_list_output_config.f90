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
!! @par Copyright
!! 2002-2011 by DWD and MPI-M
!! This software is provided for non-commercial use only.
!! See the LICENSE and the WARRANTY conditions.
!!
!! @par License
!! The use of ICON is hereby granted free of charge for an unlimited time,
!! provided the following rules are accepted and applied:
!! <ol>
!! <li> You may use or modify this code for your own non commercial and non
!!    violent purposes.
!! <li> The code may not be re-distributed without the consent of the authors.
!! <li> The copyright notice and statement of authorship must appear in all
!!    copies.
!! <li> You accept the warranty conditions (see WARRANTY).
!! <li> In case you intend to use the code commercially, we oblige you to sign
!!    an according license agreement with DWD and MPI-M.
!! </ol>
!!
!! @par Warranty
!! This code has been tested up to a certain level. Defects and weaknesses,
!! which may be included in the code, do not establish any warranties by the
!! authors.
!! The authors do not make any warranty, express or implied, or assume any
!! liability or responsibility for the use, acquisition or application of this
!! software.
!!
MODULE mo_name_list_output_config

  USE mo_kind,                  ONLY: wp
  USE mo_exception,             ONLY: finish
  USE mo_impl_constants,        ONLY: max_var_ml, max_var_pl, &
    &                                 max_var_hl, max_var_il, &
    &                                 MAX_TIME_LEVELS
  USE mo_cdi_constants,         ONLY: FILETYPE_GRB, FILETYPE_GRB2
  USE mo_var_metadata,          ONLY: t_var_metadata
  USE mo_util_string,           ONLY: toupper
  USE mo_master_control,        ONLY: is_restart_run
  USE mo_name_list_output_types,ONLY: t_output_name_list, t_output_file, &
    &                                 t_var_desc

  IMPLICIT NONE

  PUBLIC :: is_grib_output,                                  &
    &       is_output_nml_active,  is_any_output_nml_active, &
    &       is_output_file_active, is_any_output_file_active
  PUBLIC :: use_async_name_list_io
  PUBLIC :: first_output_name_list
  PUBLIC :: add_var_desc

  CHARACTER(len=*),PARAMETER,PRIVATE :: &
    &  version = '$Id$'

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
  !! @return .TRUE. if output is due for a given namelist
  FUNCTION is_output_nml_active(p_onl, sim_time, dtime, &
    &                           iadv_rcf, last_step, is_restart, var_name) RESULT(retval)
    LOGICAL                           :: retval

    TYPE(t_output_name_list), POINTER         :: p_onl      !< output name list
    REAL(wp),            INTENT(IN), OPTIONAL :: sim_time   !< elapsed simulation time
    REAL(wp),            INTENT(IN)           :: dtime      !< [s] length of a time step
    INTEGER,             INTENT(IN)           :: iadv_rcf   !< calling freq. of adv., phys.
    LOGICAL,             INTENT(IN), OPTIONAL :: last_step
    LOGICAL,             INTENT(IN), OPTIONAL :: is_restart
    CHARACTER(LEN=*),    INTENT(IN), OPTIONAL :: var_name   !< variable name
    ! local variables
    INTEGER :: ivar
!    LOGICAL :: a
    
    IF (PRESENT(last_step)) THEN
      retval = (p_onl%include_last .AND. last_step)
    ELSE
      retval = .FALSE.
    END IF

    IF (PRESENT(sim_time)) THEN
      retval = retval .OR.  &
        &      ( p_onl%next_output_time <= sim_time+REAL(iadv_rcf,wp)*dtime/2._wp )
    END IF

    ! do nothing on the first timestep during a restart
    IF (PRESENT(is_restart)) THEN
      retval = retval .AND. .NOT. is_first_timestep_during_restart(is_restart,p_onl%next_output_time)
    ENDIF

    ! if a specific variable name has been provided, loop over the
    ! variables for this output file
    IF (PRESENT(var_name)) THEN
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
    END IF

  END FUNCTION is_output_nml_active


  !-------------------------------------------------------------------------------------------------
  !>
  !! @return .TRUE. if output is due for a given namelist
  FUNCTION is_any_output_nml_active(first_output_name_list, sim_time, dtime, &
    &                               iadv_rcf, last_step, is_restart, var_name) RESULT(retval)
    LOGICAL                           :: retval

    TYPE(t_output_name_list), POINTER          :: first_output_name_list   !< head output namelist list
    REAL(wp),            INTENT(IN), OPTIONAL  :: sim_time   !< elapsed simulation time
    REAL(wp),            INTENT(IN)            :: dtime      !< [s] length of a time step
    INTEGER,             INTENT(IN)            :: iadv_rcf   !< calling freq. of adv., phys.
    LOGICAL,             INTENT(IN), OPTIONAL  :: last_step
    LOGICAL,             INTENT(IN), OPTIONAL  :: is_restart
    CHARACTER(LEN=*),    INTENT(IN), OPTIONAL  :: var_name   !< variable name
    ! local variables
    TYPE (t_output_name_list), POINTER :: p_onl

    p_onl => first_output_name_list
    retval = .FALSE.

    DO
      IF(.NOT.ASSOCIATED(p_onl)) EXIT
      IF (retval) EXIT
      retval = is_output_nml_active(p_onl, sim_time, dtime, &
        &                           iadv_rcf, last_step, is_restart, var_name=var_name)
      p_onl => p_onl%next
    END DO

  END FUNCTION is_any_output_nml_active


  !-------------------------------------------------------------------------------------------------
  !>
  !! @return .TRUE. if output is due for a given file
  !!         (and, optionally, a given logical domain ID and/or a variable name)
  !! @author  F. Prill, DWD
  FUNCTION is_output_file_active(of, sim_time, dtime,             &
    &                            iadv_rcf, last_step, is_restart, &
    &                            idom, var_name) RESULT(retval)
    LOGICAL                           :: retval

    TYPE(t_output_file), INTENT(IN), TARGET   :: of         !< output file
    REAL(wp),            INTENT(IN), OPTIONAL :: sim_time   !< elapsed simulation time
    REAL(wp),            INTENT(IN)           :: dtime      !< [s] length of a time step
    INTEGER,             INTENT(IN)           :: iadv_rcf   !< calling freq. of adv., phys.
    LOGICAL,             INTENT(IN), OPTIONAL :: last_step
    LOGICAL,             INTENT(IN), OPTIONAL :: is_restart
    INTEGER,             INTENT(IN), OPTIONAL :: idom       !< logical domain index 
    CHARACTER(LEN=*),    INTENT(IN), OPTIONAL :: var_name   !< variable name

    TYPE (t_var_metadata), POINTER :: info
    INTEGER :: iv

    ! if a specific domain index has been provided,
    ! check if it matches the output file's logical domain ID:
    IF (PRESENT(idom)) THEN
      retval = (of%log_patch_id == idom)
    ELSE
      retval = .TRUE.
    END IF

    ! check if output file is active
    IF (PRESENT(sim_time)) THEN
      retval = retval .AND. &
        &      is_output_nml_active(of%name_list, sim_time, dtime, iadv_rcf, last_step,is_restart)
    END IF
    
    ! if a specific variable name has been provided, loop over the
    ! variables for this output file
    IF (retval .AND. PRESENT(var_name)) THEN
      retval = .FALSE.
      DO iv = 1, of%num_vars
        IF (retval) EXIT
        info => of%var_desc(iv)%info
        retval = (TRIM(info%name) == TRIM(var_name))
      END DO ! iv
    END IF

    IF (PRESENT(sim_time)) THEN
      IF (sim_time < of%start_time .OR. sim_time > of%end_time) retval = .FALSE.
    END IF

    !skip the initial state during a restarted run
    IF (PRESENT(is_restart)) THEN
      retval = retval .AND. .NOT. is_first_timestep_during_restart(is_restart,of%name_list%next_output_time)
    ENDIF

  END FUNCTION is_output_file_active


  !-------------------------------------------------------------------------------------------------
  !>
  !! @return .TRUE. if output is due for any output file
  !! @author  F. Prill, DWD
  FUNCTION is_any_output_file_active(of_list, sim_time, dtime, &
    &                                iadv_rcf, last_step, is_restart, idom, var_name) RESULT(retval)
    LOGICAL  :: retval

    TYPE(t_output_file), TARGET               :: of_list(:) !< list of output files
    REAL(wp),            INTENT(IN), OPTIONAL :: sim_time   !< elapsed simulation time
    REAL(wp),            INTENT(IN)           :: dtime      !< [s] length of a time step
    INTEGER,             INTENT(IN)           :: iadv_rcf   !< calling freq. of adv., phys.
    LOGICAL,             INTENT(IN), OPTIONAL :: last_step
    LOGICAL,             INTENT(IN), OPTIONAL :: is_restart
    INTEGER,             INTENT(IN), OPTIONAL :: idom       !< logical domain index 
    CHARACTER(LEN=*),    INTENT(IN), OPTIONAL :: var_name   !< variable name

    INTEGER :: i

    retval = .FALSE.
    DO i = 1, SIZE(of_list)
      IF (retval) EXIT
      retval = is_output_file_active(of_list(i), sim_time, dtime, &
        &                            iadv_rcf, last_step, is_restart, idom=idom, var_name=var_name)
      
    END DO ! iv
  END FUNCTION is_any_output_file_active


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
        p_of%var_desc(ivar)%i_ptr => NULL()
        DO i = 1, max_time_levels
          p_of%var_desc(ivar)%tlev_rptr(i)%p => NULL()
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

  ! to get the initial setup written out to the def. output file, the first
  ! values of output bounds is set to zero. to avoid this feature in a restart
  ! run is the purpose of this function
  FUNCTION is_first_timestep_during_restart(is_restart, next_output_time) RESULT(retval)
    LOGICAL :: retval

    LOGICAL,  INTENT(IN) :: is_restart
    REAL(wp), INTENT(IN) :: next_output_time

    retval = .FALSE.
    retval = (is_restart .AND. (ABS(next_output_time) < 0.05_wp))
  END FUNCTION is_first_timestep_during_restart

END MODULE mo_name_list_output_config
