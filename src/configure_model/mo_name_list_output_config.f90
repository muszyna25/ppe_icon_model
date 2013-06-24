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

  USE mo_kind,                  ONLY: wp, i8
  USE mo_exception,             ONLY: finish
  USE mo_io_units,              ONLY: filename_max
  USE mo_impl_constants,        ONLY: max_phys_dom, max_bounds,          &
    &                                 vname_len, max_var_ml, max_var_pl, &
    &                                 max_var_hl, max_var_il, max_levels,&
    &                                 MAX_TIME_LEVELS
  USE mo_cdi_constants,         ONLY: FILETYPE_GRB, FILETYPE_GRB2
  USE mo_var_metadata,          ONLY: t_var_metadata
  USE mo_util_string,           ONLY: toupper
  USE mo_master_control,        ONLY: is_restart_run

  IMPLICIT NONE

  PUBLIC :: is_grib_output, &
    &       is_output_nml_active,  is_any_output_nml_active, &
    &       is_output_file_active, is_any_output_file_active
  PUBLIC :: use_async_name_list_io, l_output_phys_patch
  PUBLIC :: max_var_ml, max_var_pl, max_var_hl, max_bounds              
  PUBLIC :: max_levels, vname_len, t_output_name_list
  PUBLIC :: first_output_name_list, max_time_levels, t_output_file,&
    &       t_var_desc,t_rptr_5d
  PUBLIC :: add_var_desc

  CHARACTER(len=*),PARAMETER,PRIVATE :: &
    &  version = '$Id$'

  ! Flag whether async name_list I/O is used, it is set in the main program:

  LOGICAL :: use_async_name_list_io = .FALSE.
  
  ! Constant defining how many variable entries are added when resizing array:
  INTEGER, PARAMETER :: NVARS_GROW = 10

  ! The following parameter decides whether physical or logical patches are output
  ! and thus whether the domain number in output name lists pertains to physical
  ! or logical patches.

  LOGICAL, PARAMETER :: l_output_phys_patch = .TRUE. !** DO NOT CHANGE - needed for GRIB output **!

  TYPE t_output_name_list

    ! --------------------
    ! file name and format
    ! --------------------

    INTEGER                     :: filetype          ! One of CDI's FILETYPE_XXX constants
    CHARACTER(LEN=filename_max) :: output_filename   ! output filename prefix
    CHARACTER(LEN=filename_max) :: filename_format   ! output filename format (contains keywords <physdom>,<levtype> etc.)

    ! --------------------
    ! general settings
    ! --------------------

    INTEGER          :: mode                        ! 1 = forecast mode, 2 = climate mode
    INTEGER          :: dom(max_phys_dom)           ! domains for which this namelist is used, ending with -1
    INTEGER          :: output_time_unit            ! 1 = second, 2=minute, 3=hour, 4=day, 5=month, 6=year
    INTEGER          :: steps_per_file              ! Max number of output steps in one output file
    LOGICAL          :: include_last                ! Flag whether to include the last timestep in output
    LOGICAL          :: output_grid                 ! Flag whether grid information is output (in NetCDF output)

    ! post-processing times in units defined by output_time_unit: start, end, increment:
    REAL(wp)         :: output_bounds(3,max_bounds) 

    INTEGER          :: taxis_tunit   ! 1 = TUNIT_SECOND, 2 = TUNIT_MINUTE, 3 TUNIT_HOUR ... (see cdi.inc)

    ! --------------------
    ! ready file handling
    ! --------------------

    LOGICAL                     :: lwrite_ready     ! Flag. TRUE if a "ready file" (sentinel file) should be written
    CHARACTER(LEN=filename_max) :: ready_directory  ! output directory for ready files

    ! --------------------
    ! variable lists
    ! --------------------

    CHARACTER(LEN=vname_len)  :: ml_varlist(max_var_ml)   ! name of model level fields
    CHARACTER(LEN=vname_len)  :: pl_varlist(max_var_pl)   ! name of pressure level fields
    CHARACTER(LEN=vname_len)  :: hl_varlist(max_var_hl)   ! name of height level fields
    CHARACTER(LEN=vname_len)  :: il_varlist(max_var_hl)   ! name of isentropic level fields

    ! --------------------
    ! horizontal interpol.
    ! --------------------

    INTEGER  :: remap                 ! interpolate horizontally, 0: none, 1: to regular lat-lon grid, 2: to Gaussian grids, (3:...)
    LOGICAL  :: remap_internal        ! do interpolations online in the model or external (including triggering)
    INTEGER  :: lonlat_id             ! if remap=1: index of lon-lat-grid in global list "lonlat_grid_list"

    ! --------------------
    ! vertical interpol.
    ! --------------------

    REAL(wp) :: p_levels(max_levels)  ! pressure levels [hPa]
    REAL(wp) :: h_levels(max_levels)  ! height levels
    REAL(wp) :: i_levels(max_levels)  ! isentropic levels

    ! -------------------------------------
    ! Internal members, not read from input
    ! -------------------------------------

    INTEGER  :: cur_bounds_triple     ! current output_bounds triple in use
    REAL(wp) :: next_output_time      ! next output time (in seconds simulation time)
    INTEGER  :: n_output_steps
    TYPE(t_output_name_list), POINTER :: next ! Pointer to next output_name_list

  END TYPE t_output_name_list

  ! Pointer to a linked list of output name lists:
  TYPE(t_output_name_list), POINTER :: first_output_name_list => NULL()

  !------------------------------------------------------------------------------------------------

  ! Unfortunately, Fortran does not allow arrays of pointers, so we
  ! have to define extra types
  TYPE t_rptr_5d
    REAL(wp), POINTER :: p(:,:,:,:,:)
  END TYPE

  TYPE t_iptr_5d
    INTEGER,  POINTER :: p(:,:,:,:,:)
  END TYPE

  TYPE t_var_desc
    REAL(wp), POINTER :: r_ptr(:,:,:,:,:)         ! Pointer to time level independent REAL data (or NULL)
    INTEGER,  POINTER :: i_ptr(:,:,:,:,:)         ! Pointer to time level independent INTEGER data (or NULL)
    TYPE(t_rptr_5d) :: tlev_rptr(MAX_TIME_LEVELS) ! Pointers to time level dependent REAL data
    TYPE(t_iptr_5d) :: tlev_iptr(MAX_TIME_LEVELS) ! Pointers to time level dependent INTEGER data
    TYPE(t_var_metadata) :: info                  ! Info structure for variable
  END TYPE

  !------------------------------------------------------------------------------------------------
  TYPE t_output_file
    
    ! The following data must be set before opening the output file:
    CHARACTER(LEN=filename_max) :: filename_pref ! Prefix of output file name
    INTEGER                     :: output_type   ! CDI format
    INTEGER                     :: phys_patch_id ! ID of physical output patch
    INTEGER                     :: log_patch_id  ! ID of logical output patch
    REAL(wp)                    :: start_time    ! start time of model domain
    REAL(wp)                    :: end_time      ! end time of model domain
    LOGICAL                     :: initialized   ! .TRUE. if vlist setup has already been called
    INTEGER                     :: ilev_type     ! level type: level_type_ml/level_type_pl/level_type_hl/level_type_il
    INTEGER                     :: max_vars      ! maximum number of variables allocated
    INTEGER                     :: num_vars      ! number of variables in use
    TYPE(t_var_desc),ALLOCATABLE :: var_desc(:)
    TYPE(t_output_name_list), POINTER :: name_list ! Pointer to corresponding output name list

    CHARACTER(LEN=vname_len), ALLOCATABLE :: name_map(:,:) ! mapping internal names -> names in NetCDF

    INTEGER                     :: remap         ! Copy of remap from associated namelist

    INTEGER                     :: io_proc_id    ! ID of process doing I/O on this file

    !----------------------------
    ! Used for async IO only
    INTEGER(i8)                 :: my_mem_win_off
    INTEGER(i8), ALLOCATABLE    :: mem_win_off(:)
    !----------------------------

    ! The following members are set during open
    CHARACTER(LEN=filename_max) :: filename           ! Actual name of output file
    CHARACTER(LEN=filename_max) :: rdy_filename       ! Actual name of ready file (if any)
    INTEGER                     :: cdiFileId
    INTEGER                     :: cdiVlistId         ! cdi vlist handler
    INTEGER                     :: cdiCellGridID
    INTEGER                     :: cdiVertGridID
    INTEGER                     :: cdiEdgeGridID
    INTEGER                     :: cdiLonLatGridID
    INTEGER                     :: cdiZaxisID(27) ! All types of possible Zaxis ID's
    INTEGER                     :: cdiTaxisID
    INTEGER                     :: cdiTimeIndex
    INTEGER                     :: cdiInstID      ! output generating institute
    INTEGER                     :: cdi_grb2(3,2)  !< geographical position: (GRID, latitude/longitude)

  END TYPE t_output_file

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
